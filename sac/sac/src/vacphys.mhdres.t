!##############################################################################
! module vacphys.mhdres - subroutines for resistive mhd and mhdiso

!=============================================================================
subroutine getdt_res(w,ix^L)

! If resistivity is  not zero, check diffusion time limit for dt

include 'vacdef.f'

double precision:: w(ixG^T,nw),dtdiff
integer:: ix^L,idim,idirmin
save dtdiff

double precision:: current(ixG^T,7-2*ndir:3),eta(ixG^T),gradeta(ixG^T,ndim)
common/resist/current,eta,gradeta
!-----------------------------------------------------------------------------

oktest=index(teststr,'getdt')>=1
if(oktest)write(*,*)'GetDt_Res'

if(eqpar(eta_)==zero)return

if(eqpar(eta_)>zero)then
   dtdiff=dtdiffpar*minval(dx(ix^S,1:ndim))**2/eqpar(eta_)
else if(eqpar(eta_)<zero)then
   if(it==itmin)then
      call getcurrent(w,ixM^L,idirmin)
      call specialeta(w,ixM^L,idirmin)
   endif
   dtdiff=bigdouble
   do idim=1,ndim
      dtdiff=min(dtdiff,&
                 dtdiffpar/(smalldouble+maxval(eta(ix^S)/dx(ix^S,idim)**2)))
   enddo
endif
{^IFMPI call mpiallreduce(dtdiff,MPI_MIN)}

dt=min(dt,dtdiff)

if(oktest) write(*,*)'GetDt dtdiff:',dtdiff
if(oktest) write(*,*)'GetDt dt    :',dt

return
end

!=============================================================================
subroutine getcurrent(w,ix^L,idirmin)

! Calculate idirmin and the idirmin:3 components of the common current array

include 'vacdef.f'

integer, parameter:: idirmin0=7-2*ndir
double precision:: w(ixG^T,nw)
integer::          ix^L,idirmin

integer:: ixI^L,idir,jdir,kdir

! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for MHD
double precision:: current(ixG^T,7-2*ndir:3),eta(ixG^T),gradeta(ixG^T,ndim)
common/resist/current,eta,gradeta
!-----------------------------------------------------------------------------

oktest=index(teststr,'getcurrent')>=1
if(oktest)write(*,*)'GetCurrent'

ixI^L=ix^L^LADD1;

! Calculate current density within ixL: J=curl B, thus J_i=eps_ijk*d_j B_k
! Current can have components (idirmin0:3)
! Determine exact value of idirmin while doing the loop.

idirmin=4
current(ix^S,idirmin0:3)=zero

do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir
   if(lvc(idir,jdir,kdir)/=0)then
      tmp(ixI^S)=(w(ixI^S,b0_+kdir)+w(ixI^S,bg0_+kdir))
      call gradient(.true.,tmp,ix^L,jdir,tmp2)
      if(lvc(idir,jdir,kdir)==1)then
         current(ix^S,idir)=current(ix^S,idir)+tmp2(ix^S)
      else
         current(ix^S,idir)=current(ix^S,idir)-tmp2(ix^S)
      endif
      if(idir<idirmin)idirmin=idir
   endif
enddo; enddo; enddo;

if(oktest)write(*,*)'idirmin,J(idirmin:3):',idirmin,current(ixtest^D,idirmin:3)

return
end

!=============================================================================
subroutine addsource_res1(qdt,ixI^L,ix^L,iws,qtC,w,qt,wnew)

! Add resistive source to wnew within ixL if possible, otherwise shrink ixL
! Uses 3 point stencil (1 neighbour) in each direction, non-conservative

include 'vacdef.f'

integer::          ixI^L,ix^L,iws(niw_)
double precision:: qdt,qtC,qt,w(ixG^T,nw),wnew(ixG^T,nw)

integer:: ix,jx^L,hx^L,idim,idir,jdir,kdir,idirmin,iiw,iw

! Resistivity "eta" may or may not vary in time and/or space
! For ndir=2 only 3rd component of J can exist, ndir=1 is not possible for MHD
double precision:: current(ixG^T,7-2*ndir:3),eta(ixG^T),gradeta(ixG^T,ndim)
common/resist/current,eta,gradeta
!-----------------------------------------------------------------------------

oktest=index(teststr,'addsource_res')>=1
if(oktest)write(*,*)'AddSource_Res1'

! Compact resistive sources involve one extra layer only
call ensurebound(1,ixI^L,ix^L,qtC,w)

! Calculate current density within ixL and determine idirmin
call getcurrent(w,ix^L,idirmin)

! Calculate and save eta for the first time
! for eqpar(eta_)<0 call specialeta, this will also set the common gradeta
if(eqpar(eta_)>zero)then
   if(it==itmin) eta(ixG^S)=eqpar(eta_)
else
   call specialeta(w,ix^L,idirmin)
endif
if(oktest)write(*,*)'eta    :',eta(ixtest^D)
if(oktest)write(*,*)'gradeta:',gradeta(ixtest^D,1:ndim)

do idir=1,ndir

   ! Put B_idir into tmp2 and eta*Laplace B_idir into tmp
   tmp(ix^S)=zero
   tmp2(ixG^S)=(w(ixG^S,b0_+idir)+w(ixG^S,bg0_+idir))
   if(gencoord)then 
      ! Use contour integral of Grad(B_idir) along cell edges
      !!! Assumes that connected cell centers are orthogonal to interfaces
      do idim=1,ndim
         !SHIFT
         jx^L=ix^L+kr(idim,^D);
         !SHIFT MORE
         hx^L=ix^L-kr(idim,^D);
         !SHIFT BEGIN
         tmp(ix^S)=tmp(ix^S) &
           +surfaceC(ix^S,idim)*(tmp2(jx^S)-tmp2(ix^S)) &
                             /sqrt(^D&(x(jx^S,^D)-x(ix^S,^D))**2+) &
           +surfaceC(hx^S,idim)*(tmp2(hx^S)-tmp2(ix^S)) &
                             /sqrt(^D&(x(hx^S,^D)-x(ix^S,^D))**2+)
         !SHIFT END
      enddo
      tmp(ix^S)=tmp(ix^S)/dvolume(ix^S)
   else
      do idim=1,ndim
         if(typeaxial=='cylinder'.and.idim==r_)then
            ! Calculate 1/r d(r d(B_idir)/dr)/dr
            forall(ix=ixmin1:ixmax1)tmp(ix,ix^SE)=tmp(ix,ix^SE) &
              +(areaC(ix)*(tmp2(ix+1,ix^SE)-tmp2(ix,ix^SE))&
                           /(x(ix+1,ix^SE,r_)-x(ix,ix^SE,r_)) &
               +areaC(ix-1)*(tmp2(ix-1,ix^SE)-tmp2(ix,ix^SE))&
                           /(x(ix,ix^SE,r_)-x(ix-1,ix^SE,r_))) &
                             /x(ix,ix^SE,r_)/dx(ix,ix^SE,r_)
         else
            !SHIFT
            jx^L=ix^L+kr(idim,^D); 
            !SHIFT MORE
            hx^L=ix^L-kr(idim,^D);
            !SHIFT BEGIN
            tmp(ix^S)=tmp(ix^S)+&
              (tmp2(jx^S)-2*tmp2(ix^S)+tmp2(hx^S))/dx(ix^S,idim)**2 


            !SHIFT END
         endif
      enddo
   endif

   ! Multiply by eta
   tmp(ix^S)=tmp(ix^S)*eta(ix^S)

   ! Subtract grad(eta) x J = eps_ijk d_j eta J_k if eta is non-constant
   if(eqpar(eta_)<zero)then
      do jdir=1,ndim; do kdir=idirmin,3
         if(lvc(idir,jdir,kdir)/=0)then
            if(lvc(idir,jdir,kdir)==1)then
               tmp(ix^S)=tmp(ix^S)-gradeta(ix^S,jdir)*current(ix^S,kdir)
            else
               tmp(ix^S)=tmp(ix^S)+gradeta(ix^S,jdir)*current(ix^S,kdir)
            endif
         endif
      enddo; enddo
   endif

   ! Add sources related to eta*laplB-grad(eta) x J to B and e
   do iiw=1,iws(niw_); iw=iws(iiw)
      if(iw==b0_+idir)then
         ! dB_idir/dt+=tmp
         wnew(ix^S,iw)=wnew(ix^S,iw)+qdt*tmp(ix^S)
      else if(iw==e_)then
         ! de/dt+=B.tmp
         wnew(ix^S,iw)=wnew(ix^S,iw)+qdt*tmp(ix^S)*(w(ix^S,b0_+idir)+w(ix^S,bg0_+idir))
      endif
   end do  ! iiw
enddo ! idir

! de/dt+=eta*J**2
do iiw=1,iws(niw_); iw=iws(iiw)
   if(iw==e_)then
      tmp(ix^S)=zero
      do idir=idirmin,3
         tmp(ix^S)=tmp(ix^S)+current(ix^S,idir)**2
      enddo
      wnew(ix^S,iw)=wnew(ix^S,iw)+qdt*eta(ix^S)*tmp(ix^S)
   endif
enddo

return
end

!=============================================================================
subroutine addsource_res2(qdt,ixI^L,ixO^L,iws,qtC,w,qt,wnew)

! Add resistive source to wnew within ixO if possible, otherwise shrink ixO
! Uses 5 point stencil (2 neighbours) in each direction, conservative

include 'vacdef.f'

integer::          ixI^L,ixO^L,iws(niw_)
double precision:: qdt,qtC,qt,w(ixG^T,nw),wnew(ixG^T,nw)

integer:: ix^L,idir,jdir,kdir,idirmin,iiw,iw

! Resistivity may or may not vary in time and/or space
! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for MHD
double precision:: current(ixG^T,7-2*ndir:3),eta(ixG^T),gradeta(ixG^T,ndim)
common/resist/current,eta,gradeta
!-----------------------------------------------------------------------------

oktest=index(teststr,'addsource_res')>=1
if(oktest)write(*,*)'AddSource_Res2'

! Calculating resistive sources involves second derivatives, two extra layers
if(oktest)write(*,*)'calling ensurebound in Addsource_Res2'
call ensurebound(2,ixI^L,ixO^L,qtC,w)
if(oktest)write(*,*)'end calling ensurebound in Addsource_Res2'
ix^L=ixO^L^LADD1;

! Calculate current density within ixL: J=curl B, thus J_i=eps_ijk*d_j B_k
! Determine exact value of idirmin while doing the loop.
call getcurrent(w,ix^L,idirmin)

! Calculate and save eta for the first time
! for eqpar(eta_)<0 call specialeta
if(eqpar(eta_)>zero)then
   if(it==itmin) eta(ixG^S)=eqpar(eta_)
else
   call specialeta(w,ix^L,idirmin)
endif
if(oktest)write(*,*)'eta:',eta(ixtest^D)

! Calculate sources from resistivity
do iiw=1,iws(niw_); iw=iws(iiw)
   select case(iw)
   case(b^C_)
      ! dB/dt= -curl(J*eta), thus B_i=B_i-eps_ijk d_j Jeta_k
      idir=iw-b0_
      do jdir=1,ndim; do kdir=idirmin,3;
         if(lvc(idir,jdir,kdir)/=0)then
            tmp(ix^S)=current(ix^S,kdir)*eta(ix^S)*qdt
            call gradient(.true.,tmp,ixO^L,jdir,tmp2)
            if(lvc(idir,jdir,kdir)==1)then
               wnew(ixO^S,iw)=wnew(ixO^S,iw)-tmp2(ixO^S)
            else
               wnew(ixO^S,iw)=wnew(ixO^S,iw)+tmp2(ixO^S)
            endif
        endif
      enddo; enddo
   case(e_)
      ! de/dt+= div(B x Jeta), thus e-= dt*eps_ijk d_i B_j Jeta_k
      do idir=1,ndim; do jdir=1,ndir; do kdir=idirmin,3
         if(lvc(idir,jdir,kdir)/=0)then
            tmp(ix^S)=(w(ix^S,b0_+jdir)+w(ix^S,bg0_+jdir))*current(ix^S,kdir)*eta(ix^S)*qdt
            call gradient(.false.,tmp,ixO^L,idir,tmp2)
            if(lvc(idir,jdir,kdir)==1)then
               wnew(ixO^S,ee_)=wnew(ixO^S,ee_)+tmp2(ixO^S)
            else
               wnew(ixO^S,ee_)=wnew(ixO^S,ee_)-tmp2(ixO^S)
            endif
        endif
      enddo; enddo; enddo
   end select ! iw
end do        ! iiw

return
end

!=============================================================================
! end module vacphys.mhdres
!##############################################################################


subroutine gradient(bebe,q,ix^L,idim,gradq)
include 'vacdef.f'
integer:: ix^L,idim
double precision:: q(ixG^T),gradq(ixG^T)
integer:: hx^L,kx^L
logical:: bebe
!-----------------------------------------------------------------------------

hx^L=ix^L-kr(idim,^D);
kx^L=ix^L+kr(idim,^D);
gradq(ix^S)=-(q(kx^S)-q(hx^S))/dx(ix^S,idim)/two

return
end
