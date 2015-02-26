!##############################################################################
! module vacphys.mhd0 - common subroutines for mhd and mhdiso

INCLUDE:vacphys.mhdroe.t^IFTVD
INCLUDE:vacproc.projectb.t^NOONED^IFPOISSON
INCLUDE:vacproc.constrainb.t^NOONED^IFCT
INCLUDE:vacphys.mhdres.t^IFRES
!=============================================================================
subroutine physini

! Tell VAC which variables are vectors, set default entropy coefficients

include 'vacdef.f'
integer:: il
!-----------------------------------------------------------------------------

iw_vector(1)=m0_; iw_vector(2)=b0_

! The values of the constants are taken from Ryu & Jones ApJ 442, 228
do il=1,nw
   select case(il)
   case(fastRW_,fastLW_,slowRW_,slowLW_)
      entropycoef(il)=0.2
   case(alfvRW_,alfvLW_)
      entropycoef(il)=0.4
   case default
      entropycoef(il)= -one
   end select
end do

return
end

!=============================================================================
subroutine process(count,idim^LIM,w)

! Process w before it is advected in directions idim^LIM, or before save
! count=1 and 2 for first and second (half step) processing during advection
! count=ifile+2 for saving results into the file indexed by ifile

include 'vacdef.f'

integer:: count,idim^LIM
double precision:: w(ixG^T,nw)

logical:: oktime
double precision:: cputime,time1,timeproc
data timeproc /0.D0/

! The processing should eliminate divergence of B.
!-----------------------------------------------------------------------------

oktest=index(teststr,'process')>=1
oktime=index(teststr,'timeproc')>=1

if(oktest)write(*,*)'Process it,idim^LIM,count',it,idim^LIM,count

if(oktime)time1=cputime()

{^NOONED
if(count==0)then
   if(divbconstrain)then
      {^IFCT call constrainb(w)
      if(.false.)}call die('CT module is OFF: setvac -on=ct; make vac')
   endif
else
   ! Use the projection scheme 
   {^IFPOISSON call projectb(w)
   if(.false.)}call die('Poisson module is OFF: setvac -on=poisson;make vac')
endif
}

if(oktime)then
   time1=cputime()-time1
   timeproc=timeproc+time1
   write(*,*)'Time.Proc:',time1,timeproc
endif

return
end

!=============================================================================
subroutine getdt(w,ix^L)

! If resistivity is  not zero, check diffusion time limit for dt

include 'vacdef.f'

double precision:: w(ixG^T,nw)
integer:: ix^L
!-----------------------------------------------------------------------------

oktest=index(teststr,'getdt')>=1
if(oktest)write(*,*)'GetDt'

if(eqpar(eta_)==zero)return

{^IFRES call getdt_res(w,ix^L)}
{^NORES write(*,*)'Error: Resistive MHD module is OFF'
call die('Recompile with setvac -on=resist or set eqpar(eta_)=0')}

return
end

!=============================================================================
subroutine getdivb(w,ixO^L,divb)

! Calculate div B within ixO

include 'vacdef.f'

integer::          ixO^L,ix^L,idim
double precision:: w(ixG^T,nw),divb(ixG^T)
!-----------------------------------------------------------------------------

oktest=index(teststr,'getdivb')>=1
if(oktest)write(*,*)'getdivb ixO=',ixO^L

if(fourthorder)then
   ix^L=ixO^L^LADD2;
else
   ix^L=ixO^L^LADD1;
endif
divb(ixO^S)=zero
do idim=1,ndim
   tmp(ix^S)=w(ix^S,b0_+idim)+w(ix^S,bg0_+idim)

      call gradient4(.false.,tmp,ixO^L,idim,tmp2)

   divb(ixO^S)=divb(ixO^S)+tmp2(ixO^S)
enddo

if(oktest)then
   write(*,*)'divb:',divb(ixtest^D)
!   write(*,*)'bx=',w(ixtest1-1:ixtest1+1,ixtest2,b1_)
!   write(*,*)'by=',w(ixtest1,ixtest2-1:ixtest2+1,b2_)
!   write(*,*)'x =',x(ixtest1-1:ixtest1+1,ixtest2,1)
!   write(*,*)'y =',x(ixtest1,ixtest2-1:ixtest2+1,2)
!   write(*,*)'dx=',dx(ixtest1,ixtest2,1)
!   write(*,*)'dy=',dx(ixtest1,ixtest2,2)
endif

return
end

!=============================================================================
subroutine getflux(w,ix^L,iw,idim,f,transport)

! Calculate non-transport flux f_idim[iw] within ix^L.
! Set transport=.true. if a transport flux should be added

include 'vacdef.f'

integer::          ix^L,iw,idim
double precision:: w(ixG^T,nw),f(ixG^T), fb(ixG^T)
logical::          transport
!-----------------------------------------------------------------------------

oktest= index(teststr,'getflux')>=1
if(oktest.and.iw==iwtest)write(*,*)'Getflux idim,w:',&
   idim,w(ixtest^D,iwtest)

transport=.true.

     
select case(iw)
   case(rho_)
      f(ix^S)=w(ix^S,rhob_)*w(ix^S,m0_+idim)/(w(ix^S,rho_)+w(ix^S,rhob_))
      
   {case(m^C_)
      if(idim==^C)then
            
          call getptotal(w,ix^L,f)
	  call getptotal_bg(w,ix^L,fb)
          fb(ix^S)=0.d0
	  
          f(ix^S)=f(ix^S)+fb(ix^S)-(w(ix^S,b^C_)*w(ix^S,bg0_+idim)+w(ix^S,b0_+idim)*w(ix^S,bg^C_))-&
                           w(ix^S,b^C_)*w(ix^S,b0_+idim) !-&
			 !  w(ix^S,bg0_+idim)*w(ix^S,bg^C_)    !remove for perturbed
      else
          f(ix^S)=-(w(ix^S,b^C_)*w(ix^S,bg0_+idim)+w(ix^S,b0_+idim)*w(ix^S,bg^C_))-&
                    w(ix^S,b^C_)*w(ix^S,b0_+idim) !-&
		 !   -w(ix^S,bg0_+idim)*w(ix^S,bg^C_)  !remove for perturbed
 
      endif \}

   case(e_)
  
      call getptotal(w,ix^L,f)
      call getptotal_bg(w,ix^L,fb)      
      fb(ix^S)=0.d0      

      f(ix^S)=(w(ix^S,m0_+idim)*(f(ix^S)+fb(ix^S))-&
               w(ix^S,b0_+idim)*( ^C&(w(ix^S,bg^C_))*w(ix^S,m^C_)+ )-&
               w(ix^S,bg0_+idim)*( ^C&(w(ix^S,b^C_))*w(ix^S,m^C_)+ ))/(w(ix^S,rho_)+w(ix^S,rhob_))+&	      
               w(ix^S,eb_)*w(ix^S,m0_+idim)/(w(ix^S,rho_)+w(ix^S,rhob_))-&
	       w(ix^S,b0_+idim)*( ^C&(w(ix^S,b^C_))*w(ix^S,m^C_)+ )/(w(ix^S,rho_)+w(ix^S,rhob_)) !&
              !  -w(ix^S,bg0_+idim)*( ^C&(w(ix^S,bg^C_))*w(ix^S,m^C_)+ )/(w(ix^S,rho_)+w(ix^S,rhob_))  ! remove for perturbed
	       
	      


   {case(b^C_)
      if(idim==^C) then
         f(ix^S)= zero
         transport=.false.
      else
      
         f(ix^S)= -w(ix^S,m^C_)/(w(ix^S,rho_)+w(ix^S,rhob_))*(w(ix^S,b0_+idim)+w(ix^S,bg0_+idim))+ &
                  w(ix^S,m0_+idim)/(w(ix^S,rho_)+w(ix^S,rhob_))*w(ix^S,bg^C_)
		  
      endif  \}

   case default
      call die('Error in getflux: unknown flow variable')
end select

if(oktest.and.iw==iwtest)write(*,*)'transport,f:',&
   transport,f(ixtest^D)

return
end

!=============================================================================
subroutine addsource(qdt,ixI^L,ixO^L,iws,qtC,w,qt,wnew)

! Add sources from resistivity and Powell solver

include 'vacdef.f'

integer::          ixI^L,ixO^L,iws(niw_)
double precision:: qdt,qtC,qt,w(ixG^T,nw),wnew(ixG^T,nw)
!-----------------------------------------------------------------------------

oktest=index(teststr,'addsource')>=1
if(oktest)write(*,*)'Addsource, compactres,divbfix:',compactres,divbfix
if(oktest)write(*,*)'Before adding source:',wnew(ixtest^D,iwtest)

! Sources for resistivity in eqs. for e, B1, B2 and B3
if(abs(eqpar(eta_))>smalldouble)then
   {^IFRES
   if(compactres)then
      call addsource_res1(qdt,ixI^L,ixO^L,iws,qtC,w,qt,wnew)
   else
      call addsource_res2(qdt,ixI^L,ixO^L,iws,qtC,w,qt,wnew)
   endif
   if(oktest)write(*,*)'With resistive source:',wnew(ixtest^D,iwtest)
   }
   {^NORES write(*,*)'Error: Resistive MHD module is OFF'
   call die('Recompile with setvac -on=resist or set eqpar(eta_)=0')}
endif


! Sources related to div B in the Powell solver
{^NOONED if(divbfix) call addsource_divb(qdt,ixI^L,ixO^L,iws,qtC,w,qt,wnew)}

if(oktest)write(*,*)'After adding source:',wnew(ixtest^D,iwtest)

return
end

!=============================================================================
subroutine addsource_divb(qdt,ixI^L,ixO^L,iws,qtC,w,qt,wnew)

! Add Powell's divB related sources to wnew within ixO if possible, 
! otherwise shrink ixO

include 'vacdef.f'

integer::          ixI^L,ixO^L,iws(niw_),iiw,iw
double precision:: qdt,qtC,qt,w(ixG^T,nw),wnew(ixG^T,nw)
double precision:: divb(ixG^T)
!-----------------------------------------------------------------------------

! Calculating div B involves first derivatives
call ensurebound(1,ixI^L,ixO^L,qtC,w)

! We calculate now div B
call getdivb(w,ixO^L,divb)
divb(ixO^S)=qdt*divb(ixO^S)

do iiw=1,iws(niw_); iw=iws(iiw)
   select case(iw)
      {case(m^C_)
         wnew(ixO^S,iw)=wnew(ixO^S,iw)-(w(ixO^S,b^C_)+w(ixO^S,bg^C_))*divb(ixO^S)
      }
      {case(b^C_)
         wnew(ixO^S,iw)=wnew(ixO^S,iw)-w(ixO^S,m^C_)/(w(ixO^S,rho_)+w(ixO^S,rhob_))*divb(ixO^S)
      }
      case(e_)
         wnew(ixO^S,iw)=wnew(ixO^S,iw)-&
             (^C&w(ixO^S,m^C_)*(w(ixO^S,b^C_)+w(ixO^S,bg^C_))+ )/(w(ixO^S,rho_)+w(ixO^S,rhob_))*divb(ixO^S)
   end select
end do

return
end


!=============================================================================
! end module vacphys.mhd0
!##############################################################################
