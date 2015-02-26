!#############################################################################
! module vacgrid
! Subroutines for boundaries, grid, divergence of flux, gradients
! Also the limiter functions for TVD, TVDLF, TVDMU schemes

{INCLUDE:vacgrid.gencoord.t ^IFGEN}
!INCLUDE:vacgrid.setnozzle.t
!=============================================================================
subroutine boundsetup

! Variables describing boundaries
!
!    iB=1..nB                                    - boundary ID
!    ixBmin(idim,iB):ixBmax(idim,iD)             - location of boundary
!    idimB(iB)                                   - direction orthogonal to the 
!                                                  boundary layer
!    upperB(iB)                                  - .true. if boundary is at a
!                                                  high index of the phys. grid
!    typeB(iw,iB)                                - boundary type string

include 'vacdef.f'

integer:: ix^L,iB,jB,iw,idim,idm,ixG^LIM(ndim),ixM^LIM(ndim)
!-----------------------------------------------------------------------------

oktest=index(teststr,'boundsetup')>=1
if(oktest)write(*,*)'BoundSetup'

! If ixBmax(1,1)==0, the user did not specify boundary regions. Setup default.
! A schematic picture for ndim=2, where the numbers are iB for each region, and
! *-s are the extra layers for qx. 
!
!                     ************
!                     *1444444424*    ixGmax2
!                     *4144444442*
!                     *11      22*    ixMmax2
!                     *11      22*
!                     *11      22*    ixMmin2
!                     *1333333323*
!                     *3133333332*    ixGmin2
!                     ************
!                      i i    i i
!                      x x    x x
!                      G M    M G
!                      m m    m m
!                      i i    a a
!                      n n    x x
!                      1 1    1 1

^D&ixG^LIM(^D)=ixG^DL;
^D&ixM^LIM(^D)=ixM^DL;

if(ixBmax(1,1)==0)then
   nB=2*ndim
   do iB=1,nB
      idim=(iB+1)/2
      idimB(iB) = idim
      upperB(iB)= 2*idim==iB
      do idm=1,ndim
         ixB^LIM(idm,iB)=ixG^LIM(idm);
      end do
      if(upperB(iB))then
         ixBmin(idim,iB)=ixMmax(idim)+1
         ixBmax(idim,iB)=ixGmax(idim)
      else
         ixBmin(idim,iB)=ixGmin(idim)
         ixBmax(idim,iB)=ixMmin(idim)-1
      end if
   end do
else
   ! Check if boundary regions are inside grid and in between mesh and grid
   do iB=1,nB
      ix^L=ixB^LIM(^D,iB);
      {^IFMPI
      ! Convert global limits into local grid limits 
      ix^L=ix^L-ixpemin^D+1;
      ! Cut off extra cells at MPI boundaries
      ^D&if(ipe^D>0)ixmin^D=max(ixmin^D,ixGmin^D)\
      ^D&if(ipe^D<npe^D-1)ixmax^D=min(ixmax^D,ixGmax^D)\
      ! Change the global ixB^LIM(^D,iB) to the local index limits
      ^D&ixBmin(^D,iB)=ixmin^D;
      ^D&ixBmax(^D,iB)=ixmax^D;
      ! Check if local boundary region is empty
      if(ixmin^D>ixmax^D|.or.)then
          ! The limits do not contain cells, change the boundary type to 'mpi'
          typeB(1:nw,iB)='mpi'
          CYCLE
      endif
      |\}
      if(ixmin^D<ixGmin^D.or.ixmax^D>ixGmax^D|.or.) then
         write(*,*)'Error for boundary region iB,ixL=',iB,ix^L
         call die('Error in BoundSetup: Boundary region is outside grid')
      endif
      select case(idimB(iB))
      {case(^D)
         if(upperB(iB))then
            if(ixmin^D-1/=ixMmax^D.or.ixmax^D/=ixGmax^D)write(*,*)&
               'Warning in BoundSetup: Boundary does not fit, iB:',iB
         else
            if(ixmax^D+1/=ixMmin^D.or.ixmin^D/=ixGmin^D)write(*,*)&
               'Warning in BoundSetup: Boundary does not fit, iB:',iB
         endif \}
      end select
   end do
end if

! Identify the periodic pairs if they are not defined in boundlist
! Check type, direction, orientation, and size before matching pairs
do iB=1,nB
   if(typeB(1,iB)=='periodic'.and.ipairB(iB)==0)then
      do jB=iB+1,nB
         if(typeB(1,jB)=='periodic'.and.ipairB(jB)==0.and.&
         idimB(iB)==idimB(jB).and.(upperB(iB).neqv.upperB(jB)).and.&
         {ixBmax(^D,iB)-ixBmin(^D,iB)==ixBmax(^D,jB)-ixBmin(^D,jB)|.and.})then
            ipairB(iB)=jB; ipairB(jB)=iB
         endif
      end do
      if(ipairB(iB)==0)call die('Error in BoundSetup: No periodic pair')
   end if
end do

{^IFMPI
do iB=1,nB
   ! Change boundary type if processor is not at the edge or if this is a
   ! periodic boundary and there are more than 1 processors in this direction
   select case(idimB(iB))
   {case(^D)
      if(typeB(1,iB)=='periodic' .and. npe^D>1) then
          typeB(1:nw,iB)='mpiperiod'
      else if(upperB(iB) .and. ipe^D<npe^D-1 .or. &
              .not.upperB(iB) .and. ipe^D>0) then
          typeB(1:nw,iB)='mpi'
      endif \}
   end select
end do
}


! symm0 means zero orthogonal flux via the boundary 
do iw=1,nw
   do iB=1,nB
     if(typeB(iw,iB)=='symm0') nofluxB(iw,idimB(iB))=.true.
   enddo
enddo

if(oktest)then
   do iB=1,nB
      {write(unitterm,*)'ipe=',ipe ^IFMPI}
      write(unitterm,*)'iB,idimB,upperB,typeB:',iB,idimB(iB),upperB(iB),&
         ' ',typeB(1,iB)
      write(unitterm,*)'ixBmin:',(ixBmin(idm,iB),idm=1,ndim)
      write(unitterm,*)'ixBmax:',(ixBmax(idm,iB),idm=1,ndim)
   end do
   do iB=1,nB
      if(typeB(1,iB)=='periodic')write(*,*)'Pairs:',iB,ipairB(iB)
   end do
end if

return
end

!=============================================================================
subroutine ensurebound(dix,ixI^L,ixO^L,qt,w)

! Check if there is enough information for calculating derivatives.
! The requirement is that ixI is wider than ixO by dix. 
! Adjust ixI and ixO. Call getboundary if needed.

include 'vacdef.f'

integer:: dix,ixI^L,ixO^L
double precision:: qt,w(ixG^T,nw)
!-----------------------------------------------------------------------------

oktest=index(teststr,'ensurebound')>0
if(oktest)write(*,*)'EnsureBound dix,ixI,ixO:',dix,',',ixI^L,',',ixO^L

! Check wether ixO+dix is within the grid
if(ixG^L^LTixO^L^LADDdix|.or.|.or.)then
   ixO^L=ixM^L;
endif
! Check whether ixI is wider than ixO by at least dix otherwise getboundary
if(ixI^L^LTixO^L^LADDdix|.or.|.or.)then
   ixI^L=ixG^L;
   call getboundary(qt,1,nw,1,ndim,w)
   if(oktest)write(*,*)'calling getboundary'
end if

if(oktest)write(*,*)'Final       dix,ixI,ixO:',dix,',',ixI^L,',',ixO^L

return
end

!=============================================================================
subroutine getboundary(qt,iw^LIM,idim^LIM,w)

include 'vacdef.f'

integer:: iw^LIM,idim^LIM
double precision:: qt,w(ixG^T,1:nw)
integer:: ix,ix^D,ixe,ixf,ix^L,ixpair^L,idim,iw,iB
integer:: iwv,jdim
double precision:: coeffnormal,coefftransv
logical:: initialized

integer:: ixee
integer:: ixb^L

data initialized/.false./
!-----------------------------------------------------------------------------


ixb^L=ixG^LL;

oktest=index(teststr,'getboundary')>=1
if(oktest)write(*,*)'GetBoundary it,step:',it,step
if(oktest)write(*,*)'GetBoundary wold:',w(ixtest^D,iwtest)

if(extraB)call specialbound(qt,ixM^L,0,0,w)
if(smallfix)call keeppositive(ixM^L,w)

{^IFMPI
! Get boundaries from other PE-s
call mpibound(nw,w)
}

iB=0
do
   iB=iB+1
   if(oktest)write(*,*)'iB  :',iB
   if(iB>nB) exit
   idim=idimB(iB)
   ! Only boundary segments in the required direction(s) are filled in
   if(idimmin>idim.or.idimmax<idim) cycle

   ix^L=ixB^LIM(^D,iB);

   ! The possibly shifted coordinates parallel to the boundary layer 
   ! are defined by the PAIR of the periodic/overlapping boundary.
   ! Put the location of the source of information into ixpairL.
   if(ipairB(iB)>0)then
      ixpair^L=ixB^LIM(^D,ipairB(iB));
      select case(idim)
      {case(^D)
         if(upperB(iB))then
            ixpair^LIM^D=ixpair^LIM^D+dixB^LIM^D;
         else
            ixpair^LIM^D=ixpair^LIM^D-dixB^LIM^D;
         endif
      \}
      end select
   endif

   do iw= iw^LIM
      if(oktest)write(*,*)'  iw:',iw
      if(oktest)write(*,*)'typeB(iw,iB):',typeB(iw,iB)
      select case (typeB(iw,iB))



      case('contCD4')
         select case(idim)
      {   case(^D)

            if(upperB(iB))then

      	    
               ixe=ixmin^D-2
	       ixee=ixmin^D-3

               !HPF$ INDEPENDENT
               ix= ixmin^D
	          w(ix^D%ix^S,iw)=w(ixe^D%ix^S,iw)

               ix= ixmax^D
		  w(ix^D%ix^S,iw)=w(ixee^D%ix^S,iw)
	       
            else
               ixe=ixmax^D+3
	       ixee=ixmax^D+2

               !HPF$ INDEPENDENT
               ix= ixmin^D
                  w(ix^D%ix^S,iw)=w(ixe^D%ix^S,iw)

               ix= ixmax^D
		  w(ix^D%ix^S,iw)=w(ixee^D%ix^S,iw)
	       
            endif


          }
         end select

      case('zero')
         select case(idim)
      {   case(^D)

            if(upperB(iB))then
              
	      call primitive(ixG^LL,w)
	      
               ixe=ixmin^D-2
	       ixee=ixmin^D-3

	       ixbmax^D=ixmax^D
	       ixbmin^D=ixee
	       
	                     
               !HPF$ INDEPENDENT
               ix= ixmin^D
	       
                  w(ix^D%ix^S,iw)=w(ixe^D%ix^S,iw)

               ix= ixmax^D
		  w(ix^D%ix^S,iw)=w(ixee^D%ix^S,iw)
	       
	      call conserve(ixG^LL,w)	         
	       
           else
	       
	      call primitive(ixG^LL,w)
               
               ixe=ixmax^D+3
	       ixee=ixmax^D+2
	       
       
               !HPF$ INDEPENDENT
               ix= ixmin^D
                  w(ix^D%ix^S,iw)=w(ixe^D%ix^S,iw)

               ix= ixmax^D
		  w(ix^D%ix^S,iw)=w(ixee^D%ix^S,iw)
	

	      call conserve(ixG^LL,w)	         	       
	       
            endif
	
   	

          }
	  
	  
	  
         end select
     


      case('cont','fixed')
         ! For 'cont' copy w at the edge into the whole boundary region.
         ! Fot 'fixed' copy w first, later use values stored in fixB.
         ! For fullgridini=T store the w values read from the file in fixB.
         select case(idim)
         {case(^D)
            if(upperB(iB))then
               ixe=ixmin^D-1
            else
               ixe=ixmax^D+1
            endif
            if(fixedB(iw,iB))then
               !HPF$ INDEPENDENT
               {do ix^DD=ixmin^DD,ixmax^DD\}
                  w(ix^DD,iw)=fixB^D(ix^D-ixe^D%ix^DD,iw)
               {enddo^DD&\} 
            else if(typeB(iw,iB)=='cont' .or. .not.fullgridini) then
               !HPF$ INDEPENDENT
                do ix= ix^DL
                  w(ix^D%ix^S,iw)=w(ixe^D%ix^S,iw)
               end do
            end if \}
         end select
      case ('cont1','fixed1','grad1')
         ! First order extrapolation from edge and inner edge to the boundary
         ! 'cont1'  extrapolates in every time step, can be numericly unstable.
         ! 'fixed1' extrapolates first, stores VALUES into fixB, then restores.
         ! 'grad1'  extrapolates first, stores DIFFERENCES into fixB, then 
         !          adds the stored differences to the current edge value.
         select case(idim)
         {case(^D)
            if(upperB(iB))then
               ixe=ixmin^D-1; ixf=ixe-1
            else
               ixe=ixmax^D+1; ixf=ixe+1
            endif
            if(fixedB(iw,iB))then
               if(typeB(iw,iB)=='grad1')then
                  !HPF$ INDEPENDENT
                  {do ix^DD=ixmin^DD,ixmax^DD\}
                     w(ix^DD,iw)=fixB^D(ix^D-ixe^D%ix^DD,iw)+w(ixe^D%ix^DD,iw)
                  {enddo^DD&\}
               else
                  !HPF$ INDEPENDENT
                  {do ix^DD=ixmin^DD,ixmax^DD\}
                     w(ix^DD,iw)=fixB^D(ix^D-ixe^D%ix^DD,iw)
                  {enddo^DD&\}
               endif
            else if(typeB(iw,iB)=='cont1'.or. .not.fullgridini)then !HPF_
            !HPF_ endif
            !HPF_ if(.not.fixedB(iw,iB).and.&
            !HPF_    (typeB(iw,iB)=='cont1'.or. .not.fullgridini))then
               !HPF$ INDEPENDENT
               do ix= ix^DL
                  w(ix^D%ix^S,iw)=&
                     (abs(ix-ixe)+1)*w(ixe^D%ix^S,iw)- &
                      abs(ix-ixe)   *w(ixf^D%ix^S,iw)  
               end do
            end if \}
         end select
      case('periodic')
         ! Update boundary by translation of w by width of mesh (and shift)
         w(ix^S,iw)=w(ixpair^S,iw)
      case('symm','symm0','asymm')
         ! Reflect w into the boundary region, multiply by -1 for "asymm"
         ! In generalized coordinates take into account the other vector
         ! components for vector variables. The symmetry of the transverse 
         ! component is based on typeB for the jdim=idim+1 -th component.
         ! ixe is used for the reflection, normal vectors are taken at ixf+1/2
         if(gencoord.and.vectoriw(iw)>=0)then
            ! Determine direction of vector component, symmetry coefficients
            ! for normal and transverse vector components
            iwv=vectoriw(iw); jdim=idim+1-(idim/ndim)*ndim 
            coeffnormal=1; if(typeB(iwv+idim,iB)=='asymm') coeffnormal=-1
            coefftransv=1; if(typeB(iwv+jdim,iB)=='asymm') coefftransv=-1
         endif
         select case(idim)
         {case(^D)
            if(upperB(iB))then
               ixe=2*ixmin^D-1; ixf=ixmin^D-1
            else
               ixe=2*ixmax^D+1; ixf=ixmax^D
            endif
            if(gencoord.and.vectoriw(iw)>=0)then
               !HPF$ INDEPENDENT
               do ix= ix^DL
                  w(ix^D%ix^S,iw)=zero
                  do jdim=1,ndim
                     w(ix^D%ix^S,iw)=w(ix^D%ix^S,iw)+&
                        normalC(ixf^D%ix^S,idim,jdim)*w(ixe-ix^D%ix^S,iwv+jdim)
                  end do
                  w(ix^D%ix^S,iw)=w(ix^D%ix^S,iw)*&
                     normalC(ixf^D%ix^S,idim,iw-iwv)*(coeffnormal-coefftransv)&
                     +w(ixe-ix^D%ix^S,iw)*coefftransv
               end do
            else
               !HPF$ INDEPENDENT
               do ix= ix^DL
                  w(ix^D%ix^S,iw)=w(ixe-ix^D%ix^S,iw)
               end do
               if(typeB(iw,iB)=='asymm') w(ix^S,iw)=-w(ix^S,iw)
            endif  
         \}
         end select
      case('special')
            ! Skip special now, we do it after normal boundary type variables
            !HPF_ if(.false.)write(*,*)'Avoiding xlhpf90 compiler bug'
      {^IFMPI
      case('mpi','mpiperiod')
            ! This boundary is handled by MPI \}
      case default
         write(uniterr,*)'Error in GetBoundary: boundary type(', &
            iw,iB,')=',typeB(iw,iB),' is not implemented!'
         call die(' ')
      end select ! typeB(iw,iB)
   end do ! next iw
   ! Do special boundaries
   do iw= iw^LIM
      if(oktest)write(*,*)'special, iw:',iw
      if(typeB(iw,iB).eq.'special')call specialbound(qt,ix^L,iw,iB,w)
   end do ! next iw
end do ! next iB

! Fixed boundaries (fixed,fixed1) or fixed gradients (grad1) are stored into
! fixB after all boundaries have been updated.
! This needs to be done in the very first time step only.
if(.not.initialized)then
   initialized=.true.
   do iB= 1,nB
      ix^L=ixB^LIM(^D,iB);
      do iw= iw^LIM
         if( (typeB(iw,iB)=='fixed'.or.typeB(iw,iB)=='fixed1'.or.&
              typeB(iw,iB)=='grad1') .and. .not.fixedB(iw,iB))then
            fixedB(iw,iB)=.true.
            select case(idimB(iB))
               {case(^D)
                  if(upperB(iB))then
                     ixe=ixmin^D-1
                  else
                     ixe=ixmax^D+1
                  endif
                  if(typeB(iw,iB)=='grad1')then
                     !HPF$ INDEPENDENT
                     {do ix^DD= ixmin^DD,ixmax^DD\}
                      fixB^D(ix^D-ixe^D%ix^DD,iw)=w(ix^DD,iw)-w(ixe^D%ix^DD,iw)
                     {enddo^DD&\}
                  else
                     !HPF$ INDEPENDENT
                     {do ix^DD= ixmin^DD,ixmax^DD\}
                      fixB^D(ix^D-ixe^D%ix^DD,iw)=w(ix^DD,iw)
                     {enddo^DD&\}
                  endif
               \}
            end select
         end if 
      end do ! iw
   end do    ! iB
end if

if(oktest)write(*,*)'GetBoundary wnew:',w(ixtest^D,iwtest)

return 
end

!=============================================================================
subroutine setnoflux(iw,idim,ix^L,fRC,ixR^L,fLC,ixL^L)

! Set flux in direction idim to zero for variable iw if typeB is 'symm0'
! in a boundary region

include 'vacdef.f'

integer:: iw,idim,ix^L,ixL^L,ixR^L
double precision:: fRC(ixG^T), fLC(ixG^T)
integer:: iB,ixe,ixB^L

!-----------------------------------------------------------------------------

oktest=index(teststr,'setnoflux')>=1

do iB=1,nB
   if(typeB(iw,iB)=='symm0'.and.idimB(iB)==idim)then
      ixB^L=ixB^LIM(^D,iB);
      ! Calculate edge index and set the appropriate flux to zero
      select case(idim)
      {case(^D)
      if(upperB(iB))then
         ixe=ixBmin^D-1+ixRmin^D-ixmin^D
         if(oktest)write(*,*)'Setnoflux it,idim,iw,ixe:', &
             it,idim,iw,ixe,fRC(ixe^D%ixtest^DD)
         fRC(ixe^D%ixB^S)=zero
      else
         ixe=ixBmax^D+1+ixLmin^D-ixmin^D
         if(oktest)write(*,*)'Setnoflux it,idim,iw,ixe:',&
             it,idim,iw,ixe,fLC(ixe^D%ixtest^DD)
         fLC(ixe^D%ixB^S)=zero
      endif
      \}
      end select
   endif
enddo

return
end

!=============================================================================
subroutine gridsetup1

! Cartesian or polar grid. Determine x at the boundaries.
! Determine often needed combinations of x, such as dx or dvolume.
! Determine variables for axial symmetry
!
! ixe          - edge coordinate of the grid touching the boundary region
! ixf          - coordinate inside of ixe
! qx           - x with an extended index range for calculation of dx

include 'vacdef.f'

integer:: ix^L,hx^L,jx^L
integer:: ix,ixe,ixf,idim,jdim
double precision:: qx(IXG^LL^LADD1:,ndim)
double precision:: r(IXGLO1-1:IXGHI1+1),rC(IXGLO1-1:IXGHI1+1)

!-----------------------------------------------------------------------------

oktest=index(teststr,'gridsetup')>=1
if(oktest)write(*,*)'GridSetup1'

! Calculate qx in the boundary layers by linearly extrapolating x from
! the touching edge of the grid (ixe) and the inner edge (ixf).

{^IFMPI
! Fill in ghost cells for x from MPI neighbors
call mpibound(ndim,x)
}
qx(ixG^LL^LADD1:,1:ndim)=zero
qx(ixG^S,1:ndim) = x(ixG^S,1:ndim)
do idim=1,ndim
   ix^L=ixG^L^LADD1;
   select case(idim)
      {case(^D)
         ! First do the upper layers
         ixmax^D=ixGmax^D+1; ixmin^D=ixMmax^D+1
         if(fullgridini .or.mpiupperB(^D)^IFMPI) ixmin^D=ixGmax^D+1
         ixe=ixmin^D-1; ixf=ixe-1
         !!!forall replaced by do for sake of sequentializing (Adaptor)
	 do jdim=1,ndim
            do ix= ix^DL
               qx(ix^D%ix^S,jdim)=(1+abs(ixe-ix))*qx(ixe^D%ix^S,jdim)- &
                                     abs(ixe-ix) *qx(ixf^D%ix^S,jdim)
            end do
         end do
         ! Next do the lower layers
         ixmin^D=ixGmin^D-1; ixmax^D=ixMmin^D-1
         if(fullgridini .or.mpilowerB(^D)^IFMPI) ixmax^D=ixGmin^D-1
         ixe=ixmax^D+1; ixf=ixe+1
	 do jdim=1,ndim
            do ix= ix^DL
                qx(ix^D%ix^S,jdim)=(1+abs(ixe-ix))*qx(ixe^D%ix^S,jdim)- &
                                      abs(ixe-ix) *qx(ixf^D%ix^S,jdim)
            end do
         end do \}
   end select
enddo

x(ixG^S,1:ndim)=qx(ixG^S,1:ndim)

! Loop with ^D instead of idim to avoid an SP2 xlphf90 compiler error
! if qx is distributed. But it should not be
!{
!^D&jx^L=ixG^L+kr(^D,^DD); hx^L=ixG^L-kr(^D,^DD);
!   dx(ixG^S,^D)=half*(qx(jx^S,^D)-qx(hx^S,^D))
!   if(oktest)write(*,*)'dx,qxj,qxh:',dx(ixtest^DD,^D),&
!      qx(ixtest^DD+kr(^D,^DD),^D),qx(ixtest^DD-kr(^D,^DD),^D)
!\}

do idim=1,ndim
   jx^L=ixG^L+kr(idim,^D); hx^L=ixG^L-kr(idim,^D);
   dx(ixG^S,idim)=half*(qx(jx^S,idim)-qx(hx^S,idim))
end do

! Calculate geometrical factors for axial symmetry based on Boris FCT.
! Fluxes are multiplied by areaC. The cell volume is proportional to areadx.
! Gradient type source terms are multiplied by areaside = darea/areadx.

if(oktest)write(*,*)'Start calculating geometrical factors'



if(oktest)write(*,*)'Start calculating cell volumes'

! Calculate volume of cells and total volume of mesh
if(typeaxial=='slab')then
   dvolume(ixG^S)= ^D&dx(ixG^S,^D)*
else
   forall(ix= ixG^LIM1:) dvolume(ix,ixG^SE)= ^D&areadx(ix)^%1dx(ix,ixG^SE,^D)*
endif
volume=sum(dvolume(ixM^S))
{^IFMPI
! Add up volumes
call mpiallreduce(volume,MPI_SUM)
! Correct volumes in 2nd ghost cells from neighboring processors
call mpibound(1,dvolume)
}

! For polar grid dx_phi=r*dphi. 
if(polargrid)dx(ixG^S,pphi_)=x(ixG^S,r_)*dx(ixG^S,pphi_)

if(oktest)write(*,*)'Finish GridSetup1'

return 
end

!=============================================================================
 
 subroutine gradient4(realgrad,q,ix^L,idim,gradq)
 include 'vacdef.f'
 logical:: realgrad
 integer:: ix^L,idim
 double precision:: q(ixG^T),gradq(ixG^T)
 integer:: kx^L,jx^L,hx^L,gx^L
 integer:: minx1^D,maxx1^D,k
 !-----------------------------------------------------------------------------
 
 !SHIFT
 kx^L=ix^L+2*kr(idim,^D);
 !SHIFT MORE
 jx^L=ix^L+kr(idim,^D);
 !SHIFT MORE
 hx^L=ix^L-kr(idim,^D);
 !SHIFT MORE
 gx^L=ix^L-2*kr(idim,^D);
 
 !SHIFT BEGIN
 gradq(ix^S)=-(q(kx^S)-8.D0*(q(jx^S)-q(hx^S))-q(gx^S))/dx(ix^S,idim)/12.D0
 !SHIFT END
 
 minx1^D=ixmin^D+kr(idim,^D);
 maxx1^D=ixmax^D-kr(idim,^D);
 
 do k=0,1  !left-right bc
 
 if (typeB(1,2*idim-1+k) .ne. 'mpi') then
 if (upperB(2*idim-1+k)) then
 select case(idim)
 {   case(^D)
 gradq(ixmax^D^D%ix^S)=0.d0
 gradq(maxx1^D^D%ix^S)=0.d0
 }
 end select
 else
 select case(idim)
 {   case(^D)
 gradq(ixmin^D^D%ix^S)=0.d0
 gradq(minx1^D^D%ix^S)=0.d0
 }
 end select
 endif
 endif
 enddo
 
 return
 end


!=============================================================================
subroutine laplace4(q,ix^L,laplaceq)

!!! This subroutine should not use tmp or tmp2

! Calculate 4th order laplace of q within ixL in Cartesian direction idim

!!! We assume uniform Cartesian grid in slab symmetry for now

include 'vacdef.f'

integer:: ix^L
double precision:: q(ixG^T),laplaceq(ixG^T)

integer:: idim,kx^L,jx^L,hx^L,gx^L
!-----------------------------------------------------------------------------

if(gencoord)call die('Error: laplace4 does not work for gen.coords!')
if(typeaxial/='slab')&
   call die('Error: laplace4 does not work for axial symmetry yet!')

oktest= index(teststr,'lpalace')>=1

laplaceq(ix^S)=zero

do idim=1,ndim
   !SHIFT
   kx^L=ix^L+2*kr(idim,^D); 
   !SHIFT MORE
   jx^L=ix^L+kr(idim,^D); 
   !SHIFT MORE
   hx^L=ix^L-kr(idim,^D);
   !SHIFT MORE
   gx^L=ix^L-2*kr(idim,^D);

   !SHIFT BEGIN
   laplaceq(ix^S)=laplaceq(ix^S)+&
     (q(kx^S)+q(gx^S)+30*q(ix^S)-16*(q(jx^S)+q(hx^S)))/dx(ix^S,idim)**2/12
   !SHIFT END

   if(oktest)write(*,*)'idim,q(kx,jx,ix,hx,gx):',idim,&
           q(ixtest^D+2*kr(idim,^D)),q(ixtest^D+kr(idim,^D)),q(ixtest^D),&
           q(ixtest^D-kr(idim,^D)),q(ixtest^D-2*kr(idim,^D))
enddo

if(oktest)write(*,*)'laplaceq:',laplaceq(ixtest^D)

return
end

!=============================================================================
! end module vacgrid
!##############################################################################



