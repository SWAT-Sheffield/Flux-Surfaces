!#############################################################################
! module vaccd
! Centered difference scheme
!=============================================================================

!=============================================================================
subroutine centdiff4(qdt,ixI^L,ixO^L,iws,idim^LIM,qtC,wCT,qt,w)

! Advance the iws flow variables from t to t+qdt within ixO^L by 
! fourth order centered  differencing in space the dw/dt+dF_i(w)/dx_i=S 
! type equation. 
! wCT contains the time centered variables at time qtC for flux and source.
! w is the old value at qt on input and the new value at qt+qdt on output.

include 'vacdef.f'

double precision:: qdt,qtC,qt,wCT(ixG^T,nw),w(ixG^T,nw)
integer:: ixI^L,ixO^L,iws(niw_),idim^LIM
logical :: transport

double precision:: v(ixG^T),f(ixG^T), fb(ixG^T)
integer:: iiw,iw,ix^L,idim,idir
!-----------------------------------------------------------------------------


! Two extra layers are needed in each direction for which fluxes are added.
ix^L=ixO^L;
do idim= idim^LIM
   ix^L=ix^L^LADD2*kr(idim,^D);
enddo
if(ixI^L^LTix^L|.or.|.or.) call die( &
   'Error in CentDiff4: Non-conforming input limits')

! Add fluxes to w
do idim= idim^LIM
   ix^L=ixO^L^LADD2*kr(idim,^D);

   call getv(wCT,ix^L,idim,v)

   do iiw=1,iws(niw_); iw=iws(iiw)
!   print*,'iiw', iiw,idim,idir
      ! Get non-transported flux
      call getflux(wCT,ix^L,iw,idim,f,transport)

      ! Add transport flux
      if(transport)f(ix^S)=f(ix^S)+v(ix^S)*wCT(ix^S,iw)

      ! Add divergence of flux
      call gradient4(.false.,f,ixO^L,idim,tmp)
      w(ix^S,iw)=w(ix^S,iw)-qdt*tmp(ix^S)

   select case(iw)

    case(e_)

         call gradient4(.false.,v,ixO^L,idim,tmp)   
         call getptotal_bg(w,ix^L,fb)
         
         w(ix^S,iw)=w(ix^S,iw)-qdt*tmp(ix^S)*fb(ix^S)
	 
	do idir= idim^LIM 
             call gradient4(.false.,v,ixO^L,idir,tmp)   
             w(ix^S,iw)=w(ix^S,iw)+qdt*w(ix^S,bg0_+idir)*w(ix^S,bg0_+idim)*tmp(ix^S)
        enddo

   endselect        


   end do    !next iw
end do       !next idim


if(sourceunsplit) &
   call addsource2(qdt*(idimmax-idimmin+one)/ndim, &
      ixI^L,ixO^L,iws,qtC,wCT,qt,w)

return
end

!=============================================================================
! end module vaccd
!#############################################################################
