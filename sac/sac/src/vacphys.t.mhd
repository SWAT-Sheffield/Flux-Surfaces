!##############################################################################
! module vacphys - mhd

INCLUDE:vacphys.mhd0.t

!=============================================================================

subroutine keeppositive(ix^L,w)

! Keep pressure and density positive (following D.Ryu)

include 'vacdef.f'

integer::          ix^L
double precision:: w(ixG^T,nw)
logical:: toosmallp
!-----------------------------------------------------------------------------


   ! Where rho is small use vacuum state: rho=vacuumrho, v=0, p=smallp, same B
   where((w(ix^S,rho_)+w(ix^S,rhob_))<smallrho)
      ^C&w(ix^S,m^C_)=zero;
!!!      ^C&w(ix^S,m^C_)=w(ix^S,m^C_)/w(ix^S,rho_)*vacuumrho;
      w(ix^S,rho_)=vacuumrho-w(ix^S,rhob_)
      w(ix^S,e_)=smallp/(eqpar(gamma_)-one)+half*(^C&w(ix^S,b^C_)**2+)-w(ix^S,eb_)
   endwhere


! Calculate pressure without clipping toosmall values (.false.)
call getpthermal(w,ix^L,tmp)

toosmallp=any(tmp(ix^S)<max(zero,smallp))

if(toosmallp)then
   nerror(toosmallp_)=nerror(toosmallp_)+1
   if(nerror(toosmallp_)==1)then
      write(*,'(a,i2,a,i7)')&
         'Too small pressure (code=',toosmallp_,') at it=',it
      write(*,*)'Value < smallp: ',minval(tmp(ix^S)),smallp
!     write(*,*)'Location: ',minloc(tmp(ix^S)) !F77_  
      {write(*,*)'Processor:',ipe ^IFMPI}
   endif
   if(smallp>zero)&
      w(ix^S,e_)=max(tmp(ix^S),smallp)/(eqpar(gamma_)-1)+&
         half*((^C&w(ix^S,m^C_)**2+)/w(ix^S,rho_)+(^C&w(ix^S,b^C_)**2+))-w(ix^S,eb_)
endif

return
end

!=============================================================================
! end module vacphys - mhd
!##############################################################################
