!##############################################################################
! module vacphys0 - mhd

!=============================================================================
subroutine conserve(ix^L,w)

! Transform primitive variables into conservative ones

include 'vacdef.f'

integer:: ix^L
double precision:: w(ixG^T,nw)
!-----------------------------------------------------------------------------


! Calculate total energy from pressure, kinetic and magnetic energy

w(ix^S,e_)=w(ix^S,p_)/(eqpar(gamma_)-1)+&
   half*((w(ix^S,rho_)+w(ix^S,rhob_))*(^C&w(ix^S,v^C_)**2+)+(^C&(w(ix^S,b^C_))**2+))+&
   (^C&(w(ix^S,b^C_)*w(ix^S,bg^C_))+)


! Convert velocity to momentum
^C&w(ix^S,m^C_)=(w(ix^S,rho_)+w(ix^S,rhob_))*w(ix^S,v^C_);


return
end

!=============================================================================
subroutine primitive(ix^L,w)

! Transform conservative variables into primitive ones

include 'vacdef.f'

integer:: ix^L
double precision:: w(ixG^T,nw)
!-----------------------------------------------------------------------------


! Calculate pressure

call getpthermal(w,ix^L,tmp)

w(ix^S,p_)=tmp(ix^S)

! Convert momentum to velocity
^C&w(ix^S,v^C_)=w(ix^S,m^C_)/(w(ix^S,rho_)+w(ix^S,rhob_));

return
end

!=============================================================================
subroutine getv(w,ix^L,idim,v)

! Calculate v_idim=m_idim/rho within ix

include 'vacdef.f'

integer:: ix^L,idim
double precision:: w(ixG^T,nw),v(ixG^T)
!-----------------------------------------------------------------------------

oktest=index(teststr,'getv')>=1
if(oktest)write(*,*)'GetV w:',w(ixtest^D,iwtest)

v(ix^S)=w(ix^S,m0_+idim)/(w(ix^S,rho_)+w(ix^S,rhob_))

if(oktest)write(*,*)'GetV v:',v(ixtest^D)

return 
end


!=============================================================================
subroutine getcmax(new_cmax,w,ix^L,idim,cmax)

! Calculate cmax_idim=cfast_i+abs(v_idim) within ix^L
! where cfast_i=sqrt(0.5*(cf**2+sqrt(cf**4-4*cs**2*b_i**2/rho)))
! and cf**2=b**2/rho+cs**2/rho is the square of the speed of the fast wave 
! perpendicular to the magnetic field, and cs is the sound speed.

include 'vacdef.f'

logical:: new_cmax
integer:: ix^L,idim
double precision:: w(ixG^T,nw),cmax(ixG^T)
double precision:: csound2(ixG^T),cfast2(ixG^T)
save csound2,cfast2
!-----------------------------------------------------------------------------
oktest=index(teststr,'getcmax')>=1

!Direction independent part of getcmax:
if(new_cmax)then
   new_cmax=.false.
   call getcsound2(w,ix^L,csound2)
   if(oktest)write(*,*)'csound2:',csound2(ixtest^D)
   cfast2(ix^S)=(^C&(w(ix^S,b^C_)+w(ix^S,bg^C_))**2+ )/(w(ix^S,rho_)+w(ix^S,rhob_))+csound2(ix^S)
end if
if(oktest)write(*,*)'cfast2:',cfast2(ixtest^D)

cmax(ix^S)=sqrt(half*(cfast2(ix^S)+ &
   sqrt(cfast2(ix^S)**2-4*csound2(ix^S)* &
   ((w(ix^S,b0_+idim)+w(ix^S,bg0_+idim))**2)/(w(ix^S,rho_)+w(ix^S,rhob_))))) &
   +abs(w(ix^S,m0_+idim)/(w(ix^S,rho_)+w(ix^S,rhob_)))

if(oktest) write(*,*)'cmax:',cmax(ixtest^D)


return 
end

!=============================================================================
subroutine getcsound2prim(w,ix^L,csound2)

! Calculate the square of the thermal sound speed csound2 within ix^L
! from the primitive variables in w.
! csound2=gamma*p/rho

include 'vacdef.f'

double precision:: w(ixG^T,nw),csound2(ixG^T)
integer:: ix^L
!-----------------------------------------------------------------------------

if(eqpar(gamma_)<=zero)&
   call die('Correct Getcsound2prim for NONIDEAL gas in vacphys.t.mhd')

csound2(ix^S)=eqpar(gamma_)*(w(ix^S,p_)+(eqpar(gamma_)-one)*(w(ix^S,eb_)-&
                 half*( ^C&(w(ix^S,bg^C_))**2+ )))/(w(ix^S,rho_)+w(ix^S,rhob_))

return 
end

!=============================================================================
subroutine getcsound2(w,ix^L,csound2)

! Calculate the square of the thermal sound speed csound2 within ix^L.
! csound2=gamma*p/rho

include 'vacdef.f'

double precision:: w(ixG^T,nw),csound2(ixG^T)
integer:: ix^L
!-----------------------------------------------------------------------------

if(eqpar(gamma_)<=zero)&
   call die('Correct Getcsound2 for NONIDEAL gas in vacphys.t.mhd')

oktest=index(teststr,'getcsound2')>=1
if(oktest) write(*,*)'Getcsound2'

call getpthermal(w,ix^L,csound2)
if(oktest) write(*,*)'p(ixtest)=',csound2(ixtest^D)
csound2(ix^S)=eqpar(gamma_)*(csound2(ix^S)+(eqpar(gamma_)-one)*(w(ix^S,eb_)-&
                 half*( ^C&(w(ix^S,bg^C_))**2+ )))/(w(ix^S,rho_)+w(ix^S,rhob_))

return 
end

!=============================================================================
subroutine getpthermal(w,ix^L,p)

!!! This subroutine should not use tmp,tmp2


include 'vacdef.f'

double precision:: w(ixG^T,nw),p(ixG^T)
integer:: ix^L
!-----------------------------------------------------------------------------


p(ix^S)=half*( ^C&w(ix^S,m^C_)**2+ )/(w(ix^S,rho_)+w(ix^S,rhob_))

p(ix^S)=p(ix^S)+ half*( ^C&(w(ix^S,b^C_)**2)+ )+( ^C&(w(ix^S,b^C_)*w(ix^S,bg^C_))+ )

p(ix^S)=(eqpar(gamma_)-one)*(w(ix^S,e_)-p(ix^S))


return 
end

!=============================================================================
subroutine getptotal(w,ix^L,p)

include 'vacdef.f'

double precision::  w(ixG^T,nw),p(ixG^T),gamma
integer:: ix^L
!-----------------------------------------------------------------------------

if(eqpar(gamma_)<=zero)&
   call die('Correct GetPtotal for NONIDEAL gas in vacphys.t.mhd')

gamma=eqpar(gamma_)

p(ix^S)=(gamma-two)*(( ^C&(w(ix^S,b^C_)*w(ix^S,bg^C_))+ )+ half*( ^C&(w(ix^S,b^C_))**2.d0+ ))

p(ix^S)=(gamma-one)*(w(ix^S,e_)-half*( ^C&w(ix^S,m^C_)**2.d0+ )/(w(ix^S,rho_)+w(ix^S,rhob_)))-p(ix^S)


return 
end

!=============================================================================
subroutine getptotal_bg(w,ix^L,p)

include 'vacdef.f'

double precision::  w(ixG^T,nw),p(ixG^T),gamma
integer:: ix^L
!-----------------------------------------------------------------------------

if(eqpar(gamma_)<=zero)&
   call die('Correct GetPtotal for NONIDEAL gas in vacphys.t.mhd')

gamma=eqpar(gamma_)

p(ix^S)=(eqpar(gamma_)-one)*w(ix^S,eb_)-half*(eqpar(gamma_)-two)*( ^C&(w(ix^S,bg^C_)**2.d0)+ )    

return 
end

