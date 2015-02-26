  !##############################################################################
  ! module vacusr - vert
  
  
  INCLUDE:vacusr.gravity.t
  INCLUDE:vacusr.viscosity.t

  !=============================================================================
  SUBROUTINE specialini(ix^L,w)

    INCLUDE 'vacdef.f'

    INTEGER:: ix^L
    DOUBLE PRECISION:: w(ixG^T,1:nw)

    RETURN
  END SUBROUTINE specialini


  !=============================================================================
  SUBROUTINE specialsource(qdt,ixI^L,ixO^L,iws,qtC,wCT,qt,w)

    INCLUDE 'vacdef.f'

    INTEGER:: ixI^L,ixO^L,iws(niw_)
    DOUBLE PRECISION:: qdt,qtC,qt,wCT(ixG^T,nw),w(ixG^T,nw) !Intent(IN)

    INTEGER:: iw
    INTEGER:: ix_1,ix_2,ix_3

    DOUBLE PRECISION:: s_period, xc1, xc2, xc3, xxmax, yymax, zzmax
    DOUBLE PRECISION:: xc1Mm, xc2Mm, xc3Mm
    DOUBLE PRECISION:: xx, yy, zz
    DOUBLE PRECISION:: vvx(ixG^T), vvy(ixG^T), vvz(ixG^T)
    DOUBLE PRECISION:: AA
    DOUBLE PRECISION:: delta_z, delta_x, delta_y, exp_x, exp_y, exp_z, exp_xyz, tdep
    !-----------------------------------------------------------------------------

    eqpar(eta_)=0.d0
    eqpar(nu_)=1.0d0
    
    CALL addsource_grav(qdt,ixI^L,ixO^L,iws,qtC,wCT,qt,w)
    
    IF(ABS(eqpar(nu_))>smalldouble)&
         CALL addsource_visc(qdt,ixI^L,ixO^L,iws,qtC,wCT,qt,w)

    !Force Everything as Zero for good measure
    vvx(ixG^T) = 0.d0
    vvy(ixG^T) = 0.d0
    vvz(ixG^T) = 0.d0

    !Define Centre of driver (in Mm)
    xc1Mm=0.1 !Mm        z axis
    xc2Mm=1.0 !Mm        x axis
    xc3Mm=1.0 !Mm        y axis

    !Convert to meters
    xc1 = xc1Mm * 1.0d6  !m        z axis
    xc2 = xc2Mm * 1.0d6  !m        x axis
    xc3 = xc3Mm * 1.0d6  !m        y axis

    !Define size of Domain
    xxmax = 2.d6!MAXVAL(x(1,:,1,2)) !2.d6
    yymax = 2.d6!MAXVAL(x(1,1,:,3)) !2.d6
    zzmax = 1.6d6!MAXVAL(x(:,1,1,1)) !2.d6

    !Define FWHM of gaussian driver profile
    delta_z = 0.05d6
    delta_x = 0.1d6
    delta_y = 0.1d6

    !Set Amplitude of Driver
    !AA=80000.d0 ! ~2km/s hor/vert
    AA = 10.d0

    !Period of Driver
    s_period = 30.d0

    tdep=SIN(qt*2.d0*pi/s_period)!*EXP(-qt/(s_period)) 
 
    DO ix_1 = ixImin1,ixImax1
       DO ix_2 = ixImin2,ixImax2
          DO ix_3 = ixImin3,ixImax3
             
             xx = x(ix_1,ix_2,ix_3,2) - xc2
             yy = x(ix_1,ix_2,ix_3,3) - xc3
             zz = x(ix_1,ix_2,ix_3,1) - xc1  
             
             exp_z = EXP(-zz**2.d0 / (delta_z**2.d0))
             exp_x = EXP(-xx**2.d0 / (delta_x**2.d0))
             exp_y = EXP(-yy**2.d0 / (delta_y**2.d0))
             
             exp_xyz = exp_x * exp_y * exp_z
             
             vvz(ix_1,ix_2,ix_3) = AA * exp_xyz * tdep  
              
          ENDDO
       ENDDO
    ENDDO
    
    DO ix_1 = ixImin1,ixImax1
       DO ix_2 = ixImin2,ixImax2
          DO ix_3 =ixImin3,ixImax3
             
            w(ix_1,ix_2,ix_3,m1_)=w(ix_1,ix_2,ix_3,m1_)+(w(ix_1,ix_2,ix_3,rho_)+w(ix_1,ix_2,ix_3,rhob_))*vvz(ix_1,ix_2,ix_3)*qdt
            
             w(ix_1,ix_2,ix_3,e_) = w(ix_1,ix_2,ix_3,e_) + &
                  (w(ix_1,ix_2,ix_3,rho_) + w(ix_1,ix_2,ix_3,rhob_)) * &
(vvx(ix_1,ix_2,ix_3)**2.d0 + vvy(ix_1,ix_2,ix_3)**2.d0 + vvz(ix_1,ix_2,ix_3)**2.d0) * qdt / 2.d0
          ENDDO
       ENDDO
    ENDDO

    {^IFMPI IF (ipe.EQ.0)} WRITE(*,*) '***time=',qt

  END SUBROUTINE specialsource

!=============================================================================
SUBROUTINE specialbound(qt,ix^L,iw,iB,w)
  INCLUDE 'vacdef.f'

  INTEGER:: ix_1,ix_2

  INTEGER:: iw^LIM,idim^LIM
  DOUBLE PRECISION:: qt,w(ixG^T,1:nw)
  INTEGER:: ix,ix^D,ixe,ixf,ix^L,ixpair^L,idim,iw,iB

  CALL die('not defined')
  RETURN
END SUBROUTINE specialbound

!=============================================================================
SUBROUTINE getdt_special(w,ix^L)

  ! If the Coriolis force is made very strong it may require time step limiting,
  ! but this is not implemented here.

  INCLUDE 'vacdef.f'
  DOUBLE PRECISION:: w(ixG^T,nw)
  INTEGER:: ix^L
  !----------------------------------------------------------------------------

  !call getdt_diff(w,ix^L)

  IF(ABS(eqpar(nu_))>smalldouble)&
       CALL getdt_visc(w,ix^L)

  CALL getdt_grav(w,ix^L)

  RETURN
END SUBROUTINE getdt_special

!=============================================================================
SUBROUTINE specialeta(w,ix^L,idirmin)

  INCLUDE 'vacdef.f'

  DOUBLE PRECISION:: w(ixG^T,nw)
  INTEGER:: ix^L,idirmin
  !---------------------------------------------------------------------------

  CALL die('specialeta is not defined')
END SUBROUTINE specialeta

!=============================================================================
SUBROUTINE savefilelog_special
  INCLUDE 'vacdef.f'

  CALL die("Save file Log Special is not defined")

END SUBROUTINE savefilelog_special
