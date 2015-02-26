!##############################################################################
! include vacpar - mhd

! For MHD: density,momentum,energy,magnetic_field + adiabatic_index,resistivity

CHARACTER*3,PARAMETER:: typephys='mhd'           ! VACPHYS module name
CHARACTER*9,PARAMETER:: eqparname='gamma eta'    ! Equation parameter names

INTEGER,PARAMETER:: rho_=1,m0_=rho_,m^C_=m0_+^C,e_=m^NC_+1! flow variables
INTEGER,PARAMETER:: ee_=e_, b0_=e_,b^C_=b0_+^C            ! flow variables
INTEGER,PARAMETER:: eb_=b^NC_+1, rhob_=eb_+1
INTEGER,PARAMETER:: bg0_=rhob_, bg^C_=bg0_+^C
INTEGER,PARAMETER:: nw=bg^NC_                              ! No. flow variables

INTEGER,PARAMETER:: v0_=m0_, v^C_=m^C_, p_=e_, pp_=ee_    ! Primitive variables

INTEGER,PARAMETER:: fastRW_=1,fastLW_=2,slowRW_=3,slowLW_=4 ! Characteristic
INTEGER,PARAMETER:: entroW_=5,diverW_=6,alfvRW_=7,alfvLW_=8 ! waves
INTEGER,PARAMETER:: extremeRW_=fastRW_,extremeLW_=fastLW_   ! Potential extrema

INTEGER,PARAMETER:: mr_=m0_+r_,mphi_=m0_+phi_,mz_=m0_+z_  ! Polar var. names
INTEGER,PARAMETER:: br_=b0_+r_,bphi_=b0_+phi_,bz_=b0_+z_

INTEGER,PARAMETER:: nvector=2                             ! No. vector vars

INTEGER,PARAMETER:: gamma_=1,eta_=2,neqpar=2              ! equation params
INTEGER,PARAMETER:: divbcoeff_=1,divbconst_=2,divbbound_=3,nprocpar=3
                                                          ! processing params

! end include vacpar - mhd
!##############################################################################
