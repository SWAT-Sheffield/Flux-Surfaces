!------------------------------------------------------------------------------
SUBROUTINE savefilelog_special(qunit,w,ix^L)
  
  ! This is a save log file routine to calculate and save out Vpar Vperp and Vaz
  ! It mimics savefileout_bin to mantain compatibility with usual readin routines
  
  INCLUDE 'vacdef.f'

  INTEGER:: qunit,ix^L,idim,iw,ndimout
  DOUBLE PRECISION :: w(ixG^T,nw),vpar(ixG^T),vperp(ixG^T),vaz(ixG^T)
  DOUBLE PRECISION :: Bz(ixG^T), By(ixG^T), Bx(ixG^T), Vx(ixG^T), Vy(ixG^T), Vz(ixG^T),DelBz(ixG^T), DelBy(ixG^T), DelBx(ixG^T),rhot(ixG^T),DelBxB(ixG^T,3),B(ixG^T)
  DOUBLE PRECISION :: Velout(ixG^T,ndim)
  LOGICAL:: fileopen
  
  INTEGER :: ix_1,ix_2,ix_3
  CHARACTER*10:: itstring
  CHARACTER*3 :: x1
  CALL die("this dont work")
  !===========================================================================
  !Calculate the specials to save out
  !===========================================================================
  !Total density
  rhot = w(ixG^T,rhob_) + w(ixG^T,rho_)

  !Make velocity Components
  vz = w(ixG^T,m1_) / rhot(ixG^T)
  vx = w(ixG^T,m2_) / rhot(ixG^T)
  vy = w(ixG^T,m3_) / rhot(ixG^T)
  
  !Total Magnetic Field
  Bz = w(ixG^T,b1_) + w(ixG^T,bg1_)
  Bx = w(ixG^T,b2_) + w(ixG^T,bg2_)
  By = w(ixG^T,b3_) + w(ixG^T,bg3_)

  B = SQRT(Bz**2 + Bx**2 + By**2)
  
  !Calculate DelB
  DelBz = B(ixG^T) / dx(ixG^T,1)
  DelBx = B(ixG^T) / dx(ixG^T,2)
  DelBy = B(ixG^T) / dx(ixG^T,3)



  DO ix_1=ixGlo1,ixGhi1
     DO ix_2=ixGlo2,ixGhi2
        DO ix_3=ixGlo3,ixGhi3

           DelBxB(ix_1,ix_2,ix_3,:) = cross(DelBz(ix_1,ix_2,ix_3),DelBx(ix_1,ix_2,ix_3),DelBy(ix_1,ix_2,ix_3),Bz(ix_1,ix_2,ix_3),Bx(ix_1,ix_2,ix_3),By(ix_1,ix_2,ix_3))
           
           !Parallel velocity is V.B
           Vpar(ix_1,ix_2,ix_3) = Bz(ix_1,ix_2,ix_3) * Vz(ix_1,ix_2,ix_3) + Bx(ix_1,ix_2,ix_3) * Vx(ix_1,ix_2,ix_3) + By(ix_1,ix_2,ix_3) * Vy(ix_1,ix_2,ix_3)
           
           !Perp Velocity is V.DelB
           Vperp(ix_1,ix_2,ix_3) = DelBz(ix_1,ix_2,ix_3) * Vz(ix_1,ix_2,ix_3) + DelBx(ix_1,ix_2,ix_3) * Vx(ix_1,ix_2,ix_3) + DelBy(ix_1,ix_2,ix_3) * Vy(ix_1,ix_2,ix_3)
           
           !Azimuthal Velocity is DelBxB.V
           Vperp(ix_1,ix_2,ix_3) = DelBxB(ix_1,ix_2,ix_3,1) * Vz(ix_1,ix_2,ix_3) + DelBxB(ix_1,ix_2,ix_3,2) * Vx(ix_1,ix_2,ix_3) + DelBxB(ix_1,ix_2,ix_3,3) * Vy(ix_1,ix_2,ix_3)

        END DO
     END DO
  END DO

  VelOut(ixG^T,1) = Vpar
  VelOut(ixG^T,2) = Vperp
  VelOut(ixG^T,2) = Vaz
  !==========================================================================
  filenameout = filename(filelog_)
  PRINT*, filenameout
  WRITE (x1,'I3.3') ipe
  filenameout = filenameout(1:INDEX(filenameout,'.')-1)//'_np020204_'//TRIM(x1)//'.log'
  PRINT*, filenameout

  INQUIRE(qunit,opened=fileopen)
  IF(.NOT.fileopen)&
       OPEN(qunit,file=filenameout,status='unknown',form='unformatted')
  
  IF(gencoord)THEN
     ndimout= -ndim
  ELSE
     ndimout= ndim
  ENDIF
  
  WRITE(qunit)fileheadout
  WRITE(qunit)it,t,ndimout,neqpar+nspecialpar,nw
  WRITE(qunit) ixmax^D-ixmin^D+1
  WRITE(qunit)eqpar
  WRITE(qunit) "Vpar","Vperp","Vaz"
  WRITE(qunit)(dx(ix^S,idim),idim=1,ndim)
  

  

  DO iw=1,ndim
     WRITE(qunit)VelOut(ix^S,iw)
  END DO
  


  CALL flushunit(qunit)

CONTAINS
  FUNCTION cross(a1,a2,a3,b1,b2,b3)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: a1,a2,a3,b1,b2,b3
    DOUBLE PRECISION :: cross(3)

    cross(1) = a2*b3 - a3*b2
    cross(2) = a1*b3 - a3*b1
    cross(3) = a1*b2 - a2*b1
  END FUNCTION cross

END SUBROUTINE savefilelog_special
