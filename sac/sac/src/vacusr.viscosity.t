!==============================================================================
SUBROUTINE addsource_visc(qdt,ixI^L,ixO^L,iws,qtC,w,qt,wnew)

  ! Add viscosity source to wnew within ixO 

  INCLUDE 'vacdef.f'

  INTEGER::          ixI^L,ixO^L,iws(niw_)
  DOUBLE PRECISION:: qdt,qtC,qt,w(ixG^T,nw),wnew(ixG^T,nw)

  INTEGER:: ix,ix^L,idim,idir,jdir,iiw,iw

  !already declared in vacusr.f
  !double precision:: tmp2(ixG^T)
  DOUBLE PRECISION:: nushk(ixG^T,ndim)



  DOUBLE PRECISION:: tmprhoL(ixG^T), tmprhoR(ixG^T), tmprhoC(ixG^T)
  DOUBLE PRECISION:: tmpVL(ixG^T), tmpVR(ixG^T), tmpVC(ixG^T)
  DOUBLE PRECISION:: tmpBL(ixG^T), tmpBR(ixG^T), tmpBC(ixG^T)

  DOUBLE PRECISION:: tmpL(ixG^T),tmpR(ixG^T), tmpC(ixG^T)

  DOUBLE PRECISION:: nuL(ixG^T),nuR(ixG^T)

  INTEGER:: jx^L,hx^L, hxO^L

  DOUBLE PRECISION:: c_ene,c_shk

  INTEGER:: i,j,k,l,m,ii0,ii1,t00

  DOUBLE PRECISION:: sB

  !-----------------------------------------------------------------------------

  ! Calculating viscosity sources 
  ! involves second derivatives, two extra layers
  CALL ensurebound(2,ixI^L,ixO^L,qtC,w)
  ix^L=ixO^L^LADD1;

  !sehr wichtig
  CALL setnushk(w,ix^L,nushk)

  DO idim=1,ndim
     tmp(ixI^S)=w(ixI^S,rho_)
     CALL setnu(w,rho_,idim,ixO^L,nuR,nuL)      
     CALL gradient1L(tmp,ix^L,idim,tmp2)
     tmpL(ixI^S)=(nuL(ixI^S)+nushk(ixI^S,idim))*tmp2(ixI^S)	     
     CALL gradient1R(tmp,ix^L,idim,tmp2)
     tmpR(ixI^S)=(nuR(ixI^S)+nushk(ixI^S,idim))*tmp2(ixI^S)
     wnew(ixI^S,rho_)=wnew(ixI^S,rho_)+(tmpR(ixI^S)-tmpL(ixI^S))/dx(ixI^S,idim)*qdt
  ENDDO


  DO idim=1,ndim
     tmp(ixI^S)=w(ixI^S,e_)-half*((^C&w(ixI^S,b^C_)**2+)+(^C&w(ixI^S,m^C_)**2+)/(w(ixI^S,rho_)+w(ixI^S,rhob_)))
     CALL setnu(w,173,idim,ixO^L,nuR,nuL)      
     CALL gradient1L(tmp,ix^L,idim,tmp2)
     tmpL(ixI^S)=(nuL(ixI^S)+nushk(ixI^S,idim))*tmp2(ixI^S)      
     CALL gradient1R(tmp,ix^L,idim,tmp2)
     tmpR(ixI^S)=(nuR(ixI^S)+nushk(ixI^S,idim))*tmp2(ixI^S)
     wnew(ixI^S,e_)=wnew(ixI^S,e_)+(tmpR(ixI^S)-tmpL(ixI^S))/dx(ixI^S,idim)*qdt
  ENDDO




  tmprhoC(ixI^S)=w(ixI^S,rho_)+w(ixI^S,rhob_)



  DO k=1,ndim
     jx^L=ix^L+kr(k,^D); 
     hx^L=ix^L-kr(k,^D);
     tmprhoL(ix^S)=((w(ix^S,rho_)+w(ix^S,rhob_))+(w(hx^S,rho_)+w(hx^S,rhob_)))/two
     tmprhoR(ix^S)=((w(jx^S,rho_)+w(jx^S,rhob_))+(w(ix^S,rho_)+w(ix^S,rhob_)))/two

     DO l=1,ndim
	CALL setnu(w,l+m0_,k,ixO^L,nuR,nuL)      
	tmp(ixI^S)=w(ixI^S,m0_+l)/(w(ixI^S,rho_)+w(ixI^S,rhob_))


        DO ii1=0,1
           IF (ii1 .EQ. 0) THEN
              i=k
              ii0=l
           ELSE
              i=l
              ii0=k
           ENDIF



           IF (i .EQ. k) THEN 
              tmpVL(ix^S)=(w(ix^S,m0_+ii0)+w(hx^S,m0_+ii0))/two
              tmpVR(ix^S)=(w(jx^S,m0_+ii0)+w(ix^S,m0_+ii0))/two

              CALL gradient1L(tmp,ix^L,k,tmp2)
              tmpL(ixI^S)=(nuL(ixI^S)+nushk(ixI^S,k))*tmp2(ixI^S)
              CALL gradient1R(tmp,ix^L,k,tmp2)
              tmpR(ixI^S)=(nuR(ixI^S)+nushk(ixI^S,k))*tmp2(ixI^S) 

              tmp2(ixI^S)=(tmprhoR(ixI^S)*tmpR(ixI^S)-tmprhoL(ixI^S)*tmpL(ixI^S))/dx(ixI^S,k)/two

              wnew(ixI^S,m0_+ii0)=wnew(ixI^S,m0_+ii0)+tmp2(ixI^S)*qdt

              tmp2(ixI^S)=(tmpVR(ixI^S)*tmpR(ixI^S)-tmpVL(ixI^S)*tmpL(ixI^S))/dx(ixI^S,k)/two

              wnew(ixI^S,e_)=wnew(ixI^S,e_)+tmp2(ixI^S)*qdt
           ENDIF




           IF (i .NE. k) THEN
              CALL gradient1(tmp,ix^L,k,tmp2)
              tmp2(ixI^S)=tmp2(ixI^S)*(nuL(ixI^S)+nuR(ixI^S)+two*nushk(ixI^S,k))/two/two

              tmp(ixI^S)=tmprhoC(ixI^S)*tmp2(ixI^S)
              CALL gradient1(tmp,ix^L,i,tmpC)

              wnew(ixI^S,m0_+ii0)=wnew(ixI^S,m0_+ii0)+tmpC(ixI^S)*qdt

              tmp(ixI^S)=w(ixI^S,m0_+ii0)*tmp2(ixI^S)
              CALL gradient1(tmp,ix^L,i,tmpC)

              wnew(ixI^S,e_)=wnew(ixI^S,e_)+tmpC(ixI^S)*qdt
           ENDIF

        ENDDO
     ENDDO
  ENDDO





  DO k=1,ndim
     DO l=1,ndim

        IF (k .NE. l) THEN

           CALL setnu(w,b0_+l,k,ixO^L,nuR,nuL)

           DO ii1=0,1

              IF (ii1 .EQ. 0) THEN
                 ii0=k
                 m=l
                 sB=-1.d0
                 j=k
              ENDIF

              IF (ii1 .EQ. 1) THEN 
                 ii0=l    !ii0 is index B
                 m=k      !first derivative
                 sB=1.d0  !sign B
                 j=l      !first B in energy
              ENDIF



              !print*,'k,l,m,j,ii0,ii1=',k,l,m,j,ii0,ii1



              IF (m .EQ. k) THEN

                 jx^L=ix^L+kr(m,^D); 
                 hx^L=ix^L-kr(m,^D);
                 tmpBL(ix^S)=(w(ix^S,b0_+j)+w(hx^S,b0_+j))/two
                 tmpBR(ix^S)=(w(jx^S,b0_+j)+w(ix^S,b0_+j))/two

                 tmp(ixI^S)=w(ixI^S,b0_+l)

                 CALL gradient1L(tmp,ix^L,k,tmp2)
                 tmpL(ixI^S)=(nuL(ixI^S))*tmp2(ixI^S)
                 CALL gradient1R(tmp,ix^L,k,tmp2)
                 tmpR(ixI^S)=(nuR(ixI^S))*tmp2(ixI^S) 

                 wnew(ixI^S,b0_+ii0)=wnew(ixI^S,b0_+ii0)+sB*(tmpR(ixI^S)-tmpL(ixI^S))/dx(ixI^S,k)*qdt

                 wnew(ixI^S,e_)=wnew(ixI^S,e_)+sB*(tmpR(ixI^S)*tmpBR(ixI^S)-tmpL(ixI^S)*tmpBL(ixI^S))/dx(ixI^S,k)*qdt


              ENDIF



              IF (m .NE. k) THEN

                 tmp(ixI^S)=w(ixI^S,b0_+l)

                 CALL gradient1(tmp,ix^L,k,tmp2)

                 tmp2(ixI^S)=tmp2(ixI^S)*(nuL(ixI^S)+nuR(ixI^S))/two

                 CALL gradient1(tmp2,ix^L,m,tmpC)

                 wnew(ixI^S,b0_+ii0)=wnew(ixI^S,b0_+ii0)+sB*tmpC(ixI^S)*qdt

                 tmp2(ixI^S)=tmp2(ixI^S)*w(ixI^S,b0_+j)

                 CALL gradient1(tmp2,ix^L,m,tmpC)

                 wnew(ixI^S,e_)=wnew(ixI^S,e_)+sB*tmpC(ixI^S)*qdt

              ENDIF


           ENDDO
        ENDIF
     ENDDO
  ENDDO




  RETURN
END SUBROUTINE addsource_visc

!=============================================================================
SUBROUTINE setnu(w,iw,idim,ix^L,nuR,nuL)

  ! Set the viscosity coefficient nu within ixO based on w(ixI). 

  INCLUDE 'vacdef.f'

  INTEGER:: ixi^L
  DOUBLE PRECISION:: w(ixG^T,nw)
  DOUBLE PRECISION:: d1R(^SIDEADO),d1L(^SIDEADO)
  DOUBLE PRECISION:: d3R(^SIDEADO),d3L(^SIDEADO)
  DOUBLE PRECISION:: md3R(ixG^T),md3L(ixG^T)
  DOUBLE PRECISION:: md1R(ixG^T),md1L(ixG^T)
  DOUBLE PRECISION:: nuR(ixG^T),nuL(ixG^T)

  DOUBLE PRECISION:: c_tot, c_hyp,cmax(ixG^T), tmp_nu(ixG^T)
  INTEGER:: ix^L,idim, iw
  INTEGER:: kx^L,jx^L,hx^L,gx^L,ixFF^L,jxFF^L,hxFF^L
  INTEGER:: ix_1,ix_2,ix_3

  INTEGER:: ixF^LL,ixF^L,ixY^LL

  LOGICAL:: new_cmax

  DOUBLE PRECISION:: tmp_nuI(^SIDEADD)

  INTEGER:: k,iwc

  INTEGER:: ix,ixe

  {^IFMPI 

  INTEGER :: nmpirequest, mpirequests(2)
  INTEGER :: mpistatus(MPI_STATUS_SIZE,2)
  COMMON /mpirecv/ nmpirequest,mpirequests,mpistatus


  INTEGER:: hpe,jpe

  DOUBLE PRECISION:: tgtbufferR^D(1^D%^LM)
  DOUBLE PRECISION:: tgtbufferL^D(1^D%^LM)
  DOUBLE PRECISION:: srcbufferR^D(1^D%^LM)
  DOUBLE PRECISION:: srcbufferL^D(1^D%^LM)

  INTEGER:: n

  }

  !----------------------------------------------------------------------------

  new_cmax=.TRUE.

  CALL getcmax(new_cmax,w,ix^L,idim,cmax)
  c_tot=MAXVAL(cmax(ix^S))

  {^IFMPI CALL mpiallreduce(c_tot,MPI_MAX)}

  !---------------------------------------------
  ! Set HyperVis coefficients here:
  !---------------------------------------------

  c_hyp=0.4d0 ! 1.4d0 ! 0.6

  IF (iw.EQ.b^D_|.OR.) c_hyp=0.04d0 ! 2d0

  IF (iw .EQ. rho_) c_hyp=0.04d0 !5d0

  IF (iw .EQ. 173) c_hyp=0.04d0 !2d0


  !---------------------------------------------


  IF (iw .NE. 173) THEN     
     tmp_nu(ixG^T)=w(ixG^T,iw)
     IF (iw.EQ.m^D_|.OR.) tmp_nu(ixG^T)=w(ixG^T,iw)/(w(ixG^T,rho_)+w(ixG^T,rhob_))
  ENDIF

  IF (iw .EQ. 173) tmp_nu(ixG^T)=w(ixG^T,e_)-half*((^C&w(ixG^T,b^C_)**2+)+(^C&w(ixG^T,m^C_)**2+)/(w(ixG^T,rho_)+w(ixG^T,rhob_)))


  ixY^LL=ix^L^LADD2;

  ixF^LL=ixY^LL+1;

  tmp_nuI(ixF^T)=tmp_nu(ixY^T)


  {^IFMPI

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)

  n = ^D&(ixFhi^D-ixFlo^D+1)*   

  SELECT CASE(idim)
     { CASE(^D)

     n=n/(ixFhi^D-ixFlo^D+1)

     }
  END SELECT



  SELECT CASE(idim)
     {   CASE(^D)

     IF(npe^D>1)THEN

        nmpirequest =0
        mpirequests(1:2) = MPI_REQUEST_NULL


        !source
        srcbufferL^D(1^D%ixF^T)=tmp_nuI(ixFlo^D+4^D%ixF^T) !left, lower

        srcbufferR^D(1^D%ixF^T)=tmp_nuI(ixFhi^D-4^D%ixF^T) !right, upper

        CALL mpineighbors(^D,hpe,jpe)

        !Patched for intel compiler by Stuart Mumford July 2013 added (1^D%:^) which expands to (1,:,:) etc.
        IF (mpiupperB(^D)) nmpirequest=nmpirequest+1
        IF (mpiupperB(^D)) CALL MPI_IRECV(tgtbufferR^D(1^D%:^),n,MPI_DOUBLE_PRECISION, jpe,10*jpe+0,MPI_COMM_WORLD, mpirequests(nmpirequest),ierrmpi)

        IF (mpilowerB(^D)) nmpirequest=nmpirequest+1
        IF (mpilowerB(^D)) CALL MPI_IRECV(tgtbufferL^D(1^D%:^),n,MPI_DOUBLE_PRECISION, hpe,10*hpe+1,MPI_COMM_WORLD, mpirequests(nmpirequest),ierrmpi)

        CALL MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)

        IF (mpiupperB(^D)) CALL MPI_RSEND(srcbufferR^D(1^D%:^),n,MPI_DOUBLE_PRECISION, jpe,10*ipe+1,MPI_COMM_WORLD,ierrmpi)

        IF (mpilowerB(^D)) CALL MPI_RSEND(srcbufferL^D(1^D%:^),n,MPI_DOUBLE_PRECISION, hpe,10*ipe+0,MPI_COMM_WORLD,ierrmpi)

        CALL MPI_WAITALL(nmpirequest,mpirequests,mpistatus,ierrmpi)

        !target
        tmp_nuI(ixFhi^D+1^D%ixF^T)=tgtbufferR^D(1^D%ixF^T) !right, upper R

        tmp_nuI(ixFlo^D-1^D%ixF^T)=tgtbufferL^D(1^D%ixF^T) !left, lower  L


     ENDIF
     }
  END SELECT

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)
  }


  IF (iw .EQ. 173) THEN 
     iwc=e_ 
  ELSE 
     iwc=iw
  ENDIF

  DO k=0,1  !left-right bc

     IF (typeB(iwc,2*idim-1+k) .NE. 'mpi') THEN
        IF (upperB(2*idim-1+k)) THEN

           SELECT CASE(idim)
              {   CASE(^D)
              tmp_nuI(ixFhi^D+1^D%ixF^T)=tmp_nuI(ixFhi^D-5^D%ixF^T)
              }
           END SELECT

        ELSE

           SELECT CASE(idim)
              {   CASE(^D)
              tmp_nuI(ixFlo^D-1^D%ixF^T)=tmp_nuI(ixFlo^D+5^D%ixF^T)
              }
           END SELECT

        ENDIF
     ENDIF

  ENDDO

  ixF^L=ixF^LL^LSUB1; 

  kx^L=ixF^L+2*kr(idim,^D);  !5:66
  jx^L=ixF^L+kr(idim,^D);  !4:65
  hx^L=ixF^L-kr(idim,^D);  !2:63
  gx^L=ixF^L-2*kr(idim,^D);  !1:62

  ixFF^L=ixF^LL;   !2:65
  jxFF^L=ixF^LL+kr(idim,^D);  !3:66
  hxFF^L=ixF^LL-kr(idim,^D);  !1:64

  d3R(ixF^S)=ABS(3.d0*(tmp_nuI(jx^S)-tmp_nuI(ixF^S))-(tmp_nuI(kx^S)-tmp_nuI(hx^S))) !3:64
  d1R(ixFF^S)=ABS(tmp_nuI(jxFF^S)-tmp_nuI(ixFF^S)) !2:65

  {DO ix_^D=ixmin^D,ixmax^D\}    !3:62  +1=4:63

  md3R(ix_^D)=MAXVAL(d3R(ix_^D+1-kr(idim,^D):ix_^D+1+kr(idim,^D)))
  md1R(ix_^D)=MAXVAL(d1R(ix_^D+1-2*kr(idim,^D):ix_^D+1+2*kr(idim,^D)))

  {ENDDO\}

  WHERE (md1R(ix^S).GT.0.d0)
     nuR(ix^S)=c_tot*c_hyp*md3R(ix^S)/md1R(ix^S)*dx(ix^S,idim)
  ELSEWHERE 
     nuR(ix^S)=0.d0
  END WHERE

  maxviscoef=MAX(MAXVAL(nuR(ix^S)), maxviscoef)


  !************

  d3L(ixF^S)=ABS(3.d0*(tmp_nuI(ixF^S)-tmp_nuI(hx^S))-(tmp_nuI(jx^S)-tmp_nuI(gx^S)))
  d1L(ixFF^S)=ABS(tmp_nuI(ixFF^S)-tmp_nuI(hxFF^S))    

  {DO ix_^D=ixmin^D,ixmax^D\}

  md3L(ix_^D)=MAXVAL(d3L(ix_^D+1-kr(idim,^D):ix_^D+1+kr(idim,^D)))
  md1L(ix_^D)=MAXVAL(d1L(ix_^D+1-2*kr(idim,^D):ix_^D+1+2*kr(idim,^D)))

  {ENDDO\}

  WHERE (md1L(ix^S).GT.0.d0)
     nuL(ix^S)=c_tot*c_hyp*md3L(ix^S)/md1L(ix^S)*dx(ix^S,idim)
  ELSEWHERE 
     nuL(ix^S)=0.d0  
  END WHERE

  maxviscoef=MAX(MAXVAL(nuL(ix^S)), maxviscoef)

  {^IFMPI CALL mpiallreduce(maxviscoef,MPI_MAX)}

  RETURN
END SUBROUTINE setnu


!=============================================================================
!=============================================================================
SUBROUTINE setnushk(w,ix^L,nushk)

  INCLUDE 'vacdef.f'

  !double precision:: w(ixG^T,nw),tmp2(ixG^T),nushk(ixG^T,ndim)
  DOUBLE PRECISION:: w(ixG^T,nw),nushk(ixG^T,ndim)

  DOUBLE PRECISION:: c_shk

  DOUBLE PRECISION:: tmp3(ixG^T)

  INTEGER:: ix^L,idim, iw,i

  INTEGER:: ix_1,ix_2

  DO idim=1,ndim
     nushk(ix^S,idim)=0.d0
  ENDDO


  !--------------------------------------------------
  ! Comment this out and NOT USE A DAMN GOTO!
  !--------------------------------------------------
!!$c_shk=0.5d0
!!$
!!$tmp3(ix^S)=0.d0
!!$
!!$!**************************BEGIN shock viscosity*******************************
!!$      do idim=1,ndim
!!$         tmp(ix^S)=w(ix^S,m0_+idim)/(w(ix^S,rho_)+w(ix^S,rhob_))
!!$         call gradient1(tmp,ix^L,idim,tmp2)
!!$         tmp3(ix^S)=tmp3(ix^S)+tmp2(ix^S)
!!$       enddo
!!$      do idim=1,ndim
!!$        nushk(ix^S,idim)=tmp3(ix^S)*(dx(ix^S,idim)**2.d0)*c_shk
!!$	WHERE (tmp3(ix^S) .ge. 0.d0)
!!$!	  nushk(ix^S,idim)=0.d0
!!$	END WHERE
!!$	nushk(ix^S,idim)=abs(nushk(ix^S,idim))
!!$      enddo
!!$!****************************END shock viscosity*******************************


  RETURN
END SUBROUTINE setnushk



!=============================================================================
SUBROUTINE getdt_visc(w,ix^L)

  ! Check diffusion time limit for dt < dtdiffpar * dx**2 / (nu/rho)

  ! Based on Hirsch volume 2, p.631, eq.23.2.17

  INCLUDE 'vacdef.f'

  DOUBLE PRECISION:: w(ixG^T,nw),dtdiff_visc
  INTEGER:: ix^L,idim, ix_1,ix_2

  INTEGER:: aa

  ! For spatially varying nu you need a common nu array
  DOUBLE PRECISION::tmpdt(ixG^T), nuL(ixG^T),nuR(ixG^T), nushk(ixG^T,ndim)
  COMMON/visc/nuL
  COMMON/visc/nuR
  !-----------------------------------------------------------------------------

  CALL setnushk(w,ix^L,nushk)

  dtdiffpar=0.25d0

  DO idim=1,ndim
     tmpdt(ix^S)=(maxviscoef+nushk(ix^S,idim))       !/(w(ix^S,rho_)+w(ix^S,rhob_))   ! ~1/dt
     dtdiff_visc=dtdiffpar/MAXVAL(tmpdt(ix^S)/(dx(ix^S,idim)**2))
     {^IFMPI CALL mpiallreduce(dtdiff_visc,MPI_MIN)}
     dt=MIN(dt,dtdiff_visc)
  END DO

  maxviscoef=0.d0

  RETURN
END SUBROUTINE getdt_visc


!***** 2-point central finite difference gradient******

SUBROUTINE gradient1(q,ix^L,idim,gradq)
  INCLUDE 'vacdef.f'
  INTEGER:: ix^L,idim
  DOUBLE PRECISION:: q(ixG^T),gradq(ixG^T)
  INTEGER:: hx^L,kx^L
  INTEGER:: minx1^D,maxx1^D,k
  !-----------------------------------------------------------------------------

  hx^L=ix^L-kr(idim,^D);
  kx^L=ix^L+kr(idim,^D);
  gradq(ix^S)=(q(kx^S)-q(hx^S))/dx(ix^S,idim)/two

  minx1^D=ixmin^D+kr(idim,^D);
  maxx1^D=ixmax^D-kr(idim,^D);

  DO k=0,1  !left-right bc
     IF (typeB(1,2*idim-1+k) .NE. 'periodic') THEN
        IF (typeB(1,2*idim-1+k) .NE. 'mpi') THEN
           IF (upperB(2*idim-1+k)) THEN
              SELECT CASE(idim)
                 {   CASE(^D)
                 gradq(ixmax^D^D%ix^S)=0.d0
                 gradq(maxx1^D^D%ix^S)=0.d0
                 }
              END SELECT
           ELSE
              SELECT CASE(idim)
                 {   CASE(^D)
                 gradq(ixmin^D^D%ix^S)=0.d0
                 gradq(minx1^D^D%ix^S)=0.d0
                 }
              END SELECT
           ENDIF
        ENDIF
     ENDIF
  ENDDO


  RETURN
END SUBROUTINE gradient1

!=============================================================================


!*****left upwind forward 2-point non-central finite difference gradient******

SUBROUTINE gradient1L(q,ix^L,idim,gradq)
  INCLUDE 'vacdef.f'
  INTEGER:: ix^L,idim
  DOUBLE PRECISION:: q(ixG^T),gradq(ixG^T)
  INTEGER:: hx^L
  INTEGER:: minx1^D,maxx1^D,k
  !-----------------------------------------------------------------------------

  hx^L=ix^L-kr(idim,^D);
  gradq(ix^S)=(q(ix^S)-q(hx^S))/dx(ix^S,idim)

  minx1^D=ixmin^D+kr(idim,^D);
  maxx1^D=ixmax^D-kr(idim,^D);

  DO k=0,1  !left-right bc
     IF (typeB(1,2*idim-1+k) .NE. 'periodic') THEN
        IF (typeB(1,2*idim-1+k) .NE. 'mpi') THEN
           IF (upperB(2*idim-1+k)) THEN
              SELECT CASE(idim)
                 {   CASE(^D)
                 gradq(ixmax^D^D%ix^S)=0.d0
                 gradq(maxx1^D^D%ix^S)=0.d0
                 }
              END SELECT
           ELSE
              SELECT CASE(idim)
                 {   CASE(^D)
                 gradq(ixmin^D^D%ix^S)=0.d0
                 gradq(minx1^D^D%ix^S)=0.d0
                 }
              END SELECT
           ENDIF
        ENDIF
     ENDIF
  ENDDO


  RETURN
END SUBROUTINE gradient1L

!=============================================================================

!*****right upwind forward 2-point non-central finite difference gradient*****

SUBROUTINE gradient1R(q,ix^L,idim,gradq)
  INCLUDE 'vacdef.f'
  INTEGER:: ix^L,idim
  DOUBLE PRECISION:: q(ixG^T),gradq(ixG^T)
  INTEGER:: hx^L
  INTEGER:: minx1^D,maxx1^D,k
  !-----------------------------------------------------------------------------

  hx^L=ix^L+kr(idim,^D);
  gradq(ix^S)=(q(hx^S)-q(ix^S))/dx(ix^S,idim)

  minx1^D=ixmin^D+kr(idim,^D);
  maxx1^D=ixmax^D-kr(idim,^D);

  DO k=0,1  !left-right bc
     IF (typeB(1,2*idim-1+k) .NE. 'periodic') THEN
        IF (typeB(1,2*idim-1+k) .NE. 'mpi') THEN
           IF (upperB(2*idim-1+k)) THEN
              SELECT CASE(idim)
                 {   CASE(^D)
                 gradq(ixmax^D^D%ix^S)=0.d0
                 gradq(maxx1^D^D%ix^S)=0.d0
                 }
              END SELECT
           ELSE
              SELECT CASE(idim)
                 {   CASE(^D)
                 gradq(ixmin^D^D%ix^S)=0.d0
                 gradq(minx1^D^D%ix^S)=0.d0
                 }
              END SELECT
           ENDIF
        ENDIF
     ENDIF
  ENDDO


  RETURN
END SUBROUTINE gradient1R
