!##############################################################################
! module vac

!=============================================================================
program vac

! Versatile Advection Code, (c) Gabor Toth. Started on Nov 8, 1994. 

include 'vacdef.f'                !declare common variables and parameters

integer:: ifile,ierrcode,iw
double precision:: w(ixG^T,nw),wnrm2,dtold,time0,time1

! functions
logical:: timetofinish,timetosave
double precision:: cputime
!-----------------------------------------------------------------------------
{call mpiinit ^IFMPI}

verbose=.true. .and.ipe==0^IFMPI
if(verbose)then
   write(*,'(a)')'VAC 4.52 configured to'
   write(*,'(a)')'  -d=33 -phi=0 -z=0 -g=68,68,36 -p=mhd -u=Slog'
   write(*,'(a)')'  -on=cd,rk,mpi'
   write(*,'(a)')'  -off=mc,fct,tvdlf,tvd,impl,poisson,ct,gencoord,resist'
   {^IFMPI write(*,'(a,i3,a)')'Running on ',npe,' processors'}
endif

if(ipe==0)^IFMPI  time0=cputime()

call physini                  ! Initialize physics dependent variables
call readparameters(w)        ! Read filenames and parameters for advection
                              ! Read initial data, set ixM,ixG,gencoord

{^NOGEN 
if(gencoord)then
   write(*,*) 'Error: input file contains general grid'
   write(*,*) 'Recompile vac after setvac -on=gencoord is set.'
endif
}
call boundsetup               ! Initialize boundary data
! Initialize grid geometry
if(gencoord)then
   {^IFGEN call gridsetup2}
   {^NOGEN call die('Error: gencoord module is off')}
else
   call gridsetup1
endif

call startup                  ! Initialize it, t, headers etc.

if(verbose)write(*,'(a,f10.3,a)')'Start Advance  ',cputime()-time0,' sec'

call getboundary(t,1,nw,1,ndim,w)
do
   do ifile=1,nfile
      if(timetosave(ifile)) call savefile(ifile,w)
   end do

   ! Determine time step
   if(dtpar>zero)then
      dt=dtpar
   else
      if(courantpar>zero)call getdt_courant(w,ixM^L)
      call getdt(w,ixM^L)
      call getdt_special(w,ixM^L)

      if(dtcantgrow.and.it>itmin)dt=min(dt,dtold)
      dtold=dt
   endif
   if(dtmrpc>zero)dt=min(dt,dtmrpc)

   if (timetofinish(time0)) exit

   ! For slowsteps == 1, use dtpar in the first time step ONLY
   if(slowsteps==1.and.it==itmin)dtpar=-one

   ! For slowsteps > 1, reduce dt for the first few steps 
   if(slowsteps>it-itmin+1) dt=dt*(one-(one-(it-itmin+one)/slowsteps)**2)

   if(tmaxexact)dt=min(dt,tmax-t+smalldouble)

   ! Store w into wold for residual calculations and 
   ! for TVD limiting based on the previous time step.
   wold(ixG^S,1:nw)=w(ixG^S,1:nw)

   ! Advance w (except variables with typefull='nul') by dt in the full grid
   call advance(iw_full,w)

   if(residmin>zero .or. residmax<bigdouble)then
      ! calculate true residual ||w_n+1-w_n|| for steady-state calculations
      residual=zero
      do iw=1,nw
         wnrm2=sum(w(ixM^S,iw)**2)
         {^IFMPI call mpiallreduce(wnrm2,MPI_SUM)}
         if(wnrm2<smalldouble)wnrm2=one
         residual = residual + sum((w(ixM^S,iw)-wold(ixM^S,iw))**2)/wnrm2
      enddo
      {^IFMPI call mpiallreduce(residual,MPI_SUM)}
      residual=sqrt(residual/nw)
   endif  

   it=it+1
   t=t+dt
   !print*, "****t=", t
end do

time1=cputime()-time0


do ifile=1,nfile
   if(itsavelast(ifile)<it)call savefile(ifile,w)
   close(unitini+ifile)
enddo

if(verbose)write(*,'(a,f10.3,a)')'Finish Advance ',time1,' sec'

if(ipe==0)then^IFMPI

if(dt<dtmin)write(unitterm,*)'Warning: dt<dtmin !!!'
if(time1>cputimemax)write(unitterm,*)'Warning: cputimemax exceeded !!!'
do ierrcode=1,nerrcode
   if(nerror(ierrcode)>0)then
      write(*,"(a,i2,a,i5,a)")'Error (code=',&
      ierrcode,') occurred ',nerror(ierrcode),' times !!!'
      select case(ierrcode)
      case(toosmallp_)
         write(*,"(a)")'Error description: Pressure below psmall'
      case(toosmallr_)
         write(*,"(a)")'Error description: Density below rhosmall'
      case(couranterr_)
         write(*,"(a)")'Error description: Courant number above 1'
      case(poissonerr_)
         write(*,"(a)")'Error description: Poisson solver failed'
      end select
   endif
end do
endif^IFMPI

if(verbose)then
   if(implpar>zero)then
      write(*,*)'Number of explicit evaluations:',nexpl
      write(*,*)'Number of Newton iterations   :',nnewton
      write(*,*)'Number of linear iterations   :',niter
      write(*,*)'Number of MatVecs             :',nmatvec
   endif

   if(residmin>zero .or. residmax<bigdouble)then
      write(*,*)'Number of time steps          :',it-itmin
      write(*,*)'Residual and residmin         :',residual,residmin
   endif

   write(*,'(a,f10.3,a)')'Finished VAC   ',cputime()-time0,' sec'
endif

{call mpifinalize ^IFMPI}

end

!=============================================================================
subroutine startup

include 'vacdef.f'
integer:: ifile,iw,ivector,idim,qnvector
!-----------------------------------------------------------------------------

! Initialize dtmrpc which will be calculated by MRPC
dtmrpc=-one

! Initialize dtcourant, which will be calculated by TVD, TVD-MUSCL or TVDLF
do idim=1,ndim
   dtcourant(idim)=bigdouble
enddo

! If dtpar is set, and not only for the first time step (slowsteps/=1)
! then set courantpar<0, so that dtcourant is not calculated at all
if(dtpar>zero.and.slowsteps/=1)courantpar= -one

itmin=it
do ifile=1,nfile
   tsavelast(ifile)=t
   itsavelast(ifile)=it
end do
isaveout=0

! Initialize vectoriw based on iw_vector. vectoriw=-1 for scalar variables,
do iw=1,nw
   vectoriw(iw)=-1
end do
! It points to the 0-th component (iw_vector=m0_,b0_,...) for vector variables.
! Only the first ndim components of the vector variables are rotated in
! generalized coordinates. 
! qnvector is only used to avoid compiler warning when nvector=0

qnvector=nvector
do ivector=1,qnvector
   do idim=1,ndim
      vectoriw(iw_vector(ivector)+idim)=iw_vector(ivector)
   end do
end do

{^IFPHI
   !For 3D polar (r,z,phi) grid m_phi behaves as a scalar
   !and we take advantage of that to conserve angular momentum
   if(angmomfix.and.polargrid.and.gencoord)vectoriw(mphi_)=-1
}

! Initial value for residual and counters
residual=bigdouble
nexpl=0
nnewton=0
nmatvec=0
niter=0

{^IFCT
! fstore should be zero for fluxCD, fluxCT, ryuCT schemes
fstore(ixG^S,1:ndim)=zero
}

return
end

!=============================================================================
subroutine advance(iws,w)

! w(iws,t) -> w(iws,t+qdt) based on typedimsplit and typesourcesplit
!
! Add split sources and fluxes with unsplit sources

include 'vacdef.f'

integer:: iws(niw_)
double precision:: w(ixG^T,nw), w1(ixG^T,nw)
!-----------------------------------------------------------------------------

! Add split sources berforehand if this is required
if(sourcesplit)then
   w1(ixG^S,1:nw)=w(ixG^S,1:nw)
   select case(typesourcesplit)
      case('sf')
         call addsource2(dt,ixG^L,ixM^L,iws,t,w1,t,w)
      case('sfs')
         call addsource2(dt/2,ixG^L,ixM^L,iws,t,w1,t,w)
      case('ssf')
         call addsource2(dt/2,ixG^L,ixG^L,iws,t,w,t,w1)
         call addsource2(dt,ixG^L,ixM^L,iws,t,w1,t,w)
      case('ssfss')
         call addsource2(dt/4,ixG^L,ixG^L,iws,t,w,t,w1)
         call addsource2(dt/2,ixG^L,ixM^L,iws,t,w1,t,w)
      case default
         call die('Error: Unknown typesourcesplit='//typesourcesplit)
   end select
   call getboundary(t,1,nw,1,ndim,w)
endif

! Add fluxes and unsplit sources explicitly or implicitly
if(typeimpl1=='nul')then
   call advance_expl(typefull1,ixG^L,iws,w1,w)
else
   {^IFIMPL
   call advance_impl(ixG^L,w1,w)
   if(.false.)} call die('IMPL module is switched off')
endif

! Add split sources afterwards if this is required
if(sourcesplit)then
   select case(typesourcesplit)
   case('sfs')
      w1(ixG^S,1:nw)=w(ixG^S,1:nw)
      call addsource2(dt/2,ixG^L,ixM^L,iws,t+dt,w1,t+dt,w)
      call getboundary(t+dt,1,nw,1,ndim,w)
   case('ssfss')
      w1(ixG^S,1:nw)=w(ixG^S,1:nw)
      call addsource2(dt/4,ixG^L,ixG^L,iws,t+dt,w ,t+dt,w1)
      call addsource2(dt/2,ixG^L,ixM^L,iws,t+dt,w1,t+dt, w)
      call getboundary(t+dt,1,nw,1,ndim,w)
   end select
endif

return
end

!=============================================================================
subroutine advance_expl(method,ix^L,iws,w1,w)

! w(t) -> w(t+qdt) within ix^L based on typedimsplit, typesourcesplit, nproc
!
! Add fluxes and unsplit sources, possibly with dimensional splitting
! Boundaries should be kept updated by addsource2 and advect
!
! w1 can be ised freely.

include 'vacdef.f'

character*^LENTYPE :: method
integer:: ix^L,iws(niw_)
double precision:: w(ixG^T,nw),w1(ixG^T,nw)

logical:: firstsweep,lastsweep
common/first/firstsweep,lastsweep
!-----------------------------------------------------------------------------

oktest=index(teststr,'advance')>=1
if(oktest)write(*,*)'Advance method,it,w:',method,' ',it,w(ixtest^D,iwtest)

if(ix^L/=ixG^L|.or.|.or.)&
   call die('Error in Advance: No subgrids implemented yet...')

nexpl=nexpl+1
firstsweep=.true.
if(dimsplit)then
   if((it/2)*2.eq.it .or. typedimsplit=='xy')then
      !If typedimsplit='xy', always do the sweeps in order of increasing idim,
      !otherwise for even parity of "it" only, and reverse order for odd. 
      do idimsplit=1,ndim
         lastsweep= idimsplit==ndim
         call advect(method,ix^L,iws,idimsplit,idimsplit,w1,w)
      enddo
   else
      ! If the parity of "it" is odd and typedimsplit=xyyx, do sweeps backwards
      do idimsplit=ndim,1,-1
         lastsweep= idimsplit==1
         call advect(method,ix^L,iws,idimsplit,idimsplit,w1,w)
      enddo
   endif
else
   ! Add fluxes from all directions at once
   lastsweep= .true.
   call advect(method,ix^L,iws,1,ndim,w1,w)
endif

if(typefilter1/='nul')then
   ! We use updated w for the filter fluxes
   w1(ix^S,1:nw)=w(ix^S,1:nw)

   ! Filter according to typefilter1
   select case(typefilter1)
   {^IFTVD
   case('tvd1')
      ! Call tvdlimit with tvd1 (2nd order Lax-Wendroff terms off)
      call tvdlimit(typefilter1,dt,ix^L,ix^L^LSUB2,iw_filter,1,ndim,w1,t+dt,w)
   }
   {^IFTVDLF
   case('tvdlf','tvdmu','tvdlf1','tvdmu1')
      ! Call tvdmusclf with filter method and physical fluxes off
      call tvdmusclf(.false.,typefilter1,dt,ix^L,ix^L^LSUB2,iw_filter,1,ndim,&
                     t+dt,w1,t+dt,w)
   }
   case default
      call die('Error in Advance: typefilter='&
         //typefilter1//' is unknown or module is switched off!')
   endselect
   call getboundary(t+dt,1,nw,1,ndim,w)
endif

call process(0,1,ndim,w)

if(oktest)write(*,*)'Advance new w:',w(ixtest^D,iwtest)

return
end

!=============================================================================
subroutine advect(method,ix^L,iws,idim^LIM,w1,w) 
! Process w if nproc/=0:   		call process
! Add fluxes and unsplit sources in 
! directions idim=idimmin..idimmax:	call advect1
!
! Depending on typeadvance and implpar call advect1 several times

include 'vacdef.f'

character*^LENTYPE :: method
integer:: ix^L,iws(niw_),idim^LIM
double precision:: w1(ixG^T,nw),w(ixG^T,nw)

! For most Runge-Kutta type schemes one more full array is needed
! For classical RK4 another array is needed
{^ANDIFRK 
double precision:: w2(ixG^T,nw),w3(ixG^T,nw)
}

!!!MEMORY Needed for typeadvance='adams2' only
{^IFIMPL{^ANDIFRK save w2}}

logical:: firstsweep,lastsweep
common/first/firstsweep,lastsweep
!-----------------------------------------------------------------------------

oktest=index(teststr,'advect')>=1
if(oktest)write(*,*)'Advect method w:',method,' ',w(ixtest^D,iwtest)

! For negative "nproc(1)" call process, if positive check whether this is the
! first sweep and if "it-itmin" is an integer multiple of "nproc(1)" 
! (the frequency of processing before the whole timestep)
! Processing is done in advance_impl for implicit methods
if(nproc(1)/=0.and.implpar<=zero)then
   if(nproc(1)<0.or.(firstsweep.and.it-itmin==((it-itmin)/nproc(1))*nproc(1)))&
      call process(1,idim^LIM,w)
end if

! Typically use "method" and at least one extra variable w1
w1(ix^S,1:nw)=w(ix^S,1:nw)

istep=0
select case(typeadvance)
case('onestep')
   call advect1(method,dt,ix^L,iws,idim^LIM,t,w1,t,w)
{^IFIMPL{^ANDIFRK
case('adams2')
   ! w=w+dt*R+dt/2*[R-(dt/dtold)*R_n-1]
   ! Use w1=w+dt*R and w2=R_n-1/dtold
   if(it==itmin)then
      call advect1(method,dt,ix^L,iws,idim^LIM,t,w1,t,w)
      w2(ix^S,1:nw)=(w(ix^S,1:nw)-w1(ix^S,1:nw))/dt**2
   else
      call advect1(method,dt,ix^L,iws,idim^LIM,t,w,t,w1)
      w1(ix^S,1:nw)=w1(ix^S,1:nw)-w(ix^S,1:nw)
      w(ix^S,1:nw)=w(ix^S,1:nw)+w1(ix^S,1:nw)&
                   +half*(w1(ix^S,1:nw)-dt**2*w2(ix^S,1:nw))
      w2(ix^S,1:nw)=w1(ix^S,1:nw)/dt**2
   endif
}}
case('twostep')
   ! do predictor step with typepred method to calculate w1 from w, then
   ! full step with method. Fluxes and unsplit sources are taken at w1.
   call advect1(typepred1,dt/2,ix^L,iws,idim^LIM,t     ,w,t,w1)
   call advect1(method   ,dt  ,ix^L,iws,idim^LIM,t+dt/2,w1,t,w)
{^ANDIFRK
case('threestep')
   ! Shu-s third order method based on eq 2.15b of Yee II with signs corrected 
   call advect1(method,dt      ,ix^L,iws,idim^LIM,t     ,w ,t     ,w1)
   w2(ix^S,1:nw)=3*quarter*w(ix^S,1:nw)+quarter*w1(ix^S,1:nw)
   call advect1(method,dt/4    ,ix^L,iws,idim^LIM,t+dt  ,w1,t+dt/4,w2)
   w(ix^S,1:nw)=w(ix^S,1:nw)/3+w2(ix^S,1:nw)*(two/3)
   call advect1(method,dt*two/3,ix^L,iws,idim^LIM,t+dt/2,w2,t+dt/3,w )
}
{^ANDIFRK
case('fourstep')
   ! Classical four step Runge-Kutta
   ! w1=w+Dt/2*k1
   call advect1(method,dt/2    ,ix^L,iws,idim^LIM,t     ,w ,t,w1)
   ! w2=w+Dt/2*k2
   w2(ix^S,1:nw)=w(ix^S,1:nw)
   call advect1(method,dt/2    ,ix^L,iws,idim^LIM,t+dt/2,w1,t,w2)
   ! w3=w+dt*k3
   w3(ix^S,1:nw)=w(ix^S,1:nw)   
   call advect1(method,dt      ,ix^L,iws,idim^LIM,t+dt/2,w2,t,w3)
   ! w1=(w1+2*w2+w3)/3=Dt*(k1+2*k2+2*k3)/6
   w1(ix^S,1:nw)=(w1(ix^S,1:nw)+2*w2(ix^S,1:nw)+w3(ix^S,1:nw)-4*w(ix^S,1:nw))/3
   ! w=w+Dt*k4/6
   call advect1(method,dt/6    ,ix^L,iws,idim^LIM,t+dt  ,w3,t,w)
   ! w=w+Dt*(k1+2*k2+2*k3+k4)/6
   w(ix^S,1:nw)=w(ix^S,1:nw)+w1(ix^S,1:nw)
}
{^ANDIFRK
case('sterck')
   ! H. Sterck has this fourstep time integration, w2 is needed
   call advect1(method,dt*0.12,ix^L,iws,idim^LIM,t        ,w ,t,w1)
   w2(ix^S,1:nw)=w(ix^S,1:nw)
   call advect1(method,dt/4   ,ix^L,iws,idim^LIM,t+dt*0.12,w1,t,w2)
   w1(ix^S,1:nw)=w(ix^S,1:nw)
   call advect1(method,dt/2   ,ix^L,iws,idim^LIM,t+dt/4   ,w2,t,w1)
   call advect1(method,dt     ,ix^L,iws,idim^LIM,t+dt/2   ,w1,t,w )
}
{^ANDIFRK
case('jameson')
   ! Based on eq.2.15 of Yee II
   call advect1(method,dt/4,ix^L,iws,idim^LIM,t     ,w ,t,w1)
   w2(ix^S,1:nw)=w(ix^S,1:nw)
   call advect1(method,dt/3,ix^L,iws,idim^LIM,t+dt/4,w1,t,w2)
   w1(ix^S,1:nw)=w(ix^S,1:nw)
   call advect1(method,dt/2,ix^L,iws,idim^LIM,t+dt/3,w2,t,w1)
   call advect1(method,dt  ,ix^L,iws,idim^LIM,t+dt/2,w1,t,w)
}
case default
   write(*,*)'typeadvance=',typeadvance
   write(*,*)'Error in Advect: Unknown time integration method or RK is off'
   call die('Correct typeadvance or: cd src; setvac -on=rk; make vac')
end select

if(oktest)write(*,*)'Advect final w:',w(ixtest^D,iwtest)

firstsweep=.false.

return
end

!=============================================================================
subroutine advect1(method,qdt,ixI^L,iws,idim^LIM,qtC,wCT,qt,w)

! Process if not first advection and nproc<0 is set
! Advect w to w+qdt*dF(wCT)_idim/dx_idim+qdt*((idimmax-idimmin+1)/ndim)*S(wCT)
! getboundaries

include 'vacdef.f'

character*^LENTYPE :: method
integer:: ixI^L,ixO^L,iws(niw_),idim^LIM,idim
double precision:: qdt,qtC,qt,wCT(ixG^T,nw),w(ixG^T,nw)

logical:: firstsweep,lastsweep
common/first/firstsweep,lastsweep
!-----------------------------------------------------------------------------

istep=istep+1

if(index(teststr,'saveadvect1')>=1) call savefile(fileout_,wCT)

oktest=index(teststr,'advect1')>=1
if(oktest)write(*,*)'Advect1 istep,wCT,w:',&
   istep,wCT(ixtest^D,iwtest),w(ixtest^D,iwtest)

! In the first step wCT=w thus wCT is already processed if there is processing.
! Otherwise for negative "nproc(2)" call process, if positive check whether 
! this is the first sweep and if "it-itmin" is an integer multiple of 
! "nproc(2)" (the frequency of processing before intermediate steps)
! No processing here for implicit methods
if(istep>1.and.nproc(2)/=0.and.implpar<=zero)then
   if(nproc(2)<0.or.(firstsweep.and.it-itmin==((it-itmin)/nproc(2))*nproc(2)))&
      call process(2,idim^LIM,w)
end if

! Shrink ixO^L in all directions by 2
ixO^L=ixI^L^LSUB2;

select case(method)
{^ANDIFCD

case('cd4')
   call centdiff4(qdt,ixI^L,ixO^L,iws,idim^LIM,qtC,wCT,qt,w)
}

case('source')
   if(sourceunsplit)call addsource2(qdt,ixI^L,ixO^L,iws,qtC,wCT,qt,w)
case('nul')
   ! There is nothing to do
   !HPF_ if(.false.)write(*,*)'This avoids an xlhpf compiler bug'
case default
   write(*,*)'Error in Advect1:',method,' is unknown or switched off!'
   call die('Error in Advect1:'//method//' is unknown or switched off!')
end select

call getboundary(qt+qdt,1,nw,1,ndim,w)

if(oktest)write(*,*)'Advect1 final w:',w(ixtest^D,iwtest)

return
end

!=============================================================================
subroutine addsource2(qdt,ixII^L,ixOO^L,iws,qtC,wCT,qt,w)

! Add source within ixOO for iws: w=w+qdt*S[wCT]

include 'vacdef.f'

integer:: ixI^L,ixO^L,ixII^L,ixOO^L,iws(niw_)
double precision:: qdt,qtC,qt,wCT(ixG^T,nw),w(ixG^T,nw)
!-----------------------------------------------------------------------------

oktest=index(teststr,'addsource')>=1
if(oktest)write(*,*)'Add Source qdt,wCT,w:',&
   qdt,wCT(ixtest^D,iwtest),w(ixtest^D,iwtest)

! AddSource and SpecialSource may shrink ixO or expand ixI for derivatives 
ixI^L=ixII^L; ixO^L=ixOO^L;

call specialsource(qdt,ixI^L,ixO^L,iws,qtC,wCT,qt,w)

call     addsource(qdt,ixI^L,ixO^L,iws,qtC,wCT,qt,w)

if(oktest)write(*,*)'wnew:',w(ixtest^D,iwtest)

! If AddSource/SpecialSource shrunk ixO, getboundary is needed.
if(ixO^L^LTixOO^L|.or.|.or.)then
   call getboundary(qt+qdt,1,nw,1,ndim,w)
   if(oktest)write(*,*)'wnew after getboundary:',w(ixtest^D,iwtest)
end if

return
end

!=============================================================================
logical function timetofinish(time0)

! Finish when it or t reached its maximum expected value, or dt is too small,
! or residual is small enough. Other conditions may be included.

include 'vacdef.f'

double precision:: time0, cputime
logical:: okfinish
!-----------------------------------------------------------------------------

okfinish = it>=itmax .or. t>=tmax .or. dt<dtmin .or. &
          (it>itmin.and.(residual<residmin .or. residual>residmax))

if(cputimemax < bigdouble .and. .not.okfinish) &
    okfinish= cputimemax <= cputime()-time0

timetofinish=okfinish

return
end

!=============================================================================
logical function timetosave(ifile)

! Save times are defined by either tsave(isavet(ifile),ifile) or 
! itsave(isaveit(ifile),ifile) or dtsave(ifile) or ditsave(ifile)
! Other conditions may be included.

include 'vacdef.f'

integer:: ifile
logical:: oksave
!-----------------------------------------------------------------------------

oksave=.false.
if(t>=tsave(isavet(ifile),ifile))then
   oksave=.true.
   isavet(ifile)=isavet(ifile)+1
end if
if(it==itsave(isaveit(ifile),ifile))then
   oksave=.true.
   isaveit(ifile)=isaveit(ifile)+1
end if
if(it==itsavelast(ifile)+ditsave(ifile)) oksave=.true.
if(t >=tsavelast(ifile) +dtsave(ifile) ) oksave=.true.
if(oksave)then
   tsavelast(ifile) =t
   itsavelast(ifile)=it
end if
timetosave=oksave

return
end

!=============================================================================
subroutine getdt_courant(w,ix^L)

! Ensure that the courant conditions is met
! Calculate the time for the  maximum propagation speed cmax_i to cross dx_i
! in each i directions then take minimum for all grid points in the mesh and 
! for all i directions, finally multiply by courantpar.
!
! If TVD or TDLF provides dtcourant(idim) we take the minimum of those.
! In case of generalized coordinates dtcourant(idim) is correct due to the
! rotations while the value calculated here does not use a rotation.

include 'vacdef.f'

double precision:: w(ixG^T,nw),cmax(ixG^T),courantmax,dtold
integer:: ix^L,idim
logical:: new_cmax
!-----------------------------------------------------------------------------

oktest=index(teststr,'getdt')>=1

if(oktest) write(*,*)'getdt_courant'

dtold=dt
dt=bigdouble
courantmax=zero
new_cmax=.true.
do idim=1,ndim
   if(dtcourant(idim)<bigdouble)then
      ! If dtcourant(idim) is calculated, use it
!!!      if(it==itmin+1)write(*,*)'second order correction in dt_courant!!!'
!!!      dt=min(dt,dtcourant(idim),dtcourant(idim)**2/dtold,1.1*dtold)
      dt=min(dt,dtcourant(idim))
      if(oktest) write(*,*)'idim,dtcourant(idim)',idim,dtcourant(idim)
   else
      ! dx>0, but cmax>=0 may actually be 0, thus we calculate 
      ! max(cmax/dx) rather than min(dx/cmax).

      call getcmax(new_cmax,w,ix^L,idim,cmax)
      courantmax=max(courantmax,maxval(cmax(ix^S)/dx(ix^S,idim)))

      if(gencoord.and.it==itmin+1.and.verbose)write(*,*)&
          'Warning in GetDtCourant: for gencoord approx. only',&
          ', better use TVD-type methods'
      if(oktest) write(*,*)'idim,cmax:',idim,cmax(ixtest^D)
      if(oktest) write(*,*)'max(c/dx)',maxval(cmax(ix^S)/dx(ix^S,idim))
   endif
end do
{^IFMPI call mpiallreduce(courantmax,MPI_MAX)}
if(index(teststr,'dtdecline')<1)then
   do idim=1,ndim
      dtcourant(idim)=bigdouble
   enddo
endif
if(courantmax>smalldouble) dt=min(dt,courantpar/courantmax)

if(oktest) write(*,*)'GetDtCourant dt=',dt

return 
end

!=============================================================================
double precision function cputime()

! Return cputime in seconds as a double precision number.
! For g77 compiler replace F77_ with F77_ everywhere in this function
! so that f90tof77 does not touch the system_clock function.

integer:: clock,clockrate,count_max !HPF_ !F77_
!F77_ real:: etime,total,tarray(2)
!F77_ external etime
!HPF_ real:: timef
!-----------------------------------------------------------------------------

cputime=-1.D0                             ! No timing

call system_clock(clock,clockrate,count_max) !HPF_ !F77_
cputime=clock*(1.D0/clockrate)               !HPF_ !F77_
!F77_ total = etime(tarray)
!F77_ cputime=tarray(1)
!HPF_ cputime=timef()/1.0D3

!cputime=second()   ! Cray CF77 or F90 (total CPU time for more CPU-s)
!cputime=secondr()  ! Cray CF77 or F90 (elapsed time for more CPU-s)

return
end
!=============================================================================
! end module vac
!##############################################################################
