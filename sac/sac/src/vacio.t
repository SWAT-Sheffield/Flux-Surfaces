!##############################################################################
! module vacio

!=============================================================================
subroutine readparameters(w)

  ! This subroutine sets or reads all default parameters from par/DEFAULT,
  ! then it reads the par/PROBLEM parameter file through standard input,
  ! and the initial data from data/PROBLEM.ini as soon as the filename is read 
  ! from the parameter file.

  include 'vacdef.f'

  double precision:: w(ixG^T,nw)

  character*^LENTYPE :: typepred(nw),typefull(nw),typeimpl(nw),typefilter(nw)
  double precision:: muscleta
  integer:: i,j,k,iw,idim,iB,ifile,isave
  logical:: implmrpc,globalixtest^IFMPI

  ! The use of NAMELIST is not recommended by Fortran 90. It could be replaced
  ! by some string manipulations, but that is difficult in Fortran 77/Adaptor.
  ! I use NAMELIST, since it is simple and convenient.

  namelist /testlist/   teststr,ixtest1,ixtest2,ixtest3,&
       iwtest,idimtest,ipetest^IFMPI
  namelist /filelist/   filenameini,filename,varnames{,^IFMPI npe^D}, &
       typefileini,typefileout,typefilelog,&
       snapshotini,snapshotout,fullgridini,fullgridout,dixB^L
  namelist /savelist/   tsave,itsave,dtsave,ditsave
  namelist /stoplist/   itmax,tmax,tmaxexact,dtmin,residmin,residmax,t,it,&
       cputimemax
  namelist /methodlist/ wnames,fileheadout,eqpar,&
       typeadvance,typefull,typepred,typeimpl,typefilter,&
       typelimiter,typeentropy,entropycoef,artcomp,&
       typelimited,typefct,typetvd,typeaxial,&
       typepoisson,typeconstrain,constraincoef,&
       useprimitive,muscleta,musclomega,&
       acmwidth,acmnolim,acmcoef,acmexpo,fourthorder,&
       implmrpc,&
       dimsplit,typedimsplit,sourcesplit,typesourcesplit,&
       sourceunsplit,&
       divbfix,divbwave,divbconstrain,angmomfix,compactres,&
       smallfix,smallp,smallpcoeff,smallrho,smallrhocoeff,&
       vacuumrho,nproc,procpar
  namelist /boundlist/  nB,typeB,ixB^LIM,idimB,upperB,extraB,typeBscalar,ipairB
  namelist /paramlist/  courantpar,dtpar,dtdiffpar,dtcantgrow,slowsteps,&
       implmrpcpar,&
       implpar,impldiffpar,implerror,implrelax,impldwlimit,&
       implrestart,implrestart2,impliter,impliternr,&
       typeimplinit,typeimpliter,typeimplmat,&
       implnewton,implconserv,implcentered,&
       implnewmat,implpred,impl3level,impljacfast,implsource
  !-----------------------------------------------------------------------------

!!! Set new scalars for sake of unaltered par/DEFAULT
  constraincoef=one
  cputimemax=bigdouble
!!!

  ! Set default values for arrays (except the ones read from the inifile)

  do ifile=1,nfile
     do isave=1,nsavehi
        tsave(isave,ifile)=bigdouble   ! t  of saves into the output files 
        itsave(isave,ifile)=biginteger ! it of saves into the output files 
     end do
     dtsave(ifile)=bigdouble           ! time between saves
     ditsave(ifile)=biginteger         ! timesteps between saves
     isavet(ifile)=1                   ! index for saves by t
     isaveit(ifile)=1                  ! index for saves by it
  end do

  do iw=1,nw
     typepred(iw)='default'     ! Predictor scheme (will be adjusted later)
     typefull(iw)='tvdlf'       ! Full step scheme
     typeimpl(iw)='nul'         ! Implicit step scheme
     typefilter(iw)='nul'       ! Filter scheme
     typelimiter(iw)='minmod'   ! Limiter type for flow variables/characteristics
     typeentropy(iw)='nul'      ! Entropy fix type
     artcomp(iw)=.false.        ! No artificial compression for Harten type TVD
  end do

  do iw=1,nw
     acmcoef(iw)=-one       ! Coefficients (0,1) for the dissipative fluxes
  enddo                     ! negative value means no coefficient is used

  do i=1,nfile+2            ! Elements define processing for the fullstep,
     nproc(i)=0             ! halfstep and the nfile output files. If the value
  end do                    ! is 0, no processing. For nproc(1) and nproc(2)
  ! the value defines the proc. frequency. Negative
  ! value results in a call at every sweep. Positive
  ! value N results in a call at every N-th step before
  ! the first sweep. For nproc(ifile+2) the nonzero 
  ! values cause processing for that file.
  do i=1,nprocpar
     procpar(i)=-one        ! Parameters for processing
  end do
  smallp=-one               ! Default small pressure, redefined in getpthermal
  smallrho=-one             ! Default small density, redefined in keeppositive
  vacuumrho=-one            ! Density representing vacuum

  nB=2*ndim                 ! If nB is not specified by the user, gridsetup 
  do iB=1,nhiB              ! will create 2*ndim boundary regions by default.
     do iw=1,nw             
        typeB(iw,iB) ='cont'  ! Default boundary type
        fixedB(iw,iB)=.false. ! Fixed boundaries are not extrapolated into yet
     end do
     ipairB(iB)=0           ! periodic pair is unknown, but can be set or guessed
  end do
  do iw=1,nw
     do idim=1,ndim
        nofluxB(iw,idim)=.false.  ! No zero flux condition for variables
     enddo
  end do
  dixB^L=2;                 ! Default width of boundary regions
  ixBmax(1,1)=0             ! An impossible value if user specifies boundaries.

  ! Read scalar parameters from the par/DEFAULT file

  unitpar=unitini-1
  open(unitpar,file='par/DEFAULT',status='old')

  read(unitpar,testlist)
  read(unitpar,filelist)
  read(unitpar,savelist)
  read(unitpar,stoplist)
  read(unitpar,methodlist)
  read(unitpar,boundlist)
  read(unitpar,paramlist)

  close(unitpar)
  ! end defaults

  ! Initialize Kronecker delta, and Levi-Civita tensor
  do i=1,3
     do j=1,3
        if(i==j)then
           kr(i,j)=1
       else
           kr(i,j)=0
        endif
        do k=1,3
           if(i==j.or.j==k.or.k==i)then
              lvc(i,j,k)=0
           else if(i+1==j.or.i-2==j)then
              lvc(i,j,k)=1
           else
              lvc(i,j,k)=-1
           endif
        enddo
     enddo
  enddo

  ! Initialize error conunters and equation parameters
  do i=1,nerrcode
     nerror(i)=0
  end do
  do i=1,neqpar+nspecialpar
     eqpar(i)=zero
  end do

  ! read from STDIN
  unitpar=unitstdin
  {^IFMPI
  ! MPI reads from a file
  unitpar=unitini-1
  open(unitpar,file='vac.par',status='old')
  }

  ! Start reading parameters from standard input, i.e. "< par/PROBLEM"
  {ipetest = -1 ^IFMPI}
  read(unitpar,testlist)
  {^IFMPI
  ! ixtest^D is given for the full grid if ipetest was not set explicitly
  globalixtest = ipetest < 0
  ipetest      = max(ipetest,0)
  ! Erase test string for other processors unless ipetest >= npe (test all PEs)
  if(ipetest<npe.and.ipetest/=ipe)teststr=''
  }
  oktest=index(teststr,'readparameters')>=1
  if(oktest) write(unitterm,*)'ReadParameters'
  if(oktest) write(unitterm,testlist)

  varnames='default'
  read(unitpar,filelist)
  {^IFMPI
  ! Extract and check the directional processor numbers and indexes
  ! and concat the PE number to the input and output filenames
  call mpisetnpeDipeD(filenameini)
  call mpisetnpeDipeD(filename(fileout_))
  }

  if (oktest) then 
     {^IFMPI write(unitterm,*)'npe^D=',npe^D}
     write(unitterm,*)filenameini
     do ifile=1,nfile
        write(unitterm,*)filename(ifile)
     enddo
     write(unitterm,*)'Type of ini/out and log files:',&
          typefileini,typefileout,typefilelog
     if(varnames/='default')write(unitterm,*)'Varnames:',varnames
     if(snapshotini>0)write(unitterm,*)'Snapshotini:',snapshotini
     if(snapshotout>0)write(unitterm,*)'Snapshotout:',snapshotout
     write(unitterm,*)'Fullgridini,out:',fullgridini,fullgridout
  endif

  call readfileini(w)

  ! Default for output header line
  fileheadout=fileheadini

  {^IFMPI
  ! Reset global test cell indexes to local ones
  if(globalixtest)call mpiix(ixtest^D,ipetest)
  }

  read(unitpar,savelist)
  do ifile=1,nfile
     if(dtsave(ifile)<bigdouble/2.and.oktest) &
          write(unitterm,'(" DTSAVE  for file",i2," =",g10.5)') &
          ifile,dtsave(ifile)
     if(ditsave(ifile)<biginteger.and.oktest) &
          write(unitterm,'(" DITSAVE for file",i2," =",i10)') &
          ifile,ditsave(ifile)
     if(tsave(1,ifile)==bigdouble.and.itsave(1,ifile)==biginteger.and. &
          dtsave(ifile)==bigdouble.and.ditsave(ifile)==biginteger.and.verbose) &
          write(uniterr,*)'Warning in ReadParameters: ', &
          'No save condition for file ',ifile
  enddo

  read(unitpar,stoplist)
  if(oktest)then
     if(tmax<bigdouble)         write(unitterm,*) 'TMAX= ',tmax
     if(tmaxexact.and.oktest)   write(unitterm,*) 'TMAXEXACT=',tmaxexact
     if(itmax<biginteger)       write(unitterm,*) 'ITMAX=',itmax
     if(dtmin>smalldouble)      write(unitterm,*) 'DTMIN=',dtmin
     if(residmin>smalldouble)   write(unitterm,*) 'RESIDMIN=',residmin
     if(residmax<bigdouble)     write(unitterm,*) 'RESIDMAX=',residmax
     if(cputimemax<bigdouble) write(unitterm,*) 'CPUTIMEMAX=',cputimemax
  endif

  if(tmax==bigdouble.and.itmax==biginteger) write(uniterr,*) &
       'Warning in ReadParameters: Neither tmax nor itmax are given!'

  read(unitpar,methodlist)


  typefull1='nul'
  typepred1='nul'
  typeimpl1='nul'
  typefilter1='nul'
  iw_full(niw_)=0
  iw_impl(niw_)=0
  iw_semi(niw_)=0
  iw_filter(niw_)=0
  do iw=1,nw
     if(typefull(iw)/='nul')then
        typefull1=typefull(iw)
        iw_full(niw_)=iw_full(niw_)+1
        iw_full(iw_full(niw_))=iw
     end if
     if(typepred(iw)/='nul')then
        typepred1=typepred(iw)
     endif
     if(typeimpl(iw)/='nul')then
        typeimpl1=typeimpl(iw)
        iw_impl(niw_)=iw_impl(niw_)+1
        iw_impl(iw_impl(niw_))=iw
     end if
     if(typefull(iw)/='nul'.and.typeimpl(iw)=='nul')then
        iw_semi(niw_)=iw_semi(niw_)+1
        iw_semi(iw_semi(niw_))=iw
     end if
     if(typefilter(iw)/='nul')then
        typefilter1=typefilter(iw)
        iw_filter(niw_)=iw_filter(niw_)+1
        iw_filter(iw_filter(niw_))=iw
     end if
  end do

  if(typefull1=='source'.and.dimsplit.and.sourceunsplit) &
       write(*,*)'Warning: dimensional splitting for source terms!'

  if(typepred1=='default')then
     select case(typefull1)
     case('fct')
        typepred1='fct'
     case('tvdlf','tvdmu')
        typepred1='hancock'
     case('source')
        typepred1='source'
     case('tvdlf1','tvdmu1','tvd','tvd1','tvdmc','cd','cd4','mc','nul')
        typepred1='nul'
     case default
        call die('No default predictor for full step='//typefull1)
     end select
  endif

  muscleta1=(one-muscleta)/4
  muscleta2=(one+muscleta)/4

  polargrid=.false.
  if(typeaxial=='cylinder')then
     if(phi_>ndir)call die(&
          'Error: typeaxial=cylinder with phi>ndir is impossible!')
     if(phi_<=ndim.and.phi_>=2)then
        polargrid=.true.
        write(*,*) 'Using polar coordinates...'
        if(ndir==3.and.(z_<2.or.z_>ndir.or.z_==phi_))call die( &
             'z direction is not set correctly! Use setup -z=? and make vac')
        if(gencoord)write(*,*)'generalized coordinates in the r,z plane...'
     endif
  endif

  if(angmomfix)then
     write(*,*)'Angular momentum conservation is on (angmomfix=T) for iw=',mphi_
     if(typeaxial/='cylinder') call die(&
          'angmomfix works in cylindrical symmetry only!'//&
          ' Set typeaxial in par file')
     if(phi_<2.or.phi_>ndir) call die(&
          'phi direction does not satisfy 1<phi<=ndir! Use setvac -phi!')
     if(gencoord)then
        write(*,*)'angmomfix=T in generalized coordinates...'
        if(polargrid.and.phi_/=ndir) call die(&
             'phi=ndir is required for angmomfix in gen. coordinates!')
     endif
  endif

  ! Artificial compression requires Harten type TVD
  do iw=1,nw
     if(artcomp(iw))typetvd='harten'
  enddo

  ! Harmonize the parameters for dimensional splitting and source splitting
  if(typedimsplit   =='default'.and.     dimsplit)   typedimsplit='xyyx'
  if(typedimsplit   =='default'.and..not.dimsplit)   typedimsplit='unsplit'
  if(typesourcesplit=='default'.and.     sourcesplit)typesourcesplit='sfs'
  if(typesourcesplit=='default'.and..not.sourcesplit)typesourcesplit='unsplit'
  dimsplit   = typedimsplit   /='unsplit'
  sourcesplit= typesourcesplit/='unsplit'
  if(sourcesplit)sourceunsplit=.false.

  {^NOONED
  if(.not.divbconstrain.or.b0_==0)typeconstrain='nul'
  if(typeconstrain=='nul')divbconstrain=.false.

  if(typeconstrain/='nul')then
     divbfix=.false.
     !???divbwave=.false.???
     nproc(1:nfile+2)=0
     if(typeaxial/='slab')write(*,*)&
          'constrainb=T keeps div B = 0 ',&
          'only approximately in cylindrical symm.'
  endif

  if(divbfix.and.(typephys=='mhd'.or.typephys=='mhdiso').and.&
       (.not.sourcesplit).and.(.not.sourceunsplit))&
       call die('divbfix=T requires unsplitsource=T or splitsource=T !')

  if(typefilter1=='tvdmu'.or.typefilter1=='tvdlf')typelimited='predictor'
  if((.not.dimsplit.or.typeimpl1=='tvdlf1'.or.typeimpl1=='tvdmu1')&
       .and.typelimited=='original')typelimited='previous'
  }

  read(unitpar,boundlist)
  if(nB>nhiB) call die('Error in ReadParameters: Too many boundary regions')
  if(ixBmax(1,1)==0.and.oktest)then
     do iB=1,nB
        write(unitterm,'(" TYPEB(",i2,")       = ",100a10)') &
             iB,(typeB(iw,iB),iw=1,nw)
     end do
     write(unitterm,*)' EXTRAB = ',extraB
  end if
  if(ixBmax(1,1)/=0.and.oktest) write(unitterm,boundlist)
  if(dixB^L>dixB^LLIM|.or.|.or.) call die( &
       'Error in ReadParameters: adjust dixBlo,dixBhi')

  ! Initialize implicit parameters as needed
  implsource=sourceunsplit
  if(typeimpl1=='nul')then
     ! Explicit time integration
     implpar=-one
     if(implmrpc)write(*,*)'Warning: implmrpc=T but typeimpl is not set!'
     typeimplinit = 'unused'
     typeimplmat  = 'unused'
  else if(implmrpc)then
     ! MR-PC time integration
     implpar=one
     typeimpliter='vac_mrpc'
     typeimplmat ='free'
     impliter=       1
     implrestart=    5
     implerror=      1.D-3
     if(residmin>0)then
        ! MR-PC for steady state
        impl3level=     .false.
        typeimplinit=   'nul'
     else
        ! MR-PC for time accurate
        impl3level=     .true.
        typeimplinit=   'explicit2'
     endif
  else if(typeimpl1=='source')then
     ! Semi-implicit time integration for sources

     if(.not.sourceunsplit)call die('Implicit sources must be unsplit!')

     if(residmin>zero)then
        write(*,*)'Warning: Semi-implicit sources for steady state!'
        implerror=residmin
     else
        implerror=1.D-5
     endif

     ! Explicit fluxes, implicit sources --> Trapezoidal scheme
     impl3level=.false.
     implpar=half
     typeimplinit='explicit'

     ! Sources only --> Matrix free approach
     typeimplmat='free'

  else
     ! Implicit time integration
     if(residmin>zero)then
        !Steady state calculation --> Backward Euler scheme
        impl3level=.false.
        implpar=one
        typeimplinit='nul'
        implerror=residmin
        if(iw_semi(niw_)>0)write(*,*)&
             'Warning: some variables are explicit for steady state!'
     else
        ! Implicit fluxes (and sources) --> BDF2 scheme
        impl3level=.true.
        ! First step is Backward Euler, but it could be trapezoidal !!!
        implpar=one
        typeimplinit='nul'
        implerror=1.D-3
     endif
     {^IFONED 
     typeimpliter='tridiag'
     typeimplmat='with'
     }
     {^NOONED
     typeimpliter='vac_bicg'
     typeimplmat='prec'
     }
     if((typeimpl1=='tvdlf1'.or.typeimpl1=='cd').and..not.sourceunsplit)then
        impljacfast=.true.
     endif
  endif

  read(unitpar,paramlist)
  {^IFMPI
  close(unitpar)
  }

  if(typeimpl1/='nul')then
     ! Check and/or set parameters for implicit calculations

     if(implpar<=zero)call die('Error: implpar<=0 for implicit calculation!')
     if(implpar>one)call die('Error: implpar>1 for implicit calculation!')

     if(implmrpc.and.typeimpliter/='vac_mrpc')&
          call die('Error: implmrpc=T requires typeimpliter=vac_mrpc')

     if(implpred)then
        implconserv=.false.
        typeadvance='onestep'
        implpar=one
     endif

     if(impl3level)then
        implpred=.false.
        implconserv=.false.
     endif

     if(typeimpliter=='tridiag'.and.ndim==1)typeimplmat='with'
     if(typeimplmat=='prec'.and.ndim==1)typeimpliter='tridiag'

     ! steady state warnings
     if(residmin>zero)then
        if(implpar/=one)write(*,*)'Warning: implpar<1 for steady state!'
        if(impl3level)write(*,*)'Warning: 3 level for steady state!'
     endif

     if(impljacfast.and.typeimpl1/='tvdlf1'.and.typeimpl1/='cd')&
          call die('Error: impljacfast=T works with typeimpl=tvdlf1 or cd!')

     if(typeimpliter=='vac_mrpc')then
        if(implnewton)call die('Error: MRPC with Newton-Raphson iteration!')
        if(residmin>zero.and.typeimplinit/='nul')&
             write(*,*)'Warning: MRPC for steady state',&
             ' should not use typeimplinit=',typeimplinit
     endif

  endif

  ! If all predictor steps are 'nul' then 'twostep' method reduces to 'onestep'
  if(typeadvance=='twostep'.and.typepred1=='nul')typeadvance='onestep'

  if(typefull1=='tvd'.and.typeadvance/='onestep')&
       call die('tvd method should only be used as a onestep method')

  if(typefull1=='tvd' .and. .not.dimsplit)&
       write(*,*)'Warning: One step TVD without dimensional splitting !!!'

  select case(typeadvance)
  case('onestep','adams2')
     nstep=1
  case('twostep')
     nstep=2
  case('threestep')
     nstep=3
  case('fourstep','sterck','jameson')
     nstep=4
  case default
     call die('Unknown typeadvance='//typeadvance)
  end select

  if(oktest) write(unitterm,methodlist)

  if(oktest) write(unitterm,paramlist)

  call setheaderstrings

  return 
end subroutine readparameters

!=============================================================================
subroutine readfileini(w)

  ! Reads from unitini named filenameini in typefilini format.
  !
  ! The file may contain more than one snapshots, in which case the last set is 
  ! read. The compatibility of initial data with internal parameters is checked.

  include 'vacdef.f'

  double precision:: w(ixG^T,nw)

  logical:: fileexist
  character*91:: fhead
  !-----------------------------------------------------------------------------

  if(typefileini=='auto')then
     inquire(FILE=filenameini,EXIST=fileexist)
     if(.not.fileexist) call die('Stop: file does not exist, filenameini='//&
          filenameini)
     open(unitini,FILE=filenameini,STATUS='old')
     read(unitini,'(a91)')fhead
     close(unitini)

     if(ichar(fhead(1:1))/=0.and.ichar(fhead(2:2))/=0.and. &
          ichar(fhead(3:3))/=0.and.ichar(fhead(4:4))/=0)then
        typefileini='ascii'
     else if(ichar(fhead(89:89))==0.and.ichar(fhead(90:90))==0.and. &
          (ichar(fhead(88:88))==24.or.ichar(fhead(91:91))==24))then
        typefileini='binary'
     else
        typefileini='special'
     endif
     if(verbose)then
        write(*,*)'Auto typefileini=',typefileini
        if(typefileini=='special'.and. &
             ichar(fhead(89:89))==0.and.ichar(fhead(90:90))==0.and. &
             (ichar(fhead(88:88))==20.or.ichar(fhead(91:91))==20)) &
             write(*,*)'Looks like a real*4 file'
     endif
  endif
  if(typefileout=='auto')then
     typefileout=typefileini
     if(verbose)write(*,*)'Auto typefileout=',typefileout
  endif

  select case(typefileini)
  case('ascii')
     call readfileini_asc(w)
  case('binary')
     call readfileini_bin(w)
  case default
     call die('Error in VAC: Unknown typefileini='//typefileini)
  end select

  return
end subroutine readfileini

!=============================================================================
subroutine readfileini_asc(w)

  ! Reads from unitini, filenameini in ASCII format.
  !
  ! The file may contain more than one snapshots, in which case the last set is 
  ! read. The compatibility of the initial data with internal parameters checked.
  !
  ! Variables in the order they are read from the file:
  !
  !   fileheadini - a header identifying the input file
  !   it,t        - the initial timestep and time
  !   ndimini     - dimensionality of grid,   Test: ==ndim
  !   neqparini   - number of eq. parameters, Test: <=neqpar+nspecialpar
  !   nwini       - number of flow variables, Test: ==nw
  !   nx          - the grid dimensions,      Test: <=ixGhi-dixBmax-ixMmin+1
  !   eqpar       - equation parameters from filenameini
  !   varnamesini - names of the coordinates, variables, equation parameters
  !                 eg. 'x y rho mx my e bx by  gamma eta'
  !   x           - the (initial) coordinates
  !   w           - the initial flow variables

  include 'vacdef.f'

  double precision:: w(ixG^T,nw)

  logical:: fileexist
  integer:: ios                           ! 0 if not EOF, -1 if EOF, >0 if error
  integer:: ndimini,neqparini,neqparin,nwini,nwin ! values describing input data
  integer:: ix^L,ix^D,idim,iw,ieqpar,snapshot
  double precision:: eqparextra,wextra
  character*^LENNAME :: varnamesini
  !-----------------------------------------------------------------------------

  oktest=index(teststr,'readfileini')>=1

  if(oktest) write(unitterm,*)'ReadFileIni'

  inquire(file=filenameini,exist=fileexist)
  if(.not.fileexist) call die('Stop: file does not exist, filenameini='//&
       filenameini)
  open(unitini,file=filenameini,status='old')

  snapshot=0
  do
     read(unitini,'(a)',iostat=ios,end=100)fileheadini
     if(ios<0)exit                ! Cycle until the last recorded state
     if(oktest) write(unitterm,*)'fileheadini=',fileheadini(1:30)
     read(unitini,*,iostat=ios)it,t,ndimini,neqparini,nwini
     if(oktest) write(unitterm, &
          "('it=',i7,' t=',g10.3,' ndim=',i3,' neqpar=',i3,' nw=',i3)")&
          it,t,ndimini,neqparini,nwini
     gencoord= ndimini<0
     call checkNdimNeqparNw(ndimini,neqparini,nwini,neqparin,nwin)
     read(unitini,*,iostat=ios)nx
     if(oktest) write(unitterm,"('nx =',3i4)")nx
     call setixGixMix(ix^L)
     read(unitini,*,iostat=ios)(eqpar(ieqpar),ieqpar=1,neqparin),&
          (eqparextra,ieqpar=neqparin+1,neqparini)
     if(oktest) write(unitterm,*)eqpar
     read(unitini,'(a)',iostat=ios)varnamesini
     if(varnames=='default')varnames=varnamesini
     if(oktest) write(unitterm,*)varnames

     {do ix^DB=ixmin^DB,ixmax^DB;}
     read(unitini,*,iostat=ios)(x(ix^D,idim),idim=1,ndim),&
          (w(ix^D,iw),iw=1,nwin),(wextra,iw=nwin+1,nwini)
     ^D&enddo;
     if(ios/=0)then
        write(uniterr,*)'Stop: iostat=',ios
        call die('Error in reading file')
     end if
     snapshot=snapshot+1
     if(snapshot==snapshotini)exit
  end do

100 continue

  close(unitini)

  if(oktest) write(*,*)'x,w:',&
       x(ixtest^D,idimtest),w(ixtest^D,iwtest)
  if(oktest) write(*,*)'x,w:',&
       x(ixtest^D,idimtest),w(ixtest^D,1:nw)

  return
end subroutine readfileini_asc

!=============================================================================
subroutine readfileini_bin(w)

  ! Reads from unitini,filenameini in binary format.
  !
  ! The file may contain more than one snapshots, in which case the last set is 
  ! read unless snapshotini is set. 
  ! The compatibility of initial data with internal parameters is checked.

  include 'vacdef.f'

  double precision:: w(ixG^T,nw)

  logical:: fileexist
  integer:: ios                           ! 0 if not EOF, -1 if EOF, >0 if error
  integer:: ndimini,neqparini,neqparin,nwini,nwin ! values describing input data
  integer:: ix^L,idim,iw,ieqpar,snapshot
  double precision:: eqparextra
  character*79 :: varnamesini
  !-----------------------------------------------------------------------------

  oktest=index(teststr,'readfileini')>=1

  if(oktest) write(unitterm,*)'ReadFileIni'

  inquire(file=filenameini,exist=fileexist)
  if(.not.fileexist) call die('Stop: file does not exist, filenameini='//&
       filenameini)
  open(unitini,file=filenameini,status='old',form='unformatted')
 ! PRINT*, ipe, filenameini, unitini
  snapshot=0
  do
     read(unitini,iostat=ios,end=100)fileheadini
     if(ios<0)exit                ! Cycle until the last recorded state
     if(oktest) write(unitterm,*)'fileheadini=',fileheadini(1:30)
     read(unitini,iostat=ios)it,t,ndimini,neqparini,nwini
     if(oktest) write(unitterm, &
          "('it=',i7,' t=',g10.3,' ndim=',i3,' neqpar=',i3,' nw=',i3)")&
          it,t,ndimini,neqparini,nwini
     gencoord= ndimini<0
     call checkNdimNeqparNw(ndimini,neqparini,nwini,neqparin,nwin)
     read(unitini,iostat=ios)nx
     if(oktest) write(unitterm,"('nx =',3i4)")nx
     call setixGixMix(ix^L)
     read(unitini,iostat=ios)(eqpar(ieqpar),ieqpar=1,neqparin),&
          (eqparextra,ieqpar=neqparin+1,neqparini)
     if(oktest) write(unitterm,*)eqpar
     read(unitini,iostat=ios)varnamesini
     if(varnames=='default')varnames=varnamesini
     if(oktest) write(unitterm,*)varnames

     read(unitini,iostat=ios)(x(ix^S,idim),idim=1,ndim)
     ! To conform savefileout_bin we use loop for iw
     do iw=1,nwin
        read(unitini,iostat=ios)w(ix^S,iw)
     end do
     if(ios/=0)then
        write(uniterr,*)'Error in ReadFileIni: iostat=',ios
        call die('Error in reading file')
     end if
     snapshot=snapshot+1
     if(snapshot==snapshotini)exit
  end do

100 continue

  close(unitini)

  if(oktest) write(*,*)'x,w:',&
       x(ixtest^D,idimtest),w(ixtest^D,iwtest)
  if(oktest) write(*,*)'x,w:',&
       x(ixtest^D,idimtest),w(ixtest^D,1:nw)

  return
end subroutine readfileini_bin

!=============================================================================
subroutine checkNdimNeqparNw(ndimini,neqparini,nwini,neqparin,nwin)

  include 'vacdef.f'
  integer:: ndimini,neqparini,nwini,neqparin,nwin
  !-----------------------------------------------------------------------------

  if(ndim/=abs(ndimini))then
     write(*,*)'Error in ReadFileini: ndimini=',ndimini
     call die('Incompatible dimensionalities')
  endif

  if(neqpar+nspecialpar/=neqparini)write(*,"(a,i3,a,i3)")&
       'Warning in ReadFileini: number of eq.params=',neqpar,&
       ' /= neqparini=',neqparini

  if(nw/=nwini)write(*,"(a,i3,a,i3)")&
       'Warning in ReadFileini: number of variables nw=',nw,&
       ' /= nwini=',nwini

  if((neqpar+nspecialpar/=neqparini.or.nw/=nwini).and.varnames=='default')&
       call die('Define varnames (in &filelist for VAC, in 3rd line for VACINI)!')

  ! The number of equation parameters and variables to be read
  neqparin=min(neqparini,neqpar+nspecialpar)
  nwin=min(nwini,nw)

  return
end subroutine checkNdimNeqparNw


!=============================================================================
subroutine setixGixMix(ix^L)

  include 'vacdef.f'
  integer:: ix^L,qnx^IFMPI
  !-----------------------------------------------------------------------------

  ixGmin^D=ixGlo^D;
  ixMmin^D=ixGmin^D+dixBmin^D;

  ! Shave off ghost cells from nx
  if(fullgridini)then
     {^DLOOP
     if(ipe^D==0)^IFMPI       nx(^D)=nx(^D)-dixBmin^D
     if(ipe^D==npe^D-1)^IFMPI nx(^D)=nx(^D)-dixBmax^D
     \}
  endif

  ! Calculate mesh and grid sizes
  ixMmax^D=ixMmin^D+nx(^D)-1;
  ixGmax^D=ixMmax^D+dixBmax^D;

  ! Set the index range for this grid
  if(fullgridini)then
     ix^L=ixG^L;
     {^IFMPI
     ! Set index range to mesh value if the boundary is not an outer bounary
     ^D&if(ipe^D>0)ixmin^D=ixMmin^D\
     ^D&if(ipe^D<npe^D-1)ixmax^D=ixMmax^D\
     }
  else
     ix^L=ixM^L;
  endif

  if(ixGmax^D>ixGhi^D|.or.)then
     write(uniterr,*)'Stop: nxhi=',ixGhi^D-dixBmax^D-ixMmin^D+1
     write(uniterr,*)'Stop: ixGhi=',ixGhi^D
     write(uniterr,*)'Stop: ixGmax=',ixGmax^D
     call die('Error in SetixGixMix')
  end if

  nx^D=nx(^D);

  {^IFMPI
  ! set global grid size by adding up nx and 
  ! dividing by the number of processors in the orthogonal plane
  {^DLOOP
  call MPI_allreduce(nx^D,nxall^D,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierrmpi)
  nxall^D=nxall^D*npe^D/npe;
  \}

  ! Distribute global grid onto processor cube
  call mpigridsetup
  }

  return
end subroutine setixGixMix

!=============================================================================
subroutine setheaderstrings

  ! Check and/or put physics and equation parameter names into file header

  include 'vacdef.f'

  integer:: i
  character*^LENTYPE :: physics
  !-----------------------------------------------------------------------------

  ! Check physics or add _typephysNDIMNDIR

  write(physics,'(a,i1,i1)')typephys,^ND,^NC

  i=index(fileheadout,'_')
  if(i>=1)then
     if(physics/=fileheadout(i+1:i+^LENTYPE))then
        write(*,*)'This code is configured to ',physics
        call die('Error:  physics in file is '//fileheadout(i+1:i+^LENTYPE))
     endif
  else
     i=^LENNAME-len(typephys)-3
     do
        if(fileheadout(i:i)/=' ' .or. i==1)exit
        i=i-1
     enddo
     fileheadout=fileheadout(1:i)//'_'//physics
     write(*,*)'Warning: incomplete input headline.',&
          ' Added to output headline _',physics
  endif

  ! Check for equation parameter names in varnames, add them if missing

  if(varnames/='default' .and. index(varnames,eqparname)<=0)then
     i=^LENNAME-len(eqparname)-3
     do
        if(varnames(i:i)/=' ' .or. i==1)exit
        i=i-1
     enddo
     varnames=varnames(1:i)//'   '//eqparname
  endif

  ! Check for special equation parameter names in varnames, add them if missing

  if(varnames/='default' .and. index(varnames,specialparname)<=0)then
     i=^LENNAME-len(specialparname)-3
     do
        if(varnames(i:i)/=' ' .or. i==1)exit
        i=i-1
     enddo
     varnames=varnames(1:i)//'   '//specialparname
  endif

  return
end subroutine setheaderstrings

!=============================================================================
subroutine savefile(ifile,w)

  include 'vacdef.f'
  integer:: ifile,ix^L
  double precision:: w(ixG^T,nw)
  character*10:: itstring
  !-----------------------------------------------------------------------------

  !if(nproc(ifile+2)>0) call process(ifile+2,1,ndim,w)

  ! In most cases the mesh should be saved
  ix^L=ixM^L;

  select case(ifile)
  case(fileout_)
     ! Produce the output file name
     filenameout=filename(fileout_)
     if(snapshotout>0.and.isaveout>0)then
        if(isaveout==snapshotout*(isaveout/snapshotout))then
           close(unitini+ifile)
           write(itstring,'(i10)')isaveout+1
           filenameout=filenameout(1:index(filenameout,' ')-1)//'_'// &
                itstring(10-int(alog10(isaveout+1.5)):10) 
        endif
     endif
     isaveout=isaveout+1
     if(fullgridout)then
        ix^L=ixG^L;
        {^IFMPI
        ! Set index range to mesh value if not an outer boundary
        ^D&if(ipe^D>0)ixmin^D=ixMmin^D\
        ^D&if(ipe^D<npe^D-1)ixmax^D=ixMmax^D\
        }
     end if
     select case(typefileout)
     case('ascii')
        call savefileout_asc(unitini+ifile,w,ix^L)
     case('binary')
        call savefileout_bin(unitini+ifile,w,ix^L)
     case default
        call die('Error in SaveFile: Unknown typefileout:'//typefileout)
     end select
  case(filelog_)
     select case(typefilelog)
     case('default')
        call savefilelog_default(unitini+ifile,w,ix^L)
     CASE('special')
        CALL savefilelog_special(unitini+ifile,w,ix^L)
     CASE default
        call die('Error in SaveFile: Unknown typefilelog:'//typefilelog)
     end select
  case default
     write(*,*) 'No save method is defined for ifile=',ifile
     call die(' ')
  end select

  return 
end subroutine savefile

!=============================================================================
subroutine savefileout_asc(qunit,w,ix^L)

  ! This version saves into filename(fileout_) ASCII data at every save time
  ! in full accordance with the ReadFileini subroutine, except that the first
  ! line is fileheadout and not fileheadini.

  include 'vacdef.f'

  integer:: qunit,ix^L,ix^D,iw,idim,ndimout
  double precision:: w(ixG^T,nw),qw(nw)
  logical:: fileopen
  !-----------------------------------------------------------------------------

  inquire(qunit,opened=fileopen)

  if(.not.fileopen)open(qunit,file=filenameout,status='unknown')

  if(gencoord)then
     ndimout= -ndim
  else
     ndimout= ndim
  endif

  write(qunit,"(a)")fileheadout
  write(qunit,"(i7,1pe13.5,3i3)")it,t,ndimout,neqpar+nspecialpar,nw
  write(qunit,"(3i4)") ixmax^D-ixmin^D+1
  write(qunit,"(100(1pe13.5))")eqpar
  write(qunit,"(a)")varnames
  {do ix^DB= ixmin^DB,ixmax^DB \}
  ! Values with magnitude less than smalldouble are written as 0d0
  where(abs(w(ix^D,1:nw))>5.0d-16)
     qw(1:nw)=w(ix^D,1:nw)
  elsewhere
     qw(1:nw)=0d0
  end where
  write(qunit,"(100(1pe18.10))")x(ix^D,1:ndim),qw(1:nw)
enddo^D&;

call flushunit(qunit)

return 
end subroutine

!=============================================================================
subroutine savefileout_bin(qunit,w,ix^L)

  ! This version saves into filename(fileout_) binary data at every save time
  ! in full accordance with the ReadFileini subroutine, except that the first
  ! line is fileheadout and not fileheadini.

include 'vacdef.f'

integer:: qunit,ix^L,idim,iw,ndimout
double precision:: w(ixG^T,nw)
logical:: fileopen

!**************** slice
integer:: s_ixmax^D, prom_ixmax^D, s_ixmin^D, prom_ixmin^D  
!**************** endslice


!-----------------------------------------------------------------------------

inquire(qunit,opened=fileopen)
if(.not.fileopen)&
    open(qunit,file=filenameout,status='unknown',form='unformatted')

if(gencoord)then
  ndimout= -ndim
else
  ndimout= ndim
endif

write(qunit)fileheadout
write(qunit)it,t,ndimout,neqpar+nspecialpar,nw
write(qunit) ixmax^D-ixmin^D+1
write(qunit)eqpar
write(qunit)varnames
write(qunit)(x(ix^S,idim),idim=1,ndim)

! write(qunit)w(ix^S,1:nw) produces segmentation fault on Alpha, thus loop

do iw=1,nw
  write(qunit)w(ix^S,iw)
end do

call flushunit(qunit)


!**************** slice *********************************

!inquire(qunit,opened=fileopen)
!if(.not.fileopen)&
!   open(qunit,file=filenameout,status='unknown',form='unformatted')

!if(gencoord)then
!   ndimout= -ndim
!else
!   ndimout= ndim
!endif

!prom_ixmax^D=ixmax^D
!prom_ixmin^D=ixmin^D
!ixmax1=2
!ixmin1=0

!write(qunit)fileheadout
!write(qunit)it,t,ndimout,neqpar+nspecialpar,nw
!write(qunit) ixmax^D-ixmin^D+1
!write(qunit)eqpar
!write(qunit)varnames
!write(qunit)(x(ix^S,idim),idim=1,ndim)

! write(qunit)w(ix^S,1:nw) produces segmentation fault on Alpha, thus loop
!
!do iw=1,nw
!   write(qunit)w(ix^S,iw)
!end do

!call flushunit(qunit)


!************* end slice ********************************** 






return 
end subroutine savefileout_bin

!=============================================================================
subroutine savefilelog_default(qunit,w,ix^L)

  ! This version saves into filename(filelog_) the following formatted data:
  !
  !   fileheadout
  !   STRING_DESCRIBING_COLUMNS
  !   it t dt wmean(1) wmean(2) ... wmean(nw) residual
  !   it t dt wmean(1) wmean(2) ... wmean(nw) residual
  !   it t dt wmean(1) wmean(2) ... wmean(nw) residual
  !   etc.
  !
  ! at every save time. wmean is the volume averaged w, residual is saved 
  ! if residmin>0 is set in the parfile.

include 'vacdef.f'

integer:: qunit,ix^L
double precision:: w(ixG^T,nw)
integer:: iw
logical:: fileopen
double precision:: wmean(nw)
!-----------------------------------------------------------------------------

if(ipe==0)then^IFMPI

  inquire(qunit,opened=fileopen)
  if(.not.fileopen)then
     open(qunit,file=filename(filelog_),status='unknown')
     write(qunit,'(a)')fileheadout
     if(residmin>zero.or.residmax<bigdouble)then
        write(qunit,'(a15,a55,a9)')'it   t   dt    ',wnames,' residual'
     else
        write(qunit,'(a15,a64)')   'it   t   dt    ',wnames
     endif
  endif

endif^IFMPI

do iw=1,nw 
  wmean(iw)=sum(dvolume(ix^S)*w(ix^S,iw))/volume
  {^IFMPI call mpireduce(wmean(iw),MPI_SUM)}
end do

if(ipe==0)then^IFMPI

  if(residmin>zero.or.residmax<bigdouble)then
     write(qunit,'(i7,100(1pe13.5))')it,t,dt,wmean,residual
  else
     write(qunit,'(i7,100(1pe13.5))')it,t,dt,wmean
  endif
  call flushunit(qunit)

endif^IFMPI

return 
end subroutine savefilelog_default

!=============================================================================
subroutine die(message)

character(*) :: message
!-----------------------------------------------------------------------------
{^IFMPI call mpistop(message)}
write(*,*)message

stop
end subroutine die
!=============================================================================
subroutine flushunit(qunit)

  !use F90_UNIX_IO,only : flush ! F90=f95 (NAG)
implicit none

integer::qunit, ierror

!call flush(qunit,ierror) ! F90=f95 (NAG)
!call flush(qunit)   ! OS=Linux, SunOS, UNICOS, T3E, Fujitsu
!call flush_(qunit)  ! OS=AIX, F90=xlf

! no flush on Linux IA64 with Intel compiler

return
end subroutine flushunit
!=============================================================================
! end module vacio
!##############################################################################



