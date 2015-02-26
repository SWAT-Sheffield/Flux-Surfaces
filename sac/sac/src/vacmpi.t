!=============================================================================
subroutine mpiinit

! Initialize MPI variables

include 'vacdef.f'
!----------------------------------------------------------------------------
call MPI_INIT(ierrmpi)
call MPI_COMM_RANK (MPI_COMM_WORLD, ipe, ierrmpi)
call MPI_COMM_SIZE (MPI_COMM_WORLD, npe, ierrmpi)

! unset values for directional processor numbers
npe^D=-1;
! default value for test processor
ipetest=0

return
end

!==============================================================================
subroutine mpifinalize

include 'vacdef.f'

call MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)
call MPI_FINALIZE(ierrmpi)

return
end

!==============================================================================
subroutine ipe2ipeD(qipe,qipe^D)

! Convert serial processor index to directional processor indexes

include 'vacdef.f'

integer:: qipe^D, qipe
!-----------------------------------------------------------------------------
qipe1 = qipe - npe1*(qipe/npe1)
{qipe2 = qipe/npe1 - npe2*(qipe/(npe1*npe2)) ^NOONED}
{qipe3 = qipe/(npe1*npe2)                    ^IFTHREED}

return
end

!==============================================================================
subroutine ipeD2ipe(qipe^D,qipe)

! Convert directional processor indexes to serial processor index

include 'vacdef.f'

integer:: qipe^D, qipe
!-----------------------------------------------------------------------------
qipe = qipe1 {^NOONED + npe1*qipe2} {^IFTHREED + npe1*npe2*qipe3}

return
end

!==============================================================================
subroutine mpisetnpeDipeD(name)

! Set directional processor numbers and indexes based on a filename.
! The filename contains _np followed by np^D written with 2 digit integers.
! For example _np0203 means np1=2, np2=3 for 2D.

include 'vacdef.f'
character*^LENNAME :: name, nametail
integer:: i,qnpe^D
logical:: npeDknown,npeDinname
!-----------------------------------------------------------------------------

oktest = index(teststr,'mpisetnpeDipeD')>0
if(oktest)write(*,*)'mpisetnpeDipeD ipe,name=',ipe,name

! Check if npe^D is already known
npeDknown  = npe1>0
if(npedknown .and. npe^D* /= npe)then
   write(*,*)'npe=',npe,' /= product of npe^D=',npe^D
   call mpistop('ERROR in setnpeDipeD')
endif

! Check if npe^D is given in the name
i=index(name,'_np')+3
npeDinname = i>3

if(.not.(npeDknown.or.npeDinname))call mpistop( &
   'ERROR in setnpeDipeD: npeD is neither known nor given in name='//name)

if(npeDinname)then
   ! read npe^D from name
   read(name(i:i+5),'(3i2)') qnpe^D
   i=i+2*^ND
   nametail=name(i:^LENNAME)
endif

if( npeDknown .and. npeDinname )then
   ! Check agreement
   if( qnpe^D/=npe^D|.or. )then
      write(*,*)'npe^D=',npe^D,' /= qnpe^D=',qnpe^D,' read from filename=',name
      call mpistop('ERROR in mpisetnpeDipeD')
   endif
endif

if(npeDinname .and. .not.npeDknown)then
   ! set npe^D based on name
   npe^D=qnpe^D;
   if( npe^D* /= npe)then
      write(*,*)'npe=',npe,' /= product of npe^D=',npe^D,&
         ' read from filename=',name
      call mpistop('ERROR in setnpeDipeD')
   endif
endif

! Get directional processor indexes
call ipe2ipeD(ipe,ipe^D)

if(npeDknown .and. .not.npeDinname)then
   ! insert npe^D into name
   i=index(name,'.')
   nametail=name(i:^LENNAME)
   write(name(i:^LENNAME),"('_np',3i2.2)") npe^D
   i = i+3+2*^ND
endif

! insert ipe number into the filename
write(name(i:^LENNAME),"('_',i3.3,a)") ipe,nametail(1:^LENNAME-i-4)

! Set logicals about MPI boundaries for this processor
{^DLOOP
mpiupperB(^D)=ipe^D<npe^D-1
mpilowerB(^D)=ipe^D>0 \}

if(oktest)write(*,*)'mpisetnpeDipeD: ipe,npeD,ipeD,name=',ipe,npe^D,ipe^D,name

return
end

!==============================================================================
subroutine mpineighbors(idir,hpe,jpe)

! Find the hpe and jpe processors on the left and right side of this processor 
! in direction idir. The processor cube is taken to be periodic in every
! direction.

include 'vacdef.f'

integer :: idir,hpe,jpe,hpe^D,jpe^D
!-----------------------------------------------------------------------------
hpe^D=ipe^D-kr(^D,idir);
jpe^D=ipe^D+kr(^D,idir);
{^DLOOP
if(hpe^D<0)hpe^D=npe^D-1
if(jpe^D>=npe^D)jpe^D=0\}

call ipeD2ipe(hpe^D,hpe)
call ipeD2ipe(jpe^D,jpe)

return
end
!==============================================================================
subroutine mpigridsetup

! Distribute a grid of size nxall^D onto PE-s arranged in a cube of size npe^D

include 'vacdef.f'
!-----------------------------------------------------------------------------
!!!write(*,*)'nxall,npe=',nxall^D,npe^D

! Grid sizes on the processors (except for the last ones in some direction)
! This formula optimizes the load balance when nx^D is not a multiple of npe^D
nxpe^D=(nxall^D-1)/npe^D+1; 

! Global grid indexes of the first grid point stored on this PE
ixPEmin^D=ipe^D*nxpe^D+1;

! The last processors in a direction may have smaller grid sizes than nxpe
{^DLOOP 
if(ipe^D < npe^D-1)then
   nx^D = nxpe^D
else
   nx^D = nxall^D - ixpemin^D + 1
endif
\}

! Global grid index of the last grid point stored on this PE
ixPEmax^D=ixPEmin^D+nx^D-1;

return
end

!=============================================================================
subroutine mpireduce(a,mpifunc)

! reduce input for one PE 0 using mpifunc

include 'mpif.h'

double precision :: a, alocal
integer          :: mpifunc, ierrmpi
!----------------------------------------------------------------------------
alocal = a
call MPI_REDUCE(alocal,a,1,MPI_DOUBLE_PRECISION,mpifunc,&
                0,MPI_COMM_WORLD,ierrmpi)

return
end

!==============================================================================
subroutine mpiallreduce(a,mpifunc)

! reduce input onto all PE-s using mpifunc

include 'mpif.h'

double precision :: a, alocal
integer          :: mpifunc, ierrmpi
!-----------------------------------------------------------------------------
alocal = a
call MPI_ALLREDUCE(alocal,a,1,MPI_DOUBLE_PRECISION,mpifunc,&
                   MPI_COMM_WORLD,ierrmpi)

return
end

!==============================================================================
subroutine mpiix(ix^D,jpe)

! Convert ix^D physical cell index on the global grid to local indexes 
! and set the processor number jpe to the processor that contains the cell

include 'vacdef.f'
integer :: ix^D, jpe, jpe^D
!-----------------------------------------------------------------------------

! Directional processor indexes
jpe^D=(ix^D-ixMmin^D)/nxpe^D;

! Conversion to local index
ix^D=ix^D-jpe^D*nxpe^D;

! Get MPI processor index
call ipeD2ipe(jpe^D,jpe)

return
end

!==============================================================================
subroutine mpiixlimits(ix^L)

! Convert global index limits to local index limits for this PE

include 'vacdef.f'
integer :: ix^L
!-----------------------------------------------------------------------------
{^DLOOP
if(ixmin^D > ixPEmax^D)then
   ixmin^D = nx^D
   ixmax^D = nx^D-1
elseif(ixmax^D < ixPEmin^D)then
   ixmax^D = 0
   ixmin^D = 1
else
   ixmin^D = max(ixmin^D,ixPEmin^D) - ixPEmin^D + 1
   ixmax^D = min(ixmax^D,ixPEmax^D) - ixPEmin^D + 1
endif
\}

return
end
!==============================================================================

subroutine mpistop(message)

! Stop MPI run in an orderly fashion

include 'vacdef.f'

character(*) :: message
integer :: nerrmpi

!------------------------------------------------------------------------------
write(*,*)'ERROR for processor',ipe,':'
write(*,*)message
call MPI_abort(MPI_COMM_WORLD, nerrmpi, ierrmpi)

stop
end

!==============================================================================
subroutine mpibound(nvar,var)

! Fill in ghost cells of var(ixG,nvar) from other processors

include 'vacdef.f'

integer :: nvar
double precision :: var(ixG^T,nvar)

! processor indexes for left and right neighbors
integer :: hpe,jpe
! index limits for the left and right side mesh and ghost cells 
integer :: ixLM^L, ixRM^L, ixLG^L, ixRG^L
logical :: periodic

! There can be at most 2 receives in any direction for each PE
integer :: nmpirequest, mpirequests(2)
integer :: mpistatus(MPI_STATUS_SIZE,2)
common /mpirecv/ nmpirequest,mpirequests,mpistatus
!-----------------------------------------------------------------------------
oktest=index(teststr,'mpibound')>0
if(oktest)write(*,*)'mpibound ipe,nvar,varold=',&
    ipe,nvar,var(ixtest^D,min(nvar,iwtest))

{^DLOOP
if(npe^D>1)then
   nmpirequest =0
   mpirequests(1:2) = MPI_REQUEST_NULL

   periodic=typeB(1,2*^D)=='mpiperiod'

   ! Left and right side ghost cell regions (target)
   ixLG^L=ixG^L; ixLGmax^D=ixMmin^D-1;
   ixRG^L=ixG^L; ixRGmin^D=ixMmax^D+1;

   ! Left and right side mesh cell regions (source)
   ixLM^L=ixG^L; ixLMmin^D=ixMmin^D; ixLMmax^D=ixMmin^D+dixBmin^D-1;
   ixRM^L=ixG^L; ixRMmax^D=ixMmax^D; ixRMmin^D=ixMmax^D-dixBmax^D+1;

   ! Obtain left and right neighbor processors for this direction
   call mpineighbors(^D,hpe,jpe)

   if(oktest)then
      write(*,*)'mpibound ipe,idir=',ipe,^D
      write(*,*)'mpibound ipe,ixLG=',ipe,ixLG^L
      write(*,*)'mpibound ipe,ixRG=',ipe,ixRG^L
      write(*,*)'mpibound ipe,ixLM=',ipe,ixLM^L
      write(*,*)'mpibound ipe,ixRM=',ipe,ixRM^L
      write(*,*)'mpibound ipe,hpe,jpe=',ipe,hpe,jpe
   endif

   ! receive right (2) boundary from left neighbor hpe
   if(mpilowerB(^D) .or. periodic)call mpirecvbuffer(nvar,ixRM^L,hpe,2)
   ! receive left (1) boundary from right neighbor jpe
   if(mpiupperB(^D) .or. periodic)call mpirecvbuffer(nvar,ixLM^L,jpe,1)
   ! Wait for all receives to be posted
   call MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)
   ! Ready send left (1) boundary to left neighbor hpe
   if(mpilowerB(^D) .or. periodic)call mpisend(nvar,var,ixLM^L,hpe,1)
   ! Ready send right (2) boundary to right neighbor
   if(mpiupperB(^D) .or. periodic)call mpisend(nvar,var,ixRM^L,jpe,2)
   ! Wait for messages to arrive
   call MPI_WAITALL(nmpirequest,mpirequests,mpistatus,ierrmpi)
   ! Copy buffer received from right (2) physical cells into left ghost cells
   if(mpilowerB(^D) .or. periodic)call mpibuffer2var(2,nvar,var,ixLG^L)
   ! Copy buffer received from left (1) physical cells into right ghost cells
   if(mpiupperB(^D) .or. periodic)call mpibuffer2var(1,nvar,var,ixRG^L)
endif
\}

if(oktest)write(*,*)'mpibound ipe,varnew=',ipe,var(ixtest^D,min(nvar,iwtest))

return
end

!==============================================================================
subroutine mpisend(nvar,var,ix^L,qipe,iside)

! Send var(ix^L,1:nvar) to processor qipe.
! jside is 0 for min and 1 for max side of the grid for the sending PE

include 'vacdef.f'

integer :: nvar
double precision :: var(ixG^T,nvar)
integer :: ix^L, qipe, iside, n, ix^D, ivar
!----------------------------------------------------------------------------
oktest = index(teststr,'mpisend')>0

n=0
do ivar=1,nvar
   {do ix^DB=ixmin^DB,ixmax^DB;}
      n=n+1
      sendbuffer(n)=var(ix^D,ivar)
   {enddo^DLOOP\}
end do

if(oktest)then
   write(*,*)'mpisend ipe-->qipe,iside,itag',ipe,qipe,iside,10*ipe+iside
   write(*,*)'mpisend ipe,ix^L,var=',ipe,ix^L,var(ixtest^D,min(iwtest,nvar))
endif

call MPI_RSEND(sendbuffer(1),n,MPI_DOUBLE_PRECISION,qipe,10*ipe+iside,&
               MPI_COMM_WORLD,ierrmpi)

return
end

!==============================================================================
subroutine mpirecvbuffer(nvar,ix^L,qipe,iside)

! receive buffer for a ghost cell region of size ix^L sent from processor qipe
! and sent from side iside of the grid

include 'vacdef.f'

integer:: nvar, ix^L, qipe, iside, n

integer :: nmpirequest, mpirequests(2)
integer :: mpistatus(MPI_STATUS_SIZE,2)
common /mpirecv/ nmpirequest,mpirequests,mpistatus
!----------------------------------------------------------------------------

oktest = index(teststr,'mpirecv')>0

n = nvar* ^D&(ixmax^D-ixmin^D+1)*

if(oktest)write(*,*)'mpirecv ipe<--qipe,iside,itag,n',&
   ipe,qipe,iside,10*qipe+iside,n

nmpirequest = nmpirequest + 1
call MPI_IRECV(recvbuffer(1,iside),n,MPI_DOUBLE_PRECISION,qipe,10*qipe+iside,&
               MPI_COMM_WORLD,mpirequests(nmpirequest),ierrmpi)

return
end

!==============================================================================
subroutine mpibuffer2var(iside,nvar,var,ix^L)

! Copy mpibuffer(:,iside) into var(ix^L,1:nvar)
include 'vacdef.f'

integer :: nvar
double precision:: var(ixG^T,nvar)
integer:: ix^L,iside,n,ix^D,ivar
!-----------------------------------------------------------------------------
oktest = index(teststr,'buffer2var')>0

n=0
do ivar=1,nvar
   {do ix^DB=ixmin^DB,ixmax^DB;}
      n=n+1
      var(ix^D,ivar)=recvbuffer(n,iside)
   {enddo^DLOOP\}
end do

if(oktest)write(*,*)'buffer2var: ipe,iside,ix^L,var',&
   ipe,iside,ix^L,var(ixtest^D,min(iwtest,nvar))

return
end
