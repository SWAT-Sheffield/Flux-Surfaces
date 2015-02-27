!##############################################################################
! include vacdef
IMPLICIT NONE

!HPF$ PROCESSORS PP(NUMBER_OF_PROCESSORS())

! DEFINITIONS OF GLOBAL PARAMETERS AND VARIABLES
! Parameters:

! Indices for cylindrical coordinates FOR TESTS, negative value when not used:
INTEGER,PARAMETER:: r_=1, phi_=^PHI, z_=^Z

! Indices for cylindrical coordinates FOR INDEXING, always positive
INTEGER,PARAMETER:: pphi_=^PPHI, zz_=^ZZ

include 'vacpar.f'

INTEGER,PARAMETER:: ixGlo^D=1
! The next line is edited by SETVAC
INTEGER,PARAMETER:: ixGhi1=68,ixGhi2=68,ixGhi3=36,ixGhimin=36,ixGhimax=68
INTEGER,PARAMETER:: ndim=^ND, ndir=^NC

INTEGER,PARAMETER:: dixBlo=2,dixBhi=2

{^IFPOISSON INTEGER,PARAMETER:: nwrk=7} !Size of work array for VACPOISSON
INTEGER,PARAMETER:: nhi=nw*{ixGhi^D*} ! Maximum number of unknowns for VACIMPL
{^IFIMPL 
! Size of work array for VACIMPL
INTEGER,PARAMETER:: nwork=(7+nw*(2*ndim+1))*{ixGhi^D*}*nw 
! Max. number of implicit variables
INTEGER,PARAMETER:: nwimplhi=nw 
} 

INTEGER,PARAMETER:: nhiB=10           ! maximum No. boundary sections

INTEGER,PARAMETER:: nsavehi=100       ! maximum No. saves into outputfiles
                                      ! defined by arrays of tsave or itsave

INTEGER,PARAMETER:: niw_=nw+1         !Indexname for size of iw index array

INTEGER,PARAMETER:: filelog_=1,fileout_=2,nfile=2 ! outputfiles

INTEGER,PARAMETER:: unitstdin=5,unitterm=6,uniterr=6,unitini=10 ! Unit names. 
                                  ! Outputfiles use unitini+1..initini+nfile
                                  ! Default parfiles uses unitini-1

INTEGER,PARAMETER:: biginteger=10000000

DOUBLE PRECISION,PARAMETER:: pi= 3.1415926535897932384626433832795
DOUBLE PRECISION,PARAMETER:: smalldouble=1.D-99, bigdouble=1.D+99
DOUBLE PRECISION,PARAMETER:: zero=0D0,one=1D0,two=2D0,half=0.5D0,quarter=0.25D0

INTEGER,PARAMETER:: toosmallp_=1,toosmallr_=2,couranterr_=3,poissonerr_=4
INTEGER,PARAMETER:: nerrcode=4

include 'vacusrpar.f'

{^IFMPI 
! Buffer size for sending and receiving ghost cells.
! Both w and x are sent with MPI.
INTEGER :: nmpibuffer, maxndimnw
! The next line is edited by SETVAC
PARAMETER (maxndimnw=nw)
PARAMETER (nmpibuffer=dixBhi*maxndimnw*{ixGhi^D*}/ixGhimin)
! MPI header file
include 'mpif.h'
}

!-- Common variables:

{^IFMPI COMMON, INTEGER:: ipe, ipe^D, npe, npe^D, nxall^D, nxpe^D, ierrmpi}
{^IFMPI COMMON, INTEGER:: ixPEmin^D, ixPEmax^D}
{^IFMPI COMMON, LOGICAL:: mpiupperB(ndim),mpilowerB(ndim)}
{^IFMPI COMMON, DOUBLE PRECISION:: sendbuffer(nmpibuffer)}
{^IFMPI COMMON, DOUBLE PRECISION:: recvbuffer(nmpibuffer,2)}

! Unit for reading input parameters.
COMMON, INTEGER:: unitpar

! Logical to set verbosity. For MPI parallel run only PE 0 is verbose
COMMON, LOGICAL:: verbose

! General temporary arrays, any subroutine call may change them 
! except for subroutines which say the opposite in their header
COMMON, DOUBLE PRECISION:: tmp(ixG^T),tmp2(ixG^T)

! Number of errors during calculation
COMMON, INTEGER:: nerror(nerrcode)

!Kronecker delta and Levi-Civita tensors
COMMON, INTEGER:: kr(3,3),lvc(3,3,3)

!Grid parameters
COMMON, INTEGER:: ixM^L,ixG^L,nx^D,nx(ndim)
COMMON, INTEGER:: dixB^L
! x and dx are local for HPF
COMMON, DOUBLE PRECISION:: x(IXG^T,ndim),dx(IXG^T,ndim)
COMMON, DOUBLE PRECISION:: volume,dvolume(IXG^T)
COMMON, DOUBLE PRECISION:: area(IXGLO1:IXGHI1),areaC(IXGLO1:IXGHI1)
COMMON, DOUBLE PRECISION:: areadx(IXGLO1:IXGHI1),areaside(IXGLO1:IXGHI1)

! Variables for generalized coordinates and polargrid
COMMON, LOGICAL::          gencoord, polargrid
{^IFGEN COMMON, DOUBLE PRECISION:: surfaceC(IXG^T,ndim),normalC(IXG^T,ndim,ndim)}
{^NOGEN COMMON, DOUBLE PRECISION:: surfaceC(2^D&,ndim), normalC(2^D&,ndim,ndim)}

!Boundary region parameters
COMMON, DOUBLE PRECISION:: fixB^D(-dixBlo:dixBhi^D%ixGLO^DD:ixGHI^DD,nw)
COMMON, INTEGER:: nB,ixB^LIM(ndim,nhiB),idimB(nhiB),ipairB(nhiB)
COMMON, LOGICAL:: upperB(nhiB),fixedB(nw,nhiB),nofluxB(nw,ndim),extraB
COMMON, CHARACTER*^LENTYPE :: typeB(nw,nhiB),typeBscalar(nhiB)

!Equation and method parameters
COMMON, DOUBLE PRECISION:: eqpar(neqpar+nspecialpar),procpar(nprocpar)

! Time step control parameters
COMMON, DOUBLE PRECISION:: courantpar,dtpar,dtdiffpar,dtcourant(ndim),dtmrpc
COMMON, LOGICAL:: dtcantgrow
COMMON, INTEGER:: slowsteps

! Parameters for the implicit techniques
{^IFPOISSON COMMON, DOUBLE PRECISION:: wrk(ixG^T,nwrk) } 
{^IFIMPL COMMON, DOUBLE PRECISION:: work(nwork) }
COMMON, INTEGER:: nwimpl,nimpl
COMMON, DOUBLE PRECISION:: implpar,impldiffpar,implerror,implrelax,impldwlimit
COMMON, INTEGER:: implrestart,implrestart2,impliter,impliternr,implmrpcpar
COMMON, CHARACTER*^LENTYPE :: typeimplinit,typeimpliter,typeimplmat
COMMON, LOGICAL:: implconserv,implnewton,implcentered,implnewmat
COMMON, LOGICAL:: implpred,impl3level,impljacfast,implsource

!Method switches
COMMON, INTEGER:: iw_full(niw_),iw_semi(niw_),iw_impl(niw_),iw_filter(niw_)
COMMON, INTEGER:: iw_vector(nvector+1),vectoriw(nw)
          ! The upper bound+1 in iw_vector avoids F77 warnings when nvector=0
COMMON, CHARACTER*^LENTYPE :: typefull1,typepred1,typeimpl1,typefilter1
COMMON, CHARACTER*^LENTYPE :: typelimited,typefct,typetvd,typeaxial
COMMON, CHARACTER*^LENTYPE :: typepoisson, typeconstrain
COMMON, CHARACTER*^LENTYPE :: typelimiter(nw),typeentropy(nw)
COMMON, CHARACTER*^LENTYPE :: typeadvance, typedimsplit, typesourcesplit
COMMON, LOGICAL:: dimsplit,sourcesplit,sourceunsplit,artcomp(nw),useprimitive
COMMON, LOGICAL:: divbfix,divbwave,divbconstrain,angmomfix,compactres,smallfix
COMMON, INTEGER:: idimsplit
COMMON, INTEGER:: nproc(nfile+2)
COMMON, DOUBLE PRECISION:: entropycoef(nw),constraincoef
COMMON, DOUBLE PRECISION:: smallp,smallpcoeff,smallrho,smallrhocoeff,vacuumrho
COMMON, DOUBLE PRECISION:: muscleta1,muscleta2,musclomega,acmcoef(nw),acmexpo
COMMON, LOGICAL:: acmnolim, fourthorder
COMMON, INTEGER:: acmwidth

!Previous time step and residuals
COMMON, DOUBLE PRECISION:: wold(ixG^T,nw),residual,residmin,residmax

! Flux storage for flux-CT and flux-CD methods !!! for MHD only !!! 
{^IFCT COMMON, DOUBLE PRECISION:: fstore(ixG^T,ndim) }

!Time parameters
COMMON, INTEGER:: step,istep,nstep,it,itmin,itmax,nexpl,nnewton,niter,nmatvec
COMMON, DOUBLE PRECISION:: t,tmax,dt,dtmin,cputimemax
COMMON, LOGICAL:: tmaxexact
COMMON, DOUBLE PRECISION:: tsave(nsavehi,nfile),tsavelast(nfile),dtsave(nfile)
COMMON, INTEGER:: itsave(nsavehi,nfile),itsavelast(nfile),ditsave(nfile)
COMMON, INTEGER:: isavet(nfile),isaveit(nfile)

!File parameters
COMMON, CHARACTER*^LENNAME :: filenameini,filenameout,filename(nfile)
COMMON, CHARACTER*^LENNAME :: fileheadini,fileheadout,varnames,wnames
COMMON, CHARACTER*^LENTYPE :: typefileini,typefileout,typefilelog
COMMON, LOGICAL::             fullgridini,fullgridout
COMMON, INTEGER::             snapshotini,snapshotout,isaveout

!Test parameters
COMMON, CHARACTER*^LENNAME :: teststr
COMMON, INTEGER:: ixtest1,ixtest2,ixtest3,iwtest,idimtest,ipetest^IFMPI
LOGICAL:: oktest    !This is a local variable for all subroutines and functions

COMMON, DOUBLE PRECISION:: maxviscoef

! end include vacdef
!##############################################################################
