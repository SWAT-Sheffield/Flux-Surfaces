!##############################################################################
! module vacini

!INCLUDE:vacnul.process.t
!=============================================================================
program vacini

include 'vacdef.f'

integer:: iw,ieqpar,ix^L
character*20 :: typeini
double precision:: w(ixG^T,1:nw),wpar(nw)
logical:: lastiw
!-----------------------------------------------------------------------------
{call mpiinit ^IFMPI}

verbose=.true. .and.ipe==0^IFMPI
if(verbose)then
   write(*,'(a)')'VACINI 4.52 configured to'
   write(*,'(a)')'  -d=33 -phi=0 -z=0 -g=68,68,36 -p=mhd -u=Slog'
   write(*,'(a)')'  -on=cd,rk,mpi'
   write(*,'(a)')'  -off=mc,fct,tvdlf,tvd,impl,poisson,ct,gencoord,resist'
   {^IFMPI write(*,'(a,i3,a)')'Running on ',npe,' processors'}
endif

! Some default values
t=zero; it=0; 
typefileout='ascii'; typefileini='auto'
snapshotini=0
fullgridini=.false.; fullgridout=.false.
gencoord=   .false.
! There are no ghost cells in VACINI except when "readmesh" is used.
dixB^L=0; ixMmin^D=ixGlo^D; 
! Test cell
ixtest^D=ixMmin^D;
! Read parameters from STDIN
unitpar=unitstdin

{^IFMPI
! MPI reads from a file 
unitpar=unitini-1
open(unitpar,file='vacini.par',status='old')
}

if(verbose)write(*,*)'Filename for new initial file:'
read(unitpar,'(a)')filename(fileout_)
{^IFMPI 
! Extract and check the directional processor numbers and indexes
! and concat the PE number to the output filename
call mpisetnpeDipeD(filename(fileout_))
}
if(verbose)write(*,*)'Fileheader:'
read(unitpar,'(a)')fileheadout
if(verbose)write(*,*)'Variable names, e.g. "x y rho m1 m2":'
read(unitpar,'(a)')varnames

call setheaderstrings

if(verbose)then
   write(*,*)'Select action by typing one of the following words: '
   write(*,*)'   test'
   write(*,*)'   typefileini,snapshotini,read,readmesh,readnext,',&
                    'typefileout,write,save'
   write(*,*)'   domain,grid,sheargrid,shiftgrid,polargrid,',&
                    'roundgrid,rotategrid'
   write(*,*)'   transpose,regrid,stretchgrid,stretchgridcent'
   write(*,*)'   polarvar,spherevar,cartesian,rotatevar,',&
                    'setvar,perturbvar,conserve,primitive,multiply,divide'
   write(*,*)'   uniform,shocktube,wave,wave1,special'
   write(*,*)'   eqpar,gencoord'
endif

do
   if(verbose)write(*,*)'Action:'
   read(unitpar,'(a)')typeini
   if(verbose)write(*,*)'> ',typeini
   select case(typeini)
      case('verbose')
          if(verbose)write(*,*)'Verbose:'
          read(unitpar,*)verbose
          verbose=verbose.and.ipe==0^IFMPI
      case('test')
          if(verbose)write(*,*)'Teststring:'
          read(unitpar,'(a)')teststr
          if(verbose)write(*,*)'ixtest, idimtest and iwtest:'
          read(unitpar,*) ixtest^D,idimtest,iwtest
          {^IFMPI 
          call mpiix(ixtest^D,ipetest)}
      case('typefileini')
          if(verbose)write(*,*)'Type of old initial file: ascii/binary/special'
          read(unitpar,'(a)')typefileini
      case('snapshotini')
          if(verbose)write(*,*)'Number of snapshot to be read:'
          read(unitpar,*)snapshotini
      case('read','readmesh')
          if(verbose)write(*,*)'Filename for old initial file:'
          read(unitpar,'(a)')filenameini
          {call mpisetnpeDipeD(filenameini) ^IFMPI}
          if(typeini=='readmesh')then
              if(verbose)write(*,*)&
                  'Specify boundary width for old initial file:'
              read(unitpar,*) dixB^L
              fullgridini=.true.
          endif
          call readfileini(w)
      case('readnext')
          snapshotini=snapshotini+1
          call readfileini(w)
      case('typefileout')
          if(verbose)write(*,*)'Type of new initial file: ascii/binary/special'
          read(unitpar,'(a)')typefileout
      case('write')
         call savefile(fileout_,w)
      case('save')
         call savefile(fileout_,w)
         close(unitini+fileout_)
         exit
      case('grid')
         call ini_grid(.true.,ixM^L)
      case('domain')
         call ini_grid(.false.,ixM^L)
      case('shiftgrid')
         {^IFMPI call die('shiftgrid not implemented for MPI')}
         call shiftgrid(ixM^L,w)
      case('sheargrid')
         {^IFMPI call die('sheargrid not implemented for MPI')}
         call sheargrid(ixM^L,w)
      case('polargrid')
         {^NOONED call makepolargrid(ixM^L)
         gencoord=.true.
         if(.false.)}call die('Polar grid is meaningless in 1D')
      case('spheregrid')
         {^IFTHREED call spheregrid(ixM^L)
         gencoord=.true.
         if(.false.)}call die('Spherical grid is meaningful in 3D only')
      case('roundgrid')
         {^NOONED call roundgrid(ixM^L)
         gencoord=.true. 
         if(.false.)}call die('Round grid is meaningless in 1D')
      case('rotategrid')
         call rotatevar(ixM^L,1,ndim,x)
         gencoord=.true.
      case('transpose')
         {^IFMPI call die('transposexy not implemented for MPI')}
         {^IFTWOD call transposexy(ixM^L,w)
         if(.false.)}call die('Transpose is implemented in 2D only.')
      case('regrid')
         {^IFMPI call die('regrid not implemented for MPI')}
         call regrid(ixM^L,w)
      case('stretchgrid')
         {^IFMPI call die('stretchgrid not yet implemented for MPI')}
         call stretchgrid(.true.,ixM^L)
      case('stretchgridcent')
         {^IFMPI call die('stretchgridcent not yet implemented for MPI')}
         call stretchgrid(.false.,ixM^L)
      case('polarvar')
         {^NOONED call polarvar(ixM^L,w)
         if(.false.)}call die('Polar variables are meaningless in 1D')
      case('spherevar')
         {^IFTHREED call spherevar(ixM^L,w)
         if(.false.)}call die('Spherical variables are meaningful '// &
                                'in 3D only')
      case('cartesian')
         {^NOONED call cartesian(ixM^L,w)
         if(.false.)} call die('Polar variables are meaningless in 1D')
      case('rotatevar')
         call rotatevar(ixM^L,1,nw,w)
      case('setvar')
         if(verbose)write(*,*)'Give ix limits:'
         read(unitpar,*) ix^L
         {^IFMPI 
         call mpiixlimits(ix^L)}
         do
            if(verbose)write(*,*) &
                'Variable index, variable value, lastiw (T/F)?'
            read(unitpar,*)iw,wpar(iw),lastiw
            w(ix^S,iw)=wpar(iw)
            if(lastiw)exit
         enddo
      case('perturbvar')
         call perturbvar(ixM^L,w)
      case('multiply')
         if(verbose)write(*,*)'Give multiplying factors for each variable:'
         read(unitpar,*) wpar(1:nw)
         do iw=1,nw
            w(ixM^S,iw)=w(ixM^S,iw)*wpar(iw)
         end do
      case('divide')
         if(verbose)write(*,*)'Give dividing factors for each variable:'
         read(unitpar,*) wpar(1:nw)
         do iw=1,nw
            w(ixM^S,iw)=w(ixM^S,iw)/wpar(iw)
         end do
      case('conserv','conserve')
         call conserve(ixM^L,w)
      case('primitive')
         call primitive(ixM^L,w)
      case('uniform')
         if(verbose)write(*,*)'Give values for each variable:'
         read(unitpar,*) wpar(1:nw)
         do iw=1,nw
            w(ixM^S,iw)=wpar(iw)
         end do
      case('shocktube')
         call ini_shocktube(ixM^L,w)
      case('wave')
         call ini_wave(ixM^L,w)
      case('wave1')
         call wave1(ixM^L,w)
      case('special')
         call specialini(ixM^L,w)
      case('eqpar')
         if(verbose)write(*,*)'Equation params:',neqpar+nspecialpar
         read(unitpar,*)(eqpar(ieqpar),ieqpar=1,neqpar+nspecialpar)
      case('gencoord')
         if(verbose)write(*,*)'Generalized coordinates (T/F):'
         read(unitpar,*)gencoord
      case default
         call die('Error in VACIni: no such action')
   end select
end do

{^IFMPI
close(unitpar)
call mpifinalize
}

end

!=============================================================================
subroutine ini_grid(coord,ix^L)

! Setup a uniform grid. When coord is .true., the user provides the coordinates
! for the centers of the grid, otherwise the boundaries of the computational
! domaines, thus the centers start at xmin+dx/2, and end at xmax-dx/2

include 'vacdef.f'

logical:: coord
integer:: ix^L,ix^D,idim
double precision:: dx^D,xmax(ndim),xmin(ndim)
!-----------------------------------------------------------------------------

if(verbose)write(*,'(a,3i6)')'Size of mesh. Max: ',ixGhi^D
read(unitpar,*) nx^D
if(verbose)write(*,'(a,3i6)')'Size of mesh: ',nx^D
if(coord)then
   if(verbose)write(*,*)'Coordinates of cell centers at the edges'
else
   if(verbose)write(*,*)'Boundaries of the computational domain'
endif
if(verbose)write(*,*)'xmin coordinates:'
read(unitpar,*)(xmin(idim),idim=1,ndim)
if(verbose)write(*,*)'xmax coordinates:'
read(unitpar,*)(xmax(idim),idim=1,ndim)

! Calculate cell sizes and modify coordinates for 'domain' action
if(coord)then
   dx^D=(xmax(^D)-xmin(^D))/(nx^D-1);
else
   dx^D=(xmax(^D)-xmin(^D))/nx^D;
   {^DLOOP
   xmax(^D)=xmax(^D)-dx^D/2
   xmin(^D)=xmin(^D)+dx^D/2
   }
endif

{^IFMPI
! Distribute global grid onto processor cube
nxall^D=nx^D;
call mpigridsetup
}

ixmax^D=ixmin^D+nx^D-1;
if(ixmax^D>ixGhi^D|.or.) call die('Error in IniGrid: Too big grid')

{^IFMPI
! Set coordinate limits for this PE
{^DLOOP
xmin(^D) = xmin(^D) + dx^D*(ixPEmin^D-1)
xmax(^D) = xmin(^D) + dx^D*(nx^D-1)
\}
}

{forall(ix^DD=ixmin^DD:ixmax^DD) x(ix^DD,^D)= &
   ((ix^D-ixmin^D)*xmax(^D)+ &
    (ixmax^D-ix^D)*xmin(^D)) /(ixmax^D-ixmin^D) \}

return
end

!=============================================================================
{^NOONED
subroutine makepolargrid(ix^L)

include 'vacdef.f'

integer:: ix^L
double precision:: pi2
!-----------------------------------------------------------------------------

if(verbose)write(*,*)&
  'First coordinate is interpreted as radius, second as angle/2pi.'

pi2=8*atan(one)

tmp(ix^S)=x(ix^S,1)
x(ix^S,1)=tmp(ix^S)*cos(x(ix^S,2)*pi2)
x(ix^S,2)=tmp(ix^S)*sin(x(ix^S,2)*pi2)

return
end
}
!=============================================================================
{^IFTHREED
subroutine spheregrid(ix^L)

include 'vacdef.f'

integer:: ix^L
double precision:: pi2
!-----------------------------------------------------------------------------

if(verbose)write(*,*)'R, PHI/2Pi [0,1], THETA/2Pi [-0.25,0.25] --> X,Y,Z'

pi2=8*atan(one)

tmp(ix^S)=x(ix^S,1)
x(ix^S,1)=tmp(ix^S)*cos(x(ix^S,3)*pi2)*cos(x(ix^S,2)*pi2)
x(ix^S,2)=tmp(ix^S)*cos(x(ix^S,3)*pi2)*sin(x(ix^S,2)*pi2)
x(ix^S,3)=tmp(ix^S)*sin(x(ix^S,3)*pi2)

return
end
}
!=============================================================================
{^NOONED 
subroutine roundgrid(ix^L)

! Calculate the shrink factor to shrink a rectangle to an ellipse. For a sguare
!   1               in directions parallel to x and y
!   1-r+r*sqrt(0.5) in diagonal directions, where r is the radial distance
!                      normalized to 1. 

include 'vacdef.f'

integer:: ix^L
double precision:: dist1(ixG^T),dist2(ixG^T),weight(ixG^T)
double precision:: xcent^D,rounded,squared
!-----------------------------------------------------------------------------

if(verbose)write(*,*)'Coordinates of center point:'
read(unitpar,*) xcent^D
if(verbose)then
   write(*,*)'Center of rounded grid:',xcent^D
   write(*,*)'If the rectangle is mapped to the (-1,-1,1,1) square'
   write(*,*)'give the distances to be rounded, and to be squared:'
endif
read(unitpar,*)rounded,squared

! Normalized distances in the L1 and L2 norms

dist1(ix^S)=max(^D&abs((x(ix^S,^D)-xcent^D)/(x(ixmax^DD,^D)-xcent^D)))
dist2(ix^S)=sqrt(^D&((x(ix^S,^D)-xcent^D)/(x(ixmax^DD,^D)-xcent^D))**2+)

! The weight increases from 0 to 1 for distance between 0 and rounded, and
! it drops back to 0 in the distance range rounded and squared. 
where(dist1(ix^S)<rounded)
   weight(ix^S)=dist1(ix^S)/rounded
elsewhere
   weight(ix^S)=(squared-dist1(ix^S))/(squared-rounded)
endwhere
weight(ix^S)=min(one,max(zero,weight(ix^S)))

! Shrink all coordinates by the factor 1 + weight*(dist1/dist2 - 1)
! Where weight is 0 there is no distortion, where weight is 1 the grid is round

weight(ix^S)=one + weight(ix^S)*(dist1(ix^S)/dist2(ix^S) - one)

^D&x(ix^S,^D)=weight(ix^S)*(x(ix^S,^D)-xcent^D)+xcent^D;

return
end
}
!=============================================================================
{^NOONED
subroutine polarvar(ix^L,w)

! Given r=x(:,1) and phi=2*Pi*x(:,2) 
! rotate the v_r,v_phi vector components to v_x,v_y

include 'vacdef.f'

integer:: ix^L,iw,jw
double precision:: w(ixG^T,nw),wi(ixG^T),wj(ixG^T)
double precision:: cosphi(ixG^T),sinphi(ixG^T),pi2
!-----------------------------------------------------------------------------

if(verbose)write(*,*)'Indices of var_r and var_phi:'
read(unitpar,*)iw,jw

pi2=8*atan(one)
cosphi(ix^S)=cos(x(ix^S,2)*pi2)
sinphi(ix^S)=sin(x(ix^S,2)*pi2)
wi(ix^S)=w(ix^S,iw)
wj(ix^S)=w(ix^S,jw)

w(ix^S,iw)=cosphi(ix^S)*wi(ix^S)-sinphi(ix^S)*wj(ix^S)
w(ix^S,jw)=cosphi(ix^S)*wj(ix^S)+sinphi(ix^S)*wi(ix^S)

return
end
}
!=============================================================================
{^IFTHREED
subroutine spherevar(ix^L,w)

! Given r=x(:,1), phi=2*Pi*x(:,2) [0,1], and theta=2*Pi*x(:,3) [-0.25,0.25]
! rotate the v_r,v_phi,v_theta vector components to v_x,v_y,v_z

include 'vacdef.f'

integer:: ix^L,iw,jw,kw
double precision:: w(ixG^T,nw),pi2
double precision, dimension(ixG^T):: cosphi,sinphi,costheta,sintheta,wi,wj,wk
!-----------------------------------------------------------------------------

if(verbose)write(*,*)'Indices of var_r, var_phi, var_theta:'
read(unitpar,*)iw,jw,kw

pi2=8*atan(one)
sinphi(ix^S)=sin(x(ix^S,2)*pi2)
cosphi(ix^S)=cos(x(ix^S,2)*pi2)
sintheta(ix^S)=sin(x(ix^S,3)*pi2)
costheta(ix^S)=cos(x(ix^S,3)*pi2)
wi(ix^S)=w(ix^S,iw)
wj(ix^S)=w(ix^S,jw)
wk(ix^S)=w(ix^S,kw)

w(ix^S,iw)=costheta(ix^S)*cosphi(ix^S)*wi(ix^S)-sinphi(ix^S)*wj(ix^S)&
          -sintheta(ix^S)*cosphi(ix^S)*wk(ix^S)

w(ix^S,jw)=costheta(ix^S)*sinphi(ix^S)*wi(ix^S)+cosphi(ix^S)*wj(ix^S)&
          -sintheta(ix^S)*sinphi(ix^S)*wk(ix^S)

w(ix^S,kw)=sintheta(ix^S)*wi(ix^S)+costheta(ix^S)*wk(ix^S)

return
end
}
!=============================================================================
{^NOONED
subroutine cartesian(ix^L,w)

! Given r=x(:,1) and phi=x(:,2) on a polar grid
! rotate the v_x,v_y vector components to v_r,v_phi

include 'vacdef.f'

integer:: ix^L,iw,jw
double precision:: w(ixG^T,nw),wi(ixG^T),wj(ixG^T)
!-----------------------------------------------------------------------------

if(verbose)write(*,*)'Indices of var_x and var_y:'
read(unitpar,*)iw,jw

wi(ix^S)=w(ix^S,iw)
wj(ix^S)=w(ix^S,jw)

w(ix^S,iw)=cos(x(ix^S,2))*wi(ix^S)+sin(x(ix^S,2))*wj(ix^S)
w(ix^S,jw)=cos(x(ix^S,2))*wj(ix^S)-sin(x(ix^S,2))*wi(ix^S)

return
end
}
!=============================================================================
subroutine shiftgrid(ix^L,wnew)

! Shifts grid relative to the variables. Padding is done by the border values.

include 'vacdef.f'

integer:: ix^L,ix,ixinside,dix,idim,iw
double precision:: w(ixG^T,nw),wnew(ixG^T,nw)
!-----------------------------------------------------------------------------

if(verbose)write(*,*)'Give direction and size of shift: idim, dix'
read(unitpar,*)idim,dix

w(ix^S,1:nw)=wnew(ix^S,1:nw)
select case(idim)
  {case(^D)
     do iw=1,nw
        do ix= ixmin^D,ixmax^D
           ixinside=min(ixmax^D,max(ixmin^D,ix-dix))
           wnew(ix^D%ix^S,iw)=w(ixinside^D%ix^S,iw)
        end do
     end do \}
  case default
    call die('Error in ShiftGrid: Unknown direction')
end select

return
end

!=============================================================================
subroutine sheargrid(ix^L,wnew)

! Shears grid relative to the variables. Padding is done by the border values.

include 'vacdef.f'

integer:: ix^L,ix,ix^D,dix,idim1,idim2
double precision:: w(ixG^T,nw),wnew(ixG^T,nw),angle
!-----------------------------------------------------------------------------

if(verbose)write(*,*)'Give directions and angle of shear: idim1,idim2,angle'
read(unitpar,*)idim1,idim2,angle
if(idim1==idim2)call die('Error in ShearGrid: idim1==idim2')
angle=angle*atan(one)/45

w(ix^S,1:nw)=wnew(ix^S,1:nw)
{do ix^D= ix^DL \}
   select case(idim2)
     {case(^D)
        dix=nint(-tan(angle)*(ix^D-ixmin^D)) \}
     case default
        call die('Error in ShearGrid: Unknown 2nd direction')
   end select
   select case(idim1)
      {case(^D)
        ix=min(ixmax^D,max(ixmin^D,ix^D-dix))
        wnew(ix^DD,1:nw)=w(ix^D%ix^DD,1:nw) \}
      case default
         call die('Error in ShearGrid: Unknown 1st direction')
   end select
{enddo^D&\}

return
end

!=============================================================================
subroutine multiply(exponent,ix^L,w)

! Multiply all variables by some function of the first coordinate

include 'vacdef.f'

integer:: exponent,ix^L,ix
double precision:: w(ixG^T,nw),r(ixG^LLIM1:)
!-----------------------------------------------------------------------------

if(verbose)write(*,*)'Define axial symmetry: cylinder/sphere/nozzle'
read(unitpar,'(a)')typeaxial

r(ix^LIM1:)=x(ix^LIM1:,ixmin^DE,1)

select case(typeaxial)
   case('cylinder')
      area(ix^LIM1:)=r(ix^LIM1:)
   case('sphere')
      area(ix^LIM1:)=r(ix^LIM1:)**2
   case('nozzle')
      if(verbose)write(*,*)"Warning in Multiply: This is Yee's nozzle problem!"
      area(ix^LIM1:) =1.398+0.347*(1-2/(exp(1.6*r(ix^LIM1:) -8)+1))
   case default
      call die('Error in Multiply: Unknown axial symmetry type!')
end select

forall(ix= ix^LIM1:) w(ix,ix^SE,1:nw)=w(ix,ix^SE,1:nw)*area(ix)**exponent

return
end

!=============================================================================
subroutine rotatevar(ix^L,iw^LIM,w)

! Rotate a pair of vector variables around some axis

include 'vacdef.f'

integer:: ix^L,iw^LIM,iw1,iw2
double precision:: w(ixG^T,iw^LIM:),w1(ixG^T),w2(ixG^T),angle
!-----------------------------------------------------------------------------

if(verbose)write(*,*)'Give indices of 2 variables and angle in degrees'
read(unitpar,*)iw1,iw2,angle
angle=angle*atan(one)/45
if(iw1==iw2.or.iw1<1.or.iw1>nw.or.iw2<1.or.iw2>nw)call die( &
   'Error in RotateVar: Incorrect iw1 and iw2')

w1(ix^S)=w(ix^S,iw1)
w2(ix^S)=w(ix^S,iw2)
w(ix^S,iw1)=cos(angle)*w1(ix^S)-sin(angle)*w2(ix^S)
w(ix^S,iw2)=sin(angle)*w1(ix^S)+cos(angle)*w2(ix^S)

return
end

{^IFTWOD
!=============================================================================
subroutine transposexy(ix^L,w)

! Transpose the first two coordinates of the grid (x), the variables (w), 
! and exchange the vector components as required

include 'vacdef.f'

integer:: ix^L,ixold^L,ix^D,iw,ivect,idim,qnvector
double precision:: w(ixG^T,1:nw)
!-----------------------------------------------------------------------------
ixold^L=ix^L;
ix^LIM1=ixold^LIM2; 
ix^LIM2=ixold^LIM1;
! Transpose x
do idim=1,ndim
   !!!For sake of f90tof77 x(ix^S,idim)=transpose(x(ixold^S,idim)) is replaced:
   tmp(ix^S)=x(ix^S,idim)
   {do ix^D=ixmin^D,ixmax^D\}
      x(ix^D,idim)=tmp(ix^DB)
   {enddo^D&\}
enddo

! Swap the X and Y coordinates
tmp(ix^S)=x(ix^S,1)
x(ix^S,1)=x(ix^S,2)
x(ix^S,2)=tmp(ix^S)

! Transpose w
do iw=1,nw
   !!!For sake of f90tof77 w(ix^S,iw)=transpose(w(ixold^S,iw)) is replaced by
   tmp(ix^S)=w(ix^S,iw)
   {do ix^D=ixmin^D,ixmax^D\}
      w(ix^D,iw)=tmp(ix^DB)
   {enddo^D&\}
enddo

! Swap the vector variables
! qnvector is only used to avoid compiler warning when nvector=0
qnvector=nvector
do ivect=1,qnvector
   if(verbose)write(*,"(a,i1,a)")&
     'Index of first component of vector variable #',ivect,':'
   read(unitpar,*)iw
   if(iw>=nw.or.iw<1)call die('Error in TransposeXY: Incorrect iw.')
   tmp(ix^S)=w(ix^S,iw)
   w(ix^S,iw)=w(ix^S,iw+1)
   w(ix^S,iw+1)=tmp(ix^S)
end do

return
end
}
!=============================================================================
subroutine regrid(ix^L,w)

! Change the number of grid points and extrapolate and interpolate the
! original cell positions x and averaged values w.

include 'vacdef.f'

integer:: ix^L,nix^D,ixnew^L,iw,idim
double precision:: w(ixG^T,1:nw),q(ixG^T)
!-----------------------------------------------------------------------------

if(verbose)write(*,*)'Define number of gridpoints in each direction:'
read(unitpar,*) nix^D

ixnewmin^D=ixmin^D;
ixnewmax^D=ixmin^D+nix^D-1;

do idim=1,ndim
   q(ix^S)=x(ix^S,idim)
   call regrid1(ix^L,ixnew^L,q)
   x(ixnew^S,idim)=q(ixnew^S)
enddo
do iw=1,nw
   q(ix^S)=w(ix^S,iw)
   call regrid1(ix^L,ixnew^L,q)
   w(ixnew^S,iw)=q(ixnew^S)
enddo

ixmax^D=ixnewmax^D;

return
end

!=============================================================================
subroutine regrid1(ixI^L,ixO^L,q)

! Interpolate q from grid determined by ixI to ixO. Use the distances measured
! between the Cartesian gridpoints, i.e. distances in generalized coordinates.
! The DOMAIN ixImin-0.5..ixImax+0.5 is rediscretized by ixOmax-ixOmin+1 points.

include 'vacdef.f'

integer:: ixI^L,ixO^L,ixI^D,ixO^D,dixI^D
double precision:: q(ixG^T),qnew(ixG^T),dxO^D,xO^D,coeff^D(0:1)
!-----------------------------------------------------------------------------

! Grid spacing of the output grid stretched onto the integer input grid
dxO^D=(ixImax^D-ixImin^D+one)/(ixOmax^D-ixOmin^D+one);

qnew(ixO^S)=zero
{do ixO^D=ixOmin^D,ixOmax^D\}

   ! Location of the output grid point
   xO^D=ixImin^D-half+(ixO^D-ixOmin^D+half)*dxO^D;

   ! Index of the input grid point to the left of xO within ixImin..ixImax-1
   ixI^D=min(ixImax^D-1,max(ixImin^D,int(xO^D)));

   ! Calculate bilinear interpolation/extrapolation coefficients
   coeff^D(1)=xO^D-ixI^D;
   coeff^D(0)=1-coeff^D(1);

   ! Interpolate q into qnew
   {do dixI^D=0,1\}
      qnew(ixO^D)=qnew(ixO^D)+(coeff^D(dixI^D)*)*q(ixI^D+dixI^D)
   {enddo^D&\}

{enddo^D&\}

q(ixO^S)=qnew(ixO^S);

return
end

!=============================================================================
subroutine stretchgrid(qdomain,ix^L)

! Stretch the grid logarithmically in direction idim segment by segment
! The original computational domain size is preserved if qdomain is true,
! and the first and last grid center locations are preserved if it is false.

include 'vacdef.f'

logical:: qdomain
integer:: ix^L,ix,ixL,ixR,idim,iseg,nseg
integer,parameter:: qixhi=10000
double precision:: qxL,qxR,qdxL,qdxR,qdxsum,qdx(qixhi)
!-----------------------------------------------------------------------------

if(verbose)write(*,*)'Direction of stretch and number of segments'
read(unitpar,*)idim,nseg

if(idim>ndim.or.idim<1)call die('Error in StretchGrid: Invalid direction')
if(nseg<1)call die('Error in StretchGrid: Invalid number of segments')

if(qdomain)then
   qxL=1.5*x(ixmin^D,idim)-0.5*x(ixmin^D+1,idim)
   qxR=1.5*x(ixmax^D,idim)-0.5*x(ixmax^D-1,idim)
   if(verbose)write(*,*)'Old domain from',qxL,' to',qxR
else
   qxL=x(ixmin^D,idim)
   qxR=x(ixmax^D,idim)
   if(verbose)write(*,*)'Old centers from',qxL,' to',qxR
endif

if(qxL>=qxR)&
   call die('Error in StretchGrid: qxR<qxL, use "domain" or "grid" action')

if(verbose)write(*,*)'xsecond-xfirst is scaled to 1 for segment 1'
select case(idim)
   {case(^D)
      if(ixmax^D>qixhi)call die( &
           'Error in StretchGrid: Too big grid, change qixhi')
      ixL=ixmin^D+1
      qdxL=one
      qdx(ixL)=qdxL
      if(qdomain)then
          qdxsum=1.5D0
      else
          qdxsum=1.0D0
      endif
      do iseg=1,nseg
         if(verbose)write(*,*)'xlast-xprev for segment:',iseg
         read(unitpar,*)qdxR
         if(iseg<nseg)then
            if(verbose)write(*,*)'Cell number of end of segment:',iseg
            read(unitpar,*)ixR
            if(ixR<=ixL.or.ixR>=ixmax^D) &
               call die('Error in StretchGrid: bad segment position')
         else
            if(verbose)write(*,*)'Last segment ends at right boundary'
            ixR=ixmax^D
         end if
         do ix=ixL+1,ixR
            qdx(ix)=qdxL*(qdxR/qdxL)**(dble(ix-ixL)/dble(ixR-ixL))
            qdxsum=qdxsum+qdx(ix)
         enddo
         qdxL=qdxR
         ixL=ixR
      end do
      if(qdomain)qdxsum=qdxsum+half*qdx(ixR)
      qdx(ixmin^D+1:ixmax^D)=qdx(ixmin^D+1:ixmax^D)*(qxR-qxL)/qdxsum
      if(qdomain)x(ixmin^D^D%ix^S,idim)=qxL+half*qdx(ixmin^D+1)
      do ix=ixmin^D+1,ixmax^D
         x(ix^D%ix^S,idim)=x(ix-1^D%ix^S,idim)+qdx(ix)
      enddo
   \}
end select

if(qdomain)then
   qxL=1.5*x(ixmin^D,idim)-0.5*x(ixmin^D+1,idim)
   qxR=1.5*x(ixmax^D,idim)-0.5*x(ixmax^D-1,idim)
   if(verbose)write(*,*)'New domain from',qxL,' to',qxR
else
   qxL=x(ixmin^D,idim)
   qxR=x(ixmax^D,idim)
   if(verbose)write(*,*)'New centers from',qxL,' to',qxR
endif

return
end

!=============================================================================
subroutine perturbvar(ix^L,w)

! Perturb a variable within limits ixP by adding the product of sine waves
! in each direction. The phases of the waves are relative to ixPmin.

include 'vacdef.f'

integer:: ix^L,ixP^L,idim,iw
double precision:: w(ixG^T,nw),dw,wavenum(ndim),phase(ndim)
!-----------------------------------------------------------------------------

if(verbose)write(*,*)'Variable index, amplitude and',&
   ' ixP^DL '
read(unitpar,*)iw,dw,ixP^DL
{^IFMPI 
call mpiixlimits(ixP^L)}
if(verbose)write(*,*)'wavenumber and phase for each idim:'
read(unitpar,*)(wavenum(idim),phase(idim),idim=1,ndim)

w(ixP^S,iw)=w(ixP^S,iw)+dw* &
   {sin(phase(^D)+(x(ixP^S,^D)-x(ixPmin^DD,^D))*wavenum(^D))*}

return
end

!=============================================================================
subroutine ini_shocktube(ix^L,w)

! The shocktube is divided into nseg segments in the chosen idim direction.
! Linear interpolation in segments with given left and right states.

include 'vacdef.f'

integer:: ix^L
double precision:: w(ixG^T,nw)
integer:: ix,ixL,ixR,idim,iw,iseg,nseg,ieqpar
double precision:: wL(nw),wR(nw)
{^IFMPI logical:: inside}
!-----------------------------------------------------------------------------

if(verbose)write(*,*)'Normal of slab symmetry'
read(unitpar,*)idim
if(verbose)write(*,*)'Number of segments'
read(unitpar,*)nseg
if(verbose)write(*,*)'All variables at minimal position:'
read(unitpar,*)wL(1:nw)
select case(idim)
   {case(^D)
      ixL=ixmin^D
      {ixL=ixL-ixPEmin^D+1 ^IFMPI}
      do iseg=1,nseg
         if(verbose)write(*,*)'All variables at end of segment:',iseg
         read(unitpar,*)wR(1:nw)
         if(iseg<nseg)then
            if(verbose)write(*,*)'Cell number of end of segment:',iseg
            read(unitpar,*)ixR
            {ixR=ixR-ixPEmin^D+1 ^IFMPI}
            if(ixR<=ixL) &
               call die('Error in ini_shocktube: bad segment position')
         else
            if(verbose)write(*,*)'Last segment ends at right boundary'
            {ixR=ixmax^D ^NOMPI}
            {ixR=nxall^D-ixPEmin^D+1 ^IFMPI}
         end if
         do iw=1,nw
            do ix=ixL,ixR
              if(ix>=ixmin^D.and.ix<=ixmax^D) &
                  w(ix^D%ix^S,iw)=((ixR-ix)*wL(iw)+(ix-ixL)*wR(iw))/(ixR-ixL)

              !w(ix^D%ix^S,iw)=((x(ixR^D%ix^S,^D)-x(ix^D%ix^S,^D))*wL(iw)+ &
              !                 (x(ix^D%ix^S,^D)-x(ixL^D%ix^S,^D))*wR(iw))&
              !                /(x(ixR^D%ix^S,^D)-x(ixL^D%ix^S,^D)) 
            end do
         end do
         ixL=ixR
         wL(1:nw)=wR(1:nw)
      end do \}
   case default
      call die('Error in Ini_ShockTube: Unknown dimension')
end select

if(verbose)write(*,*)'Eqpar:'
read(unitpar,*)(eqpar(ieqpar),ieqpar=1,neqpar+nspecialpar)

return
end

!=============================================================================
subroutine ini_wave(ix^L,w)

! Sum of sine waves in each direction and each variable

include 'vacdef.f'

integer:: ix^L,ix^D,idim,iw,ieqpar
double precision:: w(ixG^T,nw),wmean(nw),wavenum(ndim),ampl(ndim),phase(ndim)
logical:: nextiw
!-----------------------------------------------------------------------------

if(verbose)write(*,*)'Mean values:'
read(unitpar,*)wmean(1:nw)
do iw=1,nw
   w(ix^S,iw)=wmean(iw)
end do
do iw=1,nw
   if(verbose)write(*,*)'Variable iw:',iw
   do
      if(verbose)write(*,*)'Amplitude,wavenum,phase for each dim,nextiw=T/F:'
      read(unitpar,*)(ampl(idim),wavenum(idim),phase(idim),idim=1,ndim),nextiw
      do idim=1,ndim
          w(ix^S,iw)=w(ix^S,iw)+ampl(idim)*sin(phase(idim)+&
            x(ix^S,idim)*wavenum(idim))
      enddo
      if(nextiw)exit
   enddo
enddo

if(verbose)write(*,*)'Eqpar:'
read(unitpar,*)(eqpar(ieqpar),ieqpar=1,neqpar+nspecialpar)

return
end

!=============================================================================
subroutine wave1(ix^L,w)

! Sine waves with arbitrary wave vectors and shifts using rationalized angle
! units (1.0 = full circle)

include 'vacdef.f'

integer:: ix^L,idim,iw
double precision:: w(ixG^T,nw),wmean,ampl,wavenum(ndim),phase,pi2
logical:: lastiw
!-----------------------------------------------------------------------------

pi2=8*atan(one)
if(verbose)&
   write(*,*)'w(iw)=wmean+dw*sin(2*pi*[x*kx+y*ky+phase]) (Note the 2*pi!)'
do
   if(verbose)write(*,*)'Give iw,wmean,dw,k(idim),phase,lastiw (quit with T):'
   read(unitpar,*)iw,wmean,ampl,(wavenum(idim),idim=1,ndim),phase,lastiw
   w(ix^S,iw)=wmean+ampl*sin(pi2*({wavenum(^D)*x(ix^S,^D)+}+phase))
   if(lastiw)exit
enddo

return
end
!=============================================================================
! Some interface routines for subroutines often used in the VACUSR module
! to keep the compiler happy 
!=============================================================================
subroutine gradient(realgrad,q,ix^L,idir,gradq)
logical:: realgrad
integer:: ix^L,idir
double precision:: q(*),gradq(*)
call die('Error: VACINI cannot call gradient !')
end

!=============================================================================
subroutine laplace4(q,ix^L,laplaceq)
integer:: ix^L
double precision:: q(*),laplaceq(*)
call die('Error: VACINI cannot call laplace4 !')
end

!=============================================================================
subroutine ensurebound(dix,ixI^L,ixO^L,qt,w)
integer:: dix,ixI^L,ixO^L
double precision:: qt,w(*)
call die('Error: VACINI cannot call ensurebound !')
end
!=============================================================================
! end module vacini
!##############################################################################


