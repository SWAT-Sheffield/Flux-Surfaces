!##############################################################################
! include vacusrpar - gravity
INTEGER,PARAMETER:: grav0_=neqpar+1
INTEGER,PARAMETER:: grav^D_=neqpar+^D+1
INTEGER,PARAMETER:: nu_=neqpar+4
INTEGER,PARAMETER:: nspecialpar=^ND+2

!CHARACTER*2 ,PARAMETER:: specialparname='nu'

{^IFONED   CHARACTER*5 ,PARAMETER:: specialparname='grav1'}
{^IFTWOD   CHARACTER*11,PARAMETER:: specialparname='grav1 grav2'}
{^IFTHREED CHARACTER*17,PARAMETER:: specialparname='grav1 grav2 grav3'}



! end include vacusrpar - gravity
!##############################################################################
