Effects of Expansion Rate on MHD Wave Generation from a Logarithmic Spiral Photospheric Driver
==============================================================================================

This is a super speedy paper!!

Technical Notes
---------------

###Compiler Configurations

SAC can be compiled under most FORTRAN compilers, tested ones and their flags
are listed below:

#### GNU FORTRAN
compiler: gfortran
flags: -ffree-form

#### INTEL FORTRAN
compiler:
flags: -free -mcmodel=medium -O3 [-xAVX]

#### PGI FORTRAN
compiler: pgf90
flags: -Mfreeform -w=all

Notes:
------

* The configure script does not write wnames to vac.par in a form that fortran will read, therefore wnames are not in the .out file correctly and are re created from the config file in the GDF translation step.