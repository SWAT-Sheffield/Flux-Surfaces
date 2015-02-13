Flux Surface Analysis
=====================

This repository contains the core code used to run parameter studies
using SAC. Primarily, it comprises: the SAC code, as a git
submodule in the `sac` directory; A configration system and various
helper scripts in `scripts` to manage data and parameters; and the
flux surface analysis code and descriptive IPython notebooks in the
`analysis` directory.

The first stop on using this repository should be the top level
`configure.py` and `run.py` scripts which can be used to set up and
execute various parts of the code.

Initial Conditions
------------------

This code uses initial condition data in binary form, the
mathematics in http://arxiv.org/abs/1305.4788 could be used to
recreate something with very similar properties, the implementation
of those algorithms is under development and will be included in
(pysac)[https://github.com/SWAT-Sheffield/pysac].

You can download the binary version of the initial conditions from
[FigShare](http://figshare.com/articles/SAC_Magnetohydrostatic_Background_Conditions/1308563)
by running `./run.py download ini`.

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
