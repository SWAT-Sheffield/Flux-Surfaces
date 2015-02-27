############################################################################
#
#                    VAC Makefile
#
############################################################################
#
# For the first time read carefully everything and do selections. 
# Later you will probably need to change the compiler flags only.
# You may delete compiler flags for other machines once everything works.
#
# Select compiling options by adding and removing the leading # characters.
#
############################################################################

# The Makefile should always use the sh shell.
SHELL = /bin/sh

# The location of the main VAC directory. If src/ is a symbolic link,
# you need to define it explicitly, e.g. VACDIR=/u/smith/vac
# otherwise the default VACDIR=.. definition is fine
VACDIR=..

# This is here for safety if you called "make" without arguments.
nothing	:
	@echo make what?

############################################################################
# FLAGS FOR F90, F77, AND HPF COMPILERS
#
# The variables F90, F90FLAGS, F77, F77FLAGS, HPF, and HPFFLAGS are set 
# for the respective compilers if those compilers exist and are to be used.
# Uncomment and select the definitions for your machine and REMOVE 
# (or at least comment out) all definitions meant for other machines.

##### Nag F95 compiler (on a Linux box) #########################

#F77=f95
F77=pgf77
F77FLG=-132 -dusty -Mnobounds

#F90=f95
F90=pgf90
F90FLG=-Mfreeform -Mnobounds

#
# Select one:
#
F77FLAGS= $(F77FLG) -O4 -w=all #        final compile
#F77FLAGS= $(F77FLG) -O1 -C -w=all #    developing
#F77FLAGS= $(F77FLG) -g -C #            debugging

#
# Select one:
#
F90FLAGS= $(F90FLG) -O4 -w=all #        optimized
#F90FLAGS= $(F90FLG) -O1 -C -w=all #    develop
#F90FLAGS= $(F90FLG) -O0 -gline -C #            debugging

##### GNU Fortran compiler: g77 (e.g. for Linux) ######
#
# Replace F77_ with G77_ in the function cputime in src/vac.t to allow timing!
#
#F77=g77
#F77FLG=-ffixed-line-length-none 		#no maximum line length
# F77FLG=-ffixed-line-length-none -fno-automatic #no maximum line length,static

# The -fno-automatic flag may be needed to avoid memory allocation problems 
#
# Select one:
#
#F77FLAGS= $(F77FLG) -O2 # 			optimized 
# F77FLAGS= $(F77FLG) -funroll-loops -O2 #	more optimized 
# F77FLAGS= $(F77FLG) -m486 -funroll-loops -O2#	optimized for i486 or i586
# F77FLAGS= $(F77FLG) -O0 #			no optimization (save memory)
# F77FLAGS= $(F77FLG) -O1 #			some optimization (save memory)
# F77FLAGS= $(F77FLG) -Wuninitialized #		developing
# F77FLAGS= $(F77FLG) -Wall -g # 		debugging, being "picky"

# VAST/F90 compiler with g77
#F90=f90
#F90FLG=	-Wv,-xa # 				Free format
# F90FLAGS= -O2 $(F90FLG)  #
#F90FLAGS= -keep -O2 $(F90FLG)  #

###### SGI ######
# Typical definitions for the SGI F77 compiler
#
#F77=f77
#F77FLG= -u -extend_source -static
#
# Using the F90 compiler instead of F77:
#
#F77=f90
#F77FLG= -u -extend_source -woff136

# Typical definitions for the SGI F90 compiler
#
#F90=f90
#F90FLG= -u -freeform -extend_source -woff136

# FOR AUTO PARALLEL : setenv OMP_NUM_THREADS 4
# and use the flags below

#F77=f90
#F77FLG= -extend_source -pfa keep -static -woff136

#F90=f90
#F90FLG= -u -freeform -extend_source -pfa keep -static -woff136

# AUTO PARALLEL on SGI origin 3800 (TERAS) with
# automatic parallelization option (apo); setenv OMP_NUM_THREADS 4
# Note: type of processor to be set with -Ofast can be obtained using
#       command hinv
#F90=f90
#F90FLG= -u -freeform -extend_source
#F90FLAGS=  $(F90FLG) -Ofast=ip35 -apo  
#F77=f90
#F77FLG= -u -fixedform -extend_source
#F77FLAGS= $(F77FLG) -Ofast=ip35 -apo  

#
# Select one:
#
#F77FLAGS= $(F77FLG) -O2 #		final compile
# F77FLAGS= $(F77FLG) -O2 -p #		profiling
# F77FLAGS= $(F77FLG) #			developing
# F77FLAGS= $(F77FLG) -g -C -TENV:X=0 #		debugging

#
# Select one:
#
#F90FLAGS=  $(F90FLG) -O2 #                            optimized
#F90FLAGS=  $(F90FLG) -Ofast=ip27 #                      optimized

#
# for F77 libraries
#
#F77=f90
#F77FLG= -u -fixedform -extend_source
#
# Select one:
#
#F77FLAGS= $(F77FLG) -O2 #                     final compile
#F77FLAGS= $(F77FLG) -Ofast=ip27 #               final compile

###### DEC ######
# Typical definitions for DEC Fortran 77 and 90 compilers
#
# F90= f90
# F90FLG= -u -free # -convert big_endian # convert binary for other machines
#
# Select one:
#
# F90FLAGS=  $(F90FLG) -O4 #				optimized
# F90FLAGS=  $(F90FLG) -O4 -fast -transform_loops #	fast and optimized

# F90FLAGS=  $(F90FLG) -O5 -fast -transform_loops # -feedback vacall.feed 
#very optimized, use it with make vacall, feedback file can be produced by
#pixie vacall
#vacall.pixie < par/parfile
#prof -pixie -feedback src/vacall.feed

# F90FLAGS=  $(F90FLG) -O4 -automatic #	small compile (does not seem to work)
# F90FLAGS=  $(F90FLG) -O4 -p #		profiling
# F90FLAGS=  $(F90FLG) -C #		developing
# F90FLAGS=  $(F90FLG) -C -g #		debugging
#
# F77= f77
# F77FLG= -u -extend_source
#
# Select one:
#
# F77FLAGS= $(F77FLG) -O4 #		final compile
# F77FLAGS= $(F77FLG) -O4 -p #		profiling
# F77FLAGS= $(F77FLG) #			developing
# F77FLAGS= $(F77FLG) -g -C #		debugging

###### SUN ######
# Typical definitions for the SUN f77 compiler
#
# F77= f77
# F77FLG= -u -e -w
#
# Select one:
#
# F77FLAGS= $(F77FLG) -fast #		final compile
# F77FLAGS= $(F77FLG) -fast -p #	profiling
# F77FLAGS= $(F77FLG) #			developing
# F77FLAGS= $(F77FLG) -g -C #		debugging

###### HP ######
# Typical definitions for HP Fortran 90 compiler (according to G. Bodo)
# Select "F90TOF90= | forall2do" below
#
#F90= f90
#F90FLG= +implicit_none +source=free
#
# Select one:
#
# F90FLAGS=  $(F90FLG) #                Development
# F90FLAGS=  $(F90FLG) -O #             Optimized
#
# Typical definitions for HP Fortran  77 compiler
#
#F77= f77
#F77FLG= -u +es
#
# Select one:
#
#F77FLAGS= $(F77FLG) #                  Development
#F77FLAGS= $(F77FLG) -O #               Optimized

###### IBM RS/6000 or IBM SP ######
# Typical definitions for the AIX compilers
# Select F90TOF90= | forall2do below   !
# Remove FLUSH statements from vacio.t !
# 
# HPF = pghpf
# HPFFLG= -Mfreeform
# Select one:
# HPFFLAGS= $(HPFFLG) -O3 #	final compile (-O4 has a bug)
# HPFFLAGS= $(HPFFLG) #		development
#
# HPF = xlhpf90
# HPFFLG= -u
# HPFFLAGS= $(HPFFLG) -bmaxdata:0x80000000 -qintsize=8 -q64 -O3
# HPFFLAGS= $(HPFFLG) -bmaxdata:0x70000000 -qhpf=overlap=4 -O3

# F90 = xlf90
# F90FLG= -u 
# Select one:
# F90FLAGS= $(F90FLG) -O3 #	final compile (-O4 has a bug in SELECT CASE)
# F90FLAGS= $(F90FLG) #		development
#
# F77 = f77
# F77FLG= -u -qfixed=80 # 
# Select one:
# F77FLAGS= $(F77FLG) -O2 #	final compile
# F77FLAGS= $(F77FLG) #		development

##### Cray C90 ######
# Select F90TOF90= | forall2do below.
# Note that BLAS and LAPACK libraries are linked by default.
#
# F90= f90
# F90FLG= -dp -f free #			Double precision=real*8, free format
#
# Select one:
#
# F90FLAGS=  $(F90FLG) -O3 #		final compile on a single CPU
# F90FLAGS=  $(F90FLG) -O scalar3,vector3,task3,inline3 # optimized multitask
# F90FLAGS=  $(F90FLG) -O1 -G1 -l prof#	profiling, read "man prof" on Cray
# F90FLAGS=  $(F90FLG) -O0 #		developing
# F90FLAGS=  $(F90FLG) -O0 -m 0 -en -g#	debugging
#
# F77= cf77
# F77FLG= -N80 -dp #			80 columns, double precision=real*8
#
# Select one:
#
# F77FLAGS= $(F77FLG) -O1 #		final compile
# F77FLAGS= $(F77FLG) -O1 -G -l prof#	profiling, read "man prof" on Cray
# F77FLAGS= $(F77FLG) -O0 #		developing
# F77FLAGS= $(F77FLG) -O0 -g #		debugging (Adaptor does not allow -C)

##### Cray T3E ######
# Select F90TOHPF= | f90tohpf -distr=2 -pghpf below 
# for a better speed and to avoid a PGHPF 2.3-1 compiler bug.

# HPF = pghpf
# HPFFLG= -dp -Mfreeform -Mdclchk

# Select one:
# HPFFLAGS= $(HPFFLG) -O3 #	final compile
# HPFFLAGS= $(HPFFLG) #		development

# Select F90TOF90= | forall2do below.
#
# F90= f90
# F90FLG= -dp -eI -f free #		Double precision=real*8, free format
#
# Select one:
#
# F90FLAGS=  $(F90FLG) -O3 #		final compile
# F90FLAGS=  $(F90FLG) -O3 -eA -lapp #	final compile with Apprentice enabled
# F90FLAGS=  $(F90FLG) -O1 -G1 -l prof#	profiling, read "man prof" on Cray
# F90FLAGS=  $(F90FLG) -O0 -R abcn #	developing
# F90FLAGS=  $(F90FLG) -O0 -R abcn -m 0 -en -g#	debugging
#
# F77= f90 
# F77FLG= -dp -N80 #			Double precision=real*8, fixed format
#
# Select one:
#
# F77FLAGS= $(F77FLG) -O3 #		final compile
# F77FLAGS= $(F77FLG) -O1 -G -l prof#	profiling, read "man prof" on Cray
# F77FLAGS= $(F77FLG) -O0 -R abcn #	developing
# F77FLAGS= $(F77FLG) -O0 -R abcn -m 0 -en -g#	debugging 

############################################################################
# COMPILATION DEFAULTS FOR THE GENERIC "make vac" AND "make vacini" COMMANDS
#
# Select these if your default compiler for making VAC and VACINI is F90:
#
FOR=$(F90)
FORFLG=$(F90FLAGS)
PREFOR=$(F90TOF90)

# Select these if your default compiler for making VAC and VACINI is F77:
#
# FOR=$(F77)
# FORFLG=$(F77FLAGS)
# PREFOR=$(F90TOF77)

# Select these if your default compiler for making VAC and VACINI is GCC:
# (Assumes that an F2C translator has been used!)
#
# FOR=$(GCC)
# FORFLG=$(GCCFLAGS)
# PREFOR=$(F90TOF77)

# Select these if your default compiler for making VAC and VACINI is HPF.
# However, the VACIMPL module should be off. Select VACPP = vacpp.pl -O below.
#
# FOR=$(HPF)
# FORFLG=$(HPFFLAGS)
# PREFOR=$(F90TOHPF)

# DEFAULT COMPILER FOR THE FORTRAN 77 LIBRARIES AND CODES
# Select these if you want them to be compiled by a Fortran 77 compiler.
# If you compile VAC and VACINI with Fortran 77, you may want to set
# FORFLGLIB explicitly (e.g. to -O4 or -fast) independent of F77FLAGS.
#
 FORLIB=$(F77)
 FORFLGLIB=$(F77FLAGS)

# Select these if you use a C compiler with Fortran to C translator
# FORLIB=$(GCC)
# FORFLGLIB=$(GCCFLAGS)

# If you want to use another compiler (e.g. f90) make sure that the flags
# are set  corectly for fixed format Fortran 77. An example could be this:
#
# FORLIB=$(F90)
# FORFLGLIB=-u -O -fixed

############################################################################
# PREPROCESSOR OPTIONS
#
#  Some Fortran compilers cannot handle long lines in the source code.
#  The compiler flags -e or -extended_source usually solve this problem. 
#  Use the -maxlen flag for the precompiler ONLY IF IT IS REALLY NEEDED. 
#  Then the precompilers will limit the line length to maxlen characters.
#  The resulting preprocessing commands should be one of the following:
#
#  vacpp.pl -maxlen=72 FILENAME.t > FILENAME.f # short F90
# or
#  vacpp.pl -maxlen=72 FILENAME.t | forall2do -maxlen=72 > FILENAME.f
#							  # short F90 no FORALL
# or
#  vacpp.pl FILENAME.t | f90tof77 -maxlen=72 > FILENAME.f # short F77
# or
#  vacpp.pl -O -maxlen=72 FILENAME.t | f90tohpf > FILENAME.f # short HPF  

############### VACPreProcessor ################
#For the dependency

VACPP_=vacpp.pl

# ... select one for the compilation rules:
#
VACPP =vacpp.pl #			
#VACPP =vacpp.pl -O #			Optimized shifts (always used by HPF)
#VACPP =vacpp.pl -maxlen=72 #		Reduced maximum line length
#VACPP =vacpp.pl -O -maxlen=72 #	Optimization and reduction 

### Do you use a F90 compiler that does not know FORALL? (e.g. C90,SP) ###

F90TOF90 = #				No
#F90TOF90 = | forall2do #		Yes and use default maximum line length
#F90TOF90 = | forall2do -maxlen=72 #	Yes and use different maximum length

############ Maximum line length for Fortran 77 compiler ###################

F90TOF77= | f90tof77
#F90TOF77= | f90tof77 -maxlen=72 

##################### HPF Preprocessor ##############################

# For the PGHPF 2.3-1 compiler and earlier versions use the -pghpf flag, 
# which drops "ONTO PP" in the directives and thus avoids a compiler bug.
# Distributing the second dimension -distr=2 is preferred on the Cray T3E.
# If more than one dimensions are distributed and the -pghpf flag is not set, 
# the !HPF$ PROCESSORS directive in src/vacdef.t should be set correctly. 

# Select one and modify flags if needed

F90TOHPF= | f90tohpf #			Distribute 1st dimension only    (HPF)
# F90TOHPF= | f90tohpf -distr=2 -pghpf #Distribute 2nd dimension only  (PGHPF)
# F90TOHPF= | f90tohpf -distr=23 #	Distribute 2nd and 3rd dimension (HPF)
# F90TOHPF= | f90tohpf -distrstatic #   Distribute static arrays too     (HPF)

##############################################################################
# FORTRAN 77 LIBRARIES

###### BLAS AND LAPACK #####
#
# BLAS   (Basic Linear Algorithms Library) and
# LAPACK (Linear Algebra PACKage) are installed at most places. 
# VAC comes with the necessary BLAS and LAPACK routines, 
# but you may want to use your local precompiled libraries.

# Select one and modify directory and file names if needed:

#BLASLIB =        # 	BLAS is linked by default (e.g. on the Cray C90)
#BLASLIB = -lblas # 	Use local BLAS library (e.g. /usr/local/lib/libblas.a)
BLASLIB = blas$O # 	Use BLAS provided with VAC

#LAPACKLIB =      # 	LAPACK is linked by default (e.g. on the Cray C90)
#LAPACKLIB = -llapack #	Use local LAPACK library
LAPACKLIB = lapack$O #	Use LAPACK provided with VAC

#### General dependency for executable targets

# THe definitions for LIBS and LIBS_ are commented out by setvac 
# unless the implicit module is on (setvac -on=impl).
# Libraries (and dependencies) 
# to be linked to VAC if VACIMPL module is on:

#LIBS= conjgrad$O gmresr$O tridiag$O $(LAPACKLIB) $(BLASLIB)
#LIBS_ = conjgrad$O gmresr$O tridiag$O lapack$O blas$O 

# Select one:

#MAKE_ = Makefile # 	Recompile if Makefile changed
MAKE_ = #		Do not recompile if Makefile changed

############################################################################
# EXTENSION NAMES FOR SOURCE AND OBJECT FILES
#
# $F is the suffix for the generic VAC and VACINI precompiled Fortran files
# $(F7) is the suffix for the Fortran 77 codes, usually .f
# $O is the suffix for the object files, usually .o
F=.f
F7=.f
O=.o
# CHANGE THEM ONLY IF YOUR MACHINE DOES NOT ACCEPT THESE, BUT FIRST DO 
# "make clean" to remove any precompiled Fortran files with wrong extensions.
# Check the manual pages for your F90 and/or F77 compilers to see
# what extensions they would accept.
#
# IF YOU CHANGED THE EXTENSION FOR THE VAC AND VACINI FORTRAN FILES,
# i.e. the F variable, then type "make extensions"
# to correct filenames in the precompiler scripts.

extensions:
	$(VACDIR)/Fix/extensionfix -stem=vacdef,vacpar,vacusrpar -ext=$(F) \
		vacpp.pl savemem

# IF YOU CHANGED THE EXTENSION FOR FORTRAN 77 CODES, 
# i.e. the F7 variable, then set OLDF7 to the old suffix (including the dot)
# and type "make extension77" to move the Fortran 77 source files.
OLDF7=.f
extension77:
	mv convertdata$(OLDF7)  convertdata$(F7)
	mv convertavsdx$(OLDF7) convertavsdx$(F7)
	mv convertnul$(OLDF7)   convertnul$(F7)
	mv convertusr$(OLDF7)   convertusr$(F7)
	mv convertdef$(OLDF7)   convertdef$(F7)
	mv gmresr$(OLDF7) gmresr$(F7)
	mv tridiag$(OLDF7) tridiag$(F7)
	mv conjgrad$(OLDF7) conjgrad$(F7)
	mv lapack$(OLDF7) lapack$(F7)
	mv blas$(OLDF7) blas$(F7)
	$(VACDIR)/Fix/extensionfix -stem=convertdef -ext=$(F7) \
		convertdata$(F7) convertavsdx$(F7) convertusr$(F7)
	$(VACDIR)/Fix/extensionfix -stem=vacpp -ext=$(F7) setvac

########### Conversion of VAC data file formats
#
# Conversion to AVS and DX formats requires the LOVER library (C) S.Belien).
# See Man/convert.html for installation procedure and usage.

# Select one
AVSDX=convertavsdx$(F7) #	Yes, I want to convert to AVS and DX formats
#AVSDX=convertnul$(F7) # 	No need to convert to AVS and DX formats

# Select one according to the above choice
AVSDXLIB=-L$(VACDIR)/AVSDX/lib -llover # 	Yes
#AVSDXLIB= # 					No

convertdata: convertdata$(F7) $(AVSDX) convertusr$(F7) convertdef$(F7) Makefile
	$(F77) $(F77FLAGS) -o $(VACDIR)/convertdata \
		convertdata$(F7) convertusr$(F7) $(AVSDX) $(AVSDXLIB)

# On the Cray C90 replace "$(F77) $(F77FLAGS)" by "$(F90) -dp" above
# to avoid a cf77 compiler restriction. On the T3E use "$(F90) -sdefault32"

##############################################################################
######                                                                 #######
######   YOU PROBABLY DO NOT NEED TO CHANGE ANYTHING BELOW THIS LINE   #######
######   BUT IT MAY BE WORTHWHILE TO READ IT                           #######
######                                                                 #######
##############################################################################

################# Definitions for source files #############################

INCLUDES = vacdef$F vacpar$F vacusrpar$F $(MAKE_)

VAC_FOR = vac$F vacio$F vacgrid$F vacphys0$F vacphys$F vacusr$F
VAC_OBJ = vac$O vacio$O vacgrid$O vacphys0$O vacphys$O vacusr$O

# The VACCD, VACMC, VACFCT, VACTVDLF, VACTVD, VACIMPL, VACPOISSON, VACMPI 
# modules are removed or added in the following two lines by "setvac". 
# DO NOT TOUCH THESE TWO LINES:
VACFOR = $(VAC_FOR) vaccd$F vacmpi$F
VACOBJ = $(VAC_OBJ) vaccd$O vacmpi$O

# The VACMPI module is removed or added the following two lines by "setvac".
# DO NOT TOUCH THESE TWO LINES:
VACINIFOR = vacini$F vacio$F vacphys0$F vacusr$F vacmpi$F
VACINIOBJ = vacini$O vacio$O vacphys0$O vacusr$O vacmpi$O

ROEOBJ = roetest$O vacphys0$O vacphys$O

################### Translation and dependency rules #######################

PREPROC= $(VACPP_) forall2do f90tof77 f90tohpf

.SUFFIXES:
.SUFFIXES: .t $F  $O

$(VACFOR) vacini$F vacdef$F vacpar$F vacusrpar$F		: $(PREPROC)
roetest$F 							: $(PREPROC)

$(VACOBJ) vacini$O vacall$O vaciniall$O vacsmall$O roetest$O	: $(INCLUDES)

# General precompilation rule
.t$(F):
	$(VACPP) $< $(PREFOR) > $@

# General Fortran compilation rule
$(F)$(O):
	$(FOR) $(FORFLG) -c $<

# Compilation rule for Fortran 77 libraries and Fortran vacpp
$(LIBS_) nothing$O:
	$(FORLIB) $(FORFLGLIB) -c $*$(F7)

# Dependencies for libraries
blas$O:		blas$(F7)
lapack$O:	lapack$(F7)
conjgrad$O:	conjgrad$(F7)
gmresr$O:	gmresr$(F7)
tridiag$O:	tridiag$(F7)


########### Extra dependencies for some files with "INCLUDE:" statements

vacphys$F: vacphys.t vacphys.mhd0.t vacphys.mhdres.t 
	$(VACPP) $< $(PREFOR) > $@

vacgrid$F: vacgrid.t
	$(VACPP) $< $(PREFOR) > $@

################# Makefile targets #########################################

# The "vacini" and "vac" targets use the default Fortran compiler
vacini	:$(VACINIOBJ) $(MAKE_)
	$(FOR) $(FORFLG) -o $(VACDIR)/vacini $(VACINIOBJ)

vac	:$(VACOBJ) $(LIBS_) $(MAKE_)
	$(FOR) $(FORFLG) -o $(VACDIR)/vac $(VACOBJ) $(LIBS)

# The "vacini90" and "vac90" targets always use the Fortran 90 compiler
vacini90:
	@make PREFOR='$(F90TOF90)' FOR='$(F90)' FORFLG='$(F90FLAGS)' vacini

vac90:
	@make PREFOR='$(F90TOF90)' FOR='$(F90)' FORFLG='$(F90FLAGS)' vac

# The "vacini77" and "vac77" targets always use the Fortran 77 compiler
vacini77:
	@make PREFOR='$(F90TOF77)' FOR='$(F77)' FORFLG='$(F77FLAGS)' vacini

vac77:
	@make PREFOR='$(F90TOF77)' FOR='$(F77)' FORFLG='$(F77FLAGS)' vac

# The "vacinihpf" and "vachpf" targets always use the HPF compiler
# and vacpp.pl -O optimization for shifts
vacinihpf: 
	@make PREFOR='$(F90TOHPF)' FOR='$(HPF)' FORFLG='$(HPFFLAGS)' vacini

vachpf: 
	@make VACPP='vacpp.pl -O' \
		PREFOR='$(F90TOHPF)' FOR='$(HPF)' FORFLG='$(HPFFLAGS)' vac

########### Targets using single source file, useful for debugging

vaciniall$F: $(VACINIFOR)
	cat $(VACINIFOR) > vaciniall$F

vaciniall:vaciniall$O $(MAKE_) 
	$(FOR) $(FORFLG) -o $(VACDIR)/vaciniall vaciniall$O

vacall$F: $(VACFOR)
	cat $(VACFOR) > vacall$F

vacall	:vacall$O $(LIBS_) $(MAKE_) 
	$(FOR) $(FORFLG) -o $(VACDIR)/vacall vacall$O $(LIBS)

# "vacsmall" uses less memory than "vac" by overlapping
# the local arrays in subroutines using COMMON blocks.
# However it takes some time to compile, and the VACUSR modules may not
# be compatible with the requirements of "savemem" (see Man/srclimits.html)

VACALLF90:
	@if test $(FOR) = $(F77); then				\
	   rm -f vacall$F $(INCLUDES);				\
	   make PREFOR='$(F90TOF90)' vacall$F $(INCLUDES);	\
	fi

vacsmall$F: VACALLF90 vacall$F $(INCLUDES)
	@if test $(FOR) = $(F77); then				\
	    echo 'savemem vacall$F > vacsmallF90';		\
	    savemem vacall$F > vacsmallF90;			\
	    echo 'cat vacsmallF90 $(F90TOF77) > vacsmall$F';	\
	    cat vacsmallF90 $(F90TOF77) > vacsmall$F;		\
	    echo 'rm vacsmallF90 $(INCLUDES)';			\
	    rm vacsmallF90 $(INCLUDES);				\
	    make $(INCLUDES);					\
	else							\
	    echo 'savemem vacall$F > vacsmall$F';		\
	    savemem vacall$F > vacsmall$F;			\
	fi

vacsmall: vacsmall$O $(LIBS_) $(MAKE_) 
	$(FOR) $(FORFLG) -o $(VACDIR)/vacsmall vacsmall$O $(LIBS)

########### Convert binary formats on a DEC Alpha for other machines

ieee2dec: converti.f90 converto.f90 Makefile
	$(F90) -u -convert big_endian    -c converti.f90
	$(F90) -u -convert little_endian -c converto.f90
	$(F90) -o $(VACDIR)/ieee2dec converti.o converto.o

dec2ieee: converti.f90 converto.f90 Makefile
	$(F90) -u -convert little_endian -c converti.f90
	$(F90) -u -convert big_endian    -c converto.f90
	$(F90) -o $(VACDIR)/dec2ieee converti.o converto.o

########### Roe solver can be tested by roetest

roetest	:$(ROEOBJ) 
	$(FOR) $(FORFLG) -o $(VACDIR)/roetest $(ROEOBJ)

###### Removing object files, precompiled Fortran files, and symbolic links

clean	:
	rm -f vac*$F roetest$F Vvac*$F vac*.d *$O *~

cleanall: clean
	rm -f vacphys.t vacphys0.t vacpar.t vacusrpar.t vacusr.t 
