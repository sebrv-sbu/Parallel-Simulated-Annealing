# CODE VERSION ######################################################

VERSION = 3.3

# executables to make

#TSPEXECS  = tsp_sa tsplibconvert tspdistance calc_ave_error_bar curve_fit 
TSPEXECS  = initialize tsp_sa.mpi calc_ave_error_bar curve_fit printscore
# Removing tsplibconvert - SEB RV AUG 5 2024
# Adding prinstcore and initialize - SEB RV AUG 8 2024
# FLAGS FOR -v FOR ALL EXECUTABLES ##################################
# this passes user and host name, compiler and version to the com-
# pile so it can be printed in the -v message of each executable

USRFLAG  = -DUSR=\"$(USER)\"
HOSTFLAG = -DMACHINE=\"$(HOST)\"
COMPFLAG = -DCOMPILER=\"gcc\"
FLAGFLAG = -DFLAGS=\"optimized\"
VERSFLAG = -DVERS=\"$(VERSION)\"

VFLAGS = $(USRFLAG) $(HOSTFLAG) $(COMPFLAG) $(FLAGFLAG) $(VERSFLAG) 

# find out about which architecture we're on and set compiler 
# accordingly

ifeq ($(OSTYPE),osf1)
	CC = cc
endif

HOST := $(shell uname -n)

MPIFLAGS = $(CCFLAGS) -DMPI
DEBUGFLAGS = $(DEBUGFLAGS) -DMPI
PROFILEFLAGS = $(PROFILEFLAGS) -DMPI
MPICC = mpicc
#TSPEXECS = tsp_sa tsp_sa.mpi tspdistance tsplibconvert calc_ave_error_bar curve_fit
TSPEXECS = tsp_sa.mpi calc_ave_error_bar curve_fit printscore
# Removing tsplibconvert - SEB RV AUG 5 2024
# adding printscore and initialize - SEB RV AUG 8 2024
#nompi	TSPEXECS = tsp_sa tsplibconvert tspdistance calc_ave_error_bar curve_fit

# find the compiler and set corresponding flags
CC = mpicc
MPI = on
#ifeq ($(CC),cc)
#	CCFLAGS = -std1 -fast -DALPHA_DU -DNDEBUG 
# 	DEBUGFLAGS = -std1 -O0 -g -ieee
# 	PROFILEFLAGS = -std1 -O2 -g1 -pg
#	LIBS = -lm  -ldxml
#	FLIBS = $(LIBS)
#	KCC = /bin/kcc
#	KFLAGS = -ckapargs=' -scalaropt=3 -optimize=5 -roundoff=3 -arl=4 '
# uncomment 2 lines below if you don't want kcc
#	KCC = $(CC)
#	KFLAGS = $(CCFLAGS)
#endif


#ifeq ($(CC),icc)
	#CCFLAGS = -O1 -DNDEBUG 
#	CCFLAGS = -O3 -xW -tpp7 -ipo 
#	PRECFLAGS = -mp -prec_div -pc80 
#	DEBUGFLAGS = -g  -inline_debug_info -O3  -xW -tpp7
#	PROFILEFLAGS = -prof_dir profile_data -prof_gen -O2  -xW -tpp7
#	USEPROFFLAGS = -prof_use  -prof_dir profile_data  -O3 -xW -tpp7 -ipo -opt_report
#	LIBS = -limf 
#	FLIBS = -lm -static (already commented out)
#	KCC = $(CC)
#	KFLAGS = $(CCFLAGS)
#	export ICC = "yes"
#endif

#ifeq ($(CC),gcc)
CCFLAGS = -O3 -DNDEBUG -fcommon 
DEBUGFLAGS = -Wall -g -fcommon -fsanitize=address -fno-omit-frame-pointer
PROFILEFLAGS = -g -pg -O3
#	PROFILEFLAGS = -g -pg
DOSTATIC=no
ifeq ($(DOSTATIC),yes)
	LIBS = -lm -static
else
	LIBS = -lm
endif
#	FLIBS = -lm -static
KCC = $(CC)
KFLAGS = $(CCFLAGS)
#endif

# debugging?
#DEBUG=on
ifdef DEBUG
	CCFLAGS = $(DEBUGFLAGS) -fcommon
	FLAGFLAG = -DFLAGS=\"debugging\"
else
	DEBUG = "Off"
endif

ifdef PROFILE
	CCFLAGS = $(PROFILEFLAGS)
	FLAGFLAG = -DFLAGS=\"profiling\"
	KCC = $(CC)
	KFLAGS =
else
	PROFILE = "Off"
endif
ifdef USEPROFILE
	CCFLAGS = $(USEPROFFLAGS)
endif

# export all variables that Makefiles in subdirs need
OMPI_CC=gcc #this may be needed 
export INCLUDES = -I. -I../lam -I/usr/local/include
export CFLAGS = $(CCFLAGS) $(INCLUDES)
export VFLAGS
export CC
export KCC
export MPICC
export KFLAGS
export LIBS
export MPIFLAGS
export TSPEXECS
export MPI
export OMPI_CC
#define targets

tsp: lsa
	@cd tsp && $(MAKE)

deps: 
	cd tsp && $(MAKE) -f basic.mk Makefile

lsa:
	@cd lam && $(MAKE)

clean:
	rm -f core* *.o
	rm -f */core* */*.o

veryclean:
	rm -f */core* */*.o 
	rm -f tsp/tsp_sa tsp/tsp_sa.mpi
	rm -f tsp/tsplibconvert
	rm -f tsp/tspdistance
	rm -f tsp/calc_ave_error_bar
	rm -f tsp/curve_fit
	rm -f tsp/printscore
	rm -f lam/gen_deviates
	rm -f tsp/Makefile

help:
	@echo "make: this is the Makefile for fly and/or TSP LSA code"
	@echo "      always 'make deps' first after a 'make veryclean'"
	@echo ""
	@echo "      the following targets are available:"
	@echo "      lsa:       make object files in the lam directory only"
	@echo "      tsp:       compile the TSP code (which is in 'tsp')"
	@echo "      clean:     gets rid of cores and object files"
	@echo "      veryclean: gets rid of executables and dependencies too"
	@echo ""
	@echo "      your current settings are:"   
	@echo "      compiler:  $(CC)"
	@echo "      flags:     $(CFLAGS)"
	@echo "      libs:      $(LIBS)"
	@echo "      debugging: $(DEBUG) (1=yes/empty=no)"
	@echo "      profiling: $(PROFILE)"
	@echo ""

