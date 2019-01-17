#CAMB Makefile

#Set FISHER=Y to compile bispectrum fisher matrix code
FISHER=

#Set FORUTILSPATH to the path where the libforutils.a file can be found.
#The OUTPUT_DIR will be appended.

ifneq "$(wildcard ../forutils)" ""
FORUTILSPATH ?= $(shell pwd)/../forutils
else
ifneq "$(wildcard ../../forutils)" ""
FORUTILSPATH ?= $(shell pwd)/../../forutils
else
ifneq "$(wildcard ../../../forutils)" ""
FORUTILSPATH ?= $(shell pwd)/../../../forutils
endif
endif
endif

ifeq ($(FORUTILSPATH),)
$(error Use  git clone --recurse-submodules, install forutils from https://github.com/cmbant/forutils, or set FORUTILSPATH variable)
endif

#native optimization does not work on Mac gfortran or heterogeneous clusters
CLUSTER_SAFE ?= 0
ifneq ($(CLUSTER_SAFE), 0)
NONNATIVE = 1
endif

#Will detect ifort/gfortran or edit for your compiler
ifneq ($(COMPILER),gfortran)
ifortErr = $(shell which ifort >/dev/null 2>&1; echo $$?)
else
ifortErr = 1
endif

ifeq "$(ifortErr)" "0"

#Intel compiler
F90C     = ifort

ifortVer_major = $(shell ifort -v 2>&1 | cut -d " " -f 3 | cut -d. -f 1)
ifeq ($(shell test $(ifortVer_major) -gt 15; echo $$?),0)
COMMON_FFLAGS = -fpp -qopenmp
else
COMMON_FFLAGS = -fpp -openmp
endif
COMMON_FFLAGS += -gen-dep=$$*.d

FFLAGS = -fp-model precise -W0 -WB $(COMMON_FFLAGS)
DEBUGFLAGS =  -g -check all -check noarg_temp_created -traceback -fpe0 $(COMMON_FFLAGS)

ifeq ($(shell uname -s),Darwin)
SFFLAGS = -dynamiclib -fpic
else
SFFLAGS = -shared -fpic
endif

ifdef NONNATIVE
FFLAGS+=-O3 -ipo -axCORE-AVX2
else
FFLAGS+=-fast
endif

## This is flag is passed to the Fortran compiler allowing it to link C++ if required (not usually):
F90CRLINK = -cxxlib
ifneq "$(ifortVer_major)" "14"
F90CRLINK += -qopt-report=0 -qopt-report-phase=vec
else
F90CRLINK += -vec_report0
endif
MODOUT = -module $(OUTPUT_DIR)
AR     = xiar
SMODOUT = -module $(DLL_DIR)
ifneq ($(FISHER),)
FFLAGS += -mkl
endif

else
gfortErr = $(shell which gfortran >/dev/null; echo $$?)
ifeq "$(gfortErr)" "0"
#Gfortran compiler (version 6+):
COMPILER = gfortran
F90C     = gfortran
COMMON_FFLAGS = -MMD -cpp -ffree-line-length-none -fmax-errors=4
# Using -ffast-math causes differences between Debug and Release configurations.
FFLAGS = -O3 -fopenmp $(COMMON_FFLAGS)
DEBUGFLAGS = -g -fbacktrace -ffpe-trap=invalid,overflow,zero -fbounds-check $(COMMON_FFLAGS)
SFFLAGS =  -shared -fPIC
MODOUT =  -J$(OUTPUT_DIR)
SMODOUT = -J$(DLL_DIR)

ifneq ($(FISHER),)
F90CRLINK += -lblas -llapack
endif

ifeq ($(shell uname -s),Darwin)
NONNATIVE = 1
endif
ifndef NONNATIVE
#Note this seems to make code slightly slower in some cases, use CLUSTER_SAFE=1 to test without
FFLAGS+=-march=native
endif
endif
endif

IFLAG = -I

#G95 compiler
#F90C   = g95
#FFLAGS = -O2

#SGI, -mp toggles multi-processor. Use -O2 if -Ofast gives problems.
#F90C     = f90
#FFLAGS  = -Ofast -mp

#Digital/Compaq fortran, -omp toggles multi-processor
#F90C    = f90
#FFLAGS  = -omp -O4 -arch host -math_library fast -tune host -fpe1

#Absoft ProFortran, single processor:
#F90C     = f95
#FFLAGS = -O2 -cpu:athlon -s -lU77 -w -YEXT_NAMES="LCS" -YEXT_SFX="_"

#NAGF95, single processor:
#F90C     = f95
#FFLAGS = -DNAGF95 -O3

#PGF90
#F90C = pgf90
#FFLAGS = -O2 -DESCAPEBACKSLASH -Mpreprocess

#Sun V880
#F90C = mpf90
#FFLAGS =  -O4 -openmp -ftrap=%none -dalign

#Sun parallel enterprise:
#F90C     = f95
#FFLAGS =  -O2 -xarch=native64 -openmp -ftrap=%none
#try removing -openmp if get bus errors. -03, -04 etc are dodgy.

#IBM XL Fortran, multi-processor (run gmake)
#F90C     = xlf90_r
#FFLAGS  = -DESCAPEBACKSLASH -DIBMXL -qsmp=omp -qsuffix=f=f90:cpp=F90 -O3 -qstrict -qarch=pwr3 -qtune=pwr3

#Settings for building camb_fits
#Location of FITSIO and name of library
FITSDIR       ?= /usr/local/lib
FITSLIB       = cfitsio
#Location of HEALPIX for building camb_fits
HEALPIXDIR    ?= /usr/local/healpix

ifneq ($(FISHER),)
# Its dependencies are all meet by the libutils.a which always added.
FFLAGS += -DFISHER
endif

DEBUGFLAGS ?= $(FFLAGS)
Debug: FFLAGS = $(DEBUGFLAGS)

include ./Makefile_main
