#CAMB Makefile

#Set FISHER=Y to compile bispectrum fisher matrix code
FISHER=

#Will detect ifort/gfortran or edit for your compiler
ifortErr = $(shell which ifort >/dev/null 2>&1; echo $$?)
ifeq "$(ifortErr)" "0"

#Intel compiler
F90C     = ifort
ifortVer_major = $(shell ifort -v 2>&1 | cut -d " " -f 3 | cut -d. -f 1)
FFLAGS = -openmp -fast -W0 -WB -fpp 
#FFLAGS = -openmp -fast -W0 -WB -fpp2 -vec_report0
DEBUGFLAGS =-openmp -g -check all -check noarg_temp_created -traceback -fpp -fpe0
# Activate dependency generation by the compiler
F90DEPFLAGS = -gen-dep=$$*.d
## This is flag is passed to the Fortran compiler allowing it to link C++ if required (not usually):
F90CRLINK = -cxxlib
ifeq "$(ifortVer_major)" "15"
F90DRLINK += -qopt-report=1 -qopt-report-phase=vec
else
ifeq "$(ifortVer_major)" "14"
F90CRLINK += -vec_report0
else
error "Unsupported version of ifort"
endif
endif
MODOUT = -module $(OUTPUT_DIR)
AR     = xiar
ifneq ($(FISHER),)
FFLAGS += -mkl
endif

else
gfortErr = $(shell which gfortran >/dev/null; echo $$?)
ifeq "$(gfortErr)" "0"

#Gfortran compiler:
#The options here work in v4.6+
F90C     = gfortran
COMMON_FFLAGS = -cpp -ffree-line-length-none -fmax-errors=4
FFLAGS = -O3 -fopenmp -ffast-math $(COMMON_FFLAGS)
DEBUGFLAGS = -g -fbacktrace -ffpe-trap=invalid,overflow,zero -fbounds-check $(COMMON_FFLAGS)
# Activate dependency generation by the compiler
F90DEPFLAGS = -MMD
MODOUT =  -J$(OUTPUT_DIR)

ifneq ($(FISHER),)
F90CRLINK += -lblas -llapack
endif
ifneq ($(shell uname -s),Darwin)
#native optimization does not work on Mac
O3FLAGS += -march=native
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
