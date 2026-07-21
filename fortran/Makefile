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

#Will detect ifort/gfortran/flang or edit for your compiler
ifneq ($(COMPILER),gfortran)
ifneq ($(COMPILER),flang)
ifortErr = $(shell which ifort >/dev/null 2>&1; echo $$?)
else
ifortErr = 1
endif
else
ifortErr = 1
endif

ifeq "$(ifortErr)" "0"

#Intel compiler
F90C     = ifort

ifortVer_major = $(shell ifort -v 2>&1 | cut -d " " -f 3 | cut -d. -f 1)
compiler_ver = $(shell ifort -v 2>&1)
#other tests will be done by forutils
ifeq ($(shell test $(ifortVer_major) -lt 18; echo $$?),0)
$(warning ** ifort version lower than 18.0.1 may not work with the Python wrapper **)
endif


ifeq ($(shell test $(ifortVer_major) -gt 15; echo $$?),0)
COMMON_FFLAGS = -fpp -qopenmp
else
COMMON_FFLAGS = -fpp -openmp
endif
COMMON_FFLAGS += -gen-dep=$$*.d

FFLAGS = -fp-model precise -W0 -WB $(COMMON_FFLAGS)
DEBUGFLAGS =  -g -check all -check noarg_temp_created -traceback -fpe0 $(COMMON_FFLAGS)

ifeq ($(shell uname -s),Darwin)
SFFLAGS = -dynamiclib #-fpic
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
ifeq ($(COMPILER),flang)
gfortErr = 1
endif
ifeq "$(gfortErr)" "0"
#Gfortran compiler (version 6+):
compiler_ver = $(shell gfortran -dumpversion 2>&1)
COMPILER = gfortran
F90C     = gfortran
COMMON_FFLAGS = -MMD -cpp -ffree-line-length-none -fmax-errors=4 -fopenmp
# Using -ffast-math causes differences between Debug and Release configurations.
FFLAGS = -O3 $(COMMON_FFLAGS)
DEBUGFLAGS = -g -fbacktrace -ffpe-trap=invalid,overflow,zero -fbounds-check $(COMMON_FFLAGS)
ifeq ($(shell uname -s),Darwin)
SFFLAGS = -dynamiclib #-fpic
ifneq ($(shell echo $(compiler_ver) | awk '{print ($$1 >= 14.0)}'),0)
LIBGFORTRAN_A := $(shell $(F90C) -print-file-name=libgfortran.a)
ifneq ($(shell [ -f "$(LIBGFORTRAN_A)" ] && echo yes),)
SFFLAGS += -static-libgfortran -static-libquadmath -static-libgcc
endif
endif
else
SFFLAGS = -shared -fpic
endif



MODOUT =  -J$(OUTPUT_DIR)
SMODOUT = -J$(DLL_DIR)
FISHER_LINK = -lblas -llapack

ifeq ($(shell uname -s),Darwin)
NONNATIVE = 1
endif
ifndef NONNATIVE
#Note this seems to make code slightly slower in some cases, use CLUSTER_SAFE=1 to test without
FFLAGS+=-march=native
endif

else
#LLVM flang (experimental): used if explicitly requested with COMPILER=flang,
#or as a last-resort fallback when neither ifort nor gfortran is found.
#Set F90C to override which flang executable is used.
FLANGCC := $(shell command -v flang-new 2>/dev/null || command -v flang 2>/dev/null || \
	ls /usr/bin/flang-new-[0-9]* /usr/bin/flang-[0-9]* /usr/lib/llvm-*/bin/flang 2>/dev/null | sort -V | tail -n1)
ifeq ($(FLANGCC),)
$(error No supported Fortran compiler found (checked ifort, gfortran, flang). Install one, or set COMPILER/F90C)
endif

compiler_ver = flang-$(shell $(FLANGCC) --version 2>&1 | head -n1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -n1)
COMPILER = flang
F90C     ?= $(FLANGCC)
#SFFLAGS/MODOUT are only meaningful at link time, but (as with the other compilers) get passed to every
#per-file compile too; flang's clang-based driver warns about that, gfortran/ifort don't, hence -Qunused-arguments.
COMMON_FFLAGS = -cpp -ffree-form -fopenmp -fPIC -Qunused-arguments
FFLAGS = -O3 $(COMMON_FFLAGS)
DEBUGFLAGS = -g -O0 $(COMMON_FFLAGS)

ifeq ($(shell uname -s),Darwin)
SFFLAGS = -dynamiclib -fPIC
else
SFFLAGS = -shared -fPIC
endif

MODOUT  = -module-dir $(OUTPUT_DIR)
SMODOUT = -module-dir $(DLL_DIR)
FISHER_LINK = -lblas -llapack

ifeq ($(shell uname -s),Darwin)
NONNATIVE = 1
endif
ifndef NONNATIVE
FFLAGS+=-march=native
endif
endif
endif

IFLAG = -I

ifneq ($(FISHER),)
# Its dependencies are all meet by the libutils.a which always added.
FFLAGS += -DFISHER
endif

DEBUGFLAGS ?= $(FFLAGS)
Debug: FFLAGS = $(DEBUGFLAGS)

include ./Makefile_main
