#CAMB sources makefile
#Edit for your compiler

#Intel ifort , -openmp toggles mutli-processor:
F90C     = ifort
FFLAGS = -openmp -ip -O2 -vec_report0 -W0 -WB -fpp2

#Sun, single processor:
#F90C     = f90
#FFLAGS = -O2

# Intel 9 on IA-64 (eg. COSMOS)
# (do "module load icomp90" before compiling)
#F90C = ifort
#FFLAGS = -openmp -fpp2 -w -O3 -ip -mP2OPT_hlo_prefetch=F

#Intel ifc, add -openmp for multi-processor (some have bugs):
#F90C     = ifc
#FFLAGS = -O2 -Vaxlib -ip -W0 -WB -quiet -fpp2
#some systems can can also add e.g. -tpp7 -xW

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
#FFLAGS = -DNAGF95 -O3 -mismatch

#Gfortran compiler
#F90C = gfc
#FFLAGS = -O2 

#G95 compiler
#F90C   = g95
#FFLAGS = -O2

#Sun V880
#F90C = mpf90
#FFLAGS =  -O4 -openmp -ftrap=%none -dalign -DMPI

#Sun parallel enterprise:
#F90C     = f95
#FFLAGS =  -O2 -xarch=native64 -openmp -ftrap=%none
#try removing -openmp if get bus errors. -03, -04 etc are dodgy.

#IBM XL Fortran, multi-processor (run gmake)
#F90C     = xlf90_r
#FFLAGS  = -qsmp=omp -qsuffix=f=f90:cpp=F90 -O3 -qstrict -qarch=pwr3 -qtune=pwr3


#Files containing evolution equations initial power spectrum module
EQUATIONS     = equations
POWERSPECTRUM = power_tilt
#Module doing non-linear scaling
NONLINEAR     = halofit
#Instead use perturbation theory (good for high redshift/21cm approx)
#NONLINEAR     = nonlinear_PT
REIONIZATION = reionization
RECOMBINATION = recfast

#Driver program
DRIVER        = inidriver.F90

CAMBLIB       = libcamb_sources.a

#Shouldn't need to change anything else...

F90FLAGS      = $(FFLAGS)
FC            = $(F90C)

CAMBOBJ       = constants.o utils.o subroutines.o inifile.o $(POWERSPECTRUM).o $(RECOMBINATION).o $(REIONIZATION).o modules.o \
	bessels.o $(EQUATIONS).o $(NONLINEAR).o lensing.o cmbmain.o camb.o



default: camb

all: camb $(CAMBLIB)


subroutines.o: constants.o utils.o
$(POWERSPECTRUM).o: subroutines.o  inifile.o
$(RECOMBINATION).o: subroutines.o inifile.o
$(REIONIZATION).o: constants.o inifile.o
modules.o: $(REIONIZATION).o $(POWERSPECTRUM).o $(RECOMBINATION).o
bessels.o: modules.o
$(EQUATIONS).o: bessels.o
$(NONLINEAR).o:  modules.o
lensing.o: bessels.o
cmbmain.o: lensing.o $(NONLINEAR).o $(EQUATIONS).o
camb.o: cmbmain.o


camb: $(CAMBOBJ) $(DRIVER)
	$(F90C) $(F90FLAGS) $(CAMBOBJ) $(DRIVER) -o $@

$(CAMBLIB): $(CAMBOBJ)
	ar -r $@ $?

%.o: %.f90
	$(F90C) $(F90FLAGS) -c $*.f90


%.o: %.F90
	$(F90C) $(F90FLAGS) -c $*.F90

clean:
	-rm -f *.o *.a *.d core *.mod


