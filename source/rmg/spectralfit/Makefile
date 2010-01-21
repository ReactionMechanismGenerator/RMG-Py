################################################################################
#
#	Makefile for compiling spectralfit Fortran module
#
################################################################################

F90FLAGS=-fPIC -r8

all: spectralfit.so 

spectralfit.so: spectralfit.f90 calc_freq_code.o dqed.o
	f2py --fcompiler=$(F2PY_FCOMPILER) --f90flags="$(F90FLAGS)" calc_freq_code.o dqed.o -c spectralfit.f90 -m spectralfit $(LDFLAGS)
	mv spectralfit.so ../
	
calc_freq_code.o: calc_freq_code.f90
	$(F90) $(F90FLAGS) -c calc_freq_code.f90 -o calc_freq_code.o

dqed.o: dqed.f90
	$(F90) $(F90FLAGS) -c dqed.f90 -o dqed.o

clean:
	rm -f *.o *.mod ../spectralfit.so 

