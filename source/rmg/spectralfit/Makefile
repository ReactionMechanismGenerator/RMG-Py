################################################################################
#
#	Makefile for compiling spectralfit Fortran module
#
################################################################################

F90=g95

LDFLAGS=-lblas -llapack

################################################################################

all: spectralfit.so 

spectralfit.so: spectralfit.f90 calc_freq_code.o dqed.o
	f2py --fcompiler=$(F90) calc_freq_code.o dqed.o -c spectralfit.f90 -m spectralfit $(LDFLAGS)
	cp spectralfit.so ../
	
calc_freq_code.o: calc_freq_code.f90
	$(F90) -c calc_freq_code.f90 -o calc_freq_code.o

dqed.o: dqed.f90
	$(F90) -c dqed.f90 -o dqed.o

clean:
	rm -f *.o *.mod spectralfit.so ../spectralfit.so 

