################################################################################
#
#	Makefile for compiling spectral Fortran modules
#
################################################################################

all: _modes.so _fit.so

fit.so: fit.f90 calc_freq_code.o dqed.o
	f2py --fcompiler=$(F2PY_FCOMPILER) --f90flags="$(F90FLAGS)" calc_freq_code.o dqed.o -c fit.f90 -m fit $(LDFLAGS)

_modes.so: math.o _modes.f90
	f2py --fcompiler=$(F2PY_FCOMPILER) --f90flags="$(F90FLAGS)" math.o -c _modes.f90 -m _modes $(LDFLAGS)

_fit.so: _fit.o math.o _modes.o cases.o dqed.o
	f2py --fcompiler=$(F2PY_FCOMPILER) --f90flags="$(F90FLAGS)" math.o _modes.o cases.o dqed.o -c _fit.f90 -m _fit $(LDFLAGS)

%.o: %.f90
	$(F90) $(F90FLAGS) -c $< -o $@

clean:
	rm -f *.o *.mod *.so

