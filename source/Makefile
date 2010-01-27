################################################################################
#
#	Makefile to build Cython and Fortran modules into compiled libraries
#
################################################################################

#F90=gfortran
#F2PY_FCOMPILER=gnu95
# Alternatively:
F90=g95
F2PY_FCOMPILER=g95
F90FLAGS=-fPIC -r8 -ftrace=full -fbounds-check -Wall

################################################################################

all: cython spectral unirxn
	python setup.py build_ext --inplace

cython: rmg/graph.py rmg/chem.py rmg/thermo.py 

spectral: rmg/spectral/fit.so

unirxn: rmg/unirxn/states.so rmg/unirxn/mastereqn.so rmg/unirxn/msc.so rmg/unirxn/rs.so

rmg/unirxn/%.so: rmg/unirxn/%.f90
	make -C rmg/unirxn F90=$(F90) F2PY_FCOMPILER=$(F2PY_FCOMPILER) F90FLAGS="$(F90FLAGS)"

rmg/spectral/%.so: rmg/spectral/%.f90
	make -C rmg/spectral F90=$(F90) F2PY_FCOMPILER=$(F2PY_FCOMPILER) F90FLAGS="$(F90FLAGS)" LDFLAGS=$(LDFLAGS)

clean:
	make -C rmg/spectral clean
	make -C rmg/unirxn clean
	rm -rf rmg/*.so rmg/*.c build/
	
