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

LDFLAGS=-lblas -llapack

################################################################################

all: cython spectral unirxn rmg/spectralfit.so
	python setup.py build_ext --inplace

cython: rmg/graph.py rmg/chem.py rmg/thermo.py 

spectral: rmg/spectralfit.so

unirxn: rmg/unirxn/states.so rmg/unirxn/mastereqn.so rmg/unirxn/msc.so rmg/unirxn/rs.so

rmg/unirxn/%.so: rmg/unirxn/%.f90
	make -C rmg/unirxn F90=$(F90) F2PY_FCOMPILER=$(F2PY_FCOMPILER) LDFLAGS=$(LDFLAGS)

rmg/spectralfit.so: rmg/spectralfit/*.f90
	make -C rmg/spectralfit F90=$(F90) F2PY_FCOMPILER=$(F2PY_FCOMPILER) LDFLAGS=$(LDFLAGS)

clean:
	make -C rmg/spectralfit clean
	make -C rmg/unirxn clean
	rm -rf rmg/*.so rmg/*.c build/
	