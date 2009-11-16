################################################################################
#
#	Makefile to build Cython and Fortran modules into compiled libraries
#
################################################################################

all: rmg/graph.so rmg/chem.so rmg/thermo.so rmg/spectralfit.so
	python setup.py build_ext --inplace

rmg/graph.so: rmg/graph.c

rmg/graph.c: rmg/graph.py rmg/graph.pxd rmg/chem.pxd
	cython -a rmg/graph.py

rmg/chem.so: rmg/chem.c

rmg/chem.c: rmg/chem.py rmg/chem.pxd
	cython -a rmg/chem.py

rmg/thermo.so: rmg/thermo.c 

rmg/thermo.c: rmg/thermo.py rmg/thermo.pxd
	cython -a rmg/thermo.py

rmg/spectralfit.so: rmg/spectralfit/spectralfit.f90 rmg/spectralfit/calc_freq_code.f90 rmg/spectralfit/dqed.f90
	make -C rmg/spectralfit

clean:
	rm -rf rmg/*.so rmg/*.c build/