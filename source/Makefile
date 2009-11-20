################################################################################
#
#	Makefile to build Cython and Fortran modules into compiled libraries
#
################################################################################

all: rmg/graph.py rmg/chem.py rmg/thermo.py rmg/spectralfit.so
	python setup.py build_ext --inplace

rmg/spectralfit.so: rmg/spectralfit/spectralfit.f90 rmg/spectralfit/calc_freq_code.f90 rmg/spectralfit/dqed.f90
	make -C rmg/spectralfit

clean:
	rm -rf rmg/*.so rmg/*.c build/