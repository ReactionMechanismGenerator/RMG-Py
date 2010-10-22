################################################################################
#
#	Makefile to build Cython and Fortran modules into compiled libraries
#
################################################################################

-include make.inc

all: 
	python setup.py build_ext $(CYTHON_FLAGS)
	echo "To remove the cython files without removing the compiled fortran: make clean-cython"

clean:
	python setup.py clean $(CLEAN_FLAGS)
	rm -rf rmg/spectral/*.so rmg/unirxn/*.so rmg/*.so rmg/*.c build
	rm -rf rmg/*.pyc rmg/thermo/*.pyc rmg/kinetics/*.pyc rmg/system/*.pyc rmg/spectral/*.pyc rmg/unirxn/*.pyc

clean-cython:
	rm -rf rmg/*.so rmg/*.c rmg/thermo/*.so rmg/thermo/*.c
