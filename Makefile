################################################################################
#
#	Makefile to build Cython and Fortran modules into compiled libraries
#
################################################################################

all: 
	python setup.py build_ext --inplace
	echo "To remove the cython files without removing the compiled fortran: make clean-cython"

clean:
	python setup.py clean
	rm -rf rmg/spectral/*.so rmg/unirxn/*.so rmg/*.so rmg/*.c build
	rm -rf rmg/*.pyc rmg/thermo/*.pyc rmg/kinetics/*.pyc rmg/system/*.pyc rmg/spectral/*.pyc rmg/unirxn/*.pyc

clean-cython:
	rm -rf rmg/*.so rmg/*.c rmg/thermo/*.so rmg/thermo/*.c
