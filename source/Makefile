################################################################################
#
#	Makefile to build Cython and Fortran modules into compiled libraries
#
################################################################################

all: 
	python setup.py build_ext --inplace

clean:
	python setup.py clean
	rm -rf rmg/spectral/*.so rmg/unirxn/*.so rmg/*.so rmg/*.c build/
	rm -rf rmg/*.pyc  rmg/system/*.pyc rmg/spectral/*.pyc rmg/unirxn/*.pyc