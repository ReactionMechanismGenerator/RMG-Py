################################################################################
#
#	Makefile to build Cython and Fortran modules into compiled libraries
#
################################################################################

all: 
	python setup.py build_ext --inplace

clean:
	python setup.py clean
	rm -rf *.so *.pyc build
