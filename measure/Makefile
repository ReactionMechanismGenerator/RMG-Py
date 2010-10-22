################################################################################
#
#   Makefile for MEASURE
#
################################################################################

-include make.inc

all: cython

install:
	python setup.py install

cython:
	python setup.py build_ext --inplace $(CYTHON_FLAGS)

clean:
	python setup.py clean $(CLEAN_FLAGS)

cleanall: clean
	rm -f measure/*.so measure/*.pyc
