################################################################################
#
#   Makefile for ChemPy
#
################################################################################

-include make.inc

all: cython

cython:
	python setup.py build_ext $(CYTHON_FLAGS)

clean:
	python setup.py clean $(CLEAN_FLAGS)

cleanall: clean
	rm -f chempy/*.so chempy/*.pyc chempy/ext/*.so chempy/ext/*.pyc
