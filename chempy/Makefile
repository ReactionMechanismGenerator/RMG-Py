################################################################################
#
#   Makefile for ChemPy
#
################################################################################

CYTHON_FLAGS=--build-lib . --build-temp build --pyrex-c-in-temp

CLEAN_FLAGS=--build-temp build

-include make.inc

all: cython

cython:
	python setup.py build_ext $(CYTHON_FLAGS)

clean:
	python setup.py clean $(CLEAN_FLAGS)

cleanall: clean
	rm -f chempy/*.so chempy/*.pyc chempy/ext/*.so chempy/ext/*.pyc
