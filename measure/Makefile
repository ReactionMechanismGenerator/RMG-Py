################################################################################
#
#   Makefile for MEASURE
#
################################################################################

CYTHON_FLAGS=--build-lib . --build-temp build --pyrex-c-in-temp

CLEAN_FLAGS=--build-temp build

-include make.inc

all: cython

install:
	python setup.py install

cython:
	python setup.py build_ext $(CYTHON_FLAGS)

clean:
	python setup.py clean $(CLEAN_FLAGS)

cleanall: clean
	rm -f measure/*.so measure/*.pyc
