################################################################################
#
#   Makefile for RMG Py
#
################################################################################

all: chem measure solver statesfit

minimal: measure solver

chem:
	python setup.py build_ext chem --build-lib . --build-temp build/chem --pyrex-c-in-temp

measure:
	python setup.py build_ext measure --build-lib . --build-temp build/measure --pyrex-c-in-temp

solver:
	python setup.py build_ext solver --build-lib . --build-temp build/solver --pyrex-c-in-temp

statesfit:
	make -C rmgpy/statesfit

clean:
	python setup.py clean --build-temp build
	make -C rmgpy/statesfit clean
	rm -rf build/
	find . -name *.so -exec rm -f '{}' \;
	find . -name *.pyc -exec rm -f '{}' \;
