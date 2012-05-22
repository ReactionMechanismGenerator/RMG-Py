################################################################################
#
#   Makefile for RMG Py
#
################################################################################

.PHONY : all minimal main measure solver clean decython

all: main measure solver

minimal: measure solver

main:
	python setup.py build_ext main --build-lib . --build-temp build --pyrex-c-in-temp

measure:
	python setup.py build_ext measure --build-lib . --build-temp build --pyrex-c-in-temp

solver:
	python setup.py build_ext solver --build-lib . --build-temp build --pyrex-c-in-temp

clean:
	python setup.py clean --build-temp build
	rm -rf build/
	find . -name *.so -exec rm -f '{}' \;
	find . -name *.pyc -exec rm -f '{}' \;

decython:
	# de-cythonize all but the 'minimal'. Helpful for debugging in "pure python" mode.
	find . -name *.so ! \( -name _statmech.so -o -name quantity.so -o -regex '.*rmgpy/measure/.*' -o -regex '.*rmgpy/solver/.*' \) -exec rm -f '{}' \;
	find . -name *.pyc -exec rm -f '{}' \;
