################################################################################
#
#   Makefile for RMG Py
#
################################################################################

.PHONY : all minimal main measure solver clean

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
