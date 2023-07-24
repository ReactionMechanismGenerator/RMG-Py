################################################################################
#
#   Makefile for RMG Py
#
################################################################################

.PHONY : all minimal main solver check pycheck arkane clean install decython documentation test

all: pycheck main solver check

minimal:
	python setup.py build_ext minimal --inplace --build-temp .

main:
	python setup.py build_ext main --inplace --build-temp .

solver:
	@ python utilities.py check-pydas
	python setup.py build_ext solver --inplace --build-temp .

arkane:
	python setup.py build_ext arkane --inplace --build-temp .

check:
	@ python utilities.py check-dependencies

pycheck:
	@ python utilities.py check-python

documentation:
	$(MAKE) -C documentation html
	@ echo "Start at: documentation/build/html/index.html"

clean:
	@ python utilities.py clean

clean-solver:
	@ python utilities.py clean-solver

install:
	@ python utilities.py check-pydas
	python setup.py install
	
decython:
	# de-cythonize all but the 'minimal'. Helpful for debugging in "pure python" mode.
	find . -name *.so ! \( -name _statmech.so -o -name quantity.so -o -regex '.*rmgpy/solver/.*' \) -exec rm -f '{}' \;
	find . -name *.pyc -exec rm -f '{}' \;

test-all:
	python-jl -m pytest

test test-unittests:
	python-jl -m pytest -m "not functional not database"

test-functional:
	python-jl -m pytest -m "functional"

test-database:
	python-jl -m pytest -m "database"

eg0: all
	mkdir -p testing/eg0
	rm -rf testing/eg0/*
	cp examples/rmg/superminimal/input.py testing/eg0/input.py
	@ echo "Running eg0: superminimal (H2 oxidation) example"
	python rmg.py testing/eg0/input.py

eg1: all
	mkdir -p testing/eg1
	rm -rf testing/eg1/*
	cp examples/rmg/minimal/input.py testing/eg1/input.py
	coverage erase
	@ echo "Running eg1: minimal (ethane pyrolysis) example with coverage tracking AND profiling"
	coverage run rmg.py -p testing/eg1/input.py
	coverage report
	coverage html

eg2: all
	mkdir -p testing/eg2
	rm -rf testing/eg2/*
	cp examples/rmg/1,3-hexadiene/input.py testing/eg2/input.py
	coverage erase
	@ echo "Running eg2: 1,3-hexadiene example with profiling"
	python rmg.py -p testing/eg2/input.py

eg3: all
	mkdir -p testing/eg3
	rm -rf testing/eg3/*
	cp examples/rmg/liquid_phase/input.py testing/eg3/input.py
	coverage erase
	@ echo "Running eg3: liquid_phase example with profiling"
	python rmg.py -p testing/eg3/input.py

eg5: all
	mkdir -p testing/eg5
	rm -rf testing/eg5/*
	cp examples/rmg/heptane-eg5/input.py testing/eg5/input.py
	@ echo "Running eg5: heptane example"
	python rmg.py testing/eg5/input.py

eg6: all
	mkdir -p testing/eg6
	rm -rf testing/eg6/*
	cp examples/rmg/ethane-oxidation/input.py testing/eg6/input.py
	@ echo "Running eg6: ethane-oxidation example"
	python rmg.py testing/eg6/input.py

eg7: all
	mkdir -p testing/eg7
	rm -rf testing/eg7/*
	cp examples/rmg/gri_mech_rxn_lib/input.py testing/eg7/input.py
	@ echo "Running eg7: gri_mech_rxn_lib example"
	python rmg.py testing/eg7/input.py
	
scoop: all
	mkdir -p testing/scoop
	rm -rf testing/scoop/*
	cp examples/rmg/minimal/input.py testing/scoop/input.py
	coverage erase
	@ echo "Running minimal example with SCOOP"
	python -m scoop -n 2 rmg.py -v testing/scoop/input.py

eg4: all
	mkdir -p testing/eg4
	rm -rf testing/eg4/*
	cp examples/thermoEstimator/input.py testing/eg4/input.py
	@ echo "Running thermo data estimator example. This tests QM."
	python scripts/thermoEstimator.py testing/eg4/input.py
