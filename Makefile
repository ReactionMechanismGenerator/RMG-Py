################################################################################
#
#   Makefile for RMG Py
#
################################################################################

DASPK=$(shell python -c 'import pydas.daspk; print pydas.daspk.__file__')
DASSL=$(shell python -c 'import pydas.dassl; print pydas.dassl.__file__')
RDKIT_VERSION=$(shell python -c 'import rdkit; print rdkit.__version__')

.PHONY : all minimal main solver cantherm clean decython documentation QM mopac_travis

all: main solver QM

noQM: main solver

minimal:
	python setup.py build_ext minimal --build-lib . --build-temp build --pyrex-c-in-temp

main:
	@ echo "Checking you have PyDQED..."
	@ python -c 'import pydqed; print pydqed.__file__'
	python setup.py build_ext main --build-lib . --build-temp build --pyrex-c-in-temp

solver:

ifneq ($(DASPK),)
	@ echo "DASPK solver found. Compiling with DASPK and sensitivity analysis capability..."
	@ echo "DEF DASPK = 1" > rmgpy/solver/settings.pxi 
else ifneq ($(DASSL),)
	@ echo "DASSL solver found. Compiling with DASSL.  Sensitivity analysis capabilities are off..."
	@ echo "DEF DASPK = 0" > rmgpy/solver/settings.pxi
else
	@ echo 'No PyDAS solvers found.  Please check if you have the latest version of PyDAS.'
	@ python -c 'import pydas.dassl' 
endif
	python setup.py build_ext solver --build-lib . --build-temp build --pyrex-c-in-temp

cantherm:
	python setup.py build_ext cantherm --build-lib . --build-temp build --pyrex-c-in-temp

bin/symmetry:
	mkdir -p bin
	$(MAKE) -C external/symmetry install

QM: bin/symmetry
	@ echo "Checking you have rdkit..."
	@ python -c 'import rdkit; print rdkit.__file__'
	@ echo "Checking rdkit version..."
ifneq ($(RDKIT_VERSION),)
	@ echo "Found rdkit version $(RDKIT_VERSION)"
else
	$(error RDKit version out of date, please install RDKit version 2015.03.1 or later with InChI support);
endif
	@ echo "Checking rdkit has InChI support..."
	@ python -c 'from rdkit import Chem; assert Chem.inchi.INCHI_AVAILABLE, "RDKit installed without InChI Support. Please install with InChI."'

documentation:
	$(MAKE) -C documentation html
	@ echo "Start at: documentation/build/html/index.html"

clean:
	python setup.py clean --build-temp build
	rm -rf build/
	find . -name '*.so' -exec rm -f '{}' \;
	find . -name '*.pyc' -exec rm -f '{}' \;
	$(MAKE) -C external/symmetry clean
	rm -f bin/symmetry

clean-solver:
	rm -r build/pyrex/rmgpy/solver/
	rm -r build/build/pyrex/rmgpy/solver/
	find rmgpy/solver/ -name '*.so' -exec rm -f '{}' \;
	find rmgpy/solver/ -name '*.pyc' -exec rm -f '{}' \;

decython:
	# de-cythonize all but the 'minimal'. Helpful for debugging in "pure python" mode.
	find . -name *.so ! \( -name _statmech.so -o -name quantity.so -o -regex '.*rmgpy/solver/.*' \) -exec rm -f '{}' \;
	find . -name *.pyc -exec rm -f '{}' \;

test:
	mkdir -p testing/coverage
	rm -rf testing/coverage/*
	nosetests --nocapture --nologcapture --all-modules --verbose --with-coverage --cover-inclusive --cover-package=rmgpy --cover-erase --cover-html --cover-html-dir=testing/coverage --exe rmgpy

test-database:
	nosetests -v -d testing/databaseTest.py	

eg1: noQM
	mkdir -p testing/minimal
	rm -rf testing/minimal/*
	cp examples/rmg/minimal/input.py testing/minimal/input.py
	coverage erase
	@ echo "Running minimal example with coverage tracking AND profiling"
	coverage run rmg.py -p testing/minimal/input.py
	coverage report
	coverage html
eg2: all
	mkdir -p testing/hexadiene
	rm -rf testing/hexadiene/*
	cp examples/rmg/1,3-hexadiene/input.py testing/hexadiene/input.py
	coverage erase
	@ echo "Running 1,3-hexadiene example with coverage tracking AND profiling"
	coverage run rmg.py -p testing/hexadiene/input.py
	coverage report
	coverage html

eg3: all
	mkdir -p testing/liquid_phase
	rm -rf testing/liquid_phase/*
	cp examples/rmg/liquid_phase/input.py testing/liquid_phase/input.py
	coverage erase
	@ echo "Running liquid_phase example with coverage tracking AND profiling"
	coverage run rmg.py -p testing/liquid_phase/input.py
	coverage report
	coverage html

eg5: all
	mkdir -p testing/heptane-eg5
	rm -rf testing/heptane-eg5/*
	cp examples/rmg/heptane-eg5/input.py testing/heptane-eg5/input.py
	@ echo "Running heptane-eg5 example."
	python rmg.py -q testing/heptane-eg5/input.py

eg6: all
	mkdir -p testing/ethane-oxidation
	rm -rf testing/ethane-oxidation/*
	cp examples/rmg/ethane-oxidation/input.py testing/ethane-oxidation/input.py
	@ echo "Running ethane-oxidation example."
	python rmg.py -q testing/ethane-oxidation/input.py

eg7: all
	mkdir -p testing/gri_mech_rxn_lib
	rm -rf testing/gri_mech_rxn_lib/*
	cp examples/rmg/gri_mech_rxn_lib/input.py testing/gri_mech_rxn_lib/input.py
	@ echo "Running gri_mech_rxn_lib example."
	python rmg.py testing/gri_mech_rxn_lib/input.py
	
######### 
# Section for setting up MOPAC calculations on the Travis-CI.org server
ifeq ($(TRAVIS),true)
ifneq ($(TRAVIS_SECURE_ENV_VARS),true)
SKIP_MOPAC=true
endif
endif
mopac_travis:
ifeq ($(TRAVIS),true)
ifneq ($(TRAVIS_SECURE_ENV_VARS),true)
	@echo "Don't have MOPAC licence key on this Travis build so can't test QM"
else
	@echo "Installing MOPAC key"
	@yes Yes | mopac $(MOPACKEY)
endif
else
	@#echo "Not in Travis build, no need to run this target"
endif
# End of MOPAC / TRAVIS stuff
#######

eg4: all mopac_travis
ifeq ($(SKIP_MOPAC),true)
	@echo "Skipping eg4 (without failing) because can't run MOPAC"
else
	mkdir -p testing/thermoEstimator
	rm -rf testing/thermoEstimator/*
	cp examples/thermoEstimator/input.py testing/thermoEstimator/input.py
	@ echo "Running thermo data estimator example. This tests QM."
	python scripts/thermoEstimator.py testing/thermoEstimator/input.py
endif
