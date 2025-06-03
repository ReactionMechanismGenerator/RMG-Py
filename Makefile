################################################################################
#
#   Makefile for RMG Py
#
################################################################################

CC=gcc
CXX=g++

.PHONY : all check clean install decython documentation test q2dtor

all: check install check

check:
	@ python utilities.py check-dependencies
	@ python utilities.py check-pydas

documentation:
	$(MAKE) -C documentation html
	@ echo "Start at: documentation/build/html/index.html"

clean:
	@ python utilities.py clean
	python -m pip uninstall --yes reactionmechanismgenerator || true  # can fail if RMG not installed at all

clean-solver:
	@ python utilities.py clean-solver

install:
	@ python utilities.py check-pydas
	python -m pip install -vv -e .

q2dtor:
	@ echo -e "\nInstalling Q2DTor...\n"
	@ echo -e "Q2DTor is a software for calculating the partition functions and themodynamic properties\
	of molecular systems with two or more torsional modes developed by David Ferro Costas (david.ferro@usc.es)\
	 and Antonio Fernandez Ramos (qf.ramos@usc.es) at the Universidade de Santiago de Compostela. Arkane can\
	  integrate Q2DTor to compute the quantum mechanical partition function of 2D rotors.  \n\nFor use of Q2DTor\
 and HinderedRotor2D within Arkane please cite:  \n\nD. Ferro-Costas, M. N. D. S.Cordeiro, D. G. Truhlar, A.\
		  Fernández-Ramos, Comput. Phys. Commun. 232, 190-205, 2018.\n"
	@ read -p "Press ENTER to continue" dummy
	@ mkdir -p external
	@ git clone https://github.com/cathedralpkg/Q2DTor external/Q2DTor

decython:
	# de-cythonize all but the 'minimal'. Helpful for debugging in "pure python" mode.
	find . -name *.so ! \( -name _statmech.so -o -name quantity.so -o -regex '.*rmgpy/solver/.*' \) -exec rm -f '{}' \;
	find . -name *.pyc -exec rm -f '{}' \;

test-all:
	python -m pytest

test test-unittests:
	python -m pytest -m "not functional and not database"

test-functional:
	python -m pytest -m "functional"

test-database:
	python -m pytest -m "database"

eg0: all
	mkdir -p testing/eg0
	rm -rf testing/eg0/*
	cp examples/rmg/superminimal/input.py testing/eg0/input.py
	@ echo "Running eg0: superminimal (H2 oxidation) example"
	RMG_DISABLE_JULIA=1 python rmg.py testing/eg0/input.py

eg1: all
	mkdir -p testing/eg1
	rm -rf testing/eg1/*
	cp examples/rmg/minimal/input.py testing/eg1/input.py
	coverage erase
	@ echo "Running eg1: minimal (ethane pyrolysis) example with coverage tracking AND profiling"
	RMG_DISABLE_JULIA=1 coverage run rmg.py -p testing/eg1/input.py
	coverage report
	coverage html

eg2: all
	mkdir -p testing/eg2
	rm -rf testing/eg2/*
	cp examples/rmg/1,3-hexadiene/input.py testing/eg2/input.py
	coverage erase
	@ echo "Running eg2: 1,3-hexadiene example with profiling"
	RMG_DISABLE_JULIA=1 python rmg.py -p testing/eg2/input.py

eg3: all
	mkdir -p testing/eg3
	rm -rf testing/eg3/*
	cp examples/rmg/liquid_phase/input.py testing/eg3/input.py
	coverage erase
	@ echo "Running eg3: liquid_phase example with profiling"
	RMG_DISABLE_JULIA=1 python rmg.py -p testing/eg3/input.py

eg5: all
	mkdir -p testing/eg5
	rm -rf testing/eg5/*
	cp examples/rmg/heptane-eg5/input.py testing/eg5/input.py
	@ echo "Running eg5: heptane example"
	RMG_DISABLE_JULIA=1 python rmg.py testing/eg5/input.py

eg6: all
	mkdir -p testing/eg6
	rm -rf testing/eg6/*
	cp examples/rmg/ethane-oxidation/input.py testing/eg6/input.py
	@ echo "Running eg6: ethane-oxidation example"
	RMG_DISABLE_JULIA=1 python rmg.py testing/eg6/input.py

eg7: all
	mkdir -p testing/eg7
	rm -rf testing/eg7/*
	cp examples/rmg/gri_mech_rxn_lib/input.py testing/eg7/input.py
	@ echo "Running eg7: gri_mech_rxn_lib example"
	RMG_DISABLE_JULIA=1 python rmg.py testing/eg7/input.py
	
scoop: all
	mkdir -p testing/scoop
	rm -rf testing/scoop/*
	cp examples/rmg/minimal/input.py testing/scoop/input.py
	coverage erase
	@ echo "Running minimal example with SCOOP"
	RMG_DISABLE_JULIA=1 python -m scoop -n 2 rmg.py -v testing/scoop/input.py

eg4: all
	mkdir -p testing/eg4
	rm -rf testing/eg4/*
	cp examples/thermoEstimator/input.py testing/eg4/input.py
	@ echo "Running thermo data estimator example. This tests QM."
	python scripts/thermoEstimator.py testing/eg4/input.py

# RMS reactor examples (require Julia)
eg8: all
	mkdir -p testing/eg8
	rm -rf testing/eg8/*
	cp examples/rmg/rms_constant_V/input.py testing/eg8/input.py
	@ echo "Running RMS constantVIdealGasReactor example (requires Julia)"
	python rmg.py testing/eg8/input.py

eg9: all
	mkdir -p testing/eg9
	rm -rf testing/eg9/*
	cp examples/rmg/nox_transitory_edge/input.py testing/eg9/input.py
	@ echo "Running RMS constantTPIdealGasReactor example (requires Julia)"
	python rmg.py testing/eg9/input.py

eg10: all
	mkdir -p testing/eg10
	rm -rf testing/eg10/*
	cp examples/rmg/liquid_cat/input.py testing/eg10/input.py
	@ echo "Running RMS liquidSurfaceReactor example (requires Julia)"
	python rmg.py testing/eg10/input.py
