################################################################################
#
#   Makefile for RMG Py
#
################################################################################

.PHONY : all minimal main solver check pycheck arkane clean install decython documentation mopac_travis

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

q2dtor:
	@ echo -e "\nInstalling Q2DTor...\n"
	@ echo -e "Q2DTor is a software for calculating the partition functions and themodynamic properties\
	of molecular systems with two or more torsional modes developed by David Ferro Costas (david.ferro@usc.es)\
	 and Antonio Fernandez Ramos (qf.ramos@usc.es) at the Universidade de Santiago de Compostela. Arkane can\
	  integrate Q2DTor to compute the quantum mechanical partition function of 2D rotors.  \n\nFor use of Q2DTor\
 and HinderedRotor2D within Arkane please cite:  \n\nD. Ferro-Costas, M. N. D. S.Cordeiro, D. G. Truhlar, A.\
		  Fernández-Ramos, Comput. Phys. Commun. 232, 190-205, 2018.\n"
	@ read -p "Press ENTER to continue" dummy
	@ git clone https://github.com/mjohnson541/Q2DTor.git external/Q2DTor --branch arkanepy3
	
decython:
	# de-cythonize all but the 'minimal'. Helpful for debugging in "pure python" mode.
	find . -name *.so ! \( -name _statmech.so -o -name quantity.so -o -regex '.*rmgpy/solver/.*' \) -exec rm -f '{}' \;
	find . -name *.pyc -exec rm -f '{}' \;

test-all:
ifneq ($(OS),Windows_NT)
	mkdir -p testing/coverage
	rm -rf testing/coverage/*
endif
	nosetests --nocapture --nologcapture --all-modules --verbose --with-coverage --cover-inclusive --cover-erase --cover-html --cover-html-dir=testing/coverage --exe rmgpy arkane

test test-unittests:
ifneq ($(OS),Windows_NT)
	mkdir -p testing/coverage
	rm -rf testing/coverage/*
endif
	nosetests --nocapture --nologcapture --all-modules -A 'not functional' --verbose --with-coverage --cover-inclusive --cover-erase --cover-html --cover-html-dir=testing/coverage --exe rmgpy arkane

test-functional:
ifneq ($(OS),Windows_NT)
	mkdir -p testing/coverage
	rm -rf testing/coverage/*
endif
	nosetests --nologcapture --all-modules -A 'functional' --verbose --exe rmgpy arkane

test-database:
	nosetests --nocapture --nologcapture --verbose --detailed-errors testing/databaseTest.py

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
	@ echo "Running eg2: 1,3-hexadiene example with coverage tracking AND profiling"
	coverage run rmg.py -p testing/eg2/input.py
	coverage report
	coverage html

eg3: all
	mkdir -p testing/eg3
	rm -rf testing/eg3/*
	cp examples/rmg/liquid_phase/input.py testing/eg3/input.py
	coverage erase
	@ echo "Running eg3: liquid_phase example with coverage tracking AND profiling"
	coverage run rmg.py -p testing/eg3/input.py
	coverage report
	coverage html

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
	mkdir -p testing/eg4
	rm -rf testing/eg4/*
	cp examples/thermoEstimator/input.py testing/eg4/input.py
	@ echo "Running thermo data estimator example. This tests QM."
	python scripts/thermoEstimator.py testing/eg4/input.py
endif
