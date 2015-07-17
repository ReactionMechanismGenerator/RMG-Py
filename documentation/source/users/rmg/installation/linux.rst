.. _linux:

******************
Linux Installation
******************

RMG-Py and all of its dependencies may be easily installed through a short series of Terminal commands.
The instructions listed below have been confirmed on a fresh Ubuntu 12.04 installation and should generally apply to other distributions.

* Install compilers and libraries: ::

	sudo apt-get install git g++ gfortran python-dev liblapack-dev
	sudo apt-get install python-openbabel python-setuptools python-pip

* After creating a `Github account <https://github.com/signup/free>`_, generate your public key: ::

	cd ~; ssh-keygen		# press enter to save to the default directory
					# create a password if desired
	cat .ssh/id_rsa.pub

  Copy this public key to your `Github profile <https://github.com/settings/ssh>`_.

* Install dependencies: ::

	sudo apt-get install libpng-dev libfreetype6-dev graphviz mencoder
	
	sudo pip install numpy		# install NumPy before other packages
	
	sudo pip install scipy cython nose matplotlib quantities guppy sphinx psutil xlwt 
	
	cd ~
	git clone https://github.com/ReactionMechanismGenerator/PyDAS.git
	git clone https://github.com/ReactionMechanismGenerator/PyDQED.git
	cd PyDAS; make F77=gfortran; sudo make install; cd ..
	cd PyDQED; make F77=gfortran; sudo make install; cd ..

* Install RDKit

  Full installation instructions: http://code.google.com/p/rdkit/wiki/GettingStarted
  Be sure to **build it with InChI support.** Here's a synopsis: ::
  
	cd ~
	sudo apt-get install flex bison build-essential python-numpy cmake python-dev sqlite3
	sudo apt-get install libsqlite3-dev libboost-dev libboost-python-dev libboost-regex-dev
	git clone https://github.com/rdkit/rdkit.git
	cd rdkit
	export RDBASE=`pwd`
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$RDBASE/lib
	export PYTHONPATH=$PYTHONPATH:$RDBASE
	cd External/INCHI-API
	./download-inchi.sh
	cd ../../
	mkdir build
	cd build
	cmake .. -DRDK_BUILD_INCHI_SUPPORT=ON
	make
	make install
	
  You'll need various environment variables set (you may want to add these to your `.bash_profile` file), eg.::
  
  	export RDBASE=$HOME/rdkit # CHECK THIS (maybe you put RDKit somewhere else)
  	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$RDBASE/lib
  	export PYTHONPATH=$PYTHONPATH:$RDBASE  # (or some other way to make sure it's on your Python path)

* Install RMG-Py: ::

	cd ~
	git clone https://github.com/ReactionMechanismGenerator/RMG-database.git
	git clone https://github.com/ReactionMechanismGenerator/RMG-Py.git
	sudo pip install -r RMG-Py/requirements.txt
	cd RMG-Py
	make

* Run an example: ::

	python rmg.py examples/rmg/minimal/input.py

  Verify your installation by opening the resulting output.html file under the "examples/rmg/minimal" directory.

  You can also use the Makefile targets to test and run examples: ::
  
	make test
	make eg1
	make eg2

.. _compile_sensitivity:

* Compiling RMG-Py with Sensitivity Analysis: 

  Running sensitivity analysis in RMG-Py requires the prerequisite DASPK solver and DASPK compiled wrapper in PyDAS.  
  To do so first compile daspk in PyDAS and agree to download the daspk31.tgz file when prompted. ::

    cd PyDAS/
    make
    make install

  Then compile RMG-Py normally.  It will automatically be compiled with sensitivity analysis if DASPK is found. ::

    cd RMG-Py
    make clean-solver
    make
    
  Note that using this option will allow RMG to both run with and without sensitivity.  
