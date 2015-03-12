.. _macos:

********************
MacOS X Installation
********************

There are a number of dependencies for RMG-Py. This page will guide you through installing them.
You will need the Command Line Tools for XCode. We highly recommend the `Homebrew <http://mxcl.github.com/homebrew/>`_ package manager.
The following instructions assume that you have `installed Homebrew and its requirements <https://github.com/mxcl/homebrew/wiki/installation>`_.

You will also need gfortran, Python, Numpy and Scipy. We typically install them following `these instructions
<http://www.thisisthegreenroom.com/2011/installing-python-numpy-scipy-matplotlib-and-ipython-on-lion/>`_
but other methods like the `Scipy Superpack <http://fonnesbeck.github.com/ScipySuperpack/>`_ may work as well.

* Install git (you may also like some graphical interfaces like `mxcl's GitX <https://github.com/mxcl/gitx/downloads>`_ or `GitHub for Mac <http://mac.github.com/>`_)::

	brew update
	brew install git

* Make a directory to put everything in::

	mkdir ~/Code

* Get the RMG-Py source code from GitHub::

	cd ~/Code
	git clone git@github.com:GreenGroup/RMG-database.git
	git clone git@github.com:GreenGroup/RMG-Py.git

* Install the Python dependencies listed in the :file:`RMG-Py/requirements.txt` file using `pip` (do ``easy_install pip`` if you don't already have it)::

	pip install -r RMG-Py/requirements.txt

* Get and build `PyDQED <https://github.com/GreenGroup/PyDQED>`_::

	cd ~/Code
	git clone https://github.com/GreenGroup/PyDQED.git
	cd PyDQED
	make

  It will build the fortran files, but fail to link using clang, so force it to use gcc for the linking with::

	LIBRARY_PATH=$(dirname $(gfortran -print-file-name=libgfortran.a)) LDSHARED='gcc -bundle -undefined dynamic_lookup -arch x86_64' python setup.py install

* Get and build `PyDAS <https://github.com/GreenGroup/PyDAS>`_::

	cd ~/Code
	git clone https://github.com/GreenGroup/PyDAS.git
	cd PyDAS
	make

  It will build the fortran files, but fail to link using clang, so force it to use gcc for the linking with::

	LIBRARY_PATH=$(dirname $(gfortran -print-file-name=libgfortran.a)) LDSHARED='gcc -bundle -undefined dynamic_lookup -arch x86_64' python setup.py install

* Get the RMG-Database and get and build RMG-Py::

	cd ~/Code
	git clone git@github.com:GreenGroup/RMG-database.git
	git clone git@github.com:GreenGroup/RMG-Py.git
	cd RMG-Py
	make

* Install RDKit::

	brew uninstall boost
	brew install boost --build-from-source
	brew tap edc/homebrew-rdkit
	brew install rdkit --with-inchi
	
  You'll need various environment variables set, eg.::
	
	export RDBASE=$HOME/rdkit # CHECK THIS (maybe you put RDKit somewhere else)
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$RDBASE/lib
	export PYTHONPATH=$PYTHONPATH:$RDBASE

* Run an example: ::

	cd ~/Code/RMG-Py/
	python rmg.py examples/rmg/minimal/input.py

  Verify your installation by opening the resulting output.html file under the "examples/rmg/minimal" directory.
  
  You can also use the Makefile targets to test and run examples: ::

	make test
	make eg1
	make eg2
