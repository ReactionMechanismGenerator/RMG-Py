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

* Install `Open Babel <http://openbabel.org/>`_ and its Python bindings. For now, the easiest way to get version 2.3.1 and the Python bindings is::

	brew install https://raw.github.com/rwest/homebrew/open-babel-new/Library/Formula/eigen2.rb
	brew install https://raw.github.com/rwest/homebrew/open-babel-new/Library/Formula/open-babel.rb


* Make a directory to put everything in::

	mkdir ~/Code

* Get the RMG-Py source code from GitHub::

	cd ~/Code
	git clone git@github.com:GreenGroup/RMG-database.git

* Install the Python dependencies listed in the :file:`RMG-Py/requirements.txt` file using `pip` (do ``easy_install pip`` if you don't already have it)::

	pip install -r RMG-Py/requirements.txt

* Get and build `PyDQED <https://github.com/jwallen/PyDQED>`_::

	cd ~/Code
	git clone https://github.com/jwallen/PyDQED.git
	cd PyDQED
	make

  It will build the fortran files, but fail to link using clang, so force it to use gcc for the linking with::

	LIBRARY_PATH=/usr/local/lib/gcc LDSHARED='gcc -bundle -undefined dynamic_lookup -arch x86_64' python setup.py install

* Get and build `PyDAS <https://github.com/jwallen/PyDAS>`_::

	cd ~/Code
	git clone https://github.com/jwallen/PyDAS.git
	cd PyDQED
	make

  It will build the fortran files, but fail to link using clang, so force it to use gcc for the linking with::

	LIBRARY_PATH=/usr/local/lib/gcc LDSHARED='gcc -bundle -undefined dynamic_lookup -arch x86_64' python setup.py install

* Get the RMG-Database and get and build RMG-Py::

	cd ~/Code
	git clone git@github.com:GreenGroup/RMG-database.git
	git clone git@github.com:GreenGroup/RMG-Py.git
	cd RMG-Py
	make

* Install RDKit, *if* (and only if) you want to use the new 3D geometry branch::

	brew uninstall boost
	brew install boost --build-from-source
	brew tap rwest/homebrew-rdkit
	brew install rdkit --with-inchi

* Run an example: ::

	cd ~/Code/RMG-Py/
	python rmg.py examples/rmg/minimal/input.py

  Verify your installation by opening the resulting output.html file under the "examples/rmg/minimal" directory.
