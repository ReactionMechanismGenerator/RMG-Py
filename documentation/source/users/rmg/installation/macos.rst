.. _macos:

********************
MacOS X Installation
********************

.. warning::
    This installation method is no longer maintained, and is unlikely to work as written.
    Please refer to :ref:`installation` for more up-to-date instructions.

There are a number of dependencies for RMG-Py. This page will guide you through installing them.
You will need the Command Line Tools for XCode. If you are not using Anaconda to install RMG-Py,
we highly recommend the `Homebrew <https://brew.sh>`_ package manager.
The following instructions assume that you have `installed Homebrew and its requirements <https://brew.sh>`_.
We recommend using a `Virtual Environment <https://docs.python-guide.org/dev/virtualenvs/>`_ for your Python packages,
but this is optional (without it you may need to add `sudo` before some commands to solve permission errors).

You will also need gfortran, Python, Numpy and Scipy. We typically install them using 
homebrew-python (which used to be at ``https://github.com/Homebrew/homebrew-python``)  but other methods may work as well.

* For example::

	brew tap homebrew/python
	brew install numpy
	brew install scipy
	brew install matplotlib --with-cairo --with-ghostscript --with-ticl-tk --with-pyqt --with-pygtk --withgtk3

* Install git if you don't already have it (you may also like some graphical interfaces like `GitX <https://github.com/gitx/gitx>`_ or `GitUp <https://gitup.co>`_ or dozens more)::

	brew update
	brew install git

* Optional (but recommended for Nitrogen-chemistry nomenclature): install `OpenBabel <http://openbabel.org>`_::

	brew install open-babel --with-python --HEAD

* Install `RDKit <https://www.rdkit.org>`_::

	brew tap rdkit/rdkit
	brew install rdkit --with-inchi
	brew link --overwrite rdkit

  You'll need to set an environment variable to use it, eg. put this in your `~/.bash_profile` file::
	
	export RDBASE=/usr/local/share/RDKit

* The following dependencies are also required for core RMG functions and must be installed from source before building RMG:

  **pyrdl:** RingDecomposerLib, used for ring perception. Download from https://github.com/rareylab/RingDecomposerLib. Requires CMAKE to compile.

  **lpsolve:** Mixed integer linear programming solver. Download from https://sourceforge.net/projects/lpsolve/. Python extension also required.

* Make a directory to put everything in::

	mkdir ~/Code

* Get the RMG-Py source code and the RMG-database from GitHub::

	cd ~/Code
	git clone https://github.com/ReactionMechanismGenerator/RMG-database.git
	git clone https://github.com/ReactionMechanismGenerator/RMG-Py.git

* Install the Python dependencies listed in the :file:`RMG-Py/requirements.txt` file using `pip` (do ``easy_install pip`` if you don't already have it)::

	pip install -r RMG-Py/requirements.txt

* Get and build `PyDQED <https://github.com/ReactionMechanismGenerator/PyDQED>`_::

	cd ~/Code
	git clone https://github.com/ReactionMechanismGenerator/PyDQED.git
	cd PyDQED
	export LIBRARY_PATH=$(dirname $(gfortran -print-file-name=libgfortran.a))
	make
	make install

* Get and build `PyDAS <https://github.com/ReactionMechanismGenerator/PyDAS>`_::

	cd ~/Code
	git clone https://github.com/ReactionMechanismGenerator/PyDAS.git
	cd PyDAS
	export LIBRARY_PATH=$(dirname $(gfortran -print-file-name=libgfortran.a))
	make
	make install

* Build RMG-Py::

	cd ~/Code/RMG-Py
	make -j4

* Run an example: ::

	cd ~/Code/RMG-Py/
	python rmg.py examples/rmg/minimal/input.py

  Verify your installation by opening the resulting output.html file under the "examples/rmg/minimal" directory.
  
  You can also use the Makefile targets to test and run examples: ::

	cd ~/Code/RMG-Py/
	make test
	make eg1
	make eg2

To run with on-the-fly Quantum Mechanics calculations, you will also need to install
`MOPAC <https://openmopac.net/downloads.html>`_ or `Gaussian <https://gaussian.com>`_, then run `make QM`.
