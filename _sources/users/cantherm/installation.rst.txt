************
Installation
************

Dependencies
============

Python versions 2.6+ are recommended for CanTherm.  Currently Python 3.x versions
are incompatible with CanTherm.

CanTherm relies on a number of Python packages for certain functionality. The
following lists the required Python packages that are not part of the Python
standard library:

* `NumPy <http://www.numpy.org>`_. Provides efficient array and matrix 
  operations.

* `SciPy <http://www.scipy.org>`_. Provides efficient linear algebra functions
  and special functions.
  
* `Cython v0.19+ <http://www.cython.org>`_. C-extensions for Python

* `RDKit <http://www.rdkit.org>`_. Cheminformatics libraries and functions

* `Cairo <http://cairographics.org/>`_. Cairo graphics rendering for Python for drawing reaction networks

* Quantities. For converting between different scientific units

* Argparse. For parsing input arguments when running scripts

You can install these dependencies on Linux in the following fashion. ::

	sudo apt-get install python-numpy python-scipy python-cairo
	sudo pip install cython>=0.19 quantities argparse

It is recommended that RDKit be installed manually with InChI capabilities on.


Installing CanTherm
===================

Once you have obtained the required dependencies, CanTherm can be installed by 
first obtaining the RMG-Py source code by either downloading into a directory using git::

	git clone git@github.com:ReactionMechanismGenerator/RMG-Py.git
	
or by downloading the zip file of the current RMG-Py master source code found 
`here <https://github.com/ReactionMechanismGenerator/RMG-Py/archive/master.zip>`_ and unzipping into
the appropriate directory.

Inside the root package directory for ``RMG-Py``, execute the following make command::

	make cantherm
	
The appropriate cythonization and compilation steps will now begin.  When they are completed,
you can run a test example by going into the examples folder to run a sample Cantherm job 
to verify the installation::

	cd examples/cantherm/networks/acetyl+O2
	python ../../../../cantherm.py input.py

This will allow you to use the cantherm.py as a Python script anytime you point to it. 
