************
Installation
************

Dependencies
============

CanTherm relies on a number of Python packages for certain functionality. The
following lists the required Python packages that are not part of the Python
standard library:

* `NumPy <http://www.numpy.org>`_. Provides efficient array and matrix 
  operations.

* `SciPy <http://www.scipy.org>`_. Provides efficient linear algebra functions
  and special functions.

CanTherm also depends on the Green group's base library for Python projects, 
named `ChemPy <http://github.com/jwallen/ChemPy>`_.

Installing CanTherm
===================

Once you have obtained the required dependencies, CanTherm can be installed by
executing the following command from within the root package directory::

    $ python setup.py install

This will move all of the necessary files to the appropriate Python installation
location on disk. You can also use CanTherm without installing if you prefer. 
This is done by appending the full path to the package's root directory to your 
``PYTHONPATH`` environment variable. 
