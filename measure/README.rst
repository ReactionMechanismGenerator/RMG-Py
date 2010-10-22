*********************************************************************
MEASURE - Master Equation Automatic Solver for Unimolecular REactions
*********************************************************************

**MEASURE** is a free, open-source Python package for modeling of pressure-
dependent unimolecular reaction networks using a master equation approach.

Installation
============

MEASURE depends on several other packages in order to provide its full
functional capabilities. Below is a summary of those dependencies. Of these,
Python and NumPy are required, while the rest are optional depending on
what functionality you wish to use. All dependencies are free and open
source:

* `Python <http://www.python.org/>`_ (versions 2.5.x and 2.6.x are known to work)

* `NumPy <http://numpy.scipy.org/>`_ (version 1.3.0 or later is recommended)

* `Cython <http://www.cython.org/>`_ (version 0.12.1 or later is recommended)

* C compiler

MEASURE builds on the ChemPy package.

Compilation with Cython
-----------------------

Working with the master equation for anything but the smallest networks can be
time consuming. For this reason, MEASURE has been designed to be compiled into
C-like extensions using Cython. This compilation is not required, but is 
recommended due to the significant speed boost that comes with it. To compile 
the extensions, invoke the following from within the MEASURE root directory::

    $ python setup.py build_ext --inplace

You can also use the provided Makefile to do this::

    $ make

You can set extra options by creating a file in the MEASURE root directory
named make.inc.
