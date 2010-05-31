***************************************
ChemPy - A chemistry toolkit for Python
***************************************

**ChemPy** is a free, open-source Python toolkit for chemistry, chemical
engineering, and materials science applications.

Installation
============

Below is a summary of the dependencies needed to install and use ChemPy:

* Python (versions 2.5.x and 2.6.x are known to work)

* NumPy (version 1.3.0 or later is recommended)

* Cython (version 0.12.1 or later is recommended)

* C and Fortran compilers

All of the ChemPy modules have been designed to be compiled into C extensions
using Cython. This compilation is not required, but is strongly recommended due
to the enormous speed boost that comes with it. To compile the extensions,
invoke the following from within the ChemPy root directory::

    $ python setup.py build_ext --inplace

You can also use the provided Makefile to do this::

    $ make

You can set extra options by creating a file in the ChemPy root directory
named make.inc.
