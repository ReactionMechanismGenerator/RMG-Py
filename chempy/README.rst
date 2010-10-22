***************************************
ChemPy - A chemistry toolkit for Python
***************************************

**ChemPy** is a free, open-source Python toolkit for chemistry, chemical
engineering, and materials science applications.

Installation
============

ChemPy depends on several other packages in order to provide its full
functional capabilities. Below is a summary of those dependencies. Of these,
Python and NumPy are required, while the rest are optional depending on
what functionality you wish to use. All dependencies are free and open
source:

* `Python <http://www.python.org/>`_ (versions 2.5.x and 2.6.x are known to work)

* `NumPy <http://numpy.scipy.org/>`_ (version 1.3.0 or later is recommended)

* `SciPy <http://www.scipy.org/>`_ (version 0.7.0 or later is recommended)

* `Cython <http://www.cython.org/>`_ (version 0.12.1 or later is recommended)

* `OpenBabel <http://openbabel.org/>`_ (version 2.2.0 or later is recommended)

* `Cairo <http://cairographics.org/>`_ (version 1.8.0 or later is recommended)

* C and Fortran compilers

Compilation with Cython
-----------------------

Almost all of the ChemPy modules have been designed to be compiled into C 
extensions using Cython. This compilation is not required, but is strongly 
recommended due to the enormous speed boost that comes with it. To compile the
extensions, invoke the following from within the ChemPy root directory::

    $ python setup.py build_ext --inplace

You can also use the provided Makefile to do this::

    $ make

You can set extra options by creating a file in the ChemPy root directory
named make.inc.
