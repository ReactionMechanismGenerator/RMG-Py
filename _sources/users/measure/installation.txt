************
Installation
************

Dependencies
============

Like most Python programs, MEASURE relies on a number of Python packages for
certain functionality. The following lists the required Python packages 
that are not part of the Python standard library:

* `NumPy <http://www.numpy.org>`_. Provides efficient array and matrix 
  operations.

* `SciPy <http://www.scipy.org>`_. Provides efficient linear algebra functions
  and special functions.

* `Quantities <http://packages.python.org/quantities/>`_. Provides unit
  conversions of physical quantities.

* `Cython <http://www.cython.org/>`_. Compiles certain computationally-intensive
  modules to C code, often resulting in a dramatic speed increase, mostly due
  to static typing of variables.

The following lists the optional Python packages that are not part of the 
Python standard library:

* `Cairo <http://cairographics.org/>`_. Provides 2D drawing support, used for
  automatic drawing of potential energy surfaces.

Installing MEASURE
==================

.. note::

    MEASURE can also be used directly from the web! Visit 
    http://rmg.mit.edu/measure/ for more information.

MEASURE is installed automatically as part of the installation process for 
RMG-Py. You can obtain a copy of RMG-Py from 
http://github.com/GreenGroup/RMG-Py, either by downloading the current snapshot
or by cloning the repository using `git <http://git-scm.com/>`_. The latter is
recommended so that you can easily keep up to date with new features and 
bugfixes.

Once you have obtained the required dependencies, RMG-Py and MEASURE can be 
installed by executing the following command from within the root package 
directory::

    $ python setup.py install

This will move all of the necessary files to the appropriate Python installation
location on disk. If Cython is installed, this will also compile the appropriate
modules prior to installation.

You can also use RMG-Py and MEASURE without installing if you prefer. This is 
done by appending the full path to the package's root directory to your 
``PYTHONPATH`` environment variable. To compile the Cython modules in this 
case, use the command ::

    $ python setup.py build_ext --inplace

You can remove the temporary build files - but *not* the compiled libraries! -
via the command ::

    $ python setup.py clean

If you have the ``make`` build automation tool installed, you can also use the
provided ``Makefile``, via ::

    $ make install

to install, ::

    $ make

to build in-place, or ::

    $ make clean

to remove all temporary build files *and* compiled files.
