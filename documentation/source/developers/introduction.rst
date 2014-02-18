************
Introduction
************

The :mod:`rmgpy` package contains all of the functionality of RMG.

RMG-Py naming conventions in brief:

* module names are in lowercase 
* Class and Exception names start with capital letters (e.g., CamelCase)
* all instances of Classes are mixed-case and begin with lower-case letters (e.g., camelCase)


 `Python style guide <http://www.python.org/dev/peps/pep-0008/>`_. In general, refer to this for guidance.


Dependencies
============

This package builds on a number of Python packages (in addition to those in the
Python standard library):

* `Cython <http://www.cython.org/>`_. Provides a means to compile annotated 
  Python modules to C, combining the rapid development of Python with near-C
  execution speeds.

* `NumPy <http://numpy.scipy.org/>`_. Provides efficient matrix algebra.

* `SciPy <http://www.scipy.org/>`_. Extends NumPy with a variety of mathematics 
  tools useful in scientific computing.

* `OpenBabel <http://openbabel.org/>`_. Provides functionality for converting
  between a variety of chemical formats.

* `Cairo <http://cairographics.org/>`_. Provides functionality for generation
  of 2D graphics figures.

Cython Support
==============

All of the contained modules are written in pure `Python <http://www.python.org/>`_.
However, since certain operations can be computationally expensive -- either
due to the operation itself or simply the frequency by which the operation is
performed -- most of the modules can also be compiled using 
`Cython <http://www.cython.org/>`_ for a signficant speed boost. Cython is a
modification of the Python programming language with syntax to enable 
compilation to efficient C code, particularly through the ability to statically
declare the type of variables at compile time. Most Python syntax can be
compiled via Cython as-is (especially with recent versions of Cython), though
you won't see much of a speed increase without declaring static types in your
bottleneck code.

.. note:: 
    Our use of Cython inherently ties RMG Py to the CPython implementation
    of the Python language. Jython (Java), IronPython (Common Language
    Runtime), etc. cannot be used with compiled Cython modules at this time.

Cython is one of several attempts to combine the ease and flexibility of Python
with the efficiency of a compiled language. As of this writing, a summary and
simple comparison of several such approaches is available from 
`<http://www.scipy.org/PerformancePython>`_. That page suggests that
Cython [1]_ is among the best at approaching the efficiency of a pure C
implementation.

.. [1] The authors actually used Pyrex, a similar project that Cython was 
        forked from, in their comparison. Cython has advanced significantly
        since that document was written, and today likely performs even more 
        efficiently than Pyrex.

