************
Introduction
************

The :mod:`rmgpy.chem` module provides a library of classes and functions for
chemistry, chemical engineering, and materials science applications. 

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


Modules
=======

The :mod:`rmgpy.chem` package contains a wide variety of classes and functions,
grouped by purpose into the following modules. The remainder of this document
will discuss each of these modules and their contents in detail:

=============================== ================================================
Module                          Description
=============================== ================================================
:mod:`rmgpy.chem.constants`     Declaration of common physical constants (in SI units)
:mod:`rmgpy.chem.element`       Representation of chemical elements
:mod:`rmgpy.chem.geometry`      Operations on 3D molecular geometries
:mod:`rmgpy.chem.graph`         Efficient implementation of a graph data type and graph isomorphism algorithm
:mod:`rmgpy.chem.kinetics`      Representation of kinetics models for chemical reaction rate coefficients
:mod:`rmgpy.chem.molecule`      Representation of chemical molecules and their constituent atoms and bonds
:mod:`rmgpy.chem.pattern`       Representation of molecular substructure patterns
:mod:`rmgpy.chem.reaction`      Representation of chemical reactions
:mod:`rmgpy.chem.species`       Representation of chemical species
:mod:`rmgpy.chem.states`        Representation of internal and external molecular degrees of freedom
:mod:`rmgpy.chem.thermo`        Representation of thermodynamics models for chemical species
=============================== ================================================

In addition to the above, there are a few extra modules providing extended
functionality:

======================================= ========================================
Module                                  Description
======================================= ========================================
:mod:`rmgpy.chem.ext.thermo_converter`  Conversion between various thermodynamics models
:mod:`rmgpy.chem.ext.molecule_draw`     2D depiction of structural formula of chemical molecules
======================================= ========================================



