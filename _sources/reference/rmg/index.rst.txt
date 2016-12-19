************************************************
Reaction mechanism generation (:mod:`rmgpy.rmg`)
************************************************

.. module:: rmgpy.rmg

The :mod:`rmgpy.rmg` subpackage contains the main functionality for using
RMG-Py to automatically generate detailed reaction mechanisms.



Reaction models
===============

.. currentmodule:: rmgpy.rmg.model

=============================== ================================================
Class                           Description
=============================== ================================================
:class:`Species`                A chemical species, with RMG-specific functionality
:class:`CoreEdgeReactionModel`  A reaction model comprised of core and edge species and reactions
=============================== ================================================



Input
=====

.. currentmodule:: rmgpy.rmg.input

=============================== ================================================
Function                        Description
=============================== ================================================
:func:`readInputFile`           Load an RMG job input file
:func:`saveInputFile`           Save an RMG job input file
=============================== ================================================



Output
======

.. currentmodule:: rmgpy.rmg.output

=============================== ================================================
Function                        Description
=============================== ================================================
:func:`saveOutputHTML`          Save the results of an RMG job to an HTML file
:func:`saveDiffHTML`            Save a comparison of two reaction mechanisms to an HTML file
=============================== ================================================



Job classes
===========

.. currentmodule:: rmgpy.rmg.main

=============================== ================================================
Class                           Description
=============================== ================================================
:class:`RMG`                    Main class for RMG jobs
=============================== ================================================



Pressure dependence
===================

.. currentmodule:: rmgpy.rmg.pdep

=============================== ================================================
Class                           Description
=============================== ================================================
:class:`PDepReaction`           A pressure-dependent "net" reaction
:class:`PDepNetwork`            A pressure-dependent unimolecular reaction network, with RMG-specific functionality
=============================== ================================================



Exceptions
==========

.. currentmodule:: rmgpy.rmg

=============================== ================================================
Exception                       Description
=============================== ================================================
:exc:`InputError`               Raised when an error occurs while working with an RMG input file
:exc:`OutputError`              Raised when an error occurs while saving an RMG output file
:exc:`PressureDependenceError`  Raised when an error occurs while computing pressure-dependent rate coefficients
=============================== ================================================



.. toctree::
    :hidden:
    
    coreedgereactionmodel
    input
    main
    output
    pdepnetwork
    pdepreaction
    species
    
