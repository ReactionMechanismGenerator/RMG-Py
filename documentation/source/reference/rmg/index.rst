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
:class:`CoreEdgeReactionModel`  A reaction model comprised of core and edge species and reactions
=============================== ================================================



Input
=====

.. currentmodule:: rmgpy.rmg.input

=============================== ================================================
Function                        Description
=============================== ================================================
:func:`read_input_file`         Load an RMG job input file
:func:`save_input_file`         Save an RMG job input file
=============================== ================================================



Output
======

.. currentmodule:: rmgpy.rmg.output

=============================== ================================================
Function                        Description
=============================== ================================================
:func:`save_output_html`        Save the results of an RMG job to an HTML file
:func:`save_diff_html`          Save a comparison of two reaction mechanisms to an HTML file
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



.. toctree::
    :hidden:
    
    coreedgereactionmodel
    input
    main
    output
    pdepnetwork
    pdepreaction
