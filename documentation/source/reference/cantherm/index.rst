********************************
CanTherm (:mod:`rmgpy.cantherm`)
********************************

.. module:: rmgpy.cantherm

The :mod:`rmgpy.cantherm` subpackage contains the main functionality for
CanTherm, a tool for computing thermodynamic and kinetic properties of chemical
species and reactions.



Reading Gaussian log files
==========================

.. currentmodule:: rmgpy.cantherm.gaussian

=============================== ================================================
Class                           Description
=============================== ================================================
:class:`GaussianLog`            Extract chemical parameters from Gaussian log files
=============================== ================================================



Reading Q-Chem log files
========================

.. currentmodule:: rmgpy.cantherm.qchem

=============================== ================================================
Class                           Description
=============================== ================================================
:class:`QchemLog`               Extract chemical parameters from Q-Chem log files
=============================== ================================================



Reading Molpro log files
========================

.. currentmodule:: rmgpy.cantherm.molpro

=============================== ================================================
Class                           Description
=============================== ================================================
:class:`MolproLog`              Extract chemical parameters from Molpro log files
=============================== ================================================



Input
=====

.. currentmodule:: rmgpy.cantherm.input

=============================== ================================================
Function                        Description
=============================== ================================================
:func:`loadInputFile`           Load a CanTherm job input file
=============================== ================================================



Job classes
===========

.. currentmodule:: rmgpy.cantherm

=============================== ================================================
Class                           Description
=============================== ================================================
:class:`CanTherm`               Main class for CanTherm jobs
:class:`StatMechJob`            Compute the molecular degrees of freedom for a molecular conformation
:class:`ThermoJob`              Compute the thermodynamic properties of a species
:class:`KineticsJob`            Compute the high pressure-limit rate coefficient for a reaction using transition state theory
:class:`PressureDependenceJob`  Compute the phenomenological pressure-dependent rate coefficients :math:`k(T,P)` for a unimolecular reaction network
=============================== ================================================



.. toctree::
    :hidden:
    
    gaussianlog
    qchemlog
    molprolog
    input
    kinetics
    main
    output
    pdep
    statmech
    thermo
    
