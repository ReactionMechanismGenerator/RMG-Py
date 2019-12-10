****************************
Arkane (:mod:`arkane`)
****************************

.. module:: arkane

The :mod:`arkane` subpackage contains the main functionality for
Arkane, a tool for computing thermodynamic and kinetic properties of chemical
species and reactions.



Reading electronic structure software log files
===============================================

.. currentmodule:: arkane.ess

=============================== ================================================
Class                           Description
=============================== ================================================
:class:`Log`                    Base class for generic log files
:class:`GaussianLog`            Extract chemical parameters from Gaussian log files
:class:`MolproLog`              Extract chemical parameters from Molpro log files
:class:`OrcaLog`                Extract chemical parameters from Orca log files
:class:`QChemLog`               Extract chemical parameters from Q-Chem log files
:class:`TeraChemLog`            Extract chemical parameters from TeraChem log files
=============================== ================================================



Input
=====

.. currentmodule:: arkane.input

=============================== ================================================
Function                        Description
=============================== ================================================
:func:`load_input_file`         Load an Arkane job input file
=============================== ================================================



Output
======

.. currentmodule:: arkane.output

=============================== ================================================
Function                        Description
=============================== ================================================
:class:`PrettifyVisitor`        Custom Abstract Syntax Tree (AST) visitor class
:func:`prettify`                Pretty formatting for a Python syntax string
:func:`get_str_xyz`             Pretty formatting for XYZ coordinates
:func:`save_thermo_lib`         Save an RMG thermo library
:func:`save_kinetics_lib`       Save an RMG kinetics library
=============================== ================================================



Job classes
===========

.. currentmodule:: arkane

=============================== ================================================
Class                           Description
=============================== ================================================
:class:`Arkane`                 Main class for Arkane jobs
:class:`StatMechJob`            Compute the molecular degrees of freedom for a molecular conformation
:class:`ThermoJob`              Compute the thermodynamic properties of a species
:class:`KineticsJob`            Compute the high pressure-limit rate coefficient for a reaction using transition state theory
:class:`PressureDependenceJob`  Compute the phenomenological pressure-dependent rate coefficients :math:`k(T,P)` for a unimolecular reaction network
=============================== ================================================



Sensitivity analysis
====================

.. currentmodule:: arkane.sensitivity

=============================== ================================================
Class                           Description
=============================== ================================================
:class:`KineticsSensitivity`    Perform sensitivity analysis for a kinetics job
:class:`PDepSensitivity`        Perform sensitivity analysis for a pressure dependence job
=============================== ================================================



.. toctree::
    :hidden:
    
    log
    gaussianlog
    molprolog
    orcalog
    qchemlog
    terachemlog
    input
    kinetics
    main
    output
    pdep
    statmech
    thermo
    sensitivity
    
