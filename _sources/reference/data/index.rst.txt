.. _dataapireference:

****************************
Database (:mod:`rmgpy.data`)
****************************

.. module:: rmgpy.data

General classes
===============

.. currentmodule:: rmgpy.data.base

=========================== ====================================================
Class/Function              Description
=========================== ====================================================
:class:`Entry`              An entry in a database
:class:`Database`           A database of entries
:class:`LogicNode`          A node in a database that represents a logical collection of entries
:class:`LogicAnd`           A logical collection of entries, where all entries in the collection must match
:class:`LogicOr`            A logical collection of entries, where any entry in the collection can match
:func:`makeLogicNode`       Create a :class:`LogicNode` based on a string representation
=========================== ====================================================



Thermodynamics database
=======================

.. currentmodule:: rmgpy.data.thermo

=========================== ====================================================
Class                       Description
=========================== ====================================================
:class:`ThermoDepository`   A depository of all thermodynamics parameters for one or more species
:class:`ThermoLibrary`      A library of curated thermodynamics parameters for one or more species
:class:`ThermoGroups`       A representation of a portion of a database for implementing the Benson group additivity method
:class:`ThermoDatabase`     An entire thermodynamics database, including depositories, libraries, and groups
=========================== ====================================================



Kinetics database
=================

.. currentmodule:: rmgpy.data.kinetics

=========================== ====================================================
Class                       Description
=========================== ====================================================
:class:`DepositoryReaction` A reaction with kinetics determined from querying a kinetics depository
:class:`LibraryReaction`    A reaction with kinetics determined from querying a kinetics library
:class:`TemplateReaction`   A reaction with kinetics determined from querying a kinetics group additivity or rate rules method
:class:`ReactionRecipe`     A sequence of actions that represent the process of a chemical reaction
--------------------------- ----------------------------------------------------
:class:`KineticsDepository` A depository of all kinetics parameters for one or more reactions
:class:`KineticsLibrary`    A library of curated kinetics parameters for one or more reactions
:class:`KineticsGroups`     A set of group additivity values for a reaction family, organized in a tree
:class:`KineticsRules`      A set of rate rules for a reaction family
:class:`KineticsFamily`     A kinetics database for one reaction family, including depositories, libraries, groups, and rules
:class:`KineticsDatabase`   A kinetics database for all reaction families, including depositories, libraries, groups, and rules
=========================== ====================================================



Statistical mechanics database
==============================

.. currentmodule:: rmgpy.data.statmech

=========================== ====================================================
Class                       Description
=========================== ====================================================
:class:`GroupFrequencies`   A set of characteristic frequencies for a group in the frequency database
:class:`StatmechDepository` A depository of all statistical mechanics parameters for one or more species
:class:`StatmechLibrary`    A library of curated statistical mechanics parameters for one or more species
:class:`StatmechGroups`     A set of characteristic frequencies for various functional groups, organized in a tree
:class:`StatmechDatabase`   An entire statistical mechanics database, including depositories, libraries, and groups
=========================== ====================================================


Statistical mechanics fitting
=============================

.. currentmodule:: rmgpy.data.statmechfit

=================================== ============================================
Class/Function                      Description
=================================== ============================================
:class:`DirectFit`                  DQED class for fitting a small number of vibrational frequencies and hindered rotors
:class:`PseudoFit`                  DQED class for fitting a large number of vibrational frequencies and hindered rotors by assuming degeneracies for both
:class:`PseudoRotorFit`             DQED class for fitting a moderate number of vibrational frequencies and hindered rotors by assuming degeneracies for hindered rotors only
:func:`fitStatmechDirect`           Directly fit a small number of vibrational frequencies and hindered rotors
:func:`fitStatmechPseudo`           Fit a large number of vibrational frequencies and hindered rotors by assuming degeneracies for both
:func:`fitStatmechPseudoRotors`     Fit a moderate number of vibrational frequencies and hindered rotors by assuming degeneracies for hindered rotors only
:func:`fitStatmechToHeatCapacity`   Fit vibrational and torsional degrees of freedom to heat capacity data
=================================== ============================================



Exceptions
==========

.. currentmodule:: rmgpy.data

=================================== ============================================
Exception                           Description
=================================== ============================================
:exc:`DatabaseError`                Raised when an error occurs while working with the database
:exc:`InvalidActionError`           Raised when an error occurs while applying a reaction recipe
:exc:`UndeterminableKineticsError`  Raised when the kinetics of a given reaction cannot be determined
:exc:`StatmechFitError`             Raised when an error occurs while fitting internal degrees of freedom to heat capacity data
=================================== ============================================



.. toctree::
    :hidden:
    
    database
    depositoryreaction
    entry
    groupfrequencies
    kineticsdatabase
    kineticsdepository
    kineticsfamily
    kineticsgroups
    kineticslibrary
    kineticsrules
    libraryreaction
    logicnode
    reactionrecipe
    statmechdatabase
    statmechdepository
    statmechfit
    statmechgroups
    statmechlibrary
    templatereaction
    thermodatabase
    thermodepository
    thermogroups
    thermolibrary
    
