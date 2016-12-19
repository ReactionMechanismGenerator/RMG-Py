*********************************************
Statistical mechanics (:mod:`rmgpy.statmech`)
*********************************************

.. module:: rmgpy.statmech

The :mod:`rmgpy.statmech` subpackage contains classes that represent various
statistical mechanical models of molecular degrees of freedom. These models
enable the computation of macroscopic parameters (e.g. thermodynamics, kinetics,
etc.) from microscopic parameters.

A molecular system consisting of :math:`N` atoms is described by :math:`3N`
molecular degrees of freedom. Three of these modes involve translation of the
system as a whole. Another three of these modes involve rotation of the system
as a whole, unless the system is linear (e.g. diatomics), for which there are
only two rotational modes. The remaining :math:`3N-6` (or :math:`3N-5` if 
linear) modes involve internal motion of the atoms within the system. Many of
these modes are well-described as harmonic oscillations, while others are
better modeled as torsional rotations around a bond within the system.

Molecular degrees of freedom are mathematically represented using the
Schrodinger equation :math:`\hat{H} \Psi = E \Psi`. By solving the
Schrodinger equation, we can determine the available energy states of the
molecular system, which enables computation of macroscopic parameters. 
Depending on the temperature of interest, some modes (e.g. vibrations) require 
a quantum mechanical treatment, while others (e.g. translation, rotation) can 
be described using a classical solution.



Translational degrees of freedom
================================

.. currentmodule:: rmgpy.statmech

=============================== ================================================
Class                           Description
=============================== ================================================
:class:`IdealGasTranslation`    A model of three-dimensional translation of an ideal gas
=============================== ================================================



Rotational degrees of freedom
=============================

.. currentmodule:: rmgpy.statmech

=========================== ====================================================
Class                       Description
=========================== ====================================================
:class:`LinearRotor`        A model of two-dimensional rigid rotation of a linear molecule
:class:`NonlinearRotor`     A model of three-dimensional rigid rotation of a nonlinear molecule
:class:`KRotor`             A model of one-dimensional rigid rotation of a K-rotor
:class:`SphericalTopRotor`  A model of three-dimensional rigid rotation of a spherical top molecule
=========================== ====================================================


Vibrational degrees of freedom
==============================

.. currentmodule:: rmgpy.statmech

=========================== ====================================================
Class                       Description
=========================== ====================================================
:class:`HarmonicOscillator` A model of a set of one-dimensional harmonic oscillators
=========================== ====================================================


Torsional degrees of freedom
============================

.. currentmodule:: rmgpy.statmech

=========================== ====================================================
Class                       Description
=========================== ====================================================
:class:`HinderedRotor`      A model of a one-dimensional hindered rotation
=========================== ====================================================


The Schrodinger equation
========================

.. currentmodule:: rmgpy.statmech.schrodinger

=============================== ================================================
Class                           Description
=============================== ================================================
:func:`getPartitionFunction`    Calculate the partition function at a given temperature from energy levels and degeneracies
:func:`getHeatCapacity`         Calculate the dimensionless heat capacity at a given temperature from energy levels and degeneracies
:func:`getEnthalpy`             Calculate the enthalpy at a given temperature from energy levels and degeneracies
:func:`getEntropy`              Calculate the entropy at a given temperature from energy levels and degeneracies
:func:`getSumOfStates`          Calculate the sum of states for a given energy domain from energy levels and degeneracies
:func:`getDensityOfStates`      Calculate the density of states for a given energy domain from energy levels and degeneracies
=============================== ================================================



Convolution
===========

.. currentmodule:: rmgpy.statmech.schrodinger

======================= ========================================================
Class                   Description
======================= ========================================================
:func:`convolve`        Return the convolution of two arrays
:func:`convolveBS`      Convolve a degree of freedom into a density or sum of states using the Beyer-Swinehart (BS) direct count algorithm
:func:`convolveBSSR`    Convolve a degree of freedom into a density or sum of states using the Beyer-Swinehart-Stein-Rabinovitch (BSSR) direct count algorithm
======================= ========================================================



Molecular conformers
====================

.. currentmodule:: rmgpy.statmech

======================= ========================================================
Class                   Description
======================= ========================================================
:class:`Conformer`      A model of a molecular conformation
======================= ========================================================



.. toctree::
    :hidden:
    
    idealgastranslation
    linearrotor
    nonlinearrotor
    krotor
    sphericaltoprotor
    harmonicoscillator
    hinderedrotor
    schrodinger
    conformer
