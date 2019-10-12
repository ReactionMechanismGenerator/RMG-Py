***************************************
Pressure dependence (:mod:`rmgpy.pdep`)
***************************************

.. module:: rmgpy.pdep

The :mod:`rmgpy.pdep` subpackage provides functionality for calcuating the
pressure-dependent rate coefficients :math:`k(T,P)` for unimolecular reaction
networks.

A unimolecular reaction network is defined by a set of chemically reactive 
molecular configurations - local minima on a potential energy surface - divided
into unimolecular isomers and bimolecular reactants or products. In our 
vernacular, reactants can associate to form an isomer, while such association
is neglected for products. These configurations are connected by chemical 
reactions to form a network; these are referred to as *path* reactions. The
system also consists of an excess of inert gas M, representing a thermal bath;
this allows for neglecting all collisions other than those between an isomer
and the bath gas. 

An isomer molecule at sufficiently high internal energy can be transformed by a
number of possible events:

* The isomer molecule can collide with any other molecule, resulting in an
  increase or decrease in energy

* The isomer molecule can isomerize to an adjacent isomer at the same energy

* The isomer molecule can dissociate into any directly connected bimolecular
  reactant or product channel

It is this competition between collision and reaction events that gives rise
to pressure-dependent kinetics.


Collision events
================

.. currentmodule:: rmgpy.pdep

=============================== ================================================
Class                           Description
=============================== ================================================
:class:`SingleExponentialDown`  A collisional energy transfer model based on the single exponential down model
=============================== ================================================



Reaction events
===============

.. currentmodule:: rmgpy.pdep

================================================= ================================
Function                                          Description
================================================= ================================
:func:`calculate_microcanonical_rate_coefficient` Return the microcanonical rate coefficient :math:`k(E)` for a reaction
:func:`apply_rrkm_theory`                         Use RRKM theory to compute :math:`k(E)` for a reaction
:func:`apply_inverse_laplace_transform_method`    Use the inverse Laplace transform method to compute :math:`k(E)` for a reaction
================================================= ================================



Pressure-dependent reaction networks
====================================

.. currentmodule:: rmgpy.pdep

======================= ========================================================
Class                   Description
======================= ========================================================
:class:`Configuration`  A molecular configuration on a potential energy surface
:class:`Network`        A collisional energy transfer model based on the single exponential down model
======================= ========================================================


The master equation
===================

.. currentmodule:: rmgpy.pdep.me

=============================== ================================================
Function                        Description
=============================== ================================================
:func:`generate_full_me_matrix` Return the full master equation matrix for a network
=============================== ================================================



Master equation reduction methods
=================================

.. currentmodule:: rmgpy.pdep

=========================================================== ====================
Function                                                    Description
=========================================================== ====================
:func:`msc.apply_modified_strong_collision_method`          Reduce the master equation to phenomenological rate coefficients :math:`k(T,P)` using the modified strong collision method
:func:`rs.apply_reservoir_state_method`                     Reduce the master equation to phenomenological rate coefficients :math:`k(T,P)` using the reservoir state method
:func:`cse.apply_chemically_significant_eigenvalues_method` Reduce the master equation to phenomenological rate coefficients :math:`k(T,P)` using the chemically-significant eigenvalues method
=========================================================== ====================



.. toctree::
    :hidden:
    
    singleexponentialdown
    reaction
    configuration
    network
    mastereqn
    methods
