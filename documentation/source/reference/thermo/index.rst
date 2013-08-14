************************************
Thermodynamics (:mod:`rmgpy.thermo`)
************************************

.. module:: rmgpy.thermo

The :mod:`rmgpy.thermo` subpackage contains classes that represent various
thermodynamic models of heat capacity.

Heat capacity models
====================

.. currentmodule:: rmgpy.thermo

======================= ========================================================
Class                   Description
======================= ========================================================
:class:`ThermoData`     A heat capacity model based on a set of discrete heat capacity points
:class:`Wilhoit`        A heat capacity model based on the Wilhoit polynomial
:class:`NASA`           A heat capacity model based on a set of NASA polynomials
:class:`NASAPolynomial` A heat capacity model based on a single NASA polynomial
======================= ========================================================

.. toctree::
    :hidden:
    
    thermodata
    wilhoit
    nasa
    nasapolynomial
