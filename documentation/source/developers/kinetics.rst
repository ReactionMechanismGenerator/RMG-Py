**********************************************
:mod:`rmgpy.kinetics` --- Kinetics Models
**********************************************

.. automodule:: rmgpy.kinetics

Base Class for Kinetics Models
==============================

.. autoclass:: rmgpy.kinetics.KineticsError
    :members:

.. autoclass:: rmgpy.kinetics.KineticsModel
    :members:

Pressure-Independent Kinetics Models
====================================

.. autoclass:: rmgpy.kinetics.KineticsData
    :members:

.. autoclass:: rmgpy.kinetics.Arrhenius
    :members:

.. autoclass:: rmgpy.kinetics.ArrheniusEP
    :members:

Pressure-Dependent Kinetics Models
==================================

Single-Well Models
------------------

.. autoclass:: rmgpy.kinetics.ThirdBody
    :members:

.. autoclass:: rmgpy.kinetics.Lindemann
    :members:

.. autoclass:: rmgpy.kinetics.Troe
    :members:

.. autoclass:: rmgpy.kinetics.PDepArrhenius
    :members:

Multi-Well Models
-----------------

.. autoclass:: rmgpy.kinetics.Chebyshev
    :members:

Other Kinetics Models
=====================

.. autoclass:: rmgpy.kinetics.MultiKinetics
    :members:
