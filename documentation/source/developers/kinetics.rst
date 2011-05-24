**********************************************
:mod:`rmgpy.chem.kinetics` --- Kinetics Models
**********************************************

.. automodule:: rmgpy.chem.kinetics

Base Class for Kinetics Models
==============================

.. autoclass:: rmgpy.chem.kinetics.KineticsError
    :members:

.. autoclass:: rmgpy.chem.kinetics.KineticsModel
    :members:

Pressure-Independent Kinetics Models
====================================

.. autoclass:: rmgpy.chem.kinetics.Arrhenius
    :members:

.. autoclass:: rmgpy.chem.kinetics.ArrheniusEP
    :members:

.. autoclass:: rmgpy.chem.kinetics.MultiArrhenius
    :members:

Pressure-Dependent Kinetics Models
==================================

Single-Well Models
------------------

.. autoclass:: rmgpy.chem.kinetics.ThirdBody
    :members:

.. autoclass:: rmgpy.chem.kinetics.Lindemann
    :members:

.. autoclass:: rmgpy.chem.kinetics.Troe
    :members:

.. autoclass:: rmgpy.chem.kinetics.PDepArrhenius
    :members:

Multi-Well Models
-----------------

.. autoclass:: rmgpy.chem.kinetics.Chebyshev
    :members:

