***************************************************
:mod:`rmgpy.chem.constants` --- Numerical Constants
***************************************************

.. automodule:: rmgpy.chem.constants

Physical Constants
==================

The constants available are listed below. All values were taken from
`NIST <http://physics.nist.gov/cuu/Constants/index.html>`_.

.. autodata:: rmgpy.chem.constants.Na

.. autodata:: rmgpy.chem.constants.kB

.. autodata:: rmgpy.chem.constants.R

.. autodata:: rmgpy.chem.constants.h

.. autodata:: rmgpy.chem.constants.c

.. autodata:: rmgpy.chem.constants.pi

Quantity Objects
================

Many of the physical quantities encountered throughout the ``chem`` package
have associated units and/or uncertainties. To represent these quantities in
a consistent manner, the :class:`Quantity` class has been provided.
:class:`Quantity` objects are primarily used with the classes in the ``chem``
package to store the physical quantities behind-the-scenes. Generally you 
should not need to create :class:`Quantity` objects by hand.

.. autoclass:: rmgpy.chem.constants.Quantity
    :members:

.. autofunction:: rmgpy.chem.constants.getConversionFactorToSI

.. autofunction:: rmgpy.chem.constants.getConversionFactorFromSI
