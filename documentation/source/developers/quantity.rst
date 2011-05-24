***********************************************************
:mod:`rmgpy.quantity` --- Physical Quantities and Constants
***********************************************************

.. automodule:: rmgpy.quantity

Physical Constants
==================

The constants available are listed below. All values were taken from
`NIST <http://physics.nist.gov/cuu/Constants/index.html>`_.

.. autoclass:: rmgpy.quantity.Constants
    :members:

Quantity Objects
================

Many of the physical quantities encountered throughout the ``chem`` package
have associated units and/or uncertainties. To represent these quantities in
a consistent manner, the :class:`Quantity` class has been provided.
:class:`Quantity` objects are primarily used with the classes in the ``chem``
package to store the physical quantities behind-the-scenes. Generally you 
should not need to create :class:`Quantity` objects by hand.

.. autoclass:: rmgpy.quantity.QuantityError
    :members:

.. autoclass:: rmgpy.quantity.Quantity
    :members:
