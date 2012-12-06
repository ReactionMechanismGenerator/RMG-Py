****************
Quantity objects
****************

A physical quantity is defined by a numerical value and a unit of measurement.

The :mod:`rmgpy.quantity` module contains classes and methods for working with
physical quantities. Physical quantities are represented by either the
:class:`ScalarQuantity <rmgpy.quantity.ScalarQuantity>` or
:class:`ArrayQuantity <rmgpy.quantity.ArrayQuantity>` class depending on
whether a scalar or vector (or tensor) value is used. The 
:func:`Quantity <rmgpy.quantity.Quantity>` function automatically chooses the
appropriate class based on the input value. In both cases, the value of a
physical quantity is available from the ``value`` attribute, and the units from
the ``units`` attribute. 

For efficient computation, the value is stored internally in the SI equivalent
units. The SI value can be accessed directly using the ``value_si`` attribute.
Usually it is good practice to read the ``value_si`` attribute into a local
variable and then use it for computations, especially if it is referred to
multiple times in the calculation.

Physical quantities also allow for storing of uncertainty values for both
scalars and arrays. The ``uncertaintyType`` attribute indicates whether the
given uncertainties are additive (``"+|-"``) or multiplicative (``"*|/"``),
and the ``uncertainty`` attribute contains the stored uncertainties. For
additive uncertainties these are stored in the given units (not the SI
equivalent), since they are generally not needed for efficient computations.
For multiplicative uncertainties, the uncertainty values are by definition
dimensionless.

Generating quantities
=====================

.. autofunction:: rmgpy.quantity.Quantity

Scalar quantities
=================

.. autoclass:: rmgpy.quantity.ScalarQuantity

Array quantities
================

.. autoclass:: rmgpy.quantity.ArrayQuantity

