#!/usr/bin/env python
# encoding: utf-8

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2009-2011 by the RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

"""
This module contains classes and methods for working with physical quantities,
particularly the :class:`Quantity` class for representing physical quantities
and the :class:`Constants` class for defining relevant physical constants.
"""

import math
import numpy
import quantities as pq

# Explicity set the default units to SI
pq.set_default_units('si')

# These units are not defined by the quantities package, but occur frequently
# in data handled by RMG, so we define them manually
pq.UnitQuantity('kilocalories', pq.cal*1e3, symbol='kcal')
pq.UnitQuantity('kilojoules', pq.J*1e3, symbol='kJ')
pq.UnitQuantity('kilomoles', pq.mol*1e3, symbol='kmol')
pq.UnitQuantity('molecule', pq.mol/6.02214179e23, symbol='molecule')
pq.UnitQuantity('molecules', pq.mol/6.02214179e23, symbol='molecules')

################################################################################

def getConversionFactorToSI(units):
    """
    Get the conversion factor for converting a quantity in a given set of
    `units` to the SI equivalent units.
    """
    factor = float(pq.Quantity(1.0, units).simplified)
    # Exception: don't convert wavenumbers (cm^-1) to m^-1
    if units == 'cm^-1': factor = 1.0
    return factor

def getConversionFactorFromSI(units):
    """
    Get the conversion factor for converting a quantity to a given set of
    `units` from the SI equivalent units.
    """
    factor = float(pq.Quantity(1.0, units).simplified)
    # Exception: don't convert wavenumbers (cm^-1) to m^-1
    if units == 'cm^-1': factor = 1.0
    return 1.0 / factor

class Quantity:
    """
    A single numeric quantity or an array of quantities, with optional units
    and uncertainty. The attributes are:

    =================== =================== ====================================
    Attribute           Type                Description
    =================== =================== ====================================
    `value`             ``double``          The numeric value of the quantity in SI units
    `units`             ``str``             The units the value was specified in
    `uncertainty`       ``double``          The numeric uncertainty in the value
    `uncertaintyType`   ``str``             The type of uncertainty: ``'+|-'`` for additive, ``'*|/'`` for multiplicative
    `values`            ``numpy.ndarray``   The numeric values of the quantity in SI units
    `uncertainties`     ``numpy.ndarray``   The numeric uncertainty of the values in SI units
    =================== =================== ====================================

    Only one of `value` and `values` is meaningful at a time. If `values` is
    not ``None``, then the class assumes that it represents an array of
    quantities and uses the `values` attribute. If `values` is `None`, the
    class assumed that it represents a single quantity and uses the `value`
    attribute. You can test for this using the :meth:`isArray` method.

    In order to preserve the efficiency that Cythonization provides, the actual
    value of the quantity is stored in SI units in the `value` (for a single
    quantity) or `values` (for an array of quantities) attributes. The methods
    :meth:`getConversionFactorToSI` and :meth:`getConversionFactorFromSI` have
    been provided to facilitate unit conversions when reading and writing the
    quantity data. These methods use the ``quantities`` package, and are
    generally not optimized for speed.

    When providing an array of quantities, you have the option of providing
    no uncertainty data, a single uncertainty value to use for all values in
    the array, or an array of uncertainties corresponding to each value. When
    providing a single quantity, you can either provide no uncertainty or a
    single uncertainty value for that quantity. The uncertainty data can be
    specified as either additive or multiplicative, and must be symmetric.
    """

    def __init__(self, args=None):
        # Set default attribute values
        self.value = 0.0
        self.values = None
        self.units = ''
        self.uncertaintyType = ''
        self.uncertainty = 0.0
        self.uncertainties = None

        # Unpack arguments
        # If given, the order should be: value(s), units[, uncertaintyType], uncertainty
        units = ''; uncertaintyType = ''; uncertainty = None
        if args is None:
            return
        elif isinstance(args, Quantity):
            self.value = args.value
            self.values = args.values
            self.units = args.units
            self.uncertaintyType = args.uncertaintyType
            self.uncertainty = args.uncertainty
            self.uncertainties = args.uncertainties
            return
        elif isinstance(args, numpy.ndarray):
            # We've been given just an array of numbers
            value = args
        elif isinstance(args, float) or isinstance(args, int):
            # We've been given just a single number
            value = args
        elif isinstance(args, list) or isinstance(args, tuple):
            if len(args) == 1:
                value = args
            elif len(args) == 2:
                value, units = args
            elif len(args) == 3:
                value, units, uncertainty = args
            elif len(args) == 4:
                value, units, uncertaintyType, uncertainty = args
            else:
                value = args
            if not isinstance(units, str) or not isinstance(uncertaintyType, str):
                value = args; units = ''; uncertaintyType = ''; uncertainty = None

        # Process value parameter
        if isinstance(value, list) or isinstance(value, tuple):
            self.value = 0.0
            self.values = numpy.array(value, numpy.float64)
        elif isinstance(value, numpy.ndarray):
            self.value = 0.0
            self.values = value
        elif isinstance(value, float) or isinstance(value, int):
            self.value = value
            self.values = None
        else:
            raise ValueError('Unexpected type "%s" for value parameter.' % (value.__class__))

        # Process units and uncertainty type parameters
        self.units = units
        self.uncertaintyType = uncertaintyType

        # Process uncertainty parameter
        if isinstance(uncertainty, list) or isinstance(uncertainty, tuple):
            self.uncertainty = 0.0
            self.uncertainties = numpy.array(uncertainty, numpy.float64)
        elif isinstance(uncertainty, numpy.ndarray):
            self.uncertainty = 0.0
            self.uncertainties = uncertainty
        elif isinstance(uncertainty, float) or isinstance(uncertainty, int):
            self.uncertainty = uncertainty
            self.uncertainties = None
        elif uncertainty is not None:
            raise ValueError('Unexpected type "%s" for uncertainty parameter.' % (uncertainty.__class__))

        # Having multiple uncertainties for a single value is nonsensical
        # Add assertion to ensure we don't ever try to do this
        assert self.values is not None or self.uncertainties is None, "Attempted to create Quantity with one value but multiple uncertainties."
        # If multiple values and multiple uncertainties are specified, they
        # should be of equal length
        # Add assertion to ensure this is the case
        if self.values is not None and self.uncertainties is not None:
            assert len(self.values) == len(self.uncertainties), "Provided multiple uncertainties of different length than multiple values."

        # Get the output (SI) units corresponding to the input units
        factor = getConversionFactorToSI(self.units)

        # Convert the value and uncertainty to SI units
        self.value *= factor
        if self.values is not None: self.values *= factor
        if not self.isUncertaintyAdditive(): factor = 1.0
        self.uncertainty *= factor
        if self.uncertainties is not None: self.uncertainties *= factor

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        d = {}
        d['value'] = self.value
        d['units'] = self.units
        d['uncertaintyType'] = self.uncertaintyType
        d['uncertainty'] = self.uncertainty
        d['values'] = self.values
        d['uncertainties'] = self.uncertainties
        return (Quantity, tuple([None]), d)

    def __setstate__(self, d):
        """
        A helper function used when unpickling an object.
        """
        self.value = d['value']
        self.units = d['units']
        self.uncertaintyType = d['uncertaintyType']
        self.uncertainty = d['uncertainty']
        self.values = d['values']
        self.uncertainties = d['uncertainties']

    def __str__(self):
        """
        Return a string representation of the object.
        """
        factor = getConversionFactorFromSI(self.units)
        string = ''
        if not self.isArray():
            string += '%g' % (self.value * factor)
        else:
            string += '[%s]' % (','.join(['%g' % (v * factor) for v in self.values]))
        if self.uncertaintyType != '':
            string += ' %s ' % (self.uncertaintyType)
            if not self.isUncertaintyAdditive(): factor = 1.0
            if self.uncertainties is None:
                string += '%g' % (self.uncertainty * factor)
            else:
                string += '[%s]' % (','.join(['%g' % (u * factor) for u in self.uncertainties]))
        string += ' %s' % (self.units)
        return string

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        object.
        """
        factor = getConversionFactorFromSI(self.units)
        string = ''
        if not self.isArray():
            string += '%g' % (self.value * factor)
        else:
            string += '[%s]' % (','.join(['%g' % (v * factor) for v in self.values]))
        if self.units != '' or self.uncertaintyType != '':
            string += ',"%s"' % (self.units)
            if self.uncertaintyType != '' and (self.uncertainty != 0 or self.uncertainties is not None):
                string += ',"%s"' % (self.uncertaintyType)
                if not self.isUncertaintyAdditive(): factor = 1.0
                if self.uncertainties is None:
                    string += ',%g' % (self.uncertainty * factor)
                else:
                    string += ',[%s]' % (','.join(['%g' % (u * factor) for u in self.uncertainties]))
            string = '(' + string + ')'
        return string

    def __add__(self, other):
        """
        Add two objects together and return the sum as a new object.
        """
        if not isinstance(other, Quantity): raise ValueError('Unexpected type "%s" for other parameter.' % (other.__class__))
        assert self.units == other.units
        assert self.uncertaintyType in ['', '+|-']
        assert other.uncertaintyType in ['', '+|-']

        new = Quantity()
        if self.isArray() and self.uncertainties is not None:
            new.values = self.values + other.values
            new.units = self.units
            new.uncertaintyType = self.uncertaintyType
            new.uncertainties = numpy.sqrt(self.uncertainties**2 + other.uncertainties**2)
        elif self.isArray() and self.uncertainties is None:
            new.values = self.values + other.values
            new.units = self.units
            new.uncertaintyType = self.uncertaintyType
            new.uncertainty = numpy.sqrt(self.uncertainty**2 + other.uncertainty**2)
        else:
            new.value = self.value + other.value
            new.units = self.units
            new.uncertaintyType = self.uncertaintyType
            new.uncertainty = numpy.sqrt(self.uncertainty**2 + other.uncertainty**2)
        return new
        
    def isArray(self):
        """
        Return ``True`` if this quantity contains an array of values or 
        ``False`` if it contains only a single value. If the former, you should
        use the `values` attribute to get the values (in SI units) of the
        quantities as a NumPy array. If the latter, you should use the `value` 
        attribute to get the value of the quantity as a ``float``.
        """
        return self.values is not None

    def isUncertaintyAdditive(self):
        """
        Return ``True`` if the uncertainty is specified in additive format
        and ``False`` otherwise.
        """
        return self.uncertaintyType == '+|-'

    def isUncertaintyMultiplicative(self):
        """
        Return ``True`` if the uncertainty is specified in multiplicative format
        and ``False`` otherwise.
        """
        return self.uncertaintyType == '*|/'

################################################################################

class Constants:
    """
    A class defining several physical constants:

    =============== =========== ================================================
    Attribute       Type        Description
    =============== =========== ================================================
    `Na`            ``double``  The Avogadro constant :math:`N_\\mathrm{A}`, in :math:`\\mathrm{mol^{-1}}`
    `kB`            ``double``  The Boltzmann constant :math:`k_\\mathrm{B}`, in :math:`\\mathrm{J/K}`
    `R`             ``double``  The gas law constant :math:`R`, in :math:`\\mathrm{J/mol \\cdot K}`
    `h`             ``double``  The Planck constant :math:`h`, in :math:`\\mathrm{J \\cdot s}`
    `c`             ``double``  The speed of light in a vacuum :math:`c`, in :math:`\\mathrm{m/s}`
    `pi`            ``double``  The mathematical constant :math:`\\pi = 3.14159...`
    =============== =========== ================================================
    
    """
    
    def __init__(self):
        self.Na = 6.02214179e23
        self.kB = 1.3806504e-23
        self.R = 8.314472
        self.h = 6.62606896e-34
        self.c = 299792458
        self.pi = float(math.pi)

# An instance of the Constants class providing easy access to the physical constants
constants = Constants()
