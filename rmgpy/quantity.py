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
import cython

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

class QuantityError(Exception):
    """
    An exception to be raised when an error occurs while working with physical
    quantities in RMG. Pass a string describing the circumstances of the
    exceptional behavior.
    """
    pass

################################################################################

class Quantity:
    """
    A representation of a physical quantity, with optional units and 
    uncertainty information. The attributes are:

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
    attribute. You can test for this using the :meth:`isArray()` method.

    In order to preserve the efficiency that Cythonization provides, the actual
    value of the quantity is stored in SI units in the `value` (for a single
    quantity) or `values` (for an array of quantities) attributes. The methods
    :meth:`getConversionFactorToSI()` and :meth:`getConversionFactorFromSI()` 
    have been provided to facilitate unit conversions when reading and writing 
    the quantity data. These methods use the ``quantities`` package, and are
    generally not optimized for speed.

    When providing an array of quantities, you have the option of providing
    no uncertainty data, a single uncertainty value to use for all values in
    the array, or an array of uncertainties corresponding to each value. When
    providing a single quantity, you can either provide no uncertainty or a
    single uncertainty value for that quantity. The uncertainty data can be
    specified as either additive or multiplicative, and must be symmetric.
    """
    
    # A dict of conversion factors (to SI) for each of the frequent units
    conversionFactors = {}
    
    def __init__(self, *args):
        """
        Create a new :class:`Quantity` object. The `args` parameter can contain
        from zero to four parameters. If no parameters are provided, the 
        attributes are initialized to default values. The first parameter should
        always be the value(s) of the quantity. The second parameter, if given, 
        should always be the units. The third and fourth parameters, if given,
        should be the uncertainty type and uncertainty values. If the
        uncertainty type is omitted but the uncertainty values are given, the
        type is assumed to be additive.
        """
        value = 0.0; units = ''; uncertaintyType = ''; uncertainty = 0.0
        
        # Unpack args if necessary
        if isinstance(args, tuple) and len(args) == 1 and isinstance(args[0], tuple):
            args = args[0]
            
        # Process args    
        if len(args) == 0:
            # No parameters were given, so initialize with default values
            self.value = 0.0
            self.values = None
            self.units = ''
            self.uncertaintyType = ''
            self.uncertainty = 0.0
            self.uncertainties = None
            return
        elif len(args) == 1:
            # If one parameter is given, it should be a single value or 
            # array of values; it could also be another Quantity object to
            # make a copy of
            other = args[0]
            if isinstance(other, Quantity):
                # We were given another Quantity object, so make a (shallow) copy of it
                self.value = other.value
                self.values = other.values
                self.units = other.units
                self.uncertaintyType = other.uncertaintyType
                self.uncertainty = other.uncertainty
                self.uncertainties = other.uncertainties
                return
            elif isinstance(other, list) or isinstance(other, tuple) or isinstance(other, numpy.ndarray):
                # We've been given just an array of values
                value = other
            elif isinstance(other, float) or isinstance(other, int):
                # We've been given just a single value
                value = other
        elif len(args) == 2:
            # If two parameters are given, it should be a single value or 
            # array of values, plus units
            value, units = args
        elif len(args) == 3:
            # If three parameters are given, it should be a single value or 
            # array of values, plus units and uncertainty
            value, units, uncertainty = args; 
            # Assume the intended uncertainty type is additive
            uncertaintyType = '+|-'
        elif len(args) == 4:
            # If three parameters are given, it should be a single value or 
            # array of values, plus units and uncertainty
            value, units, uncertaintyType, uncertainty = args
        else:
            raise QuantityError('Invalid parameters {0!r} passed to init method of Quantity object.'.format(args))
        
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
            raise QuantityError('Unexpected type "{0}" for value parameter.'.format(value.__class__))

        # Process units and uncertainty type parameters
        self.units = units
        self.uncertaintyType = uncertaintyType
        if uncertaintyType not in ['', '+|-', '*|/']:
            raise QuantityError('Unexpected uncertainty type "{0}"; valid values are "+|-" and "*|/".'.format(uncertaintyType))
            
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
            raise QuantityError('Unexpected type "{0}" for uncertainty parameter.'.format(uncertainty.__class__))
    
        # Having multiple uncertainties for a single value is nonsensical
        # Add a check to ensure we don't ever try to do this
        if self.values is None and self.uncertainties is not None:
            raise QuantityError('Attempted to create Quantity with one value but multiple uncertainties.')
        
        # If multiple values and multiple uncertainties are specified, they
        # should be of equal length
        # Add a check to ensure we don't ever try to do this
        if self.values is not None and self.uncertainties is not None and len(self.values) != len(self.uncertainties):
            raise QuantityError('Provided multiple uncertainties of different length than multiple values.')

        # Get the output (SI) units corresponding to the input units
        factor = self.getConversionFactorToSI()

        # Convert the value and uncertainty to SI units
        self.value *= factor
        if self.values is not None: self.values *= factor
        if not self.isUncertaintyAdditive(): factor = 1.0
        self.uncertainty *= factor
        if self.uncertainties is not None: self.uncertainties *= factor
        
    def __reduce__(self):
        """
        A helper function used when pickling a Quantity object.
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
        A helper function used when unpickling a Quantity object.
        """
        self.value = d['value']
        self.units = d['units']
        self.uncertaintyType = d['uncertaintyType']
        self.uncertainty = d['uncertainty']
        self.values = d['values']
        self.uncertainties = d['uncertainties']

    def __str__(self):
        """
        Return a human-readable string representation of the Quantity object.
        """
        factor = self.getConversionFactorFromSI()
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
        Quantity object.
        """
        factor = self.getConversionFactorFromSI()
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

    def getConversionFactorToSI(self):
        """
        Get the conversion factor for converting a quantity in a given set of
        `units` to the SI equivalent units.
        """
        cython.declare(factor=cython.double)
        try:
            # Process several common units manually for speed
            factor = Quantity.conversionFactors[self.units]
        except KeyError:
            # Fall back to (slow!) quantities package for less common units
            factor = float(pq.Quantity(1.0, self.units).simplified)
            # Cache the conversion factor so we don't ever need to use
            # quantities to compute it again
            Quantity.conversionFactors[self.units] = factor
        return factor

    def getConversionFactorFromSI(self):
        """
        Get the conversion factor for converting a quantity to a given set of
        `units` from the SI equivalent units.
        """
        cython.declare(factor=cython.double)
        factor = self.getConversionFactorToSI()
        return 1.0 / factor

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
        Return ``True`` if the uncertainty is specified in multiplicative 
        format and ``False`` otherwise.
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
