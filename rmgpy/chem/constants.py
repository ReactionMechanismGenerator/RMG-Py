#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
#   ChemPy - A chemistry toolkit for Python
#
#   Copyright (c) 2010 by Joshua W. Allen (jwallen@mit.edu)
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
This module contains a number of physical constants to be made available
throughout RMG Py. RMG Py uses SI units throughout; accordingly, all of the
constants in this module are stored in combinations of meters, seconds,
kilograms, moles, etc.

The constants available are listed below. All values were taken from
`NIST <http://physics.nist.gov/cuu/Constants/index.html>`_.

"""

import math
import cython
import numpy
import quantities as pq

pq.UnitQuantity('kilocalories', pq.cal*1e3, symbol='kcal')
pq.UnitQuantity('kilojoules', pq.J*1e3, symbol='kJ')
pq.UnitQuantity('kilomoles', pq.mol*1e3, symbol='kmol')

################################################################################

#: The Avogadro constant :math:`N_\mathrm{a}`, in :math:`\mathrm{mol^{-1}}`.
Na = 6.02214179e23

#: The Boltzmann constant :math:`k_\mathrm{B}`, in :math:`\mathrm{J/K}`.
kB = 1.3806504e-23

#: The gas law constant :math:`R`, in :math:`\mathrm{J/mol \cdot K}`.
R = 8.314472

#: The Planck constant :math:`h`, in :math:`\mathrm{J \cdot s}`.
h = 6.62606896e-34

#: The speed of light in a vacuum :math:`c`, in :math:`\mathrm{m/s}`.
c = 299792458

#: The mathematical constant :math:`\pi = 3.14159...`
pi = float(math.pi)

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
    A single numeric quantity, with optional units and uncertainty. The
    attributes are:

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

    In order for this class to be efficient for mathematical operations, we
    must minimize the number of unit conversions that are performed. To this
    end, all numeric values and uncertainties are stored in SI units, rather
    than the units they are specified as.
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
        elif isinstance(args, numpy.ndarray):
            # We've been given just an array of numbers
            value = args
        elif isinstance(args, float) or isinstance(args, int):
            # We've been given just a single number
            value = args
        elif isinstance(args, list) or isinstance(args, tuple):
            if len(args) == 1:
                value = args[0]
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

    def __str__(self):
        """
        Return a string representation of the object.
        """
        factor = getConversionFactorFromSI(self.units)
        string = ''
        if self.values is None:
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
        if self.values is None:
            string += '%g' % (self.value * factor)
        else:
            string += '[%s]' % (','.join(['%g' % (v * factor) for v in self.values]))
        if self.units != '' or self.uncertaintyType != '':
            string += ',"%s"' % (self.units)
            if self.uncertaintyType != '':
                string += ',"%s"' % (self.uncertaintyType)
                if not self.isUncertaintyAdditive(): factor = 1.0
                if self.uncertainties is None:
                    string += ',%g' % (self.uncertainty * factor)
                else:
                    string += ',[%s]' % (','.join(['%g' % (u * factor) for u in self.uncertainties]))
            string = '(' + string + ')'
        return string

    def __eq__(self, other):
        """
        Return ``True`` if two quantity objects are equal or ``False`` if not.
        """
        if not isinstance(other, Quantity): return False
        if self.values is not None and self.uncertainties is not None:
            if other.values is None or other.uncertainties is None:
                return False
            else:
                return (self.value == other.value and self.units == other.units and
                    self.uncertaintyType == other.uncertaintyType and self.uncertainty == other.uncertainty and
                    (self.values == other.values).all() and (self.uncertainties == other.uncertainties).all())
        elif self.values is not None and self.uncertainties is None:
            if other.values is None or other.uncertainties is not None:
                return False
            else:
                return (self.value == other.value and self.units == other.units and
                    self.uncertaintyType == other.uncertaintyType and self.uncertainty == other.uncertainty and
                    (self.values == other.values).all())
        else:
            if other.values is not None or other.uncertainties is not None:
                return False
            else:
                return (self.value == other.value and self.units == other.units and
                    self.uncertaintyType == other.uncertaintyType and self.uncertainty == other.uncertainty)

    def __add__(self, other):
        """
        Add two objects together and return the sum as a new object.
        """
        if not isinstance(other, Quantity): raise ValueError('Unexpected type "%s" for other parameter.' % (other.__class__))
        assert self.units == other.units
        assert self.uncertaintyType == other.uncertaintyType == '+|-'

        new = Quantity()
        if self.values is not None and self.uncertainties is not None:
            new.values = self.values + other.values
            new.units = self.units
            new.uncertaintyType = self.uncertaintyType
            new.uncertainties = numpy.sqrt(self.uncertainties**2 + other.uncertainties**2)
        elif self.values is not None and self.uncertainties is None:
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
