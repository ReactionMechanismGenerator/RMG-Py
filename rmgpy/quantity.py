#!/usr/bin/env python
# encoding: utf-8

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2017 Prof. William H. Green (whgreen@mit.edu), 
#   Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)
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
particularly the :class:`Quantity` class for representing physical quantities.
"""

import numpy
import quantities as pq

import rmgpy.constants as constants
from rmgpy.exceptions import QuantityError
import re
from logging import warn

################################################################################

# Explicitly set the default units to SI
pq.set_default_units('si')

# These units are not defined by the quantities package, but occur frequently
# in data handled by RMG, so we define them manually
pq.UnitQuantity('kilocalories', pq.cal*1e3, symbol='kcal')
pq.UnitQuantity('kilojoules', pq.J*1e3, symbol='kJ')
pq.UnitQuantity('kilomoles', pq.mol*1e3, symbol='kmol')
pq.UnitQuantity('molecule', pq.mol/6.02214179e23, symbol='molecule')
pq.UnitQuantity('molecules', pq.mol/6.02214179e23, symbol='molecules')
pq.UnitQuantity('debye', 1.0/(constants.c*1e21)*pq.C*pq.m, symbol='De')

################################################################################

# Units that should not be used in RMG-Py:
NOT_IMPLEMENTED_UNITS = [
                        'degC',
                        'C',
                        'degF',
                        'F',
                        'degR',
                        'R'
                        ]

################################################################################

class Units(object):
    """
    The :class:`Units` class provides a representation of the units of a
    physical quantity. The attributes are:

    =================== ========================================================
    Attribute           Description
    =================== ========================================================
    `units`             A string representation of the units
    =================== ========================================================

    Functions that return the conversion factors to and from SI units are
    provided.
    """
    
    # A dict of conversion factors (to SI) for each of the frequent units
    # Here we also define that cm^-1 is not to be converted to m^-1 (or Hz, J, K, etc.)
    conversionFactors = {'cm^-1': 1.0}
    
    def __init__(self, units=''):
        if units in NOT_IMPLEMENTED_UNITS:
            raise NotImplementedError(
        'The units {} are not yet supported. Please choose SI units.'.format(units)
                                     )
        self.units = units

    def getConversionFactorToSI(self):
        """
        Return the conversion factor for converting a quantity in a given set
        of`units` to the SI equivalent units.
        """
        try:
            # Process several common units manually for speed
            factor = Units.conversionFactors[self.units]
        except KeyError:
            # Fall back to (slow!) quantities package for less common units
            factor = float(pq.Quantity(1.0, self.units).simplified)
            # Cache the conversion factor so we don't ever need to use
            # quantities to compute it again
            Units.conversionFactors[self.units] = factor
        return factor

    def getConversionFactorFromSI(self):
        """
        Return the conversion factor for converting a quantity to a given set
        of `units` from the SI equivalent units.
        """
        return 1.0 / self.getConversionFactorToSI()

    def si_units(self):
        """
        Return the corresponding si units for dimensionality checks.
        """
        return pq.Quantity(1.0, self.units).simplified.units

    def return_compounded_units(self, other):
        """
        Return the compound unit upon multiplication of two ScalarQuantities
        """
        compound_unit = pq.CompoundUnit(self.units + "*" + other.units).name
        return Units.simplify_units(compound_unit[1:-1])  # [1:-1] removes extraneous parenthesis

    def return_divided_units(self, other):
        """
        Return the compound unit upon division of two ScalarQuantities
        """
        divided_unit = pq.CompoundUnit(self.units + "/(" + other.units + ")").name
        return Units.simplify_units(divided_unit[1:-1])  # [1:-1] removes extraneous parenthesis

# Helper functions for simplifying units

    @staticmethod
    def simplify_units(unit_string):
        """
        Return a simplified version of a unit string to remove duplicate entries for a given unit. This will simplify
        things like m**3/m**2 to m, but will not simplify things like N*m to J or things like ft*m.
        """

        original_unit_string = unit_string  # save original string for checking at the end

        # Replace all '**' exponential operators with '^' so regex does not confuse with '*'
        unit_string = re.sub('\*\*', '^', unit_string) + '*'  # extra '*' for pattern matching of last unit

        # Simplify groups inside of parenthesis first
        start, end = Units.find_unit_pattern('\([^\)]*\)', unit_string)  # get indices of the embedded parenthesis
        loop_count = 0
        while start is not None:
            paren_unit = unit_string[start + 1:end - 1] + '*'
            paren_dict = Units.return_unit_dictionary(paren_unit)  # dictionary of units with their respective powers

            # If this term was being divided, invert the powers for every unit and switch to multiplication
            if (unit_string[start - 1] == '/') and (start-1 > 0):
                unit_string = unit_string[:start - 1] + '*' + unit_string[start:]  # switch to multiplication
                for unit in paren_dict.iterkeys():
                    paren_dict[unit] *= -1.0

            # If this term was being raised to a power, adjust the powers accordingly
            if unit_string[end] == '^':
                # Find the exponent
                end_string = unit_string[end:]
                power_start, power_end = Units.find_unit_pattern('\^[^\*,\(,\),/,\^]*', end_string)
                overall_power = float(end_string[power_start + 1:power_end])
                end = end + power_end  # this will allow us to get rid of the overall power later
                for unit in paren_dict.iterkeys():
                    paren_dict[unit] *= overall_power

            # Replace the part in parenthesis with the above simplified string to remove parenthesis
            unit_string = unit_string[:start] + Units.process_unit_dict(paren_dict) + unit_string[end:]

            # Move on to next parenthesis (if any)
            start, end = Units.find_unit_pattern('\([^\)]*\)', unit_string)

            if loop_count > 100:
                warn('Unable to simplify units {}. Maximum iteration for parenthesis exceeded. Original unit string '
                     'returned'.format(original_unit_string))
                return original_unit_string

        # To simplify the string, convert the unit string into a unit-power dictionary
        new_unit_dict = Units.return_unit_dictionary(unit_string)
        new_unit_string = Units.process_unit_dict(new_unit_dict)

        # Check that the dimensions of the unit string was not change (if so this was unsuccessful)
        if pq.CompoundUnit(new_unit_string) == pq.CompoundUnit(original_unit_string):
            return new_unit_string
        else:
            warn('unit string {} could not be simplified. Original unit string returned'.format(original_unit_string))
            return original_unit_string

    @staticmethod
    def return_unit_dictionary(unit_string):
        """
        Simplify a unit string (that does not contain parenthesis) by converting into a dictionary with units as keys and
        simplified powers as values
        """

        # Convert any division groups into multiplication to the inverse power
        start, end = Units.find_unit_pattern('/[^\*, /]*\*', unit_string)  # Get the indices of the division group
        loop_count = 0
        while start is not None:
            dividing_unit = unit_string[start:end].split('^')  # Separate units and their powers
            if len(dividing_unit) == 2:
                power = float(dividing_unit[1][:-1])
                unit_string = unit_string[:start] + '*' + dividing_unit[0][1:] + '^' + str(power * -1.0) + \
                                                    dividing_unit[1][-1] + unit_string[end:]
            elif len(dividing_unit) == 1:  # Then the unit was to the first power
                power = '-1'
                unit_string = unit_string[:start] + '*' + dividing_unit[0][1:-1] + '^' + power + \
                                                    dividing_unit[0][-1] + unit_string[end:]

            # Move on to next division group, if any
            start, end = Units.find_unit_pattern('/[^\*, /]*[\*, /]', unit_string)

            if loop_count > 100:
                warn('Unable to return unit-power dictionary. Maximum loop count exceeded')
                return {}  # This will be dealt with in simplify_units

        # At this point, we have only units (with powers) being multiplied. Split into groups of units^power
        unit_parts = unit_string[:-1].split('*')  # [:-1] gets rid of extraneous '*' that we added
        new_unit_dict = {}

        for parts in unit_parts:
            # Split into units and powers
            part_fragment = parts.split('^')
            if len(part_fragment) == 2:
                unit = part_fragment[0]
                power = part_fragment[1]
            elif len(part_fragment) == 1:  # Then this unit was to the first power
                unit = part_fragment[0]
                power = 1

            # Add units and powers to the unit-power dictionary and update along the way
            if unit in new_unit_dict.keys():
                new_unit_dict[unit] += float(power)
            else:
                new_unit_dict[unit] = float(power)

        return new_unit_dict

    @staticmethod
    def process_unit_dict(new_unit_dict):
        """
        Convert a unit-power dictionary {unit1:power1, ...} into a sting with units raised to their respective powers
        """
        new_unit_string = ''
        for unit, power in new_unit_dict.iteritems():
            if power != 0:
                if power == 1:  # No need for a power
                    new_unit_string += unit + '*'
                else:
                    new_unit_string += unit + '^' + str(power) + '*'
        if new_unit_string == '':
            return new_unit_string
        else:
            return new_unit_string[:-1]  # [:-1] gets rid of extraneous '*' at end

    @staticmethod
    def find_unit_pattern(pattern, string):
        """
        Return the start and end indices of a sub-string in 'string' that matches the regex pattern
        """

        try:
            start = re.search(pattern, string).start()
            end = re.search(pattern, string).end()
        except AttributeError:  # No pattern found, return None for verification
            start = None
            end = None

        return start, end

################################################################################


class ScalarQuantity(Units):
    """
    The :class:`ScalarQuantity` class provides a representation of a scalar
    physical quantity, with optional units and uncertainty information. The
    attributes are:

    =================== ========================================================
    Attribute           Description
    =================== ========================================================
    `value`             The numeric value of the quantity in the given units
    `units`             The units the value was specified in
    `uncertainty`       The numeric uncertainty in the value in the given units (unitless if multiplicative)
    `uncertaintyType`   The type of uncertainty: ``'+|-'`` for additive, ``'*|/'`` for multiplicative
    `value_si`          The numeric value of the quantity in the corresponding SI units
    `uncertainty_si`    The numeric value of the uncertainty in the corresponding SI units (unitless if multiplicative)
    =================== ========================================================

    It is often more convenient to perform computations using SI units instead
    of the given units of the quantity. For this reason, the SI equivalent of
    the `value` attribute can be directly accessed using the `value_si`
    attribute. This value is cached on the :class:`ScalarQuantity` object for
    speed.
    """

    def __init__(self, value, units='', uncertainty=None, uncertaintyType='+|-'):
        Units.__init__(self, units)
        self.value = value
        self.uncertaintyType = uncertaintyType
        self.uncertainty = float(uncertainty) if uncertainty is not None else 0.0

    def __reduce__(self):
        """
        Return a tuple of information used to pickle the scalar quantity.
        """
        return (ScalarQuantity, (self.value, self.units, self.uncertainty, self.uncertaintyType))

    def __str__(self):
        """
        Return a string representation of the scalar quantity.
        """
        result = '{0:g}'.format(self.value)
        if self.uncertainty != 0.0:
            result += ' {0} {1:g}'.format(self.uncertaintyType, self.uncertainty)
        if self.units != '':
            result += ' {0}'.format(self.units)
        return result

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        scalar quantity.
        """
        if self.units == '' and self.uncertainty == 0.0:
            return '{0:g}'.format(self.value)
        else:
            result = '({0:g},{1!r}'.format(self.value, self.units)
            if self.uncertainty != 0.0:
                result += ',{0!r},{1:g}'.format(self.uncertaintyType, self.uncertainty)
            result += ')'
            return result
    
    def __add__(self, other):
        if self.has_same_dimensions_as(other):
            sum_si = float(self.value_si + other.value_si)
            sum_uncertainty = self.additive_uncertainty(other)

            return ScalarQuantity(sum_si*self.getConversionFactorFromSI(), self.units, uncertainty=sum_uncertainty,
                                  uncertaintyType='+|-')

        else:
            raise QuantityError('Cannot add items with units of "{0}" and "{1}"'.format(self.units, other.units))

    def __sub__(self, other):
        if self.has_same_dimensions_as(other):
            result_si = float(self.value_si - other.value_si)
            result_uncertainty = self.additive_uncertainty(other)

            return ScalarQuantity(result_si*self.getConversionFactorFromSI(), self.units,
                                  uncertainty=result_uncertainty, uncertaintyType='+|-')

        else:
            raise QuantityError('Cannot subtract items with units of "{0}" and "{1}"'.format(self.units, other.units))

    def __mul__(self, other):
        product = self.value * other.value
        product_uncertainty = self.multiplicative_uncertainty(other)
        return ScalarQuantity(product, self.return_compounded_units(other), uncertainty=product_uncertainty,
                              uncertaintyType='*|/')

    def __div__(self, other):
        result = float(self.value) / float(other.value)
        result_uncertainty = self.multiplicative_uncertainty(other)
        return ScalarQuantity(result, self.return_divided_units(other), uncertainty=result_uncertainty,
                              uncertaintyType='*|/')

    def __richcmp__(self, other, cmp_type):
        if self.has_same_dimensions_as(other):
            if cmp_type == 0:
                return self.value_si < other.value_si
            if cmp_type == 1:
                return self.value_si <= other.value_si
            if cmp_type == 2:
                return self.value_si == other.value_si
            if cmp_type == 3:
                return self.value_si != other.value_si
            if cmp_type == 4:
                return self.value_si > other.value_si
            if cmp_type == 5:
                return self.value_si >= other.value_si
            else:
                raise Exception('Unsupported operation between ScalarQuantity types')
        else:
            raise QuantityError('Cannot compare Scalar Quantities with dimensions of {0} and {1}. Original units of '
                                'the objects are {2} and {3}'.format(self.si_units(), other.si_units(),
                                                                     self.units, other.units))

    def has_same_dimensions_as(self, other):
        """
        Returns True if two ScalarQuantity objects have the same dimensions
        """
        return self.si_units() == other.si_units()

    def additive_uncertainty(self, other):
        """
        Returns the value of the resultant uncertainty for addition and subtraction operations
        """

        if self.isUncertaintyAdditive():
            uncertainty_a = self.uncertainty_si
        else:
            uncertainty_a = self.uncertainty * self.value_si

        if other.isUncertaintyAdditive():
            uncertainty_b = other.uncertainty_si
        else:
            uncertainty_b = other.uncertainty * other.value_si

        return ((uncertainty_a ** 2 + uncertainty_b ** 2) ** 0.5) * self.getConversionFactorFromSI()

    def multiplicative_uncertainty(self, other):
        """
        Returns the value of the resultant uncertainty for multiplication and division operations
        """

        if self.isUncertaintyMultiplicative():
            uncertainty_a = self.uncertainty
        else:
            uncertainty_a = self.uncertainty_si/self.value_si

        if other.isUncertaintyMultiplicative():
            uncertainty_b = other.uncertainty
        else:
            uncertainty_b = other.uncertainty_si/other.value_si

        return (uncertainty_a ** 2 + uncertainty_b ** 2)**0.5

    def copy(self):
        """
        Return a copy of the quantity.
        """
        return ScalarQuantity(self.value, self.units, self.uncertainty, self.uncertaintyType)

    def getValue(self):
        """
        The numeric value of the quantity, in the given units
        """
        return self.value_si * self.getConversionFactorFromSI()
    def setValue(self, v):
        self.value_si = float(v) * self.getConversionFactorToSI()
    value = property(getValue, setValue)
    
    def getUncertainty(self):
        """
        The numeric value of the uncertainty, in the given units if additive, or no units if multiplicative.
        """
        if self.isUncertaintyAdditive():
            return self.uncertainty_si * self.getConversionFactorFromSI()
        else:
            return self.uncertainty_si
    def setUncertainty(self, v):
        if self.isUncertaintyAdditive():
            self.uncertainty_si = float(v) * self.getConversionFactorToSI()
        else:
            self.uncertainty_si = float(v)
    uncertainty = property(getUncertainty, setUncertainty)
    
    def getUncertaintyType(self):
        """
        The type of uncertainty: ``'+|-'`` for additive, ``'*|/'`` for multiplicative
        """
        return self._uncertaintyType
    def setUncertaintyType(self, v):
        """
        Check the uncertainty type is valid, then set it.
        """
        if v not in ['+|-','*|/']:
            raise QuantityError('Unexpected uncertainty type "{0}"; valid values are "+|-" and "*|/".'.format(v))
        self._uncertaintyType = v

    uncertaintyType = property(getUncertaintyType, setUncertaintyType)
    
    
    
    def equals(self, quantity):
        """
        Return ``True`` if the everything in a quantity object matches
        the parameters in this object.  If there are lists of values or uncertainties,
        each item in the list must be matching and in the same order.
        Otherwise, return ``False``
        (Originally intended to return warning if units capitalization was
        different, however, Quantity object only parses units matching in case, so
        this will not be a problem.)
        """

        def approx_equal(x, y, atol = .01):
            """
            Returns true if two float/double values are approximately equal
            within a relative error of 1% or under a user specific absolute tolerance.
            """
            return abs(x-y) <= 1e-2*abs(x) or abs(x-y) <= 1e-2*abs(y) or abs(x-y) <= atol

        if isinstance(quantity, ScalarQuantity):
            if (self.uncertaintyType == quantity.uncertaintyType and
                approx_equal(self.uncertainty * self.getConversionFactorToSI(), quantity.uncertainty * quantity.getConversionFactorToSI()) and
                self.units == quantity.units):

                if self.units == "kcal/mol":
                    # set absolute tolerance to .01 kcal/mol = 42 J/mol
                    atol = 42
                else:
                    # for other units, set it to .01
                    atol = .01

                if not approx_equal(self.value_si, quantity.value_si, atol):
                    return False

                return True

        return False

    def isUncertaintyAdditive(self):
        """
        Return ``True`` if the uncertainty is specified in additive format
        and ``False`` otherwise.
        """
        return self._uncertaintyType == '+|-'

    def isUncertaintyMultiplicative(self):
        """
        Return ``True`` if the uncertainty is specified in multiplicative 
        format and ``False`` otherwise.
        """
        return self._uncertaintyType == '*|/'

################################################################################

class ArrayQuantity(Units):
    """
    The :class:`ArrayQuantity` class provides a representation of an array of
    physical quantity values, with optional units and uncertainty information. 
    The attributes are:

    =================== ========================================================
    Attribute           Description
    =================== ========================================================
    `value`             The numeric value of the quantity in the given units
    `units`             The units the value was specified in
    `uncertainty`       The numeric uncertainty in the value (unitless if multiplicative)
    `uncertaintyType`   The type of uncertainty: ``'+|-'`` for additive, ``'*|/'`` for multiplicative
    `value_si`          The numeric value of the quantity in the corresponding SI units
    `uncertainty_si`    The numeric value of the uncertainty in the corresponding SI units (unitless if multiplicative)
    =================== ========================================================

    It is often more convenient to perform computations using SI units instead
    of the given units of the quantity. For this reason, the SI equivalent of
    the `value` attribute can be directly accessed using the `value_si`
    attribute. This value is cached on the :class:`ArrayQuantity` object for 
    speed.
    """

    def __init__(self, value, units='', uncertainty=None, uncertaintyType='+|-'):
        Units.__init__(self, units)
        self.value = value
        self.uncertaintyType = uncertaintyType
        if uncertainty is None:
            self.uncertainty = numpy.zeros_like(self.value)
        elif isinstance(uncertainty, (int,float)):
            self.uncertainty = numpy.ones_like(self.value) * uncertainty
        else:
            uncertainty = numpy.array(uncertainty)
            if uncertainty.ndim != self.value.ndim:
                raise QuantityError('The given uncertainty has {0:d} dimensions, while the given value has {1:d} dimensions.'.format(uncertainty.ndim, self.value.ndim))
            for i in range(self.value.ndim):
                if self.value.shape[i] != uncertainty.shape[i]:
                    raise QuantityError('Dimension {0:d} has {1:d} elements for the given value, but {2:d} elements for the given uncertainty.'.format(i, self.value.shape[i], uncertainty.shape[i]))
            else:
                self.uncertainty = uncertainty
    
    def __reduce__(self):
        """
        Return a tuple of information used to pickle the array quantity.
        """
        return (ArrayQuantity, (self.value, self.units, self.uncertainty, self.uncertaintyType))

    def __str__(self):
        """
        Return a string representation of the array quantity.
        """
        if self.value.ndim == 1:
            value = '[{0}]'.format(','.join(['{0:g}'.format(float(v)) for v in self.value]))
        elif self.value.ndim == 2:
            value = []
            for i in range(self.value.shape[0]):
                value.append('[{0}]'.format(','.join(['{0:g}'.format(float(self.value[i,j])) for j in range(self.value.shape[1])])))
            value = '[{0}]'.format(','.join(value))

        if self.uncertainty.ndim == 1:
            uncertainty = '[{0}]'.format(','.join(['{0:g}'.format(float(v)) for v in self.uncertainty]))
        elif self.uncertainty.ndim == 2:
            uncertainty = []
            for i in range(self.uncertainty.shape[0]):
                uncertainty.append('[{0}]'.format(','.join(['{0:g}'.format(float(self.uncertainty[i,j])) for j in range(self.uncertainty.shape[1])])))
            uncertainty = '[{0}]'.format(','.join(uncertainty))
        
        result = '{0}'.format(value)
        if any(self.uncertainty != 0.0):
            result += ' {0} {1}'.format(self.uncertaintyType, uncertainty)
        if self.units != '':
            result += ' {0}'.format(self.units)
        return result

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        array quantity.
        """
        if self.value.ndim == 1:
            value = '[{0}]'.format(','.join(['{0:g}'.format(float(v)) for v in self.value]))
        elif self.value.ndim == 2:
            value = []
            for i in range(self.value.shape[0]):
                value.append('[{0}]'.format(','.join(['{0:g}'.format(float(self.value[i,j])) for j in range(self.value.shape[1])])))
            value = '[{0}]'.format(','.join(value))

        if self.uncertainty.ndim == 1:
            uncertainty = '[{0}]'.format(','.join(['{0:g}'.format(float(v)) for v in self.uncertainty]))
        elif self.uncertainty.ndim == 2:
            uncertainty = []
            for i in range(self.uncertainty.shape[0]):
                uncertainty.append('[{0}]'.format(','.join(['{0:g}'.format(float(self.uncertainty[i,j])) for j in range(self.uncertainty.shape[1])])))
            uncertainty = '[{0}]'.format(','.join(uncertainty))

        if self.units == '' and not numpy.any(self.uncertainty != 0.0):
            return '{0}'.format(value)
        else:
            result = '({0},{1!r}'.format(value, self.units)
            if numpy.any(self.uncertainty != 0.0):
                result += ',{0!r},{1}'.format(self.uncertaintyType, uncertainty)
            result += ')'
            return result

    def copy(self):
        """
        Return a copy of the quantity.
        """
        return ArrayQuantity(self.value.copy(), self.units, self.uncertainty.copy(), self.uncertaintyType)

    def getValue(self):
        """
        The numeric value of the array quantity, in the given units.
        """
        return self.value_si * self.getConversionFactorFromSI()
    def setValue(self, v):
        self.value_si = numpy.array(v) * self.getConversionFactorToSI()
    value = property(getValue, setValue)
    
    def getUncertainty(self):
        """
        The numeric value of the uncertainty, in the given units if additive, or no units if multiplicative.
        """
        if self.isUncertaintyAdditive():
            return self.uncertainty_si * self.getConversionFactorFromSI()
        else:
            return self.uncertainty_si
    def setUncertainty(self, v):
        if self.isUncertaintyAdditive():
            self.uncertainty_si = numpy.array(v) * self.getConversionFactorToSI()
        else:
            self.uncertainty_si = numpy.array(v)
    uncertainty = property(getUncertainty, setUncertainty)
    
    def getUncertaintyType(self):
        """
        The type of uncertainty: ``'+|-'`` for additive, ``'*|/'`` for multiplicative
        """
        return self._uncertaintyType
    def setUncertaintyType(self, v):
        """
        Check the uncertainty type is valid, then set it.
        
        If you set the uncertainty then change the type, we have no idea what to do with 
        the units. This ensures you set the type first.
        """
        if v not in ['+|-','*|/']:
            raise QuantityError('Unexpected uncertainty type "{0}"; valid values are "+|-" and "*|/".'.format(v))
        self._uncertaintyType = v

    uncertaintyType = property(getUncertaintyType, setUncertaintyType)

    def equals(self, quantity):
        """
        Return ``True`` if the everything in a quantity object matches
        the parameters in this object.  If there are lists of values or uncertainties,
        each item in the list must be matching and in the same order.
        Otherwise, return ``False``
        (Originally intended to return warning if units capitalization was
        different, however, Quantity object only parses units matching in case, so
        this will not be a problem.)
        """

        def approx_equal(x, y, atol = .01):
            """
            Returns true if two float/double values are approximately equal
            within a relative error of 1% or under a user specific absolute tolerance.
            """
            return abs(x-y) <= 1e-2*abs(x) or abs(x-y) <= 1e-2*abs(y) or abs(x-y) <= atol

        if isinstance(quantity, ArrayQuantity):
            if (self.uncertaintyType == quantity.uncertaintyType and self.units == quantity.units):

                if self.units == "kcal/mol":
                    # set absolute tolerance to .01 kcal/mol = 42 J/mol
                    atol = 42
                else:
                    # for other units, set it to .01
                    atol = .01

                if self.value.ndim != quantity.value.ndim:
                    return False
                for i in range(self.value.ndim):
                    if self.value.shape[i] != quantity.value.shape[i]:
                        return False
                for v1, v2 in zip(self.value.flat, quantity.value.flat):
                    if not approx_equal(v1, v2, atol):
                        return False
                
                if self.uncertainty.ndim != quantity.uncertainty.ndim:
                    return False
                for i in range(self.uncertainty.ndim):
                    if self.uncertainty.shape[i] != quantity.uncertainty.shape[i]:
                        return False
                for v1, v2 in zip(self.uncertainty.flat, quantity.uncertainty.flat):
                    if not approx_equal(v1, v2, atol):
                        return False

                return True

        return False

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

def Quantity(*args, **kwargs):
    """
    Create a :class:`ScalarQuantity` or :class:`ArrayQuantity` object for a
    given physical quantity. The physical quantity can be specified in several
    ways:
    
    * A scalar-like or array-like value (for a dimensionless quantity)
    
    * An array of arguments (including keyword arguments) giving some or all of
      the `value`, `units`, `uncertainty`, and/or `uncertaintyType`.
    
    * A tuple of the form ``(value,)``, ``(value,units)``, 
      ``(value,units,uncertainty)``, or 
      ``(value,units,uncertaintyType,uncertainty)``
    
    * An existing :class:`ScalarQuantity` or :class:`ArrayQuantity` object, for
      which a copy is made
    
    """
    # Initialize attributes
    value = None
    units = ''
    uncertaintyType = '+|-'
    uncertainty = None
    
    if len(args) == 1 and len(kwargs) == 0 and args[0] is None:
        return None
    
    # Unpack args if necessary
    if isinstance(args, tuple) and len(args) == 1 and isinstance(args[0], tuple):
        args = args[0]
        
    # Process args    
    Nargs = len(args)
    if Nargs == 1 and isinstance(args[0], (ScalarQuantity,ArrayQuantity)):
        # We were given another quantity object, so make a (shallow) copy of it
        other = args[0]
        value = other.value
        units = other.units
        uncertaintyType = other.uncertaintyType
        uncertainty = other.uncertainty
    elif Nargs == 1:
        # If one parameter is given, it should be a single value
        value, = args
    elif Nargs == 2:
        # If two parameters are given, it should be a value and units
        value, units = args
    elif Nargs == 3:
        # If three parameters are given, it should be a value, units and uncertainty
        value, units, uncertainty = args
    elif Nargs == 4:
        # If four parameters are given, it should be a value, units, uncertainty type, and uncertainty
        value, units, uncertaintyType, uncertainty = args
    elif Nargs != 0:
        raise QuantityError('Invalid parameters {0!r} passed to ArrayQuantity.__init__() method.'.format(args))
    
    # Process kwargs
    for k, v in kwargs.items():
        if k == 'value':
            if len(args) >= 1:
                raise QuantityError('Multiple values for argument {0} passed to ArrayQuantity.__init__() method.'.format(k))
            else:
                value = v
        elif k == 'units':
            if len(args) >= 2:
                raise QuantityError('Multiple values for argument {0} passed to ArrayQuantity.__init__() method.'.format(k))
            else:
                units = v
        elif k == 'uncertainty':
            if len(args) >= 3:
                raise QuantityError('Multiple values for argument {0} passed to ArrayQuantity.__init__() method.'.format(k))
            else:
                uncertainty = v
        elif k == 'uncertaintyType':
            if len(args) >= 4:
                raise QuantityError('Multiple values for argument {0} passed to ArrayQuantity.__init__() method.'.format(k))
            else:
                uncertaintyType = v
        else:
            raise QuantityError('Invalid keyword argument {0} passed to ArrayQuantity.__init__() method.'.format(k))
    
    # Process units and uncertainty type parameters
    if uncertaintyType not in ['+|-', '*|/']:
        raise QuantityError('Unexpected uncertainty type "{0}"; valid values are "+|-" and "*|/".'.format(uncertaintyType))

    if isinstance(value, (list,tuple,numpy.ndarray)):
        return ArrayQuantity(value, units, uncertainty, uncertaintyType)
    
    try:
        value = float(value)
    except TypeError:
        return ArrayQuantity(value, units, uncertainty, uncertaintyType)
    
    uncertainty = 0.0 if uncertainty is None else float(uncertainty)
    return ScalarQuantity(value, units, uncertainty, uncertaintyType)

################################################################################

class UnitType:
    """
    The :class:`UnitType` class represents a factory for producing
    :class:`ScalarQuantity` or :class:`ArrayQuantity` objects of a given unit
    type, e.g. time, volume, etc.
    """

    def __init__(self, units, commonUnits=None, extraDimensionality=None):
        self.units = units
        self.dimensionality = pq.Quantity(1.0, units).simplified.dimensionality
        self.commonUnits = commonUnits or []
        self.extraDimensionality = {}
        if extraDimensionality:
            for unit, factor in extraDimensionality.items():
                self.extraDimensionality[pq.Quantity(1.0, unit).simplified.dimensionality] = factor
        
    def __call__(self, *args, **kwargs):
        # Make a ScalarQuantity or ArrayQuantity object out of the given parameter
        quantity = Quantity(*args, **kwargs)
        if quantity is None:
            return quantity
        
        units = quantity.units
        
        # If the units are in the common units, then we can do the conversion
        # very quickly and avoid the slow calls to the quantities package
        if units == self.units or units in self.commonUnits:
            return quantity
        
        # Check that the units are consistent with this unit type
        # This uses the quantities package (slow!)
        units = pq.Quantity(1.0, units)
        dimensionality = units.simplified.dimensionality
        if dimensionality == self.dimensionality:
            pass
        elif dimensionality in self.extraDimensionality:
            quantity.value_si *= self.extraDimensionality[dimensionality]
            quantity.units = self.units
        else:
            raise QuantityError('Invalid units {0!r}. Try common units: {1}'.format(quantity.units,self.commonUnits))
        
        # Return the Quantity or ArrayQuantity object object
        return quantity

Acceleration = UnitType('m/s^2')

Area = UnitType('m^2')

Concentration = UnitType('mol/m^3')

Dimensionless = UnitType('')

DipoleMoment = UnitType('C*m', extraDimensionality={
    'De': 1.0 / (1.0e21 * constants.c), 
})

"We have to allow 'energies' to be created in units of Kelvins, because Chemkin does so"
Energy = Enthalpy = FreeEnergy = UnitType('J/mol',
    commonUnits=['kJ/mol', 'cal/mol', 'kcal/mol'],
    extraDimensionality={'K': constants.R },
)


Entropy = HeatCapacity = UnitType('J/(mol*K)', commonUnits=['kJ/(mol*K)', 'cal/(mol*K)', 'kcal/(mol*K)'])

Flux = UnitType('mol/(m^2*s)')

Frequency = UnitType('cm^-1', extraDimensionality={
    's^-1': 1.0 / (constants.c * 100.),
    'Hz': 1.0 / (constants.c * 100.),
    'J': 1.0 / (constants.h * constants.c * 100.),
    'K': constants.kB / (constants.h * constants.c * 100.), 
})

Force = UnitType('N')

Inertia = UnitType('kg*m^2')

Length = UnitType('m')

Mass = UnitType('amu', extraDimensionality={'kg/mol': 1000.*constants.amu})

Momentum = UnitType('kg*m/s^2')

Power = UnitType('W')

Pressure = UnitType('Pa', commonUnits=['bar', 'atm', 'torr', 'psi', 'mbar'])

Temperature = UnitType('K', commonUnits=[])

Time = UnitType('s')

Velocity = UnitType('m/s')

Volume = UnitType('m^3')

# Polarizability = UnitType('C*m^2*V^-1')
"""
What's called Polarizability in the transport properties is in fact a polarizability volume,
which is related by $4*\pi*\epsilon_0$ where $\epsilon_0$ is the permittivity of free space.
Rather than mess around with conversions, I suggest we just use "Volume" as the units for 
what we call 'polarizability'. Chemkin expects it in Angstrom^3. We'll store it in m^3.
"""

# RateCoefficient is handled as a special case since it can take various
# units depending on the reaction order
RATECOEFFICIENT_CONVERSION_FACTORS = {
    (1.0/pq.s).dimensionality: 1.0,              
    (pq.m**3/pq.s).dimensionality: 1.0,              
    (pq.m**6/pq.s).dimensionality: 1.0,              
    (pq.m**9/pq.s).dimensionality: 1.0,              
    (pq.m**3/(pq.mol*pq.s)).dimensionality: 1.0,              
    (pq.m**6/(pq.mol**2*pq.s)).dimensionality: 1.0,              
    (pq.m**9/(pq.mol**3*pq.s)).dimensionality: 1.0,              
}
RATECOEFFICIENT_COMMON_UNITS = ['s^-1', 'm^3/(mol*s)', 'cm^3/(mol*s)', 'm^3/(molecule*s)', 'cm^3/(molecule*s)']
def RateCoefficient(*args, **kwargs):
    # Make a ScalarQuantity or ArrayQuantity object out of the given parameter
    quantity = Quantity(*args, **kwargs)
    if quantity is None:
        return quantity
    
    units = quantity.units
        
    # If the units are in the common units, then we can do the conversion
    # very quickly and avoid the slow calls to the quantities package
    if units in RATECOEFFICIENT_COMMON_UNITS:
        return quantity

    dimensionality = pq.Quantity(1.0, quantity.units).simplified.dimensionality
    try:
        factor = RATECOEFFICIENT_CONVERSION_FACTORS[dimensionality]
        quantity.value_si *= factor
    except KeyError:
        raise QuantityError('Invalid units {0!r}. Common units are {1}'.format(quantity.units,RATECOEFFICIENT_COMMON_UNITS))

    # Return the Quantity or ArrayQuantity object object
    return quantity
