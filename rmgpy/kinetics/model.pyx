# cython: embedsignature=True, cdivision=True

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2009 Prof. William H. Green (whgreen@mit.edu) and the
#   RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the "Software"),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

"""
This module contains base classes that represent various rate coefficient
models.
"""

import numpy
import quantities as pq

import rmgpy.quantity as quantity

################################################################################

cpdef str getRateCoefficientUnitsFromReactionOrder(order):
    """
    Given a reaction `order`, return the corresponding SI units of the rate
    coefficient. These are the units that rate coefficients are stored in
    internally, as well as the units of the rate coefficient obtained using
    the ``simplified`` attribute of a :class:`Quantity` object that represents
    a rate coefficient. Raises a :class:`ValueError` if the units could not be
    determined.
    """
    if order == 0: 
        kunits = 'mol/(m^3*s)'
    elif order == 1:
        kunits = 's^-1'
    elif order == 2:
        kunits = 'm^3/(mol*s)'
    elif order == 3:
        kunits = 'm^6/(mol^2*s)'
    elif order == 4:
        kunits = 'm^9/(mol^3*s)'
    else:
        raise ValueError('Invalid reaction order {0}.'.format(order))
    return kunits

cpdef int getReactionOrderFromRateCoefficientUnits(kunits) except -1:
    """
    Given a set of rate coefficient units `kunits`, return the corresponding
    reaction order. Raises a :class:`ValueError` if the reaction order could 
    not be determined.
    """
    import quantities as pq
    dimensionality = pq.Quantity(1.0, kunits).simplified.dimensionality
    order = -1
    if len(dimensionality) == 3 and pq.mol in dimensionality and pq.m in dimensionality and pq.s in dimensionality:
        if dimensionality[pq.s] == -1 and dimensionality[pq.m] == -3 and dimensionality[pq.mol] == 1:
            order = 0
        elif dimensionality[pq.s] == -1 and dimensionality[pq.m] == 3 and dimensionality[pq.mol] == -1:
            order = 2
        elif dimensionality[pq.s] == -1 and dimensionality[pq.m] == 6 and dimensionality[pq.mol] == -2:
            order = 3
        elif dimensionality[pq.s] == -1 and dimensionality[pq.m] == 9 and dimensionality[pq.mol] == -3:
            order = 4
    elif len(dimensionality) == 1 and pq.s in dimensionality:
        if dimensionality[pq.s] == -1:
            order = 1
    if order == -1:
        raise ValueError('Invalid rate coefficient units "{0}".'.format(str(dimensionality)))
    return order

################################################################################

cdef class KineticsModel:
    """
    A base class for chemical kinetics models, containing several attributes
    common to all models:
    
    =============== ============================================================
    Attribute       Description
    =============== ============================================================
    `Tmin`          The minimum temperature at which the model is valid, or zero if unknown or undefined
    `Tmax`          The maximum temperature at which the model is valid, or zero if unknown or undefined
    `comment`       Information about the model (e.g. its source)
    =============== ============================================================

    """
    
    def __init__(self, Tmin=None, Tmax=None, comment=''):
        self.Tmin = Tmin
        self.Tmax = Tmax
        self.comment = comment
        
    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        KineticsModel object.
        """
        return 'KineticsModel(Tmin={0!r}, Tmax={1!r}, comment="""{2}""")'.format(self.Tmin, self.Tmax, self.comment)

    def __reduce__(self):
        """
        A helper function used when pickling a KineticsModel object.
        """
        return (KineticsModel, (self.Tmin, self.Tmax, self.comment))

    property Tmin:
        """The minimum temperature at which the model is valid, or ``None`` if not defined."""
        def __get__(self):
            return self._Tmin
        def __set__(self, value):
            self._Tmin = quantity.Temperature(value)

    property Tmax:
        """The maximum temperature at which the model is valid, or ``None`` if not defined."""
        def __get__(self):
            return self._Tmax
        def __set__(self, value):
            self._Tmax = quantity.Temperature(value)

    cpdef bint isPressureDependent(self) except -2:
        """
        Return ``False`` since, by default, all objects derived from KineticsModel
        represent pressure-independent kinetics.
        """
        return False

    cpdef bint isTemperatureValid(self, double T) except -2:
        """
        Return ``True`` if the temperature `T` in K is within the valid
        temperature range of the kinetic data, or ``False`` if not. If
        the minimum and maximum temperature are not defined, ``True`` is 
        returned.
        """
        return (self.Tmin is None or self._Tmin.value_si <= T) and (self.Tmax is None or T <= self._Tmax.value_si)

    cpdef double getRateCoefficient(self, double T, double P=0.0) except -1:
        """
        Return the value of the rate coefficient :math:`k(T)` in units of cm^3,
        mol, and s at the specified temperature `T` in K. This method must be
        overloaded in the derived class.
        """
        raise NotImplementedError('Unexpected call to KineticsModel.getRateCoefficient(); you should be using a class derived from KineticsModel.')

################################################################################

cdef class PDepKineticsModel(KineticsModel):
    """
    A base class for chemical kinetics models that depend on both temperature
    and pressure, containing several attributes common to all such models:
    
    =============== ============================================================
    Attribute       Description
    =============== ============================================================
    `Tmin`          The minimum temperature at which the model is valid, or zero if unknown or undefined
    `Tmax`          The maximum temperature at which the model is valid, or zero if unknown or undefined
    `Pmin`          The minimum pressure at which the model is valid, or zero if unknown or undefined
    `Pmax`          The maximum pressure at which the model is valid, or zero if unknown or undefined
    `efficiencies`  A dict associating chemical species with associated efficiencies
    `highPlimit`    The high-pressure limit kinetics (optional)
    `comment`       Information about the model (e.g. its source)
    =============== ============================================================

    """
    
    def __init__(self, Tmin=None, Tmax=None, Pmin=None, Pmax=None, efficiencies=None, highPlimit=None, comment=''):
        KineticsModel.__init__(self, Tmin, Tmax, comment)
        self.Pmin = Pmin
        self.Pmax = Pmax
        self.efficiencies = efficiencies or {}
        self.highPlimit = highPlimit
        
    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        PDepKineticsModel object.
        """
        return 'PDepKineticsModel(Tmin={0!r}, Tmax={1!r}, Pmin={2!r}, Pmax={3!r}, efficiencies={4!r}, highPlimit={5!r}, comment="""{6}""")'.format(self.Tmin, self.Tmax, self.Pmin, self.Pmax, self.efficiencies, self.highPlimit, self.comment)

    def __reduce__(self):
        """
        A helper function used when pickling a PDepKineticsModel object.
        """
        return (PDepKineticsModel, (self.Tmin, self.Tmax, self.Pmin, self.Pmax, self.efficiencies, self.highPlimit, self.comment))

    property Pmin:
        """The minimum pressure at which the model is valid, or ``None`` if not defined."""
        def __get__(self):
            return self._Pmin
        def __set__(self, value):
            self._Pmin = quantity.Pressure(value)

    property Pmax:
        """The maximum pressure at which the model is valid, or ``None`` if not defined."""
        def __get__(self):
            return self._Pmax
        def __set__(self, value):
            self._Pmax = quantity.Pressure(value)

    cpdef bint isPressureDependent(self) except -2:
        """
        Return ``True`` since all objects derived from PDepKineticsModel
        represent pressure-dependent kinetics.
        """
        return True

    cpdef bint isPressureValid(self, double P) except -2:
        """
        Return ``True`` if the pressure `P` in Pa is within the valid
        pressure range of the kinetic data, or ``False`` if not. If
        the minimum and maximum pressure are not defined, ``True`` is 
        returned.
        """
        return (self.Pmin is None or self._Pmin.value_si <= P) and (self.Pmax is None or P <= self._Pmax.value_si)
    
    cpdef double getEffectivePressure(self, double P, list species, numpy.ndarray fractions) except -1:
        """
        Return the effective pressure in Pa for a system at a given pressure
        `P` in bar composed of the given list of `species` with the given
        `fractions`.
        """
        cdef numpy.ndarray[numpy.float64_t,ndim=1] _fractions
        cdef double Peff, frac, eff, total_frac, eff_frac
        cdef int i
        
        assert len(species) == len(fractions)
        
        _fractions = fractions
        
        # We probably have fewer efficiencies than we do fractions, so 
        # iterating over the species with efficiencies is faster
        Peff = 0.0
        eff_frac = 0.0
        for spec, eff in self.efficiencies.items():
            try:
                i = species.index(spec)
            except ValueError:
                # Species not in list of fractions, so assume fraction of zero
                # and skip to the next species
                continue
            
            frac = _fractions[i]
            Peff += eff * frac
            eff_frac += frac
        
        # For the species with no efficiency data, assume an efficiency of 
        # unity and add to the calculation of the effective pressure
        total_frac = numpy.sum(_fractions)
        Peff += (total_frac - eff_frac) * 1.0
        
        # Don't forget to include the actual pressure and scale by the total
        # fraction (in case fractions is not normalized)
        Peff *= P / total_frac
        
        return Peff
        
    cpdef double getRateCoefficient(self, double T, double P=0.0) except -1:
        """
        Return the value of the rate coefficient :math:`k(T)` in units of cm^3,
        mol, and s at the specified temperature `T` in K and pressure `P` in
        Pa. If you wish to consider collision efficiencies, then you should
        first use :meth:`getEffectivePressure()` to compute the effective
        pressure, and pass that value as the pressure to this method. This
        method must be overloaded in the derived class.
        """
        raise NotImplementedError('Unexpected call to PDepKineticsModel.getRateCoefficient(); you should be using a class derived from PDepKineticsModel.')

################################################################################

cdef class TunnelingModel:
    """
    A base class for models of quantum mechanical tunneling through a reaction
    barrier. The attributes are:
    
    =============== ============================================================
    Attribute       Description
    =============== ============================================================
    `frequency`     The negative frequency along the reaction coordinate
    =============== ============================================================

    """

    def __init__(self, frequency=None):
        self.frequency = frequency

    property frequency:
        """The negative frequency along the reaction coordinate."""
        def __get__(self):
            return self._frequency
        def __set__(self, value):
            self._frequency = quantity.Frequency(value)

    cpdef double calculateTunnelingFactor(self, double T) except -100000000:
        """
        Calculate and return the value of the tunneling correction for
        the reaction at the temperature `T` in K.
        """
        raise NotImplementedError('Unexpected call to Tunneling.calculateTunnelingFactor(); you should be using a class derived from TunnelingModel.')

    cpdef numpy.ndarray calculateTunnelingFunction(self, numpy.ndarray Elist):
        """
        Calculate and return the value of the tunneling correction for
        the reaction at the energies `Elist` in J/mol.
        """
        raise NotImplementedError('Unexpected call to Tunneling.calculateTunnelingFunction(); you should be using a class derived from TunnelingModel.')
