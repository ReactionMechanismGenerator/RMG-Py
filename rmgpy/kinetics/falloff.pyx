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
This module contains classes representing pressure-dependent kinetics models
of "standard" falloff.
"""

import numpy
from libc.math cimport exp, log, log10

cimport rmgpy.constants as constants
import rmgpy.quantity as quantity

################################################################################

cdef class ThirdBody(PDepKineticsModel):
    """
    A kinetic model of a phenomenological rate coefficient :math:`k(T, P)`
    using third-body kinetics. The attributes are:

    =================== ========================================================
    Attribute           Description
    =================== ========================================================
    `arrheniusLow`         The Arrhenius kinetics at the low-pressure limit
    `Tmin`              The minimum temperature at which the model is valid, or zero if unknown or undefined
    `Tmax`              The maximum temperature at which the model is valid, or zero if unknown or undefined
    `Pmin`              The minimum pressure at which the model is valid, or zero if unknown or undefined
    `Pmax`              The maximum pressure at which the model is valid, or zero if unknown or undefined
    `efficiencies`      A dict associating chemical species with associated efficiencies
    `comment`           Information about the model (e.g. its source)
    =================== ========================================================
    
    """

    def __init__(self, arrheniusLow=None, Tmin=None, Tmax=None, Pmin=None, Pmax=None, efficiencies=None, comment=''):
        PDepKineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, efficiencies=efficiencies, comment=comment)
        self.arrheniusLow = arrheniusLow
        
    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        ThirdBody object.
        """
        string = 'ThirdBody(arrheniusLow={0!r}'.format(self.arrheniusLow)
        if self.efficiencies is not None: string += ', efficiencies={0!r}'.format(self.efficiencies)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.Pmin is not None: string += ', Pmin={0!r}'.format(self.Pmin)
        if self.Pmax is not None: string += ', Pmax={0!r}'.format(self.Pmax)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling a ThirdBody object.
        """
        return (ThirdBody, (self.arrheniusLow, self.Tmin, self.Tmax, self.Pmin, self.Pmax, self.efficiencies, self.comment))

    cpdef double getRateCoefficient(self, double T, double P=0.0) except -1:
        """
        Return the value of the rate coefficient :math:`k(T)` in units of m^3,
        mol, and s at the specified temperature `T` in K and pressure `P` in
        Pa. If you wish to consider collision efficiencies, then you should
        first use :meth:`getEffectivePressure()` to compute the effective
        pressure, and pass that value as the pressure to this method.
        """
        cdef double C, k0
                
        C = P / constants.R / T     # bath gas concentration in mol/m^3
        k0 = self.arrheniusLow.getRateCoefficient(T)
        
        return k0 * C

    cpdef bint isIdenticalTo(self, KineticsModel otherKinetics) except -2:
        """
        Checks to see if kinetics matches that of other kinetics and returns ``True``
        if coeffs, kunits, Tmin,
        """
        if not isinstance(otherKinetics, Troe):
            return False
        if not KineticsModel.isIdenticalTo(self,otherKinetics):
            return False
        if not self.arrheniusLow.isIdenticalTo(otherKinetics.arrheniusLow):
            return False
        
        return True
    
    cpdef changeRate(self, double factor):
        """
        Changes kinetics rate by a multiple ``factor``.
        """
        self.arrheniusLow.changeRate(factor)

################################################################################

cdef class Lindemann(PDepKineticsModel):
    """
    A kinetic model of a phenomenological rate coefficient :math:`k(T, P)`
    using the Lindemann formulation. The attributes are:

    =================== ========================================================
    Attribute           Description
    =================== ========================================================
    `arrheniusHigh`     The Arrhenius kinetics at the high-pressure limit
    `arrheniusLow`      The Arrhenius kinetics at the low-pressure limit
    `Tmin`              The minimum temperature at which the model is valid, or zero if unknown or undefined
    `Tmax`              The maximum temperature at which the model is valid, or zero if unknown or undefined
    `Pmin`              The minimum pressure at which the model is valid, or zero if unknown or undefined
    `Pmax`              The maximum pressure at which the model is valid, or zero if unknown or undefined
    `efficiencies`      A dict associating chemical species with associated efficiencies
    `comment`           Information about the model (e.g. its source)
    =================== ========================================================
    
    """

    def __init__(self, arrheniusHigh=None, arrheniusLow=None, Tmin=None, Tmax=None, Pmin=None, Pmax=None, efficiencies=None, comment=''):
        PDepKineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, efficiencies=efficiencies, comment=comment)
        self.arrheniusHigh = arrheniusHigh
        self.arrheniusLow = arrheniusLow
        
    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        Lindemann object.
        """
        string = 'Lindemann(arrheniusHigh={0!r}, arrheniusLow={1!r}'.format(self.arrheniusHigh, self.arrheniusLow)
        if self.efficiencies is not None: string += ', efficiencies={0!r}'.format(self.efficiencies)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.Pmin is not None: string += ', Pmin={0!r}'.format(self.Pmin)
        if self.Pmax is not None: string += ', Pmax={0!r}'.format(self.Pmax)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling a Lindemann object.
        """
        return (Lindemann, (self.arrheniusHigh, self.arrheniusLow, self.Tmin, self.Tmax, self.Pmin, self.Pmax, self.efficiencies, self.comment))

    cpdef double getRateCoefficient(self, double T, double P=0.0) except -1:
        """
        Return the value of the rate coefficient :math:`k(T)` in units of m^3,
        mol, and s at the specified temperature `T` in K and pressure `P` in
        Pa. If you wish to consider collision efficiencies, then you should
        first use :meth:`getEffectivePressure()` to compute the effective
        pressure, and pass that value as the pressure to this method.
        """
        cdef double C, k0, kinf, Pr
                
        C = P / constants.R / T     # bath gas concentration in mol/m^3
        k0 = self.arrheniusLow.getRateCoefficient(T)
        kinf = self.arrheniusHigh.getRateCoefficient(T)
        Pr = k0 * C / kinf
        
        return kinf * (Pr / (1 + Pr))

    cpdef bint isIdenticalTo(self, KineticsModel otherKinetics) except -2:
        """
        Checks to see if kinetics matches that of other kinetics and returns ``True``
        if coeffs, kunits, Tmin,
        """
        if not isinstance(otherKinetics, Troe):
            return False
        if not KineticsModel.isIdenticalTo(self,otherKinetics):
            return False
        if not self.arrheniusLow.isIdenticalTo(otherKinetics.arrheniusLow):
            return False
        if not self.arrheniusHigh.isIdenticalTo(otherKinetics.arrheniusHigh):
            return False
        
        return True
    
    cpdef changeRate(self, double factor):
        """
        Changes kinetics rate by a multiple ``factor``.
        """
        self.arrheniusLow.changeRate(factor)
        self.arrheniusHigh.changeRate(factor)

################################################################################

cdef class Troe(PDepKineticsModel):
    """
    A kinetic model of a phenomenological rate coefficient :math:`k(T, P)`
    using the Troe formulation. The attributes are:

    =================== ========================================================
    Attribute           Description
    =================== ========================================================
    `arrheniusHigh`     The Arrhenius kinetics at the high-pressure limit
    `arrheniusLow`      The Arrhenius kinetics at the low-pressure limit
    `alpha`             The :math:`\\alpha` parameter
    `T1`                The :math:`T_1` parameter
    `T2`                The :math:`T_2` parameter
    `T3`                The :math:`T_3` parameter
    `Tmin`              The minimum temperature at which the model is valid, or zero if unknown or undefined
    `Tmax`              The maximum temperature at which the model is valid, or zero if unknown or undefined
    `Pmin`              The minimum pressure at which the model is valid, or zero if unknown or undefined
    `Pmax`              The maximum pressure at which the model is valid, or zero if unknown or undefined
    `efficiencies`      A dict associating chemical species with associated efficiencies
    `comment`           Information about the model (e.g. its source)
    =================== ========================================================
    
    """

    def __init__(self, arrheniusHigh=None, arrheniusLow=None, alpha=0.0, T3=None, T1=None, T2=None, Tmin=None, Tmax=None, Pmin=None, Pmax=None, efficiencies=None, comment=''):
        PDepKineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, efficiencies=efficiencies, comment=comment)
        self.arrheniusHigh = arrheniusHigh
        self.arrheniusLow = arrheniusLow
        self.alpha = alpha
        self.T3 = T3
        self.T1 = T1
        self.T2 = T2
        
    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        Troe object.
        """
        string = 'Troe(arrheniusHigh={0!r}, arrheniusLow={1!r}'.format(self.arrheniusHigh, self.arrheniusLow)
        if self.alpha != 0.0: string += ', alpha={0:g}'.format(self.alpha)
        if self.T3 is not None: string += ', T3={0!r}'.format(self.T3)
        if self.T1 is not None: string += ', T1={0!r}'.format(self.T1)
        if self.T2 is not None: string += ', T2={0!r}'.format(self.T2)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.Pmin is not None: string += ', Pmin={0!r}'.format(self.Pmin)
        if self.Pmax is not None: string += ', Pmax={0!r}'.format(self.Pmax)
        if self.efficiencies is not None: string += ', efficiencies={0!r}'.format(self.efficiencies)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling a Troe object.
        """
        return (Troe, (self.arrheniusHigh, self.arrheniusLow, self.alpha, self.T3, self.T1, self.T2, self.Tmin, self.Tmax, self.Pmin, self.Pmax, self.efficiencies, self.comment))

    property T1:
        """The Troe :math:`T_1` parameter."""
        def __get__(self):
            return self._T1
        def __set__(self, value):
            self._T1 = quantity.Temperature(value)

    property T2:
        """The Troe :math:`T_2` parameter."""
        def __get__(self):
            return self._T2
        def __set__(self, value):
            self._T2 = quantity.Temperature(value)

    property T3:
        """The Troe :math:`T_3` parameter."""
        def __get__(self):
            return self._T3
        def __set__(self, value):
            self._T3 = quantity.Temperature(value)

    cpdef double getRateCoefficient(self, double T, double P=0.0) except -1:
        """
        Return the value of the rate coefficient :math:`k(T)` in units of m^3,
        mol, and s at the specified temperature `T` in K and pressure `P` in
        Pa. If you wish to consider collision efficiencies, then you should
        first use :meth:`getEffectivePressure()` to compute the effective
        pressure, and pass that value as the pressure to this method.
        """
        cdef double C, k0, kinf, Pr
        cdef double d, n, c, Fcent, F
        cdef double alpha, T1, T2, T3
        
        C = P / constants.R / T     # bath gas concentration in mol/m^3
        k0 = self.arrheniusLow.getRateCoefficient(T)
        kinf = self.arrheniusHigh.getRateCoefficient(T)
        Pr = k0 * C / kinf
        
        alpha = self.alpha
        T1 = self._T1.value_si if self._T1 is not None else 0.0
        T2 = self._T2.value_si if self._T2 is not None else 0.0
        T3 = self._T3.value_si if self._T3 is not None else 0.0
        
        if T1 == 0 and T3 == 0:
            F = 1.0
        else:
            Fcent = (1 - alpha) * exp(-T / T3) + alpha * exp(-T / T1)
            if T2 != 0.0: Fcent += exp(-T2 / T)
            d = 0.14
            n = 0.75 - 1.27 * log10(Fcent)
            c = -0.4 - 0.67 * log10(Fcent)
            F = 10.0**(log10(Fcent)/(1 + ((log10(Pr) + c)/(n - d * (log10(Pr))))**2))

        return kinf * (Pr / (1 + Pr)) * F

    cpdef bint isIdenticalTo(self, KineticsModel otherKinetics) except -2:
        """
        Checks to see if kinetics matches that of other kinetics and returns ``True``
        if coeffs, kunits, Tmin,
        """
        if not isinstance(otherKinetics, Troe):
            return False
        if not KineticsModel.isIdenticalTo(self,otherKinetics):
            return False
        if not self.arrheniusLow.isIdenticalTo(otherKinetics.arrheniusLow):
            return False
        if not self.arrheniusHigh.isIdenticalTo(otherKinetics.arrheniusHigh):
            return False
        if not self.alpha == otherKinetics.alpha:
            return False
        if not self.T1.equals(otherKinetics.T1):
            return False
        if not self.T3.equals(otherKinetics.T3):
            return False
        if self.T2 is not None and not self.T2.equals(otherKinetics.T2):
            return False
        
        return True
    
    cpdef changeRate(self, double factor):
        """
        Changes kinetics rate by a multiple ``factor``.
        """
        self.arrheniusLow.changeRate(factor)
        self.arrheniusHigh.changeRate(factor)
