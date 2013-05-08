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

import numpy
import cython
from libc.math cimport sqrt, log

cimport rmgpy.constants as constants
import rmgpy.quantity as quantity

################################################################################

cdef class ThermoData(HeatCapacityModel):
    """
    A heat capacity model based on a set of discrete heat capacity data points.
    The attributes are:
    
    =============== ============================================================
    Attribute       Description
    =============== ============================================================
    `Tdata`         An array of temperatures at which the heat capacity is known
    `Cpdata`        An array of heat capacities at the given temperatures
    `H298`          The standard enthalpy of formation at 298 K
    `S298`          The standard entropy at 298 K
    `Cp0`           The heat capacity at zero temperature
    `CpInf`         The heat capacity at infinite temperature
    `Tmin`          The minimum temperature at which the model is valid, or zero if unknown or undefined
    `Tmax`          The maximum temperature at which the model is valid, or zero if unknown or undefined
    `E0`            The energy at zero Kelvin (including zero point energy)
    `comment`       Information about the model (e.g. its source)
    =============== ============================================================
    
    """

    def __init__(self, Tdata=None, Cpdata=None, H298=None, S298=None, Cp0=None, CpInf=None, Tmin=None, Tmax=None, E0=None, comment=''):
        HeatCapacityModel.__init__(self, Tmin=Tmin, Tmax=Tmax, E0=E0, comment=comment)
        self.H298 = H298
        self.S298 = S298
        self.Tdata = Tdata
        self.Cpdata = Cpdata
        self.Cp0 = Cp0
        self.CpInf = CpInf
    
    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        ThermoData object.
        """
        string = 'ThermoData(Tdata={0!r}, Cpdata={1!r}, H298={2!r}, S298={3!r}'.format(self.Tdata, self.Cpdata, self.H298, self.S298)
        if self.Cp0 is not None: string += ', Cp0={0!r}'.format(self.Cp0)
        if self.CpInf is not None: string += ', CpInf={0!r}'.format(self.CpInf)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.E0 is not None: string += ', E0={0!r}'.format(self.E0)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling a ThermoData object.
        """
        return (ThermoData, (self.Tdata, self.Cpdata, self.H298, self.S298, self.Cp0, self.CpInf, self.Tmin, self.Tmax, self.E0, self.comment))

    property Tdata:
        """An array of temperatures at which the heat capacity is known."""
        def __get__(self):
            return self._Tdata
        def __set__(self, value):
            self._Tdata = quantity.Temperature(value)

    property Cpdata:
        """An array of heat capacities at the given temperatures."""
        def __get__(self):
            return self._Cpdata
        def __set__(self, value):
            self._Cpdata = quantity.HeatCapacity(value)

    property H298:
        """The standard enthalpy of formation at 298 K."""
        def __get__(self):
            return self._H298
        def __set__(self, value):
            self._H298 = quantity.Enthalpy(value)

    property S298:
        """The standard entropy of formation at 298 K."""
        def __get__(self):
            return self._S298
        def __set__(self, value):
            self._S298 = quantity.Entropy(value)

    property Cp0:
        """The heat capacity at zero temperature."""
        def __get__(self):
            return self._Cp0
        def __set__(self, value):
            self._Cp0 = quantity.HeatCapacity(value)

    property CpInf:
        """The heat capacity at infinite temperature."""
        def __get__(self):
            return self._CpInf
        def __set__(self, value):
            self._CpInf = quantity.HeatCapacity(value)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef double getHeatCapacity(self, double T) except -1000000000:
        """
        Return the constant-pressure heat capacity in J/mol*K at the specified
        temperature `T` in K.
        """
        cdef numpy.ndarray[numpy.float64_t,ndim=1] Tdata, Cpdata
        cdef double Cp0, CpInf
        cdef double Tlow, Thigh, Cplow, Cphigh
        cdef double Cp
        cdef int i, N
        
        Tdata = self._Tdata.value_si
        Cpdata = self._Cpdata.value_si
        Cp0 = self._Cp0.value_si if self._Cp0 is not None else 0.0
        CpInf = self._CpInf.value_si if self._CpInf is not None else 0.0
        N = Cpdata.shape[0]
        
        # If T is outside the range of Cp data, make sure we have a value of
        # Cp0 or CpInf we need for reasonable extrapolation
        # Allow a small window of safe extrapolation past each endpoint
        # (ostensibly so we can extrapolate from 300 K to 298 K)
        if T < Tdata[0] - 5 and Cp0 == 0.0:
            raise ValueError('Unable to compute heat capacity at {0:g} K using ThermoData model; please supply a value for Cp0.'.format(T))
        elif T > Tdata[N-1] + 5 and CpInf == 0.0:
            raise ValueError('Unable to compute heat capacity at {0:g} K using ThermoData model; please supply a value for CpInf.'.format(T))

        if T < Tdata[0] and Cp0 != 0.0:
            # Extrapolate towards zero temperature using the slope at the lowest temperature
            # However, if the computed value is less than Cp0, use Cp0 instead
            Tlow = Tdata[0]; Thigh = Tdata[1]
            Cplow = Cpdata[0]; Cphigh = Cpdata[1]
            Cp = Cplow + (T - Tlow) / (Thigh - Tlow) * (Cphigh - Cplow)
            if Cp < Cp0: Cp = Cp0
        elif T > Tdata[N-1] and CpInf != 0.0:
            # Extrapolate towards infinite temperature using the slope at 1500 K
            # However, if the computed value is greater than CpInf, use CpInf instead
            Tlow = Tdata[N-2]; Thigh = Tdata[N-1]
            Cplow = Cpdata[N-2]; Cphigh = Cpdata[N-1]
            Cp = Cplow + (T - Tlow) / (Thigh - Tlow) * (Cphigh - Cplow)
            if Cp > CpInf: Cp = CpInf
        else:
            for i in range(N-1):
                Tlow = Tdata[i]; Thigh = Tdata[i+1]
                Cplow = Cpdata[i]; Cphigh = Cpdata[i+1]
                if Tlow <= T and T <= Thigh:
                    Cp = (Cphigh - Cplow) * ((T - Tlow) / (Thigh - Tlow)) + Cplow
                    break
                
        return Cp
    
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef double getEnthalpy(self, double T) except 1000000000:
        """
        Return the enthalpy in J/mol at the specified temperature `T` in K.
        """
        cdef numpy.ndarray[numpy.float64_t,ndim=1] Tdata, Cpdata
        cdef double Cp0, CpInf
        cdef double Tlow, Thigh, Cplow, Cphigh
        cdef double H, slope, intercept, T0
        cdef int i, N

        Tdata = self._Tdata.value_si
        Cpdata = self._Cpdata.value_si
        Cp0 = self._Cp0.value_si if self._Cp0 is not None else 0.0
        CpInf = self._CpInf.value_si if self._CpInf is not None else 0.0
        N = Cpdata.shape[0]

        H = self._H298.value_si
        
        # If T is outside the range of Cp data, make sure we have a value of
        # Cp0 or CpInf we need for reasonable extrapolation
        # Allow a small window of safe extrapolation past each endpoint
        # (ostensibly so we can extrapolate from 300 K to 298 K)
        if T < Tdata[0] - 5 and Cp0 == 0.0:
            raise ValueError('Unable to compute enthalpy at {0:g} K using ThermoData model; please supply a value for Cp0.'.format(T))
        elif T > Tdata[N-1] + 5 and CpInf == 0.0:
            raise ValueError('Unable to compute enthalpy at {0:g} K using ThermoData model; please supply a value for CpInf.'.format(T))

        # Correct the enthalpy from 298 K to the temperature of the lowest heat capacity point
        assert Tdata[0] >= 298
        Tlow = Tdata[0]; Thigh = Tdata[1]
        Cplow = Cpdata[0]; Cphigh = Cpdata[1]
        slope = (Cphigh - Cplow) / (Thigh - Tlow)
        intercept = (Cplow * Thigh - Cphigh * Tlow) / (Thigh - Tlow)
        H -= 0.5 * slope * (298*298 - Tlow*Tlow) + intercept * (298 - Tlow)
        
        if T < Tdata[0]:
            Tlow = Tdata[0]; Thigh = Tdata[1]
            Cplow = Cpdata[0]; Cphigh = Cpdata[1]
            slope = (Cphigh - Cplow) / (Thigh - Tlow)
            intercept = (Cplow * Thigh - Cphigh * Tlow) / (Thigh - Tlow)
            T0 = (Cp0 - Tlow) / slope + Tlow
            if T > T0 or slope <= 0 or T0 > Tlow:
                H += 0.5 * slope * (T*T - Tlow*Tlow) + intercept * (T - Tlow)
            else:
                H += 0.5 * slope * (T0*T0 - Tlow*Tlow) + intercept * (T0 - Tlow) + Cp0 * (T0 - T)
        
        for i in range(N-1):
            Tlow = Tdata[i]; Thigh = Tdata[i+1]
            Cplow = Cpdata[i]; Cphigh = Cpdata[i+1]
            if T > Tlow:
                slope = (Cphigh - Cplow) / (Thigh - Tlow)
                intercept = (Cplow * Thigh - Cphigh * Tlow) / (Thigh - Tlow)
                if T <= Thigh: 
                    H += 0.5 * slope * (T*T - Tlow*Tlow) + intercept * (T - Tlow)
                    break
                else:
                    H += 0.5 * slope * (Thigh*Thigh - Tlow*Tlow) + intercept * (Thigh - Tlow)

        if T > Tdata[N-1]:
            Tlow = Tdata[N-2]; Thigh = Tdata[N-1]
            Cplow = Cpdata[N-2]; Cphigh = Cpdata[N-1]
            slope = (Cphigh - Cplow) / (Thigh - Tlow)
            intercept = (Cplow * Thigh - Cphigh * Tlow) / (Thigh - Tlow)
            T0 = (CpInf - Cphigh) / slope + Thigh
            if T <= T0 or slope <= 0 or T0 < Thigh:
                H += 0.5 * slope * (T*T - Thigh*Thigh) + intercept * (T - Thigh)
            else:
                H += 0.5 * slope * (T0*T0 - Thigh*Thigh) + intercept * (T0 - Thigh) + CpInf * (T - T0)
        
        return H
        
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef double getEntropy(self, double T) except -1000000000:
        """
        Return the entropy in J/mol*K at the specified temperature `T` in K.
        """
        cdef numpy.ndarray[numpy.float64_t,ndim=1] Tdata, Cpdata
        cdef double Cp0, CpInf
        cdef double Tlow, Thigh, Cplow, Cphigh
        cdef double S, slope, intercept, T0
        cdef int i, N

        Tdata = self._Tdata.value_si
        Cpdata = self._Cpdata.value_si
        Cp0 = self._Cp0.value_si if self._Cp0 is not None else 0.0
        CpInf = self._CpInf.value_si if self._CpInf is not None else 0.0
        N = Cpdata.shape[0]

        S = self._S298.value_si
         
        # Correct the entropy from 298 K to the temperature of the lowest heat capacity point
        assert Tdata[0] > 298
        Tlow = Tdata[0]; Thigh = Tdata[1]
        Cplow = Cpdata[0]; Cphigh = Cpdata[1]
        slope = (Cphigh - Cplow) / (Thigh - Tlow)
        intercept = (Cplow * Thigh - Cphigh * Tlow) / (Thigh - Tlow)
        S -= slope * (298 - Tlow) + intercept * log(298 / Tlow)
        
        # If T is outside the range of Cp data, make sure we have a value of
        # Cp0 or CpInf we need for reasonable extrapolation
        # Allow a small window of safe extrapolation past each endpoint
        # (ostensibly so we can extrapolate from 300 K to 298 K)
        if T < Tdata[0] - 5 and Cp0 == 0.0:
            raise ValueError('Unable to compute entropy at {0:g} K using ThermoData model; please supply a value for Cp0.'.format(T))
        elif T > Tdata[N-1] + 5 and CpInf == 0.0:
            raise ValueError('Unable to compute entropy at {0:g} K using ThermoData model; please supply a value for CpInf.'.format(T))

        if T < Tdata[0]:
            Tlow = Tdata[0]; Thigh = Tdata[1]
            Cplow = Cpdata[0]; Cphigh = Cpdata[1]
            slope = (Cphigh - Cplow) / (Thigh - Tlow)
            intercept = (Cplow * Thigh - Cphigh * Tlow) / (Thigh - Tlow)
            T0 = (Cp0 - Tlow) / slope + Tlow
            if T > T0 or slope <= 0 or T0 >= Tlow:
                S += slope * (T - Tlow) + intercept * log(T / Tlow)
            else:
                S += slope * (T0 - Tlow) + intercept * log(T0 / Tlow) + Cp0 * log(T0 / T)
        
        for i in range(N-1):
            Tlow = Tdata[i]; Thigh = Tdata[i+1]
            Cplow = Cpdata[i]; Cphigh = Cpdata[i+1]
            if T > Tlow:
                slope = (Cphigh - Cplow) / (Thigh - Tlow)
                intercept = (Cplow * Thigh - Cphigh * Tlow) / (Thigh - Tlow)
                if T <= Thigh:
                    S += slope * (T - Tlow) + intercept * log(T / Tlow)
                else:
                    S += slope * (Thigh - Tlow) + intercept * log(Thigh / Tlow)

        if T > Tdata[N-1]:
            Tlow = Tdata[N-2]; Thigh = Tdata[N-1]
            Cplow = Cpdata[N-2]; Cphigh = Cpdata[N-1]
            slope = (Cphigh - Cplow) / (Thigh - Tlow)
            intercept = (Cplow * Thigh - Cphigh * Tlow) / (Thigh - Tlow)
            if slope > 0:
                T0 = (CpInf - Cphigh) / slope + Thigh
                if T <= T0:
                    S += slope * (T - Thigh) + intercept * log(T / Thigh)
                else:
                    S += slope * (T0 - Thigh) + intercept * log(T0 / Thigh) + CpInf * log(T / T0)

        return S
    
    cpdef double getFreeEnergy(self, double T) except 1000000000:
        """
        Return the Gibbs free energy in J/mol at the specified temperature
        `T` in K.
        """
        cdef numpy.ndarray[numpy.float64_t,ndim=1] Tdata
        cdef double Cp0, CpInf
        cdef int N

        Tdata = self._Tdata.value_si
        Cp0 = self._Cp0.value_si if self._Cp0 is not None else 0.0
        CpInf = self._CpInf.value_si if self._CpInf is not None else 0.0
        N = Tdata.shape[0]

        # If T is outside the range of Cp data, make sure we have a value of
        # Cp0 or CpInf we need for reasonable extrapolation
        # Allow a small window of safe extrapolation past each endpoint
        # (ostensibly so we can extrapolate from 300 K to 298 K)
        if T < Tdata[0] - 5 and Cp0 == 0.0:
            raise ValueError('Unable to compute Gibbs free energy at {0:g} K using ThermoData model; please supply a value for Cp0.'.format(T))
        elif T > Tdata[N-1] + 5 and CpInf == 0.0:
            raise ValueError('Unable to compute Gibbs free energy at {0:g} K using ThermoData model; please supply a value for CpInf.'.format(T))

        return self.getEnthalpy(T) - T * self.getEntropy(T)

    cpdef Wilhoit toWilhoit(self):
        """
        Convert the Benson model to a Wilhoit model. For the conversion to
        succeed, you must have set the `Cp0` and `CpInf` attributes of the
        Benson model.
        """
        if self.Cp0 == 0.0 or self.CpInf == 0.0:
            raise Exception('Cannot convert Benson model to Wilhoit model; first specify Cp0 and CpInf.')
        from rmgpy.thermo.wilhoit import Wilhoit
        
        Tdata = self._Tdata.value_si
        Cpdata = self._Cpdata.value_si
        H298 = self.getEnthalpy(298)
        S298 = self.getEntropy(298)
        Cp0 = self._Cp0.value_si
        CpInf = self._CpInf.value_si
        
        return Wilhoit().fitToData(Tdata, Cpdata, Cp0, CpInf, H298, S298)

    cpdef NASA toNASA(self, double Tmin, double Tmax, double Tint, bint fixedTint=False, bint weighting=True, int continuity=3):
        """
        Convert the object to a :class:`NASA` object. You must specify the
        minimum and maximum temperatures of the fit `Tmin` and `Tmax` in K, as
        well as the intermediate temperature `Tint` in K to use as the bridge
        between the two fitted polynomials. The remaining parameters can be
        used to modify the fitting algorithm used:
        
        * `fixedTint` - ``False`` to allow `Tint` to vary in order to improve the fit, or ``True`` to keep it fixed
    
        * `weighting` - ``True`` to weight the fit by :math:`T^{-1}` to emphasize good fit at lower temperatures, or ``False`` to not use weighting
    
        * `continuity` - The number of continuity constraints to enforce at `Tint`:
    
            - 0: no constraints on continuity of :math:`C_\\mathrm{p}(T)` at `Tint`
    
            - 1: constrain :math:`C_\\mathrm{p}(T)` to be continous at `Tint`
    
            - 2: constrain :math:`C_\\mathrm{p}(T)` and :math:`\\frac{d C_\\mathrm{p}}{dT}` to be continuous at `Tint`
    
            - 3: constrain :math:`C_\\mathrm{p}(T)`, :math:`\\frac{d C_\\mathrm{p}}{dT}`, and :math:`\\frac{d^2 C_\\mathrm{p}}{dT^2}` to be continuous at `Tint`
    
            - 4: constrain :math:`C_\\mathrm{p}(T)`, :math:`\\frac{d C_\\mathrm{p}}{dT}`, :math:`\\frac{d^2 C_\\mathrm{p}}{dT^2}`, and :math:`\\frac{d^3 C_\\mathrm{p}}{dT^3}` to be continuous at `Tint`
    
            - 5: constrain :math:`C_\\mathrm{p}(T)`, :math:`\\frac{d C_\\mathrm{p}}{dT}`, :math:`\\frac{d^2 C_\\mathrm{p}}{dT^2}`, :math:`\\frac{d^3 C_\\mathrm{p}}{dT^3}`, and :math:`\\frac{d^4 C_\\mathrm{p}}{dT^4}` to be continuous at `Tint`
            
        Note that values of `continuity` of 5 or higher effectively constrain all
        the coefficients to be equal and should be equivalent to fitting only one
        polynomial (rather than two).
    
        Returns the fitted :class:`NASA` object containing the two fitted
        :class:`NASAPolynomial` objects.
        """
        return self.toWilhoit().toNASA(Tmin, Tmax, Tint, fixedTint, weighting, continuity)
