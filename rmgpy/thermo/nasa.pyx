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

from libc.math cimport log

cimport rmgpy.constants as constants
import rmgpy.quantity as quantity

################################################################################

cdef class NASAPolynomial(HeatCapacityModel):
    """
    A heat capacity model based on the NASA polynomial. Both the 
    seven-coefficient and nine-coefficient variations are supported.
    The attributes are:
    
    =============== ============================================================
    Attribute       Description
    =============== ============================================================
    `coeffs`        The seven or nine NASA polynomial coefficients
    `Tmin`          The minimum temperature in K at which the model is valid, or zero if unknown or undefined
    `Tmax`          The maximum temperature in K at which the model is valid, or zero if unknown or undefined
    `E0`            The energy at zero Kelvin (including zero point energy)
    `comment`       Information about the model (e.g. its source)
    =============== ============================================================

    """
    
    def __init__(self, coeffs=None, Tmin=None, Tmax=None, E0=None, comment=''):
        HeatCapacityModel.__init__(self, Tmin=Tmin, Tmax=Tmax, E0=E0, comment=comment)
        self.coeffs = coeffs
        
    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        object.
        """
        string = 'NASAPolynomial('
        if self.cm2 == 0 and self.cm1 == 0:
            string += 'coeffs=[{0:g},{1:g},{2:g},{3:g},{4:g},{5:g},{6:g}]'.format(self.c0, self.c1, self.c2, self.c3, self.c4, self.c5, self.c6)
        else:
            string += 'coeffs=[{0:g},{1:g},{2:g},{3:g},{4:g},{5:g},{6:g},{7:g},{8:g}]'.format(self.cm2, self.cm1, self.c0, self.c1, self.c2, self.c3, self.c4, self.c5, self.c6)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.E0 is not None: string += ', E0={0!r}'.format(self.E0)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (NASAPolynomial, ([self.cm2, self.cm1, self.c0, self.c1, self.c2, self.c3, self.c4, self.c5, self.c6], self.Tmin, self.Tmax, self.E0, self.comment))

    property coeffs:
        """The set of seven or nine NASA polynomial coefficients."""
        def __get__(self):
            if self.cm2 == 0 and self.cm1 == 0:
                return numpy.array([self.c0, self.c1, self.c2, self.c3, self.c4, self.c5, self.c6])
            else:
                return numpy.array([self.cm2, self.cm1, self.c0, self.c1, self.c2, self.c3, self.c4, self.c5, self.c6])
        def __set__(self, value):
            if value is None:
                value = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            if len(value) == 7:
                self.c0, self.c1, self.c2, self.c3, self.c4, self.c5, self.c6 = value
                self.cm2 = 0; self.cm1 = 0
            elif len(value) == 9:
                self.cm2, self.cm1, self.c0, self.c1, self.c2, self.c3, self.c4, self.c5, self.c6 = value
            else:
                raise ValueError('Invalid number of NASA polynomial coefficients; expected 7 or 9, got {0:d}.'.format(len(value)))
    
    cpdef double getHeatCapacity(self, double T) except -1000000000:
        """
        Return the constant-pressure heat capacity in J/mol*K at the specified 
        temperature `T` in K.
        """
        return ((self.cm2 / T + self.cm1) / T + self.c0 + T*(self.c1 + T*(self.c2 + T*(self.c3 + self.c4*T)))) * constants.R
    
    cpdef double getEnthalpy(self, double T) except 1000000000:
        """
        Return the enthalpy in J/mol at the specified temperature `T` in K.
        """
        cdef double T2 = T * T
        cdef double T4 = T2 * T2
        return ((-self.cm2 / T + self.cm1 * log(T)) / T + self.c0 + self.c1*T/2. + self.c2*T2/3. + self.c3*T2*T/4. + self.c4*T4/5. + self.c5/T) * constants.R * T
    
    cpdef double getEntropy(self, double T) except -1000000000:
        """
        Return the entropy in J/mol*K at the specified temperature `T` in K.
        """
        cdef double T2 = T * T
        cdef double T4 = T2 * T2
        return ((-self.cm2 / T / 2. - self.cm1) / T + self.c0*log(T) + self.c1*T + self.c2*T2/2. + self.c3*T2*T/3. + self.c4*T4/4. + self.c6) * constants.R
    
    cpdef double getFreeEnergy(self, double T) except 1000000000:
        """
        Return the Gibbs free energy in J/mol at the specified temperature `T`
        in K.
        """
        return self.getEnthalpy(T) - T * self.getEntropy(T)
    
    cpdef changeBaseEnthalpy(self, double deltaH):
        """
        Add deltaH in J/mol to the base enthalpy of formation H298.
        """
        self.c5 += deltaH / constants.R

    cdef double integral2_T0(self, double T):
        """
        Return the value of the dimensionless integral
        
        .. math:: \\int \\left[ \\frac{C_\\mathrm{p}^\\mathrm{NASA}(T')}{R} \\right]^2 \\ dT'
        
        evaluated at the given temperature `T` in K.
        """
        cdef double c0, c1, c2, c3, c4
        cdef double T2, T4, T8, result
        c0, c1, c2, c3, c4 = self.c0, self.c1, self.c2, self.c3, self.c4
        T2 = T * T
        T4 = T2 * T2
        T8 = T4 * T4
        result = (
            c0*c0*T + c0*c1*T2 + 2./3.*c0*c2*T2*T + 0.5*c0*c3*T4 + 0.4*c0*c4*T4*T +
            c1*c1*T2*T/3. + 0.5*c1*c2*T4 + 0.4*c1*c3*T4*T + c1*c4*T4*T2/3. +
            0.2*c2*c2*T4*T + c2*c3*T4*T2/3. + 2./7.*c2*c4*T4*T2*T +
            c3*c3*T4*T2*T/7. + 0.25*c3*c4*T8 +
            c4*c4*T8*T/9.
        )
        return result
    
    cdef double integral2_TM1(self, double T):
        """
        Return the value of the dimensionless integral
        
        .. math:: \\int \\left[ \\frac{C_\\mathrm{p}^\\mathrm{NASA}(T')}{R} \\right]^2 (T')^{-1} \\ dT'
        
        evaluated at the given temperature `T` in K.
        """
        cdef double c0, c1, c2, c3, c4
        cdef double T2, T4, logT, result
        c0, c1, c2, c3, c4 = self.c0, self.c1, self.c2, self.c3, self.c4
        T2 = T * T
        T4 = T2 * T2
        logT = log(T)
        result = (
            c0*c0*logT + 2*c0*c1*T + c0*c2*T2 + 2./3.*c0*c3*T2*T + 0.5*c0*c4*T4 +
            0.5*c1*c1*T2 + 2./3.*c1*c2*T2*T + 0.5*c1*c3*T4 + 0.4*c1*c4*T4*T +
            0.25*c2*c2*T4 + 0.4*c2*c3*T4*T + c2*c4*T4*T2/3. +
            c3*c3*T4*T2/6. + 2./7.*c3*c4*T4*T2*T +
            c4*c4*T4*T4/8.
        )
        return result

################################################################################

cdef class NASA(HeatCapacityModel):
    """
    A heat capacity model based on a set of one, two, or three 
    :class:`NASAPolynomial` objects. The attributes are:
    
    =============== ============================================================
    Attribute       Description
    =============== ============================================================
    `polynomials`   The list of NASA polynomials to use in this model
    `Tmin`          The minimum temperature in K at which the model is valid, or zero if unknown or undefined
    `Tmax`          The maximum temperature in K at which the model is valid, or zero if unknown or undefined
    `E0`            The energy at zero Kelvin (including zero point energy)
    `comment`       Information about the model (e.g. its source)
    =============== ============================================================

    """
    
    def __init__(self, polynomials=None, Tmin=None, Tmax=None, E0=None, comment=''):
        HeatCapacityModel.__init__(self, Tmin=Tmin, Tmax=Tmax, E0=E0, comment=comment)
        self.polynomials = polynomials
    
    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        object.
        """
        polys = self.polynomials
        string = 'NASA(polynomials={0!r}'.format(polys if polys else None)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.E0 is not None: string += ', E0={0!r}'.format(self.E0)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (NASA, (self.polynomials, self.Tmin, self.Tmax, self.E0, self.comment))

    property polynomials:
        """The set of one, two, or three NASA polynomials."""
        def __get__(self):
            polys = []
            if self.poly1 is not None: polys.append(self.poly1)
            if self.poly2 is not None: polys.append(self.poly2)
            if self.poly3 is not None: polys.append(self.poly3)
            return polys
        def __set__(self, value):
            self.poly1 = None
            self.poly2 = None
            self.poly3 = None
            if value is not None:
                if len(value) == 1:
                    self.poly1 = value[0]
                elif len(value) == 2:
                    self.poly1 = value[0]
                    self.poly2 = value[1]
                elif len(value) == 3:
                    self.poly1 = value[0]
                    self.poly2 = value[1]
                    self.poly3 = value[2]
                elif len(value) > 3:
                    raise ValueError('Only one, two, or three NASA polynomials can be stored in a single NASA object.')

    cpdef NASAPolynomial selectPolynomial(self, double T):
        if self.poly1 is not None and self.poly1.isTemperatureValid(T):
            return self.poly1
        elif self.poly2 is not None and self.poly2.isTemperatureValid(T):
            return self.poly2
        elif self.poly3 is not None and self.poly3.isTemperatureValid(T):
            return self.poly3
        else:
            raise ValueError('No valid NASA polynomial at temperature {0:g} K.'.format(T))
    
    cpdef double getHeatCapacity(self, double T) except -1000000000:
        """
        Return the dimensionless constant-pressure heat capacity
        :math:`C_\\mathrm{p}(T)/R` at the specified temperature `T` in K.
        """
        return self.selectPolynomial(T).getHeatCapacity(T)
    
    cpdef double getEnthalpy(self, double T) except 1000000000:
        """
        Return the dimensionless enthalpy :math:`H(T)/RT` at the specified
        temperature `T` in K.
        """
        return self.selectPolynomial(T).getEnthalpy(T)
    
    cpdef double getEntropy(self, double T) except -1000000000:
        """
        Return the dimensionless entropy :math:`S(T)/R` at the specified
        temperature `T` in K.
        """
        return self.selectPolynomial(T).getEntropy(T)
    
    cpdef double getFreeEnergy(self, double T) except 1000000000:
        """
        Return the dimensionless Gibbs free energy :math:`G(T)/RT` at the
        specified temperature `T` in K.
        """
        return self.selectPolynomial(T).getFreeEnergy(T)

    cpdef ThermoData toThermoData(self, double Cp0=0.0, double CpInf=0.0):
        """
        Convert the Wilhoit model to a :class:`ThermoData` object.
        """
        from rmgpy.thermo.thermodata import ThermoData
        
        Tdata = [300,400,500,600,800,1000,1500]
        Cpdata = [self.getHeatCapacity(T) for T in Tdata]
        
        return ThermoData(
            Tdata = (Tdata,"K"),
            Cpdata = (Cpdata,"J/(mol*K)"),
            H298 = (self.getEnthalpy(298)*0.001,"kJ/mol"),
            S298 = (self.getEntropy(298),"J/(mol*K)"),
            Cp0 = (Cp0,"J/(mol*K)"),
            CpInf = (CpInf,"J/(mol*K)"),
            E0 = self.E0,
        )

    @cython.boundscheck(False)
    cpdef Wilhoit toWilhoit(self, double Cp0, double CpInf):
        """
        Convert a :class:`MultiNASA` object `multiNASA` to a :class:`Wilhoit` 
        object. You must specify the linearity of the molecule `linear`, the number
        of vibrational modes `Nfreq`, and the number of hindered rotor modes
        `Nrotors` so the algorithm can determine the appropriate heat capacity
        limits at zero and infinite temperature.
        """
        cdef double Tmin, Tmax, dT, H298, S298
        cdef numpy.ndarray[numpy.float64_t, ndim=1] Tdata, Cpdata
        cdef int i
        
        from rmgpy.thermo.wilhoit import Wilhoit
        
        Tmin = self.Tmin.value_si
        Tmax = self.Tmax.value_si
        dT = min(50.0, (Tmax - Tmin) / 100.)
        
        Tdata = numpy.arange(Tmin, Tmax, dT)
        Cpdata = numpy.zeros_like(Tdata)
        
        for i in range(Tdata.shape[0]):
            Cpdata[i] = self.getHeatCapacity(Tdata[i])
        H298 = self.getEnthalpy(298)
        S298 = self.getEntropy(298)
        
        return Wilhoit().fitToData(Tdata, Cpdata, Cp0, CpInf, H298, S298)
    
    cpdef NASA changeBaseEnthalpy(self, double deltaH):
        """
        Add deltaH in J/mol to the base enthalpy of formation H298 and return the
        modified NASA object.  
        """
        for poly in self.polynomials:
            poly.changeBaseEnthalpy(deltaH)
        return self
    
    def toCantera(self):
        """
        Return the cantera equivalent NasaPoly2 object from this NASA object.
        """
        
        from cantera import NasaPoly2
        import copy

        cdef numpy.ndarray[numpy.float64_t, ndim=1] coeffs
        
        polys = self.polynomials
        assert len(polys) == 2, "Cantera NasaPoly2 objects only accept 2 polynomials"
        assert len(polys[0].coeffs) == 7 and len(polys[1].coeffs) == 7, "Cantera NasaPoly2 polynomials can only contain 7 coefficients."
        
        # In RMG's NASA object, the first polynoamial is low temperature, and the second is 
        # high temperature
        coeffs = numpy.zeros(15)
        coeffs[0] = polys[0].Tmax.value_si # mid point temperature between two polynomials
        coeffs[1:8] = copy.deepcopy(polys[1].coeffs) # 7 coefficients of the high temperature polynomial
        coeffs[8:15] = copy.deepcopy(polys[0].coeffs) # 7 coefficients of the low temperature polynomial

        # initialize cantera.NasaPoly2(T_low, T_high, P_ref, coeffs)
        return NasaPoly2(polys[0].Tmin.value_si, polys[1].Tmax.value_si, 10000.0, coeffs)