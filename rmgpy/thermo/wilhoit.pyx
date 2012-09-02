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

cdef class Wilhoit(HeatCapacityModel):
    """
    A heat capacity model based on the Wilhoit equation. The attributes are:
    
    =============== ============================================================
    Attribute       Description
    =============== ============================================================
    `Cp0`           The heat capacity at zero temperature
    `CpInf`         The heat capacity at infinite temperature
    `a0`            The zeroth-order Wilhoit polynomial coefficient
    `a1`            The first-order Wilhoit polynomial coefficient
    `a2`            The second-order Wilhoit polynomial coefficient
    `a3`            The third-order Wilhoit polynomial coefficient
    `H0`            The integration constant for enthalpy
    `S0`            The integration constant for entropy
    `B`             The Wilhoit scaled temperature coefficient in K
    `Tmin`          The minimum temperature in K at which the model is valid, or zero if unknown or undefined
    `Tmax`          The maximum temperature in K at which the model is valid, or zero if unknown or undefined
    `comment`       Information about the model (e.g. its source)
    =============== ============================================================
    
    """

    def __init__(self, Cp0=None, CpInf=None, a0=0.0, a1=0.0, a2=0.0, a3=0.0, H0=None, S0=None, B=None, Tmin=None, Tmax=None, comment=''):
        HeatCapacityModel.__init__(self, Tmin=Tmin, Tmax=Tmax, comment=comment)
        self.Cp0 = Cp0
        self.CpInf = CpInf
        self.B = B
        self.a0 = a0
        self.a1 = a1
        self.a2 = a2
        self.a3 = a3
        self.H0 = H0
        self.S0 = S0
    
    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        Wilhoit object.
        """
        string = 'Wilhoit(Cp0={0!r}, CpInf={1!r}, a0={2:g}, a1={3:g}, a2={4:g}, a3={5:g}, H0={6!r}, S0={7!r}, B={8!r}'.format(
            self.Cp0, self.CpInf, self.a0, self.a1, self.a2, self.a3, self.H0, self.S0, self.B)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling a Wilhoit object.
        """
        return (Wilhoit, (self.Cp0, self.CpInf, self.a0, self.a1, self.a2, self.a3, self.H0, self.S0, self.B, self.Tmin, self.Tmax, self.comment))

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

    property B:
        """The Wilhoit scaled temperature coefficient."""
        def __get__(self):
            return self._B
        def __set__(self, value):
            self._B = quantity.Temperature(value)

    property H0:
        """The integration constant for enthalpy."""
        def __get__(self):
            return self._H0
        def __set__(self, value):
            self._H0 = quantity.Enthalpy(value)

    property S0:
        """The integration constant for entropy."""
        def __get__(self):
            return self._S0
        def __set__(self, value):
            self._S0 = quantity.Entropy(value)

    cpdef double getHeatCapacity(self, double T) except -1000000000:
        """
        Return the constant-pressure heat capacity in J/mol*K at the specified
        temperature `T` in K.
        """
        cdef double Cp0, CpInf, B, a0, a1, a2, a3
        cdef double y
        Cp0, CpInf, B, a0, a1, a2, a3 = self._Cp0.value_si, self._CpInf.value_si, self._B.value_si, self.a0, self.a1, self.a2, self.a3
        y = T / (T + B)
        return Cp0 + (CpInf - Cp0) * y * y * (
            1 + (y - 1) * (a0 + y * (a1 + y * (a2 + y * a3))) 
        )
            
    cpdef double getEnthalpy(self, double T) except 1000000000:
        """
        Return the enthalpy in J/mol at the specified temperature `T` in K.
        """
        cdef double Cp0, CpInf, B, a0, a1, a2, a3
        cdef double y
        Cp0, CpInf, B, a0, a1, a2, a3 = self._Cp0.value_si, self._CpInf.value_si, self._B.value_si, self.a0, self.a1, self.a2, self.a3
        y = T / (T + B)
        return self._H0.value_si + Cp0 * T - (CpInf - Cp0) * T * (
            y * y * ((3 * a0 + a1 + a2 + a3) / 6. + 
                (4 * a1 + a2 + a3) * y / 12. + 
                (5 * a2 + a3) * y * y / 20. + 
                a3 * y * y * y / 5.) + 
            (2 + a0 + a1 + a2 + a3) * (y / 2. - 1 + (1.0 / y - 1.) * log(B + T))
        )
    
    cpdef double getEntropy(self, double T) except -1000000000:
        """
        Return the entropy in J/mol*K at the specified temperature `T` in K.
        """
        cdef double Cp0, CpInf, B, a0, a1, a2, a3
        cdef double y, logT, logy
        Cp0, CpInf, B, a0, a1, a2, a3 = self._Cp0.value_si, self._CpInf.value_si, self._B.value_si, self.a0, self.a1, self.a2, self.a3
        y = T / (T + B)
        logT = log(T)
        logy = log(y)
        return self._S0.value_si + CpInf * logT - (CpInf - Cp0) * (
            logy + y * (1 + y * (a0 / 2. + y * (a1 / 3. + y * (a2 / 4. + y * a3 / 5.))))
        )
    
    cpdef double getFreeEnergy(self, double T) except 1000000000:
        """
        Return the Gibbs free energy in J/mol at the specified temperature `T`
        in K.
        """
        return self.getEnthalpy(T) - self.getEntropy(T)
    
    cpdef Wilhoit copy(self):
        """
        Return a copy of the Wilhoit object.
        """
        return Wilhoit(
            self.Cp0, self.CpInf, 
            self.a0, self.a1, self.a2, self.a3, 
            self.H0, self.S0, self.B, 
            Tmin=self.Tmin, Tmax=self.Tmax, comment=self.comment,
        )
    
    def __residual(self, B, Tdata, Cpdata, Cp0, CpInf, H298, S298):
        # The residual corresponding to the fitToData() method
        # Parameters are the same as for that method
        cdef double res = 0.0, diff
        cdef int i
        self.fitToDataForConstantB(Tdata, Cpdata, Cp0, CpInf, H298, S298, B)
        # Objective function is linear least-squares
        for i in range(Cpdata.shape[0]):
            diff = self.getHeatCapacity(Tdata[i]) - Cpdata[i]
            res += diff * diff
        return res
    
    def fitToData(self, 
                  numpy.ndarray[numpy.float64_t, ndim=1] Tdata, 
                  numpy.ndarray[numpy.float64_t, ndim=1] Cpdata, 
                  double Cp0, double CpInf,
                  double H298, double S298, double B0=500.0):
        """
        Fit a Wilhoit model to the data points provided, allowing the 
        characteristic temperature `B` to vary so as to improve the fit. This
        procedure requires an optimization, using the ``fminbound`` function
        in the ``scipy.optimize`` module. The data consists of a set
        of heat capacity points `Cpdata` in J/mol*K at a given set of 
        temperatures `Tdata` in K, along with the enthalpy `H298` in kJ/mol and
        entropy `S298` in J/mol*K at 298 K. The linearity of the molecule, 
        number of vibrational frequencies, and number of internal rotors 
        (`linear`, `Nfreq`, and `Nrotors`, respectively) is used to set the 
        limits at zero and infinite temperature.
        """
        self.B = (B0,"K")
        import scipy.optimize
        scipy.optimize.fminbound(self.__residual, 300.0, 3000.0, args=(Tdata, Cpdata, Cp0, CpInf, H298, S298))
        return self
    
    def fitToDataForConstantB(self, 
                              numpy.ndarray[numpy.float64_t, ndim=1] Tdata, 
                              numpy.ndarray[numpy.float64_t, ndim=1] Cpdata, 
                              double Cp0, double CpInf,
                              double H298, double S298, double B):
        """
        Fit a Wilhoit model to the data points provided using a specified value
        of the characteristic temperature `B`. The data consists of a set
        of dimensionless heat capacity points `Cpdata` at a given set of 
        temperatures `Tdata` in K, along with the dimensionless heat capacity
        at zero and infinite temperature, the dimensionless enthalpy `H298` at 
        298 K, and the dimensionless entropy `S298` at 298 K. 
        """
        cdef numpy.ndarray[numpy.float64_t, ndim=1] b, x
        cdef numpy.ndarray[numpy.float64_t, ndim=2] A
        cdef double y
        cdef int i, j
        
        self.Cp0 = (Cp0,"J/(mol*K)")
        self.CpInf = (CpInf,"J/(mol*K)")
            
        if Cp0 == CpInf:
            # The heat capacity is constant at all temperatures 
            # (i.e. probably a monatomic species)
            self.B = (B,"K")
            self.a0 = 0.0
            self.a1 = 0.0
            self.a2 = 0.0
            self.a3 = 0.0
    
        else:
            # Polyatomic species
                
            # What remains is to fit the polynomial coefficients (a0, a1, a2, a3)
            # This can be done directly - no iteration required
            A = numpy.empty((Cpdata.shape[0],4), numpy.float64)
            b = numpy.empty(Cpdata.shape[0], numpy.float64)
            for i in range(Cpdata.shape[0]):
                y = Tdata[i] / (Tdata[i] + B)
                for j in range(4):
                    A[i,j] = (y*y*y - y*y) * y**j
                b[i] = ((Cpdata[i] - Cp0) / (CpInf - Cp0) - y*y)
            x, residues, rank, s = numpy.linalg.lstsq(A, b)
            
            self.B = (float(B),"K")
            self.a0 = float(x[0])
            self.a1 = float(x[1])
            self.a2 = float(x[2])
            self.a3 = float(x[3])

        self.H0 = (0.0,"kJ/mol")
        self.S0 = (0.0,"J/(mol*K)")
        self._H0.value_si = H298 - self.getEnthalpy(298)
        self._S0.value_si = S298 - self.getEntropy(298)

        return self
