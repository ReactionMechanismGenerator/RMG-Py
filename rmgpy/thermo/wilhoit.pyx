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
import scipy.linalg

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
    `H0`            The integration constant for enthalpy (not H at T=0)
    `S0`            The integration constant for entropy (not S at T=0)
    `E0`            The energy at zero Kelvin (including zero point energy)
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
        """The integration constant for enthalpy.
        
        NB. this is not equal to the enthlapy at 0 Kelvin, which you can access via E0"""
        def __get__(self):
            return self._H0
        def __set__(self, value):
            self._H0 = quantity.Enthalpy(value)

    property E0:
        """The ground state energy (J/mol) at zero Kelvin, including zero point energy.
        
        For the Wilhoit class, this is calculated as the Enthalpy at 0.001 Kelvin."""
        def __get__(self):
            cdef double E0
            E0 = self.getEnthalpy(0.001) # in J/mol
            return quantity.Enthalpy(E0 * 0.001, "kJ/mol")
        def __set__(self, value):
            assert value is None, "You should not be setting E0 on a Wilhoit object - it is determined from the Enthalpy at 0.001 Kelvin."

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
        return self.getEnthalpy(T) - T * self.getEntropy(T)
    
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
    
    @cython.boundscheck(False)
    @cython.wraparound(False)
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
    
    @cython.boundscheck(False)
    @cython.wraparound(False)
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
    
    @cython.boundscheck(False)
    @cython.wraparound(False)
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

    cdef double integral_T0(self, double T):
        """
        Return the value of the dimensionless integral
        
        .. math:: \\int \\frac{C_\\mathrm{p}^\\mathrm{Wilhoit}(T')}{R} \\ dT'
        
        evaluated at the given temperature `T` in kK. The implementation
        differs from that given in the Yelvington thesis for enthalpy by a 
        parameter-dependent (but temperature-independent) constant.
        """
        cdef double Cp0, CpInf, B, a0, a1, a2, a3
        cdef double y, y2, logBplusT, result
        Cp0, CpInf, B, a0, a1, a2, a3 = self._Cp0.value_si, self._CpInf.value_si, self._B.value_si, self.a0, self.a1, self.a2, self.a3
        y = T / (T + B)
        y2 = y * y
        logBplusT = log(B + T)
        result = Cp0*T - (CpInf-Cp0)*T*(y2*((3*a0 + a1 + a2 + a3)/6. + (4*a1 + a2 + a3)*y/12. + (5*a2 + a3)*y2/20. + a3*y2*y/5.) + (2 + a0 + a1 + a2 + a3)*( y/2. - 1 + (1/y-1)*logBplusT))
        return result
    
    cdef double integral_TM1(self, double T):
        """
        Return the value of the dimensionless integral
        
        .. math:: \\int \\frac{C_\\mathrm{p}^\\mathrm{Wilhoit}(T')}{R} (T')^{-1} \\ dT'
        
        evaluated at the given temperature `T` in kK. The implementation
        differs from that given in the Yelvington thesis for entropy by a 
        parameter-dependent (but temperature-independent) constant.
        """
        cdef double Cp0, CpInf, B, a0, a1, a2, a3
        cdef double y, logy, logT, result
        Cp0, CpInf, B, a0, a1, a2, a3 = self._Cp0.value_si, self._CpInf.value_si, self._B.value_si, self.a0, self.a1, self.a2, self.a3
        y = T / (T + B)
        logy = log(y)
        logT = log(T)
        result = CpInf*logT-(CpInf-Cp0)*(logy+y*(1+y*(a0/2+y*(a1/3 + y*(a2/4 + y*a3/5)))))
        return result
    
    cdef double integral_T1(self, double T):
        """
        Return the value of the dimensionless integral
        
        .. math:: \\int \\frac{C_\\mathrm{p}^\\mathrm{Wilhoit}(T')}{R} T' \\ dT'
        
        evaluated at the given temperature `T` in kK.
        """
        cdef double Cp0, CpInf, B, a0, a1, a2, a3
        cdef double logBplusT, result
        Cp0, CpInf, B, a0, a1, a2, a3 = self._Cp0.value_si, self._CpInf.value_si, self._B.value_si, self.a0, self.a1, self.a2, self.a3
        logBplusT = log(B + T)
        result = ( (2 + a0 + a1 + a2 + a3)*B*(Cp0 - CpInf)*T + (CpInf*T**2)/2. + (a3*B**7*(-Cp0 + CpInf))/(5.*(B + T)**5) + ((a2 + 6*a3)*B**6*(Cp0 - CpInf))/(4.*(B + T)**4) -
            ((a1 + 5*(a2 + 3*a3))*B**5*(Cp0 - CpInf))/(3.*(B + T)**3) + ((a0 + 4*a1 + 10*(a2 + 2*a3))*B**4*(Cp0 - CpInf))/(2.*(B + T)**2) -
            ((1 + 3*a0 + 6*a1 + 10*a2 + 15*a3)*B**3*(Cp0 - CpInf))/(B + T) - (3 + 3*a0 + 4*a1 + 5*a2 + 6*a3)*B**2*(Cp0 - CpInf)*logBplusT)
        return result
    
    cdef double integral_T2(self, double T):
        """
        Return the value of the dimensionless integral
        
        .. math:: \\int \\frac{C_\\mathrm{p}^\\mathrm{Wilhoit}(T')}{R} (T')^2 \\ dT'
        
        evaluated at the given temperature `T` in kK.
        """
        cdef double Cp0, CpInf, B, a0, a1, a2, a3
        cdef double logBplusT, result
        Cp0, CpInf, B, a0, a1, a2, a3 = self._Cp0.value_si, self._CpInf.value_si, self._B.value_si, self.a0, self.a1, self.a2, self.a3
        logBplusT = log(B + T)
        result = ( -((3 + 3*a0 + 4*a1 + 5*a2 + 6*a3)*B**2*(Cp0 - CpInf)*T) + ((2 + a0 + a1 + a2 + a3)*B*(Cp0 - CpInf)*T**2)/2. + (CpInf*T**3)/3. + (a3*B**8*(Cp0 - CpInf))/(5.*(B + T)**5) -
            ((a2 + 7*a3)*B**7*(Cp0 - CpInf))/(4.*(B + T)**4) + ((a1 + 6*a2 + 21*a3)*B**6*(Cp0 - CpInf))/(3.*(B + T)**3) - ((a0 + 5*(a1 + 3*a2 + 7*a3))*B**5*(Cp0 - CpInf))/(2.*(B + T)**2) +
            ((1 + 4*a0 + 10*a1 + 20*a2 + 35*a3)*B**4*(Cp0 - CpInf))/(B + T) + (4 + 6*a0 + 10*a1 + 15*a2 + 21*a3)*B**3*(Cp0 - CpInf)*logBplusT)
        return result
    
    cdef double integral_T3(self, double T):
        """
        Return the value of the dimensionless integral
        
        .. math:: \\int \\frac{C_\\mathrm{p}^\\mathrm{Wilhoit}(T')}{R} (T')^3 \\ dT'
        
        evaluated at the given temperature `T` in kK.
        """
        cdef double Cp0, CpInf, B, a0, a1, a2, a3
        cdef double logBplusT, result
        Cp0, CpInf, B, a0, a1, a2, a3 = self._Cp0.value_si, self._CpInf.value_si, self._B.value_si, self.a0, self.a1, self.a2, self.a3
        logBplusT = log(B + T)
        result = ( (4 + 6*a0 + 10*a1 + 15*a2 + 21*a3)*B**3*(Cp0 - CpInf)*T + ((3 + 3*a0 + 4*a1 + 5*a2 + 6*a3)*B**2*(-Cp0 + CpInf)*T**2)/2. + ((2 + a0 + a1 + a2 + a3)*B*(Cp0 - CpInf)*T**3)/3. +
            (CpInf*T**4)/4. + (a3*B**9*(-Cp0 + CpInf))/(5.*(B + T)**5) + ((a2 + 8*a3)*B**8*(Cp0 - CpInf))/(4.*(B + T)**4) - ((a1 + 7*(a2 + 4*a3))*B**7*(Cp0 - CpInf))/(3.*(B + T)**3) +
            ((a0 + 6*a1 + 21*a2 + 56*a3)*B**6*(Cp0 - CpInf))/(2.*(B + T)**2) - ((1 + 5*a0 + 15*a1 + 35*a2 + 70*a3)*B**5*(Cp0 - CpInf))/(B + T) -
            (5 + 10*a0 + 20*a1 + 35*a2 + 56*a3)*B**4*(Cp0 - CpInf)*logBplusT)
        return result
    
    cdef double integral_T4(self, double T):
        """
        Return the value of the dimensionless integral
        
        .. math:: \\int \\frac{C_\\mathrm{p}^\\mathrm{Wilhoit}(T')}{R} (T')^4 \\ dT'
        
        evaluated at the given temperature `T` in kK.
        """
        cdef double Cp0, CpInf, B, a0, a1, a2, a3
        cdef double logBplusT, result
        Cp0, CpInf, B, a0, a1, a2, a3 = self._Cp0.value_si, self._CpInf.value_si, self._B.value_si, self.a0, self.a1, self.a2, self.a3
        logBplusT = log(B + T)
        result = ( -((5 + 10*a0 + 20*a1 + 35*a2 + 56*a3)*B**4*(Cp0 - CpInf)*T) + ((4 + 6*a0 + 10*a1 + 15*a2 + 21*a3)*B**3*(Cp0 - CpInf)*T**2)/2. +
            ((3 + 3*a0 + 4*a1 + 5*a2 + 6*a3)*B**2*(-Cp0 + CpInf)*T**3)/3. + ((2 + a0 + a1 + a2 + a3)*B*(Cp0 - CpInf)*T**4)/4. + (CpInf*T**5)/5. + (a3*B**10*(Cp0 - CpInf))/(5.*(B + T)**5) -
            ((a2 + 9*a3)*B**9*(Cp0 - CpInf))/(4.*(B + T)**4) + ((a1 + 8*a2 + 36*a3)*B**8*(Cp0 - CpInf))/(3.*(B + T)**3) - ((a0 + 7*(a1 + 4*(a2 + 3*a3)))*B**7*(Cp0 - CpInf))/(2.*(B + T)**2) +
            ((1 + 6*a0 + 21*a1 + 56*a2 + 126*a3)*B**6*(Cp0 - CpInf))/(B + T) + (6 + 15*a0 + 35*a1 + 70*a2 + 126*a3)*B**5*(Cp0 - CpInf)*logBplusT)
        return result
    
    cdef double integral2_T0(self, double T):
        """
        Return the value of the dimensionless integral
        
        .. math:: \\int \\left[ \\frac{C_\\mathrm{p}^\\mathrm{Wilhoit}(T')}{R} \\right]^2 \\ dT'
        
        evaluated at the given temperature `T` in kK.
        """
        cdef double Cp0, CpInf, B, a0, a1, a2, a3
        cdef double logBplusT, result
        Cp0, CpInf, B, a0, a1, a2, a3 = self._Cp0.value_si, self._CpInf.value_si, self._B.value_si, self.a0, self.a1, self.a2, self.a3
        logBplusT = log(B + T)
        result = (CpInf**2*T - (a3**2*B**12*(Cp0 - CpInf)**2)/(11.*(B + T)**11) + (a3*(a2 + 5*a3)*B**11*(Cp0 - CpInf)**2)/(5.*(B + T)**10) -
            ((a2**2 + 18*a2*a3 + a3*(2*a1 + 45*a3))*B**10*(Cp0 - CpInf)**2)/(9.*(B + T)**9) + ((4*a2**2 + 36*a2*a3 + a1*(a2 + 8*a3) + a3*(a0 + 60*a3))*B**9*(Cp0 - CpInf)**2)/(4.*(B + T)**8) -
            ((a1**2 + 14*a1*(a2 + 4*a3) + 2*(14*a2**2 + a3 + 84*a2*a3 + 105*a3**2 + a0*(a2 + 7*a3)))*B**8*(Cp0 - CpInf)**2)/(7.*(B + T)**7) +
            ((3*a1**2 + a2 + 28*a2**2 + 7*a3 + 126*a2*a3 + 126*a3**2 + 7*a1*(3*a2 + 8*a3) + a0*(a1 + 6*a2 + 21*a3))*B**7*(Cp0 - CpInf)**2)/(3.*(B + T)**6) -
            (B**6*(Cp0 - CpInf)*(a0**2*(Cp0 - CpInf) + 15*a1**2*(Cp0 - CpInf) + 10*a0*(a1 + 3*a2 + 7*a3)*(Cp0 - CpInf) + 2*a1*(1 + 35*a2 + 70*a3)*(Cp0 - CpInf) +
             2*(35*a2**2*(Cp0 - CpInf) + 6*a2*(1 + 21*a3)*(Cp0 - CpInf) + a3*(5*(4 + 21*a3)*Cp0 - 21*(CpInf + 5*a3*CpInf)))))/(5.*(B + T)**5) +
            (B**5*(Cp0 - CpInf)*(14*a2*Cp0 + 28*a2**2*Cp0 + 30*a3*Cp0 + 84*a2*a3*Cp0 + 60*a3**2*Cp0 + 2*a0**2*(Cp0 - CpInf) + 10*a1**2*(Cp0 - CpInf) +
             a0*(1 + 10*a1 + 20*a2 + 35*a3)*(Cp0 - CpInf) + a1*(5 + 35*a2 + 56*a3)*(Cp0 - CpInf) - 15*a2*CpInf - 28*a2**2*CpInf - 35*a3*CpInf - 84*a2*a3*CpInf - 60*a3**2*CpInf))/
             (2.*(B + T)**4) - (B**4*(Cp0 - CpInf)*((1 + 6*a0**2 + 15*a1**2 + 32*a2 + 28*a2**2 + 50*a3 + 72*a2*a3 + 45*a3**2 + 2*a1*(9 + 21*a2 + 28*a3) + a0*(8 + 20*a1 + 30*a2 + 42*a3))*Cp0 -
             (1 + 6*a0**2 + 15*a1**2 + 40*a2 + 28*a2**2 + 70*a3 + 72*a2*a3 + 45*a3**2 + a0*(8 + 20*a1 + 30*a2 + 42*a3) + a1*(20 + 42*a2 + 56*a3))*CpInf))/(3.*(B + T)**3) +
            (B**3*(Cp0 - CpInf)*((2 + 2*a0**2 + 3*a1**2 + 9*a2 + 4*a2**2 + 11*a3 + 9*a2*a3 + 5*a3**2 + a0*(5 + 5*a1 + 6*a2 + 7*a3) + a1*(7 + 7*a2 + 8*a3))*Cp0 -
             (2 + 2*a0**2 + 3*a1**2 + 15*a2 + 4*a2**2 + 21*a3 + 9*a2*a3 + 5*a3**2 + a0*(6 + 5*a1 + 6*a2 + 7*a3) + a1*(10 + 7*a2 + 8*a3))*CpInf))/(B + T)**2 -
            (B**2*((2 + a0 + a1 + a2 + a3)**2*Cp0**2 - 2*(5 + a0**2 + a1**2 + 8*a2 + a2**2 + 9*a3 + 2*a2*a3 + a3**2 + 2*a0*(3 + a1 + a2 + a3) + a1*(7 + 2*a2 + 2*a3))*Cp0*CpInf +
             (6 + a0**2 + a1**2 + 12*a2 + a2**2 + 14*a3 + 2*a2*a3 + a3**2 + 2*a1*(5 + a2 + a3) + 2*a0*(4 + a1 + a2 + a3))*CpInf**2))/(B + T) +
            2*(2 + a0 + a1 + a2 + a3)*B*(Cp0 - CpInf)*CpInf*logBplusT)
        return result
    
    cdef double integral2_TM1(self, double T):
        """
        Return the value of the dimensionless integral
        
        .. math:: \\int \\left[ \\frac{C_\\mathrm{p}^\\mathrm{Wilhoit}(T')}{R} \\right]^2 (T')^{-1} \\ dT'
        
        evaluated at the given temperature `T` in kK.
        """
        cdef double Cp0, CpInf, B, a0, a1, a2, a3
        cdef double logBplusT, result
        Cp0, CpInf, B, a0, a1, a2, a3 = self._Cp0.value_si, self._CpInf.value_si, self._B.value_si, self.a0, self.a1, self.a2, self.a3
        logBplusT = log(B + T); logT = log(T)
        result = ( (a3**2*B**11*(Cp0 - CpInf)**2)/(11.*(B + T)**11) - (a3*(2*a2 + 9*a3)*B**10*(Cp0 - CpInf)**2)/(10.*(B + T)**10) +
            ((a2**2 + 16*a2*a3 + 2*a3*(a1 + 18*a3))*B**9*(Cp0 - CpInf)**2)/(9.*(B + T)**9) -
            ((7*a2**2 + 56*a2*a3 + 2*a1*(a2 + 7*a3) + 2*a3*(a0 + 42*a3))*B**8*(Cp0 - CpInf)**2)/(8.*(B + T)**8) +
            ((a1**2 + 21*a2**2 + 2*a3 + 112*a2*a3 + 126*a3**2 + 2*a0*(a2 + 6*a3) + 6*a1*(2*a2 + 7*a3))*B**7*(Cp0 - CpInf)**2)/(7.*(B + T)**7) -
            ((5*a1**2 + 2*a2 + 30*a1*a2 + 35*a2**2 + 12*a3 + 70*a1*a3 + 140*a2*a3 + 126*a3**2 + 2*a0*(a1 + 5*(a2 + 3*a3)))*B**6*(Cp0 - CpInf)**2)/(6.*(B + T)**6) +
            (B**5*(Cp0 - CpInf)*(10*a2*Cp0 + 35*a2**2*Cp0 + 28*a3*Cp0 + 112*a2*a3*Cp0 + 84*a3**2*Cp0 + a0**2*(Cp0 - CpInf) + 10*a1**2*(Cp0 - CpInf) + 2*a1*(1 + 20*a2 + 35*a3)*(Cp0 - CpInf) +
            4*a0*(2*a1 + 5*(a2 + 2*a3))*(Cp0 - CpInf) - 10*a2*CpInf - 35*a2**2*CpInf - 30*a3*CpInf - 112*a2*a3*CpInf - 84*a3**2*CpInf))/(5.*(B + T)**5) -
            (B**4*(Cp0 - CpInf)*(18*a2*Cp0 + 21*a2**2*Cp0 + 32*a3*Cp0 + 56*a2*a3*Cp0 + 36*a3**2*Cp0 + 3*a0**2*(Cp0 - CpInf) + 10*a1**2*(Cp0 - CpInf) +
            2*a0*(1 + 6*a1 + 10*a2 + 15*a3)*(Cp0 - CpInf) + 2*a1*(4 + 15*a2 + 21*a3)*(Cp0 - CpInf) - 20*a2*CpInf - 21*a2**2*CpInf - 40*a3*CpInf - 56*a2*a3*CpInf - 36*a3**2*CpInf))/
            (4.*(B + T)**4) + (B**3*(Cp0 - CpInf)*((1 + 3*a0**2 + 5*a1**2 + 14*a2 + 7*a2**2 + 18*a3 + 16*a2*a3 + 9*a3**2 + 2*a0*(3 + 4*a1 + 5*a2 + 6*a3) + 2*a1*(5 + 6*a2 + 7*a3))*Cp0 -
            (1 + 3*a0**2 + 5*a1**2 + 20*a2 + 7*a2**2 + 30*a3 + 16*a2*a3 + 9*a3**2 + 2*a0*(3 + 4*a1 + 5*a2 + 6*a3) + 2*a1*(6 + 6*a2 + 7*a3))*CpInf))/(3.*(B + T)**3) -
            (B**2*((3 + a0**2 + a1**2 + 4*a2 + a2**2 + 4*a3 + 2*a2*a3 + a3**2 + 2*a1*(2 + a2 + a3) + 2*a0*(2 + a1 + a2 + a3))*Cp0**2 -
            2*(3 + a0**2 + a1**2 + 7*a2 + a2**2 + 8*a3 + 2*a2*a3 + a3**2 + 2*a1*(3 + a2 + a3) + a0*(5 + 2*a1 + 2*a2 + 2*a3))*Cp0*CpInf +
            (3 + a0**2 + a1**2 + 10*a2 + a2**2 + 12*a3 + 2*a2*a3 + a3**2 + 2*a1*(4 + a2 + a3) + 2*a0*(3 + a1 + a2 + a3))*CpInf**2))/(2.*(B + T)**2) +
            (B*(Cp0 - CpInf)*(Cp0 - (3 + 2*a0 + 2*a1 + 2*a2 + 2*a3)*CpInf))/(B + T) + Cp0**2*logT + (-Cp0**2 + CpInf**2)*logBplusT)
        return result

    cpdef ThermoData toThermoData(self):
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
            Cp0 = self.Cp0,
            CpInf = self.CpInf,
            E0 = self.E0,
        )
    
    cpdef NASA toNASA(self, double Tmin, double Tmax, double Tint, bint fixedTint=False, bint weighting=True, int continuity=3):
        """
        Convert the Wilhoit object to a :class:`NASA` object. You must specify
        the minimum and maximum temperatures of the fit `Tmin` and `Tmax` in K,
        as well as the intermediate temperature `Tint` in K to use as the bridge
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
        cdef Wilhoit wilhoit_scaled
        cdef NASAPolynomial nasa_low, nasa_high
        cdef double iseUnw, rmsUnw, iseWei, rmsWei, T
        cdef str rmsStr
        
        from rmgpy.thermo.nasa import NASA, NASAPolynomial

        # Scale the temperatures to kK
        Tmin /= 1000.
        Tint /= 1000.
        Tmax /= 1000.
        
        # Make copy of Wilhoit data so we don't modify the original
        wilhoit_scaled = self.copy() 
        # Rescale Wilhoit parameters
        wilhoit_scaled._Cp0.value_si /= constants.R
        wilhoit_scaled._CpInf.value_si /= constants.R
        wilhoit_scaled._B.value_si /= 1000.
        
        # If we are using fixed Tint, do not allow Tint to float
        if fixedTint:
            nasa_low, nasa_high = Wilhoit_to_NASA(wilhoit_scaled, Tmin, Tmax, Tint, weighting, continuity)
        else:
            nasa_low, nasa_high, Tint = Wilhoit_to_NASA_TintOpt(wilhoit_scaled, Tmin, Tmax, weighting, continuity)
        iseUnw = Wilhoit_to_NASA_TintOpt_objFun(Tint, wilhoit_scaled, Tmin, Tmax, 0, continuity) #the scaled, unweighted ISE (integral of squared error)
        rmsUnw = sqrt(iseUnw/(Tmax-Tmin))
        rmsStr = 'Unweighted RMS error = %.3f*R; '%(rmsUnw)
        if (weighting == 1):
            iseWei = Wilhoit_to_NASA_TintOpt_objFun(Tint, wilhoit_scaled, Tmin, Tmax, weighting, continuity) #the scaled, weighted ISE
            rmsWei = sqrt(iseWei/log(Tmax/Tmin))
            rmsStr = 'Weighted RMS error = %.3f*R; '%(rmsWei)+rmsStr
    
        # Print a warning if the RMS fit is worse that 0.25*R
        #if (rmsUnw > 0.25 or rmsWei > 0.25):
        #    print("Poor Wilhoit-to-NASA fit quality: RMS error = {0:.3f}*R".format(rmsWei if weighting == 1 else rmsUnw))
    
        # Restore to conventional units of K for Tint and units based on K rather than kK in NASA polynomial coefficients
        Tint *= 1000.
        Tmin *= 1000.
        Tmax *= 1000.
        nasa_low.c1 *= 1.0e-3
        nasa_low.c2 *= 1.0e-6
        nasa_low.c3 *= 1.0e-9
        nasa_low.c4 *= 1.0e-12
        nasa_high.c1 *= 1.0e-3
        nasa_high.c2 *= 1.0e-6
        nasa_high.c3 *= 1.0e-9
        nasa_high.c4 *= 1.0e-12
        
        # output comment
        # comment = 'NASA function fitted to Wilhoit function with B = {0:g} K. {1}\n{2}'.format(self.B.value_si, rmsStr, self.comment)
    
        # For the low polynomial, we want the results to match the Wilhoit value at 298 K
        nasa_low.c5 = (self.getEnthalpy(298) - nasa_low.getEnthalpy(298)) / constants.R
        nasa_low.c6 = (self.getEntropy(298) - nasa_low.getEntropy(298)) / constants.R
        
        # For the high polynomial, we want the results to match the low polynomial value at tint
        nasa_high.c5 = (nasa_low.getEnthalpy(Tint) - nasa_high.getEnthalpy(Tint)) / constants.R
        nasa_high.c6 = (nasa_low.getEntropy(Tint) - nasa_high.getEntropy(Tint)) / constants.R
    
        nasa = NASA(
            polynomials = [nasa_low, nasa_high],
            Tmin = nasa_low.Tmin,
            Tmax = nasa_high.Tmax,
            E0 = self.E0,
            comment = self.comment,
        )
    
        return nasa

################################################################################

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef Wilhoit_to_NASA(Wilhoit wilhoit, double Tmin, double Tmax, double Tint, bint weighting, int contCons):
    """
    Convert a Wilhoit polynomial to a pair of NASA polynomials.
    
    :param wilhoit: The Wilhoit polynomial to convert, with dimensionless heat
                    capacity limits and scaled temperature coefficient in kK
    :param Tmin:    The minimum temperature of the low-temperature NASA 
                    polynomial, in kK
    :param Tmax:    The maximum temperature of the high-temperature NASA 
                    polynomial, in kK
    :param Tint:    The intermediate temperature dividing the low-temperature
                    and high-temperature NASA polynomials, in kK
    :param weighting: ``True`` to weight the fit by inverse temperature, ``False`` to apply no weighting
    :param contCons:  The number of continuity constraints to apply to the 
                      fitted NASA polynomials at `Tint`:
        0: no constraints on continuity of Cp(T) at Tint
        1: constrain Cp to be continuous at Tint
        2: constrain Cp and dCp/dT to be continuous at Tint
        3 (default): constrain Cp, dCp/dT, and d2Cp/dT2 to be continuous at Tint
        4: constrain Cp, dCp/dT, d2Cp/dT2, and d3Cp/dT3 to be continuous at Tint
        5: constrain Cp, dCp/dT, d2Cp/dT2, d3Cp/dT3, and d4Cp/dT4 to be continuous at Tint; note: this effectively constrains all the coefficients to be equal and should be equivalent to fitting only one polynomial (rather than two)
        
        note: 5th (and higher) derivatives of NASA Cp(T) are zero and hence will automatically be continuous at Tint by the form of the Cp(T) function
                     
    :result: The pair of NASA polynomials with scaled parameters
    """
    cdef numpy.ndarray[numpy.float64_t, ndim=2] A
    cdef numpy.ndarray[numpy.float64_t, ndim=1] b, x
    cdef double w0min, w1min, w2min, w3min, w4min, wM1min
    cdef double w0int, w1int, w2int, w3int, w4int, wM1int
    cdef double w0max, w1max, w2max, w3max, w4max, wM1max
    cdef NASA nasa
    cdef int i, j
    
    #construct (typically 13*13) symmetric A matrix (in A*x = b); other elements will be zero
    A = numpy.zeros([10+contCons,10+contCons])
    b = numpy.zeros([10+contCons])

    if weighting:
        A[0,0] = 2*log(Tint/Tmin)
        A[0,1] = 2*(Tint - Tmin)
        A[0,2] = Tint*Tint - Tmin*Tmin
        A[0,3] = 2.*(Tint*Tint*Tint - Tmin*Tmin*Tmin)/3
        A[0,4] = (Tint*Tint*Tint*Tint - Tmin*Tmin*Tmin*Tmin)/2
        A[1,4] = 2.*(Tint*Tint*Tint*Tint*Tint - Tmin*Tmin*Tmin*Tmin*Tmin)/5
        A[2,4] = (Tint*Tint*Tint*Tint*Tint*Tint - Tmin*Tmin*Tmin*Tmin*Tmin*Tmin)/3
        A[3,4] = 2.*(Tint*Tint*Tint*Tint*Tint*Tint*Tint - Tmin*Tmin*Tmin*Tmin*Tmin*Tmin*Tmin)/7
        A[4,4] = (Tint*Tint*Tint*Tint*Tint*Tint*Tint*Tint - Tmin*Tmin*Tmin*Tmin*Tmin*Tmin*Tmin*Tmin)/4
    else:
        A[0,0] = 2*(Tint - Tmin)
        A[0,1] = Tint*Tint - Tmin*Tmin
        A[0,2] = 2.*(Tint*Tint*Tint - Tmin*Tmin*Tmin)/3
        A[0,3] = (Tint*Tint*Tint*Tint - Tmin*Tmin*Tmin*Tmin)/2
        A[0,4] = 2.*(Tint*Tint*Tint*Tint*Tint - Tmin*Tmin*Tmin*Tmin*Tmin)/5
        A[1,4] = (Tint*Tint*Tint*Tint*Tint*Tint - Tmin*Tmin*Tmin*Tmin*Tmin*Tmin)/3
        A[2,4] = 2.*(Tint*Tint*Tint*Tint*Tint*Tint*Tint - Tmin*Tmin*Tmin*Tmin*Tmin*Tmin*Tmin)/7
        A[3,4] = (Tint*Tint*Tint*Tint*Tint*Tint*Tint*Tint - Tmin*Tmin*Tmin*Tmin*Tmin*Tmin*Tmin*Tmin)/4
        A[4,4] = 2.*(Tint*Tint*Tint*Tint*Tint*Tint*Tint*Tint*Tint - Tmin*Tmin*Tmin*Tmin*Tmin*Tmin*Tmin*Tmin*Tmin)/9
    A[1,1] = A[0,2]
    A[1,2] = A[0,3]
    A[1,3] = A[0,4]
    A[2,2] = A[0,4]
    A[2,3] = A[1,4]
    A[3,3] = A[2,4]

    if weighting:
        A[5,5] = 2*log(Tmax/Tint)
        A[5,6] = 2*(Tmax - Tint)
        A[5,7] = Tmax*Tmax - Tint*Tint
        A[5,8] = 2.*(Tmax*Tmax*Tmax - Tint*Tint*Tint)/3
        A[5,9] = (Tmax*Tmax*Tmax*Tmax - Tint*Tint*Tint*Tint)/2
        A[6,9] = 2.*(Tmax*Tmax*Tmax*Tmax*Tmax - Tint*Tint*Tint*Tint*Tint)/5
        A[7,9] = (Tmax*Tmax*Tmax*Tmax*Tmax*Tmax - Tint*Tint*Tint*Tint*Tint*Tint)/3
        A[8,9] = 2.*(Tmax*Tmax*Tmax*Tmax*Tmax*Tmax*Tmax - Tint*Tint*Tint*Tint*Tint*Tint*Tint)/7
        A[9,9] = (Tmax*Tmax*Tmax*Tmax*Tmax*Tmax*Tmax*Tmax - Tint*Tint*Tint*Tint*Tint*Tint*Tint*Tint)/4
    else:
        A[5,5] = 2*(Tmax - Tint)
        A[5,6] = Tmax*Tmax - Tint*Tint
        A[5,7] = 2.*(Tmax*Tmax*Tmax - Tint*Tint*Tint)/3
        A[5,8] = (Tmax*Tmax*Tmax*Tmax - Tint*Tint*Tint*Tint)/2
        A[5,9] = 2.*(Tmax*Tmax*Tmax*Tmax*Tmax - Tint*Tint*Tint*Tint*Tint)/5
        A[6,9] = (Tmax*Tmax*Tmax*Tmax*Tmax*Tmax - Tint*Tint*Tint*Tint*Tint*Tint)/3
        A[7,9] = 2.*(Tmax*Tmax*Tmax*Tmax*Tmax*Tmax*Tmax - Tint*Tint*Tint*Tint*Tint*Tint*Tint)/7
        A[8,9] = (Tmax*Tmax*Tmax*Tmax*Tmax*Tmax*Tmax*Tmax - Tint*Tint*Tint*Tint*Tint*Tint*Tint*Tint)/4
        A[9,9] = 2.*(Tmax*Tmax*Tmax*Tmax*Tmax*Tmax*Tmax*Tmax*Tmax - Tint*Tint*Tint*Tint*Tint*Tint*Tint*Tint*Tint)/9
    A[6,6] = A[5,7]
    A[6,7] = A[5,8]
    A[6,8] = A[5,9]
    A[7,7] = A[5,9]
    A[7,8] = A[6,9]
    A[8,8] = A[7,9]

    if(contCons > 0):#set non-zero elements in the 11th column for Cp(T) continuity contraint
        A[0,10] = 1.
        A[1,10] = Tint
        A[2,10] = Tint*Tint
        A[3,10] = A[2,10]*Tint
        A[4,10] = A[3,10]*Tint
        A[5,10] = -A[0,10]
        A[6,10] = -A[1,10]
        A[7,10] = -A[2,10]
        A[8,10] = -A[3,10]
        A[9,10] = -A[4,10]
        if(contCons > 1): #set non-zero elements in the 12th column for dCp/dT continuity constraint
            A[1,11] = 1.
            A[2,11] = 2*Tint
            A[3,11] = 3*A[2,10]
            A[4,11] = 4*A[3,10]
            A[6,11] = -A[1,11]
            A[7,11] = -A[2,11]
            A[8,11] = -A[3,11]
            A[9,11] = -A[4,11]
            if(contCons > 2): #set non-zero elements in the 13th column for d2Cp/dT2 continuity constraint
                A[2,12] = 2.
                A[3,12] = 6*Tint
                A[4,12] = 12*A[2,10]
                A[7,12] = -A[2,12]
                A[8,12] = -A[3,12]
                A[9,12] = -A[4,12]
                if(contCons > 3): #set non-zero elements in the 14th column for d3Cp/dT3 continuity constraint
                    A[3,13] = 6
                    A[4,13] = 24*Tint
                    A[8,13] = -A[3,13]
                    A[9,13] = -A[4,13]
                    if(contCons > 4): #set non-zero elements in the 15th column for d4Cp/dT4 continuity constraint
                        A[4,14] = 24
                        A[9,14] = -A[4,14]

    # make the matrix symmetric
    for i in range(1,10+contCons):
        for j in range(0, i):
            A[i,j] = A[j,i]

    #construct b vector
    w0int = wilhoit.integral_T0(Tint)
    w1int = wilhoit.integral_T1(Tint)
    w2int = wilhoit.integral_T2(Tint)
    w3int = wilhoit.integral_T3(Tint)
    w0min = wilhoit.integral_T0(Tmin)
    w1min = wilhoit.integral_T1(Tmin)
    w2min = wilhoit.integral_T2(Tmin)
    w3min = wilhoit.integral_T3(Tmin)
    w0max = wilhoit.integral_T0(Tmax)
    w1max = wilhoit.integral_T1(Tmax)
    w2max = wilhoit.integral_T2(Tmax)
    w3max = wilhoit.integral_T3(Tmax)
    if weighting:
        wM1int = wilhoit.integral_TM1(Tint)
        wM1min = wilhoit.integral_TM1(Tmin)
        wM1max = wilhoit.integral_TM1(Tmax)
    else:
        w4int = wilhoit.integral_T4(Tint)
        w4min = wilhoit.integral_T4(Tmin)
        w4max = wilhoit.integral_T4(Tmax)

    if weighting:
        b[0] = 2*(wM1int - wM1min)
        b[1] = 2*(w0int - w0min)
        b[2] = 2*(w1int - w1min)
        b[3] = 2*(w2int - w2min)
        b[4] = 2*(w3int - w3min)
        b[5] = 2*(wM1max - wM1int)
        b[6] = 2*(w0max - w0int)
        b[7] = 2*(w1max - w1int)
        b[8] = 2*(w2max - w2int)
        b[9] = 2*(w3max - w3int)
    else:
        b[0] = 2*(w0int - w0min)
        b[1] = 2*(w1int - w1min)
        b[2] = 2*(w2int - w2min)
        b[3] = 2*(w3int - w3min)
        b[4] = 2*(w4int - w4min)
        b[5] = 2*(w0max - w0int)
        b[6] = 2*(w1max - w1int)
        b[7] = 2*(w2max - w2int)
        b[8] = 2*(w3max - w3int)
        b[9] = 2*(w4max - w4int)

    # solve A*x=b for x (note that factor of 2 in b vector and 10*10 submatrix of A
    # matrix is not required; not including it should give same result, except
    # Lagrange multipliers will differ by a factor of two)
    x = scipy.linalg.solve(A,b,overwrite_a=1,overwrite_b=1)

    nasa_low = NASAPolynomial(
        [x[0], x[1], x[2], x[3], x[4], 0.0, 0.0],
        Tmin = (Tmin * 1000.,"K"),
        Tmax = (Tint * 1000.,"K"),
    )
    nasa_high = NASAPolynomial(
        [x[5], x[6], x[7], x[8], x[9], 0.0, 0.0],
        Tmin = (Tint * 1000.,"K"),
        Tmax = (Tmax * 1000.,"K"),
    )

    return nasa_low, nasa_high

cpdef Wilhoit_to_NASA_TintOpt(Wilhoit wilhoit, double Tmin, double Tmax, bint weighting, int contCons):
    """
    Convert a Wilhoit polynomial to a pair of NASA polynomials, using an
    optimization algorithm to choose the best value of the intermediate
    temperature. The parameters are the same as for the :func:`Wilhoit_to_NASA`
    function.
    """
    import scipy.optimize
    Tint = scipy.optimize.fminbound(Wilhoit_to_NASA_TintOpt_objFun, Tmin, Tmax, args=(wilhoit, Tmin, Tmax, weighting, contCons))
    Tint = float(Tint) # fminbound returns a numpy.ndarray object
    #note that we have not used any guess when using this minimization routine
    #2. determine the bi parameters based on the optimized Tint (alternatively, maybe we could have Wilhoit_to_NASA_TintOpt_objFun also return these parameters, along with the objective function, which would avoid an extra calculation)
    nasa_low, nasa_high = Wilhoit_to_NASA(wilhoit, Tmin, Tmax, Tint, weighting, contCons)
    return nasa_low, nasa_high, Tint

cpdef double Wilhoit_to_NASA_TintOpt_objFun(double Tint, Wilhoit wilhoit, double Tmin, double Tmax, bint weighting, int contCons):
    """
    Evaluate the objective function used to convert a Wilhoit polynomial to a
    pair of NASA polynomials. The parameters are the same as for the 
    :func:`Wilhoit_to_NASA` function.
    """
    if (weighting == 1):
        result = Wilhoit_to_NASA_TintOpt_objFun_W(Tint, wilhoit, Tmin, Tmax, contCons)
    else:
        result = Wilhoit_to_NASA_TintOpt_objFun_NW(Tint, wilhoit, Tmin, Tmax, contCons)

    # numerical errors could accumulate to give a slightly negative result
    # this is unphysical (it's the integral of a *squared* error) so we
    # set it to zero to avoid later problems when we try find the square root.
    if result < 0:
        #print("Negative ISE of {0:g} reset to zero.".format(result))
        result = 0

    return result

cpdef double Wilhoit_to_NASA_TintOpt_objFun_NW(double Tint, Wilhoit wilhoit, double Tmin, double Tmax, int contCons):
    """
    Evaluate the unweighted objective function used to convert a Wilhoit 
    polynomial to a pair of NASA polynomials. The parameters are the same as 
    for the :func:`Wilhoit_to_NASA` function.
    """
    cdef NASAPolynomial nasa_low, nasa_high
    cdef double b1, b2, b3, b4, b5, b6, b7, b8, b9, b10
    cdef double qM1, q0, q1, q2, q3, result

    nasa_low, nasa_high = Wilhoit_to_NASA(wilhoit,Tmin,Tmax,Tint, 0, contCons)
    b1 = nasa_low.c0
    b2 = nasa_low.c1
    b3 = nasa_low.c2
    b4 = nasa_low.c3
    b5 = nasa_low.c4
    b6 = nasa_high.c0
    b7 = nasa_high.c1
    b8 = nasa_high.c2
    b9 = nasa_high.c3
    b10 = nasa_high.c4

    q0 = wilhoit.integral_T0(Tint)
    q1 = wilhoit.integral_T1(Tint)
    q2 = wilhoit.integral_T2(Tint)
    q3 = wilhoit.integral_T3(Tint)
    q4 = wilhoit.integral_T4(Tint)
    result = (wilhoit.integral2_T0(Tmax) - wilhoit.integral2_T0(Tmin) +
                 nasa_low.integral2_T0(Tint) - nasa_low.integral2_T0(Tmin) +
                 nasa_high.integral2_T0(Tmax) - nasa_high.integral2_T0(Tint)
                 - 2* (b6*(wilhoit.integral_T0(Tmax)-q0)+b1*(q0-wilhoit.integral_T0(Tmin))
                 +b7*(wilhoit.integral_T1(Tmax) - q1) +b2*(q1 - wilhoit.integral_T1(Tmin))
                 +b8*(wilhoit.integral_T2(Tmax) - q2) +b3*(q2 - wilhoit.integral_T2(Tmin))
                 +b9*(wilhoit.integral_T3(Tmax) - q3) +b4*(q3 - wilhoit.integral_T3(Tmin))
                 +b10*(wilhoit.integral_T4(Tmax) - q4)+b5*(q4 - wilhoit.integral_T4(Tmin))))

    return result

cpdef double Wilhoit_to_NASA_TintOpt_objFun_W(double Tint, Wilhoit wilhoit, double Tmin, double Tmax, int contCons):
    """
    Evaluate the weighted objective function used to convert a Wilhoit 
    polynomial to a pair of NASA polynomials. The parameters are the same as 
    for the :func:`Wilhoit_to_NASA` function. The weighting is by inverse
    temperature, to bias the fit towards the lower temperatures, where the
    heat capacity is changing more rapidly.
    
    If the fit is close to perfect, the result may be slightly negative due to 
    numerical errors in evaluating this integral.
    """
    cdef NASAPolynomial nasa_low, nasa_high
    cdef double b1, b2, b3, b4, b5, b6, b7, b8, b9, b10
    cdef double qM1, q0, q1, q2, q3, result
    
    nasa_low, nasa_high = Wilhoit_to_NASA(wilhoit,Tmin,Tmax,Tint, 1, contCons)
    b1 = nasa_low.c0
    b2 = nasa_low.c1
    b3 = nasa_low.c2
    b4 = nasa_low.c3
    b5 = nasa_low.c4
    b6 = nasa_high.c0
    b7 = nasa_high.c1
    b8 = nasa_high.c2
    b9 = nasa_high.c3
    b10 = nasa_high.c4

    qM1 = wilhoit.integral_TM1(Tint)
    q0 = wilhoit.integral_T0(Tint)
    q1 = wilhoit.integral_T1(Tint)
    q2 = wilhoit.integral_T2(Tint)
    q3 = wilhoit.integral_T3(Tint)
    result = (wilhoit.integral2_TM1(Tmax) - wilhoit.integral2_TM1(Tmin) +
                 nasa_low.integral2_TM1(Tint) - nasa_low.integral2_TM1(Tmin) +
                 nasa_high.integral2_TM1(Tmax) - nasa_high.integral2_TM1(Tint)
                 - 2* (b6*(wilhoit.integral_TM1(Tmax)-qM1)+b1*(qM1 - wilhoit.integral_TM1(Tmin))
                 +b7*(wilhoit.integral_T0(Tmax)-q0)+b2*(q0 - wilhoit.integral_T0(Tmin))
                 +b8*(wilhoit.integral_T1(Tmax)-q1)+b3*(q1 - wilhoit.integral_T1(Tmin))
                 +b9*(wilhoit.integral_T2(Tmax)-q2)+b4*(q2 - wilhoit.integral_T2(Tmin))
                 +b10*(wilhoit.integral_T3(Tmax)-q3)+b5*(q3 - wilhoit.integral_T3(Tmin))))

    return result
