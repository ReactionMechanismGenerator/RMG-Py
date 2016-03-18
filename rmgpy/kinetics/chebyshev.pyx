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
from libc.math cimport exp, log, sqrt, log10

cimport rmgpy.constants as constants
import rmgpy.quantity as quantity

################################################################################

cdef class Chebyshev(PDepKineticsModel):
    """
    A model of a phenomenological rate coefficient :math:`k(T,P)` using a
    set of Chebyshev polynomials in temperature and pressure. The attributes
    are:
    
    =============== ============================================================
    Attribute       Description
    =============== ============================================================
    `coeffs`        Matrix of Chebyshev coefficients, such that the resulting :math:`k(T,P)` has units of cm^3, mol, s
    `kunits`        The units of the rate coefficient
    `degreeT`       The number of terms in the inverse temperature direction
    `degreeP`       The number of terms in the log pressure direction
    `Tmin`          The minimum temperature at which the model is valid, or zero if unknown or undefined
    `Tmax`          The maximum temperature at which the model is valid, or zero if unknown or undefined
    `Pmin`          The minimum pressure at which the model is valid, or zero if unknown or undefined
    `Pmax`          The maximum pressure at which the model is valid, or zero if unknown or undefined
    `comment`       Information about the model (e.g. its source)
    =============== ============================================================
    
    """

    def __init__(self, coeffs=None, kunits='', highPlimit=None, Tmin=None, Tmax=None, Pmin=None, Pmax=None, comment=''):
        PDepKineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, highPlimit=highPlimit, comment=comment)
        self.coeffs = coeffs
        self.kunits = kunits
        if self.coeffs is not None:
            self.degreeT = self._coeffs.value_si.shape[0]
            self.degreeP = self._coeffs.value_si.shape[1]
            factor = quantity.RateCoefficient(1.0, kunits).value_si
            self.coeffs.value_si[0,0] += log10(factor)
        
    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        Chebyshev object.
        """
        coeffs = self.coeffs.copy()
        factor = quantity.RateCoefficient(1.0, self.kunits).value_si
        coeffs.value_si[0,0] -= log10(factor)
        string = 'Chebyshev(coeffs={0!r}, kunits={1!r}'.format(coeffs, self.kunits)
        if self.highPlimit is not None: string += ', highPlimit={0!r}'.format(self.highPlimit)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.Pmin is not None: string += ', Pmin={0!r}'.format(self.Pmin)
        if self.Pmax is not None: string += ', Pmax={0!r}'.format(self.Pmax)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling a Chebyshev object.
        """
        coeffs = self.coeffs.copy()
        factor = quantity.RateCoefficient(1.0, self.kunits).value_si
        coeffs.value_si[0,0] -= log10(factor)
        return (Chebyshev, (coeffs, self.kunits, self.highPlimit, self.Tmin, self.Tmax, self.Pmin, self.Pmax, self.comment))
    
    property coeffs:
        """The Chebyshev coefficients."""
        def __get__(self):
            return self._coeffs
        def __set__(self, value):
            self._coeffs = quantity.Dimensionless(value)
    
    cdef double chebyshev(self, int n, double x):
        """
        Return the value of the nth-order Chebyshev polynomial at the given
        value of `x`.
        """
        cdef double T0, T1, T
        cdef int i
        if n == 0:
            return 1
        elif n == 1:
            return x
        else:
            T0 = 1
            T1 = x
            for i in range(1,n):
                T = 2 * x * T1 - T0
                T0 = T1
                T1 = T
            return T

    cdef double getReducedTemperature(self, double T) except -1000:
        """
        Return the reduced temperature corresponding to the given temperature
        `T` in K. This maps the inverse of the temperature onto the domain 
        [-1, 1] using the `Tmin` and `Tmax` attributes as the limits.
        """
        cdef double Tmin, Tmax
        Tmin = self._Tmin.value_si
        Tmax = self._Tmax.value_si
        return (2.0/T - 1.0/Tmin - 1.0/Tmax) / (1.0/Tmax - 1.0/Tmin)
    
    cdef double getReducedPressure(self, double P) except -1000:
        """
        Return the reduced pressure corresponding to the given pressure
        `P` in Pa. This maps the logarithm of the pressure onto the domain 
        [-1, 1] using the `Pmin` and `Pmax` attributes as the limits.
        """
        cdef double Pmin, Pmax
        Pmin = self._Pmin.value_si
        Pmax = self._Pmax.value_si
        return (2.0*log10(P) - log10(Pmin) - log10(Pmax)) / (log10(Pmax) - log10(Pmin))
    
    cpdef double getRateCoefficient(self, double T, double P=0) except -1:
        """
        Return the rate coefficient in the appropriate combination of m^3, 
        mol, and s at temperature `T` in K and pressure `P` in Pa by 
        evaluating the Chebyshev expression.
        """
        cdef numpy.ndarray[numpy.float64_t,ndim=2] coeffs
        cdef double Tred, Pred, k
        cdef int i, j, t, p
        
        if P == 0:
            raise ValueError('No pressure specified to pressure-dependent Chebyshev.getRateCoefficient().')

        coeffs = self._coeffs.value_si

        k = 0.0
        Tred = self.getReducedTemperature(T)
        Pred = self.getReducedPressure(P)
        for t in range(self.degreeT):
            for p in range(self.degreeP):
                k += coeffs[t,p] * self.chebyshev(t, Tred) * self.chebyshev(p, Pred)
        return 10.0**k

    cpdef fitToData(self, numpy.ndarray Tlist, numpy.ndarray Plist, numpy.ndarray K,
        str kunits, int degreeT, int degreeP, double Tmin, double Tmax, double Pmin, double Pmax):
        """
        Fit a Chebyshev kinetic model to a set of rate coefficients `K`, which
        is a matrix corresponding to the temperatures `Tlist` in K and pressures
        `Plist` in Pa. `degreeT` and `degreeP` are the degree of the 
        polynomials in temperature and pressure, while `Tmin`, `Tmax`, `Pmin`,
        and `Pmax` set the edges of the valid temperature and pressure ranges
        in K and bar, respectively.
        """
        cdef int nT = len(Tlist), nP = len(Plist)
        cdef list Tred, Pred
        cdef int t1, p1, t2, p2
        cdef double T, P

        # Set temperature and pressure ranges
        self.Tmin = (Tmin,"K")
        self.Tmax = (Tmax,"K")
        self.Pmin = (Pmin*1e-5,"bar")
        self.Pmax = (Pmax*1e-5,"bar")

        # Calculate reduced temperatures and pressures
        Tred = [self.getReducedTemperature(T) for T in Tlist]
        Pred = [self.getReducedPressure(P) for P in Plist]

        K = quantity.RateCoefficient(K,kunits).value_si

        # Create matrix and vector for coefficient fit (linear least-squares)
        A = numpy.zeros((nT*nP, degreeT*degreeP), numpy.float64)
        b = numpy.zeros((nT*nP), numpy.float64)
        for t1, T in enumerate(Tred):
            for p1, P in enumerate(Pred):
                for t2 in range(degreeT):
                    for p2 in range(degreeP):
                        A[p1*nT+t1, p2*degreeT+t2] = self.chebyshev(t2, T) * self.chebyshev(p2, P)
                b[p1*nT+t1] = log10(K[t1,p1])

        # Do linear least-squares fit to get coefficients
        x, residues, rank, s = numpy.linalg.lstsq(A, b)

        # Extract coefficients
        coeffs = numpy.zeros((degreeT,degreeP), numpy.float64)
        for t2 in range(degreeT):
            for p2 in range(degreeP):
                coeffs[t2,p2] = x[p2*degreeT+t2]
        self.coeffs = coeffs

        self.degreeT = degreeT
        self.degreeP = degreeP
        
        self.kunits = kunits
        
        return self

    cpdef bint isIdenticalTo(self, KineticsModel otherKinetics) except -2:
        """
        Checks to see if kinetics matches that of other kinetics and returns ``True``
        if coeffs, kunits, Tmin,
        """
        if not isinstance(otherKinetics, Chebyshev):
            return False
        if not KineticsModel.isIdenticalTo(self,otherKinetics):
            return False
        if self.degreeT != otherKinetics.degreeT or self.degreeP != otherKinetics.degreeP:
            return False
        if self.kunits != otherKinetics.kunits or (self.coeffs != otherKinetics.coeffs).any():
            return False

        return True

    cpdef changeRate(self, double factor):
        """ 
        Changes kinetics rates by a multiple ``factor``.
        """
        self.coeffs.value_si[0,0] += log10(factor)

    def setCanteraKinetics(self, ctReaction, speciesList=[]):
        """
        Sets the kinetics parameters for a Cantera ChebyshevReaction() object
        Uses set_parameters(self,Tmin,Tmax,Pmin,Pmax,coeffs)
        where T's are in units of K, P's in units of Pa, and coeffs is 2D array of (nTemperature, nPressure).
        """
        import cantera as ct
        import copy
        assert isinstance(ctReaction, ct.ChebyshevReaction), "Must be a Cantera ChebyshevReaction object"

        Tmin = self.Tmin.value_si
        Tmax = self.Tmax.value_si
        Pmin = self.Pmin.value_si
        Pmax = self.Pmax.value_si
        coeffs = copy.deepcopy(self.coeffs.value_si)
        
        # The first coefficient must be adjusted to match Cantera units
        rateUnitsDimensionality = {'1/s':0,
                                   's^-1':0, 
                               'm^3/(mol*s)':1,
                               'm^6/(mol^2*s)':2,
                               'cm^3/(mol*s)':1,
                               'cm^6/(mol^2*s)':2,
                               'm^3/(molecule*s)': 1, 
                               'm^6/(molecule^2*s)': 2,
                               'cm^3/(molecule*s)': 1,
                               'cm^6/(molecule^2*s)': 2,
                               }
        try:
            factor = 1000**rateUnitsDimensionality[self.kunits]
            coeffs[0,0] += log10(factor)
        except:
            raise Exception('Chebyshev units {0} not found among accepted units for converting to Cantera Chebyshev object.'.format(self.kunits))
        ctReaction.set_parameters(Tmin, Tmax, Pmin, Pmax, coeffs)