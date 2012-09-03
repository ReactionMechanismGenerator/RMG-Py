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

cdef class Arrhenius(KineticsModel):
    """
    A kinetics model based on the (modified) Arrhenius equation. The attributes
    are:

    =============== =============================================================
    Attribute       Description
    =============== =============================================================
    `A`             The preexponential factor
    `T0`            The reference temperature
    `n`             The temperature exponent
    `Ea`            The activation energy
    `Tmin`          The minimum temperature at which the model is valid, or zero if unknown or undefined
    `Tmax`          The maximum temperature at which the model is valid, or zero if unknown or undefined
    `comment`       Information about the model (e.g. its source)
    =============== =============================================================
    
    """
    
    def __init__(self, A=None, n=0.0, Ea=None, T0=(1.0,"K"), Tmin=None, Tmax=None, comment=''):
        KineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, comment=comment)
        self.A = A
        self.n = n
        self.Ea = Ea
        self.T0 = T0
        
    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        Arrhenius object.
        """
        string = 'Arrhenius(A={0!r}, n={1:g}, Ea={2!r}, T0={3!r}'.format(self.A, self.n, self.Ea, self.T0)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling an Arrhenius object.
        """
        return (Arrhenius, (self.A, self.n, self.Ea, self.T0, self.Tmin, self.Tmax, self.comment))

    property A:
        """The preexponential factor."""
        def __get__(self):
            return self._A
        def __set__(self, value):
            self._A = quantity.RateCoefficient(value)

    property Ea:
        """The activation energy."""
        def __get__(self):
            return self._Ea
        def __set__(self, value):
            self._Ea = quantity.Energy(value)

    property T0:
        """The reference temperature."""
        def __get__(self):
            return self._T0
        def __set__(self, value):
            self._T0 = quantity.Temperature(value)

    cpdef double getRateCoefficient(self, double T, double P=0.0) except -1:
        """
        Return the rate coefficient in the appropriate combination of m^3, 
        mol, and s at temperature `T` in K. 
        """
        cdef double A, n, Ea, T0
        A = self._A.value_si
        n = self.n
        Ea = self._Ea.value_si
        T0 = self._T0.value_si
        return A * (T / T0)**n * exp(-Ea / (constants.R * T))

    cpdef changeT0(self, double T0):
        """
        Changes the reference temperature used in the exponent to `T0` in K, 
        and adjusts the preexponential factor accordingly.
        """
        self._A.value_si /= (self._T0.value_si / T0)**self.n
        self._T0.value_si = T0

    cpdef fitToData(self, numpy.ndarray Tlist, numpy.ndarray klist, str kunits, double T0=1, numpy.ndarray weights=None, bint threeParams=True):
        """
        Fit the Arrhenius parameters to a set of rate coefficient data `klist`
        in units of `kunits` corresponding to a set of temperatures `Tlist` in 
        K. A linear least-squares fit is used, which guarantees that the 
        resulting parameters provide the best possible approximation to the 
        data.
        """
        import numpy.linalg
        import scipy.stats
        if threeParams:
            A = numpy.zeros((len(Tlist),3), numpy.float64)
            A[:,0] = numpy.ones_like(Tlist)
            A[:,1] = numpy.log(Tlist / T0)
            A[:,2] = -1.0 / constants.R / Tlist
        else:
            A = numpy.zeros((len(Tlist),2), numpy.float64)
            A[:,0] = numpy.ones_like(Tlist)
            A[:,1] = -1.0 / constants.R / Tlist
        b = numpy.log(klist)
        if weights is not None:
            for n in range(b.size):
                A[n,:] *= weights[n]
                b[n] *= weights[n]
        x, residues, rank, s = numpy.linalg.lstsq(A,b)
        
        # Determine covarianace matrix to obtain parameter uncertainties
        count = klist.size
        cov = residues[0] / (count - 3) * numpy.linalg.inv(numpy.dot(A.T, A))
        t = scipy.stats.t.ppf(0.975, count - 3)
            
        if not threeParams:
            x = numpy.array([x[0], 0, x[1]])
            cov = numpy.array([[cov[0,0], 0, cov[0,1]], [0,0,0], [cov[1,0], 0, cov[1,1]]])
        
        self.A = (exp(x[0]),kunits)
        self.n = x[1]
        self.Ea = (x[2] * 0.001,"kJ/mol")
        self.T0 = (T0,"K")
        self.Tmin = (numpy.min(Tlist),"K")
        self.Tmax = (numpy.max(Tlist),"K")
        self.comment = 'Fitted to {0:d} data points; dA = *|/ {1:g}, dn = +|- {2:g}, dEa = +|- {3:g} kJ/mol'.format(
            len(Tlist),
            exp(sqrt(cov[0,0])),
            sqrt(cov[1,1]),
            sqrt(cov[2,2]) * 0.001,
        )
        
        return self
