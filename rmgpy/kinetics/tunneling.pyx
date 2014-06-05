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
This module contains classes representing various models of tunneling through
a reaction barrier.
"""

import numpy

import cython
from libc.math cimport abs, exp, sqrt, cosh

cimport rmgpy.constants as constants
import rmgpy.quantity as quantity

################################################################################

cdef class Wigner(TunnelingModel):
    """
    A tunneling model based on the Wigner formula. The attributes are:

    =============== =============================================================
    Attribute       Description
    =============== =============================================================
    `frequency`     The imaginary frequency of the transition state
    =============== =============================================================
    
    """
    
    def __init__(self, frequency):
        TunnelingModel.__init__(self, frequency)
    
    def __repr__(self):
        """
        Return a string representation of the tunneling model.
        """
        return 'Wigner(frequency={0!r})'.format(self.frequency)

    def __reduce__(self):
        """
        A helper function used when pickling an Wigner object.
        """
        return (Wigner, (self.frequency,))
    
    cpdef double calculateTunnelingFactor(self, double T) except -100000000:
        """
        Calculate and return the value of the Wigner tunneling correction for
        the reaction at the temperature `T` in K.
        """
        cdef double factor, frequency
        frequency = abs(self._frequency.value_si) * constants.c * 100.0
        factor = constants.h * frequency / (constants.kB * T)
        return 1.0 + factor * factor / 24.0
    
    cpdef numpy.ndarray calculateTunnelingFunction(self, numpy.ndarray Elist):
        """
        Raises :class:`NotImplementedError`, as the Wigner tunneling model
        does not have a well-defined energy-dependent tunneling function. 
        """
        raise NotImplementedError('The Wigner tunneling correction does not have a well-defined k(E) function.')

################################################################################

cdef class Eckart(TunnelingModel):
    """
    A tunneling model based on the Eckart model. The attributes are:

    =============== =============================================================
    Attribute       Description
    =============== =============================================================
    `frequency`     The imaginary frequency of the transition state
    `E0_reac`       The ground-state energy of the reactants
    `E0_TS`         The ground-state energy of the transition state
    `E0_prod`       The ground-state energy of the products
    =============== =============================================================
    
    If `E0_prod` is not given, it is assumed to be the same as the reactants;
    this results in the so-called "symmetric" Eckart model. Providing 
    `E0_prod`, and thereby using the "asymmetric" Eckart model, is the
    recommended approach.
    """
    
    def __init__(self, frequency, E0_reac, E0_TS, E0_prod=None):
        TunnelingModel.__init__(self, frequency)
        self.E0_reac = E0_reac
        self.E0_TS = E0_TS
        self.E0_prod = E0_prod
    
    def __repr__(self):
        """
        Return a string representation of the tunneling model.
        """
        return 'Eckart(frequency={0!r}, E0_reac={1!r}, E0_TS={2!r}, E0_prod={3!r})'.format(self.frequency, self.E0_reac, self.E0_TS, self.E0_prod)

    def __reduce__(self):
        """
        A helper function used when pickling an Eckart object.
        """
        return (Eckart, (self.frequency, self.E0_reac, self.E0_TS, self.E0_prod))
    
    property E0_reac:
        """The ground-state energy of the reactants."""
        def __get__(self):
            return self._E0_reac
        def __set__(self, value):
            self._E0_reac = quantity.Energy(value)

    property E0_TS:
        """The ground-state energy of the transition state."""
        def __get__(self):
            return self._E0_TS
        def __set__(self, value):
            self._E0_TS = quantity.Energy(value)

    property E0_prod:
        """The ground-state energy of the products."""
        def __get__(self):
            return self._E0_prod
        def __set__(self, value):
            self._E0_prod = quantity.Energy(value)

    cpdef double calculateTunnelingFactor(self, double T) except -100000000:
        """
        Calculate and return the value of the Eckart tunneling correction for
        the reaction at the temperature `T` in K.
        """
        cdef double E0_reac, E0_prod, E0_TS
        cdef double E0, dE, beta, dV1, dV2
        cdef numpy.ndarray Elist, kappaE
        
        beta = 1. / (constants.R * T)        # [=] mol/J
        
        E0_reac = self._E0_reac.value_si
        E0_TS = self._E0_TS.value_si
        E0_prod = self._E0_prod.value_si
        
        # Calculate intermediate constants
        if E0_reac > E0_prod:
            E0 = E0_reac
            dV1 = E0_TS - E0_reac
            dV2 = E0_TS - E0_prod
        else:
            E0 = E0_prod
            dV1 = E0_TS - E0_prod
            dV2 = E0_TS - E0_reac

        if dV1 < 0 or dV2 < 0:
            raise ValueError('One or both of the barrier heights of {0:g} and {1:g} kJ/mol encountered in Eckart method are invalid.'.format(dV1 / 1000., dV2 / 1000.)) 

        # Ensure that dV1 is smaller than dV2
        assert dV1 <= dV2

        # Evaluate microcanonical tunneling function kappa(E)
        dE = 100.
        Elist = numpy.arange(E0, E0 + 2. * (E0_TS - E0) + 40. * constants.R * T, dE)
        kappaE = self.calculateTunnelingFunction(Elist)
        
        # Integrate to get kappa(T)
        kappa = exp(dV1 * beta) * numpy.sum(kappaE * numpy.exp(-beta * (Elist - E0))) * dE * beta
        
        # Return the calculated Eckart correction
        return kappa
    
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef numpy.ndarray calculateTunnelingFunction(self, numpy.ndarray Elist):
        """
        Calculate and return the value of the Eckart tunneling function for
        the reaction at the energies `Elist` in J/mol.
        """
        cdef double frequency, E0_reac, E0_prod, E0_TS
        cdef numpy.ndarray[numpy.float64_t,ndim=1] kappa, _Elist
        cdef double E0, dV1, dV2, alpha1, alpha2, E, xi, twopia, twopib, twopid
        cdef int r, r0
        
        frequency = abs(self._frequency.value_si) * constants.h * constants.c * 100. * constants.Na
        E0_reac = self._E0_reac.value_si
        E0_TS = self._E0_TS.value_si
        E0_prod = self._E0_prod.value_si

        _Elist = Elist
        
        # Calculate intermediate constants
        if E0_reac > E0_prod:
            E0 = E0_reac
            dV1 = E0_TS - E0_reac
            dV2 = E0_TS - E0_prod
        else:
            E0 = E0_prod
            dV1 = E0_TS - E0_prod
            dV2 = E0_TS - E0_reac
        
        # Ensure that dV1 is smaller than dV2
        assert dV1 <= dV2
        
        alpha1 = 2 * constants.pi * dV1 / frequency
        alpha2 = 2 * constants.pi * dV2 / frequency

        kappa = numpy.zeros_like(Elist)
        for r0 in range(_Elist.shape[0]):
            if _Elist[r0] >= E0:
                break
        
        for r in range(r0, _Elist.shape[0]):
            E = _Elist[r]
            
            xi = (E - E0) / dV1
            # 2 * pi * a
            twopia = 2.*sqrt(alpha1*xi)/(1./sqrt(alpha1)+1./sqrt(alpha2))
            # 2 * pi * b
            twopib = 2.*sqrt(abs((xi-1.)*alpha1+alpha2))/(1./sqrt(alpha1)+1/sqrt(alpha2))
            # 2 * pi * d
            twopid = 2.*sqrt(abs(alpha1*alpha2-4*constants.pi*constants.pi/16.))
            
            # We use different approximate versions of the integrand to avoid
            # domain errors when evaluating cosh(x) for large x
            # If all of 2*pi*a, 2*pi*b, and 2*pi*d are sufficiently small,
            # compute as normal
            if twopia < 200. and twopib < 200. and twopid < 200.:
                kappa[r] = 1 - (cosh(twopia-twopib)+cosh(twopid)) / (cosh(twopia+twopib)+cosh(twopid))
            # If one of the following is true, then we can eliminate most of the
            # exponential terms after writing out the definition of cosh and
            # dividing all terms by exp(2*pi*d)
            elif twopia-twopib-twopid > 10 or twopib-twopia-twopid > 10 or twopia+twopib-twopid > 10:
                kappa[r] = 1 - exp(-2*twopia) - exp(-2*twopib) - exp(-twopia-twopib+twopid) - exp(-twopia-twopib-twopid)
            # Otherwise expand each cosh(x) in terms of its exponentials and divide
            # all terms by exp(2*pi*d) before evaluating
            else:
                kappa[r] = 1 - (exp(twopia-twopib-twopid) + exp(-twopia+twopib-twopid) + 1 + exp(-2*twopid)) / (exp(twopia+twopib-twopid) + exp(-twopia-twopib-twopid) + 1 + exp(-2*twopid))
        
        return kappa
