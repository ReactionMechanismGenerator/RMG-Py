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
This module contains classes and functions for working with collision models.
"""

import numpy

cimport rmgpy.constants as constants
import rmgpy.quantity as quantity
from libc.math cimport exp, sqrt

################################################################################

def CollisionError(Exception):
    pass

################################################################################

cdef class LennardJones:
    """
    A set of Lennard-Jones collision parameters. The attributes are:

    =================== ========================================================
    Attribute           Description
    =================== ========================================================
    `sigma`             Distance at which the inter-particle potential is minimum
    `epsilon`           Depth of the potential well
    =================== ========================================================
    
    """

    def __init__(self, sigma=None, epsilon=None):
        self.sigma = sigma
        self.epsilon = epsilon

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        object.
        """
        return 'LennardJones(sigma={0!r}, epsilon={1!r})'.format(self.sigma, self.epsilon)

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (LennardJones, (self.sigma, self.epsilon))
    
    property sigma:
        """The distance at which the inter-particle potential is minimum."""
        def __get__(self):
            return self._sigma
        def __set__(self, value):
            self._sigma = quantity.Length(value)

    property epsilon:
        """The depth of the potential well."""
        def __get__(self):
            return self._epsilon
        def __set__(self, value):
            try:
                self._epsilon = quantity.Temperature(value)
                self._epsilon.value_si *= constants.R
                self._epsilon.units = 'kJ/mol'
            except quantity.QuantityError:
                self._epsilon = quantity.Energy(value)

    cpdef double getCollisionFrequency(self, double T, double M, double mu) except -1:
        """
        Return the value of the Lennard-Jones collision frequency in Hz at the
        given temperature `T` in K for colliders with the given concentration
        `M` in mol/m^3 and reduced mass `mu` in amu.
        """
        cdef double sigma, epsilon
        cdef double Tred, omega22
        sigma = self._sigma.value_si
        epsilon = self._epsilon.value_si
        M *= constants.Na       # mol/m^3 -> molecules/m^3
        Tred = constants.R * T / epsilon
        omega22 = 1.16145 * Tred**(-0.14874) + 0.52487 * exp(-0.77320 * Tred) + 2.16178 * exp(-2.43787 * Tred)
        mu *= constants.amu
        return omega22 * sqrt(8 * constants.kB * T / constants.pi / mu) * constants.pi * sigma * sigma * M
