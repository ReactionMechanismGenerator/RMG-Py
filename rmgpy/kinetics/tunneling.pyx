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
