#!/usr/bin/env python
# encoding: utf-8

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2009-2011 by the RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

"""
This module contains classes and methods for working with statistical mechanics
models for various molecular degrees of freedom. All such models derive from 
the :class:`Mode` base class, and include:

* :class:`Translation` - 3D translational motion in an ideal gas

* :class:`RigidRotor` - 2D (linear) or 3D (nonlinear) external rotational motion, modeled as a rigid rotor

* :class:`HarmonicOscillator` - A set of independent 1D vibrational motions, modeled as harmonic oscillators

* :class:`HinderedRotor` - 1D internal (hindered) torsional rotation using a simple cosine or Fourier series potential

A list of molecular degrees of freedom can be stored in a :class:`StatesModel`
object.
"""

################################################################################

import rmgpy.constants as constants
cimport rmgpy.constants as constants

import cython
import numpy
cimport numpy

@cython.boundscheck(False)
def convolve(numpy.ndarray[numpy.float64_t,ndim=1] rho1,
    numpy.ndarray[numpy.float64_t,ndim=1] rho2,
    numpy.ndarray[numpy.float64_t,ndim=1] Elist):
    """
    Convolutes two density of states arrays `rho1` and `rho2` with corresponding
    energies `Elist` together using the equation

    .. math:: \\rho(E) = \\int_0^E \\rho_1(x) \\rho_2(E-x) \\, dx

    The units of the parameters do not matter so long as they are consistent.
    """

    cdef numpy.ndarray[numpy.float64_t,ndim=1] rho
    cdef bint found1, found2
    cdef double dE
    cdef int nE, i, j
    
    rho = numpy.zeros_like(Elist)

    found1 = rho1.any(); found2 = rho2.any()
    if not found1 and not found2:
        pass
    elif found1 and not found2:
        rho = rho1
    elif not found1 and found2:
        rho = rho2
    else:
        dE = Elist[1] - Elist[0]
        nE = Elist.shape[0]
        for i in range(nE):
            for j in range(i+1):
                rho[i] += rho2[i-j] * rho1[j] * dE

    return rho

@cython.boundscheck(False)
def convolveBS(numpy.ndarray[numpy.float64_t,ndim=1] frequencies,
    numpy.ndarray[numpy.float64_t,ndim=1] Elist,
    numpy.ndarray[numpy.float64_t,ndim=1] rho0=None):
    
    cdef numpy.ndarray[numpy.float64_t,ndim=1] rho
    cdef double freq, dE
    cdef int nE, n, dn, s, Nfreq
    
    if rho0 is not None:
        rho = rho0
    else:
        rho = numpy.zeros_like(Elist)
        rho[0] = 1.0
    
    dE = Elist[1] - Elist[0]
    nE = Elist.shape[0]
    Nfreq = frequencies.shape[0]
    for s in range(Nfreq):
        freq = frequencies[s]
        dn = int(freq * constants.h * constants.c * 100 * constants.Na / dE)
        for n in range(dn+1, nE):
            rho[n] = rho[n] + rho[n-dn]
    return rho
