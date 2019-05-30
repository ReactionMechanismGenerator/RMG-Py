# cython: embedsignature=True, cdivision=True

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2019 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

"""
The :mod:`rmgpy.statmech.schrodinger` module contains functionality for
working with the Schrodinger equation and its solution. In particular, it
contains functions for using the energy levels and corresponding degeneracies
obtained from solving the Schrodinger equation to compute various thermodynamic
and statistical mechanical properties, such as heat capacity, enthalpy,
entropy, partition function, and the sum and density of states.
"""

import math
import numpy
cimport cython
from libc.math cimport exp, log

cimport rmgpy.constants as constants

################################################################################


def unitDegeneracy(int n):
    return 1


def unit_qk_energy(int n):
    return 1


cpdef double getPartitionFunction(double T, energy, degeneracy=unitDegeneracy, int n0=0,
                                  qk_energy=unit_qk_energy, int nmax=10000, double tol=1e-12) except -1:
    """
    Return the value of the partition function :math:`Q(T)` at a given
    temperature `T` in K. The solution to the Schrodinger equation is given
    using functions `energy` and `degeneracy` that accept as argument a quantum
    number and return the corresponding energy in J/mol and degeneracy of that
    level. The quantum number always begins at `n0` and increases by ones. You
    can also change the relative tolerance `tol` and the maximum allowed value
    of the quantum number `nmax`.
    """
    cdef int n
    cdef double q, dq, e_n, qk
    cdef int g_n
    cdef double beta = 1. / (constants.R * T)  # in mol/J

    q = 0.0
    for n in range(n0, nmax):
        e_n, qk, g_n = energy(n), qk_energy(n), degeneracy(n)
            
        dq = g_n * exp(-beta * e_n) * qk
        q += dq

        if dq < tol * q:
            break

    return q


cpdef double getHeatCapacity(double T, energy, degeneracy=unitDegeneracy, int n0=0,
                             qk_energy=unit_qk_energy, int nmax=10000, double tol=1e-12) except -100000000:
    """
    Return the value of the dimensionless heat capacity :math:`C_\\mathrm{v}(T)/R`
    at a given temperature `T` in K. The solution to the Schrodinger equation
    is given using functions `energy` and `degeneracy` that accept as argument
    a quantum number and  return the corresponding energy in J/mol and
    degeneracy of that level. The quantum number always begins at `n0` and
    increases by ones. You can also change the relative tolerance `tol` and the
    maximum allowed value of the quantum number `nmax`.
    """    
    cdef int n
    cdef double q, dq, sum_e, dsum_e, sum_e2, dsum_e2, e_n, qk
    cdef int g_n
    cdef double beta = 1. / (constants.R * T)  # in mol/J
    
    q, sum_e, sum_e2 = 0.0, 0.0, 0.0
    for n in range(n0, nmax):
        e_n, qk, g_n = energy(n), qk_energy(n), degeneracy(n)
            
        dq = g_n * exp(-beta * e_n) * qk
        dsum_e = e_n * qk ** 2 * dq
        dsum_e2 = e_n * qk ** 2 * dsum_e
        q += dq
        sum_e += dsum_e
        sum_e2 += dsum_e2

        if dq < tol * q and dsum_e < tol * sum_e and dsum_e2 < tol * sum_e2:
            break

    return beta ** 2 * (sum_e2 / q - sum_e ** 2 / (q ** 2))


cpdef double getEnthalpy(double T, energy, degeneracy=unitDegeneracy, int n0=0,
                         qk_energy=unit_qk_energy, int nmax=10000, double tol=1e-12) except 100000000:
    """
    Return the value of the dimensionless enthalpy :math:`H(T)/RT` at a given
    temperature `T` in K. The solution to the Schrodinger equation is given
    using functions `energy` and `degeneracy` that accept as argument a quantum
    number and return the corresponding energy in J/mol and degeneracy of that
    level. The quantum number always begins at `n0` and increases by ones. You
    can also change the relative tolerance `tol` and the maximum allowed value
    of the quantum number `nmax`.
    """
    cdef int n
    cdef double q, dq, sum_e, dsum_e, e_n
    cdef int g_n
    cdef double beta = 1. / (constants.R * T)  # in mol/J
    
    q, sum_e = 0.0, 0.0
    for n in range(n0, nmax):
        e_n, qk, g_n = energy(n), qk_energy(n), degeneracy(n)
            
        dq = g_n * exp(-beta * e_n) * qk
        dsum_e = e_n * qk ** 2 * dq
        q += dq
        sum_e += dsum_e

        if dq < tol * q and dsum_e < tol * sum_e:
            break

    return beta * sum_e / q


cpdef double getEntropy(double T, energy, degeneracy=unitDegeneracy, int n0=0,
                        qk_energy=unit_qk_energy, int nmax=10000, double tol=1e-12) except -100000000:
    """
    Return the value of the dimensionless entropy :math:`S(T)/R` at a given
    temperature `T` in K. The solution to the Schrodinger equation is given
    using functions `energy` and `degeneracy` that accept as argument a quantum
    number and return the corresponding energy in J/mol and degeneracy of that
    level. The quantum number always begins at `n0` and increases by ones. You
    can also change the relative tolerance `tol` and the maximum allowed value
    of the quantum number `nmax`.
    """
    cdef int n
    cdef double q, dq, sum_e, dsum_e, e_n
    cdef int g_n
    cdef double beta = 1. / (constants.R * T)  # in mol/J
    
    q, sum_e = 0.0, 0.0
    for n in range(n0, nmax):
        e_n, qk, g_n = energy(n), qk_energy(n), degeneracy(n)
        
        dq = g_n * exp(-beta * e_n) * qk
        dsum_e = e_n * qk ** 2 * dq
        q += dq
        sum_e += dsum_e
        
        if dq < tol * q and dsum_e < tol * sum_e:
            break

    return log(q) + beta * sum_e / q

################################################################################


cpdef numpy.ndarray getSumOfStates(numpy.ndarray Elist, energy, degeneracy=unitDegeneracy, int n0=0,
                                   numpy.ndarray sumStates0=None, qk_energy=unit_qk_energy):
    """
    Return the values of the sum of states :math:`N(E)` for a given set of
    energies `Elist` in J/mol above the ground state using an initial sum of
    states `sumStates0`. The solution to the Schrodinger equation is given
    using functions `energy` and `degeneracy` that accept as argument a quantum
    number and return the corresponding energy in J/mol and degeneracy of that
    level. The quantum number always begins at `n0` and increases by ones.
    """
    if sumStates0 is None:
        sumStates0 = numpy.ones_like(Elist)
    return convolveBSSR(Elist, sumStates0, energy, degeneracy, n0, qk_energy)


cpdef numpy.ndarray getDensityOfStates(numpy.ndarray Elist, energy, degeneracy=unitDegeneracy, int n0=0,
                                       numpy.ndarray densStates0=None, qk_energy=unit_qk_energy):
    """
    Return the values of the dimensionless density of states
    :math:`\\rho(E) \\ dE` for a given set of energies `Elist` in J/mol above
    the ground state using an initial density of states `densStates0`. The
    solution to the Schrodinger equation is given using functions `energy` and
    `degeneracy` that accept as argument a quantum number and return the
    corresponding energy in J/mol and degeneracy of that level. The quantum 
    number always begins at `n0` and increases by ones.
    """
    if densStates0 is None:
        densStates0 = numpy.zeros_like(Elist)
        densStates0[0] = 1.0
    return convolveBSSR(Elist, densStates0, energy, degeneracy, n0, qk_energy)

################################################################################


@cython.boundscheck(False)
@cython.wraparound(False)
def convolve(numpy.ndarray[numpy.float64_t, ndim=1] rho1, numpy.ndarray[numpy.float64_t, ndim=1] rho2):
    """
    Return the convolution of two arrays `rho1` and `rho2`.
    """
    cdef numpy.ndarray[numpy.float64_t, ndim=1] rho
    cdef int i, j, n_e
    
    if rho1.shape[0] != rho2.shape[0]:
        raise ValueError('Attempted to convolve an array of length {0:d} with an array of length {1:d}.'.format(
            len(rho1), len(rho2)))
    
    n_e = rho1.shape[0]
    rho = numpy.zeros_like(rho1)
    
    for i in range(n_e):
        for j in range(i+1):
            rho[i] += rho2[i-j] * rho1[j]

    return rho


@cython.boundscheck(False)
@cython.wraparound(False)
def convolveBS(numpy.ndarray[numpy.float64_t, ndim=1] Elist, numpy.ndarray[numpy.float64_t, ndim=1] rho0,
               energy, degeneracy=1):
    """
    Convolve a molecular degree of freedom into a density or sum of states
    using the Beyer-Swinehart (BS) direct count algorithm. This algorithm is
    suitable for unevenly-spaced energy levels in the array of energy grains
    `Elist` (in J/mol), but assumes the solution of the Schrodinger equation
    gives evenly-spaced energy levels with spacing `energy` in kJ/mol and
    degeneracy `degeneracy`.
    """
    cdef int n = 0
    cdef int i, j, n_e = Elist.shape[0], g_n
    cdef numpy.ndarray[numpy.float64_t, ndim=1] rho
    
    rho = rho0.copy()
    
    for i in range(n_e):
        for j in range(i, n_e):
            if Elist[j] - Elist[i] >= 0.9999 * energy:
                rho[j] += degeneracy * rho[i]
                break
        
    return rho


@cython.boundscheck(False)
@cython.wraparound(False)
def convolveBSSR(numpy.ndarray[numpy.float64_t, ndim=1] Elist, numpy.ndarray[numpy.float64_t, ndim=1] rho0,
                 energy, degeneracy=unitDegeneracy, int n0=0, qk_energy=unit_qk_energy):
    """
    Convolve a molecular degree of freedom into a density or sum of states
    using the Beyer-Swinehart-Stein-Rabinovitch (BSSR) direct count algorithm.
    This algorithm is suitable for unevenly-spaced energy levels in both the
    array of energy grains `Elist` (in J/mol) and the energy levels
    corresponding to the solution of the Schrodinger equation.
    """
    cdef int n = n0
    cdef double e_max = numpy.max(Elist), e_n
    cdef int i, j, n_e = Elist.shape[0], g_n
    cdef numpy.ndarray[numpy.float64_t, ndim=1] rho
    
    rho = numpy.zeros_like(rho0)
    
    e_n, qk, g_n = energy(n), qk_energy(n), degeneracy(n)
    while e_n < e_max:

        for i in range(n_e):
            for j in range(i, n_e):
                if Elist[j] - Elist[i] >= 0.9999 * e_n * qk:
                    rho[j] += g_n * rho0[i]
                    break

        n += 1
        e_n, qk, g_n = energy(n), qk_energy(n), degeneracy(n)
        
    return rho
