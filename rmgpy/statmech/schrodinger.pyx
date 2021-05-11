# cython: embedsignature=True, cdivision=True

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2021 Prof. William H. Green (whgreen@mit.edu),           #
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

cimport cython
import numpy as np
cimport numpy as np
from libc.math cimport exp, log

cimport rmgpy.constants as constants

################################################################################

def unit_degeneracy(n):
    return 1

cpdef double get_partition_function(double T, energy, degeneracy=unit_degeneracy, int n0=0, int nmax=10000,
                                    double tol=1e-12) except -1:
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
    cdef double Q, dQ, E_n
    cdef int g_n
    cdef double beta = 1. / (constants.R * T)  # [=] mol/J

    Q = 0.0
    for n in range(n0, nmax):
        E_n = energy(n)
        g_n = degeneracy(n)

        dQ = g_n * exp(-beta * E_n)
        Q += dQ

        if dQ < tol * Q:
            break

    return Q

cpdef double get_heat_capacity(double T, energy, degeneracy=unit_degeneracy, int n0=0, int nmax=10000,
                               double tol=1e-12) except -100000000:
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
    cdef double Q, dQ, sumE, dsumE, sumE2, dsumE2, E_n
    cdef int g_n
    cdef double beta = 1. / (constants.R * T)  # [=] mol/J

    Q = 0.0
    sumE = 0.0
    sumE2 = 0.0
    for n in range(n0, nmax):
        E_n = energy(n)
        g_n = degeneracy(n)

        dQ = g_n * exp(-beta * E_n)
        dsumE = E_n * dQ
        dsumE2 = E_n * dsumE
        Q += dQ
        sumE += dsumE
        sumE2 += dsumE2

        if dQ < tol * Q and dsumE < tol * sumE and dsumE2 < tol * sumE2:
            break

    return beta * beta * (sumE2 / Q - sumE * sumE / (Q * Q))

cpdef double get_enthalpy(double T, energy, degeneracy=unit_degeneracy, int n0=0, int nmax=10000,
                          double tol=1e-12) except 100000000:
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
    cdef double Q, dQ, sumE, dsumE, E_n
    cdef int g_n
    cdef double beta = 1. / (constants.R * T)  # [=] mol/J

    Q = 0.0
    sumE = 0.0
    for n in range(n0, nmax):
        E_n = energy(n)
        g_n = degeneracy(n)

        dQ = g_n * exp(-beta * E_n)
        dsumE = E_n * dQ
        Q += dQ
        sumE += dsumE

        if dQ < tol * Q and dsumE < tol * sumE:
            break

    return beta * sumE / Q

cpdef double get_entropy(double T, energy, degeneracy=unit_degeneracy, int n0=0, int nmax=10000,
                         double tol=1e-12) except -100000000:
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
    cdef double Q, dQ, sumE, dsumE, E_n
    cdef int g_n
    cdef double beta = 1. / (constants.R * T)  # [=] mol/J

    Q = 0.0
    sumE = 0.0
    for n in range(n0, nmax):
        E_n = energy(n)
        g_n = degeneracy(n)

        dQ = g_n * exp(-beta * E_n)
        dsumE = E_n * dQ
        Q += dQ
        sumE += dsumE

        if dQ < tol * Q and dsumE < tol * sumE:
            break

    return log(Q) + beta * sumE / Q

################################################################################

cpdef np.ndarray get_sum_of_states(np.ndarray e_list, energy, degeneracy=unit_degeneracy, int n0=0,
                                   np.ndarray sum_states_0=None):
    """
    Return the values of the sum of states :math:`N(E)` for a given set of
    energies `e_list` in J/mol above the ground state using an initial sum of
    states `sum_states_0`. The solution to the Schrodinger equation is given
    using functions `energy` and `degeneracy` that accept as argument a quantum
    number and return the corresponding energy in J/mol and degeneracy of that
    level. The quantum number always begins at `n0` and increases by ones.
    """
    if sum_states_0 is None:
        sum_states_0 = np.ones_like(e_list)
    return convolve_bssr(e_list, sum_states_0, energy, degeneracy, n0)

cpdef np.ndarray get_density_of_states(np.ndarray e_list, energy, degeneracy=unit_degeneracy, int n0=0,
                                       np.ndarray dens_states_0=None):
    """
    Return the values of the dimensionless density of states
    :math:`\\rho(E) \\ dE` for a given set of energies `e_list` in J/mol above
    the ground state using an initial density of states `dens_states_0`. The
    solution to the Schrodinger equation is given using functions `energy` and
    `degeneracy` that accept as argument a quantum number and return the
    corresponding energy in J/mol and degeneracy of that level. The quantum 
    number always begins at `n0` and increases by ones.
    """
    if dens_states_0 is None:
        dens_states_0 = np.zeros_like(e_list)
        dens_states_0[0] = 1.0
    return convolve_bssr(e_list, dens_states_0, energy, degeneracy, n0)

################################################################################

@cython.boundscheck(False)
@cython.wraparound(False)
def convolve(np.ndarray[np.float64_t, ndim=1] rho1, np.ndarray[np.float64_t, ndim=1] rho2):
    """
    Return the convolution of two arrays `rho1` and `rho2`.
    """
    cdef np.ndarray[np.float64_t, ndim=1] rho
    cdef int i, j, nE

    if rho1.shape[0] != rho2.shape[0]:
        raise ValueError('Attempted to convolve an array of length {0:d} with an array of '
                         'length {1:d}.'.format(len(rho1), len(rho2)))

    nE = rho1.shape[0]
    rho = np.zeros_like(rho1)

    for i in range(nE):
        for j in range(i + 1):
            rho[i] += rho2[i - j] * rho1[j]

    return rho

@cython.boundscheck(False)
@cython.wraparound(False)
def convolve_bs(np.ndarray[np.float64_t, ndim=1] e_list,
                np.ndarray[np.float64_t, ndim=1] rho0,
                double energy, int degeneracy=1):
    """
    Convolve a molecular degree of freedom into a density or sum of states
    using the Beyer-Swinehart (BS) direct count algorithm. This algorithm is
    suitable for unevenly-spaced energy levels in the array of energy grains
    `e_list` (in J/mol), but assumes the solution of the Schrodinger equation
    gives evenly-spaced energy levels with spacing `energy` in kJ/mol and
    degeneracy `degeneracy`.
    """
    cdef int n = 0
    cdef double Emax = np.max(e_list), E_n
    cdef int i, j, nE = e_list.shape[0], g_n
    cdef np.ndarray[np.float64_t, ndim=1] rho

    rho = rho0.copy()

    for i in range(nE):
        for j in range(i, nE):
            if e_list[j] - e_list[i] >= 0.9999 * energy:
                rho[j] += degeneracy * rho[i]
                break

    return rho

@cython.boundscheck(False)
@cython.wraparound(False)
def convolve_bssr(np.ndarray[np.float64_t, ndim=1] e_list,
                  np.ndarray[np.float64_t, ndim=1] rho0,
                  energy, degeneracy=unit_degeneracy, int n0=0):
    """
    Convolve a molecular degree of freedom into a density or sum of states
    using the Beyer-Swinehart-Stein-Rabinovitch (BSSR) direct count algorithm.
    This algorithm is suitable for unevenly-spaced energy levels in both the
    array of energy grains `e_list` (in J/mol) and the energy levels
    corresponding to the solution of the Schrodinger equation.
    """
    cdef int n = n0
    cdef double Emax = np.max(e_list), E_n
    cdef int i, j, nE = e_list.shape[0], g_n
    cdef np.ndarray[np.float64_t, ndim=1] rho

    rho = np.zeros_like(rho0)

    E_n = energy(n)
    g_n = degeneracy(n)
    while E_n < Emax:

        for i in range(nE):
            for j in range(i, nE):
                if e_list[j] - e_list[i] >= 0.9999 * E_n:
                    rho[j] += g_n * rho0[i]
                    break

        n += 1
        E_n = energy(n)
        g_n = degeneracy(n)

    return rho
