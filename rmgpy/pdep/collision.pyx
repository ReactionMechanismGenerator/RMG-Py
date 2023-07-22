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
This module contains classes and functions for working with collision models.
"""

import logging

cimport cython
import numpy as np
cimport numpy as np
from libc.math cimport exp

cimport rmgpy.constants as constants
import rmgpy.quantity as quantity
from rmgpy.exceptions import CollisionError

################################################################################


cdef class SingleExponentialDown(RMGObject):
    """
    A representation of a single exponential down model of collisional energy
    transfer. The attributes are:
    
    =================== ========================================================
    Attribute           Description
    =================== ========================================================
    `alpha0`            The average energy transferred in a deactivating collision at the reference temperature
    `T0`                The reference temperature
    `n`                 The temperature exponent
    =================== ========================================================
    
    Based around the collisional energy transfer probability function
    
    .. math:: P(E, E^\\prime) = C(E^\\prime) \\exp \\left( -\\frac{E^\\prime - E}{\\alpha} \\right) \\hspace{40pt} E < E^\\prime
    
    where the parameter :math:`\\alpha = \\left< \\Delta E_\\mathrm{d} \\right>`
    represents the average energy transferred in a deactivating collision. This
    is the most commonly-used collision model, simply because it only has one
    parameter to determine. The parameter :math:`\\alpha` is specified using the
    equation
    
    .. math:: \\alpha = \\alpha_0 \\left( \\frac{T}{T_0} \\right)^n
    
    where :math:`\\alpha_0` is the value of :math:`\\alpha` at temperature
    :math:`T_0` in K. Set the exponent :math:`n` to zero to obtain a
    temperature-independent value for :math:`\\alpha`.
    """

    def __init__(self, alpha0=None, T0=None, n=0.0):
        self.alpha0 = alpha0
        self.T0 = T0
        self.n = n

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        object.
        """
        return 'SingleExponentialDown(alpha0={0!r}, T0={1!r}, n={2:g})'.format(self.alpha0, self.T0, self.n)

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return SingleExponentialDown, (self.alpha0, self.T0, self.n)

    property alpha0:
        """The average energy transferred in a deactivating collision at the reference temperature."""
        def __get__(self):
            return self._alpha0
        def __set__(self, value):
            if value is not None:
                try:
                    self._alpha0 = quantity.Frequency(value)
                    self._alpha0.value_si *= constants.h * constants.c * 100. * constants.Na
                    self._alpha0.units = 'kJ/mol'
                except quantity.QuantityError:
                    self._alpha0 = quantity.Energy(value)

    property T0:
        """The reference temperature."""
        def __get__(self):
            return self._t0
        def __set__(self, value):
            self._t0 = quantity.Temperature(value)

    cpdef double get_alpha(self, double T) except -1000000000:
        """
        Return the value of the :math:`\\alpha` parameter - the average energy
        transferred in a deactivating collision - in J/mol at temperature `T`
        in K.
        """
        cdef double alpha0, t0
        alpha0 = self._alpha0.value_si
        if self._t0 is None:
            return alpha0
        else:
            t0 = self._t0.value_si
            return alpha0 * (T / t0) ** self.n

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def generate_collision_matrix(self, double T,
                                  np.ndarray[float_t,ndim=2] dens_states,
                                  np.ndarray[float_t,ndim=1] e_list,
                                  np.ndarray[np.int_t,ndim=1] j_list=None):
        """
        Generate and return the collision matrix
        :math:`\\matrix{M}_\\mathrm{coll} / \\omega = \\matrix{P} - \\matrix{I}`
        corresponding to this collision model for a given set of energies
        `e_list` in J/mol, temperature `T` in K, and isomer density of states
        `dens_states`.
        """

        cdef double alpha, beta
        cdef double c, left, right
        cdef int n_grains, start, i, r, s, u, v
        cdef np.ndarray[float_t,ndim=1] rho
        cdef np.ndarray[float_t,ndim=2] phi, p0
        cdef np.ndarray[float_t,ndim=4] p

        n_grains = e_list.shape[0]
        n_j = j_list.shape[0] if j_list is not None else 1
        p = np.zeros((n_grains, n_j, n_grains, n_j), float)
        p0 = np.zeros((n_grains, n_grains), float)

        alpha = 1.0 / self.get_alpha(T)
        beta = 1.0 / (constants.R * T)
        
        if n_j > 1:
            rho = np.zeros(n_grains)
            for r in range(n_grains):
                rho[r] = np.sum((2 * j_list + 1) * dens_states[r, :])
        else:
            rho = dens_states[:, 0]
        
        for start in range(n_grains):
            if rho[start] > 0:
                break

        # Determine unnormalized entries in collisional transfer probability matrix
        for r in range(start, n_grains):
            for s in range(start, r + 1):
                p0[s, r] = exp(-(e_list[r] - e_list[s]) * alpha)
            for s in range(r+1,n_grains):
                p0[s, r] = exp(-(e_list[s] - e_list[r]) * alpha) * rho[s] / rho[r] * exp(-(e_list[s] - e_list[r]) * beta)
        
        # Normalize using detailed balance
        # This method is much more robust, and corresponds to:
        #    [ 1 1 1 1 ...]
        #    [ 1 2 2 2 ...]
        #    [ 1 2 3 3 ...]
        #    [ 1 2 3 4 ...]
        for r in range(start, n_grains):
            left, right = 0.0, 0.0
            for s in range(start, r):
                left += p0[s, r]
            for s in range(r, n_grains):
                right += p0[s, r]
            c = (1 - left) / right
            # Check for normalization consistency (i.e. all numbers are positive)
            if c < 0:
                raise CollisionError('Encountered negative normalization coefficient while normalizing '
                                     'collisional transfer probabilities matrix.')
            for s in range(r + 1, n_grains):
                p0[r, s] *= c
                p0[s, r] *= c
            p0[r, r] = p0[r, r] * c - 1
        # This method is described by Pilling and Holbrook, and corresponds to:
        #    [ ... 4 3 2 1 ]
        #    [ ... 3 3 2 1 ]
        #    [ ... 2 2 2 1 ]
        #    [ ... 1 1 1 1 ]
        # for r in range(n_grains, start, -1):
        #     left = 0.0; right = 0.0
        #     for s in range(start, r): left += p0[s,r]
        #     for s in range(r, n_grains): right += p0[s,r]
        #     c = (1 - right) / left
        #     # Check for normalization consistency (i.e. all numbers are positive)
        #     if c < 0: raise CollisionError('Encountered negative normalization coefficient while normalizing '
        #                                    'collisional transfer probabilities matrix.')
        #     for s in range(r-1):
        #         p0[r,s] *= c
        #         p0[s,r] *= c
        #     p0[r,r] = p0[r,r] * c - 1

        # If solving the 2D master equation, compute P(E,J,E',J') from P(E,E')
        # by assuming that the J distribution after the collision is independent
        # of that before the collision (the strong collision approximation in J)
        if n_j > 1:
            phi = np.zeros_like(dens_states)
            for s in range(n_j):
                phi[:,s] = (2 * j_list[s] + 1) * dens_states[:, s]
            for r in range(start, n_grains):
                phi[r,:] /= rho[r]
            for r in range(start, n_grains):
                for s in range(n_j):
                    for u in range(start, n_grains):
                        for v in range(n_j):
                            p[r, s, u, v] = p0[r, u] * phi[r, s]
        else:
            p[:,0,:,0] = p0
            
        return p

    def calculate_collision_efficiency(self,
                                       double T,
                                       np.ndarray[float_t,ndim=1] e_list,
                                       np.ndarray[np.int_t,ndim=1] j_list,
                                       np.ndarray[float_t,ndim=2] dens_states,
                                       double E0, double e_reac):
        """
        Calculate an efficiency factor for collisions, particularly useful for the
        modified strong collision method. The collisions involve the given 
        `species` with density of states `dens_states` corresponding to energies
        e_list` in J/mol, ground-state energy `E0` in kJ/mol, and first
        reactive energy `e_reac` in kJ/mol. The collisions occur at temperature `T`
        in K and are described by the average energy transferred in a deactivating
        collision `d_e_down` in kJ/mol. The algorithm here is implemented as
        described by A. Y. Chang, J. W. Bozzelli, and A. M. Dean.
        *Z. Phys. Chem.* **214**, p. 1533-1568 (2000).
        `doi: 10.1524/zpch.2000.214.11.1533 <https://doi.org/10.1524/zpch.2000.214.11.1533>`_
    
        """
    
        cdef double d_e_down, d_e, fe, fe_num, fe_den, delta1, delta2, delta_n, delta, value, beta
        cdef double gas_constant = constants.R
        cdef int n_grains, n_j, r
    
        # Ensure that the barrier height is sufficiently above the ground state
        # Otherwise invalid efficiencies are observed
        if e_reac - E0 < 100:
            e_reac = E0 + 100
    
        d_e_down = self.get_alpha(T)
    
        n_grains = len(e_list)
        n_j = 1 if j_list is None else len(j_list)
        d_e = e_list[1] - e_list[0]
        
        fe_num, fe_den, delta1, delta2, delta_n, delta = 0, 0, 0, 0, 0, 1
    
        for r in range(n_grains):
            value = 0.0
            for s in range(n_j):
                value += dens_states[r, s] * (2 * j_list[s] + 1) * exp(-e_list[r] / (gas_constant * T))
            if e_list[r] > e_reac:
                fe_num += value
                if fe_den == 0:
                    fe_den = value * gas_constant * T / d_e
        if fe_den == 0:
            return 1.0
        fe = fe_num / fe_den
    
        # Chang, Bozzelli, and Dean recommend "freezing out" fe at values greater
        # than 1e6 to avoid issues of roundoff error
        # They claim that the collision efficiency isn't too temperature-dependent
        # in this regime, so it's an okay approximation to use
        if fe > 1e6:
            fe = 1e6
        
        for r in range(n_grains):
            value = 0.0
            for s in range(n_j):
                value += dens_states[r, s] * (2 * j_list[s] + 1) * exp(-e_list[r] / (gas_constant * T))
            # Delta
            if e_list[r] < e_reac:
                delta1 += value
                delta2 += value * exp(-(e_reac - e_list[r]) / (fe * gas_constant * T))
            delta_n += value
    
        delta1 /= delta_n
        delta2 /= delta_n
    
        delta = delta1 - (fe * gas_constant * T) / (d_e_down + fe * gas_constant * T) * delta2
    
        beta = (d_e_down / (d_e_down + fe * gas_constant * T)) ** 2 / delta
    
        if beta > 1:
            logging.debug('Collision efficiency {0:.3f} calculated at {1:g} K is greater than unity, '
                          'so it will be set to unity.'.format(beta, T))
            beta = 1
        if beta < 0:
            raise CollisionError('Invalid collision efficiency {0:.3f} calculated at {1:g} K.'.format(beta, T))
        
        return beta
