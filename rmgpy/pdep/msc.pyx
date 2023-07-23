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
Contains functionality for computing pressure-dependent phenomenological
rate coefficients :math:`k(T,P)` using the modified strong collision method.
"""

import logging

import numpy as np
cimport numpy as np
from libc.math cimport exp

import rmgpy.constants as constants
from rmgpy.exceptions import ModifiedStrongCollisionError

################################################################################

cpdef apply_modified_strong_collision_method(network, str efficiency_model='default'):
    """A method for applying the Modified Strong Collision approach for solving the master equation."""
    cdef np.ndarray[np.int_t,ndim=1] j_list
    cdef np.ndarray[np.float64_t,ndim=1] e_list, coll_freq, coll_eff, d_e_down, E0, e_reac
    cdef np.ndarray[np.float64_t,ndim=2] a_mat, b, k, x
    cdef np.ndarray[np.float64_t,ndim=3] dens_states
    cdef np.ndarray[np.float64_t,ndim=4] k_ij, g_nj, f_im, pa
    cdef double temperature, val, beta
    cdef int n_isom, n_reac, n_prod, n_grains, n_j
    cdef int i, j, n, r, s, start, src

    temperature = network.T
    e_list = network.e_list
    j_list = network.j_list
    dens_states = network.dens_states
    coll_freq = network.coll_freq
    k_ij = network.Kij
    f_im = network.Fim
    g_nj = network.Gnj
    E0 = network.E0
    n_isom = network.n_isom
    n_reac = network.n_reac
    n_prod = network.n_prod
    n_grains = network.n_grains
    n_j = network.n_j
    
    if np.isnan(dens_states.sum()):
        raise AttributeError('Network {0} has NaN in the density of states. '
                             'This will prevent adequate solution to the network'.format(network.label))

    k = np.zeros((n_isom + n_reac + n_prod, n_isom + n_reac + n_prod), float)
    pa = np.zeros((n_isom, n_isom + n_reac, n_grains, n_j), float)

    beta = 1. / (constants.R * temperature)  # [=] mol/kJ
    
    # Determine the starting grain for the calculation
    e_reac = np.zeros(n_isom)
    start = n_grains
    for i in range(n_isom):
        for r in range(n_grains):
            if (k_ij[:, i, r, 0] > 0).any() or (g_nj[:, i, r, 0] > 0).any():
                if start > r: start = r
                e_reac[i] = e_list[r]
                break
        else:
            raise ModifiedStrongCollisionError('Unable to determine starting grain; check active-state energies.')
    
    d_e_down = np.zeros(n_isom)
    for i in range(n_isom):
        d_e_down[i] = network.isomers[i].species[0].energy_transfer_model.get_alpha(temperature)
    
    # Compute collision efficiencies
    coll_eff = np.ones(n_isom)
    if efficiency_model == 'default':
        for i in range(n_isom):
            coll_eff[i] = network.isomers[i].species[0].energy_transfer_model.calculate_collision_efficiency(
                temperature, e_list, j_list, dens_states[i, :, :], E0[i], e_reac[i])
    elif efficiency_model == 'none':
        pass
    else:
        raise ValueError('Unknown efficiency model "{0}".'.format(efficiency_model))
    
    # Zero LHS matrix and RHS vectors
    a_mat = np.zeros((n_isom, n_isom), float)
    b = np.zeros((n_isom, n_isom + n_reac), float)

    # Iterate over the grains, calculating the PSSA concentrations
    for r in range(start, n_grains):
        for s in range(n_j):
            # Populate LHS matrix
            # Collisional deactivation
            for i in range(n_isom):
                a_mat[i, i] = -coll_freq[i] * coll_eff[i]
            # Isomerization reactions
            for i in range(n_isom):
                for j in range(i):
                    a_mat[i, j] = k_ij[i, j, r, s]
                    a_mat[j, j] -= k_ij[i, j, r, s]
                    a_mat[j, i] = k_ij[j, i, r, s]
                    a_mat[i, i] -= k_ij[j, i, r, s]
            # Dissociation reactions
            for n in range(n_reac + n_prod):
                for j in range(n_isom):
                    a_mat[j, j] -= g_nj[n, j, r, s]
    
            # Populate RHS vectors, one per isomer and reactant
            for i in range(n_isom):
                # Thermal activation via collisions
                b[i, i] = coll_freq[i] * coll_eff[i] * dens_states[i, r, s] \
                          * (2 * j_list[s] + 1) * exp(-e_list[r] * beta)
                if np.isnan(b[i, i]):
                    logging.warning('Non-number generated for grain {0} for isomer {1}'.format(r,network.isomers[i]))
            for n in range(n_isom, n_isom + n_reac):
                # Chemical activation via association reaction
                for j in range(n_isom):
                    b[j, n] = f_im[j, n - n_isom, r, s] * dens_states[n, r, s] * (2 * j_list[s] + 1) \
                              * exp(-e_list[r] * beta)
                    if np.isnan(b[j, n]):
                        logging.warning('Non-number generated for grain {0} for isomer {1} and isomer/reactant '
                                        '{2}'.format(r, network.isomers[j], network.reactants[ - n_isom]))
                        logging.debug(str([f_im[j, n - n_isom, r, s], dens_states[n, r, s],
                                           (2 * j_list[s] + 1), exp(-e_list[r] * beta), e_list[r], beta]))
            # Solve for steady-state population
            x = -np.linalg.solve(a_mat, b)
            for n in range(n_isom + n_reac):
                for i in range(n_isom):
                    pa[i, n, r, s] = x[i, n]
                        
    # Check that our populations are all positive
    if not (pa >= 0).all():
        for reactant_index in range(len(pa[:, 0, 0, 0])):
            for isomer_index in range(len(pa[0, :, 0, 0])):
                populations =pa[reactant_index,isomer_index, :, :]
                if not (populations >= 0).all():
                    logging.debug('A negative concentration was encountered for reactant/isomer {0} and isomer {1} '
                                  'with matrix\n{2}'.format(network.isomers + network.reactants,
                                                            network.reactants,populations))
        raise ModifiedStrongCollisionError('A negative steady-state concentration was encountered.')

    # Compute rate coefficients from PSSA concentrations
    for src in range(n_isom + n_reac):
        # Calculate stabilization rates (i.e.) R + R' --> Ai or M --> Ai
        for i in range(n_isom):
            if i != src:
                val = coll_freq[i] * coll_eff[i] * np.sum(pa[i,src, :, :])
                k[i, src] += val
                k[src, src] -= val
        # Calculate dissociation rates (i.e.) R + R' --> Bn + Cn or M --> Bn + Cn
        for n in range(n_reac+n_prod):
            for j in range(n_isom):
                if n + n_isom != src:
                    val = np.sum(g_nj[n, j, :, :] * pa[j, src, :, :])
                    k[n + n_isom, src] += val
                    k[src, src] -= val
    # To complete pa we need the Boltzmann distribution at low energies
    for i in range(n_isom):
        for r in range(n_grains):
            for s in range(n_j):
                if pa[i, i, r, s] == 0:
                    pa[i,i, r, s] = dens_states[i, r, s] * (2 * j_list[s] + 1) * exp(-e_list[r] * beta)
                
    # Return the matrix of k(T,P) values and the pseudo-steady population distributions
    return k, pa
