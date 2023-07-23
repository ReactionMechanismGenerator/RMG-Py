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
rate coefficients :math:`k(T,P)` using the reservoir state method.
"""

import scipy.linalg

import numpy as np
cimport numpy as np

import rmgpy.constants as constants
from rmgpy.exceptions import ReservoirStateError

################################################################################


cpdef apply_reservoir_state_method(network):
    """A method for applying the Reservoir State approach for solving the master equation."""
    cdef np.ndarray[np.int_t,ndim=1] j_list
    cdef np.ndarray[np.int_t,ndim=2] n_res, n_act
    cdef np.ndarray[np.int_t,ndim=3] indices
    cdef np.ndarray[np.float64_t,ndim=1] e_list
    cdef np.ndarray[np.float64_t,ndim=2] active_state_mat, source_vectors, pss_active_state, k
    cdef np.ndarray[np.float64_t,ndim=3] dens_states, eq_dist
    cdef np.ndarray[np.float64_t,ndim=4] k_ij, g_nj, f_im, pa
    cdef np.ndarray[np.float64_t,ndim=5] m_coll
    cdef list ind
    cdef double temperature, tol, y, beta
    cdef int n_isom, n_reac, n_prod, n_grains, n_j, bandwidth, halfbandwidth, width, width0
    cdef int i, j, n, r, s, u, v, row, iter

    temperature = network.T
    e_list = network.e_list
    j_list = network.j_list
    dens_states = network.dens_states
    m_coll = network.Mcoll
    k_ij = network.Kij
    f_im = network.Fim
    g_nj = network.Gnj
    n_isom = network.n_isom
    n_reac = network.n_reac
    n_prod = network.n_prod
    n_grains = network.n_grains
    n_j = network.n_j

    beta = 1. / (constants.R * temperature)  # [=] mol/kJ

    k = np.zeros((n_isom + n_reac + n_prod, n_isom + n_reac + n_prod), float)
    pa = np.zeros((n_isom, n_isom + n_reac, n_grains, n_j), float)

    # Determine the reservoir cutoff grain for each isomer
    # Start by simply placing it at the lowest reactive grain
    n_res = np.zeros((n_isom, n_j), int)
    for i in range(n_isom):
        for s in range(n_j):
            for r in range(n_grains):
                if dens_states[i, r, s] != 0 and ((k_ij[:, i, r, s] > 0).any() or (g_nj[:, i, r, s] > 0).any()):
                    # We need at least one reservoir grain for the RS method to be successful
                    if r == 0 or dens_states[i, r - 1, s] == 0:
                        n_res[i, s] = r + 1
                    else:
                        n_res[i, s] = r
                    break
    n_act = n_grains - n_res

    # Determine equilibrium distributions
    eq_dist = np.zeros((n_isom + n_reac, n_grains, n_j), float)
    for i in range(n_isom + n_reac):
        for s in range(n_j):
            eq_dist[i, :, s] = dens_states[i, :, s] * (2 * j_list[s] + 1) * np.exp(-e_list * beta)

    # Determine pseudo-steady state populations of active state
    row = 0
    indices = -np.ones((n_isom, n_grains, n_j), int)
    for r in range(n_grains):
        for s in range(n_j):
            for i in range(n_isom):
                if r >= n_res[i, s]:
                    indices[i, r, s] = row
                    row += 1

    # Choose the half-bandwidth using the deepest isomer well
    width = 0
    tol = 1e-12
    for i in range(n_isom):
        for s in range(n_j):
            r = n_res[i, s]
            if m_coll[i, r, s, r, s] == 0: continue
            ratio = np.abs(m_coll[i, :, s, r, s] / m_coll[i, r, s, r, s])
            ind = [j for j,y in enumerate(ratio) if y > tol]
            if len(ind) > 0:
                width0 = max(r - min(ind), max(ind) - r)
                if width0 > width:
                    width = width0
    if width == 0:
        raise ReservoirStateError('Unable to determine half-bandwidth for active-state matrix; '
                                  'the wells may be too shallow to use the RS method.')
    halfbandwidth = (width + 1) * n_isom * n_j - n_isom
    bandwidth = 2 * halfbandwidth + 1

    # Populate active-state matrix and source vectors
    active_state_mat = np.zeros((bandwidth, np.sum(n_act)), float)
    source_vectors = np.zeros((np.sum(n_act), n_isom + n_reac), float)
    # Collisional terms
    for i in range(n_isom):
        for u in range(n_j):
            for v in range(n_j):
                for r in range(n_res[i, u], n_grains):
                    for s in range(max(n_res[i, v], r - width), min(n_grains, r + width + 1)):
                        active_state_mat[halfbandwidth + indices[i, r, u] - indices[i, s, v], indices[i, s, v]] = \
                            m_coll[i, r, u, s, v]
                    source_vectors[indices[i, r, u], i] = np.sum(m_coll[i, r, u, 0: n_res[i, u], v] *
                                                                 eq_dist[i, 0: n_res[i, u], v])

    # Isomerization terms
    for i in range(n_isom):
        for j in range(i):
            for u in range(n_j):
                for r in range(max(n_res[i, u], n_res[j, u]), n_grains):
                    active_state_mat[halfbandwidth + indices[j, r, u] - indices[i, r, u], indices[i, r, u]] = \
                        k_ij[j, i, r, u]
                    active_state_mat[halfbandwidth, indices[i, r, u]] -= k_ij[j, i, r, u]
                    active_state_mat[halfbandwidth + indices[i, r, u] - indices[j, r, u], indices[j, r, u]] = \
                        k_ij[i, j, r, u]
                    active_state_mat[halfbandwidth, indices[j, r, u]] -= k_ij[i, j, r, u]
    # Dissociation/association terms
    for i in range(n_isom):
        for n in range(n_reac + n_prod):
            for u in range(n_j):
                for r in range(n_res[i, u], n_grains):
                    active_state_mat[halfbandwidth, indices[i, r, u]] -= g_nj[n, i, r, u]
        for n in range(n_reac):
            for u in range(n_j):
                for r in range(n_res[i, u], n_grains):
                    source_vectors[indices[i, r, u], n + n_isom] = f_im[i, n, r, u] * eq_dist[n + n_isom, r, u]

    # Solve for pseudo-steady state populations of active state
    pss_active_state = scipy.linalg.solve_banded((halfbandwidth, halfbandwidth), active_state_mat, -source_vectors,
                                                 overwrite_ab=True, overwrite_b=True)
    for i in range(n_isom):
        for u in range(n_j):
            for r in range(n_res[i, u], n_grains):
                for n in range(n_isom + n_reac):
                    pa[i, n, r, u] = pss_active_state[indices[i, r, u], n]

    # Double-check to ensure that we have all positive populations
    if not (pa >= 0).all():
        raise ReservoirStateError('A negative steady-state population was encountered.')

    # Put the reservoir populations into pa as well
    for i in range(n_isom):
        for u in range(n_j):
            for r in range(n_res[i, u]):
                pa[i, i, r, u] = eq_dist[i, r, u]

    # Determine the phenomenological rate coefficients using the general procedure
    # This should be exactly the same as the procedure below, which is more
    # specific to the RS method
    # Previously it was noted that this more general approach was more robust;
    # however, more recently it seems that this is no longer the case
    # k = computeRateCoefficients(m_coll, k_ij, f_im, g_nj, pa, n_isom, n_reac, n_prod)

    # Determine the phenomenological rate coefficients
    k = np.zeros((n_isom+n_reac+n_prod, n_isom+n_reac+n_prod), float)
    # Rows relating to isomers
    for i in range(n_isom):
        for u in range(n_j):
            for v in range(n_j):
                # Collisional rearrangement within the reservoir of isomer i
                k[i, i] = k[i, i] + np.sum(np.dot(m_coll[i, 0: n_res[i, u], u, 0: n_res[i, v], v],
                                                  eq_dist[i, 0: n_res[i, v], v]))
                # Isomerization from isomer j to isomer i
                for j in range(n_isom):
                    k[i, j] = k[i, j] + np.sum(np.dot(m_coll[i, 0: n_res[i, u], u, n_res[i, v]: n_grains, v],
                                                      pa[i, j, n_res[i, v]: n_grains, v]))
                # Association from reactant n to isomer i
                for n in range(n_isom, n_isom + n_reac):
                    k[i, n] = k[i, n] + np.sum(np.dot(m_coll[i, 0: n_res[i, u], u, n_res[i, v]: n_grains, v],
                                                      pa[i, n, n_res[i, v]: n_grains, v]))
    # Rows relating to reactants
    for n in range(n_reac):
        # Association loss
        for i in range(n_isom):
            k[n_isom + n, n_isom + n] = k[n_isom + n, n_isom + n] - np.sum(f_im[i, n, :, :] * eq_dist[n + n_isom, :, :])
        # Reaction from isomer or reactant j to reactant n
        for j in range(n_isom + n_reac):
            for i in range(n_isom):
                for u in range(n_j):
                    k[n_isom + n, j] = k[n_isom + n, j] + np.sum(g_nj[n, i, n_res[i, u]: n_grains, u]
                                                                 * pa[i, j, n_res[i, u]: n_grains, u])
    # Rows relating to products
    for n in range(n_reac, n_reac + n_prod):
        # Reaction from isomer or reactant j to product n
        for j in range(n_isom + n_reac):
            for i in range(n_isom):
                for u in range(n_j):
                    k[n_isom + n, j] = k[n_isom + n, j] + np.sum(g_nj[n, i, n_res[i, u]: n_grains, u]
                                                                 * pa[i, j, n_res[i, u]: n_grains, u])

    # Ensure matrix is conservative
    for n in range(n_isom+n_reac):
        k[n, n] = k[n, n] - np.sum(k[:, n])

    # Return the matrix of k(T,P) values and the pseudo-steady population distributions
    return k, pa
