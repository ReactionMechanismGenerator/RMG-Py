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
Contains functionality for generating the master equation matrix for a given
pressure-dependent reaction network.
"""

import numpy as np
cimport numpy as np
from libc.math cimport exp

import rmgpy.constants as constants

################################################################################


cpdef generate_full_me_matrix(network, bint products=True, bint exclude_association=False, bint neglect_high_energy_collisions=False, high_energy_rate_tol=0.01):
    """
    Generate the full master equation matrix for the network.
    """

    cdef np.ndarray[np.int_t,ndim=1] j_list
    cdef np.ndarray[np.int_t,ndim=3] indices
    cdef np.ndarray[np.float64_t,ndim=1] e_list
    cdef np.ndarray[np.float64_t,ndim=2] me_mat
    cdef np.ndarray[np.float64_t,ndim=3] dens_states
    cdef np.ndarray[np.float64_t,ndim=4] k_ij, g_nj, f_im
    cdef np.ndarray[np.float64_t,ndim=5] m_coll
    cdef double temperature, pressure, beta, val, ktot
    cdef int n_isom, n_reac, n_prod, n_grains, n_j
    cdef int i, n, r, s, u, v, ind1, ind2, ind

    temperature = network.T
    # pressure = network.P  # not used in this module
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
    
    beta = 1. / (constants.R * temperature)
    
    # Construct accounting matrix
    indices = -np.ones((n_isom, n_grains, n_j), np.int)
    n_rows = 0
    for r in range(n_grains):
        for s in range(n_j):
            for i in range(n_isom):
                if dens_states[i, r, s] > 0:
                    indices[i, r, s] = n_rows
                    n_rows += 1
    n_rows += n_reac
    if products:
        n_rows += n_prod
    
    # Construct full ME matrix
    me_mat = np.zeros([n_rows, n_rows], np.float64)
    
    # Collision terms
    for i in range(n_isom):
        for r in range(n_grains):
            for s in range(n_j):
                if indices[i, r, s] > -1:
                    for u in range(r, n_grains):
                        for v in range(s, n_j):
                            me_mat[indices[i, r, s], indices[i, u, v]] = m_coll[i, r, s, u, v]
                            me_mat[indices[i, u, v], indices[i, r, s]] = m_coll[i, u, v, r, s]
    
    # Isomerization terms
    for i in range(n_isom):
        for j in range(i):
            if k_ij[i, j, n_grains - 1, 0] > 0 or k_ij[j, i, n_grains - 1, 0] > 0:
                for r in range(n_grains):
                    for s in range(n_j):
                        u, v = indices[i, r, s], indices[j, r, s]
                        if u > -1 and v > -1:
                            me_mat[v, u] = k_ij[j, i, r, s]
                            me_mat[u, u] -= k_ij[j, i, r, s]
                            me_mat[u, v] = k_ij[i, j, r, s]
                            me_mat[v, v] -= k_ij[i, j, r, s]
    
    # Association/dissociation terms
    for i in range(n_isom):
        for n in range(n_reac + n_prod):
            if g_nj[n, i, n_grains - 1, 0] > 0:
                for r in range(n_grains):
                    for s in range(n_j):
                        u = indices[i, r, s]
                        v = n_rows - n_reac - n_prod + n if products else n_rows - n_reac + n
                        if u > -1:
                            me_mat[u, u] -= g_nj[n, i, r, s]
                            if n < n_reac or products:
                                me_mat[v, u] = g_nj[n, i, r, s]
                            if n < n_reac and not exclude_association:
                                val = f_im[i, n, r, s] * dens_states[n + n_isom, r, s] \
                                      * (2 * j_list[s] + 1) * exp(-e_list[r] * beta)
                                me_mat[u,v] = val
                                me_mat[v,v] -= val

    if neglect_high_energy_collisions:
        for i in range(n_isom):
            for Eind in range(1,n_grains):
                for Jind in range(n_j):
                    ind1 = indices[i, Eind, Jind]
                    if ind1 > -1: #we're handling the loss, but what about production?
                        w = 0.0 #loss due to collisions
                        ktot = 0.0
                        for eind in range(n_grains):
                            ind = indices[i, eind, Jind]
                            if ind > -1 and ind != ind1:
                                w += me_mat[ind,ind1]
                        for j in range(n_isom):
                            if i != j:
                                ind = indices[j,Eind,Jind]
                                ktot += me_mat[ind,ind1]
                        for k in range(n_reac+n_prod):
                            ind = -k-1
                            ktot += me_mat[ind,ind1]
                        if ktot*high_energy_rate_tol > w:
                            for eind in range(n_grains):
                                if ind != ind1:
                                    ind = indices[i, eind, Jind]
                                    me_mat[ind,ind] -= me_mat[ind1,ind]
                                    me_mat[ind1,ind] = 0.0
                                    me_mat[ind,ind1] = 0.0

                            me_mat[ind1,ind1] = -ktot

    if exclude_association:
        return me_mat[:-(n_reac+n_prod),:-(n_reac+n_prod)], indices
    else:
        return me_mat, indices

def states_to_configurations(network,indices,state,exclude_association=False):
    """
    sum full master equation state into total species concentrations
    """
    if exclude_association:
        xs = np.zeros(network.n_isom)
    else:
        xs = np.zeros(network.n_isom+network.n_reac+network.n_prod)
    for i in range(network.n_isom):
        cum = np.float64(0.0)
        for r in range(network.n_grains):
            for s in range(network.n_j):
                index = indices[i,r,s]
                if index == -1:
                    continue
                cum += state[index]
        xs[i] += cum
    if not exclude_association:
        for i in range(network.n_reac+network.n_prod):
            xs[i+network.n_isom] += state[-network.n_reac-network.n_prod+i]
    return xs
