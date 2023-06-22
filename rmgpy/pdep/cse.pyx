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
rate coefficients :math:`k(T,P)` using the chemically-significant eigenvalues
method.
"""

import logging
import scipy.sparse as sparse
import numpy as np
cimport numpy as np
import scipy.linalg
import scipy.sparse as sparse
from libc.math cimport exp, sqrt

import rmgpy.constants as constants
from rmgpy.exceptions import ChemicallySignificantEigenvaluesError
from rmgpy.pdep.me import generate_full_me_matrix, states_to_configurations

################################################################################


def apply_chemically_significant_eigenvalues_method(network, list lumping_order=None, bint neglect_high_energy_collisions=False, double high_energy_rate_tol=0.01):
    """A method for applying the Chemically Significant Eigenvalues approach for solving the master equation."""
    cdef np.ndarray[np.int_t, ndim=1] j_list
    cdef np.ndarray[np.int_t, ndim=3] indices
    cdef np.ndarray[np.float64_t, ndim=1] e_list, s_mat, s_mat_inv, omega0, omega, eq_ratios
    cdef np.ndarray[np.float64_t, ndim=2] me_mat, k, eigen_vectors0, eigen_vectors, z_mat, z_mat_inv, y, x
    cdef np.ndarray[np.float64_t, ndim=3] dens_states
    cdef np.ndarray[np.float64_t, ndim=4] g_nj, pa
    cdef list lumping, unlumping
    cdef double temperature, pressure, ym_b
    cdef int n_isom, n_reac, n_prod, n_grains, n_j, n_chem, n_cse, n_rows
    cdef int i, n, r, s, index

    temperature = network.T
    pressure = network.P
    e_list = network.e_list
    j_list = network.j_list
    dens_states = network.dens_states
    g_nj = network.Gnj
    eq_ratios = network.eq_ratios
    n_isom = network.n_isom
    n_reac = network.n_reac
    n_prod = network.n_prod
    n_grains = network.n_grains
    n_j = network.n_j

    n_grains = len(e_list)
    n_chem = n_isom + n_reac

    ym_b = 1.0e-6 * pressure / (constants.R * temperature)

    # Generate the full master equation matrix
    me_mat, indices = generate_full_me_matrix(network, products=False,
                                              neglect_high_energy_collisions=neglect_high_energy_collisions,
                                              high_energy_rate_tol=high_energy_rate_tol)
    n_rows = me_mat.shape[0]
    me_mat[:, n_rows-n_reac:] *= ym_b

    # Generate symmetrization matrix and its inverse
    s_mat = np.zeros(n_rows, np.float64)
    s_mat_inv = np.zeros_like(s_mat)
    for i in range(n_isom):
        for r in range(n_grains):
            for s in range(n_j):
                index = indices[i, r, s]
                if index > -1:
                    s_mat[index] = sqrt(dens_states[i, r, s] * (2 * j_list[s] + 1) \
                                        * exp(-e_list[r] / (constants.R * temperature)) * eq_ratios[i])
                    s_mat_inv[index] = 1.0 / s_mat[index]
    for n in range(n_reac):
        index = n_rows - n_reac + n
        s_mat[index] = sqrt(eq_ratios[n + n_isom] / ym_b)
        s_mat_inv[index] = 1.0 / s_mat[index]

    # Symmetrize master equation matrix: me_mat = s_mat * m_symm * s_mat_inv
    # Since s_mat and s_mat_inv are diagonal we can do this very efficiently
    for r in range(n_rows):
        for s in range(n_rows):
            me_mat[r, s] = s_mat_inv[r] * me_mat[r, s] * s_mat[s]

    # DEBUG: Check that the matrix has been properly symmetrized
    properly_symmetrized = True
    for r in range(n_rows):
        for s in range(r):
            if me_mat[r, s] != 0:
                if abs(me_mat[r, s] - me_mat[s, r]) > 0.01 * min(me_mat[r, s], me_mat[s, r]) \
                    and max(me_mat[r, s], me_mat[s, r]) > 1e-200:
                        print(r, s, me_mat[r, s], me_mat[s, r])
                        properly_symmetrized = False
    if not properly_symmetrized:
        raise ChemicallySignificantEigenvaluesError('Master equation matrix not properly symmetrized.')

    # Get eigenvalues and eigenvectors
    # We only need the slowest n_chem + 1 eigenmodes, so only compute those
    try:
        # omega0, eigen_vectors0 = scipy.linalg.eigh(me_mat, eigvals=(n_rows - n_chem - 1, n_rows - 1),
        #                                            overwrite_a=True, overwrite_b=True)
        omega0, eigen_vectors0 = scipy.linalg.eigh(me_mat, overwrite_a=True, overwrite_b=True)
    except np.linalg.LinAlgError:
        raise ChemicallySignificantEigenvaluesError('Eigenvalue calculation failed to converge.')

    # We can't assume that eigh returns them in sorted order
    ind = omega0.argsort()

    # Count the number of distinct eigenvalues
    n_cse = 0
    for i in range(n_chem):
        if abs(omega0[ind[-n_chem - 1]] / omega0[ind[-1 - i]]) > 3.0:
            n_cse += 1

    k = np.zeros((n_isom + n_reac + n_prod, n_isom + n_reac + n_prod), np.float64)
    pa = np.zeros((n_isom, n_isom + n_reac, n_grains, n_j), np.float64)

    # Check that we have the correct number of distinct eigenvalues and that
    # there is a zero eigenvalue if there should be (i.e. no product channels)
    # If not, print an error and return
    if n_cse != n_chem:
        logging.error('Could only identify {0:d} distinct eigenvalues, when {1:d} are required.'.format(n_cse, n_chem))
        logging.info('Last IERE = {0:g}    First CSE = {1:g}    Ratio = {2:g}'.format(
                      omega0[ind[-n_chem-1]], omega0[ind[-n_chem]], omega0[ind[-n_chem-1]] / omega0[ind[-n_chem]]))
        if lumping_order is None or len(lumping_order) < n_chem - n_cse:
            # If we don't have a lumping order, then don't try to recover from this situation
            return k, pa
        lumping = lumping_order[0:n_chem - n_cse]
        unlumping = [i for i in range(n_chem) if i not in lumping]

    # elif n_prod == 0 and abs(omega0[ind[-1]]) > 1e-3:
    #    logging.error('Could not identify zero eigenvalue.')
    #    logging.info('Zero CSE = {0:g}    Last CSE = {1:g}    Ratio = {2:g}'.format(
    #        omega0[ind[-1]], omega0[ind[-2]], omega0[ind[-1]] / omega0[ind[-2]]))
    else:
        lumping = []
        unlumping = list(range(n_chem))

    # Extract the chemically-significant eigenvalues and eigenvectors
    omega = omega0.take(ind[-n_cse:])
    eigen_vectors = eigen_vectors0.take(ind[-n_cse:], axis=1)

    # Unsymmetrize the eigenvectors
    for j in range(n_cse):
        eigen_vectors[:, j] *= s_mat

    # Use the "long-time" method to extract the k(T,P) values
    # This method is more numerically robust
    # It also doesn't require finagling with various initial conditions
    # Source: Robertson, Pilling, Jitariu, and Hillier, Phys. Chem. Chem. Phys 9, p. 4085-4097 (2007).
    z_mat = np.zeros((n_cse, n_cse), np.float64)
    z_mat_inv = np.zeros((n_cse, n_cse), np.float64)
    y = np.zeros((n_prod, n_cse), np.float64)
    for j in range(n_cse):
        for i in range(n_isom):
            if i in lumping:
                continue
            i1 = unlumping.index(i)
            for r in range(n_grains):
                for s in range(n_j):
                    index = indices[i, r, s]
                    if index > -1:
                        z_mat[i1,j] += eigen_vectors[index, j]
    for j in range(n_cse):
        for i in range(n_isom):
            for r in range(n_grains):
                for s in range(n_j):
                    index = indices[i, r, s]
                    if index > -1:
                        for n in range(n_prod):
                            y[n, j] += g_nj[n_reac + n, i, r, s] * eigen_vectors[index, j]
    for j in range(n_cse):
        for n in range(n_reac):
            index = n_rows - n_reac + n
            z_mat[n_isom + n - len(lumping), j] += eigen_vectors[index, j]

    z_mat_inv = np.linalg.inv(z_mat)

    for i, m in enumerate(unlumping):
        for j, n in enumerate(unlumping):
            k[m, n] = np.sum(z_mat[i, :] * omega * z_mat_inv[:, j])
    for n in range(n_prod):
        for j, m in enumerate(unlumping):
            k[n_isom + n_reac + n, m] = np.sum(y[n, :] * z_mat_inv[:, j])

    # Compute pa
    x = np.dot(eigen_vectors, z_mat_inv)
    for src, src1 in enumerate(unlumping):
        for i in range(n_isom):
            if i in lumping:
                continue
            for r in range(n_grains):
                for s in range(n_j):
                    index = indices[i, r, s]
                    if index > -1:
                        pa[i, src1, r, s] = x[index, src]

    pa[:, n_isom:, :, :] /= ym_b
    k[:, n_isom:] /= ym_b

    # Return the matrix of k(T,P) values and the pseudo-steady population distributions
    return k, pa

def get_rate_coefficients_CSE_Advanced(network, T, P, neglect_high_energy_collisions=False, high_energy_rate_tol=0.01):
    """
    CSE using the Georgievskii et al. 2013 method https://doi.org/10.1021/jp4060704
    """
    if network.T != T or network.P != P:
        network.set_conditions(T,P)
        network.calculate_equilibrium_ratios()

    lumping_order = None

    temperature = network.T
    pressure = network.P
    beta = 1.0/(constants.R*temperature)
    e_list = network.e_list
    j_list = network.j_list
    dens_states = network.dens_states
    g_nj = network.Gnj
    f_im = network.Fim
    eq_ratios = network.eq_ratios
    n_isom = network.n_isom
    n_reac = network.n_reac
    n_prod = network.n_prod
    n_grains = network.n_grains
    n_j = network.n_j
    n_channels = n_isom+n_reac

    channels = network.isomers+network.reactants

    n_grains = len(e_list)
    n_chem = n_isom #+ n_reac


    # Generate the full master equation matrix
    me_mat, indices = generate_full_me_matrix(network, products=False, exclude_association=True,
                                                  neglect_high_energy_collisions=neglect_high_energy_collisions,
                                                  high_energy_rate_tol=high_energy_rate_tol)
    n_rows = me_mat.shape[0]

    # Generate symmetrization matrix and its inverse
    s_mat = np.zeros(n_rows, np.float64)
    s_mat_inv = np.zeros_like(s_mat)
    for i in range(n_isom):
        for r in range(n_grains):
            for s in range(n_j):
                index = indices[i, r, s]
                if index > -1:
                    s_mat[index] = sqrt(dens_states[i, r, s] * (2 * j_list[s] + 1) \
                                            * exp(-e_list[r] / (constants.R * temperature)) * eq_ratios[i])
                    s_mat_inv[index] = 1.0 / s_mat[index]

    # Symmetrize master equation matrix: me_mat = s_mat * m_symm * s_mat_inv
    # Since s_mat and s_mat_inv are diagonal we can do this very efficiently
    for r in range(n_rows):
        for s in range(n_rows):
            me_mat[r, s] = s_mat_inv[r] * me_mat[r, s] * s_mat[s]

    # DEBUG: Check that the matrix has been properly symmetrized
    properly_symmetrized = True
    for r in range(n_rows):
        for s in range(r):
            if me_mat[r, s] != 0:
                if abs(me_mat[r, s] - me_mat[s, r]) > 0.01 * min(me_mat[r, s], me_mat[s, r]) \
                    and max(me_mat[r, s], me_mat[s, r]) > 1e-200:
                        properly_symmetrized = False
    if not properly_symmetrized:
        raise ChemicallySignificantEigenvaluesError('Master equation matrix not properly symmetrized.')

    try:
        omega0, eigen_vectors0 = scipy.linalg.eigh(me_mat, overwrite_a=True, overwrite_b=True)
    except np.linalg.LinAlgError:
        raise ChemicallySignificantEigenvaluesError('Eigenvalue calculation failed to converge.')

    ind = omega0.argsort()[::-1]

    n_cse = n_chem
    eigratio = abs(omega0[ind[n_chem]]/omega0[ind[n_chem-1]])
    if eigratio < 3.0:
        logging.warning("eigratio of {} was less than 3.0".format(eigratio))

    k = np.zeros((n_isom + n_reac + n_prod, n_isom + n_reac + n_prod), np.float64)

    # Extract the chemically-significant eigenvalues and eigenvectors
    omega = np.real(omega0.take(ind[:]))

    cse_omega = np.real(omega0.take(ind[:n_cse]))
    eigen_vectors = np.real(eigen_vectors0.take(ind[:], axis=1))

    for j in range(eigen_vectors.shape[1]):
        eigen_vectors[:, j] *= s_mat
        eigen_vectors[:, j] /= np.linalg.norm(eigen_vectors[:, j])

    cse_eigen_vectors = eigen_vectors[:,:n_cse]

    bi = np.zeros((n_reac,n_grains,n_j))
    biunnormed = np.zeros((n_reac,n_grains,n_j))
    for n in range(n_reac):
        for eind in range(n_grains):
            for jind in range(n_j):
                bi[n,eind,jind] = dens_states[n+n_isom,eind,jind]*(2 * j_list[jind] + 1) * exp(-e_list[eind] * beta)

    biunnormed[:,:,:] = bi[:,:,:]

    for n in range(n_reac):
        bi[n,:,:] /= np.sum(bi[n,:,:]) #normalize


    Pv = np.zeros((me_mat.shape[0],n_reac))
    k2v = np.zeros((me_mat.shape[0],n_reac))
    for n in range(n_reac):
        for i in range(n_isom):
            for eind in range(n_grains):
                for jind in range(n_j):
                    index = indices[i,eind,jind]
                    if index > -1:
                        Pv[index,n] += f_im[i,n,eind,jind]*bi[n,eind,jind]
                        k2v[index,n] += g_nj[n,i,eind,jind]


    pvlmbd = np.linalg.solve(eigen_vectors,Pv)
    Gvlmbd = np.matmul(eigen_vectors.T, k2v)

    W = np.zeros((n_isom,n_isom))
    for k in range(n_chem):
        W[:,k] = states_to_configurations(network,indices,cse_eigen_vectors[:,k],exclude_association=True)

    Winv = np.linalg.inv(W)


    transformed_omega = np.matmul(W, np.matmul(sparse.diags(omega[:n_chem]), Winv))

    kmat = np.zeros((n_channels,n_channels))

    for i in range(n_isom):
        for j in range(n_isom):
            kmat[i,j] = transformed_omega[i,j]

    for i,react in enumerate(channels):
        for j,prod in enumerate(channels):
            if i == j:
                continue
            if len(react.species) == 2 and len(prod.species) == 2:
                kmat[j,i] = 0.0
                for k,lmbd in enumerate(omega[n_cse:]):
                    kmat[j,i] += -pvlmbd[k+n_cse,i-n_isom]*Gvlmbd[k+n_cse,j-n_isom]/lmbd

            elif len(react.species) == 2:
                kmat[j,i] = 0.0
                kmat[j,i] = np.matmul(W, pvlmbd[:n_cse,i-n_isom])[j]

            elif len(prod.species) == 2:
                kmat[j,i] = 0.0
                kmat[j,i] = np.matmul(Gvlmbd[:n_cse,j-n_isom], Winv[:,i])

    return kmat

def apply_chemically_significant_eigenvalues_method_georgievskii(network, neglect_high_energy_collisions=False, high_energy_rate_tol=0.01):
    return get_rate_coefficients_CSE_Advanced(network, network.T, network.P, neglect_high_energy_collisions=neglect_high_energy_collisions,
                                              high_energy_rate_tol=high_energy_rate_tol)
