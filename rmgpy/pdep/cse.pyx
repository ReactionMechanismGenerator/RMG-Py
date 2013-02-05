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
Contains functionality for computing pressure-dependent phenomenological
rate coefficients :math:`k(T,P)` using the chemically-significant eigenvalues
method.
"""

import numpy
cimport numpy
import logging
import scipy.linalg

from libc.math cimport exp, log, sqrt

import rmgpy.constants as constants

from rmgpy.pdep.me import generateFullMEMatrix

################################################################################

class ChemicallySignificantEigenvaluesError(Exception):
    """
    An exception raised when the reservoir state method is unsuccessful for
    any reason. Pass a string describing the cause of the exceptional behavior.
    """
    pass

################################################################################

def applyChemicallySignificantEigenvaluesMethod(network, list lumpingOrder=None):

    cdef numpy.ndarray[numpy.int_t,ndim=1] Jlist
    cdef numpy.ndarray[numpy.int_t,ndim=3] indices
    cdef numpy.ndarray[numpy.float64_t,ndim=1] Elist, S, Sinv, W0, W, eqRatios
    cdef numpy.ndarray[numpy.float64_t,ndim=2] M, K, V0, V, Z, Zinv, Y, X
    cdef numpy.ndarray[numpy.float64_t,ndim=3] densStates
    cdef numpy.ndarray[numpy.float64_t,ndim=4] Kij, Gnj, Fim, pa
    cdef list lumping, unlumping
    cdef double T, P, ymB
    cdef int Nisom, Nreac, Nprod, Ngrains, NJ, Nchem, Ncse, Nrows
    cdef int i, n, r, s, index

    T = network.T
    P = network.P
    Elist = network.Elist
    Jlist = network.Jlist
    densStates = network.densStates
    collFreq = network.collFreq
    Kij = network.Kij
    Fim = network.Fim
    Gnj = network.Gnj
    E0 = network.E0
    eqRatios = network.eqRatios
    Nisom = network.Nisom
    Nreac = network.Nreac
    Nprod = network.Nprod
    Ngrains = network.Ngrains
    NJ = network.NJ
    
    Ngrains = len(Elist)
    Nchem = Nisom + Nreac
    
    ymB = 1.0e-6 * P / constants.R / T
    
    # Generate the full master equation matrix
    M, indices = generateFullMEMatrix(network, products=False)
    Nrows = M.shape[0]
    M[:,Nrows-Nreac:] *= ymB
    
    # Generate symmetrization matrix and its inverse
    S = numpy.zeros(Nrows, numpy.float64)
    Sinv = numpy.zeros_like(S)
    for i in range(Nisom):
        for r in range(Ngrains):
            for s in range(NJ):
                index = indices[i,r,s]
                if index > -1:
                    S[index] = sqrt(densStates[i,r,s] * (2*Jlist[s]+1) * exp(-Elist[r] / constants.R / T) * eqRatios[i])
                    Sinv[index] = 1.0 / S[index]
    for n in range(Nreac):
        index = Nrows - Nreac + n
        S[index] = sqrt(eqRatios[n+Nisom] / ymB)
        Sinv[index] = 1.0 / S[index]

    # Symmetrize master equation matrix: M = S * Msymm * Sinv
    # Since S and Sinv are diagonal we can do this very efficiently
    for r in range(Nrows):
        for s in range(Nrows):
            M[r,s] = Sinv[r] * M[r,s] * S[s]

    # DEBUG: Check that the matrix has been properly symmetrized
    properlySymmetrized = True
    for r in range(Nrows):
        for s in range(r):
            if M[r,s] != 0:
                if abs(M[r,s] - M[s,r]) > 0.01 * M[r,s]:
                    if M[r,s] > 1e-200 or M[s,r] > 1e-200:
                        print r, s, M[r,s], M[s,r]
                        properlySymmetrized = False
    if not properlySymmetrized:
        raise ChemicallySignificantEigenvaluesError('Master equation matrix not properly symmetrized.')

    # Get eigenvalues and eigenvectors
    # We only need the slowest Nchem + 1 eigenmodes, so only compute those
    try:
        #W0, V0 = scipy.linalg.eigh(M, eigvals=(Nrows-Nchem-1,Nrows-1), overwrite_a=True, overwrite_b=True)
        W0, V0 = scipy.linalg.eigh(M, overwrite_a=True, overwrite_b=True)
    except numpy.linalg.LinAlgError:
        raise ChemicallySignificantEigenvaluesError('Eigenvalue calculation failed to converge.')
    
    # We can't assume that eigh returns them in sorted order
    ind = W0.argsort()

    # Count the number of distinct eigenvalues
    Ncse = 0
    for i in range(Nchem):
        if abs(W0[ind[-Nchem-1]] / W0[ind[-1-i]]) > 3.0: Ncse += 1

    K = numpy.zeros((Nisom+Nreac+Nprod, Nisom+Nreac+Nprod), numpy.float64)
    pa = numpy.zeros((Nisom, Nisom+Nreac, Ngrains, NJ), numpy.float64)


    
    
    # Check that we have the correct number of distinct eigenvalues and that
    # there is a zero eigenvalue if there should be (i.e. no product channels)
    # If not, print an error and return
    if Ncse != Nchem:
        logging.error('Could only identify {0:d} distinct eigenvalues, when {1:d} are required.'.format(Ncse, Nchem))
        logging.info('Last IERE = {0:g}    First CSE = {1:g}    Ratio = {2:g}'.format(W0[ind[-Nchem-1]], W0[ind[-Nchem]], W0[ind[-Nchem-1]] / W0[ind[-Nchem]]))
        if lumpingOrder is None or len(lumpingOrder) < Nchem - Ncse:
            # If we don't have a lumping order, then don't try to recover from
            # this situation
            return K, pa
        lumping = lumpingOrder[0:Nchem-Ncse]
        unlumping = [i for i in range(Nchem) if i not in lumping]
    
    #elif Nprod == 0 and abs(W0[ind[-1]]) > 1e-3:
    #    logging.error('Could not identify zero eigenvalue.')
    #    logging.info('Zero CSE = {0:g}    Last CSE = {1:g}    Ratio = {2:g}'.format(W0[ind[-1]], W0[ind[-2]], W0[ind[-1]] / W0[ind[-2]]))
    else:
        lumping = []
        unlumping = range(Nchem)

    # Extract the chemically-significant eigenvalues and eigenvectors
    W = W0.take(ind[-Ncse:])
    V = V0.take(ind[-Ncse:], axis=1)
    
    # Unsymmetrize the eigenvectors
    for j in range(Ncse):
        V[:,j] *= S

    # Use the "long-time" method to extract the k(T,P) values
    # This method is more numerically robust
    # It also doesn't require finagling with various initial conditions
    # Source: Robertson, Pilling, Jitariu, and Hillier, Phys. Chem. Chem. Phys 9, p. 4085-4097 (2007).
    Z = numpy.zeros((Ncse,Ncse), numpy.float64)
    Zinv = numpy.zeros((Ncse,Ncse), numpy.float64)
    Y = numpy.zeros((Nprod,Ncse), numpy.float64)
    for j in range(Ncse):
        for i in range(Nisom):
            if i in lumping: continue
            i1 = unlumping.index(i)
            for r in range(Ngrains):
                for s in range(NJ):
                    index = indices[i,r,s]
                    if index > -1: 
                        Z[i1,j] += V[index,j]
    for j in range(Ncse):
        for i in range(Nisom):
            for r in range(Ngrains):
                for s in range(NJ):
                    index = indices[i,r,s]
                    if index > -1: 
                        for n in range(Nprod):
                            Y[n,j] += Gnj[Nreac+n,i,r,s] * V[index,j]
    for j in range(Ncse):
        for n in range(Nreac):
            index = Nrows - Nreac + n
            Z[Nisom+n-len(lumping),j] += V[index,j]
    
    Zinv = numpy.linalg.inv(Z)
    
    for i, m in enumerate(unlumping):
        for j, n in enumerate(unlumping):
            K[m,n] = numpy.sum(Z[i,:] * W * Zinv[:,j])
    for n in range(Nprod):
        for j, m in enumerate(unlumping):
            K[Nisom+Nreac+n,m] = numpy.sum(Y[n,:] * Zinv[:,j])
        
    # Compute pa
    X = numpy.dot(V, Zinv)
    for src, src1 in enumerate(unlumping):
        for i in range(Nisom):
            if i in lumping: continue
            for r in range(Ngrains):
                for s in range(NJ):
                    index = indices[i,r,s]
                    if index > -1:
                        pa[i,src1,r,s] = X[index,src]
    
    pa[:,Nisom:,:,:] /= ymB
    K[:,Nisom:] /= ymB

    # Return the matrix of k(T,P) values and the pseudo-steady population distributions
    return K, pa
