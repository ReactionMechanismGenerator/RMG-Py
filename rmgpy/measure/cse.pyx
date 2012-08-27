################################################################################
#
#   MEASURE - Master Equation Automatic Solver for Unimolecular REactions
#
#   Copyright (c) 2010 by Joshua W. Allen (jwallen@mit.edu)
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
This module provides an implementation of the chemically-significant eigenvalues
method for reducing a master equation model of unimolecular reaction networks to
a set of phenomenological rate coefficients :math:`k(T,P)`.
"""

import math
import numpy
cimport numpy
import numpy.linalg
import logging

cimport rmgpy.constants as constants

from me import generateFullMEMatrix

################################################################################

class ChemicallySignificantEigenvaluesError(Exception):
    """
    An exception raised when the reservoir state method is unsuccessful for
    any reason. Pass a string describing the cause of the exceptional behavior.
    """
    pass

################################################################################

def applyChemicallySignificantEigenvaluesMethod(double T, double P,
    numpy.ndarray[numpy.float64_t,ndim=1] Elist,
    numpy.ndarray[numpy.float64_t,ndim=2] densStates,
    numpy.ndarray[numpy.float64_t,ndim=3] Mcoll,
    numpy.ndarray[numpy.float64_t,ndim=3] Kij,
    numpy.ndarray[numpy.float64_t,ndim=3] Fim,
    numpy.ndarray[numpy.float64_t,ndim=3] Gnj,
    numpy.ndarray[numpy.float64_t,ndim=1] eqRatios,
    int Nisom, int Nreac, int Nprod):
    """
    Use the chemically-significant eigenvalues method to reduce the master
    equation model to a set of phenomenological rate coefficients :math:`k(T,P)`
    and a set of time-independent population vectors :math:`\\vector{u}_{ij}`
    and :math:`\\vector{v}_{im}`. Inputs are the temperature `T` in K; pressure
    `P` in Pa; list of energy grains `Elist` in J/mol; dimensionless densities
    of states for each isomer and reactant channel `densStates`; collision
    matrix `Mcoll` for each isomer; isomerization, association, and dissociation
    microcanonical rate coefficients `Kij`, `Fim`, and `Gnj`, respectively;
    energies of the first reactive grain for each isomer `Ereac` in J/mol;
    and the numbers of isomers, reactant channels, and product channels `Nisom`,
    `Nreac`, and `Nprod`, respectively.
    """

    cdef int Ngrains, Nchem, Nrows, Ncse, index, i, n, r, s
    cdef numpy.ndarray[numpy.int_t,ndim=2] indices
    cdef numpy.ndarray[numpy.float64_t,ndim=1] S, Sinv, W0, W
    cdef numpy.ndarray[numpy.float64_t,ndim=2] M, K, V0, V, Z, Zinv, Y, X
    cdef numpy.ndarray[numpy.float64_t,ndim=3] pa
    cdef double ymB
    
    Ngrains = len(Elist)
    Nchem = Nisom + Nreac
    
    # The concentration of the major species for all bimolecular channels
    # This is the species whose concentration is constant (pseudo-first order
    # approximation)
    # Tweaking this number may improve the numerics somewhat (not guaranteed)
    ymB = 1.0
    
    # Construct accounting matrix
    Nrows = 0
    # Row is grain number, column is well number, value is index into matrix
    indices = -numpy.ones((Ngrains,Nisom), numpy.int)
    for r in range(Ngrains):
        for i in range(Nisom):
            if densStates[i,r] > 0:
                indices[r,i] = Nrows
                Nrows += 1
    # Reactant wells are appended to matrix, one row each
    # Product wells are NOT added because we are neglecting reassociation
    Nrows += Nreac

#    # Multiply association rates by Boltzmann distribution to get association fluxes
#    Fim = Fim.copy()
#    for i in range(Nisom):
#        for n in range(Nreac):
#            Fim[i,n,:] *= densStates[n+Nisom,:] * numpy.exp(-Elist[:] / constants.R / T)

    M = generateFullMEMatrix(Mcoll, Kij, Fim * ymB, Gnj, Ngrains, Nisom, Nreac, Nprod, indices)

    # Generate symmetrization matrix and its inverse
    S = numpy.zeros(Nrows, numpy.float64)
    Sinv = numpy.zeros_like(S)
    for i in range(Nisom):
        for r in range(Ngrains):
            index = indices[r,i]
            if index > -1:
                S[index] = math.sqrt(densStates[i,r] * math.exp(-Elist[r] / constants.R / T) * eqRatios[i])
                Sinv[index] = 1.0 / S[index]
    for n in range(Nreac):
        index = Nrows - Nreac + n
        S[index] = math.sqrt(eqRatios[n+Nisom] / ymB)
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
                if abs(M[r,s] - M[s,r]) / M[r,s] > 0.01:
                    if M[r,s] > 1e-200 or M[s,r] > 1e-200:
                        print r, s, M[r,s], M[s,r]
                        properlySymmetrized = False
    if not properlySymmetrized:
        raise ChemicallySignificantEigenvaluesError('Master equation matrix not properly symmetrized.')

    # Get eigenvalues and eigenvectors
    try:
        W0, V0 = numpy.linalg.eigh(M)
    except numpy.linalg.LinAlgError:
        raise ChemicallySignificantEigenvaluesError('Eigenvalue calculation failed to converge.')
    
    # We can't assume that eigh returns them in sorted order
    ind = W0.argsort()

    # Count the number of distinct eigenvalues
    Ncse = 0
    for i in range(Nchem):
        if abs(W0[ind[-Nchem-1]] / W0[ind[-1-i]]) > 3.0: Ncse += 1

    K = numpy.zeros((Nisom+Nreac+Nprod, Nisom+Nreac+Nprod), numpy.float64)
    pa = numpy.zeros((Ngrains, Nisom, Nisom+Nreac), numpy.float64)

    # Check that we have the correct number of distinct eigenvalues and that
    # there is a zero eigenvalue if there should be (i.e. no product channels)
    # If not, print an error and return
    if Ncse != Nchem:
        logging.error('Could only identify {0:d} distinct eigenvalues, when {1:d} are required.'.format(Ncse, Nchem))
        logging.info('Last IERE = {0:g}    First CSE = {1:g}    Ratio = {2:g}'.format(W0[ind[-Nchem-1]], W0[ind[-Nchem]], W0[ind[-Nchem-1]] / W0[ind[-Nchem]]))
    #elif Nprod == 0 and abs(W0[ind[-1]]) > 1e-3:
    #    logging.error('Could not identify zero eigenvalue.')
    #    logging.info('Zero CSE = {0:g}    Last CSE = {1:g}    Ratio = {2:g}'.format(W0[ind[-1]], W0[ind[-2]], W0[ind[-1]] / W0[ind[-2]]))
    else:

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
        Y = numpy.zeros((Nprod,Ncse), numpy.float64)
        for j in range(Ncse):
            for i in range(Nisom):
                for r in range(Ngrains):
                    index = indices[r,i]
                    if index > -1: 
                        Z[i,j] += V[index,j]
                        for n in range(Nprod):
                            Y[n,j] += Gnj[Nreac+n,i,r] * V[index,j]
            for n in range(Nreac):
                index = Nrows - Nreac + n
                Z[Nisom+n,j] += V[index,j]
        Zinv = numpy.linalg.inv(Z)
        for i in range(Ncse):
            for j in range(Ncse):
                K[i,j] = numpy.sum(Z[i,:] * W * Zinv[:,j])
        for n in range(Nprod):
            for j in range(Ncse):
                K[Nisom+Nreac+n,j] = numpy.sum(Y[n,:] * Zinv[:,j])
        
        # Compute pa
        X = numpy.dot(V, Zinv)
        for src in range(Nisom+Nreac):
            for i in range(Nisom):
                for r in range(Ngrains):
                    index = indices[r,i]
                    if index > -1:
                        pa[r,i,src] = X[index,src]
     
    K[:,Nisom:] = K[:,Nisom:] / ymB
        
#    import pylab
#    pylab.semilogy(Elist / 4184, pa[:,0,:])
#    pylab.show()
#    quit()

    # Return the matrix of k(T,P) values and the pseudo-steady population distributions
    return K, pa

################################################################################

