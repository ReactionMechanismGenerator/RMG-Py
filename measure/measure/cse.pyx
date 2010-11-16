#!/usr/bin/env python
# -*- coding: utf-8 -*-

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

import chempy.constants as constants

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
    cdef numpy.ndarray[numpy.float64_t,ndim=1] S, Sinv, W0, W, C
    cdef numpy.ndarray[numpy.float64_t,ndim=2] eqDist, M, K, V0, V, X, Xinv, dXij
    cdef numpy.ndarray[numpy.float64_t,ndim=3] pa,

    Ngrains = len(Elist)
    Nchem = Nisom + Nreac
    
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

    M = generateFullMEMatrix(Mcoll, Kij, Fim, Gnj, Ngrains, Nisom, Nreac, Nprod, indices)

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
        S[index] = math.sqrt(1.0 * eqRatios[n+Nisom])
        Sinv[index] = 1.0 / S[index]

    # Symmetrize master equation matrix: M = S * Msymm * Sinv
    # Since S and Sinv are diagonal we can do this very efficiently
    for r in range(Nrows):
        for s in range(Nrows):
            M[r,s] = Sinv[r] * M[r,s] * S[s]

    # DEBUG: Check that the matrix has been properly symmetrized
    for r in range(Nrows):
        for s in range(r):
            if M[r,s] != 0:
                if abs(M[r,s] - M[s,r]) / M[r,s] > 0.01:
                    print r, s, M[r,s], M[s,r]
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
        logging.error('Could only identify %i distinct eigenvalues, when %i are required.' % (Ncse, Nchem))
        logging.info('Last IERE = %g    First CSE = %g    Ratio = %g' % (W0[ind[-Nchem-1]], W0[ind[-Nchem]], W0[ind[-Nchem-1]] / W0[ind[-Nchem]]))
    #elif Nprod == 0 and abs(W0[ind[-1]]) > 1e-3:
    #    logging.error('Could not identify zero eigenvalue.')
    #    logging.info('Zero CSE = %g    Last CSE = %g    Ratio = %g' % (W0[ind[-1]], W0[ind[-2]], W0[ind[-1]] / W0[ind[-2]]))
    else:

        # Extract the chemically-significant eigenvalues and eigenvectors
        W = W0.take(ind[-Ncse:])
        V = V0.take(ind[-Ncse:], axis=1)
        
        # Unsymmetrize the eigenvectors and their inverses
        X = V.copy()
        Xinv = V.copy().transpose()
        for j in range(Ncse):
            X[:,j] *= S
            Xinv[j,:] *= Sinv
        
        # Generate set of initial condition vectors, one per isomer and reactant
        # Each initial condition vector corresponds to a case where a Boltzmann
        # distribution of only one isomer or reactant is present
        # The Boltzmann distribution occurs as a result of the fast relaxation of
        # internal energy modes, which is separable from the chemical modes
        eqDist = numpy.zeros((Nrows, Nisom+Nreac), numpy.float64)
        for i in range(Nisom):
            for r in range(Ngrains):
                index = indices[r,i]
                if index > -1: eqDist[index,i] = densStates[i,r] * math.exp(-Elist[r] / constants.R / T)
        for n in range(Nreac):
            index = Nrows - Nreac + n
            eqDist[index,n+Nisom] = 1.0

        # Calculate the phenomenological rate constants
        # This version follows the notation of Miller and Klippenstein
        # Iterate over isomers and reactants to determine k(T,P) with each as the
        # reactant
        for src in range(Nisom+Nreac):
            # dXij contains the change in isomer/reactant/product i as a result
            # of the jth eigenmode, using isomer/reactant n as the starting
            # point
            dXij = numpy.zeros((Nisom+Nreac+Nprod, Ncse), numpy.float64)
            for j in range(Ncse):

                C = X[:,j] * numpy.sum(Xinv[j,:] * eqDist[:,src])

                for i in range(Nisom):
                    for r in range(Ngrains):
                        index = indices[r,i]
                        if index > -1:
                            # Isomers
                            dXij[i,j] += C[index]
                            pa[r,i,src] += C[index]
                            # Product channels
                            for n in range(Nreac, Nreac+Nprod):
                                dXij[Nisom+n,j] += C[index] * Gnj[n,i,r] / W[j]
                for n in range(Nreac):
                    # Reactant channels
                    index = Nrows - Nreac + n
                    dXij[Nisom+n,j] += C[index]

                # Convert the dXij information for eigenmode j into phenomenological
                # rate coefficients
                K[:,src] += W[j] * dXij[:,j]

#        # Use the "long-time" method to extract the k(T,P) values
#        # This method is supposedly more numerically robust
#        # It also doesn't require finagling with various initial conditions
#        cdef numpy.ndarray[numpy.float64_t,ndim=2] Z, Zinv
#        Z = numpy.zeros((Ncse,Ncse), numpy.float64)
#        for j in range(Ncse):
#            for i in range(Nisom):
#                for r in range(Ngrains):
#                    index = indices[r,i]
#                    if index > -1: Z[i,j] += X[index,j]
#            for n in range(Nreac):
#                index = Nrows - Nreac + n
#                Z[Nisom+n,j] += X[index,j]
#        Zinv = numpy.linalg.inv(Z)
#        for i in range(Ncse):
#            for j in range(Ncse):
#                K[i,j] = numpy.sum(Z[i,:] * W * Zinv[:,j])
#        # Note that we still need k(T,P) values for reactions to product channels!
#        # For now we keep those from the initial rate method
        
#    import pylab
#    pylab.semilogy(Elist / 4184, pa[:,0,:])
#    pylab.show()
#    quit()

    # Return the matrix of k(T,P) values and the pseudo-steady population distributions
    return K, pa

################################################################################

