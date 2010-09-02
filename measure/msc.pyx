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
This module provides an implementation of the modified strong collsion method 
for reducing a master equation model of unimolecular reaction networks to a set
of phenomenological rate coefficients :math:`k(T,P)`.
"""

import math
import numpy
cimport numpy
import cython

import chempy.constants as constants

################################################################################

class ModifiedStrongCollisionError(Exception): 
    """
    An exception raised when the modified strong collsion method is unsuccessful
    for any reason. Pass a string describing the cause of the exceptional 
    behavior.
    """
    pass

################################################################################

@cython.boundscheck(False)
def applyModifiedStrongCollisionMethod(double T, double P,
    numpy.ndarray[numpy.float64_t,ndim=1] Elist,
    numpy.ndarray[numpy.float64_t,ndim=2] densStates,
    numpy.ndarray[numpy.float64_t,ndim=1] collFreq,
    numpy.ndarray[numpy.float64_t,ndim=3] Kij,
    numpy.ndarray[numpy.float64_t,ndim=3] Fim,
    numpy.ndarray[numpy.float64_t,ndim=3] Gnj,
    numpy.ndarray[numpy.float64_t,ndim=1] Ereac,
    int Nisom, int Nreac, int Nprod):
    """
    Use the modified strong collsion method to reduce the master equation model
    to a set of phenomenological rate coefficients :math:`k(T,P)` and a set of
    time-independent population vectors :math:`\\vector{u}_{ij}` and
    :math:`\\vector{v}_{im}`. Inputs are the temperature `T` in K; pressure `P`
    in Pa; list of energy grains `Elist` in J/mol; dimensionless densities of 
    states for each isomer and reactant channel `densStates`; modified collision
    frequencies `collFreq` for each isomer in s^-1; isomerization, association, 
    and dissociation microcanonical rate coefficients `Kij`, `Fim`, and `Gnj`, 
    respectively; energies of the first reactive grain for each isomer `Ereac`
    in J/mol; and the numbers of isomers, reactant channels, and product 
    channels `Nisom`, `Nreac`, and `Nprod`, respectively.
    """
    
    cdef int Ngrains, start, i, j, n, r, src
    cdef double E, Emin
    cdef numpy.ndarray[numpy.float64_t,ndim=2] A, b, K
    cdef numpy.ndarray[numpy.float64_t,ndim=3] pa

    Ngrains = len(Elist)

    K = numpy.zeros((Nisom+Nreac+Nprod, Nisom+Nreac+Nprod), numpy.float64)
    pa = numpy.zeros((Ngrains,Nisom,Nisom+Nreac), numpy.float64)

    # Determine the starting grain for the calculation based on the
    # active-state cutoff energy
    Emin = numpy.min(Ereac); start = -1
    for i in range(len(Elist)):
        if Elist[i] > Emin:
            start = i
            break
    if start < 0:
        raise ModifiedStrongCollisionError('Unable to determine starting grain; check active-state energies.')

    # Iterate over the grains, calculating the PSSA concentrations
    for r in range(start, Ngrains):

        # Zero LHS matrix and RHS vectors
        A = numpy.zeros((Nisom,Nisom), numpy.float64)
        b = numpy.zeros((Nisom,Nisom+Nreac), numpy.float64)

        # Populate LHS matrix
        # Collisional deactivation
        for i in range(Nisom):
            A[i,i] -= collFreq[i]
        # Isomerization reactions
        for i in range(Nisom):
            for j in range(i):
                A[i,j] = Kij[i,j,r]
                A[j,j] -= Kij[i,j,r]
                A[j,i] = Kij[j,i,r]
                A[i,i] -= Kij[j,i,r]
        # Dissociation reactions
        for n in range(Nreac+Nprod):
            for j in range(Nisom):
                A[j,j] -= Gnj[n,j,r]

        # Populate RHS vectors, one per isomer and reactant
        for i in range(Nisom):
            # Thermal activation via collisions
            b[i,i] = collFreq[i] * densStates[i,r] * math.exp(-Elist[r] / constants.R / T)
        for n in range(Nisom, Nisom+Nreac):
            # Chemical activation via association reaction
            b[:,n] = Fim[:,n-Nisom,r] #* (densStates[n,r] * math.exp(-Elist[r] / constants.R / T))

        # Solve for steady-state population
        pa[r,:,:] = -numpy.linalg.solve(A, b)
        
            
    # Check that our populations are all positive
    if not (pa >= 0).all():
        raise ModifiedStrongCollisionError('A negative steady-state concentration was encountered.')

    # Compute rate coefficients from PSSA concentrations
    for src in range(Nisom+Nreac):
        # Calculate stabilization rates (i.e.) R + R' --> Ai or M --> Ai
        for i in range(Nisom):
            if i != src:
                val = collFreq[i] * numpy.sum(pa[:,i,src])
                K[i,src] += val
                K[src,src] -= val
        # Calculate dissociation rates (i.e.) R + R' --> Bn + Cn or M --> Bn + Cn
        for n in range(Nreac+Nprod):
            for j in range(Nisom):
                if n + Nisom != src:
                    val = numpy.sum(Gnj[n,j,:] * pa[:,j,src])
                    K[n+Nisom,src] += val
                    K[src,src] -= val
    
    # To complete pa we need the Boltzmann distribution at low energies
    for i in range(Nisom):
        for r in range(Ngrains):
            if pa[r,i,i] == 0: pa[r,i,i] = densStates[i,r] * math.exp(-Elist[r] / constants.R / T)

    # Return the matrix of k(T,P) values and the pseudo-steady population distributions
    return K, pa
    