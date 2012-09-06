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
rate coefficients :math:`k(T,P)` using the modified strong collision method.
"""

import numpy
cimport numpy

from libc.math cimport exp, log, sqrt

import rmgpy.constants as constants

################################################################################

class ModifiedStrongCollisionError(Exception): 
    """
    An exception raised when the modified strong collsion method is unsuccessful
    for any reason. Pass a string describing the cause of the exceptional 
    behavior.
    """
    pass

################################################################################

cpdef applyModifiedStrongCollisionMethod(network, str efficiencyModel='default'):

    cdef numpy.ndarray[numpy.int_t,ndim=1] Jlist
    cdef numpy.ndarray[numpy.float64_t,ndim=1] Elist, collFreq, collEff, dEdown, E0, Ereac
    cdef numpy.ndarray[numpy.float64_t,ndim=2] A, b, K, x
    cdef numpy.ndarray[numpy.float64_t,ndim=3] densStates
    cdef numpy.ndarray[numpy.float64_t,ndim=4] Kij, Gnj, Fim, pa
    cdef double T, P, E, Emin, val, beta
    cdef int Nisom, Nreac, Nprod, Ngrains, NJ
    cdef int i, j, n, r, s, start, src

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
    Nisom = network.Nisom
    Nreac = network.Nreac
    Nprod = network.Nprod
    Ngrains = network.Ngrains
    NJ = network.NJ
    
    K = numpy.zeros((Nisom+Nreac+Nprod, Nisom+Nreac+Nprod), numpy.float64)
    pa = numpy.zeros((Nisom,Nisom+Nreac,Ngrains,NJ), numpy.float64)

    beta = 1. / (constants.R * T)        # [=] mol/kJ
    
    # Determine the starting grain for the calculation
    Ereac = numpy.zeros(Nisom)
    start = Ngrains
    for i in range(Nisom):
        for r in range(Ngrains):
            if (Kij[:,i,r,0] > 0).any() or (Gnj[:,i,r,0] > 0).any():
                if start > r: start = r
                Ereac[i] = Elist[r]
                break
        else:
            raise ModifiedStrongCollisionError('Unable to determine starting grain; check active-state energies.')
    
    dEdown = numpy.zeros(Nisom)
    for i in range(Nisom):
        dEdown[i] = network.isomers[i].species[0].energyTransferModel.getAlpha(T)
    
    # Compute collision efficiencies
    collEff = numpy.ones(Nisom)
    if efficiencyModel == 'default':
        for i in range(Nisom):
            collEff[i] = network.isomers[i].species[0].energyTransferModel.calculateCollisionEfficiency(T, Elist, Jlist, densStates[i,:,:], E0[i], Ereac[i])
    elif efficiencyModel == 'none':
        pass
    else:
        raise ValueError('Unknown efficiency model "{0}".'.format(efficiencyModel))
    
    # Zero LHS matrix and RHS vectors
    A = numpy.zeros((Nisom,Nisom), numpy.float64)
    b = numpy.zeros((Nisom,Nisom+Nreac), numpy.float64)

    # Iterate over the grains, calculating the PSSA concentrations
    for r in range(start, Ngrains):
        for s in range(NJ):
            # Populate LHS matrix
            # Collisional deactivation
            for i in range(Nisom):
                A[i,i] = -collFreq[i] * collEff[i]
            # Isomerization reactions
            for i in range(Nisom):
                for j in range(i):
                    A[i,j] = Kij[i,j,r,s]
                    A[j,j] -= Kij[i,j,r,s]
                    A[j,i] = Kij[j,i,r,s]
                    A[i,i] -= Kij[j,i,r,s]
            # Dissociation reactions
            for n in range(Nreac+Nprod):
                for j in range(Nisom):
                    A[j,j] -= Gnj[n,j,r,s]
    
            # Populate RHS vectors, one per isomer and reactant
            for i in range(Nisom):
                # Thermal activation via collisions
                b[i,i] = collFreq[i] * collEff[i] * densStates[i,r,s] * (2*Jlist[s]+1) * exp(-Elist[r] * beta)
            for n in range(Nisom, Nisom+Nreac):
                # Chemical activation via association reaction
                for j in range(Nisom):
                    b[j,n] = Fim[j,n-Nisom,r,s] * densStates[n,r,s] * (2*Jlist[s]+1) * exp(-Elist[r] * beta)
                    
            # Solve for steady-state population
            x = -numpy.linalg.solve(A, b)
            for n in range(Nisom+Nreac):
                for i in range(Nisom):
                    pa[i,n,r,s] = x[i,n]
                        
    # Check that our populations are all positive
    if not (pa >= 0).all():
        raise ModifiedStrongCollisionError('A negative steady-state concentration was encountered.')

    # Compute rate coefficients from PSSA concentrations
    for src in range(Nisom+Nreac):
        # Calculate stabilization rates (i.e.) R + R' --> Ai or M --> Ai
        for i in range(Nisom):
            if i != src:
                val = collFreq[i] * collEff[i] * numpy.sum(pa[i,src,:,:])
                K[i,src] += val
                K[src,src] -= val
        # Calculate dissociation rates (i.e.) R + R' --> Bn + Cn or M --> Bn + Cn
        for n in range(Nreac+Nprod):
            for j in range(Nisom):
                if n + Nisom != src:
                    val = numpy.sum(Gnj[n,j,:,:] * pa[j,src,:,:])
                    K[n+Nisom,src] += val
                    K[src,src] -= val
    # To complete pa we need the Boltzmann distribution at low energies
    for i in range(Nisom):
        for r in range(Ngrains):
            for s in range(NJ):
                if pa[i,i,r,s] == 0: pa[i,i,r,s] = densStates[i,r,s] * (2*Jlist[s]+1) * exp(-Elist[r] * beta)
                
    # Return the matrix of k(T,P) values and the pseudo-steady population distributions
    return K, pa
