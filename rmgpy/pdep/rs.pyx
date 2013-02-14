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
rate coefficients :math:`k(T,P)` using the reservoir state method.
"""

import numpy
cimport numpy
import scipy.linalg

from libc.math cimport exp, log, sqrt

import rmgpy.constants as constants

################################################################################

class ReservoirStateError(Exception): 
    """
    An exception raised when the reservoir state method is unsuccessful for
    any reason. Pass a string describing the cause of the exceptional behavior.
    """
    pass

################################################################################

cpdef applyReservoirStateMethod(network):

    cdef numpy.ndarray[numpy.int_t,ndim=1] Jlist
    cdef numpy.ndarray[numpy.int_t,ndim=2] Nres, Nact
    cdef numpy.ndarray[numpy.int_t,ndim=3] indices
    cdef numpy.ndarray[numpy.float64_t,ndim=1] Elist, collFreq, E0
    cdef numpy.ndarray[numpy.float64_t,ndim=2] L, Z, X, K
    cdef numpy.ndarray[numpy.float64_t,ndim=3] densStates, eqDist
    cdef numpy.ndarray[numpy.float64_t,ndim=4] Kij, Gnj, Fim, pa
    cdef numpy.ndarray[numpy.float64_t,ndim=5] Mcoll
    cdef list ind
    cdef double T, P, E, tol, y, dfactor, beta
    cdef int Nisom, Nreac, Nprod, Ngrains, NJ, bandwidth, halfbandwidth, width, width0
    cdef int i, j, n, r, s, u, v, row, iter

    T = network.T
    P = network.P
    Elist = network.Elist
    Jlist = network.Jlist
    densStates = network.densStates
    collFreq = network.collFreq
    Mcoll = network.Mcoll
    Kij = network.Kij
    Fim = network.Fim
    Gnj = network.Gnj
    E0 = network.E0
    Nisom = network.Nisom
    Nreac = network.Nreac
    Nprod = network.Nprod
    Ngrains = network.Ngrains
    NJ = network.NJ
    
    beta = 1. / (constants.R * T)        # [=] mol/kJ

    K = numpy.zeros((Nisom+Nreac+Nprod, Nisom+Nreac+Nprod), numpy.float64)
    pa = numpy.zeros((Nisom,Nisom+Nreac,Ngrains,NJ), numpy.float64)

    # Determine the reservoir cutoff grain for each isomer
    # Start by simply placing it at the lowest reactive grain
    Nres = numpy.zeros((Nisom,NJ), numpy.int)
    for i in range(Nisom):
        for s in range(NJ):
            for r in range(Ngrains):
                if densStates[i,r,s] != 0 and ((Kij[:,i,r,s] > 0).any() or (Gnj[:,i,r,s] > 0).any()):
                    # We need at least one reservoir grain for the RS method to be successful
                    if r == 0 or densStates[i,r-1,s] == 0:
                        Nres[i,s] = r + 1
                    else:
                        Nres[i,s] = r
                    break
    Nact = Ngrains - Nres
    
    # Determine equilibrium distributions
    eqDist = numpy.zeros((Nisom+Nreac,Ngrains,NJ), numpy.float64)
    for i in range(Nisom+Nreac):
        for s in range(NJ):
            eqDist[i,:,s] = densStates[i,:,s] * (2*Jlist[s]+1) * numpy.exp(-Elist * beta)

    # Determine pseudo-steady state populations of active state
    row = 0
    indices = -numpy.ones((Nisom,Ngrains,NJ), numpy.int)
    for r in range(Ngrains):
        for s in range(NJ):
            for i in range(Nisom):
                if r >= Nres[i,s]:
                    indices[i,r,s] = row
                    row += 1
    
    # Choose the half-bandwidth using the deepest isomer well
    width = 0
    tol = 1e-12
    for i in range(Nisom):
        for s in range(NJ):
            r = Nres[i,s]
            if Mcoll[i,r,s,r,s] == 0: continue
            ratio = numpy.abs(Mcoll[i,:,s,r,s] / Mcoll[i,r,s,r,s])
            ind = [j for j,y in enumerate(ratio) if y > tol]
            if len(ind) > 0:
                width0 = max(r - min(ind), max(ind) - r)
                if width0 > width:
                    width = width0
    if width == 0:
        raise ReservoirStateError('Unable to determine half-bandwidth for active-state matrix; the wells may be too shallow to use the RS method.')
    halfbandwidth = (width + 1) * Nisom * NJ - Nisom
    bandwidth = 2 * halfbandwidth + 1
    
    # Populate active-state matrix and source vectors
    L = numpy.zeros((bandwidth,numpy.sum(Nact)), numpy.float64)
    Z = numpy.zeros((numpy.sum(Nact),Nisom+Nreac), numpy.float64)
    # Collisional terms
    for i in range(Nisom):
        for u in range(NJ):
            for v in range(NJ):
                for r in range(Nres[i,u], Ngrains):
                    for s in range(max(Nres[i,v], r-width), min(Ngrains, r+width+1)):
                        L[halfbandwidth + indices[i,r,u] - indices[i,s,v], indices[i,s,v]] = Mcoll[i,r,u,s,v]
                    Z[indices[i,r,u],i] = numpy.sum(Mcoll[i,r,u,0:Nres[i,u],v] * eqDist[i,0:Nres[i,u],v])

    # Isomerization terms
    for i in range(Nisom):
        for j in range(i):
            for u in range(NJ):
                for r in range(max(Nres[i,u], Nres[j,u]), Ngrains):
                    L[halfbandwidth + indices[j,r,u] - indices[i,r,u], indices[i,r,u]] = Kij[j,i,r,u]
                    L[halfbandwidth, indices[i,r,u]] -= Kij[j,i,r,u]
                    L[halfbandwidth + indices[i,r,u] - indices[j,r,u], indices[j,r,u]] = Kij[i,j,r,u]
                    L[halfbandwidth, indices[j,r,u]] -= Kij[i,j,r,u]
    # Dissociation/association terms
    for i in range(Nisom):
        for n in range(Nreac+Nprod):
            for u in range(NJ):
                for r in range(Nres[i,u], Ngrains):
                    L[halfbandwidth, indices[i,r,u]] -= Gnj[n,i,r,u]
        for n in range(Nreac):
            for u in range(NJ):
                for r in range(Nres[i,u], Ngrains):
                    Z[indices[i,r,u], n+Nisom] = Fim[i,n,r,u] * eqDist[n+Nisom,r,u]
    
    # Solve for pseudo-steady state populations of active state
    X = scipy.linalg.solve_banded((halfbandwidth,halfbandwidth), L, -Z, overwrite_ab=True, overwrite_b=True)
    for i in range(Nisom):
        for u in range(NJ):
            for r in range(Nres[i,u], Ngrains):
                for n in range(Nisom+Nreac):
                    pa[i,n,r,u] = X[indices[i,r,u], n]
    
    # Double-check to ensure that we have all positive populations
    if not (pa >= 0).all():
        raise ReservoirStateError('A negative steady-state population was encountered.')

    # Put the reservoir populations into pa as well
    for i in range(Nisom):
        for u in range(NJ):
            for r in range(Nres[i,u]):
                pa[i,i,r,u] = eqDist[i,r,u]

    # Determine the phenomenological rate coefficients using the general procedure
    # This should be exactly the same as the procedure below, which is more
    # specific to the RS method
    # Previously it was noted that this more general approach was more robust;
    # however, more recently it seems that this is no longer the case
    #K = computeRateCoefficients(Mcoll, Kij, Fim, Gnj, pa, Nisom, Nreac, Nprod)

    # Determine the phenomenological rate coefficients
    K = numpy.zeros((Nisom+Nreac+Nprod, Nisom+Nreac+Nprod), numpy.float64)
    # Rows relating to isomers
    for i in range(Nisom):
        for u in range(NJ):
            for v in range(NJ):
                # Collisional rearrangement within the reservoir of isomer i
                K[i,i] = K[i,i] + numpy.sum(numpy.dot(Mcoll[i,0:Nres[i,u],u,0:Nres[i,v],v], eqDist[i,0:Nres[i,v],v]))
                # Isomerization from isomer j to isomer i
                for j in range(Nisom):
                    K[i,j] = K[i,j] + numpy.sum(numpy.dot(Mcoll[i,0:Nres[i,u],u,Nres[i,v]:Ngrains,v], pa[i,j,Nres[i,v]:Ngrains,v]))
                # Association from reactant n to isomer i
                for n in range(Nisom, Nisom+Nreac):
                    K[i,n] = K[i,n] + numpy.sum(numpy.dot(Mcoll[i,0:Nres[i,u],u,Nres[i,v]:Ngrains,v], pa[i,n,Nres[i,v]:Ngrains,v]))
    # Rows relating to reactants
    for n in range(Nreac):
        # Association loss
        for i in range(Nisom):
            K[Nisom+n,Nisom+n] = K[Nisom+n,Nisom+n] - numpy.sum(Fim[i,n,:,:] * eqDist[n+Nisom,:,:])
        # Reaction from isomer or reactant j to reactant n
        for j in range(Nisom+Nreac):
            for i in range(Nisom):
                for u in range(NJ):
                    K[Nisom+n,j] = K[Nisom+n,j] + numpy.sum(Gnj[n,i,Nres[i,u]:Ngrains,u] * pa[i,j,Nres[i,u]:Ngrains,u])
    # Rows relating to products
    for n in range(Nreac, Nreac+Nprod):
        # Reaction from isomer or reactant j to product n
        for j in range(Nisom+Nreac):
            for i in range(Nisom):
                for u in range(NJ):
                    K[Nisom+n,j] = K[Nisom+n,j] + numpy.sum(Gnj[n,i,Nres[i,u]:Ngrains,u] * pa[i,j,Nres[i,u]:Ngrains,u])

    # Ensure matrix is conservative
    for n in range(Nisom+Nreac):
        K[n,n] = K[n,n] - numpy.sum(K[:,n])

    # Return the matrix of k(T,P) values and the pseudo-steady population distributions
    return K, pa
