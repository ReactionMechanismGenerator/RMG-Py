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
This module provides an implementation of the reservoir state method for
reducing a master equation model of unimolecular reaction networks to a set of
phenomenological rate coefficients :math:`k(T,P)`.
"""

import math
import numpy
cimport numpy
import scipy.linalg
import cython

cimport rmgpy.constants as constants
from me cimport computeRateCoefficients

################################################################################

class ReservoirStateError(Exception): 
    """
    An exception raised when the reservoir state method is unsuccessful for
    any reason. Pass a string describing the cause of the exceptional behavior.
    """
    pass

################################################################################

@cython.boundscheck(False)
def applyReservoirStateMethod(double T, double P,
    numpy.ndarray[numpy.float64_t,ndim=1] Elist,
    numpy.ndarray[numpy.float64_t,ndim=2] densStates,
    numpy.ndarray[numpy.float64_t,ndim=3] Mcoll,
    numpy.ndarray[numpy.float64_t,ndim=3] Kij,
    numpy.ndarray[numpy.float64_t,ndim=3] Fim,
    numpy.ndarray[numpy.float64_t,ndim=3] Gnj,
    numpy.ndarray[numpy.float64_t,ndim=1] Ereac,
    int Nisom, int Nreac, int Nprod):
    """
    Use the reservoir state method to reduce the master equation model to a
    set of phenomenological rate coefficients :math:`k(T,P)` and a set of
    time-independent population vectors :math:`\\vector{u}_{ij}` and
    :math:`\\vector{v}_{im}`. Inputs are the temperature `T` in K; pressure `P`
    in Pa; list of energy grains `Elist` in J/mol; dimensionless densities of 
    states for each isomer and reactant channel `densStates`; collision
    matrix `Mcoll` for each isomer; isomerization, association, and dissociation
    microcanonical rate coefficients `Kij`, `Fim`, and `Gnj`, respectively;
    energies of the first reactive grain for each isomer `Ereac` in J/mol;
    and the numbers of isomers, reactant channels, and product channels `Nisom`,
    `Nreac`, and `Nprod`, respectively. The method involves a significant linear
    solve, which is accelerated by taking advantage of the bandedness of the
    active-state matrix. The nonreactive grains are placed in the reservoir,
    while the reactive grains are placed in the active-state.
    """

    cdef int Ngrains, bandwidth, halfbandwidth, row, i, j, n, r, s
    cdef double E, tol
    cdef list ind
    cdef numpy.ndarray[numpy.int_t,ndim=1] Nres, Nact
    cdef numpy.ndarray[numpy.int_t,ndim=2] indices
    cdef numpy.ndarray[numpy.float64_t,ndim=1] ratio
    cdef numpy.ndarray[numpy.float64_t,ndim=2] eqDist, L, Z, X, K
    cdef numpy.ndarray[numpy.float64_t,ndim=3] pa

    Ngrains = len(Elist)

    # Determine the reservoir cutoff grain for each isomer
    # Start by simply placing it at the lowest reactive grain
    Nres = numpy.zeros(Nisom, numpy.int)
    for i in range(Nisom):
        for r in range(Ngrains):
            if densStates[i,r] != 0 and Elist[r] > Ereac[i] and Nres[i] == 0:
                # We need at least one reservoir grain for the RS method to be successful
                if r == 0 or densStates[i,r-1] == 0:
                    Nres[i] = r + 1
                else:
                    Nres[i] = r
                break
    Nact = Ngrains - Nres
    
    # Determine equilibrium distributions
    eqDist = numpy.zeros((Nisom+Nreac,Ngrains), numpy.float64)
    for i in range(Nisom+Nreac):
        eqDist[i,:] = densStates[i,:] * numpy.exp(-Elist / constants.R / T)

    # Determine pseudo-steady state populations of active state
    row = 0
    indices = -numpy.ones((Ngrains,Nisom), numpy.int)
    for r in range(Ngrains):
        for i in range(Nisom):
            if r >= Nres[i]:
                indices[r,i] = row
                row += 1
    
    # Choose the half-bandwidth
    r = int(Ngrains / 2)
    tol = 1e-12
    ratio = numpy.abs(Mcoll[0,:,r] / Mcoll[0,r,r])
    ind = [i for i,x in enumerate(ratio) if x > tol]
    halfbandwidth = max(r - min(ind), max(ind) - r) * Nisom
    bandwidth = 2 * halfbandwidth + 1
    
    # Populate active-state matrix and source vectors
    L = numpy.zeros((bandwidth,numpy.sum(Nact)), numpy.float64)
    Z = numpy.zeros((numpy.sum(Nact),Nisom+Nreac), numpy.float64)
    # Collisional terms
    for i in range(Nisom):
        for r in range(Nres[i], Ngrains):
            for s in range(max(Nres[i], r-halfbandwidth/Nisom), min(Ngrains, r+halfbandwidth/Nisom)):
                L[halfbandwidth + indices[r,i] - indices[s,i], indices[s,i]] = Mcoll[i,r,s]
            Z[indices[r,i],i] = numpy.sum(Mcoll[i,r,0:Nres[i]] * eqDist[i,0:Nres[i]])
    # Isomerization terms
    for i in range(Nisom):
        for j in range(i):
            for r in range(max(Nres[i], Nres[j]), Ngrains):
                L[halfbandwidth + indices[r,j] - indices[r,i], indices[r,i]] = Kij[j,i,r]
                L[halfbandwidth, indices[r,i]] -= Kij[j,i,r]
                L[halfbandwidth + indices[r,i] - indices[r,j], indices[r,j]] = Kij[i,j,r]
                L[halfbandwidth, indices[r,j]] -= Kij[i,j,r]
    # Dissociation/association terms
    for i in range(Nisom):
        for n in range(Nreac+Nprod):
            for r in range(Nres[i], Ngrains):
                L[halfbandwidth, indices[r,i]] -= Gnj[n,i,r]
        for n in range(Nreac):
            for r in range(Nres[i], Ngrains):
                Z[indices[r,i], n+Nisom] = Fim[i,n,r] #* eqDist[n+Nisom,r]
        
    # Solve for pseudo-steady state populations of active state
    X = scipy.linalg.solve_banded((halfbandwidth,halfbandwidth), L, -Z, overwrite_ab=True, overwrite_b=True)
    pa = numpy.zeros((Ngrains,Nisom,Nisom+Nreac), numpy.float64)
    for i in range(Nisom):
        for r in range(Nres[i], Ngrains):
            for n in range(Nisom+Nreac):
                pa[r,i,n] = X[indices[r,i], n]
    
    # Double-check to ensure that we have all positive populations
    if not (pa >= 0).all():
        raise ReservoirStateError('A negative steady-state population was encountered.')

    # Put the reservoir populations into pa as well
    for i in range(Nisom):
        for r in range(Nres[i]):
            pa[r,i,i] = eqDist[i,r]

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
        # Collisional rearrangement within the reservoir of isomer i
        K[i,i] = K[i,i] + numpy.sum(numpy.dot(Mcoll[i,0:Nres[i],0:Nres[i]], eqDist[i,0:Nres[i]]))
        # Isomerization from isomer j to isomer i
        for j in range(Nisom):
            K[i,j] = K[i,j] + numpy.sum(numpy.dot(Mcoll[i,0:Nres[i],Nres[i]:Ngrains], pa[Nres[i]:Ngrains,i,j]))
        # Association from reactant n to isomer i
        for n in range(Nisom, Nisom+Nreac):
            K[i,n] = K[i,n] + numpy.sum(numpy.dot(Mcoll[i,0:Nres[i],Nres[i]:Ngrains], pa[Nres[i]:Ngrains,i,n]))
    # Rows relating to reactants
    for n in range(Nreac):
        # Association loss
        for i in range(Nisom):
            K[Nisom+n,Nisom+n] = K[Nisom+n,Nisom+n] - numpy.sum(Fim[i,n,:])# * eqDist[n+Nisom,:])
        # Reaction from isomer or reactant j to reactant n
        for j in range(Nisom+Nreac):
            for i in range(Nisom):
                K[Nisom+n,j] = K[Nisom+n,j] + numpy.sum(Gnj[n,i,Nres[i]:Ngrains] * pa[Nres[i]:Ngrains,i,j])
    # Rows relating to products
    for n in range(Nreac, Nreac+Nprod):
        # Reaction from isomer or reactant j to product n
        for j in range(Nisom+Nreac):
            for i in range(Nisom):
                K[Nisom+n,j] = K[Nisom+n,j] + numpy.sum(Gnj[n,i,Nres[i]:Ngrains] * pa[Nres[i]:Ngrains,i,j])

    # Ensure matrix is conservative
    for n in range(Nisom+Nreac):
        K[n,n] = K[n,n] - numpy.sum(K[:,n])

    # Return the matrix of k(T,P) values and the pseudo-steady population distributions
    return K, pa
