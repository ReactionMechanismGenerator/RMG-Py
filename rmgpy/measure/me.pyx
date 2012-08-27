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
This module provides functions for working with the full master equation,
including a method for extracting phenomenological rate coefficients
:math:`k(T,P)` using direct simulation and branching ratios.
"""

cdef extern from "math.h":
    cdef double exp(double x)

import numpy
cimport numpy
import cython

cimport rmgpy.constants as constants

from simulate import solveFullME

################################################################################

cpdef computeRateCoefficients(
    numpy.ndarray[numpy.float64_t,ndim=3] Mcoll,
    numpy.ndarray[numpy.float64_t,ndim=3] Kij,
    numpy.ndarray[numpy.float64_t,ndim=3] Fim,
    numpy.ndarray[numpy.float64_t,ndim=3] Gnj,
    numpy.ndarray[numpy.float64_t,ndim=3] p0,
    int Nisom, int Nreac, int Nprod):
    """
    Using the time-independent population vectors `p0` returned from one of
    the approximate methods and the components of the full master equation
    matrix (`Mcoll`, `Kij`, `Fim`, and `Gnj`), compute and return the
    corresponding phenomenological rate coefficients :math:`k(T,P)`. This
    should produce rates that are exactly the same as those returned by
    each method (which provides a way to check that the method is working
    as expected).
    """

    cdef int i, j, m, n, r, Ngrains
    cdef numpy.ndarray[numpy.float64_t,ndim=2] K
    
    Ngrains = p0.shape[0]
    K = numpy.zeros((Nisom+Nreac+Nprod,Nisom+Nreac+Nprod), numpy.float64)

    for i in range(Nisom):
        # Isomerization
        for j in range(Nisom):
            if i != j:
                K[i,j] = 0.0  # numpy.sum(numpy.dot(Mcoll[i,:,:], p0[:,i,j]))
                for r in range(Ngrains):
                    for l in range(Nisom):
                        K[i,j] -= Kij[l,i,r] * p0[r,i,j]
                    for n in range(Nreac+Nprod):
                        K[i,j] -= Gnj[n,i,r] * p0[r,i,j]
                    for l in range(Nisom):
                        K[i,j] += Kij[i,l,r] * p0[r,l,j]
        # Association
        for m in range(Nisom, Nisom+Nreac):
            K[i,m] = 0.0  # numpy.sum(numpy.dot(Mcoll[i,:,:], p0[:,i,m]))
            for r in range(Ngrains):
                for l in range(Nisom):
                    K[i,m] -= Kij[l,i,r] * p0[r,i,m]
                for n in range(Nreac+Nprod):
                    K[i,m] -= Gnj[n,i,r] * p0[r,i,m]
                for l in range(Nisom):
                    K[i,m] += Kij[i,l,r] * p0[r,l,m]
                K[i,m] += Fim[i,m-Nisom,r]

    for n in range(Nisom, Nisom+Nreac+Nprod):
        # Dissociation
        for j in range(Nisom):
            for r in range(Ngrains):
                for l in range(Nisom):
                    K[n,j] += Gnj[n-Nisom,l,r] * p0[r,l,j]
        # Bimolecular
        for m in range(Nisom, Nisom+Nreac):
            if m != n:
                for r in range(Ngrains):
                    for l in range(Nisom):
                        K[n,m] += Gnj[n-Nisom,l,r] * p0[r,l,m]

    for i in range(Nisom+Nreac):
        K[i,i] = -numpy.sum(K[:,i])

    return K

################################################################################

@cython.boundscheck(True)
def applyBranchingRatiosMethod(double T, double P,
    numpy.ndarray[numpy.float64_t,ndim=1] Elist,
    numpy.ndarray[numpy.float64_t,ndim=2] densStates,
    numpy.ndarray[numpy.float64_t,ndim=3] Mcoll,
    numpy.ndarray[numpy.float64_t,ndim=3] Kij,
    numpy.ndarray[numpy.float64_t,ndim=3] Fim,
    numpy.ndarray[numpy.float64_t,ndim=3] Gnj,
    numpy.ndarray[numpy.float64_t,ndim=1] Ereac,
    int Nisom, int Nreac, int Nprod):
    """
    Use the branching ratios method to reduce the master equation model
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

    cdef int Ngrains, Nrows, index, i, n, r, s
    cdef numpy.ndarray[numpy.int_t,ndim=1] Nres
    cdef numpy.ndarray[numpy.int_t,ndim=2] indices
    cdef numpy.ndarray[numpy.float64_t,ndim=1] x0, t
    cdef numpy.ndarray[numpy.float64_t,ndim=2] M0, M1, M, K, x, eqDist
    cdef numpy.ndarray[numpy.float64_t,ndim=3] pa, p
    cdef double ktot

    Ngrains = len(Elist)

    K = numpy.zeros((Nisom+Nreac+Nprod, Nisom+Nreac+Nprod), numpy.float64)
    pa = numpy.zeros((Ngrains,Nisom,Nisom+Nreac), numpy.float64)

    # Determine the reservoir cutoff grain for each isomer
    # Start by simply placing it at the lowest reactive grain
    Nres = numpy.zeros(Nisom, numpy.int)
    for i in range(Nisom):
        for r in range(Ngrains):
            if Elist[r] > Ereac[i] and Nres[i] == 0:
                Nres[i] = r
                break

    # Determine equilibrium distributions
    eqDist = numpy.zeros((Nisom+Nreac,Ngrains), numpy.float64)
    for i in range(Nisom+Nreac):
        eqDist[i,:] = densStates[i,:] * numpy.exp(-Elist / constants.R / T)

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
    Nrows += Nreac + Nprod

    M0 = generateFullMEMatrix(Mcoll, Kij, Fim, Gnj, Ngrains, Nisom, Nreac, Nprod, indices, products=True)

    tlist = numpy.array([10.**n for n in range(-15,1)], numpy.float64)

    # Thermal activation: isomer as source term
    for i in range(Nisom):
        M = M0.copy()
        # Forbid collisions causing reactivation from nonreactive grain to reactive grain
        # except for source isomer
        for j in range(Nisom):
            if i != j:
                for r in range(Ngrains):
                    for s in range(Ngrains):
                        if Ereac[j] < Elist[r] and Ereac[j] > Elist[s] and indices[r,j] > -1 and indices[s,j] > -1 and r != s:
                            M[indices[s,j],indices[s,j]] += M[indices[r,j],indices[s,j]]
                            M[indices[r,j],indices[s,j]] = 0.0
        # Make deactivation of source isomer from reactive grain to nonreactive grain irreversible
        for r in range(Ngrains):
            for s in range(Ngrains):
                if Ereac[i] > Elist[r] and Ereac[i] < Elist[s] and indices[r,i] > -1 and indices[s,i] > -1 and r != s:
                    M[indices[r,i],indices[s,i]] = 0.0
        # Forbid reassociation of reactant channels
        for n in range(Nreac):
            M[:,Nrows-Nreac-Nprod+n] = numpy.zeros((Nrows), numpy.float64)

        x0 = numpy.zeros((Nisom+Nreac+Nprod), numpy.float64)
        x0[i] = 1.0
        t, p, x = solveFullME(T, P, Elist, tlist, x0, M, indices, densStates, Nisom, Nreac, Nprod)
        
        # Get total rate of thermal activation from this isomer
        # Not sure that this works under all conditions
        ktot = numpy.sum(numpy.dot(Mcoll[i,Nres[i]:Ngrains,0:Nres[i]], eqDist[i,0:Nres[i]]))
        
        # Multiply total rate of activation by branching ratios to get k(T,P) values
        for j in range(Nisom+Nreac+Nprod):
            if i != j:
                K[j,i] = ktot * x[-1,j]
        K[i,i] = -numpy.sum(K[:,i])

    # Chemical activation: reactant channel as source term
    M1 = M0.copy()
    # Forbid collisions causing reactivation from nonreactive grain to reactive grain
    for i in range(Nisom):
        for r in range(Ngrains):
            for s in range(Ngrains):
                if Ereac[i] < Elist[r] and Ereac[i] > Elist[s] and indices[r,i] > -1 and indices[s,i] > -1 and r != s:
                    M1[indices[s,i],indices[s,i]] += M1[indices[r,i],indices[s,i]]
                    M1[indices[r,i],indices[s,i]] = 0.0
    
    for n in range(Nreac):
        M = M1.copy()
        # Make dissociation back to reactants irreversible by removing dissociation terms from reactant channel row
        M[Nrows-Nreac-Nprod+n,0:Nrows-Nreac-Nprod] = numpy.zeros((1,Nrows-Nreac-Nprod), numpy.float64)

        x0 = numpy.zeros((Nisom+Nreac+Nprod), numpy.float64)
        x0[Nisom+n] = 1.0
        t, p, x = solveFullME(T, P, Elist, tlist, x0, M, indices, densStates, Nisom, Nreac, Nprod)

        # Get total rate of chemical activation from this reactant channel to
        # any isomer
        ktot = 0.0
        for i in range(Nisom):
            if Gnj[n,i,-1] > 0 or Fim[i,n,-1] > 0:
                ktot += numpy.sum(Fim[i,n,:])
                
        # Multiply total rate of activation by branching ratios to get k(T,P) values
        for i in range(Nisom+Nreac+Nprod):
            if i != n+Nisom:
                K[i,n+Nisom] = ktot * x[-1,i]
        K[n+Nisom,n+Nisom] = -numpy.sum(K[:,n+Nisom])

    # Return the matrix of k(T,P) values and the pseudo-steady population distributions
    return K, pa

################################################################################

def generateFullMEMatrix(
    numpy.ndarray[numpy.float64_t,ndim=3] Mcoll,
    numpy.ndarray[numpy.float64_t,ndim=3] Kij,
    numpy.ndarray[numpy.float64_t,ndim=3] Fim,
    numpy.ndarray[numpy.float64_t,ndim=3] Gnj,
    int Ngrains, int Nisom, int Nreac, int Nprod,
    numpy.ndarray[numpy.int_t,ndim=2] indices,
    bint products=False,
):

    cdef int Nrows, i, j, n, r, s, u, v
    cdef numpy.ndarray[numpy.float64_t,ndim=2] M

    Nrows = numpy.max(indices) + 1 + Nreac

    if products: Nrows += Nprod
    M = numpy.zeros((Nrows, Nrows), numpy.float64)
    if products: Nrows -= Nprod

    # Add collision terms
    for i in range(Nisom):
        for r in range(Ngrains):
            if indices[r,i] > -1:
                for s in range(Ngrains):
                    if indices[s,i] > -1:
                        M[indices[r,i], indices[s,i]] = Mcoll[i,r,s]

    # Add isomerization terms
    for i in range(Nisom):
        for j in range(i):
            if Kij[i,j,-1] > 0 or Kij[j,i,-1] > 0:
                for r in range(Ngrains):
                    u = indices[r,i]; v = indices[r,j]
                    if u > -1 and v > -1:
                        M[v,u] = Kij[j,i,r]
                        M[u,u] -= Kij[j,i,r]
                        M[u,v] = Kij[i,j,r]
                        M[v,v] -= Kij[i,j,r]

    # Add dissociation/association terms
    for i in range(Nisom):
        for n in range(Nreac):
            if Gnj[n,i,-1] > 0 or Fim[i,n,-1] > 0:
                for r in range(Ngrains):
                    u = indices[r,i]; v = Nrows - Nreac + n
                    if u > -1:
                        M[v,u] = Gnj[n,i,r]
                        M[u,u] -= Gnj[n,i,r]
                        M[u,v] = Fim[i,n,r]
                        M[v,v] -= Fim[i,n,r]
    for i in range(Nisom):
        for n in range(Nreac, Nreac+Nprod):
            if Gnj[n,i,-1] > 0:
                for r in range(Ngrains):
                    u = indices[r,i]; v = Nrows - Nreac + n
                    if u > -1:
                        M[u,u] -= Gnj[n,i,r]
                        if products: M[v,u] = Gnj[n,i,r]
                        
    return M
