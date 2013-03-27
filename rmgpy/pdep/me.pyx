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
Contains functionality for generating the master equation matrix for a given
pressure-dependent reaction network.
"""

import numpy
cimport numpy

from libc.math cimport exp

import rmgpy.constants as constants

################################################################################

cpdef generateFullMEMatrix(network, bint products=True):
    """
    Generate the full master equation matrix for the network.
    """
    
    cdef numpy.ndarray[numpy.int_t,ndim=1] Jlist
    cdef numpy.ndarray[numpy.int_t,ndim=3] indices
    cdef numpy.ndarray[numpy.float64_t,ndim=1] Elist
    cdef numpy.ndarray[numpy.float64_t,ndim=2] M
    cdef numpy.ndarray[numpy.float64_t,ndim=3] densStates
    cdef numpy.ndarray[numpy.float64_t,ndim=4] Kij, Gnj, Fim
    cdef numpy.ndarray[numpy.float64_t,ndim=5] Mcoll
    cdef double T, P, beta, val
    cdef int Nisom, Nreac, Nprod, Ngrains, NJ
    cdef int i, n, r, s, u, v

    T = network.T
    P = network.P
    Elist = network.Elist
    Jlist = network.Jlist
    densStates = network.densStates
    Mcoll = network.Mcoll
    Kij = network.Kij
    Fim = network.Fim
    Gnj = network.Gnj
    Nisom = network.Nisom
    Nreac = network.Nreac
    Nprod = network.Nprod
    Ngrains = network.Ngrains
    NJ = network.NJ
    
    beta = 1. / (constants.R * T)
    
    # Construct accounting matrix
    indices = -numpy.ones((Nisom,Ngrains,NJ), numpy.int)
    Nrows = 0
    for r in range(Ngrains):
        for s in range(NJ):
            for i in range(Nisom):
                if densStates[i,r,s] > 0:
                    indices[i,r,s] = Nrows
                    Nrows += 1
    Nrows += Nreac
    if products:
        Nrows += Nprod
    
    # Construct full ME matrix
    M = numpy.zeros([Nrows,Nrows], numpy.float64)
    
    # Collision terms
    for i in range(Nisom):
        for r in range(Ngrains):
            for s in range(NJ):
                if indices[i,r,s] > -1:
                    for u in range(r, Ngrains):
                        for v in range(s, NJ):
                            M[indices[i,r,s], indices[i,u,v]] = Mcoll[i,r,s,u,v]
                            M[indices[i,u,v], indices[i,r,s]] = Mcoll[i,u,v,r,s]
    
    # Isomerization terms
    for i in range(Nisom):
        for j in range(i):
            if Kij[i,j,Ngrains-1,0] > 0 or Kij[j,i,Ngrains-1,0] > 0:
                for r in range(Ngrains):
                    for s in range(NJ):
                        u = indices[i,r,s]; v = indices[j,r,s]
                        if u > -1 and v > -1:
                            M[v,u] = Kij[j,i,r,s]
                            M[u,u] -= Kij[j,i,r,s]
                            M[u,v] = Kij[i,j,r,s]
                            M[v,v] -= Kij[i,j,r,s]
    
    # Association/dissociation terms
    for i in range(Nisom):
        for n in range(Nreac+Nprod):
            if Gnj[n,i,Ngrains-1,0] > 0:
                for r in range(Ngrains):
                    for s in range(NJ):
                        u = indices[i,r,s]
                        if products: 
                            v = Nrows - Nreac - Nprod + n
                        else:
                            v = Nrows - Nreac + n
                        if u > -1:
                            M[u,u] -= Gnj[n,i,r,s]
                            if n < Nreac or products:
                                M[v,u] = Gnj[n,i,r,s]
                            if n < Nreac:
                                val = Fim[i,n,r,s] * densStates[n+Nisom,r,s] * (2*Jlist[s]+1) * exp(-Elist[r] * beta)
                                M[u,v] = val
                                M[v,v] -= val

    return M, indices
