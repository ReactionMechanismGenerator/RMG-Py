#!/usr/bin/python
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
Contains the :class:`Network` class, which defines an internal representation
of a unimolecular reaction network. This provides a convienent means of
keeping track of things while we are performing the master equation calculation.
"""

import cython
import numpy
cimport numpy

################################################################################

@cython.boundscheck(False)
def mapDensitiesOfStates(numpy.ndarray[numpy.float64_t,ndim=1] Elist,
    numpy.ndarray[numpy.float64_t,ndim=1] Elist0,
    numpy.ndarray[numpy.float64_t,ndim=2] densStates0,
    double T,
    numpy.ndarray[numpy.float64_t,ndim=1] E0,
    int Nisom, int Nreac, int Nprod):
    """
    Map the overall densities of states to the current energy grains.
    Semi-logarithmic interpolation will be used if the grain sizes of 
    `Elist0` and `Elist` do not match; this should not be a significant
    source of error as long as the grain sizes are sufficiently small.
    """
    
    cdef numpy.ndarray[numpy.float64_t,ndim=2] densStates
    cdef int i, r, s, Ngrains, Ngrains0
    
    Ngrains = Elist.shape[0]
    Ngrains0 = Elist0.shape[0]
    
    densStates = numpy.zeros((Nisom+Nreac+Nprod, Ngrains), numpy.float64)
    for i in range(Nisom+Nreac+Nprod):
        for r in range(Ngrains):
            if Elist[r] >= E0[i]:
                for s in range(Ngrains0):
                    if E0[i] + Elist0[s] >= Elist[r]:
                        if densStates0[i,s-1] > 0 and densStates0[i,s] > 0:
                            densStates[i,r] = densStates0[i,s] * (densStates0[i,s-1] / densStates0[i,s]) ** ((Elist[r] - E0[i] - Elist0[s]) / (Elist0[s-1] - Elist0[s]))
                        else:
                            densStates[i,r] = densStates0[i,s] + (densStates0[i,s-1] - densStates0[i,s]) * (Elist[r] - E0[i] - Elist0[s]) / (Elist0[s-1] - Elist0[s])
                        break
                else:
                    if i < Nisom:
                        raise Exception('Unable to determine density of states for isomer {0} at {1:g} K.'.format(i,T))
                    elif i < Nisom+Nreac:
                        raise Exception('Unable to determine density of states for reactant channel {0} at {1:g} K.'.format(i-Nisom,T))
                    else:
                        raise Exception('Unable to determine density of states for product channel {0} at {1:g} K.'.format(i-Nisom-Nreac,T))
    
    return densStates
