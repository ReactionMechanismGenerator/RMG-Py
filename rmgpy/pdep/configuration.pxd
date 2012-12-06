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

cimport numpy

################################################################################

cdef class Configuration:

    cdef public list species
    cdef public numpy.ndarray Elist
    cdef public numpy.ndarray densStates
    cdef public numpy.ndarray sumStates
    cdef public bint activeJRotor
    cdef public bint activeKRotor

    cpdef cleanup(self)

    cpdef bint isUnimolecular(self) except -2
    
    cpdef bint isBimolecular(self) except -2
    
    cpdef bint isTransitionState(self) except -2
    
    cpdef bint hasStatMech(self) except -2
    
    cpdef bint hasThermo(self) except -2
    
    cpdef double getHeatCapacity(self, double T) except -100000000

    cpdef double getEnthalpy(self, double T) except 100000000

    cpdef double getEntropy(self, double T) except -100000000

    cpdef double getFreeEnergy(self, double T) except 100000000
    
    cpdef double calculateCollisionFrequency(self, double T, double P, dict bathGas) except -1
        
    cpdef numpy.ndarray generateCollisionMatrix(self, double T, numpy.ndarray densStates, numpy.ndarray Elist, numpy.ndarray Jlist=?)
    
    cpdef calculateDensityOfStates(self, numpy.ndarray Elist, bint activeJRotor=?, bint activeKRotor=?, bint rmgmode=?)
