################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2009-2011 by the RMG Team (rmg_dev@mit.edu)
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

cimport numpy

from rmgpy.statmech.mode cimport Mode
from rmgpy.quantity cimport ScalarQuantity, ArrayQuantity

################################################################################

cdef class Vibration(Mode):
    pass

################################################################################

cdef class HarmonicOscillator(Vibration):

    cdef public ArrayQuantity _frequencies

    cpdef double getPartitionFunction(self, double T) except -1

    cpdef double getHeatCapacity(self, double T) except -100000000

    cpdef double getEnthalpy(self, double T) except 100000000

    cpdef double getEntropy(self, double T) except -100000000

    cpdef numpy.ndarray getSumOfStates(self, numpy.ndarray Elist, numpy.ndarray sumStates0=?)
    
    cpdef numpy.ndarray getDensityOfStates(self, numpy.ndarray Elist, numpy.ndarray densStates0=?)
