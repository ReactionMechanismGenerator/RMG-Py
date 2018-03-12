###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2018 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

cimport numpy

from rmgpy.statmech.mode cimport Mode
from rmgpy.quantity cimport ScalarQuantity, ArrayQuantity

################################################################################

cdef class Rotation(Mode):

    cdef public int symmetry

################################################################################

cdef class LinearRotor(Rotation):

    cdef public ScalarQuantity _inertia

    cpdef double getLevelEnergy(self, int J) except -1
    
    cpdef int getLevelDegeneracy(self, int J) except -1
    
    cpdef double getPartitionFunction(self, double T) except -1
        
    cpdef double getHeatCapacity(self, double T) except -100000000

    cpdef double getEnthalpy(self, double T) except 100000000

    cpdef double getEntropy(self, double T) except -100000000

    cpdef numpy.ndarray getSumOfStates(self, numpy.ndarray Elist, numpy.ndarray sumStates0=?)
    
    cpdef numpy.ndarray getDensityOfStates(self, numpy.ndarray Elist, numpy.ndarray densStates0=?)

################################################################################

cdef class NonlinearRotor(Rotation):

    cdef public ArrayQuantity _inertia

    cdef numpy.ndarray getRotationalConstantEnergy(self)

    cpdef double getPartitionFunction(self, double T) except -1

    cpdef double getHeatCapacity(self, double T) except -100000000

    cpdef double getEnthalpy(self, double T) except 100000000

    cpdef double getEntropy(self, double T) except -100000000

    cpdef numpy.ndarray getSumOfStates(self, numpy.ndarray Elist, numpy.ndarray sumStates0=?)

    cpdef numpy.ndarray getDensityOfStates(self, numpy.ndarray Elist, numpy.ndarray densStates0=?)

################################################################################

cdef class KRotor(Rotation):

    cdef public ScalarQuantity _inertia

    cpdef double getLevelEnergy(self, int J) except -1

    cpdef int getLevelDegeneracy(self, int J) except -1

    cpdef double getPartitionFunction(self, double T) except -1

    cpdef double getHeatCapacity(self, double T) except -100000000

    cpdef double getEnthalpy(self, double T) except 100000000

    cpdef double getEntropy(self, double T) except -100000000

    cpdef numpy.ndarray getSumOfStates(self, numpy.ndarray Elist, numpy.ndarray sumStates0=?)

    cpdef numpy.ndarray getDensityOfStates(self, numpy.ndarray Elist, numpy.ndarray densStates0=?)

################################################################################

cdef class SphericalTopRotor(Rotation):

    cdef public ScalarQuantity _inertia

    cpdef double getLevelEnergy(self, int J) except -1
    
    cpdef int getLevelDegeneracy(self, int J) except -1
    
    cpdef double getPartitionFunction(self, double T) except -1
    
    cpdef double getHeatCapacity(self, double T) except -100000000

    cpdef double getEnthalpy(self, double T) except 100000000

    cpdef double getEntropy(self, double T) except -100000000

    cpdef numpy.ndarray getSumOfStates(self, numpy.ndarray Elist, numpy.ndarray sumStates0=?)
    
    cpdef numpy.ndarray getDensityOfStates(self, numpy.ndarray Elist, numpy.ndarray densStates0=?)
