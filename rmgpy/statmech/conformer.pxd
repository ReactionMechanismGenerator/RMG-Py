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

from rmgpy.quantity cimport ScalarQuantity, ArrayQuantity

################################################################################

cdef class Conformer:

    cdef public ScalarQuantity _E0
    cdef public list modes
    cdef public int spinMultiplicity
    cdef public int opticalIsomers
    cdef public ArrayQuantity _number
    cdef public ArrayQuantity _mass
    cdef public ArrayQuantity _coordinates

    cpdef double getPartitionFunction(self, double T) except -1

    cpdef double getHeatCapacity(self, double T) except -100000000

    cpdef double getEnthalpy(self, double T) except 100000000

    cpdef double getEntropy(self, double T) except -100000000

    cpdef double getFreeEnergy(self, double T) except 100000000

    cpdef numpy.ndarray getSumOfStates(self, numpy.ndarray Elist)

    cpdef numpy.ndarray getDensityOfStates(self, numpy.ndarray Elist)

    cpdef double getTotalMass(self, atoms=?) except -1

    cpdef numpy.ndarray getCenterOfMass(self, atoms=?)

    cpdef numpy.ndarray getMomentOfInertiaTensor(self)

    cpdef getPrincipalMomentsOfInertia(self)

    cpdef double getInternalReducedMomentOfInertia(self, pivots, top1) except -1

    cpdef getSymmetricTopRotors(self)

    cpdef list getActiveModes(self, bint activeJRotor=?, bint activeKRotor=?)
    
    cpdef getNumberDegreesOfFreedom(self)

################################################################################

cpdef double phi(double beta, int k, double E, logQ) except -10000000
