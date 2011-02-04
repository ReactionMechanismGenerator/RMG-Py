################################################################################
#
#   ChemPy - A chemistry toolkit for Python
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

cimport numpy

cdef class Mode:

    cpdef numpy.ndarray getPartitionFunctions(self, numpy.ndarray Tlist)

    cpdef numpy.ndarray getHeatCapacities(self, numpy.ndarray Tlist)

    cpdef numpy.ndarray getEnthalpies(self, numpy.ndarray Tlist)

    cpdef numpy.ndarray getEntropies(self, numpy.ndarray Tlist)

################################################################################

cdef class Translation(Mode):
    
    cdef public double mass
    
    cpdef double getPartitionFunction(self, double T)
    
    cpdef double getHeatCapacity(self, double T)
    
    cpdef double getEnthalpy(self, double T)
    
    cpdef double getEntropy(self, double T)
    
    cpdef numpy.ndarray getDensityOfStates(self, numpy.ndarray Elist)

################################################################################

cdef class RigidRotor(Mode):
    
    cdef public list inertia
    cdef public bint linear
    cdef public int symmetry

    cpdef double getPartitionFunction(self, double T)

    cpdef double getHeatCapacity(self, double T)

    cpdef double getEnthalpy(self, double T)

    cpdef double getEntropy(self, double T)

    cpdef numpy.ndarray getDensityOfStates(self, numpy.ndarray Elist)

################################################################################

cdef class HinderedRotor(Mode):
    
    cdef public double inertia
    cdef public double barrier
    cdef public int symmetry
    cdef public numpy.ndarray fourier
    cdef numpy.ndarray energies

    cpdef double getPartitionFunction(self, double T)

    cpdef double getHeatCapacity(self, double T)

    cpdef double getEnthalpy(self, double T)

    cpdef double getEntropy(self, double T)

    cpdef numpy.ndarray getDensityOfStates(self, numpy.ndarray Elist)

    cpdef numpy.ndarray getPotential(self, numpy.ndarray phi)

    cpdef double getFrequency(self)

cdef double besseli0(double x)
cdef double besseli1(double x)
cdef double cellipk(double x)

################################################################################

cdef class HarmonicOscillator(Mode):
    
    cdef public list frequencies
    
    cpdef double getPartitionFunction(self, double T)

    cpdef double getHeatCapacity(self, double T)

    cpdef double getEnthalpy(self, double T)

    cpdef double getEntropy(self, double T)

    cpdef numpy.ndarray getDensityOfStates(self, numpy.ndarray Elist, numpy.ndarray rho0=?)

################################################################################

cdef class StatesModel:
    
    cdef public list modes
    cdef public int spinMultiplicity

    cpdef double getPartitionFunction(self, double T)

    cpdef double getHeatCapacity(self, double T)

    cpdef double getEnthalpy(self, double T)

    cpdef double getEntropy(self, double T)

    cpdef numpy.ndarray getDensityOfStates(self, numpy.ndarray Elist)

    cpdef numpy.ndarray getSumOfStates(self, numpy.ndarray Elist)

    cpdef numpy.ndarray getDensityOfStatesILT(self, numpy.ndarray Elist, int order=?)

    cpdef numpy.ndarray getPartitionFunctions(self, numpy.ndarray Tlist)

    cpdef numpy.ndarray getHeatCapacities(self, numpy.ndarray Tlist)

    cpdef numpy.ndarray getEnthalpies(self, numpy.ndarray Tlist)

    cpdef numpy.ndarray getEntropies(self, numpy.ndarray Tlist)
    
################################################################################

cpdef numpy.ndarray convolve(numpy.ndarray rho1, numpy.ndarray rho2, numpy.ndarray Elist)

