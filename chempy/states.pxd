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
	pass

################################################################################

cdef class Translation(Mode):
	
	cdef public double mass
	cdef public double volume
	cdef public int dimension
	
	cpdef getPartitionFunction(self, numpy.ndarray Tlist)
	
	cpdef getHeatCapacity(self, numpy.ndarray Tlist)
	
	cpdef getEnthalpy(self, numpy.ndarray Tlist)
	
	cpdef getEntropy(self, numpy.ndarray Tlist)
	
	cpdef getDensityOfStates(self, numpy.ndarray Elist)

################################################################################

cdef class RigidRotor(Mode):
	
	cdef public list frequencies
	cdef public bint linear
	
	cpdef getPartitionFunction(RigidRotor self, Tlist)
	
	cpdef getHeatCapacity(RigidRotor self, Tlist)
	
	cpdef getDensityOfStates(RigidRotor self, Elist)

################################################################################

cdef class HinderedRotor(Mode):
	
	cdef public double frequency
	cdef public double barrier
	
	cpdef getPartitionFunction(HinderedRotor self, Tlist)
	
	cpdef getHeatCapacity(HinderedRotor self, Tlist)
	
	cpdef getDensityOfStates(HinderedRotor self, Elist)

cdef double besseli0(double x)
cdef double besseli1(double x)
cdef double cellipk(double x)

################################################################################

cdef class HarmonicOscillator(Mode):
	
	cdef public double frequency
	
	cpdef getPartitionFunction(HarmonicOscillator self, Tlist)
	
	cpdef getHeatCapacity(HarmonicOscillator self, Tlist)
	
	cpdef getDensityOfStates(HarmonicOscillator self, Elist)

################################################################################

cdef class StatesModel:
	
	cdef public list modes
	
	cpdef getHeatCapacity(StatesModel self, Tlist)
	
	cpdef getPartitionFunction(StatesModel self, Tlist)
	
	cpdef getDensityOfStates(StatesModel self, Elist, bint linear)
