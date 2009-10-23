################################################################################
#
#	RMG - Reaction Mechanism Generator
#
#	Copyright (c) 2002-2009 Prof. William H. Green (whgreen@mit.edu) and the
#	RMG Team (rmg_dev@mit.edu)
#
#	Permission is hereby granted, free of charge, to any person obtaining a
#	copy of this software and associated documentation files (the 'Software'),
#	to deal in the Software without restriction, including without limitation
#	the rights to use, copy, modify, merge, publish, distribute, sublicense,
#	and/or sell copies of the Software, and to permit persons to whom the
#	Software is furnished to do so, subject to the following conditions:
#
#	The above copyright notice and this permission notice shall be included in
#	all copies or substantial portions of the Software.
#
#	THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#	FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#	DEALINGS IN THE SOFTWARE.
#
################################################################################

cdef extern from "dictobject.h":
	ctypedef class __builtin__.dict [object PyDictObject]:
		pass

################################################################################

cdef class ThermoData:
	
	cdef public object Trange # to be removed
	cdef public str comment
	cdef public float Tmin
	cdef public float Tmax
	
	cpdef bint isTemperatureValid(ThermoData self, float T)

################################################################################

cdef class ThermoNASAPolynomial(ThermoData):
	
	cdef public float c0, c1, c2, c3, c4, c5, c6 
	
	cpdef float getHeatCapacity(ThermoNASAPolynomial self, float T)
	
	cpdef float getEnthalpy(ThermoNASAPolynomial self, float T)
	
	cpdef float getEntropy(ThermoNASAPolynomial self, float T)
	
	cpdef float getFreeEnergy(ThermoNASAPolynomial self, float T)	

################################################################################

