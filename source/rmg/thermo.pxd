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


cdef extern from "math.h":
	cdef float log(float theta)
################################################################################

cdef class ThermoData:
	
	cdef public object Trange # to be removed
	cdef public str comment
	cdef public float Tmin
	cdef public float Tmax
	
	cpdef bint isTemperatureValid(ThermoData self, float T) except -2

################################################################################

cdef class ThermoGAData(ThermoData):

	cdef public float H298, S298
	cdef public list Cp
	cdef public str index
	
	# can't cpdef special methods like __add__ :-(
	#cpdef ThermoGAData __add__(ThermoGAData self, ThermoGAData other)
	
	cpdef float getHeatCapacity(ThermoGAData self, float T)
	
	cpdef float getEnthalpy(ThermoGAData self, float T)

################################################################################

cdef class ThermoWilhoitData(ThermoData):

	cdef public float cp0
	cdef public float cpInf
	cdef public float B
	cdef public float a0
	cdef public float a1
	cdef public float a2
	cdef public float a3
	cdef public float H0
	cdef public float S0

	cpdef float getHeatCapacity(ThermoWilhoitData self, float T)

	cpdef float getEnthalpy(ThermoWilhoitData self, float T)

	cpdef float getEntropy(ThermoWilhoitData self, float T)

	cpdef float getFreeEnergy(ThermoWilhoitData self, float T)

	cpdef float integral_T0(ThermoWilhoitData self, float t)

	cpdef float integral_TM1(ThermoWilhoitData self, float t)

	cpdef float integral_T1(ThermoWilhoitData self, float t)

	cpdef float integral_T2(ThermoWilhoitData self, float t)

	cpdef float integral_T3(ThermoWilhoitData self, float t)

	cpdef float integral_T4(ThermoWilhoitData self, float t)

	cpdef float integral2_T0(ThermoWilhoitData self, float t)

	cpdef float integral2_TM1(ThermoWilhoitData self, float t)

################################################################################

cdef class ThermoNASAPolynomial(ThermoData):
	
	cdef public float c0, c1, c2, c3, c4, c5, c6 
	
	cpdef float getHeatCapacity(ThermoNASAPolynomial self, float T)
	
	cpdef float getEnthalpy(ThermoNASAPolynomial self, float T)
	
	cpdef float getEntropy(ThermoNASAPolynomial self, float T)
	
	cpdef float getFreeEnergy(ThermoNASAPolynomial self, float T)
	
	cpdef float integral2_T0(ThermoNASAPolynomial self, float t)
	
	cpdef float integral2_TM1(ThermoNASAPolynomial self, float t)
	

################################################################################

