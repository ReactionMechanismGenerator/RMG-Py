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

cdef extern from "math.h":
	cdef double acos(double x)
	cdef double cos(double x)
	cdef double exp(double x)
	cdef double log(double x)
	cdef double log10(double x)
	cdef double pow(double base, double exponent)

################################################################################

cdef class KineticsModel:
	
	cdef public double Tmin
	cdef public double Tmax
	cdef public double Pmin
	cdef public double Pmax
	cdef public int numReactants
	cdef public str comment
	
	cpdef bint isTemperatureValid(KineticsModel self, double T) except -2
	
	cpdef bint isPressureValid(KineticsModel self, double P) except -2

################################################################################

cdef class ArrheniusModel(KineticsModel):
	
	cdef public double A
	cdef public double Ea
	cdef public double n
	
	cpdef double getRateConstant(ArrheniusModel self, double T)

################################################################################

cdef class ArrheniusEPModel(KineticsModel):
	
	cdef public double A
	cdef public double E0
	cdef public double n
	cdef public double alpha
	
	cpdef double getActivationEnergy(ArrheniusEPModel self, double dHrxn)
	
	cpdef double getRateConstant(ArrheniusEPModel self, double T, double dHrxn)

################################################################################

cdef class PDepArrheniusModel(KineticsModel):
	
	cdef public list pressures
	cdef public list arrhenius
	
	cpdef tuple __getAdjacentExpressions(PDepArrheniusModel self, double P)
	
	cpdef getRateConstant(PDepArrheniusModel self, double T, double P)

################################################################################

cdef class ChebyshevModel(KineticsModel):
	
	cdef public object coeffs
	cdef public int degreeT
	cdef public int degreeP
	
	cpdef double __chebyshev(ChebyshevModel self, double n, double x)
	
	cpdef double __getReducedTemperature(ChebyshevModel self, double T)
	
	cpdef double __getReducedPressure(ChebyshevModel self, double P)
	
	cpdef double getRateConstant(ChebyshevModel self, double T, double P)
	