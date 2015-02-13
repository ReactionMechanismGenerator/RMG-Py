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

from rmgpy.thermo.model cimport HeatCapacityModel
from rmgpy.thermo.thermodata cimport ThermoData
from rmgpy.thermo.nasa cimport NASA, NASAPolynomial
from rmgpy.quantity cimport ScalarQuantity, ArrayQuantity

################################################################################

cdef class Wilhoit(HeatCapacityModel):
    
    cdef public ScalarQuantity _Cp0, _CpInf, _B, _H0, _S0
    cdef public double a0, a1, a2, a3
    
    cpdef double getHeatCapacity(self, double T) except -1000000000

    cpdef double getEnthalpy(self, double T) except 1000000000

    cpdef double getEntropy(self, double T) except -1000000000

    cpdef double getFreeEnergy(self, double T) except 1000000000
    
    cpdef Wilhoit copy(self)
    
    cdef double integral_T0(self, double T)
    
    cdef double integral_TM1(self, double T)
    
    cdef double integral_T1(self, double T)
    
    cdef double integral_T2(self, double T)
    
    cdef double integral_T3(self, double T)
    
    cdef double integral_T4(self, double T)
    
    cdef double integral2_T0(self, double T)
    
    cdef double integral2_TM1(self, double T)
    
    cpdef ThermoData toThermoData(self)

    cpdef NASA toNASA(self, double Tmin, double Tmax, double Tint, bint fixedTint=?, bint weighting=?, int continuity=?)

################################################################################

cpdef Wilhoit_to_NASA(Wilhoit wilhoit, double Tmin, double Tmax, double Tint, bint weighting, int contCons)

cpdef Wilhoit_to_NASA_TintOpt(Wilhoit wilhoit, double Tmin, double Tmax, bint weighting, int contCons)

cpdef double Wilhoit_to_NASA_TintOpt_objFun(double Tint, Wilhoit wilhoit, double Tmin, double Tmax, bint weighting, int contCons)

cpdef double Wilhoit_to_NASA_TintOpt_objFun_NW(double Tint, Wilhoit wilhoit, double Tmin, double Tmax, int contCons)

cpdef double Wilhoit_to_NASA_TintOpt_objFun_W(double Tint, Wilhoit wilhoit, double Tmin, double Tmax, int contCons)
