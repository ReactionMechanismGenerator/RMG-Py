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
from rmgpy.thermo.wilhoit cimport Wilhoit

################################################################################

cdef class NASAPolynomial(HeatCapacityModel):

    cdef public double cm2, cm1, c0, c1, c2, c3, c4, c5, c6
    
    cpdef double getHeatCapacity(self, double T) except -100000000

    cpdef double getEnthalpy(self, double T) except 100000000

    cpdef double getEntropy(self, double T) except -100000000

    cpdef double getFreeEnergy(self, double T) except 100000000    
    
    cpdef changeBaseEnthalpy(self, double deltaH)

    cdef double integral2_T0(self, double T)
    
    cdef double integral2_TM1(self, double T)
    
cdef class NASA(HeatCapacityModel):

    cdef public NASAPolynomial poly1, poly2, poly3
    
    cpdef NASAPolynomial selectPolynomial(self, double T)

    cpdef double getHeatCapacity(self, double T) except -100000000

    cpdef double getEnthalpy(self, double T) except 100000000

    cpdef double getEntropy(self, double T) except -100000000

    cpdef double getFreeEnergy(self, double T) except 100000000

    cpdef ThermoData toThermoData(self, double Cp0=?, double CpInf=?)

    cpdef Wilhoit toWilhoit(self, double Cp0, double CpInf)
    
    cpdef NASA changeBaseEnthalpy(self, double deltaH)
