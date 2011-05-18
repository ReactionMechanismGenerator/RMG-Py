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

from constants cimport Quantity

################################################################################

cdef class ThermoModel:
    
    cdef public Quantity Tmin, Tmax
    cdef public str comment
    
    cpdef bint isTemperatureValid(ThermoModel self, double T) except -2

#    cpdef double getHeatCapacity(self, double T)
#
#    cpdef double getEnthalpy(self, double T)
#
#    cpdef double getEntropy(self, double T)
#
#    cpdef double getFreeEnergy(self, double T)

    cpdef numpy.ndarray getHeatCapacities(self, numpy.ndarray Tlist)

    cpdef numpy.ndarray getEnthalpies(self, numpy.ndarray Tlist)

    cpdef numpy.ndarray getEntropies(self, numpy.ndarray Tlist)

    cpdef numpy.ndarray getFreeEnergies(self, numpy.ndarray Tlist)
    
################################################################################

cdef class ThermoData(ThermoModel):

    cdef public Quantity Tdata, Cpdata, H298, S298
    
    cpdef double getHeatCapacity(self, double T)

    cpdef double getEnthalpy(self, double T)

    cpdef double getEntropy(self, double T)

    cpdef double getFreeEnergy(self, double T)

################################################################################

cdef class Wilhoit(ThermoModel):
    
    cdef public Quantity cp0, cpInf, B, a0, a1, a2, a3, H0, S0
    
    cpdef double getHeatCapacity(self, double T)

    cpdef double getEnthalpy(self, double T)

    cpdef double getEntropy(self, double T)

    cpdef double getFreeEnergy(self, double T)
    
    cpdef double __residual(self, double B, numpy.ndarray Tlist, numpy.ndarray Cplist, 
        bint linear, int nFreq, int nRotors, double H298, double S298)
    
    cpdef Wilhoit fitToData(self, numpy.ndarray Tlist, numpy.ndarray Cplist,
        bint linear, int nFreq, int nRotors, double H298, double S298, double B0=?)
    
    cpdef Wilhoit fitToDataForConstantB(self, numpy.ndarray Tlist, numpy.ndarray Cplist,
        bint linear, int nFreq, int nRotors, double B, double H298, double S298)
    
################################################################################

cdef class NASA(ThermoModel):
    
    cdef public double cm2, cm1, c0, c1, c2, c3, c4, c5, c6
    
    cpdef double getHeatCapacity(self, double T)

    cpdef double getEnthalpy(self, double T)

    cpdef double getEntropy(self, double T)

    cpdef double getFreeEnergy(self, double T)
    
################################################################################

cdef class MultiNASA(ThermoModel):
    
    cdef public list polynomials
    
    cpdef double getHeatCapacity(self, double T)

    cpdef double getEnthalpy(self, double T)

    cpdef double getEntropy(self, double T)

    cpdef double getFreeEnergy(self, double T)
    
    cpdef NASA __selectPolynomialForTemperature(self, double T)
