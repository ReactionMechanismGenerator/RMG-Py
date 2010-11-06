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

################################################################################

cdef class ThermoModel:
    
    cdef public double Tmin
    cdef public double Tmax
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

cdef class ThermoGAModel(ThermoModel):

    cdef public numpy.ndarray Tdata, Cpdata, dCp
    cdef public double H298, S298, dH, dS
    
    cpdef double getHeatCapacity(self, double T)

    cpdef double getEnthalpy(self, double T)

    cpdef double getEntropy(self, double T)

    cpdef double getFreeEnergy(self, double T)

################################################################################

cdef class WilhoitModel(ThermoModel):
    
    cdef public double cp0
    cdef public double cpInf
    cdef public double B
    cdef public double a0
    cdef public double a1
    cdef public double a2
    cdef public double a3
    cdef public double H0
    cdef public double S0
    
    cpdef double getHeatCapacity(self, double T)

    cpdef double getEnthalpy(self, double T)

    cpdef double getEntropy(self, double T)

    cpdef double getFreeEnergy(self, double T)
    
    cpdef double __residual(self, double B, numpy.ndarray Tlist, numpy.ndarray Cplist, 
        bint linear, int nFreq, int nRotors, double H298, double S298)
    
    cpdef WilhoitModel fitToData(self, numpy.ndarray Tlist, numpy.ndarray Cplist,
        bint linear, int nFreq, int nRotors, double H298, double S298, double B0=?)
    
    cpdef WilhoitModel fitToDataForConstantB(self, numpy.ndarray Tlist, numpy.ndarray Cplist,
        bint linear, int nFreq, int nRotors, double B, double H298, double S298)
    
################################################################################

cdef class NASAPolynomial(ThermoModel):
    
    cdef public double cm2, cm1, c0, c1, c2, c3, c4, c5, c6
    
    cpdef double getHeatCapacity(self, double T)

    cpdef double getEnthalpy(self, double T)

    cpdef double getEntropy(self, double T)

    cpdef double getFreeEnergy(self, double T)
    
################################################################################

cdef class NASAModel(ThermoModel):
    
    cdef public list polynomials
    
    cpdef double getHeatCapacity(self, double T)

    cpdef double getEnthalpy(self, double T)

    cpdef double getEntropy(self, double T)

    cpdef double getFreeEnergy(self, double T)
    
    cpdef NASAPolynomial __selectPolynomialForTemperature(self, double T)
