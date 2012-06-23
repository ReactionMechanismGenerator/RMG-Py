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

from quantity cimport Quantity, constants

cdef extern from "math.h":
    double log(double)

################################################################################

cdef class ThermoModel:
    
    cdef public Quantity Tmin, Tmax
    cdef public str comment
    
    cpdef bint isTemperatureValid(ThermoModel self, double T) except -2

    cpdef double getHeatCapacity(self, double T)

    cpdef double getEnthalpy(self, double T)

    cpdef double getEntropy(self, double T)

    cpdef double getFreeEnergy(self, double T)

    cpdef numpy.ndarray getHeatCapacities(self, numpy.ndarray Tlist)

    cpdef numpy.ndarray getEnthalpies(self, numpy.ndarray Tlist)

    cpdef numpy.ndarray getEntropies(self, numpy.ndarray Tlist)

    cpdef numpy.ndarray getFreeEnergies(self, numpy.ndarray Tlist)
    
    cpdef bint isSimilarTo(self, ThermoModel other) except -2

    cpdef bint isIdenticalTo(self, ThermoModel other) except -2

################################################################################

cdef class ThermoData(ThermoModel):

    cdef public Quantity Tdata, Cpdata, H298, S298
    
    cpdef double getHeatCapacity(self, double T)

    cpdef double getEnthalpy(self, double T)

    cpdef double getEntropy(self, double T)

    cpdef double getFreeEnergy(self, double T)

################################################################################

cdef class Wilhoit(ThermoModel):
    
    cdef public Quantity cp0, cpInf, B, H0, S0
    cdef public double a0, a1, a2, a3
    
    cpdef double getHeatCapacity(self, double T)

    cpdef double getEnthalpy(self, double T)

    cpdef double getEntropy(self, double T)

    cpdef double getFreeEnergy(self, double T)
    
    cpdef double __residual(self, double B, numpy.ndarray Tlist, numpy.ndarray Cplist, 
        bint linear, int nFreq, int nRotors, double H298, double S298)
    
    cpdef Wilhoit fitToData(self, numpy.ndarray Tlist, numpy.ndarray Cplist,
        bint linear, int nFreq, int nRotors, double H298, double S298, double B0=?, double Bmin=?, double Bmax=?)
    
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

################################################################################

cpdef MultiNASA convertWilhoitToNASA(Wilhoit wilhoit, double Tmin, double Tmax, double Tint, bint fixedTint=?, bint weighting=?, int continuity=?)

cpdef Wilhoit2NASA(Wilhoit wilhoit, double tmin, double tmax, double tint, bint weighting, int contCons)

cpdef Wilhoit2NASA_TintOpt(Wilhoit wilhoit, double tmin, double tmax, bint weighting, int contCons)

cpdef TintOpt_objFun(tint, Wilhoit wilhoit, double tmin, double tmax, bint weighting, int contCons)

cpdef TintOpt_objFun_NW(double tint, Wilhoit wilhoit, double tmin, double tmax, int contCons)

cpdef TintOpt_objFun_W(double tint, Wilhoit wilhoit, double tmin, double tmax, int contCons)

cpdef double Wilhoit_integral_T0(Wilhoit wilhoit, double t)

cpdef double Wilhoit_integral_TM1(Wilhoit wilhoit, double t)

cpdef double Wilhoit_integral_T1(Wilhoit wilhoit, double t)

cpdef double Wilhoit_integral_T2(Wilhoit wilhoit, double t)

cpdef double Wilhoit_integral_T3(Wilhoit wilhoit, double t)

cpdef double Wilhoit_integral_T4(Wilhoit wilhoit, double t)

cpdef double Wilhoit_integral2_T0(Wilhoit wilhoit, double t)

cpdef double Wilhoit_integral2_TM1(Wilhoit wilhoit, double t)

cpdef double Wilhoit_Dintegral_T0(Wilhoit wilhoit, double t1, double t2)

cpdef double Wilhoit_Dintegral_TM1(Wilhoit wilhoit, double t1, double t2)

cpdef double Wilhoit_Dintegral_T1(Wilhoit wilhoit, double t1, double t2)

cpdef double Wilhoit_Dintegral_T2(Wilhoit wilhoit, double t1, double t2)

cpdef double Wilhoit_Dintegral_T3(Wilhoit wilhoit, double t1, double t2)

cpdef double Wilhoit_Dintegral_T4(Wilhoit wilhoit, double t1, double t2)

cpdef double Wilhoit_Dintegral2_T0(Wilhoit wilhoit, double t1, double t2)

cpdef double Wilhoit_Dintegral2_TM1(Wilhoit wilhoit, double t1, double t2)

cpdef double NASA_integral2_T0(NASA polynomial, double T)

cpdef double NASA_integral2_TM1(NASA polynomial, double T)

cpdef double NASA_Dintegral2_T0(NASA polynomial, double Ta, double Tb)

cpdef double NASA_Dintegral2_TM1(NASA polynomial, double Ta, double Tb)

