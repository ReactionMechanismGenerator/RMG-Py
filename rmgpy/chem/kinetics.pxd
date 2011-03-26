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

from constants cimport Quantity

cdef extern from "math.h":
    cdef double acos(double x)
    cdef double cos(double x)
    cdef double exp(double x)
    cdef double log(double x)
    cdef double log10(double x)
    cdef double pow(double base, double exponent)

################################################################################

cdef class KineticsModel:
    
    cdef public Quantity Tmin, Tmax, Pmin, Pmax
    cdef public str comment
    
    cpdef bint isTemperatureValid(self, double T) except -2
    
    cpdef bint isPressureValid(self, double P) except -2

    cpdef bint isPressureDependent(self)
    
    cpdef numpy.ndarray getRateCoefficients(self, numpy.ndarray Tlist)

################################################################################

cdef class Arrhenius(KineticsModel):
    
    cdef public Quantity A, T0, Ea, n
    
    cpdef bint isPressureDependent(self)

    cpdef double getRateCoefficient(self, double T, double P=?)

    cpdef changeT0(self, double T0)

    cpdef fitToData(self, Quantity Tlist, Quantity klist, double T0=?)

################################################################################

cdef class ArrheniusEP(KineticsModel):
    
    cdef public Quantity A, E0, n, alpha

    cpdef bint isPressureDependent(self)

    cpdef double getActivationEnergy(self, double dHrxn)
    
    cpdef double getRateCoefficient(self, double T, double dHrxn)

################################################################################

cdef class MultiArrhenius(KineticsModel):

    cdef public list arrheniusList

    cpdef bint isPressureDependent(self)

    cpdef double getRateCoefficient(self, double T, double P=?)

################################################################################

cdef class PDepArrhenius(KineticsModel):
    
    cdef public Quantity pressures
    cdef public list arrhenius
    
    cpdef bint isPressureDependent(self)

    cpdef tuple __getAdjacentExpressions(self, double P)
    
    cpdef double getRateCoefficient(self, double T, double P)

    cpdef fitToData(self, numpy.ndarray Tlist, numpy.ndarray Plist, numpy.ndarray K, double T0=?)

################################################################################

cdef class Chebyshev(KineticsModel):
    
    cdef public object coeffs
    cdef public int degreeT
    cdef public int degreeP
    
    cpdef bint isPressureDependent(self)

    cpdef double __chebyshev(self, double n, double x)
    
    cpdef double __getReducedTemperature(self, double T)
    
    cpdef double __getReducedPressure(self, double P)
    
    cpdef double getRateCoefficient(self, double T, double P)

    cpdef fitToData(self, numpy.ndarray Tlist, numpy.ndarray Plist, numpy.ndarray K,
        int degreeT, int degreeP, double Tmin, double Tmax, double Pmin, double Pmax)

################################################################################

cdef class ThirdBody(KineticsModel):

    cdef public Arrhenius arrheniusHigh
    cdef public dict efficiencies
    
    cpdef bint isPressureDependent(self)

    cpdef getColliderEfficiency(self, collider)

    cpdef double getRateCoefficient(self, double T, double P, collider=?)

################################################################################

cdef class Lindemann(ThirdBody):

    cdef public Arrhenius arrheniusLow
    
    cpdef double getRateCoefficient(self, double T, double P, collider=?)

################################################################################

cdef class Troe(Lindemann):

    cdef public Quantity alpha, T1, T2, T3
    
    cpdef double getRateCoefficient(self, double T, double P, collider=?)
