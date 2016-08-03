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

from rmgpy.quantity cimport ScalarQuantity, ArrayQuantity

################################################################################

cpdef str getRateCoefficientUnitsFromReactionOrder(order)

cpdef int getReactionOrderFromRateCoefficientUnits(kunits) except -1

################################################################################

cdef class KineticsModel:
    
    cdef public ScalarQuantity _Tmin, _Tmax
    cdef public ScalarQuantity _Pmin, _Pmax
    cdef public str comment
    
    cpdef bint isPressureDependent(self) except -2
    
    cpdef bint isTemperatureValid(self, double T) except -2

    cpdef double getRateCoefficient(self, double T, double P=?) except -1
    
    cpdef toHTML(self)

    cpdef bint isSimilarTo(self, KineticsModel otherKinetics) except -2

    cpdef bint isIdenticalTo(self, KineticsModel otherKinetics) except -2
    
    cpdef double discrepancy(self, KineticsModel otherKinetics) except -2

cdef class PDepKineticsModel(KineticsModel):
    
    cdef public dict efficiencies
    cdef public KineticsModel highPlimit
    
    cpdef bint isPressureDependent(self) except -2
    
    cpdef bint isPressureValid(self, double P) except -2

    cpdef double getEffectivePressure(self, double P, list species, numpy.ndarray fractions) except -1
    
    cpdef numpy.ndarray getEffectiveColliderEfficiencies(self, list species)

    cpdef double getRateCoefficient(self, double T, double P=?) except -1

    cpdef toHTML(self)

    cpdef bint isSimilarTo(self, KineticsModel otherKinetics) except -2

    cpdef bint isIdenticalTo(self, KineticsModel otherKinetics) except -2

################################################################################

cdef class TunnelingModel:

    cdef public ScalarQuantity _frequency

    cpdef double calculateTunnelingFactor(self, double T) except -100000000

    cpdef numpy.ndarray calculateTunnelingFunction(self, numpy.ndarray Elist)
