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

from rmgpy.kinetics.model cimport KineticsModel, PDepKineticsModel
from rmgpy.quantity cimport ScalarQuantity, ArrayQuantity

################################################################################

cdef class Arrhenius(KineticsModel):
    
    cdef public ScalarQuantity _A
    cdef public ScalarQuantity _n
    cdef public ScalarQuantity _Ea
    cdef public ScalarQuantity _T0
    
    cpdef double getRateCoefficient(self, double T, double P=?) except -1

    cpdef changeT0(self, double T0)

    cpdef fitToData(self, numpy.ndarray Tlist, numpy.ndarray klist, str kunits, double T0=?, numpy.ndarray weights=?, bint threeParams=?)

    cpdef bint isIdenticalTo(self, KineticsModel otherKinetics) except -2

################################################################################

cdef class ArrheniusEP(KineticsModel):
    
    cdef public ScalarQuantity _A
    cdef public ScalarQuantity _n
    cdef public ScalarQuantity _alpha
    cdef public ScalarQuantity _E0
    
    cpdef double getRateCoefficient(self, double T, double dHrxn=?) except -1

    cpdef double getActivationEnergy(self, double dHrxn) except -1
    
    cpdef Arrhenius toArrhenius(self, double dHrxn)

    cpdef bint isIdenticalTo(self, KineticsModel otherKinetics) except -2
    
################################################################################

cdef class PDepArrhenius(PDepKineticsModel):
    
    cdef public ArrayQuantity _pressures
    cdef public list arrhenius
    
    cdef getAdjacentExpressions(self, double P)
    
    cpdef double getRateCoefficient(self, double T, double P=?) except -1
    
    cpdef fitToData(self, numpy.ndarray Tlist, numpy.ndarray Plist, numpy.ndarray K, str kunits, double T0=?)

    cpdef bint isIdenticalTo(self, KineticsModel otherKinetics) except -2

################################################################################

cdef class MultiArrhenius(KineticsModel):
    
    cdef public list arrhenius
    
    cpdef double getRateCoefficient(self, double T, double P=?) except -1

    cpdef bint isIdenticalTo(self, KineticsModel otherKinetics) except -2
    
    cpdef Arrhenius toArrhenius(self, double Tmin=?, double Tmax=?)

################################################################################

cdef class MultiPDepArrhenius(PDepKineticsModel):
    
    cdef public list arrhenius
    
    cpdef double getRateCoefficient(self, double T, double P=?) except -1

    cpdef bint isIdenticalTo(self, KineticsModel otherKinetics) except -2
    