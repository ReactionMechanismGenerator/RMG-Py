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
from rmgpy.thermo.wilhoit cimport Wilhoit
from rmgpy.thermo.nasa cimport NASA
from rmgpy.quantity cimport ScalarQuantity, ArrayQuantity

################################################################################

cdef class ThermoData(HeatCapacityModel):
    
    cdef public ScalarQuantity _H298, _S298
    cdef public ArrayQuantity _Tdata, _Cpdata
    cdef public ScalarQuantity _Cp0, _CpInf
    
    cpdef double getHeatCapacity(self, double T) except -1000000000

    cpdef double getEnthalpy(self, double T) except 1000000000

    cpdef double getEntropy(self, double T) except -1000000000

    cpdef double getFreeEnergy(self, double T) except 1000000000

    cpdef Wilhoit toWilhoit(self)

    cpdef NASA toNASA(self, double Tmin, double Tmax, double Tint, bint fixedTint=?, bint weighting=?, int continuity=?)
