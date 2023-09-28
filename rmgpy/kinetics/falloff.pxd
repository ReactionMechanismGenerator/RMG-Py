###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2023 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

from rmgpy.kinetics.model cimport KineticsModel, PDepKineticsModel
from rmgpy.kinetics.arrhenius cimport Arrhenius
from rmgpy.quantity cimport ScalarQuantity, ArrayQuantity

################################################################################

cdef class ThirdBody(PDepKineticsModel):
    
    cdef public Arrhenius arrheniusLow
    
    cpdef double get_rate_coefficient(self, double T, double P=?) except -1

    cpdef bint is_identical_to(self, KineticsModel other_kinetics) except -2
    
    cpdef change_rate(self, double factor)

################################################################################

cdef class Lindemann(PDepKineticsModel):
    
    cdef public Arrhenius arrheniusHigh
    cdef public Arrhenius arrheniusLow
    
    cpdef double get_rate_coefficient(self, double T, double P=?) except -1

    cpdef bint is_identical_to(self, KineticsModel other_kinetics) except -2
    
    cpdef change_rate(self, double factor)

################################################################################

cdef class Troe(PDepKineticsModel):
    
    cdef public Arrhenius arrheniusHigh
    cdef public Arrhenius arrheniusLow
    cdef public double alpha
    cdef public ScalarQuantity _T1, _T2, _T3
    
    cpdef double get_rate_coefficient(self, double T, double P=?) except -1

    cpdef bint is_identical_to(self, KineticsModel other_kinetics) except -2
    
    cpdef change_rate(self, double factor)
