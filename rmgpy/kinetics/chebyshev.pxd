###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2018 Prof. William H. Green (whgreen@mit.edu),           #
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

cimport numpy

from rmgpy.kinetics.model cimport KineticsModel, PDepKineticsModel
from rmgpy.quantity cimport ScalarQuantity, ArrayQuantity

################################################################################

cdef class Chebyshev(PDepKineticsModel):
    
    cdef public ArrayQuantity _coeffs
    cdef public int degreeT
    cdef public int degreeP
    cdef public str kunits
    
    cdef double chebyshev(self, int n, double x)
    
    cdef double getReducedTemperature(self, double T) except -1000
    
    cdef double getReducedPressure(self, double P) except -1000
    
    cpdef double getRateCoefficient(self, double T, double P=?) except -1

    cpdef fitToData(self, numpy.ndarray Tlist, numpy.ndarray Plist, numpy.ndarray K, str kunits,
        int degreeT, int degreeP, double Tmin, double Tmax, double Pmin, double Pmax)

    cpdef bint isIdenticalTo(self, KineticsModel otherKinetics) except -2
    
    cpdef changeRate(self, double factor)
