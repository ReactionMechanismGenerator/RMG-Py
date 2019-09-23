###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2019 Prof. William H. Green (whgreen@mit.edu),           #
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

cimport numpy as np

from rmgpy.rmgobject cimport RMGObject


cpdef list NOT_IMPLEMENTED_UNITS

################################################################################

cdef class Units(RMGObject):

    cdef public str units

    cpdef double get_conversion_factor_to_si(self) except -1

    cpdef double get_conversion_factor_from_si(self) except -1

    cpdef double get_conversion_factor_from_si_to_cm_mol_s(self) except -1
    
################################################################################

cdef class ScalarQuantity(Units):

    cdef public double value_si
    cdef        str _uncertaintyType
    cdef public double uncertainty_si

    cpdef dict as_dict(self)

    cpdef str getUncertaintyType(self)
    cpdef     setUncertaintyType(self, str v)
    
    cpdef bint is_uncertainty_additive(self) except -2

    cpdef bint is_uncertainty_multiplicative(self) except -2
    
    cpdef ScalarQuantity copy(self)

################################################################################

cdef class ArrayQuantity(Units):

    cdef public np.ndarray value_si
    cdef public str _uncertaintyType
    cdef public np.ndarray uncertainty_si

    cpdef dict as_dict(self)

    cpdef str getUncertaintyType(self)
    cpdef     setUncertaintyType(self, str v)

    cpdef bint is_uncertainty_additive(self) except -2

    cpdef bint is_uncertainty_multiplicative(self) except -2
    
    cpdef ArrayQuantity copy(self)
