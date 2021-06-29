###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2021 Prof. William H. Green (whgreen@mit.edu),           #
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

from rmgpy.kinetics.model cimport KineticsModel, PDepKineticsModel
from rmgpy.quantity cimport ScalarQuantity, ArrayQuantity

################################################################################

cdef class Arrhenius(KineticsModel):
    
    cdef public ScalarQuantity _A
    cdef public ScalarQuantity _n
    cdef public ScalarQuantity _Ea
    cdef public ScalarQuantity _T0
    
    cpdef double get_rate_coefficient(self, double T, double P=?) except -1

    cpdef change_t0(self, double T0)

    cpdef fit_to_data(self, np.ndarray Tlist, np.ndarray klist, str kunits, double T0=?, np.ndarray weights=?, bint three_params=?)

    cpdef bint is_identical_to(self, KineticsModel other_kinetics) except -2
    
    cpdef change_rate(self, double factor)

    cpdef ArrheniusEP to_arrhenius_ep(self, double alpha=?, double dHrxn=?)

################################################################################

cdef class ArrheniusEP(KineticsModel):
    
    cdef public ScalarQuantity _A
    cdef public ScalarQuantity _n
    cdef public ScalarQuantity _alpha
    cdef public ScalarQuantity _E0
    
    cpdef double get_rate_coefficient(self, double T, double dHrxn=?) except -1

    cpdef double get_activation_energy(self, double dHrxn) except -1
    
    cpdef Arrhenius to_arrhenius(self, double dHrxn)

    cpdef bint is_identical_to(self, KineticsModel other_kinetics) except -2
    
    cpdef change_rate(self, double factor)
################################################################################

cdef class ArrheniusBM(KineticsModel):
    
    cdef public ScalarQuantity _A
    cdef public ScalarQuantity _n
    cdef public ScalarQuantity _w0
    cdef public ScalarQuantity _E0
    
    cpdef double get_rate_coefficient(self, double T, double dHrxn=?) except -1

    cpdef double get_activation_energy(self, double dHrxn) except -1
    
    cpdef Arrhenius to_arrhenius(self, double dHrxn)

    cpdef bint is_identical_to(self, KineticsModel other_kinetics) except -2
    
    cpdef change_rate(self, double factor)
################################################################################

cdef class PDepArrhenius(PDepKineticsModel):
    
    cdef public ArrayQuantity _pressures
    cdef public list arrhenius
    
    cdef get_adjacent_expressions(self, double P)
    
    cpdef double get_rate_coefficient(self, double T, double P=?) except -1
    
    cpdef fit_to_data(self, np.ndarray Tlist, np.ndarray Plist, np.ndarray K, str kunits, double T0=?)

    cpdef bint is_identical_to(self, KineticsModel other_kinetics) except -2
    
    cpdef change_rate(self, double factor)

################################################################################

cdef class MultiArrhenius(KineticsModel):
    
    cdef public list arrhenius
    
    cpdef double get_rate_coefficient(self, double T, double P=?) except -1

    cpdef bint is_identical_to(self, KineticsModel other_kinetics) except -2
    
    cpdef Arrhenius to_arrhenius(self, double Tmin=?, double Tmax=?)
    
    cpdef change_rate(self, double factor)

################################################################################

cdef class MultiPDepArrhenius(PDepKineticsModel):
    
    cdef public list arrhenius
    
    cpdef double get_rate_coefficient(self, double T, double P=?) except -1

    cpdef bint is_identical_to(self, KineticsModel other_kinetics) except -2
    
    cpdef change_rate(self, double factor)
