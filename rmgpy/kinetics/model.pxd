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

from rmgpy.quantity cimport ScalarQuantity, ArrayQuantity
from rmgpy.kinetics.uncertainties cimport RateUncertainty
from rmgpy.data.solvation import SoluteData
################################################################################

cpdef str get_rate_coefficient_units_from_reaction_order(n_gas, n_surf=?)

cpdef int get_reaction_order_from_rate_coefficient_units(kunits) except -1

################################################################################

cdef class KineticsModel:
    
    cdef public ScalarQuantity _Tmin, _Tmax
    cdef public ScalarQuantity _Pmin, _Pmax
    cdef public RateUncertainty uncertainty
    cdef public object solute

    cdef public str comment
    
    cpdef bint is_pressure_dependent(self) except -2
    
    cpdef bint is_temperature_valid(self, double T) except -2

    cpdef double get_rate_coefficient(self, double T, double P=?) except -1
    
    cpdef to_html(self)

    cpdef bint is_similar_to(self, KineticsModel other_kinetics) except -2

    cpdef bint is_identical_to(self, KineticsModel other_kinetics) except -2
    
    cpdef double discrepancy(self, KineticsModel other_kinetics) except -2
    

cdef class PDepKineticsModel(KineticsModel):
    
    cdef public dict efficiencies
    cdef public KineticsModel highPlimit
    
    cpdef bint is_pressure_dependent(self) except -2
    
    cpdef bint is_pressure_valid(self, double P) except -2

    cpdef double get_effective_pressure(self, double P, list species, np.ndarray fractions) except -1
    
    cpdef np.ndarray get_effective_collider_efficiencies(self, list species)

    cpdef double get_rate_coefficient(self, double T, double P=?) except -1

    cpdef to_html(self)

    cpdef bint is_similar_to(self, KineticsModel other_kinetics) except -2

    cpdef bint is_identical_to(self, KineticsModel other_kinetics) except -2

################################################################################

cdef class TunnelingModel:

    cdef public ScalarQuantity _frequency

    cpdef double calculate_tunneling_factor(self, double T) except -100000000

    cpdef np.ndarray calculate_tunneling_function(self, np.ndarray Elist)
