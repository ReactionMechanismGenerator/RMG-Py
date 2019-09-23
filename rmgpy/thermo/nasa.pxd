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

from rmgpy.thermo.model cimport HeatCapacityModel
from rmgpy.thermo.thermodata cimport ThermoData
from rmgpy.thermo.wilhoit cimport Wilhoit

################################################################################

cdef class NASAPolynomial(HeatCapacityModel):

    cdef public double cm2, cm1, c0, c1, c2, c3, c4, c5, c6
    
    cpdef double get_heat_capacity(self, double T) except -1000000000

    cpdef double get_enthalpy(self, double T) except 1000000000

    cpdef double get_entropy(self, double T) except -1000000000

    cpdef double get_free_energy(self, double T) except 1000000000    
    
    cpdef change_base_enthalpy(self, double deltaH)

    cpdef change_base_entropy(self, double deltaS)

    cdef double integral2_T0(self, double T)
    
    cdef double integral2_TM1(self, double T)
    
cdef class NASA(HeatCapacityModel):

    cdef public NASAPolynomial poly1, poly2, poly3
    
    cpdef NASAPolynomial select_polynomial(self, double T)

    cpdef dict as_dict(self)

    cpdef double get_heat_capacity(self, double T) except -1000000000

    cpdef double get_enthalpy(self, double T) except 1000000000

    cpdef double get_entropy(self, double T) except -1000000000

    cpdef double get_free_energy(self, double T) except 1000000000

    cpdef ThermoData to_thermo_data(self)

    cpdef Wilhoit to_wilhoit(self)
    
    cpdef NASA change_base_enthalpy(self, double deltaH)

    cpdef NASA change_base_entropy(self, double deltaS)
