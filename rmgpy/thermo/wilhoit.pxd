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

cimport numpy as np

from rmgpy.thermo.model cimport HeatCapacityModel
from rmgpy.thermo.thermodata cimport ThermoData
from rmgpy.thermo.nasa cimport NASA, NASAPolynomial
from rmgpy.quantity cimport ScalarQuantity, ArrayQuantity

################################################################################

cdef class Wilhoit(HeatCapacityModel):
    
    cdef public ScalarQuantity _B, _H0, _S0
    cdef public double a0, a1, a2, a3
    cdef public dict _thermo_coverage_dependence

    cpdef dict as_dict(self)

    cpdef make_object(self, dict data, dict class_dict)
    
    cpdef double get_heat_capacity(self, double T) except -1000000000

    cpdef double get_enthalpy(self, double T) except 1000000000

    cpdef double get_entropy(self, double T) except -1000000000

    cpdef double get_free_energy(self, double T) except 1000000000
    
    cpdef Wilhoit copy(self)
    
    cdef double integral_T0(self, double T)
    
    cdef double integral_TM1(self, double T)
    
    cdef double integral_T1(self, double T)
    
    cdef double integral_T2(self, double T)
    
    cdef double integral_T3(self, double T)
    
    cdef double integral_T4(self, double T)
    
    cdef double integral2_T0(self, double T)
    
    cdef double integral2_TM1(self, double T)
    
    cpdef ThermoData to_thermo_data(self)

    cpdef NASA to_nasa(self, double Tmin, double Tmax, double Tint, bint fixedTint=?, bint weighting=?, int continuity=?)

################################################################################

cpdef wilhoit_to_nasa(Wilhoit wilhoit, double Tmin, double Tmax, double Tint, bint weighting, int cont_cons)

cpdef wilhoit_to_nasa_t_int_opt(Wilhoit wilhoit, double Tmin, double Tmax, bint weighting, int cont_cons)

cpdef double wilhoit_to_nasa_t_int_opt_obj_fun(double Tint, Wilhoit wilhoit, double Tmin, double Tmax, bint weighting, int cont_cons)

cpdef double wilhoit_to_nasa_t_int_opt_obj_fun_nw(double Tint, Wilhoit wilhoit, double Tmin, double Tmax, int cont_cons)

cpdef double wilhoit_to_nasa_t_int_opt_obj_fun_w(double Tint, Wilhoit wilhoit, double Tmin, double Tmax, int cont_cons)
