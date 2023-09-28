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

from rmgpy.statmech.mode cimport Mode
from rmgpy.quantity cimport ScalarQuantity, ArrayQuantity

################################################################################

cdef class Rotation(Mode):

    cdef public int symmetry

    cpdef make_object(self, dict data, dict class_dict)

################################################################################

cdef class LinearRotor(Rotation):

    cdef public ScalarQuantity _inertia

    cpdef double get_level_energy(self, int J) except -1
    
    cpdef int get_level_degeneracy(self, int J) except -1
    
    cpdef double get_partition_function(self, double T) except -1
        
    cpdef double get_heat_capacity(self, double T) except -100000000

    cpdef double get_enthalpy(self, double T) except 100000000

    cpdef double get_entropy(self, double T) except -100000000

    cpdef np.ndarray get_sum_of_states(self, np.ndarray e_list, np.ndarray sum_states_0=?)
    
    cpdef np.ndarray get_density_of_states(self, np.ndarray e_list, np.ndarray dens_states_0=?)

################################################################################

cdef class NonlinearRotor(Rotation):

    cdef public ArrayQuantity _inertia

    cdef np.ndarray get_rotational_constant_energy(self)

    cpdef double get_partition_function(self, double T) except -1

    cpdef double get_heat_capacity(self, double T) except -100000000

    cpdef double get_enthalpy(self, double T) except 100000000

    cpdef double get_entropy(self, double T) except -100000000

    cpdef np.ndarray get_sum_of_states(self, np.ndarray e_list, np.ndarray sum_states_0=?)

    cpdef np.ndarray get_density_of_states(self, np.ndarray e_list, np.ndarray dens_states_0=?)

################################################################################

cdef class KRotor(Rotation):

    cdef public ScalarQuantity _inertia

    cpdef double get_level_energy(self, int J) except -1

    cpdef int get_level_degeneracy(self, int J) except -1

    cpdef double get_partition_function(self, double T) except -1

    cpdef double get_heat_capacity(self, double T) except -100000000

    cpdef double get_enthalpy(self, double T) except 100000000

    cpdef double get_entropy(self, double T) except -100000000

    cpdef np.ndarray get_sum_of_states(self, np.ndarray e_list, np.ndarray sum_states_0=?)

    cpdef np.ndarray get_density_of_states(self, np.ndarray e_list, np.ndarray dens_states_0=?)

################################################################################

cdef class SphericalTopRotor(Rotation):

    cdef public ScalarQuantity _inertia

    cpdef double get_level_energy(self, int J) except -1
    
    cpdef int get_level_degeneracy(self, int J) except -1
    
    cpdef double get_partition_function(self, double T) except -1
    
    cpdef double get_heat_capacity(self, double T) except -100000000

    cpdef double get_enthalpy(self, double T) except 100000000

    cpdef double get_entropy(self, double T) except -100000000

    cpdef np.ndarray get_sum_of_states(self, np.ndarray e_list, np.ndarray sum_states_0=?)
    
    cpdef np.ndarray get_density_of_states(self, np.ndarray e_list, np.ndarray dens_states_0=?)
