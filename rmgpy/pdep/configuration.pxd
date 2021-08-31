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

################################################################################

cdef class Configuration(object):

    cdef public list species
    cdef public np.ndarray e_list
    cdef public np.ndarray dens_states
    cdef public np.ndarray sum_states
    cdef public bint active_j_rotor
    cdef public bint active_k_rotor
    cdef public float energy_correction

    cpdef cleanup(self)

    cpdef bint is_unimolecular(self) except -2
    
    cpdef bint is_bimolecular(self) except -2

    cpdef bint is_termolecular(self) except -2

    cpdef bint is_transition_state(self) except -2
    
    cpdef bint has_statmech(self) except -2
    
    cpdef bint has_thermo(self) except -2
    
    cpdef double get_heat_capacity(self, double T) except -100000000

    cpdef double get_enthalpy(self, double T) except 100000000

    cpdef double get_entropy(self, double T) except -100000000

    cpdef double get_free_energy(self, double T) except 100000000
    
    cpdef double calculate_collision_frequency(self, double T, double P, dict bath_gas) except -1
        
    cpdef np.ndarray generate_collision_matrix(self, double T, np.ndarray dens_states,
                                             np.ndarray e_list, np.ndarray j_list=?)
    
    cpdef calculate_density_of_states(self, np.ndarray e_list, bint active_j_rotor=?, bint active_k_rotor=?, bint rmgmode=?)
