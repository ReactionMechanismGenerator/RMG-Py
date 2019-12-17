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

from rmgpy.quantity cimport ScalarQuantity, ArrayQuantity
from rmgpy.statmech.mode cimport Mode

################################################################################

cdef class Torsion(Mode):
    cpdef make_object(self, dict data, dict class_dict)

################################################################################

cdef class HinderedRotor(Torsion):

    cdef public ScalarQuantity _inertia
    cdef public ArrayQuantity _fourier
    cdef public ScalarQuantity _barrier
    cdef public int symmetry
    cdef public bint semiclassical
    cdef public double frequency
    cdef public np.ndarray energies
    
    cpdef double get_frequency(self) except -1

    cpdef double get_level_energy(self, int J) except -1
    
    cpdef int get_level_degeneracy(self, int J) except -1
    
    cpdef np.ndarray solve_schrodinger_equation(self, int n_basis=?)
    
    cpdef np.ndarray get_hamiltonian(self, int n_basis)
    
    cdef double get_rotational_constant_energy(self)
    
    cpdef double get_potential(self, double phi) except -100000000

    cpdef double get_partition_function(self, double T) except -1
        
    cpdef double get_heat_capacity(self, double T) except -100000000

    cpdef double get_enthalpy(self, double T) except 100000000

    cpdef double get_entropy(self, double T) except -100000000

    cpdef np.ndarray get_sum_of_states(self, np.ndarray e_list, np.ndarray sum_states_0=?)
    
    cpdef np.ndarray get_density_of_states(self, np.ndarray e_list, np.ndarray dens_states_0=?)

    cpdef fit_fourier_potential_to_data(self, np.ndarray angle, np.ndarray V)
    
    cpdef fit_cosine_potential_to_data(self, np.ndarray angle, np.ndarray V)

cdef class FreeRotor(Torsion):
    
    cdef public ScalarQuantity _inertia
    cdef public int symmetry
    
    cdef double get_rotational_constant_energy(self)
    
    cpdef double get_partition_function(self, double T) except -1
        
    cpdef double get_heat_capacity(self, double T) except -100000000

    cpdef double get_enthalpy(self, double T) except 100000000

    cpdef double get_entropy(self, double T) except -100000000