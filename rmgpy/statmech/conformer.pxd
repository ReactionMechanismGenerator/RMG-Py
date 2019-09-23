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
from rmgpy.quantity cimport ScalarQuantity, ArrayQuantity

################################################################################

cdef class Conformer(RMGObject):

    cdef public ScalarQuantity _E0
    cdef public list modes
    cdef public int spin_multiplicity
    cdef public int optical_isomers
    cdef public ArrayQuantity _number
    cdef public ArrayQuantity _mass
    cdef public ArrayQuantity _coordinates

    cpdef double get_partition_function(self, double T) except -1

    cpdef double get_heat_capacity(self, double T) except -100000000

    cpdef double get_enthalpy(self, double T) except 100000000

    cpdef double get_entropy(self, double T) except -100000000

    cpdef double get_free_energy(self, double T) except 100000000

    cpdef np.ndarray get_sum_of_states(self, np.ndarray e_list)

    cpdef np.ndarray get_density_of_states(self, np.ndarray e_list)

    cpdef double get_total_mass(self, atoms=?) except -1

    cpdef np.ndarray get_center_of_mass(self, atoms=?)

    cpdef np.ndarray get_moment_of_inertia_tensor(self)

    cpdef get_principal_moments_of_inertia(self)

    cpdef double get_internal_reduced_moment_of_inertia(self, pivots, top1, option=?) except -1

    cpdef get_symmetric_top_rotors(self)

    cpdef list get_active_modes(self, bint active_j_rotor=?, bint active_k_rotor=?)
    
    cpdef get_number_degrees_of_freedom(self)

################################################################################

cpdef double phi(double beta, int k, double E, logQ) except -10000000
