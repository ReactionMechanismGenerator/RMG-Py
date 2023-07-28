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

cimport rmgpy.constants as constants
from rmgpy.quantity cimport ScalarQuantity, ArrayQuantity
from rmgpy.thermo.model cimport HeatCapacityModel
from rmgpy.statmech.conformer cimport Conformer
from rmgpy.kinetics.model cimport TunnelingModel
from rmgpy.molecule.molecule cimport Atom, Bond, Molecule
from rmgpy.molecule.graph cimport Graph

################################################################################

cdef class Species:
    
    cdef public int index
    cdef public str label
    cdef public object thermo
    cdef public Conformer conformer
    cdef public object transport_data
    cdef public list molecule
    cdef public ScalarQuantity _molecular_weight
    cdef public bint reactive
    cdef public object energy_transfer_model
    cdef public dict props
    cdef public str aug_inchi
    cdef public float symmetry_number
    cdef public bint is_solvent
    cdef public int creation_iteration
    cdef public bint explicitly_allowed
    cdef public object liquid_volumetric_mass_transfer_coefficient_data
    cdef public object henry_law_constant_data
    cdef str _fingerprint
    cdef str _inchi
    cdef str _smiles

    cpdef generate_resonance_structures(self, bint keep_isomorphic=?, bint filter_structures=?, bint save_order=?)

    cpdef get_net_charge(self)

    cpdef bint is_isomorphic(self, other, bint generate_initial_map=?, bint save_order=?, bint strict=?) except -2

    cpdef bint is_identical(self, other, bint strict=?) except -2

    cpdef bint is_structure_in_list(self, list species_list) except -2
    
    cpdef from_adjacency_list(self, adjlist, bint raise_atomtype_exception=?, bint raise_charge_exception=?)

    cpdef from_smiles(self, smiles)
    
    cpdef to_adjacency_list(self)
    
    cpdef bint contains_surface_site(self) except -2

    cpdef bint is_surface_site(self) except -2

    cpdef bint is_electron(self) except -2

    cpdef bint is_proton(self) except -2

    cpdef bint has_statmech(self) except -2

    cpdef bint has_thermo(self) except -2

    cpdef double get_partition_function(self, double T) except -1

    cpdef double get_heat_capacity(self, double T) except -100000000

    cpdef double get_enthalpy(self, double T) except 100000000

    cpdef double get_entropy(self, double T) except -100000000

    cpdef double get_free_energy(self, double T) except 100000000

    cpdef np.ndarray get_sum_of_states(self, np.ndarray e_list)

    cpdef np.ndarray get_density_of_states(self, np.ndarray e_list)

    cpdef double calculate_cp0(self) except -1

    cpdef double calculate_cpinf(self) except -1

    cpdef bint has_reactive_molecule(self) except -1

    cpdef Species copy(self, bint deep=?)

    cpdef set_structure(self, str structure)

    cpdef object get_henry_law_constant_data(self, list Ts=?)

    cpdef object get_liquid_volumetric_mass_transfer_coefficient_data(self, list Ts=?)
    
################################################################################

cdef class TransitionState:
    
    cdef public str label
    cdef public Conformer conformer
    cdef public ScalarQuantity _frequency
    cdef public int degeneracy
    cdef public TunnelingModel tunneling

    cpdef double get_partition_function(self, double T) except -1

    cpdef double get_heat_capacity(self, double T) except -100000000

    cpdef double get_enthalpy(self, double T) except 100000000

    cpdef double get_entropy(self, double T) except -100000000

    cpdef double get_free_energy(self, double T) except 100000000

    cpdef np.ndarray get_sum_of_states(self, np.ndarray e_list)

    cpdef np.ndarray get_density_of_states(self, np.ndarray e_list)
    
    cpdef double calculate_tunneling_factor(self, double T) except -1
    
    cpdef np.ndarray calculate_tunneling_function(self, np.ndarray e_list)
