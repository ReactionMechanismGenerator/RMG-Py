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

cimport rmgpy.constants as constants
from rmgpy.species cimport Species, TransitionState
from rmgpy.molecule.molecule cimport Atom, Molecule
from rmgpy.molecule.graph cimport Vertex, Graph
from rmgpy.molecule.element cimport Element
from rmgpy.kinetics.model cimport KineticsModel
from rmgpy.kinetics.arrhenius cimport Arrhenius
from rmgpy.kinetics.surface cimport SurfaceArrhenius, StickingCoefficient

cimport numpy as np

################################################################################

cdef class Reaction:
    cdef public int index
    cdef public str label
    cdef public list reactants
    cdef public list products
    cdef public Species specific_collider
    cdef public bint reversible
    cdef public TransitionState transition_state
    cdef public KineticsModel kinetics
    cdef public Arrhenius network_kinetics
    cdef public SurfaceArrhenius
    cdef public bint duplicate
    cdef public float _degeneracy
    cdef public list pairs
    cdef public bint allow_pdep_route
    cdef public bint elementary_high_p
    cdef public str comment
    cdef public dict k_effective_cache
    cdef public bint is_forward
    cdef public bint allow_max_rate_violation
    cdef public object rank

    cpdef bint is_isomerization(self)

    cpdef bint is_dissociation(self)

    cpdef bint is_association(self)

    cpdef bint is_unimolecular(self)

    cpdef bint is_surface_reaction(self)

    cpdef bint has_template(self, list reactants, list products)

    cpdef bint matches_species(self, list reactants, list products=?)

    cpdef bint is_isomorphic(self, Reaction other, bint either_direction=?, bint check_identical=?,
                             bint check_only_label=?, bint check_template_rxn_products=?, bint generate_initial_map=?,
                             bint strict=?, bint save_order=?) except -2

    cpdef double get_enthalpy_of_reaction(self, double T)

    cpdef double get_entropy_of_reaction(self, double T)

    cpdef double get_free_energy_of_reaction(self, double T)

    cpdef double get_equilibrium_constant(self, double T, str type=?, double surface_site_density=?)

    cpdef np.ndarray get_enthalpies_of_reaction(self, np.ndarray Tlist)

    cpdef np.ndarray get_entropies_of_reaction(self, np.ndarray Tlist)

    cpdef np.ndarray get_free_energies_of_reaction(self, np.ndarray Tlist)

    cpdef np.ndarray get_equilibrium_constants(self, np.ndarray Tlist, str type=?)

    cpdef int get_stoichiometric_coefficient(self, Species spec)

    cpdef double get_rate_coefficient(self, double T, double P=?, double surface_site_density=?)

    cpdef double get_surface_rate_coefficient(self, double T, double surface_site_density) except -2

    cpdef fix_barrier_height(self, bint force_positive=?)

    cpdef reverse_arrhenius_rate(self, Arrhenius k_forward, str reverse_units, Tmin=?, Tmax=?)

    cpdef reverse_surface_arrhenius_rate(self, SurfaceArrhenius k_forward, str reverse_units, Tmin=?, Tmax=?)

    cpdef reverse_sticking_coeff_rate(self, StickingCoefficient k_forward, str reverse_units, double surface_site_density, Tmin=?, Tmax=?)

    cpdef generate_reverse_rate_coefficient(self, bint network_kinetics=?, Tmin=?, Tmax=?, double surface_site_density=?)

    cpdef np.ndarray calculate_tst_rate_coefficients(self, np.ndarray Tlist)

    cpdef double calculate_tst_rate_coefficient(self, double T) except -2

    cpdef bint can_tst(self) except -2

    cpdef calculate_microcanonical_rate_coefficient(self, np.ndarray e_list, np.ndarray j_list,
                                                    np.ndarray reac_dens_states, np.ndarray prod_dens_states=?,
                                                    double T=?)

    cpdef bint is_balanced(self)

    cpdef generate_pairs(self)

    cpdef copy(self)

    cpdef ensure_species(self, bint reactant_resonance=?, bint product_resonance=?, bint save_order=?)

    cpdef list check_collision_limit_violation(self, float t_min, float t_max, float p_min, float p_max)

    cpdef calculate_coll_limit(self, float temp, bint reverse=?)

    cpdef get_reduced_mass(self, bint reverse=?)

    cpdef get_mean_sigma_and_epsilon(self, bint reverse=?)

cpdef bint same_species_lists(list list1, list list2, bint check_identical=?, bint only_check_label=?,
                              bint generate_initial_map=?, bint strict=?, bint save_order=?) except -2
