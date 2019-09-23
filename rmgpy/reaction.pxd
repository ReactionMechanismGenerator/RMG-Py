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

cimport rmgpy.constants as constants
from rmgpy.species cimport Species, TransitionState
from rmgpy.molecule.molecule cimport Atom, Molecule
from rmgpy.molecule.element cimport Element
from rmgpy.kinetics.model cimport KineticsModel
from rmgpy.kinetics.arrhenius cimport Arrhenius

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
    
    cpdef bint isIsomerization(self)

    cpdef bint isDissociation(self)

    cpdef bint isAssociation(self)
    
    cpdef bint isUnimolecular(self)

    cpdef bint isSurfaceReaction(self)

    cpdef bint hasTemplate(self, list reactants, list products)
    
    cpdef bint matchesSpecies(self, list reactants, list products=?)

    cpdef bint is_isomorphic(self, Reaction other, bint either_direction=?, bint check_identical=?, bint check_only_label=?,
                            bint check_template_rxn_products=?, bint generate_initial_map=?, bint strict=?) except -2

    cpdef double getEnthalpyOfReaction(self, double T)

    cpdef double getEntropyOfReaction(self, double T)

    cpdef double getFreeEnergyOfReaction(self, double T)

    cpdef double getEquilibriumConstant(self, double T, str type=?)

    cpdef np.ndarray getEnthalpiesOfReaction(self, np.ndarray Tlist)

    cpdef np.ndarray getEntropiesOfReaction(self, np.ndarray Tlist)

    cpdef np.ndarray getFreeEnergiesOfReaction(self, np.ndarray Tlist)

    cpdef np.ndarray getEquilibriumConstants(self, np.ndarray Tlist, str type=?)

    cpdef int getStoichiometricCoefficient(self, Species spec)

    cpdef double get_rate_coefficient(self, double T, double P=?)

    cpdef double getSurfaceRateCoefficient(self, double T, double surfaceSiteDensity) except -2

    cpdef fixBarrierHeight(self, bint forcePositive=?)

    cpdef reverseThisArrheniusRate(self, Arrhenius kForward, str reverseUnits, Tmin=?, Tmax=?)

    cpdef generateReverseRateCoefficient(self, bint network_kinetics=?, Tmin=?, Tmax=?)

    cpdef np.ndarray calculateTSTRateCoefficients(self, np.ndarray Tlist)

    cpdef double calculateTSTRateCoefficient(self, double T) except -2

    cpdef bint canTST(self) except -2

    cpdef calculateMicrocanonicalRateCoefficient(self, np.ndarray Elist, np.ndarray Jlist, np.ndarray reacDensStates, np.ndarray prodDensStates=?, double T=?)

    cpdef bint isBalanced(self)
    
    cpdef generatePairs(self)
    
    cpdef copy(self)

    cpdef ensure_species(self, bint reactant_resonance=?, bint product_resonance=?)

    cpdef list check_collision_limit_violation(self, float t_min, float t_max, float p_min, float p_max)

    cpdef calculate_coll_limit(self, float temp, bint reverse=?)

    cpdef get_reduced_mass(self, bint reverse=?)

    cpdef get_mean_sigma_and_epsilon(self, bint reverse=?)

cpdef bint same_species_lists(list list1, list list2, bint check_identical=?, bint only_check_label=?, bint generate_initial_map=?, bint strict=?) except -2
