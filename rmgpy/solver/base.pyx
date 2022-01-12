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

"""
Contains the :class:`ReactionSystem` class, a base class for all RMG reaction
systems.
"""

cdef extern from "math.h":
    double sqrt(double)

import csv
import itertools
import logging

import cython
import numpy as np
cimport numpy as np
from cpython cimport bool

include "settings.pxi"
if DASPK == 1:
    from pydas.daspk cimport DASPK as DASx
    from pydas.daspk import DASPKError as DASxError
else:
    from pydas.dassl cimport DASSL as DASx
    from pydas.daspk import DASSLError as DASxError

import rmgpy.constants as constants
cimport rmgpy.constants as constants
from rmgpy.chemkin import get_species_identifier
from rmgpy.reaction import Reaction
from rmgpy.quantity import Quantity
from rmgpy.species import Species
from rmgpy.solver.termination import TerminationTime, TerminationConversion, TerminationRateRatio
################################################################################

cdef class ReactionSystem(DASx):
    """
    A base class for all RMG reaction systems.
    """

    def __init__(self, termination=None, sensitive_species=None, sensitivity_threshold=1e-3):
        DASx.__init__(self)

        # reactor state variables:
        self.t0 = 0.0
        self.y0 = None
        self.dydt0 = None

        #  variables that determine the dimensions of arrays and matrices:
        self.num_core_species = -1
        self.num_core_reactions = -1
        self.num_edge_species = -1
        self.num_edge_reactions = -1
        self.num_pdep_networks = -1

        """
        The number of differential equations that will 
        be solved simulatenously by the solver.

        """
        self.neq = -1

        # variables that store stoichiometry data
        """
        species_index is a dictionary with species as keys, and 
        values equal to the index integer that is generated
        when enumerating the core and edge species .
        """
        self.species_index = {}

        """
        reaction_index is a dictionary with reactions as keys, and 
        values equal to the index integer that is generated
        when enumerating the core and edge reactions .
        """
        self.reaction_index = {}

        """
        A matrix for the reactants and products.

        The matrix has dimensions n x 3, with :
        - n the sum of the number of core and edge reactions,
        - 3 the maximum number of molecules allowed in either the reactant or
            product side of a reaction.

        The matrix element (j,l), with
        - j the index of the reaction and 
        - l the index of the lth reactant/product in the reaction contains the index
            of the molecule in the species_index dictionary.
        """
        self.reactant_indices = None
        self.product_indices = None

        self.network_indices = None

        # matrices that cache kinetic and rate data
        self.kf = None  # forward rate coefficients
        self.kb = None  # reverse rate coefficients
        self.Keq = None  # equilibrium constants
        self.network_leak_coefficients = None
        self.jacobian_matrix = None

        self.core_species_concentrations = None

        # The reaction and species rates at the current time (in mol/m^3*s)
        self.core_species_rates = None
        self.core_reaction_rates = None
        self.core_species_production_rates = None
        self.core_species_consumption_rates = None
        self.edge_species_rates = None
        self.edge_reaction_rates = None

        self.network_leak_rates = None

        #surface indices
        self.surface_species_indices = None
        self.surface_reaction_indices = None
        self.valid_layering_indices = None

        # variables that cache maximum rate (ratio) data
        self.max_edge_species_rate_ratios = None
        self.max_network_leak_rate_ratios = None

        #for managing prunable edge species
        self.prunable_species = []
        self.prunable_networks = []
        self.prunable_species_indices = None
        self.prunable_network_indices = None

        # sensitivity variables
        self.sensmethod = 2  # sensmethod = 1 for staggered corrector sensitivities, 0 (simultaneous corrector), 2 (staggered direct)
        self.sensitivity_coefficients = None
        self.sensitive_species = sensitive_species
        self.sensitivity_threshold = sensitivity_threshold
        self.senpar = None

        # tolerance settings

        """
        Arrays containing the absolute and relative tolerances.
        The dimension of the arrays correspond to the number of 
        differential equations in the system.
        """
        self.atol_array = None
        self.rtol_array = None

        self.termination = termination or []

        # Flag to indicate whether or not reactions with 3 reactants are present 
        self.trimolecular = False

        # reaction filtration, unimolecular_threshold is a vector with length of number of core species
        # bimolecular_threshold is a square matrix with length of number of core species
        # trimolecular_threshold is a rank-3 tensor with length of number of core species
        # A value of 1 in the matrix indicates the species is above the threshold to react or participate in those reactions
        self.unimolecular_threshold = None
        self.bimolecular_threshold = None
        self.trimolecular_threshold = None

        # number of allowed model resurrections
        self.retry = 0
        self.max_retries = 5

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (self.__class__, (self.termination,))

    cpdef initialize_model(self, list core_species, list core_reactions, list edge_species, list edge_reactions,
                           list surface_species=None, list surface_reactions=None, list pdep_networks=None,
                           atol=1e-16, rtol=1e-8, sensitivity=False, sens_atol=1e-6, sens_rtol=1e-4,
                           filter_reactions=False, dict conditions=None):
        """
        Initialize a simulation of the reaction system using the provided
        kinetic model. You will probably want to create your own version of this
        method in the derived class; don't forget to also call the base class
        version, too.
        """
        if surface_species is None:
            surface_species = []
        if surface_reactions is None:
            surface_reactions = []

        if conditions:
            is_conc = hasattr(self, 'initial_concentrations')
            is_surf = hasattr(self, 'initial_gas_mole_fractions')
            # ToDo: I think this block is incompatible with surface.pyx catalyst reactors
            keys = conditions.keys()
            if 'T' in keys and hasattr(self, 'T'):
                self.T = Quantity(conditions['T'], 'K')
            if 'P' in keys and hasattr(self, 'P'):
                self.P = Quantity(conditions['P'], 'Pa')
            for k in keys:
                if is_conc:
                    if k in self.initial_concentrations.keys():
                        self.initial_concentrations[k] = conditions[k]  #already in SI units
                elif is_surf:
                    if k in self.initial_gas_mole_fractions.keys():
                        self.initial_gas_mole_fractions[k] = conditions[k]  #already in SI units
                    if k in self.initial_surface_coverages.keys():
                        self.initial_surface_coverages[k] = conditions[k]  #already in SI units
                else:
                    if k in self.initial_mole_fractions.keys():
                        self.initial_mole_fractions[k] = conditions[k]

        self.num_core_species = len(core_species)
        self.num_core_reactions = len(core_reactions)
        self.num_edge_species = len(edge_species)
        self.num_edge_reactions = len(edge_reactions)

        self.initiate_tolerances(atol, rtol, sensitivity, sens_atol, sens_rtol)

        pdep_networks = pdep_networks or []
        self.num_pdep_networks = len(pdep_networks)

        self.kf = np.zeros((self.num_core_reactions + self.num_edge_reactions), np.float64)
        self.kb = np.zeros_like(self.kf)
        self.Keq = np.zeros_like(self.kf)

        self.generate_species_indices(core_species, edge_species)
        self.generate_reaction_indices(core_reactions, edge_reactions)
        self.generate_reactant_product_indices(core_reactions, edge_reactions)

        self.core_species_concentrations = np.zeros((self.num_core_species), np.float64)
        self.core_species_production_rates = np.zeros((self.num_core_species), np.float64)
        self.core_species_consumption_rates = np.zeros((self.num_core_species), np.float64)
        self.core_reaction_rates = np.zeros((self.num_core_reactions), np.float64)
        self.edge_reaction_rates = np.zeros((self.num_edge_reactions), np.float64)
        self.core_species_rates = np.zeros((self.num_core_species), np.float64)
        self.edge_species_rates = np.zeros((self.num_edge_species), np.float64)
        self.network_leak_rates = np.zeros((self.num_pdep_networks), np.float64)
        self.max_network_leak_rate_ratios = np.zeros((len(self.prunable_networks)), np.float64)
        self.sensitivity_coefficients = np.zeros((self.num_core_species, self.num_core_reactions), np.float64)
        self.unimolecular_threshold = np.zeros((self.num_core_species), bool)
        self.bimolecular_threshold = np.zeros((self.num_core_species, self.num_core_species), bool)
        if self.trimolecular:
            self.trimolecular_threshold = np.zeros((self.num_core_species, self.num_core_species, self.num_core_species),
                                                     bool)

        surface_species, surface_reactions = self.initialize_surface(core_species, core_reactions, surface_species,
                                                                   surface_reactions)

        self.set_prunable_indices(edge_species, pdep_networks)

    def initialize_solver(self):
        DASx.initialize(self, self.t0, self.y0, self.dydt0, self.senpar, self.atol_array, self.rtol_array)

    def reset_max_edge_species_rate_ratios(self):
        """
        This function sets max_edge_species_rate_ratios back to zero
        for pruning of ranged reactors it is important to avoid doing this
        every initialization
        """
        self.max_edge_species_rate_ratios = np.zeros((len(self.prunable_species)), np.float64)

    def set_prunable_indices(self, edge_species, pdep_networks):
        cdef object spc
        cdef list temp
        temp = []
        for i, spc in enumerate(self.prunable_species):
            try:
                temp.append(edge_species.index(spc))
            except ValueError:
                self.max_edge_species_rate_ratios[i] = np.inf  #avoid pruning of species that have been moved to core

        self.prunable_species_indices = np.array(temp)

        temp = []
        for i, spc in enumerate(self.prunable_networks):
            try:
                temp.append(pdep_networks.index(spc))
            except:
                self.max_network_leak_rate_ratios[i] = np.inf  #avoid pruning of lost networks

        self.prunable_network_indices = np.array(temp)

    @cython.boundscheck(False)
    cpdef initialize_surface(self, list core_species, list core_reactions, list surface_species, list surface_reactions):
        """
        removes surface_species and surface_reactions from  until they are self consistent: 
            1) every reaction has one species in the surface
            2) every species participates in a surface reaction
        """
        cdef np.ndarray[np.int_t, ndim=2] product_indices, reactant_indices
        cdef list surface_species_indices, surface_reaction_indices, remove_inds
        cdef set possible_species_indices
        cdef int i, j
        cdef bool not_in_surface
        cdef object obj

        logging.debug('Initializing surface...')

        product_indices = self.product_indices
        reactant_indices = self.reactant_indices

        surface_species_indices = []
        surface_reaction_indices = []
        possible_species_indices = set()
        remove_inds = []

        for obj in surface_species:
            surface_species_indices.append(core_species.index(obj))

        for obj in surface_reactions:
            surface_reaction_indices.append(core_reactions.index(obj))

        for i in surface_reaction_indices:  #remove surface reactions whose species have been moved to the bulk core
            not_in_surface = True
            for j in product_indices[i]:
                possible_species_indices.add(j)
                if j in surface_species_indices:
                    not_in_surface = False
            for j in reactant_indices[i]:
                possible_species_indices.add(j)
                if j in surface_species_indices:
                    not_in_surface = False
            if not_in_surface:
                logging.info('removing disconnected reaction from surface: {0}'.format(str(core_reactions[i])))
                surface_reactions.remove(core_reactions[i])
                surface_reaction_indices.remove(i)

        possible_species_indices -= {-1}  #remove the -1 indexes that indicate there is no third/second reactant/product

        for i in surface_species_indices:  #remove species without reactions in the surface
            if not (i in possible_species_indices):
                logging.info('removing disconnected species from surface: {0}'.format(core_species[i].label))
                remove_inds.append(i)

        for i in remove_inds:
            surface_species.remove(core_species[i])
            surface_species_indices.remove(i)

        self.surface_species_indices = np.array(surface_species_indices, dtype=np.int)
        self.surface_reaction_indices = np.array(surface_reaction_indices, dtype=np.int)

        self.valid_layering_indices = self.get_layering_indices()

        surface_species = [core_species[i] for i in surface_species_indices]
        surface_reactions = [core_reactions[i] for i in surface_reaction_indices]

        logging.debug('Surface initialization complete')

        return surface_species, surface_reactions

    def initiate_tolerances(self, atol=1e-16, rtol=1e-8, sensitivity=False, sens_atol=1e-6, sens_rtol=1e-4):
        """
        Computes the number of differential equations and initializes the tolerance arrays.        
        """
        # Compute number of equations    
        if sensitivity:
            # Set DASPK sensitivity analysis to ON
            self.sensitivity = True
            # Compute number of variables
            self.neq = self.num_core_species * (self.num_core_reactions + self.num_core_species + 1)

            self.atol_array = np.ones(self.neq, np.float64) * sens_atol
            self.atol_array[:self.num_core_species] = atol

            self.rtol_array = np.ones(self.neq, np.float64) * sens_rtol
            self.rtol_array[:self.num_core_species] = rtol

            self.senpar = np.zeros(self.num_core_reactions + self.num_core_species, np.float64)

        else:
            self.neq = self.num_core_species

            self.atol_array = np.ones(self.neq, np.float64) * atol
            self.rtol_array = np.ones(self.neq, np.float64) * rtol

            self.senpar = np.zeros(self.num_core_reactions, np.float64)

    def get_species_index(self, spc):
        """
        Retrieves the index that is associated with the parameter species
        from the species index dictionary.
        """
        return self.species_index[spc]

    def generate_reactant_product_indices(self, core_reactions, edge_reactions):
        """
        Creates a matrix for the reactants and products.
        """

        self.reactant_indices = -np.ones((self.num_core_reactions + self.num_edge_reactions, 3), np.int)
        self.product_indices = -np.ones_like(self.reactant_indices)

        for rxn in itertools.chain(core_reactions, edge_reactions):
            j = self.reaction_index[rxn]
            for l, spec in enumerate(rxn.reactants):
                i = self.get_species_index(spec)
                self.reactant_indices[j, l] = i
            for l, spec in enumerate(rxn.products):
                i = self.get_species_index(spec)
                self.product_indices[j, l] = i

    def generate_species_indices(self, core_species, edge_species):
        """
        Assign an index to each species (core first, then edge) and 
        store the (species, index) pair in a dictionary.
        """

        for index, spec in enumerate(itertools.chain(core_species, edge_species)):
            self.species_index[spec] = index

    def generate_reaction_indices(self, core_reactions, edge_reactions):
        """
        Assign an index to each reaction (core first, then edge) and 
        store the (reaction, index) pair in a dictionary.
        """

        for index, rxn in enumerate(itertools.chain(core_reactions, edge_reactions)):
            self.reaction_index[rxn] = index

    def set_initial_conditions(self):
        """
        Sets the common initial conditions of the rate equations that 
        represent the reaction system.

        - Sets the initial time of the reaction system to 0
        - Initializes the species moles to a n x 1 array with zeros
        """

        self.t0 = 0.0

        self.y0 = np.zeros(self.neq, np.float64)

    def set_initial_reaction_thresholds(self):

        # Set unimolecular and bimolecular thresholds as true for any concentrations greater than 0
        num_core_species = len(self.core_species_concentrations)
        for i in range(num_core_species):
            if self.core_species_concentrations[i] > 0:
                self.unimolecular_threshold[i] = True
        for i in range(num_core_species):
            for j in range(i, num_core_species):
                if self.core_species_concentrations[i] > 0 and self.core_species_concentrations[j] > 0:
                    self.bimolecular_threshold[i, j] = True
        if self.trimolecular:
            for i in range(num_core_species):
                for j in range(i, num_core_species):
                    for k in range(j, num_core_species):
                        if (self.core_species_concentrations[i] > 0
                                and self.core_species_concentrations[j] > 0
                                and self.core_species_concentrations[k] > 0):
                            self.trimolecular_threshold[i, j, k] = True

    def set_initial_derivative(self):
        """
        Sets the derivative of the species moles with respect to the independent variable (time)
        equal to the residual.
        """
        self.dydt0 = - self.residual(self.t0, self.y0, np.zeros(self.neq, np.float64), self.senpar)[0]

    def compute_network_variables(self, pdep_networks=None):
        """
        Initialize the arrays containing network information:

        - NetworkLeakCoefficients is a n x 1 array with 
            n the number of pressure-dependent networks.
        - NetworkIndices is a n x 3 matrix with 
            n the number of pressure-dependent networks and 
            3 the maximum number of molecules allowed in either the reactant or
            product side of a reaction. 
        """

        pdep_networks = pdep_networks or []

        self.network_indices = -np.ones((self.num_pdep_networks, 3), np.int)
        self.network_leak_coefficients = np.zeros((self.num_pdep_networks), np.float64)

        for j, network in enumerate(pdep_networks):
            self.network_leak_coefficients[j] = network.get_leak_coefficient(self.T.value_si, self.P.value_si)
            for k, spec in enumerate(network.source):
                try:
                    i = self.get_species_index(spec)
                except KeyError:
                    self.network_indices[j, :] = [-1,-1,-1]
                    break
                if i >= self.num_core_species: #an edge species is in source
                    self.network_indices[j, :] = [-1,-1,-1]
                    break
                self.network_indices[j, k] = i

    @cython.boundscheck(False)
    cpdef get_layering_indices(self):
        """
        determines the edge reaction indices that indicate reactions that are valid for movement from
        edge to surface based on the layering constraint
        """

        cdef np.ndarray[np.int_t, ndim=1] surf_species_indices
        cdef int i, j, k, index, num_core_species, num_core_reactions, num_edge_reactions
        cdef list valid_indices

        surf_species_indices = self.surface_species_indices
        num_core_species = self.num_core_species
        num_core_reactions = self.num_core_reactions
        num_edge_reactions = self.num_edge_reactions
        product_indices = self.product_indices
        reactant_indices = self.reactant_indices

        valid_indices = []

        for index in range(num_edge_reactions):
            for j in product_indices[index + num_core_reactions]:
                if j in surf_species_indices or j >= num_core_species:
                    break
            else:
                valid_indices.append(index)
                continue
            for j in reactant_indices[index + num_core_reactions]:
                if j in surf_species_indices or j >= num_core_species:
                    break
            else:
                valid_indices.append(index)

        return np.array(valid_indices)

    cpdef add_reactions_to_surface(self, list new_surface_reactions, list new_surface_reaction_inds, list surface_species,
                                   list surface_reactions, list edge_species):
        """
        moves new surface reactions to the surface
        done after the while loop before the simulate call ends
        """
        cdef object srxn
        cdef int k, sind, num_core_species, i, num_core_reactions
        cdef np.ndarray[np.int_t, ndim=2] product_indices, reactant_indices

        product_indices = self.product_indices
        reactant_indices = self.reactant_indices
        num_core_species = self.num_core_species
        num_core_reactions = self.num_core_reactions

        for k in range(len(new_surface_reactions)):
            srxn = new_surface_reactions[k]
            sind = new_surface_reaction_inds[k]
            surface_reactions.append(srxn)  #add to surface trackers

            for i in product_indices[sind + num_core_reactions]:
                if i >= num_core_species:
                    surface_species.append(edge_species[i - num_core_species])
            for i in reactant_indices[sind + num_core_reactions]:
                if i >= num_core_species:
                    surface_species.append(edge_species[i - num_core_species])

        return surface_species, surface_reactions

    @cython.boundscheck(False)
    cpdef simulate(self, list core_species, list core_reactions, list edge_species,
                   list edge_reactions, list surface_species, list surface_reactions,
                   list pdep_networks=None, bool prune=False, bool sensitivity=False, list sens_worksheet=None,
                   object model_settings=None, object simulator_settings=None, dict conditions=None):
        """
        Simulate the reaction system with the provided reaction model,
        consisting of lists of core species, core reactions, edge species, and
        edge reactions. As the simulation proceeds the system is monitored for
        validity. If the model becomes invalid (e.g. due to an excessively
        large edge flux), the simulation is interrupted and the object causing
        the model to be invalid is returned. If the simulation completes to
        the desired termination criteria and the model remains valid throughout,
        ``None`` is returned.
        """
        cdef double tol_keep_in_edge, tol_move_to_core, tol_move_edge_reaction_to_core, tol_interrupt_simulation
        cdef double tol_move_edge_reaction_to_core_interrupt, tol_move_edge_reaction_to_surface
        cdef double tol_move_surface_species_to_core, tol_move_surface_reaction_to_core
        cdef double tol_move_edge_reaction_to_surface_interrupt, b_num
        cdef bool ignore_overall_flux_criterion, filter_reactions
        cdef double atol, rtol, sens_atol, sens_rtol
        cdef dict species_index
        cdef list row, temp_surface_objects
        cdef list sorted_inds, temp_new_objects, temp_new_object_inds, temp_new_object_vals, tempInds
        cdef int index, spc_index, max_species_index, max_network_index
        cdef int num_core_species, num_edge_species, num_pdep_networks, num_core_reactions
        cdef double step_time, char_rate, max_species_rate, max_network_rate, maxEdgeReactionAccum, stdan
        cdef np.ndarray[np.float64_t, ndim=1] y0  # Vector containing the number of moles of each species
        cdef np.ndarray[np.float64_t, ndim=1] core_species_rates, edge_species_rates, network_leak_rates
        cdef np.ndarray[np.float64_t, ndim=1] core_species_production_rates, core_species_consumption_rates, total_div_accum_nums
        cdef np.ndarray[np.float64_t, ndim=1] max_edge_species_rate_ratios, max_network_leak_rate_ratios
        cdef bint terminated
        cdef object max_species, max_network
        cdef int i, j, k
        cdef np.float64_t conversion
        cdef np.ndarray[np.float64_t, ndim=1] surface_species_production, surface_species_consumption, branching_nums
        cdef np.ndarray[np.float64_t, ndim=1] surface_total_div_accum_nums, surface_species_rate_ratios
        cdef np.ndarray[np.float64_t, ndim=1] forward_rate_coefficients, core_species_concentrations
        cdef double prev_time, total_moles, c, volume, RTP, max_char_rate, br, rr
        cdef double unimolecular_threshold_val, bimolecular_threshold_val, trimolecular_threshold_val
        cdef bool useDynamicsTemp, first_time, use_dynamics, terminate_at_max_objects, schanged, invalid_objects_print_boolean
        cdef np.ndarray[np.float64_t, ndim=1] edge_reaction_rates
        cdef double reaction_rate, production, consumption
        cdef np.ndarray[np.int_t, ndim=1] surface_species_indices, surface_reaction_indices
        # cython declations for sensitivity analysis
        cdef np.ndarray[np.int_t, ndim=1] sens_species_indices, reactant_side, product_side
        cdef np.ndarray[np.float64_t, ndim=1] mole_sens, dVdk, norm_sens
        cdef list time_array, norm_sens_array, new_surface_reactions, new_surface_reaction_inds, new_objects, new_object_inds



        zero_production = False
        zero_consumption = False
        pdep_networks = pdep_networks or []

        schanged = False

        num_core_species = len(core_species)
        num_edge_species = len(edge_species)
        num_pdep_networks = len(pdep_networks)
        num_core_reactions = len(core_reactions)

        assert set(core_reactions) >= set(surface_reactions), 'given surface reactions are not a subset of core reactions'
        assert set(core_species) >= set(surface_species), 'given surface species are not a subset of core species'

        tol_keep_in_edge = model_settings.tol_keep_in_edge if prune else 0
        tol_move_to_core = model_settings.tol_move_to_core
        tol_move_edge_reaction_to_core = model_settings.tol_move_edge_rxn_to_core
        tol_interrupt_simulation = model_settings.tol_interrupt_simulation
        tol_move_edge_reaction_to_core_interrupt = model_settings.tol_move_edge_rxn_to_core
        tol_move_edge_reaction_to_surface = model_settings.tol_move_edge_rxn_to_surface
        tol_move_surface_species_to_core = model_settings.tol_move_surface_spc_to_core
        tol_move_surface_reaction_to_core = model_settings.tol_move_surface_rxn_to_core
        tol_move_edge_reaction_to_surface_interrupt = model_settings.tol_move_edge_rxn_to_surface_interrupt
        ignore_overall_flux_criterion = model_settings.ignore_overall_flux_criterion
        atol = simulator_settings.atol
        rtol = simulator_settings.rtol
        sens_atol = simulator_settings.sens_atol
        sens_rtol = simulator_settings.sens_rtol
        filter_reactions = model_settings.filter_reactions
        max_num_objs_per_iter = model_settings.max_num_objects_per_iter

        if model_settings.tol_branch_rxn_to_core != 0.0:
            branch_factor = 1.0 / model_settings.tol_branch_rxn_to_core
            br_max = model_settings.branching_ratio_max
            branching_index = model_settings.branching_index
        else:
            branch_factor = 0.0

        #if not pruning always terminate at max objects, otherwise only do so if terminate_at_max_objects=True
        terminate_at_max_objects = True if not prune else model_settings.terminate_at_max_objects

        dynamics_time_scale = model_settings.dynamics_time_scale

        use_dynamics = not (tol_move_edge_reaction_to_core == np.inf and tol_move_edge_reaction_to_surface == np.inf)

        species_index = {}
        for index, spec in enumerate(core_species):
            species_index[spec] = index

        self.initialize_model(core_species, core_reactions,
                              edge_species, edge_reactions,
                              surface_species, surface_reactions,
                              pdep_networks, atol, rtol, sensitivity,
                              sens_atol, sens_rtol,
                              filter_reactions, conditions)

        prunable_species_indices = self.prunable_species_indices
        prunable_network_indices = self.prunable_network_indices

        surface_species_indices = self.surface_species_indices
        surface_reaction_indices = self.surface_reaction_indices

        # the product of the ratios between accumulation numbers
        # with and without a given reaction for products and reactants
        total_div_accum_nums = None
        branching_nums = None

        invalid_objects = []
        new_surface_reactions = []
        new_surface_reaction_inds = []
        terminated = False
        max_species_index = -1
        max_species = None
        max_species_rate = 0.0
        max_network_index = -1
        max_network = None
        max_network_rate = 0.0
        iteration = 0
        conversion = 0.0
        max_char_rate = 0.0

        max_edge_species_rate_ratios = self.max_edge_species_rate_ratios
        max_network_leak_rate_ratios = self.max_network_leak_rate_ratios
        forward_rate_coefficients = self.kf
        unimolecular_threshold = self.unimolecular_threshold
        bimolecular_threshold = self.bimolecular_threshold
        trimolecular_threshold = self.trimolecular_threshold

        # Copy the initial conditions to use in evaluating conversions
        y0 = self.y.copy()

        # a list with the time, Volume, number of moles of core species
        self.snapshots = []

        if sensitivity:
            time_array = []
            norm_sens_array = [[] for spec in self.sensitive_species]
            RTP = constants.R * self.T.value_si / self.P.value_si
            # identify sensitive species indices
            sens_species_indices = np.array([species_index[spec] for spec in self.sensitive_species],
                                               np.int)  # index within core_species list of the sensitive species

        step_time = 1e-12
        prev_time = self.t

        first_time = True

        invalid_objects_print_boolean = True  
        while not terminated:
            # Integrate forward in time by one time step

            if not first_time:
                try:
                    self.step(step_time)
                    if np.isnan(self.y).any():
                        raise DASxError("nans in moles")
                except DASxError as e:
                    self.retry += 1
                    logging.error("Trying to step from time {0} to {1} resulted in a solver (DASPK) error: "
                                  "{2!s}".format(prev_time, step_time, e))

                    if self.retry > self.max_retries:
                        raise RuntimeError(
                            'Reached the maximum permitted number of retries due to solver errors ({0}), ending simulation'.format(self.max_retries))


                    logging.info('Resurrecting Model...')

                    conversion = 0.0
                    for term in self.termination:
                        if isinstance(term, TerminationConversion):
                            index = species_index[term.species]
                            conversion = 1 - (y_core_species[index] / y0[index])

                    if invalid_objects == []:
                        
                        #species flux criterion
                        if len(edge_species_rate_ratios) > 0:
                            ind = np.argmax(edge_species_rate_ratios)
                            obj = edge_species[ind]
                            logging.info('At time {0:10.4e} s, species {1} at rate ratio {2} was added to model core '
                                            'in model resurrection process.'.format(self.t, obj, edge_species_rates[ind]))
                            invalid_objects.append(obj)


                        if total_div_accum_nums and len(total_div_accum_nums) > 0:  #if dynamics data available
                            ind = np.argmax(total_div_accum_nums)
                            obj = edge_reactions[ind]
                            logging.info('At time {0:10.4e} s, Reaction {1} at dynamics number {2} was added to model core '
                                            'in model resurrection process.'.format(self.t, obj, total_div_accum_nums[ind]))
                            invalid_objects.append(obj)


                        if pdep_networks != [] and network_leak_rate_ratios != []:
                            ind = np.argmax(network_leak_rate_ratios)
                            obj = pdep_networks[ind]
                            logging.info('At time {0:10.4e} s, PDepNetwork #{1:d} at network leak rate {2} '
                                         'was sent for exploring during model resurrection process'
                                         ''.format(self.t, obj.index, network_leak_rate_ratios[ind]))
                            invalid_objects.append(obj)
                            
                        logging.info('Retry # {0} of {1}'.format(
                            self.retry, self.max_retries))

                    if invalid_objects != []:
                        return False, True, invalid_objects, surface_species, surface_reactions, self.t, conversion
                    else:
                        logging.error('Model Resurrection has failed')
                        logging.error("Core species names: {!r}".format([get_species_identifier(s) for s in core_species]))
                        logging.error("Core species moles: {!r}".format(self.y[:num_core_species]))
                        logging.error("Volume: {!r}".format(self.V))
                        logging.error("Core species net rates: {!r}".format(self.core_species_rates))
                        logging.error("Edge species net rates: {!r}".format(self.edge_species_rates))
                        logging.error("Network leak rates: {!r}".format(self.network_leak_rates))
                        raise ValueError('invalid_objects could not be filled during resurrection process')

            y_core_species = self.y[:num_core_species]
            total_moles = np.sum(y_core_species)
            if sensitivity:
                time_array.append(self.t)
                mole_sens = self.y[num_core_species:]
                volume = self.V

                dVdk = np.zeros(num_core_reactions + num_core_species, np.float64)
                if not self.constant_volume:
                    for j in range(num_core_reactions + num_core_species):
                        dVdk[j] = np.sum(mole_sens[j * num_core_species:(j + 1) * num_core_species]) * RTP  # Contains [ dV_dk and dV_dG ]
                for i in range(len(self.sensitive_species)):
                    norm_sens = np.zeros(num_core_reactions + num_core_species, np.float64)
                    c = self.core_species_concentrations[sens_species_indices[i]]
                    if c != 0:
                        for j in range(num_core_reactions):
                            norm_sens[j] = 1/volume * (mole_sens[j*num_core_species+sens_species_indices[i]]-c*dVdk[j]) * forward_rate_coefficients[j]/c
                        for j in range(num_core_reactions,num_core_reactions+num_core_species):
                            # no normalization against dG, conversion to kcal/mol units
                            norm_sens[j] = 1/volume * (mole_sens[j*num_core_species+sens_species_indices[i]]-c*dVdk[j]) / c * 4184
                    norm_sens_array[i].append(norm_sens)

            snapshot = [self.t, self.V]
            snapshot.extend(y_core_species)
            self.snapshots.append(snapshot)

            # Get the characteristic flux
            char_rate = sqrt(np.sum(self.core_species_rates * self.core_species_rates))

            if char_rate > max_char_rate:
                max_char_rate = char_rate

            core_species_rates = np.abs(self.core_species_rates)
            edge_reaction_rates = self.edge_reaction_rates
            core_species_consumption_rates = self.core_species_consumption_rates
            core_species_production_rates = self.core_species_production_rates
            edge_species_rates = np.abs(self.edge_species_rates)
            network_leak_rates = np.abs(self.network_leak_rates)
            core_species_rate_ratios = np.abs(self.core_species_rates / char_rate)
            edge_species_rate_ratios = np.abs(self.edge_species_rates / char_rate)
            network_leak_rate_ratios = np.abs(self.network_leak_rates / char_rate)
            num_edge_reactions = self.num_edge_reactions
            core_reaction_rates = self.core_reaction_rates
            product_indices = self.product_indices
            reactant_indices = self.reactant_indices
            core_species_concentrations = self.core_species_concentrations

            # Update the maximum species rate and maximum network leak rate arrays
            for i, index in enumerate(prunable_species_indices):
                if max_edge_species_rate_ratios[i] < edge_species_rate_ratios[index]:
                    max_edge_species_rate_ratios[i] = edge_species_rate_ratios[index]
            for i, index in enumerate(prunable_network_indices):
                if max_network_leak_rate_ratios[i] < network_leak_rate_ratios[index]:
                    max_network_leak_rate_ratios[i] = network_leak_rate_ratios[index]

            if char_rate == 0 and len(edge_species_rates) > 0:  # this deals with the case when there is no flux in the system
                max_species_index = np.argmax(edge_species_rates)
                max_species = edge_species[max_species_index]
                max_species_rate = edge_species_rates[max_species_index]
                logging.info('At time {0:10.4e} s, species {1} was added to model core to avoid '
                             'singularity'.format(self.t, max_species))
                invalid_objects.append(max_species)
                break

            if branch_factor != 0.0 and not first_time:
                ######################################################
                # Calculation of branching numbers for edge reactions#
                ######################################################
                branching_nums = np.zeros(num_edge_reactions)
                for index in range(num_edge_reactions):
                    reaction_rate = edge_reaction_rates[index]

                    if reaction_rate > 0:
                        reactant_side = self.reactant_indices[index + num_core_reactions, :]
                        product_side = self.product_indices[index + num_core_reactions, :]
                    else:
                        reactant_side = self.product_indices[index + num_core_reactions, :]
                        product_side = self.reactant_indices[index + num_core_reactions, :]

                    mults = []
                    for i in product_side:
                        if i == -1:
                            continue
                        elif i < num_core_species:
                            mults.append(core_species[i].molecule[0].multiplicity)
                        else:
                            mults.append(edge_species[i - num_core_species].molecule[0].multiplicity)

                    if max(mults) > 2:
                        continue

                    for spc_index in reactant_side:
                        if spc_index != -1 and spc_index < num_core_species:
                            if core_species[spc_index].molecule[0].multiplicity != 2:
                                continue
                            consumption = core_species_consumption_rates[spc_index]
                            if consumption != 0:  #if consumption = 0 ignore species
                                br = reaction_rate / consumption
                                rr = core_species_rate_ratios[spc_index]
                                if br > br_max:
                                    br = br_max
                                b_num = branch_factor * br * rr ** branching_index
                                if b_num > branching_nums[index]:
                                    branching_nums[index] = b_num

            if use_dynamics and not first_time and self.t >= dynamics_time_scale:
                #######################################################
                # Calculation of dynamics criterion for edge reactions#
                #######################################################

                total_div_accum_nums = np.ones(num_edge_reactions)
                for index in range(num_edge_reactions):
                    reaction_rate = edge_reaction_rates[index]
                    for spc_index in self.reactant_indices[index + num_core_reactions, :]:
                        if spc_index != -1 and spc_index < num_core_species:
                            consumption = core_species_consumption_rates[spc_index]
                            if consumption != 0:  #if consumption = 0 ignore species
                                total_div_accum_nums[index] *= (reaction_rate + consumption) / consumption

                    for spc_index in self.product_indices[index + num_core_reactions, :]:
                        if spc_index != -1 and spc_index < num_core_species:
                            production = core_species_production_rates[spc_index]
                            if production != 0:  #if production = 0 ignore species
                                total_div_accum_nums[index] *= (reaction_rate + production) / production

                total_div_ln_accum_nums = np.log(total_div_accum_nums)

                surface_species_indices = self.surface_species_indices
                surface_reaction_indices = self.surface_reaction_indices

                ##########################################################
                # Calculation of dynamics criterion for surface reactions#
                ##########################################################

                surface_total_div_accum_nums = np.ones(len(surface_reaction_indices))

                for i in range(len(surface_reaction_indices)):
                    index = surface_reaction_indices[i]
                    reaction_rate = core_reaction_rates[index]
                    for spc_index in reactant_indices[index, :]:
                        if spc_index != -1 and spc_index < num_core_species and not (spc_index in surface_species_indices):
                            consumption = core_species_consumption_rates[spc_index]
                            if consumption != 0:  #if consumption=0 ignore species
                                if abs(abs(consumption) - abs(reaction_rate)) < atol:
                                    surface_total_div_accum_nums[i] = np.inf
                                elif reaction_rate > 0:
                                    surface_total_div_accum_nums[i] *= consumption / (consumption - reaction_rate)
                                else:
                                    surface_total_div_accum_nums[i] *= (consumption - reaction_rate) / consumption

                    for spc_index in product_indices[index, :]:
                        if spc_index != -1 and spc_index < num_core_species and not (spc_index in surface_species_indices):
                            production = core_species_production_rates[spc_index]
                            if production != 0:  #if production = 0 ignore species
                                if abs(abs(production) - abs(reaction_rate)) < atol:
                                    surface_total_div_accum_nums[i] = np.inf
                                elif reaction_rate > 0:
                                    surface_total_div_accum_nums[i] *= production / (production - reaction_rate)
                                else:
                                    surface_total_div_accum_nums[i] *= (production - reaction_rate) / production

                surface_total_div_accum_nums = np.log(surface_total_div_accum_nums)

                ###############################################
                # Move objects from surface to core on-the-fly#
                ###############################################

                surface_objects = []
                surface_object_indices = []

                #Determination of reactions moving from surface to core on-the-fly

                for ind, stdan in enumerate(surface_total_div_accum_nums):
                    if stdan > tol_move_surface_reaction_to_core:
                        sind = surface_reaction_indices[ind]
                        surface_object_indices.append(sind)
                        surface_objects.append(core_reactions[sind])

                #Determination of species moving from surface to core on-the-fly

                surface_species_rate_ratios = np.zeros(len(surface_species_indices))
                surface_species_production = core_species_production_rates[surface_species_indices]
                surface_species_consumption = core_species_consumption_rates[surface_species_indices]

                for i in range(len(surface_species_indices)):
                    rr = max(abs(surface_species_production[i]), abs(surface_species_consumption[i])) / char_rate
                    surface_species_rate_ratios[i] = rr

                    if rr > tol_move_surface_species_to_core:
                        sind = surface_species_indices[i]
                        surface_object_indices.append(sind)
                        surface_objects.append(core_species[sind])

                #process objects to be moved from surface to core

                if len(surface_objects) > 0:
                    schanged = True
                    for ind, obj in enumerate(surface_objects):
                        if isinstance(obj, Reaction):
                            logging.info('Moving reaction {0} from surface to core'.format(obj))
                            surface_reactions.remove(obj)
                        elif isinstance(obj, Species):
                            logging.info('Moving species {0} from surface to core'.format(obj))
                            surface_species.remove(obj)
                        else:
                            raise ValueError
                    surface_species, surface_reactions = self.initialize_surface(core_species, core_reactions,
                                                                               surface_species, surface_reactions)
                    logging.info('Surface now has {0} Species and {1} Reactions'
                                 ''.format(len(self.surface_species_indices), len(self.surface_reaction_indices)))

            if filter_reactions:
                # Calculate thresholds for reactions
                (unimolecular_threshold_rate_constant,
                 bimolecular_threshold_rate_constant,
                 trimolecular_threshold_rate_constant) = self.get_threshold_rate_constants(model_settings)
                unimolecular_threshold_val = tol_move_to_core * char_rate / unimolecular_threshold_rate_constant
                bimolecular_threshold_val = tol_move_to_core * char_rate / bimolecular_threshold_rate_constant
                trimolecular_threshold_val = tol_move_to_core * char_rate / trimolecular_threshold_rate_constant
                for i in range(num_core_species):
                    if not unimolecular_threshold[i]:
                        # Check if core species concentration has gone above threshold for unimolecular reaction
                        if core_species_concentrations[i] > unimolecular_threshold_val:
                            unimolecular_threshold[i] = True
                for i in range(num_core_species):
                    for j in range(i, num_core_species):
                        if not bimolecular_threshold[i, j]:
                            if core_species_concentrations[i] * core_species_concentrations[j] > bimolecular_threshold_val:
                                bimolecular_threshold[i, j] = True
                if self.trimolecular:
                    for i in range(num_core_species):
                        for j in range(i, num_core_species):
                            for k in range(j, num_core_species):
                                if not trimolecular_threshold[i, j, k]:
                                    if (core_species_concentrations[i] *
                                            core_species_concentrations[j] *
                                            core_species_concentrations[k]
                                            > trimolecular_threshold_val):
                                        trimolecular_threshold[i, j, k] = True

            ###############################################################################
            # Movement from edge to core or surface processing and interrupt determination#
            ###############################################################################

            new_object_inds = []
            new_objects = []
            new_object_vals = []
            new_object_type = []

            temp_new_objects = []
            temp_new_object_inds = []
            temp_new_object_vals = []
            temp_new_object_type = []

            new_surface_rxn_inds = []
            interrupt = False

            #movement of species to core based on rate ratios

            if not ignore_overall_flux_criterion:
                for ind, obj in enumerate(edge_species):
                    rr = edge_species_rate_ratios[ind]
                    if rr > tol_move_to_core:
                        if not (obj in new_objects or obj in invalid_objects):
                            temp_new_objects.append(edge_species[ind])
                            temp_new_object_inds.append(ind)
                            temp_new_object_vals.append(rr)
                            temp_new_object_type.append('rr')
                    if rr > tol_interrupt_simulation:
                        logging.info('At time {0:10.4e} s, species {1} at {2} exceeded the minimum rate for simulation '
                                     'interruption of {3}'.format(self.t, obj, rr, tol_interrupt_simulation))
                        interrupt = True

                sorted_inds = np.argsort(np.array(temp_new_object_vals)).tolist()[::-1]

                new_objects.extend([temp_new_objects[q] for q in sorted_inds])
                new_object_inds.extend([temp_new_object_inds[q] for q in sorted_inds])
                new_object_vals.extend([temp_new_object_vals[q] for q in sorted_inds])
                new_object_type.extend([temp_new_object_type[q] for q in sorted_inds])

                temp_new_objects = []
                temp_new_object_inds = []
                temp_new_object_vals = []
                temp_new_object_type = []

            if branch_factor != 0.0 and not first_time:
                #movement of reactions to core based on branching number
                for ind, obj in enumerate(edge_reactions):
                    b_num = branching_nums[ind]
                    if b_num > 1:
                        if not (obj in new_objects or obj in invalid_objects):
                            temp_new_objects.append(edge_reactions[ind])
                            temp_new_object_inds.append(ind)
                            temp_new_object_vals.append(b_num)
                            temp_new_object_type.append('branching')

                sorted_inds = np.argsort(np.array(temp_new_object_vals)).tolist()[::-1]

                new_objects.extend([temp_new_objects[q] for q in sorted_inds])
                new_object_inds.extend([temp_new_object_inds[q] for q in sorted_inds])
                new_object_vals.extend([temp_new_object_vals[q] for q in sorted_inds])
                new_object_type.extend([temp_new_object_type[q] for q in sorted_inds])

                temp_new_objects = []
                temp_new_object_inds = []
                temp_new_object_vals = []
                temp_new_object_type = []

            if use_dynamics and not first_time and self.t >= dynamics_time_scale:
                #movement of reactions to core/surface based on dynamics number  
                valid_layering_indices = self.valid_layering_indices
                temp_surface_objects = []

                for ind, obj in enumerate(edge_reactions):
                    dlnaccum = total_div_ln_accum_nums[ind]
                    if dlnaccum > tol_move_edge_reaction_to_core:
                        if not (obj in new_objects or obj in invalid_objects):
                            temp_new_objects.append(edge_reactions[ind])
                            temp_new_object_inds.append(ind)
                            temp_new_object_vals.append(dlnaccum)
                            temp_new_object_type.append('dyn')
                    elif dlnaccum > tol_move_edge_reaction_to_surface and ind in valid_layering_indices:
                        if not (obj in new_objects or obj in invalid_objects):
                            temp_new_objects.append(edge_reactions[ind])
                            temp_new_object_inds.append(ind)
                            temp_new_object_vals.append(dlnaccum)
                            temp_surface_objects.append(edge_reactions[ind])
                            temp_new_object_type.append('dyn')
                    if dlnaccum > tol_move_edge_reaction_to_core_interrupt:
                        logging.info('At time {0:10.4e} s, Reaction {1} at {2} exceeded the minimum difference in '
                                     'total log(accumulation number) for simulation interruption of {3}'
                                     ''.format(self.t, obj,dlnaccum,tol_move_edge_reaction_to_core_interrupt))
                        interrupt = True

                sorted_inds = np.argsort(np.array(temp_new_object_vals)).tolist()[::-1]

                new_objects.extend([temp_new_objects[q] for q in sorted_inds])
                new_object_inds.extend([temp_new_object_inds[q] for q in sorted_inds])
                new_object_vals.extend([temp_new_object_vals[q] for q in sorted_inds])
                new_object_type.extend([temp_new_object_type[q] for q in sorted_inds])

                new_surface_rxn_inds = [new_objects.index(obj) for obj in temp_surface_objects]

                temp_new_objects = []
                temp_new_object_inds = []
                temp_new_object_vals = []
                temp_new_object_type = []

            #Determination of pdep_networks in need of exploring

            if pdep_networks:
                for ind, obj in enumerate(pdep_networks):
                    lr = network_leak_rate_ratios[ind]
                    if lr > tol_move_to_core:
                        if not (obj in new_objects or obj in invalid_objects):
                            temp_new_objects.append(pdep_networks[ind])
                            temp_new_object_inds.append(ind)
                            temp_new_object_vals.append(lr)
                            temp_new_object_type.append('pdep')
                    if lr > tol_interrupt_simulation:
                        logging.info('At time {0:10.4e} s, PDepNetwork #{1:d} at {2} exceeded the minimum rate for '
                                     'simulation interruption of {3}'.format(self.t, obj.index, lr,
                                                                             tol_interrupt_simulation))
                        interrupt = True

                sorted_inds = np.argsort(np.array(temp_new_object_vals)).tolist()[::-1]

                new_objects.extend([temp_new_objects[q] for q in sorted_inds])
                new_object_inds.extend([temp_new_object_inds[q] for q in sorted_inds])
                new_object_vals.extend([temp_new_object_vals[q] for q in sorted_inds])
                new_object_type.extend([temp_new_object_type[q] for q in sorted_inds])

                temp_new_objects = []
                temp_new_object_inds = []
                temp_new_object_vals = []
                temp_new_object_type = []

            ###########################
            #Overall Object Processing#
            ###########################

            #remove excess objects
            if len(invalid_objects) + len(new_objects) > max_num_objs_per_iter:
                if invalid_objects_print_boolean:
                    logging.info('Exceeded max number of objects...removing excess objects')
                    invalid_objects_print_boolean = False
                num = max_num_objs_per_iter - len(invalid_objects)
                new_objects = new_objects[:num]
                new_object_inds = new_object_inds[:num]
                new_object_vals = new_object_vals[:num]
                new_object_type = new_object_type[:num]

            if terminate_at_max_objects and len(invalid_objects) + len(new_objects) >= max_num_objs_per_iter:
                logging.info('Reached max number of objects...preparing to terminate')
                interrupt = True

            if new_objects != []:
                for i, obj in enumerate(new_objects):  # log everything and add reactions to surface
                    val = new_object_vals[i]
                    ind = new_object_inds[i]
                    if isinstance(obj, Species):
                        logging.info('At time {0:10.4e} s, species {1} at rate ratio {2} exceeded the minimum rate for '
                                     'moving to model core of {3}'.format(self.t, obj, val, tol_move_to_core))
                    elif isinstance(obj, Reaction):
                        if new_object_type[i] == 'dyn':
                            if i in new_surface_rxn_inds:
                                logging.info('At time {0:10.4e} s, Reaction {1} at {2} exceeded the minimum difference '
                                             'in total log(accumulation number) for moving to model surface of {3}'
                                             ''.format(self.t, obj, val, tol_move_edge_reaction_to_surface))
                                new_surface_reactions.append(obj)
                                new_surface_reaction_inds.append(ind)
                            else:
                                logging.info('At time {0:10.4e} s, Reaction {1} at {2} exceeded the minimum difference '
                                             'in total log(accumulation number) for moving to model core of {3}'
                                             ''.format(self.t, obj, val, tol_move_edge_reaction_to_core))
                        elif new_object_type[i] == 'branching':
                            logging.info('At time {0:10.4e} s, Reaction {1} at a branching number of {2} exceeded the '
                                         'threshold of 1 for moving to model core'.format(self.t, obj, val))
                    else:
                        logging.info('At time {0:10.4e} s, PDepNetwork #{1:d} at {2} exceeded the minimum rate for '
                                     'exploring of {3}'.format(self.t, obj.index, val, tol_move_to_core))

                invalid_objects += new_objects

            if schanged:  # reinitialize surface
                surface_species, surface_reactions = self.initialize_surface(core_species, core_reactions,
                                                                             surface_species, surface_reactions)
                schanged = False

            if first_time:  # turn off first_time
                first_time = False

            if interrupt:  # breaks while loop terminating iterations
                logging.info('terminating simulation due to interrupt...')
                break

            # Finish simulation if any of the termination criteria are satisfied
            for term in self.termination:
                if isinstance(term, TerminationTime):
                    if self.t > term.time.value_si:
                        terminated = True
                        logging.info('At time {0:10.4e} s, reached target termination time.'.format(term.time.value_si))
                        self.log_conversions(species_index, y0)
                        break
                elif isinstance(term, TerminationConversion):
                    index = species_index[term.species]
                    conversion = 1 - (y_core_species[index] / y0[index])
                    if 1 - (y_core_species[index] / y0[index]) > term.conversion:
                        terminated = True
                        logging.info('At time {0:10.4e} s, reached target termination conversion: {1:f} of '
                                     '{2}'.format(self.t,term.conversion,term.species))
                        self.log_conversions(species_index, y0)
                        break
                elif isinstance(term, TerminationRateRatio):
                    if max_char_rate != 0.0 and char_rate / max_char_rate < term.ratio:
                        terminated = True
                        logging.info('At time {0:10.4e} s, reached target termination RateRatio: '
                                     '{1}'.format(self.t,char_rate/max_char_rate))
                        self.log_conversions(species_index, y0)

            # Increment destination step time if necessary
            if self.t >= 0.9999 * step_time:
                step_time *= 10.0

        # Change surface species and reactions based on what will be added to the surface
        surface_species, surface_reactions = self.add_reactions_to_surface(new_surface_reactions,
                                                                           new_surface_reaction_inds,
                                                                           surface_species,
                                                                           surface_reactions,
                                                                           edge_species)

        # notify reaction system listeners
        self.notify()

        if sensitivity:
            for i in range(len(self.sensitive_species)):
                with open(sens_worksheet[i], 'w') as outfile:
                    worksheet = csv.writer(outfile)
                    reactions_above_threshold = []
                    for j in range(num_core_reactions + num_core_species):
                        for k in range(len(time_array)):
                            if abs(norm_sens_array[i][k][j]) > self.sensitivity_threshold:
                                reactions_above_threshold.append(j)
                                break
                    species_name = get_species_identifier(self.sensitive_species[i])
                    headers = ['Time (s)']
                    headers.extend(['dln[{0}]/dln[k{1}]: {2}'.format(species_name, j + 1, core_reactions[j].to_chemkin(kinetics=False)) if j < num_core_reactions
                                    else 'dln[{0}]/dG[{1}]'.format(species_name, get_species_identifier(core_species[j - num_core_reactions])) for j in reactions_above_threshold])
                    worksheet.writerow(headers)

                    for k in range(len(time_array)):
                        row = [time_array[k]]
                        row.extend([norm_sens_array[i][k][j] for j in reactions_above_threshold])
                        worksheet.writerow(row)

        self.max_edge_species_rate_ratios = max_edge_species_rate_ratios
        self.max_network_leak_rate_ratios = max_network_leak_rate_ratios

        self.unimolecular_threshold = unimolecular_threshold
        self.bimolecular_threshold = bimolecular_threshold
        self.trimolecular_threshold = trimolecular_threshold

        # Return the invalid object (if the simulation was invalid) or None
        # (if the simulation was valid)
        return terminated, False, invalid_objects, surface_species, surface_reactions, self.t, conversion

    cpdef log_rates(self, double char_rate, object species, double species_rate, double max_dif_ln_accum_num, object network,
                    double network_rate):
        """
        Log information about the current maximum species and network rates.
        """
        logging.info('    Characteristic rate: {0:10.4e} mol/m^3*s'.format(char_rate))
        logging.info('    Dynamics number:  {0:10.4e}'.format(max_dif_ln_accum_num))
        if char_rate == 0.0:
            logging.info('    {0} rate: {1:10.4e} mol/m^3*s'.format(species, species_rate))
            if network is not None:
                logging.info('    PDepNetwork #{0:d} leak rate: {1:10.4e} mol/m^3*s'.format(network.index, network_rate))
        else:
            logging.info('    {0} rate: {1:10.4e} mol/m^3*s ({2:.4g})'.format(
                species, species_rate, species_rate / char_rate))
            if network is not None:
                logging.info('    PDepNetwork #{0:d} leak rate: {1:10.4e} mol/m^3*s ({2:.4g})'.format(
                    network.index, network_rate, network_rate / char_rate))

    cpdef log_conversions(self, species_index, y0):
        """
        Log information about the current conversion values.
        """
        for term in self.termination:
            if isinstance(term, TerminationConversion):
                index = species_index[term.species]
                X = 1 - (self.y[index] / y0[index])
                logging.info('    {0} conversion: {1:<10.4g}'.format(term.species, X))

    @cython.boundscheck(False)
    def compute_rate_derivative(self):
        """
        Returns derivative vector df/dk_j where dy/dt = f(y, t, k) and
        k_j is the rate parameter for the jth core reaction.
        """
        cdef np.ndarray[np.int_t, ndim=2] ir, ip
        cdef np.ndarray[np.float64_t, ndim=1] kf, kr, C, deriv
        cdef np.ndarray[np.float64_t, ndim=2] rate_deriv
        cdef double fderiv, rderiv, flux, V
        cdef int j, num_core_reactions, num_core_species

        cdef double RT_inverse, gderiv

        ir = self.reactant_indices
        ip = self.product_indices

        kf = self.kf
        kr = self.kb

        num_core_reactions = len(self.core_reaction_rates)
        num_core_species = len(self.core_species_concentrations)

        # Use stored volume, since this function is only called from residual function. 
        RT_inverse = 1 / (constants.R * self.T.value_si)
        V = self.V

        C = self.core_species_concentrations

        rate_deriv = np.zeros((num_core_species, num_core_reactions + num_core_species), np.float64)

        for j in range(num_core_reactions):
            if ir[j, 1] == -1:  # only one reactant
                fderiv = C[ir[j, 0]]
            elif ir[j, 2] == -1:  # only two reactants
                fderiv = C[ir[j, 0]] * C[ir[j, 1]]
            else:  # three reactants
                fderiv = C[ir[j, 0]] * C[ir[j, 1]] * C[ir[j, 2]]

            if ip[j, 1] == -1:  # only one reactant
                rderiv = kr[j] / kf[j] * C[ip[j, 0]]
            elif ip[j, 2] == -1:  # only two reactants
                rderiv = kr[j] / kf[j] * C[ip[j, 0]] * C[ip[j, 1]]
            else:  # three reactants
                rderiv = kr[j] / kf[j] * C[ip[j, 0]] * C[ip[j, 1]] * C[ip[j, 2]]

            flux = fderiv - rderiv
            gderiv = rderiv * kf[j] * RT_inverse

            deriv = np.zeros(num_core_species, np.float64)  # derivative for reaction j with respect to dG_species i

            deriv[ir[j, 0]] += gderiv
            if ir[j, 1] != -1:  # only two reactants
                deriv[ir[j, 1]] += gderiv
                if ir[j, 2] != -1:  # three reactants
                    deriv[ir[j, 2]] += gderiv

            deriv[ip[j, 0]] -= gderiv
            if ip[j, 1] != -1:  # only two reactants
                deriv[ip[j, 1]] -= gderiv
                if ip[j, 2] != -1:  # three reactants
                    deriv[ip[j, 2]] -= gderiv

            rate_deriv[ir[j, 0], j] -= flux
            rate_deriv[ir[j, 0], num_core_reactions:num_core_reactions + num_core_species] -= deriv
            if ir[j, 1] != -1:
                rate_deriv[ir[j, 1], j] -= flux
                rate_deriv[ir[j, 1], num_core_reactions:num_core_reactions + num_core_species] -= deriv
                if ir[j, 2] != -1:
                    rate_deriv[ir[j, 2], j] -= flux
                    rate_deriv[ir[j, 2], num_core_reactions:num_core_reactions + num_core_species] -= deriv

            rate_deriv[ip[j, 0], j] += flux
            rate_deriv[ip[j, 0], num_core_reactions:num_core_reactions + num_core_species] += deriv
            if ip[j, 1] != -1:
                rate_deriv[ip[j, 1], j] += flux
                rate_deriv[ip[j, 1], num_core_reactions:num_core_reactions + num_core_species] += deriv
                if ip[j, 2] != -1:
                    rate_deriv[ip[j, 2], j] += flux
                    rate_deriv[ip[j, 2], num_core_reactions:num_core_reactions + num_core_species] += deriv

        rate_deriv = V * rate_deriv

        return rate_deriv
