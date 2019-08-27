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
from rmgpy.chemkin import getSpeciesIdentifier
from rmgpy.reaction import Reaction
from rmgpy.quantity import Quantity
from rmgpy.species import Species

################################################################################

cdef class ReactionSystem(DASx):
    """
    A base class for all RMG reaction systems.
    """

    def __init__(self, termination=None, sensitiveSpecies=None, sensitivityThreshold=1e-3):
        DASx.__init__(self)

        # reactor state variables:
        self.t0 = 0.0
        self.y0 = None
        self.dydt0 = None

        #  variables that determine the dimensions of arrays and matrices:
        self.numCoreSpecies = -1
        self.numCoreReactions = -1
        self.numEdgeSpecies = -1
        self.numEdgeReactions = -1
        self.numPdepNetworks = -1

        """
        The number of differential equations that will 
        be solved simulatenously by the solver.

        """
        self.neq = -1

        # variables that store stoichiometry data
        """
        speciesIndex is a dictionary with species as keys, and 
        values equal to the index integer that is generated
        when enumerating the core and edge species .
        """
        self.speciesIndex = {}

        """
        reactionIndex is a dictionary with reactions as keys, and 
        values equal to the index integer that is generated
        when enumerating the core and edge reactions .
        """
        self.reactionIndex = {}

        """
        A matrix for the reactants and products.

        The matrix has dimensions n x 3, with :
        - n the sum of the number of core and edge reactions,
        - 3 the maximum number of molecules allowed in either the reactant or
            product side of a reaction.

        The matrix element (j,l), with
        - j the index of the reaction and 
        - l the index of the lth reactant/product in the reaction contains the index
            of the molecule in the speciesIndex dictionary.
        """
        self.reactantIndices = None
        self.productIndices = None

        self.networkIndices = None

        # matrices that cache kinetic and rate data
        self.kf = None  # forward rate coefficients
        self.kb = None  # reverse rate coefficients
        self.Keq = None  # equilibrium constants
        self.networkLeakCoefficients = None
        self.jacobianMatrix = None

        self.coreSpeciesConcentrations = None

        # The reaction and species rates at the current time (in mol/m^3*s)
        self.coreSpeciesRates = None
        self.coreReactionRates = None
        self.coreSpeciesProductionRates = None
        self.coreSpeciesConsumptionRates = None
        self.edgeSpeciesRates = None
        self.edgeReactionRates = None

        self.networkLeakRates = None

        #surface indices
        self.surfaceSpeciesIndices = None
        self.surfaceReactionIndices = None
        self.validLayeringIndices = None

        # variables that cache maximum rate (ratio) data
        self.maxEdgeSpeciesRateRatios = None
        self.maxNetworkLeakRateRatios = None

        #for managing prunable edge species
        self.prunableSpecies = []
        self.prunableNetworks = []
        self.prunableSpeciesIndices = None
        self.prunableNetworkIndices = None

        # sensitivity variables
        self.sensmethod = 2  # sensmethod = 1 for staggered corrector sensitivities, 0 (simultaneous corrector), 2 (staggered direct)
        self.sensitivityCoefficients = None
        self.sensitiveSpecies = sensitiveSpecies
        self.sensitivityThreshold = sensitivityThreshold
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

        # reaction filtration, unimolecularThreshold is a vector with length of number of core species
        # bimolecularThreshold is a square matrix with length of number of core species
        # trimolecularThreshold is a rank-3 tensor with length of number of core species
        # A value of 1 in the matrix indicates the species is above the threshold to react or participate in those reactions
        self.unimolecularThreshold = None
        self.bimolecularThreshold = None
        self.trimolecularThreshold = None

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (self.__class__, (self.termination,))

    cpdef initializeModel(self, list coreSpecies, list coreReactions, list edgeSpecies, list edgeReactions,
                          list surfaceSpecies=None, list surfaceReactions=None, list pdepNetworks=None,
                          atol=1e-16, rtol=1e-8, sensitivity=False, sens_atol=1e-6, sens_rtol=1e-4,
                          filterReactions=False, dict conditions=None):
        """
        Initialize a simulation of the reaction system using the provided
        kinetic model. You will probably want to create your own version of this
        method in the derived class; don't forget to also call the base class
        version, too.
        """
        if surfaceSpecies is None:
            surfaceSpecies = []
        if surfaceReactions is None:
            surfaceReactions = []

        if conditions:
            is_conc = hasattr(self, 'initialConcentrations')
            is_surf = hasattr(self, 'initialGasMoleFractions')
            # ToDo: I think this block is incompatible with surface.pyx catalyst reactors
            keys = conditions.keys()
            if 'T' in keys and hasattr(self, 'T'):
                self.T = Quantity(conditions['T'], 'K')
            if 'P' in keys and hasattr(self, 'P'):
                self.P = Quantity(conditions['P'], 'Pa')
            for k in keys:
                if is_conc:
                    if k in self.initialConcentrations.keys():
                        self.initialConcentrations[k] = conditions[k]  #already in SI units
                elif is_surf:
                    if k in self.initialGasMoleFractions.keys():
                        self.initialGasMoleFractions[k] = conditions[k]  #already in SI units
                    if k in self.initialSurfaceCoverages.keys():
                        self.initialSurfaceCoverages[k] = conditions[k]  #already in SI units
                else:
                    if k in self.initialMoleFractions.keys():
                        self.initialMoleFractions[k] = conditions[k]

        self.numCoreSpecies = len(coreSpecies)
        self.numCoreReactions = len(coreReactions)
        self.numEdgeSpecies = len(edgeSpecies)
        self.numEdgeReactions = len(edgeReactions)

        self.initiate_tolerances(atol, rtol, sensitivity, sens_atol, sens_rtol)

        pdepNetworks = pdepNetworks or []
        self.numPdepNetworks = len(pdepNetworks)

        self.kf = np.zeros((self.numCoreReactions + self.numEdgeReactions), np.float64)
        self.kb = np.zeros_like(self.kf)
        self.Keq = np.zeros_like(self.kf)

        self.generate_species_indices(coreSpecies, edgeSpecies)
        self.generate_reaction_indices(coreReactions, edgeReactions)
        self.generate_reactant_product_indices(coreReactions, edgeReactions)

        self.coreSpeciesConcentrations = np.zeros((self.numCoreSpecies), np.float64)
        self.coreSpeciesProductionRates = np.zeros((self.numCoreSpecies), np.float64)
        self.coreSpeciesConsumptionRates = np.zeros((self.numCoreSpecies), np.float64)
        self.coreReactionRates = np.zeros((self.numCoreReactions), np.float64)
        self.edgeReactionRates = np.zeros((self.numEdgeReactions), np.float64)
        self.coreSpeciesRates = np.zeros((self.numCoreSpecies), np.float64)
        self.edgeSpeciesRates = np.zeros((self.numEdgeSpecies), np.float64)
        self.networkLeakRates = np.zeros((self.numPdepNetworks), np.float64)
        self.maxNetworkLeakRateRatios = np.zeros((len(self.prunableNetworks)), np.float64)
        self.sensitivityCoefficients = np.zeros((self.numCoreSpecies, self.numCoreReactions), np.float64)
        self.unimolecularThreshold = np.zeros((self.numCoreSpecies), bool)
        self.bimolecularThreshold = np.zeros((self.numCoreSpecies, self.numCoreSpecies), bool)
        if self.trimolecular:
            self.trimolecularThreshold = np.zeros((self.numCoreSpecies, self.numCoreSpecies, self.numCoreSpecies),
                                                     bool)

        surfaceSpecies, surfaceReactions = self.initialize_surface(coreSpecies, coreReactions, surfaceSpecies,
                                                                   surfaceReactions)

        self.set_prunable_indices(edgeSpecies, pdepNetworks)

    def initialize_solver(self):
        DASx.initialize(self, self.t0, self.y0, self.dydt0, self.senpar, self.atol_array, self.rtol_array)

    def reset_max_edge_species_rate_ratios(self):
        """
        This function sets maxEdgeSpeciesRateRatios back to zero
        for pruning of ranged reactors it is important to avoid doing this
        every initialization
        """
        self.maxEdgeSpeciesRateRatios = np.zeros((len(self.prunableSpecies)), np.float64)

    def set_prunable_indices(self, edgeSpecies, pdepNetworks):
        cdef object spc
        cdef list temp
        temp = []
        for i, spc in enumerate(self.prunableSpecies):
            try:
                temp.append(edgeSpecies.index(spc))
            except ValueError:
                self.maxEdgeSpeciesRateRatios[i] = np.inf  #avoid pruning of species that have been moved to core

        self.prunableSpeciesIndices = np.array(temp)

        temp = []
        for i, spc in enumerate(self.prunableNetworks):
            try:
                temp.append(pdepNetworks.index(spc))
            except:
                self.maxNetworkLeakRateRatios[i] = np.inf  #avoid pruning of lost networks

        self.prunableNetworkIndices = np.array(temp)

    @cython.boundscheck(False)
    cpdef initialize_surface(self, list coreSpecies, list coreReactions, list surfaceSpecies, list surfaceReactions):
        """
        removes surfaceSpecies and surfaceReactions from  until they are self consistent: 
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

        product_indices = self.productIndices
        reactant_indices = self.reactantIndices

        surface_species_indices = []
        surface_reaction_indices = []
        possible_species_indices = set()
        remove_inds = []

        for obj in surfaceSpecies:
            surface_species_indices.append(coreSpecies.index(obj))

        for obj in surfaceReactions:
            surface_reaction_indices.append(coreReactions.index(obj))

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
                logging.info('removing disconnected reaction from surface: {0}'.format(str(coreReactions[i])))
                surfaceReactions.remove(coreReactions[i])
                surface_reaction_indices.remove(i)

        possible_species_indices -= {-1}  #remove the -1 indexes that indicate there is no third/second reactant/product

        for i in surface_species_indices:  #remove species without reactions in the surface
            if not (i in possible_species_indices):
                logging.info('removing disconnected species from surface: {0}'.format(coreSpecies[i].label))
                remove_inds.append(i)

        for i in remove_inds:
            surfaceSpecies.remove(coreSpecies[i])
            surface_species_indices.remove(i)

        self.surfaceSpeciesIndices = np.array(surface_species_indices, dtype=np.int)
        self.surfaceReactionIndices = np.array(surface_reaction_indices, dtype=np.int)

        self.validLayeringIndices = self.getLayeringIndices()

        surfaceSpecies = [coreSpecies[i] for i in surface_species_indices]
        surfaceReactions = [coreReactions[i] for i in surface_reaction_indices]

        logging.debug('Surface initialization complete')

        return surfaceSpecies, surfaceReactions

    def initiate_tolerances(self, atol=1e-16, rtol=1e-8, sensitivity=False, sens_atol=1e-6, sens_rtol=1e-4):
        """
        Computes the number of differential equations and initializes the tolerance arrays.        
        """
        # Compute number of equations    
        if sensitivity:
            # Set DASPK sensitivity analysis to ON
            self.sensitivity = True
            # Compute number of variables
            self.neq = self.numCoreSpecies * (self.numCoreReactions + self.numCoreSpecies + 1)

            self.atol_array = np.ones(self.neq, np.float64) * sens_atol
            self.atol_array[:self.numCoreSpecies] = atol

            self.rtol_array = np.ones(self.neq, np.float64) * sens_rtol
            self.rtol_array[:self.numCoreSpecies] = rtol

            self.senpar = np.zeros(self.numCoreReactions + self.numCoreSpecies, np.float64)

        else:
            self.neq = self.numCoreSpecies

            self.atol_array = np.ones(self.neq, np.float64) * atol
            self.rtol_array = np.ones(self.neq, np.float64) * rtol

            self.senpar = np.zeros(self.numCoreReactions, np.float64)

    def get_species_index(self, spc):
        """
        Retrieves the index that is associated with the parameter species
        from the species index dictionary.
        """
        return self.speciesIndex[spc]

    def generate_reactant_product_indices(self, coreReactions, edgeReactions):
        """
        Creates a matrix for the reactants and products.
        """

        self.reactantIndices = -np.ones((self.numCoreReactions + self.numEdgeReactions, 3), np.int)
        self.productIndices = -np.ones_like(self.reactantIndices)

        for rxn in itertools.chain(coreReactions, edgeReactions):
            j = self.reactionIndex[rxn]
            for l, spec in enumerate(rxn.reactants):
                i = self.get_species_index(spec)
                self.reactantIndices[j, l] = i
            for l, spec in enumerate(rxn.products):
                i = self.get_species_index(spec)
                self.productIndices[j, l] = i

    def generate_species_indices(self, coreSpecies, edgeSpecies):
        """
        Assign an index to each species (core first, then edge) and 
        store the (species, index) pair in a dictionary.
        """

        for index, spec in enumerate(itertools.chain(coreSpecies, edgeSpecies)):
            self.speciesIndex[spec] = index

    def generate_reaction_indices(self, coreReactions, edgeReactions):
        """
        Assign an index to each reaction (core first, then edge) and 
        store the (reaction, index) pair in a dictionary.
        """

        for index, rxn in enumerate(itertools.chain(coreReactions, edgeReactions)):
            self.reactionIndex[rxn] = index

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
        num_core_species = len(self.coreSpeciesConcentrations)
        for i in range(num_core_species):
            if self.coreSpeciesConcentrations[i] > 0:
                self.unimolecularThreshold[i] = True
        for i in range(num_core_species):
            for j in range(i, num_core_species):
                if self.coreSpeciesConcentrations[i] > 0 and self.coreSpeciesConcentrations[j] > 0:
                    self.bimolecularThreshold[i, j] = True
        if self.trimolecular:
            for i in range(num_core_species):
                for j in range(i, num_core_species):
                    for k in range(j, num_core_species):
                        if (self.coreSpeciesConcentrations[i] > 0
                                and self.coreSpeciesConcentrations[j] > 0
                                and self.coreSpeciesConcentrations[k] > 0):
                            self.trimolecularThreshold[i, j, k] = True

    def set_initial_derivative(self):
        """
        Sets the derivative of the species moles with respect to the independent variable (time)
        equal to the residual.
        """
        self.dydt0 = - self.residual(self.t0, self.y0, np.zeros(self.neq, np.float64), self.senpar)[0]

    def compute_network_variables(self, pdepNetworks=None):
        """
        Initialize the arrays containing network information:

        - NetworkLeakCoefficients is a n x 1 array with 
            n the number of pressure-dependent networks.
        - NetworkIndices is a n x 3 matrix with 
            n the number of pressure-dependent networks and 
            3 the maximum number of molecules allowed in either the reactant or
            product side of a reaction. 
        """

        pdepNetworks = pdepNetworks or []

        self.networkIndices = -np.ones((self.numPdepNetworks, 3), np.int)
        self.networkLeakCoefficients = np.zeros((self.numPdepNetworks), np.float64)

        for j, network in enumerate(pdepNetworks):
            self.networkLeakCoefficients[j] = network.getLeakCoefficient(self.T.value_si, self.P.value_si)
            for l, spec in enumerate(network.source):
                try:
                    i = self.get_species_index(spec)
                except KeyError:
                    # spec is probably in the edge, hence is not a key in the speciesIndex dictionary.
                    # Sicne networkIndices is only used to identify the number of reactants, set the
                    # corresponding value to be different than `-1`.
                    i = -2
                self.networkIndices[j, l] = i

    @cython.boundscheck(False)
    cpdef getLayeringIndices(self):
        """
        determines the edge reaction indices that indicate reactions that are valid for movement from
        edge to surface based on the layering constraint
        """

        cdef np.ndarray[np.int_t, ndim=1] surf_species_indices
        cdef int i, j, k, index, num_core_species, num_core_reactions, num_edge_reactions
        cdef list valid_indices

        surf_species_indices = self.surfaceSpeciesIndices
        num_core_species = self.numCoreSpecies
        num_core_reactions = self.numCoreReactions
        num_edge_reactions = self.numEdgeReactions
        product_indices = self.productIndices
        reactant_indices = self.reactantIndices

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

    cpdef addReactionsToSurface(self, list newSurfaceReactions, list newSurfaceReactionInds, list surfaceSpecies,
                                list surfaceReactions, list edgeSpecies):
        """
        moves new surface reactions to the surface
        done after the while loop before the simulate call ends
        """
        cdef object srxn
        cdef int k, sind, num_core_species, i, num_core_reactions
        cdef np.ndarray[np.int_t, ndim=2] product_indices, reactant_indices

        product_indices = self.productIndices
        reactant_indices = self.reactantIndices
        num_core_species = self.numCoreSpecies
        num_core_reactions = self.numCoreReactions

        for k in range(len(newSurfaceReactions)):
            srxn = newSurfaceReactions[k]
            sind = newSurfaceReactionInds[k]
            surfaceReactions.append(srxn)  #add to surface trackers

            for i in product_indices[sind + num_core_reactions]:
                if i >= num_core_species:
                    surfaceSpecies.append(edgeSpecies[i - num_core_species])
            for i in reactant_indices[sind + num_core_reactions]:
                if i >= num_core_species:
                    surfaceSpecies.append(edgeSpecies[i - num_core_species])

        return surfaceSpecies, surfaceReactions

    @cython.boundscheck(False)
    cpdef simulate(self, list coreSpecies, list coreReactions, list edgeSpecies,
                   list edgeReactions, list surfaceSpecies, list surfaceReactions,
                   list pdepNetworks=None, bool prune=False, bool sensitivity=False, list sensWorksheet=None,
                   object modelSettings=None, object simulatorSettings=None, dict conditions=None):
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
        cdef bool useDynamicsTemp, first_time, use_dynamics, terminate_at_max_objects, schanged
        cdef np.ndarray[np.float64_t, ndim=1] edge_reaction_rates
        cdef double reaction_rate, production, consumption
        cdef np.ndarray[np.int_t, ndim=1] surface_species_indices, surface_reaction_indices
        # cython declations for sensitivity analysis
        cdef np.ndarray[np.int_t, ndim=1] sens_species_indices, reactant_side, product_side
        cdef np.ndarray[np.float64_t, ndim=1] mole_sens, dVdk, norm_sens
        cdef list time_array, norm_sens_array, new_surface_reactions, new_surface_reaction_inds, new_objects, new_object_inds

        zero_production = False
        zero_consumption = False
        pdepNetworks = pdepNetworks or []

        schanged = False

        num_core_species = len(coreSpecies)
        num_edge_species = len(edgeSpecies)
        num_pdep_networks = len(pdepNetworks)
        num_core_reactions = len(coreReactions)

        assert set(coreReactions) >= set(surfaceReactions), 'given surface reactions are not a subset of core reactions'
        assert set(coreSpecies) >= set(surfaceSpecies), 'given surface species are not a subset of core species'

        tol_keep_in_edge = modelSettings.fluxToleranceKeepInEdge if prune else 0
        tol_move_to_core = modelSettings.fluxToleranceMoveToCore
        tol_move_edge_reaction_to_core = modelSettings.toleranceMoveEdgeReactionToCore
        tol_interrupt_simulation = modelSettings.fluxToleranceInterrupt
        tol_move_edge_reaction_to_core_interrupt = modelSettings.toleranceMoveEdgeReactionToCore
        tol_move_edge_reaction_to_surface = modelSettings.toleranceMoveEdgeReactionToSurface
        tol_move_surface_species_to_core = modelSettings.toleranceMoveSurfaceSpeciesToCore
        tol_move_surface_reaction_to_core = modelSettings.toleranceMoveSurfaceReactionToCore
        tol_move_edge_reaction_to_surface_interrupt = modelSettings.toleranceMoveEdgeReactionToSurfaceInterrupt
        ignore_overall_flux_criterion = modelSettings.ignoreOverallFluxCriterion
        atol = simulatorSettings.atol
        rtol = simulatorSettings.rtol
        sens_atol = simulatorSettings.sens_atol
        sens_rtol = simulatorSettings.sens_rtol
        filter_reactions = modelSettings.filterReactions
        max_num_objs_per_iter = modelSettings.maxNumObjsPerIter

        if modelSettings.toleranceBranchReactionToCore != 0.0:
            branch_factor = 1.0 / modelSettings.toleranceBranchReactionToCore
            br_max = modelSettings.branchingRatioMax
            branching_index = modelSettings.branchingIndex
        else:
            branch_factor = 0.0

        #if not pruning always terminate at max objects, otherwise only do so if terminate_at_max_objects=True
        terminate_at_max_objects = True if not prune else modelSettings.terminateAtMaxObjects

        dynamics_time_scale = modelSettings.dynamicsTimeScale

        use_dynamics = not (tol_move_edge_reaction_to_core == np.inf and tol_move_edge_reaction_to_surface == np.inf)

        species_index = {}
        for index, spec in enumerate(coreSpecies):
            species_index[spec] = index

        self.initializeModel(coreSpecies, coreReactions, edgeSpecies, edgeReactions, surfaceSpecies, surfaceReactions,
                             pdepNetworks, atol, rtol, sensitivity,
                             sens_atol, sens_rtol,
                             filter_reactions, conditions)

        prunable_species_indices = self.prunableSpeciesIndices
        prunable_network_indices = self.prunableNetworkIndices

        surface_species_indices = self.surfaceSpeciesIndices
        surface_reaction_indices = self.surfaceReactionIndices

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

        max_edge_species_rate_ratios = self.maxEdgeSpeciesRateRatios
        max_network_leak_rate_ratios = self.maxNetworkLeakRateRatios
        forward_rate_coefficients = self.kf
        unimolecular_threshold = self.unimolecularThreshold
        bimolecular_threshold = self.bimolecularThreshold
        trimolecular_threshold = self.trimolecularThreshold

        # Copy the initial conditions to use in evaluating conversions
        y0 = self.y.copy()

        # a list with the time, Volume, number of moles of core species
        self.snapshots = []

        if sensitivity:
            time_array = []
            norm_sens_array = [[] for spec in self.sensitiveSpecies]
            RTP = constants.R * self.T.value_si / self.P.value_si
            # identify sensitive species indices
            sens_species_indices = np.array([species_index[spec] for spec in self.sensitiveSpecies],
                                               np.int)  # index within coreSpecies list of the sensitive species

        step_time = 1e-12
        prev_time = self.t

        first_time = True

        while not terminated:
            # Integrate forward in time by one time step

            if not first_time:
                try:
                    self.step(step_time)
                    if np.isnan(self.y).any():
                        raise DASxError("nans in moles")
                except DASxError as e:
                    logging.error("Trying to step from time {} to {} resulted in a solver (DASPK) error: "
                                  "{}".format(prev_time, step_time, e.message))

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
                            obj = edgeSpecies[ind]
                            logging.info('At time {0:10.4e} s, species {1} at rate ratio {2} was added to model core '
                                         'in model resurrection process'.format(self.t, obj,edge_species_rates[ind]))
                            invalid_objects.append(obj)

                        if total_div_accum_nums and len(total_div_accum_nums) > 0:  #if dynamics data available
                            ind = np.argmax(total_div_accum_nums)
                            obj = edgeReactions[ind]
                            logging.info('At time {0:10.4e} s, Reaction {1} at dynamics number {2} was added to model core '
                                         'in model resurrection process'.format(self.t, obj,total_div_accum_nums[ind]))
                            invalid_objects.append(obj)

                        if pdepNetworks != [] and network_leak_rate_ratios != []:
                            ind = np.argmax(network_leak_rate_ratios)
                            obj = pdepNetworks[ind]
                            logging.info('At time {0:10.4e} s, PDepNetwork #{1:d} at network leak rate {2} '
                                         'was sent for exploring during model resurrection process'
                                         ''.format(self.t, obj.index, network_leak_rate_ratios[ind]))
                            invalid_objects.append(obj)

                    if invalid_objects != []:
                        return False, True, invalid_objects, surfaceSpecies, surfaceReactions, self.t, conversion
                    else:
                        logging.error('Model Resurrection has failed')
                        logging.error("Core species names: {!r}".format([getSpeciesIdentifier(s) for s in coreSpecies]))
                        logging.error("Core species moles: {!r}".format(self.y[:num_core_species]))
                        logging.error("Volume: {!r}".format(self.V))
                        logging.error("Core species net rates: {!r}".format(self.coreSpeciesRates))
                        logging.error("Edge species net rates: {!r}".format(self.edgeSpeciesRates))
                        logging.error("Network leak rates: {!r}".format(self.networkLeakRates))
                        raise ValueError('invalid_objects could not be filled during resurrection process')

            y_core_species = self.y[:num_core_species]
            total_moles = np.sum(y_core_species)
            if sensitivity:
                time_array.append(self.t)
                mole_sens = self.y[num_core_species:]
                volume = self.V

                dVdk = np.zeros(num_core_reactions + num_core_species, np.float64)
                if not self.constantVolume:
                    for j in range(num_core_reactions + num_core_species):
                        dVdk[j] = np.sum(mole_sens[j * num_core_species:(j + 1) * num_core_species]) * RTP  # Contains [ dV_dk and dV_dG ]
                for i in range(len(self.sensitiveSpecies)):
                    norm_sens = np.zeros(num_core_reactions + num_core_species, np.float64)
                    c = self.coreSpeciesConcentrations[sens_species_indices[i]]
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
            char_rate = sqrt(np.sum(self.coreSpeciesRates * self.coreSpeciesRates))

            if char_rate > max_char_rate:
                max_char_rate = char_rate

            core_species_rates = np.abs(self.coreSpeciesRates)
            edge_reaction_rates = self.edgeReactionRates
            core_species_consumption_rates = self.coreSpeciesConsumptionRates
            core_species_production_rates = self.coreSpeciesProductionRates
            edge_species_rates = np.abs(self.edgeSpeciesRates)
            network_leak_rates = np.abs(self.networkLeakRates)
            core_species_rate_ratios = np.abs(self.coreSpeciesRates / char_rate)
            edge_species_rate_ratios = np.abs(self.edgeSpeciesRates / char_rate)
            network_leak_rate_ratios = np.abs(self.networkLeakRates / char_rate)
            num_edge_reactions = self.numEdgeReactions
            core_reaction_rates = self.coreReactionRates
            product_indices = self.productIndices
            reactant_indices = self.reactantIndices
            core_species_concentrations = self.coreSpeciesConcentrations

            # Update the maximum species rate and maximum network leak rate arrays
            for i, index in enumerate(prunable_species_indices):
                if max_edge_species_rate_ratios[i] < edge_species_rate_ratios[index]:
                    max_edge_species_rate_ratios[i] = edge_species_rate_ratios[index]
            for i, index in enumerate(prunable_network_indices):
                if max_network_leak_rate_ratios[i] < network_leak_rate_ratios[index]:
                    max_network_leak_rate_ratios[i] = network_leak_rate_ratios[index]

            if char_rate == 0 and len(edge_species_rates) > 0:  # this deals with the case when there is no flux in the system
                max_species_index = np.argmax(edge_species_rates)
                max_species = edgeSpecies[max_species_index]
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
                        reactant_side = self.reactantIndices[index + num_core_reactions, :]
                        product_side = self.productIndices[index + num_core_reactions, :]
                    else:
                        reactant_side = self.productIndices[index + num_core_reactions, :]
                        product_side = self.reactantIndices[index + num_core_reactions, :]

                    mults = []
                    for i in product_side:
                        if i == -1:
                            continue
                        elif i < num_core_species:
                            mults.append(coreSpecies[i].molecule[0].multiplicity)
                        else:
                            mults.append(edgeSpecies[i - num_core_species].molecule[0].multiplicity)

                    if max(mults) > 2:
                        continue

                    for spc_index in reactant_side:
                        if spc_index != -1 and spc_index < num_core_species:
                            if coreSpecies[spc_index].molecule[0].multiplicity != 2:
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
                    for spc_index in self.reactantIndices[index + num_core_reactions, :]:
                        if spc_index != -1 and spc_index < num_core_species:
                            consumption = core_species_consumption_rates[spc_index]
                            if consumption != 0:  #if consumption = 0 ignore species
                                total_div_accum_nums[index] *= (reaction_rate + consumption) / consumption

                    for spc_index in self.productIndices[index + num_core_reactions, :]:
                        if spc_index != -1 and spc_index < num_core_species:
                            production = core_species_production_rates[spc_index]
                            if production != 0:  #if production = 0 ignore species
                                total_div_accum_nums[index] *= (reaction_rate + production) / production

                total_div_ln_accum_nums = np.log(total_div_accum_nums)

                surface_species_indices = self.surfaceSpeciesIndices
                surface_reaction_indices = self.surfaceReactionIndices

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
                        surface_objects.append(coreReactions[sind])

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
                        surface_objects.append(coreSpecies[sind])

                #process objects to be moved from surface to core

                if len(surface_objects) > 0:
                    schanged = True
                    for ind, obj in enumerate(surface_objects):
                        if isinstance(obj, Reaction):
                            logging.info('Moving reaction {0} from surface to core'.format(obj))
                            surfaceReactions.remove(obj)
                        elif isinstance(obj, Species):
                            logging.info('Moving species {0} from surface to core'.format(obj))
                            surfaceSpecies.remove(obj)
                        else:
                            raise ValueError
                    surfaceSpecies, surfaceReactions = self.initialize_surface(coreSpecies, coreReactions,
                                                                               surfaceSpecies, surfaceReactions)
                    logging.info('Surface now has {0} Species and {1} Reactions'
                                 ''.format(len(self.surfaceSpeciesIndices), len(self.surfaceReactionIndices)))

            if filter_reactions:
                # Calculate thresholds for reactions
                (unimolecular_threshold_rate_constant,
                 bimolecular_threshold_rate_constant,
                 trimolecular_threshold_rate_constant) = self.get_threshold_rate_constants(modelSettings)
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
                for ind, obj in enumerate(edgeSpecies):
                    rr = edge_species_rate_ratios[ind]
                    if rr > tol_move_to_core:
                        if not (obj in new_objects or obj in invalid_objects):
                            temp_new_objects.append(edgeSpecies[ind])
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
                for ind, obj in enumerate(edgeReactions):
                    b_num = branching_nums[ind]
                    if b_num > 1:
                        if not (obj in new_objects or obj in invalid_objects):
                            temp_new_objects.append(edgeReactions[ind])
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
                valid_layering_indices = self.validLayeringIndices
                temp_surface_objects = []

                for ind, obj in enumerate(edgeReactions):
                    dlnaccum = total_div_ln_accum_nums[ind]
                    if dlnaccum > tol_move_edge_reaction_to_core:
                        if not (obj in new_objects or obj in invalid_objects):
                            temp_new_objects.append(edgeReactions[ind])
                            temp_new_object_inds.append(ind)
                            temp_new_object_vals.append(dlnaccum)
                            temp_new_object_type.append('dyn')
                    elif dlnaccum > tol_move_edge_reaction_to_surface and ind in valid_layering_indices:
                        if not (obj in new_objects or obj in invalid_objects):
                            temp_new_objects.append(edgeReactions[ind])
                            temp_new_object_inds.append(ind)
                            temp_new_object_vals.append(dlnaccum)
                            temp_surface_objects.append(edgeReactions[ind])
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

            #Determination of pdepNetworks in need of exploring

            if pdepNetworks:
                for ind, obj in enumerate(pdepNetworks):
                    lr = network_leak_rate_ratios[ind]
                    if lr > tol_move_to_core:
                        if not (obj in new_objects or obj in invalid_objects):
                            temp_new_objects.append(pdepNetworks[ind])
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
                logging.info('Exceeded max number of objects...removing excess objects')
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
                surfaceSpecies, surfaceReactions = self.initialize_surface(coreSpecies, coreReactions, surfaceSpecies,
                                                                           surfaceReactions)
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
                        self.logConversions(species_index, y0)
                        break
                elif isinstance(term, TerminationConversion):
                    index = species_index[term.species]
                    conversion = 1 - (y_core_species[index] / y0[index])
                    if 1 - (y_core_species[index] / y0[index]) > term.conversion:
                        terminated = True
                        logging.info('At time {0:10.4e} s, reached target termination conversion: {1:f} of '
                                     '{2}'.format(self.t,term.conversion,term.species))
                        self.logConversions(species_index, y0)
                        break
                elif isinstance(term, TerminationRateRatio):
                    if max_char_rate != 0.0 and char_rate / max_char_rate < term.ratio:
                        terminated = True
                        logging.info('At time {0:10.4e} s, reached target termination RateRatio: '
                                     '{1}'.format(self.t,char_rate/max_char_rate))
                        self.logConversions(species_index, y0)

            # Increment destination step time if necessary
            if self.t >= 0.9999 * step_time:
                step_time *= 10.0

        # Change surface species and reactions based on what will be added to the surface
        surfaceSpecies, surfaceReactions = self.addReactionsToSurface(new_surface_reactions, new_surface_reaction_inds,
                                                                      surfaceSpecies, surfaceReactions, edgeSpecies)

        # notify reaction system listeners
        self.notify()

        if sensitivity:
            for i in range(len(self.sensitiveSpecies)):
                with open(sensWorksheet[i], 'w') as outfile:
                    worksheet = csv.writer(outfile)
                    reactions_above_threshold = []
                    for j in range(num_core_reactions + num_core_species):
                        for k in range(len(time_array)):
                            if abs(norm_sens_array[i][k][j]) > self.sensitivityThreshold:
                                reactions_above_threshold.append(j)
                                break
                    species_name = getSpeciesIdentifier(self.sensitiveSpecies[i])
                    headers = ['Time (s)']
                    headers.extend(['dln[{0}]/dln[k{1}]: {2}'.format(species_name, j+1, coreReactions[j].toChemkin(kinetics=False)) if j < num_core_reactions
                                    else 'dln[{0}]/dG[{1}]'.format(species_name, getSpeciesIdentifier(coreSpecies[j-num_core_reactions])) for j in reactions_above_threshold])
                    worksheet.writerow(headers)

                    for k in range(len(time_array)):
                        row = [time_array[k]]
                        row.extend([norm_sens_array[i][k][j] for j in reactions_above_threshold])
                        worksheet.writerow(row)

        self.maxEdgeSpeciesRateRatios = max_edge_species_rate_ratios
        self.maxNetworkLeakRateRatios = max_network_leak_rate_ratios

        self.unimolecularThreshold = unimolecular_threshold
        self.bimolecularThreshold = bimolecular_threshold
        self.trimolecularThreshold = trimolecular_threshold

        # Return the invalid object (if the simulation was invalid) or None
        # (if the simulation was valid)
        return terminated, False, invalid_objects, surfaceSpecies, surfaceReactions, self.t, conversion

    cpdef logRates(self, double charRate, object species, double speciesRate, double maxDifLnAccumNum, object network,
                   double networkRate):
        """
        Log information about the current maximum species and network rates.
        """
        logging.info('    Characteristic rate: {0:10.4e} mol/m^3*s'.format(charRate))
        logging.info('    Dynamics number:  {0:10.4e}'.format(maxDifLnAccumNum))
        if charRate == 0.0:
            logging.info('    {0} rate: {1:10.4e} mol/m^3*s'.format(species, speciesRate))
            if network is not None:
                logging.info('    PDepNetwork #{0:d} leak rate: {1:10.4e} mol/m^3*s'.format(network.index, networkRate))
        else:
            logging.info('    {0} rate: {1:10.4e} mol/m^3*s ({2:.4g})'.format(
                species, speciesRate, speciesRate / charRate))
            if network is not None:
                logging.info('    PDepNetwork #{0:d} leak rate: {1:10.4e} mol/m^3*s ({2:.4g})'.format(
                    network.index, networkRate, networkRate / charRate))

    cpdef logConversions(self, speciesIndex, y0):
        """
        Log information about the current conversion values.
        """
        for term in self.termination:
            if isinstance(term, TerminationConversion):
                index = speciesIndex[term.species]
                X = 1 - (self.y[index] / y0[index])
                logging.info('    {0} conversion: {1:<10.4g}'.format(term.species, X))

    @cython.boundscheck(False)
    def computeRateDerivative(self):
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

        ir = self.reactantIndices
        ip = self.productIndices

        kf = self.kf
        kr = self.kb

        num_core_reactions = len(self.coreReactionRates)
        num_core_species = len(self.coreSpeciesConcentrations)

        # Use stored volume, since this function is only called from residual function. 
        RT_inverse = 1 / (constants.R * self.T.value_si)
        V = self.V

        C = self.coreSpeciesConcentrations

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


################################################################################

class TerminationTime:
    """
    Represent a time at which the simulation should be terminated. This class
    has one attribute: the termination `time` in seconds.
    """

    def __init__(self, time=(0.0, 's')):
        self.time = Quantity(time)


################################################################################

class TerminationConversion:
    """
    Represent a conversion at which the simulation should be terminated. This
    class has two attributes: the `species` to monitor and the fractional
    `conversion` at which to terminate.
    """

    def __init__(self, spec=None, conv=0.0):
        self.species = spec
        self.conversion = conv


class TerminationRateRatio:
    """
    Represent a fraction of the maximum characteristic rate of the simulation
    at which the simulation should be terminated.  This class has one attribute
    the ratio between the current and maximum characteristic rates at which
    to terminate
    """

    def __init__(self, ratio=0.01):
        self.ratio = ratio
