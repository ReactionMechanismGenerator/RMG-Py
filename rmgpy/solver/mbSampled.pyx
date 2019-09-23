################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2017 Prof. William H. Green (whgreen@mit.edu), 
#   Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

"""
Contains the :class:`MBSampledReactor` class, providing a reaction system
consisting of a homogeneous, isothermal, isobaric batch reactor that is being 
sampled by a molecular beam (MB).
"""

import itertools
import logging

cimport cython
import numpy as np
cimport numpy as np

import rmgpy.constants as constants
cimport rmgpy.constants as constants
from rmgpy.quantity import Quantity, RateCoefficient
from rmgpy.quantity cimport ScalarQuantity
from rmgpy.solver.base cimport ReactionSystem


cdef class MBSampledReactor(ReactionSystem):
    """
    A reaction system consisting of a homogeneous, isothermal, isobaric batch
    reactor that is being sample by a molecular beam. 
    The sampling process is modeled as a unimolecular reaction.
    These assumptions allow for a number of optimizations that enable
    this solver to complete very rapidly, even for large kinetic models.

    This is currently only intended for use with the ``simulate.py``script,
    and cannot be used for a standard RMG job.
    """

    cdef public ScalarQuantity T
    cdef public ScalarQuantity P
    cdef public ScalarQuantity k_sampling
    cdef public double V
    cdef public bint constantVolume
    cdef public dict initialMoleFractions
    cdef public list constantSpeciesList

    # collider variables

    """
    pdepColliderKinetics:
    an array that contains a reference to the kinetics object of the reaction
    that has pressure dependent kinetics.
    """
    cdef public list pdepColliderKinetics

    """
    colliderEfficiencies:
    an array consisting of array elements, each element corresponding to a reaction.
    Each element is an array with each position in the array corresponding to the collider efficiency
    of the core species. The collider efficiency is set to 1 if the species was not found in the list
    of colliders.
    """
    cdef public np.ndarray colliderEfficiencies

    """
    pdepColliderReactionIndices: 
    array that contains the indices of those reactions that 
    have pressure dependent kinetics. E.g. [4, 10, 2, 123]
    """
    cdef public np.ndarray pdepColliderReactionIndices

    """
    pdepSpecificColliderKinetics:
    an array that contains a reference to the kinetics object of the reaction
    that has pressure dependent kinetics with a specific species as a third body collider.
    """
    cdef public list pdepSpecificColliderKinetics

    """
    specificColliderSpecies:
    a list that contains object references to species which are specific third body colliders
    in the respective reactions in pdepSpecificColliderReactionIndices.
    """
    cdef public list specificColliderSpecies

    """
    pdepSpecificColliderReactionIndices:
    an array that contains the indices of reactions that have
    a specifcCollider attribyte. E.g. [16, 155, 90]
    """
    cdef public np.ndarray pdepSpecificColliderReactionIndices

    def __init__(self, T, P, initialMoleFractions, k_sampling, constantSpeciesList, termination, sensitiveSpecies=None,
                 sensitivityThreshold=1e-3):
        ReactionSystem.__init__(self, termination, sensitiveSpecies, sensitivityThreshold)
        self.T = Quantity(T)
        self.P = Quantity(P)
        self.initialMoleFractions = initialMoleFractions
        self.k_sampling = RateCoefficient(k_sampling)
        self.constantSpeciesList = constantSpeciesList

        self.V = 0  # will be set in initializeModel
        self.constantVolume = False

        self.pdepColliderReactionIndices = None
        self.pdepColliderKinetics = None
        self.colliderEfficiencies = None
        self.pdepSpecificColliderReactionIndices = None
        self.pdepSpecificColliderKinetics = None
        self.specificColliderSpecies = None

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (self.__class__,
                (self.T, self.P, self.initialMoleFractions, self.termination))

    def convertInitialKeysToSpeciesObjects(self, speciesDict):
        """
        Convert the initial_mole_fractions dictionary from species names into species objects,
        using the given dictionary of species.
        """
        initial_mole_fractions = {}
        for label, moleFrac in self.initialMoleFractions.items():
            initial_mole_fractions[speciesDict[label]] = moleFrac
        self.initialMoleFractions = initial_mole_fractions

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

        # First call the base class version of the method
        # This initializes the attributes declared in the base class
        ReactionSystem.initializeModel(self, coreSpecies=coreSpecies, coreReactions=coreReactions,
                                       edgeSpecies=edgeSpecies, edgeReactions=edgeReactions,
                                       surfaceSpecies=surfaceSpecies, surfaceReactions=surfaceReactions,
                                       pdepNetworks=pdepNetworks, atol=atol, rtol=rtol, sensitivity=sensitivity,
                                       sens_atol=sens_atol, sens_rtol=sens_rtol)

        # Set initial conditions
        self.set_initial_conditions()

        # Compute reaction thresholds if reaction filtering is turned on
        if filterReactions:
            ReactionSystem.set_initial_reaction_thresholds(self)

        self.set_colliders(coreReactions, edgeReactions, coreSpecies)

        ReactionSystem.compute_network_variables(self, pdepNetworks)

        # Generate forward and reverse rate coefficients k(T,P)
        self.generate_rate_coefficients(coreReactions, edgeReactions)

        ReactionSystem.set_initial_derivative(self)
        # Initialize the model
        ReactionSystem.initialize_solver(self)

    def calculate_effective_pressure(self, rxn):
        """
        Computes the effective pressure for a reaction as:
            Peff = P * sum(yi * effi / sum(y))
        with:
            - P the pressure of the reactor,
            - y the array of initial moles of the core species
        or as:
            Peff = P * y_specificCollider / sum(y)
        if a secificCollider is mentioned
        """

        y0_core_species = self.y0[:self.numCoreSpecies]
        sum_core_species = np.sum(y0_core_species)

        j = self.reactionIndex[rxn]
        for i in range(self.pdepColliderReactionIndices.shape[0]):
            if j == self.pdepColliderReactionIndices[i]:
                # Calculate effective pressure
                if rxn.specific_collider is None:
                    Peff = self.P.value_si * np.sum(self.colliderEfficiencies[i] * y0_core_species / sum_core_species)
                else:
                    logging.debug("Calculating Peff using {0} as a specific_collider".format(rxn.specific_collider))
                    Peff = self.P.value_si * self.y0[self.speciesIndex[rxn.specific_collider]] / sum_core_species
                return Peff
        return self.P.value_si

    def generate_rate_coefficients(self, coreReactions, edgeReactions):
        """
        Populates the forward rate coefficients (kf), reverse rate coefficients (kb)
        and equilibrium constants (Keq) arrays with the values computed at the temperature
        and (effective) pressure of the reaction system.
        """

        for rxn in itertools.chain(coreReactions, edgeReactions):
            j = self.reactionIndex[rxn]
            Peff = self.calculate_effective_pressure(rxn)
            self.kf[j] = rxn.getRateCoefficient(self.T.value_si, Peff)

            if rxn.reversible:
                self.Keq[j] = rxn.getEquilibriumConstant(self.T.value_si)
                self.kb[j] = self.kf[j] / self.Keq[j]

    def set_colliders(self, coreReactions, edgeReactions, coreSpecies):
        """
        Store collider efficiencies and reaction indices for pdep reactions that have collider efficiencies,
        and store specific collider indices
        """
        pdep_collider_reaction_indices = []
        self.pdepColliderKinetics = []
        collider_efficiencies = []
        pdep_specific_collider_reaction_indices = []
        self.pdepSpecificColliderKinetics = []
        self.specificColliderSpecies = []

        for rxn in itertools.chain(coreReactions, edgeReactions):
            if rxn.kinetics.isPressureDependent():
                if rxn.kinetics.efficiencies:
                    j = self.reactionIndex[rxn]
                    pdep_collider_reaction_indices.append(j)
                    self.pdepColliderKinetics.append(rxn.kinetics)
                    collider_efficiencies.append(rxn.kinetics.getEffectiveColliderEfficiencies(coreSpecies))
                if rxn.specific_collider:
                    pdep_specific_collider_reaction_indices.append(self.reactionIndex[rxn])
                    self.pdepSpecificColliderKinetics.append(rxn.kinetics)
                    self.specificColliderSpecies.append(rxn.specific_collider)

        self.pdepColliderReactionIndices = np.array(pdep_collider_reaction_indices, np.int)
        self.colliderEfficiencies = np.array(collider_efficiencies, np.float64)
        self.pdepSpecificColliderReactionIndices = np.array(pdep_specific_collider_reaction_indices, np.int)

    def set_initial_conditions(self):
        """
        Sets the initial conditions of the rate equations that represent the
        current reactor model.

        The volume is set to the value derived from the ideal gas law, using the
        user-defined pressure, temperature, and the number of moles of initial species.

        The species moles array (y0) is set to the values stored in the
        initial mole fractions dictionary.

        The initial species concentration is computed and stored in the
        coreSpeciesConcentrations array.

        """

        ReactionSystem.set_initial_conditions(self)

        for spec, moleFrac in self.initialMoleFractions.items():
            i = self.get_species_index(spec)
            self.y0[i] = moleFrac

        # Use ideal gas law to compute volume
        self.V = constants.R * self.T.value_si * np.sum(self.y0[:self.numCoreSpecies]) / self.P.value_si  # volume in m^3
        for j in range(self.numCoreSpecies):
            self.coreSpeciesConcentrations[j] = self.y0[j] / self.V

    @cython.boundscheck(False)
    def residual(self, double t, np.ndarray[np.float64_t, ndim=1] y, np.ndarray[np.float64_t, ndim=1] dydt,
                 np.ndarray[np.float64_t, ndim=1] senpar = np.zeros(1, np.float64)):
        """
        Return the residual function for the governing DAE system for the
        simple reaction system.
        """
        cdef np.ndarray[np.int_t, ndim=2] ir, ip, inet
        cdef np.ndarray[np.float64_t, ndim=1] res, kf, kr, knet, delta, equilibrium_constants
        cdef int num_core_species, num_core_reactions, num_edge_species, num_edge_reactions, num_pdep_networks
        cdef int i, j, z, first, second, third, real_species_index
        cdef double k, V, reaction_rate, rev_reaction_rate, T, P, Peff, core_species_rate
        cdef np.ndarray[np.float64_t, ndim=1] core_species_concentrations, core_species_rates, core_reaction_rates
        cdef np.ndarray[np.float64_t, ndim=1] edge_species_rates, edge_reaction_rates, network_leak_rates
        cdef np.ndarray[np.float64_t, ndim=1] core_species_consumption_rates, core_species_production_rates
        cdef np.ndarray[np.float64_t, ndim=1] C, y_core_species
        cdef np.ndarray[np.float64_t, ndim=2] jacobian, dgdk, collider_efficiencies
        cdef np.ndarray[np.int_t, ndim=1] pdep_collider_reaction_indices, pdep_specific_collider_reaction_indices
        cdef list pdep_collider_kinetics, pdep_specific_collider_kinetics

        ir = self.reactantIndices
        ip = self.productIndices

        num_core_species = len(self.coreSpeciesRates)
        num_core_reactions = len(self.coreReactionRates)
        num_edge_species = len(self.edgeSpeciesRates)
        num_edge_reactions = len(self.edgeReactionRates)
        num_pdep_networks = len(self.networkLeakRates)
        kf = self.kf
        kr = self.kb

        y_core_species = y[:num_core_species]

        # Recalculate any forward and reverse rate coefficients that involve pdep collision efficiencies
        if self.pdepColliderReactionIndices.shape[0] != 0:
            T = self.T.value_si
            P = self.P.value_si
            equilibrium_constants = self.Keq
            pdep_collider_reaction_indices = self.pdepColliderReactionIndices
            pdep_collider_kinetics = self.pdepColliderKinetics
            collider_efficiencies = self.colliderEfficiencies
            for i in range(pdep_collider_reaction_indices.shape[0]):
                # Calculate effective pressure
                Peff = P * np.sum(collider_efficiencies[i] * y_core_species / np.sum(y_core_species))
                j = pdep_collider_reaction_indices[i]
                kf[j] = pdep_collider_kinetics[i].getRateCoefficient(T, Peff)
                kr[j] = kf[j] / equilibrium_constants[j]
        if self.pdepSpecificColliderReactionIndices.shape[0] != 0:
            T = self.T.value_si
            P = self.P.value_si
            equilibrium_constants = self.Keq
            pdep_specific_collider_reaction_indices = self.pdepSpecificColliderReactionIndices
            pdep_specific_collider_kinetics = self.pdepSpecificColliderKinetics
            specific_collider_species = self.specificColliderSpecies
            for i in range(pdep_specific_collider_reaction_indices.shape[0]):
                # Calculate effective pressure
                Peff = P * y[self.speciesIndex[specific_collider_species[i]]] / np.sum(y_core_species)
                j = pdep_specific_collider_reaction_indices[i]
                kf[j] = pdep_specific_collider_kinetics[i].getRateCoefficient(T, Peff)
                kr[j] = kf[j] / equilibrium_constants[j]

        inet = self.networkIndices
        knet = self.networkLeakCoefficients

        res = np.zeros(num_core_species, np.float64)

        core_species_concentrations = np.zeros_like(self.coreSpeciesConcentrations)
        core_species_rates = np.zeros_like(self.coreSpeciesRates)
        core_reaction_rates = np.zeros_like(self.coreReactionRates)
        core_species_consumption_rates = np.zeros_like(self.coreSpeciesConsumptionRates)
        core_species_production_rates = np.zeros_like(self.coreSpeciesProductionRates)
        edge_species_rates = np.zeros_like(self.edgeSpeciesRates)
        edge_reaction_rates = np.zeros_like(self.edgeReactionRates)
        network_leak_rates = np.zeros_like(self.networkLeakRates)

        C = np.zeros_like(self.coreSpeciesConcentrations)

        # Use ideal gas law to compute volume
        V = constants.R * self.T.value_si * np.sum(y_core_species) / self.P.value_si
        self.V = V

        for j in range(num_core_species):
            C[j] = y[j] / V
            core_species_concentrations[j] = C[j]

        for j in range(ir.shape[0]):
            k = kf[j]
            if ir[j, 0] >= num_core_species or ir[j, 1] >= num_core_species or ir[j, 2] >= num_core_species:
                f_reaction_rate = 0.0
            elif ir[j, 1] == -1:  # only one reactant
                f_reaction_rate = k * C[ir[j, 0]]
            elif ir[j, 2] == -1:  # only two reactants
                f_reaction_rate = k * C[ir[j, 0]] * C[ir[j, 1]]
            else:  # three reactants
                f_reaction_rate = k * C[ir[j, 0]] * C[ir[j, 1]] * C[ir[j, 2]]
            k = kr[j]
            if ip[j, 0] >= num_core_species or ip[j, 1] >= num_core_species or ip[j, 2] >= num_core_species:
                rev_reaction_rate = 0.0
            elif ip[j, 1] == -1:  # only one reactant
                rev_reaction_rate = k * C[ip[j, 0]]
            elif ip[j, 2] == -1:  # only two reactants
                rev_reaction_rate = k * C[ip[j, 0]] * C[ip[j, 1]]
            else:  # three reactants
                rev_reaction_rate = k * C[ip[j, 0]] * C[ip[j, 1]] * C[ip[j, 2]]

            reaction_rate = f_reaction_rate - rev_reaction_rate

            # Set the reaction and species rates
            if j < num_core_reactions:
                # The reaction is a core reaction
                core_reaction_rates[j] = reaction_rate

                # Add/substract the total reaction rate from each species rate
                # Since it's a core reaction we know that all of its reactants
                # and products are core species
                first = ir[j, 0]
                core_species_rates[first] -= reaction_rate
                core_species_consumption_rates[first] += f_reaction_rate
                core_species_production_rates[first] += rev_reaction_rate
                second = ir[j, 1]
                if second != -1:
                    core_species_rates[second] -= reaction_rate
                    core_species_consumption_rates[second] += f_reaction_rate
                    core_species_production_rates[second] += rev_reaction_rate
                    third = ir[j, 2]
                    if third != -1:
                        core_species_rates[third] -= reaction_rate
                        core_species_consumption_rates[third] += f_reaction_rate
                        core_species_production_rates[third] += rev_reaction_rate
                first = ip[j, 0]
                core_species_rates[first] += reaction_rate
                core_species_production_rates[first] += f_reaction_rate
                core_species_consumption_rates[first] += rev_reaction_rate
                second = ip[j, 1]
                if second != -1:
                    core_species_rates[second] += reaction_rate
                    core_species_production_rates[second] += f_reaction_rate
                    core_species_consumption_rates[second] += rev_reaction_rate
                    third = ip[j, 2]
                    if third != -1:
                        core_species_rates[third] += reaction_rate
                        core_species_production_rates[third] += f_reaction_rate
                        core_species_consumption_rates[third] += rev_reaction_rate

            else:
                # The reaction is an edge reaction
                edge_reaction_rates[j - num_core_reactions] = reaction_rate

                # Add/substract the total reaction rate from each species rate
                # Since it's an edge reaction its reactants and products could
                # be either core or edge species
                # We're only interested in the edge species
                first = ir[j, 0]
                if first >= num_core_species: edge_species_rates[first - num_core_species] -= reaction_rate
                second = ir[j, 1]
                if second != -1:
                    if second >= num_core_species: edge_species_rates[second - num_core_species] -= reaction_rate
                    third = ir[j, 2]
                    if third != -1:
                        if third >= num_core_species: edge_species_rates[third - num_core_species] -= reaction_rate
                first = ip[j, 0]
                if first >= num_core_species: edge_species_rates[first - num_core_species] += reaction_rate
                second = ip[j, 1]
                if second != -1:
                    if second >= num_core_species: edge_species_rates[second - num_core_species] += reaction_rate
                    third = ip[j, 2]
                    if third != -1:
                        if third >= num_core_species: edge_species_rates[third - num_core_species] += reaction_rate

        for j in range(inet.shape[0]):
            k = knet[j]
            if inet[j, 1] == -1:  # only one reactant
                reaction_rate = k * C[inet[j, 0]]
            elif inet[j, 2] == -1:  # only two reactants
                reaction_rate = k * C[inet[j, 0]] * C[inet[j, 1]]
            else:  # three reactants
                reaction_rate = k * C[inet[j, 0]] * C[inet[j, 1]] * C[inet[j, 2]]
            network_leak_rates[j] = reaction_rate

        # Add unimolecular reactions describing the transport of sampled molecules to the measurement zone
        # (i.e., ionization region of a Mass Spectrometer, intersection with probing laser, etc.)
        for j, core_species_rate in enumerate(core_species_rates):
            species = list(self.speciesIndex.keys())[list(self.speciesIndex.values()).index(j)]
            if '_obs' in species.label:
                for species2 in self.speciesIndex.keys():
                    if species2.label == species.label.replace('_obs', ''):
                        real_species_index = self.speciesIndex[species2]
                        if real_species_index > len(C) - 1:
                            continue
                        else:
                            core_species_rates[j] = self.k_sampling.value_si * C[real_species_index] - self.k_sampling.value_si * C[j]
                            core_species_consumption_rates[j] = self.k_sampling.value_si * C[j]
                            core_species_production_rates[j] = self.k_sampling.value_si * C[real_species_index]
                            break

        self.coreSpeciesConcentrations = core_species_concentrations
        self.coreSpeciesRates = core_species_rates
        self.coreSpeciesProductionRates = core_species_production_rates
        self.coreSpeciesConsumptionRates = core_species_consumption_rates
        self.coreReactionRates = core_reaction_rates
        self.edgeSpeciesRates = edge_species_rates
        self.edgeReactionRates = edge_reaction_rates
        self.networkLeakRates = network_leak_rates

        res = core_species_rates * V

        delta = res - dydt

        # Return DELTA, IRES.  IRES is set to 1 in order to tell DASPK to evaluate the sensitivity residuals
        return delta, 1
