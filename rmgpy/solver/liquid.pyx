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
Contains the :class:`LiquidReactor` class, providing a reaction system
consisting of a homogeneous, isothermal, isobaric batch reactor.
"""

import itertools

cimport cython
import numpy as np
cimport numpy as np

import rmgpy.constants as constants
cimport rmgpy.constants as constants
from rmgpy.quantity import Quantity
from rmgpy.quantity cimport ScalarQuantity, ArrayQuantity
from rmgpy.solver.base cimport ReactionSystem


cdef class LiquidReactor(ReactionSystem):
    """
    A reaction system consisting of a homogeneous, isothermal, constant volume batch
    reactor. These assumptions allow for a number of optimizations that enable
    this solver to complete very rapidly, even for large kinetic models.
    """

    cdef public ScalarQuantity T
    cdef public ScalarQuantity P
    cdef public double V
    cdef public bint constantVolume
    cdef public double viscosity
    cdef public list constSPCNames
    cdef public list constSPCIndices
    cdef public dict initialConcentrations
    cdef public list Trange
    cdef public int nSims
    cdef public dict sensConditions

    def __init__(self, T, initialConcentrations, nSims=1, termination=None, sensitiveSpecies=None,
                 sensitivityThreshold=1e-3, sensConditions=None, constSPCNames=None):

        ReactionSystem.__init__(self, termination, sensitiveSpecies, sensitivityThreshold)

        if type(T) != list:
            self.T = Quantity(T)
        else:
            self.Trange = [Quantity(t) for t in T]

        self.P = Quantity(100000., 'kPa')  # Arbitrary high pressure (1000 Bar) to get reactions in the high-pressure limit!
        self.initialConcentrations = initialConcentrations  # should be passed in SI
        self.V = 0  # will be set from initialConcentrations in initializeModel
        self.constantVolume = True
        self.viscosity = 0  # in Pa*s

        #Constant concentration attributes
        self.constSPCIndices = None
        self.constSPCNames = constSPCNames  #store index of constant species
        self.sensConditions = sensConditions
        self.nSims = nSims

    def convertInitialKeysToSpeciesObjects(self, speciesDict):
        """
        Convert the initial_concentrations dictionary from species names into species objects,
        using the given dictionary of species.
        """
        initial_concentrations = {}
        for label, moleFrac in self.initialConcentrations.items():
            if label == 'T':
                continue
            initial_concentrations[speciesDict[label]] = moleFrac
        self.initialConcentrations = initial_concentrations

        conditions = {}
        if self.sensConditions is not None:
            for label, value in self.sensConditions.items():
                if label == 'T':
                    conditions[label] = value
                else:
                    conditions[speciesDict[label]] = value
        self.sensConditions = conditions

    def get_constSPCIndices(self, coreSpecies):
        """Allow to identify constant Species position in solver"""
        for spc in self.constSPCNames:
            if self.constSPCIndices is None:  # initialize once the list if constant SPC declared
                self.constSPCIndices = []
            for item in coreSpecies:
                # Need to identify the species object corresponding to the the string written in the input file
                if item.label == spc:
                    self.constSPCIndices.append(coreSpecies.index(item))

    cpdef initializeModel(self, list coreSpecies, list coreReactions, list edgeSpecies, list edgeReactions,
                          list surfaceSpecies=None, list surfaceReactions=None, list pdepNetworks=None,
                          atol=1e-16, rtol=1e-8, sensitivity=False, sens_atol=1e-6, sens_rtol=1e-4,
                          filterReactions=False, dict conditions=None):
        """
        Initialize a simulation of the liquid reactor using the provided kinetic
        model.
        """
        if surfaceSpecies is None:
            surfaceSpecies = []
        if surfaceReactions is None:
            surfaceReactions = []

        # First call the base class version of the method
        # This initializes the attributes declared in the base class
        ReactionSystem.initializeModel(self, coreSpecies, coreReactions, edgeSpecies, edgeReactions,
                                       surfaceSpecies=surfaceSpecies, surfaceReactions=surfaceReactions,
                                       pdepNetworks=pdepNetworks, atol=atol, rtol=rtol,
                                       sensitivity=sensitivity, sens_atol=sens_atol, sens_rtol=sens_rtol,
                                       filterReactions=filterReactions, conditions=conditions)

        # Set initial conditions
        self.set_initial_conditions()

        # Compute reaction thresholds if reaction filtering is turned on
        if filterReactions:
            ReactionSystem.set_initial_reaction_thresholds(self)

        # Generate forward and reverse rate coefficients k(T,P)
        self.generate_rate_coefficients(coreReactions, edgeReactions)

        ReactionSystem.compute_network_variables(self, pdepNetworks)

        ReactionSystem.set_initial_derivative(self)

        # Initialize the model
        ReactionSystem.initialize_solver(self)

    def generate_rate_coefficients(self, coreReactions, edgeReactions):
        """
        Populates the forwardRateCoefficients, reverseRateCoefficients and equilibriumConstants
        arrays with the values computed at the temperature and (effective) pressure of the 
        reacion system.
        """

        for rxn in itertools.chain(coreReactions, edgeReactions):
            j = self.reactionIndex[rxn]
            self.kf[j] = rxn.get_rate_coefficient(self.T.value_si, self.P.value_si)
            if rxn.reversible:
                self.Keq[j] = rxn.get_equilibrium_constant(self.T.value_si)
                self.kb[j] = self.kf[j] / self.Keq[j]

    def get_threshold_rate_constants(self, modelSettings):
        """
        Get the threshold rate constants for reaction filtering.

        modelSettings is not used here, but is needed so that the method
        matches the one in simpleReactor.
        """
        # Set the maximum unimolecular rate to be kB*T/h
        unimolecular_threshold_rate_constant = 2.08366122e10 * self.T.value_si
        # Set the maximum bi/trimolecular rates based on the Smoluchowski and Stokes-Einstein equations
        bimolecular_threshold_rate_constant = 22.2 * self.T.value_si / self.viscosity
        trimolecular_threshold_rate_constant = 0.11 * self.T.value_si / self.viscosity
        return (unimolecular_threshold_rate_constant,
                bimolecular_threshold_rate_constant,
                trimolecular_threshold_rate_constant)

    def set_initial_conditions(self):
        """
        Sets the initial conditions of the rate equations that represent the 
        current reactor model.

        The volume is set to the value in m3 required to contain 
        one mole total of core species at start.

        The coreSpeciesConcentrations array is set to the values stored in the
        initial concentrations dictionary.

        The initial number of moles of a species j is computed and stored in the
        y0 instance attribute.

        """
        ReactionSystem.set_initial_conditions(self)

        for spec, conc in self.initialConcentrations.items():
            i = self.get_species_index(spec)
            self.coreSpeciesConcentrations[i] = conc

        V = 1.0 / np.sum(self.coreSpeciesConcentrations)
        self.V = V

        for j in range(self.numCoreSpecies):
            self.y0[j] = self.coreSpeciesConcentrations[j] * V

    @cython.boundscheck(False)
    def residual(self, double t, np.ndarray[np.float64_t, ndim=1] y, np.ndarray[np.float64_t, ndim=1] dydt,
                 np.ndarray[np.float64_t, ndim=1] senpar = np.zeros(1, np.float64)):

        """
        Return the residual function for the governing DAE system for the
        liquid reaction system.
        """
        cdef np.ndarray[np.int_t, ndim=2] ir, ip, inet
        cdef np.ndarray[np.float64_t, ndim=1] res, kf, kr, knet, delta, equilibrium_constants
        cdef int num_core_species, num_core_reactions, num_edge_species, num_edge_reactions, num_pdep_networks
        cdef int i, j, z, first, second, third
        cdef double k, V, reaction_rate
        cdef np.ndarray[np.float64_t,ndim=1] core_species_concentrations, core_species_rates, core_reaction_rates
        cdef np.ndarray[np.float64_t,ndim=1] edge_species_rates, edge_reaction_rates, network_leak_rates
        cdef np.ndarray[np.float64_t,ndim=1] core_species_consumption_rates, core_species_production_rates
        cdef np.ndarray[np.float64_t, ndim=1] C
        cdef np.ndarray[np.float64_t, ndim=2] jacobian, dgdk

        ir = self.reactantIndices
        ip = self.productIndices
        equilibrium_constants = self.Keq

        kf = self.kf
        kr = self.kb

        inet = self.networkIndices
        knet = self.networkLeakCoefficients

        num_core_species = len(self.coreSpeciesRates)
        num_core_reactions = len(self.coreReactionRates)
        num_edge_species = len(self.edgeSpeciesRates)
        num_edge_reactions = len(self.edgeReactionRates)
        num_pdep_networks = len(self.networkLeakRates)

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
        V = self.V  # constant volume reactor

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

        # chatelak: Same as in Java, core species rate = 0 if declared as constant
        if self.constSPCIndices is not None:
            for spc_index in self.constSPCIndices:
                core_species_rates[spc_index] = 0

        self.coreSpeciesConcentrations = core_species_concentrations
        self.coreSpeciesRates = core_species_rates
        self.coreReactionRates = core_reaction_rates
        self.coreSpeciesProductionRates = core_species_production_rates
        self.coreSpeciesConsumptionRates = core_species_consumption_rates
        self.edgeSpeciesRates = edge_species_rates
        self.edgeReactionRates = edge_reaction_rates
        self.networkLeakRates = network_leak_rates

        res = core_species_rates * V

        if self.sensitivity:
            delta = np.zeros(len(y), np.float64)
            delta[:num_core_species] = res
            if self.jacobianMatrix is None:
                jacobian = self.jacobian(t, y, dydt, 0, senpar)
            else:
                jacobian = self.jacobianMatrix
            dgdk = ReactionSystem.computeRateDerivative(self)
            for j in range(num_core_reactions + num_core_species):
                for i in range(num_core_species):
                    for z in range(num_core_species):
                        delta[(j + 1) * num_core_species + i] += jacobian[i, z] * y[(j + 1) * num_core_species + z]
                    delta[(j + 1) * num_core_species + i] += dgdk[i, j]

        else:
            delta = res
        delta = delta - dydt

        # Return DELTA, IRES.  IRES is set to 1 in order to tell DASPK to evaluate the sensitivity residuals
        return delta, 1

    @cython.boundscheck(False)
    def jacobian(self, double t, np.ndarray[np.float64_t, ndim=1] y, np.ndarray[np.float64_t, ndim=1] dydt,
                 double cj, np.ndarray[np.float64_t, ndim=1] senpar = np.zeros(1, np.float64)):
        """
        Return the analytical Jacobian for the reaction system.
        """
        cdef np.ndarray[np.int_t, ndim=2] ir, ip
        cdef np.ndarray[np.float64_t, ndim=1] kf, kr, C
        cdef np.ndarray[np.float64_t, ndim=2] pd
        cdef int num_core_reactions, num_core_species, i, j
        cdef double k, V, Ctot, deriv, corr

        ir = self.reactantIndices
        ip = self.productIndices

        kf = self.kf
        kr = self.kb
        num_core_reactions = len(self.coreReactionRates)
        num_core_species = len(self.coreSpeciesConcentrations)

        pd = -cj * np.identity(num_core_species, np.float64)

        V = self.V  # volume is constant

        C = np.zeros_like(self.coreSpeciesConcentrations)
        for j in range(num_core_species):
            C[j] = y[j] / V

        for j in range(num_core_reactions):

            k = kf[j]
            if ir[j, 1] == -1:  # only one reactant
                deriv = k
                pd[ir[j, 0], ir[j, 0]] -= deriv

                pd[ip[j, 0], ir[j, 0]] += deriv
                if ip[j, 1] != -1:
                    pd[ip[j, 1], ir[j, 0]] += deriv
                    if ip[j, 2] != -1:
                        pd[ip[j, 2], ir[j, 0]] += deriv


            elif ir[j, 2] == -1:  # only two reactants
                if ir[j, 0] == ir[j, 1]:  # reactants are the same
                    deriv = 2 * k * C[ir[j, 0]]
                    pd[ir[j, 0], ir[j, 0]] -= 2 * deriv

                    pd[ip[j, 0], ir[j, 0]] += deriv
                    if ip[j, 1] != -1:
                        pd[ip[j, 1], ir[j, 0]] += deriv
                        if ip[j, 2] != -1:
                            pd[ip[j, 2], ir[j, 0]] += deriv

                else:
                    # Derivative with respect to reactant 1
                    deriv = k * C[ir[j, 1]]
                    pd[ir[j, 0], ir[j, 0]] -= deriv
                    pd[ir[j, 1], ir[j, 0]] -= deriv

                    pd[ip[j, 0], ir[j, 0]] += deriv
                    if ip[j, 1] != -1:
                        pd[ip[j, 1], ir[j, 0]] += deriv
                        if ip[j, 2] != -1:
                            pd[ip[j, 2], ir[j, 0]] += deriv

                    # Derivative with respect to reactant 2
                    deriv = k * C[ir[j, 0]]
                    pd[ir[j, 0], ir[j, 1]] -= deriv
                    pd[ir[j, 1], ir[j, 1]] -= deriv

                    pd[ip[j, 0], ir[j, 1]] += deriv
                    if ip[j, 1] != -1:
                        pd[ip[j, 1], ir[j, 1]] += deriv
                        if ip[j, 2] != -1:
                            pd[ip[j, 2], ir[j, 1]] += deriv


            else:  # three reactants!! (really?)
                if (ir[j, 0] == ir[j, 1] & ir[j, 0] == ir[j, 2]):
                    deriv = 3 * k * C[ir[j, 0]] * C[ir[j, 0]]
                    pd[ir[j, 0], ir[j, 0]] -= 3 * deriv

                    pd[ip[j, 0], ir[j, 0]] += deriv
                    if ip[j, 1] != -1:
                        pd[ip[j, 1], ir[j, 0]] += deriv
                        if ip[j, 2] != -1:
                            pd[ip[j, 2], ir[j, 0]] += deriv

                elif ir[j, 0] == ir[j, 1]:
                    # derivative with respect to reactant 1
                    deriv = 2 * k * C[ir[j, 0]] * C[ir[j, 2]]
                    pd[ir[j, 0], ir[j, 0]] -= 2 * deriv
                    pd[ir[j, 2], ir[j, 0]] -= deriv

                    pd[ip[j, 0], ir[j, 0]] += deriv
                    if ip[j, 1] != -1:
                        pd[ip[j, 1], ir[j, 0]] += deriv
                        if ip[j, 2] != -1:
                            pd[ip[j, 2], ir[j, 0]] += deriv

                    # derivative with respect to reactant 3
                    deriv = k * C[ir[j, 0]] * C[ir[j, 0]]
                    pd[ir[j, 0], ir[j, 2]] -= 2 * deriv
                    pd[ir[j, 2], ir[j, 2]] -= deriv

                    pd[ip[j, 0], ir[j, 2]] += deriv
                    if ip[j, 1] != -1:
                        pd[ip[j, 1], ir[j, 2]] += deriv
                        if ip[j, 2] != -1:
                            pd[ip[j, 2], ir[j, 2]] += deriv


                elif ir[j, 1] == ir[j, 2]:
                    # derivative with respect to reactant 1
                    deriv = k * C[ir[j, 1]] * C[ir[j, 1]]
                    pd[ir[j, 0], ir[j, 0]] -= deriv
                    pd[ir[j, 1], ir[j, 0]] -= 2 * deriv

                    pd[ip[j, 0], ir[j, 0]] += deriv
                    if ip[j, 1] != -1:
                        pd[ip[j, 1], ir[j, 0]] += deriv
                        if ip[j, 2] != -1:
                            pd[ip[j, 2], ir[j, 0]] += deriv
                    # derivative with respect to reactant 2
                    deriv = 2 * k * C[ir[j, 0]] * C[ir[j, 1]]
                    pd[ir[j, 0], ir[j, 1]] -= deriv
                    pd[ir[j, 1], ir[j, 1]] -= 2 * deriv

                    pd[ip[j, 0], ir[j, 1]] += deriv
                    if ip[j, 1] != -1:
                        pd[ip[j, 1], ir[j, 1]] += deriv
                        if ip[j, 2] != -1:
                            pd[ip[j, 2], ir[j, 1]] += deriv

                elif ir[j, 0] == ir[j, 2]:
                    # derivative with respect to reactant 1
                    deriv = 2 * k * C[ir[j, 0]] * C[ir[j, 1]]
                    pd[ir[j, 0], ir[j, 0]] -= 2 * deriv
                    pd[ir[j, 1], ir[j, 0]] -= deriv

                    pd[ip[j, 0], ir[j, 0]] += deriv
                    if ip[j, 1] != -1:
                        pd[ip[j, 1], ir[j, 0]] += deriv
                        if ip[j, 2] != -1:
                            pd[ip[j, 2], ir[j, 0]] += deriv
                    # derivative with respect to reactant 2
                    deriv = k * C[ir[j, 0]] * C[ir[j, 0]]
                    pd[ir[j, 0], ir[j, 1]] -= 2 * deriv
                    pd[ir[j, 1], ir[j, 1]] -= deriv

                    pd[ip[j, 0], ir[j, 1]] += deriv
                    if ip[j, 1] != -1:
                        pd[ip[j, 1], ir[j, 1]] += deriv
                        if ip[j, 2] != -1:
                            pd[ip[j, 2], ir[j, 1]] += deriv

                else:
                    # derivative with respect to reactant 1
                    deriv = k * C[ir[j, 1]] * C[ir[j, 2]]
                    pd[ir[j, 0], ir[j, 0]] -= deriv
                    pd[ir[j, 1], ir[j, 0]] -= deriv
                    pd[ir[j, 2], ir[j, 0]] -= deriv

                    pd[ip[j, 0], ir[j, 0]] += deriv
                    if ip[j, 1] != -1:
                        pd[ip[j, 1], ir[j, 0]] += deriv
                        if ip[j, 2] != -1:
                            pd[ip[j, 2], ir[j, 0]] += deriv

                    # derivative with respect to reactant 2
                    deriv = k * C[ir[j, 0]] * C[ir[j, 2]]
                    pd[ir[j, 0], ir[j, 1]] -= deriv
                    pd[ir[j, 1], ir[j, 1]] -= deriv
                    pd[ir[j, 2], ir[j, 1]] -= deriv

                    pd[ip[j, 0], ir[j, 1]] += deriv
                    if ip[j, 1] != -1:
                        pd[ip[j, 1], ir[j, 1]] += deriv
                        if ip[j, 2] != -1:
                            pd[ip[j, 2], ir[j, 1]] += deriv

                    # derivative with respect to reactant 3
                    deriv = k * C[ir[j, 0]] * C[ir[j, 1]]
                    pd[ir[j, 0], ir[j, 2]] -= deriv
                    pd[ir[j, 1], ir[j, 2]] -= deriv
                    pd[ir[j, 2], ir[j, 2]] -= deriv

                    pd[ip[j, 0], ir[j, 2]] += deriv
                    if ip[j, 1] != -1:
                        pd[ip[j, 1], ir[j, 2]] += deriv
                        if ip[j, 2] != -1:
                            pd[ip[j, 2], ir[j, 2]] += deriv

            k = kr[j]
            if ip[j, 1] == -1:  # only one product
                deriv = k
                pd[ip[j, 0], ip[j, 0]] -= deriv

                pd[ir[j, 0], ip[j, 0]] += deriv
                if ir[j, 1] != -1:
                    pd[ir[j, 1], ip[j, 0]] += deriv
                    if ir[j, 2] != -1:
                        pd[ir[j, 2], ip[j, 0]] += deriv


            elif ip[j, 2] == -1:  # only two products
                if ip[j, 0] == ip[j, 1]:
                    deriv = 2 * k * C[ip[j, 0]]
                    pd[ip[j, 0], ip[j, 0]] -= 2 * deriv

                    pd[ir[j, 0], ip[j, 0]] += deriv
                    if ir[j, 1] != -1:
                        pd[ir[j, 1], ip[j, 0]] += deriv
                        if ir[j, 2] != -1:
                            pd[ir[j, 2], ip[j, 0]] += deriv

                else:
                    # Derivative with respect to product 1
                    deriv = k * C[ip[j, 1]]
                    pd[ip[j, 0], ip[j, 0]] -= deriv
                    pd[ip[j, 1], ip[j, 0]] -= deriv

                    pd[ir[j, 0], ip[j, 0]] += deriv
                    if ir[j, 1] != -1:
                        pd[ir[j, 1], ip[j, 0]] += deriv
                        if ir[j, 2] != -1:
                            pd[ir[j, 2], ip[j, 0]] += deriv

                    # Derivative with respect to product 2
                    deriv = k * C[ip[j, 0]]
                    pd[ip[j, 0], ip[j, 1]] -= deriv
                    pd[ip[j, 1], ip[j, 1]] -= deriv

                    pd[ir[j, 0], ip[j, 1]] += deriv
                    if ir[j, 1] != -1:
                        pd[ir[j, 1], ip[j, 1]] += deriv
                        if ir[j, 2] != -1:
                            pd[ir[j, 2], ip[j, 1]] += deriv


            else:  # three products
                if (ip[j, 0] == ip[j, 1] & ip[j, 0] == ip[j, 2]):
                    deriv = 3 * k * C[ip[j, 0]] * C[ip[j, 0]]
                    pd[ip[j, 0], ip[j, 0]] -= 3 * deriv

                    pd[ir[j, 0], ip[j, 0]] += deriv
                    if ir[j, 1] != -1:
                        pd[ir[j, 1], ip[j, 0]] += deriv
                        if ir[j, 2] != -1:
                            pd[ir[j, 2], ip[j, 0]] += deriv

                elif ip[j, 0] == ip[j, 1]:
                    # derivative with respect to product 1
                    deriv = 2 * k * C[ip[j, 0]] * C[ip[j, 2]]
                    pd[ip[j, 0], ip[j, 0]] -= 2 * deriv
                    pd[ip[j, 2], ip[j, 0]] -= deriv

                    pd[ir[j, 0], ip[j, 0]] += deriv
                    if ir[j, 1] != -1:
                        pd[ir[j, 1], ip[j, 0]] += deriv
                        if ir[j, 2] != -1:
                            pd[ir[j, 2], ip[j, 0]] += deriv
                    # derivative with respect to product 3
                    deriv = k * C[ip[j, 0]] * C[ip[j, 0]]
                    pd[ip[j, 0], ip[j, 2]] -= 2 * deriv
                    pd[ip[j, 2], ip[j, 2]] -= deriv

                    pd[ir[j, 0], ip[j, 2]] += deriv
                    if ir[j, 1] != -1:
                        pd[ir[j, 1], ip[j, 2]] += deriv
                        if ir[j, 2] != -1:
                            pd[ir[j, 2], ip[j, 2]] += deriv


                elif ip[j, 1] == ip[j, 2]:
                    # derivative with respect to product 1
                    deriv = k * C[ip[j, 1]] * C[ip[j, 1]]
                    pd[ip[j, 0], ip[j, 0]] -= deriv
                    pd[ip[j, 1], ip[j, 0]] -= 2 * deriv

                    pd[ir[j, 0], ip[j, 0]] += deriv
                    if ir[j, 1] != -1:
                        pd[ir[j, 1], ip[j, 0]] += deriv
                        if ir[j, 2] != -1:
                            pd[ir[j, 2], ip[j, 0]] += deriv

                    # derivative with respect to product 2
                    deriv = 2 * k * C[ip[j, 0]] * C[ip[j, 1]]
                    pd[ip[j, 0], ip[j, 1]] -= deriv
                    pd[ip[j, 1], ip[j, 1]] -= 2 * deriv

                    pd[ir[j, 0], ip[j, 1]] += deriv
                    if ir[j, 1] != -1:
                        pd[ir[j, 1], ip[j, 1]] += deriv
                        if ir[j, 2] != -1:
                            pd[ir[j, 2], ip[j, 1]] += deriv


                elif ip[j, 0] == ip[j, 2]:
                    # derivative with respect to product 1
                    deriv = 2 * k * C[ip[j, 0]] * C[ip[j, 1]]
                    pd[ip[j, 0], ip[j, 0]] -= 2 * deriv
                    pd[ip[j, 1], ip[j, 0]] -= deriv

                    pd[ir[j, 0], ip[j, 0]] += deriv
                    if ir[j, 1] != -1:
                        pd[ir[j, 1], ip[j, 0]] += deriv
                        if ir[j, 2] != -1:
                            pd[ir[j, 2], ip[j, 0]] += deriv
                    # derivative with respect to product 2
                    deriv = k * C[ip[j, 0]] * C[ip[j, 0]]
                    pd[ip[j, 0], ip[j, 1]] -= 2 * deriv
                    pd[ip[j, 1], ip[j, 1]] -= deriv

                    pd[ir[j, 0], ip[j, 1]] += deriv
                    if ir[j, 1] != -1:
                        pd[ir[j, 1], ip[j, 1]] += deriv
                        if ir[j, 2] != -1:
                            pd[ir[j, 2], ip[j, 1]] += deriv

                else:
                    # derivative with respect to product 1
                    deriv = k * C[ip[j, 1]] * C[ip[j, 2]]
                    pd[ip[j, 0], ip[j, 0]] -= deriv
                    pd[ip[j, 1], ip[j, 0]] -= deriv
                    pd[ip[j, 2], ip[j, 0]] -= deriv

                    pd[ir[j, 0], ip[j, 0]] += deriv
                    if ir[j, 1] != -1:
                        pd[ir[j, 1], ip[j, 0]] += deriv
                        if ir[j, 2] != -1:
                            pd[ir[j, 2], ip[j, 0]] += deriv

                    # derivative with respect to product 2
                    deriv = k * C[ip[j, 0]] * C[ip[j, 2]]
                    pd[ip[j, 0], ip[j, 1]] -= deriv
                    pd[ip[j, 1], ip[j, 1]] -= deriv
                    pd[ip[j, 2], ip[j, 1]] -= deriv

                    pd[ir[j, 0], ip[j, 1]] += deriv
                    if ir[j, 1] != -1:
                        pd[ir[j, 1], ip[j, 1]] += deriv
                        if ir[j, 2] != -1:
                            pd[ir[j, 2], ip[j, 1]] += deriv

                    # derivative with respect to product 3
                    deriv = k * C[ip[j, 0]] * C[ip[j, 1]]
                    pd[ip[j, 0], ip[j, 2]] -= deriv
                    pd[ip[j, 1], ip[j, 2]] -= deriv
                    pd[ip[j, 2], ip[j, 2]] -= deriv

                    pd[ir[j, 0], ip[j, 2]] += deriv
                    if ir[j, 1] != -1:
                        pd[ir[j, 1], ip[j, 2]] += deriv
                        if ir[j, 2] != -1:
                            pd[ir[j, 2], ip[j, 2]] += deriv

        self.jacobianMatrix = pd + cj * np.identity(num_core_species, np.float64)
        return pd
