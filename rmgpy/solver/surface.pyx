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
Contains the :class:`SimpleReactor` class, providing a reaction system
consisting of a homogeneous, isothermal, isobaric batch reactor.
"""

import itertools
import logging

cimport cython
import numpy as np
cimport numpy as np

import rmgpy.constants as constants
cimport rmgpy.constants as constants
from rmgpy.quantity import Quantity
from rmgpy.quantity cimport ScalarQuantity
from rmgpy.solver.base cimport ReactionSystem


cdef class SurfaceReactor(ReactionSystem):
    """
    A reaction system consisting of a heterogeneous, isothermal, constant volume batch
    reactor. 
    """

    cdef public ScalarQuantity T
    cdef public ScalarQuantity initialP
    cdef public double V
    cdef public bint constantVolume

    cdef public list Trange
    cdef public list Prange
    cdef public int nSims
    cdef public dict sensConditions

    cdef public dict initialGasMoleFractions
    cdef public dict initialSurfaceCoverages
    cdef public ScalarQuantity surfaceVolumeRatio
    cdef public ScalarQuantity surfaceSiteDensity
    cdef public np.ndarray reactionsOnSurface  # (catalyst surface, not core/edge surface)
    cdef public np.ndarray speciesOnSurface  # (catalyst surface, not core/edge surface)

    def __init__(self,
                 T,
                 initialP,
                 initialGasMoleFractions,
                 initialSurfaceCoverages,
                 surfaceVolumeRatio,
                 surfaceSiteDensity,
                 nSims=None,
                 termination=None,
                 sensitiveSpecies=None,
                 sensitivityThreshold=1e-3,
                 sensConditions=None,
                 ):
        ReactionSystem.__init__(self,
                                termination,
                                sensitiveSpecies,
                                sensitivityThreshold)

        if isinstance(T, list):
            self.Trange = [Quantity(t) for t in T]
        else:
            self.T = Quantity(T)
        if isinstance(initialP, list):
            raise NotImplementedError("Can't do ranges of initial pressures for surface reactors yet")
        else:
            self.initialP = Quantity(initialP)
        self.initialGasMoleFractions = initialGasMoleFractions
        self.initialSurfaceCoverages = initialSurfaceCoverages
        self.surfaceVolumeRatio = Quantity(surfaceVolumeRatio)
        self.surfaceSiteDensity = Quantity(surfaceSiteDensity)
        self.V = 0  # will be set from ideal gas law in initializeModel
        self.constantVolume = True
        self.sensConditions = sensConditions
        self.nSims = nSims

    def convertInitialKeysToSpeciesObjects(self, speciesDict):
        """
        Convert the initial_gas_mole_fractions and initial_surface_coverages dictionaries
        from species names into species objects,
        using the given dictionary of species.
        """
        initial_gas_mole_fractions = {}
        for label, moleFrac in self.initialGasMoleFractions.items():
            initial_gas_mole_fractions[speciesDict[label]] = moleFrac
        self.initialGasMoleFractions = initial_gas_mole_fractions
        initial_surface_coverages = {}
        for label, surfaceCoverage in self.initialSurfaceCoverages.items():
            initial_surface_coverages[speciesDict[label]] = surfaceCoverage
        self.initialSurfaceCoverages = initial_surface_coverages

    cpdef initializeModel(self,
                          list coreSpecies,
                          list coreReactions,
                          list edgeSpecies,
                          list edgeReactions,
                          list surfaceSpecies=[],
                          list surfaceReactions=[],
                          list pdepNetworks=None,
                          atol=1e-16,
                          rtol=1e-8,
                          sensitivity=False,
                          sens_atol=1e-6,
                          sens_rtol=1e-4,
                          filterReactions=False,
                          dict conditions=None,
                          ):
        """
        Initialize a simulation of the simple reactor using the provided kinetic
        model.
        """

        # First call the base class version of the method
        # This initializes the attributes declared in the base class
        ReactionSystem.initializeModel(self,
                                       coreSpecies=coreSpecies,
                                       coreReactions=coreReactions,
                                       edgeSpecies=edgeSpecies,
                                       edgeReactions=edgeReactions,
                                       surfaceSpecies=surfaceSpecies,
                                       surfaceReactions=surfaceReactions,
                                       pdepNetworks=pdepNetworks,
                                       atol=atol,
                                       rtol=rtol,
                                       sensitivity=sensitivity,
                                       sens_atol=sens_atol,
                                       sens_rtol=sens_rtol,
                                       filterReactions=filterReactions,
                                       conditions=conditions,
                                       )
        cdef np.ndarray[np.int_t, ndim=1] species_on_surface, reactions_on_surface
        cdef int index
        #: 1 if it's on a surface, 0 if it's in the gas phase
        reactions_on_surface = np.zeros((self.numCoreReactions + self.numEdgeReactions), np.int)
        species_on_surface = np.zeros((self.numCoreSpecies), np.int)
        for spec, index in self.speciesIndex.items():
            if index >= self.numCoreSpecies:
                continue
            if spec.containsSurfaceSite():
                species_on_surface[index] = 1
        for rxn, index in self.reactionIndex.items():
            if rxn.isSurfaceReaction():
                reactions_on_surface[index] = 1
        self.speciesOnSurface = species_on_surface
        self.reactionsOnSurface = reactions_on_surface

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
        Populates the kf, kb and equilibriumConstants
        arrays with the values computed at the temperature and (effective) pressure of the 
        reaction system.
        """

        cdef double P, surface_volume_ratio_si

        surface_volume_ratio_si = self.surfaceVolumeRatio.value_si
        # ToDo: Pressure should come from ideal gas law?
        P = self.initialP.value_si

        warned = False
        for rxn in itertools.chain(coreReactions, edgeReactions):
            j = self.reactionIndex[rxn]

            # ToDo: getRateCoefficient should also depend on surface coverages vector

            if rxn.isSurfaceReaction():
                """
                Be careful! From here on kf and kb will now be in Volume units,
                even for surface reactions (which you may expect to be in Area units).
                This is to avoid repeatedly multiplying a bunch of things inside every 
                loop of the ODE solver.
                """
                self.kf[j] = (surface_volume_ratio_si *
                              rxn.getSurfaceRateCoefficient(self.T.value_si,
                                                            self.surfaceSiteDensity.value_si
                                                            ))
            else:
                if not warned and rxn.kinetics.is_pressure_dependent():
                    logging.warning("Pressure may be varying, but using initial pressure to evaluate k(T,P) expressions!")
                    warned = True
                self.kf[j] = rxn.getRateCoefficient(self.T.value_si, P)
            if rxn.reversible:
                # ToDo: getEquilibriumConstant should be coverage dependent
                self.Keq[j] = rxn.getEquilibriumConstant(self.T.value_si)
                self.kb[j] = self.kf[j] / self.Keq[j]

    def log_initial_conditions(self, number=None):
        """
        Log to the console some information about this reaction system.
        
        Should correspond to the calculations done in set_initial_conditions.
        """
        logging.info("\nSurface reaction system {}".format(number if number is not None else ""))
        logging.info("Gas phase mole fractions:")
        total_gas_moles = 0
        for spec, moleFrac in self.initialGasMoleFractions.items():
            logging.info("  {0:20s} {1:.5g}".format(spec, moleFrac))
            total_gas_moles += moleFrac
        logging.info("Total gas phase:          {:.3g} moles".format(total_gas_moles))
        logging.info("Pressure:                 {:.3g} Pa".format(self.initialP.value_si))
        logging.info("Temperature:              {} K".format(self.T.value_si))
        V = constants.R * self.T.value_si * total_gas_moles / self.initialP.value_si
        logging.info("Reactor volume:           {:.3g} m3".format(V))
        surface_volume_ratio_si = self.surfaceVolumeRatio.value_si  # 1/m
        logging.info("Surface/volume ratio:     {:.3g} m2/m3".format(surface_volume_ratio_si))
        logging.info("Surface site density:     {:.3g} mol/m2".format(self.surfaceSiteDensity.value_si))
        total_surface_sites = V * surface_volume_ratio_si * self.surfaceSiteDensity.value_si  # total surface sites in reactor
        logging.info("Surface sites in reactor: {:.3g} moles".format(total_surface_sites))
        logging.info("Initial surface coverages (and amounts):")
        for spec, coverage in self.initialSurfaceCoverages.items():
            logging.info("  {:18s} {:.5g} = {:.5g} moles".format(spec, coverage, total_surface_sites * coverage))

    def set_initial_conditions(self):
        """
        Sets the initial conditions of the rate equations that represent the 
        current reactor model.

        The volume is set to the value in m3 required to contain 
        one mole total of gas phase core species at start.

        The total surface sites are calculated from surfaceVolumeRatio and surfaceSiteDensity
        allowing initialSurfaceCoverages to determine the number of moles of surface species.
        The number of moles of gas phase species is taken from initialGasMoleFractions.
        
        The coreSpeciesConcentrations array is then determined, in mol/m3 for gas phase
        and mol/m2 for surface species.

        The initial number of moles of a species j in the reactor is computed and stored in the
        y0 instance attribute.

        """
        ### 
        ### WARNING -- When updating this method, please be sure 
        ###            to also update the log_initial_conditions above
        ###            which unfortunately is maintained separately.
        ### 
        cdef double V, P, surface_volume_ratio_si

        ReactionSystem.set_initial_conditions(self)
        # self.y0 tracks number of moles of each core species

        # add only the gas phase species first
        self.y0 *= 0.
        for spec, moleFrac in self.initialGasMoleFractions.items():
            i = self.get_species_index(spec)
            self.y0[i] = moleFrac  # moles in reactor

        # Use ideal gas law to compute reactor volume
        V = constants.R * self.T.value_si * np.sum(self.y0[:self.numCoreSpecies]) / self.initialP.value_si
        self.V = V  # volume in m^3  (assume reactor volume is gas phase volume, i.e catalyst takes no space)

        surface_volume_ratio_si = self.surfaceVolumeRatio.value_si  # 1/m
        total_surface_sites = V * surface_volume_ratio_si * self.surfaceSiteDensity.value_si  # total surface sites in reactor

        for spec, coverage in self.initialSurfaceCoverages.items():
            i = self.get_species_index(spec)
            self.y0[i] = total_surface_sites * coverage  # moles in reactor

        for j, isSurfaceSpecies in enumerate(self.speciesOnSurface):  # should only go up to core species
            if isSurfaceSpecies:
                self.coreSpeciesConcentrations[j] = self.y0[j] / V / surface_volume_ratio_si  # moles per m2 of surface
            else:
                self.coreSpeciesConcentrations[j] = self.y0[j] / V  # moles per m3 of gas

    def compute_network_variables(self, pdepNetworks=None):
        # ToDo: this should allow pressure to vary?
        # for now, just call the base class version.
        ReactionSystem.compute_network_variables(self, pdepNetworks)

    def get_threshold_rate_constants(self, modelSettings):
        """
        Get the threshold rate constants for reaction filtering.
        """
        raise NotImplementedError("filterReactions=True for SurfaceReactor")
        # Set the maximum unimolecular rate to be kB*T/h
        unimolecular_threshold_rate_constant = 2.08366122e10 * self.T.value_si
        # Set the maximum bi/trimolecular rate by using the user-defined rate constant threshold
        bimolecular_threshold_rate_constant = modelSettings.filterThreshold
        # Maximum trimolecular rate constants are approximately three
        # orders of magnitude smaller (accounting for the unit
        # conversion from m^3/mol/s to m^6/mol^2/s) based on
        # extending the Smoluchowski equation to three molecules
        trimolecular_threshold_rate_constant = modelSettings.filterThreshold / 1e3
        return (unimolecular_threshold_rate_constant,
                bimolecular_threshold_rate_constant,
                trimolecular_threshold_rate_constant)

    @cython.boundscheck(False)
    def residual(self,
                 double t,
                 np.ndarray[np.float64_t, ndim=1] N,
                 np.ndarray[np.float64_t, ndim=1] dNdt,
                 np.ndarray[np.float64_t, ndim=1] senpar = np.zeros(1, np.float64)
                 ):

        """
        Return the residual function for the governing DAE system for the
        simple reaction system.
        """
        cdef np.ndarray[np.int_t, ndim=2] ir, ip, inet
        cdef np.ndarray[np.int_t, ndim=1] reactions_on_surface, species_on_surface
        cdef np.ndarray[np.float64_t, ndim=1] res, kf, kr, knet, delta, equilibrium_constants
        cdef int num_core_species, num_core_reactions, num_edge_species, num_edge_reactions, num_pdep_networks
        cdef int i, j, z, first, second, third
        cdef double k, V, reaction_rate, surface_volume_ratio_si
        cdef np.ndarray[np.float64_t, ndim=1] core_species_concentrations, core_species_rates, core_reaction_rates
        cdef np.ndarray[np.float64_t, ndim=1] edge_species_rates, edge_reaction_rates, network_leak_rates
        cdef np.ndarray[np.float64_t, ndim=1] C
        cdef np.ndarray[np.float64_t, ndim=2] jacobian, dgdk

        ir = self.reactantIndices
        ip = self.productIndices
        equilibrium_constants = self.Keq

        kf = self.kf  # are already 'per m3 of reactor' even for surface reactions
        kr = self.kb  # are already 'per m3 of reactor' even for surface reactions

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
        edge_species_rates = np.zeros_like(self.edgeSpeciesRates)
        edge_reaction_rates = np.zeros_like(self.edgeReactionRates)
        network_leak_rates = np.zeros_like(self.networkLeakRates)

        reactions_on_surface = self.reactionsOnSurface
        species_on_surface = self.speciesOnSurface
        surface_volume_ratio_si = self.surfaceVolumeRatio.value_si

        C = np.zeros_like(self.coreSpeciesConcentrations)
        V = self.V  # constant volume reactor

        for j in range(num_core_species):
            if species_on_surface[j]:
                C[j] = (N[j] / V) / surface_volume_ratio_si
            else:
                C[j] = N[j] / V
            #: surface species are in mol/m2, gas phase are in mol/m3
            core_species_concentrations[j] = C[j]

        for j in range(ir.shape[0]):
            k = kf[j]
            if ir[j, 0] >= num_core_species or ir[j, 1] >= num_core_species or ir[j, 2] >= num_core_species:
                reaction_rate = 0.0
            elif ir[j, 1] == -1:  # only one reactant
                reaction_rate = k * C[ir[j, 0]]
            elif ir[j, 2] == -1:  # only two reactants
                reaction_rate = k * C[ir[j, 0]] * C[ir[j, 1]]
            else:  # three reactants!! (really?)
                reaction_rate = k * C[ir[j, 0]] * C[ir[j, 1]] * C[ir[j, 2]]
            k = kr[j]
            if ip[j, 0] >= num_core_species or ip[j, 1] >= num_core_species or ip[j, 2] >= num_core_species:
                pass
            elif ip[j, 1] == -1:  # only one reactant
                reaction_rate -= k * C[ip[j, 0]]
            elif ip[j, 2] == -1:  # only two reactants
                reaction_rate -= k * C[ip[j, 0]] * C[ip[j, 1]]
            else:  # three reactants!! (really?)
                reaction_rate -= k * C[ip[j, 0]] * C[ip[j, 1]] * C[ip[j, 2]]

            "reaction_rate is now in mol/m3/s"
            # Set the reaction and species rates
            if j < num_core_reactions:
                # The reaction is a core reaction
                core_reaction_rates[j] = reaction_rate

                # Add/subtract the total reaction rate from each species rate
                # Since it's a core reaction we know that all of its reactants
                # and products are core species
                first = ir[j, 0]
                core_species_rates[first] -= reaction_rate
                second = ir[j, 1]
                if second != -1:
                    core_species_rates[second] -= reaction_rate
                    third = ir[j, 2]
                    if third != -1:
                        core_species_rates[third] -= reaction_rate
                first = ip[j, 0]
                core_species_rates[first] += reaction_rate
                second = ip[j, 1]
                if second != -1:
                    core_species_rates[second] += reaction_rate
                    third = ip[j, 2]
                    if third != -1:
                        core_species_rates[third] += reaction_rate

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
            else:  # three reactants!! (really?)
                reaction_rate = k * C[inet[j, 0]] * C[inet[j, 1]] * C[inet[j, 2]]
            network_leak_rates[j] = reaction_rate

        self.coreSpeciesConcentrations = core_species_concentrations
        self.coreSpeciesRates = core_species_rates
        self.coreReactionRates = core_reaction_rates
        self.edgeSpeciesRates = edge_species_rates
        self.edgeReactionRates = edge_reaction_rates
        self.networkLeakRates = network_leak_rates

        res = core_species_rates * V
        # mol/s

        if self.sensitivity and False:
            delta = np.zeros(len(N), np.float64)
            delta[:num_core_species] = res
            if self.jacobianMatrix is None:
                jacobian = self.jacobian(t, N, dNdt, 0, senpar)
            else:
                jacobian = self.jacobianMatrix
            dgdk = ReactionSystem.computeRateDerivative(self)
            for j in range(num_core_reactions + num_core_species):
                for i in range(num_core_species):
                    for z in range(num_core_species):
                        delta[(j + 1) * num_core_species + i] += jacobian[i, z] * N[(j + 1) * num_core_species + z]
                    delta[(j + 1) * num_core_species + i] += dgdk[i, j]
        else:
            delta = res
        delta = delta - dNdt

        # Return DELTA, IRES.  IRES is set to 1 in order to tell DASPK to evaluate the sensitivity residuals
        return delta, 1
