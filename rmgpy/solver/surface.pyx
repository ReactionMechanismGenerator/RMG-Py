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
    cdef public ScalarQuantity P_initial
    cdef public double V
    cdef public bint constant_volume

    cdef public list Trange
    cdef public list Prange
    cdef public int n_sims
    cdef public dict sens_conditions

    cdef public dict initial_gas_mole_fractions
    cdef public dict initial_surface_coverages
    cdef public ScalarQuantity surface_volume_ratio
    cdef public ScalarQuantity surface_site_density
    cdef public np.ndarray reactions_on_surface  # (catalyst surface, not core/edge surface)
    cdef public np.ndarray species_on_surface  # (catalyst surface, not core/edge surface)

    def __init__(self,
                 T,
                 P_initial,
                 initial_gas_mole_fractions,
                 initial_surface_coverages,
                 surface_volume_ratio,
                 surface_site_density,
                 n_sims=None,
                 termination=None,
                 sensitive_species=None,
                 sensitivity_threshold=1e-3,
                 sens_conditions=None,
                 ):
        ReactionSystem.__init__(self,
                                termination,
                                sensitive_species,
                                sensitivity_threshold)

        if isinstance(T, list):
            self.Trange = [Quantity(t) for t in T]
        else:
            self.T = Quantity(T)
        if isinstance(P_initial, list):
            raise NotImplementedError("Can't do ranges of initial pressures for surface reactors yet")
        else:
            self.P_initial = Quantity(P_initial)
        self.initial_gas_mole_fractions = initial_gas_mole_fractions
        self.initial_surface_coverages = initial_surface_coverages
        self.surface_volume_ratio = Quantity(surface_volume_ratio)
        self.surface_site_density = Quantity(surface_site_density)
        self.V = 0  # will be set from ideal gas law in initialize_model
        self.constant_volume = True
        self.sens_conditions = sens_conditions
        self.n_sims = n_sims

    def convert_initial_keys_to_species_objects(self, species_dict):
        """
        Convert the initial_gas_mole_fractions and initial_surface_coverages dictionaries
        from species names into species objects,
        using the given dictionary of species.
        """
        initial_gas_mole_fractions = {}
        for label, moleFrac in self.initial_gas_mole_fractions.items():
            initial_gas_mole_fractions[species_dict[label]] = moleFrac
        self.initial_gas_mole_fractions = initial_gas_mole_fractions
        initial_surface_coverages = {}
        for label, surfaceCoverage in self.initial_surface_coverages.items():
            initial_surface_coverages[species_dict[label]] = surfaceCoverage
        self.initial_surface_coverages = initial_surface_coverages

    cpdef initialize_model(self,
                          list core_species,
                          list core_reactions,
                          list edge_species,
                          list edge_reactions,
                          list surface_species=[],
                          list surface_reactions=[],
                          list pdep_networks=None,
                          atol=1e-16,
                          rtol=1e-8,
                          sensitivity=False,
                          sens_atol=1e-6,
                          sens_rtol=1e-4,
                          filter_reactions=False,
                          dict conditions=None,
                          ):
        """
        Initialize a simulation of the simple reactor using the provided kinetic
        model.
        """

        # First call the base class version of the method
        # This initializes the attributes declared in the base class
        ReactionSystem.initialize_model(self,
                                       core_species=core_species,
                                       core_reactions=core_reactions,
                                       edge_species=edge_species,
                                       edge_reactions=edge_reactions,
                                       surface_species=surface_species,
                                       surface_reactions=surface_reactions,
                                       pdep_networks=pdep_networks,
                                       atol=atol,
                                       rtol=rtol,
                                       sensitivity=sensitivity,
                                       sens_atol=sens_atol,
                                       sens_rtol=sens_rtol,
                                       filter_reactions=filter_reactions,
                                       conditions=conditions,
                                       )
        cdef np.ndarray[np.int_t, ndim=1] species_on_surface, reactions_on_surface
        cdef int index
        #: 1 if it's on a surface, 0 if it's in the gas phase
        reactions_on_surface = np.zeros((self.num_core_reactions + self.num_edge_reactions), np.int)
        species_on_surface = np.zeros((self.num_core_species), np.int)
        for spec, index in self.species_index.items():
            if index >= self.num_core_species:
                continue
            if spec.contains_surface_site():
                species_on_surface[index] = 1
        for rxn, index in self.reaction_index.items():
            if rxn.is_surface_reaction():
                reactions_on_surface[index] = 1
        self.species_on_surface = species_on_surface
        self.reactions_on_surface = reactions_on_surface

        # Set initial conditions
        self.set_initial_conditions()

        # Compute reaction thresholds if reaction filtering is turned on
        if filter_reactions:
            ReactionSystem.set_initial_reaction_thresholds(self)

        # Generate forward and reverse rate coefficients k(T,P)
        self.generate_rate_coefficients(core_reactions, edge_reactions)

        ReactionSystem.compute_network_variables(self, pdep_networks)

        ReactionSystem.set_initial_derivative(self)

        # Initialize the model
        ReactionSystem.initialize_solver(self)

    def generate_rate_coefficients(self, core_reactions, edge_reactions):
        """
        Populates the kf, kb and equilibriumConstants
        arrays with the values computed at the temperature and (effective) pressure of the 
        reaction system.
        """

        cdef double P, surface_volume_ratio_si

        surface_volume_ratio_si = self.surface_volume_ratio.value_si
        # ToDo: Pressure should come from ideal gas law?
        P = self.P_initial.value_si

        warned = False
        for rxn in itertools.chain(core_reactions, edge_reactions):
            j = self.reaction_index[rxn]

            # ToDo: get_rate_coefficient should also depend on surface coverages vector

            if rxn.is_surface_reaction():
                """
                Be careful! From here on kf and kb will now be in Volume units,
                even for surface reactions (which you may expect to be in Area units).
                This is to avoid repeatedly multiplying a bunch of things inside every 
                loop of the ODE solver.
                """
                self.kf[j] = (surface_volume_ratio_si *
                              rxn.get_surface_rate_coefficient(self.T.value_si,
                                                               self.surface_site_density.value_si
                                                               ))
            else:
                if not warned and rxn.kinetics.is_pressure_dependent():
                    logging.warning("Pressure may be varying, but using initial pressure to evaluate k(T,P) expressions!")
                    warned = True
                self.kf[j] = rxn.get_rate_coefficient(self.T.value_si, P)
            if rxn.reversible:
                # ToDo: get_equilibrium_constant should be coverage dependent
                self.Keq[j] = rxn.get_equilibrium_constant(self.T.value_si)
                self.kb[j] = self.kf[j] / self.Keq[j]

    def log_initial_conditions(self, number=None):
        """
        Log to the console some information about this reaction system.
        
        Should correspond to the calculations done in set_initial_conditions.
        """
        logging.info("\nSurface reaction system {}".format(number if number is not None else ""))
        logging.info("Gas phase mole fractions:")
        total_gas_moles = 0
        for spec, moleFrac in self.initial_gas_mole_fractions.items():
            logging.info("  {0:20s} {1:.5g}".format(spec, moleFrac))
            total_gas_moles += moleFrac
        logging.info("Total gas phase:          {:.3g} moles".format(total_gas_moles))
        logging.info("Pressure:                 {:.3g} Pa".format(self.P_initial.value_si))
        logging.info("Temperature:              {} K".format(self.T.value_si))
        V = constants.R * self.T.value_si * total_gas_moles / self.P_initial.value_si
        logging.info("Reactor volume:           {:.3g} m3".format(V))
        surface_volume_ratio_si = self.surface_volume_ratio.value_si  # 1/m
        logging.info("Surface/volume ratio:     {:.3g} m2/m3".format(surface_volume_ratio_si))
        logging.info("Surface site density:     {:.3g} mol/m2".format(self.surface_site_density.value_si))
        total_surface_sites = V * surface_volume_ratio_si * self.surface_site_density.value_si  # total surface sites in reactor
        logging.info("Surface sites in reactor: {:.3g} moles".format(total_surface_sites))
        logging.info("Initial surface coverages (and amounts):")
        for spec, coverage in self.initial_surface_coverages.items():
            logging.info("  {:18s} {:.5g} = {:.5g} moles".format(spec, coverage, total_surface_sites * coverage))

    def set_initial_conditions(self):
        """
        Sets the initial conditions of the rate equations that represent the 
        current reactor model.

        The volume is set to the value in m3 required to contain 
        one mole total of gas phase core species at start.

        The total surface sites are calculated from surface_volume_ratio and surface_site_density
        allowing initial_surface_coverages to determine the number of moles of surface species.
        The number of moles of gas phase species is taken from initial_gas_mole_fractions.
        
        The core_species_concentrations array is then determined, in mol/m3 for gas phase
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
        for spec, moleFrac in self.initial_gas_mole_fractions.items():
            i = self.get_species_index(spec)
            self.y0[i] = moleFrac  # moles in reactor

        # Use ideal gas law to compute reactor volume
        V = constants.R * self.T.value_si * np.sum(self.y0[:self.num_core_species]) / self.P_initial.value_si
        self.V = V  # volume in m^3  (assume reactor volume is gas phase volume, i.e catalyst takes no space)

        surface_volume_ratio_si = self.surface_volume_ratio.value_si  # 1/m
        total_surface_sites = V * surface_volume_ratio_si * self.surface_site_density.value_si  # total surface sites in reactor

        for spec, coverage in self.initial_surface_coverages.items():
            i = self.get_species_index(spec)
            self.y0[i] = total_surface_sites * coverage  # moles in reactor

        for j, isSurfaceSpecies in enumerate(self.species_on_surface):  # should only go up to core species
            if isSurfaceSpecies:
                self.core_species_concentrations[j] = self.y0[j] / V / surface_volume_ratio_si  # moles per m2 of surface
            else:
                self.core_species_concentrations[j] = self.y0[j] / V  # moles per m3 of gas

    def compute_network_variables(self, pdep_networks=None):
        # ToDo: this should allow pressure to vary?
        # for now, just call the base class version.
        ReactionSystem.compute_network_variables(self, pdep_networks)

    def get_threshold_rate_constants(self, model_settings):
        """
        Get the threshold rate constants for reaction filtering.
        """
        raise NotImplementedError("filter_reactions=True for SurfaceReactor")
        # Set the maximum unimolecular rate to be kB*T/h
        unimolecular_threshold_rate_constant = 2.08366122e10 * self.T.value_si
        # Set the maximum bi/trimolecular rate by using the user-defined rate constant threshold
        bimolecular_threshold_rate_constant = model_settings.filter_threshold
        # Maximum trimolecular rate constants are approximately three
        # orders of magnitude smaller (accounting for the unit
        # conversion from m^3/mol/s to m^6/mol^2/s) based on
        # extending the Smoluchowski equation to three molecules
        trimolecular_threshold_rate_constant = model_settings.filter_threshold / 1e3
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

        ir = self.reactant_indices
        ip = self.product_indices
        equilibrium_constants = self.Keq

        kf = self.kf  # are already 'per m3 of reactor' even for surface reactions
        kr = self.kb  # are already 'per m3 of reactor' even for surface reactions

        inet = self.network_indices
        knet = self.network_leak_coefficients

        num_core_species = len(self.core_species_rates)
        num_core_reactions = len(self.core_reaction_rates)
        num_edge_species = len(self.edge_species_rates)
        num_edge_reactions = len(self.edge_reaction_rates)
        num_pdep_networks = len(self.network_leak_rates)

        res = np.zeros(num_core_species, np.float64)

        core_species_concentrations = np.zeros_like(self.core_species_concentrations)
        core_species_rates = np.zeros_like(self.core_species_rates)
        core_reaction_rates = np.zeros_like(self.core_reaction_rates)
        edge_species_rates = np.zeros_like(self.edge_species_rates)
        edge_reaction_rates = np.zeros_like(self.edge_reaction_rates)
        network_leak_rates = np.zeros_like(self.network_leak_rates)

        reactions_on_surface = self.reactions_on_surface
        species_on_surface = self.species_on_surface
        surface_volume_ratio_si = self.surface_volume_ratio.value_si

        C = np.zeros_like(self.core_species_concentrations)
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
            if inet[j, 0] != -1: #all source species are in the core
                k = knet[j]
                if inet[j, 1] == -1:  # only one reactant
                    reaction_rate = k * C[inet[j, 0]]
                elif inet[j, 2] == -1:  # only two reactants
                    reaction_rate = k * C[inet[j, 0]] * C[inet[j, 1]]
                else:  # three reactants
                    reaction_rate = k * C[inet[j, 0]] * C[inet[j, 1]] * C[inet[j, 2]]
                    network_leak_rates[j] = reaction_rate
            else:
                network_leak_rates[j] = 0.0

        self.core_species_concentrations = core_species_concentrations
        self.core_species_rates = core_species_rates
        self.core_reaction_rates = core_reaction_rates
        self.edge_species_rates = edge_species_rates
        self.edge_reaction_rates = edge_reaction_rates
        self.network_leak_rates = network_leak_rates

        res = core_species_rates * V
        # mol/s

        if self.sensitivity and False:
            delta = np.zeros(len(N), np.float64)
            delta[:num_core_species] = res
            if self.jacobian_matrix is None:
                jacobian = self.jacobian(t, N, dNdt, 0, senpar)
            else:
                jacobian = self.jacobian_matrix
            dgdk = ReactionSystem.compute_rate_derivative(self)
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
