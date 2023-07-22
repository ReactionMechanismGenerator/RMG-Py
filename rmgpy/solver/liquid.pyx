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
    cdef public bint constant_volume
    cdef public double viscosity
    cdef public list const_spc_names
    cdef public list const_spc_indices
    cdef public dict initial_concentrations
    cdef public list Trange
    cdef public int n_sims
    cdef public dict sens_conditions

    def __init__(self, T, initial_concentrations, n_sims=1, termination=None, sensitive_species=None,
                 sensitivity_threshold=1e-3, sens_conditions=None, const_spc_names=None):

        ReactionSystem.__init__(self, termination, sensitive_species, sensitivity_threshold)

        if type(T) != list:
            self.T = Quantity(T)
        else:
            self.Trange = [Quantity(t) for t in T]

        self.P = Quantity(100000., 'kPa')  # Arbitrary high pressure (1000 Bar) to get reactions in the high-pressure limit!
        self.initial_concentrations = initial_concentrations  # should be passed in SI
        self.V = 0  # will be set from initial_concentrations in initialize_model
        self.constant_volume = True
        self.viscosity = 0  # in Pa*s

        #Constant concentration attributes
        self.const_spc_indices = None
        self.const_spc_names = const_spc_names  #store index of constant species
        self.sens_conditions = sens_conditions
        self.n_sims = n_sims

    def convert_initial_keys_to_species_objects(self, species_dict):
        """
        Convert the initial_concentrations dictionary from species names into species objects,
        using the given dictionary of species.
        """
        initial_concentrations = {}
        for label, moleFrac in self.initial_concentrations.items():
            if label == 'T':
                continue
            initial_concentrations[species_dict[label]] = moleFrac
        self.initial_concentrations = initial_concentrations

        conditions = {}
        if self.sens_conditions is not None:
            for label, value in self.sens_conditions.items():
                if label == 'T':
                    conditions[label] = value
                else:
                    conditions[species_dict[label]] = value
        self.sens_conditions = conditions

    def get_const_spc_indices(self, core_species):
        """Allow to identify constant Species position in solver"""
        for spc in self.const_spc_names:
            if self.const_spc_indices is None:  # initialize once the list if constant SPC declared
                self.const_spc_indices = []
            for item in core_species:
                # Need to identify the species object corresponding to the the string written in the input file
                if item.label == spc:
                    self.const_spc_indices.append(core_species.index(item))

    cpdef initialize_model(self, list core_species, list core_reactions, list edge_species, list edge_reactions,
                          list surface_species=None, list surface_reactions=None, list pdep_networks=None,
                          atol=1e-16, rtol=1e-8, sensitivity=False, sens_atol=1e-6, sens_rtol=1e-4,
                          filter_reactions=False, dict conditions=None):
        """
        Initialize a simulation of the liquid reactor using the provided kinetic
        model.
        """
        if surface_species is None:
            surface_species = []
        if surface_reactions is None:
            surface_reactions = []

        # First call the base class version of the method
        # This initializes the attributes declared in the base class
        ReactionSystem.initialize_model(self, core_species, core_reactions, edge_species, edge_reactions,
                                       surface_species=surface_species, surface_reactions=surface_reactions,
                                       pdep_networks=pdep_networks, atol=atol, rtol=rtol,
                                       sensitivity=sensitivity, sens_atol=sens_atol, sens_rtol=sens_rtol,
                                       filter_reactions=filter_reactions, conditions=conditions)

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
        Populates the forwardRateCoefficients, reverseRateCoefficients and equilibriumConstants
        arrays with the values computed at the temperature and (effective) pressure of the 
        reacion system.
        """

        for rxn in itertools.chain(core_reactions, edge_reactions):
            j = self.reaction_index[rxn]
            self.kf[j] = rxn.get_rate_coefficient(self.T.value_si, self.P.value_si)
            if rxn.reversible:
                self.Keq[j] = rxn.get_equilibrium_constant(self.T.value_si)
                self.kb[j] = self.kf[j] / self.Keq[j]

    def get_threshold_rate_constants(self, model_settings):
        """
        Get the threshold rate constants for reaction filtering.

        model_settings is not used here, but is needed so that the method
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

        The core_species_concentrations array is set to the values stored in the
        initial concentrations dictionary.

        The initial number of moles of a species j is computed and stored in the
        y0 instance attribute.

        """
        ReactionSystem.set_initial_conditions(self)

        for spec, conc in self.initial_concentrations.items():
            i = self.get_species_index(spec)
            self.core_species_concentrations[i] = conc

        V = 1.0 / np.sum(self.core_species_concentrations)
        self.V = V

        for j in range(self.num_core_species):
            self.y0[j] = self.core_species_concentrations[j] * V

    @cython.boundscheck(False)
    def residual(self, double t, np.ndarray[float_t, ndim=1] y, np.ndarray[float_t, ndim=1] dydt,
                 np.ndarray[float_t, ndim=1] senpar = np.zeros(1, float)):

        """
        Return the residual function for the governing DAE system for the
        liquid reaction system.
        """
        cdef np.ndarray[np.int_t, ndim=2] ir, ip, inet
        cdef np.ndarray[float_t, ndim=1] res, kf, kr, knet, delta, equilibrium_constants
        cdef Py_ssize_t num_core_species, num_core_reactions, num_edge_species, num_edge_reactions, num_pdep_networks
        cdef Py_ssize_t i, j, z, first, second, third
        cdef double k, V, reaction_rate
        cdef np.ndarray[float_t,ndim=1] core_species_concentrations, core_species_rates, core_reaction_rates
        cdef np.ndarray[float_t,ndim=1] edge_species_rates, edge_reaction_rates, network_leak_rates
        cdef np.ndarray[float_t,ndim=1] core_species_consumption_rates, core_species_production_rates
        cdef np.ndarray[float_t, ndim=1] C
        cdef np.ndarray[float_t, ndim=2] jacobian, dgdk

        ir = self.reactant_indices
        ip = self.product_indices
        equilibrium_constants = self.Keq

        kf = self.kf
        kr = self.kb

        inet = self.network_indices
        knet = self.network_leak_coefficients

        num_core_species = len(self.core_species_rates)
        num_core_reactions = len(self.core_reaction_rates)
        num_edge_species = len(self.edge_species_rates)
        num_edge_reactions = len(self.edge_reaction_rates)
        num_pdep_networks = len(self.network_leak_rates)

        res = np.zeros(num_core_species, float)

        core_species_concentrations = np.zeros_like(self.core_species_concentrations)
        core_species_rates = np.zeros_like(self.core_species_rates)
        core_reaction_rates = np.zeros_like(self.core_reaction_rates)
        core_species_consumption_rates = np.zeros_like(self.core_species_consumption_rates)
        core_species_production_rates = np.zeros_like(self.core_species_production_rates)
        edge_species_rates = np.zeros_like(self.edge_species_rates)
        edge_reaction_rates = np.zeros_like(self.edge_reaction_rates)
        network_leak_rates = np.zeros_like(self.network_leak_rates)

        C = np.zeros_like(self.core_species_concentrations)
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

        # chatelak: Same as in Java, core species rate = 0 if declared as constant
        if self.const_spc_indices is not None:
            for spc_index in self.const_spc_indices:
                core_species_rates[spc_index] = 0

        self.core_species_concentrations = core_species_concentrations
        self.core_species_rates = core_species_rates
        self.core_reaction_rates = core_reaction_rates
        self.core_species_production_rates = core_species_production_rates
        self.core_species_consumption_rates = core_species_consumption_rates
        self.edge_species_rates = edge_species_rates
        self.edge_reaction_rates = edge_reaction_rates
        self.network_leak_rates = network_leak_rates

        res = core_species_rates * V

        if self.sensitivity:
            delta = np.zeros(len(y), float)
            delta[:num_core_species] = res
            if self.jacobian_matrix is None:
                jacobian = self.jacobian(t, y, dydt, 0, senpar)
            else:
                jacobian = self.jacobian_matrix
            dgdk = ReactionSystem.compute_rate_derivative(self)
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
    def jacobian(self, double t, np.ndarray[float_t, ndim=1] y, np.ndarray[float_t, ndim=1] dydt,
                 double cj, np.ndarray[float_t, ndim=1] senpar = np.zeros(1, float)):
        """
        Return the analytical Jacobian for the reaction system.
        """
        cdef np.ndarray[np.int_t, ndim=2] ir, ip
        cdef np.ndarray[float_t, ndim=1] kf, kr, C
        cdef np.ndarray[float_t, ndim=2] pd
        cdef Py_ssize_t num_core_reactions, num_core_species, i, j
        cdef double k, V, Ctot, deriv, corr

        ir = self.reactant_indices
        ip = self.product_indices

        kf = self.kf
        kr = self.kb
        num_core_reactions = len(self.core_reaction_rates)
        num_core_species = len(self.core_species_concentrations)

        pd = -cj * np.identity(num_core_species, float)

        V = self.V  # volume is constant

        C = np.zeros_like(self.core_species_concentrations)
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

        self.jacobian_matrix = pd + cj * np.identity(num_core_species, float)
        return pd
