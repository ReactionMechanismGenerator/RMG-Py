#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2023 Prof. William H. Green (whgreen@mit.edu),           #
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

import itertools
import logging
import numpy as np
from typing import Dict, List, Optional, Union

import rmgpy.constants as constants
from rmgpy.quantity import Quantity
from rmgpy.solver.base import ReactionSystem, TerminationConversion, TerminationRateRatio, TerminationTime
from rmgpy.solver.polymer import HybridPolymerSystem, MassTransferConfig, PolymerPoolConfig
from rmgpy.species import Species


def _base_label(label):
    """Polymer base label = everything before the first '(' (strips the
    proxy multiplicity / index suffix, e.g. 'PS(2)' -> 'PS')."""
    return label.partition("(")[0]


class HybridPolymerReactor(ReactionSystem):
    """
    A biphasic reactor input specification for polymer pyrolysis and degradation simulations.

    This reactor models two distinct phases: a gas phase and a polymer melt phase.
    It couples discrete chemical species (gas and explicit oligomers) with a statistical
    method-of-moments representation for long polymer chains.

    Args:
        temperature (Quantity): The initial temperature of the reactor (e.g. '300 K').
        pressure (Quantity): The initial pressure of the reactor (e.g. '1 bar').
        initialMoles (dict): A dictionary {Species: float} representing the initial composition
                             of the GAS phase. Values are interpreted as MOLES, not mole fractions.
        polymerPhase (PolymerPhase): Configuration object containing all properties of the polymer
                                     melt phase (density, initial moments, pools, mass transfer).
        terminationConversion (dict, optional): A dictionary {Species: float} or {str: float}
                                                specifying the fractional conversion at which to
                                                terminate the simulation (0.0 to 1.0).
        terminationTime (Quantity, optional): The maximum time to simulate (e.g. '10 s').
        terminationRateRatio (float, optional): The minimum ratio of production rate to consumption rate
                                               for all core species before terminating the simulation.
        sensitivity (list, optional): A list of Species objects or reaction labels to calculate
                                      sensitivities for.
        sensitivityThreshold (float, optional): The cutoff threshold for sensitivity analysis.
                                                Default is 1e-3.
        sens_conditions (dict, optional): A dictionary specifying conditions for sensitivity
                                            analysis (e.g. {'T': 800, 'P': 1e5}).
        constant_gas_volume (bool, optional): If True, the gas phase volume remains fixed at its
                                              initial value calculated from T, P, and initial gas moles.
                                              If False (default), the gas volume expands/contracts
                                              isobarically to maintain constant Pressure.
        allow_unpaired_reference_state (bool, optional): If True, bypass the build-time refusal on
                                              reversible reactions whose thermo reference-state term
                                              is unpaired (U > 3 decades); the census is still
                                              logged. The deck author asserts their thermo handles
                                              the melt reference state. See the invariant section
                                              of docs/multi_pool_design.md.
    """
    def __init__(self,
                 temperature,
                 pressure,
                 initialMoles: Dict[Species, float],
                 polymerPhase: 'PolymerPhase',
                 terminationConversion: Optional[Dict[Union[Species, str], float]] = None,
                 terminationTime: Optional[Union[Quantity, float]] = None,
                 terminationRateRatio=None,
                 sensitivity: Optional[List[Union[Species, str]]] = None,
                 sensitivityThreshold: float = 1e-3,
                 sens_conditions: Optional[Dict[str, Union[float, Quantity]]] = None,
                 constant_gas_volume: bool = False,
                 allow_unpaired_reference_state: bool = False,
                 n_sims=1,
                 ):
        ReactionSystem.__init__(self)

        if not hasattr(self, 'listeners'):
            self.listeners = []

        if type(temperature) != list:
            self.T = Quantity(temperature)
            self.Trange = [self.T, self.T]
            self.temperature = Quantity(temperature)
        else:
            self.Trange = [Quantity(t) for t in temperature]
            self.temperature = [Quantity(t) for t in temperature]

        if type(pressure) != list:
            self.P = Quantity(pressure)
            self.Prange = [self.P, self.P]
            self.pressure = Quantity(pressure)
        else:
            self.Prange = [Quantity(p) for p in pressure]
            self.pressure = [Quantity(p) for p in pressure]

        self.initialMoles = initialMoles
        self.initial_mole_fractions = initialMoles
        self.polymerPhase = polymerPhase
        self.terminationConversion = terminationConversion
        self.terminationTime = terminationTime
        self.terminationRateRatio = terminationRateRatio
        self.sensitivity = sensitivity
        self.sensitive_species = list()
        self.sensitivityThreshold = sensitivityThreshold
        self.sens_conditions = sens_conditions
        self.constant_gas_volume = constant_gas_volume
        self.allow_unpaired_reference_state = allow_unpaired_reference_state
        self.n_sims = n_sims
        self.const_spc_names = []
        self.solver = None

    def initialize_model(self, core_species, core_reactions, edge_species, edge_reactions,
                          surface_species=None, surface_reactions=None, pdep_networks=None,
                          atol=1e-16, rtol=1e-8, sensitivity=False, sens_atol=1e-6, sens_rtol=1e-4,
                          filter_reactions=False, conditions=None, **kwargs):
        """
        Standard RMG hook:
        1. Update conditions (T, P) on the blueprint.
        2. Create the fast Cython solver.
        3. Delegate numerical setup to the solver.
        """
        # 1. Update the Blueprint (Reactor) settings based on RMG 'conditions'
        ReactionSystem.initialize_model(self, core_species=core_species, core_reactions=core_reactions,
                                        edge_species=edge_species, edge_reactions=edge_reactions,
                                        surface_species=surface_species, surface_reactions=surface_reactions,
                                        pdep_networks=pdep_networks, atol=atol, rtol=rtol, sensitivity=sensitivity,
                                        sens_atol=sens_atol, sens_rtol=sens_rtol, filter_reactions=filter_reactions,
                                        conditions=conditions)

        # 2. Build the Solver Engine (using the updated T/P from step 1)
        self.solver = self.to_solver_object(core_species, core_reactions, edge_species, edge_reactions)
        self.solver.listeners = self.listeners

        # 3. Delegate ALL numerical initialization to the Solver Engine
        result = self.solver.initialize_model(core_species=core_species, core_reactions=core_reactions,
                                            edge_species=edge_species, edge_reactions=edge_reactions,
                                            surface_species=surface_species, surface_reactions=surface_reactions,
                                            pdep_networks=pdep_networks, atol=atol, rtol=rtol, sensitivity=sensitivity,
                                            sens_atol=sens_atol, sens_rtol=sens_rtol, filter_reactions=filter_reactions,
                                            conditions=conditions)

        # 4. Sync threshold arrays back from solver so RMG main can read them
        self.unimolecular_threshold = self.solver.unimolecular_threshold
        self.bimolecular_threshold = self.solver.bimolecular_threshold
        if hasattr(self.solver, 'trimolecular_threshold') and self.solver.trimolecular_threshold is not None:
            self.trimolecular_threshold = self.solver.trimolecular_threshold

        return result

    def simulate(self, core_species, core_reactions, edge_species, edge_reactions, **kwargs):
        """
        Run the simulation using the underlying solver.
        pass the command to the numerical solver created in initialize()
        """
        n_core = len(core_species)
        n_rxn = len(core_reactions)
        if (self.solver is None or
                self.solver.num_core_species != n_core or
                self.solver.num_core_reactions != n_rxn):
            self.initialize_model(core_species, core_reactions, edge_species, edge_reactions, **kwargs)
        result = self.solver.simulate(core_species, core_reactions, edge_species, edge_reactions, **kwargs)

        # Sync threshold arrays back from solver so RMG main can read them
        self.unimolecular_threshold = self.solver.unimolecular_threshold
        self.bimolecular_threshold = self.solver.bimolecular_threshold
        if hasattr(self.solver, 'trimolecular_threshold') and self.solver.trimolecular_threshold is not None:
            self.trimolecular_threshold = self.solver.trimolecular_threshold

        return result

    def convert_initial_keys_to_species_objects(self, species_dict):
        """
        Convert species labels into species objects across the reactor
        and the associated polymer phase.
        """
        # 1. Convert Reactor-level initialMoles
        # Matches your signature: initialMoles = {Species/str: float}
        new_initial_moles = {}
        for key, value in self.initialMoles.items():
            if isinstance(key, str):
                new_initial_moles[species_dict[key]] = value
            else:
                new_initial_moles[key] = value
        self.initialMoles = new_initial_moles

        # 2. Convert Phase-level initial_explicit
        # The polymerPhase was initialized with labels; swap them for objects now.
        if hasattr(self.polymerPhase, 'initial_explicit'):
            new_explicit = {}
            for key, value in self.polymerPhase.initial_explicit.items():
                if isinstance(key, str):
                    new_explicit[species_dict[key]] = value
                else:
                    new_explicit[key] = value
            self.polymerPhase.initial_explicit = new_explicit

        # 3. Convert Pool-level species (Monomers and Mu_species)
        if hasattr(self.polymerPhase, 'pools'):
            for pool in self.polymerPhase.pools:
                # Resolve Monomer
                if isinstance(pool.monomer, str):
                    pool.monomer = species_dict[pool.monomer]

                # Resolve Mu_species (the 3 placeholder species for moments)
                if pool.mu_species:
                    new_mu = []
                    for spc in pool.mu_species:
                        if isinstance(spc, str):
                            new_mu.append(species_dict[spc])
                        else:
                            new_mu.append(spc)
                    pool.mu_species = new_mu

                # Resolve Explicit Map (Oligomers)
                if pool.explicit_map:
                    new_map = {}
                    for dp, spc in pool.explicit_map.items():
                        if isinstance(spc, str):
                            new_map[dp] = species_dict[spc]
                        else:
                            new_map[dp] = spc
                    pool.explicit_map = new_map

        # 4. Handle Sensitivity (if applicable)
        if self.sensitivity:
            new_sens = []
            for spec in self.sensitivity:
                if isinstance(spec, str):
                    new_sens.append(species_dict[spec])
                else:
                    new_sens.append(spec)
            self.sensitive_species = new_sens

    def to_solver_object(self, core_species, core_reactions, edge_species, edge_reactions):
        """
        Convert this Input settings object into a runnable Solver engine.
        """

        # 0. Create efficient lookup map (Performance: O(1) vs O(N))
        spc_map = {spc: i for i, spc in enumerate(core_species)}

        # Validate initialMoles keys
        unknown_initials = [spc for spc in self.initialMoles if spc not in spc_map]
        if unknown_initials:
            raise ValueError(f"Initial moles specified for species not in core: {unknown_initials}")

        # 1. Calculate Polymer Volume (Mass / Density)
        # Note: Validates that explicit species don't exceed total distribution mass.
        V_poly = self.polymerPhase.calculate_volume()

        # 2. Identify Phases (Robust Masking)
        # Determines which core species are Gas vs Polymer
        gas_mask = self.polymerPhase.get_gas_mask(core_species)
        if len(gas_mask) != len(core_species):
            # Emergency resize if get_gas_mask returns a shorter list
            new_mask = np.zeros(len(core_species), dtype=bool)
            new_mask[:len(gas_mask)] = gas_mask
            gas_mask = new_mask

        # Item 17 (spec 2026-06-12 SS3(a), A5-2) stage-1 prospective seed: the
        # SAME config-keyed classifier over chain(core, edge) -- same
        # function, same phase object, longer list; no second classifier
        # exists to drift. A5-2: this build-frozen seed went STALE on the
        # engine-reuse path (the engine outlives the edge list that produced
        # the seed), so it is now only a one-time first-build hint; the
        # AUTHORITATIVE production source is the classifier handle passed below
        # (prospective_classifier=self.polymerPhase.get_gas_mask), which the
        # solver re-runs over the LIVE chain(core, edge) on every
        # initialize_model. The solver's stage-2 pass + rider R1 (core-prefix
        # parity raise) + R1-EDGE (edge-suffix provenance raise) complete and
        # verify the construction at initialize_model time.
        prospective_seed = self.polymerPhase.get_gas_mask(
            list(core_species) + list(edge_species))

        # 3. Calculate Initial Gas Volume (Headspace)
        # Logic: Sum MOLES of species that are actually in the gas phase.
        total_gas_moles = 0.0
        for spc, moles in self.initialMoles.items():
            # We already validated spc is in spc_map
            idx = spc_map[spc]
            if gas_mask[idx]:
                total_gas_moles += moles

        V_gas0 = None
        if total_gas_moles > 0:
            V_gas0 = (total_gas_moles * constants.R * self.temperature.value_si) / self.pressure.value_si

        # Enforce consistency for constant volume constraint
        if self.constant_gas_volume and (V_gas0 is None or V_gas0 <= 0):
            raise ValueError("HybridPolymerReactor with constant_gas_volume=True requires positive initial gas moles "
                "to define the headspace volume.")

        # 4. Construct Termination Objects (Handle Species vs Label)
        termination = list()
        if self.terminationTime is not None:
            termination.append(TerminationTime(self.terminationTime))

        if self.terminationConversion:
            for spec_key, conv in self.terminationConversion.items():
                spec_obj = None

                # Case A: Input is a string label
                if isinstance(spec_key, str):
                    matches = [s for s in core_species if getattr(s, "label", None) == spec_key]
                    if len(matches) == 0:
                        raise ValueError(f"TerminationConversion label '{spec_key}' not found in core species.")
                    if len(matches) > 1:
                        raise ValueError(f"TerminationConversion label '{spec_key}' is ambiguous (matches multiple species).")
                    spec_obj = matches[0]

                # Case B: Input is a Species object
                else:
                    spec_obj = spec_key
                    if spec_obj not in spc_map:
                        raise ValueError(f"TerminationConversion species '{spec_obj}' is not in the core species list.")

                termination.append(TerminationConversion(spec_obj, conv))

        if self.terminationRateRatio is not None:
            termination.append(TerminationRateRatio(self.terminationRateRatio))

        # 5. Convert Input Objects -> Solver Configs
        # Pass spc_map to avoid re-searching the list
        pool_configs = [p.to_config(spc_map) for p in self.polymerPhase.pools]
        # Auto-register scission/spawn daughter pools (proxy_reaction_reality_rules.md
        # Layer 2): daughters are registered as core Polymer species + _muN dummies by
        # _register_polymer, but polymerPhase.pools never grows at runtime. Derive
        # their pool configs from the live core species so their stamped
        # SCISSION_FRAGMENT/MIGRATION flux resolves instead of demoting to UNRESOLVED.
        static_pool_labels = {p.label for p in self.polymerPhase.pools}
        pool_configs += derive_daughter_pool_configs(core_species, spc_map, static_pool_labels)
        mt_configs = [mt.to_config(spc_map) for mt in self.polymerPhase.mass_transfer]

        # Validate indices consistent with get_gas_mask()
        for mt_cfg in mt_configs:
            gi = mt_cfg.gas_index
            pi = mt_cfg.poly_index

            if not gas_mask[gi]:
                raise ValueError("MassTransfer error: gas_species mapped to non-gas by get_gas_mask().")
            if gas_mask[pi]:
                raise ValueError("MassTransfer error: poly_species mapped to gas by get_gas_mask().")

        # 5.5. Extract Pressure-Dependent Collider Data
        # We mirror SimpleReactor logic to ensure gas-phase P-dep reactions
        # use efficiency-weighted effective pressure.
        pdep_indices, pdep_kinetics, efficiencies = [], [], []

        for i, rxn in enumerate(core_reactions):
            if rxn.kinetics is not None and rxn.kinetics.is_pressure_dependent():
                if hasattr(rxn.kinetics, 'efficiencies') and rxn.kinetics.efficiencies:
                    pdep_indices.append(i)
                    pdep_kinetics.append(rxn.kinetics)
                    # Calculate efficiencies for all core species
                    efficiencies.append(rxn.kinetics.get_effective_collider_efficiencies(core_species))

        pdep_collision_indices = np.array(pdep_indices, int)
        collider_eff = np.array(efficiencies, float)

        poly_labels = set()
        for spc in self.polymerPhase.initial_explicit.keys():
            poly_labels.add(spc.label)
        for pool in self.polymerPhase.pools:
            poly_labels.add(pool.proxy_species.label)
            for mu in pool.mu_species:
                poly_labels.add(mu.label)
            if pool.explicit_map:
                for dp, spc in pool.explicit_map.items():
                    poly_labels.add(spc.label)
        for mt in self.polymerPhase.mass_transfer:
            poly_labels.add(mt.poly_species.label)

        # 6. Instantiate the Numerical Engine
        # Note: We pass 'initialMoles' to the solver's 'initial_mole_fractions' argument
        # to satisfy the base class signature, but the solver correctly interprets them as moles.
        solver = HybridPolymerSystem(
            T=self.temperature.value_si,
            P=self.pressure.value_si,
            initial_mole_fractions=self.initialMoles,  # Passed as moles
            V_poly=V_poly,
            polymer_pools=pool_configs,
            mass_transfer=mt_configs,
            polymer_species_labels=poly_labels,
            gas_species_mask=gas_mask,
            prospective_gas_mask=prospective_seed,
            # A5-2: the authoritative production source for stage 1 -- the
            # solver re-runs this over the LIVE chain(core, edge) on every
            # initialize_model (base.pyx:simulate already calls
            # initialize_model with the live edge), so the frozen seed above
            # can never go stale and silently demote the build to the
            # edge-defaults-GAS fallback. Production builds are held to the
            # stage-1 contract (allow_default_prospective_edge left default
            # False -> R1-EDGE raises if the fallback ever fires here).
            prospective_classifier=self.polymerPhase.get_gas_mask,
            # Spec 2026-06-29: the live edge-daughter condensed-mask predicate.
            # Bound method re-run over chain(core, edge) every initialize_model
            # (callable, never frozen -> no stale-set hazard).
            prospective_condensed_edge_daughter_classifier=self.polymerPhase.get_condensed_edge_daughter_bases,
            constant_gas_volume=self.constant_gas_volume,
            V_gas0=V_gas0,
            initial_polymer_moments=self.polymerPhase.initial_moments,
            initial_explicit_species=self.polymerPhase.initial_explicit,
            termination=termination,
            sensitive_species=self.sensitive_species,
            sensitivity_threshold=self.sensitivityThreshold,
            sens_conditions=self.sens_conditions,
            const_spc_names=self.const_spc_names,
            pdep_collision_reaction_indices=pdep_collision_indices,
            pdep_collider_kinetics=pdep_kinetics,
            collider_efficiencies=collider_eff,
            allow_unpaired_reference_state=self.allow_unpaired_reference_state,
        )

        solver.V = (V_gas0 if V_gas0 is not None else 0.0) + V_poly

        species_to_pool = np.full(len(core_species), -1, dtype=np.int32)

        for p_idx, pool in enumerate(self.polymerPhase.pools):
            proxy = getattr(pool, "proxy_species", None)
            if proxy is None:
                continue
            # Find the RMG Species object that acts as the proxy for this pool
            try:
                s_idx = core_species.index(proxy)
                species_to_pool[s_idx] = p_idx
            except ValueError:
                continue

        solver.species_to_pool_indices = species_to_pool

        return solver

    def set_colliders(self, core_reactions, edge_reactions, core_species):
        """
        Store collider efficiencies and reaction indices for pdep reactions that have collider efficiencies,
        and store specific collider indices
        """
        pdep_collider_reaction_indices = []
        self.pdep_collider_kinetics = []
        collider_efficiencies = []
        pdep_specific_collider_reaction_indices = []
        self.pdep_specific_collider_kinetics = []
        self.specific_collider_species = []

        for rxn in itertools.chain(core_reactions, edge_reactions):
            if rxn.kinetics.is_pressure_dependent():
                if rxn.kinetics.efficiencies:
                    j = self.reaction_index[rxn]
                    pdep_collider_reaction_indices.append(j)
                    self.pdep_collider_kinetics.append(rxn.kinetics)
                    collider_efficiencies.append(rxn.kinetics.get_effective_collider_efficiencies(core_species))
                if rxn.specific_collider:
                    pdep_specific_collider_reaction_indices.append(self.reaction_index[rxn])
                    self.pdep_specific_collider_kinetics.append(rxn.kinetics)
                    self.specific_collider_species.append(rxn.specific_collider)

        self.pdep_collision_reaction_indices = np.array(pdep_collider_reaction_indices, int)
        self.collider_efficiencies = np.array(collider_efficiencies, float)
        self.pdep_specific_collider_reaction_indices = np.array(pdep_specific_collider_reaction_indices, int)

    def generate_rate_coefficients(self, core_reactions, edge_reactions):
        """
        Populates the forward rate coefficients (kf), reverse rate coefficients (kb)
        and equilibrium constants (Keq) arrays with the values computed at the temperature
        and (effective) pressure of the reaction system.
        """
        for rxn in itertools.chain(core_reactions, edge_reactions):
            j = self.reaction_index[rxn]
            Peff = self.calculate_effective_pressure(rxn)
            self.kf[j] = rxn.get_rate_coefficient(self.T.value_si, Peff)

            if rxn.reversible:
                self.Keq[j] = rxn.get_equilibrium_constant(self.T.value_si)
                self.kb[j] = self.kf[j] / self.Keq[j]
            else:
                self.kb[j] = 0.0
                self.Keq[j] = np.inf

    def calculate_effective_pressure(self, rxn):
        """
        Computes the effective pressure for a reaction as:

        .. math:: P_{eff} = P * \\sum_i \\frac{y_i * eff_i}{\\sum_j y_j}

        with:
            - P the pressure of the reactor,
            - y the array of initial moles of the core species

        or as:

        .. math:: P_{eff} = \\frac{P * y_{specific_collider}}{\\sum_j y_j}

        if a specific_collider is mentioned.
        """

        y0_core_species = self.y0[:self.num_core_species]
        sum_core_species = np.sum(y0_core_species)

        j = self.reaction_index[rxn]
        for i in range(self.pdep_collision_reaction_indices.shape[0]):
            if j == self.pdep_collision_reaction_indices[i]:
                # Calculate effective pressure
                if rxn.specific_collider is None:
                    Peff = self.P.value_si * np.sum(self.collider_efficiencies[i] * y0_core_species / sum_core_species)
                else:
                    Peff = self.P.value_si * self.y0[self.species_index[rxn.specific_collider]] / sum_core_species
                return Peff
        return self.P.value_si


class PolymerPhase(object):
    """
    Input container for polymer phase properties.
    """

    def __init__(self,
                 density,
                 initial_moments,
                 initial_explicit,
                 pools,
                 mass_transfer=None,
                 ):
        self.density = density
        self.initial_moments = initial_moments
        self.initial_explicit = initial_explicit
        self.pools = pools
        self.mass_transfer = mass_transfer or list()

    def calculate_volume(self):
        """
        Calculates V_poly = Mass_total / Density.

        Mass Logic:
        Mass_total = Mass(Explicit Species) + Mass(Tails)

        Note: If initial_explicit contains non-polymer species (e.g. dissolved gases),
        their mass contributes to the total phase volume. This assumes 'density'
        refers to the mixture density.
        """
        total_mass_kg = 0.0

        # 1. Add Mass of Explicit Species (Directly)
        for species, moles in self.initial_explicit.items():
            # species.molecular_weight.value_si is in kg/molecule
            total_mass_kg += float(moles) * species.molecular_weight.value_si * constants.Na

        # 2. Add Mass of Tails
        for pool in self.pools:
            label = pool.label
            if label not in self.initial_moments:
                continue

            # Validate Moment Array Shape
            moments = self.initial_moments[label]
            if len(moments) < 2:
                raise ValueError(f"Pool '{label}': initial_moments must provide at least (mu0, mu1).")

            mu0, mu1 = moments[0], moments[1]

            if not pool.monomer:
                if mu1 > 1e-9:
                    raise ValueError(f"Pool '{label}' has moments but no monomer defined.")
                continue

            monomer_mw = pool.monomer.get_molecular_weight()  # kg/mol

            # Calculate explicit contribution to Mu1
            explicit_mu1 = 0.0
            if pool.explicit_map:
                for n_dp, spc in pool.explicit_map.items():
                    if spc in self.initial_explicit:
                        explicit_mu1 += float(n_dp) * float(self.initial_explicit[spc])

            # Sanity Check: Explicit mass cannot exceed Total mass
            tail_mu1 = mu1 - explicit_mu1
            if tail_mu1 < -1e-12:
                raise ValueError(
                    f"Polymer pool '{label}': Explicit mass (mu1={explicit_mu1:.3e}) exceeds "
                    f"Total defined moments (mu1={mu1:.3e}). Check inputs."
                )

            total_mass_kg += max(0.0, tail_mu1) * monomer_mw

        rho_kg_m3 = self.density.value_si
        if rho_kg_m3 <= 0.0:
            raise ValueError(f"Polymer density must be positive, got {rho_kg_m3}.")

        return total_mass_kg / rho_kg_m3

    def get_gas_mask(self, core_species) -> np.ndarray:
        """
        Returns boolean array (True=Gas, False=Polymer).
        Uses ID checks with Label fallback for robustness against species copying.
        Warns if duplicate labels prevent reliable fallback.
        """
        poly_ids = set()
        poly_labels = set()

        def register(spc):
            if spc:
                poly_ids.add(id(spc))
                if hasattr(spc, 'label') and spc.label:
                    poly_labels.add(spc.label)

        # A. Explicit Initials — register only species that are genuinely
        #    polymer-phase. initial_explicit also holds gas-phase solvents and
        #    dissolved gases (e.g. N2).
        #
        #    The is_polymer_proxy tag MUST NOT be used as the discriminator
        #    here: family.py:1657 blanket-stamps it onto EVERY structure of a
        #    proxy-touching reaction, and model.py:486 propagates the tag to the
        #    Species. A gas solvent that merely reacts with the proxy therefore
        #    acquires is_polymer_proxy=True during reaction generation -- AFTER
        #    the build-time gas_species_mask was computed but BEFORE the solver
        #    rebuilds the prospective mask. Keying section A on that mutable tag
        #    makes the SAME core species (e.g. N2) classify GAS at build and
        #    CONDENSED at simulate, which is exactly the divergence RIDER R1
        #    (solver/polymer.pyx) raises on.
        #
        #    is_moment_dummy is set once at species creation (model.py:387) and
        #    is never propagated, so it is the stable polymer-phase marker.
        #    Genuine proxies / monomers / explicit oligomers / mu-species are
        #    registered structurally in sections B/C below and condensed by the
        #    pool-label override pass, independent of the mutable tag.
        for spc in self.initial_explicit.keys():
            if getattr(spc, 'is_moment_dummy', False):
                register(spc)

        # B. Pool Definitions
        for pool in self.pools:
            register(pool.monomer)
            if pool.explicit_map:
                for spc in pool.explicit_map.values():
                    register(spc)
            if pool.mu_species:
                for spc in pool.mu_species:
                    register(spc)

        # C. Mass Transfer
        for mt in self.mass_transfer:
            register(mt.poly_species)

        # Check for Label Ambiguity in Core Species
        core_labels = [getattr(s, "label", None) for s in core_species]
        core_labels = [lab for lab in core_labels if lab]
        label_fallback_safe = (len(core_labels) == len(set(core_labels)))

        mask = np.ones(len(core_species), dtype=bool)
        for i, spc in enumerate(core_species):
            # Check ID match
            if id(spc) in poly_ids:
                mask[i] = False
            # Check Label match (only if safe)
            elif label_fallback_safe and (spc.label and spc.label in poly_labels):
                mask[i] = False

        return mask

    def get_condensed_edge_daughter_bases(self, combined_species) -> set:
        """Return the set of BASE labels of qualifying polymer daughters in
        ``combined_species`` (chain(core, edge)). A species qualifies iff it is
        a Polymer, is_polymer_proxy, its base is NOT a static-deck pool label,
        and a complete {base}_mu0/_mu1/_mu2 moment-dummy triplet (each with
        is_moment_dummy=True) is present in combined_species.

        Spec 2026-06-29. Recomputed by the solver over the LIVE combined list on
        every initialize_model (callable, never a frozen set) -- mirrors
        prospective_classifier so it cannot go stale on the engine-reuse path.
        The solver restricts APPLICATION to edge slots; this predicate is
        base-level. is_moment_dummy is the stable marker; is_polymer_proxy is
        over-stamped (family.py:1657) and never gates alone."""
        from rmgpy.polymer import Polymer  # local: avoids an import cycle

        static_pool_labels = {p.label for p in self.pools}

        # Index the confirmed moment-dummy universe: (base, k) present iff a
        # species labelled "{base}_mu{k}" carries is_moment_dummy=True.
        present = {}
        for spc in combined_species:
            label = getattr(spc, "label", None)
            if not label:
                continue
            for k in (0, 1, 2):
                suffix = "_mu%d" % k
                if label.endswith(suffix) and getattr(spc, "is_moment_dummy", False):
                    raw = label[:-len(suffix)]
                    present[(_base_label(raw), k)] = True

        qualifying = set()
        for spc in combined_species:
            if not isinstance(spc, Polymer):
                continue
            if not getattr(spc, "is_polymer_proxy", False):
                continue
            base = _base_label(spc.label)
            if base in static_pool_labels:
                continue
            if all(present.get((base, k), False) for k in (0, 1, 2)):
                qualifying.add(base)
        return qualifying


class PolymerPhaseBlueprint(object):
    """Temporary container for input file settings."""
    def __init__(self,
                 label: str,
                 species: List[str],
                 solvent: str,
                 density: Union[float, tuple] = (1000.0, 'kg/m^3'),
                 mass_transfer: Optional[List['MassTransfer']] = None,
                 ):
        self.label = label
        self.species_labels = species
        self.solvent_label = solvent
        self.density = density
        self.mass_transfer = mass_transfer or list()


def polymer_phase(label: str,
                  species: List[str],
                  solvent: str,
                  density: Union[float, tuple] = (1000.0, 'kg/m^3'),
                  mass_transfer: Optional[List['MassTransfer']] = None):
    """
    Input file helper to define the polymer phase contents.
    Returns a blueprint that hybrid_polymer_reactor will compile into a real PolymerPhase.

    Args:
        mass_transfer (Optional[List[MassTransfer]]): MassTransfer objects describing
            gas <-> polymer-melt transport for dissolved species. Defaults to an empty list.
    """
    return PolymerPhaseBlueprint(label, species, solvent, density, mass_transfer)


def compile_polymer_phase(blueprint: Union[PolymerPhaseBlueprint, PolymerPhase],
                          initial_moles: Dict[Species, float],
                          species_dict: Dict[str, Species]) -> PolymerPhase:
    """
    Converts a Blueprint + Initial Conditions into a fully realized PolymerPhase object.
    Calculates moments, generates pools, and maps species.
    """
    # If the user somehow passed a ready-made object, return it.
    if not isinstance(blueprint, PolymerPhaseBlueprint):
        return blueprint

    # A. Resolve Density
    rho = Quantity(blueprint.density)

    # B. Resolve Phase Species
    phase_species_set = set()
    for label in blueprint.species_labels:
        if label not in species_dict:
            raise ValueError(f"PolymerPhase species '{label}' not defined in species block.")
        phase_species_set.add(species_dict[label])

    # B2. Validate Mass Transfer definitions
    mass_transfer = blueprint.mass_transfer or list()
    for mt in mass_transfer:
        if not isinstance(mt, MassTransfer):
            raise ValueError(
                f"PolymerPhase mass_transfer must contain MassTransfer objects, got {type(mt).__name__}.")

    # C. Compile State
    initial_moments = {}
    initial_explicit = {}
    pools = []

    for spc, moles in initial_moles.items():
        # Skip if this species isn't in the polymer phase list
        if spc not in phase_species_set:
            continue

        # Check for Polymer Type
        # We assume the Species object was enriched with polymer attributes in Phase 1
        if hasattr(spc, 'Mn') and hasattr(spc, 'monomer'):
            # --- CALCULATE MOMENTS ---
            if hasattr(spc, 'monomer_mw_g_mol'):
                monomer_mw = spc.monomer_mw_g_mol / 1000.0
            else:
                # Fallback: Calculate from graph
                monomer_mw = spc.monomer.get_molecular_weight().value_si

            # Create distinct dummy species for moments
            mu0_spc = species_dict[f"{spc.label}_mu0"]
            mu1_spc = species_dict[f"{spc.label}_mu1"]
            mu2_spc = species_dict[f"{spc.label}_mu2"]
            for m_spc in [mu0_spc, mu1_spc, mu2_spc]:
                if m_spc.label not in species_dict:
                    species_dict[m_spc.label] = m_spc

            mu0 = moles
            mn_kg = spc.Mn / 1000.0
            mw_kg = spc.Mw / 1000.0

            dp_n = mn_kg / monomer_mw
            dp_w = mw_kg / monomer_mw

            mu1 = mu0 * dp_n
            mu2 = mu1 * dp_w

            initial_moments[spc.label] = (mu0, mu1, mu2)

            # --- Single source of truth for the t=0 moment state ---
            # initialMoles[proxy] is authoritative for mu0: it is what the
            # solver actually integrates as y0. The Polymer object's .moments,
            # by contrast, are what write_polymer_pools_sidecar serializes
            # (it reads Polymer.moments). If we leave .moments at its
            # distribution-derived value (initial_mass/Mn), the sidecar would
            # report a t=0 state the solver never simulated. Reconcile here so
            # the report MIRRORS the engine: the sidecar reports exactly the
            # moments the solver integrates.
            spc.moments = np.array([mu0, mu1, mu2], dtype=np.float64)

            # Warn if the pool's stated loading (initial_mass/Mn) implies a
            # different chain count (mu0) than initialMoles. Both are valid
            # ways to state the t=0 amount, but they must agree; if they don't,
            # the solver silently honors initialMoles and the initial_mass +
            # Mn/Mw in the deck misrepresent the simulated state.
            initial_mass_g = getattr(spc, 'initial_mass_g', None)
            if initial_mass_g and spc.Mn:
                mu0_from_mass = initial_mass_g / spc.Mn
                if abs(mu0_from_mass - mu0) > 1e-6 * max(abs(mu0_from_mass), abs(mu0)):
                    logging.warning(
                        "Polymer pool '%s': initial_mass/Mn implies mu0=%g chains, but "
                        "initialMoles gives mu0=%g. These disagree; the solver integrates "
                        "the initialMoles value (mu0=%g) and the sidecar mirrors it. "
                        "Reconcile initial_mass/Mn/Mw with initialMoles in the input deck.",
                        spc.label, mu0_from_mass, mu0, mu0)

            # Create Pool Config
            # Note: Using the species itself as the placeholder for moments (Integration Test Hack)
            # In production, this needs distinct dummy species.
            pool = PolymerPool(
                label=spc.label,
                xs=spc.cutoff,
                monomer=spc.monomer,
                explicit_map={},
                mu_species=[mu0_spc, mu1_spc, mu2_spc],
                k_scission=spc.k_scission,
                k_unzip=spc.k_unzip,
                proxy_species=spc,
                monomer_product=getattr(spc, 'monomer_product_species', None),
            )
            pools.append(pool)
        else:
            # Standard Solvents / Dissolved Gases
            initial_explicit[spc] = moles

    # D. Instantiate Real Object
    return PolymerPhase(
        density=rho,
        initial_moments=initial_moments,
        initial_explicit=initial_explicit,
        pools=pools,
        mass_transfer=mass_transfer
    )


class PolymerPool(object):
    """
    Input class for defining a polymer pool configuration.

    This class configures the hybrid Method of Moments (HMOM) representation for a
    specific polymer type (e.g., Polyethylene, Polystyrene). It defines the boundary
    between explicit oligomers and the statistical tail, as well as the kinetic
    parameters driving the distribution dynamics.

    Args:
        label (str): A unique name for this polymer pool (e.g., 'PE', 'PS').
                     Used for logging and identification.
        xs (int): The hybrid cutoff index. Chains with length n <= xs are treated as
                  explicit chemical species. Chains with n > xs are tracked statistically
                  via moments.
        monomer (Species): The RMG Species object representing the monomer unit.
                           Used to calculate molecular weights and mass balances.
        explicit_map (Dict[int, Species]): A dictionary mapping degree of polymerization (DP)
                                           to explicit Species objects.
                                           Format: {1: Monomer, 2: Dimer, ..., xs: Oligomer_xs}.
        mu_species (List[Species]): A list of exactly three Species objects representing the
                                    statistical moments [Mu0, Mu1, Mu2]. These are placeholder
                                    species used by the solver to track the moment values.
        k_scission (float, optional): The random scission rate coefficient [1/s].
                                      Defaults to 0.0.
        k_unzip (float, optional): The chain-end scission (unzipping) rate coefficient [1/s].
                                   This parameter drives the physical flux ('handshake')
                                   from the statistical tail into the explicit oligomers.
                                   Defaults to 0.0.
    """
    def __init__(self,
                 label: str,
                 xs: int,
                 monomer: Species,
                 explicit_map: Dict[int, Species],
                 mu_species: List[Species],
                 k_scission: float = 0.0,
                 k_unzip: float = 0.0,
                 proxy_species: Optional[Species] = None,
                 monomer_product: Optional[Species] = None,
                 ):
        self.label = label
        self.xs = xs
        self.monomer = monomer
        self.explicit_map = explicit_map
        self.mu_species = mu_species
        self.k_scission = k_scission
        self.k_unzip = k_unzip
        self.proxy_species = proxy_species
        # Real reactive species (e.g. styrene) that chain-unzip releases into the
        # melt. Distinct from `monomer` (the repeat-unit Molecule used for mass).
        self.monomer_product = monomer_product

    def to_config(self, spc_map):
        """
        Converts Input Object -> Solver Config (resolving indices using pre-built map).
        """
        # 1. Resolve Explicit Map Indices
        explicit_indices = dict()
        if self.explicit_map:
            for dp, spc in self.explicit_map.items():
                if dp > self.xs:
                    raise ValueError(f"Pool '{self.label}': explicit_map contains DP={dp} > xs={self.xs}.")

                if spc in spc_map:
                    explicit_indices[dp] = spc_map[spc]
                else:
                    raise ValueError(f"Pool {self.label}: Explicit species for DP={dp} ({spc}) not in core species.")

        # 2. Resolve Moment Indices
        if len(self.mu_species) != 3:
            raise ValueError(f"Pool {self.label}: mu_species must contain exactly 3 species objects.")

        try:
            mu_idxs = tuple(spc_map[s] for s in self.mu_species)
        except KeyError as e:
            raise ValueError(f"Pool {self.label}: Moment species {e} missing from core list.")

        # 3. Resolve Monomer Index
        # The released-monomer target must be a real core Species. `self.monomer`
        # is the repeat-unit Molecule (never a key in the Species-keyed spc_map),
        # so resolve from `monomer_product` when provided.
        monomer_idx = None
        if self.monomer_product is not None:
            monomer_idx = spc_map.get(self.monomer_product)
            if monomer_idx is None:
                raise ValueError(
                    f"Pool {self.label}: monomer_product {self.monomer_product} not in core species; "
                    f"cannot wire unzip-to-monomer release.")

        # HARD ERROR: a negative k_unzip is not a valid rate constant. Every
        # downstream unzip consumer is gated on k_unzip > 0, so a negative
        # value would silently become an inert channel instead of failing --
        # and a wired monomer_product must not dodge the sign check.
        if self.k_unzip < 0.0:
            raise ValueError(
                f"Pool {self.label}: k_unzip={self.k_unzip:g} is negative -- a "
                f"negative k_unzip is not a valid rate constant. Set k_unzip >= 0.")

        # HARD ERROR: k_unzip > 0 with no resolvable monomer emission target.
        # The solver drains the condensed moments unconditionally when
        # k_unzip > 0 (polymer.pyx: dmu1_dt -= k_unzip*mu0) but only emits the
        # released monomer when monomer_poly_index is not None -- so this
        # config shape would silently un-conserve mass (the drained repeat
        # units go nowhere). Refuse at configuration time, before the solver.
        if self.k_unzip > 0.0 and monomer_idx is None:
            raise ValueError(
                f"Pool {self.label}: k_unzip={self.k_unzip:g} > 0 but no monomer_product "
                f"is defined for this pool. The unzip channel would drain the condensed "
                f"moments (mu1/mu2) with no released-monomer emission target, leaving "
                f"mass silently un-conserved. Define monomer_product (the real monomer "
                f"species released on unzip) or set k_unzip=0.")

        # 4. Monomer (repeat-unit) MW [g/mol] for the spawn-gate snapshot AND the
        #    reference-state tripwire chain_window (spec 2026-06-10 §3, same idiom
        #    as Polymer.monomer_mw_g_mol). self.monomer is normally a Molecule (the
        #    polymer() input helper / Polymer._validate_monomer), but may be a
        #    Species when resolved from species_dict; handle BOTH. A Species carries
        #    a .molecule list; a Molecule answers get_molecular_weight() directly.
        #    Reading only the Species idiom left monomer_mw_g_mol=0 for Molecule
        #    monomers, collapsing chain_window to the slack and leaking small gas
        #    fragments into the melt reference-state sum. Best-effort: 0.0 (-> the
        #    gate defers) when no resolvable structure.
        monomer_mw_g_mol = 0.0
        mol_list = getattr(self.monomer, "molecule", None)
        try:
            if mol_list:
                monomer_mw_g_mol = mol_list[0].get_molecular_weight() * 1000.0
            elif self.monomer is not None and hasattr(self.monomer, "get_molecular_weight"):
                monomer_mw_g_mol = self.monomer.get_molecular_weight() * 1000.0
        except Exception:
            monomer_mw_g_mol = 0.0

        return PolymerPoolConfig(
            label=self.label,
            xs=self.xs,
            explicit_dp_to_species_index=explicit_indices,
            mu_indices=mu_idxs,
            monomer_poly_index=monomer_idx,
            monomer_mw_g_mol=monomer_mw_g_mol,
            k_scission=self.k_scission,
            k_unzip=self.k_unzip
        )


def derive_daughter_pool_configs(core_species, spc_map, existing_pool_labels):
    """
    Build a :class:`PolymerPoolConfig` for each daughter Polymer registered as a
    core species but not covered by a static deck pool.

    A scission tail or spawn-intent (``<parent>_d{n}``) daughter is registered by
    ``CoreEdgeReactionModel._register_polymer`` as a core Polymer species plus its
    own ``_mu0`` / ``_mu1`` / ``_mu2`` moment-dummy species -- but ``pool_configs``
    is built only from the static deck list ``polymerPhase.pools``. Without a pool
    config the daughter's species map to ``-1`` in ``species_to_pool_indices``, so
    a stamped SCISSION_FRAGMENT / MIGRATION reaction targeting it demotes to
    UNRESOLVED (``polymer.pyx`` "could not resolve their solver pool(s)"). Deriving
    the config from the registered core species closes that gap uniformly for both
    daughter paths; the solver binds the pool by label on the next rebuild.

    Daughters spawn honest-empty (no explicit oligomers); only label + xs +
    moment indices are needed for the solver to resolve the pool.

    Labels use the solver's binding convention: RMG appends a ``(N)`` index to
    registered species labels (a proxy displays as ``PS(2)`` while its dummies
    stay the clean ``PS_mu0``), and ``_apply_pool_phase_overrides`` binds on
    ``label.partition('(')[0] == pool.label``. So the derived pool label and the
    moment-dummy lookup both use the ``(``-stripped base label.
    """
    from rmgpy.polymer import Polymer

    def base(label):
        return label.partition('(')[0]

    static = set(existing_pool_labels)
    # Base-label -> core index, using the same '('-stripping the solver binds with.
    base_to_index = {}
    for spc in core_species:
        base_to_index.setdefault(base(spc.label), spc_map[spc])

    configs = []
    seen = set()
    for spc in core_species:
        if not isinstance(spc, Polymer):
            continue
        b = base(spc.label)
        if b in static or b in seen:
            # Static deck pool's proxy (never re-derived) or a base already
            # configured (the proxy plus its index-suffixed copies).
            continue
        mu_indices = tuple(base_to_index.get(f"{b}_mu{k}") for k in (0, 1, 2))
        if any(i is None for i in mu_indices):
            # Daughter registered without its full moment-dummy triplet in the
            # core map; skip rather than build an unresolvable pool.
            continue
        seen.add(b)
        # Carry the daughter's own monomer MW [g/mol] (same repeat unit as the
        # parent; the daughter Polymer computed it in __init__). Omitting it left
        # the config at 0.0, which drags max(monomer_mw over pools) -> 0 and
        # collapses the reference-state tripwire chain_window to the bare slack,
        # leaking small gas scission fragments into the melt sum.
        configs.append(PolymerPoolConfig(
            label=b,
            xs=spc.cutoff,
            explicit_dp_to_species_index={},
            mu_indices=mu_indices,
            monomer_mw_g_mol=float(getattr(spc, "monomer_mw_g_mol", 0.0) or 0.0),
        ))
    return configs


class MassTransfer(object):
    """
    Input class for Mass Transfer definition.

    Defines the transport of a specific chemical species between the gas phase headspace
    and the polymer melt phase. The flux is driven by the concentration difference relative
    to equilibrium: J = kLa * (C_poly - K * C_gas).

    Args:
        gas_species (Species): The RMG Species object representing the component in the gas phase.
        poly_species (Species): The RMG Species object representing the component dissolved in the polymer phase.
        K (Union[float, Quantity]): The partition coefficient (Equilibrium Constant), defined as
                                    K = C_poly_eq / C_gas_eq. Dimensionless.
        kLa (Union[float, Quantity]): The volumetric mass transfer coefficient [1/s].
    """

    def __init__(self,
                 gas_species: Species,
                 poly_species: Species,
                 K: Union[float, Quantity],
                 kLa: Union[float, Quantity],
                 ):
        self.gas_species = gas_species
        self.poly_species = poly_species
        self.K = K
        self.kLa = kLa

    def to_config(self, spc_map):
        """
        Converts Input Object -> Solver Config (resolving indices using pre-built map).
        """
        if self.gas_species not in spc_map:
            raise ValueError(f"MassTransfer gas species '{self.gas_species}' not found in core species.")

        if self.poly_species not in spc_map:
            raise ValueError(f"MassTransfer polymer species '{self.poly_species}' not found in core species.")

        gas_index = spc_map[self.gas_species]
        poly_index = spc_map[self.poly_species]

        # Handle Quantities if present
        K_val = self.K.value_si if hasattr(self.K, 'value_si') else float(self.K)
        kLa_val = self.kLa.value_si if hasattr(self.kLa, 'value_si') else float(self.kLa)

        # Enforce Physical Bounds
        if K_val <= 0.0:
            raise ValueError(f"MassTransfer K (partition coeff) must be > 0, got {K_val}.")
        if kLa_val < 0.0:
            raise ValueError(f"MassTransfer kLa must be >= 0, got {kLa_val}.")

        return MassTransferConfig(
            gas_index=gas_index,
            poly_index=poly_index,
            K=K_val,
            kLa=kLa_val,
        )
