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

"""
Biphasic hybrid reactor solver for polymer pyrolysis/degradation in RMG-Py.

Phase Policy: "Produce-then-Transfer"
1. Reactions are strictly phase-pure.
   - Gas reactants -> Gas products.
   - Polymer reactants -> Polymer products.
   - Cross-phase core products disqualify a reaction (rate=0).
   - Edge products (not in solver state) allowed for diagnostics.
2. Mass Transfer is explicit via kLa driving force.

Hybrid Representation:
- Explicit Oligomers: Species tracking chains with degree of polymerization (DP) <= xs.
- Moment Tail: Statistical moments (mu0, mu1, mu2) for chains with DP > xs.
- Hybrid Handshake: Dynamic flux from Tail -> Explicit at boundary xs.

Units:
- Time: s
- Residual: mol/s (dn/dt - dydt)
- Moments (Input/State): Moles of Moments (Mu_k = mu_k * V_poly)
- Tail kinetics interface: tail_kinetics receives μk/V_poly [mol/m^3] and returns d(μk/V_poly)/dt [mol/m^3/s]
  (and small-species sources in mol/m^3/s). The solver multiplies by V_poly when accumulating into dn_dt [mol/s].
- Diagnostics: Core rates in mol/m^3_phase/s (or just 'rate' units)
"""

from __future__ import annotations

import itertools
import logging
import math
from dataclasses import dataclass
from typing import Callable, Dict, List, Optional, Tuple

cimport cython
import numpy as np
cimport numpy as np

import rmgpy.constants as constants
from rmgpy.quantity import Quantity
from rmgpy.solver.base import ReactionSystem


SMALL_EPS = 1e-30
LN_EXP_OVERFLOW_GUARD = 700.0
TAIL_CONC_MIN = 1e-9  # Minimum concentration (mol/m^3) to actuate handshake


# ======================================================================================
# DATA CONFIGURATION
# ======================================================================================

@dataclass(frozen=True)
class PolymerPoolConfig:
    """
    Configuration for one polymer pool.
    Indices must refer to valid CORE species indices.
    """
    label: str
    xs: int
    explicit_dp_to_species_index: Dict[int, int]
    mu_indices: Tuple[int, int, int]
    monomer_poly_index: Optional[int] = None

    # Kinetics Parameters
    k_scission: float = 0.0

    # Recession Rate [1/s]
    # REQUIRED for the Hybrid Handshake (Tail -> Explicit flux).
    # This parameter sets the timescale for both the default recession chemistry
    # (n -> n-1) AND the physical handshake flux across the boundary xs.
    k_unzip: float = 0.0

    # Custom Kinetics Hook
    # f(T, P, mu0, mu1, mu2, mu3) -> (dmu0_dt, dmu1_dt, dmu2_dt, small_species_sources)
    # Units:
    #   mu* are intensive moments [mol/m^3_poly]
    #   dmu*/dt are [mol/m^3_poly/s]
    #   small_species_sources: Dict[core_index -> mol/m^3_poly/s]
    # The solver multiplies these volumetric rates by V_poly to update extensive state [mol] / [mol/s].
    # Contract: If provided, this handles "Chemical Event" derivatives.
    # It does NOT handle the "Physical Handshake" flux.
    tail_kinetics: Optional[Callable[[float, float, float, float, float, float],
                                     Tuple[float, float, float, Dict[int, float]]]] = None

    def __repr__(self):
        return f"<PolymerPool '{self.label}' xs={self.xs}>"


@dataclass(frozen=True)
class MassTransferConfig:
    """
    Explicit mass transfer between a gas-phase species and its dissolved counterpart.
    """
    gas_index: int
    poly_index: int
    K: float    # K = C_poly_eq / C_gas_eq
    kLa: float  # [1/s]


# ======================================================================================
# HELPERS
# ======================================================================================

def _safe_mu3_from_mu012(mu0: float, mu1: float, mu2: float) -> float:
    """Log-Lagrange closure."""
    if mu0 <= SMALL_EPS or mu1 <= SMALL_EPS or mu2 <= SMALL_EPS:
        return 0.0
    with np.errstate(divide="raise", invalid="raise", over="raise"):
        try:
            ln_mu3 = 3.0 * np.log(mu2) - 3.0 * np.log(mu1) + np.log(mu0)
            if ln_mu3 > LN_EXP_OVERFLOW_GUARD:
                return float("inf")
            return float(np.exp(ln_mu3))
        except (FloatingPointError, ValueError):
            return 0.0

def _explicit_moment_contributions(V_poly: float,
                                   y: np.ndarray,
                                   explicit_dp_to_species_index: Dict[int, int],
                                   ) -> Tuple[float, float, float]:
    """Calculate explicit oligomer contributions to moments."""
    if V_poly <= 0: return 0.0, 0.0, 0.0
    mu0, mu1, mu2 = 0.0, 0.0, 0.0

    for n, idx in explicit_dp_to_species_index.items():
        if n <= 0: continue
        if not (0 <= idx < len(y)):
            raise IndexError(f"Explicit species index {idx} out of bounds.")

        Nn = max(0.0, y[idx]) / V_poly
        mu0 += Nn
        mu1 += n * Nn
        mu2 += (n * n) * Nn
    return mu0, mu1, mu2


def _gamma_params_from_mu012(mu0: float, mu1: float, mu2: float) -> Optional[Tuple[float, float]]:
    if mu0 <= SMALL_EPS or mu1 <= SMALL_EPS or mu2 <= SMALL_EPS:
        return None
    mean = mu1 / mu0
    if mean <= SMALL_EPS: return None
    pdi = (mu2 * mu0) / (mu1 * mu1)

    # Monodisperse fallback trigger
    if not np.isfinite(pdi) or pdi <= 1.0 + 1e-6:
        return None

    k = 1.0 / (pdi - 1.0)
    theta = mean / k
    if k <= 0.0 or theta <= 0.0 or not np.isfinite(k) or not np.isfinite(theta):
        return None
    return float(k), float(theta)


def _gamma_prob_conditional_hybrid(m: int, xs: int, k: float, theta: float) -> float:
    if m <= 0:
        return 0.0
    try:
        from scipy.special import gammainc
        def cdf(x):
            return float(gammainc(k, x/theta))
        F_b = cdf(m + 0.5)
        F_a = cdf(max(0.0, m - 0.5))
        F_cut = cdf(xs + 0.5)

        bin_prob = max(0.0, F_b - F_a)
        tail_prob = max(0.0, 1.0 - F_cut)
        if tail_prob <= 1e-12: return 0.0
        return bin_prob / tail_prob
    except ImportError:
        return _discrete_gamma_fallback(m, xs, k, theta)


def _discrete_gamma_fallback(target: int, xs: int, k: float, theta: float) -> float:
    if target <= xs:
        return 0.0
    def logpdf(n):
        x = float(n)
        if x <= 0: return -1e9
        return (k - 1.0) * math.log(x) - (x / theta) - k * math.log(theta) - math.lgamma(k)

    mean = k * theta
    sigma = math.sqrt(k * theta * theta)
    nmax = int(max(xs + 1, mean + 6.0 * sigma + 20.0))
    nmax = min(nmax, xs + 2000)

    logs = []
    target_idx = -1
    for i, n in enumerate(range(xs + 1, nmax + 1)):
        val = logpdf(n)
        logs.append(val)
        if n == target:
            target_idx = i

    if not logs:
        return 0.0
    max_log = max(logs)
    denom = sum(math.exp(v - max_log) for v in logs)
    if denom <= 0.0 or target_idx == -1:
        return 0.0
    return math.exp(logs[target_idx] - max_log) / denom


# ======================================================================================
# SOLVER CLASS
# ======================================================================================

class HybridPolymerSystem(ReactionSystem):
    """
    Biphasic Hybrid Reactor Solver.

    Residual returns mol/s.
    Enforces 'Produce-then-Transfer' policy for interphase transport.
    """

    def __init__(
        self,
        T,
        P,
        initial_mole_fractions: Dict,
        V_poly: float,
        polymer_pools: List[PolymerPoolConfig],
        mass_transfer: Optional[List[MassTransferConfig]] = None,
        polymer_species_labels=None,
        gas_species_mask: Optional[np.ndarray] = None,
        constant_gas_volume: bool = False,
        V_gas0: Optional[float] = None,
        initial_polymer_moments: Optional[Dict[str, Tuple[float, float, float]]] = None,
        initial_explicit_species: Optional[Dict[str, Dict[int, float]]] = None,
        termination=None,
        sensitive_species=None,
        sensitivity_threshold: float = 1e-3,
        sens_conditions=None,
        const_spc_names: Optional[List[str]] = None,
        pdep_collision_reaction_indices: Optional[np.ndarray] = None,
        pdep_collider_kinetics: Optional[List] = None,
        collider_efficiencies: Optional[np.ndarray] = None,
    ):
        super().__init__(termination=termination,
                         sensitive_species=sensitive_species,
                         sensitivity_threshold=sensitivity_threshold)
        cdef np.ndarray species_to_pool_indices

        self.T = Quantity(T)
        self.P = Quantity(P)
        self.initial_mole_fractions = initial_mole_fractions

        self.V_poly = float(V_poly)
        if self.V_poly <= 0.0:
            raise ValueError(f"V_poly must be > 0, got {self.V_poly}")

        self.constant_gas_volume = bool(constant_gas_volume)
        self.V_gas0 = float(V_gas0) if V_gas0 is not None else None
        self.V_gas = 0.0
        self.V = 0.0

        self.polymer_pools = list(polymer_pools)
        self.mass_transfer = list(mass_transfer) if mass_transfer else []
        self.gas_species_mask = gas_species_mask
        self.const_spc_names = const_spc_names or []
        self.const_spc_indices = None
        self.sens_conditions = sens_conditions

        self.pdep_collision_reaction_indices = pdep_collision_reaction_indices if pdep_collision_reaction_indices is not None else np.array([], int)
        self.pdep_collider_kinetics = pdep_collider_kinetics if pdep_collider_kinetics is not None else []
        self.collider_efficiencies = collider_efficiencies if collider_efficiencies is not None else np.array([[]], float)
        self.initial_polymer_moments = initial_polymer_moments or {}
        self.initial_explicit_species = initial_explicit_species or {}
        self.polymer_species_labels = set(polymer_species_labels) if polymer_species_labels else set()

        self.jacobian_matrix = None

        self._scratch_C_gas = None
        self._scratch_C_poly = None
        self._scratch_dn_dt = None

    def validate_configuration(self):
        """Strict validation of indices and masks."""
        n_core = self.num_core_species

        if self.gas_species_mask is None:
            raise ValueError("gas_species_mask must be set before validation.")

        if len(self.gas_species_mask) != n_core:
            raise ValueError(f"Gas mask len {len(self.gas_species_mask)} != num_core {n_core}")

        moment_indices = set()

        for pool in self.polymer_pools:
            for idx in pool.mu_indices:
                if not (0 <= idx < n_core):
                    raise ValueError(f"Pool {pool.label} moment index {idx} out of range.")
                if self.gas_species_mask[idx]:
                    raise ValueError(f"Pool {pool.label} moment index {idx} masked as GAS.")
                moment_indices.add(idx)

            for dp, idx in pool.explicit_dp_to_species_index.items():
                if not (0 <= idx < n_core):
                    raise ValueError(f"Pool {pool.label} explicit DP={dp} index {idx} out of range.")
                if self.gas_species_mask[idx]:
                    raise ValueError(f"Pool {pool.label} explicit DP={dp} index {idx} masked as GAS.")

            if pool.monomer_poly_index is not None:
                if not (0 <= pool.monomer_poly_index < n_core):
                    raise ValueError(f"Pool {pool.label} monomer index out of range.")
                if self.gas_species_mask[pool.monomer_poly_index]:
                    raise ValueError(f"Pool {pool.label} monomer index masked as GAS.")

        for mt in self.mass_transfer:
            if not (0 <= mt.poly_index < n_core):
                raise ValueError(f"Mass transfer poly_index {mt.poly_index} out of range.")
            if not (0 <= mt.gas_index < n_core):
                raise ValueError(f"Mass transfer gas_index {mt.gas_index} out of range.")

            if self.gas_species_mask[mt.poly_index]:
                raise ValueError(f"Mass transfer poly_index {mt.poly_index} is masked as GAS.")
            if not self.gas_species_mask[mt.gas_index]:
                raise ValueError(f"Mass transfer gas_index {mt.gas_index} is NOT masked as GAS.")

        # Moment isolation check
        ir = self.reactant_indices
        ip = self.product_indices
        all_rxn_indices = {i for i in itertools.chain(ir.flatten(), ip.flatten()) if i >= 0}
        overlap = moment_indices.intersection(all_rxn_indices)
        if overlap:
            raise ValueError(f"Configuration Error: Moment indices {overlap} appear in reaction stoichiometry. "
                             f"Moments must evolve only via tail_kinetics.")

    def get_const_spc_indices(self, core_species):
        """Identify indices of constant species."""
        if not self.const_spc_names:
            return
        self.const_spc_indices = list()
        for name in self.const_spc_names:
            for spc in core_species:
                if spc.label == name:
                    self.const_spc_indices.append(core_species.index(spc))
                    break

    def convert_initial_keys_to_species_objects(self, species_dict):
        """Convert initial mole fraction keys to species objects."""
        initial_mole_fractions = dict()
        for label, moleFrac in self.initial_mole_fractions.items():
            initial_mole_fractions[species_dict[label]] = moleFrac
        self.initial_mole_fractions = initial_mole_fractions

        conditions = dict()
        if self.sens_conditions is not None:
            for label, value in self.sens_conditions.items():
                if label in ("T", "P"):
                    conditions[label] = value
                else:
                    conditions[species_dict[label]] = value
        self.sens_conditions = conditions

    def initialize_model(self, core_species, core_reactions, edge_species, edge_reactions,
                      surface_species=None, surface_reactions=None, pdep_networks=None,
                      atol=1e-16, rtol=1e-8, sensitivity=False, sens_atol=1e-6, sens_rtol=1e-4,
                      filter_reactions=False, conditions=None, **kwargs):
        """
        Initialize the polymer hybrid reactor model.
        """
        ReactionSystem.initialize_model(
            self, core_species, core_reactions, edge_species, edge_reactions,
            surface_species, surface_reactions, pdep_networks,
            atol, rtol, sensitivity, sens_atol, sens_rtol, filter_reactions, conditions
        )

        cdef int n_core = self.num_core_species
        if n_core <= 0:
            raise ValueError(f"Solver received an empty core species list (n_core={n_core}).")

        if len(core_species) != n_core:
            raise ValueError(
                f"Core species length mismatch: len(core_species)={len(core_species)} != num_core_species={n_core}. "
                "Refusing to proceed with inconsistent allocations."
            )

        if self.gas_species_mask is None:
            mask = np.ones(n_core, dtype=bool)
            poly_labels = getattr(self, "polymer_species_labels", set())
            for i in range(n_core):
                if core_species[i].label in poly_labels:
                    mask[i] = False
            self.gas_species_mask = mask
        else:
            if len(self.gas_species_mask) != n_core:
                raise ValueError(
                    f"Provided gas_species_mask length {len(self.gas_species_mask)} != n_core {n_core}. "
                    "Mask must match core size exactly."
                )
            self.gas_species_mask = np.asarray(self.gas_species_mask, dtype=bool)

        self.get_const_spc_indices(core_species)
        self.set_initial_conditions()

        self._scratch_C_gas = np.zeros(n_core, float)
        self._scratch_C_poly = np.zeros(n_core, float)
        self._scratch_dn_dt = np.zeros(n_core, float)

        if filter_reactions:
            ReactionSystem.set_initial_reaction_thresholds(self)
        self.generate_rate_coefficients(core_reactions, edge_reactions)
        ReactionSystem.compute_network_variables(self, pdep_networks)

        ReactionSystem.set_initial_derivative(self)
        ReactionSystem.initialize_solver(self)


    def set_initial_conditions(self):
        """
        Sets initial state y0.

        WARNING: The dictionary `initial_mole_fractions` is named by the base class contract,
        but for this reactor type, the values are interpreted as MOLES of gas species.
        If `constant_gas_volume=False`, `V_gas` is inferred from the total gas moles via the
        ideal gas law at (T, P).
        """
        cdef int n_core = self.num_core_species

        ReactionSystem.set_initial_conditions(self)

        if self.gas_species_mask is None or self.gas_species_mask.shape[0] != n_core:
            raise ValueError(f"CRITICAL DIMENSION MISMATCH: gas_species_mask size "
                             f"({None if self.gas_species_mask is None else self.gas_species_mask.shape[0]}) "
                             f"does not match num_core_species ({n_core}).")

        # 1. Gas Species (Interpret values as MOLES)
        for spec, val_moles in self.initial_mole_fractions.items():
            i = self.get_species_index(spec)
            self.y0[i] = val_moles

        # 2. Polymer Explicit Species
        for pool in self.polymer_pools:
            if pool.label in self.initial_explicit_species:
                dp_map = self.initial_explicit_species[pool.label]
                for dp, moles in dp_map.items():
                    if dp in pool.explicit_dp_to_species_index:
                        idx = pool.explicit_dp_to_species_index[dp]
                        self.y0[idx] = max(0.0, moles)
            elif pool.label not in self.initial_polymer_moments:
                logging.warning(f"Polymer pool '{pool.label}' has no initial conditions provided (Explicit or Moments).")

        # 3. Polymer Moments
        for pool in self.polymer_pools:
            if pool.label in self.initial_polymer_moments:
                moms = self.initial_polymer_moments[pool.label]
                for k in range(3):
                    idx = pool.mu_indices[k]
                    self.y0[idx] = max(0.0, moms[k])

        # 4. Initialize Gas Volume
        mask = self.gas_species_mask
        if self.constant_gas_volume:
            if self.V_gas0 is None or self.V_gas0 <= 0.0:
                raise ValueError("constant_gas_volume=True requires V_gas0 > 0")
            self.V_gas = self.V_gas0
        else:
            # Sum moles of gas species only
            n_gas0 = float(np.sum(self.y0[:n_core][mask]))
            self.V_gas = constants.R * self.T.value_si * n_gas0 / self.P.value_si if n_gas0 > 0 else 1.0
        self.V = self.V_gas + self.V_poly

        # 5. Set Concentrations
        for j in range(n_core):
            if mask[j]:
                self.core_species_concentrations[j] = self.y0[j] / self.V_gas
            else:
                self.core_species_concentrations[j] = self.y0[j] / self.V_poly

        # 6. Moment Consistency Check
        tol = 1e-12
        for pool in self.polymer_pools:
            idx_mu0, idx_mu1, idx_mu2 = pool.mu_indices

            mu0_tot = max(0.0, self.y0[idx_mu0]) / self.V_poly
            mu1_tot = max(0.0, self.y0[idx_mu1]) / self.V_poly
            mu2_tot = max(0.0, self.y0[idx_mu2]) / self.V_poly

            mu0_exp, mu1_exp, mu2_exp = _explicit_moment_contributions(
                V_poly=self.V_poly,
                y=self.y0,
                explicit_dp_to_species_index=pool.explicit_dp_to_species_index)

            r0 = mu0_tot - mu0_exp
            r1 = mu1_tot - mu1_exp
            r2 = mu2_tot - mu2_exp

            if r1 < -tol:
                logging.warning(f"Polymer pool {pool.label}: Explicit mass exceeds Total. Clamping tail to 0.")

            self.y0[idx_mu0] = max(0.0, r0) * self.V_poly
            self.y0[idx_mu1] = max(0.0, r1) * self.V_poly
            self.y0[idx_mu2] = max(0.0, r2) * self.V_poly

    def generate_rate_coefficients(self, core_reactions, edge_reactions):
        for rxn in itertools.chain(core_reactions, edge_reactions):
            j = self.reaction_index[rxn]
            self.kf[j] = rxn.get_rate_coefficient(self.T.value_si, self.P.value_si)
            if rxn.reversible:
                self.Keq[j] = rxn.get_equilibrium_constant(self.T.value_si)
                self.kb[j] = self.kf[j] / self.Keq[j]
            else:
                self.kb[j] = 0.0
                self.Keq[j] = np.inf

    def get_threshold_rate_constants(self, model_settings):
        """
        Get the theoretical maximum rate constants for reaction filtering.
        These are used to normalize edge reaction fluxes.
        """
        # Unimolecular limit: kB * T / h (Transition State Theory limit)
        unimolecular_threshold_rate_constant = (constants.kB / constants.h) * self.T.value_si

        # Bimolecular limit: User-defined (usually 1e8 - 1e12 m^3/mol/s)
        bimolecular_threshold_rate_constant = model_settings.filter_threshold

        # Trimolecular limit: Adjusted for volume units
        trimolecular_threshold_rate_constant = model_settings.filter_threshold / 1e3

        return (unimolecular_threshold_rate_constant,
                bimolecular_threshold_rate_constant,
                trimolecular_threshold_rate_constant)

    @cython.boundscheck(False)
    def residual(self, double t, np.ndarray[np.float64_t, ndim=1] y,
                 np.ndarray[np.float64_t, ndim=1] dydt,
                 np.ndarray[np.float64_t, ndim=1] senpar = np.zeros(1, float)):
        """
        Compute the residual (dn/dt - dydt) at time t for state y.
        """
        cdef int n_core = len(self.core_species_rates)
        cdef int n_rxn = len(self.core_reaction_rates)
        cdef int p0_idx, p1_idx, p2_idx, prod_p_idx, p_slot, p_idx_tmp
        cdef double mu1_0, mu1_1, rf_hijack

        if self.gas_species_mask.shape[0] != n_core:
            raise ValueError(f"State/Mask mismatch: y={n_core}, mask={self.gas_species_mask.shape[0]}")

        # Disable caching
        self.jacobian_matrix = None

        ir = self.reactant_indices
        ip = self.product_indices
        inet = self.network_indices
        knet = self.network_leak_coefficients

        n_core = len(self.core_species_rates)
        n_rxn = len(self.core_reaction_rates)
        n_edge = len(self.edge_species_rates)

        # 1. Lazy Allocation (Safety check)
        if self._scratch_C_gas is None or len(self._scratch_C_gas) != n_core:
            self._scratch_C_gas = np.zeros(n_core, float)
            self._scratch_C_poly = np.zeros(n_core, float)
            self._scratch_dn_dt = np.zeros(n_core, float)

        # 2. Update Volumes
        if self.constant_gas_volume:
            V_gas = self.V_gas0
        else:
            n_gas = float(np.sum(y[:n_core][self.gas_species_mask]))
            V_gas = constants.R * self.T.value_si * n_gas / self.P.value_si if n_gas > 0 else 1.0
        self.V_gas = V_gas
        V_poly = self.V_poly
        self.V = self.V_gas + self.V_poly

        # 3. Concentrations (Vectorized)
        C_gas = self._scratch_C_gas
        C_poly = self._scratch_C_poly
        C_gas[:] = 0.0
        C_poly[:] = 0.0

        mask = self.gas_species_mask
        C_gas[mask] = np.maximum(0.0, y[mask]) / V_gas
        C_poly[~mask] = np.maximum(0.0, y[~mask]) / V_poly

        # Sync diagnostic concentrations
        self.core_species_concentrations[mask] = C_gas[mask]
        self.core_species_concentrations[~mask] = C_poly[~mask]

        # 4. Clear Accumulators
        dn_dt = self._scratch_dn_dt
        dn_dt[:] = 0.0

        self.core_reaction_rates[:] = 0.0
        self.edge_reaction_rates[:] = 0.0
        self.edge_species_rates[:] = 0.0
        self.core_species_consumption_rates[:] = 0.0
        self.core_species_production_rates[:] = 0.0
        self.network_leak_rates[:] = 0.0

        # 4.5 Recalculate P-dependent rates in the GAS phase
        if self.pdep_collision_reaction_indices.shape[0] != 0:
            T_val = self.T.value_si
            P_val = self.P.value_si

            y_gas = y[:n_core] * mask
            sum_y_gas = np.sum(y_gas)

            if sum_y_gas > 1e-20:
                for i in range(self.pdep_collision_reaction_indices.shape[0]):
                    j = self.pdep_collision_reaction_indices[i]

                    Peff = P_val * np.sum(self.collider_efficiencies[i] * y_gas / sum_y_gas)

                    self.kf[j] = self.pdep_collider_kinetics[i].get_rate_coefficient(T_val, Peff)
                    if self.Keq[j] != np.inf:
                        self.kb[j] = self.kf[j] / self.Keq[j]
            else:
                pass

        # 5. Standard Kinetics (Optimized Loop)
        for r_idx in range(ir.shape[0]):
            r0, r1, r2 = ir[r_idx, 0], ir[r_idx, 1], ir[r_idx, 2]

            if r0 >= n_core: continue

            # Determine phase from first reactant (all must match)
            is_r0_gas = mask[r0]

            # Check other reactants match r0's phase
            if r1 != -1 and (r1 >= n_core or mask[r1] != is_r0_gas):
                continue
            if r2 != -1 and (r2 >= n_core or mask[r2] != is_r0_gas):
                continue

            if is_r0_gas:
                phase = 'gas'
                V_rxn = V_gas
            else:
                phase = 'poly'
                V_rxn = V_poly

            # Product Checks
            p0, p1, p2 = ip[r_idx, 0], ip[r_idx, 1], ip[r_idx, 2]

            prods_phase_ok = True
            has_edge_prod = False
            has_any_prod = False

            # Unrolled Phase/Edge Check
            if p0 != -1:
                has_any_prod = True
                if p0 < n_core:
                    if mask[p0] != is_r0_gas:
                        prods_phase_ok = False
                else:
                    has_edge_prod = True

            if prods_phase_ok and p1 != -1:
                has_any_prod = True
                if p1 < n_core:
                    if mask[p1] != is_r0_gas:
                        prods_phase_ok = False
                else:
                    has_edge_prod = True

            if prods_phase_ok and p2 != -1:
                has_any_prod = True
                if p2 < n_core:
                    if mask[p2] != is_r0_gas:
                        prods_phase_ok = False
                else:
                    has_edge_prod = True

            if not prods_phase_ok:
                continue

            # Rate Calculation
            kf = self.kf[r_idx]
            kb = self.kb[r_idx]

            def _C(idx):
                if idx == -1:
                    return 1.0
                return C_gas[idx] if phase == 'gas' else C_poly[idx]

            if r1 == -1:
                rf = kf * _C(r0)
            elif r2 == -1:
                rf = kf * _C(r0) * _C(r1)
            else:
                rf = kf * _C(r0) * _C(r1) * _C(r2)

            # Reverse Rate (Unrolled, Hole-Safe)
            rr = 0.0
            if has_any_prod and not has_edge_prod:
                rr = kb
                if p0 != -1: rr *= _C(p0)
                if p1 != -1: rr *= _C(p1)
                if p2 != -1: rr *= _C(p2)

            rate = rf - rr

            # 1. Map Reactants to Polymer Pools
            p0_pool_idx = self.species_to_pool_indices[r0]
            p1_pool_idx = -1 if r1 == -1 else self.species_to_pool_indices[r1]

            # 2. Kinetic Hijack: Scale Rate by Polymer Site Density (mu1)
            if p0_pool_idx != -1 or p1_pool_idx != -1:
                if p0_pool_idx != -1:
                    # Reactant 0 is the polymer
                    conc_factor = max(0.0, y[self.polymer_pools[p0_pool_idx].mu_indices[1]]) / V_poly
                    rate = kf * conc_factor
                    if r1 != -1: rate *= _C(r1)
                    if r2 != -1: rate *= _C(r2)
                else:
                    # Reactant 1 is the polymer
                    conc_factor = max(0.0, y[self.polymer_pools[p1_pool_idx].mu_indices[1]]) / V_poly
                    rate = kf * _C(r0) * conc_factor
                    if r2 != -1: rate *= _C(r2)

            r_mol_s = rate * V_rxn

            # Update diagnostic rates for RMG's pruning logic
            if r_idx < n_rxn:
                self.core_reaction_rates[r_idx] = rate
            else:
                self.edge_reaction_rates[r_idx - n_rxn] = rate

            # 3. Apply Fluxes to Reactants (Redirect to Moments if needed)
            if p0_pool_idx != -1:
                dn_dt[self.polymer_pools[p0_pool_idx].mu_indices[1]] -= r_mol_s
            else:
                dn_dt[r0] -= r_mol_s

            if r1 != -1:
                if p1_pool_idx != -1:
                    dn_dt[self.polymer_pools[p1_pool_idx].mu_indices[1]] -= r_mol_s
                else:
                    dn_dt[r1] -= r_mol_s
            if r2 != -1:
                dn_dt[r2] -= r_mol_s

            # 4. Apply Fluxes to Products (Redirect to Moments if needed)
            for p_slot in range(3):
                p_idx_tmp = ip[r_idx, p_slot]
                if p_idx_tmp == -1: continue

                if p_idx_tmp < n_core:
                    prod_p_idx = self.species_to_pool_indices[p_idx_tmp]
                    if prod_p_idx != -1:
                        # Increment Monomer Units moment for the product pool
                        dn_dt[self.polymer_pools[prod_p_idx].mu_indices[1]] += r_mol_s
                    else:
                        dn_dt[p_idx_tmp] += r_mol_s
                else:
                    # Track edge species for model enlargement
                    self.edge_species_rates[p_idx_tmp - n_core] += rate

        # 6. Network Leaks
        for j in range(inet.shape[0]):
            s0 = inet[j, 0]
            if s0 == -1 or s0 >= n_core or not mask[s0]:
                continue

            s1, s2 = inet[j, 1], inet[j, 2]
            if (s1 != -1 and not mask[s1]) or (s2 != -1 and not mask[s2]):
                continue

            k = knet[j]
            rate = k * C_gas[s0]
            if s1 != -1: rate *= C_gas[s1]
            if s2 != -1: rate *= C_gas[s2]

            self.network_leak_rates[j] = rate
            dn_dt[s0] -= rate * V_gas

        # 7. Polymer Tail & Handshake
        for pool in self.polymer_pools:
            idx_mu0, idx_mu1, idx_mu2 = pool.mu_indices
            xs = pool.xs

            mu0_mol = max(0.0, y[idx_mu0])
            mu1_mol = max(0.0, y[idx_mu1])
            mu2_mol = max(0.0, y[idx_mu2])

            mu0 = mu0_mol / V_poly
            mu1 = mu1_mol / V_poly
            mu2 = mu2_mol / V_poly
            mu3 = _safe_mu3_from_mu012(mu0, mu1, mu2)

            dmu0_dt = 0.0
            dmu1_dt = 0.0
            dmu2_dt = 0.0
            small_src = dict()

            if pool.tail_kinetics:
                dmu0_dt, dmu1_dt, dmu2_dt, small_src = pool.tail_kinetics(
                    self.T.value_si, self.P.value_si, mu0, mu1, mu2, mu3)
            else:
                if pool.k_scission > 0:
                    dmu0_dt += pool.k_scission * mu1
                    if np.isfinite(mu3):
                         dmu2_dt += pool.k_scission * mu1 * (mu3 / max(mu1, SMALL_EPS) - mu2)

                if pool.k_unzip > 0:
                    r_events = pool.k_unzip * mu0
                    dmu1_dt -= r_events
                    dmu2_dt -= pool.k_unzip * (2.0 * mu1 - mu0)
                    if pool.monomer_poly_index is not None:
                        small_src[pool.monomer_poly_index] = r_events

            # Hybrid Handshake
            tail_mean = mu1 / mu0 if mu0 > SMALL_EPS else 0.0
            valid_tail = (mu0 > TAIL_CONC_MIN) and (tail_mean > xs + 1e-9)

            if valid_tail and pool.k_unzip > 0:
                explicit_xs_idx = pool.explicit_dp_to_species_index.get(xs, None)
                if explicit_xs_idx is not None:
                    params = _gamma_params_from_mu012(mu0, mu1, mu2)
                    p_cond = 0.0

                    if params:
                        k_shape, theta = params
                        p_cond = _gamma_prob_conditional_hybrid(xs + 1, xs, k_shape, theta)
                    else:
                        if tail_mean <= xs + 1.0:
                            p_cond = 0.0
                        elif tail_mean >= xs + 2.0:
                            p_cond = 0.0
                        else:
                            p_cond = 1.0 - abs(tail_mean - (xs + 1.5)) / 0.5

                    p_cond = min(1.0, max(0.0, p_cond))

                    # Population (conc) at boundary
                    N_boundary = mu0 * p_cond

                    # Flux Clamp
                    N_boundary = min(N_boundary, mu0)
                    if xs > 0:
                        N_boundary = min(N_boundary, mu1 / xs)
                        N_boundary = min(N_boundary, mu2 / (xs * xs))

                    F = pool.k_unzip * N_boundary

                    if F > 0.0:
                        dn_dt[explicit_xs_idx] += F * V_poly
                        dmu0_dt -= F
                        dmu1_dt -= xs * F
                        dmu2_dt -= (xs * xs) * F

            dn_dt[idx_mu0] += dmu0_dt * V_poly
            dn_dt[idx_mu1] += dmu1_dt * V_poly
            dn_dt[idx_mu2] += dmu2_dt * V_poly

            for spc_idx, r_vol in small_src.items():
                if 0 <= spc_idx < n_core:
                    dn_dt[spc_idx] += r_vol * V_poly

        # 8. Mass Transfer
        for mt in self.mass_transfer:
            ig = mt.gas_index
            ipoly = mt.poly_index

            Cg = C_gas[ig]
            Cp = C_poly[ipoly]

            J = mt.kLa * (Cp - mt.K * Cg)
            dn = J * V_poly

            dn_dt[ig] += dn
            dn_dt[ipoly] -= dn

        # 9. Constants
        if self.const_spc_indices is not None:
            for i in self.const_spc_indices:
                dn_dt[i] = 0.0

        # 10. Diagnostics
        self.core_species_rates[mask] = dn_dt[mask] / V_gas
        self.core_species_rates[~mask] = dn_dt[~mask] / V_poly

        return (dn_dt - dydt), 1
