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
# Attribution trust floor multiplier (item #14a, spec 2026-06-11 §2/§3):
# trust = max(SMALL_EPS, ATTRIBUTION_TRUST_K * atol_mu0). K = 100 buys ~two
# decades of separation between the integrator's own error budget (the atol
# on the pool's mu0 slot) and what the spawn-gate snapshot will TRUST as an
# E[n] denominator — a margin decision on record (like the thermo bounds),
# not a default. Tolerance-anchored on purpose: any absolute mole constant
# is wrong across deck scales (a ~10 mg TGA sample at Mn ~1e4-1e5 g/mol puts
# a whole physical deck's mu0 at ~1e-10-1e-9 mol). SMALL_EPS itself is
# untouched and keeps its job as the solver's REALIZABILITY clamp; ONLY the
# gate snapshot's attribution distrusts the band between them — the gate
# deliberately distrusts what the solver still legally computes (two
# constants, two jobs).
ATTRIBUTION_TRUST_K = 100.0
LN_EXP_OVERFLOW_GUARD = 700.0
TAIL_CONC_MIN = 1e-9  # Minimum concentration (mol/m^3) to actuate handshake

# Pool moment-flux archetypes. Mirror of rmgpy.polymer.PolymerFluxArchetype
# (not imported to avoid a solver->polymer module cycle); equality is pinned
# by test_flux_archetype_constants_match_enum.
FLUX_NONE = 0
FLUX_SAME_POOL = 1
FLUX_MIGRATION = 2
FLUX_SCISSION_FRAGMENT = 3
FLUX_UNRESOLVED = 4
FLUX_DISCRETE_CHIP = 5


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

    # Monomer (repeat-unit) molecular weight [g/mol]. Consumed by
    # spawn_gate_flux_snapshot (E[n]*MW event-mass calibration, spec
    # 2026-06-10-mass-flux-spawn-gate-design.md §3). Default 0.0 zeroes the
    # pool's snapshot contribution (honest degradation: the spawn gate
    # defers; nothing is fabricated).
    monomer_mw_g_mol: float = 0.0

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
    """Log-Lagrange closure μ3 = μ0·(μ2/μ1)³, with a realizability guard.

    A valid chain-length distribution (k ≥ 1) always satisfies μ1 ≥ μ0. If
    solver noise or a bad source term pushes the state out of the realizable
    cone, the (μ2/μ1)³ extrapolation can explode into a DASSL singularity, so
    we return 0.0 (no closure contribution) rather than amplifying garbage.
    """
    if mu0 <= SMALL_EPS or mu1 <= SMALL_EPS or mu2 <= SMALL_EPS:
        return 0.0
    if mu1 < mu0:  # unrealizable: more chains than monomer units
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
# THERMO REFERENCE-STATE TRIPWIRE (spec 2026-06-11)
# ======================================================================================

# Standard-state constants for the unpaired reference-state magnitude U.
# P0 is the NASA-polynomial gas standard pressure; C0 = 1 mol/m^3 is the
# standard concentration that makes log10(P0/(R*T*C0)) dimensionless -- at
# 1000 K the term is 1.08 decades. Omit C0 in a re-derivation in other units
# and the formula merely LOOKS wrong (spec §5.1).
REFERENCE_STATE_P0 = 1.0e5      # Pa
REFERENCE_STATE_C0 = 1.0        # mol/m^3
# Bounds (spec §6 -- chemistry decisions on record, not defaults): census
# above the measured benign ceiling (0.33 decades + small-margin items);
# refuse with a wide margin below the pathological floor (~10 decades),
# grounded in bimodality plus the §5.2 sign argument. The 0.5-3.0 band is
# deliberately wide and deliberately EMPTY on today's measurements; anything
# that ever lands in it is news AND the trigger to compute rotation into U
# properly (spec §5.2 escalation, on record in multi_pool_design.md §5.2).
REFERENCE_STATE_CENSUS_DECADES = 0.5
REFERENCE_STATE_REFUSE_DECADES = 3.0
# Chain-scale window slack on top of the largest pool monomer MW. ONE window
# (max pool monomer MW + slack) serves BOTH the physically-melt class
# tag-branch gate (spec §5.1, C3-amended) and the §5.3 provenance
# counterparty predicate -- one definition, two uses, so the class and the
# sensor cannot drift apart.
REFERENCE_STATE_MW_SLACK_G_MOL = 10.0


def _sackur_tetrode_s_trans(mw_kg_mol, T):
    """Exact Sackur-Tetrode translational entropy [J/mol/K] of an ideal gas
    of molar mass ``mw_kg_mol`` [kg/mol] at ``T`` [K] and P = REFERENCE_STATE_P0."""
    m = mw_kg_mol / constants.Na
    return constants.R * (
        math.log((2.0 * constants.pi * m * constants.kB * T / constants.h ** 2) ** 1.5
                 * constants.kB * T / REFERENCE_STATE_P0) + 2.5)


def _unpaired_reference_decades(reactant_melt_mws, product_melt_mws, T):
    """Unpaired reference-state magnitude U [decades] for one reversible
    reaction (spec 2026-06-11 §5.1):

        U = |sum_prod,melt S_trans - sum_react,melt S_trans| / (R ln10)
            + |dn_melt| * log10(P0 / (R*T*C0))

    Inputs are the molar masses [kg/mol] of the PHYSICALLY-MELT participants
    (condensed-mask OR is_polymer_proxy-tagged with chain-scale MW -- spec
    §5.1, C3-amended) on each side, multiplicity
    included; gas participants are excluded because their gas reference
    state is simply CORRECT. Paired same-mass chains cancel in the signed
    sum exactly as they do in the real thermo -- no pairing optimization is
    needed. Translational-only BY DESIGN: the measured U distribution is
    bimodal with a four-decade gap (benign <= 0.33, pathological >= ~10), so
    the omitted bounded rotational term (+1.6-3.1 decades, same sign where
    it matters) cannot flip a decision against the §6 thresholds. Escalation
    trigger (spec §5.2): a census entry in the 1-5 decade band means rotation
    must be computed into U properly."""
    s_sum = 0.0
    for mw in product_melt_mws:
        s_sum += _sackur_tetrode_s_trans(mw, T)
    for mw in reactant_melt_mws:
        s_sum -= _sackur_tetrode_s_trans(mw, T)
    dn_melt = len(product_melt_mws) - len(reactant_melt_mws)
    return (abs(s_sum) / (constants.R * math.log(10.0))
            + abs(dn_melt) * math.log10(
                REFERENCE_STATE_P0 / (constants.R * T * REFERENCE_STATE_C0)))


def _assert_chain_scale_melt_member(label, mw_kg_mol, gas_classified, window_kg_mol):
    """Cannot-happen leak guard (spec §5.1, C3 amendment 2026-06-11), called
    for every species entering the melt sum. A gas-classified species can be
    physically-melt ONLY via the tag branch, whose chain-scale MW conjunct
    (MW >= max pool monomer + slack) is part of the CLASS DEFINITION -- so a
    gas-classified member below the window cannot reach the sum through the
    amended gate. A tagged below-window species being EXCLUDED by the gate is
    expected and silent (the family.py:1657 over-tagging fingerprint, H2 on
    every proxy-touching reaction); the raise here fires only if such a
    species REACHES the sum (a future refactor recomputing membership without
    the conjunct), converting the mistake into a loud CLASSIFICATION error
    instead of a silently-large U and a misattributed reference-state
    refusal. Condensed-branch members (gas_classified=False; pool-configured
    by input) are exempt from the gate and from this guard."""
    if gas_classified and mw_kg_mol < window_kg_mol:
        raise ValueError(
            "THERMO REFERENCE-STATE TRIPWIRE (classification leak): the "
            "is_polymer_proxy tag includes a non-chain species in the melt "
            f"sum ({label}, MW = {mw_kg_mol * 1000.0:.2f} g/mol < chain-scale "
            f"window {window_kg_mol * 1000.0:.2f} g/mol); physically-melt "
            "class definition violated -- this is a classification leak, NOT "
            "a thermo problem; do not respond by touching reference states. "
            "See the proxy-tag propagation chain (family.py) and the "
            "invariant section of docs/multi_pool_design.md.")


def _thermo_provenance(spc):
    """Classify a species' thermo source from its comment: 'library', 'gav',
    or None (no thermo / unrecognized -- never warned on). Substrings pinned
    against rmgpy/data/thermo.py: 'Thermo library: <label>' (:1818) and
    'Thermo group additivity estimation: ...' (:232/:2237; gav_keywords
    :2845). 'Thermo library' is checked FIRST: an HBI radical estimated from
    a library parent ('Thermo library: X + radical(Y)') classifies 'library',
    which is correct -- a library-parent HBI resolves through the same parent
    and PRESERVES the structural cancellation (spec §2)."""
    comment = getattr(getattr(spc, "thermo", None), "comment", "") or ""
    if "Thermo library" in comment:
        return "library"
    if "group additivity" in comment:
        return "gav"
    return None


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
        prospective_gas_mask: Optional[np.ndarray] = None,
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
        debug_check_realizability: bool = False,
        allow_unpaired_reference_state: bool = False,
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
        # Item 17 (spec 2026-06-12 SS3(a)/(b)): stage-1 seed for the
        # prospective mask over chain(core, edge), normally computed by
        # to_solver_object with the SAME config-keyed classifier as
        # gas_species_mask (one classifier, longer list). Stored as a seed;
        # initialize_model rebuilds self.prospective_gas_mask on EVERY call.
        # RIDER R3: prospective_gas_mask is GATE-INPUT ONLY -- nothing but
        # the residual product gates (and riders R1/R2, which verify and
        # report them) may read it; gas_species_mask stays the sole source
        # of truth for every other phase behavior (concentrations, volumes,
        # validation, the thermo tripwire's melt class, the artifact mask
        # harvest).
        self._prospective_gas_mask_seed = prospective_gas_mask
        self.prospective_gas_mask = None
        self.const_spc_names = const_spc_names or []
        self.const_spc_indices = None
        self.sens_conditions = sens_conditions

        self.pdep_collision_reaction_indices = pdep_collision_reaction_indices if pdep_collision_reaction_indices is not None else np.array([], int)
        self.pdep_collider_kinetics = pdep_collider_kinetics if pdep_collider_kinetics is not None else []
        self.collider_efficiencies = collider_efficiencies if collider_efficiencies is not None else np.array([[]], float)
        self.initial_polymer_moments = initial_polymer_moments or {}

        # Opt-in moment-realizability diagnostic (off by default; zero cost when
        # off). When on, the residual logs ONCE per pool if its moment state
        # leaves the realizable cone (μ1≥μ0≥0, μ0·μ2≥μ1²) — a fast way to localize
        # a bad moment source term. It only logs (never raises): the μ3 closure
        # already degrades gracefully on out-of-cone states.
        self.debug_check_realizability = bool(debug_check_realizability)
        self._realizability_warned = set()

        # Thermo reference-state tripwire (spec 2026-06-11 §7). The override
        # is the deck author's CONSCIOUS assertion that their thermo handles
        # unpaired reference states -- the eventual per-deck switch-on point
        # for full melt consistency. Census/diagnostics are repopulated by
        # _reference_state_tripwire on every initialize_model rebuild:
        # reference_state_census = [(str(rxn), U)] for U > census bound;
        # reference_state_max_decades = max U over reversible core reactions.
        self.allow_unpaired_reference_state = bool(allow_unpaired_reference_state)
        self.reference_state_census = []
        self.reference_state_max_decades = 0.0
        self.initial_explicit_species = initial_explicit_species or {}
        self.polymer_species_labels = set(polymer_species_labels) if polymer_species_labels else set()

        self.jacobian_matrix = None

        self._scratch_C_gas = None
        self._scratch_C_poly = None
        self._scratch_dn_dt = None

        self.pool_mu0_indices = np.full(len(self.polymer_pools), -1, dtype=np.int32)
        self.pool_mu1_indices = np.full(len(polymer_pools), -1, dtype=np.int32)

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

    def _reference_state_tripwire(self, core_species, core_reactions):
        """Build-time thermo reference-state tripwire (spec 2026-06-11).

        Invariant (docs/multi_pool_design.md §5.2): per species, per reaction
        side, the thermo reference state must match the phase residence.
        Today every species carries the gas reference UNIFORMLY, so
        boundary-crossing reactions between same-mass chains cancel exactly
        -- and the cancellation is STRUCTURAL (HBI saturates the radical onto
        the same GAV parent the proxy runs on). The danger is the PARTIAL
        fix: a melt correction applied to the condensed set alone injects up
        to ~11.6 decades into every boundary-crossing Keq. This pass
        measures, per reversible CORE reaction, the unpaired reference-state
        magnitude U over the physically-melt participants (condensed
        gas_species_mask OR is_polymer_proxy-tagged with chain-scale MW --
        spec §5.1, C3-amended), logs a census above
        REFERENCE_STATE_CENSUS_DECADES, refuses above
        REFERENCE_STATE_REFUSE_DECADES (unless
        allow_unpaired_reference_state -- the census still logs), and warns
        on mixed library-vs-GAV thermo provenance among chain-scale
        counterparties (the spec-§5.3 decoupling fingerprint that U, being
        MW-only, is structurally blind to). Cost: ~|core reactions|
        Sackur-Tetrode evaluations per rebuild -- nil.
        """
        cdef int n_core = self.num_core_species
        T = self.T.value_si
        self.reference_state_census = []
        self.reference_state_max_decades = 0.0

        mask = self.gas_species_mask

        # ONE chain-scale window (spec §5.1 C3 amendment + §5.3): largest
        # configured pool monomer MW + slack, in kg/mol. Shared by BOTH the
        # physically-melt class gate (tag branch, below) and the provenance
        # counterparty predicate (i) -- one definition, two uses, so the
        # class and the sensor cannot drift apart. Predicate (ii) (sharing
        # the proxy's saturated parent) was REJECTED as the cheap-at-init
        # test: it needs graph saturation + isomorphism per gas species per
        # rebuild; (i) is a float compare on data already in hand.
        chain_window_kg = (max((float(p.monomer_mw_g_mol) for p in self.polymer_pools),
                               default=0.0)
                           + REFERENCE_STATE_MW_SLACK_G_MOL) / 1000.0

        is_melt = [False] * n_core
        mws = [0.0] * n_core
        for i in range(n_core):
            spc = core_species[i]
            # Input contract: consumer-world species carry no structure
            # (molecule == []); MW arrives via the species-level
            # molecular_weight quantity, populated by the runner from the
            # artifact's composition. Prefer it (value_si is kg/molecule);
            # the structure branch below is a defensive fallback (normally
            # unreachable: the molecular_weight property lazily derives from
            # the structure when one exists -- species.py:278-282 -- so the
            # quantity path above short-circuits it).
            mol_list = getattr(spc, "molecule", None)
            mw_q = getattr(spc, "molecular_weight", None)
            if mw_q is not None and mw_q.value_si > 0.0:
                mws[i] = mw_q.value_si * constants.Na
            elif mol_list:
                mws[i] = mol_list[0].get_molecular_weight()
            else:
                mws[i] = 0.0
            # Physically-melt CLASS (spec §5.1, C3-AMENDED): the condensed
            # branch (pool-configured by input) unconditionally; the tag
            # branch only at chain-scale MW. The MW conjunct is part of the
            # class DEFINITION, not bolted onto the provenance set:
            # family.py:1657 blanket-tags every structure of a
            # proxy-touching reaction (including H2), so a raw tag-read
            # would be correct only by spawn-pass ordering -- the structural
            # gate cannot be broken by a lifecycle reordering the way a
            # tag-read can. A stale tag on a below-window species simply
            # FAILS the conjunct and is excluded: expected and silent (its
            # gas reference state is CORRECT).
            is_melt[i] = ((not mask[i])
                          or (bool(getattr(spc, "is_polymer_proxy", False))
                              and mws[i] >= chain_window_kg))

        offenders = []
        for rxn in core_reactions:
            if not getattr(rxn, "reversible", False):
                continue
            j = self.reaction_index[rxn]
            r_idx = [int(k) for k in self.reactant_indices[j] if 0 <= k < n_core]
            p_idx = [int(k) for k in self.product_indices[j] if 0 <= k < n_core]
            melt_r = [k for k in r_idx if is_melt[k]]
            melt_p = [k for k in p_idx if is_melt[k]]
            if not melt_r and not melt_p:
                continue  # all-gas reaction: gas reference uniformly correct

            # Leak self-assertion (spec §5.1 C3 amendment): cannot-happen
            # guard on every species entering the melt sum. A gas-classified
            # member can only be here via the tag branch, whose MW conjunct
            # is enforced in the gate above; if a below-window species ever
            # reaches this point (a future refactor recomputing membership
            # without the conjunct), raise the CLASSIFICATION error loudly
            # instead of computing a large U and misattributing it to thermo.
            for k in melt_r + melt_p:
                # Log-domain guard (input contract): a melt participant with
                # no resolvable molar mass would send mw <= 0 into
                # _sackur_tetrode_s_trans's math.log. A melt chain with no
                # structure is the tripwire's MAIN case (consumer world), not
                # an edge case -- never skip it silently, and never let a raw
                # 'math domain error' escape.
                if mws[k] <= 0.0:
                    raise ValueError(
                        "THERMO REFERENCE-STATE TRIPWIRE (input contract): no "
                        f"molar mass available for melt participant "
                        f"'{core_species[k].label}' (molecule list empty and "
                        "species-level molecular_weight unset); consumer-world "
                        "species must carry molecular_weight populated from "
                        "the artifact's composition/mw fields -- this is an "
                        "input-contract violation, NOT a thermo problem; do "
                        "not respond by touching reference states.")
                _assert_chain_scale_melt_member(
                    core_species[k].label, mws[k], bool(mask[k]),
                    chain_window_kg)

            u = _unpaired_reference_decades(
                [mws[k] for k in melt_r], [mws[k] for k in melt_p], T)
            if u > self.reference_state_max_decades:
                self.reference_state_max_decades = u
            if u > REFERENCE_STATE_CENSUS_DECADES:
                self.reference_state_census.append((str(rxn), u))
                logging.warning(
                    "THERMO REFERENCE-STATE CENSUS: reaction %s: U = %.2f "
                    "decades at T = %.1f K (above the %.1f-decade census "
                    "bound; see the invariant section of "
                    "docs/multi_pool_design.md).",
                    rxn, u, T, REFERENCE_STATE_CENSUS_DECADES)
            if u > REFERENCE_STATE_REFUSE_DECADES:
                offenders.append((rxn, u))

            # Provenance sensor (spec §5.3): counterparty set = melt
            # participants + gas participants inside the shared chain-scale
            # MW window of a melt participant. Deliberately NARROW: small
            # co-reactants (H2, CH4, every abstraction partner) legitimately
            # take library thermo while chains take GAV -- sweeping them in
            # would warn on every healthy deck (alarm fatigue re-arming the
            # exact landmine this sensor guards).
            melt_set = set(melt_r) | set(melt_p)
            counterparties = set(melt_set)
            for k in set(r_idx) | set(p_idx):
                if k in melt_set:
                    continue
                if any(abs(mws[k] - mws[m]) <= chain_window_kg for m in melt_set):
                    counterparties.add(k)
            examples = {}
            for k in sorted(counterparties):
                prov = _thermo_provenance(core_species[k])
                if prov is not None and prov not in examples:
                    examples[prov] = core_species[k].label
            if "library" in examples and "gav" in examples:
                logging.warning(
                    "THERMO REFERENCE-STATE PROVENANCE: mixed thermo "
                    "provenance among chain-scale counterparties in reaction "
                    "%s (library: %s; group additivity: %s) -- the structural "
                    "cancellation that keeps gas-reference thermo safe on "
                    "melt chains may be broken for this pair; see the "
                    "invariant section of docs/multi_pool_design.md.",
                    rxn, examples["library"], examples["gav"])

        if offenders:
            if self.allow_unpaired_reference_state:
                logging.warning(
                    "THERMO REFERENCE-STATE TRIPWIRE: "
                    "allow_unpaired_reference_state=True: bypassing the "
                    "reference-state refusal for %d reaction(s); max U = "
                    "%.2f decades. The deck author asserts the thermo "
                    "handles the melt reference state.",
                    len(offenders), self.reference_state_max_decades)
            else:
                lines = "\n".join(
                    f"  {rxn}: U = {u:.2f} decades" for rxn, u in offenders)
                raise ValueError(
                    "THERMO REFERENCE-STATE TRIPWIRE: this deck has "
                    f"{len(offenders)} reversible core reaction(s) with an "
                    "unpaired reference-state term (U > "
                    f"{REFERENCE_STATE_REFUSE_DECADES} decades) at "
                    f"T = {T:.1f} K:\n{lines}\n"
                    "Do NOT fix this by applying a melt correction to the "
                    "condensed set alone -- that injects the full mismatch "
                    "into every boundary-crossing Keq (decapping kb by up to "
                    "~12 orders of magnitude). See the thermo reference-state "
                    "invariant section of docs/multi_pool_design.md. If the "
                    "deck's thermo genuinely handles the melt reference "
                    "state, set allow_unpaired_reference_state=True on the "
                    "reactor.")

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

    def _apply_pool_phase_overrides(self, mask_arr, species_list,
                                    record_indices):
        """Stage 2 of the two-stage phase classifier (item 17, spec
        2026-06-12 SS3(a)): the per-pool config-label override pass, factored
        so ONE code path serves both masks. Run with
        (gas_species_mask, core_species, record_indices=True) -- behavior
        bit-identical to the historical inline pass, including the index
        bookkeeping (species_to_pool_indices / is_pool_proxy /
        pool_mu1_indices / pool_mu0_indices and the mu1 fallback) -- and
        with (prospective_gas_mask, chain(core, edge),
        record_indices=False), where ONLY mask writes happen (the index
        arrays are core-sized; an edge match must never write them)."""
        n_core = self.num_core_species
        n_total = len(species_list)
        for pool_i, pool in enumerate(self.polymer_pools):
            for i in range(n_total):
                base_label = species_list[i].label.partition('(')[0]
                if base_label == pool.label:
                    mask_arr[i] = False
                    if record_indices:
                        self.species_to_pool_indices[i] = pool_i
                        self.is_pool_proxy[i] = 1
                    break
            mu1_target_label = f"{pool.label}_mu1"
            for i in range(n_total):
                # Handle RMG renaming: "PS_mu1(2)" -> "PS_mu1"
                base_label = species_list[i].label.partition('(')[0]
                if base_label == mu1_target_label:
                    if record_indices:
                        self.pool_mu1_indices[pool_i] = i
                    mask_arr[i] = False  # Ensure it's poly phase

            mu0_target_label = f"{pool.label}_mu0"
            for i in range(n_total):
                base_label = species_list[i].label.partition('(')[0]
                if base_label == mu0_target_label:
                    if record_indices:
                        self.pool_mu0_indices[pool_i] = i
                    mask_arr[i] = False

            # map explicit oligomers
            for dp, idx in pool.explicit_dp_to_species_index.items():
                if 0 <= idx < n_core:
                    if record_indices:
                        self.species_to_pool_indices[idx] = pool_i
                    mask_arr[idx] = False

            # map the moment indices
            for idx in pool.mu_indices:
                if 0 <= idx < n_core:
                    if record_indices:
                        self.species_to_pool_indices[idx] = pool_i
                    mask_arr[idx] = False

            # map monomer-in-poly if used
            if pool.monomer_poly_index is not None and 0 <= pool.monomer_poly_index < n_core:
                if record_indices:
                    self.species_to_pool_indices[pool.monomer_poly_index] = pool_i
                mask_arr[pool.monomer_poly_index] = False

            if record_indices and self.pool_mu1_indices[pool_i] == -1:
                # Fallback: Try to use the config index if label lookup failed
                # (This catches cases where renaming didn't happen as expected)
                cfg_idx = pool.mu_indices[1]
                if 0 <= cfg_idx < n_core:
                    self.pool_mu1_indices[pool_i] = cfg_idx
                else:
                    print(f"WARNING: Could not locate mu1 species for pool {pool.label}. Polymer chemistry will be inert.")

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
        cdef int n_rxn = len(core_reactions) + len(edge_reactions)
        # Flag end-group (terminal) reactions so their proxy rate scales by
        # chain-end density (mu0) instead of monomer-unit density (mu1). Iterate
        # in the SAME order that builds kf/kb and the ir/ip arrays
        # (generate_rate_coefficients below) so the index matches r_idx.
        self.is_end_group_reaction = np.zeros(n_rxn, dtype=np.int8)
        for i, rxn in enumerate(itertools.chain(core_reactions, edge_reactions)):
            if getattr(rxn, "is_end_group_reaction", False):
                self.is_end_group_reaction[i] = 1

        # Per-reaction pool moment-flux archetype (spec 2026-06-09). Same
        # chain(core, edge) order as is_end_group_reaction so indices match
        # r_idx in the residual.
        self.reaction_flux_archetype = np.zeros(n_rxn, dtype=np.int8)
        self.reaction_src_pool = np.full(n_rxn, -1, dtype=np.int32)
        self.reaction_dst_pool = np.full(n_rxn, -1, dtype=np.int32)
        # Stamped chip repeat-unit counts (spec 2026-06-10 §4.3); same
        # chain(core, edge) order so the index matches r_idx in the residual.
        self.reaction_chip_units = np.zeros(n_rxn, dtype=np.int32)
        for i, rxn in enumerate(itertools.chain(core_reactions, edge_reactions)):
            self.reaction_flux_archetype[i] = int(getattr(rxn, "polymer_flux_archetype", 0))
            self.reaction_chip_units[i] = int(getattr(rxn, "polymer_chip_units", 0))
        # Item 17 (spec 2026-06-12 SS3(e)): generation-time stamps, captured
        # AFTER the stamp-read loop and BEFORE the init demotion pass
        # (:NONE->UNRESOLVED + unresolvable stamped shapes) mutates the
        # array in place. For reactions flip-demoted at GENERATION time
        # (demote_flipped_polymer_archetype mutates the object itself),
        # "pre-demotion" deliberately means pre-SOLVER-demotion: the
        # captured value is legitimately UNRESOLVED and the census reports
        # it as such.
        self.reaction_pre_demotion_archetype = self.reaction_flux_archetype.copy()

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

        # Build mapping: core species index -> polymer pool index (or -1)
        self.species_to_pool_indices = np.full(n_core, -1, dtype=np.int32)
        self.is_pool_proxy = np.zeros(n_core, dtype=np.int8)
        self.pool_mu1_indices.fill(-1)
        self._apply_pool_phase_overrides(self.gas_species_mask, core_species,
                                         record_indices=True)

        # --- prospective_gas_mask (item 17, spec 2026-06-12 SS3(a)) -------
        # A SECOND array over chain(core, edge), never a resize of
        # gas_species_mask (SS3(b): the core size is load-bearing -- six hard
        # raises, ~10 length-coupled consumers and two silent fallbacks key
        # on it). Stage 1 = the blueprint classifier seed (combined-list
        # get_gas_mask, passed by to_solver_object); stage 2 = the SAME
        # per-pool override pass run on the combined list (label-keyed
        # writes scan core+edge; index-keyed writes land identically in the
        # core prefix). NEVER polymer-identity shortcuts: probed INCONSISTENT
        # with the post-promotion mask (spec SS2).
        n_edge_spc = len(edge_species) if edge_species is not None else 0
        combined_species = list(core_species) + (list(edge_species)
                                                 if edge_species else [])
        seed = self._prospective_gas_mask_seed
        if seed is not None and len(seed) == n_core + n_edge_spc:
            self.prospective_gas_mask = np.asarray(seed, dtype=bool).copy()
        else:
            if seed is not None:
                # Engine-reuse path (polymer_input.py:176-180): a
                # constructor-era seed can be sized for a previous build's
                # edge. Loud fallback, never a crash; R1 below still proves
                # the core prefix every build. At HEAD config the fallback
                # is verdict-identical to a fresh stage-1 seed: every
                # species stage 1 classifies condensed is input-configured
                # into core (premise A3), and stage 2 re-applies pool labels
                # on the CURRENT combined list.
                logging.warning(
                    "PROSPECTIVE-MASK SEED STALE: stage-1 seed length %d != "
                    "n_core + n_edge = %d; rebuilding prospective mask from "
                    "the fallback (core mask + edge defaults GAS).",
                    len(seed), n_core + n_edge_spc)
            # Fallback path (direct solver construction without the
            # blueprint seed -- tests, polymer_moments_runner). Stated
            # premise, not an implication (amendment A3): edge species
            # default prospectively-GAS here, which would diverge from
            # production stage 1 only for an edge species stage 1 classifies
            # condensed -- unreachable today, because get_gas_mask's registry
            # draws only from input-configured objects and all of those live
            # in core_species (the mass-transfer raises index the CORE-sized
            # mask). Fixtures express a prospectively-condensed EDGE species
            # through configured pool labels (stage 2) -- exactly how
            # production daughters become condensed under item 16.
            self.prospective_gas_mask = np.concatenate([
                np.asarray(self.gas_species_mask, dtype=bool).copy(),
                np.ones(n_edge_spc, dtype=bool)])
        self._apply_pool_phase_overrides(self.prospective_gas_mask,
                                         combined_species,
                                         record_indices=False)

        # RIDER R1 -- core-prefix parity tripwire (spec SS3(d)). The
        # architecture's central claim ("the prospective mask is the real
        # mask evaluated early") made self-verifying: the gates may use
        # prospective[p] uniformly for all p precisely BECAUSE this prefix
        # is proven equal, every build. Raise, never warn. Re-runs on every
        # rebuild, so spawned daughters and config changes are re-checked.
        if not np.array_equal(
                np.asarray(self.prospective_gas_mask[:n_core], dtype=bool),
                np.asarray(self.gas_species_mask, dtype=bool)):
            diverging = [
                f"index {i} ({core_species[i].label}): core mask says "
                f"{'GAS' if self.gas_species_mask[i] else 'CONDENSED'}, "
                f"prospective says "
                f"{'GAS' if self.prospective_gas_mask[i] else 'CONDENSED'}"
                for i in range(n_core)
                if bool(self.prospective_gas_mask[i])
                != bool(self.gas_species_mask[i])]
            raise ValueError(
                "PROSPECTIVE-MASK TRIPWIRE: prospective_gas_mask core prefix "
                "diverges from gas_species_mask: " + "; ".join(diverging)
                + ". First suspect: duplicate-label fallback disablement in "
                "the combined get_gas_mask call (a duplicate label present "
                "only in the edge disables the label fallback for the "
                "combined list -- polymer_input.py:609-612).")

        # Resolve per-reaction source/target pools from species indices (no
        # label matching), and remap unstamped proxy-touching reactions
        # (archetype NONE, e.g. unpickled) to UNRESOLVED so legacy mu1 flux
        # applies instead of silently dropping pool moment flux.
        ir_arr = self.reactant_indices
        ip_arr = self.product_indices
        n_unstamped = 0
        n_unresolvable = 0
        for i in range(n_rxn):
            src = -1
            for slot in range(3):
                ridx = ir_arr[i, slot]
                if ridx != -1 and ridx < n_core and self.species_to_pool_indices[ridx] != -1:
                    src = self.species_to_pool_indices[ridx]
                    break
            dst = -1
            for slot in range(3):
                pidx = ip_arr[i, slot]
                if pidx != -1 and pidx < n_core:
                    pool_j = self.species_to_pool_indices[pidx]
                    if pool_j != -1:
                        if pool_j != src:
                            dst = pool_j  # prefer the cross-pool product
                            break
                        # pool_j == src here: keep the same-pool fold-back only
                        # as a fallback while no cross-pool product was found.
                        if dst == -1:
                            dst = pool_j
            self.reaction_src_pool[i] = src
            self.reaction_dst_pool[i] = dst
            if self.reaction_flux_archetype[i] == FLUX_NONE and (src != -1 or dst != -1):
                self.reaction_flux_archetype[i] = FLUX_UNRESOLVED
                n_unstamped += 1
            if ((self.reaction_flux_archetype[i] in (FLUX_MIGRATION, FLUX_SCISSION_FRAGMENT)
                    and (src == -1 or dst == -1))
                    or (self.reaction_flux_archetype[i] == FLUX_DISCRETE_CHIP
                        and src == -1)):
                # A stamped archetype needs its pool(s) resolved in the
                # solver: MIGRATION/SCISSION_FRAGMENT need BOTH src and dst
                # (e.g. scission daughters are registered as core Polymer
                # species but have no pool config yet); DISCRETE_CHIP needs
                # only src (no dst: complement folds back to the same pool
                # and the chip is a plain gas species). Demote to the legacy
                # mu1-only transfer so the parent drain is never silently
                # zeroed (mass would otherwise be duplicated).
                # src == dst is deliberately NOT demoted: for that shape
                # (fold-back proxy product + non-pool daughter) the dispatch
                # skip and the legacy transfer produce identical pool flux
                # (reactant -r and fold-back +r cancel on the same mu1), so
                # demotion would change nothing.
                self.reaction_flux_archetype[i] = FLUX_UNRESOLVED
                n_unresolvable += 1
        if n_unresolvable:
            logging.warning(
                "%d reactions stamped MIGRATION/SCISSION_FRAGMENT/DISCRETE_CHIP "
                "could not resolve their solver pool(s); demoted to legacy "
                "mu1-only moment flux (UNRESOLVED).", n_unresolvable)
        if n_unstamped:
            logging.warning(
                "%d proxy-touching reactions arrived without a polymer_flux_archetype "
                "stamp; applying legacy mu1-only pool moment flux for them.",
                n_unstamped)

        # Enforce the moment-isolation invariant and pool/mass-transfer index
        # sanity now that gas_species_mask and the reaction index tables are
        # populated. Moments must evolve only via the tail block, never through
        # generic reaction stoichiometry.
        # RIDER R2 static half (spec SS3(e)): once per rebuild by
        # construction -- the enumeration runs exactly once per
        # initialize_model; no keying set needed.
        self._static_phase_gate_census(core_species, core_reactions)

        self.validate_configuration()

        # Thermo reference-state tripwire (spec 2026-06-11 §7): runs AFTER
        # the archetype demotion pass and validate_configuration (masks,
        # pool membership and archetypes are final here) and BEFORE
        # generate_rate_coefficients computes any kb from Keq -- refusal
        # gates SOLVING. Re-runs on every rebuild, so spawned daughters are
        # checked automatically.
        self._reference_state_tripwire(core_species, core_reactions)

        self._scratch_C_gas = np.zeros(n_core, float)
        self._scratch_C_poly = np.zeros(n_core, float)
        self._scratch_dn_dt = np.zeros(n_core, float)
        self._scratch_proxy_activity = np.zeros(n_core, float)

        # RIDER R2 dynamic half (item 17, spec 2026-06-12 SS3(e)): the edge
        # counterfactual -- what enlargement WOULD have seen absent the
        # gates. Ungated rows mirror their real writes; gate-zeroed edge
        # rows carry their counterfactual here and ONLY here. Warn-once
        # keying (gate_code, edge reaction index) clears per engine rebuild
        # = per RMG iteration: a persistent gated channel re-announces once
        # per iteration, deliberately (correct-but-loud; measured ~14
        # lines/iteration on the reference EPDM deck).
        self.edge_reaction_gate_code = np.zeros(self.num_edge_reactions,
                                                dtype=np.int8)
        self.edge_reaction_rates_ungated = np.zeros(self.num_edge_reactions,
                                                    float)
        self.edge_species_rates_ungated = np.zeros(self.num_edge_species,
                                                   float)
        self._phase_gate_census_emitted = set()

        self.get_const_spc_indices(core_species)
        self.set_initial_conditions()

        if filter_reactions:
            ReactionSystem.set_initial_reaction_thresholds(self)
        self.generate_rate_coefficients(core_reactions, edge_reactions)
        ReactionSystem.compute_network_variables(self, pdep_networks)

        ReactionSystem.set_initial_derivative(self)
        ReactionSystem.initialize_solver(self)

        self.diagnose_polymer_mapping(core_species)




    def diagnose_polymer_mapping(self, core_species):
        w = 90
        print(f"\n{'=' * w}")
        print(f"{'POLYMER SOLVER DIAGNOSTIC':^{w}}")
        print(f"{'=' * w}")

        # --- Per-pool sections ---
        for pool_i, pool in enumerate(self.polymer_pools):
            print(f"\n--- Pool {pool_i}: '{pool.label}'  (xs={pool.xs}, "
                  f"k_scission={pool.k_scission:.2e}, k_unzip={pool.k_unzip:.2e}) ---")
            print(f"  {'Role':<16} {'Species Label':<30} {'Index':<7} {'y0':>13}")
            print(f"  {'-'*16} {'-'*30} {'-'*7} {'-'*13}")

            # Moments
            for k, mu_label in enumerate(["mu0", "mu1", "mu2"]):
                idx = pool.mu_indices[k]
                lbl = core_species[idx].label if 0 <= idx < len(core_species) else "???"
                val = self.y0[idx] if 0 <= idx < len(self.y0) else float('nan')
                print(f"  {'Moment ('+mu_label+')':<16} {lbl:<30} {idx:<7} {val:>13.4e}")

            # Proxy species
            for i in range(len(core_species)):
                if self.species_to_pool_indices[i] == pool_i and self.is_pool_proxy[i] and i not in pool.mu_indices:
                    print(f"  {'Proxy':<16} {core_species[i].label:<30} {i:<7} {self.y0[i]:>13.4e}")

            # Explicit oligomers
            if pool.explicit_dp_to_species_index:
                for dp in sorted(pool.explicit_dp_to_species_index):
                    idx = pool.explicit_dp_to_species_index[dp]
                    lbl = core_species[idx].label if 0 <= idx < len(core_species) else "???"
                    val = self.y0[idx] if 0 <= idx < len(self.y0) else float('nan')
                    print(f"  {'Explicit DP='+str(dp):<16} {lbl:<30} {idx:<7} {val:>13.4e}")

            # Monomer
            if pool.monomer_poly_index is not None:
                idx = pool.monomer_poly_index
                lbl = core_species[idx].label if 0 <= idx < len(core_species) else "???"
                val = self.y0[idx] if 0 <= idx < len(self.y0) else float('nan')
                print(f"  {'Monomer(poly)':<16} {lbl:<30} {idx:<7} {val:>13.4e}")

            # Custom kinetics
            if pool.tail_kinetics is not None:
                print(f"  Custom tail_kinetics: YES")

        # --- Unassigned non-gas species (catch anything missed) ---
        unassigned = []
        for i in range(len(core_species)):
            if self.species_to_pool_indices[i] == -1 and not self.gas_species_mask[i]:
                unassigned.append(i)
        if unassigned:
            print(f"\n--- Unassigned non-gas species (not in any pool) ---")
            print(f"  {'Species Label':<30} {'Index':<7} {'y0':>13}")
            print(f"  {'-'*30} {'-'*7} {'-'*13}")
            for i in unassigned:
                print(f"  {core_species[i].label:<30} {i:<7} {self.y0[i]:>13.4e}")

        # --- Mass transfer ---
        if self.mass_transfer:
            print(f"\n--- Mass Transfer ---")
            for mt in self.mass_transfer:
                gas_lbl = core_species[mt.gas_index].label if 0 <= mt.gas_index < len(core_species) else "???"
                poly_lbl = core_species[mt.poly_index].label if 0 <= mt.poly_index < len(core_species) else "???"
                print(f"  {gas_lbl} (gas, idx={mt.gas_index}) <-> {poly_lbl} (poly, idx={mt.poly_index})  "
                      f"kLa={mt.kLa:.2e}")

        print(f"\n{'=' * w}\n")







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

    def _phase_gate_flux_census(self, core_species, edge_species,
                                edge_reactions, char_rate, tol_move_to_core):
        """RIDER R2 dynamic half (item 17, spec 2026-06-12 SS3(e)): emit one
        census line per (gate, edge reaction) per engine rebuild when the
        UNGATED ratio would have cleared the enlargement bar -- the
        species-level quantity base.pyx:1046 actually tests. Lives in
        simulate() via the base hook because tol_move_to_core is a simulate
        local and char_rate exists only per accepted snapshot; reads the
        most recent residual evaluation's arrays with EXACTLY the staleness
        of the enlargement read it audits (amendment A2 -- a feature).
        String assembly happens here (python level, once per emission),
        never in the residual."""
        if char_rate == 0.0:
            return  # the base.pyx singularity path owns the no-flux case
        gate_codes = getattr(self, "edge_reaction_gate_code", None)
        if gate_codes is None:
            return
        n_core = self.num_core_species
        n_core_rxn = self.num_core_reactions
        ip = self.product_indices
        for j in range(gate_codes.shape[0]):
            code = int(gate_codes[j])
            if code == 0:
                continue
            key = (code, j)
            if key in self._phase_gate_census_emitted:
                continue
            best_ratio = 0.0
            best_label = "<none>"
            condensed_labels = []
            for slot in range(3):
                p = ip[n_core_rxn + j, slot]
                if p == -1:
                    continue
                if not self.prospective_gas_mask[p]:
                    condensed_labels.append(
                        edge_species[p - n_core].label if p >= n_core
                        else core_species[p].label)
                if p >= n_core:
                    ratio = abs(self.edge_species_rates_ungated[p - n_core]) \
                        / char_rate
                    if ratio > best_ratio:
                        best_ratio = ratio
                        best_label = edge_species[p - n_core].label
            if best_ratio <= tol_move_to_core:
                continue
            self._phase_gate_census_emitted.add(key)
            rxn = edge_reactions[j]
            idx = n_core_rxn + j
            decisive = (", ".join(condensed_labels) if code == 1
                        else "no prospectively-condensed product")
            fam = (getattr(rxn, "family", None)
                   or getattr(rxn, "label", "") or type(rxn).__name__)
            logging.warning(
                "PHASE-GATE FLUX CENSUS: reaction %s gate=%s "
                "ungated_ratio=%.6e via edge species %s "
                "(tol_move_to_core=%.3e); decisive=%s; archetype "
                "pre-demotion=%d post-demotion=%d family=%s",
                rxn, "A" if code == 1 else "B", best_ratio, best_label,
                tol_move_to_core, decisive,
                int(self.reaction_pre_demotion_archetype[idx]),
                int(self.reaction_flux_archetype[idx]), fam)

    def _static_phase_gate_census(self, core_species, core_reactions):
        """RIDER R2 static half (item 17, spec 2026-06-12 SS3(e), amendment
        A1): enumerate CORE reactions whose phase-gate verdict is zero and
        announce each once per rebuild. Covers reactions that arrive
        core-gated without crossing the edge on their own flux: the third
        route (independent species promotion -- each participant promoted on
        other channels' flux) and legacy/restart cores. Gate verdicts for
        core rows are STATIC -- masks + event type only, no rates -- so this
        runs at init with zero residual cost; the kinetics loop's bare
        continue for core-gated rows survives untouched. The gate
        classification below MUST mirror residual section 6 (the
        prospective-mask gates) -- keep in sync."""
        n_core = self.num_core_species
        ir = self.reactant_indices
        ip = self.product_indices
        mask = self.gas_species_mask
        pmask = self.prospective_gas_mask
        for i in range(len(core_reactions)):
            r0, r1, r2 = ir[i, 0], ir[i, 1], ir[i, 2]
            # Defensive parity with the residual's reactant skips (core
            # reactions have all-core participants by construction).
            if r0 == -1 or r0 >= n_core:
                continue
            if (r1 != -1 and r1 >= n_core) or (r2 != -1 and r2 >= n_core):
                continue
            is_poly_event = ((not mask[r0])
                             or (r1 != -1 and not mask[r1])
                             or (r2 != -1 and not mask[r2]))
            has_condensed_prod = False
            condensed_labels = []
            for slot in range(3):
                p = ip[i, slot]
                if p == -1:
                    continue
                if not pmask[p]:
                    has_condensed_prod = True
                    if p < n_core:
                        condensed_labels.append(core_species[p].label)
            if (not is_poly_event) and has_condensed_prod:
                code = 1
            elif is_poly_event and not has_condensed_prod:
                code = 2
            else:
                continue
            rxn = core_reactions[i]
            decisive = (", ".join(condensed_labels) if code == 1
                        else "no prospectively-condensed product")
            fam = (getattr(rxn, "family", None)
                   or getattr(rxn, "label", "") or type(rxn).__name__)
            logging.warning(
                "PHASE-GATE FLUX CENSUS: reaction %s gate=%s static (core, "
                "init-time); decisive=%s; archetype pre-demotion=%d "
                "post-demotion=%d family=%s",
                rxn, "A" if code == 1 else "B", decisive,
                int(self.reaction_pre_demotion_archetype[i]),
                int(self.reaction_flux_archetype[i]), fam)

    def spawn_gate_flux_snapshot(self, motif_counts_by_pool=None):
        """Engine half of the mass-flux spawn gate.

        Spec: docs/superpowers/specs/2026-06-10-mass-flux-spawn-gate-design.md
        §3/§4.1 (AMENDED). The engine cannot attribute representatives (it
        has no ledger), so the split is: the engine reads arrays, the
        python-side gate does ledger-dependent attribution. Returns a
        3-tuple ``(gross, pool_stats, proxy_event_mass_total)``:

        1. ``gross``: dict core-species label ->
           ``max(0, core_species_production_rates[i])`` for ALL core species
           (ordinary species have real entries since change (a), spec §4.6).
           Labels are ``Species.label`` verbatim (the same labels the ledger
           records at absorption). Duplicate core labels RAISE ValueError —
           uniqueness is NOT enforced on the standard non-RMS path
           (model.py's dedup is gated on edge.phase_system), and silent
           overwrite would misattribute gate flux; the main.py caller turns
           the raise into warn + None snapshot -> the gate defers.
        2. ``pool_stats``: dict pool label -> ``(E_n, monomer_mw_g_mol)``
           with ``E_n = y[mu1]/y[mu0] if y[mu0] > trust else 0.0`` where
           ``trust = max(SMALL_EPS, ATTRIBUTION_TRUST_K * atol_mu0)`` is the
           ATTRIBUTION TRUST FLOOR (item #14a): mu0 at or below the
           integrator's own noise floor cannot substantiate a mass
           attribution, so the band (SMALL_EPS, trust] zeroes E_n (the band
           defers) and logs a ``SPAWN-GATE ATTRIBUTION CENSUS`` warning once
           per snapshot per pool. atol_mu0 is the pool's mu0-slot entry when
           the tolerance is a vector, the scalar atol otherwise. mu0
           exhausted (<= SMALL_EPS) still defers SILENTLY (unchanged). The
           solver's own E[n] consumers and the residual keep SMALL_EPS
           realizability semantics — errs toward deferral either way.
           ``motif_counts_by_pool`` (optional dict pool label -> int) lets
           the main.py stash report in the census line how many ledger
           motifs currently attribute to the zeroed pool (the engine has no
           ledger; the caller does).
        3. ``proxy_event_mass_total``: float — sum of
           ``gross_j * E_n[pool(j)] * mw[pool(j)]`` over CANONICAL pool
           proxies (``species_to_pool_indices[j] != -1 and is_pool_proxy[j]``
           — the engine CAN attribute those).

        GROSS production, never net dn_dt: canonical proxies have
        dn_dt ~= 0 BY DESIGN (the archetype apportionment reroutes their
        flux to pool moments) and ordinary species net to ~0 at steady
        state. E[n] is read LIVE from the state vector (never
        recorded-and-stale). All terms are polymer-phase volumetric, so
        V_poly cancels in any fraction of these numbers. Pure read of
        already-maintained state — no bookkeeping here beyond change (a)'s
        residual writes.
        """
        gross = {}
        pool_stats = {}
        proxy_event_mass_total = 0.0
        stp = getattr(self, "species_to_pool_indices", None)
        prod = getattr(self, "core_species_production_rates", None)
        y = getattr(self, "y", None)
        if stp is None or prod is None or y is None:
            return gross, pool_stats, proxy_event_mass_total
        n_core = self.num_core_species
        # Live per-pool stats, gated by the ATTRIBUTION TRUST FLOOR
        # (item #14a, spec 2026-06-11 §2/§3): per-pool mu0-slot atol when
        # the tolerance is a vector (the floor is PER-POOL — a deck mixing
        # large and trace pools gets per-pool trust, a feature), scalar atol
        # otherwise. atol_array unavailable (snapshot before
        # initialize_model) -> atol_mu0 = 0 -> the floor degenerates to
        # SMALL_EPS (pre-floor behavior, honest).
        n_pools = len(self.polymer_pools)
        e_n_by_pool = [0.0] * n_pools
        atol_arr = getattr(self, "atol_array", None)
        for p in range(n_pools):
            i0 = self.pool_mu0_indices[p]
            i1 = self.pool_mu1_indices[p]
            if i0 < 0 or i1 < 0:
                continue
            mu0 = y[i0]
            mu1 = y[i1]
            if atol_arr is None:
                atol_mu0 = 0.0
            elif np.ndim(atol_arr) == 0:
                atol_mu0 = float(atol_arr)
            elif i0 < len(atol_arr):
                atol_mu0 = float(atol_arr[i0])
            else:
                atol_mu0 = float(np.max(atol_arr))
            trust = max(SMALL_EPS, ATTRIBUTION_TRUST_K * atol_mu0)
            if mu0 > trust:
                e_n_by_pool[p] = mu1 / mu0
            elif mu0 > SMALL_EPS:
                # The distrust band (SMALL_EPS, trust]: integrator-noise-
                # scale mu0. Zero-out is the honest measurement ("this pool
                # cannot substantiate any mass attribution this iteration"),
                # drags the window statistic down, and preserves the
                # subsystem asymmetry: under-attribution defers, never
                # falsely spawns. The census line is the greppable reason a
                # deck's spawns mysteriously defer AND the standing sensor
                # for the floor ever biting legitimate chemistry (the
                # rejected-absolute-constant failure mode, made observable).
                n_motifs = 0
                if motif_counts_by_pool:
                    n_motifs = int(motif_counts_by_pool.get(
                        self.polymer_pools[p].label, 0))
                logging.warning(
                    "SPAWN-GATE ATTRIBUTION CENSUS: pool %s mu0=%.6e is in "
                    "the distrust band (SMALL_EPS=%.0e < mu0 <= trust "
                    "floor=%.6e = max(SMALL_EPS, %.0f * atol_mu0=%.6e)); "
                    "E[n] attribution zeroed for this snapshot; ledger "
                    "motifs attributed to this pool: %d.",
                    self.polymer_pools[p].label, mu0, SMALL_EPS, trust,
                    ATTRIBUTION_TRUST_K, atol_mu0, n_motifs)
            # mu0 <= SMALL_EPS: exhausted/empty -> 0.0 silently (unchanged).
        for p in range(n_pools):
            pool = self.polymer_pools[p]
            mw = float(getattr(pool, "monomer_mw_g_mol", 0.0) or 0.0)
            pool_stats[pool.label] = (e_n_by_pool[p], mw)
        # index -> label for CORE species (species_index covers core+edge).
        labels = {}
        for spc, idx in self.species_index.items():
            if idx < n_core:
                labels[idx] = spc.label
        for i in range(n_core):
            label = labels.get(i)
            if label is None:
                continue
            g = max(0.0, float(prod[i]))
            if label in gross:
                # Label uniqueness is NOT an enforced invariant on the
                # standard (non-RMS) path: model.py's dedup loop is gated on
                # edge.phase_system ("!!! Not maintained when
                # require_rms=False?"). A duplicate key would silently
                # misattribute gate flux. Raise instead: the main.py stash
                # wraps this call, logs a warning and leaves the snapshot
                # None -> the gate defers honestly (spec §4.5).
                raise ValueError(
                    f"duplicate core species label {label!r} in "
                    f"spawn_gate_flux_snapshot(); spawn-gate attribution "
                    f"would be unreliable")
            gross[label] = g
            p = stp[i]
            if p >= 0 and self.is_pool_proxy[i]:
                e_n, mw = pool_stats[self.polymer_pools[p].label]
                proxy_event_mass_total += g * e_n * mw
        return gross, pool_stats, proxy_event_mass_total

    def _chain_bundle(self, int pool_idx, y, double V_poly, bint end_group):
        """
        Expected (chains, units, units^2) carried by ONE picked chain of pool
        ``pool_idx``: (1, E[k], E[k^2]). end_group=True -> uniform chain pick
        (rate ~ mu0): (1, mu1/mu0, mu2/mu0). Otherwise length-biased pick
        (rate ~ mu1): (1, mu2/mu1, mu3/mu1) with the guarded mu3 closure.
        Returns (b0, b1, b2, mu2_ok); b0 == 0.0 means the pool is too empty to
        move a chain (denominator below SMALL_EPS) -- caller skips the term.
        mu2_ok False means apply b0/b1 but skip the mu2 component (mu3 = inf).
        """
        idx0, idx1, idx2 = self.polymer_pools[pool_idx].mu_indices
        mu0 = max(0.0, y[idx0]) / V_poly
        mu1 = max(0.0, y[idx1]) / V_poly
        mu2 = max(0.0, y[idx2]) / V_poly
        if end_group:
            if mu0 <= SMALL_EPS:
                return 0.0, 0.0, 0.0, False
            return 1.0, mu1 / mu0, mu2 / mu0, True
        if mu1 <= SMALL_EPS:
            return 0.0, 0.0, 0.0, False
        mu3 = _safe_mu3_from_mu012(mu0, mu1, mu2)
        if np.isfinite(mu3):
            return 1.0, mu2 / mu1, mu3 / mu1, True
        return 1.0, mu2 / mu1, 0.0, False

    @cython.boundscheck(False)
    def residual(self, double t, np.ndarray[np.float64_t, ndim=1] y,
                 np.ndarray[np.float64_t, ndim=1] dydt,
                 np.ndarray[np.float64_t, ndim=1] senpar = np.zeros(1, float)):
        """
        Compute the residual (dn/dt - dydt) at time t for state y.
        """
        cdef int i
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

        # Item 17 (spec 2026-06-12 SS3(b)): combined-length raise mirroring
        # the gas_species_mask raise above -- the product gates below index
        # prospective_gas_mask with raw product slots (core AND edge).
        if (self.prospective_gas_mask is None
                or self.prospective_gas_mask.shape[0] != n_core + n_edge):
            raise ValueError(
                f"State/Prospective-mask mismatch: n_core+n_edge="
                f"{n_core + n_edge}, prospective="
                f"{None if self.prospective_gas_mask is None else self.prospective_gas_mask.shape[0]}")

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
        pmask = self.prospective_gas_mask  # gate-input only (rider R3)
        C_gas[mask] = np.maximum(0.0, y[mask]) / V_gas
        C_poly[~mask] = np.maximum(0.0, y[~mask]) / V_poly

        # Sync diagnostic concentrations
        self.core_species_concentrations[mask] = C_gas[mask]
        self.core_species_concentrations[~mask] = C_poly[~mask]

        # 4. Clear Accumulators
        dn_dt = self._scratch_dn_dt
        dn_dt[:] = 0.0

        proxy_activity = self._scratch_proxy_activity
        proxy_activity[:] = 0.0

        self.core_reaction_rates[:] = 0.0
        self.edge_reaction_rates[:] = 0.0
        self.edge_species_rates[:] = 0.0
        self.edge_reaction_gate_code[:] = 0
        self.edge_reaction_rates_ungated[:] = 0.0
        self.edge_species_rates_ungated[:] = 0.0
        self.core_species_consumption_rates[:] = 0.0
        self.core_species_production_rates[:] = 0.0
        self.network_leak_rates[:] = 0.0

        # 5. Recalculate P-dependent rates in the GAS phase
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

        # 6. Standard Kinetics (Optimized Loop)
        for r_idx in range(ir.shape[0]):
            r0, r1, r2 = ir[r_idx, 0], ir[r_idx, 1], ir[r_idx, 2]

            if r0 == -1 or r0 >= n_core:
                continue

            r0_is_gas = mask[r0]
            r1_is_gas = True
            r2_is_gas = True

            if r1 != -1:
                if r1 >= n_core:
                    continue
                r1_is_gas = mask[r1]

            if r2 != -1:
                if r2 >= n_core:
                    continue
                r2_is_gas = mask[r2]

            is_poly_event = (not r0_is_gas) or (r1 != -1 and not r1_is_gas) or (r2 != -1 and not r2_is_gas)

            if is_poly_event:
                phase = 'poly'
                V_rxn = V_poly
            else:
                phase = 'gas'
                V_rxn = V_gas

            # Product Checks
            p0, p1, p2 = ip[r_idx, 0], ip[r_idx, 1], ip[r_idx, 2]

            prods_phase_ok = True
            has_edge_prod = False
            has_any_prod = False
            has_condensed_prod = False
            gate_code = 0  # 0 = open, 1 = Gate A, 2 = Gate B

            # PROSPECTIVE-MASK GATES (item 17, spec 2026-06-12 SS3(c)). The
            # phase verdict for EVERY product, core or edge alike, is
            # prospective_gas_mask[p]; for p < n_core it is bit-identical to
            # mask[p] by rider R1 (core-prefix parity raise at init), so
            # CORE-product behavior is unchanged by construction. The old
            # has_edge_prod PHASE bypass (maskless edge product -> exempt)
            # is dead: it is what let enlargement promote on flux the
            # post-promotion model zeroes (the umbrella invariant's mask
            # projection). has_edge_prod survives ONLY for the reverse-rate
            # concentration-availability hole below (:rr block) -- an edge
            # product has no state in y, so rr is UNCOMPUTABLE there; that
            # is simple.pyx parity (Z6), not a phase verdict.
            def _check_prod(p):
                nonlocal has_edge_prod, has_any_prod, has_condensed_prod
                if p == -1:
                    return
                has_any_prod = True
                if p >= n_core:
                    has_edge_prod = True
                if not pmask[p]:
                    has_condensed_prod = True

            _check_prod(p0)
            _check_prod(p1)
            _check_prod(p2)

            if (not is_poly_event) and has_condensed_prod:
                # Gate A: a gas event with ANY prospectively-condensed
                # product (core or edge) is zeroed.
                prods_phase_ok = False
                gate_code = 1
            elif is_poly_event and (not has_condensed_prod):
                # Gate B: a poly event with NO prospectively-condensed
                # product (core or edge) is zeroed.
                prods_phase_ok = False
                gate_code = 2

            gated = False
            if not prods_phase_ok:
                if r_idx >= n_rxn:
                    # RIDER R2 dynamic half: record the gate verdict and fall
                    # through the EXISTING rate computation (pool mapping,
                    # _C, site scaling, throttle; the rr hole applies as
                    # usual) with every REAL write suppressed -- the zeros in
                    # edge_reaction_rates/edge_species_rates ARE the
                    # consistency point; the counterfactual lands only in
                    # the *_ungated arrays.
                    self.edge_reaction_gate_code[r_idx - n_rxn] = gate_code
                    gated = True
                else:
                    # Core-gated rows keep the bare continue: zero residual
                    # cost. Their loudness is the STATIC census at
                    # initialize_model (spec SS3(e) static half).
                    continue

            # 1. Map Reactants to Polymer Pools  (MOVED UP before _C)
            p0_pool_idx = self.species_to_pool_indices[r0]
            p1_pool_idx = -1 if r1 == -1 else self.species_to_pool_indices[r1]
            p2_pool_idx = -1 if r2 == -1 else self.species_to_pool_indices[r2]  # optional but consistent

            # Rate Calculation
            kf = self.kf[r_idx]
            kb = self.kb[r_idx]

            def _C(idx):
                if idx == -1:
                    return 1.0
                if self.species_to_pool_indices[idx] != -1:
                    return 1.0
                return C_gas[idx] if mask[idx] else C_poly[idx]

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

            # 2. Polymer scaling
            if p0_pool_idx != -1 or p1_pool_idx != -1 or p2_pool_idx != -1:
                if p0_pool_idx != -1 or p1_pool_idx != -1 or p2_pool_idx != -1:
                    # Determine the target pool index
                    target_pool_idx = p0_pool_idx
                    if target_pool_idx == -1: target_pool_idx = p1_pool_idx
                    if target_pool_idx == -1: target_pool_idx = p2_pool_idx

                    # Default to Mu1 (Site Density)
                    moment_idx = self.polymer_pools[target_pool_idx].mu_indices[1]

                    # Check for specific "End Group" physics (Tier 1 / Initiation)
                    if self.is_end_group_reaction[r_idx]:
                        moment_idx = self.polymer_pools[target_pool_idx].mu_indices[0]  # Scale by Mu0 (Chain Density)

                    # Calculate the site/chain concentration
                    # Note: We use the *mapped* index from the pool, ensuring we grab the real state variable
                    site = max(0.0, y[moment_idx]) / V_poly

                    # Exhaustion throttle (spec 2026-06-10 s5 AMENDMENT): a
                    # mu0-scaled DISCRETE_CHIP drains mu1 but never mu0, so
                    # its rate would not self-limit -- mu1 would run linearly
                    # negative past exhaustion while chips keep being created.
                    # Counting inequality (not a heuristic): every chain with
                    # n >= a holds at least a units, so donatable ends are
                    # bounded by mu1/a always. Throttling the SITE (not
                    # write-gating the mu1 leg, which would duplicate mass)
                    # gives dmu1/dt >= -kf*mu1 globally; it multiplies BOTH
                    # directions by design. a = 0 chips (drain nothing) and
                    # mu1-scaled chips (already self-limiting) are exempt.
                    # (mirrored in get_reaction_rates' hijack block -- keep in sync)
                    if (self.reaction_flux_archetype[r_idx] == FLUX_DISCRETE_CHIP
                            and self.is_end_group_reaction[r_idx]
                            and self.reaction_chip_units[r_idx] > 0):
                        mu_idx = self.polymer_pools[target_pool_idx].mu_indices
                        site = min(
                            max(0.0, y[mu_idx[0]]),
                            max(0.0, y[mu_idx[1]]) / float(self.reaction_chip_units[r_idx]),
                        ) / V_poly

                    rf *= site
                    rr *= site

            # Net rate (volumetric, in the phase volume chosen earlier)
            rate = rf - rr
            r_mol_s = rate * V_rxn
            abs_flux = abs(r_mol_s)

            if r_idx < n_rxn:
                self.core_reaction_rates[r_idx] = rate
            elif gated:
                # counterfactual only -- the real entry stays 0 (the
                # consistency point enlargement reads).
                self.edge_reaction_rates_ungated[r_idx - n_rxn] = rate
            else:
                self.edge_reaction_rates[r_idx - n_rxn] = rate
                self.edge_reaction_rates_ungated[r_idx - n_rxn] = rate

            core_rxn = r_idx < n_rxn

            # 3. Apply Fluxes to Reactants (core reactions only -- edge
            #    reactions are diagnostic-only, matching simple.pyx).
            #    proxy_activity stays UNGATED on purpose: it feeds proxy
            #    similarity/spawn diagnostics, which want edge flux too --
            #    except gate-zeroed edge rows, suppressed just below (T9).
            #    Pool MOMENT flux is handled once per reaction in section 5.
            if self.is_pool_proxy[r0]:
                if not gated:
                    # Ghost-flux suppression (spec SS3(e)/T9): a gate-zeroed
                    # edge row's |flux| is judged unphysical under core
                    # semantics and must not feed spawn/similarity
                    # diagnostics.
                    proxy_activity[r0] += abs_flux
                if core_rxn:
                    self.core_species_consumption_rates[r0] += rf
                    self.core_species_production_rates[r0] += rr
            elif core_rxn:
                dn_dt[r0] -= r_mol_s
                # Change (a) (spec 2026-06-10 mass-flux-spawn-gate s3.1/s4.6):
                # gross bookkeeping for ORDINARY core species too (simple.pyx
                # parity), so absorbed-representative production has a real
                # record for spawn_gate_flux_snapshot(). Cost-gated on the
                # EPDM deck (<= ~5% wall-clock; fallback: ledger-tracked-only
                # writes via a flag array).
                self.core_species_consumption_rates[r0] += rf
                self.core_species_production_rates[r0] += rr

            if r1 != -1:
                if self.is_pool_proxy[r1]:
                    if not gated:
                        proxy_activity[r1] += abs_flux
                    if core_rxn:
                        self.core_species_consumption_rates[r1] += rf
                        self.core_species_production_rates[r1] += rr
                elif core_rxn:
                    dn_dt[r1] -= r_mol_s
                    self.core_species_consumption_rates[r1] += rf
                    self.core_species_production_rates[r1] += rr

            if r2 != -1 and core_rxn:
                # r2 is always ordinary today (mirror the proxy logic above
                # if polymer is ever allowed in r2); change (a) gross writes
                # apply as for any ordinary reactant.
                dn_dt[r2] -= r_mol_s
                self.core_species_consumption_rates[r2] += rf
                self.core_species_production_rates[r2] += rr

            # 4. Apply Fluxes to Products (same core_rxn gate as section 3)
            for p_slot in range(3):
                p_idx_tmp = ip[r_idx, p_slot]
                if p_idx_tmp == -1:
                    continue

                if p_idx_tmp < n_core:
                    if self.is_pool_proxy[p_idx_tmp]:
                        if not gated:
                            proxy_activity[p_idx_tmp] += abs_flux
                        if core_rxn:
                            self.core_species_production_rates[p_idx_tmp] += rf
                            self.core_species_consumption_rates[p_idx_tmp] += rr
                    elif core_rxn:
                        # feature polymer species (e.g., PS_rad, scission oligomer, etc.) stays explicit
                        dn_dt[p_idx_tmp] += r_mol_s
                        # Change (a): gross bookkeeping for ordinary products
                        # (simple.pyx parity; see the reactant-side comment).
                        self.core_species_production_rates[p_idx_tmp] += rf
                        self.core_species_consumption_rates[p_idx_tmp] += rr
                else:
                    self.edge_species_rates_ungated[p_idx_tmp - n_core] += rate
                    if not gated:
                        self.edge_species_rates[p_idx_tmp - n_core] += rate

            # 5. Pool moment flux -- archetype dispatch (core reactions only).
            #    Spec: docs/superpowers/specs/2026-06-09-proxy-moment-flux-
            #    apportionment-design.md. SAME_POOL and NONE apply nothing:
            #    fold-back flux is net-zero by construction, so skipping avoids
            #    roundoff and closure calls.
            if core_rxn:
                arch = self.reaction_flux_archetype[r_idx]
                if arch == FLUX_MIGRATION:
                    src = self.reaction_src_pool[r_idx]
                    dst = self.reaction_dst_pool[r_idx]
                    # -1 cannot reach here (init demotes unresolved pools to
                    # UNRESOLVED); the checks are defensive only.
                    if src != -1 and dst != -1 and src != dst:
                        # Per-direction bundles: forward moves A-statistics
                        # chains A->B, reverse moves B-statistics chains B->A.
                        # Each direction is guarded by its OWN source pool.
                        for ev_rate, from_pool, to_pool in (
                                (rf, src, dst), (rr, dst, src)):
                            if ev_rate <= 0.0:
                                continue
                            ev_mol = ev_rate * V_rxn
                            b0, b1, b2, mu2_ok = self._chain_bundle(
                                from_pool, y, V_poly,
                                self.is_end_group_reaction[r_idx])
                            if b0 == 0.0:
                                continue
                            f_idx = self.polymer_pools[from_pool].mu_indices
                            t_idx = self.polymer_pools[to_pool].mu_indices
                            dn_dt[f_idx[0]] -= ev_mol * b0
                            dn_dt[f_idx[1]] -= ev_mol * b1
                            dn_dt[t_idx[0]] += ev_mol * b0
                            dn_dt[t_idx[1]] += ev_mol * b1
                            if mu2_ok:
                                dn_dt[f_idx[2]] -= ev_mol * b2
                                dn_dt[t_idx[2]] += ev_mol * b2
                elif arch == FLUX_SCISSION_FRAGMENT:
                    src = self.reaction_src_pool[r_idx]
                    dst = self.reaction_dst_pool[r_idx]
                    # -1 cannot reach here (init demotes unresolved pools).
                    # src == dst (fold-back + non-pool daughter) is skipped:
                    # outcome-identical to the legacy transfer (see the init
                    # demotion comment).
                    if src != -1 and dst != -1 and src != dst:
                        s_idx = self.polymer_pools[src].mu_indices
                        d_idx = self.polymer_pools[dst].mu_indices
                        # Length-biased parent-chain statistics; mirrors
                        # _chain_bundle's non-end-group branch (kept inline
                        # because the scission factors 1/2, 2/3, 1/3 differ
                        # per moment and per side).
                        mu0_p = max(0.0, y[s_idx[0]]) / V_poly
                        mu1_p = max(0.0, y[s_idx[1]]) / V_poly
                        mu2_p = max(0.0, y[s_idx[2]]) / V_poly
                        ok = mu1_p > SMALL_EPS
                        if ok and r_mol_s < 0.0:
                            # Net reverse = coupling bookkeeping; it depletes
                            # the DAUGHTER, so guard its moments too. (Stated
                            # approximation: parent statistics, sign-flipped.)
                            if (max(0.0, y[d_idx[0]]) / V_poly <= SMALL_EPS or
                                    max(0.0, y[d_idx[1]]) / V_poly <= SMALL_EPS):
                                ok = False
                        if ok:
                            # Complement-stays-in-parent accounting: parent
                            # mu0 net 0 (chain broke, complement remains);
                            # fragment (uniform cut of a length-biased chain:
                            # E[a] = E[n]/2, E[a^2] = E[n^2]/3) moves to the
                            # daughter. mu1 conserves exactly per reaction.
                            e_n = mu2_p / mu1_p
                            dn_dt[s_idx[1]] -= r_mol_s * e_n / 2.0
                            dn_dt[d_idx[0]] += r_mol_s
                            dn_dt[d_idx[1]] += r_mol_s * e_n / 2.0
                            mu3_p = _safe_mu3_from_mu012(mu0_p, mu1_p, mu2_p)
                            if np.isfinite(mu3_p):
                                e_n2 = mu3_p / mu1_p
                                dn_dt[s_idx[2]] -= r_mol_s * (2.0 / 3.0) * e_n2
                                dn_dt[d_idx[2]] += r_mol_s * e_n2 / 3.0
                elif arch == FLUX_DISCRETE_CHIP:
                    # Discrete chip (spec 2026-06-10 §5): an end-anchored cut
                    # ejects a stamped a-unit chip to the gas phase; the
                    # complement folds back into the SAME pool, so mu0 is
                    # unchanged and there is no dst pool. Closure-free in
                    # effect: only E[n] (an exact ratio of tracked moments)
                    # is consumed -- b2/mu2_ok are ignored here, so no mu2_ok
                    # branch is needed. (For mu1-scaled chips _chain_bundle's
                    # length-biased path still computes the mu3 closure
                    # internally; its output is unused by this branch.) The
                    # bundle and the rate scaling key on the same stored flag;
                    # in the throttled exhaustion regime the rate measure is
                    # donor-limited while the pick stays population-based --
                    # accepted approximation, see spec s5 amendment. The chip
                    # species itself gains/loses moles through the standard
                    # gas dn_dt path (section 4) -- no special handling here.
                    src = self.reaction_src_pool[r_idx]
                    # src == -1 cannot reach here (init demotes); defensive.
                    if src != -1:
                        chip_a = float(self.reaction_chip_units[r_idx])
                        b0, b1, _b2, _mu2_ok = self._chain_bundle(
                            src, y, V_poly, self.is_end_group_reaction[r_idx])
                        if b0 != 0.0:
                            s_idx = self.polymer_pools[src].mu_indices
                            chip_e_n = b1
                            if rf > 0.0:
                                # Forward (ejection): Delta n = -a,
                                # Delta n^2 = -(2aE[n] - a^2). Clamp the mu2
                                # decrement at >= 0: 2aE[n] - a^2 < 0 is
                                # impossible per-chain (n >= a always) but
                                # reachable in expectation when the pool mean
                                # length has decayed toward chip size; by
                                # then the moment description is marginal.
                                rf_mol = rf * V_rxn
                                dn_dt[s_idx[1]] -= rf_mol * chip_a
                                chip_dmu2 = 2.0 * chip_a * chip_e_n - chip_a * chip_a
                                if chip_dmu2 > 0.0:
                                    dn_dt[s_idx[2]] -= rf_mol * chip_dmu2
                            if rr > 0.0:
                                # Reverse (re-attachment) -- EXACT extension
                                # form, NOT the forward sign-flip:
                                # (n+a)^2 - n^2 = +(2aE[n] + a^2). E[n] on the
                                # same single pool (the re-formed chain
                                # extends a parent-statistics chain by a).
                                # Unconditionally positive: no clamp.
                                rr_mol = rr * V_rxn
                                dn_dt[s_idx[1]] += rr_mol * chip_a
                                dn_dt[s_idx[2]] += rr_mol * (
                                    2.0 * chip_a * chip_e_n + chip_a * chip_a)
                elif arch == FLUX_UNRESOLVED:
                    # Legacy mu1-only transfer (pre-apportionment behavior),
                    # replicated exactly: -r per reactant proxy, +r per
                    # product proxy. NOTE: mu0-scaled shapes here share the
                    # chip exhaustion structure but are deliberately NOT
                    # throttled (bit-exact legacy contract; see
                    # docs/multi_pool_design.md limitation 14 / spec A4).
                    if self.is_pool_proxy[r0]:
                        dn_dt[self.polymer_pools[p0_pool_idx].mu_indices[1]] -= r_mol_s
                    if r1 != -1 and self.is_pool_proxy[r1]:
                        dn_dt[self.polymer_pools[p1_pool_idx].mu_indices[1]] -= r_mol_s
                    for p_slot in range(3):
                        p_idx_tmp = ip[r_idx, p_slot]
                        if (p_idx_tmp != -1 and p_idx_tmp < n_core
                                and self.is_pool_proxy[p_idx_tmp]):
                            pool_idx = self.species_to_pool_indices[p_idx_tmp]
                            dn_dt[self.polymer_pools[pool_idx].mu_indices[1]] += r_mol_s

        # 7. Network Leaks
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

        # 8. Polymer Tail & Handshake
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

            if self.debug_check_realizability and pool.label not in self._realizability_warned:
                # Use raw (pre-clamp) moles so negatives are caught too. All terms
                # are extensive, so the inequalities are scale-invariant in V_poly.
                raw0, raw1, raw2 = y[idx_mu0], y[idx_mu1], y[idx_mu2]
                if (raw0 < -1e-9 or raw1 < -1e-9 or raw2 < -1e-9
                        or raw1 + 1e-9 < raw0
                        or raw0 * raw2 + 1e-9 < raw1 * raw1):
                    self._realizability_warned.add(pool.label)
                    logging.warning(
                        "Polymer pool '%s' moment state left the realizable cone "
                        "(mu0=%.6g, mu1=%.6g, mu2=%.6g; require mu1>=mu0>=0 and "
                        "mu0*mu2>=mu1^2). The mu3 closure is guarded, but a moment "
                        "source term is likely wrong.", pool.label, raw0, raw1, raw2)

            dmu0_dt = 0.0
            dmu1_dt = 0.0
            dmu2_dt = 0.0
            small_src = dict()

            if pool.tail_kinetics:
                dmu0_dt, dmu1_dt, dmu2_dt, small_src = pool.tail_kinetics(
                    self.T.value_si, self.P.value_si, mu0, mu1, mu2, mu3)
            else:
                if pool.k_scission > 0:
                    # Random backbone scission, discrete-bond (Ziff-McGrady)
                    # convention: a chain of length k has (k-1) breakable bonds,
                    # so the distribution's breakable-bond count is (mu1 - mu0):
                    #   dμ0/dt = k_s·(μ1 − μ0)     (one new chain per bond broken)
                    #   dμ1/dt = 0                 (monomer units / mass conserved)
                    #   dμ2/dt = (k_s/3)·(μ1 − μ3)
                    # See docs/multi_pool_design.md §5. The μ0 term MUST carry the
                    # −μ0 depletion part: with the bare +k_s·μ1 form μ0 grows past
                    # μ1, the state leaves the realizable cone (μ1 ≥ μ0 always for a
                    # k≥1 distribution), and the μ3 closure blows up to a DASSL
                    # singularity. The (μ1 − μ0) form self-limits (rate → 0 as the
                    # pool reaches all-monomer) and structurally keeps μ0 ≤ μ1.
                    dmu0_dt += pool.k_scission * (mu1 - mu0)
                    if np.isfinite(mu3):
                        dmu2_dt += pool.k_scission * (mu1 - mu3) / 3.0

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

        # 9. Mass Transfer
        for mt in self.mass_transfer:
            ig = mt.gas_index
            ipoly = mt.poly_index

            Cg = C_gas[ig]
            Cp = C_poly[ipoly]

            J = mt.kLa * (Cp - mt.K * Cg)
            dn = J * V_poly

            dn_dt[ig] += dn
            dn_dt[ipoly] -= dn

        # 10. Constants
        if self.const_spc_indices is not None:
            for i in self.const_spc_indices:
                dn_dt[i] = 0.0

        # 11. Diagnostics
        self.core_species_rates[mask] = dn_dt[mask] / V_gas
        self.core_species_rates[~mask] = dn_dt[~mask] / V_poly

        for i in range(n_core):
            if self.is_pool_proxy[i]:
                self.core_species_rates[i] = proxy_activity[i] / V_poly

        return (dn_dt - dydt), 1

    def get_reaction_rates(self, y_in):
        """
        Return the net rates of all core and edge reactions at state y.
        """
        # 1. Cast Input to Memoryview
        cdef np.ndarray[np.float64_t, ndim=1] y_arr = np.asarray(y_in, dtype=np.float64)
        cdef double[:] y = y_arr

        cdef int n_core = self.num_core_species
        cdef int n_rxn = self.num_core_reactions + self.num_edge_reactions

        # Output Array
        cdef np.ndarray[np.float64_t, ndim=1] rates_arr = np.zeros(n_rxn, float)
        cdef double[:] reaction_rates = rates_arr

        # Scratch Arrays
        cdef np.ndarray[np.float64_t, ndim=1] Cg_arr = np.zeros(n_core, float)
        cdef double[:] C_gas = Cg_arr

        cdef np.ndarray[np.float64_t, ndim=1] Cp_arr = np.zeros(n_core, float)
        cdef double[:] C_poly = Cp_arr

        # Declarations
        cdef int r_idx, r0, r1, r2, p0_pool_idx
        cdef int i  # <--- Moved declaration to top
        cdef int p0, p1, p2
        cdef double kf, kb, rf, rr, site, V_rxn
        cdef double V_gas, V_poly
        cdef double n_gas

        # 2. Volumes
        V_poly = self.V_poly
        if self.constant_gas_volume:
            V_gas = self.V_gas0
        else:
            n_gas = 0.0
            for i in range(n_core):
                if self.gas_species_mask[i]:
                    n_gas += y[i]

            V_gas = constants.R * self.T.value_si * n_gas / self.P.value_si if n_gas > 0 else 1.0

        # 3. Concentrations
        for i in range(n_core):
            if self.gas_species_mask[i]:
                C_gas[i] = max(0.0, y[i]) / V_gas
            else:
                C_poly[i] = max(0.0, y[i]) / V_poly

        ir = self.reactant_indices
        ip = self.product_indices

        # 4. Kinetics Loop
        for r_idx in range(n_rxn):
            r0 = ir[r_idx, 0]
            r1 = ir[r_idx, 1]
            r2 = ir[r_idx, 2]

            # Determine Phase for Volume Basis
            is_poly = False
            if r0 != -1 and not self.gas_species_mask[r0]:
                is_poly = True
            elif r1 != -1 and not self.gas_species_mask[r1]:
                is_poly = True
            elif r2 != -1 and not self.gas_species_mask[r2]:
                is_poly = True

            V_rxn = V_poly if is_poly else V_gas

            kf = self.kf[r_idx]
            kb = self.kb[r_idx]

            # Forward Rate
            rf = kf
            if self.species_to_pool_indices[r0] != -1:
                rf *= 1.0
            elif self.gas_species_mask[r0]:
                rf *= C_gas[r0]
            else:
                rf *= C_poly[r0]

            if r1 != -1:
                if self.species_to_pool_indices[r1] != -1:
                    rf *= 1.0
                elif self.gas_species_mask[r1]:
                    rf *= C_gas[r1]
                else:
                    rf *= C_poly[r1]

            if r2 != -1:
                if self.species_to_pool_indices[r2] != -1:
                    rf *= 1.0
                elif self.gas_species_mask[r2]:
                    rf *= C_gas[r2]
                else:
                    rf *= C_poly[r2]

            # Reverse Rate
            rr = 0.0
            if kb > 0:
                p0 = ip[r_idx, 0]
                p1 = ip[r_idx, 1]
                p2 = ip[r_idx, 2]
                rr = kb

                if p0 != -1:
                    if self.species_to_pool_indices[p0] != -1:
                        rr *= 1.0
                    elif self.gas_species_mask[p0]:
                        rr *= C_gas[p0]
                    else:
                        rr *= C_poly[p0]

                if p1 != -1:
                    if self.species_to_pool_indices[p1] != -1:
                        rr *= 1.0
                    elif self.gas_species_mask[p1]:
                        rr *= C_gas[p1]
                    else:
                        rr *= C_poly[p1]

                if p2 != -1:
                    if self.species_to_pool_indices[p2] != -1:
                        rr *= 1.0
                    elif self.gas_species_mask[p2]:
                        rr *= C_gas[p2]
                    else:
                        rr *= C_poly[p2]

            # [THE HIJACK] Scale by Moment if Reactant is Proxy
            p0_pool_idx = self.species_to_pool_indices[r0] if r0 != -1 else -1
            if p0_pool_idx != -1:
                moment_idx = self.polymer_pools[p0_pool_idx].mu_indices[1]

                if self.pool_mu1_indices[p0_pool_idx] != -1:
                    moment_idx = self.pool_mu1_indices[p0_pool_idx]

                if self.is_end_group_reaction[r_idx]:
                    moment_idx = self.pool_mu0_indices[p0_pool_idx]

                site = max(0.0, y[moment_idx]) / V_poly

                # Exhaustion throttle -- parity with the residual's exhaustion
                # throttle (spec 2026-06-10 s5 AMENDMENT): the diagnostic rate
                # must not report a flux the residual refuses to deliver for a
                # mu0-scaled DISCRETE_CHIP past unit exhaustion. Same condition
                # and formula as the residual's section-2 site scaling; here
                # moment_idx is already the resolved mu0 index (end-group
                # branch above), and the mu1 index is resolved with the same
                # pool_mu1_indices override as the default branch.
                if (self.reaction_flux_archetype[r_idx] == FLUX_DISCRETE_CHIP
                        and self.is_end_group_reaction[r_idx]
                        and self.reaction_chip_units[r_idx] > 0):
                    mu1_idx = self.polymer_pools[p0_pool_idx].mu_indices[1]
                    if self.pool_mu1_indices[p0_pool_idx] != -1:
                        mu1_idx = self.pool_mu1_indices[p0_pool_idx]
                    site = min(
                        max(0.0, y[moment_idx]),
                        max(0.0, y[mu1_idx]) / float(self.reaction_chip_units[r_idx]),
                    ) / V_poly

                rf *= site
                rr *= site

            reaction_rates[r_idx] = (rf - rr)

        return rates_arr

    def get_edge_reaction_rates(self, core_y, edge_rates):
        """
        Calculate ONLY edge reaction rates using the hijacked logic.
        """
        # Call full calculation
        full_rates = self.get_reaction_rates(core_y)

        # Copy to output buffer
        cdef int n_core_rxn = self.num_core_reactions
        cdef int n_edge_rxn = self.num_edge_reactions
        cdef double[:] edge_view = edge_rates  # Cast input to view
        cdef double[:] full_view = full_rates

        cdef int i
        for i in range(n_edge_rxn):
            edge_view[i] = full_view[n_core_rxn + i]
