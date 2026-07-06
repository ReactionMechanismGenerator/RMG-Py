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

import copy
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
# Exhaustion-tail census/tripwire floor multiplier (adjudicated round 81,
# the PP run-7 IDID=-7 wall; re-adjudicated round 86: DIAGNOSTIC ONLY --
# the r81 RHS read-projection is reverted, the RHS always computes from the
# raw state): per-state floor = max(SMALL_EPS, K * atol[state]).
# A pool is counted EXHAUSTED for the census only when |mu0|, |mu1| AND
# |mu2| ALL sit below their floors (never mu0 alone -- tiny mu0 with
# nontrivial mu1 is a few long chains, not an empty pool). K = 100 is the
# same two-decade tolerance-anchored margin as ATTRIBUTION_TRUST_K (a
# separate constant on purpose: the spawn-gate trust band and the solver
# exhaustion band are two jobs that may diverge). Any moment more negative
# than -floor is a HARD error (integrator corruption), never
# max(..., SMALL_EPS). Generic solver infrastructure OUTSIDE the 2.8 kernel
# contract -- the sidecar recipe strings are deliberately untouched.
EXHAUSTION_FLOOR_K = 100.0
LN_EXP_OVERFLOW_GUARD = 700.0
TAIL_CONC_MIN = 1e-9  # Minimum concentration (mol/m^3) to actuate handshake

# Gas constant for the radical_qssa_unzip Arrhenius evaluations, J/(mol K).
# Pinned to 8.314 by the M1/M2 channel contract (Ea is stored in J/mol by
# convention, not dimensionally enforced) -- deliberately NOT constants.R
# so the channel's k(T) is bit-reproducible against the documented law.
QSSA_R_GAS = 8.314

# +inf sentinel for the RHS runtime degenerate-rate guard: `x < QSSA_INF` is
# a single comparison that rejects BOTH inf and NaN (NaN fails every
# comparison), keeping the per-RHS-call hot path to a plain compare chain.
QSSA_INF = float("inf")


def _raise_degenerate_qssa_rates(pool_label, T, ki, kdp, kt, ktr=None,
                                 kia=None, ktrec=None, ktdisp=None):
    """Cold path for the radical_qssa_unzip runtime degenerate-rate guard.

    M1 validates config FINITENESS, but Arrhenius EVALUATION at the solver
    temperature can still degenerate: kt(T) can underflow to exactly 0.0
    (R_ss = sqrt(f*ki*B/kt) then divides by zero into inf) and a large A or n
    can overflow A*T**n to inf (NaN downstream). Fail loud, naming the pool,
    the temperature and every offending constant -- never silently zero or
    clamp the rate. Only called off the hot path, when the cheap comparison
    chain in the RHS has already found a degenerate value."""
    problems = []
    for name, val in (("initiation ki", ki), ("depropagation kdp", kdp),
                      ("termination kt", kt), ("transfer ktr", ktr),
                      ("initiation_allyl kia", kia),
                      ("termination_recombination ktrec", ktrec),
                      ("termination_disproportionation ktdisp", ktdisp)):
        if val is None:
            continue
        if not math.isfinite(val):
            problems.append(f"{name}(T)={val!r} is non-finite "
                            f"(Arrhenius overflow at evaluation)")
        elif name == "termination kt" and val <= 0.0:
            problems.append(
                f"{name}(T)={val!r} underflowed to <= 0 "
                f"(R_ss = sqrt(f*ki*B/kt) would divide by zero)")
    raise ValueError(
        f"Pool {pool_label}: radical_qssa_unzip rate evaluation is "
        f"degenerate at T={T:g} K: " + "; ".join(problems) + ". Config-time "
        f"validation (M1) pins finiteness of A/n/Ea but cannot catch "
        f"runtime overflow/underflow of k(T) = A*T**n*exp(-Ea/(R*T)); fix "
        f"the offending Arrhenius block for the reactor temperature range.")

# Pool moment-flux archetypes. Mirror of rmgpy.polymer.PolymerFluxArchetype
# (not imported to avoid a solver->polymer module cycle); equality is pinned
# by test_flux_archetype_constants_match_enum.
FLUX_NONE = 0
FLUX_SAME_POOL = 1
FLUX_MIGRATION = 2
FLUX_SCISSION_FRAGMENT = 3
FLUX_UNRESOLVED = 4
FLUX_DISCRETE_CHIP = 5
FLUX_VOLATILE_EJECTION = 6


# ======================================================================================
# DATA CONFIGURATION
# ======================================================================================

# Pinned basis string for the radical QSSA unzip channel (forward-compat pin):
# the M2 rate law will count initiation sites from the pool's backbone-bond
# count (mu1 - mu0). Any other basis is rejected until a rate law that
# consumes it actually ships, so a deck cannot silently request semantics
# that do not exist yet.
RADICAL_QSSA_UNZIP_BASIS = "backbone_bonds_mu1_minus_mu0"

_RADICAL_QSSA_ARRHENIUS_KEYS = ("A", "n", "Ea")
_RADICAL_QSSA_MANDATORY_BLOCKS = ("initiation", "depropagation", "termination")

# Weak-link allyl/U-state vocabulary (schema-2.2 milestone i: config +
# validation ONLY -- no RHS reads it yet, and the sidecar schema bump is a
# later milestone). ALL-OR-NOTHING as a group, and mutually exclusive with
# the legacy SUMMED 'termination' block: the U-state (unsaturated tail ends)
# is PRODUCED by the disproportionation branch specifically, so a summed kt
# cannot source U -- the split blocks are structurally required.
_RADICAL_QSSA_WEAKLINK_TRIPLET_BLOCKS = (
    "initiation_allyl", "termination_recombination",
    "termination_disproportionation")
_RADICAL_QSSA_WEAKLINK_KEYS = _RADICAL_QSSA_WEAKLINK_TRIPLET_BLOCKS + (
    "unsaturated_tail_ends_initial",)

_RADICAL_QSSA_ALLOWED_KEYS = frozenset(
    _RADICAL_QSSA_MANDATORY_BLOCKS
    + ("transfer", "efficiency", "monomer_yield", "basis")
    + _RADICAL_QSSA_WEAKLINK_KEYS)


def _validate_qssa_arrhenius_triplet(pool_label, block_name, triplet):
    """Validate one radical_qssa_unzip Arrhenius block {A, n, Ea}; return it
    normalized to plain floats.

    Rules (from day one, unlike the legacy k_unzip scalar): every value must
    be a real number and FINITE (NaN/inf rejected explicitly), A > 0,
    Ea >= 0, n any finite float. Raises ValueError naming the pool and the
    offending block/key.
    """
    if not isinstance(triplet, dict):
        raise ValueError(
            f"Pool {pool_label}: radical_qssa_unzip block '{block_name}' must "
            f"be a dict {{A, n, Ea}}, got {type(triplet).__name__}.")
    missing = [k for k in _RADICAL_QSSA_ARRHENIUS_KEYS if k not in triplet]
    extra = sorted(set(triplet) - set(_RADICAL_QSSA_ARRHENIUS_KEYS))
    if missing or extra:
        raise ValueError(
            f"Pool {pool_label}: radical_qssa_unzip block '{block_name}' must "
            f"have exactly the keys {{A, n, Ea}}; missing {missing or 'none'}, "
            f"unknown {extra or 'none'}.")
    out = {}
    for key in _RADICAL_QSSA_ARRHENIUS_KEYS:
        val = triplet[key]
        if isinstance(val, bool) or not isinstance(val, (int, float)):
            raise ValueError(
                f"Pool {pool_label}: radical_qssa_unzip {block_name}.{key}="
                f"{val!r} must be a number.")
        val = float(val)
        if not math.isfinite(val):
            raise ValueError(
                f"Pool {pool_label}: radical_qssa_unzip {block_name}.{key}="
                f"{val!r} is not finite (NaN/inf are rejected).")
        out[key] = val
    if out["A"] <= 0.0:
        raise ValueError(
            f"Pool {pool_label}: radical_qssa_unzip {block_name}.A={out['A']:g} "
            f"must be > 0 (a non-positive pre-exponential is not a valid rate "
            f"constant).")
    if out["Ea"] < 0.0:
        raise ValueError(
            f"Pool {pool_label}: radical_qssa_unzip {block_name}.Ea="
            f"{out['Ea']:g} must be >= 0 [J/mol].")
    return out


def validate_radical_qssa_unzip(pool_label, channel):
    """Validate a radical_qssa_unzip channel config dict and return it
    normalized (defaults filled, all numerics coerced to float).

    This is the single source of truth for the channel's FIELD rules, shared
    by the deck helper (rmgpy/rmg/input.py, re-raised as InputError),
    PolymerPool.to_config (rmgpy/rmg/polymer_input.py) and the solver's own
    validate_configuration (last line of defense: a directly-constructed
    PolymerPoolConfig bypasses both upstream layers). The per-layer CROSS
    invariants (resolvable monomer routing; mutual exclusion with k_unzip > 0)
    stay with the callers, whose field names differ.

    Channel contract (M1; the QSSA rate law lands in M2 -- until then the
    stored config is inert, nothing in the residual reads it):

    - ``initiation``, ``depropagation``, ``termination``: mandatory Arrhenius
      triplets ``{A, n, Ea}``. SI units BY CONVENTION, not dimensionally
      enforced: A [s^-1] for the unimolecular blocks (initiation,
      depropagation), A [m^3 mol^-1 s^-1] for the bimolecular termination
      block; Ea [J/mol]; n dimensionless. All values finite; A > 0; Ea >= 0.
    - ``transfer``: optional Arrhenius triplet (default None), same
      finite/positivity rules. Units A [s^-1]: the M2 rate law is
      PSEUDO-FIRST-ORDER in the active end R (ktr multiplies R directly in
      the balance ktr*R + 2*kt*R^2 = G_R), so a literature bimolecular
      k_tr [L mol^-1 s^-1 / m^3 mol^-1 s^-1] must be premultiplied by the
      relevant substrate concentration [mol/m^3, SI] BEFORE entering this
      config -- do not drop that concentration factor silently.
    - ``efficiency`` (f_i) and ``monomer_yield`` (y_m): floats in (0, 1],
      default 1.0.
    - ``basis``: must equal RADICAL_QSSA_UNZIP_BASIS (default; forward-compat
      pin).

    Weak-link allyl/U-state vocabulary (schema-2.2 milestone i; OPTIONAL as
    a GROUP, all-or-nothing within it):

    - ``initiation_allyl``: Arrhenius triplet (A [s^-1]), same rules as
      ``initiation`` -- the weak-link (allylic) initiation channel driven by
      the pool's unsaturated-tail-ends state U.
    - ``termination_recombination`` and ``termination_disproportionation``:
      Arrhenius triplets (A [m^3 mol^-1 s^-1]), same rules as the legacy
      summed ``termination``.
    - ``unsaturated_tail_ends_initial``: finite float >= 0, the initial U
      amount [mol] -- the SAME amount basis as mu0 (consumers divide by
      V_poly).
    - MUTUAL EXCLUSION: if ANY weak-link key is present, ALL four must be
      present and the legacy summed ``termination`` must be ABSENT. U is
      produced by the disproportionation branch specifically, so a summed
      kt cannot source U. A channel with NONE of the weak-link keys
      normalizes EXACTLY as before this vocabulary existed (legacy freeze).

    Raises ValueError naming the pool on any violation.
    """
    if not isinstance(channel, dict):
        raise ValueError(
            f"Pool {pool_label}: radical_qssa_unzip must be a dict, got "
            f"{type(channel).__name__}.")
    unknown = sorted(set(channel) - _RADICAL_QSSA_ALLOWED_KEYS)
    if unknown:
        raise ValueError(
            f"Pool {pool_label}: radical_qssa_unzip has unknown key(s) "
            f"{unknown}; allowed keys are {sorted(_RADICAL_QSSA_ALLOWED_KEYS)}.")

    weaklink_present = [k for k in _RADICAL_QSSA_WEAKLINK_KEYS
                        if k in channel]
    if weaklink_present:
        if "termination" in channel:
            raise ValueError(
                f"Pool {pool_label}: radical_qssa_unzip carries the legacy "
                f"SUMMED 'termination' block together with weak-link key(s) "
                f"{weaklink_present}. They are mutually exclusive: the "
                f"weak-link U-state channel requires the SPLIT termination "
                f"blocks (termination_recombination + "
                f"termination_disproportionation) because U production is "
                f"sourced by the disproportionation branch specifically -- "
                f"a summed termination cannot source U. Remove 'termination' "
                f"(split it into the two branches) or remove all weak-link "
                f"keys.")
        missing_weaklink = [k for k in _RADICAL_QSSA_WEAKLINK_KEYS
                            if k not in channel]
        if missing_weaklink:
            raise ValueError(
                f"Pool {pool_label}: radical_qssa_unzip weak-link channel is "
                f"incomplete: key(s) {weaklink_present} present but "
                f"{missing_weaklink} missing. The weak-link vocabulary is "
                f"all-or-nothing: initiation_allyl, "
                f"termination_recombination, termination_disproportionation "
                f"and unsaturated_tail_ends_initial must ALL be present "
                f"(and the legacy summed 'termination' absent).")
        mandatory_blocks = (("initiation", "depropagation")
                            + _RADICAL_QSSA_WEAKLINK_TRIPLET_BLOCKS)
    else:
        mandatory_blocks = _RADICAL_QSSA_MANDATORY_BLOCKS

    normalized = {}
    for block in mandatory_blocks:
        if block not in channel:
            raise ValueError(
                f"Pool {pool_label}: radical_qssa_unzip is missing the "
                f"mandatory Arrhenius block '{block}' (required blocks: "
                f"{list(mandatory_blocks)}).")
        normalized[block] = _validate_qssa_arrhenius_triplet(
            pool_label, block, channel[block])

    if weaklink_present:
        u0 = channel["unsaturated_tail_ends_initial"]
        if isinstance(u0, bool) or not isinstance(u0, (int, float)):
            raise ValueError(
                f"Pool {pool_label}: radical_qssa_unzip "
                f"unsaturated_tail_ends_initial={u0!r} must be a number "
                f"[mol; same amount basis as mu0].")
        u0 = float(u0)
        if not math.isfinite(u0):
            raise ValueError(
                f"Pool {pool_label}: radical_qssa_unzip "
                f"unsaturated_tail_ends_initial={u0!r} is not finite "
                f"(NaN/inf are rejected).")
        if u0 < 0.0:
            raise ValueError(
                f"Pool {pool_label}: radical_qssa_unzip "
                f"unsaturated_tail_ends_initial={u0:g} must be >= 0 "
                f"[mol; same amount basis as mu0].")
        normalized["unsaturated_tail_ends_initial"] = u0

    transfer = channel.get("transfer", None)
    normalized["transfer"] = (
        None if transfer is None
        else _validate_qssa_arrhenius_triplet(pool_label, "transfer", transfer))

    for name in ("efficiency", "monomer_yield"):
        val = channel.get(name, 1.0)
        if isinstance(val, bool) or not isinstance(val, (int, float)):
            raise ValueError(
                f"Pool {pool_label}: radical_qssa_unzip {name}={val!r} must "
                f"be a number.")
        val = float(val)
        if not math.isfinite(val) or not (0.0 < val <= 1.0):
            raise ValueError(
                f"Pool {pool_label}: radical_qssa_unzip {name}={val!r} must "
                f"be a finite float in (0, 1].")
        normalized[name] = val

    basis = channel.get("basis", RADICAL_QSSA_UNZIP_BASIS)
    if basis != RADICAL_QSSA_UNZIP_BASIS:
        raise ValueError(
            f"Pool {pool_label}: radical_qssa_unzip basis={basis!r} is not "
            f"supported; the only allowed basis is "
            f"'{RADICAL_QSSA_UNZIP_BASIS}' (forward-compat pin: other bases "
            f"are rejected until a rate law that consumes them ships).")
    normalized["basis"] = basis

    return normalized


# Ratified label convention (radical-homolysis conduit, adjudicated round 66;
# wording per round 67 ruling a) for the two end-radical daughter pools
# spawned by a k_homolysis parent. The contract is POSITIONAL: 'primary' in
# the suffix means the open-*1 end radical, 'secondary' the open-*2 end
# radical of the backbone C-C cut. The names are NOT a claim of
# primary/secondary radical character -- that reading holds only for PP under
# its head-to-tail repeat convention, is orientation-dependent, and is false
# for other polymers.
K_HOMOLYSIS_DAUGHTER_SUFFIXES = ("_rad_primary_end", "_rad_secondary_end")


def _validate_kernel_arrhenius_triplet(pool_label, param_name, triplet):
    """Shared field validator for pool-level kernel Arrhenius triplets
    {A, n, Ea} (k_homolysis, k_depropagation); returns the triplet
    normalized to plain floats. Error text stays byte-compatible with the
    historical k_homolysis messages (parameter name interpolated)."""
    if not isinstance(triplet, dict):
        raise ValueError(
            f"Pool {pool_label}: {param_name} must be a dict {{A, n, Ea}} "
            f"(Arrhenius triplet, SI convention: A [s^-1], Ea [J/mol]), got "
            f"{type(triplet).__name__}. A bare scalar is rejected by design "
            f"(round 66): the kernel must evaluate k(T) at the solver's "
            f"runtime temperature.")
    missing = [k for k in _RADICAL_QSSA_ARRHENIUS_KEYS if k not in triplet]
    extra = sorted(set(triplet) - set(_RADICAL_QSSA_ARRHENIUS_KEYS))
    if missing or extra:
        raise ValueError(
            f"Pool {pool_label}: {param_name} must have exactly the keys "
            f"{{A, n, Ea}}; missing {missing or 'none'}, unknown "
            f"{extra or 'none'}.")
    out = {}
    for key in _RADICAL_QSSA_ARRHENIUS_KEYS:
        val = triplet[key]
        if isinstance(val, bool) or not isinstance(val, (int, float)):
            raise ValueError(
                f"Pool {pool_label}: {param_name} {key}={val!r} must be a "
                f"number.")
        val = float(val)
        if not math.isfinite(val):
            raise ValueError(
                f"Pool {pool_label}: {param_name} {key}={val!r} is not finite "
                f"(NaN/inf are rejected).")
        out[key] = val
    if out["A"] <= 0.0:
        raise ValueError(
            f"Pool {pool_label}: {param_name} A={out['A']:g} must be > 0 (a "
            f"non-positive pre-exponential is not a valid rate constant).")
    if out["Ea"] < 0.0:
        raise ValueError(
            f"Pool {pool_label}: {param_name} Ea={out['Ea']:g} must be >= 0 "
            f"[J/mol].")
    return out


def validate_k_homolysis(pool_label, triplet):
    """Validate a k_homolysis Arrhenius triplet {A, n, Ea} and return it
    normalized to plain floats (single source of truth for the field rules,
    shared by the deck helper rmgpy/rmg/input.py -- re-raised as InputError
    -- PolymerPool.to_config and the solver's own validate_configuration).

    Same rules and SI convention as the radical_qssa_unzip Arrhenius blocks
    (deliberately NOT a bare scalar -- round 66: 'a scalar precomputed at
    1100 K will poison any ramp/TA replay'; the solver evaluates
    k(T) = A*T^n*exp(-Ea/(R_gas*T)) at its runtime T): exactly the keys
    {A, n, Ea}; every value a finite real number; A > 0 [s^-1]; Ea >= 0
    [J/mol]; n dimensionless. Units are convention, not dimensionally
    enforced (matches the QSSA blocks' posture).
    """
    return _validate_kernel_arrhenius_triplet(pool_label, "k_homolysis",
                                              triplet)


def validate_k_depropagation(pool_label, triplet):
    """Validate a k_depropagation Arrhenius triplet {A, n, Ea} and return it
    normalized to plain floats (single source of truth for the field rules,
    shared by the deck helper rmgpy/rmg/input.py -- re-raised as InputError
    -- derive_daughter_pool_configs and the solver's own
    validate_configuration).

    End-radical DEPROPAGATION kernel (adjudicated round 74 SS2, the run-6
    no-outlet wall fix): per-chain-end unzip frequency
    k_dep(T) = A*T^n*exp(-Ea/(R_gas*T)) evaluated at the solver's RUNTIME
    temperature (same not-a-bare-scalar posture as k_homolysis, round 66).
    A [s^-1], Ea [J/mol], n dimensionless; units are convention, not
    dimensionally enforced."""
    return _validate_kernel_arrhenius_triplet(pool_label, "k_depropagation",
                                              triplet)


# Side-group homolysis initiation kernel (FR1-K1, adjudicated adversarial
# round 70): pool-level side-group X-loss (e.g. aliphatic C-Br -> Br.(gas) +
# mid-chain backbone radical; the chain length is UNCHANGED -- this is NOT a
# chain cut, so it is deliberately a distinct kernel from k_homolysis, and
# the S2 H-loss conduit is NOT reused: the ruling requires its own machinery
# and fingerprint class). One X-loss feature pool per (parent, channel),
# labeled '{parent}_sidegrp_{sanitized channel label}'.
SIDE_GROUP_DAUGHTER_INFIX = "_sidegrp_"

_SIDE_GROUP_CHANNEL_KEYS = ("label", "A", "n", "Ea", "site_selector",
                            "sites_per_unit", "gas_product")

# Structural site-selector vocabulary (round-72 P1 adjudication): every
# side_group_homolysis channel MUST name the carbon environment of the X
# atom it removes -- the kinetics label alone must never pick the site.
# Implemented as explicit in-repo atom-environment predicates over RMG's
# own molecular graph (adjudication-equivalent alternative to site_smarts:
# full SMARTS would require an RMG<->RDKit atom-index round-trip whose
# ordering is not contract-stable -- a misordered mapping would silently
# select the wrong atom, the exact bug class being fixed). The three
# classes PARTITION the carbon-bound X sites (mutually exclusive, jointly
# exhaustive), so two channels can only resolve to the SAME atom set by
# duplicating a (selector, element) pair -- which the validator hard-fails
# as a double-carry.
SIDE_GROUP_SITE_SELECTORS = ("aryl", "benzylic", "aliphatic")


def side_group_site_atom_indices(monomer, element_symbol, site_selector):
    """Deterministic structural site selection for the side_group_homolysis
    kernel (round-72 P1): the atom indices, in ``monomer.atoms`` order (the
    canonical deterministic order; the producer removes the FIRST match),
    of the removable side-group atoms of ``element_symbol`` whose carbon
    environment matches ``site_selector``.

    Removable side-group X atom (the same predicate the producer's removal
    always used): element match, unlabeled (never a stitch atom),
    closed-shell, terminal (exactly one bond), single-bonded. The selector
    then classifies the carbon the X hangs off (an X on a non-carbon
    neighbor matches NO selector -- the v1 vocabulary is carbon-environment
    only, consistent with the kernel's C-X homolysis law):

    - 'aryl':      the neighbor carbon is in an aromatic ring
                   (Molecule.get_aromatic_rings -- robust to kekulized vs
                   aromatic bond representations).
    - 'benzylic':  the neighbor carbon is NOT aromatic but is bonded to an
                   aromatic-ring carbon.
    - 'aliphatic': the neighbor carbon is neither aromatic nor benzylic.
    """
    aromatic_atoms = set()
    for ring in monomer.get_aromatic_rings()[0]:
        aromatic_atoms.update(ring)
    matches = []
    for idx, atom in enumerate(monomer.atoms):
        if (atom.symbol != element_symbol or atom.label
                or atom.radical_electrons != 0 or len(atom.bonds) != 1
                or not next(iter(atom.bonds.values())).is_single()):
            continue
        neighbor = next(iter(atom.bonds.keys()))
        if not neighbor.is_carbon():
            continue
        if neighbor in aromatic_atoms:
            cls = "aryl"
        elif any(nb in aromatic_atoms and nb.is_carbon()
                 for nb in neighbor.bonds.keys()):
            cls = "benzylic"
        else:
            cls = "aliphatic"
        if cls == site_selector:
            matches.append(idx)
    return tuple(matches)


def sanitize_side_group_channel_label(channel_label):
    """Deterministic sanitizer for a side_group_homolysis channel label used
    inside pool labels and fingerprints: every character outside [A-Za-z0-9_]
    becomes '_' (e.g. 'aliphatic_C-Br' -> 'aliphatic_C_Br'). Single source of
    truth shared by the validator (sanitized-collision duplicate check), the
    producer (rmgpy.polymer.Polymer.generate_side_loss_daughters) and the
    flattener's destination-pool resolution."""
    return "".join(c if (c.isalnum() or c == "_") else "_"
                   for c in channel_label)


def side_group_daughter_pool_label(pool_label, channel_label):
    """The ratified label of the X-loss feature pool spawned for one
    (parent, channel) pair: '{parent}_sidegrp_{sanitized channel label}'."""
    return (f"{pool_label}{SIDE_GROUP_DAUGHTER_INFIX}"
            f"{sanitize_side_group_channel_label(channel_label)}")


def _side_group_gas_mw_g_mol(gas_product):
    """Molar mass [g/mol] of a channel's gas_product SMILES (M_X of the
    ejected X radical). Lazy Molecule import (cycle-avoidance idiom)."""
    from rmgpy.molecule import Molecule
    return Molecule().from_smiles(gas_product).get_molecular_weight() * 1000.0


def validate_side_group_homolysis(pool_label, channels, monomer=None):
    """Validate a side_group_homolysis channel LIST and return it normalized
    (single source of truth for the field rules, shared by the deck helper
    rmgpy/rmg/input.py -- re-raised as InputError -- PolymerPool.to_config,
    the producer and the solver's validate_configuration; mirrors
    validate_k_homolysis).

    Each channel is a dict with EXACTLY the keys
    {label, A, n, Ea, site_selector, sites_per_unit, gas_product}:
    - label: non-empty str naming the bond class (e.g. 'aliphatic_C-Br').
      Channels stay SEPARATE, never lumped (round-70 ruling); duplicate
      labels -- raw or after sanitization -- are rejected.
    - A > 0 [s^-1 per site], n dimensionless, Ea >= 0 [J/mol]: same SI
      convention and finite-rejection posture as validate_k_homolysis (the
      solver evaluates k(T) = A*T^n*exp(-Ea/(R_gas*T)) at its runtime T).
    - site_selector: REQUIRED structural site selector (round-72 P1), one
      of SIDE_GROUP_SITE_SELECTORS ('aryl' | 'benzylic' | 'aliphatic') --
      names the carbon environment of the X atom the channel removes. The
      kinetics label alone must NEVER pick the site: on a mixed-site
      monomer a channel labeled 'aliphatic_C-Br' would otherwise remove
      the first X in atom order (an ARYL one), minting the wrong -- or the
      same -- defect structure under a different rate.
    - sites_per_unit > 0: side-group X sites per monomer repeat unit
      (e.g. aryl C-Br = 4, benzylic C-Br = 2). CHECKED, not trusted, when
      ``monomer`` is given (below).
    - gas_product: SMILES of the ejected X radical. Must parse and be a
      MONO-RADICAL; v1 supports only MONOATOMIC radicals (e.g. '[Br]') --
      the deterministic feature-unit derivation removes a single X atom.

    STRUCTURAL layer (round-72 P1): when ``monomer`` -- the pool's parsed
    repeat-unit Molecule -- is given, each channel's selector is resolved
    against it via side_group_site_atom_indices and the law is enforced:
    (a) the selector must match >= 1 removable X atom of the channel's
    gas_product element (no match -> hard fail naming the channel);
    (b) sites_per_unit must EQUAL the selector's structural match count
    (mismatch -> hard fail with both numbers); (c) no two channels may
    resolve to the SAME atom set (that is a double-carry the distinct
    labels were hiding -> hard fail naming both channels). The layers that
    hold the monomer structure (deck helper, PolymerPool.to_config, the
    producer) pass it; the solver's directly-constructed PolymerPoolConfig
    backstop carries no structure and validates shape + vocabulary only.
    """
    if isinstance(channels, dict) or not isinstance(channels, (list, tuple)):
        raise ValueError(
            f"Pool {pool_label}: side_group_homolysis must be a list of "
            f"channel dicts (one entry per bond class, e.g. "
            f"[{{'label': 'aliphatic_C-Br', 'A': ..., 'n': ..., 'Ea': ..., "
            f"'sites_per_unit': ..., 'gas_product': '[Br]'}}]), got "
            f"{type(channels).__name__}.")
    if len(channels) == 0:
        raise ValueError(
            f"Pool {pool_label}: side_group_homolysis must not be an empty "
            f"list -- omit the parameter entirely to disable the kernel "
            f"(an empty list would be a silently inert channel).")
    normalized = []
    seen_sanitized = {}
    for pos, ch in enumerate(channels):
        if not isinstance(ch, dict):
            raise ValueError(
                f"Pool {pool_label}: side_group_homolysis entry {pos} must "
                f"be a channel dict, got {type(ch).__name__}.")
        missing = [k for k in _SIDE_GROUP_CHANNEL_KEYS if k not in ch]
        extra = sorted(set(ch) - set(_SIDE_GROUP_CHANNEL_KEYS))
        if missing or extra:
            raise ValueError(
                f"Pool {pool_label}: side_group_homolysis entry {pos} "
                f"(label={ch.get('label', '<unset>')!r}) must have exactly "
                f"the keys {{label, A, n, Ea, site_selector, sites_per_unit, "
                f"gas_product}}; missing {missing or 'none'}, unknown "
                f"{extra or 'none'}.")
        label = ch["label"]
        if not isinstance(label, str) or not label.strip():
            raise ValueError(
                f"Pool {pool_label}: side_group_homolysis entry {pos} label "
                f"must be a non-empty string naming the bond class, got "
                f"{label!r}.")
        out = {"label": label}
        sel = ch["site_selector"]
        if not isinstance(sel, str) or sel not in SIDE_GROUP_SITE_SELECTORS:
            raise ValueError(
                f"Pool {pool_label}: side_group_homolysis channel '{label}' "
                f"site_selector={sel!r} must be one of "
                f"{SIDE_GROUP_SITE_SELECTORS} -- the REQUIRED structural "
                f"site selector pins WHICH X atom the channel removes "
                f"(round-72: a kinetics label alone must never pick the "
                f"site).")
        out["site_selector"] = sel
        for key in ("A", "n", "Ea", "sites_per_unit"):
            val = ch[key]
            if isinstance(val, bool) or not isinstance(val, (int, float)):
                raise ValueError(
                    f"Pool {pool_label}: side_group_homolysis channel "
                    f"'{label}' {key}={val!r} must be a number.")
            val = float(val)
            if not math.isfinite(val):
                raise ValueError(
                    f"Pool {pool_label}: side_group_homolysis channel "
                    f"'{label}' {key}={val!r} is not finite (NaN/inf are "
                    f"rejected).")
            out[key] = val
        if out["A"] <= 0.0:
            raise ValueError(
                f"Pool {pool_label}: side_group_homolysis channel '{label}' "
                f"A={out['A']:g} must be > 0 (a non-positive pre-exponential "
                f"is not a valid rate constant).")
        if out["Ea"] < 0.0:
            raise ValueError(
                f"Pool {pool_label}: side_group_homolysis channel '{label}' "
                f"Ea={out['Ea']:g} must be >= 0 [J/mol].")
        if out["sites_per_unit"] <= 0.0:
            raise ValueError(
                f"Pool {pool_label}: side_group_homolysis channel '{label}' "
                f"sites_per_unit={out['sites_per_unit']:g} must be > 0 "
                f"(side-group X sites per monomer repeat unit).")
        gp = ch["gas_product"]
        if not isinstance(gp, str) or not gp.strip():
            raise ValueError(
                f"Pool {pool_label}: side_group_homolysis channel '{label}' "
                f"gas_product must be a SMILES string, got {gp!r}.")
        from rmgpy.molecule import Molecule
        try:
            gmol = Molecule().from_smiles(gp)
        except Exception as e:
            raise ValueError(
                f"Pool {pool_label}: side_group_homolysis channel '{label}' "
                f"gas_product={gp!r} does not parse as SMILES ({e}).")
        if gmol.get_radical_count() != 1:
            raise ValueError(
                f"Pool {pool_label}: side_group_homolysis channel '{label}' "
                f"gas_product={gp!r} must be a mono-radical (the homolysis "
                f"co-product is the ejected X radical); got radical count "
                f"{gmol.get_radical_count()}.")
        if len(gmol.atoms) != 1:
            raise ValueError(
                f"Pool {pool_label}: side_group_homolysis channel '{label}' "
                f"gas_product={gp!r} has {len(gmol.atoms)} atoms -- v1 "
                f"supports only monoatomic radical gas products (e.g. "
                f"'[Br]'): the deterministic feature-unit derivation removes "
                f"a single X atom from the repeat unit.")
        out["gas_product"] = gp
        san = sanitize_side_group_channel_label(label)
        if san in seen_sanitized:
            raise ValueError(
                f"Pool {pool_label}: side_group_homolysis has duplicate "
                f"channel labels: {seen_sanitized[san]!r} and {label!r} "
                f"collide (sanitized '{san}'). Two decks of the same channel "
                f"would double-carry the same bond class; channels must be "
                f"unique.")
        seen_sanitized[san] = label
        normalized.append(out)

    # STRUCTURAL layer (round-72 P1) -- only when the caller holds the
    # parsed repeat-unit structure; see the docstring for the three laws.
    if monomer is not None:
        from rmgpy.molecule import Molecule
        seen_sites = {}
        for out in normalized:
            label = out["label"]
            sym = Molecule().from_smiles(out["gas_product"]).atoms[0].symbol
            idxs = side_group_site_atom_indices(monomer, sym,
                                                out["site_selector"])
            if not idxs:
                raise ValueError(
                    f"Pool {pool_label}: side_group_homolysis channel "
                    f"'{label}' site_selector='{out['site_selector']}' "
                    f"matches NO removable side-group {sym} atom (terminal, "
                    f"single-bonded, closed-shell, unlabeled, on a matching "
                    f"carbon environment) in the repeat unit -- the "
                    f"selector must structurally pin which {sym} the "
                    f"channel removes; a non-matching selector would mint "
                    f"the wrong defect (round-72 P1).")
            if float(len(idxs)) != out["sites_per_unit"]:
                raise ValueError(
                    f"Pool {pool_label}: side_group_homolysis channel "
                    f"'{label}' sites_per_unit={out['sites_per_unit']:g} "
                    f"contradicts the structural match count: "
                    f"site_selector='{out['site_selector']}' matches "
                    f"{len(idxs)} {sym} site(s) in the repeat unit. "
                    f"sites_per_unit is CHECKED against the monomer, never "
                    f"trusted (round-72 P1); set it to {len(idxs)} or fix "
                    f"the selector.")
            key = (sym, frozenset(idxs))
            if key in seen_sites:
                raise ValueError(
                    f"Pool {pool_label}: side_group_homolysis channels "
                    f"'{seen_sites[key]}' and '{label}' resolve to the SAME "
                    f"{sym} atom set {sorted(idxs)} in the repeat unit -- "
                    f"two rate channels on one structural site double-carry "
                    f"the loss (the distinct labels were hiding it; "
                    f"round-72 P1). Merge them or pick disjoint selectors.")
            seen_sites[key] = label
    return normalized


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

    # Monomer (repeat-unit) heavy-atom (non-H) count. The SECOND axis of the
    # r89 dual-axis polymer-sized melt gate in _reference_state_tripwire
    # (mirror of rmgpy.polymer._discrete_is_polymer_sized: mass AND heavy
    # structure, both computable, both at/above
    # _IMPOSTOR_DISCRETE_MONOMER_UNITS monomer-equivalents against >= 1
    # pool). Populated by every config producer from the monomer STRUCTURE
    # (rmgpy/rmg/polymer_input.py to_pool_config /
    # derive_daughter_pool_configs; consumer world reconstructs it from the
    # artifact's monomer_adj_list) -- per the r89 data-flow constraint the
    # heavy axis is never approximated from mass. Default 0 = axis
    # UNCOMPUTABLE (honest degradation: tag-branch candidates classify
    # conservative-gas and the tripwire emits an axis-undecidable census
    # warning; nothing is classified blind).
    monomer_heavy_atoms: int = 0

    # Kinetics Parameters
    k_scission: float = 0.0

    # Recession Rate [1/s]
    # REQUIRED for the Hybrid Handshake (Tail -> Explicit flux).
    # This parameter sets the timescale for both the default recession chemistry
    # (n -> n-1) AND the physical handshake flux across the boundary xs.
    k_unzip: float = 0.0

    # Radical QSSA unzip channel (M1: config + validation ONLY -- nothing in
    # the residual reads this yet; the QSSA rate law lands in M2). Normalized
    # dict per validate_radical_qssa_unzip: {initiation, depropagation,
    # termination: {A, n, Ea}, transfer: {A, n, Ea} | None, efficiency,
    # monomer_yield: float in (0, 1], basis: RADICAL_QSSA_UNZIP_BASIS}.
    # Units BY CONVENTION (documented, not dimensionally enforced):
    # initiation A [s^-1], depropagation A [s^-1], termination A
    # [m^3 mol^-1 s^-1] (the only bimolecular block), transfer A [s^-1]
    # (PSEUDO-first-order: ktr multiplies R directly in the M2 rate law --
    # premultiply a literature bimolecular k_tr by the substrate
    # concentration [mol/m^3] before configuring it); Ea [J/mol].
    # Mutually exclusive with k_unzip > 0
    # (double-count guard) and requires monomer_poly_index (the channel
    # reuses the pool's existing monomer routing -- no new routing field).
    radical_qssa_unzip: Optional[Dict[str, object]] = None

    # Radical-homolysis initiation kernel (Stage 1, adjudicated round 66).
    # Normalized Arrhenius triplet {A, n, Ea} per validate_k_homolysis
    # (SI convention: A [s^-1], Ea [J/mol]; NOT a bare scalar -- the RHS
    # evaluates k(T) = A*T^n*exp(-Ea/(R_gas*T)) at the runtime T). Event
    # rate R = k(T)*max(mu1-mu0, 0) [backbone bonds]; the parent pool is
    # debited the bond-weighted breaking chain and BOTH end-radical daughter
    # pools ({label}_rad_primary_end / {label}_rad_secondary_end, created at
    # model setup by the producer) are credited the uniform-cut fragments.
    # One-way (no reverse leg); releases NO gas species. Mutually exclusive
    # with k_scission > 0 (same random backbone-break event; double-count),
    # with radical_qssa_unzip (QSSA initiation IS random backbone
    # homolysis; double-count), and (round 67) with k_unzip > 0 (legacy
    # closed-chain monomer-loss channel vs radical-end pools feeding
    # explicit beta-scission/unzip chemistry; double-carried
    # depolymerization).
    k_homolysis: Optional[Dict[str, float]] = None

    # End-radical DEPROPAGATION kernel (adjudicated round 74 SS2, the run-6
    # no-outlet wall fix): normalized Arrhenius triplet {A, n, Ea} per
    # validate_k_depropagation (SI convention: A [s^-1], Ea [J/mol]; NOT a
    # bare scalar -- the RHS evaluates k_dep(T) = A*T^n*exp(-Ea/(R_gas*T))
    # at the runtime T) or None. Consumption channel for RADICAL-END pools
    # (one active radical end per chain): per unzip event one monomer
    # volatile is released to monomer_poly_index (gas amount basis).
    #   R    = k_dep(T) * mu0 * g     [g = smooth exhaustion gate, == 1 in
    #                                  the realizable region mean DP >= 1]
    #   dmu1 = -R;  gas = +R (the SAME float -- exact mass invariant)
    #   dmu2 = -k_dep*(2*mu1 - mu0)   [smooth positive-part on the excess]
    #   dmu0 = -k_dep*N1, N1 = mu0 * P(DP=1) from the EXISTING
    #          discrete/gamma closure (never a permanent dmu0 = 0 -- r74:
    #          that stalls at a one-repeat-per-chain residue)
    # Mutually exclusive on one pool with k_unzip > 0 (legacy
    # phenomenological form of the SAME chain-end monomer-release event),
    # with radical_qssa_unzip (its depropagation block IS lumped chain-end
    # depropagation) and with k_homolysis (multi-generation homolysis of
    # radical-ended chains is DEFERRED, r74 SS3; and a closed-chain
    # initiation pool has no radical end to depropagate). Requires a
    # resolvable monomer_poly_index (released units must land somewhere).
    k_depropagation: Optional[Dict[str, float]] = None

    # Side-group homolysis initiation kernel (FR1-K1, adjudicated round 70):
    # normalized channel LIST per validate_side_group_homolysis (each entry
    # exactly {label, A, n, Ea, site_selector, sites_per_unit, gas_product};
    # site_selector is the round-72 REQUIRED structural site selector) or
    # None. Per
    # channel with s = sites_per_unit and k(T) = A*T^n*exp(-Ea/(R_gas*T))
    # per site: total sites = s*mu1, so the event rate is R = k(T)*s*mu1 and
    # the reacting chain is picked with probability ~ its length n
    # (site-weighted). Parent debit dmu_j -= k*s*mu_{j+1} (j = 0, 1, 2; mu3
    # from the log-Lagrange closure); the '{label}_sidegrp_{channel}' X-loss
    # feature pool is credited EXACTLY the parent debit (the chain transfers
    # INTACT -- no chain cut); the gas X radical is credited +R. Mutually
    # exclusive with k_unzip > 0 and with radical_qssa_unzip (double-carry);
    # MAY coexist with k_homolysis / k_scission (different bonds: side-group
    # C-X vs backbone C-C). v1 LIMITATION: the kernel acts on the PARENT
    # pool only -- feature pools carry no further side-group loss (they
    # saturate as terminal X-loss sinks; no multi-loss cascade).
    side_group_homolysis: Optional[List[Dict[str, object]]] = None
    # Resolved core index of each channel's gas_product species, aligned
    # with side_group_homolysis (PolymerPool.to_config resolves it; a
    # kernel-enabled pool without resolved indices is a hard config error --
    # the emitted X radical would silently vanish).
    side_group_gas_indices: Optional[Tuple[int, ...]] = None

    # Exact per-chain mass defect [g/mol] of an X-loss feature pool (FR1-K1
    # mass contract, the round-70 #1 P1 trap): the feature pool keeps the
    # parent's monomer_mw_g_mol (the chain transferred INTACT), but every
    # feature chain lost exactly ONE X atom (v1: no multi-loss cascade), so
    # the pool's EXACT condensed mass is
    #     mass [g] = mu1*monomer_mw_g_mol - mu0*chain_mass_defect_g_mol
    # (see condensed_mass_g). The producer pins it to M_X of the spawning
    # channel's gas_product; _flatten_side_group_state hard-errors on any
    # kernel destination whose defect does not pin M_X exactly. 0.0 on
    # ordinary pools (the accessor then reduces to mu1*MW).
    chain_mass_defect_g_mol: float = 0.0

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

    def condensed_mass_g(self, mu0_mol, mu1_mol):
        """EXACT condensed mass [g] of this pool at the given extensive
        moments [mol] (FR1-K1 mass contract): mu1*monomer_mw_g_mol minus the
        per-chain X-loss defect mu0*chain_mass_defect_g_mol. For ordinary
        pools (defect 0) this is the plain moment-derived mass; for X-loss
        feature pools it carries the one-lost-X-per-chain correction so that
        d(condensed mass)/dt + d(gas X mass)/dt = 0 holds exactly on the
        side-group kernel's contributions (pinned by test). Also works on
        rates: feeding d(mu)/dt returns d(mass)/dt."""
        return (mu1_mol * self.monomer_mw_g_mol
                - mu0_mol * self.chain_mass_defect_g_mol)

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


# End-radical DEPROPAGATION kernel (adjudicated round 74 SS2) smooth-gate
# width, dimensionless in mean DP. sp(x) = x^3/(x^2 + W^2) is a C2 smooth
# positive-part: exactly 0 for x <= 0, asymptotically x - W^2/x for x >> W.
# The gate g = 1 - sp(1 - mean) is therefore EXACTLY 1 in the realizable
# region mean >= 1 (the healthy law is exact, no perturbation) and rolls
# off C2-smoothly below it (r74 SS5: no hard max(...,0) cliff at
# exhaustion -- the DASPK grind / IDID=-7 class). Residual monomer release
# from a pathological mu1 = 0, mu0 > 0 noise state is bounded by
# k_dep*mu0*W^2 (documented honest degradation; the chain count keeps
# draining at the full -k_dep*mu0 there, so the state cannot become a
# stiff no-outlet grind).
KDEP_GATE_WIDTH = 1.0e-2


def _smooth_pos(x: float, w: float) -> float:
    """C2 smooth positive-part x^3/(x^2 + w^2): 0 for x <= 0, ~x for
    x >> w, with value/first/second derivative all continuous at 0
    (strictly smoother than the max(x, 0) clamps elsewhere)."""
    if x <= 0.0:
        return 0.0
    return x * x * x / (x * x + w * w)


def _deprop_dp1_fraction(mu0: float, mu1: float, mu2: float) -> float:
    """DP=1 chain fraction p1 (N1 = mu0 * p1) for the end-radical
    depropagation kernel (r74 SS2: 'Estimate N1 using the existing
    discrete/gamma-style moment closure').

    Gamma leg -- the EXISTING closure machinery: shape/scale from
    _gamma_params_from_mu012 (k = 1/(PDI-1), theta = mean/k), DP=1 bin
    probability from _gamma_prob_conditional_hybrid(1, 0, k, theta), i.e.
    [F(1.5) - F(0.5)] / [1 - F(0.5)] on the half-integer-bin discretization
    conditioned on n >= 1 (the same convention the hybrid handshake uses).

    Smooth terminal boundary floor (disclosed per r74): the gamma leg is
    unavailable for degenerate/monodisperse shapes (PDI <= 1 + 1e-6 returns
    no params), so p1 is floored by a C1 smoothstep 1 -> 0 over mean DP in
    [1, 2]. The floor is what guarantees NO stall: any state at or below
    one repeat unit per chain has p1 = 1, so dmu0 = -k_dep*mu0 terminates
    the residue instead of freezing it (r74: no permanent dmu0 = 0).

    Returns a value in [0, 1]; smooth (max of two continuous legs) in the
    moment state."""
    if mu0 <= SMALL_EPS:
        return 0.0
    mean = mu1 / mu0
    t = mean - 1.0
    if t <= 0.0:
        p_floor = 1.0
    elif t >= 1.0:
        p_floor = 0.0
    else:
        p_floor = 1.0 - (3.0 * t * t - 2.0 * t * t * t)
    p_gamma = 0.0
    params = _gamma_params_from_mu012(mu0, mu1, mu2)
    if params:
        k_shape, theta = params
        p_gamma = _gamma_prob_conditional_hybrid(1, 0, k_shape, theta)
    p1 = p_gamma if p_gamma > p_floor else p_floor
    if p1 < 0.0:
        return 0.0
    if p1 > 1.0:
        return 1.0
    return p1


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


def _dual_axis_polymer_sized(mw_g_mol, heavy_count, pool_axes, units):
    """r89 dual-axis polymer-sized test at the solver seam -- the EXACT
    mirror of rmgpy.polymer._discrete_is_polymer_sized measured against the
    configured pools, on pre-plumbed numbers instead of structures (the
    solver's PolymerPoolConfig sidecars and consumer-world species carry no
    Molecule):

    * mass axis:  ``mw_g_mol`` [g/mol] vs ``units * pool monomer MW``;
    * heavy axis: ``heavy_count`` (non-H atoms; -1 = unknown, e.g. a
      consumer-world species whose artifact carried no composition) vs
      ``units * pool monomer_heavy_atoms`` (0 = unknown denominator).

    Returns ``(sized, missing)``. ``sized`` is True iff at least one pool
    decides True on BOTH computable axes (any-pool, matching
    _discrete_is_chain_scale_proxy_derived). ``missing`` is None when every
    consulted pool DECIDED (True somewhere, or computable evidence decided
    False everywhere); otherwise it names the uncomputable axis/axes of a
    pool whose computable evidence did not already decide False -- the r85
    P2(a) undecidable case: an uncomputable axis NEVER defers to the other
    (a single-axis melt call is exactly the run-9 false positive in reverse),
    the answer degrades to conservative-gas, and the caller must announce it
    (mirroring rmgpy.polymer._warn_impostor_axis_undecidable semantics)."""
    missing = None
    for monomer_mw, monomer_heavy in pool_axes:
        mass_computable = monomer_mw > 0.0 and mw_g_mol > 0.0
        mass_sized = mass_computable and mw_g_mol >= units * monomer_mw
        heavy_computable = monomer_heavy > 0 and heavy_count >= 0
        heavy_sized = heavy_computable and heavy_count >= units * monomer_heavy
        if mass_computable and heavy_computable:
            if mass_sized and heavy_sized:
                return True, None
            continue  # both axes computable: decided False for this pool
        # At least one axis uncomputable for this pool: undecidable UNLESS a
        # computable axis already decided False (then it is a decided False).
        if not any([mass_computable and not mass_sized,
                    heavy_computable and not heavy_sized]):
            missing = ("mass and structure"
                       if not mass_computable and not heavy_computable
                       else ("mass" if not mass_computable else "structure"))
    return False, missing


def _assert_chain_scale_melt_member(label, mw_kg_mol, heavy_count,
                                    gas_classified, pool_axes, units):
    """Cannot-happen leak guard (spec §5.1 C3 amendment 2026-06-11, r89
    dual-axis form), called for every species entering the melt sum. A
    gas-classified species can be physically-melt ONLY via the tag branch,
    whose dual-axis polymer-sized conjunct (_dual_axis_polymer_sized: mass
    AND heavy axes, both computable, both >= units monomer-equivalents
    against >= 1 pool) is part of the CLASS DEFINITION -- so a gas-classified
    member below the threshold (or with an undecidable axis) cannot reach
    the sum through the amended gate. A tagged below-threshold species being
    EXCLUDED by the gate is expected and silent (the family.py:1657
    over-tagging fingerprint, H2 on every proxy-touching reaction; the PP
    run-9 DP-2 hexadiene at 2.0 monomer-equivalents); the raise here fires
    only if such a species REACHES the sum (a future refactor recomputing
    membership without the conjunct), converting the mistake into a loud
    CLASSIFICATION error instead of a silently-large U and a misattributed
    reference-state refusal. Condensed-branch members (gas_classified=False;
    pool-configured by input) are exempt from the gate and from this
    guard."""
    if not gas_classified:
        return
    sized, _missing = _dual_axis_polymer_sized(
        mw_kg_mol * 1000.0, heavy_count, pool_axes, units)
    if not sized:
        raise ValueError(
            "THERMO REFERENCE-STATE TRIPWIRE (classification leak): the "
            "is_polymer_proxy tag includes a non-chain species in the melt "
            f"sum ({label}, MW = {mw_kg_mol * 1000.0:.2f} g/mol, "
            f"heavy atoms = {heavy_count}, not polymer-sized at "
            f"{units} monomer-equivalents on both axes against any pool); "
            "physically-melt class definition violated -- this is a "
            "classification leak, NOT a thermo problem; do not respond by "
            "touching reference states. See the proxy-tag propagation chain "
            "(family.py) and the invariant section of "
            "docs/multi_pool_design.md.")


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
        prospective_classifier=None,
        prospective_condensed_edge_daughter_classifier=None,
        allow_default_prospective_edge: bool = False,
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
        allow_unstamped_proxy_rows: bool = False,
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
        # Item 17 A5-2 (spec 2026-06-12 SS3(a)/(d)): the §3(a) staging fix.
        # _prospective_classifier is a CALLABLE f(species_list) -> bool array
        # (the bound polymerPhase.get_gas_mask), plumbed onto the engine by
        # to_solver_object so initialize_model can re-run stage-1 over the
        # LIVE chain(core, edge) on EVERY call (base.pyx:simulate already
        # re-runs initialize_model with the live edge, so no reactor-level
        # trigger is needed). Used for CLASSIFICATION ONLY -- R3-clean, never a
        # phase verdict for any core behavior. None on direct-test/runner
        # construction. _allow_default_prospective_edge is the opt-in flag that
        # marks a build as legitimately fallback-permitted (direct-test/runner
        # with no blueprint phase object); R1-EDGE raises iff a build is NOT so
        # flagged AND its prospective edge suffix was default-filled.
        # _prospective_edge_provenance is the per-build int8 marker over the
        # edge entries (1 = stage-1-classified, 0 = default-filled).
        self._prospective_classifier = prospective_classifier
        # Spec 2026-06-29: callable f(combined_species) -> set of qualifying
        # daughter BASE labels (the bound PolymerPhase.get_condensed_edge_daughter_bases),
        # plumbed by to_solver_object. Re-run over the LIVE chain(core, edge) in
        # initialize_model on EVERY call (never a frozen set -> cannot go stale
        # on the engine-reuse path). None on direct-test/runner construction.
        self.prospective_condensed_edge_daughter_classifier = \
            prospective_condensed_edge_daughter_classifier
        self._allow_default_prospective_edge = bool(allow_default_prospective_edge)
        # r71 FIX 4 escape hatch (test/runner-only, mirroring
        # allow_default_prospective_edge): permits the legacy mu1-only
        # NONE -> UNRESOLVED demotion for unstamped LIVE proxy rows instead
        # of the production hard-fail at initialize_model. A production
        # build (rmgpy.rmg construction) leaves this default-False.
        self._allow_unstamped_proxy_rows = bool(allow_unstamped_proxy_rows)
        self._prospective_edge_provenance = None
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
        # Weak-link U-state (milestones ii+iii): number of appended per-pool
        # U slots and their dU/dt scratch. 0/None until initiate_tolerances /
        # validate_configuration lay the state out; legacy-only systems keep
        # num_qssa_u == 0 and the pre-milestone neq == n_core layout.
        self.num_qssa_u = 0
        self._scratch_du_dt = None

        self.pool_mu0_indices = np.full(len(self.polymer_pools), -1, dtype=np.int32)
        self.pool_mu1_indices = np.full(len(polymer_pools), -1, dtype=np.int32)

        # Exhaustion-tail conditioning state (r81 B). None/empty until
        # initialize_model lays out the state vector and the atol-derived
        # per-moment floors; the census set is per-rebuild (cleared there).
        self._pool_exhausted = None
        self._pool_mu_floors = None
        self._exhaustion_census_emitted = set()

    def initiate_tolerances(self, atol=1e-16, rtol=1e-8, sensitivity=False,
                            sens_atol=1e-6, sens_rtol=1e-4):
        """Extend the base ODE layout with one U slot per weak-link pool
        (milestones ii+iii).

        Layout contract: y = [core species amounts (mol) | U slots (mol)],
        one appended slot per pool whose radical_qssa_unzip channel carries
        the weak-link vocabulary, in pool order. Pools with a legacy-only
        channel (or none) get NO slot: the legacy state layout is UNCHANGED
        (neq == n_core, pinned bitwise by the golden RHS test). Slot
        positions are bound in _flatten_radical_qssa_state (qssa_u_slot);
        the count here keys on the same "initiation_allyl" key presence, and
        the flattening pass cross-checks the two counts.

        The DASPK sensitivity layout interleaves per-species blocks
        (neq = n_core*(n_rxn+n_core+1)); a trailing U slot does not fit it,
        so sensitivity + weak-link is refused loudly rather than silently
        corrupting the sensitivity state vector.
        """
        ReactionSystem.initiate_tolerances(self, atol, rtol, sensitivity,
                                           sens_atol, sens_rtol)
        n_weak = 0
        for pool in self.polymer_pools:
            q = pool.radical_qssa_unzip
            if isinstance(q, dict) and "initiation_allyl" in q:
                n_weak += 1
        self.num_qssa_u = n_weak
        if n_weak == 0:
            return
        if sensitivity:
            raise ValueError(
                f"{n_weak} weak-link allyl/U-state pool(s) configured with "
                f"sensitivity analysis enabled: the DASPK sensitivity state "
                f"layout has no room for the appended U slots. Disable "
                f"sensitivity or remove the weak-link channel(s).")
        self.neq += n_weak
        self.atol_array = np.concatenate(
            [self.atol_array, np.full(n_weak, atol, dtype=float)])
        self.rtol_array = np.concatenate(
            [self.rtol_array, np.full(n_weak, rtol, dtype=float)])

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
                # INVERTED contract (incident 2026-07-03, design B-prime): the
                # released-monomer target is the deck's gas volatile and the
                # unzip/QSSA release deposits into the GAS amount basis. A
                # condensed-masked target would re-create the reference-state
                # conflation the tripwire refuses (one Species carrying two
                # phase residences).
                if not self.gas_species_mask[pool.monomer_poly_index]:
                    raise ValueError(
                        f"Pool {pool.label} monomer index masked as CONDENSED. "
                        f"The unzip/QSSA release target (monomer_product) is "
                        f"emitted to the GAS phase; it must be classified gas.")

            # Unzip-channel invariants. This is the LAST line of defense: the
            # deck helper (rmg/input.py), PolymerPool.to_config
            # (rmg/polymer_input.py) and the artifact runner
            # (tools/polymer_moments_runner.py) all guard the same shapes, but
            # a directly-constructed PolymerPoolConfig bypasses every one of
            # them. The residual drains the condensed moments unconditionally
            # when k_unzip > 0 (dmu1_dt -= k_unzip*mu0; dmu2_dt -=
            # k_unzip*(2*mu1 - mu0)) while the released-monomer emission is
            # gated on monomer_poly_index is not None -- so k_unzip > 0 with
            # no monomer index makes mass leave the condensed phase
            # un-conserved (the drained repeat units go nowhere).
            # Finiteness first: NaN passes BOTH the `< 0` and `> 0` gates as
            # False, so a non-finite k_unzip/k_scission would make the
            # channel SILENTLY INERT (a laundered no-op) while inf poisons
            # the residual. Mirror the QSSA triplet validator's posture.
            for _rate_name, _rate_val in (("k_unzip", pool.k_unzip),
                                          ("k_scission", pool.k_scission)):
                if not math.isfinite(_rate_val):
                    raise ValueError(
                        f"Pool {pool.label}: {_rate_name}={_rate_val!r} is not "
                        f"finite -- NaN/inf is not a valid rate constant (NaN "
                        f"would silently disable the channel; inf would poison "
                        f"the residual). Set {_rate_name} to a finite value "
                        f">= 0.")
            if pool.k_unzip < 0.0:
                raise ValueError(
                    f"Pool {pool.label}: k_unzip={pool.k_unzip:g} is negative -- "
                    f"a negative k_unzip is not a valid rate constant (every unzip "
                    f"consumer is gated on k_unzip > 0, so it would silently "
                    f"masquerade as a frozen pool). Set k_unzip >= 0.")
            if pool.k_unzip > 0.0 and pool.monomer_poly_index is None:
                raise ValueError(
                    f"Pool {pool.label}: k_unzip={pool.k_unzip:g} > 0 but "
                    f"monomer_poly_index is None. The unzip channel drains this "
                    f"pool's condensed moments unconditionally, and without a "
                    f"released-monomer emission target the drained mass leaves "
                    f"the condensed phase silently un-conserved. Wire "
                    f"monomer_poly_index (the released monomer's core index) or "
                    f"set k_unzip=0.")

            # Radical-QSSA unzip channel invariants (M1: the channel is inert
            # -- nothing in the residual reads it until the M2 rate law -- but
            # the solver is the LAST line of defense, same rationale as the
            # k_unzip guards above: a directly-constructed PolymerPoolConfig
            # bypasses the deck helper and PolymerPool.to_config, so the
            # solver re-validates the channel shape itself).
            if pool.radical_qssa_unzip is not None:
                # CONSUME the validator's normalized return and store it back
                # (review round 21, finding 1): a directly-constructed
                # PolymerPoolConfig with a minimal-but-valid channel would
                # otherwise pass validation but stay UN-normalized in storage
                # (missing efficiency/monomer_yield/transfer/basis keys), so
                # any future q["efficiency"] read would KeyError. The config
                # is a frozen dataclass, so the store-back uses
                # object.__setattr__ (the standard frozen-dataclass idiom; no
                # prior precedent in this codebase). deepcopy makes the stored
                # structure defensively independent of every caller-provided
                # sub-dict (finding 2's aliasing half): post-hoc mutation of
                # the caller's dict cannot reach the stored config.
                normalized_qssa = validate_radical_qssa_unzip(
                    pool.label, pool.radical_qssa_unzip)
                object.__setattr__(pool, 'radical_qssa_unzip',
                                   copy.deepcopy(normalized_qssa))
                # (The milestone-i anti-silent-no-op guard that refused
                # weak-link configs here was REMOVED by milestones ii+iii:
                # the U-state slot and the allyl initiation RHS below now
                # consume the vocabulary for real.)
                if pool.monomer_poly_index is None:
                    raise ValueError(
                        f"Pool {pool.label}: radical_qssa_unzip is configured "
                        f"but monomer_poly_index is None. The QSSA unzip "
                        f"channel releases monomer through the pool's existing "
                        f"monomer routing; without an emission target the "
                        f"depropagated repeat units would leave the condensed "
                        f"phase silently un-conserved. Wire monomer_poly_index "
                        f"(the released monomer's core index) or remove "
                        f"radical_qssa_unzip.")
                if pool.k_unzip > 0.0:
                    raise ValueError(
                        f"Pool {pool.label}: radical_qssa_unzip is configured "
                        f"AND k_unzip={pool.k_unzip:g} > 0. These are two "
                        f"representations of the SAME chain-end depropagation "
                        f"channel and are mutually exclusive on a pool "
                        f"(enabling both would double-count the unzip flux). "
                        f"Set k_unzip=0 or remove radical_qssa_unzip.")

            # Radical-homolysis kernel invariants (Stage 1, adjudicated round
            # 66). LAST line of defense, same rationale as above: the deck
            # helper and PolymerPool.to_config guard these shapes, but a
            # directly-constructed PolymerPoolConfig bypasses both. The
            # normalized triplet is stored back (frozen-dataclass store-back
            # idiom, mirroring the QSSA channel above) so the flattener can
            # trust its shape.
            if pool.k_homolysis is not None:
                normalized_khom = validate_k_homolysis(
                    pool.label, pool.k_homolysis)
                object.__setattr__(pool, 'k_homolysis',
                                   copy.deepcopy(normalized_khom))
                if pool.k_scission > 0.0:
                    raise ValueError(
                        f"Pool {pool.label}: k_homolysis is configured AND "
                        f"k_scission={pool.k_scission:g} > 0. Both "
                        f"parameterize the SAME random backbone-break event "
                        f"(homolysis IS random scission, with the products "
                        f"routed to the end-radical pools) and are mutually "
                        f"exclusive on a pool -- enabling both would "
                        f"double-count chain initiation. Set k_scission=0 or "
                        f"remove k_homolysis.")
                if pool.radical_qssa_unzip is not None:
                    raise ValueError(
                        f"Pool {pool.label}: k_homolysis is configured AND "
                        f"radical_qssa_unzip is configured. The QSSA "
                        f"channel's initiation block IS random backbone "
                        f"homolysis, so the two are mutually exclusive on a "
                        f"pool -- enabling both would double-count "
                        f"initiation. Remove one of them.")
                if pool.k_unzip > 0.0:
                    raise ValueError(
                        f"Pool {pool.label}: k_homolysis is configured AND "
                        f"k_unzip={pool.k_unzip:g} > 0. Legacy k_unzip is a "
                        f"phenomenological closed-chain monomer-loss "
                        f"channel, while k_homolysis creates radical-end "
                        f"pools that feed explicit beta-scission/unzip "
                        f"chemistry; the two are mutually exclusive on a "
                        f"pool -- enabling both would double-carry "
                        f"depolymerization. Set k_unzip=0 or remove "
                        f"k_homolysis.")

            # End-radical depropagation kernel invariants (adjudicated round
            # 74 SS2). LAST line of defense, same rationale as the kernels
            # above: the deck helper (parent-declared, daughter-applied) and
            # derive_daughter_pool_configs guard these shapes, but a
            # directly-constructed PolymerPoolConfig bypasses both. The
            # normalized triplet is stored back (frozen-dataclass store-back
            # idiom) so the flattener can trust its shape.
            if pool.k_depropagation is not None:
                normalized_kdep = validate_k_depropagation(
                    pool.label, pool.k_depropagation)
                object.__setattr__(pool, 'k_depropagation',
                                   copy.deepcopy(normalized_kdep))
                if pool.monomer_poly_index is None:
                    raise ValueError(
                        f"Pool {pool.label}: k_depropagation is configured "
                        f"but monomer_poly_index is None. The kernel "
                        f"releases ONE monomer volatile per unzip event; "
                        f"without a resolvable gas emission target the "
                        f"released units would leave the condensed phase "
                        f"silently un-conserved. Wire monomer_poly_index "
                        f"(the released monomer's core index) or remove "
                        f"k_depropagation.")
                if pool.k_unzip > 0.0:
                    raise ValueError(
                        f"Pool {pool.label}: k_depropagation is configured "
                        f"AND k_unzip={pool.k_unzip:g} > 0. Legacy k_unzip "
                        f"is the phenomenological scalar form of the SAME "
                        f"chain-end monomer-release event (and its law "
                        f"stalls at a permanent dmu0 = 0), so the two are "
                        f"mutually exclusive on a pool -- enabling both "
                        f"would double-carry depropagation. Set k_unzip=0 "
                        f"or remove k_depropagation.")
                if pool.radical_qssa_unzip is not None:
                    raise ValueError(
                        f"Pool {pool.label}: k_depropagation is configured "
                        f"AND radical_qssa_unzip is configured. The QSSA "
                        f"channel's depropagation block IS lumped chain-end "
                        f"depropagation, so the two are mutually exclusive "
                        f"on a pool -- enabling both would double-count the "
                        f"unzip flux. Remove one of them.")
                if pool.k_homolysis is not None:
                    raise ValueError(
                        f"Pool {pool.label}: k_depropagation is configured "
                        f"AND k_homolysis is configured. Multi-generation "
                        f"homolysis of radical-ended chains is DEFERRED "
                        f"(adjudicated round 74 SS3: internal homolysis of a "
                        f"radical-ended chain creates radical-state pools "
                        f"this producer cannot generate), and a closed-chain "
                        f"initiation pool has no radical end to "
                        f"depropagate; the two kernels are mutually "
                        f"exclusive on one pool. Remove one of them.")

            # Side-group homolysis kernel invariants (FR1-K1, adjudicated
            # round 70). LAST line of defense, same rationale as the two
            # kernels above: deck helper and PolymerPool.to_config guard
            # these shapes, but a directly-constructed PolymerPoolConfig
            # bypasses both. NOTE the adjudicated coexistence: k_homolysis /
            # k_scission act on backbone C-C bonds while this kernel acts on
            # side-group C-X bonds -- no exclusion between them.
            if (not math.isfinite(pool.chain_mass_defect_g_mol)
                    or pool.chain_mass_defect_g_mol < 0.0):
                raise ValueError(
                    f"Pool {pool.label}: chain_mass_defect_g_mol="
                    f"{pool.chain_mass_defect_g_mol!r} must be a finite "
                    f"value >= 0 [g/mol] (the X-loss feature-pool mass "
                    f"contract; garbage here silently mints/destroys "
                    f"condensed mass).")
            if pool.side_group_homolysis is not None:
                # Shape + selector-vocabulary backstop only: the config
                # carries no monomer structure, so the round-72 STRUCTURAL
                # selector layer (match >= 1, sites_per_unit == match
                # count, disjoint atom sets) runs where the parsed monomer
                # lives -- deck helper, PolymerPool.to_config, producer.
                normalized_sgh = validate_side_group_homolysis(
                    pool.label, pool.side_group_homolysis)
                object.__setattr__(pool, 'side_group_homolysis',
                                   copy.deepcopy(normalized_sgh))
                if pool.k_unzip > 0.0:
                    raise ValueError(
                        f"Pool {pool.label}: side_group_homolysis is "
                        f"configured AND k_unzip={pool.k_unzip:g} > 0. "
                        f"Legacy k_unzip is a phenomenological closed-chain "
                        f"monomer-loss channel, while side_group_homolysis "
                        f"creates radical-defect feature pools that feed "
                        f"explicit degradation chemistry; the two are "
                        f"mutually exclusive on a pool -- enabling both "
                        f"would double-carry degradation. Set k_unzip=0 or "
                        f"remove side_group_homolysis.")
                if pool.radical_qssa_unzip is not None:
                    raise ValueError(
                        f"Pool {pool.label}: side_group_homolysis is "
                        f"configured AND radical_qssa_unzip is configured. "
                        f"The QSSA channel's initiation block is random "
                        f"backbone homolysis feeding the same lumped "
                        f"depropagation, so the two initiation carriers are "
                        f"mutually exclusive on a pool -- enabling both "
                        f"would double-carry initiation. Remove one of "
                        f"them.")
                gi = pool.side_group_gas_indices
                if gi is None or len(gi) != len(normalized_sgh):
                    raise ValueError(
                        f"Pool {pool.label}: side_group_homolysis is "
                        f"configured but the gas-product core indices are "
                        f"not resolved (side_group_gas_indices={gi!r} for "
                        f"{len(normalized_sgh)} channel(s)). The kernel "
                        f"emits one X radical per event; without a resolved "
                        f"gas destination the ejected X would silently "
                        f"vanish (un-conserved mass).")
                for ci, g_idx in enumerate(gi):
                    if not (0 <= g_idx < n_core):
                        raise ValueError(
                            f"Pool {pool.label}: side_group_homolysis "
                            f"channel '{normalized_sgh[ci]['label']}' gas "
                            f"index {g_idx} out of range.")
                    if not self.gas_species_mask[g_idx]:
                        raise ValueError(
                            f"Pool {pool.label}: side_group_homolysis "
                            f"channel '{normalized_sgh[ci]['label']}' "
                            f"gas_product index {g_idx} is NOT masked as "
                            f"GAS -- the ejected X radical is a gas-phase "
                            f"species; emitting it into a condensed slot "
                            f"would corrupt the phase bookkeeping.")

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

        # Flatten the (now validated + normalized) radical_qssa_unzip channels
        # into solver-owned per-pool flat arrays. Runs AFTER the per-pool loop
        # above so every stored channel dict is guaranteed normalized.
        self._flatten_radical_qssa_state()

        # Flatten the (now validated + normalized) k_homolysis kernels and
        # resolve each enabled parent's end-radical daughter pools. Same
        # timing contract as the QSSA flattener above.
        self._flatten_homolysis_state()

        # Flatten the (now validated + normalized) side_group_homolysis
        # channels and resolve each enabled parent's X-loss feature pools +
        # gas destinations (FR1-K1). Same timing contract.
        self._flatten_side_group_state()

        # Flatten the (now validated + normalized) k_depropagation kernels
        # and resolve each enabled pool's released-monomer gas destination
        # (r74 SS2). Same timing contract.
        self._flatten_depropagation_state()

    def _flatten_depropagation_state(self):
        """Flatten each pool's validated+normalized k_depropagation triplet
        into solver-owned per-pool flat arrays (r74 SS2 end-radical
        depropagation kernel).

        Contract (mirrors _flatten_homolysis_state):
        - Populated ONLY here, on the initialize_model ->
          validate_configuration path, from the normalized dict just stored.
        - The RHS reads ONLY these arrays, never the dict.

        Layout (index = pool position in self.polymer_pools):
        - kdep_enabled[i] (int8): 1 iff the pool has k_depropagation.
        - kdep_A/n/Ea[i]: the Arrhenius triplet (A [s^-1], Ea [J/mol]).
        - kdep_gas[i] (int32): core index of the released monomer volatile
          (== pool.monomer_poly_index, guaranteed non-None and gas-masked by
          validate_configuration); -1 on kernel-free pools.
        """
        n_pools = len(self.polymer_pools)
        self.kdep_enabled = np.zeros(n_pools, dtype=np.int8)
        self.kdep_A = np.zeros(n_pools, dtype=float)
        self.kdep_n = np.zeros(n_pools, dtype=float)
        self.kdep_Ea = np.zeros(n_pools, dtype=float)
        self.kdep_gas = np.full(n_pools, -1, dtype=np.int32)
        for i, pool in enumerate(self.polymer_pools):
            trip = pool.k_depropagation
            if trip is None:
                continue
            self.kdep_enabled[i] = 1
            self.kdep_A[i] = trip["A"]
            self.kdep_n[i] = trip["n"]
            self.kdep_Ea[i] = trip["Ea"]
            self.kdep_gas[i] = pool.monomer_poly_index

    def _flatten_homolysis_state(self):
        """Flatten each pool's validated+normalized k_homolysis triplet into
        solver-owned per-pool flat arrays, and resolve the two end-radical
        daughter pools' moment slots (Stage 1, adjudicated round 66).

        Contract (mirrors _flatten_radical_qssa_state):
        - Populated ONLY here, on the initialize_model ->
          validate_configuration path, from the normalized dict just stored.
        - The RHS reads ONLY these arrays, never the dict.

        Layout (index = pool position in self.polymer_pools):
        - khom_enabled[i] (int8): 1 iff the pool has k_homolysis configured.
        - khom_A/n/Ea[i]: the Arrhenius triplet (A [s^-1], Ea [J/mol]).
        - khom_prim_mu0/mu1/mu2[i] and khom_sec_mu0/mu1/mu2[i] (int32): the
          absolute core indices of the {label}_rad_primary_end /
          {label}_rad_secondary_end daughter pools' moment slots; -1 on
          kernel-free pools.

        HARD ERROR when a kernel-enabled pool's daughter pools are missing
        from the configured pool list: the producer creates them at model
        setup, so an unresolvable daughter here means the kernel's fragment
        credits would silently vanish (mass laundered out of the melt).
        """
        n_pools = len(self.polymer_pools)
        # warn-once registry for the runtime moment-shape guard (mirrors
        # self._realizability_warned): one loud line per pool per rebuild.
        self._khom_guard_warned = set()
        self.khom_enabled = np.zeros(n_pools, dtype=np.int8)
        self.khom_A = np.zeros(n_pools, dtype=float)
        self.khom_n = np.zeros(n_pools, dtype=float)
        self.khom_Ea = np.zeros(n_pools, dtype=float)
        self.khom_prim_mu0 = np.full(n_pools, -1, dtype=np.int32)
        self.khom_prim_mu1 = np.full(n_pools, -1, dtype=np.int32)
        self.khom_prim_mu2 = np.full(n_pools, -1, dtype=np.int32)
        self.khom_sec_mu0 = np.full(n_pools, -1, dtype=np.int32)
        self.khom_sec_mu1 = np.full(n_pools, -1, dtype=np.int32)
        self.khom_sec_mu2 = np.full(n_pools, -1, dtype=np.int32)

        label_to_pool = {p.label: p for p in self.polymer_pools}
        for i, pool in enumerate(self.polymer_pools):
            trip = pool.k_homolysis
            if trip is None:
                continue
            missing = [f"{pool.label}{suffix}"
                       for suffix in K_HOMOLYSIS_DAUGHTER_SUFFIXES
                       if f"{pool.label}{suffix}" not in label_to_pool]
            if missing:
                raise ValueError(
                    f"Pool {pool.label}: k_homolysis is enabled but the "
                    f"end-radical daughter pool(s) {missing} are not among "
                    f"the configured solver pools (expected "
                    f"'{pool.label}{K_HOMOLYSIS_DAUGHTER_SUFFIXES[0]}' and "
                    f"'{pool.label}{K_HOMOLYSIS_DAUGHTER_SUFFIXES[1]}'). The "
                    f"producer spawns both at model setup; without them the "
                    f"kernel's fragment credits would silently vanish.")
            self.khom_enabled[i] = 1
            self.khom_A[i] = trip["A"]
            self.khom_n[i] = trip["n"]
            self.khom_Ea[i] = trip["Ea"]
            prim = label_to_pool[
                f"{pool.label}{K_HOMOLYSIS_DAUGHTER_SUFFIXES[0]}"]
            sec = label_to_pool[
                f"{pool.label}{K_HOMOLYSIS_DAUGHTER_SUFFIXES[1]}"]
            (self.khom_prim_mu0[i], self.khom_prim_mu1[i],
             self.khom_prim_mu2[i]) = prim.mu_indices
            (self.khom_sec_mu0[i], self.khom_sec_mu1[i],
             self.khom_sec_mu2[i]) = sec.mu_indices

    def _flatten_side_group_state(self):
        """Flatten each pool's validated+normalized side_group_homolysis
        channel list into solver-owned per-channel row arrays, resolve each
        channel's X-loss feature pool + gas destination, and ENFORCE the
        exact mass-defect contract (FR1-K1, adjudicated round 70).

        Contract (mirrors _flatten_homolysis_state):
        - Populated ONLY here, on the initialize_model ->
          validate_configuration path, from the normalized channel lists.
        - The RHS reads ONLY these arrays, never the pool dicts.

        Layout:
        - sgh_enabled[i] (int8, per pool): 1 iff pool i has channels.
        - sgh_row_start[i]/sgh_row_end[i] (int32, per pool): this pool's
          contiguous slice of the row arrays below (one row per channel).
        - sgh_A/n/Ea/sites[row]: the channel's Arrhenius triplet
          (A [s^-1 per site], Ea [J/mol]) and sites_per_unit.
        - sgh_gas[row] (int32): core index of the ejected X radical.
        - sgh_dst_mu0/mu1/mu2[row] (int32): the X-loss feature pool's
          moment slots.
        - sgh_channel_labels[row]: channel label (guard warn-once keys and
          error messages).

        HARD ERRORS (never silent):
        - a channel whose feature pool is missing from the configured pools
          (the transferred chains would silently vanish);
        - a feature pool whose chain_mass_defect_g_mol does not pin the
          channel's gas_product molar mass M_X (the round-70 #1 P1 trap:
          same-monomer_mw feature pools double-count X unless the defect is
          explicit and exact);
        - a feature pool whose monomer_mw_g_mol diverges from the parent's
          (the chain transfers intact -- a diverging repeat-unit mass
          fabricates/destroys condensed mass).
        """
        n_pools = len(self.polymer_pools)
        # warn-once registry for the runtime out-of-domain guard, keyed
        # (pool label, channel label); one loud line per key per rebuild.
        self._sgh_guard_warned = set()
        self.sgh_enabled = np.zeros(n_pools, dtype=np.int8)
        self.sgh_row_start = np.zeros(n_pools, dtype=np.int32)
        self.sgh_row_end = np.zeros(n_pools, dtype=np.int32)
        rows_A = []
        rows_n = []
        rows_Ea = []
        rows_sites = []
        rows_gas = []
        rows_d0 = []
        rows_d1 = []
        rows_d2 = []
        self.sgh_channel_labels = []

        label_to_pool = {p.label: p for p in self.polymer_pools}
        for i, pool in enumerate(self.polymer_pools):
            self.sgh_row_start[i] = len(rows_A)
            self.sgh_row_end[i] = len(rows_A)
            channels = pool.side_group_homolysis
            if not channels:
                continue
            for ci, ch in enumerate(channels):
                d_label = side_group_daughter_pool_label(
                    pool.label, ch["label"])
                dst = label_to_pool.get(d_label)
                if dst is None:
                    raise ValueError(
                        f"Pool {pool.label}: side_group_homolysis channel "
                        f"'{ch['label']}' is enabled but its X-loss feature "
                        f"pool '{d_label}' is not among the configured "
                        f"solver pools. The producer spawns it at model "
                        f"setup (rmgpy/rmg/input.py); without it the "
                        f"kernel's transferred chains would silently vanish "
                        f"(mass laundered out of the melt).")
                # Exact mass-defect contract (round-70 P1 trap): the feature
                # pool keeps the parent's monomer_mw, so its condensed mass
                # is mu1*MW - mu0*M_X -- the defect MUST pin the channel's
                # gas_product molar mass, and the repeat-unit mass must not
                # diverge from the parent's.
                mw_x = _side_group_gas_mw_g_mol(ch["gas_product"])
                defect = float(dst.chain_mass_defect_g_mol)
                if not (defect > 0.0) or abs(defect - mw_x) > 1.0e-6 * mw_x:
                    raise ValueError(
                        f"Pool {pool.label}: X-loss feature pool "
                        f"'{d_label}' has chain_mass_defect_g_mol="
                        f"{defect!r} but channel '{ch['label']}' ejects "
                        f"gas_product={ch['gas_product']!r} with M_X="
                        f"{mw_x:g} g/mol. The mass contract (feature-pool "
                        f"condensed mass = mu1*MW - mu0*M_X) requires the "
                        f"defect to pin M_X exactly -- anything else mints "
                        f"or destroys condensed mass while gas X appears.")
                mw_parent = float(pool.monomer_mw_g_mol)
                mw_dst = float(dst.monomer_mw_g_mol)
                if abs(mw_dst - mw_parent) > 1.0e-9 * max(abs(mw_parent),
                                                          1.0):
                    raise ValueError(
                        f"Pool {pool.label}: X-loss feature pool "
                        f"'{d_label}' monomer_mw_g_mol={mw_dst:g} diverges "
                        f"from the parent's {mw_parent:g}. The chain "
                        f"transfers INTACT (same repeat unit; the X loss is "
                        f"carried by chain_mass_defect_g_mol, never by the "
                        f"monomer mass), so a diverging monomer_mw "
                        f"fabricates/destroys condensed mass.")
                rows_A.append(ch["A"])
                rows_n.append(ch["n"])
                rows_Ea.append(ch["Ea"])
                rows_sites.append(ch["sites_per_unit"])
                rows_gas.append(int(pool.side_group_gas_indices[ci]))
                d0, d1, d2 = dst.mu_indices
                rows_d0.append(d0)
                rows_d1.append(d1)
                rows_d2.append(d2)
                self.sgh_channel_labels.append(ch["label"])
            self.sgh_enabled[i] = 1
            self.sgh_row_end[i] = len(rows_A)

        self.sgh_A = np.asarray(rows_A, dtype=float)
        self.sgh_n = np.asarray(rows_n, dtype=float)
        self.sgh_Ea = np.asarray(rows_Ea, dtype=float)
        self.sgh_sites = np.asarray(rows_sites, dtype=float)
        self.sgh_gas = np.asarray(rows_gas, dtype=np.int32)
        self.sgh_dst_mu0 = np.asarray(rows_d0, dtype=np.int32)
        self.sgh_dst_mu1 = np.asarray(rows_d1, dtype=np.int32)
        self.sgh_dst_mu2 = np.asarray(rows_d2, dtype=np.int32)

    def _flatten_radical_qssa_state(self):
        """Flatten each pool's validated+normalized radical_qssa_unzip channel
        into solver-owned per-pool flat numpy arrays (review round 21,
        finding 1 -- the M2-prep half; finding 2's mutation-bypass closure).

        Contract:
        - Populated ONLY here, on the initialize_model -> validate_configuration
          path, from the normalized dict validate_configuration just stored.
        - The M2 QSSA rate law will read ONLY these arrays, NEVER the dict.
          Once flattening has run, mutating the (mutable) nested dict behind
          the frozen PolymerPoolConfig is harmless to solver dynamics: the
          arrays are the solver's sole source of truth for the channel.
        - M1: NOTHING in the residual/RHS reads these yet (pinned by the
          bitwise residual-equality test); they are inert state.

        Layout (index = pool position in self.polymer_pools):
        - qssa_enabled[i] (int8): 1 iff the pool has a channel configured.
        - qssa_ki_A/n/Ea[i]: initiation Arrhenius triplet (A [s^-1]).
        - qssa_kdp_A/n/Ea[i]: depropagation triplet (A [s^-1]).
        - qssa_kt_A/n/Ea[i]: termination triplet (A [m^3 mol^-1 s^-1]; the
          only bimolecular block).
        - qssa_efficiency[i], qssa_monomer_yield[i]: f_i, y_m in (0, 1]
          (default 1.0 on channel-absent pools -- the inert value).
        - qssa_has_transfer[i] (int8) + qssa_ktr_A/n/Ea[i]: optional transfer
          triplet (A [s^-1], pseudo-first-order: ktr multiplies R directly
          in the M2 rate law); zeros when absent.
        Disabled rows keep zero Arrhenius slots: every M2 consumer must gate
        on qssa_enabled, never on A != 0.

        Weak-link U-state extension (milestones ii+iii):
        - qssa_weaklink[i] (int8): 1 iff the channel carries the weak-link
          vocabulary. On weak-link rows the legacy summed qssa_kt_* slots
          STAY ZERO (there is no summed 'termination' block); the RHS must
          gate on qssa_weaklink, never on kt_A != 0.
        - qssa_kia_A/n/Ea[i]: initiation_allyl triplet (A [s^-1]).
        - qssa_ktrec_A/n/Ea[i] / qssa_ktdisp_A/n/Ea[i]: split termination
          triplets (A [m^3 mol^-1 s^-1]); the RHS uses kt_total =
          kt_rec(T) + kt_disp(T) wherever the legacy summed kt appeared,
          and ONLY kt_disp sources dU/dt.
        - qssa_u0[i]: unsaturated_tail_ends_initial [mol] (state IC, set
          into y0 by set_initial_conditions after the census trap there).
        - qssa_u_slot[i] (int32): absolute ODE index of the pool's U slot
          (n_core + running weak-link ordinal, matching the count
          initiate_tolerances used to extend neq), or -1.
        """
        n_pools = len(self.polymer_pools)
        self.qssa_enabled = np.zeros(n_pools, dtype=np.int8)
        self.qssa_ki_A = np.zeros(n_pools, dtype=float)
        self.qssa_ki_n = np.zeros(n_pools, dtype=float)
        self.qssa_ki_Ea = np.zeros(n_pools, dtype=float)
        self.qssa_kdp_A = np.zeros(n_pools, dtype=float)
        self.qssa_kdp_n = np.zeros(n_pools, dtype=float)
        self.qssa_kdp_Ea = np.zeros(n_pools, dtype=float)
        self.qssa_kt_A = np.zeros(n_pools, dtype=float)
        self.qssa_kt_n = np.zeros(n_pools, dtype=float)
        self.qssa_kt_Ea = np.zeros(n_pools, dtype=float)
        self.qssa_efficiency = np.ones(n_pools, dtype=float)
        self.qssa_monomer_yield = np.ones(n_pools, dtype=float)
        self.qssa_has_transfer = np.zeros(n_pools, dtype=np.int8)
        self.qssa_ktr_A = np.zeros(n_pools, dtype=float)
        self.qssa_ktr_n = np.zeros(n_pools, dtype=float)
        self.qssa_ktr_Ea = np.zeros(n_pools, dtype=float)
        self.qssa_weaklink = np.zeros(n_pools, dtype=np.int8)
        self.qssa_kia_A = np.zeros(n_pools, dtype=float)
        self.qssa_kia_n = np.zeros(n_pools, dtype=float)
        self.qssa_kia_Ea = np.zeros(n_pools, dtype=float)
        self.qssa_ktrec_A = np.zeros(n_pools, dtype=float)
        self.qssa_ktrec_n = np.zeros(n_pools, dtype=float)
        self.qssa_ktrec_Ea = np.zeros(n_pools, dtype=float)
        self.qssa_ktdisp_A = np.zeros(n_pools, dtype=float)
        self.qssa_ktdisp_n = np.zeros(n_pools, dtype=float)
        self.qssa_ktdisp_Ea = np.zeros(n_pools, dtype=float)
        self.qssa_u0 = np.zeros(n_pools, dtype=float)
        self.qssa_u_slot = np.full(n_pools, -1, dtype=np.int32)

        u_count = 0
        for i, pool in enumerate(self.polymer_pools):
            q = pool.radical_qssa_unzip
            if q is None:
                continue
            self.qssa_enabled[i] = 1
            self.qssa_ki_A[i] = q["initiation"]["A"]
            self.qssa_ki_n[i] = q["initiation"]["n"]
            self.qssa_ki_Ea[i] = q["initiation"]["Ea"]
            self.qssa_kdp_A[i] = q["depropagation"]["A"]
            self.qssa_kdp_n[i] = q["depropagation"]["n"]
            self.qssa_kdp_Ea[i] = q["depropagation"]["Ea"]
            if "initiation_allyl" in q:
                # Weak-link channel: split termination blocks (validator
                # invariant: no legacy summed 'termination' key exists on
                # this channel, so qssa_kt_* stays 0 on this row).
                self.qssa_weaklink[i] = 1
                self.qssa_kia_A[i] = q["initiation_allyl"]["A"]
                self.qssa_kia_n[i] = q["initiation_allyl"]["n"]
                self.qssa_kia_Ea[i] = q["initiation_allyl"]["Ea"]
                self.qssa_ktrec_A[i] = q["termination_recombination"]["A"]
                self.qssa_ktrec_n[i] = q["termination_recombination"]["n"]
                self.qssa_ktrec_Ea[i] = q["termination_recombination"]["Ea"]
                self.qssa_ktdisp_A[i] = q["termination_disproportionation"]["A"]
                self.qssa_ktdisp_n[i] = q["termination_disproportionation"]["n"]
                self.qssa_ktdisp_Ea[i] = q["termination_disproportionation"]["Ea"]
                self.qssa_u0[i] = q["unsaturated_tail_ends_initial"]
                self.qssa_u_slot[i] = self.num_core_species + u_count
                u_count += 1
            else:
                self.qssa_kt_A[i] = q["termination"]["A"]
                self.qssa_kt_n[i] = q["termination"]["n"]
                self.qssa_kt_Ea[i] = q["termination"]["Ea"]
            self.qssa_efficiency[i] = q["efficiency"]
            self.qssa_monomer_yield[i] = q["monomer_yield"]
            if q["transfer"] is not None:
                self.qssa_has_transfer[i] = 1
                self.qssa_ktr_A[i] = q["transfer"]["A"]
                self.qssa_ktr_n[i] = q["transfer"]["n"]
                self.qssa_ktr_Ea[i] = q["transfer"]["Ea"]

        # Layout cross-check: initiate_tolerances extended neq by ITS count
        # of weak-link pools (key presence on the raw configs, before
        # validation); the flattening above counted the normalized ones. A
        # mismatch means the pool configs changed between the two passes --
        # the U slots would alias core species. Fail loud, never integrate
        # a misaligned state vector.
        if u_count != self.num_qssa_u:
            raise ValueError(
                f"Weak-link U-state layout drift: initiate_tolerances "
                f"allocated {self.num_qssa_u} U slot(s) but "
                f"validate_configuration found {u_count} weak-link pool(s). "
                f"The pool configs changed between ODE layout and "
                f"validation; rebuild the reactor (initialize_model) instead "
                f"of mutating pool configs in place.")

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
        gas_species_mask OR is_polymer_proxy-tagged, unvetoed AND dual-axis
        polymer-sized -- spec §5.1 C3-amended, r89 dual-axis), logs a census
        above
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

        # Chain-scale window (spec §5.3): largest configured pool monomer MW
        # + slack, in kg/mol. Since the r89 dual-axis amendment this window
        # serves ONLY the provenance counterparty predicate (i) below -- the
        # physically-melt class gate (tag branch) uses the dual-axis
        # polymer-sized test instead (the window form melt-classified DP-2
        # gas volatiles: PP run-9 1,5-hexadiene, 82.15 g/mol > propene
        # window 52.1 but only 2.0 monomer-equivalents). Predicate (ii)
        # (sharing the proxy's saturated parent) was REJECTED as the
        # cheap-at-init test: it needs graph saturation + isomorphism per
        # gas species per rebuild; (i) is a float compare on data already in
        # hand. FR1-K2 mass audit note (round-72 P2): deliberately NOT
        # defect-aware. This is a PER-REPEAT-UNIT MW scale window, not a
        # condensed-mass computation; chain_mass_defect_g_mol is a
        # PER-CHAIN correction (one lost X per chain), and X-loss feature
        # pools keep the parent's monomer_mw by design (intact-chain
        # transfer), so the per-unit window is exactly right as is.
        chain_window_kg = (max((float(p.monomer_mw_g_mol) for p in self.polymer_pools),
                               default=0.0)
                           + REFERENCE_STATE_MW_SLACK_G_MOL) / 1000.0

        # r89 dual-axis gate inputs. The threshold is imported at call time
        # (documented solver <-> polymer module cycle, same posture as
        # rmgpy.polymer's function-local solver imports) so the gate and
        # rmgpy.polymer._discrete_is_polymer_sized share ONE constant; the
        # per-pool axes come from the flattened configs (monomer_mw_g_mol +
        # monomer_heavy_atoms), plumbed producer-side from the monomer
        # STRUCTURE -- the heavy axis is never approximated from mass.
        from rmgpy.polymer import _IMPOSTOR_DISCRETE_MONOMER_UNITS as _units
        pool_axes = [(float(getattr(p, "monomer_mw_g_mol", 0.0) or 0.0),
                      int(getattr(p, "monomer_heavy_atoms", 0) or 0))
                     for p in self.polymer_pools]

        is_melt = [False] * n_core
        mws = [0.0] * n_core
        heavies = [-1] * n_core   # -1 = unknown (no structure, no stamp)
        # Backstop census: polymer-sized (dual-axis) members whose melt
        # classification the durable gas veto SUPPRESSED. create_reacted_copy
        # returns None both for a genuine gas volatile and for a wing-match
        # FAILURE on a real chain-scale fragment; the veto trusts that None as
        # authoritative "gas", so a mis-vetoed genuine chain would be silently
        # dropped from the melt sum instead of loudly refused. Recording +
        # warning here keeps that case visible for a human without regressing
        # the alpha-methylstyrene fix (the build no longer refuses). r89: the
        # census keys on the SAME dual-axis chain-scale notion as the gate --
        # a genuinely-volatile vetoed species (AMS, 1.13 monomer-equivalents)
        # is no longer announced; only chain-scale suppressions are.
        veto_suppressed_chain_scale = []
        # r89 undecidable-axis census (mirror of rmgpy.polymer
        # ._warn_impostor_axis_undecidable): a tag-branch candidate whose
        # polymer-sized verdict could not be established on BOTH axes is
        # NEVER classified blind -- it degrades to conservative-gas -- but
        # the degradation is announced, per species, so a silent axis-data
        # gap cannot quietly unpair a genuine chain.
        axis_undecidable = []
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
            # Heavy-atom (non-H) count for the r89 dual-axis gate: from the
            # structure when one exists; otherwise from the species-level
            # props stamp the consumer runner populates from the artifact's
            # composition (real composition data both ways -- the r89
            # data-flow constraint forbids deriving this axis from mass).
            # -1 = unknown: the heavy axis is then uncomputable and the
            # dual-axis verdict degrades to conservative-gas + census.
            _props = getattr(spc, "props", None)
            if mol_list:
                _m0 = mol_list[0]
                heavies[i] = _m0.get_num_atoms() - _m0.get_num_atoms('H')
            elif isinstance(_props, dict):
                # HARDCODED literal must equal rmgpy.polymer
                # .POLYMER_HEAVY_ATOM_COUNT_KEY -- pinned by
                # test_heavy_atom_key_literal_matches_solver_gate
                # (rename both).
                _hv = _props.get("polymer_heavy_atom_count")
                if _hv is not None:
                    heavies[i] = int(_hv)
            # Physically-melt CLASS (spec §5.1 C3-AMENDED, r89 dual-axis):
            # the condensed branch (pool-configured by input)
            # unconditionally; the tag branch only when DUAL-AXIS
            # POLYMER-SIZED (mass AND heavy axes both computable and both
            # >= _IMPOSTOR_DISCRETE_MONOMER_UNITS monomer-equivalents
            # against >= 1 pool -- _dual_axis_polymer_sized, the exact
            # mirror of rmgpy.polymer._discrete_is_polymer_sized). The size
            # conjunct is part of the class DEFINITION, not bolted onto the
            # provenance set: family.py:1657 blanket-tags every structure
            # of a proxy-touching reaction (including H2), so a raw
            # tag-read would be correct only by spawn-pass ordering -- the
            # structural gate cannot be broken by a lifecycle reordering
            # the way a tag-read can. A stale tag on a below-threshold
            # species simply FAILS the conjunct and is excluded: expected
            # and silent (its gas reference state is CORRECT -- PP run-9's
            # 1,5-hexadiene, a genuine DP-2 gas volatile at 2.0
            # monomer-equivalents, is exactly this class; the pre-r89 MW
            # window wrongly melt-classified it and refused an all-gas
            # allyl recombination at U ~ 11 decades).
            #
            # DURABLE GAS-VOLATILE VETO (rmgpy.polymer): the tag branch keys off
            # is_polymer_proxy, which is a monotonic multi-writer sticky cache
            # with no authoritative "gas" clear point -- a genuine discrete gas
            # volatile (e.g. alpha-methylstyrene) that got proxy-contaminated
            # and happens to clear the MW window would be wrongly counted as a
            # melt participant (the reference-state-tripwire leak). The polymer
            # handshake / chip demotion stamps such a species with a POSITIVE,
            # durable veto in props (POLYMER_REFERENCE_STATE_GAS_VETO_KEY) that
            # the proxy stamping machinery never touches; honor it here. C3
            # mask-lagged chains never receive the veto, so they stay melt.
            _gas_veto = False
            if isinstance(_props, dict):
                # HARDCODED literal must equal rmgpy.polymer
                # .POLYMER_REFERENCE_STATE_GAS_VETO_KEY -- pinned by
                # test_gas_veto_key_literal_matches_solver_gate (rename both).
                _gas_veto = bool(_props.get("polymer_reference_state_gas_veto", False))
            _proxy_i = bool(getattr(spc, "is_polymer_proxy", False))
            # r89 gate order: veto precedence is unchanged (a vetoed species
            # NEVER classifies melt via the tag branch, whatever its size);
            # the dual-axis verdict is evaluated for every gas-masked
            # tagged candidate because both censuses key on it.
            _sized = False
            _missing = None
            if mask[i] and _proxy_i:
                _sized, _missing = _dual_axis_polymer_sized(
                    mws[i] * 1000.0, heavies[i], pool_axes, _units)
                if not _sized and _missing is not None:
                    axis_undecidable.append(
                        (getattr(spc, "label", "?"), _missing))
            is_melt[i] = ((not mask[i])
                          or (_proxy_i and not _gas_veto and _sized))
            # The veto is only a backstop concern when it suppressed an
            # otherwise-melt POLYMER-SIZED tag-branch member (gas-masked,
            # proxy, dual-axis sized). A below-threshold vetoed species
            # would fail the size conjunct anyway -- its exclusion is
            # correct and silent (r89: alpha-methylstyrene, 1.13
            # monomer-equivalents, is no longer announced).
            if _gas_veto and _proxy_i and mask[i] and _sized:
                veto_suppressed_chain_scale.append(
                    (getattr(spc, "label", "?"), mws[i]))

        # Preserve the reference-state backstop for the ambiguous handshake
        # false-None case (code-review IMPORTANT #1): a genuine chain-scale
        # fragment the wing-matcher missed is now gas-vetoed rather than
        # refused, so surface it for a human. NOT a refusal -- census only.
        self.gas_veto_census = list(veto_suppressed_chain_scale)
        if veto_suppressed_chain_scale:
            logging.warning(
                "THERMO REFERENCE-STATE GAS VETO: %d chain-scale "
                "(dual-axis polymer-sized at >= %.1f monomer-equivalents "
                "of mass AND heavy atoms) product(s) were excluded from the "
                "melt reference-state by the durable gas-volatile veto: %s. "
                "If any is actually a polymer chain, it is a handshake "
                "create_reacted_copy false-None (wing-match failure), NOT a "
                "thermo problem -- investigate the handshake, do not touch "
                "reference states.",
                len(veto_suppressed_chain_scale), _units,
                ", ".join("%s (%.1f g/mol)" % (lbl, mw * 1000.0)
                          for lbl, mw in veto_suppressed_chain_scale))

        # r89 undecidable-axis census: conservative-gas degradations are
        # announced, never silent (mirror of rmgpy.polymer
        # ._warn_impostor_axis_undecidable -- never classify blind; the
        # solver tripwire's condensed branch and the r71 hard-fail stay the
        # loud backstops).
        self.reference_state_axis_undecidable = list(axis_undecidable)
        if axis_undecidable:
            logging.warning(
                "THERMO REFERENCE-STATE AXIS UNDECIDABLE (r89 census): %d "
                "proxy-tagged gas-masked candidate(s) had an uncomputable "
                "polymer-sized axis, so the dual-axis melt conjunct could "
                "not be established on BOTH axes; classifying "
                "conservative-gas (never melt-classify blind): %s. If one "
                "is a genuine chain, plumb its missing axis data (species "
                "composition / pool monomer structure), do not touch "
                "reference states.",
                len(axis_undecidable),
                ", ".join("%s (missing: %s)" % (lbl, miss)
                          for lbl, miss in axis_undecidable))

        offenders = []
        for rxn in core_reactions:
            if not getattr(rxn, "reversible", False):
                continue
            if getattr(rxn, "polymer_refused", False):
                # Refused rows (item-18 stamp-but-keep + the PP v1
                # gas-association refusal, adjudicated round 63) have their
                # WHOLE flux suppressed in the residual (reaction_refused),
                # so their reference-state pairing is moot: no Keq of theirs
                # ever moves a mole. Skip them here -- keyed on the refused
                # stamp ONLY, no other relaxation (the SAME row unrefused
                # still trips; pinned by
                # test_refused_row_with_unpaired_reference_state_passes_initialize).
                continue
            j = self.reaction_index[rxn]
            r_idx = [int(k) for k in self.reactant_indices[j] if 0 <= k < n_core]
            p_idx = [int(k) for k in self.product_indices[j] if 0 <= k < n_core]
            melt_r = [k for k in r_idx if is_melt[k]]
            melt_p = [k for k in p_idx if is_melt[k]]
            if not melt_r and not melt_p:
                continue  # all-gas reaction: gas reference uniformly correct

            # Leak self-assertion (spec §5.1 C3 amendment, r89 dual-axis):
            # cannot-happen guard on every species entering the melt sum. A
            # gas-classified member can only be here via the tag branch,
            # whose dual-axis polymer-sized conjunct is enforced in the gate
            # above; if a below-threshold (or axis-undecidable) species ever
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
                    core_species[k].label, mws[k], heavies[k], bool(mask[k]),
                    pool_axes, _units)

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
        # ONE base-label convention (function-local import: documented
        # solver<->polymer cycle note at the module head): trailing-index
        # strip ONLY -- the former first-'(' truncation could alias a
        # SMILES-labelled species (structural parentheses) onto a pool label.
        from rmgpy.polymer import strip_rmg_index_suffix
        n_core = self.num_core_species
        n_total = len(species_list)
        for pool_i, pool in enumerate(self.polymer_pools):
            for i in range(n_total):
                base_label = strip_rmg_index_suffix(species_list[i].label)
                if base_label == pool.label:
                    mask_arr[i] = False
                    if record_indices:
                        self.species_to_pool_indices[i] = pool_i
                        self.is_pool_proxy[i] = 1
                    break
            mu1_target_label = f"{pool.label}_mu1"
            for i in range(n_total):
                # Handle RMG renaming: "PS_mu1(2)" -> "PS_mu1"
                base_label = strip_rmg_index_suffix(species_list[i].label)
                if base_label == mu1_target_label:
                    if record_indices:
                        self.pool_mu1_indices[pool_i] = i
                    mask_arr[i] = False  # Ensure it's poly phase

            mu0_target_label = f"{pool.label}_mu0"
            for i in range(n_total):
                base_label = strip_rmg_index_suffix(species_list[i].label)
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

            # monomer_poly_index (the unzip/QSSA release target) is
            # deliberately NOT mapped or condensed here (incident 2026-07-03,
            # design B-prime): the release target is the mechanism's gas
            # volatile (e.g. styrene for PS), and force-condensing it made
            # every reversible gas core reaction producing it carry an
            # ~11-decade unpaired reference-state term (correctly refused by
            # the tripwire). The release deposits into the GAS amount basis
            # (small_src -> dn_dt, mol/s); validate_configuration now REQUIRES
            # the monomer index to be masked GAS.

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

        # r71 FIX 2 (PP run-5 stall) -- rebuild RESTAMPING at the last honest
        # chokepoint. The generation-time refusal stamps are ad-hoc object
        # attributes that have PROVEN losable between stamping and this
        # rebuild (canonical dedup in check_for_existing_reaction; the
        # __reduce__ serialization contract -- rmgpy/reaction.py documents
        # "the solver re-derives them at initialize_model", and THIS is that
        # re-derivation). Idempotently re-run the shape predicate over every
        # row (core AND edge) BEFORE reaction_refused is captured and BEFORE
        # the thermo reference-state tripwire scans core rows: this loop sees
        # the final post-promotion/post-dedup list regardless of construction
        # path. stamp_gas_association_refusal early-returns on already-refused
        # rows, so an upstream qssa-invalid census reason is never overwritten
        # (precedence preserved) and stamped rows are unchanged. The other
        # refusal sibling (stamp_polymer_flux_archetype's UNRESOLVED
        # FEATURE-radical branch) is NOT re-runnable here -- it needs the
        # generation-time handshake context (resolved polymer_reactants, chip
        # surgery) -- which is exactly why FIX 4 below hard-fails any
        # remaining unclassified live proxy row instead of trusting it.
        # Function-local import: avoids a solver->polymer module cycle.
        # pool_registry (r87 shape B): the reference-state-split
        # isomerization row has no Polymer participant, so the restamp
        # supplies the registered pools collected from the SPECIES lists --
        # real Polymer objects with monomer structures (the r85 dual-axis
        # size gate stays computable), NOT the structureless
        # PolymerPoolConfig sidecars. Collected once per rebuild.
        from rmgpy.polymer import collect_polymer_pool_registry, stamp_gas_association_refusal
        restamp_pool_registry = collect_polymer_pool_registry(
            core_species or [], edge_species or [])
        for rxn in itertools.chain(core_reactions, edge_reactions):
            stamp_gas_association_refusal(rxn, pool_registry=restamp_pool_registry)

        # item 18: capture the Task-3 refusal stamp in the SAME chain(core, edge)
        # order so self.reaction_refused[r_idx] aligns with the residual's r_idx
        # exactly like is_end_group_reaction. A refused FEATURE-radical reaction's
        # whole flux is suppressed in the residual so it fabricates no MW-weighted
        # backbone mass (stamp-but-keep: it stays in the edge/model).
        self.reaction_refused = np.zeros(n_rxn, dtype=np.int8)
        for i, rxn in enumerate(itertools.chain(core_reactions, edge_reactions)):
            if getattr(rxn, "polymer_refused", False):
                self.reaction_refused[i] = 1

        # item 18 (T3): correct-but-loud census of refused FEATURE-radical
        # reactions, tagged by radical class (eliminating vs accumulating) and
        # refuse reason (conduit-deferred vs qssa-invalid). Built once per
        # rebuild over the SAME chain(core, edge) order as the capture above so
        # reaction_refused[i] aligns. Changes NO flux/state -- it only reports
        # (mirrors the thermo reference_state_census posture). Helpers are
        # function-local imports to avoid a solver->polymer module cycle.
        from rmgpy.polymer import _warn_once_refused, _reaction_census_label, _warn_once_double_count
        self.refused_reaction_census = []
        for i, rxn in enumerate(itertools.chain(core_reactions, edge_reactions)):
            if not self.reaction_refused[i]:
                continue
            accumulating = bool(getattr(rxn, "polymer_refused_accumulating", False))
            entry = {
                "reaction": _reaction_census_label(rxn),
                "radical_class": "accumulating" if accumulating else "eliminating",
                "reason": "qssa-invalid" if accumulating else "conduit-deferred",
            }
            self.refused_reaction_census.append(entry)
            _warn_once_refused(entry)

        # Per-reaction pool moment-flux archetype (spec 2026-06-09). Same
        # chain(core, edge) order as is_end_group_reaction so indices match
        # r_idx in the residual.
        self.reaction_flux_archetype = np.zeros(n_rxn, dtype=np.int8)
        self.reaction_src_pool = np.full(n_rxn, -1, dtype=np.int32)
        self.reaction_dst_pool = np.full(n_rxn, -1, dtype=np.int32)
        # Stamped chip repeat-unit counts (spec 2026-06-10 §4.3); same
        # chain(core, edge) order so the index matches r_idx in the residual.
        self.reaction_chip_units = np.zeros(n_rxn, dtype=np.int32)
        # Stamped volatile-ejection unit debit a = Sigma(volatile MW) /
        # source-pool monomer MW (fractional, e.g. ~1.135 for
        # alpha-methylstyrene off a styrene pool). Same chain(core, edge)
        # order so the index matches r_idx in the residual.
        self.reaction_eject_units = np.zeros(n_rxn, dtype=np.float64)
        for i, rxn in enumerate(itertools.chain(core_reactions, edge_reactions)):
            self.reaction_flux_archetype[i] = int(getattr(rxn, "polymer_flux_archetype", 0))
            self.reaction_chip_units[i] = int(getattr(rxn, "polymer_chip_units", 0))
            self.reaction_eject_units[i] = float(getattr(rxn, "polymer_eject_units", 0.0))
        # Item 17 (spec 2026-06-12 SS3(e)): generation-time stamps, captured
        # AFTER the stamp-read loop and BEFORE the init demotion pass
        # (:NONE->UNRESOLVED + unresolvable stamped shapes) mutates the
        # array in place. For reactions kinetics-FLIPPED at GENERATION time
        # (restamp_flipped_polymer_archetype, r92: restamp-or-refuse mutates
        # the object itself), "pre-demotion" deliberately means
        # pre-SOLVER-demotion: the captured value is the restamped
        # flipped-direction archetype (or UNRESOLVED for flip-REFUSED rows,
        # which are flux-dead via reaction_refused) and the census reports
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

        # Characteristic-flux include-mask (fix #2a, 2026-06-27). The base
        # enlargement char_rate is an L2 norm over core_species_rates; the
        # moment-coordinate positions (PS_mu0/_mu1/_mu2) carry moment
        # derivatives (the residual writes dn_dt[mu_indices] = dmu*_dt), not
        # molar species fluxes, and contaminate that norm: the lumped k_unzip
        # channel scales char_rate ~linearly in k_unzip and buries family
        # chemistry below tol_move_to_core, leaving the core empty.
        #
        # Authority = the pool mu_indices, NOT is_moment_dummy alone. Species
        # carry an is_moment_dummy flag (set in model.py when the dummies are
        # created), but Species.copy() (species.py:784) does NOT preserve it,
        # so a copied core dummy reads False and an is_moment_dummy-only mask
        # would be silently live-inert. mu_indices are integer core positions
        # the solver binds itself from the pool configs -- copy-proof and
        # exactly the positions the residual treats as moment coordinates. We
        # union with the flag too (cheap; covers any flagged dummy outside an
        # active pool's mu_indices). Everything real -- the proxy, the solvent,
        # the monomer-release routing species -- stays True: we drop bookkeeping
        # coordinates, NOT real chemistry (2a, not 2b). Rebuilt every
        # initialize_model (per RMG iteration) so it tracks the live core.
        include_mask = np.array(
            [not getattr(core_species[i], "is_moment_dummy", False)
             for i in range(n_core)], dtype=bool)
        for pool in self.polymer_pools:
            for idx in pool.mu_indices:
                if 0 <= idx < n_core:
                    include_mask[idx] = False
        self._char_rate_include_mask = include_mask

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
        # A5-2 (spec 2026-06-12 SS3(a) REDESIGNED): stage-1 precedence, in
        # order. The edge-provenance marker records HOW the edge suffix was
        # produced (1 = stage-1-classified, 0 = default-filled) so R1-EDGE can
        # raise on a production fallback R1's core-prefix check is blind to.
        edge_provenance = np.ones(n_edge_spc, dtype=np.int8)
        if self._prospective_classifier is not None:
            # (1) PRODUCTION PATH: rebuild stage 1 against the LIVE
            # chain(core, edge) via the plumbed classifier handle
            # (polymerPhase.get_gas_mask over the CURRENT lists -- no frozen
            # seed to go stale). The edge suffix is genuinely classified, so
            # provenance stays all-1.
            stage1 = np.asarray(
                self._prospective_classifier(combined_species), dtype=bool)
            if stage1.shape[0] != n_core + n_edge_spc:
                raise ValueError(
                    "PROSPECTIVE-MASK CLASSIFIER: prospective_classifier "
                    "returned length %d for a chain(core, edge) of length %d."
                    % (stage1.shape[0], n_core + n_edge_spc))
            self.prospective_gas_mask = stage1.copy()
        elif seed is not None and len(seed) == n_core + n_edge_spc:
            # (2) DIRECT-TEST SEED PATH: a full doctored seed supplied at
            # construction (e.g. R1's T4). The seed IS a stage-1 product by
            # contract, so provenance is stage-1-classified.
            self.prospective_gas_mask = np.asarray(seed, dtype=bool).copy()
        else:
            # (3) FALLBACK -- genuine last resort, direct-test/runner only.
            # Edge defaults prospectively-GAS; the edge suffix provenance is
            # default-filled (0), so R1-EDGE raises unless the build is flagged
            # allow_default_prospective_edge. The PROSPECTIVE-MASK SEED STALE
            # warning is retained ONLY for back-compat when a stale seed was
            # present (a production build never reaches here -- classifier
            # present => branch 1).
            if seed is not None:
                logging.warning(
                    "PROSPECTIVE-MASK SEED STALE: stage-1 seed length %d != "
                    "n_core + n_edge = %d; rebuilding prospective mask from "
                    "the fallback (core mask + edge defaults GAS).",
                    len(seed), n_core + n_edge_spc)
            # Stated premise, not an implication (amendment A3): edge species
            # default prospectively-GAS here. Fixtures express a
            # prospectively-condensed EDGE species through configured pool
            # labels (stage 2) -- exactly how production daughters become
            # condensed under item 16.
            self.prospective_gas_mask = np.concatenate([
                np.asarray(self.gas_species_mask, dtype=bool).copy(),
                np.ones(n_edge_spc, dtype=bool)])
            edge_provenance = np.zeros(n_edge_spc, dtype=np.int8)
        self._prospective_edge_provenance = edge_provenance
        self._apply_pool_phase_overrides(self.prospective_gas_mask,
                                         combined_species,
                                         record_indices=False)

        # Edge-daughter prospective condensed-mask (spec 2026-06-29). After the
        # stage-2 override, BEFORE R1: flip qualifying EDGE daughters CONDENSED
        # so Gate B (is_poly_event AND not has_condensed_prod) does not zero
        # their real scission flux. The classifier is re-run over the LIVE
        # combined list here (callable, never a frozen set). EDGE slots only --
        # R1 (next) re-verifies the core prefix is untouched. LOUD on failure:
        # a raising classifier must NOT silently default the daughter to GAS
        # (that would re-hide Gate B).
        if self.prospective_condensed_edge_daughter_classifier is not None:
            try:
                condensed_bases = \
                    self.prospective_condensed_edge_daughter_classifier(combined_species)
            except Exception:
                logging.error(
                    "EDGE-DAUGHTER CONDENSED-MASK: classifier raised while "
                    "building prospective_gas_mask (n_core=%d, n_edge=%d); not "
                    "defaulting to GAS.", n_core, n_edge_spc)
                raise
            # LOCKSTEP with the classifier's base-label convention
            # (rmgpy.polymer.strip_rmg_index_suffix, via polymer_input's
            # _base_label): trailing-index strip ONLY. The former first-'('
            # truncation aliased SMILES-labelled species onto shared bases
            # (PP run-5 tripwire crash).
            from rmgpy.polymer import strip_rmg_index_suffix
            for i in range(n_core, n_core + n_edge_spc):
                if strip_rmg_index_suffix(combined_species[i].label) in condensed_bases:
                    self.prospective_gas_mask[i] = False

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
                + ". Known suspects: (a) duplicate-label fallback "
                "disablement in the combined get_gas_mask call (a duplicate "
                "label present only in the edge disables label_fallback_safe "
                "for the combined list -- PolymerPhase.get_gas_mask, "
                "polymer_input.py); (b) any base-label producer/consumer "
                "drifting from rmgpy.polymer.strip_rmg_index_suffix (the PP "
                "run-5 crash was first-'(' truncation aliasing SMILES "
                "labels onto shared bases -- fixed 2026-07-05).")

        # RIDER R1-EDGE -- edge-suffix provenance guard (item 17 A5-2, spec
        # SS3(d); DISTINCT from R1 by construction). R1's core-prefix check
        # CANNOT see the §3(a) staging defect: when a build falls through to
        # edge-defaults-GAS, the core PREFIX still matches gas_species_mask
        # verbatim (the fallback copies it), so R1 passes green while the edge
        # SUFFIX is silently default-filled. A PRODUCTION build (not flagged
        # allow_default_prospective_edge) whose edge suffix came from defaults
        # means the live-edge stage-1 rebuild never fired -- the latent no-op
        # A5-2 exists to kill. Raise, never warn. The legitimate last-resort
        # fallback (direct-test/runner, no blueprint phase) is opt-in via the
        # flag. Runs at the same point R1 does, every rebuild.
        if (not self._allow_default_prospective_edge
                and self._prospective_edge_provenance is not None
                and self._prospective_edge_provenance.shape[0] > 0
                and not np.all(self._prospective_edge_provenance)):
            n_default = int(np.count_nonzero(
                self._prospective_edge_provenance == 0))
            raise ValueError(
                "PROSPECTIVE-MASK PROVENANCE: fallback fired during a "
                "production enlargement: %d edge entries, %d default-filled; "
                "a production build must classify the live edge via stage-1 "
                "(prospective_classifier over chain(core, edge)) -- the "
                "rebuild trigger missed an edge-list change, or no classifier "
                "handle was plumbed onto the engine. Direct-test/runner "
                "construction that legitimately needs the fallback must set "
                "allow_default_prospective_edge=True."
                % (self._prospective_edge_provenance.shape[0], n_default))

        # Resolve per-reaction source/target pools from species indices (no
        # label matching). r71 FIX 4 (PP run-5 stall): an unstamped
        # proxy-touching row (archetype NONE after the FIX 2 restamp) that is
        # NOT refused has no classification the solver can conserve mass
        # under -- the former silent NONE -> UNRESOLVED remap (legacy
        # mu1-only flux) is exactly what ran live against run-5's exhausted
        # pool. Such rows now HARD-FAIL unless the build explicitly opts into
        # the legacy fallback (allow_unstamped_proxy_rows=True --
        # direct-test/runner construction only, mirroring
        # allow_default_prospective_edge). Refused rows stay archetype NONE:
        # they are flux-dead via reaction_refused, and routing them onto the
        # legacy path would misreport an adjudicated row as legacy-live.
        ir_arr = self.reactant_indices
        ip_arr = self.product_indices
        n_unstamped = 0
        n_unresolvable = 0
        unstamped_live_rows = []
        combined_rxns = list(itertools.chain(core_reactions, edge_reactions))
        # r91 (PP run-10, spar-adjudicated): spawned-pool demotion refusal
        # bookkeeping. A stamped pool-coupled row whose required pool
        # endpoint resolves to -1 BECAUSE the endpoint is a mid-run-SPAWNED,
        # config-less pool (a Polymer species present in the species lists
        # but absent from self.polymer_pools) is REFUSED conduit-deferred
        # (reaction_refused -> zero whole-row flux, the same mechanics as
        # every other refused row) instead of demoted to the legacy mu1-only
        # transfer: the legacy path applied an unapportioned mu1-only drain
        # (run-10 drove pool polypropylene_mod to mu0=+0.0047 / mu1=-5.6e-5
        # at simulate-leg 21; the r81 negative-beyond-floor tripwire killed
        # the run). The refusal is re-derived EVERY rebuild and deliberately
        # NOT stamped on the reaction object: once item 16 configures the
        # spawned pool, the SAME row resolves and runs fully apportioned.
        # Unresolved endpoints NOT attributable to a spawned config-less
        # Polymer keep the legacy demotion unchanged, censused separately as
        # a loud anomaly (should be impossible for configured pools; the
        # constructible routes are a non-Polymer daughter species, an
        # edge-only proxy, or a duplicate core base label).
        from rmgpy.polymer import Polymer as _Polymer
        from rmgpy.polymer import strip_rmg_index_suffix as _strip_idx
        _configured_pool_labels = {p.label for p in self.polymer_pools}
        _arch_names = {FLUX_MIGRATION: "MIGRATION",
                       FLUX_SCISSION_FRAGMENT: "SCISSION_FRAGMENT",
                       FLUX_DISCRETE_CHIP: "DISCRETE_CHIP",
                       FLUX_VOLATILE_EJECTION: "VOLATILE_EJECTION"}
        spawned_refusals = []      # (row_idx, rxn, archetype name, labels)
        demotion_anomalies = []    # (row_idx, rxn, archetype name)
        self.spawned_pool_refusal_census = []
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
                if self.reaction_refused[i]:
                    # Refused row (adjudicated, flux-dead via
                    # reaction_refused): stays archetype NONE -- never
                    # routed onto the legacy mu1-only path (r71 FIX 4).
                    pass
                elif self._allow_unstamped_proxy_rows:
                    self.reaction_flux_archetype[i] = FLUX_UNRESOLVED
                    n_unstamped += 1
                else:
                    unstamped_live_rows.append(
                        "%s row %d: %s" % (
                            "core" if i < len(core_reactions) else "edge",
                            i, str(combined_rxns[i])))
            arch_i = self.reaction_flux_archetype[i]
            if ((arch_i in (FLUX_MIGRATION, FLUX_SCISSION_FRAGMENT, FLUX_VOLATILE_EJECTION)
                    and (src == -1 or dst == -1))
                    or (arch_i == FLUX_DISCRETE_CHIP
                        and src == -1)):
                # A stamped archetype needs its pool(s) resolved in the
                # solver: MIGRATION/SCISSION_FRAGMENT/VOLATILE_EJECTION need
                # BOTH src and dst
                # (e.g. scission daughters are registered as core Polymer
                # species but have no pool config yet); DISCRETE_CHIP needs
                # only src (no dst: complement folds back to the same pool
                # and the chip is a plain gas species).
                # src == dst is deliberately NOT touched: for that shape
                # (fold-back proxy product + non-pool daughter) the dispatch
                # skip and the legacy transfer produce identical pool flux
                # (reactant -r and fold-back +r cancel on the same mu1).
                #
                # r91: attribute each MISSING required endpoint. If every
                # missing endpoint is explained by a spawned config-less
                # Polymer participant on that side, REFUSE the row
                # (conduit-deferred, zero flux -- archetype cleared to NONE
                # so no dispatch/census path sees a live stamped shape; the
                # generation stamp survives in
                # reaction_pre_demotion_archetype). Otherwise keep the
                # legacy mu1-only demotion (parent drain never silently
                # zeroed) and census the anomaly loudly.
                rxn_i = combined_rxns[i]
                spawned_labels = set()
                attributable = True
                sides = []
                if src == -1:
                    sides.append(rxn_i.reactants or [])
                if dst == -1 and arch_i != FLUX_DISCRETE_CHIP:
                    sides.append(rxn_i.products or [])
                for side in sides:
                    side_labels = [
                        _strip_idx(s.label) for s in side
                        if isinstance(s, _Polymer)
                        and _strip_idx(s.label) not in _configured_pool_labels]
                    if side_labels:
                        spawned_labels.update(side_labels)
                    else:
                        attributable = False
                if attributable and spawned_labels:
                    self.reaction_refused[i] = 1
                    self.reaction_flux_archetype[i] = FLUX_NONE
                    spawned_refusals.append(
                        (i, rxn_i, _arch_names[arch_i],
                         sorted(spawned_labels)))
                else:
                    self.reaction_flux_archetype[i] = FLUX_UNRESOLVED
                    n_unresolvable += 1
                    demotion_anomalies.append((i, rxn_i, _arch_names[arch_i]))
        if spawned_refusals:
            self.spawned_pool_refusal_census = [
                {"reaction": _reaction_census_label(r_rxn),
                 "archetype": r_arch,
                 "pools": r_pools,
                 "reason": "conduit-deferred"}
                for (_, r_rxn, r_arch, r_pools) in spawned_refusals]
            logging.warning(
                "SPAWNED-POOL DEMOTION REFUSAL: %d stamped pool-coupled "
                "row(s) targeting unconfigured spawned pools refused "
                "conduit-deferred instead of demoted to legacy mu1-only; "
                "item-16 pending; archetypes=%s, pools=%s, first_rows=%s",
                len(spawned_refusals),
                ",".join(sorted({a for (_, _, a, _) in spawned_refusals})),
                ",".join(sorted({lbl for (_, _, _, ls) in spawned_refusals
                                 for lbl in ls})),
                "; ".join("%s row %d: %s" % (
                    "core" if r_i < len(core_reactions) else "edge",
                    r_i, str(r_rxn))
                    for (r_i, r_rxn, _, _) in spawned_refusals[:3]))
        if n_unresolvable:
            logging.warning(
                "%d reactions stamped MIGRATION/SCISSION_FRAGMENT/DISCRETE_CHIP "
                "could not resolve their solver pool(s); demoted to legacy "
                "mu1-only moment flux (UNRESOLVED).", n_unresolvable)
            logging.warning(
                "CONFIGURED-POOL UNRESOLVED DEMOTION ANOMALY: %d stamped "
                "pool-coupled row(s) had an unresolved required pool "
                "endpoint NOT attributable to an unconfigured spawned pool "
                "(should be impossible for configured pools; constructible "
                "routes: non-Polymer daughter species, edge-only proxy, "
                "duplicate core base label); legacy mu1-only demotion kept "
                "unchanged. rows=%s",
                len(demotion_anomalies),
                "; ".join("%s row %d [%s]: %s" % (
                    "core" if a_i < len(core_reactions) else "edge",
                    a_i, a_arch, str(a_rxn))
                    for (a_i, a_rxn, a_arch) in demotion_anomalies[:5]))
        if unstamped_live_rows:
            raise ValueError(
                "UNSTAMPED PROXY ROW(S): %d proxy-touching reaction(s) "
                "arrived at the solver rebuild with no polymer_flux_archetype "
                "stamp and no polymer_refused adjudication -- the legacy "
                "mu1-only fallback would apply silent unclassified flux to "
                "the pool (the PP run-5 stall channel, adjudicated r71). "
                "Offending row(s): %s. A production build must classify "
                "every live proxy row at generation "
                "(stamp_polymer_flux_archetype) or refuse it "
                "(stamp_gas_association_refusal). Direct-test/runner "
                "construction that deliberately exercises the legacy "
                "mu1-only fallback must set allow_unstamped_proxy_rows=True."
                % (len(unstamped_live_rows), "; ".join(unstamped_live_rows)))
        if n_unstamped:
            logging.warning(
                "%d proxy-touching reactions arrived without a polymer_flux_archetype "
                "stamp; applying legacy mu1-only pool moment flux for them "
                "(allow_unstamped_proxy_rows escape hatch -- direct-test/"
                "runner construction only).",
                n_unstamped)

        # item 18 (T2): double-count tripwire. A pool carrying BOTH a surviving
        # explicit beta-scission/chip reaction sourced from it (archetype
        # SCISSION_FRAGMENT or DISCRETE_CHIP) AND a nonzero phenomenological
        # k_scission/k_unzip double-counts chain degradation: the explicit
        # chemistry and the phenomenological stand-in both model the same chain
        # break. Census it (correct-but-loud, warn-once). Placement is
        # load-bearing: AFTER the demotion loop completes, so it scans the FINAL
        # archetypes -- a reaction demoted to UNRESOLVED above is no longer
        # explicit scission and must NOT trip the tripwire. Diagnostic-only:
        # changes NO flux/state. Default severity is census-loud (warn only; NO
        # refuse-cliff -- severity is calibrated later by item 19).
        self.double_count_census = []
        scission_src_pools = set()
        for r_dc in range(n_rxn):
            if self.reaction_flux_archetype[r_dc] in (FLUX_SCISSION_FRAGMENT, FLUX_DISCRETE_CHIP):
                if self.reaction_src_pool[r_dc] >= 0:
                    scission_src_pools.add(self.reaction_src_pool[r_dc])
        for p_idx, pool in enumerate(self.polymer_pools):
            if (pool.k_scission > 0.0 or pool.k_unzip > 0.0) and p_idx in scission_src_pools:
                entry = {"pool": pool.label, "k_scission": pool.k_scission, "k_unzip": pool.k_unzip}
                self.double_count_census.append(entry)
                _warn_once_double_count(entry)

        # Enforce the moment-isolation invariant and pool/mass-transfer index
        # sanity now that gas_species_mask and the reaction index tables are
        # populated. Moments must evolve only via the tail block, never through
        # generic reaction stoichiometry.
        # RIDER R2 static half (spec SS3(e)): once per rebuild by
        # construction -- the enumeration runs exactly once per
        # initialize_model; no keying set needed.
        self._static_phase_gate_census(core_species, core_reactions)

        self.validate_configuration()

        # radical_qssa_unzip double-count census (M2). QSSA initiation is
        # backbone homolysis, so generated-chemistry SCISSION_FRAGMENT /
        # VOLATILE_EJECTION reactions sourced from the same pool may
        # represent overlapping physics. WARN-ONCE census mirroring the
        # item-18 tripwire above -- warn, NEVER refuse (generated scission
        # may cover different bonds; hard exclusion would be wrong). The
        # k_unzip co-presence stays a hard error (M1 mutual exclusion).
        # Placement is load-bearing twice over: AFTER the demotion loop (it
        # scans FINAL archetypes, like the item-18 census) and AFTER
        # validate_configuration (qssa_enabled -- the flattened solver-owned
        # gate, the ONLY channel signal the solver trusts -- exists here).
        from rmgpy.polymer import warn_once_qssa_double_count
        self.qssa_double_count_census = []
        qssa_overlap_src_pools = set()
        for r_dc in range(self.reaction_flux_archetype.shape[0]):
            if self.reaction_flux_archetype[r_dc] in (FLUX_SCISSION_FRAGMENT, FLUX_VOLATILE_EJECTION):
                if self.reaction_src_pool[r_dc] >= 0:
                    qssa_overlap_src_pools.add(self.reaction_src_pool[r_dc])
        for p_idx, pool in enumerate(self.polymer_pools):
            if not self.qssa_enabled[p_idx]:
                continue
            if p_idx in qssa_overlap_src_pools:
                entry = {"pool": pool.label,
                         "overlap": "generated_scission_ve"}
                self.qssa_double_count_census.append(entry)
                warn_once_qssa_double_count(entry)
            # QSSA + the pool's own k_scission: BOTH are random backbone
            # homolysis -- the most direct initiation double-count of all
            # (unlike generated scission, k_scission cannot claim to cover
            # different bonds by construction; only the user knows whether
            # the two were parameterized for disjoint physics). Same
            # warn-once helper, distinct census key.
            if pool.k_scission > 0.0:
                entry = {"pool": pool.label, "overlap": "k_scission",
                         "k_scission": pool.k_scission}
                self.qssa_double_count_census.append(entry)
                warn_once_qssa_double_count(entry)

        # k_homolysis supersession census (Stage 1, adjudicated round 66):
        # refused explicit homolysis/association rows co-existing with the
        # kernel are FINE -- they carry zero flux (reaction_refused) -- but
        # the census must state the pairing EXPLICITLY: a warn-level line
        # names the refused conduit-deferred rows as superseded by the
        # kernel. Placement mirrors the QSSA census above (after
        # validate_configuration: khom_enabled -- the flattened solver-owned
        # gate -- exists here; after the refused census: superseded rows are
        # final). Diagnostic-only: NO flux/state change, warn NEVER refuse.
        from rmgpy.polymer import warn_once_homolysis_supersession
        self.homolysis_supersession_census = []
        deferred_rows = [c["reaction"] for c in self.refused_reaction_census
                         if c.get("reason") == "conduit-deferred"]
        if deferred_rows and np.any(self.khom_enabled):
            for p_idx, pool in enumerate(self.polymer_pools):
                if not self.khom_enabled[p_idx]:
                    continue
                entry = {"pool": pool.label,
                         "superseded_rows": list(deferred_rows)}
                self.homolysis_supersession_census.append(entry)
                warn_once_homolysis_supersession(entry)

        # side_group_homolysis supersession census (FR1-K1): refused
        # explicit gas-radical<->condensed rows REMAIN refused next to the
        # live kernel (they carry zero flux -- reaction_refused), and the
        # census pairs them explicitly, same contract as the k_homolysis
        # census above. Diagnostic-only: NO flux/state change, warn NEVER
        # refuse. Same placement rationale (after validate_configuration:
        # sgh_enabled exists; after the refused census: rows are final).
        from rmgpy.polymer import warn_once_side_group_supersession
        self.side_group_supersession_census = []
        if deferred_rows and np.any(self.sgh_enabled):
            for p_idx, pool in enumerate(self.polymer_pools):
                if not self.sgh_enabled[p_idx]:
                    continue
                entry = {"pool": pool.label,
                         "superseded_rows": list(deferred_rows)}
                self.side_group_supersession_census.append(entry)
                warn_once_side_group_supersession(entry)

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
        self._scratch_du_dt = np.zeros(self.num_qssa_u, float)
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

        # Exhaustion-tail conditioning (r81 B): per-pool per-moment floors
        # derived from the solver ABSOLUTE tolerances -- floor_k =
        # max(SMALL_EPS, EXHAUSTION_FLOOR_K * atol[state]) on the pool's own
        # mu0/mu1/mu2 state slots (mole basis, same basis atol governs).
        # Census keying is per-rebuild, mirroring the sibling warn-once sets.
        n_pools_exh = len(self.polymer_pools)
        self._pool_exhausted = np.zeros(n_pools_exh, dtype=np.uint8)
        self._exhaustion_census_emitted = set()
        self._pool_mu_floors = np.full((n_pools_exh, 3), SMALL_EPS,
                                       dtype=float)
        atol_arr_exh = getattr(self, "atol_array", None)
        for p in range(n_pools_exh):
            mu_idx_exh = self.polymer_pools[p].mu_indices
            for k in range(3):
                s_idx = int(mu_idx_exh[k])
                if atol_arr_exh is None:
                    a_exh = 0.0
                elif np.ndim(atol_arr_exh) == 0:
                    a_exh = float(atol_arr_exh)
                elif s_idx < len(atol_arr_exh):
                    a_exh = float(atol_arr_exh[s_idx])
                else:
                    a_exh = float(np.max(atol_arr_exh))
                self._pool_mu_floors[p, k] = max(
                    SMALL_EPS, EXHAUSTION_FLOOR_K * a_exh)

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
            # k_homolysis kernel + end-radical daughter pools must be
            # traceable in run logs (round 66: radical inventories); the
            # daughters themselves enumerate here as ordinary pools, with
            # their role tagged so a reader can find their mu0/mu1 growth.
            khom = pool.k_homolysis
            if khom is not None:
                khom_str = (f"A={khom['A']:.2e},n={khom['n']:g},"
                            f"Ea={khom['Ea']:.4g}J/mol")
            else:
                khom_str = "off"
            # side_group_homolysis channels (FR1-K1) must be traceable in
            # run logs too: the parent names every channel; feature pools
            # are tagged; the v1 saturation limitation is disclosed.
            sgh_channels = pool.side_group_homolysis
            if sgh_channels:
                sgh_str = ("[" + ",".join(ch["label"] for ch in sgh_channels)
                           + "]")
            else:
                sgh_str = "off"
            # k_depropagation kernel (r74 SS2): the radical-end consumption
            # channel must be traceable in run logs too.
            kdep = pool.k_depropagation
            if kdep is not None:
                kdep_str = (f"A={kdep['A']:.2e},n={kdep['n']:g},"
                            f"Ea={kdep['Ea']:.4g}J/mol")
            else:
                kdep_str = "off"
            role_tag = ""
            if pool.label.endswith(K_HOMOLYSIS_DAUGHTER_SUFFIXES):
                role_tag = "  [end-radical daughter]"
            elif SIDE_GROUP_DAUGHTER_INFIX in pool.label:
                role_tag = "  [side-group X-loss daughter]"
            print(f"\n--- Pool {pool_i}: '{pool.label}'  (xs={pool.xs}, "
                  f"k_scission={pool.k_scission:.2e}, "
                  f"k_unzip={pool.k_unzip:.2e}, "
                  f"k_homolysis={khom_str}, "
                  f"k_depropagation={kdep_str}, "
                  f"side_group_homolysis={sgh_str}){role_tag} ---")
            if sgh_channels:
                print(f"  side_group_homolysis v1 saturation: the kernel "
                      f"acts on THIS parent pool only; its X-loss feature "
                      f"pools carry no further side-group loss (no "
                      f"multi-loss cascade -- they saturate as terminal "
                      f"X-loss sinks).")
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

        # --- Phase / Volume summary (DEBUG) ---
        n_gas_species = int(np.sum(self.gas_species_mask))
        n_poly_species = len(core_species) - n_gas_species
        n_gas_moles = float(np.sum(self.y0[:len(core_species)][self.gas_species_mask]))
        n_poly_moles = float(np.sum(self.y0[:len(core_species)][~self.gas_species_mask]))
        n_gas_via_residual = float(np.sum(self.y0[:len(core_species)] * self.gas_species_mask))
        print(f"\n--- Phase / Volume summary (DEBUG) ---")
        print(f"  T.value_si = {self.T.value_si}  P.value_si = {self.P.value_si}  "
              f"R = {constants.R}  constant_gas_volume = {self.constant_gas_volume}  V_gas0 = {self.V_gas0}")
        print(f"  n_gas (mask*y0 sum) = {n_gas_via_residual:.6e} mol   "
              f"ideal V_gas = {constants.R*self.T.value_si*n_gas_via_residual/self.P.value_si:.6e} m^3")
        print(f"  V_gas = {self.V_gas:.6e} m^3   V_poly = {self.V_poly:.6e} m^3   V = {self.V:.6e} m^3")
        print(f"  gas species: {n_gas_species} (sum y0 = {n_gas_moles:.6e} mol)   "
              f"poly species: {n_poly_species} (sum y0 = {n_poly_moles:.6e} mol)")
        print(f"  {'GAS Species':<30} {'Index':<7} {'y0':>13}")
        for i in range(len(core_species)):
            if self.gas_species_mask[i] and self.y0[i] != 0.0:
                print(f"  {core_species[i].label:<30} {i:<7} {self.y0[i]:>13.4e}")
        # largest-magnitude initial gas concentration (collapse check)
        if self.V_gas > 0:
            cmax = 0.0
            cmax_lbl = "-"
            for i in range(len(core_species)):
                if self.gas_species_mask[i]:
                    c = abs(self.y0[i]) / self.V_gas
                    if c > cmax:
                        cmax = c
                        cmax_lbl = core_species[i].label
            print(f"  max |C_gas0| = {cmax:.4e} mol/m^3  ({cmax_lbl})")

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

        # 7. Weak-link U-state slots (milestones ii+iii): U(0) =
        # unsaturated_tail_ends_initial [mol -- SAME amount basis as mu0],
        # behind the CENSUS TRAP (adversarial-review requirement): each
        # chain carries at most 2 tail ends, so U0 must fit the TAIL-
        # DISTRIBUTION chain-end capacity 2*mu0_tail (the tail mu0 written
        # back in step 6). TAIL-ONLY basis (review r37): U is a tail-
        # distribution state -- B and the RHS capacity throttle count the
        # tail moments only -- so the census must use the SAME basis;
        # explicit-species ends must not back a U0 that would then evolve
        # against tail-only capacity (mixed semantics). A mol vs mol/L (or
        # other per-volume) typo would otherwise become a silent hidden
        # initiation source. U0 = 0 always passes (mu0 may be 0 for an
        # empty pool).
        u_slots = getattr(self, "qssa_u_slot", None)
        if u_slots is not None:
            for i, pool in enumerate(self.polymer_pools):
                slot = int(u_slots[i])
                if slot < 0:
                    continue
                u0 = float(self.qssa_u0[i])
                mu0_amount = max(0.0, self.y0[pool.mu_indices[0]])
                capacity = 2.0 * mu0_amount
                if u0 > capacity:
                    raise ValueError(
                        f"Pool {pool.label}: radical_qssa_unzip "
                        f"unsaturated_tail_ends_initial={u0:g} mol exceeds "
                        f"the pool's tail-distribution chain-end capacity "
                        f"2*mu0 = {capacity:g} mol (initial TAIL mu0 = "
                        f"{mu0_amount:g} mol; each chain carries at most 2 "
                        f"tail ends, and U is a TAIL-distribution state -- "
                        f"explicit-species chain ends do not back U, "
                        f"matching the RHS capacity throttle basis). U is "
                        f"on the SAME amount basis as mu0 [mol] -- a mol "
                        f"vs mol/L (or other per-volume) typo here would "
                        f"silently become a hidden initiation source. Fix "
                        f"the amount basis or the value.")
                self.y0[slot] = u0

    def generate_rate_coefficients(self, core_reactions, edge_reactions):
        for rxn in itertools.chain(core_reactions, edge_reactions):
            j = self.reaction_index[rxn]
            kf = rxn.get_rate_coefficient(self.T.value_si, self.P.value_si)
            # Guard against non-finite forward coefficients (e.g. extreme PLOG/Troe
            # extrapolation) so a single bad reaction cannot poison the residual.
            self.kf[j] = kf if np.isfinite(kf) else 0.0
            if rxn.reversible:
                Keq = rxn.get_equilibrium_constant(self.T.value_si)
                # When Keq under/overflows (thermo at extreme T for large fragments),
                # kb = kf/Keq -> inf, and inf * 0 (zero product conc) = NaN, which
                # crashes DASPK ("nans in moles"). Treat such degenerate reverse
                # directions as irreversible rather than emitting a non-finite kb.
                if (not np.isfinite(Keq)) or Keq <= 0.0:
                    self.Keq[j] = np.inf
                    self.kb[j] = 0.0
                else:
                    kb = self.kf[j] / Keq
                    self.Keq[j] = Keq
                    self.kb[j] = kb if np.isfinite(kb) else 0.0
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
        2. ``pool_stats``: dict pool label ->
           ``(E_n, monomer_mw_g_mol, chain_mass_defect_g_mol)``
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
           ``gross_j * max(0, condensed_mass_g(1, E_n[pool(j)]))``
           = ``gross_j * max(0, E_n*mw - chain_mass_defect)`` over
           CANONICAL pool proxies
           (``species_to_pool_indices[j] != -1 and is_pool_proxy[j]``
           — the engine CAN attribute those). The defect term is the
           FR1-K2 mass-consumer audit (round-72 P2): an X-loss feature
           pool's chains each lost exactly one X, so the EXACT per-chain
           event mass is the pool's normative condensed_mass_g at
           (mu0=1, mu1=E_n) — the raw E_n*mw would overstate every event
           by M_X. Ordinary pools (defect 0) reduce to the legacy
           product bit-identically; the max(0, .) clamp errs toward
           deferral.

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
            # r81 (B) negative-moment rule: the pre-r81 mu0 <= SMALL_EPS
            # branch zeroed SILENTLY and would swallow a negative. Now a
            # mu0 beyond -floor is a hard error (integrator corruption);
            # inside [-floor, SMALL_EPS] it projects to E[n] = 0 WITH the
            # exhaustion census below. Floor from the centralized
            # atol-derived exhaustion floors when available
            # (post-initialize_model); SMALL_EPS otherwise (pre-init
            # snapshot -> honest degenerate, mirroring atol_mu0 = 0).
            floors_exh = self._pool_mu_floors
            floor0 = (float(floors_exh[p, 0]) if floors_exh is not None
                      and p < len(floors_exh) else SMALL_EPS)
            if mu0 < -floor0:
                raise ValueError(
                    f"Polymer pool '{self.polymer_pools[p].label}' mu0="
                    f"{mu0!r} mol is negative beyond the exhaustion floor "
                    f"-{floor0:.6e} mol (= -max(SMALL_EPS, "
                    f"{EXHAUSTION_FLOOR_K:.0f}*atol[state])) at the "
                    f"spawn-gate snapshot. A negative moment beyond the "
                    f"integrator's own error budget is integrator "
                    f"corruption, not exhaustion -- refusing to attribute "
                    f"from it (r81 negative-moment rule: hard error, never "
                    f"a silent 0.0).")
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
            else:
                # mu0 in [-floor0, SMALL_EPS]: exhausted/empty -> E[n] = 0.0
                # WITH the exhaustion census (r81 B5, once per pool per
                # rebuild) -- never silently: a swallowed in-band negative
                # would hide the onset of integrator corruption.
                key_exh = (self.polymer_pools[p].label, "spawn-gate")
                if key_exh not in self._exhaustion_census_emitted:
                    self._exhaustion_census_emitted.add(key_exh)
                    logging.warning(
                        "POOL EXHAUSTION CENSUS: pool %s spawn-gate "
                        "mu0=%.6e mol is within the exhaustion band "
                        "[-%.6e, %.0e]; E[n] attribution projected to 0.0 "
                        "for this snapshot (floor = max(SMALL_EPS, "
                        "%.0f*atol[mu0]); negative beyond -floor is a "
                        "hard error).",
                        self.polymer_pools[p].label, mu0, floor0, SMALL_EPS,
                        EXHAUSTION_FLOOR_K)
        for p in range(n_pools):
            pool = self.polymer_pools[p]
            mw = float(getattr(pool, "monomer_mw_g_mol", 0.0) or 0.0)
            defect = float(getattr(pool, "chain_mass_defect_g_mol", 0.0)
                           or 0.0)
            pool_stats[pool.label] = (e_n_by_pool[p], mw, defect)
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
                pool = self.polymer_pools[p]
                e_n = pool_stats[pool.label][0]
                # EXACT per-chain event mass (FR1-K2 mass audit, round-72
                # P2): the pool's normative condensed_mass_g at
                # (mu0=1, mu1=E_n) = E_n*mw - chain_mass_defect. Ordinary
                # pools (defect 0) reduce to the legacy E_n*mw
                # bit-identically; max(0, .) errs toward deferral.
                proxy_event_mass_total += g * max(
                    0.0, pool.condensed_mass_g(1.0, e_n))
        return gross, pool_stats, proxy_event_mass_total

    def _update_pool_exhaustion(self, y):
        """Centralized pool-exhaustion CENSUS/TRIPWIRE driver (r81 B, the
        PP run-7 IDID=-7 wall; re-adjudicated r86: the r81 RHS
        read-projection is reverted -- ``self._pool_exhausted`` is kept
        for DIAGNOSTIC LOGGING ONLY and has ZERO RHS effect; the RHS
        always computes from the raw state vector). Called once per
        residual evaluation; refreshes ``self._pool_exhausted`` from the
        RAW state vector ``y``:

        * a pool is counted EXHAUSTED (census once per pool per rebuild)
          only when ``abs(mu0)``, ``abs(mu1)`` AND ``abs(mu2)`` are ALL
          below their per-state floors
          ``max(SMALL_EPS, EXHAUSTION_FLOOR_K * atol[state])`` (never mu0
          alone: tiny mu0 with nontrivial mu1 is a few long chains);
        * any moment more negative than ``-floor_k`` is a HARD error
          (integrator corruption, never ``max(..., SMALL_EPS)``) -- this
          tripwire is unchanged from r81;
        * a negative inside ``[-floor_k, +floor_k]`` is ANNOUNCED through
          the exhaustion census instead of silence (the RHS's
          pre-existing local ``max(0.0, .)`` reads are what clamp it).

        ``y`` itself is never mutated. mu3 < mu2 OUTSIDE the exhaustion
        band deliberately stays a real guard event (the existing warn-once
        kernel guards are untouched)."""
        cdef int p, n_pools
        cdef double f0, f1, f2, raw0, raw1, raw2
        n_pools = len(self.polymer_pools)
        for p in range(n_pools):
            idx0, idx1, idx2 = self.polymer_pools[p].mu_indices
            f0 = self._pool_mu_floors[p, 0]
            f1 = self._pool_mu_floors[p, 1]
            f2 = self._pool_mu_floors[p, 2]
            raw0 = y[idx0]
            raw1 = y[idx1]
            raw2 = y[idx2]
            if raw0 < -f0 or raw1 < -f1 or raw2 < -f2:
                pool = self.polymer_pools[p]
                raise ValueError(
                    f"Polymer pool '{pool.label}' moment state went "
                    f"negative beyond the exhaustion floor: mu0={raw0!r}, "
                    f"mu1={raw1!r}, mu2={raw2!r} mol vs floors "
                    f"-{f0:.6e}/-{f1:.6e}/-{f2:.6e} mol "
                    f"(= -max(SMALL_EPS, {EXHAUSTION_FLOOR_K:.0f}*"
                    f"atol[state])). A negative moment beyond the "
                    f"integrator's own error budget is integrator "
                    f"corruption, not exhaustion -- refusing to integrate "
                    f"it (r81 negative-moment rule: hard error, never "
                    f"max(..., SMALL_EPS)).")
            if (abs(raw0) <= f0 and abs(raw1) <= f1 and abs(raw2) <= f2):
                self._pool_exhausted[p] = 1
                self._emit_pool_exhaustion_census(
                    p, raw0, raw1, raw2,
                    "exhausted (census only -- r86: RHS computes live "
                    "from the raw state)")
            else:
                self._pool_exhausted[p] = 0
                if raw0 < 0.0 or raw1 < 0.0 or raw2 < 0.0:
                    self._emit_pool_exhaustion_census(
                        p, raw0, raw1, raw2,
                        "in-band negative moment (census only; pool NOT "
                        "exhausted)")

    def _emit_pool_exhaustion_census(self, int p, double raw0, double raw1,
                                     double raw2, reason):
        """Once-per-pool-per-rebuild exhaustion census (r81 B3,
        re-adjudicated r86: log-only diagnostics -- nothing is suppressed):
        label, raw moments, floors, and which kernels / pool-coupled
        reaction rows draw on the pool while it sits in the census band.
        Cold path -- string and row-scan work happen at most once per
        (pool, reason) per rebuild."""
        pool = self.polymer_pools[p]
        key = (pool.label, reason)
        if key in self._exhaustion_census_emitted:
            return
        self._exhaustion_census_emitted.add(key)
        f0 = self._pool_mu_floors[p, 0]
        f1 = self._pool_mu_floors[p, 1]
        f2 = self._pool_mu_floors[p, 2]
        kernels = []
        if getattr(pool, "tail_kinetics", None):
            kernels.append("tail_kinetics")
        if float(getattr(pool, "k_scission", 0.0) or 0.0) > 0.0:
            kernels.append("k_scission")
        if float(getattr(pool, "k_unzip", 0.0) or 0.0) > 0.0:
            kernels.append("k_unzip")
        khom = getattr(self, "khom_enabled", None)
        if khom is not None and len(khom) > p and khom[p]:
            kernels.append("k_homolysis")
        sgh = getattr(self, "sgh_enabled", None)
        if sgh is not None and len(sgh) > p and sgh[p]:
            kernels.append("side_group_homolysis")
        qssa = getattr(self, "qssa_enabled", None)
        if qssa is not None and len(qssa) > p and qssa[p]:
            kernels.append("radical_qssa_unzip")
        kdep = getattr(self, "kdep_enabled", None)
        if kdep is not None and len(kdep) > p and kdep[p]:
            kernels.append("k_depropagation")
        coupled_rows = []
        n_rows = 0
        ir = getattr(self, "reactant_indices", None)
        ip = getattr(self, "product_indices", None)
        stp = getattr(self, "species_to_pool_indices", None)
        if ir is not None and ip is not None and stp is not None:
            n_core = self.num_core_species
            n_rows = ir.shape[0]
            for j in range(n_rows):
                hit = False
                for slot in range(3):
                    a = ir[j, slot]
                    b = ip[j, slot]
                    if ((a != -1 and a < n_core and stp[a] == p)
                            or (b != -1 and b < n_core and stp[b] == p)):
                        hit = True
                        break
                if hit:
                    coupled_rows.append(j)
        logging.warning(
            "POOL EXHAUSTION CENSUS: pool %s %s (raw mu0=%.6e, mu1=%.6e, "
            "mu2=%.6e mol; floors %.6e/%.6e/%.6e mol = max(SMALL_EPS, "
            "%.0f*atol[state])); pool kernels drawing on this state: %s; "
            "pool-coupled reaction rows: %d of %d (core+edge)%s.",
            pool.label, reason, raw0, raw1, raw2, f0, f1, f2,
            EXHAUSTION_FLOOR_K, kernels if kernels else "none",
            len(coupled_rows), n_rows,
            " (first rows: %s)" % coupled_rows[:8]
            if coupled_rows else "")

    def get_total_polymer_condensed_mass_g(self, y=None):
        """r86 terminationPolymerConversion metric (defect-adjusted):

            M_pool = max(0, mu1_p*monomer_mw_g_mol_p
                            - mu0_p*chain_mass_defect_g_mol_p)
            M_poly = sum over ALL solver polymer pools of M_pool

        summed over EVERY pool the solver carries at this rebuild --
        configured, spawned, homolysis-daughter, deprop-daughter AND
        side-group feature pools (defect pools apply the schema-2.7
        mass-defect formula through PolymerPoolConfig.condensed_mass_g, so
        an FR1 Br-loss pool never falsely retains Br mass). Explicit-DP
        oligomer carriers hold genuine polymer repeat-unit mass OUTSIDE
        the mu state slots (the init consistency check subtracts their
        share from mu1, r86 probe 3), so their contribution is added back
        via _explicit_moment_contributions. Only each pool's FINAL mass
        contribution is clamped at zero (numerical fuzz); individual
        moments are never clamped. Returns grams; reads self.y when no
        state vector is passed."""
        cdef int p, i0, i1
        cdef double total, mu0_mol, mu1_mol, V_poly
        if y is None:
            y = self.y
        if self.polymer_pools is None or len(self.polymer_pools) == 0:
            return None
        V_poly = self.V_poly
        total = 0.0
        for p in range(len(self.polymer_pools)):
            pool = self.polymer_pools[p]
            i0 = self.pool_mu0_indices[p] if self.pool_mu0_indices is not None else -1
            i1 = self.pool_mu1_indices[p] if self.pool_mu1_indices is not None else -1
            if i0 < 0 or i1 < 0:
                # Pre-initialize_model fallback: the config-time layout.
                i0, i1 = pool.mu_indices[0], pool.mu_indices[1]
            mu0_mol = y[i0]
            mu1_mol = y[i1]
            if pool.explicit_dp_to_species_index:
                e0, e1, _e2 = _explicit_moment_contributions(
                    V_poly, y, pool.explicit_dp_to_species_index)
                mu0_mol += e0 * V_poly
                mu1_mol += e1 * V_poly
            total += max(0.0, pool.condensed_mass_g(mu0_mol, mu1_mol))
        return total

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
        if (self._scratch_du_dt is None
                or len(self._scratch_du_dt) != self.num_qssa_u):
            self._scratch_du_dt = np.zeros(self.num_qssa_u, float)

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
        # NOTE the [:n_core] slice: y may carry trailing weak-link U slots
        # beyond the core-species block; the mask covers only the prefix.
        C_gas[mask] = np.maximum(0.0, y[:n_core][mask]) / V_gas
        C_poly[~mask] = np.maximum(0.0, y[:n_core][~mask]) / V_poly

        # Sync diagnostic concentrations
        self.core_species_concentrations[mask] = C_gas[mask]
        self.core_species_concentrations[~mask] = C_poly[~mask]

        # 4. Clear Accumulators
        dn_dt = self._scratch_dn_dt
        dn_dt[:] = 0.0

        du_dt = self._scratch_du_dt
        du_dt[:] = 0.0

        proxy_activity = self._scratch_proxy_activity
        proxy_activity[:] = 0.0

        # r81 (B) exhaustion census/tripwire, re-adjudicated r86: refresh
        # the per-pool exhausted flags from the RAW state ONCE per residual
        # evaluation (hard-errors on any moment beyond -floor; census once
        # per pool per rebuild). The flags are DIAGNOSTIC ONLY -- no RHS
        # read below consults them (the r81 read-projection is reverted);
        # site factors, bundles and kernels always compute from raw y.
        self._update_pool_exhaustion(y)

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

                    kf_j = self.pdep_collider_kinetics[i].get_rate_coefficient(T_val, Peff)
                    self.kf[j] = kf_j if np.isfinite(kf_j) else 0.0
                    if self.Keq[j] != np.inf and self.Keq[j] > 0.0:
                        kb_j = self.kf[j] / self.Keq[j]
                        self.kb[j] = kb_j if np.isfinite(kb_j) else 0.0
            else:
                pass

        # 6. Standard Kinetics (Optimized Loop)
        for r_idx in range(ir.shape[0]):
            r0, r1, r2 = ir[r_idx, 0], ir[r_idx, 1], ir[r_idx, 2]

            if r0 == -1:
                continue

            # simple.pyx parity (:443-459, :509-515): a row with an EDGE
            # REACTANT is still evaluated -- the edge species has no state in
            # y, so its concentration is ZERO and rf vanishes; the discovery
            # flux for reverse-stored decomposition channels (daughter as
            # edge reactant, rr = kb*[core products]) flows through the
            # reactant-side edge accumulator in section 4b. Phase verdict for
            # an edge reactant comes from prospective_gas_mask (the only mask
            # covering edge slots; its core prefix is bit-identical to
            # gas_species_mask by rider R1).
            r0_is_gas = mask[r0] if r0 < n_core else pmask[r0]
            r1_is_gas = True
            r2_is_gas = True

            if r1 != -1:
                r1_is_gas = mask[r1] if r1 < n_core else pmask[r1]

            if r2 != -1:
                r2_is_gas = mask[r2] if r2 < n_core else pmask[r2]

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
            #    Edge reactants map to -1: species_to_pool_indices covers the
            #    core prefix only, and an edge species is never a pool proxy.
            p0_pool_idx = -1 if r0 >= n_core else self.species_to_pool_indices[r0]
            p1_pool_idx = -1 if (r1 == -1 or r1 >= n_core) else self.species_to_pool_indices[r1]
            p2_pool_idx = -1 if (r2 == -1 or r2 >= n_core) else self.species_to_pool_indices[r2]  # optional but consistent

            # Rate Calculation
            kf = self.kf[r_idx]
            kb = self.kb[r_idx]

            def _C(idx):
                if idx == -1:
                    return 1.0
                if idx >= n_core:
                    # Edge species: no state in y -> zero concentration, so
                    # rf vanishes for edge-reactant rows (simple.pyx parity;
                    # NOT a guard -- definitional, no silent zeroing).
                    return 0.0
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
                    elif (self.reaction_flux_archetype[r_idx] == FLUX_VOLATILE_EJECTION
                            and self.is_end_group_reaction[r_idx]
                            and self.reaction_src_pool[r_idx] != -1
                            and self.reaction_src_pool[r_idx] == self.reaction_dst_pool[r_idx]
                            and self.reaction_eject_units[r_idx] > 0.0):
                        # Same-pool a>0 VE shares the DISCRETE_CHIP exhaustion
                        # structure: mu0-scaled drain of mu1 that never touches
                        # mu0, so it would run mu1 linearly negative past
                        # exhaustion. Throttle site = min(mu0, mu1/a). Guard
                        # a>0: a<0 GROWS the chain (no exhaustion) and mu1/a<0
                        # would be a spurious negative site. (mirrored in
                        # get_reaction_rates' hijack block -- keep in sync)
                        mu_idx = self.polymer_pools[target_pool_idx].mu_indices
                        site = min(
                            max(0.0, y[mu_idx[0]]),
                            max(0.0, y[mu_idx[1]]) / float(self.reaction_eject_units[r_idx]),
                        ) / V_poly

                    rf *= site
                    # Direction-specific source availability (run-5 DASPK
                    # IDID=-7 forensics, adjudicated Part C): each direction's
                    # event flux must scale with the moments of the pool
                    # ACTUALLY BEING DEBITED in that direction. The REVERSE
                    # leg of a cross-pool exchange debits the dst pool
                    # (section-5 legs run (rf, src->dst), (rr, dst->src)), so
                    # its site factor comes from the dst pool's OWN moments --
                    # with the forward-reactant site a near-empty spawned pool
                    # was drained at a rate set by the healthy source pool
                    # (mu2 negative in ~1e-25 s, corrector divergence,
                    # h -> ~4e-21, IDID=-7). Same moment order as the forward
                    # site (mu1; mu0 for end-group rows), so the outflux -> 0
                    # CONTINUOUSLY (linearly) as the debited pool empties: no
                    # dt-dependent clamp, no activation threshold. Same-pool
                    # rows (incl. the throttled DISCRETE_CHIP / same-pool-VE
                    # shapes above, whose throttle multiplies BOTH directions
                    # by design) and rows with no product-side pool (dst ==
                    # -1: the reverse leg debits gas/explicit species only,
                    # whose availability is already in _C) keep the shared
                    # reactant-pool site. (mirrored in get_reaction_rates'
                    # hijack block -- keep in sync)
                    dst_pool_idx = self.reaction_dst_pool[r_idx]
                    if dst_pool_idx != -1 and dst_pool_idx != target_pool_idx:
                        moment_idx = self.polymer_pools[dst_pool_idx].mu_indices[1]
                        if self.is_end_group_reaction[r_idx]:
                            moment_idx = self.polymer_pools[dst_pool_idx].mu_indices[0]
                        rr *= max(0.0, y[moment_idx]) / V_poly
                    else:
                        rr *= site
            elif has_any_prod and not has_edge_prod:
                # PRODUCT-side proxy (review round 51): a reverse-stored
                # association row (e.g. edge daughter + edge daughter <=>
                # pool proxy, or gas radical + gas radical <=> pool proxy)
                # carries the proxy on the PRODUCT side, so rr's _C(proxy)
                # placeholder (1.0) must be replaced by the pool site
                # density, mirroring the reactant-side hijack above --
                # otherwise reverse homolysis discovery flux is unscaled.
                # Products are all core here (has_edge_prod is False), so
                # species_to_pool_indices indexing is safe. rf is untouched:
                # no reactant is a proxy on this branch (elif), so the
                # forward direction has no proxy concentration to replace.
                prod_pool_idx = -1
                for prod_p_idx in (p0, p1, p2):
                    if prod_p_idx != -1 and self.species_to_pool_indices[prod_p_idx] != -1:
                        prod_pool_idx = self.species_to_pool_indices[prod_p_idx]
                        break
                if prod_pool_idx != -1:
                    # Default to Mu1 (site density); end-group physics
                    # scales by Mu0 (chain density) -- same selection as
                    # the reactant-side block.
                    moment_idx = self.polymer_pools[prod_pool_idx].mu_indices[1]
                    if self.is_end_group_reaction[r_idx]:
                        moment_idx = self.polymer_pools[prod_pool_idx].mu_indices[0]
                    rr *= max(0.0, y[moment_idx]) / V_poly

            # Net rate (volumetric, in the phase volume chosen earlier)
            rate = rf - rr
            r_mol_s = rate * V_rxn
            abs_flux = abs(r_mol_s)

            if r_idx < n_rxn:
                self.core_reaction_rates[r_idx] = rate
            elif gated or self.reaction_refused[r_idx]:
                # counterfactual only -- the real entry stays 0 (the
                # consistency point enlargement reads). r71 FIX 3: refused
                # edge rows are routed here too -- their real
                # edge_reaction_rates entry must be zero (enlargement
                # input), the ungated entry stays a strictly-diagnostic
                # counterfactual (census only, never flux).
                self.edge_reaction_rates_ungated[r_idx - n_rxn] = rate
            else:
                self.edge_reaction_rates[r_idx - n_rxn] = rate
                self.edge_reaction_rates_ungated[r_idx - n_rxn] = rate

            core_rxn = r_idx < n_rxn

            if self.reaction_refused[r_idx]:
                # item 18 + r71 FIX 3: a refused reaction must contribute NO
                # flux ANYWHERE (stamp-but-keep): skip all reactant/product/
                # moment writes so no MW-weighted backbone mass is fabricated
                # in core -- AND no edge promotion flux is delivered. The
                # pre-r71 core-only gate left refused EDGE rows feeding
                # edge_species_rates / proxy_activity, which is exactly how
                # PP run-5's dead chain-scale radicals were promoted to core
                # (forensics: run5/RMG.log:910-955), where the adjudication
                # was lost and legacy mu1 flux collapsed DASPK. Refused rows
                # must not move species to core; the per-reaction ungated
                # counterfactual above is the only surviving record
                # (diagnostic, never flux).
                continue

            # 3. Apply Fluxes to Reactants (core reactions only -- edge
            #    reactions are diagnostic-only, matching simple.pyx).
            #    proxy_activity stays UNGATED on purpose: it feeds proxy
            #    similarity/spawn diagnostics, which want edge flux too --
            #    except gate-zeroed edge rows, suppressed just below (T9).
            #    Pool MOMENT flux is handled once per reaction in section 5.
            if r0 >= n_core:
                # Edge reactant: no core state, never a proxy; its
                # reactant-side edge flux is accumulated in section 4b
                # (diagnostic-only).
                pass
            elif self.is_pool_proxy[r0]:
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
                if r1 >= n_core:
                    pass  # edge reactant: see the r0 branch above
                elif self.is_pool_proxy[r1]:
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

            # 4b. Reactant-side edge accumulation (simple.pyx :509-515 sign
            #     parity): edge reactants get -= rate, so a reverse-stored
            #     decomposition row (rf = 0, net rate = -rr) registers
            #     POSITIVE daughter production for enlargement. Diagnostic-
            #     only: no dn_dt / core bookkeeping (edge species have no
            #     state), exactly like the product-side branch above.
            for p_slot in range(3):
                p_idx_tmp = ir[r_idx, p_slot]
                if p_idx_tmp == -1:
                    continue
                if p_idx_tmp >= n_core:
                    self.edge_species_rates_ungated[p_idx_tmp - n_core] -= rate
                    if not gated:
                        self.edge_species_rates[p_idx_tmp - n_core] -= rate

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
                elif arch == FLUX_VOLATILE_EJECTION:
                    # Volatile ejection (spec 2026-06-2x): a MIGRATION mirror
                    # A(pool) <=> volatile_gas + B(other pool). The moved chain
                    # loses `a` = reaction_eject_units backbone units to the
                    # discrete gas volatile as it lands in the destination
                    # pool, so the to-leg carries the SAME b0 chains but with
                    # E[n] shifted by sa (sa = -a forward: chain shrinks
                    # landing in dst; +a reverse: chain re-grows landing in
                    # src). from-leg loses the FULL bundle (identical to
                    # MIGRATION). mu2 shift is the exact quadratic
                    # (n + sa)^2 - n^2 = 2*sa*n + sa^2; sa*sa == a*a in both
                    # directions -> (b2 - 2a*b1 + a^2*b0) forward,
                    # (b2 + 2a*b1 + a^2*b0) reverse. The gas volatile itself is
                    # credited by the standard section-4 net-rate product path
                    # (NO gas moles added inside this branch).
                    src = self.reaction_src_pool[r_idx]
                    dst = self.reaction_dst_pool[r_idx]
                    a = float(self.reaction_eject_units[r_idx])
                    # -1 cannot reach here (init demotes unresolved pools);
                    # the checks are defensive only (mirror MIGRATION).
                    if src != -1 and dst != -1 and src != dst:
                        for ev_rate, from_pool, to_pool, sa in (
                                (rf, src, dst, -a), (rr, dst, src, +a)):
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
                            # from-leg loses the FULL bundle
                            dn_dt[f_idx[0]] -= ev_mol * b0
                            dn_dt[f_idx[1]] -= ev_mol * b1
                            # to-leg gains the a-SHIFTED bundle
                            dn_dt[t_idx[0]] += ev_mol * b0
                            dn_dt[t_idx[1]] += ev_mol * (b1 + sa * b0)
                            if mu2_ok:
                                dn_dt[f_idx[2]] -= ev_mol * b2
                                dn_dt[t_idx[2]] += ev_mol * (
                                    b2 + 2.0 * sa * b1 + a * a * b0)
                    elif src != -1 and src == dst:
                        # SAME-POOL volatile ejection (spec round-13): unzip /
                        # depropagation (a>0 sheds mass) or monomer/radical
                        # addition (a<0 gains mass) that lands back in the SAME
                        # pool. A dedicated single-pool chip-style write, NOT
                        # the cross-pool per-direction cancellation above: the
                        # cancellation can push mu2 the wrong way near
                        # exhaustion, so we write mu2 with the DISCRETE_CHIP
                        # clamp (>0) instead. `a` (signed) was resolved above.
                        # mu0 is untouched (no chain created/destroyed); the gas
                        # volatile flows through the standard section-4 path.
                        b0, b1, _b2, _mu2_ok = self._chain_bundle(
                            src, y, V_poly, self.is_end_group_reaction[r_idx])
                        if b0 != 0.0:
                            s_idx = self.polymer_pools[src].mu_indices
                            e_n = b1
                            if rf > 0.0:
                                # Forward: Delta n = -a. Signed a: for a<0 this
                                # ADDS mass (chain growth) -- correct. Clamp the
                                # mu2 decrement at >0 (identical to
                                # DISCRETE_CHIP): 2a*E[n]-a^2 < 0 is unphysical
                                # per-chain for a>0 but reachable in expectation
                                # near chip-size mean length; for a<0 it is
                                # always <0, so the mu2 growth term is dropped
                                # here (the exact +extension lives on reverse).
                                rf_mol = rf * V_rxn
                                dn_dt[s_idx[1]] -= rf_mol * a
                                mu2_dec = 2.0 * a * e_n - a * a
                                if mu2_dec > 0.0:
                                    dn_dt[s_idx[2]] -= rf_mol * mu2_dec
                            if rr > 0.0:
                                # Reverse: EXACT extension form (n+a)^2 - n^2 =
                                # +(2a*E[n] + a^2); unconditionally positive for
                                # a>0 (no clamp), mirrors DISCRETE_CHIP reverse.
                                rr_mol = rr * V_rxn
                                dn_dt[s_idx[1]] += rr_mol * a
                                dn_dt[s_idx[2]] += rr_mol * (
                                    2.0 * a * e_n + a * a)
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
        for pool_i, pool in enumerate(self.polymer_pools):
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
                        # Released monomer is emitted to the GAS species
                        # amount basis (incident 2026-07-03, design B-prime):
                        # small_src [mol/(m^3_poly s)] lands as
                        # dn_dt += r*V_poly [mol/s] on the gas-masked
                        # monomer_poly_index. Mass conservation: one gas
                        # monomer mole per drained mu1 repeat unit.
                        small_src[pool.monomer_poly_index] = r_events

            # Radical-homolysis initiation kernel (Stage 1, adjudicated round
            # 66). Independent of tail_kinetics (a custom tail closure does
            # not describe chain initiation); mutually exclusive with
            # k_scission > 0 and with the QSSA channel (validate_configuration
            # hard-errors on both shapes). Reads ONLY the flattened khom_*
            # arrays from _flatten_homolysis_state -- never the pool dict.
            #
            # Event: random backbone C-C homolysis at the RUNTIME temperature
            # (round 66: 'a scalar precomputed at 1100 K will poison any
            # ramp/TA replay'):
            #   k(T) = A * T^n * exp(-Ea/(R_gas*T))
            #   R    = k(T) * max(mu1 - mu0, 0)   [a chain of length n has
            #                                      n-1 breakable bonds]
            # Bond-weighted breaking-chain moments (E[n] = B1, E[n^2] = B2):
            #   B1 = (mu2 - mu1)/(mu1 - mu0);  B2 = (mu3 - mu2)/(mu1 - mu0)
            # with mu3 from the established log-Lagrange closure (computed
            # once per pool above). Parent debit:
            #   dmu0 -= R;  dmu1 -= R*B1;  dmu2 -= R*B2
            # EACH end-radical daughter pool (uniform cut at bond j of an
            # n-chain: E[j] = n/2, E[j^2] = n(2n-1)/6):
            #   dmu0 += R;  dmu1 += R*B1/2;  dmu2 += R*(2*B2 - B1)/6
            # STABLE FORMS (adjudicated round 67, P2): dividing by
            # (mu1 - mu0) and re-multiplying by R cancels catastrophically
            # near DP -> 1 exhaustion (mu1 ~= mu0), so every product-side
            # term is computed in its algebraically identical direct form:
            #   R*B1         = k*(mu2 - mu1)
            #   R*B2         = k*(mu3 - mu2)
            #   R*(2B2-B1)/6 = k*(2*mu3 - 3*mu2 + mu1)/6
            # (R = k*max(mu1-mu0, 0) itself stays the mu0 event rate).
            # Totals reproduce the legacy Ziff-McGrady random-scission moment
            # sources EXACTLY: dmu0_tot = +R = k*(mu1-mu0); dmu1_tot = 0
            # (mass conserved, machine precision); dmu2_tot = -R*(B1+B2)/3
            # = k*(mu1-mu3)/3. NO gas species term (homolysis releases no
            # volatiles) and NO reverse leg (one-way; round 66
            # §Reversibility -- recombination arrives via the discovered
            # chemistry conduit, not this kernel).
            if self.khom_enabled[pool_i]:
                B_hom = mu1 - mu0
                if B_hom > 0.0 and mu0 > 0.0:
                    T_hom = self.T.value_si
                    k_hom = (self.khom_A[pool_i]
                             * T_hom ** self.khom_n[pool_i]
                             * math.exp(-self.khom_Ea[pool_i]
                                        / (QSSA_R_GAS * T_hom)))
                    if not (0.0 < k_hom < QSSA_INF):
                        raise ValueError(
                            f"Pool {pool.label}: k_homolysis(T={T_hom:g} K) "
                            f"evaluated to a degenerate rate constant "
                            f"({k_hom!r}) -- A={self.khom_A[pool_i]:g}, "
                            f"n={self.khom_n[pool_i]:g}, "
                            f"Ea={self.khom_Ea[pool_i]:g} J/mol. Refusing "
                            f"to integrate a poisoned kernel.")
                    mu3_ok = np.isfinite(mu3)
                    # Direct (stable) moment differences: d1 = R*B1/k,
                    # d2 = R*B2/k, c2 = 6*(daughter mu2 credit)/k.
                    d1_hom = mu2 - mu1
                    d2_hom = (mu3 - mu2) if mu3_ok else -1.0
                    c2_hom = (2.0 * mu3 - 3.0 * mu2 + mu1) if mu3_ok else -1.0
                    if ((not mu3_ok) or d1_hom < 0.0 or d2_hom < 0.0
                            or c2_hom < 0.0):
                        # Loud guard (warn once per pool per rebuild,
                        # mirroring the realizability idiom above): the
                        # moment state is outside the kernel's domain
                        # (mu2 < mu1, a degenerate/overflowed mu3 closure,
                        # or -- round 67 P1 -- a nonrealizable combination
                        # where B1, B2 >= 0 but the daughter mu2 credit
                        # factor 2*B2 - B1 = c2/(mu1-mu0) is NEGATIVE, which
                        # would credit negative mu2 to both end-radical
                        # daughters). Contribute ZERO flux this step rather
                        # than integrate garbage fragment moments.
                        if pool.label not in self._khom_guard_warned:
                            self._khom_guard_warned.add(pool.label)
                            logging.warning(
                                "Polymer pool '%s': k_homolysis kernel "
                                "skipped -- nonfinite/negative moment "
                                "combination (mu0=%.6g, mu1=%.6g, mu2=%.6g, "
                                "mu3=%.6g; mu2-mu1=%.6g, mu3-mu2=%.6g, "
                                "daughter-mu2 factor 2*mu3-3*mu2+mu1=%.6g). "
                                "The kernel contributes zero flux until the "
                                "state re-enters its realizable domain.",
                                pool.label, mu0, mu1, mu2, mu3,
                                d1_hom, d2_hom, c2_hom)
                    else:
                        r_hom = k_hom * B_hom
                        dmu0_dt -= r_hom
                        dmu1_dt -= k_hom * d1_hom
                        dmu2_dt -= k_hom * d2_hom
                        # Fragment credits go straight to the daughter
                        # pools' moment slots (they are core species; the
                        # daughter's own accumulator flush is additive).
                        frag_mu0 = r_hom * V_poly
                        frag_mu1 = (k_hom * d1_hom / 2.0) * V_poly
                        frag_mu2 = (k_hom * c2_hom / 6.0) * V_poly
                        dn_dt[self.khom_prim_mu0[pool_i]] += frag_mu0
                        dn_dt[self.khom_prim_mu1[pool_i]] += frag_mu1
                        dn_dt[self.khom_prim_mu2[pool_i]] += frag_mu2
                        dn_dt[self.khom_sec_mu0[pool_i]] += frag_mu0
                        dn_dt[self.khom_sec_mu1[pool_i]] += frag_mu1
                        dn_dt[self.khom_sec_mu2[pool_i]] += frag_mu2

            # Side-group homolysis initiation kernel (FR1-K1, adjudicated
            # round 70). Reads ONLY the flattened sgh_* row arrays from
            # _flatten_side_group_state -- never the pool dicts. MAY coexist
            # with the k_homolysis kernel above (different bonds: side-group
            # C-X here vs backbone C-C there; no double-carry).
            #
            # Event, per channel with s = sites_per_unit and
            # k(T) = A*T^n*exp(-Ea/(R_gas*T)) per site, at the RUNTIME T:
            # every repeat unit carries s X sites, so a chain of length n
            # loses X at rate k*s*n and the reacting chain is picked with
            # probability ~ n (site-weighted). The chain transfers INTACT
            # (no chain cut) to the channel's X-loss feature pool:
            #   event rate      R = k*s*mu1
            #   parent debit    dmu_j -= k*s*mu_{j+1}   (j = 0, 1, 2;
            #                   equivalently R*E[n^j] with picked-chain
            #                   moments E[n^j] = mu_{j+1}/mu1 -- computed in
            #                   the STABLE direct product forms, never a
            #                   ratio re-multiplied by R)
            #   feature credit  EXACTLY the parent debit (total mu0/mu1/mu2
            #                   across the two pools conserved bitwise)
            #   gas X credit    +R of gas_product (one X radical per event)
            # mu3 from the established log-Lagrange closure (computed once
            # per pool above). MASS CONTRACT: the transferred chain lost one
            # X atom while the feature pool keeps the parent's monomer_mw,
            # so the feature pool's EXACT condensed mass is
            # mu1*MW - mu0*M_X (chain_mass_defect_g_mol = M_X, enforced at
            # flatten); on this kernel's contributions
            # d(condensed mass)/dt = -R*M_X = -d(gas X mass)/dt exactly.
            # v1 LIMITATION (adjudicated): the kernel acts on the PARENT
            # pool only -- feature pools carry no side_group_homolysis of
            # their own (no multi-loss cascade; they saturate as terminal
            # X-loss sinks).
            if self.sgh_enabled[pool_i] and mu1 > 0.0 and mu0 > 0.0:
                T_sgh = self.T.value_si
                mu3_ok_sgh = np.isfinite(mu3)
                for sgh_row in range(self.sgh_row_start[pool_i],
                                     self.sgh_row_end[pool_i]):
                    k_sgh = (self.sgh_A[sgh_row]
                             * T_sgh ** self.sgh_n[sgh_row]
                             * math.exp(-self.sgh_Ea[sgh_row]
                                        / (QSSA_R_GAS * T_sgh)))
                    if not (0.0 < k_sgh < QSSA_INF):
                        raise ValueError(
                            f"Pool {pool.label}: side_group_homolysis "
                            f"channel "
                            f"'{self.sgh_channel_labels[sgh_row]}' "
                            f"k(T={T_sgh:g} K) evaluated to a degenerate "
                            f"rate constant ({k_sgh!r}) -- "
                            f"A={self.sgh_A[sgh_row]:g}, "
                            f"n={self.sgh_n[sgh_row]:g}, "
                            f"Ea={self.sgh_Ea[sgh_row]:g} J/mol. Refusing "
                            f"to integrate a poisoned kernel.")
                    if (not mu3_ok_sgh) or mu2 < mu1 or mu3 < mu2:
                        # Loud out-of-domain guard (warn once per
                        # pool/channel per rebuild, sibling-kernel
                        # convention): the debit consumes mu2 and the mu3
                        # closure directly, so a nonrealizable shape
                        # (mu2 < mu1 -- impossible for chains with n >= 1;
                        # nonfinite mu3; mu3 < mu2) would transfer negative
                        # or garbage moment bundles. Contribute ZERO flux
                        # this step. (mu1 <= 0 / mu0 <= 0 is the SILENT
                        # in-domain zero above: no sites, no events.)
                        key = (pool.label,
                               self.sgh_channel_labels[sgh_row])
                        if key not in self._sgh_guard_warned:
                            self._sgh_guard_warned.add(key)
                            logging.warning(
                                "Polymer pool '%s': side_group_homolysis "
                                "channel '%s' skipped -- "
                                "nonfinite/negative moment combination "
                                "(mu0=%.6g, mu1=%.6g, mu2=%.6g, mu3=%.6g; "
                                "require mu2 >= mu1 and finite "
                                "mu3 >= mu2). The kernel contributes zero "
                                "flux until the state re-enters its "
                                "realizable domain.",
                                pool.label,
                                self.sgh_channel_labels[sgh_row],
                                mu0, mu1, mu2, mu3)
                        continue
                    ks_sgh = k_sgh * self.sgh_sites[sgh_row]
                    r0_sgh = ks_sgh * mu1     # R: events == chains == gas X
                    r1_sgh = ks_sgh * mu2     # R*E[n], stable direct form
                    r2_sgh = ks_sgh * mu3     # R*E[n^2], stable direct form
                    dmu0_dt -= r0_sgh
                    dmu1_dt -= r1_sgh
                    dmu2_dt -= r2_sgh
                    # Feature-pool credits go straight to the destination
                    # moment slots -- the SAME r*V_poly products as the
                    # parent debit flush below, so the two-pool totals
                    # cancel bitwise.
                    dn_dt[self.sgh_dst_mu0[sgh_row]] += r0_sgh * V_poly
                    dn_dt[self.sgh_dst_mu1[sgh_row]] += r1_sgh * V_poly
                    dn_dt[self.sgh_dst_mu2[sgh_row]] += r2_sgh * V_poly
                    # Gas X radical: +R through the same small_src ->
                    # dn_dt * V_poly path as the unzip/QSSA releases
                    # (additive: channels may share one gas species).
                    g_idx_sgh = int(self.sgh_gas[sgh_row])
                    small_src[g_idx_sgh] = (
                        small_src.get(g_idx_sgh, 0.0) + r0_sgh)

            # radical_qssa_unzip channel (M2 rate law). Reads ONLY the
            # flattened solver-owned qssa_* arrays from M1 -- NEVER the pool
            # dict -- and gates on qssa_enabled, never on A != 0. Mutually
            # exclusive with k_unzip > 0 (M1 hard error), independent of
            # tail_kinetics (a custom tail closure does not describe
            # radical depropagation).
            #
            # Rate law (isothermal QSSA on the active chain-END radical R):
            #   B    = max(mu1 - mu0, 0)         breakable backbone bonds
            #                                    (mol/m^3 of condensed volume:
            #                                    a chain of length k has k-1
            #                                    bonds, so sum = mu1 - mu0)
            #   G_R  = 2 f ki(T) B               radical generation (homolysis
            #                                    makes TWO end radicals/bond)
            #   no transfer:  0 = G_R - 2 kt R^2          -> R = sqrt(f ki B / kt)
            #   transfer:     0 = G_R - ktr R - 2 kt R^2  -> R = (sqrt(ktr^2
            #                     + 8 kt G_R) - ktr) / (4 kt)
            #     (transfer = first-order active-END loss: H-abstraction turns
            #     an unzipping end radical into a mid-chain radical)
            #   r_mono = y_m kdp(T) R            monomer release, mol/(m^3 s)
            # Arrhenius k(T) = A T^n exp(-Ea/(R_gas T)), Ea in J/mol (M1 SI
            # pin). Units: ki/kdp in s^-1; kt in m^3/(mol s) (bimolecular);
            # ktr in s^-1 -- PSEUDO-FIRST-ORDER, because ktr multiplies R
            # DIRECTLY in the balance above (a literature bimolecular k_tr
            # [L/(mol s) or m^3/(mol s)] must be premultiplied by the
            # relevant substrate concentration [mol/m^3, SI] BEFORE entering
            # the config). The mol/m^3 concentration basis established at
            # C_poly = y/V_poly above means NO unit conversion is applied
            # here.
            r_qssa = 0.0
            if self.qssa_enabled[pool_i]:
                B_qssa = max(mu1 - mu0, 0.0)
                if self.qssa_weaklink[pool_i]:
                    # Weak-link allyl/U-state channel (milestones ii+iii;
                    # rate law per review r36). Extensions over the legacy
                    # algebra (everything else identical -- the bridge pin
                    # proves it):
                    #   kt_total = kt_rec(T) + kt_disp(T)   replaces kt
                    #   u_active = min(max(U,0)/V_poly, B)  [mol/m^3]
                    #     ACTIVE-SITE CLAMP (r36 P1-2): the weak channel
                    #     fissions a backbone bond ALLYLIC to the
                    #     unsaturated end; ends beyond the remaining
                    #     backbone-bond count have nothing to fission. At
                    #     B = 0 the whole channel is INERT (gate below) --
                    #     no monomer drain that could push mu1 < mu0 past
                    #     the DP->1 self-termination floor, no U sink, no
                    #     production.
                    #   G_R = 2 f ki(T) B + 1 f ki_allyl(T) u_active
                    #     (nu = 1: ONE unzipping radical per weak-link
                    #     fission; the allylic co-fragment does not unzip)
                    #   dU/dt = + kt_disp R_ss^2 * max(0, 1 - U/(2*mu0))
                    #             * V_poly
                    #           - f ki_allyl(T) u_active V_poly
                    #     Production: halved radical-disappearance
                    #     convention (kt_disp*R^2 IS the disproportionation
                    #     EVENT rate, 1 unsaturated end per event) -- ONLY
                    #     disproportionation sources U; recombination and
                    #     random initiation do not. NO f on production
                    #     (r36): R_ss already carries the escape
                    #     efficiency; termination is not a caged event.
                    #     CAPACITY THROTTLE (r36 P1-3): the linear factor
                    #     max(0, 1 - U/(2*mu0_amount)) is EXACTLY zero at
                    #     the chain-end capacity 2*mu0 (each chain has at
                    #     most 2 tail ends), so U(t) cannot integrate past
                    #     it -- the t=0 census guards the IC, this guards
                    #     the trajectory. mu0 here is the TAIL mu0 amount
                    #     [mol], the same distribution B is counted on; at
                    #     mu0 = 0 the throttle is 0 (no ends, no
                    #     capacity).
                    #     Sink: f-SYMMETRIC with the G_R term (r36 P1-1) --
                    #     f is the escaped-radical-pair efficiency, a
                    #     caged recombination restores the allylic bond,
                    #     so the caged fraction consumes no U.
                    #   Volume convention matches the moment ODEs below
                    #   exactly: concentration-rate * V_poly -> amount-rate
                    #   (dn_dt[idx_mu1] += dmu1_dt * V_poly).
                    u_amount = max(0.0, y[self.qssa_u_slot[pool_i]])
                    if B_qssa > 0.0:
                        u_active = u_amount / V_poly
                        if u_active > B_qssa:
                            u_active = B_qssa
                        T_qssa = self.T.value_si
                        RT_qssa = QSSA_R_GAS * T_qssa
                        ki_qssa = (self.qssa_ki_A[pool_i]
                                   * T_qssa ** self.qssa_ki_n[pool_i]
                                   * math.exp(-self.qssa_ki_Ea[pool_i] / RT_qssa))
                        kdp_qssa = (self.qssa_kdp_A[pool_i]
                                    * T_qssa ** self.qssa_kdp_n[pool_i]
                                    * math.exp(-self.qssa_kdp_Ea[pool_i] / RT_qssa))
                        kia_qssa = (self.qssa_kia_A[pool_i]
                                    * T_qssa ** self.qssa_kia_n[pool_i]
                                    * math.exp(-self.qssa_kia_Ea[pool_i] / RT_qssa))
                        ktrec_qssa = (self.qssa_ktrec_A[pool_i]
                                      * T_qssa ** self.qssa_ktrec_n[pool_i]
                                      * math.exp(-self.qssa_ktrec_Ea[pool_i] / RT_qssa))
                        ktdisp_qssa = (self.qssa_ktdisp_A[pool_i]
                                       * T_qssa ** self.qssa_ktdisp_n[pool_i]
                                       * math.exp(-self.qssa_ktdisp_Ea[pool_i] / RT_qssa))
                        kt_qssa = ktrec_qssa + ktdisp_qssa
                        # Same degenerate-rate posture as the legacy branch:
                        # kt_total underflowing to 0 divides R_ss by zero;
                        # any inf/NaN evaluation poisons the residual.
                        if not (ki_qssa < QSSA_INF and kdp_qssa < QSSA_INF
                                and kia_qssa < QSSA_INF
                                and ktrec_qssa < QSSA_INF
                                and ktdisp_qssa < QSSA_INF
                                and 0.0 < kt_qssa < QSSA_INF):
                            _raise_degenerate_qssa_rates(
                                pool.label, T_qssa, ki_qssa, kdp_qssa,
                                kt_qssa, kia=kia_qssa, ktrec=ktrec_qssa,
                                ktdisp=ktdisp_qssa)
                        fkiB = self.qssa_efficiency[pool_i] * ki_qssa * B_qssa
                        G_R_qssa = (2.0 * fkiB
                                    + self.qssa_efficiency[pool_i]
                                    * kia_qssa * u_active)
                        if self.qssa_has_transfer[pool_i]:
                            ktr_qssa = (self.qssa_ktr_A[pool_i]
                                        * T_qssa ** self.qssa_ktr_n[pool_i]
                                        * math.exp(-self.qssa_ktr_Ea[pool_i] / RT_qssa))
                            if not ktr_qssa < QSSA_INF:
                                _raise_degenerate_qssa_rates(
                                    pool.label, T_qssa, ki_qssa, kdp_qssa,
                                    kt_qssa, ktr_qssa, kia=kia_qssa,
                                    ktrec=ktrec_qssa, ktdisp=ktdisp_qssa)
                            R_ss = ((math.sqrt(ktr_qssa * ktr_qssa
                                               + 8.0 * kt_qssa * G_R_qssa)
                                     - ktr_qssa) / (4.0 * kt_qssa))
                        else:
                            R_ss = math.sqrt(G_R_qssa / (2.0 * kt_qssa))
                        r_qssa = (self.qssa_monomer_yield[pool_i]
                                  * kdp_qssa * R_ss)
                        # dU/dt (amount basis, mol/s): capacity-throttled
                        # production by the disproportionation branch ONLY
                        # (no f -- R_ss carries the escape efficiency);
                        # f-symmetric sink by allyl fission of the ACTIVE
                        # unsaturated ends. See the law block above.
                        cap_qssa = 2.0 * mu0_mol  # tail chain-end capacity
                        if cap_qssa > 0.0:
                            throttle_qssa = 1.0 - u_amount / cap_qssa
                            if throttle_qssa < 0.0:
                                throttle_qssa = 0.0
                        else:
                            throttle_qssa = 0.0
                        du_dt[self.qssa_u_slot[pool_i] - n_core] += (
                            ktdisp_qssa * R_ss * R_ss * throttle_qssa
                            - self.qssa_efficiency[pool_i]
                            * kia_qssa * u_active) * V_poly
                elif B_qssa > 0.0:
                    T_qssa = self.T.value_si
                    RT_qssa = QSSA_R_GAS * T_qssa
                    ki_qssa = (self.qssa_ki_A[pool_i]
                               * T_qssa ** self.qssa_ki_n[pool_i]
                               * math.exp(-self.qssa_ki_Ea[pool_i] / RT_qssa))
                    kdp_qssa = (self.qssa_kdp_A[pool_i]
                                * T_qssa ** self.qssa_kdp_n[pool_i]
                                * math.exp(-self.qssa_kdp_Ea[pool_i] / RT_qssa))
                    kt_qssa = (self.qssa_kt_A[pool_i]
                               * T_qssa ** self.qssa_kt_n[pool_i]
                               * math.exp(-self.qssa_kt_Ea[pool_i] / RT_qssa))
                    # Runtime degenerate-rate guard (fail-loud, cheap hot
                    # path). M1 pins config finiteness, but EVALUATION can
                    # still degenerate at solver T: kt underflows to 0.0
                    # (sqrt(f*ki*B/kt) divides by zero -> inf) or a large
                    # A/n overflows to inf/NaN. `x < QSSA_INF` rejects inf
                    # AND NaN in one compare; ki/kdp underflow to 0.0 is
                    # fine (the rate smoothly reaches 0), kt is not.
                    if not (ki_qssa < QSSA_INF and kdp_qssa < QSSA_INF
                            and 0.0 < kt_qssa < QSSA_INF):
                        _raise_degenerate_qssa_rates(
                            pool.label, T_qssa, ki_qssa, kdp_qssa, kt_qssa)
                    fkiB = self.qssa_efficiency[pool_i] * ki_qssa * B_qssa
                    G_R_qssa = 2.0 * fkiB  # two end radicals per homolysis
                    if self.qssa_has_transfer[pool_i]:
                        ktr_qssa = (self.qssa_ktr_A[pool_i]
                                    * T_qssa ** self.qssa_ktr_n[pool_i]
                                    * math.exp(-self.qssa_ktr_Ea[pool_i] / RT_qssa))
                        if not ktr_qssa < QSSA_INF:
                            _raise_degenerate_qssa_rates(
                                pool.label, T_qssa, ki_qssa, kdp_qssa,
                                kt_qssa, ktr_qssa)
                        R_ss = ((math.sqrt(ktr_qssa * ktr_qssa
                                           + 8.0 * kt_qssa * G_R_qssa)
                                 - ktr_qssa) / (4.0 * kt_qssa))
                    else:
                        # == sqrt(G_R / (2 kt)), simplified
                        R_ss = math.sqrt(fkiB / kt_qssa)
                    r_qssa = self.qssa_monomer_yield[pool_i] * kdp_qssa * R_ss
                if r_qssa > 0.0:
                    # Chain-END monomer release signature: mu0 untouched (no
                    # chain created/destroyed), mu1 drains one unit per
                    # release, mu2 drains (2 E[n] - 1) per release with the
                    # same-pool-VE clamp (>0 only: the drain must never make
                    # mu2 increase; mu0 ~ 0 guarded by the eps clamp).
                    # monomer_yield already scales r_qssa, so the moment
                    # drain and the gas emission below scale TOGETHER --
                    # scaling only one side would fabricate/destroy mass.
                    dmu1_dt -= r_qssa
                    qssa_mu2_dec = 2.0 * (mu1 / max(mu0, SMALL_EPS)) - 1.0
                    if qssa_mu2_dec > 0.0:
                        dmu2_dt -= r_qssa * qssa_mu2_dec
                    # monomer_poly_index is non-None whenever enabled (M1
                    # invariant); emission flows through the SAME small_src
                    # -> dn_dt * V_poly path as the k_unzip channel, i.e.
                    # to the GAS species amount basis (incident 2026-07-03,
                    # design B-prime).
                    small_src[pool.monomer_poly_index] = (
                        small_src.get(pool.monomer_poly_index, 0.0) + r_qssa)

            # End-radical DEPROPAGATION kernel (adjudicated round 74 SS2, the
            # run-6 no-outlet wall fix). Reads ONLY the flattened kdep_*
            # arrays from _flatten_depropagation_state -- never the pool
            # dict. Gates on kdep_enabled; mutually exclusive with
            # k_unzip > 0 / radical_qssa_unzip / k_homolysis on the same
            # pool (validate_configuration hard-errors), so at most one
            # chain-end release arm fires per pool.
            #
            # Law, per radical-end pool (ONE active radical end per chain),
            # k_dep(T) = A*T^n*exp(-Ea/(R_gas*T)) at the RUNTIME T:
            #   R    = k_dep * mu0 * g   unzip events == monomer release
            #   gas  = +R at kdep_gas    (the SAME float as the mu1 drain:
            #                             d(condensed) + d(gas monomer) = 0
            #                             EXACTLY under MW multiplication)
            #   dmu1 = -R
            #   dmu2 = -k_dep*(2*mu1 - mu0), computed as
            #          -k_dep*mu0*(g + 2*sp(mean - 1)): a DP=n chain loses
            #          (2n - 1) from mu2 per event; the smooth positive-part
            #          sp keeps the excess term nonnegative without a
            #          max(...,0) cliff (identical to the nominal law up to
            #          O(W^2/(mean-1)) in the realizable region).
            #   dmu0 = -k_dep * N1, N1 = mu0 * p1 from the EXISTING
            #          discrete/gamma closure + smooth terminal boundary
            #          floor (_deprop_dp1_fraction): DP=1 chains release
            #          their last unit and TERMINATE -- never a permanent
            #          dmu0 = 0 (r74: that stalls at a one-repeat-per-chain
            #          residue).
            # SMOOTH exhaustion gate (r74 SS5, designed in from the start):
            # g = 1 - sp(1 - mean) is EXACTLY 1 for mean >= 1 and rolls off
            # C2-smoothly below; dmu0 is deliberately NOT gated by g, so a
            # noise state mean < 1 drains chains FASTER than units
            # (self-healing back toward the realizable cone, no stiff
            # no-outlet grind).
            if self.kdep_enabled[pool_i] and mu0 > 0.0:
                T_kdep = self.T.value_si
                k_dep = (self.kdep_A[pool_i]
                         * T_kdep ** self.kdep_n[pool_i]
                         * math.exp(-self.kdep_Ea[pool_i]
                                    / (QSSA_R_GAS * T_kdep)))
                if not (0.0 < k_dep < QSSA_INF):
                    raise ValueError(
                        f"Pool {pool.label}: k_depropagation(T={T_kdep:g} K) "
                        f"evaluated to a degenerate rate constant "
                        f"({k_dep!r}) -- A={self.kdep_A[pool_i]:g}, "
                        f"n={self.kdep_n[pool_i]:g}, "
                        f"Ea={self.kdep_Ea[pool_i]:g} J/mol. Refusing to "
                        f"integrate a poisoned kernel.")
                mean_kdep = mu1 / mu0
                g_kdep = 1.0 - _smooth_pos(1.0 - mean_kdep, KDEP_GATE_WIDTH)
                r_kdep = k_dep * mu0 * g_kdep
                if r_kdep > 0.0:
                    dmu1_dt -= r_kdep
                    dmu2_dt -= k_dep * mu0 * (
                        g_kdep + 2.0 * _smooth_pos(mean_kdep - 1.0,
                                                   KDEP_GATE_WIDTH))
                    g_idx_kdep = int(self.kdep_gas[pool_i])
                    small_src[g_idx_kdep] = (
                        small_src.get(g_idx_kdep, 0.0) + r_kdep)
                dmu0_dt -= k_dep * mu0 * _deprop_dp1_fraction(mu0, mu1, mu2)

            # Hybrid Handshake
            tail_mean = mu1 / mu0 if mu0 > SMALL_EPS else 0.0
            valid_tail = (mu0 > TAIL_CONC_MIN) and (tail_mean > xs + 1e-9)

            # Per-chain unzip frequency feeding the handshake: the legacy
            # k_unzip IS that frequency; the QSSA equivalent is
            # k = r_mono / mu0 (release events per chain per second), so
            # QSSA pools do not strand low-DP condensed residue. The two
            # channels are mutually exclusive (M1), so at most one arm fires.
            k_chain_handshake = 0.0
            if pool.k_unzip > 0:
                k_chain_handshake = pool.k_unzip
            elif self.qssa_enabled[pool_i] and r_qssa > 0.0:
                k_chain_handshake = r_qssa / max(mu0, SMALL_EPS)

            if valid_tail and k_chain_handshake > 0.0:
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

                    F = k_chain_handshake * N_boundary

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

        # --- DEBUG: report first non-finite / blow-up in residual ---
        if not getattr(self, '_dbg_reported', False):
            if (not np.isfinite(V_gas)) or (not np.all(np.isfinite(dn_dt))):
                self._dbg_reported = True
                bad = np.where(~np.isfinite(dn_dt))[0]
                print(f"\n[POLY-DBG] NON-FINITE in residual at t={t:.3e}")
                print(f"[POLY-DBG]   V_gas={V_gas:.6e} V_poly={V_poly:.6e} "
                      f"n_gas={float(np.sum(y[:n_core][mask])):.6e}")
                cg = C_gas[mask]
                print(f"[POLY-DBG]   max|C_gas|={float(np.max(np.abs(cg))) if cg.size else 0.0:.4e} "
                      f"max|dn_dt finite|={float(np.max(np.abs(dn_dt[np.isfinite(dn_dt)]))) if np.any(np.isfinite(dn_dt)) else 0.0:.4e}")
                print(f"[POLY-DBG]   non-finite dn_dt indices (first 10): {list(bad[:10])}")
                for bi in bad[:10]:
                    print(f"[POLY-DBG]     idx {bi}: y={y[bi]:.4e} gas={bool(mask[bi])}")

        if self.num_qssa_u > 0:
            # Weak-link layout: y = [core species | U slots]; the residual
            # vector follows the same layout. The legacy (num_qssa_u == 0)
            # return below is bit-for-bit the pre-milestone expression.
            return (np.concatenate((dn_dt, du_dt)) - dydt), 1
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

            # Determine Phase for Volume Basis. Edge participants (review
            # round 51: edge-reactant rows are evaluated, not skipped) are
            # judged by prospective_gas_mask -- the only mask covering edge
            # slots; its core prefix is bit-identical to gas_species_mask by
            # rider R1. Raw indexing of the core-sized arrays with an edge
            # slot would be an IndexError.
            is_poly = False
            if r0 != -1 and not (self.gas_species_mask[r0] if r0 < n_core
                                 else self.prospective_gas_mask[r0]):
                is_poly = True
            elif r1 != -1 and not (self.gas_species_mask[r1] if r1 < n_core
                                   else self.prospective_gas_mask[r1]):
                is_poly = True
            elif r2 != -1 and not (self.gas_species_mask[r2] if r2 < n_core
                                   else self.prospective_gas_mask[r2]):
                is_poly = True

            V_rxn = V_poly if is_poly else V_gas

            kf = self.kf[r_idx]
            kb = self.kb[r_idx]

            # Forward Rate. An EDGE reactant has no state in y -> zero
            # concentration, so rf vanishes (simple.pyx :443-450 parity).
            rf = kf
            if r0 >= n_core:
                rf *= 0.0
            elif self.species_to_pool_indices[r0] != -1:
                rf *= 1.0
            elif self.gas_species_mask[r0]:
                rf *= C_gas[r0]
            else:
                rf *= C_poly[r0]

            if r1 != -1:
                if r1 >= n_core:
                    rf *= 0.0
                elif self.species_to_pool_indices[r1] != -1:
                    rf *= 1.0
                elif self.gas_species_mask[r1]:
                    rf *= C_gas[r1]
                else:
                    rf *= C_poly[r1]

            if r2 != -1:
                if r2 >= n_core:
                    rf *= 0.0
                elif self.species_to_pool_indices[r2] != -1:
                    rf *= 1.0
                elif self.gas_species_mask[r2]:
                    rf *= C_gas[r2]
                else:
                    rf *= C_poly[r2]

            # Reverse Rate. An EDGE product has no state in y, so rr is
            # UNCOMPUTABLE and stays 0 (simple.pyx :452-453 parity -- the
            # concentration-availability hole, same as the residual).
            p0 = ip[r_idx, 0]
            p1 = ip[r_idx, 1]
            p2 = ip[r_idx, 2]
            rr = 0.0
            if kb > 0 and p0 < n_core and p1 < n_core and p2 < n_core:
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

            # [THE HIJACK] Scale by Moment if Reactant is Proxy (edge
            # reactants are never proxies: pool proxies are core species)
            p0_pool_idx = self.species_to_pool_indices[r0] if (r0 != -1 and r0 < n_core) else -1
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
                elif (self.reaction_flux_archetype[r_idx] == FLUX_VOLATILE_EJECTION
                        and self.is_end_group_reaction[r_idx]
                        and self.reaction_src_pool[r_idx] != -1
                        and self.reaction_src_pool[r_idx] == self.reaction_dst_pool[r_idx]
                        and self.reaction_eject_units[r_idx] > 0.0):
                    # Same-pool a>0 VE exhaustion throttle -- parity with the
                    # residual's section-2 site scaling (keep in sync). site =
                    # min(mu0, mu1/a); guard a>0 (a<0 grows, no throttle).
                    mu1_idx = self.polymer_pools[p0_pool_idx].mu_indices[1]
                    if self.pool_mu1_indices[p0_pool_idx] != -1:
                        mu1_idx = self.pool_mu1_indices[p0_pool_idx]
                    site = min(
                        max(0.0, y[moment_idx]),
                        max(0.0, y[mu1_idx]) / float(self.reaction_eject_units[r_idx]),
                    ) / V_poly

                rf *= site
                # Direction-specific source availability -- keep in sync with
                # the residual's section-2 scaling (adjudicated Part C): the
                # reverse leg of a cross-pool exchange debits the dst pool,
                # so rr's site comes from the dst pool's own moments (mu1;
                # mu0 for end-group rows; same pool_mu1_indices override as
                # the forward site). Same-pool and no-dst rows keep the
                # shared (possibly throttled) reactant-pool site.
                dst_pool_idx = self.reaction_dst_pool[r_idx]
                if dst_pool_idx != -1 and dst_pool_idx != p0_pool_idx:
                    moment_idx = self.polymer_pools[dst_pool_idx].mu_indices[1]
                    if self.pool_mu1_indices[dst_pool_idx] != -1:
                        moment_idx = self.pool_mu1_indices[dst_pool_idx]
                    if self.is_end_group_reaction[r_idx]:
                        moment_idx = self.pool_mu0_indices[dst_pool_idx]
                    rr *= max(0.0, y[moment_idx]) / V_poly
                else:
                    rr *= site
            elif rr != 0.0:
                # PRODUCT-side proxy mirror (review round 51; keep in sync
                # with the residual's product-side scaling): a reverse-stored
                # association row carries the proxy on the PRODUCT side, so
                # rr's 1.0 placeholder is replaced by the pool site density.
                # rr != 0.0 implies all products are core (the rr hole above),
                # so species_to_pool_indices indexing is safe. rf untouched:
                # no reactant is a proxy on this branch (elif).
                prod_pool_idx = -1
                if p0 != -1 and self.species_to_pool_indices[p0] != -1:
                    prod_pool_idx = self.species_to_pool_indices[p0]
                elif p1 != -1 and self.species_to_pool_indices[p1] != -1:
                    prod_pool_idx = self.species_to_pool_indices[p1]
                elif p2 != -1 and self.species_to_pool_indices[p2] != -1:
                    prod_pool_idx = self.species_to_pool_indices[p2]
                if prod_pool_idx != -1:
                    moment_idx = self.polymer_pools[prod_pool_idx].mu_indices[1]
                    if self.pool_mu1_indices[prod_pool_idx] != -1:
                        moment_idx = self.pool_mu1_indices[prod_pool_idx]
                    if self.is_end_group_reaction[r_idx]:
                        moment_idx = self.pool_mu0_indices[prod_pool_idx]
                    rr *= max(0.0, y[moment_idx]) / V_poly

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
