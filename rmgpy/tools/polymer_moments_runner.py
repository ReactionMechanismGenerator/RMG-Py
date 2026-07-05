#!/usr/bin/env python3
"""CLI reference runner for the polymer moments artifact (design spec §8).

Drives HybridPolymerSystem — the ORACLE — from consumer-world inputs only:
  polymer_pools.json (schema 2.0) + chem.yaml + piecewise-isothermal
  T-profile + time grid + mass-transfer spec  ->  moment-trajectory CSV.

It is the cross-validation oracle for the numpy/Cantera consumer (TA), not a
reimplementation. Normative contract: docs/polymer_moments_format.md.

Restart pattern per temperature segment (proven on the analytic two-segment
check; see test/rmgpy/tools/polymerMomentsRunnerTest.py):
  rs.T = Quantity((T_k, 'K')); rs.generate_rate_coefficients(core_rxns, []);
  rs.t0 = t_start; rs.y0 = y_carry; rs.set_initial_derivative();
  rs.initialize_solver(); then rs.advance(t) per grid point.

Caveats on oracle equivalence
------------------------------
Oracle equivalence to the *generating RMG run* is:

* **Exact** for the proxy/pool subset: every artifact-listed reaction is
  restamped and its ``kinetics.reversible`` restores the generating solver's
  reversibility flag, so these reactions run identically to the original.

* **Exact** for the time-grid machinery, moment channels, and mass-transfer.

* **Exact for ordinary gas chemistry (post-fix chem.yaml):** chem.yaml now
  records reversibility in the equation arrow (``rmgpy/cantera.py
  get_reaction_equation`` emits ``=>`` for irreversible reactions, ``<=>``
  otherwise), so reactions with NO artifact entry (i.e. not proxy-touching)
  round-trip their reversibility faithfully from the arrow in every consumer —
  Cantera, the reference runner, and the TA numpy consumer alike.  The
  artifact's ``kinetics.reversible`` restoration for listed entries now
  agrees with the arrow and is kept as belt-and-braces.

  **Legacy caveat:** chem.yaml files written by RMG *before* this exporter fix
  carry an all-``<=>`` arrow and no separate reversibility marker, so an
  originally-irreversible ordinary (non-artifact-listed) reaction in such a
  file runs *reversible* in every consumer.  A warning is emitted at runtime
  when a non-artifact-listed reaction is present so old files are flagged.
"""

import argparse
import contextlib
import csv
import io
import json
import logging
import math
import sys

import numpy as np
import yaml

from rmgpy.kinetics import Arrhenius
from rmgpy.molecule.element import get_element
from rmgpy.quantity import Quantity
from rmgpy.reaction import Reaction
from rmgpy.solver.polymer import (
    HybridPolymerSystem,
    MassTransferConfig,
    PolymerPoolConfig,
    validate_radical_qssa_unzip,
)
from rmgpy.species import Species
from rmgpy.thermo import NASA, NASAPolynomial

ARCHETYPE_INTS = {
    "same_pool/1": 1,
    "migration/1": 2,
    "scission_fragment/1": 3,
    "legacy_mu1/1": 4,
    "discrete_chip/1": 5,
    "volatile_ejection/1": 6,
}


def _validated_eject_units(e):
    """Validate + return a volatile_ejection/1 row's SIGNED eject_units.

    The emitter writes exactly ONE VE params sub-shape:
    ``params = {"eject_units": float}`` (rmgpy/polymer.py
    compile_polymer_reaction_entries; a = net non-polymer mass /
    source-monomer MW, SIGNED -- a < 0 is net mass gain from the gas
    co-reactants, see the signed-VE contract in rmgpy/polymer.py). Reject
    anything else with an actionable error, never KeyError and never a
    silent 0.0 default: defaulting would launder the atom-transfer debit
    away (the moved chain lands un-shrunk while the gas volatile still
    appears -- fabricated mass), and adapting chip-style vocabulary would
    guess semantics."""
    eid = e.get("id")
    params = e.get("params")
    if not isinstance(params, dict) or set(params) != {"eject_units"}:
        raise ValueError(
            f"reactions[] entry {eid!r} (volatile_ejection/1) must carry "
            f"params = {{'eject_units': <signed float>}} exactly -- the "
            f"only VE params shape the emitter writes -- got {params!r}. "
            f"Fix the artifact; defaulting would silently zero the "
            f"atom-transfer debit.")
    a = params["eject_units"]
    if isinstance(a, bool) or not isinstance(a, (int, float)) \
            or not math.isfinite(float(a)):
        raise ValueError(
            f"reactions[] entry {eid!r} (volatile_ejection/1) has "
            f"eject_units={a!r}; it must be a finite SIGNED number (the "
            f"source-monomer-equivalents transferred to the gas "
            f"co-participants, rmgpy/polymer.py "
            f"compute_volatile_ejection_units). Fix the artifact.")
    return float(a)


def _species_from_yaml(entry):
    name = entry["name"]
    th = entry.get("thermo", {})
    thermo = None
    if th.get("model") == "NASA7":
        tranges = th["temperature-ranges"]
        rows = th["data"]
        polys = []
        for i, coeffs in enumerate(rows):
            polys.append(NASAPolynomial(coeffs=[float(c) for c in coeffs],
                                        Tmin=(float(tranges[i]), "K"),
                                        Tmax=(float(tranges[i + 1]), "K")))
        thermo = NASA(polynomials=polys,
                      Tmin=(float(tranges[0]), "K"),
                      Tmax=(float(tranges[-1]), "K"))
    spc = Species(label=name, thermo=thermo)
    spc.molecule = []  # consumer-world species: label + thermo only
    # Structure does not cross the artifact boundary, so MW must come from
    # artifact fields too: the chem.yaml elemental composition. The solver's
    # thermo reference-state tripwire reads this species-level quantity.
    mw_kg_mol = 0.0
    for symbol, count in (entry.get("composition") or {}).items():
        if symbol == "E":
            continue  # electron count: negligible mass; not a get_element symbol
        mw_kg_mol += get_element(symbol).mass * float(count)  # .mass is kg/mol
    if mw_kg_mol > 0.0:
        spc.molecular_weight = (mw_kg_mol, "kg/mol")
    return spc


def _parse_equation(eq):
    if "(+ M)" in eq or "(+M)" in eq or " + M " in eq:
        raise NotImplementedError(f"third-body reactions unsupported in v1: {eq}")
    if "<=>" in eq:
        lhs, rhs = eq.split("<=>")
        reversible = True
    elif "=>" in eq:
        lhs, rhs = eq.split("=>")
        reversible = False
    else:
        raise ValueError(f"cannot parse equation: {eq}")
    reactants = [tok.strip() for tok in lhs.split(" + ")]
    products = [tok.strip() for tok in rhs.split(" + ")]
    return [r for r in reactants if r], [p for p in products if p], reversible


_A_UNITS_BY_ORDER = {1: "s^-1", 2: "m^3/(mol*s)", 3: "m^6/(mol^2*s)"}


def load_chem_yaml(path):
    """chem.yaml -> ([Species], [Reaction]) preserving the yaml reactions
    order (== the artifact's cantera.index space). Plain Arrhenius only."""
    with open(path) as fh:
        data = yaml.safe_load(fh)
    species = [_species_from_yaml(e) for e in data.get("species", [])]
    by_name = {s.label: s for s in species}
    reactions = []
    for entry in data.get("reactions", []):
        if "rate-constant" not in entry or "type" in entry:
            raise NotImplementedError(
                f"only elementary Arrhenius reactions are supported in v1; "
                f"offending entry: {entry.get('equation')}")
        r_names, p_names, reversible = _parse_equation(entry["equation"])
        rc = entry["rate-constant"]
        order = len(r_names)
        kin = Arrhenius(A=(float(rc["A"]), _A_UNITS_BY_ORDER[order]),
                        n=float(rc["b"]), Ea=(float(rc["Ea"]), "J/mol"),
                        T0=(1.0, "K"))
        rxn = Reaction(reactants=[by_name[n] for n in r_names],
                       products=[by_name[n] for n in p_names],
                       kinetics=kin, reversible=reversible,
                       duplicate=bool(entry.get("duplicate", False)))
        reactions.append(rxn)
    return species, reactions


def _restamp_and_extend(artifact, species, reactions):
    """Restore the solver stamps from artifact entries; reconstruct dropped
    (cantera-null) reactions from their recorded stoichiometry + kinetics."""
    by_name = {s.label: s for s in species}
    all_reactions = list(reactions)
    restamped_indices = set()
    for e in artifact["reactions"]:
        arch_name = e["archetype"]
        if arch_name not in ARCHETYPE_INTS:
            raise ValueError(
                f"reactions[] entry {e.get('id')!r} carries unknown "
                f"archetype {arch_name!r}; this loader's term-type "
                f"vocabulary is CLOSED to {sorted(ARCHETYPE_INTS)} "
                f"(docs/polymer_moments_format.md §3). An unknown archetype "
                f"is flux this consumer cannot reproduce -- upgrade the "
                f"loader or regenerate the sidecar with a matching RMG-Py "
                f"polymer branch.")
        arch = ARCHETYPE_INTS[arch_name]
        if e["cantera"] is not None:
            rxn = reactions[e["cantera"]["index"]]
            restamped_indices.add(e["cantera"]["index"])
        else:
            kin = e["kinetics"]
            if kin is None:
                raise ValueError(f"cantera-null entry {e['id']} carries no kinetics")
            order = len(e["reactants"])
            rxn = Reaction(
                reactants=[by_name[n] for n in e["reactants"]],
                products=[by_name[n] for n in e["products"]],
                kinetics=Arrhenius(A=(kin["A"], _A_UNITS_BY_ORDER[order]),
                                   n=kin["n"], Ea=(kin["Ea"], "J/mol"),
                                   T0=(1.0, "K")),
                reversible=bool(kin["reversible"]))
            all_reactions.append(rxn)
        rxn.polymer_flux_archetype = arch
        rxn.is_end_group_reaction = (e["scaling"] == "mu0")
        rxn.polymer_chip_units = int(e.get("params", {}).get("a", 0))
        # VOLATILE_EJECTION rows carry the SIGNED atom-transfer stamp
        # (params.eject_units; the only VE params sub-shape the emitter
        # writes). The oracle reads it back as Reaction.polymer_eject_units
        # (polymer.pyx:1551): cross-pool, the from-leg drops the full
        # interior/end-group bundle and the to-leg lands a-shifted
        # ((b0, b1 + sa*b0, b2 + 2*sa*b1 + a^2*b0), sa = -a forward / +a
        # reverse); same-pool (src == dst fold-back), the chip-style signed
        # mu1/mu2 single-pool write applies. src/dst pools are re-resolved
        # from the species at initialize_model like every other archetype
        # (the artifact's src_pool/dst_pool fields stay documentation).
        # Restore exactly; validate, never default (see
        # _validated_eject_units).
        if arch_name == "volatile_ejection/1":
            rxn.polymer_eject_units = _validated_eject_units(e)
        # Physically-melt classification, like MW, must cross the artifact
        # boundary: the reference-state tripwire's tag branch reads
        # is_polymer_proxy, which generation world stamps by blanket-tagging
        # every participant of a proxy-touching reaction (family.py:1657).
        # Reproduce that fingerprint from the artifact's proxy_* fields; the
        # tripwire's chain-scale MW conjunct (part of the class DEFINITION)
        # filters out small co-participants (H2, H...) identically to
        # generation world.
        #
        # MEMBERSHIP DIVERGENCE (known, fail-loud). This reconstruction is
        # NOT equivalent to generation-world membership:
        #   restorable    -- tags whose source is an ENTRY-LISTED reaction
        #                    (participants of artifact reactions[] entries
        #                    with non-empty proxy_* fields), i.e. exactly the
        #                    fingerprint above.
        #   NOT restorable -- three divergence classes whose tag source never
        #                    reaches the artifact: (1) explicit-oligomer
        #                    reactions (oligomer participants are not pool
        #                    proxies, so no entry is compiled), (2) spawned
        #                    UNCONFIGURED daughter proxies (registry pools
        #                    without solver configs; their reactions are not
        #                    entry-listed under configured_pools), and
        #                    (3) edge-reaction tag sources (only CORE
        #                    reactions compile to entries).
        # Direction is FAIL-LOUD: a chain-scale species that loses its tag
        # here is treated as gas, unpairs its condensed counterparty, and the
        # tripwire REFUSES (U ~ 11.4 decades) on chemistry the generation-
        # world tripwire scored 0.0 -- a spurious refusal, never a silent
        # acceptance. (Silent pair-down would require a generation run that
        # used allow_unpaired_reference_state AND a non-family pool-touching
        # reaction: contrived.) The robust alternative -- a per-species
        # is_polymer_proxy flag carried IN the artifact -- is a schema
        # addition awaiting a design decision; until then, diagnose a
        # consumer-only refusal as melt-classification divergence (see
        # docs/polymer_moments_format.md §8), not as a thermo problem.
        if e.get("proxy_reactants") or e.get("proxy_products"):
            for s in rxn.reactants + rxn.products:
                s.is_polymer_proxy = True
        # Refused-row marker (schema 2.4, format doc §12): restore the
        # generating solver's stamp so initialize_model rebuilds
        # reaction_refused and the residual suppresses the row's WHOLE flux
        # (moment and species alike) exactly like the generating run. The
        # accumulating class is recovered from the census reason (the same
        # bijection polymer.pyx:1526-1530 serialized it with). The row's
        # Cantera entry (when retained) needs no zeroing here -- the runner
        # drives the oracle directly, and the oracle's suppression is
        # whole-reaction.
        if _validate_refused_entry(e):
            rxn.polymer_refused = True
            rxn.polymer_refused_accumulating = (
                e["refused_reason"] == "qssa-invalid")
        if e["kinetics"] is not None:
            # Belt-and-braces reversibility for listed entries: post-fix
            # chem.yaml records reversibility in the equation arrow
            # (rmgpy/cantera.py get_reaction_equation emits '=>'/'<=>'), so this
            # now AGREES with the parsed arrow. The artifact's
            # kinetics.reversible records what the generating solver used and is
            # retained as a defence-in-depth restore (and recovers reversibility
            # from PRE-FIX all-'<=>' chem.yaml files for listed entries).
            rxn.reversible = bool(e["kinetics"]["reversible"])
    # Non-artifact (ordinary) reactions get their reversibility from the arrow
    # in post-fix chem.yaml. Warn only so that PRE-FIX all-'<=>' files — where
    # such a reaction would silently run reversible — are flagged.
    n_unrestamped = len(reactions) - len(restamped_indices)
    if n_unrestamped > 0:
        logging.warning(
            "polymer_moments_runner: %d chem.yaml reaction(s) have no artifact "
            "entry; their reversibility comes from the equation arrow "
            "('=>'/'<=>'). Post-fix chem.yaml round-trips this faithfully, but a "
            "PRE-FIX file (written before the cantera '=>'-for-irreversible fix) "
            "records every equation as '<=>', so an originally-irreversible "
            "ordinary reaction in such a file would run REVERSIBLE here (see "
            "module docstring).",
            n_unrestamped,
        )
    return all_reactions


# Pinned A-units per radical_qssa_unzip Arrhenius block. Deliberately pinned
# HERE, independently of the emitter's RADICAL_QSSA_SIDECAR_A_UNITS
# (rmgpy/polymer.py): the loader guards the artifact BOUNDARY, so a sidecar
# claiming any other unit system -- e.g. ktr in m^3/(mol*s), i.e. a
# bimolecular transfer constant that was never premultiplied by the substrate
# concentration -- must ERROR, never be silently converted or accepted.
_QSSA_PINNED_A_UNITS = {
    "initiation": "s^-1",
    "depropagation": "s^-1",
    "termination": "m^3/(mol*s)",
    "transfer": "s^-1",
    # Weak-link vocabulary (schema 2.2): allylic initiation is unimolecular
    # in the active unsaturated end; the split termination triplets carry
    # the bimolecular units of the summed block they replace.
    "initiation_allyl": "s^-1",
    "termination_recombination": "m^3/(mol*s)",
    "termination_disproportionation": "m^3/(mol*s)",
}
# Pinned units note for the U-state initial condition (schema 2.2). U is
# STATE, not a rate constant: same amount basis as mu0 [mol]; the consumer
# divides by V_poly. Pinned byte-for-byte against the emitter's
# RADICAL_QSSA_SIDECAR_U0_UNITS -- a sidecar claiming any other basis note
# (e.g. a bare "mol", or a per-volume unit) must ERROR, never be adapted.
_QSSA_U0_PINNED_UNITS = (
    "mol — tail-distribution state; consumer divides by V_poly")
# The weak-link allyl/U-state vocabulary (schema 2.2): all-or-nothing as a
# group, mutually exclusive with the legacy summed 'termination' (enforced
# by the shared validator; the loader only routes the keys through).
_QSSA_WEAKLINK_KEYS = frozenset(
    ("initiation_allyl", "termination_recombination",
     "termination_disproportionation", "unsaturated_tail_ends_initial"))
# Every key the sidecar block may carry: the channel-config vocabulary
# (validate_radical_qssa_unzip) plus the sidecar-only envelope fields.
_QSSA_SIDECAR_KEYS = frozenset(_QSSA_PINNED_A_UNITS) | frozenset(
    ("enabled", "basis", "efficiency", "monomer_yield", "provenance",
     "recipe", "unsaturated_tail_ends_initial"))

# Normative machine-readable QSSA recipe, pinned HERE independently of the
# emitter's RADICAL_QSSA_SIDECAR_RECIPE (rmgpy/polymer.py) -- same idiom as
# the units pin above: the loader guards the artifact BOUNDARY, so a sidecar
# claiming a different rate algebra (or omitting the recipe entirely) must
# ERROR, never be adapted to. Each string matches the implemented RHS in
# rmgpy/solver/polymer.pyx (M2 rate law; SMALL_EPS = 1e-30 at polymer.pyx:71).
_QSSA_PINNED_RECIPE = {
    "bond_basis": ("B = max(mu1 - mu0, 0) on concentration moments "
                   "(mol/m^3 condensed)"),
    "rate_no_transfer": ("r_mono = monomer_yield * kdp * "
                         "sqrt(efficiency * ki * B / kt)"),
    "rate_with_transfer": ("r_mono = monomer_yield * kdp * "
                           "(sqrt(ktr^2 + 8*kt*(2*efficiency*ki*B)) "
                           "- ktr) / (4*kt)"),
    "moment_signature": ("dmu0 = 0; dmu1 -= r_mono; dmu2 -= r_mono * "
                         "max(2*mu1/max(mu0, small_eps) - 1, 0)"),
    "small_eps": 1e-30,
    "volume_note": ("kt is bimolecular: rates depend on condensed "
                    "volume V_poly; consumers MUST evaluate on "
                    "concentration moments mu_k = n_k / V_poly and "
                    "convert emitted rate back with *V_poly"),
}

# Normative machine-readable recipe for the WEAK-LINK allyl/U-state variant
# (schema 2.2), pinned HERE independently of the emitter's
# RADICAL_QSSA_SIDECAR_RECIPE_WEAKLINK (rmgpy/polymer.py) -- same idiom as
# the legacy pin above. Each string matches the implemented weak-link RHS
# in rmgpy/solver/polymer.pyx (weak-link branch of the pool loop; U0 census
# in set_initial_conditions).
_QSSA_PINNED_RECIPE_WEAKLINK = {
    "bond_basis": ("B = max(mu1 - mu0, 0) on concentration moments "
                   "(mol/m^3 condensed)"),
    "channel_gate": ("channel gates on B > 0: at B = 0 it is INERT even at "
                     "U > 0 (no monomer drain, no U production, no U sink; "
                     "the DP->1 self-termination floor survives)"),
    "u_state": ("U is a per-pool amount state [mol], tail-distribution "
                "basis, MASSLESS (it never enters condensed mass or any "
                "Mn/Mw/PDI consumer); daughter pools spawn with U0 = 0 "
                "(constants inherit, state resets)"),
    "u_active": ("u_active = min(max(U, 0)/V_poly, B) (active-site clamp, "
                 "mol/m^3)"),
    "radical_generation": ("G_R = 2*efficiency*ki*B + "
                           "1*efficiency*ki_allyl*u_active (nu = 1: one "
                           "unzipping radical per weak-link fission; the "
                           "allylic co-fragment does not unzip)"),
    "kt_split": ("kt_total = kt_rec + kt_disp replaces the legacy summed kt "
                 "everywhere (halved radical-disappearance convention "
                 "unchanged)"),
    "rate_no_transfer": ("r_mono = monomer_yield * kdp * "
                         "sqrt(G_R / (2*kt_total))"),
    "rate_with_transfer": ("r_mono = monomer_yield * kdp * "
                           "(sqrt(ktr^2 + 8*kt_total*G_R) - ktr) / "
                           "(4*kt_total)"),
    "du_dt": ("dU/dt = kt_disp*R_ss^2*max(0, 1 - U/(2*mu0))*V_poly - "
              "efficiency*ki_allyl*u_active*V_poly; production is the "
              "disproportionation EVENT rate with NO efficiency factor "
              "(R_ss already carries the escape efficiency) under a linear "
              "capacity throttle exactly zero at the TAIL chain-end "
              "capacity 2*mu0 [mol], where mu0 is the INSTANTANEOUS "
              "tail-distribution mu0(t) amount read from the current "
              "state at each RHS evaluation (NOT the frozen initial mu0 "
              "of the u0_census bound); at mu0(t) <= 0 the throttle is "
              "EXACTLY 0 by explicit branch (no ends, no capacity, no U "
              "production; no division is performed and no small_eps "
              "floor is applied); the sink is efficiency-SYMMETRIC "
              "with the G_R allyl term (a caged recombination restores "
              "the allylic bond)"),
    "u0_census": ("U0 = unsaturated_tail_ends_initial [mol] is rejected at "
                  "solver init when U0 > 2*mu0_tail evaluated on the "
                  "INITIAL (t = 0) tail mu0 (TAIL-only chain-end capacity, "
                  "amount basis); the loader passes the value through"),
    "u_transport": ("U is NEVER advected or transferred by any inter-pool "
                    "or ejection flux archetype (migration, "
                    "scission_fragment, discrete_chip and "
                    "volatile_ejection move pool moments and species "
                    "amounts only); the ONLY writers of U are this "
                    "recipe's du_dt law and the t = 0 initial condition, "
                    "and daughter pools spawn with U0 = 0 (constants "
                    "inherit, state resets)"),
    "moment_signature": ("dmu0 = 0; dmu1 -= r_mono; dmu2 -= r_mono * "
                         "max(2*mu1/max(mu0, small_eps) - 1, 0)"),
    "small_eps": 1e-30,
    "volume_note": ("kt is bimolecular: rates depend on condensed "
                    "volume V_poly; consumers MUST evaluate on "
                    "concentration moments mu_k = n_k / V_poly and "
                    "convert emitted rate back with *V_poly"),
}

# Normative explicit-DP handshake recipe (schema 2.3), pinned HERE
# independently of the emitter's EXPLICIT_DP_SIDECAR_RECIPE
# (rmgpy/polymer.py) -- same boundary-guard idiom as the QSSA pins: a
# sidecar claiming a different handshake algebra (or omitting the recipe)
# must ERROR, never be adapted to. Each string matches the implemented
# oracle law in rmgpy/solver/polymer.pyx (hybrid handshake; gamma helpers
# _gamma_params_from_mu012 / _gamma_prob_conditional_hybrid) and the
# generation-side t=0 accounting (PolymerPhase.calculate_volume /
# set_initial_conditions).
_EXPLICIT_DP_PINNED_RECIPE = {
    "tail_split": ("declared pools[].moments are TOTAL-INCLUSIVE (explicit "
                   "chains counted in); the solver seeds the explicit "
                   "species' initial_moles as species amounts (clamped "
                   ">= 0), then subtracts each mapped DP's contribution "
                   "(N_dp, dp*N_dp, dp^2*N_dp on concentration moments) "
                   "from the seeded mu0/mu1/mu2, clamped >= 0 with a clamp "
                   "warning (set_initial_conditions step 6) -- the "
                   "integrated tail moments are total - explicit; the "
                   "generation-side V_poly mass split applies the same "
                   "subtraction and hard-errors when explicit mu1 exceeds "
                   "declared mu1 beyond -1e-12"),
    "boundary_flux": ("gated on mu0 > 1e-9 mol/m^3 AND mu1/mu0 > xs + 1e-9; "
                      "p_cond = P(DP = xs+1 | DP > xs) from the gamma "
                      "distribution moment-matched to (mu0, mu1, mu2) "
                      "(k = 1/(PDI - 1), theta = mean/k; half-integer bins: "
                      "[F((xs+1.5)/theta) - F((xs+0.5)/theta)] / "
                      "[1 - F((xs+0.5)/theta)] with F the regularized lower "
                      "incomplete gamma); triangular fallback on tail_mean "
                      "in (xs+1, xs+2) peaking 1.0 at xs+1.5 when the gamma "
                      "is unrealizable (any moment <= 1e-30, PDI <= 1+1e-6, "
                      "or non-finite params); p_cond clamped to [0, 1]; "
                      "N_boundary = min(mu0*p_cond, mu0, mu1/xs, mu2/xs^2) "
                      "[mol/m^3]; F_flux = k_chain * N_boundary; "
                      "dn(species[xs])/dt += F_flux*V_poly; dmu0 -= F_flux; "
                      "dmu1 -= xs*F_flux; dmu2 -= xs^2*F_flux"),
    "k_chain": ("k_unzip when k_unzip > 0, else r_qssa/max(mu0, 1e-30) when "
                "the radical_qssa_unzip channel is active (mutually "
                "exclusive upstream; at most one arm fires)"),
    "transport": ("one-way: statistical moment tail -> explicit real "
                  "condensed species at DP == handshake_target_dp; no "
                  "reverse flux"),
}
# The pool block's own revision token, pinned byte-for-byte against the
# emitter's EXPLICIT_DP_BLOCK_RECIPE_REVISION.
_EXPLICIT_DP_PINNED_REVISION = "2026-07-04-explicit-dp"
# Closed key vocabulary of the pool-level explicit_dp block (schema 2.3).
_EXPLICIT_DP_BLOCK_KEYS = frozenset(
    ("enabled", "species", "initial_moles", "handshake_target_dp",
     "recipe_revision", "recipe"))

# --- Pinned homolysis_initiation contract (schema 2.6, Stage 2 of the
# radical-homolysis initiation arc, adjudicated rounds 66/67). Pinned HERE
# independently of the emitter (rmgpy/polymer.py HOMOLYSIS_*), the same
# boundary-guard idiom as the QSSA/explicit-DP pins: the loader guards the
# artifact BOUNDARY -- a sidecar claiming a different kernel, recipe, or
# unit system must ERROR, never be adapted to. Each recipe string matches
# the implemented oracle law (rmgpy/solver/polymer.pyx, the khom_* RHS
# section: stable direct product forms, round-67 P2).
_HOMOLYSIS_PINNED_KERNEL = "radical_homolysis_initiation/1"
_HOMOLYSIS_PINNED_REVISION = "2026-07-05-radical-homolysis"
_HOMOLYSIS_BLOCK_KEYS = frozenset(
    ("enabled", "kinetics", "open_site_1_radical_pool",
     "open_site_2_radical_pool", "kernel", "recipe_revision", "recipe"))
_HOMOLYSIS_KINETICS_KEYS = frozenset(("A", "n", "Ea", "units"))
# Pinned units: the kernel triplet is SI (A [s^-1] -- unimolecular event
# rate per breakable bond mole; Ea [J/mol]), the radical_qssa_unzip
# convention. Any other claim is malformed.
_HOMOLYSIS_PINNED_UNITS = {"A": "s^-1", "Ea": "J/mol"}
# Ratified daughter-label convention (round 66; solver
# K_HOMOLYSIS_DAUGHTER_SUFFIXES): POSITIONAL open-*1 / open-*2 termini.
_HOMOLYSIS_DAUGHTER_SUFFIXES = ("_rad_primary_end", "_rad_secondary_end")
# The spawn provenance the Stage-1 producer stamps on both daughters
# (Polymer.generate_end_radical_daughters).
_HOMOLYSIS_SPAWN_SOURCE = "k_homolysis_end_radical"
_HOMOLYSIS_PINNED_RECIPE = {
    "event_rate": ("R = k(T)*max(mu1 - mu0, 0) [mol/(m^3 s)]; "
                   "k(T) = A*T^n*exp(-Ea/(R_gas*T)) evaluated at the "
                   "RUNTIME temperature (round 66: never a precomputed "
                   "scalar); a chain of length n has n-1 breakable bonds"),
    "parent_debit": ("dmu0 -= R; dmu1 -= R*B1 computed as k*(mu2 - mu1); "
                     "dmu2 -= R*B2 computed as k*(mu3 - mu2); mu3 from the "
                     "log_lagrange/1 closure"),
    "daughter_credit": ("EACH of the two end-radical daughter pools: "
                        "dmu0 += R; dmu1 += k*(mu2 - mu1)/2; "
                        "dmu2 += k*(2*mu3 - 3*mu2 + mu1)/6 -- STABLE "
                        "direct forms (round-67 P2): never B1/B2 ratios "
                        "re-multiplied by R (catastrophic cancellation "
                        "near DP -> 1 exhaustion)"),
    "totals": ("dmu0_total = +R = k(T)*max(mu1 - mu0, 0); dmu1_total = 0 "
               "(mass conserved, machine precision); dmu2_total = "
               "k*(mu1 - mu3)/3 = -k*(mu3 - mu1)/3 (the legacy "
               "Ziff-McGrady random-scission second-moment source)"),
    "out_of_domain": ("zero flux (warn once per pool per rebuild) when "
                      "mu2 < mu1, mu3 is nonfinite, or "
                      "2*mu3 - 3*mu2 + mu1 < 0 (B1, B2 >= 0 does NOT imply "
                      "a nonnegative daughter mu2 credit; round-67 P1)"),
    "reversibility": ("one-way, NO gas release: homolysis releases no "
                      "volatiles, and recombination arrives via the "
                      "discovered chemistry conduit, never this kernel"),
}

# --- Pinned side_group_homolysis contract (schema 2.7, FR1-K2, the
# sibling of the 2.6 homolysis_initiation surface; adjudicated rounds
# 72/73). Pinned HERE independently of the emitter (rmgpy/polymer.py
# SIDE_GROUP_*), the same boundary-guard idiom: the loader guards the
# artifact BOUNDARY -- a sidecar claiming a different kernel, recipe,
# unit system, selector vocabulary or feature-pool label convention must
# ERROR, never be adapted to. Each recipe string matches the implemented
# oracle law (rmgpy/solver/polymer.pyx, the sgh_* RHS section) and the
# NORMATIVE mass formula (PolymerPoolConfig.condensed_mass_g).
_SIDE_GROUP_PINNED_KERNEL = "side_group_homolysis/1"
_SIDE_GROUP_PINNED_REVISION = "2026-07-06-side-group-homolysis"
_SIDE_GROUP_BLOCK_KEYS = frozenset(
    ("enabled", "channels", "kernel", "recipe_revision", "recipe"))
_SIDE_GROUP_CHANNEL_KEYS = frozenset(
    ("label", "kinetics", "site_selector", "sites_per_unit",
     "site_atom_indices", "gas_product", "gas_species", "gas_mw_g_mol",
     "feature_pool"))
_SIDE_GROUP_KINETICS_KEYS = frozenset(("A", "n", "Ea", "units"))
# Pinned units: the kernel triplet is SI PER SITE (the RHS multiplies by
# sites_per_unit*mu1) -- deliberately distinct from the sibling kernel's
# plain "s^-1". Any other claim is malformed.
_SIDE_GROUP_PINNED_UNITS = {"A": "s^-1 per site", "Ea": "J/mol"}
# The CLOSED round-72 structural selector vocabulary (solver
# SIDE_GROUP_SITE_SELECTORS): the three classes PARTITION the carbon-bound
# X sites of a repeat unit. An unknown selector is a site this consumer
# cannot classify -- reject, never adapt.
_SIDE_GROUP_SELECTORS = ("aryl", "benzylic", "aliphatic")
# Ratified X-loss feature-pool label convention (solver
# side_group_daughter_pool_label): '{parent}_sidegrp_{sanitized label}',
# with the sanitizer mapping every character outside [A-Za-z0-9_] to '_'.
_SIDE_GROUP_DAUGHTER_INFIX = "_sidegrp_"
# The spawn provenance the FR1-K1 producer stamps on every X-loss feature
# daughter (Polymer.generate_side_loss_daughters): source + the SPAWNING
# CHANNEL's label.
_SIDE_GROUP_SPAWN_SOURCE = "side_group_homolysis"
_SIDE_GROUP_PINNED_RECIPE = {
    "event_rate": ("per channel: R = k(T)*s*mu1 [mol/(m^3 s)]; "
                   "k(T) = A*T^n*exp(-Ea/(R_gas*T)) per site evaluated at "
                   "the RUNTIME temperature; s = sites_per_unit -- every "
                   "repeat unit carries s X sites, so a chain of length n "
                   "reacts at k*s*n and the reacting chain is picked with "
                   "probability ~ n (site-weighted)"),
    "parent_debit": ("dmu_j -= k*s*mu_{j+1} for j = 0, 1, 2; mu3 from the "
                     "log_lagrange/1 closure -- STABLE direct product "
                     "forms, never picked-chain ratios re-multiplied "
                     "by R"),
    "feature_credit": ("the channel's X-loss feature pool "
                       "('{parent}_sidegrp_{sanitized channel label}') is "
                       "credited EXACTLY the parent debit -- the chain "
                       "transfers INTACT (no chain cut), two-pool "
                       "mu0/mu1/mu2 totals conserved bitwise, arriving-"
                       "flux mean length = the parent's length-biased "
                       "mean mu2/mu1"),
    "gas_credit": ("dn(gas_product) += R via the small-source -> "
                   "dn_dt*V_poly path (bit-identical to the feature mu0 "
                   "credit); channels sharing one gas species are "
                   "additive"),
    "out_of_domain": ("zero flux (warn once per (pool, channel) per "
                      "rebuild) when mu2 < mu1, mu3 is nonfinite, or "
                      "mu3 < mu2; mu1 <= 0 / mu0 <= 0 is the silent "
                      "in-domain zero (no sites, no events); a degenerate "
                      "k(T) raises"),
    "mass": ("condensed_mass_g = mu1*monomer_mw_g_mol - "
             "mu0*chain_mass_defect_g_mol -- NORMATIVE feature-pool mass "
             "formula: the feature pool keeps the parent's monomer_mw "
             "(intact-chain transfer) and every feature chain lost "
             "exactly ONE X (v1), so chain_mass_defect_g_mol = M_X of the "
             "spawning channel's gas_product and d(condensed mass)/dt + "
             "d(gas X mass)/dt = 0 exactly on the kernel's "
             "contributions"),
    "saturation": ("v1 LIMITATION: the kernel acts on the PARENT pool "
                   "only; feature pools carry no side_group_homolysis of "
                   "their own and saturate as terminal X-loss sinks (no "
                   "multi-loss cascade)"),
}


def _side_group_sanitize(channel_label):
    """The ratified channel-label sanitizer (pinned here independently of
    the solver's sanitize_side_group_channel_label, boundary-guard idiom):
    every character outside [A-Za-z0-9_] becomes '_'."""
    return "".join(c if (c.isalnum() or c == "_") else "_"
                   for c in channel_label)


# The QSSA channel vocabulary entered the sidecar schema at 2.1 (channel-
# vocabulary growth = minor bump); the emitter stamps >= 2.1 whenever it
# writes the block, so a 2.0 artifact carrying one is malformed. The
# weak-link allyl/U-state sub-vocabulary entered at 2.2, and the pool-level
# explicit_dp block at 2.3 -- same rule.
_QSSA_MIN_SCHEMA_MINOR = 1
_WEAKLINK_MIN_SCHEMA_MINOR = 2
_EXPLICIT_DP_MIN_SCHEMA_MINOR = 3
_REFUSED_MIN_SCHEMA_MINOR = 4
_SPAWNED_MIN_SCHEMA_MINOR = 5
_HOMOLYSIS_MIN_SCHEMA_MINOR = 6
_SIDE_GROUP_MIN_SCHEMA_MINOR = 7
# The CLOSED refused_reason vocabulary (format doc §12): the emitter derives
# the reason bijectively from the accumulating stamp, so exactly these two
# strings can exist. _restamp_and_extend reconstructs the accumulating class
# FROM the reason (== "qssa-invalid", else non-accumulating), so an unknown
# reason would silently change solver semantics -- reject, never adapt.
REFUSED_REASONS = frozenset({"conduit-deferred", "qssa-invalid"})
# Maximum schema minor this loader implements. Weak-link milestone iv
# POLICY CHANGE (was minor-permissive): a newer-minor artifact may carry
# vocabulary OUTSIDE the channel blocks (new conventions, new pool fields)
# that the unknown-key guards here never inspect, so an older loader must
# fail loud on it instead of loading additively. Raised to 3 with the
# explicit-DP handshake block (stage B). Raised to 4 with the refused-row
# marker (format doc §12): this loader restores the marker onto the
# reconstructed reactions (Reaction.polymer_refused) so the oracle's
# reaction_refused suppression zeroes the row's whole flux exactly like the
# generating run. Raised to 5 with the spawned-pool closure (format doc
# §13): conventions.spawned_pools is classification vocabulary whose
# runtime effect rides conventions.condensed_species -- which this loader
# already honors verbatim for its phase mask -- while pool moment blocks
# stay keyed on configured_pools (spawned pools remain solver-inert, §2),
# so 2.5 acceptance is truthful. Raised to 6 with the pool-level
# homolysis_initiation block (Stage 2, rounds 66/67): this loader validates
# the block strictly (_check_homolysis_initiation /
# _parse_homolysis_initiation_block) and wires the kernel triplet into the
# reconstructed pool configs, so the rebuilt oracle's RHS carries the
# generating solver's homolysis flux and its supersession census re-runs
# (refused rows stay refused/zero-flux) -- 2.6 acceptance is truthful.
# Raised to 7 with the pool-level side_group_homolysis block + the X-loss
# feature-pool chain_mass_defect_g_mol mass contract (FR1-K2, rounds
# 72/73): this loader validates both strictly
# (_check_side_group_homolysis / _parse_side_group_homolysis_block --
# including selector/feature closure FROM SERIALIZED DATA ALONE via the
# per-channel site_atom_indices, never re-deriving from a monomer graph
# it does not have) and wires the channel list + gas routing + defect
# into the reconstructed pool configs, so the rebuilt oracle's RHS
# carries the generating solver's side-group flux, its exact mass
# contract re-enforces at initialize_model (_flatten_side_group_state)
# and its supersession census re-runs -- 2.7 acceptance is truthful;
# 2.8+ stays rejected.
_MAX_KNOWN_SCHEMA_MINOR = 7


def _check_schema_version_known(artifact):
    """Reject any artifact whose schema_version is not one this loader
    implements (2.0 <= 2.x <= 2.3)."""
    ver = str(artifact.get("schema_version", ""))
    parts = ver.split(".")
    ok = (len(parts) == 2 and parts[0] == "2" and parts[1].isdigit()
          and int(parts[1]) <= _MAX_KNOWN_SCHEMA_MINOR)
    if not ok:
        raise ValueError(
            f"artifact schema_version {ver!r} is not implemented by this "
            f"loader (known: 2.0 .. 2.{_MAX_KNOWN_SCHEMA_MINOR}). A newer "
            f"minor may carry vocabulary outside the channel blocks that "
            f"this loader would silently ignore -- upgrade the loader or "
            f"regenerate the sidecar with a matching RMG-Py polymer "
            f"branch.")


def _check_qssa_schema_version(artifact):
    """Reject an artifact carrying any channels.radical_qssa_unzip block
    under a schema_version below 2.1, or the weak-link allyl/U-state
    sub-vocabulary below 2.2 (or a non-2.x version either way).

    Scans ALL pool entries, configured or not: the vocabulary appearing
    anywhere means the artifact claims the corresponding shape. Artifacts
    with no QSSA block anywhere are untouched by this check (the envelope
    gate is _check_schema_version_known)."""
    blocks = [
        (p.get("channels") or {}).get("radical_qssa_unzip")
        for p in artifact.get("pools", [])
        if isinstance(p, dict)
        and "radical_qssa_unzip" in (p.get("channels") or {})]
    if not blocks:
        return
    ver = str(artifact.get("schema_version", ""))
    parts = ver.split(".")
    minor = (int(parts[1]) if len(parts) == 2 and parts[0] == "2"
             and parts[1].isdigit() else -1)
    if minor < _QSSA_MIN_SCHEMA_MINOR:
        raise ValueError(
            f"artifact schema_version {ver!r} cannot carry a "
            f"channels.radical_qssa_unzip block: the QSSA channel vocabulary "
            f"was introduced in schema 2.1 (channel-vocabulary growth is a "
            f"minor bump), and the emitter stamps >= 2.1 whenever it writes "
            f"the block. This artifact is malformed -- regenerate the "
            f"sidecar with a current RMG-Py polymer branch.")
    weaklink = sorted(set().union(*(
        _QSSA_WEAKLINK_KEYS & set(b) for b in blocks
        if isinstance(b, dict))))
    if weaklink and minor < _WEAKLINK_MIN_SCHEMA_MINOR:
        raise ValueError(
            f"artifact schema_version {ver!r} cannot carry the weak-link "
            f"allyl/U-state vocabulary {weaklink}: it was introduced in "
            f"schema 2.2, and the emitter stamps 2.2 whenever it writes it. "
            f"A 2.1-stamped artifact must not smuggle U keys (schema/"
            f"vocabulary consistency; never loaded permissively) -- "
            f"regenerate the sidecar with a current RMG-Py polymer branch.")


# Recipe-revision gate for the monomer-gas contract (P1-B, incident
# 2026-07-03). Revision tokens that implement the GAS routed-monomer
# semantics (docs/polymer_moments_format.md, revision note
# 2026-07-03-monomer-gas): the unzip/QSSA release target is a gas-phase
# species, never listed in conventions.condensed_species or
# pools[].phase_species. Pre-revision tokens implemented the CONDENSED
# routed-monomer semantics this loader no longer supports.
_MONOMER_GAS_RECIPE_REVISIONS = frozenset({
    "2026-07-03-monomer-gas",
    "2026-07-03-qssa-monomer-gas",
    "2026-07-03-weaklink-u-monomer-gas",
    # Explicit-DP re-stamps (schema 2.3, stage B): same 3-way channel-family
    # split, all implementing the GAS routed-monomer contract (the
    # -monomer-gas suffix is retained; the explicit-DP handshake adds
    # condensed-side algebra and does not touch the routed monomer's phase).
    "2026-07-04-explicit-dp-monomer-gas",
    "2026-07-04-explicit-dp-qssa-monomer-gas",
    "2026-07-04-explicit-dp-weaklink-u-monomer-gas",
})
_PRE_MONOMER_GAS_RECIPE_REVISIONS = frozenset({
    "2026-06-10",
    "2026-07-02",
    "2026-07-03-weaklink-u",
})


def _check_explicit_dp_schema_version(artifact):
    """Reject an artifact carrying any pool-level explicit_dp block under a
    schema_version below 2.3 (or a non-2.x version). Scans ALL pool entries,
    configured or not, mirroring _check_qssa_schema_version: the vocabulary
    appearing anywhere means the artifact claims the 2.3 shape."""
    carriers = [p.get("label") for p in artifact.get("pools", [])
                if isinstance(p, dict) and "explicit_dp" in p]
    if not carriers:
        return
    ver = str(artifact.get("schema_version", ""))
    parts = ver.split(".")
    minor = (int(parts[1]) if len(parts) == 2 and parts[0] == "2"
             and parts[1].isdigit() else -1)
    if minor < _EXPLICIT_DP_MIN_SCHEMA_MINOR:
        raise ValueError(
            f"artifact schema_version {ver!r} cannot carry a pool-level "
            f"explicit_dp block (pools {carriers}): the explicit-DP "
            f"handshake vocabulary was introduced in schema 2.3, and the "
            f"emitter stamps 2.3 whenever it writes it. This artifact is "
            f"malformed -- regenerate the sidecar with a current RMG-Py "
            f"polymer branch.")


def _check_refused_schema_version(artifact):
    """Reject an artifact carrying the refused-row vocabulary (a "refused"
    or "refused_reason" key on any reactions[] entry) under a schema_version
    below 2.4 (or a non-2.x version). Mirrors _check_explicit_dp_schema_version:
    the vocabulary appearing anywhere means the artifact claims the 2.4
    shape (format doc §12) -- an older parse would drop the marker silently
    and integrate flux the generating solver zeroed."""
    carriers = [e.get("id") for e in artifact.get("reactions", [])
                if isinstance(e, dict)
                and ("refused" in e or "refused_reason" in e)]
    if not carriers:
        return
    ver = str(artifact.get("schema_version", ""))
    parts = ver.split(".")
    minor = (int(parts[1]) if len(parts) == 2 and parts[0] == "2"
             and parts[1].isdigit() else -1)
    if minor < _REFUSED_MIN_SCHEMA_MINOR:
        raise ValueError(
            f"artifact schema_version {ver!r} cannot carry the refused-row "
            f"marker (entries {carriers}): the refused vocabulary was "
            f"introduced in schema 2.4, and the emitter stamps 2.4 whenever "
            f"it writes it. This artifact is malformed -- regenerate the "
            f"sidecar with a current RMG-Py polymer branch.")


def _check_spawned_pools_schema_version(artifact):
    """Reject an artifact carrying the spawned-pool closure vocabulary (a
    ``conventions.spawned_pools`` KEY -- presence, not truthiness) under a
    schema_version below 2.5 (or a non-2.x version), and shape-guard the
    key itself. Mirrors _check_refused_schema_version on the version axis
    and _validate_refused_entry on the shape axis: the emitter writes the
    key ONLY as a non-empty list of unique, non-empty pool-label strings
    and stamps 2.5 whenever it does (format doc §13), so an empty list, a
    non-list (a bare string would be consumed as an iterable of
    CHARACTERS), non-string / empty entries, or duplicate labels are all
    malformed -- reject, never adapt. Also rejects an overlap with
    configured_pools: the closure is the complement of the configured set
    by construction (a label in both would be simultaneously
    solver-configured and solver-inert -- the emitter never produces
    it)."""
    conv = artifact.get("conventions") or {}
    if "spawned_pools" not in conv:
        return
    spawned = conv["spawned_pools"]
    ver = str(artifact.get("schema_version", ""))
    parts = ver.split(".")
    minor = (int(parts[1]) if len(parts) == 2 and parts[0] == "2"
             and parts[1].isdigit() else -1)
    if minor < _SPAWNED_MIN_SCHEMA_MINOR:
        raise ValueError(
            f"artifact schema_version {ver!r} cannot carry "
            f"conventions.spawned_pools ({spawned!r}): the spawned-pool "
            f"closure vocabulary was introduced in schema 2.5, and the "
            f"emitter stamps 2.5 whenever it writes it. This artifact is "
            f"malformed -- regenerate the sidecar with a current RMG-Py "
            f"polymer branch.")
    if not isinstance(spawned, (list, tuple)):
        raise ValueError(
            f"conventions.spawned_pools is {spawned!r} "
            f"({type(spawned).__name__}), not a list: the emitter writes a "
            f"non-empty list of pool-label strings (format doc §13). This "
            f"artifact is malformed -- fix/regenerate it at the source; "
            f"this loader never adapts it.")
    if not spawned:
        raise ValueError(
            "conventions.spawned_pools is an empty list: the emitter "
            "writes the key ONLY when the closure is non-empty "
            "(presence-based 2.5 stamping, format doc §13), so an empty "
            "list is malformed -- fix/regenerate it at the source; this "
            "loader never adapts it.")
    for lbl in spawned:
        if not isinstance(lbl, str) or not lbl:
            raise ValueError(
                f"conventions.spawned_pools entry {lbl!r} is not a "
                f"non-empty string: the closure lists registry pool LABELS "
                f"(format doc §13). This artifact is malformed -- "
                f"fix/regenerate it at the source; this loader never "
                f"adapts it.")
    dupes = sorted({lbl for lbl in spawned if spawned.count(lbl) > 1})
    if dupes:
        raise ValueError(
            f"conventions.spawned_pools carries duplicate label(s) "
            f"{dupes}: the closure is a set of registry pool labels "
            f"(format doc §13) and the emitter never repeats one. This "
            f"artifact is malformed -- fix/regenerate it at the source; "
            f"this loader never adapts it.")
    overlap = sorted(set(spawned) & set(conv.get("configured_pools") or []))
    if overlap:
        raise ValueError(
            f"conventions.spawned_pools overlaps configured_pools on "
            f"{overlap}: the spawned closure is the complement of the "
            f"configured set (format doc §13; a pool cannot be both "
            f"solver-configured and runtime-spawned/inert). This artifact "
            f"is malformed -- fix/regenerate it at the source; this loader "
            f"never adapts it.")


def _parse_homolysis_initiation_block(lab, pool_entry):
    """Parse + validate a pool entry's pool-level homolysis_initiation
    block (schema 2.6). Returns the validated Arrhenius triplet
    ``{"A": float, "n": float, "Ea": float}`` (SI: A [s^-1], Ea [J/mol]),
    or ``None`` when the block is absent.

    Boundary rules mirror _parse_explicit_dp_block: closed key vocabulary,
    boolean ``enabled`` (present-disabled REJECTED -- the emitter never
    writes disabled blocks), key-presence + SHAPE validation (never
    truthiness: JSON ``true`` is not a rate constant, the r63 lesson), the
    pinned kernel name (round-67 ruling (c): an unknown kernel is flux this
    consumer cannot reproduce), pinned units (ERROR on mismatch, never
    convert), and the normative ``recipe`` + block ``recipe_revision``
    pinned by exact match (reject, never adapt)."""
    raw = pool_entry.get("homolysis_initiation")
    if raw is None:
        return None
    if not isinstance(raw, dict):
        raise ValueError(
            f"Pool {lab!r}: homolysis_initiation must be a dict, got "
            f"{type(raw).__name__}. Fix the artifact.")
    unknown = sorted(set(raw) - _HOMOLYSIS_BLOCK_KEYS)
    if unknown:
        raise ValueError(
            f"Pool {lab!r}: homolysis_initiation has unknown key(s) "
            f"{unknown}; allowed keys are {sorted(_HOMOLYSIS_BLOCK_KEYS)}. "
            f"Fix the artifact (unknown vocabulary is never dropped "
            f"permissively).")
    missing = sorted(_HOMOLYSIS_BLOCK_KEYS - set(raw))
    if missing:
        raise ValueError(
            f"Pool {lab!r}: homolysis_initiation is missing key(s) "
            f"{missing}. The emitter always writes the full block "
            f"(kinetics, both open-site daughter fields, kernel, "
            f"recipe_revision, recipe) -- regenerate the sidecar.")
    enabled = raw["enabled"]
    if not isinstance(enabled, bool):
        raise ValueError(
            f"Pool {lab!r}: homolysis_initiation must carry a boolean "
            f"'enabled' field, got {enabled!r}. Fix the artifact.")
    if not enabled:
        raise ValueError(
            f"Pool {lab!r}: homolysis_initiation carries enabled=false. "
            f"The emitter never writes disabled blocks: a disabled kernel "
            f"must be ABSENT from the sidecar, not present-disabled. Fix "
            f"the artifact (remove the block).")
    kernel = raw["kernel"]
    if kernel != _HOMOLYSIS_PINNED_KERNEL:
        raise ValueError(
            f"Pool {lab!r}: homolysis_initiation names kernel {kernel!r}; "
            f"this loader implements exactly "
            f"{_HOMOLYSIS_PINNED_KERNEL!r}. An unknown kernel is flux this "
            f"consumer cannot reproduce (round-67 ruling (c) supersession "
            f"contract) -- upgrade the loader or regenerate the sidecar.")
    if raw["recipe_revision"] != _HOMOLYSIS_PINNED_REVISION:
        raise ValueError(
            f"Pool {lab!r}: homolysis_initiation recipe_revision must "
            f"equal {_HOMOLYSIS_PINNED_REVISION!r} exactly, got "
            f"{raw['recipe_revision']!r}. An artifact claiming a different "
            f"kernel recipe must be fixed at the source; this loader "
            f"validates, never adapts.")
    recipe = raw["recipe"]
    if not isinstance(recipe, dict):
        raise ValueError(
            f"Pool {lab!r}: homolysis_initiation must carry the normative "
            f"'recipe' dict (schema 2.6), got {recipe!r} -- regenerate the "
            f"sidecar.")
    unknown_recipe = sorted(set(recipe) - set(_HOMOLYSIS_PINNED_RECIPE))
    if unknown_recipe:
        raise ValueError(
            f"Pool {lab!r}: homolysis_initiation recipe has unknown "
            f"key(s) {unknown_recipe}; allowed keys are "
            f"{sorted(_HOMOLYSIS_PINNED_RECIPE)}. Fix the artifact.")
    for key, pinned in _HOMOLYSIS_PINNED_RECIPE.items():
        if key not in recipe or recipe[key] != pinned:
            raise ValueError(
                f"Pool {lab!r}: homolysis_initiation recipe[{key!r}] must "
                f"equal the pinned normative recipe exactly; got "
                f"{recipe.get(key)!r}, expected {pinned!r}. An artifact "
                f"claiming a different kernel algebra must be fixed at the "
                f"source; this loader validates, never adapts.")
    kin = raw["kinetics"]
    if not isinstance(kin, dict):
        raise ValueError(
            f"Pool {lab!r}: homolysis_initiation kinetics must be a dict "
            f"{{A, n, Ea, units}}, got {type(kin).__name__}. Fix the "
            f"artifact.")
    unknown_kin = sorted(set(kin) - _HOMOLYSIS_KINETICS_KEYS)
    missing_kin = sorted(_HOMOLYSIS_KINETICS_KEYS - set(kin))
    if unknown_kin or missing_kin:
        raise ValueError(
            f"Pool {lab!r}: homolysis_initiation kinetics must carry "
            f"exactly the keys {sorted(_HOMOLYSIS_KINETICS_KEYS)} "
            f"(unknown: {unknown_kin}, missing: {missing_kin}). Fix the "
            f"artifact.")
    if kin["units"] != dict(_HOMOLYSIS_PINNED_UNITS):
        raise ValueError(
            f"Pool {lab!r}: homolysis_initiation kinetics units must be "
            f"exactly {_HOMOLYSIS_PINNED_UNITS!r} (SI, the "
            f"radical_qssa_unzip convention), got {kin['units']!r}. A "
            f"sidecar claiming any other unit system must ERROR, never be "
            f"silently converted.")
    out = {}
    for key in ("A", "n", "Ea"):
        val = kin[key]
        if isinstance(val, bool) or not isinstance(val, (int, float)) or \
                not math.isfinite(float(val)):
            raise ValueError(
                f"Pool {lab!r}: homolysis_initiation kinetics {key}="
                f"{val!r} must be a finite number (shape validation, "
                f"never truthiness). Fix the artifact.")
        out[key] = float(val)
    if out["A"] <= 0.0:
        raise ValueError(
            f"Pool {lab!r}: homolysis_initiation kinetics A={out['A']:g} "
            f"must be > 0 (a zero/negative kernel must be ABSENT, not "
            f"present-inert). Fix the artifact.")
    if out["Ea"] < 0.0:
        raise ValueError(
            f"Pool {lab!r}: homolysis_initiation kinetics "
            f"Ea={out['Ea']:g} must be >= 0 [J/mol]. Fix the artifact.")
    for field, suffix in zip(("open_site_1_radical_pool",
                              "open_site_2_radical_pool"),
                             _HOMOLYSIS_DAUGHTER_SUFFIXES):
        d_label = raw[field]
        if not isinstance(d_label, str) or d_label != f"{lab}{suffix}":
            raise ValueError(
                f"Pool {lab!r}: homolysis_initiation {field}={d_label!r} "
                f"must be the ratified daughter pool label "
                f"{lab + suffix!r} (round-66 POSITIONAL open-site "
                f"convention). Fix the artifact.")
    return out


def _check_homolysis_initiation(artifact):
    """Vocabulary/version cross-check + closure guard for the schema-2.6
    pool-level homolysis_initiation block. Mirrors
    _check_explicit_dp_schema_version on the version axis (the block
    anywhere under a below-2.6 stamp is malformed) and adds the daughter
    closure the emitter guarantees (round-67 §Stage 2 Scope, tightened by
    the round-68 adjudication): the kernel-carrying pool itself must be in
    conventions.configured_pools (build_system_from_artifact only builds
    configured pools -- a block on an unconfigured carrier would be a
    silently dropped kernel), and each kernel-carrying pool's two
    end-radical daughter pools must

    * be present in pools[] (their moment slots receive the kernel's
      fragment credits),
    * be in conventions.configured_pools (the daughters are eagerly
      solver-configured by design, polymer.pyx _flatten_homolysis_state
      hard-errors otherwise; build_system_from_artifact never builds a
      pool outside configured_pools, so any other classification has no
      solver home for the kernel's credits) and NOT in
      conventions.spawned_pools (schema 2.5 defines that list as the
      configured set's complement -- membership there contradicts the
      eager-configured daughter design),
    * be condensed per the condensed closure (phase_species non-empty and
      fully inside conventions.condensed_species -- the item-16
      mass-balance hazard otherwise), and
    * carry the Stage-1 spawn provenance
      (spawn_event_metadata.source == 'k_homolysis_end_radical'): the
      daughters are producer-spawned pools, never user pools."""
    carriers = [p for p in artifact.get("pools", [])
                if isinstance(p, dict) and "homolysis_initiation" in p]
    if not carriers:
        return
    ver = str(artifact.get("schema_version", ""))
    parts = ver.split(".")
    minor = (int(parts[1]) if len(parts) == 2 and parts[0] == "2"
             and parts[1].isdigit() else -1)
    if minor < _HOMOLYSIS_MIN_SCHEMA_MINOR:
        raise ValueError(
            f"artifact schema_version {ver!r} cannot carry a pool-level "
            f"homolysis_initiation block (pools "
            f"{[p.get('label') for p in carriers]}): the radical-homolysis "
            f"kernel vocabulary was introduced in schema 2.6, and the "
            f"emitter stamps 2.6 whenever it writes it. This artifact is "
            f"malformed -- regenerate the sidecar with a current RMG-Py "
            f"polymer branch.")
    conv = artifact.get("conventions") or {}
    configured = set(conv.get("configured_pools") or [])
    spawned_raw = conv.get("spawned_pools")
    spawned = (set(spawned_raw)
               if isinstance(spawned_raw, (list, tuple)) else set())
    condensed = set(conv.get("condensed_species") or [])
    by_label = {p.get("label"): p for p in artifact.get("pools", [])
                if isinstance(p, dict)}
    for carrier in carriers:
        lab = carrier.get("label")
        # Round-68 P1: the carrier itself must be solver-configured.
        # build_system_from_artifact only constructs PolymerPoolConfigs
        # for conventions.configured_pools, so a valid-looking block on an
        # unconfigured pool would be SKIPPED -- a silently dropped kernel.
        if lab not in configured:
            raise ValueError(
                f"Pool {lab!r}: carries a homolysis_initiation block but "
                f"is not listed in conventions.configured_pools. The "
                f"loader only builds configured pools, so the kernel "
                f"would be silently dropped (never integrated). The "
                f"producer only serializes the block on solver-configured "
                f"pools; this artifact is hand-edited/corrupted. Fix the "
                f"artifact.")
        _parse_homolysis_initiation_block(lab, carrier)  # shape guard
        block = carrier["homolysis_initiation"]
        for field in ("open_site_1_radical_pool",
                      "open_site_2_radical_pool"):
            d_label = block[field]
            d_entry = by_label.get(d_label)
            if d_entry is None:
                raise ValueError(
                    f"Pool {lab!r}: homolysis_initiation daughter pool "
                    f"{d_label!r} ({field}) is missing from the artifact's "
                    f"pools[] -- the kernel's fragment credits would have "
                    f"no moment slots. This artifact is malformed; "
                    f"regenerate the sidecar.")
            # Round-68 P1: daughters are eagerly solver-CONFIGURED by
            # design; spawned_pools is the configured set's complement
            # (schema 2.5), so membership there contradicts that design.
            # Check the spawned conjunct first so the contradiction is
            # named even when the 2.5 overlap guard was bypassed.
            if d_label in spawned:
                raise ValueError(
                    f"Pool {lab!r}: homolysis_initiation daughter pool "
                    f"{d_label!r} is classified in "
                    f"conventions.spawned_pools. Schema 2.5 defines that "
                    f"list as the configured set's complement, and "
                    f"homolysis daughters are eagerly solver-configured "
                    f"by design (polymer.pyx _flatten_homolysis_state "
                    f"hard-errors otherwise) -- a spawned-classified "
                    f"daughter is never built by the loader and the "
                    f"kernel's credits would have no solver home. Fix "
                    f"the artifact.")
            if d_label not in configured:
                raise ValueError(
                    f"Pool {lab!r}: homolysis_initiation daughter pool "
                    f"{d_label!r} is not listed in "
                    f"conventions.configured_pools. The daughters are "
                    f"eagerly solver-configured by the producer "
                    f"(polymer.pyx _flatten_homolysis_state hard-errors "
                    f"otherwise), and the loader only builds configured "
                    f"pools -- an unconfigured daughter has no solver "
                    f"home for the kernel's credits. Fix the artifact.")
            members = d_entry.get("phase_species") or []
            not_condensed = sorted(m for m in members if m not in condensed)
            if not members or not_condensed:
                raise ValueError(
                    f"Pool {lab!r}: homolysis_initiation daughter pool "
                    f"{d_label!r} is not condensed per the condensed "
                    f"closure (phase_species={members}, not condensed="
                    f"{not_condensed}); its moment slots would be "
                    f"phase-classified GAS (the item-16 mass-balance "
                    f"hazard). Fix the artifact.")
            meta = d_entry.get("spawn_event_metadata")
            if not isinstance(meta, dict) or \
                    meta.get("source") != _HOMOLYSIS_SPAWN_SOURCE:
                raise ValueError(
                    f"Pool {lab!r}: homolysis_initiation daughter pool "
                    f"{d_label!r} lacks the Stage-1 spawn provenance "
                    f"(expected spawn_event_metadata.source == "
                    f"{_HOMOLYSIS_SPAWN_SOURCE!r}, got {meta!r}). The "
                    f"daughters are producer-spawned pools, never user "
                    f"pools -- a row claiming other provenance is "
                    f"hand-edited. Fix the artifact.")


def _parse_side_group_homolysis_block(lab, pool_entry):
    """Parse + validate a pool entry's pool-level side_group_homolysis
    block (schema 2.7, FR1-K2). Returns the validated channel list -- one
    dict per channel with EXACTLY the serialized vocabulary (label,
    kinetics triplet flattened to A/n/Ea floats, site_selector,
    sites_per_unit, site_atom_indices, gas_product, gas_species,
    gas_mw_g_mol, feature_pool) -- or ``None`` when the block is absent.

    Boundary rules mirror _parse_homolysis_initiation_block: closed key
    vocabulary at every level, boolean ``enabled`` (present-disabled
    REJECTED), key-presence + SHAPE validation (never truthiness), pinned
    kernel/units/recipe/recipe_revision by exact match (reject, never
    adapt). Structural selector closure is validated FROM SERIALIZED DATA
    ALONE (round-73: K2 must not inherit the solver-backstop structural
    gap -- this loader has no monomer graph to re-derive from):

    * ``site_selector`` must come from the CLOSED vocabulary
      {'aryl','benzylic','aliphatic'};
    * ``site_atom_indices`` -- the emitter-resolved selector match
      indices in monomer_adj_list atom order -- must be a non-empty list
      of unique non-negative ints whose LENGTH EQUALS ``sites_per_unit``
      (the serialized rendering of the round-72 'sites_per_unit is
      CHECKED, never trusted' law);
    * no two channels may resolve to the SAME (gas element, atom set) --
      the double-carry distinct labels were hiding;
    * ``gas_product`` must parse to a monoatomic mono-radical whose molar
      mass matches the serialized ``gas_mw_g_mol`` pin (the value the
      feature pool's chain_mass_defect_g_mol must also pin);
    * ``feature_pool`` must be the ratified
      '{parent}_sidegrp_{sanitized label}' daughter label;
    * duplicate channel labels (raw or sanitized-collision) reject."""
    raw = pool_entry.get("side_group_homolysis")
    if raw is None:
        return None
    if not isinstance(raw, dict):
        raise ValueError(
            f"Pool {lab!r}: side_group_homolysis must be a dict, got "
            f"{type(raw).__name__}. Fix the artifact.")
    unknown = sorted(set(raw) - _SIDE_GROUP_BLOCK_KEYS)
    if unknown:
        raise ValueError(
            f"Pool {lab!r}: side_group_homolysis has unknown key(s) "
            f"{unknown}; allowed keys are "
            f"{sorted(_SIDE_GROUP_BLOCK_KEYS)}. Fix the artifact "
            f"(unknown vocabulary is never dropped permissively).")
    missing = sorted(_SIDE_GROUP_BLOCK_KEYS - set(raw))
    if missing:
        raise ValueError(
            f"Pool {lab!r}: side_group_homolysis is missing key(s) "
            f"{missing}. The emitter always writes the full block "
            f"(channels, kernel, recipe_revision, recipe) -- regenerate "
            f"the sidecar.")
    enabled = raw["enabled"]
    if not isinstance(enabled, bool):
        raise ValueError(
            f"Pool {lab!r}: side_group_homolysis must carry a boolean "
            f"'enabled' field, got {enabled!r}. Fix the artifact.")
    if not enabled:
        raise ValueError(
            f"Pool {lab!r}: side_group_homolysis carries enabled=false. "
            f"The emitter never writes disabled blocks: a disabled kernel "
            f"must be ABSENT from the sidecar, not present-disabled. Fix "
            f"the artifact (remove the block).")
    kernel = raw["kernel"]
    if kernel != _SIDE_GROUP_PINNED_KERNEL:
        raise ValueError(
            f"Pool {lab!r}: side_group_homolysis names kernel {kernel!r}; "
            f"this loader implements exactly "
            f"{_SIDE_GROUP_PINNED_KERNEL!r}. An unknown kernel is flux "
            f"this consumer cannot reproduce (the supersession contract) "
            f"-- upgrade the loader or regenerate the sidecar.")
    if raw["recipe_revision"] != _SIDE_GROUP_PINNED_REVISION:
        raise ValueError(
            f"Pool {lab!r}: side_group_homolysis recipe_revision must "
            f"equal {_SIDE_GROUP_PINNED_REVISION!r} exactly, got "
            f"{raw['recipe_revision']!r}. An artifact claiming a "
            f"different kernel recipe must be fixed at the source; this "
            f"loader validates, never adapts.")
    recipe = raw["recipe"]
    if not isinstance(recipe, dict):
        raise ValueError(
            f"Pool {lab!r}: side_group_homolysis must carry the normative "
            f"'recipe' dict (schema 2.7), got {recipe!r} -- regenerate "
            f"the sidecar.")
    unknown_recipe = sorted(set(recipe) - set(_SIDE_GROUP_PINNED_RECIPE))
    if unknown_recipe:
        raise ValueError(
            f"Pool {lab!r}: side_group_homolysis recipe has unknown "
            f"key(s) {unknown_recipe}; allowed keys are "
            f"{sorted(_SIDE_GROUP_PINNED_RECIPE)}. Fix the artifact.")
    for key, pinned in _SIDE_GROUP_PINNED_RECIPE.items():
        if key not in recipe or recipe[key] != pinned:
            raise ValueError(
                f"Pool {lab!r}: side_group_homolysis recipe[{key!r}] must "
                f"equal the pinned normative recipe exactly (including "
                f"the NORMATIVE mass formula, recipe['mass']); got "
                f"{recipe.get(key)!r}, expected {pinned!r}. An artifact "
                f"claiming a different kernel algebra or mass contract "
                f"must be fixed at the source; this loader validates, "
                f"never adapts.")
    channels = raw["channels"]
    if not isinstance(channels, list) or len(channels) == 0:
        raise ValueError(
            f"Pool {lab!r}: side_group_homolysis channels must be a "
            f"non-empty list of channel dicts, got {channels!r}. A "
            f"kernel with no channels must be ABSENT from the sidecar. "
            f"Fix the artifact.")
    # r75 P1-2: the serialized site_atom_indices are 0-based positions in
    # the carrier's monomer_adj_list atom order; without that text the
    # indices cannot be bounds-anchored AT ALL. The atom count is parsed
    # from the TEXT alone (one numbered atom line per atom -- the first
    # token is the 1-based atom number), never by re-deriving a molecule
    # graph (the round-73 posture).
    adj = pool_entry.get("monomer_adj_list")
    n_atoms = 0
    if isinstance(adj, str):
        n_atoms = sum(1 for line in adj.splitlines()
                      if line.split() and line.split()[0].isdigit())
    if n_atoms <= 0:
        raise ValueError(
            f"Pool {lab!r}: carries a side_group_homolysis block but its "
            f"monomer_adj_list is missing/empty (names {n_atoms} atoms) "
            f"-- site_atom_indices are indices in monomer_adj_list atom "
            f"order and cannot be bounds-checked without it (round-75 "
            f"P1). Fix the artifact (regenerate the sidecar).")
    # Lazy import (the boundary guard needs to parse the gas_product
    # SMILES to re-pin M_X and the gas element -- a one-atom molecule,
    # never a monomer graph).
    from rmgpy.molecule import Molecule
    out = []
    seen_sanitized = {}
    seen_sites = {}
    for pos, ch in enumerate(channels):
        if not isinstance(ch, dict):
            raise ValueError(
                f"Pool {lab!r}: side_group_homolysis channels[{pos}] must "
                f"be a channel dict, got {type(ch).__name__}. Fix the "
                f"artifact.")
        unknown_ch = sorted(set(ch) - _SIDE_GROUP_CHANNEL_KEYS)
        missing_ch = sorted(_SIDE_GROUP_CHANNEL_KEYS - set(ch))
        if unknown_ch or missing_ch:
            raise ValueError(
                f"Pool {lab!r}: side_group_homolysis channels[{pos}] "
                f"(label={ch.get('label', '<unset>')!r}) must carry "
                f"exactly the keys {sorted(_SIDE_GROUP_CHANNEL_KEYS)} "
                f"(unknown: {unknown_ch}, missing: {missing_ch}). Fix "
                f"the artifact.")
        label = ch["label"]
        if not isinstance(label, str) or not label.strip():
            raise ValueError(
                f"Pool {lab!r}: side_group_homolysis channels[{pos}] "
                f"label must be a non-empty string, got {label!r}. Fix "
                f"the artifact.")
        san = _side_group_sanitize(label)
        if san in seen_sanitized:
            raise ValueError(
                f"Pool {lab!r}: side_group_homolysis has duplicate "
                f"channel labels: {seen_sanitized[san]!r} and {label!r} "
                f"collide (sanitized {san!r}). Two channels of the same "
                f"bond class would double-carry the loss. Fix the "
                f"artifact.")
        seen_sanitized[san] = label
        kin = ch["kinetics"]
        if not isinstance(kin, dict):
            raise ValueError(
                f"Pool {lab!r}: side_group_homolysis channel {label!r} "
                f"kinetics must be a dict {{A, n, Ea, units}}, got "
                f"{type(kin).__name__}. Fix the artifact.")
        unknown_kin = sorted(set(kin) - _SIDE_GROUP_KINETICS_KEYS)
        missing_kin = sorted(_SIDE_GROUP_KINETICS_KEYS - set(kin))
        if unknown_kin or missing_kin:
            raise ValueError(
                f"Pool {lab!r}: side_group_homolysis channel {label!r} "
                f"kinetics must carry exactly the keys "
                f"{sorted(_SIDE_GROUP_KINETICS_KEYS)} (unknown: "
                f"{unknown_kin}, missing: {missing_kin}). Fix the "
                f"artifact.")
        if kin["units"] != dict(_SIDE_GROUP_PINNED_UNITS):
            raise ValueError(
                f"Pool {lab!r}: side_group_homolysis channel {label!r} "
                f"kinetics units must be exactly "
                f"{_SIDE_GROUP_PINNED_UNITS!r} (SI, PER SITE -- the RHS "
                f"multiplies by sites_per_unit*mu1), got {kin['units']!r}. "
                f"A sidecar claiming any other unit system must ERROR, "
                f"never be silently converted.")
        parsed = {"label": label}
        for key in ("A", "n", "Ea"):
            val = kin[key]
            if isinstance(val, bool) or not isinstance(val, (int, float)) \
                    or not math.isfinite(float(val)):
                raise ValueError(
                    f"Pool {lab!r}: side_group_homolysis channel "
                    f"{label!r} kinetics {key}={val!r} must be a finite "
                    f"number (shape validation, never truthiness). Fix "
                    f"the artifact.")
            parsed[key] = float(val)
        if parsed["A"] <= 0.0:
            raise ValueError(
                f"Pool {lab!r}: side_group_homolysis channel {label!r} "
                f"kinetics A={parsed['A']:g} must be > 0 (a zero/negative "
                f"channel must be ABSENT, not present-inert). Fix the "
                f"artifact.")
        if parsed["Ea"] < 0.0:
            raise ValueError(
                f"Pool {lab!r}: side_group_homolysis channel {label!r} "
                f"kinetics Ea={parsed['Ea']:g} must be >= 0 [J/mol]. Fix "
                f"the artifact.")
        sel = ch["site_selector"]
        if not isinstance(sel, str) or sel not in _SIDE_GROUP_SELECTORS:
            raise ValueError(
                f"Pool {lab!r}: side_group_homolysis channel {label!r} "
                f"site_selector={sel!r} must be one of "
                f"{_SIDE_GROUP_SELECTORS} -- the round-72 CLOSED "
                f"structural selector vocabulary. An unknown selector is "
                f"a site this consumer cannot classify; reject, never "
                f"adapt.")
        parsed["site_selector"] = sel
        spu = ch["sites_per_unit"]
        if isinstance(spu, bool) or not isinstance(spu, (int, float)) or \
                not math.isfinite(float(spu)) or float(spu) <= 0.0:
            raise ValueError(
                f"Pool {lab!r}: side_group_homolysis channel {label!r} "
                f"sites_per_unit={spu!r} must be a finite number > 0. "
                f"Fix the artifact.")
        parsed["sites_per_unit"] = float(spu)
        idxs = ch["site_atom_indices"]
        if (not isinstance(idxs, list) or len(idxs) == 0
                or any(isinstance(i, bool) or not isinstance(i, int)
                       or i < 0 for i in idxs)
                or len(set(idxs)) != len(idxs)):
            raise ValueError(
                f"Pool {lab!r}: side_group_homolysis channel {label!r} "
                f"site_atom_indices={idxs!r} must be a non-empty list of "
                f"unique non-negative ints (the emitter-resolved selector "
                f"match indices in monomer_adj_list atom order -- the "
                f"loader-side structural closure, round-73). Fix the "
                f"artifact.")
        if float(len(idxs)) != parsed["sites_per_unit"]:
            raise ValueError(
                f"Pool {lab!r}: side_group_homolysis channel {label!r} "
                f"sites_per_unit={parsed['sites_per_unit']:g} contradicts "
                f"the serialized structural match count: "
                f"site_atom_indices={sorted(idxs)} names {len(idxs)} "
                f"site(s). sites_per_unit is CHECKED against the "
                f"serialized selector resolution, never trusted "
                f"(round-72). Fix the artifact.")
        parsed["site_atom_indices"] = [int(i) for i in idxs]
        oob = sorted(i for i in parsed["site_atom_indices"]
                     if i >= n_atoms)
        if oob:
            raise ValueError(
                f"Pool {lab!r}: side_group_homolysis channel {label!r} "
                f"site_atom_indices {oob} are out of range for the "
                f"pool's monomer_adj_list, which names {n_atoms} atom(s) "
                f"(the indices are 0-based positions in its atom order; "
                f"round-75 P1). Fix the artifact.")
        gp = ch["gas_product"]
        if not isinstance(gp, str) or not gp.strip():
            raise ValueError(
                f"Pool {lab!r}: side_group_homolysis channel {label!r} "
                f"gas_product must be a SMILES string, got {gp!r}. Fix "
                f"the artifact.")
        try:
            gmol = Molecule().from_smiles(gp)
        except Exception as e:
            raise ValueError(
                f"Pool {lab!r}: side_group_homolysis channel {label!r} "
                f"gas_product={gp!r} does not parse as SMILES ({e}). Fix "
                f"the artifact.")
        if len(gmol.atoms) != 1 or gmol.get_radical_count() != 1:
            raise ValueError(
                f"Pool {lab!r}: side_group_homolysis channel {label!r} "
                f"gas_product={gp!r} must be a monoatomic mono-radical "
                f"(v1, e.g. '[Br]'); got {len(gmol.atoms)} atom(s), "
                f"radical count {gmol.get_radical_count()}. Fix the "
                f"artifact.")
        parsed["gas_product"] = gp
        gas_sym = gmol.atoms[0].symbol
        mw_x = gmol.get_molecular_weight() * 1000.0
        gw = ch["gas_mw_g_mol"]
        if (isinstance(gw, bool) or not isinstance(gw, (int, float))
                or not math.isfinite(float(gw))
                or abs(float(gw) - mw_x) > 1.0e-6 * mw_x):
            raise ValueError(
                f"Pool {lab!r}: side_group_homolysis channel {label!r} "
                f"gas_mw_g_mol={gw!r} does not pin the gas_product's "
                f"molar mass (M_X({gp}) = {mw_x:g} g/mol). This is the "
                f"value the feature pool's chain_mass_defect_g_mol must "
                f"carry (the NORMATIVE mass formula); a diverging pin "
                f"mints/destroys condensed mass. Fix the artifact.")
        parsed["gas_mw_g_mol"] = float(gw)
        # Double-carry guard, from serialized data alone: EMPTY pairwise
        # intersection between channels' atom sets per gas element
        # (round-72 P1: distinct labels hiding one structural site;
        # round-75 P1-1: OVERLAPPING sets -- subset/superset -- hide it
        # just as well while each satisfies len == sites_per_unit).
        new_set = set(parsed["site_atom_indices"])
        for prev_label, prev_set in seen_sites.get(gas_sym, []):
            shared = sorted(prev_set & new_set)
            if shared:
                raise ValueError(
                    f"Pool {lab!r}: side_group_homolysis channels "
                    f"{prev_label!r} and {label!r} resolve to the SAME "
                    f"{gas_sym} atom set or overlap on atom indices "
                    f"{shared} -- two rate channels claiming one "
                    f"structural site double-carry the loss (rounds "
                    f"72/75 P1: disjointness is EMPTY pairwise "
                    f"intersection, not merely non-identical sets). Fix "
                    f"the artifact.")
        seen_sites.setdefault(gas_sym, []).append((label, new_set))
        gs = ch["gas_species"]
        if not isinstance(gs, str) or not gs.strip():
            raise ValueError(
                f"Pool {lab!r}: side_group_homolysis channel {label!r} "
                f"gas_species={gs!r} must be a non-empty artifact species "
                f"label -- the routing target of the kernel's +R gas "
                f"credit (the ejected X would silently vanish otherwise). "
                f"Fix the artifact.")
        parsed["gas_species"] = gs
        fp = ch["feature_pool"]
        expected_fp = f"{lab}{_SIDE_GROUP_DAUGHTER_INFIX}{san}"
        if not isinstance(fp, str) or fp != expected_fp:
            raise ValueError(
                f"Pool {lab!r}: side_group_homolysis channel {label!r} "
                f"feature_pool={fp!r} must be the ratified X-loss feature "
                f"pool label {expected_fp!r} "
                f"('{{parent}}_sidegrp_{{sanitized channel label}}'). Fix "
                f"the artifact.")
        parsed["feature_pool"] = fp
        out.append(parsed)
    return out


def _check_side_group_homolysis(artifact):
    """Vocabulary/version cross-check + closure guard for the schema-2.7
    side_group_homolysis vocabulary (FR1-K2), the strict mirror of the
    producer's _assert_side_group_serialization_closure (the r68
    producer/consumer mirror-property): BOTH the pool-level block and the
    X-loss feature-pool chain_mass_defect_g_mol field anywhere under a
    below-2.7 stamp are malformed, and every conjunct is checkable from
    serialized data alone. Per carrier: the pool itself must be in
    conventions.configured_pools; per channel, the named feature pool
    must

    * be present in pools[] (its moment slots receive the kernel's
      intact-chain credits),
    * be in conventions.configured_pools and NOT in
      conventions.spawned_pools (feature daughters are eagerly
      solver-configured by design; polymer.pyx _flatten_side_group_state
      hard-errors otherwise),
    * be condensed per the condensed closure (phase_species non-empty and
      fully inside conventions.condensed_species),
    * carry the FR1-K1 spawn provenance (spawn_event_metadata ==
      {'source': 'side_group_homolysis', 'channel': <channel label>}),
    * pin the carrier's monomer_mw_g_mol EXACTLY (intact-chain transfer:
      the X loss is carried by the defect, never the repeat-unit mass),
    * pin chain_mass_defect_g_mol == the channel's gas_mw_g_mol (the
      NORMATIVE mass formula's M_X).

    Reverse closure: a pool entry carrying the side_group_homolysis spawn
    PROVENANCE that NO carrier channel claims is an orphan (the carrier's
    block was lost) -- reject, never adapt. A defect field on any pool is
    shape-validated (finite, > 0), but a defect WITHOUT the provenance is
    legal: Polymer.copy() carries it to downstream daughters, where the
    normative mass formula stays exact with no live channel."""
    pools = [p for p in artifact.get("pools", []) if isinstance(p, dict)]
    carriers = [p for p in pools if "side_group_homolysis" in p]
    defect_pools = [p for p in pools if "chain_mass_defect_g_mol" in p]
    if not carriers and not defect_pools:
        return
    ver = str(artifact.get("schema_version", ""))
    parts = ver.split(".")
    minor = (int(parts[1]) if len(parts) == 2 and parts[0] == "2"
             and parts[1].isdigit() else -1)
    if minor < _SIDE_GROUP_MIN_SCHEMA_MINOR:
        present = []
        if carriers:
            present.append(
                f"side_group_homolysis block on pools "
                f"{[p.get('label') for p in carriers]}")
        if defect_pools:
            present.append(
                f"chain_mass_defect_g_mol on pools "
                f"{[p.get('label') for p in defect_pools]}")
        raise ValueError(
            f"artifact schema_version {ver!r} cannot carry the "
            f"side-group homolysis vocabulary ({'; '.join(present)}): it "
            f"was introduced in schema 2.7, and the emitter stamps 2.7 "
            f"whenever it writes either the block or the X-loss mass "
            f"contract. This artifact is malformed -- regenerate the "
            f"sidecar with a current RMG-Py polymer branch.")
    conv = artifact.get("conventions") or {}
    configured = set(conv.get("configured_pools") or [])
    spawned_raw = conv.get("spawned_pools")
    spawned = (set(spawned_raw)
               if isinstance(spawned_raw, (list, tuple)) else set())
    condensed = set(conv.get("condensed_species") or [])
    by_label = {p.get("label"): p for p in pools}
    claimed = {}
    for carrier in carriers:
        lab = carrier.get("label")
        if lab not in configured:
            raise ValueError(
                f"Pool {lab!r}: carries a side_group_homolysis block but "
                f"is not listed in conventions.configured_pools. The "
                f"loader only builds configured pools, so the kernel "
                f"would be silently dropped (never integrated). The "
                f"producer only serializes the block on solver-configured "
                f"pools; this artifact is hand-edited/corrupted. Fix the "
                f"artifact.")
        if "chain_mass_defect_g_mol" in carrier:
            raise ValueError(
                f"Pool {lab!r}: carries BOTH a side_group_homolysis "
                f"block and a chain_mass_defect_g_mol field "
                f"({carrier.get('chain_mass_defect_g_mol')!r}). v1 "
                f"saturation: the parent pool owns the kernel and X-loss "
                f"feature pools own the defect, NEVER both -- a defected "
                f"carrier claims its chains already lost an X while the "
                f"live kernel debits them for losing another (no "
                f"multi-loss cascade in v1; round-75 P1). Fix the "
                f"artifact.")
        channels = _parse_side_group_homolysis_block(lab, carrier)
        carrier_mw = carrier.get("monomer_mw_g_mol")
        if (isinstance(carrier_mw, bool)
                or not isinstance(carrier_mw, (int, float))
                or not math.isfinite(float(carrier_mw))
                or float(carrier_mw) <= 0.0):
            raise ValueError(
                f"Pool {lab!r}: carries a side_group_homolysis block but "
                f"monomer_mw_g_mol={carrier_mw!r} is not a finite value "
                f"> 0 -- the NORMATIVE mass formula (condensed_mass_g = "
                f"mu1*monomer_mw_g_mol - mu0*chain_mass_defect_g_mol) is "
                f"unevaluatable without it. Fix the artifact.")
        carrier_mw = float(carrier_mw)
        for ch in channels:
            d_label = ch["feature_pool"]
            claimed[d_label] = (lab, ch["label"])
            d_entry = by_label.get(d_label)
            if d_entry is None:
                raise ValueError(
                    f"Pool {lab!r}: side_group_homolysis channel "
                    f"{ch['label']!r} feature pool {d_label!r} is missing "
                    f"from the artifact's pools[] -- the kernel's "
                    f"intact-chain credits would have no moment slots. "
                    f"This artifact is malformed; regenerate the "
                    f"sidecar.")
            if d_label in spawned:
                raise ValueError(
                    f"Pool {lab!r}: side_group_homolysis feature pool "
                    f"{d_label!r} is classified in "
                    f"conventions.spawned_pools. Schema 2.5 defines that "
                    f"list as the configured set's complement, and X-loss "
                    f"feature daughters are eagerly solver-configured by "
                    f"design (polymer.pyx _flatten_side_group_state "
                    f"hard-errors otherwise) -- a spawned-classified "
                    f"feature pool is never built by the loader and the "
                    f"kernel's credits would have no solver home. Fix "
                    f"the artifact.")
            if d_label not in configured:
                raise ValueError(
                    f"Pool {lab!r}: side_group_homolysis feature pool "
                    f"{d_label!r} is not listed in "
                    f"conventions.configured_pools. Feature daughters are "
                    f"eagerly solver-configured by the producer "
                    f"(polymer.pyx _flatten_side_group_state hard-errors "
                    f"otherwise), and the loader only builds configured "
                    f"pools -- an unconfigured feature pool has no solver "
                    f"home for the kernel's credits. Fix the artifact.")
            members = d_entry.get("phase_species") or []
            not_condensed = sorted(m for m in members
                                   if m not in condensed)
            if not members or not_condensed:
                raise ValueError(
                    f"Pool {lab!r}: side_group_homolysis feature pool "
                    f"{d_label!r} is not condensed per the condensed "
                    f"closure (phase_species={members}, not condensed="
                    f"{not_condensed}); its moment slots would be "
                    f"phase-classified GAS (the item-16 mass-balance "
                    f"hazard). Fix the artifact.")
            meta = d_entry.get("spawn_event_metadata")
            if (not isinstance(meta, dict)
                    or meta.get("source") != _SIDE_GROUP_SPAWN_SOURCE
                    or meta.get("channel") != ch["label"]):
                raise ValueError(
                    f"Pool {lab!r}: side_group_homolysis feature pool "
                    f"{d_label!r} lacks the FR1-K1 spawn provenance "
                    f"(expected spawn_event_metadata == {{'source': "
                    f"{_SIDE_GROUP_SPAWN_SOURCE!r}, 'channel': "
                    f"{ch['label']!r}}}, got {meta!r}). Feature daughters "
                    f"are producer-spawned pools, never user pools -- a "
                    f"row claiming other provenance is hand-edited. Fix "
                    f"the artifact.")
            d_mw = d_entry.get("monomer_mw_g_mol")
            if (isinstance(d_mw, bool)
                    or not isinstance(d_mw, (int, float))
                    or not math.isfinite(float(d_mw))
                    or abs(float(d_mw) - carrier_mw) >
                    1.0e-9 * max(abs(carrier_mw), 1.0)):
                raise ValueError(
                    f"Pool {lab!r}: side_group_homolysis feature pool "
                    f"{d_label!r} monomer_mw_g_mol={d_mw!r} does not pin "
                    f"the carrier's ({carrier_mw!r}). The chain transfers "
                    f"INTACT -- the X loss is carried by "
                    f"chain_mass_defect_g_mol, never by the repeat-unit "
                    f"mass; a diverging monomer_mw fabricates/destroys "
                    f"condensed mass. Fix the artifact.")
            defect = d_entry.get("chain_mass_defect_g_mol")
            mw_x = ch["gas_mw_g_mol"]
            if (isinstance(defect, bool)
                    or not isinstance(defect, (int, float))
                    or not math.isfinite(float(defect))
                    or not float(defect) > 0.0
                    or abs(float(defect) - mw_x) > 1.0e-6 * mw_x):
                raise ValueError(
                    f"Pool {lab!r}: side_group_homolysis feature pool "
                    f"{d_label!r} has chain_mass_defect_g_mol={defect!r} "
                    f"but channel {ch['label']!r} ejects "
                    f"gas_product={ch['gas_product']!r} with "
                    f"M_X={mw_x:g} g/mol. The NORMATIVE mass formula "
                    f"(condensed_mass_g = mu1*monomer_mw_g_mol - "
                    f"mu0*chain_mass_defect_g_mol) requires the defect to "
                    f"pin M_X exactly -- anything else mints or destroys "
                    f"condensed mass while gas X appears. Fix the "
                    f"artifact.")
    for p in pools:
        lbl = p.get("label")
        # Any serialized defect must at least be a well-formed exact
        # per-chain mass loss (finite, > 0) -- the emitter only writes
        # positive defects, and the reconstructed config feeds it straight
        # into the normative condensed_mass_g accessor.
        if "chain_mass_defect_g_mol" in p:
            defect = p.get("chain_mass_defect_g_mol")
            if (isinstance(defect, bool)
                    or not isinstance(defect, (int, float))
                    or not math.isfinite(float(defect))
                    or not float(defect) > 0.0):
                raise ValueError(
                    f"Pool {lbl!r}: chain_mass_defect_g_mol={defect!r} "
                    f"must be a finite value > 0 (the exact per-chain "
                    f"X-loss mass; the emitter never writes any other "
                    f"shape). Fix the artifact.")
        # Reverse closure (orphan guard), keyed on the spawn PROVENANCE
        # only: a pool stamped as a producer-spawned X-loss feature
        # daughter that NO channel claims is corrupted (its label is a
        # deterministic function of (parent, channel), so unclaimed
        # provenance means the carrier's block was lost). A
        # defect-carrying pool WITHOUT the provenance is deliberately
        # legal: Polymer.copy() carries chain_mass_defect_g_mol to
        # downstream daughters (e.g. an S2 _mod child of a feature pool
        # -- its chains each lost one X too), and the NORMATIVE mass
        # formula stays exact there with no live channel.
        meta = p.get("spawn_event_metadata")
        if (isinstance(meta, dict) and meta.get("source") ==
                _SIDE_GROUP_SPAWN_SOURCE and lbl not in claimed):
            raise ValueError(
                f"Pool {lbl!r}: carries the side_group_homolysis spawn "
                f"provenance (a producer-spawned X-loss feature daughter) "
                f"but NO side_group_homolysis channel in the artifact "
                f"claims it as its feature pool. The feature-pool label "
                f"is a deterministic function of (parent, channel), so an "
                f"unclaimed provenance pin means the carrier's block was "
                f"lost -- the kernel's flux into this pool cannot be "
                f"reproduced. Fix the artifact.")


def _validate_refused_entry(e):
    """Shape-guard one reactions[] entry's refused vocabulary (format doc
    §12; reject, never adapt): the emitter writes ``refused: true`` plus a
    non-empty ``refused_reason`` on POOL-MAPPED rows, or NOTHING (absent,
    not false). Returns True iff the row is a valid refused row."""
    if "refused" not in e and "refused_reason" not in e:
        return False
    eid = e.get("id")
    if e.get("refused") is not True:
        raise ValueError(
            f"reactions[] entry {eid!r} carries malformed refused vocabulary "
            f"(refused={e.get('refused')!r}): the emitter writes refused: "
            f"true + refused_reason on refused rows and NOTHING otherwise "
            f"(absent, not false). A refused_reason without refused: true "
            f"is equally malformed. Fix the artifact.")
    reason = e.get("refused_reason")
    if not isinstance(reason, str) or not reason:
        raise ValueError(
            f"reactions[] entry {eid!r} has refused: true but no valid "
            f"refused_reason ({reason!r}): the emitter always writes the "
            f"stamp-time census reason (a non-empty string, e.g. "
            f"'conduit-deferred' / 'qssa-invalid'). Fix the artifact.")
    if reason not in REFUSED_REASONS:
        raise ValueError(
            f"reactions[] entry {eid!r} has unknown refused_reason "
            f"{reason!r}: the schema-2.4 vocabulary is CLOSED to "
            f"{sorted(REFUSED_REASONS)} (format doc §12). The accumulating "
            f"class is reconstructed from this reason, so an unknown value "
            f"would silently change solver semantics. Fix the artifact.")
    if not (e.get("proxy_reactants") or e.get("proxy_products")):
        raise ValueError(
            f"reactions[] entry {eid!r} has refused: true but no pool-mapped "
            f"participant (empty proxy_reactants AND proxy_products): the "
            f"refused marker is legal ONLY on pool-mapped rows -- honoring "
            f"it here would silently zero ordinary gas chemistry. Fix the "
            f"artifact.")
    return True


# The artifact-level recipe-revision tokens that carry the explicit-DP
# handshake algebra (the 3-way channel-family split, rmgpy/polymer.py
# POLYMER_RATE_RECIPE_REVISION_EXPLICIT_DP*). The pairing with the pool
# vocabulary is exact in BOTH directions (P2 gates, stage-B adversarial
# review): a token without any explicit_dp block claims handshake algebra
# the pools do not carry, and a block without the token claims rates from a
# recipe with no handshake algebra. The producer can emit neither shape;
# hand-edited artifacts must fail loud.
_EXPLICIT_DP_RECIPE_REVISIONS = frozenset({
    "2026-07-04-explicit-dp-monomer-gas",
    "2026-07-04-explicit-dp-qssa-monomer-gas",
    "2026-07-04-explicit-dp-weaklink-u-monomer-gas",
})


def _check_explicit_dp_recipe_revision(artifact):
    """Enforce the explicit-DP token/vocabulary pairing in both directions
    (see _EXPLICIT_DP_RECIPE_REVISIONS)."""
    rev = (artifact.get("conventions") or {}).get("recipe_revision")
    carriers = [p.get("label") for p in artifact.get("pools", [])
                if isinstance(p, dict) and "explicit_dp" in p]
    if rev in _EXPLICIT_DP_RECIPE_REVISIONS and not carriers:
        raise ValueError(
            f"artifact recipe_revision {rev!r} declares the explicit-dp "
            f"handshake algebra but NO pool carries an explicit_dp block. "
            f"The emitter stamps an explicit-dp token exactly when a block "
            f"is present; this artifact is hand-edited/corrupted -- "
            f"regenerate the sidecar.")
    if carriers and rev not in _EXPLICIT_DP_RECIPE_REVISIONS:
        raise ValueError(
            f"pools {carriers} carry explicit_dp blocks but the artifact "
            f"recipe_revision is {rev!r} (not an explicit-dp token: "
            f"{sorted(_EXPLICIT_DP_RECIPE_REVISIONS)}). The rates would "
            f"claim a recipe with no handshake algebra; reject, never "
            f"integrate permissively -- regenerate the sidecar.")


def _parse_explicit_dp_block(lab, pool_entry):
    """Parse + validate a pool entry's pool-level explicit_dp block (schema
    2.3). Returns ``(dp_to_label, initial_moles)`` — ``{int dp: str label}``
    and ``{int dp: float moles}`` — or ``None`` when the block is absent.

    Boundary rules mirror the QSSA channel parser: closed key vocabulary,
    boolean ``enabled`` (present-disabled REJECTED — the emitter never
    writes disabled blocks), the normative ``recipe`` and the block
    ``recipe_revision`` pinned by exact match (reject, never adapt), and
    the handshake invariants the emitter guarantees (exactly one species
    entry in v1, at DP == handshake_target_dp == the pool's cutoff;
    initial_moles keyed identically, finite and >= 0)."""
    raw = pool_entry.get("explicit_dp")
    if raw is None:
        return None
    if not isinstance(raw, dict):
        raise ValueError(
            f"Pool {lab!r}: explicit_dp must be a dict, got "
            f"{type(raw).__name__}. Fix the artifact.")
    unknown = sorted(set(raw) - _EXPLICIT_DP_BLOCK_KEYS)
    if unknown:
        raise ValueError(
            f"Pool {lab!r}: explicit_dp has unknown key(s) {unknown}; "
            f"allowed keys are {sorted(_EXPLICIT_DP_BLOCK_KEYS)}. Fix the "
            f"artifact (unknown vocabulary is never dropped permissively).")
    missing = sorted(_EXPLICIT_DP_BLOCK_KEYS - set(raw))
    if missing:
        raise ValueError(
            f"Pool {lab!r}: explicit_dp is missing key(s) {missing}. The "
            f"emitter always writes the full block -- regenerate the "
            f"sidecar.")
    enabled = raw["enabled"]
    if not isinstance(enabled, bool):
        raise ValueError(
            f"Pool {lab!r}: explicit_dp must carry a boolean 'enabled' "
            f"field, got {enabled!r}. Fix the artifact.")
    if not enabled:
        raise ValueError(
            f"Pool {lab!r}: explicit_dp carries enabled=false. The emitter "
            f"never writes disabled blocks: a disabled handshake must be "
            f"absent from the sidecar, not present-disabled. Fix the "
            f"artifact (remove the block).")
    if raw["recipe_revision"] != _EXPLICIT_DP_PINNED_REVISION:
        raise ValueError(
            f"Pool {lab!r}: explicit_dp recipe_revision must equal "
            f"{_EXPLICIT_DP_PINNED_REVISION!r} exactly, got "
            f"{raw['recipe_revision']!r}. An artifact claiming a different "
            f"handshake recipe must be fixed at the source; this loader "
            f"validates, never adapts.")
    recipe = raw["recipe"]
    if not isinstance(recipe, dict):
        raise ValueError(
            f"Pool {lab!r}: explicit_dp must carry the normative 'recipe' "
            f"dict (schema 2.3), got {recipe!r} -- regenerate the sidecar.")
    unknown_recipe = sorted(set(recipe) - set(_EXPLICIT_DP_PINNED_RECIPE))
    if unknown_recipe:
        raise ValueError(
            f"Pool {lab!r}: explicit_dp recipe has unknown key(s) "
            f"{unknown_recipe}; allowed keys are "
            f"{sorted(_EXPLICIT_DP_PINNED_RECIPE)}. Fix the artifact.")
    for key, pinned in _EXPLICIT_DP_PINNED_RECIPE.items():
        if key not in recipe or recipe[key] != pinned:
            raise ValueError(
                f"Pool {lab!r}: explicit_dp recipe[{key!r}] must equal the "
                f"pinned normative recipe exactly; got {recipe.get(key)!r}, "
                f"expected {pinned!r}. An artifact claiming a different "
                f"handshake algebra must be fixed at the source; this "
                f"loader validates, never adapts.")
    species = raw["species"]
    moles = raw["initial_moles"]
    if not isinstance(species, dict) or not species:
        raise ValueError(
            f"Pool {lab!r}: explicit_dp species must be a non-empty "
            f"{{dp: label}} dict, got {species!r}. Fix the artifact.")
    if not isinstance(moles, dict) or set(moles) != set(species):
        raise ValueError(
            f"Pool {lab!r}: explicit_dp initial_moles must be keyed "
            f"identically to species (got keys {sorted(moles) if isinstance(moles, dict) else moles!r} "
            f"vs {sorted(species)}). Fix the artifact.")
    target = raw["handshake_target_dp"]
    cutoff = pool_entry.get("cutoff")
    if not isinstance(target, int) or isinstance(target, bool) or (
            cutoff is not None and int(cutoff) != target):
        raise ValueError(
            f"Pool {lab!r}: explicit_dp handshake_target_dp must be the "
            f"pool's cutoff as an integer (cutoff={cutoff!r}), got "
            f"{target!r}. Fix the artifact.")
    if sorted(species) != [str(target)]:
        raise ValueError(
            f"Pool {lab!r}: explicit_dp species must carry exactly ONE "
            f"entry at DP == handshake_target_dp ({target}) in v1, got keys "
            f"{sorted(species)}. Fix the artifact.")
    dp_to_label = {}
    initial_moles = {}
    for dp_key, label in species.items():
        if not isinstance(label, str) or not label:
            raise ValueError(
                f"Pool {lab!r}: explicit_dp species[{dp_key!r}] must be a "
                f"chem.yaml species-label string, got {label!r}. Fix the "
                f"artifact.")
        m = moles[dp_key]
        if not isinstance(m, (int, float)) or isinstance(m, bool) or \
                not math.isfinite(float(m)) or float(m) < 0.0:
            raise ValueError(
                f"Pool {lab!r}: explicit_dp initial_moles[{dp_key!r}] must "
                f"be a finite number >= 0, got {m!r}. Fix the artifact.")
        dp_to_label[int(dp_key)] = label
        initial_moles[int(dp_key)] = float(m)
    return dp_to_label, initial_moles


def _check_recipe_revision_monomer_phase(artifact):
    """Revision-gate the routed-monomer phase semantics (P1-B).

    Only artifacts that route released monomer (any configured pool with
    a non-null ``monomer_routing``) carry the semantics at all; artifacts
    without routing pass untouched, whatever their revision.

    * NEW (monomer-gas) recipe_revision: the routed monomer must NOT
      appear in ``conventions.condensed_species`` or the pool's
      ``phase_species`` -- such an artifact is internally contradictory
      (gas contract stamped, condensed membership listed) and is
      REJECTED rather than mis-phased.
    * OLD (pre-monomer-gas) or unknown recipe_revision: HARD REFUSAL
      (decision (a)). Legacy acceptance would require re-condensing the
      routed monomer on the live solver path -- the exact
      reference-state conflation revision 2026-07-03-monomer-gas
      removed, and the solver now validates the release target GAS
      (validate_configuration raises on a condensed monomer index). The
      gate is revision-keyed, NOT semantics-sniffed, so re-freezing an
      old artifact with gas-looking lists cannot launder it past this
      check."""
    conv = artifact.get("conventions") or {}
    rev = conv.get("recipe_revision")
    configured = set(conv.get("configured_pools") or [])
    condensed = set(conv.get("condensed_species") or [])
    # Non-string monomer_routing values are NOT gated here: they fall
    # through to the existing per-pool type guard in
    # build_system_from_artifact, whose message names the pool and type.
    routed = [(p.get("label"), p.get("monomer_routing"), p)
              for p in artifact.get("pools", [])
              if isinstance(p, dict) and p.get("label") in configured
              and isinstance(p.get("monomer_routing"), str)
              and p.get("monomer_routing")]
    if not routed:
        return  # no routed monomer: no monomer-phase semantics in play
    if rev in _MONOMER_GAS_RECIPE_REVISIONS:
        for lab, routing, p in routed:
            phase_species = p.get("phase_species") or []
            if routing in condensed or routing in phase_species:
                raise ValueError(
                    f"Pool {lab!r}: recipe_revision {rev!r} declares the "
                    f"monomer-gas contract (the routed monomer is a "
                    f"GAS-phase species), but the monomer_routing target "
                    f"{routing!r} is still listed in "
                    f"conventions.condensed_species / pools[].phase_species. "
                    f"The artifact is internally contradictory; regenerate "
                    f"it with current code (monomer-gas contract, "
                    f"docs/polymer_moments_format.md revision note "
                    f"2026-07-03-monomer-gas).")
        return
    kind = ("predates" if rev in _PRE_MONOMER_GAS_RECIPE_REVISIONS
            else "is unknown to")
    routed_desc = ", ".join(f"pool {lab!r} -> {routing!r}"
                            for lab, routing, _ in routed)
    raise ValueError(
        f"artifact recipe_revision {rev!r} {kind} the monomer-gas contract "
        f"(2026-07-03: the unzip/QSSA release target is a GAS-phase "
        f"species) and the artifact routes released monomer "
        f"({routed_desc}). This loader implements ONLY the gas-monomer "
        f"semantics: legacy acceptance would re-condense the routed "
        f"monomer -- the exact reference-state conflation the revision "
        f"removed -- and the live solver refuses a condensed release "
        f"target. Regenerate the artifact with current code (monomer-gas "
        f"contract; recipe_revision one of "
        f"{sorted(_MONOMER_GAS_RECIPE_REVISIONS)}).")


def _parse_radical_qssa_channel(lab, pool_entry):
    """Parse + validate a pool entry's channels.radical_qssa_unzip block.

    Returns the validated channel config dict (the exact shape
    PolymerPoolConfig.radical_qssa_unzip expects, so the M1 flattening picks
    it up), or None when the block is absent. A present block with
    enabled == false is REJECTED, not skipped: the serializer never emits
    disabled blocks (a disabled channel is absent), so present-disabled can
    only mean a hand-edited/corrupted artifact whose constants nothing
    would validate.

    Validation is layered: the sidecar-envelope rules (boolean ``enabled``,
    unknown keys, per-block ``units`` pinned to the emitter's convention,
    the normative ``recipe`` pinned by exact match) live here, at the
    artifact boundary; the FIELD rules (finite triplets, A > 0, Ea >= 0,
    (0,1] fractions, pinned ``basis``) are delegated to
    validate_radical_qssa_unzip -- the shared single source of truth, not
    duplicated. Raises ValueError naming the pool on any violation.
    """
    raw = (pool_entry.get("channels") or {}).get("radical_qssa_unzip")
    if raw is None:
        return None
    if not isinstance(raw, dict):
        raise ValueError(
            f"Pool {lab!r}: channels.radical_qssa_unzip must be a dict, got "
            f"{type(raw).__name__}. Fix the artifact.")
    unknown = sorted(set(raw) - _QSSA_SIDECAR_KEYS)
    if unknown:
        raise ValueError(
            f"Pool {lab!r}: channels.radical_qssa_unzip has unknown key(s) "
            f"{unknown}; allowed keys are {sorted(_QSSA_SIDECAR_KEYS)}. "
            f"Fix the artifact (unknown sub-vocabulary is never dropped "
            f"permissively).")
    enabled = raw.get("enabled")
    if not isinstance(enabled, bool):
        raise ValueError(
            f"Pool {lab!r}: channels.radical_qssa_unzip must carry a boolean "
            f"'enabled' field, got {enabled!r}. Fix the artifact.")
    if not enabled:
        raise ValueError(
            f"Pool {lab!r}: channels.radical_qssa_unzip carries "
            f"enabled=false. The emitter never writes disabled blocks: a "
            f"disabled channel must be absent from the sidecar, not "
            f"present-disabled (whose constants nothing validates). Fix the "
            f"artifact (remove the block).")

    # The weak-link allyl/U-state sub-vocabulary (schema 2.2) selects the
    # weak-link recipe pin; the group's all-or-nothing rule and the summed-
    # termination exclusion are the shared validator's job below.
    weaklink = bool(_QSSA_WEAKLINK_KEYS & set(raw))
    pinned_recipe = (_QSSA_PINNED_RECIPE_WEAKLINK if weaklink
                     else _QSSA_PINNED_RECIPE)

    # Normative recipe pin (schema 2.1/2.2): every field must match the
    # loader's own copy EXACTLY -- reject, never adapt (units-pin idiom). A
    # QSSA block without a recipe is malformed: the emitter always writes
    # one, and a weak-link block carrying only the LEGACY recipe (or vice
    # versa) claims an algebra this loader does not implement for it.
    recipe = raw.get("recipe")
    if not isinstance(recipe, dict):
        raise ValueError(
            f"Pool {lab!r}: channels.radical_qssa_unzip must carry the "
            f"normative 'recipe' block (schema 2.1), got "
            f"{recipe!r}. A QSSA block without its machine-readable recipe "
            f"is malformed -- regenerate the sidecar.")
    unknown_recipe = sorted(set(recipe) - set(pinned_recipe))
    if unknown_recipe:
        raise ValueError(
            f"Pool {lab!r}: channels.radical_qssa_unzip recipe has unknown "
            f"key(s) {unknown_recipe}; allowed keys are "
            f"{sorted(pinned_recipe)}. Fix the artifact (unknown "
            f"recipe vocabulary is never dropped permissively).")
    for key, pinned in pinned_recipe.items():
        if key not in recipe or recipe[key] != pinned:
            raise ValueError(
                f"Pool {lab!r}: channels.radical_qssa_unzip recipe[{key!r}] "
                f"must equal the pinned normative recipe exactly; got "
                f"{recipe.get(key)!r}, expected {pinned!r}. An artifact "
                f"claiming a different rate algebra must be fixed at the "
                f"source; this loader validates, never adapts.")

    channel = {}
    for block_name, pinned_a in _QSSA_PINNED_A_UNITS.items():
        trip = raw.get(block_name)
        if trip is None:
            # transfer: null is the valid channel-absent shape (and a
            # weak-link block carries termination: null -- the split
            # triplets replace it, so the key is simply not forwarded); a
            # missing mandatory block is diagnosed by the shared validator
            # below.
            if block_name == "transfer":
                channel["transfer"] = None
            continue
        if not isinstance(trip, dict):
            # let the shared validator produce its canonical triplet error
            channel[block_name] = trip
            continue
        units = trip.get("units")
        expected = {"A": pinned_a, "Ea": "J/mol"}
        if units != expected:
            note = ""
            if block_name == "transfer":
                note = (" ktr is PSEUDO-first-order: premultiply a "
                        "bimolecular literature k_tr by the substrate "
                        "concentration [mol/m^3] upstream, do not relabel "
                        "the units.")
            raise ValueError(
                f"Pool {lab!r}: channels.radical_qssa_unzip {block_name} "
                f"units must be exactly {expected} (pinned convention), got "
                f"{units!r}. A sidecar claiming different units must be "
                f"fixed at the source; this loader never converts.{note}")
        channel[block_name] = {k: v for k, v in trip.items() if k != "units"}
    if "unsaturated_tail_ends_initial" in raw:
        # U0 is STATE, serialized {value, units-note}. The shape and the
        # units note are pinned here at the boundary; the VALUE rules
        # (number, finite, >= 0) are the shared validator's, and the
        # capacity census (U0 <= 2*mu0_tail) is the SOLVER's at init --
        # the loader passes the value through.
        u0 = raw["unsaturated_tail_ends_initial"]
        if not isinstance(u0, dict) or set(u0) != {"value", "units"}:
            raise ValueError(
                f"Pool {lab!r}: channels.radical_qssa_unzip "
                f"unsaturated_tail_ends_initial must be a "
                f"{{value, units}} dict (schema 2.2), got {u0!r}. Fix the "
                f"artifact.")
        if u0["units"] != _QSSA_U0_PINNED_UNITS:
            raise ValueError(
                f"Pool {lab!r}: channels.radical_qssa_unzip "
                f"unsaturated_tail_ends_initial units must be exactly "
                f"{_QSSA_U0_PINNED_UNITS!r} (pinned convention: same "
                f"amount basis as mu0), got {u0['units']!r}. A sidecar "
                f"claiming a different basis must be fixed at the source; "
                f"this loader never converts.")
        channel["unsaturated_tail_ends_initial"] = u0["value"]
    for name in ("efficiency", "monomer_yield", "basis"):
        if name in raw:
            channel[name] = raw[name]
    # shared single source of truth for field rules (basis pin included)
    return validate_radical_qssa_unzip(lab, channel)


def build_system_from_artifact(artifact, species, reactions,
                               T0, P, V_poly, initial_moles,
                               mass_transfer_spec, initial_moments=None):
    """Assemble the HybridPolymerSystem oracle from consumer-world inputs.

    Returns (system, core_species, all_reactions) — all_reactions includes
    the cantera-null reconstructions and is needed by run_segments for the
    per-segment generate_rate_coefficients call (HybridPolymerSystem is a
    cdef class; do not hang extra attributes on it). The system is fully
    initialized at T0 (initialize_model runs through initialize_solver,
    rmgpy/solver/polymer.pyx:601-610)."""
    # Envelope gate first (schema minor must be one this loader implements),
    # then the QSSA-vocabulary/version cross-check: a 2.0 artifact carrying a
    # radical_qssa_unzip block (or a 2.1 artifact carrying the weak-link
    # sub-vocabulary) is malformed regardless of pool configuration.
    _check_schema_version_known(artifact)
    _check_qssa_schema_version(artifact)
    # Explicit-DP vocabulary/version cross-check (schema 2.3): a 2.0-2.2
    # artifact carrying a pool-level explicit_dp block is malformed.
    _check_explicit_dp_schema_version(artifact)
    # Refused-row vocabulary/version cross-check (schema 2.4): a 2.0-2.3
    # artifact carrying a refused marker on any reactions[] row is
    # malformed (the emitter stamps 2.4 whenever it writes one).
    _check_refused_schema_version(artifact)
    # Spawned-pool closure vocabulary/version cross-check (schema 2.5): a
    # 2.0-2.4 artifact carrying conventions.spawned_pools is malformed, and
    # so is any overlap with configured_pools (the closure is the
    # configured set's complement by construction).
    _check_spawned_pools_schema_version(artifact)
    # Homolysis-initiation vocabulary/version cross-check + daughter
    # closure (schema 2.6, r68-tightened): a 2.0-2.5 artifact carrying the
    # pool-level block is malformed, and a 2.6 block whose carrier is
    # unconfigured or whose end-radical daughters are missing from
    # pools[], not solver-configured (or spawned-classified),
    # un-condensed, or provenance-stripped is rejected before any pool
    # config is built.
    _check_homolysis_initiation(artifact)
    # Side-group homolysis vocabulary/version cross-check + feature-pool
    # closure (schema 2.7, FR1-K2): a 2.0-2.6 artifact carrying the
    # pool-level block OR the X-loss chain_mass_defect_g_mol field is
    # malformed; a 2.7 block whose carrier is unconfigured or whose
    # feature pools are missing from pools[], not solver-configured (or
    # spawned-classified), un-condensed, provenance-stripped, or whose
    # monomer_mw/defect pins break the NORMATIVE mass formula is rejected
    # before any pool config is built -- as is any orphan X-loss pool no
    # channel claims.
    _check_side_group_homolysis(artifact)
    # ... and the explicit-dp token/vocabulary pairing is exact in both
    # directions (P2 gates): token without block, block without token --
    # both hand-edited shapes fail loud.
    _check_explicit_dp_recipe_revision(artifact)
    # Monomer-phase semantics are revision-gated (P1-B): reject
    # contradictory NEW-revision artifacts, hard-refuse pre-monomer-gas
    # ones that route monomer (see _check_recipe_revision_monomer_phase).
    _check_recipe_revision_monomer_phase(artifact)

    core = list(species)
    idx = {s.label: i for i, s in enumerate(core)}
    all_reactions = _restamp_and_extend(artifact, core, reactions)

    conv = artifact["conventions"]
    condensed = set(conv["condensed_species"])
    mask = np.array([s.label not in condensed for s in core], dtype=bool)

    pools = []
    moments0 = {}
    initial_explicit_by_pool = {}
    for p in artifact["pools"]:
        lab = p["label"]
        if lab not in conv["configured_pools"]:
            continue
        mu_idx = tuple(idx[f"{lab}_mu{k}"] for k in range(3))
        routing = p.get("monomer_routing")
        k_unzip = float(p["channels"]["unzip"]["A"])
        k_scission = float(p["channels"]["scission"]["A"])
        # HARD ERROR: a non-finite unzip/scission A is not a valid rate
        # constant. NaN passes BOTH the `< 0` and `> 0` gates as False, so a
        # hand-edited artifact could make the channel SILENTLY INERT (a
        # laundered no-op); inf poisons the residual. Mirror the negative-A
        # guard below.
        for _rate_name, _channel_name, _rate_val in (
                ("k_unzip", "unzip", k_unzip),
                ("k_scission", "scission", k_scission)):
            if not math.isfinite(_rate_val):
                raise ValueError(
                    f"Pool {lab!r}: artifact declares {_channel_name} "
                    f"A={_rate_val!r} -- a non-finite {_rate_name} is not a "
                    f"valid rate constant (NaN would silently disable the "
                    f"channel; inf would poison the residual). Fix the "
                    f"artifact's {_channel_name} channel A (finite, >= 0).")
        # HARD ERROR: a negative unzip A is not a valid rate constant. Every
        # downstream k_unzip consumer is gated on k_unzip > 0, so a negative
        # value would silently become an inert channel instead of failing.
        if k_unzip < 0.0:
            raise ValueError(
                f"Pool {lab!r}: artifact declares unzip A={k_unzip:g} -- a negative "
                f"k_unzip is not a valid rate constant. Fix the artifact's unzip "
                f"channel A (set it >= 0).")
        # Routing schema: monomer_routing is a species-label string (see
        # polymer.py _serialize_pool_for_sidecar). A hand-edited artifact with
        # any other type used to die with a raw TypeError from idx[routing]
        # (e.g. unhashable dict); an unknown label died with a raw KeyError.
        # Both must be actionable ValueErrors naming the pool.
        if routing is not None and not isinstance(routing, str):
            raise ValueError(
                f"Pool {lab!r}: monomer_routing must be a species-label string "
                f"(the released monomer's label in the deck), got "
                f"{type(routing).__name__}: {routing!r}. Fix the artifact's "
                f"monomer_routing.")
        # HARD ERROR: same guard as generation-side PolymerPool.to_config. A
        # config with k_unzip > 0 and monomer_poly_index=None makes the solver
        # drain the condensed moments with no monomer emission -- silently
        # un-conserved mass. A post-guard generation run can never write such
        # an artifact; refuse legacy/hand-edited ones at assembly time.
        if k_unzip > 0.0 and not routing:
            raise ValueError(
                f"Pool {lab!r}: artifact declares k_unzip={k_unzip:g} > 0 but no "
                f"monomer_routing target. The unzip channel would drain the condensed "
                f"moments with no released-monomer emission, leaving mass silently "
                f"un-conserved. Fix the artifact's monomer_routing or set the unzip "
                f"channel A to 0.")
        # Radical QSSA unzip channel (milestone 3): parse + validate the
        # sidecar block (units pin, shared field validator) and enforce the
        # cross-invariants at the artifact boundary, mirroring the k_unzip
        # guards above and the generation-side layers.
        qssa_cfg = _parse_radical_qssa_channel(lab, p)
        if qssa_cfg is not None:
            if k_unzip > 0.0:
                raise ValueError(
                    f"Pool {lab!r}: artifact declares an enabled "
                    f"radical_qssa_unzip channel AND unzip A={k_unzip:g} > 0 "
                    f"(k_unzip). These are two representations of the SAME "
                    f"chain-end depropagation channel and are mutually "
                    f"exclusive on a pool (both would double-count the unzip "
                    f"flux). Fix the artifact: zero the unzip channel or "
                    f"disable radical_qssa_unzip.")
            if not routing:
                raise ValueError(
                    f"Pool {lab!r}: artifact declares an enabled "
                    f"radical_qssa_unzip channel but no monomer_routing "
                    f"target. The QSSA channel releases depropagated monomer "
                    f"through the pool's monomer routing; without a target "
                    f"the condensed mass would leave silently un-conserved. "
                    f"Fix the artifact's monomer_routing.")
        monomer_idx = None
        if routing:
            monomer_idx = idx.get(routing)
            if monomer_idx is None:
                raise ValueError(
                    f"Pool {lab!r}: monomer_routing target {routing!r} is not in "
                    f"the deck's species list; cannot wire unzip-to-monomer "
                    f"release. Fix the artifact's monomer_routing to name a "
                    f"species present in chem.yaml.")
        # Explicit-DP handshake block (schema 2.3): resolve the species
        # labels against chem.yaml and wire the pool's explicit map + the
        # t=0 loadings, so the oracle's hybrid handshake deposits
        # boundary-crossing chains exactly as the generation run did. A
        # label the deck cannot resolve is a hard error (the handshake
        # would silently strand its flux otherwise).
        explicit_map = {}
        parsed_explicit = _parse_explicit_dp_block(lab, p)
        if parsed_explicit is not None:
            dp_to_label, dp_moles = parsed_explicit
            for dp, spc_label in dp_to_label.items():
                spc_idx = idx.get(spc_label)
                if spc_idx is None:
                    raise ValueError(
                        f"Pool {lab!r}: explicit_dp species {spc_label!r} "
                        f"(DP {dp}) is not in the deck's species list; "
                        f"cannot wire the hybrid handshake target. Fix the "
                        f"artifact's explicit_dp species to name a species "
                        f"present in chem.yaml.")
                explicit_map[dp] = spc_idx
            if any(v > 0.0 for v in dp_moles.values()):
                initial_explicit_by_pool[lab] = dict(dp_moles)
        # Radical-homolysis initiation kernel (schema 2.6): parse +
        # validate the pool-level block (pinned kernel/units/recipe;
        # artifact-level closure already enforced by
        # _check_homolysis_initiation above) and wire the triplet into the
        # reconstructed config, so the rebuilt oracle's RHS carries the
        # generating solver's kernel flux instead of silently integrating
        # a melt with no initiation (the flat/false-trajectory failure
        # class). The solver's own validate_configuration +
        # _flatten_homolysis_state re-enforce the daughter-configured
        # invariant and the k_scission/QSSA/k_unzip mutual exclusions
        # natively on initialize_model.
        khom_cfg = _parse_homolysis_initiation_block(lab, p)
        # Side-group homolysis kernel (schema 2.7, FR1-K2): parse +
        # validate the pool-level block (pinned kernel/units/recipe;
        # data-alone selector closure; artifact-level feature-pool closure
        # already enforced by _check_side_group_homolysis above), enforce
        # the generation-side mutual exclusions at the artifact boundary,
        # resolve each channel's gas_species routing against chem.yaml,
        # and wire the channel list + gas indices into the reconstructed
        # config -- so the rebuilt oracle's RHS carries the generating
        # solver's side-group flux and its EXACT mass contract
        # (chain_mass_defect_g_mol, wired below) re-enforces natively at
        # initialize_model (_flatten_side_group_state).
        sgh_channels = None
        sgh_gas_indices = None
        sgh_parsed = _parse_side_group_homolysis_block(lab, p)
        if sgh_parsed is not None:
            if k_unzip > 0.0:
                raise ValueError(
                    f"Pool {lab!r}: artifact declares an enabled "
                    f"side_group_homolysis block AND unzip A="
                    f"{k_unzip:g} > 0 (k_unzip). Legacy k_unzip is a "
                    f"phenomenological closed-chain monomer-loss channel, "
                    f"while side_group_homolysis creates radical-defect "
                    f"feature pools feeding explicit degradation "
                    f"chemistry; the two are mutually exclusive on a pool "
                    f"(generation-side layers enforce the same). Fix the "
                    f"artifact.")
            if qssa_cfg is not None:
                raise ValueError(
                    f"Pool {lab!r}: artifact declares an enabled "
                    f"side_group_homolysis block AND an enabled "
                    f"radical_qssa_unzip channel. Two lumped initiation "
                    f"carriers on one pool are mutually exclusive "
                    f"(generation-side layers enforce the same). Fix the "
                    f"artifact.")
            sgh_channels = [
                {k: ch[k] for k in ("label", "A", "n", "Ea",
                                    "site_selector", "sites_per_unit",
                                    "gas_product")}
                for ch in sgh_parsed]
            gas_idx = []
            for ch in sgh_parsed:
                g_idx = idx.get(ch["gas_species"])
                if g_idx is None:
                    raise ValueError(
                        f"Pool {lab!r}: side_group_homolysis channel "
                        f"{ch['label']!r} gas_species "
                        f"{ch['gas_species']!r} is not in the deck's "
                        f"species list; cannot wire the kernel's gas "
                        f"X-radical release (the ejected X would silently "
                        f"vanish -- un-conserved mass). Fix the artifact "
                        f"to name a species present in chem.yaml.")
                gas_idx.append(g_idx)
            sgh_gas_indices = tuple(gas_idx)
        pools.append(PolymerPoolConfig(
            label=lab, xs=int(p.get("cutoff") or 0),
            explicit_dp_to_species_index=explicit_map,
            mu_indices=mu_idx,
            monomer_poly_index=monomer_idx,
            # The tripwire's ONE chain-scale window (and the spawn-gate
            # snapshot) consume this; without it the consumer-world window
            # collapses to the bare slack and the tag-branch class drifts
            # from generation world.
            monomer_mw_g_mol=float(p.get("monomer_mw_g_mol") or 0.0),
            k_scission=k_scission,
            k_unzip=k_unzip,
            radical_qssa_unzip=qssa_cfg,
            k_homolysis=khom_cfg,
            side_group_homolysis=sgh_channels,
            side_group_gas_indices=sgh_gas_indices,
            # X-loss feature-pool EXACT mass contract (schema 2.7): the
            # solver's _flatten_side_group_state re-checks defect == M_X
            # per channel destination at initialize_model, and every
            # defect-aware mass consumer reads it via condensed_mass_g.
            chain_mass_defect_g_mol=float(
                p.get("chain_mass_defect_g_mol") or 0.0),
        ))
        if initial_moments and lab in initial_moments:
            moments0[lab] = tuple(initial_moments[lab])
        elif p.get("moments") is not None:
            moments0[lab] = tuple(p["moments"])

    for m in (mass_transfer_spec or []):
        for role, key in (("gas", m["gas"]), ("poly", m["poly"])):
            if key not in idx:
                raise ValueError(
                    f"--mass-transfer {role} label {key!r} not found in "
                    f"chem.yaml species")
    mt_configs = [MassTransferConfig(gas_index=idx[m["gas"]],
                                     poly_index=idx[m["poly"]],
                                     K=float(m["K"]), kLa=float(m["kLa"]))
                  for m in (mass_transfer_spec or [])]

    init_moles = {}
    for label, v in (initial_moles or {}).items():
        if label not in idx:
            raise ValueError(
                f"--initial-moles label {label!r} not found in chem.yaml species")
        init_moles[core[idx[label]]] = float(v)

    rs = HybridPolymerSystem(
        T=(T0, "K"), P=(P, "Pa"),
        initial_mole_fractions=init_moles,  # interpreted as moles
        V_poly=float(V_poly),
        polymer_pools=pools, mass_transfer=mt_configs,
        gas_species_mask=mask, constant_gas_volume=False,
        initial_polymer_moments=moments0,
        # Stage-A solver contract {pool_label: {dp: moles}} recovered from
        # the artifact's explicit_dp.initial_moles (set_initial_conditions
        # step 2 seeds it verbatim, matching the generation run).
        initial_explicit_species=initial_explicit_by_pool,
        termination=[],
        # Item 17 A5-2: the runner is a direct (no-blueprint-phase) construction
        # -- a legitimate last-resort prospective-mask fallback. Flag it so the
        # R1-EDGE provenance guard does not raise on a default-filled edge
        # suffix. (The runner passes edge_species=[] today, so its edge suffix
        # is empty and provenance is vacuously clean even without the flag, but
        # set it explicitly for clarity and to cover non-empty edges.)
        allow_default_prospective_edge=True,
        # r71 FIX 4 escape (documented runner posture): the oracle replays
        # historical decks faithfully, including legacy-schema decks whose
        # pool-mapped rows carry no archetype in the sidecar (pre-archetype
        # emitters). Those rows ran legacy mu1-only flux in the generating
        # solver, and the oracle must reproduce that flux, not refuse it.
        # Production RMG builds (rmgpy/rmg/polymer_input.py) leave the flag
        # default-False and hard-fail on unstamped live proxy rows.
        allow_unstamped_proxy_rows=True)
    with contextlib.redirect_stdout(io.StringIO()):  # mute the mapping banner
        rs.initialize_model(core, all_reactions, [], [])
    return rs, core, all_reactions


def run_segments(rs, core, artifact, all_reactions, segments,
                 n_points_per_segment=50):
    """Piecewise-isothermal integration. ``segments`` = [(t_end, T_K), ...]
    with strictly increasing t_end. Returns rows:
    [t, T, <pool>_mu0.., <pool>_mu1.., <pool>_mu2.. per configured pool,
     n(monomer_routing) per routed pool]."""
    conv = artifact["conventions"]
    idx = {s.label: i for i, s in enumerate(core)}
    pool_labels = list(conv["configured_pools"])
    mu_cols = [(lab, tuple(idx[f"{lab}_mu{k}"] for k in range(3)))
               for lab in pool_labels]
    routed = [(p["label"], idx[p["monomer_routing"]])
              for p in artifact["pools"]
              if p["label"] in pool_labels and p.get("monomer_routing")]

    t_ends = [t for t, _ in segments]
    if any(t <= 0 for t in t_ends) or t_ends != sorted(t_ends) or len(t_ends) != len(set(t_ends)):
        raise ValueError(
            f"t-profile t_end values must be strictly increasing and positive; "
            f"got {t_ends}")

    rows = []
    t_start = 0.0
    y_carry = None
    for seg_i, (t_end, T_k) in enumerate(segments):
        if seg_i > 0:
            # the proven restart pattern (see module docstring)
            rs.T = Quantity((T_k, "K"))
            rs.generate_rate_coefficients(all_reactions, [])
            rs.t0 = t_start
            rs.y0 = y_carry.copy()
            rs.set_initial_derivative()
            rs.initialize_solver()
        for t in np.linspace(t_start, t_end, n_points_per_segment + 1)[1:]:
            rs.advance(t)
            y = np.asarray(rs.y)
            row = [float(t), float(T_k)]
            for _lab, (i0, i1, i2) in mu_cols:
                row.extend([float(y[i0]), float(y[i1]), float(y[i2])])
            for _lab, ri in routed:
                row.append(float(y[ri]))
            rows.append(row)
        y_carry = np.asarray(rs.y).copy()
        t_start = t_end
    return rows


def _csv_header(artifact):
    conv = artifact["conventions"]
    header = ["t_s", "T_K"]
    for lab in conv["configured_pools"]:
        header.extend([f"{lab}_mu0_mol", f"{lab}_mu1_mol", f"{lab}_mu2_mol"])
    for p in artifact["pools"]:
        if p["label"] in conv["configured_pools"] and p.get("monomer_routing"):
            header.append(f"n_{p['monomer_routing']}_mol")
    return header


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Polymer moments artifact reference runner (oracle). "
                    "See docs/polymer_moments_format.md.")
    parser.add_argument("--artifact", required=True,
                        help="polymer_pools.json (schema 2.0/2.1/2.2)")
    parser.add_argument("--chem", required=True, help="chem.yaml from the same RMG run")
    parser.add_argument("--t-profile", required=True,
                        help="JSON: [{\"t_end\": s, \"T\": K}, ...] piecewise-isothermal")
    parser.add_argument("--n-points", type=int, default=50,
                        help="output points per segment (default 50)")
    parser.add_argument("--pressure", type=float, default=1.0e5, help="Pa")
    parser.add_argument("--v-poly", type=float, required=True, help="m^3")
    parser.add_argument("--initial-moles", required=True,
                        help="JSON: {chem.yaml label: mol}")
    parser.add_argument("--initial-moments", default=None,
                        help="JSON: {pool label: [mu0, mu1, mu2] mol}; "
                             "default = artifact pools[].moments")
    parser.add_argument("--mass-transfer", default=None,
                        help="JSON: [{gas, poly, K, kLa}] (labels; operating "
                             "condition, not artifact content)")
    parser.add_argument("--output", required=True, help="CSV path")
    args = parser.parse_args(argv)

    with open(args.artifact) as fh:
        artifact = json.load(fh)
    if not str(artifact.get("schema_version", "")).startswith("2."):
        sys.exit(f"artifact schema_version {artifact.get('schema_version')!r} "
                 "is not 2.x — regenerate with a current RMG-Py polymer branch")
    with open(args.t_profile) as fh:
        profile = [(float(seg["t_end"]), float(seg["T"])) for seg in json.load(fh)]
    with open(args.initial_moles) as fh:
        initial_moles = json.load(fh)
    initial_moments = None
    if args.initial_moments:
        with open(args.initial_moments) as fh:
            initial_moments = json.load(fh)
    mass_transfer_spec = []
    if args.mass_transfer:
        with open(args.mass_transfer) as fh:
            mass_transfer_spec = json.load(fh)

    species, reactions = load_chem_yaml(args.chem)
    rs, core, all_reactions = build_system_from_artifact(
        artifact, species, reactions,
        T0=profile[0][1], P=args.pressure, V_poly=args.v_poly,
        initial_moles=initial_moles, mass_transfer_spec=mass_transfer_spec,
        initial_moments=initial_moments)
    rows = run_segments(rs, core, artifact, all_reactions, profile,
                        n_points_per_segment=args.n_points)

    with open(args.output, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(_csv_header(artifact))
        writer.writerows(rows)
    print(f"wrote {len(rows)} rows to {args.output}")


if __name__ == "__main__":
    main()
