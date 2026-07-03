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
}


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
        arch = ARCHETYPE_INTS[e["archetype"]]
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

# The QSSA channel vocabulary entered the sidecar schema at 2.1 (channel-
# vocabulary growth = minor bump); the emitter stamps >= 2.1 whenever it
# writes the block, so a 2.0 artifact carrying one is malformed. The
# weak-link allyl/U-state sub-vocabulary entered at 2.2, same rule.
_QSSA_MIN_SCHEMA_MINOR = 1
_WEAKLINK_MIN_SCHEMA_MINOR = 2
# Maximum schema minor this loader implements. Weak-link milestone iv
# POLICY CHANGE (was minor-permissive): a 2.3+ artifact may carry
# vocabulary OUTSIDE the channel blocks (new conventions, new pool fields)
# that the unknown-key guards here never inspect, so an older loader must
# fail loud on it instead of loading additively.
_MAX_KNOWN_SCHEMA_MINOR = 2


def _check_schema_version_known(artifact):
    """Reject any artifact whose schema_version is not one this loader
    implements (2.0 <= 2.x <= 2.2)."""
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
})
_PRE_MONOMER_GAS_RECIPE_REVISIONS = frozenset({
    "2026-06-10",
    "2026-07-02",
    "2026-07-03-weaklink-u",
})


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
        pools.append(PolymerPoolConfig(
            label=lab, xs=int(p.get("cutoff") or 0),
            explicit_dp_to_species_index={},
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
        initial_polymer_moments=moments0, termination=[],
        # Item 17 A5-2: the runner is a direct (no-blueprint-phase) construction
        # -- a legitimate last-resort prospective-mask fallback. Flag it so the
        # R1-EDGE provenance guard does not raise on a default-filled edge
        # suffix. (The runner passes edge_species=[] today, so its edge suffix
        # is empty and provenance is vacuously clean even without the flag, but
        # set it explicitly for clarity and to cover non-empty edges.)
        allow_default_prospective_edge=True)
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
