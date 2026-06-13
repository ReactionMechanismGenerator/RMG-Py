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
        pools.append(PolymerPoolConfig(
            label=lab, xs=int(p.get("cutoff") or 0),
            explicit_dp_to_species_index={},
            mu_indices=mu_idx,
            monomer_poly_index=idx[routing] if routing else None,
            # The tripwire's ONE chain-scale window (and the spawn-gate
            # snapshot) consume this; without it the consumer-world window
            # collapses to the bare slack and the tag-branch class drifts
            # from generation world.
            monomer_mw_g_mol=float(p.get("monomer_mw_g_mol") or 0.0),
            k_scission=float(p["channels"]["scission"]["A"]),
            k_unzip=float(p["channels"]["unzip"]["A"]),
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
    parser.add_argument("--artifact", required=True, help="polymer_pools.json (schema 2.0)")
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
