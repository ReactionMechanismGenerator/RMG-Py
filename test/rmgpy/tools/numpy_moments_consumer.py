"""Numpy-only reference consumer for the polymer moments artifact.

TEST / TA-REFERENCE TWIN -- NOT A PRODUCTION LOADER (r42 P2). This module
lives under test/ and exists to verify the format doc's trajectory
semantics independently of the compiled oracle. It INTENTIONALLY DIVERGES
from the production reference runner on one point: it SIMULATES
moment_credit_conduit/1 rows (the DESIGN §2.2 step-6 bundle) that
rmgpy/tools/polymer_moments_runner.py hard-REFUSES to replay pending the
M18.4 solver dispatch arm (the runner's compiled oracle forbids any
consumer-side moment write). That divergence is deliberate -- it is the
only normative trajectory verification of the conduit bundle law -- and
it MUST NOT leak into production dispatch: never route a production
replay through this module, and never relax the runner's fail-closed
refusal to match this twin.

Implements docs/polymer_moments_format.md §4-§7 with numpy + stdlib ONLY.
THIS MODULE MUST NOT IMPORT rmgpy (that is the artifact's entire point);
test_consumer_module_is_rmgpy_free enforces it.

Rate constants come from each entry's own `kinetics` (every entry in the
test decks is emitted cantera-null), so no Cantera is needed here either;
Keq for reversible entries is computed from caller-supplied NASA7 data via
the documented recipe.
"""

import numpy as np

# MUST equal the oracle's rmgpy.constants.R (CODATA-2006), NOT the 2018-SI
# exact value 8.31446261815324: a mismatched R shifts every gas-phase / Keq /
# mass-transfer term by a constant 1.13e-6 (format doc §4 definitions).
R = 8.314472          # J/(mol K)
P0 = 1.0e5            # Pa, Keq reference pressure
SMALL_EPS = 1e-30
LN_EXP_OVERFLOW_GUARD = 700.0
# Near-exhaustion bundle limiter band (round-27 P1-A + round-29 N2 C1
# soft-min; keep in sync with the generating solver's
# BUNDLE_LIMITER_E_LO/_E_HI/_SOFTMIN_P, rmgpy/solver/polymer.pyx):
# dimensionless accepted-state floor distances E = softmin_p(mu_k/f_k) with
# f_k = max(SMALL_EPS, 100*atol). E >= E_HI: S_free == S_base exactly (bulk
# untouched); E <= E_LO: soft bundle cap fully active (softmin_p over the
# per-moment caps and S_base); between: C1 smoothstep blend. S_free then
# passes the INDEPENDENT round-30 cone-margin drain gate (its own M band
# below) to give S_eff.
BUNDLE_LIMITER_E_LO = 1.0e2
BUNDLE_LIMITER_E_HI = 1.0e4
BUNDLE_LIMITER_SOFTMIN_P = 8.0
# Cone-margin drain gate band (round-30 N1; keep in sync with the
# generating solver's CONE_MARGIN_M_LO/_M_HI): dimensionless margin
# distance M = Q10/f with Q10 = mu1 - mu0 and f the r81 floor. For
# cone-shrinking debits (b1 > b0 = 1): Q10 <= 0 returns 0 REGARDLESS of
# E; M >= M_HI returns S_free exactly; M <= M_LO applies
# softmin_p(S_free, S_cone), S_cone = Q10/(V_poly*(b1 - 1)); between: C1
# smoothstep blend.
CONE_MARGIN_M_LO = 1.0e2
CONE_MARGIN_M_HI = 1.0e4


def softmin_p(terms, p=BUNDLE_LIMITER_SOFTMIN_P):
    """Reciprocal p-norm soft-min over positive terms (round-29 N2):
    softmin_p(x_i) = (sum x_i^(-p))^(-1/p), evaluated in the overflow-safe
    factored form m*(sum (m/x_i)^p)^(-1/p) with m the hard min
    (algebraically identical). Any non-positive term zeroes the result
    (softmin_p <= min <= every term)."""
    m = min(terms)
    if m <= 0.0:
        return 0.0
    return m * sum((m / x) ** p for x in terms) ** (-1.0 / p)

# The only channel kinds this consumer dispatches. Any other key inside a
# pool's ``channels`` (e.g. radical_qssa_unzip) must fail at construction,
# not be dropped permissively -- silently integrating without an enabled
# channel produces a flat/false trajectory. Mirrors the TA loader's guard
# (~/Code/TA/ta/mechanism.py SUPPORTED_CHANNELS).
SUPPORTED_CHANNELS = frozenset({"scission", "unzip"})

# The CLOSED refused_reason vocabulary (schema 2.4, format doc §12). Local
# copy -- this module never imports rmgpy -- of the reference loader's
# REFUSED_REASONS pin: a refused row's WHOLE flux (moment and species
# alike) was zeroed by the generating solver, so this consumer must skip
# the row entirely, never convert its kinetics into flux (round-67 ruling
# (c): refused explicit homolysis rows superseded by a live k_homolysis
# kernel stay zero). Unknown reasons are rejected at construction, never
# adapted.
REFUSED_REASONS = frozenset({"conduit-deferred", "qssa-invalid",
                             "qssa-unassessable"})

# ---------------------------------------------------------------------------
# M18.3 (T7 commit 2): envelope gate + CLOSED archetype vocabulary.
# Before this fix the module had NO schema_version reference and an if/elif
# archetype dispatch with no else -- a one-row artifact with an unknown
# archetype under an arbitrary stamp integrated its step-5 species dispatch
# while crediting zero pool moments: silent condensed-mass fabrication,
# PROVEN by TestSilentMassFabricationProof (commit 1 of the ratified
# test-first rider) at ~2 mol of gas fabricated over the standard 0.2 s
# window with every pool moment bit-identical.
# ---------------------------------------------------------------------------

# Accepted envelope, IDENTICAL to the reference runner's
# (rmgpy/tools/polymer_moments_runner.py _check_schema_version_known):
# 2.0 .. 2.8 and 3.0 .. 3.1. Local pins -- this module never imports rmgpy.
KNOWN_SCHEMA_2_MAX_MINOR = 8
KNOWN_SCHEMA_3_MAX_MINOR = 1
# Comparable ordinal ladder (the runner's _schema_minor): 2.x -> x,
# 3.y -> 10 + y (the 2.x line is CLOSED at 2.9). Malformed -> -1.
CONDUIT_MIN_SCHEMA_ORDINAL = 11

# The CLOSED archetype vocabulary this consumer dispatches (mirror of the
# runner's ARCHETYPE_INTS). An unknown archetype is flux this consumer
# cannot reproduce -- reject at construction, never integrate around it.
KNOWN_ARCHETYPES = frozenset({
    "same_pool/1", "migration/1", "scission_fragment/1", "legacy_mu1/1",
    "discrete_chip/1", "volatile_ejection/1", "moment_credit_conduit/1",
})


def _schema_ordinal(ver):
    parts = str(ver).split(".")
    if len(parts) == 2 and parts[1].isdigit():
        if parts[0] == "2":
            return int(parts[1])
        if parts[0] == "3":
            return 10 + int(parts[1])
    return -1


def _check_schema_version(artifact):
    """Envelope gate (M18.3): accept exactly the reference runner's set --
    2.0..2.8 and 3.0..3.1 -- so the two reference consumers agree on the
    refusal surface. A newer minor may carry vocabulary this consumer
    would silently ignore; refuse loud, never load permissively."""
    ver = str(artifact.get("schema_version", ""))
    parts = ver.split(".")
    ok = (len(parts) == 2 and parts[1].isdigit()
          and ((parts[0] == "2"
                and int(parts[1]) <= KNOWN_SCHEMA_2_MAX_MINOR)
               or (parts[0] == "3"
                   and int(parts[1]) <= KNOWN_SCHEMA_3_MAX_MINOR)))
    if not ok:
        raise ValueError(
            f"artifact schema_version {ver!r} is not implemented by this "
            f"consumer (known: 2.0 .. 2.{KNOWN_SCHEMA_2_MAX_MINOR}, "
            f"3.0 .. 3.{KNOWN_SCHEMA_3_MAX_MINOR}). A newer minor may "
            f"carry vocabulary this consumer would silently ignore -- "
            f"upgrade the consumer or regenerate the sidecar.")


def _validated_conduit(e, pools, configured, condensed):
    """rmgpy-free mirror of the reference loader's load-bearing
    moment_credit_conduit/1 reject rules (DESIGN §2.1; the full §2.4
    chem.yaml cross-pins live in the reference runner, which holds the
    composition data this consumer deliberately does not). Returns the
    validated (u, dst) pair."""
    eid = e.get("id")

    def _bad(msg):
        raise ValueError(
            f"reactions[] entry {eid!r} (moment_credit_conduit/1) {msg} "
            f"Fix the artifact.")

    if e.get("cantera") is None:
        _bad("carries cantera: null -- the conduit contract requires the "
             "row to exist in chem.yaml at its index.")
    equation = str((e.get("cantera") or {}).get("equation", ""))
    if "<=>" in equation or "=>" not in equation:
        _bad(f"must export the irreversible arrow '=>' (got {equation!r}).")
    kin = e.get("kinetics")
    if not isinstance(kin, dict) or kin.get("reversible") is not False:
        _bad("must carry kinetics.reversible: false.")
    if e.get("src_pool") is not None:
        _bad(f"must carry src_pool: null (got {e.get('src_pool')!r}).")
    dst = e.get("dst_pool")
    if not dst or dst not in pools or dst not in configured:
        _bad(f"must name a serialized AND configured destination pool "
             f"(dst_pool={dst!r}).")
    if e.get("proxy_reactants") != []:
        _bad(f"must carry proxy_reactants: [] (got "
             f"{e.get('proxy_reactants')!r}).")
    if not e.get("proxy_products"):
        _bad("must carry non-empty proxy_products.")
    if "refused" in e or "refused_reason" in e:
        _bad("also carries refused vocabulary -- mutually exclusive.")
    params = e.get("params")
    required = {"admission_direction", "chain_units", "gas_products",
                "gas_units", "candidate_key"}
    allowed = required | {"candidate_key_note"}
    if not isinstance(params, dict) or not required <= set(params) \
            or not set(params) <= allowed:
        _bad(f"must carry the closed §2.1 params vocabulary; got "
             f"{sorted(params) if isinstance(params, dict) else params!r}.")
    if params.get("admission_direction") != "chain_to_pool":
        _bad(f"admission_direction is CLOSED to 'chain_to_pool' (got "
             f"{params.get('admission_direction')!r}).")
    u = params.get("chain_units")
    if isinstance(u, bool) or not isinstance(u, (int, float)) \
            or not np.isfinite(float(u)) or float(u) < 1.0:
        _bad(f"params.chain_units={u!r} must be a finite float >= 1.0 "
             f"(landing cone, §2.3).")
    gps = params.get("gas_products")
    if not isinstance(gps, list) or len(gps) != 1 \
            or not isinstance(gps[0], dict) or gps[0].get("stoich") != 1:
        _bad(f"must carry EXACTLY ONE gas product with stoich 1 "
             f"[r39-P5]; got {gps!r}.")
    gas_label = str(gps[0].get("species"))
    if gas_label not in (e.get("products") or []):
        _bad(f"gas product {gas_label!r} is not among the row's products.")
    if gas_label in condensed:
        _bad(f"gas product {gas_label!r} is condensed-classified.")
    if not any(r in condensed for r in (e.get("reactants") or [])):
        _bad("no reactant is melt-classified: the conduit event is "
             "condensed by contract (§2.2).")
    # candidate_key semantic pin [r42 P1-5] (mirror of the reference
    # runner's check): non-empty AND recomputes from the row's OWN
    # serialized identity (sides sorted + '+'-joined, ordered
    # lexicographically around '<>'); the §4.4 flux-census accounting is
    # partitioned on it.
    key = params.get("candidate_key")
    if not isinstance(key, str) or not key.strip():
        _bad(f"params.candidate_key={key!r} must be a non-empty string "
             f"(r42 P1-5).")
    side_a = "+".join(sorted(str(t) for t in (e.get("reactants") or [])))
    side_b = "+".join(sorted(str(t) for t in (e.get("products") or [])))
    lo, hi = sorted((side_a, side_b))
    key_recomputed = f"{lo}<>{hi}"
    if key != key_recomputed:
        _bad(f"params.candidate_key={key!r} does not recompute from the "
             f"row's own serialized reactants/products (expected "
             f"{key_recomputed!r}; r42 P1-5).")
    return float(u), dst


def _validated_refused(e):
    """Strict shape guard for one reactions[] entry's refused vocabulary
    (format doc §12; the reference loader's _validate_refused_entry,
    mirrored rmgpy-free): the emitter writes ``refused: true`` plus a
    known ``refused_reason`` on POOL-MAPPED rows, or NOTHING (absent, not
    false). Returns True iff the row is a valid refused row."""
    if "refused" not in e and "refused_reason" not in e:
        return False
    eid = e.get("id")
    if e.get("refused") is not True:
        raise ValueError(
            f"reactions[] entry {eid!r} carries malformed refused "
            f"vocabulary (refused={e.get('refused')!r}): the emitter "
            f"writes refused: true + refused_reason on refused rows and "
            f"NOTHING otherwise (absent, not false). Fix the artifact.")
    reason = e.get("refused_reason")
    if not isinstance(reason, str) or reason not in REFUSED_REASONS:
        raise ValueError(
            f"reactions[] entry {eid!r} has invalid refused_reason "
            f"{reason!r}: the schema-2.4 vocabulary is CLOSED to "
            f"{sorted(REFUSED_REASONS)} (format doc §12). Fix the "
            f"artifact.")
    if not (e.get("proxy_reactants") or e.get("proxy_products")):
        raise ValueError(
            f"reactions[] entry {eid!r} has refused: true but no "
            f"pool-mapped participant: the refused marker is legal ONLY "
            f"on pool-mapped rows -- honoring it here would silently zero "
            f"ordinary gas chemistry. Fix the artifact.")
    return True


def _validated_eject_units(e):
    """Validate + return a volatile_ejection/1 row's SIGNED eject_units
    (rmgpy-free mirror of the reference loader's _validated_eject_units,
    rmgpy/tools/polymer_moments_runner.py): the emitter writes exactly ONE
    VE params sub-shape, ``params = {"eject_units": <signed float>}``.
    Reject anything else loudly, never KeyError and never a silent 0.0
    default -- defaulting would launder the atom-transfer debit away (the
    moved chain lands un-shrunk while the gas volatile still appears:
    fabricated mass)."""
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
            or not np.isfinite(float(a)):
        raise ValueError(
            f"reactions[] entry {eid!r} (volatile_ejection/1) has "
            f"eject_units={a!r}; it must be a finite SIGNED number (the "
            f"source-monomer-equivalents transferred to the gas "
            f"co-participants). Fix the artifact.")
    return float(a)


def safe_mu3(mu0, mu1, mu2):
    """log_lagrange/1 closure with realizability guard (format doc §6)."""
    if mu0 <= SMALL_EPS or mu1 <= SMALL_EPS or mu2 <= SMALL_EPS:
        return 0.0
    if mu1 < mu0:
        return 0.0
    ln_mu3 = 3.0 * np.log(mu2) - 3.0 * np.log(mu1) + np.log(mu0)
    if ln_mu3 > LN_EXP_OVERFLOW_GUARD:
        return float("inf")
    return float(np.exp(ln_mu3))


def nasa_g_over_rt(coeffs, T):
    """G/(R*T) from one NASA7 row [a0..a4, a5, a6]."""
    a0, a1, a2, a3, a4, a5, a6 = coeffs
    h_rt = (a0 + a1 * T / 2.0 + a2 * T**2 / 3.0 + a3 * T**3 / 4.0
            + a4 * T**4 / 5.0 + a5 / T)
    s_r = (a0 * np.log(T) + a1 * T + a2 * T**2 / 2.0 + a3 * T**3 / 3.0
           + a4 * T**4 / 4.0 + a6)
    return h_rt - s_r


def _base(label):
    """Strip ONLY a trailing RMG index suffix ``(<int>)`` from a label
    ('PS(2)' -> 'PS', 'C[CH]CC(C)C(6)' -> 'C[CH]CC(C)C'). Local duplicate
    of the ONE canonical convention, rmgpy.polymer.strip_rmg_index_suffix
    (this module deliberately never imports rmgpy): truncating at the
    FIRST '(' mangles SMILES-derived labels, whose parentheses are
    branching syntax; a parenthesised group of bare digits is never valid
    SMILES."""
    if label and label.endswith(")"):
        head, sep, tail = label[:-1].rpartition("(")
        if sep and tail.isdigit():
            return head
    return label


class ArtifactConsumer:
    """Integrates the artifact's pool moments + species ODEs.

    Parameters
    ----------
    artifact : dict           parsed polymer_pools.json (schema 2.0)
    species_order : [str]     chem.yaml labels defining the y-vector layout
                              (same order as the oracle's core species)
    P : float                 pressure [Pa]
    V_poly : float            condensed-phase volume [m^3]
    mass_transfer : [dict]    [{"gas": label, "poly": label, "K": f, "kLa": f}]
    nasa : {label: {"Tmid": f, "low": [7], "high": [7]}}, optional
                              NASA7 data for Keq of reversible entries
    atol : float, optional    generating solver's absolute tolerance (mole
                              basis). Sets the near-exhaustion bundle
                              limiter's per-moment floor
                              f_k = max(SMALL_EPS, 100*atol), mirroring the
                              solver's _pool_mu_floors (scalar-atol form).
                              Must match the oracle's atol for parity.
    """

    def __init__(self, artifact, species_order, P, V_poly,
                 mass_transfer=None, nasa=None, atol=1e-16):
        # M18.3 envelope gate FIRST: refuse unknown schema versions loudly
        # (this consumer used to accept ANY stamp silently -- the proven
        # fabrication precondition).
        _check_schema_version(artifact)
        self.P = float(P)
        self.V_poly = float(V_poly)
        self.nasa = nasa or {}
        # Near-exhaustion bundle limiter floor (round-27 P1-A; keep in sync
        # with the generating solver's initialize_model floors
        # max(SMALL_EPS, EXHAUSTION_FLOOR_K*atol[state]) -- this consumer
        # takes the scalar-atol form, uniform across moment slots).
        self.mu_floor = max(SMALL_EPS, 100.0 * float(atol))
        self.idx = {lab: i for i, lab in enumerate(species_order)}
        n = len(species_order)

        conv = artifact["conventions"]
        condensed = set(conv["condensed_species"])
        self.gas_mask = np.array([lab not in condensed for lab in species_order],
                                 dtype=bool)
        self.configured = list(conv["configured_pools"])

        # pools: label -> dict(mu indices, channels, monomer routing index)
        self.pools = {}
        for p in artifact["pools"]:
            lab = p["label"]
            # Schema-2.6 radical-homolysis initiation kernel: NOT
            # implemented by this consumer. Fail at construction, never
            # drop permissively -- integrating a kernel-enabled melt
            # without its initiation flux produces a flat/false
            # trajectory, and the refused explicit rows it supersedes
            # carry zero flux by contract (round-67 ruling (c)), so
            # nothing here may stand in for the kernel.
            if "homolysis_initiation" in p:
                raise ValueError(
                    f"pool {lab!r}: unsupported pool-level "
                    f"homolysis_initiation block (schema 2.6 "
                    f"radical-homolysis initiation kernel); this consumer "
                    f"implements channels {sorted(SUPPORTED_CHANNELS)} "
                    f"only. Integrating without the kernel would silently "
                    f"produce a wrong trajectory, and refused rows must "
                    f"never be converted into its flux.")
            # Schema-2.7 side-group homolysis kernel + its X-loss
            # feature-pool exact mass contract: NOT implemented by this
            # consumer. Fail at construction, never drop permissively --
            # integrating a kernel-enabled melt without the channel (or a
            # defect-carrying pool without the normative mass formula
            # condensed_mass_g = mu1*monomer_mw_g_mol -
            # mu0*chain_mass_defect_g_mol) silently mints condensed mass
            # while gas X appears (the round-70 P1 trap), and the refused
            # explicit rows the kernel supersedes carry zero flux by
            # contract, so nothing here may stand in for it.
            # Schema-2.8 end-radical depropagation kernel: NOT implemented
            # by this consumer. Fail at construction, never drop
            # permissively -- integrating a kernel-enabled melt without
            # the channel leaves the radical-end pools outlet-free (the
            # run-6 no-outlet wall) and drops the gas monomer source the
            # generating solver emitted (un-conserved mass), silently
            # producing a wrong trajectory.
            if "end_radical_depropagation" in p:
                raise ValueError(
                    f"pool {lab!r}: unsupported pool-level "
                    f"end_radical_depropagation block (schema 2.8 "
                    f"end-radical depropagation kernel); this consumer "
                    f"implements channels {sorted(SUPPORTED_CHANNELS)} "
                    f"only. Integrating without the kernel would silently "
                    f"produce a wrong trajectory (outlet-free radical-end "
                    f"pools, missing gas monomer source).")
            if "side_group_homolysis" in p:
                raise ValueError(
                    f"pool {lab!r}: unsupported schema-2.7 side-group "
                    f"homolysis vocabulary (pool-level "
                    f"side_group_homolysis block); this consumer "
                    f"implements channels "
                    f"{sorted(SUPPORTED_CHANNELS)} only. Integrating "
                    f"without the kernel would silently produce a "
                    f"wrong trajectory, and refused rows must never be "
                    f"converted into its flux.")
            # Pool-level chain_mass_defect_g_mol (additive optional field,
            # P1-2 atom-transfer mass defect): per-chain mass [g/mol]
            # already shed to the gas (e.g. one abstracted H per chain on
            # an H-loss _mod daughter). This consumer implements the SAME
            # defect-aware booking as the generating solver: the VE moment
            # shift uses a_moment = a_mass - (defect_dst - defect_src)/
            # monomer_mw_src (exactly 0 for H-loss conduits -- chains
            # transfer with UNCHANGED mu1), and the condensed-mass closure
            # is mu1*monomer_mw_g_mol - mu0*chain_mass_defect_g_mol.
            defect = p.get("chain_mass_defect_g_mol", 0.0)
            if isinstance(defect, bool) or not isinstance(defect, (int, float)) \
                    or not np.isfinite(float(defect)) or float(defect) < 0.0:
                raise ValueError(
                    f"pool {lab!r}: chain_mass_defect_g_mol={defect!r} "
                    f"must be a finite non-negative per-chain g/mol mass "
                    f"defect. Fix the artifact.")
            if lab not in self.configured:
                continue
            unknown = sorted(set(p["channels"]) - SUPPORTED_CHANNELS)
            if unknown:
                raise ValueError(
                    f"pool {lab!r}: unsupported channel key(s) {unknown}; "
                    f"this consumer implements {sorted(SUPPORTED_CHANNELS)} "
                    f"only. Integrating without an enabled channel would "
                    f"silently produce a wrong trajectory.")
            mu_idx = tuple(self.idx[f"{lab}_mu{k}"] for k in range(3))
            routing = p.get("monomer_routing")
            self.pools[lab] = {
                "mu": mu_idx,
                "k_s": p["channels"]["scission"]["A"],
                "k_u": p["channels"]["unzip"]["A"],
                "routing": self.idx[routing] if routing else None,
                "mw": float(p.get("monomer_mw_g_mol", 0.0) or 0.0),
                "defect": float(defect),
            }

        # reactions: precompute index forms
        self.entries = []
        pool_labels = {p.get("label") for p in artifact.get("pools", [])
                       if isinstance(p, dict)}
        for e in artifact["reactions"]:
            # M18.3 CLOSED archetype dispatch (mirror of the reference
            # runner's raise): an unknown archetype is flux this consumer
            # cannot reproduce -- integrating its step-5 species dispatch
            # anyway fabricates gas mass with no condensed debit (the
            # proven pre-fix behavior). Hard else-refusal, never skip.
            if e["archetype"] not in KNOWN_ARCHETYPES:
                raise ValueError(
                    f"reactions[] entry {e.get('id')!r} carries unknown "
                    f"archetype {e['archetype']!r}; this consumer's "
                    f"term-type vocabulary is CLOSED to "
                    f"{sorted(KNOWN_ARCHETYPES)} (docs/"
                    f"polymer_moments_format.md §3). An unknown archetype "
                    f"is flux this consumer cannot reproduce -- upgrade "
                    f"the consumer or regenerate the sidecar.")
            conduit_u = None
            conduit_dst = None
            if e["archetype"] == "moment_credit_conduit/1":
                if (_schema_ordinal(artifact.get("schema_version"))
                        < CONDUIT_MIN_SCHEMA_ORDINAL):
                    raise ValueError(
                        f"artifact schema_version "
                        f"{artifact.get('schema_version')!r} cannot carry "
                        f"moment_credit_conduit/1 rows: the conduit "
                        f"vocabulary was introduced in schema 3.1 "
                        f"(presence-elected). This artifact is malformed.")
                conduit_u, conduit_dst = _validated_conduit(
                    e, pool_labels, set(self.configured), condensed)
            kin = e["kinetics"]
            assert kin is not None, f"entry {e['id']} has no kinetics"
            ridx = [self.idx[s] for s in e["reactants"]]
            pidx = [self.idx[s] for s in e["products"]]
            pool_mapped = ({self.idx[s] for s in e["proxy_reactants"]}
                           | {self.idx[s] for s in e["proxy_products"]})
            self.entries.append({
                # Refused-row marker (schema 2.4, format doc §12):
                # validated strictly at construction; a refused entry's
                # WHOLE flux is skipped in rhs (moment and species alike),
                # exactly the generating solver's reaction_refused
                # suppression.
                "refused": _validated_refused(e),
                "A": kin["A"], "n": kin["n"], "Ea": kin["Ea"],
                "reversible": kin["reversible"],
                "ridx": ridx, "pidx": pidx, "pool_mapped": pool_mapped,
                "r_labels": list(e["reactants"]), "p_labels": list(e["products"]),
                "proxy_r_pools": [_base(s) for s in e["proxy_reactants"]],
                "proxy_p_pools": [_base(s) for s in e["proxy_products"]],
                "scaling": e["scaling"],
                "src": e["src_pool"], "dst": e["dst_pool"],
                "arch": e["archetype"],
                "a": int(e.get("params", {}).get("a", 0)),
                # volatile_ejection/1 rows carry the SIGNED atom-transfer
                # stamp; validated strictly (never defaulted -- see
                # _validated_eject_units).
                "eject": (_validated_eject_units(e)
                          if e["archetype"] == "volatile_ejection/1"
                          else 0.0),
                # moment_credit_conduit/1 (M18.3, DESIGN §2.2): the
                # validated credit u and destination pool.
                "conduit_u": conduit_u,
                "conduit_dst": conduit_dst,
            })
            # P1-2 atom-transfer mass defect (keep in sync with the
            # generating solver's reaction_eject_units_moment): the artifact
            # eject_units stamp is the MASS defect a_mass; the MOMENT shift
            # for a cross-pool VE row is
            #   a_moment = a_mass - (defect_dst - defect_src)/monomer_mw_src
            # (exactly 0 for H-loss conduits: chains transfer with UNCHANGED
            # mu1 while the gas gains the H mass through the ordinary step-5
            # species path; byte-identical a_moment == a_mass when the pools
            # carry equal defects, e.g. true monomer/chip ejection rows).
            ent = self.entries[-1]
            a_mom = float(ent["eject"])
            if (ent["arch"] == "volatile_ejection/1"
                    and ent["src"] is not None and ent["dst"] is not None
                    and ent["src"] != ent["dst"]
                    and ent["src"] in self.pools
                    and ent["dst"] in self.pools):
                d_src = self.pools[ent["src"]]["defect"]
                d_dst = self.pools[ent["dst"]]["defect"]
                if d_dst != d_src:
                    mw_src = self.pools[ent["src"]]["mw"]
                    if mw_src <= 0.0:
                        raise ValueError(
                            f"entry {e['id']}: pools {ent['src']!r} -> "
                            f"{ent['dst']!r} carry differing "
                            f"chain_mass_defect_g_mol but the source pool "
                            f"has no positive monomer_mw_g_mol; cannot "
                            f"convert the mass defect to a moment shift.")
                    a_mom = a_mom - (d_dst - d_src) / mw_src
                    if abs(a_mom) <= 1.0e-9 * max(1.0, abs(ent["eject"])):
                        a_mom = 0.0
            ent["eject_moment"] = a_mom

        self.mass_transfer = []
        for mt in (mass_transfer or []):
            self.mass_transfer.append((self.idx[mt["gas"]], self.idx[mt["poly"]],
                                       float(mt["K"]), float(mt["kLa"])))

    # ----- helpers ------------------------------------------------------

    def _keq(self, entry, T):
        def g(label):
            d = self.nasa[_label_lookup(self.nasa, label)]
            row = d["low"] if T <= d["Tmid"] else d["high"]
            return nasa_g_over_rt(row, T)  # G/(R T)
        dg_rt = sum(g(s) for s in entry["p_labels"]) - sum(g(s) for s in entry["r_labels"])
        # dn counts ALL species as written (condensed/proxy included) —
        # format doc §4 step 1.
        dn = len(entry["p_labels"]) - len(entry["r_labels"])
        return (P0 / (R * T)) ** dn * np.exp(-dg_rt)

    def _chain_bundle(self, pool, y, end_group):
        i0, i1, i2 = self.pools[pool]["mu"]
        mu0 = max(0.0, y[i0]) / self.V_poly
        mu1 = max(0.0, y[i1]) / self.V_poly
        mu2 = max(0.0, y[i2]) / self.V_poly
        if end_group:
            if mu0 <= SMALL_EPS:
                return 0.0, 0.0, 0.0, False
            return 1.0, mu1 / mu0, mu2 / mu0, True
        if mu1 <= SMALL_EPS:
            return 0.0, 0.0, 0.0, False
        mu3 = safe_mu3(mu0, mu1, mu2)
        if np.isfinite(mu3):
            return 1.0, mu2 / mu1, mu3 / mu1, True
        return 1.0, mu2 / mu1, 0.0, False

    def _bundle_cap(self, pool, y, end_group):
        """P1-1 directional bundle throttle cap (keep in sync with the
        generating solver's _bundle_availability_cap): the largest
        event-site density [mol/m^3] the debited pool sustains so every
        moment's drain vanishes linearly near exhaustion --
        S_cap = softmin_p over positive bundle terms b_k of
        mu_k/(V_poly*b_k) (round-29 N2 soft-min). 0.0 for an empty bundle
        (gas and moment flux throttle off together). The round-30
        cone-margin drain gate is NOT part of this cap -- it is an
        independent second gate in _bundle_limited_site."""
        b0, b1, b2, _ok = self._chain_bundle(pool, y, end_group)
        if b0 <= 0.0:
            return 0.0
        i0, i1, i2 = self.pools[pool]["mu"]
        terms = [max(0.0, y[i0]) / (self.V_poly * b0)]
        if b1 > 0.0:
            terms.append(max(0.0, y[i1]) / (self.V_poly * b1))
        if b2 > 0.0:
            terms.append(max(0.0, y[i2]) / (self.V_poly * b2))
        return softmin_p(terms)

    def _floor_distance(self, pool, y):
        """E: the debited pool's accepted-state floor distance,
        E = softmin_p(mu0/f0, mu1/f1, mu2/f2) with
        f_k = max(SMALL_EPS, 100*atol) (round-29 N2 soft moment-min; keep
        in sync with the generating solver's _pool_floor_distance)."""
        i0, i1, i2 = self.pools[pool]["mu"]
        return softmin_p([max(0.0, y[i0]), max(0.0, y[i1]),
                          max(0.0, y[i2])]) / self.mu_floor

    def _bundle_limited_site(self, pool, y, end_group, s_base):
        """Near-exhaustion bundle limiter (tail-only smoothstep; C1
        soft-min, cone-margin drain guard), round-27 P1-A throttle +
        round-29 N2 soft-min + round-30 N1 independent cone gate (keep in
        sync with the generating solver's _bundle_limited_site). Stage 1
        (E band):
            S_free = w*S_base + (1-w)*softmin_p(S_base, S_cap),
            w = 3e^2 - 2e^3 on e = clamp((E - E_lo)/(E_hi - E_lo), 0, 1);
            E >= E_hi returns s_base EXACTLY (bulk law untouched).
        Stage 2 (M band, INDEPENDENT of E; only cone-shrinking debits,
        b1 > b0 = 1): Q10 = mu1 - mu0 <= 0 -> 0 regardless of E;
        M = Q10/floor >= M_hi -> S_free exactly; M <= M_lo ->
        softmin_p(S_free, S_cone) with S_cone = Q10/(V_poly*(b1 - 1));
        between: C1 v-smoothstep blend."""
        # stage 1: exhaustion tail limiter
        e_dist = self._floor_distance(pool, y)
        if e_dist >= BUNDLE_LIMITER_E_HI:
            s_free = s_base
        else:
            cap = self._bundle_cap(pool, y, end_group)
            if cap <= 0.0 or s_base <= 0.0:
                cap = 0.0
            else:
                cap = softmin_p([s_base, cap])
            if e_dist <= BUNDLE_LIMITER_E_LO:
                s_free = cap
            else:
                e_n = ((e_dist - BUNDLE_LIMITER_E_LO)
                       / (BUNDLE_LIMITER_E_HI - BUNDLE_LIMITER_E_LO))
                w = e_n * e_n * (3.0 - 2.0 * e_n)
                s_free = w * s_base + (1.0 - w) * cap
        # stage 2: cone-margin drain gate (independent of E)
        i0, i1, i2 = self.pools[pool]["mu"]
        y0c = max(0.0, y[i0])
        y1c = max(0.0, y[i1])
        y2c = max(0.0, y[i2])
        if end_group:
            if y0c / self.V_poly <= SMALL_EPS:
                return s_free       # empty pool: stage 1 owns it
            b1c = y1c / y0c
        else:
            if y1c / self.V_poly <= SMALL_EPS:
                return s_free
            b1c = y2c / y1c
        if b1c <= 1.0:              # b0 = 1: debit does not shrink Q10
            return s_free
        q10 = y1c - y0c
        if q10 <= 0.0:
            return 0.0              # out of cone: drain OFF, whatever E
        m_dist = q10 / self.mu_floor
        if m_dist >= CONE_MARGIN_M_HI:
            return s_free           # margin safely bulk: gate inactive
        s_cone = q10 / (self.V_poly * (b1c - 1.0))
        if s_free <= 0.0:
            return s_free
        cap = softmin_p([s_free, s_cone])
        if m_dist <= CONE_MARGIN_M_LO:
            return cap
        v_n = ((m_dist - CONE_MARGIN_M_LO)
               / (CONE_MARGIN_M_HI - CONE_MARGIN_M_LO))
        v = v_n * v_n * (3.0 - 2.0 * v_n)
        return v * s_free + (1.0 - v) * cap

    def condensed_mass_g(self, y):
        """Defect-aware condensed-mass closure over all configured pools
        (P1-2; keep in sync with the generating solver's
        get_total_polymer_condensed_mass_g):
        sum_p max(0, mu1_p*monomer_mw_g_mol_p -
        mu0_p*chain_mass_defect_g_mol_p) [g]."""
        total = 0.0
        for pool in self.pools.values():
            i0, i1, _ = pool["mu"]
            total += max(0.0, max(0.0, y[i1]) * pool["mw"]
                         - max(0.0, y[i0]) * pool["defect"])
        return total

    # ----- RHS (format doc §4-§7) ----------------------------------------

    def rhs(self, y, T):
        dn = np.zeros_like(y)
        n_gas = float(np.sum(np.clip(y[self.gas_mask], 0.0, None)))
        V_gas = n_gas * R * T / self.P if n_gas > 0 else 1.0
        Vp = self.V_poly
        C = np.where(self.gas_mask,
                     np.clip(y, 0.0, None) / V_gas,
                     np.clip(y, 0.0, None) / Vp)

        for e in self.entries:
            if e["refused"]:
                # Stamp-but-keep (format doc §12): the row stays listed
                # (consumers may need it to zero a mapped Cantera
                # reaction) but its ENTIRE flux is suppressed.
                continue
            kf = e["A"] * T ** e["n"] * np.exp(-e["Ea"] / (R * T))
            kb = kf / self._keq(e, T) if e["reversible"] else 0.0

            # step 2: phase + gate
            is_poly_event = any(not self.gas_mask[i] for i in e["ridx"])
            V_rxn = Vp if is_poly_event else V_gas
            has_poly_prod = any(not self.gas_mask[i] for i in e["pidx"])
            if is_poly_event and not has_poly_prod:
                continue
            if (not is_poly_event) and has_poly_prod:
                continue

            # step 3: concentration products
            rf = kf
            for i in e["ridx"]:
                rf *= 1.0 if i in e["pool_mapped"] else C[i]
            rr = kb
            for i in e["pidx"]:
                rr *= 1.0 if i in e["pool_mapped"] else C[i]

            # step 4: site scaling
            if e["src"] is not None:
                i0, i1, _ = self.pools[e["src"]]["mu"]
                is_ve = e["arch"] == "volatile_ejection/1"
                if e["arch"] == "discrete_chip/1" and e["scaling"] == "mu0" and e["a"] > 0:
                    site = min(max(0.0, y[i0]), max(0.0, y[i1]) / e["a"]) / Vp
                elif (is_ve and e["scaling"] == "mu0"
                        and e["eject_moment"] > 0.0):
                    # a>0 END-GROUP VE exhaustion throttle -- parity with
                    # the generating solver's section-2 site scaling
                    # (polymer.pyx; keep in sync): site = min(mu0, mu1/a);
                    # a<0 GROWS the chain (exempt, no spurious negative
                    # site from mu1/a). Ruling round 20 item C: applies to
                    # same-pool AND cross-pool rows (the forward leg debits
                    # src either way). P1-2: the throttle rides the MOMENT
                    # shift a_moment, not the mass stamp.
                    site = min(max(0.0, y[i0]),
                               max(0.0, y[i1]) / e["eject_moment"]) / Vp
                else:
                    mi = i0 if e["scaling"] == "mu0" else i1
                    site = max(0.0, y[mi]) / Vp
                # P1-1 directional bundle throttle, round-27 P1-A tail-only
                # form -- parity with the generating solver's section-2
                # scaling (keep in sync): each debited direction of a
                # cross-pool VE/MIGRATION leg runs at the near-exhaustion
                # bundle limiter (tail-only smoothstep; C1 soft-min,
                # cone-margin drain guard)
                # S_eff = w*S_base + (1-w)*softmin_p(S_base, S_cap) of the
                # pool it debits, BEFORE rf - rr, so species flux and
                # moment flux never diverge; S_eff == S_base exactly in
                # bulk.
                cross_bt = (e["arch"] in ("volatile_ejection/1",
                                          "migration/1")
                            and e["dst"] is not None
                            and e["dst"] != e["src"])
                if cross_bt:
                    site = self._bundle_limited_site(
                        e["src"], y, e["scaling"] == "mu0", site)
                rf *= site
                if cross_bt:
                    # Direction-specific source availability for cross-pool
                    # VE/MIGRATION (parity with the solver's adjudicated
                    # Part C scaling): the reverse leg debits the DST pool,
                    # so its site factor comes from the dst pool's own
                    # moments (same moment order as the forward site). Item
                    # C mirror direction: for a<0 the reverse leg sheds
                    # |a|/event from dst, so its end-group availability
                    # carries the same min(mu0, mu1/|a|) throttle. P1-1
                    # (round-27 P1-A tail-only form): plus the same
                    # near-exhaustion bundle limiter (tail-only smoothstep;
                    # C1 soft-min, cone-margin drain guard) from dst's own
                    # moments, floor distance and cone margin.
                    d0, d1, _ = self.pools[e["dst"]]["mu"]
                    if (is_ve and e["scaling"] == "mu0"
                            and e["eject_moment"] < 0.0):
                        site_rev = min(
                            max(0.0, y[d0]),
                            max(0.0, y[d1]) / (-e["eject_moment"])) / Vp
                    else:
                        di = d0 if e["scaling"] == "mu0" else d1
                        site_rev = max(0.0, y[di]) / Vp
                    site_rev = self._bundle_limited_site(
                        e["dst"], y, e["scaling"] == "mu0", site_rev)
                    rr *= site_rev
                else:
                    rr *= site

            r_mol = (rf - rr) * V_rxn

            # step 5: stoichiometric flux for non-pool-mapped species
            for i in e["ridx"]:
                if i not in e["pool_mapped"]:
                    dn[i] -= r_mol
            for i in e["pidx"]:
                if i not in e["pool_mapped"]:
                    dn[i] += r_mol

            # step 6: archetype bundles
            arch = e["arch"]
            if arch == "migration/1":
                src, dst = e["src"], e["dst"]
                if src and dst and src != dst:
                    for ev, frm, to in ((rf, src, dst), (rr, dst, src)):
                        if ev <= 0.0:
                            continue
                        ev_mol = ev * V_rxn
                        b0, b1, b2, ok = self._chain_bundle(frm, y, e["scaling"] == "mu0")
                        if b0 == 0.0:
                            continue
                        f = self.pools[frm]["mu"]
                        t = self.pools[to]["mu"]
                        dn[f[0]] -= ev_mol * b0
                        dn[f[1]] -= ev_mol * b1
                        dn[t[0]] += ev_mol * b0
                        dn[t[1]] += ev_mol * b1
                        if ok:
                            dn[f[2]] -= ev_mol * b2
                            dn[t[2]] += ev_mol * b2
            elif arch == "scission_fragment/1":
                src, dst = e["src"], e["dst"]
                if src and dst and src != dst:
                    s = self.pools[src]["mu"]
                    d = self.pools[dst]["mu"]
                    mu0p = max(0.0, y[s[0]]) / Vp
                    mu1p = max(0.0, y[s[1]]) / Vp
                    mu2p = max(0.0, y[s[2]]) / Vp
                    ok = mu1p > SMALL_EPS
                    if ok and r_mol < 0.0:
                        if (max(0.0, y[d[0]]) / Vp <= SMALL_EPS
                                or max(0.0, y[d[1]]) / Vp <= SMALL_EPS):
                            ok = False
                    if ok:
                        e_n = mu2p / mu1p
                        dn[s[1]] -= r_mol * e_n / 2.0
                        dn[d[0]] += r_mol
                        dn[d[1]] += r_mol * e_n / 2.0
                        mu3p = safe_mu3(mu0p, mu1p, mu2p)
                        if np.isfinite(mu3p):
                            e_n2 = mu3p / mu1p
                            dn[s[2]] -= r_mol * (2.0 / 3.0) * e_n2
                            dn[d[2]] += r_mol * e_n2 / 3.0
            elif arch == "discrete_chip/1":
                src = e["src"]
                if src:
                    a = float(e["a"])
                    b0, b1, _b2, _ok = self._chain_bundle(src, y, e["scaling"] == "mu0")
                    if b0 != 0.0:
                        s = self.pools[src]["mu"]
                        e_n = b1
                        if rf > 0.0:
                            rf_mol = rf * V_rxn
                            dn[s[1]] -= rf_mol * a
                            dmu2 = 2.0 * a * e_n - a * a
                            if dmu2 > 0.0:
                                dn[s[2]] -= rf_mol * dmu2
                        if rr > 0.0:
                            rr_mol = rr * V_rxn
                            dn[s[1]] += rr_mol * a
                            dn[s[2]] += rr_mol * (2.0 * a * e_n + a * a)
            elif arch == "volatile_ejection/1":
                # Mirrors the generating solver's VE dispatch (polymer.pyx
                # section 5; keep in sync). The gas volatile itself flows
                # through the standard step-5 species path -- NO gas moles
                # are written here (single-count mass booking, ruling
                # round 20 increment 5).
                src, dst = e["src"], e["dst"]
                a = float(e["eject_moment"])
                if src and dst and src != dst:
                    # Cross-pool: from-leg loses the FULL bundle; to-leg
                    # gains the a-SHIFTED bundle (sa = -a forward, +a
                    # reverse; mu2 shift is the exact quadratic
                    # 2*sa*E[n] + sa^2, sa*sa == a*a both directions).
                    for ev, frm, to, sa in ((rf, src, dst, -a),
                                            (rr, dst, src, +a)):
                        if ev <= 0.0:
                            continue
                        ev_mol = ev * V_rxn
                        b0, b1, b2, ok = self._chain_bundle(
                            frm, y, e["scaling"] == "mu0")
                        if b0 == 0.0:
                            continue
                        f = self.pools[frm]["mu"]
                        t = self.pools[to]["mu"]
                        dn[f[0]] -= ev_mol * b0
                        dn[f[1]] -= ev_mol * b1
                        dn[t[0]] += ev_mol * b0
                        dn[t[1]] += ev_mol * (b1 + sa * b0)
                        if ok:
                            dn[f[2]] -= ev_mol * b2
                            dn[t[2]] += ev_mol * (
                                b2 + 2.0 * sa * b1 + a * a * b0)
                elif src and src == dst:
                    # Same-pool (src == dst fold-back): chip-style signed
                    # single-pool write; mu0 untouched. Forward mu2
                    # decrement carries the DISCRETE_CHIP `> 0` clamp --
                    # for a < 0 it is ALWAYS negative, so the forward mu2
                    # growth term is dropped (the exact +extension lives
                    # on the reverse leg), the documented a<0 convention.
                    b0, b1, _b2, _ok = self._chain_bundle(
                        src, y, e["scaling"] == "mu0")
                    if b0 != 0.0:
                        s = self.pools[src]["mu"]
                        e_n = b1
                        if rf > 0.0:
                            rf_mol = rf * V_rxn
                            dn[s[1]] -= rf_mol * a
                            dmu2 = 2.0 * a * e_n - a * a
                            if dmu2 > 0.0:
                                dn[s[2]] -= rf_mol * dmu2
                        if rr > 0.0:
                            rr_mol = rr * V_rxn
                            dn[s[1]] += rr_mol * a
                            dn[s[2]] += rr_mol * (2.0 * a * e_n + a * a)
            elif arch == "legacy_mu1/1":
                for pool in e["proxy_r_pools"]:
                    dn[self.pools[pool]["mu"][1]] -= r_mol
                for pool in e["proxy_p_pools"]:
                    dn[self.pools[pool]["mu"][1]] += r_mol
            elif arch == "moment_credit_conduit/1":
                # M18.3 step-6 bundle (DESIGN §2.2): forward-only r = rf
                # (kb = 0 structurally -- the row is irreversible by
                # contract, and any nonzero reverse contribution is
                # impossible-by-construction); point-mass credit
                # (1, u, u^2) * ev_mol to the destination pool ONLY. NO
                # site scaling (src_pool is null, so step 4 never ran), NO
                # exhaustion-tail limiter, NO cone-margin drain gate: the
                # conduit never debits pool moments, and its consumed
                # species is a real amount clamped >= 0 by the ordinary
                # step-3 concentration clamp. The credited bundle with
                # u >= 1 is inside the realizability cone by construction.
                if rf > 0.0:
                    ev_mol = rf * V_rxn
                    u = e["conduit_u"]
                    d = self.pools[e["conduit_dst"]]["mu"]
                    dn[d[0]] += ev_mol
                    dn[d[1]] += ev_mol * u
                    dn[d[2]] += ev_mol * u * u
            # same_pool/1: no moment flux

        # step 7a: channels (format doc §5)
        for pool in self.pools.values():
            i0, i1, i2 = pool["mu"]
            mu0 = max(0.0, y[i0]) / Vp
            mu1 = max(0.0, y[i1]) / Vp
            mu2 = max(0.0, y[i2]) / Vp
            dmu0 = dmu1 = dmu2 = 0.0
            if pool["k_s"] > 0.0:
                mu3 = safe_mu3(mu0, mu1, mu2)
                dmu0 += pool["k_s"] * (mu1 - mu0)
                if np.isfinite(mu3):
                    dmu2 += pool["k_s"] * (mu1 - mu3) / 3.0
            if pool["k_u"] > 0.0:
                r_ev = pool["k_u"] * mu0
                dmu1 -= r_ev
                dmu2 -= pool["k_u"] * (2.0 * mu1 - mu0)
                if pool["routing"] is not None:
                    dn[pool["routing"]] += r_ev * Vp
            dn[i0] += dmu0 * Vp
            dn[i1] += dmu1 * Vp
            dn[i2] += dmu2 * Vp

        # step 7b: mass transfer (format doc §7)
        for ig, ip, K, kLa in self.mass_transfer:
            J = kLa * (C[ip] - K * C[ig])
            dq = J * Vp
            dn[ig] += dq
            dn[ip] -= dq

        return dn

    def integrate_euler(self, y0, T, dt, n_steps, record_every=1):
        """Fixed-step forward Euler; returns (times, trajectory[n_rec, n_spc])."""
        y = np.array(y0, dtype=float)
        times, traj = [0.0], [y.copy()]
        for k in range(n_steps):
            y = y + dt * self.rhs(y, T)
            if (k + 1) % record_every == 0:
                times.append((k + 1) * dt)
                traj.append(y.copy())
        return np.array(times), np.array(traj)


def _label_lookup(nasa, label):
    """NASA data may be keyed by full chem.yaml label or its base form."""
    if label in nasa:
        return label
    base = _base(label)
    if base in nasa:
        return base
    raise KeyError(f"no NASA thermo for {label!r}")
