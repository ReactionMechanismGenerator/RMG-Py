"""Census-only refusal classifier for the moment-credit conduit
(item-20, increment M18.2; adjudicated rounds 32/36).

Classifies refused reaction rows (r93 general-branch refusals and upstream
FEATURE-RADICAL refusals) into admission buckets for the future surrogate
archetype ``moment_credit_conduit/1``. First increment (round-32
adjudication): shapes A and B only, admission is IRREVERSIBLE-ONLY, gas
products must satisfy the gas-phase MW threshold (<= 1.5 x monomer MW).

CENSUS-ONLY CONTRACT (M18.2): nothing in this module changes behavior.
Refused rows stay refused, no row is stamped, no flux is applied, nothing
is emitted into the sidecar (no schema/format-doc change), and the solver
is untouched. The ONLY visible effect is APPEND-ONLY annotation of the two
existing refusal census warning lines in :mod:`rmgpy.polymer` (the r93
general branch and the FEATURE-RADICAL census) plus this module's own
census summary/divergence warnings.

The module has two layers:

* a PURE CORE (no RMG imports, no I/O, deterministic): record-dict in,
  classification dict out. The record contract is::

      record = {
          "census": "r93_general" | "feature_radical",
          "reaction": str,                # as written in the census line
          "reversible": bool,             # "<=>" vs "=>"
          "reactants": [species, ...],
          "products":  [species, ...],
          "log_reason": str,              # refusal text from the census line
      }
      species = {
          "token": str,        # e.g. "CH3(8)" or "<smiles>(8158)"
          "label": str,        # token without the trailing (index)
          "index": int|None,
          "formula": str|None,
          "mw": float|None,          # g/mol
          "heavy_atoms": int|None,
          "link_marker": bool,       # SMILES carried CL/CR link pseudo-atoms
      }

* a THIN ADAPTER (``record_from_reaction``) that builds records from live
  RMG objects, REUSING the in-repo predicates
  (:func:`rmgpy.polymer._discrete_is_chain_scale_proxy_derived`,
  :func:`rmgpy.polymer._discrete_resolves_to_pool_state`) instead of
  re-deriving MW/heavy-atom thresholds, and reading the gas MW threshold
  from the ROW'S OWN pool participants (GAS_MW_FACTOR x min monomer MW
  over the row's pools), never from a module constant. RMG imports happen
  lazily inside the adapter so the core stays pure.

Buckets:
    ADMISSIBLE_A / ADMISSIBLE_B   shape A/B rows passing all admission gates
    DEFERRED_A / DEFERRED_B       shape A/B rows failing a gate (reason set)
    DEFERRED_C .. DEFERRED_F      shapes deferred by design in this increment
    FEATURE_RADICAL               upstream feature-radical census (deferred)
    UNCLASSIFIED                  shape not in the A-F vocabulary (reason set)

Shape vocabulary (roles per side, direction-normalized so the pool-crediting
side is the product side):
    A: CHAIN            <=> POOL + DISC
    B: CHAIN + DISC     <=> POOL + DISC
    C: CHAIN            <=> POOL + POOL
    D: CHAIN + CHAIN    <=> POOL + POOL
    E: CHAIN + POOL     <=> POOL + POOL     (POOL reactant = chip resolvable
                                             to pool state, e.g. trimer_rad*)
    F: CHAIN + DISC     <=> POOL            (gas association, no gas product)

Note on D/E vs the round-32 shape strings: round 32 wrote D as "chain-disc +
disc <=> pool + pool" and E as "pool + disc <=> pool + pool". Empirically
(poly_102_conduit3) every D/E reactant pair is two chain-SCALE species; the
operative D/E split is whether one of them is pool-state-resolvable (the
trimer_rad* chips -> E, 54 rows) or not (-> D, 270 rows). The tables here
encode the empirical split, which reproduces the adjudicated counts exactly.

Round-36 landing P1s carried by this module:

(a) STABLE CANDIDATE KEY + OVERLAP PRECEDENCE. 8,561 rows appear in BOTH
    the r93 and feature-radical censuses. :func:`conduit_candidate_key`
    gives every row a deterministic, orientation- and census-independent
    identity; the process-wide registry (:func:`register_candidate`)
    records which censuses have seen each key. FEATURE-RADICAL refusal
    WINS as the upstream blocker until admission exists: an overlapped
    candidate's ``effective_bucket`` is ``FEATURE_RADICAL`` regardless of
    its r93-side classification (kept as ``shadow_bucket`` so M18.3 can
    see what it WOULD be) -- M18.3 must consult the registry and never
    double-book a key.

(b) APPEND-ONLY WARNING EDITS. This module only ever produces SUFFIX
    strings (:func:`census_suffix`) appended after the existing warning
    text; existing headers / refusal tokens / reason strings are never
    renamed or restructured (downstream greps are assumed).

(c) LOUD CENSUS. :func:`census_summary` ALWAYS surfaces the UNCLASSIFIED
    count, including when it is zero; every appended suffix carries the
    running ``unclassified=<n>`` counter. Label-test vs isomorphism-test
    divergence in the adapter is census-logged (its own warning token,
    ``CONDUIT CLASSIFIER DIVERGENCE``) and recorded on the result --
    the in-repo isomorphism verdict is used, never a silent override of
    either side.
"""

import logging
import math
import threading
from collections import Counter
from dataclasses import dataclass
from typing import Optional, Tuple

# ---------------------------------------------------------------------------
# Configuration (round-32 adjudicated constants)
# ---------------------------------------------------------------------------

#: Repeat-unit (monomer) of the phenol-formaldehyde proxy. The deck comment
#: calls it "C8H9O" but the deck's own SMILES is C9H10O — the SMILES is trusted.
MONOMER_FORMULA = "C9H10O"
MONOMER_MW = 134.178          # g/mol  (9*12.011 + 10*1.008 + 15.999)
MONOMER_HEAVY_ATOMS = 10

#: Chain-scale test (mirrors rmgpy/polymer.py refusal logic): a discrete is
#: chain-scale iff it reaches >= 2.5 monomer-equivalents on BOTH the mass and
#: heavy-atom axes.
CHAIN_SCALE_FACTOR = 2.5
CHAIN_SCALE_MW = CHAIN_SCALE_FACTOR * MONOMER_MW              # 335.445
CHAIN_SCALE_HEAVY = CHAIN_SCALE_FACTOR * MONOMER_HEAVY_ATOMS  # 25.0

#: Gas-phase species threshold: MW <= ~1.5 x monomer MW. NOTE: this module
#: default exists for the pure core / offline census only; the in-repo
#: adapter derives the threshold from the ROW'S OWN pool participants
#: (see :func:`gas_mw_threshold_for_pools`).
GAS_MW_FACTOR = 1.5
GAS_MW_THRESHOLD = GAS_MW_FACTOR * MONOMER_MW                 # 201.267 g/mol

#: Near-threshold band for census reporting (+-10% of the threshold).
NEAR_THRESHOLD_BAND = 0.10

#: Known pool-label prefix. Every polymer pool in poly_102_conduit3 (configured
#: pool + spawned mod pools) carries this prefix; spawned pools resolve to it.
POOL_LABEL_PREFIXES = ("phenol_formaldehyde",)

#: Condensed chip species that RESOLVE to a pool state (sidecar
#: condensed_species list, polymer_pools.json schema 2.5). The live machinery
#: decides this by isomorphism (_discrete_resolves_to_pool_state); census-side
#: we reconstruct it by label membership. These take the POOL role in shape
#: classification (this is what separates shape E from shape D).
POOL_STATE_RESOLVABLE_LABELS = ("trimer_rad33", "trimer_rad38", "trimer_rad44")

# Species roles
POOL = "POOL"
CHAIN = "CHAIN"
DISC = "DISC"
UNKNOWN = "UNKNOWN"

# Refusal-reason labels (for anything non-admissible)
REASON_SHAPE_DEFERRED = "shape-deferred-first-increment"      # C/D/E/F
REASON_GAS_MW_FAIL = "gas-product-exceeds-mw-threshold"       # A/B, product too heavy
REASON_FEATURE_RADICAL = "feature-radical-upstream-refusal"
REASON_DIRECTION = "direction-inadmissible-irreversible-source"
REASON_UNRESOLVED = "unresolved-species-data"
REASON_UNKNOWN_SHAPE = "shape-outside-vocabulary"

_SHAPE_TABLE = {
    (("CHAIN",), ("DISC", "POOL")): "A",
    (("CHAIN", "DISC"), ("DISC", "POOL")): "B",
    (("CHAIN",), ("POOL", "POOL")): "C",
    (("CHAIN", "CHAIN"), ("POOL", "POOL")): "D",
    (("CHAIN", "POOL"), ("POOL", "POOL")): "E",
    (("CHAIN", "DISC"), ("POOL",)): "F",
}

#: Overlap precedence rule (round-36 P1(a)): the upstream blocker wins.
PRECEDENCE_RULE = "feature_radical"


# ---------------------------------------------------------------------------
# Species-level classification (pure core)
# ---------------------------------------------------------------------------

def _role_label(species):
    """Label used for ROLE decisions. The optional ``label_for_roles``
    override lets the in-repo adapter fold an isomorphism verdict into the
    census label test without mutating the display label (round-36 P1(c):
    the isomorphism test is the in-repo source of truth; the divergence is
    census-logged, never silent)."""
    lab = species.get("label_for_roles")
    if lab is None:
        lab = species.get("label")
    return lab or ""


def is_pool(species):
    """True if the species is a polymer pool participant (by label)."""
    label = _role_label(species)
    return any(label.startswith(p) for p in POOL_LABEL_PREFIXES)


def is_pool_state_resolvable(species):
    """True for condensed chips that resolve to a pool state (by label)."""
    return _role_label(species) in POOL_STATE_RESOLVABLE_LABELS


def is_chain_scale(species):
    """Chain-scale discrete: >= 2.5 monomer-equiv on mass AND heavy atoms."""
    mw = species.get("mw")
    heavy = species.get("heavy_atoms")
    if mw is None or heavy is None:
        return False
    return mw >= CHAIN_SCALE_MW and heavy >= CHAIN_SCALE_HEAVY


def species_role(species):
    if is_pool(species) or is_pool_state_resolvable(species):
        return POOL
    if species.get("mw") is None:
        return UNKNOWN
    if is_chain_scale(species):
        return CHAIN
    return DISC


def gas_mw_ok(species, threshold=GAS_MW_THRESHOLD):
    mw = species.get("mw")
    return mw is not None and mw <= threshold


def near_threshold(species, threshold=GAS_MW_THRESHOLD, band=NEAR_THRESHOLD_BAND):
    mw = species.get("mw")
    if mw is None:
        return False
    return abs(mw - threshold) <= band * threshold


# ---------------------------------------------------------------------------
# Row-level classification (pure core)
# ---------------------------------------------------------------------------

def _side_roles(side):
    return tuple(sorted(species_role(s) for s in side))


def _resolve_shape(record):
    """Return (shape, admitted_direction, chain_side, pool_side).

    Direction is normalized so the pool-crediting side is the product side:
    the side holding the CHAIN (or, for chain-free shape E, the single-pool
    side) is the consumed side. admitted_direction is "forward" when the
    normalized orientation matches the row as written, "reverse" otherwise.
    """
    r_roles = _side_roles(record["reactants"])
    p_roles = _side_roles(record["products"])
    if (r_roles, p_roles) in _SHAPE_TABLE:
        return _SHAPE_TABLE[(r_roles, p_roles)], "forward", record["reactants"], record["products"]
    if (p_roles, r_roles) in _SHAPE_TABLE:
        return _SHAPE_TABLE[(p_roles, r_roles)], "reverse", record["products"], record["reactants"]
    return None, None, None, None


def _destination_pool(consumed_side, produced_side):
    """Pool the chain-scale species resolves to: the unique pool participant
    on the produced (pool-crediting) side. None when ambiguous (2-pool rows)."""
    pools = [s["label"] for s in produced_side if is_pool(s)]
    if len(pools) == 1:
        return pools[0]
    return None


def classify_record(record, gas_mw_threshold=GAS_MW_THRESHOLD):
    """Classify one refusal record. Returns a new dict (input not mutated)."""
    all_species = list(record["reactants"]) + list(record["products"])
    has_unresolved = any(s.get("mw") is None for s in all_species)
    has_link_marker = any(s.get("link_marker") for s in all_species)

    shape, direction, consumed, produced = _resolve_shape(record)

    result = {
        "census": record["census"],
        "reaction": record["reaction"],
        "candidate_key": conduit_candidate_key(record),
        "shape": shape,
        "bucket": None,
        "admissible": False,
        "admitted_direction": None,
        "irreversible_only": True,          # admission mode is always one-way
        "gas_mw_threshold": gas_mw_threshold,
        "gas_products": [],                 # non-pool species on produced side
        "gas_products_mw_ok": None,
        "gas_products_near_threshold": [],
        "gas_reactants_over_threshold": [], # informational flag (not a gate)
        "destination_pool": None,
        "refusal_reason": None,
        "flags": [],
    }
    if has_link_marker:
        result["flags"].append("link-marker-species")

    # Feature-radical census: entirely deferred in this increment.
    if record["census"] == "feature_radical":
        result["bucket"] = "FEATURE_RADICAL"
        result["refusal_reason"] = REASON_FEATURE_RADICAL
        return result

    if has_unresolved:
        result["bucket"] = "UNCLASSIFIED"
        result["refusal_reason"] = REASON_UNRESOLVED
        return result

    if shape is None:
        result["bucket"] = "UNCLASSIFIED"
        result["refusal_reason"] = REASON_UNKNOWN_SHAPE
        return result

    # Direction: admitted direction consumes the chain and credits the pool.
    result["admitted_direction"] = direction
    result["destination_pool"] = _destination_pool(consumed, produced)

    gas_products = [s for s in produced if not is_pool(s)]
    result["gas_products"] = [s["token"] for s in gas_products]
    result["gas_products_mw_ok"] = all(gas_mw_ok(s, gas_mw_threshold) for s in gas_products)
    result["gas_products_near_threshold"] = [
        s["token"] for s in gas_products if near_threshold(s, gas_mw_threshold)
    ]
    result["gas_reactants_over_threshold"] = [
        s["token"] for s in consumed
        if not is_pool(s) and not is_chain_scale(s) and not gas_mw_ok(s, gas_mw_threshold)
    ]

    if shape not in ("A", "B"):
        result["bucket"] = f"DEFERRED_{shape}"
        result["refusal_reason"] = REASON_SHAPE_DEFERRED
        return result

    # Shape A/B admission gates.
    # Gate 1 (direction): admission is irreversible-only in the admitted
    # direction. A reversible source row is admitted one-way. An irreversible
    # source row whose written direction opposes the admitted direction cannot
    # be admitted at all.
    if not record["reversible"] and direction == "reverse":
        result["bucket"] = f"DEFERRED_{shape}"
        result["refusal_reason"] = REASON_DIRECTION
        return result

    # Gate 2 (gas MW): every discrete product entering the gas phase must sit
    # at or below the gas-phase threshold.
    if not result["gas_products_mw_ok"]:
        result["bucket"] = f"DEFERRED_{shape}"
        result["refusal_reason"] = REASON_GAS_MW_FAIL
        return result

    result["bucket"] = f"ADMISSIBLE_{shape}"
    result["admissible"] = True
    return result


def classify_all(records, gas_mw_threshold=GAS_MW_THRESHOLD):
    """Classify an iterable of records. Returns a list of result dicts."""
    return [classify_record(r, gas_mw_threshold) for r in records]


def bucket_counts(results):
    """Counter of bucket labels over classification results. ALWAYS carries
    an UNCLASSIFIED entry (round-36 P1(c) loud-census rule), zero or not."""
    counts = Counter(r["bucket"] for r in results)
    counts.setdefault("UNCLASSIFIED", 0)
    return counts


def refusal_reason_counts(results):
    """Counter of refusal reasons over non-admissible results."""
    return Counter(r["refusal_reason"] for r in results if not r["admissible"])


# ---------------------------------------------------------------------------
# Round-36 P1(a): stable candidate key + overlap/precedence registry
# ---------------------------------------------------------------------------

def candidate_key_from_label(reaction_label):
    """Deterministic, orientation- and census-independent identity of a row
    from its census equation string.

    Sides are split on '<=>' / '=>', tokens on ' + ', each side's tokens are
    sorted, and the two sides are ordered lexicographically -- so the SAME
    row seen written in either orientation, from either census, in any run,
    maps to ONE key. The arrow (reversibility) is deliberately NOT part of
    the key: it is a property of the row, not of its identity, and the
    overlap ledger must join r93/feature-radical sightings of one row."""
    label = (reaction_label or "").strip()
    for arrow in ("<=>", "=>", "<="):
        if arrow in label:
            left, right = label.split(arrow, 1)
            break
    else:
        left, right = label, ""
    side_a = "+".join(sorted(t.strip() for t in left.split(" + ") if t.strip()))
    side_b = "+".join(sorted(t.strip() for t in right.split(" + ") if t.strip()))
    lo, hi = sorted((side_a, side_b))
    return f"{lo}<>{hi}"


def conduit_candidate_key(record):
    """Candidate key for a record dict (delegates to the label form so the
    string-only feature-radical hook and the full r93 adapter hook produce
    IDENTICAL keys for the same row)."""
    return candidate_key_from_label(record.get("reaction") or "")


class _CensusRegistry:
    """Process-wide, thread-safe candidate ledger (census-only bookkeeping).

    One entry per candidate key:
        {"censuses": set, "bucket_by_census": dict,
         "effective_bucket": str, "shadow_bucket": str|None,
         "precedence": "feature_radical"|None}

    FEATURE-RADICAL WINS (round-36 P1(a)): once a key has been seen in the
    feature-radical census, its effective bucket is FEATURE_RADICAL (the
    upstream blocker) no matter what the r93 classifier says; the r93
    verdict is preserved as ``shadow_bucket`` for M18.3. M18.3 must call
    :meth:`lookup` before booking a candidate and skip keys already booked.
    """

    def __init__(self):
        self._lock = threading.Lock()
        self._entries = {}
        self.unclassified_total = 0

    def register(self, key, census, bucket, epoch=None):
        with self._lock:
            entry = self._entries.setdefault(
                key, {"censuses": set(), "bucket_by_census": {},
                      "effective_bucket": None, "shadow_bucket": None,
                      "precedence": None, "epochs": set()})
            first_sighting = census not in entry["censuses"]
            entry["censuses"].add(census)
            # M18.3 ledger epochs (DESIGN §3.3): the engine rebuild-signature
            # id (or RMG iteration int) under which this sighting was
            # recorded. Bookkeeping only -- FR stickiness stays set-membership
            # of "feature_radical" in censuses, across ALL epochs of the run.
            if epoch is not None:
                entry["epochs"].add(epoch)
            entry["bucket_by_census"][census] = bucket
            if bucket == "UNCLASSIFIED" and first_sighting:
                self.unclassified_total += 1
            if "feature_radical" in entry["censuses"]:
                entry["effective_bucket"] = "FEATURE_RADICAL"
                entry["shadow_bucket"] = entry["bucket_by_census"].get(
                    "r93_general")
                entry["precedence"] = (
                    PRECEDENCE_RULE if len(entry["censuses"]) > 1 else None)
            else:
                entry["effective_bucket"] = bucket
                entry["shadow_bucket"] = None
                entry["precedence"] = None
            return dict(entry, censuses=set(entry["censuses"]),
                        epochs=set(entry["epochs"]))

    def lookup(self, key):
        with self._lock:
            entry = self._entries.get(key)
            if entry is None:
                return None
            return dict(entry, censuses=set(entry["censuses"]),
                        epochs=set(entry["epochs"]))

    def counts(self):
        with self._lock:
            eff = Counter(e["effective_bucket"] for e in self._entries.values())
            eff.setdefault("UNCLASSIFIED", 0)
            overlap = sum(1 for e in self._entries.values()
                          if len(e["censuses"]) > 1)
            return eff, overlap, len(self._entries)

    def reset(self):
        with self._lock:
            self._entries.clear()
            self.unclassified_total = 0


#: The process-wide registry instance.
CENSUS_REGISTRY = _CensusRegistry()


def register_candidate(key, census, bucket, epoch=None):
    """Record a candidate sighting; returns the resolved registry entry."""
    return CENSUS_REGISTRY.register(key, census, bucket, epoch=epoch)


def lookup_candidate(key):
    """M18.3 entry point: the resolved entry for a key, or None."""
    return CENSUS_REGISTRY.lookup(key)


def reset_conduit_state():
    """Run-boundary HARD reset (DESIGN §3.3, called once from RMG
    initialization -- rmgpy/rmg/main.py RMG.initialize): clears the
    process-wide candidate ledger AND the warn-once sets that feed it,
    together ("reset both or neither"). The FEATURE-RADICAL census is
    warn-once per (reaction, reason) (rmgpy/polymer.py
    _refused_census_warned); resetting the ledger without that set would
    silently starve the FR side of re-sightings in the new run, so a
    fresh run could admit a candidate its own censuses never got to
    re-block. Also clears the run-level conduit flux accumulator (§4.4).
    Import-safe: the rmgpy.polymer sets are cleared lazily so this module
    keeps its no-module-level-RMG-imports contract."""
    CENSUS_REGISTRY.reset()
    _CONDUIT_FLUX_TOTALS.clear()
    try:
        import rmgpy.polymer as _polymer
        _polymer._refused_census_warned.clear()
        # Conduit keys of the archetype warn-once set (none exist in
        # M18.3 -- no conduit census line is keyed there yet -- but the
        # both-or-neither rule is cheap to keep exhaustive).
        _polymer._flux_archetype_warned.difference_update(
            {k for k in _polymer._flux_archetype_warned
             if "conduit" in str(k).lower()})
    except ImportError:  # pragma: no cover - pure-core usage
        pass


def reset_census_registry():
    """Historical test hook, kept as an ALIAS of :func:`reset_conduit_state`
    (M18.3 W1.5): the ledger and its warn-once feeder sets reset together
    or not at all (DESIGN §3.3)."""
    reset_conduit_state()


def census_summary():
    """One-line loud census summary. ALWAYS names UNCLASSIFIED, zero or not."""
    eff, overlap, total = CENSUS_REGISTRY.counts()
    parts = " ".join(f"{b}={eff[b]}" for b in sorted(eff))
    return (f"conduit-census/1 summary: candidates={total} overlap={overlap} "
            f"{parts}")


def log_census_summary():
    """Census-log the summary (warning level, like the refusal censuses)."""
    logging.warning("MOMENT-CREDIT CONDUIT CENSUS (M18.2, census-only): %s",
                    census_summary())


def census_suffix(result, entry):
    """APPEND-ONLY structured suffix for an existing refusal warning line
    (round-36 P1(b)): never a replacement of any existing header/token --
    callers append it verbatim after their unchanged message text."""
    overlap = "+".join(sorted(entry["censuses"]))
    precedence = entry["precedence"] or "none"
    shadow = (f" shadow_bucket={entry['shadow_bucket']}"
              if entry["shadow_bucket"] else "")
    return (f" [conduit-census/1 key={result['candidate_key']} "
            f"bucket={entry['effective_bucket']}{shadow} "
            f"shape={result['shape']} admissible={result['admissible']} "
            f"censuses={overlap} precedence={precedence} "
            f"unclassified={CENSUS_REGISTRY.unclassified_total}]")


# ---------------------------------------------------------------------------
# Thin adapter over live RMG objects (lazy imports; census-only)
# ---------------------------------------------------------------------------

def gas_mw_threshold_for_pools(row_pools):
    """Gas-phase MW threshold from the ROW'S OWN pool participants:
    GAS_MW_FACTOR x the smallest available monomer MW among the row's pools
    (conservative for multi-pool rows: the tightest pool gates). Falls back
    to the module default when no pool carries a usable monomer MW."""
    mws = []
    for poly in row_pools or []:
        mw = float(getattr(poly, "monomer_mw_g_mol", 0.0) or 0.0)
        if mw > 0.0:
            mws.append(mw)
    if not mws:
        return GAS_MW_THRESHOLD
    return GAS_MW_FACTOR * min(mws)


def _species_entry(species, row_pools):
    """Species record from a live RMG object, reusing the in-repo predicates.

    Round-36 P1(c) divergence rule: the census-side LABEL test for
    pool-state resolvability is compared against the in-repo ISOMORPHISM
    test (:func:`rmgpy.polymer._discrete_resolves_to_pool_state`); on
    divergence the mismatch is census-logged (token: CONDUIT CLASSIFIER
    DIVERGENCE) and flagged on the entry, and the ISOMORPHISM verdict is
    used -- neither side is silently overridden.
    """
    from rmgpy.molecule import Molecule
    from rmgpy.polymer import (Polymer, _discrete_resolves_to_pool_state,
                               _heavy_atom_count)

    label = getattr(species, "label", "") or ""
    index = getattr(species, "index", None)
    if isinstance(index, int) and index >= 0:
        token = f"{label}({index})"
    else:
        token = label
        index = None
    entry = {"token": token, "label": label, "index": index,
             "formula": None, "mw": None, "heavy_atoms": None,
             "link_marker": False, "_divergence": False}

    if isinstance(species, Polymer):
        # Pool participants classify by label prefix in the core; carry the
        # monomer identity for completeness.
        entry["mw"] = float(getattr(species, "monomer_mw_g_mol", 0.0) or 0.0) or None
        return entry, species

    mol_list = getattr(species, "molecule", None)
    mol = mol_list[0] if mol_list else (
        species if isinstance(species, Molecule) else None)
    if mol is not None:
        try:
            entry["mw"] = mol.get_molecular_weight() * 1000.0
            entry["heavy_atoms"] = _heavy_atom_count(mol)
            entry["formula"] = mol.get_formula()
        except (ValueError, AttributeError):
            pass

    # Divergence check: census label test vs in-repo isomorphism test.
    label_says = label in POOL_STATE_RESOLVABLE_LABELS
    try:
        iso_says = bool(_discrete_resolves_to_pool_state(species, row_pools))
    except Exception:  # defensive: census-only code must not raise
        iso_says = False
    if label_says != iso_says:
        entry["_divergence"] = True
        logging.warning(
            "CONDUIT CLASSIFIER DIVERGENCE (M18.2 census-only): species %s "
            "pool-state resolvability disagrees between the census label "
            "test (%s) and the in-repo isomorphism test (%s); using the "
            "isomorphism verdict. This is a finding, not an override -- "
            "record it for the next adversarial round.",
            token, label_says, iso_says)
    if iso_says and not label_says:
        # Make the core's role assignment follow the isomorphism verdict
        # without touching POOL_STATE_RESOLVABLE_LABELS: mark via label
        # override field the core checks first.
        entry["_iso_pool_state"] = True
    if label_says and not iso_says:
        entry["_iso_pool_state"] = False
    return entry, None


def record_from_reaction(forward, row_pools, census="r93_general"):
    """Build a classifier record from a live (refused) RMG reaction.

    ``row_pools`` is the row's own Polymer participants (exactly what the
    r93 branch computes). The gas MW threshold consumers should pass to
    :func:`classify_record` comes from :func:`gas_mw_threshold_for_pools`.
    """
    from rmgpy.polymer import _reaction_census_label

    reactants = list(getattr(forward, "reactants", None) or [])
    products = list(getattr(forward, "products", None) or [])
    r_entries, p_entries = [], []
    divergent = False
    for side_in, side_out in ((reactants, r_entries), (products, p_entries)):
        for s in side_in:
            entry, _poly = _species_entry(s, row_pools)
            divergent = divergent or entry.pop("_divergence", False)
            side_out.append(entry)
    record = {
        "census": census,
        "reaction": _reaction_census_label(forward),
        "reversible": bool(getattr(forward, "reversible", True)),
        "reactants": r_entries,
        "products": p_entries,
        "log_reason": "",
        "label_isomorphism_divergence": divergent,
    }
    return record


def _apply_iso_overrides(record):
    """Fold the adapter's isomorphism verdicts into the core's label test:
    entries carrying _iso_pool_state get their label swapped onto/off the
    resolvable list surrogate (handled by giving them a role hint)."""
    for side in (record["reactants"], record["products"]):
        for s in side:
            hint = s.pop("_iso_pool_state", None)
            if hint is True and s["label"] not in POOL_STATE_RESOLVABLE_LABELS:
                s["label_for_roles"] = POOL_STATE_RESOLVABLE_LABELS[0]
            elif hint is False and s["label"] in POOL_STATE_RESOLVABLE_LABELS:
                s["label_for_roles"] = ""
    return record


def annotate_refused_row(forward, row_pools, census="r93_general",
                         verdict=None, epoch=None):
    """CENSUS-ONLY annotation hook for the refusal warning sites in
    :mod:`rmgpy.polymer`. Classifies the refused row, registers it in the
    candidate ledger, and returns the APPEND-ONLY suffix for the existing
    warning line. NEVER raises (a census annotation must not be able to
    change generation behavior): on any internal error it census-logs the
    failure loudly and returns an empty suffix.

    M18.3: when the caller evaluated admission (the r93 stamp site), it
    passes the :class:`AdmissionVerdict` so the suffix carries the would-be
    verdict while :data:`CONDUIT_ADMISSION_ENABLED` is False
    (``would_admit=1 ...`` / ``deny=<reason>``, append-only)."""
    try:
        record = _apply_iso_overrides(
            record_from_reaction(forward, row_pools, census=census))
        threshold = gas_mw_threshold_for_pools(row_pools)
        result = classify_record(record, gas_mw_threshold=threshold)
        if record.get("label_isomorphism_divergence"):
            result["flags"].append("label-isomorphism-divergence")
        entry = register_candidate(result["candidate_key"], census,
                                   result["bucket"], epoch=epoch)
        # r42 P1-4(b) ordering pin: the caller's admission verdict was
        # evaluated BEFORE this row's sighting was registered, so its G1
        # ledger consult can miss a feature-radical sighting that landed in
        # between (the FR census is warn-once and runs in either order
        # relative to the r93 stamp site). Re-check against the
        # POST-registration entry so the emitted line is deterministic
        # w.r.t. ledger state that includes the current row: an
        # FR-overlapped key NEVER prints would_admit=1.
        if (verdict is not None and verdict.admitted
                and "feature_radical" in entry["censuses"]):
            verdict = _deny(verdict.candidate_key or result["candidate_key"],
                            "feature-radical-overlap")
        return census_suffix(result, entry) + admission_census_suffix(verdict)
    except Exception as exc:  # pragma: no cover - defensive fail-open
        logging.warning(
            "MOMENT-CREDIT CONDUIT CENSUS (M18.2): annotation failed for a "
            "refused row (%s: %s); the refusal itself is unaffected "
            "(census-only code path).", type(exc).__name__, exc)
        return ""


# ---------------------------------------------------------------------------
# M18.3: admission layer (FAIL-CLOSED DEAD CODE behind the census)
# ---------------------------------------------------------------------------

#: MASTER admission switch. M18.3 lands the admission gates + stamping arm as
#: fail-closed DEAD CODE behind the census: every row still refuses exactly
#: as M18.2 left it while the census logs the would-be verdict (WOULD-ADMIT /
#: deny reason) for sizing. Flipping this to True is M18.4, gated on THREE
#: open items (BUILD_SPEC / OPEN_QUESTIONS):
#:   OQ-1  sizing of the deferred populations (FR-overlap +
#:         direction-requires-flip-rewrite denial counts),
#:   OQ-2  the .reversible in-place-mutation safety audit,
#:   OQ-10 TA's cert-path answer (sidecar rows vs classic-Cantera-only).
CONDUIT_ADMISSION_ENABLED = False

#: In-band caveat serialized beside every candidate_key (DESIGN §2.1):
#: keys are built from label(index) census strings whose indices do not
#: survive regeneration.
CANDIDATE_KEY_NOTE = ("run-scoped provenance only; species indices are not "
                      "stable across regenerations -- never join across "
                      "artifacts")

#: Closed deny-reason vocabulary (BUILD_SPEC W1.3; used by census + tests).
#: admission-evaluation-error:* is the G7 family (suffixed with the
#: exception type name).
ADMISSION_DENY_REASONS = frozenset({
    "classifier-not-admissible", "classifier-divergence",
    "feature-radical-overlap", "direction-inadmissible",
    "direction-requires-flip-rewrite", "gas-product-count",
    "gas-mw-threshold-unresolvable", "gas-mw-over-threshold",
    "landing-cone-violation", "destination-unresolvable",
    "chain-not-condensed", "not-balanced", "kinetics-not-exportable",
    "kinetics-not-yet-assigned",
})

#: PROVISIONAL deny reasons (G6 re-adjudication defect fix): the r93 stamp
#: site (rmgpy/rmg/model.py make_new_reaction, BEFORE check_existing) runs
#: before make_new_reaction assigns kinetics, so family-generated
#: TemplateReactions always arrive at G6 with ``kinetics is None`` -- that
#: is a sequencing artifact, not a property of the row's eventual kinetics.
#: Such rows deny provisionally (a provisional row NEVER prints
#: would_admit=1) and get G6 re-evaluated against the now-final kinetics by
#: :func:`rmgpy.polymer.readjudicate_conduit_admission` after the kinetics
#: conversion / barrier-correction block (and across the canonical-dedup
#: merge). Every other deny reason is FINAL.
PROVISIONAL_DENY_REASONS = frozenset({"kinetics-not-yet-assigned"})


@dataclass(frozen=True)
class AdmissionVerdict:
    """Result of :func:`evaluate_conduit_admission`. Default-DENY: every
    field starts at its refusing value; admission is built ONLY through
    affirmative gate passes."""
    admitted: bool = False
    deny_reason: Optional[str] = None
    chain_units: Optional[float] = None          # u (DESIGN §2.3)
    gas_product: Optional[Tuple[str, float]] = None  # (species_token, mw)
    gas_units: Optional[float] = None            # a = mw_gas / M(dst)
    dst_pool: Optional[str] = None
    candidate_key: str = ""
    needs_irreversible_rewrite: bool = False     # True for aligned-reversible


def _deny(key, reason):
    return AdmissionVerdict(admitted=False, deny_reason=reason,
                            candidate_key=key)


def evaluate_conduit_admission(forward, row_pools):
    """Admission gates G0-G7 (DESIGN §3.2) over a live RMG reaction at the
    r93 general-branch stamp site. FAIL-CLOSED: the default is DENY; any
    exception anywhere denies (G7) -- admission code never "recovers into
    admit". M18.3: the caller's admit arm is dead behind
    :data:`CONDUIT_ADMISSION_ENABLED`; this evaluator only feeds the census.
    """
    key = ""
    try:
        from rmgpy.polymer import (Polymer,
                                   _discrete_is_chain_scale_proxy_derived,
                                   strip_rmg_index_suffix)

        record = _apply_iso_overrides(
            record_from_reaction(forward, row_pools, census="r93_general"))
        key = conduit_candidate_key(record)

        # G0 (part 1) -- classifier divergence: a label/isomorphism
        # divergence is a finding, never an admission basis.
        if record.get("label_isomorphism_divergence"):
            return _deny(key, "classifier-divergence")

        # G1 -- ledger: an FR sighting in ANY epoch of this process blocks
        # admission for the whole run (sticky, fail-closed).
        entry = lookup_candidate(key)
        if entry is not None and "feature_radical" in entry["censuses"]:
            return _deny(key, "feature-radical-overlap")

        # Unresolved species data (no MW on a non-pool participant): the
        # classifier buckets such rows UNCLASSIFIED; uncertainty never
        # admits.
        r_side, p_side = record["reactants"], record["products"]
        if any(species_role(s) == UNKNOWN
               for s in list(r_side) + list(p_side)):
            return _deny(key, "classifier-not-admissible")

        # G2 -- direction/orientation [r39-P1]: the shape must resolve with
        # a UNIQUE chain_to_pool admitted direction (chain strictly on the
        # consumed side, pool participation strictly on the credited side).
        r_roles = [species_role(s) for s in r_side]
        p_roles = [species_role(s) for s in p_side]
        if (POOL in p_roles and POOL not in r_roles
                and CHAIN in r_roles and CHAIN not in p_roles):
            direction = "forward"
            consumed, produced = r_side, p_side
            consumed_objs = list(getattr(forward, "reactants", None) or [])
            produced_objs = list(getattr(forward, "products", None) or [])
        elif (POOL in r_roles and POOL not in p_roles
                and CHAIN in p_roles and CHAIN not in r_roles):
            direction = "reverse"
            consumed, produced = p_side, r_side
            consumed_objs = list(getattr(forward, "products", None) or [])
            produced_objs = list(getattr(forward, "reactants", None) or [])
        else:
            return _deny(key, "classifier-not-admissible")
        if direction == "reverse":
            # The admitted direction is the written-REVERSE orientation:
            # (i) an irreversible source written against it can never run
            # it; (ii) a reversible source would need FABRICATED reverse
            # kinetics (an Arrhenius fit of kf/Keq(T)) -- real rewriting
            # with real fit error -- deferred until that rewrite is
            # implemented and adjudicated (v2). Counted separately so the
            # deferred population can be sized (OQ-1).
            if record["reversible"]:
                return _deny(key, "direction-requires-flip-rewrite")
            return _deny(key, "direction-inadmissible")

        # G3 -- gas product: EXACTLY ONE non-pool species on the credited
        # side, with stoichiometric multiplicity 1 [r39-P5] (tokens counted
        # WITH repeats: a stoich-2 product appears twice).
        gas_entries = [s for s in produced if species_role(s) != POOL]
        if len(gas_entries) != 1:
            return _deny(key, "gas-product-count")
        gas = gas_entries[0]

        # G3 threshold -- derived from the ROW'S OWN pools; the census-time
        # module-default fallback (gas_mw_threshold_for_pools) is FORBIDDEN
        # on the admission path: no usable pool MW -> DENY.
        pool_mws = []
        for poly in row_pools or []:
            mw = float(getattr(poly, "monomer_mw_g_mol", 0.0) or 0.0)
            if mw > 0.0:
                pool_mws.append(mw)
        if not pool_mws:
            return _deny(key, "gas-mw-threshold-unresolvable")
        threshold = GAS_MW_FACTOR * min(pool_mws)
        gas_mw = gas.get("mw")
        if gas_mw is None:
            return _deny(key, "gas-mw-threshold-unresolvable")
        if gas_mw > threshold:
            return _deny(key, "gas-mw-over-threshold")

        # G0 (part 2) -- the classifier's own verdict, with the admission
        # threshold: nothing broader than shapes A/B is admitted.
        result = classify_record(record, gas_mw_threshold=threshold)
        if result["bucket"] not in ("ADMISSIBLE_A", "ADMISSIBLE_B"):
            return _deny(key, "classifier-not-admissible")

        # G5 -- destination: exactly one pool participant on the credited
        # side, resolving to a solver-configured pool with a usable monomer
        # MW. v1 keeps chip-resolved (pool-state) destinations deferred:
        # only a real Polymer participant is a resolvable destination.
        dst_polys = [s for s in produced_objs if isinstance(s, Polymer)]
        if len(dst_polys) != 1:
            return _deny(key, "destination-unresolvable")
        dst_poly = dst_polys[0]
        monomer_mw = float(getattr(dst_poly, "monomer_mw_g_mol", 0.0) or 0.0)
        if monomer_mw <= 0.0:
            return _deny(key, "destination-unresolvable")
        dst_label = strip_rmg_index_suffix(
            str(getattr(dst_poly, "label", "")))
        if not dst_label:
            return _deny(key, "destination-unresolvable")

        # G5 (chain phase): the chain species must be melt-classified
        # (chain-scale proxy-derived) so the consumer's step-2 phase gate
        # passes the event (V_rxn = V_poly).
        chain_positions = [i for i, s in enumerate(consumed)
                           if species_role(s) == CHAIN]
        if len(chain_positions) != 1:
            return _deny(key, "classifier-not-admissible")
        chain_obj = consumed_objs[chain_positions[0]]
        if not _discrete_is_chain_scale_proxy_derived(chain_obj, row_pools):
            return _deny(key, "chain-not-condensed")

        # G4 -- landing cone (DESIGN §2.3, exact MWs; closed >= 1.0
        # semantics, equality-boundary rider T6 / OQ-4):
        #   u_raw = (MW(chain) + sum MW(disc reactants)) / M
        #   a     = mw_gas / M
        #   u     = u_raw - a + d/M
        # guard: u >= 1.0  (per-event surplus mu1 - mu0 = u - 1 >= 0).
        consumed_mws = [s.get("mw") for s in consumed
                        if species_role(s) != POOL]
        if any(mw is None for mw in consumed_mws):
            return _deny(key, "classifier-not-admissible")
        defect = float(
            getattr(dst_poly, "chain_mass_defect_g_mol", 0.0) or 0.0)
        u_raw = sum(float(mw) for mw in consumed_mws) / monomer_mw
        a = float(gas_mw) / monomer_mw
        u = u_raw - a + defect / monomer_mw
        if not math.isfinite(u) or u < 1.0:
            return _deny(key, "landing-cone-violation")

        # G6 -- exportability [P1-5]: the chem.yaml export is load-bearing;
        # only balanced rows with plain-Arrhenius (directly serializable,
        # rmgpy/cantera.py) kinetics are admissible, which makes the
        # emission-time cantera:null outcome unreachable except by bugs
        # (which the serializer tripwire then catches by RAISING).
        if not forward.is_balanced():
            return _deny(key, "not-balanced")
        if forward.kinetics is None:
            # G6 PROVISIONAL deny (re-adjudication defect fix): the only
            # call chain into this evaluator (stamp_gas_association_refusal
            # <- make_new_reaction) runs BEFORE kinetics assignment, so a
            # family-generated row's kinetics is ALWAYS None here --
            # "kinetics-not-exportable" would misdescribe the row's eventual
            # kinetics. Deny provisionally (would_admit stays 0, fail-closed)
            # and let readjudicate_conduit_admission re-run this gate once
            # the final kinetics exists (PROVISIONAL_DENY_REASONS).
            return _deny(key, "kinetics-not-yet-assigned")
        from rmgpy.kinetics.arrhenius import Arrhenius as _Arrhenius
        if type(forward.kinetics) is not _Arrhenius:
            # FINAL, fail-closed last word: genuinely non-Arrhenius kinetics
            # (Chebyshev/MultiArrhenius/PDep/...; strict type -- subclasses
            # included) is not directly serializable and never admits.
            return _deny(key, "kinetics-not-exportable")

        # All gates passed affirmatively -> ADMIT (dead behind the flag).
        return AdmissionVerdict(
            admitted=True,
            deny_reason=None,
            chain_units=float(u),
            gas_product=(str(gas.get("token") or ""), float(gas_mw)),
            gas_units=float(a),
            dst_pool=dst_label,
            candidate_key=key,
            needs_irreversible_rewrite=bool(record["reversible"]),
        )
    except Exception as exc:  # G7 -- anything else: DENY, loudly.
        logging.warning(
            "MOMENT-CREDIT CONDUIT ADMISSION (M18.3): evaluation failed for "
            "a candidate row (%s: %s); DENYING fail-closed "
            "(admission-evaluation-error). Admission code never recovers "
            "into admit.", type(exc).__name__, exc)
        return _deny(key, f"admission-evaluation-error:{type(exc).__name__}")


def admission_census_suffix(verdict):
    """APPEND-ONLY census tokens carrying the would-be admission verdict
    while :data:`CONDUIT_ADMISSION_ENABLED` is False (BUILD_SPEC W1.6):
    an ADMISSIBLE verdict appends ``would_admit=1 deny=None
    rewrite=<bool>``; a denied verdict appends ``deny=<reason>``. The
    tokens ride their own bracketed group AFTER the M18.2
    ``[conduit-census/1 ...]`` suffix, so the existing header/tokens (and
    the line's closing bracket) stay byte-identical (round-36 P1(b))."""
    if verdict is None:
        return ""
    if verdict.admitted:
        return (f" [conduit-admission/1 would_admit=1 deny=None "
                f"rewrite={verdict.needs_irreversible_rewrite}]")
    return (f" [conduit-admission/1 would_admit=0 "
            f"deny={verdict.deny_reason}]")


# ---------------------------------------------------------------------------
# M18.3: run-level conduit flux accumulator API (DESIGN §4.4)
# ---------------------------------------------------------------------------

#: Run-level admitted-conduit gas-mass accumulator,
#: {candidate_key: {"grams": float, "revoked": bool}}. ALWAYS EMPTY in
#: M18.3: the writer that populates it is the solver's conduit dispatch arm
#: (M18.4); the artifact-side census writer + both loaders land now so the
#: block's contract is pinned before any real number exists.
_CONDUIT_FLUX_TOTALS = {}


def get_conduit_flux_totals():
    """Snapshot of the run-level conduit gas-mass accumulator:
    ``{candidate_key: {"grams": <float g>, "revoked": <bool>}}``.
    Cumulative over the generating run, additive across rebuild epochs; a
    REVOKED row's accumulated mass stays counted (it happened). Empty until
    the M18.4 solver dispatch populates it."""
    return {k: dict(v) for k, v in _CONDUIT_FLUX_TOTALS.items()}


def annotate_feature_radical(reaction_label):
    """String-only annotation hook for the FEATURE-RADICAL census line
    (emitted via :func:`rmgpy.polymer._warn_once_refused`; its caller lives
    in the solver, which M18.2 must not edit, so this hook works from the
    census label alone -- sufficient because every feature-radical row
    buckets FEATURE_RADICAL regardless of shape). NEVER raises."""
    try:
        key = candidate_key_from_label(reaction_label)
        entry = register_candidate(key, "feature_radical", "FEATURE_RADICAL")
        result = {"candidate_key": key, "shape": None, "admissible": False}
        # r42 P1-4(a): EVERY census line carries the admission verdict
        # tokens (BUILD_SPEC W1.6 promised them on every line, and the FR
        # line lacked them). An FR sighting is the upstream blocker: once
        # registered (the line above), any admission evaluation of this
        # key G1-denies feature-radical-overlap -- printed here so the FR
        # census line and the ledger can never disagree.
        return (census_suffix(result, entry)
                + admission_census_suffix(
                    _deny(key, "feature-radical-overlap")))
    except Exception as exc:  # pragma: no cover - defensive fail-open
        logging.warning(
            "MOMENT-CREDIT CONDUIT CENSUS (M18.2): feature-radical "
            "annotation failed (%s: %s); the refusal itself is unaffected "
            "(census-only code path).", type(exc).__name__, exc)
        return ""
