"""Census-only refusal classifier for the moment-credit conduit
(item-20, increment M18.2; adjudicated rounds 32/36).

Classifies refused reaction rows (r93 general-branch refusals and, for
the narrow slice of upstream FEATURE-RADICAL refusals genuinely gated on
QSSA validity -- see :data:`GENUINE_FEATURE_RADICAL_REASONS`, currently
'qssa-invalid' and 'qssa-unassessable' -- FEATURE-RADICAL refusals too)
into admission buckets for the future surrogate archetype
``moment_credit_conduit/1``. First increment (round-32
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
          "census": "r93_general" | "feature_radical" | "conduit_echo",
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
    CONDUIT_ECHO                  a refusal echoing through the FR warn-once
                                   hook for a NON-genuine reason (round-50
                                   FR-census scoping; e.g. 'conduit-deferred'):
                                   NOT an upstream feature-radical blocker,
                                   just the same candidate refusing again
                                   through the shared warn-once site. G1
                                   (:func:`evaluate_conduit_admission`)
                                   ignores it entirely -- it must never
                                   deny-by-overlap its own (or any co-keyed)
                                   candidate the way real 'feature_radical'
                                   membership does.
    UNCLASSIFIED                  shape not in the A-F vocabulary (reason set)

Precedence when one candidate key is seen across multiple censuses
(overlap ledger, round-36 P1(a) / round-50 scoping): 'feature_radical'
(the genuine upstream blocker) ranks highest and wins over everything;
'conduit_echo' ranks LOWEST -- any other census on the key (including
r93_general) beats it for ``effective_bucket``, since an echo sighting
carries no classification information of its own about the row.

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

(a) STABLE CANDIDATE KEY + OVERLAP PRECEDENCE. A non-trivial share of rows
    appear in BOTH the r93 and feature-radical censuses (the exact overlap
    count is run-dependent and tracked live by :func:`census_summary`'s
    ``overlap=`` token, not pinned here). :func:`conduit_candidate_key`
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

    SAME-CENSUS RESIGHT RESOLUTION (round-50 P1, determinism finding):
    when the SAME (candidate, census) pair is sighted more than once with
    DIFFERENT buckets (e.g. a row re-classified across rebuild epochs),
    that census's contribution to ``effective_bucket`` is no longer
    last-write-wins (which made the final ledger state depend on sighting
    order). Instead the registry keeps the FULL SET of buckets ever
    sighted for that (candidate, census) pair and resolves it
    deterministically via :data:`BUCKET_DECLARATION_ORDER`: the bucket
    appearing LATEST in that fixed vocabulary order wins ("most
    conservative wins" -- ADMISSIBLE_* loses to DEFERRED_*, which loses to
    FEATURE_RADICAL/CONDUIT_ECHO, which loses to UNCLASSIFIED; see
    :func:`_most_conservative_bucket`). This makes the resolved
    per-census bucket, and therefore ``effective_bucket``, a pure function
    of the SET of sightings -- order- and epoch-boundary-independent, as
    this module's determinism contract requires. A (candidate, census)
    pair that ever accumulates more than one distinct bucket is a
    "resight divergence"; the count of such pairs is surfaced, loudly, as
    an appended ``resight_divergence=<n>`` token on :func:`census_summary`
    (never inserted into the per-row :func:`census_suffix`, which stays
    unchanged). ``unclassified_total`` is likewise computed fresh from
    the current registry state on every access (never incrementally
    accumulated during ingestion), so it can never desync from the
    sightings that produced it.
"""

import atexit
import base64
import hashlib
import logging
import math
import re
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

# ---------------------------------------------------------------------------
# round-53 adjudication: DECLARED vs ACTIVE pool-state label oracle
# ---------------------------------------------------------------------------
#
# :data:`POOL_STATE_RESOLVABLE_LABELS` above is the DECLARED set: a
# deck-curated, human-reviewable snapshot (historically drifted -- an old
# "conduit3" deck's trimer labels turned out, on a real conduit4 run, to be
# 108/108 core discrete trimer RADICALS, never pool species: the label
# oracle and the in-repo isomorphism oracle disagreed on every single
# sighting). The census layer must never trust DECLARED membership as
# ground truth on its own; every declared label is validated against the
# in-repo isomorphism oracle (:func:`rmgpy.polymer._discrete_resolves_to_pool_state`)
# before it is trusted at runtime, producing the ACTIVE set (declared
# members that pass validation). The oracle may only validate-or-DROP a
# DECLARED label -- it never adds a label the deck did not declare.
#
# Validation reasons (closed vocabulary):
LABEL_VALIDATION_ISO_PASS = "iso-pass"              # kept active
LABEL_VALIDATION_ISO_FAIL = "iso-fail"              # declared, not isomorphic to any pool
LABEL_VALIDATION_MISSING_SPECIES = "missing-species"  # declared label never sighted
LABEL_VALIDATION_AMBIGUOUS = "ambiguous"            # matches more than one distinct pool


def _validate_declared_label(species, row_pools):
    """Validate a DECLARED pool-state label against a row's registered
    pools using the in-repo isomorphism oracle
    (:func:`rmgpy.polymer._discrete_resolves_to_pool_state`) as the
    pass/fail gate. Returns ``(status, reason, pool_or_None)``.

    The oracle itself only reports pass/fail; the extra per-pool walk
    (``poly.is_isomorphic(species)``, the exact same test the oracle
    applies internally) recovers WHICH pool matched, for the LABEL
    VALIDATION line's ``pool=`` field, and detects a label mapping onto
    more than one distinct pool proxy (``ambiguous``)."""
    from rmgpy.polymer import _discrete_resolves_to_pool_state
    if not _discrete_resolves_to_pool_state(species, row_pools):
        return "dropped", LABEL_VALIDATION_ISO_FAIL, None
    matches = []
    for poly in row_pools or []:
        try:
            if poly.is_isomorphic(species):
                pool_label = getattr(poly, "label", None) or ""
                if pool_label not in matches:
                    matches.append(pool_label)
        except (ValueError, AttributeError):
            continue
    if len(matches) > 1:
        return "dropped", LABEL_VALIDATION_AMBIGUOUS, None
    if len(matches) == 1:
        return "active", LABEL_VALIDATION_ISO_PASS, matches[0]
    # Defensive: the oracle said True (e.g. species IS a Polymer, the
    # oracle's isinstance short-circuit) but the re-walk found no named
    # pool match. Treat as a genuine pass with no named pool rather than
    # silently dropping a label the oracle itself validated.
    return "active", LABEL_VALIDATION_ISO_PASS, None


class _LabelOracleState:
    """Process-wide, thread-safe DECLARED -> ACTIVE pool-state label
    validation ledger (round-53).

    Species/pools become known to this census layer per-refused-row (see
    :func:`record_from_reaction`'s ``row_pools`` parameter) rather than
    all up front at a single init point, so validation is LAZY: each
    declared label is validated ONCE, the first time a species carrying
    that label is sighted (:meth:`note_sighting`) -- never re-validated,
    so a later mid-run drift shows up as a runtime DIVERGENCE (see
    :func:`_species_entry`), not a silent re-validation. A declared label
    that is never sighted during the run is resolved as
    ``missing-species`` when the health line is finalized
    (:meth:`finalize`) -- there is no species to validate it against, so
    it can never earn ``iso-pass``.

    LIFECYCLE (round-55 P1-2/P1-3 redesign): finalization is an explicit
    END-OF-LIFECYCLE event, never a read-path side effect.
    :func:`census_summary` is READ-ONLY; the boundary emitters are the
    epoch reset (:meth:`reset` -- production-called via
    :func:`reset_conduit_state` from rmgpy/rmg/main.py RMG.initialize)
    and the module's process-exit ``atexit`` hook, each of which closes
    the CURRENT lifecycle (finalize + one-shot health line) before any
    state is cleared. Exactly one health line exists per lifecycle -- a
    lifecycle opened by a reset or by any declared-label sighting --
    even when zero census rows fired in it; the virgin import->first-reset
    interval is not a lifecycle and emits nothing. A sighting landing
    AFTER true finalization keeps the finalized verdict and, when that
    verdict is ``dropped``, surfaces a loud versioned
    ``CONDUIT CLASSIFIER ORACLE ANOMALY/1`` line (once per label per
    lifecycle) instead of raising: this module's census-only contract is
    that its code paths NEVER raise into generation code, and the
    defensive fail-open handlers upstream (annotate_refused_row / G7)
    would swallow a raise into a generic annotation failure -- strictly
    quieter than the anomaly line.
    """

    def __init__(self):
        self._lock = threading.Lock()
        self._validated = {}   # label -> (status, reason, pool)
        self._finalized = False
        self._health_line = None
        self._open = False               # a lifecycle is in progress
        self._post_final_warned = set()  # anomaly warn-once, per lifecycle
        # round-60 P2-2: warn-once dedup for the
        # concurrent-validation-disagreement anomaly (see note_sighting),
        # keyed by label, per lifecycle. Cleared in reset() alongside the
        # other lifecycle-scoped warn-once state below.
        self._disagreement_warned = set()

    def reset(self):
        # round-55 P1-2/P1-3: the reset boundary CLOSES the outgoing
        # lifecycle -- finalize (resolving never-sighted declared labels
        # as missing-species) and emit its one-shot health line BEFORE
        # clearing -- then opens a fresh one.
        self.close_lifecycle()
        with self._lock:
            self._validated = {}
            self._finalized = False
            self._health_line = None
            self._post_final_warned = set()
            self._disagreement_warned = set()
            self._open = True

    def close_lifecycle(self):
        """Boundary emitter (round-55 P1-2; round-56 F1): finalize-and-emit
        the one-shot health line for the CURRENT lifecycle if one is open
        (opened by a reset or by any declared-label sighting); a no-op on
        virgin process state, so the process's first reset does not
        report a lifecycle that never existed.

        round-56 F1(b): idempotent per lifecycle -- a second close on an
        already-closed (finalized) lifecycle returns the cached health line
        WITHOUT re-emitting anything. This lets BOTH the initialize-time
        reset boundary (:meth:`reset`) and the new production end-of-run
        hook (:func:`close_conduit_lifecycle`, wired at rmgpy/rmg/main.py
        RMG.finish) call it safely without ever double-emitting. finalize()
        is itself idempotent via ``_finalized``; this explicit short-circuit
        makes the no-double-emit contract local to close_lifecycle and
        independent of finalize internals."""
        with self._lock:
            if self._finalized:
                # Already closed this lifecycle -- no-op, cached line only.
                return self._health_line
            active = self._open or bool(self._validated)
        if not active:
            return None
        return self.finalize()

    def note_sighting(self, label, species, row_pools):
        """Validate a DECLARED label the first time it is sighted attached
        to a live species; a no-op for already-validated or non-declared
        labels. NEVER raises (census-only code path).

        round-55 P1-3: a sighting arriving AFTER true finalization keeps
        the finalized verdict; if that verdict is ``dropped`` (chiefly
        ``missing-species`` -- a late sighting directly contradicts it)
        the contradiction is surfaced as a loud versioned ANOMALY/1 line,
        once per label per lifecycle. See the class docstring for why the
        anomaly line is chosen over raising.

        round-56 F2 (finalize/note_sighting race, approach ii): the label
        validation (:func:`_validate_declared_label`) runs OUTSIDE the lock
        -- it is a cheap in-repo isomorphism check, but re-acquiring the
        lock inside it would deadlock the deterministic interleave test and,
        in production, a concurrent :meth:`finalize` can still slip into the
        gap between the cached check and the setdefault. Approach (ii) is
        used instead of holding the lock across validation: after the
        setdefault we RE-CHECK the finalized state under the same lock, and
        if finalize() won the race and marked this label dropped/
        missing-species in the gap, we route through the SAME loud
        post-finalization ANOMALY/1 path a normally-ordered late sighting
        takes -- never silently returning the dropped verdict with no
        anomaly. (Approach ii chosen over "hold the lock across validation"
        precisely because the validation callable can re-enter the oracle:
        the interleave test drives finalize() from inside it, which a held
        non-reentrant lock would deadlock.)"""
        if label not in POOL_STATE_RESOLVABLE_LABELS:
            return None
        with self._lock:
            self._open = True
            cached = self._validated.get(label)
            warn_anomaly = (
                self._finalized and cached is not None
                and cached[0] == "dropped"
                and label not in self._post_final_warned)
            if warn_anomaly:
                self._post_final_warned.add(label)
        if warn_anomaly:
            logging.warning(
                "CONDUIT CLASSIFIER ORACLE ANOMALY/1 label=%s "
                "event=post-finalization-sighting status=dropped reason=%s "
                "action=keep-finalized-verdict", label, cached[1])
        if cached is not None:
            return cached
        try:
            computed = _validate_declared_label(species, row_pools)
        except Exception:  # pragma: no cover - defensive fail-closed
            computed = ("dropped", LABEL_VALIDATION_ISO_FAIL, None)
        with self._lock:
            result = self._validated.setdefault(label, computed)
            won = result is computed
            # F2 race close: finalize() interleaved and populated this label
            # as dropped while we validated outside the lock. Do not return
            # that dropped verdict silently -- flag the post-finalization
            # anomaly here, under the lock, exactly once per label.
            race_anomaly = (
                not won and self._finalized and result[0] == "dropped"
                and label not in self._post_final_warned)
            if race_anomaly:
                self._post_final_warned.add(label)
            # round-60 P2-2: dedupe the concurrent-validation-disagreement
            # anomaly below to once per (label, lifecycle) -- the check and
            # the add happen together under self._lock (same lock/reset
            # discipline as every other warn-once set in this class) so a
            # losing race on a hot mistyped/flaky label cannot spam this
            # line once per losing thread.
            disagreement_anomaly = (
                not won and not race_anomaly and result != computed
                and label not in self._disagreement_warned)
            if disagreement_anomaly:
                self._disagreement_warned.add(label)
        if won:
            # This call won the (label not yet validated) race -- log it.
            status, reason, pool = result
            logging.warning(
                "CONDUIT CLASSIFIER LABEL VALIDATION/1 label=%s status=%s "
                "reason=%s pool=%s", label, status, reason, pool or "none")
        elif race_anomaly:
            logging.warning(
                "CONDUIT CLASSIFIER ORACLE ANOMALY/1 label=%s "
                "event=post-finalization-sighting status=dropped reason=%s "
                "action=keep-finalized-verdict", label, result[1])
        elif disagreement_anomaly:
            # P1-B (round-58): two concurrent FIRST sightings of the SAME
            # label both validate outside the lock; setdefault() picks a
            # winner deterministically (first-write-wins), but until now a
            # disagreement between the winner's cached verdict and this
            # (losing) call's freshly-computed verdict vanished silently --
            # the loser's result was simply discarded. Keep the cached
            # (first-stored) verdict -- first-write-wins is still the
            # resolution rule -- but make the disagreement loud instead of
            # silent. round-60 P2-2: rate-limited to once per (label,
            # lifecycle) -- see disagreement_anomaly above -- since every
            # losing thread in a race would otherwise re-log this line.
            logging.warning(
                "CONDUIT CLASSIFIER ORACLE ANOMALY/1 label=%s "
                "reason=concurrent-validation-disagreement kept=%s "
                "discarded=%s", label, result, computed)
        return result

    def active_labels(self):
        with self._lock:
            return frozenset(l for l, (status, _reason, _pool)
                             in self._validated.items() if status == "active")

    def finalize(self):
        """Resolve any never-sighted DECLARED label as ``missing-species``
        and emit the one-shot health line (idempotent: repeat calls --
        e.g. the reset boundary and the atexit hook racing at process
        end, or an explicit :func:`finalize_label_oracle` before either --
        return the already-computed line without re-emitting anything).

        round-55 P1-4: the health line says ``validation=lazy-first-sight``
        -- validation IS lazy (first-sight per label, finalize-at-boundary
        for never-sighted labels); the token it replaces
        (``validation=init``) claimed an init-time validation that never
        existed. The rename is legal under the append-only token rule
        because this line class was born in the same unpushed commit."""
        with self._lock:
            if self._finalized:
                return self._health_line
            newly_missing = []
            for label in POOL_STATE_RESOLVABLE_LABELS:
                if label not in self._validated:
                    self._validated[label] = (
                        "dropped", LABEL_VALIDATION_MISSING_SPECIES, None)
                    newly_missing.append(label)
            configured = len(POOL_STATE_RESOLVABLE_LABELS)
            active = sum(1 for status, _reason, _pool
                        in self._validated.values() if status == "active")
            dropped = configured - active
            line = (f"CONDUIT CLASSIFIER ORACLE HEALTH/1 configured={configured} "
                    f"active={active} dropped={dropped} source=deck "
                    f"validation=lazy-first-sight")
            self._health_line = line
            self._finalized = True
        for label in newly_missing:
            logging.warning(
                "CONDUIT CLASSIFIER LABEL VALIDATION/1 label=%s status=dropped "
                "reason=missing-species pool=none", label)
        logging.warning(line)
        return line


#: Process-wide DECLARED -> ACTIVE label oracle instance.
_LABEL_ORACLE = _LabelOracleState()


def finalize_label_oracle():
    """Explicit end-of-lifecycle finalization for the DECLARED/ACTIVE
    label oracle (round-55 P1-3 redesign): resolve every never-sighted
    declared label as ``missing-species`` and emit the one-shot
    ``CONDUIT CLASSIFIER ORACLE HEALTH/1`` line; returns that line.
    Idempotent within a lifecycle. Production never needs to call this
    directly -- the epoch-reset boundary (:func:`reset_conduit_state`,
    called from rmgpy/rmg/main.py RMG.initialize), the production
    end-of-run hook (:func:`close_conduit_lifecycle`, wired at RMG.finish),
    and the process-exit ``atexit`` guard below all fire it automatically,
    which is what keeps :func:`census_summary` READ-ONLY."""
    return _LABEL_ORACLE.finalize()


def close_conduit_lifecycle():
    """PRODUCTION end-of-run close for the label-oracle lifecycle (round-56
    F1a): finalize the CURRENT lifecycle and emit its one-shot
    ``CONDUIT CLASSIFIER ORACLE HEALTH/1`` line if one is open, else no-op.
    Wired at the deterministic end-of-generation point (rmgpy/rmg/main.py
    RMG.finish, the tail of RMG.execute) so this run's health line lands in
    THIS run's log rather than being deferred to the next run's initialize
    or left to the fragile process-exit path. Idempotent
    (:meth:`_LabelOracleState.close_lifecycle`): safe to call alongside the
    reset boundary and the atexit guard without double-emitting.

    Design §1.2 'Close-hook extension' (round-68 amendment 6): ALSO closes
    the CORE EPOCH provider's lifecycle, emitting one ``CONDUIT EPOCH
    MAP/1`` line per distinct epoch advanced this lifecycle
    (:meth:`_EpochProvider.close_lifecycle`) -- idempotent and never-raising
    exactly like the oracle close above. The return value stays the
    oracle's health line (unchanged shape/contract for existing callers);
    the epoch-provider close is invoked for its emission side effect.

    Design §2.5 (commit-3, round-72 P3 fold-in #2): ALSO closes the
    standing-admit registry's lifecycle, emitting one
    ``CONDUIT STANDING ADMIT HEALTH/1`` line. Each of the three surfaces
    (oracle, epoch provider, standing-admit registry) is now wrapped in ITS
    OWN never-raise guard: :meth:`_LabelOracleState.close_lifecycle` and
    :meth:`_EpochProvider.close_lifecycle` and
    :meth:`_StandingAdmitRegistry.close_lifecycle` are already internally
    never-raising by contract, but this function no longer trusts that
    blindly -- a defensive ``try/except`` per surface means one surface
    raising unexpectedly (e.g. a monkeypatched/broken close in a test, or a
    future regression) can never skip the other two surfaces' emission.
    Each guard logs a rate-limited anomaly line and continues; NEVER
    raises overall."""
    try:
        oracle_line = _LABEL_ORACLE.close_lifecycle()
    except Exception as exc:  # pragma: no cover - defensive fail-closed
        _raw_name = type(exc).__name__
        _sanitized_name = _sanitize_dynamic_token(_raw_name)
        logging.warning(
            "CONDUIT LIFECYCLE CLOSE ANOMALY/1 surface=oracle reason=%s%s "
            "action=continue-remaining-surfaces",
            _sanitized_name, _raw_sha_suffix(_raw_name, _sanitized_name))
        oracle_line = None
    try:
        _EPOCH_PROVIDER.close_lifecycle()
    except Exception as exc:  # pragma: no cover - defensive fail-closed
        _raw_name = type(exc).__name__
        _sanitized_name = _sanitize_dynamic_token(_raw_name)
        logging.warning(
            "CONDUIT LIFECYCLE CLOSE ANOMALY/1 surface=epoch-provider "
            "reason=%s%s action=continue-remaining-surfaces",
            _sanitized_name, _raw_sha_suffix(_raw_name, _sanitized_name))
    try:
        _STANDING_ADMITS.close_lifecycle()
    except Exception as exc:  # pragma: no cover - defensive fail-closed
        _raw_name = type(exc).__name__
        _sanitized_name = _sanitize_dynamic_token(_raw_name)
        logging.warning(
            "CONDUIT LIFECYCLE CLOSE ANOMALY/1 surface=standing-admits "
            "reason=%s%s action=continue-remaining-surfaces",
            _sanitized_name, _raw_sha_suffix(_raw_name, _sanitized_name))
    return oracle_line


# round-56 F1(c)(d): atexit is KEPT, but strictly as a LAST-RESORT guard for
# processes that crash or exit early WITHOUT reaching the RMG.finish hook
# (:func:`close_conduit_lifecycle`). For a normal run the finish hook is the
# primary close point and it runs first; because close_lifecycle() is
# idempotent (no-op once finalized), the atexit guard then emits nothing --
# so the health line is deterministic and never doubled.
#
# Double-registration hazard (importlib.reload re-executes this module body
# in the SAME module namespace, which would register a second callback):
# guard the register() call with a flag stored on the persistent ``atexit``
# module itself (NOT a module-level name, which reload would reset). The
# registered function references ``close_conduit_lifecycle`` BY NAME, so a
# reload that rebinds the module global still closes the CURRENT singletons
# through the same (in-place-updated) module __dict__.
#
# round-70 P2-d: this last-resort path used to call ONLY
# ``_LABEL_ORACLE.close_lifecycle()``, skipping the CORE EPOCH provider
# entirely -- a process that crashed/exited before RMG.finish reached
# _conduit_lifecycle_close would lose that run's EPOCH MAP/HEALTH lines even
# though the label-oracle health line still landed. Route through the SAME
# :func:`close_conduit_lifecycle` the production finish hook uses so no
# surface is skipped; it is idempotent and never-raising exactly like the
# oracle-only call it replaces.
def _close_conduit_lifecycle_atexit():  # pragma: no cover - process teardown
    close_conduit_lifecycle()


if not getattr(atexit, "_conduit_lifecycle_atexit_registered", False):
    atexit.register(_close_conduit_lifecycle_atexit)
    atexit._conduit_lifecycle_atexit_registered = True


def _canonical_sig_sha16(signature):
    """Token-safe digest of a raw core-topology signature (design §1.1
    'Token-safe digest encoding'): the raw signature is a nested
    ``((label,index)...),((index,label)...))`` tuple (`polymer.py:7703-7707`)
    -- not itself log-token-safe -- so only its sha256 (truncated to the
    first 16 hex chars) ever leaves this module. ``repr()`` of a tuple of
    hashable, order-stable primitives is deterministic across processes for
    the same signature, matching the module's determinism contract (design
    §1.2 'Determinism'). NEVER raises: an unrepresentable/exotic signature
    still yields a stable, charset-safe fallback digest rather than
    propagating an exception into the census/generation path."""
    try:
        payload = repr(signature).encode("utf-8", errors="backslashreplace")
    except Exception:  # pragma: no cover - defensive fail-closed
        payload = repr(None).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()[:16]


#: round-70 P2-e: warn-once dedup state for the advance-after-close anomaly
#: line in :meth:`_EpochProvider.advance` -- a single lifecycle key (there
#: is only one meaningful event here, unlike the ``(census, bucket)``-keyed
#: sets below) so a hot post-close caller cannot spam the log once per call.
#: Checked/set under :data:`_WARN_ONCE_LOCK` via :func:`_warn_once` (both
#: defined later in this module; referenced here only at call time, after
#: the whole module has finished loading). Cleared once per lifecycle by
#: :func:`reset_conduit_state`, mirroring :data:`_bucket_anomaly_warned`.
_epoch_advance_after_close_warned = set()

#: round-71 P1: warn-once dedup state for the sighting-on-burned-epoch
#: anomaly line in :meth:`_EpochProvider.note_sighted` -- keyed by the
#: (already-burned) ordinal token, since more than one such token could in
#: principle exist per lifecycle. Should never fire in production (a burn
#: only follows a raise, which aborts the rebuild before any further
#: sighting could register), but stays rate-limited defensively like every
#: other anomaly line in this class. Checked/set under
#: :data:`_WARN_ONCE_LOCK` via :func:`_warn_once`. Cleared once per
#: lifecycle by :func:`reset_conduit_state`, mirroring
#: :data:`_epoch_advance_after_close_warned`.
_epoch_sighting_on_burned_warned = set()


class _EpochProvider:
    """Process-wide CORE EPOCH provider singleton (design §1.1-§1.4),
    mirroring the :class:`_CensusRegistry` (`CENSUS_REGISTRY`) and
    :class:`_LabelOracleState` (`_LABEL_ORACLE`) singleton patterns above:
    one module-global instance, one :class:`threading.Lock` guarding all
    state, an idempotent-per-lifecycle close hook, and a lifecycle reset
    that the module's "reset both or neither" boundary
    (:func:`reset_conduit_state`) drives alongside the other singletons.

    A CORE EPOCH is a run-local ATTEMPTED-rebuild ordinal bound to a
    core-topology signature (design §1.1): ``current_epoch()`` returns a
    short ``e{N}`` token (never the raw signature -- see
    :func:`_canonical_sig_sha16`), or the pre-first-advance sentinel
    ``"epre"``. Epochs are monotone in TIME, not in signature-space: a
    signature that reverts to an earlier value (A->B->A) still earns a NEW
    ordinal (design §1.2).

    Every public method NEVER raises into the census/generation path
    (module contract, design §1.4 'Constraints fixed'): on an internal
    error each method logs a versioned ``CONDUIT EPOCH PROVIDER ANOMALY/1``
    line (mirroring the ``CONDUIT CLASSIFIER ORACLE ANOMALY/1`` pattern
    above) and degrades to the last-known-good state rather than
    propagating the exception."""

    def __init__(self):
        self._lock = threading.Lock()
        self._ordinal = -1                    # -1 before first advance
        self._current_sig = None              # last-advanced raw signature
        self._sig_by_ordinal = {}             # {ordinal: sha256_hex[:16]}
        self._advanced_this_lifecycle = []     # [(ordinal, sha16), ...]
        # round-70 P1-a: burned accounting is a SET of burned ordinal
        # TOKENS (e.g. "e3") rather than a blind incrementing counter, so
        # burning the same ordinal twice (idempotent double-report) is a
        # no-op instead of double-counting. A rebuild failure that follows
        # a DEDUP no-op advance (no new ordinal was created) can never
        # burn anything -- it increments the separate ``_failed_attempts``
        # counter instead (see :meth:`note_conduit_rebuild_failed`).
        self._burned_ordinals = set()
        self._failed_attempts = 0
        # round-71 P1: SET of ordinal TOKENS a census producer has
        # resolved epoch=None (or an explicit epoch=) through THIS
        # provider for (:meth:`note_sighted`, driven from
        # :meth:`_CensusRegistry.register`). A token in this set was
        # ACTUALLY SEEN carrying a real census sighting -- the burn
        # decision in :meth:`note_conduit_rebuild_failed` consults it so
        # "burned" never again means "created AND rebuild raised" when the
        # ordinal in fact went on to be sighted first (the finding: the
        # polymer initialize path registers sightings early, then a LATER
        # tripwire in the same initialize can still raise).
        self._sighted_ordinals = set()
        # round-71 P1: rebuild failures whose token WAS sighted before the
        # failure hook ran -- kept separate from ``_burned_ordinals`` so
        # "burned" stays honestly "attempted, never sighted".
        # round-72 P3 (commit-3 fold-in): a SET of tokens, not a bare
        # counter -- a duplicate failure report for the SAME sighted
        # ordinal (e.g. two independent tripwires in the same rebuild both
        # calling :meth:`note_conduit_rebuild_failed` with the same
        # ``token``) must not double-count. Count is ``len(...)`` via the
        # :attr:`_failed_after_sighting` property below.
        self._failed_after_sighting_tokens = set()
        self._closed = False                   # EPOCH MAP emitted this lifecycle
        self._epoch_map_lines = None           # cached lines once closed
        # round-71 P2: set on FIRST ACTIVITY -- the first advance()
        # attempt (any outcome: created, dedup no-op, or post-close
        # degrade), the first note_sighted() resolution, or the first
        # note_conduit_rebuild_failed() report. A provider that was NEVER
        # opened has nothing to report: :meth:`close_lifecycle` emits
        # NOTHING (no MAP lines -- there are none -- and no HEALTH line
        # either) instead of a fake ``epochs=0 ... last=epre`` line, e.g.
        # for a run that fails before RMG.initialize ever opens a
        # lifecycle.
        self._opened = False

    @property
    def _failed_after_sighting(self):
        """round-72 P3: token-idempotent count -- ``len()`` of the token
        set, so a duplicate failure report against the same already-sighted
        ordinal is not double-counted. Read-only; attribute-access syntax
        is unchanged from the old bare-counter so every existing caller/
        test keeps working."""
        return len(self._failed_after_sighting_tokens)

    def current_epoch(self):
        """``f"e{ordinal}"`` if an advance has fired, else the sentinel
        ``"epre"`` (design §1.1/§1.4). Never raises; never returns
        ``None``."""
        try:
            with self._lock:
                ordinal = self._ordinal
        except Exception:  # pragma: no cover - defensive fail-closed
            return "epre"
        return f"e{ordinal}" if ordinal >= 0 else "epre"

    def advance(self, signature):
        """Advance to a new core epoch for ``signature`` (design §1.2).
        Identical-consecutive-signature is a no-op (returns the current
        token); any other signature -- including a revert to a prior value
        -- increments the ordinal and returns the new token. Records the
        signature's sha16 digest in the audit map and appends to this
        lifecycle's EPOCH MAP backlog. NEVER raises: on internal failure,
        degrades to the last-known-good token via :meth:`current_epoch`.

        Returns ``(token, created)`` (round-70 P1-a/P1-b): ``created`` is
        True only when a NEW ordinal was actually minted this call, so
        callers (chiefly ``RMG._advance_conduit_epoch``) can tell a genuine
        new epoch apart from a dedup no-op or a post-close degrade -- the
        distinction :meth:`note_conduit_rebuild_failed` needs to decide
        whether a subsequent rebuild failure may burn an ordinal.

        round-70 P2-e (advance-after-close): once :meth:`close_lifecycle`
        has fired for the current lifecycle, a further ``advance()`` must
        NOT mint a new ordinal -- those epochs would never be emitted by
        an EPOCH MAP line (close already ran). It instead returns the
        CURRENT token unchanged with ``created=False`` and logs a
        rate-limited ``CONDUIT EPOCH PROVIDER ANOMALY/1
        event=advance-after-close`` line (once per lifecycle)."""
        try:
            with self._lock:
                self._opened = True  # round-71 P2: first advance attempt
                if self._closed:
                    ordinal = self._ordinal
                    token = f"e{ordinal}" if ordinal >= 0 else "epre"
                    after_close = True
                else:
                    after_close = False
                    if signature == self._current_sig:
                        ordinal = self._ordinal
                        token = f"e{ordinal}" if ordinal >= 0 else "epre"
                        return (token, False)
                    sig_sha = _canonical_sig_sha16(signature)
                    self._ordinal += 1
                    self._current_sig = signature
                    self._sig_by_ordinal[self._ordinal] = sig_sha
                    self._advanced_this_lifecycle.append(
                        (self._ordinal, sig_sha))
                    return (f"e{self._ordinal}", True)
            # after_close is always True to reach here (every branch above
            # either returns or falls through only on after_close).
            _warn_once(
                _epoch_advance_after_close_warned, "advance-after-close",
                "CONDUIT EPOCH PROVIDER ANOMALY/1 event=advance-after-close "
                "reason=advance-after-close action=return-current-token")
            return (token, False)
        except Exception as exc:  # pragma: no cover - defensive fail-closed
            _raw_name = type(exc).__name__
            _sanitized_name = _sanitize_dynamic_token(_raw_name)
            logging.warning(
                "CONDUIT EPOCH PROVIDER ANOMALY/1 event=advance-failed "
                "reason=%s%s action=degrade-to-current-epoch",
                _sanitized_name, _raw_sha_suffix(_raw_name, _sanitized_name))
            return (self.current_epoch(), False)

    def note_sighted(self, token):
        """Mark ordinal ``token`` (e.g. ``"e3"``) SIGHTED (round-71 P1
        finding): called from :meth:`_CensusRegistry.register` every time
        a census producer (``register_candidate`` /
        ``annotate_refused_row`` / ``annotate_feature_radical``) resolves
        an ``epoch`` for a sighting through this provider -- both the
        ``epoch=None`` fallback (:func:`current_epoch`) and an
        explicitly-passed ``epoch=`` token are marked. This is cheap (one
        set add under the lock) and counts as lifecycle "activity" (round-71
        P2, see :meth:`close_lifecycle`). The pre-advance ``"epre"``
        sentinel and a missing token are never marked -- they name no real
        ordinal.

        Decision (round-71 P1, explicit-epoch case): marking is
        unconditional -- it does not first check that ``token`` is a
        currently-known/created ordinal. The ONLY consumer of the sighted
        set is the burn-vs-``failed_after_sighting`` decision below, which
        itself only ever tests an actually provider-created token (the one
        ``advance()`` returned with ``created=True``); marking a synthetic
        or test-only epoch string (e.g. ``"sig-1"``) sighted has no
        observable downstream effect, so the cheaper unconditional add was
        chosen over an extra "is this a real ordinal" membership check.

        Ordering subtlety (round-71 P1): the burn decision in
        :meth:`note_conduit_rebuild_failed` is taken at failure-hook time
        against whatever this set contains AT THAT MOMENT -- a sighting
        registered here BEFORE the failure hook runs (the exact scenario
        the finding describes: the polymer initialize path registers a
        sighting early via ``stamp_gas_association_refusal``, then a LATER
        tripwire in that same initialize raises) correctly flips the
        outcome to ``failed_after_sighting`` instead of burning. A sighting
        landing here AFTER its token has already been burned should never
        happen in production (a burn only follows a raise, which aborts
        the rebuild before any further sighting can be registered against
        it); if it somehow does, this does NOT un-burn the token or
        mutate burn history -- it logs a rate-limited
        ``CONDUIT EPOCH PROVIDER ANOMALY/1 event=sighting-on-burned-epoch``
        line, once per token per lifecycle, via the same warn-once
        machinery as the advance-after-close anomaly. NEVER raises."""
        try:
            if not token or token == "epre":
                return
            with self._lock:
                self._opened = True  # round-71 P2: sighting is activity
                already_burned = token in self._burned_ordinals
                self._sighted_ordinals.add(token)
            if already_burned:
                _warn_once(
                    _epoch_sighting_on_burned_warned, token,
                    "CONDUIT EPOCH PROVIDER ANOMALY/1 "
                    "event=sighting-on-burned-epoch epoch=%s "
                    "reason=sighting-on-burned-epoch "
                    "action=keep-burned-history", token)
        except Exception:  # pragma: no cover - defensive fail-closed
            pass

    def note_conduit_rebuild_failed(self, token=None, created=False):
        """Mark the just-advanced ordinal as BURNED (design §1.2/§1.4
        amendment 7; round-70 P1-a honest accounting; round-71 P1
        sighting-aware refinement): the advance fired AND created a new
        ordinal (``created=True``, ``token`` naming it), but the rebuild it
        labeled then raised. If that ordinal was NEVER sighted (no census
        producer resolved a sighting against it before this call), it is
        burned -- "attempted, never sighted" stays true. The burned set is
        idempotent -- burning the same ``token`` twice is a no-op, never a
        double count.

        round-71 P1: when ``created`` is True AND ``token`` WAS already
        sighted (see :meth:`note_sighted`) -- the polymer initialize path
        registers a census sighting early, then a LATER tripwire in that
        same initialize still raises -- the ordinal is NOT burned (burning
        it would corrupt "burned = attempted, never sighted" by
        simultaneously claiming the epoch was both sighted and never
        sighted). Instead the separate ``_failed_after_sighting`` counter
        increments, an honest third outcome distinct from both "burned"
        and "failed_attempts". No re-evaluation happens later: the
        sighted/not-sighted decision is taken exactly once, at the moment
        this hook runs, against whatever :attr:`_sighted_ordinals`
        contains at that instant.

        When ``created`` is False (a DEDUP no-op advance never minted an
        ordinal to burn) or ``token`` is missing (the advance itself
        failed internally -- round-70 P1-b, no reliable ordinal-created
        signal exists), nothing is burned or counted as
        failed-after-sighting; instead the separate ``_failed_attempts``
        counter increments, so a rebuild failure can never corrupt
        burned-epoch accounting by attributing itself to an ordinal that
        was never actually (re)created. NEVER raises."""
        try:
            with self._lock:
                self._opened = True  # round-71 P2: a failure is activity
                if created and token:
                    if token in self._sighted_ordinals:
                        # round-72 P3: idempotent by token -- a duplicate
                        # report for the same already-sighted ordinal is a
                        # set no-op, mirroring the burned-ordinal-set
                        # idempotency a few lines below.
                        self._failed_after_sighting_tokens.add(token)
                    else:
                        self._burned_ordinals.add(token)
                else:
                    self._failed_attempts += 1
        except Exception:  # pragma: no cover - defensive fail-closed
            pass

    def close_lifecycle(self):
        """Boundary emitter (design §1.2 'Close-hook extension', mirroring
        :meth:`_LabelOracleState.close_lifecycle`): emit one
        ``CONDUIT EPOCH MAP/1`` line per distinct epoch advanced THIS
        lifecycle, in ordinal order, followed by exactly one
        ``CONDUIT EPOCH HEALTH/1`` summary line (round-70 P1/P2-c; round-71
        P1 adds the ``failed_after_sighting=`` field) -- ``epochs=`` the
        count of ordinals actually created this lifecycle, ``burned=`` the
        size of the burned-ordinal set, ``failed_attempts=`` rebuild
        failures that could not be attributed to a created ordinal,
        ``failed_after_sighting=`` rebuild failures whose created ordinal
        WAS already sighted (so honestly NOT burned -- see
        :meth:`note_conduit_rebuild_failed`), ``last=`` the current token.
        Field order is fixed and documented here; it must not change
        without updating every format pin in the test suite. All cached
        together. Idempotent per lifecycle -- a second call before the next
        :meth:`reset` returns the cached lines WITHOUT re-emitting
        anything, so both the production end-of-run hook
        (:func:`close_conduit_lifecycle`) and the lifecycle-reset boundary
        (:meth:`reset`) can call this safely.

        round-71 P2 (virgin-lifecycle guard): if NOTHING has happened this
        lifecycle -- no advance attempt, no sighting, no recorded failure
        (:attr:`_opened` is still False) -- this emits NOTHING: no MAP
        lines (there genuinely are none) and no HEALTH line either,
        instead of the previous fake ``epochs=0 burned=0
        failed_attempts=0 failed_after_sighting=0 last=epre`` line a
        never-opened provider used to print (e.g. for a run that fails
        before RMG.initialize ever calls :func:`reset_conduit_state`).
        This is a TRUE no-op -- it does NOT set ``self._closed`` -- so a
        later real close (once the lifecycle has since seen activity)
        still fires normally; this mirrors
        :meth:`_LabelOracleState.close_lifecycle`'s existing
        open-lifecycle guard for the exact same reason.

        NEVER raises."""
        try:
            with self._lock:
                if self._closed:
                    return list(self._epoch_map_lines or [])
                if not self._opened:
                    # round-71 P2: virgin lifecycle -- stay a true no-op.
                    return []
                advanced = list(self._advanced_this_lifecycle)
                burned_count = len(self._burned_ordinals)
                failed_attempts = self._failed_attempts
                failed_after_sighting = self._failed_after_sighting
                ordinal = self._ordinal
                self._closed = True
            last = f"e{ordinal}" if ordinal >= 0 else "epre"
            lines = []
            for o, sig_sha in advanced:
                parent = f"e{o - 1}" if o > 0 else "-"
                lines.append(
                    f"CONDUIT EPOCH MAP/1 epoch=e{o} "
                    f"sig_sha={sig_sha} parent={parent}")
            lines.append(
                f"CONDUIT EPOCH HEALTH/1 epochs={len(advanced)} "
                f"burned={burned_count} failed_attempts={failed_attempts} "
                f"failed_after_sighting={failed_after_sighting} "
                f"last={last}")
            with self._lock:
                self._epoch_map_lines = lines
            for line in lines:
                logging.warning(line)
            return list(lines)
        except Exception as exc:  # pragma: no cover - defensive fail-closed
            _raw_name = type(exc).__name__
            _sanitized_name = _sanitize_dynamic_token(_raw_name)
            logging.warning(
                "CONDUIT EPOCH PROVIDER ANOMALY/1 event=close-lifecycle-failed "
                "reason=%s%s action=skip-map-emission",
                _sanitized_name, _raw_sha_suffix(_raw_name, _sanitized_name))
            return []

    def reset(self):
        """Lifecycle reset (design §1.2 'Lifecycle reset', the same "reset
        both or neither" discipline :func:`reset_conduit_state` applies to
        every other singleton): close the OUTGOING lifecycle first (so its
        ``CONDUIT EPOCH MAP/1`` lines are emitted from pre-clear state --
        same ordering the label oracle uses), then clear all state so the
        next run restarts at ``epre`` -> ``e0``. NEVER raises."""
        self.close_lifecycle()
        try:
            with self._lock:
                self._ordinal = -1
                self._current_sig = None
                self._sig_by_ordinal = {}
                self._advanced_this_lifecycle = []
                self._burned_ordinals = set()
                self._failed_attempts = 0
                self._sighted_ordinals = set()
                self._failed_after_sighting_tokens = set()
                self._closed = False
                self._epoch_map_lines = None
                self._opened = False
        except Exception:  # pragma: no cover - defensive fail-closed
            pass

    def attempted_and_burned(self):
        """Read-only snapshot of ``(epochs_seen, burned)`` -- the count of
        ATTEMPTED core-epoch advances (``ordinal + 1``, 0 if the provider
        is still virgin) and the size of the burned-ordinal set. Consumed
        by :class:`_StandingAdmitRegistry`'s own health line (design §2.5)
        so the two health lines read consistently without that registry
        reaching into this provider's private attributes directly. NEVER
        raises."""
        try:
            with self._lock:
                return (self._ordinal + 1, len(self._burned_ordinals))
        except Exception:  # pragma: no cover - defensive fail-closed
            return (0, 0)


#: Process-wide CORE EPOCH provider instance (design §1.1).
_EPOCH_PROVIDER = _EpochProvider()


def current_epoch():
    """Module-level accessor for :meth:`_EpochProvider.current_epoch`
    (design §1.2)."""
    return _EPOCH_PROVIDER.current_epoch()


def advance_conduit_epoch(signature):
    """Module-level accessor for :meth:`_EpochProvider.advance` (design
    §1.2). Production callers (rmgpy/rmg/main.py) call this before every
    RMG-owned polymer rebuild/simulate with
    ``core_topology_signature(core_species, core_reactions)``.

    round-70 P1-a/P1-b: returns ``(token, created)`` -- ``created`` is True
    only when a NEW ordinal was actually minted, so callers can tell a
    genuine advance apart from a dedup no-op (or a post-close degrade) and
    pass that outcome on to :func:`note_conduit_rebuild_failed` honestly."""
    return _EPOCH_PROVIDER.advance(signature)


def note_conduit_rebuild_failed(token=None, created=False):
    """Module-level accessor for
    :meth:`_EpochProvider.note_conduit_rebuild_failed` (design §1.2/§1.4
    amendment 7; round-70 P1-a/P1-b honest accounting). Production callers
    (rmgpy/rmg/main.py) call this from the ``except`` around each guarded
    advance site when the labeled rebuild then raises, passing the
    ``(token, created)`` outcome of the immediately-preceding
    :func:`advance_conduit_epoch` call at that same site: only a
    ``created=True`` advance's ordinal may ever be burned; everything else
    (a dedup no-op, or an advance that itself failed internally) increments
    ``failed_attempts`` instead."""
    return _EPOCH_PROVIDER.note_conduit_rebuild_failed(
        token=token, created=created)


def active_pool_state_labels():
    """The current ACTIVE (validated) subset of
    :data:`POOL_STATE_RESOLVABLE_LABELS` (round-53 DECLARED/ACTIVE split).
    Never ground truth on its own before validation -- see
    :class:`_LabelOracleState`."""
    return _LABEL_ORACLE.active_labels()


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
#: P2a fix (classify_record/conduit_echo contract): reason attached to a
#: record whose ``census`` is 'conduit_echo' -- a candidate echoing back
#: through the same warn-once hook a genuine feature_radical sighting
#: uses, for a non-genuine reason (see GENUINE_FEATURE_RADICAL_REASONS).
#: Distinct from REASON_FEATURE_RADICAL: an echo sighting is never itself
#: an upstream blocker.
REASON_CONDUIT_ECHO = "conduit-echo-non-genuine-refusal"

_SHAPE_TABLE = {
    (("CHAIN",), ("DISC", "POOL")): "A",
    (("CHAIN", "DISC"), ("DISC", "POOL")): "B",
    (("CHAIN",), ("POOL", "POOL")): "C",
    (("CHAIN", "CHAIN"), ("POOL", "POOL")): "D",
    (("CHAIN", "POOL"), ("POOL", "POOL")): "E",
    (("CHAIN", "DISC"), ("POOL",)): "F",
}

#: Overlap precedence rule (round-36 P1(a)): the upstream blocker wins.
#: Full order (highest to lowest), see also
#: :data:`CONDUIT_ECHO_CENSUS`/:meth:`_CensusRegistry._effective_bucket_below_fr`
#: (round-50 determinism finding): 'feature_radical' > every other census
#: > 'conduit_echo'.
PRECEDENCE_RULE = "feature_radical"

#: round-50 P1 same-census resight resolution: the FULL bucket vocabulary,
#: in a fixed total order from LEAST conservative (most willing to admit)
#: to MOST conservative (most willing to keep a row un-admitted/deferred).
#: round-56 F4: this is an ENFORCED closed vocabulary -- :func:`register_candidate`
#: coerces any bucket name outside it to UNCLASSIFIED (with a versioned
#: anomaly line), so it is no longer merely an assumed/documented contract.
#: When a single (candidate, census) pair is sighted with more than one
#: DISTINCT bucket across multiple registrations, :func:`_most_conservative_bucket`
#: picks the one appearing LATEST here -- never the last-write, so the
#: resolved bucket (and therefore ``effective_bucket``) does not depend on
#: sighting order. See the module docstring's "SAME-CENSUS RESIGHT
#: RESOLUTION" paragraph for the full rationale.
BUCKET_DECLARATION_ORDER = (
    "ADMISSIBLE_A", "ADMISSIBLE_B",
    "DEFERRED_A", "DEFERRED_B", "DEFERRED_C", "DEFERRED_D", "DEFERRED_E",
    "DEFERRED_F",
    "FEATURE_RADICAL", "CONDUIT_ECHO",
    "UNCLASSIFIED",
)
_BUCKET_CONSERVATISM_RANK = {b: i for i, b in enumerate(BUCKET_DECLARATION_ORDER)}


def _most_conservative_bucket(buckets):
    """Resolve a set of buckets sighted under the SAME census for the SAME
    candidate to a single deterministic bucket, per the fixed total order
    :data:`BUCKET_DECLARATION_ORDER`: the bucket appearing LATEST in that
    order wins. A set of one element returns that element unchanged (the
    common case, and the only case for a census sighted just once).

    An unrecognized bucket string is now ENFORCED out at registration
    (round-56 F4): :func:`register_candidate` coerces any name outside the
    closed :data:`BUCKET_DECLARATION_ORDER` vocabulary to UNCLASSIFIED and
    emits a versioned anomaly line, so an unknown string can no longer reach
    this function through the public API. The defensive after-everything
    rank below is therefore pure defense-in-depth: were an unknown name to
    arrive anyway it sorts after everything (never losing to a known,
    less-conservative bucket), with lexicographic tie-breaking.

    round-55 M2 (determinism defect fix): known buckets have pairwise
    distinct ranks, but every UNKNOWN bucket string shares the defensive
    after-everything rank -- under the old single-component key, max() over
    two unknown names returned whichever the set yielded first, and str
    hash randomization makes set iteration order vary across processes.
    Rank ties therefore break lexicographically (largest name wins), so
    the result is a pure function of the SET, never of iteration order."""
    return max(buckets,
              key=lambda b: (_BUCKET_CONSERVATISM_RANK.get(
                  b, len(BUCKET_DECLARATION_ORDER)), str(b)))


# ---------------------------------------------------------------------------
# Species-level classification (pure core)
# ---------------------------------------------------------------------------

def is_pool(species):
    """True if the species is a polymer pool participant (by label)."""
    label = species.get("label") or ""
    return any(label.startswith(p) for p in POOL_LABEL_PREFIXES)


def is_pool_state_resolvable(species):
    """True for condensed chips that resolve to a pool state (by label).

    round-53: consults the VALIDATED ACTIVE label set
    (:func:`active_pool_state_labels`), never the raw DECLARED set
    (:data:`POOL_STATE_RESOLVABLE_LABELS`) -- a declared label is not
    ground truth until the in-repo isomorphism oracle has validated it
    (see the module's "DECLARED vs ACTIVE pool-state label oracle"
    section)."""
    label = species.get("label") or ""
    return label in active_pool_state_labels()


def is_chain_scale(species):
    """Chain-scale discrete: >= 2.5 monomer-equiv on mass AND heavy atoms."""
    mw = species.get("mw")
    heavy = species.get("heavy_atoms")
    if mw is None or heavy is None:
        return False
    return mw >= CHAIN_SCALE_MW and heavy >= CHAIN_SCALE_HEAVY


def species_role(species):
    """round-53 crash fix: ``pool_role_override`` (set by
    :func:`_apply_iso_overrides` from a live per-species isomorphism
    verdict) is an explicit BOOLEAN and is consulted ahead of the
    ACTIVE-set label test -- it can assign POOL even with an EMPTY
    ACTIVE set (an iso=True discrete whose label was never declared, or
    was dropped), and it can equally veto a stale ACTIVE-set label match
    (iso=False mid-run drift) by falling through to the size-based role
    test instead. This replaces the old surrogate-label mechanism
    (``label_for_roles = POOL_STATE_RESOLVABLE_LABELS[0]``), which raised
    IndexError whenever the resolvable-label tuple was empty."""
    # round-55 P2(b), reviewed and kept: is_pool() is consulted BEFORE
    # pool_role_override, deliberately. The pool-label-prefix test is the
    # pure core's stand-in for isinstance(..., Polymer) -- a STRUCTURAL
    # identity of the participant, not an oracle verdict -- so the
    # oracle-layer override (which exists solely to arbitrate
    # DECLARED-label vs isomorphism disagreements on discrete chips) never
    # gets to veto it. The in-repo adapter cannot even produce that
    # conflict: _species_entry returns Polymer participants before the
    # oracle/divergence machinery runs, so no pool participant ever
    # carries an override; only a hand-built core record could combine
    # the two, and pool identity should still win there.
    if is_pool(species):
        return POOL
    override = species.get("pool_role_override")
    if override is True:
        return POOL
    if override is False:
        pass  # explicit non-pool verdict wins over a stale ACTIVE match
    elif is_pool_state_resolvable(species):
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

    # P2a fix: conduit_echo census, documented at the top of this module
    # as a valid ``record["census"]`` value, was previously unhandled here
    # and fell through to r93-style shape classification (or UNCLASSIFIED)
    # instead of the CONDUIT_ECHO bucket the docstring promises. Mirrors
    # the feature_radical branch above: a pure-core caller passing a
    # conduit_echo record gets CONDUIT_ECHO unconditionally, never a shape
    # verdict an echo sighting carries no information to support.
    if record["census"] == "conduit_echo":
        result["bucket"] = "CONDUIT_ECHO"
        result["refusal_reason"] = REASON_CONDUIT_ECHO
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

    CONDUIT-ECHO RANKS LOWEST (round-50 determinism finding): below
    feature_radical, the censuses are NOT last-write-wins -- 'conduit_echo'
    is an ordered LOWEST-precedence rank, so any other census sighted on
    the same key (chiefly 'r93_general') always beats it for
    ``effective_bucket``, in EITHER registration order. Without this, a key
    seen r93_general-then-conduit_echo and one seen
    conduit_echo-then-r93_general would resolve to different
    effective_bucket values for what is otherwise an identical sighting
    set -- an order dependency across rebuild epochs that violates this
    module's determinism contract. See :meth:`_effective_bucket_below_fr`.

    SAME-CENSUS RESIGHT RESOLUTION (round-50 P1, determinism finding): the
    per-census value in ``bucket_by_census`` is not last-write-wins either.
    Internally each census tracks the FULL SET of buckets ever sighted for
    it on this candidate (``bucket_sightings_by_census``); the exposed
    ``bucket_by_census[census]`` is that set resolved via
    :func:`_most_conservative_bucket`, so it -- and everything downstream
    of it, including ``effective_bucket`` -- is a pure function of the
    sighting set, independent of registration order or which rebuild
    epoch a sighting landed in. See the module docstring's "SAME-CENSUS
    RESIGHT RESOLUTION" paragraph.
    """

    def __init__(self):
        self._lock = threading.Lock()
        self._entries = {}

    def register(self, key, census, bucket, epoch=None):
        with self._lock:
            entry = self._entries.setdefault(
                key, {"censuses": set(), "bucket_by_census": {},
                      "bucket_sightings_by_census": {},
                      "effective_bucket": None, "shadow_bucket": None,
                      "precedence": None, "epochs": set()})
            entry["censuses"].add(census)
            # M18.3 ledger epochs (DESIGN §3.3): the engine rebuild-signature
            # id (or RMG iteration int) under which this sighting was
            # recorded. Bookkeeping only -- FR stickiness stays set-membership
            # of "feature_radical" in censuses, across ALL epochs of the run.
            #
            # m18.4 §1.3: DELIBERATE CONTRACT CHANGE -- epoch=None used to
            # leave entry["epochs"] untouched (skip). It now falls back to
            # the process-wide CORE EPOCH provider's current token
            # (:func:`current_epoch`): every production caller passes
            # epoch=None today (design §INVERSION-1), so this single
            # fallback wires ALL of them (annotate_refused_row,
            # annotate_feature_radical) with zero edits to polymer.py. An
            # explicit epoch= argument still overrides the fallback.
            if epoch is None:
                epoch = current_epoch()
            entry["epochs"].add(epoch)
            # round-71 P1 finding: mark this token SIGHTED on the process-
            # wide CORE EPOCH provider -- whether ``epoch`` came from the
            # fallback above or was passed explicitly (see
            # :meth:`_EpochProvider.note_sighted` for the explicit-token
            # decision). This is what lets a LATER rebuild failure on the
            # SAME epoch tell "attempted, never sighted" (burn) apart from
            # "attempted, but this row was sighted first" (NOT a burn --
            # see :meth:`_EpochProvider.note_conduit_rebuild_failed`).
            # Never raises; cheap (one set add under the provider's own
            # lock, a different lock than this one -- no cross-singleton
            # deadlock risk).
            _EPOCH_PROVIDER.note_sighted(epoch)
            # round-73 P2 (virgin-guard false negative, design §2.5): every
            # registration -- admit AND deny -- is census activity for the
            # standing-admit registry's own virgin-lifecycle guard, not
            # only a would-admit insert (see
            # _StandingAdmitRegistry.note_lifecycle_activity). Forward
            # reference to _STANDING_ADMITS is safe: this method only ever
            # runs after full module import, by which point the singleton
            # below is defined (same late-binding already relied on for
            # _EPOCH_PROVIDER.note_sighted above).
            _STANDING_ADMITS.note_lifecycle_activity()
            sightings = entry["bucket_sightings_by_census"].setdefault(
                census, set())
            sightings.add(bucket)
            # round-50 P1: resolve the (possibly multi-valued) sighting set
            # for this census deterministically -- never last-write-wins.
            entry["bucket_by_census"][census] = _most_conservative_bucket(
                sightings)
            if "feature_radical" in entry["censuses"]:
                entry["effective_bucket"] = "FEATURE_RADICAL"
                entry["shadow_bucket"] = entry["bucket_by_census"].get(
                    "r93_general")
                entry["precedence"] = (
                    PRECEDENCE_RULE if len(entry["censuses"]) > 1 else None)
            else:
                entry["effective_bucket"] = self._effective_bucket_below_fr(
                    entry)
                entry["shadow_bucket"] = None
                entry["precedence"] = None
            return self._copy_entry(entry)

    @staticmethod
    def _effective_bucket_below_fr(entry):
        """Ordered precedence among the non-feature_radical censuses seen
        on ``entry`` (round-50 determinism finding): 'conduit_echo' is an
        explicit LOWEST-precedence rank -- any OTHER census sighted on this
        key (today, only 'r93_general') always wins for
        ``effective_bucket``, regardless of which was registered first or
        last. A key seen ONLY via conduit_echo falls back to its own echo
        bucket. Expressed as a rank over ``entry["censuses"]`` rather than
        a last-write special case so a future non-echo census slots in
        without new branching."""
        non_echo = entry["censuses"] - {CONDUIT_ECHO_CENSUS}
        if non_echo:
            # Deterministic tie-break by census name; today exactly one
            # non-echo, non-FR census exists ('r93_general').
            winner = sorted(non_echo)[0]
            return entry["bucket_by_census"][winner]
        return entry["bucket_by_census"][CONDUIT_ECHO_CENSUS]

    @staticmethod
    def _copy_entry(entry):
        """Defensive copy: callers must never be able to mutate internal
        registry state through a returned entry. EVERY nested mutable
        container is duplicated: the top-level dict, the censuses/epochs
        sets, the ``bucket_by_census`` dict (round-55 P1-1: this one was
        previously returned LIVE, so register/lookup callers could reach
        registry internals through the "copy"), and each per-census
        sighting set. All remaining leaf values are strings or None
        (immutable), so the returned view is fully detached."""
        return dict(
            entry,
            censuses=set(entry["censuses"]),
            epochs=set(entry["epochs"]),
            bucket_by_census=dict(entry["bucket_by_census"]),
            bucket_sightings_by_census={
                c: set(b) for c, b in entry["bucket_sightings_by_census"].items()
            },
        )

    def lookup(self, key):
        with self._lock:
            entry = self._entries.get(key)
            if entry is None:
                return None
            return self._copy_entry(entry)

    def counts(self):
        """Return ``(eff, overlap, total, resight_divergence)``, ALL
        recomputed fresh from the current final registry state (round-50
        P1 fix #5): never incrementally accumulated during ingestion, so
        none of these can desync from the sightings that produced them,
        and all are permutation-of-sightings independent."""
        with self._lock:
            eff = Counter(e["effective_bucket"] for e in self._entries.values())
            eff.setdefault("UNCLASSIFIED", 0)
            overlap = sum(1 for e in self._entries.values()
                          if len(e["censuses"]) > 1)
            resight_divergence = sum(
                1
                for e in self._entries.values()
                for sightings in e["bucket_sightings_by_census"].values()
                if len(sightings) > 1
            )
            return eff, overlap, len(self._entries), resight_divergence

    @property
    def unclassified_total(self):
        """Candidates whose current ``effective_bucket`` is UNCLASSIFIED,
        computed fresh from final registry state on every access (round-50
        P1 fix #5) -- never an incrementally-maintained counter, so it can
        never desync from a same-census resighting that changes a
        candidate's resolved bucket after its first sighting."""
        with self._lock:
            return sum(1 for e in self._entries.values()
                      if e["effective_bucket"] == "UNCLASSIFIED")

    def reset(self):
        with self._lock:
            self._entries.clear()


#: The process-wide registry instance.
CENSUS_REGISTRY = _CensusRegistry()

#: round-60 P2-3: single module-level lock guarding every check-then-add
#: sequence against the module-level warn-once sets below
#: (_bucket_anomaly_warned, _admission_token_anomaly_warned). These
#: previously did an UNLOCKED "if key not in set: add; log" -- two
#: concurrent threads could both pass the membership check before either
#: added the key, so both would log (double emission). reset_conduit_state()
#: clears these same sets under this lock too (:func:`_warn_once_clear`), so
#: a lifecycle reset can never race a concurrent check-then-add in
#: :func:`_warn_once` and leave the set in a state where the loser's key
#: silently vanishes without ever logging.
_WARN_ONCE_LOCK = threading.Lock()


def _warn_once(warned_set, key, message, *args):
    """Thread-safe check-then-add against a module-level warn-once set
    (round-60 P2-3): adds ``key`` to ``warned_set`` and logs ``message`` %
    ``args`` at most once per key, with the check-then-add itself guarded by
    :data:`_WARN_ONCE_LOCK` so concurrent callers can never both win the
    race. Returns True iff this call was the one that logged."""
    with _WARN_ONCE_LOCK:
        if key in warned_set:
            return False
        warned_set.add(key)
    logging.warning(message, *args)
    return True


def _warn_once_clear(*warned_sets):
    """Clear one or more warn-once sets atomically under
    :data:`_WARN_ONCE_LOCK` (round-60 P2-3), so a lifecycle reset can never
    interleave with a concurrent :func:`_warn_once` check-then-add."""
    with _WARN_ONCE_LOCK:
        for s in warned_sets:
            s.clear()


#: round-58 P2-B: warn-once dedup state for the unknown-bucket anomaly line
#: below, keyed by ``(census, bucket)``. Mirrors the module-level warn-once
#: set pattern used throughout :mod:`rmgpy.polymer` (e.g.
#: ``_refused_census_warned``, ``_flux_archetype_warned``): a hot
#: coercion producer (a persistently-mistyped bucket name feeding
#: :func:`register_candidate` on every registration) would otherwise spam
#: one log line PER REGISTRATION. The underlying COERCION below still
#: applies on EVERY call -- only the LOG LINE is deduplicated. Cleared by
#: :func:`reset_conduit_state` so a new lifecycle re-arms the line. Checked
#: and set under :data:`_WARN_ONCE_LOCK` via :func:`_warn_once` (round-60
#: P2-3).
_bucket_anomaly_warned = set()


def register_candidate(key, census, bucket, epoch=None):
    """Record a candidate sighting; returns the resolved registry entry.

    round-56 F4 (unknown-bucket domination fix): this public entry point is
    where the bucket vocabulary is ENFORCED. :data:`BUCKET_DECLARATION_ORDER`
    is a closed vocabulary; a bucket name outside it (e.g. a typo) would --
    absent this check -- reach :func:`_most_conservative_bucket`, where its
    defensive after-everything rank sorts it ABOVE every known bucket
    (including UNCLASSIFIED) and would silently dominate the correct
    conservative classification. Instead, an out-of-vocabulary bucket is
    COERCED to UNCLASSIFIED and a loud versioned anomaly line is emitted, so
    no unknown string ever reaches the registry. Never raises (census-only
    contract).

    round-58 P2-B: the coercion itself runs on EVERY call (registry state
    must stay correct every time); the anomaly LOG LINE is deduplicated to
    once per ``(census, bucket)`` pair per lifecycle (:data:`_bucket_anomaly_warned`)
    so a hot mistyped-bucket producer cannot spam the log once per
    registration."""
    if bucket not in _BUCKET_CONSERVATISM_RANK:
        _warn_once(
            _bucket_anomaly_warned, (census, bucket),
            "CONDUIT CENSUS BUCKET ANOMALY/1 census=%s bucket=%s "
            "action=coerced-unclassified", census, bucket)
        bucket = "UNCLASSIFIED"
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
    re-block. Also clears the run-level conduit flux accumulator (§4.4) and
    the round-53 DECLARED/ACTIVE label-oracle ledger (:data:`_LABEL_ORACLE`)
    so a fresh run re-validates every declared label from scratch.

    round-55 P1-2: this reset is also the label oracle's LIFECYCLE
    boundary -- before its ledger is cleared, the outgoing lifecycle is
    finalized (never-sighted declared labels resolve as missing-species)
    and its one-shot ``CONDUIT CLASSIFIER ORACLE HEALTH/1`` line is
    emitted, exactly once per epoch-reset cycle, even when zero census
    rows fired in it (the virgin import->first-reset interval is not a
    lifecycle and emits nothing; the process-exit atexit hook closes the
    final lifecycle of the process).

    Import-safe: the rmgpy.polymer sets are cleared lazily so this module
    keeps its no-module-level-RMG-imports contract.

    round-73 P2 fix (reset ordering zeroed epoch counters in standing
    health): ALL THREE outgoing-lifecycle closes (oracle, epoch, standing)
    must run and emit with LIVE counts BEFORE any of the three singletons'
    state is cleared -- the standing registry's own close reads
    ``_EPOCH_PROVIDER.attempted_and_burned()`` (design §2.5), so if the
    epoch provider's state were cleared first (as a naive
    ``_EPOCH_PROVIDER.reset()`` before ``_STANDING_ADMITS.reset()`` would
    do), a lifecycle closed via this reset would always report
    ``epochs_seen=0 burned=0`` regardless of real prior activity. Fixed by
    literally reusing :func:`close_conduit_lifecycle` (the same
    three-surface, each-under-its-own-never-raise-guard close helper) as
    the CLOSE-ONLY phase here, before any clearing happens; each
    singleton's own ``.reset()`` below then finds its lifecycle already
    closed (its internal ``close_lifecycle()`` call is idempotent, a
    no-op) and only CLEARS state -- so reset and an explicit
    :func:`close_conduit_lifecycle` call can never drift apart in
    behavior.
    """
    # Close all three lifecycles FIRST, with state still live, so the
    # standing registry's health line reads real (pre-clear) epoch
    # counts -- see the docstring note above. Each surface is already
    # isolated under its own never-raise guard inside
    # close_conduit_lifecycle() itself.
    close_conduit_lifecycle()
    # Now clear everything together ("reset both or neither", m18.4
    # §1.2/design §2.5). Each singleton's .reset() calls its own
    # close_lifecycle() again first -- a no-op here since the close above
    # already ran -- then clears state.
    _LABEL_ORACLE.reset()
    _EPOCH_PROVIDER.reset()
    CENSUS_REGISTRY.reset()
    _STANDING_ADMITS.reset()
    _CONDUIT_FLUX_TOTALS.clear()
    # round-58 P2-B: the unknown-bucket and F3 token-anomaly warn-once
    # dedup sets are part of the same "reset both or neither" contract --
    # a new lifecycle must re-arm these lines, not inherit an already-fired
    # dedup key from the previous run. round-60 P2-3: cleared atomically
    # under _WARN_ONCE_LOCK (via _warn_once_clear) so this reset can never
    # race a concurrent _warn_once() check-then-add on either set.
    # round-70 P2-e: the advance-after-close anomaly dedup set joins the
    # same both-or-neither contract, for the same reason. round-71 P1: the
    # sighting-on-burned-epoch anomaly dedup set joins it too.
    _warn_once_clear(_bucket_anomaly_warned, _admission_token_anomaly_warned,
                     _epoch_advance_after_close_warned,
                     _epoch_sighting_on_burned_warned)
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
    """One-line loud census summary. ALWAYS names UNCLASSIFIED, zero or not.

    round-50 P1: APPENDS a ``resight_divergence=<n>`` token after the
    existing tokens -- the count of (candidate, census) pairs that were
    ever sighted with more than one distinct bucket (see the module
    docstring's "SAME-CENSUS RESIGHT RESOLUTION" paragraph). This is a
    pure append: every existing token/ordering in this line is unchanged.

    round-55 P1-3: READ-ONLY. This used to be the "summary time"
    finalization hook for the lazy DECLARED/ACTIVE label oracle, which
    made any EARLY (diagnostic) summary call permanently resolve every
    not-yet-sighted declared label as ``missing-species`` -- poisoning
    legitimate later sightings with the cached dropped verdict. A summary
    now reports current state without finalizing anything; finalization
    lives at the explicit end-of-lifecycle boundaries (the epoch reset
    :func:`reset_conduit_state` and the process-exit atexit hook; see
    :func:`finalize_label_oracle` and :class:`_LabelOracleState`)."""
    eff, overlap, total, resight_divergence = CENSUS_REGISTRY.counts()
    parts = " ".join(f"{b}={eff[b]}" for b in sorted(eff))
    return (f"conduit-census/1 summary: candidates={total} overlap={overlap} "
            f"{parts} resight_divergence={resight_divergence}")


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

    Round-36 P1(c) divergence rule, reshaped by round-53: the census-side
    LABEL test for pool-state resolvability is now the VALIDATED ACTIVE
    set (never the raw DECLARED set -- see the module's "DECLARED vs
    ACTIVE pool-state label oracle" section), compared against the
    in-repo ISOMORPHISM test
    (:func:`rmgpy.polymer._discrete_resolves_to_pool_state`) run fresh on
    THIS row's own species/pools. A DECLARED label is validated once, at
    its first sighting here (:meth:`_LabelOracleState.note_sighting`), so
    that same-row validation and the divergence check can never disagree
    on the very sighting that establishes the label's ACTIVE-set
    membership; a genuine divergence now means the ACTIVE-set verdict and
    THIS row's structural verdict disagree -- either real mid-run drift,
    or a species that is genuinely isomorphic under a label that was
    never declared (or was dropped). On divergence the mismatch is
    census-logged (versioned token: CONDUIT CLASSIFIER DIVERGENCE/1) and
    flagged on the entry, and the ISOMORPHISM verdict is used -- neither
    side is silently overridden.
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

    # round-53: validate a DECLARED label at its first sighting (lazy --
    # species/pools become known to this census layer per-row, not at a
    # single init point), then compare the resulting ACTIVE-set membership
    # against THIS row's own structural isomorphism verdict.
    _LABEL_ORACLE.note_sighting(label, species, row_pools)
    label_says = label in _LABEL_ORACLE.active_labels()
    try:
        iso_says = bool(_discrete_resolves_to_pool_state(species, row_pools))
    except Exception:  # defensive: census-only code must not raise
        iso_says = False
    if label_says != iso_says:
        entry["_divergence"] = True
        logging.warning(
            "CONDUIT CLASSIFIER DIVERGENCE/1 species=%s label=%s "
            "label_says=%d iso_says=%d active_label=%d action=use-iso",
            token, label, int(label_says), int(iso_says), int(label_says))
        # Make the core's role assignment follow the isomorphism verdict via
        # the explicit boolean override field (round-53 crash fix; see
        # :func:`_apply_iso_overrides` / :func:`species_role`) -- never a
        # surrogate label swap.
        entry["_iso_pool_state"] = iso_says
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
    """Fold the adapter's per-species isomorphism verdict into the core's
    role test via an explicit BOOLEAN override field (round-53 crash fix):
    entries carrying an ``_iso_pool_state`` hint (set by
    :func:`_species_entry` exactly on divergence) get ``pool_role_override``
    set to that hint, consumed by :func:`species_role` ahead of the
    ACTIVE-set label test.

    This replaces the old surrogate-label mechanism
    (``label_for_roles = POOL_STATE_RESOLVABLE_LABELS[0]``), which raised
    IndexError whenever the resolvable-label tuple was empty -- exactly
    the state an all-dropped DECLARED set now reaches routinely. Species
    with no hint (the overwhelming common case: no divergence) are left
    untouched."""
    for side in (record["reactants"], record["products"]):
        for s in side:
            hint = s.pop("_iso_pool_state", None)
            if hint is not None:
                s["pool_role_override"] = hint
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
        # design §2.2/§2.3 (commit 3, ruling D): populate the standing-
        # admit registry from the FINAL verdict (post FR-overlap
        # downgrade above), so a would-admit row that just lost its
        # admission to an FR overlap is never registered as standing.
        register_standing_admit_would_admit(forward, row_pools, record,
                                            result, verdict)
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
#: exception type name). 'echo-not-evaluated' (round-50 finding, contract
#: fix) is the label-only :func:`annotate_feature_radical` hook's own
#: reason: a pure conduit_echo sighting never carries a live reaction
#: object, so G0/G2-G7 never run there -- fail-closed default-deny stays,
#: with an HONEST reason (never fabricate "feature-radical-overlap" for a
#: key G1 does not actually block). 'admission-token-anomaly' (round-56 F3)
#: is the fail-closed fallback emitted when a reserved token (currently
#: 'echo-not-evaluated') is misused: :func:`admission_census_suffix`
#: deny-by-defaults with this reason instead of propagating the misused
#: token, and logs a versioned CONDUIT ADMISSION TOKEN ANOMALY/1 line.
ADMISSION_DENY_REASONS = frozenset({
    "classifier-not-admissible", "classifier-divergence",
    "feature-radical-overlap", "direction-inadmissible",
    "direction-requires-flip-rewrite", "gas-product-count",
    "gas-mw-threshold-unresolvable", "gas-mw-over-threshold",
    "landing-cone-violation", "destination-unresolvable",
    "chain-not-condensed", "not-balanced", "kinetics-not-exportable",
    "kinetics-not-yet-assigned", "echo-not-evaluated",
    "admission-token-anomaly",
})

#: round-58 P2-B: warn-once dedup state for the F3 token-anomaly lines in
#: :func:`admission_census_suffix`, keyed by ``(candidate_key, reason)``.
#: Same rationale/pattern as :data:`_bucket_anomaly_warned`: the fail-closed
#: substitution runs on EVERY misuse, only the log line is deduplicated.
#: Cleared by :func:`reset_conduit_state` so a new lifecycle re-arms it.
_admission_token_anomaly_warned = set()

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


def _sanitize_dynamic_token(raw, max_len=64):
    """Shared charset-safe sanitizer for any dynamic fragment (chiefly
    ``type(exc).__name__``) entering a ``/1`` protocol key=value line
    (round-60 P2-1 origin, round-73 P3 extraction): collapses every
    character outside the strict ``[A-Za-z0-9_.-]+`` charset every static
    token in this module already uses to a single ``_``, falling back to
    the literal ``"unsanitizable"`` if that would otherwise leave nothing,
    then caps the result at ``max_len`` characters (round-66 P2: an
    unbounded sanitized token is itself a log/line-length DoS vector).
    ``max_len=None`` skips the cap -- used by :func:`admission_census_suffix`,
    which needs the FULL sanitized string to compute its own
    dirty-vs-truncated ``action=`` classification before applying its own
    64-char cap.

    This is the single shared CORE both :func:`admission_census_suffix`'s
    richer ``admission-evaluation-error:*`` audit trail (which additionally
    records raw_sha/raw_b64/dirty-vs-truncated ``action=`` for that
    specific deny-reason family) and every plainer dynamic-exception-name
    site in this module (the standing-admit/epoch/lifecycle-close anomaly
    lines) route through -- so no value can ever reach a ``/1`` line
    carrying embedded whitespace, ``=``, or any other charset-breaking
    character, regardless of source. NEVER raises."""
    try:
        sanitized = re.sub(r"[^A-Za-z0-9_.-]", "_", str(raw)) or "unsanitizable"
        return sanitized if max_len is None else sanitized[:max_len]
    except Exception:  # pragma: no cover - defensive fail-closed
        return "unsanitizable"


def _raw_sha_suffix(raw, sanitized):
    """Minimal forensic audit trail for the plain sanitized-exception-name
    anomaly sites in this module (round-74 P3): when two distinct hostile
    exception CLASS names collapse to the SAME sanitized token under
    :func:`_sanitize_dynamic_token` (e.g. ``"Bad Name"`` and ``"Bad=Name"``
    both -> ``"Bad_Name"``), the ``/1`` line otherwise loses forensic
    distinguishability between them. This deliberately does NOT build the
    full raw_sha/raw_b64/action machinery :func:`admission_census_suffix`
    carries for its ``admission-evaluation-error:*`` family (round-64) --
    that machinery exists for a different threat model (attacker/data-
    influenced raw suffixes). This case is bounded, code-controlled
    exception class names, so it stays minimal: a `` raw_sha=<12 hex>``
    fragment (the first 12 hex chars of the sha256 of the raw string,
    encoded with ``errors="backslashreplace"`` exactly like
    :func:`admission_census_suffix`'s own raw digest), computed and
    returned ONLY when sanitization actually changed the string
    (``sanitized != str(raw)``) -- the clean/common path (no sanitization
    needed) emits nothing extra, avoiding noise. NEVER raises."""
    try:
        raw_str = str(raw)
        if sanitized == raw_str:
            return ""
        digest = hashlib.sha256(
            raw_str.encode("utf-8", errors="backslashreplace")).hexdigest()
        return f" raw_sha={digest[:12]}"
    except Exception:  # pragma: no cover - defensive fail-closed
        return ""


def admission_census_suffix(verdict):
    """APPEND-ONLY census tokens carrying the would-be admission verdict
    while :data:`CONDUIT_ADMISSION_ENABLED` is False (BUILD_SPEC W1.6):
    an ADMISSIBLE verdict appends ``would_admit=1 deny=None stage=final
    rewrite=<bool>``; a denied verdict appends ``deny=<reason>
    stage=<provisional|final>``. The tokens ride their own bracketed group
    AFTER the M18.2 ``[conduit-census/1 ...]`` suffix, so the existing
    header/tokens (and the line's closing bracket) stay byte-identical
    (round-36 P1(b)).

    ``stage=`` (round-49): the G6 re-adjudication hook makes the same
    candidate emit TWO census lines -- one at stamp time (provisional
    ``kinetics-not-yet-assigned``) and one after kinetics are final.
    Without a machine-readable discriminator, grep tallies double-count
    and naive consumers retain the provisional deny as a real denial.
    ``stage=provisional`` marks exactly the verdicts in
    :data:`PROVISIONAL_DENY_REASONS` (subject to re-adjudication); every
    other line -- every admit and every one-shot deny -- is the final
    word for its emission and says so with ``stage=final``. Count/filter
    with ``stage=final``.

    P2b fix (echo-token reservation, fail-closed): this is the ONE
    chokepoint every caller (:func:`annotate_refused_row` and
    :func:`annotate_feature_radical`) routes a verdict through before it
    is serialized, so it is where ``echo-not-evaluated`` is reserved.
    That reason means "this key was never evaluated because the sighting
    is a pure conduit_echo with no live reaction object" -- it is a lie
    if the SAME candidate key also carries any non-echo census
    membership (that key WAS genuinely evaluated, under a different
    census), and it also describes a RECORDED echo sighting, so a key with
    NO ledger entry at all (never registered under ANY census) misuses it
    too. Either misuse is a caller bug (the only legitimate producer,
    :func:`annotate_feature_radical`, registers an echo-only sighting
    before building the verdict).

    round-56 F3 (never-raise-contract alignment): both misuses previously
    RAISED :class:`ValueError`. That exception, escaping into the
    census-only production wrappers, was swallowed into a generic
    ``annotation-failed`` line -- losing the specific invariant info AND
    violating this module's never-raise contract (census/bookkeeping code
    must never raise into generation paths). Each is now surfaced as a
    loud VERSIONED anomaly line (``CONDUIT ADMISSION TOKEN ANOMALY/1``)
    naming the key, the token, the specific reason, and the census
    memberships, and the suffix falls back to a fail-closed deny that does
    NOT propagate the reserved ``echo-not-evaluated`` token (deny-by-
    default, but never emit the token that was misused).

    round-58 P2-B: both F3 anomaly lines below are rate-limited to once per
    ``(candidate_key, reason)`` per lifecycle (:data:`_admission_token_anomaly_warned`),
    mirroring the same warn-once dedup pattern used for the unknown-bucket
    anomaly in :func:`register_candidate` -- the fail-closed substitution
    itself still happens on EVERY call, only the log line is deduplicated.

    round-58 P2-C: a denied verdict's ``deny_reason`` is validated against
    the closed :data:`ADMISSION_DENY_REASONS` vocabulary (plus the
    dynamically-suffixed ``admission-evaluation-error:*`` G7 family) before
    it is serialized, so a future typo'd reason can never leak as a
    structured token into the census output: an out-of-vocabulary reason is
    surfaced with a loud ``CONDUIT ADMISSION TOKEN ANOMALY/1`` line and
    substituted with ``admission-token-anomaly`` -- the same conservative,
    in-vocabulary fail-closed token this function already uses for the
    echo-token misuses above.

    round-60 P2-1 (token-injection hole): the dynamic
    ``admission-evaluation-error:*`` family (G7's fail-closed catch-all,
    suffixed with the failing exception's TYPE NAME) deliberately bypasses
    the closed-vocabulary check above -- it is not a fixed literal. That
    made it the one path that serialized an attacker/data-influenced string
    into the census line WITHOUT going through any charset validation: a
    dynamically-constructed exception class with a crafted ``__name__``
    could embed spaces, ``=``, or brackets and forge what looks like
    additional ``key=value`` tokens in the line. This chokepoint now
    sanitizes that dynamic suffix -- same strict charset every other
    already-safe token value in this module uses -- before it can ever be
    serialized (see the sanitization block below), so no value can ever
    produce a token with embedded whitespace or ``=`` regardless of input.

    round-60 P2-3: the anomaly lines in this function (echo-token misuse,
    dynamic-token sanitization, unknown-deny-reason) are all rate-limited
    via :func:`_warn_once`, which guards the check-then-add against
    :data:`_admission_token_anomaly_warned` with :data:`_WARN_ONCE_LOCK` so
    concurrent callers can never both win the race and double-emit."""
    if verdict is None:
        return ""
    if verdict.deny_reason == "echo-not-evaluated":
        entry = lookup_candidate(verdict.candidate_key)
        if entry is None:
            # F3: unregistered-key echo token -- loud versioned anomaly,
            # fail-closed deny WITHOUT the reserved token.
            _warn_once(
                _admission_token_anomaly_warned,
                (verdict.candidate_key, "unregistered-key"),
                "CONDUIT ADMISSION TOKEN ANOMALY/1 key=%s "
                "token=echo-not-evaluated reason=unregistered-key "
                "action=fail-closed-deny", verdict.candidate_key)
            verdict = _deny(verdict.candidate_key, "admission-token-anomaly")
        else:
            non_echo_censuses = entry["censuses"] - {CONDUIT_ECHO_CENSUS}
            if non_echo_censuses:
                # F3: mixed-membership echo token -- kept EQUALLY loud and
                # versioned (not downgraded), fail-closed deny WITHOUT the
                # reserved token.
                _warn_once(
                    _admission_token_anomaly_warned,
                    (verdict.candidate_key, "mixed-census-membership"),
                    "CONDUIT ADMISSION TOKEN ANOMALY/1 key=%s "
                    "token=echo-not-evaluated "
                    "reason=mixed-census-membership censuses=%s "
                    "action=fail-closed-deny",
                    verdict.candidate_key,
                    "+".join(sorted(non_echo_censuses)))
                verdict = _deny(verdict.candidate_key,
                                "admission-token-anomaly")
    if (verdict.deny_reason or "").startswith("admission-evaluation-error:"):
        # P2-1: sanitize the dynamic G7 suffix at the serialization
        # chokepoint. Only a strict [A-Za-z0-9_.-]+ charset survives --
        # matching the charset every static deny-reason literal in
        # ADMISSION_DENY_REASONS already satisfies -- so this token can
        # never carry embedded whitespace, '=', or brackets that a naive
        # downstream parser could mistake for other structured tokens.
        raw_suffix = verdict.deny_reason[len("admission-evaluation-error:"):]
        # round-73 P3: routed through the shared _sanitize_dynamic_token
        # core (same charset regex, same fallback) rather than a private
        # copy of the regex -- see that function's docstring.
        sanitized_full = _sanitize_dynamic_token(raw_suffix, max_len=None)
        was_dirty = sanitized_full != raw_suffix
        # round-66 P2 (Finding 1): the sanitized suffix above was still
        # UNBOUNDED -- a huge raw_suffix (whether charset-dirty or already
        # charset-clean) produced an equally huge sanitized string that was
        # both logged verbatim in the anomaly line's `sanitized=` field AND
        # emitted as the final census deny token, so the length-cap half of
        # the round-65 DoS fix was never actually wired up. Cap the
        # sanitized form at 64 chars and use the IDENTICAL truncated string
        # for both the census token and the `sanitized=` field -- they must
        # never diverge. Truncating an already charset-clean string
        # trivially still matches [A-Za-z0-9_.-]+ (tested).
        was_long = len(sanitized_full) > 64
        sanitized = sanitized_full[:64] if was_long else sanitized_full
        if was_dirty or was_long:
            # round-64 P2 (sanitized deny reasons silently collapse distinct
            # causes): the census TOKEN stays exactly the sanitized value
            # above (cardinality/vocabulary unchanged), but distinct raw
            # suffixes that happen to sanitize to the same token (e.g.
            # "bad token" / "bad=token" / "bad_token") would otherwise be
            # indistinguishable, and the raw value was never recorded
            # anywhere auditable.
            #
            # round-65 P2 (Findings 1+2+3, joint fix): repr(raw_suffix) could
            # embed a literal space inside its quotes, breaking the one-line
            # whitespace-parseable key=value contract (F1); logging the full
            # raw with no cap is a log/memory DoS vector for a huge raw
            # suffix (F2); and a 12-hex (48-bit) dedup key risked suppressing
            # a genuinely distinct second raw on a 48-bit collision (F3).
            # Fixed by: hashing the FULL raw (never-raise via
            # errors="backslashreplace" so a malformed/surrogate-bearing raw
            # can never raise inside this census-only, never-raise path);
            # keying dedup on the FULL sha256 hexdigest (not the truncated
            # display form); and replacing the repr'd raw with two
            # whitespace-free bounded fields -- raw_len (full character
            # length) and raw_b64 (base64url, unpadded, of only the first 96
            # encoded bytes, truncated BEFORE base64 so the field itself is
            # always bounded -- round-66 P3 (Finding 3) switched this from
            # padded standard base64 to unpadded base64url so the field can
            # never contain a `=` that could confuse a naive
            # split-on-first-`=` key=value parser). The DISPLAYED raw_sha
            # stays the first 12 hex chars of the full digest, unchanged in
            # shape from round-64.
            #
            # round-66 P2 (Finding 1, continued): the anomaly now ALSO
            # fires when the only problem is over-length (charset-clean but
            # >64 chars) -- previously such a raw sailed through with no
            # anomaly line at all, leaving the operator with no way to
            # audit that a token had been silently truncated. `action=`
            # distinguishes the cases from a closed vocabulary:
            # `sanitized-truncated` (both problems), `truncated`
            # (over-length only), `sanitized` (charset-dirty only,
            # unchanged from round-65). The raw_sha (12-hex display of the
            # FULL sha256, still the dedup key) already disambiguates
            # collisions between distinct raws that truncate to the same
            # token, so no hash is appended to the token itself -- that
            # would change census cardinality/vocabulary shape, which
            # remains forbidden.
            raw_bytes = raw_suffix.encode("utf-8", errors="backslashreplace")
            raw_digest = hashlib.sha256(raw_bytes).hexdigest()
            raw_sha = raw_digest[:12]
            raw_len = len(raw_suffix)
            raw_b64 = base64.urlsafe_b64encode(
                raw_bytes[:96]).rstrip(b"=").decode("ascii")
            if was_dirty and was_long:
                action = "sanitized-truncated"
            elif was_long:
                action = "truncated"
            else:
                action = "sanitized"
            _warn_once(
                _admission_token_anomaly_warned,
                (verdict.candidate_key, raw_digest),
                "CONDUIT ADMISSION TOKEN ANOMALY/1 key=%s "
                "reason=unsanitized-error-token raw_sha=%s raw_len=%s "
                "raw_b64=%s sanitized=%s action=%s",
                verdict.candidate_key, raw_sha, raw_len, raw_b64, sanitized,
                action)
            verdict = _deny(verdict.candidate_key,
                            f"admission-evaluation-error:{sanitized}")
    if verdict.admitted:
        return (f" [conduit-admission/1 would_admit=1 deny=None "
                f"stage=final rewrite={verdict.needs_irreversible_rewrite}]")
    if (verdict.deny_reason not in ADMISSION_DENY_REASONS
            and not (verdict.deny_reason or "").startswith(
                "admission-evaluation-error:")):
        # P2-C: an unknown/typo'd deny reason must never leak as a
        # structured token -- surface it loudly and fail-closed with the
        # existing conservative in-vocabulary fallback. round-60 P3-1:
        # rate-limited (once per (candidate_key, reason) per lifecycle) so
        # a persistently-misbehaving producer cannot spam this line once
        # per call, mirroring every other TOKEN ANOMALY line in this
        # function.
        _warn_once(
            _admission_token_anomaly_warned,
            (verdict.candidate_key, "unknown-deny-reason", verdict.deny_reason),
            "CONDUIT ADMISSION TOKEN ANOMALY/1 key=%s "
            "reason=unknown-deny-reason value=%s action=fail-closed-deny",
            verdict.candidate_key, verdict.deny_reason)
        verdict = _deny(verdict.candidate_key, "admission-token-anomaly")
    stage = ("provisional" if verdict.deny_reason in PROVISIONAL_DENY_REASONS
             else "final")
    return (f" [conduit-admission/1 would_admit=0 "
            f"deny={verdict.deny_reason} stage={stage}]")


# ---------------------------------------------------------------------------
# M18.4 commit 3: standing-admit registry (design §2, POPULATION + HEALTH
# ONLY -- the revocation sweep is adjudicated amendment 3, explicitly OUT OF
# SCOPE for this arc; CONDUIT_ADMISSION_ENABLED stays False throughout, so
# every entry this registry ever holds is a WOULD-ADMIT row, never a real
# admit -- see ruling D / the round-53 "exists and reads zero" precedent).
# ---------------------------------------------------------------------------

def canonical_species_id(species):
    """design §2.4.2 (amendment 5): a species identity that survives index
    reshuffle/regeneration, unlike the run-scoped ``label(index)`` tokens
    :func:`candidate_key_from_label` builds from. Recipe, tried in order,
    NEVER raising and NEVER mutating the live ``species`` object:

      1. InChI-on-a-COPY: ``species.copy(deep=True)`` then
         ``get_augmented_inchi()`` on the copy. A deep copy isolates the
         copy's own ``.molecule``/``.aug_inchi``/resonance-structure cache
         from the live object -- ``generate_resonance_structures()`` (which
         InChI generation triggers) mutates the species it runs on, and
         this recipe must never do that to a species still live in the
         reaction being classified. Returns ``f"inchi:{aug_inchi}"``.
      2. Label-stripped, atom-canonical adjacency-list fallback: for a
         species whose InChI generation fails (e.g. resonance/valence
         edge cases some polymer proxy structures hit), take the first
         live ``Molecule``, ``.copy(deep=True)`` it (same live-mutation
         isolation), ``.clear_labeled_atoms()`` it (round-73 P2 fix:
         ``to_adjacency_list(label="", ...)`` only suppresses the
         MOLECULE header -- per-atom ``.label`` values still serialize,
         and ``Atom.copy()`` preserves them, so two copies of the SAME
         graph differing only in reaction/template atom labels would
         otherwise resolve to different, non-canonical ids), then
         ``.sort_atoms()`` it (deterministic canonical atom order for a
         given graph), then serialize with ``to_adjacency_list(label="",
         ...)`` (now genuinely label-independent, molecule header AND
         every atom). Returns ``f"adjlist:{sha256(...)[:16]}"``.
      3. Terminal sentinel: if BOTH of the above fail (or there is no
         usable ``Molecule`` at all -- e.g. an empty ``.molecule`` list),
         return a deterministic-per-failure ``f"csid-unresolved:{sha8}"``
         (the sha8 is over the failure's exception repr, so two DIFFERENT
         failure modes get visibly different sentinels rather than
         colliding on one opaque string; two IDENTICAL failures collide on
         purpose -- that is the point of a sentinel).

    NEVER raises."""
    try:
        copy = species.copy(deep=True)
        aug_inchi = copy.get_augmented_inchi()
        if aug_inchi:
            return f"inchi:{aug_inchi}"
    except Exception:
        pass
    try:
        mol_list = getattr(species, "molecule", None)
        if not mol_list:
            raise ValueError("no live Molecule to canonicalize")
        mol_copy = mol_list[0].copy(deep=True)
        mol_copy.clear_labeled_atoms()
        mol_copy.sort_atoms()
        adjlist = mol_copy.to_adjacency_list(label="", remove_h=False)
        digest = hashlib.sha256(
            adjlist.encode("utf-8", errors="backslashreplace")).hexdigest()
        return f"adjlist:{digest[:16]}"
    except Exception as exc:
        digest = hashlib.sha256(
            repr(exc).encode("utf-8", errors="backslashreplace")).hexdigest()
        return f"csid-unresolved:{digest[:8]}"


def _rxn_signature_orientation(record):
    """Read-only resolution of the G2 admitted-direction basis, DUPLICATED
    deliberately from :func:`evaluate_conduit_admission`'s own G2 block
    rather than shared with it: that evaluator is a frozen, heavily-tested
    fail-closed admission gate this arc must not touch or add a new
    caller/import edge to, and this helper is pure identity bookkeeping
    for the standing-admit registry that must never be able to influence
    an admission decision. Returns ``(direction, produced_side_entries)``
    -- ``"forward"``/``"reverse"`` and the record-side entry list for the
    credited side -- or ``(None, None)`` if no unique chain-to-pool
    direction resolves (mirrors G2's own "classifier-not-admissible"
    non-resolution). NEVER raises."""
    try:
        r_side, p_side = record["reactants"], record["products"]
        r_roles = [species_role(s) for s in r_side]
        p_roles = [species_role(s) for s in p_side]
        if (POOL in p_roles and POOL not in r_roles
                and CHAIN in r_roles and CHAIN not in p_roles):
            return "forward", p_side
        if (POOL in r_roles and POOL not in p_roles
                and CHAIN in p_roles and CHAIN not in r_roles):
            return "reverse", r_side
    except Exception:  # pragma: no cover - defensive fail-closed
        pass
    return None, None


def _canonical_reaction_signature(forward, record, result, verdict):
    """design §2.4.1 (amendment 4): the canonical REACTION signature for
    the standing-admit registry -- SIDE-SEPARATED (never side-merged the
    way the legacy :func:`candidate_key_from_label` sorts/merges sides),
    multiplicity-preserving multisets of :func:`canonical_species_id` over
    the RAW live ``forward.reactants``/``forward.products`` (so
    ``A+B->C`` and ``C->A+B`` never alias), plus the destination-pool
    identity (``verdict.dst_pool`` -- already the canonical
    index-stripped pool label G5 resolved), the gas-product identity
    (``canonical_species_id`` of the live gas species object, NOT the
    unstable ``verdict.gas_product[0]`` token), and an explicit
    ``orientation`` tuple (admitted direction, classifier shape, the G2
    irreversible-rewrite basis) so the identity captures the SAME
    direction basis G2 gates admission on -- an aligned and a
    direction-flipped sighting of otherwise-identical species never
    collide.

    Returns ``(rxn_sig_hash, ok)``. ``ok`` is False (with a deterministic
    ``rxnsig-unresolved:*`` fallback hash) if the row's orientation cannot
    be resolved here -- this should not happen for a
    ``verdict.admitted=True`` row (G2 already proved a direction resolves
    for it), but this population path never trusts that invariant
    blindly; a resolution failure is a no-op upstream (see
    :func:`register_standing_admit_would_admit`), never a corrupt insert.
    NEVER raises."""
    try:
        direction, produced_side = _rxn_signature_orientation(record)
        if direction is None:
            return "rxnsig-unresolved:no-orientation", False
        produced_objs = (
            list(getattr(forward, "products", None) or [])
            if direction == "forward"
            else list(getattr(forward, "reactants", None) or []))
        gas_idx = next((i for i, s in enumerate(produced_side)
                        if species_role(s) != POOL), None)
        gas_obj = (produced_objs[gas_idx]
                  if gas_idx is not None and gas_idx < len(produced_objs)
                  else None)
        gas_id = (canonical_species_id(gas_obj) if gas_obj is not None
                 else "gas-unresolved")
        r_objs = list(getattr(forward, "reactants", None) or [])
        p_objs = list(getattr(forward, "products", None) or [])
        r_multiset = tuple(sorted(
            Counter(canonical_species_id(s) for s in r_objs).items()))
        p_multiset = tuple(sorted(
            Counter(canonical_species_id(s) for s in p_objs).items()))
        orientation = (direction, str(result.get("shape")),
                      bool(verdict.needs_irreversible_rewrite))
        sig = (r_multiset, p_multiset, str(verdict.dst_pool), gas_id,
              orientation)
        digest = hashlib.sha256(repr(sig).encode("utf-8")).hexdigest()
        return digest[:16], True
    except Exception as exc:  # pragma: no cover - defensive fail-closed
        digest = hashlib.sha256(
            repr(exc).encode("utf-8", errors="backslashreplace")).hexdigest()
        return f"rxnsig-unresolved:{digest[:8]}", False


class _StandingAdmitRegistry:
    """Process-wide standing-admit ledger (design §2.1-§2.5). This arc
    (commit 3) implements POPULATION + HEALTH only -- adjudicated
    amendment 3 explicitly defers the revocation sweep (no entry's
    ``status`` ever becomes ``"revoked"`` here; that field/value is
    reserved for the sweep's own commit). Mirrors
    :class:`_CensusRegistry`/:class:`_EpochProvider`: one module-global
    instance, one lock, never-raise, idempotent-per-lifecycle health
    emission, reset alongside the other singletons.

    Entry, keyed by the canonical ``rxn_sig_hash`` (design §2.4 -- the
    PRIMARY key; unlike ``candidate_key`` it survives an index reshuffle
    across regeneration):
      ``rxn_sig_hash``, ``candidate_key`` (run-scoped, carried for
      cross-referencing the census/admission lines), ``status``
      (``"would-admit"`` this arc; ``"admitted"``/``"revoked"`` reserved
      for later commits), ``admit_epoch`` (write-once, the CORE EPOCH at
      first insert), ``last_seen_epoch`` (mutable, refreshed every
      re-sighting), ``snapshot`` (write-once: a plain-dict, no live-object
      copy of the direction/reversible/bucket basis the verdict was built
      on, for a future sweep to diff against), ``live_ref`` (mutable
      ``(forward, row_pools)`` tuple, refreshed every re-sighting -- the
      DEFERRED sweep's future evaluator input, populated now per ruling
      G-Q5 so it starts from real data once it lands).

    Ruling G-Q7 (write-once, double-insert protection): re-registering an
    already-known ``rxn_sig_hash`` (a readjudication re-sighting of the
    same admitted structure) refreshes ``live_ref``/``last_seen_epoch``
    ONLY; ``admit_epoch``/``snapshot`` are set exactly once, at first
    insert, and never overwritten.
    """

    def __init__(self):
        self._lock = threading.Lock()
        self._entries = {}
        self._hash_cache = {}  # memo_key -> rxn_sig_hash (run-scoped memo;
        # memo_key = (candidate_key, direction, shape,
        # needs_irreversible_rewrite) -- round-73 P1, see cached_hash()
        self._closed = False
        self._health_line = None
        self._opened = False

    def note_lifecycle_activity(self):
        """round-73 P2 fix (virgin-guard false negative): open this
        lifecycle on ANY census registration activity -- admit or deny --
        not only a would-admit insert. Previously :attr:`_opened` flipped
        only inside :meth:`register_would_admit`, so a real lifecycle that
        only ever produced DENIED rows (never a would-admit) closed
        silently -- no ``CONDUIT STANDING ADMIT HEALTH/1`` line at all,
        indistinguishable from a broken emission path, violating the
        "exists-and-reads-zero" discipline the rest of this module follows
        (see :meth:`close_lifecycle`). Called from the SAME choke point
        that marks every registration's epoch sighting
        (:meth:`_CensusRegistry.register`), so any run with census
        activity -- admit or deny -- opens this registry; a truly virgin
        process (zero registrations of any kind) still closes silently.
        Cheap: one flag set under the same lock as the rest of the
        registry. NEVER raises."""
        try:
            with self._lock:
                self._opened = True
        except Exception:  # pragma: no cover - defensive fail-closed
            pass

    def register_would_admit(self, rxn_sig_hash, candidate_key, snapshot,
                             live_ref):
        """Insert-or-refresh (design §2.3): a NEW ``rxn_sig_hash`` inserts
        with ``admit_epoch``/``snapshot`` from THIS call (write-once); an
        EXISTING one refreshes ``live_ref``/``last_seen_epoch`` only
        (ruling G-Q7). NEVER raises."""
        try:
            epoch = current_epoch()
            with self._lock:
                self._opened = True
                entry = self._entries.get(rxn_sig_hash)
                if entry is None:
                    self._entries[rxn_sig_hash] = {
                        "rxn_sig_hash": rxn_sig_hash,
                        "candidate_key": candidate_key,
                        "status": "would-admit",
                        "admit_epoch": epoch,
                        "last_seen_epoch": epoch,
                        "snapshot": snapshot,
                        "live_ref": live_ref,
                    }
                else:
                    entry["last_seen_epoch"] = epoch
                    entry["live_ref"] = live_ref
                    entry["candidate_key"] = candidate_key
        except Exception as exc:  # pragma: no cover - defensive fail-closed
            _raw_name = type(exc).__name__
            _sanitized_name = _sanitize_dynamic_token(_raw_name)
            logging.warning(
                "CONDUIT STANDING ADMIT ANOMALY/1 event=register-failed "
                "reason=%s%s action=skip-registration",
                _sanitized_name, _raw_sha_suffix(_raw_name, _sanitized_name))

    def cached_hash(self, memo_key):
        """run-scoped memo (design §2.4.2 'compute once, memoize'): a
        cache of already-computed ``rxn_sig_hash`` values, so a repeated
        re-sighting of the same live reaction across epochs
        (readjudication) skips recomputing InChI/adjacency-list identities
        entirely.

        round-73 P1 fix (hash-cache aliasing): ``memo_key`` is
        ``(candidate_key, direction, shape, needs_irreversible_rewrite)``,
        NOT bare ``candidate_key``. ``candidate_key`` deliberately omits
        arrow/reversibility (:func:`conduit_candidate_key`), but
        :func:`_canonical_reaction_signature` folds
        ``needs_irreversible_rewrite`` (and the resolved direction/shape)
        into the ``orientation`` tuple it hashes -- a bare-candidate_key
        cache would silently alias a row sighted once as ``=>`` and once
        as ``<=>`` onto the FIRST sighting's ``rxn_sig_hash``, corrupting
        the second sighting's write-once ``admit_epoch``/``snapshot``
        under the wrong identity. Widening the memo key to include every
        orientation-affecting input the signature depends on (all cheap
        to read off the verdict/record/result before the full signature
        is computed) keeps the memo correct without giving up caching.
        NEVER raises."""
        try:
            with self._lock:
                return self._hash_cache.get(memo_key)
        except Exception:  # pragma: no cover - defensive fail-closed
            return None

    def cache_hash(self, memo_key, rxn_sig_hash):
        try:
            with self._lock:
                self._hash_cache[memo_key] = rxn_sig_hash
        except Exception:  # pragma: no cover - defensive fail-closed
            pass

    def lookup(self, rxn_sig_hash):
        try:
            with self._lock:
                entry = self._entries.get(rxn_sig_hash)
                return dict(entry) if entry is not None else None
        except Exception:  # pragma: no cover - defensive fail-closed
            return None

    def close_lifecycle(self):
        """Boundary emitter (design §2.5): exactly ONE
        ``CONDUIT STANDING ADMIT HEALTH/1`` line, recomputed fresh from
        final registry state (never incrementally accumulated), emitted
        once an OPENED lifecycle closes -- present even at zero for an
        opened lifecycle (the round-53 exists-and-reads-zero witness:
        ``admitted=0`` here is a PROOF that
        :data:`CONDUIT_ADMISSION_ENABLED` is False, not a policy choice
        made by this method) -- exactly once per lifecycle, idempotent
        exactly like :meth:`_EpochProvider.close_lifecycle`.

        Virgin-lifecycle guard (design §2.5, "consistent with the epoch
        provider's ``_opened`` pattern"): if NOTHING has happened this
        lifecycle -- no ``register_would_admit`` call at all
        (:attr:`_opened` still False, e.g. a run that fails before any
        row is ever classified) -- this emits NOTHING and returns
        ``None``, a TRUE no-op that does NOT set :attr:`_closed`, so a
        later real close (once the lifecycle has since seen activity)
        still fires normally -- exactly
        :meth:`_EpochProvider.close_lifecycle`'s own virgin guard.

        ``epochs_seen=``/``burned=`` are read from the epoch provider
        (:meth:`_EpochProvider.attempted_and_burned`) so the two health
        lines report consistent epoch accounting. ``orphaned=`` and the
        real ``revoked=``/``admitted=`` counts stay reserved at 0/0/0 this
        arc (amendment 3: the sweep that would ever populate them is out
        of scope). NEVER raises."""
        try:
            with self._lock:
                if self._closed:
                    return self._health_line
                if not self._opened:
                    # Virgin lifecycle -- stay a true no-op (round-71 P2
                    # pattern, ported per design §2.5).
                    return None
                entries = list(self._entries.values())
                self._closed = True
            would_admit = sum(
                1 for e in entries if e["status"] == "would-admit")
            admitted = sum(1 for e in entries if e["status"] == "admitted")
            revoked = sum(1 for e in entries if e["status"] == "revoked")
            orphaned = 0  # sweep out of scope this arc (amendment 3)
            frozen_grams = 0.0  # pre-flip: no real admit ever contributes
            epochs_seen, burned = _EPOCH_PROVIDER.attempted_and_burned()
            line = (
                "CONDUIT STANDING ADMIT HEALTH/1 "
                f"epochs_seen={epochs_seen} burned={burned} "
                f"would_admit={would_admit} admitted={admitted} "
                f"revoked={revoked} orphaned={orphaned} "
                f"frozen_grams={frozen_grams}")
            with self._lock:
                self._health_line = line
            logging.warning(line)
            return line
        except Exception as exc:  # pragma: no cover - defensive fail-closed
            _raw_name = type(exc).__name__
            _sanitized_name = _sanitize_dynamic_token(_raw_name)
            logging.warning(
                "CONDUIT STANDING ADMIT ANOMALY/1 "
                "event=close-lifecycle-failed reason=%s%s "
                "action=skip-health-emission",
                _sanitized_name, _raw_sha_suffix(_raw_name, _sanitized_name))
            return None

    def reset(self):
        """Lifecycle reset (design §2.5, the same "reset both or neither"
        discipline every other singleton follows): close the OUTGOING
        lifecycle first (its health line is computed from pre-clear
        state), then clear all state. NEVER raises."""
        self.close_lifecycle()
        try:
            with self._lock:
                self._entries = {}
                self._hash_cache = {}
                self._closed = False
                self._health_line = None
                self._opened = False
        except Exception:  # pragma: no cover - defensive fail-closed
            pass


#: Process-wide standing-admit registry instance (design §2.1).
_STANDING_ADMITS = _StandingAdmitRegistry()


def register_standing_admit_would_admit(forward, row_pools, record, result,
                                        verdict):
    """design §2.2/§2.3, ruling D (adjudicated): populate the standing-
    admit registry from a WOULD-ADMIT verdict, pre-flip. Called from
    :func:`annotate_refused_row` right AFTER the r42 P1-4(b) post-
    registration FR re-check, so a verdict the FR-overlap check just
    downgraded to denied is never populated as a would-admit row.

    ``status`` is always ``"would-admit"`` this arc -- a real admit is
    structurally impossible while :data:`CONDUIT_ADMISSION_ENABLED` is
    False; the health line's ``admitted=0`` is therefore a PROOF (round-53
    "exists and reads zero"), not a policy decision made here.

    A no-op (never raises) when ``verdict`` is None or not admitted, or
    when the canonical reaction signature cannot be resolved -- a
    malformed/unresolvable row's identity is never trustworthy enough to
    stand as an admit entry, so it is skipped rather than inserted under a
    fallback sentinel key that could collide across unrelated rows."""
    try:
        if verdict is None or not verdict.admitted:
            return
        candidate_key = (verdict.candidate_key or result.get("candidate_key")
                        or "")
        # round-73 P1 fix (hash-cache aliasing): the memo key must include
        # every orientation-affecting input _canonical_reaction_signature
        # folds into its ``orientation`` tuple that ``candidate_key`` itself
        # omits (candidate_key deliberately drops arrow/reversibility --
        # see conduit_candidate_key). direction/shape/rewrite-basis are all
        # cheap to read off the record/result/verdict already in hand here,
        # so widening the key (option a) is preferred over dropping the
        # cache outright: it keeps the memo correct for readjudication
        # re-sightings within a single orientation while still telling a
        # `=>` sighting and a `<=>` sighting of the same candidate_key
        # apart, rather than aliasing the second onto the first's hash.
        direction, _produced_side = _rxn_signature_orientation(record)
        memo_key = (candidate_key, direction, str(result.get("shape")),
                   bool(verdict.needs_irreversible_rewrite))
        rxn_sig_hash = _STANDING_ADMITS.cached_hash(memo_key)
        if rxn_sig_hash is None:
            rxn_sig_hash, ok = _canonical_reaction_signature(
                forward, record, result, verdict)
            if not ok:
                return
            _STANDING_ADMITS.cache_hash(memo_key, rxn_sig_hash)
        snapshot = {
            "reversible": bool(record.get("reversible")),
            "shape": result.get("shape"),
            "direction": direction,
            "needs_irreversible_rewrite": bool(
                verdict.needs_irreversible_rewrite),
            "dst_pool": verdict.dst_pool,
            "gas_product_token": (verdict.gas_product[0]
                                  if verdict.gas_product else None),
            "bucket": result.get("bucket"),
        }
        _STANDING_ADMITS.register_would_admit(
            rxn_sig_hash, candidate_key, snapshot, (forward, row_pools))
    except Exception as exc:  # pragma: no cover - defensive fail-closed
        _raw_name = type(exc).__name__
        _sanitized_name = _sanitize_dynamic_token(_raw_name)
        logging.warning(
            "CONDUIT STANDING ADMIT ANOMALY/1 event=population-failed "
            "reason=%s%s action=skip-registration",
            _sanitized_name, _raw_sha_suffix(_raw_name, _sanitized_name))


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


#: round-50 FR-census scoping (ratified, adversarial review round 50):
#: ONLY these refusal reasons are genuine upstream feature-radical blockers.
#: Every other reason seen at the warn-once hook (chiefly
#: 'conduit-deferred', the accumulating=False branch -- exactly the
#: population moment_credit_conduit/1 (M18.4) will eventually try to
#: admit) is a candidate ECHOING through the same warn-once site, not a
#: genuine FR sighting, and must NOT poison the 'feature_radical' census:
#: doing so previously made a row's own refusal echo permanently G1-deny
#: itself (and any co-keyed row) as "feature-radical-overlap" -- the
#: downgrade trap this split fixes.
GENUINE_FEATURE_RADICAL_REASONS = frozenset({"qssa-invalid",
                                             "qssa-unassessable"})

#: Census name for a non-genuine (echoed) refusal sighting recorded via
#: :func:`annotate_feature_radical`. Deliberately a DIFFERENT ledger census
#: than 'feature_radical' -- G1 (:func:`evaluate_conduit_admission`) and the
#: r42 post-registration downgrade (:func:`annotate_refused_row`) key on
#: literal membership in the 'feature_radical' census set, so a row seen
#: ONLY here can never trigger FR-overlap denial of itself or of any other
#: row sharing its candidate key.
CONDUIT_ECHO_CENSUS = "conduit_echo"


def annotate_feature_radical(reaction_label, reason=None):
    """String-only annotation hook for the FEATURE-RADICAL census line
    (emitted via :func:`rmgpy.polymer._warn_once_refused`; its caller lives
    in the solver, which M18.2 must not edit, so this hook works from the
    census label alone). round-50 FR-census scoping: ``reason`` (the row's
    ``polymer_refused_reason``) decides which ledger census the sighting
    registers into -- only :data:`GENUINE_FEATURE_RADICAL_REASONS` land in
    the real 'feature_radical' census that G1 and the r42 downgrade consult;
    every other reason (esp. 'conduit-deferred') lands in
    :data:`CONDUIT_ECHO_CENSUS`, which no FR consumer reads. ``reason=None``
    (unknown/legacy caller) is treated as non-genuine -- fail-closed against
    accidentally re-widening the FR census. NEVER raises."""
    try:
        key = candidate_key_from_label(reaction_label)
        is_genuine = reason in GENUINE_FEATURE_RADICAL_REASONS
        census = "feature_radical" if is_genuine else CONDUIT_ECHO_CENSUS
        bucket = "FEATURE_RADICAL" if is_genuine else "CONDUIT_ECHO"
        entry = register_candidate(key, census, bucket)
        result = {"candidate_key": key, "shape": None, "admissible": False}
        suffix = census_suffix(result, entry)
        # r42 P1-4(a) + round-50 echo-token fix: EVERY census line carries
        # the admission verdict tokens (BUILD_SPEC W1.6 promised them on
        # every line; the FR line lacked them, and later a pure echo line
        # lacked them too). This label-only hook never has the live
        # reaction object, so it can only speak to the ONE gate it
        # genuinely evaluates here -- G1, the FR ledger membership check --
        # never to G0/G2-G7 (which need the live row and only run at the
        # r93 stamp site). When THIS key genuinely carries 'feature_radical'
        # membership (this sighting or an earlier one under a genuine
        # reason), G1 denies feature-radical-overlap -- printed here so the
        # FR census line and the ledger can never disagree. Otherwise (a
        # pure echo sighting with no feature_radical co-membership) G1 does
        # not block it, but nothing else was evaluated here either: stay
        # fail-closed (never fabricate an admit) with an honest reason that
        # says so, instead of a fabricated "feature-radical-overlap" deny a
        # pure echo never earned.
        if "feature_radical" in entry["censuses"]:
            verdict = _deny(key, "feature-radical-overlap")
        else:
            verdict = _deny(key, "echo-not-evaluated")
        suffix += admission_census_suffix(verdict)
        return suffix
    except Exception as exc:  # pragma: no cover - defensive fail-open
        logging.warning(
            "MOMENT-CREDIT CONDUIT CENSUS (M18.2): feature-radical "
            "annotation failed (%s: %s); the refusal itself is unaffected "
            "(census-only code path).", type(exc).__name__, exc)
        return ""
