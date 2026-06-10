# Mass-Flux Spawn Gate (multi-pool §4.4) — Design

**Date:** 2026-06-10
**Branch:** `polymer`
**Status:** Approved design (brainstormed with Alon; supersedes the §4.4
sketch in `docs/multi_pool_design.md`, which this change amends).

## 1. Problem

`_estimate_relative_flux` (`rmgpy/polymer.py:2399`) returns `0.5`
unconditionally. With the default `mass_flux_threshold=0.01`, the Phase-D
mass-flux spawn gate in `process_polymer_candidates_multipool` never blocks
anything; pool spawning is gated only by similarity-merge and the
`max_pools` cap. `MassFluxAccumulator` is implemented and unit-tested but
never instantiated in the live path. One test
(`test_mass_flux_below_threshold_defers_spawn`) "passes" only by setting the
threshold above the hardcoded constant — a fake gate with a test laundering
it. This change replaces the stub with a real, honest gate.

## 2. Decisions (with rationale, in dialogue order)

1. **Purpose: noise suppression, not slot prioritization.** The gate is a
   minimum-significance bar keeping trivial/transient motifs from spawning
   pools (a spurious pool's real cost is the reaction-generation cascade of
   its proxy variants, not its 3 ODEs). It is explicitly NOT a ranking
   mechanism: a threshold cannot prioritize slots, and an important motif
   arriving after the cap fills is a different feature (ranking/eviction),
   out of scope. That feature becomes a design item if and when a real deck
   hits `max_pools` with a significant motif locked out.
2. **Re-evaluation is arrival-driven.** Deferred candidates are never
   re-presented (`new_species_list` contains only new species), so the gate
   re-checks a motif only when a NEW candidate carrying it arrives.
   Rationale: arrivals are themselves the relevant signal — a channel that
   matters keeps generating distinct motif-carrying products — and a
   deferred motif's mass still flows through the parent pool's accounting
   (absorbed as proxy variants), so the cost of the known blind spot (flux
   grows while arrivals are quiet) is statistical distinction, not mass
   conservation. The iteration-boundary ledger re-check (full §4.4
   semantics) is the documented upgrade; **upgrade trigger:** a real deck
   shows a deferred motif's flux growing while its arrivals are quiet.
3. **Recorded quantity: snapshot-attributed gross mass-flux fraction**
   (formula in §3). Rejected alternatives: pure arrival-count bar (arrivals
   decide WHEN to check; the gate's question — "would distinct accounting
   change the chemistry?" — is mass-driven), and net species rates (proxy
   variants have `dn_dt ≈ 0` BY DESIGN, the apportionment routes their flux
   to pool moments; a net-rate gate would be born dead and silently never
   spawn — see the tripwire test, §7.1).
4. **First sighting can never spawn.** Zero prior representatives → zero
   numerator → defer. Earliest spawn is the second arrival, ≥1 iteration
   later. Verified consequence scope: scission-daughter pools (the Polymer
   handshake in `make_new_reaction`, model.py:625-656) and input-file pools
   do NOT route through this gate, so the EPDM end-to-end baseline is
   expected to be a no-op (verified by test, §7.8, not assumed). Three unit
   tests that assert spawn-on-first-sighting are deliberately re-baselined
   (§7.4).

## 3. The fraction

For proxy species `i` mapped to pool `p(i)` (via the engine's
`species_to_pool_indices`):

```
g_i = max(0, core_species_production_rates[i]) * E[n]_p(i) * monomer_MW_p(i)

E[n]_p = y[mu1]_p / y[mu0]_p   if y[mu0]_p > SMALL_EPS, else 0.0

fraction(motif) = sum_{i in representatives(motif)} g_i
                  / sum_{all proxy i} g_i
```

- **Gross, not net.** The numerator reads
  `core_species_production_rates` — the gross production array maintained
  for proxy species at `rmgpy/solver/polymer.pyx:1082-1117` precisely so
  diagnostics survive the moment rerouting. Never `dn_dt`-derived net rates
  (see decision 3).
- **E[n] calibration.** One mole of representative production is one mole of
  *events*; the mass entering the motif class per event is a chain's worth
  (~E[n]·monomer_MW), not the proxy fragment's ~3·monomer_MW. Using proxy MW
  would miscalibrate the threshold by E[n]/3 (a 30%-of-mass channel at
  E[n]=60 would read 1.5%) and drift with pool statistics. E[n] is read LIVE
  from the engine state at snapshot time (never recorded-and-stale).
- **ε-guard errs toward deferral.** When `mu0 <= SMALL_EPS`, E[n] is 0 and
  `g_i = 0`: the motif defers. The inflating direction (mu1 large while
  mu0 ≈ 0) is effectively unreachable inside the realizability cone. The
  exhaustion test (§7.7) asserts the deferral direction, not just
  finiteness. `SMALL_EPS` is the EXISTING solver constant
  (`rmgpy/solver/polymer.pyx:70`), and the guard is the exact idiom already
  used for `tail_mean` at polymer.pyx:1343 — no new constant.
- **Denominator = total polymer-derived gross event-mass** (same form summed
  over ALL proxy species). Immune to the two failures of a net-dµ1 form:
  transfer chemistry double-counting (parent loss + daughter gain ≈ 2×) and
  SAME_POOL chemistry netting to zero (a deck doing furious modification
  chemistry with slow net degradation would otherwise have a near-zero
  denominator and wave trace motifs through). Fraction ∈ [0,1] by
  construction: "this motif's share of polymer-derived gross production."
- **Units.** `core_species_production_rates` entries are volumetric
  (mol/m³-phase/s); every term in both sums is polymer-phase, so V_poly
  cancels. No volume plumbing.

## 4. Components and data flow

### 4.1 Engine snapshot (`rmgpy/solver/polymer.pyx` — `make` before tests)

New method on `HybridPolymerSystem`, `spawn_gate_flux_snapshot()`:
returns `(per_species, total)` where `per_species` maps core-species label →
`g_i` for every species with `species_to_pool_indices[i] != -1` and
`is_pool_proxy[i]`, computed from `core_species_production_rates`,
`pool_mu0_indices`/`pool_mu1_indices`, live `y`, and each pool's
`monomer_MW`; `total = sum(per_species.values())`. Pure read of
already-maintained state; no new bookkeeping inside the residual.

### 4.2 main.py stash

Immediately after each `simulate()`, main.py computes the snapshot on the
ENGINE — `system.solver`, not the `HybridPolymerReactor` blueprint (the
established blueprint-vs-engine gotcha) — and stashes it on the reaction
model as `reaction_model.polymer_flux_snapshot`. Skipped (left `None`) for
non-polymer systems. This stash + the ledger are exactly the shared
infrastructure the §2.2 upgrade would reuse.

### 4.3 Motif ledger (attribute on `CoreEdgeReactionModel`)

`reaction_model.polymer_motif_ledger`: a list of entries
`{motif: Group, representatives: [species labels], accumulator records,
last_recorded_iteration: int, spawned: bool}`.

- **Lookup is by Group isomorphism** (the same matching `similarity_merge`
  uses), not a canonical string key — sidesteps Group canonicalization; the
  ledger is `max_pools`-scale, O(n) isomorphism lookup is fine.
  `MassFluxAccumulator` keeps its API; its keys are per-entry opaque ids
  assigned at first sighting.
- **Lifecycle:** entries are created at a motif's first arrival; when a
  motif spawns, its entry is marked `spawned=True` and later arrivals of
  that motif never re-run the gate (they classify against the new pool in
  Phase A anyway; the flag removes the dead re-evaluation path and makes
  the lifecycle explicit).
- **Not serialized — restart consequence (correct-but-loud):** the ledger is
  in-memory state on the reaction model; an RMG restart resets windows and
  deferred motifs re-earn their bar. Same philosophy as unstamped-reaction
  demotion: graceful, conservative, logged (the gate logs at INFO when a
  motif is deferred, including fraction, bar, and window occupancy).

### 4.4 The gate (`_estimate_relative_flux` replaced)

At a candidate's arrival in Phase D, with the candidate's discovered motif:

1. Ledger lookup (isomorphism). Miss → create entry.
2. If `entry.spawned`: skip the gate (should be unreachable — Phase A
   matches the spawned pool first; assert-log if hit).
3. **One record per motif per RMG iteration:** if
   `entry.last_recorded_iteration < current_iteration`, compute
   `fraction(motif)` from the stashed snapshot over `entry.representatives`
   and `record()` it; otherwise do NOT record — same-iteration burst
   arrivals re-check the gate against the existing window only.
   Representatives absent from the snapshot (absorbed this iteration, not
   yet simulated) contribute 0 to the numerator — stated, not incidental. (Multiple
   same-motif candidates routinely land in one `enlarge()` pass and would
   otherwise fill the window with triplicates of ONE simulation's
   measurement, silently converting "windowed over N iterations" into
   "windowed over N arrivals".)
4. Gate statistic: **`flux(key) / window`** — sum of in-window records
   divided by the FIXED window length N (zero-filled semantics). A
   single-snapshot spike must be N× the bar to clear it; a channel
   persisting at fraction f for N iterations reads f.
5. Spawn if statistic ≥ `mass_flux_threshold`; on spawn, mark
   `spawned=True` and stamp `SpawnIntent.mass_flux_at_spawn` with the real
   statistic.
6. Either way, append the arriving candidate's label to
   `entry.representatives` (it is absorbed as a proxy variant and will
   carry rates in the next snapshot).

The RMG iteration counter is plumbed alongside the snapshot stash (the
accumulator's `record(..., iteration=...)` argument already exists).

### 4.5 Honest degradation

No snapshot stashed (`polymer_flux_snapshot is None`: iteration 0, or
`reaction_model is None` in the unit-test path) → fraction 0.0 → defer,
logged. The `0.5` constant and the `reaction_model is None → 0.5` special
case are DELETED. Tests that exercise Phase E supply a fabricated snapshot
and a pre-populated ledger; no production code path fakes a number.

## 5. `triggering_moles` — named-consumer limitation (in scope: the TODO; out of scope: the fix)

`SpawnIntent.triggering_moles` is populated from
`getattr(cand, "amount", 1.0)` — `amount` is never assigned anywhere, so the
value is always the placeholder 1.0. Its consumer is
`drain_spawn_intents` (`rmgpy/polymer.py:3123`): it seeds the daughter
pool's initial moments (`mu0 = N`, `mu1 = N·DP`, `mu2 = N·DP²`). Every
mid-run spawned pool is therefore currently seeded with mu0 = 1.0 mol of
chains — a physically enormous fiction. Removing the field is not possible
(it has a real consumer); fixing the seeding is a separate design item
(candidate source: the triggering reaction's snapshot flux × dt-scale; needs
its own chemistry decision). This change: (a) replaces the bare `getattr`
with an explicit named-consumer TODO comment at the construction site
stating "placeholder 1.0 mol seeds daughter moments in
drain_spawn_intents — see <running-log item>"; (b) records the item in the
polymer running log. No behavior change to seeding here.

## 6. Doc amendment (`docs/multi_pool_design.md`, same change)

§4.4 is rewritten to state: the arrival-driven rule (decision 2 wording),
the gross-production E[n]-calibrated fraction (§3 formula verbatim incl.
the ε-clamp direction sentence), the one-record-per-iteration rule and
sum/N window statistic, the first-sighting-defers consequence, honest
degradation (§4.5), the restart consequence (§4.3), and the §2.2 upgrade
trigger verbatim. The §8 "NOT YET ACTIVE" limitation note is removed.

## 7. Tests

1. **Tripwire (gross-array pin):** engine-integrated test — after a solve
   with apportioned (non-UNRESOLVED) reactions, a proxy-variant
   representative must yield a NONZERO numerator via
   `spawn_gate_flux_snapshot()`. Dies if the numerator is ever rewired to
   net rates (whose proxy entries are ≈0 by design).
2. **Second-sighting spawns:** novel motif, first arrival defers; fabricated
   snapshot giving above-bar fraction; arrivals on later iterations spawn
   once sum/N clears the bar.
3. **Trace motif stays deferred** across many arrivals (below-bar fraction).
4. **Re-baselined first-sighting tests:** `test_novel_motif_spawns_one_pool`,
   `test_iteration_boundary_pass_spawns_and_registers`,
   `test_novel_product_spawns_and_registers_and_serializes` become
   second-sighting tests (pre-populated ledger + snapshot).
   `test_mass_flux_below_threshold_defers_spawn` (the gamed test) is
   DELETED, replaced by #2/#3.
5. **Same-iteration burst:** three same-iteration arrivals of one motif at
   an above-bar fraction must still defer (window holds ONE entry → sum/N
   below bar), and the ledger must hold all three representatives.
6. **E[n] calibration:** a migration-shaped representative with parent-pool
   E[n]=60 must contribute ~60·monomer_MW per event-mole (≈20× the
   proxy-MW accounting) — pins decision 3's calibration.
7. **ε-clamp direction:** mu0 < ε with tiny mu1 → fraction underestimates →
   defers (asserts direction, not just finiteness).
8. **EPDM end-to-end no-op:** rerun the EPDM fixture; spawn behavior,
   species/reaction counts, and sidecar content unchanged (its only daughter
   pool is the Path-C scission handshake, which bypasses this gate).
9. **Ledger isomorphism lookup:** same motif Group with permuted atom
   ordering must hit one ledger entry (pins the isomorphism-not-string-key
   choice).
10. **Spawned-flag lifecycle:** after spawn, a later same-motif arrival does
    not re-run the gate (assert via the entry's record count).
11. **Window statistic units:** accumulator sum/N math, zero-filled
    semantics, eviction at window edge (extends the existing accumulator
    tests).

## 8. Scope boundary

No iteration-boundary re-check hook (documented upgrade only), no pool
ranking/eviction, no change to Path C (scission handshake) or Path A (input
pools), no change to `max_pools`, no change to `discover_repeat_motif`
(item #9), no daughter-moment seeding fix (§5 item), no artifact/sidecar
schema change (`mass_flux_at_spawn` already exists in spawn metadata and
simply becomes real). Defaults retained: `mass_flux_threshold=0.01`,
`window=3` — both now mean what they say.

## 9. Affected files

- `rmgpy/solver/polymer.pyx` — `spawn_gate_flux_snapshot()` (run `make`).
- `rmgpy/rmg/main.py` — snapshot + iteration stash after `simulate()`.
- `rmgpy/rmg/model.py` — ledger attribute init; pass iteration/snapshot into
  the spawn pass.
- `rmgpy/polymer.py` — real `_estimate_relative_flux` (or its replacement;
  the gate logic per §4.4), ledger helpers, `MassFluxAccumulator` wiring,
  `triggering_moles` TODO.
- `docs/multi_pool_design.md` — §4.4 rewrite, §8 note removal.
- `test/rmgpy/polymerMultiPoolTest.py` (+ a solver-side test file for the
  tripwire) — §7.
- Running log memory — `triggering_moles` seeding item.
