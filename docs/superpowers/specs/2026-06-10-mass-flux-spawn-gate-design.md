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
4. **First sighting can never spawn — and the floor is steeper than one
   iteration.** At first arrival (iteration k) the entry has zero
   representatives → record 0, defer. The k-arrival only carries a
   simulated rate from iteration k+1 onward, so the first REAL fraction
   enters the window at the second arrival's iteration (k+1 at the
   earliest). With the sum/N statistic the exact floor is: **earliest spawn
   is the second arrival's iteration, and only if that single record is
   itself ≥ N× the bar; an exactly-at-bar motif needs real records filling
   the window, spawning at its N-th recording iteration** (k+3 for N=3 with
   arrivals every iteration; ~1.5× the bar spawns at k+2). Re-checks happen
   only at arrivals (decision 2), so arrival gaps stretch the floor
   further. Spawn-timing tests encode this arithmetic — NOT the looser
   "second sighting spawns" — so the bar is never "fixed" by loosening.
   Verified consequence scope: scission-daughter pools (the Polymer
   handshake in `make_new_reaction`, model.py:625-656) and input-file pools
   do NOT route through this gate, so the EPDM end-to-end baseline is
   expected to be a no-op (verified by test, §7.8, not assumed). Three unit
   tests that assert spawn-on-first-sighting are deliberately re-baselined
   (§7.4).

## 3. The fraction

For a species `i` attributed to pool `P`:

```
g_i(P) = max(0, core_species_production_rates[i]) * E[n]_P * monomer_MW_P

E[n]_P = y[mu1]_P / y[mu0]_P   if y[mu0]_P > SMALL_EPS, else 0.0

numerator(motif) = sum_{(i, P_i) in representatives(motif)} g_i(P_i)

denominator = sum_{canonical pool proxies j} g_j(pool(j))
            + sum_{DEDUPED representatives (i, P_i) across the ledger} g_i(P_i)

fraction(motif) = numerator(motif) / denominator
```

**Attribution rules (AMENDED 2026-06-10 after plan-author probe — see §3.1):**

- **Canonical pool proxies** (exactly one per configured pool,
  `is_pool_proxy`, label-matched at polymer.pyx:478-485) attribute to their
  own pool via the engine's `species_to_pool_indices`.
- **Representatives are ORDINARY absorbed species** — the engine has no
  pool mapping for them. Each representative is recorded in the ledger as a
  `(species_label, parent_pool_label)` pair, where the parent pool is the
  pool that absorbed the candidate (for Phase-D deferred candidates: the
  pool that would parent its spawn intent — currently `pool_registry[0]`;
  if multi-parent attribution ever lands, this follows it). E[n] and
  monomer_MW are read LIVE from that pool's stats in the same snapshot
  (live freshness kept; the wrong-mapping hole closed).
- **Multi-motif double-counting is a stated decision, not emergent:** E[n]
  is the absorbing pool's chain statistics, used as the event-mass
  calibration for THAT motif's attribution; a species absorbed into pool P
  contributes its gross production to EVERY motif-entry that lists it as a
  representative, with P's E[n]. A species legitimately carrying two motifs
  (decoration + backbone) double-counts across the two gate evaluations —
  accepted, because the two motifs are competing for DIFFERENT pool slots
  and each should see the mass.
- **The denominator dedups:** a species appearing in multiple motif entries
  is counted ONCE in the denominator. Including all canonical proxies plus
  all deduped ledger representatives makes each motif's numerator a subset
  of denominator terms, so each fraction ∈ [0,1] (the pinned property);
  the SUM of fractions across motifs may exceed 1 — accepted, per the
  double-counting decision above.

### 3.1 Why representatives need change (a) — the born-dead hole

The original spec sourced representative flux via `species_to_pool_indices`
and the proxy-only gross arrays. Plan-author probing proved that is
identically zero in production: `is_pool_proxy` marks only the ONE
canonical proxy per pool, and the gross arrays are maintained ONLY inside
the `is_pool_proxy` branches (polymer.pyx:1082-1117) — absorbed
representatives are ordinary species with `dn_dt` writes only. A gate built
on that returns zero for every motif under all production conditions: not
conservative — broken, with green tests (the §7 fabricated-snapshot suite
cannot see it). **Change (a):** the polymer residual maintains gross
production/consumption for ALL core species (mirroring `simple.pyx`'s
per-species bookkeeping), so representative production — which is real;
unlike canonical proxies their flux is NOT rerouted to moments — has a
production-side record. Cost-gated per §4.6.

- **Gross, not net.** Both sums read `core_species_production_rates`.
  For canonical proxies the gross arrays exist precisely so diagnostics
  survive the moment rerouting (polymer.pyx:1082-1117); for ordinary
  species (representatives) they are maintained by change (a) (§3.1/§4.6).
  Never `dn_dt`-derived net rates (see decision 3) — net is rerouted-to-≈0
  for proxies and consumption-cancelled at steady state for ordinary
  species.
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
- **Denominator = total polymer-derived gross event-mass** (canonical
  proxies + deduped ledger representatives, per the formula block above).
  Immune to the two failures of a net-dµ1 form:
  transfer chemistry double-counting (parent loss + daughter gain ≈ 2×) and
  SAME_POOL chemistry netting to zero (a deck doing furious modification
  chemistry with slow net degradation would otherwise have a near-zero
  denominator and wave trace motifs through). Each motif's fraction ∈ [0,1]
  by construction: "this motif's share of polymer-derived gross
  production."
- **Units.** `core_species_production_rates` entries are volumetric
  (mol/m³-phase/s); every term in both sums is polymer-phase, so V_poly
  cancels. No volume plumbing.

## 4. Components and data flow

### 4.1 Engine snapshot (`rmgpy/solver/polymer.pyx` — `make` before tests)

New method on `HybridPolymerSystem`, `spawn_gate_flux_snapshot()` —
AMENDED: the engine cannot attribute representatives (it has no ledger), so
the engine/python split is: the engine reads arrays, the gate does
ledger-dependent attribution. Returns a 3-tuple:

1. `gross`: dict core-species label → `max(0, core_species_production_rates[i])`
   for ALL core species (label via the species↔index mapping the engine
   already holds; requires change (a), §4.6, for ordinary species to be
   nonzero).
2. `pool_stats`: dict pool label → `(E_n, monomer_mw_g_mol)` with the
   SMALL_EPS guard applied to E_n (per-pool, from
   `pool_mu0_indices`/`pool_mu1_indices` and live `y`).
3. `proxy_event_mass_total`: float — `Σ g_j(pool(j))` over canonical
   proxies (the engine CAN attribute those, via `species_to_pool_indices`).

The gate (python side) computes representative `g_i(P_i)` as
`gross[label] · E_n[parent_pool] · mw[parent_pool]`, the numerator per
motif, and `denominator = proxy_event_mass_total + Σ deduped representative
g_i`. Pool monomer MW reaches the engine via a new
`PolymerPoolConfig.monomer_mw_g_mol` field (default 0.0 → g=0 → defer,
honest degradation), plumbed in `polymer_input.py::to_config` with the same
`get_molecular_weight()*1000` idiom as `Polymer.monomer_mw_g_mol`
(polymer.py:261). Beyond change (a), no new bookkeeping inside the
residual.

### 4.2 main.py stash

Immediately after each `simulate()`, main.py computes the snapshot on the
ENGINE — `system.solver`, not the `HybridPolymerReactor` blueprint (the
established blueprint-vs-engine gotcha) — and stashes it on the reaction
model as `reaction_model.polymer_flux_snapshot`. Skipped (left `None`) for
non-polymer systems. This stash + the ledger are exactly the shared
infrastructure the §2.2 upgrade would reuse.

### 4.3 Motif ledger (attribute on `CoreEdgeReactionModel`)

`reaction_model.polymer_motif_ledger`: a list of entries
`{motif: Group, representatives: [(species_label, parent_pool_label)],
accumulator records, last_recorded_iteration: int, spawned: bool}` — the
parent-pool label recorded at absorption per §3's attribution rules.

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
   `fraction(motif)` per §3 — numerator from `gross` × the recorded
   parent-pool's `pool_stats`, denominator =
   `proxy_event_mass_total` + deduped representative event-mass across the
   ledger — and `record()` it; otherwise do NOT record — same-iteration
   burst arrivals re-check the gate against the existing window only.
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
6. Either way, append `(candidate_label, parent_pool_label)` to
   `entry.representatives` (the candidate is absorbed as an ordinary
   explicit species; with change (a) its gross production is recorded from
   the next simulation onward).

The RMG iteration counter is plumbed alongside the snapshot stash (the
accumulator's `record(..., iteration=...)` argument already exists).

### 4.5 Honest degradation

No snapshot stashed (`polymer_flux_snapshot is None`: iteration 0, or
`reaction_model is None` in the unit-test path) → fraction 0.0 → defer,
logged. The `0.5` constant and the `reaction_model is None → 0.5` special
case are DELETED. Tests that exercise Phase E supply a fabricated snapshot
and a pre-populated ledger; no production code path fakes a number.

### 4.6 Change (a): gross arrays for all core species — cost-gated

The polymer residual's sections 3/4 (reactant/product flux application)
gain gross production/consumption writes for ordinary (non-proxy) core
species, mirroring `simple.pyx`'s per-species bookkeeping (reactants:
`consumption += rf, production += rr`; products: `production += rf,
consumption += rr`; core reactions only, same gate as the existing
branches). This is the first change since the apportionment work that adds
unconditional per-species writes back into the residual hot loop, so it is
COST-GATED, decided on measured data, not assumption:

- **Measure before committing:** EPDM deck wall-clock, before vs after
  (same machine, same session, ≥2 runs each). The simple.pyx precedent
  says this is almost certainly fine — but simple.pyx's residual does not
  also carry the moment dispatch.
- **Acceptance is a DUAL gate:** (i) slowdown within run-to-run noise
  (≤ ~5% wall-clock), AND (ii) a CONSERVATION gate — change (a) adds
  diagnostic writes only, so the EPDM deck's final pool moments (and
  Σµ1·MW across pools) must be IDENTICAL (rtol 1e-12) between the before
  and after timing runs. Speed was never the risk on this change; silent
  mass drift is — a conservation regression in the hot loop would surface
  as slow drift the per-task functional tests cannot catch.
- **Fallback if the wall-clock bites:** maintain gross writes only for
  LEDGER-TRACKED species (a per-species flag array set from the ledger at
  `initialize_model`) — narrower, uglier, bounded. The spec prefers the
  unconditional form for simplicity and simple.pyx parity. The
  conservation gate applies to the fallback too.

## 5. `triggering_moles` — named-consumer limitation (in scope: the TODO; out of scope: the fix)

`SpawnIntent.triggering_moles` is populated from
`getattr(cand, "amount", 1.0)` — `amount` is never assigned anywhere, so the
value is always the placeholder 1.0. Its consumer is
`drain_spawn_intents` (`rmgpy/polymer.py:3123`): it seeds the daughter
pool's initial moments (`mu0 = N`, `mu1 = N·DP`, `mu2 = N·DP²`). Every
mid-run GATE-PATH spawned pool is therefore currently seeded with
mu0 = 1.0 mol of chains — a physically enormous fiction, the same defect
class this change exists to kill, and larger in physical consequence than
the 0.5 gate (a wrongly-spawned pool wastes compute; a wrongly-seeded pool
injects fictional mass into the conservation accounting).

**Scope of the defect (verified):** Path-B (gate-path) spawns only. Path-C
scission daughters are constructed with `moments=None, initial_mass=0.0`,
registered via `_register_polymer` (which seeds no moments), and never
receive engine pool config — their reactions demote to UNRESOLVED. The
EPDM baseline therefore carries no mu0-fiction (§7.8 verifies correctness,
not bug-for-bug).

**Sequencing (explicit):** the gate shipped here is correct in isolation;
a gate-path spawned pool's STATE remains known-fictional pending the
seeding item. Nobody should trust mid-run gate-path pool masses until both
land. The running-log item is filed as the next physics item — not a
someday-TODO — with candidate source: the triggering reaction's snapshot
flux × a dt-scale (needs its own chemistry decision).

Removing the field is not possible (real consumer). This change: (a)
replaces the bare `getattr` with an explicit named-consumer TODO comment at
the construction site stating "placeholder 1.0 mol seeds daughter moments
in drain_spawn_intents — see <running-log item>"; (b) records the item in
the polymer running log as above. No behavior change to seeding here.

## 6. Doc amendment (`docs/multi_pool_design.md`, same change)

§4.4 is rewritten to state: the arrival-driven rule (decision 2 wording),
the gross-production E[n]-calibrated fraction (§3 formula verbatim incl.
the ε-clamp direction sentence), the one-record-per-iteration rule and
sum/N window statistic, the spawn-floor arithmetic (decision 4 verbatim:
second arrival's iteration only at ≥ N× the bar, N-th recording iteration
at the bar), honest degradation (§4.5), the restart consequence (§4.3),
and the §2.2 upgrade trigger verbatim. The §8 "NOT YET ACTIVE" limitation note is removed.

## 7. Tests

**Standing convention (the deliverable of the 2026-06-10 review pass,
recorded in the running log):** every laundered-quantity fix ships a test
that FAILS in the live path until the quantity is real. Fabricated-input
tests prove formula correctness; they cannot prove the live path sources
the quantity — twice now (the 0.5 gate, the born-dead representative flux)
a green fabricated-input suite hid an identically-zero production value.
This convention applies to the remaining spawn/seeding items (#14
`triggering_moles` seeding) without per-spec rediscovery.

1. **Integrated tripwire (live path, two halves, RED-FIRST).** An
   integrated solve on a real `HybridPolymerSystem` — NOT a fabricated
   snapshot. Fixture: a pool with apportioned (non-UNRESOLVED) reactions
   where at least one product is an ordinary (non-canonical-proxy) core
   species standing in for an absorbed representative (representative
   status is a python/ledger concept; to the solver it is any ordinary
   species produced by pool-touching chemistry — no multi-iteration RMG run
   needed). After the solve:
   - **Numerator half (the regression that would have caught this class):**
     that ordinary species has a NONZERO `gross` entry in
     `spawn_gate_flux_snapshot()`, equal to
     `max(0, core_species_production_rates[i])` recomputed independently
     from the engine arrays, and its `g_i` under its parent pool's
     `pool_stats` equals `gross · E[n] · monomer_MW`. **This half MUST be
     written and confirmed RED against current HEAD before change (a) is
     implemented** — ordinary species have no gross writes today. That red
     run is the proof the fix is real, and a gate to executing the rest of
     the plan. **Pin the reason, not just the redness:** a LIVENESS
     assertion (`dn_dt[i] > 0` for the same species) must precede the
     gross assertion inside the test body, so the red can only mean
     "chemically alive but no gross record" — never "fixture dead".
     xfail(strict) cannot distinguish those; the inner assertion can. The
     red traceback is inspected BY HAND (controller, not implementer
     say-so) before any production code lands.
   - **Denominator half:** the canonical proxy's net `dn_dt` contribution
     is ≈ 0 (the apportionment reroutes proxy flux to pool moments) while
     its gross entry is nonzero — the assertion that is only true of the
     gross array and dies if the denominator path is ever rewired to net
     rates.
   The polymer.pyx line citations in §3 are informational; the test asserts
   the mechanism, not the address.
2. **Spawn-timing per the §2.4 floor:** novel motif, first arrival defers;
   fabricated snapshots; assert BOTH branches of the floor arithmetic — a
   ≥ N×-bar motif spawns at its second arrival's iteration, an
   exactly-at-bar motif defers until its N-th recording iteration. (Not the
   looser "second sighting spawns".)
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
   species/reaction counts, and sidecar content unchanged. Its only
   daughter pool is the Path-C scission handshake, which bypasses this gate
   AND the `triggering_moles` seeding (verified, §5) — so this asserts
   correctness, not a bug-for-bug match.
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

- `rmgpy/solver/polymer.pyx` — `spawn_gate_flux_snapshot()`; change (a)
  gross writes for ordinary species in residual sections 3/4 (cost-gated,
  §4.6); `PolymerPoolConfig.monomer_mw_g_mol` field (run `make`).
- `rmgpy/solver/polymer_input.py` (or wherever `to_config` lives — plan
  pins the exact path) — monomer MW plumb into `PolymerPoolConfig`.
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
