# Dynamic Multi-Pool Polymer Spawning — Design Doc

**Branch:** `polymer`
**Status:** Design complete, implementation in progress
**Companion doc:** `~/Code/TA/PRD.md` (consumer of the output sidecar)
**Last updated:** 2026-06-07 (reconciled with implementation: scission/unzip moment ODEs specified in §5; the §7 in-place CVODES-resize design was NOT built — spawning uses per-iteration solver rebuild; §10 limitations expanded)

---

## 1. Problem

The current polymer feature on this branch supports exactly **one** [μ₀, μ₁, μ₂] moment vector per `polymer()` declaration. For a phenolic-resin pyrolysis (the carbon-phenol use case), the network produces structurally distinct chain populations during the run — a parent novolak chain and a char-precursor chain (cross-linked, pendant-CH₂ stripped, aromatic-fused). One moment vector cannot represent both. Without a second pool, the model collapses everything that is structurally "polymer" into the parent's distribution, which loses physical meaning.

We need **dynamic multi-pool spawning**: detect the emergence of a structurally distinct chain population during the run and spawn a new pool with its own [μ₀, μ₁, μ₂].

## 2. Resolved design decisions

| # | Decision | Resolution |
|---|---|---|
| 1 | Spawn trigger | Auto-detect from product structure (no user pre-declaration) |
| 2 | Detection algorithm | Hybrid: classify-against-all-pools first; novel-monomer discovery fallback |
| 3 | Subgraph isomorphism backend | RMG-native `find_subgraph_isomorphisms`, fingerprint-gated by element formula. RDKit only as documented fallback if benchmarks show RMG's pure-Python iso is intractable for the auto-detect inner loop. |
| 4 | Spawn timing | Between RMG iterations only. Queue spawn intents during iteration; drain → register daughter `_muN` species; solver is rebuilt on the next iteration (no in-place reinit — see §7) |
| 5 | Conservation | Event-spawn instantiation from triggering product (B.μₖ = N · DPᵏ); ongoing transfer reactions use Schulz-Flory closure |
| 6 | Runaway guard | (a) similarity-merge to existing pool first, (b) mass-flux gate ≥1% over N iters, (c) `max_pools` cap (default 5, configurable) |
| 7 | Daughter end_groups | Inherit parent's. Char-precursor edge differences are a known phase-1 limitation. |
| 8 | Output contract | Pseudo-species `<pool_label>_mu0/1/2` in chemkin and cantera + sidecar `polymer_pools.json` (single source of truth for pool metadata) |
| 9 | Acceptance | Synthetic unit test + carbon-phenol integration test |

These match the PRD §6.2.

## 3. Code surface

| File | Role | Change scope |
|---|---|---|
| `rmgpy/polymer.py` | `Polymer` class, `classify_structure()`, `_analyze_wing_matches()`, `process_polymer_candidates()` | Extend: classify-against-all-pools loop; add novel-monomer discovery; spawn-intent generation |
| `rmgpy/rmg/polymer_input.py` | `polymer()`, `polymer_phase()`, `hybridPolymerReactor()` input directives | Add optional `max_pools=5` and `mass_flux_threshold=0.01` knobs to `polymer()`. Pool registry (currently 1) becomes a list. |
| `rmgpy/rmg/model.py` | Core-model iteration; `is_polymer_proxy` propagation | Drain spawn-intent queue between iterations; instantiate new `Polymer` objects; register with reaction_model |
| `rmgpy/solver/polymer.pyx` | `HybridPolymerSystem`; state-vector layout; `mu_indices` per pool | Per-pool moment ODEs (scission/unzip, §5) + closure helpers. Moment indices resolved by label in `initialize_model`; new pools picked up on solver rebuild (no in-place resize). |
| `rmgpy/data/kinetics/family.py` | Reaction generation; polymer-proxy contagion | Pass *all* live pools to `process_polymer_candidates` (not just one) so contagion classifies against the full pool registry |
| `rmgpy/yml.py` (Cantera writer) | Cantera YAML emission | Loop over pool registry to emit `<label>_muN` pseudo-species; emit sidecar `polymer_pools.json` |
| `rmgpy/chemkin.py` | Chemkin writer | Same pool loop for chem.inp emission |
| `test/rmgpy/polymerTest.py` | Polymer unit tests | New `TestMultiPoolSpawn` class (synthetic) |
| `test/rmgpy/solver/solverPolymerTest.py` | Solver tests | Carbon-phenol integration test |
| `examples/rmg/polystyrene/input.py` | Reference example | Verify still works with multi-pool code path (single-pool case ≡ multi-pool with N=1) |

## 4. Algorithm

### 4.1 Spawn-intent generation (in `process_polymer_candidates`)

```
INPUT: candidates: list[Species],       # newly-generated polymer-proxy products
       pool_registry: list[Polymer],    # all currently live pools
       reaction_model

OUTPUT: classified_candidates,          # mapping product -> (pool, polymer_class)
        spawn_intents: list[SpawnIntent]

1. spawn_intents = []
2. for product in candidates:
3.     # Phase A: classify against every existing pool, pick best non-GAS/UNKNOWN
4.     best = None
5.     for pool in pool_registry:
6.         klass, details = classify_structure(product, pool)
7.         if klass not in (GAS, UNKNOWN):
8.             best = (pool, klass, details)
9.             break  # first non-trivial classification wins; matches existing single-pool semantics
10.    if best:
11.        classified_candidates[product] = best
12.        record_mass_flux(best[0], product)   # for future spawn gating
13.        continue
14.    # Phase B: novel-monomer discovery
15.    motif = discover_repeat_motif(product)   # see 4.2
16.    if motif is None:
17.        classify as GAS                       # not a polymer chain after all
18.        continue
19.    # Phase C: similarity-merge — does this motif belong to an existing pool?
20.    merged_pool = similarity_merge(motif, pool_registry)
21.    if merged_pool:
22.        merged_pool.feature_monomers.append(motif)   # extend pool
23.        classified_candidates[product] = (merged_pool, FEATURE, {...})
24.        continue
25.    # Phase D: would-be-novel-pool — gate on mass-flux + cap
26.    if not mass_flux_gate(motif, threshold=0.01):
27.        defer (do not spawn yet; tag product as polymer-proxy without pool assignment)
28.        continue
29.    if len(pool_registry) >= max_pools:
30.        log warning; defer
31.        continue
32.    # Phase E: queue spawn intent
33.    spawn_intents.append(SpawnIntent(
34.        parent_pool=best_parent_for_motif(motif, candidates),
35.        monomer=motif,
36.        end_groups=parent.end_groups,        # phase-1 inheritance
37.        triggering_product=product,
38.        triggering_dp=count_motif_in(product, motif),
39.        triggering_moles=stoichiometric_amount_of(product),
40.    ))
41. return classified_candidates, spawn_intents
```

### 4.2 `discover_repeat_motif(product)` — novel-monomer discovery

```
INPUT:  product: Molecule
OUTPUT: motif: Group | None

Approach: enumerate connected sub-graphs of product, find the largest
that occurs ≥2 times via subgraph-isomorphism, return its Group.

1. if product.atoms < 2 * polymer_branch_min_motif_atoms:    # = 4 by default
2.     return None
3. candidates = enumerate_connected_subgraphs(product, min_atoms=4, max_atoms=product.atoms // 2)
4. # Order candidates largest-first to prefer chemically-meaningful (longer) motifs
5. candidates.sort(key=size, reverse=True)
6. for sub in candidates:
7.     # Cheap pre-filter: element-formula gate
8.     if not formula_appears_twice_in(product, sub):
9.         continue
10.    # Full subgraph iso check — count disjoint occurrences using Maximum Set Packing
11.    # (same algorithm polymer.py:_analyze_wing_matches already uses)
12.    occurrences = count_disjoint_subgraph_isomorphisms(product, sub)
13.    if occurrences >= 2:
14.        return Group.from_molecule_subgraph(sub)   # relax labels and radicals
15. return None
```

The expensive step is `enumerate_connected_subgraphs`. Two practical optimizations:
- Stop once the first candidate satisfies the ≥2 condition (largest-first ordering).
- Cache the molecule's element-formula histogram; cheap rejection of subgraphs whose formula doesn't appear twice.

If profiling shows this inner loop dominates an RMG iteration, fall back to RDKit's MCS / connected-subgraph utilities. Document the regression and fix in a follow-up.

### 4.3 `similarity_merge(motif, pool_registry)`

```
1. for pool in pool_registry:
2.     for known_motif in [pool.monomer] + pool.feature_monomers:
3.         if motif.is_isomorphic(known_motif, strict=False):
4.             return pool
5.         # Loose match: graph-edit distance ≤ 1 (single atom or bond difference)
6.         if graph_edit_distance(motif, known_motif) <= 1:
7.             return pool
8. return None
```

Strict isomorphism handles the trivial duplicate case. Edit-distance ≤ 1 catches single-atom modifications (e.g., a hydroxyl loss yielding a near-identical aromatic motif) — those should *not* spawn a new pool, they should extend the existing pool's feature_monomer list.

### 4.4 Mass-flux gate

**ACTIVE** (spec `docs/superpowers/specs/2026-06-10-mass-flux-spawn-gate-design.md`,
as AMENDED 2026-06-10). The gate is a minimum-significance bar keeping
trivial/transient motifs from spawning pools — noise suppression, NOT slot
ranking: an important motif arriving after `max_pools` fills is a different
feature (ranking/eviction), out of scope until a real deck hits the cap with
a significant motif locked out.

**Arrival-driven re-evaluation.** Deferred candidates are never re-presented
(`new_species_list` contains only new species), so the gate re-checks a motif
only when a NEW candidate carrying it arrives. Arrivals are themselves the
relevant signal — a channel that matters keeps generating distinct
motif-carrying products — and a deferred motif's mass still flows through the
parent pool's accounting (absorbed as proxy variants), so the cost of the
known blind spot (flux grows while arrivals are quiet) is statistical
distinction, not mass conservation. The iteration-boundary ledger re-check is
the documented upgrade; **upgrade trigger:** a real deck shows a deferred
motif's flux growing while its arrivals are quiet.

**Recorded quantity** — snapshot-attributed gross mass-flux fraction. The
engine snapshot (`HybridPolymerSystem.spawn_gate_flux_snapshot()`) is a
3-tuple `(gross, pool_stats, proxy_event_mass_total)`: `gross` maps EVERY
core-species label to `max(0, core_species_production_rates[i])` (ordinary
species have real gross records — the residual's sections 3/4 maintain
production/consumption for ALL core species, simple.pyx parity, since
change (a)); `pool_stats` maps pool label to `(E[n], monomer_MW)`;
`proxy_event_mass_total` is the engine-attributed canonical-proxy sum.
Representatives are ORDINARY absorbed species, recorded in the ledger as
`(species_label, parent_pool_label)` pairs — the parent pool that absorbed
the candidate, recorded at absorption (currently `pool_registry[0]` for
gate-path candidates):

```
g_i(P) = max(0, core_species_production_rates[i]) * E[n]_P * monomer_MW_P

E[n]_P = y[mu1]_P / y[mu0]_P   if y[mu0]_P > SMALL_EPS, else 0.0

numerator(motif) = sum_{(i, P_i) in representatives(motif)} g_i(P_i)

denominator = proxy_event_mass_total
            + sum_{DEDUPED representatives (i, P_i) across the ledger} g_i(P_i)

fraction(motif) = numerator(motif) / denominator
```

Gross, never net: canonical proxies have `dn_dt ≈ 0` BY DESIGN (the archetype
apportionment reroutes their flux to pool moments) and ordinary species net
to ≈0 at steady state; the gross arrays exist precisely so diagnostics
survive that rerouting. E[n] is read LIVE from engine state at snapshot time
(never recorded-and-stale). When `mu0 <= SMALL_EPS`, `E[n]` is 0 and
`g_i = 0`: the ε-clamp errs toward deferral — the inflating direction (`mu1`
large while `mu0 ≈ 0`) is effectively unreachable inside the realizability
cone. **The denominator dedups:** a species appearing in multiple motif
entries counts ONCE in the denominator, so each motif's numerator is a
subset of denominator terms and each fraction is in [0, 1] by construction.
**Multi-motif double-counting is a stated decision, not emergent:** a
species legitimately carrying two motifs contributes its gross production to
EVERY motif-entry numerator that lists it (the motifs compete for DIFFERENT
pool slots), so the SUM of fractions across motifs may exceed 1 — accepted.
All terms are polymer-phase volumetric: V_poly cancels, no volume plumbing.

**Ledger and window.** Motifs live in `reaction_model.polymer_motif_ledger`
(`MotifLedgerEntry`); lookup is by Group **isomorphism** (the
`similarity_merge` matching idiom), never a canonical string key. One record
per motif per RMG iteration: same-iteration burst arrivals re-check the gate
against the existing window only — otherwise "windowed over N iterations"
silently becomes "windowed over N arrivals". Representatives absent from the
snapshot (absorbed this iteration, not yet simulated) contribute 0 to the
numerator — stated, not incidental. The gate statistic is the window sum
divided by the FIXED window length N (zero-filled): a single-snapshot spike
must be N× the bar to clear it; a channel persisting at fraction f for N
iterations reads f.

**Spawn-floor arithmetic.** First sighting can never spawn (zero
representatives → record 0, defer). The earliest spawn is the second
arrival's iteration, and only if that single record is itself ≥ N× the bar;
an exactly-at-bar motif needs real records filling the window, spawning at
its N-th recording iteration (k+3 for N=3 with arrivals every iteration;
~1.5× the bar spawns at k+2). Re-checks happen only at arrivals, so arrival
gaps stretch the floor further.

**Honest degradation.** No snapshot stashed (iteration 0, or no polymer
reaction system) → fraction 0.0 → defer, logged at INFO with the statistic,
the bar, and the window occupancy. No production code path fakes a number.

**Restart consequence (correct-but-loud).** The ledger is in-memory state on
the reaction model — not serialized. An RMG restart resets windows and
deferred motifs re-earn their bar: graceful, conservative, logged (same
philosophy as unstamped-reaction demotion).

Scission-daughter pools (the Path-C handshake in `make_new_reaction`) and
input-file pools do NOT route through this gate.

### 4.5 Spawn drain (in `model.py` between iterations)

```
1. for intent in queued_spawn_intents:
2.     new_pool = Polymer(
3.         label=auto_label(intent.parent_pool, intent.monomer),  # e.g., "phenol_formaldehyde_d1"
4.         monomer=intent.monomer,
5.         feature_monomer=None,
6.         end_groups=intent.end_groups,
7.         cutoff=intent.parent_pool.cutoff,
8.     )
9.     pool_registry.append(new_pool)
10.    # Allocate state-vector slots
11.    mu0_idx, mu1_idx, mu2_idx = state_vector_allocate(3)
12.    new_pool.mu_indices = (mu0_idx, mu1_idx, mu2_idx)
13.    # Initialize moments from triggering event
14.    state[mu0_idx] = intent.triggering_moles
15.    state[mu1_idx] = intent.triggering_moles * intent.triggering_dp
16.    state[mu2_idx] = intent.triggering_moles * intent.triggering_dp ** 2
17.    # Record sidecar metadata
18.    new_pool.spawn_metadata = {
19.        "spawn_iteration": current_iteration,
20.        "triggering_product_smiles": intent.triggering_product.smiles,
21.        "triggering_dp": intent.triggering_dp,
22.        "mass_flux_at_spawn": current_motif_flux(intent.monomer),
23.    }
24. # NO in-place CVODES resize (see §7). The new pool registers its
25. #   _mu0/_mu1/_mu2 dummy core species; the next RMG iteration rebuilds
26. #   the reaction system and initialize_model() resolves the new pool's
27. #   moment indices from the (now larger) core-species list by label.
28. queued_spawn_intents.clear()
```

Chip events never queue a spawn intent: `surge_chip_products` replaces the
SCISSION daughter at stamping time, upstream of this iteration-boundary hook.

## 5. Conservation invariants

For every reaction `A_chain -> B_chain + volatile`, where `A_chain` is in pool A and `B_chain` is in pool B:

**Moment effect on A (sink):**
```
dA.μ₀/dt -= r
dA.μ₁/dt -= r · ⟨DP_A⟩          where ⟨DP_A⟩ = A.μ₁ / A.μ₀
dA.μ₂/dt -= r · ⟨DP²_A⟩         where ⟨DP²_A⟩ = (Schulz-Flory closure)
```

**Schulz-Flory closure for ⟨DP²⟩:**
```
⟨DP²⟩ = (2·μ₀·μ₂ - μ₁²) / μ₀² + (μ₁/μ₀)²        (for moments-of-distribution)
       = closure_func(μ₀, μ₁, μ₂)                  in code
```

**Moment effect on B (source):**
- If B is being created in this event: B.μₖ = N · DPᵏ_product (event-spawn rule).
- If B is gaining flux from an A→B transfer: dB.μ₀/dt += r; dB.μ₁/dt += r · ⟨DP_A⟩; dB.μ₂/dt += r · ⟨DP²_A⟩.

**Volatile sink:**
- The volatile product is a regular (non-pool) species with its own ODE. Its mass flux is `r · m_volatile` (g/s), accounted for by the standard species ODE.

**Per-pool degradation source terms (`k_scission`, `k_unzip`) — added 2026-06-07.**
These are the lumped homopolymer-degradation knobs applied directly to each pool's
moments (distinct from family-generated reactions). They MUST be implemented exactly
as below; the implementation lives in `solver/polymer.pyx` and is unit-tested.

Random backbone scission, discrete-bond (Ziff–McGrady) convention — a chain of
length k has (k−1) breakable bonds, so the distribution's breakable-bond count is
(μ1−μ0):
```
dμ0/dt = +k_scission · (μ1 − μ0)     # one new chain per bond broken; SELF-LIMITING
dμ1/dt =  0                          # mass (monomer units) conserved
dμ2/dt =  (k_scission / 3) · (μ1 − μ3)
```
The `(μ1 − μ0)` form is mandatory: the naive `+k_scission·μ1` lets μ0 grow past μ1,
the state leaves the realizable cone (μ1 ≥ μ0 always holds for a k≥1 distribution),
and the μ3 closure then explodes to a DASSL singularity. With `(μ1 − μ0)` the rate
→ 0 as the pool reaches all-monomer, structurally preserving μ0 ≤ μ1.

Chain-end depropagation (unzip) — each chain loses one monomer unit per event:
```
dμ0/dt =  0
dμ1/dt = -k_unzip · μ0               # released monomer added to the gas-phase species
dμ2/dt = -k_unzip · (2μ1 − μ0)
```

Safety net: both the solver closure (`_safe_mu3_from_mu012`) and the post-processing
closure (`get_closing_moment`) return a bounded 0.0 when `μ1 < μ0`, so a transient
out-of-cone state cannot amplify into a singularity. An opt-in
`debug_check_realizability` flag on `HybridPolymerSystem` logs (once per pool) when a
moment state violates μ1≥μ0≥0 or μ0·μ2≥μ1².

**Mass-conservation invariant (test assertion):**

Define total polymer-derived mass:
```
M_polymer(t) = Σ_pools μ₁(pool, t) · m_monomer(pool)  +  Σ_explicit_oligomers n(spc, t) · m(spc)
```

Volatile mass:
```
M_volatile(t) = Σ_volatiles n(spc, t) · m(spc)
```

Invariant: `d(M_polymer + M_volatile)/dt = 0` to numerical tolerance (1e-9 for unit test, 1e-6 for integration test).

### 5.1 Discrete-reaction moment apportionment (flux archetypes)

Every proxy-touching reaction is stamped at generation time with a flux
archetype (`rmgpy.polymer.PolymerFluxArchetype`, carried on
`Reaction.polymer_flux_archetype` and resolved to solver arrays at
`initialize_model`). The residual dispatches pool moment flux on it:

| Archetype | Pool moment flux per net event rate r |
|---|---|
| SAME_POOL | none (fold-back is net-zero by construction) |
| MIGRATION | whole-chain bundle (1, E[k], E[k²]) per direction: forward rf moves source-pool statistics, reverse rr moves target-pool statistics. µ1-scaled reactions pick length-biased chains (E[k]=µ2/µ1, E[k²]=µ3/µ1, guarded closure); µ0-scaled end-group reactions pick uniformly (µ1/µ0, µ2/µ0). |
| SCISSION_FRAGMENT | complement-stays-in-parent: parent (0, −r·µ2/(2µ1), −r·(2/3)·µ3/µ1), daughter (+r, +r·µ2/(2µ1), +r·µ3/(3µ1)). Conserves µ1 exactly and adds one chain per event, independent of whether head/tail/both reactions were generated. |
| DISCRETE_CHIP | end-anchored cut ejects a stamped `a`-unit discrete chip (`Reaction.polymer_chip_units`, `a = round(chip_MW/monomer_MW)`, a = 0 legal); the complement folds back into the SAME pool. Per direction: forward (0, −rf·a, −rf·max(0, 2a·E[n]−a²)) with the µ2 decrement clamped at ≥ 0; reverse uses the EXACT extension form (0, +rr·a, +rr·(2a·E[n]+a²)) — plus a², never clamped. E[n] is scaling-aware (µ1/µ0 when µ0-scaled, µ2/µ1 when µ1-scaled — the rate scaling is the sampling measure). Closure-free (no µ3). The chip species flows through the standard gas path; exact invariant: d/dt[Σ_pools µ1 + Σ_chips a_i·n_i] = 0. µ0-scaled chips with a > 0 additionally throttle their site scaling in exhaustion — site = min(max(0, µ0), max(0, µ1)/a)/V_poly — a counting-inequality bound (every chain with n ≥ a holds at least a units, so donatable ends ≤ µ1/a always), which makes µ1 decay at worst exponentially instead of running linearly negative once the pool's units exhaust. |
| UNRESOLVED | legacy µ1-only transfer + one-time warning. Fallback for unstamped reactions (e.g. restored from pickles) AND for stamped MIGRATION/SCISSION_FRAGMENT whose src/dst pool (or DISCRETE_CHIP whose src pool) cannot be resolved in the solver (e.g. scission daughters registered as species before their pool exists) — demoting preserves the parent drain instead of silently zeroing it. |

Edge reactions are diagnostic-only (matching `simple.pyx`): their fluxes never
enter `dn_dt` or the consumption/production accumulators.

**Input hygiene:** the phenomenological `k_scission` coexists with
family-generated scission reactions as a chemically distinct (thermal random)
channel. A `k_scission` fitted to bulk degradation data that already includes
radical chemistry double-counts — specify the thermal-only rate.

Full derivations and decision rationale:
`docs/superpowers/specs/2026-06-09-proxy-moment-flux-apportionment-design.md`.

**Discreteness gate + chip surgery (2026-06-10).** `classify_structure`
rejects backbone impostors in the wing_count ≥ 2 branch (one-sided heavy-atom
bound: keep iff heavy ≥ proxy_heavy − round(0.35·proxy_heavy); reason
`backbone_impostor`). Flag-true (µ0-scaled) scission shapes are rewritten at
stamping time by `surge_chip_products` into [discrete chip,
CHIP-stamped fold-back] — no `_scission` pool is ever spawned for a chip
event (the surgery replaces the daughter before the spawn-candidate pass).
Infeasible surgery stamps UNRESOLVED, never SCISSION_FRAGMENT.

**Liveness (partially live, not dormant):** the flag fires today only for
shapes with an END_MOD product alongside the scission piece — those route
DISCRETE_CHIP immediately. Single-step end-anchored cuts (eliminations,
retro-ene) are mis-scaled µ1 today and route SCISSION_FRAGMENT until the
end-anchor detector item lands (tripwire warning
"probable mis-scaled end-anchored cut" gives the census). CHIP archetypes in
logs are expected, not a bug.

**Dormant knob:** `discrete_dp_threshold` (per pool, default 4) is parsed and
stored but has NO behavioral use under the fixed trimer proxy: the backbone
gate is proxy-relative, scission routing keys on `is_end_group_reaction`, and
the conditional DP backstop (reclassify-toward-chip when proxy repeat-count >
threshold AND piece DP < threshold) only activates for longer proxies; NO backstop code exists today — the condition is specified for the future longer-proxy work and its dormancy is pinned by `test_backstop_dormant_under_trimer_proxy`.

## 6. Sidecar JSON schema

> **SUPERSEDED for schema ≥ 2.0:** the sidecar is now schema 2.0 — envelope,
> pool additions, compiled `reactions[]` flux terms, and `conventions` are
> normatively specified in `docs/polymer_moments_format.md`. The block below
> documents the original 1.0 field set, which 2.0 preserves verbatim (note:
> the per-pool `current_moments` object sketched here was never emitted; the
> 2.0 field is `moments`, a `[mu0, mu1, mu2]` list).

Path: `<output_dir>/polymer_pools.json` alongside `chem.yaml`.

```jsonc
{
  "schema_version": "1.0",
  "generated_at": "2026-05-10T14:23:00Z",
  "rmg_iteration": 12,
  "pools": [
    {
      "label": "phenol_formaldehyde",
      "monomer_smiles": "[CH2]c1c(O)c([CH2])c(C)cc1",
      "monomer_adj_list": "1 C u1 ...",      // RMG adjacency list, exact pattern
      "feature_monomers_smiles": [],
      "end_groups": ["[H]", "[H]"],
      "cutoff": 3,
      "parent_pool": null,                     // root
      "spawn_iteration": 0,
      "spawn_event_metadata": {"source": "input"},
      "mu_indices": {"mu0_idx": 12, "mu1_idx": 13, "mu2_idx": 14},
      "current_moments": {"mu0": 1.2e-3, "mu1": 0.84, "mu2": 612.4}   // optional snapshot
    },
    {
      "label": "phenol_formaldehyde_d1",
      "monomer_smiles": "...",
      "monomer_adj_list": "...",
      "feature_monomers_smiles": [],
      "end_groups": ["[H]", "[H]"],
      "cutoff": 3,
      "parent_pool": "phenol_formaldehyde",
      "spawn_iteration": 7,
      "spawn_event_metadata": {
        "triggering_product_smiles": "...",
        "triggering_reaction_index": 314,
        "triggering_dp": 4,
        "triggering_moles": 1.42e-5,
        "mass_flux_at_spawn": 0.0123
      },
      "mu_indices": {"mu0_idx": 41, "mu1_idx": 42, "mu2_idx": 43},
      "current_moments": {"mu0": 8.7e-5, "mu1": 0.034, "mu2": 18.2}
    }
  ]
}
```

`monomer_adj_list` is the canonical machine-readable pattern. `monomer_smiles` is informational. The TA loader uses `mu_indices` to identify and skip pseudo-species in the cantera Solution, and uses pool semantics to auto-classify condensed vs. volatile.

## 7. State-vector lifecycle

```
RMG iteration N:
  1. Generate new core species via family.py
  2. process_polymer_candidates() classifies products against existing pools,
     produces spawn_intents
  3. Build / extend reaction set
  4. Run model.simulate() with current state-vector size
       — solver runs to ODE termination criteria
  5. End of iteration: if spawn_intents non-empty:
       a. Drain intents, append new pools to registry
       b. Register each new pool's _mu0/_mu1/_mu2 dummy core species
       c. Seed daughter moments from the triggering event (B.μk = N·DPᵏ)
  6. Iteration N+1: RMG rebuilds the reaction system; the larger core-species
     list yields a larger y0 and a freshly-constructed CVODES. The new pool's
     moment indices are resolved by label in initialize_model().
```

**Implementation note (reconciled 2026-06-07):** there is NO in-place state-vector
resize and no `CVodeReInit` — `reinit_solver`/`solver_reinit` do not exist. Spawning
relies on RMG's existing per-iteration teardown: the spawned pool's `_muN` dummy
species enter the core list, and the next `initialize_model()` builds a fresh solver
of the correct size. This is simpler than an in-place resize and avoids the
size-mismatch/segfault hazard the original design contemplated.

## 8. Test strategy

### 8.1 Synthetic unit test (fast — gates PRs)

`test/rmgpy/polymerTest.py::TestMultiPoolSpawn`

Construct a synthetic polymer + a hard-coded set of products that includes:
- A baseline product (should classify into parent pool, no spawn)
- A char-precursor-like product (should trigger novel-monomer discovery + spawn)
- A near-duplicate of the parent's monomer (should similarity-merge, no spawn)
- A small molecule (should classify as GAS, no spawn)

Assertions:
- `pool_registry` has exactly 2 pools after processing
- New pool's `monomer` SMILES matches expected
- Mass conservation: total moles atoms (polymer chains + volatiles) equals input total to <1e-9
- Sidecar JSON dumps and parses round-trip
- Schulz-Flory closure: assert `⟨DP²⟩ ≥ ⟨DP⟩²` (Cauchy-Schwarz; sanity check, always holds for valid distributions)

### 8.2 Carbon-phenol integration test (slow — gates releases)

`test/rmgpy/solver/solverPolymerTest.py::TestCarbonPhenolMultiPool`

Run a small carbon-phenol model with limited families and constraints (so it finishes in ~5 min on CI). Assertions:
- `polymer_pools.json` is emitted with `schema_version` set
- `len(pools) >= 2`
- At least one pool has `parent_pool == "phenol_formaldehyde"`
- Mass-balance invariant holds across the run within 1e-6
- No NaN in any μ trajectory
- `<pool_label>_muN` pseudo-species appear in cantera YAML output

### 8.3 Regression

The existing `examples/rmg/polystyrene/input.py` must still produce the same output as before this work, since polystyrene with PS-only is the single-pool case (≡ multi-pool with N=1, no spawn intents triggered). Add a snapshot test if not already present.

## 9. Implementation order with TDD checkpoints

1. **Synthetic unit test (red)** — `TestMultiPoolSpawn` written; all assertions fail because nothing implements them yet. Commit.
2. **Pool registry refactor** — `polymer_input.py` builds a list of pools. Existing single-pool behavior preserved by registry of length 1. All existing tests pass. Commit.
3. **Auto-label generator** — `auto_label(parent, motif) -> str` implemented + tested. Commit.
4. **`discover_repeat_motif`** — implement; unit test against simple chains (e.g., 3× propylene oxide). Commit.
5. **`similarity_merge`** — implement; unit test on near-duplicate motifs. Commit.
6. **Mass-flux gate** — implement rolling accumulator; unit test. Commit.
7. **Spawn-intent generation in `process_polymer_candidates`** — wire phases A–E. Synthetic unit test starts passing partial assertions. Commit.
8. **Solver picks up new pools on rebuild** — daughter `_muN` species register; `initialize_model` resolves their indices by label on the next iteration (no in-place resize/CVODES reinit — that approach was dropped). Synthetic test asserts moment continuity across spawn. Commit.
9. **Schulz-Flory closure helpers** — implement; unit test against analytic Schulz-Flory chains. Commit.
10. **Sidecar JSON writer** — `yml.py` and/or `chemkin.py` extension. Synthetic test asserts JSON round-trip. Commit.
11. **Family.py contagion update** — pass full pool registry. Existing tests must still pass. Commit.
12. **Integration test (red → green)** — `TestCarbonPhenolMultiPool` written; iterate on real run until green. Commit.
13. **Polystyrene regression** — add snapshot if missing; verify unchanged. Commit.
14. **Doc update** — point `documentation/source/...` to this design doc; update `examples/rmg/polystyrene/input.py` comments if needed. Commit.

Each numbered step is its own commit. The synthetic unit test from step 1 acts as the running progress meter.

## 10. Known limitations (phase 1)

1. Daughter pools inherit parent end_groups verbatim. Char-precursor termination chemistry differs from novolak; this is acknowledged. Phase 2 work item: structural inference of daughter end_groups.
2. Novel-motif discovery's connected-subgraph enumeration is O(2^n) worst case. We order largest-first and exit on first ≥2-occurrence hit, but pathological products (highly symmetric, large) may be slow. Profile and switch to RDKit MCS only if needed; document.
3. Schulz-Flory closure is exact only for free-radical chain-growth polymerization. Phenolic resin is step-growth — Flory's `2-1/μ₀` form applies. Closure functions documented inline; alternative closures (log-normal) are pluggable but not phase-1.
4. `max_pools` cap is hard. After cap reached, would-be-spawned products are still tagged `is_polymer_proxy=True` and tracked as explicit oligomers; their moments are not separately accounted. Acceptable for phase 1; phase 2 may add cap-relaxation under user override.
5. The cantera writer emits `<label>_muN` as pseudo-species. Cantera does not "understand" them — they are decorative in YAML and meaningful only via the sidecar. TA's loader filters them out before constructing a Cantera `Solution`.
6. **Crosslink / chain–chain coupling is not modeled.** A reaction whose product is a crosslink (>2 intact wings, i.e. two chains joined) is *rejected* (`PolymerCrosslinkError` → the reaction is discarded) rather than represented, because the single-distribution-per-pool moment model has no coupling term. This conserves mass but drops char-forming coupling chemistry — a phase-2 item.
7. **Gas/polymer classification is purely topological** (wing count). There is no volatility/MW check, so a small volatile that happens to contain an [end-group + monomer] subgraph can be retained as a polymer fragment. Acceptable for phase 1; revisit with a size/boiling-point gate. *Partially addressed 2026-06-10:* the discreteness gate rejects backbone impostors (heavy-atom bound) and DISCRETE_CHIP tracks end-cut fragments discretely; true volatility (Tb vs reactor T) remains a mass-transfer-layer item.
8. **The proxy is a fixed single-buffer trimer** (`[Head]–O–O–O–[Tail]`), with no length knob — `cutoff`/`x_s` set the solver's explicit/tail boundary, NOT the proxy length. Interior chemistry and thermo are therefore taken from a center monomer that sits one monomer from a real chain end; for short repeat units this may let end effects leak in.
9. **Cantera export drops element-unbalanced polymer-proxy reactions** (the size-changing `*_scission_tail/_scission_head` daughters) so `ct.Solution` doesn't reject the mechanism via `checkBalance`. Their mass is carried by the moment solver + `polymer_pools.json` sidecar, not by Cantera — the TA consumer must rely on the sidecar for that flux.
10. **Scission daughters seed as EMPTY pools** (`moments=None`, `initial_mass=0` → zero moments, halved Mn/Mw) and fill via reaction flux. They do NOT inherit the parent's distribution (that would duplicate mass).
11. **Chip fold-back ignores the new chain end created by the cut** (a
    potential unzip initiator) — the same approximation END_MOD fold-back
    already makes everywhere; end-specific kinetics belong to k_unzip/END_MOD
    chemistry, not a spawned pool. (A1)
12. **End-group mass drift on chips** (chip carries the old cap, parent image
    keeps both) is sub-monomer and consistent with existing µ-accounting;
    shape-(a)'s Cantera overbalance is this same drift made visible at the
    species-balance level — dropped from the export and counted, like
    registered scissions. (A2)
13. **True volatility (Tb vs reactor T) is out of scope** — mass-transfer
    layer, phase 2; see also the end-anchor detector follow-up in §11. (A3)
14. **The legacy UNRESOLVED path has the same latent exhaustion structure**
    for µ0-scaled reactions (rate ∝ µ0, drains µ1, µ0 untouched) as an
    unthrottled chip — documented, NOT fixed: the UNRESOLVED contract is
    bit-exact legacy behavior, and the end-anchor detector (§11) migrates
    exactly these reactions to DISCRETE_CHIP, where the exhaustion throttle
    protects them. The unprotected population shrinks by design. (A4)

## 11. Open items for follow-up

- **End-anchor detector** (scoped 2026-06-10, deliberately deferred): extend
  `is_end_group_reaction` determination to template/site anchoring (reacted
  site in the cap or cap-adjacent unit, from template-labeled atoms before
  `clear_labeled_atoms`). Chemistry-visible: flips affected reactions
  µ1 → µ0 scaling (rate × E[n], 2–4 orders of magnitude) — needs its own
  EPDM delta analysis, reaction by reaction. Acceptance criteria: (a) the
  tripwire census on real decks drops to ~zero; (b) the EPDM delta is
  explained reaction-by-reaction. The seam is engineered so the detector
  touches only detection + re-scaling: routing, bundles, and solver tests
  downstream of the flag are already correct (CHIP activates with zero
  solver changes).
- Benchmark `discover_repeat_motif` on the carbon-phenol case; decide on RDKit fallback yes/no with data.
- Decide on rate-rule heuristics for inter-pool transfer reactions when family.py doesn't naturally produce them.
- TA loader implementation (TA repo's Phase 2) — proceed in parallel against this contract.
