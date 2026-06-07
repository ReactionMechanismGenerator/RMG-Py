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

```
state: per-motif accumulator updated during reaction generation
       motif_flux[motif_signature] = sum of |product moles produced via reactions yielding this motif|
                                   over the trailing N RMG iterations (default N=3)

mass_flux_gate(motif, threshold):
    return motif_flux[motif.signature] / total_polymer_derived_mass >= threshold
```

The signature is a canonical SMARTS or RMG `Group` adjacency-list hash. Trailing-N rolling window prevents a single transient peak from triggering a spawn.

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

## 6. Sidecar JSON schema

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
7. **Gas/polymer classification is purely topological** (wing count). There is no volatility/MW check, so a small volatile that happens to contain an [end-group + monomer] subgraph can be retained as a polymer fragment. Acceptable for phase 1; revisit with a size/boiling-point gate.
8. **The proxy is a fixed single-buffer trimer** (`[Head]–O–O–O–[Tail]`), with no length knob — `cutoff`/`x_s` set the solver's explicit/tail boundary, NOT the proxy length. Interior chemistry and thermo are therefore taken from a center monomer that sits one monomer from a real chain end; for short repeat units this may let end effects leak in.
9. **Cantera export drops element-unbalanced polymer-proxy reactions** (the size-changing `*_scission_tail/_scission_head` daughters) so `ct.Solution` doesn't reject the mechanism via `checkBalance`. Their mass is carried by the moment solver + `polymer_pools.json` sidecar, not by Cantera — the TA consumer must rely on the sidecar for that flux.
10. **Scission daughters seed as EMPTY pools** (`moments=None`, `initial_mass=0` → zero moments, halved Mn/Mw) and fill via reaction flux. They do NOT inherit the parent's distribution (that would duplicate mass).

> NOTE (§4.4 mass-flux gate): NOT YET ACTIVE. `_estimate_relative_flux` is a `0.5`
> stub, so spawning is currently gated only by `max_pools` + similarity-merge. The
> rolling-window `MassFluxAccumulator` is implemented but not wired into the live path.

## 11. Open items for follow-up

- Benchmark `discover_repeat_motif` on the carbon-phenol case; decide on RDKit fallback yes/no with data.
- Decide on rate-rule heuristics for inter-pool transfer reactions when family.py doesn't naturally produce them.
- TA loader implementation (TA repo's Phase 2) — proceed in parallel against this contract.
