# Design — proxy reaction moment-flux apportionment (µ0/µ1/µ2)

**Date:** 2026-06-09
**Branch:** `polymer`
**Status:** approved design, pending implementation plan
**Open item resolved:** #4 from `~/handoffs/handoff-polymer-branch-open-design-items.md`
("proxy reaction flux updates only µ1, not µ0/µ2").
**Includes:** edge-reaction `dn_dt` gating fix (pre-existing divergence from
`simple.pyx`, discovered while scoping this item).

## Problem

When a reaction touches a pool-proxy species, the polymer hybrid solver
(`rmgpy/solver/polymer.pyx`, RHS residual ~lines 942–981) applies the whole
molar event flux `r_mol_s` to that pool's **µ1 only**: reactant side `-r`,
product side `+r`. µ0 and µ2 are never touched by discrete reactions — they
evolve only through the dedicated `k_scission`/`k_unzip` pool ODEs.

Consequences by reaction shape:

- **Same-pool round trips** (BASELINE/FEATURE/END_MOD products folding back
  into the parent pool): `-r` and `+r` cancel on the same µ1 slot. Correct
  today by accident of symmetry.
- **Cross-pool reactions** (product classifies into a different pool —
  `_scission_head`/`_scission_tail` daughters, or a modification that
  reclassifies the chain into another pool): exactly 1 monomer-unit-mole of µ1
  moves per event-mole, and **zero chains (µ0) and zero µ2** move. A daughter
  pool accumulates µ1 against a frozen near-zero µ0 — its Mn drifts toward
  infinity — and the parent never loses a chain. Mn/Mw/PDI (the model's
  headline outputs) desynchronize.
- **Edge reactions** (all-core reactants, ≥1 edge product): `polymer.pyx`
  applies their `dn_dt` and consumption/production fluxes unconditionally,
  whereas standard RMG (`simple.pyx:464-502`) treats edge reactions as
  diagnostic-only. Edge fluxes currently perturb the integrated core state,
  including pushing µ1 flux into pools.

## Decision summary

Resolved during brainstorming (2026-06-09):

1. **Event semantics — whole-chain migration.** When a discrete event makes a
   chain change pool identity, one event moves a *whole chain's* statistics,
   not one monomer unit. The event rate stays scaled by the picking moment
   (µ1 site density, or µ0 for end-group reactions): a longer chain is
   proportionally more likely to be picked, and the picked chain's statistics
   are correspondingly length-biased.
2. **Scission ownership — full chain-split bookkeeping with
   "complement stays in parent" accounting.** Family-generated scission
   reactions carry real moment flux (Option A). Each registered scission
   reaction has only ONE polymer fragment product (the reason they are
   Cantera-unbalanced), so the accounting treats the un-represented complement
   fragment as remaining in the parent pool. Each reaction then independently
   conserves monomer units and adds exactly one chain per event — no
   double-counting if RMG generates both head and tail reactions for the same
   physical cut, no leak if it generates only one. The phenomenological
   `k_scission` ODE coexists as a chemically distinct channel (thermal random
   scission vs. radical-mediated site-specific scission). **Input-hygiene
   caveat for docs:** a user who fitted `k_scission` to bulk data that already
   includes radical chemistry double-counts.
3. **Mechanism — generation-time archetype stamping (Approach 1).** Per-
   reaction flux archetype stamped at reaction generation (where the graph
   classification lives), carried on the `Reaction` object, resolved to solver
   arrays at `initialize_model`, dispatched in the RHS. Rejected: solver-side
   inference (cannot distinguish MIGRATION into a daughter pool from
   SCISSION_FRAGMENT — only `_reacted_class` knows; label-based inference was
   already rejected during the µ0-scaling work) and post-hoc moment repair
   (erases the PDI evolution the moment model exists to track).
4. **Edge fix included.** Gate `dn_dt` and `core_species_consumption/
   production_rates` writes on `r_idx < n_rxn`, matching `simple.pyx`.

## Per-archetype moment flux (the chemistry)

Let `r` = net event flux (mol events/s, `r_mol_s`; may be negative). Chain-pick
statistics from the **source** pool, by the reaction's scaling moment:

| Scaling | E[1] | E[k] | E[k²] |
|---|---|---|---|
| µ1-scaled (site-picked → length-biased chain) | 1 | µ2/µ1 | µ3/µ1 |
| µ0-scaled (`is_end_group_reaction`, uniform chain) | 1 | µ1/µ0 | µ2/µ0 |

µ3 comes from the existing guarded closure `_safe_mu3_from_mu012`.

### SAME_POOL

BASELINE/FEATURE/END_MOD product folds back into the reactant's pool. Net
moment flux is exactly zero — **explicitly skip** pool moment writes (today's
`-r`/`+r` on µ1 cancels anyway; skipping avoids roundoff and closure calls).
Gas co-reactant/product fluxes and proxy diagnostics unchanged.

### MIGRATION

Intact-chain product classifies into pool B ≠ reactant pool A. The bundle is
applied **per direction**, not on the net rate: the forward term moves
A-statistics chains A→B and the reverse term moves B-statistics chains B→A.
With b_A, b_B = (1, E[k], E[k²]) evaluated from each pool and rf, rr the
forward/reverse molar event rates (already computed separately in the
residual, `polymer.pyx:895-908`):

    A: −rf·b_A + rr·b_B        B: +rf·b_A − rr·b_B

Conservation is exact per direction. This costs one extra bundle evaluation
only when rr > 0 and removes any bundle-source-by-sign approximation — a
near-equilibrium cross-pool reaction correctly exchanges long chains one way
and short chains back even when the net µ0 flux is ~0.

Realizability: b satisfies E[k²] ≥ E[k]² (Cauchy–Schwarz), so adding it to a
realizable pool stays in the cone; the closure guard and
`debug_check_realizability` remain the safety net.

### SCISSION_FRAGMENT

Product polymer stamped `PolymerClass.SCISSION`, mapping to daughter pool B.
µ1-scaled pick (length-biased parent chain n: E[n] = µ2/µ1, E[n²] = µ3/µ1),
uniform cut position (fragment a uniform along the chain: E[a] = E[n]/2,
E[a²] = E[n²]/3). The complement (length n − a) stays in the parent.

| | Δµ0 | Δµ1 | Δµ2 |
|---|---|---|---|
| Parent A | 0 | −r·µ2/(2µ1) | −r·(2/3)·µ3/µ1 |
| Daughter B | +r | +r·µ2/(2µ1) | +r·µ3/(3µ1) |

Consistency checks: µ1 exactly conserved; total µ0 grows by +r per event (one
cut = one new chain); total µ2 changes by −r·µ3/(3µ1) < 0 (scission narrows
the distribution). The accounting is per-reaction-independent: head-only,
tail-only, or both-generated cases all conserve.

**Stated approximation:** for r < 0 (net reverse = chain coupling, physically
negligible under degradation conditions and rejected as CROSSLINK at
generation anyway), the same parent-statistics formulas are sign-flipped.
SCISSION_FRAGMENT uses the net rate (unlike MIGRATION's per-direction split)
because the reverse statistics of a coupling event are not representable.

**µ0-scaled scission is excluded.** The table above assumes a µ1-scaled pick
and a uniform cut position. An end-initiated scission
(`is_end_group_reaction` ∧ SCISSION product) violates both assumptions
(uniform-by-chain pick, cut near the chain end); silently applying the table
would systematically overestimate daughter Mn. Generation routes this
combination to **UNRESOLVED** with a warning. If the warning fires on real
decks, a second bundle (fixed small fragment length ℓ: daughter
(+r, +rℓ, +rℓ²), parent (0, −rℓ, −r·(2ℓ·µ1/µ0 − ℓ²))) is the follow-up.

### UNRESOLVED / legacy fallback

If generation cannot stamp a single unambiguous archetype (more than one
polymer product with a cross-pool member, end-initiated scission, |R| > 1
with a cross-pool product), or a proxy-touching reaction arrives unstamped
(e.g. restart from an older pickle), the solver applies **today's legacy
µ1-only transfer** and logs a one-time warning. Real flux is never silently
zeroed.

The unstamped-arrival remap (touches-a-proxy ∧ archetype == NONE →
UNRESOLVED) happens once in `initialize_model` while filling the arrays —
the RHS hot loop never branches on it, and the one-time warning is trivial
there.

Follow-up if the |R| == 2 warning fires often on real decks (e.g. radical
transfer between chains in different pools): `create_reacted_copy` lineage
can in principle resolve which reactant each product came from; out of scope
for this phase.

## Metadata channel and data flow (4 layers, mirroring the µ0-scaling work)

```
create_reacted_copy stamps product._reacted_class      (already exists)
        │
        ▼
rmgpy/polymer.py: classify_reaction_flux_archetype(reactants, products) -> int
        │   Let R = set of reactant polymer pools; a product polymer is
        │   "cross-pool" iff its pool ∉ R.
        │   any product SCISSION ∧ NOT is_end_group     -> SCISSION_FRAGMENT (3)
        │   any product SCISSION ∧ is_end_group         -> UNRESOLVED (4)
        │       (end-initiated scission: bundle table assumptions invalid)
        │   exactly 1 polymer product total, it is
        │   cross-pool, AND |R| == 1                    -> MIGRATION (2)
        │   no cross-pool products                      -> SAME_POOL (1)
        │       (covers multi-proxy reactions like inter-chain H-abstraction,
        │        where each product folds back into its own reactant pool)
        │   anything else (>1 polymer product with a
        │   cross-pool member, >1 cross-pool, or
        │   ambiguous source with |R| > 1)              -> UNRESOLVED (4) + warn once
        ▼
rmgpy/reaction.pxd/.py: cdef public int polymer_flux_archetype  (default 0=NONE)
        │   set in rmgpy/rmg/model.py make_new_reaction after each handshake
        │   (both directions), beside the existing is_end_group_reaction stamp;
        │   serialization mirrors is_end_group_reaction
        ▼
rmgpy/solver/polymer.pyx initialize_model:
        reaction_flux_archetype[n]  (int8, chain(core, edge) order)
        reaction_src_pool[n], reaction_dst_pool[n]
            resolved from species_to_pool_indices (NO label matching);
            for MIGRATION/SCISSION_FRAGMENT src = the single reactant
            proxy's pool, dst = the cross-pool/daughter product's pool
        ▼
RHS residual: replace the three µ1-only writes (~943/952/971) with ONE
per-reaction archetype dispatch applying the bundles above. Per-slot gas
fluxes and proxy diagnostics (proxy_activity, consumption/production rates)
keep their current per-slot structure.
```

Archetype enum values: 0 NONE (no proxy involvement), 1 SAME_POOL,
2 MIGRATION, 3 SCISSION_FRAGMENT, 4 UNRESOLVED.

The reaction-rate diagnostic block (~line 1269) computes rates only, not
`dn_dt`; it keeps its existing µ0/µ1 scaling and needs no apportionment
change (verify during implementation).

## Guards

Each guard follows the **direction of depletion**, not just the nominal
source pool:

- MIGRATION: the rf term is guarded by A's denominator moments, the rr term
  by B's — each direction is skipped independently if its own source pool's
  denominator (µ1, or µ0 for end-group bundles) is below `SMALL_EPS`·V_poly
  (no chains to move; prevents 0/0).
- SCISSION_FRAGMENT: for r ≥ 0 guard the parent's µ1; for r < 0 the daughter
  is the depleted pool (Δµ0 = −|r| on B), so additionally require the
  daughter's µ0 and µ1 above threshold — otherwise an empty daughter is
  driven µ0-negative through a guard that only looked at the parent.
- µ3 not finite from `_safe_mu3_from_mu012`: skip **only the µ2 component**,
  still apply µ0/µ1 components (mirrors `polymer.pyx:1052`; keeps
  `test_residual_stays_finite_for_extreme_moment_states` green).
- All moment reads `max(0, y[idx])`-clamped, as elsewhere in the residual.
- The parent-side scission µ2 drain uses closure-estimated µ3; a closure
  overestimate drains parent µ2 too fast, so the `max(0,·)` clamp and
  `debug_check_realizability` carry real weight there (targeted test below).

## Jacobian

Verified 2026-06-09: `HybridPolymerSystem` defines no `jacobian` attribute
(`hasattr` is False; the analytic Jacobian path belongs to `SimpleReactor`
only), so pydas DASSL uses a numeric finite-difference Jacobian — nothing to
update for the new cross-pool µ0/µ2 coupling. Expect some added stiffness
from the µ3-closure derivative near the realizability boundary; the numeric
Jacobian absorbs it, but integration-step behavior on the EPDM deck is part
of end-to-end verification.

## Edge-reaction gating fix

Gate `dn_dt` writes **and** `core_species_consumption/production_rates`
accumulation on `r_idx < n_rxn`, matching `simple.pyx:464-502` where edge
reactions contribute only to `edge_reaction_rates`/`edge_species_rates`.
`edge_reaction_rates`, `edge_species_rates`, and `proxy_activity` keep their
current behavior. This also closes the "edge polymer product" hole: core
reactions cannot have edge products, so bundles only ever apply between
resolved pools.

## Testing

- **polymer.py unit:** `classify_reaction_flux_archetype` — fold-back →
  SAME_POOL; cross-pool FEATURE → MIGRATION; scission product →
  SCISSION_FRAGMENT; two cross-pool products → UNRESOLVED + warning.
- **reactionTest:** `polymer_flux_archetype` field default + assignment
  (mirror `test_is_end_group_reaction_default_and_kwarg`).
- **modelTest/polymerTest:** `make_new_reaction` stamps the archetype (mirror
  the µ0-scaling handshake-flag tests).
- **solverPolymerTest:**
  - Migration: two pools + one migration reaction; assert dµ0/dµ1/dµ2 against
    the bundle formulas analytically and exact A-loss == B-gain.
  - Scission: parent + daughter; assert the SCISSION_FRAGMENT table; total µ2
    decreases.
  - Same-pool: net-zero pool moments (regression of today's behavior).
  - Guards: empty source pool → zero bundle, no NaN; µ3 = inf → µ0/µ1 applied,
    µ2 skipped (extend the extreme-state finiteness sweep).
  - Edge gating: an edge reaction with core reactants → zero `dn_dt` for its
    reactants while `edge_species_rates` is nonzero.
  - Realizability trajectory: forward-Euler mini-trajectory with mixed
    migration + scission stays in the realizable cone (mirror
    `test_scission_moment_trajectory_stays_realizable`), with a **global
    conservation assertion** (Σ over pools of µ1 + gas-phase monomer moles =
    constant).
  - Monodisperse limit: a PDI→1 pool has closed-form scission moment answers
    — sharp analytic check on the bundle/closure interaction.
  - Closure-overestimate clamp: parent µ2 driven down by an inflated µ3
    estimate stays non-negative and warns via `debug_check_realizability`.
  - Unstamped remap: a proxy-touching reaction arriving with archetype NONE
    is remapped to UNRESOLVED at `initialize_model` (legacy µ1 flux applied,
    one-time warning logged).
- **End-to-end:** EPDM deck re-run (`~/runs/RMG/epdm_v0_2026-06-06b`) —
  scission-ODE dominated, expect `MODEL GENERATION COMPLETED`, 26 species /
  28 reactions unchanged.
- `make` after every `.pyx`/`.pxd` edit before running tests.

## Out of scope

- Reverse-rate site scaling for proxy reactions (`rr *= site` uses the
  reactant pool's density; pre-existing approximation, unchanged — the
  MIGRATION per-direction bundles consume rf/rr as computed).
- The `k_scission`/`k_unzip` pool ODEs (already correct; untouched).
- Item #5 (gas/polymer volatility gate), #8 (mass-flux spawn gate), #9
  (motif discovery) from the handoff.
- Representing scission complements as explicit products (would require
  balanced two-polymer-product reactions; the complement-stays accounting
  makes this unnecessary for moment bookkeeping).
