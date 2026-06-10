# Discreteness Gate + DISCRETE_CHIP Archetype — Design Spec

**Date:** 2026-06-10
**Branch:** `polymer`
**Resolves:** open item #5 (gas/polymer split is topological only) and, for the
reachable shapes, residual (e) from the moment-flux apportionment work
(end-initiated scission accounting).
**Builds on:** `docs/superpowers/specs/2026-06-09-proxy-moment-flux-apportionment-design.md`
(archetype stamping, `_chain_bundle`, complement-stays convention).

---

## 1. Problem and reframing

`classify_structure` (`rmgpy/polymer.py:1599-1664`) decides gas vs polymer
purely topologically (wing count). Two failure modes:

1. **Backbone impostors.** A small molecule that accidentally contains both
   [end-group + monomer] wing subgraphs classifies as an intact-backbone
   polymer image and its mass is captured by the pool. (Verified live case:
   bibenzyl double-wing-matches the PS proxy today. The original item-#5
   example — propane/butane vs PE — turned out on probing to classify via the
   wing_count == 1 scission path instead, where the §3.2 routing governs it;
   see §10-V6.)
2. **Small fragments as pools.** A scission piece of literal DP ~1-2 folded
   into pool bookkeeping adds µ0 = 1 with tiny µ1, crashing Mn; monomer/dimer
   yields are headline pyrolysis outputs that should be explicit species.

**Reframing (decided):** the gate decides *discreteness* — is this product
representable by moment statistics, or should it be a discrete tracked
species — structurally and temperature-independently, at generation time.
Actual volatility (Tb estimate vs reactor T) is a property of the
**mass-transfer layer** and is explicitly out of scope (phase 2). A discrete
chip is discrete whether it is volatile or condensed; the moment accounting
below is correct either way. Hence the archetype name DISCRETE_CHIP, not
EVAPORATIVE_FRAGMENT — so nobody later "fixes" a condensed chip back into the
pool path.

**Core principle (recurs throughout):** literal size of a probe is a property
of the *representation*, not of any real species. The baseline proxy is a
deliberately small trimer surrogate (`_stitch_trimer`, polymer.py:763-791:
head + 3 repeat units + tail). A size rule is valid only where the
representation can actually express the distinction the rule tests.

## 2. Decisions (settled in brainstorm, 2026-06-10)

| # | Decision | Rationale |
|---|----------|-----------|
| D1 | Gate criterion = discreteness, not volatility; Tb/evaporation deferred to mass-transfer layer | T-independent, generation-time safe; "styrene trimer is condensed at 298 K but a classic volatile at 700 K" is a phase question, not a representability question |
| D2 | Backbone branch (wing_count ≥ 2): proxy-relative size gate, NOT the DP knob | A wing_count ≥ 2 candidate is an image of the proxy representation; its literal DP (~3) is a truncation artifact. The only real question is "genuine chain image vs accidental double-wing match", which is inherently proxy-relative |
| D3 | Scission branch (wing_count == 1): discriminator is **`is_end_group_reaction`** (per-chain vs per-site rate scaling), NOT piece-size/buffer counting | Piece counts on a 3-unit proxy carry no positional information: a u1–u2 cut (the image of *every* interior backbone bond) yields 1- and 2-unit pieces, both artifacts. No buffer value classifies the trimer correctly. The rate scaling IS the sampling measure, so keying routing on it makes rate and statistics consistent by construction |
| D4 | Chip complement folds back to the parent pool (no `_scission` pool spawn) | A chip event shortens the chain by a known small `a`; the remainder's statistics are ≈ the parent's. Spawning a pool indistinguishable from the parent is pure proliferation, and double-represents the remainder (the bundle already encodes complement-stays) |
| D5 | Chip µ2 term uses scaling-aware E[n]: uniform µ1/µ0 when µ0-scaled, length-biased µ2/µ1 when µ1-scaled | The rate scaling is the sampling measure; pairing a µ0 rate with length-biased E[n] is an internally inconsistent estimator that systematically over-drains parent µ2 |
| D6 | Detector gap accepted: DISCRETE_CHIP keys on the flag **as implemented** (END_MOD product stamps); extending detection to template/site anchoring is a separate scoped follow-up item | Flipping µ1→µ0 scaling changes rates by E[n] (2–4 orders of magnitude) — a chemistry-visible change that demands its own EPDM delta analysis and would destroy this item's "EPDM unchanged" regression gate. A mis-scaled end cut today gets µ1 rate AND uniform-cut statistics: wrong, but coherently wrong, attributable to the pre-existing scaling bug |
| D7 | `discrete_dp_threshold` knob: defined in pool config (default 4), **dormant** under the fixed trimer proxy; no C40 heavy-atom ceiling | Backbone branch doesn't use it (D2); the flag dichotomy does the scission work (D3); the backstop (D8) only activates for longer proxies. Defined-but-documented beats undefined intent; the ceiling has no reachable trigger |
| D8 | DP backstop is conditional: applies only when proxy repeat-count > threshold | "Center-touching piece with DP < threshold → exact-a accounting" presumes the proxy can express "small but interior". The trimer cannot (every piece has DP ≤ 3 < 4); unconditional, the backstop would route ALL mid-chain scission to chip and kill µ0 growth/Mn halving entirely |

## 3. Classification changes (`rmgpy/polymer.py`)

### 3.1 Backbone impostor gate (wing_count ≥ 2 branch)

In `classify_structure`, inside the `wing_count >= 2` branch, **before** the
BASELINE isomorphism check (an exact-isomorphic baseline trivially passes
anyway, but gating first keeps one decision point):

```
keep iff heavy_atoms(candidate) >= heavy_atoms(baseline_proxy) - tol
tol = round(0.35 * heavy_atoms(baseline_proxy))
```

- Below the bound → `PolymerClass.GAS`, reason `"backbone_impostor"`.
- **One-sided on purpose:** legitimate images can be larger than the proxy
  (FEATURE grafting a peroxy/side group) or modestly smaller (H loss,
  side-group elimination). Heavy atoms, not MW: H-insensitive, and a lost
  side group is a known heavy-atom delta.
- The exact tolerance is not delicate: an impostor that matched both wings
  sits far below 0.65× of the proxy (bibenzyl: 14 heavy atoms vs the
  PS-proxy bound of 16; see §10-V6 — a head-wing+tail-wing directly-joined
  molecule is mathematically never below the bound, so the gate targets
  accidental overlap matches like this one).
- **No ceiling**, contingent on verification (see §10-V3) that
  polymer+polymer coupling cannot reach this branch (CROSSLINK rejection
  upstream raises `PolymerCrosslinkError` in `create_reacted_copy` before
  classification of coupling shapes). If a reachable ≥2×-proxy path is found
  during implementation, add a ceiling at that point with the same
  proxy-relative form.

PS arithmetic pinning D2: a genuine PS proxy image with H caps is
~314 g/mol → literal DP = round(314/104) = 3 < 4; a DP knob here would
gas-classify the pool's own bookkeeping device.

### 3.2 Scission-shaped routing (wing_count == 1) — discriminator swap

`classify_structure` itself is **unchanged** for wing_count == 1 (still
returns `SCISSION`); routing happens at reaction stamping time (§4), because
the discriminator `is_end_group_reaction` is a reaction-level fact not
available per-product:

- scission-shaped ∧ flag **true** → **DISCRETE_CHIP** (product surgery §4.2)
- scission-shaped ∧ flag **false** (µ1-scaled) → **SCISSION_FRAGMENT**,
  *regardless of literal piece size* (artifact, per D3)

Piece structure survives in exactly two roles:
1. **Supplying `a`** (the chip's monomer-unit count) — for end-anchored
   chemistry the represented piece is exact: u1 really is one unit from the
   end on the real chain too.
2. **Diagnostics** (the tripwire, §4.4) — never routing.

### 3.3 New `PolymerClass.CHIP` member

Stamped (as `_reacted_class`) on the fold-back parent copy produced by chip
product surgery (§4.2), so `classify_reaction_flux_archetype` can recognize
the shape. Not produced by `classify_structure` itself.

### 3.4 a = 0 chips are legal

Bare end-cap ejection (e.g. CH3 loss) gives `a = 0`: the archetype fires with
zero µ1/µ2 drain — net pool effect ≈ SAME_POOL, chain count unchanged. No
special-casing beyond the formulas in §5 (they degrade to zero naturally).

## 4. Generation stamping (`rmgpy/polymer.py`, `rmgpy/rmg/model.py`, `rmgpy/reaction.pxd/.py`)

### 4.1 Archetype enum + classifier branch

- `PolymerFluxArchetype` gains `DISCRETE_CHIP = 5`; solver mirror constant
  `FLUX_DISCRETE_CHIP = 5` (extend the constants-match-enum test).
- `classify_reaction_flux_archetype`: a CHIP-stamped product polymer →
  `DISCRETE_CHIP`. **This check must precede the SCISSION branch:** after the
  (b) surgery (§4.2) the product list has no END_MOD member, so the
  classifier's internal `is_end_group_reaction(products)` call would return
  False — the CHIP short-circuit keeps the classifier from ever reaching that
  recompute for surged shapes.
- **The "end-initiated scission → UNRESOLVED" branch does not survive in its
  old form.** Flag-true scission shapes now resolve to DISCRETE_CHIP;
  flag-false ones go to SCISSION_FRAGMENT. UNRESOLVED remains only for
  invariant-violation shapes (§4.2 surgery failure). This spec supersedes the
  2026-06-09 spec's routing for that one shape.

### 4.2 Product surgery + represented-product invariant

At stamping time in `make_new_reaction` (after `_handshake_structures` and
the `is_end_group_reaction` flag computation, both already in place), a
flag-true reaction with a SCISSION-stamped product is one of **two
sub-shapes**, both driven to the same end state
**[discrete chip, CHIP-stamped fold-back]**:

- **(a) Represented polymer = macro daughter** (no END_MOD product; the
  SCISSION-stamped piece is the large remainder): replace the daughter
  Polymer product with `parent.copy(deep=True)` stamped
  `_reacted_class = PolymerClass.CHIP` (fold-back, mirroring the END_MOD
  fold-back mechanism), and derive/register the chip as a discrete
  (non-Polymer) product (its structure is derivable from the proxy graph).
- **(b) SCISSION piece = chip, END_MOD fold-back present** — the only
  flag-true shape live today (§4.5): the small SCISSION-stamped piece IS the
  chip and the macro fold-back already exists as the END_MOD product.
  Demote the chip to a discrete Molecule (undo its handshake conversion) and
  **re-stamp the existing END_MOD fold-back as CHIP**. Applying (a)'s
  procedure here would replace the chip with a second fold-back — losing the
  chip and double-folding the parent.

Chip identification within a sub-shape: the smaller piece (for an
end-anchored cut the pieces are end+ε and rest) — size used for
*identification* only, not positional routing; the position is already known
from the flag.

If the generation shape makes the surgery infeasible (chip unrepresentable,
ambiguous piece identification, neither sub-shape's preconditions met) →
stamp **UNRESOLVED + warn-once**, never SCISSION_FRAGMENT (which would apply
uniform-cut statistics to a near-end cut AND leave the chip mass
unaccounted).

**Flag-stability riders (b):** the re-stamp happens *after*
`is_end_group_reaction` is computed and stored on the reaction, so the stored
flag survives — and **nothing downstream may recompute the flag from product
stamps** (the re-stamp removes the END_MOD member, so a recompute would flip
it). The solver reads the stored field; the classifier must check
CHIP-stamped products **before** its SCISSION branch (§4.1), so its internal
`is_end_group_reaction(products)` call is never reached for surged shapes.

**Never-queue (verified, §10-V4):** spawn intents are created by the
iteration-boundary hook `_apply_multipool_spawn_pass` (model.py:389) on
candidate species *after* reaction generation. Surgery at stamping time
replaces the daughter before the candidates pass, so no `_scission` spawn
intent is ever queued for a chip event — no cancel hook needed.

### 4.3 `Reaction.polymer_chip_units`

New `cdef public int polymer_chip_units` (default 0) following the existing
pattern exactly: declared in `reaction.pxd`, kwarg in `reaction.py`,
**NOT serialized in `__reduce__`** (like `is_end_group_reaction` and
`polymer_flux_archetype`; unstamped arrivals demote at solver init).
Stamped in model.py during surgery:

```
a = round(chip_MW_g_mol / parent_monomer_MW_g_mol)   # parent.monomer_mw_g_mol
```

MW ratio, not subgraph matching — robust to radical sites, unsaturation and
H-transfer deltas; degrades gracefully and is monotone. Copolymers: per-pool
monomer MW is the number-average comonomer MW with an explicit config
override.

### 4.4 Tripwire diagnostic (the legitimate remnant of structural detection)

At stamping, when a scission-shaped **µ1-scaled** reaction's represented
piece is end-confined (piece ⊆ wing match **+ at most 1 repeat unit** —
"⊆ wing only" would miss cap+1-unit pieces, plausibly the most common
single-step end cuts, biasing the census low), log warn-once:

> "probable mis-scaled end-anchored cut; routed SCISSION_FRAGMENT pending
> detector item"

Structure is used for *diagnostics only*, never routing — no rate/statistics
inconsistency. The warn-once counter gives a measured census on real decks of
how much chemistry waits on the detector item (§9), which sets that item's
priority.

### 4.5 Liveness — partially live, not fully dormant

Precisely: the flag fires **today** for shapes with an END_MOD-stamped
product alongside the scission piece; with the discriminator swap those route
to DISCRETE_CHIP immediately (replacing their current UNRESOLVED dead-end) —
correctly, since flag-true reactions are µ0-scaled and the bundle matches.
What stays dormant is the single-step end-anchored cut population
(eliminations, retro-ene, concerted H-shift+scission), which is mis-scaled µ1
today and routes to SCISSION_FRAGMENT until the detector item lands. CHIP
archetypes appearing in logs are expected, not a bug. (The canonical unzip
case is unaffected by the gap: END_MOD fold-back means an end-activated chain
never exists as a separate core species, so its β-scission is never
family-generated — that chemistry lives in the `k_unzip` channel by design.)

## 5. Solver dispatch (`rmgpy/solver/polymer.pyx`)

New archetype branch in the section-5 dispatch, per direction
(`(rf, src), (rr, …)` — note the chip shape has no dst pool; the reverse leg
re-forms the parent chain from chip + folded parent, see below):

```
b0, b1, _, _ = _chain_bundle(src, y, V_poly, end_group=is_end_group[r_idx])
if b0 == 0: skip          # empty-pool sentinel, inherited
E_n = b1                  # µ1/µ0 (uniform) or µ2/µ1 (length-biased)
a   = chip_units[r_idx]

parent: dmu0 += 0
        dmu1 += -r * a
        dmu2 += -r * max(0.0, 2.0*a*E_n - a*a)
```

- **Closure-free — a genuine robustness win.** Only E[n] (first moment under
  the pick) is needed, never E[n²]: no `_safe_mu3_from_mu012` call, no
  `mu2_ok` branch, both E[n] forms are exact ratios of tracked moments.
- **One new guard (forward leg only):** `2a·E[n] − a²` goes negative when
  `a ≥ 2·E[n]` — impossible per-chain (n ≥ a always) but reachable in
  expectation for a degenerate pool whose mean length has decayed toward chip
  size. Clamp the µ2 decrement at ≥ 0 (skip the write); by that regime the
  pool is near-exhausted and the moment description is marginal anyway. The
  reverse term below is unconditionally positive — no clamp.
- **Exhaustion throttle (AMENDMENT 2026-06-10, adversarial-review finding):**
  DISCRETE_CHIP is the only stamped archetype that decouples the drained
  moment (µ1) from the rate-scaling moment (µ0 — which chip events
  deliberately never drain, per D4). Unthrottled, the µ0-scaled rate is
  constant in µ1: after the pool's units exhaust, µ1 runs linearly negative
  while chip moles keep being created — silent physical mass creation and a
  permanent cone exit (µ1 < 0 < µ0). Write-gating the µ1 leg alone would
  recreate mass duplication (chips produced at full rate, no drain), so the
  **rate itself throttles**, in the section-2 site scaling, for
  `archetype == DISCRETE_CHIP ∧ scaling == mu0 ∧ a > 0`:

  ```
  site = min(max(0, µ0), max(0, µ1)/a) / V_poly
  ```

  This is a **counting inequality, not a heuristic**: every chain with
  n ≥ a contributes at least a units to µ1, so #{chains long enough to
  donate} ≤ µ1/a always — the throttle is "eligible ends = min(all ends,
  an upper bound on donatable ends)". Throttled regime:
  dµ1/dt = −kf·(µ1/a)·a = −kf·µ1 exactly; unthrottled regime has
  µ0 < µ1/a so the bound dµ1/dt ≥ −kf·µ1 holds globally — µ1 decays at
  worst exponentially and never crosses zero. Because the throttle lives in
  the site scaling it **multiplies both directions** (the reverse
  re-attachment leg also throttles in exhaustion — acceptable: a
  near-zero-unit pool is a degenerate regime, and this is consistent with
  the existing src-density-both-directions convention; stated here so the
  asymmetry is not later read as a bug) and conservation stays exact at
  every instant (chip production and µ1 drain scale together; the
  Σµ1 + Σa·n invariant is unaffected). Interaction with the µ2 clamp: in
  the throttled regime E[n] = µ1/µ0 < a, so once E[n] < a/2 the clamp
  engages and µ2 writes stop while µ1 decays — the exhaustion test asserts
  the FULL cone (µ1 ≥ 0 AND µ0·µ2 ≥ µ1²), not just µ1 non-negativity.
  a = 0 chips are exempt (they drain nothing); µ1-scaled chips already
  self-limit exponentially and are untouched. A hard cutoff at µ1 < a·µ0
  was rejected: discontinuous, and it fires on the *mean*, freezing a
  polydisperse pool whose tail still has donatable chains, while the
  throttle drains it at the parcel-limited rate.
  **Accepted approximation (re-review finding D-4):** in the throttled
  regime the effective rate measure becomes donor-limited (chains with
  n ≥ a, conditional mean ≥ a) while the µ2 bundle still uses the
  population pick E[n] = µ1/µ0 < a — so µ2 is under-drained (the clamp
  engages) and PDI biases high in exhaustion. The bias is one-sided INTO
  the cone (µ2 too high over-satisfies µ0·µ2 ≥ µ1²) and touches no
  conserved quantity; the exact conditional mean is unrecoverable from
  three moments, so no better estimator exists at this closure order.
- The chip itself gains/loses moles through the **standard gas-species
  dn_dt path** — no special handling.
- **Reverse leg — NOT the sign-flip of the forward.** Extending a chain by a:
  Δ(n²) = (n+a)² − n² = +(2an + a²), so
  dmu1 += +rr·a, dmu2 += **+rr·(2a·E[n] + a²)** — plus a², not minus.
  E[n] evaluated on the same (single) pool — the re-formed chain extends a
  parent-statistics chain by a. Unlike SCISSION's r < 0 case, where the
  sign-flip was a declared approximation, the exact form here is free — use
  it.
- **Demotion:** stamped DISCRETE_CHIP with unresolved src pool → UNRESOLVED
  at `initialize_model`, extending the existing unresolvable-pool demotion
  path (same aggregate warning).
- **Mis-scaling caveat evaporates for this archetype** (unlike MIGRATION,
  where it's carried): routing and bundle both key on the same flag, so they
  cannot disagree with the implemented rate by construction.

**Conservation invariant (exact, closed system, chip reactions only):**

```
d/dt [ Σ_pools µ1  +  Σ_chips a_i · n_i ]  =  0
```

chip moles weighted by their stamped unit count — *not* raw chip moles.

## 6. Config

- `discrete_dp_threshold` per pool, **default 4** (monomer through trimer
  explicit), parsed from the input file. **Dormant under the fixed trimer
  proxy** — documented as such in `multi_pool_design.md`: the backbone branch
  is proxy-relative (D2), the scission branch keys on the flag (D3), and the
  backstop below requires proxy repeat-count > threshold.
- **Conditional DP backstop (D8), for future longer proxies:** when proxy
  repeat-count > threshold, a µ1-scaled cut producing a piece with literal
  DP < threshold → exact-a accounting (DISCRETE_CHIP-style bundle with the
  length-biased E[n], consistent with its µ1 scaling) — the representation
  genuinely resolves the cut position there, and exact-a beats uniform-cut
  statistics. This is a reclassification **toward** chip, never toward pool.
- No C40 heavy-atom ceiling (D7).

## 7. Cantera / balance audit

Chip-reaction balance is **shape-dependent** (§4.2 sub-shapes):

- **Shape (a)** fold-back is an *unmodified* `parent.copy`, so the reaction
  reads proxy → proxy + chip: **over-balanced by the chip mass**,
  `is_balanced()` fails. (a)-shapes join the unbalanced-proxy population and
  must be handled exactly like registered scissions (dropped from the Cantera
  export, logged/counted). This overbalance *is* the cap-mass drift A2
  accepts — the two statements agree by construction (cross-ref §10-A2).
- **Shape (b)** fold-back is a structurally modified END_MOD image, so
  proxy → modified-proxy + chip can be atom-exact: (b)-shapes can genuinely
  pass `is_balanced()` and **must survive** the unbalanced-proxy filter.

Audit anything keyed on "proxy scission ⇒ unbalanced":
- `rmgpy/cantera.py` `generate_cantera_data`'s filter: keeps balanced
  (b)-shapes, drops unbalanced (a)-shapes — verify both with tests (§8).
- Any balance validators/diagnostics that special-case proxy reactions.

## 8. Tests

Classification / stamping (`test/rmgpy/polymerTest.py`):
1. **Tolerance-pinning pair:** an engineered small molecule containing both
   wing substructures without the backbone (easiest with methyl-terminated
   wings) must gas-classify (`backbone_impostor`); a FEATURE-modified proxy
   image with +1 side group and −2 H must pass. Pins both sides of the
   tolerance.
2. **Discriminator regression (the tripwire for piece counting):** a
   µ1-scaled (flag-false) trimer cut whose represented piece is 1 unit must
   classify **SCISSION_FRAGMENT** — guards against anyone reintroducing
   piece-size routing.
3. Sub-shape (a) (macro daughter represented, no END_MOD product) →
   DISCRETE_CHIP: daughter replaced by CHIP fold-back, chip derived and
   registered discrete, `polymer_chip_units` stamped from MW ratio.
3b. **Sub-shape (b)** (SCISSION piece = chip + END_MOD fold-back — the live
   shape): chip demoted to discrete Molecule, END_MOD fold-back re-stamped
   CHIP, stored `is_end_group_reaction` still True after surgery (pins the
   flag-stability rider), archetype = DISCRETE_CHIP.
4. Surgery-infeasible shape → UNRESOLVED + warn-once (and never
   SCISSION_FRAGMENT).
5. Tripwire warn-once fires exactly once for a wing-confined µ1-scaled piece.
6. a = 0 chip (bare cap ejection) stamps cleanly.
7. Fold-back never-queue: after a chip-event stamping pass, no `_scission_*`
   Polymer is registered and no spawn intent exists.

Reaction field (`test/rmgpy/reactionTest.py`):
8. `polymer_chip_units` default/kwarg; absent from `__reduce__` round-trip.

Solver (`test/rmgpy/solver/solverPolymerTest.py`):
9. Constants-match-enum extended to DISCRETE_CHIP.
10. Monodisperse closed-form: pool of N chains, length L, chip a → dµ1 =
    −r·a, dµ2 = −r·(2aL − a²) for both picks (uniform pick E[n] = L =
    length-biased on monodisperse — also asserts the two picks agree there).
11. µ0-scaled vs µ1-scaled E[n] ratio test on a polydisperse pool (µ2/µ1 ≠
    µ1/µ0 — asserts the dispatch reads the right bundle per flag).
12. Clamp regime: degenerate pool with a ≥ 2·E[n] → µ2 write skipped, RHS
    finite.
12b. **Exhaustion trajectory (throttle amendment):** µ0-scaled chip reaction
    integrated past unit exhaustion — µ1 decays exponentially toward 0,
    never crosses zero, chip production rate → 0 in lockstep, and the FULL
    cone holds throughout (µ1 ≥ 0 AND µ0·µ2 ≥ µ1², not just µ1
    non-negativity — the µ2 channel is where a residual throttle/clamp
    inconsistency would hide). Also asserts Σµ1 + a·n_chip constant.
13. Conservation trajectory: Σ_pools µ1 + Σ a_i·n_i constant over a chip-only
    forward-Euler trajectory; cone preserved.
14. Reverse-leg test (rs.kb override trick, kf = 0): expects the exact
    extension form dµ1 = +rr·a, dµ2 = +rr·(2a·E[n] + a²) — NOT the forward
    sign-flip; also asserts no clamp on the reverse term.
15. Demotion: stamped DISCRETE_CHIP with src = −1 → UNRESOLVED.
16. Backstop dormancy: threshold ≥ proxy repeat-count → mid-cut still
    SCISSION_FRAGMENT (ties D8 to code).

Cantera (`test/rmgpy/canteraTest.py`):
17. A balanced (b)-shape chip reaction survives `generate_cantera_data`
    filtering.
17b. An over-balanced (a)-shape chip reaction is dropped and counted exactly
    like a registered unbalanced scission.

Back-compat: unstamped synthetic `Reaction()` objects continue to demote to
UNRESOLVED → legacy flux; all existing suites must stay green; EPDM
end-to-end unchanged (26 species / 28 reactions) — EPDM has no flag-true
scission shapes, so this is a true no-op gate for it.

## 9. Follow-up item (out of scope, scoped now)

**End-anchor detector:** extend `is_end_group_reaction` determination to
template/site anchoring (the reacted site on the proxy lies in the cap or
cap-adjacent unit, from template-labeled atoms before `clear_labeled_atoms`).
Chemistry-visible: flips affected reactions µ1 → µ0 scaling (rate × E[n],
2–4 orders of magnitude) — needs its own EPDM delta analysis,
reaction-by-reaction. **Acceptance criteria (written now per D6):**
(a) the §4.4 tripwire count on real decks drops to ~zero;
(b) the EPDM delta is explained reaction-by-reaction.
The seam is engineered so the detector item touches only detection and
re-scaling: routing, bundles, and tests downstream of the flag are already
correct — CHIP activates with zero solver changes.

## 10. Stated approximations + verified mechanism facts

Approximations (one line each goes to `multi_pool_design.md`):
- **A1:** Fold-back ignores the new chain end created by the cut
  (radical/unsaturated, a potential unzip initiator) — the same approximation
  END_MOD fold-back already makes everywhere; end-specific kinetics belong to
  k_unzip/END_MOD chemistry, not a spawned pool.
- **A2:** End-group mass drift on the chip (chip carries the old cap, parent
  image keeps both) is sub-monomer and consistent with existing µ-accounting
  (µ1 counts monomer units; cap masses are not tracked in moments). Shape
  (a)'s Cantera overbalance (§7) is this same drift made visible at the
  species-balance level — the two statements agree by construction.
- **A3:** Out of scope: true volatility (Tb vs reactor T) — mass-transfer
  layer, phase 2; the detector item (§9).
- **A4 (added with the exhaustion-throttle amendment):** the legacy
  UNRESOLVED path has the same latent structure for µ0-scaled reactions
  (rate ∝ µ0, drains µ1, µ0 untouched) — documented, NOT fixed. The
  UNRESOLVED contract is bit-exact legacy behavior (that is what keeps the
  pre-existing test population untouched), and the detector item (§9)
  migrates exactly these reactions to DISCRETE_CHIP, where the throttle
  protects them. The unprotected population shrinks by design; changing it
  now would silently break the back-compat property.

Verified against code during design (2026-06-10):
- **V1:** `is_end_group_reaction` = `any(p._reacted_class == END_MOD ...)`
  — product-stamp-based, no template anchoring (basis of D6).
- **V2:** Baseline proxy = head + 3 repeat units + tail
  (`_stitch_trimer`, polymer.py:763-791) — basis of D2/D3/D8.
- **V3 (to re-verify at implementation):** polymer+polymer coupling cannot
  deliver a ≥2×-proxy candidate to the wing_count ≥ 2 branch
  (`PolymerCrosslinkError` raised upstream). If false → add ceiling (§3.1).
- **V4:** spawn intents queue in `_apply_multipool_spawn_pass`
  (model.py:389), an iteration-boundary hook downstream of
  `make_new_reaction` → never-queue mechanism (§4.2) is sound.
- **V5:** the "polymer reactant with no polymer product → UNRESOLVED" guard
  already exists (polymer.py:1453-1460, tested) — tier 1 of the classifier
  interaction needed no new work.
- **V6 (probed 2026-06-10, plan-authoring):** propane/butane/pentane against
  CH3-capped PE classify SCISSION (wing_count 1), not 2-wing — the original
  item-#5 impostor example was wrong. The smallest 2-wing alkane (hexane)
  sits above the PE bound, and a head-wing+tail-wing directly-joined molecule
  is never below the 0.35 bound. The verified live impostor is **bibenzyl vs
  the PS proxy** (double-wing-matches today, classifies UNKNOWN; 14 heavy <
  bound 16) — the §8 test 1 pair is built on it, with probed pass-side
  candidates (27-heavy +CH3/−2H image; 19-heavy lost-phenyl image).
