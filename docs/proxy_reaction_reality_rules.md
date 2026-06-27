# Proxy Reaction Reality Rules

**Status:** design spec (2026-06-27). Governs which proxy-touching reactions are
physically meaningful for the moment-represented polymer, in which direction, and how a
kept reaction's flux maps onto polymer moments. Companion to
[`multi_pool_design.md`](multi_pool_design.md) (esp. §4 spawning, §5.1 flux archetypes,
§5.4 item-18 refuse) and [`polymer_moments_format.md`](polymer_moments_format.md).

This document has **two layers**, which must be kept conceptually separate:

- **Layer 1 — Reality rules:** *which* proxy reactions to keep, and *in which direction*.
- **Layer 2 — Representation:** *how* a kept reaction's flux maps onto polymer moments
  (feature/daughter pools + the moment kernels of §5.1).

Conflating "which reactions" with "what moment recipe" is how this subsystem gets muddled;
the rules below decide membership/orientation, the kernels decide bookkeeping.

---

## 0. Governing principle — what the proxy *is*

The trimer proxy (e.g. polystyrene `PS(2)`, C25H28) is **not a species whose amount the
solver tracks** — its concentration is *pinned* (~1.0). It is a **lumped local stand-in
for an interior segment of a long chain** whose real state lives in the moments
(μ0 = moles of chains, μ1 = moles of repeat units, μ2 = second moment / spread). A proxy
reaction's only legitimate job is to **transduce a local site-event into a change in
(μ0, μ1, μ2) + released small molecules, weighted by how many such sites the real chain
actually has** (`site_scaling`). The proxy being *consumed or formed as a species* is
meaningless — what matters is the reaction's effect on the moments.

From this, one discriminator governs everything:

> **Suppress only what is *structurally invalid* for a real long chain. Never suppress
> what is merely *slow* — that is RMG's rate-based flux analysis to decide.**

---

## Layer 1 — Reality rules (membership and orientation)

### R1 — Real site, weighted by its moment
The reactive center must be a site that exists on the real chain: an **interior**
backbone C–C / C–H, a **pendant** group (e.g. the phenyl on each PS unit), or a **genuine
chain end / cap**. Caps are *allowed* — the proxy's `[CH3]`/`[H]` termini stand in for the
chain's real end groups, which do react (end-group initiation, end scission, unzip). Site
type sets the abundance weighting: interior sites scale with **μ1** (abundant), end/cap
sites scale with **μ0** (rare; ratio ≈ μ1/μ0 = DPn). So cap reactivity is principally
allowed and comes out ~DPn× rarer through `site_scaling`, automatically.

### R2 — Local, never whole-proxy
On a real chain the proxy's two ends are ~DPn units apart and never meet. Any reaction that
needs **both ends / the finite span** of the proxy to interact — transannular ring closure
across the proxy, retro-Diels–Alder spanning the molecule, intra-disproportionation between
the two caps, end-to-end H-transfer — is a **finite-size artifact** and is **structurally
invalid**. Suppress it. (Reactions local to a *single* unit — including that unit's own
pendant ring — are fine; locality, not ring-involvement, is the test.)

### R3 — Canonicalize to the physical direction; orient-and-keep
A reaction family's canonical direction is **arbitrary** and is often the reverse of the
physics: RMG may enumerate a real chain scission `proxy → A + B` written backwards as
`A + B → proxy` (recombination), because the C25 fragments happen to recombine into exactly
the proxy structure (mapped back by isomorphism). **Never filter on the surface direction.**
Canonicalize each proxy-touching reaction to its physical direction first; if the physical
direction is proxy-*consuming* (degradation / local modification) and it satisfies R1–R2 and
R4, **keep and orient it**. Suppress only reactions that remain genuine proxy-*assembly*
*after* orientation (no degrading reading exists) — these are the only Rule-3 suppressions.

> **Why this rule exists (failure it prevents):** a syntactic test ("proxy appears only as a
> product ⇒ spurious") silently deletes real scission written in reverse, and tempts a
> fallback to a single crude lumped rate that discards the families' bond-resolved chemistry.
> Orient first; preserve the physics.

### R4 — Represent, don't refuse (see Layer 2)
Every structurally-valid reaction's effect must be expressible as **flux between polymer
populations + small molecules, conserving backbone repeat units** (Σ μ1 invariant over the
discrete subset; + monomer for unzip). When the product is a population the current basis
cannot represent — e.g. a *same-length interior backbone radical* from mid-chain
H-abstraction — **do not refuse the flux (the item-18 stopgap, §5.4) and do not leak the
whole chain to the gas (which fabricates backbone mass).** Instead **spawn the lumped
feature pool** that represents it and route flux into it (Layer 2). Refusal is an
implementation limitation, not a reality criterion.

### R5 — Relevance is RMG's job
Do **not** hand-suppress chemistry on "this is probably slow" grounds. Give RMG every
structurally-valid reaction; its thermo/kinetics estimation pays the correct energetic
penalty (e.g. aromaticity loss for a pendant-phenyl dearomatization → correctly tiny flux),
and its move-to-core flux analysis leaves negligible reactions in the edge. Hand-suppression
substitutes intuition for the rate analysis that is RMG's entire purpose. The only cost of
keeping slow chemistry is **model size**, which is handled orthogonally by a polymer-scoped
constraints block (`generatePolymerConstraints`), not by the reality filter.

### The double-count contract
Once oriented family scission/unzip reactions carry mechanistic degradation (R3), they must
not be **double-counted** with the explicit phenomenological `k_scission` / `k_unzip`
channels. Contract: **families carry mechanistic degradation; the lumped channels are
explicit user overrides for processes the families do not capture; no event is counted
twice.** The existing double-count tripwire and the §5.1 input-hygiene hazard (a reaction
with BOTH a family archetype AND a nonzero phenomenological channel) enforce this.

### Suppression scope (small, by construction)
After R1–R5, hand-suppression is limited to: **(R2)** finite-size / whole-proxy events, and
**(R3)** proxy-*assembly* reactions with no degrading reading after orientation. Everything
else is kept and handed to Layer 2 + the rate analysis.

---

## Layer 2 — Representation: feature/daughter pools + moment kernels

### Feature pools (the item-20 successor to item-18 refuse)
A **feature pool** is a lumped sub-population of chains characterized by a local structural
feature (a radical, a modified unit) — its own `[μ0, μ1, μ2]`, **born at zero
concentration**, with flux routed *in* (formation) and later *out* (consumption). This
generalizes the existing daughter-pool spawning (§4.2–4.3 `discover_repeat_motif` /
`similarity_merge` / `feature_monomers`) from "structurally distinct chain types" to
"chains bearing a reactive feature", and **supersedes the item-18 refuse machinery (§5.4)**:
mid-chain H-abstraction routes `parent → feature` (conserving μ1, since the chain did not
change length), and the feature's downstream β-scission routes `feature → product pools +
small molecules`.

### Feature-pool equivalencing key (the tractability hinge)
Feature identity is the **canonicalized local feature-type signature only** — *not* exact
chain length and *not* position along the chain. Example (locked): *"a chain bearing one
mid-chain secondary benzylic radical"* is **one** feature pool regardless of where on the
chain the radical sits. Chain length is carried by the pool's moments; radical position is
handled statistically by the kernel (below). This is the same lump-by-motif discipline that
makes the proxy itself finite, applied one level up — without it, feature pools proliferate
into a pool explosion (the finite-size problem re-emerging). Reuse `similarity_merge`
(strict isomorphism + edit-distance ≤ 1 on the feature motif) for the key.

### Lifecycle and solver registration
Spawn at zero concentration → register in the **solver pool config** (`species_to_pool_indices`
/ `polymer_pools`) so MIGRATION/SCISSION archetypes can *resolve* their source/target pools.
The 2 "UNRESOLVED" reactions observed in runs failed precisely because a spawned daughter
pool was never registered — **the same registration gap feature pools must close.** R4 and
the daughter-pool-registration fix are therefore *one* piece of work.

### QSSA: eliminate vs. accumulate
A radical feature pool is usually a quasi-steady-state intermediate (formed and consumed
fast, low concentration). The eliminating-vs-accumulating classifier (§5.4) decides whether
to **QSSA-eliminate** it algebraically (flux-in ≈ flux-out, no state variable) or
**dynamically integrate** it as a real population. "Move flux in, then out" is compatible
with both; default to QSSA-eliminate for fast-β-scission radicals, accumulate only when the
feature is long-lived.

### The moment kernel — already exists; the bug is unstamped reactions
The product **size distribution** of a kept reaction comes from the §5.1 archetype kernels,
**not** from the proxy's fixed C25 size. A naïve reading ("scission of C25 always gives
C12 + C13") is wrong; the real chain has a length distribution and the break is at a
distributed position, so products are a distribution. The `SCISSION_FRAGMENT` archetype
already encodes this (parent `−r·µ2/(2µ1)`, daughter `+r·µ2/(2µ1)`, "conserves µ1 exactly";
`E[k]=µ2/µ1` is the length-biased mean) — i.e. the random-scission / Ziff–McGrady moment
apportionment is **implemented**. The reaction supplies the **rate and bond-class**; the
kernel + parent/feature moments supply the **product size distribution**.

**So the real defect is upstream:** stamping fires only when a `Polymer` is a *reactant*
(`make_new_reaction`), but the proxy is overwhelmingly a *product* (recombination written
backwards), so 858/860 reactions reach the solver **unstamped** and fall back to crude
"legacy mu1-only" flux instead of the kernel. R3 orientation is the fix: orient first, then
stamp the physical direction with the correct archetype so it uses the §5.1 kernel.

---

## Worked examples

| Reaction | Verdict | Why / recipe |
|---|---|---|
| `proxy + •CH3 → proxy•(interior) + CH4` (mid-chain H-abstraction) | **keep → feature pool** | Same-length radical; R4 spawns "mid-chain secondary benzylic radical" feature pool, flux `parent→feature`, μ1 conserved. NOT refused. |
| `proxy•(interior) → shorter chain + unsaturated end` (β-scission) | **keep, kernel** | SCISSION_FRAGMENT archetype; product moments from §5.1 position kernel, not fixed C25 split. |
| `A + B → proxy` (recombination as enumerated) | **orient → keep as `proxy → A+B`** | R3: physical direction is scission; orient and stamp it. Do not suppress by surface direction. |
| `A + B → proxy` with no degrading reading | **suppress** | R3: genuine proxy-assembly; structurally meaningless for a pinned stand-in. |
| ring closure spanning the whole C25 proxy | **suppress** | R2: finite-size artifact; impossible on a real long chain. |
| pendant-phenyl dearomatization (one unit's own ring) | **keep** | R5: local + representable; energetically unfavorable → RMG rates leave it in the edge. Cost is size, not correctness. |

---

## Code map (as of 2026-06-27 — verify before editing)

- Archetype stamp: `rmgpy/polymer.py:~1682` `stamp_polymer_flux_archetype` →
  `classify_reaction_flux_archetype:~1523`; called in `rmgpy/rmg/model.py:~663` / `~705`,
  **gated on a `Polymer` reactant** (`if polymer_reactants:`) — the R3 gap.
- Solver read + warnings: `rmgpy/solver/polymer.pyx:~910` (read
  `polymer_flux_archetype`), `~1102-1133` ("arrived without a stamp" / "could not resolve
  their solver pool(s)").
- Refuse machinery (item 18, to be superseded by Layer 2): `rmgpy/solver/polymer.pyx:~2024-2032`;
  census at `~886-898`.
- Pool registration: `rmgpy/rmg/polymer_input.py` `get_gas_mask` / pool config;
  `rmgpy/solver/polymer.pyx:~790-834` `_apply_pool_phase_overrides` (`species_to_pool_indices`).
- Double-count tripwire: `rmgpy/solver/polymer.pyx:~1138`.

## Known omissions / open items
- **Cyclization / crosslinking** (indane/indene fused-ring formation) is **logged as an
  accepted omission** for the depolymerization moment model; revisit only if a deck needs it.
- **Feature-pool proliferation guard**: the equivalencing key bounds it, but a cap +
  ranking/eviction (cf. §4.3 note) may be needed if a real deck spawns many features.
- **`is_polymer_proxy` sticky-contamination** (family.py blanket-stamping; see
  `multi_pool_design.md` and the `get_gas_mask` fix that keys phase off `is_moment_dummy`)
  must not leak the polymer-scoped constraints or pool identity onto ordinary coreactants.

## Cross-references
- [`multi_pool_design.md`](multi_pool_design.md): §4.2–4.3 (motif discovery / similarity_merge),
  §5.1 (flux archetypes + kernels), §5.4 (item-18 refuse — superseded by Layer 2), §5.3
  (prospective mask).
- [`polymer_moments_format.md`](polymer_moments_format.md): moment basis, unzip invariant,
  `monomer_routing`, mass-transfer.
