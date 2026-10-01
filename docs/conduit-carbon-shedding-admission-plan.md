# Conduit carbon-shedding admission — design plan

Status: **DESIGN, not implemented.** Sequenced, adversarially-reviewed plan for admitting
genuine pool→small-carbon scission daughters into the polymer moment balance so a novolac
(phenol-formaldehyde) mechanism sheds mass by elementary chemistry, mass-conservingly.

Prerequisite already landed: the `feature_radical` tag-narrowing hygiene fix
(`annotate_feature_radical` / `_warn_once_refused` in `rmgpy/polymer_conduit.py` +
`rmgpy/polymer.py`), which un-blocks the self-referential admission gate G1 for
conduit-deferred rows. This plan builds on it. Everything below lands behind the existing
`polymerConduitAdmission` opt-in (default-off); decks that do not opt in stay byte-identical.

## Problem

A novolac deck's model sheds far less mass than measured (~11% modelled vs ~47% measured).
During RMG generation the moment/pool kernel refuses a large set of pool-coupled reactions via
a "feature-radical" census, zeroing their flux (stamp-but-keep). The refused set includes
genuine scission channels that would produce small carbon daughters (cresol/xylenol-sized
aromatics, and demethylation fragments). Because that flux is zeroed and the degradation is
left in the intact pool moments (the residual skips the row *before* the pool-moment section),
the elementary carbon fragments never carry flux and are never generated — the pool merely
accumulates radical isomers.

## Current contract (established by characterization; see the reviewed findings)

- **Flux-zeroing** is the whole-row residual skip in `rmgpy/solver/polymer.pyx` (keyed on
  `reaction_refused[r_idx]`), which returns before the pool-moment section, so scission mass
  stays in-pool by *omission* of a debit — not by an active recycle.
- **OFF (default):** the reaction's flux is zeroed; the daughter species is still generated
  ("stamp-but-keep") but carries no flux; scission mass remains in the pool moments.
- **ON (`polymerConduitAdmission=True`):** does not, by itself, admit these rows. The
  admission *policy* must change — the direction and mass-debit machinery below are required.

## Classification of the refused set (audited, r105: 634 rows)

- **12 genuine pool→small-carbon scission** rows: 4 aromatic volatiles (C7 cresyl, C8
  xylenol-yl ×2, C9) + 8 C1 ejections (4 CH2(S) carbene + 4 CH3 methyl = demethylation).
- **622 bookkeeping** rows (correctly deferred): 288 pool+pool disproportionation, 205
  isomerization relabel, 129 bath-gas H-abstraction — all daughters pool-sized.
- The two populations are gap-separated by daughter carbon count ({1,7,8,9} then pool-sized),
  so a structural predicate can separate them — **but** each scission row is written in the
  *recombination* direction (`small_frag + large_frag <=> pool`, pool on the product side).
  The mass-shed is the **reverse** (pool → small volatile + large residue).

## Why this sequence (the adversarial review is the spec)

Three findings gate everything:

1. Narrowing `feature_radical` (landed) un-blocks G1, but the rows then re-deny downstream
   (G2 `classifier-not-admissible` / G3 `gas-product-count`) — admission is not reachable
   until the reverse (scission) orientation is made runnable. **Flip-rewrite is prerequisite
   #1, sequenced first.**
2. There is **no** existing per-row same-pool pool-debit path. The arch-7 conduit is
   credit-only (pool-*growing*); arch-6 `VOLATILE_EJECTION` is a cross-pool migration mirror
   (`rmgpy/solver/polymer.pyx`, requires `src != dst` — two pools), and the novolac deck has
   one pool. **A same-pool mass-conserving debit is prerequisite #2**, and it must use the
   length-biased bundle math (b0/b1/b2 + an `a`-unit shift), **not** a Dirac (F, F·u, F·u²)
   form.
3. The predicate must be defined against the moment-pool **distribution** (Mn/Mw), not a
   discrete proxy, or it overfits.

## Phase 0 — premise re-probes (read-only; each can invert the plan)

- **0.1 Reversibility of the 12:** confirm per-row that the scission direction has reverse
  kinetics (mechanism is ~99% reversible / `keepIrreversible=False`, but any
  irreversible-as-written row's scission direction can never run and must **not** be admitted —
  a fluxless daughter is a no-op at best, a fabricated source at worst).
- **0.2 Same-pool debit primitive:** read the `VOLATILE_EJECTION` same-pool leg and the
  `DISCRETE_CHIP` branch in `rmgpy/solver/polymer.pyx`; pick the primitive that conserves
  mu0/mu1/mu2 for "pool consumes one chain → gas fragment (`a` units) + residue that stays in
  the same pool." Do not hand-roll a Dirac debit.
- **0.3 Reference-state size for the predicate:** express "C_fragment < C_pool" against the
  distribution via the monomer-equivalent / reference-state MW (with the existing slack),
  not a discrete proxy count; confirm the 12-vs-622 separation survives.
- **0.4 k_scission overlap:** the phenomenological random-scission channel (`k_scission`) is
  active on the novolac pool; it conserves mu1 (sheds no gas today) but books the same C–C
  homolysis population dynamics. Decide the reconciliation rule now (zero `k_scission` on a
  pool that gains elementary scission, or prove disjoint) — it is a Phase-2 acceptance gate.

## Phase 1 — prerequisite #1: reverse-kinetics flip-rewrite (direction machinery)

- Where: the direction gate in `rmgpy/polymer_conduit.py` (today a reversible-needs-reverse
  row lands on `direction-requires-flip-rewrite`, deferred).
- What: for a reversible recombination-written row whose scission orientation is otherwise
  admissible, produce a runnable scission-direction reaction — either synthesize reverse
  Arrhenius from `kf` and `Keq(T)` from the participants' thermo, or emit the flipped reaction
  and use RMG's existing reverse-rate path. Deterministic; inert when admission is off.
- Hardest sub-problem: `Keq` needs pool thermo + reference-state MW, and the pool is a
  distribution, so "Keq of a distribution scission" must be defined consistently with the
  moment closure (probes 0.1/0.3 feed this).
- Acceptance: the scission rate matches `kf`/`Keq(T)` within tolerance across the reactor T
  range; zero effect with admission off; irreversible-as-written rows remain inadmissible
  (not fabricated).

## Phase 2 — prerequisite #2: same-pool mass-conserving debit archetype

- Where: `rmgpy/solver/polymer.pyx` — a new/extended archetype for "pool → gas fragment(s) +
  same-pool residue" (distinct from cross-pool `VOLATILE_EJECTION` and the credit-only conduit).
- What: on the scission direction, debit the source pool with the length-biased bundle
  (b0/b1/b2) and an `a` = ejected-fragment backbone-units shift (mirror the `VOLATILE_EJECTION`
  same-pool leg); credit the gas fragment by ordinary net-rate stoichiometry; the residue
  remains in the same pool.
- **Mass invariant (acceptance gate):** for each admitted row, pool moment mu1 debit ==
  mass(gas fragment); residue mass == mu1_before − mass(gas); the mu2 shift is the exact
  quadratic `(n+sa)² − n²`. No Dirac `F·u²` shortcut.
- **Double-count reconciliation (enforce, not warn):** apply the 0.4 rule; upgrade the
  double-count tripwire from warn-only to a real guard for this path.

## Phase 3 — admission predicate/policy (depends on 1+2)

- Where: `evaluate_conduit_admission` in `rmgpy/polymer_conduit.py` — a scission-admission
  classifier alongside the existing shape buckets, evaluated after the (now-narrowed)
  feature_radical gate and the direction gate.
- Predicate: the row is a pool scission whose non-pool side has ≥1 carbon fragment strictly
  smaller than the pool's reference-state size (0.3 definition), with the complementary
  fragment(s) accounting for the balance (mass-conserving split), passing the existing gas_mw
  + balance gates.
- **Cutoff policy decision (explicit, chemistry-justified, not tuned):** start
  aromatic-volatiles-only (C7–C9 stable molecular volatiles; the 4 clean rows). Treat C1
  carbene/methyl ejection (the 8 demethylation rows) as a separate, later toggle — CH2(S)/CH3
  are transient radicals, not stable volatiles, and are the more arguable half. Neither cutoff
  is chosen to hit a mass-loss number.

## Phase 4 — acceptance gates + regression (every phase gates on these)

- **Per-row moment-debit mass audit (hard gate):** automated assertion — for every admitted
  scission row, pool moment debit == gas + residue mass within tolerance (solver invariant +
  unit tests with synthetic pools).
- **OFF byte-identical regression:** decks without the opt-in are byte-identical.
- **622-defer regression:** the bookkeeping rows must still defer — no new admissions among them.
- Determinism + reset/epoch behavior unchanged.

## Phase 5 — full validation run (schedule separately)

Native novolac generation (r105/r106-class deck), dual-stream captured, measuring: resin
mass-loss % (vs ~47% measured / ~11% current), admitted-scission-row count (expect ~4
aromatic, or ~12 with C1), per-row moment-debit audit pass, and the 622-defer regression.
Pair with the OFF byte-identical smoke.

## Guardrail (carried into every phase)

If admitting scission still undershoots ~47%, that is a **finding about non-refused physics**
(k_scission/k_unzip magnitudes, reference-state MW, missing families — truly-small fragments
already pass unrefused), **not** license to widen admission. Do not touch admission-widening
code to chase the number; do not reshape the predicate/cutoff to make the TGA number come out
right.

## Adversarial-review checkpoints (each phase, before merge)

- Phase 1: attack the Keq-of-a-distribution derivation and irreversible-row handling.
- Phase 2: attack mu0/mu1/mu2 conservation on unequal fragment splits; verify against the
  `VOLATILE_EJECTION` same-pool math, not a Dirac form; verify k_scission reconciliation is
  enforced.
- Phase 3: attack predicate overfit against the distribution pool and on a non-novolac deck.

Each merge is user-gated.
