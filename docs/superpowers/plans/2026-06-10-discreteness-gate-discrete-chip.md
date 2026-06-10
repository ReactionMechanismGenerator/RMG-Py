# Discreteness Gate + DISCRETE_CHIP Archetype Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement the discreteness gate (backbone impostor rejection) and the DISCRETE_CHIP flux archetype per docs/superpowers/specs/2026-06-10-discreteness-gate-discrete-chip-design.md.

**Architecture:** Generation-time chip product surgery (`surge_chip_products` + `stamp_polymer_flux_archetype` in `rmgpy/polymer.py`, called from both `make_new_reaction` stamping branches) rewrites flag-true scission shapes into [discrete chip Molecule, CHIP-stamped fold-back Polymer] and stamps `Reaction.polymer_chip_units`; the solver resolves a `reaction_chip_units` array at `initialize_model` and dispatches a closure-free per-direction DISCRETE_CHIP branch in the RHS. The heavy-atom discreteness gate in `classify_structure` rejects backbone impostors before they ever become pool images. Spec: `docs/superpowers/specs/2026-06-10-discreteness-gate-discrete-chip-design.md` (READ IT FIRST — it has the decision table D1–D8 and the math).

**Tech Stack:** Python 3.9 (conda rmg_env), Cython (.pyx/.pxd via make), pytest.

---

**Working directory:** `/home/alon/Code/RMG-Py` (all paths below relative to it).

**Conventions (environment facts — do not skip):**
- Run tests with `~/anaconda3/envs/rmg_env/bin/python -m pytest <file>::<test> -v` (or `-q` for suites).
- After ANY edit to a `.pyx` or `.pxd` file, run `make` from the repo root BEFORE running tests, or tests silently run stale compiled code.
- cdef class fields MUST be declared in the `.pxd` (`Reaction` is a cdef class; `HybridPolymerSystem` in `polymer.pyx` is a plain Python class — its new attributes need NO pxd entry). `Reaction.__reduce__` deliberately omits the polymer fields; unstamped arrivals demote at solver init.
- Commit style: `topic: imperative subject` ≤60 chars, body only if needed, trailer `Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>`. Stage explicit paths, never `git add -A`. One logical change per commit. Never push.
- Warn-once registries: tests touching `_warn_unresolved_archetype` must clear `rmgpy.polymer._flux_archetype_warned` in setup; tests touching the new tripwire must clear `rmgpy.polymer._chip_tripwire_warned`.
- Synthetic solver tests: an unstamped `Reaction()` auto-demotes to UNRESOLVED → legacy flux; stamp `rxn.polymer_flux_archetype` AND `rxn.polymer_chip_units` explicitly to exercise the new paths. The `rs.kb[0]` override drives reverse rates without thermo (set `rs.kf[0] = 0.0`, inject product/chip moles into `y` directly).

**Key chemistry (from the spec, restated for the implementer):**

A DISCRETE_CHIP event ejects a discrete chip of `a` repeat units (`a = round(chip_MW / monomer_MW)`, stamped at generation time; `a = 0` is legal) from a chain of the SOURCE pool; the complement folds back into the same pool. Per direction, with `E[n] = b1` from `_chain_bundle(src, …, end_group=flag)` (µ1/µ0 uniform when µ0-scaled, µ2/µ1 length-biased when µ1-scaled — the rate scaling IS the sampling measure, D5):

| Leg | dµ0 | dµ1 | dµ2 |
|---|---|---|---|
| forward (rf, ejection) | 0 | −rf·a | −rf·max(0, 2a·E[n] − a²) (clamp forward only) |
| reverse (rr, re-attachment) | 0 | +rr·a | +rr·(2a·E[n] + a²) — EXACT extension form, NOT the sign-flip; never clamped |

Closure-free: only E[n] is needed (exact moment ratio), never µ3. The chip species itself flows through the standard gas `dn_dt` path. Conservation invariant (closed system, chip reactions only): `d/dt [ Σ_pools µ1 + Σ_chips a_i·n_i ] = 0` (chip moles weighted by stamped unit count).

**Empirically verified fixtures (probed against current code, 2026-06-10 — trust these):**
- PS proxy (`monomer='[CH2][CH]c1ccccc1'` or the labeled adj used in tests, end_groups `['[CH3]','[H]']`): baseline proxy = 25 heavy atoms, MW 328.5 g/mol, monomer MW 104.15 g/mol; gate tol = round(0.35·25) = 9, bound = 16.
- Bibenzyl `c1ccccc1CCc1ccccc1` (14 heavy) genuinely matches BOTH PS wings today (`wing_count == 2`, classifies UNKNOWN/`unclassified_intact_backbone`) → the canonical backbone impostor; 14 < 16 → gate rejects it.
- Modified-proxy pass candidates: `CC(CC(C)=C(CC(C)c1ccccc1)c1ccccc1)c1ccccc1` (27 heavy, +CH3 side group −2H) and `CC(CCCC(C)c1ccccc1)c1ccccc1` (19 heavy, lost center phenyl) — both wing_count 2, both ≥ 16 → pass.
- The live chip shape (b): `_handshake_structures([stitch(head_wing, CH3*2), radicalized_proxy], [ps])` yields `[PS_scission_tail (SCISSION), PS (END_MOD)]`, `is_end_group_reaction(products) is True`; the fragment is `CCC(C)C1=CC=CC=C1`, MW 134.2 → `a = round(134.2/104.15) = 1`.
- Scission Polymers do NOT currently retain their source fragment (`_source_molecule` does not exist) — Task 4 adds it.

---

## File Structure

| File | Change |
|---|---|
| `rmgpy/polymer.py` | Discreteness gate in `classify_structure`; `PolymerClass.CHIP`; `PolymerFluxArchetype.DISCRETE_CHIP = 5`; CHIP classifier branch + reworded unsurged fallback; `_source_molecule` stamping in `create_reacted_copy`; tripwire (`_chip_tripwire_warned` + `_warn_probable_end_cut`); `surge_chip_products`; `stamp_polymer_flux_archetype`; `discrete_dp_threshold` field (+ `copy()`) |
| `rmgpy/reaction.pxd` | `cdef public int polymer_chip_units` |
| `rmgpy/reaction.py` | `polymer_chip_units=0` kwarg + assignment (NOT in `__reduce__`) |
| `rmgpy/rmg/model.py` | Both `make_new_reaction` stamping branches call `stamp_polymer_flux_archetype` |
| `rmgpy/solver/polymer.pyx` | `FLUX_DISCRETE_CHIP = 5`; `reaction_chip_units` int32 array; demotion extension (src == −1); DISCRETE_CHIP dispatch branch |
| `rmgpy/rmg/input.py` | `discrete_dp_threshold` kwarg on `polymer()` (dormant) |
| `docs/multi_pool_design.md` | §5.1 DISCRETE_CHIP row, gate/liveness/knob notes, A1–A3, detector follow-up |
| Tests | `test/rmgpy/polymerTest.py`, `test/rmgpy/reactionTest.py`, `test/rmgpy/solver/solverPolymerTest.py`, `test/rmgpy/canteraTest.py` |

---

### Task 1: Backbone impostor gate in `classify_structure` (spec §3.1)

**Files:**
- Modify: `rmgpy/polymer.py` — `classify_structure`, inside the `wing_count >= 2` branch (~line 1623), after the CROSSLINK check (~1627) and before the BASELINE isomorphism check (~1632)
- Test: `test/rmgpy/polymerTest.py` — class `TestPolymerClassification` (~line 1735; fixture provides `self.p`, the PS polymer with CH3/H caps; `Molecule`, `Species`, `polymer` are imported at the top of the file)

- [ ] **Step 1: Write the failing tests**

Add to `TestPolymerClassification`, right after `test_branch_gas_too_few_atoms` (~line 1807):

```python
    def test_discreteness_gate_rejects_backbone_impostor(self):
        """
        Discreteness gate (spec 2026-06-10 §3.1): a wing_count >= 2 candidate
        far smaller than the baseline proxy is a backbone impostor, not a
        chain image. Bibenzyl (14 heavy atoms) genuinely matches both PS wing
        subgraphs today but sits below the one-sided bound
        proxy_heavy - round(0.35*proxy_heavy) = 25 - 9 = 16 -> GAS.
        """
        impostor = Molecule(smiles="c1ccccc1CCc1ccccc1")  # bibenzyl, 14 heavy
        p_class, details = polymer.classify_structure(
            Species(molecule=[impostor]), self.p)
        assert p_class == polymer.PolymerClass.GAS
        assert details["reason"] == "backbone_impostor"

    def test_discreteness_gate_keeps_modified_proxy_images(self):
        """
        Tolerance pin, pass side (one-sided on purpose): legitimately LARGER
        images (+1 side group, -2 H => 27 heavy) and modestly SMALLER images
        (lost center phenyl => 19 heavy >= bound 16) must NOT gas-classify as
        impostors.
        """
        bigger = Molecule(smiles="CC(CC(C)=C(CC(C)c1ccccc1)c1ccccc1)c1ccccc1")
        p_class, details = polymer.classify_structure(
            Species(molecule=[bigger]), self.p)
        assert details["reason"] != "backbone_impostor"
        assert p_class != polymer.PolymerClass.GAS

        smaller = Molecule(smiles="CC(CCCC(C)c1ccccc1)c1ccccc1")  # 19 heavy
        p_class, details = polymer.classify_structure(
            Species(molecule=[smaller]), self.p)
        assert details["reason"] != "backbone_impostor"
        assert p_class != polymer.PolymerClass.GAS
```

- [ ] **Step 2: Run to verify failure**

Run: `~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/polymerTest.py -q -k "discreteness_gate"`
Expected: `test_discreteness_gate_rejects_backbone_impostor` FAILS (bibenzyl currently classifies UNKNOWN, reason `unclassified_intact_backbone`); the keeps-test PASSES already (it is a regression pin for the pass side — that is fine).

- [ ] **Step 3: Implement the gate**

In `rmgpy/polymer.py` `classify_structure`, the `wing_count >= 2` branch currently reads (lines ~1622-1633):

```python
    # BRANCH A: INTACT BACKBONE (2 or more wings) ---
    if wing_count >= 2:
        base_details.update(wing_match_details)

        # 1. Crosslinking Check (>2 wings means bi-molecular polymer combination)
        if wing_count > 2:
            return PolymerClass.CROSSLINK, {**base_details, "reason": "more_than_two_wings_found"}

        # 2. Baseline Check (Is it the exact unreacted proxy?)
        # Evaluates: [X-O]-O-[O-Y]
        if original_polymer.baseline_proxy.is_isomorphic(species):
            return PolymerClass.BASELINE, {**base_details, "reason": "unreacted_proxy"}
```

Insert between the crosslink check and the baseline check:

```python
        # 1.5 Discreteness gate: backbone impostor rejection (spec 2026-06-10
        # docs/superpowers/specs/2026-06-10-discreteness-gate-discrete-chip-design.md
        # §3.1). A genuine 2-wing candidate is an image of the baseline proxy
        # REPRESENTATION; a small molecule that merely contains both wing
        # subgraphs (e.g. bibenzyl against a PS proxy) sits far below it in
        # heavy atoms. One-sided on purpose: legitimate images can be larger
        # (FEATURE side groups) or modestly smaller (H loss, side-group
        # elimination). Heavy atoms, not MW: H-insensitive, and a lost side
        # group is a known heavy-atom delta.
        # No upper ceiling (spec §10-V3, verified 2026-06-10): polymer+polymer
        # coupling cannot deliver a >=2x-proxy candidate here -- a coupling
        # product carries >2 wings and returns CROSSLINK above, and
        # create_reacted_copy raises PolymerCrosslinkError for it upstream
        # (polymer.py create_reacted_copy crosslink guard), so make_new_reaction
        # discards the whole reaction before classification of coupled shapes.
        proxy_heavy = sum(
            1 for a in original_polymer.baseline_proxy.molecule[0].atoms
            if not a.is_hydrogen())
        cand_heavy = sum(
            1 for a in species.molecule[0].atoms if not a.is_hydrogen())
        if cand_heavy < proxy_heavy - round(0.35 * proxy_heavy):
            return PolymerClass.GAS, {**base_details, "reason": "backbone_impostor"}
```

- [ ] **Step 4: Verify spec §10-V3 (the no-ceiling contingency) and record it**

This is a verification step, not new code. Confirm both facts the comment above asserts:
1. In `classify_structure`, the `wing_count > 2 → CROSSLINK` return (step 1) precedes the gate insertion point — a coupling candidate with >2 detectable wings can never reach the gate. (Read the diff you just made.)
2. `Polymer.create_reacted_copy` (rmgpy/polymer.py, ~lines 816-828) raises `PolymerCrosslinkError` for CROSSLINK-classified products, and `make_new_reaction` (rmgpy/rmg/model.py, ~lines 636-641 and 691-696) catches it and discards the reaction.

If you find a reachable path that delivers a ≥2×-proxy candidate to the `wing_count == 2` gate (e.g. a coupling product whose extra wings are destroyed by the coupling bond), STOP and add a ceiling with the same proxy-relative form (`cand_heavy > 2*proxy_heavy + tol → GAS`) plus a test; otherwise leave the comment as the record. (Analysis note: a coupling product with exactly 2 surviving wings is structurally a longer single backbone — indistinguishable from a legitimate image and bounded near 2× proxy, which is why no ceiling is required.)

- [ ] **Step 5: Run to verify pass + no classification regressions**

Run: `~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/polymerTest.py -q`
Expected: all passed (+2 new). The existing classification fixtures are all full proxy images (≥19 heavy vs bound 16) and PE candidates below 6 heavy never reach `wing_count >= 2`, so no existing test should trip the gate. If one does, inspect whether it asserted a small two-wing molecule staying polymer-classified (that was the bug this gate fixes — update it citing the spec); anything else, STOP.

- [ ] **Step 6: Run the multipool suite (classify_structure callers)**

Run: `~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/polymerMultiPoolTest.py -q`
Expected: all passed.

- [ ] **Step 7: Commit**

```bash
git add rmgpy/polymer.py test/rmgpy/polymerTest.py
git commit -m "polymer: gate backbone impostors out of pool classification"
```

---

### Task 2: `PolymerClass.CHIP`, `DISCRETE_CHIP` enum, classifier CHIP branch (spec §3.3, §4.1)

**Files:**
- Modify: `rmgpy/polymer.py` — `PolymerClass` (~line 1363), `PolymerFluxArchetype` (~line 1398), `classify_reaction_flux_archetype` (~line 1426)
- Test: `test/rmgpy/polymerTest.py` — class `TestPolymer`, after `test_classify_reaction_flux_archetype` (~line 833; fixture provides `self.polymer_1`)

- [ ] **Step 1: Write the failing tests**

Add to `TestPolymer` right after `test_classify_reaction_flux_archetype`:

```python
    def test_classify_chip_product_returns_discrete_chip(self):
        """
        A CHIP-stamped product polymer (left by surge_chip_products' fold-back)
        short-circuits to DISCRETE_CHIP BEFORE the SCISSION branch. Order
        matters (spec 2026-06-10 §4.1): after the (b)-surgery the product list
        has no END_MOD member, so the SCISSION branch's internal
        is_end_group_reaction(products) recompute would return False and
        misroute; the CHIP check must win even alongside a SCISSION member.
        """
        from rmgpy.polymer import PolymerFluxArchetype, classify_reaction_flux_archetype

        p = self.polymer_1
        chip_gas = Molecule(smiles="CC")
        fold = p.copy()
        fold._reacted_class = PolymerClass.CHIP
        assert (classify_reaction_flux_archetype([p], [chip_gas, fold])
                == PolymerFluxArchetype.DISCRETE_CHIP)

        # Defensive order pin: CHIP beats a co-present SCISSION stamp.
        sc = p.copy()
        sc.label = f"{p.label}_scission_head"
        sc._reacted_class = PolymerClass.SCISSION
        assert (classify_reaction_flux_archetype([p], [sc, fold])
                == PolymerFluxArchetype.DISCRETE_CHIP)

    def test_flag_false_one_unit_piece_routes_scission_fragment(self):
        """
        Discriminator regression (spec 2026-06-10 test 2, decision D3): a
        mu1-scaled (flag-false) cut whose represented piece is ONE repeat unit
        must still route SCISSION_FRAGMENT. On a 3-unit proxy a u1-u2 cut (the
        image of EVERY interior backbone bond) yields 1- and 2-unit pieces --
        literal piece size is a representation artifact and must never be a
        routing input. Guards against reintroducing piece-size routing.
        """
        import rmgpy.polymer as polymer_mod
        polymer_mod._flux_archetype_warned.clear()
        from rmgpy.polymer import PolymerFluxArchetype, classify_reaction_flux_archetype

        p = self.polymer_1
        gas = Molecule(smiles="CC")
        piece = p.copy()
        piece.label = f"{p.label}_scission_tail"
        piece._reacted_class = PolymerClass.SCISSION
        # Represented piece: cap + ONE styrene unit (10 heavy atoms).
        piece._source_molecule = Molecule(smiles="CCC(C)c1ccccc1")
        assert (classify_reaction_flux_archetype([p], [piece, gas])
                == PolymerFluxArchetype.SCISSION_FRAGMENT)
```

- [ ] **Step 2: Run to verify failure**

Run: `~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/polymerTest.py -q -k "chip_product_returns or one_unit_piece"`
Expected: `test_classify_chip_product_returns_discrete_chip` FAILS with `AttributeError: CHIP` (enum member missing). `test_flag_false_one_unit_piece_routes_scission_fragment` PASSES already (it pins current flag-false routing; `_source_molecule` is just an ad-hoc attribute until Task 4 makes it real). If the second test fails, STOP — the classifier has a real bug.

- [ ] **Step 3: Implement**

(a) In `PolymerClass` (rmgpy/polymer.py ~1378-1380), extend the chain-breaking block:

```python
    # --- Chain Breaking/Linking States ---
    SCISSION = 'SCISSION'    # Only 1 wing found; chain broke (e.g., [X-O]-W)
    CROSSLINK = 'CROSSLINK'  # >2 wings found; bi-molecular polymer recombination
    CHIP = 'CHIP'            # Fold-back parent copy left by DISCRETE_CHIP product
                             # surgery (surge_chip_products, spec 2026-06-10 §4.2).
                             # Never produced by classify_structure itself.
```

(b) In `PolymerFluxArchetype` (~1405-1409), add after `UNRESOLVED = 4`:

```python
    UNRESOLVED = 4         # ambiguous/unstamped: solver applies legacy mu1 flux
    DISCRETE_CHIP = 5      # end-anchored cut ejects a stamped a-unit discrete
                           # chip; complement folds back to the SAME pool
                           # (spec 2026-06-10). Mirror: solver FLUX_DISCRETE_CHIP.
```

(c) In `classify_reaction_flux_archetype` (~1426), insert the CHIP branch directly after the NONE early-return (`if not reactant_pools and not product_polymers: return ... NONE`, ~line 1437) and BEFORE the SCISSION `any(...)` check (~1439):

```python
    if any(getattr(p, '_reacted_class', None) == PolymerClass.CHIP
           for p in product_polymers):
        # Chip product surgery (surge_chip_products) already rewrote this
        # product list to [discrete chip, CHIP-stamped fold-back]. This check
        # MUST precede the SCISSION branch: after the (b)-surgery there is no
        # END_MOD member left, so is_end_group_reaction(products) would
        # recompute False and misroute. The solver reads the STORED Reaction
        # flag; nothing downstream may recompute it from product stamps.
        return PolymerFluxArchetype.DISCRETE_CHIP
```

(d) Reword the existing SCISSION ∧ end-group branch (~1441-1449) as the unsurged/invariant-violation fallback — replace its comment and warn reason (the return value is unchanged):

```python
        if is_end_group_reaction(products):
            # Unsurged end-initiated scission: chip product surgery
            # (surge_chip_products, spec 2026-06-10 §4.2) was either not
            # attempted or infeasible for this shape. Surged shapes never
            # reach here (the CHIP branch above short-circuits). UNRESOLVED,
            # never SCISSION_FRAGMENT: uniform-cut statistics near a chain
            # end are wrong AND the chip mass would go unaccounted.
            _warn_unresolved_archetype(
                "unsurged end-initiated scission",
                tuple(sorted(p.label for p in product_polymers)))
            return PolymerFluxArchetype.UNRESOLVED
        return PolymerFluxArchetype.SCISSION_FRAGMENT
```

Note: the existing `test_classify_reaction_flux_archetype` asserts the UNRESOLVED return and counts warn-once entries by registry size, not by reason string — the reason rewording does not break it. Verify in step 4.

- [ ] **Step 4: Run to verify pass**

Run: `~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/polymerTest.py -q`
Expected: all passed (+2 new, all pre-existing classifier tests green).

- [ ] **Step 5: Commit**

```bash
git add rmgpy/polymer.py test/rmgpy/polymerTest.py
git commit -m "polymer: add CHIP class and DISCRETE_CHIP archetype"
```

---

### Task 3: `Reaction.polymer_chip_units` field (spec §4.3, test 8)

**Files:**
- Modify: `rmgpy/reaction.pxd` (after line 67, `cdef public int polymer_flux_archetype`)
- Modify: `rmgpy/reaction.py` (`__init__` kwargs ~line 130, assignment ~line 162)
- Test: `test/rmgpy/reactionTest.py` — class `TestReaction`, after `test_polymer_flux_archetype_default_and_kwarg` (~line 796)

- [ ] **Step 1: Write the failing test**

Add to `TestReaction` right after `test_polymer_flux_archetype_default_and_kwarg`:

```python
    def test_polymer_chip_units_default_and_kwarg(self):
        """
        Reaction carries ``polymer_chip_units`` (default 0): the repeat-unit
        count of the discrete chip ejected by a DISCRETE_CHIP reaction,
        stamped by chip product surgery at generation time. Like
        is_end_group_reaction / polymer_flux_archetype it is deliberately NOT
        serialized in __reduce__ -- unstamped arrivals demote at solver init.
        """
        import pickle

        assert Reaction().polymer_chip_units == 0
        assert Reaction(polymer_chip_units=2).polymer_chip_units == 2

        round_trip = pickle.loads(pickle.dumps(Reaction(polymer_chip_units=2)))
        assert round_trip.polymer_chip_units == 0
```

- [ ] **Step 2: Run to verify failure**

Run: `~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/reactionTest.py -q -k test_polymer_chip_units_default_and_kwarg`
Expected: FAIL — `AttributeError: ... has no attribute 'polymer_chip_units'` (cdef class: undeclared attributes don't exist).

- [ ] **Step 3: Implement**

In `rmgpy/reaction.pxd`, after line 67 (`cdef public int polymer_flux_archetype`):

```pxd
    cdef public int polymer_chip_units
```

In `rmgpy/reaction.py` `__init__`, extend the kwargs (lines 129-131):

```python
                 is_end_group_reaction=False,
                 polymer_flux_archetype=0,
                 polymer_chip_units=0,
                 ):
```

and after `self.polymer_flux_archetype = polymer_flux_archetype` (line 162):

```python
        # Repeat-unit count of the discrete chip for DISCRETE_CHIP reactions
        # (a = round(chip_MW / monomer_MW), a == 0 legal), stamped by chip
        # product surgery in rmgpy.polymer.stamp_polymer_flux_archetype via
        # make_new_reaction. Like the two fields above it is NOT serialized in
        # __reduce__; the solver demotes unstamped arrivals at initialize_model.
        self.polymer_chip_units = polymer_chip_units
```

Do NOT modify `__reduce__` (line 209) — the omission is the point.

- [ ] **Step 4: Rebuild Cython (pxd changed)**

Run: `make`
Expected: compiles without error.

- [ ] **Step 5: Run to verify pass**

Run: `~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/reactionTest.py -q -k "polymer_chip_units or polymer_flux_archetype or is_end_group_reaction_default"`
Expected: 3 passed.

- [ ] **Step 6: Commit**

```bash
git add rmgpy/reaction.pxd rmgpy/reaction.py test/rmgpy/reactionTest.py
git commit -m "reaction: declare polymer_chip_units field"
```

---

### Task 4: Chip product surgery + tripwire in `rmgpy/polymer.py` (spec §4.2, §4.4; tests 3, 3b, 4, 5, 6)

**Files:**
- Modify: `rmgpy/polymer.py` — `create_reacted_copy` (~line 854), new module functions after `classify_reaction_flux_archetype` (~line 1470), tripwire inside `classify_reaction_flux_archetype`
- Test: `test/rmgpy/polymerTest.py` — `TestPolymer` (tripwire) and `TestHandshakeStructures` (~line 2672; provides `self.ps`, `self._handshake`, and the module provides `_methyl_radical_adj` (~2639) and `radicalize_head_end_group` (~2347))

- [ ] **Step 1: Write the failing tests**

(a) Add to `TestHandshakeStructures` after `test_handshake_products_classify_flux_archetype` (~line 2778):

```python
    def test_surge_chip_sub_shape_b_live_end_mod_fold_back(self):
        """
        Spec test 3b -- the only flag-true shape live today: products =
        [SCISSION piece, END_MOD fold-back]. Surgery demotes the chip back to
        a discrete Molecule (undoing its handshake conversion), re-stamps the
        END_MOD fold-back CHIP, and returns a = round(134.2/104.15) = 1.
        Flag-stability rider: the recompute over surged products flips to
        False (END_MOD member gone) -- which is exactly why nothing downstream
        may recompute the flag from product stamps.
        """
        from rmgpy.polymer import (PolymerFluxArchetype, is_end_group_reaction,
                                   classify_reaction_flux_archetype,
                                   surge_chip_products)

        head_wing = self.ps._stitch_wing("head")
        methyl_star2 = Molecule().from_adjacency_list(_methyl_radical_adj("*2"))
        frag = polymer.stitch_molecules_by_labeled_atoms(head_wing, methyl_star2)
        assert frag is not None
        end_mod = self.ps.baseline_proxy.molecule[0].copy(deep=True)
        radicalize_head_end_group(self.ps, end_mod)

        products = [frag.copy(deep=True), end_mod]
        self._handshake(products, [self.ps])
        assert products[0]._reacted_class == PolymerClass.SCISSION
        assert products[1]._reacted_class == PolymerClass.END_MOD
        assert is_end_group_reaction(products) is True  # the stored flag's value

        a = surge_chip_products(products, self.ps)

        assert a == 1
        assert isinstance(products[0], Molecule)          # chip demoted
        assert not isinstance(products[0], Polymer)
        assert products[0].get_formula() == frag.get_formula()
        assert products[1]._reacted_class == PolymerClass.CHIP
        # Recompute now flips -- pins the no-recompute rule.
        assert is_end_group_reaction(products) is False
        assert (classify_reaction_flux_archetype([self.ps], products)
                == PolymerFluxArchetype.DISCRETE_CHIP)

    def test_surge_chip_sub_shape_a_macro_daughter(self):
        """
        Spec test 3 -- sub-shape (a) (dormant today, live when the end-anchor
        detector lands): the SCISSION-stamped Polymer is the MACRO daughter
        and the chip is the single discrete co-product. Surgery replaces the
        daughter with parent.copy(deep=True) stamped CHIP; the chip stays
        as-is; a stamps from the chip's MW ratio.
        """
        from rmgpy.polymer import (PolymerFluxArchetype,
                                   classify_reaction_flux_archetype,
                                   surge_chip_products)

        head_wing = self.ps._stitch_wing("head")
        methyl_star2 = Molecule().from_adjacency_list(_methyl_radical_adj("*2"))
        frag = polymer.stitch_molecules_by_labeled_atoms(head_wing, methyl_star2)
        prods = [frag]
        self._handshake(prods, [self.ps])
        daughter = prods[0]
        assert isinstance(daughter, Polymer)
        assert daughter._reacted_class == PolymerClass.SCISSION

        chip = Molecule(smiles="C=Cc1ccccc1")   # styrene, 104.15 g/mol -> a = 1
        products = [daughter, chip]
        a = surge_chip_products(products, self.ps)

        assert a == 1
        fold = products[0]
        assert isinstance(fold, Polymer)
        assert fold.label == self.ps.label                # PARENT pool fold-back
        assert fold._reacted_class == PolymerClass.CHIP
        assert products[1] is chip                        # chip untouched, discrete
        assert (classify_reaction_flux_archetype([self.ps], products)
                == PolymerFluxArchetype.DISCRETE_CHIP)

    def test_surge_chip_infeasible_stamps_unresolved_never_scission(self):
        """
        Spec test 4: surgery-infeasible flag-true scission shapes stamp
        UNRESOLVED + warn-once via stamp_polymer_flux_archetype -- NEVER
        SCISSION_FRAGMENT (uniform-cut statistics near an end + unaccounted
        chip mass). Two infeasible shapes: (b) without a demotable source
        molecule, and (a) without a discrete co-product.
        """
        import rmgpy.polymer as polymer_mod
        polymer_mod._flux_archetype_warned.clear()
        from rmgpy.reaction import Reaction
        from rmgpy.polymer import PolymerFluxArchetype, stamp_polymer_flux_archetype

        # Infeasible (b): SCISSION chip with no _source_molecule.
        sc = self.ps.copy()
        sc.label = "PS_scission_tail"
        sc._reacted_class = PolymerClass.SCISSION
        end = self.ps.copy()
        end._reacted_class = PolymerClass.END_MOD
        rxn = Reaction(reactants=[self.ps], products=[sc, end],
                       is_end_group_reaction=True)
        stamp_polymer_flux_archetype(rxn, [self.ps], [self.ps])
        assert rxn.polymer_flux_archetype == int(PolymerFluxArchetype.UNRESOLVED)
        assert rxn.polymer_chip_units == 0

        # Infeasible (a): flag-true scission shape, no discrete co-product.
        sc2 = self.ps.copy()
        sc2.label = "PS_scission_head"
        sc2._reacted_class = PolymerClass.SCISSION
        rxn2 = Reaction(reactants=[self.ps], products=[sc2],
                        is_end_group_reaction=True)
        stamp_polymer_flux_archetype(rxn2, [self.ps], [self.ps])
        assert rxn2.polymer_flux_archetype == int(PolymerFluxArchetype.UNRESOLVED)

        # Warn-once: repeating an already-warned shape adds no registry entry.
        n = len(polymer_mod._flux_archetype_warned)
        stamp_polymer_flux_archetype(rxn2, [self.ps], [self.ps])
        assert len(polymer_mod._flux_archetype_warned) == n

    def test_surge_chip_a_zero_bare_cap_ejection(self):
        """
        Spec test 6: a = 0 chips are legal (bare end-cap ejection, e.g. CH3
        loss): surgery succeeds and returns 0 (NOT None) -- the archetype
        fires with zero mu1/mu2 drain, net pool effect ~ SAME_POOL.
        """
        from rmgpy.polymer import (PolymerFluxArchetype,
                                   classify_reaction_flux_archetype,
                                   surge_chip_products)

        sc = self.ps.copy()
        sc.label = "PS_scission_tail"
        sc._reacted_class = PolymerClass.SCISSION
        sc._source_molecule = Molecule(smiles="C")   # CH4 cap image, 16 g/mol
        end = self.ps.copy()
        end._reacted_class = PolymerClass.END_MOD
        products = [sc, end]

        a = surge_chip_products(products, self.ps)

        assert a == 0
        assert a is not None                         # 0 != infeasible
        assert isinstance(products[0], Molecule)
        assert products[1]._reacted_class == PolymerClass.CHIP
        assert (classify_reaction_flux_archetype([self.ps], products)
                == PolymerFluxArchetype.DISCRETE_CHIP)
```

(b) Add to `TestPolymer` after `test_flag_false_one_unit_piece_routes_scission_fragment` (Task 2):

```python
    def test_tripwire_warns_once_for_wing_confined_mu1_piece(self, caplog):
        """
        Spec test 5 (§4.4): a mu1-scaled (flag-false) scission whose
        represented piece is end-confined (piece <= wing + at most 1 repeat
        unit by heavy atoms) logs the probable-mis-scaled-end-cut warning
        exactly once. Diagnostics only -- routing stays SCISSION_FRAGMENT.
        The census sets the priority of the end-anchor detector follow-up.
        """
        import logging as _logging
        import rmgpy.polymer as polymer_mod
        polymer_mod._flux_archetype_warned.clear()
        polymer_mod._chip_tripwire_warned.clear()
        from rmgpy.polymer import PolymerFluxArchetype, classify_reaction_flux_archetype

        p = self.polymer_1
        gas = Molecule(smiles="CC")
        piece = p.copy()
        piece.label = f"{p.label}_scission_tail"
        piece._reacted_class = PolymerClass.SCISSION
        # cap (1 heavy) + 1 unit (8 heavy) + stitch CH3 = 10 heavy
        #   <= max_cap_heavy(1) + 2*monomer_heavy(16) = 17 -> end-confined.
        piece._source_molecule = Molecule(smiles="CCC(C)c1ccccc1")

        with caplog.at_level(_logging.WARNING):
            arch = classify_reaction_flux_archetype([p], [piece, gas])
        assert arch == PolymerFluxArchetype.SCISSION_FRAGMENT  # never routes
        hits = [r for r in caplog.records
                if "probable mis-scaled end-anchored cut" in r.getMessage()]
        assert len(hits) == 1

        with caplog.at_level(_logging.WARNING):
            classify_reaction_flux_archetype([p], [piece, gas])  # repeat
        hits = [r for r in caplog.records
                if "probable mis-scaled end-anchored cut" in r.getMessage()]
        assert len(hits) == 1                                   # warn-once
```

- [ ] **Step 2: Run to verify failure**

Run: `~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/polymerTest.py -q -k "surge_chip or tripwire"`
Expected: 5 FAIL — `ImportError: cannot import name 'surge_chip_products'` (and `_chip_tripwire_warned` AttributeError for the tripwire test).

- [ ] **Step 3: Implement — source-molecule stamp in `create_reacted_copy`**

In `rmgpy/polymer.py` `create_reacted_copy`, the verdict stamp currently reads (~lines 851-854):

```python
        # Stamp the classification verdict so the polymer handshake can flag
        # END_MOD reactions for chain-end (mu0) scaling in the solver. Read by
        # is_end_group_reaction(products); a transient generation-time marker.
        new_poly._reacted_class = klass
```

Add directly after it:

```python
        # Keep the sanitized reacted fragment so chip product surgery
        # (surge_chip_products, spec 2026-06-10 §4.2) can demote a SCISSION
        # chip back to a discrete Molecule and size chips by MW. Transient
        # generation-time marker like _reacted_class (deliberately not carried
        # by Polymer.copy()).
        new_poly._source_molecule = probe
```

(`probe` is the labels-cleared deep copy already computed at the top of `create_reacted_copy`.)

- [ ] **Step 4: Implement — tripwire registry + diagnostic**

In `rmgpy/polymer.py`, directly after `_warn_unresolved_archetype` (~line 1423), add:

```python
_chip_tripwire_warned = set()


def _warn_probable_end_cut(detail) -> None:
    """
    Diagnostics-only census (spec 2026-06-10 §4.4): warn once per distinct
    piece when a mu1-scaled scission's represented piece is end-confined.
    Never affects routing. The accumulated count on real decks measures how
    much chemistry waits on the end-anchor detector follow-up item.
    """
    if detail not in _chip_tripwire_warned:
        _chip_tripwire_warned.add(detail)
        logging.warning(
            "Polymer scission piece %s is end-confined (wing + <=1 repeat "
            "unit) but the reaction is mu1-scaled: probable mis-scaled "
            "end-anchored cut; routed SCISSION_FRAGMENT pending the "
            "end-anchor detector item.", detail)
```

Then, in `classify_reaction_flux_archetype`, replace the bare `return PolymerFluxArchetype.SCISSION_FRAGMENT` (end of the SCISSION branch, after the Task 2 rewording) with:

```python
        # Tripwire diagnostic (spec 2026-06-10 §4.4): structure is used for
        # DIAGNOSTICS only, never routing. "wing + at most 1 repeat unit"
        # (not "wing only") so cap+1-unit pieces -- plausibly the most common
        # single-step end cuts -- are counted. Heavy-atom bound:
        # max(cap heavy) + 2 * monomer heavy (wing = cap + 1 unit, plus 1).
        parent = next((r for r in reactants if isinstance(r, Polymer)), None)
        piece = next(
            (p for p in product_polymers
             if getattr(p, '_reacted_class', None) == PolymerClass.SCISSION),
            None)
        src_mol = getattr(piece, '_source_molecule', None) if piece is not None else None
        if parent is not None and src_mol is not None:
            piece_heavy = sum(1 for a in src_mol.atoms if not a.is_hydrogen())
            mon_heavy = sum(1 for a in parent.monomer.atoms if not a.is_hydrogen())
            cap_heavy = max(
                sum(1 for a in eg.atoms if not a.is_hydrogen())
                for eg in parent.end_groups)
            if piece_heavy <= cap_heavy + 2 * mon_heavy:
                _warn_probable_end_cut(piece.label)
        return PolymerFluxArchetype.SCISSION_FRAGMENT
```

(Polymers without `_source_molecule` — e.g. fabricated test objects or pre-Task-4 pickles — skip the diagnostic silently.)

- [ ] **Step 5: Implement — `surge_chip_products` + `stamp_polymer_flux_archetype`**

Add to `rmgpy/polymer.py` directly after `classify_reaction_flux_archetype` (~line 1470, before `MatchMapping = ...`):

```python
def surge_chip_products(products, parent: 'Polymer') -> Optional[int]:
    """
    Chip product surgery (spec 2026-06-10 §4.2): rewrite a flag-true
    (end-group-scaled) scission-shaped product list IN PLACE into the
    canonical DISCRETE_CHIP end state [discrete chip Molecule, CHIP-stamped
    fold-back Polymer], and return the chip's repeat-unit count
    ``a = round(chip_MW / parent.monomer_mw_g_mol)``. ``a == 0`` is legal
    (bare cap ejection) and distinct from the infeasible return ``None``.

    Sub-shapes (chip identification = the smaller piece; the cut POSITION is
    already known from the reaction's stored is_end_group_reaction flag):

    (b) END_MOD fold-back present -- the only flag-true shape live today:
        the SCISSION-stamped Polymer is the chip. Demote it back to its
        sanitized source Molecule (undoing the handshake conversion) and
        re-stamp the existing END_MOD fold-back as CHIP. Applying (a) here
        would replace the chip with a second fold-back -- losing the chip
        and double-folding the parent.
    (a) no END_MOD product (dormant until the end-anchor detector item):
        the SCISSION-stamped Polymer is the MACRO daughter; the chip is the
        single non-Polymer co-product (already discrete -- the handshake left
        it a Molecule). Replace the daughter with ``parent.copy(deep=True)``
        stamped ``PolymerClass.CHIP``.

    Returns ``None`` WITHOUT modifying ``products`` when the shape is not a
    feasible chip shape (no/ambiguous scission piece, chip unrepresentable or
    not the smaller piece, missing source molecule). The caller
    (stamp_polymer_flux_archetype) then stamps UNRESOLVED -- never
    SCISSION_FRAGMENT.
    """
    product_polymers = [p for p in products if isinstance(p, Polymer)]
    scissions = [p for p in product_polymers
                 if getattr(p, '_reacted_class', None) == PolymerClass.SCISSION]
    end_mods = [p for p in product_polymers
                if getattr(p, '_reacted_class', None) == PolymerClass.END_MOD]
    if len(scissions) != 1:
        return None  # not a chip shape, or ambiguous piece identification
    daughter = scissions[0]
    proxy_mw = parent.baseline_proxy.molecule[0].get_molecular_weight()

    if end_mods:
        # --- sub-shape (b): SCISSION piece = chip, END_MOD = fold-back ---
        if len(end_mods) != 1:
            return None
        chip_src = getattr(daughter, '_source_molecule', None)
        if chip_src is None:
            return None  # cannot demote: handshake source unavailable
        if chip_src.get_molecular_weight() >= proxy_mw:
            return None  # "chip" is not the smaller piece
        chip_mol = chip_src.copy(deep=True)
        chip_mol.clear_labeled_atoms()
        for i, p in enumerate(products):
            if p is daughter:
                products[i] = chip_mol
                break
        end_mods[0]._reacted_class = PolymerClass.CHIP
        chip_mw_g = chip_mol.get_molecular_weight() * 1000.0
        return int(round(chip_mw_g / parent.monomer_mw_g_mol))

    # --- sub-shape (a): SCISSION piece = macro daughter, chip discrete ---
    chips = [p for p in products
             if not isinstance(p, Polymer) and isinstance(p, Molecule)]
    if len(chips) != 1:
        return None  # chip absent or ambiguous
    chip_mol = chips[0]
    daughter_src = getattr(daughter, '_source_molecule', None)
    ref_mw = (daughter_src.get_molecular_weight()
              if daughter_src is not None else proxy_mw)
    if chip_mol.get_molecular_weight() >= ref_mw:
        return None  # the discrete co-product is not the smaller piece
    fold = parent.copy(deep=True)
    fold._reacted_class = PolymerClass.CHIP
    for i, p in enumerate(products):
        if p is daughter:
            products[i] = fold
            break
    chip_mw_g = chip_mol.get_molecular_weight() * 1000.0
    return int(round(chip_mw_g / parent.monomer_mw_g_mol))


def stamp_polymer_flux_archetype(forward, reactants, polymer_reactants) -> None:
    """
    Stamp ``forward.polymer_flux_archetype`` (and ``polymer_chip_units``)
    AFTER the handshake and AFTER ``forward.is_end_group_reaction`` is stored.
    Called from both make_new_reaction stamping branches (rmgpy/rmg/model.py).

    Flag-true shapes run chip product surgery FIRST so the classifier sees
    the surged product list (its CHIP branch precedes the SCISSION
    recompute). An infeasible flag-true scission shape stamps UNRESOLVED +
    warn-once -- never SCISSION_FRAGMENT (spec 2026-06-10 §4.2).
    """
    chip_a = None
    if forward.is_end_group_reaction and len(polymer_reactants) == 1:
        chip_a = surge_chip_products(forward.products, polymer_reactants[0])
    if chip_a is not None:
        forward.polymer_chip_units = chip_a
    elif forward.is_end_group_reaction and any(
            getattr(p, '_reacted_class', None) == PolymerClass.SCISSION
            for p in forward.products if isinstance(p, Polymer)):
        _warn_unresolved_archetype(
            "infeasible chip surgery",
            tuple(sorted(getattr(p, 'label', '?') for p in forward.products
                         if isinstance(p, Polymer))))
        forward.polymer_flux_archetype = int(PolymerFluxArchetype.UNRESOLVED)
        return
    forward.polymer_flux_archetype = int(
        classify_reaction_flux_archetype(reactants, forward.products))
```

(`Optional` is already imported in `rmgpy/polymer.py`'s typing imports; `logging` likewise.)

- [ ] **Step 6: Run to verify pass**

Run: `~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/polymerTest.py -q -k "surge_chip or tripwire"`
Expected: 5 passed.

- [ ] **Step 7: Run the whole polymer suite**

Run: `~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/polymerTest.py test/rmgpy/polymerMultiPoolTest.py -q`
Expected: all passed. (The `_source_molecule` stamp is additive; if any handshake/registration test fails, inspect whether it asserted an exact attribute set on Polymers — none should.)

- [ ] **Step 8: Commit**

```bash
git add rmgpy/polymer.py test/rmgpy/polymerTest.py
git commit -m "polymer: chip product surgery and end-cut tripwire"
```

---

### Task 5: Wire the surgery into `make_new_reaction` (spec §4.2 never-queue, §4.5; test 7)

**Files:**
- Modify: `rmgpy/rmg/model.py` — import line 47 and both stamping branches (~lines 644-646 and 685-687)
- Test: `test/rmgpy/polymerTest.py` — class `TestMakeNewReactionPolymer` (~line 3141; provides `self.model` (a `CoreEdgeReactionModel`) and `self.ps` registered via `_register_polymer`)

- [ ] **Step 1: Write the failing test**

Add to `TestMakeNewReactionPolymer` after `test_make_new_reaction_pairs_regenerated` (~line 3243):

```python
    def _make_chip_template_reaction(self):
        """
        Proxy -> [cap+unit fragment, END_MOD image]: handshakes into the live
        chip shape (b) (SCISSION piece + END_MOD fold-back), per the probed
        recipe in the 2026-06-10 spec work. Fragment MW 134.2 -> a = 1.
        """
        from rmgpy.data.kinetics.family import TemplateReaction
        from rmgpy.kinetics import Arrhenius

        proxy_mol = self.ps.baseline_proxy.molecule[0].copy(deep=True)
        head_wing = self.ps._stitch_wing("head")
        methyl_star2 = Molecule().from_adjacency_list(_methyl_radical_adj("*2"))
        frag = polymer.stitch_molecules_by_labeled_atoms(head_wing, methyl_star2)
        end_mod = self.ps.baseline_proxy.molecule[0].copy(deep=True)
        radicalize_head_end_group(self.ps, end_mod)
        return TemplateReaction(
            reactants=[proxy_mol],
            products=[frag, end_mod],
            family='R_Recombination',
            is_forward=True,
            kinetics=Arrhenius(A=(1e13, 's^-1'), n=0.0, Ea=(50.0, 'kcal/mol')),
            pairs=[(proxy_mol, frag), (proxy_mol, end_mod)],
        )

    def test_chip_event_stamps_discrete_chip_and_never_queues_pool(self):
        """
        Spec test 7 (+3b's flag-survival rider at the model level): a chip
        event through make_new_reaction stamps DISCRETE_CHIP with
        polymer_chip_units = 1, keeps the STORED is_end_group_reaction True
        (the surgery removed the END_MOD member, so a recompute would flip
        it -- nothing recomputes), registers NO _scission_* daughter, and the
        iteration-boundary spawn pass finds nothing to spawn (never-queue:
        surgery replaced the daughter before the candidates pass).
        """
        from rmgpy.polymer import PolymerFluxArchetype

        rxn = self._make_chip_template_reaction()
        result_rxn, is_new = self.model.make_new_reaction(
            rxn, check_existing=False, generate_thermo=False,
            generate_kinetics=False,
        )
        assert result_rxn is not None

        # Stamps: archetype, chip units, stored-flag survival.
        assert result_rxn.is_end_group_reaction is True
        assert (result_rxn.polymer_flux_archetype
                == int(PolymerFluxArchetype.DISCRETE_CHIP))
        assert result_rxn.polymer_chip_units == 1

        # Products: a discrete chip + the PS fold-back; no scission daughter.
        labels = [getattr(p, 'label', '') for p in result_rxn.products]
        assert not any('_scission' in lbl for lbl in labels)
        assert any(isinstance(p, Polymer) and p.label == 'PS'
                   for p in result_rxn.products)

        # Never-queue: no _scission_* Polymer registered...
        assert not any('_scission' in s.label for s in self.model.new_species_list)
        # ...and the spawn pass has nothing to drain (no daughter pools appear).
        self.model._apply_multipool_spawn_pass(self.model.new_species_list)
        assert not any('_scission' in s.label for s in self.model.new_species_list)
        pools = [s for s in self.model.new_species_list if isinstance(s, Polymer)]
        assert all(p.label == 'PS' for p in pools)
```

- [ ] **Step 2: Run to verify failure**

Run: `~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/polymerTest.py -q -k test_chip_event_stamps_discrete_chip`
Expected: FAIL — today the products handshake to [SCISSION, END_MOD] and the classifier's unsurged fallback stamps UNRESOLVED (4), so the `polymer_flux_archetype == 5` assertion fails (and `polymer_chip_units == 0`). The failure message should show archetype 4 — that confirms the pre-wiring routing.

- [ ] **Step 3: Implement the wiring**

In `rmgpy/rmg/model.py` line 47, change the import from:

```python
from rmgpy.polymer import Polymer, PolymerCrosslinkError, PolymerFluxArchetype, is_end_group_reaction, classify_reaction_flux_archetype
```

to:

```python
from rmgpy.polymer import Polymer, PolymerCrosslinkError, PolymerFluxArchetype, is_end_group_reaction, stamp_polymer_flux_archetype
```

Branch 1 (`forward.family and forward.is_forward`, lines ~644-646) currently:

```python
                forward.is_end_group_reaction = is_end_group_reaction(forward.products)
                forward.polymer_flux_archetype = int(
                    classify_reaction_flux_archetype(reactants, forward.products))
```

becomes:

```python
                forward.is_end_group_reaction = is_end_group_reaction(forward.products)
                # Chip product surgery + archetype stamping (spec 2026-06-10
                # §4.2): MUST run after the stored is_end_group_reaction flag
                # and before make_new_species registers the products, so the
                # SCISSION daughter is replaced before any spawn-candidate or
                # registration pass can see it (never-queue).
                stamp_polymer_flux_archetype(forward, reactants, polymer_reactants)
```

Branch 2 (lines ~685-687) currently:

```python
                    forward.is_end_group_reaction = is_end_group_reaction(forward.products)
                    forward.polymer_flux_archetype = int(
                        classify_reaction_flux_archetype(reactants, forward.products))
```

becomes:

```python
                    forward.is_end_group_reaction = is_end_group_reaction(forward.products)
                    # Same surgery + stamping as the is_forward branch.
                    stamp_polymer_flux_archetype(forward, reactants, polymer_reactants)
```

Flag-stability note (spec §4.2): the stored flag is written BEFORE the surgery and nothing downstream recomputes it from product stamps — `demote_flipped_polymer_archetype` (model.py:2478) only resets the archetype field on in-place flips (DISCRETE_CHIP included, which is correct: a flip reverses the roles the stamp encodes), and the solver reads the stored fields only.

- [ ] **Step 4: Run to verify pass**

Run: `~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/polymerTest.py -q`
Expected: all passed (+1 new; the Task 3-of-prior-plan handshake/classifier contract tests still pass because `stamp_polymer_flux_archetype` delegates to the same classifier).

- [ ] **Step 5: Run the model suite**

Run: `~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/rmg/modelTest.py test/rmgpy/polymerMultiPoolTest.py -q`
Expected: all passed.

- [ ] **Step 6: Commit**

```bash
git add rmgpy/rmg/model.py test/rmgpy/polymerTest.py
git commit -m "model: run chip surgery in make_new_reaction stamping"
```

---

### Task 6: Solver constant, `reaction_chip_units` array, demotion (spec §5 demotion; tests 9, 15)

**Files:**
- Modify: `rmgpy/solver/polymer.pyx` — constants (~line 77), `initialize_model` array fill (~lines 436-443) and demotion (~lines 555-573)
- Test: `test/rmgpy/solver/solverPolymerTest.py` — `TestHybridPolymerReactor` (extend `test_flux_archetype_constants_match_enum` ~line 572; new test after `test_stamped_scission_without_daughter_pool_demotes_to_legacy` ~line 669)

- [ ] **Step 1: Write the failing tests**

(a) Extend `test_flux_archetype_constants_match_enum` (~line 572) by appending one line inside it:

```python
        assert sp.FLUX_DISCRETE_CHIP == int(PolymerFluxArchetype.DISCRETE_CHIP) == 5
```

(b) Add a new test after `test_stamped_scission_without_daughter_pool_demotes_to_legacy` (~line 669):

```python
    def test_stamped_chip_without_src_pool_demotes_to_unresolved(self):
        """
        Spec test 15: a reaction stamped DISCRETE_CHIP whose reactant does not
        resolve to a solver pool (src == -1) demotes to UNRESOLVED at
        initialize_model -- the same aggregate unresolvable-pool demotion path
        as MIGRATION/SCISSION_FRAGMENT (chip needs only src; there is no dst:
        the complement folds back and the chip is a gas species).
        """
        Proxy = _spc("CCCC", "poly")
        Mu0 = _spc("CO", "poly_mu0")
        Mu1 = _spc("C=O", "poly_mu1")
        Mu2 = _spc("C#N", "poly_mu2")
        A = _spc("C", "A_gas")
        B = _spc("[CH3]", "B_gas")
        core = [Proxy, Mu0, Mu1, Mu2, A, B]
        mask = np.array([False] * 4 + [True, True], dtype=bool)

        rxn = Reaction(reactants=[A], products=[B], **_KIN)   # no pool reactant
        rxn.polymer_flux_archetype = 5
        rxn.polymer_chip_units = 2

        pool = PolymerPoolConfig(
            label="poly", xs=2, explicit_dp_to_species_index={},
            mu_indices=(1, 2, 3), monomer_poly_index=None,
            k_scission=0.0, k_unzip=0.0, tail_kinetics=None,
        )
        rs = HybridPolymerSystem(
            T=800.0, P=1.0e5, initial_mole_fractions={A: 1.0}, V_poly=1.0,
            polymer_pools=[pool], mass_transfer=[],
            gas_species_mask=mask.copy(), constant_gas_volume=False,
            initial_polymer_moments={"poly": (1.0, 5.0, 30.0)}, termination=[],
        )
        rs.initialize_model(core, [rxn], [], [])

        import rmgpy.solver.polymer as sp
        assert rs.reaction_src_pool[0] == -1
        assert rs.reaction_flux_archetype[0] == sp.FLUX_UNRESOLVED   # demoted
        assert rs.reaction_chip_units[0] == 2    # array filled regardless
```

- [ ] **Step 2: Run to verify failure**

Run: `~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/solver/solverPolymerTest.py -q -k "constants_match_enum or chip_without_src_pool"`
Expected: 2 FAIL — `AttributeError: ... 'FLUX_DISCRETE_CHIP'` and `'reaction_chip_units'`.

- [ ] **Step 3: Implement**

(a) In `rmgpy/solver/polymer.pyx`, the constants block (lines 74-81) gains one line:

```python
# Pool moment-flux archetypes. Mirror of rmgpy.polymer.PolymerFluxArchetype
# (not imported to avoid a solver->polymer module cycle); equality is pinned
# by test_flux_archetype_constants_match_enum.
FLUX_NONE = 0
FLUX_SAME_POOL = 1
FLUX_MIGRATION = 2
FLUX_SCISSION_FRAGMENT = 3
FLUX_UNRESOLVED = 4
FLUX_DISCRETE_CHIP = 5
```

(b) In `initialize_model`, the archetype fill loop (lines ~439-443) currently:

```python
        self.reaction_flux_archetype = np.zeros(n_rxn, dtype=np.int8)
        self.reaction_src_pool = np.full(n_rxn, -1, dtype=np.int32)
        self.reaction_dst_pool = np.full(n_rxn, -1, dtype=np.int32)
        for i, rxn in enumerate(itertools.chain(core_reactions, edge_reactions)):
            self.reaction_flux_archetype[i] = int(getattr(rxn, "polymer_flux_archetype", 0))
```

becomes:

```python
        self.reaction_flux_archetype = np.zeros(n_rxn, dtype=np.int8)
        self.reaction_src_pool = np.full(n_rxn, -1, dtype=np.int32)
        self.reaction_dst_pool = np.full(n_rxn, -1, dtype=np.int32)
        # Stamped chip repeat-unit counts (spec 2026-06-10 §4.3); same
        # chain(core, edge) order so the index matches r_idx in the residual.
        self.reaction_chip_units = np.zeros(n_rxn, dtype=np.int32)
        for i, rxn in enumerate(itertools.chain(core_reactions, edge_reactions)):
            self.reaction_flux_archetype[i] = int(getattr(rxn, "polymer_flux_archetype", 0))
            self.reaction_chip_units[i] = int(getattr(rxn, "polymer_chip_units", 0))
```

(c) Extend the unresolvable-pool demotion (lines ~555-568). Currently:

```python
            if (self.reaction_flux_archetype[i] in (FLUX_MIGRATION, FLUX_SCISSION_FRAGMENT)
                    and (src == -1 or dst == -1)):
```

becomes (keep the entire existing comment block inside the `if` body unchanged, and append the chip clause to the condition + one comment line):

```python
            if ((self.reaction_flux_archetype[i] in (FLUX_MIGRATION, FLUX_SCISSION_FRAGMENT)
                    and (src == -1 or dst == -1))
                    or (self.reaction_flux_archetype[i] == FLUX_DISCRETE_CHIP
                        and src == -1)):
```

and inside the body, after the existing comment, add:

```python
                # DISCRETE_CHIP needs only src (no dst: complement folds back
                # to the same pool and the chip is a plain gas species).
```

(d) Update the aggregate warning text (lines ~569-573) from `"stamped MIGRATION/SCISSION_FRAGMENT could not resolve both solver pools (src/dst)"` to:

```python
        if n_unresolvable:
            logging.warning(
                "%d reactions stamped MIGRATION/SCISSION_FRAGMENT/DISCRETE_CHIP "
                "could not resolve their solver pool(s); demoted to legacy "
                "mu1-only moment flux (UNRESOLVED).", n_unresolvable)
```

- [ ] **Step 4: Rebuild and run**

Run: `make`
Then: `~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/solver/solverPolymerTest.py -q`
Expected: all passed (existing suite green, +1 new test, 1 extended).

- [ ] **Step 5: Commit**

```bash
git add rmgpy/solver/polymer.pyx test/rmgpy/solver/solverPolymerTest.py
git commit -m "solver: resolve chip units array and demotion at init"
```

---

### Task 7: DISCRETE_CHIP dispatch branch in the residual (spec §5; tests 10, 11, 12, 14)

**Files:**
- Modify: `rmgpy/solver/polymer.pyx` — section-5 archetype dispatch (~lines 1090-1178), new branch between `FLUX_SCISSION_FRAGMENT` and `FLUX_UNRESOLVED`
- Test: `test/rmgpy/solver/solverPolymerTest.py` — `TestHybridPolymerReactor`, after `test_bundle_guard_mu3_overflow_skips_mu2_component_only` (~line 1115); reuses the module fixtures `_two_pool_species` / `_two_pool_rs` / `_KIN` (lines 48-80) and `constants` (imported at line 33)

- [ ] **Step 1: Write the failing tests**

```python
    def test_discrete_chip_monodisperse_closed_form_both_picks(self):
        """
        Spec test 10: monodisperse pool (mu_j = N*L^j) -> E[n] = L under BOTH
        picks (uniform mu1/mu0 == length-biased mu2/mu1 == L on monodisperse),
        so the chip drain has the same closed form for either flag value:
        dmu0 = 0, dmu1 = -r*a, dmu2 = -r*(2aL - a^2). Only the rate's site
        scaling differs (mu0 vs mu1). The chip species flows through the
        standard gas path (+r).
        """
        N, L, a = 2.0, 10.0, 3
        for end_group in (False, True):
            sp, core, mask = _two_pool_species()
            rxn = Reaction(reactants=[sp["A"]], products=[sp["A"], sp["G"]], **_KIN)
            rxn.polymer_flux_archetype = 5
            rxn.polymer_chip_units = a
            rxn.is_end_group_reaction = end_group
            rs = _two_pool_rs(rxn, core, mask, (N, N * L, N * L * L), (0.1, 0.2, 0.5))

            dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]

            kf = rxn.get_rate_coefficient(800.0, 1.0e5)
            r = kf * (N if end_group else N * L)          # site scaling per flag
            assert np.isclose(dn_dt[1], 0.0, atol=1e-14), end_group   # mu0
            assert np.isclose(dn_dt[2], -r * a), end_group            # mu1
            assert np.isclose(dn_dt[3], -r * (2.0 * a * L - a * a)), end_group
            assert np.allclose(dn_dt[5:8], 0.0, atol=1e-14)           # pool B idle
            assert np.isclose(dn_dt[8], +r)                           # chip gas

    def test_discrete_chip_uses_scaling_consistent_e_n(self):
        """
        Spec test 11 (decision D5): polydisperse pool (1, 5, 30) separates the
        picks -- uniform E[n] = mu1/mu0 = 5, length-biased E[n] = mu2/mu1 = 6.
        The mu2 drain must use the E[n] matching the reaction's rate-scaling
        flag; pairing a mu0 rate with length-biased E[n] (or vice versa)
        fails this test loudly.
        """
        a = 1
        for end_group, e_n, site in ((False, 6.0, 5.0), (True, 5.0, 1.0)):
            sp, core, mask = _two_pool_species()
            rxn = Reaction(reactants=[sp["A"]], products=[sp["A"], sp["G"]], **_KIN)
            rxn.polymer_flux_archetype = 5
            rxn.polymer_chip_units = a
            rxn.is_end_group_reaction = end_group
            rs = _two_pool_rs(rxn, core, mask, (1.0, 5.0, 30.0), (0.1, 0.2, 0.5))

            dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]

            kf = rxn.get_rate_coefficient(800.0, 1.0e5)
            r = kf * site
            assert np.isclose(dn_dt[2], -r * a), end_group
            assert np.isclose(dn_dt[3], -r * (2.0 * a * e_n - a * a)), end_group

    def test_discrete_chip_clamp_regime_skips_mu2_write(self):
        """
        Spec test 12: a >= 2*E[n] makes 2aE[n] - a^2 <= 0 -- impossible
        per-chain (n >= a) but reachable in expectation for a pool whose mean
        length decayed toward chip size. The forward mu2 decrement is clamped
        (write skipped), mu1 still drains, RHS finite.
        """
        a = 13   # length-biased E[n] = 6 -> 2*13*6 - 169 = -13 < 0
        sp, core, mask = _two_pool_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["A"], sp["G"]], **_KIN)
        rxn.polymer_flux_archetype = 5
        rxn.polymer_chip_units = a
        rs = _two_pool_rs(rxn, core, mask, (1.0, 5.0, 30.0), (0.1, 0.2, 0.5))

        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]

        kf = rxn.get_rate_coefficient(800.0, 1.0e5)
        r = kf * 5.0                                  # mu1-scaled
        assert np.all(np.isfinite(dn_dt))
        assert np.isclose(dn_dt[2], -r * a)           # mu1 still drains
        assert dn_dt[3] == 0.0                        # mu2 write skipped
        assert np.isclose(dn_dt[1], 0.0, atol=1e-14)  # mu0 untouched

    def test_discrete_chip_reverse_leg_exact_extension_form(self):
        """
        Spec test 14: the reverse leg is the EXACT extension form
        Delta(n^2) = (n+a)^2 - n^2 = +(2aE[n] + a^2): dmu1 = +rr*a,
        dmu2 = +rr*(2aE[n] + a^2) -- PLUS a^2, not the forward sign-flip
        (which would subtract it) -- and never clamps, even at a >= 2*E[n]
        where the forward leg would. Driven via the rs.kb override with
        kf = 0 and injected chip moles.
        """
        a = 13                                # forward would clamp here
        sp, core, mask = _two_pool_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["A"], sp["G"]], **_KIN)
        rxn.polymer_flux_archetype = 5
        rxn.polymer_chip_units = a
        rxn.is_end_group_reaction = True      # uniform pick: E[n] = mu1/mu0 = 5
        rs = _two_pool_rs(rxn, core, mask, (1.0, 5.0, 30.0), (0.1, 0.2, 0.5))
        rs.kf[0] = 0.0                        # silence the forward leg
        rs.kb[0] = 0.4                        # drive the reverse leg directly

        y = rs.y.copy()
        y[8] = 1.0                            # inject chip (G) moles

        dn_dt = rs.residual(0.0, y, np.zeros_like(y))[0]

        # G is the only gas species with moles -> C_G = P/(R*T) exactly.
        c_g = 1.0e5 / (constants.R * 800.0)
        site = 1.0                            # mu0 of pool A (flag True)
        rr = 0.4 * 1.0 * c_g * site           # kb * C(fold-back proxy)=1 * C(G) * site
        e_n = 5.0
        assert np.isclose(dn_dt[1], 0.0, atol=1e-14)              # mu0 untouched
        assert np.isclose(dn_dt[2], +rr * a)                      # exact form
        assert np.isclose(dn_dt[3], +rr * (2.0 * a * e_n + a * a))  # +a^2, no clamp
        assert np.isclose(dn_dt[8], -rr)      # chip consumed via standard path
```

- [ ] **Step 2: Run to verify failure**

Run: `~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/solver/solverPolymerTest.py -q -k "discrete_chip"`
Expected: 4 FAIL. Today a stamped archetype 5 falls through the dispatch chain without matching any branch (no chip branch, and 5 != FLUX_UNRESOLVED), so all pool-moment derivatives are zero — e.g. test 10 sees `dn_dt[2] == 0` instead of `-r*a`. If instead the residual CRASHES, STOP and check Task 6's array wiring.

- [ ] **Step 3: Implement the dispatch branch**

In `rmgpy/solver/polymer.pyx`, in the section-5 dispatch inside `residual`, insert between the end of the `elif arch == FLUX_SCISSION_FRAGMENT:` block (ends ~line 1164) and `elif arch == FLUX_UNRESOLVED:` (~line 1165):

```python
                elif arch == FLUX_DISCRETE_CHIP:
                    # Discrete chip (spec 2026-06-10 §5): an end-anchored cut
                    # ejects a stamped a-unit chip to the gas phase; the
                    # complement folds back into the SAME pool, so mu0 is
                    # unchanged and there is no dst pool. Closure-free: only
                    # E[n] (an exact ratio of tracked moments) is needed --
                    # no _safe_mu3_from_mu012 call, no mu2_ok branch. The
                    # bundle and the rate scaling key on the same stored flag,
                    # so they cannot disagree by construction. The chip
                    # species itself gains/loses moles through the standard
                    # gas dn_dt path (section 4) -- no special handling here.
                    src = self.reaction_src_pool[r_idx]
                    # src == -1 cannot reach here (init demotes); defensive.
                    if src != -1:
                        chip_a = float(self.reaction_chip_units[r_idx])
                        b0, b1, _b2, _mu2_ok = self._chain_bundle(
                            src, y, V_poly, self.is_end_group_reaction[r_idx])
                        if b0 != 0.0:
                            s_idx = self.polymer_pools[src].mu_indices
                            chip_e_n = b1
                            if rf > 0.0:
                                # Forward (ejection): Delta n = -a,
                                # Delta n^2 = -(2aE[n] - a^2). Clamp the mu2
                                # decrement at >= 0: 2aE[n] - a^2 < 0 is
                                # impossible per-chain (n >= a always) but
                                # reachable in expectation when the pool mean
                                # length has decayed toward chip size; by
                                # then the moment description is marginal.
                                rf_mol = rf * V_rxn
                                dn_dt[s_idx[1]] -= rf_mol * chip_a
                                chip_dmu2 = 2.0 * chip_a * chip_e_n - chip_a * chip_a
                                if chip_dmu2 > 0.0:
                                    dn_dt[s_idx[2]] -= rf_mol * chip_dmu2
                            if rr > 0.0:
                                # Reverse (re-attachment) -- EXACT extension
                                # form, NOT the forward sign-flip:
                                # (n+a)^2 - n^2 = +(2aE[n] + a^2). E[n] on the
                                # same single pool (the re-formed chain
                                # extends a parent-statistics chain by a).
                                # Unconditionally positive: no clamp.
                                rr_mol = rr * V_rxn
                                dn_dt[s_idx[1]] += rr_mol * chip_a
                                dn_dt[s_idx[2]] += rr_mol * (
                                    2.0 * chip_a * chip_e_n + chip_a * chip_a)
```

Notes for the implementer:
- `a = 0` chips degrade to zero flux naturally — no special case (spec §3.4).
- `b0 == 0.0` is the inherited empty-pool sentinel from `_chain_bundle` — skip.
- The `_chain_bundle` call computes the µ3 closure internally for the length-biased pick; its `b2`/`mu2_ok` outputs are deliberately ignored (`_b2`, `_mu2_ok`) — the chip formulas never use them, which is the closure-free property the spec claims. Do NOT add a `mu2_ok` gate here.
- Do not reuse the names `e_n`/`mu3_p` from the SCISSION branch — the `chip_*` names avoid any cross-branch aliasing confusion in this long function.

- [ ] **Step 4: Rebuild and run the new tests**

Run: `make && ~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/solver/solverPolymerTest.py -q -k "discrete_chip"`
Expected: 4 passed.

- [ ] **Step 5: Run the full solver suite (no regressions)**

Run: `~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/solver/solverPolymerTest.py -q`
Expected: all passed.

- [ ] **Step 6: Commit**

```bash
git add rmgpy/solver/polymer.pyx test/rmgpy/solver/solverPolymerTest.py
git commit -m "solver: dispatch DISCRETE_CHIP moment flux per direction"
```

---

### Task 8: Conservation trajectory, dormant knob, Cantera audit, docs (spec §6, §7, §8 tests 13/16/17/17b, §10 A1–A3)

**Files:**
- Test: `test/rmgpy/solver/solverPolymerTest.py` (test 13), `test/rmgpy/polymerTest.py` (test 16 + knob), `test/rmgpy/canteraTest.py` (tests 17/17b, class `TestCanteraWriter`)
- Modify: `rmgpy/polymer.py` (`__init__` + `copy()` for the knob), `rmgpy/rmg/input.py` (`polymer()` kwarg), `docs/multi_pool_design.md`

- [ ] **Step 1: Write the conservation trajectory test (spec test 13)**

Add to `TestHybridPolymerReactor` in `test/rmgpy/solver/solverPolymerTest.py` after the Task 7 tests:

```python
    def test_discrete_chip_trajectory_conserves_units_and_cone(self):
        """
        Spec test 13 -- the exact conservation invariant (closed system,
        chip reactions only): d/dt [ Sigma_pools mu1 + Sigma_chips a_i*n_i ]
        = 0, chip moles weighted by the STAMPED unit count (not raw moles).
        Over a chip-only forward-Euler trajectory: the weighted invariant is
        constant to roundoff, mu0 never changes (chain count preserved), and
        the pool stays in the realizable cone (mu1 >= mu0 >= 0).
        """
        sp, core, mask = _two_pool_species()
        a = 3
        rxn = Reaction(reactants=[sp["A"]], products=[sp["A"], sp["G"]], **_KIN)
        rxn.polymer_flux_archetype = 5
        rxn.polymer_chip_units = a
        rxn.is_end_group_reaction = True
        rs = _two_pool_rs(rxn, core, mask, (1.0, 50.0, 3000.0), (0.0, 0.0, 0.0))

        y = rs.y.copy()
        invariant0 = y[2] + a * y[8]          # mu1_A + a * n_chip
        # r = kf*mu0 = 2/s -> dmu1 = -6/s; 200 * 1e-4 s drains ~0.12 of 50:
        # far from depletion, forward Euler is safe.
        dt = 1e-4
        for _ in range(200):
            dn_dt = rs.residual(0.0, y, np.zeros_like(y))[0]
            assert np.all(np.isfinite(dn_dt))
            y = y + dt * dn_dt
            assert y[1] >= -1e-12             # mu0 >= 0
            assert y[2] - y[1] >= -1e-9       # cone: mu1 >= mu0
        assert np.isclose(y[2] + a * y[8], invariant0, rtol=1e-12)
        assert np.isclose(y[1], 1.0)          # chain count unchanged
        assert y[8] > 0.0                     # chips accumulated
```

Run: `~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/solver/solverPolymerTest.py -q -k trajectory_conserves_units`
Expected: PASS (Task 7's branch makes the invariant exact by construction; if it fails, that is a real Task 7 bug — STOP and fix there, do not tune this test).

- [ ] **Step 2: Write the failing knob test, then implement `discrete_dp_threshold`**

Add to `TestPolymer` in `test/rmgpy/polymerTest.py` (after the Task 4 tripwire test):

```python
    def test_discrete_dp_threshold_config_field(self):
        """
        discrete_dp_threshold (spec 2026-06-10 §6, decisions D7/D8): per-pool
        config knob, default 4 (monomer through trimer explicit), DORMANT
        under the fixed trimer proxy -- no behavioral use yet. Stored on the
        Polymer, survives copy(), overridable per pool.
        """
        assert self.polymer_1.discrete_dp_threshold == 4
        assert self.polymer_1.copy().discrete_dp_threshold == 4
        p = Polymer(label='PE_thresh', monomer='[CH2][CH2]',
                    end_groups=['[H]', '[H]'], cutoff=3,
                    Mn=1000.0, Mw=2500.0, initial_mass=1.0,
                    discrete_dp_threshold=6)
        assert p.discrete_dp_threshold == 6
```

Run to verify FAIL (`AttributeError: ... 'discrete_dp_threshold'`), then implement:

(a) `rmgpy/polymer.py` `__init__` (~lines 235-238) — extend the kwargs.pop block:

```python
        # Pop polymer-specific kwargs before passing the rest to Species.__init__
        # — Species does not accept k_unzip/k_scission and would raise TypeError.
        self.k_unzip = kwargs.pop('k_unzip', 0.0)
        self.k_scission = kwargs.pop('k_scission', 0.0)
        # Discreteness threshold (spec 2026-06-10 §6, D7/D8): chains with
        # literal DP < threshold are candidates for discrete tracking. Default
        # 4 = monomer..trimer explicit. DORMANT under the fixed trimer proxy:
        # the backbone gate is proxy-relative (D2), scission routing keys on
        # is_end_group_reaction (D3), and the conditional DP backstop (D8)
        # activates only when the proxy repeat-count exceeds this threshold
        # (a 3-unit proxy never does). Defined-but-documented beats undefined
        # intent; no behavioral use yet.
        self.discrete_dp_threshold = kwargs.pop('discrete_dp_threshold', 4)
```

(b) `rmgpy/polymer.py` `copy()` (~line 491, beside the k_scission/k_unzip carry-over):

```python
        other.k_scission = self.k_scission
        other.k_unzip = self.k_unzip
        other.discrete_dp_threshold = getattr(self, 'discrete_dp_threshold', 4)
```

(c) `rmgpy/rmg/input.py` `polymer()` (~lines 248-292) — add the kwarg and pass it through:

```python
def polymer(label: str,
            monomer: Union[Molecule, str],
            end_groups: Optional[List[Union[str, Molecule]]],
            cutoff: int = 4,
            Mn: Optional[float] = None,
            Mw: Optional[float] = None,
            moments: Optional[List[float]] = None,
            initial_mass: float = 1.0,
            k_scission: float = 0.0,
            k_unzip: float = 0.0,
            discrete_dp_threshold: int = 4,
            ):
```

docstring addition (after the `k_unzip` entry):

```python
        discrete_dp_threshold (int): Chains with literal DP below this are
            candidates for discrete (explicit-species) tracking instead of
            pool moments. Default 4. DORMANT under the fixed trimer proxy
            (see docs/multi_pool_design.md §5.1) -- reserved for the
            conditional DP backstop with longer proxies.
```

and in the `Polymer(...)` construction:

```python
        k_scission=k_scission,
        k_unzip=k_unzip,
        discrete_dp_threshold=discrete_dp_threshold,
    )
```

Run: `~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/polymerTest.py -q -k discrete_dp_threshold`
Expected: 1 passed.

- [ ] **Step 3: Write the backstop-dormancy test (spec test 16 — no production code)**

Add to `TestPolymer` after the knob test:

```python
    def test_backstop_dormant_under_trimer_proxy(self):
        """
        Spec test 16 / decision D8: the conditional DP backstop ("mu1-scaled
        cut producing literal DP < threshold -> exact-a accounting") applies
        ONLY when the proxy repeat-count exceeds discrete_dp_threshold. The
        fixed trimer proxy (3 units) never exceeds the default threshold (4),
        so a mu1-scaled mid-cut keeps routing SCISSION_FRAGMENT regardless of
        literal piece DP. No backstop code exists yet -- this pins the
        routing so adding it later cannot silently activate on trimer decks
        (unconditional, it would route ALL mid-chain scission to chip and
        kill mu0 growth/Mn halving).
        """
        from rmgpy.polymer import PolymerFluxArchetype, classify_reaction_flux_archetype

        p = self.polymer_1
        assert p.discrete_dp_threshold == 4          # >= 3 proxy repeat units
        gas = Molecule(smiles="CC")
        piece = p.copy()
        piece.label = f"{p.label}_scission_tail"
        piece._reacted_class = PolymerClass.SCISSION
        piece._source_molecule = Molecule(smiles="CCC(C)c1ccccc1")  # DP ~1 piece
        assert (classify_reaction_flux_archetype([p], [piece, gas])
                == PolymerFluxArchetype.SCISSION_FRAGMENT)
```

Run: `~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/polymerTest.py -q -k backstop_dormant`
Expected: PASS immediately (routing pin; the test exists so a future backstop implementation must consciously satisfy the proxy-repeat-count > threshold condition).

- [ ] **Step 4: Write the Cantera audit tests (spec tests 17/17b — no production code expected)**

Add to `TestCanteraWriter` in `test/rmgpy/canteraTest.py`, after `test_generate_cantera_data_drops_unbalanced_polymer_proxy_reaction` (~line 306):

```python
    def test_generate_cantera_data_keeps_balanced_chip_reaction(self):
        """
        Spec test 17: a sub-shape-(b) chip reaction's fold-back is a
        structurally MODIFIED END_MOD image, so proxy -> modified-proxy + chip
        can be atom-exact. Balanced (b)-shapes must SURVIVE the
        unbalanced-proxy filter (it keys on is_balanced(), not on "is a proxy
        scission").
        """
        parent = self._create_dummy_species("parent", "CCC", index=1)
        parent.is_polymer_proxy = True
        fold_b = self._create_dummy_species("parent_mod", "[CH2]CC", index=2)
        fold_b.is_polymer_proxy = True
        h = self._create_dummy_species("H", "[H]", index=3)

        chip_rxn_b = Reaction(
            reactants=[parent], products=[fold_b, h],
            kinetics=Arrhenius(A=(1e13, "s^-1"), n=0, Ea=(300, "kJ/mol"), T0=(1, "K")),
        )
        chip_rxn_b.polymer_flux_archetype = 5
        chip_rxn_b.polymer_chip_units = 0
        assert chip_rxn_b.is_balanced()   # sanity: (b)-shapes can balance

        species = [parent, fold_b, h]
        data = generate_cantera_data(species, [chip_rxn_b], is_plasma=False)
        equations = [d['equation'] for d in data.get('reactions', [])]
        assert len(equations) == 1
        assert any('parent_mod' in eq for eq in equations)

    def test_generate_cantera_data_drops_overbalanced_chip_reaction(self, caplog):
        """
        Spec test 17b: a sub-shape-(a) chip reaction's fold-back is an
        UNMODIFIED parent copy, so it reads proxy -> proxy + chip:
        over-balanced by the chip mass (the spec-A2 cap-mass drift made
        visible at the species-balance level). It must be dropped and counted
        exactly like a registered unbalanced scission.
        """
        import logging as _logging

        parent = self._create_dummy_species("parent", "CCC", index=1)
        parent.is_polymer_proxy = True
        chip = self._create_dummy_species("chip", "C", index=2)

        chip_rxn_a = Reaction(
            reactants=[parent], products=[parent, chip],
            kinetics=Arrhenius(A=(1e13, "s^-1"), n=0, Ea=(300, "kJ/mol"), T0=(1, "K")),
        )
        chip_rxn_a.polymer_flux_archetype = 5
        chip_rxn_a.polymer_chip_units = 1
        assert not chip_rxn_a.is_balanced()   # over-balanced by the chip

        r = self._create_dummy_species("R", "[CH2]O", index=3)
        p = self._create_dummy_species("P", "C[O]", index=4)
        normal = Reaction(
            reactants=[r], products=[p],
            kinetics=Arrhenius(A=(1e10, "s^-1"), n=0, Ea=(0, "J/mol"), T0=(1, "K")),
        )

        species = [parent, chip, r, p]
        with caplog.at_level(_logging.INFO):
            data = generate_cantera_data(species, [chip_rxn_a, normal], is_plasma=False)
        equations = [d['equation'] for d in data.get('reactions', [])]
        assert "R(3) <=> P(4)" in equations
        assert len(equations) == 1            # chip (a)-shape dropped
        assert any("dropped 1 element-unbalanced" in rec.getMessage()
                   for rec in caplog.records) # counted like a scission drop
```

Run: `~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/canteraTest.py -q -k "chip_reaction"`
Expected: 2 passed WITHOUT touching `rmgpy/cantera.py` — the existing filter (`generate_cantera_data`, ~lines 216-245: drop iff `_involves_polymer_proxy(rxn) and not rxn.is_balanced()`) already implements exactly the spec-§7 behavior. If either fails, the filter keys on something other than balance — STOP, read `rmgpy/cantera.py:200-245`, and fix the FILTER (not the test) to be balance-keyed. Also audit (read-only) that no other balance validator special-cases proxy reactions: `grep -n "is_balanced" rmgpy/cantera.py rmgpy/rmg/*.py` — anything keyed on "proxy scission ⇒ unbalanced" must be reported in the final summary.

- [ ] **Step 5: Documentation — `docs/multi_pool_design.md`**

Make these edits:

(a) In the §5.1 table, append a row after SCISSION_FRAGMENT and extend the UNRESOLVED row:

```markdown
| DISCRETE_CHIP | end-anchored cut ejects a stamped `a`-unit discrete chip (`Reaction.polymer_chip_units`, `a = round(chip_MW/monomer_MW)`, a = 0 legal); the complement folds back into the SAME pool. Per direction: forward (0, −rf·a, −rf·max(0, 2a·E[n]−a²)) with the µ2 decrement clamped at ≥ 0; reverse uses the EXACT extension form (0, +rr·a, +rr·(2a·E[n]+a²)) — plus a², never clamped. E[n] is scaling-aware (µ1/µ0 when µ0-scaled, µ2/µ1 when µ1-scaled — the rate scaling is the sampling measure). Closure-free (no µ3). The chip species flows through the standard gas path; exact invariant: d/dt[Σ_pools µ1 + Σ_chips a_i·n_i] = 0. |
```

In the UNRESOLVED row, change "stamped MIGRATION/SCISSION_FRAGMENT whose src/dst pool cannot be resolved" to "stamped MIGRATION/SCISSION_FRAGMENT whose src/dst pool (or DISCRETE_CHIP whose src pool) cannot be resolved".

(b) After the §5.1 table paragraphs, add:

```markdown
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
threshold AND piece DP < threshold) only activates for longer proxies.
```

(c) In §4.5 (Spawn drain), append one line:

```markdown
Chip events never queue a spawn intent: `surge_chip_products` replaces the
SCISSION daughter at stamping time, upstream of this iteration-boundary hook.
```

(d) In §10 (Known limitations), amend item 7 by appending: "*Partially addressed 2026-06-10:* the discreteness gate rejects backbone impostors (heavy-atom bound) and DISCRETE_CHIP tracks end-cut fragments discretely; true volatility (Tb vs reactor T) remains a mass-transfer-layer item." Then append three new numbered items:

```markdown
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
```

(e) In §11 (Open items), add:

```markdown
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
```

- [ ] **Step 6: Commit Task 8 code + tests, then docs**

```bash
git add rmgpy/polymer.py rmgpy/rmg/input.py test/rmgpy/polymerTest.py test/rmgpy/solver/solverPolymerTest.py test/rmgpy/canteraTest.py
git commit -m "polymer: chip conservation tests and dormant dp knob"
git add docs/multi_pool_design.md
git commit -m "docs: discreteness gate and chip archetype notes"
```

- [ ] **Step 7: Full verification — all six affected suites**

Run:
```bash
~/anaconda3/envs/rmg_env/bin/python -m pytest \
  test/rmgpy/polymerTest.py \
  test/rmgpy/solver/solverPolymerTest.py \
  test/rmgpy/canteraTest.py \
  test/rmgpy/polymerMultiPoolTest.py \
  test/rmgpy/reactionTest.py \
  test/rmgpy/rmg/modelTest.py -q
```
Expected: 0 failures (this plan adds ~17 tests across the suites). Investigate ANY failure before proceeding — only a test that directly asserted the OLD end-initiated-scission UNRESOLVED warn reason string ("end-initiated scission") may be updated to the new reason ("unsurged end-initiated scission"), and Task 2 already checked there is none.

- [ ] **Step 8: EPDM end-to-end no-op gate**

Run:
```bash
cd ~/runs/RMG/epdm_v0_2026-06-06b && ~/anaconda3/envs/rmg_env/bin/python ~/Code/RMG-Py/rmg.py input.py
```
Expected: `MODEL GENERATION COMPLETED`, **26 core species / 28 reactions** — EPDM has no flag-true scission shapes, so this work is a true no-op for it (spec §8). Also scan the log: no `backbone_impostor` reclassification of previously-polymer species (would change the species count), no flood of the new warn-once messages. Then `cd ~/Code/RMG-Py`.

- [ ] **Step 9: Update the auto-memory running log**

Edit `~/.claude/projects/-home-alon-Code-RMG-Py/memory/project_polymer_branch_review.md`: mark open item #5 (gas/polymer split is topological only) as IMPLEMENTED with the commit range and the headline test names, mirroring the style of the existing FIXED entries.

---

## Verification summary (spec → task map)

| Spec requirement | Task |
|---|---|
| §3.1 backbone impostor gate (one-sided heavy-atom bound, `backbone_impostor` reason, no ceiling + V3 verification) — test 1 | 1 |
| §3.3 `PolymerClass.CHIP`; §4.1 `DISCRETE_CHIP = 5`, CHIP-before-SCISSION classifier branch, reworded unsurged fallback — test 2 | 2 |
| §4.3 `Reaction.polymer_chip_units` (pxd + kwarg, NOT in `__reduce__`) — test 8 | 3 |
| §4.2 surgery both sub-shapes, smaller-piece identification, MW-ratio `a`, a=0 legal, infeasible→UNRESOLVED never SCISSION_FRAGMENT; §4.4 tripwire — tests 3, 3b, 4, 5, 6 | 4 |
| §4.2 never-queue + §4.5 flag stability at the model level — test 7 (+3b rider) | 5 |
| §4.1 mirror constant; chip-units array in chain(core,edge) order; §5 demotion (DISCRETE_CHIP src==−1) — tests 9, 15 | 6 |
| §5 dispatch: forward clamp, reverse exact extension form, scaling-aware E[n], closure-free, b0 sentinel — tests 10, 11, 12, 14 | 7 |
| §5 conservation invariant trajectory — test 13; §6 `discrete_dp_threshold` knob (dormant) + D8 backstop-dormancy — test 16; §7 Cantera audit — tests 17, 17b; §10 A1–A3 + §9 detector item → docs; EPDM no-op gate | 8 |
