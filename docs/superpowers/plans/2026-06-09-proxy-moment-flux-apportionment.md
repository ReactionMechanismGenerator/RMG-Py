# Proxy Moment-Flux Apportionment Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Discrete proxy reactions apportion their molar event flux across all three pool moments (µ0/µ1/µ2) by per-reaction archetype, instead of pushing everything into µ1; edge reactions stop perturbing the integrated core state.

**Architecture:** Generation-time archetype stamping (`Reaction.polymer_flux_archetype`, mirroring the `is_end_group_reaction` 4-layer wiring from the µ0-scaling work) → solver int8/int32 arrays resolved at `initialize_model` → archetype dispatch in the RHS residual applying chain-statistics "bundles". Spec: `docs/superpowers/specs/2026-06-09-proxy-moment-flux-apportionment-design.md` (READ IT FIRST — it has the math derivations and decision rationale).

**Tech Stack:** RMG-Py `polymer` branch, Python 3.9 (`~/anaconda3/envs/rmg_env/bin/python`), Cython (`.pyx`/`.pxd` — **run `make` from the repo root after ANY edit to a `.pyx` or `.pxd` file, or tests will silently run stale compiled code**), pytest.

**Working directory:** `/home/alon/Code/RMG-Py` (all paths below relative to it).

**Conventions:**
- Run tests with `python -m pytest <file>::<class>::<test> -q` using the rmg_env python.
- Commit style: lowercase `topic: imperative subject` ≤60 chars, body only when the why is non-obvious. No `git add -A`; stage explicit paths. Never push.
- The four polymer suites: `test/rmgpy/polymerTest.py test/rmgpy/solver/solverPolymerTest.py test/rmgpy/canteraTest.py test/rmgpy/polymerMultiPoolTest.py` — currently **191 passed**; this plan adds ~17 tests, every run must end with 0 failures.

**Key chemistry (from the spec, restated for the implementer):**

A "bundle" b = (b0, b1, b2) is the expected (chains, monomer units, units²) carried by ONE picked chain:
- µ1-scaled reaction (rate ∝ site density → the picked chain is length-biased): b = (1, µ2/µ1, µ3/µ1), µ3 from the existing `_safe_mu3_from_mu012` closure.
- µ0-scaled reaction (`is_end_group_reaction` → uniform chain pick): b = (1, µ1/µ0, µ2/µ0).

| Archetype | int | Pool moment flux |
|---|---|---|
| NONE | 0 | none (no proxy in reaction) |
| SAME_POOL | 1 | none (product folds back; −r/+r on the same pool cancels exactly, so skip) |
| MIGRATION | 2 | per direction: A −rf·b_A, B +rf·b_A (forward); A +rr·b_B, B −rr·b_B (reverse) |
| SCISSION_FRAGMENT | 3 | parent: (0, −r·µ2/(2µ1), −r·(2/3)·µ3/µ1); daughter: (+r, +r·µ2/(2µ1), +r·µ3/(3µ1)) — "complement stays in parent" |
| UNRESOLVED | 4 | legacy µ1-only transfer (today's behavior) + one-time warning |

---

## File Structure

| File | Change |
|---|---|
| `rmgpy/polymer.py` | Add `PolymerFluxArchetype` IntEnum + `classify_reaction_flux_archetype()` helper + warn-once registry |
| `rmgpy/reaction.pxd` | Declare `cdef public int polymer_flux_archetype` |
| `rmgpy/reaction.py` | `polymer_flux_archetype=0` kwarg + assignment (NOT added to `__reduce__` — mirrors `is_end_group_reaction`; the init-time remap covers unpickled reactions) |
| `rmgpy/rmg/model.py` | Stamp the archetype in `make_new_reaction` (both branches, beside the `is_end_group_reaction` stamps at lines 644 and 683) |
| `rmgpy/solver/polymer.pyx` | Module FLUX_* constants; `initialize_model`: archetype/src/dst arrays + NONE→UNRESOLVED remap; `_chain_bundle()` method; residual: edge gating + archetype dispatch |
| `docs/multi_pool_design.md` | §5 subsection documenting the apportionment + k_scission input-hygiene caveat |
| Tests | `test/rmgpy/polymerTest.py`, `test/rmgpy/reactionTest.py`, `test/rmgpy/solver/solverPolymerTest.py` |

---

### Task 1: Archetype enum + classifier helper (`rmgpy/polymer.py`)

**Files:**
- Modify: `rmgpy/polymer.py` (add after `is_end_group_reaction`, ~line 1395)
- Test: `test/rmgpy/polymerTest.py` (add after `test_is_end_group_reaction_helper`, ~line 760, same class — it provides `self.polymer_1`)

- [ ] **Step 1: Write the failing tests**

Add to `test/rmgpy/polymerTest.py` right after `test_is_end_group_reaction_helper` (the class whose setup provides `self.polymer_1`; `PolymerClass`, `Polymer`, `Molecule` are already imported at the top of the file):

```python
    def test_classify_reaction_flux_archetype(self):
        """
        classify_reaction_flux_archetype(reactants, products) drives the solver's
        per-reaction pool moment apportionment (spec 2026-06-09):
        SCISSION product -> SCISSION_FRAGMENT; single cross-pool polymer product
        -> MIGRATION; fold-backs -> SAME_POOL; no polymers -> NONE; ambiguous or
        end-initiated-scission shapes -> UNRESOLVED (legacy mu1 flux + warning).
        """
        from rmgpy.polymer import PolymerFluxArchetype, classify_reaction_flux_archetype

        p = self.polymer_1
        gas = Molecule(smiles="CC")

        # No polymer on either side -> NONE
        assert classify_reaction_flux_archetype([gas], [gas]) == PolymerFluxArchetype.NONE

        # Fold-back (same label as the reactant pool) -> SAME_POOL
        fold = p.copy()
        fold._reacted_class = PolymerClass.FEATURE
        assert (classify_reaction_flux_archetype([p], [fold, gas])
                == PolymerFluxArchetype.SAME_POOL)

        # Single cross-pool polymer product -> MIGRATION
        other = p.copy()
        other.label = "other_pool"
        other._reacted_class = PolymerClass.FEATURE
        assert (classify_reaction_flux_archetype([p], [other, gas])
                == PolymerFluxArchetype.MIGRATION)

        # SCISSION-stamped product -> SCISSION_FRAGMENT
        sc = p.copy()
        sc.label = f"{p.label}_scission_head"
        sc._reacted_class = PolymerClass.SCISSION
        assert (classify_reaction_flux_archetype([p], [sc, gas])
                == PolymerFluxArchetype.SCISSION_FRAGMENT)

        # End-initiated scission (SCISSION + END_MOD product) -> UNRESOLVED:
        # the uniform-cut bundle assumptions don't hold near a chain end.
        end = p.copy()
        end._reacted_class = PolymerClass.END_MOD
        assert (classify_reaction_flux_archetype([p], [sc, end])
                == PolymerFluxArchetype.UNRESOLVED)

        # Two polymer products with a cross-pool member -> UNRESOLVED
        assert (classify_reaction_flux_archetype([p], [other, fold])
                == PolymerFluxArchetype.UNRESOLVED)

        # Inter-chain (two reactant pools, each product folds back) -> SAME_POOL
        q = p.copy()
        q.label = "second_pool"
        fold_q = q.copy()
        fold_q._reacted_class = PolymerClass.FEATURE
        assert (classify_reaction_flux_archetype([p, q], [fold, fold_q])
                == PolymerFluxArchetype.SAME_POOL)

        # Cross-pool product with TWO reactant pools (ambiguous source) -> UNRESOLVED
        assert (classify_reaction_flux_archetype([p, q], [other])
                == PolymerFluxArchetype.UNRESOLVED)
```

- [ ] **Step 2: Run to verify failure**

Run: `python -m pytest test/rmgpy/polymerTest.py -q -k test_classify_reaction_flux_archetype`
Expected: FAIL with `ImportError: cannot import name 'PolymerFluxArchetype'`

- [ ] **Step 3: Implement**

In `rmgpy/polymer.py`, change line 106 from `from enum import Enum` to `from enum import Enum, IntEnum`. Then add immediately after the `is_end_group_reaction` function (~line 1395):

```python
class PolymerFluxArchetype(IntEnum):
    """
    Per-reaction pool moment-flux archetype, stamped at generation time on
    ``Reaction.polymer_flux_archetype`` and dispatched by the polymer hybrid
    solver's residual. See
    docs/superpowers/specs/2026-06-09-proxy-moment-flux-apportionment-design.md.
    """
    NONE = 0               # no proxy involvement
    SAME_POOL = 1          # product folds back into the reactant pool (net-zero)
    MIGRATION = 2          # whole chain migrates to a different pool
    SCISSION_FRAGMENT = 3  # chain cut; fragment to daughter, complement stays
    UNRESOLVED = 4         # ambiguous/unstamped: solver applies legacy mu1 flux


_flux_archetype_warned = set()


def _warn_unresolved_archetype(reason, detail):
    """Log each distinct UNRESOLVED-archetype cause once (not per reaction)."""
    key = (reason, detail)
    if key not in _flux_archetype_warned:
        _flux_archetype_warned.add(key)
        logging.warning(
            "Polymer flux archetype UNRESOLVED (%s): %s -- the solver will "
            "apply legacy mu1-only moment flux for this reaction shape.",
            reason, detail)


def classify_reaction_flux_archetype(reactants, products) -> 'PolymerFluxArchetype':
    """
    Classify a reaction's pool moment-flux archetype from its (handshaked)
    reactant and product lists. Product Polymers carry ``_reacted_class``
    stamped by :meth:`Polymer.create_reacted_copy`; pool identity is the
    Polymer label (the same key the solver's ``initialize_model`` uses to
    map species to pools).
    """
    reactant_pools = {r.label for r in reactants if isinstance(r, Polymer)}
    product_polymers = [p for p in products if isinstance(p, Polymer)]
    if not reactant_pools and not product_polymers:
        return PolymerFluxArchetype.NONE

    if any(getattr(p, '_reacted_class', None) == PolymerClass.SCISSION
           for p in product_polymers):
        if is_end_group_reaction(products):
            # End-initiated scission: the SCISSION_FRAGMENT bundle assumes a
            # length-biased chain pick and a uniform cut position; an
            # end-group-scaled scission violates both (uniform-by-chain pick,
            # cut near the chain end) and would overestimate daughter Mn.
            _warn_unresolved_archetype(
                "end-initiated scission",
                tuple(sorted(p.label for p in product_polymers)))
            return PolymerFluxArchetype.UNRESOLVED
        return PolymerFluxArchetype.SCISSION_FRAGMENT

    cross_pool = [p for p in product_polymers if p.label not in reactant_pools]
    if not cross_pool:
        return PolymerFluxArchetype.SAME_POOL
    if (len(product_polymers) == 1 and len(cross_pool) == 1
            and len(reactant_pools) == 1):
        return PolymerFluxArchetype.MIGRATION
    _warn_unresolved_archetype(
        "ambiguous cross-pool shape",
        tuple(sorted(p.label for p in product_polymers)))
    return PolymerFluxArchetype.UNRESOLVED
```

Note: `logging` is already imported in `rmgpy/polymer.py`.

- [ ] **Step 4: Run to verify pass**

Run: `python -m pytest test/rmgpy/polymerTest.py -q -k test_classify_reaction_flux_archetype`
Expected: 1 passed

- [ ] **Step 5: Run the whole polymer unit suite (no regressions)**

Run: `python -m pytest test/rmgpy/polymerTest.py -q`
Expected: all passed (was 135+ before; +1)

- [ ] **Step 6: Commit**

```bash
git add rmgpy/polymer.py test/rmgpy/polymerTest.py
git commit -m "polymer: add reaction flux archetype classifier"
```

---

### Task 2: `Reaction.polymer_flux_archetype` field

**Files:**
- Modify: `rmgpy/reaction.pxd:66` (after `cdef public bint is_end_group_reaction`)
- Modify: `rmgpy/reaction.py:129-154` (`__init__`)
- Test: `test/rmgpy/reactionTest.py:787` (after `test_is_end_group_reaction_default_and_kwarg`)

- [ ] **Step 1: Write the failing test**

Add to `test/rmgpy/reactionTest.py` right after `test_is_end_group_reaction_default_and_kwarg` (same class):

```python
    def test_polymer_flux_archetype_default_and_kwarg(self):
        """
        Reaction carries a ``polymer_flux_archetype`` int (default 0 = NONE)
        stamped by the polymer handshake; the hybrid solver dispatches pool
        moment flux on it (values mirror rmgpy.polymer.PolymerFluxArchetype).
        """
        assert Reaction().polymer_flux_archetype == 0
        assert Reaction(polymer_flux_archetype=3).polymer_flux_archetype == 3
```

- [ ] **Step 2: Run to verify failure**

Run: `python -m pytest test/rmgpy/reactionTest.py -q -k test_polymer_flux_archetype_default_and_kwarg`
Expected: FAIL — `AttributeError: ... has no attribute 'polymer_flux_archetype'` (cdef class: undeclared attributes don't exist)

- [ ] **Step 3: Implement**

In `rmgpy/reaction.pxd`, after line 66 (`cdef public bint is_end_group_reaction`):

```pxd
    cdef public int polymer_flux_archetype
```

In `rmgpy/reaction.py` `__init__`: add the kwarg after `is_end_group_reaction=False,` (line 129):

```python
                 is_end_group_reaction=False,
                 polymer_flux_archetype=0,
```

and after the `self.is_end_group_reaction = is_end_group_reaction` assignment (line 154):

```python
        # Pool moment-flux archetype (int values of
        # rmgpy.polymer.PolymerFluxArchetype), stamped by the polymer handshake
        # in rmgpy.rmg.model.make_new_reaction. 0 = NONE. Like
        # is_end_group_reaction it is NOT serialized in __reduce__; the solver
        # remaps unstamped proxy-touching reactions to UNRESOLVED at
        # initialize_model (legacy mu1 flux), so restarts stay correct-but-loud.
        self.polymer_flux_archetype = polymer_flux_archetype
```

Do NOT modify `__reduce__` (line 201) — `is_end_group_reaction` is deliberately absent there too.

- [ ] **Step 4: Rebuild Cython (pxd changed)**

Run: `make`
Expected: compiles without error (rebuilds `rmgpy/reaction` and dependents)

- [ ] **Step 5: Run to verify pass**

Run: `python -m pytest test/rmgpy/reactionTest.py -q -k "test_polymer_flux_archetype_default_and_kwarg or test_is_end_group_reaction_default_and_kwarg"`
Expected: 2 passed

- [ ] **Step 6: Commit**

```bash
git add rmgpy/reaction.pxd rmgpy/reaction.py test/rmgpy/reactionTest.py
git commit -m "reaction: declare polymer_flux_archetype field"
```

---

### Task 3: Stamp the archetype in `make_new_reaction`

**Files:**
- Modify: `rmgpy/rmg/model.py:47` (import) and the two handshake blocks (lines ~644 and ~683)
- Test: `test/rmgpy/polymerTest.py` (the handshake test class with `self.ps` and `self._handshake`, after `test_handshake_end_mod_flags_end_group_reaction` ~line 2654)

- [ ] **Step 1: Write the failing test**

This mirrors how the µ0-scaling stamp was tested: handshake products, then assert the classifier (the exact expression `make_new_reaction` evaluates) returns the right archetype. Add to the handshake test class in `test/rmgpy/polymerTest.py` (it provides `self.ps`, `self._handshake`, `_methyl_radical_adj`, and `polymer.stitch_molecules_by_labeled_atoms` — see `test_handshake_head_scission_fragment_returns_scission_tail_polymer` in the same class for the fragment-construction pattern):

```python
    def test_handshake_products_classify_flux_archetype(self):
        """
        After the handshake, classify_reaction_flux_archetype (the expression
        make_new_reaction stamps onto Reaction.polymer_flux_archetype) returns
        SAME_POOL for fold-back products and SCISSION_FRAGMENT for scission
        fragments.
        """
        from rmgpy.polymer import PolymerFluxArchetype, classify_reaction_flux_archetype

        # Baseline fold-back -> SAME_POOL
        base = [self.ps.baseline_proxy.molecule[0].copy(deep=True)]
        self._handshake(base, [self.ps])
        assert isinstance(base[0], Polymer)
        assert (classify_reaction_flux_archetype([self.ps], base)
                == PolymerFluxArchetype.SAME_POOL)

        # Head-wing-only fragment -> scission_tail Polymer -> SCISSION_FRAGMENT
        head_wing = self.ps._stitch_wing("head")
        methyl_star2 = Molecule().from_adjacency_list(_methyl_radical_adj("*2"))
        frag = polymer.stitch_molecules_by_labeled_atoms(head_wing, methyl_star2)
        assert frag is not None
        prods = [frag]
        self._handshake(prods, [self.ps])
        assert isinstance(prods[0], Polymer)
        assert (classify_reaction_flux_archetype([self.ps], prods)
                == PolymerFluxArchetype.SCISSION_FRAGMENT)
```

- [ ] **Step 2: Run to verify it fails meaningfully**

Run: `python -m pytest test/rmgpy/polymerTest.py -q -k test_handshake_products_classify_flux_archetype`
Expected: PASS already (the classifier exists from Task 1 — this test pins the handshake-level contract). If it FAILS, the classifier or `_reacted_class` stamping has a real bug: STOP and fix before proceeding.

- [ ] **Step 3: Implement the stamping**

In `rmgpy/rmg/model.py` line 47, extend the import:

```python
from rmgpy.polymer import Polymer, PolymerCrosslinkError, is_end_group_reaction, classify_reaction_flux_archetype
```

At line ~644 (the `forward.family and forward.is_forward` branch), directly after `forward.is_end_group_reaction = is_end_group_reaction(forward.products)`:

```python
                forward.polymer_flux_archetype = int(
                    classify_reaction_flux_archetype(reactants, forward.products))
```

At line ~683 (the other branch), directly after the same stamp there:

```python
                    forward.polymer_flux_archetype = int(
                        classify_reaction_flux_archetype(reactants, forward.products))
```

- [ ] **Step 4: Run suites**

Run: `python -m pytest test/rmgpy/polymerTest.py test/rmgpy/polymerMultiPoolTest.py -q`
Expected: all passed (no count change beyond Task 1's +1 and this +1)

- [ ] **Step 5: Commit**

```bash
git add rmgpy/rmg/model.py test/rmgpy/polymerTest.py
git commit -m "model: stamp polymer_flux_archetype at handshake"
```

---

### Task 4: Solver arrays + NONE→UNRESOLVED remap (`initialize_model`)

**Files:**
- Modify: `rmgpy/solver/polymer.pyx` — module constants near line 70, and `initialize_model` (~lines 422-475)
- Test: `test/rmgpy/solver/solverPolymerTest.py`

- [ ] **Step 1: Write the failing tests**

Add to `test/rmgpy/solver/solverPolymerTest.py` inside `TestHybridPolymerReactor`:

```python
    def test_flux_archetype_constants_match_enum(self):
        """The solver's mirror constants must equal PolymerFluxArchetype."""
        from rmgpy.polymer import PolymerFluxArchetype
        import rmgpy.solver.polymer as sp
        assert sp.FLUX_NONE == int(PolymerFluxArchetype.NONE) == 0
        assert sp.FLUX_SAME_POOL == int(PolymerFluxArchetype.SAME_POOL) == 1
        assert sp.FLUX_MIGRATION == int(PolymerFluxArchetype.MIGRATION) == 2
        assert sp.FLUX_SCISSION_FRAGMENT == int(PolymerFluxArchetype.SCISSION_FRAGMENT) == 3
        assert sp.FLUX_UNRESOLVED == int(PolymerFluxArchetype.UNRESOLVED) == 4

    def test_unstamped_proxy_reaction_remaps_to_unresolved(self):
        """
        A proxy-touching reaction arriving with the default archetype 0 (NONE)
        — e.g. restored from a pickle, which does not serialize the stamp —
        must be remapped to UNRESOLVED at initialize_model so the solver
        applies legacy mu1 flux instead of silently dropping pool moment flux.
        Pure-gas reactions must stay NONE.
        """
        Proxy = _spc("CCCC", "poly")
        Mu0 = _spc("CO", "poly_mu0")
        Mu1 = _spc("C=O", "poly_mu1")
        Mu2 = _spc("C#N", "poly_mu2")
        A = _spc("C", "A")
        B = _spc("[CH3]", "B")

        core_species = [Proxy, Mu0, Mu1, Mu2, A, B]
        gas_species_mask = np.array([False, False, False, False, True, True], dtype=bool)

        kin = Arrhenius(A=(2.0, "1/s"), n=0.0, Ea=(0.0, "kcal/mol"), T0=(298.15, "K"))
        proxy_rxn = Reaction(reactants=[Proxy], products=[B], kinetics=kin, reversible=False)
        gas_rxn = Reaction(reactants=[A], products=[B], kinetics=kin, reversible=False)
        assert proxy_rxn.polymer_flux_archetype == 0  # unstamped

        pool = PolymerPoolConfig(
            label="poly", xs=2, explicit_dp_to_species_index={},
            mu_indices=(1, 2, 3), monomer_poly_index=None,
            k_scission=0.0, k_unzip=0.0, tail_kinetics=None,
        )
        rs = HybridPolymerSystem(
            T=800.0, P=1.0e5, initial_mole_fractions={A: 1.0}, V_poly=1.0,
            polymer_pools=[pool], mass_transfer=[],
            gas_species_mask=gas_species_mask.copy(), constant_gas_volume=False,
            initial_polymer_moments={"poly": (1.0, 5.0, 30.0)}, termination=[],
        )
        rs.initialize_model(core_species, [proxy_rxn, gas_rxn], [], [])

        assert rs.reaction_flux_archetype[0] == 4   # UNRESOLVED
        assert rs.reaction_src_pool[0] == 0
        assert rs.reaction_dst_pool[0] == -1        # gas-only products
        assert rs.reaction_flux_archetype[1] == 0   # pure-gas stays NONE
        assert rs.reaction_src_pool[1] == -1
```

- [ ] **Step 2: Run to verify failure**

Run: `python -m pytest test/rmgpy/solver/solverPolymerTest.py -q -k "test_flux_archetype_constants_match_enum or test_unstamped_proxy_reaction_remaps_to_unresolved"`
Expected: FAIL — `AttributeError` (no `FLUX_NONE` / no `reaction_flux_archetype`)

- [ ] **Step 3: Implement**

In `rmgpy/solver/polymer.pyx`, near the module constants (after `SMALL_EPS = 1e-30`, ~line 70):

```python
# Pool moment-flux archetypes. Mirror of rmgpy.polymer.PolymerFluxArchetype
# (not imported to avoid a solver->polymer module cycle); equality is pinned
# by test_flux_archetype_constants_match_enum.
FLUX_NONE = 0
FLUX_SAME_POOL = 1
FLUX_MIGRATION = 2
FLUX_SCISSION_FRAGMENT = 3
FLUX_UNRESOLVED = 4
```

In `initialize_model`, directly after the `is_end_group_reaction` fill loop (lines 422-425), add the stamp transfer:

```python
        # Per-reaction pool moment-flux archetype (spec 2026-06-09). Same
        # chain(core, edge) order as is_end_group_reaction so indices match
        # r_idx in the residual.
        self.reaction_flux_archetype = np.zeros(n_rxn, dtype=np.int8)
        self.reaction_src_pool = np.full(n_rxn, -1, dtype=np.int32)
        self.reaction_dst_pool = np.full(n_rxn, -1, dtype=np.int32)
        for i, rxn in enumerate(itertools.chain(core_reactions, edge_reactions)):
            self.reaction_flux_archetype[i] = int(getattr(rxn, "polymer_flux_archetype", 0))
```

(`n_rxn` here is the initialize_model local = core+edge total, line 417 — NOT the residual's core-only `n_rxn`.)

Then, AFTER the pool-mapping loop that fills `species_to_pool_indices` / `is_pool_proxy` (the `for pool_i, pool in enumerate(self.polymer_pools):` block starting at line 454 — insert after that whole loop ends), add src/dst resolution + remap:

```python
        # Resolve per-reaction source/target pools from species indices (no
        # label matching), and remap unstamped proxy-touching reactions
        # (archetype NONE, e.g. unpickled) to UNRESOLVED so legacy mu1 flux
        # applies instead of silently dropping pool moment flux.
        ir_arr = self.reactant_indices
        ip_arr = self.product_indices
        n_unstamped = 0
        for i in range(n_rxn):
            src = -1
            for slot in range(3):
                ridx = ir_arr[i, slot]
                if ridx != -1 and ridx < n_core and self.species_to_pool_indices[ridx] != -1:
                    src = self.species_to_pool_indices[ridx]
                    break
            dst = -1
            for slot in range(3):
                pidx = ip_arr[i, slot]
                if pidx != -1 and pidx < n_core:
                    pool_j = self.species_to_pool_indices[pidx]
                    if pool_j != -1:
                        if pool_j != src:
                            dst = pool_j  # prefer the cross-pool product
                            break
                        if dst == -1:
                            dst = pool_j  # fold-back fallback
            self.reaction_src_pool[i] = src
            self.reaction_dst_pool[i] = dst
            if self.reaction_flux_archetype[i] == FLUX_NONE and (src != -1 or dst != -1):
                self.reaction_flux_archetype[i] = FLUX_UNRESOLVED
                n_unstamped += 1
        if n_unstamped:
            logging.warning(
                "%d proxy-touching reactions arrived without a polymer_flux_archetype "
                "stamp; applying legacy mu1-only pool moment flux for them.",
                n_unstamped)
```

(`itertools` and `logging` are already imported in `polymer.pyx` at lines 55-56.)

- [ ] **Step 4: Rebuild and run**

Run: `make`
Then: `python -m pytest test/rmgpy/solver/solverPolymerTest.py -q`
Expected: all passed (existing 10+ tests green, +2 new)

- [ ] **Step 5: Commit**

```bash
git add rmgpy/solver/polymer.pyx test/rmgpy/solver/solverPolymerTest.py
git commit -m "solver: resolve polymer flux archetype arrays at init"
```

---

### Task 5: Gate edge-reaction fluxes out of the core ODE

**Files:**
- Modify: `rmgpy/solver/polymer.pyx` residual, lines ~942-981
- Test: `test/rmgpy/solver/solverPolymerTest.py`

Background: `simple.pyx:464-502` applies species fluxes ONLY for core reactions (`j < num_core_reactions`); edge reactions feed just the diagnostic `edge_reaction_rates`/`edge_species_rates`. The polymer residual currently applies `dn_dt` and consumption/production for edge reactions too — a divergence this task removes.

- [ ] **Step 1: Write the failing test**

```python
    def test_edge_reaction_fluxes_are_diagnostic_only(self):
        """
        Edge reactions (all-core reactants, edge product) must not perturb the
        integrated core state: dn_dt and consumption/production stay zero, while
        edge_reaction_rates / edge_species_rates carry the diagnostic flux.
        Matches simple.pyx semantics (core ODE = core reactions only).
        """
        A = _spc("C", "A")        # core gas reactant
        E = _spc("[CH3]", "E")    # edge product

        kin = Arrhenius(A=(2.0, "1/s"), n=0.0, Ea=(0.0, "kcal/mol"), T0=(298.15, "K"))
        edge_rxn = Reaction(reactants=[A], products=[E], kinetics=kin, reversible=False)

        rs = HybridPolymerSystem(
            T=800.0, P=1.0e5, initial_mole_fractions={A: 1.0}, V_poly=1.0,
            polymer_pools=[], mass_transfer=[],
            gas_species_mask=np.array([True], dtype=bool), constant_gas_volume=False,
            termination=[],
        )
        rs.initialize_model([A], [], [E], [edge_rxn])

        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]

        assert dn_dt[0] == 0.0                          # A untouched by edge rxn
        assert rs.core_species_consumption_rates[0] == 0.0
        assert rs.edge_reaction_rates[0] > 0.0          # diagnostics still flow
        assert rs.edge_species_rates[0] > 0.0
```

- [ ] **Step 2: Run to verify failure**

Run: `python -m pytest test/rmgpy/solver/solverPolymerTest.py -q -k test_edge_reaction_fluxes_are_diagnostic_only`
Expected: FAIL — `dn_dt[0]` is negative today (edge flux leaks into the ODE)

- [ ] **Step 3: Implement**

In the residual, replace the flux blocks at lines ~942-981 with a `core_rxn`-gated version. The proxy µ1 writes stay 1:1 for now (Task 6 replaces them with the dispatch); only the gating is added here:

```python
            core_rxn = r_idx < n_rxn

            # 3. Apply Fluxes to Reactants (core reactions only -- edge
            #    reactions are diagnostic-only, matching simple.pyx)
            if self.is_pool_proxy[r0]:
                proxy_activity[r0] += abs_flux
                if core_rxn:
                    dn_dt[self.polymer_pools[p0_pool_idx].mu_indices[1]] -= r_mol_s
                    self.core_species_consumption_rates[r0] += rf
                    self.core_species_production_rates[r0] += rr
            elif core_rxn:
                dn_dt[r0] -= r_mol_s

            if r1 != -1:
                if self.is_pool_proxy[r1]:
                    proxy_activity[r1] += abs_flux
                    if core_rxn:
                        dn_dt[self.polymer_pools[p1_pool_idx].mu_indices[1]] -= r_mol_s
                        self.core_species_consumption_rates[r1] += rf
                        self.core_species_production_rates[r1] += rr
                elif core_rxn:
                    dn_dt[r1] -= r_mol_s

            if r2 != -1 and core_rxn:
                # if you ever allow polymer in r2, mirror the same logic; otherwise keep as-is
                dn_dt[r2] -= r_mol_s

            # 4. Apply Fluxes to Products
            for p_slot in range(3):
                p_idx_tmp = ip[r_idx, p_slot]
                if p_idx_tmp == -1:
                    continue

                if p_idx_tmp < n_core:
                    if self.is_pool_proxy[p_idx_tmp]:
                        proxy_activity[p_idx_tmp] += abs_flux
                        if core_rxn:
                            pool_idx = self.species_to_pool_indices[p_idx_tmp]
                            dn_dt[self.polymer_pools[pool_idx].mu_indices[1]] += r_mol_s
                            self.core_species_production_rates[p_idx_tmp] += rf
                            self.core_species_consumption_rates[p_idx_tmp] += rr
                    elif core_rxn:
                        # feature polymer species (e.g., PS_rad, scission oligomer, etc.) stays explicit
                        dn_dt[p_idx_tmp] += r_mol_s
                else:
                    self.edge_species_rates[p_idx_tmp - n_core] += rate
```

(Note: in the residual, `n_rxn` is the CORE reaction count — line 743.)

- [ ] **Step 4: Rebuild and run the full solver suite**

Run: `make && python -m pytest test/rmgpy/solver/solverPolymerTest.py -q`
Expected: all passed. If any pre-existing test fails, inspect whether it asserted edge-flux leakage into `dn_dt` (encoding the old divergence) — if so update it citing simple.pyx:464-502 semantics; if it asserted anything else, STOP, the gating broke something real.

- [ ] **Step 5: Commit**

```bash
git add rmgpy/solver/polymer.pyx test/rmgpy/solver/solverPolymerTest.py
git commit -m "solver: gate edge-reaction fluxes out of core dn_dt"
```

---

### Task 6: Archetype dispatch in the residual (the core change)

**Files:**
- Modify: `rmgpy/solver/polymer.pyx` — new `_chain_bundle` method on `HybridPolymerSystem` + replace the proxy µ1 writes from Task 5 with the dispatch
- Test: `test/rmgpy/solver/solverPolymerTest.py`

- [ ] **Step 1: Write the failing tests**

A shared two-pool fixture is used by several tests. The reaction must be built from the SAME species objects as the core list (`initialize_model` maps reaction slots by species object identity), so the fixture exposes the species dict and takes the reaction as an argument. Add as module-level helpers after `_spc` in `test/rmgpy/solver/solverPolymerTest.py`:

```python
def _two_pool_species():
    """Species + mask for a two-pool system: A (mu 1-3), B (mu 5-7), gas G at 8."""
    sp = {
        "A": _spc("CCCC", "A"),
        "A_mu0": _spc("CO", "A_mu0"), "A_mu1": _spc("C=O", "A_mu1"), "A_mu2": _spc("C#N", "A_mu2"),
        "B": _spc("CCCCC", "B"),
        "B_mu0": _spc("CCO", "B_mu0"), "B_mu1": _spc("CC=O", "B_mu1"), "B_mu2": _spc("CC#N", "B_mu2"),
        "G": _spc("[CH3]", "G"),
    }
    core = [sp["A"], sp["A_mu0"], sp["A_mu1"], sp["A_mu2"],
            sp["B"], sp["B_mu0"], sp["B_mu1"], sp["B_mu2"], sp["G"]]
    mask = np.array([False] * 8 + [True], dtype=bool)
    return sp, core, mask


def _two_pool_rs(rxn, core, mask, mom_a, mom_b):
    pool_a = PolymerPoolConfig(label="A", xs=2, explicit_dp_to_species_index={},
                               mu_indices=(1, 2, 3), monomer_poly_index=None,
                               k_scission=0.0, k_unzip=0.0, tail_kinetics=None)
    pool_b = PolymerPoolConfig(label="B", xs=2, explicit_dp_to_species_index={},
                               mu_indices=(5, 6, 7), monomer_poly_index=None,
                               k_scission=0.0, k_unzip=0.0, tail_kinetics=None)
    rs = HybridPolymerSystem(
        T=800.0, P=1.0e5, initial_mole_fractions={core[8]: 0.0}, V_poly=1.0,
        polymer_pools=[pool_a, pool_b], mass_transfer=[],
        gas_species_mask=mask.copy(), constant_gas_volume=False,
        initial_polymer_moments={"A": mom_a, "B": mom_b}, termination=[],
    )
    rs.initialize_model(core, [rxn], [], [])
    return rs

_KIN = dict(kinetics=Arrhenius(A=(2.0, "1/s"), n=0.0, Ea=(0.0, "kcal/mol"), T0=(298.15, "K")),
            reversible=False)
```

Now the tests (in `TestHybridPolymerReactor`):

```python
    def test_migration_moves_whole_chain_bundle(self):
        """
        MIGRATION (archetype 2): one event moves a whole length-biased chain
        from pool A to pool B: bundle (1, mu2/mu1, mu3/mu1) with the
        log-Lagrange closure mu3 = mu0*(mu2/mu1)^3. A loses exactly what B
        gains (conservation by construction).
        """
        sp, core, mask = _two_pool_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["B"]], **_KIN)
        rxn.polymer_flux_archetype = 2
        mu0a, mu1a, mu2a = 1.0, 5.0, 30.0
        rs = _two_pool_rs(rxn, core, mask, (mu0a, mu1a, mu2a), (2.0, 4.0, 10.0))

        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]

        kf = rxn.get_rate_coefficient(800.0, 1.0e5)
        r = kf * mu1a                       # site-scaled by A's mu1, V_poly=1
        mu3a = mu0a * (mu2a / mu1a) ** 3    # 216.0
        b1 = mu2a / mu1a                    # 6.0
        b2 = mu3a / mu1a                    # 43.2

        assert np.isclose(dn_dt[1], -r * 1.0)        # A mu0
        assert np.isclose(dn_dt[2], -r * b1)         # A mu1
        assert np.isclose(dn_dt[3], -r * b2)         # A mu2
        assert np.isclose(dn_dt[5], +r * 1.0)        # B mu0
        assert np.isclose(dn_dt[6], +r * b1)         # B mu1
        assert np.isclose(dn_dt[7], +r * b2)         # B mu2
        # exact pairwise conservation
        assert np.isclose(dn_dt[1] + dn_dt[5], 0.0, atol=1e-14)
        assert np.isclose(dn_dt[2] + dn_dt[6], 0.0, atol=1e-14)
        assert np.isclose(dn_dt[3] + dn_dt[7], 0.0, atol=1e-14)

    def test_migration_reverse_leg_uses_target_pool_statistics(self):
        """
        Per-direction MIGRATION bundles: the reverse (rr) leg must move
        B-statistics chains B->A, not A-statistics. Pool stats are chosen
        distinguishable (b_A=(1,6,43.2) vs b_B=(1,2.5,7.8125)) so a
        wrong-source bundle fails loudly. Also pins the rf/rr volume-factor
        identity: at b0=1 on both legs, the net mu0 flux must equal the
        legacy net molar rate (rf-rr)*V_rxn.
        """
        sp, core, mask = _two_pool_species()
        # reversible=False so generate_rate_coefficients needs no thermo
        # (kb=0); the reverse leg is then driven by overriding rs.kb directly.
        rxn = Reaction(reactants=[sp["A"]], products=[sp["B"]], **_KIN)
        rxn.polymer_flux_archetype = 2
        mu_a = (1.0, 5.0, 30.0)
        mu_b = (2.0, 4.0, 10.0)
        rs = _two_pool_rs(rxn, core, mask, mu_a, mu_b)
        rs.kb[0] = 0.6

        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]

        kf = rxn.get_rate_coefficient(800.0, 1.0e5)
        rf = kf * mu_a[1]               # 2.0 * 5 = 10 (site-scaled by A mu1)
        rr = 0.6 * mu_a[1]              # kb * C(proxyB)=1, then *= site -> 3.0
        mu3_a = mu_a[0] * (mu_a[2] / mu_a[1]) ** 3   # 216.0
        mu3_b = mu_b[0] * (mu_b[2] / mu_b[1]) ** 3   # 31.25
        bA1, bA2 = mu_a[2] / mu_a[1], mu3_a / mu_a[1]    # 6.0, 43.2
        bB1, bB2 = mu_b[2] / mu_b[1], mu3_b / mu_b[1]    # 2.5, 7.8125

        assert np.isclose(dn_dt[1], -(rf - rr))              # == legacy net (b0=1)
        assert np.isclose(dn_dt[2], -rf * bA1 + rr * bB1)    # -52.5
        assert np.isclose(dn_dt[3], -rf * bA2 + rr * bB2)    # -408.5625
        assert np.isclose(dn_dt[5], +(rf - rr))
        assert np.isclose(dn_dt[6], +rf * bA1 - rr * bB1)
        assert np.isclose(dn_dt[7], +rf * bA2 - rr * bB2)

    def test_end_group_migration_uses_uniform_chain_bundle(self):
        """
        A mu0-scaled (is_end_group_reaction) MIGRATION picks chains uniformly:
        bundle (1, mu1/mu0, mu2/mu0), and the rate itself scales by mu0.
        """
        sp, core, mask = _two_pool_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["B"]], **_KIN)
        rxn.polymer_flux_archetype = 2
        rxn.is_end_group_reaction = True
        mu0a, mu1a, mu2a = 1.0, 5.0, 30.0
        rs = _two_pool_rs(rxn, core, mask, (mu0a, mu1a, mu2a), (2.0, 4.0, 10.0))

        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]

        kf = rxn.get_rate_coefficient(800.0, 1.0e5)
        r = kf * mu0a                    # site-scaled by mu0, not mu1
        assert np.isclose(dn_dt[1], -r * 1.0)
        assert np.isclose(dn_dt[2], -r * mu1a / mu0a)        # -10
        assert np.isclose(dn_dt[3], -r * mu2a / mu0a)        # -60
        assert np.isclose(dn_dt[5], +r * 1.0)
        assert np.isclose(dn_dt[6], +r * mu1a / mu0a)
        assert np.isclose(dn_dt[7], +r * mu2a / mu0a)

    def test_scission_fragment_complement_stays_in_parent(self):
        """
        SCISSION_FRAGMENT (archetype 3), complement-stays accounting: parent
        (0, -r*mu2/(2 mu1), -r*(2/3)*mu3/mu1); daughter (+r, +r*mu2/(2 mu1),
        +r*mu3/(3 mu1)). mu1 conserves exactly; total mu0 +r per event; total
        mu2 drops by r*mu3/(3 mu1).
        """
        sp, core, mask = _two_pool_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["B"], sp["G"]], **_KIN)
        rxn.polymer_flux_archetype = 3
        mu0a, mu1a, mu2a = 1.0, 5.0, 30.0
        rs = _two_pool_rs(rxn, core, mask, (mu0a, mu1a, mu2a), (0.1, 0.2, 0.5))

        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]

        kf = rxn.get_rate_coefficient(800.0, 1.0e5)
        r = kf * mu1a
        mu3a = mu0a * (mu2a / mu1a) ** 3    # 216.0
        e_n = mu2a / mu1a                   # 6.0
        e_n2 = mu3a / mu1a                  # 43.2

        assert np.isclose(dn_dt[1], 0.0, atol=1e-14)             # parent mu0
        assert np.isclose(dn_dt[2], -r * e_n / 2.0)              # parent mu1
        assert np.isclose(dn_dt[3], -r * (2.0 / 3.0) * e_n2)     # parent mu2
        assert np.isclose(dn_dt[5], +r)                          # daughter mu0
        assert np.isclose(dn_dt[6], +r * e_n / 2.0)              # daughter mu1
        assert np.isclose(dn_dt[7], +r * e_n2 / 3.0)             # daughter mu2
        assert np.isclose(dn_dt[2] + dn_dt[6], 0.0, atol=1e-14)  # mu1 conserved
        assert np.isclose(dn_dt[3] + dn_dt[7], -r * e_n2 / 3.0)  # mu2 destroyed
        assert np.isclose(dn_dt[8], +r)                          # gas co-product G

    def test_same_pool_reaction_leaves_moments_untouched(self):
        """
        SAME_POOL (archetype 1): fold-back means net-zero pool moment flux —
        the dispatch skips pool writes entirely. Gas co-species still flow.
        """
        sp, core, mask = _two_pool_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["A"], sp["G"]], **_KIN)
        rxn.polymer_flux_archetype = 1
        rs = _two_pool_rs(rxn, core, mask, (1.0, 5.0, 30.0), (2.0, 4.0, 10.0))

        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]

        kf = rxn.get_rate_coefficient(800.0, 1.0e5)
        assert np.allclose(dn_dt[1:4], 0.0, atol=1e-14)   # A moments untouched
        assert np.allclose(dn_dt[5:8], 0.0, atol=1e-14)   # B untouched
        assert np.isclose(dn_dt[8], kf * 5.0)             # gas G produced

    def test_unresolved_applies_legacy_mu1_flux(self):
        """UNRESOLVED (archetype 4) replicates the pre-apportionment behavior:
        whole event flux on mu1 only (reactant pool -r, product pool +r)."""
        sp, core, mask = _two_pool_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["B"]], **_KIN)
        rxn.polymer_flux_archetype = 4
        rs = _two_pool_rs(rxn, core, mask, (1.0, 5.0, 30.0), (2.0, 4.0, 10.0))

        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]

        kf = rxn.get_rate_coefficient(800.0, 1.0e5)
        r = kf * 5.0
        assert np.isclose(dn_dt[2], -r)                   # A mu1 only
        assert np.isclose(dn_dt[6], +r)                   # B mu1 only
        assert np.allclose([dn_dt[1], dn_dt[3], dn_dt[5], dn_dt[7]], 0.0, atol=1e-14)

    def test_bundle_guard_empty_source_pool_no_nan(self):
        """
        An empty source pool (all moments zero) must produce zero flux and no
        NaN: the site scaling already zeroes the rate, and the bundle guard
        prevents the 0/0 in mu2/mu1.
        """
        sp, core, mask = _two_pool_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["B"]], **_KIN)
        rxn.polymer_flux_archetype = 2
        rs = _two_pool_rs(rxn, core, mask, (0.0, 0.0, 0.0), (2.0, 4.0, 10.0))

        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
        assert np.all(np.isfinite(dn_dt))
        assert np.allclose(dn_dt[1:8], 0.0, atol=1e-14)

    def test_bundle_guard_mu3_overflow_skips_mu2_component_only(self):
        """
        When the mu3 closure overflows to inf (huge mu2), MIGRATION still
        applies the mu0/mu1 bundle components and skips ONLY mu2 (mirrors the
        scission-ODE finiteness gate at polymer.pyx ~1052).
        """
        sp, core, mask = _two_pool_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["B"]], **_KIN)
        rxn.polymer_flux_archetype = 2
        mu0a, mu1a, mu2a = 1.0, 5.0, 1.0e120
        rs = _two_pool_rs(rxn, core, mask, (mu0a, mu1a, mu2a), (2.0, 4.0, 10.0))

        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
        kf = rxn.get_rate_coefficient(800.0, 1.0e5)
        r = kf * mu1a
        assert np.all(np.isfinite(dn_dt))
        assert np.isclose(dn_dt[1], -r)                       # mu0 applied
        assert np.isclose(dn_dt[2], -r * mu2a / mu1a)         # mu1 applied
        assert dn_dt[3] == 0.0                                # mu2 skipped
        assert dn_dt[7] == 0.0
```

- [ ] **Step 2: Run to verify failures**

Run: `python -m pytest test/rmgpy/solver/solverPolymerTest.py -q -k "migration or scission_fragment_complement or same_pool_reaction_leaves or unresolved_applies or bundle_guard"`
Expected: 8 FAIL (today everything lands µ1-only: migration test sees `dn_dt[1] == 0`, etc.)

- [ ] **Step 3: Implement `_chain_bundle`**

Add as a method on `HybridPolymerSystem` (plain `def`, near `residual`):

```python
    def _chain_bundle(self, int pool_idx, y, double V_poly, bint end_group):
        """
        Expected (chains, units, units^2) carried by ONE picked chain of pool
        ``pool_idx``: (1, E[k], E[k^2]). end_group=True -> uniform chain pick
        (rate ~ mu0): (1, mu1/mu0, mu2/mu0). Otherwise length-biased pick
        (rate ~ mu1): (1, mu2/mu1, mu3/mu1) with the guarded mu3 closure.
        Returns (b0, b1, b2, mu2_ok); b0 == 0.0 means the pool is too empty to
        move a chain (denominator below SMALL_EPS) -- caller skips the term.
        mu2_ok False means apply b0/b1 but skip the mu2 component (mu3 = inf).
        """
        idx0, idx1, idx2 = self.polymer_pools[pool_idx].mu_indices
        mu0 = max(0.0, y[idx0]) / V_poly
        mu1 = max(0.0, y[idx1]) / V_poly
        mu2 = max(0.0, y[idx2]) / V_poly
        if end_group:
            if mu0 <= SMALL_EPS:
                return 0.0, 0.0, 0.0, False
            return 1.0, mu1 / mu0, mu2 / mu0, True
        if mu1 <= SMALL_EPS:
            return 0.0, 0.0, 0.0, False
        mu3 = _safe_mu3_from_mu012(mu0, mu1, mu2)
        if np.isfinite(mu3):
            return 1.0, mu2 / mu1, mu3 / mu1, True
        return 1.0, mu2 / mu1, 0.0, False
```

- [ ] **Step 4: Implement the dispatch**

In the residual, replace the three gated proxy µ1 writes from Task 5 (the `dn_dt[...mu_indices[1]...]` lines inside the reactant/product blocks — remove ONLY those `dn_dt` lines; `proxy_activity` and consumption/production stay) and append the dispatch block after the product loop.

**Do NOT touch the pool-index computations** at lines ~880-882 (`p0_pool_idx = self.species_to_pool_indices[r0]` etc., marked "MOVED UP before _C") — they feed the site-scaling block AND the UNRESOLVED legacy branch below. Removing them breaks both.

The reactant blocks become:

```python
            if self.is_pool_proxy[r0]:
                proxy_activity[r0] += abs_flux
                if core_rxn:
                    self.core_species_consumption_rates[r0] += rf
                    self.core_species_production_rates[r0] += rr
            elif core_rxn:
                dn_dt[r0] -= r_mol_s

            if r1 != -1:
                if self.is_pool_proxy[r1]:
                    proxy_activity[r1] += abs_flux
                    if core_rxn:
                        self.core_species_consumption_rates[r1] += rf
                        self.core_species_production_rates[r1] += rr
                elif core_rxn:
                    dn_dt[r1] -= r_mol_s

            if r2 != -1 and core_rxn:
                # if you ever allow polymer in r2, mirror the same logic; otherwise keep as-is
                dn_dt[r2] -= r_mol_s
```

The product loop becomes:

```python
            for p_slot in range(3):
                p_idx_tmp = ip[r_idx, p_slot]
                if p_idx_tmp == -1:
                    continue

                if p_idx_tmp < n_core:
                    if self.is_pool_proxy[p_idx_tmp]:
                        proxy_activity[p_idx_tmp] += abs_flux
                        if core_rxn:
                            self.core_species_production_rates[p_idx_tmp] += rf
                            self.core_species_consumption_rates[p_idx_tmp] += rr
                    elif core_rxn:
                        # feature polymer species (e.g., PS_rad, scission oligomer, etc.) stays explicit
                        dn_dt[p_idx_tmp] += r_mol_s
                else:
                    self.edge_species_rates[p_idx_tmp - n_core] += rate
```

Then append the dispatch (same indentation level, still inside the reaction loop):

```python
            # 5. Pool moment flux -- archetype dispatch (core reactions only).
            #    Spec: docs/superpowers/specs/2026-06-09-proxy-moment-flux-
            #    apportionment-design.md. SAME_POOL and NONE apply nothing:
            #    fold-back flux is net-zero by construction, so skipping avoids
            #    roundoff and closure calls.
            if core_rxn:
                arch = self.reaction_flux_archetype[r_idx]
                if arch == FLUX_MIGRATION:
                    src = self.reaction_src_pool[r_idx]
                    dst = self.reaction_dst_pool[r_idx]
                    if src != -1 and dst != -1 and src != dst:
                        # Per-direction bundles: forward moves A-statistics
                        # chains A->B, reverse moves B-statistics chains B->A.
                        # Each direction is guarded by its OWN source pool.
                        for ev_rate, from_pool, to_pool in (
                                (rf, src, dst), (rr, dst, src)):
                            if ev_rate <= 0.0:
                                continue
                            ev_mol = ev_rate * V_rxn
                            b0, b1, b2, mu2_ok = self._chain_bundle(
                                from_pool, y, V_poly,
                                self.is_end_group_reaction[r_idx])
                            if b0 == 0.0:
                                continue
                            f_idx = self.polymer_pools[from_pool].mu_indices
                            t_idx = self.polymer_pools[to_pool].mu_indices
                            dn_dt[f_idx[0]] -= ev_mol * b0
                            dn_dt[f_idx[1]] -= ev_mol * b1
                            dn_dt[t_idx[0]] += ev_mol * b0
                            dn_dt[t_idx[1]] += ev_mol * b1
                            if mu2_ok:
                                dn_dt[f_idx[2]] -= ev_mol * b2
                                dn_dt[t_idx[2]] += ev_mol * b2
                elif arch == FLUX_SCISSION_FRAGMENT:
                    src = self.reaction_src_pool[r_idx]
                    dst = self.reaction_dst_pool[r_idx]
                    if src != -1 and dst != -1 and src != dst:
                        s_idx = self.polymer_pools[src].mu_indices
                        d_idx = self.polymer_pools[dst].mu_indices
                        mu0_p = max(0.0, y[s_idx[0]]) / V_poly
                        mu1_p = max(0.0, y[s_idx[1]]) / V_poly
                        mu2_p = max(0.0, y[s_idx[2]]) / V_poly
                        ok = mu1_p > SMALL_EPS
                        if ok and r_mol_s < 0.0:
                            # Net reverse = coupling bookkeeping; it depletes
                            # the DAUGHTER, so guard its moments too. (Stated
                            # approximation: parent statistics, sign-flipped.)
                            if (max(0.0, y[d_idx[0]]) / V_poly <= SMALL_EPS or
                                    max(0.0, y[d_idx[1]]) / V_poly <= SMALL_EPS):
                                ok = False
                        if ok:
                            # Complement-stays-in-parent accounting: parent
                            # mu0 net 0 (chain broke, complement remains);
                            # fragment (uniform cut of a length-biased chain:
                            # E[a] = E[n]/2, E[a^2] = E[n^2]/3) moves to the
                            # daughter. mu1 conserves exactly per reaction.
                            e_n = mu2_p / mu1_p
                            dn_dt[s_idx[1]] -= r_mol_s * e_n / 2.0
                            dn_dt[d_idx[0]] += r_mol_s
                            dn_dt[d_idx[1]] += r_mol_s * e_n / 2.0
                            mu3_p = _safe_mu3_from_mu012(mu0_p, mu1_p, mu2_p)
                            if np.isfinite(mu3_p):
                                e_n2 = mu3_p / mu1_p
                                dn_dt[s_idx[2]] -= r_mol_s * (2.0 / 3.0) * e_n2
                                dn_dt[d_idx[2]] += r_mol_s * e_n2 / 3.0
                elif arch == FLUX_UNRESOLVED:
                    # Legacy mu1-only transfer (pre-apportionment behavior),
                    # replicated exactly: -r per reactant proxy, +r per
                    # product proxy.
                    if self.is_pool_proxy[r0]:
                        dn_dt[self.polymer_pools[p0_pool_idx].mu_indices[1]] -= r_mol_s
                    if r1 != -1 and self.is_pool_proxy[r1]:
                        dn_dt[self.polymer_pools[p1_pool_idx].mu_indices[1]] -= r_mol_s
                    for p_slot in range(3):
                        p_idx_tmp = ip[r_idx, p_slot]
                        if (p_idx_tmp != -1 and p_idx_tmp < n_core
                                and self.is_pool_proxy[p_idx_tmp]):
                            pool_idx = self.species_to_pool_indices[p_idx_tmp]
                            dn_dt[self.polymer_pools[pool_idx].mu_indices[1]] += r_mol_s
```

- [ ] **Step 5: Rebuild and run the new tests**

Run: `make && python -m pytest test/rmgpy/solver/solverPolymerTest.py -q -k "migration or scission_fragment_complement or same_pool_reaction_leaves or unresolved_applies or bundle_guard"`
Expected: 8 passed

- [ ] **Step 6: Verify the rate-diagnostic block needs no change**

Read `rmgpy/solver/polymer.pyx` around line 1269 (`[THE HIJACK]` in `get_reaction_rates`): it computes reaction *rates* only (µ0/µ1 site scaling, no `dn_dt` writes), so the apportionment does not touch it. Confirm it is unchanged by your edits (`git diff rmgpy/solver/polymer.pyx` should show nothing in `get_reaction_rates`).

- [ ] **Step 6b: Verify the Jacobian path is numeric**

The dispatch adds cross-pool µ0/µ2 couplings and the nonlinear µ3 closure to the RHS; an analytic Jacobian would now be wrong. Confirm DASPK uses a numeric finite-difference Jacobian for this system (pydas sets the analytic-Jacobian flag iff the instance has a `jacobian` attribute):

```bash
~/anaconda3/envs/rmg_env/bin/python -c "
from rmgpy.solver.polymer import HybridPolymerSystem
assert not hasattr(HybridPolymerSystem, 'jacobian'), 'BLOCKER: analytic jacobian present'
print('numeric Jacobian confirmed')"
```

Expected: `numeric Jacobian confirmed` (verified 2026-06-09 — the analytic `jacobian()` belongs to `SimpleReactor` only). If this ever fails, STOP: the analytic Jacobian must be updated or disabled before the dispatch can land.

- [ ] **Step 7: Run the full four-suite set (regression)**

Run: `python -m pytest test/rmgpy/polymerTest.py test/rmgpy/solver/solverPolymerTest.py test/rmgpy/canteraTest.py test/rmgpy/polymerMultiPoolTest.py -q`
Expected: ALL passed (synthetic decks built with raw `Reaction()` objects arrive unstamped → remap to UNRESOLVED → legacy flux, so pre-existing tests keep their old numbers). Investigate ANY failure before proceeding — only tests that directly asserted µ1-only flux on a *stamped* reaction may be updated, and there should be none.

- [ ] **Step 8: Commit**

```bash
git add rmgpy/solver/polymer.pyx test/rmgpy/solver/solverPolymerTest.py
git commit -m "solver: apportion pool moment flux by reaction archetype"
```

---

### Task 7: Trajectory, conservation, monodisperse, and clamp tests

**Files:**
- Test: `test/rmgpy/solver/solverPolymerTest.py` (uses the Task 6 fixtures)

These pin the spec's verification requirements beyond single RHS evaluations. No production code should change in this task — if a test fails, that's a real bug in Task 6: STOP and fix there.

- [ ] **Step 1: Write the monodisperse-limit test**

For a monodisperse pool (µ0, µ1, µ2) = (N, N·k, N·k²) the log-Lagrange closure is EXACT (µ3 = N·k³), so the scission table has sharp closed forms: E[n] = k, E[n²] = k².

```python
    def test_scission_monodisperse_limit_closed_form(self):
        """
        Monodisperse pool (PDI=1, mu_j = N*k^j): the mu3 closure is exact
        (mu3 = N*k^3), so SCISSION_FRAGMENT has closed-form derivatives:
        parent (0, -r*k/2, -r*(2/3)*k^2), daughter (+r, +r*k/2, +r*k^2/3).
        """
        sp, core, mask = _two_pool_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["B"], sp["G"]], **_KIN)
        rxn.polymer_flux_archetype = 3
        N, k = 2.0, 5.0
        rs = _two_pool_rs(rxn, core, mask, (N, N * k, N * k * k), (0.1, 0.2, 0.5))

        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
        kf = rxn.get_rate_coefficient(800.0, 1.0e5)
        r = kf * (N * k)                                  # site-scaled by mu1

        assert np.isclose(dn_dt[1], 0.0, atol=1e-14)
        assert np.isclose(dn_dt[2], -r * k / 2.0)
        assert np.isclose(dn_dt[3], -r * (2.0 / 3.0) * k * k)
        assert np.isclose(dn_dt[5], +r)
        assert np.isclose(dn_dt[6], +r * k / 2.0)
        assert np.isclose(dn_dt[7], +r * k * k / 3.0)
```

- [ ] **Step 2: Write the conservation + realizability trajectory test**

```python
    def test_apportionment_trajectory_conserves_mu1_and_stays_realizable(self):
        """
        Forward-Euler trajectory with a SCISSION_FRAGMENT reaction A -> B + G:
        total monomer units (mu1_A + mu1_B) stay constant (the gas co-product
        G tracks events, not units), total chain count mu0_A + mu0_B grows,
        and both pools stay in the realizable cone (mu1 >= mu0 >= 0).
        """
        sp, core, mask = _two_pool_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["B"], sp["G"]], **_KIN)
        rxn.polymer_flux_archetype = 3
        rs = _two_pool_rs(rxn, core, mask, (1.0, 50.0, 3000.0), (0.0, 0.0, 0.0))

        y = rs.y.copy()
        mu1_total0 = y[2] + y[6]
        mu0_total0 = y[1] + y[5]
        # Step sizing: parent mu1 drains at ~kf*mu2/2 = 3000/s initially, so
        # keep t_total = 2e-3 s (~12% of parent mu1) to stay far from the
        # depletion overshoot regime that forward Euler handles poorly.
        dt = 1e-5
        for _ in range(200):
            dn_dt = rs.residual(0.0, y, np.zeros_like(y))[0]
            assert np.all(np.isfinite(dn_dt))
            y = y + dt * dn_dt
            for i0, i1 in ((1, 2), (5, 6)):
                assert y[i1] >= -1e-9                  # mu1 >= 0
                assert y[i1] - y[i0] >= -1e-6          # mu1 >= mu0 (cone)
        assert np.isclose(y[2] + y[6], mu1_total0, rtol=1e-9)   # units conserved
        assert y[1] + y[5] > mu0_total0                          # chains created
```

- [ ] **Step 3: Write the closure-overestimate clamp test**

```python
    def test_scission_mu2_overdrain_stays_finite_and_warns(self, caplog):
        """
        A high-PDI parent makes the closure-estimated mu3 huge; the resulting
        parent mu2 drain can overshoot in an explicit step. The residual must
        stay finite when evaluated at the resulting negative-mu2 state (reads
        are max(0,.)-clamped), and debug_check_realizability must log the cone
        violation rather than raise.
        """
        import logging as _logging
        sp, core, mask = _two_pool_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["B"], sp["G"]], **_KIN)
        rxn.polymer_flux_archetype = 3
        rs = _two_pool_rs(rxn, core, mask, (1.0, 5.0, 1000.0), (0.1, 0.2, 0.5))
        rs.debug_check_realizability = True

        y = rs.y.copy()
        dn_dt = rs.residual(0.0, y, np.zeros_like(y))[0]
        assert np.all(np.isfinite(dn_dt))
        # Overshoot mu2 negative with one huge explicit step, then re-evaluate.
        y2 = y + 1e3 * dn_dt
        assert y2[3] < 0.0, "test setup: expected a mu2 overshoot"
        with caplog.at_level(_logging.WARNING):
            dn_dt2 = rs.residual(0.0, y2, np.zeros_like(y2))[0]
        assert np.all(np.isfinite(dn_dt2))
        assert any("realizable cone" in rec.message for rec in caplog.records)
```

(The asserted substring matches the actual `debug_check_realizability` message — "Polymer pool '%s' moment state left the realizable cone ..." at `polymer.pyx:1023-1027`, verified 2026-06-09.)

- [ ] **Step 4: Run the three tests**

Run: `python -m pytest test/rmgpy/solver/solverPolymerTest.py -q -k "monodisperse_limit or trajectory_conserves or mu2_overdrain"`
Expected: 3 passed. (If `test setup` assertions fail, tune the magic numbers — dt, steps, the 1e3 step factor — NOT the production code.)

- [ ] **Step 5: Commit**

```bash
git add test/rmgpy/solver/solverPolymerTest.py
git commit -m "solver: trajectory tests for moment apportionment"
```

---

### Task 8: Docs, full verification, end-to-end

**Files:**
- Modify: `docs/multi_pool_design.md` §5 (the section that specifies the scission/unzip moment ODEs)

- [ ] **Step 1: Document the apportionment in the design doc**

Append to `docs/multi_pool_design.md` §5 (after the scission/unzip ODE specification):

```markdown
### 5.x Discrete-reaction moment apportionment (flux archetypes)

Every proxy-touching reaction is stamped at generation time with a flux
archetype (`rmgpy.polymer.PolymerFluxArchetype`, carried on
`Reaction.polymer_flux_archetype` and resolved to solver arrays at
`initialize_model`). The residual dispatches pool moment flux on it:

| Archetype | Pool moment flux per net event rate r |
|---|---|
| SAME_POOL | none (fold-back is net-zero by construction) |
| MIGRATION | whole-chain bundle (1, E[k], E[k²]) per direction: forward rf moves source-pool statistics, reverse rr moves target-pool statistics. µ1-scaled reactions pick length-biased chains (E[k]=µ2/µ1, E[k²]=µ3/µ1, guarded closure); µ0-scaled end-group reactions pick uniformly (µ1/µ0, µ2/µ0). |
| SCISSION_FRAGMENT | complement-stays-in-parent: parent (0, −r·µ2/(2µ1), −r·(2/3)·µ3/µ1), daughter (+r, +r·µ2/(2µ1), +r·µ3/(3µ1)). Conserves µ1 exactly and adds one chain per event, independent of whether head/tail/both reactions were generated. |
| UNRESOLVED | legacy µ1-only transfer + one-time warning (also the fallback for unstamped reactions, e.g. restored from pickles). |

Edge reactions are diagnostic-only (matching `simple.pyx`): their fluxes never
enter `dn_dt` or the consumption/production accumulators.

**Input hygiene:** the phenomenological `k_scission` coexists with
family-generated scission reactions as a chemically distinct (thermal random)
channel. A `k_scission` fitted to bulk degradation data that already includes
radical chemistry double-counts — specify the thermal-only rate.

Full derivations and decision rationale:
`docs/superpowers/specs/2026-06-09-proxy-moment-flux-apportionment-design.md`.
```

(Adapt the subsection number `5.x` to the actual next free number in the file.)

- [ ] **Step 2: Full four-suite run**

Run: `python -m pytest test/rmgpy/polymerTest.py test/rmgpy/solver/solverPolymerTest.py test/rmgpy/canteraTest.py test/rmgpy/polymerMultiPoolTest.py -q`
Expected: ~208 passed (191 + ~17 new), 0 failures

- [ ] **Step 3: Reaction + model suites (the other touched layers)**

Run: `python -m pytest test/rmgpy/reactionTest.py test/rmgpy/rmg/modelTest.py -q`
Expected: all passed

- [ ] **Step 4: EPDM end-to-end (no regression)**

Run: `cd ~/runs/RMG/epdm_v0_2026-06-06b && ~/anaconda3/envs/rmg_env/bin/python ~/Code/RMG-Py/rmg.py input.py`
Expected: `MODEL GENERATION COMPLETED`, 26 core species / 28 reactions (deck is scission-ODE dominated; the apportionment must not change its outcome). Also check the log for unexpected floods of the UNRESOLVED warning — a handful is fine (restart artifacts), per-reaction spam is not.

- [ ] **Step 5: Commit docs**

```bash
cd ~/Code/RMG-Py
git add docs/multi_pool_design.md
git commit -m "docs: moment apportionment in multi_pool_design"
```

- [ ] **Step 6: Update the auto-memory running log**

Edit `~/.claude/projects/-home-alon-Code-RMG-Py/memory/project_polymer_branch_review.md`: mark open item #4 as IMPLEMENTED with the commit range and test names, mirroring the style of the existing FIXED entries.

---

## Verification summary (spec → task map)

| Spec requirement | Task |
|---|---|
| Archetype classifier + UNRESOLVED rules (incl. end-initiated scission, multi-product, |R|>1) | 1 |
| `Reaction.polymer_flux_archetype` declared field, no `__reduce__` change | 2 |
| Stamp at both `make_new_reaction` branches | 3 |
| Solver arrays in chain(core,edge) order, src/dst from species indices, NONE→UNRESOLVED remap in `initialize_model` | 4 |
| Edge-reaction gating (`dn_dt` + consumption/production) | 5 |
| SAME_POOL skip; MIGRATION per-direction bundles (incl. reverse-leg B-statistics + rf/rr volume identity); µ0-scaled end-group bundle; SCISSION complement-stays table; UNRESOLVED legacy; guards (empty pool, µ3 finite, daughter-depletion on r<0); numeric-Jacobian check | 6 |
| Monodisperse limit, global µ1 conservation, realizability, clamp+warn | 7 |
| Docs (apportionment + k_scission caveat), EPDM end-to-end, full suites | 8 |
