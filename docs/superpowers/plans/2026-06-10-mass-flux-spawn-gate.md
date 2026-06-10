# Mass-Flux Spawn Gate Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the unconditional `0.5` stub in `_estimate_relative_flux` with the real, honest mass-flux spawn gate approved in `docs/superpowers/specs/2026-06-10-mass-flux-spawn-gate-design.md` (the normative source for every behavior below).

**REVISED 2026-06-10 against the AMENDED spec (commit `cded166bc`).** The original plan's representative-flux sourcing was proven born-dead in production (`is_pool_proxy` marks one canonical proxy per pool; gross arrays are maintained only in `is_pool_proxy` branches; representatives are ordinary species never in `species_to_pool_indices` — numerator identically zero in every real run, all fabricated-snapshot tests green). This revision implements the amendment: change (a) (gross writes for ALL core species, cost-gated, spec §4.6), the 3-tuple snapshot (spec §4.1), `(species_label, parent_pool_label)` representative pairs with a deduped denominator (spec §3/§4.3), and the RED-FIRST integrated tripwire (spec §7.1) as the FIRST executable task — its confirmed-red run on current HEAD gates the rest of the plan.

**Standing convention (spec §7 preamble — read this before implementing anything):**

> Every laundered-quantity fix ships a test that FAILS in the live path until the quantity is real. Fabricated-input tests prove formula correctness; they cannot prove the live path sources the quantity — twice now (the 0.5 gate, the born-dead representative flux) a green fabricated-input suite hid an identically-zero production value. This convention applies to the remaining spawn/seeding items (#14 `triggering_moles` seeding) without per-spec rediscovery.

**Architecture:** Change (a) makes the polymer residual maintain gross production/consumption for ALL core species (simple.pyx parity; today only `is_pool_proxy` species get them). A new engine method `HybridPolymerSystem.spawn_gate_flux_snapshot()` returns a 3-tuple — `gross` (label → max(0, production) for every core species), `pool_stats` (pool label → (E[n], monomer_MW) with the SMALL_EPS guard), and `proxy_event_mass_total` (the engine-attributed canonical-proxy event-mass sum). `main.py` stashes that snapshot on the reaction model after each `simulate()`. The python-side gate does the ledger-dependent attribution: representatives are `(species_label, parent_pool_label)` pairs recorded at absorption; `g_i = gross[label] · E[n][parent_pool] · MW[parent_pool]`; the denominator is `proxy_event_mass_total` + the DEDUPED representative event-mass across the whole ledger (each motif's fraction ∈ [0,1]; multi-motif double-counting across numerators is a stated decision). A Group-isomorphism motif ledger plus the existing `MassFluxAccumulator` (one record per motif per iteration, sum/N zero-filled statistic) on `CoreEdgeReactionModel` drives the Phase-D spawn decision. No snapshot → fraction 0.0 → defer, logged: no production code path fakes a number.

**Tech Stack:** Python 3 / Cython (`rmgpy/solver/polymer.pyx` — run `make` after editing), pytest (`~/anaconda3/envs/rmg_env/bin/python -m pytest`), RMG-Py `polymer` branch.

---

## Verified baseline facts (re-probed 2026-06-10, post-amendment)

- Full 10-file regression suite: **380 passed, 2 skipped** (43 s) — re-run and confirmed.
- `test/rmgpy/polymerMultiPoolTest.py`: 26 tests; `test/rmgpy/solver/solverPolymerTest.py`: 41 tests.
- `SMALL_EPS = 1e-30` at `rmgpy/solver/polymer.pyx:70`; the `tail_mean` guard idiom is `mu1 / mu0 if mu0 > SMALL_EPS else 0.0`.
- **Gross arrays are proxy-only today (the born-dead hole, spec §3.1):** `core_species_consumption_rates`/`core_species_production_rates` are zeroed at polymer.pyx:898-899 each residual call and written ONLY inside the `is_pool_proxy` branches of residual sections 3/4 (reactant slots r0 at :1082-1086, r1 at :1090-1095; product slots at :1109-1114). Ordinary species get `dn_dt` writes only (:1087-1088, :1096-1097, :1115-1117); the r2 reactant slot (:1099-1101) is always ordinary today and carries a "mirror the same logic" comment.
- **simple.pyx parity pattern** (rmgpy/solver/simple.pyx:472-498): for each reactant slot `consumption += f_reaction_rate; production += rev_reaction_rate`; for each product slot `production += f_reaction_rate; consumption += rev_reaction_rate`; core reactions only.
- `PolymerPoolConfig` (polymer.pyx:89-124, frozen dataclass) has **no monomer-MW field**; `monomer_poly_index: Optional[int] = None` is the last field before the `# Kinetics Parameters` comment.
- `rmgpy/rmg/polymer_input.py` has TWO `to_config` methods: `PolymerPool.to_config` at :779 is the ONLY one constructing `PolymerPoolConfig` (at :807); `MassTransfer.to_config` at :845 constructs `MassTransferConfig` and is untouched by this change. In `PolymerPool.to_config`, `self.monomer` is a resolved `Species`; the `Polymer.monomer_mw_g_mol` idiom is `self.monomer.get_molecular_weight() * 1000.0` (rmgpy/polymer.py:261).
- Engine label mapping: `self.species_index` (Species → index, core first then edge; filled by `generate_species_indices`, rmgpy/solver/base.pyx:419-426). `is_pool_proxy[i] = 1` only where `core_species[i].label.partition('(')[0] == pool.label` (polymer.pyx:478-485) — exactly ONE canonical proxy per pool in these fixtures; moment dummies are pool-mapped but NOT `is_pool_proxy`. `pool_mu0_indices`/`pool_mu1_indices` are per-pool index arrays filled in `initialize_model` (:491/:498). `self.y` comes from the pydas `DASx` base and is in-place writable (`rs.y[:] = ...`).
- Phase rules permit the tripwire fixture: a polymer-phase ordinary species (mask False, unmapped) as product of a poly-event reaction passes `prods_phase_ok` (polymer.pyx:953-989); proxy concentrations are 1.0 and the rate is site-scaled by mu1/V_poly for non-end-group pool reactions (:1024-1040); `PolymerFluxArchetype.SAME_POOL == 1` (polymer.py:1429, `FLUX_SAME_POOL` polymer.pyx:78) is apportioned/non-UNRESOLVED and applies no moment flux.
- `HybridPolymerSystem` is a plain Python class in the .pyx (polymer.pyx:253; no .pxd entry) — new methods/attrs need no declaration.
- Iteration counter: `self.reaction_model.iteration_num`, incremented at `rmgpy/rmg/main.py:912`; initialized in `CoreEdgeReactionModel.__init__` (`rmgpy/rmg/model.py:229`).
- Post-simulate stash point: immediately before `self.rmg_memories[index].add_t_conv_N(t, x, len(obj))` (`rmgpy/rmg/main.py:1022`, unique anchor). Blueprint/engine gotcha: the runnable engine is `system.solver` (idiom at main.py:2120: `engine = getattr(system, "solver", None) or system`).
- Spawn pass: `CoreEdgeReactionModel._apply_multipool_spawn_pass` (model.py:389-426), called from `enlarge` at model.py:954; its `process_polymer_candidates_multipool` call (:415-419) does NOT pass `iteration` or `flux_accumulator` (both kwargs already exist in the signature, polymer.py:2419-2428; the `iteration=` at model.py:424 belongs to `apply_spawn_intents`).
- `triggering_moles=float(getattr(cand, "amount", 1.0))` is at polymer.py:2526; its consumer is `drain_spawn_intents` seeding `mu_k = N·DP^k` at polymer.py:3123-3125. Running-log item #14 already records the defect.
- `MassFluxAccumulator.record(motif_key, mass, iteration)` evicts entries older than `window` relative to the recording iteration (polymer.py:2245-2255); `flux` sums the window (:2257).
- The doc note to remove ("NOT YET ACTIVE") is at `docs/multi_pool_design.md:481-483` — at the end of **§10 Known limitations**, not §8 as the spec says (reported; remove it where it lives). §4.4 to rewrite is at :149-161.
- EPDM deck: `~/runs/RMG/epdm_v0_2026-06-06b/` (contains `input.py`); baseline log lines `MODEL GENERATION COMPLETED` and `The final model core has 26 species and 28 reactions`; sidecar key path `d["conventions"]["configured_pools"] == ["epdm"]`.
- Group-isomorphism check verified live: permuted-atom-order C-C-O group adjlists are `is_isomorphic` → True; `discover_repeat_motif` on the phenolic trimer SMILES returns a `Group`.
- `parent_polymer` fixture label is `"PE"` (polymerMultiPoolTest.py:64-75); the `TestReactionModelIntegration` parent is `"PS"`.

## File-structure map

| File | Change | Responsibility |
|---|---|---|
| `rmgpy/solver/polymer.pyx` | modify | change (a): gross production/consumption writes for ordinary core species in residual sections 3/4 (cost-gated, spec §4.6); `PolymerPoolConfig.monomer_mw_g_mol` field; `spawn_gate_flux_snapshot()` 3-tuple engine method |
| `rmgpy/rmg/polymer_input.py` | modify | plumb monomer MW (g/mol) into `PolymerPoolConfig` in `PolymerPool.to_config()` (:779/:807 — the only PolymerPoolConfig constructor; MassTransfer.to_config at :845 untouched) |
| `rmgpy/polymer.py` | modify | `MassFluxAccumulator.gate_statistic`/`window_occupancy`; `MotifLedgerEntry` (pair representatives), `_ledger_lookup`, `_snapshot_event_mass`, `_spawn_gate_fraction` (deduped denominator), `_evaluate_spawn_gate`; Phase-D/E gate flip; DELETE `_estimate_relative_flux`; `triggering_moles` TODO |
| `rmgpy/rmg/model.py` | modify | ledger/accumulator/snapshot attrs on `CoreEdgeReactionModel`; pass `iteration` into the spawn pass |
| `rmgpy/rmg/main.py` | modify | stash the 3-tuple snapshot + iteration on the reaction model after each `simulate()` |
| `docs/multi_pool_design.md` | modify | §4.4 rewrite (amended active-gate semantics); remove the §10 "NOT YET ACTIVE" note |
| `test/rmgpy/solver/solverPolymerTest.py` | modify | gate fixture helpers; `TestIntegratedSpawnGateTripwire` (2 tests, RED-FIRST); `TestSpawnGateFluxSnapshot` (2 tests); ε-direction end-to-end test |
| `test/rmgpy/polymerMultiPoolTest.py` | modify | 2 accumulator tests; ledger tests (isomorphism + dedup); gate scaffolding + 6 gate tests; 3 re-baselined first-sighting tests; gamed test DELETED |
| `~/runs/RMG/epdm_gate_cost_2026-06-10/` | create (out-of-repo) | cost-gate timing runs: `before_1`, `before_2` (Task 2, pre-.pyx-change), `after_1`, `after_2` (Task 3) |
| `~/.claude/projects/-home-alon-Code-RMG-Py/memory/project_polymer_branch_review.md` | modify (final task) | mark running-log item #8 implemented |

**Test arithmetic:** 380 baseline − 1 deleted (`test_mass_flux_below_threshold_defers_spawn`) + 15 new (2 integrated tripwire + 2 snapshot + 1 ε-direction = 5 solver-side; 2 accumulator + 2 ledger + 6 gate behavior = 10 multipool-side) = **394 passed, 2 skipped** at the end. Per file: solverPolymerTest 41 → 46; polymerMultiPoolTest 26 → 35.

All pytest commands run from the repo root `/home/alon/Code/RMG-Py`.

---

## Task 1: Integrated tripwire — RED-FIRST checkpoint (spec §7.1; gates the whole plan)

**Files:**
- `test/rmgpy/solver/solverPolymerTest.py` (import + helpers after the `_KIN` block at ~line 81; test class at end of file)

This task touches NO production code. It writes both halves of the spec-§7.1 integrated tripwire, runs them against CURRENT HEAD, and **the numerator half MUST FAIL at the gross-array assertion** — the proof that ordinary species have no gross writes today (the born-dead hole, spec §3.1). That red run is an explicit checkpoint: if it does NOT fail there, STOP the plan and report (the amendment's premise would be wrong). The tests are then committed `xfail(strict=True)` so the suite stays green while preserving the red proof; Task 3 removes the markers as part of going green.

- [ ] **Add the import.** In `test/rmgpy/solver/solverPolymerTest.py`, Edit:

  old_string:
```python
import numpy as np
import pytest
```
  new_string:
```python
import dataclasses

import numpy as np
import pytest
```

- [ ] **Add the gate fixture helpers.** Edit (insert immediately after the `_KIN` block, before `class TestHybridPolymerReactor:`):

  old_string:
```python
_KIN = dict(kinetics=Arrhenius(A=(2.0, "1/s"), n=0.0, Ea=(0.0, "kcal/mol"), T0=(298.15, "K")),
            reversible=False)
```
  new_string:
```python
_KIN = dict(kinetics=Arrhenius(A=(2.0, "1/s"), n=0.0, Ea=(0.0, "kcal/mol"), T0=(298.15, "K")),
            reversible=False)


def _gate_pool_config(monomer_mw_g_mol=28.0):
    """PolymerPoolConfig for the spawn-gate fixtures.

    ``monomer_mw_g_mol`` is added by the mass-flux-spawn-gate change (spec
    2026-06-10 §4.1). The field-presence guard below is deliberate RED-FIRST
    scaffolding, NOT a compatibility shim: it lets the Task-1 integrated
    tripwire run on pre-change HEAD and die at the GROSS-ARRAY assertion
    (the born-dead mechanism, spec §3.1) instead of at fixture construction.
    Once the field lands, the guard always takes the kwargs branch.
    """
    kwargs = dict(label="A", xs=2, explicit_dp_to_species_index={},
                  mu_indices=(1, 2, 3), monomer_poly_index=None)
    if any(f.name == "monomer_mw_g_mol"
           for f in dataclasses.fields(PolymerPoolConfig)):
        kwargs["monomer_mw_g_mol"] = monomer_mw_g_mol
    return PolymerPoolConfig(**kwargs)


def _one_pool_gate_species(rep_smiles="CCO"):
    """Species + mask for the one-pool spawn-gate fixture: canonical proxy A
    at 0, A_mu0/1/2 at 1-3, ordinary POLYMER-PHASE species R at 4 (stands in
    for an absorbed representative — representative status is a python/ledger
    concept; to the solver it is ANY ordinary species produced by
    pool-touching chemistry), gas G at 5."""
    sp = {
        "A": _spc("CCCC", "A"),
        "A_mu0": _spc("CO", "A_mu0"), "A_mu1": _spc("C=O", "A_mu1"), "A_mu2": _spc("C#N", "A_mu2"),
        "R": _spc(rep_smiles, "R"),
        "G": _spc("[CH3]", "G"),
    }
    core = [sp["A"], sp["A_mu0"], sp["A_mu1"], sp["A_mu2"], sp["R"], sp["G"]]
    mask = np.array([False] * 5 + [True], dtype=bool)
    return sp, core, mask


def _one_pool_gate_rs(rxn, core, mask, moments, monomer_mw_g_mol=28.0):
    rs = HybridPolymerSystem(
        T=800.0, P=1.0e5, initial_mole_fractions={core[5]: 0.0}, V_poly=1.0,
        polymer_pools=[_gate_pool_config(monomer_mw_g_mol)], mass_transfer=[],
        gas_species_mask=mask.copy(), constant_gas_volume=False,
        initial_polymer_moments={"A": moments}, termination=[],
    )
    rs.initialize_model(core, [rxn], [], [])
    return rs
```

- [ ] **Write the tripwire (both halves, NO xfail markers yet).** Append this class at the very end of the file:

```python
class TestIntegratedSpawnGateTripwire:
    """Spec §7.1 — INTEGRATED live-path tripwire, two halves, RED-FIRST.

    A real HybridPolymerSystem integrated solve (the apportionment-trajectory
    forward-Euler idiom), NOT a fabricated snapshot. The SAME_POOL reaction
    A -> A + R is apportioned (non-UNRESOLVED) pool-touching chemistry whose
    product R is an ordinary (non-canonical-proxy) core species. The
    numerator half MUST fail on pre-change-(a) HEAD: gross arrays are
    maintained only in is_pool_proxy branches today (spec §3.1) — that red
    run is the proof the fix is real and the gate to executing the rest of
    the plan. The polymer.pyx line citations in spec §3 are informational;
    these tests assert the MECHANISM, not the address.
    """

    @staticmethod
    def _solve():
        sp, core, mask = _one_pool_gate_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["A"], sp["R"]], **_KIN)
        rxn.polymer_flux_archetype = 1  # SAME_POOL: apportioned, non-UNRESOLVED
        rs = _one_pool_gate_rs(rxn, core, mask, (1.0, 5.0, 30.0), monomer_mw_g_mol=28.0)
        # Short integrated solve (the trajectory-test idiom): evolve the
        # engine state, then evaluate the residual AT the evolved state so
        # the gross arrays hold the last evaluation — exactly what
        # spawn_gate_flux_snapshot() reads after a production simulate().
        y = rs.y.copy()
        dt = 1e-4
        for _ in range(50):
            dn_dt = rs.residual(0.0, y, np.zeros_like(y))[0]
            assert np.all(np.isfinite(dn_dt))
            y = y + dt * dn_dt
        rs.y[:] = y
        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
        return sp, rs, dn_dt

    def test_numerator_half_ordinary_species_gross_is_real(self):
        """Numerator half (the regression that would have caught the
        born-dead class): the ordinary product R has a NONZERO gross entry
        in the snapshot, equal to max(0, core_species_production_rates[R])
        recomputed independently from the engine arrays, and
        g_R = gross * E[n] * monomer_MW under its parent pool's pool_stats.
        """
        sp, rs, dn_dt = self._solve()
        i_r = 4
        prod_r = float(rs.core_species_production_rates[i_r])
        # THE red assertion: fails on pre-change-(a) HEAD because ordinary
        # core species get dn_dt writes only — no gross record exists.
        assert prod_r > 0.0, (
            "ordinary core species R has no gross production record: the "
            "gross arrays are proxy-only (spec §3.1 born-dead hole; "
            "change (a) pending)"
        )
        # Independent recompute: irreversible reaction, V_poly = V_rxn = 1
        # -> R's gross production must equal its net dn_dt exactly.
        assert prod_r == pytest.approx(float(dn_dt[i_r]), rel=1e-12)

        gross, pool_stats, proxy_total = rs.spawn_gate_flux_snapshot()
        assert gross["R"] == pytest.approx(max(0.0, prod_r), rel=1e-12)
        e_n, mw = pool_stats["A"]
        assert e_n == pytest.approx(float(rs.y[2]) / float(rs.y[1]), rel=1e-12)
        assert mw == pytest.approx(28.0)
        g_r = gross["R"] * e_n * mw
        assert g_r == pytest.approx(
            max(0.0, prod_r) * (float(rs.y[2]) / float(rs.y[1])) * 28.0, rel=1e-12)
        assert g_r > 0.0

    def test_denominator_half_proxy_net_rerouted_gross_nonzero(self):
        """Denominator half: the canonical proxy's net dn_dt contribution is
        ~= 0 (the apportionment reroutes proxy flux to pool moments) while
        its gross entry is nonzero — an assertion that is only true of the
        GROSS array and dies if the denominator path is ever rewired to net
        rates."""
        sp, rs, dn_dt = self._solve()
        assert dn_dt[0] == pytest.approx(0.0, abs=1e-14)
        gross, pool_stats, proxy_total = rs.spawn_gate_flux_snapshot()
        assert gross["A"] > 0.0
        e_n, mw = pool_stats["A"]
        assert proxy_total == pytest.approx(gross["A"] * e_n * mw, rel=1e-12)
        assert proxy_total > 0.0
```

- [ ] **CHECKPOINT — run against current HEAD and confirm the red proof:**

```
~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/solver/solverPolymerTest.py::TestIntegratedSpawnGateTripwire -v
```

  Expected: **2 failed.**
  - `test_numerator_half_ordinary_species_gross_is_real` MUST fail at the `assert prod_r > 0.0` line with the "gross arrays are proxy-only" message (prod_r is 0.0 on HEAD — ordinary species have no gross writes). **This specific failure is the gate to the rest of the plan.** If it fails anywhere else (fixture error, phase disqualification) — fix the fixture until the failure is the gross assertion. If it PASSES — STOP THE PLAN and report: the amendment's premise (spec §3.1) would be contradicted by the live code.
  - `test_denominator_half_proxy_net_rerouted_gross_nonzero` fails with `AttributeError: ... 'spawn_gate_flux_snapshot'` (the method does not exist until Task 3); its `dn_dt[0] ≈ 0` assertion passes before that.

- [ ] **Mark both tests xfail (strict) so the suite stays green while preserving the red proof.** Two Edits:

  old_string:
```python
    def test_numerator_half_ordinary_species_gross_is_real(self):
```
  new_string:
```python
    @pytest.mark.xfail(strict=True, reason="change (a) pending: gross arrays proxy-only")
    def test_numerator_half_ordinary_species_gross_is_real(self):
```

  old_string:
```python
    def test_denominator_half_proxy_net_rerouted_gross_nonzero(self):
```
  new_string:
```python
    @pytest.mark.xfail(strict=True, reason="spawn_gate_flux_snapshot pending: Task 3 implements it with change (a)")
    def test_denominator_half_proxy_net_rerouted_gross_nonzero(self):
```

  (`strict=True` means these tests ERROR the suite if they ever unexpectedly pass — the red proof cannot silently rot. Task 3 removes both markers.)

- [ ] **Run the whole file — verify green-with-xfails:**

```
~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/solver/solverPolymerTest.py -v
```

  Expected: **41 passed, 2 xfailed**.

- [ ] **Commit:**

```
git add test/rmgpy/solver/solverPolymerTest.py
git commit -m "test: pin integrated spawn-gate tripwire red on HEAD (spec s7.1)

Both halves of the live-path tripwire on a real HybridPolymerSystem solve:
numerator (ordinary product R must have a real gross record + g = gross*
E[n]*MW under pool_stats) confirmed FAILING at the gross-array assertion
on current HEAD - the born-dead proof (spec s3.1): gross arrays are
maintained only in is_pool_proxy branches. Denominator half pins proxy net
dn_dt ~= 0 while gross is nonzero. Committed xfail(strict=True); the
change-(a) task removes the markers as it goes green.

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

## Task 2: EPDM wall-clock BASELINE timing (BEFORE any .pyx change — cost gate, spec §4.6)

**Files:** none in-repo. Creates `~/runs/RMG/epdm_gate_cost_2026-06-10/before_1` and `before_2` (out-of-repo scratch).

Change (a) is the first change since the apportionment work that adds unconditional per-species writes back into the residual hot loop, so it is COST-GATED on measured data. The baseline MUST be timed on the current tree (Task 1 committed only an xfail'd test — zero runtime impact), BEFORE Task 3 touches the .pyx. Each run takes a few minutes; runs go in fresh per-run directories (RMG restarts from existing output otherwise, which would contaminate timing).

- [ ] **Confirm the tree carries no .pyx modifications** (the baseline must be the pre-change residual):

```
git status --porcelain -- rmgpy/solver/ && git log --oneline -1
```

  Expected: no output from the status line (clean solver tree); HEAD is the Task-1 commit.

- [ ] **Prepare the timing directories:**

```
mkdir -p ~/runs/RMG/epdm_gate_cost_2026-06-10/before_1 ~/runs/RMG/epdm_gate_cost_2026-06-10/before_2
cp ~/runs/RMG/epdm_v0_2026-06-06b/input.py ~/runs/RMG/epdm_gate_cost_2026-06-10/before_1/
cp ~/runs/RMG/epdm_v0_2026-06-06b/input.py ~/runs/RMG/epdm_gate_cost_2026-06-10/before_2/
```

- [ ] **Baseline run 1:**

```
cd ~/runs/RMG/epdm_gate_cost_2026-06-10/before_1 && /usr/bin/time -v ~/anaconda3/envs/rmg_env/bin/python ~/Code/RMG-Py/rmg.py input.py > run.log 2> time.log; grep -E "Elapsed \(wall clock\)|Maximum resident" time.log; grep -c "MODEL GENERATION COMPLETED" RMG.log
```

  Expected: the wall-clock line (record it), and `1` (the run completed — a timing of a crashed run is invalid).

- [ ] **Baseline run 2:**

```
cd ~/runs/RMG/epdm_gate_cost_2026-06-10/before_2 && /usr/bin/time -v ~/anaconda3/envs/rmg_env/bin/python ~/Code/RMG-Py/rmg.py input.py > run.log 2> time.log; grep -E "Elapsed \(wall clock\)|Maximum resident" time.log; grep -c "MODEL GENERATION COMPLETED" RMG.log
```

  Expected: as above.

- [ ] **Record the two baseline wall-clock times** (write them into this plan file under this checkbox, or in the worker's report — Task 3's cost decision needs them):

```
BEFORE_1 = <m:ss.ss from before_1/time.log>
BEFORE_2 = <m:ss.ss from before_2/time.log>
```

  Sanity: the two runs should agree within run-to-run noise (~±5%). If they differ wildly, the machine is loaded — rerun until two consistent baselines exist; the cost gate is meaningless against a noisy baseline.

No commit (nothing in-repo changed).

---

## Task 3: Change (a) + monomer-MW plumb + 3-tuple snapshot — un-xfail, `make`, AFTER timing, cost decision

**Files:**
- `test/rmgpy/solver/solverPolymerTest.py` (remove the two xfail markers; add `TestSpawnGateFluxSnapshot`)
- `rmgpy/solver/polymer.pyx` (sections 3/4 gross writes; `monomer_mw_g_mol` field; `spawn_gate_flux_snapshot()` before `_chain_bundle` at :805)
- `rmgpy/rmg/polymer_input.py` (`PolymerPool.to_config`, :779-816)

This task touches a `.pyx` file: **`make` MUST be run from the repo root after the edits and BEFORE running tests.** The AFTER-timing runs happen at the end of this task, with ONLY this task's residual change on top of the Task-2 baseline (the gate flip and main.py stash land later and do not touch the residual) — so the timing isolates change (a)'s cost exactly.

- [ ] **Remove both xfail markers** (the tripwire now becomes this task's red phase). Two Edits in `test/rmgpy/solver/solverPolymerTest.py`:

  old_string:
```python
    @pytest.mark.xfail(strict=True, reason="change (a) pending: gross arrays proxy-only")
    def test_numerator_half_ordinary_species_gross_is_real(self):
```
  new_string:
```python
    def test_numerator_half_ordinary_species_gross_is_real(self):
```

  old_string:
```python
    @pytest.mark.xfail(strict=True, reason="spawn_gate_flux_snapshot pending: Task 3 implements it with change (a)")
    def test_denominator_half_proxy_net_rerouted_gross_nonzero(self):
```
  new_string:
```python
    def test_denominator_half_proxy_net_rerouted_gross_nonzero(self):
```

- [ ] **Write the failing snapshot unit tests.** Append this class at the end of `test/rmgpy/solver/solverPolymerTest.py` (after `TestIntegratedSpawnGateTripwire`):

```python
class TestSpawnGateFluxSnapshot:
    """spawn_gate_flux_snapshot() unit pins (spec §4.1): 3-tuple shape,
    all-core gross coverage, engine-attributed proxy event-mass total, and
    the E[n]*MW calibration (spec §7.6 / decision 3)."""

    def test_snapshot_three_tuple_covers_all_core_species(self):
        sp, core, mask = _one_pool_gate_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["A"], sp["R"]], **_KIN)
        rxn.polymer_flux_archetype = 1
        rs = _one_pool_gate_rs(rxn, core, mask, (1.0, 5.0, 30.0), monomer_mw_g_mol=28.0)

        rs.residual(0.0, rs.y, np.zeros_like(rs.y))
        gross, pool_stats, proxy_total = rs.spawn_gate_flux_snapshot()

        # gross: EVERY core species by label, max(0, production) — moment
        # dummies and the untouched gas species carry explicit zeros.
        assert set(gross.keys()) == {"A", "A_mu0", "A_mu1", "A_mu2", "R", "G"}
        assert gross["A_mu0"] == 0.0 and gross["A_mu1"] == 0.0 and gross["A_mu2"] == 0.0
        assert gross["G"] == 0.0
        assert gross["R"] == pytest.approx(
            max(0.0, float(rs.core_species_production_rates[4])), rel=1e-12)
        assert gross["R"] > 0.0
        assert gross["A"] > 0.0
        # pool_stats: pool label -> (E[n], monomer MW), live E[n] = mu1/mu0.
        assert set(pool_stats.keys()) == {"A"}
        e_n, mw = pool_stats["A"]
        assert e_n == pytest.approx(5.0)
        assert mw == pytest.approx(28.0)
        # proxy_event_mass_total: engine-attributed CANONICAL PROXIES only
        # (species_to_pool_indices + is_pool_proxy); attributing the
        # ordinary R is the python ledger's job (spec §4.1).
        assert proxy_total == pytest.approx(gross["A"] * 5.0 * 28.0, rel=1e-12)

    def test_snapshot_e_n_calibration_dominates_fragment_mw(self):
        """Spec §7.6 (decision 3): one mole of representative production is
        one mole of EVENTS; the mass entering the motif class per event is a
        chain's worth (E[n]*monomer_MW), not the representative fragment's
        own MW. With parent-pool E[n]=60 and a ~3-monomer-sized
        representative (hexane vs C2H4 repeat unit) the calibrated
        event-mass must read ~20x a fragment-MW accounting."""
        sp, core, mask = _one_pool_gate_species(rep_smiles="CCCCCC")
        rxn = Reaction(reactants=[sp["A"]], products=[sp["A"], sp["R"]], **_KIN)
        rxn.polymer_flux_archetype = 1
        # E[n] = 60; realizable: mu0*mu2 = 3700 >= mu1^2 = 3600.
        rs = _one_pool_gate_rs(rxn, core, mask, (1.0, 60.0, 3700.0), monomer_mw_g_mol=28.0)

        rs.residual(0.0, rs.y, np.zeros_like(rs.y))
        gross, pool_stats, _ = rs.spawn_gate_flux_snapshot()

        prod = float(rs.core_species_production_rates[4])
        assert prod > 0.0
        e_n, mw = pool_stats["A"]
        g_r = gross["R"] * e_n * mw  # the gate's representative g_i (spec §3)
        assert g_r == pytest.approx(prod * 60.0 * 28.0, rel=1e-12)
        rep_mw_g_mol = sp["R"].molecule[0].get_molecular_weight() * 1000.0  # ~86.18
        ratio = g_r / (prod * rep_mw_g_mol)
        assert ratio == pytest.approx(60.0 * 28.0 / rep_mw_g_mol, rel=1e-12)
        assert 15.0 < ratio < 25.0, "E[n]-calibrated mass must read ~20x the fragment-MW accounting"
```

- [ ] **Run — verify the red phase:**

```
~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/solver/solverPolymerTest.py::TestIntegratedSpawnGateTripwire test/rmgpy/solver/solverPolymerTest.py::TestSpawnGateFluxSnapshot -v
```

  Expected: **4 failed** — the numerator half again at the gross assertion (the red proof, now un-xfailed), the other three with `AttributeError: ... 'spawn_gate_flux_snapshot'` / missing field behavior.

- [ ] **Change (a), reactant slots.** In `rmgpy/solver/polymer.pyx`, Edit:

  old_string:
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
  new_string:
```python
            if self.is_pool_proxy[r0]:
                proxy_activity[r0] += abs_flux
                if core_rxn:
                    self.core_species_consumption_rates[r0] += rf
                    self.core_species_production_rates[r0] += rr
            elif core_rxn:
                dn_dt[r0] -= r_mol_s
                # Change (a) (spec 2026-06-10 mass-flux-spawn-gate s3.1/s4.6):
                # gross bookkeeping for ORDINARY core species too (simple.pyx
                # parity), so absorbed-representative production has a real
                # record for spawn_gate_flux_snapshot(). Cost-gated on the
                # EPDM deck (<= ~5% wall-clock; fallback: ledger-tracked-only
                # writes via a flag array).
                self.core_species_consumption_rates[r0] += rf
                self.core_species_production_rates[r0] += rr

            if r1 != -1:
                if self.is_pool_proxy[r1]:
                    proxy_activity[r1] += abs_flux
                    if core_rxn:
                        self.core_species_consumption_rates[r1] += rf
                        self.core_species_production_rates[r1] += rr
                elif core_rxn:
                    dn_dt[r1] -= r_mol_s
                    self.core_species_consumption_rates[r1] += rf
                    self.core_species_production_rates[r1] += rr

            if r2 != -1 and core_rxn:
                # r2 is always ordinary today (mirror the proxy logic above
                # if polymer is ever allowed in r2); change (a) gross writes
                # apply as for any ordinary reactant.
                dn_dt[r2] -= r_mol_s
                self.core_species_consumption_rates[r2] += rf
                self.core_species_production_rates[r2] += rr
```

- [ ] **Change (a), product slots.** Edit:

  old_string:
```python
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
```
  new_string:
```python
                if p_idx_tmp < n_core:
                    if self.is_pool_proxy[p_idx_tmp]:
                        proxy_activity[p_idx_tmp] += abs_flux
                        if core_rxn:
                            self.core_species_production_rates[p_idx_tmp] += rf
                            self.core_species_consumption_rates[p_idx_tmp] += rr
                    elif core_rxn:
                        # feature polymer species (e.g., PS_rad, scission oligomer, etc.) stays explicit
                        dn_dt[p_idx_tmp] += r_mol_s
                        # Change (a): gross bookkeeping for ordinary products
                        # (simple.pyx parity; see the reactant-side comment).
                        self.core_species_production_rates[p_idx_tmp] += rf
                        self.core_species_consumption_rates[p_idx_tmp] += rr
                else:
```

- [ ] **Add the monomer-MW field.** Edit:

  old_string:
```python
    monomer_poly_index: Optional[int] = None

    # Kinetics Parameters
```
  new_string:
```python
    monomer_poly_index: Optional[int] = None

    # Monomer (repeat-unit) molecular weight [g/mol]. Consumed by
    # spawn_gate_flux_snapshot (E[n]*MW event-mass calibration, spec
    # 2026-06-10-mass-flux-spawn-gate-design.md §3). Default 0.0 zeroes the
    # pool's snapshot contribution (honest degradation: the spawn gate
    # defers; nothing is fabricated).
    monomer_mw_g_mol: float = 0.0

    # Kinetics Parameters
```

- [ ] **Add the 3-tuple engine method.** Edit (insert before `_chain_bundle`):

  old_string:
```python
    def _chain_bundle(self, int pool_idx, y, double V_poly, bint end_group):
```
  new_string:
```python
    def spawn_gate_flux_snapshot(self):
        """Engine half of the mass-flux spawn gate.

        Spec: docs/superpowers/specs/2026-06-10-mass-flux-spawn-gate-design.md
        §3/§4.1 (AMENDED). The engine cannot attribute representatives (it
        has no ledger), so the split is: the engine reads arrays, the
        python-side gate does ledger-dependent attribution. Returns a
        3-tuple ``(gross, pool_stats, proxy_event_mass_total)``:

        1. ``gross``: dict core-species label ->
           ``max(0, core_species_production_rates[i])`` for ALL core species
           (ordinary species have real entries since change (a), spec §4.6).
           Labels are ``Species.label`` verbatim (the same labels the ledger
           records at absorption); a duplicate label would overwrite — core
           labels are unique in practice (make_new_species suffixes).
        2. ``pool_stats``: dict pool label -> ``(E_n, monomer_mw_g_mol)``
           with ``E_n = y[mu1]/y[mu0] if y[mu0] > SMALL_EPS else 0.0`` (the
           tail_mean guard idiom; errs toward deferral — mu0 exhausted ->
           g_i = 0 -> the gate defers).
        3. ``proxy_event_mass_total``: float — sum of
           ``gross_j * E_n[pool(j)] * mw[pool(j)]`` over CANONICAL pool
           proxies (``species_to_pool_indices[j] != -1 and is_pool_proxy[j]``
           — the engine CAN attribute those).

        GROSS production, never net dn_dt: canonical proxies have
        dn_dt ~= 0 BY DESIGN (the archetype apportionment reroutes their
        flux to pool moments) and ordinary species net to ~0 at steady
        state. E[n] is read LIVE from the state vector (never
        recorded-and-stale). All terms are polymer-phase volumetric, so
        V_poly cancels in any fraction of these numbers. Pure read of
        already-maintained state — no bookkeeping here beyond change (a)'s
        residual writes.
        """
        gross = {}
        pool_stats = {}
        proxy_event_mass_total = 0.0
        stp = getattr(self, "species_to_pool_indices", None)
        prod = getattr(self, "core_species_production_rates", None)
        y = getattr(self, "y", None)
        if stp is None or prod is None or y is None:
            return gross, pool_stats, proxy_event_mass_total
        n_core = self.num_core_species
        # Live per-pool stats.
        n_pools = len(self.polymer_pools)
        e_n_by_pool = [0.0] * n_pools
        for p in range(n_pools):
            i0 = self.pool_mu0_indices[p]
            i1 = self.pool_mu1_indices[p]
            if i0 < 0 or i1 < 0:
                continue
            mu0 = y[i0]
            mu1 = y[i1]
            e_n_by_pool[p] = mu1 / mu0 if mu0 > SMALL_EPS else 0.0
        for p in range(n_pools):
            pool = self.polymer_pools[p]
            mw = float(getattr(pool, "monomer_mw_g_mol", 0.0) or 0.0)
            pool_stats[pool.label] = (e_n_by_pool[p], mw)
        # index -> label for CORE species (species_index covers core+edge).
        labels = {}
        for spc, idx in self.species_index.items():
            if idx < n_core:
                labels[idx] = spc.label
        for i in range(n_core):
            label = labels.get(i)
            if label is None:
                continue
            g = max(0.0, float(prod[i]))
            gross[label] = g
            p = stp[i]
            if p >= 0 and self.is_pool_proxy[i]:
                e_n, mw = pool_stats[self.polymer_pools[p].label]
                proxy_event_mass_total += g * e_n * mw
        return gross, pool_stats, proxy_event_mass_total

    def _chain_bundle(self, int pool_idx, y, double V_poly, bint end_group):
```

- [ ] **Plumb the MW from the input blueprint.** In `rmgpy/rmg/polymer_input.py` (`PolymerPool.to_config` at :779 — the ONLY constructor of `PolymerPoolConfig`; `MassTransfer.to_config` at :845 builds `MassTransferConfig` and needs nothing), Edit:

  old_string:
```python
        # 3. Resolve Monomer Index
        monomer_idx = spc_map.get(self.monomer)

        return PolymerPoolConfig(
            label=self.label,
            xs=self.xs,
            explicit_dp_to_species_index=explicit_indices,
            mu_indices=mu_idxs,
            monomer_poly_index=monomer_idx,
            k_scission=self.k_scission,
            k_unzip=self.k_unzip
        )
```
  new_string:
```python
        # 3. Resolve Monomer Index
        monomer_idx = spc_map.get(self.monomer)

        # 4. Monomer (repeat-unit) MW [g/mol] for the spawn-gate snapshot
        #    (spec 2026-06-10 §3, same idiom as Polymer.monomer_mw_g_mol).
        #    Best-effort: 0.0 (-> the gate defers) when the monomer Species
        #    carries no resolvable structure.
        monomer_mw_g_mol = 0.0
        mol_list = getattr(self.monomer, "molecule", None)
        if mol_list:
            try:
                monomer_mw_g_mol = mol_list[0].get_molecular_weight() * 1000.0
            except Exception:
                monomer_mw_g_mol = 0.0

        return PolymerPoolConfig(
            label=self.label,
            xs=self.xs,
            explicit_dp_to_species_index=explicit_indices,
            mu_indices=mu_idxs,
            monomer_poly_index=monomer_idx,
            monomer_mw_g_mol=monomer_mw_g_mol,
            k_scission=self.k_scission,
            k_unzip=self.k_unzip
        )
```

- [ ] **Compile** (MANDATORY before tests — .pyx changed):

```
cd /home/alon/Code/RMG-Py && make
```

  Expected: compiles cleanly (Cython regenerates `rmgpy/solver/polymer.c`, no errors).

- [ ] **Run the solver file — verify green:**

```
~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/solver/solverPolymerTest.py -v
```

  Expected: **45 passed** (41 baseline + 2 tripwire now real-green + 2 snapshot), 0 xfailed, 0 failures.

- [ ] **Regression sanity** (gross arrays are diagnostics consumed elsewhere — chemkin/cantera/artifact paths must not notice):

```
~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/polymerMultiPoolTest.py test/rmgpy/rmg/modelTest.py test/rmgpy/polymerTest.py test/rmgpy/polymerArtifactTest.py -q
```

  Expected: all pass at their Task-1 baselines (26 + modelTest + polymerTest + artifact, unchanged), 0 failures.

- [ ] **AFTER timing, run 1 and 2** (same machine, same session as Task 2):

```
mkdir -p ~/runs/RMG/epdm_gate_cost_2026-06-10/after_1 ~/runs/RMG/epdm_gate_cost_2026-06-10/after_2
cp ~/runs/RMG/epdm_v0_2026-06-06b/input.py ~/runs/RMG/epdm_gate_cost_2026-06-10/after_1/
cp ~/runs/RMG/epdm_v0_2026-06-06b/input.py ~/runs/RMG/epdm_gate_cost_2026-06-10/after_2/
cd ~/runs/RMG/epdm_gate_cost_2026-06-10/after_1 && /usr/bin/time -v ~/anaconda3/envs/rmg_env/bin/python ~/Code/RMG-Py/rmg.py input.py > run.log 2> time.log; grep -E "Elapsed \(wall clock\)" time.log; grep -c "MODEL GENERATION COMPLETED" RMG.log
cd ~/runs/RMG/epdm_gate_cost_2026-06-10/after_2 && /usr/bin/time -v ~/anaconda3/envs/rmg_env/bin/python ~/Code/RMG-Py/rmg.py input.py > run.log 2> time.log; grep -E "Elapsed \(wall clock\)" time.log; grep -c "MODEL GENERATION COMPLETED" RMG.log
```

  Expected: two wall-clock lines (record as AFTER_1/AFTER_2) and `1` from each completion grep.

- [ ] **COST DECISION (spec §4.6 — decided on measured data, not assumption):**
  - Compute `mean(AFTER_1, AFTER_2) / mean(BEFORE_1, BEFORE_2)`.
  - **Accept** if ≤ 1.05 (slowdown within ~5% / run-to-run noise): proceed to the commit step, and record all four numbers + the ratio in the commit message.
  - **Fallback decision point if it bites (> 1.05):** STOP — do NOT commit the unconditional form. The documented fallback is gross writes only for LEDGER-TRACKED species: a per-species `int8` flag array on the engine (proxies always set; OR'd with the labels present in `reaction_model.polymer_motif_ledger` representatives, plumbed through `initialize_model`), gating the new section-3/4 writes. Implement it, `make`, re-run the solver tests (the tripwire's R must then be made ledger-tracked in the fixture via the plumb), re-time, and report the numbers to Alon before continuing the plan. The spec prefers the unconditional form for simplicity and simple.pyx parity — switch only on measured evidence.

- [ ] **Commit** (fill in the measured numbers):

```
git add rmgpy/solver/polymer.pyx rmgpy/rmg/polymer_input.py test/rmgpy/solver/solverPolymerTest.py
git commit -m "solver: gross flux records for all core species + spawn-gate snapshot

Change (a) (spec 2026-06-10 s3.1/s4.6): residual sections 3/4 now maintain
gross production/consumption for ordinary core species (simple.pyx parity,
core reactions only) - absorbed-representative production has a real
record. spawn_gate_flux_snapshot() returns the amended 3-tuple (gross for
ALL core species, pool_stats label->(E[n],MW) with the SMALL_EPS guard,
engine-attributed proxy_event_mass_total). PolymerPoolConfig.monomer_mw_g_mol
plumbed from PolymerPool.to_config. Integrated tripwire un-xfailed: red on
pre-change HEAD at the gross assertion, green now.

Cost gate (EPDM deck, spec s4.6): before <BEFORE_1>/<BEFORE_2>, after
<AFTER_1>/<AFTER_2>, ratio <RATIO> - within the <=~5% acceptance.

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

## Task 4: `MassFluxAccumulator.gate_statistic` / `window_occupancy`

**Files:**
- `rmgpy/polymer.py` (append methods after `flux`, line ~2262)
- `test/rmgpy/polymerMultiPoolTest.py` (extend `TestMassFluxAccumulator`, ~line 149)

- [ ] **Write the failing tests.** In `test/rmgpy/polymerMultiPoolTest.py`, Edit (append two methods to `TestMassFluxAccumulator`):

  old_string:
```python
    def test_unknown_motif_returns_zero(self):
        from rmgpy.polymer import MassFluxAccumulator

        acc = MassFluxAccumulator(window=3)
        assert acc.flux("never_recorded") == 0.0
```
  new_string:
```python
    def test_unknown_motif_returns_zero(self):
        from rmgpy.polymer import MassFluxAccumulator

        acc = MassFluxAccumulator(window=3)
        assert acc.flux("never_recorded") == 0.0

    def test_gate_statistic_divides_by_fixed_window_length(self):
        """Spec §7.11: sum/N with zero-filled semantics — one record divides
        by the FIXED window N, not the record count, so a single-snapshot
        spike must be N x the bar to clear the gate."""
        from rmgpy.polymer import MassFluxAccumulator

        acc = MassFluxAccumulator(window=3)
        acc.record("m", 0.06, iteration=4)
        assert acc.gate_statistic("m") == pytest.approx(0.02)
        acc.record("m", 0.03, iteration=5)
        assert acc.gate_statistic("m") == pytest.approx(0.03)
        assert acc.gate_statistic("never_recorded") == 0.0

    def test_gate_statistic_evicts_at_window_edge(self):
        """Spec §7.11: records older than window iterations relative to the
        latest recording iteration are evicted from the statistic."""
        from rmgpy.polymer import MassFluxAccumulator

        acc = MassFluxAccumulator(window=3)
        acc.record("m", 0.09, iteration=1)
        acc.record("m", 0.03, iteration=2)
        acc.record("m", 0.03, iteration=4)   # cutoff 4-3+1=2: evicts iteration 1
        assert acc.window_occupancy("m") == 2
        assert acc.gate_statistic("m") == pytest.approx(0.02)
```

- [ ] **Run — verify failure:**

```
~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/polymerMultiPoolTest.py::TestMassFluxAccumulator -v
```

  Expected: 3 pass (existing), 2 fail with `AttributeError: 'MassFluxAccumulator' object has no attribute 'gate_statistic'`.

- [ ] **Implement.** In `rmgpy/polymer.py`, Edit:

  old_string:
```python
    def flux(self, motif_key: str) -> float:
        """Sum of masses currently in the rolling window for ``motif_key``."""
        if motif_key not in self._records:
            return 0.0
        return sum(m for (_, m) in self._records[motif_key])
```
  new_string:
```python
    def flux(self, motif_key: str) -> float:
        """Sum of masses currently in the rolling window for ``motif_key``."""
        if motif_key not in self._records:
            return 0.0
        return sum(m for (_, m) in self._records[motif_key])

    def gate_statistic(self, motif_key: str) -> float:
        """Window sum divided by the FIXED window length (zero-filled
        semantics, spec 2026-06-10 §4.4 step 4): a single-snapshot spike must
        be ``window``x the bar to clear the gate; a channel persisting at
        fraction f for ``window`` iterations reads f."""
        return self.flux(motif_key) / float(self.window)

    def window_occupancy(self, motif_key: str) -> int:
        """Number of records currently in the rolling window."""
        return len(self._records.get(motif_key, []))
```

- [ ] **Run — verify pass:**

```
~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/polymerMultiPoolTest.py -v
```

  Expected: **28 passed** (26 baseline + 2 new).

- [ ] **Commit:**

```
git add rmgpy/polymer.py test/rmgpy/polymerMultiPoolTest.py
git commit -m "polymer: add windowed gate statistic to MassFluxAccumulator

sum/N zero-filled semantics (spec 2026-06-10 s4.4 step 4) + window
occupancy for the gate's deferral log line.

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

## Task 5: Motif ledger — pair representatives, deduped denominator, gate evaluator

**Files:**
- `rmgpy/polymer.py` (insert machinery before `def _bfs_grow_heavy_subset(` at line ~2264; the `_estimate_relative_flux` stub stays untouched until Task 6)
- `test/rmgpy/polymerMultiPoolTest.py` (new test class)

- [ ] **Write the failing tests.** In `test/rmgpy/polymerMultiPoolTest.py`, Edit (insert a new class before `class TestProcessPolymerCandidatesMultiPool:` at ~line 183):

  old_string:
```python
class TestProcessPolymerCandidatesMultiPool:
    """Multi-pool aware product classification + spawn-intent generation."""
```
  new_string:
```python
class TestMotifLedger:
    """Group-isomorphism motif ledger + amended fraction math
    (spec 2026-06-10 §3/§4.3)."""

    def test_ledger_lookup_is_group_isomorphism_not_string_key(self):
        """Spec §7.9: the same motif Group with permuted atom ordering must
        hit ONE ledger entry — pins the isomorphism-not-string-key choice."""
        from rmgpy.molecule.group import Group
        from rmgpy.polymer import MotifLedgerEntry, _ledger_lookup

        g1 = Group().from_adjacency_list("""
1 C u0 {2,S}
2 C u0 {1,S} {3,S}
3 O u0 {2,S}
""")
        # Same motif, permuted atom ordering.
        g2 = Group().from_adjacency_list("""
1 O u0 {2,S}
2 C u0 {1,S} {3,S}
3 C u0 {2,S}
""")
        entry = MotifLedgerEntry(motif=g1, accumulator_key="motif-0")
        ledger = [entry]
        assert _ledger_lookup(ledger, g2) is entry, (
            "permuted atom ordering must resolve to the same ledger entry "
            "(Group isomorphism, not a canonical string key)"
        )
        assert _ledger_lookup(ledger, Group().from_adjacency_list("""
1 C u0 {2,S}
2 C u0 {1,S}
""")) is None, "a structurally different motif must miss"

    def test_denominator_dedups_shared_representative_numerators_do_not(self):
        """Spec §3 (amended): a species appearing in multiple motif entries
        counts ONCE in the denominator (deduped), while EVERY motif-entry
        numerator that lists it sees its mass — the stated multi-motif
        double-counting decision (the two motifs compete for DIFFERENT pool
        slots). Each fraction stays in [0,1] (numerators are subsets of
        denominator terms); the SUM across motifs may exceed 1 — accepted."""
        from rmgpy.molecule.group import Group
        from rmgpy.polymer import MotifLedgerEntry, _spawn_gate_fraction

        g1 = Group().from_adjacency_list("""
1 C u0 {2,S}
2 C u0 {1,S} {3,S}
3 O u0 {2,S}
""")
        g2 = Group().from_adjacency_list("""
1 C u0 {2,S}
2 C u0 {1,S}
""")
        e1 = MotifLedgerEntry(motif=g1, accumulator_key="motif-0",
                              representatives=[("X", "PE")])
        e2 = MotifLedgerEntry(motif=g2, accumulator_key="motif-1",
                              representatives=[("X", "PE"), ("Y", "PE")])
        ledger = [e1, e2]
        # pool_stats PE: E[n]=2, MW=5 -> g_X = 0.3*10 = 3.0, g_Y = 0.1*10 = 1.0;
        # engine-attributed canonical-proxy total = 6.0.
        snapshot = ({"X": 0.3, "Y": 0.1}, {"PE": (2.0, 5.0)}, 6.0)

        # Denominator = proxies (6) + DEDUPED reps (g_X + g_Y = 4) = 10 —
        # X appears in two entries but counts once.
        assert _spawn_gate_fraction(e1, ledger, snapshot) == pytest.approx(3.0 / 10.0)
        assert _spawn_gate_fraction(e2, ledger, snapshot) == pytest.approx(4.0 / 10.0)

        # A representative whose recorded parent pool has no stats defers
        # (g = 0, the deferral direction); a missing snapshot defers.
        e3 = MotifLedgerEntry(motif=g1, accumulator_key="motif-2",
                              representatives=[("Z", "NOPOOL")])
        assert _spawn_gate_fraction(e3, ledger + [e3], snapshot) == pytest.approx(0.0)
        assert _spawn_gate_fraction(e1, ledger, None) == 0.0


class TestProcessPolymerCandidatesMultiPool:
    """Multi-pool aware product classification + spawn-intent generation."""
```

- [ ] **Run — verify failure:**

```
~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/polymerMultiPoolTest.py::TestMotifLedger -v
```

  Expected: 2 failures — `ImportError: cannot import name 'MotifLedgerEntry'`.

- [ ] **Implement the machinery.** In `rmgpy/polymer.py`, Edit (insert before `_bfs_grow_heavy_subset`; this also pre-builds the gate evaluator that Task 6 wires into Phase D):

  old_string:
```python
def _bfs_grow_heavy_subset(
```
  new_string:
```python
@dataclass
class MotifLedgerEntry:
    """Per-motif spawn-gate ledger entry (design doc §4.4; spec 2026-06-10
    §4.3, AMENDED).

    Lives in ``reaction_model.polymer_motif_ledger`` — in-memory ONLY (an RMG
    restart resets windows and deferred motifs re-earn their bar:
    correct-but-loud, same philosophy as unstamped-reaction demotion).
    Lookup is by Group isomorphism (:func:`_ledger_lookup`), never a canonical
    string key. ``representatives`` are ``(species_label, parent_pool_label)``
    pairs — the parent pool recorded at absorption (for Phase-D deferred
    candidates: the pool that would parent the spawn intent, currently
    ``pool_registry[0]``; if multi-parent attribution ever lands, this
    follows it). E[n]/monomer_MW for a representative are read LIVE from that
    pool's stats in the snapshot — live freshness kept, the wrong-mapping
    hole closed. ``accumulator_key`` is the opaque per-entry id used with
    :class:`MassFluxAccumulator`.
    """

    motif: Group
    accumulator_key: str
    representatives: List[Tuple[str, str]] = field(default_factory=list)
    last_recorded_iteration: int = -1
    spawned: bool = False


def _ledger_lookup(
    ledger: List[MotifLedgerEntry],
    motif: Group,
) -> Optional[MotifLedgerEntry]:
    """Find the entry whose motif Group is isomorphic to ``motif``.

    Same matching idiom as :func:`similarity_merge` — sidesteps Group
    canonicalization. The ledger is ``max_pools``-scale, so the O(n)
    isomorphism scan is fine.
    """
    for entry in ledger:
        try:
            if motif.is_isomorphic(entry.motif):
                return entry
        except (NotImplementedError, AttributeError, ValueError):
            continue
    return None


def _snapshot_event_mass(
    snapshot: Tuple[Dict[str, float], Dict[str, Tuple[float, float]], float],
    species_label: str,
    parent_pool_label: str,
) -> float:
    """g_i for one representative per spec §3 (amended):
    ``gross[label] * E[n]_parent * monomer_MW_parent``, with the E[n]/MW pair
    read from the snapshot's ``pool_stats`` for the parent pool RECORDED AT
    ABSORPTION. Labels absent from ``gross`` (absorbed this iteration, not
    yet simulated) and parent pools absent from ``pool_stats`` contribute 0
    — stated, not incidental; both err toward deferral."""
    gross, pool_stats, _ = snapshot
    e_n, mw = pool_stats.get(parent_pool_label, (0.0, 0.0))
    return gross.get(species_label, 0.0) * e_n * mw


def _spawn_gate_fraction(
    entry: MotifLedgerEntry,
    ledger: List[MotifLedgerEntry],
    snapshot: Optional[Tuple[Dict[str, float], Dict[str, Tuple[float, float]], float]],
) -> float:
    """fraction(motif) per spec §3 (AMENDED), from a stashed engine snapshot.

    ``snapshot`` is the 3-tuple from
    ``HybridPolymerSystem.spawn_gate_flux_snapshot()``:
    ``(gross, pool_stats, proxy_event_mass_total)``.

        numerator   = sum of g_i over THIS entry's (label, parent_pool) pairs
        denominator = proxy_event_mass_total
                      + sum of g_i over DEDUPED representatives across the
                        WHOLE ledger (a species in multiple motif entries
                        counts ONCE here, keyed by species label)

    Numerators are subsets of denominator terms, so each motif's fraction is
    in [0,1]; the SUM across motifs may exceed 1 (the stated multi-motif
    double-counting decision — competing for different pool slots). No
    snapshot (iteration 0, or no polymer reaction system) or an empty
    denominator -> 0.0: honest degradation, the gate defers; no production
    code path fakes a number.
    """
    if not snapshot:
        return 0.0
    try:
        _, _, proxy_event_mass_total = snapshot
    except (TypeError, ValueError):
        return 0.0
    numerator = sum(
        _snapshot_event_mass(snapshot, lbl, pool_lbl)
        for (lbl, pool_lbl) in entry.representatives
    )
    deduped: Dict[str, str] = {}
    for e in ledger:
        for (lbl, pool_lbl) in e.representatives:
            deduped.setdefault(lbl, pool_lbl)
    representative_total = sum(
        _snapshot_event_mass(snapshot, lbl, pool_lbl)
        for lbl, pool_lbl in deduped.items()
    )
    denominator = float(proxy_event_mass_total) + representative_total
    if denominator <= 0.0:
        return 0.0
    return numerator / denominator


def _evaluate_spawn_gate(
    cand: 'Species',
    motif: Group,
    reaction_model: Any,
    iteration: int,
    flux_accumulator: Optional[MassFluxAccumulator],
    mass_flux_threshold: float,
) -> Tuple[bool, float, Optional[MotifLedgerEntry]]:
    """Mass-flux spawn gate (design doc §4.4; spec 2026-06-10 §4.4).

    Arrival-driven: runs only when a NEW candidate carrying ``motif`` reaches
    Phase D (deferred candidates are never re-presented). Records at most ONE
    snapshot-attributed gross mass-flux fraction per motif per RMG iteration
    — same-iteration burst arrivals re-check the existing window only,
    otherwise "windowed over N iterations" silently becomes "windowed over N
    arrivals". The gate statistic is the window sum divided by the FIXED
    window length (zero-filled): a single-snapshot spike must be window x the
    bar to clear it.

    Returns ``(spawn, statistic, entry)``. ``entry`` is ``None`` when there
    is no reaction model / accumulator to hold gate state (bare unit-test
    path) — then the candidate defers with statistic 0.0.
    """
    cand_label = getattr(cand, "label", "") or repr(cand)
    if reaction_model is None or flux_accumulator is None:
        logging.info(
            "Polymer spawn gate: no reaction-model gate state available; "
            "deferring spawn for %s (statistic 0.0 < bar %.4g).",
            cand_label, mass_flux_threshold,
        )
        return False, 0.0, None

    ledger = getattr(reaction_model, "polymer_motif_ledger", None)
    if ledger is None:
        ledger = []
        reaction_model.polymer_motif_ledger = ledger

    entry = _ledger_lookup(ledger, motif)
    if entry is None:
        entry = MotifLedgerEntry(motif=motif, accumulator_key=f"motif-{len(ledger)}")
        ledger.append(entry)

    if entry.spawned:
        # Should be unreachable: Phase A classifies arrivals against the
        # spawned pool first. Assert-log; never re-run the gate.
        logging.warning(
            "Polymer spawn gate: arrival %s hit already-spawned ledger entry "
            "%s — Phase A should have classified it against the spawned pool.",
            cand_label, entry.accumulator_key,
        )
        return False, 0.0, entry

    snapshot = getattr(reaction_model, "polymer_flux_snapshot", None)
    if entry.last_recorded_iteration < iteration:
        fraction = _spawn_gate_fraction(entry, ledger, snapshot)
        flux_accumulator.record(entry.accumulator_key, fraction, iteration)
        entry.last_recorded_iteration = iteration

    statistic = flux_accumulator.gate_statistic(entry.accumulator_key)
    spawn = statistic >= mass_flux_threshold
    if not spawn:
        logging.info(
            "Polymer spawn gate: deferring spawn for %s — statistic %.4g < "
            "bar %.4g (window %d/%d records%s).",
            cand_label, statistic, mass_flux_threshold,
            flux_accumulator.window_occupancy(entry.accumulator_key),
            flux_accumulator.window,
            "" if snapshot is not None else "; no flux snapshot stashed",
        )
    return spawn, statistic, entry


def _bfs_grow_heavy_subset(
```

- [ ] **Run — verify pass (and no regressions):**

```
~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/polymerMultiPoolTest.py -v
```

  Expected: **30 passed** (28 + 2 new).

- [ ] **Commit:**

```
git add rmgpy/polymer.py test/rmgpy/polymerMultiPoolTest.py
git commit -m "polymer: motif ledger with pair representatives + deduped denominator

MotifLedgerEntry holds (species_label, parent_pool_label) pairs recorded at
absorption (spec 2026-06-10 s4.3 AMENDED); _snapshot_event_mass reads
E[n]/MW live from the snapshot's pool_stats for the recorded parent pool;
_spawn_gate_fraction implements the amended s3 formula - numerator per
motif, denominator = proxy_event_mass_total + DEDUPED representative
event-mass across the whole ledger (fraction in [0,1]; multi-motif
double-counting across numerators is the stated decision).
_evaluate_spawn_gate: one record per motif per iteration, spawned-flag
assert-log skip. Not yet wired into Phase D.

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

## Task 6: Flip Phase D/E to the real gate; delete the stub and the gamed test; re-baseline the first-sighting tests; reaction-model wiring

**Files:**
- `rmgpy/polymer.py` (`_estimate_relative_flux` deleted, lines ~2399-2416; Phase D/E in `process_polymer_candidates_multipool` ~2502-2530; accumulator fallback at ~2450)
- `rmgpy/rmg/model.py` (import line 47; `__init__` after `self.iteration_num = 0` line 229; `_apply_multipool_spawn_pass` ~415)
- `test/rmgpy/polymerMultiPoolTest.py` (scaffolding after the `parent_polymer` fixture; `TestSpawnGateBehavior` new class; `test_mass_flux_below_threshold_defers_spawn` DELETED ~251; `test_novel_motif_spawns_one_pool` ~228, `test_iteration_boundary_pass_spawns_and_registers` ~365, `test_novel_product_spawns_and_registers_and_serializes` ~420 re-baselined)

This is the atomic flip: the gate goes live, so the three first-sighting tests and the gamed test MUST change in the same commit to keep the suite green. Fabricated snapshots below are the amended 3-tuple `(gross, pool_stats, proxy_event_mass_total)` — per spec §4.5 fabrication is a TEST device only; production code never fakes a number.

- [ ] **Add the test scaffolding.** In `test/rmgpy/polymerMultiPoolTest.py`, Edit (insert after the `parent_polymer` fixture, before the `# discover_repeat_motif` banner):

  old_string:
```python
# ---------------------------------------------------------------------------
# discover_repeat_motif
# ---------------------------------------------------------------------------
```
  new_string:
```python
# ---------------------------------------------------------------------------
# Spawn-gate test scaffolding (spec 2026-06-10-mass-flux-spawn-gate-design.md)
# ---------------------------------------------------------------------------

class _GateModel:
    """Bare reaction-model stand-in carrying the spawn-gate state the gate
    reads off the real CoreEdgeReactionModel (spec §4.3): the motif ledger,
    the stashed 3-tuple flux snapshot, and the shared accumulator."""

    def __init__(self, window=3):
        from rmgpy.polymer import MassFluxAccumulator

        self.polymer_motif_ledger = []
        self.polymer_flux_accumulator = MassFluxAccumulator(window=window)
        self.polymer_flux_snapshot = None


# pool_stats with E[n]=1, MW=1 make a representative's g equal its gross
# entry — fabricated-snapshot arithmetic reads off directly.
_PE_STATS = {"PE": (1.0, 1.0)}


def _phenolic(label):
    """A fresh phenolic-trimer arrival (novel motif w.r.t. the PE parent)."""
    from rmgpy.species import Species

    s = Species(smiles="Oc1ccc(Cc2ccc(Cc3ccc(O)cc3)cc2)cc1")
    s.label = label
    return s


def _gate_pass(model, parent, cand, iteration, threshold=0.01):
    """One Phase A-E pass for a single candidate against one parent pool."""
    from rmgpy.polymer import process_polymer_candidates_multipool

    return process_polymer_candidates_multipool(
        candidates=[cand],
        reaction_model=model,
        pool_registry=[parent],
        mass_flux_threshold=threshold,
        iteration=iteration,
    )


# ---------------------------------------------------------------------------
# discover_repeat_motif
# ---------------------------------------------------------------------------
```

- [ ] **Add the gate-behavior tests.** In `test/rmgpy/polymerMultiPoolTest.py`, Edit (insert the new class before `class TestMotifLedger:` added in Task 5):

  old_string:
```python
class TestMotifLedger:
    """Group-isomorphism motif ledger + amended fraction math
    (spec 2026-06-10 §3/§4.3)."""
```
  new_string:
```python
class TestSpawnGateBehavior:
    """The live Phase-D mass-flux gate (spec 2026-06-10 §4.4/§7.2-7.5/7.10).

    Snapshots are fabricated 3-tuples per spec §4.5: tests exercising Phase E
    supply a fabricated snapshot and a pre-populated ledger; production code
    never fakes a number. With _PE_STATS (E[n]=1, MW=1) a representative's g
    equals its gross entry, so fraction = gross / (proxy_total + rep_total).
    """

    def test_spike_spawns_at_second_arrival_iteration(self, parent_polymer):
        """Spec §7.2 branch 1 (floor arithmetic, decision 4): first sighting
        can never spawn; a single record >= window x bar clears the gate at
        the SECOND arrival's iteration."""
        model = _GateModel(window=3)

        # First arrival (iteration 1): zero representatives -> record 0, defer.
        model.polymer_flux_snapshot = ({}, dict(_PE_STATS), 1.0)
        _, intents = _gate_pass(model, parent_polymer, _phenolic("phen_1"), iteration=1)
        assert intents == [], "First sighting can never spawn"
        entry = model.polymer_motif_ledger[0]
        assert entry.representatives == [("phen_1", "PE")]
        assert entry.spawned is False

        # Second arrival (iteration 2): the first arrival now carries flux.
        # fraction = 0.5/(0.5+0.5) = 0.5 >= 3 x the 0.01 bar -> statistic
        # 0.5/3 clears the gate.
        model.polymer_flux_snapshot = ({"phen_1": 0.5}, dict(_PE_STATS), 0.5)
        _, intents = _gate_pass(model, parent_polymer, _phenolic("phen_2"), iteration=2)
        assert len(intents) == 1, "a >= window x bar motif spawns at its 2nd arrival's iteration"
        assert intents[0].mass_flux_at_spawn == pytest.approx(0.5 / 3.0)
        assert entry.spawned is True
        assert entry.representatives == [("phen_1", "PE"), ("phen_2", "PE")]

    def test_exactly_at_bar_defers_until_nth_recording_iteration(self, parent_polymer):
        """Spec §7.2 branch 2: an exactly-at-bar motif needs real records
        filling the window — spawn at the window-th REAL recording iteration
        (k+3 for window=3 with arrivals every iteration), NOT at second
        sighting."""
        model = _GateModel(window=3)
        bar = 0.01

        # Arrival 1 (iteration 1): no representatives -> record 0.0, defer.
        model.polymer_flux_snapshot = ({}, dict(_PE_STATS), 1.0)
        _, intents = _gate_pass(model, parent_polymer, _phenolic("phen_1"), iteration=1, threshold=bar)
        assert intents == []

        # The motif persists at EXACTLY the bar from here on:
        # fraction = 0.01/(0.99+0.01) = 0.01.
        model.polymer_flux_snapshot = ({"phen_1": 0.01}, dict(_PE_STATS), 0.99)
        _, intents = _gate_pass(model, parent_polymer, _phenolic("phen_2"), iteration=2, threshold=bar)
        assert intents == [], "window sum 0.01 / 3 is below the bar"
        _, intents = _gate_pass(model, parent_polymer, _phenolic("phen_3"), iteration=3, threshold=bar)
        assert intents == [], "window sum 0.02 / 3 is below the bar"

        # Arrival 4 (iteration 4): window holds iterations 2,3,4 (the
        # iteration-1 zero evicted) -> sum 0.03 / 3 == bar -> spawn.
        _, intents = _gate_pass(model, parent_polymer, _phenolic("phen_4"), iteration=4, threshold=bar)
        assert len(intents) == 1, "exactly-at-bar spawns at its window-th real recording iteration"
        assert intents[0].mass_flux_at_spawn == pytest.approx(bar)

    def test_trace_motif_stays_deferred_across_many_arrivals(self, parent_polymer):
        """Spec §7.3: a persistent below-bar motif (0.1% of gross event-mass
        vs the 1% bar) never spawns, however many arrivals it produces."""
        model = _GateModel(window=3)
        model.polymer_flux_snapshot = ({}, dict(_PE_STATS), 1.0)
        _, intents = _gate_pass(model, parent_polymer, _phenolic("phen_0"), iteration=0)
        assert intents == []
        for it in range(1, 8):
            # fraction = 0.001/(0.999+0.001) = 0.001 at every arrival.
            model.polymer_flux_snapshot = ({"phen_0": 0.001}, dict(_PE_STATS), 0.999)
            _, intents = _gate_pass(model, parent_polymer, _phenolic(f"phen_{it}"), iteration=it)
            assert intents == [], f"trace motif must stay deferred (iteration {it})"
        entry = model.polymer_motif_ledger[0]
        assert entry.spawned is False
        assert len(entry.representatives) == 8
        assert all(pool == "PE" for (_, pool) in entry.representatives)

    def test_same_iteration_burst_records_once_and_defers(self, parent_polymer):
        """Spec §7.5: three same-iteration arrivals at an ABOVE-bar fraction
        (0.02 > 0.01) still defer — the window holds ONE record, sum/N =
        0.02/3 < bar — and the ledger must hold all three representatives."""
        from rmgpy.polymer import MotifLedgerEntry, discover_repeat_motif

        model = _GateModel(window=3)
        seed = _phenolic("phen_seed")
        motif = discover_repeat_motif(seed.molecule[0])
        assert motif is not None
        model.polymer_motif_ledger.append(MotifLedgerEntry(
            motif=motif, accumulator_key="motif-0",
            representatives=[("phen_seed", "PE")],
        ))
        # fraction = 0.02/(0.98+0.02) = 0.02.
        model.polymer_flux_snapshot = ({"phen_seed": 0.02}, dict(_PE_STATS), 0.98)

        for n in range(3):
            _, intents = _gate_pass(model, parent_polymer, _phenolic(f"phen_burst_{n}"), iteration=5)
            assert intents == [], "above-bar single record must still defer (sum/N)"

        assert len(model.polymer_motif_ledger) == 1
        entry = model.polymer_motif_ledger[0]
        assert model.polymer_flux_accumulator.window_occupancy("motif-0") == 1, (
            "three same-iteration arrivals must produce ONE window record"
        )
        assert entry.representatives == [
            ("phen_seed", "PE"), ("phen_burst_0", "PE"),
            ("phen_burst_1", "PE"), ("phen_burst_2", "PE")]

    def test_spawned_flag_skips_gate_without_new_records(self, parent_polymer):
        """Spec §7.10: after a spawn, a later same-motif arrival does not
        re-run the gate (record count unchanged, no new intent)."""
        model = _GateModel(window=3)

        model.polymer_flux_snapshot = ({}, dict(_PE_STATS), 1.0)
        _gate_pass(model, parent_polymer, _phenolic("phen_1"), iteration=1)
        model.polymer_flux_snapshot = ({"phen_1": 0.5}, dict(_PE_STATS), 0.5)
        _, intents = _gate_pass(model, parent_polymer, _phenolic("phen_2"), iteration=2)
        assert len(intents) == 1
        entry = model.polymer_motif_ledger[0]
        assert entry.spawned is True
        occupancy_before = model.polymer_flux_accumulator.window_occupancy(entry.accumulator_key)

        # In production Phase A would classify this arrival against the new
        # pool; here the registry was not extended, so it reaches the gate
        # and must hit the spawned-flag (assert-log) skip.
        _, intents = _gate_pass(model, parent_polymer, _phenolic("phen_3"), iteration=3)
        assert intents == [], "a spawned motif must never re-spawn"
        assert model.polymer_flux_accumulator.window_occupancy(entry.accumulator_key) == occupancy_before, (
            "the spawned-flag skip must not record (gate not re-run)"
        )

    def test_no_snapshot_defers_honestly(self, parent_polymer):
        """Spec §4.5/§7: no stashed snapshot (iteration 0 path) -> fraction
        0.0 -> defer; the entry still records the zero (zero-filled window),
        and the candidate is absorbed as a proxy variant."""
        model = _GateModel(window=3)
        assert model.polymer_flux_snapshot is None

        processed, intents = _gate_pass(model, parent_polymer, _phenolic("phen_1"), iteration=0)
        assert intents == [], "no snapshot must defer, never fabricate"
        assert len(processed) == 1, "deferred candidate is absorbed as a proxy variant"
        assert getattr(processed[0], "is_polymer_proxy", False) is True
        entry = model.polymer_motif_ledger[0]
        assert model.polymer_flux_accumulator.window_occupancy(entry.accumulator_key) == 1
        assert model.polymer_flux_accumulator.flux(entry.accumulator_key) == 0.0


class TestMotifLedger:
    """Group-isomorphism motif ledger + amended fraction math
    (spec 2026-06-10 §3/§4.3)."""
```

- [ ] **Delete the gamed test.** In `test/rmgpy/polymerMultiPoolTest.py`, Edit (remove the whole method):

  old_string:
```python
    def test_mass_flux_below_threshold_defers_spawn(self, parent_polymer):
        """Even a novel motif should not spawn if its mass flux is below the gate."""
        from rmgpy.polymer import process_polymer_candidates_multipool
        from rmgpy.species import Species

        # Single tiny-flux candidate — should not trigger spawn.
        phenolic = Species(smiles="Oc1ccc(Cc2ccc(Cc3ccc(O)cc3)cc2)cc1")
        processed, intents = process_polymer_candidates_multipool(
            candidates=[phenolic],
            reaction_model=None,
            pool_registry=[parent_polymer],
            mass_flux_threshold=0.99,  # impossible threshold
        )
        assert intents == [], (
            "Mass flux below threshold must defer spawning"
        )

    def test_max_pools_cap_blocks_additional_spawn(self, parent_polymer):
```
  new_string:
```python
    def test_max_pools_cap_blocks_additional_spawn(self, parent_polymer):
```

  (The gamed test "passed" only because the threshold was set above the hardcoded 0.5; its honest replacements are `test_exactly_at_bar_defers_until_nth_recording_iteration` and `test_trace_motif_stays_deferred_across_many_arrivals`.)

- [ ] **Re-baseline `test_novel_motif_spawns_one_pool` to a second sighting.** Edit, replacing the whole method:

  old_string:
```python
    def test_novel_motif_spawns_one_pool(self, parent_polymer):
        from rmgpy.polymer import process_polymer_candidates_multipool
        from rmgpy.species import Species

        # A phenolic 3-mer is structurally distinct from PE — must spawn.
        # NB: a real phenolic 3-mer is a complex adj-list; we use a SMILES proxy
        #     and let the implementation handle motif discovery. The test will
        #     fail concretely with a 'cannot construct chain' message until the
        #     production code knows how to handle this case.
        phenolic_chain = Species(smiles="Oc1ccc(Cc2ccc(Cc3ccc(O)cc3)cc2)cc1")
        processed, intents = process_polymer_candidates_multipool(
            candidates=[phenolic_chain],
            reaction_model=None,
            pool_registry=[parent_polymer],
        )
        assert len(intents) == 1, (
            f"Phenolic novel motif should spawn one pool; got {len(intents)}"
        )
        intent = intents[0]
        assert intent.parent_pool is parent_polymer
        assert intent.monomer is not None
        assert intent.triggering_dp >= 2
```
  new_string:
```python
    def test_novel_motif_spawns_one_pool(self, parent_polymer):
        """Second sighting (spec §7.4 re-baseline): first sighting can never
        spawn, so the ledger is pre-populated with a prior representative
        (recorded with its parent pool, spec §4.3) carrying 50% of gross
        event-mass in a fabricated 3-tuple snapshot — statistic 0.5/3 clears
        the default 0.01 bar at this arrival."""
        from rmgpy.polymer import (MotifLedgerEntry, discover_repeat_motif,
                                   process_polymer_candidates_multipool)
        from rmgpy.species import Species

        # A phenolic 3-mer is structurally distinct from PE — must spawn.
        phenolic_chain = Species(smiles="Oc1ccc(Cc2ccc(Cc3ccc(O)cc3)cc2)cc1")
        phenolic_chain.label = "phenolic_2nd"

        model = _GateModel(window=3)
        motif = discover_repeat_motif(phenolic_chain.molecule[0])
        assert motif is not None
        model.polymer_motif_ledger.append(MotifLedgerEntry(
            motif=motif, accumulator_key="motif-0",
            representatives=[("phenolic_1st", "PE")],
        ))
        # fraction = 0.5/(0.5+0.5) = 0.5.
        model.polymer_flux_snapshot = ({"phenolic_1st": 0.5}, {"PE": (1.0, 1.0)}, 0.5)

        processed, intents = process_polymer_candidates_multipool(
            candidates=[phenolic_chain],
            reaction_model=model,
            pool_registry=[parent_polymer],
            iteration=1,
        )
        assert len(intents) == 1, (
            f"Phenolic novel motif should spawn one pool; got {len(intents)}"
        )
        intent = intents[0]
        assert intent.parent_pool is parent_polymer
        assert intent.monomer is not None
        assert intent.triggering_dp >= 2
        assert intent.mass_flux_at_spawn == pytest.approx(0.5 / 3.0), (
            "mass_flux_at_spawn must carry the REAL gate statistic"
        )
```

- [ ] **Re-baseline `test_iteration_boundary_pass_spawns_and_registers`.** Edit, replacing the section of the method between the parent registration and the spawn-pass call:

  old_string:
```python
        # Synthetic phenolic-trimer candidate; tagged polymer-proxy so the
        # integration helper picks it up.
        phenolic = Species(smiles="Oc1ccc(Cc2ccc(Cc3ccc(O)cc3)cc2)cc1")
        phenolic.is_polymer_proxy = True

        model._apply_multipool_spawn_pass([phenolic])
```
  new_string:
```python
        # Synthetic phenolic-trimer candidate; tagged polymer-proxy so the
        # integration helper picks it up.
        phenolic = Species(smiles="Oc1ccc(Cc2ccc(Cc3ccc(O)cc3)cc2)cc1")
        phenolic.label = "phenolic_2nd"
        phenolic.is_polymer_proxy = True

        # Second sighting (spec §7.4 re-baseline): pre-populate the model's
        # motif ledger and stash a fabricated 3-tuple snapshot in which the
        # first arrival's representative (parent pool PS, recorded at
        # absorption) carries 50% of the gross event-mass -> statistic 0.5/3
        # clears the default 0.01 bar at this arrival.
        from rmgpy.polymer import MotifLedgerEntry, discover_repeat_motif
        motif = discover_repeat_motif(phenolic.molecule[0])
        assert motif is not None
        model.polymer_motif_ledger.append(MotifLedgerEntry(
            motif=motif, accumulator_key="motif-0",
            representatives=[("phenolic_1st", "PS")],
        ))
        model.polymer_flux_snapshot = ({"phenolic_1st": 0.5}, {"PS": (1.0, 1.0)}, 0.5)

        model._apply_multipool_spawn_pass([phenolic])
```

- [ ] **Re-baseline `test_novel_product_spawns_and_registers_and_serializes`.** Edit:

  old_string:
```python
        # A candidate that is structurally novel relative to the PE parent
        # (a phenolic 3-mer surrogate) — should trigger a spawn intent.
        novel = Species(smiles="Oc1ccc(Cc2ccc(Cc3ccc(O)cc3)cc2)cc1")

        processed, intents = process_polymer_candidates_multipool(
            candidates=[novel],
            reaction_model=None,
            pool_registry=[parent_polymer],
        )
        assert len(intents) == 1, "Novel product must produce one spawn intent"
```
  new_string:
```python
        # A candidate that is structurally novel relative to the PE parent
        # (a phenolic 3-mer surrogate). Second sighting (spec §7.4
        # re-baseline): the gate needs a prior representative carrying
        # snapshot flux — first sighting can never spawn.
        novel = Species(smiles="Oc1ccc(Cc2ccc(Cc3ccc(O)cc3)cc2)cc1")
        novel.label = "phenolic_2nd"

        from rmgpy.polymer import MotifLedgerEntry, discover_repeat_motif
        gate_model = _GateModel(window=3)
        motif = discover_repeat_motif(novel.molecule[0])
        assert motif is not None
        gate_model.polymer_motif_ledger.append(MotifLedgerEntry(
            motif=motif, accumulator_key="motif-0",
            representatives=[("phenolic_1st", "PE")],
        ))
        gate_model.polymer_flux_snapshot = ({"phenolic_1st": 0.5}, {"PE": (1.0, 1.0)}, 0.5)

        processed, intents = process_polymer_candidates_multipool(
            candidates=[novel],
            reaction_model=gate_model,
            pool_registry=[parent_polymer],
            iteration=1,
        )
        assert len(intents) == 1, "Novel product must produce one spawn intent"
        assert intents[0].mass_flux_at_spawn == pytest.approx(0.5 / 3.0)
```

- [ ] **Run the test file — verify the new/changed tests FAIL** (gate not flipped yet; the new gate tests find the stub still spawning/deferring wrongly):

```
~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/polymerMultiPoolTest.py -v
```

  Expected: failures in `TestSpawnGateBehavior` (e.g. `test_no_snapshot_defers_honestly` gets an intent because the stub returns 0.5 ≥ 0.01; `test_spike_spawns_at_second_arrival_iteration` fails on the empty-ledger / representatives assertions) and in the re-baselined tests (`mass_flux_at_spawn == 0.5`, not 0.5/3). This is the red phase.

- [ ] **Flip Phase D/E.** In `rmgpy/polymer.py`, three Edits.

  Edit 1 — delete the stub (the whole function):

  old_string:
```python
def _estimate_relative_flux(
    candidate: 'Species',
    pool_registry: List['Polymer'],
    reaction_model: Any,
) -> float:
    """Estimate the fraction of polymer-derived mass flowing into ``candidate``.

    Phase-1 simplification: when ``reaction_model`` is ``None`` (unit-test
    path) return 0.5, so the spawn gate is exercised by the threshold knob
    without needing a full reaction-rate integrator. The real implementation,
    used during an RMG run, will integrate reaction rates and species
    molecular weights over the trailing window via :class:`MassFluxAccumulator`
    (design doc §4.4).
    """
    if reaction_model is None:
        return 0.5
    # TODO(multi-pool §4.4): real flux calculation against reaction_model
    return 0.5


def process_polymer_candidates_multipool(
```
  new_string:
```python
def process_polymer_candidates_multipool(
```

  Edit 2 — accumulator fallback at the top of the function body:

  old_string:
```python
    processed: List['Species'] = []
    spawn_intents: List[SpawnIntent] = []
```
  new_string:
```python
    # The shared accumulator lives on the reaction model (alongside the
    # motif ledger and the stashed snapshot, spec 2026-06-10 §4.3); an
    # explicitly-passed accumulator (test injection) wins.
    if flux_accumulator is None and reaction_model is not None:
        flux_accumulator = getattr(reaction_model, "polymer_flux_accumulator", None)

    processed: List['Species'] = []
    spawn_intents: List[SpawnIntent] = []
```

  Edit 3 — Phase D/E:

  old_string:
```python
        # Phase D: gates — relative flux and max_pools cap.
        relative_flux = _estimate_relative_flux(cand, pool_registry, reaction_model)
        if relative_flux < mass_flux_threshold:
            _tag_polymer_proxy(cand, is_proxy=True)
            processed.append(cand)
            continue
        if len(pool_registry) >= max_pools:
            _tag_polymer_proxy(cand, is_proxy=True)
            processed.append(cand)
            continue

        # Phase E: queue the spawn intent.
        triggering_dp = count_disjoint_subgraph_isomorphisms(mol, motif)
        parent_for_intent = pool_registry[0] if pool_registry else None
        if parent_for_intent is None:
            continue
        spawn_intents.append(
            SpawnIntent(
                parent_pool=parent_for_intent,
                monomer=motif,
                end_groups=list(parent_for_intent.end_groups),
                triggering_product=cand,
                triggering_dp=triggering_dp,
                triggering_moles=float(getattr(cand, "amount", 1.0)),
                mass_flux_at_spawn=relative_flux,
            )
        )
        _tag_polymer_proxy(cand, is_proxy=True)
        processed.append(cand)
```
  new_string:
```python
        # Phase D: gates — mass-flux spawn gate (design doc §4.4, spec
        # 2026-06-10) and max_pools cap.
        spawn_ok, gate_statistic, ledger_entry = _evaluate_spawn_gate(
            cand=cand,
            motif=motif,
            reaction_model=reaction_model,
            iteration=iteration,
            flux_accumulator=flux_accumulator,
            mass_flux_threshold=mass_flux_threshold,
        )
        # Either way, the arriving candidate becomes a representative of its
        # motif (spec §4.4 step 6): it is absorbed as a proxy variant and is
        # the handle a future snapshot attributes flux to. Recorded as a
        # (species_label, parent_pool_label) pair — the parent pool that
        # would parent its spawn intent (currently pool_registry[0]; if
        # multi-parent attribution ever lands, this follows it), per the
        # spec-§3 attribution rules.
        cand_label = getattr(cand, "label", "") or ""
        parent_pool_label = pool_registry[0].label if pool_registry else ""
        if (ledger_entry is not None and cand_label
                and all(lbl != cand_label for (lbl, _) in ledger_entry.representatives)):
            ledger_entry.representatives.append((cand_label, parent_pool_label))
        if not spawn_ok:
            _tag_polymer_proxy(cand, is_proxy=True)
            processed.append(cand)
            continue
        if len(pool_registry) >= max_pools:
            _tag_polymer_proxy(cand, is_proxy=True)
            processed.append(cand)
            continue

        # Phase E: queue the spawn intent.
        triggering_dp = count_disjoint_subgraph_isomorphisms(mol, motif)
        parent_for_intent = pool_registry[0] if pool_registry else None
        if parent_for_intent is None:
            continue
        spawn_intents.append(
            SpawnIntent(
                parent_pool=parent_for_intent,
                monomer=motif,
                end_groups=list(parent_for_intent.end_groups),
                triggering_product=cand,
                triggering_dp=triggering_dp,
                triggering_moles=float(getattr(cand, "amount", 1.0)),
                mass_flux_at_spawn=gate_statistic,
            )
        )
        if ledger_entry is not None:
            ledger_entry.spawned = True
        _tag_polymer_proxy(cand, is_proxy=True)
        processed.append(cand)
```

- [ ] **Wire the reaction model.** In `rmgpy/rmg/model.py`, three Edits.

  Edit 1 — import:

  old_string:
```python
from rmgpy.polymer import Polymer, PolymerCrosslinkError, PolymerFluxArchetype, is_end_group_reaction, stamp_polymer_flux_archetype
```
  new_string:
```python
from rmgpy.polymer import MassFluxAccumulator, Polymer, PolymerCrosslinkError, PolymerFluxArchetype, is_end_group_reaction, stamp_polymer_flux_archetype
```

  Edit 2 — `__init__` attributes:

  old_string:
```python
        self.save_edge_species = False
        self.iteration_num = 0
```
  new_string:
```python
        self.save_edge_species = False
        self.iteration_num = 0
        # Mass-flux spawn-gate state (multi-pool §4.4, spec 2026-06-10): the
        # motif ledger + trailing-window accumulator live on the reaction
        # model. In-memory ONLY — an RMG restart resets windows and deferred
        # motifs re-earn their bar (correct-but-loud, same philosophy as
        # unstamped-reaction demotion). main.py stashes the 3-tuple
        # polymer_flux_snapshot (gross, pool_stats, proxy_event_mass_total)
        # plus the iteration it was taken at right after each polymer
        # simulate(); it stays None for non-polymer systems.
        self.polymer_motif_ledger = []
        self.polymer_flux_accumulator = MassFluxAccumulator()
        self.polymer_flux_snapshot = None
        self.polymer_flux_snapshot_iteration = -1
```

  Edit 3 — pass the iteration into the spawn pass:

  old_string:
```python
        _, intents = process_polymer_candidates_multipool(
            candidates=proxy,
            reaction_model=self,
            pool_registry=pool_registry,
        )
```
  new_string:
```python
        _, intents = process_polymer_candidates_multipool(
            candidates=proxy,
            reaction_model=self,
            pool_registry=pool_registry,
            iteration=getattr(self, "iteration_num", 0),
        )
```

- [ ] **Run — verify green:**

```
~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/polymerMultiPoolTest.py test/rmgpy/rmg/modelTest.py test/rmgpy/polymerTest.py -v
```

  Expected: **all pass**; `polymerMultiPoolTest.py` collects **35 tests** (30 from Task 5 + 6 new − 1 deleted), `modelTest.py` and `polymerTest.py` unchanged from baseline, 0 failures.

- [ ] **Commit:**

```
git add rmgpy/polymer.py rmgpy/rmg/model.py test/rmgpy/polymerMultiPoolTest.py
git commit -m "polymer: replace 0.5 flux stub with real mass-flux spawn gate

Phase D now runs _evaluate_spawn_gate: Group-isomorphism ledger lookup,
(label, parent_pool) pair representatives recorded at absorption, amended
s3 fraction (deduped denominator on proxy_event_mass_total + ledger reps),
one record per motif per iteration, sum/N zero-filled statistic, honest
deferral when no snapshot is stashed. The 0.5 constant, the
reaction_model-is-None special case, and the gamed
test_mass_flux_below_threshold_defers_spawn are DELETED; the three
first-sighting tests are re-baselined to second sightings per the
spawn-floor arithmetic (spec 2026-06-10 s2.4/s7.2-7.5/7.10).

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

## Task 7: main.py — stash the 3-tuple snapshot after each `simulate()`

**Files:**
- `rmgpy/rmg/main.py` (insert before line ~1022, the unique `self.rmg_memories[index].add_t_conv_N(t, x, len(obj))` call)

No new pytest here: the stash is main-loop plumbing exercised end-to-end by the EPDM verification in Task 11 (a polymer deck stashes a real snapshot; non-polymer decks leave it `None`).

- [ ] **Insert the stash.** In `rmgpy/rmg/main.py`, Edit:

  old_string:
```python
                        self.rmg_memories[index].add_t_conv_N(t, x, len(obj))
```
  new_string:
```python
                        # Mass-flux spawn-gate snapshot (multi-pool §4.4, spec
                        # 2026-06-10): read the 3-tuple (gross for all core
                        # species, pool_stats, proxy_event_mass_total) off the
                        # ENGINE — `system.solver`, never the
                        # HybridPolymerReactor blueprint (the established
                        # blueprint-vs-engine gotcha) — and stash on the
                        # reaction model for the Phase-D gate. Stays None for
                        # non-polymer systems (honest degradation: the gate
                        # defers). This stash + the motif ledger are the
                        # shared infrastructure the spec-§2.2 iteration-
                        # boundary re-check upgrade would reuse.
                        engine = getattr(reaction_system, "solver", None) or reaction_system
                        if callable(getattr(engine, "spawn_gate_flux_snapshot", None)):
                            try:
                                self.reaction_model.polymer_flux_snapshot = engine.spawn_gate_flux_snapshot()
                                self.reaction_model.polymer_flux_snapshot_iteration = self.reaction_model.iteration_num
                            except Exception as exc:
                                self.reaction_model.polymer_flux_snapshot = None
                                logging.warning(
                                    "Polymer spawn-gate snapshot failed (all spawns will defer): %s", exc)

                        self.rmg_memories[index].add_t_conv_N(t, x, len(obj))
```

- [ ] **Sanity-check no import/regression breakage:**

```
~/anaconda3/envs/rmg_env/bin/python -c "import rmgpy.rmg.main" && ~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/polymerMultiPoolTest.py test/rmgpy/rmg/modelTest.py -q
```

  Expected: import succeeds; all tests pass (35 + modelTest baseline), 0 failures.

- [ ] **Commit:**

```
git add rmgpy/rmg/main.py
git commit -m "main: stash spawn-gate flux snapshot after each simulate

Computed on the engine (system.solver, not the blueprint) and parked on
reaction_model.polymer_flux_snapshot (+ the iteration it was taken at)
for the Phase-D gate; None for non-polymer systems.

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

## Task 8: ε-direction end-to-end — exhausted pool defers (never inflates)

**Files:**
- `test/rmgpy/solver/solverPolymerTest.py` (append a test to `TestSpawnGateFluxSnapshot`)

This is a regression pin across the engine→gate seam (spec §7.7); it is expected to pass immediately given Tasks 3-6 — write it, confirm it passes, and keep it as the permanent direction tripwire.

- [ ] **Write the test.** Append to `TestSpawnGateFluxSnapshot` (end of `test/rmgpy/solver/solverPolymerTest.py`):

```python
    def test_snapshot_mu0_exhaustion_defers_not_inflates(self):
        """Spec §7.7: mu0 <= SMALL_EPS with tiny mu1 -> pool_stats E[n]
        clamps to 0 -> g_i = 0 -> the gate DEFERS. Asserts the deferral
        DIRECTION, not just finiteness: the naive mu1/mu0 would explode
        toward +inf and wave the motif through. Note the amended split:
        gross itself stays nonzero (it is the raw production record); the
        zeroing lives in pool_stats.
        """
        from rmgpy.polymer import (MassFluxAccumulator, MotifLedgerEntry,
                                   Polymer, discover_repeat_motif,
                                   process_polymer_candidates_multipool)

        sp, core, mask = _one_pool_gate_species()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["A"], sp["R"]], **_KIN)
        rxn.polymer_flux_archetype = 1
        # mu0 exhausted (0 <= SMALL_EPS), mu1 tiny but nonzero.
        rs = _one_pool_gate_rs(rxn, core, mask, (0.0, 1e-25, 1e-20))

        rs.residual(0.0, rs.y, np.zeros_like(rs.y))
        snapshot = rs.spawn_gate_flux_snapshot()
        gross, pool_stats, proxy_total = snapshot

        assert gross["A"] > 0.0, (
            "gross production is nonzero — only the E[n] clamp zeroes g_i"
        )
        assert pool_stats["A"][0] == 0.0, "E[n] must clamp to 0 under SMALL_EPS"
        assert pool_stats["A"][1] == pytest.approx(28.0)
        assert proxy_total == 0.0

        # Feed the engine snapshot to the gate: the exhausted pool must DEFER.
        parent = Polymer(label="PE", monomer="[CH2][CH2]", end_groups=["[H]", "[H]"],
                         cutoff=3, Mn=1000.0, Mw=2500.0, initial_mass=1.0)
        cand = Species(smiles="Oc1ccc(Cc2ccc(Cc3ccc(O)cc3)cc2)cc1")
        cand.label = "phenolic_arrival"
        motif = discover_repeat_motif(cand.molecule[0])
        assert motif is not None

        class _Model:
            pass

        model = _Model()
        model.polymer_motif_ledger = [MotifLedgerEntry(
            motif=motif, accumulator_key="motif-0",
            representatives=[("R", "A")],  # parent pool "A" recorded at absorption
        )]
        model.polymer_flux_snapshot = snapshot

        _, intents = process_polymer_candidates_multipool(
            candidates=[cand],
            reaction_model=model,
            pool_registry=[parent],
            iteration=2,
            flux_accumulator=MassFluxAccumulator(window=3),
        )
        assert intents == [], (
            "a mu0-exhausted pool must defer the spawn (epsilon errs toward deferral)"
        )
```

- [ ] **Run — verify pass:**

```
~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/solver/solverPolymerTest.py -v
```

  Expected: **46 passed** (45 + 1 new). If the new test FAILS, that is a real gate/snapshot bug — debug before proceeding (do not weaken the assertion to finiteness).

- [ ] **Commit:**

```
git add test/rmgpy/solver/solverPolymerTest.py
git commit -m "solver: pin mu0-exhaustion deferral direction end-to-end

pool_stats clamps E[n] to 0 under SMALL_EPS while gross stays a real
record, proxy_event_mass_total zeroes, and the gate defers
(spec 2026-06-10 s7.7 - direction, not just finiteness).

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

## Task 9: `triggering_moles` named-consumer TODO

**Files:**
- `rmgpy/polymer.py` (inside the Phase-E `SpawnIntent(...)` construction)

No behavior change to seeding (spec §5: in scope is the TODO; out of scope is the fix).

- [ ] **Edit.** In `rmgpy/polymer.py`:

  old_string:
```python
                triggering_moles=float(getattr(cand, "amount", 1.0)),
```
  new_string:
```python
                # TODO(polymer running-log item #14): `amount` is never
                # assigned anywhere, so this is ALWAYS the placeholder
                # 1.0 mol. Named consumer: drain_spawn_intents seeds the
                # daughter pool's initial moments from it (mu_k = N * DP^k,
                # see N/DP in drain_spawn_intents) — every gate-path spawned
                # pool currently starts with a fictional mu0 = 1.0 mol of
                # chains. Honest seeding (candidate source: the triggering
                # reaction's snapshot flux x a dt-scale — needs its own
                # chemistry decision) is the NEXT physics item; do not trust
                # mid-run gate-path pool masses until it lands.
                triggering_moles=float(getattr(cand, "amount", 1.0)),
```

- [ ] **Run — verify no breakage:**

```
~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/polymerMultiPoolTest.py -q
```

  Expected: **35 passed**.

- [ ] **Commit:**

```
git add rmgpy/polymer.py
git commit -m "polymer: name the triggering_moles placeholder consumer (item #14)

The 1.0 mol placeholder seeds daughter moments in drain_spawn_intents;
the honest-seeding fix is the next physics item (spec 2026-06-10 s5 -
TODO only, no seeding behavior change).

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

## Task 10: Doc amendment — `docs/multi_pool_design.md` §4.4 rewrite + stale note removal

**Files:**
- `docs/multi_pool_design.md` (§4.4 at lines 149-161; the "NOT YET ACTIVE" note at lines 481-483 — NOTE: the spec calls this the "§8 note", but it actually sits at the end of §10 "Known limitations"; remove it where it lives)

- [ ] **Rewrite §4.4.** Edit:

  old_string:

`````markdown
### 4.4 Mass-flux gate

```
state: per-motif accumulator updated during reaction generation
       motif_flux[motif_signature] = sum of |product moles produced via reactions yielding this motif|
                                   over the trailing N RMG iterations (default N=3)

mass_flux_gate(motif, threshold):
    return motif_flux[motif.signature] / total_polymer_derived_mass >= threshold
```

The signature is a canonical SMARTS or RMG `Group` adjacency-list hash. Trailing-N rolling window prevents a single transient peak from triggering a spawn.
`````

  new_string:

`````markdown
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
`````

- [ ] **Remove the stale note.** Edit:

  old_string:
```markdown
> NOTE (§4.4 mass-flux gate): NOT YET ACTIVE. `_estimate_relative_flux` is a `0.5`
> stub, so spawning is currently gated only by `max_pools` + similarity-merge. The
> rolling-window `MassFluxAccumulator` is implemented but not wired into the live path.

## 11. Open items for follow-up
```
  new_string:
```markdown
## 11. Open items for follow-up
```

- [ ] **Verify the doc edits applied** (both anchors gone):

```
grep -c "NOT YET ACTIVE" docs/multi_pool_design.md && echo "FAIL: note still present" || echo "OK: note removed"
grep -q "Spawn-floor arithmetic" docs/multi_pool_design.md && echo "OK: 4.4 rewritten"
grep -q "proxy_event_mass_total" docs/multi_pool_design.md && echo "OK: amended formula in doc"
```

  Expected: `OK: note removed`, `OK: 4.4 rewritten`, `OK: amended formula in doc`.

- [ ] **Commit:**

```
git add docs/multi_pool_design.md
git commit -m "docs: rewrite multi-pool 4.4 for the live amended mass-flux gate

Arrival-driven rule, amended gross E[n]-calibrated fraction ((label,
parent_pool) representatives, deduped denominator on proxy_event_mass_total
+ ledger reps, stated double-counting decision, epsilon-clamp direction),
one-record-per-iteration + sum/N window, spawn-floor arithmetic, honest
degradation, restart consequence, and the 2.2 upgrade trigger; the stale
NOT-YET-ACTIVE limitation note is removed.

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

## Task 11: Final verification — full suite + EPDM end-to-end no-op + running log

**Files:**
- none in-repo (verification only)
- `~/.claude/projects/-home-alon-Code-RMG-Py/memory/project_polymer_branch_review.md` (running-log update)

- [ ] **Full regression suite:**

```
~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/polymerTest.py test/rmgpy/solver/solverPolymerTest.py test/rmgpy/canteraTest.py test/rmgpy/polymerMultiPoolTest.py test/rmgpy/reactionTest.py test/rmgpy/rmg/modelTest.py test/rmgpy/chemkinTest.py test/rmgpy/polymerArtifactTest.py test/rmgpy/tools/polymerMomentsConsumerTest.py test/rmgpy/tools/polymerMomentsRunnerTest.py -q
```

  Expected: **394 passed, 2 skipped** (baseline 380 − 1 deleted gamed test + 15 new: 2 integrated tripwire + 2 snapshot + 1 ε-direction solver-side; 2 accumulator + 2 ledger + 6 gate behavior multipool-side). No xfails remain.

- [ ] **EPDM end-to-end no-op (spec §7.8 — manual verification, not pytest).** The Task-3 `after_*` timing runs predate the gate flip and main.py stash, so a FRESH run at final HEAD is required. Copy the deck to a scratch dir (never dirty the fixture at `~/runs/RMG/epdm_v0_2026-06-06b/`):

```
mkdir -p ~/runs/RMG/epdm_gate_check_2026-06-10
cp ~/runs/RMG/epdm_v0_2026-06-06b/input.py ~/runs/RMG/epdm_gate_check_2026-06-10/
cd ~/runs/RMG/epdm_gate_check_2026-06-10 && ~/anaconda3/envs/rmg_env/bin/python ~/Code/RMG-Py/rmg.py input.py
```

  Expected: the run finishes with `MODEL GENERATION COMPLETED` (a few minutes). Then verify the no-op:

```
grep -E "MODEL GENERATION COMPLETED|The final model core has" ~/runs/RMG/epdm_gate_check_2026-06-10/RMG.log
~/anaconda3/envs/rmg_env/bin/python -c "import json; d = json.load(open('/home/alon/runs/RMG/epdm_gate_check_2026-06-10/chemkin/polymer_pools.json')); print(d['conventions']['configured_pools'])"
```

  Expected output:
  - `MODEL GENERATION COMPLETED`
  - `The final model core has 26 species and 28 reactions`
  - `['epdm']`

  Rationale (spec decision 4 / §7.8): EPDM's only daughter pool is the Path-C scission handshake, which bypasses this gate AND the `triggering_moles` seeding — so identical counts assert correctness, not a bug-for-bug match. If counts differ, STOP and investigate (the gate must not have touched Path A/C).

- [ ] **Update the running log.** In `~/.claude/projects/-home-alon-Code-RMG-Py/memory/project_polymer_branch_review.md`, find open item **8** (`_estimate_relative_flux` is a `0.5` stub) and append this sentence to the end of that item's text (keeping everything already there; substitute the four measured wall-clock values and ratio from Tasks 2/3):

```
**IMPLEMENTED 2026-06-10 per the AMENDED spec** (plan `docs/superpowers/plans/2026-06-10-mass-flux-spawn-gate.md`, revised post-`cded166bc`): real gate live — change (a) gross writes for ALL core species in the residual (simple.pyx parity; COST GATE measured on the EPDM deck: before <BEFORE_1>/<BEFORE_2>, after <AFTER_1>/<AFTER_2>, ratio <RATIO> ≤ 1.05 accepted), engine `spawn_gate_flux_snapshot()` 3-tuple (gross/pool_stats/proxy_event_mass_total; + `PolymerPoolConfig.monomer_mw_g_mol` plumbed from `PolymerPool.to_config`), main.py snapshot stash, Group-isomorphism motif ledger with (label, parent_pool_label) pair representatives + DEDUPED denominator on CoreEdgeReactionModel, 0.5 stub and gamed test DELETED, 3 first-sighting tests re-baselined. The INTEGRATED tripwire (spec §7.1) was written FIRST and confirmed RED on pre-change HEAD at the gross-array assertion (the born-dead proof), committed xfail(strict=True), then un-xfailed green with change (a) — the §7 standing convention (laundered-quantity fixes ship a live-path test confirmed red first) applied in full. Suite 394 passed / 2 skipped; EPDM no-op verified at final HEAD (26 sp / 28 rxn, configured_pools ['epdm']). The born-dead representative-flux blind spot is CLOSED (representatives now have real gross records); remaining open follow-ups: #14 triggering_moles seeding (TODO comment placed at the construction site), §2.2 iteration-boundary re-check upgrade (trigger documented).
```

- [ ] **No push.** Work stays on the `polymer` branch as local commits.
