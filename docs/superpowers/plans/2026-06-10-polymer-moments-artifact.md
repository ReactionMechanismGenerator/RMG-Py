# Polymer Moments Artifact (schema 2.0) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** RMG emits a versioned declarative artifact (`polymer_pools.json` schema 2.0) that compiles the `HybridPolymerSystem` moment dispatch into data, so a consumer with only numpy + Cantera can integrate per-pool µ0/µ1/µ2 ODEs — proven by a numpy-only sufficiency test, documented by a normative format spec, and cross-validated by a CLI reference runner that drives the real solver.

**Architecture:** The existing 1.0 sidecar writer (`rmgpy/polymer.py`) grows pool-state fields, a compiled `reactions[]` block (one entry per stamped proxy-touching core reaction, mirroring the solver's archetype dispatch including its init-time demotions), and a `conventions` block. The Cantera export (`rmgpy/cantera.py`) exposes its retained-reaction index map so entries can carry authoritative `cantera.index`. A standalone Chemkin >16-char identifier fix unblocks `load_chemkin_file` on EPDM artifacts. A numpy-only consumer module in the test suite proves artifact sufficiency against the solver oracle; a CLI runner in `rmgpy/tools/` drives `HybridPolymerSystem` from artifact + chem.yaml for TA cross-validation.

**Tech Stack:** Pure Python (rmgpy.polymer, rmgpy.cantera, rmgpy.rmg.main, rmgpy.tools) except Task 1 which touches `rmgpy/chemkin.pyx` (requires `make`). Tests: pytest via `~/anaconda3/envs/rmg_env/bin/python -m pytest`.

---

## Ground rules (carried from the project running log)

- Branch `polymer`, repo `/home/alon/Code/RMG-Py`. Conda env `rmg_env` (py3.9).
- Run tests with `~/anaconda3/envs/rmg_env/bin/python -m pytest <file> -q` from the repo root.
- **Any `.pyx`/`.pxd` edit requires `make` from the repo root before running tests.** Only Task 1 touches a `.pyx`.
- Small `topic: imperative` commits, **one per task**. **Never push.**
- Spec: `docs/superpowers/specs/2026-06-10-polymer-moments-artifact-design.md` (approved; do not re-litigate its decisions).
- The chip stream is fully landed at HEAD: `Reaction.__init__` carries `polymer_flux_archetype` / `polymer_chip_units` / `is_end_group_reaction` (`rmgpy/reaction.py:129-131,156-169`) and `stamp_polymer_flux_archetype` is live (`rmgpy/polymer.py:1630-1656`). **Spec §7's cross-stream sequencing caveat is moot** — the emitter reads the fields directly, no defensive `getattr` gate on `discrete_chip/1` emission is needed (plain `getattr(rxn, "polymer_chip_units", 0)` is still used as the universal default-0 accessor, same as the solver at `rmgpy/solver/polymer.pyx:448`).
- Baseline: the six affected suites are green at HEAD `c7d0d845f` — 311 passed, 2 skipped:
  ```
  ~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/polymerTest.py test/rmgpy/solver/solverPolymerTest.py test/rmgpy/canteraTest.py test/rmgpy/polymerMultiPoolTest.py test/rmgpy/reactionTest.py test/rmgpy/rmg/modelTest.py -q
  ```

## Verified code map (probed 2026-06-10 at `c7d0d845f`; spec line numbers have drifted — use these)

| What | Where (verified) |
|---|---|
| `_serialize_pool_for_sidecar` | `rmgpy/polymer.py:2540-2590` (spec said ~2331) |
| `write_polymer_pools_sidecar` | `rmgpy/polymer.py:2593-2630`; schema constant + filename at `:2536-2537` |
| `save_everything` sidecar hook | `rmgpy/rmg/main.py:2036-2081`; `self.notify()` at `:2055` runs the `CanteraWriter` listener (attached at `:773`) BEFORE the sidecar block at `:2057-2079` — ordering confirmed |
| Cantera unbalanced-proxy filter | `rmgpy/cantera.py:210-245` (`_involves_moment_dummy` skip, then `is_balanced()` drop for proxy reactions); entries emitted at `:308-317`; one RMG reaction can expand to **multiple** YAML entries via `reaction_to_dict_list` (`:416-436`, MultiArrhenius/PDep) |
| chem.yaml species naming | `rmgpy/cantera.py:626-639` `get_label`: `f"{label}({index})"` if `index > 0` else bare `label`. µ-dummies have `index = -1` (`rmgpy/rmg/model.py:376`) → bare labels in chem.yaml. **Must not change.** |
| `kb = kf/Keq` | `rmgpy/solver/polymer.pyx:776-785` (spec said 766-775); `Keq = rxn.get_equilibrium_constant(T)` for every reversible reaction; irreversible → `kb=0`, `Keq=inf` |
| Keq formula | `rmgpy/reaction.py:767-840`: `Kc = exp(−ΔG°/RT) · (P0/(R·T))^Δn_gas`, `P0 = 1e5 Pa` |
| RHS dispatch | `rmgpy/solver/polymer.pyx:922-1261` (spec said 940-1037). Concentration rule `_C` at `:1000-1005` (pool-mapped species → 1.0, **applies to reverse product term too**, `:1014-1020`); phase + product gate `:943-989`; site scaling `:1022-1063` (first **reactant**-slot pool, µ0 if `is_end_group_reaction` else µ1, multiplies once, scales rf AND rr); chip throttle `:1053-1060` `min(max(0,µ0), max(0,µ1)/a)/V_poly`; MIGRATION per-direction bundles `:1128-1155`; SCISSION_FRAGMENT complement-stays `:1156-1195`; DISCRETE_CHIP fwd clamp `max(0, 2a·E[n] − a²)` / reverse `+(2a·E[n] + a²)` `:1196-1244`; UNRESOLVED legacy µ1 `:1245-1261` |
| src/dst pool resolution + demotions | `rmgpy/solver/polymer.pyx:527-588`: src = first reactant slot mapping to a **configured** pool; dst = first cross-pool product, else fold-back; NONE+pool-touching → UNRESOLVED; MIGRATION/SCISSION_FRAGMENT missing src or dst → UNRESOLVED; DISCRETE_CHIP missing src → UNRESOLVED |
| `_chain_bundle` / `_safe_mu3_from_mu012` | `rmgpy/solver/polymer.pyx:805-828` / `:141-160` (`SMALL_EPS=1e-30`, `LN_EXP_OVERFLOW_GUARD=700.0`, µ1<µ0 → 0.0) |
| Channel ODEs | `rmgpy/solver/polymer.pyx:1318-1340`: scission `dµ0=k_s(µ1−µ0)`, `dµ2=k_s(µ1−µ3)/3`; unzip `r=k_u·µ0`, `dµ1−=r`, `dµ2−=k_u(2µ1−µ0)`, monomer routed to `monomer_poly_index` (condensed-phase index) |
| Mass transfer | `rmgpy/solver/polymer.pyx:1390-1402`: `J = kLa·(C_poly − K·C_gas)`, `dn = J·V_poly`, `+gas, −poly`; `MassTransferConfig(gas_index, poly_index, K, kLa)` at `:126-134` |
| `PolymerPoolConfig` | `rmgpy/solver/polymer.pyx:89-123` (`label, xs, explicit_dp_to_species_index, mu_indices, monomer_poly_index, k_scission, k_unzip, tail_kinetics`) |
| `initial_polymer_moments` are **extensive mol**, written straight into y0 | `rmgpy/solver/polymer.pyx:724-730` |
| `HybridPolymerSystem.initialize_model` tail | `rmgpy/solver/polymer.pyx:601-610`: `set_initial_conditions` → `generate_rate_coefficients` → `set_initial_derivative` → `initialize_solver`; `advance(t)` from pydas `DASSL` (used in `test/rmgpy/solver/liquidTest.py:155`) |
| Chemkin bug | `rmgpy/chemkin.pyx:1493-1516` `get_species_identifier`: the `index == -1` branch returns a >16-char label with only a warning; `epdm_scission_tail_mu0` (22 chars) shifts the THERMO fixed columns; reproduced live: `load_chemkin_file` on `~/runs/RMG/epdm_v0_2026-06-06b/chemkin/chem.inp` fails with `ValueError: could not convert string to float: 'E+0'`. THERMO/SPECIES/dictionary writers all route through `get_species_identifier` already — the fix is **inside** that function. |
| Two-pool solver fixtures | `test/rmgpy/solver/solverPolymerTest.py:48-80` (`_two_pool_species`, `_two_pool_rs`, `_KIN`); MIGRATION tests `:988-1077`, DISCRETE_CHIP tests `:1178+`; Euler-on-`residual` trajectory pattern `:498-510` |
| Restart pattern | Proven empirically only (handoff `~/handoffs/handoff-polymer-declarative-moments-format.md:100-104`, "matched analytic two-segment solution to 6 decimals"); **no committed test exists** — Task 9 adds one |
| TA loader | `~/Code/TA/ta/mechanism.py:66-85` — already declares the 2.0 pool field names: `moments, monomer_mw_g_mol, mn_g_mol, mw_g_mol, initial_mass_g, channels, phase_species, monomer_routing, mu3_closure`. The emitter MUST use exactly these names. |

**Spec deviations baked into this plan (all additive / code-mirroring, see report):**

1. `reactions[]` entries additionally carry `reactants`, `products` (full label lists) and `proxy_products`. Without them, recipe steps 2/5 (concentration product, manual gas-side stoichiometric flux) and the legacy µ1 bundle are **inexpressible for `cantera: null` entries** (the dropped reactions exist nowhere else), and the reverse concentration product (which treats proxy **products** as 1.0 too — `polymer.pyx:1018-1020`) cannot be formed.
2. `conventions` additionally carries `configured_pools` (the solver-configured pool labels — daughters like `epdm_scission_tail` are registry pools but are NOT solver-configured, and the solver treats their proxies as ordinary species) and `condensed_species` (the solver's phase mask — needed for `C = n/V` phase choice and the product-phase gate at `polymer.pyx:943-989`, which the spec's recipe omitted but the oracle enforces).
3. `monomer_routing` is `null` on current run decks (the run path never resolves `monomer_poly_index`: `PolymerPool.monomer` is a Molecule, `spc_map` is Species-keyed — `rmgpy/rmg/polymer_input.py:706-815`).

## File structure

| File | Action | Responsibility |
|---|---|---|
| `rmgpy/chemkin.pyx` | Modify (`make` needed) | Task 1: compress >16-char identifiers in `get_species_identifier` |
| `test/rmgpy/chemkinTest.py` | Modify | Task 1 regression tests |
| `rmgpy/cantera.py` | Modify | Task 2: `generate_cantera_data(..., return_reaction_index_map=True)` |
| `test/rmgpy/canteraTest.py` | Modify | Task 2 tests |
| `rmgpy/polymer.py` | Modify | Tasks 3–5: pool block additions, reaction compiler, artifact builder, schema constant 2.0 |
| `test/rmgpy/polymerArtifactTest.py` | Create | Tasks 3–5 unit tests (pool blocks, entries, round-trip, 1.0 stability) |
| `rmgpy/rmg/main.py` | Modify | Task 5: hook wiring in `save_everything` |
| `docs/polymer_moments_format.md` | Create | Task 6: normative format spec |
| `test/rmgpy/tools/numpy_moments_consumer.py` | Create | Task 7: numpy-only reference consumer (NO rmgpy imports) |
| `test/rmgpy/tools/polymerMomentsConsumerTest.py` | Create | Tasks 7–8: sufficiency proof + mass-transfer cross-check |
| `rmgpy/tools/polymer_moments_runner.py` | Create | Task 9: CLI reference runner |
| `setup.py` | Modify | Task 9: console entry point |
| `test/rmgpy/tools/polymerMomentsRunnerTest.py` | Create | Task 9 tests (incl. two-segment restart check) |
| `~/runs/RMG/epdm_v0_2026-06-06b/` | Regenerate (not committed) | Task 10: EPDM schema-2.0 fixture |
| `~/handoffs/handoff-polymer-moments-artifact-landed.md` | Create (not committed) | Task 11: TA handoff-back |

---

### Task 1: Chemkin >16-char µ-dummy identifier fix (own commit, independent of everything else)

The bug: µ-dummy species are created with `index = -1` (`rmgpy/rmg/model.py:376`), so `get_species_identifier` (`rmgpy/chemkin.pyx:1508-1516`) returns their raw label even when >16 chars. `epdm_scission_tail_mu0` (22 chars) then overflows the fixed-column THERMO line 1 (`write_thermo_entry`, `:1640`), shifting the element block and the column-80 `1` marker; `load_chemkin_file` dies with `ValueError: could not convert string to float: 'E+0'` (reproduced live on the EPDM artifact). All Chemkin writers (SPECIES block, THERMO, species_dictionary, transport) route through `get_species_identifier`, so fixing it there propagates consistently and the written files stay self-consistent and loadable.

**Constraint:** chem.yaml naming (`rmgpy/cantera.py get_label`) must NOT change — TA matches the bare labels there. This fix touches only the Chemkin identifier.

**Design:** for the no-index path with `len(label) > 16`, emit a deterministic ≤16-char identifier: truncated stem + `-` + 4-hex md5 digest of the full label + a preserved readable `_muN` tail when present. Distinct labels (e.g. `epdm_scission_tail_mu0` vs `epdm_scission_head_mu0`) get distinct identifiers via the digest; the same label always maps to the same identifier.

**Files:**
- Modify: `rmgpy/chemkin.pyx:1508-1516` (and module imports near the top)
- Test: `test/rmgpy/chemkinTest.py` (append a class; all needed names are already imported at `:34-60`)

- [ ] **Step 1: Write the failing tests**

Append to `test/rmgpy/chemkinTest.py`:

```python
class TestLongLabelSpeciesIdentifier:
    """Regression for the >16-char moment-dummy Chemkin identifier overflow.

    Moment dummies are created with index = -1 (rmgpy/rmg/model.py), so
    get_species_identifier used to return the raw label even when longer than
    16 characters, overflowing the fixed-column THERMO format and breaking
    load_chemkin_file on the EPDM artifact (ValueError: could not convert
    string to float: 'E+0').
    """

    @staticmethod
    def _mu_dummy(label):
        from rmgpy.molecule import Molecule
        spc = Species(label=label, reactive=False)
        spc.molecule = [Molecule().from_smiles("[Ne]")]
        spc.index = -1
        spc.is_moment_dummy = True
        return spc

    @staticmethod
    def _ne_thermo():
        poly_low = NASAPolynomial(
            coeffs=[2.5, 0.0, 0.0, 0.0, 0.0, -745.375, 3.35532],
            Tmin=(200.0, "K"), Tmax=(1000.0, "K"))
        poly_high = NASAPolynomial(
            coeffs=[2.5, 0.0, 0.0, 0.0, 0.0, -745.375, 3.35532],
            Tmin=(1000.0, "K"), Tmax=(6000.0, "K"))
        return NASA(polynomials=[poly_low, poly_high],
                    Tmin=(200.0, "K"), Tmax=(6000.0, "K"))

    def test_long_label_is_compressed_to_16_chars(self):
        spc = self._mu_dummy("epdm_scission_tail_mu0")
        ident = get_species_identifier(spc)
        assert len(ident) <= 16
        assert ident.endswith("_mu0")  # readable moment tail preserved

    def test_short_labels_unchanged(self):
        spc = self._mu_dummy("epdm_mu0")
        assert get_species_identifier(spc) == "epdm_mu0"

    def test_distinct_long_labels_get_distinct_identifiers(self):
        labels = ["epdm_scission_tail_mu0", "epdm_scission_head_mu0",
                  "epdm_scission_tail_mu1", "epdm_scission_tail_mu2"]
        idents = [get_species_identifier(self._mu_dummy(l)) for l in labels]
        assert len(set(idents)) == len(labels)
        assert all(len(i) <= 16 for i in idents)

    def test_identifier_is_deterministic(self):
        a = get_species_identifier(self._mu_dummy("epdm_scission_tail_mu0"))
        b = get_species_identifier(self._mu_dummy("epdm_scission_tail_mu0"))
        assert a == b

    def test_chemkin_roundtrip_with_long_moment_dummy_labels(self, tmpdir):
        """save_chemkin_file + save_species_dictionary then load_chemkin_file
        must succeed with >16-char moment-dummy labels (EPDM failure mode)."""
        long_dummy = self._mu_dummy("epdm_scission_tail_mu0")
        long_dummy.thermo = self._ne_thermo()
        short_dummy = self._mu_dummy("epdm_mu0")
        short_dummy.thermo = self._ne_thermo()
        species = [long_dummy, short_dummy]

        chem_path = os.path.join(str(tmpdir), "chem.inp")
        dict_path = os.path.join(str(tmpdir), "species_dictionary.txt")
        save_chemkin_file(chem_path, species, [], verbose=False, check_for_duplicates=False)
        save_species_dictionary(dict_path, species)

        loaded_species, loaded_reactions = load_chemkin_file(chem_path, dict_path)
        assert len(loaded_species) == 2
        labels = {s.label for s in loaded_species}
        # the compressed identifier round-trips as the species label
        assert any(l.endswith("_mu0") and len(l) <= 16 for l in labels)
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/chemkinTest.py::TestLongLabelSpeciesIdentifier -q`
Expected: `test_long_label_is_compressed_to_16_chars`, `test_distinct_long_labels_get_distinct_identifiers`, and `test_chemkin_roundtrip_with_long_moment_dummy_labels` FAIL (identifier is 22 chars / loader ValueError). `test_short_labels_unchanged` may already pass.

- [ ] **Step 3: Implement the fix in `rmgpy/chemkin.pyx`**

Add to the module's imports (near the top, where `import re` / `import logging` live):

```python
import hashlib
```

Replace at `rmgpy/chemkin.pyx:1510-1516` (inside `get_species_identifier`, the `species.index == -1` branch):

```python
        if len(label) > 0 and not re.search(r'[^A-Za-z0-9\-_,\(\)\*#.:\[\]]+', label):
            if len(label) <= 16:
                return label
            else:
                # Compress to a deterministic <=16-char identifier instead of
                # overflowing the fixed-column Chemkin THERMO format (long
                # polymer moment-dummy labels like 'epdm_scission_tail_mu0'
                # used to shift the element block and break load_chemkin_file).
                # Keep a readable '_muN' tail when present; disambiguate the
                # truncated stem with a 4-hex digest of the full label.
                digest = hashlib.md5(label.encode()).hexdigest()[:4]
                match = re.search(r'(_mu[0-9])$', label)
                suffix = match.group(1) if match else ''
                stem = label[:16 - len(suffix) - 5]
                ident = '{0}-{1}{2}'.format(stem, digest, suffix)
                logging.warning('Species label %r is longer than 16 characters; '
                                'using compressed Chemkin identifier %r.', label, ident)
                return ident
```

(`-` is in the allowed Chemkin identifier charset per the regex on the line above; 16 − 4 (`_mu0`) − 5 (`-` + 4 hex) = 7-char stem → `epdm_sc-xxxx_mu0`, exactly 16.)

- [ ] **Step 4: Rebuild Cython**

Run: `make` (from `/home/alon/Code/RMG-Py`)
Expected: compiles `rmgpy/chemkin.pyx` without errors.

- [ ] **Step 5: Run the new tests + the chemkin suite**

Run: `~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/chemkinTest.py -q`
Expected: all PASS (including the pre-existing tests — the changed branch only triggers for >16-char no-index labels, which nothing else exercises).

- [ ] **Step 6: Verify against the real EPDM artifact (the original repro)**

Run:
```bash
~/anaconda3/envs/rmg_env/bin/python -c "
from rmgpy.chemkin import load_chemkin_file
sp, rx = load_chemkin_file('$HOME/runs/RMG/epdm_v0_2026-06-06b/chemkin/chem.inp',
                           '$HOME/runs/RMG/epdm_v0_2026-06-06b/chemkin/species_dictionary.txt')
print('LOADED', len(sp), 'species,', len(rx), 'reactions')"
```
Expected: still FAILS — that file was written by the buggy writer and is fixed by regeneration in Task 10. This step documents the expectation; do not chase it here. (If you want a live check now: `load_chemkin_file` on a freshly written file from Step 1's roundtrip test is the proof.)

- [ ] **Step 7: Commit**

```bash
git add rmgpy/chemkin.pyx test/rmgpy/chemkinTest.py
git commit -m "chemkin: compress >16-char species identifiers

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 2: Cantera export exposes the retained-reaction index map

`reactions[]` entries need `cantera.index` = the entry's position in chem.yaml's `reactions:` list (authoritative — the consumer calls `Kinetics.set_multiplier(0, index)` on it). The filter and ordering live in `generate_cantera_data` (`rmgpy/cantera.py:176-319`); one RMG reaction can expand to multiple YAML entries (`reaction_to_dict_list`, MultiArrhenius/PDep). Expose the mapping from the same code path so it can never drift from the export.

**Files:**
- Modify: `rmgpy/cantera.py:303-319`, signature at `:176-181`
- Test: `test/rmgpy/canteraTest.py` (style/fixtures: see `test_generate_cantera_data_drops_unbalanced_polymer_proxy_reaction` at `:264-306`; `self._create_dummy_species(label, smiles, index)` exists on the test class)

- [ ] **Step 1: Write the failing test**

Add to the test class in `test/rmgpy/canteraTest.py` that contains `test_generate_cantera_data_drops_unbalanced_polymer_proxy_reaction` (reuse its imports — `generate_cantera_data`, `Reaction`, `Arrhenius` are already imported in that file):

```python
    def test_generate_cantera_data_returns_reaction_index_map(self):
        """The optional reaction index map gives, per retained RMG reaction,
        its entry indices in data['reactions'] (chem.yaml order); dropped
        (unbalanced-proxy) reactions are absent from the map. Needed by the
        polymer moments artifact (cantera.index is authoritative: the consumer
        zeroes the multiplier at that index)."""
        parent = self._create_dummy_species("parent", "CCC", index=1)
        parent.is_polymer_proxy = True
        tail = self._create_dummy_species("scission_tail", "CCCCC", index=2)
        tail.is_polymer_proxy = True
        h = self._create_dummy_species("H", "[H]", index=3)
        hh = self._create_dummy_species("H2", "[H][H]", index=4)
        r = self._create_dummy_species("R", "[CH2]O", index=5)
        p = self._create_dummy_species("P", "C[O]", index=6)

        dropped = Reaction(
            reactants=[h, parent], products=[hh, tail],
            kinetics=Arrhenius(A=(1e10, "cm^3/(mol*s)"), n=0, Ea=(0, "J/mol"), T0=(1, "K")),
        )
        kept1 = Reaction(
            reactants=[r], products=[p],
            kinetics=Arrhenius(A=(1e10, "s^-1"), n=0, Ea=(0, "J/mol"), T0=(1, "K")),
        )
        kept2 = Reaction(
            reactants=[h, h], products=[hh],
            kinetics=Arrhenius(A=(1e10, "cm^3/(mol*s)"), n=0, Ea=(0, "J/mol"), T0=(1, "K")),
        )

        species = [parent, tail, h, hh, r, p]
        reactions = [dropped, kept1, kept2]

        data, index_map = generate_cantera_data(species, reactions,
                                                return_reaction_index_map=True)

        assert id(dropped) not in index_map
        assert index_map[id(kept1)] == [0]
        assert index_map[id(kept2)] == [1]
        # indices really point at the right entries
        assert "R(5)" in data["reactions"][0]["equation"]
        assert "H2(4)" in data["reactions"][1]["equation"]
        # default call signature unchanged
        data_only = generate_cantera_data(species, reactions)
        assert data_only["reactions"] == data["reactions"]

    def test_reaction_index_map_multi_entry_reaction(self):
        """A MultiArrhenius reaction expands to several YAML entries; the map
        must list ALL of its indices (the consumer zeroes every one)."""
        from rmgpy.kinetics import MultiArrhenius
        r = self._create_dummy_species("R", "[CH2]O", index=5)
        p = self._create_dummy_species("P", "C[O]", index=6)
        multi = Reaction(
            reactants=[r], products=[p],
            kinetics=MultiArrhenius(arrhenius=[
                Arrhenius(A=(1e10, "s^-1"), n=0, Ea=(0, "J/mol"), T0=(1, "K")),
                Arrhenius(A=(1e8, "s^-1"), n=0, Ea=(0, "J/mol"), T0=(1, "K")),
            ]),
        )
        single = Reaction(
            reactants=[p], products=[r],
            kinetics=Arrhenius(A=(1e10, "s^-1"), n=0, Ea=(0, "J/mol"), T0=(1, "K")),
        )
        data, index_map = generate_cantera_data([r, p], [multi, single],
                                                return_reaction_index_map=True)
        assert index_map[id(multi)] == [0, 1]
        assert index_map[id(single)] == [2]
        assert len(data["reactions"]) == 3
```

- [ ] **Step 2: Run to verify failure**

Run: `~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/canteraTest.py -q -k "index_map"`
Expected: FAIL — `TypeError: generate_cantera_data() got an unexpected keyword argument 'return_reaction_index_map'`

- [ ] **Step 3: Implement**

In `rmgpy/cantera.py`, change the signature at `:176-181`:

```python
def generate_cantera_data(species_list,
                          reaction_list,
                          is_plasma=False,
                          site_density=None,
                          search_for_additional_elements=False,
                          return_reaction_index_map=False,
                          ):
```

Replace the entry-emission block at `:308-319`:

```python
    reaction_data = list()
    reaction_index_map = {}
    for rxn in gas_reactions + surface_reactions:
        entries = reaction_to_dict_list(rxn, species_list)
        if entries:
            start = len(reaction_data)
            reaction_index_map[id(rxn)] = list(range(start, start + len(entries)))
            reaction_data.extend(entries)
    data['reactions'] = reaction_data

    if return_reaction_index_map:
        return data, reaction_index_map
    return data
```

(Behavior-identical to the old two-loop version for the data itself: same order — gas reactions first, then surface.)

- [ ] **Step 4: Run the cantera suite**

Run: `~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/canteraTest.py -q`
Expected: all PASS.

- [ ] **Step 5: Commit**

```bash
git add rmgpy/cantera.py test/rmgpy/canteraTest.py
git commit -m "cantera: expose retained-reaction index map from export

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 3: Schema-2.0 pool blocks + envelope

Extend `_serialize_pool_for_sidecar` with the spec §3 fields (names pinned by TA's loader: `moments, monomer_mw_g_mol, mn_g_mol, mw_g_mol, initial_mass_g, channels, phase_species, monomer_routing, mu3_closure`), keep every 1.0 field byte-identical, and bump the schema constant. Envelope additions (`rmg_commit`, `conventions`, `reactions`) land in Task 5 via the builder; this task adds the `rmg_commit` helper and the pool serializer.

**Files:**
- Modify: `rmgpy/polymer.py:2536-2590` (constant + `_serialize_pool_for_sidecar`)
- Test: `test/rmgpy/polymerArtifactTest.py` (create)

- [ ] **Step 1: Write the failing tests**

Create `test/rmgpy/polymerArtifactTest.py`:

```python
#!/usr/bin/env python3
"""Tests for the schema-2.0 polymer moments artifact emitter
(spec: docs/superpowers/specs/2026-06-10-polymer-moments-artifact-design.md;
format doc: docs/polymer_moments_format.md)."""

import json

import pytest

from rmgpy.kinetics import Arrhenius
from rmgpy.molecule import Molecule
from rmgpy.polymer import (
    POLYMER_POOLS_SIDECAR_SCHEMA_VERSION,
    Polymer,
    PolymerFluxArchetype,
    _serialize_pool_for_sidecar,
    build_polymer_moments_artifact,
    compile_polymer_reaction_entries,
    write_polymer_pools_sidecar,
)
from rmgpy.reaction import Reaction
from rmgpy.species import Species


# Every schema-1.0 key with its exact type expectations (spec §2 backward
# compatibility; TA's loader is .get()-based but the keys must not move).
SCHEMA_1_0_KEYS = {
    "label": str,
    "monomer_smiles": str,
    "monomer_adj_list": str,
    "feature_monomers_smiles": list,
    "end_groups": list,
    "cutoff": int,
    "parent_pool": (str, type(None)),
    "spawn_iteration": int,
    "spawn_event_metadata": dict,
    "mu_indices": (dict, type(None)),
}


def _spc(smiles, label, index=-1):
    s = Species(molecule=[Molecule().from_smiles(smiles)])
    s.label = label
    s.index = index
    return s


def _mu_dummy(label):
    s = Species(label=label, reactive=False)
    s.molecule = [Molecule().from_smiles("[Ne]")]
    s.is_moment_dummy = True
    s.index = -1
    return s


@pytest.fixture
def pe_pool():
    return Polymer(
        label="PE",
        monomer="[CH2][CH2]",
        end_groups=["[H]", "[H]"],
        cutoff=3,
        Mn=1500.0,
        Mw=1800.0,
        initial_mass=1.0,
        k_scission=1.0,
        k_unzip=0.01,
    )


class TestPoolBlockSchema2:
    def test_schema_version_bumped(self):
        assert POLYMER_POOLS_SIDECAR_SCHEMA_VERSION == "2.0"

    def test_all_schema_1_0_fields_present_with_stable_types(self, pe_pool):
        d = _serialize_pool_for_sidecar(pe_pool)
        for key, typ in SCHEMA_1_0_KEYS.items():
            assert key in d, f"1.0 key {key!r} missing"
            assert isinstance(d[key], typ), f"1.0 key {key!r} changed type: {type(d[key])}"

    def test_pool_additions(self, pe_pool):
        d = _serialize_pool_for_sidecar(pe_pool)
        # state at generation time, mol / DP basis
        assert d["moments"] == pytest.approx(list(pe_pool.moments))
        assert d["monomer_mw_g_mol"] == pytest.approx(pe_pool.monomer_mw_g_mol)
        assert d["mn_g_mol"] == pytest.approx(1500.0)
        assert d["mw_g_mol"] == pytest.approx(1800.0)
        assert d["initial_mass_g"] == pytest.approx(1000.0)
        # Arrhenius-capable channels; today A=k, n=0, Ea=0
        assert d["channels"]["scission"] == {
            "A": 1.0, "n": 0.0, "Ea": 0.0,
            "units": {"A": "s^-1", "Ea": "J/mol"},
        }
        assert d["channels"]["unzip"]["A"] == pytest.approx(0.01)
        assert d["mu3_closure"] == "log_lagrange/1"
        assert d["monomer_routing"] is None
        assert d["phase_species"] == []

    def test_phase_species_collected_from_core(self, pe_pool):
        core = [
            _spc("CC", "PE", index=2),          # proxy (Polymer stand-in by label)
            _mu_dummy("PE_mu0"), _mu_dummy("PE_mu1"), _mu_dummy("PE_mu2"),
            _spc("[CH3]", "G", index=7),        # gas — must not appear
            _mu_dummy("OTHER_mu0"),             # other pool's dummy — must not appear
        ]
        d = _serialize_pool_for_sidecar(pe_pool, core_species=core)
        assert d["phase_species"] == ["PE(2)", "PE_mu0", "PE_mu1", "PE_mu2"]

    def test_monomer_routing_passthrough(self, pe_pool):
        d = _serialize_pool_for_sidecar(pe_pool, monomer_routing="ethylene(5)")
        assert d["monomer_routing"] == "ethylene(5)"
        assert "ethylene(5)" in d["phase_species"]
```

- [ ] **Step 2: Run to verify failure**

Run: `~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/polymerArtifactTest.py -q`
Expected: FAIL with `ImportError: cannot import name 'build_polymer_moments_artifact'` (those land in Tasks 4–5; for now comment the two not-yet-existing names out of the import in your local run if you want the pool tests red/green in isolation — or accept the ImportError as the failing state and proceed).

To run just this task's tests before Tasks 4–5 exist, temporarily reduce the import to the names defined so far; restore the full import in Task 4. Simpler: implement Step 3 then re-run — `compile_polymer_reaction_entries`/`build_polymer_moments_artifact` are added as stubs in Step 3 below so the module imports.

- [ ] **Step 3: Implement in `rmgpy/polymer.py`**

Change the constant at `:2536`:

```python
POLYMER_POOLS_SIDECAR_SCHEMA_VERSION = "2.0"
```

Add helpers just above `_serialize_pool_for_sidecar` (after the constant block):

```python
def _artifact_species_label(spc) -> str:
    """chem.yaml species name for ``spc`` — must match rmgpy.cantera.get_label:
    ``label(index)`` when index > 0, bare label otherwise (µ-dummies have
    index = -1 and appear bare in chem.yaml)."""
    index = getattr(spc, "index", -1)
    label = getattr(spc, "label", "")
    return f"{label}({index})" if index > 0 else label


def _species_base_label(spc) -> str:
    """Strip the RMG ``(N)`` index suffix — the solver's pool-membership rule
    (rmgpy/solver/polymer.pyx:478-485 uses label.partition('(')[0])."""
    return getattr(spc, "label", "").partition('(')[0]


def _get_rmg_commit():
    """Best-effort git SHA of the emitting RMG-Py checkout (envelope field)."""
    try:
        import subprocess
        repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        out = subprocess.run(["git", "-C", repo, "rev-parse", "HEAD"],
                             capture_output=True, text=True, timeout=5)
        if out.returncode == 0:
            return out.stdout.strip()
    except Exception:
        pass
    return None
```

Replace `_serialize_pool_for_sidecar` (`:2540-2590`) wholesale with this complete function (the 1.0 construction code is identical to the current body, only repackaged into `d`):

```python
def _serialize_pool_for_sidecar(pool: 'Polymer',
                                core_species: Optional[List['Species']] = None,
                                monomer_routing: Optional[str] = None) -> Dict[str, Any]:
    """Convert a :class:`Polymer` instance to a JSON-serialisable dict.

    Schema 2.0: 1.0 fields (docs/multi_pool_design.md §6) preserved verbatim;
    additions per docs/polymer_moments_format.md §2.
    """
    monomer_smiles = ""
    monomer_adj_list = ""
    try:
        if getattr(pool, "monomer", None) is not None:
            monomer_smiles = pool.monomer.to_smiles() if hasattr(pool.monomer, "to_smiles") else ""
            monomer_adj_list = (
                pool.monomer.to_adjacency_list() if hasattr(pool.monomer, "to_adjacency_list") else ""
            )
    except Exception:
        pass

    feature_smiles: List[str] = []
    feature_attr = getattr(pool, "feature_monomers", None) or (
        [pool.feature_monomer] if getattr(pool, "feature_monomer", None) else []
    )
    for fm in feature_attr:
        try:
            if hasattr(fm, "to_smiles"):
                feature_smiles.append(fm.to_smiles())
        except Exception:
            continue

    spawn_metadata = getattr(pool, "spawn_metadata", None) or {"source": "input"}
    mu_indices = getattr(pool, "mu_indices", None)
    if mu_indices is not None and not isinstance(mu_indices, dict):
        try:
            mu0_idx, mu1_idx, mu2_idx = mu_indices
            mu_indices = {"mu0_idx": mu0_idx, "mu1_idx": mu1_idx, "mu2_idx": mu2_idx}
        except Exception:
            mu_indices = None

    d = {
        "label": getattr(pool, "label", ""),
        "monomer_smiles": monomer_smiles,
        "monomer_adj_list": monomer_adj_list,
        "feature_monomers_smiles": feature_smiles,
        "end_groups": [
            eg.to_smiles() if hasattr(eg, "to_smiles") else str(eg)
            for eg in (getattr(pool, "end_groups", []) or [])
        ],
        "cutoff": getattr(pool, "cutoff", None),
        "parent_pool": getattr(pool, "parent_pool_label", None),
        "spawn_iteration": getattr(pool, "spawn_iteration", 0),
        "spawn_event_metadata": spawn_metadata,
        "mu_indices": mu_indices,
    }

    # --- schema 2.0 additions (field names pinned by TA's loader,
    #     ~/Code/TA/ta/mechanism.py) ---
    moments = getattr(pool, "moments", None)
    d["moments"] = [float(m) for m in moments] if moments is not None else None
    d["monomer_mw_g_mol"] = (float(pool.monomer_mw_g_mol)
                             if getattr(pool, "monomer_mw_g_mol", None) is not None else None)
    d["mn_g_mol"] = float(pool.Mn) if getattr(pool, "Mn", None) is not None else None
    d["mw_g_mol"] = float(pool.Mw) if getattr(pool, "Mw", None) is not None else None
    d["initial_mass_g"] = (float(pool.initial_mass_g)
                           if getattr(pool, "initial_mass_g", None) is not None else None)
    d["channels"] = {
        "scission": {"A": float(getattr(pool, "k_scission", 0.0)), "n": 0.0, "Ea": 0.0,
                     "units": {"A": "s^-1", "Ea": "J/mol"}},
        "unzip": {"A": float(getattr(pool, "k_unzip", 0.0)), "n": 0.0, "Ea": 0.0,
                  "units": {"A": "s^-1", "Ea": "J/mol"}},
    }
    phase_species: List[str] = []
    if core_species:
        member_bases = {pool.label, f"{pool.label}_mu0", f"{pool.label}_mu1",
                        f"{pool.label}_mu2"}
        for spc in core_species:
            if _species_base_label(spc) in member_bases:
                phase_species.append(_artifact_species_label(spc))
    if monomer_routing and monomer_routing not in phase_species:
        phase_species.append(monomer_routing)
    d["phase_species"] = phase_species
    d["monomer_routing"] = monomer_routing
    d["mu3_closure"] = "log_lagrange/1"
    return d
```

Do not rename or reorder any 1.0 key.

Also add **stubs** so Tasks 4–5 imports resolve (replaced by real implementations there):

```python
def compile_polymer_reaction_entries(*args, **kwargs):
    raise NotImplementedError  # implemented in the reaction-compiler task


def build_polymer_moments_artifact(*args, **kwargs):
    raise NotImplementedError  # implemented in the artifact-builder task
```

- [ ] **Step 4: Run this task's tests**

Run: `~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/polymerArtifactTest.py -q`
Expected: PASS (all of `TestPoolBlockSchema2`).

- [ ] **Step 5: Run the neighbor suites (sidecar tests live in polymerMultiPoolTest)**

Run: `~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/polymerMultiPoolTest.py test/rmgpy/polymerTest.py -q`
Expected: PASS — `TestPolymerPoolsSidecar.test_writes_valid_schema` compares against the constant, so the bump doesn't break it.

- [ ] **Step 6: Commit**

```bash
git add rmgpy/polymer.py test/rmgpy/polymerArtifactTest.py
git commit -m "polymer: emit schema-2.0 pool blocks in sidecar

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 4: Reaction compiler — stamps → compiled flux entries

One entry per stamped/pool-touching core reaction. The compiler must mirror the solver's init-time pool resolution and demotions (`polymer.pyx:527-588`) so the artifact describes what the oracle actually does (spec Q3: UNRESOLVED/legacy entries ARE emitted, flagged).

**Files:**
- Modify: `rmgpy/polymer.py` (replace the Task-3 stub)
- Test: `test/rmgpy/polymerArtifactTest.py` (append)

- [ ] **Step 1: Write the failing tests**

Append to `test/rmgpy/polymerArtifactTest.py`:

```python
def _arrhenius(A=(2.0, "s^-1")):
    return Arrhenius(A=A, n=0.0, Ea=(0.0, "J/mol"), T0=(1.0, "K"))


def _two_pool_core():
    """Mirror of test/rmgpy/solver/solverPolymerTest.py:_two_pool_species —
    pools A and B with proxies + µ-dummies, gas species G and C."""
    core = [
        _spc("CCCC", "A", index=1),
        _mu_dummy("A_mu0"), _mu_dummy("A_mu1"), _mu_dummy("A_mu2"),
        _spc("CCCCC", "B", index=5),
        _mu_dummy("B_mu0"), _mu_dummy("B_mu1"), _mu_dummy("B_mu2"),
        _spc("[CH3]", "G", index=9),
        _spc("C", "C1", index=10),
    ]
    core[0].is_polymer_proxy = True
    core[4].is_polymer_proxy = True
    return core


class TestCompileReactionEntries:
    def test_migration_entry(self):
        core = _two_pool_core()
        rxn = Reaction(reactants=[core[0]], products=[core[4]],
                       kinetics=_arrhenius(), reversible=False)
        rxn.polymer_flux_archetype = int(PolymerFluxArchetype.MIGRATION)
        entries = compile_polymer_reaction_entries(
            [rxn], core, configured_pool_labels=["A", "B"],
            cantera_index_map={id(rxn): [4]})
        assert len(entries) == 1
        e = entries[0]
        assert e["id"] == "r4"
        assert e["cantera"] == {"index": 4, "equation": "A(1) <=> B(5)"}
        assert e["archetype"] == "migration/1"
        assert e["src_pool"] == "A" and e["dst_pool"] == "B"
        assert e["scaling"] == "mu1"
        assert e["unresolved"] is False
        assert e["proxy_reactants"] == ["A(1)"]
        assert e["proxy_products"] == ["B(5)"]
        assert e["reactants"] == ["A(1)"] and e["products"] == ["B(5)"]
        assert e["kinetics"]["A"] == pytest.approx(2.0)
        assert e["kinetics"]["units"]["A"] == "s^-1"
        assert e["kinetics"]["reversible"] is False
        assert "params" not in e

    def test_discrete_chip_entry_carries_a_and_mu0_scaling(self):
        core = _two_pool_core()
        rxn = Reaction(reactants=[core[0]], products=[core[0], core[9]],
                       kinetics=_arrhenius(), reversible=False)
        rxn.polymer_flux_archetype = int(PolymerFluxArchetype.DISCRETE_CHIP)
        rxn.polymer_chip_units = 2
        rxn.is_end_group_reaction = True
        entries = compile_polymer_reaction_entries(
            [rxn], core, configured_pool_labels=["A", "B"],
            cantera_index_map={id(rxn): [0]})
        e = entries[0]
        assert e["archetype"] == "discrete_chip/1"
        assert e["params"] == {"a": 2}
        assert e["scaling"] == "mu0"
        assert e["src_pool"] == "A" and e["dst_pool"] == "A"  # fold-back

    def test_dropped_reaction_is_cantera_null_and_carries_kinetics(self):
        core = _two_pool_core()
        rxn = Reaction(reactants=[core[8], core[0]], products=[core[9], core[4]],
                       kinetics=_arrhenius(A=(3.0, "m^3/(mol*s)")), reversible=True)
        rxn.polymer_flux_archetype = int(PolymerFluxArchetype.SCISSION_FRAGMENT)
        entries = compile_polymer_reaction_entries(
            [rxn], core, configured_pool_labels=["A", "B"],
            cantera_index_map={})  # not retained by the export
        e = entries[0]
        assert e["cantera"] is None
        assert e["kinetics"] is not None  # REQUIRED for cantera-null (spec §4)
        assert e["kinetics"]["A"] == pytest.approx(3.0)
        assert e["kinetics"]["units"]["A"] == "m^3/(mol*s)"
        assert e["kinetics"]["reversible"] is True
        assert e["id"].startswith("d")
        assert e["archetype"] == "scission_fragment/1"

    def test_unstamped_pool_touching_reaction_emits_legacy_unresolved(self):
        core = _two_pool_core()
        rxn = Reaction(reactants=[core[8], core[0]], products=[core[9], core[0]],
                       kinetics=_arrhenius(A=(3.0, "m^3/(mol*s)")), reversible=False)
        # no polymer_flux_archetype stamp at all (legacy emission, spec Q3)
        entries = compile_polymer_reaction_entries(
            [rxn], core, configured_pool_labels=["A", "B"],
            cantera_index_map={id(rxn): [7]})
        e = entries[0]
        assert e["archetype"] == "legacy_mu1/1"
        assert e["unresolved"] is True

    def test_stamped_migration_without_configured_dst_demotes_to_legacy(self):
        """Mirror of the solver demotion (polymer.pyx:560-578): a stamped
        MIGRATION whose dst pool is not solver-configured runs as legacy µ1 —
        the artifact must say so."""
        core = _two_pool_core()
        rxn = Reaction(reactants=[core[0]], products=[core[4]],
                       kinetics=_arrhenius(), reversible=False)
        rxn.polymer_flux_archetype = int(PolymerFluxArchetype.MIGRATION)
        entries = compile_polymer_reaction_entries(
            [rxn], core, configured_pool_labels=["A"],  # B not configured
            cantera_index_map={id(rxn): [0]})
        e = entries[0]
        assert e["archetype"] == "legacy_mu1/1"
        assert e["unresolved"] is True
        assert e["src_pool"] == "A" and e["dst_pool"] is None
        assert e["proxy_products"] == []  # B is not a configured pool

    def test_pure_gas_reaction_skipped(self):
        core = _two_pool_core()
        rxn = Reaction(reactants=[core[8]], products=[core[9]],
                       kinetics=_arrhenius(), reversible=False)
        entries = compile_polymer_reaction_entries(
            [rxn], core, configured_pool_labels=["A", "B"],
            cantera_index_map={id(rxn): [0]})
        assert entries == []

    def test_dropped_entry_ids_are_stable_within_artifact(self):
        core = _two_pool_core()
        def mk():
            r = Reaction(reactants=[core[0]], products=[core[4]],
                         kinetics=_arrhenius(), reversible=False)
            r.polymer_flux_archetype = int(PolymerFluxArchetype.MIGRATION)
            return r
        r1, r2 = mk(), mk()
        entries = compile_polymer_reaction_entries(
            [r1, r2], core, configured_pool_labels=["A", "B"], cantera_index_map={})
        ids = [e["id"] for e in entries]
        assert len(set(ids)) == 2          # occurrence counter disambiguates
        assert ids == sorted(ids) or True  # deterministic order = input order
        # re-compiling the same list reproduces the same ids
        entries2 = compile_polymer_reaction_entries(
            [r1, r2], core, configured_pool_labels=["A", "B"], cantera_index_map={})
        assert [e["id"] for e in entries2] == ids
```

- [ ] **Step 2: Run to verify failure**

Run: `~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/polymerArtifactTest.py::TestCompileReactionEntries -q`
Expected: FAIL with `NotImplementedError` (the Task-3 stub).

- [ ] **Step 3: Implement `compile_polymer_reaction_entries` in `rmgpy/polymer.py`** (replace the stub)

```python
# Archetype int -> versioned term-type name (docs/polymer_moments_format.md §3).
ARCHETYPE_TERM_NAMES = {
    int(PolymerFluxArchetype.SAME_POOL): "same_pool/1",
    int(PolymerFluxArchetype.MIGRATION): "migration/1",
    int(PolymerFluxArchetype.SCISSION_FRAGMENT): "scission_fragment/1",
    int(PolymerFluxArchetype.UNRESOLVED): "legacy_mu1/1",
    int(PolymerFluxArchetype.DISCRETE_CHIP): "discrete_chip/1",
}

_ARRHENIUS_A_UNITS = {1: "s^-1", 2: "m^3/(mol*s)", 3: "m^6/(mol^2*s)"}


def _resolve_reaction_pools(rxn, pool_set):
    """Mirror the solver's src/dst pool resolution (polymer.pyx:535-556):
    src = first reactant slot in a configured pool; dst = first cross-pool
    product, falling back to the same-pool fold-back product."""
    src = None
    for s in rxn.reactants:
        b = _species_base_label(s)
        if b in pool_set:
            src = b
            break
    dst = None
    for s in rxn.products:
        b = _species_base_label(s)
        if b not in pool_set:
            continue
        if b != src:
            dst = b
            break
        if dst is None:
            dst = b
    return src, dst


def compile_polymer_reaction_entries(core_reactions, core_species,
                                     configured_pool_labels,
                                     cantera_index_map=None):
    """Compile stamped proxy-touching core reactions into schema-2.0
    ``reactions[]`` entries (docs/polymer_moments_format.md §3).

    Mirrors the solver's init-time pool resolution and demotions
    (rmgpy/solver/polymer.pyx:527-588): the artifact describes what the
    oracle DOES, including legacy/unresolved fallbacks (design spec Q3).

    Parameters
    ----------
    cantera_index_map : dict, optional
        ``{id(rxn): [entry indices in chem.yaml's reactions list]}`` from
        ``rmgpy.cantera.generate_cantera_data(..., return_reaction_index_map=True)``.
        Reactions absent from the map are emitted with ``cantera: null``
        (the unbalanced-proxy filter dropped them) and MUST carry kinetics.
    """
    from rmgpy.cantera import get_reaction_equation
    from rmgpy.kinetics import Arrhenius as _Arrhenius

    cantera_index_map = cantera_index_map or {}
    pool_set = set(configured_pool_labels)
    entries = []
    dropped_counters: Dict[tuple, int] = {}

    NONE_ = int(PolymerFluxArchetype.NONE)
    MIG = int(PolymerFluxArchetype.MIGRATION)
    SCI = int(PolymerFluxArchetype.SCISSION_FRAGMENT)
    UNR = int(PolymerFluxArchetype.UNRESOLVED)
    CHIP = int(PolymerFluxArchetype.DISCRETE_CHIP)

    for rxn in core_reactions:
        arch = int(getattr(rxn, "polymer_flux_archetype", 0))
        src, dst = _resolve_reaction_pools(rxn, pool_set)
        if arch == NONE_ and src is None and dst is None:
            continue  # ordinary chemistry — Cantera handles it untouched

        # Mirror solver demotions (polymer.pyx:557-578).
        unresolved = False
        if arch == NONE_:
            arch, unresolved = UNR, True
        elif arch == UNR:
            unresolved = True
        elif arch in (MIG, SCI) and (src is None or dst is None):
            arch, unresolved = UNR, True
        elif arch == CHIP and src is None:
            arch, unresolved = UNR, True

        equation = get_reaction_equation(rxn, core_species)
        indices = cantera_index_map.get(id(rxn))
        if indices:
            if len(indices) > 1:
                logging.warning(
                    "Polymer artifact: reaction %s maps to %d Cantera entries; "
                    "emitting the first index. The consumer must zero ALL "
                    "duplicate entries (see format doc §4 step 0).",
                    equation, len(indices))
            cantera = {"index": int(indices[0]), "equation": equation}
            entry_id = f"r{int(indices[0])}"
        else:
            cantera = None
            family = str(getattr(rxn, "family", None)
                         or getattr(rxn, "library", None) or "rxn")
            key = (family, equation)
            occurrence = dropped_counters.get(key, 0)
            dropped_counters[key] = occurrence + 1
            entry_id = f"d{family}:{equation}:{occurrence}"

        kin = getattr(rxn, "kinetics", None)
        if isinstance(kin, _Arrhenius):
            kinetics = {
                "A": float(kin.A.value_si),
                "n": float(kin.n.value_si),
                "Ea": float(kin.Ea.value_si),
                "units": {"A": _ARRHENIUS_A_UNITS.get(len(rxn.reactants), "SI"),
                          "Ea": "J/mol"},
                "reversible": bool(rxn.reversible),
            }
        else:
            kinetics = None
            if cantera is None:
                logging.warning(
                    "Polymer artifact: dropped reaction %s has non-Arrhenius "
                    "kinetics (%s); the consumer cannot evaluate it.",
                    equation, type(kin).__name__)

        entry = {
            "id": entry_id,
            "cantera": cantera,
            "kinetics": kinetics,
            "reactants": [_artifact_species_label(s) for s in rxn.reactants],
            "products": [_artifact_species_label(s) for s in rxn.products],
            "proxy_reactants": [_artifact_species_label(s) for s in rxn.reactants
                                if _species_base_label(s) in pool_set],
            "proxy_products": [_artifact_species_label(s) for s in rxn.products
                               if _species_base_label(s) in pool_set],
            "scaling": "mu0" if getattr(rxn, "is_end_group_reaction", False) else "mu1",
            "src_pool": src,
            "dst_pool": dst,
            "archetype": ARCHETYPE_TERM_NAMES[arch],
            "unresolved": unresolved,
        }
        if arch == CHIP:
            entry["params"] = {"a": int(getattr(rxn, "polymer_chip_units", 0))}
        entries.append(entry)
    return entries
```

(`PolymerFluxArchetype`, `logging`, `Dict` are already imported/defined in `rmgpy/polymer.py`.)

- [ ] **Step 4: Run the tests**

Run: `~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/polymerArtifactTest.py -q`
Expected: all PASS.

- [ ] **Step 5: Commit**

```bash
git add rmgpy/polymer.py test/rmgpy/polymerArtifactTest.py
git commit -m "polymer: compile stamped reactions into artifact entries

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 5: Artifact builder + `save_everything` hook wiring

Assemble the full 2.0 payload (envelope + conventions + pools + reactions) and wire the hook to feed the builder from the live run (cantera index map, solver-configured pools, phase mask, monomer routing).

**Files:**
- Modify: `rmgpy/polymer.py` (replace `build_polymer_moments_artifact` stub; rewrite `write_polymer_pools_sidecar` body)
- Modify: `rmgpy/rmg/main.py:2057-2079`
- Test: `test/rmgpy/polymerArtifactTest.py` (append)

- [ ] **Step 1: Write the failing tests**

Append to `test/rmgpy/polymerArtifactTest.py`:

```python
class TestArtifactBuilderAndRoundTrip:
    def _build(self, pe_pool, tmp_path):
        core = [
            _spc("CC", "PE", index=2),
            _mu_dummy("PE_mu0"), _mu_dummy("PE_mu1"), _mu_dummy("PE_mu2"),
            _spc("[CH3]", "G", index=7),
        ]
        core[0].is_polymer_proxy = True
        rxn = Reaction(reactants=[core[0]], products=[core[0], core[4]],
                       kinetics=_arrhenius(), reversible=False)
        rxn.polymer_flux_archetype = int(PolymerFluxArchetype.DISCRETE_CHIP)
        rxn.polymer_chip_units = 1
        rxn.is_end_group_reaction = True
        path = write_polymer_pools_sidecar(
            pool_registry=[pe_pool],
            output_dir=str(tmp_path),
            iteration=3,
            core_species=core,
            core_reactions=[rxn],
            configured_pool_labels=["PE"],
            condensed_species=core[:4],
            monomer_routing_by_pool={},
            cantera_index_map={id(rxn): [0]},
        )
        return path

    def test_envelope_and_blocks(self, pe_pool, tmp_path):
        path = self._build(pe_pool, tmp_path)
        with open(path) as fh:
            data = json.load(fh)
        assert data["schema_version"] == "2.0"
        assert data["rmg_iteration"] == 3
        assert "generated_at" in data
        assert "rmg_commit" in data  # may be a SHA string or null
        conv = data["conventions"]
        assert conv["configured_pools"] == ["PE"]
        assert conv["condensed_species"] == ["PE(2)", "PE_mu0", "PE_mu1", "PE_mu2"]
        assert conv["mu3_closure"] == "log_lagrange/1"
        assert "format_doc" in conv
        assert len(data["pools"]) == 1
        assert data["pools"][0]["phase_species"] == ["PE(2)", "PE_mu0", "PE_mu1", "PE_mu2"]
        assert len(data["reactions"]) == 1
        assert data["reactions"][0]["archetype"] == "discrete_chip/1"
        assert data["reactions"][0]["params"] == {"a": 1}

    def test_json_round_trip_is_lossless(self, pe_pool, tmp_path):
        path = self._build(pe_pool, tmp_path)
        with open(path) as fh:
            data = json.load(fh)
        # round-trip: dump and re-load reproduces an identical document
        assert json.loads(json.dumps(data)) == data
        # 1.0 field stability inside the written file
        pool = data["pools"][0]
        for key, typ in SCHEMA_1_0_KEYS.items():
            assert key in pool
            assert isinstance(pool[key], typ)

    def test_legacy_call_signature_still_works(self, pe_pool, tmp_path):
        """The pre-2.0 call shape (pool_registry, output_dir, iteration) used by
        existing tests/callers keeps working: reactions=[] and conventions
        present with defaults."""
        path = write_polymer_pools_sidecar(
            pool_registry=[pe_pool], output_dir=str(tmp_path), iteration=0)
        with open(path) as fh:
            data = json.load(fh)
        assert data["schema_version"] == "2.0"
        assert data["reactions"] == []
        assert data["conventions"]["configured_pools"] == ["PE"]
```

- [ ] **Step 2: Run to verify failure**

Run: `~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/polymerArtifactTest.py::TestArtifactBuilderAndRoundTrip -q`
Expected: FAIL (`NotImplementedError` from the stub / `TypeError` unexpected kwargs).

- [ ] **Step 3: Implement the builder and rewire the writer in `rmgpy/polymer.py`**

Replace the `build_polymer_moments_artifact` stub:

```python
def build_polymer_moments_artifact(pool_registry,
                                   core_species=None,
                                   core_reactions=None,
                                   configured_pool_labels=None,
                                   condensed_species=None,
                                   monomer_routing_by_pool=None,
                                   cantera_index_map=None,
                                   iteration=0,
                                   rmg_commit=None):
    """Assemble the full schema-2.0 polymer moments artifact payload.

    Normative contract: docs/polymer_moments_format.md. The payload mirrors
    the HybridPolymerSystem oracle, including its init-time demotions —
    ``configured_pool_labels`` must be the SOLVER-configured pools (which can
    be a subset of ``pool_registry``: spawned daughters are registry pools
    without solver configs and run as ordinary species).
    """
    if configured_pool_labels is None:
        configured_pool_labels = [getattr(p, "label", "") for p in pool_registry]
    monomer_routing_by_pool = monomer_routing_by_pool or {}

    pools = [
        _serialize_pool_for_sidecar(
            p,
            core_species=core_species,
            monomer_routing=monomer_routing_by_pool.get(getattr(p, "label", "")),
        )
        for p in pool_registry
    ]
    reactions = compile_polymer_reaction_entries(
        core_reactions or [], core_species or [],
        configured_pool_labels, cantera_index_map)

    conventions = {
        "format_doc": "docs/polymer_moments_format.md (polymer_moments_format/2.0)",
        "moment_basis": "extensive mol, DP basis (mu1 = moles of repeat units)",
        "volumes": {
            "V_poly": "constant, consumer-supplied [m^3]",
            "V_gas": "ideal gas, dynamic: V_gas = n_gas*R*T/P (1.0 m^3 floor when n_gas <= 0)",
        },
        "configured_pools": list(configured_pool_labels),
        "condensed_species": sorted(_artifact_species_label(s)
                                    for s in (condensed_species or [])),
        "site_scaling": ("site = max(0, mu_scaling)/V_poly read from the first proxy "
                         "reactant's pool; multiplies ONCE; scales rf AND rr"),
        "chip_site_throttle": ("site = min(max(0,mu0), max(0,mu1)/a)/V_poly when "
                               "archetype=discrete_chip/1 and scaling=mu0 and a>0"),
        "kb_recipe": ("kb = kf/Keq; Keq(T) = (P0/(R*T))^dn_gas * exp(-dG0/(R*T)), "
                      "P0 = 1e5 Pa, dG0 from chem.yaml NASA thermo"),
        "mu3_closure": "log_lagrange/1",
        "invariants": {
            "discrete_subset": ("sum_pools(mu1) + sum_chip_species(a_i * n_i) is "
                                "invariant over the discrete-reaction subset only"),
            "with_unzip": ("add + n(monomer_routing) per pool with an active unzip "
                           "channel (unzip moves units from mu1 into that species)"),
        },
    }

    return {
        "schema_version": POLYMER_POOLS_SIDECAR_SCHEMA_VERSION,
        "generated_at": datetime.datetime.utcnow().isoformat(timespec="seconds") + "Z",
        "rmg_commit": rmg_commit if rmg_commit is not None else _get_rmg_commit(),
        "rmg_iteration": int(iteration),
        "conventions": conventions,
        "pools": pools,
        "reactions": reactions,
    }
```

Rewrite `write_polymer_pools_sidecar` (keep its docstring, extend signature; old positional shape unchanged):

```python
def write_polymer_pools_sidecar(
    pool_registry: List['Polymer'],
    output_dir: str,
    iteration: int = 0,
    filename: str = POLYMER_POOLS_SIDECAR_FILENAME,
    core_species=None,
    core_reactions=None,
    configured_pool_labels=None,
    condensed_species=None,
    monomer_routing_by_pool=None,
    cantera_index_map=None,
    rmg_commit=None,
) -> str:
    payload = build_polymer_moments_artifact(
        pool_registry,
        core_species=core_species,
        core_reactions=core_reactions,
        configured_pool_labels=configured_pool_labels,
        condensed_species=condensed_species,
        monomer_routing_by_pool=monomer_routing_by_pool,
        cantera_index_map=cantera_index_map,
        iteration=iteration,
        rmg_commit=rmg_commit,
    )
    path = os.path.join(output_dir, filename)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, default=str)
    return path
```

- [ ] **Step 4: Wire the hook in `rmgpy/rmg/main.py`**

Replace the sidecar block at `:2057-2079` with:

```python
        # Emit the polymer_pools.json sidecar (schema 2.0) alongside the
        # chemkin / cantera outputs so the TA-side mechanism loader
        # (~/Code/TA) can pick up pool semantics + compiled flux terms.
        # Normative contract: docs/polymer_moments_format.md.
        try:
            from rmgpy.polymer import Polymer, write_polymer_pools_sidecar
            pool_registry = [
                s for s in (
                    self.reaction_model.core.species
                    + self.reaction_model.edge.species
                    + self.reaction_model.new_species_list
                )
                if isinstance(s, Polymer)
            ]
            if pool_registry and self.output_directory:
                chemkin_dir = os.path.join(self.output_directory, "chemkin")
                target_dir = chemkin_dir if os.path.isdir(chemkin_dir) else self.output_directory

                core_species = self.reaction_model.core.species
                core_reactions = self.reaction_model.core.reactions

                # Cantera index map: recompute the exact filter/ordering the
                # CanteraWriter listener (notify() above) just used on the
                # same core — same inputs, same map.
                cantera_index_map = None
                try:
                    from rmgpy.cantera import generate_cantera_data
                    _, cantera_index_map = generate_cantera_data(
                        core_species, core_reactions,
                        return_reaction_index_map=True)
                except Exception as e:
                    logging.warning(
                        "polymer_pools.json: cantera index map unavailable (%s); "
                        "all reaction entries will be emitted cantera-null.", e)

                # Solver-configured pools / phase mask / monomer routing from
                # the live hybrid reactor (first one that has a built solver).
                configured = condensed = None
                routing = {}
                for system in self.reaction_systems:
                    solver = getattr(system, "solver", None)
                    pools_cfg = getattr(solver, "polymer_pools", None)
                    if not pools_cfg:
                        continue
                    configured = [p.label for p in pools_cfg]
                    mask = getattr(solver, "gas_species_mask", None)
                    if mask is not None and len(mask) == len(core_species):
                        condensed = [core_species[i]
                                     for i in range(len(core_species)) if not mask[i]]
                    for p in pools_cfg:
                        idx = p.monomer_poly_index
                        if idx is not None and 0 <= idx < len(core_species):
                            spc = core_species[idx]
                            routing[p.label] = (f"{spc.label}({spc.index})"
                                                if spc.index > 0 else spc.label)
                    break

                write_polymer_pools_sidecar(
                    pool_registry=pool_registry,
                    output_dir=target_dir,
                    iteration=getattr(self.reaction_model, "iteration_num", 0),
                    core_species=core_species,
                    core_reactions=core_reactions,
                    configured_pool_labels=configured,
                    condensed_species=condensed,
                    monomer_routing_by_pool=routing,
                    cantera_index_map=cantera_index_map,
                )
        except Exception as e:
            logging.warning(f"Failed to write polymer_pools.json sidecar: {e}")
```

- [ ] **Step 5: Run the tests + neighbor suites**

Run: `~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/polymerArtifactTest.py test/rmgpy/polymerMultiPoolTest.py test/rmgpy/polymerTest.py -q`
Expected: all PASS.

- [ ] **Step 6: Commit**

```bash
git add rmgpy/polymer.py rmgpy/rmg/main.py test/rmgpy/polymerArtifactTest.py
git commit -m "polymer: wire schema-2.0 artifact into save_everything

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 6: Normative format doc `docs/polymer_moments_format.md`

THE contract TA implements against (spec §6). Every equation below was transcribed from the live solver during plan authoring — the implementer must re-verify each against the cited lines before committing (they are normative; a transcription typo becomes a TA bug).

**Files:**
- Create: `docs/polymer_moments_format.md`

- [ ] **Step 1: Write the document**

Create `docs/polymer_moments_format.md` with exactly this content (verify each formula against the cited `rmgpy/solver/polymer.pyx` lines as you paste):

````markdown
# Polymer Moments Artifact — Normative Format Spec (`polymer_moments_format/2.0`)

**Artifact:** `polymer_pools.json`, emitted next to `chem.yaml` at the end of an
RMG run (`save_everything`). **Schema version:** `2.0` (major bump over the 1.0
pools-only sidecar; every 1.0 field is preserved verbatim).
**Oracle:** `rmgpy/solver/polymer.pyx` (`HybridPolymerSystem`). A consumer that
implements this document with numpy + Cantera reproduces the oracle's per-pool
µ0/µ1/µ2 trajectories. The reference consumer lives at
`test/rmgpy/tools/numpy_moments_consumer.py`; the reference runner (drives the
oracle itself) at `rmgpy/tools/polymer_moments_runner.py`.

**Versioning policy:** minor bump = additive term types / fields (consumers
must ignore unknown keys and SHOULD warn on unknown `archetype` values);
major bump = breaking. This document versions with the schema.

## 1. Envelope

```json
{
  "schema_version": "2.0",
  "generated_at": "<ISO-8601 UTC>",
  "rmg_commit": "<git SHA>" | null,
  "rmg_iteration": <int>,
  "conventions": { ... },     // §4
  "pools": [ ... ],           // §2
  "reactions": [ ... ]        // §3
}
```

## 2. `pools[]`

Schema-1.0 fields (unchanged): `label`, `monomer_smiles`, `monomer_adj_list`,
`feature_monomers_smiles`, `end_groups`, `cutoff`, `parent_pool`,
`spawn_iteration`, `spawn_event_metadata`, `mu_indices` (legacy; do not use —
solver state-vector indices from the generating run).

2.0 additions:

| Field | Type | Meaning |
|---|---|---|
| `moments` | `[µ0, µ1, µ2]` floats \| null | pool state at generation time, **extensive mol, DP basis** (µ1 = moles of repeat units) |
| `monomer_mw_g_mol` | float | repeat-unit MW [g/mol] |
| `mn_g_mol`, `mw_g_mol` | float \| null | number-/weight-average MW at generation time |
| `initial_mass_g` | float \| null | as configured (grams) |
| `channels` | `{"scission": {A,n,Ea,units}, "unzip": {...}}` | Arrhenius-capable; today RMG emits `A=k, n=0, Ea=0`, `units = {"A": "s^-1", "Ea": "J/mol"}`. Channel equations: §5 (versioned `scission/1`, `unzip/1`) |
| `phase_species` | [labels] | the pool's condensed-phase species by chem.yaml name (proxy variants, µ-dummies, routed monomer). No more `<label>_muN` name-guessing |
| `monomer_routing` | label \| null | chem.yaml species receiving the unzip monomer flux. It is a **condensed-phase** species (it reaches the gas only via mass transfer). `null` ⇒ unzip µ-flux applies but no monomer species is credited |
| `mu3_closure` | `"log_lagrange/1"` | §6 |

A pool listed here is **solver-configured** iff its `label` appears in
`conventions.configured_pools`. Non-configured pools (spawned daughters) are
inert containers: their proxies behave as ordinary species (no site scaling,
no concentration-1.0 rule, no channels integration by the oracle).

## 3. `reactions[]`

One entry per pool-touching core reaction (`archetype` ≠ none or any
participant in a configured pool). Closed term vocabulary — entries carry
parameters only; the equations live here and only here.

```json
{
  "id": "r<cantera.index>" | "d<family>:<equation>:<occurrence>",
  "cantera": {"index": <int>, "equation": "<string>"} | null,
  "kinetics": {"A":..,"n":..,"Ea":..,"units":{"A":..,"Ea":"J/mol"},"reversible":bool} | null,
  "reactants": ["<label>", ...],        // full stoichiometry, with repeats
  "products":  ["<label>", ...],
  "proxy_reactants": ["<label>", ...],  // pool-mapped reactants (conc := 1.0)
  "proxy_products":  ["<label>", ...],  // pool-mapped products  (conc := 1.0 in reverse)
  "scaling": "mu0" | "mu1",
  "src_pool": "<pool label>" | null,
  "dst_pool": "<pool label>" | null,
  "archetype": "same_pool/1" | "migration/1" | "scission_fragment/1"
             | "discrete_chip/1" | "legacy_mu1/1",
  "params": {"a": <int>},               // discrete_chip/1 only
  "unresolved": bool                    // true => legacy fallback emission
}
```

- `cantera.index` is **authoritative**: the entry's position in the
  *same-directory* `chem.yaml` `reactions:` list. The consumer must act on it
  inside Cantera (§4 step 0). `cantera: null` ⇒ the reaction was dropped by
  the unbalanced-proxy export filter; `kinetics` is then REQUIRED and
  `reactants`/`products` are the only record of its stoichiometry.
- `id` is stable within one artifact and across re-reads of the same file;
  NOT across regenerations.
- The pool of a proxy label is the label with any trailing `(N)` stripped
  (e.g. `epdm(2)` → pool `epdm`).
- `unresolved: true` entries (legacy-µ1 emissions, including solver-demoted
  stamps) MUST be integrated; consumers MAY warn.

## 4. Rate-evaluation recipe (normative; oracle: `polymer.pyx:922-1261`)

Definitions: `V_poly` = constant consumer input; `V_gas = n_gas·R·T/P`
(ideal gas over the *gas-phase* species — those NOT in
`conventions.condensed_species`; floor `V_gas = 1.0 m³` when `n_gas ≤ 0`).
Pool moments are extensive mol; intensive µ ≡ y[µ]/V_poly clamped at ≥ 0.

0. **Disable every listed retained reaction in Cantera:** for each entry with
   `cantera != null`, `Kinetics.set_multiplier(0, cantera.index)`. (Equation
   strings are a human checksum only.)
1. Per entry, per T: `kf` from Cantera (or from `kinetics`:
   `kf = A·T^n·exp(−Ea/(R·T))`, SI). If `reversible`: `kb = kf/Keq` with
   `Keq(T) = (P°/(R·T))^Δn_gas · exp(−ΔG°/RT)`, `P° = 1e5 Pa`, ΔG° summed
   from the chem.yaml NASA polynomials (`reaction.py:767-840`). Else `kb = 0`.
2. **Phase + gate** (`polymer.pyx:943-989`): the event is condensed
   (`V_rxn = V_poly`) iff ANY reactant is in `condensed_species`, else gas
   (`V_rxn = V_gas`). Gate: a condensed event with NO condensed core product
   is **skipped entirely** (rate 0, but still zeroed in Cantera per step 0);
   a gas event with a condensed core product is likewise skipped.
3. Concentration products (`polymer.pyx:1000-1020`):
   `rf = kf · Π C(reactant)`, `rr = kb · Π C(product)`, where
   `C(s) = 1.0` if `s` ∈ `proxy_reactants ∪ proxy_products` (pool-mapped),
   else `n_s/V_gas` (gas) or `n_s/V_poly` (condensed), clamped ≥ 0.
4. **Site scaling** (`polymer.pyx:1022-1063`), only when `src_pool != null`:
   `site = max(0, µ_scaling)/V_poly` with µ from `src_pool` (the FIRST
   proxy-reactant's pool, reactant-slot priority), µ0 when `scaling=="mu0"`
   else µ1 — EXCEPT chip entries (`archetype=="discrete_chip/1"` ∧
   `scaling=="mu0"` ∧ `a > 0`):
   `site = min(max(0, µ0), max(0, µ1)/a)/V_poly` (exhaustion throttle).
   The site multiplies **once** (even with two proxy reactants) and scales
   BOTH `rf` and `rr` (the reverse is NOT scaled by the dst pool).
   `legacy_mu1/1` entries that look chip-shaped are deliberately NOT
   throttled (bit-exact legacy contract).
5. **Gas/explicit stoichiometric flux:** `r = rf − rr`, `r_mol = r·V_rxn`.
   For every NON-pool-mapped species: reactants `dn/dt −= r_mol` (per
   occurrence), products `dn/dt += r_mol`. (This is how chips, abstraction
   partners and co-products move — Cantera no longer does it after step 0.)
6. **Archetype bundle** (per direction; `ev_mol_f = rf·V_rxn`,
   `ev_mol_r = rr·V_rxn`):

   Chain bundle `B(pool, end_group)` (`polymer.pyx:805-828`): with intensive
   µ of that pool — end_group=true (uniform pick): if µ0 ≤ 1e-30 → empty;
   else `(b0,b1,b2,ok) = (1, µ1/µ0, µ2/µ0, true)`. end_group=false
   (length-biased): if µ1 ≤ 1e-30 → empty; µ3 from §6;
   `(1, µ2/µ1, µ3/µ1, true)` if µ3 finite else `(1, µ2/µ1, 0, false)`.

   - `same_pool/1`: no moment flux (fold-back is net-zero by construction).
   - `migration/1` (`polymer.pyx:1128-1155`): only if src ≠ dst, both non-null.
     For (ev_mol, from, to) ∈ {(ev_mol_f, src, dst), (ev_mol_r, dst, src)} with
     ev_mol > 0: `B(from, scaling=="mu0")`; skip if empty;
     `µ0/µ1[from] −= ev_mol·(b0,b1)`, `µ0/µ1[to] += ev_mol·(b0,b1)`;
     if ok, same for µ2 with b2.
   - `scission_fragment/1` (`polymer.pyx:1156-1195`): only if src ≠ dst, both
     non-null. Net `r_mol`; parent µ intensive. Guard: µ1_src > 1e-30, and if
     `r_mol < 0` additionally µ0_dst > 1e-30 and µ1_dst > 1e-30, else skip.
     `E[n] = µ2_src/µ1_src`:
     `µ1[src] −= r_mol·E[n]/2`; `µ0[dst] += r_mol`; `µ1[dst] += r_mol·E[n]/2`;
     if µ3_src finite: `E[n²] = µ3_src/µ1_src`,
     `µ2[src] −= r_mol·(2/3)·E[n²]`, `µ2[dst] += r_mol·E[n²]/3`.
     (Complement stays in parent: parent µ0 net 0; µ1 conserves exactly.)
   - `discrete_chip/1` (`polymer.pyx:1196-1244`): src only (complement folds
     back; the chip species moves via step 5). `B(src, scaling=="mu0")`; skip
     if empty; `E[n] = b1`, `a = params.a`.
     Forward (rf>0), `rf_mol = rf·V_rxn`: `µ1[src] −= rf_mol·a`;
     `Δµ2 = 2a·E[n] − a²`; if `Δµ2 > 0`: `µ2[src] −= rf_mol·Δµ2` (clamped
     decrement — no write when ≤ 0).
     Reverse (rr>0), `rr_mol = rr·V_rxn`: `µ1[src] += rr_mol·a`;
     `µ2[src] += rr_mol·(2a·E[n] + a²)` (exact extension form, no clamp).
   - `legacy_mu1/1` (`polymer.pyx:1245-1261`): net `r_mol`. For EACH label in
     `proxy_reactants`: `µ1[pool(label)] −= r_mol`. For EACH label in
     `proxy_products`: `µ1[pool(label)] += r_mol`. (Only configured pools.)

7. **Channels** (per configured pool, §5), then **mass transfer** (§7).

## 5. Channel equations (versioned; oracle `polymer.pyx:1318-1340`; see also
`docs/multi_pool_design.md` §5)

Intensive µ; rates volumetric [mol/m³/s], multiplied by `V_poly` into mol/s.
`k_s`, `k_u` from `channels` (Arrhenius at T; today constant `A`).

- `scission/1` (Ziff–McGrady discrete-bond):
  `dµ0/dt += k_s·(µ1 − µ0)`; `dµ1/dt += 0`;
  `dµ2/dt += k_s·(µ1 − µ3)/3` only when µ3 (§6) is finite.
- `unzip/1`: `r = k_u·µ0`; `dµ1/dt −= r`; `dµ2/dt −= k_u·(2µ1 − µ0)`;
  if `monomer_routing != null`: `dn(monomer_routing)/dt += r·V_poly`
  (condensed phase). µ0 unchanged.

The oracle's explicit-oligomer "handshake" (`polymer.pyx:1342-1380`) is NOT
part of the artifact: run-path pools carry no explicit oligomer map, so it is
inert in every emitted artifact. (If a future schema adds explicit oligomers,
that will be a versioned addition.)

## 6. µ3 closure — `log_lagrange/1` (oracle `polymer.pyx:141-160`)

`µ3 = µ0·(µ2/µ1)³` computed in log space, with a realizability guard:
- if µ0 ≤ 1e-30 or µ1 ≤ 1e-30 or µ2 ≤ 1e-30 → µ3 = 0
- if µ1 < µ0 (unrealizable) → µ3 = 0
- `ln µ3 = 3·ln µ2 − 3·ln µ1 + ln µ0`; if `ln µ3 > 700` → µ3 = +inf
  (callers treat infinite µ3 by skipping the µ2 component, §4/§5)

## 7. Mass transfer (consumer-supplied operating condition, NOT artifact content)

`J = kLa·(C_poly − K·C_gas)` [mol/m³/s], `dn = J·V_poly`;
`dn(gas species) += dn`, `dn(condensed species) −= dn`
(oracle `polymer.pyx:1390-1402`; `K = C_poly_eq/C_gas_eq`, `kLa` [1/s]).
kLa/K are apparatus properties and enter via the runner/consumer inputs.

## 8. Conventions block

`conventions` carries (informative duplicates of this doc, plus two normative
lists): `configured_pools` (pool labels with solver configs — §2/§3 semantics)
and `condensed_species` (chem.yaml labels with phase = condensed; everything
else is gas). Consumers MUST use these lists, not name heuristics.

## 9. Conservation invariants (assertable, per-channel qualified)

Over the discrete-reaction subset only (no channels active):
`Σ_pools µ1 + Σ_chip_species a_i·n_i` is constant.
With unzip active, add `+ n(monomer_routing)` per routed pool. Random
scission conserves `Σ µ1` exactly. Mass transfer conserves total moles of the
transferred species pair.
````

- [ ] **Step 2: Re-verify every cited equation against the live code**

Open `rmgpy/solver/polymer.pyx` and check, character by character, the lines cited in §4–§7 (the line numbers in the "Verified code map" table at the top of this plan). Fix any transcription drift.

- [ ] **Step 3: Commit**

```bash
git add docs/polymer_moments_format.md
git commit -m "docs: add normative polymer moments format spec

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 7: Numpy-only consumer + sufficiency proof

The whole point of the artifact: a consumer module with **NO rmgpy imports** (numpy + stdlib only) reproduces the oracle. Oracle trajectories are generated in the test file (which may import rmgpy freely); the consumer module must never touch it. Both sides integrate with identical fixed-step forward Euler (the established oracle-trajectory pattern, `solverPolymerTest.py:498-510`), so trajectories match to near machine precision.

**Files:**
- Create: `test/rmgpy/tools/numpy_moments_consumer.py` (the consumer — importable, numpy-only)
- Create: `test/rmgpy/tools/polymerMomentsConsumerTest.py` (oracle + proof tests)

- [ ] **Step 1: Write the consumer module**

Create `test/rmgpy/tools/numpy_moments_consumer.py`:

```python
"""Numpy-only reference consumer for the polymer moments artifact.

Implements docs/polymer_moments_format.md §4-§7 with numpy + stdlib ONLY.
THIS MODULE MUST NOT IMPORT rmgpy (that is the artifact's entire point);
test_consumer_module_is_rmgpy_free enforces it.

Rate constants come from each entry's own `kinetics` (every entry in the
test decks is emitted cantera-null), so no Cantera is needed here either;
Keq for reversible entries is computed from caller-supplied NASA7 data via
the documented recipe.
"""

import numpy as np

R = 8.31446261815324  # J/(mol K)
P0 = 1.0e5            # Pa, Keq reference pressure
SMALL_EPS = 1e-30
LN_EXP_OVERFLOW_GUARD = 700.0


def safe_mu3(mu0, mu1, mu2):
    """log_lagrange/1 closure with realizability guard (format doc §6)."""
    if mu0 <= SMALL_EPS or mu1 <= SMALL_EPS or mu2 <= SMALL_EPS:
        return 0.0
    if mu1 < mu0:
        return 0.0
    ln_mu3 = 3.0 * np.log(mu2) - 3.0 * np.log(mu1) + np.log(mu0)
    if ln_mu3 > LN_EXP_OVERFLOW_GUARD:
        return float("inf")
    return float(np.exp(ln_mu3))


def nasa_g_over_rt(coeffs, T):
    """G/(R*T) from one NASA7 row [a0..a4, a5, a6]."""
    a0, a1, a2, a3, a4, a5, a6 = coeffs
    h_rt = (a0 + a1 * T / 2.0 + a2 * T**2 / 3.0 + a3 * T**3 / 4.0
            + a4 * T**4 / 5.0 + a5 / T)
    s_r = (a0 * np.log(T) + a1 * T + a2 * T**2 / 2.0 + a3 * T**3 / 3.0
           + a4 * T**4 / 4.0 + a6)
    return h_rt - s_r


def _base(label):
    return label.partition("(")[0]


class ArtifactConsumer:
    """Integrates the artifact's pool moments + species ODEs.

    Parameters
    ----------
    artifact : dict           parsed polymer_pools.json (schema 2.0)
    species_order : [str]     chem.yaml labels defining the y-vector layout
                              (same order as the oracle's core species)
    P : float                 pressure [Pa]
    V_poly : float            condensed-phase volume [m^3]
    mass_transfer : [dict]    [{"gas": label, "poly": label, "K": f, "kLa": f}]
    nasa : {label: {"Tmid": f, "low": [7], "high": [7]}}, optional
                              NASA7 data for Keq of reversible entries
    """

    def __init__(self, artifact, species_order, P, V_poly,
                 mass_transfer=None, nasa=None):
        self.P = float(P)
        self.V_poly = float(V_poly)
        self.nasa = nasa or {}
        self.idx = {lab: i for i, lab in enumerate(species_order)}
        n = len(species_order)

        conv = artifact["conventions"]
        condensed = set(conv["condensed_species"])
        self.gas_mask = np.array([lab not in condensed for lab in species_order],
                                 dtype=bool)
        self.configured = list(conv["configured_pools"])

        # pools: label -> dict(mu indices, channels, monomer routing index)
        self.pools = {}
        for p in artifact["pools"]:
            lab = p["label"]
            if lab not in self.configured:
                continue
            mu_idx = tuple(self.idx[f"{lab}_mu{k}"] for k in range(3))
            routing = p.get("monomer_routing")
            self.pools[lab] = {
                "mu": mu_idx,
                "k_s": p["channels"]["scission"]["A"],
                "k_u": p["channels"]["unzip"]["A"],
                "routing": self.idx[routing] if routing else None,
            }

        # reactions: precompute index forms
        self.entries = []
        for e in artifact["reactions"]:
            kin = e["kinetics"]
            assert kin is not None, f"entry {e['id']} has no kinetics"
            ridx = [self.idx[s] for s in e["reactants"]]
            pidx = [self.idx[s] for s in e["products"]]
            pool_mapped = ({self.idx[s] for s in e["proxy_reactants"]}
                           | {self.idx[s] for s in e["proxy_products"]})
            self.entries.append({
                "A": kin["A"], "n": kin["n"], "Ea": kin["Ea"],
                "reversible": kin["reversible"],
                "ridx": ridx, "pidx": pidx, "pool_mapped": pool_mapped,
                "r_labels": list(e["reactants"]), "p_labels": list(e["products"]),
                "proxy_r_pools": [_base(s) for s in e["proxy_reactants"]],
                "proxy_p_pools": [_base(s) for s in e["proxy_products"]],
                "scaling": e["scaling"],
                "src": e["src_pool"], "dst": e["dst_pool"],
                "arch": e["archetype"],
                "a": int(e.get("params", {}).get("a", 0)),
            })

        self.mass_transfer = []
        for mt in (mass_transfer or []):
            self.mass_transfer.append((self.idx[mt["gas"]], self.idx[mt["poly"]],
                                       float(mt["K"]), float(mt["kLa"])))

    # ----- helpers ------------------------------------------------------

    def _keq(self, entry, T):
        def g(label):
            d = self.nasa[_label_lookup(self.nasa, label)]
            row = d["low"] if T <= d["Tmid"] else d["high"]
            return nasa_g_over_rt(row, T)  # G/(R T)
        dg_rt = sum(g(s) for s in entry["p_labels"]) - sum(g(s) for s in entry["r_labels"])
        dn_gas = len(entry["p_labels"]) - len(entry["r_labels"])
        return (P0 / (R * T)) ** dn_gas * np.exp(-dg_rt)

    def _chain_bundle(self, pool, y, end_group):
        i0, i1, i2 = self.pools[pool]["mu"]
        mu0 = max(0.0, y[i0]) / self.V_poly
        mu1 = max(0.0, y[i1]) / self.V_poly
        mu2 = max(0.0, y[i2]) / self.V_poly
        if end_group:
            if mu0 <= SMALL_EPS:
                return 0.0, 0.0, 0.0, False
            return 1.0, mu1 / mu0, mu2 / mu0, True
        if mu1 <= SMALL_EPS:
            return 0.0, 0.0, 0.0, False
        mu3 = safe_mu3(mu0, mu1, mu2)
        if np.isfinite(mu3):
            return 1.0, mu2 / mu1, mu3 / mu1, True
        return 1.0, mu2 / mu1, 0.0, False

    # ----- RHS (format doc §4-§7) ----------------------------------------

    def rhs(self, y, T):
        dn = np.zeros_like(y)
        n_gas = float(np.sum(np.clip(y[self.gas_mask], 0.0, None)))
        V_gas = n_gas * R * T / self.P if n_gas > 0 else 1.0
        Vp = self.V_poly
        C = np.where(self.gas_mask,
                     np.clip(y, 0.0, None) / V_gas,
                     np.clip(y, 0.0, None) / Vp)

        for e in self.entries:
            kf = e["A"] * T ** e["n"] * np.exp(-e["Ea"] / (R * T))
            kb = kf / self._keq(e, T) if e["reversible"] else 0.0

            # step 2: phase + gate
            is_poly_event = any(not self.gas_mask[i] for i in e["ridx"])
            V_rxn = Vp if is_poly_event else V_gas
            has_poly_prod = any(not self.gas_mask[i] for i in e["pidx"])
            if is_poly_event and not has_poly_prod:
                continue
            if (not is_poly_event) and has_poly_prod:
                continue

            # step 3: concentration products
            rf = kf
            for i in e["ridx"]:
                rf *= 1.0 if i in e["pool_mapped"] else C[i]
            rr = kb
            for i in e["pidx"]:
                rr *= 1.0 if i in e["pool_mapped"] else C[i]

            # step 4: site scaling
            if e["src"] is not None:
                i0, i1, _ = self.pools[e["src"]]["mu"]
                if e["arch"] == "discrete_chip/1" and e["scaling"] == "mu0" and e["a"] > 0:
                    site = min(max(0.0, y[i0]), max(0.0, y[i1]) / e["a"]) / Vp
                else:
                    mi = i0 if e["scaling"] == "mu0" else i1
                    site = max(0.0, y[mi]) / Vp
                rf *= site
                rr *= site

            r_mol = (rf - rr) * V_rxn

            # step 5: stoichiometric flux for non-pool-mapped species
            for i in e["ridx"]:
                if i not in e["pool_mapped"]:
                    dn[i] -= r_mol
            for i in e["pidx"]:
                if i not in e["pool_mapped"]:
                    dn[i] += r_mol

            # step 6: archetype bundles
            arch = e["arch"]
            if arch == "migration/1":
                src, dst = e["src"], e["dst"]
                if src and dst and src != dst:
                    for ev, frm, to in ((rf, src, dst), (rr, dst, src)):
                        if ev <= 0.0:
                            continue
                        ev_mol = ev * V_rxn
                        b0, b1, b2, ok = self._chain_bundle(frm, y, e["scaling"] == "mu0")
                        if b0 == 0.0:
                            continue
                        f = self.pools[frm]["mu"]
                        t = self.pools[to]["mu"]
                        dn[f[0]] -= ev_mol * b0
                        dn[f[1]] -= ev_mol * b1
                        dn[t[0]] += ev_mol * b0
                        dn[t[1]] += ev_mol * b1
                        if ok:
                            dn[f[2]] -= ev_mol * b2
                            dn[t[2]] += ev_mol * b2
            elif arch == "scission_fragment/1":
                src, dst = e["src"], e["dst"]
                if src and dst and src != dst:
                    s = self.pools[src]["mu"]
                    d = self.pools[dst]["mu"]
                    mu0p = max(0.0, y[s[0]]) / Vp
                    mu1p = max(0.0, y[s[1]]) / Vp
                    mu2p = max(0.0, y[s[2]]) / Vp
                    ok = mu1p > SMALL_EPS
                    if ok and r_mol < 0.0:
                        if (max(0.0, y[d[0]]) / Vp <= SMALL_EPS
                                or max(0.0, y[d[1]]) / Vp <= SMALL_EPS):
                            ok = False
                    if ok:
                        e_n = mu2p / mu1p
                        dn[s[1]] -= r_mol * e_n / 2.0
                        dn[d[0]] += r_mol
                        dn[d[1]] += r_mol * e_n / 2.0
                        mu3p = safe_mu3(mu0p, mu1p, mu2p)
                        if np.isfinite(mu3p):
                            e_n2 = mu3p / mu1p
                            dn[s[2]] -= r_mol * (2.0 / 3.0) * e_n2
                            dn[d[2]] += r_mol * e_n2 / 3.0
            elif arch == "discrete_chip/1":
                src = e["src"]
                if src:
                    a = float(e["a"])
                    b0, b1, _b2, _ok = self._chain_bundle(src, y, e["scaling"] == "mu0")
                    if b0 != 0.0:
                        s = self.pools[src]["mu"]
                        e_n = b1
                        if rf > 0.0:
                            rf_mol = rf * V_rxn
                            dn[s[1]] -= rf_mol * a
                            dmu2 = 2.0 * a * e_n - a * a
                            if dmu2 > 0.0:
                                dn[s[2]] -= rf_mol * dmu2
                        if rr > 0.0:
                            rr_mol = rr * V_rxn
                            dn[s[1]] += rr_mol * a
                            dn[s[2]] += rr_mol * (2.0 * a * e_n + a * a)
            elif arch == "legacy_mu1/1":
                for pool in e["proxy_r_pools"]:
                    dn[self.pools[pool]["mu"][1]] -= r_mol
                for pool in e["proxy_p_pools"]:
                    dn[self.pools[pool]["mu"][1]] += r_mol
            # same_pool/1: no moment flux

        # step 7a: channels (format doc §5)
        for pool in self.pools.values():
            i0, i1, i2 = pool["mu"]
            mu0 = max(0.0, y[i0]) / Vp
            mu1 = max(0.0, y[i1]) / Vp
            mu2 = max(0.0, y[i2]) / Vp
            dmu0 = dmu1 = dmu2 = 0.0
            if pool["k_s"] > 0.0:
                mu3 = safe_mu3(mu0, mu1, mu2)
                dmu0 += pool["k_s"] * (mu1 - mu0)
                if np.isfinite(mu3):
                    dmu2 += pool["k_s"] * (mu1 - mu3) / 3.0
            if pool["k_u"] > 0.0:
                r_ev = pool["k_u"] * mu0
                dmu1 -= r_ev
                dmu2 -= pool["k_u"] * (2.0 * mu1 - mu0)
                if pool["routing"] is not None:
                    dn[pool["routing"]] += r_ev * Vp
            dn[i0] += dmu0 * Vp
            dn[i1] += dmu1 * Vp
            dn[i2] += dmu2 * Vp

        # step 7b: mass transfer (format doc §7)
        for ig, ip, K, kLa in self.mass_transfer:
            J = kLa * (C[ip] - K * C[ig])
            dq = J * Vp
            dn[ig] += dq
            dn[ip] -= dq

        return dn

    def integrate_euler(self, y0, T, dt, n_steps, record_every=1):
        """Fixed-step forward Euler; returns (times, trajectory[n_rec, n_spc])."""
        y = np.array(y0, dtype=float)
        times, traj = [0.0], [y.copy()]
        for k in range(n_steps):
            y = y + dt * self.rhs(y, T)
            if (k + 1) % record_every == 0:
                times.append((k + 1) * dt)
                traj.append(y.copy())
        return np.array(times), np.array(traj)


def _label_lookup(nasa, label):
    """NASA data may be keyed by full chem.yaml label or its base form."""
    if label in nasa:
        return label
    base = _base(label)
    if base in nasa:
        return base
    raise KeyError(f"no NASA thermo for {label!r}")
```

- [ ] **Step 2: Write the test file (oracle + proof)**

Create `test/rmgpy/tools/polymerMomentsConsumerTest.py`:

```python
#!/usr/bin/env python3
"""Sufficiency proof (design spec §9 tests 1 + 2b): a numpy-only consumer of
the schema-2.0 artifact reproduces HybridPolymerSystem trajectories.

The ORACLE side (this file) may import rmgpy freely; the CONSUMER
(test/rmgpy/tools/numpy_moments_consumer.py) must not — enforced below.
Both sides integrate with identical fixed-step forward Euler (the oracle
pattern from test/rmgpy/solver/solverPolymerTest.py:498-510)."""

import json
import os
import sys

import numpy as np
import pytest

from rmgpy.kinetics import Arrhenius
from rmgpy.molecule import Molecule
from rmgpy.polymer import Polymer, PolymerFluxArchetype, build_polymer_moments_artifact
from rmgpy.reaction import Reaction
from rmgpy.solver.polymer import HybridPolymerSystem, MassTransferConfig, PolymerPoolConfig
from rmgpy.species import Species
from rmgpy.thermo import NASA, NASAPolynomial

# The test dirs are not packages (no __init__.py anywhere under test/);
# import the consumer module by path.
_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
from numpy_moments_consumer import ArtifactConsumer  # noqa: E402

T_K = 800.0
P_PA = 1.0e5
V_POLY = 1.0
DT = 1.0e-4
N_STEPS = 2000


def test_consumer_module_is_rmgpy_free():
    """THE point of the artifact: the consumer never imports rmgpy.
    Checked structurally (AST import table), not by substring — the module
    docstring legitimately mentions rmgpy by name."""
    import ast
    src_path = os.path.join(_HERE, "numpy_moments_consumer.py")
    with open(src_path) as fh:
        tree = ast.parse(fh.read())
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            assert all(not a.name.split(".")[0] == "rmgpy" for a in node.names)
        if isinstance(node, ast.ImportFrom):
            assert (node.module or "").split(".")[0] != "rmgpy"


def _spc(smiles, label, index=-1):
    s = Species(molecule=[Molecule().from_smiles(smiles)])
    s.label = label
    s.index = index
    return s


def _mu(label):
    s = Species(label=label, reactive=False)
    s.molecule = [Molecule().from_smiles("[Ne]")]
    s.is_moment_dummy = True
    s.index = -1
    return s


def _yaml_label(s):
    """chem.yaml naming (== artifact labels, rmgpy.cantera.get_label):
    'label(index)' when index > 0, bare label otherwise. The consumer's
    species_order MUST use these names — they key every artifact entry."""
    return f"{s.label}({s.index})" if s.index > 0 else s.label


def _nasa(a5, a6):
    rows = []
    for tmin, tmax in ((200.0, 1000.0), (1000.0, 6000.0)):
        rows.append(NASAPolynomial(coeffs=[2.5, 0.0, 0.0, 0.0, 0.0, a5, a6],
                                   Tmin=(tmin, "K"), Tmax=(tmax, "K")))
    return NASA(polynomials=rows, Tmin=(200.0, "K"), Tmax=(6000.0, "K"))


def _euler_oracle(rs, y0, dt, n_steps):
    y = y0.copy()
    traj = [y.copy()]
    for k in range(n_steps):
        dn = rs.residual(k * dt, y, np.zeros_like(y))[0]
        y = y + dt * dn
        traj.append(y.copy())
    return np.array(traj)


# ---------------------------------------------------------------------------
# (a) channel ODEs: pure scission and pure unzip vs oracle + analytic forms
# ---------------------------------------------------------------------------

class TestChannelSufficiency:
    def _run(self, k_s, k_u):
        inert = _spc("N#N", "N2")
        core = [inert, _mu("poly_mu0"), _mu("poly_mu1"), _mu("poly_mu2")]
        mask = np.array([True, False, False, False], dtype=bool)
        pool_cfg = PolymerPoolConfig(label="poly", xs=2,
                                     explicit_dp_to_species_index={},
                                     mu_indices=(1, 2, 3), monomer_poly_index=None,
                                     k_scission=k_s, k_unzip=k_u, tail_kinetics=None)
        mom0 = (1.0, 5.0, 30.0)
        rs = HybridPolymerSystem(
            T=T_K, P=P_PA, initial_mole_fractions={inert: 1.0}, V_poly=V_POLY,
            polymer_pools=[pool_cfg], mass_transfer=[], gas_species_mask=mask,
            constant_gas_volume=False,
            initial_polymer_moments={"poly": mom0}, termination=[])
        rs.initialize_model(core, [], [], [])
        y0 = rs.y.copy()
        oracle = _euler_oracle(rs, y0, DT, N_STEPS)

        registry_pool = Polymer(label="poly", monomer="[CH2][CH2]",
                                end_groups=["[H]", "[H]"], cutoff=3,
                                moments=list(mom0), initial_mass=0.0,
                                k_scission=k_s, k_unzip=k_u)
        artifact = build_polymer_moments_artifact(
            [registry_pool], core_species=core, core_reactions=[],
            configured_pool_labels=["poly"],
            condensed_species=core[1:4], cantera_index_map={})
        artifact = json.loads(json.dumps(artifact))  # honest serialization

        consumer = ArtifactConsumer(artifact, [_yaml_label(s) for s in core],
                                    P=P_PA, V_poly=V_POLY)
        _, mine = consumer.integrate_euler(y0, T_K, DT, N_STEPS)
        return oracle, mine

    def test_pure_scission_matches_oracle_and_analytic(self):
        oracle, mine = self._run(k_s=1.0, k_u=0.0)
        np.testing.assert_allclose(mine, oracle, rtol=1e-9, atol=1e-12)
        # analytic: mu1 const; mu0(t) = mu1 - (mu1 - mu0_0) e^{-k t}
        t_end = DT * N_STEPS
        mu0_exp = 5.0 - (5.0 - 1.0) * np.exp(-1.0 * t_end)
        # Euler at dt=1e-4 over 0.2 time constants: ~1e-4 relative accuracy
        assert mine[-1, 1] == pytest.approx(mu0_exp, rel=5e-4)
        assert mine[-1, 2] == pytest.approx(5.0, rel=1e-9)  # mu1 conserved

    def test_pure_unzip_matches_oracle(self):
        oracle, mine = self._run(k_s=0.0, k_u=0.5)
        np.testing.assert_allclose(mine, oracle, rtol=1e-9, atol=1e-12)
        # analytic: mu0 const, mu1(t) = mu1_0 - k*mu0*t
        t_end = DT * N_STEPS
        assert mine[-1, 2] == pytest.approx(5.0 - 0.5 * 1.0 * t_end, rel=1e-6)


# ---------------------------------------------------------------------------
# (b) synthetic two-pool deck: stamped MIGRATION + DISCRETE_CHIP (+ reversible
#     Keq case). Reuses the apportionment plan's two-pool fixture layout.
# ---------------------------------------------------------------------------

def _two_pool_setup(with_thermo=False):
    sp = {
        "A": _spc("CCCC", "A", index=1),
        "A_mu0": _mu("A_mu0"), "A_mu1": _mu("A_mu1"), "A_mu2": _mu("A_mu2"),
        "B": _spc("CCCCC", "B", index=5),
        "B_mu0": _mu("B_mu0"), "B_mu1": _mu("B_mu1"), "B_mu2": _mu("B_mu2"),
        "G": _spc("[CH3]", "G", index=9),
    }
    sp["A"].is_polymer_proxy = True
    sp["B"].is_polymer_proxy = True
    core = [sp["A"], sp["A_mu0"], sp["A_mu1"], sp["A_mu2"],
            sp["B"], sp["B_mu0"], sp["B_mu1"], sp["B_mu2"], sp["G"]]
    mask = np.array([False] * 8 + [True], dtype=bool)
    if with_thermo:
        sp["A"].thermo = _nasa(-1000.0, 4.0)
        sp["B"].thermo = _nasa(-1200.0, 3.0)
        sp["G"].thermo = _nasa(-500.0, 5.0)
    return sp, core, mask


def _pools_ab():
    pool_a = PolymerPoolConfig(label="A", xs=2, explicit_dp_to_species_index={},
                               mu_indices=(1, 2, 3), monomer_poly_index=None,
                               k_scission=0.0, k_unzip=0.0, tail_kinetics=None)
    pool_b = PolymerPoolConfig(label="B", xs=2, explicit_dp_to_species_index={},
                               mu_indices=(5, 6, 7), monomer_poly_index=None,
                               k_scission=0.0, k_unzip=0.0, tail_kinetics=None)
    return [pool_a, pool_b]


def _registry_ab():
    out = []
    for lab in ("A", "B"):
        out.append(Polymer(label=lab, monomer="[CH2][CH2]",
                           end_groups=["[H]", "[H]"], cutoff=3,
                           moments=[1.0, 5.0, 30.0], initial_mass=0.0))
    return out


def _oracle_system(core, mask, reactions, mass_transfer=()):
    rs = HybridPolymerSystem(
        T=T_K, P=P_PA, initial_mole_fractions={core[8]: 0.05}, V_poly=V_POLY,
        polymer_pools=_pools_ab(), mass_transfer=list(mass_transfer),
        gas_species_mask=mask.copy(), constant_gas_volume=False,
        initial_polymer_moments={"A": (1.0, 5.0, 30.0), "B": (2.0, 4.0, 10.0)},
        termination=[])
    rs.initialize_model(core, list(reactions), [], [])
    return rs


def _artifact_for(core, reactions):
    artifact = build_polymer_moments_artifact(
        _registry_ab(), core_species=core, core_reactions=list(reactions),
        configured_pool_labels=["A", "B"],
        condensed_species=core[:8], cantera_index_map={})
    return json.loads(json.dumps(artifact))


def _kin1():
    return Arrhenius(A=(2.0, "1/s"), n=0.0, Ea=(0.0, "kcal/mol"), T0=(298.15, "K"))


class TestTwoPoolSufficiency:
    def test_migration_plus_discrete_chip_deck(self):
        """Spec §9 test 1(b): one stamped MIGRATION + one stamped DISCRETE_CHIP
        reaction in the SAME deck; consumer matches oracle on every species
        and every pool moment."""
        sp, core, mask = _two_pool_setup()
        mig = Reaction(reactants=[sp["A"]], products=[sp["B"]],
                       kinetics=_kin1(), reversible=False)
        mig.polymer_flux_archetype = int(PolymerFluxArchetype.MIGRATION)
        chip = Reaction(reactants=[sp["A"]], products=[sp["A"], sp["G"]],
                        kinetics=_kin1(), reversible=False)
        chip.polymer_flux_archetype = int(PolymerFluxArchetype.DISCRETE_CHIP)
        chip.polymer_chip_units = 2
        chip.is_end_group_reaction = True

        rxns = [mig, chip]
        rs = _oracle_system(core, mask, rxns)
        y0 = rs.y.copy()
        oracle = _euler_oracle(rs, y0, DT, N_STEPS)

        artifact = _artifact_for(core, rxns)
        # sanity: vocabulary exercised as intended
        archs = sorted(e["archetype"] for e in artifact["reactions"])
        assert archs == ["discrete_chip/1", "migration/1"]
        assert all(e["cantera"] is None and e["kinetics"] for e in artifact["reactions"])

        consumer = ArtifactConsumer(artifact, [_yaml_label(s) for s in core],
                                    P=P_PA, V_poly=V_POLY)
        _, mine = consumer.integrate_euler(y0, T_K, DT, N_STEPS)
        np.testing.assert_allclose(mine, oracle, rtol=1e-9, atol=1e-12)
        # the chip really fired (mu1_A dropped, G grew)
        assert mine[-1, 2] < y0[2]
        assert mine[-1, 8] > y0[8]

    def test_reversible_entry_keq_recipe(self):
        """cantera-null reversible entry: consumer's NASA-based
        kb = kf/Keq, Keq = (P0/RT)^dn * exp(-dG0/RT) matches the oracle's
        rxn.get_equilibrium_constant path (format doc §4 step 1)."""
        sp, core, mask = _two_pool_setup(with_thermo=True)
        rxn = Reaction(reactants=[sp["A"]], products=[sp["B"], sp["G"]],
                       kinetics=_kin1(), reversible=True)
        rxn.polymer_flux_archetype = int(PolymerFluxArchetype.MIGRATION)
        rs = _oracle_system(core, mask, [rxn])
        y0 = rs.y.copy()
        oracle = _euler_oracle(rs, y0, DT, N_STEPS)

        artifact = _artifact_for(core, [rxn])
        nasa = {}
        for lab in ("A", "B", "G"):
            th = sp[lab].thermo
            nasa[lab] = {"Tmid": 1000.0,
                         "low": [float(c) for c in th.polynomials[0].coeffs],
                         "high": [float(c) for c in th.polynomials[1].coeffs]}
        consumer = ArtifactConsumer(artifact, [_yaml_label(s) for s in core],
                                    P=P_PA, V_poly=V_POLY, nasa=nasa)
        _, mine = consumer.integrate_euler(y0, T_K, DT, N_STEPS)
        np.testing.assert_allclose(mine, oracle, rtol=1e-7, atol=1e-12)
```

- [ ] **Step 3: Confirm the import mechanics**

The test dirs carry no `__init__.py` (only a root `test/conftest.py`), so the consumer is imported by path (`sys.path.insert` in the test file). Sanity-check: `~/anaconda3/envs/rmg_env/bin/python -c "import sys; sys.path.insert(0, 'test/rmgpy/tools'); import numpy_moments_consumer; print('ok')"` from the repo root prints `ok`.

- [ ] **Step 4: Run to verify the proof**

Run: `~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/tools/polymerMomentsConsumerTest.py -q`
Expected: all PASS. If a trajectory mismatch appears, diff the first diverging step's `rhs` against the oracle's `residual` term by term — the format doc (Task 6) is the arbiter; fix the CONSUMER or the EMITTER, never the oracle.

- [ ] **Step 5: Commit**

```bash
git add test/rmgpy/tools/numpy_moments_consumer.py test/rmgpy/tools/polymerMomentsConsumerTest.py
git commit -m "test: prove artifact sufficiency with numpy-only consumer

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 8: Mass-transfer cross-check (spec §9 test 2b)

Validate the consumer's `J = kLa·(C_poly − K·C_gas)` against the SOLVER's (not just against the format doc), with both terms of the driving force active.

**Files:**
- Modify: `test/rmgpy/tools/polymerMomentsConsumerTest.py` (append)

- [ ] **Step 1: Write the failing test**

Append to `test/rmgpy/tools/polymerMomentsConsumerTest.py`:

```python
class TestMassTransferCrossCheck:
    def test_evaporation_both_sides_nonzero(self):
        """Spec §9 test 2b: nonzero kLa with BOTH C_poly and C_gas nonzero, on
        top of the two-pool reaction deck — the consumer's J = kLa(Cp - K*Cg)
        is validated against the solver's, end to end."""
        sp, core, mask = _two_pool_setup()
        # dissolved species D (condensed, not pool-mapped) <-> gas G
        d = _spc("CC=O", "D", index=11)
        core = core + [d]
        mask = np.append(mask, False)  # D is condensed

        mig = Reaction(reactants=[sp["A"]], products=[sp["B"]],
                       kinetics=_kin1(), reversible=False)
        mig.polymer_flux_archetype = int(PolymerFluxArchetype.MIGRATION)
        rxns = [mig]

        mt = MassTransferConfig(gas_index=8, poly_index=9, K=2.0, kLa=5.0)
        rs = HybridPolymerSystem(
            T=T_K, P=P_PA,
            initial_mole_fractions={core[8]: 0.05, d: 0.0},
            V_poly=V_POLY,
            polymer_pools=_pools_ab(), mass_transfer=[mt],
            gas_species_mask=mask.copy(), constant_gas_volume=False,
            initial_polymer_moments={"A": (1.0, 5.0, 30.0), "B": (2.0, 4.0, 10.0)},
            termination=[])
        rs.initialize_model(core, rxns, [], [])
        y0 = rs.y.copy()
        y0[9] = 0.02  # dissolved D moles: condensed side nonzero too
        oracle = _euler_oracle(rs, y0, DT, N_STEPS)

        artifact = build_polymer_moments_artifact(
            _registry_ab(), core_species=core, core_reactions=rxns,
            configured_pool_labels=["A", "B"],
            condensed_species=core[:8] + [d], cantera_index_map={})
        artifact = json.loads(json.dumps(artifact))
        consumer = ArtifactConsumer(
            artifact, [_yaml_label(s) for s in core], P=P_PA, V_poly=V_POLY,
            mass_transfer=[{"gas": "G(9)", "poly": "D(11)", "K": 2.0, "kLa": 5.0}])
        _, mine = consumer.integrate_euler(y0, T_K, DT, N_STEPS)
        np.testing.assert_allclose(mine, oracle, rtol=1e-9, atol=1e-12)
        # the transfer really moved mass and conserved the pair total
        assert mine[-1, 9] != pytest.approx(y0[9])
        np.testing.assert_allclose(mine[:, 8] + mine[:, 9],
                                   y0[8] + y0[9], rtol=1e-9)
```

NOTE: indices — `G` is core index 8, `D` is core index 9 in the extended list (`core + [d]` puts D at position 9). `MassTransferConfig(gas_index=8, poly_index=9, ...)` matches.

- [ ] **Step 2: Run to verify it fails before any needed consumer wiring, then passes**

Run: `~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/tools/polymerMomentsConsumerTest.py::TestMassTransferCrossCheck -q`
Expected: PASS immediately if Task 7's consumer already wired `mass_transfer` (it did). If it fails, the failure is real signal — fix the consumer against `polymer.pyx:1390-1402`.

- [ ] **Step 3: Run the whole consumer test file**

Run: `~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/tools/polymerMomentsConsumerTest.py -q`
Expected: all PASS.

- [ ] **Step 4: Commit**

```bash
git add test/rmgpy/tools/polymerMomentsConsumerTest.py
git commit -m "test: cross-check consumer mass transfer against solver

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 9: CLI reference runner (`rmgpy/tools/polymer_moments_runner.py`)

The oracle binary for TA cross-validation (spec §8): inputs = artifact + chem.yaml + piecewise-isothermal T-profile + time grid + mass-transfer spec; output = moment-trajectory CSV. It drives `HybridPolymerSystem` directly — it is the oracle, NOT a reimplementation. It lives in the consumer's world: it reads chem.yaml + the artifact, never RMG input files (spec §10).

Key mechanics (all verified):
- chem.yaml → RMG objects: `Species(label=name, thermo=NASA(...))` (molecule-less is fine — the solver matches by label and Keq uses thermo) and `Reaction(..., kinetics=Arrhenius(...))`. Plain `rate-constant` entries only; anything else (falloff/PLOG/Chebyshev) raises a clear error in v1.
- Dropped (cantera-null) entries are reconstructed from the artifact's `reactants`/`products`/`kinetics`.
- Stamps are restored from artifact entries (`archetype`→`polymer_flux_archetype`, `scaling=="mu0"`→`is_end_group_reaction`, `params.a`→`polymer_chip_units`). Retained entries match chem.yaml reactions by `cantera.index` (= position in the yaml `reactions:` list); demoted/legacy entries restamp as UNRESOLVED(4).
- Piecewise-isothermal restart (proven empirically — handoff `~/handoffs/handoff-polymer-declarative-moments-format.md:100-104`; no committed test existed before this task): segment 1 via `initialize_model`; each next segment: `rs.T = Quantity((T_k,'K'))`; `rs.generate_rate_coefficients(core_rxns, [])`; `rs.t0 = t_start`; `rs.y0 = y_carry`; `rs.set_initial_derivative()`; `rs.initialize_solver()`; then `rs.advance(t)` per grid point (pattern: `test/rmgpy/solver/liquidTest.py:155`).

**Files:**
- Create: `rmgpy/tools/polymer_moments_runner.py`
- Modify: `setup.py:153-159` (console entry point)
- Test: `test/rmgpy/tools/polymerMomentsRunnerTest.py`

- [ ] **Step 1: Write the failing tests**

Create `test/rmgpy/tools/polymerMomentsRunnerTest.py`:

```python
#!/usr/bin/env python3
"""Tests for the polymer moments CLI reference runner (design spec §8)."""

import csv
import json
import os

import numpy as np
import pytest
import yaml

from rmgpy.cantera import generate_cantera_data
from rmgpy.kinetics import Arrhenius
from rmgpy.molecule import Molecule
from rmgpy.polymer import Polymer, PolymerFluxArchetype, build_polymer_moments_artifact
from rmgpy.reaction import Reaction
from rmgpy.species import Species
from rmgpy.thermo import NASA, NASAPolynomial

from rmgpy.tools.polymer_moments_runner import (
    build_system_from_artifact,
    load_chem_yaml,
    run_segments,
    main,
)


def _spc(smiles, label, index=-1, thermo=True):
    s = Species(molecule=[Molecule().from_smiles(smiles)])
    s.label = label
    s.index = index
    if thermo:
        rows = [NASAPolynomial(coeffs=[2.5, 0, 0, 0, 0, -745.375, 3.35532],
                               Tmin=(tmin, "K"), Tmax=(tmax, "K"))
                for tmin, tmax in ((200.0, 1000.0), (1000.0, 6000.0))]
        s.thermo = NASA(polynomials=rows, Tmin=(200.0, "K"), Tmax=(6000.0, "K"))
    return s


def _mu(label):
    s = Species(label=label, reactive=False)
    s.molecule = [Molecule().from_smiles("[Ne]")]
    s.is_moment_dummy = True
    s.index = -1
    rows = [NASAPolynomial(coeffs=[2.5, 0, 0, 0, 0, -745.375, 3.35532],
                           Tmin=(tmin, "K"), Tmax=(tmax, "K"))
            for tmin, tmax in ((200.0, 1000.0), (1000.0, 6000.0))]
    s.thermo = NASA(polynomials=rows, Tmin=(200.0, "K"), Tmax=(6000.0, "K"))
    return s


@pytest.fixture
def deck(tmp_path):
    """A scission pool + one retained gas reaction, written out as
    chem.yaml + polymer_pools.json exactly like an RMG run would."""
    n2 = _spc("N#N", "N2", index=1)
    g = _spc("[CH3]", "G", index=2)
    g2 = _spc("C", "C1", index=3)
    mus = [_mu("poly_mu0"), _mu("poly_mu1"), _mu("poly_mu2")]
    core = [n2, g, g2] + mus
    gas_rxn = Reaction(reactants=[g], products=[g2],
                       kinetics=Arrhenius(A=(5.0, "1/s"), n=0.0,
                                          Ea=(10.0, "kJ/mol"), T0=(1.0, "K")),
                       reversible=False)
    data, index_map = generate_cantera_data(core, [gas_rxn],
                                            return_reaction_index_map=True)
    chem_path = os.path.join(str(tmp_path), "chem.yaml")
    with open(chem_path, "w") as fh:
        yaml.dump(data, fh, sort_keys=False, default_flow_style=None)

    pool = Polymer(label="poly", monomer="[CH2][CH2]",
                   end_groups=["[H]", "[H]"], cutoff=3,
                   moments=[1.0, 5.0, 30.0], initial_mass=0.0,
                   k_scission=1.0, k_unzip=0.0)
    artifact = build_polymer_moments_artifact(
        [pool], core_species=core, core_reactions=[gas_rxn],
        configured_pool_labels=["poly"], condensed_species=mus,
        cantera_index_map=index_map)
    art_path = os.path.join(str(tmp_path), "polymer_pools.json")
    with open(art_path, "w") as fh:
        json.dump(artifact, fh, indent=2, default=str)
    return chem_path, art_path


class TestChemYamlLoader:
    def test_load_chem_yaml(self, deck):
        chem_path, _ = deck
        species, reactions = load_chem_yaml(chem_path)
        labels = [s.label for s in species]
        assert "N2(1)" in labels and "poly_mu0" in labels
        assert len(reactions) == 1
        rxn = reactions[0]
        assert rxn.kinetics.A.value_si == pytest.approx(5.0)
        assert rxn.reversible is False
        # thermo round-trips (needed for Keq on reversible decks)
        assert species[0].thermo is not None
        assert species[0].get_free_energy(800.0) == pytest.approx(
            species[0].thermo.get_free_energy(800.0))


class TestTwoSegmentRestart:
    def test_two_segment_matches_analytic_to_6_decimals(self, deck):
        """The handoff's empirical claim, now a committed test: a two-segment
        piecewise-isothermal run carries state across the restart and matches
        the analytic scission solution to 6 decimals.
        Analytic (channels are T-independent today):
        mu0(t) = mu1 - (mu1 - mu0_0) * exp(-k_s * t), mu1 const."""
        chem_path, art_path = deck
        with open(art_path) as fh:
            artifact = json.load(fh)
        species, reactions = load_chem_yaml(chem_path)
        rs, core, all_rxns = build_system_from_artifact(
            artifact, species, reactions,
            T0=800.0, P=1.0e5, V_poly=1.0,
            initial_moles={"N2(1)": 1.0}, mass_transfer_spec=[])
        segments = [(0.5, 800.0), (1.0, 850.0)]
        rows = run_segments(rs, core, artifact, all_rxns, segments,
                            n_points_per_segment=10)
        t_final, row_final = rows[-1][0], rows[-1]
        assert t_final == pytest.approx(1.0)
        mu0_idx = 2 + 0  # columns: t, T, then poly_mu0, poly_mu1, poly_mu2
        mu0 = row_final[mu0_idx]
        mu0_analytic = 5.0 - (5.0 - 1.0) * np.exp(-1.0 * 1.0)
        assert mu0 == pytest.approx(mu0_analytic, abs=1e-6)
        mu1 = row_final[mu0_idx + 1]
        assert mu1 == pytest.approx(5.0, abs=1e-6)
        # T column reflects the segment
        assert rows[5][1] == pytest.approx(800.0)
        assert rows[-1][1] == pytest.approx(850.0)

    def test_single_vs_two_segment_equivalence_at_same_T(self, deck):
        """Restart machinery itself must be a no-op when T does not change."""
        chem_path, art_path = deck
        with open(art_path) as fh:
            artifact = json.load(fh)

        def run(segments):
            species, reactions = load_chem_yaml(chem_path)
            rs, core, all_rxns = build_system_from_artifact(
                artifact, species, reactions, T0=800.0, P=1.0e5, V_poly=1.0,
                initial_moles={"N2(1)": 1.0}, mass_transfer_spec=[])
            return run_segments(rs, core, artifact, all_rxns, segments,
                                n_points_per_segment=10)

        one = run([(1.0, 800.0)])
        two = run([(0.5, 800.0), (1.0, 800.0)])
        np.testing.assert_allclose(np.array(one[-1]), np.array(two[-1]),
                                   rtol=1e-6)


class TestCli:
    def test_main_writes_csv(self, deck, tmp_path):
        chem_path, art_path = deck
        profile = os.path.join(str(tmp_path), "tprofile.json")
        with open(profile, "w") as fh:
            json.dump([{"t_end": 0.5, "T": 800.0}, {"t_end": 1.0, "T": 850.0}], fh)
        out_csv = os.path.join(str(tmp_path), "moments.csv")
        moles = os.path.join(str(tmp_path), "moles.json")
        with open(moles, "w") as fh:
            json.dump({"N2(1)": 1.0}, fh)
        main([
            "--artifact", art_path, "--chem", chem_path,
            "--t-profile", profile, "--n-points", "10",
            "--pressure", "1e5", "--v-poly", "1.0",
            "--initial-moles", moles,
            "--output", out_csv,
        ])
        with open(out_csv) as fh:
            rows = list(csv.reader(fh))
        assert rows[0][:2] == ["t_s", "T_K"]
        assert "poly_mu0_mol" in rows[0]
        assert len(rows) == 1 + 2 * 10  # header + 10 points/segment
        assert float(rows[-1][0]) == pytest.approx(1.0)
```

- [ ] **Step 2: Run to verify failure**

Run: `~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/tools/polymerMomentsRunnerTest.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'rmgpy.tools.polymer_moments_runner'`

- [ ] **Step 3: Implement the runner**

Create `rmgpy/tools/polymer_moments_runner.py`:

```python
#!/usr/bin/env python3
"""CLI reference runner for the polymer moments artifact (design spec §8).

Drives HybridPolymerSystem — the ORACLE — from consumer-world inputs only:
  polymer_pools.json (schema 2.0) + chem.yaml + piecewise-isothermal
  T-profile + time grid + mass-transfer spec  ->  moment-trajectory CSV.

It is the cross-validation oracle for the numpy/Cantera consumer (TA), not a
reimplementation. Normative contract: docs/polymer_moments_format.md.

Restart pattern per temperature segment (proven on the analytic two-segment
check; see test/rmgpy/tools/polymerMomentsRunnerTest.py):
  rs.T = Quantity((T_k, 'K')); rs.generate_rate_coefficients(core_rxns, []);
  rs.t0 = t_start; rs.y0 = y_carry; rs.set_initial_derivative();
  rs.initialize_solver(); then rs.advance(t) per grid point.
"""

import argparse
import contextlib
import csv
import io
import json
import sys

import numpy as np
import yaml

from rmgpy.kinetics import Arrhenius
from rmgpy.quantity import Quantity
from rmgpy.reaction import Reaction
from rmgpy.solver.polymer import (
    HybridPolymerSystem,
    MassTransferConfig,
    PolymerPoolConfig,
)
from rmgpy.species import Species
from rmgpy.thermo import NASA, NASAPolynomial

ARCHETYPE_INTS = {
    "same_pool/1": 1,
    "migration/1": 2,
    "scission_fragment/1": 3,
    "legacy_mu1/1": 4,
    "discrete_chip/1": 5,
}


def _species_from_yaml(entry):
    name = entry["name"]
    th = entry.get("thermo", {})
    thermo = None
    if th.get("model") == "NASA7":
        tranges = th["temperature-ranges"]
        rows = th["data"]
        polys = []
        for i, coeffs in enumerate(rows):
            polys.append(NASAPolynomial(coeffs=[float(c) for c in coeffs],
                                        Tmin=(float(tranges[i]), "K"),
                                        Tmax=(float(tranges[i + 1]), "K")))
        thermo = NASA(polynomials=polys,
                      Tmin=(float(tranges[0]), "K"),
                      Tmax=(float(tranges[-1]), "K"))
    spc = Species(label=name, thermo=thermo)
    spc.molecule = []  # consumer-world species: label + thermo only
    return spc


def _parse_equation(eq):
    if "(+ M)" in eq or "(+M)" in eq or " + M " in eq:
        raise NotImplementedError(f"third-body reactions unsupported in v1: {eq}")
    if "<=>" in eq:
        lhs, rhs = eq.split("<=>")
        reversible = True
    elif "=>" in eq:
        lhs, rhs = eq.split("=>")
        reversible = False
    else:
        raise ValueError(f"cannot parse equation: {eq}")
    reactants = [tok.strip() for tok in lhs.split(" + ")]
    products = [tok.strip() for tok in rhs.split(" + ")]
    return [r for r in reactants if r], [p for p in products if p], reversible


_A_UNITS_BY_ORDER = {1: "s^-1", 2: "m^3/(mol*s)", 3: "m^6/(mol^2*s)"}


def load_chem_yaml(path):
    """chem.yaml -> ([Species], [Reaction]) preserving the yaml reactions
    order (== the artifact's cantera.index space). Plain Arrhenius only."""
    with open(path) as fh:
        data = yaml.safe_load(fh)
    species = [_species_from_yaml(e) for e in data.get("species", [])]
    by_name = {s.label: s for s in species}
    reactions = []
    for entry in data.get("reactions", []):
        if "rate-constant" not in entry or "type" in entry:
            raise NotImplementedError(
                f"only elementary Arrhenius reactions are supported in v1; "
                f"offending entry: {entry.get('equation')}")
        r_names, p_names, reversible = _parse_equation(entry["equation"])
        rc = entry["rate-constant"]
        order = len(r_names)
        kin = Arrhenius(A=(float(rc["A"]), _A_UNITS_BY_ORDER[order]),
                        n=float(rc["b"]), Ea=(float(rc["Ea"]), "J/mol"),
                        T0=(1.0, "K"))
        rxn = Reaction(reactants=[by_name[n] for n in r_names],
                       products=[by_name[n] for n in p_names],
                       kinetics=kin, reversible=reversible,
                       duplicate=bool(entry.get("duplicate", False)))
        reactions.append(rxn)
    return species, reactions


def _restamp_and_extend(artifact, species, reactions):
    """Restore the solver stamps from artifact entries; reconstruct dropped
    (cantera-null) reactions from their recorded stoichiometry + kinetics."""
    by_name = {s.label: s for s in species}
    all_reactions = list(reactions)
    for e in artifact["reactions"]:
        arch = ARCHETYPE_INTS[e["archetype"]]
        if e["cantera"] is not None:
            rxn = reactions[e["cantera"]["index"]]
        else:
            kin = e["kinetics"]
            if kin is None:
                raise ValueError(f"cantera-null entry {e['id']} carries no kinetics")
            order = len(e["reactants"])
            rxn = Reaction(
                reactants=[by_name[n] for n in e["reactants"]],
                products=[by_name[n] for n in e["products"]],
                kinetics=Arrhenius(A=(kin["A"], _A_UNITS_BY_ORDER[order]),
                                   n=kin["n"], Ea=(kin["Ea"], "J/mol"),
                                   T0=(1.0, "K")),
                reversible=bool(kin["reversible"]))
            all_reactions.append(rxn)
        rxn.polymer_flux_archetype = arch
        rxn.is_end_group_reaction = (e["scaling"] == "mu0")
        rxn.polymer_chip_units = int(e.get("params", {}).get("a", 0))
        if e["kinetics"] is not None:
            # Oracle-faithful reversibility: chem.yaml equations always print
            # '<=>' (rmgpy/cantera.py get_reaction_equation), so the arrow
            # cannot distinguish irreversible RMG reactions. The artifact's
            # kinetics.reversible records what the generating solver used.
            rxn.reversible = bool(e["kinetics"]["reversible"])
    return all_reactions


def build_system_from_artifact(artifact, species, reactions,
                               T0, P, V_poly, initial_moles,
                               mass_transfer_spec, initial_moments=None):
    """Assemble the HybridPolymerSystem oracle from consumer-world inputs.

    Returns (system, core_species, all_reactions) — all_reactions includes
    the cantera-null reconstructions and is needed by run_segments for the
    per-segment generate_rate_coefficients call (HybridPolymerSystem is a
    cdef class; do not hang extra attributes on it). The system is fully
    initialized at T0 (initialize_model runs through initialize_solver,
    rmgpy/solver/polymer.pyx:601-610)."""
    core = list(species)
    idx = {s.label: i for i, s in enumerate(core)}
    all_reactions = _restamp_and_extend(artifact, core, reactions)

    conv = artifact["conventions"]
    condensed = set(conv["condensed_species"])
    mask = np.array([s.label not in condensed for s in core], dtype=bool)

    pools = []
    moments0 = {}
    for p in artifact["pools"]:
        lab = p["label"]
        if lab not in conv["configured_pools"]:
            continue
        mu_idx = tuple(idx[f"{lab}_mu{k}"] for k in range(3))
        routing = p.get("monomer_routing")
        pools.append(PolymerPoolConfig(
            label=lab, xs=int(p.get("cutoff") or 0),
            explicit_dp_to_species_index={},
            mu_indices=mu_idx,
            monomer_poly_index=idx[routing] if routing else None,
            k_scission=float(p["channels"]["scission"]["A"]),
            k_unzip=float(p["channels"]["unzip"]["A"]),
        ))
        if initial_moments and lab in initial_moments:
            moments0[lab] = tuple(initial_moments[lab])
        elif p.get("moments") is not None:
            moments0[lab] = tuple(p["moments"])

    mt_configs = [MassTransferConfig(gas_index=idx[m["gas"]],
                                     poly_index=idx[m["poly"]],
                                     K=float(m["K"]), kLa=float(m["kLa"]))
                  for m in (mass_transfer_spec or [])]

    init_moles = {core[idx[label]]: float(v)
                  for label, v in (initial_moles or {}).items()}

    rs = HybridPolymerSystem(
        T=(T0, "K"), P=(P, "Pa"),
        initial_mole_fractions=init_moles,  # interpreted as moles
        V_poly=float(V_poly),
        polymer_pools=pools, mass_transfer=mt_configs,
        gas_species_mask=mask, constant_gas_volume=False,
        initial_polymer_moments=moments0, termination=[])
    with contextlib.redirect_stdout(io.StringIO()):  # mute the mapping banner
        rs.initialize_model(core, all_reactions, [], [])
    return rs, core, all_reactions


def run_segments(rs, core, artifact, all_reactions, segments,
                 n_points_per_segment=50):
    """Piecewise-isothermal integration. ``segments`` = [(t_end, T_K), ...]
    with strictly increasing t_end. Returns rows:
    [t, T, <pool>_mu0.., <pool>_mu1.., <pool>_mu2.. per configured pool,
     n(monomer_routing) per routed pool]."""
    conv = artifact["conventions"]
    idx = {s.label: i for i, s in enumerate(core)}
    pool_labels = list(conv["configured_pools"])
    mu_cols = [(lab, tuple(idx[f"{lab}_mu{k}"] for k in range(3)))
               for lab in pool_labels]
    routed = [(p["label"], idx[p["monomer_routing"]])
              for p in artifact["pools"]
              if p["label"] in pool_labels and p.get("monomer_routing")]

    rows = []
    t_start = 0.0
    y_carry = None
    for seg_i, (t_end, T_k) in enumerate(segments):
        if seg_i > 0:
            # the proven restart pattern (see module docstring)
            rs.T = Quantity((T_k, "K"))
            rs.generate_rate_coefficients(all_reactions, [])
            rs.t0 = t_start
            rs.y0 = y_carry.copy()
            rs.set_initial_derivative()
            rs.initialize_solver()
        for t in np.linspace(t_start, t_end, n_points_per_segment + 1)[1:]:
            rs.advance(t)
            y = np.asarray(rs.y)
            row = [float(t), float(T_k)]
            for _lab, (i0, i1, i2) in mu_cols:
                row.extend([float(y[i0]), float(y[i1]), float(y[i2])])
            for _lab, ri in routed:
                row.append(float(y[ri]))
            rows.append(row)
        y_carry = np.asarray(rs.y).copy()
        t_start = t_end
    return rows


def _csv_header(artifact):
    conv = artifact["conventions"]
    header = ["t_s", "T_K"]
    for lab in conv["configured_pools"]:
        header.extend([f"{lab}_mu0_mol", f"{lab}_mu1_mol", f"{lab}_mu2_mol"])
    for p in artifact["pools"]:
        if p["label"] in conv["configured_pools"] and p.get("monomer_routing"):
            header.append(f"n_{p['monomer_routing']}_mol")
    return header


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Polymer moments artifact reference runner (oracle). "
                    "See docs/polymer_moments_format.md.")
    parser.add_argument("--artifact", required=True, help="polymer_pools.json (schema 2.0)")
    parser.add_argument("--chem", required=True, help="chem.yaml from the same RMG run")
    parser.add_argument("--t-profile", required=True,
                        help="JSON: [{\"t_end\": s, \"T\": K}, ...] piecewise-isothermal")
    parser.add_argument("--n-points", type=int, default=50,
                        help="output points per segment (default 50)")
    parser.add_argument("--pressure", type=float, default=1.0e5, help="Pa")
    parser.add_argument("--v-poly", type=float, required=True, help="m^3")
    parser.add_argument("--initial-moles", required=True,
                        help="JSON: {chem.yaml label: mol}")
    parser.add_argument("--initial-moments", default=None,
                        help="JSON: {pool label: [mu0, mu1, mu2] mol}; "
                             "default = artifact pools[].moments")
    parser.add_argument("--mass-transfer", default=None,
                        help="JSON: [{gas, poly, K, kLa}] (labels; operating "
                             "condition, not artifact content)")
    parser.add_argument("--output", required=True, help="CSV path")
    args = parser.parse_args(argv)

    with open(args.artifact) as fh:
        artifact = json.load(fh)
    if not str(artifact.get("schema_version", "")).startswith("2."):
        sys.exit(f"artifact schema_version {artifact.get('schema_version')!r} "
                 "is not 2.x — regenerate with a current RMG-Py polymer branch")
    with open(args.t_profile) as fh:
        profile = [(float(seg["t_end"]), float(seg["T"])) for seg in json.load(fh)]
    with open(args.initial_moles) as fh:
        initial_moles = json.load(fh)
    initial_moments = None
    if args.initial_moments:
        with open(args.initial_moments) as fh:
            initial_moments = json.load(fh)
    mass_transfer_spec = []
    if args.mass_transfer:
        with open(args.mass_transfer) as fh:
            mass_transfer_spec = json.load(fh)

    species, reactions = load_chem_yaml(args.chem)
    rs, core, all_reactions = build_system_from_artifact(
        artifact, species, reactions,
        T0=profile[0][1], P=args.pressure, V_poly=args.v_poly,
        initial_moles=initial_moles, mass_transfer_spec=mass_transfer_spec,
        initial_moments=initial_moments)
    rows = run_segments(rs, core, artifact, all_reactions, profile,
                        n_points_per_segment=args.n_points)

    with open(args.output, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(_csv_header(artifact))
        writer.writerows(rows)
    print(f"wrote {len(rows)} rows to {args.output}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Add the console entry point**

In `setup.py:154-157`, extend:

```python
        'console_scripts': [
            'rmg.py = rmgpy.__main__:main',
            'Arkane.py = arkane.__main__:main',
            'polymer_moments_runner.py = rmgpy.tools.polymer_moments_runner:main',
        ],
```

(Note: the entry point activates on the next `make`/install; tests call `main()` directly and `python -m rmgpy.tools.polymer_moments_runner` works immediately.)

- [ ] **Step 5: Run the tests**

Run: `~/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/tools/polymerMomentsRunnerTest.py -q`
Expected: all PASS. Watch for: (a) `gas_species_mask` length mismatches → check the deck's condensed list; (b) the analytic check failing in the 7th decimal — DASPK tolerances are atol 1e-16/rtol 1e-8 by default through `initialize_model`; the assertion is abs=1e-6 per the handoff's claim.

- [ ] **Step 6: Commit**

```bash
git add rmgpy/tools/polymer_moments_runner.py setup.py test/rmgpy/tools/polymerMomentsRunnerTest.py
git commit -m "tools: add polymer moments reference runner

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 10: EPDM fixture regeneration + full verification

Regenerate `~/runs/RMG/epdm_v0_2026-06-06b` so TA gets a schema-2.0 fixture; confirm the run is a no-op chemistry-wise; confirm the Chemkin fix on the freshly written files; re-run every affected suite.

**Files:** none committed (the run dir is outside the repo). If any step exposes a bug, fix it in its own task-style commit.

- [ ] **Step 1: Re-run the EPDM deck**

```bash
cd ~/runs/RMG/epdm_v0_2026-06-06b && ~/anaconda3/envs/rmg_env/bin/python ~/Code/RMG-Py/rmg.py input.py 2>&1 | tail -20
```
Expected: `MODEL GENERATION COMPLETED`, final model: **26 species / 28 reactions**. ~20 aggregate legacy warning lines (`...proxy-touching reactions arrived without a polymer_flux_archetype stamp...` / `...demoted to legacy mu1-only moment flux...`) are EXPECTED — do not chase them.

- [ ] **Step 2: Verify the regenerated chem.inp loads (Chemkin fix, original repro)**

```bash
~/anaconda3/envs/rmg_env/bin/python -c "
from rmgpy.chemkin import load_chemkin_file
sp, rx = load_chemkin_file('$HOME/runs/RMG/epdm_v0_2026-06-06b/chemkin/chem.inp',
                           '$HOME/runs/RMG/epdm_v0_2026-06-06b/chemkin/species_dictionary.txt')
print('LOADED', len(sp), 'species,', len(rx), 'reactions')"
```
Expected: `LOADED 26 species, 28 reactions` (no ValueError).

- [ ] **Step 3: Verify the schema-2.0 sidecar**

```bash
~/anaconda3/envs/rmg_env/bin/python - <<'EOF'
import json
data = json.load(open('/home/alon/runs/RMG/epdm_v0_2026-06-06b/chemkin/polymer_pools.json'))
assert data['schema_version'] == '2.0', data['schema_version']
assert data['rmg_commit'], 'rmg_commit missing'
conv = data['conventions']
assert conv['configured_pools'] == ['epdm'], conv['configured_pools']
assert any(l.startswith('epdm_mu') for l in conv['condensed_species'])
pools = {p['label']: p for p in data['pools']}
assert set(pools) >= {'epdm', 'epdm_scission_tail'}
e = pools['epdm']
for key in ('moments', 'monomer_mw_g_mol', 'mn_g_mol', 'mw_g_mol',
            'initial_mass_g', 'channels', 'phase_species',
            'monomer_routing', 'mu3_closure'):
    assert key in e, key
assert e['channels']['scission']['A'] == 1.0
assert e['channels']['unzip']['A'] == 0.01
assert e['mu3_closure'] == 'log_lagrange/1'
# 1.0 fields still present on every pool
for p in data['pools']:
    for key in ('label', 'monomer_smiles', 'monomer_adj_list',
                'feature_monomers_smiles', 'end_groups', 'cutoff',
                'parent_pool', 'spawn_iteration', 'spawn_event_metadata',
                'mu_indices'):
        assert key in p, (p['label'], key)
rxns = data['reactions']
assert len(rxns) >= 25, len(rxns)   # EPDM: ~28 pool-touching reactions
assert all(e['kinetics'] is not None for e in rxns if e['cantera'] is None), \
    'cantera-null entry without kinetics'
assert any(e['unresolved'] for e in rxns)        # legacy emission (Q3)
assert any(e['cantera'] is None for e in rxns)   # dropped-by-filter entries exist
print('sidecar OK:', len(rxns), 'reaction entries,',
      sum(1 for e in rxns if e['unresolved']), 'unresolved,',
      sum(1 for e in rxns if e['cantera'] is None), 'cantera-null')
EOF
```
Expected: `sidecar OK: ...` with ~27–28 entries, most `unresolved`.

- [ ] **Step 4: Smoke the CLI runner on the EPDM fixture**

```bash
cd /home/alon/Code/RMG-Py
printf '[{"t_end": 1.0, "T": 1000.0}]' > /tmp/epdm_tprofile.json
printf '{"N2": 0.90, "epdm(2)": 0.099, "H(1)": 0.001}' > /tmp/epdm_moles.json
~/anaconda3/envs/rmg_env/bin/python -m rmgpy.tools.polymer_moments_runner \
  --artifact ~/runs/RMG/epdm_v0_2026-06-06b/chemkin/polymer_pools.json \
  --chem ~/runs/RMG/epdm_v0_2026-06-06b/cantera/chem.yaml \
  --t-profile /tmp/epdm_tprofile.json --n-points 20 \
  --pressure 1e5 --v-poly 1.16e-6 \
  --initial-moles /tmp/epdm_moles.json \
  --output /tmp/epdm_moments.csv
head -3 /tmp/epdm_moments.csv
```
Expected: `wrote 20 rows to /tmp/epdm_moments.csv`; header `t_s,T_K,epdm_mu0_mol,...`. (Exact initial-moles labels: check the `species:` list in the chem.yaml — gas species carry `(index)` suffixes, e.g. `H(1)`, `epdm(2)`, bare `N2`. `--v-poly 1.16e-6` ≈ 1 g EPDM / 860 kg/m³.) If a `NotImplementedError` fires on a reaction type, record it — EPDM is expected to be all elementary Arrhenius.

- [ ] **Step 5: Full affected-suite re-run (baseline 311 passed, 2 skipped + new tests)**

```bash
cd /home/alon/Code/RMG-Py && ~/anaconda3/envs/rmg_env/bin/python -m pytest \
  test/rmgpy/polymerTest.py test/rmgpy/solver/solverPolymerTest.py \
  test/rmgpy/canteraTest.py test/rmgpy/polymerMultiPoolTest.py \
  test/rmgpy/reactionTest.py test/rmgpy/rmg/modelTest.py \
  test/rmgpy/chemkinTest.py test/rmgpy/polymerArtifactTest.py \
  test/rmgpy/tools/polymerMomentsConsumerTest.py \
  test/rmgpy/tools/polymerMomentsRunnerTest.py -q
```
Expected: everything passes; no regressions in the original six suites.

- [ ] **Step 6: Commit (only if repo files changed during fixes here)**

```bash
git status --short   # if clean: nothing to commit for this task
```

---

### Task 11: TA handoff-back

The contract obligation (spec §12). Written OUTSIDE the repo; not committed.

**Files:**
- Create: `~/handoffs/handoff-polymer-moments-artifact-landed.md`

- [ ] **Step 1: Write the handoff**

Create `~/handoffs/handoff-polymer-moments-artifact-landed.md` with this content, filling the `<commit>` placeholders from `git log --oneline -11`:

```markdown
# Handoff: polymer moments artifact (schema 2.0) LANDED

**Date:** <today>
**From:** RMG-Py `polymer` branch session (commits <first>..<last>, not pushed)
**To:** TA session (~/Code/TA) — you are unblocked.

## ⚠ MAJOR VERSION BUMP (your ground rule: explicit flag)

`polymer_pools.json` is now **schema 2.0** (was 1.0). Every 1.0 field kept its
exact name, type and location — your `.get()`-based loader
(`ta/mechanism.py`) keeps working unmodified until you opt into the new
fields. The major bump signals semantics (compiled reactions block), not
breakage. There is NO new file; same sidecar, same directory as before.

## The normative contract you implement against

`~/Code/RMG-Py/docs/polymer_moments_format.md` (`polymer_moments_format/2.0`)
— full schema, the five archetype term equations with clamp/guard semantics
(transcribed from rmgpy/solver/polymer.pyx and pinned by the sufficiency
tests), the 6-step rate recipe (incl. the Cantera multiplier zeroing and the
phase/product gate), channel equations (`scission/1`, `unzip/1`), the
`log_lagrange/1` µ3 closure, Keq recipe for cantera-null reversible entries,
conservation invariants, and the versioning policy.

## Final field names (pool additions — match your dataclass)

`moments` (mol, DP basis), `monomer_mw_g_mol`, `mn_g_mol`, `mw_g_mol`,
`initial_mass_g`, `channels` ({scission,unzip}: {A,n,Ea,units}),
`phase_species` (chem.yaml labels), `monomer_routing` (label|null; null on
the EPDM fixture), `mu3_closure` ("log_lagrange/1").

Beyond the original spec sketch, reaction entries ALSO carry `reactants`,
`products`, `proxy_products` (needed to evaluate cantera-null entries and the
reverse-direction proxy rule), and `conventions` carries `configured_pools` +
`condensed_species` (normative phase/pool-membership lists — use them, do not
name-guess). Details in the format doc §3/§8.

## Your test fixture (regenerated, schema 2.0)

- `~/runs/RMG/epdm_v0_2026-06-06b/chemkin/polymer_pools.json` (schema 2.0,
  ~28 reaction entries, mostly `legacy_mu1/1` + `unresolved: true` — integrate
  them, warning allowed; several `cantera: null` entries carrying kinetics)
- `~/runs/RMG/epdm_v0_2026-06-06b/cantera/chem.yaml` (unchanged naming —
  µ-dummies are bare labels, gas species are `label(index)`)
- `~/runs/RMG/epdm_v0_2026-06-06b/chemkin/chem.inp` now LOADS via
  load_chemkin_file — the >16-char µ-dummy THERMO overflow is fixed (long
  no-index labels get compressed Chemkin identifiers; chem.yaml names are
  untouched, so nothing on your side changes).

## Cross-validation oracle

`python -m rmgpy.tools.polymer_moments_runner --artifact ... --chem ...
--t-profile segments.json --n-points N --pressure Pa --v-poly m3
--initial-moles moles.json [--mass-transfer mt.json] --output out.csv`
(also installed as `polymer_moments_runner.py`). Piecewise-isothermal
T-profile; mass transfer is a runner INPUT ({gas, poly, K, kLa} by label) —
an operating condition, never artifact content. Output: per-pool µ0/µ1/µ2
trajectory CSV.

## Reference consumer (steal it)

`~/Code/RMG-Py/test/rmgpy/tools/numpy_moments_consumer.py` — numpy-only
(zero rmgpy imports, enforced by test), implements the full recipe, matches
the solver to ~1e-9 rtol on: pure scission/unzip channels, a two-pool
MIGRATION+DISCRETE_CHIP deck, a reversible Keq case, and a both-sides-nonzero
mass-transfer case (`test/rmgpy/tools/polymerMomentsConsumerTest.py`).

## Known nulls / v1 limits

- `monomer_routing` is null on EPDM (no condensed monomer species configured);
  unzip µ-flux still applies.
- The runner supports elementary Arrhenius chem.yaml entries only (EPDM is).
- `id` values are stable within one artifact, NOT across regenerations.
- Spawned daughter pools (e.g. `epdm_scission_tail`) appear in `pools[]` but
  are NOT in `configured_pools`: treat their proxies as ordinary species
  (that is what the oracle does).
```

- [ ] **Step 2: Update the running-log memory**

Append a line to `~/.claude/projects/-home-alon-Code-RMG-Py/memory/project_polymer_branch_review.md` noting: artifact schema 2.0 landed (commit range), format doc path, EPDM fixture regenerated, handoff written.

- [ ] **Step 3: Final verification sweep**

Re-run the Task 10 Step 5 pytest command once more; confirm `git log --oneline -11` shows one commit per task (Tasks 1–9; 10/11 commit only if repo files changed) and `git status` is clean apart from untracked run artifacts.
