# Declarative Polymer Moment-Kinetics Artifact — Design Spec (DRAFT, awaiting user review)

**Date:** 2026-06-10
**Branch:** `polymer`
**Origin:** feature request from the TA session
(`/home/alon/handoffs/handoff-polymer-declarative-moments-format.md`); TA
development is design-blocked on this schema.
**Goal:** RMG emits a versioned declarative artifact compiling the
`HybridPolymerSystem` moment dispatch into data, so a consumer with only
numpy + Cantera integrates per-pool µ0/µ1/µ2 ODEs coupled to the gas solution
— no rmgpy import. The same boundary pattern as chem.yaml ⇄ Cantera:
RMG's classification/stamping logic keeps churning; the artifact contract is
a small closed vocabulary.

## 1. Decisions (settled with the user, 2026-06-10)

| # | Decision | Rationale |
|---|----------|-----------|
| Q1 | **`polymer_pools.json` schema 2.0** — single sidecar, major version bump (not a new file) | User's call. One artifact to plumb; TA's loader is `.get()`-based with no version check, so the bump is mechanically safe; the handoff-back to TA must flag the major bump explicitly (its ground rule) |
| Q2 | **Archetype-level term vocabulary** — closed set mirroring the solver dispatch one-to-one | Bundle equations (incl. clamps) are normative in the format doc; per-reaction entries carry parameters only. New archetype = new term type = minor version bump. Explicit term algebra rejected: bigger surface, clamp semantics need names anyway |
| Q3 | **UNRESOLVED/legacy-µ1 reactions ARE emitted**, as `legacy_mu1/1` with `unresolved: true` | The acceptance criterion is trajectory match with the solver; EPDM has ~27 such reactions. Consumer may warn but must integrate them |
| Q4 | **CLI reference runner ships in v1.0** | User's call. Artifact + chem.yaml + piecewise-isothermal T-profile in, moment-trajectory CSV out, driving `HybridPolymerSystem` directly — an oracle binary for TA cross-validation, not a reimplementation |

## 2. Envelope (schema 2.0)

```json
{
  "schema_version": "2.0",
  "generated_at": "<ISO-8601 UTC>",
  "rmg_commit": "<git SHA of the emitting RMG-Py checkout>",
  "rmg_iteration": <int>,
  "conventions": { ... },          // §5 — rate-evaluation recipe, bases, units
  "pools": [ ... ],                // §3 — extended per-pool blocks
  "reactions": [ ... ]             // §4 — compiled per-reaction flux terms
}
```

**Backward compatibility:** every schema-1.0 field keeps its exact name,
type, and location (`label`, `monomer_smiles`, `monomer_adj_list`,
`feature_monomers_smiles`, `end_groups`, `cutoff`, `parent_pool`,
`spawn_iteration`, `spawn_event_metadata`, `mu_indices`). TA's current loader
(`~/Code/TA/ta/mechanism.py:66-80`, `.get()`-based, ignores unknown keys)
continues to work unmodified until it opts into the new fields. The major
bump signals semantics, not breakage.

## 3. `pools[]` additions (per pool, all absent from 1.0)

| Field | Type | Meaning |
|---|---|---|
| `moments` | [µ0, µ1, µ2] floats | pool state at generation time, **mol, DP basis** (basis/units stated in `conventions`) |
| `monomer_mw_g_mol` | float | repeat-unit MW |
| `mn_g_mol`, `mw_g_mol` | float | number-/weight-average MW at generation time |
| `initial_mass_g` | float \| null | as configured |
| `channels` | `{scission: {A, n, Ea, units}, unzip: {...}}` | Arrhenius-capable; today RMG fills `A=k, n=0, Ea=0`. Channel *equations* are normative by versioned name (`scission/1`, `unzip/1`) referencing `docs/multi_pool_design.md` §5 — not duplicated in the artifact |
| `phase_species` | [labels] | the pool's condensed-phase species: proxy variants, µ-dummies, explicit oligomers, condensed monomer — explicit, so TA stops name-guessing `<label>_muN` |
| `monomer_routing` | label \| null | the chem.yaml species receiving unzip monomer flux (replaces the internal `monomer_poly_index` — labels, not solver state-vector indices) |
| `mu3_closure` | `"log_lagrange/1"` | versioned closure name; equation + realizability guard normative in the format doc |

## 4. `reactions[]` — compiled flux terms

One entry per stamped proxy-touching core reaction:

```json
{
  "id": "<stable string id>",
  "cantera": {"index": <int>, "equation": "<string>"} | null,
  "kinetics": {"A": ..., "n": ..., "Ea": ..., "units": ..., "reversible": bool},
  "proxy_reactants": ["<species label>", ...],
  "scaling": "mu0" | "mu1",
  "src_pool": "<pool label>",
  "dst_pool": "<pool label>" | null,
  "archetype": "same_pool/1" | "migration/1" | "scission_fragment/1"
             | "discrete_chip/1" | "legacy_mu1/1",
  "params": {"a": <int>},          // discrete_chip only
  "unresolved": bool               // true => legacy_mu1 fallback emission (Q3)
}
```

- `cantera.index` is authoritative (the reaction's position in the emitted
  chem.yaml); `equation` is a human-readable checksum. `cantera: null` means
  the reaction was dropped by the unbalanced-proxy filter — then `kinetics`
  is REQUIRED, since the consumer cannot get its rate from Cantera.
- Bundle equations per term type (including the chip forward clamp
  `max(0, 2a·E[n] − a²)` and the unclamped reverse `+(2a·E[n] + a²)`, the
  MIGRATION per-direction src/dst bundles, the SCISSION_FRAGMENT
  complement-stays forms, and all denominators/µ3 guards) are normative in
  the format doc (§6), transcribed from `rmgpy/solver/polymer.pyx` — the
  artifact carries parameters only.

## 5. `conventions` block — the rate-evaluation recipe

The critical non-obvious semantics (verified against
`polymer.pyx:987-1032`): **retained chem.yaml proxy reactions cannot be
evaluated by Cantera as-is.** The recipe the consumer must implement:

1. Take `kf`/`kb` (rate constants) from Cantera at the current T, P (or from
   the entry's own `kinetics` when `cantera` is null).
2. Form the concentration product over reactants, **treating every species
   in `proxy_reactants` as concentration 1.0**; gas species use `C = n/V_gas`,
   condensed species `C = n/V_poly`.
3. Multiply both directions by the site density
   `max(0, µ_scaling)/V_poly` — µ1, or µ0 when `scaling == "mu0"`.
4. Event rate per direction: `ev_mol = rate · V_rxn` (V_rxn = V_poly for
   pool events); apply the archetype bundle per direction — forward moves
   src-pool statistics, reverse moves dst-pool statistics.

Also stated here: moment basis (mol, DP), volume conventions
(V_gas ideal-gas dynamic, V_poly constant consumer input), and the
conservation invariants the consumer can assert (Σµ1 bookkeeping incl. the
chip-weighted sum `Σ_pools µ1 + Σ_chips a_i·n_i`).

## 6. Format spec doc — `docs/polymer_moments_format.md`

THE normative contract TA implements against: the full 2.0 schema, the
term-type equations with clamp/guard semantics, the rate recipe, channel
equations by reference to `multi_pool_design.md` §5, conservation
invariants, and the versioning policy (minor = additive term types/fields;
major = breaking). The schema churns with RMG; the format doc versions with
the schema.

## 7. Emitter

- Extends `_serialize_pool_for_sidecar` / `write_polymer_pools_sidecar`
  (`rmgpy/polymer.py:2331-2421`) for the pool additions.
- New reaction-compiler function reading the stamps
  (`polymer_flux_archetype`, `is_end_group_reaction`,
  `polymer_chip_units`) off live core reactions at the existing
  `save_everything` hook (`rmgpy/rmg/main.py:2057-2079`).
- The Cantera index mapping comes from the export step, which already runs
  before the sidecar write in the same hook (ordering verified).

## 8. CLI reference runner (v1.0)

`rmgpy/tools/` module + console entry point. Inputs: artifact + chem.yaml +
piecewise-isothermal T-profile + time grid. Output: moment-trajectory CSV.
Drives `HybridPolymerSystem` directly (the proven restart pattern: set
`rs.T`, rebuild pool configs, `generate_rate_coefficients`, set `t0/y0`,
`set_initial_derivative()`, `initialize_solver()` — carries state across
segments to 6 decimals on the analytic two-segment check). It is the oracle,
not a reimplementation.

## 9. Acceptance / tests

1. **Sufficiency proof:** a consumer-side mini-integrator in the test suite
   using ONLY numpy + the artifact (no rmgpy objects) reproduces the
   solver's trajectories on pure-scission and pure-unzip analytic cases.
2. **EPDM fixture:** regenerate `~/runs/RMG/epdm_v0_2026-06-06b` artifacts
   with schema 2.0 as the downstream TA test fixture.
3. Emitter unit tests: stamps → reaction entries mapping (incl. dropped
   reactions carrying kinetics, UNRESOLVED flagged emission).
4. Schema round-trip + 1.0-field stability test (every 1.0 key present,
   unchanged type).
5. **Chemkin >16-char µ-dummy label bug fixed** (handoff: "fix
   regardless") + regression test — long moment-dummy names (e.g.
   `epdm_scission_tail_mu0`, 22 chars) overflow the fixed-column THERMO
   format and break `load_chemkin_file`; route through
   `get_species_identifier` or equivalent. Separate commit.

## 10. Decided without asking (flag if wrong)

- The runner consumes chem.yaml + artifact, NOT RMG input files — it lives
  in the consumer's world.
- `monomer_routing` is a species label, not a state-vector index; solver
  allocation details (`mu_indices`) remain only for 1.0 compatibility.

## 11. Out of scope (v1.1+)

Per-event ΔH for DTA heat-release reconstruction (schema reserves nothing;
additive when it comes); continuous non-isothermal T-profiles in the runner;
alternative closures (log-normal); emitting the artifact mid-run per
iteration (end-of-run only, matching the current sidecar hook).

## 12. Handoff-back obligation

When this lands: write the handoff for the TA session (final schema, format
doc path, EPDM fixture paths, the major-bump flag per Q1) or note it in the
running-log memory — TA resumes from that.
