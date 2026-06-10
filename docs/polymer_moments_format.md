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
- `kinetics` is only guaranteed non-null when the underlying RMG rate is plain
  Arrhenius; a dropped reaction with composite kinetics (MultiArrhenius/PDep)
  emits both `cantera: null` AND `kinetics: null` and cannot be evaluated by
  the consumer (RMG logs a warning at emission). The emitted `A` is normalized
  to `T0 = 1 K` (`A·T^n·exp(−Ea/RT)` evaluates directly).
- `id` is stable within one artifact and across re-reads of the same file;
  NOT across regenerations.
- The pool of a proxy label is the label with any trailing `(N)` stripped
  (e.g. `epdm(2)` → pool `epdm`).
- `unresolved: true` entries (legacy-µ1 emissions, including solver-demoted
  stamps) MUST be integrated; consumers MAY warn.
- For retained entries with composite kinetics, `kinetics` may be `null`; the
  consumer then takes rates from Cantera. In post-fix chem.yaml the equation
  arrow records reversibility (`=>` irreversible, `<=>` reversible), so the
  consumer reads reversibility from the arrow — composite-kinetics retained
  entries no longer need a treat-as-reversible assumption. (In PRE-FIX files,
  written before the cantera `=>`-for-irreversible fix, every equation is
  `<=>`, so such an entry must still be treated as reversible.)
- **Caveat (ordinary reactions):** post-fix chem.yaml records reversibility in
  the equation arrow (`=>` irreversible, `<=>` reversible), so reactions
  WITHOUT an artifact entry round-trip their reversibility faithfully in every
  consumer of the file (Cantera, the reference runner, the numpy consumer).
  *Legacy:* chem.yaml files written by RMG BEFORE this fix print every equation
  as `<=>` with no separate reversibility marker, so an originally-irreversible
  ordinary reaction in such a file runs reversible in every consumer.
  Artifact-listed entries are exempt either way — their `kinetics.reversible`
  restores the generating solver's behavior (and now agrees with the arrow in
  post-fix files).

## 4. Rate-evaluation recipe (normative; oracle: `polymer.pyx:922-1261`)

Definitions: `V_poly` = constant consumer input; `V_gas = n_gas·R·T/P`
(ideal gas over the *gas-phase* species — those NOT in
`conventions.condensed_species`; floor `V_gas = 1.0 m³` when `n_gas ≤ 0`).
Pool moments are extensive mol; intensive µ ≡ y[µ]/V_poly clamped at ≥ 0.
**`R = 8.314472 J/(mol·K)`** — the oracle's `rmgpy.constants.R`
(CODATA-2006), used in `V_gas`, `kf`, `Keq`, and mass transfer alike. Do
NOT substitute the 2018-SI exact value (8.31446261815324): the 1.13e-6
relative offset shifts every gas-phase term off the oracle.

0. **Disable every listed retained reaction in Cantera:** for each entry with
   `cantera != null`, `Kinetics.set_multiplier(0, cantera.index)`. (Equation
   strings are a human checksum only.)
1. Per entry, per T: `kf` from Cantera (or from `kinetics`:
   `kf = A·T^n·exp(−Ea/(R·T))`, SI). If `reversible`: `kb = kf/Keq` with
   `Keq(T) = (P°/(R·T))^Δn · exp(−ΔG°/RT)`, `P° = 1e5 Pa`, ΔG° summed
   from the chem.yaml NASA polynomials (`reaction.py:767-840`). Else `kb = 0`.
   **`Δn` counts ALL species as written in `reactants`/`products`** (with
   repeats), condensed and proxy participants included — the oracle's
   `get_equilibrium_constant` (`reaction.py:805-821`) excludes only *surface
   sites*, which never occur in this artifact. Do NOT compute the exponent
   over the `condensed_species`-complement: that silently mis-scales `kb` by
   `(P°/RT)^±k` for reversible reactions with condensed participants.
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
   `site = max(0, y[µ_scaling])/V_poly` — NOTE: `y[µ]` here is the EXTENSIVE
   moment (raw moles, before the preamble's intensive division), so the single
   `/V_poly` shown is the whole conversion; do not divide twice. µ read from
   `src_pool` (the FIRST proxy-reactant's pool, reactant-slot priority), µ0
   when `scaling=="mu0"` else µ1 — EXCEPT chip entries
   (`archetype=="discrete_chip/1"` ∧ `scaling=="mu0"` ∧ `a > 0`):
   `site = min(max(0, y[µ0]), max(0, y[µ1])/a)/V_poly` (exhaustion throttle).
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
