# Design — wire `is_end_group_reaction` (END_MOD → μ0 scaling)

**Date:** 2026-06-07
**Branch:** `polymer`
**Status:** approved design, pending implementation plan
**Open item resolved:** #6 from `~/handoffs/handoff-polymer-branch-open-decisions.md`
("`is_end_group_reaction` dead code" — decision: WIRE IT UP, not remove).

## Problem

The polymer hybrid solver (`rmgpy/solver/polymer.pyx`) scales a proxy
reaction's rate by a pool moment when a reactant or product is a pool proxy.
Today every proxy reaction scales by **μ1** (monomer-unit / site density). The
solver contains scaffolding to instead scale *end-group* reactions by **μ0**
(chain-end density):

- `self.is_end_group_reaction` is allocated as zeros at `initialize_model`
  (line ~418) and **never assigned 1** anywhere.
- Two consumer branches read it and would switch the scaling moment from μ1 to
  μ0: the RHS residual block (~line 915) and the reaction-rate block
  (~line 1269).

Because the flag is never set, both branches are dead and the μ0 path is
unreachable.

### Why this is wrong chemistry

A reaction that modifies a chain **end** occurs at a rate proportional to the
number of chain ends, i.e. μ0. A reaction on a random backbone repeat unit is
proportional to the number of such units, i.e. μ1. Scaling an end-group
reaction by μ1 over-counts its rate by ~DP (degree of polymerization, μ1/μ0),
which for a long chain is a large error. The scaffolding anticipates the
correct physics; it is simply inert.

## Decision summary

Resolved during brainstorming (2026-06-07):

1. **Scope — which classes scale by μ0:** `END_MOD` only. `BASELINE`/`FEATURE`
   keep μ1; `SCISSION` stays in the dedicated `k_scission` moment ODEs;
   `CROSSLINK` is rejected upstream; `DISCARD`/`GAS`/`UNKNOWN` are not
   pool-scaled. The single existing classifier verdict maps cleanly to the
   single existing flag.
2. **Channel — how the verdict reaches the solver:** a declared per-reaction
   attribute `Reaction.is_end_group_reaction`, set at reaction-generation time
   from the product classification, read by the solver when filling its `int8`
   array. (Rejected alternatives: solver re-classifies at init — expensive,
   duplicates logic, needs the parent-pool map inside the solver; set passed to
   the constructor — extra arg + identity matching for no gain over an
   attribute the reaction already carries through its life.)

## Data flow

```
classify_structure (already runs in create_reacted_copy)
        │  klass == END_MOD ?
        ▼
Polymer._reacted_class            (transient, generation-time)
        │  carried via _handshake_structures into forward.products
        ▼
Reaction.is_end_group_reaction    (declared cdef public bint, default False)
        │  read in initialize_model, same order as kf/kb and ir/ip
        ▼
solver self.is_end_group_reaction[r_idx]   (int8)
        │
        ▼
existing μ0 branches fire (polymer.pyx ~915 RHS, ~1269 reaction-rate)
```

## Integration points

### 1. `rmgpy/reaction.pxd` + `rmgpy/reaction.py`

`Reaction` is a Cython `cdef class` (extension type, see `reaction.pxd:43`),
so a dynamically-set attribute is rejected — the field must be declared.

- `reaction.pxd`: add `cdef public bint is_end_group_reaction`
- `reaction.py __init__`: add parameter `is_end_group_reaction=False` and
  `self.is_end_group_reaction = is_end_group_reaction`.
- Subclasses (`TemplateReaction`) inherit the field. Requires `make` (recompile).

### 2. `rmgpy/polymer.py` — `create_reacted_copy`

`klass` is already computed at line ~821 (used today for the `CROSSLINK`
guard). Retain it on the returned polymer so the verdict survives to
reaction-build time:

```python
new_poly._reacted_class = klass
```

`Polymer` is a plain Python class with `__dict__`; `_reacted_class` is a
transient generation-time marker, not serialized.

### 3. `rmgpy/rmg/model.py` — `make_new_reaction`

After each `_handshake_structures(forward.products, polymer_reactants)` call
(both the `is_forward` and the reverse branch), the products were mutated
in place into `Polymer` objects carrying `_reacted_class`. Set the per-reaction
flag from them:

```python
forward.is_end_group_reaction = any(
    getattr(p, '_reacted_class', None) == PolymerClass.END_MOD
    for p in forward.products if isinstance(p, Polymer)
)
```

Add `PolymerClass` to the existing `from rmgpy.polymer import ...` line.

### 4. `rmgpy/solver/polymer.pyx` — `initialize_model`

Keep the zero-allocation, then fill from the reactions. Iterate
`itertools.chain(core_reactions, edge_reactions)` — the **same order** that
populates `kf/kb` and the `ir/ip` arrays (see `generate_rate_coefficients`,
line ~683) so `r_idx` matches:

```python
self.is_end_group_reaction = np.zeros(n_rxn, dtype=np.int8)
for i, rxn in enumerate(itertools.chain(core_reactions, edge_reactions)):
    if getattr(rxn, 'is_end_group_reaction', False):
        self.is_end_group_reaction[i] = 1
```

The two consumer branches already map to valid moment positions
(`mu_indices[0]` at ~916; `pool_mu0_indices`, filled at ~467, at ~1270), so no
further solver change is needed.

## Behavior & invariants

- **Only END_MOD-tagged proxy reactions change behavior.** Everything else
  defaults `False` → μ1, identical to current output. Gas–gas, pdep, and any
  reaction that never went through the polymer handshake keep μ1; the `getattr`
  default keeps the solver robust if a reaction lacks the attribute.
- **One flag per reaction, direction-independent.** "This reaction acts at a
  chain end" holds whether the proxy is consumed (reactant block ~1269) or
  produced (product block ~915); both read the same flag.
- **No double-count with scission.** `SCISSION` products never set the END_MOD
  flag; random scission remains entirely in the `k_scission` moment ODEs.

## Testing (TDD — tests written first)

1. **Solver unit** (`test/rmgpy/solver/solverPolymerTest.py`): build a
   `HybridPolymerSystem` with μ0 ≠ μ1 and two proxy reactions — one with
   `is_end_group_reaction=True`, one `False`. Assert the flagged reaction's
   rate scales with μ0 and the unflagged with μ1 (ratio check on the returned
   reaction rates / residual contribution).
2. **Classification→flag unit** (`test/rmgpy/polymerTest.py` or model-level):
   a synthetic proxy reaction whose product classifies `END_MOD` ⇒
   `reaction.is_end_group_reaction is True`; a `FEATURE` product ⇒ `False`.
3. **Regression:** `make`, then the full polymer suite
   (`polymerTest.py solverPolymerTest.py canteraTest.py
   polymerMultiPoolTest.py`) stays green (currently 184 passed).

## Cost / limitations

- One new `bint` field on the `Reaction` extension type + a recompile.
- `_reacted_class` is a private transient attribute; a reaction reconstructed
  outside the handshake path (e.g. library reactions) defaults to μ1 —
  acceptable for phase 1.
- Scope is `END_MOD` only. An explicit "initiation" reaction set (terminal
  radical creation not captured by END_MOD topology) is out of scope; revisit
  if a deck needs it.
