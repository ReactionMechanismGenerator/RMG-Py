#!/usr/bin/env python3
###############################################################################
# Tests for dynamic multi-pool polymer spawning.
#
# This is the TDD red-phase specification — every test in this module
# targets a public API that does not exist yet. As implementation lands
# (per docs/multi_pool_design.md §9), tests turn green incrementally.
#
# Companion to:
#   - docs/multi_pool_design.md (RMG-Py-side design)
#   - ~/Code/TA/PRD.md          (TA-side consumer)
###############################################################################

"""
Unit tests for dynamic multi-pool polymer spawning in :mod:`rmgpy.polymer`.

Targets the new APIs:
    * ``rmgpy.polymer.discover_repeat_motif``
    * ``rmgpy.polymer.similarity_merge``
    * ``rmgpy.polymer.MassFluxAccumulator``
    * ``rmgpy.polymer.SpawnIntent``
    * ``rmgpy.polymer.process_polymer_candidates_multipool``

All currently raise ``ImportError`` / ``AttributeError`` until §9 of the
design doc is implemented.
"""

import pytest

from rmgpy.molecule import Molecule
from rmgpy.polymer import Polymer


# ---------------------------------------------------------------------------
# Adjacency-list fixtures
# ---------------------------------------------------------------------------

# 3-mer of polyethylene: -CH2-CH2-CH2-CH2-CH2-CH2- (ethylene repeat ×3)
PE_3MER_ADJ = """
1  C u0 p0 c0 {2,S} {7,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {3,S} {10,S} {11,S}
3  C u0 p0 c0 {2,S} {4,S} {12,S} {13,S}
4  C u0 p0 c0 {3,S} {5,S} {14,S} {15,S}
5  C u0 p0 c0 {4,S} {6,S} {16,S} {17,S}
6  C u0 p0 c0 {5,S} {18,S} {19,S} {20,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
"""


@pytest.fixture
def parent_polymer():
    """A simple PE-like parent polymer pool."""
    return Polymer(
        label="PE",
        monomer="[CH2][CH2]",
        end_groups=["[H]", "[H]"],
        cutoff=3,
        Mn=1000.0,
        Mw=2500.0,
        initial_mass=1.0,
    )


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

class TestDiscoverRepeatMotif:
    """Auto-detect a repeat motif within a single molecule."""

    def test_methane_returns_none(self):
        from rmgpy.polymer import discover_repeat_motif

        mol = Molecule(smiles="C")
        assert discover_repeat_motif(mol) is None, (
            "A single-carbon molecule has no internal polymer motif"
        )

    def test_pe_3mer_returns_ch2_motif(self):
        from rmgpy.polymer import discover_repeat_motif

        mol = Molecule().from_adjacency_list(PE_3MER_ADJ)
        motif = discover_repeat_motif(mol)
        assert motif is not None, "PE 3-mer must yield a discoverable motif"
        # The motif is a Group built from heavy atoms only; one ethylene
        # repeat unit has 2 atoms.
        assert len(motif.atoms) == 2, (
            f"PE motif should be ethylene-sized; got {len(motif.atoms)} GroupAtoms"
        )

    def test_motif_count_meets_cutoff(self):
        """The discovered motif must occur at least twice in the source molecule."""
        from rmgpy.polymer import (
            count_disjoint_subgraph_isomorphisms,
            discover_repeat_motif,
        )

        mol = Molecule().from_adjacency_list(PE_3MER_ADJ)
        motif = discover_repeat_motif(mol)
        n_occ = count_disjoint_subgraph_isomorphisms(mol, motif)
        assert n_occ >= 2, "Discovered motif must repeat at least twice"


# ---------------------------------------------------------------------------
# similarity_merge
# ---------------------------------------------------------------------------

class TestSimilarityMerge:
    """A near-duplicate motif should merge into an existing pool, not spawn."""

    def test_exact_duplicate_merges(self, parent_polymer):
        from rmgpy.polymer import similarity_merge

        # Use the parent's own backbone group as the candidate motif: must merge.
        candidate = parent_polymer.backbone_group
        merged = similarity_merge(candidate, [parent_polymer])
        assert merged is parent_polymer, (
            "Exact backbone match must merge to the parent pool"
        )

    def test_distinct_motif_returns_none(self, parent_polymer):
        from rmgpy.polymer import similarity_merge

        # A phenolic motif is structurally distinct from PE.
        phenolic = Molecule(smiles="[CH2]c1c(O)c([CH2])c(C)cc1")
        merged = similarity_merge(phenolic, [parent_polymer])
        assert merged is None, (
            "A structurally distinct motif must not merge into the PE pool"
        )


# ---------------------------------------------------------------------------
# MassFluxAccumulator (rolling window)
# ---------------------------------------------------------------------------

class TestMassFluxAccumulator:
    """Trailing-N rolling window for spawn gating."""

    def test_records_and_sums_within_window(self):
        from rmgpy.polymer import MassFluxAccumulator

        acc = MassFluxAccumulator(window=3)
        acc.record("motif_a", 0.005, iteration=0)
        acc.record("motif_a", 0.008, iteration=1)
        acc.record("motif_a", 0.002, iteration=2)
        assert acc.flux("motif_a") == pytest.approx(0.015)

    def test_drops_entries_outside_window(self):
        from rmgpy.polymer import MassFluxAccumulator

        acc = MassFluxAccumulator(window=3)
        acc.record("motif_a", 0.005, iteration=0)
        acc.record("motif_a", 0.008, iteration=1)
        acc.record("motif_a", 0.001, iteration=5)
        assert acc.flux("motif_a") == pytest.approx(0.001), (
            "Iteration 5 with window=3 must evict iteration 0 and 1 entries"
        )

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


# ---------------------------------------------------------------------------
# process_polymer_candidates_multipool — end-to-end behavior
# ---------------------------------------------------------------------------

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

    def test_baseline_products_yield_no_spawn(self, parent_polymer):
        from rmgpy.polymer import process_polymer_candidates_multipool
        from rmgpy.species import Species

        # An exact PE 3-mer matches the parent's baseline classification.
        baseline = Species().from_adjacency_list(PE_3MER_ADJ)
        processed, intents = process_polymer_candidates_multipool(
            candidates=[baseline],
            reaction_model=None,
            pool_registry=[parent_polymer],
        )
        assert intents == [], (
            "Baseline products should classify into the existing pool, not spawn"
        )
        assert len(processed) == 1

    def test_unknown_against_all_pools_stays_polymer_not_gas(self, parent_polymer, monkeypatch):
        """An UNKNOWN classification (intact backbone, no clean pool match) must
        be kept as a polymer proxy, NOT leaked to the gas phase via motif
        discovery. This unifies the multipool semantics with the single-pool
        process_polymer_candidates (classification != GAS -> polymer proxy).
        """
        import rmgpy.polymer as polymer_mod
        from rmgpy.polymer import process_polymer_candidates_multipool, PolymerClass
        from rmgpy.species import Species

        cand = Species(smiles="CCCCCC")  # structure irrelevant; class is forced below
        # Force every pool comparison to return UNKNOWN (intact backbone, unmatched).
        monkeypatch.setattr(polymer_mod, "classify_structure",
                            lambda c, p, **kw: (PolymerClass.UNKNOWN, {}))

        processed, intents = process_polymer_candidates_multipool(
            candidates=[cand],
            reaction_model=None,
            pool_registry=[parent_polymer],
        )
        assert intents == [], "UNKNOWN must not spawn a new pool"
        assert cand in processed, "UNKNOWN intact-backbone product must be kept"
        assert getattr(cand, "is_polymer_proxy", False) is True, (
            "UNKNOWN must stay in the polymer phase, not be tagged gas"
        )

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

    def test_max_pools_cap_blocks_additional_spawn(self, parent_polymer):
        from rmgpy.polymer import process_polymer_candidates_multipool
        from rmgpy.species import Species

        phenolic = Species(smiles="Oc1ccc(Cc2ccc(Cc3ccc(O)cc3)cc2)cc1")
        processed, intents = process_polymer_candidates_multipool(
            candidates=[phenolic],
            reaction_model=None,
            pool_registry=[parent_polymer],
            max_pools=1,  # already at cap
        )
        assert intents == [], "max_pools=1 must block any new spawn"


# ---------------------------------------------------------------------------
# SpawnIntent dataclass shape
# ---------------------------------------------------------------------------

class _FakeReactionModel:
    """Minimal stand-in for CoreEdgeReactionModel.

    Records every make_new_species call so the test can assert that
    register_spawned_pools delegated correctly. The real model would
    proceed to install moment dummies via _register_polymer — that path
    is covered by the existing TestPolymerRegistration suite, so we
    only check the delegation here.
    """

    def __init__(self):
        self.registered: List = []

    def make_new_species(self, obj, **kwargs):
        self.registered.append(obj)
        return obj, True


class TestRegisterSpawnedPools:
    """Daughter pools landed by drain_spawn_intents must be registered with
    the reaction model so the next RMG iteration's solver sees them as
    core species and gives them their μ-dummies."""

    def test_each_new_pool_registered_with_reaction_model(self, parent_polymer):
        from rmgpy.polymer import (
            SpawnIntent,
            drain_spawn_intents,
            register_spawned_pools,
        )

        intent = SpawnIntent(
            parent_pool=parent_polymer,
            monomer=parent_polymer.backbone_group,
            end_groups=["[H]", "[H]"],
            triggering_dp=3,
            triggering_moles=1.0e-5,
        )
        new_pools = drain_spawn_intents([intent], iteration=4)
        fake = _FakeReactionModel()

        register_spawned_pools(fake, new_pools)

        # All daughter Polymers passed to the reaction model.
        assert fake.registered == new_pools

    def test_apply_spawn_intents_drains_and_registers(self, parent_polymer):
        """End-to-end hook that callers in model.py will use:
        intents go in, registered Polymers come out, and the reaction
        model has seen them. Single entry point for the iteration-boundary
        glue.
        """
        from rmgpy.polymer import SpawnIntent, apply_spawn_intents

        intent = SpawnIntent(
            parent_pool=parent_polymer,
            monomer=parent_polymer.backbone_group,
            end_groups=["[H]", "[H]"],
            triggering_dp=3,
            triggering_moles=1.0e-5,
        )
        fake = _FakeReactionModel()

        spawned = apply_spawn_intents(
            fake, [intent], iteration=5,
            existing_pools=[parent_polymer],
        )

        assert len(spawned) == 1
        assert fake.registered == spawned
        # Daughter's μ-indices land past the parent's default (0,1,2) slots.
        assert spawned[0].mu_indices == (3, 4, 5)


class TestReactionModelIntegration:
    """Wiring between the multi-pool classifier and the live
    CoreEdgeReactionModel — the hook the RMG main loop calls at the
    iteration boundary.
    """

    def test_iteration_boundary_pass_spawns_and_registers(self):
        """With a parent Polymer registered and a novel polymer-proxy
        candidate, calling _apply_multipool_spawn_pass on the reaction
        model produces a daughter Polymer in new_species_list, with
        its own _mu0/_mu1/_mu2 dummy core species.
        """
        from rmgpy.rmg.model import CoreEdgeReactionModel
        from rmgpy.species import Species

        model = CoreEdgeReactionModel()
        parent = Polymer(
            label="PS",
            monomer="[CH2][CH]c1ccccc1",
            end_groups=["[CH3]", "[H]"],
            cutoff=3,
            Mn=5000.0,
            Mw=6000.0,
            initial_mass=1.0,
        )
        model._register_polymer(parent, generate_thermo=False)
        parent.mu_indices = (0, 1, 2)

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

        # A daughter Polymer should be registered, with auto-attached dummies.
        daughter_polys = [
            s for s in model.new_species_list
            if isinstance(s, Polymer) and s is not parent
        ]
        assert len(daughter_polys) == 1, (
            f"expected one daughter Polymer, got {len(daughter_polys)}: "
            f"{[p.label for p in daughter_polys]}"
        )
        daughter = daughter_polys[0]
        dummy_labels = {f"{daughter.label}_mu0", f"{daughter.label}_mu1",
                        f"{daughter.label}_mu2"}
        present = {s.label for s in model.new_species_list}
        assert dummy_labels.issubset(present), (
            f"Daughter moment dummies missing: {dummy_labels - present}"
        )


class TestMultiPoolPipelineEndToEnd:
    """Composition test: candidates → classify → drain → register → sidecar.

    Exercises every multi-pool public API in a single flow without a real
    RMG run. The carbon-phenol full-RMG integration test (when wired into
    the main loop) is heavier and lives in solverPolymerTest.py.
    """

    def test_novel_product_spawns_and_registers_and_serializes(
        self, parent_polymer, tmp_path
    ):
        import json
        from rmgpy.polymer import (
            apply_spawn_intents,
            process_polymer_candidates_multipool,
            write_polymer_pools_sidecar,
        )
        from rmgpy.species import Species

        # The parent pool occupies state-vector slots 0..2 by convention.
        parent_polymer.mu_indices = (0, 1, 2)

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

        fake = _FakeReactionModel()
        spawned = apply_spawn_intents(
            fake, intents, iteration=7, existing_pools=[parent_polymer],
        )
        assert len(spawned) == 1
        assert fake.registered == spawned, (
            "Daughter Polymer must be handed to the reaction model"
        )

        # Sidecar reflects the full registry — parent first, daughter second.
        registry = [parent_polymer] + spawned
        path = write_polymer_pools_sidecar(
            pool_registry=registry, output_dir=str(tmp_path), iteration=7,
        )
        with open(path) as fh:
            data = json.load(fh)

        labels = [p["label"] for p in data["pools"]]
        assert labels[0] == parent_polymer.label
        assert labels[1].startswith(parent_polymer.label + "_")
        # Daughter records its parent lineage and spawn iteration.
        daughter = data["pools"][1]
        assert daughter["parent_pool"] == parent_polymer.label
        assert daughter["spawn_iteration"] == 7
        assert daughter["mu_indices"]["mu0_idx"] == 3


class TestSchulzFloryClosure:
    """Closure helper relating μ₂ to (μ₀, μ₁) for a Schulz-Flory distribution.

    Used by the inter-pool transfer-reaction moment effects (design doc §5).
    """

    def test_closed_form_on_known_distribution(self):
        """For Flory-Schulz P(n) = p^(n-1)·(1-p), ⟨n⟩ = 5, N = 1:
        analytic μ₂ = N·(1+p)/(1-p)² = 45.
        Closed-form check: μ₂ = 2·μ₁²/μ₀ − μ₁ = 2·25/1 − 5 = 45.
        """
        from rmgpy.polymer import schulz_flory_mu2

        assert schulz_flory_mu2(mu0=1.0, mu1=5.0) == pytest.approx(45.0)


class TestPolymerPoolsSidecar:
    """Sidecar JSON writer (``polymer_pools.json``) — design doc §6."""

    def test_writes_valid_schema(self, parent_polymer, tmp_path):
        import json
        from rmgpy.polymer import (
            POLYMER_POOLS_SIDECAR_SCHEMA_VERSION,
            write_polymer_pools_sidecar,
        )

        path = write_polymer_pools_sidecar(
            pool_registry=[parent_polymer],
            output_dir=str(tmp_path),
            iteration=0,
        )
        with open(path, "r", encoding="utf-8") as fh:
            data = json.load(fh)

        assert data["schema_version"] == POLYMER_POOLS_SIDECAR_SCHEMA_VERSION
        assert data["rmg_iteration"] == 0
        assert isinstance(data["pools"], list) and len(data["pools"]) == 1
        pool = data["pools"][0]
        assert pool["label"] == "PE"
        assert pool["end_groups"] == ["[H]", "[H]"]
        assert pool["cutoff"] == 3
        # Root pool has no parent.
        assert pool["parent_pool"] is None

    def test_serialises_multiple_pools(self, parent_polymer, tmp_path):
        import json
        from rmgpy.polymer import write_polymer_pools_sidecar

        # Synthesize a second pool by hand to exercise multi-pool serialisation.
        second = Polymer(
            label="PE_d1",
            monomer="[CH2][CH2]",
            end_groups=["[H]", "[H]"],
            cutoff=3,
            Mn=800.0,
            Mw=1500.0,
            initial_mass=0.5,
        )
        second.parent_pool_label = parent_polymer.label
        second.spawn_iteration = 7
        second.spawn_metadata = {"triggering_dp": 4, "mass_flux_at_spawn": 0.012}
        second.mu_indices = (41, 42, 43)

        path = write_polymer_pools_sidecar(
            pool_registry=[parent_polymer, second],
            output_dir=str(tmp_path),
            iteration=7,
        )
        data = json.loads(open(path).read())
        assert len(data["pools"]) == 2
        d1 = data["pools"][1]
        assert d1["parent_pool"] == "PE"
        assert d1["spawn_iteration"] == 7
        assert d1["mu_indices"] == {"mu0_idx": 41, "mu1_idx": 42, "mu2_idx": 43}
        assert d1["spawn_event_metadata"]["triggering_dp"] == 4


class TestDrainSpawnIntents:
    """Pure-Python drain of queued SpawnIntents into new Polymer pools.

    This is the iteration-boundary hook (design doc §4.5) running on the
    Python side. The Cython solver-reinit half follows separately and
    consumes the Polymer objects produced here.
    """

    def test_drains_one_intent_into_one_pool(self, parent_polymer):
        from rmgpy.polymer import SpawnIntent, drain_spawn_intents

        intent = SpawnIntent(
            parent_pool=parent_polymer,
            monomer=parent_polymer.backbone_group,
            end_groups=["[H]", "[H]"],
            triggering_product=None,
            triggering_dp=4,
            triggering_moles=1.42e-5,
            mass_flux_at_spawn=0.05,
        )

        new_pools = drain_spawn_intents([intent], iteration=7)

        assert len(new_pools) == 1
        spawned = new_pools[0]
        # New pool's structural metadata mirrors the intent + parent.
        assert spawned.parent_pool_label == parent_polymer.label
        assert spawned.spawn_iteration == 7
        assert spawned.end_groups_str == ["[H]", "[H]"]
        # Label is namespaced under the parent so the sidecar can show lineage.
        assert spawned.label.startswith(parent_polymer.label + "_")

    def test_stamps_spawn_event_metadata_for_sidecar(self, parent_polymer):
        """``write_polymer_pools_sidecar`` reads ``pool.spawn_metadata``;
        drain must populate it from the intent so the JSON record is complete.
        """
        from rmgpy.polymer import SpawnIntent, drain_spawn_intents

        intent = SpawnIntent(
            parent_pool=parent_polymer,
            monomer=parent_polymer.backbone_group,
            end_groups=["[H]", "[H]"],
            triggering_dp=4,
            triggering_moles=1.42e-5,
            mass_flux_at_spawn=0.07,
        )
        spawned = drain_spawn_intents([intent], iteration=11)[0]

        meta = spawned.spawn_metadata
        assert meta["triggering_dp"] == 4
        assert meta["triggering_moles"] == pytest.approx(1.42e-5)
        assert meta["mass_flux_at_spawn"] == pytest.approx(0.07)

    def test_allocates_contiguous_mu_indices_after_existing(self, parent_polymer):
        """Daughter pools must receive μ-indices that don't collide with the
        indices already held by live pools. The Cython solver uses these
        as offsets into the resized state vector.
        """
        from rmgpy.polymer import SpawnIntent, drain_spawn_intents

        # Parent occupies the first three slots.
        parent_polymer.mu_indices = (0, 1, 2)
        intent = SpawnIntent(
            parent_pool=parent_polymer,
            monomer=parent_polymer.backbone_group,
            end_groups=["[H]", "[H]"],
            triggering_dp=3,
            triggering_moles=1.0e-5,
        )

        spawned = drain_spawn_intents(
            [intent], iteration=0, existing_pools=[parent_polymer]
        )[0]

        assert spawned.mu_indices == (3, 4, 5)

    def test_initialises_daughter_moments_from_event(self, parent_polymer):
        """B.μ_k = N·DP^k from the triggering event (design doc §5)."""
        from rmgpy.polymer import SpawnIntent, drain_spawn_intents

        N = 1.42e-5
        DP = 4
        intent = SpawnIntent(
            parent_pool=parent_polymer,
            monomer=parent_polymer.backbone_group,
            end_groups=["[H]", "[H]"],
            triggering_dp=DP,
            triggering_moles=N,
        )
        spawned = drain_spawn_intents([intent], iteration=0)[0]

        assert spawned.moments[0] == pytest.approx(N)
        assert spawned.moments[1] == pytest.approx(N * DP)
        assert spawned.moments[2] == pytest.approx(N * DP * DP)

    def test_subsequent_calls_avoid_label_collisions(self, parent_polymer):
        """Two drain calls for the same parent must produce distinct labels.

        The pool registry grows across iterations; without an
        ``existing_pools`` hint the second call would reuse ``<parent>_d1``
        and clash with a live pool.
        """
        from rmgpy.polymer import SpawnIntent, drain_spawn_intents

        intent = SpawnIntent(
            parent_pool=parent_polymer,
            monomer=parent_polymer.backbone_group,
            end_groups=["[H]", "[H]"],
            triggering_dp=3,
            triggering_moles=1.0e-5,
        )

        first = drain_spawn_intents([intent], iteration=7,
                                    existing_pools=[parent_polymer])
        second = drain_spawn_intents([intent], iteration=8,
                                     existing_pools=[parent_polymer] + first)
        assert first[0].label != second[0].label, (
            "Cross-call drains must avoid label collisions on the same parent"
        )


class TestSpawnIntentShape:
    """The SpawnIntent dataclass must carry the contract fields TA expects."""

    def test_required_fields_present(self, parent_polymer):
        from rmgpy.polymer import SpawnIntent

        intent = SpawnIntent(
            parent_pool=parent_polymer,
            monomer=parent_polymer.backbone_group,
            end_groups=["[H]", "[H]"],
            triggering_product=None,
            triggering_dp=4,
            triggering_moles=1.42e-5,
        )
        # Mirrors the polymer_pools.json schema (see design doc §6).
        assert intent.parent_pool is parent_polymer
        assert intent.triggering_dp == 4
        assert intent.triggering_moles == pytest.approx(1.42e-5)
        assert intent.end_groups == ["[H]", "[H]"]
