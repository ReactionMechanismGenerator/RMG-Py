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
        # Heuristic: the motif's heavy-atom count should equal one ethylene unit (2 C).
        heavy_atoms = [a for a in motif.atoms if a.element.symbol != "H"]
        assert len(heavy_atoms) == 2, (
            f"PE motif should be ethylene-sized; got {len(heavy_atoms)} heavy atoms"
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


# ---------------------------------------------------------------------------
# process_polymer_candidates_multipool — end-to-end behavior
# ---------------------------------------------------------------------------

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
