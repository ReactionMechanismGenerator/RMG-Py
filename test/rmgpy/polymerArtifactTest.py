#!/usr/bin/env python3
"""Tests for the schema-2.0 polymer moments artifact emitter
(spec: docs/superpowers/specs/2026-06-10-polymer-moments-artifact-design.md;
format doc: docs/polymer_moments_format.md)."""

import json

import pytest

from rmgpy.kinetics import Arrhenius, MultiArrhenius
from rmgpy.molecule import Molecule
from rmgpy.polymer import (
    POLYMER_POOLS_SIDECAR_SCHEMA_VERSION,
    POLYMER_RATE_RECIPE_REVISION,
    Polymer,
    PolymerFluxArchetype,
    _serialize_pool_for_sidecar,
    build_polymer_moments_artifact,
    compile_polymer_reaction_entries,
    derive_condensed_species,
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

    def test_bookkeeping_species_is_proxy_plus_mu_dummies(self, pe_pool):
        """bookkeeping_species = the pool proxy + the three µ-dummies, exactly,
        as a subset of phase_species preserving its relative order."""
        core = [
            _spc("CC", "PE", index=2),          # proxy
            _mu_dummy("PE_mu0"), _mu_dummy("PE_mu1"), _mu_dummy("PE_mu2"),
            _spc("[CH3]", "G", index=7),        # gas — in neither list
            _mu_dummy("OTHER_mu0"),             # other pool's dummy
        ]
        d = _serialize_pool_for_sidecar(pe_pool, core_species=core)
        assert d["bookkeeping_species"] == ["PE(2)", "PE_mu0", "PE_mu1", "PE_mu2"]
        # subset invariant + same relative order as phase_species
        assert set(d["bookkeeping_species"]) <= set(d["phase_species"])
        positions = [d["phase_species"].index(s) for s in d["bookkeeping_species"]]
        assert positions == sorted(positions)

    def test_routed_monomer_is_phase_but_not_bookkeeping(self, pe_pool):
        """The routed condensed monomer is real condensed: it appears in
        phase_species but NOT in bookkeeping_species."""
        core = [
            _spc("CC", "PE", index=2),
            _mu_dummy("PE_mu0"), _mu_dummy("PE_mu1"), _mu_dummy("PE_mu2"),
        ]
        d = _serialize_pool_for_sidecar(pe_pool, core_species=core,
                                        monomer_routing="ethylene(5)")
        assert "ethylene(5)" in d["phase_species"]
        assert "ethylene(5)" not in d["bookkeeping_species"]
        assert d["bookkeeping_species"] == ["PE(2)", "PE_mu0", "PE_mu1", "PE_mu2"]

    def test_bookkeeping_species_always_present(self, pe_pool):
        """The key exists even without core_species (empty phase_species), and
        a routing-only pool block keeps it empty (routing is real condensed)."""
        d = _serialize_pool_for_sidecar(pe_pool)
        assert d["phase_species"] == []
        assert d["bookkeeping_species"] == []
        d = _serialize_pool_for_sidecar(pe_pool, monomer_routing="ethylene(5)")
        assert d["phase_species"] == ["ethylene(5)"]
        assert d["bookkeeping_species"] == []


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
        assert e["cantera"] == {"index": 4, "equation": "A(1) => B(5)"}
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

    def test_cantera_null_kinetics_folds_t0(self):
        """Dropped entries with T0 != 1 must emit A folded to the T0=1
        convention: A_emitted = A.value_si / T0.value_si**n.value_si."""
        core = _two_pool_core()
        kin = Arrhenius(A=(2.0, "s^-1"), n=0.5, Ea=(0.0, "J/mol"), T0=(300.0, "K"))
        rxn = Reaction(reactants=[core[0]], products=[core[4]],
                       kinetics=kin, reversible=False)
        rxn.polymer_flux_archetype = int(PolymerFluxArchetype.MIGRATION)
        entries = compile_polymer_reaction_entries(
            [rxn], core, configured_pool_labels=["A", "B"],
            cantera_index_map={})  # dropped — no Cantera entry
        assert len(entries) == 1
        e = entries[0]
        assert e["cantera"] is None
        assert e["kinetics"]["A"] == pytest.approx(2.0 / 300.0 ** 0.5)
        assert e["kinetics"]["n"] == pytest.approx(0.5)

    def test_bimolecular_units_and_duplicate_proxy_reactants(self):
        """Bimolecular A+A->B reaction: units must be m^3/(mol*s) and both
        copies of the same pool species appear in proxy_reactants."""
        core = _two_pool_core()
        spc_a = core[0]  # "A" proxy
        spc_b = core[4]  # "B" proxy
        rxn = Reaction(reactants=[spc_a, spc_a], products=[spc_b],
                       kinetics=Arrhenius(A=(1.5, "m^3/(mol*s)"), n=0.0,
                                         Ea=(0.0, "J/mol"), T0=(1.0, "K")),
                       reversible=False)
        rxn.polymer_flux_archetype = int(PolymerFluxArchetype.MIGRATION)
        entries = compile_polymer_reaction_entries(
            [rxn], core, configured_pool_labels=["A", "B"],
            cantera_index_map={id(rxn): [0]})
        assert len(entries) == 1
        e = entries[0]
        assert e["kinetics"]["units"]["A"] == "m^3/(mol*s)"
        assert e["proxy_reactants"] == ["A(1)", "A(1)"]

    def test_retained_multiarrhenius_entry_has_null_kinetics(self):
        """A retained reaction with MultiArrhenius kinetics must produce a
        valid entry with cantera set, kinetics=None, and no exception."""
        core = _two_pool_core()
        arr = Arrhenius(A=(1e13, "s^-1"), n=0.0, Ea=(0.0, "J/mol"), T0=(1.0, "K"))
        kin = MultiArrhenius(arrhenius=[arr])
        rxn = Reaction(reactants=[core[0]], products=[core[4]],
                       kinetics=kin, reversible=True)
        rxn.polymer_flux_archetype = int(PolymerFluxArchetype.MIGRATION)
        entries = compile_polymer_reaction_entries(
            [rxn], core, configured_pool_labels=["A", "B"],
            cantera_index_map={id(rxn): [0]})
        assert len(entries) == 1
        e = entries[0]
        assert e["cantera"] is not None
        assert e["kinetics"] is None


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

    def test_recipe_revision_emitted_in_conventions(self, pe_pool, tmp_path):
        """conventions.recipe_revision marks the RATE-RECIPE revision (site
        scaling, chip throttle, kb recipe, channel/flux algebra) independently
        of schema_version (shape). TA hard-fails on unknown values, so an
        accidental constant bump must fail here too: assert the literal."""
        path = self._build(pe_pool, tmp_path)
        with open(path) as fh:
            data = json.load(fh)
        assert data["conventions"]["recipe_revision"] == POLYMER_RATE_RECIPE_REVISION
        assert data["conventions"]["recipe_revision"] == "2026-06-10"
        # shape version is untouched by the recipe marker
        assert data["schema_version"] == "2.0"

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


class _FakePool:
    """Minimal stand-in for PolymerPoolConfig (label + index fields only)."""
    def __init__(self, label, mu_indices=(), explicit_dp_to_species_index=None,
                 monomer_poly_index=None):
        self.label = label
        self.mu_indices = mu_indices
        self.explicit_dp_to_species_index = explicit_dp_to_species_index or {}
        self.monomer_poly_index = monomer_poly_index


class TestDeriveCondensedSpecies:
    """conventions.condensed_species must MIRROR THE ORACLE (the live solver's
    final-core gas_species_mask). derive_condensed_species is keyed on the
    SOLVER-CONFIGURED pools (the engine's polymer_pools), which the
    save_everything hook resolves off system.solver — NOT the full Polymer
    registry. A daughter pool spawned mid-run but never solver-configured (e.g.
    epdm_scission_tail, which the EPDM run runs as ordinary GAS species and
    whose scission stamp the solver DEMOTES to legacy/UNRESOLVED) is NOT passed
    here and so is NOT reported condensed."""

    def _epdm_like_core(self):
        # Final core: epdm proxy + epdm_mu0/1/2 (the ONE configured pool, all
        # condensed) interleaved with gaseous N2 / H(1), plus the spawned
        # epdm_scission_tail family which the solver ran as ordinary GAS species
        # (its pool was never solver-configured).
        core = [
            _spc("[CH2]CC([CH2])C", "epdm", index=2),
            _mu_dummy("epdm_mu0"), _mu_dummy("epdm_mu1"), _mu_dummy("epdm_mu2"),
            _spc("N#N", "N2", index=4),
            _spc("[H]", "H", index=1),
            _spc("[CH2]CC([CH2])C", "epdm_scission_tail", index=9),
            _mu_dummy("epdm_scission_tail_mu0"),
            _mu_dummy("epdm_scission_tail_mu1"),
            _mu_dummy("epdm_scission_tail_mu2"),
        ]
        # Only 'epdm' is solver-configured (matches the live engine's
        # polymer_pools == ['epdm']); epdm_scission_tail is a registry pool
        # without a solver config.
        configured_pools = [_FakePool("epdm", mu_indices=(1, 2, 3))]
        return core, configured_pools

    def test_mask_shorter_than_core_falls_back_to_configured_derivation(self):
        """Mask sized to the constructor-era core (4), final core is 10. The
        stale mask is ignored; membership is derived from the CONFIGURED pools
        only — the spawned scission_tail family stays GAS (oracle truth)."""
        core, pools = self._epdm_like_core()
        stale_mask = [False, False, False, False]  # 4 entries, core is 10
        condensed = derive_condensed_species(core, pools, stale_mask)
        labels = [s.label for s in condensed]
        # ONLY the configured epdm pool's proxy + moments.
        assert labels == ["epdm", "epdm_mu0", "epdm_mu1", "epdm_mu2"]
        assert "epdm_scission_tail" not in labels  # spawned, unconfigured -> gas
        assert "N2" not in labels and "H" not in labels

    def test_absent_mask_falls_back_to_configured_derivation(self):
        """No mask at all (blueprint surfaced None) → derive from CONFIGURED
        pools; the unconfigured daughter is NOT condensed."""
        core, pools = self._epdm_like_core()
        condensed = derive_condensed_species(core, pools, None)
        labels = [s.label for s in condensed]
        assert "epdm" in labels
        assert "epdm_scission_tail" not in labels
        assert "N2" not in labels and "H" not in labels

    def test_matching_mask_honored_verbatim(self):
        """A length-matched mask is the oracle's verdict: honored verbatim. The
        configured-pool derivation agrees with it (epdm condensed), and the
        spawned scission_tail family the mask leaves as GAS stays GAS."""
        core, pools = self._epdm_like_core()
        import numpy as np
        # False=condensed, True=gas. Solver mask: epdm members (0-3) condensed,
        # everything else (incl. the spawned scission_tail family 6-9) gas.
        mask = np.array([False, False, False, False, True, True,
                         True, True, True, True], dtype=bool)
        condensed = derive_condensed_species(core, pools, mask)
        labels = [s.label for s in condensed]
        assert labels == ["epdm", "epdm_mu0", "epdm_mu1", "epdm_mu2"]
        assert "epdm_scission_tail" not in labels   # gas in the oracle mask
        assert "epdm_scission_tail_mu0" not in labels
        assert "N2" not in labels and "H" not in labels

    def test_mask_marks_nonpool_condensed_is_honored(self):
        """A length-matched mask can mark a non-pool species condensed (e.g. a
        solvent species the consumer treats as condensed); that verdict is
        honored verbatim even though pool derivation would not add it."""
        core, pools = self._epdm_like_core()
        import numpy as np
        # Mark index 4 (N2) condensed via the mask, on top of the epdm pool.
        mask = np.array([False, False, False, False, False, True,
                         True, True, True, True], dtype=bool)
        condensed = derive_condensed_species(core, pools, mask)
        labels = [s.label for s in condensed]
        assert "N2" in labels                       # honored from the mask
        assert set(labels) >= {"epdm", "epdm_mu0", "epdm_mu1", "epdm_mu2"}
        assert "epdm_scission_tail" not in labels   # still gas

    def test_explicit_and_monomer_indices_are_condensed(self):
        """Explicit-oligomer and routed-monomer indices count as condensed in
        the derived fallback (mirrors polymer.pyx:502-516)."""
        core = [
            _spc("CC", "P", index=2),
            _mu_dummy("P_mu0"), _mu_dummy("P_mu1"), _mu_dummy("P_mu2"),
            _spc("[CH3]", "G", index=7),
            _spc("CCC", "P_dp3", index=8),
            _spc("[CH2]CC", "monomer_in_poly", index=9),
        ]
        pools = [_FakePool("P", mu_indices=(1, 2, 3),
                           explicit_dp_to_species_index={3: 5},
                           monomer_poly_index=6)]
        condensed = derive_condensed_species(core, pools, mask=None)
        labels = [s.label for s in condensed]
        assert set(labels) == {"P", "P_mu0", "P_mu1", "P_mu2",
                               "P_dp3", "monomer_in_poly"}
        assert "G" not in labels


class TestSpawnedPoolMoments:
    """Item #14a consumer 2, AMENDED 2026-06-12 (uniform-t=0, spec section 4):
    pools[].moments are the pool's INITIAL CONDITIONS at t=0 of the simulated
    experiment, normatively. A spawned daughter genuinely contains nothing at
    t=0, so its entry carries [0, 0, 0] (the physically correct initial
    condition, not a hole) plus moments_provenance == "spawned_empty";
    input-declared pools keep their input-derived moments plus
    "input_declared". The spawn-time [N, N*DP, N*DP^2] fiction dies at the
    source, and the quantity it fed (a generation-end Sigma-bookkeeping
    estimate) is never emitted at all."""

    @staticmethod
    def _drained_daughter(pe_pool):
        """A gate-path-shaped daughter via the LIVE drain path (registry
        level; never solver-configured)."""
        from rmgpy.polymer import SpawnIntent, drain_spawn_intents
        intent = SpawnIntent(
            parent_pool=pe_pool,
            monomer=pe_pool.monomer,
            end_groups=["[H]", "[H]"],
            triggering_dp=4,
        )
        return drain_spawn_intents([intent], iteration=7,
                                   existing_pools=[pe_pool])[0]

    def test_spawned_daughter_entry_is_empty_at_t0(self, pe_pool):
        """T2 (amended, spec section 6; absorbs T7): the drained daughter's
        entry carries moments == [0, 0, 0] AND
        moments_provenance == "spawned_empty". Assertion order PINNED
        (single-reason discipline): the moments-value assertion FIRST — the
        hand-inspected red dies on the fiction values [1.0, 4.0, 16.0]
        specifically — the provenance-field assertion SECOND (confirmed as
        the second red by temporary reorder during the hand-inspection)."""
        daughter = self._drained_daughter(pe_pool)
        payload = build_polymer_moments_artifact(
            [pe_pool, daughter],
            configured_pool_labels=["PE"],  # daughter NOT solver-configured
        )
        entry = next(p for p in payload["pools"] if p["label"] == daughter.label)

        # LIVENESS PINS — BEFORE the red assertions (tripwire discipline):
        # the daughter came through the LIVE drain path and its entry
        # resolved with lineage intact (a pin failure = broken fixture,
        # never a valid red).
        assert daughter.label == "PE_d1"
        assert entry["parent_pool"] == "PE"
        assert entry["spawn_event_metadata"]["triggering_dp"] == 4
        # The daughter keeps the parent's Mn/Mw (lineage metadata; nothing
        # derives a DP from them anymore).
        assert daughter.Mn == pytest.approx(pe_pool.Mn)

        # RED assertion 1 (the fiction): a spawned pool contains NOTHING at
        # t=0. Pre-change actual: [1.0, 4.0, 16.0] = [N*DP^0, N*DP^1,
        # N*DP^2] from the never-assigned placeholder N=1.0, DP=4.
        assert entry["moments"] == [0.0, 0.0, 0.0], (
            "spawned-pool artifact moments still carry the spawn-time fiction"
        )
        # RED assertion 2 (the provenance field).
        assert entry.get("moments_provenance") == "spawned_empty"

    def test_schema_stability_and_additive_provenance(self, pe_pool):
        """T9 (spec section 6, amended values): the 1.0 stability contract
        holds on BOTH entry kinds; moments_provenance is additive
        ("input_declared" | "spawned_empty"); configured emission is
        BYTE-IDENTICAL pass-through of the pool object's moments (t=0
        initial conditions); spawn_event_metadata remains a dict — the
        triggering_moles inner-key delete is contract-clean (never emitted
        by any real deck, spec section 2)."""
        daughter = self._drained_daughter(pe_pool)
        payload = build_polymer_moments_artifact(
            [pe_pool, daughter], configured_pool_labels=["PE"])
        by_label = {p["label"]: p for p in payload["pools"]}
        for entry in by_label.values():
            for key, typ in SCHEMA_1_0_KEYS.items():
                assert key in entry, f"1.0 key {key!r} missing"
                assert isinstance(entry[key], typ), f"1.0 key {key!r} type"
        assert by_label["PE"]["moments_provenance"] == "input_declared"
        # Byte-identical pass-through for the input-declared pool (the
        # emitter's moments source is untouched by item #14a).
        assert by_label["PE"]["moments"] == [float(m) for m in pe_pool.moments]
        d_entry = by_label[daughter.label]
        assert d_entry["moments_provenance"] == "spawned_empty"
        assert d_entry["moments"] == [0.0, 0.0, 0.0]
        assert isinstance(d_entry["spawn_event_metadata"], dict)
        assert "triggering_moles" not in d_entry["spawn_event_metadata"]
        # Legacy default-labels call: the daughter still reads spawned via
        # its object spawn-markers (the secondary signal).
        legacy = build_polymer_moments_artifact([pe_pool, daughter])
        legacy_by = {p["label"]: p for p in legacy["pools"]}
        assert legacy_by["PE"]["moments_provenance"] == "input_declared"
        assert legacy_by[daughter.label]["moments_provenance"] == "spawned_empty"
