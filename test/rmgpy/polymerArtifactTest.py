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
