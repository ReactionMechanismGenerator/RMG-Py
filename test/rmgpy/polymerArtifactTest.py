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
