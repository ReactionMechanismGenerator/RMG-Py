#!/usr/bin/env python3
"""Tests for the polymer moments CLI reference runner (design spec §8)."""

import csv
import json
import logging
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
    _restamp_and_extend,
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


QSSA_RAW_CFG = {
    "initiation": {"A": 1.0e13, "n": 0.0, "Ea": 3.0e5},
    "depropagation": {"A": 1.0e14, "n": 0.5, "Ea": 9.0e4},
    "termination": {"A": 1.0e8, "n": 0.0, "Ea": 1.0e4},
    "transfer": {"A": 2.0e3, "n": 0.0, "Ea": 5.0e4},
    "efficiency": 0.8,
    "monomer_yield": 0.9,
}


@pytest.fixture
def qssa_deck(tmp_path):
    """A radical-QSSA pool (k_unzip == 0) with wired monomer routing, written
    out as chem.yaml + polymer_pools.json exactly like an RMG run would."""
    n2 = _spc("N#N", "N2", index=1)
    sty = _spc("C=Cc1ccccc1", "styrene", index=2)
    mus = [_mu("PS_mu0"), _mu("PS_mu1"), _mu("PS_mu2")]
    core = [n2, sty] + mus
    data, index_map = generate_cantera_data(core, [],
                                            return_reaction_index_map=True)
    chem_path = os.path.join(str(tmp_path), "chem.yaml")
    with open(chem_path, "w") as fh:
        yaml.dump(data, fh, sort_keys=False, default_flow_style=None)

    pool = Polymer(label="PS", monomer="[CH2][CH](c1ccccc1)",
                   end_groups=["[H]", "[H]"], cutoff=3,
                   moments=[1.0, 50.0, 3000.0], initial_mass=0.0,
                   k_scission=0.0, k_unzip=0.0,
                   radical_qssa_unzip=dict(QSSA_RAW_CFG))
    artifact = build_polymer_moments_artifact(
        [pool], core_species=core, core_reactions=[],
        configured_pool_labels=["PS"], condensed_species=mus + [sty],
        monomer_routing_by_pool={"PS": "styrene(2)"},
        cantera_index_map=index_map)
    art_path = os.path.join(str(tmp_path), "polymer_pools.json")
    with open(art_path, "w") as fh:
        json.dump(artifact, fh, indent=2, default=str)
    return chem_path, art_path


def _build_qssa(qssa_deck, artifact=None):
    chem_path, art_path = qssa_deck
    if artifact is None:
        with open(art_path) as fh:
            artifact = json.load(fh)
    species, reactions = load_chem_yaml(chem_path)
    return build_system_from_artifact(
        artifact, species, reactions, T0=650.0, P=1.0e5, V_poly=1.0,
        initial_moles={"N2(1)": 1.0}, mass_transfer_spec=[])


def _load_artifact(qssa_deck):
    with open(qssa_deck[1]) as fh:
        return json.load(fh)


class TestRadicalQssaArtifactLoader:
    """build_system_from_artifact must parse + validate the sidecar's
    radical_qssa_unzip block (milestone 3): shared-validator field rules,
    pinned units (ERROR on mismatch, never convert), enabled-without-routing
    hard error, and the M1 flattening picking up the parsed config."""

    def test_round_trip_flattened_arrays_match_input_constants(self, qssa_deck):
        rs, core, _ = _build_qssa(qssa_deck)
        assert rs.qssa_enabled[0] == 1
        assert rs.qssa_ki_A[0] == pytest.approx(1.0e13)
        assert rs.qssa_ki_n[0] == pytest.approx(0.0)
        assert rs.qssa_ki_Ea[0] == pytest.approx(3.0e5)
        assert rs.qssa_kdp_A[0] == pytest.approx(1.0e14)
        assert rs.qssa_kdp_n[0] == pytest.approx(0.5)
        assert rs.qssa_kdp_Ea[0] == pytest.approx(9.0e4)
        assert rs.qssa_kt_A[0] == pytest.approx(1.0e8)
        assert rs.qssa_kt_n[0] == pytest.approx(0.0)
        assert rs.qssa_kt_Ea[0] == pytest.approx(1.0e4)
        assert rs.qssa_has_transfer[0] == 1
        assert rs.qssa_ktr_A[0] == pytest.approx(2.0e3)
        assert rs.qssa_ktr_Ea[0] == pytest.approx(5.0e4)
        assert rs.qssa_efficiency[0] == pytest.approx(0.8)
        assert rs.qssa_monomer_yield[0] == pytest.approx(0.9)
        # routing survived: the pool config's monomer index points at styrene
        pool_cfg = rs.polymer_pools[0]
        assert pool_cfg.monomer_poly_index is not None
        assert core[pool_cfg.monomer_poly_index].label == "styrene(2)"
        assert pool_cfg.k_unzip == 0.0

    def test_legacy_sidecar_without_block_loads_unchanged(self, deck):
        """A pre-QSSA sidecar (no radical_qssa_unzip key) must load exactly
        as before, with the channel disabled in the flattened state."""
        chem_path, art_path = deck
        with open(art_path) as fh:
            artifact = json.load(fh)
        assert "radical_qssa_unzip" not in artifact["pools"][0]["channels"]
        species, reactions = load_chem_yaml(chem_path)
        rs, _, _ = build_system_from_artifact(
            artifact, species, reactions, T0=800.0, P=1.0e5, V_poly=1.0,
            initial_moles={"N2(1)": 1.0}, mass_transfer_spec=[])
        assert rs.qssa_enabled[0] == 0
        assert rs.polymer_pools[0].radical_qssa_unzip is None

    def test_rejects_disabled_block_as_malformed(self, qssa_deck):
        """The serializer never emits enabled=false: a disabled channel is
        ABSENT from the sidecar. A present-disabled block would carry
        constants nobody validates (enabled=false + garbage used to load
        silently), so it is rejected outright (round-23)."""
        artifact = _load_artifact(qssa_deck)
        artifact["pools"][0]["channels"]["radical_qssa_unzip"]["enabled"] = False
        with pytest.raises(ValueError, match=r"PS.*present-disabled"):
            _build_qssa(qssa_deck, artifact)

    def test_rejects_enabled_without_routing(self, qssa_deck):
        artifact = _load_artifact(qssa_deck)
        artifact["pools"][0]["monomer_routing"] = None
        with pytest.raises(ValueError,
                           match=r"PS.*radical_qssa_unzip.*monomer_routing"):
            _build_qssa(qssa_deck, artifact)

    def test_rejects_wrong_units_string_no_conversion(self, qssa_deck):
        """A sidecar claiming a different unit system must ERROR, not
        convert: termination pinned bimolecular m^3/(mol*s)."""
        artifact = _load_artifact(qssa_deck)
        block = artifact["pools"][0]["channels"]["radical_qssa_unzip"]
        block["termination"]["units"]["A"] = "s^-1"
        with pytest.raises(ValueError, match=r"PS.*termination.*units"):
            _build_qssa(qssa_deck, artifact)

    def test_rejects_bimolecular_transfer_units(self, qssa_deck):
        """ktr is PSEUDO-first-order [s^-1]; a sidecar claiming bimolecular
        transfer units must ERROR (premultiplication happened upstream or
        not at all -- never converted silently here)."""
        artifact = _load_artifact(qssa_deck)
        block = artifact["pools"][0]["channels"]["radical_qssa_unzip"]
        block["transfer"]["units"]["A"] = "m^3/(mol*s)"
        with pytest.raises(ValueError, match=r"PS.*transfer.*units"):
            _build_qssa(qssa_deck, artifact)

    def test_rejects_missing_units_block(self, qssa_deck):
        artifact = _load_artifact(qssa_deck)
        block = artifact["pools"][0]["channels"]["radical_qssa_unzip"]
        del block["initiation"]["units"]
        with pytest.raises(ValueError, match=r"PS.*initiation.*units"):
            _build_qssa(qssa_deck, artifact)

    def test_rejects_unknown_basis(self, qssa_deck):
        artifact = _load_artifact(qssa_deck)
        block = artifact["pools"][0]["channels"]["radical_qssa_unzip"]
        block["basis"] = "chain_ends_mu0"
        with pytest.raises(ValueError, match=r"PS.*basis"):
            _build_qssa(qssa_deck, artifact)

    def test_rejects_qssa_with_nonzero_unzip(self, qssa_deck):
        """Mutual exclusion mirrors the generation-side guard: both channel
        representations on one pool would double-count the unzip flux."""
        artifact = _load_artifact(qssa_deck)
        artifact["pools"][0]["channels"]["unzip"]["A"] = 0.5
        with pytest.raises(ValueError,
                           match=r"PS.*radical_qssa_unzip.*k_unzip"):
            _build_qssa(qssa_deck, artifact)

    def test_rejects_missing_or_non_bool_enabled(self, qssa_deck):
        artifact = _load_artifact(qssa_deck)
        block = artifact["pools"][0]["channels"]["radical_qssa_unzip"]
        del block["enabled"]
        with pytest.raises(ValueError, match=r"PS.*enabled"):
            _build_qssa(qssa_deck, artifact)
        artifact = _load_artifact(qssa_deck)
        artifact["pools"][0]["channels"]["radical_qssa_unzip"]["enabled"] = 1
        with pytest.raises(ValueError, match=r"PS.*enabled"):
            _build_qssa(qssa_deck, artifact)

    def test_rejects_unknown_key_in_block(self, qssa_deck):
        artifact = _load_artifact(qssa_deck)
        block = artifact["pools"][0]["channels"]["radical_qssa_unzip"]
        block["efficency"] = 0.5  # typo'd efficiency must not be dropped
        with pytest.raises(ValueError, match=r"PS.*efficency"):
            _build_qssa(qssa_deck, artifact)


class TestQssaSchemaVersionGate:
    """Version contract (round-23): the QSSA channel vocabulary was
    introduced in schema 2.1 (channel-vocabulary growth = minor bump).
    Acceptance matrix: 2.0 no-QSSA -> load; 2.0 with-QSSA -> reject
    (malformed); 2.1 either way -> load; future 2.x minor -> load (the
    envelope gate has always been minor-permissive: main() pins
    startswith('2.')); non-2.x major carrying QSSA -> reject."""

    def test_qssa_artifact_is_stamped_2_1_and_accepted(self, qssa_deck):
        artifact = _load_artifact(qssa_deck)
        assert artifact["schema_version"] == "2.1"
        rs, _, _ = _build_qssa(qssa_deck, artifact)
        assert rs.qssa_enabled[0] == 1

    def test_legacy_deck_still_stamped_2_0_and_accepted(self, deck):
        chem_path, art_path = deck
        with open(art_path) as fh:
            artifact = json.load(fh)
        assert artifact["schema_version"] == "2.0"
        species, reactions = load_chem_yaml(chem_path)
        rs, _, _ = build_system_from_artifact(
            artifact, species, reactions, T0=800.0, P=1.0e5, V_poly=1.0,
            initial_moles={"N2(1)": 1.0}, mass_transfer_spec=[])
        assert rs.qssa_enabled[0] == 0

    def test_accepts_2_1_without_qssa(self, deck):
        chem_path, art_path = deck
        with open(art_path) as fh:
            artifact = json.load(fh)
        artifact["schema_version"] = "2.1"
        species, reactions = load_chem_yaml(chem_path)
        rs, _, _ = build_system_from_artifact(
            artifact, species, reactions, T0=800.0, P=1.0e5, V_poly=1.0,
            initial_moles={"N2(1)": 1.0}, mass_transfer_spec=[])
        assert rs.qssa_enabled[0] == 0

    def test_rejects_qssa_block_in_2_0_artifact(self, qssa_deck):
        """A 2.0 artifact containing the QSSA vocabulary is malformed --
        the emitter stamps 2.1 whenever it writes the block."""
        artifact = _load_artifact(qssa_deck)
        artifact["schema_version"] = "2.0"
        with pytest.raises(ValueError, match=r"schema_version.*2\.1"):
            _build_qssa(qssa_deck, artifact)

    def test_rejects_qssa_block_in_non_2x_artifact(self, qssa_deck):
        artifact = _load_artifact(qssa_deck)
        artifact["schema_version"] = "3.0"
        with pytest.raises(ValueError, match=r"schema_version"):
            _build_qssa(qssa_deck, artifact)

    def test_accepts_future_minor_with_qssa(self, qssa_deck):
        artifact = _load_artifact(qssa_deck)
        artifact["schema_version"] = "2.7"
        rs, _, _ = _build_qssa(qssa_deck, artifact)
        assert rs.qssa_enabled[0] == 1


class TestQssaNormativeRecipe:
    """The QSSA block's machine-readable recipe sub-block (round-23): the
    loader validates every field by EXACT match against its own pinned copy
    (same idiom as the units pin -- reject, never adapt) and rejects a QSSA
    block that lacks the recipe entirely."""

    PINNED_RECIPE = {
        "bond_basis": ("B = max(mu1 - mu0, 0) on concentration moments "
                       "(mol/m^3 condensed)"),
        "rate_no_transfer": ("r_mono = monomer_yield * kdp * "
                             "sqrt(efficiency * ki * B / kt)"),
        "rate_with_transfer": ("r_mono = monomer_yield * kdp * "
                               "(sqrt(ktr^2 + 8*kt*(2*efficiency*ki*B)) "
                               "- ktr) / (4*kt)"),
        "moment_signature": ("dmu0 = 0; dmu1 -= r_mono; dmu2 -= r_mono * "
                             "max(2*mu1/max(mu0, small_eps) - 1, 0)"),
        "small_eps": 1e-30,
        "volume_note": ("kt is bimolecular: rates depend on condensed "
                        "volume V_poly; consumers MUST evaluate on "
                        "concentration moments mu_k = n_k / V_poly and "
                        "convert emitted rate back with *V_poly"),
    }

    def test_recipe_round_trips_and_validates(self, qssa_deck):
        artifact = _load_artifact(qssa_deck)
        block = artifact["pools"][0]["channels"]["radical_qssa_unzip"]
        assert block["recipe"] == self.PINNED_RECIPE
        rs, _, _ = _build_qssa(qssa_deck, artifact)
        assert rs.qssa_enabled[0] == 1

    def test_rejects_missing_recipe(self, qssa_deck):
        """QSSA block without recipe = malformed (never emitted that way)."""
        artifact = _load_artifact(qssa_deck)
        del artifact["pools"][0]["channels"]["radical_qssa_unzip"]["recipe"]
        with pytest.raises(ValueError, match=r"PS.*recipe"):
            _build_qssa(qssa_deck, artifact)

    def test_rejects_recipe_field_mismatch_never_adapts(self, qssa_deck):
        artifact = _load_artifact(qssa_deck)
        recipe = artifact["pools"][0]["channels"]["radical_qssa_unzip"]["recipe"]
        recipe["rate_no_transfer"] = "r_mono = kdp * sqrt(ki * B / kt)"
        with pytest.raises(ValueError, match=r"PS.*rate_no_transfer"):
            _build_qssa(qssa_deck, artifact)

    def test_rejects_recipe_missing_field(self, qssa_deck):
        artifact = _load_artifact(qssa_deck)
        recipe = artifact["pools"][0]["channels"]["radical_qssa_unzip"]["recipe"]
        del recipe["moment_signature"]
        with pytest.raises(ValueError, match=r"PS.*moment_signature"):
            _build_qssa(qssa_deck, artifact)

    def test_rejects_recipe_small_eps_drift(self, qssa_deck):
        artifact = _load_artifact(qssa_deck)
        recipe = artifact["pools"][0]["channels"]["radical_qssa_unzip"]["recipe"]
        recipe["small_eps"] = 1e-29
        with pytest.raises(ValueError, match=r"PS.*small_eps"):
            _build_qssa(qssa_deck, artifact)

    def test_rejects_unknown_recipe_key(self, qssa_deck):
        artifact = _load_artifact(qssa_deck)
        recipe = artifact["pools"][0]["channels"]["radical_qssa_unzip"]["recipe"]
        recipe["extra_note"] = "x"
        with pytest.raises(ValueError, match=r"PS.*extra_note"):
            _build_qssa(qssa_deck, artifact)


class TestChemYamlLoader:
    def test_load_chem_yaml(self, deck):
        chem_path, art_path = deck
        species, reactions = load_chem_yaml(chem_path)
        labels = [s.label for s in species]
        assert "N2(1)" in labels and "poly_mu0" in labels
        assert len(reactions) == 1
        rxn = reactions[0]
        assert rxn.kinetics.A.value_si == pytest.approx(5.0)
        # chem.yaml now records reversibility in the equation arrow
        # (rmgpy/cantera.py get_reaction_equation emits '=>' for irreversible
        # reactions, '<=>' otherwise), so the raw loader recovers RMG
        # irreversibility straight from the arrow — this fixture's gas_rxn is
        # reversible=False, so it round-trips reversible=False without any
        # artifact help. (_parse_equation in the runner handles '=>'.) The
        # artifact's kinetics.reversible restoration remains as belt-and-braces
        # for proxy-touching entries and now agrees with the arrow.
        assert rxn.reversible is False
        with open(art_path) as fh:
            artifact = json.load(fh)
        assert artifact["reactions"] == []  # G->C1 is ordinary chemistry
        _restamp_and_extend(artifact, species, reactions)
        assert rxn.reversible is False  # unchanged: arrow already irreversible
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


class TestUnzipRoutingGuard:
    def test_build_system_rejects_unzip_pool_without_monomer_routing(self, deck):
        """An artifact pool with k_unzip > 0 but no monomer_routing target would
        reach the solver as PolymerPoolConfig(k_unzip>0, monomer_poly_index=None)
        -- the exact shape that drains condensed moments with no gas emission
        (silently un-conserved mass). The consumer-world assembler must refuse
        it at configuration time, same as generation-side to_config."""
        chem_path, art_path = deck
        with open(art_path) as fh:
            artifact = json.load(fh)
        pool_entry = artifact["pools"][0]
        # Premise: the fixture's pool has no monomer routing target.
        assert not pool_entry.get("monomer_routing")
        pool_entry["channels"]["unzip"]["A"] = 1.0
        species, reactions = load_chem_yaml(chem_path)
        with pytest.raises(ValueError, match=r"poly.*k_unzip.*un-conserved"):
            build_system_from_artifact(
                artifact, species, reactions, T0=800.0, P=1.0e5, V_poly=1.0,
                initial_moles={"N2(1)": 1.0}, mass_transfer_spec=[])

    def test_build_system_rejects_negative_unzip_rate(self, deck):
        """A negative unzip A in the artifact is not a valid rate constant. All
        downstream k_unzip consumers are gated on k_unzip > 0, so a negative
        value would silently become an inert channel -- the assembler must
        refuse it with a clear ValueError naming the pool."""
        chem_path, art_path = deck
        with open(art_path) as fh:
            artifact = json.load(fh)
        artifact["pools"][0]["channels"]["unzip"]["A"] = -1.0
        species, reactions = load_chem_yaml(chem_path)
        with pytest.raises(ValueError,
                           match=r"poly.*k_unzip.*not a valid rate constant"):
            build_system_from_artifact(
                artifact, species, reactions, T0=800.0, P=1.0e5, V_poly=1.0,
                initial_moles={"N2(1)": 1.0}, mass_transfer_spec=[])

    @pytest.mark.parametrize("channel,field", [("unzip", "k_unzip"),
                                               ("scission", "k_scission")])
    @pytest.mark.parametrize("bad", [float("nan"), float("inf"),
                                     float("-inf")],
                             ids=["nan", "+inf", "-inf"])
    def test_build_system_rejects_non_finite_channel_rate(self, deck, channel,
                                                          field, bad):
        """NaN passes BOTH the `< 0` and `> 0` gates as False, so a
        hand-edited artifact with a non-finite unzip/scission A would make
        the channel SILENTLY INERT (NaN) or poison the residual (inf). The
        assembler must refuse it with a clear ValueError naming the pool,
        mirroring the negative-A guard."""
        chem_path, art_path = deck
        with open(art_path) as fh:
            artifact = json.load(fh)
        artifact["pools"][0]["channels"][channel]["A"] = bad
        species, reactions = load_chem_yaml(chem_path)
        with pytest.raises(ValueError,
                           match=rf"poly.*{field}.*not a valid rate constant"):
            build_system_from_artifact(
                artifact, species, reactions, T0=800.0, P=1.0e5, V_poly=1.0,
                initial_moles={"N2(1)": 1.0}, mass_transfer_spec=[])

    def test_build_system_rejects_non_string_monomer_routing(self, deck):
        """monomer_routing is a species-label string in the artifact schema
        (polymer.py _serialize_pool_for_sidecar). A hand-edited artifact with a
        non-string routing value (e.g. a dict) used to die with a raw TypeError
        from idx[routing] (unhashable) -- the assembler must instead raise an
        actionable ValueError naming the pool and the bad routing value."""
        chem_path, art_path = deck
        with open(art_path) as fh:
            artifact = json.load(fh)
        pool_entry = artifact["pools"][0]
        pool_entry["monomer_routing"] = {"target": "C1(3)"}
        pool_entry["channels"]["unzip"]["A"] = 1.0
        species, reactions = load_chem_yaml(chem_path)
        with pytest.raises(ValueError, match=r"poly.*monomer_routing"):
            build_system_from_artifact(
                artifact, species, reactions, T0=800.0, P=1.0e5, V_poly=1.0,
                initial_moles={"N2(1)": 1.0}, mass_transfer_spec=[])

    def test_build_system_rejects_unknown_monomer_routing_target(self, deck):
        """A monomer_routing label absent from the deck's species list used to
        die with a raw KeyError from idx[routing] -- the assembler must raise a
        clear ValueError naming the pool and the unresolvable target."""
        chem_path, art_path = deck
        with open(art_path) as fh:
            artifact = json.load(fh)
        pool_entry = artifact["pools"][0]
        pool_entry["monomer_routing"] = "styrene(99)"
        pool_entry["channels"]["unzip"]["A"] = 1.0
        species, reactions = load_chem_yaml(chem_path)
        with pytest.raises(ValueError, match=r"poly.*styrene\(99\)"):
            build_system_from_artifact(
                artifact, species, reactions, T0=800.0, P=1.0e5, V_poly=1.0,
                initial_moles={"N2(1)": 1.0}, mass_transfer_spec=[])


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


class TestValidation:
    """F2 + F3: input validation guards."""

    def test_run_segments_rejects_decreasing_t_end(self, deck):
        """Backward-in-time profile must raise ValueError, not silently integrate."""
        chem_path, art_path = deck
        with open(art_path) as fh:
            artifact = json.load(fh)
        species, reactions = load_chem_yaml(chem_path)
        rs, core, all_rxns = build_system_from_artifact(
            artifact, species, reactions,
            T0=800.0, P=1.0e5, V_poly=1.0,
            initial_moles={"N2(1)": 1.0}, mass_transfer_spec=[])
        with pytest.raises(ValueError, match="strictly increasing"):
            run_segments(rs, core, artifact, all_rxns,
                         [(1.0, 800.0), (0.5, 850.0)])

    def test_run_segments_rejects_duplicate_t_end(self, deck):
        """Duplicate t_end values (zero-length segment) must raise ValueError."""
        chem_path, art_path = deck
        with open(art_path) as fh:
            artifact = json.load(fh)
        species, reactions = load_chem_yaml(chem_path)
        rs, core, all_rxns = build_system_from_artifact(
            artifact, species, reactions,
            T0=800.0, P=1.0e5, V_poly=1.0,
            initial_moles={"N2(1)": 1.0}, mass_transfer_spec=[])
        with pytest.raises(ValueError, match="strictly increasing"):
            run_segments(rs, core, artifact, all_rxns,
                         [(0.5, 800.0), (0.5, 850.0)])

    def test_build_system_rejects_unknown_initial_moles_label(self, deck):
        """Unknown label in initial_moles must raise ValueError, not KeyError."""
        chem_path, art_path = deck
        with open(art_path) as fh:
            artifact = json.load(fh)
        species, reactions = load_chem_yaml(chem_path)
        with pytest.raises(ValueError, match="not found"):
            build_system_from_artifact(
                artifact, species, reactions,
                T0=800.0, P=1.0e5, V_poly=1.0,
                initial_moles={"NOSUCHSPECIES": 1.0}, mass_transfer_spec=[])


class TestReferenceStateTripwireConsumerWorld:
    """The thermo reference-state tripwire (polymer.pyx _reference_state_tripwire)
    derives per-species MW from spc.molecule[0]. Consumer-world species built by
    the runner's _species_from_yaml carry molecule == [] (label + NASA thermo
    only), so every MW comes out 0.0 and any reversible melt-touching reaction
    sends mw=0 into _sackur_tetrode_s_trans -> math.log(0) -> ValueError. The
    base ``deck`` fixture dodges this (its sole reaction is irreversible and
    all-gas), so this deck adds the missing shape: a balanced REVERSIBLE
    H-abstraction touching condensed chain species."""

    @pytest.fixture
    def reversible_melt_deck(self, tmp_path):
        """The ``deck`` fixture plus two condensed chain-scale species (PR =
        pentadecyl C15H31, PH = pentadecane C15H32) on a balanced reversible
        H-abstraction PR + CH4 <=> PH + CH3. Benign by construction: the melt
        participants sit one on each side at (near-)equal chain mass, so the
        unpaired reference-state magnitude U is ~0.003 decades — far below the
        0.5-decade census bound. A correct tripwire must stay SILENT here."""
        n2 = _spc("N#N", "N2", index=1)
        g = _spc("[CH3]", "G", index=2)
        g2 = _spc("C", "C1", index=3)
        pr = _spc("[CH2]CCCCCCCCCCCCCC", "PR")  # C15H31, condensed chain radical
        ph = _spc("CCCCCCCCCCCCCCC", "PH")      # C15H32, condensed chain parent
        mus = [_mu("poly_mu0"), _mu("poly_mu1"), _mu("poly_mu2")]
        core = [n2, g, g2, pr, ph] + mus
        gas_rxn = Reaction(reactants=[g], products=[g2],
                           kinetics=Arrhenius(A=(5.0, "1/s"), n=0.0,
                                              Ea=(10.0, "kJ/mol"), T0=(1.0, "K")),
                           reversible=False)
        melt_rxn = Reaction(reactants=[pr, g2], products=[ph, g],
                            kinetics=Arrhenius(A=(1.0e3, "m^3/(mol*s)"), n=0.0,
                                               Ea=(20.0, "kJ/mol"), T0=(1.0, "K")),
                            reversible=True)
        data, index_map = generate_cantera_data(core, [gas_rxn, melt_rxn],
                                                return_reaction_index_map=True)
        chem_path = os.path.join(str(tmp_path), "chem.yaml")
        with open(chem_path, "w") as fh:
            yaml.dump(data, fh, sort_keys=False, default_flow_style=None)

        pool = Polymer(label="poly", monomer="[CH2][CH2]",
                       end_groups=["[H]", "[H]"], cutoff=3,
                       moments=[1.0, 5.0, 30.0], initial_mass=0.0,
                       k_scission=1.0, k_unzip=0.0)
        artifact = build_polymer_moments_artifact(
            [pool], core_species=core, core_reactions=[gas_rxn, melt_rxn],
            configured_pool_labels=["poly"],
            condensed_species=mus + [pr, ph],
            cantera_index_map=index_map)
        art_path = os.path.join(str(tmp_path), "polymer_pools.json")
        with open(art_path, "w") as fh:
            json.dump(artifact, fh, indent=2, default=str)
        return chem_path, art_path

    def test_tripwire_runs_silently_on_consumer_world_species(
            self, reversible_melt_deck, caplog):
        chem_path, art_path = reversible_melt_deck
        with open(art_path) as fh:
            artifact = json.load(fh)
        species, reactions = load_chem_yaml(chem_path)

        # --- liveness pins (must hold BEFORE the act; a dead fixture cannot
        # produce this red) ---
        # 1. consumer-world species: structure does not cross the artifact
        #    boundary, so every species carries molecule == [].
        assert all(s.molecule == [] for s in species)
        # 2. the melt-touching reaction round-trips reversible=True from the
        #    '<=>' arrow.
        rev = [r for r in reactions if r.reversible]
        assert len(rev) == 1
        melt_rxn = rev[0]
        assert sorted(s.label for s in melt_rxn.reactants) == ["C1(3)", "PR"]
        assert sorted(s.label for s in melt_rxn.products) == ["G(2)", "PH"]
        # 3. both chain species are in the artifact's condensed set, which is
        #    exactly what build_system_from_artifact turns into the
        #    gas_species_mask (mask False == condensed/melt class).
        condensed = set(artifact["conventions"]["condensed_species"])
        assert {"PR", "PH"} <= condensed

        # --- act: the runner's system-build path runs initialize_model, which
        # runs _reference_state_tripwire. This deck is benign (melt chain mass
        # paired across the reaction), so the tripwire must complete WITHOUT
        # raising and stay silent — silence, not absence.
        with caplog.at_level(logging.WARNING):
            rs, core, all_rxns = build_system_from_artifact(
                artifact, species, reactions,
                T0=800.0, P=1.0e5, V_poly=1.0,
                initial_moles={"N2(1)": 1.0}, mass_transfer_spec=[])

        assert not any("THERMO REFERENCE-STATE" in rec.getMessage()
                       for rec in caplog.records)
        # the tripwire ran on this rebuild and saw the condensed species: the
        # paired chain masses keep U far below the 0.5-decade census bound
        # (REFERENCE_STATE_CENSUS_DECADES).
        idx = {s.label: i for i, s in enumerate(core)}
        assert not rs.gas_species_mask[idx["PR"]]
        assert not rs.gas_species_mask[idx["PH"]]
        assert rs.reference_state_max_decades < 0.5

    @pytest.fixture
    def entry_listed_melt_deck(self, tmp_path):
        """A deck whose artifact HAS a ``reactions[]`` entry (unlike
        ``reversible_melt_deck``, whose melt reaction has no pool participant
        and so compiles to zero entries). The EPDM-shaped reversible capture

            PR + C2H5 <=> poly        (C15H31 + C2H5 -> C17H36, balanced)

        has the POOL PROXY itself as a product, so the artifact compiler
        (rmgpy/polymer.py compile_polymer_reaction_entries) emits an entry,
        and the runner's _restamp_and_extend tag restoration + window plumbing
        are both LIVE on this deck:

        * PR (C15H31, 211.4 g/mol) is a GAS chain-scale counterparty: it is in
          the entry's reactants but NOT in proxy_reactants/proxy_products, so
          its physically-melt membership in the consumer world exists ONLY via
          the runner's blanket tag restoration. Paired against the condensed
          proxy (240.5 g/mol) it keeps U ~ 0.08 decades -> SILENT.
        * C2H5 (29.1 g/mol) sits in the bare-slack..window MW band: above the
          10 g/mol bare slack the chain window degrades to if the
          monomer_mw_g_mol plumb is deleted, below the real window
          (pool monomer C2H4 28.05 + 10 = 38.05 g/mol). With the plumb intact
          the tag-branch MW conjunct excludes it; without it, C2H5 enters the
          melt sum and U blows up to ~10 decades -> refusal."""
        n2 = _spc("N#N", "N2", index=1)
        et = _spc("C[CH2]", "C2H5")              # ethyl, 29.06 g/mol: band species
        pr = _spc("[CH2]CCCCCCCCCCCCCC", "PR")   # C15H31, GAS chain radical
        proxy = _spc("CCCCCCCCCCCCCCCCC", "poly")  # C17H36 pool-proxy species
        mus = [_mu("poly_mu0"), _mu("poly_mu1"), _mu("poly_mu2")]
        core = [n2, et, pr, proxy] + mus
        melt_rxn = Reaction(reactants=[pr, et], products=[proxy],
                            kinetics=Arrhenius(A=(1.0e3, "m^3/(mol*s)"), n=0.0,
                                               Ea=(20.0, "kJ/mol"), T0=(1.0, "K")),
                            reversible=True)
        data, index_map = generate_cantera_data(core, [melt_rxn],
                                                return_reaction_index_map=True)
        chem_path = os.path.join(str(tmp_path), "chem.yaml")
        with open(chem_path, "w") as fh:
            yaml.dump(data, fh, sort_keys=False, default_flow_style=None)

        pool = Polymer(label="poly", monomer="[CH2][CH2]",
                       end_groups=["[H]", "[H]"], cutoff=3,
                       moments=[1.0, 5.0, 30.0], initial_mass=0.0,
                       k_scission=1.0, k_unzip=0.0)
        artifact = build_polymer_moments_artifact(
            [pool], core_species=core, core_reactions=[melt_rxn],
            configured_pool_labels=["poly"],
            condensed_species=mus + [proxy],
            cantera_index_map=index_map)
        art_path = os.path.join(str(tmp_path), "polymer_pools.json")
        with open(art_path, "w") as fh:
            json.dump(artifact, fh, indent=2, default=str)
        return chem_path, art_path

    def test_entry_listed_tag_and_window_restorations_keep_tripwire_silent(
            self, entry_listed_melt_deck, caplog):
        """Pins BOTH consumer-world restorations on a live entry-listed deck:
        (1) _restamp_and_extend's is_polymer_proxy tag restoration gives the
        gas chain counterparty PR its melt membership (neuter it -> PR drops
        out of the melt sum, the condensed proxy is unpaired, U ~ 11.4 ->
        refusal); (2) the pools[].monomer_mw_g_mol -> PolymerPoolConfig plumb
        keeps the ONE chain-scale window at monomer MW + slack (delete it ->
        the window degrades to the bare 10 g/mol slack, the tagged band
        species C2H5 enters the melt sum, U ~ 10 -> refusal)."""
        chem_path, art_path = entry_listed_melt_deck
        with open(art_path) as fh:
            artifact = json.load(fh)
        species, reactions = load_chem_yaml(chem_path)
        by_label = {s.label: s for s in species}

        def mw_g_mol(label):
            # species-level quantity value_si is kg/molecule (species.py)
            return by_label[label].molecular_weight.value_si * 6.02214179e23 * 1000.0

        # --- liveness pins (must hold BEFORE the act) ---
        # 1. the melt reaction IS entry-listed (proxy product 'poly' resolves
        #    a configured pool), unlike the zero-entry reversible_melt_deck.
        assert len(artifact["reactions"]) >= 1
        entry = artifact["reactions"][0]
        assert entry["cantera"]["equation"] == "PR + C2H5 <=> poly"
        # 2. PR appears ONLY via reactants, never in the proxy_* lists, and is
        #    NOT in the condensed set: its melt membership in the consumer
        #    world exists ONLY via the runner's tag restoration.
        assert "PR" in entry["reactants"]
        assert "PR" not in entry["proxy_reactants"] + entry["proxy_products"]
        assert entry["proxy_products"] == ["poly"]
        condensed = set(artifact["conventions"]["condensed_species"])
        assert "PR" not in condensed and "poly" in condensed
        # 3. window band pins, numeric and self-checking: the REAL chain
        #    window is pool monomer MW + 10 g/mol slack (polymer.pyx
        #    REFERENCE_STATE_MW_SLACK_G_MOL); the bare-slack degraded window
        #    is 10 g/mol. PR must clear the real window (chain-scale);
        #    C2H5 must sit strictly INSIDE the bare-slack..window band.
        monomer_mw = artifact["pools"][0]["monomer_mw_g_mol"]
        assert monomer_mw == pytest.approx(28.05, abs=0.1)
        real_window = monomer_mw + 10.0
        assert mw_g_mol("PR") > real_window
        assert 10.0 < mw_g_mol("C2H5") < real_window

        # --- act: build runs initialize_model -> _reference_state_tripwire.
        # PR's restored tag pairs it against the condensed proxy (U ~ 0.08);
        # the plumbed window excludes C2H5. The build must complete SILENTLY.
        with caplog.at_level(logging.WARNING):
            rs, core, all_rxns = build_system_from_artifact(
                artifact, species, reactions,
                T0=800.0, P=1.0e5, V_poly=1.0,
                initial_moles={"N2(1)": 1.0}, mass_transfer_spec=[])

        assert not any("THERMO REFERENCE-STATE" in rec.getMessage()
                       for rec in caplog.records)
        idx = {s.label: i for i, s in enumerate(core)}
        # PR stayed GAS in the mask: its melt-class membership came from the
        # restored tag, not the condensed branch.
        assert rs.gas_species_mask[idx["PR"]]
        assert not rs.gas_species_mask[idx["poly"]]
        # the tripwire DID sum melt participants on this reaction (max U is
        # the PR-vs-proxy translational mismatch, nonzero but far below the
        # 0.5-decade census bound).
        assert 0.0 < rs.reference_state_max_decades < 0.5

    def test_species_from_yaml_populates_mw_from_composition(self):
        """The runner sources MW the way every consumer-world computation
        sources data: from artifact fields. C15H32 hand-computes to
        15*12.011 + 32*1.00794 = 212.42 g/mol."""
        from rmgpy.tools.polymer_moments_runner import _species_from_yaml
        spc = _species_from_yaml({"name": "PH", "composition": {"C": 15, "H": 32}})
        assert spc.molecule == []
        assert spc.molecular_weight is not None
        mw_g_mol = spc.molecular_weight.value_si * 6.02214179e23 * 1000.0
        assert mw_g_mol == pytest.approx(212.41, abs=0.1)

    def test_unresolvable_melt_mw_raises_contract_violation(self):
        """Defense-in-depth log-domain guard: a melt participant whose MW
        cannot be resolved (molecule == [] AND species-level
        molecular_weight unset) on a reversible reaction must raise the
        tripwire's input-contract ValueError naming the species -- never a
        bare 'math domain error', and never a silent skip (a melt chain
        with no structure is the tripwire's main case)."""
        from rmgpy.solver.polymer import HybridPolymerSystem
        rows = [NASAPolynomial(coeffs=[2.5, 0, 0, 0, 0, -745.375, 3.35532],
                               Tmin=(tmin, "K"), Tmax=(tmax, "K"))
                for tmin, tmax in ((200.0, 1000.0), (1000.0, 6000.0))]
        ghost = Species(label="GHOST")
        ghost.molecule = []  # no structure, no molecular_weight: unresolvable
        ghost.thermo = NASA(polynomials=rows, Tmin=(200.0, "K"), Tmax=(6000.0, "K"))
        gas = _spc("[CH3]", "G")
        core = [ghost, gas]
        mask = np.array([False, True], dtype=bool)  # GHOST is condensed/melt
        rxn = Reaction(reactants=[ghost], products=[gas],
                       kinetics=Arrhenius(A=(1.0, "1/s"), n=0.0,
                                          Ea=(10.0, "kJ/mol"), T0=(1.0, "K")),
                       reversible=True)
        rs = HybridPolymerSystem(
            T=800.0, P=1.0e5, initial_mole_fractions={gas: 0.0}, V_poly=1.0,
            polymer_pools=[], mass_transfer=[],
            gas_species_mask=mask.copy(), constant_gas_volume=False,
            initial_polymer_moments={}, termination=[])
        with pytest.raises(ValueError, match="THERMO REFERENCE-STATE TRIPWIRE") as exc:
            rs.initialize_model(core, [rxn], [], [])
        assert "GHOST" in str(exc.value)
        assert "math domain error" not in str(exc.value)
