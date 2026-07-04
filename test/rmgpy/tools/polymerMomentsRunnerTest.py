#!/usr/bin/env python3
"""Tests for the polymer moments CLI reference runner (design spec §8)."""

import csv
import json
import logging
import os
from copy import deepcopy as _dc

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
        configured_pool_labels=["PS"],
        # styrene(2) (the monomer_routing target) is GAS since recipe
        # revision 2026-07-03-monomer-gas -- only the mu dummies are
        # condensed.
        condensed_species=mus,
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

    def test_rejects_weaklink_vocabulary_in_2_1_stamped_artifact(self, qssa_deck):
        """SCHEMA/VOCABULARY CONSISTENCY (weak-link milestone iv, replaces
        the milestone-i unknown-key pin WITHOUT dropping its concern): the
        weak-link allyl/U-state vocabulary entered the sidecar at schema
        2.2, and the emitter stamps 2.2 whenever it writes it. A 2.1-stamped
        artifact must not smuggle U keys -- reject, never load the channel
        permissively (the milestone-i "no laundering under 2.1" pin)."""
        artifact = _load_artifact(qssa_deck)
        assert artifact["schema_version"] == "2.1"
        block = artifact["pools"][0]["channels"]["radical_qssa_unzip"]
        block["initiation_allyl"] = {
            "A": 2.0e14, "n": 0.0, "Ea": 2.4e5,
            "units": {"A": "s^-1", "Ea": "J/mol"}}
        block["unsaturated_tail_ends_initial"] = {
            "value": 0.02, "units": WEAKLINK_U0_UNITS}
        with pytest.raises(ValueError, match=r"schema_version.*2\.2"):
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

    def test_rejects_unknown_future_minor(self, qssa_deck, deck):
        """Weak-link milestone iv POLICY CHANGE (was minor-permissive): the
        loader pins the maximum schema minor it implements (2.5 since the
        spawned-pool closure; was 2.4 at the refused-row marker, 2.3 at the
        explicit-DP handshake block). A newer-minor artifact may carry
        vocabulary outside the channel blocks that the unknown-key guards
        never see (new conventions, new pool fields), so an older loader
        must fail loud instead of loading additively."""
        artifact = _load_artifact(qssa_deck)
        artifact["schema_version"] = "2.6"
        with pytest.raises(ValueError, match=r"schema_version.*2\.6"):
            _build_qssa(qssa_deck, artifact)
        artifact = _load_artifact(qssa_deck)
        artifact["schema_version"] = "2.8"
        with pytest.raises(ValueError, match=r"schema_version"):
            _build_qssa(qssa_deck, artifact)
        # no-QSSA legacy artifact: same envelope pin
        chem_path, art_path = deck
        with open(art_path) as fh:
            legacy = json.load(fh)
        legacy["schema_version"] = "2.6"
        species, reactions = load_chem_yaml(chem_path)
        with pytest.raises(ValueError, match=r"schema_version.*2\.6"):
            build_system_from_artifact(
                legacy, species, reactions, T0=800.0, P=1.0e5, V_poly=1.0,
                initial_moles={"N2(1)": 1.0}, mass_transfer_spec=[])

    def test_accepts_2_5_without_spawned(self, deck):
        """2.5 is now an implemented minor: a 2.5 stamp with no
        conventions.spawned_pools anywhere loads (mirror of
        test_accepts_2_4_without_refused)."""
        chem_path, art_path = deck
        with open(art_path) as fh:
            artifact = json.load(fh)
        artifact["schema_version"] = "2.5"
        species, reactions = load_chem_yaml(chem_path)
        rs, _, _ = build_system_from_artifact(
            artifact, species, reactions, T0=800.0, P=1.0e5, V_poly=1.0,
            initial_moles={"N2(1)": 1.0}, mass_transfer_spec=[])
        assert len(rs.polymer_pools) == 1

    def test_rejects_spawned_pools_in_2_4_stamped_artifact(self, deck):
        """Vocabulary/version cross-check (the refused-row precedent): a
        below-2.5 artifact carrying conventions.spawned_pools is malformed
        -- the emitter stamps 2.5 whenever it writes the key."""
        chem_path, art_path = deck
        with open(art_path) as fh:
            artifact = json.load(fh)
        artifact["schema_version"] = "2.4"
        artifact["conventions"]["spawned_pools"] = ["poly_d1"]
        species, reactions = load_chem_yaml(chem_path)
        with pytest.raises(ValueError, match=r"spawned_pools.*2\.5"):
            build_system_from_artifact(
                artifact, species, reactions, T0=800.0, P=1.0e5, V_poly=1.0,
                initial_moles={"N2(1)": 1.0}, mass_transfer_spec=[])

    def test_rejects_spawned_configured_overlap(self, deck):
        """The closure is the configured set's complement by construction: a
        label in BOTH lists is simultaneously solver-configured and
        solver-inert -- malformed, reject-never-adapt."""
        chem_path, art_path = deck
        with open(art_path) as fh:
            artifact = json.load(fh)
        artifact["schema_version"] = "2.5"
        pool_label = artifact["conventions"]["configured_pools"][0]
        artifact["conventions"]["spawned_pools"] = [pool_label]
        species, reactions = load_chem_yaml(chem_path)
        with pytest.raises(ValueError, match=r"overlaps configured_pools"):
            build_system_from_artifact(
                artifact, species, reactions, T0=800.0, P=1.0e5, V_poly=1.0,
                initial_moles={"N2(1)": 1.0}, mass_transfer_spec=[])

    def _spawned_shape_case(self, deck, schema_version, spawned, match):
        """Load a fresh deck, plant conventions.spawned_pools = ``spawned``
        under ``schema_version``, and assert the loader rejects it."""
        chem_path, art_path = deck
        with open(art_path) as fh:
            artifact = json.load(fh)
        artifact["schema_version"] = schema_version
        artifact["conventions"]["spawned_pools"] = spawned
        species, reactions = load_chem_yaml(chem_path)
        with pytest.raises(ValueError, match=match):
            build_system_from_artifact(
                artifact, species, reactions, T0=800.0, P=1.0e5, V_poly=1.0,
                initial_moles={"N2(1)": 1.0}, mass_transfer_spec=[])

    def test_rejects_empty_spawned_pools_in_2_4_stamped_artifact(self, deck):
        """KEY PRESENCE, not truthiness, is the vocabulary signal: the
        emitter stamps 2.5 whenever it writes conventions.spawned_pools at
        all, so even an EMPTY list under a 2.4 stamp is malformed (a
        truthiness gate silently waved it through)."""
        self._spawned_shape_case(deck, "2.4", [], r"spawned_pools.*2\.5")

    def test_rejects_empty_spawned_pools_in_2_5_stamped_artifact(self, deck):
        """The emitter writes the key ONLY when the closure is non-empty
        (presence-based 2.5 stamping: spawned-free artifacts stay
        byte-identical to their older stamps), so an empty list under 2.5
        is equally malformed -- reject, never adapt."""
        self._spawned_shape_case(deck, "2.5", [], r"spawned_pools.*empty")

    def test_rejects_non_list_spawned_pools(self, deck):
        """A bare string ('poly_mod') is truthy and iterable, so a
        truthiness gate consumed it as an iterable of CHARACTERS; the
        emitter writes a list of pool labels -- anything else is
        malformed."""
        self._spawned_shape_case(deck, "2.5", "poly_mod",
                                 r"spawned_pools.*list")

    def test_rejects_malformed_spawned_pool_entries(self, deck):
        """Closure entries are registry pool LABELS: non-string entries and
        empty strings are malformed regardless of position."""
        for bad in ([42], [None], [""], ["poly_d1", ""], ["poly_d1", 7]):
            self._spawned_shape_case(deck, "2.5", bad,
                                     r"spawned_pools.*non-empty string")

    def test_rejects_duplicate_spawned_pool_labels(self, deck):
        """The closure is a SET of registry pool labels (identity-deduped
        at collection, labels disambiguated at registration): the emitter
        never repeats one, so a duplicate is malformed."""
        self._spawned_shape_case(deck, "2.5", ["poly_d1", "poly_d1"],
                                 r"spawned_pools.*duplicate")

    def test_accepts_2_4_without_refused(self, deck):
        """2.4 is now an implemented minor: a 2.4 stamp with no refused row
        anywhere loads (mirror of test_accepts_2_3_without_explicit_dp)."""
        chem_path, art_path = deck
        with open(art_path) as fh:
            artifact = json.load(fh)
        artifact["schema_version"] = "2.4"
        species, reactions = load_chem_yaml(chem_path)
        rs, _, _ = build_system_from_artifact(
            artifact, species, reactions, T0=800.0, P=1.0e5, V_poly=1.0,
            initial_moles={"N2(1)": 1.0}, mass_transfer_spec=[])
        assert int(rs.reaction_refused.sum()) == 0

    def test_accepts_2_3_without_explicit_dp(self, deck):
        """2.3 is now an implemented minor: a 2.3 stamp with no explicit_dp
        block anywhere loads (mirror of test_accepts_2_1_without_qssa)."""
        chem_path, art_path = deck
        with open(art_path) as fh:
            artifact = json.load(fh)
        artifact["schema_version"] = "2.3"
        species, reactions = load_chem_yaml(chem_path)
        rs, _, _ = build_system_from_artifact(
            artifact, species, reactions, T0=800.0, P=1.0e5, V_poly=1.0,
            initial_moles={"N2(1)": 1.0}, mass_transfer_spec=[])
        assert rs.polymer_pools[0].explicit_dp_to_species_index == {}


def _refused_deck(tmp_path, mark_refused=True):
    """A deck whose artifact carries one pool-mapped refused row (item 18
    stamp-but-keep: poly(2) + H -> R + H2, the chain radical R condensed-
    classified past Gate B exactly like production).
    ``mark_refused=False`` builds the same chemistry WITHOUT the refused
    stamp -- the pre-2.4 artifact shape.

    NOTE: the runner's volatile_ejection/1 vocabulary gap this note used to
    record is CLOSED -- consumer-side VE coverage (load, numeric bundle
    pins, refused-VE zero) lives in TestVolatileEjectionConsumer below; the
    solver-level refused-vs-VE regression pin remains in solverPolymerTest.py
    test_refused_zero_while_volatile_ejection_contributes."""
    n2 = _spc("N#N", "N2", index=1)
    poly = _spc("CCCC", "poly", index=2)
    rad = _spc("[CH2]CCC", "R", index=3)
    h = _spc("[H]", "H", index=4)
    h2 = _spc("[H][H]", "H2", index=5)
    mus = [_mu("poly_mu0"), _mu("poly_mu1"), _mu("poly_mu2")]
    core = [n2, poly, rad, h, h2] + mus
    poly.is_polymer_proxy = True

    refused_rxn = Reaction(
        reactants=[poly, h], products=[rad, h2],
        kinetics=Arrhenius(A=(2.0, "m^3/(mol*s)"), n=0.0, Ea=(0.0, "J/mol"),
                           T0=(1.0, "K")),
        reversible=False)
    refused_rxn.polymer_flux_archetype = int(PolymerFluxArchetype.UNRESOLVED)
    if mark_refused:
        refused_rxn.polymer_refused = True
        refused_rxn.polymer_refused_accumulating = False
    rxns = [refused_rxn]

    data, index_map = generate_cantera_data(core, [refused_rxn],
                                            return_reaction_index_map=True)
    chem_path = os.path.join(str(tmp_path), "chem.yaml")
    with open(chem_path, "w") as fh:
        yaml.dump(data, fh, sort_keys=False, default_flow_style=None)

    pool = Polymer(label="poly", monomer="[CH2][CH2]",
                   end_groups=["[H]", "[H]"], cutoff=3,
                   moments=[1.0, 5.0, 30.0], initial_mass=0.0,
                   k_scission=0.0, k_unzip=0.0)
    artifact = build_polymer_moments_artifact(
        [pool], core_species=core, core_reactions=rxns,
        configured_pool_labels=["poly"],
        # R past Gate B: the leaked chain radical is condensed-classified,
        # exactly the production refused shape (all-gas products would be
        # phase-gate-disqualified and never reach the refused suppression).
        condensed_species=[poly, rad] + mus,
        cantera_index_map=index_map)
    art_path = os.path.join(str(tmp_path), "polymer_pools.json")
    with open(art_path, "w") as fh:
        json.dump(artifact, fh, indent=2, default=str)
    return chem_path, art_path


def _build_refused(chem_path, artifact):
    species, reactions = load_chem_yaml(chem_path)
    return build_system_from_artifact(
        artifact, species, reactions, T0=800.0, P=1.0e5, V_poly=1.0,
        initial_moles={"N2(1)": 1.0, "H(4)": 1.0}, mass_transfer_spec=[])


class TestRefusedRowConsumer:
    """Schema 2.4 refused-row contract, consumer side (reference loader):
    honor the marker (zero the row's WHOLE flux, moment and species alike --
    the generating solver's reaction_refused suppression restored across the
    artifact boundary), guard it (pool-mapped rows only, stamped >= 2.4,
    emitter never writes refused: false), and keep non-refused rows
    integrating."""

    def _artifact(self, tmp_path, **kw):
        chem_path, art_path = _refused_deck(tmp_path, **kw)
        with open(art_path) as fh:
            return chem_path, json.load(fh)

    def test_refused_row_emitted_marked_and_stamped_2_4(self, tmp_path):
        _, artifact = self._artifact(tmp_path)
        assert artifact["schema_version"] == "2.4"
        marked = [e for e in artifact["reactions"] if e.get("refused")]
        assert len(marked) == 1
        assert marked[0]["refused"] is True
        assert marked[0]["refused_reason"] == "conduit-deferred"
        assert marked[0]["proxy_reactants"] == ["poly(2)"]

    def test_refused_row_restamped_and_contributes_exactly_zero(
            self, tmp_path):
        """Pin 5 (RMG half, numbers not strings): the consumer-world oracle
        rebuilt from a 2.4 artifact suppresses the refused row's ENTIRE
        residual contribution -- with the refused row as the only reaction,
        every residual entry is exactly zero."""
        chem_path, artifact = self._artifact(tmp_path)
        rs, core, _ = _build_refused(chem_path, artifact)
        assert int(rs.reaction_refused.sum()) == 1
        dn = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
        assert np.all(dn == 0.0), f"refused row leaked flux: {dn}"

    def test_stripped_marker_reproduces_pre_2_4_over_integration(
            self, tmp_path):
        """Pin 3 (RMG half): the SAME row without the marker -- a pre-2.4
        artifact shape -- integrates the flux the generating solver zeroed
        (mu1 drained, chain radical R fabricated). This is the erratum the
        format doc records: pre-2.4 sidecars with refused RMG rows
        over-integrate in consumers."""
        chem_path, art_path = _refused_deck(tmp_path, mark_refused=False)
        with open(art_path) as fh:
            artifact = json.load(fh)
        assert artifact["schema_version"] == "2.0"  # pre-2.4 shape
        assert all("refused" not in e for e in artifact["reactions"])
        rs, core, _ = _build_refused(chem_path, artifact)
        assert int(rs.reaction_refused.sum()) == 0
        dn = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
        labels = [s.label for s in core]
        i = {lab: k for k, lab in enumerate(labels)}
        # the pre-2.4 over-integration, numerically: mu1 drains while the
        # chain radical R gains -- MW-weighted backbone mass fabrication
        assert dn[i["poly_mu1"]] < 0.0
        assert dn[i["R(3)"]] > 0.0

    def test_rejects_refused_marker_in_2_3_stamped_artifact(self, tmp_path):
        """Vocabulary/version cross-check (the 2.1-QSSA / 2.3-explicit-dp
        precedent): a below-2.4 artifact carrying the refused marker is
        malformed -- the emitter stamps 2.4 whenever it writes one."""
        chem_path, artifact = self._artifact(tmp_path)
        artifact["schema_version"] = "2.3"
        with pytest.raises(ValueError, match=r"schema_version.*2\.4"):
            _build_refused(chem_path, artifact)

    def test_rejects_refused_on_non_pool_mapped_row(self, tmp_path):
        """Refused is legal ONLY on pool-mapped rows: a hand-edited row with
        the marker but no proxy participant is corrupt (the suppression
        would silently zero ordinary gas chemistry)."""
        chem_path, artifact = self._artifact(tmp_path)
        (row,) = [e for e in artifact["reactions"] if e.get("refused")]
        row["proxy_reactants"] = []
        row["proxy_products"] = []
        with pytest.raises(ValueError, match=r"refused.*pool-mapped|pool-mapped.*refused"):
            _build_refused(chem_path, artifact)

    def test_rejects_refused_false_and_reason_shapes(self, tmp_path):
        """The emitter writes refused: true + a non-empty refused_reason or
        NOTHING (absent, not false). Present-but-false, a missing/empty
        reason, and a reason without the marker are all malformed -- reject,
        never adapt."""
        chem_path, artifact = self._artifact(tmp_path)
        (row,) = [e for e in artifact["reactions"] if e.get("refused")]

        bad = json.loads(json.dumps(artifact))
        (brow,) = [e for e in bad["reactions"] if "refused" in e]
        brow["refused"] = False
        with pytest.raises(ValueError, match=r"refused"):
            _build_refused(chem_path, bad)

        bad = json.loads(json.dumps(artifact))
        (brow,) = [e for e in bad["reactions"] if "refused" in e]
        del brow["refused_reason"]
        with pytest.raises(ValueError, match=r"refused_reason"):
            _build_refused(chem_path, bad)

        bad = json.loads(json.dumps(artifact))
        for e in bad["reactions"]:
            if "refused" in e:
                del e["refused"]
        # reason without marker (artifact still stamped 2.4 -- shape alone)
        with pytest.raises(ValueError, match=r"refused_reason"):
            _build_refused(chem_path, bad)

    def test_rejects_unknown_refused_reason(self, tmp_path):
        """The refused_reason vocabulary is CLOSED (format doc §12): exactly
        'conduit-deferred' / 'qssa-invalid'. _restamp_and_extend reconstructs
        the accumulating class FROM the reason, so an unknown string would
        silently coerce to non-accumulating semantics -- reject at load,
        naming the allowed vocabulary, before the coercion runs."""
        chem_path, artifact = self._artifact(tmp_path)
        (row,) = [e for e in artifact["reactions"] if e.get("refused")]
        row["refused_reason"] = "bogus"
        with pytest.raises(ValueError,
                           match=r"conduit-deferred.*qssa-invalid"):
            _build_refused(chem_path, artifact)

    def test_refused_reason_bijection_restores_accumulating_class(
            self, tmp_path):
        """The emitter bijection round-trips: a valid 'qssa-invalid' reason
        restores the accumulating stamp on restamp, 'conduit-deferred'
        restores non-accumulating. (A 'bogus' reason never reaches this
        coercion -- test_rejects_unknown_refused_reason pins the earlier
        rejection.)"""
        for reason, accumulating in (("qssa-invalid", True),
                                     ("conduit-deferred", False)):
            chem_path, artifact = self._artifact(tmp_path)
            (row,) = [e for e in artifact["reactions"] if e.get("refused")]
            row["refused_reason"] = reason
            species, reactions = load_chem_yaml(chem_path)
            _restamp_and_extend(artifact, species, reactions)
            (rxn,) = [r for r in reactions
                      if getattr(r, "polymer_refused", False)]
            assert rxn.polymer_refused_accumulating is accumulating, (
                f"reason {reason!r} must restore accumulating="
                f"{accumulating}")


def _ve_deck(tmp_path, eject_units=2.0, mark_refused=False):
    """A cross-pool VOLATILE_EJECTION deck: polyA(2) -> vol(4) + polyB(3)
    (interior/mu1-scaled; the moved chain loses ``eject_units`` backbone
    units to the discrete gas volatile as it lands in polyB), written out as
    chem.yaml + polymer_pools.json exactly like an RMG run would. The VE row
    is proxy-touching, so chem.yaml drops it and the artifact entry carries
    cantera: null + its own kinetics -- the same on-disk shape production PS
    runs emit (dRetroene rows: cross-pool, scaling mu1,
    params={eject_units})."""
    n2 = _spc("N#N", "N2", index=1)
    pa = _spc("CCCC", "polyA", index=2)
    pb = _spc("CCC", "polyB", index=3)
    vol = _spc("C=C", "vol", index=4)
    mus = ([_mu(f"polyA_mu{k}") for k in range(3)]
           + [_mu(f"polyB_mu{k}") for k in range(3)])
    core = [n2, pa, pb, vol] + mus
    pa.is_polymer_proxy = True
    pb.is_polymer_proxy = True

    ve_rxn = Reaction(
        reactants=[pa], products=[vol, pb],
        kinetics=Arrhenius(A=(2.5, "1/s"), n=0.0, Ea=(0.0, "J/mol"),
                           T0=(1.0, "K")),
        reversible=False)
    ve_rxn.polymer_flux_archetype = int(
        PolymerFluxArchetype.VOLATILE_EJECTION)
    ve_rxn.polymer_eject_units = eject_units
    if mark_refused:
        ve_rxn.polymer_refused = True
        ve_rxn.polymer_refused_accumulating = False

    data, index_map = generate_cantera_data(core, [ve_rxn],
                                            return_reaction_index_map=True)
    chem_path = os.path.join(str(tmp_path), "chem.yaml")
    with open(chem_path, "w") as fh:
        yaml.dump(data, fh, sort_keys=False, default_flow_style=None)

    pool_a = Polymer(label="polyA", monomer="[CH2][CH2]",
                     end_groups=["[H]", "[H]"], cutoff=3,
                     moments=[1.0, 10.0, 120.0], initial_mass=0.0,
                     k_scission=0.0, k_unzip=0.0)
    pool_b = Polymer(label="polyB", monomer="[CH2][CH2]",
                     end_groups=["[H]", "[H]"], cutoff=3,
                     moments=[1.0, 5.0, 30.0], initial_mass=0.0,
                     k_scission=0.0, k_unzip=0.0)
    artifact = build_polymer_moments_artifact(
        [pool_a, pool_b], core_species=core, core_reactions=[ve_rxn],
        configured_pool_labels=["polyA", "polyB"],
        condensed_species=[pa, pb] + mus,
        cantera_index_map=index_map)
    art_path = os.path.join(str(tmp_path), "polymer_pools.json")
    with open(art_path, "w") as fh:
        json.dump(artifact, fh, indent=2, default=str)
    return chem_path, art_path


def _build_ve(chem_path, artifact):
    species, reactions = load_chem_yaml(chem_path)
    return build_system_from_artifact(
        artifact, species, reactions, T0=800.0, P=1.0e5, V_poly=1.0,
        initial_moles={"N2(1)": 1.0}, mass_transfer_spec=[])


class TestVolatileEjectionConsumer:
    """volatile_ejection/1 vocabulary in the reference runner (pre-S2
    blocker, round 59): the emitter has written VE rows since the signed-VE
    spec (rmgpy/polymer.py ARCHETYPE_TERM_NAMES + params.eject_units) and
    shipped PS artifacts carry them, but ARCHETYPE_INTS had no entry, so any
    sidecar with a VE row crashed the runner with
    KeyError: 'volatile_ejection/1'. The consumer must restore BOTH stamps
    the oracle reads (polymer_flux_archetype = 6 and the SIGNED
    polymer_eject_units; src/dst pools are re-resolved from the species like
    every other archetype) and reject malformed eject vocabulary with
    actionable errors, never a silent 0.0 default."""

    def _artifact(self, tmp_path, **kw):
        chem_path, art_path = _ve_deck(tmp_path, **kw)
        with open(art_path) as fh:
            return chem_path, json.load(fh)

    def test_emitter_writes_ve_row_shape(self, tmp_path):
        """Premise pin: the emitter's one on-disk VE shape (matches the
        shipped PS artifacts, e.g. b2fix-psrun ve_1782931902): cross-pool
        row, cantera null (proxy-touching rows are dropped from chem.yaml),
        scaling mu1, params carrying exactly eject_units."""
        _, artifact = self._artifact(tmp_path)
        (row,) = artifact["reactions"]
        assert row["archetype"] == "volatile_ejection/1"
        assert row["cantera"] is None
        assert row["kinetics"] is not None
        assert row["scaling"] == "mu1"
        assert row["src_pool"] == "polyA"
        assert row["dst_pool"] == "polyB"
        assert row["params"] == {"eject_units": 2.0}
        assert row["proxy_reactants"] == ["polyA(2)"]
        assert row["proxy_products"] == ["polyB(3)"]
        assert row["unresolved"] is False

    def test_ve_row_loads_and_restores_solver_stamps(self, tmp_path):
        """Pre-fix this died with KeyError: 'volatile_ejection/1' (the
        ARCHETYPE_INTS lookup in _restamp_and_extend). Post-fix the oracle
        carries the VE archetype int, the eject-units stamp, and the
        re-resolved src/dst pools."""
        chem_path, artifact = self._artifact(tmp_path)
        rs, _, _ = _build_ve(chem_path, artifact)
        assert int(rs.reaction_flux_archetype[0]) == 6
        assert float(rs.reaction_eject_units[0]) == pytest.approx(2.0)
        assert int(rs.reaction_src_pool[0]) == 0
        assert int(rs.reaction_dst_pool[0]) == 1

    def test_signed_negative_eject_units_round_trips(self, tmp_path):
        """eject_units is SIGNED (polymer.py signed-VE contract: a < 0 is
        net mass GAIN from the gas co-reactants); the consumer must restore
        the sign, never clamp."""
        chem_path, artifact = self._artifact(tmp_path, eject_units=-1.5)
        rs, _, _ = _build_ve(chem_path, artifact)
        assert float(rs.reaction_eject_units[0]) == pytest.approx(-1.5)

    def test_ve_forward_bundle_and_gas_pins(self, tmp_path):
        """Numeric behavioral pins through the real runner path (numbers,
        not strings). Forward cross-pool VE with mu(polyA) = [1, 10, 120]
        mol, V_poly = 1 m^3, kf = 2.5 s^-1, a = eject_units = 2,
        hand-computed against the solver VE law (polymer.pyx VE branch +
        _chain_bundle):

          event rate  ev = kf * mu1 = 25 mol/s   (interior: mu1-scaled site)
          length-biased bundle: b0 = 1, b1 = mu2/mu1 = 12,
              mu3 = mu0*(mu2/mu1)^3 = 1728 (log-Lagrange closure),
              b2 = mu3/mu1 = 172.8
          src debit (full bundle):  -ev*(b0, b1, b2) = (-25, -300, -4320)
          dst credit (a-shifted, sa = -a): +ev*(b0, b1 - a,
              b2 - 2*a*b1 + a^2) = (+25, +250, +3220)
          gas volatile (standard net-rate product path): +ev = +25
        """
        chem_path, artifact = self._artifact(tmp_path)
        rs, core, _ = _build_ve(chem_path, artifact)
        dn = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
        i = {s.label: k for k, s in enumerate(core)}
        ev = 2.5 * 10.0
        assert dn[i["polyA_mu0"]] == pytest.approx(-ev, rel=1e-12)
        assert dn[i["polyA_mu1"]] == pytest.approx(-ev * 12.0, rel=1e-12)
        assert dn[i["polyA_mu2"]] == pytest.approx(-ev * 172.8, rel=1e-12)
        assert dn[i["polyB_mu0"]] == pytest.approx(+ev, rel=1e-12)
        assert dn[i["polyB_mu1"]] == pytest.approx(+ev * (12.0 - 2.0),
                                                   rel=1e-12)
        assert dn[i["polyB_mu2"]] == pytest.approx(
            +ev * (172.8 - 2.0 * 2.0 * 12.0 + 2.0 ** 2), rel=1e-12)
        assert dn[i["vol(4)"]] == pytest.approx(+ev, rel=1e-12)

    def test_refused_ve_row_contributes_exactly_zero(self, tmp_path):
        """A refused VE row (schema 2.4 stamp-but-keep) must suppress its
        WHOLE flux in the consumer oracle too: with the refused VE row as
        the only reaction, every residual entry is exactly zero."""
        chem_path, artifact = self._artifact(tmp_path, mark_refused=True)
        (row,) = artifact["reactions"]
        assert artifact["schema_version"] == "2.4"
        assert row["refused"] is True
        assert row["archetype"] == "volatile_ejection/1"
        rs, _, _ = _build_ve(chem_path, artifact)
        assert int(rs.reaction_refused.sum()) == 1
        dn = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
        assert np.all(dn == 0.0), f"refused VE row leaked flux: {dn}"

    def test_rejects_malformed_eject_params_never_keyerror(self, tmp_path):
        """The emitter writes params = {'eject_units': float} on every VE
        row (rmgpy/polymer.py compile_polymer_reaction_entries) -- the ONLY
        VE params sub-shape that exists on disk. Anything else is a
        hand-edited/corrupted sidecar: reject with an actionable ValueError
        naming eject_units -- never KeyError, never a silent 0.0 default
        (which would launder the atom-transfer debit away: the chain lands
        un-shrunk while the gas volatile still appears)."""
        cases = [
            ("missing params", lambda r: r.pop("params")),
            ("empty params", lambda r: r.__setitem__("params", {})),
            ("chip vocabulary", lambda r: r.__setitem__("params", {"a": 2})),
            ("extra key", lambda r: r["params"].__setitem__("a", 2)),
            ("string value", lambda r: r["params"].__setitem__(
                "eject_units", "2.0")),
            ("boolean value", lambda r: r["params"].__setitem__(
                "eject_units", True)),
            ("non-finite", lambda r: r["params"].__setitem__(
                "eject_units", float("nan"))),
        ]
        for name, mutate in cases:
            chem_path, artifact = self._artifact(tmp_path)
            (row,) = artifact["reactions"]
            mutate(row)
            with pytest.raises(ValueError, match=r"eject_units"):
                _build_ve(chem_path, artifact)
            del artifact  # each case mutates a fresh load

    def test_rejects_unknown_archetype_with_actionable_error(self, tmp_path):
        """The archetype term-type vocabulary is CLOSED: an unknown name
        means flux this consumer cannot reproduce -- ValueError naming the
        offender and the known vocabulary, not the bare KeyError that was
        the pre-fix failure mode."""
        chem_path, artifact = self._artifact(tmp_path)
        (row,) = artifact["reactions"]
        row["archetype"] = "volatile_ejection/2"
        with pytest.raises(ValueError,
                           match=r"unknown archetype.*volatile_ejection/2"):
            _build_ve(chem_path, artifact)


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


# --- Weak-link allyl/U-state channel (sidecar schema 2.2, milestone iv) ---

# Raw deck-shaped weak-link config: split termination REPLACES the legacy
# summed block; the four weak-link keys are all-or-nothing.
WEAKLINK_RAW_CFG = {
    "initiation": {"A": 1.0e13, "n": 0.0, "Ea": 3.0e5},
    "depropagation": {"A": 1.0e14, "n": 0.5, "Ea": 9.0e4},
    "transfer": {"A": 2.0e3, "n": 0.0, "Ea": 5.0e4},
    "initiation_allyl": {"A": 2.0e14, "n": 0.0, "Ea": 2.4e5},
    "termination_recombination": {"A": 6.0e7, "n": 0.0, "Ea": 8.0e3},
    "termination_disproportionation": {"A": 4.0e7, "n": 0.0, "Ea": 1.2e4},
    "unsaturated_tail_ends_initial": 0.02,
    "efficiency": 0.8,
    "monomer_yield": 0.9,
}

# Pinned U0 units note -- must match the emitter byte-for-byte.
WEAKLINK_U0_UNITS = "mol — tail-distribution state; consumer divides by V_poly"


def _weaklink_deck(tmp_path, u0=None):
    """A weak-link QSSA pool with wired monomer routing, written out as
    chem.yaml + polymer_pools.json exactly like an RMG run would."""
    n2 = _spc("N#N", "N2", index=1)
    sty = _spc("C=Cc1ccccc1", "styrene", index=2)
    mus = [_mu("PS_mu0"), _mu("PS_mu1"), _mu("PS_mu2")]
    core = [n2, sty] + mus
    data, index_map = generate_cantera_data(core, [],
                                            return_reaction_index_map=True)
    chem_path = os.path.join(str(tmp_path), "chem.yaml")
    with open(chem_path, "w") as fh:
        yaml.dump(data, fh, sort_keys=False, default_flow_style=None)

    cfg = dict(WEAKLINK_RAW_CFG)
    if u0 is not None:
        cfg["unsaturated_tail_ends_initial"] = u0
    pool = Polymer(label="PS", monomer="[CH2][CH](c1ccccc1)",
                   end_groups=["[H]", "[H]"], cutoff=3,
                   moments=[1.0, 50.0, 3000.0], initial_mass=0.0,
                   k_scission=0.0, k_unzip=0.0,
                   radical_qssa_unzip=cfg)
    artifact = build_polymer_moments_artifact(
        [pool], core_species=core, core_reactions=[],
        configured_pool_labels=["PS"],
        # styrene(2) (the monomer_routing target) is GAS since recipe
        # revision 2026-07-03-monomer-gas -- only the mu dummies are
        # condensed.
        condensed_species=mus,
        monomer_routing_by_pool={"PS": "styrene(2)"},
        cantera_index_map=index_map)
    art_path = os.path.join(str(tmp_path), "polymer_pools.json")
    with open(art_path, "w") as fh:
        json.dump(artifact, fh, indent=2, default=str)
    return chem_path, art_path


@pytest.fixture
def weaklink_deck(tmp_path):
    return _weaklink_deck(tmp_path)


class TestWeakLinkArtifactLoader:
    """Sidecar schema 2.2 loader (milestone iv): the weak-link vocabulary
    round-trips config -> sidecar -> load -> config exactly, the same
    single-source validator covers the loaded keys, and the solver-side U0
    census fires end-to-end through the runner path (no laundering
    window)."""

    def test_weaklink_artifact_stamps_2_2_and_loads(self, weaklink_deck):
        artifact = _load_artifact(weaklink_deck)
        assert artifact["schema_version"] == "2.2"
        assert artifact["conventions"]["format_doc"] == (
            "docs/polymer_moments_format.md (polymer_moments_format/2.2)")
        assert artifact["conventions"]["recipe_revision"] == \
            "2026-07-03-weaklink-u-monomer-gas"
        rs, core, _ = _build_qssa(weaklink_deck, artifact)
        assert rs.qssa_enabled[0] == 1
        assert rs.qssa_weaklink[0] == 1
        assert rs.qssa_kia_A[0] == pytest.approx(2.0e14)
        assert rs.qssa_kia_Ea[0] == pytest.approx(2.4e5)
        assert rs.qssa_ktrec_A[0] == pytest.approx(6.0e7)
        assert rs.qssa_ktdisp_A[0] == pytest.approx(4.0e7)
        assert rs.qssa_u0[0] == pytest.approx(0.02)
        # the U slot exists: one trailing ODE slot beyond the core block
        assert rs.neq == len(core) + 1

    def test_weaklink_round_trip_config_identity(self, weaklink_deck):
        """config -> sidecar -> load -> config must be IDENTICAL to the
        shared validator's normalization of the deck config (dict
        equality, not approx): the loader feeds exactly the shape
        validate_radical_qssa_unzip accepts."""
        from rmgpy.solver.polymer import validate_radical_qssa_unzip

        rs, _, _ = _build_qssa(weaklink_deck)
        expected = validate_radical_qssa_unzip("PS", dict(WEAKLINK_RAW_CFG))
        assert rs.polymer_pools[0].radical_qssa_unzip == expected

    def test_nan_in_loaded_initiation_allyl_rejected(self, weaklink_deck):
        """Boundary guard coverage for the new keys: a NaN smuggled into a
        loaded weak-link triplet dies in the single-source validator, not
        silently disabling (or poisoning) the channel."""
        artifact = _load_artifact(weaklink_deck)
        block = artifact["pools"][0]["channels"]["radical_qssa_unzip"]
        block["initiation_allyl"]["A"] = float("nan")
        with pytest.raises(ValueError,
                           match=r"PS.*initiation_allyl\.A.*finite"):
            _build_qssa(weaklink_deck, artifact)

    def test_wrong_units_on_weaklink_triplet_rejected(self, weaklink_deck):
        """The split termination triplets carry the bimolecular pin of the
        block they replace; a sidecar claiming s^-1 must ERROR, never be
        converted."""
        artifact = _load_artifact(weaklink_deck)
        block = artifact["pools"][0]["channels"]["radical_qssa_unzip"]
        block["termination_disproportionation"]["units"]["A"] = "s^-1"
        with pytest.raises(
                ValueError,
                match=r"PS.*termination_disproportionation.*units"):
            _build_qssa(weaklink_deck, artifact)

    def test_wrong_u0_units_note_rejected(self, weaklink_deck):
        """U0's units note is pinned byte-for-byte: a bare 'mol' (or any
        other string) means the amount-basis contract is not the one this
        loader implements."""
        artifact = _load_artifact(weaklink_deck)
        block = artifact["pools"][0]["channels"]["radical_qssa_unzip"]
        block["unsaturated_tail_ends_initial"]["units"] = "mol"
        with pytest.raises(
                ValueError,
                match=r"PS.*unsaturated_tail_ends_initial.*units"):
            _build_qssa(weaklink_deck, artifact)

    def test_weaklink_recipe_pinned_exactly(self, weaklink_deck):
        """The machine-pinned weak-link law: any mutation or omission of a
        recipe field is rejected (reject, never adapt), and the legacy
        recipe alone is NOT accepted for a weak-link block."""
        artifact = _load_artifact(weaklink_deck)
        recipe = artifact["pools"][0]["channels"]["radical_qssa_unzip"]["recipe"]
        recipe["du_dt"] = "dU/dt = kt_disp*R_ss^2*V_poly"
        with pytest.raises(ValueError, match=r"PS.*du_dt"):
            _build_qssa(weaklink_deck, artifact)

        artifact = _load_artifact(weaklink_deck)
        recipe = artifact["pools"][0]["channels"]["radical_qssa_unzip"]["recipe"]
        del recipe["radical_generation"]
        with pytest.raises(ValueError, match=r"PS.*radical_generation"):
            _build_qssa(weaklink_deck, artifact)

    def test_u0_census_fires_end_to_end_through_runner_path(self, tmp_path):
        """Test 7 ruling: the U0 > 2*mu0_tail census lives in the SOLVER
        (set_initial_conditions); the loader passes the value through and
        the solver rejects at init. Prove artifact-in -> error-out through
        the full runner path so no laundering window exists: mu0 = 1.0 mol
        -> tail capacity 2.0 mol; U0 = 5.0 mol must die at build time."""
        deck = _weaklink_deck(tmp_path, u0=5.0)
        with pytest.raises(
                ValueError,
                match=r"PS.*unsaturated_tail_ends_initial.*chain-end "
                      r"capacity"):
            _build_qssa(deck)


class TestSchema21FixtureRegression:
    """Two FROZEN 2.1 artifacts (channel constants and shape derived from
    the kdpswap PS-rerun baseline sidecar; deliberately NOT
    emitter-produced, so this regression is independent of emitter
    changes), split by recipe revision -- P1-B, the anti-laundering pin:

    * FROZEN_2_1_LEGACY_ARTIFACT is the TRUE pre-monomer-gas freeze
      (recipe_revision 2026-07-02, routed monomer styrene(2) CONDENSED;
      byte-content recovered from the b35d54047 freeze). The loader must
      HARD-REFUSE it with a regenerate message: this loader implements
      only the gas-monomer semantics, and legacy acceptance would
      re-condense the routed monomer -- the exact reference-state
      conflation revision 2026-07-03-monomer-gas removed.
    * FROZEN_2_1_MONOMER_GAS_ARTIFACT carries the NEW revision
      (2026-07-03-qssa-monomer-gas, routed monomer GAS) and pins the
      QSSA-parse regression the old single fixture used to pin.

    A fixture carrying the OLD revision token with NEW gas semantics (or
    vice versa) is semantically contradictory and is pinned REJECTED
    below -- re-freezing across the revision boundary is the laundering
    vector this split exists to forbid."""

    FROZEN_2_1_LEGACY_ARTIFACT = {
        "schema_version": "2.1",
        "rmg_commit": "kdpswap-baseline-derived",
        "rmg_iteration": 0,
        "conventions": {
            "format_doc": ("docs/polymer_moments_format.md "
                           "(polymer_moments_format/2.1)"),
            "recipe_revision": "2026-07-02",
            "configured_pools": ["PS"],
            "condensed_species": ["PS_mu0", "PS_mu1", "PS_mu2",
                                  "styrene(2)"],
        },
        "pools": [{
            "label": "PS",
            "cutoff": 3,
            "moments": [0.01, 0.4800819207743936, 27.65743807853174],
            "monomer_mw_g_mol": 104.14962093690001,
            "channels": {
                "scission": {"A": 0.0, "n": 0.0, "Ea": 0.0,
                             "units": {"A": "s^-1", "Ea": "J/mol"}},
                "unzip": {"A": 0.0, "n": 0.0, "Ea": 0.0,
                          "units": {"A": "s^-1", "Ea": "J/mol"}},
                "radical_qssa_unzip": {
                    "enabled": True,
                    "basis": "backbone_bonds_mu1_minus_mu0",
                    "efficiency": 1.0,
                    "monomer_yield": 1.0,
                    "initiation": {"A": 1.0e15, "n": 0.0, "Ea": 281600.0,
                                   "units": {"A": "s^-1", "Ea": "J/mol"}},
                    "depropagation": {"A": 3.1e12, "n": 0.0, "Ea": 100000.0,
                                      "units": {"A": "s^-1",
                                                "Ea": "J/mol"}},
                    "termination": {"A": 57750000.0, "n": 0.0, "Ea": 9600.0,
                                    "units": {"A": "m^3/(mol*s)",
                                              "Ea": "J/mol"}},
                    "transfer": None,
                    "recipe": {
                        "bond_basis": ("B = max(mu1 - mu0, 0) on "
                                       "concentration moments (mol/m^3 "
                                       "condensed)"),
                        "rate_no_transfer": (
                            "r_mono = monomer_yield * kdp * "
                            "sqrt(efficiency * ki * B / kt)"),
                        "rate_with_transfer": (
                            "r_mono = monomer_yield * kdp * "
                            "(sqrt(ktr^2 + 8*kt*(2*efficiency*ki*B)) "
                            "- ktr) / (4*kt)"),
                        "moment_signature": (
                            "dmu0 = 0; dmu1 -= r_mono; dmu2 -= r_mono * "
                            "max(2*mu1/max(mu0, small_eps) - 1, 0)"),
                        "small_eps": 1e-30,
                        "volume_note": (
                            "kt is bimolecular: rates depend on condensed "
                            "volume V_poly; consumers MUST evaluate on "
                            "concentration moments mu_k = n_k / V_poly and "
                            "convert emitted rate back with *V_poly"),
                    },
                    "provenance": {
                        "radical_balance": (
                            "G_R = 2*f*ki*B; loss = ktr*R + 2*kt*R^2; "
                            "Rss no-transfer = sqrt(f*ki*B/kt)"),
                    },
                },
            },
            "phase_species": [],
            "bookkeeping_species": [],
            "monomer_routing": "styrene(2)",
        }],
        "reactions": [],
    }

    # The monomer-gas twin differs from the legacy freeze by EXACTLY the
    # 2026-07-03-monomer-gas revision delta (recipe_revision token + the
    # routed monomer leaving the condensed list) -- derived here so the
    # delta is explicit and nothing else can drift between the two.
    FROZEN_2_1_MONOMER_GAS_ARTIFACT = _dc(FROZEN_2_1_LEGACY_ARTIFACT)
    FROZEN_2_1_MONOMER_GAS_ARTIFACT["conventions"]["recipe_revision"] = \
        "2026-07-03-qssa-monomer-gas"
    FROZEN_2_1_MONOMER_GAS_ARTIFACT["conventions"]["condensed_species"] = \
        ["PS_mu0", "PS_mu1", "PS_mu2"]

    def test_frozen_monomer_gas_artifact_loads_with_identical_parsed_output(
            self, qssa_deck):
        from copy import deepcopy

        from rmgpy.solver.polymer import validate_radical_qssa_unzip

        chem_path, _ = qssa_deck
        species, reactions = load_chem_yaml(chem_path)
        rs, core, _ = build_system_from_artifact(
            deepcopy(self.FROZEN_2_1_MONOMER_GAS_ARTIFACT), species,
            reactions, T0=650.0, P=1.0e5, V_poly=1.0,
            initial_moles={"N2(1)": 1.0}, mass_transfer_spec=[])
        expected = validate_radical_qssa_unzip("PS", {
            "initiation": {"A": 1.0e15, "n": 0.0, "Ea": 281600.0},
            "depropagation": {"A": 3.1e12, "n": 0.0, "Ea": 100000.0},
            "termination": {"A": 57750000.0, "n": 0.0, "Ea": 9600.0},
            "transfer": None,
        })
        assert rs.polymer_pools[0].radical_qssa_unzip == expected
        assert rs.qssa_enabled[0] == 1
        assert rs.qssa_weaklink[0] == 0
        assert rs.qssa_kt_A[0] == pytest.approx(57750000.0)
        # legacy layout: no trailing U slot
        assert rs.neq == len(core)

    def _load(self, artifact, qssa_deck):
        from copy import deepcopy
        chem_path, _ = qssa_deck
        species, reactions = load_chem_yaml(chem_path)
        return build_system_from_artifact(
            deepcopy(artifact), species, reactions,
            T0=650.0, P=1.0e5, V_poly=1.0,
            initial_moles={"N2(1)": 1.0}, mass_transfer_spec=[])

    def test_frozen_legacy_artifact_hard_refused_with_regenerate_message(
            self, qssa_deck):
        """RED pin (P1-B direction 1): the TRUE legacy artifact (old
        revision, old condensed-monomer semantics) is hard-refused AT THE
        LOADER with an actionable regenerate message -- never a deep
        solver error, never silent legacy acceptance (that would
        re-condense the routed monomer, the exact defect revision
        2026-07-03-monomer-gas removed)."""
        with pytest.raises(ValueError,
                           match=r"recipe_revision.*2026-07-02.*[Rr]egenerate"
                                 r".*monomer-gas"):
            self._load(self.FROZEN_2_1_LEGACY_ARTIFACT, qssa_deck)

    def test_old_revision_with_gas_semantics_still_refused(self, qssa_deck):
        """RED pin (P1-B, the laundering vector itself): an artifact
        carrying the OLD revision token but re-frozen NEW gas semantics
        (routed monomer absent from the condensed list -- exactly the
        contradictory re-freeze the review forbids) must ALSO be refused:
        the gate is revision-keyed, not semantics-sniffed, so a re-freeze
        cannot launder an old artifact past it."""
        from copy import deepcopy
        laundered = deepcopy(self.FROZEN_2_1_LEGACY_ARTIFACT)
        laundered["conventions"]["condensed_species"] = \
            ["PS_mu0", "PS_mu1", "PS_mu2"]
        with pytest.raises(ValueError,
                           match=r"recipe_revision.*2026-07-02.*[Rr]egenerate"
                                 r".*monomer-gas"):
            self._load(laundered, qssa_deck)

    def test_new_revision_with_condensed_routed_monomer_rejected(
            self, qssa_deck):
        """RED pin (P1-B direction 2): an artifact stamping a NEW
        monomer-gas recipe_revision whose condensed list still contains
        the routed monomer is internally contradictory and must be
        REJECTED at the loader with a clear message (not mis-phased, not
        left to a deep solver error)."""
        from copy import deepcopy
        contradictory = deepcopy(self.FROZEN_2_1_MONOMER_GAS_ARTIFACT)
        contradictory["conventions"]["condensed_species"] = \
            ["PS_mu0", "PS_mu1", "PS_mu2", "styrene(2)"]
        with pytest.raises(ValueError,
                           match=r"internally contradictory.*regenerate"):
            self._load(contradictory, qssa_deck)

    def test_new_revision_with_routed_monomer_in_phase_species_rejected(
            self, qssa_deck):
        """Same contradiction through the pools[].phase_species list (the
        other membership list the monomer-gas revision removed the routed
        monomer from)."""
        from copy import deepcopy
        contradictory = deepcopy(self.FROZEN_2_1_MONOMER_GAS_ARTIFACT)
        contradictory["pools"][0]["phase_species"] = ["styrene(2)"]
        with pytest.raises(ValueError,
                           match=r"internally contradictory.*regenerate"):
            self._load(contradictory, qssa_deck)

    def test_old_revision_without_routed_monomer_still_loads(
            self, qssa_deck):
        """Scope pin: the hard refusal is keyed on routed-monomer
        semantics being in play. An OLD-revision artifact that routes no
        monomer (no monomer_routing, unzip/QSSA silent) has no
        monomer-phase semantics to mis-read and still loads."""
        from copy import deepcopy
        legacy_unrouted = deepcopy(self.FROZEN_2_1_LEGACY_ARTIFACT)
        legacy_unrouted["pools"][0]["monomer_routing"] = None
        del legacy_unrouted["pools"][0]["channels"]["radical_qssa_unzip"]
        rs, core, _ = self._load(legacy_unrouted, qssa_deck)
        assert rs.qssa_enabled[0] == 0


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


@pytest.fixture
def explicit_dp_deck(tmp_path):
    """A stage-A explicit-DP pool (deck flag explicit_dp=True attached the
    capped DP=cutoff oligomer as a real core species), written out as
    chem.yaml + polymer_pools.json exactly like an RMG run would."""
    n2 = _spc("N#N", "N2", index=1)
    dp3 = _spc("CCC", "poly_dp3", index=2)
    mus = [_mu("poly_mu0"), _mu("poly_mu1"), _mu("poly_mu2")]
    core = [n2, dp3] + mus
    data, index_map = generate_cantera_data(core, [],
                                            return_reaction_index_map=True)
    chem_path = os.path.join(str(tmp_path), "chem.yaml")
    with open(chem_path, "w") as fh:
        yaml.dump(data, fh, sort_keys=False, default_flow_style=None)

    pool = Polymer(label="poly", monomer="[CH2][CH2]",
                   end_groups=["[H]", "[H]"], cutoff=3,
                   moments=[1.0, 5.0, 30.0], initial_mass=0.0,
                   k_scission=1.0, k_unzip=0.0)
    pool.explicit_dp = True
    pool.explicit_dp_species = dp3
    artifact = build_polymer_moments_artifact(
        [pool], core_species=core, core_reactions=[],
        configured_pool_labels=["poly"], condensed_species=mus + [dp3],
        cantera_index_map=index_map,
        initial_explicit_by_pool={"poly": {3: 0.25}})
    art_path = os.path.join(str(tmp_path), "polymer_pools.json")
    with open(art_path, "w") as fh:
        json.dump(artifact, fh, indent=2, default=str)
    return chem_path, art_path


def _build_explicit(explicit_dp_deck, artifact=None):
    chem_path, art_path = explicit_dp_deck
    if artifact is None:
        with open(art_path) as fh:
            artifact = json.load(fh)
    species, reactions = load_chem_yaml(chem_path)
    return build_system_from_artifact(
        artifact, species, reactions, T0=800.0, P=1.0e5, V_poly=1.0,
        initial_moles={"N2(1)": 1.0}, mass_transfer_spec=[])


class TestExplicitDpArtifactLoader:
    """Schema-2.3 explicit_dp pool block (stage B): the loader parses and
    validates the block at the artifact boundary (pinned recipe/revision,
    closed key vocabulary, handshake invariants) and wires the oracle's
    explicit map + t=0 loadings from it."""

    def test_round_trip_wires_explicit_map_and_initial_moles(
            self, explicit_dp_deck):
        with open(explicit_dp_deck[1]) as fh:
            artifact = json.load(fh)
        assert artifact["schema_version"] == "2.3"
        assert artifact["conventions"]["recipe_revision"] == \
            "2026-07-04-explicit-dp-monomer-gas"
        block = artifact["pools"][0]["explicit_dp"]
        assert block["species"] == {"3": "poly_dp3(2)"}
        assert block["initial_moles"] == {"3": 0.25}
        rs, core, _ = _build_explicit(explicit_dp_deck, artifact)
        labels = [s.label for s in core]
        dp_idx = labels.index("poly_dp3(2)")
        assert rs.polymer_pools[0].explicit_dp_to_species_index == {3: dp_idx}
        # t=0 loading seeded as species amount (set_initial_conditions
        # step 2, clamped >= 0)
        assert rs.y0[dp_idx] == pytest.approx(0.25)
        # BEHAVIORAL tail-split pin (adversarial-review P1): declared
        # pools[].moments are TOTAL-INCLUSIVE and set_initial_conditions
        # step 6 subtracts each mapped DP's contribution (N, dp*N, dp^2*N)
        # from the seeded mu0/mu1/mu2, clamped >= 0. Declared [1, 5, 30]
        # with 0.25 mol at DP=3, V_poly=1 -> seeded mu amounts
        # [1-0.25, 5-3*0.25, 30-9*0.25] = [0.75, 4.25, 27.75]. Verbatim
        # (no-subtraction) seeding would leave [1, 5, 30] — the false claim
        # the original tail_split contract text made; this pin fails RED on
        # any emitter/oracle drift back to that reading.
        mu_seeded = [float(rs.y0[i]) for i in rs.polymer_pools[0].mu_indices]
        assert mu_seeded == pytest.approx([0.75, 4.25, 27.75])

    def test_legacy_artifact_without_block_loads_unchanged(self, deck):
        chem_path, art_path = deck
        with open(art_path) as fh:
            artifact = json.load(fh)
        assert "explicit_dp" not in json.dumps(artifact)
        species, reactions = load_chem_yaml(chem_path)
        rs, _, _ = build_system_from_artifact(
            artifact, species, reactions, T0=800.0, P=1.0e5, V_poly=1.0,
            initial_moles={"N2(1)": 1.0}, mass_transfer_spec=[])
        assert rs.polymer_pools[0].explicit_dp_to_species_index == {}

    def test_rejects_block_in_2_2_stamped_artifact(self, explicit_dp_deck):
        with open(explicit_dp_deck[1]) as fh:
            artifact = json.load(fh)
        artifact["schema_version"] = "2.2"
        with pytest.raises(ValueError, match=r"explicit_dp.*2\.3|2\.3.*explicit_dp"):
            _build_explicit(explicit_dp_deck, artifact)

    def test_rejects_disabled_block_as_malformed(self, explicit_dp_deck):
        with open(explicit_dp_deck[1]) as fh:
            artifact = json.load(fh)
        artifact["pools"][0]["explicit_dp"]["enabled"] = False
        with pytest.raises(ValueError, match="enabled=false"):
            _build_explicit(explicit_dp_deck, artifact)

    def test_rejects_tampered_recipe_no_adaptation(self, explicit_dp_deck):
        with open(explicit_dp_deck[1]) as fh:
            artifact = json.load(fh)
        artifact["pools"][0]["explicit_dp"]["recipe"]["k_chain"] = "k_unzip"
        with pytest.raises(ValueError, match="pinned normative recipe"):
            _build_explicit(explicit_dp_deck, artifact)
        with open(explicit_dp_deck[1]) as fh:
            artifact = json.load(fh)
        artifact["pools"][0]["explicit_dp"]["recipe_revision"] = "2026-01-01"
        with pytest.raises(ValueError, match="recipe_revision"):
            _build_explicit(explicit_dp_deck, artifact)

    def test_rejects_unknown_block_key(self, explicit_dp_deck):
        with open(explicit_dp_deck[1]) as fh:
            artifact = json.load(fh)
        artifact["pools"][0]["explicit_dp"]["dp_ladder"] = [2, 3]
        with pytest.raises(ValueError, match="unknown key"):
            _build_explicit(explicit_dp_deck, artifact)

    def test_rejects_unresolvable_species_label(self, explicit_dp_deck):
        with open(explicit_dp_deck[1]) as fh:
            artifact = json.load(fh)
        artifact["pools"][0]["explicit_dp"]["species"]["3"] = "GHOST(99)"
        with pytest.raises(ValueError, match="GHOST"):
            _build_explicit(explicit_dp_deck, artifact)

    def test_rejects_target_dp_cutoff_mismatch(self, explicit_dp_deck):
        with open(explicit_dp_deck[1]) as fh:
            artifact = json.load(fh)
        artifact["pools"][0]["explicit_dp"]["handshake_target_dp"] = 2
        with pytest.raises(ValueError, match="handshake_target_dp"):
            _build_explicit(explicit_dp_deck, artifact)

    def test_handshake_residual_deposits_and_drains_through_runner_path(
            self, tmp_path):
        """BEHAVIORAL handshake pin through the runner path (adversarial
        review: string pins alone gave false confidence). A runner-built
        system with an active k_unzip arm and an explicit-DP target must,
        at the RHS level, deposit F*V_poly into the explicit species and
        drain the moments by exactly (F, xs*F, xs^2*F) on top of the unzip
        channel's own (0, r, k_u*(2*mu1-mu0)) drain."""
        n2 = _spc("N#N", "N2", index=1)
        dp3 = _spc("CCC", "poly_dp3", index=2)
        mono = _spc("C=C", "C2H4", index=3)
        mus = [_mu("poly_mu0"), _mu("poly_mu1"), _mu("poly_mu2")]
        core = [n2, dp3, mono] + mus
        data, index_map = generate_cantera_data(
            core, [], return_reaction_index_map=True)
        chem_path = os.path.join(str(tmp_path), "chem.yaml")
        with open(chem_path, "w") as fh:
            yaml.dump(data, fh, sort_keys=False, default_flow_style=None)
        pool = Polymer(label="poly", monomer="[CH2][CH2]",
                       end_groups=["[H]", "[H]"], cutoff=3,
                       moments=[1.0, 5.0, 30.0], initial_mass=0.0,
                       k_scission=0.0, k_unzip=0.1)
        pool.explicit_dp = True
        pool.explicit_dp_species = dp3
        artifact = build_polymer_moments_artifact(
            [pool], core_species=core, core_reactions=[],
            configured_pool_labels=["poly"],
            condensed_species=mus + [dp3],
            monomer_routing_by_pool={"poly": "C2H4(3)"},
            cantera_index_map=index_map)
        species, reactions = load_chem_yaml(chem_path)
        rs, core2, _ = build_system_from_artifact(
            artifact, species, reactions, T0=800.0, P=1.0e5, V_poly=1.0,
            initial_moles={"N2(1)": 1.0}, mass_transfer_spec=[])
        labels = [s.label for s in core2]
        dp_idx = labels.index("poly_dp3(2)")
        mono_idx = labels.index("C2H4(3)")
        p = rs.polymer_pools[0]
        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
        # V_poly = 1, mu(conc) = (1, 5, 30), xs = 3, k_u = 0.1:
        F = float(dn_dt[dp_idx])
        assert F > 0.0  # boundary flux is live
        # only the handshake touches mu0 (k_scission = 0)
        assert -float(dn_dt[p.mu_indices[0]]) == pytest.approx(F)
        # unzip drain r = k_u*mu0 = 0.1 stacks on the handshake's xs*F
        assert -float(dn_dt[p.mu_indices[1]]) == pytest.approx(3.0 * F + 0.1)
        # unzip mu2 drain k_u*(2*mu1 - mu0) = 0.9 stacks on xs^2*F
        assert -float(dn_dt[p.mu_indices[2]]) == pytest.approx(9.0 * F + 0.9)
        # the unzip release deposits r into the GAS monomer, not the target
        assert float(dn_dt[mono_idx]) == pytest.approx(0.1)

    def test_rejects_explicit_dp_token_without_any_block(
            self, explicit_dp_deck):
        """P2 pairing gate, direction (a): a conventions recipe_revision
        carrying an explicit-dp token claims handshake algebra the pools do
        not carry — a hand-edited artifact must fail loud (the producer
        stamps the token iff a block is present)."""
        with open(explicit_dp_deck[1]) as fh:
            artifact = json.load(fh)
        del artifact["pools"][0]["explicit_dp"]
        assert "explicit_dp" not in json.dumps(artifact["pools"])
        with pytest.raises(ValueError, match="explicit-dp|explicit_dp"):
            _build_explicit(explicit_dp_deck, artifact)

    def test_rejects_block_without_explicit_dp_token(self, explicit_dp_deck):
        """P2 pairing gate, direction (b): an explicit_dp block under a
        non-explicit-dp recipe_revision claims rates from a recipe with no
        handshake algebra — reject, never integrate permissively."""
        with open(explicit_dp_deck[1]) as fh:
            artifact = json.load(fh)
        artifact["conventions"]["recipe_revision"] = "2026-07-03-monomer-gas"
        with pytest.raises(ValueError, match="recipe_revision"):
            _build_explicit(explicit_dp_deck, artifact)


def _s2_conduit_ve_deck(tmp_path):
    """Stage-S2 conduit deck: the live PP tertiary H-abstraction row
    (H + polypropylene -> H2 + polypropylene_mod) ROUTED by the real
    generation machinery (per-product verdicts -> _handshake_structures ->
    stamp_polymer_flux_archetype), then written out as chem.yaml +
    polymer_pools.json in the proven _ve_deck idiom (plain proxy Species
    stand-ins in core; pool objects in the artifact). Returns
    (chem_path, art_path, monomer_mw, mw_h, mw_h2, eject_units)."""
    from rmgpy.data.kinetics.family import _handshake_structures
    from rmgpy.polymer import (compute_h_loss_feature_verdicts,
                               is_end_group_reaction,
                               stamp_polymer_flux_archetype)

    parent = Polymer(label="polypropylene", monomer="[CH2][CH](C)",
                     end_groups=["[H]", "[H]"], cutoff=3,
                     moments=[1.0, 40.0, 2000.0], initial_mass=0.0,
                     k_scission=0.0, k_unzip=0.0)
    h_spc = Species(label="H", molecule=[Molecule().from_smiles("[H]")])
    h2_spc = Species(label="H2", molecule=[Molecule().from_smiles("[H][H]")])
    routed = Reaction(reactants=[h_spc, parent],
                      products=[h2_spc,
                                Molecule().from_smiles("CCC[C](C)CC(C)C")],
                      reversible=False)
    polymer_reactants = [parent]
    verdicts = compute_h_loss_feature_verdicts(
        routed.reactants, routed.products, polymer_reactants)
    _handshake_structures(routed.products, polymer_reactants,
                          h_loss_verdicts=verdicts)
    routed.is_end_group_reaction = is_end_group_reaction(routed.products)
    stamp_polymer_flux_archetype(routed, routed.reactants, polymer_reactants)
    daughter = routed.products[1]
    assert isinstance(daughter, Polymer) and daughter.label == \
        "polypropylene_mod", (
        f"S2 routing must land the tertiary daughter in the feature pool, "
        f"got {getattr(daughter, 'label', type(daughter))!r}")
    assert routed.polymer_flux_archetype == int(
        PolymerFluxArchetype.VOLATILE_EJECTION)
    eject_units = float(routed.polymer_eject_units)

    # deck stand-ins (same idiom as _ve_deck): proxy Species in core, pool
    # objects (parent declared + routed born-at-zero daughter) in the artifact
    n2 = _spc("N#N", "N2", index=1)
    hatom = _spc("[H]", "H", index=2)
    h2 = _spc("[H][H]", "H2", index=3)
    pa = _spc("CCCC(C)CC(C)C", "polypropylene", index=4)
    pb = Species(molecule=[daughter.get_proxy_species().molecule[0].copy(
        deep=True)])
    pb.label = daughter.label
    pb.index = 5
    rows = [NASAPolynomial(coeffs=[2.5, 0, 0, 0, 0, -745.375, 3.35532],
                           Tmin=(tmin, "K"), Tmax=(tmax, "K"))
            for tmin, tmax in ((200.0, 1000.0), (1000.0, 6000.0))]
    pb.thermo = NASA(polynomials=rows, Tmin=(200.0, "K"), Tmax=(6000.0, "K"))
    mus = ([_mu(f"polypropylene_mu{k}") for k in range(3)]
           + [_mu(f"polypropylene_mod_mu{k}") for k in range(3)])
    core = [n2, hatom, h2, pa, pb] + mus
    pa.is_polymer_proxy = True
    pb.is_polymer_proxy = True

    ve_rxn = Reaction(
        reactants=[hatom, pa], products=[h2, pb],
        kinetics=Arrhenius(A=(1.0e3, "m^3/(mol*s)"), n=0.0,
                           Ea=(0.0, "J/mol"), T0=(1.0, "K")),
        reversible=False)
    ve_rxn.polymer_flux_archetype = routed.polymer_flux_archetype
    ve_rxn.polymer_eject_units = eject_units

    data, index_map = generate_cantera_data(core, [ve_rxn],
                                            return_reaction_index_map=True)
    chem_path = os.path.join(str(tmp_path), "chem.yaml")
    with open(chem_path, "w") as fh:
        yaml.dump(data, fh, sort_keys=False, default_flow_style=None)
    artifact = build_polymer_moments_artifact(
        [parent, daughter], core_species=core, core_reactions=[ve_rxn],
        configured_pool_labels=["polypropylene", "polypropylene_mod"],
        condensed_species=[pa, pb] + mus,
        cantera_index_map=index_map)
    art_path = os.path.join(str(tmp_path), "polymer_pools.json")
    with open(art_path, "w") as fh:
        json.dump(artifact, fh, indent=2, default=str)
    mw_h = hatom.molecule[0].get_molecular_weight() * 1000.0
    mw_h2 = h2.molecule[0].get_molecular_weight() * 1000.0
    return (chem_path, art_path, parent.monomer_mw_g_mol, mw_h, mw_h2,
            eject_units)


class TestS2ConduitVEMassInvariant:
    """S2 pin 6: the ROUTED H-abstraction conduit row conserves mass at the
    RHS level -- pool source debit + feature-pool credit + gas credit nets
    to zero including the single transferred H -- through the real emitter +
    reference-runner path (the db6f46e37 VE numeric-pin pattern)."""

    def _build(self, tmp_path):
        deck = _s2_conduit_ve_deck(tmp_path)
        chem_path, art_path = deck[0], deck[1]
        with open(art_path) as fh:
            artifact = json.load(fh)
        species, reactions = load_chem_yaml(chem_path)
        rs, core, _ = build_system_from_artifact(
            artifact, species, reactions, T0=800.0, P=1.0e5, V_poly=1.0,
            initial_moles={"N2(1)": 1.0, "H(2)": 0.5},
            mass_transfer_spec=[])
        return deck, artifact, rs, core

    def test_routed_row_emits_configured_cross_pool_ve(self, tmp_path):
        """The routed row's sidecar shape: volatile_ejection/1, src the
        parent pool, dst the (configured) feature pool, NOT unresolved,
        signed atom-transfer eject_units = +MW(H)/monomer_MW; the feature
        pool itself is spawned_empty at [0,0,0] with the parent monomer
        MW."""
        deck, artifact, _, _ = self._build(tmp_path)
        _, _, monomer_mw, mw_h, _, eject_units = deck
        assert eject_units == pytest.approx(mw_h / monomer_mw, rel=1e-9)
        (row,) = artifact["reactions"]
        assert row["archetype"] == "volatile_ejection/1"
        assert row["src_pool"] == "polypropylene"
        assert row["dst_pool"] == "polypropylene_mod"
        assert row["unresolved"] is False
        assert row["scaling"] == "mu1"
        assert row["params"] == {"eject_units": pytest.approx(eject_units)}
        assert "refused" not in row
        pools = {p["label"]: p for p in artifact["pools"]}
        mod = pools["polypropylene_mod"]
        assert mod["moments"] == [0.0, 0.0, 0.0]
        assert mod["moments_provenance"] == "spawned_empty"
        assert mod["parent_pool"] == "polypropylene"
        assert mod["monomer_mw_g_mol"] == pytest.approx(monomer_mw)

    def test_rhs_mass_conservation_including_transferred_h(self, tmp_path):
        """Numeric pin through the real solver residual: with a = MW(H)/
        monomer_MW, per event the source pool loses a full length-biased
        bundle, the feature pool gains the a-shifted bundle, and the gas
        side nets +MW(H2) - MW(H) = +MW(H) -- total mass rate exactly
        zero."""
        deck, _, rs, core = self._build(tmp_path)
        _, _, monomer_mw, mw_h, mw_h2, a = deck
        dn = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
        i = {s.label: k for k, s in enumerate(core)}
        ev = float(dn[i["H2(3)"]])
        assert ev > 0.0, "the routed VE row must carry live flux"
        assert dn[i["H(2)"]] == pytest.approx(-ev, rel=1e-12)
        # length-biased bundle off mu = [1, 40, 2000]: b0=1, b1=mu2/mu1=50
        b1 = 2000.0 / 40.0
        assert dn[i["polypropylene_mu0"]] == pytest.approx(-ev, rel=1e-9)
        assert dn[i["polypropylene_mu1"]] == pytest.approx(-ev * b1,
                                                           rel=1e-9)
        assert dn[i["polypropylene_mod_mu0"]] == pytest.approx(+ev, rel=1e-9)
        assert dn[i["polypropylene_mod_mu1"]] == pytest.approx(
            +ev * (b1 - a), rel=1e-9)
        chain_mass_rate = monomer_mw * (dn[i["polypropylene_mu1"]]
                                        + dn[i["polypropylene_mod_mu1"]])
        gas_mass_rate = mw_h2 * dn[i["H2(3)"]] + mw_h * dn[i["H(2)"]]
        assert chain_mass_rate == pytest.approx(-ev * mw_h, rel=1e-9)
        assert gas_mass_rate == pytest.approx(+ev * mw_h, rel=1e-9)
        assert chain_mass_rate + gas_mass_rate == pytest.approx(
            0.0, abs=abs(chain_mass_rate) * 1e-9)
