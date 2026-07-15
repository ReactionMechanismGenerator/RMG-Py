#!/usr/bin/env python3
"""Tests for the polymer moments CLI reference runner (design spec §8)."""

import csv
import json
import logging
import os
import sys
import time
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
    _csv_header,
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


class TestStaleTopologyRejection:
    """Round-27 P1-C enforcement: the reference runner must REJECT an
    artifact stamped conventions.stale_topology: true (its engine-derived
    surfaces describe the PRE-rebuild model and may lie about liveness)
    unless the explicit allow_stale=True / --allow-stale debug flag is
    passed, with a loud error naming the flag."""

    def test_stale_artifact_rejected_without_flag(self, deck):
        chem_path, art_path = deck
        with open(art_path) as fh:
            artifact = json.load(fh)
        artifact.setdefault("conventions", {})["stale_topology"] = True
        species, reactions = load_chem_yaml(chem_path)
        with pytest.raises(ValueError, match=r"stale_topology.*--allow-stale"
                                             r"|--allow-stale"):
            build_system_from_artifact(
                artifact, species, reactions, T0=800.0, P=1.0e5, V_poly=1.0,
                initial_moles={"N2(1)": 1.0}, mass_transfer_spec=[])

    def test_stale_artifact_accepted_with_flag(self, deck):
        chem_path, art_path = deck
        with open(art_path) as fh:
            artifact = json.load(fh)
        artifact.setdefault("conventions", {})["stale_topology"] = True
        species, reactions = load_chem_yaml(chem_path)
        rs, core, _ = build_system_from_artifact(
            artifact, species, reactions, T0=800.0, P=1.0e5, V_poly=1.0,
            initial_moles={"N2(1)": 1.0}, mass_transfer_spec=[],
            allow_stale=True)
        assert rs is not None and len(core) > 0

    def test_fresh_artifact_unaffected_by_default(self, deck):
        """No marker => the gate is a no-op (fresh artifacts load exactly
        as before, flag or no flag)."""
        chem_path, art_path = deck
        with open(art_path) as fh:
            artifact = json.load(fh)
        assert "stale_topology" not in artifact.get("conventions", {})
        species, reactions = load_chem_yaml(chem_path)
        rs, core, _ = build_system_from_artifact(
            artifact, species, reactions, T0=800.0, P=1.0e5, V_poly=1.0,
            initial_moles={"N2(1)": 1.0}, mass_transfer_spec=[])
        assert rs is not None and len(core) > 0


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
        # 3.0 is now an IMPLEMENTED major (SGH kernel-v2) and _schema_minor
        # maps it above every 2.x _MIN_SCHEMA_MINOR, so a 3.0 artifact
        # legitimately carries the 2.x QSSA vocabulary. Use an unimplemented
        # major (4.0) to keep the original intent: QSSA vocabulary cannot
        # ride a schema this loader does not implement.
        artifact = _load_artifact(qssa_deck)
        artifact["schema_version"] = "4.0"
        with pytest.raises(ValueError, match=r"schema_version"):
            _build_qssa(qssa_deck, artifact)

    def test_rejects_unknown_future_minor(self, qssa_deck, deck):
        """Weak-link milestone iv POLICY CHANGE (was minor-permissive): the
        loader pins the maximum schema minor it implements (2.7 since the
        side_group_homolysis block; was 2.6 at the homolysis_initiation
        block, 2.5 at the spawned-pool closure, 2.4 at the refused-row
        marker, 2.3 at the explicit-DP handshake block). A newer-minor
        artifact may carry vocabulary outside the channel blocks that the
        unknown-key guards never see (new conventions, new pool fields),
        so an older loader must fail loud instead of loading additively.
        Per-bump precedent: each schema bump moves this probe one minor up
        (2.8 -> 2.9 at the end-radical depropagation block; 2.7 -> 2.8 at
        the side-group block)."""
        artifact = _load_artifact(qssa_deck)
        artifact["schema_version"] = "2.9"
        with pytest.raises(ValueError, match=r"schema_version.*2\.9"):
            _build_qssa(qssa_deck, artifact)
        artifact = _load_artifact(qssa_deck)
        artifact["schema_version"] = "2.10"
        with pytest.raises(ValueError, match=r"schema_version"):
            _build_qssa(qssa_deck, artifact)
        # no-QSSA legacy artifact: same envelope pin
        chem_path, art_path = deck
        with open(art_path) as fh:
            legacy = json.load(fh)
        legacy["schema_version"] = "2.9"
        species, reactions = load_chem_yaml(chem_path)
        with pytest.raises(ValueError, match=r"schema_version.*2\.9"):
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
        restores non-accumulating, and 'qssa-unassessable' (round-20
        increment 7: an accumulating refusal whose radical/consumers were
        not visible at the rebuild census) restores accumulating. The
        reason attr itself is restored too (single-source stamp: a re-emit
        reads Reaction.polymer_refused_reason and must round-trip
        byte-identically). (A 'bogus' reason never reaches this coercion
        -- test_rejects_unknown_refused_reason pins the earlier
        rejection.)"""
        for reason, accumulating in (("qssa-invalid", True),
                                     ("conduit-deferred", False),
                                     ("qssa-unassessable", True)):
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
            assert rxn.polymer_refused_reason == reason, (
                f"single-source stamp: the restored reason attr must "
                f"round-trip {reason!r}")


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

          length-biased bundle: b0 = 1, b1 = mu2/mu1 = 12,
              mu3 = mu0*(mu2/mu1)^3 = 1728 (log-Lagrange closure),
              b2 = mu3/mu1 = 172.8
          event rate: the pool is HEALTHY (mu0 = 1 mol sits ~1e14 r81
              floors above exhaustion), so the near-exhaustion bundle
              limiter (tail-only smoothstep, round-27 P1-A) is exactly
              inactive and the row runs at the plain adjudicated site law
              S_base = mu1/V_poly: ev = kf * mu1 = 2.5 * 10 (the P1-1
              bundle cap min(10, 1, 10/12, 120/172.8) applies only inside
              the depletion band)
          src debit (full bundle):  -ev*(b0, b1, b2)
          dst credit (a-shifted, sa = -a): +ev*(b0, b1 - a,
              b2 - 2*a*b1 + a^2)
          gas volatile (standard net-rate product path, the SAME rate --
              gas flux and moment flux never diverge): +ev
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

            PR + C2H5 <=> poly        (C22H45 + C2H5 -> C24H50, balanced)

        has the POOL PROXY itself as a product, so the artifact compiler
        (rmgpy/polymer.py compile_polymer_reaction_entries) emits an entry,
        and the runner's _restamp_and_extend tag restoration + window plumbing
        are both LIVE on this deck:

        * PR (C22H45, 309.6 g/mol) is a GAS chain-scale counterparty (clears
          the r95 absolute floor): it is in the entry's reactants but NOT in
          proxy_reactants/proxy_products, so its physically-melt membership in
          the consumer world exists ONLY via the runner's blanket tag
          restoration. Paired against the condensed proxy (338.6 g/mol) it
          keeps U ~ 0.06 decades -> SILENT.
        * C2H5 (29.1 g/mol) sits in the bare-slack..window MW band: above the
          10 g/mol bare slack the chain window degrades to if the
          monomer_mw_g_mol plumb is deleted, below the real window
          (pool monomer C2H4 28.05 + 10 = 38.05 g/mol). With the plumb intact
          the tag-branch MW conjunct excludes it; without it, C2H5 enters the
          melt sum and U blows up to ~10 decades -> refusal."""
        n2 = _spc("N#N", "N2", index=1)
        et = _spc("C[CH2]", "C2H5")              # ethyl, 29.06 g/mol: band species
        # C22 GAS chain radical (309.60 g/mol / 22 heavy): a genuine
        # multi-repeat chain clearing the r95 absolute floor, so its restored
        # tag makes it melt. PR + C2H5 -> C24H50 proxy (balanced).
        pr = _spc("[CH2]CCCCCCCCCCCCCCCCCCCCC", "PR")   # C22H45, GAS chain radical
        proxy = _spc("CCCCCCCCCCCCCCCCCCCCCCCC", "poly")  # C24H50 pool-proxy species
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
        parent pool, dst the engine-configured feature pool, NOT
        unresolved, signed atom-transfer eject_units = +MW(H)/monomer_MW;
        the feature pool itself is spawned_empty at [0,0,0] with the
        parent monomer MW -- and (item-16 P1 split) it is PUBLISHED via
        conventions.spawned_pools, never configured_pools."""
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
        # Item-16 P1 split: the engine-configured daughter is published
        # through conventions.spawned_pools, NEVER configured_pools (the
        # CKMG consumer hard-flags dead configured pools), while the row
        # above stays live -- build_system_from_artifact still builds the
        # pool (live-coupled, mu triplet mechanism-resident), which the
        # RHS mass-conservation pin below exercises numerically.
        conv = artifact["conventions"]
        assert "polypropylene_mod" not in conv["configured_pools"]
        assert conv["spawned_pools"] == ["polypropylene_mod"]
        assert artifact["schema_version"] == "2.5"

    def test_rhs_mass_conservation_including_transferred_h(self, tmp_path):
        """Numeric pin through the real solver residual, defect-aware
        (P1-2): a = MW(H)/monomer_MW is a MASS defect, not a DP decrement.
        The H-loss feature pool carries chain_mass_defect_g_mol = MW(H)
        (emitted by the real serializer, restored by the runner), so the
        row's MOMENT shift is a_moment = a - defect/monomer_MW = 0: per
        event the source pool loses a full length-biased bundle and the
        feature pool gains the SAME bundle UNSHIFTED (booking the a
        decrement landed monomeric chains below the realizability cone
        mu1 >= mu0 by construction -- the regen-#2 seed defect). The
        condensed-mass closure books the loss through the defect ledger
        (mu1*MW - mu0*defect), and the gas side nets
        +MW(H2) - MW(H) = +MW(H) -- total mass rate exactly zero."""
        deck, artifact, rs, core = self._build(tmp_path)
        _, _, monomer_mw, mw_h, mw_h2, a = deck
        # the serialized feature pool carries the exact per-chain defect
        (mod_pool,) = [p for p in artifact["pools"]
                       if p["label"] == "polypropylene_mod"]
        defect = mod_pool["chain_mass_defect_g_mol"]
        assert defect == pytest.approx(mw_h, rel=1e-9)
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
        # chains land with UNCHANGED mu1 (a_moment == 0), NOT b1 - a
        assert dn[i["polypropylene_mod_mu1"]] == pytest.approx(
            +ev * b1, rel=1e-9)
        assert dn[i["polypropylene_mod_mu1"]] != pytest.approx(
            +ev * (b1 - a), rel=1e-9)
        chain_mass_rate = (monomer_mw * (dn[i["polypropylene_mu1"]]
                                         + dn[i["polypropylene_mod_mu1"]])
                           - defect * dn[i["polypropylene_mod_mu0"]])
        gas_mass_rate = mw_h2 * dn[i["H2(3)"]] + mw_h * dn[i["H(2)"]]
        assert chain_mass_rate == pytest.approx(-ev * mw_h, rel=1e-9)
        assert gas_mass_rate == pytest.approx(+ev * mw_h, rel=1e-9)
        assert chain_mass_rate + gas_mass_rate == pytest.approx(
            0.0, abs=abs(chain_mass_rate) * 1e-9)

    def test_replay_csv_reports_buildable_spawned_pool_state(self, tmp_path):
        """Item-16 P1 replay-output pin: the CSV surface must not hide
        integrated spawned-pool state. The live-VE artifact's buildable
        spawned pool (polypropylene_mod: live cross-pool row + mechanism-
        resident mu triplet) gets its own mu0/mu1/mu2 columns AFTER the
        configured pool's (same <label>_muN_mol naming), the header matches
        the run_segments row shape exactly, and the emitted values ARE the
        integrated solver state (which the VE conduit feeds, so mu0 is
        live-nonzero, never a padded 0 column)."""
        deck = _s2_conduit_ve_deck(tmp_path)
        chem_path, art_path = deck[0], deck[1]
        with open(art_path) as fh:
            artifact = json.load(fh)
        species, reactions = load_chem_yaml(chem_path)
        rs, core, all_rxns = build_system_from_artifact(
            artifact, species, reactions, T0=800.0, P=1.0e5, V_poly=1.0,
            initial_moles={"N2(1)": 1.0, "H(2)": 0.5},
            mass_transfer_spec=[])
        header = _csv_header(artifact, core)
        assert header[:2] == ["t_s", "T_K"]
        assert header[2:5] == ["polypropylene_mu0_mol",
                               "polypropylene_mu1_mol",
                               "polypropylene_mu2_mol"]
        assert header[5:8] == ["polypropylene_mod_mu0_mol",
                               "polypropylene_mod_mu1_mol",
                               "polypropylene_mod_mu2_mol"]
        # t_end pinned at 1 ns: the born-empty _mod pool makes DASPK's
        # first-step heuristic stall for larger first targets (pre-existing
        # integrator quirk, unrelated to the CSV surface); ~2e-4 mol of VE
        # flux lands in the daughter by then -- ample to pin live columns.
        rows = run_segments(rs, core, artifact, all_rxns,
                            [(1.0e-9, 800.0)], n_points_per_segment=3)
        last = rows[-1]
        assert len(last) == len(header)
        i = {s.label: k for k, s in enumerate(core)}
        y = np.asarray(rs.y)
        for k in range(3):
            col = header.index(f"polypropylene_mod_mu{k}_mol")
            assert last[col] == float(y[i[f"polypropylene_mod_mu{k}"]])
        assert last[header.index("polypropylene_mod_mu0_mol")] > 0.0


# ---------------------------------------------------------------------------
# Radical-homolysis initiation block, consumer side (schema 2.6, Stage 2 of
# the adjudicated round-66/67 arc): strict block guard (key presence + shape,
# never truthiness), version cross-check, closure membership (daughters in
# pools[] / classified / condensed), kernel + recipe pins, and the
# machine-checkable kernel-supersession census (refused rows stay refused
# and carry ZERO flux while the kernel is live).
# ---------------------------------------------------------------------------


def _homolysis_deck(tmp_path, with_refused=False):
    """A kernel-enabled parent pool + its two eagerly-configured end-radical
    daughter pools (the live registration shape: polymer.pyx
    _flatten_homolysis_state hard-errors on unconfigured daughters), written
    out as chem.yaml + polymer_pools.json exactly like an RMG run would.
    ``with_refused=True`` adds one pool-mapped conduit-deferred refused row
    (the supersession pairing of round-66/67 ruling (c))."""
    n2 = _spc("N#N", "N2", index=1)
    pool = Polymer(label="PP", monomer="[CH2][CH](C)",
                   end_groups=["[H]", "[H]"], cutoff=3,
                   moments=[1.0, 50.0, 3000.0], initial_mass=0.0,
                   k_homolysis={"A": 1.0e13, "n": 0.5, "Ea": 1.2e5})
    prim, sec = pool.generate_end_radical_daughters()
    mus = []
    for base in ("PP", prim.label, sec.label):
        mus += [_mu(f"{base}_mu{k}") for k in range(3)]
    core = [n2]
    rxns = []
    condensed = list(mus)
    if with_refused:
        proxy = _spc("CCC(C)CC(C)C", "PP", index=2)
        rad = _spc("[CH2]C(C)CC(C)C", "R", index=3)
        h = _spc("[H]", "H", index=4)
        h2 = _spc("[H][H]", "H2", index=5)
        proxy.is_polymer_proxy = True
        refused_rxn = Reaction(
            reactants=[proxy, h], products=[rad, h2],
            kinetics=Arrhenius(A=(2.0, "m^3/(mol*s)"), n=0.0,
                               Ea=(0.0, "J/mol"), T0=(1.0, "K")),
            reversible=False)
        refused_rxn.polymer_flux_archetype = int(
            PolymerFluxArchetype.UNRESOLVED)
        refused_rxn.polymer_refused = True
        refused_rxn.polymer_refused_accumulating = False
        core += [proxy, rad, h, h2]
        rxns = [refused_rxn]
        condensed = [proxy, rad] + mus
    core += mus
    data, index_map = generate_cantera_data(core, rxns,
                                            return_reaction_index_map=True)
    chem_path = os.path.join(str(tmp_path), "chem.yaml")
    with open(chem_path, "w") as fh:
        yaml.dump(data, fh, sort_keys=False, default_flow_style=None)
    artifact = build_polymer_moments_artifact(
        [pool, prim, sec], core_species=core, core_reactions=rxns,
        configured_pool_labels=["PP", prim.label, sec.label],
        condensed_species=condensed, cantera_index_map=index_map)
    art_path = os.path.join(str(tmp_path), "polymer_pools.json")
    with open(art_path, "w") as fh:
        json.dump(artifact, fh, indent=2, default=str)
    return chem_path, art_path


@pytest.fixture
def homolysis_deck(tmp_path):
    return _homolysis_deck(tmp_path)


def _build_homolysis(deck, artifact=None, initial_moles=None):
    chem_path, art_path = deck
    if artifact is None:
        with open(art_path) as fh:
            artifact = json.load(fh)
    species, reactions = load_chem_yaml(chem_path)
    return build_system_from_artifact(
        artifact, species, reactions, T0=800.0, P=1.0e5, V_poly=1.0,
        initial_moles=initial_moles or {"N2(1)": 1.0},
        mass_transfer_spec=[])


def _homolysis_artifact(deck):
    with open(deck[1]) as fh:
        return json.load(fh)


class TestHomolysisInitiationConsumer:
    """Schema 2.6 homolysis_initiation contract, consumer side."""

    def test_round_trip_kernel_flattened_and_daughters_resolved(
            self, homolysis_deck):
        """RED pin (full round trip): a real kernel-enabled PolymerPhase
        pool serialized by the emitter loads GREEN, wires the kernel into
        the reconstructed oracle (khom_* flattened arrays), and resolves
        both end-radical daughter pools' moment slots."""
        artifact = _homolysis_artifact(homolysis_deck)
        assert artifact["schema_version"] == "2.6"
        entry = next(p for p in artifact["pools"] if p["label"] == "PP")
        block = entry.get("homolysis_initiation")
        assert isinstance(block, dict)
        assert block["kernel"] == "radical_homolysis_initiation/1"
        assert block["recipe_revision"] == "2026-07-05-radical-homolysis"

        rs, core, _ = _build_homolysis(homolysis_deck, artifact)
        assert len(rs.polymer_pools) == 3
        assert rs.khom_enabled[0] == 1
        assert rs.khom_A[0] == pytest.approx(1.0e13)
        assert rs.khom_n[0] == pytest.approx(0.5)
        assert rs.khom_Ea[0] == pytest.approx(1.2e5)
        assert rs.khom_prim_mu0[0] >= 0
        assert rs.khom_sec_mu0[0] >= 0
        labels = [s.label for s in core]
        assert labels[rs.khom_prim_mu0[0]] == "PP_rad_primary_end_mu0"
        assert labels[rs.khom_sec_mu0[0]] == "PP_rad_secondary_end_mu0"

    def test_accepts_2_6_without_homolysis(self, deck):
        """2.6 is now an implemented minor: a 2.6 stamp with no
        homolysis_initiation block anywhere loads (mirror of
        test_accepts_2_5_without_spawned)."""
        chem_path, art_path = deck
        with open(art_path) as fh:
            artifact = json.load(fh)
        artifact["schema_version"] = "2.6"
        species, reactions = load_chem_yaml(chem_path)
        rs, _, _ = build_system_from_artifact(
            artifact, species, reactions, T0=800.0, P=1.0e5, V_poly=1.0,
            initial_moles={"N2(1)": 1.0}, mass_transfer_spec=[])
        assert len(rs.polymer_pools) == 1
        assert rs.khom_enabled[0] == 0

    def test_rejects_block_in_2_5_stamped_artifact(self, homolysis_deck):
        """Vocabulary/version cross-check (the 2.4-refused / 2.5-spawned
        precedent): a below-2.6 artifact carrying a homolysis_initiation
        block is malformed -- the emitter stamps 2.6 whenever it writes
        one."""
        artifact = _homolysis_artifact(homolysis_deck)
        artifact["schema_version"] = "2.5"
        with pytest.raises(ValueError, match=r"homolysis_initiation.*2\.6"):
            _build_homolysis(homolysis_deck, artifact)

    def _shape_case(self, deck, mutate, match):
        """Load a fresh homolysis artifact, apply ``mutate(artifact)``, and
        assert the loader rejects with an actionable error."""
        artifact = _homolysis_artifact(deck)
        mutate(artifact)
        with pytest.raises(ValueError, match=match):
            _build_homolysis(deck, artifact)

    @staticmethod
    def _block(artifact):
        return next(p for p in artifact["pools"]
                    if p["label"] == "PP")["homolysis_initiation"]

    def test_rejects_daughter_missing_from_pools(self, homolysis_deck):
        def cut(a):
            a["pools"] = [p for p in a["pools"]
                          if p["label"] != "PP_rad_primary_end"]
        self._shape_case(homolysis_deck, cut,
                         r"PP_rad_primary_end.*pools")

    def test_rejects_unclassified_daughter(self, homolysis_deck):
        """A daughter in NEITHER conventions.configured_pools NOR
        conventions.spawned_pools is unclassified in the configured/spawned
        closure -- its moment credits would have no solver home."""
        def cut(a):
            a["conventions"]["configured_pools"] = [
                lbl for lbl in a["conventions"]["configured_pools"]
                if lbl != "PP_rad_secondary_end"]
        self._shape_case(homolysis_deck, cut,
                         r"PP_rad_secondary_end.*(configured|spawned|classif)")

    def test_rejects_spawned_only_daughter(self, homolysis_deck):
        """r68 P1: a daughter classified ONLY in conventions.spawned_pools
        must be rejected BEFORE solver construction --
        build_system_from_artifact constructs PolymerPoolConfigs from
        configured_pools alone, so a spawned-only daughter is never built
        and the kernel's fragment credits have no solver home. Schema 2.5
        defines spawned_pools as the configured set's complement, so
        membership there contradicts the eager-configured daughter
        design."""
        def cut(a):
            conv = a["conventions"]
            conv["configured_pools"] = [
                lbl for lbl in conv["configured_pools"]
                if lbl != "PP_rad_primary_end"]
            conv["spawned_pools"] = ["PP_rad_primary_end"]
        self._shape_case(homolysis_deck, cut,
                         r"PP_rad_primary_end.*spawned_pools")

    @pytest.mark.parametrize("also_spawn", [False, True],
                             ids=["dropped", "respawned"])
    def test_rejects_unconfigured_carrier(self, homolysis_deck, also_spawn):
        """r68 P1: the pool CARRYING the homolysis_initiation block must
        itself be in conventions.configured_pools --
        build_system_from_artifact skips unconfigured pools, so a block on
        an unconfigured carrier would be a silently dropped kernel."""
        def cut(a):
            conv = a["conventions"]
            conv["configured_pools"] = [
                lbl for lbl in conv["configured_pools"] if lbl != "PP"]
            if also_spawn:
                conv["spawned_pools"] = ["PP"]
        self._shape_case(homolysis_deck, cut,
                         r"'PP'.*configured_pools")

    def test_rejects_uncondensed_daughter(self, homolysis_deck):
        def cut(a):
            a["conventions"]["condensed_species"] = [
                lbl for lbl in a["conventions"]["condensed_species"]
                if not lbl.startswith("PP_rad_primary_end")]
        self._shape_case(homolysis_deck, cut,
                         r"PP_rad_primary_end.*condensed")

    def test_rejects_missing_kinetics_key(self, homolysis_deck):
        def cut(a):
            del self._block(a)["kinetics"]["Ea"]
        self._shape_case(homolysis_deck, cut, r"kinetics")

    def test_rejects_wrong_kinetics_units(self, homolysis_deck):
        def cut(a):
            self._block(a)["kinetics"]["units"]["A"] = "1/min"
        self._shape_case(homolysis_deck, cut, r"units")

    def test_rejects_non_positive_A(self, homolysis_deck):
        def cut(a):
            self._block(a)["kinetics"]["A"] = 0.0
        self._shape_case(homolysis_deck, cut, r"A")

    def test_rejects_boolean_kinetics_value(self, homolysis_deck):
        """Key-presence + SHAPE validation, not truthiness (the r63 lesson):
        JSON true is not a rate constant."""
        def cut(a):
            self._block(a)["kinetics"]["A"] = True
        self._shape_case(homolysis_deck, cut, r"A")

    def test_rejects_present_disabled_block(self, homolysis_deck):
        def cut(a):
            self._block(a)["enabled"] = False
        self._shape_case(homolysis_deck, cut, r"present-disabled|enabled")

    def test_rejects_unknown_block_key(self, homolysis_deck):
        def cut(a):
            self._block(a)["gas_release"] = "CH4"
        self._shape_case(homolysis_deck, cut, r"unknown key")

    def test_rejects_missing_kernel_field(self, homolysis_deck):
        def cut(a):
            del self._block(a)["kernel"]
        self._shape_case(homolysis_deck, cut, r"kernel|missing")

    def test_rejects_unknown_kernel_name(self, homolysis_deck):
        """Ruling (c): the kernel field is the machine-checkable
        supersession contract -- an unknown kernel is flux this consumer
        cannot reproduce."""
        def cut(a):
            self._block(a)["kernel"] = "beta_scission_zip/9"
        self._shape_case(homolysis_deck, cut, r"kernel")

    def test_rejects_missing_recipe_revision(self, homolysis_deck):
        def cut(a):
            del self._block(a)["recipe_revision"]
        self._shape_case(homolysis_deck, cut, r"recipe_revision|missing")

    def test_rejects_unknown_recipe_revision(self, homolysis_deck):
        def cut(a):
            self._block(a)["recipe_revision"] = "2026-01-01-bogus"
        self._shape_case(homolysis_deck, cut, r"recipe_revision")

    def test_rejects_tampered_recipe(self, homolysis_deck):
        def cut(a):
            block = self._block(a)
            key = sorted(block["recipe"])[0]
            block["recipe"][key] = "dmu1 += R (mass fabricated)"
        self._shape_case(homolysis_deck, cut, r"recipe")

    def test_rejects_daughter_without_spawn_provenance(self, homolysis_deck):
        """The daughters are producer-spawned pools; a daughter row claiming
        user/input provenance is a hand-edited artifact."""
        def cut(a):
            entry = next(p for p in a["pools"]
                         if p["label"] == "PP_rad_primary_end")
            entry["spawn_event_metadata"] = {"source": "input"}
        self._shape_case(homolysis_deck, cut, r"provenance|spawn")

    def test_supersession_census_refused_rows_zero_flux(self, tmp_path):
        """Ruling (c), machine-checkable: with the kernel live, refused
        conduit-deferred rows stay marked refused on the reconstructed
        reactions, the rebuilt oracle's homolysis supersession census names
        them, and they contribute EXACTLY zero flux -- while the kernel
        itself credits the daughter pools (so the zero is suppression, not
        a dead system)."""
        deck = _homolysis_deck(tmp_path, with_refused=True)
        artifact = _homolysis_artifact(deck)
        assert artifact["schema_version"] == "2.6"
        marked = [e for e in artifact["reactions"] if e.get("refused")]
        assert len(marked) == 1
        assert marked[0]["refused_reason"] == "conduit-deferred"

        rs, core, _ = _build_homolysis(
            deck, artifact, initial_moles={"N2(1)": 1.0, "H(4)": 1.0})
        # refused marker survived the boundary
        assert int(rs.reaction_refused.sum()) == 1
        # machine-checkable census: the kernel-enabled pool names the
        # superseded refused rows
        census = list(rs.homolysis_supersession_census)
        assert census and census[0]["pool"] == "PP"
        assert census[0]["superseded_rows"]
        # zero flux from the refused row: its non-pool participants gain
        # nothing...
        labels = [s.label for s in core]
        i = {lab: k for k, lab in enumerate(labels)}
        dn = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
        assert dn[i["R(3)"]] == 0.0
        assert dn[i["H2(5)"]] == 0.0
        # ...while the live kernel credits both end-radical daughters
        assert dn[i["PP_rad_primary_end_mu0"]] > 0.0
        assert dn[i["PP_rad_secondary_end_mu0"]] > 0.0


# ---------------------------------------------------------------------------
# Side-group homolysis block, consumer side (schema 2.7, FR1-K2, adjudicated
# rounds 72/73): strict block guard (key presence + shape, never truthiness),
# vocabulary/version cross-check for BOTH the pool-level block and the X-loss
# feature-pool chain_mass_defect_g_mol field, selector/feature closure
# validated FROM SERIALIZED DATA ALONE (site_atom_indices -- the loader has
# no monomer graph to re-derive from), kernel + recipe + NORMATIVE mass
# formula pins, round-trip wiring into the reconstructed PolymerPoolConfig,
# and the machine-checkable kernel-supersession census (refused rows stay
# refused and carry ZERO flux while the kernel is live).
# ---------------------------------------------------------------------------


def _side_group_deck(tmp_path, with_refused=False):
    """A single SGH kernel-v2 carrier pool (PVBr, side_group_homolysis/2),
    written out as chem.yaml + polymer_pools.json exactly like an RMG run
    would. SGH kernel-v2 (schema 3.0) spawns NO X-loss feature pool --
    multi-loss X depletion is tracked by an auxiliary per-(pool, channel) Z
    inventory on the carrier -- so there is no daughter pool to configure.
    ``with_refused=True`` adds one pool-mapped conduit-deferred refused row
    (the FR1-K1 supersession pairing)."""
    n2 = _spc("N#N", "N2", index=1)
    br = _spc("[Br]", "Br", index=2)
    pool = Polymer(label="PVBr", monomer="[CH2][CH]Br",
                   end_groups=["[H]", "[H]"], cutoff=3,
                   moments=[1.0, 50.0, 3000.0], initial_mass=0.0,
                   side_group_homolysis=[dict(
                       label="aliphatic_C-Br", A=1.0e13, n=0.5, Ea=1.2e5,
                       site_selector="aliphatic", sites_per_unit=1.0,
                       gas_product="[Br]")])
    pool.side_group_gas_species = [br]
    mus = [_mu(f"PVBr_mu{k}") for k in range(3)]
    core = [n2, br]
    rxns = []
    condensed = list(mus)
    if with_refused:
        proxy = _spc("CC(Br)CC(Br)C", "PVBr", index=3)
        radp = _spc("C[CH]CC(Br)C", "RP", index=4)
        hbr = _spc("Br", "HBr", index=5)
        proxy.is_polymer_proxy = True
        refused_rxn = Reaction(
            reactants=[proxy, br], products=[radp, hbr],
            kinetics=Arrhenius(A=(2.0, "m^3/(mol*s)"), n=0.0,
                               Ea=(0.0, "J/mol"), T0=(1.0, "K")),
            reversible=False)
        refused_rxn.polymer_flux_archetype = int(
            PolymerFluxArchetype.UNRESOLVED)
        refused_rxn.polymer_refused = True
        refused_rxn.polymer_refused_accumulating = False
        core += [proxy, radp, hbr]
        rxns = [refused_rxn]
        condensed = [proxy, radp] + mus
    core += mus
    data, index_map = generate_cantera_data(core, rxns,
                                            return_reaction_index_map=True)
    chem_path = os.path.join(str(tmp_path), "chem.yaml")
    with open(chem_path, "w") as fh:
        yaml.dump(data, fh, sort_keys=False, default_flow_style=None)
    artifact = build_polymer_moments_artifact(
        [pool], core_species=core, core_reactions=rxns,
        configured_pool_labels=["PVBr"],
        condensed_species=condensed, cantera_index_map=index_map)
    art_path = os.path.join(str(tmp_path), "polymer_pools.json")
    with open(art_path, "w") as fh:
        json.dump(artifact, fh, indent=2, default=str)
    return chem_path, art_path


@pytest.fixture
def side_group_deck(tmp_path):
    return _side_group_deck(tmp_path)


def _build_side_group(deck, artifact=None, initial_moles=None):
    chem_path, art_path = deck
    if artifact is None:
        with open(art_path) as fh:
            artifact = json.load(fh)
    species, reactions = load_chem_yaml(chem_path)
    return build_system_from_artifact(
        artifact, species, reactions, T0=800.0, P=1.0e5, V_poly=1.0,
        initial_moles=initial_moles or {"N2(1)": 1.0},
        mass_transfer_spec=[])


def _side_group_artifact(deck):
    with open(deck[1]) as fh:
        return json.load(fh)


class TestSideGroupHomolysisConsumer:
    """Consumer/loader side of the SGH kernel-v2 contract (schema 3.0,
    side_group_homolysis/2, explicit Br-inventory depletion). v2 DELETES
    the v1 feature/daughter pool: the block carries no per-channel
    feature_pool and the carrier no chain_mass_defect_g_mol; a legacy /1
    block is hard-rejected LOUDLY (P1-A firewall), never reinterpreted."""

    def test_accepts_side_group_homolysis_v2_artifact(self, tmp_path):
        """GOLDEN: the emitter's SGH kernel-v2 artifact (schema 3.0,
        side_group_homolysis/2) loads GREEN and the rebuilt oracle carries
        the kernel -- sgh_* row arrays populated from the block, the gas
        credit routed to the chem.yaml Br species. v2 spawns NO feature
        pool: the block carries no per-channel feature_pool and the carrier
        no chain_mass_defect_g_mol."""
        deck = _side_group_deck(tmp_path)
        artifact = _side_group_artifact(deck)
        assert artifact["schema_version"] == "3.0"
        entry = next(p for p in artifact["pools"] if p["label"] == "PVBr")
        block = entry.get("side_group_homolysis")
        assert block is not None
        assert block["kernel"] == "side_group_homolysis/2"
        assert block["recipe_revision"] == \
            "2026-07-08-side-group-inventory-depletion"
        assert block["recipe"]["mass"] == \
            "mu1*MW - sum_c max(0, sites_c*mu1 - Z_c)*M_X_c"
        assert block["recipe"]["state"] == (
            "Z_c(0)=sites_c*mu1(0); dZ_c/dt=-k*Z_c/V_poly (conc basis); "
            "pool moments unchanged")
        # v2 spawns no feature pool: no per-channel feature_pool key, and
        # the carrier carries no X-loss chain_mass_defect_g_mol.
        assert all("feature_pool" not in ch for ch in block["channels"])
        assert "chain_mass_defect_g_mol" not in entry
        rs, core, _ = _build_side_group(deck, artifact)
        assert len(rs.polymer_pools) == 1
        assert rs.sgh_enabled[0] == 1
        assert rs.sgh_A[0] == pytest.approx(1.0e13)
        assert rs.sgh_n[0] == pytest.approx(0.5)
        assert rs.sgh_Ea[0] == pytest.approx(1.2e5)
        assert rs.sgh_sites[0] == pytest.approx(1.0)
        assert rs.sgh_mw_X[0] == pytest.approx(79.904, rel=1e-3)
        labels = [s.label for s in core]
        assert labels[rs.sgh_gas[0]] == "Br(2)"

    def test_rejects_legacy_side_group_homolysis_v1_artifact(
            self, side_group_deck):
        """P1-A firewall: an OLD side_group_homolysis/1 block must be HARD-
        REJECTED loudly, never silently reinterpreted as /2 (the state-
        vector shape and mass contract both changed). Overwrite a valid v2
        block to the v1 kernel + v1 recipe_revision + a per-channel
        feature_pool key (each of which independently trips the legacy
        guard)."""
        artifact = _side_group_artifact(side_group_deck)
        block = self._block(artifact)
        block["kernel"] = "side_group_homolysis/1"
        block["recipe_revision"] = "2026-07-06-side-group-homolysis"
        block["channels"][0]["feature_pool"] = \
            "PVBr_sidegrp_aliphatic_C_Br"
        with pytest.raises(ValueError, match=r"side_group_homolysis/1"):
            _build_side_group(side_group_deck, artifact)

    def test_rejects_v1_feature_pool_spawn_provenance(self, side_group_deck):
        """A pool carrying the v1 X-loss feature-pool spawn provenance
        (spawn_event_metadata.source == 'side_group_homolysis') is a legacy
        feature-pool artifact -- v2 spawns none -- and is hard-rejected."""
        from rmgpy.tools.polymer_moments_runner import (
            _check_side_group_homolysis)
        artifact = _side_group_artifact(side_group_deck)
        artifact["pools"].append({
            "label": "PVBr_sidegrp_aliphatic_C_Br",
            "spawn_event_metadata": {
                "source": "side_group_homolysis",
                "channel": "aliphatic_C-Br"}})
        with pytest.raises(ValueError, match=r"spawns NO feature pool"):
            _check_side_group_homolysis(artifact)

    def test_rejects_explicit_dp_on_side_group_carrier(self, side_group_deck):
        """Gap D (P1-B, artifact boundary): an SGH kernel-v2 carrier that also
        carries a pool-level explicit_dp block is a silent mass-corruption hole
        -- the condensed-mass tally would mix tail-only Z with explicit-
        augmented mu1 and corrupt Br mass/speciation -- and must be HARD-
        REJECTED. The producer never serializes both on one pool (the deck and
        solver-build both refuse it); a v2 carrier carrying explicit_dp is
        hand-edited/corrupted."""
        from rmgpy.tools.polymer_moments_runner import (
            _check_side_group_homolysis)
        artifact = _side_group_artifact(side_group_deck)
        entry = next(p for p in artifact["pools"] if p["label"] == "PVBr")
        entry["explicit_dp"] = {"enabled": True}
        with pytest.raises(ValueError,
                           match=r"PVBr.*side_group_homolysis.*explicit_dp"):
            _check_side_group_homolysis(artifact)

    def test_accepts_2_7_without_side_group(self, deck):
        """A 2.7 stamp alone (no side-group vocabulary) loads fine: the
        stamp gates vocabulary, absence keeps legacy behavior."""
        chem_path, art_path = deck
        with open(art_path) as fh:
            artifact = json.load(fh)
        artifact["schema_version"] = "2.7"
        species, reactions = load_chem_yaml(chem_path)
        rs, _, _ = build_system_from_artifact(
            artifact, species, reactions, T0=800.0, P=1.0e5, V_poly=1.0,
            initial_moles={"N2(1)": 1.0}, mass_transfer_spec=[])
        assert len(rs.polymer_pools) == 1
        assert rs.sgh_enabled[0] == 0

    def test_rejects_block_in_2x_stamped_artifact(self, side_group_deck):
        """Vocabulary/version cross-check: the SGH kernel-v2 block is the
        MAJOR bump schema 3.0 (state-vector shape + mass contract both
        change), so a v2 block under any 2.x stamp is malformed -- the
        emitter stamps 3.0 whenever it writes one."""
        artifact = _side_group_artifact(side_group_deck)
        artifact["schema_version"] = "2.6"
        with pytest.raises(ValueError,
                           match=r"side_group_homolysis.*3\.0"):
            _build_side_group(side_group_deck, artifact)

    def _shape_case(self, deck, mutate, match):
        """Load the side-group artifact, apply ``mutate(artifact)``, and
        require the loader to reject with an actionable error."""
        artifact = _side_group_artifact(deck)
        mutate(artifact)
        with pytest.raises(ValueError, match=match):
            _build_side_group(deck, artifact)

    @staticmethod
    def _block(artifact):
        return next(p for p in artifact["pools"]
                    if p["label"] == "PVBr")["side_group_homolysis"]

    @classmethod
    def _channel(cls, artifact):
        return cls._block(artifact)["channels"][0]

    def test_rejects_missing_kinetics_key(self, side_group_deck):
        def cut(a):
            del self._channel(a)["kinetics"]["Ea"]
        self._shape_case(side_group_deck, cut, r"kinetics")

    def test_rejects_wrong_kinetics_units(self, side_group_deck):
        """The per-site A unit is pinned exactly: the sibling kernel's
        plain 's^-1' claims a DIFFERENT rate law and must reject."""
        def cut(a):
            self._channel(a)["kinetics"]["units"]["A"] = "s^-1"
        self._shape_case(side_group_deck, cut, r"units")

    def test_rejects_non_positive_A(self, side_group_deck):
        def cut(a):
            self._channel(a)["kinetics"]["A"] = 0.0
        self._shape_case(side_group_deck, cut, r"A")

    def test_rejects_boolean_kinetics_value(self, side_group_deck):
        """Key-presence + SHAPE validation, not truthiness: JSON true is
        not a rate constant."""
        def cut(a):
            self._channel(a)["kinetics"]["A"] = True
        self._shape_case(side_group_deck, cut, r"finite number")

    def test_rejects_present_disabled_block(self, side_group_deck):
        def cut(a):
            self._block(a)["enabled"] = False
        self._shape_case(side_group_deck, cut, r"present-disabled|ABSENT")

    def test_rejects_unknown_block_key(self, side_group_deck):
        def cut(a):
            self._block(a)["gas_release"] = "CH4"
        self._shape_case(side_group_deck, cut, r"unknown key")

    def test_rejects_missing_block_key(self, side_group_deck):
        def cut(a):
            del self._block(a)["kernel"]
        self._shape_case(side_group_deck, cut, r"missing key")

    def test_rejects_unknown_kernel(self, side_group_deck):
        """An unknown kernel is flux this consumer cannot reproduce."""
        def cut(a):
            self._block(a)["kernel"] = "beta_scission_zip/9"
        self._shape_case(side_group_deck, cut, r"kernel")

    def test_rejects_unknown_recipe_revision(self, side_group_deck):
        def cut(a):
            self._block(a)["recipe_revision"] = "2026-01-01-bogus"
        self._shape_case(side_group_deck, cut, r"recipe_revision")

    def test_rejects_tampered_recipe(self, side_group_deck):
        """The recipe strings -- including the NORMATIVE mass formula --
        are pinned independently of the emitter, exact match."""
        def cut(a):
            self._block(a)["recipe"]["mass"] = "condensed_mass_g = mu1*MW"
        self._shape_case(side_group_deck, cut, r"mass")

    def test_rejects_unknown_channel_key(self, side_group_deck):
        def cut(a):
            self._channel(a)["site_smarts"] = "[Br]"
        self._shape_case(side_group_deck, cut, r"exactly the keys")

    def test_rejects_missing_channel_key(self, side_group_deck):
        """site_selector is REQUIRED (round-72 P1: the kinetics label
        alone must never pick the site)."""
        def cut(a):
            del self._channel(a)["site_selector"]
        self._shape_case(side_group_deck, cut, r"exactly the keys")

    def test_rejects_unknown_site_selector(self, side_group_deck):
        """The selector vocabulary is CLOSED: an unknown selector is a
        site this consumer cannot classify."""
        def cut(a):
            self._channel(a)["site_selector"] = "allylic"
        self._shape_case(side_group_deck, cut, r"site_selector")

    def test_rejects_sites_per_unit_match_count_mismatch(
            self, side_group_deck):
        """sites_per_unit is CHECKED against the SERIALIZED selector
        resolution (site_atom_indices), never trusted -- the round-72 law,
        validated from data alone (round-73: no solver-backstop gap)."""
        def cut(a):
            self._channel(a)["sites_per_unit"] = 2.0
        self._shape_case(side_group_deck, cut, r"contradicts")

    def test_rejects_malformed_site_atom_indices(self, side_group_deck):
        def cut(a):
            self._channel(a)["site_atom_indices"] = [2, 2]
        self._shape_case(side_group_deck, cut, r"site_atom_indices")

    def test_rejects_same_atom_set_channels(self, side_group_deck):
        """Two rate channels resolving to ONE structural site double-carry
        the loss (round-72 P1) -- checkable from the serialized atom sets
        alone."""
        def cut(a):
            block = self._block(a)
            twin = _dc(block["channels"][0])
            twin["label"] = "other_C-Br"
            block["channels"].append(twin)
        self._shape_case(side_group_deck, cut, r"SAME.*atom set")

    def _overlap_case(self, deck, base_indices, twin_indices):
        """Hand-build a two-channel v2 block whose channels' atom sets
        OVERLAP without being identical (round-75 P1-1) -- each channel
        still satisfies len(site_atom_indices) == sites_per_unit -- and
        hand it straight to the runner's artifact guard (v2 spawns no
        feature pool, so nothing else needs planting)."""
        from rmgpy.tools.polymer_moments_runner import (
            _check_side_group_homolysis)
        artifact = _side_group_artifact(deck)
        block = self._block(artifact)
        base = block["channels"][0]
        base["site_atom_indices"] = sorted(base_indices)
        base["sites_per_unit"] = float(len(base_indices))
        twin = _dc(base)
        twin["label"] = "other_C-Br"
        twin["site_atom_indices"] = sorted(twin_indices)
        twin["sites_per_unit"] = float(len(twin_indices))
        block["channels"].append(twin)
        with pytest.raises(ValueError, match=r"overlap"):
            _check_side_group_homolysis(artifact)

    def test_rejects_superset_overlap_atom_set_channels(
            self, side_group_deck):
        """Round-75 P1-1: the same-set guard keyed on identical
        frozensets lets a twin claiming a strict SUPERSET of the base
        channel's Br site pass -- atom 2 is double-carried while both
        channels satisfy the count law. ANY non-empty pairwise
        intersection per gas element must reject."""
        self._overlap_case(side_group_deck, [2], [2, 4])

    def test_rejects_subset_overlap_atom_set_channels(
            self, side_group_deck):
        """Round-75 P1-1, subset direction: the twin claims a strict
        SUBSET of the base channel's atom set -- same double-carry, same
        rejection."""
        self._overlap_case(side_group_deck, [2, 4], [2])

    def test_rejects_out_of_range_site_atom_indices(self, side_group_deck):
        """Round-75 P1-2: the serialized indices are 0-based positions in
        the carrier's monomer_adj_list atom order -- an index past the
        serialized atom count is structurally meaningless and must
        reject, not pass on non-negativity alone."""
        def cut(a):
            self._channel(a)["site_atom_indices"] = [999]
        self._shape_case(side_group_deck, cut, r"out of range")

    def test_rejects_carrier_without_monomer_adj_list(
            self, side_group_deck):
        """Round-75 P1-2: without the carrier's monomer_adj_list the
        indices cannot be bounds-anchored at all -- a kernel-carrying
        pool REQUIRES the serialized structure text."""
        def cut(a):
            next(p for p in a["pools"]
                 if p["label"] == "PVBr")["monomer_adj_list"] = ""
        self._shape_case(side_group_deck, cut, r"monomer_adj_list")

    def test_rejects_carrier_block_plus_defect(self, side_group_deck):
        """A v2 carrier carries NO chain_mass_defect_g_mol: its pool
        moments are unchanged by the kernel and multi-loss X depletion is
        tracked on the carrier's auxiliary Z inventory, not a feature-pool
        defect. A carrier with both the block AND a defect field is
        malformed."""
        def cut(a):
            next(p for p in a["pools"] if p["label"] == "PVBr")[
                "chain_mass_defect_g_mol"] = 79.904
        self._shape_case(side_group_deck, cut,
                         r"chain_mass_defect_g_mol")

    def test_rejects_duplicate_channel_labels(self, side_group_deck):
        def cut(a):
            block = self._block(a)
            block["channels"].append(_dc(block["channels"][0]))
        self._shape_case(side_group_deck, cut, r"duplicate")

    def test_rejects_gas_mw_mismatch(self, side_group_deck):
        """gas_mw_g_mol must pin M_X of the gas_product -- it is the value
        the feature pool's defect must carry."""
        def cut(a):
            self._channel(a)["gas_mw_g_mol"] = 5.0
        self._shape_case(side_group_deck, cut, r"gas_mw_g_mol")

    def test_rejects_non_monoatomic_gas_product(self, side_group_deck):
        def cut(a):
            self._channel(a)["gas_product"] = "[CH3]"
        self._shape_case(side_group_deck, cut, r"monoatomic|mono-radical")

    def test_rejects_unresolvable_gas_species(self, side_group_deck):
        """The gas routing target must exist in chem.yaml -- the ejected X
        would silently vanish otherwise (un-conserved mass)."""
        def cut(a):
            self._channel(a)["gas_species"] = "Xe(99)"
        self._shape_case(side_group_deck, cut,
                         r"not in the deck's species list")

    def test_rejects_unconfigured_carrier(self, side_group_deck):
        """A block on an unconfigured carrier is a silently dropped
        kernel."""
        def cut(a):
            conv = a["conventions"]
            conv["configured_pools"] = [
                lbl for lbl in conv["configured_pools"] if lbl != "PVBr"]
        self._shape_case(side_group_deck, cut,
                         r"'PVBr'.*configured_pools")

    def test_rejects_malformed_defect_value(self, side_group_deck):
        """Any serialized chain_mass_defect_g_mol -- on any pool -- must be
        a finite value > 0 (shape validation, never truthiness). v2 spawns
        no feature pool, so the field is planted on a copy-carried-defect
        pool shape riding alongside the carrier."""
        def cut(a):
            extra = _dc(next(p for p in a["pools"]
                             if p["label"] == "PVBr"))
            extra["label"] = "PVBr_defect_mod"
            extra.pop("side_group_homolysis", None)
            extra["chain_mass_defect_g_mol"] = True
            a["pools"].append(extra)
        self._shape_case(side_group_deck, cut, r"finite value")

    def test_rejects_k_unzip_coexistence(self, side_group_deck):
        """Generation-side mutual exclusion re-enforced at the artifact
        boundary: side_group_homolysis + k_unzip > 0 double-carries
        degradation."""
        def cut(a):
            entry = next(p for p in a["pools"] if p["label"] == "PVBr")
            entry["channels"]["unzip"]["A"] = 2.0
            entry["monomer_routing"] = "Br(2)"
        self._shape_case(side_group_deck, cut, r"mutually exclusive")

    def test_kernel_supersession_census_and_refused_marker_survives(
            self, tmp_path):
        """FR1-K1 supersession, machine-checkable across the artifact
        boundary: with the SGH kernel-v2 live, refused conduit-deferred
        rows stay marked refused on the reconstructed reactions and the
        rebuilt oracle's side_group_supersession_census names them (so the
        refused rows they supersede are accounted for, not a dead system).
        v2 spawns no feature pool -- the census keys on the kernel-enabled
        carrier alone."""
        deck = _side_group_deck(tmp_path, with_refused=True)
        artifact = _side_group_artifact(deck)
        assert artifact["schema_version"] == "3.0"
        marked = [e for e in artifact["reactions"] if e.get("refused")]
        assert len(marked) == 1
        assert marked[0]["refused_reason"] == "conduit-deferred"

        rs, core, _ = _build_side_group(
            deck, artifact, initial_moles={"N2(1)": 1.0, "Br(2)": 0.5})
        # refused marker survived the boundary
        assert int(rs.reaction_refused.sum()) == 1
        # the SGH kernel is live ...
        assert rs.sgh_enabled[0] == 1
        # ... and its supersession census names the superseded refused rows
        census = list(rs.side_group_supersession_census)
        assert census and census[0]["pool"] == "PVBr"
        assert census[0]["superseded_rows"]


# ---------------------------------------------------------------------------
# End-radical depropagation block, consumer side (schema 2.8, r74 SS2 kernel;
# r78 serialization rulings): strict block guard (closed keys, pinned
# kernel/units/recipe/gate_width), version cross-check, daughter/parent
# closure (carrier is a provenance-pinned end-radical daughter of a
# homolysis carrier; sibling symmetry), the solver's validate_configuration
# exclusion mirror (k_unzip / radical_qssa_unzip / k_homolysis -- plus the
# r78-adjudicated k_scission rejection), the gas routing/MW cross-pins, and
# truthful 2.8 acceptance (the kernel triplet is wired into the rebuilt
# oracle's kdep_* flattened arrays).
# ---------------------------------------------------------------------------


def _deprop_deck(tmp_path):
    """A parent pool declaring k_homolysis + k_depropagation + monomer
    routing, plus its two eagerly-configured end-radical daughter pools
    (which carry the live kernel), written out as chem.yaml +
    polymer_pools.json exactly like an RMG run would."""
    n2 = _spc("N#N", "N2", index=1)
    c3h6 = _spc("C=CC", "C3H6", index=5)
    pool = Polymer(label="PP", monomer="[CH2][CH](C)",
                   end_groups=["[H]", "[H]"], cutoff=3,
                   moments=[1.0, 50.0, 3000.0], initial_mass=0.0,
                   k_homolysis={"A": 1.0e13, "n": 0.5, "Ea": 1.2e5},
                   k_depropagation={"A": 9.4e14, "n": 0.0, "Ea": 117152.0})
    pool.monomer_product_species = c3h6
    prim, sec = pool.generate_end_radical_daughters()
    mus = []
    for base in ("PP", prim.label, sec.label):
        mus += [_mu(f"{base}_mu{k}") for k in range(3)]
    core = [n2, c3h6] + mus
    data, index_map = generate_cantera_data(core, [],
                                            return_reaction_index_map=True)
    chem_path = os.path.join(str(tmp_path), "chem.yaml")
    with open(chem_path, "w") as fh:
        yaml.dump(data, fh, sort_keys=False, default_flow_style=None)
    artifact = build_polymer_moments_artifact(
        [pool, prim, sec], core_species=core, core_reactions=[],
        configured_pool_labels=["PP", prim.label, sec.label],
        condensed_species=mus, cantera_index_map=index_map)
    art_path = os.path.join(str(tmp_path), "polymer_pools.json")
    with open(art_path, "w") as fh:
        json.dump(artifact, fh, indent=2, default=str)
    return chem_path, art_path


@pytest.fixture
def deprop_deck(tmp_path):
    return _deprop_deck(tmp_path)


def _deprop_artifact(deck):
    with open(deck[1]) as fh:
        return json.load(fh)


def _build_deprop(deck, artifact=None, initial_moles=None):
    chem_path, art_path = deck
    if artifact is None:
        with open(art_path) as fh:
            artifact = json.load(fh)
    species, reactions = load_chem_yaml(chem_path)
    return build_system_from_artifact(
        artifact, species, reactions, T0=800.0, P=1.0e5, V_poly=1.0,
        initial_moles=initial_moles or {"N2(1)": 1.0},
        mass_transfer_spec=[])


class TestEndRadicalDepropagationConsumer:
    """Schema 2.8 end_radical_depropagation contract, consumer side."""

    def test_round_trip_kernel_flattened_on_daughters(self, deprop_deck):
        """RED pin (full round trip, truthful 2.8 acceptance): a real
        kernel-enabled deck serialized by the emitter loads GREEN, wires
        the triplet into BOTH daughter pools' kdep_* flattened arrays with
        the gas routing resolved -- and leaves the parent kernel-free."""
        artifact = _deprop_artifact(deprop_deck)
        assert artifact["schema_version"] == "2.8"
        for lbl in ("PP_rad_primary_end", "PP_rad_secondary_end"):
            entry = next(p for p in artifact["pools"] if p["label"] == lbl)
            block = entry.get("end_radical_depropagation")
            assert isinstance(block, dict)
            assert block["kernel"] == "end_radical_depropagation/1"
            assert block["recipe_revision"] == \
                "2026-07-06-end-radical-depropagation"
            assert block["gate_width"] == 1.0e-2

        rs, core, _ = _build_deprop(deprop_deck, artifact)
        assert len(rs.polymer_pools) == 3
        by_label = {p.label: i for i, p in enumerate(rs.polymer_pools)}
        i_gas = [s.label for s in core].index("C3H6(5)")
        assert rs.kdep_enabled[by_label["PP"]] == 0
        for lbl in ("PP_rad_primary_end", "PP_rad_secondary_end"):
            k = by_label[lbl]
            assert rs.kdep_enabled[k] == 1
            assert rs.kdep_A[k] == pytest.approx(9.4e14)
            assert rs.kdep_n[k] == pytest.approx(0.0)
            assert rs.kdep_Ea[k] == pytest.approx(117152.0)
            assert int(rs.kdep_gas[k]) == i_gas

    def test_accepts_2_8_without_deprop(self, deck):
        """2.8 is now an implemented minor: a 2.8 stamp with no
        end_radical_depropagation block anywhere loads (mirror of
        test_accepts_2_6_without_homolysis)."""
        chem_path, art_path = deck
        with open(art_path) as fh:
            artifact = json.load(fh)
        artifact["schema_version"] = "2.8"
        species, reactions = load_chem_yaml(chem_path)
        rs, _, _ = build_system_from_artifact(
            artifact, species, reactions, T0=800.0, P=1.0e5, V_poly=1.0,
            initial_moles={"N2(1)": 1.0}, mass_transfer_spec=[])
        assert len(rs.polymer_pools) == 1
        assert rs.kdep_enabled[0] == 0

    def test_rejects_block_in_2_7_stamped_artifact(self, deprop_deck):
        """Vocabulary/version cross-check (the 2.6/2.7 precedent): a
        below-2.8 artifact carrying an end_radical_depropagation block is
        malformed -- the emitter stamps 2.8 whenever it writes one."""
        artifact = _deprop_artifact(deprop_deck)
        artifact["schema_version"] = "2.7"
        with pytest.raises(ValueError,
                           match=r"end_radical_depropagation.*2\.8"):
            _build_deprop(deprop_deck, artifact)

    def _shape_case(self, deck, mutate, match):
        artifact = _deprop_artifact(deck)
        mutate(artifact)
        with pytest.raises(ValueError, match=match):
            _build_deprop(deck, artifact)

    @staticmethod
    def _entry(artifact, label="PP_rad_primary_end"):
        return next(p for p in artifact["pools"] if p["label"] == label)

    @classmethod
    def _block(cls, artifact, label="PP_rad_primary_end"):
        return cls._entry(artifact, label)["end_radical_depropagation"]

    def test_rejects_missing_kinetics_key(self, deprop_deck):
        def cut(a):
            del self._block(a)["kinetics"]["Ea"]
        self._shape_case(deprop_deck, cut, r"kinetics")

    def test_rejects_wrong_kinetics_units(self, deprop_deck):
        def cut(a):
            self._block(a)["kinetics"]["units"]["A"] = "1/min"
        self._shape_case(deprop_deck, cut, r"units")

    def test_rejects_non_positive_A(self, deprop_deck):
        def cut(a):
            self._block(a)["kinetics"]["A"] = 0.0
        self._shape_case(deprop_deck, cut, r"A")

    def test_rejects_boolean_kinetics_value(self, deprop_deck):
        def cut(a):
            self._block(a)["kinetics"]["A"] = True
        self._shape_case(deprop_deck, cut, r"A")

    def test_rejects_present_disabled_block(self, deprop_deck):
        def cut(a):
            self._block(a)["enabled"] = False
        self._shape_case(deprop_deck, cut, r"present-disabled|enabled")

    def test_rejects_unknown_block_key(self, deprop_deck):
        def cut(a):
            self._block(a)["turbo"] = True
        self._shape_case(deprop_deck, cut, r"unknown key")

    def test_rejects_unknown_kernel_name(self, deprop_deck):
        def cut(a):
            self._block(a)["kernel"] = "beta_scission_zip/9"
        self._shape_case(deprop_deck, cut, r"kernel")

    def test_rejects_unknown_recipe_revision(self, deprop_deck):
        def cut(a):
            self._block(a)["recipe_revision"] = "2026-01-01-bogus"
        self._shape_case(deprop_deck, cut, r"recipe_revision")

    def test_rejects_tampered_recipe(self, deprop_deck):
        def cut(a):
            block = self._block(a)
            key = sorted(block["recipe"])[0]
            block["recipe"][key] = "dmu0 = 0 forever (stall recipe)"
        self._shape_case(deprop_deck, cut, r"recipe")

    def test_rejects_wrong_gate_width(self, deprop_deck):
        """The gate_width field is pinned BITWISE to the solver constant
        (KDEP_GATE_WIDTH = 1e-2): the RHS the artifact claims to replicate
        integrates with exactly that width, so any other value is a
        different law (the 1e-12 TA-twin contract)."""
        def cut(a):
            self._block(a)["gate_width"] = 5.0e-2
        self._shape_case(deprop_deck, cut, r"gate_width")

    def test_rejects_missing_gate_width(self, deprop_deck):
        def cut(a):
            del self._block(a)["gate_width"]
        self._shape_case(deprop_deck, cut, r"gate_width|missing")

    def test_rejects_gas_species_routing_mismatch(self, deprop_deck):
        """block.gas_species and the pool surface's monomer_routing must
        name the SAME species (one routing, cross-pinned)."""
        def cut(a):
            self._block(a)["gas_species"] = "N2(1)"
        self._shape_case(deprop_deck, cut, r"gas_species|routing")

    def test_rejects_gas_mw_mismatch(self, deprop_deck):
        """r78 mass pin: gas_mw_g_mol must equal the carrier's
        monomer_mw_g_mol (each event moves exactly one repeat unit) --
        anything else mints/destroys mass on every unzip event."""
        def cut(a):
            self._block(a)["gas_mw_g_mol"] = 2.016
        self._shape_case(deprop_deck, cut, r"mw|mass")

    def test_rejects_missing_monomer_routing(self, deprop_deck):
        def cut(a):
            self._entry(a)["monomer_routing"] = None
        self._shape_case(deprop_deck, cut, r"routing")

    def test_rejects_sibling_block_dropped(self, deprop_deck):
        """The producer copies ONE parent triplet onto BOTH daughters; an
        artifact where only one sibling carries the block is corrupted."""
        def cut(a):
            del self._entry(a, "PP_rad_secondary_end")[
                "end_radical_depropagation"]
        self._shape_case(deprop_deck, cut, r"sibling")

    def test_rejects_sibling_kinetics_divergence(self, deprop_deck):
        def cut(a):
            self._block(a, "PP_rad_secondary_end")["kinetics"]["A"] = 1.0e3
        self._shape_case(deprop_deck, cut, r"sibling")

    def test_rejects_provenance_stripped_carrier(self, deprop_deck):
        """The kernel's home is a producer-spawned end-radical daughter; a
        carrier claiming user/input provenance is hand-edited."""
        def cut(a):
            self._entry(a)["spawn_event_metadata"] = {"source": "input"}
        self._shape_case(deprop_deck, cut, r"provenance|spawn")

    def test_rejects_carrier_without_parent_homolysis(self, deprop_deck):
        """Daughter/parent closure: the carrier's parent pool entry must
        carry the homolysis_initiation block that spawned it -- a deprop
        daughter with no initiation feed is an orphan shape production
        cannot generate."""
        def cut(a):
            del self._entry(a, "PP")["homolysis_initiation"]
        self._shape_case(deprop_deck, cut, r"homolysis")

    def test_rejects_unzip_on_carrier(self, deprop_deck):
        """Solver validate_configuration exclusion mirror #1: legacy
        k_unzip is the scalar form of the SAME chain-end release event."""
        def cut(a):
            self._entry(a)["channels"]["unzip"]["A"] = 0.5
        self._shape_case(deprop_deck, cut,
                         r"PP_rad_primary_end.*unzip")

    def test_rejects_scission_on_carrier(self, deprop_deck):
        """r78 adjudication (RED-pinned): k_scission + k_depropagation on
        one pool is a direct-config-only shape production cannot generate
        (daughters are born with k_scission = 0), so no generating run
        ever integrated the combined law -- reject, never adapt."""
        def cut(a):
            self._entry(a)["channels"]["scission"]["A"] = 0.5
        self._shape_case(deprop_deck, cut,
                         r"PP_rad_primary_end.*scission")

    def test_rejects_homolysis_block_on_carrier(self, deprop_deck):
        """Solver validate_configuration exclusion mirror #3:
        multi-generation homolysis of radical-ended chains is DEFERRED
        (r74 SS3) -- a carrier entry claiming both kernels is malformed."""
        def cut(a):
            entry = self._entry(a)
            entry["homolysis_initiation"] = _dc(
                self._entry(a, "PP")["homolysis_initiation"])
        self._shape_case(deprop_deck, cut, r"homolysis")

    def test_rejects_unconfigured_carrier(self, deprop_deck):
        def cut(a):
            conv = a["conventions"]
            conv["configured_pools"] = [
                lbl for lbl in conv["configured_pools"]
                if lbl != "PP_rad_primary_end"]
        self._shape_case(deprop_deck, cut,
                         r"PP_rad_primary_end.*(configured|classif)")

    def test_rejects_defect_bearing_carrier(self, deprop_deck):
        """r79 P1 (RED-pinned): chain_mass_defect_g_mol on a deprop
        carrier mints mass at terminal DP1 events (condensed drains
        R*(MW - defect) while gas credits R*MW -> +R*defect per event);
        the side-group orphan guard deliberately legalizes copied defect
        pools (non-side-group provenance), so without THIS rejection the
        shape loads. v2 defect-chain depropagation would need a different
        mass law / gas product."""
        def cut(a):
            self._entry(a)["chain_mass_defect_g_mol"] = 79.904
        self._shape_case(deprop_deck, cut,
                         r"PP_rad_primary_end.*chain_mass_defect")

    def test_rejects_spawned_classified_carrier(self, deprop_deck):
        """r79 item 2: the spawned-classified conjunct, named even when
        the carrier is ALSO missing from configured_pools (the 2.6
        spawned-first ordering rationale). Isolated from the sibling 2.6
        daughter guard by stripping the parent's homolysis block (which
        fires first on shared shapes)."""
        def cut(a):
            del self._entry(a, "PP")["homolysis_initiation"]
            conv = a["conventions"]
            conv["configured_pools"] = [
                lbl for lbl in conv["configured_pools"]
                if lbl != "PP_rad_primary_end"]
            conv["spawned_pools"] = ["PP_rad_primary_end"]
        self._shape_case(deprop_deck, cut,
                         r"PP_rad_primary_end.*spawned_pools")

    def test_rejects_uncondensed_carrier(self, deprop_deck):
        """r79 item 2: an un-condensed carrier's moment slots would be
        phase-classified GAS (the item-16 mass-balance hazard) while the
        kernel drains them as condensed mass. Isolated from the sibling
        2.6 daughter guard by stripping the parent's homolysis block."""
        def cut(a):
            del self._entry(a, "PP")["homolysis_initiation"]
            a["conventions"]["condensed_species"] = [
                lbl for lbl in a["conventions"]["condensed_species"]
                if not lbl.startswith("PP_rad_primary_end")]
        self._shape_case(deprop_deck, cut,
                         r"PP_rad_primary_end.*condensed")

    def test_rejects_qssa_channel_on_carrier(self, deprop_deck):
        """r79 item 2 / solver validate_configuration exclusion mirror
        #2: the QSSA channel's depropagation block IS this lumped
        chain-end channel -- a carrier entry claiming both double-counts
        the unzip flux."""
        def cut(a):
            self._entry(a)["channels"]["radical_qssa_unzip"] = {
                "enabled": True}
        self._shape_case(deprop_deck, cut,
                         r"PP_rad_primary_end.*radical_qssa_unzip")


# ---------------------------------------------------------------------------
# Round-31 P2: --atol/--rtol replay parity. atol is a MODEL knob for the
# polymer kernel -- it anchors the r81 accepted-state floors
# max(SMALL_EPS, 100*atol) and with them the exhaustion band (E) and the
# cone-margin band (M) of the near-exhaustion bundle limiter (format doc
# section 4b). A replay must be able to match the generating deck's
# tolerances, and BOTH consumers (the solver-backed runner and the numpy
# oracle consumer) must anchor the SAME floors from the same atol.
# ---------------------------------------------------------------------------
_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
from numpy_moments_consumer import ArtifactConsumer  # noqa: E402


class TestAtolReplayParity:

    def _load(self, deck):
        chem_path, art_path = deck
        with open(art_path) as fh:
            artifact = json.load(fh)
        species, reactions = load_chem_yaml(chem_path)
        return chem_path, art_path, artifact, species, reactions

    def test_atol_reaches_solver_floors_and_numpy_consumer(self, deck):
        """build_system_from_artifact(atol=...) must anchor the solver's
        per-pool r81 floors at max(SMALL_EPS, 100*atol) AND the DASPK
        atol array; ArtifactConsumer(atol=...) must anchor the identical
        mu_floor -- the two consumers share one regularization envelope."""
        _, _, artifact, species, reactions = self._load(deck)
        rs, core, _ = build_system_from_artifact(
            artifact, species, reactions, T0=650.0, P=1.0e5, V_poly=1.0,
            initial_moles={"N2(1)": 1.0}, mass_transfer_spec=[],
            atol=1e-12, rtol=1e-4)
        floors = np.asarray(rs._pool_mu_floors)
        assert floors.shape == (1, 3)
        assert np.all(floors == 1.0e-10)          # max(SMALL_EPS, 100*atol)
        atol_arr = np.asarray(rs.atol_array)
        assert np.all(atol_arr == 1.0e-12)
        assert np.all(np.asarray(rs.rtol_array) == 1.0e-4)
        # the numpy oracle consumer anchors the SAME floor from the same
        # deck atol (scalar-atol form).
        consumer = ArtifactConsumer(
            artifact, [s.label for s in core], P=1.0e5, V_poly=1.0,
            atol=1e-12)
        assert consumer.mu_floor == floors[0, 0] == 1.0e-10

    def test_default_tolerances_unchanged(self, deck):
        """Backward compatibility: omitting atol/rtol keeps the historical
        runner defaults (1e-16/1e-8 -> floors 1e-14), byte-identical
        replays for existing invocations."""
        _, _, artifact, species, reactions = self._load(deck)
        rs, core, _ = build_system_from_artifact(
            artifact, species, reactions, T0=650.0, P=1.0e5, V_poly=1.0,
            initial_moles={"N2(1)": 1.0}, mass_transfer_spec=[])
        assert np.all(np.asarray(rs._pool_mu_floors) == 1.0e-14)
        assert np.all(np.asarray(rs.atol_array) == 1.0e-16)
        assert np.all(np.asarray(rs.rtol_array) == 1.0e-8)
        consumer = ArtifactConsumer(
            artifact, [s.label for s in core], P=1.0e5, V_poly=1.0)
        assert consumer.mu_floor == 1.0e-14

    def test_cli_atol_propagates_through_main(self, deck, tmp_path,
                                              monkeypatch):
        """--atol/--rtol reach build_system_from_artifact from the CLI
        (captured via a pass-through wrapper; the run itself completes,
        so the wrapper exercised the real replay path)."""
        import rmgpy.tools.polymer_moments_runner as runner_mod
        chem_path, art_path, artifact, species, reactions = self._load(deck)
        captured = {}
        real_build = runner_mod.build_system_from_artifact

        def capturing_build(*args, **kwargs):
            captured["atol"] = kwargs.get("atol")
            captured["rtol"] = kwargs.get("rtol")
            return real_build(*args, **kwargs)

        monkeypatch.setattr(runner_mod, "build_system_from_artifact",
                            capturing_build)
        profile_path = os.path.join(str(tmp_path), "profile.json")
        with open(profile_path, "w") as fh:
            json.dump([{"t_end": 1.0e-4, "T": 650.0}], fh)
        moles_path = os.path.join(str(tmp_path), "moles.json")
        with open(moles_path, "w") as fh:
            json.dump({"N2(1)": 1.0}, fh)
        out_path = os.path.join(str(tmp_path), "out.csv")
        runner_mod.main([
            "--artifact", art_path, "--chem", chem_path,
            "--t-profile", profile_path, "--n-points", "2",
            "--v-poly", "1.0", "--initial-moles", moles_path,
            "--output", out_path,
            "--atol", "1e-12", "--rtol", "1e-4"])
        assert captured == {"atol": 1e-12, "rtol": 1e-4}
        assert os.path.exists(out_path)


# ---------------------------------------------------------------------------
# M18.1 item 2 (rounds 30-31): the poly_102 regen-#3 79-species/82-reaction
# DEATH-CONFIGURATION replay, reconstructed from the run's own artifacts.
# Requires the READ-ONLY forensic run dir; skipped wherever it is absent.
# ---------------------------------------------------------------------------
_POLY102_RUN = "/home/alon/runs/RMG/poly_102_conduit3"


@pytest.mark.functional
@pytest.mark.skipif(not os.path.isdir(_POLY102_RUN),
                    reason="poly_102_conduit3 forensic run dir not present")
class TestRegen3SavedCoreReplay:
    """Replay the ACTUAL regen-#3 death configuration under the fixed law
    (round-29 N2 soft-min + round-30 N1 independent cone gate).

    Reconstruction sources (and their fidelity):
    * cantera/chem0079.yaml -- the POST-rebuild 79-species/82-reaction
      core the crashing stage-1 simulate integrated (FAITHFUL: species,
      kinetics, ordering match the crash dump's core list).
    * chemkin/polymer_pools.json -- iteration-6 sidecar, stamped
      stale_topology=true (PRE-rebuild emission): pool configs/moments
      are deck-faithful, but per-row refusal/liveness state may differ
      from what the post-rebuild solver ran (6 of 21 pool-coupled rows
      refused pre-rebuild). Loaded via allow_stale=True -- DIAGNOSTIC
      reconstruction, not a certified artifact replay.
    * RMG.log crash dump (t = 22.936 s IDID=-7 failure) -- the EXACT
      accepted-state y vector (79 entries), core label order, and total
      volume, parsed from the last "Core species names/moles" dump.
    * input.py deck -- T = 1100 K, P = 1 bar, initialMoles
      {N2: 0.9, H: 0.001}, atol = 1e-12, rtol = 1e-4 (via the round-31
      --atol/--rtol replay-parity plumbing); V_poly = 1.182975e-4 m^3
      (validated below: the reconstructed system's total volume at the
      crash state matches the dump's 0.08252302901147644 m^3 to ~3e-10
      relative).

    Known fidelity gaps (documented, not fabricated): (a) the stale
    sidecar's row liveness may not equal the post-rebuild solver's -- the
    retained-row set is the pre-rebuild one; (b) the crashing simulate's
    DASPK internal state (step history, order, Jacobian age) is not
    recoverable -- we re-initialize AT the accepted crash state; (c) the
    0 -> 22.936 s trajectory of the original run is not replayed here
    (the from-deck replay reproduces its shape but grinds for minutes in
    the t ~ 10-20 s multi-daughter decay regime; done out-of-band in the
    round-31 forensics, not in-test).

    Pins: at the EXACT logged crash state the RHS is finite, the
    out-of-cone daughter (mu = [9.956e-4, 1.497e-7, 1.819e-3],
    b1 = mu2/mu1 ~ 1.2e4) has its cone-deepening drain OFF (nonnegative
    net moment rates; FD-Jacobian entries w.r.t. its mu1 finite and
    BOUNDED -- the old law fed DASPK d(b1)/d(mu1) ~ 8e10 here), and
    DASPK integrates THROUGH the death point across the full 20-30 s
    crash window without IDID=-7 (the 100 s horizon was verified
    out-of-band: completes at ~400 s wall, daughter parks at
    mu1 ~ 4.5e-6 with the gate holding the drain off)."""

    _CRASH_T = 22.936          # s, from the IDID=-7 resurrection error
    _DUMP_V = 0.08252302901147644   # m^3, "Error: Volume:" crash dump

    def _reconstruct(self):
        with open(os.path.join(_POLY102_RUN,
                               "chemkin/polymer_pools.json")) as fh:
            artifact = json.load(fh)
        species, reactions = load_chem_yaml(
            os.path.join(_POLY102_RUN, "cantera/chem0079.yaml"))
        rs, core, all_rxns = build_system_from_artifact(
            artifact, species, reactions,
            T0=1100.0, P=1.0e5, V_poly=1.182975e-4,
            initial_moles={"N2": 0.90, "H(1)": 0.001},
            mass_transfer_spec=[], initial_moments=None,
            allow_stale=True, atol=1e-12, rtol=1e-4)
        return rs, core

    def _crash_state(self, core):
        """Parse the LAST crash dump (the 79-core one) from RMG.log and
        map it positionally onto the reconstructed core (the dump prints
        formula-based labels for a few species and compressed chemkin
        ids for the mu slots; the ORDER is the core order -- anchored
        below on the shared trailing RMG indices)."""
        import re as _re
        with open(os.path.join(_POLY102_RUN, "RMG.log")) as fh:
            log = fh.read()
        labels = [s.strip().strip("'") for s in _re.findall(
            r"Error: Core species names: \[(.*?)\]", log, _re.S)[-1].split(",")]
        moles = eval(_re.findall(
            r"Error: Core species moles: array\((\[.*?\])\)",
            log, _re.S)[-1].replace("\n", " "))
        core_labels = [s.label for s in core]
        assert len(core_labels) == len(labels) == 79
        for i, (cl, dl) in enumerate(zip(core_labels, labels)):
            mc = _re.search(r"\((\d+)\)$", cl)
            md = _re.search(r"\((\d+)\)$", dl)
            if mc and md:
                assert mc.group(1) == md.group(1), (i, cl, dl)
        return np.array(moles, dtype=float)

    def _crash_ready(self, rtol):
        """Build the death configuration and initialize AT the exact
        logged crash state (deck atol 1e-12; rtol per test)."""
        with open(os.path.join(_POLY102_RUN,
                               "chemkin/polymer_pools.json")) as fh:
            artifact = json.load(fh)
        species, reactions = load_chem_yaml(
            os.path.join(_POLY102_RUN, "cantera/chem0079.yaml"))
        rs, core, _ = build_system_from_artifact(
            artifact, species, reactions,
            T0=1100.0, P=1.0e5, V_poly=1.182975e-4,
            initial_moles={"N2": 0.90, "H(1)": 0.001},
            mass_transfer_spec=[], initial_moments=None,
            allow_stale=True, atol=1e-12, rtol=rtol)
        y = np.zeros(len(rs.y))
        y[:79] = self._crash_state(core)
        dn = np.asarray(rs.residual(self._CRASH_T, y.copy(),
                                    np.zeros_like(y))[0])
        rs.initialize(self._CRASH_T, y.copy(), dn.copy(),
                      atol=1e-12, rtol=rtol)
        return rs

    def test_crash_window_cost_bar_deck_tolerances(self):
        """CHEAP cost-bar regression at the deck tolerances (rtol=1e-4,
        atol=1e-12) under the no-F-gate law (round-35 revert).

        LINEAGE (rounds 33-35): this pin descends from the strict-xfail
        cost-bar test of rounds 33-34. The round-33 linear sub-floor
        ramp re-tripped IDID=-7 at t = 22.996 here; the round-34 wide
        log-floor ramp survived but took 1240.7 s wall to reach t = 23.0
        (~388x the round-32 law). Round-35 CONVICTED the floor-anchored
        ramp mechanism class (parked pools sit 2-4 floors inside any
        such ramp = permanent Jacobian resident; du/dmu = u'(F)/f peaks
        near the floor) and reverted the gate; the xfail flips to this
        plain green pin. The round-32 law does 22.936 -> 23.0 in ~3.4 s;
        <60 s catches any grind-class regression without flakiness."""
        rs = self._crash_ready(rtol=1e-4)
        wall = time.monotonic()
        rs.advance(23.0)
        assert time.monotonic() - wall < 60.0, "crash-window cost bar"
        yv = np.asarray(rs.y)
        assert np.all(np.isfinite(yv))
        rs._assert_pool_moments_accepted()

    def test_prestress_crash_window_rtol_1e6(self):
        """Round-35 pre-regen stress gate 5a: the exact-crash replay at
        rtol=1e-6, atol=1e-12 crosses t = 23.0 (and 24.0) clean and
        fast under the no-F-gate law (measured ~2.7 s to t = 24.0).
        (Round-37 policy note: 1e-6 fixes the minimal K2 fixture but is
        NOT full-system safe -- see the canary below -- so NO regen deck
        tolerance is currently certified.)

        DOCUMENTED CANARY (round-35 finding, P1 input for round-36): at
        rtol=1e-6 this replay then dies IDID=-7 at t = 24.639 -- where
        the SAME law at rtol=1e-4 integrates to t = 100 (389 s wall).
        Tightening rtol inverts the failure mode on the 79-dim system.
        The canary below asserts the death so the finding stays pinned;
        if a future round fixes it, the canary raises loudly here and
        this test must be updated to assert full-window health instead."""
        rs = self._crash_ready(rtol=1e-6)
        wall = time.monotonic()
        rs.advance(23.0)
        rs.advance(24.0)
        assert time.monotonic() - wall < 60.0
        yv = np.asarray(rs.y)
        assert np.all(np.isfinite(yv))
        rs._assert_pool_moments_accepted()
        # the documented canary: the 1e-6 death at t ~ 24.639
        from pydas.daspk import DASPKError
        with pytest.raises(DASPKError):
            rs.advance(25.0)
        assert 24.0 < rs.t < 25.0     # died inside the window, as logged

    @pytest.mark.xfail(
        strict=True,
        reason="round-35 pre-regen stress gate 5b, RED (round-36 "
               "finding): the from-deck 79/82 window at rtol=1e-6 "
               "(atol=1e-12; NO regen tolerance is certified -- "
               "round-37 policy) GRINDS worse "
               "than at the convicted rtol=1e-4 -- 0 -> 13 s took "
               "202.3 s wall (1e-4: 16.5 s) and 13 -> 14 took 1765.8 s "
               "(29.4 min per sim-second; 1e-4 pre-gate law: ~16 min "
               "then IDID=-7 at t = 14.2445). Worse, mod_5 is STILL "
               "dragged sub-floor at 1e-6 (accepted 1.82e-11 mol at "
               "t = 14 vs the 1e-10 floor -- the H1 signature the "
               "tolerance conviction was expected to remove). "
               "Tightening rtol does NOT resolve the multi-daughter "
               "near-floor regime on the full system; combined with the "
               "crash-state replay's 1e-6 death at t = 24.639 (see the "
               "prestress canary), the tolerance interaction is a "
               "round-36 P1. This xfail asserts a 120 s budget at the "
               "t = 13 checkpoint so it burns ~3.5 min, not the full "
               "grind; a law/tolerance combination that traverses the "
               "window flips it loudly.")
    def test_prestress_fromdeck_window_rtol_1e6(self):
        with open(os.path.join(_POLY102_RUN,
                               "chemkin/polymer_pools.json")) as fh:
            artifact = json.load(fh)
        species, reactions = load_chem_yaml(
            os.path.join(_POLY102_RUN, "cantera/chem0079.yaml"))
        rs, core, _ = build_system_from_artifact(
            artifact, species, reactions,
            T0=1100.0, P=1.0e5, V_poly=1.0 / 1050.0,
            initial_moles={"N2": 0.90, "H(1)": 0.001},
            mass_transfer_spec=[], initial_moments=None,
            allow_stale=True, atol=1e-12, rtol=1e-6)
        y0 = np.array(rs.y, dtype=float).copy()
        dn0 = rs.residual(0.0, y0, np.zeros_like(y0))[0]
        rs.initialize(0.0, y0.copy(), dn0, atol=1e-12, rtol=1e-6)
        wall = time.monotonic()
        budgets = {13.0: 120.0, 14.0: 400.0, 14.5: 500.0}
        for t, budget in budgets.items():
            rs.advance(t)
            assert time.monotonic() - wall < budget, (t, "wall budget")
            y = np.asarray(rs.y)
            assert np.all(np.isfinite(y)), t
            rs._assert_pool_moments_accepted()
        assert rs.t >= 14.5 - 1e-9    # crossed the old death locus

    def test_death_configuration_replay_through_crash_window(self):
        rs, core = self._reconstruct()
        # reconstruction shape: the death configuration
        assert len(core) == 79
        assert len(rs.kf) == 82
        assert len(rs.polymer_pools) == 6
        assert np.all(np.asarray(rs._pool_mu_floors) == 1.0e-10)
        pools = {p.label: p.mu_indices for p in rs.polymer_pools}
        mod = pools["phenol_formaldehyde_mod"]

        y = np.zeros(len(rs.y))
        y[:79] = self._crash_state(core)
        # the parked out-of-cone daughter, exactly as logged
        assert y[mod[0]] == pytest.approx(9.95590522e-04, rel=1e-8, abs=0.0)
        assert y[mod[1]] == pytest.approx(1.49748726e-07, rel=1e-8, abs=0.0)
        assert y[mod[2]] == pytest.approx(1.81855523e-03, rel=1e-8, abs=0.0)
        assert y[mod[1]] < y[mod[0]]            # out of cone

        # RHS finite at the exact crash state; volume fidelity
        dn = np.asarray(rs.residual(self._CRASH_T, y.copy(),
                                    np.zeros_like(y))[0])
        assert np.all(np.isfinite(dn))
        assert rs.V == pytest.approx(self._DUMP_V, rel=1e-8, abs=0.0)
        # the cone gate holds the killer drain OFF: the out-of-cone
        # daughter's moments are not drained deeper into the cone
        for k in range(3):
            assert dn[mod[k]] >= 0.0, (k, dn[mod[k]])
        # FD-Jacobian entries w.r.t. the daughter's mu1: finite and
        # BOUNDED (the pre-fix law carried d(b1)/d(mu1) ~ 8e10 here)
        for h in (1.0e-10, 1.0e-8):
            yp = y.copy()
            yp[mod[1]] += h
            dnp = np.asarray(rs.residual(self._CRASH_T, yp,
                                         np.zeros_like(y))[0])
            jac = (dnp - dn) / h
            assert np.all(np.isfinite(jac))
            assert np.max(np.abs(jac)) < 1.0e6

        # DASPK integrates THROUGH the death point across the crash
        # window (the run died at t = 22.936 stepping toward 100).
        rs.initialize(self._CRASH_T, y.copy(), dn.copy(),
                      atol=1e-12, rtol=1e-4)
        wall = time.monotonic()
        for t in (23.0, 23.5, 24.0, 25.0, 26.0, 28.0, 30.0):
            rs.advance(t)
            yv = np.asarray(rs.y)
            assert np.all(np.isfinite(yv)), t
            rs._assert_pool_moments_accepted()
            # accepted daughter mu1 never dragged negative
            assert yv[mod[1]] >= 0.0, (t, yv[mod[1]])
        assert rs.t >= 30.0 - 1e-9
        # generous wall bound: observed ~170 s; a stepsize collapse
        # would blow this by orders of magnitude
        assert time.monotonic() - wall < 900.0
