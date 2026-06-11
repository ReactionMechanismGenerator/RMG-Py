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

    @pytest.mark.xfail(strict=True,
                       reason="item 15 residual (iii): tripwire derives MW from "
                              "molecule[]; consumer-world species crash "
                              "_sackur_tetrode_s_trans with math domain error")
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
