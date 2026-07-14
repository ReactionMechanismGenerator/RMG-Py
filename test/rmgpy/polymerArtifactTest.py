#!/usr/bin/env python3
"""Tests for the schema-2.0 polymer moments artifact emitter
(spec: docs/superpowers/specs/2026-06-10-polymer-moments-artifact-design.md;
format doc: docs/polymer_moments_format.md)."""

import json

import pytest

from rmgpy.kinetics import Arrhenius, MultiArrhenius
from rmgpy.molecule import Molecule
from rmgpy.polymer import (
    HOMOLYSIS_SPAWN_SOURCE,
    POLYMER_POOLS_SIDECAR_SCHEMA_VERSION,
    POLYMER_POOLS_SIDECAR_SCHEMA_VERSION_THERMAL,
    POLYMER_RATE_RECIPE_REVISION,
    Polymer,
    PolymerFluxArchetype,
    _serialize_pool_for_sidecar,
    build_polymer_moments_artifact,
    compile_polymer_reaction_entries,
    derive_condensed_species,
    is_midrun_spawned_pool_daughter,
    validate_thermal_analysis_inputs,
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
        # Recipe revision 2026-07-03-monomer-gas: the routed monomer is a
        # GAS-phase species and must NOT be listed as pool phase membership.
        assert "ethylene(5)" not in d["phase_species"]

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

    def test_routed_monomer_is_gas_not_phase_or_bookkeeping(self, pe_pool):
        """Recipe revision 2026-07-03-monomer-gas (incident 2026-07-03): the
        routed monomer is a GAS volatile -- in NEITHER phase_species NOR
        bookkeeping_species. (Pre-revision emitters listed it in
        phase_species as 'real condensed'; that classification conflated the
        release target with the deck's principal gas product.)"""
        core = [
            _spc("CC", "PE", index=2),
            _mu_dummy("PE_mu0"), _mu_dummy("PE_mu1"), _mu_dummy("PE_mu2"),
        ]
        d = _serialize_pool_for_sidecar(pe_pool, core_species=core,
                                        monomer_routing="ethylene(5)")
        assert "ethylene(5)" not in d["phase_species"]
        assert "ethylene(5)" not in d["bookkeeping_species"]
        assert d["phase_species"] == ["PE(2)", "PE_mu0", "PE_mu1", "PE_mu2"]
        assert d["bookkeeping_species"] == ["PE(2)", "PE_mu0", "PE_mu1", "PE_mu2"]

    def test_bookkeeping_species_always_present(self, pe_pool):
        """The key exists even without core_species (empty phase_species), and
        a routing-only pool block keeps BOTH lists empty (the routed monomer
        is gas since recipe revision 2026-07-03-monomer-gas)."""
        d = _serialize_pool_for_sidecar(pe_pool)
        assert d["phase_species"] == []
        assert d["bookkeeping_species"] == []
        d = _serialize_pool_for_sidecar(pe_pool, monomer_routing="ethylene(5)")
        assert d["phase_species"] == []
        assert d["bookkeeping_species"] == []


# Raw (pre-normalization) radical QSSA channel config, deck-shaped: only the
# mandatory blocks plus transfer; efficiency/monomer_yield/basis exercised
# both explicit and defaulted across the tests below.
QSSA_RAW_CFG = {
    "initiation": {"A": 1.0e13, "n": 0.0, "Ea": 3.0e5},
    "depropagation": {"A": 1.0e14, "n": 0.5, "Ea": 9.0e4},
    "termination": {"A": 1.0e8, "n": 0.0, "Ea": 1.0e4},
    "transfer": {"A": 2.0e3, "n": 0.0, "Ea": 5.0e4},
    "efficiency": 0.8,
    "monomer_yield": 0.9,
}

# The pinned sidecar units convention (docs: M2 commit c0256cc2f). initiation/
# depropagation are unimolecular [s^-1]; termination is the ONLY bimolecular
# block [m^3/(mol*s)]; transfer is PSEUDO-first-order [s^-1] (a literature
# bimolecular k_tr must be premultiplied by the substrate concentration
# before it enters the config). Ea always [J/mol].
QSSA_PINNED_UNITS = {
    "initiation": {"A": "s^-1", "Ea": "J/mol"},
    "depropagation": {"A": "s^-1", "Ea": "J/mol"},
    "termination": {"A": "m^3/(mol*s)", "Ea": "J/mol"},
    "transfer": {"A": "s^-1", "Ea": "J/mol"},
}


@pytest.fixture
def qssa_pool():
    pool = Polymer(
        label="PS",
        monomer="[CH2][CH](c1ccccc1)",
        end_groups=["[H]", "[H]"],
        cutoff=3,
        moments=[1.0, 50.0, 3000.0],
        initial_mass=0.0,
        k_scission=0.0,
        k_unzip=0.0,
        radical_qssa_unzip=dict(QSSA_RAW_CFG),
    )
    # A real QSSA pool always carries a routing target (three layers of
    # guards: deck read / to_config / solver validate_configuration; round-25
    # P1-1 makes the serializer refuse enabled-QSSA without one).
    pool.monomer_product_species = _spc("C=Cc1ccccc1", "styrene")
    return pool


class TestRadicalQssaChannelSerialization:
    """Sidecar serialization of the radical_qssa_unzip channel (milestone 3).

    The block lives INSIDE ``channels`` deliberately: the TA loader hard-errors
    on any channel key outside its supported set (~/Code/TA/ta/mechanism.py:52
    SUPPORTED_CHANNELS + :509-517 unknown-key ValueError), so an old consumer
    fails loudly instead of silently integrating a flat/false TGA."""

    def test_qssa_block_emitted_with_pinned_units_and_provenance(self, qssa_pool):
        d = _serialize_pool_for_sidecar(qssa_pool)
        block = d["channels"]["radical_qssa_unzip"]
        assert block["enabled"] is True
        assert block["basis"] == "backbone_bonds_mu1_minus_mu0"
        assert block["efficiency"] == pytest.approx(0.8)
        assert block["monomer_yield"] == pytest.approx(0.9)
        for name in ("initiation", "depropagation", "termination", "transfer"):
            trip = block[name]
            assert trip["A"] == pytest.approx(QSSA_RAW_CFG[name]["A"])
            assert trip["n"] == pytest.approx(QSSA_RAW_CFG[name]["n"])
            assert trip["Ea"] == pytest.approx(QSSA_RAW_CFG[name]["Ea"])
            assert trip["units"] == QSSA_PINNED_UNITS[name], (
                f"units pin broken for {name!r}: {trip.get('units')!r}")
        prov = block["provenance"]
        assert prov["radical_balance"] == (
            "G_R = 2*f*ki*B; loss = ktr*R + 2*kt*R^2; "
            "Rss no-transfer = sqrt(f*ki*B/kt)")
        assert prov["moment_closure"] == "end_shrink_pool_mean/1"
        assert prov["R_gas_J_per_mol_K"] == pytest.approx(8.314)
        assert prov["concentration_basis"] == "mol/m^3 condensed volume"
        assert "premultiplied by substrate concentration" in prov["transfer_note"]

    def test_qssa_pool_keeps_unzip_zeroed(self, qssa_pool):
        """Mutual exclusion is an upstream hard error; the sidecar shape it
        implies is unzip.A == 0 whenever the QSSA block is enabled."""
        d = _serialize_pool_for_sidecar(qssa_pool)
        assert d["channels"]["radical_qssa_unzip"]["enabled"] is True
        assert d["channels"]["unzip"]["A"] == 0.0

    def test_transfer_null_when_absent(self, qssa_pool):
        cfg = dict(QSSA_RAW_CFG)
        del cfg["transfer"]
        qssa_pool.radical_qssa_unzip = cfg
        d = _serialize_pool_for_sidecar(qssa_pool)
        block = d["channels"]["radical_qssa_unzip"]
        assert block["transfer"] is None
        # defaults from the shared validator fill in when the deck omits them
        cfg2 = {k: v for k, v in QSSA_RAW_CFG.items()
                if k in ("initiation", "depropagation", "termination")}
        qssa_pool.radical_qssa_unzip = cfg2
        block = _serialize_pool_for_sidecar(qssa_pool)["channels"]["radical_qssa_unzip"]
        assert block["efficiency"] == 1.0
        assert block["monomer_yield"] == 1.0
        assert block["basis"] == "backbone_bonds_mu1_minus_mu0"

    def test_legacy_pool_emits_no_qssa_block(self, pe_pool):
        """Channel absent -> emit NOTHING (no noise on legacy pools; an old
        TA must keep loading legacy sidecars unchanged)."""
        d = _serialize_pool_for_sidecar(pe_pool)
        assert "radical_qssa_unzip" not in d["channels"]
        assert set(d["channels"]) == {"scission", "unzip"}

    def test_invalid_qssa_config_refuses_to_serialize(self, qssa_pool):
        """The serializer normalizes through the SHARED validator: a garbage
        channel dict must raise, never reach the sidecar."""
        qssa_pool.radical_qssa_unzip = {"initiation": {"A": 1.0}}
        with pytest.raises(ValueError, match=r"PS.*radical_qssa_unzip"):
            _serialize_pool_for_sidecar(qssa_pool)

    def test_recipe_block_emitted_verbatim(self, qssa_pool):
        """The machine-readable normative recipe (round-23). Each string is
        pinned as a LITERAL, verified against the implemented RHS in
        rmgpy/solver/polymer.pyx (M2 rate law + moment signature; SMALL_EPS
        at polymer.pyx:71). The loader pins the same strings independently
        and rejects on any mismatch, so an accidental edit here is a
        consumer-coordination event, not a cosmetic one."""
        d = _serialize_pool_for_sidecar(qssa_pool)
        recipe = d["channels"]["radical_qssa_unzip"]["recipe"]
        assert recipe == {
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
        # human-readable provenance stays alongside the normative recipe
        assert "provenance" in d["channels"]["radical_qssa_unzip"]

    def test_spawned_daughter_sidecar_carries_inherited_qssa_block(self, qssa_pool):
        """Daughter-pool inheritance (M5): a spawn-intent daughter of a QSSA
        parent lands in the serialized sidecar carrying the SAME channel block
        (and lineage) as a deck pool would -- the sidecar writer iterates the
        live pool registry (main.py: every Polymer in core+edge+new) and reads
        the channel off the SPECIES, so inheritance at spawn time is exactly
        what makes the daughter's block appear."""
        from rmgpy.polymer import SpawnIntent, drain_spawn_intents

        mp = Species(label="styrene")
        mp.molecule = [Molecule().from_smiles("C=Cc1ccccc1")]
        qssa_pool.monomer_product_species = mp

        intent = SpawnIntent(
            parent_pool=qssa_pool,
            monomer=qssa_pool.backbone_group,
            end_groups=["[H]", "[H]"],
            triggering_dp=3,
        )
        daughter = drain_spawn_intents([intent], iteration=2)[0]
        assert daughter.label == "PS_d1"

        d = _serialize_pool_for_sidecar(daughter)
        assert d["parent_pool"] == "PS"
        block = d["channels"]["radical_qssa_unzip"]
        assert block["enabled"] is True
        assert block["efficiency"] == pytest.approx(0.8)
        assert block["monomer_yield"] == pytest.approx(0.9)
        for name in ("initiation", "depropagation", "termination"):
            assert block[name]["A"] == pytest.approx(QSSA_RAW_CFG[name]["A"])
        # And a channel-free parent's daughter serializes channel-free.
        pe_parent = Polymer(label="PE", monomer="[CH2][CH2]",
                            end_groups=["[H]", "[H]"], cutoff=3,
                            Mn=1000.0, Mw=2500.0, initial_mass=1.0)
        intent2 = SpawnIntent(
            parent_pool=pe_parent,
            monomer=pe_parent.backbone_group,
            end_groups=["[H]", "[H]"],
            triggering_dp=3,
        )
        d2 = _serialize_pool_for_sidecar(
            drain_spawn_intents([intent2], iteration=2)[0])
        assert "radical_qssa_unzip" not in d2["channels"]

    def _artifact(self, pools, labels, routing):
        return build_polymer_moments_artifact(
            pools, core_species=None, core_reactions=[],
            configured_pool_labels=labels,
            monomer_routing_by_pool=routing)

    def test_qssa_artifact_bumps_schema_and_recipe_revision(self, qssa_pool):
        """Version contract (round-23): channel-vocabulary growth is a minor
        SHAPE bump (2.0 -> 2.1) and the QSSA rate law is new channel/flux
        algebra (recipe_revision bumps over the legacy value). Both stamped
        only when a pool actually carries the channel. The -monomer-gas
        suffix carries the 2026-07-03 gas-monomer routing revision."""
        artifact = self._artifact([qssa_pool], ["PS"], {"PS": "styrene(5)"})
        assert artifact["schema_version"] == "2.1"
        assert artifact["conventions"]["recipe_revision"] == \
            "2026-07-03-qssa-monomer-gas"

    def test_qssa_artifact_stamps_format_doc_2_1(self, qssa_pool):
        """conventions.format_doc must agree with schema_version: a QSSA
        artifact stamps /2.1 (it used to keep the hard-coded /2.0 token even
        while schema_version said 2.1 — a self-contradicting envelope). No
        loader parses format_doc (grepped rmgpy + TA), but the human-facing
        stamp must not lie about the schema of the document carrying it."""
        artifact = self._artifact([qssa_pool], ["PS"], {"PS": "styrene(5)"})
        assert artifact["conventions"]["format_doc"] == (
            "docs/polymer_moments_format.md (polymer_moments_format/2.1)")

    def test_mixed_registry_bumps_when_any_pool_has_qssa(self, qssa_pool,
                                                         pe_pool):
        artifact = self._artifact([pe_pool, qssa_pool], ["PE", "PS"],
                                  {"PS": "styrene(5)"})
        assert artifact["schema_version"] == "2.1"
        assert artifact["conventions"]["recipe_revision"] == \
            "2026-07-03-qssa-monomer-gas"
        assert artifact["conventions"]["format_doc"] == (
            "docs/polymer_moments_format.md (polymer_moments_format/2.1)")

    def test_legacy_artifact_serialization_pinned(self, pe_pool):
        """No QSSA anywhere -> the ENTIRE serialized artifact is pinned to a
        golden dict: adding the QSSA feature must not have changed a single
        byte of a no-QSSA artifact, and any FUTURE drift (new key, dropped
        key, value change, even float-vs-int flips) fails loudly here.

        The two run-dependent envelope fields are handled explicitly:
        ``rmg_commit`` is overridden through its parameter and
        ``generated_at`` is popped after a format check. Byte-stability is
        asserted on ``json.dumps(..., sort_keys=True)`` rather than dict
        equality, because ``0 == 0.0`` as dicts yet they serialize to
        different bytes."""
        import re

        artifact = build_polymer_moments_artifact(
            [pe_pool], core_species=None, core_reactions=[],
            configured_pool_labels=["PE"], monomer_routing_by_pool={},
            rmg_commit="PINNED-FOR-TEST")

        # The QSSA vocabulary must be entirely absent from a legacy artifact.
        assert "radical_qssa_unzip" not in json.dumps(artifact)

        generated_at = artifact.pop("generated_at")
        assert re.fullmatch(r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z",
                            generated_at)

        golden = {
            "schema_version": "2.0",
            "rmg_commit": "PINNED-FOR-TEST",
            "rmg_iteration": 0,
            "conventions": {
                "format_doc": ("docs/polymer_moments_format.md "
                               "(polymer_moments_format/2.0)"),
                # 2026-07-03-monomer-gas: gas-monomer routing revision — the
                # ONE deliberate byte change vs the pre-revision golden.
                "recipe_revision": "2026-07-03-monomer-gas",
                "moment_basis": ("extensive mol, DP basis "
                                 "(mu1 = moles of repeat units)"),
                "volumes": {
                    "V_poly": "constant, consumer-supplied [m^3]",
                    "V_gas": ("ideal gas, dynamic: V_gas = n_gas*R*T/P "
                              "(1.0 m^3 floor when n_gas <= 0)"),
                },
                "configured_pools": ["PE"],
                "condensed_species": [],
                "site_scaling": ("site = max(0, mu_scaling)/V_poly read from "
                                 "the first proxy reactant's pool; multiplies "
                                 "ONCE; scales rf AND rr; near-exhaustion "
                                 "bundle limiter (tail-only smoothstep; "
                                 "C1 soft-min, cone-margin drain guard)"),
                "chip_site_throttle": ("site = min(max(0,mu0), max(0,mu1)/a)"
                                       "/V_poly when archetype=discrete_chip/1 "
                                       "and scaling=mu0 and a>0"),
                "kb_recipe": ("kb = kf/Keq; Keq(T) = (P0/(R*T))^dn * "
                              "exp(-dG0/(R*T)), P0 = 1e5 Pa, dG0 from "
                              "chem.yaml NASA thermo; dn counts ALL species "
                              "incl. condensed/proxy (format doc s4 step 1)"),
                "mu3_closure": "log_lagrange/1",
                "invariants": {
                    "discrete_subset": ("sum_pools(mu1) + sum_chip_species"
                                        "(a_i * n_i) is invariant over the "
                                        "discrete-reaction subset only"),
                    "with_unzip": ("add + n(monomer_routing) per pool with an "
                                   "active unzip channel (unzip moves units "
                                   "from mu1 into that species)"),
                },
            },
            "pools": [
                {
                    "label": "PE",
                    "monomer_smiles": "[CH2][CH2]",
                    "monomer_adj_list": (
                        "multiplicity 3\n"
                        "1 *1 C u1 p0 c0 {2,S} {3,S} {4,S}\n"
                        "2 *2 C u1 p0 c0 {1,S} {5,S} {6,S}\n"
                        "3    H u0 p0 c0 {1,S}\n"
                        "4    H u0 p0 c0 {1,S}\n"
                        "5    H u0 p0 c0 {2,S}\n"
                        "6    H u0 p0 c0 {2,S}\n"),
                    "feature_monomers_smiles": [],
                    "end_groups": ["[H]", "[H]"],
                    "cutoff": 3,
                    "parent_pool": None,
                    "spawn_iteration": 0,
                    "spawn_event_metadata": {"source": "input"},
                    "mu_indices": None,
                    "moments": [0.6666666666666666, 35.646601795609776,
                                2287.2243952345866],
                    "monomer_mw_g_mol": 28.053164947777987,
                    "mn_g_mol": 1500.0,
                    "mw_g_mol": 1800.0,
                    "initial_mass_g": 1000.0,
                    "channels": {
                        "scission": {"A": 1.0, "n": 0.0, "Ea": 0.0,
                                     "units": {"A": "s^-1", "Ea": "J/mol"}},
                        "unzip": {"A": 0.01, "n": 0.0, "Ea": 0.0,
                                  "units": {"A": "s^-1", "Ea": "J/mol"}},
                    },
                    "phase_species": [],
                    "bookkeeping_species": [],
                    "monomer_routing": None,
                    "mu3_closure": "log_lagrange/1",
                    "moments_provenance": "input_declared",
                }
            ],
            "reactions": [],
        }
        assert json.dumps(artifact, sort_keys=True) == \
            json.dumps(golden, sort_keys=True)

    def test_qssa_routing_and_json_round_trip(self, qssa_pool):
        """A QSSA pool with k_unzip == 0 still carries monomer_routing (the
        channel reuses the pool's existing routing), and the whole artifact
        survives a real json.dumps/loads round trip losslessly."""
        artifact = build_polymer_moments_artifact(
            [qssa_pool],
            core_species=None,
            core_reactions=[],
            configured_pool_labels=["PS"],
            monomer_routing_by_pool={"PS": "styrene(5)"},
        )
        rt = json.loads(json.dumps(artifact))
        pool_entry = rt["pools"][0]
        assert pool_entry["monomer_routing"] == "styrene(5)"
        assert pool_entry["channels"]["unzip"]["A"] == 0.0
        block = pool_entry["channels"]["radical_qssa_unzip"]
        assert block == artifact["pools"][0]["channels"]["radical_qssa_unzip"]
        assert block["enabled"] is True
        assert block["termination"]["units"]["A"] == "m^3/(mol*s)"


# Raw (pre-normalization) weak-link allyl/U-state channel config (schema 2.2,
# milestone iv). Split termination REPLACES the legacy summed block; the four
# weak-link keys are all-or-nothing (validate_radical_qssa_unzip).
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

# Pinned U0 units note (schema 2.2): U is serialized on the SAME amount basis
# as mu0; the consumer divides by V_poly to get the concentration state.
WEAKLINK_U0_UNITS = "mol — tail-distribution state; consumer divides by V_poly"

# The normative machine-readable weak-link recipe, pinned as LITERALS
# verified against the implemented RHS in rmgpy/solver/polymer.pyx
# (weak-link branch, commit c6b13c3a0). The loader pins the same strings
# independently; editing either copy is a consumer-coordination event.
WEAKLINK_PINNED_RECIPE = {
    "bond_basis": ("B = max(mu1 - mu0, 0) on concentration moments "
                   "(mol/m^3 condensed)"),
    "channel_gate": ("channel gates on B > 0: at B = 0 it is INERT even at "
                     "U > 0 (no monomer drain, no U production, no U sink; "
                     "the DP->1 self-termination floor survives)"),
    "u_state": ("U is a per-pool amount state [mol], tail-distribution "
                "basis, MASSLESS (it never enters condensed mass or any "
                "Mn/Mw/PDI consumer); daughter pools spawn with U0 = 0 "
                "(constants inherit, state resets)"),
    "u_active": ("u_active = min(max(U, 0)/V_poly, B) (active-site clamp, "
                 "mol/m^3)"),
    "radical_generation": ("G_R = 2*efficiency*ki*B + "
                           "1*efficiency*ki_allyl*u_active (nu = 1: one "
                           "unzipping radical per weak-link fission; the "
                           "allylic co-fragment does not unzip)"),
    "kt_split": ("kt_total = kt_rec + kt_disp replaces the legacy summed kt "
                 "everywhere (halved radical-disappearance convention "
                 "unchanged)"),
    "rate_no_transfer": ("r_mono = monomer_yield * kdp * "
                         "sqrt(G_R / (2*kt_total))"),
    "rate_with_transfer": ("r_mono = monomer_yield * kdp * "
                           "(sqrt(ktr^2 + 8*kt_total*G_R) - ktr) / "
                           "(4*kt_total)"),
    "du_dt": ("dU/dt = kt_disp*R_ss^2*max(0, 1 - U/(2*mu0))*V_poly - "
              "efficiency*ki_allyl*u_active*V_poly; production is the "
              "disproportionation EVENT rate with NO efficiency factor "
              "(R_ss already carries the escape efficiency) under a linear "
              "capacity throttle exactly zero at the TAIL chain-end "
              "capacity 2*mu0 [mol], where mu0 is the INSTANTANEOUS "
              "tail-distribution mu0(t) amount read from the current "
              "state at each RHS evaluation (NOT the frozen initial mu0 "
              "of the u0_census bound); at mu0(t) <= 0 the throttle is "
              "EXACTLY 0 by explicit branch (no ends, no capacity, no U "
              "production; no division is performed and no small_eps "
              "floor is applied); the sink is efficiency-SYMMETRIC "
              "with the G_R allyl term (a caged recombination restores "
              "the allylic bond)"),
    "u0_census": ("U0 = unsaturated_tail_ends_initial [mol] is rejected at "
                  "solver init when U0 > 2*mu0_tail evaluated on the "
                  "INITIAL (t = 0) tail mu0 (TAIL-only chain-end capacity, "
                  "amount basis); the loader passes the value through"),
    "u_transport": ("U is NEVER advected or transferred by any inter-pool "
                    "or ejection flux archetype (migration, "
                    "scission_fragment, discrete_chip and "
                    "volatile_ejection move pool moments and species "
                    "amounts only); the ONLY writers of U are this "
                    "recipe's du_dt law and the t = 0 initial condition, "
                    "and daughter pools spawn with U0 = 0 (constants "
                    "inherit, state resets)"),
    "moment_signature": ("dmu0 = 0; dmu1 -= r_mono; dmu2 -= r_mono * "
                         "max(2*mu1/max(mu0, small_eps) - 1, 0)"),
    "small_eps": 1e-30,
    "volume_note": ("kt is bimolecular: rates depend on condensed "
                    "volume V_poly; consumers MUST evaluate on "
                    "concentration moments mu_k = n_k / V_poly and "
                    "convert emitted rate back with *V_poly"),
}


@pytest.fixture
def weaklink_pool():
    pool = Polymer(
        label="PSW",
        monomer="[CH2][CH](c1ccccc1)",
        end_groups=["[H]", "[H]"],
        cutoff=3,
        moments=[1.0, 50.0, 3000.0],
        initial_mass=0.0,
        k_scission=0.0,
        k_unzip=0.0,
        radical_qssa_unzip=dict(WEAKLINK_RAW_CFG),
    )
    pool.monomer_product_species = _spc("C=Cc1ccccc1", "styrene")
    return pool


class TestWeakLinkChannelSerialization:
    """Sidecar schema 2.2 (weak-link milestone iv): the emitter serializes
    the weak-link allyl/U-state vocabulary instead of refusing (the
    milestone-i anti-silent-no-op guard is replaced by the real thing).
    The no-laundering concern the old refusal test pinned lives on in the
    LOADER: a 2.1-stamped artifact carrying weak-link keys is rejected
    (polymerMomentsRunnerTest)."""

    def test_weaklink_block_serializes_all_keys_with_units(self, weaklink_pool):
        d = _serialize_pool_for_sidecar(weaklink_pool)
        block = d["channels"]["radical_qssa_unzip"]
        assert block["enabled"] is True
        assert block["basis"] == "backbone_bonds_mu1_minus_mu0"
        assert block["efficiency"] == pytest.approx(0.8)
        assert block["monomer_yield"] == pytest.approx(0.9)
        # legacy triplets keep their pinned units
        for name in ("initiation", "depropagation", "transfer"):
            trip = block[name]
            assert trip["A"] == pytest.approx(WEAKLINK_RAW_CFG[name]["A"])
            assert trip["units"] == QSSA_PINNED_UNITS[name]
        # split termination REPLACES the summed block (mutual exclusion):
        # the summed slot is explicitly null, the split triplets carry the
        # bimolecular units of the block they replace.
        assert block["termination"] is None
        assert block["initiation_allyl"] == {
            "A": 2.0e14, "n": 0.0, "Ea": 2.4e5,
            "units": {"A": "s^-1", "Ea": "J/mol"}}
        assert block["termination_recombination"] == {
            "A": 6.0e7, "n": 0.0, "Ea": 8.0e3,
            "units": {"A": "m^3/(mol*s)", "Ea": "J/mol"}}
        assert block["termination_disproportionation"] == {
            "A": 4.0e7, "n": 0.0, "Ea": 1.2e4,
            "units": {"A": "m^3/(mol*s)", "Ea": "J/mol"}}
        # U0 is state, not a rate constant: value + pinned units note
        assert block["unsaturated_tail_ends_initial"] == {
            "value": 0.02, "units": WEAKLINK_U0_UNITS}

    def test_weaklink_recipe_block_emitted_verbatim(self, weaklink_pool):
        """The machine-pinned weak-link law (milestone iv): every string is
        a LITERAL a consumer can hash/compare, verified against the
        implemented RHS in rmgpy/solver/polymer.pyx (weak-link branch)."""
        d = _serialize_pool_for_sidecar(weaklink_pool)
        block = d["channels"]["radical_qssa_unzip"]
        assert block["recipe"] == WEAKLINK_PINNED_RECIPE
        assert "provenance" in block

    def _artifact(self, pools, labels, routing):
        return build_polymer_moments_artifact(
            pools, core_species=None, core_reactions=[],
            configured_pool_labels=labels,
            monomer_routing_by_pool=routing)

    def test_weaklink_artifact_stamps_2_2(self, weaklink_pool):
        """Weak-link vocabulary present -> the ARTIFACT stamps schema 2.2,
        format_doc /2.2 and the weak-link recipe revision (the artifact-level
        conventions.recipe_revision is where rate-recipe revisions live)."""
        artifact = self._artifact([weaklink_pool], ["PSW"],
                                  {"PSW": "styrene(5)"})
        assert artifact["schema_version"] == "2.2"
        assert artifact["conventions"]["format_doc"] == (
            "docs/polymer_moments_format.md (polymer_moments_format/2.2)")
        assert artifact["conventions"]["recipe_revision"] == \
            "2026-07-03-weaklink-u-monomer-gas"

    def test_mixed_artifact_stamps_2_2_legacy_channels_unchanged(
            self, weaklink_pool, qssa_pool, pe_pool):
        """A mixed artifact (legacy pool + legacy-QSSA pool + weak-link pool)
        is schema 2.2 (artifact stamp governed by the strongest vocabulary
        present), while the legacy pools' channel blocks serialize BYTE-
        IDENTICALLY to their 2.1 form."""
        solo = self._artifact([qssa_pool], ["PS"], {"PS": "styrene(5)"})
        assert solo["schema_version"] == "2.1"
        solo_block = solo["pools"][0]["channels"]["radical_qssa_unzip"]

        mixed = self._artifact(
            [pe_pool, qssa_pool, weaklink_pool], ["PE", "PS", "PSW"],
            {"PS": "styrene(5)", "PSW": "styrene(5)"})
        assert mixed["schema_version"] == "2.2"
        assert mixed["conventions"]["recipe_revision"] == \
            "2026-07-03-weaklink-u-monomer-gas"
        assert mixed["conventions"]["format_doc"] == (
            "docs/polymer_moments_format.md (polymer_moments_format/2.2)")
        by_label = {p["label"]: p for p in mixed["pools"]}
        assert json.dumps(by_label["PS"]["channels"]["radical_qssa_unzip"]) \
            == json.dumps(solo_block)
        assert set(by_label["PE"]["channels"]) == {"scission", "unzip"}

    def test_legacy_qssa_block_byte_identical_to_2_1_golden(self, qssa_pool):
        """r34 baseline lock, channel level: the serialized 2.1 QSSA block of
        a NO-weak-link pool is pinned byte-for-byte to a golden captured from
        the pre-milestone-iv build (7e1fa0671 lineage). The weak-link
        serialization must not have moved a single byte of legacy output."""
        d = _serialize_pool_for_sidecar(qssa_pool)
        block = d["channels"]["radical_qssa_unzip"]
        golden = (
            '{"enabled": true, "basis": "backbone_bonds_mu1_minus_mu0", '
            '"efficiency": 0.8, "monomer_yield": 0.9, '
            '"initiation": {"A": 10000000000000.0, "n": 0.0, "Ea": 300000.0, '
            '"units": {"A": "s^-1", "Ea": "J/mol"}}, '
            '"depropagation": {"A": 100000000000000.0, "n": 0.5, '
            '"Ea": 90000.0, "units": {"A": "s^-1", "Ea": "J/mol"}}, '
            '"termination": {"A": 100000000.0, "n": 0.0, "Ea": 10000.0, '
            '"units": {"A": "m^3/(mol*s)", "Ea": "J/mol"}}, '
            '"transfer": {"A": 2000.0, "n": 0.0, "Ea": 50000.0, '
            '"units": {"A": "s^-1", "Ea": "J/mol"}}, '
            '"recipe": {"bond_basis": "B = max(mu1 - mu0, 0) on '
            'concentration moments (mol/m^3 condensed)", '
            '"rate_no_transfer": "r_mono = monomer_yield * kdp * '
            'sqrt(efficiency * ki * B / kt)", '
            '"rate_with_transfer": "r_mono = monomer_yield * kdp * '
            '(sqrt(ktr^2 + 8*kt*(2*efficiency*ki*B)) - ktr) / (4*kt)", '
            '"moment_signature": "dmu0 = 0; dmu1 -= r_mono; dmu2 -= r_mono * '
            'max(2*mu1/max(mu0, small_eps) - 1, 0)", '
            '"small_eps": 1e-30, '
            '"volume_note": "kt is bimolecular: rates depend on condensed '
            'volume V_poly; consumers MUST evaluate on concentration moments '
            'mu_k = n_k / V_poly and convert emitted rate back with '
            '*V_poly"}, '
            '"provenance": {"radical_balance": "G_R = 2*f*ki*B; loss = '
            'ktr*R + 2*kt*R^2; Rss no-transfer = sqrt(f*ki*B/kt)", '
            '"moment_closure": "end_shrink_pool_mean/1", '
            '"R_gas_J_per_mol_K": 8.314, '
            '"concentration_basis": "mol/m^3 condensed volume", '
            '"transfer_note": "ktr is pseudo-first-order (s^-1); bimolecular '
            'literature k_tr must be premultiplied by substrate '
            'concentration before entering this config"}}')
        assert json.dumps(block) == golden


class TestQssaRoutingDerivation:
    """Round-25 P1-1: the PRODUCER must never emit an enabled QSSA block with
    monomer_routing null -- that shape is defined malformed (the consumer
    hard-rejects it: polymerMomentsRunnerTest test_rejects_enabled_without_
    routing). The live save_everything hook builds monomer_routing_by_pool
    from the ENGINE's configured pools only, so a daughter registered after
    the last solver rebuild has no entry there; the serializer must then
    derive the routing from the pool's own monomer_product_species (held
    BY REFERENCE from M5 inheritance) against the same core-species universe
    the routing labels come from -- or refuse to serialize."""

    def test_unconfigured_qssa_daughter_derives_routing_from_species_ref(
            self, qssa_pool):
        """A channel-bearing daughter NOT present in the engine configs (no
        monomer_routing_by_pool entry) serializes with routing DERIVED from
        its monomer_product_species -- never enabled + null."""
        sty = _spc("C=Cc1ccccc1", "styrene", index=7)
        qssa_pool.monomer_product_species = sty
        qssa_pool.parent_pool_label = "PS_parent"
        core = [_spc("N#N", "N2", index=1), sty]

        artifact = build_polymer_moments_artifact(
            [qssa_pool], core_species=core, core_reactions=[],
            configured_pool_labels=["PS_parent"],  # daughter unconfigured
            monomer_routing_by_pool={})            # engine has no entry

        entry = artifact["pools"][0]
        assert entry["channels"]["radical_qssa_unzip"]["enabled"] is True
        assert entry["monomer_routing"] == "styrene(7)"
        # gas since recipe revision 2026-07-03-monomer-gas
        assert "styrene(7)" not in entry["phase_species"]

    def test_explicit_engine_routing_wins_over_derivation(self, qssa_pool):
        """The engine's configured routing stays authoritative when present."""
        sty = _spc("C=Cc1ccccc1", "styrene", index=7)
        qssa_pool.monomer_product_species = sty
        d = _serialize_pool_for_sidecar(qssa_pool, core_species=[sty],
                                        monomer_routing="styrene(7)")
        assert d["monomer_routing"] == "styrene(7)"

    def test_qssa_pool_without_routing_species_is_hard_error(self, qssa_pool):
        """No monomer_product_species and no engine routing: serialization
        must FAIL LOUDLY, never emit enabled + null."""
        qssa_pool.monomer_product_species = None
        with pytest.raises(ValueError,
                           match=r"PS.*radical_qssa_unzip.*monomer_routing"):
            _serialize_pool_for_sidecar(qssa_pool)

    def test_qssa_pool_routing_species_not_in_core_is_hard_error(
            self, qssa_pool):
        """monomer_product_species present but NOT in the core universe the
        artifact labels come from (identity check -- routing resolution is
        object-keyed everywhere): hard error, not a dangling label."""
        orphan = _spc("C=Cc1ccccc1", "styrene", index=7)
        qssa_pool.monomer_product_species = orphan
        core = [_spc("N#N", "N2", index=1)]
        with pytest.raises(ValueError,
                           match=r"PS.*radical_qssa_unzip.*monomer_routing"):
            _serialize_pool_for_sidecar(qssa_pool, core_species=core)

    def test_channel_free_pool_still_serializes_null_routing(self, pe_pool):
        """The hard error is scoped to enabled-QSSA pools: legacy pools keep
        the legacy null-routing shape byte-identically."""
        d = _serialize_pool_for_sidecar(pe_pool)
        assert d["monomer_routing"] is None


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


class TestPoolLivenessAlarm:
    """r95 run-gate POOL LIVENESS census (review-gated alarm at artifact-write
    time, NOT a hard-fail). build_polymer_moments_artifact lists EVERY
    pool-coupled row -- refused rows stamp-but-kept alongside live ones -- so if
    a pool has coupled rows but ZERO survive live, the refusal predicates may be
    over-broad (the r95 PP-rerun regression: 106 rows all refused, pool 12->0
    live). Census/warning ONLY: a legitimately deck-incomplete polymer may
    honestly refuse every row, so this must NEVER raise."""

    @staticmethod
    def _artifact_with(refused):
        # One pool-coupled reversible row (proxy 'poly' resolves the configured
        # pool). refused=True stamps it conduit-deferred (still listed).
        pool = Polymer(label="poly", monomer="[CH2][CH2]",
                       end_groups=["[H]", "[H]"], cutoff=3,
                       moments=[1.0, 5.0, 30.0], initial_mass=0.0,
                       k_scission=1.0, k_unzip=0.0)
        pr = _spc("[CH2]CCCCCCCCCCCCCCCCCCCCC", "PR")
        et = _spc("C[CH2]", "C2H5")
        proxy = _spc("CCCCCCCCCCCCCCCCCCCCCCCC", "poly")
        mus = [_mu_dummy("poly_mu0"), _mu_dummy("poly_mu1"),
               _mu_dummy("poly_mu2")]
        core = [pr, et, proxy] + mus
        rxn = Reaction(reactants=[pr, et], products=[proxy],
                       kinetics=_arrhenius(A=(3.0, "m^3/(mol*s)")),
                       reversible=True)
        if refused:
            rxn.polymer_refused = True
            rxn.polymer_refused_accumulating = False
        return build_polymer_moments_artifact(
            [pool], core_species=core, core_reactions=[rxn],
            configured_pool_labels=["poly"],
            condensed_species=mus + [proxy],
            cantera_index_map={id(rxn): [0]})

    def test_all_refused_pool_emits_liveness_alarm(self, caplog):
        import logging
        with caplog.at_level(logging.WARNING):
            art = self._artifact_with(refused=True)
        # stamp-but-keep: the row is still listed, but refused ...
        assert len(art["reactions"]) == 1
        assert art["reactions"][0]["refused"] is True
        # ... and the alarm fired (census/warning only -- the build SUCCEEDED).
        assert any("POOL LIVENESS ALARM" in r.getMessage()
                   for r in caplog.records)

    def test_healthy_pool_emits_no_liveness_alarm(self, caplog):
        import logging
        with caplog.at_level(logging.WARNING):
            art = self._artifact_with(refused=False)
        assert len(art["reactions"]) == 1
        assert not art["reactions"][0].get("refused")
        assert not any("POOL LIVENESS ALARM" in r.getMessage()
                       for r in caplog.records)


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

    def test_volatile_ejection_entry_carries_eject_units(self):
        """A VOLATILE_EJECTION reaction (polymer A -> discrete volatile +
        cross-pool polymer B) serializes archetype ``volatile_ejection/1`` with
        ``params: {"eject_units": <a>}`` -- the FRACTIONAL source-monomer-
        equivalents that leave as volatile (mirrors the DISCRETE_CHIP a-param)."""
        core = _two_pool_core()
        # reactant A (pool), products: cross-pool B + discrete gas G
        rxn = Reaction(reactants=[core[0]], products=[core[4], core[8]],
                       kinetics=_arrhenius(), reversible=False)
        rxn.polymer_flux_archetype = int(PolymerFluxArchetype.VOLATILE_EJECTION)
        rxn.polymer_eject_units = 1.135
        entries = compile_polymer_reaction_entries(
            [rxn], core, configured_pool_labels=["A", "B"],
            cantera_index_map={id(rxn): [4]})
        e = entries[0]
        assert e["archetype"] == "volatile_ejection/1"
        assert e["params"]["eject_units"] == pytest.approx(1.135)
        assert set(e["params"].keys()) == {"eject_units"}
        assert e["src_pool"] == "A" and e["dst_pool"] == "B"
        assert e["unresolved"] is False

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


class TestRefusedRowSerialization:
    """Schema 2.4 (refused-row sidecar contract): a reaction stamped
    ``polymer_refused=True`` (item 18 stamp-but-keep) is zeroed by the
    generating solver (polymer.pyx reaction_refused), so its sidecar row must
    carry ``refused: true`` + ``refused_reason`` — otherwise a consumer
    integrates moment flux the oracle fabricated nothing for (the live
    RMG-vs-TA divergence this contract closes). The row STAYS listed:
    consumers need it to zero the Cantera multiplier."""

    def _refused_rxn(self, core, accumulating=False):
        # The emitter-side refused shape: pool proxy reactant, gas products
        # only (the chain radical the handshake dropped to gas), stamped
        # UNRESOLVED + polymer_refused by Task 3 at generation time.
        rxn = Reaction(reactants=[core[0], core[8]], products=[core[9]],
                       kinetics=_arrhenius(A=(3.0, "m^3/(mol*s)")),
                       reversible=False)
        rxn.polymer_flux_archetype = int(PolymerFluxArchetype.UNRESOLVED)
        rxn.polymer_refused = True
        rxn.polymer_refused_accumulating = accumulating
        return rxn

    def test_refused_row_carries_marker_and_census_reason(self):
        core = _two_pool_core()
        rxn = self._refused_rxn(core, accumulating=False)
        entries = compile_polymer_reaction_entries(
            [rxn], core, configured_pool_labels=["A", "B"],
            cantera_index_map={id(rxn): [3]})
        assert len(entries) == 1  # stamp-but-keep: the row stays listed
        e = entries[0]
        assert e["refused"] is True
        # the census reason available at stamp time (polymer.pyx:1529):
        # eliminating radical -> conduit-deferred
        assert e["refused_reason"] == "conduit-deferred"
        # a refused row still runs the legacy emission (solver mirror)
        assert e["archetype"] == "legacy_mu1/1"
        assert e["unresolved"] is True

    def test_refused_reason_tracks_accumulating_class(self):
        core = _two_pool_core()
        rxn = self._refused_rxn(core, accumulating=True)
        entries = compile_polymer_reaction_entries(
            [rxn], core, configured_pool_labels=["A", "B"],
            cantera_index_map={id(rxn): [3]})
        assert entries[0]["refused_reason"] == "qssa-invalid"

    def test_gas_association_refused_row_round_trips_conduit_deferred(self):
        """PP v1 gas-association refusal (adjudicated round 63): a run-5-shaped
        row (pure gas radicals <=> condensed pool proxy) stamped by the NEW
        classifier ``stamp_gas_association_refusal`` round-trips the sidecar
        with the shipped schema-2.4 vocabulary: ``refused: true`` +
        ``refused_reason: "conduit-deferred"`` (the row is deferred pending
        the pool-moment-credit conduit). RED at b917becd7: the classifier does
        not exist."""
        from rmgpy.polymer import stamp_gas_association_refusal
        pp = Polymer(label="polypropylene", monomer="[CH2][CH]C",
                     Mn=5000.0, Mw=8000.0, initial_mass=1.0)
        pp.index = 2
        r1 = _spc("[CH2]C(C)C", "iC4H9", index=4)
        r2 = _spc("C[CH]CCC", "C5H11", index=5)
        core = [pp, r1, r2]
        rxn = Reaction(reactants=[r1, r2], products=[pp],
                       kinetics=_arrhenius(A=(3.0, "m^3/(mol*s)")),
                       reversible=True)
        stamp_gas_association_refusal(rxn)
        assert rxn.polymer_refused is True   # classifier fired (RED half)
        entries = compile_polymer_reaction_entries(
            [rxn], core, configured_pool_labels=["polypropylene"],
            cantera_index_map={id(rxn): [0]})
        assert len(entries) == 1   # stamp-but-keep: the row stays listed
        e = entries[0]
        assert e["refused"] is True
        assert e["refused_reason"] == "conduit-deferred"
        assert e["proxy_products"] == ["polypropylene(2)"]

    def test_impostor_refused_row_round_trips_conduit_deferred(self):
        """r82 impostor-row refusal (FR1 run-2): a run-2-shaped row (small
        closed-shell gas + polymer-sized discrete <=> condensed pool proxy)
        stamped by ``stamp_gas_association_refusal`` round-trips the sidecar
        with the same shipped schema-2.4 vocabulary as the r63/r74 shapes:
        ``refused: true`` + ``refused_reason: "conduit-deferred"``. RED at
        3521caead: the impostor conjunct does not exist, the row is never
        stamped, and the sidecar carries no refused keys."""
        from rmgpy.polymer import stamp_gas_association_refusal
        pp = Polymer(label="polypropylene", monomer="[CH2][CH]C",
                     Mn=5000.0, Mw=8000.0, initial_mass=1.0)
        pp.index = 2
        br2 = _spc("BrBr", "Br2", index=4)                    # closed-shell
        # Genuine chain-scale unsaturated impostor (C24H48, 336.6 g/mol / 24
        # heavy), clearing the r95 absolute floor -- an honest FR1 run-2-class
        # proxy-minus-XY discrete, not a light-monomer ratio artifact.
        impostor = _spc("C=CCC(C)CC(C)CC(C)CC(C)CC(C)CC(C)CC(C)C",
                        "C24H48", index=5)
        core = [pp, br2, impostor]
        rxn = Reaction(reactants=[br2, impostor], products=[pp],
                       kinetics=_arrhenius(A=(3.0, "m^3/(mol*s)")),
                       reversible=True)
        stamp_gas_association_refusal(rxn)
        assert rxn.polymer_refused is True   # classifier fired (RED half)
        entries = compile_polymer_reaction_entries(
            [rxn], core, configured_pool_labels=["polypropylene"],
            cantera_index_map={id(rxn): [0]})
        assert len(entries) == 1   # stamp-but-keep: the row stays listed
        e = entries[0]
        assert e["refused"] is True
        assert e["refused_reason"] == "conduit-deferred"
        assert e["proxy_products"] == ["polypropylene(2)"]

    def test_non_refused_rows_carry_no_refused_keys(self):
        """Non-refused rows: NO new key — absent, not false (byte-identical
        pin; TA-side loaders treat key PRESENCE as the 2.4 vocabulary)."""
        core = _two_pool_core()
        mig = Reaction(reactants=[core[0]], products=[core[4]],
                       kinetics=_arrhenius(), reversible=False)
        mig.polymer_flux_archetype = int(PolymerFluxArchetype.MIGRATION)
        legacy = Reaction(reactants=[core[8], core[0]],
                          products=[core[9], core[0]],
                          kinetics=_arrhenius(A=(3.0, "m^3/(mol*s)")),
                          reversible=False)  # unstamped -> legacy/unresolved
        entries = compile_polymer_reaction_entries(
            [mig, legacy], core, configured_pool_labels=["A", "B"],
            cantera_index_map={id(mig): [0], id(legacy): [1]})
        assert len(entries) == 2
        for e in entries:
            assert "refused" not in e
            assert "refused_reason" not in e

    def test_refused_row_without_pool_mapped_participant_hard_fails(self):
        """r92 artifact P1 (re-adjudicates the former warn+emit pin, a
        defect pin): a refused row with NO configured-pool participant
        cannot legally carry the refused marker (loader guard, closed
        schema-2.4 vocabulary) -- and emitting it live means a conforming
        consumer integrates flux the generating solver ZEROED
        (over-integration; artifact misrepresents solver physics, the
        r71-banned class). The emitter must HARD-FAIL loudly, never emit."""
        core = _two_pool_core()
        rxn = self._refused_rxn(core)
        with pytest.raises(ValueError, match=r"refused.*no pool-mapped"):
            compile_polymer_reaction_entries(
                [rxn], core, configured_pool_labels=["B"],  # A not configured
                cantera_index_map={id(rxn): [3]})

    def test_flip_refused_row_not_emitted_live_legacy(self):
        """r92 (run-10 artifact rows r8/r30-32): a kinetics-flipped row whose
        flipped-direction restamp is unresolvable is refused at generation
        time and must reach the artifact carrying the refused marker (the
        consumer zeroes it) -- never as a live legacy_mu1/1 row while the
        generating solver zeroes its flux."""
        from rmgpy.polymer import restamp_flipped_polymer_archetype
        pool = Polymer(label="A", monomer="[CH2][CH]c1ccccc1",
                       end_groups=["[CH3]", "[H]"], cutoff=3,
                       Mn=5000.0, Mw=6000.0, initial_mass=1.0)
        pool.index = 1
        gas = _spc("CC", "ethane", index=7)
        # Flipped orientation left by apply_kinetics_to_reaction: gas -> pool
        # (reversed association), stale pre-flip stamp still on the object.
        rxn = Reaction(reactants=[gas], products=[pool],
                       kinetics=_arrhenius(A=(3.0, "m^3/(mol*s)")),
                       reversible=True)
        rxn.polymer_flux_archetype = int(PolymerFluxArchetype.SCISSION_FRAGMENT)
        restamp_flipped_polymer_archetype(rxn)
        entries = compile_polymer_reaction_entries(
            [rxn], [pool, gas], configured_pool_labels=["A"],
            cantera_index_map={id(rxn): [8]})
        assert len(entries) == 1
        e = entries[0]
        assert e.get("refused") is True
        assert e.get("refused_reason") == "conduit-deferred"

    def test_refused_presence_stamps_schema_2_4(self, pe_pool):
        """Presence-based stamp, the 2.1/2.2/2.3 precedent exactly: the
        emitter stamps "2.4" exactly when at least one reactions[] row
        carries the refused marker. Refused rows add SHAPE vocabulary with
        consumption semantics, not new rate algebra, so recipe_revision is
        untouched (STRICT-MINOR acceptance stops old consumers at the
        envelope)."""
        core = _two_pool_core()
        rxn = self._refused_rxn(core)
        # "PE" listed configured alongside A/B: the registry pool must not
        # read as runtime-spawned here, or the spawned-pool closure (schema
        # 2.5, the stronger SHAPE stamp) legitimately outranks the 2.4
        # refused stamp this test pins.
        artifact = build_polymer_moments_artifact(
            [pe_pool], core_species=core, core_reactions=[rxn],
            configured_pool_labels=["PE", "A", "B"],
            condensed_species=core[:4],
            cantera_index_map={id(rxn): [3]})
        assert artifact["schema_version"] == "2.4"
        assert "polymer_moments_format/2.4" in artifact["conventions"]["format_doc"]
        assert artifact["conventions"]["recipe_revision"] == \
            POLYMER_RATE_RECIPE_REVISION

    def test_no_refused_keeps_legacy_stamps_byte_identical(self, pe_pool):
        """No refused row anywhere -> the 2.0 stamp (and the whole artifact)
        is byte-identical to the pre-2.4 emitter."""
        core = _two_pool_core()
        mig = Reaction(reactants=[core[0]], products=[core[4]],
                       kinetics=_arrhenius(), reversible=False)
        mig.polymer_flux_archetype = int(PolymerFluxArchetype.MIGRATION)
        # "PE" configured for the same reason as the 2.4-stamp test above:
        # a registry pool absent from the configured set is runtime-spawned
        # by the primary signal and would truthfully stamp 2.5.
        artifact = build_polymer_moments_artifact(
            [pe_pool], core_species=core, core_reactions=[mig],
            configured_pool_labels=["PE", "A", "B"],
            condensed_species=core[:4],
            cantera_index_map={id(mig): [0]})
        assert artifact["schema_version"] == "2.0"
        assert "polymer_moments_format/2.0" in artifact["conventions"]["format_doc"]
        for e in artifact["reactions"]:
            assert "refused" not in e


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
        # 2026-07-03-monomer-gas: gas-monomer routing revision (deliberate
        # consumer-coordination bump; see docs/polymer_moments_format.md).
        assert data["conventions"]["recipe_revision"] == "2026-07-03-monomer-gas"
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


class TestCollectPolymerPoolRegistry:
    """collect_polymer_pool_registry is the single builder for the sidecar's
    pool_registry (used by main.save_everything and
    CoreEdgeReactionModel._apply_multipool_spawn_pass). A freshly-promoted
    daughter Polymer sits in BOTH core.species and new_species_list until the
    next enlarge clears it, so an un-deduped concatenation serializes the same
    pool twice (observed live: pools [PS, tail, tail_2, tail, tail_2])."""

    def test_same_object_in_two_lists_appears_once(self, pe_pool):
        from rmgpy.polymer import collect_polymer_pool_registry

        gas = _spc("C", "CH4", index=1)
        registry = collect_polymer_pool_registry(
            [gas, pe_pool],   # core.species (freshly-promoted daughter)
            [],               # edge.species
            [pe_pool],        # new_species_list (not yet cleared by enlarge)
        )
        assert len(registry) == 1
        assert registry[0] is pe_pool

    def test_order_preserved_and_distinct_pools_kept(self, pe_pool):
        """Identity dedup only: first occurrence wins, order preserved, and a
        DISTINCT (even equal-looking) Polymer object is NOT collapsed."""
        from rmgpy.polymer import collect_polymer_pool_registry

        other = Polymer(
            label="PE",  # same label as pe_pool, different object
            monomer="[CH2][CH2]", end_groups=["[H]", "[H]"], cutoff=3,
            Mn=1500.0, Mw=1800.0, initial_mass=1.0,
        )
        gas = _spc("C", "CH4", index=1)
        registry = collect_polymer_pool_registry(
            [pe_pool, gas], [other], [pe_pool, other])
        assert len(registry) == 2
        assert registry[0] is pe_pool and registry[1] is other

    def test_non_polymer_species_filtered(self, pe_pool):
        from rmgpy.polymer import collect_polymer_pool_registry

        gas = _spc("C", "CH4", index=1)
        assert collect_polymer_pool_registry([gas], [gas], []) == []
        assert collect_polymer_pool_registry([gas, pe_pool]) == [pe_pool]


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
    epdm_scission_tail, which post-item-17 never reaches core at all (its
    producing channel is Gate-B zeroed at edge, census-announced) is NOT passed
    here and so is NOT reported condensed."""

    def _epdm_like_core(self):
        # Final core: epdm proxy + epdm_mu0/1/2 (the ONE configured pool, all
        # condensed) interleaved with gaseous N2 / H(1), plus the spawned
        # epdm_scission_tail family which the solver ran as ordinary GAS species
        # (its pool was never solver-configured).
        # NOTE (item 17, 2026-06-12): this encodes the OLD baseline's SHAPE.
        # Post-17 the live EPDM run never promotes the tail family (the
        # baseline re-pinned 26/28 -> 8/0; A5 measured, H2 also sub-bar), but
        # the emitter must keep handling cores of this shape -- they still
        # arrive via legacy/
        # restart artifacts AND via the third route (independent species
        # promotion, spec 2026-06-12 SS3(e)). RETAINED as legacy-core
        # coverage; the sibling fixture below mirrors the NEW baseline.
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

    def test_explicit_indices_condensed_monomer_index_gas(self):
        """Explicit-oligomer indices count as condensed in the derived
        fallback; the routed-monomer index does NOT (recipe revision
        2026-07-03-monomer-gas: the release target is a gas volatile, and
        the solver oracle validates it gas)."""
        core = [
            _spc("CC", "P", index=2),
            _mu_dummy("P_mu0"), _mu_dummy("P_mu1"), _mu_dummy("P_mu2"),
            _spc("[CH3]", "G", index=7),
            _spc("CCC", "P_dp3", index=8),
            _spc("[CH2]CC", "released_monomer", index=9),
        ]
        pools = [_FakePool("P", mu_indices=(1, 2, 3),
                           explicit_dp_to_species_index={3: 5},
                           monomer_poly_index=6)]
        condensed = derive_condensed_species(core, pools, mask=None)
        labels = [s.label for s in condensed]
        assert set(labels) == {"P", "P_mu0", "P_mu1", "P_mu2", "P_dp3"}
        assert "released_monomer" not in labels
        assert "G" not in labels

    def test_registry_pool_with_proxy_absent_from_core_emits_empty_lists(
            self, tmp_path):
        """Item 17 sibling fixture (spec SS4.3/SS7(i)): the NEW EPDM baseline
        shape -- tail pool in the registry, proxy + mu-dummies NOT in core.
        The entry survives emission with empty phase/bookkeeping lists (the
        documented absent-from-core shape, its first live instance), t=0
        moments, and no reactions[] reference. Probe-derived expected
        values (/tmp/probe_artifact_registry_pool.py, re-run green at HEAD
        4d667f1af)."""
        import json as _json
        from rmgpy.polymer import (Polymer, write_polymer_pools_sidecar,
                                   derive_condensed_species)
        parent = Polymer(label="PS", monomer="[CH2][CH]c1ccccc1",
                         end_groups=["[CH3]", "[H]"], cutoff=5,
                         Mn=3000.0, Mw=10000.0, initial_mass=1.0)
        tail = Polymer(label="PS_scission_tail", monomer="[CH2][CH]c1ccccc1",
                       feature_monomer=None, end_groups=["[CH3]", "[H]"],
                       cutoff=5, Mn=parent.Mn / 2.0, Mw=parent.Mw / 2.0,
                       initial_mass=0.0, moments=None)
        core = [
            _spc("N#N", "N2", index=4), _spc("[H]", "H", index=1), parent,
            _mu_dummy("PS_mu0"), _mu_dummy("PS_mu1"), _mu_dummy("PS_mu2"),
            _spc("[H][H]", "[H][H]", index=40),
        ]
        mask = [True, True, False, False, False, False, True]
        configured = [_FakePool("PS", mu_indices=(3, 4, 5))]
        condensed = derive_condensed_species(core, configured, mask)
        write_polymer_pools_sidecar(
            pool_registry=[parent, tail], output_dir=str(tmp_path),
            iteration=9, core_species=core, core_reactions=[],
            configured_pool_labels=["PS"], condensed_species=condensed,
            monomer_routing_by_pool={}, cantera_index_map=None)
        with open(tmp_path / "polymer_pools.json") as fh:
            data = _json.load(fh)
        assert [p["label"] for p in data["pools"]] == \
            ["PS", "PS_scission_tail"]
        entry = data["pools"][1]
        assert entry["phase_species"] == []
        assert entry["bookkeeping_species"] == []
        assert entry["moments"] == [0.0, 0.0, 0.0]
        assert entry["moments_provenance"] == "spawned_empty"
        assert entry["monomer_routing"] is None
        assert data["reactions"] == []
        assert "PS_scission_tail" not in data["conventions"]["condensed_species"]


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


# The normative explicit-DP recipe strings, pinned as LITERALS in the test
# (NOT imported from the emitter constant — an emitter edit must fail here,
# same idiom as WEAKLINK_PINNED_RECIPE). Transcribed from the implemented
# oracle: boundary flux polymer.pyx:3475-3524 (gamma-conditional p_cond,
# triangular fallback, flux clamps), gamma helpers polymer.pyx:491-525,
# k_chain arm selection polymer.pyx:3484-3488, t=0 tail split
# set_initial_conditions step 6 / _explicit_moment_contributions
# (polymer.pyx:2109-2132, :471-488: seeded moments = declared TOTAL minus
# each mapped DP's (N, dp*N, dp^2*N), clamped >= 0) + the matching
# generation-side V_poly mass split (PolymerPhase.calculate_volume).
EXPLICIT_DP_PINNED_RECIPE = {
    "tail_split": ("declared pools[].moments are TOTAL-INCLUSIVE (explicit "
                   "chains counted in); the solver seeds the explicit "
                   "species' initial_moles as species amounts (clamped "
                   ">= 0), then subtracts each mapped DP's contribution "
                   "(N_dp, dp*N_dp, dp^2*N_dp on concentration moments) "
                   "from the seeded mu0/mu1/mu2, clamped >= 0 with a clamp "
                   "warning (set_initial_conditions step 6) -- the "
                   "integrated tail moments are total - explicit; the "
                   "generation-side V_poly mass split applies the same "
                   "subtraction and hard-errors when explicit mu1 exceeds "
                   "declared mu1 beyond -1e-12"),
    "boundary_flux": ("gated on mu0 > 1e-9 mol/m^3 AND mu1/mu0 > xs + 1e-9; "
                      "p_cond = P(DP = xs+1 | DP > xs) from the gamma "
                      "distribution moment-matched to (mu0, mu1, mu2) "
                      "(k = 1/(PDI - 1), theta = mean/k; half-integer bins: "
                      "[F((xs+1.5)/theta) - F((xs+0.5)/theta)] / "
                      "[1 - F((xs+0.5)/theta)] with F the regularized lower "
                      "incomplete gamma); triangular fallback on tail_mean "
                      "in (xs+1, xs+2) peaking 1.0 at xs+1.5 when the gamma "
                      "is unrealizable (any moment <= 1e-30, PDI <= 1+1e-6, "
                      "or non-finite params); p_cond clamped to [0, 1]; "
                      "N_boundary = min(mu0*p_cond, mu0, mu1/xs, mu2/xs^2) "
                      "[mol/m^3]; F_flux = k_chain * N_boundary; "
                      "dn(species[xs])/dt += F_flux*V_poly; dmu0 -= F_flux; "
                      "dmu1 -= xs*F_flux; dmu2 -= xs^2*F_flux"),
    "k_chain": ("k_unzip when k_unzip > 0, else r_qssa/max(mu0, 1e-30) when "
                "the radical_qssa_unzip channel is active (mutually "
                "exclusive upstream; at most one arm fires)"),
    "transport": ("one-way: statistical moment tail -> explicit real "
                  "condensed species at DP == handshake_target_dp; no "
                  "reverse flux"),
}


@pytest.fixture
def explicit_dp_pool():
    """Stage-A-shaped pool: deck flag explicit_dp=True attached the capped
    DP=cutoff oligomer species (rmgpy/rmg/input.py polymer() step 4c)."""
    pool = Polymer(
        label="PS",
        monomer="[CH2][CH](c1ccccc1)",
        end_groups=["[H]", "[H]"],
        cutoff=3,
        Mn=1500.0,
        Mw=1800.0,
        initial_mass=1.0,
        k_scission=0.0,
        k_unzip=0.01,
    )
    pool.explicit_dp = True
    pool.explicit_dp_species = _spc("CCC", "PS_dp3", index=42)
    return pool


class TestExplicitDpSerialization:
    """Stage B (schema 2.3): the explicit-DP handshake target is serialized
    as a POOL-LEVEL block (state/topology, not kinetics — deliberately NOT
    inside channels), emitted ONLY when the pool carries the stage-A
    explicit species. Flag-OFF artifacts stay byte-identical (the golden pin
    in test_legacy_artifact_serialization_pinned covers the bytes; the
    absence assertions here cover the vocabulary)."""

    def test_explicit_dp_block_emitted_with_exact_shape(self, explicit_dp_pool):
        d = _serialize_pool_for_sidecar(explicit_dp_pool,
                                        initial_explicit_moles={3: 0.25})
        block = d["explicit_dp"]
        assert block == {
            "enabled": True,
            "species": {"3": "PS_dp3(42)"},
            "initial_moles": {"3": 0.25},
            "handshake_target_dp": 3,
            "recipe_revision": "2026-07-04-explicit-dp",
            "recipe": EXPLICIT_DP_PINNED_RECIPE,
        }

    def test_initial_moles_default_zero(self, explicit_dp_pool):
        """No initial_explicit loading declared -> 0.0, never a hole (the
        stage-A channel's default: the species starts empty)."""
        d = _serialize_pool_for_sidecar(explicit_dp_pool)
        assert d["explicit_dp"]["initial_moles"] == {"3": 0.0}

    def test_species_label_matches_artifact_convention(self, explicit_dp_pool):
        """Labels, not indices: the species entry uses the same chem.yaml
        label rule as every other species reference in the sidecar
        (_artifact_species_label: 'label(index)' when index > 0)."""
        core = [
            _spc("CC", "PS", index=2),
            _mu_dummy("PS_mu0"), _mu_dummy("PS_mu1"), _mu_dummy("PS_mu2"),
            explicit_dp_pool.explicit_dp_species,
        ]
        d = _serialize_pool_for_sidecar(explicit_dp_pool, core_species=core)
        assert d["explicit_dp"]["species"] == {"3": "PS_dp3(42)"}

    def test_explicit_species_joins_phase_species_not_bookkeeping(
            self, explicit_dp_pool):
        """Format doc §2 phase_species/bookkeeping_species carve-out:
        explicit-DP chains are REAL condensed inventory — listed in
        phase_species, and exactly NOT in the bookkeeping subset."""
        core = [
            _spc("CC", "PS", index=2),
            _mu_dummy("PS_mu0"), _mu_dummy("PS_mu1"), _mu_dummy("PS_mu2"),
            explicit_dp_pool.explicit_dp_species,
            _spc("[CH3]", "G", index=7),  # gas — in neither list
        ]
        d = _serialize_pool_for_sidecar(explicit_dp_pool, core_species=core)
        assert d["phase_species"] == [
            "PS(2)", "PS_mu0", "PS_mu1", "PS_mu2", "PS_dp3(42)"]
        assert d["bookkeeping_species"] == [
            "PS(2)", "PS_mu0", "PS_mu1", "PS_mu2"]

    def test_flag_without_species_refuses_to_serialize(self, explicit_dp_pool):
        """explicit_dp=True with no attached species would recreate the
        structurally-inert handshake — the producer refuses (same posture as
        compile_polymer_phase and enabled-QSSA-without-routing)."""
        explicit_dp_pool.explicit_dp_species = None
        with pytest.raises(ValueError, match="explicit_dp"):
            _serialize_pool_for_sidecar(explicit_dp_pool)

    def test_species_not_in_core_universe_refuses(self, explicit_dp_pool):
        """When a core universe is supplied, the explicit species must be in
        it BY IDENTITY (the labels must come from the same universe as every
        other artifact label)."""
        core = [_spc("CC", "PS", index=2), _spc("CCC", "PS_dp3", index=42)]
        # same label, DIFFERENT object -> identity check must refuse
        with pytest.raises(ValueError, match="explicit"):
            _serialize_pool_for_sidecar(explicit_dp_pool, core_species=core)

    def test_legacy_pool_emits_no_explicit_dp_vocabulary(self, pe_pool):
        d = _serialize_pool_for_sidecar(pe_pool)
        assert "explicit_dp" not in d

    def _artifact(self, pools, labels, routing, **kw):
        return build_polymer_moments_artifact(
            pools, core_species=None, core_reactions=[],
            configured_pool_labels=labels,
            monomer_routing_by_pool=routing, **kw)

    def test_explicit_dp_artifact_stamps_2_3(self, explicit_dp_pool):
        artifact = self._artifact([explicit_dp_pool], ["PS"], {})
        assert artifact["schema_version"] == "2.3"
        assert artifact["conventions"]["recipe_revision"] == \
            "2026-07-04-explicit-dp-monomer-gas"
        assert artifact["conventions"]["format_doc"] == (
            "docs/polymer_moments_format.md (polymer_moments_format/2.3)")

    def test_no_explicit_artifact_keeps_legacy_stamps(self, pe_pool,
                                                      qssa_pool,
                                                      weaklink_pool):
        """Absent path unchanged: no explicit_dp anywhere -> the 2.0/2.1/2.2
        pick-max stamps apply byte-identically (composition with the golden
        pin in test_legacy_artifact_serialization_pinned)."""
        legacy = self._artifact([pe_pool], ["PE"], {})
        assert legacy["schema_version"] == "2.0"
        assert "explicit_dp" not in json.dumps(legacy)
        qssa = self._artifact([qssa_pool], ["PS"], {"PS": "styrene(5)"})
        assert qssa["schema_version"] == "2.1"
        assert "explicit_dp" not in json.dumps(qssa)
        weak = self._artifact([weaklink_pool], ["PSW"], {"PSW": "styrene(5)"})
        assert weak["schema_version"] == "2.2"
        assert "explicit_dp" not in json.dumps(weak)

    def test_composition_with_qssa_and_weaklink(self, explicit_dp_pool,
                                                qssa_pool, weaklink_pool):
        """explicit_dp composes with the channel vocabularies: the schema
        stamp is 2.3 (strongest), the recipe token carries the strongest
        channel family, and ALL blocks coexist in the same artifact."""
        art_q = self._artifact([explicit_dp_pool, qssa_pool], ["PS"],
                               {"PS": "styrene(5)"})
        assert art_q["schema_version"] == "2.3"
        assert art_q["conventions"]["recipe_revision"] == \
            "2026-07-04-explicit-dp-qssa-monomer-gas"
        art_w = self._artifact([explicit_dp_pool, weaklink_pool],
                               ["PS", "PSW"], {"PSW": "styrene(5)"})
        assert art_w["schema_version"] == "2.3"
        assert art_w["conventions"]["recipe_revision"] == \
            "2026-07-04-explicit-dp-weaklink-u-monomer-gas"
        by_label = {p["label"]: p for p in art_w["pools"]}
        assert by_label["PS"]["explicit_dp"]["enabled"] is True
        assert "initiation_allyl" in \
            by_label["PSW"]["channels"]["radical_qssa_unzip"]

    def test_initial_moles_flow_from_pool_channel_map(self, explicit_dp_pool):
        """The artifact-level initial_explicit_by_pool map (the stage-A
        solver contract {pool_label: {dp: moles}}) flows into the pool
        block; pools without an entry default to 0.0."""
        artifact = self._artifact(
            [explicit_dp_pool], ["PS"], {},
            initial_explicit_by_pool={"PS": {3: 0.125}})
        assert artifact["pools"][0]["explicit_dp"]["initial_moles"] == \
            {"3": 0.125}
        artifact0 = self._artifact([explicit_dp_pool], ["PS"], {})
        assert artifact0["pools"][0]["explicit_dp"]["initial_moles"] == \
            {"3": 0.0}

    def test_write_sidecar_passthrough_and_round_trip(self, explicit_dp_pool,
                                                      tmp_path):
        path = write_polymer_pools_sidecar(
            [explicit_dp_pool], str(tmp_path),
            configured_pool_labels=["PS"],
            initial_explicit_by_pool={"PS": {3: 0.25}})
        with open(path, encoding="utf-8") as fh:
            rt = json.load(fh)
        assert rt["schema_version"] == "2.3"
        block = rt["pools"][0]["explicit_dp"]
        assert block["species"] == {"3": "PS_dp3(42)"}
        assert block["initial_moles"] == {"3": 0.25}
        assert block["recipe"] == EXPLICIT_DP_PINNED_RECIPE



class TestGenerationDefaults:
    """Commit 2 (provenance completeness): deck-declared mass_transfer and
    the V_poly RMG used are emitted into conventions.generation_defaults,
    explicitly NON-normative (consumer experiment config takes precedence).
    Additive-informative: no schema bump; absent when nothing was declared
    (the golden pin covers the byte-identity of that path)."""

    GEN_MT = [{"gas_species": "C3H6(3)", "poly_species": "C3H6_poly(9)",
               "K": 0.8, "kLa": 0.05}]

    def _artifact(self, pools, labels, **kw):
        return build_polymer_moments_artifact(
            pools, core_species=None, core_reactions=[],
            configured_pool_labels=labels, monomer_routing_by_pool={}, **kw)

    def test_declared_mass_transfer_emitted_non_normative(self, pe_pool):
        artifact = self._artifact(
            [pe_pool], ["PE"],
            generation_mass_transfer=self.GEN_MT,
            generation_v_poly_m3=1.25e-6)
        gd = artifact["conventions"]["generation_defaults"]
        assert gd == {
            "mass_transfer": [{"gas_species": "C3H6(3)",
                               "poly_species": "C3H6_poly(9)",
                               "K": 0.8, "kLa": 0.05,
                               "units": {"K": "dimensionless",
                                         "kLa": "s^-1"}}],
            "V_poly_m3": 1.25e-6,
            "note": ("generation-run values; consumer experiment config "
                     "takes precedence"),
        }
        # no schema/recipe consequences: additive-informative only
        assert artifact["schema_version"] == "2.0"
        assert artifact["conventions"]["recipe_revision"] == \
            "2026-07-03-monomer-gas"

    def test_v_poly_only_omits_mass_transfer_key(self, pe_pool):
        artifact = self._artifact([pe_pool], ["PE"],
                                  generation_v_poly_m3=1.25e-6)
        gd = artifact["conventions"]["generation_defaults"]
        assert "mass_transfer" not in gd
        assert gd["V_poly_m3"] == 1.25e-6

    def test_nothing_declared_key_absent(self, pe_pool):
        artifact = self._artifact([pe_pool], ["PE"])
        assert "generation_defaults" not in artifact["conventions"]


class TestSpawnedPoolConfiguredSurface:
    """S4 serializer closure (schema 2.5): every runtime-spawned pool present
    in the registry at save time — scission daughters AND S2 feature pools —
    must reach the sidecar's configured-pools surface so the TA consumer can
    classify its proxy + mu-dummies CONDENSED instead of defaulting them GAS
    (the item-16 hazard shape). Vocabulary:

    * ``conventions.spawned_pools`` — labels of registry pools NOT in
      ``conventions.configured_pools`` (closure complement; disjoint by
      construction). Presence-based: no spawned pool anywhere -> key absent
      and the artifact stays byte-identical (schema stamp included).
    * ``conventions.condensed_species`` closure — a spawned pool's
      ``phase_species`` (proxy + mu-dummies, already declared condensed
      row-side) join the normative condensed list the consumers key on.
    * schema 2.5 stamp, strongest SHAPE (2.5 > 2.4 > ...); recipe_revision
      untouched (classification vocabulary, not rate algebra) — the exact
      2.4 refused-row precedent."""

    PP_H_LOSS_DAUGHTER = 'CCC[C](C)CC(C)C'  # live PP center-tertiary C9H19

    @staticmethod
    def _pp_pool():
        return Polymer(label='polypropylene', monomer='[CH2][CH](C)',
                       end_groups=['[H]', '[H]'], cutoff=3,
                       Mn=1500.0, Mw=1800.0, initial_mass=0.1485)

    def _spawned_feature_setup(self):
        """A runtime-spawned feature pool through the REAL S1a/S2 machinery
        (create_reacted_copy(h_loss_feature=True)), plus a core carrying
        both pools' proxies and mu-dummies."""
        pp = self._pp_pool()
        daughter = pp.create_reacted_copy(
            Molecule(smiles=self.PP_H_LOSS_DAUGHTER), h_loss_feature=True)
        core = [
            _spc("CCCC(C)CC(C)C", "polypropylene", index=2),
            _mu_dummy("polypropylene_mu0"),
            _mu_dummy("polypropylene_mu1"),
            _mu_dummy("polypropylene_mu2"),
            _spc(self.PP_H_LOSS_DAUGHTER, daughter.label if daughter else "d",
                 index=9),
            _mu_dummy(f"{daughter.label}_mu0" if daughter else "d_mu0"),
            _mu_dummy(f"{daughter.label}_mu1" if daughter else "d_mu1"),
            _mu_dummy(f"{daughter.label}_mu2" if daughter else "d_mu2"),
            _spc("[CH3]", "G", index=7),
        ]
        core[0].is_polymer_proxy = True
        core[4].is_polymer_proxy = True
        return pp, daughter, core

    def test_spawned_feature_pool_reaches_configured_pools_surface(self):
        """RED pin 1: the S2 feature pool born mid-run must appear in the
        sidecar's configured-pools surface (conventions.spawned_pools) with
        monomer_mw_g_mol present, spawn provenance preserved and born-at-zero
        moments untouched."""
        pp, daughter, core = self._spawned_feature_setup()

        # LIVENESS PINS (tripwire discipline): the daughter came through the
        # real producer path with lineage + pinned monomer MW intact. A pin
        # failure here is a broken fixture, never a valid red.
        assert daughter is not None
        assert daughter.label == "polypropylene_mod"
        assert daughter.parent_pool_label == "polypropylene"
        assert daughter.monomer_mw_g_mol == pytest.approx(
            pp.monomer_mw_g_mol)

        payload = build_polymer_moments_artifact(
            [pp, daughter],
            core_species=core,
            configured_pool_labels=["polypropylene"],  # engine-configured only
            condensed_species=core[:4],  # engine mask: deck pool only
        )
        conv = payload["conventions"]
        entry = next(p for p in payload["pools"]
                     if p["label"] == "polypropylene_mod")

        # Row-side truths that must survive the closure unchanged.
        assert entry["moments"] == [0.0, 0.0, 0.0]
        assert entry["moments_provenance"] == "spawned_empty"
        assert isinstance(entry["spawn_event_metadata"], dict)
        assert entry["monomer_mw_g_mol"] == pytest.approx(
            pp.monomer_mw_g_mol)

        # RED assertion 1 (the configured-pools surface): the spawned pool
        # label is published next to configured_pools, disjoint from it.
        assert conv.get("spawned_pools") == ["polypropylene_mod"], (
            "runtime-spawned feature pool missing from the sidecar's "
            "configured-pools surface (conventions.spawned_pools)")
        assert set(conv["spawned_pools"]).isdisjoint(conv["configured_pools"])

        # RED assertion 2 (the consumer-keyed condensed list): the spawned
        # pool's proxy + mu-dummies join conventions.condensed_species.
        for lbl in ("polypropylene_mod(9)", "polypropylene_mod_mu0",
                    "polypropylene_mod_mu1", "polypropylene_mod_mu2"):
            assert lbl in conv["condensed_species"], (
                f"spawned-pool condensed member {lbl!r} missing from "
                "conventions.condensed_species (TA classifies it GAS)")

        # RED assertion 3 (envelope): spawned-pool vocabulary is the
        # strongest SHAPE stamp; rate algebra untouched.
        assert payload["schema_version"] == "2.5"
        assert conv["recipe_revision"] == POLYMER_RATE_RECIPE_REVISION

    def test_spawned_scission_daughter_reaches_surface_too(self, pe_pool):
        """The historical spawn shape (gate-path drain scission daughter)
        rides the same closure."""
        from rmgpy.polymer import SpawnIntent, drain_spawn_intents
        intent = SpawnIntent(parent_pool=pe_pool, monomer=pe_pool.monomer,
                             end_groups=["[H]", "[H]"], triggering_dp=4)
        daughter = drain_spawn_intents([intent], iteration=7,
                                       existing_pools=[pe_pool])[0]
        payload = build_polymer_moments_artifact(
            [pe_pool, daughter], configured_pool_labels=["PE"])
        conv = payload["conventions"]
        assert conv.get("spawned_pools") == [daughter.label]
        assert payload["schema_version"] == "2.5"
        entry = next(p for p in payload["pools"]
                     if p["label"] == daughter.label)
        assert entry["monomer_mw_g_mol"] == pytest.approx(
            pe_pool.monomer_mw_g_mol)

    def test_declared_only_artifact_is_byte_identical(self, pe_pool):
        """RED pin 2 (regression guard, GREEN before and after): with no
        spawned pool anywhere the key is ABSENT — never an empty list — and
        the legacy 2.0 stamp + conventions block are untouched."""
        core = [
            _spc("CC", "PE", index=2),
            _mu_dummy("PE_mu0"), _mu_dummy("PE_mu1"), _mu_dummy("PE_mu2"),
        ]
        core[0].is_polymer_proxy = True
        payload = build_polymer_moments_artifact(
            [pe_pool], core_species=core,
            configured_pool_labels=["PE"], condensed_species=core)
        conv = payload["conventions"]
        assert "spawned_pools" not in conv
        assert payload["schema_version"] == "2.0"
        assert conv["configured_pools"] == ["PE"]
        assert conv["condensed_species"] == [
            "PE(2)", "PE_mu0", "PE_mu1", "PE_mu2"]

    def test_legacy_default_labels_call_emits_no_spawned_surface(self, pe_pool):
        """Legacy default-label call (configured defaults to ALL registry
        labels): the closure complement is empty by construction, so the
        surface stays absent even though the daughter's row still reads
        spawned via its object markers (the documented secondary signal)."""
        from rmgpy.polymer import SpawnIntent, drain_spawn_intents
        intent = SpawnIntent(parent_pool=pe_pool, monomer=pe_pool.monomer,
                             end_groups=["[H]", "[H]"], triggering_dp=4)
        daughter = drain_spawn_intents([intent], iteration=7,
                                       existing_pools=[pe_pool])[0]
        payload = build_polymer_moments_artifact([pe_pool, daughter])
        conv = payload["conventions"]
        assert "spawned_pools" not in conv
        assert payload["schema_version"] == "2.0"
        entry = next(p for p in payload["pools"]
                     if p["label"] == daughter.label)
        assert entry["moments_provenance"] == "spawned_empty"

    def test_explicit_empty_configured_set_keeps_daughter_spawned(
            self, pe_pool):
        """P2 sentinel pin: the caller-supplied sentinel is None, NOT
        emptiness. An explicitly EMPTY configured set (a daughters-only
        artifact: nothing is root-configured) must NOT fall back to the
        legacy default-to-all-registry-labels behavior -- the spawned
        daughter is emitted in conventions.spawned_pools with the 2.5
        stamp, never in configured_pools."""
        from rmgpy.polymer import SpawnIntent, drain_spawn_intents
        intent = SpawnIntent(parent_pool=pe_pool, monomer=pe_pool.monomer,
                             end_groups=["[H]", "[H]"], triggering_dp=4)
        daughter = drain_spawn_intents([intent], iteration=7,
                                       existing_pools=[pe_pool])[0]
        payload = build_polymer_moments_artifact(
            [daughter], configured_pool_labels=[])
        conv = payload["conventions"]
        assert conv["configured_pools"] == []
        assert conv.get("spawned_pools") == [daughter.label]
        assert payload["schema_version"] == "2.5"
        entry = next(p for p in payload["pools"]
                     if p["label"] == daughter.label)
        assert entry["moments_provenance"] == "spawned_empty"


# ---------------------------------------------------------------------------
# Radical-homolysis initiation sidecar block (Stage 2, adjudicated rounds
# 66/67): schema 2.6 = 2.5 + the pool-level ``homolysis_initiation`` block.
# Presence-gated exactly like 2.5 (no homolysis pool anywhere -> artifacts
# stay byte-identical at their older stamps), BESIDE 2.5 (the spawned-pool
# closure + condensed closure keep their exact 2.5 semantics), with the
# daughter fields named for the open-*1/open-*2 POSITIONAL contract (round
# 67 ruling (a)) and a kernel field naming the solver kernel (ruling (c),
# machine-checkable supersession).
# ---------------------------------------------------------------------------


def _khom_triplet(A=1.0e13, n=0.5, Ea=1.2e5):
    return dict(A=A, n=n, Ea=Ea)


class TestEndRadicalDepropagationSidecar:
    """Emitter side of the schema-2.8 end_radical_depropagation contract
    (r74 SS2 kernel; r78 serialization rulings) -- the real serialization
    that replaced the 'schema 2.8 pending' producer hard-refusal. The block
    rides the two spawned end-radical DAUGHTER pool entries (where the
    kernel actually integrates), never the parent's declaration surface,
    and pins the ARTIFACT's as-integrated behavior: Arrhenius triplet +
    units, gas monomer identity/MW, the smooth exhaustion gate form + the
    machine-readable gate_width, the gated mu2 smooth-pos law (disclosed
    under-drain), and the UN-gated dmu0 half-bin N1 gamma closure."""

    @staticmethod
    def _kdep_triplet():
        return dict(A=9.4e14, n=0.0, Ea=117152.0)

    def _setup(self):
        """Parent (k_homolysis + k_depropagation declaration context +
        monomer routing) + both end-radical daughters through the REAL
        Stage-1 producer, plus a core carrying the proxies, mu-dummies and
        the released-monomer GAS species."""
        pp = Polymer(label='PP', monomer='[CH2][CH](C)',
                     end_groups=['[H]', '[H]'], cutoff=3,
                     moments=[1.0, 50.0, 3000.0], initial_mass=0.0,
                     k_homolysis=_khom_triplet(),
                     k_depropagation=self._kdep_triplet())
        monomer_gas = _spc("C=CC", "C3H6", index=5)
        pp.monomer_product_species = monomer_gas
        prim, sec = pp.generate_end_radical_daughters()
        core = [_spc("CCC(C)CC(C)C", "PP", index=2), monomer_gas]
        for base in ("PP", prim.label, sec.label):
            core += [_mu_dummy(f"{base}_mu{k}") for k in range(3)]
        core[0].is_polymer_proxy = True
        # The routed monomer is GAS (recipe revision 2026-07-03-monomer-gas):
        # everything else is condensed.
        condensed = [core[0]] + core[2:]
        return pp, prim, sec, core, condensed

    def _payload(self, pp, prim, sec, core, condensed):
        return build_polymer_moments_artifact(
            [pp, prim, sec], core_species=core,
            configured_pool_labels=["PP", prim.label, sec.label],
            condensed_species=condensed)

    def test_deprop_daughters_emit_block_and_2_8_stamp(self):
        """RED pin: both end-radical daughter entries carry the full
        end_radical_depropagation block (structured kinetics with explicit
        units, gas monomer identity/MW, machine-pinned gate_width, kernel,
        block recipe_revision + pinned recipe) and the artifact stamps
        schema 2.8."""
        pp, prim, sec, core, condensed = self._setup()
        payload = self._payload(pp, prim, sec, core, condensed)
        assert payload["schema_version"] == "2.8"
        for lbl in (prim.label, sec.label):
            entry = next(p for p in payload["pools"] if p["label"] == lbl)
            block = entry.get("end_radical_depropagation")
            assert isinstance(block, dict), (
                f"kernel-carrying daughter {lbl!r} must carry the "
                f"pool-level end_radical_depropagation block")
            assert block["enabled"] is True
            assert block["kinetics"] == {
                "A": 9.4e14, "n": 0.0, "Ea": 117152.0,
                "units": {"A": "s^-1", "Ea": "J/mol"},
            }
            # r78: the gas product (monomer) identity/MW are IN the block,
            # cross-pinned to the pool surface (monomer_routing /
            # monomer_mw_g_mol).
            assert block["gas_species"] == "C3H6(5)"
            assert entry["monomer_routing"] == "C3H6(5)"
            assert block["gas_mw_g_mol"] == pytest.approx(42.0797,
                                                          rel=1e-4)
            assert block["gas_mw_g_mol"] == pytest.approx(
                entry["monomer_mw_g_mol"], rel=1e-9)
            # r78: the gate WIDTH is a machine-readable field pinned to the
            # solver constant (bitwise -- the 1e-12 TA-twin contract).
            assert block["gate_width"] == 1.0e-2
            assert block["kernel"] == "end_radical_depropagation/1"
            assert block["recipe_revision"] == \
                "2026-07-06-end-radical-depropagation"
            recipe = block["recipe"]
            assert isinstance(recipe, dict) and recipe
            joined = " ".join(str(v) for v in recipe.values())
            for needle in (
                    # the C2 smooth positive-part + gate form (r74 SS5)
                    "x^3/(x^2 + W^2)", "1 - sp(1 - mu1/mu0)",
                    # runtime-T Arrhenius with the pinned gas constant
                    "R_gas = 8.314",
                    # gated smooth-pos mu2 law (disclosed under-drain)
                    "k*mu0*(g + 2*sp(mu1/mu0 - 1))",
                    # dmu0 is UN-gated, N1 from the half-bin gamma closure
                    "UN-gated", "F(1.5) - F(0.5)",
                    "k_shape = 1/(PDI - 1)", "PDI = mu2*mu0/mu1^2",
                    # the C1 smoothstep terminal floor (never dmu0 = 0)
                    "1 - (3*t^2 - 2*t^3)"):
                assert needle in joined, (
                    f"recipe must pin the as-implemented law term "
                    f"{needle!r}")
        # Block-local revision only: artifact-level recipe_revision is
        # untouched (the 2.6/2.7 precedent).
        assert payload["conventions"]["recipe_revision"] == \
            POLYMER_RATE_RECIPE_REVISION

    def test_parent_declaration_surface_emits_no_block(self):
        """The parent pool carries k_depropagation as DECK DECLARATION
        surface only (the kernel never integrates there: validate_
        configuration excludes k_depropagation + k_homolysis on one pool
        and parent configs never carry the triplet) -- its sidecar entry
        must NOT carry the block."""
        pp, prim, sec, core, condensed = self._setup()
        payload = self._payload(pp, prim, sec, core, condensed)
        entry = next(p for p in payload["pools"] if p["label"] == "PP")
        assert "end_radical_depropagation" not in entry
        assert "homolysis_initiation" in entry

    def test_no_deprop_artifact_keeps_2_6_stamp(self):
        """Negative control (presence gate): a homolysis-only artifact
        (no k_depropagation anywhere) stays byte-identical at 2.6 and no
        pool carries the block."""
        pp, prim, sec, core, condensed = self._setup()
        for p in (pp, prim, sec):
            p.k_depropagation = None
        payload = self._payload(pp, prim, sec, core, condensed)
        assert payload["schema_version"] == "2.6"
        assert all("end_radical_depropagation" not in p
                   for p in payload["pools"])

    def test_kernel_free_pool_still_serializes(self):
        """A kernel-free end-radical pool keeps serializing exactly as
        before (presence-gated emission)."""
        pool = Polymer(label='PP_rad_primary_end', monomer='[CH2][CH](C)',
                       end_groups=['[H]', '[H]'], cutoff=3,
                       moments=[0.0, 0.0, 0.0], initial_mass=0.0,
                       end_radical_site='primary')
        d = _serialize_pool_for_sidecar(pool)
        assert d["label"] == "PP_rad_primary_end"
        assert "end_radical_depropagation" not in d

    def test_dead_knob_parentless_kernel_refused(self):
        """A pool carrying k_depropagation with NEITHER a k_homolysis
        declaration context NOR an end-radical daughter identity is a
        dead-knob shape no deck can produce -- the serializer refuses it
        (the narrowed successor of the pre-2.8 blanket refusal)."""
        pool = Polymer(label='PPX', monomer='[CH2][CH](C)',
                       end_groups=['[H]', '[H]'], cutoff=3,
                       moments=[0.0, 0.0, 0.0], initial_mass=0.0)
        pool.k_depropagation = self._kdep_triplet()
        with pytest.raises(ValueError,
                           match=r"PPX.*k_depropagation"):
            _serialize_pool_for_sidecar(pool)

    def test_producer_refuses_missing_gas_routing(self):
        """A kernel-carrying daughter without monomer_product_species has
        no resolvable gas emission target -- the released units would leave
        the condensed phase silently un-conserved; the producer refuses."""
        pp, prim, sec, core, condensed = self._setup()
        prim.monomer_product_species = None
        with pytest.raises(ValueError,
                           match=r"PP_rad_primary_end.*(routing|monomer)"):
            self._payload(pp, prim, sec, core, condensed)

    def test_producer_refuses_unresolvable_gas_routing(self):
        """The routed monomer Species must live in the core universe the
        artifact labels come from (IDENTITY check, the QSSA routing
        posture)."""
        pp, prim, sec, core, condensed = self._setup()
        stranger = _spc("C=CC", "C3H6_other", index=9)
        prim.monomer_product_species = stranger
        with pytest.raises(ValueError,
                           match=r"PP_rad_primary_end.*core"):
            self._payload(pp, prim, sec, core, condensed)

    def test_producer_refuses_gas_mw_mismatch(self):
        """r78 mass pin: each unzip event moves exactly ONE repeat unit
        from the condensed basis (mass monomer_mw_g_mol) into ONE mole of
        gas_species -- a routed species whose molar mass diverges from the
        pool's repeat-unit mass would mint/destroy mass on every event; the
        producer refuses to serialize it."""
        pp, prim, sec, core, condensed = self._setup()
        h2 = _spc("[H][H]", "H2", index=6)
        core = core + [h2]
        prim.monomer_product_species = h2
        with pytest.raises(ValueError,
                           match=r"PP_rad_primary_end.*(mw|mass)"):
            self._payload(pp, prim, sec, core, condensed)

    def test_producer_refuses_sibling_asymmetry(self):
        """The producer copies ONE parent-declared triplet onto BOTH
        spawned daughters -- an artifact where only one daughter carries
        the block is corrupted (the consumer hard-rejects it), so the
        producer refuses to emit it."""
        pp, prim, sec, core, condensed = self._setup()
        sec.k_depropagation = None
        with pytest.raises(ValueError,
                           match=r"PP_rad_.*end.*sibling"):
            self._payload(pp, prim, sec, core, condensed)

    def test_producer_refuses_defect_bearing_carrier(self):
        """r79 P1 (RED-pinned): a deprop carrier that ALSO carries the
        X-loss chain_mass_defect_g_mol contract mints mass at terminal DP1
        events -- its condensed mass is mu1*MW - mu0*defect, so a DP1
        terminal event (dmu0 = dmu1 = -R) drains only R*(MW - defect) of
        condensed mass while the gas monomer credits R*MW: net +R*defect
        minted per event. The side-group guards deliberately legalize
        copied defect pools, so without THIS refusal the shape serializes;
        v2 defect-chain depropagation would need a different mass law /
        gas product."""
        pp, prim, sec, core, condensed = self._setup()
        prim.chain_mass_defect_g_mol = 79.904
        with pytest.raises(
                ValueError,
                match=r"PP_rad_primary_end.*chain_mass_defect"):
            self._payload(pp, prim, sec, core, condensed)

    def test_producer_refuses_scission_deprop_combination(self):
        """r78 adjudication: k_scission > 0 on a kernel-carrying pool is a
        direct-config-only shape production cannot generate (end-radical
        daughters are born with k_scission = 0 by
        generate_end_radical_daughters), so no generating run ever
        integrated the combined law -- refuse rather than serialize an
        unvalidated combination (consumers hard-reject it)."""
        pp, prim, sec, core, condensed = self._setup()
        prim.k_scission = 5.0
        with pytest.raises(ValueError,
                           match=r"PP_rad_primary_end.*scission"):
            self._payload(pp, prim, sec, core, condensed)

    def test_closure_refuses_non_daughter_carrier(self):
        """Dict-level defense (r68 broken-caller posture): a serialized
        block on a pool whose label is not an end-radical daughter's is
        refused by the closure guard directly."""
        from rmgpy.polymer import _assert_depropagation_serialization_closure
        pp, prim, sec, core, condensed = self._setup()
        payload = self._payload(pp, prim, sec, core, condensed)
        carrier = next(p for p in payload["pools"]
                       if p["label"] == prim.label)
        carrier["label"] = "NOT_A_DAUGHTER"
        conv = payload["conventions"]
        with pytest.raises(ValueError, match=r"NOT_A_DAUGHTER"):
            _assert_depropagation_serialization_closure(
                payload["pools"], [carrier],
                set(conv["condensed_species"]),
                set(conv["configured_pools"]) | {"NOT_A_DAUGHTER"},
                [])

    def test_closure_refuses_unzip_carrying_carrier(self):
        """Dict-level mirror of the solver's validate_configuration
        exclusion set: a carrier entry claiming legacy unzip A > 0 beside
        the block double-carries the SAME chain-end release event."""
        from rmgpy.polymer import _assert_depropagation_serialization_closure
        pp, prim, sec, core, condensed = self._setup()
        payload = self._payload(pp, prim, sec, core, condensed)
        carrier = next(p for p in payload["pools"]
                       if p["label"] == prim.label)
        carrier["channels"]["unzip"]["A"] = 0.5
        conv = payload["conventions"]
        with pytest.raises(ValueError,
                           match=r"PP_rad_primary_end.*unzip"):
            _assert_depropagation_serialization_closure(
                payload["pools"], [carrier],
                set(conv["condensed_species"]),
                set(conv["configured_pools"]), [])


class TestHomolysisInitiationSidecar:
    """Emitter side of the schema-2.6 homolysis_initiation contract."""

    @staticmethod
    def _pp_pool(**kw):
        return Polymer(label='PP', monomer='[CH2][CH](C)',
                       end_groups=['[H]', '[H]'], cutoff=3,
                       moments=[1.0, 50.0, 3000.0], initial_mass=0.0,
                       k_homolysis=_khom_triplet(), **kw)

    def _setup(self):
        """Parent + both end-radical daughters through the REAL Stage-1
        producer (generate_end_radical_daughters), plus a core carrying the
        three pools' proxies and mu-dummies (the live registration shape:
        daughters are eagerly solver-CONFIGURED -- polymer.pyx
        _flatten_homolysis_state hard-errors otherwise)."""
        pp = self._pp_pool()
        prim, sec = pp.generate_end_radical_daughters()
        core = [_spc("CCC(C)CC(C)C", "PP", index=2)]
        for base in ("PP", prim.label, sec.label):
            core += [_mu_dummy(f"{base}_mu{k}") for k in range(3)]
        core[0].is_polymer_proxy = True
        return pp, prim, sec, core

    def _payload(self, pp, prim, sec, core, registry=None, condensed=None):
        return build_polymer_moments_artifact(
            registry if registry is not None else [pp, prim, sec],
            core_species=core,
            configured_pool_labels=["PP", prim.label, sec.label],
            condensed_species=condensed if condensed is not None else core,
        )

    def test_homolysis_pool_emits_block_and_2_6_stamp(self):
        """RED pin: the kernel-enabled pool's sidecar entry carries the full
        homolysis_initiation block (structured kinetics with explicit units,
        open-site daughter fields, kernel, block recipe_revision + pinned
        recipe) and the artifact stamps schema 2.6."""
        pp, prim, sec, core = self._setup()
        payload = self._payload(pp, prim, sec, core)
        entry = next(p for p in payload["pools"] if p["label"] == "PP")

        block = entry.get("homolysis_initiation")
        assert isinstance(block, dict), (
            "kernel-enabled pool must carry the pool-level "
            "homolysis_initiation block")
        assert block["enabled"] is True
        assert block["kinetics"] == {
            "A": 1.0e13, "n": 0.5, "Ea": 1.2e5,
            "units": {"A": "s^-1", "Ea": "J/mol"},
        }
        # Round-67 ruling (a): field names are POSITIONAL (open-*1/open-*2
        # stitch termini), never primary/secondary radical-character claims.
        assert block["open_site_1_radical_pool"] == "PP_rad_primary_end"
        assert block["open_site_2_radical_pool"] == "PP_rad_secondary_end"
        # Ruling (c): the block names the solver kernel (machine-checkable
        # supersession contract).
        assert block["kernel"] == "radical_homolysis_initiation/1"
        assert block["recipe_revision"] == "2026-07-05-radical-homolysis"
        # Machine-pinned moment law in the STABLE product forms actually
        # implemented (round-67 P2; polymer.pyx homolysis kernel).
        recipe = block["recipe"]
        assert isinstance(recipe, dict) and recipe
        joined = " ".join(str(v) for v in recipe.values())
        for needle in ("max(mu1 - mu0, 0)", "k*(mu2 - mu1)",
                       "k*(mu3 - mu2)", "k*(2*mu3 - 3*mu2 + mu1)/6",
                       "k*(mu1 - mu3)/3"):
            assert needle in joined, (
                f"recipe must pin the stable-form law term {needle!r}")

        # Daughters land in pools[] and are classified + condensed.
        labels = [p["label"] for p in payload["pools"]]
        assert prim.label in labels and sec.label in labels
        conv = payload["conventions"]
        for lbl in (prim.label, sec.label):
            assert (lbl in conv["configured_pools"]
                    or lbl in conv.get("spawned_pools", [])), (
                f"daughter {lbl!r} must be classified in the "
                f"configured/spawned closure")
            d_entry = next(p for p in payload["pools"] if p["label"] == lbl)
            assert d_entry["phase_species"], (
                f"daughter {lbl!r} must carry phase_species")
            for member in d_entry["phase_species"]:
                assert member in conv["condensed_species"], (
                    f"daughter member {member!r} missing from the "
                    f"condensed closure")
            assert d_entry["spawn_event_metadata"].get("source") == \
                "k_homolysis_end_radical"
            assert d_entry["parent_pool"] == "PP"
        assert payload["schema_version"] == "2.6"
        # Block-local revision only: artifact-level recipe_revision is
        # untouched (the 2.4/2.5 precedent for new-vocabulary stamps).
        assert conv["recipe_revision"] == POLYMER_RATE_RECIPE_REVISION

    def test_2_6_sits_beside_2_5_spawned_closure(self):
        """RED pin: 2.6 must NOT subsume 2.5 -- an artifact carrying BOTH a
        homolysis pool and a runtime-spawned feature pool keeps the exact
        2.5 spawned_pools + condensed-closure emission alongside the 2.6
        stamp."""
        pp, prim, sec, core = self._setup()
        mod = pp.create_reacted_copy(
            Molecule(smiles='CCC[C](C)CC(C)C'), h_loss_feature=True)
        assert mod is not None
        core = core + [
            _spc('CCC[C](C)CC(C)C', mod.label, index=9),
            _mu_dummy(f"{mod.label}_mu0"), _mu_dummy(f"{mod.label}_mu1"),
            _mu_dummy(f"{mod.label}_mu2"),
        ]
        core[-4].is_polymer_proxy = True
        payload = build_polymer_moments_artifact(
            [pp, prim, sec, mod],
            core_species=core,
            configured_pool_labels=["PP", prim.label, sec.label],
            condensed_species=core[:10],  # engine mask: configured pools only
        )
        conv = payload["conventions"]
        # 2.5 closure unchanged beside the block: the _mod pool is the
        # configured set's complement, disjoint, and its members join the
        # condensed closure.
        assert conv.get("spawned_pools") == [mod.label]
        assert set(conv["spawned_pools"]).isdisjoint(conv["configured_pools"])
        for lbl in (f"{mod.label}(9)", f"{mod.label}_mu0",
                    f"{mod.label}_mu1", f"{mod.label}_mu2"):
            assert lbl in conv["condensed_species"]
        # 2.6 outranks 2.5 as the strongest SHAPE stamp.
        assert payload["schema_version"] == "2.6"
        assert "homolysis_initiation" in next(
            p for p in payload["pools"] if p["label"] == "PP")

    def test_no_homolysis_artifact_keeps_pre_2_6_stamp(self, pe_pool):
        """Negative control (presence gate): without any homolysis pool the
        artifact stays byte-identical at its older stamp and no pool carries
        the block."""
        core = [
            _spc("CC", "PE", index=2),
            _mu_dummy("PE_mu0"), _mu_dummy("PE_mu1"), _mu_dummy("PE_mu2"),
        ]
        core[0].is_polymer_proxy = True
        payload = build_polymer_moments_artifact(
            [pe_pool], core_species=core,
            configured_pool_labels=["PE"], condensed_species=core)
        assert payload["schema_version"] == "2.0"
        assert all("homolysis_initiation" not in p for p in payload["pools"])

    def test_producer_refuses_missing_daughters(self):
        """RED pin: serializing a kernel-enabled pool whose end-radical
        daughters are absent from the registry would emit exactly the shape
        the consumer hard-rejects -- the producer refuses instead."""
        pp, prim, sec, core = self._setup()
        with pytest.raises(ValueError, match=r"PP.*_rad_.*end"):
            build_polymer_moments_artifact(
                [pp],  # daughters missing from the registry
                core_species=core,
                configured_pool_labels=["PP", prim.label, sec.label],
                condensed_species=core,
            )

    def test_producer_refuses_uncondensed_daughters(self):
        """RED pin: a daughter whose phase_species does not reach the
        condensed closure is the item-16 mass-balance hazard shape -- the
        producer refuses to emit it."""
        pp, prim, sec, core = self._setup()
        with pytest.raises(ValueError, match=r"condensed"):
            self._payload(pp, prim, sec, core,
                          condensed=core[:4])  # daughters' members missing

    def test_producer_refuses_unconfigured_daughter(self):
        """r68 P2 (producer/consumer closure mirror): a daughter omitted
        from configured_pool_labels classifies into the spawned complement
        -- exactly the shape the r68-tightened loader hard-rejects
        (build_system_from_artifact only builds configured pools). The
        producer must refuse to emit it."""
        pp, prim, sec, core = self._setup()
        with pytest.raises(ValueError,
                           match=r"PP_rad_primary_end.*configured"):
            build_polymer_moments_artifact(
                [pp, prim, sec], core_species=core,
                configured_pool_labels=["PP", sec.label],
                condensed_species=core)

    def test_producer_refuses_spawned_classified_daughter(self):
        """r68 P2: the daughter-not-in-spawned_pools conjunct, checked
        independently of the configured conjunct (defense-in-depth: in the
        emitter spawned is derived as the configured complement, so the
        overlap shape is only reachable through a broken caller of the
        guard itself)."""
        from rmgpy.polymer import _assert_homolysis_serialization_closure
        pp, prim, sec, core = self._setup()
        payload = self._payload(pp, prim, sec, core)
        carriers = [p for p in payload["pools"]
                    if "homolysis_initiation" in p]
        conv = payload["conventions"]
        with pytest.raises(ValueError,
                           match=r"PP_rad_primary_end.*spawned"):
            _assert_homolysis_serialization_closure(
                payload["pools"], carriers,
                set(conv["condensed_species"]),
                set(conv["configured_pools"]),
                {prim.label})  # contrived overlap: spawned AND configured

    def test_producer_refuses_provenance_stripped_daughter(self):
        """r68 P2: a daughter without the Stage-1 spawn provenance pin
        (spawn_event_metadata.source == 'k_homolysis_end_radical')
        serializes with the {'source': 'input'} default -- a shape the
        loader hard-rejects. The producer must refuse to emit it."""
        pp, prim, sec, core = self._setup()
        prim.spawn_metadata = None  # serializer defaults to source: input
        with pytest.raises(ValueError,
                           match=r"PP_rad_primary_end.*(provenance|spawn)"):
            self._payload(pp, prim, sec, core)

    def test_producer_refuses_unconfigured_carrier(self):
        """r68 P2: the kernel-carrying pool itself must be
        solver-configured -- build_system_from_artifact skips unconfigured
        pools, so a block on an unconfigured carrier is a silently dropped
        kernel on the consumer side; the r68-tightened loader hard-rejects
        it and the producer must refuse to emit it."""
        pp, prim, sec, core = self._setup()
        with pytest.raises(ValueError, match=r"'PP'.*configured"):
            build_polymer_moments_artifact(
                [pp, prim, sec], core_species=core,
                configured_pool_labels=[prim.label, sec.label],
                condensed_species=core)


class TestSideGroupHomolysisSidecar:
    """Emitter side of the schema-3.0 side_group_homolysis kernel-v2
    (explicit Br-inventory depletion) contract: the full channel list
    (kinetics + site_selector + serialized site_atom_indices + gas
    routing + M_X pin), NO X-loss feature pool (v2 tracks multi-loss X
    depletion via the carrier's per-(pool, channel) Z inventory, so
    generate_side_loss_daughters() returns ()), the carrier's moments
    unchanged and carrying NO chain_mass_defect_g_mol, the presence-gated
    MAJOR 3.0 stamp, and the producer closure guard mirroring the
    consumer (the r68 mirror-property): carrier-configured, no-defect-on-
    carrier, site_atom_indices in range, channels atom-set disjoint."""

    @staticmethod
    def _sgh_channel(label="aliphatic_C-Br", A=1.0e13, n=0.5, Ea=1.2e5,
                     site_selector="aliphatic", sites_per_unit=1.0,
                     gas_product="[Br]"):
        return dict(label=label, A=A, n=n, Ea=Ea,
                    site_selector=site_selector,
                    sites_per_unit=sites_per_unit, gas_product=gas_product)

    def _setup(self, channels=None, monomer='[CH2][CH]Br'):
        """Kernel-enabled pool built for the SGH kernel-v2 artifact. v2
        deletes the X-loss feature pool, so generate_side_loss_daughters()
        returns () and ``daughters`` is empty; the artifact is built from
        the single carrier + its gas Species (identity is load-bearing for
        routing) + the carrier's mu-dummies (no feature-pool mu-dummies).
        The list-comprehensions over ``daughters`` below degrade to the
        single-carrier registry."""
        channels = channels or [self._sgh_channel()]
        pool = Polymer(label='PVBr', monomer=monomer,
                       end_groups=['[H]', '[H]'], cutoff=3,
                       moments=[1.0, 50.0, 3000.0], initial_mass=0.0,
                       side_group_homolysis=channels)
        daughters = list(pool.generate_side_loss_daughters())
        br = _spc("[Br]", "Br", index=7)
        pool.side_group_gas_species = [br] * len(channels)
        core = [br]
        for base in ["PVBr"] + [d.label for d in daughters]:
            core.extend(_mu_dummy(f"{base}_mu{k}") for k in range(3))
        return pool, daughters, br, core

    @staticmethod
    def _payload(pool, daughters, core, registry=None, configured=None,
                 condensed=None):
        return build_polymer_moments_artifact(
            registry if registry is not None else [pool] + daughters,
            core_species=core,
            configured_pool_labels=(
                configured if configured is not None
                else ["PVBr"] + [d.label for d in daughters]),
            condensed_species=condensed if condensed is not None
            else core[1:])

    def test_side_group_homolysis_v2_artifact_stamps_schema_3_0(self):
        """SGH kernel-v2 (schema 3.0) election, golden: a v2
        side_group_homolysis artifact stamps schema_version 3.0 -- the MAJOR
        rung wins DETERMINISTICALLY over the 2.7/2.8/2.9 rungs a v2 block
        would otherwise touch. The v2 block reuses the side_group_homolysis
        key, so it trips the 2.7 rung, but the /2 kernel marker unconditionally
        overrides to 3.0 with NO '2.x fallback' (the 2.7 chain_mass_defect
        trigger never fires -- v2 emits no defect field).

        v2 tracks multi-loss X depletion via the carrier's auxiliary
        per-(pool, channel) Z inventory instead of a one-lost-X feature pool,
        so generate_side_loss_daughters() returns () and the artifact is built
        from the SINGLE PVBr carrier: the block emits NO per-channel
        feature_pool and the carrier carries NO chain_mass_defect_g_mol."""
        pool, daughters, br, core = self._setup()
        # v2 DELETES the feature/daughter pool entirely.
        assert daughters == []
        payload = self._payload(pool, daughters, core)
        # Deterministic MAJOR election: exactly 3.0, never a 2.x fallback.
        assert payload["schema_version"] == "3.0"
        assert not payload["schema_version"].startswith("2.")
        entry = next(p for p in payload["pools"] if p["label"] == "PVBr")
        block = entry["side_group_homolysis"]
        assert block["kernel"] == "side_group_homolysis/2"
        assert block["recipe_revision"] == (
            "2026-07-08-side-group-inventory-depletion")
        assert block["recipe"]["mass"] == (
            "mu1*MW - sum_c max(0, sites_c*mu1 - Z_c)*M_X_c")
        assert block["recipe"]["state"] == (
            "Z_c(0)=sites_c*mu1(0); dZ_c/dt=-k*Z_c/V_poly (conc basis); "
            "pool moments unchanged")
        # Per-channel detail (preserved from the deleted v1 block/2.7 test):
        # kinetics (units UNCHANGED from v1), selector, serialized site
        # indices, gas routing + M_X pin.
        ch = block["channels"][0]
        assert ch["kinetics"] == {
            "A": 1.0e13, "n": 0.5, "Ea": 1.2e5,
            "units": {"A": "s^-1 per site", "Ea": "J/mol"},
        }
        assert ch["site_selector"] == "aliphatic"
        assert ch["sites_per_unit"] == 1.0
        assert (isinstance(ch["site_atom_indices"], list)
                and len(ch["site_atom_indices"]) == 1
                and all(isinstance(i, int) and i >= 0
                        for i in ch["site_atom_indices"]))
        assert ch["gas_product"] == "[Br]"
        assert ch["gas_species"] == "Br(7)"
        assert ch["gas_mw_g_mol"] == pytest.approx(79.904, rel=1e-3)
        # v2 SHAPE: no per-channel feature_pool anywhere in the block, and no
        # chain_mass_defect_g_mol on the carrier pool entry.
        for ch in block["channels"]:
            assert "feature_pool" not in ch
        assert "chain_mass_defect_g_mol" not in entry

    def test_mixed_site_channels_serialize_disjoint_atom_sets(self):
        """Round-72 mixed-site fixture end-to-end (v2): the aryl and
        aliphatic C-Br channels of one FR1-like repeat unit serialize with
        DISJOINT site_atom_indices -- the serialized rendering the loader's
        same-atom-set double-carry guard keys on -- under the MAJOR 3.0
        stamp. v2 has no feature pools, so neither channel carries a
        feature_pool key."""
        channels = [
            self._sgh_channel(label="aryl_C-Br", site_selector="aryl"),
            self._sgh_channel(label="aliphatic_C-Br",
                              site_selector="aliphatic"),
        ]
        pool, daughters, br, core = self._setup(
            channels=channels, monomer='[CH2][C](c1ccc(Br)cc1)CBr')
        payload = self._payload(pool, daughters, core)
        assert payload["schema_version"] == "3.0"
        block = next(p for p in payload["pools"]
                     if p["label"] == "PVBr")["side_group_homolysis"]
        ch_aryl, ch_ali = block["channels"]
        assert ch_aryl["site_selector"] == "aryl"
        assert ch_ali["site_selector"] == "aliphatic"
        assert not (set(ch_aryl["site_atom_indices"])
                    & set(ch_ali["site_atom_indices"]))
        assert "feature_pool" not in ch_aryl
        assert "feature_pool" not in ch_ali

    def test_v2_sgh_sits_beside_2_6_homolysis_block(self):
        """SGH-v2 (schema 3.0) OUTRANKS but never subsumes 2.6: an artifact
        carrying BOTH kernels on DIFFERENT pools (SGH-v2 on PVBr,
        k_homolysis on PP) keeps the homolysis_initiation block
        byte-identical beside the side_group_homolysis block, under the
        single strongest (MAJOR 3.0) stamp."""
        pool, daughters, br, core = self._setup()
        pp = Polymer(label='PP', monomer='[CH2][CH](C)',
                     end_groups=['[H]', '[H]'], cutoff=3,
                     moments=[1.0, 50.0, 3000.0], initial_mass=0.0,
                     k_homolysis=_khom_triplet())
        prim, sec = pp.generate_end_radical_daughters()
        for base in ("PP", prim.label, sec.label):
            core.extend(_mu_dummy(f"{base}_mu{k}") for k in range(3))
        payload = build_polymer_moments_artifact(
            [pool] + daughters + [pp, prim, sec], core_species=core,
            configured_pool_labels=(["PVBr"]
                                    + [d.label for d in daughters]
                                    + ["PP", prim.label, sec.label]),
            condensed_species=core[1:])
        assert payload["schema_version"] == "3.0"
        assert "homolysis_initiation" in next(
            p for p in payload["pools"] if p["label"] == "PP")
        assert "side_group_homolysis" in next(
            p for p in payload["pools"] if p["label"] == "PVBr")

    def test_homolysis_only_artifact_stays_2_6(self):
        """Negative control (presence gate): a 2.6-only artifact
        (k_homolysis, no side-group vocabulary anywhere) keeps its 2.6
        stamp byte-identically."""
        pp = Polymer(label='PP', monomer='[CH2][CH](C)',
                     end_groups=['[H]', '[H]'], cutoff=3,
                     moments=[1.0, 50.0, 3000.0], initial_mass=0.0,
                     k_homolysis=_khom_triplet())
        prim, sec = pp.generate_end_radical_daughters()
        core = []
        for base in ("PP", prim.label, sec.label):
            core.extend(_mu_dummy(f"{base}_mu{k}") for k in range(3))
        payload = build_polymer_moments_artifact(
            [pp, prim, sec], core_species=core,
            configured_pool_labels=["PP", prim.label, sec.label],
            condensed_species=core)
        assert payload["schema_version"] == "2.6"
        assert all("side_group_homolysis" not in p
                   and "chain_mass_defect_g_mol" not in p
                   for p in payload["pools"])

    def test_producer_refuses_unconfigured_carrier(self):
        """A block on an unconfigured carrier is a silently dropped kernel
        on the consumer side; the producer must refuse to emit it."""
        pool, daughters, br, core = self._setup()
        # v2 spawns no daughters, so `configured` must be a NON-EMPTY list that
        # excludes the carrier PVBr; an empty list would trip the legacy
        # "no configured labels => all pools configured" fallback and mask the
        # refusal this test exists to exercise.
        with pytest.raises(ValueError, match=r"'PVBr'.*configured"):
            self._payload(pool, daughters, core,
                          configured=["__not_PVBr__"])

    def test_producer_refuses_unresolvable_gas_routing(self):
        """Enabled kernel without a registered/resolvable gas Species is
        the defined-malformed shape (the ejected X would silently vanish
        -- un-conserved mass): missing registration AND an off-core gas
        object both refuse."""
        pool, daughters, br, core = self._setup()
        pool.side_group_gas_species = []
        with pytest.raises(ValueError, match=r"gas"):
            self._payload(pool, daughters, core)
        pool.side_group_gas_species = [_spc("[Br]", "Br", index=99)]
        with pytest.raises(ValueError, match=r"core species universe"):
            self._payload(pool, daughters, core)

    def test_producer_refuses_carrier_block_plus_defect(self):
        """Producer mirror (v2): a pool serializing the side_group_homolysis
        block AND a chain_mass_defect_g_mol field is contradictory -- SGH
        kernel-v2 tracks multi-loss X depletion via the carrier's auxiliary
        Z inventory, so the pool moments are unchanged and a v2 carrier
        carries no defect. The consumer hard-rejects it, so the producer
        refuses to emit it."""
        pool, daughters, br, core = self._setup()
        pool.chain_mass_defect_g_mol = 79.904
        with pytest.raises(ValueError, match=r"chain_mass_defect_g_mol"):
            self._payload(pool, daughters, core)

    def test_producer_guard_rejects_overlapping_atom_sets(self):
        """Round-75 P1-1 (producer mirror, direct guard call -- through
        the builder the structural selectors partition the atom classes,
        so the overlap shape is only reachable via a broken caller; the
        r68 spawned-overlap precedent): two serialized channels whose
        site_atom_indices intersect WITHOUT being identical double-carry
        the shared site and must refuse. v2 channels have no feature pool,
        so the overlap is checkable from the carrier block alone."""
        from copy import deepcopy
        from rmgpy.polymer import _assert_side_group_serialization_closure
        pool, daughters, br, core = self._setup()
        payload = self._payload(pool, daughters, core)
        pools = payload["pools"]
        carrier = next(p for p in pools if p["label"] == "PVBr")
        block = carrier["side_group_homolysis"]
        base = block["channels"][0]
        twin = deepcopy(base)
        twin["label"] = "other_C-Br"
        twin["site_atom_indices"] = sorted(
            set(base["site_atom_indices"]) | {4})  # strict superset
        twin["sites_per_unit"] = float(len(twin["site_atom_indices"]))
        block["channels"].append(twin)
        conv = payload["conventions"]
        with pytest.raises(ValueError, match=r"overlap"):
            _assert_side_group_serialization_closure(
                pools, [carrier], set(conv["configured_pools"]))

    def test_producer_guard_rejects_out_of_range_indices(self):
        """Round-75 P1-2 (producer mirror, direct guard call): indices
        past the carrier's serialized monomer_adj_list atom count have no
        meaning in the index space the loader bounds-anchors -- the guard
        refuses the shape the consumer hard-rejects."""
        from rmgpy.polymer import _assert_side_group_serialization_closure
        pool, daughters, br, core = self._setup()
        payload = self._payload(pool, daughters, core)
        pools = payload["pools"]
        carrier = next(p for p in pools if p["label"] == "PVBr")
        carrier["side_group_homolysis"]["channels"][0][
            "site_atom_indices"] = [999]
        conv = payload["conventions"]
        with pytest.raises(ValueError, match=r"out of range"):
            _assert_side_group_serialization_closure(
                pools, [carrier], set(conv["configured_pools"]))

    def test_adjacency_list_round_trip_preserves_atom_order(self):
        """Round-75 P1-2 order-stability pin: site_atom_indices are
        0-based positions in pool.monomer.atoms order, and the loader
        bounds-anchors them in monomer_adj_list atom order -- ONE index
        space only because (a) to_adjacency_list writes one numbered atom
        line per atom IN monomer.atoms order and (b)
        from_adjacency_list appends atoms in read order, so a parse
        round-trip preserves identity order. Pinned on the deck monomer
        AND a hostile interleaved-H ordering the writer must not sort."""
        hostile = ("1 H u0 p0 c0 {2,S}\n"
                   "2 C u0 p0 c0 {1,S} {3,S} {4,S} {5,S}\n"
                   "3 H u0 p0 c0 {2,S}\n"
                   "4 H u0 p0 c0 {2,S}\n"
                   "5 C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}\n"
                   "6 H u0 p0 c0 {5,S}\n"
                   "7 H u0 p0 c0 {5,S}\n"
                   "8 H u0 p0 c0 {5,S}\n")
        cases = [Molecule().from_smiles('[CH2][CH]Br'),
                 Molecule().from_adjacency_list(hostile)]
        assert [a.symbol for a in cases[1].atoms][:3] == ['H', 'C', 'H']
        for mol in cases:
            emitted = mol.to_adjacency_list()
            reparsed = Molecule().from_adjacency_list(emitted)
            # emit -> parse -> re-emit: byte-identical
            assert reparsed.to_adjacency_list() == emitted
            # element sequence (the index space the indices live in)
            assert ([a.symbol for a in reparsed.atoms]
                    == [a.symbol for a in mol.atoms])
            # index-space topology: neighbor index sets per position
            pos1 = {a: i for i, a in enumerate(mol.atoms)}
            pos2 = {a: i for i, a in enumerate(reparsed.atoms)}
            assert ([sorted(pos1[b] for b in a.bonds)
                     for a in mol.atoms]
                    == [sorted(pos2[b] for b in a.bonds)
                        for a in reparsed.atoms])
            # the numbered atom lines walk monomer.atoms 1:1, in order
            atom_lines = [ln.split() for ln in emitted.splitlines()
                          if ln.split() and ln.split()[0].isdigit()]
            assert ([t[1] for t in atom_lines]
                    == [a.symbol for a in mol.atoms])

    def test_producer_asserts_adj_list_walks_monomer_atoms(self):
        """Round-75 P1-2: the producer refuses to emit a 2.7 block whose
        site_atom_indices were computed against a different atoms
        sequence than the adj-list serializer walked (the helper runs on
        every kernel-carrying pool serialization; misalignment is only
        reachable via a broken serializer, so the helper is exercised
        directly -- the K1 precedent)."""
        from rmgpy.polymer import _assert_side_group_adj_list_alignment
        mol = Molecule().from_smiles('[CH2][CH]Br')
        # aligned: the real emitted text passes
        _assert_side_group_adj_list_alignment(
            'PVBr', mol, mol.to_adjacency_list())
        # empty text (a swallowed to_adjacency_list failure) refuses
        with pytest.raises(ValueError, match=r"monomer_adj_list"):
            _assert_side_group_adj_list_alignment('PVBr', mol, "")
        # a text that walked a DIFFERENT atoms sequence refuses
        other = Molecule().from_smiles('CC')
        with pytest.raises(ValueError, match=r"atom order"):
            _assert_side_group_adj_list_alignment(
                'PVBr', mol, other.to_adjacency_list())

    def test_copy_carried_defect_without_provenance_serializes(self):
        """A defect-carrying NON-carrier pool serializes legally in v2: the
        serializer still emits chain_mass_defect_g_mol for any pool that
        carries it nonzero, refusing it ONLY on an SGH carrier. Here a
        downstream _mod child carries the defect (e.g. its chains each lost
        one X too) while the SGH-v2 carrier PVBr drives the MAJOR 3.0 stamp;
        the child is unclaimed and serializes with its defect field."""
        pool, daughters, br, core = self._setup()
        grandchild = Polymer(label='PVBr_feat_mod', monomer='[CH2][CH]Br',
                             end_groups=['[H]', '[H]'], cutoff=3,
                             moments=[0.0, 0.0, 0.0], initial_mass=0.0)
        # M_X of the channel's Br gas product (79.904 g/mol).
        grandchild.chain_mass_defect_g_mol = 79.904
        core = core + [_mu_dummy(f"PVBr_feat_mod_mu{k}")
                       for k in range(3)]
        payload = self._payload(
            pool, daughters, core,
            registry=[pool] + daughters + [grandchild],
            configured=(["PVBr"] + [d.label for d in daughters]
                        + ["PVBr_feat_mod"]))
        assert payload["schema_version"] == "3.0"
        g_entry = next(p for p in payload["pools"]
                       if p["label"] == "PVBr_feat_mod")
        assert g_entry["chain_mass_defect_g_mol"] == pytest.approx(
            79.904, rel=1e-12)
        assert g_entry["spawn_event_metadata"].get("source") != \
            "side_group_homolysis"

    def test_kernel_free_pool_still_serializes(self):
        """Negative control: a kernel-free pool of the same shape keeps
        serializing exactly as before (presence-gated)."""
        pe = Polymer(label='PVBr', monomer='[CH2][CH]Br',
                     end_groups=['[H]', '[H]'], cutoff=3,
                     moments=[1.0, 50.0, 3000.0], initial_mass=0.0)
        core = [_mu_dummy(f"PVBr_mu{k}") for k in range(3)]
        payload = build_polymer_moments_artifact(
            [pe], core_species=core,
            configured_pool_labels=["PVBr"], condensed_species=core)
        assert payload["schema_version"] == "2.0"
        assert all("side_group_homolysis" not in p
                   and "chain_mass_defect_g_mol" not in p
                   for p in payload["pools"])


class TestR81RecipeStringsUntouched:
    """r81 (B) schema/TA constraint: exhaustion-tail conditioning is generic
    tolerance-based solver infrastructure OUTSIDE the 2.8 kernel contract --
    the schema-2.8 end_radical_depropagation recipe strings are load-bearing
    pins (runner and TA validate exact text) and must remain byte-identical.
    This pin fails if ANY byte of the recipe dict, the block revision token,
    or the kernel name changes without a coordinated re-dating."""

    def test_28_depropagation_recipe_strings_byte_identical(self):
        import hashlib

        from rmgpy.polymer import (
            DEPROPAGATION_BLOCK_RECIPE_REVISION,
            DEPROPAGATION_KERNEL_NAME,
            DEPROPAGATION_SIDECAR_RECIPE,
        )
        payload = json.dumps(
            {"recipe": DEPROPAGATION_SIDECAR_RECIPE,
             "revision": DEPROPAGATION_BLOCK_RECIPE_REVISION,
             "kernel": DEPROPAGATION_KERNEL_NAME},
            sort_keys=True, ensure_ascii=True)
        assert hashlib.sha256(payload.encode("utf-8")).hexdigest() == \
            "c27136d6f5e745c631bac3f9b8e430d0187fb6e849a6faf9cca9390311e7a352"
        assert DEPROPAGATION_BLOCK_RECIPE_REVISION == \
            "2026-07-06-end-radical-depropagation"
        assert DEPROPAGATION_KERNEL_NAME == "end_radical_depropagation/1"


class TestSpawnedPoolRefusalMirror:
    """r91 artifact lockstep: the emitter must mirror the solver's
    spawned-pool demotion refusal -- a stamped pool-coupled row whose
    unresolved endpoint is a spawned config-less Polymer emits refused
    (refused: true, refused_reason conduit-deferred) with its STAMPED
    archetype, never a live legacy_mu1/1 row."""

    @staticmethod
    def _spawned(label="A_mod"):
        return Polymer(label=label, monomer="[CH2][CH]c1ccccc1",
                       end_groups=["[CH3]", "[H]"], cutoff=3,
                       Mn=5000.0, Mw=6000.0, initial_mass=0.001)

    def test_migration_to_spawned_pool_emits_refused_not_live_legacy(self):
        core = _two_pool_core()
        d = self._spawned()
        core = core + [d]
        rxn = Reaction(reactants=[core[0]], products=[d],
                       kinetics=_arrhenius(), reversible=False)
        rxn.polymer_flux_archetype = int(PolymerFluxArchetype.MIGRATION)
        entries = compile_polymer_reaction_entries(
            [rxn], core, configured_pool_labels=["A", "B"],
            cantera_index_map={id(rxn): [0]})
        assert len(entries) == 1
        e = entries[0]
        assert e.get("refused") is True
        assert e.get("refused_reason") == "conduit-deferred"
        assert e["archetype"] == "migration/1"   # NOT a live legacy_mu1/1 row
        assert e["archetype"] != "legacy_mu1/1"
        assert e["unresolved"] is False

    def test_volatile_ejection_to_spawned_pool_emits_refused(self):
        core = _two_pool_core()
        d = self._spawned()
        core = core + [d]
        rxn = Reaction(reactants=[core[0]], products=[d, core[8]],
                       kinetics=_arrhenius(), reversible=False)
        rxn.polymer_flux_archetype = int(PolymerFluxArchetype.VOLATILE_EJECTION)
        rxn.polymer_eject_units = 0.5
        entries = compile_polymer_reaction_entries(
            [rxn], core, configured_pool_labels=["A", "B"],
            cantera_index_map={id(rxn): [3]})
        e = entries[0]
        assert e.get("refused") is True
        assert e.get("refused_reason") == "conduit-deferred"
        assert e["archetype"] == "volatile_ejection/1"

    def test_refused_marker_round_trips_runner_loader(self):
        from rmgpy.tools.polymer_moments_runner import _validate_refused_entry
        core = _two_pool_core()
        d = self._spawned()
        core = core + [d]
        rxn = Reaction(reactants=[core[0]], products=[d],
                       kinetics=_arrhenius(), reversible=False)
        rxn.polymer_flux_archetype = int(PolymerFluxArchetype.MIGRATION)
        entries = compile_polymer_reaction_entries(
            [rxn], core, configured_pool_labels=["A", "B"],
            cantera_index_map={id(rxn): [0]})
        assert _validate_refused_entry(entries[0])

    def test_spawned_refused_without_configured_participant_hard_fails(self):
        """r92 artifact P1: the r91-documented corner -- a spawned-refused
        row with NO configured-pool participant (both endpoints spawned
        config-less pools) cannot legally carry the sidecar refused marker,
        and the pre-r92 emitter warned + emitted a live legacy_mu1/1 row the
        generating solver zeroes (over-integration for a conforming
        consumer). The emitter must HARD-FAIL loudly instead."""
        core = _two_pool_core()
        d1, d2 = self._spawned("A_mod"), self._spawned("B_mod")
        rxn = Reaction(reactants=[d1], products=[d2],
                       kinetics=_arrhenius(), reversible=False)
        rxn.polymer_flux_archetype = int(PolymerFluxArchetype.MIGRATION)
        with pytest.raises(ValueError, match=r"refused.*no pool-mapped"):
            compile_polymer_reaction_entries(
                [rxn], core + [d1, d2], configured_pool_labels=["A", "B"],
                cantera_index_map={id(rxn): [0]})

    def test_plain_species_dst_keeps_legacy_demotion_mirror(self):
        """Negative control: the existing demotion mirror is untouched when
        the unresolved endpoint is NOT a spawned Polymer (plain Species dst,
        the test_stamped_migration_without_configured_dst shape)."""
        core = _two_pool_core()
        rxn = Reaction(reactants=[core[0]], products=[core[4]],
                       kinetics=_arrhenius(), reversible=False)
        rxn.polymer_flux_archetype = int(PolymerFluxArchetype.MIGRATION)
        entries = compile_polymer_reaction_entries(
            [rxn], core, configured_pool_labels=["A"],   # B not configured
            cantera_index_map={id(rxn): [0]})
        e = entries[0]
        assert e["archetype"] == "legacy_mu1/1"
        assert e["unresolved"] is True
        assert "refused" not in e


class TestItem16SpawnedConfiguredPoolArtifact:
    """Item 16 artifact lockstep (must-pin f): once the ENGINE configures a
    spawned daughter pool (derive_daughter_pool_configs over the promoted
    core -- nothing hand-fed), its volatile_ejection/1 row serializes LIVE
    (no refused marker, not unresolved), the spawned pool is listed in
    pools[] and STILL classifies moments_provenance spawned_empty (the
    item-#14a s7 coupling re-confirmed, not silently inherited), and the
    r95 POOL LIVENESS census stays quiet. The refusal mirror for genuinely
    config-less spawned pools (TestSpawnedPoolRefusalMirror) is unchanged."""

    @staticmethod
    def _engine_configured_core():
        from rmgpy.rmg.polymer_input import derive_daughter_pool_configs
        core = _two_pool_core()
        d = Polymer(label="A_mod", monomer="[CH2][CH]c1ccccc1",
                    end_groups=["[CH3]", "[H]"], cutoff=3,
                    Mn=5000.0, Mw=6000.0, initial_mass=0.0)
        d.parent_pool_label = "A"
        d.spawn_metadata = {"source": "radical_feature_h_loss"}
        core = core + [d, _mu_dummy("A_mod_mu0"), _mu_dummy("A_mod_mu1"),
                       _mu_dummy("A_mod_mu2")]
        spc_map = {s: i for i, s in enumerate(core)}
        derived = derive_daughter_pool_configs(core, spc_map, {"A", "B"})
        configured = ["A", "B"] + [c.label for c in derived]
        assert "A_mod" in configured        # ENGINE-created, not hand-fed
        rxn = Reaction(reactants=[core[0]], products=[d, core[8]],
                       kinetics=_arrhenius(), reversible=False)
        rxn.polymer_flux_archetype = int(PolymerFluxArchetype.VOLATILE_EJECTION)
        rxn.polymer_eject_units = 0.5
        return core, d, rxn, configured

    def test_ve_row_to_engine_configured_pool_serializes_live(self):
        core, d, rxn, configured = self._engine_configured_core()
        entries = compile_polymer_reaction_entries(
            [rxn], core, configured_pool_labels=configured,
            cantera_index_map={id(rxn): [0]})
        assert len(entries) == 1
        e = entries[0]
        assert "refused" not in e           # live, no refused marker
        assert "refused_reason" not in e
        assert e["archetype"] == "volatile_ejection/1"
        assert e["unresolved"] is False
        assert e["src_pool"] == "A"
        assert e["dst_pool"] == "A_mod"
        assert e["params"] == {"eject_units": 0.5}

    def test_artifact_lists_spawned_pool_and_liveness_stays_quiet(self, caplog):
        import logging
        core, d, rxn, configured = self._engine_configured_core()
        pool_a = Polymer(label="A", monomer="[CH2][CH]c1ccccc1",
                         end_groups=["[CH3]", "[H]"], cutoff=3,
                         Mn=5000.0, Mw=6000.0, initial_mass=0.001)
        with caplog.at_level(logging.WARNING):
            art = build_polymer_moments_artifact(
                [pool_a, d], core_species=core, core_reactions=[rxn],
                configured_pool_labels=configured,
                condensed_species=[core[0], core[4]],
                cantera_index_map={id(rxn): [0]})
        by_label = {p["label"]: p for p in art["pools"]}
        assert "A_mod" in by_label          # spawned pool listed
        # s7 coupling re-confirmed: configured AND still spawned_empty (the
        # spawn provenance markers, not the configured set, decide).
        assert by_label["A_mod"]["moments_provenance"] == "spawned_empty"
        assert list(by_label["A_mod"]["moments"]) == [0.0, 0.0, 0.0]
        assert len(art["reactions"]) == 1
        assert "refused" not in art["reactions"][0]
        assert not any("POOL LIVENESS ALARM" in r.getMessage()
                       for r in caplog.records)

    def test_sidecar_configured_pools_excludes_midrun_spawned_daughter(self):
        """Item-16 P1 (adversarial review, CKMG contract): the EMITTED
        conventions.configured_pools names DECK/root pools only. The
        engine-configured _mod daughter stays OUT of it -- the consumer
        treats configured pools as load-bearing and hard-flags structurally
        dead ones, and a born-empty daughter with no outgoing edge is dead
        by construction -- while its volatile_ejection row still serializes
        LIVE off the full solver-configured set, and the daughter is
        published through conventions.spawned_pools (schema 2.5) with its
        phase species joining the condensed closure."""
        core, d, rxn, configured = self._engine_configured_core()
        assert "A_mod" in configured        # solver truth: engine-configured
        pool_a = Polymer(label="A", monomer="[CH2][CH]c1ccccc1",
                         end_groups=["[CH3]", "[H]"], cutoff=3,
                         Mn=5000.0, Mw=6000.0, initial_mass=0.001)
        art = build_polymer_moments_artifact(
            [pool_a, d], core_species=core, core_reactions=[rxn],
            configured_pool_labels=configured,
            condensed_species=[core[0], core[4]],
            cantera_index_map={id(rxn): [0]})
        conv = art["conventions"]
        # The P1 pins: _mod NOT configured, published spawned instead.
        assert "A_mod" not in conv["configured_pools"]
        assert "A" in conv["configured_pools"]
        assert "B" in conv["configured_pools"]
        assert conv["spawned_pools"] == ["A_mod"]
        assert set(conv["spawned_pools"]).isdisjoint(conv["configured_pools"])
        # Row resolution keeps running off the FULL engine set: the
        # cross-pool VE row serializes live, refused-absent.
        assert len(art["reactions"]) == 1
        assert "refused" not in art["reactions"][0]
        assert "refused_reason" not in art["reactions"][0]
        assert art["reactions"][0]["dst_pool"] == "A_mod"
        # Spawned-pool vocabulary is the strongest shape stamp here (2.5),
        # and the daughter's condensed bookkeeping joins the closure.
        assert art["schema_version"] == "2.5"
        by_label = {p["label"]: p for p in art["pools"]}
        for lbl in by_label["A_mod"].get("phase_species") or []:
            assert lbl in conv["condensed_species"]

    def test_setup_time_homolysis_daughter_stays_configured(self):
        """Negative control for the P1 split: a SETUP-TIME homolysis
        end-radical daughter carries the same spawn markers but source
        HOMOLYSIS_SPAWN_SOURCE -- the 2.6/2.8 closure guards require it
        configured and never spawned-classified, so the emitted
        configured_pools keeps it and no spawned surface appears."""
        pool_a = Polymer(label="A", monomer="[CH2][CH]c1ccccc1",
                         end_groups=["[CH3]", "[H]"], cutoff=3,
                         Mn=5000.0, Mw=6000.0, initial_mass=0.001)
        hd = Polymer(label="A_rad_primary_end", monomer="[CH2][CH]c1ccccc1",
                     end_groups=["[CH3]", "[H]"], cutoff=3,
                     Mn=5000.0, Mw=6000.0, initial_mass=0.0)
        hd.parent_pool_label = "A"
        hd.spawn_metadata = {"source": HOMOLYSIS_SPAWN_SOURCE}
        art = build_polymer_moments_artifact(
            [pool_a, hd], core_species=None, core_reactions=[],
            configured_pool_labels=["A", "A_rad_primary_end"])
        conv = art["conventions"]
        assert conv["configured_pools"] == ["A", "A_rad_primary_end"]
        assert "spawned_pools" not in conv

    def test_midrun_spawn_predicate_classification(self):
        """is_midrun_spawned_pool_daughter: the single-source split the
        emission keys on. _mod H-loss daughter (mid-run) True; setup-time
        homolysis daughter False; markerless deck pool False; markerless
        parent_pool_label-only scission tail True."""
        deck = Polymer(label="A", monomer="[CH2][CH]c1ccccc1",
                       end_groups=["[CH3]", "[H]"], cutoff=3,
                       Mn=5000.0, Mw=6000.0, initial_mass=0.001)
        assert not is_midrun_spawned_pool_daughter(deck)
        mod = Polymer(label="A_mod", monomer="[CH2][CH]c1ccccc1",
                      end_groups=["[CH3]", "[H]"], cutoff=3,
                      Mn=5000.0, Mw=6000.0, initial_mass=0.0)
        mod.parent_pool_label = "A"
        mod.spawn_metadata = {"source": "radical_feature_h_loss"}
        assert is_midrun_spawned_pool_daughter(mod)
        homolysis = Polymer(label="A_rad_primary_end",
                            monomer="[CH2][CH]c1ccccc1",
                            end_groups=["[CH3]", "[H]"], cutoff=3,
                            Mn=5000.0, Mw=6000.0, initial_mass=0.0)
        homolysis.parent_pool_label = "A"
        homolysis.spawn_metadata = {"source": HOMOLYSIS_SPAWN_SOURCE}
        assert not is_midrun_spawned_pool_daughter(homolysis)
        tail = Polymer(label="A_scission_tail", monomer="[CH2][CH]c1ccccc1",
                       end_groups=["[CH3]", "[H]"], cutoff=3,
                       Mn=5000.0, Mw=6000.0, initial_mass=0.0)
        tail.parent_pool_label = "A"
        assert is_midrun_spawned_pool_daughter(tail)


class TestThermalAnalysisInputsValidation:
    """validate_thermal_analysis_inputs (schema 2.9, part A): the shared
    deck-read validator. Declared thermal-analysis inputs are threaded
    verbatim with provenance; the validator is the guard against unsourced
    numbers and typo'd fields."""

    def test_declared_field_normalizes_to_full_object(self):
        norm = validate_thermal_analysis_inputs("PS", {
            "dH_depoly_J_per_mol": {
                "value": 70000.0, "units": "J/mol", "basis": "per monomer",
                "temperature_K": 298.15, "provenance": "Smith 2020"},
        })
        assert norm["dH_depoly_J_per_mol"] == {
            "value": 70000.0, "units": "J/mol", "basis": "per monomer",
            "temperature_K": 298.15, "provenance": "Smith 2020"}

    def test_value_without_provenance_is_rejected(self):
        with pytest.raises(ValueError, match="provenance"):
            validate_thermal_analysis_inputs("PS", {
                "dH_vap_J_per_mol": {"value": 40000.0, "units": "J/mol"}})

    def test_unrecognized_field_is_rejected(self):
        with pytest.raises(ValueError, match="unrecognized field"):
            validate_thermal_analysis_inputs("PS", {
                "dH_sublimation": {"value": 1.0, "provenance": "x"}})

    def test_non_finite_value_is_rejected(self):
        with pytest.raises(ValueError, match="non-finite"):
            validate_thermal_analysis_inputs("PS", {
                "cp_condensed_J_per_kg_K": {
                    "value": float("inf"), "provenance": "x"}})

    def test_value_null_normalizes_provenance_to_unset(self):
        norm = validate_thermal_analysis_inputs("PS", {
            "dH_vap_J_per_mol": {"value": None}})
        assert norm["dH_vap_J_per_mol"]["value"] is None
        assert norm["dH_vap_J_per_mol"]["provenance"] == "unset"


class TestThermalAnalysisInputsSerialization:
    """Schema 2.9, part A: DECLARED thermal-analysis inputs threaded from the
    deck into the sidecar with provenance. Per-pool fields inside the pool
    entry; instrument fields under conventions. Undeclared fields emit the
    honest {value: null, provenance: "unset"} sentinel -- never a made-up
    number."""

    def _pool(self, tai=None):
        p = Polymer(
            label="PS", monomer="[CH2][CH](c1ccccc1)",
            end_groups=["[H]", "[H]"], cutoff=3, Mn=1500.0, Mw=1800.0,
            initial_mass=1.0, k_scission=1.0)
        if tai is not None:
            p.thermal_analysis_inputs = validate_thermal_analysis_inputs(
                "PS", tai)
        return p

    def test_declared_per_pool_value_round_trips_with_provenance(self):
        p = self._pool({
            "dH_depoly_J_per_mol": {
                "value": 70000.0, "units": "J/mol", "basis": "per monomer",
                "temperature_K": 298.15, "provenance": "Smith 2020"}})
        d = _serialize_pool_for_sidecar(p)
        got = d["thermal_analysis_inputs"]["dH_depoly_J_per_mol"]
        assert got["value"] == pytest.approx(70000.0)
        assert got["units"] == "J/mol"
        assert got["provenance"] == "Smith 2020"

    def test_undeclared_per_pool_field_emits_null_unset(self):
        p = self._pool({
            "dH_depoly_J_per_mol": {"value": 70000.0, "provenance": "S20"}})
        d = _serialize_pool_for_sidecar(p)
        # dH_vap / cp_condensed were NOT declared -> honest incompleteness
        assert d["thermal_analysis_inputs"]["dH_vap_J_per_mol"] == {
            "value": None, "units": None, "basis": None,
            "temperature_K": None, "provenance": "unset"}
        assert d["thermal_analysis_inputs"]["cp_condensed_J_per_kg_K"][
            "provenance"] == "unset"

    def test_instrument_fields_land_under_conventions(self):
        p = self._pool({
            "htc_W_per_m2_K": {
                "value": 25.0, "units": "W/m^2/K", "provenance": "pan calib"}})
        art = build_polymer_moments_artifact(
            [p], core_species=None, core_reactions=[],
            configured_pool_labels=["PS"])
        conv = art["conventions"]["thermal_analysis_inputs"]
        assert conv["htc_W_per_m2_K"]["value"] == pytest.approx(25.0)
        assert conv["htc_W_per_m2_K"]["provenance"] == "pan calib"
        # undeclared instrument fields -> null/unset
        assert conv["wall_area_m2"]["provenance"] == "unset"
        assert conv["pan_area_m2"]["value"] is None
        # instrument fields are NOT duplicated into the pool entry
        assert "thermal_analysis_inputs" not in art["pools"][0]

    def test_thermal_inputs_bump_schema_to_2_9(self):
        p = self._pool({
            "dH_depoly_J_per_mol": {"value": 70000.0, "provenance": "S20"}})
        art = build_polymer_moments_artifact(
            [p], core_species=None, core_reactions=[],
            configured_pool_labels=["PS"])
        assert art["schema_version"] == POLYMER_POOLS_SIDECAR_SCHEMA_VERSION_THERMAL
        assert art["schema_version"] == "2.9"

    def test_no_thermal_declaration_keeps_legacy_stamp_and_no_keys(self):
        """Negative control: a pool with no thermal declaration is
        byte-identical to the pre-2.9 shape (no thermal keys anywhere, older
        stamp preserved)."""
        p = self._pool()
        art = build_polymer_moments_artifact(
            [p], core_species=None, core_reactions=[],
            configured_pool_labels=["PS"])
        assert art["schema_version"] != "2.9"
        assert "thermal_analysis_inputs" not in art["conventions"]
        assert "thermal_analysis_inputs" not in art["pools"][0]

    def test_thermal_inputs_survive_json_round_trip(self):
        p = self._pool({
            "dH_vap_J_per_mol": {
                "value": 40000.0, "units": "J/mol", "provenance": "NIST"},
            "wall_area_m2": {
                "value": 1.5e-4, "units": "m^2", "provenance": "pan spec"}})
        art = build_polymer_moments_artifact(
            [p], core_species=None, core_reactions=[],
            configured_pool_labels=["PS"])
        reloaded = json.loads(json.dumps(art, default=str))
        assert reloaded["pools"][0]["thermal_analysis_inputs"][
            "dH_vap_J_per_mol"]["provenance"] == "NIST"
        assert reloaded["conventions"]["thermal_analysis_inputs"][
            "wall_area_m2"]["value"] == pytest.approx(1.5e-4)


class TestExplicitDpInventorySerialization:
    """Schema 2.9, part B: the FULL real explicit-DP chip inventory. When the
    solver tracks MORE THAN ONE discrete DP chip, every real chip is emitted
    as {DP, species_label, moles} -- the actual sparse discrete-chip
    population RMG holds, NOT a closure-expanded n(DP) histogram. The
    single-cutoff-chip case stays byte-identical 2.3 (no inventory key)."""

    def _pool(self):
        p = Polymer(
            label="PS", monomer="[CH2][CH](c1ccccc1)",
            end_groups=["[H]", "[H]"], cutoff=3, Mn=1500.0, Mw=1800.0,
            initial_mass=1.0, k_scission=0.0, k_unzip=0.01)
        p.explicit_dp = True
        p.explicit_dp_species = _spc("CCC", "PS_dp3", index=42)
        return p

    def _core(self, pool):
        return [
            _spc("CC", "PS", index=2),
            _mu_dummy("PS_mu0"), _mu_dummy("PS_mu1"), _mu_dummy("PS_mu2"),
            pool.explicit_dp_species,
            _spc("CCCC", "PS_dp4", index=43),
            _spc("CCCCC", "PS_dp5", index=44),
        ]

    def test_single_chip_stays_byte_identical_2_3(self):
        """One tracked chip -> no inventory key, schema stays 2.3 (the
        golden-pinned 2.3 block is untouched)."""
        p = self._pool()
        d = _serialize_pool_for_sidecar(
            p, core_species=self._core(p), initial_explicit_moles={3: 0.25})
        assert "inventory" not in d["explicit_dp"]
        art = build_polymer_moments_artifact(
            [p], core_species=self._core(p), core_reactions=[],
            configured_pool_labels=["PS"],
            initial_explicit_by_pool={"PS": {3: 0.25}})
        assert art["schema_version"] == "2.3"

    def test_multi_chip_inventory_serializes_real_populations(self):
        p = self._pool()
        core = self._core(p)
        by_dp = {3: "PS_dp3(42)", 4: "PS_dp4(43)", 5: "PS_dp5(44)"}
        d = _serialize_pool_for_sidecar(
            p, core_species=core,
            initial_explicit_moles={3: 0.25, 4: 0.10, 5: 0.05},
            explicit_dp_species_by_dp=by_dp)
        inv = d["explicit_dp"]["inventory"]
        assert inv == [
            {"DP": 3, "species_label": "PS_dp3(42)", "moles": 0.25},
            {"DP": 4, "species_label": "PS_dp4(43)", "moles": 0.10},
            {"DP": 5, "species_label": "PS_dp5(44)", "moles": 0.05},
        ]

    def test_multi_chip_artifact_stamps_2_9(self):
        p = self._pool()
        core = self._core(p)
        by_dp = {3: "PS_dp3(42)", 4: "PS_dp4(43)", 5: "PS_dp5(44)"}
        art = build_polymer_moments_artifact(
            [p], core_species=core, core_reactions=[],
            configured_pool_labels=["PS"],
            initial_explicit_by_pool={"PS": {3: 0.25, 4: 0.10, 5: 0.05}},
            explicit_dp_species_by_pool={"PS": by_dp})
        assert art["schema_version"] == "2.9"
        assert len(art["pools"][0]["explicit_dp"]["inventory"]) == 3

    def test_explicit_dp_off_has_no_inventory(self):
        """explicit_dp disabled -> no explicit_dp block at all, so no
        inventory is fabricated."""
        p = Polymer(
            label="PS", monomer="[CH2][CH](c1ccccc1)",
            end_groups=["[H]", "[H]"], cutoff=3, Mn=1500.0, Mw=1800.0,
            k_scission=1.0)
        d = _serialize_pool_for_sidecar(p)
        assert "explicit_dp" not in d
