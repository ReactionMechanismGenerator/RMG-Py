#!/usr/bin/env python3
"""Sufficiency proof (design spec §9 tests 1 + 2b): a numpy-only consumer of
the schema-2.0 artifact reproduces HybridPolymerSystem trajectories.

The ORACLE side (this file) may import rmgpy freely; the CONSUMER
(test/rmgpy/tools/numpy_moments_consumer.py) must not — enforced below.
Both sides integrate with identical fixed-step forward Euler (the oracle
pattern from test/rmgpy/solver/solverPolymerTest.py:498-510)."""

import json
import os
import sys

import numpy as np
import pytest

from rmgpy.kinetics import Arrhenius
from rmgpy.molecule import Molecule
from rmgpy.polymer import Polymer, PolymerFluxArchetype, build_polymer_moments_artifact
from rmgpy.reaction import Reaction
from rmgpy.solver.polymer import HybridPolymerSystem, MassTransferConfig, PolymerPoolConfig
from rmgpy.species import Species
from rmgpy.thermo import NASA, NASAPolynomial

# The test dirs are not packages (no __init__.py anywhere under test/);
# import the consumer module by path.
_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
from numpy_moments_consumer import ArtifactConsumer  # noqa: E402

T_K = 800.0
P_PA = 1.0e5
V_POLY = 1.0
DT = 1.0e-4
N_STEPS = 2000


def test_consumer_module_is_rmgpy_free():
    """THE point of the artifact: the consumer never imports rmgpy.
    Checked structurally (AST import table), not by substring — the module
    docstring legitimately mentions rmgpy by name."""
    import ast
    src_path = os.path.join(_HERE, "numpy_moments_consumer.py")
    with open(src_path) as fh:
        tree = ast.parse(fh.read())
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            assert all(not a.name.split(".")[0] == "rmgpy" for a in node.names)
        if isinstance(node, ast.ImportFrom):
            assert (node.module or "").split(".")[0] != "rmgpy"


def _spc(smiles, label, index=-1):
    s = Species(molecule=[Molecule().from_smiles(smiles)])
    s.label = label
    s.index = index
    return s


def _mu(label):
    s = Species(label=label, reactive=False)
    s.molecule = [Molecule().from_smiles("[Ne]")]
    s.is_moment_dummy = True
    s.index = -1
    return s


def _yaml_label(s):
    """chem.yaml naming (== artifact labels, rmgpy.cantera.get_label):
    'label(index)' when index > 0, bare label otherwise. The consumer's
    species_order MUST use these names — they key every artifact entry."""
    return f"{s.label}({s.index})" if s.index > 0 else s.label


def _nasa(a5, a6):
    rows = []
    for tmin, tmax in ((200.0, 1000.0), (1000.0, 6000.0)):
        rows.append(NASAPolynomial(coeffs=[2.5, 0.0, 0.0, 0.0, 0.0, a5, a6],
                                   Tmin=(tmin, "K"), Tmax=(tmax, "K")))
    return NASA(polynomials=rows, Tmin=(200.0, "K"), Tmax=(6000.0, "K"))


def _euler_oracle(rs, y0, dt, n_steps):
    y = y0.copy()
    traj = [y.copy()]
    for k in range(n_steps):
        dn = rs.residual(k * dt, y, np.zeros_like(y))[0]
        y = y + dt * dn
        traj.append(y.copy())
    return np.array(traj)


# ---------------------------------------------------------------------------
# (a) channel ODEs: pure scission and pure unzip vs oracle + analytic forms
# ---------------------------------------------------------------------------

class TestChannelSufficiency:
    def _run(self, k_s, k_u):
        inert = _spc("N#N", "N2")
        # GAS released-monomer target M (recipe revision
        # 2026-07-03-monomer-gas: the solver validates the routing target
        # gas): k_unzip > 0 requires a wired emission target on BOTH sides
        # (the solver invariant refuses monomer_poly_index=None, because the
        # drain is unconditional while the emission is index-gated; the
        # consumer mirrors it with dn[routing] += k_u*mu0*Vp) -- so the
        # oracle comparison now covers the release flux too.
        monomer = _spc("C", "M")
        core = [inert, _mu("poly_mu0"), _mu("poly_mu1"), _mu("poly_mu2"),
                monomer]
        mask = np.array([True, False, False, False, True], dtype=bool)
        pool_cfg = PolymerPoolConfig(label="poly", xs=2,
                                     explicit_dp_to_species_index={},
                                     mu_indices=(1, 2, 3), monomer_poly_index=4,
                                     k_scission=k_s, k_unzip=k_u, tail_kinetics=None)
        mom0 = (1.0, 5.0, 30.0)
        rs = HybridPolymerSystem(
            T=T_K, P=P_PA, initial_mole_fractions={inert: 1.0}, V_poly=V_POLY,
            polymer_pools=[pool_cfg], mass_transfer=[], gas_species_mask=mask,
            constant_gas_volume=False,
            initial_polymer_moments={"poly": mom0}, termination=[])
        rs.initialize_model(core, [], [], [])
        y0 = rs.y.copy()
        oracle = _euler_oracle(rs, y0, DT, N_STEPS)

        registry_pool = Polymer(label="poly", monomer="[CH2][CH2]",
                                end_groups=["[H]", "[H]"], cutoff=3,
                                moments=list(mom0), initial_mass=0.0,
                                k_scission=k_s, k_unzip=k_u)
        artifact = build_polymer_moments_artifact(
            [registry_pool], core_species=core, core_reactions=[],
            configured_pool_labels=["poly"],
            # routed monomer M is GAS (2026-07-03-monomer-gas): only the
            # mu dummies are condensed.
            condensed_species=core[1:4],
            monomer_routing_by_pool={"poly": _yaml_label(monomer)},
            cantera_index_map={})
        artifact = json.loads(json.dumps(artifact))  # honest serialization

        consumer = ArtifactConsumer(artifact, [_yaml_label(s) for s in core],
                                    P=P_PA, V_poly=V_POLY)
        _, mine = consumer.integrate_euler(y0, T_K, DT, N_STEPS)
        return oracle, mine

    def test_pure_scission_matches_oracle_and_analytic(self):
        oracle, mine = self._run(k_s=1.0, k_u=0.0)
        np.testing.assert_allclose(mine, oracle, rtol=1e-9, atol=1e-12)
        # analytic: mu1 const; mu0(t) = mu1 - (mu1 - mu0_0) e^{-k t}
        t_end = DT * N_STEPS
        mu0_exp = 5.0 - (5.0 - 1.0) * np.exp(-1.0 * t_end)
        # Euler at dt=1e-4 over 0.2 time constants: ~1e-4 relative accuracy
        assert mine[-1, 1] == pytest.approx(mu0_exp, rel=5e-4)
        assert mine[-1, 2] == pytest.approx(5.0, rel=1e-9)  # mu1 conserved

    def test_pure_unzip_matches_oracle(self):
        oracle, mine = self._run(k_s=0.0, k_u=0.5)
        np.testing.assert_allclose(mine, oracle, rtol=1e-9, atol=1e-12)
        # analytic: mu0 const, mu1(t) = mu1_0 - k*mu0*t
        t_end = DT * N_STEPS
        assert mine[-1, 2] == pytest.approx(5.0 - 0.5 * 1.0 * t_end, rel=1e-6)


class TestUnknownChannelRejection:
    def test_consumer_rejects_qssa_channel_loudly(self):
        """This consumer implements scission/unzip ONLY. An artifact whose
        pool carries a radical_qssa_unzip channel (or any unknown channel
        kind) must fail at construction, not silently integrate without the
        channel -- the silent path produces a flat/false trajectory, the
        exact failure class the sidecar schema exists to kill. Mirrors TA's
        guard (~/Code/TA/ta/mechanism.py:509-517)."""
        inert = _spc("N#N", "N2")
        monomer = _spc("C", "M")
        core = [inert, _mu("poly_mu0"), _mu("poly_mu1"), _mu("poly_mu2"),
                monomer]
        registry_pool = Polymer(
            label="poly", monomer="[CH2][CH2]",
            end_groups=["[H]", "[H]"], cutoff=3,
            moments=[1.0, 5.0, 30.0], initial_mass=0.0,
            k_scission=0.0, k_unzip=0.0,
            radical_qssa_unzip={
                "initiation": {"A": 1.0e13, "n": 0.0, "Ea": 3.0e5},
                "depropagation": {"A": 1.0e14, "n": 0.0, "Ea": 9.0e4},
                "termination": {"A": 1.0e8, "n": 0.0, "Ea": 1.0e4},
            })
        artifact = build_polymer_moments_artifact(
            [registry_pool], core_species=core, core_reactions=[],
            configured_pool_labels=["poly"],
            # routed monomer M is GAS (2026-07-03-monomer-gas): only the
            # mu dummies are condensed.
            condensed_species=core[1:4],
            monomer_routing_by_pool={"poly": _yaml_label(monomer)},
            cantera_index_map={})
        artifact = json.loads(json.dumps(artifact))
        with pytest.raises(ValueError, match=r"poly.*radical_qssa_unzip"):
            ArtifactConsumer(artifact, [_yaml_label(s) for s in core],
                             P=P_PA, V_poly=V_POLY)


# ---------------------------------------------------------------------------
# (b) synthetic two-pool deck: stamped MIGRATION + DISCRETE_CHIP (+ reversible
#     Keq case). Reuses the apportionment plan's two-pool fixture layout.
# ---------------------------------------------------------------------------

def _two_pool_setup(with_thermo=False):
    sp = {
        "A": _spc("CCCC", "A", index=1),
        "A_mu0": _mu("A_mu0"), "A_mu1": _mu("A_mu1"), "A_mu2": _mu("A_mu2"),
        "B": _spc("CCCCC", "B", index=5),
        "B_mu0": _mu("B_mu0"), "B_mu1": _mu("B_mu1"), "B_mu2": _mu("B_mu2"),
        "G": _spc("[CH3]", "G", index=9),
    }
    sp["A"].is_polymer_proxy = True
    sp["B"].is_polymer_proxy = True
    core = [sp["A"], sp["A_mu0"], sp["A_mu1"], sp["A_mu2"],
            sp["B"], sp["B_mu0"], sp["B_mu1"], sp["B_mu2"], sp["G"]]
    mask = np.array([False] * 8 + [True], dtype=bool)
    if with_thermo:
        sp["A"].thermo = _nasa(-1000.0, 4.0)
        sp["B"].thermo = _nasa(-1200.0, 3.0)
        sp["G"].thermo = _nasa(-500.0, 5.0)
    return sp, core, mask


def _pools_ab():
    pool_a = PolymerPoolConfig(label="A", xs=2, explicit_dp_to_species_index={},
                               mu_indices=(1, 2, 3), monomer_poly_index=None,
                               k_scission=0.0, k_unzip=0.0, tail_kinetics=None)
    pool_b = PolymerPoolConfig(label="B", xs=2, explicit_dp_to_species_index={},
                               mu_indices=(5, 6, 7), monomer_poly_index=None,
                               k_scission=0.0, k_unzip=0.0, tail_kinetics=None)
    return [pool_a, pool_b]


def _registry_ab():
    out = []
    for lab in ("A", "B"):
        out.append(Polymer(label=lab, monomer="[CH2][CH2]",
                           end_groups=["[H]", "[H]"], cutoff=3,
                           moments=[1.0, 5.0, 30.0], initial_mass=0.0))
    return out


def _oracle_system(core, mask, reactions, mass_transfer=()):
    rs = HybridPolymerSystem(
        T=T_K, P=P_PA, initial_mole_fractions={core[8]: 0.05}, V_poly=V_POLY,
        polymer_pools=_pools_ab(), mass_transfer=list(mass_transfer),
        gas_species_mask=mask.copy(), constant_gas_volume=False,
        initial_polymer_moments={"A": (1.0, 5.0, 30.0), "B": (2.0, 4.0, 10.0)},
        termination=[])
    rs.initialize_model(core, list(reactions), [], [])
    return rs


def _artifact_for(core, reactions):
    artifact = build_polymer_moments_artifact(
        _registry_ab(), core_species=core, core_reactions=list(reactions),
        configured_pool_labels=["A", "B"],
        condensed_species=core[:8], cantera_index_map={})
    return json.loads(json.dumps(artifact))


def _kin1():
    return Arrhenius(A=(2.0, "1/s"), n=0.0, Ea=(0.0, "kcal/mol"), T0=(298.15, "K"))


class TestTwoPoolSufficiency:
    def test_migration_plus_discrete_chip_deck(self):
        """Spec §9 test 1(b): one stamped MIGRATION + one stamped DISCRETE_CHIP
        reaction in the SAME deck; consumer matches oracle on every species
        and every pool moment."""
        sp, core, mask = _two_pool_setup()
        mig = Reaction(reactants=[sp["A"]], products=[sp["B"]],
                       kinetics=_kin1(), reversible=False)
        mig.polymer_flux_archetype = int(PolymerFluxArchetype.MIGRATION)
        chip = Reaction(reactants=[sp["A"]], products=[sp["A"], sp["G"]],
                        kinetics=_kin1(), reversible=False)
        chip.polymer_flux_archetype = int(PolymerFluxArchetype.DISCRETE_CHIP)
        chip.polymer_chip_units = 2
        chip.is_end_group_reaction = True

        rxns = [mig, chip]
        rs = _oracle_system(core, mask, rxns)
        y0 = rs.y.copy()
        oracle = _euler_oracle(rs, y0, DT, N_STEPS)

        artifact = _artifact_for(core, rxns)
        # sanity: vocabulary exercised as intended
        archs = sorted(e["archetype"] for e in artifact["reactions"])
        assert archs == ["discrete_chip/1", "migration/1"]
        assert all(e["cantera"] is None and e["kinetics"] for e in artifact["reactions"])

        consumer = ArtifactConsumer(artifact, [_yaml_label(s) for s in core],
                                    P=P_PA, V_poly=V_POLY)
        _, mine = consumer.integrate_euler(y0, T_K, DT, N_STEPS)
        np.testing.assert_allclose(mine, oracle, rtol=1e-9, atol=1e-12)
        # the chip really fired (mu1_A dropped, G grew)
        assert mine[-1, 2] < y0[2]
        assert mine[-1, 8] > y0[8]

    def test_reversible_entry_keq_recipe(self):
        """cantera-null reversible entry: consumer's NASA-based
        kb = kf/Keq, Keq = (P0/RT)^dn * exp(-dG0/RT) matches the oracle's
        rxn.get_equilibrium_constant path (format doc §4 step 1)."""
        sp, core, mask = _two_pool_setup(with_thermo=True)
        rxn = Reaction(reactants=[sp["A"]], products=[sp["B"], sp["G"]],
                       kinetics=_kin1(), reversible=True)
        rxn.polymer_flux_archetype = int(PolymerFluxArchetype.MIGRATION)
        rs = _oracle_system(core, mask, [rxn])
        y0 = rs.y.copy()
        oracle = _euler_oracle(rs, y0, DT, N_STEPS)

        artifact = _artifact_for(core, [rxn])
        nasa = {}
        for lab in ("A", "B", "G"):
            th = sp[lab].thermo
            nasa[lab] = {"Tmid": 1000.0,
                         "low": [float(c) for c in th.polynomials[0].coeffs],
                         "high": [float(c) for c in th.polynomials[1].coeffs]}
        consumer = ArtifactConsumer(artifact, [_yaml_label(s) for s in core],
                                    P=P_PA, V_poly=V_POLY, nasa=nasa)
        _, mine = consumer.integrate_euler(y0, T_K, DT, N_STEPS)
        np.testing.assert_allclose(mine, oracle, rtol=1e-7, atol=1e-12)


class TestVolatileEjectionSufficiency:
    """Ruling round 20, increment 5 (consumer parity): the reference
    consumer implements the volatile_ejection/1 branch, mirroring the
    generating solver's dispatch -- cross-pool a-shifted bundle exchange
    with direction-specific source availability, same-pool chip-style
    signed single-pool write including the a<0 forward-mu2 clamp
    convention, and the same-pool a>0 end-group exhaustion throttle.
    Before this branch existed, VE rows integrated species-only with NO
    moment flux: the gas volatile appeared while the melt kept its mass --
    exactly the fabrication class the r29 conduit rows must not re-open."""

    def test_cross_pool_ve_deck_matches_oracle(self):
        """r29 conduit shape at deck level (pool -> pool_mod + volatile,
        a = +MW(H)/monomer_MW scale): consumer matches oracle on every
        species and every pool moment."""
        a = 0.01008
        sp, core, mask = _two_pool_setup()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["B"], sp["G"]],
                       kinetics=_kin1(), reversible=False)
        rxn.polymer_flux_archetype = int(PolymerFluxArchetype.VOLATILE_EJECTION)
        rxn.polymer_eject_units = a
        rs = _oracle_system(core, mask, [rxn])
        y0 = rs.y.copy()
        oracle = _euler_oracle(rs, y0, DT, N_STEPS)

        artifact = _artifact_for(core, [rxn])
        (entry,) = artifact["reactions"]
        assert entry["archetype"] == "volatile_ejection/1"
        assert entry["params"] == {"eject_units": pytest.approx(a)}

        consumer = ArtifactConsumer(artifact, [_yaml_label(s) for s in core],
                                    P=P_PA, V_poly=V_POLY)
        _, mine = consumer.integrate_euler(y0, T_K, DT, N_STEPS)
        np.testing.assert_allclose(mine, oracle, rtol=1e-9, atol=1e-12)
        # the row really fired: volatile released, chains migrated A -> B
        assert mine[-1, 8] > y0[8]
        assert mine[-1, 1] < y0[1]
        assert mine[-1, 5] > y0[5]

    def test_cross_pool_ve_reversible_directional_scaling(self):
        """Reversible cross-pool VE: the reverse leg debits the DST pool,
        so its rate scales with the dst pool's OWN moments (adjudicated
        Part C in the solver). The consumer must mirror that directional
        source-availability scaling or a reversible deck diverges."""
        sp, core, mask = _two_pool_setup(with_thermo=True)
        rxn = Reaction(reactants=[sp["A"]], products=[sp["B"], sp["G"]],
                       kinetics=_kin1(), reversible=True)
        rxn.polymer_flux_archetype = int(PolymerFluxArchetype.VOLATILE_EJECTION)
        rxn.polymer_eject_units = 1.135
        rs = _oracle_system(core, mask, [rxn])
        y0 = rs.y.copy()
        oracle = _euler_oracle(rs, y0, DT, N_STEPS)

        artifact = _artifact_for(core, [rxn])
        nasa = {}
        for lab in ("A", "B", "G"):
            th = sp[lab].thermo
            nasa[lab] = {"Tmid": 1000.0,
                         "low": [float(c) for c in th.polynomials[0].coeffs],
                         "high": [float(c) for c in th.polynomials[1].coeffs]}
        consumer = ArtifactConsumer(artifact, [_yaml_label(s) for s in core],
                                    P=P_PA, V_poly=V_POLY, nasa=nasa)
        _, mine = consumer.integrate_euler(y0, T_K, DT, N_STEPS)
        np.testing.assert_allclose(mine, oracle, rtol=1e-7, atol=1e-12)

    def test_same_pool_ve_positive_a_matches_oracle(self):
        """Same-pool unzip shape (A -> A + G, a > 0): single-pool signed
        mu1/mu2 write; mu0 untouched; volatile through the species path."""
        sp, core, mask = _two_pool_setup()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["A"], sp["G"]],
                       kinetics=_kin1(), reversible=False)
        rxn.polymer_flux_archetype = int(PolymerFluxArchetype.VOLATILE_EJECTION)
        rxn.polymer_eject_units = 1.135
        rs = _oracle_system(core, mask, [rxn])
        y0 = rs.y.copy()
        oracle = _euler_oracle(rs, y0, DT, N_STEPS)

        artifact = _artifact_for(core, [rxn])
        consumer = ArtifactConsumer(artifact, [_yaml_label(s) for s in core],
                                    P=P_PA, V_poly=V_POLY)
        _, mine = consumer.integrate_euler(y0, T_K, DT, N_STEPS)
        np.testing.assert_allclose(mine, oracle, rtol=1e-9, atol=1e-12)
        assert mine[-1, 2] < y0[2]                    # A mu1 drained
        assert mine[-1, 1] == pytest.approx(y0[1])    # mu0 untouched
        assert mine[-1, 8] > y0[8]                    # volatile released

    def test_same_pool_ve_negative_a_clamp_matches_oracle(self):
        """Signed a < 0 (gas addition onto the same pool, G + A -> A):
        forward mu1 RISES by |a|*ev; the forward mu2 decrement
        2a*E[n] - a^2 is always negative for a < 0, so the `> 0` clamp
        SKIPS the mu2 write -- the documented a<0 convention the consumer
        must reproduce bit-for-bit."""
        kin2 = Arrhenius(A=(2.0, "m^3/(mol*s)"), n=0.0, Ea=(0.0, "kcal/mol"),
                         T0=(298.15, "K"))
        sp, core, mask = _two_pool_setup()
        rxn = Reaction(reactants=[sp["G"], sp["A"]], products=[sp["A"]],
                       kinetics=kin2, reversible=False)
        rxn.polymer_flux_archetype = int(PolymerFluxArchetype.VOLATILE_EJECTION)
        rxn.polymer_eject_units = -1.135
        rs = _oracle_system(core, mask, [rxn])
        y0 = rs.y.copy()
        oracle = _euler_oracle(rs, y0, DT, N_STEPS)

        artifact = _artifact_for(core, [rxn])
        (entry,) = artifact["reactions"]
        assert entry["params"]["eject_units"] == pytest.approx(-1.135)
        consumer = ArtifactConsumer(artifact, [_yaml_label(s) for s in core],
                                    P=P_PA, V_poly=V_POLY)
        _, mine = consumer.integrate_euler(y0, T_K, DT, N_STEPS)
        np.testing.assert_allclose(mine, oracle, rtol=1e-9, atol=1e-12)
        assert mine[-1, 2] > y0[2]     # mu1 rises (mass gained from gas)
        assert mine[-1, 8] < y0[8]     # gas consumed

    def test_same_pool_ve_end_group_throttle_parity(self):
        """Same-pool a>0 END-GROUP VE: the mu0-scaled drain carries the
        min(mu0, mu1/a) exhaustion throttle in the solver's site scaling;
        the consumer must apply the SAME throttle (a=6 with A moments
        (1, 5, 30): mu1/a = 5/6 < mu0 = 1 -- the throttle bites)."""
        sp, core, mask = _two_pool_setup()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["A"], sp["G"]],
                       kinetics=_kin1(), reversible=False)
        rxn.polymer_flux_archetype = int(PolymerFluxArchetype.VOLATILE_EJECTION)
        rxn.polymer_eject_units = 6.0
        rxn.is_end_group_reaction = True
        rs = _oracle_system(core, mask, [rxn])
        y0 = rs.y.copy()
        oracle = _euler_oracle(rs, y0, DT, N_STEPS)

        artifact = _artifact_for(core, [rxn])
        (entry,) = artifact["reactions"]
        assert entry["scaling"] == "mu0"
        consumer = ArtifactConsumer(artifact, [_yaml_label(s) for s in core],
                                    P=P_PA, V_poly=V_POLY)
        _, mine = consumer.integrate_euler(y0, T_K, DT, N_STEPS)
        np.testing.assert_allclose(mine, oracle, rtol=1e-9, atol=1e-12)

    def test_cross_pool_ve_end_group_forward_throttle_parity(self):
        """Item C consumer parity (ruling round 20): cross-pool a>0
        END-GROUP VE -- the forward (src-debiting) leg carries the same
        min(mu0, mu1/a) source-availability throttle as the solver
        (a=6, A moments (1, 5, 30): mu1/a = 5/6 < mu0 = 1, it bites)."""
        sp, core, mask = _two_pool_setup()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["B"], sp["G"]],
                       kinetics=_kin1(), reversible=False)
        rxn.polymer_flux_archetype = int(PolymerFluxArchetype.VOLATILE_EJECTION)
        rxn.polymer_eject_units = 6.0
        rxn.is_end_group_reaction = True
        rs = _oracle_system(core, mask, [rxn])
        y0 = rs.y.copy()
        oracle = _euler_oracle(rs, y0, DT, N_STEPS)

        artifact = _artifact_for(core, [rxn])
        consumer = ArtifactConsumer(artifact, [_yaml_label(s) for s in core],
                                    P=P_PA, V_poly=V_POLY)
        _, mine = consumer.integrate_euler(y0, T_K, DT, N_STEPS)
        np.testing.assert_allclose(mine, oracle, rtol=1e-9, atol=1e-12)

    def test_cross_pool_ve_end_group_reverse_throttle_parity(self):
        """Item C mirror direction, consumer parity: a<0 END-GROUP
        cross-pool VE -- the reverse leg sheds |a|/event from the dst pool
        it debits, so its availability carries min(mu0_dst, mu1_dst/|a|)
        (B moments (2, 4, 10), |a|=6: 4/6 < 2, it bites)."""
        sp, core, mask = _two_pool_setup(with_thermo=True)
        rxn = Reaction(reactants=[sp["A"]], products=[sp["B"], sp["G"]],
                       kinetics=_kin1(), reversible=True)
        rxn.polymer_flux_archetype = int(PolymerFluxArchetype.VOLATILE_EJECTION)
        rxn.polymer_eject_units = -6.0
        rxn.is_end_group_reaction = True
        rs = _oracle_system(core, mask, [rxn])
        y0 = rs.y.copy()
        oracle = _euler_oracle(rs, y0, DT, N_STEPS)

        artifact = _artifact_for(core, [rxn])
        nasa = {}
        for lab in ("A", "B", "G"):
            th = sp[lab].thermo
            nasa[lab] = {"Tmid": 1000.0,
                         "low": [float(c) for c in th.polynomials[0].coeffs],
                         "high": [float(c) for c in th.polynomials[1].coeffs]}
        consumer = ArtifactConsumer(artifact, [_yaml_label(s) for s in core],
                                    P=P_PA, V_poly=V_POLY, nasa=nasa)
        _, mine = consumer.integrate_euler(y0, T_K, DT, N_STEPS)
        np.testing.assert_allclose(mine, oracle, rtol=1e-7, atol=1e-12)

    def test_rejects_ve_row_without_eject_units(self):
        """A volatile_ejection/1 entry without the exact
        params = {'eject_units': ...} shape must fail LOUDLY at
        construction -- defaulting to 0.0 would launder the atom-transfer
        debit away (fabricated mass)."""
        sp, core, mask = _two_pool_setup()
        rxn = Reaction(reactants=[sp["A"]], products=[sp["B"], sp["G"]],
                       kinetics=_kin1(), reversible=False)
        rxn.polymer_flux_archetype = int(PolymerFluxArchetype.VOLATILE_EJECTION)
        rxn.polymer_eject_units = 1.0
        artifact = _artifact_for(core, [rxn])
        del artifact["reactions"][0]["params"]
        with pytest.raises(ValueError, match="eject_units"):
            ArtifactConsumer(artifact, [_yaml_label(s) for s in core],
                             P=P_PA, V_poly=V_POLY)


class TestSmilesBranchPoolLabel:
    def test_legacy_mu1_deck_with_smiles_branch_pool_label(self):
        """PP run-5 erratum parity: base labels strip ONLY a trailing
        '(<int>)' RMG index (rmgpy.polymer.strip_rmg_index_suffix). This
        consumer's _base() used to truncate at the FIRST '(' -- correct for
        'PS(2)' but SMILES-mangling for a pool labelled with branching
        parentheses: 'C[CH]CC(C)C(2)' -> 'C[CH]CC', so its legacy_mu1 proxy
        rows dereferenced a nonexistent pool (KeyError at rhs time).
        End-to-end: consumer must match the oracle on a deck whose pool
        label IS branched SMILES."""
        smiles = "C[CH]CC(C)C"
        inert = _spc("N#N", "N2")
        proxy = _spc(smiles, smiles, index=2)  # chem.yaml 'C[CH]CC(C)C(2)'
        proxy.is_polymer_proxy = True
        mus = [_mu(f"{smiles}_mu{k}") for k in range(3)]
        gas = _spc("[CH3]", "G", index=9)
        core = [inert, proxy] + mus + [gas]
        mask = np.array([True, False, False, False, False, True], dtype=bool)

        pool_cfg = PolymerPoolConfig(
            label=smiles, xs=2, explicit_dp_to_species_index={},
            mu_indices=(2, 3, 4), monomer_poly_index=None,
            k_scission=0.0, k_unzip=0.0, tail_kinetics=None)
        # deliberately UNSTAMPED (no polymer_flux_archetype): the artifact
        # row lands on 'legacy_mu1/1', the branch whose proxy_r/p_pools go
        # through the consumer's _base().
        rxn = Reaction(reactants=[proxy], products=[proxy, gas],
                       kinetics=_kin1(), reversible=False)

        rs = HybridPolymerSystem(
            T=T_K, P=P_PA, initial_mole_fractions={inert: 1.0}, V_poly=V_POLY,
            polymer_pools=[pool_cfg], mass_transfer=[],
            gas_species_mask=mask.copy(), constant_gas_volume=False,
            initial_polymer_moments={smiles: (1.0, 5.0, 30.0)},
            termination=[],
            # r71 FIX 4 escape: this fixture's row is deliberately unstamped
            # to exercise the legacy_mu1/1 artifact branch (direct-test
            # posture; production default False hard-fails).
            allow_unstamped_proxy_rows=True)
        rs.initialize_model(core, [rxn], [], [])
        y0 = rs.y.copy()
        oracle = _euler_oracle(rs, y0, DT, N_STEPS)

        registry_pool = Polymer(label=smiles, monomer="[CH2][CH2]",
                                end_groups=["[H]", "[H]"], cutoff=3,
                                moments=[1.0, 5.0, 30.0], initial_mass=0.0)
        artifact = build_polymer_moments_artifact(
            [registry_pool], core_species=core, core_reactions=[rxn],
            configured_pool_labels=[smiles],
            condensed_species=core[1:5], cantera_index_map={})
        artifact = json.loads(json.dumps(artifact))
        entry, = artifact["reactions"]
        assert entry["archetype"] == "legacy_mu1/1"
        assert entry["proxy_reactants"] == [f"{smiles}(2)"]

        consumer = ArtifactConsumer(artifact, [_yaml_label(s) for s in core],
                                    P=P_PA, V_poly=V_POLY)
        _, mine = consumer.integrate_euler(y0, T_K, DT, N_STEPS)
        np.testing.assert_allclose(mine, oracle, rtol=1e-9, atol=1e-12)
        assert mine[-1, 5] > y0[5]  # the gas product actually accumulated


class TestMassTransferCrossCheck:
    def test_evaporation_both_sides_nonzero(self):
        """Spec §9 test 2b: nonzero kLa with BOTH C_poly and C_gas nonzero, on
        top of the two-pool reaction deck — the consumer's J = kLa(Cp - K*Cg)
        is validated against the solver's, end to end."""
        sp, core, mask = _two_pool_setup()
        # dissolved species D (condensed, not pool-mapped) <-> gas G
        d = _spc("CC=O", "D", index=11)
        core = core + [d]
        mask = np.append(mask, False)  # D is condensed

        mig = Reaction(reactants=[sp["A"]], products=[sp["B"]],
                       kinetics=_kin1(), reversible=False)
        mig.polymer_flux_archetype = int(PolymerFluxArchetype.MIGRATION)
        rxns = [mig]

        mt = MassTransferConfig(gas_index=8, poly_index=9, K=2.0, kLa=5.0)
        rs = HybridPolymerSystem(
            T=T_K, P=P_PA,
            initial_mole_fractions={core[8]: 0.05, d: 0.0},
            V_poly=V_POLY,
            polymer_pools=_pools_ab(), mass_transfer=[mt],
            gas_species_mask=mask.copy(), constant_gas_volume=False,
            initial_polymer_moments={"A": (1.0, 5.0, 30.0), "B": (2.0, 4.0, 10.0)},
            termination=[])
        rs.initialize_model(core, rxns, [], [])
        y0 = rs.y.copy()
        y0[9] = 0.02  # dissolved D moles: condensed side nonzero too
        oracle = _euler_oracle(rs, y0, DT, N_STEPS)

        artifact = build_polymer_moments_artifact(
            _registry_ab(), core_species=core, core_reactions=rxns,
            configured_pool_labels=["A", "B"],
            condensed_species=core[:8] + [d], cantera_index_map={})
        artifact = json.loads(json.dumps(artifact))
        consumer = ArtifactConsumer(
            artifact, [_yaml_label(s) for s in core], P=P_PA, V_poly=V_POLY,
            mass_transfer=[{"gas": "G(9)", "poly": "D(11)", "K": 2.0, "kLa": 5.0}])
        _, mine = consumer.integrate_euler(y0, T_K, DT, N_STEPS)
        np.testing.assert_allclose(mine, oracle, rtol=1e-9, atol=1e-12)
        # the transfer really moved mass and conserved the pair total
        assert mine[-1, 9] != pytest.approx(y0[9])
        np.testing.assert_allclose(mine[:, 8] + mine[:, 9],
                                   y0[8] + y0[9], rtol=1e-9)


# ---------------------------------------------------------------------------
# (c) refused-row + homolysis-kernel supersession guards (schema 2.4 / 2.6,
#     Stage 2 of the adjudicated round-66/67 arc, ruling (c)): this consumer
#     must never convert a refused row's kinetics into flux, must validate
#     the CLOSED refused_reason vocabulary strictly, and must hard-reject a
#     pool-level homolysis_initiation block it does not implement (silently
#     integrating without the kernel produces a flat/false trajectory --
#     the exact failure class the sidecar schema exists to kill).
# ---------------------------------------------------------------------------


class TestRefusedRowAndHomolysisGuards:

    def _refused_artifact(self, mark_refused=True):
        """One pool-mapped row (poly(2) + H -> R + H2), optionally stamped
        refused/conduit-deferred -- the production stamp-but-keep shape."""
        n2 = _spc("N#N", "N2")
        proxy = _spc("CCCC", "poly", index=2)
        rad = _spc("[CH2]CCC", "R", index=3)
        h = _spc("[H]", "H", index=4)
        h2 = _spc("[H][H]", "H2", index=5)
        proxy.is_polymer_proxy = True
        mus = [_mu("poly_mu0"), _mu("poly_mu1"), _mu("poly_mu2")]
        core = [n2, proxy, rad, h, h2] + mus
        rxn = Reaction(
            reactants=[proxy, h], products=[rad, h2],
            kinetics=Arrhenius(A=(2.0, "m^3/(mol*s)"), n=0.0,
                               Ea=(0.0, "J/mol"), T0=(1.0, "K")),
            reversible=False)
        rxn.polymer_flux_archetype = int(PolymerFluxArchetype.UNRESOLVED)
        if mark_refused:
            rxn.polymer_refused = True
            rxn.polymer_refused_accumulating = False
        pool = Polymer(label="poly", monomer="[CH2][CH2]",
                       end_groups=["[H]", "[H]"], cutoff=3,
                       moments=[1.0, 5.0, 30.0], initial_mass=0.0,
                       k_scission=0.0, k_unzip=0.0)
        artifact = build_polymer_moments_artifact(
            [pool], core_species=core, core_reactions=[rxn],
            configured_pool_labels=["poly"],
            condensed_species=[proxy, rad] + mus,
            cantera_index_map={})
        return json.loads(json.dumps(artifact)), core

    def _consumer_and_y0(self, artifact, core):
        order = [_yaml_label(s) for s in core]
        consumer = ArtifactConsumer(artifact, order, P=P_PA, V_poly=V_POLY)
        i = {lab: k for k, lab in enumerate(order)}
        y0 = np.zeros(len(order))
        y0[i["N2"]] = 1.0
        y0[i["H(4)"]] = 2.0
        y0[i["poly_mu0"]], y0[i["poly_mu1"]], y0[i["poly_mu2"]] = 1.0, 5.0, 30.0
        return consumer, y0, i

    def test_refused_row_contributes_exactly_zero_flux(self):
        """RED pin (ruling (c)): the refused row is the ONLY reaction and
        both channels are off, so an honest consumer's RHS is exactly zero
        everywhere -- converting the row's kinetics to flux is the pre-2.4
        over-integration erratum."""
        artifact, core = self._refused_artifact(mark_refused=True)
        assert artifact["schema_version"] == "2.4"
        consumer, y0, _ = self._consumer_and_y0(artifact, core)
        dn = consumer.rhs(y0, 800.0)
        assert np.all(dn == 0.0), f"refused row converted to flux: {dn}"

    def test_unmarked_row_still_integrates(self):
        """Control: the SAME chemistry without the refused stamp integrates
        (proves the zero above is suppression, not a dead fixture)."""
        artifact, core = self._refused_artifact(mark_refused=False)
        consumer, y0, i = self._consumer_and_y0(artifact, core)
        dn = consumer.rhs(y0, 800.0)
        assert dn[i["R(3)"]] > 0.0
        assert dn[i["poly_mu1"]] < 0.0

    def test_qssa_unassessable_reason_accepted_and_zero_flux(self):
        """Round-20 increment 7: 'qssa-unassessable' joins the CLOSED
        refused_reason vocabulary (the rebuild census's spelling of an
        accumulating refusal whose radical/consumers were not visible --
        missing evidence is never 'qssa-invalid'). Suppression semantics
        identical: whole-row zero flux."""
        artifact, core = self._refused_artifact(mark_refused=True)
        (row,) = [e for e in artifact["reactions"] if e.get("refused")]
        row["refused_reason"] = "qssa-unassessable"
        consumer, y0, _ = self._consumer_and_y0(artifact, core)
        dn = consumer.rhs(y0, 800.0)
        assert np.all(dn == 0.0)

    def test_rejects_unknown_refused_reason(self):
        """The refused_reason vocabulary is CLOSED (format doc §12):
        reject at construction, never adapt."""
        artifact, core = self._refused_artifact(mark_refused=True)
        (row,) = [e for e in artifact["reactions"] if e.get("refused")]
        row["refused_reason"] = "totally-new-reason"
        order = [_yaml_label(s) for s in core]
        with pytest.raises(ValueError, match=r"refused_reason"):
            ArtifactConsumer(artifact, order, P=P_PA, V_poly=V_POLY)

    def test_rejects_refused_false_shape(self):
        """The emitter writes refused: true + reason or NOTHING (absent,
        not false) -- present-but-false is malformed."""
        artifact, core = self._refused_artifact(mark_refused=True)
        (row,) = [e for e in artifact["reactions"] if e.get("refused")]
        row["refused"] = False
        order = [_yaml_label(s) for s in core]
        with pytest.raises(ValueError, match=r"refused"):
            ArtifactConsumer(artifact, order, P=P_PA, V_poly=V_POLY)

    def test_rejects_homolysis_initiation_block_loudly(self):
        """RED pin (ruling (c) supersession contract): this consumer does
        not implement the radical-homolysis kernel, so a pool carrying the
        schema-2.6 homolysis_initiation block must fail at construction --
        the silent path integrates a melt with no initiation flux (flat/
        false trajectory). Mirrors the radical_qssa_unzip channel guard."""
        pool = Polymer(label="PP", monomer="[CH2][CH](C)",
                       end_groups=["[H]", "[H]"], cutoff=3,
                       moments=[1.0, 50.0, 3000.0], initial_mass=0.0,
                       k_homolysis={"A": 1.0e13, "n": 0.5, "Ea": 1.2e5})
        prim, sec = pool.generate_end_radical_daughters()
        core = [_spc("N#N", "N2")]
        for base in ("PP", prim.label, sec.label):
            core += [_mu(f"{base}_mu{k}") for k in range(3)]
        artifact = build_polymer_moments_artifact(
            [pool, prim, sec], core_species=core, core_reactions=[],
            configured_pool_labels=["PP", prim.label, sec.label],
            condensed_species=core[1:],
            cantera_index_map={})
        artifact = json.loads(json.dumps(artifact))
        assert artifact["schema_version"] == "2.6"
        order = [_yaml_label(s) for s in core]
        with pytest.raises(ValueError, match=r"PP.*homolysis_initiation"):
            ArtifactConsumer(artifact, order, P=P_PA, V_poly=V_POLY)

    @staticmethod
    def _side_group_artifact():
        """A real emitter artifact carrying the schema-3.0 SGH kernel-v2
        (side_group_homolysis/2) block; v2 spawns no feature pool."""
        br = _spc("[Br]", "Br", index=7)
        pool = Polymer(label="PVBr", monomer="[CH2][CH]Br",
                       end_groups=["[H]", "[H]"], cutoff=3,
                       moments=[1.0, 50.0, 3000.0], initial_mass=0.0,
                       side_group_homolysis=[dict(
                           label="aliphatic_C-Br", A=1.0e13, n=0.5,
                           Ea=1.2e5, site_selector="aliphatic",
                           sites_per_unit=1.0, gas_product="[Br]")])
        pool.side_group_gas_species = [br]
        core = ([_spc("N#N", "N2"), br]
                + [_mu(f"PVBr_mu{k}") for k in range(3)])
        artifact = build_polymer_moments_artifact(
            [pool], core_species=core, core_reactions=[],
            configured_pool_labels=["PVBr"],
            condensed_species=core[2:],
            cantera_index_map={})
        return json.loads(json.dumps(artifact)), core

    def test_rejects_side_group_homolysis_block_loudly(self):
        """RED pin (FR1-K2 supersession contract): this reference consumer
        does not implement the side-group homolysis kernel, so a pool
        carrying the schema-3.0 SGH kernel-v2 (side_group_homolysis/2)
        block must fail at construction -- the silent path integrates a
        melt with no X-loss flux."""
        artifact, core = self._side_group_artifact()
        assert artifact["schema_version"] == "3.0"
        order = [_yaml_label(s) for s in core]
        with pytest.raises(ValueError,
                           match=r"PVBr.*side_group_homolysis"):
            ArtifactConsumer(artifact, order, P=P_PA, V_poly=V_POLY)

    def test_accepts_chain_mass_defect_pool_with_normative_closure(self):
        """P1-2 flip of the former RED pin: this consumer now IMPLEMENTS
        the normative defect-aware condensed-mass closure
        (mu1*monomer_mw_g_mol - mu0*chain_mass_defect_g_mol) and the
        defect-aware VE moment booking, so a pool carrying
        chain_mass_defect_g_mol (additive optional field, e.g. an H-loss
        _mod daughter) constructs cleanly and its mass closure applies the
        per-chain defect. A malformed defect still fails loudly."""
        artifact, core = self._side_group_artifact()
        for p in artifact["pools"]:
            p.pop("side_group_homolysis", None)
        pool = next(p for p in artifact["pools"] if p["label"] == "PVBr")
        pool["chain_mass_defect_g_mol"] = 79.904
        order = [_yaml_label(s) for s in core]
        consumer = ArtifactConsumer(artifact, order, P=P_PA, V_poly=V_POLY)
        assert consumer.pools["PVBr"]["defect"] == 79.904
        y = np.zeros(len(order))
        i0, i1, _ = consumer.pools["PVBr"]["mu"]
        y[i0] = 2.0
        y[i1] = 10.0
        mw = consumer.pools["PVBr"]["mw"]
        assert mw > 0.0
        assert np.isclose(consumer.condensed_mass_g(y),
                          10.0 * mw - 2.0 * 79.904, rtol=1e-12)
        # Malformed defect (negative / non-finite / boolean) stays a loud
        # construction error -- never a silent 0.0 default.
        for bad in (-1.0, float("nan"), True):
            pool["chain_mass_defect_g_mol"] = bad
            with pytest.raises(ValueError,
                               match=r"chain_mass_defect_g_mol"):
                ArtifactConsumer(artifact, order, P=P_PA, V_poly=V_POLY)

    def test_rejects_end_radical_depropagation_block_loudly(self):
        """RED pin (schema 2.8, r74 SS2 kernel): this consumer does not
        implement the end-radical depropagation kernel, so a daughter pool
        carrying the end_radical_depropagation block must fail at
        construction -- the silent path integrates radical-end pools with
        no consumption channel (the run-6 no-outlet wall) and drops the
        gas monomer source (un-conserved mass)."""
        c3h6 = _spc("C=CC", "C3H6", index=5)
        pool = Polymer(label="PP", monomer="[CH2][CH](C)",
                       end_groups=["[H]", "[H]"], cutoff=3,
                       moments=[1.0, 50.0, 3000.0], initial_mass=0.0,
                       k_homolysis={"A": 1.0e13, "n": 0.5, "Ea": 1.2e5},
                       k_depropagation={"A": 9.4e14, "n": 0.0,
                                        "Ea": 117152.0})
        pool.monomer_product_species = c3h6
        prim, sec = pool.generate_end_radical_daughters()
        core = [_spc("N#N", "N2"), c3h6]
        for base in ("PP", prim.label, sec.label):
            core += [_mu(f"{base}_mu{k}") for k in range(3)]
        artifact = build_polymer_moments_artifact(
            [pool, prim, sec], core_species=core, core_reactions=[],
            configured_pool_labels=["PP", prim.label, sec.label],
            condensed_species=core[2:],
            cantera_index_map={})
        artifact = json.loads(json.dumps(artifact))
        assert artifact["schema_version"] == "2.8"
        # Isolate THIS guard from the parent's sibling 2.6 rejection
        # (which fires first in pool order) -- the chain_mass_defect
        # precedent: strip the other kernel's block and pin the deprop
        # rejection specifically.
        for p in artifact["pools"]:
            p.pop("homolysis_initiation", None)
        assert any("end_radical_depropagation" in p
                   for p in artifact["pools"])
        order = [_yaml_label(s) for s in core]
        with pytest.raises(
                ValueError,
                match=r"PP_rad_.*end.*end_radical_depropagation"):
            ArtifactConsumer(artifact, order, P=P_PA, V_poly=V_POLY)


# ---------------------------------------------------------------------------
# M18.3 T7, COMMIT 2 (the fix): the commit-1 fabrication proof INVERTED into
# raise-assertions, + envelope-gate tests, + the conduit bundle (DESIGN
# §2.2) with the trajectory twin of the runner's T5 replay (which the
# runner itself cannot integrate until the M18.4 solver dispatch arm --
# see polymerMomentsRunnerTest.py's M18.4 refusal pin).
# ---------------------------------------------------------------------------

class TestSilentMassFabricationProof:
    """T7 -- commit 1 of the ratified test-first rider PROVED the silent
    fabrication with assertions written to the buggy behavior (2.0 mol of
    gas fabricated over the 0.2 s window, ~30 g, with every pool moment
    bit-identical and no raise/warning). This commit INVERTS them: the
    same fixture must now refuse at construction, at BOTH gates."""

    def _unknown_archetype_artifact(self, core):
        rxn = Reaction(reactants=[core[0]], products=[core[4], core[8]],
                       kinetics=_kin1(), reversible=False)
        rxn.polymer_flux_archetype = int(PolymerFluxArchetype.MIGRATION)
        artifact = _artifact_for(core, [rxn])
        (row,) = artifact["reactions"]
        row["archetype"] = "mystery_channel/1"     # unknown to every loader
        artifact["schema_version"] = "9.9"         # unknown envelope
        return artifact

    def test_unknown_version_now_refused_at_the_envelope(self):
        sp, core, mask = _two_pool_setup()
        artifact = self._unknown_archetype_artifact(core)
        with pytest.raises(ValueError, match=r"schema_version '9\.9'"):
            ArtifactConsumer(artifact, [_yaml_label(s) for s in core],
                             P=P_PA, V_poly=V_POLY)

    def test_unknown_archetype_now_refused_at_the_closed_set(self):
        """Even under an ACCEPTED stamp, the unknown archetype refuses at
        the closed-set dispatch guard (the step-5 fabrication is
        unreachable)."""
        sp, core, mask = _two_pool_setup()
        artifact = self._unknown_archetype_artifact(core)
        artifact["schema_version"] = "2.4"         # accepted envelope
        with pytest.raises(ValueError, match=r"mystery_channel/1.*CLOSED"):
            ArtifactConsumer(artifact, [_yaml_label(s) for s in core],
                             P=P_PA, V_poly=V_POLY)


class TestConsumerEnvelopeGate:
    """T7 commit 2: the consumer's accepted set is IDENTICAL to the
    reference runner's (2.0..2.8 + 3.0..3.1) so the two reference
    consumers agree on the refusal surface."""

    def _plain_artifact(self, core):
        return _artifact_for(core, [])

    @pytest.mark.parametrize("ver", ["2.0", "2.4", "2.8", "3.0", "3.1"])
    def test_known_versions_accepted(self, ver):
        sp, core, mask = _two_pool_setup()
        artifact = self._plain_artifact(core)
        artifact["schema_version"] = ver
        ArtifactConsumer(artifact, [_yaml_label(s) for s in core],
                         P=P_PA, V_poly=V_POLY)

    @pytest.mark.parametrize("ver", ["2.9", "2.10", "3.2", "4.0", "1.0",
                                     "garbage", "3", "3.1.1"])
    def test_unknown_versions_refused_loudly(self, ver):
        sp, core, mask = _two_pool_setup()
        artifact = self._plain_artifact(core)
        artifact["schema_version"] = ver
        with pytest.raises(ValueError, match=r"not implemented"):
            ArtifactConsumer(artifact, [_yaml_label(s) for s in core],
                             P=P_PA, V_poly=V_POLY)


def _conduit_setup(chain_units=None, rate_a=0.05):
    """Consumer-world conduit deck: CHAIN => poly + CH2O, emitted by the
    REAL serializer (cantera-indexed so the row carries kinetics +
    cantera), balanced against the pool's representative molecule.
    ``chain_units`` pins u exactly by overriding the pool monomer MW
    (T6 equality boundary: 1.0)."""
    pf = Polymer(label="poly", monomer="[CH2]c1ccc(cc1)C([CH2])O",
                 end_groups=["[H]", "[H]"], cutoff=3,
                 moments=[1.0, 5.0, 30.0], initial_mass=0.0)
    pf.index = 2
    proxy_smiles = pf.molecule[0].to_smiles()
    assert proxy_smiles.startswith("C")
    chain = _spc("OCC" + proxy_smiles[1:], "CHAIN", index=3)
    gas = _spc("C=O", "CH2O", index=4)
    mus = [_mu("poly_mu0"), _mu("poly_mu1"), _mu("poly_mu2")]
    core = [pf] + mus + [chain, gas]

    chain_mw = chain.molecule[0].get_molecular_weight() * 1000.0
    gas_mw = gas.molecule[0].get_molecular_weight() * 1000.0
    if chain_units is not None:
        pf.monomer_mw_g_mol = (chain_mw - gas_mw) / float(chain_units)
    m = pf.monomer_mw_g_mol
    u = (chain_mw - gas_mw) / m

    rxn = Reaction(reactants=[chain], products=[pf, gas],
                   kinetics=Arrhenius(A=(rate_a, "1/s"), n=0.0,
                                      Ea=(0.0, "J/mol"), T0=(1.0, "K")),
                   reversible=False)
    rxn.polymer_flux_archetype = int(
        PolymerFluxArchetype.MOMENT_CREDIT_CONDUIT)
    rxn.polymer_conduit_dst_pool = "poly"
    rxn.polymer_conduit_params = {
        "admission_direction": "chain_to_pool",
        "chain_units": u,
        "gas_products": [{"species": "CH2O(4)", "stoich": 1,
                          "mw_g_mol": gas_mw}],
        "gas_units": gas_mw / m,
        "candidate_key": "CH2O(4)+poly(2)<>CHAIN(3)",
        "candidate_key_note": "run-scoped provenance only",
    }
    artifact = build_polymer_moments_artifact(
        [pf], core_species=core, core_reactions=[rxn],
        configured_pool_labels=["poly"],
        condensed_species=[pf, chain] + mus,
        cantera_index_map={id(rxn): [0]})
    artifact = json.loads(json.dumps(artifact))
    order = [_yaml_label(s) for s in core]
    return artifact, order, u


class TestConduitBundleConsumer:
    """T7 commit 2 / T5-twin / T6 replay leg: the DESIGN §2.2 bundle law
    in the numpy reference consumer (the runner refuses conduit replays
    until M18.4 -- this is the normative trajectory verification)."""

    def _y0(self, order, chain_moles=1.0e-3):
        y0 = np.zeros(len(order))
        y0[order.index("poly_mu0")] = 1.0
        y0[order.index("poly_mu1")] = 5.0
        y0[order.index("poly_mu2")] = 30.0
        y0[order.index("CHAIN(3)")] = chain_moles
        return y0

    def test_artifact_stamps_3_1_and_loads(self):
        artifact, order, u = _conduit_setup()
        assert artifact["schema_version"] == "3.1"
        ArtifactConsumer(artifact, order, P=P_PA, V_poly=V_POLY)

    def test_conduit_row_under_2_4_stamp_refused(self):
        """T1(d), numpy half: conduit vocabulary under an accepted-but-old
        stamp is malformed."""
        artifact, order, u = _conduit_setup()
        artifact["schema_version"] = "2.4"
        with pytest.raises(ValueError, match=r"3\.1"):
            ArtifactConsumer(artifact, order, P=P_PA, V_poly=V_POLY)

    def test_conduit_row_under_3_0_stamp_refused(self):
        artifact, order, u = _conduit_setup()
        artifact["schema_version"] = "3.0"
        with pytest.raises(ValueError, match=r"3\.1"):
            ArtifactConsumer(artifact, order, P=P_PA, V_poly=V_POLY)

    @pytest.mark.parametrize("mutate, pattern", [
        (lambda row: row.__setitem__("cantera", None), r"cantera: null"),
        (lambda row: row["kinetics"].__setitem__("reversible", True),
         r"reversible"),
        (lambda row: row.__setitem__("src_pool", "poly"), r"src_pool"),
        (lambda row: row.__setitem__("dst_pool", "ghost"), r"destination"),
        (lambda row: row.__setitem__("proxy_reactants", ["poly(2)"]),
         r"proxy_reactants"),
        (lambda row: row.__setitem__("proxy_products", []),
         r"proxy_products"),
        (lambda row: row["params"].__setitem__("chain_units", 0.5),
         r"landing cone|>= 1\.0"),
        (lambda row: row["params"].__setitem__("admission_direction",
                                               "pool_to_chain"),
         r"chain_to_pool"),
        (lambda row: row["params"]["gas_products"][0].__setitem__(
            "stoich", 2), r"EXACTLY ONE"),
        (lambda row: row["params"].__setitem__("extra", 1), r"closed"),
        # r42 P1-5: candidate_key semantic pin (mirror of the runner's)
        (lambda row: row["params"].__setitem__("candidate_key", ""),
         r"candidate_key.*non-empty"),
        (lambda row: row["params"].__setitem__("candidate_key",
                                               "TOTALLY+BOGUS<>KEY"),
         r"candidate_key.*does not recompute"),
    ])
    def test_conduit_reject_rules_mirror(self, mutate, pattern):
        artifact, order, u = _conduit_setup()
        (row,) = [e for e in artifact["reactions"]
                  if e["archetype"] == "moment_credit_conduit/1"]
        mutate(row)
        with pytest.raises(ValueError, match=pattern):
            ArtifactConsumer(artifact, order, P=P_PA, V_poly=V_POLY)

    def test_bundle_trajectory_matches_hand_computed_credit(self):
        """T5-twin: CHAIN (n0 = 1e-3 mol) decays first-order at
        k = 0.05 s^-1 (unimolecular, condensed V_rxn = V_poly); every
        consumed mole credits the pool (1, u, u^2) and releases one CH2O:
        mu_k(t) = mu_k(0) + u^k * n0 * (1 - e^-kt)."""
        artifact, order, u = _conduit_setup()
        consumer = ArtifactConsumer(artifact, order, P=P_PA, V_poly=V_POLY)
        n0, k = 1.0e-3, 0.05
        y0 = self._y0(order, n0)
        _, traj = consumer.integrate_euler(y0, T_K, DT, N_STEPS)
        y = traj[-1]
        t_end = DT * N_STEPS
        converted = n0 * (1.0 - np.exp(-k * t_end))
        i_mu = [order.index(f"poly_mu{j}") for j in range(3)]
        # forward-Euler discretization error is ~k*dt/2 = 2.5e-6 relative;
        # 1e-5 pins the LAW (a wrong u or a missing credit errs at O(1)).
        assert y[i_mu[0]] - 1.0 == pytest.approx(converted, rel=1e-5)
        assert y[i_mu[1]] - 5.0 == pytest.approx(u * converted, rel=1e-5)
        assert y[i_mu[2]] - 30.0 == pytest.approx(u * u * converted,
                                                  rel=1e-5)
        # species side: chain consumed, one gas released per event
        assert y[order.index("CHAIN(3)")] == pytest.approx(n0 - converted,
                                                           rel=1e-5)
        assert y[order.index("CH2O(4)")] == pytest.approx(converted,
                                                          rel=1e-5)
        # mass audit (§4.5): condensed pool mass gained == chain mass lost
        # minus gas mass released, per element balance
        m = next(p for p in artifact["pools"]
                 if p["label"] == "poly")["monomer_mw_g_mol"]
        gained_g = (y[i_mu[1]] - 5.0) * m
        chain_mw = 434.57                       # C28H34O4
        gas_mw = 30.026
        assert gained_g == pytest.approx(converted * (chain_mw - gas_mw),
                                         rel=1e-3)

    def test_forward_only_credit_saturates(self):
        """r = rf with NO reverse leg: at t >> 1/k the credit saturates at
        exactly n0 * (1, u, u^2) -- a reverse contribution would pull the
        moments back."""
        artifact, order, u = _conduit_setup(rate_a=50.0)
        consumer = ArtifactConsumer(artifact, order, P=P_PA, V_poly=V_POLY)
        n0 = 1.0e-3
        y0 = self._y0(order, n0)
        _, traj = consumer.integrate_euler(y0, T_K, DT, N_STEPS)
        y = traj[-1]
        i_mu = [order.index(f"poly_mu{j}") for j in range(3)]
        assert y[i_mu[0]] - 1.0 == pytest.approx(n0, rel=1e-3)
        assert y[i_mu[1]] - 5.0 == pytest.approx(u * n0, rel=1e-3)
        assert y[order.index("CHAIN(3)")] == pytest.approx(0.0, abs=1e-7)

    def test_equality_boundary_replay_reports_clean(self):
        """T6 (replay leg, rider; OQ-4): chain_units == 1.0 credits a
        POINT MASS ON the cone (per-event mu1 - mu0 surplus exactly zero).
        The trajectory must integrate cleanly: equal credits across
        mu0/mu1/mu2, whole-pool Q10 = mu1 - mu0 unchanged, mu3 closure
        finite, no negative moments, no spurious refusal. REPORT: the
        closed >= semantics survive the boundary -- no stiffness or
        corruption observed (admission-side leg agrees,
        polymerConduitTest.TestLandingConeEqualityBoundary)."""
        artifact, order, u = _conduit_setup(chain_units=1.0)
        assert u == pytest.approx(1.0, abs=1e-12)
        consumer = ArtifactConsumer(artifact, order, P=P_PA, V_poly=V_POLY)
        n0, k = 1.0e-3, 0.05
        y0 = self._y0(order, n0)
        _, traj = consumer.integrate_euler(y0, T_K, DT, N_STEPS)
        y = traj[-1]
        t_end = DT * N_STEPS
        converted = n0 * (1.0 - np.exp(-k * t_end))
        i_mu = [order.index(f"poly_mu{j}") for j in range(3)]
        d = [y[i_mu[j]] - (1.0, 5.0, 30.0)[j] for j in range(3)]
        assert d[0] == pytest.approx(converted, rel=1e-5)
        assert d[1] == pytest.approx(converted, rel=1e-5)   # u == 1
        assert d[2] == pytest.approx(converted, rel=1e-5)   # u^2 == 1
        # The three PER-STEP increments are bitwise equal at u == 1.0;
        # the accumulated deltas differ only by float-addition rounding
        # against the different base magnitudes (1.0 / 5.0 / 30.0) --
        # benign, and part of the T6 boundary REPORT (no stiffness).
        assert d[1] == pytest.approx(d[0], rel=1e-8)
        assert d[2] == pytest.approx(d[0], rel=1e-8)
        # whole-pool cone margin unchanged by the ON-cone credit
        assert (y[i_mu[1]] - y[i_mu[0]]) == pytest.approx(4.0, rel=1e-9)
        # trajectory clean: no negative moments anywhere
        assert np.all(traj[:, i_mu] >= 0.0)
        # mu3 closure stays finite on the final state
        from numpy_moments_consumer import safe_mu3
        mu3 = safe_mu3(y[i_mu[0]], y[i_mu[1]], y[i_mu[2]])
        assert np.isfinite(mu3)
