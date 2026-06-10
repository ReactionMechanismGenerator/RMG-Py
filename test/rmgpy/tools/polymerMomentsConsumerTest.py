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
        core = [inert, _mu("poly_mu0"), _mu("poly_mu1"), _mu("poly_mu2")]
        mask = np.array([True, False, False, False], dtype=bool)
        pool_cfg = PolymerPoolConfig(label="poly", xs=2,
                                     explicit_dp_to_species_index={},
                                     mu_indices=(1, 2, 3), monomer_poly_index=None,
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
            condensed_species=core[1:4], cantera_index_map={})
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
