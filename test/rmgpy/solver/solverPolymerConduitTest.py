#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2023 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

"""TDD tripwire for the conduit-admission moment-credit archetype (arch 7,
``PolymerFluxArchetype.MOMENT_CREDIT_CONDUIT``).

The object side of conduit admission is complete-but-dead behind
``rmgpy.polymer_conduit.CONDUIT_ADMISSION_ENABLED = False``: an admitted row
is stamped ``polymer_flux_archetype = 7`` with ``polymer_conduit_dst_pool`` and
``polymer_conduit_params`` (see ``rmgpy/polymer.py`` M18.4 admit arm). The
SOLVER (``rmgpy/solver/polymer.pyx``) has NO dispatch arm for archetype 7, so
an admitted conduit row injects ZERO pool-moment flux today.

These tests encode the CORRECTED credit-only moment law (Codex round 80):

    Credit-only Dirac point bundle -- src_pool is NULL, NO source-pool debit,
    NO extra gas dn/dt write in the arm (the standard stoichiometric species
    writes already consume the reactant and produce the gas). For an event
    molar flux ``F`` [mol/s], destination monomer MW ``M``, defect ``d``, gas
    MW ``G``, and chain-units ``u = (sum(MW_reactants) - G + d) / M``:

        dmu0_dst += F
        dmu1_dst += F * u
        dmu2_dst += F * u * u

They are RED until the ``elif arch == FLUX_MOMENT_CREDIT_CONDUIT`` arm lands in
the residual (and the solver mirror enum gains ``FLUX_MOMENT_CREDIT_CONDUIT``).
Route (a) integration test: drive one admitted arch-7 row through
``initialize_model`` + ``residual`` and read the destination pool's mu slots.
"""

import numpy as np
import pytest

from rmgpy.kinetics import Arrhenius
from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.species import Species
from rmgpy.thermo import NASA, NASAPolynomial
from rmgpy.solver.polymer import HybridPolymerSystem, PolymerPoolConfig


def _spc(smiles: str, label: str) -> Species:
    s = Species(molecule=[Molecule().from_smiles(smiles)])
    s.label = label
    return s


def _trivial_nasa():
    """Minimal valid NASA thermo (constant 2.5R heat capacity) so a REVERSIBLE
    fixture survives generate_rate_coefficients' Keq path (identical thermo on
    both sides -> finite Keq). Only the reversible-refusal negative test needs
    it; the happy path is irreversible."""
    poly = NASAPolynomial(coeffs=[2.5, 0.0, 0.0, 0.0, 0.0, -745.375, 3.35],
                          Tmin=(200, "K"), Tmax=(6000, "K"))
    return NASA(polynomials=[poly], Tmin=(200, "K"), Tmax=(6000, "K"))


# Chosen, self-consistent numbers for the moment law + closure identity. The
# residual arm uses ONLY ``u`` (chain_units); M / d / G / sum(MW_reactants) are
# fixed here purely so the mass-closure identity is exact by construction:
#   u = (sum_mw_reactants - gas_mw + defect) / dst_monomer_mw
DST_MONOMER_MW = 28.0        # M  (destination pool monomer MW, g/mol)
GAS_MW = 16.0                # G  (single gas product MW, g/mol)
DEFECT = 2.0                 # d  (chain-scale discrete defect, g/mol)
SUM_MW_REACTANTS = 114.0     # sum(MW_reactants), g/mol
CHAIN_UNITS_U = (SUM_MW_REACTANTS - GAS_MW + DEFECT) / DST_MONOMER_MW  # u = 3.5714...


def _admitted_conduit_row():
    """One admitted arch-7 conduit row + its minimal solver.

    Shape (rate-carrying, NOT disqualified by the produce-then-transfer phase
    policy): a POLY-phase chain-scale discrete ``D`` is consumed; products are
    the gas product ``Gp`` plus the destination pool proxy ``A`` (so the solver
    resolves ``reaction_dst_pool`` to pool A while ``reaction_src_pool`` stays
    -1 -- D is a plain poly species, not a pool proxy). The dst credit lands on
    pool A's mu slots (core indices 1, 2, 3).
    """
    rs, _rxn, slots, i_D, i_Gp = _build_conduit_row()
    # Endpoint resolution invariant (already GREEN on the object side): the
    # admitted arch-7 row survives with src NULL and dst resolved to pool A.
    assert int(rs.reaction_flux_archetype[0]) == 7
    assert int(rs.reaction_src_pool[0]) == -1     # credit-only: NO source pool
    assert int(rs.reaction_dst_pool[0]) == 0      # destination = pool A
    assert int(rs.reaction_refused[0]) == 0       # admitted, live
    return rs, slots, i_D, i_Gp


def _default_conduit_params():
    """The exact ``polymer_conduit_params`` dict the M18.4 admit arm stamps
    (rmgpy/polymer.py). ``chain_units`` u is pinned to the mass-closure
    identity u = (sum_MW_reactants - gas_MW + defect) / dst_monomer_MW."""
    return {
        "admission_direction": "chain_to_pool",
        "chain_units": float(CHAIN_UNITS_U),
        "gas_products": [{"species": "Gp", "stoich": 1,
                          "mw_g_mol": float(GAS_MW)}],
        "gas_units": float(GAS_MW / DST_MONOMER_MW),
        "candidate_key": "tdd-arch7-key",
    }


def _build_conduit_row(reversible=False, params="__default__", dst_pool="A",
                       attach_thermo=False, rs_kwargs=None):
    """Build + initialize a one-row arch-7 solver, WITHOUT any admitted/refused
    assertions (negative tests inspect the outcome themselves).

    Shape (rate-carrying, NOT disqualified by the produce-then-transfer phase
    policy): a POLY-phase chain-scale discrete ``D`` is consumed; products are
    the gas product ``Gp`` plus the destination pool proxy ``A`` (so the solver
    resolves ``reaction_dst_pool`` to pool A while ``reaction_src_pool`` stays
    -1 -- D is a plain poly species, not a pool proxy). The dst credit lands on
    pool A's mu slots (core indices 1, 2, 3). The destination pool carries the
    REAL ``chain_mass_defect_g_mol = DEFECT`` so the solver's condensed-mass
    law (PolymerPoolConfig.condensed_mass_g) uses the SAME defect the closure
    assertion pins.
    """
    A = _spc("CCCC", "A")            # dst pool proxy      (core 0, poly)
    A_mu0 = _spc("CO", "A_mu0")      # dst mu0             (core 1)
    A_mu1 = _spc("C=O", "A_mu1")     # dst mu1             (core 2)
    A_mu2 = _spc("C#N", "A_mu2")     # dst mu2             (core 3)
    D = _spc("CCCCCCCC", "D")        # chain-scale discrete reactant (core 4, poly)
    Gp = _spc("[CH3]", "Gp")         # single gas product  (core 5)
    core = [A, A_mu0, A_mu1, A_mu2, D, Gp]
    mask = np.array([False, False, False, False, False, True], dtype=bool)
    if attach_thermo:
        for s in core:
            s.thermo = _trivial_nasa()

    rxn = Reaction(
        reactants=[D],
        products=[Gp, A],
        kinetics=Arrhenius(A=(2.0, "1/s"), n=0.0, Ea=(0.0, "kcal/mol"),
                           T0=(298.15, "K")),
        reversible=reversible,
    )
    # Stamp exactly what the M18.4 admit arm stamps (rmgpy/polymer.py).
    rxn.polymer_flux_archetype = 7   # PolymerFluxArchetype.MOMENT_CREDIT_CONDUIT
    rxn.polymer_conduit_dst_pool = dst_pool
    rxn.polymer_conduit_params = (_default_conduit_params()
                                  if params == "__default__" else params)

    pool_a = PolymerPoolConfig(
        label="A", xs=2, explicit_dp_to_species_index={},
        mu_indices=(1, 2, 3), monomer_poly_index=None,
        k_scission=0.0, k_unzip=0.0, tail_kinetics=None,
        monomer_mw_g_mol=DST_MONOMER_MW,
        chain_mass_defect_g_mol=DEFECT,
    )
    kwargs = dict(
        T=800.0, P=1.0e5, initial_mole_fractions={D: 1.0, Gp: 0.0}, V_poly=1.0,
        polymer_pools=[pool_a], mass_transfer=[],
        gas_species_mask=mask.copy(), constant_gas_volume=False,
        initial_polymer_moments={"A": (1.0, 5.0, 30.0)}, termination=[],
        allow_unstamped_proxy_rows=True,
    )
    kwargs.update(rs_kwargs or {})
    rs = HybridPolymerSystem(**kwargs)
    rs.initialize_model(core, [rxn], [], [])
    # (mu slots, reactant D index, gas product Gp index)
    return rs, rxn, (1, 2, 3), 4, 5


class TestConduitMomentCreditArchetype:

    def test_arch7_solver_enum_constant_matches_object_enum(self):
        """The solver mirror must gain FLUX_MOMENT_CREDIT_CONDUIT == 7,
        matching PolymerFluxArchetype.MOMENT_CREDIT_CONDUIT (pinned like the
        other archetypes by test_flux_archetype_constants_match_enum). RED
        now: the solver enum stops at FLUX_VOLATILE_EJECTION = 6."""
        from rmgpy.polymer import PolymerFluxArchetype
        import rmgpy.solver.polymer as solver_mod
        assert hasattr(solver_mod, "FLUX_MOMENT_CREDIT_CONDUIT"), (
            "solver mirror missing FLUX_MOMENT_CREDIT_CONDUIT (arch 7)")
        assert (solver_mod.FLUX_MOMENT_CREDIT_CONDUIT
                == int(PolymerFluxArchetype.MOMENT_CREDIT_CONDUIT) == 7)

    def test_arch7_moment_credit_law_and_closure(self):
        """Credit-only Dirac point-bundle law for one admitted arch-7 row.

        RED now: the solver has no archetype-7 dispatch arm, so the dst pool's
        mu slots receive ZERO credit (dmu0_dst == 0 != F). GREEN once the
        `elif arch == FLUX_MOMENT_CREDIT_CONDUIT` residual arm injects
        (F, F*u, F*u*u) onto the destination pool moments with NO source debit
        and NO extra gas write.
        """
        rs, (i_mu0, i_mu1, i_mu2), i_D, i_Gp = _admitted_conduit_row()

        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]

        # Event molar flux F [mol/s]: the reactant chain-scale discrete is
        # consumed ONLY by this row (stoich 1), and the arm adds NO source
        # debit, so -dn_dt[D] == F both before and after the arm lands.
        F = -dn_dt[i_D]
        assert F > 0.0, (
            "conduit row carries zero rate; fixture is not exercising the "
            "arm (should be F = k*[D]*V_poly = 2.0 mol/s)")
        # Standard stoichiometric gas write (NOT the arm): gas gains +F.
        assert np.isclose(dn_dt[i_Gp], F, rtol=0, atol=1e-12)

        u = CHAIN_UNITS_U
        dmu0 = dn_dt[i_mu0]
        dmu1 = dn_dt[i_mu1]
        dmu2 = dn_dt[i_mu2]

        # --- Core credit-only moment law (the tripwire) --------------------
        assert np.isclose(dmu0, F, rtol=0, atol=1e-12), (
            f"dst mu0 credit {dmu0!r} != F {F!r} -- arch-7 arm absent "
            f"(zero conduit flux injected)")
        assert np.isclose(dmu1, F * u, rtol=0, atol=1e-12), (
            f"dst mu1 credit {dmu1!r} != F*u {F * u!r}")
        assert np.isclose(dmu2, F * u * u, rtol=0, atol=1e-12), (
            f"dst mu2 credit {dmu2!r} != F*u^2 {F * u * u!r}")

        # --- Dirac point-bundle consistency: mu2*mu0 == mu1^2 --------------
        assert np.isclose(dmu2 * dmu0, dmu1 * dmu1, rtol=1e-12, atol=1e-18), (
            "credited bundle is not a Dirac point mass (mu2*mu0 != mu1^2)")

        # --- Mass closure (algebraic identity): M*dmu1 - d*dmu0 ============
        lhs = DST_MONOMER_MW * dmu1 - DEFECT * dmu0
        rhs = F * (SUM_MW_REACTANTS - GAS_MW)
        assert np.isclose(lhs, rhs, rtol=1e-12, atol=1e-9), (
            f"mass not conserved by the credited bundle: "
            f"M*dmu1 - d*dmu0 = {lhs!r} != F*(sumMW - G) = {rhs!r}")

        # --- Mass closure against the SOLVER's ACTUAL mass law (FIX 1) ------
        # The dst pool carries the REAL chain_mass_defect_g_mol = DEFECT, so
        # the solver's own condensed-mass accessor uses the SAME defect the
        # closure pins. Two independent solver-side reads:
        #  (1) per-pool law applied to the moment RATES (exact, linear); and
        #  (2) a finite-difference of the total condensed-mass accessor along
        #      dn_dt -- the end-to-end quantity r86 terminationPolymerConversion
        #      actually integrates. Both must equal F*(sumMW - gas_MW): the
        #      discrete carbon leaving == pool condensed mass gained + gas out.
        pool = rs.polymer_pools[0]
        dmass_pool = pool.condensed_mass_g(dmu0, dmu1)  # = M*dmu1 - d*dmu0
        assert np.isclose(dmass_pool, rhs, rtol=1e-12, atol=1e-9), (
            f"SOLVER pool mass law {dmass_pool!r} != F*(sumMW - G) {rhs!r}")

        dt = 1.0e-4
        mass0 = rs.get_total_polymer_condensed_mass_g(rs.y)
        mass1 = rs.get_total_polymer_condensed_mass_g(rs.y + dt * dn_dt)
        dmass_dt = (mass1 - mass0) / dt
        assert np.isclose(dmass_dt, rhs, rtol=1e-7, atol=1e-6), (
            f"SOLVER d(condensed_mass)/dt {dmass_dt!r} != F*(sumMW - G) "
            f"{rhs!r}; the pool mass law (defect={DEFECT}) disagrees with the "
            f"credited bundle -- phantom condensed mass")

        # The DEFECT must be load-bearing: with d=0 the solver mass law would
        # read M*dmu1 (!= rhs), so a zero-defect pool config would FAIL the
        # checks above. This guards the FIX-1 fixture from silently regressing
        # to a defect-free false-green.
        assert not np.isclose(DST_MONOMER_MW * dmu1, rhs), (
            "DEFECT is not load-bearing (M*dmu1 already equals rhs) -- the "
            "solver-side mass check would pass even with a zero defect")

        # --- No source-pool debit: the ONLY pool moments touched are dst.
        # (Single-pool fixture, so this also pins 'credit-only'.)
        assert dmu0 >= 0.0 and dmu1 >= 0.0 and dmu2 >= 0.0

    # ---------------------------------------------------------------------
    # Negative tests (FIX 6): a mis-stamped arch-7 row must be REFUSED
    # (flux-dead) at init, never credited. Each asserts the row is marked
    # refused + demoted to FLUX_NONE and injects ZERO dst-pool moment credit.
    # ---------------------------------------------------------------------

    @staticmethod
    def _assert_refused_and_zero_credit(rs):
        assert int(rs.reaction_refused[0]) == 1, "row was not refused"
        import rmgpy.solver.polymer as solver_mod
        assert int(rs.reaction_flux_archetype[0]) == solver_mod.FLUX_NONE, (
            "refused row was not demoted to FLUX_NONE")
        dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
        # Zero credit on dst pool mu slots (1, 2, 3) -- and, since a refused
        # row hits the RHS continue guard, zero flux everywhere.
        assert dn_dt[1] == 0.0 and dn_dt[2] == 0.0 and dn_dt[3] == 0.0, (
            f"refused arch-7 row still credited dst moments: "
            f"{dn_dt[1]!r}, {dn_dt[2]!r}, {dn_dt[3]!r}")

    def test_arch7_reversible_row_is_refused(self):
        """A REVERSIBLE arch-7 stamp is refused at init: r_mol_s = rf - rr is
        signed and has no reverse moment law (admitted rows must be
        irreversible-export rewritten). Needs thermo so the Keq path in
        generate_rate_coefficients does not choke before the refusal is read.
        """
        rs, _rxn, _slots, _iD, _iGp = _build_conduit_row(
            reversible=True, attach_thermo=True,
            rs_kwargs={"allow_unpaired_reference_state": True})
        self._assert_refused_and_zero_credit(rs)

    def test_arch7_missing_chain_units_is_refused(self):
        params = _default_conduit_params()
        del params["chain_units"]
        rs, _rxn, _slots, _iD, _iGp = _build_conduit_row(params=params)
        self._assert_refused_and_zero_credit(rs)

    def test_arch7_chain_units_below_one_is_refused(self):
        params = _default_conduit_params()
        params["chain_units"] = 0.5           # u < 1.0
        rs, _rxn, _slots, _iD, _iGp = _build_conduit_row(params=params)
        self._assert_refused_and_zero_credit(rs)

    def test_arch7_nonpositive_gas_mw_is_refused(self):
        params = _default_conduit_params()
        params["gas_products"][0]["mw_g_mol"] = 0.0   # gas MW <= 0
        rs, _rxn, _slots, _iD, _iGp = _build_conduit_row(params=params)
        self._assert_refused_and_zero_credit(rs)

    def test_arch7_nonfinite_gas_mw_is_refused(self):
        """A NaN/inf gas MW must be refused too: ``mw <= 0`` alone lets
        non-finite MW slip through (NaN/inf both compare False to <= 0),
        which would store a poisoned mw for the accepted-step grams
        booking. Validation must gate on math.isfinite."""
        for bad_mw in (float("nan"), float("inf")):
            params = _default_conduit_params()
            params["gas_products"][0]["mw_g_mol"] = bad_mw
            rs, _rxn, _slots, _iD, _iGp = _build_conduit_row(params=params)
            self._assert_refused_and_zero_credit(rs)

    def test_arch7_missing_candidate_key_is_refused(self):
        params = _default_conduit_params()
        params["candidate_key"] = ""          # missing/empty key
        rs, _rxn, _slots, _iD, _iGp = _build_conduit_row(params=params)
        self._assert_refused_and_zero_credit(rs)

    def test_arch7_unresolvable_dst_pool_is_refused(self):
        """dst-less path: the stamped destination label resolves to no
        configured pool -> refuse (cannot credit anything)."""
        rs, _rxn, _slots, _iD, _iGp = _build_conduit_row(
            dst_pool="NO_SUCH_POOL")
        self._assert_refused_and_zero_credit(rs)
