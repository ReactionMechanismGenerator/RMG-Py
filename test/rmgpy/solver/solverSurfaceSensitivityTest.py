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


"""
End-to-end regression test for sensitivity analysis in a SurfaceReactor with
coverage-dependent kinetics and thermochemistry.

The sensitivity coefficients come from DASPK integrating the sensitivity
equations d(dn/dp)/dt = J dn/dp + df/dp, so they depend on both the state
Jacobian J and the parameter derivatives df/dp being consistent with the
residual. This test checks the result against brute force: re-simulating with
each parameter perturbed and taking central differences.
"""

import logging

import numpy as np
import pytest

from rmgpy.kinetics import SurfaceArrhenius
from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.solver.surface import SurfaceReactor
from rmgpy.species import Species
from rmgpy.thermo import ThermoData

T_END = 1e-5  # s; mid-transient, with theta_HX rising from 0.1 to about 0.3
HX_ADJ = "1 H u0 p0 {2,S} \n 2 X u0 p0 {1,S}"
OCX_ADJ = "1 O u0 p2 c0 {2,D} \n 2 C u0 p0 c0 {1,D} {3,D} \n 3 X u0 p0 c0 {2,D}"


def _polynomial(h1):
    return {"model": "polynomial",
            "enthalpy-coefficients": [(h1, "J/mol"), (0, "J/mol"), (0, "J/mol")],
            "entropy-coefficients": [(0, "J/(mol*K)"), (0, "J/(mol*K)"), (0, "J/(mol*K)")]}


def _simulate(thermo_coverage, kinetic_coverage, sensitivity=False, k_scale=(1.0, 1.0), dH=(0.0, 0.0)):
    """
    Integrate H2 + 2X <=> 2HX and CO + X <=> OCX to T_END and return the
    reactor. HX and OCX interact through their thermo coverage dependence, and
    each reaction's kinetics depend on the coverage of species that are not
    its reactants. ``k_scale`` scales each reaction's A factor (so kf and kb,
    with Keq fixed) and ``dH`` shifts the enthalpy, so also the free energy,
    of HX and OCX, in J/mol.
    """
    def thermo(H298, S298, coverage_dependence=None):
        return ThermoData(Tdata=([300, 400, 500, 600, 800, 1000, 1500], "K"),
                          Cpdata=([5.0] * 7, "cal/(mol*K)"),
                          H298=(H298, "kcal/mol"), S298=(S298, "cal/(mol*K)"),
                          thermo_coverage_dependence=coverage_dependence if thermo_coverage else None)

    h2 = Species(molecule=[Molecule().from_smiles("[H][H]")], thermo=thermo(0.0, 31.1))
    co = Species(molecule=[Molecule().from_smiles("[C-]#[O+]")], thermo=thermo(-26.4, 47.2))
    x = Species(molecule=[Molecule().from_adjacency_list("1 X u0 p0")], thermo=thermo(0.0, 0.0))
    hx = Species(molecule=[Molecule().from_adjacency_list(HX_ADJ)],
                 thermo=thermo(-11.3 + dH[0] / 4184, 0.4, {HX_ADJ: _polynomial(20000), OCX_ADJ: _polynomial(5000)}))
    ocx = Species(molecule=[Molecule().from_adjacency_list(OCX_ADJ)],
                  thermo=thermo(-45.0 + dH[1] / 4184, 10.0, {OCX_ADJ: _polynomial(15000), HX_ADJ: _polynomial(5000)}))
    rxn1 = Reaction(
        reactants=[h2, x, x],
        products=[hx, hx],
        kinetics=SurfaceArrhenius(
            A=(9.05e15 * k_scale[0], "cm^5/(mol^2*s)"), n=0.5, Ea=(5.0, "kJ/mol"), T0=(1.0, "K"),
            coverage_dependence={x: {"a": 1.0, "m": -0.5, "E": (0.0, "J/mol")},
                                 ocx: {"a": 0.0, "m": 0.0, "E": (10.0, "kJ/mol")}} if kinetic_coverage else None,
        ),
    )
    rxn2 = Reaction(
        reactants=[co, x],
        products=[ocx],
        kinetics=SurfaceArrhenius(
            A=(1.0e5 * k_scale[1], "m^3/(mol*s)"), n=0.0, Ea=(0.0, "kJ/mol"), T0=(1.0, "K"),
            coverage_dependence={hx: {"a": -1.0, "m": 0.0, "E": (5.0, "kJ/mol")}} if kinetic_coverage else None,
        ),
    )
    reactor = SurfaceReactor(
        600,
        1.0e5,
        n_sims=1,
        initial_gas_mole_fractions={h2: 0.5, co: 0.5},
        # Start every adsorbate above the 1e-6 coverage cutoff, where the kinetic
        # coverage factor switches on abruptly, so the integration never crosses it.
        initial_surface_coverages={x: 0.8, hx: 0.1, ocx: 0.1},
        surface_volume_ratio=(1e1, "m^-1"),
        surface_site_density=(2.72e-9, "mol/cm^2"),
        coverage_dependence=kinetic_coverage,
        thermo_coverage_dependence=thermo_coverage,
        termination=[],
    )
    reactor.initialize_model([h2, co, x, hx, ocx], [rxn1, rxn2], [], [], sensitivity=sensitivity,
                             atol=1e-20, rtol=1e-10, sens_atol=1e-12, sens_rtol=1e-8)
    # Integrate to exactly T_END: the termination-based simulate() stops on the first
    # step past the end time, which differs between runs and spoils the differences.
    reactor.advance(T_END)
    return reactor


@pytest.mark.parametrize("thermo_coverage, kinetic_coverage",
                         [(False, False), (True, False), (False, True), (True, True)])
def test_sensitivities_match_brute_force(thermo_coverage, kinetic_coverage):
    """
    d n_i / d ln(k_j) and d n_i / d G_s from the sensitivity equations must
    match central differences of re-simulations, for every species i (gas and
    surface), both reactions j and the free energies of the two adsorbates s.

    Both are compared in moles, not as logarithmic sensitivities. There is far
    more gas (0.5 mol of each) than surface sites (about 1e-5 mol), so
    d ln n / d ln k of a gas species is about 1e-5 times that of an adsorbate,
    and a tolerance relative to the largest entry would not test it at all.
    In moles, the gas and surface entries of a column are comparable.
    """
    logging.disable(logging.WARNING)  # DASPK step-count notices
    try:
        reactor = _simulate(thermo_coverage, kinetic_coverage, sensitivity=True)
        n = reactor.num_core_species
        n_rxn = reactor.num_core_reactions
        # Sensitivity block p holds d n / d p_p, for p = kf of each reaction, then G of each species
        dn_dp = reactor.y[n:].reshape(n_rxn + n, n)

        # A gas amount changes only by about 1e-9 mol per 1e-3 change in ln k, so a smaller
        # step lets the solver's error in the 0.5 mol of gas dominate the difference.
        rel = 1e-2
        for j in range(n_rxn):
            up, down = [1.0, 1.0], [1.0, 1.0]
            up[j], down[j] = 1 + rel, 1 - rel
            n_plus = _simulate(thermo_coverage, kinetic_coverage, k_scale=tuple(up)).y[:n]
            n_minus = _simulate(thermo_coverage, kinetic_coverage, k_scale=tuple(down)).y[:n]
            brute_force = (n_plus - n_minus) / (np.log(1 + rel) - np.log(1 - rel))
            from_sensitivity = dn_dp[j] * reactor.kf[j]
            assert np.allclose(from_sensitivity, brute_force, rtol=0, atol=2e-3 * np.abs(brute_force).max()), \
                f"d n / d ln k for reaction {j}"

        dG = 10.0  # J/mol
        for k, species_index in enumerate((3, 4)):  # HX, OCX
            up, down = [0.0, 0.0], [0.0, 0.0]
            up[k], down[k] = dG, -dG
            n_plus = _simulate(thermo_coverage, kinetic_coverage, dH=tuple(up)).y[:n]
            n_minus = _simulate(thermo_coverage, kinetic_coverage, dH=tuple(down)).y[:n]
            brute_force = (n_plus - n_minus) / (2 * dG)
            from_sensitivity = dn_dp[n_rxn + species_index]
            assert np.allclose(from_sensitivity, brute_force, rtol=0, atol=2e-3 * np.abs(brute_force).max()), \
                f"d n / d G for species {species_index}"
    finally:
        logging.disable(logging.NOTSET)
