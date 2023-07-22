#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2021 Prof. William H. Green (whgreen@mit.edu),           #
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


import time

import numpy as np

import rmgpy.constants as constants
from rmgpy.kinetics import SurfaceArrhenius, StickingCoefficient
from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.solver.surface import SurfaceReactor
from rmgpy.species import Species
from rmgpy.thermo import ThermoData, NASA, NASAPolynomial


class SurfaceReactorCheck:
    def test_solve_h2(self):
        """
        Test the surface batch reactor with a dissociative adsorption of H2

        Here we choose a kinetic model consisting of the dissociative adsorption reaction
        H2 + 2X <=> 2 HX
        We use a SurfaceArrhenius for the rate expression.
        """
        h2 = Species(
            molecule=[Molecule().from_smiles("[H][H]")],
            thermo=ThermoData(
                Tdata=([300, 400, 500, 600, 800, 1000, 1500], "K"),
                Cpdata=(
                    [6.955, 6.955, 6.956, 6.961, 7.003, 7.103, 7.502],
                    "cal/(mol*K)",
                ),
                H298=(0, "kcal/mol"),
                S298=(31.129, "cal/(mol*K)"),
            ),
        )
        x = Species(
            molecule=[Molecule().from_adjacency_list("1 X u0 p0")],
            thermo=ThermoData(
                Tdata=([300, 400, 500, 600, 800, 1000, 1500], "K"),
                Cpdata=([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], "cal/(mol*K)"),
                H298=(0.0, "kcal/mol"),
                S298=(0.0, "cal/(mol*K)"),
            ),
        )
        hx = Species(
            molecule=[Molecule().from_adjacency_list("1 H u0 p0 {2,S} \n 2 X u0 p0 {1,S}")],
            thermo=ThermoData(
                Tdata=([300, 400, 500, 600, 800, 1000, 1500], "K"),
                Cpdata=([1.50, 2.58, 3.40, 4.00, 4.73, 5.13, 5.57], "cal/(mol*K)"),
                H298=(-11.26, "kcal/mol"),
                S298=(0.44, "cal/(mol*K)"),
            ),
        )

        rxn1 = Reaction(
            reactants=[h2, x, x],
            products=[hx, hx],
            kinetics=SurfaceArrhenius(A=(9.05e18, "cm^5/(mol^2*s)"), n=0.5, Ea=(5.0, "kJ/mol"), T0=(1.0, "K")),
        )

        core_species = [h2, x, hx]
        edge_species = []
        core_reactions = [rxn1]
        edge_reactions = []

        T = 600
        P_initial = 1.0e5
        rxn_system = SurfaceReactor(
            T,
            P_initial,
            n_sims=1,
            initial_gas_mole_fractions={h2: 1.0},
            initial_surface_coverages={x: 1.0},
            surface_volume_ratio=(1e1, "m^-1"),
            surface_site_density=(2.72e-9, "mol/cm^2"),
            termination=[],
        )

        rxn_system.initialize_model(core_species, core_reactions, edge_species, edge_reactions)

        tlist = np.logspace(-13, -5, 81, dtype=float)

        # Integrate to get the solution at each time point
        t = []
        y = []
        reaction_rates = []
        species_rates = []
        for t1 in tlist:
            rxn_system.advance(t1)
            t.append(rxn_system.t)
            # You must make a copy of y because it is overwritten by DASSL at
            # each call to advance()
            y.append(rxn_system.y.copy())
            reaction_rates.append(rxn_system.core_reaction_rates.copy())
            species_rates.append(rxn_system.core_species_rates.copy())

        # Convert the solution vectors to np arrays
        t = np.array(t, float)
        y = np.array(y, float)
        reaction_rates = np.array(reaction_rates, float)
        species_rates = np.array(species_rates, float)
        total_sites = y[0, 1]

        # Check that we're computing the species fluxes correctly
        for i in range(t.shape[0]):
            assert abs(reaction_rates[i, 0] - -1.0 * species_rates[i, 0]) < 1e-6 * reaction_rates[i, 0]
            assert abs(reaction_rates[i, 0] - -0.5 * species_rates[i, 1]) < 1e-6 * reaction_rates[i, 0]
            assert abs(reaction_rates[i, 0] - 0.5 * species_rates[i, 2]) < 1e-6 * reaction_rates[i, 0]

        # Check that we've reached equilibrium
        assert abs(reaction_rates[-1, 0] - 0.0) < 1e-2

        # # Visualize the simulation results
        # import pylab
        # fig = pylab.figure(figsize=(6, 6))
        # pylab.subplot(2, 1, 1)
        # pylab.semilogx(t, y[:, 2] / total_sites)
        # pylab.ylabel('Surface coverage')
        # pylab.legend(['HX'], loc=4)
        # pylab.subplot(2, 1, 2)
        # pylab.semilogx(t, species_rates)
        # pylab.legend(['H2', 'X', 'HX'], loc=4)
        # pylab.xlabel('Time (s)')
        # pylab.ylabel('Rate (mol/m$^\\mathdefault{3 or 2}$*s)')
        # # fig.subplots_adjust(left=0.21, bottom=0.10, right=0.95, top=0.95, wspace=0.20, hspace=0.35)
        # pylab.tight_layout()
        # # pylab.show()
        # pylab.savefig('surfaceTestH2.pdf')
        # return

    def test_solve_ch3(self):
        """
        Test the surface batch reactor with a nondissociative adsorption of CH3

        Here we choose a kinetic model consisting of the  adsorption reaction
        CH3 + X <=>  CH3X
        We use a sticking coefficient for the rate expression.
        """

        ch3 = Species(
            molecule=[Molecule().from_smiles("[CH3]")],
            thermo=NASA(
                polynomials=[
                    NASAPolynomial(
                        coeffs=[
                            3.91547,
                            0.00184155,
                            3.48741e-06,
                            -3.32746e-09,
                            8.49953e-13,
                            16285.6,
                            0.351743,
                        ],
                        Tmin=(100, "K"),
                        Tmax=(1337.63, "K"),
                    ),
                    NASAPolynomial(
                        coeffs=[
                            3.54146,
                            0.00476786,
                            -1.82148e-06,
                            3.28876e-10,
                            -2.22545e-14,
                            16224,
                            1.66032,
                        ],
                        Tmin=(1337.63, "K"),
                        Tmax=(5000, "K"),
                    ),
                ],
                Tmin=(100, "K"),
                Tmax=(5000, "K"),
                E0=(135.382, "kJ/mol"),
                comment="""Thermo library: primaryThermoLibrary + radical(CH3)""",
            ),
            molecular_weight=(15.0345, "amu"),
        )

        x = Species(
            molecule=[Molecule().from_adjacency_list("1 X u0 p0")],
            thermo=NASA(
                polynomials=[
                    NASAPolynomial(coeffs=[0, 0, 0, 0, 0, 0, 0], Tmin=(298, "K"), Tmax=(1000, "K")),
                    NASAPolynomial(coeffs=[0, 0, 0, 0, 0, 0, 0], Tmin=(1000, "K"), Tmax=(2000, "K")),
                ],
                Tmin=(298, "K"),
                Tmax=(2000, "K"),
                E0=(-6.19426, "kJ/mol"),
                comment="""Thermo library: surfaceThermo""",
            ),
        )

        ch3x = Species(
            molecule=[
                Molecule().from_adjacency_list(
                    """1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 X u0 p0 c0 {1,S}"""
                )
            ],
            thermo=NASA(
                polynomials=[
                    NASAPolynomial(
                        coeffs=[
                            -0.552219,
                            0.026442,
                            -3.55617e-05,
                            2.60044e-08,
                            -7.52707e-12,
                            -4433.47,
                            0.692144,
                        ],
                        Tmin=(298, "K"),
                        Tmax=(1000, "K"),
                    ),
                    NASAPolynomial(
                        coeffs=[
                            3.62557,
                            0.00739512,
                            -2.43797e-06,
                            1.86159e-10,
                            3.6485e-14,
                            -5187.22,
                            -18.9668,
                        ],
                        Tmin=(1000, "K"),
                        Tmax=(2000, "K"),
                    ),
                ],
                Tmin=(298, "K"),
                Tmax=(2000, "K"),
                E0=(-39.1285, "kJ/mol"),
                comment="""Thermo library: surfaceThermoNi111""",
            ),
        )

        rxn1 = Reaction(
            reactants=[ch3, x],
            products=[ch3x],
            kinetics=StickingCoefficient(
                A=0.1,
                n=0,
                Ea=(0, "kcal/mol"),
                T0=(1, "K"),
                Tmin=(200, "K"),
                Tmax=(3000, "K"),
                comment="""Exact match found for rate rule (Adsorbate;VacantSite)""",
            )
            # kinetics=SurfaceArrhenius(A=(2.7e10, 'cm^3/(mol*s)'),
            #                           n=0.5,
            #                           Ea=(5.0, 'kJ/mol'),
            #                           T0=(1.0, 'K'))
        )
        core_species = [ch3, x, ch3x]
        edge_species = []
        core_reactions = [rxn1]
        edge_reactions = []

        T = 800.0
        P_initial = 1.0e5
        rxn_system = SurfaceReactor(
            T,
            P_initial,
            n_sims=1,
            initial_gas_mole_fractions={ch3: 1.0},
            initial_surface_coverages={x: 1.0},
            surface_volume_ratio=(1.0, "m^-1"),
            surface_site_density=(2.72e-9, "mol/cm^2"),
            termination=[],
        )
        # in chemkin, the sites are mostly occupied in about 1e-8 seconds.

        rxn_system.initialize_model(core_species, core_reactions, edge_species, edge_reactions)

        tlist = np.logspace(-13, -5, 81, dtype=float)

        print("Surface site density:", rxn_system.surface_site_density.value_si)

        print(
            "rxn1 rate coefficient",
            rxn1.get_surface_rate_coefficient(rxn_system.T.value_si, rxn_system.surface_site_density.value_si),
        )

        # Integrate to get the solution at each time point
        t = []
        y = []
        reaction_rates = []
        species_rates = []
        t.append(rxn_system.t)
        # You must make a copy of y because it is overwritten by DASSL at
        # each call to advance()
        y.append(rxn_system.y.copy())
        reaction_rates.append(rxn_system.core_reaction_rates.copy())
        species_rates.append(rxn_system.core_species_rates.copy())
        print("time: ", t)
        print("moles:", y)
        print("reaction rates:", reaction_rates)
        print("species rates:", species_rates)
        for t1 in tlist:
            rxn_system.advance(t1)
            t.append(rxn_system.t)
            # You must make a copy of y because it is overwritten by DASSL at
            # each call to advance()
            y.append(rxn_system.y.copy())
            reaction_rates.append(rxn_system.core_reaction_rates.copy())
            species_rates.append(rxn_system.core_species_rates.copy())

        # Convert the solution vectors to np arrays
        t = np.array(t, float)
        y = np.array(y, float)
        reaction_rates = np.array(reaction_rates, float)
        species_rates = np.array(species_rates, float)
        V = constants.R * rxn_system.T.value_si * np.sum(y) / rxn_system.P_initial.value_si

        # Check that we're computing the species fluxes correctly
        for i in range(t.shape[0]):
            assert abs(reaction_rates[i, 0] - -species_rates[i, 0]) < 1e-6 * reaction_rates[i, 0]
            assert abs(reaction_rates[i, 0] - -species_rates[i, 1]) < 1e-6 * reaction_rates[i, 0]
            assert abs(reaction_rates[i, 0] - species_rates[i, 2]) < 1e-6 * reaction_rates[i, 0]

        # Check that we've reached equilibrium by the end
        assert abs(reaction_rates[-1, 0] - 0.0) < 1e-2

    def test_solve_h2_coverage_dependence(self):
        """
        Test the surface batch reactor can properly apply coverage dependent parameters
        with the dissociative adsorption of H2.

        Here we choose a kinetic model consisting of the dissociative adsorption reaction
        H2 + 2X <=> 2 HX
        We use a SurfaceArrhenius for the rate expression.
        """
        h2 = Species(
            molecule=[Molecule().from_smiles("[H][H]")],
            thermo=ThermoData(
                Tdata=([300, 400, 500, 600, 800, 1000, 1500], "K"),
                Cpdata=(
                    [6.955, 6.955, 6.956, 6.961, 7.003, 7.103, 7.502],
                    "cal/(mol*K)",
                ),
                H298=(0, "kcal/mol"),
                S298=(31.129, "cal/(mol*K)"),
            ),
        )
        x = Species(
            molecule=[Molecule().from_adjacency_list("1 X u0 p0")],
            thermo=ThermoData(
                Tdata=([300, 400, 500, 600, 800, 1000, 1500], "K"),
                Cpdata=([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], "cal/(mol*K)"),
                H298=(0.0, "kcal/mol"),
                S298=(0.0, "cal/(mol*K)"),
            ),
        )
        hx = Species(
            molecule=[Molecule().from_adjacency_list("1 H u0 p0 {2,S} \n 2 X u0 p0 {1,S}")],
            thermo=ThermoData(
                Tdata=([300, 400, 500, 600, 800, 1000, 1500], "K"),
                Cpdata=([1.50, 2.58, 3.40, 4.00, 4.73, 5.13, 5.57], "cal/(mol*K)"),
                H298=(-11.26, "kcal/mol"),
                S298=(0.44, "cal/(mol*K)"),
            ),
        )

        rxn1 = Reaction(
            reactants=[h2, x, x],
            products=[hx, hx],
            kinetics=SurfaceArrhenius(
                A=(9.05e18, "cm^5/(mol^2*s)"),
                n=0.5,
                Ea=(5.0, "kJ/mol"),
                T0=(1.0, "K"),
                coverage_dependence={x: {"a": 0.0, "m": -1.0, "E": (0.0, "J/mol")}},
            ),
        )

        rxn2 = Reaction(
            reactants=[h2, x, x],
            products=[hx, hx],
            kinetics=SurfaceArrhenius(
                A=(9.05e-18, "cm^5/(mol^2*s)"),  # 1e36 times slower
                n=0.5,
                Ea=(5.0, "kJ/mol"),
                T0=(1.0, "K"),
                coverage_dependence={x: {"a": 0.0, "m": -1.0, "E": (10.0, "J/mol")}},
            ),
        )

        core_species = [h2, x, hx]
        edge_species = []
        core_reactions = [rxn1]
        edge_reactions = []

        # make it slower, for benchmarking
        for j in range(200):
            core_species.append(hx.copy())
        for j in range(1000):
            core_reactions.append(rxn2)

        T = 600
        P_initial = 1.0e5
        rxn_system = SurfaceReactor(
            T,
            P_initial,
            n_sims=1,
            initial_gas_mole_fractions={h2: 1.0},
            initial_surface_coverages={x: 1.0},
            surface_volume_ratio=(1e1, "m^-1"),
            surface_site_density=(2.72e-9, "mol/cm^2"),
            coverage_dependence=True,
            termination=[],
        )

        rxn_system.initialize_model(core_species, core_reactions, edge_species, edge_reactions)

        tlist = np.logspace(-13, -5, 81, dtype=float)

        assert isinstance(rxn1.kinetics.coverage_dependence, dict)  # check to make sure coverage_dependence is still the correct type
        for species, parameters in rxn1.kinetics.coverage_dependence.items():
            assert isinstance(species, Species)  # species should be a Species
            assert isinstance(parameters, dict)
            assert parameters["a"] is not None
            assert parameters["m"] is not None
            assert parameters["E"] is not None

        # Integrate to get the solution at each time point
        t = []
        y = []
        reaction_rates = []
        species_rates = []
        start_time = time.time()
        for t1 in tlist:
            rxn_system.advance(t1)
            t.append(rxn_system.t)
            # You must make a copy of y because it is overwritten by DASSL at
            # each call to advance()
            y.append(rxn_system.y.copy())
            reaction_rates.append(rxn_system.core_reaction_rates.copy())
            species_rates.append(rxn_system.core_species_rates.copy())
        run_time = time.time() - start_time
        print(f"Simulation took {run_time:.3e} seconds in {self.id()}")

        # Convert the solution vectors to np arrays
        t = np.array(t, float)
        y = np.array(y, float)
        reaction_rates = np.array(reaction_rates, float)
        species_rates = np.array(species_rates, float)
        total_sites = y[0, 1]

        # Check that we're computing the species fluxes correctly
        for i in range(t.shape[0]):
            assert abs(reaction_rates[i, 0] - -1.0 * species_rates[i, 0]) < 1e-6 * reaction_rates[i, 0]
            assert abs(reaction_rates[i, 0] - -0.5 * species_rates[i, 1]) < 1e-6 * reaction_rates[i, 0]
            assert abs(reaction_rates[i, 0] - 0.5 * species_rates[i, 2]) < 1e-6 * reaction_rates[i, 0]

        # Check that we've reached equilibrium
        assert abs(reaction_rates[-1, 0] - 0.0) < 1e-2

    def test_solve_ch3_coverage_dependence(self):
        """
        Test the surface batch reactor can properly apply coverage dependent parameters
        with the nondissociative adsorption of CH3

        Here we choose a kinetic model consisting of the  adsorption reaction
        CH3 + X <=>  CH3X
        We use a sticking coefficient for the rate expression.
        """

        ch3 = Species(
            molecule=[Molecule().from_smiles("[CH3]")],
            thermo=NASA(
                polynomials=[
                    NASAPolynomial(
                        coeffs=[
                            3.91547,
                            0.00184155,
                            3.48741e-06,
                            -3.32746e-09,
                            8.49953e-13,
                            16285.6,
                            0.351743,
                        ],
                        Tmin=(100, "K"),
                        Tmax=(1337.63, "K"),
                    ),
                    NASAPolynomial(
                        coeffs=[
                            3.54146,
                            0.00476786,
                            -1.82148e-06,
                            3.28876e-10,
                            -2.22545e-14,
                            16224,
                            1.66032,
                        ],
                        Tmin=(1337.63, "K"),
                        Tmax=(5000, "K"),
                    ),
                ],
                Tmin=(100, "K"),
                Tmax=(5000, "K"),
                E0=(135.382, "kJ/mol"),
                comment="""Thermo library: primaryThermoLibrary + radical(CH3)""",
            ),
            molecular_weight=(15.0345, "amu"),
        )

        x = Species(
            molecule=[Molecule().from_adjacency_list("1 X u0 p0")],
            thermo=NASA(
                polynomials=[
                    NASAPolynomial(coeffs=[0, 0, 0, 0, 0, 0, 0], Tmin=(298, "K"), Tmax=(1000, "K")),
                    NASAPolynomial(coeffs=[0, 0, 0, 0, 0, 0, 0], Tmin=(1000, "K"), Tmax=(2000, "K")),
                ],
                Tmin=(298, "K"),
                Tmax=(2000, "K"),
                E0=(-6.19426, "kJ/mol"),
                comment="""Thermo library: surfaceThermo""",
            ),
        )

        ch3x = Species(
            molecule=[
                Molecule().from_adjacency_list(
                    """1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 X u0 p0 c0 {1,S}"""
                )
            ],
            thermo=NASA(
                polynomials=[
                    NASAPolynomial(
                        coeffs=[
                            -0.552219,
                            0.026442,
                            -3.55617e-05,
                            2.60044e-08,
                            -7.52707e-12,
                            -4433.47,
                            0.692144,
                        ],
                        Tmin=(298, "K"),
                        Tmax=(1000, "K"),
                    ),
                    NASAPolynomial(
                        coeffs=[
                            3.62557,
                            0.00739512,
                            -2.43797e-06,
                            1.86159e-10,
                            3.6485e-14,
                            -5187.22,
                            -18.9668,
                        ],
                        Tmin=(1000, "K"),
                        Tmax=(2000, "K"),
                    ),
                ],
                Tmin=(298, "K"),
                Tmax=(2000, "K"),
                E0=(-39.1285, "kJ/mol"),
                comment="""Thermo library: surfaceThermoNi111""",
            ),
        )

        rxn1 = Reaction(
            reactants=[ch3, x],
            products=[ch3x],
            kinetics=StickingCoefficient(
                A=0.1,
                n=0,
                Ea=(0, "kcal/mol"),
                T0=(1, "K"),
                Tmin=(200, "K"),
                Tmax=(3000, "K"),
                coverage_dependence={x: {"a": 0.0, "m": -1.0, "E": (0.0, "J/mol")}},
                comment="""Exact match found for rate rule (Adsorbate;VacantSite)""",
            ),
        )
        core_species = [ch3, x, ch3x]
        edge_species = []
        core_reactions = [rxn1]
        edge_reactions = []

        T = 800.0
        P_initial = 1.0e5
        rxn_system = SurfaceReactor(
            T,
            P_initial,
            n_sims=1,
            initial_gas_mole_fractions={ch3: 1.0},
            initial_surface_coverages={x: 1.0},
            surface_volume_ratio=(1.0, "m^-1"),
            surface_site_density=(2.72e-9, "mol/cm^2"),
            coverage_dependence=True,
            termination=[],
        )
        # in chemkin, the sites are mostly occupied in about 1e-8 seconds.

        rxn_system.initialize_model(core_species, core_reactions, edge_species, edge_reactions)

        tlist = np.logspace(-13, -5, 81, dtype=float)

        print("Surface site density:", rxn_system.surface_site_density.value_si)

        print(
            "rxn1 rate coefficient",
            rxn1.get_surface_rate_coefficient(rxn_system.T.value_si, rxn_system.surface_site_density.value_si),
        )

        assert isinstance(rxn1.kinetics.coverage_dependence, dict)  # check to make sure coverage_dependence is still the correct type
        for species, parameters in rxn1.kinetics.coverage_dependence.items():
            assert isinstance(species, Species)  # species should be a Species
            assert isinstance(parameters, dict)
            assert parameters["a"] is not None
            assert parameters["m"] is not None
            assert parameters["E"] is not None

        # Integrate to get the solution at each time point
        t = []
        y = []
        reaction_rates = []
        species_rates = []
        t.append(rxn_system.t)
        # You must make a copy of y because it is overwritten by DASSL at
        # each call to advance()
        y.append(rxn_system.y.copy())
        reaction_rates.append(rxn_system.core_reaction_rates.copy())
        species_rates.append(rxn_system.core_species_rates.copy())
        print("time: ", t)
        print("moles:", y)
        print("reaction rates:", reaction_rates)
        print("species rates:", species_rates)
        start_time = time.time()
        for t1 in tlist:
            rxn_system.advance(t1)
            t.append(rxn_system.t)
            # You must make a copy of y because it is overwritten by DASSL at
            # each call to advance()
            y.append(rxn_system.y.copy())
            reaction_rates.append(rxn_system.core_reaction_rates.copy())
            species_rates.append(rxn_system.core_species_rates.copy())
        run_time = time.time() - start_time
        print(f"Simulation took {run_time:.3e} seconds in {self.id()}")

        # Convert the solution vectors to np arrays
        t = np.array(t, float)
        y = np.array(y, float)
        reaction_rates = np.array(reaction_rates, float)
        species_rates = np.array(species_rates, float)
        V = constants.R * rxn_system.T.value_si * np.sum(y) / rxn_system.P_initial.value_si

        # Check that we're computing the species fluxes correctly
        for i in range(t.shape[0]):
            assert abs(reaction_rates[i, 0] - -species_rates[i, 0]) < 1e-6 * reaction_rates[i, 0]
            assert abs(reaction_rates[i, 0] - -species_rates[i, 1]) < 1e-6 * reaction_rates[i, 0]
            assert abs(reaction_rates[i, 0] - species_rates[i, 2]) < 1e-6 * reaction_rates[i, 0]

        # Check that we've reached equilibrium by the end
        assert abs(reaction_rates[-1, 0] - 0.0) < 1e-2

        # Run model with Covdep off so we can test that it is actually being implemented
        rxn_system = SurfaceReactor(
            T,
            P_initial,
            n_sims=1,
            initial_gas_mole_fractions={ch3: 1.0},
            initial_surface_coverages={x: 1.0},
            surface_volume_ratio=(1.0, "m^-1"),
            surface_site_density=(2.72e-9, "mol/cm^2"),
            termination=[],
        )

        rxn_system.initialize_model(core_species, core_reactions, edge_species, edge_reactions)

        tlist = np.logspace(-13, -5, 81, dtype=float)

        # Integrate to get the solution at each time point
        t = []
        y_off = []
        species_rates_off = []
        t.append(rxn_system.t)

        # You must make a copy of y because it is overwritten by DASSL at
        # each call to advance()
        y_off.append(rxn_system.y.copy())
        species_rates_off.append(rxn_system.core_species_rates.copy())
        for t1 in tlist:
            rxn_system.advance(t1)
            t.append(rxn_system.t)
            # You must make a copy of y because it is overwritten by DASSL at
            # each call to advance()
            y_off.append(rxn_system.y.copy())
            species_rates_off.append(rxn_system.core_species_rates.copy())
        run_time = time.time() - start_time
        print(f"Simulation took {run_time:.3e} seconds in {self.id()}")

        # Convert the solution vectors to np arrays
        t = np.array(t, float)
        y_off = np.array(y_off, float)
        species_rates_off = np.array(species_rates_off, float)

        # Check that we've reached equilibrium
        assert abs(species_rates_off[-1, 0] - 0.0) < 1e-2

        # Check that coverages are different
        assert not np.allclose(y, y_off)
        assert not np.allclose(species_rates, species_rates_off)
