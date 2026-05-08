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


import time

import numpy as np

import rmgpy.constants as constants
from rmgpy.kinetics import SurfaceArrhenius, StickingCoefficient
from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.solver.surface import SurfaceReactor
from rmgpy.species import Species
from rmgpy.thermo import ThermoData, NASA, NASAPolynomial
from rmgpy.solver.termination import TerminationTime
from rmgpy.rmg.settings import ModelSettings, SimulatorSettings


def get_i_thing(thing, thing_list):
    for i in range(len(thing_list)):
        if thing.is_isomorphic(thing_list[i]):
            return i
    return -1

class SurfaceReactorTest:
    def setup_class(self):
        # Define a simple surface mechanism that we can use for testing
        # This is faster than loading in chemkin
        self.species_list = [
            Species(
                molecule=[Molecule().from_smiles("[Ar]")],
                thermo=NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,4.37967], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,4.37967], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), comment="""Thermo library: primaryThermoLibrary"""),
            ),
            Species(
                molecule=[Molecule().from_smiles("[Ne]")],
                thermo=NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,3.35532], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,3.35532], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), comment="""Thermo library: primaryThermoLibrary"""),
            ),
            Species(
                molecule=[Molecule().from_smiles("N#N")],
                thermo=NASA(polynomials=[NASAPolynomial(coeffs=[3.53101,-0.000123661,-5.02999e-07,2.43531e-09,-1.40881e-12,-1046.98,2.96747], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.95258,0.0013969,-4.92632e-07,7.8601e-11,-4.60755e-15,-923.949,5.87189], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), comment="""Thermo library: primaryThermoLibrary"""),
            ),
            Species(
                molecule=[Molecule().from_smiles("C")],
                thermo=NASA(polynomials=[NASAPolynomial(coeffs=[4.20542,-0.00535559,2.51124e-05,-2.13763e-08,5.97526e-12,-10161.9,-0.921283], Tmin=(100,'K'), Tmax=(1084.12,'K')), NASAPolynomial(coeffs=[0.908259,0.0114541,-4.57174e-06,8.29193e-10,-5.66316e-14,-9719.97,13.9931], Tmin=(1084.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), comment="""Thermo library: primaryThermoLibrary"""),
            ),
            Species(
                molecule=[Molecule().from_smiles("[O][O]")],
                thermo=NASA(polynomials=[NASAPolynomial(coeffs=[3.53732,-0.00121572,5.3162e-06,-4.89446e-09,1.45846e-12,-1038.59,4.68368], Tmin=(100,'K'), Tmax=(1074.55,'K')), NASAPolynomial(coeffs=[3.15382,0.00167804,-7.69974e-07,1.51275e-10,-1.08782e-14,-1040.82,6.16756], Tmin=(1074.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), comment="""Thermo library: primaryThermoLibrary"""),
            ),
            Species(
                molecule=[Molecule().from_smiles("O")],
                thermo=NASA(polynomials=[NASAPolynomial(coeffs=[4.05764,-0.000787929,2.90875e-06,-1.47516e-09,2.12833e-13,-30281.6,-0.311362], Tmin=(100,'K'), Tmax=(1130.23,'K')), NASAPolynomial(coeffs=[2.84325,0.00275108,-7.81028e-07,1.07243e-10,-5.79385e-15,-29958.6,5.9104], Tmin=(1130.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), comment="""Thermo library: primaryThermoLibrary"""),
            ),
            Species(
                molecule=[Molecule().from_smiles("[H][H]")],
                thermo=NASA(polynomials=[NASAPolynomial(coeffs=[3.43536,0.000212712,-2.78629e-07,3.4027e-10,-7.76039e-14,-1031.36,-3.90842], Tmin=(100,'K'), Tmax=(1959.07,'K')), NASAPolynomial(coeffs=[2.78819,0.000587616,1.59022e-07,-5.52763e-11,4.34328e-15,-596.156,0.112618], Tmin=(1959.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), comment="""Thermo library: primaryThermoLibrary"""),
            ),
            Species(
                molecule=[Molecule().from_smiles("*")],
                thermo=NASA(polynomials=[NASAPolynomial(coeffs=[0,0,0,0,0,0,0], Tmin=(298,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[0,0,0,0,0,0,0], Tmin=(1000,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), comment="""Thermo library: surfaceThermoPt111"""),
            ),
            Species(
                molecule=[Molecule().from_smiles("[H]*")],
                thermo=NASA(polynomials=[NASAPolynomial(coeffs=[-2.0757,0.0173581,-2.60921e-05,1.89282e-08,-5.38836e-12,-3166.19,8.15362], Tmin=(298,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.72248,-0.00106817,1.98654e-06,-1.12048e-09,2.09812e-13,-4218.24,-15.3207], Tmin=(1000,'K'), Tmax=(2000,'K'))], Tmin=(298,'K'), Tmax=(2000,'K'), comment="""Thermo library: surfaceThermoPt111"""),
            ),
            Species(
                molecule=[Molecule().from_smiles("[O]=*")],
                thermo=NASA(polynomials=[NASAPolynomial(coeffs=[-0.294476,0.0144163,-2.61323e-05,2.19006e-08,-6.98019e-12,-16461.9,-0.199446], Tmin=(298,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.90245,-0.000338584,6.43373e-07,-3.66327e-10,6.90094e-14,-17049.7,-15.256], Tmin=(1000,'K'), Tmax=(2000,'K'))], Tmin=(298,'K'), Tmax=(2000,'K'), comment="""Thermo library: surfaceThermoPt111"""),
            ),
            Species(
                molecule=[Molecule().from_smiles("[CH3]*")],
                thermo=NASA(polynomials=[NASAPolynomial(coeffs=[-0.0444549,0.0194368,-1.91029e-05,1.11269e-08,-2.73736e-12,-6388.04,-0.173376], Tmin=(298,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[8.65705,-0.00790308,1.401e-05,-7.40016e-09,1.31517e-12,-8635.99,-44.3353], Tmin=(1000,'K'), Tmax=(2000,'K'))], Tmin=(298,'K'), Tmax=(2000,'K'), comment="""Thermo library: surfaceThermoPt111"""),
            ),
            Species(
                molecule=[Molecule().from_smiles("[H][H].*")],
                thermo=NASA(polynomials=[NASAPolynomial(coeffs=[3.86406,0.000753456,-1.65571e-06,1.55223e-09,-4.46782e-13,-1689.28,-8.85807], Tmin=(298,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.0688,-0.000495807,6.59234e-07,-1.72598e-10,7.62965e-15,-1700.7,-9.71918], Tmin=(1000,'K'), Tmax=(2000,'K'))], Tmin=(298,'K'), Tmax=(2000,'K'), comment="""Thermo library: surfaceThermoPt111"""),
            ),
            Species(
                molecule=[Molecule().from_smiles("[OH]*")],
                thermo=NASA(polynomials=[NASAPolynomial(coeffs=[1.4236,0.015783,-2.91659e-05,2.50433e-08,-8.04088e-12,-18999.3,-3.1523], Tmin=(298,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.03574,-0.00134422,2.25916e-06,-1.08548e-09,1.77876e-13,-19634.4,-20.0345], Tmin=(1000,'K'), Tmax=(2000,'K'))], Tmin=(298,'K'), Tmax=(2000,'K'), comment="""Thermo library: surfaceThermoPt111"""),
            ),
            Species(
                molecule=[Molecule().from_smiles("O.*")],
                thermo=NASA(polynomials=[NASAPolynomial(coeffs=[2.72971,0.00871052,-1.29132e-05,1.07295e-08,-3.39434e-12,-32612.7,-6.0448], Tmin=(298,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.85496,-0.00328847,5.56991e-06,-2.73008e-09,4.55898e-13,-33304.6,-21.3518], Tmin=(1000,'K'), Tmax=(2000,'K'))], Tmin=(298,'K'), Tmax=(2000,'K'), comment="""Thermo library: surfaceThermoPt111"""),
            ),
            Species(
                molecule=[Molecule().from_smiles("[CH2]=*")],
                thermo=NASA(polynomials=[NASAPolynomial(coeffs=[-2.23007,0.0292223,-4.33155e-05,3.31428e-08,-9.96471e-12,-222.256,8.30173], Tmin=(298,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[6.8346,-0.00514926,9.15491e-06,-4.84917e-09,8.63767e-13,-2258.98,-36.2215], Tmin=(1000,'K'), Tmax=(2000,'K'))], Tmin=(298,'K'), Tmax=(2000,'K'), comment="""Thermo library: surfaceThermoPt111"""),
            ),
            Species(
                molecule=[Molecule().from_smiles("[CH]#*")],
                thermo=NASA(polynomials=[NASAPolynomial(coeffs=[-2.66805,0.0290693,-4.82654e-05,3.87589e-08,-1.19749e-11,-2918.16,9.72941], Tmin=(298,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.9043,-0.00263865,4.71729e-06,-2.51267e-09,4.49659e-13,-4464.41,-26.7108], Tmin=(1000,'K'), Tmax=(2000,'K'))], Tmin=(298,'K'), Tmax=(2000,'K'), comment="""Thermo library: surfaceThermoPt111"""),
            ),
            Species(
                molecule=[Molecule().from_smiles("[C]$*")],
                thermo=NASA(polynomials=[NASAPolynomial(coeffs=[-1.94351,0.0197767,-3.36337e-05,2.69027e-08,-8.27959e-12,7000.57,7.1747], Tmin=(298,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.81347,-0.000693952,1.30308e-06,-7.387e-10,1.38796e-13,6060.03,-15.5738], Tmin=(1000,'K'), Tmax=(2000,'K'))], Tmin=(298,'K'), Tmax=(2000,'K'), comment="""Thermo library: surfaceThermoPt111"""),
            )
        ]
        X = self.species_list[get_i_thing(Species(molecule=[Molecule().from_smiles("*")]), self.species_list)]
        O2 = self.species_list[get_i_thing(Species(molecule=[Molecule().from_smiles("[O][O]")]), self.species_list)]
        OX = self.species_list[get_i_thing(Species(molecule=[Molecule().from_smiles("O=*")]), self.species_list)]
        H2 = self.species_list[get_i_thing(Species(molecule=[Molecule().from_smiles("[H][H]")]), self.species_list)]
        HX = self.species_list[get_i_thing(Species(molecule=[Molecule().from_smiles("[H]*")]), self.species_list)]
        CH4 = self.species_list[get_i_thing(Species(molecule=[Molecule().from_smiles("C")]), self.species_list)]
        CH3X = self.species_list[get_i_thing(Species(molecule=[Molecule().from_smiles("C*")]), self.species_list)]
        OHX = self.species_list[get_i_thing(Species(molecule=[Molecule().from_smiles("[OH]*")]), self.species_list)]
        H2O = self.species_list[get_i_thing(Species(molecule=[Molecule().from_smiles("O")]), self.species_list)]
        H2OX = self.species_list[get_i_thing(Species(molecule=[Molecule().from_smiles("O.*")]), self.species_list)]
        CHX = self.species_list[get_i_thing(Species(molecule=[Molecule().from_smiles("[CH]#*")]), self.species_list)]
        CH2X = self.species_list[get_i_thing(Species(molecule=[Molecule().from_smiles("[CH2]=*")]), self.species_list)]

        self.reaction_list = [
            Reaction(
                reactants=[X, X, O2],
                products=[OX, OX],
                kinetics=StickingCoefficient(A=0.07, n=0, Ea=(0,'kcal/mol'), T0=(1,'K'))
            ),
            Reaction(
                reactants=[X, X, H2],
                products=[HX, HX],
                kinetics=StickingCoefficient(A=0.046, n=0, Ea=(0,'kcal/mol'), T0=(1,'K'))
            ),
            Reaction(
                reactants=[HX, CH3X],
                products=[X, X, CH4],
                kinetics=SurfaceArrhenius(A=(3.3e+21,'cm^2/(mol*s)'), n=0, Ea=(11.95,'kcal/mol'), T0=(1,'K'))
            ),
            Reaction(
                reactants=[X, OHX],
                products=[OX, HX],
                kinetics=SurfaceArrhenius(A=(7.39e+19,'cm^2/(mol*s)'), n=0, Ea=(18.475,'kcal/mol'), T0=(1,'K'))
            ),
            Reaction(
                reactants=[X, H2O],
                products=[H2OX],
                kinetics=StickingCoefficient(A=0.75, n=0, Ea=(0,'kcal/mol'), T0=(1,'K'))
            ),
            Reaction(
                reactants=[X, H2OX],
                products=[HX, OHX],
                kinetics=SurfaceArrhenius(A=(1.15e+19,'cm^2/(mol*s)'), n=0, Ea=(24.235,'kcal/mol'), T0=(1,'K'))
            ),
            Reaction(
                reactants=[OHX, CHX],
                products=[OX, CH2X],
                kinetics=SurfaceArrhenius(A=(4.4e+22,'cm^2/(mol*s)'), n=0.101, Ea=(10.143,'kcal/mol'), T0=(1,'K'))
            ),
            Reaction(
                reactants=[OHX, CH2X],
                products=[OX, CH3X],
                kinetics=SurfaceArrhenius(A=(1.39e+21,'cm^2/(mol*s)'), n=0.101, Ea=(4.541,'kcal/mol'), T0=(1,'K'))
            ),
            Reaction(
                reactants=[X, OX, CH4],
                products=[OHX, CH3X],
                kinetics=SurfaceArrhenius(A=(5e+18,'cm^4/(mol^2*s)'), n=0.7, Ea=(10.038,'kcal/mol'), T0=(1,'K'))
            )
        ]

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
            assert abs(reaction_rates[i, 0] - -1.0 * species_rates[i, 0]) < abs(
                1e-6 * reaction_rates[i, 0]
            )
            assert abs(reaction_rates[i, 0] - -0.5 * species_rates[i, 1]) < abs(
                1e-6 * reaction_rates[i, 0]
            )
            assert abs(reaction_rates[i, 0] - 0.5 * species_rates[i, 2]) < abs(
                1e-6 * reaction_rates[i, 0]
            )

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
            assert abs(reaction_rates[i, 0] - -species_rates[i, 0]) < abs(
                1e-6 * reaction_rates[i, 0]
            )
            assert abs(reaction_rates[i, 0] - -species_rates[i, 1]) < abs(
                1e-6 * reaction_rates[i, 0]
            )
            assert abs(reaction_rates[i, 0] - species_rates[i, 2]) < abs(
                1e-6 * reaction_rates[i, 0]
            )

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
        print(f"Simulation took {run_time:.3e} seconds")

        # Convert the solution vectors to np arrays
        t = np.array(t, float)
        y = np.array(y, float)
        reaction_rates = np.array(reaction_rates, float)
        species_rates = np.array(species_rates, float)
        total_sites = y[0, 1]

        # Check that we're computing the species fluxes correctly
        for i in range(t.shape[0]):
            assert abs(reaction_rates[i, 0] - -1.0 * species_rates[i, 0]) < abs(
                1e-6 * reaction_rates[i, 0]
            )
            assert abs(reaction_rates[i, 0] - -0.5 * species_rates[i, 1]) < abs(
                1e-6 * reaction_rates[i, 0]
            )
            assert abs(reaction_rates[i, 0] - 0.5 * species_rates[i, 2]) < abs(
                1e-6 * reaction_rates[i, 0]
            )

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
        print(f"Simulation took {run_time:.3e} seconds")

        # Convert the solution vectors to np arrays
        t = np.array(t, float)
        y = np.array(y, float)
        reaction_rates = np.array(reaction_rates, float)
        species_rates = np.array(species_rates, float)
        V = constants.R * rxn_system.T.value_si * np.sum(y) / rxn_system.P_initial.value_si

        # Check that we're computing the species fluxes correctly
        for i in range(t.shape[0]):
            assert abs(reaction_rates[i, 0] - -species_rates[i, 0]) < abs(
                1e-6 * reaction_rates[i, 0]
            )
            assert abs(reaction_rates[i, 0] - -species_rates[i, 1]) < abs(
                1e-6 * reaction_rates[i, 0]
            )
            assert abs(reaction_rates[i, 0] - species_rates[i, 2]) < abs(
                1e-6 * reaction_rates[i, 0]
            )

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
        print(f"Simulation took {run_time:.3e} seconds")

        # Convert the solution vectors to np arrays
        t = np.array(t, float)
        y_off = np.array(y_off, float)
        species_rates_off = np.array(species_rates_off, float)

        # Check that we've reached equilibrium
        assert abs(species_rates_off[-1, 0] - 0.0) < 1e-2

        # Check that coverages are different
        assert not np.allclose(y, y_off)
        assert not np.allclose(species_rates, species_rates_off)

    def test_solve_h2_thermo_coverage_dependence(self):
        """
        Test the surface batch reactor can properly apply thermodynamic coverage dependent parameters
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
                thermo_coverage_dependence={"1 H u0 p0 {2,S} \n 2 X u0 p0 {1,S}":{'model':'polynomial', 'enthalpy-coefficients':[(1,'J/mol'),(2,'J/mol'),(3,'J/mol')], "entropy-coefficients":[(1,'J/(mol*K)'),(5,'J/(mol*K)'),(3,'J/(mol*K)')]},}
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
            ),
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
            coverage_dependence=True,
            thermo_coverage_dependence=True,
            termination=[],
        )

        rxn_system.initialize_model(core_species, core_reactions, edge_species, edge_reactions)

        tlist = np.logspace(-13, -5, 81, dtype=float)

        assert isinstance(hx.thermo.thermo_coverage_dependence, dict)  # check to make sure coverage_dependence is still the correct type
        for species, parameters in hx.thermo.thermo_coverage_dependence.items():
            assert isinstance(species, str)  # species should be an ajacency list
            assert isinstance(parameters, dict)
            assert parameters["model"] is not None
            assert parameters["enthalpy-coefficients"] is not None
            assert parameters["entropy-coefficients"] is not None
        assert np.array_equal(rxn_system.stoi_matrix, np.array([[-1., -2.,  2.]]))
        thermo_coeffs = np.array([np.zeros((3,6))]*3)
        thermo_coeffs[-1][-1] = [1., 2., 3., 1., 5., 3.]
        assert np.array_equal(rxn_system.thermo_coeff_matrix, thermo_coeffs)
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
        print(f"Simulation took {run_time:.3e} seconds")

        # Convert the solution vectors to np arrays
        t = np.array(t, float)
        y = np.array(y, float)
        reaction_rates = np.array(reaction_rates, float)
        species_rates = np.array(species_rates, float)
        total_sites = y[0, 1]

        # Check that we're computing the species fluxes correctly
        for i in range(t.shape[0]):
            assert abs(reaction_rates[i, 0] - -1.0 * species_rates[i, 0]) < abs(
                1e-6 * reaction_rates[i, 0]
            )
            assert abs(reaction_rates[i, 0] - -0.5 * species_rates[i, 1]) < abs(
                1e-6 * reaction_rates[i, 0]
            )
            assert abs(reaction_rates[i, 0] - 0.5 * species_rates[i, 2]) < abs(
                1e-6 * reaction_rates[i, 0]
            )

        # Check that we've reached equilibrium
        assert abs(reaction_rates[-1, 0] - 0.0) < 1e-2

    def test_solve_ch3_thermo_coverage_dependence(self):
        """
        Test the surface batch reactor can properly apply thermodynamic coverage dependent parameters
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
                thermo_coverage_dependence={"1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S} \n 2 H u0 p0 c0 {1,S} \n 3 H u0 p0 c0 {1,S} \n 4 H u0 p0 c0 {1,S} \n 5 X u0 p0 c0 {1,S}":
                                            {'model':'polynomial', 'enthalpy-coefficients':[(1e5,'J/mol'),(2,'J/mol'),(3,'J/mol')], "entropy-coefficients":[(1,'J/(mol*K)'),(5,'J/(mol*K)'),(3,'J/(mol*K)')]},},
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
            thermo_coverage_dependence=True,
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

        assert isinstance(ch3x.thermo.thermo_coverage_dependence, dict)  # check to make sure coverage_dependence is still the correct type
        for species, parameters in ch3x.thermo.thermo_coverage_dependence.items():
            assert isinstance(species, str)  # species should be an ajacency list
            assert isinstance(parameters, dict)
            assert parameters["model"] is not None
            assert parameters["enthalpy-coefficients"] is not None
            assert parameters["entropy-coefficients"] is not None
        
        # check thermo_coverage_dependence is on
        # and thermo_coeff_matrix and stoi_matrix are correctly created
        assert rxn_system.thermo_coverage_dependence is True
        assert np.array_equal(rxn_system.stoi_matrix, np.array([[-1, -1, 1]], dtype=float))
        thermo_coeff_matrix = np.array([np.zeros((3,6))]*3)
        thermo_coeff_matrix[-1][-1] = [1e5, 2, 3, 1, 5, 3]
        assert np.array_equal(rxn_system.thermo_coeff_matrix, thermo_coeff_matrix)

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
        print(f"Simulation took {run_time:.3e} seconds")

        # Convert the solution vectors to np arrays
        t = np.array(t, float)
        y = np.array(y, float)
        reaction_rates = np.array(reaction_rates, float)
        species_rates = np.array(species_rates, float)
        V = constants.R * rxn_system.T.value_si * np.sum(y) / rxn_system.P_initial.value_si

        # Check that we're computing the species fluxes correctly
        for i in range(t.shape[0]):
            assert abs(reaction_rates[i, 0] - -species_rates[i, 0]) < abs(
                1e-6 * reaction_rates[i, 0]
            )
            assert abs(reaction_rates[i, 0] - -species_rates[i, 1]) < abs(
                1e-6 * reaction_rates[i, 0]
            )
            assert abs(reaction_rates[i, 0] - species_rates[i, 2]) < abs(
                1e-6 * reaction_rates[i, 0]
            )

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
        
        # check thermo_coverage_dependence is off
        # and thermo_coeff_matrix and stoi_matrix are not created
        assert rxn_system.thermo_coverage_dependence is False
        assert rxn_system.thermo_coeff_matrix is None
        assert rxn_system.stoi_matrix is None

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
        print(f"Simulation took {run_time:.3e} seconds")

        # Convert the solution vectors to np arrays
        t = np.array(t, float)
        y_off = np.array(y_off, float)
        species_rates_off = np.array(species_rates_off, float)

        # Check that we've reached equilibrium
        assert abs(species_rates_off[-1, 0] - 0.0) < 1e-2

        # Check that coverages are different
        assert not np.allclose(y, y_off)
        assert not np.allclose(species_rates, species_rates_off)

    def test_jacobian_with_finite_differences(self):

        # make a simple mechanism
        x_O2 = 3e-5
        x_CH4 = 1e-5
        T = (900, 'K')
        P = (1, 'atm')

        surface_volume_ratio = (1.0, "m^-1")  # TODO, try higher surface volume ratios
        surface_site_density = (2.483e-8, "kmol/m^2")  # read from Cantera yaml just to be sure these match
        termination = TerminationTime((1.0, 's'))


        Ar = self.species_list[get_i_thing(Species(smiles='[Ar]'), self.species_list)]
        CH4 = self.species_list[get_i_thing(Species(smiles='C'), self.species_list)]
        CO2  = self.species_list[get_i_thing(Species(smiles='O=C=O'), self.species_list)]
        O2 = self.species_list[get_i_thing(Species(smiles='[O][O]'), self.species_list)]
        X = self.species_list[get_i_thing(Species(smiles='*'), self.species_list)]

        initial_gas_mole_fractions = {O2: x_O2, CH4: x_CH4, Ar: 1.0 - x_O2 - x_CH4}
        initial_surface_coverages = {X: 1.0}
        sensitive_species = [CO2]

        reaction_system = SurfaceReactor(
            T,
            P,
            n_sims=1,
            initial_gas_mole_fractions=initial_gas_mole_fractions,
            initial_surface_coverages=initial_surface_coverages,
            surface_volume_ratio=surface_volume_ratio,
            surface_site_density=surface_site_density,
            termination=[termination],
            sensitive_species=sensitive_species,
        )

        # Need to list some sens_worksheets so this doesn't crash when sensitivity is turned on
        reaction_system.simulate(
            core_species=self.species_list,
            core_reactions=self.reaction_list,
            edge_species=[],
            edge_reactions=[],
            surface_species=[],
            surface_reactions=[],
            model_settings=ModelSettings(tol_move_to_core=1e5),  # tol_move_to_core isn't set by default which causes an error
            simulator_settings=SimulatorSettings(),  # defaults
            sensitivity=True,
            sens_worksheet=['temp_sensitivity.csv'],
        )

        # now comapare the jacobian to a finite difference approximation of the jacobian for the sensitive species
        t = reaction_system.t
        y = reaction_system.y
        dydt = reaction_system.dydt

        J_analytical = reaction_system.jacobian(t, y, dydt, cj=0.0)

        # turn off sensitivity so that it only handles the first n values of y and doesn't get bogged
        # down in all the partial derivatives that get stored in y for sensitivity calculations
        reaction_system.sensitivity = False

        n = reaction_system.num_core_species
        eps_machine = np.sqrt(np.finfo(float).eps)
        n_scale = 1e-10

        # Get baseline residual f(n)
        f0, _ = reaction_system.residual(t, y[:n], dydt[:n])

        # Build finite difference Jacobian by column
        J_fd = np.zeros((n, n))
        for s in range(n):
            y_perturbed = y.copy()
            delta_n_s = eps_machine * max(abs(y[s]), n_scale)  # figure out how much to perturb the value
            y_perturbed[s] += delta_n_s  # add \Delta n_s
            
            f_perturbed, _ = reaction_system.residual(t, y_perturbed[:n], dydt[:n])
            J_fd[:, s] = (f_perturbed - f0) / delta_n_s

        assert np.allclose(J_analytical, J_fd, rtol=1e-2, atol=1e-6)
