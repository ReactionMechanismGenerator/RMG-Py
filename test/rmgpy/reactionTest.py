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

"""
This module contains unit tests of the rmgpy.reaction module.
"""


import cantera as ct
import numpy
import yaml
from copy import deepcopy

import pytest

import rmgpy.constants as constants
from rmgpy.kinetics import (
    Arrhenius,
    ArrheniusEP,
    MultiArrhenius,
    PDepArrhenius,
    MultiPDepArrhenius,
    ThirdBody,
    Troe,
    Lindemann,
    Chebyshev,
    SurfaceArrhenius,
    StickingCoefficient,
)
from rmgpy.molecule import Molecule
from rmgpy.quantity import Quantity
from rmgpy.reaction import Reaction
from rmgpy.species import Species, TransitionState
from rmgpy.statmech.conformer import Conformer
from rmgpy.statmech.rotation import NonlinearRotor
from rmgpy.statmech.torsion import HinderedRotor
from rmgpy.statmech.translation import IdealGasTranslation
from rmgpy.statmech.vibration import HarmonicOscillator
from rmgpy.thermo import Wilhoit, ThermoData, NASA, NASAPolynomial


class PseudoSpecies(object):
    """
    Can be used in place of a :class:`rmg.species.Species` for isomorphism checks.

    PseudoSpecies('a') is isomorphic with PseudoSpecies('A')
    but nothing else.
    """

    def __init__(self, label):
        self.label = label

    def __repr__(self):
        return "PseudoSpecies('{0}')".format(self.label)

    def __str__(self):
        return self.label

    def is_isomorphic(self, other, generate_initial_map=False, strict=True, save_order=False):
        return self.label.lower() == other.label.lower()


class TestReactionIsomorphism:
    """
    Contains unit tests of the isomorphism testing of the  Reaction class.
    """

    def make_reaction(self, reaction_string):
        """ "
        Make a Reaction (containing PseudoSpecies) of from a string like 'Ab=CD'
        """
        reactants, products = reaction_string.split("=")
        reactants = [PseudoSpecies(i) for i in reactants]
        products = [PseudoSpecies(i) for i in products]
        return Reaction(reactants=reactants, products=products)

    def test1to1(self):
        r1 = self.make_reaction("A=B")
        assert r1.is_isomorphic(self.make_reaction("a=B"))
        assert r1.is_isomorphic(self.make_reaction("b=A"))
        assert not r1.is_isomorphic(self.make_reaction("B=a"), either_direction=False)
        assert not r1.is_isomorphic(self.make_reaction("A=C"))
        assert not r1.is_isomorphic(self.make_reaction("A=BB"))

    def test1to2(self):
        r1 = self.make_reaction("A=BC")
        assert r1.is_isomorphic(self.make_reaction("a=Bc"))
        assert r1.is_isomorphic(self.make_reaction("cb=a"))
        assert r1.is_isomorphic(self.make_reaction("a=cb"), either_direction=False)
        assert not r1.is_isomorphic(self.make_reaction("bc=a"), either_direction=False)
        assert not r1.is_isomorphic(self.make_reaction("a=c"))
        assert not r1.is_isomorphic(self.make_reaction("ab=c"))

    def test2to2(self):
        r1 = self.make_reaction("AB=CD")
        assert r1.is_isomorphic(self.make_reaction("ab=cd"))
        assert r1.is_isomorphic(self.make_reaction("ab=dc"), either_direction=False)
        assert r1.is_isomorphic(self.make_reaction("dc=ba"))
        assert not r1.is_isomorphic(self.make_reaction("cd=ab"), either_direction=False)
        assert not r1.is_isomorphic(self.make_reaction("ab=ab"))
        assert not r1.is_isomorphic(self.make_reaction("ab=cde"))

    def test2to3(self):
        r1 = self.make_reaction("AB=CDE")
        assert r1.is_isomorphic(self.make_reaction("ab=cde"))
        assert r1.is_isomorphic(self.make_reaction("ba=edc"), either_direction=False)
        assert r1.is_isomorphic(self.make_reaction("dec=ba"))
        assert not r1.is_isomorphic(self.make_reaction("cde=ab"), either_direction=False)
        assert not r1.is_isomorphic(self.make_reaction("ab=abc"))
        assert not r1.is_isomorphic(self.make_reaction("abe=cde"))

    def test2to3_using_check_only_label(self):
        r1 = self.make_reaction("AB=CDE")
        assert r1.is_isomorphic(self.make_reaction("AB=CDE"), check_only_label=True)
        assert r1.is_isomorphic(
            self.make_reaction("BA=EDC"),
            either_direction=False,
            check_only_label=True,
        )
        assert not r1.is_isomorphic(self.make_reaction("Ab=CDE"), check_only_label=True)
        assert not r1.is_isomorphic(
            self.make_reaction("BA=EDd"),
            either_direction=False,
            check_only_label=True,
        )


class TestSurfaceReaction:
    """Test surface reactions"""

    def setup_class(self):
        m_h2 = Molecule().from_smiles("[H][H]")
        m_x = Molecule().from_adjacency_list("1 X u0 p0")
        m_hx = Molecule().from_smiles("[H][*]")
        # m_hx = Molecule().from_adjacency_list("1 H u0 p0 {2,S} \n 2 X u0 p0 {1,S}")
        m_ch3 = Molecule().from_smiles("[CH3]")
        m_ch3x = Molecule().from_adjacency_list(
            """1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 X u0 p0 c0 {1,S}"""
        )

        s_h2 = Species(
            label="H2(1)",
            molecule=[m_h2],
            thermo=ThermoData(
                Tdata=([300, 400, 500, 600, 800, 1000, 1500, 2000], "K"),
                Cpdata=(
                    [6.955, 6.955, 6.956, 6.961, 7.003, 7.103, 7.502, 8.17],
                    "cal/(mol*K)",
                ),
                H298=(0, "kcal/mol"),
                S298=(31.129, "cal/(mol*K)"),
            ),
        )
        s_x = Species(
            label="X(2)",
            molecule=[m_x],
            thermo=ThermoData(
                Tdata=([300, 400, 500, 600, 800, 1000, 1500, 2000], "K"),
                Cpdata=([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], "cal/(mol*K)"),
                H298=(0.0, "kcal/mol"),
                S298=(0.0, "cal/(mol*K)"),
            ),
        )
        s_hx = Species(
            label="HX(3)",
            molecule=[m_hx],
            thermo=ThermoData(
                Tdata=([300, 400, 500, 600, 800, 1000, 1500, 2000], "K"),
                Cpdata=(
                    [1.50, 2.58, 3.40, 4.00, 4.73, 5.13, 5.57, 5.82],
                    "cal/(mol*K)",
                ),
                H298=(-11.26, "kcal/mol"),
                S298=(0.44, "cal/(mol*K)"),
            ),
        )

        s_ch3 = Species(
            label="CH3(4)",
            molecule=[m_ch3],
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

        s_ch3x = Species(
            label="CH3X(5)",
            molecule=[m_ch3x],
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

        rxn1s = Reaction(
            reactants=[s_h2, s_x, s_x],
            products=[s_hx, s_hx],
            kinetics=SurfaceArrhenius(A=(9.05e18, "cm^5/(mol^2*s)"), n=0.5, Ea=(5.0, "kJ/mol"), T0=(1.0, "K")),
        )
        self.rxn1s = rxn1s

        rxn1m = Reaction(reactants=[m_h2, m_x, m_x], products=[m_hx, m_hx])
        self.rxn1m = rxn1m

        self.rxn2sSC = Reaction(
            reactants=[s_ch3, s_x],
            products=[s_ch3x],
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
        self.rxn2sSA = Reaction(
            reactants=[s_ch3, s_x],
            products=[s_ch3x],
            kinetics=SurfaceArrhenius(
                A=(2.7e10, "cm^3/(mol*s)"),
                n=0.5,
                Ea=(5.0, "kJ/mol"),
                T0=(1.0, "K"),
                comment="""Approximate rate""",
            ),
        )

        # adding coverage dependent reaction
        self.rxn_covdep = Reaction(
            reactants=[s_h2, s_x, s_x],
            products=[s_hx, s_hx],
            kinetics=SurfaceArrhenius(
                A=(9.05e18, "cm^5/(mol^2*s)"),
                n=0.5,
                Ea=(5.0, "kJ/mol"),
                T0=(1.0, "K"),
                coverage_dependence={
                    Species().from_adjacency_list("1 X u0 p0 c0"): {
                        "E": (0.1, "J/mol"),
                        "m": -1.0,
                        "a": 1.0,
                    },
                    Species().from_adjacency_list("1 X u0 p0 c0 {2,S}\n2 H u0 p0 c0 {1,S}"): {"E": (0.2, "J/mol"), "m": -2.0, "a": 2.0},
                },
            ),
        )

    def test_is_surface_reaction_species(self):
        """Test is_surface_reaction for reaction based on Species"""
        assert self.rxn1s.is_surface_reaction()

    def test_is_surface_reaction_molecules(self):
        """Test is_surface_reaction for reaction based on Molecules"""
        assert self.rxn1m.is_surface_reaction()

    def test_methyl_adsorption_surface_arrhenius(self):
        """Test the CH3 adsorption rate given by SurfaceArrhenius"""
        T = 800
        surface_site_density = Quantity(2.72e-9, "mol/cm^2").value_si
        calculated = self.rxn2sSA.get_surface_rate_coefficient(T, surface_site_density)
        target = 1e6  # mol/m2
        assert round(abs(numpy.log10(calculated) - numpy.log10(target)), 0) == 0

    def test_methyl_adsorption_sticking_coefficient(self):
        """Test the CH3 adsorption rate given by StickingCoefficient"""

        # First, check the molecular weight is in units we expect
        assert round(abs(self.rxn2sSC.reactants[0].molecular_weight.value_si / constants.amu / 1000 - 15.0345e-3), 7) == 0  # kg/mol

        T = 800
        surface_site_density = Quantity(2.72e-9, "mol/cm^2").value_si
        calculated = self.rxn2sSC.get_surface_rate_coefficient(T, surface_site_density)
        target = 1e6  # mol/m2
        assert round(abs(numpy.log10(calculated) - numpy.log10(target)), 0) == 0

    def test_get_rate_coefficient_units_from_reaction_order(self):
        assert self.rxn1s.generate_reverse_rate_coefficient().A.units == "m^2/(mol*s)"

    def test_equilibrium_constant_surface_kc(self):
        """
        Test the Reaction.get_equilibrium_constant() method for Kc of a surface reaction.
        """
        Tlist = numpy.arange(400.0, 1600.0, 200.0, numpy.float64)
        Kclist0 = [15375.20186, 1.566753, 0.017772, 0.0013485, 0.000263180, 8.73504e-05]
        Kclist = self.rxn1s.get_equilibrium_constants(Tlist, type="Kc")
        for i in range(len(Tlist)):
            assert round(abs(Kclist[i] / Kclist0[i] - 1.0), 4) == 0

    def test_reverse_sticking_coeff_rate(self):
        """
        Test the Reaction.reverse_sticking_coeff_rate() method works for StickingCoefficient format.
        """

        original_kinetics = self.rxn2sSC.kinetics
        rxn_copy = deepcopy(self.rxn2sSC)
        reverse_kinetics = self.rxn2sSC.generate_reverse_rate_coefficient(surface_site_density=2.5e-5)

        rxn_copy.kinetics = reverse_kinetics
        # reverse reactants, products to ensure Keq is correctly computed
        rxn_copy.reactants, rxn_copy.products = rxn_copy.products, rxn_copy.reactants
        reverse_reverse_kinetics = rxn_copy.generate_reverse_rate_coefficient(surface_site_density=2.5e-5)
        rxn_copy.kinetics = reverse_reverse_kinetics

        # check that reverting the reverse yields the original
        Tlist = numpy.arange(
            original_kinetics.Tmin.value_si,
            original_kinetics.Tmax.value_si,
            200.0,
            numpy.float64,
        )
        P = 0
        for T in Tlist:
            korig = self.rxn2sSC.get_rate_coefficient(T, P, surface_site_density=2.5e-5)
            krevrev = rxn_copy.get_rate_coefficient(T, P, surface_site_density=2.5e-5)
            assert round(abs(korig / krevrev - 1.0), 0) == 0


class TestReaction:
    """
    Contains unit tests of the Reaction class.
    """

    def setup_class(self):
        """
        A method that is called prior to each unit test in this class.
        """
        ethylene = Species(
            label="C2H4",
            conformer=Conformer(
                E0=(44.7127, "kJ/mol"),
                modes=[
                    IdealGasTranslation(
                        mass=(28.0313, "amu"),
                    ),
                    NonlinearRotor(
                        inertia=(
                            [3.41526, 16.6498, 20.065],
                            "amu*angstrom^2",
                        ),
                        symmetry=4,
                    ),
                    HarmonicOscillator(
                        frequencies=(
                            [
                                828.397,
                                970.652,
                                977.223,
                                1052.93,
                                1233.55,
                                1367.56,
                                1465.09,
                                1672.25,
                                3098.46,
                                3111.7,
                                3165.79,
                                3193.54,
                            ],
                            "cm^-1",
                        ),
                    ),
                ],
                spin_multiplicity=1,
                optical_isomers=1,
            ),
        )

        hydrogen = Species(
            label="H",
            conformer=Conformer(
                E0=(211.794, "kJ/mol"),
                modes=[
                    IdealGasTranslation(
                        mass=(1.00783, "amu"),
                    ),
                ],
                spin_multiplicity=2,
                optical_isomers=1,
            ),
        )

        ethyl = Species(
            label="C2H5",
            conformer=Conformer(
                E0=(111.603, "kJ/mol"),
                modes=[
                    IdealGasTranslation(
                        mass=(29.0391, "amu"),
                    ),
                    NonlinearRotor(
                        inertia=(
                            [4.8709, 22.2353, 23.9925],
                            "amu*angstrom^2",
                        ),
                        symmetry=1,
                    ),
                    HarmonicOscillator(
                        frequencies=(
                            [
                                482.224,
                                791.876,
                                974.355,
                                1051.48,
                                1183.21,
                                1361.36,
                                1448.65,
                                1455.07,
                                1465.48,
                                2688.22,
                                2954.51,
                                3033.39,
                                3101.54,
                                3204.73,
                            ],
                            "cm^-1",
                        ),
                    ),
                    HinderedRotor(
                        inertia=(1.11481, "amu*angstrom^2"),
                        symmetry=6,
                        barrier=(0.244029, "kJ/mol"),
                        semiclassical=None,
                    ),
                ],
                spin_multiplicity=2,
                optical_isomers=1,
            ),
        )

        TS = TransitionState(
            label="TS",
            conformer=Conformer(
                E0=(266.694, "kJ/mol"),
                modes=[
                    IdealGasTranslation(
                        mass=(29.0391, "amu"),
                    ),
                    NonlinearRotor(
                        inertia=(
                            [6.78512, 22.1437, 22.2114],
                            "amu*angstrom^2",
                        ),
                        symmetry=1,
                    ),
                    HarmonicOscillator(
                        frequencies=(
                            [
                                412.75,
                                415.206,
                                821.495,
                                924.44,
                                982.714,
                                1024.16,
                                1224.21,
                                1326.36,
                                1455.06,
                                1600.35,
                                3101.46,
                                3110.55,
                                3175.34,
                                3201.88,
                            ],
                            "cm^-1",
                        ),
                    ),
                ],
                spin_multiplicity=2,
                optical_isomers=1,
            ),
            frequency=(-750.232, "cm^-1"),
        )

        self.reaction = Reaction(
            reactants=[hydrogen, ethylene],
            products=[ethyl],
            kinetics=Arrhenius(
                A=(501366000.0, "cm^3/(mol*s)"),
                n=1.637,
                Ea=(4.32508, "kJ/mol"),
                T0=(1, "K"),
                Tmin=(300, "K"),
                Tmax=(2500, "K"),
            ),
            transition_state=TS,
            degeneracy=2,
        )
        self.reaction.kinetics.comment = """
        Multiplied by reaction path degeneracy 2.0
        """

        # CC(=O)O[O]
        acetylperoxy = Species(
            label="acetylperoxy",
            molecule=[Molecule(smiles="CC(=O)O[O]")],
            thermo=Wilhoit(
                Cp0=(4.0 * constants.R, "J/(mol*K)"),
                CpInf=(21.0 * constants.R, "J/(mol*K)"),
                a0=-3.95,
                a1=9.26,
                a2=-15.6,
                a3=8.55,
                B=(500.0, "K"),
                H0=(-6.151e04, "J/mol"),
                S0=(-790.2, "J/(mol*K)"),
            ),
        )

        # C[C]=O
        acetyl = Species(
            label="acetyl",
            molecule=[Molecule(smiles="C[C]=O")],
            thermo=Wilhoit(
                Cp0=(4.0 * constants.R, "J/(mol*K)"),
                CpInf=(15.5 * constants.R, "J/(mol*K)"),
                a0=0.2541,
                a1=-0.4712,
                a2=-4.434,
                a3=2.25,
                B=(500.0, "K"),
                H0=(-1.439e05, "J/mol"),
                S0=(-524.6, "J/(mol*K)"),
            ),
        )

        # [O][O]
        oxygen = Species(
            label="oxygen",
            molecule=[Molecule(smiles="[O][O]")],
            thermo=Wilhoit(
                Cp0=(3.5 * constants.R, "J/(mol*K)"),
                CpInf=(4.5 * constants.R, "J/(mol*K)"),
                a0=-0.9324,
                a1=26.18,
                a2=-70.47,
                a3=44.12,
                B=(500.0, "K"),
                H0=(1.453e04, "J/mol"),
                S0=(-12.19, "J/(mol*K)"),
            ),
        )

        self.reaction2 = Reaction(
            reactants=[acetyl, oxygen],
            products=[acetylperoxy],
            kinetics=Arrhenius(
                A=(2.65e12, "cm^3/(mol*s)"),
                n=0.0,
                Ea=(0.0, "kJ/mol"),
                T0=(1, "K"),
                Tmin=(300, "K"),
                Tmax=(2000, "K"),
            ),
        )

        oxygen_atom = Species().from_smiles("[O]")
        so2 = Species().from_smiles("O=S=O")
        so3 = Species().from_smiles("O=S(=O)=O")

        self.reaction3 = Reaction(
            reactants=[oxygen_atom, so2],
            products=[so3],
            kinetics=Arrhenius(A=(3.7e11, "cm^3/(mol*s)"), n=0, Ea=(1689, "cal/mol"), T0=(1, "K")),
        )

        H2 = Species().from_smiles("[H][H]")
        PO3 = Species().from_smiles("[O]P(=O)=O")
        HOPO2 = Species().from_smiles("OP(=O)=O")
        H_atom = Species().from_smiles("[H]")

        self.reaction4 = Reaction(
            reactants=[H2, PO3],
            products=[HOPO2, H_atom],
            kinetics=Arrhenius(A=(2.4e7, "cm^3/(mol*s)"), n=1.38, Ea=(15.38, "kcal/mol"), T0=(1, "K")),
        )
        self.reaction4_pairs = [(PO3, HOPO2), (H2, H_atom)]

    def test_is_isomerization(self):
        """
        Test the Reaction.is_isomerization() method.
        """
        isomerization = Reaction(reactants=[Species()], products=[Species()])
        association = Reaction(reactants=[Species(), Species()], products=[Species()])
        dissociation = Reaction(reactants=[Species()], products=[Species(), Species()])
        bimolecular = Reaction(reactants=[Species(), Species()], products=[Species(), Species()])
        assert isomerization.is_isomerization()
        assert not association.is_isomerization()
        assert not dissociation.is_isomerization()
        assert not bimolecular.is_isomerization()

    def test_is_association(self):
        """
        Test the Reaction.is_association() method.
        """
        isomerization = Reaction(reactants=[Species()], products=[Species()])
        association = Reaction(reactants=[Species(), Species()], products=[Species()])
        dissociation = Reaction(reactants=[Species()], products=[Species(), Species()])
        bimolecular = Reaction(reactants=[Species(), Species()], products=[Species(), Species()])
        assert not isomerization.is_association()
        assert association.is_association()
        assert not dissociation.is_association()
        assert not bimolecular.is_association()

    def test_is_dissociation(self):
        """
        Test the Reaction.is_dissociation() method.
        """
        isomerization = Reaction(reactants=[Species()], products=[Species()])
        association = Reaction(reactants=[Species(), Species()], products=[Species()])
        dissociation = Reaction(reactants=[Species()], products=[Species(), Species()])
        bimolecular = Reaction(reactants=[Species(), Species()], products=[Species(), Species()])
        assert not isomerization.is_dissociation()
        assert not association.is_dissociation()
        assert dissociation.is_dissociation()
        assert not bimolecular.is_dissociation()

    def test_has_template(self):
        """
        Test the Reaction.has_template() method.
        """
        reactants = self.reaction.reactants[:]
        products = self.reaction.products[:]
        assert self.reaction.has_template(reactants, products)
        assert self.reaction.has_template(products, reactants)
        assert not self.reaction2.has_template(reactants, products)
        assert not self.reaction2.has_template(products, reactants)

        reactants.reverse()
        products.reverse()
        assert self.reaction.has_template(reactants, products)
        assert self.reaction.has_template(products, reactants)
        assert not self.reaction2.has_template(reactants, products)
        assert not self.reaction2.has_template(products, reactants)

        reactants = self.reaction2.reactants[:]
        products = self.reaction2.products[:]
        assert not self.reaction.has_template(reactants, products)
        assert not self.reaction.has_template(products, reactants)
        assert self.reaction2.has_template(reactants, products)
        assert self.reaction2.has_template(products, reactants)

        reactants.reverse()
        products.reverse()
        assert not self.reaction.has_template(reactants, products)
        assert not self.reaction.has_template(products, reactants)
        assert self.reaction2.has_template(reactants, products)
        assert self.reaction2.has_template(products, reactants)

    def test_enthalpy_of_reaction(self):
        """
        Test the Reaction.get_enthalpy_of_reaction() method.
        """
        Tlist = numpy.arange(200.0, 2001.0, 200.0, numpy.float64)
        Hlist0 = [
            float(v)
            for v in [
                "-146007",
                "-145886",
                "-144195",
                "-141973",
                "-139633",
                "-137341",
                "-135155",
                "-133093",
                "-131150",
                "-129316",
            ]
        ]
        Hlist = self.reaction2.get_enthalpies_of_reaction(Tlist)
        for i in range(len(Tlist)):
            assert round(abs(Hlist[i] / 1000.0 - Hlist0[i] / 1000.0), 2) == 0

    def test_entropy_of_reaction(self):
        """
        Test the Reaction.get_entropy_of_reaction() method.
        """
        Tlist = numpy.arange(200.0, 2001.0, 200.0, numpy.float64)
        Slist0 = [
            float(v)
            for v in [
                "-156.793",
                "-156.872",
                "-153.504",
                "-150.317",
                "-147.707",
                "-145.616",
                "-143.93",
                "-142.552",
                "-141.407",
                "-140.441",
            ]
        ]
        Slist = self.reaction2.get_entropies_of_reaction(Tlist)
        for i in range(len(Tlist)):
            assert round(abs(Slist[i] - Slist0[i]), 2) == 0

    def test_free_energy_of_reaction(self):
        """
        Test the Reaction.get_free_energy_of_reaction() method.
        """
        Tlist = numpy.arange(200.0, 2001.0, 200.0, numpy.float64)
        Glist0 = [
            float(v)
            for v in [
                "-114648",
                "-83137.2",
                "-52092.4",
                "-21719.3",
                "8073.53",
                "37398.1",
                "66346.8",
                "94990.6",
                "123383",
                "151565",
            ]
        ]
        Glist = self.reaction2.get_free_energies_of_reaction(Tlist)
        for i in range(len(Tlist)):
            assert round(abs(Glist[i] / 1000.0 - Glist0[i] / 1000.0), 2) == 0

    def test_equilibrium_constant_ka(self):
        """
        Test the Reaction.get_equilibrium_constant() method.
        """
        Tlist = numpy.arange(200.0, 2001.0, 200.0, numpy.float64)
        Kalist0 = [
            float(v)
            for v in [
                "8.75951e+29",
                "7.1843e+10",
                "34272.7",
                "26.1877",
                "0.378696",
                "0.0235579",
                "0.00334673",
                "0.000792389",
                "0.000262777",
                "0.000110053",
            ]
        ]
        Kalist = self.reaction2.get_equilibrium_constants(Tlist, type="Ka")
        for i in range(len(Tlist)):
            assert round(abs(Kalist[i] / Kalist0[i] - 1.0), 4) == 0

    def test_equilibrium_constant_kc(self):
        """
        Test the Reaction.get_equilibrium_constant() method.
        """
        Tlist = numpy.arange(200.0, 2001.0, 200.0, numpy.float64)
        Kclist0 = [
            float(v)
            for v in [
                "1.45661e+28",
                "2.38935e+09",
                "1709.76",
                "1.74189",
                "0.0314866",
                "0.00235045",
                "0.000389568",
                "0.000105413",
                "3.93273e-05",
                "1.83006e-05",
            ]
        ]
        Kclist = self.reaction2.get_equilibrium_constants(Tlist, type="Kc")
        for i in range(len(Tlist)):
            assert round(abs(Kclist[i] / Kclist0[i] - 1.0), 4) == 0

        rxn2_copy = self.reaction2.copy()
        rxn2_copy.reactants[0].molecule = []
        Kclist_2 = rxn2_copy.get_equilibrium_constants(Tlist, type="Kc")
        for i in range(len(Tlist)):
            assert round(abs(Kclist[i] / Kclist_2[i] - 1.0), 4) == 0

    def test_equilibrium_constant_kp(self):
        """
        Test the Reaction.get_equilibrium_constant() method.
        """
        Tlist = numpy.arange(200.0, 2001.0, 200.0, numpy.float64)
        Kplist0 = [
            float(v)
            for v in [
                "8.75951e+24",
                "718430",
                "0.342727",
                "0.000261877",
                "3.78696e-06",
                "2.35579e-07",
                "3.34673e-08",
                "7.92389e-09",
                "2.62777e-09",
                "1.10053e-09",
            ]
        ]
        Kplist = self.reaction2.get_equilibrium_constants(Tlist, type="Kp")
        for i in range(len(Tlist)):
            assert round(abs(Kplist[i] / Kplist0[i] - 1.0), 4) == 0

    def test_stoichiometric_coefficient(self):
        """
        Test the Reaction.get_stoichiometric_coefficient() method.
        """
        for reactant in self.reaction.reactants:
            assert self.reaction.get_stoichiometric_coefficient(reactant) == -1
        for product in self.reaction.products:
            assert self.reaction.get_stoichiometric_coefficient(product) == 1
        for reactant in self.reaction2.reactants:
            assert self.reaction.get_stoichiometric_coefficient(reactant) == 0
        for product in self.reaction2.products:
            assert self.reaction.get_stoichiometric_coefficient(product) == 0

    def test_rate_coefficient(self):
        """
        Test the Reaction.get_rate_coefficient() method.
        """
        Tlist = numpy.arange(200.0, 2001.0, 200.0, numpy.float64)
        P = 1e5
        for T in Tlist:
            assert round(abs(self.reaction.get_rate_coefficient(T, P) / self.reaction.kinetics.get_rate_coefficient(T) - 1.0), 6) == 0

    def test_generate_reverse_rate_coefficient(self):
        """
        Test the Reaction.generate_reverse_rate_coefficient() method.
        """
        Tlist = numpy.arange(200.0, 2001.0, 200.0, numpy.float64)
        P = 1e5
        reverse_kinetics = self.reaction2.generate_reverse_rate_coefficient()
        for T in Tlist:
            kr0 = self.reaction2.get_rate_coefficient(T, P) / self.reaction2.get_equilibrium_constant(T)
            kr = reverse_kinetics.get_rate_coefficient(T)
            assert round(abs(kr0 / kr - 1.0), 0) == 0

    def test_fix_barrier_height(self):
        """
        Test that fix_barrier_height:
            1) raises Ea to match endothermicity of reaction
            2) forces Ea to be positive if force_positive=True
            3) Evans-Polanyi kinetics are handled so that negative Ea if Ea<E0 are set to min(0,E0)
        """

        # setup
        rxn = self.reaction2.copy()
        rev_rxn = rxn.copy()
        rev_rxn.reactants = rxn.products
        rev_rxn.products = rxn.reactants

        # test that endothermicity is matched
        rxn.fix_barrier_height()
        Ea = rxn.kinetics.Ea.value_si
        assert Ea == 0.0

        rev_rxn.fix_barrier_height()
        Ea = rev_rxn.kinetics.Ea.value_si
        H0 = sum([spec.get_thermo_data().E0.value_si for spec in rxn.products]) - sum([spec.get_thermo_data().E0.value_si for spec in rxn.reactants])
        assert round(abs(Ea - -H0), 3) == 0

        # test that Ea is forced to be positive if force_positive is set to True
        Ea = Quantity((-10000.0, "J/mol"))
        rxn.kinetics.Ea = Ea
        rxn.fix_barrier_height()
        assert rxn.kinetics.Ea.value_si == Ea.value_si

        rxn.fix_barrier_height(force_positive=True)
        assert rxn.kinetics.Ea.value_si == 0.0

        # Test for ArrheniusEP handling
        # if calculated Ea < 0 and Ea < E0, Ea is set to min(0,E0)
        H298 = rxn.get_enthalpy_of_reaction(298)
        E0s = [-1000000.0, -10.0, 0.0, 10.0, 1000000.0]

        for i, E0 in enumerate(E0s):
            kinetics = ArrheniusEP(
                A=(1.0, rxn.kinetics.A.units),
                n=(0, rxn.kinetics.n.units),
                alpha=1.0,
                E0=(E0, "J/mol"),
            )
            rxn.kinetics = kinetics
            rxn.fix_barrier_height()
            Ea = rxn.kinetics.Ea.value_si
            if i < 2:
                assert Ea == E0
            elif i < 4:
                assert Ea == 0.0
            else:
                assert Ea == E0 + H298

    def test_generate_reverse_rate_coefficient_arrhenius(self):
        """
        Test the Reaction.generate_reverse_rate_coefficient() method works for the Arrhenius format.
        """
        original_kinetics = Arrhenius(
            A=(2.65e12, "cm^3/(mol*s)"),
            n=0.0,
            Ea=(0.0, "kJ/mol"),
            T0=(1, "K"),
            Tmin=(300, "K"),
            Tmax=(2000, "K"),
        )
        self.reaction2.kinetics = original_kinetics

        reverse_kinetics = self.reaction2.generate_reverse_rate_coefficient()

        self.reaction2.kinetics = reverse_kinetics
        # reverse reactants, products to ensure Keq is correctly computed
        self.reaction2.reactants, self.reaction2.products = (
            self.reaction2.products,
            self.reaction2.reactants,
        )
        reverse_reverse_kinetics = self.reaction2.generate_reverse_rate_coefficient()

        # check that reverting the reverse yields the original
        Tlist = numpy.arange(
            original_kinetics.Tmin.value_si,
            original_kinetics.Tmax.value_si,
            200.0,
            numpy.float64,
        )
        P = 1e5
        for T in Tlist:
            korig = original_kinetics.get_rate_coefficient(T, P)
            krevrev = reverse_reverse_kinetics.get_rate_coefficient(T, P)
            assert round(abs(korig / krevrev - 1.0), 0) == 0

    def test_reverse_surface_arrhenius_rate(self):
        """
        Test the Reaction.reverse_surface_arrhenius_rate() method works for SurfaceArrhenius format.
        """
        original_kinetics = SurfaceArrhenius(
            A=(1.195e12, "m^2/(mol*s)"),
            n=0.0,
            Ea=(14.989, "kcal/mol"),
            T0=(1, "K"),
            Tmin=(300, "K"),
            Tmax=(2000, "K"),
        )
        self.reaction2.kinetics = original_kinetics

        reverse_kinetics = self.reaction2.generate_reverse_rate_coefficient()

        self.reaction2.kinetics = reverse_kinetics
        # reverse reactants, products to ensure Keq is correctly computed
        self.reaction2.reactants, self.reaction2.products = (
            self.reaction2.products,
            self.reaction2.reactants,
        )
        reverse_reverse_kinetics = self.reaction2.generate_reverse_rate_coefficient()

        # check that reverting the reverse yields the original
        Tlist = numpy.arange(
            original_kinetics.Tmin.value_si,
            original_kinetics.Tmax.value_si,
            200.0,
            numpy.float64,
        )
        P = 1e5
        for T in Tlist:
            korig = original_kinetics.get_rate_coefficient(T, P)
            krevrev = reverse_reverse_kinetics.get_rate_coefficient(T, P)
            assert round(abs(korig / krevrev - 1.0), 0) == 0

    @pytest.mark.skip(reason="WIP")
    def test_generate_reverse_rate_coefficient_arrhenius_ep(self):
        """
        Test the Reaction.generate_reverse_rate_coefficient() method works for the ArrheniusEP format.
        """

        original_kinetics = ArrheniusEP(
            A=(2.65e12, "cm^3/(mol*s)"),
            n=0.0,
            alpha=0.5,
            E0=(41.84, "kJ/mol"),
            Tmin=(300, "K"),
            Tmax=(2000, "K"),
        )
        self.reaction2.kinetics = original_kinetics

        reverse_kinetics = self.reaction2.generate_reverse_rate_coefficient()

        self.reaction2.kinetics = reverse_kinetics
        # reverse reactants, products to ensure Keq is correctly computed
        self.reaction2.reactants, self.reaction2.products = (
            self.reaction2.products,
            self.reaction2.reactants,
        )
        reverse_reverse_kinetics = self.reaction2.generate_reverse_rate_coefficient()

        # check that reverting the reverse yields the original
        Tlist = numpy.arange(original_kinetics.Tmin, original_kinetics.Tmax, 200.0, numpy.float64)
        P = 1e5
        for T in Tlist:
            korig = original_kinetics.get_rate_coefficient(T, P)
            krevrev = reverse_reverse_kinetics.get_rate_coefficient(T, P)
            assert round(abs(korig / krevrev - 1.0), 0) == 0

    def test_generate_reverse_rate_coefficient_pdep_arrhenius(self):
        """
        Test the Reaction.generate_reverse_rate_coefficient() method works for the PDepArrhenius format.
        """

        arrhenius0 = Arrhenius(
            A=(1.0e6, "s^-1"),
            n=1.0,
            Ea=(10.0, "kJ/mol"),
            T0=(300.0, "K"),
            Tmin=(300.0, "K"),
            Tmax=(2000.0, "K"),
            comment="""This data is completely made up""",
        )

        arrhenius1 = Arrhenius(
            A=(1.0e12, "s^-1"),
            n=1.0,
            Ea=(20.0, "kJ/mol"),
            T0=(300.0, "K"),
            Tmin=(300.0, "K"),
            Tmax=(2000.0, "K"),
            comment="""This data is completely made up""",
        )

        pressures = numpy.array([0.1, 10.0])
        arrhenius = [arrhenius0, arrhenius1]
        Tmin = 300.0
        Tmax = 2000.0
        Pmin = 0.1
        Pmax = 10.0
        comment = """This data is completely made up"""

        original_kinetics = PDepArrhenius(
            pressures=(pressures, "bar"),
            arrhenius=arrhenius,
            Tmin=(Tmin, "K"),
            Tmax=(Tmax, "K"),
            Pmin=(Pmin, "bar"),
            Pmax=(Pmax, "bar"),
            comment=comment,
        )

        self.reaction2.kinetics = original_kinetics

        reverse_kinetics = self.reaction2.generate_reverse_rate_coefficient()

        self.reaction2.kinetics = reverse_kinetics
        # reverse reactants, products to ensure Keq is correctly computed
        self.reaction2.reactants, self.reaction2.products = (
            self.reaction2.products,
            self.reaction2.reactants,
        )
        reverse_reverse_kinetics = self.reaction2.generate_reverse_rate_coefficient()

        # check that reverting the reverse yields the original
        Tlist = numpy.arange(Tmin, Tmax, 200.0, numpy.float64)
        P = 1e5
        for T in Tlist:
            korig = original_kinetics.get_rate_coefficient(T, P)
            krevrev = reverse_reverse_kinetics.get_rate_coefficient(T, P)
            assert round(abs(korig / krevrev - 1.0), 0) == 0

    def test_generate_reverse_rate_coefficient_pdep_multi_arrhenius(self):
        """
        Test the Reaction.generate_reverse_rate_coefficient() method works for the PDepArrhenius format with MultiArrhenius rates.
        """

        arrhenius0 = MultiArrhenius(
            arrhenius=[
                Arrhenius(
                    A=(1.0e6, "s^-1"),
                    n=1.0,
                    Ea=(10.0, "kJ/mol"),
                    T0=(300.0, "K"),
                    Tmin=(300.0, "K"),
                    Tmax=(2000.0, "K"),
                ),
                Arrhenius(
                    A=(1.0e6, "s^-1"),
                    n=1.0,
                    Ea=(10.0, "kJ/mol"),
                    T0=(300.0, "K"),
                    Tmin=(300.0, "K"),
                    Tmax=(2000.0, "K"),
                ),
            ],
            comment="""This data is completely made up""",
        )

        arrhenius1 = MultiArrhenius(
            arrhenius=[
                Arrhenius(
                    A=(1.0e12, "s^-1"),
                    n=1.0,
                    Ea=(10.0, "kJ/mol"),
                    T0=(300.0, "K"),
                    Tmin=(300.0, "K"),
                    Tmax=(2000.0, "K"),
                ),
                Arrhenius(
                    A=(1.0e12, "s^-1"),
                    n=1.0,
                    Ea=(10.0, "kJ/mol"),
                    T0=(300.0, "K"),
                    Tmin=(300.0, "K"),
                    Tmax=(2000.0, "K"),
                ),
            ],
            comment="""This data is completely made up""",
        )

        pressures = numpy.array([0.1, 10.0])
        arrhenius = [arrhenius0, arrhenius1]
        Tmin = 300.0
        Tmax = 2000.0
        Pmin = 0.1
        Pmax = 10.0
        comment = """This data is completely made up"""

        original_kinetics = PDepArrhenius(
            pressures=(pressures, "bar"),
            arrhenius=arrhenius,
            Tmin=(Tmin, "K"),
            Tmax=(Tmax, "K"),
            Pmin=(Pmin, "bar"),
            Pmax=(Pmax, "bar"),
            comment=comment,
        )

        self.reaction2.kinetics = original_kinetics

        reverse_kinetics = self.reaction2.generate_reverse_rate_coefficient()

        self.reaction2.kinetics = reverse_kinetics
        # reverse reactants, products to ensure Keq is correctly computed
        self.reaction2.reactants, self.reaction2.products = (
            self.reaction2.products,
            self.reaction2.reactants,
        )
        reverse_reverse_kinetics = self.reaction2.generate_reverse_rate_coefficient()

        # check that reverting the reverse yields the original
        Tlist = numpy.arange(Tmin, Tmax, 200.0, numpy.float64)
        P = 1e5
        for T in Tlist:
            korig = original_kinetics.get_rate_coefficient(T, P)
            krevrev = reverse_reverse_kinetics.get_rate_coefficient(T, P)
            assert round(abs(korig / krevrev - 1.0), 0) == 0

    def test_generate_reverse_rate_coefficient_multi_arrhenius(self):
        """
        Test the Reaction.generate_reverse_rate_coefficient() method works for the MultiArrhenius format.
        """

        pressures = numpy.array([0.1, 10.0])
        Tmin = 300.0
        Tmax = 2000.0
        Pmin = 0.1
        Pmax = 10.0
        comment = """This data is completely made up"""

        arrhenius = [
            Arrhenius(
                A=(9.3e-14, "cm^3/(molecule*s)"),
                n=0.0,
                Ea=(4740 * constants.R * 0.001, "kJ/mol"),
                T0=(1, "K"),
                Tmin=(Tmin, "K"),
                Tmax=(Tmax, "K"),
                comment=comment,
            ),
            Arrhenius(
                A=(1.4e-9, "cm^3/(molecule*s)"),
                n=0.0,
                Ea=(11200 * constants.R * 0.001, "kJ/mol"),
                T0=(1, "K"),
                Tmin=(Tmin, "K"),
                Tmax=(Tmax, "K"),
                comment=comment,
            ),
        ]

        original_kinetics = MultiArrhenius(
            arrhenius=arrhenius,
            Tmin=(Tmin, "K"),
            Tmax=(Tmax, "K"),
            comment=comment,
        )

        self.reaction2.kinetics = original_kinetics

        reverse_kinetics = self.reaction2.generate_reverse_rate_coefficient()

        self.reaction2.kinetics = reverse_kinetics
        # reverse reactants, products to ensure Keq is correctly computed
        self.reaction2.reactants, self.reaction2.products = (
            self.reaction2.products,
            self.reaction2.reactants,
        )
        reverse_reverse_kinetics = self.reaction2.generate_reverse_rate_coefficient()

        # check that reverting the reverse yields the original
        Tlist = numpy.arange(Tmin, Tmax, 200.0, numpy.float64)
        P = 1e5
        for T in Tlist:
            korig = original_kinetics.get_rate_coefficient(T, P)
            krevrev = reverse_reverse_kinetics.get_rate_coefficient(T, P)
            assert round(abs(korig / krevrev - 1.0), 0) == 0

    def test_generate_reverse_rate_coefficient_multi_pdep_arrhenius(self):
        """
        Test the Reaction.generate_reverse_rate_coefficient() method works for the MultiPDepArrhenius format.
        """

        Tmin = 350.0
        Tmax = 1500.0
        Pmin = 1e-1
        Pmax = 1e1
        pressures = numpy.array([1e-1, 1e1])
        comment = "CH3 + C2H6 <=> CH4 + C2H5 (Baulch 2005)"
        arrhenius = [
            PDepArrhenius(
                pressures=(pressures, "bar"),
                arrhenius=[
                    Arrhenius(
                        A=(9.3e-16, "cm^3/(molecule*s)"),
                        n=0.0,
                        Ea=(4740 * constants.R * 0.001, "kJ/mol"),
                        T0=(1, "K"),
                        Tmin=(Tmin, "K"),
                        Tmax=(Tmax, "K"),
                        comment=comment,
                    ),
                    Arrhenius(
                        A=(9.3e-14, "cm^3/(molecule*s)"),
                        n=0.0,
                        Ea=(4740 * constants.R * 0.001, "kJ/mol"),
                        T0=(1, "K"),
                        Tmin=(Tmin, "K"),
                        Tmax=(Tmax, "K"),
                        comment=comment,
                    ),
                ],
                Tmin=(Tmin, "K"),
                Tmax=(Tmax, "K"),
                Pmin=(Pmin, "bar"),
                Pmax=(Pmax, "bar"),
                comment=comment,
            ),
            PDepArrhenius(
                pressures=(pressures, "bar"),
                arrhenius=[
                    Arrhenius(
                        A=(1.4e-11, "cm^3/(molecule*s)"),
                        n=0.0,
                        Ea=(11200 * constants.R * 0.001, "kJ/mol"),
                        T0=(1, "K"),
                        Tmin=(Tmin, "K"),
                        Tmax=(Tmax, "K"),
                        comment=comment,
                    ),
                    Arrhenius(
                        A=(1.4e-9, "cm^3/(molecule*s)"),
                        n=0.0,
                        Ea=(11200 * constants.R * 0.001, "kJ/mol"),
                        T0=(1, "K"),
                        Tmin=(Tmin, "K"),
                        Tmax=(Tmax, "K"),
                        comment=comment,
                    ),
                ],
                Tmin=(Tmin, "K"),
                Tmax=(Tmax, "K"),
                Pmin=(Pmin, "bar"),
                Pmax=(Pmax, "bar"),
                comment=comment,
            ),
        ]

        original_kinetics = MultiPDepArrhenius(
            arrhenius=arrhenius,
            Tmin=(Tmin, "K"),
            Tmax=(Tmax, "K"),
            Pmin=(Pmin, "bar"),
            Pmax=(Pmax, "bar"),
            comment=comment,
        )

        self.reaction2.kinetics = original_kinetics

        reverse_kinetics = self.reaction2.generate_reverse_rate_coefficient()

        self.reaction2.kinetics = reverse_kinetics
        # reverse reactants, products to ensure Keq is correctly computed
        self.reaction2.reactants, self.reaction2.products = (
            self.reaction2.products,
            self.reaction2.reactants,
        )
        reverse_reverse_kinetics = self.reaction2.generate_reverse_rate_coefficient()

        # check that reverting the reverse yields the original
        Tlist = numpy.arange(Tmin, Tmax, 200.0, numpy.float64)
        P = 1e5
        for T in Tlist:
            korig = original_kinetics.get_rate_coefficient(T, P)
            krevrev = reverse_reverse_kinetics.get_rate_coefficient(T, P)
            assert round(abs(korig / krevrev - 1.0), 0) == 0

    def test_generate_reverse_rate_coefficient_third_body(self):
        """
        Test the Reaction.generate_reverse_rate_coefficient() method works for the ThirdBody format.
        """

        arrhenius_low = Arrhenius(
            A=(2.62e33, "cm^6/(mol^2*s)"),
            n=-4.76,
            Ea=(10.21, "kJ/mol"),
            T0=(1, "K"),
        )
        efficiencies = {
            "C": 3,
            "C(=O)=O": 2,
            "CC": 3,
            "O": 6,
            "[Ar]": 0.7,
            "[C]=O": 1.5,
            "[H][H]": 2,
        }
        Tmin = 300.0
        Tmax = 2000.0
        Pmin = 0.01
        Pmax = 100.0
        comment = """H + CH3 -> CH4"""
        third_body = ThirdBody(
            arrheniusLow=arrhenius_low,
            Tmin=(Tmin, "K"),
            Tmax=(Tmax, "K"),
            Pmin=(Pmin, "bar"),
            Pmax=(Pmax, "bar"),
            efficiencies=efficiencies,
            comment=comment,
        )

        original_kinetics = third_body

        self.reaction2.kinetics = original_kinetics

        reverse_kinetics = self.reaction2.generate_reverse_rate_coefficient()

        self.reaction2.kinetics = reverse_kinetics
        # reverse reactants, products to ensure Keq is correctly computed
        self.reaction2.reactants, self.reaction2.products = (
            self.reaction2.products,
            self.reaction2.reactants,
        )
        reverse_reverse_kinetics = self.reaction2.generate_reverse_rate_coefficient()

        # check that reverting the reverse yields the original
        Tlist = numpy.arange(Tmin, Tmax, 200.0, numpy.float64)
        P = 1e5
        for T in Tlist:
            korig = original_kinetics.get_rate_coefficient(T, P)
            krevrev = reverse_reverse_kinetics.get_rate_coefficient(T, P)
            assert round(abs(korig / krevrev - 1.0), 0) == 0

    def test_generate_reverse_rate_coefficient_lindemann(self):
        """
        Test the Reaction.generate_reverse_rate_coefficient() method works for the Lindemann format.
        """

        arrhenius_high = Arrhenius(
            A=(1.39e16, "cm^3/(mol*s)"),
            n=-0.534,
            Ea=(2.243, "kJ/mol"),
            T0=(1, "K"),
        )
        arrhenius_low = Arrhenius(
            A=(2.62e33, "cm^6/(mol^2*s)"),
            n=-4.76,
            Ea=(10.21, "kJ/mol"),
            T0=(1, "K"),
        )
        efficiencies = {
            "C": 3,
            "C(=O)=O": 2,
            "CC": 3,
            "O": 6,
            "[Ar]": 0.7,
            "[C]=O": 1.5,
            "[H][H]": 2,
        }
        Tmin = 300.0
        Tmax = 2000.0
        Pmin = 0.01
        Pmax = 100.0
        comment = """H + CH3 -> CH4"""
        lindemann = Lindemann(
            arrheniusHigh=arrhenius_high,
            arrheniusLow=arrhenius_low,
            Tmin=(Tmin, "K"),
            Tmax=(Tmax, "K"),
            Pmin=(Pmin, "bar"),
            Pmax=(Pmax, "bar"),
            efficiencies=efficiencies,
            comment=comment,
        )

        original_kinetics = lindemann

        self.reaction2.kinetics = original_kinetics

        reverse_kinetics = self.reaction2.generate_reverse_rate_coefficient()

        self.reaction2.kinetics = reverse_kinetics
        # reverse reactants, products to ensure Keq is correctly computed
        self.reaction2.reactants, self.reaction2.products = (
            self.reaction2.products,
            self.reaction2.reactants,
        )
        reverse_reverse_kinetics = self.reaction2.generate_reverse_rate_coefficient()

        # check that reverting the reverse yields the original
        Tlist = numpy.arange(Tmin, Tmax, 200.0, numpy.float64)
        P = 1e5
        for T in Tlist:
            korig = original_kinetics.get_rate_coefficient(T, P)
            krevrev = reverse_reverse_kinetics.get_rate_coefficient(T, P)
            assert round(abs(korig / krevrev - 1.0), 0) == 0

    def test_generate_reverse_rate_coefficient_troe(self):
        """
        Test the Reaction.generate_reverse_rate_coefficient() method works for the Troe format.
        """

        arrhenius_high = Arrhenius(
            A=(1.39e16, "cm^3/(mol*s)"),
            n=-0.534,
            Ea=(2.243, "kJ/mol"),
            T0=(1, "K"),
        )
        arrhenius_low = Arrhenius(
            A=(2.62e33, "cm^6/(mol^2*s)"),
            n=-4.76,
            Ea=(10.21, "kJ/mol"),
            T0=(1, "K"),
        )
        alpha = 0.783
        T3 = 74
        T1 = 2941
        T2 = 6964
        efficiencies = {
            "C": 3,
            "C(=O)=O": 2,
            "CC": 3,
            "O": 6,
            "[Ar]": 0.7,
            "[C]=O": 1.5,
            "[H][H]": 2,
        }
        Tmin = 300.0
        Tmax = 2000.0
        Pmin = 0.01
        Pmax = 100.0
        comment = """H + CH3 -> CH4"""
        troe = Troe(
            arrheniusHigh=arrhenius_high,
            arrheniusLow=arrhenius_low,
            alpha=alpha,
            T3=(T3, "K"),
            T1=(T1, "K"),
            T2=(T2, "K"),
            Tmin=(Tmin, "K"),
            Tmax=(Tmax, "K"),
            Pmin=(Pmin, "bar"),
            Pmax=(Pmax, "bar"),
            efficiencies=efficiencies,
            comment=comment,
        )

        original_kinetics = troe

        self.reaction2.kinetics = original_kinetics

        reverse_kinetics = self.reaction2.generate_reverse_rate_coefficient()

        self.reaction2.kinetics = reverse_kinetics
        # reverse reactants, products to ensure Keq is correctly computed
        self.reaction2.reactants, self.reaction2.products = (
            self.reaction2.products,
            self.reaction2.reactants,
        )
        reverse_reverse_kinetics = self.reaction2.generate_reverse_rate_coefficient()

        # check that reverting the reverse yields the original
        Tlist = numpy.arange(Tmin, Tmax, 200.0, numpy.float64)
        P = 1e5
        for T in Tlist:
            korig = original_kinetics.get_rate_coefficient(T, P)
            krevrev = reverse_reverse_kinetics.get_rate_coefficient(T, P)
            assert round(abs(korig / krevrev - 1.0), 0) == 0

    def test_tst_calculation(self):
        """
        A test of the transition state theory k(T) calculation function,
        using the reaction H + C2H4 -> C2H5.
        """
        Tlist = 1000.0 / numpy.arange(0.4, 3.35, 0.01)
        klist = numpy.array([self.reaction.calculate_tst_rate_coefficient(T) for T in Tlist])
        arrhenius = Arrhenius().fit_to_data(Tlist, klist, kunits="m^3/(mol*s)")
        klist2 = numpy.array([arrhenius.get_rate_coefficient(T) for T in Tlist])

        # Check that the correct Arrhenius parameters are returned
        assert abs(arrhenius.A.value_si - 2265.2488) < 1e-2
        assert abs(arrhenius.n.value_si - 1.45419) < 1e-4
        assert abs(arrhenius.Ea.value_si - 6645.24) < 1e-2
        # Check that the fit is satisfactory (defined here as always within 5%)
        for i in range(len(Tlist)):
            assert abs(klist[i] - klist2[i]) < 5e-2 * klist[i]

    def test_pickle(self):
        """
        Test that a Reaction object can be successfully pickled and
        unpickled with no loss of information.
        """
        import pickle

        reaction = pickle.loads(pickle.dumps(self.reaction, -1))

        assert len(self.reaction.reactants) == len(reaction.reactants)
        assert len(self.reaction.products) == len(reaction.products)
        for reactant0, reactant in zip(self.reaction.reactants, reaction.reactants):
            assert round(abs(reactant0.conformer.E0.value_si / 1e6 - reactant.conformer.E0.value_si / 1e6), 2) == 0
            assert reactant0.conformer.E0.units == reactant.conformer.E0.units
        for product0, product in zip(self.reaction.products, reaction.products):
            assert round(abs(product0.conformer.E0.value_si / 1e6 - product.conformer.E0.value_si / 1e6), 2) == 0
            assert product0.conformer.E0.units == product.conformer.E0.units
        assert round(abs(self.reaction.transition_state.conformer.E0.value_si / 1e6 - reaction.transition_state.conformer.E0.value_si / 1e6), 2) == 0
        assert self.reaction.transition_state.conformer.E0.units == reaction.transition_state.conformer.E0.units
        assert round(abs(self.reaction.transition_state.frequency.value_si - reaction.transition_state.frequency.value_si), 2) == 0
        assert self.reaction.transition_state.frequency.units == reaction.transition_state.frequency.units

        assert abs(self.reaction.kinetics.A.value_si - reaction.kinetics.A.value_si) < 1e-6
        assert abs(self.reaction.kinetics.n.value_si - reaction.kinetics.n.value_si) < 1e-6
        assert abs(self.reaction.kinetics.T0.value_si - reaction.kinetics.T0.value_si) < 1e-6
        assert abs(self.reaction.kinetics.Ea.value_si - reaction.kinetics.Ea.value_si) < 1e-6
        assert self.reaction.kinetics.comment == reaction.kinetics.comment

        assert self.reaction.duplicate == reaction.duplicate
        assert self.reaction.degeneracy == reaction.degeneracy

    def test_output(self):
        """
        Test that a Reaction object can be successfully reconstructed
        from its repr() output with no loss of information.
        """
        namespace = {}
        exec("reaction = {0!r}".format(self.reaction), globals(), namespace)
        assert "reaction" in namespace
        reaction = namespace["reaction"]

        assert len(self.reaction.reactants) == len(reaction.reactants)
        assert len(self.reaction.products) == len(reaction.products)
        for reactant0, reactant in zip(self.reaction.reactants, reaction.reactants):
            assert round(abs(reactant0.conformer.E0.value_si / 1e6 - reactant.conformer.E0.value_si / 1e6), 2) == 0
            assert reactant0.conformer.E0.units == reactant.conformer.E0.units
        for product0, product in zip(self.reaction.products, reaction.products):
            assert round(abs(product0.conformer.E0.value_si / 1e6 - product.conformer.E0.value_si / 1e6), 2) == 0
            assert product0.conformer.E0.units == product.conformer.E0.units
        assert round(abs(self.reaction.transition_state.conformer.E0.value_si / 1e6 - reaction.transition_state.conformer.E0.value_si / 1e6), 2) == 0
        assert self.reaction.transition_state.conformer.E0.units == reaction.transition_state.conformer.E0.units
        assert round(abs(self.reaction.transition_state.frequency.value_si - reaction.transition_state.frequency.value_si), 2) == 0
        assert self.reaction.transition_state.frequency.units == reaction.transition_state.frequency.units

        assert abs(self.reaction.kinetics.A.value_si - reaction.kinetics.A.value_si) < 1e-6
        assert abs(self.reaction.kinetics.n.value_si - reaction.kinetics.n.value_si) < 1e-6
        assert abs(self.reaction.kinetics.T0.value_si - reaction.kinetics.T0.value_si) < 1e-6
        assert abs(self.reaction.kinetics.Ea.value_si - reaction.kinetics.Ea.value_si) < 1e-6
        assert self.reaction.kinetics.comment == reaction.kinetics.comment

        assert self.reaction.duplicate == reaction.duplicate
        assert self.reaction.degeneracy == reaction.degeneracy

    def test_degeneracy_updates_rate(self):
        """
        This method tests that a change in degeneracy will result in a modified rate constant
        """

        prefactor = self.reaction.kinetics.A.value_si
        degeneracyFactor = 2
        self.reaction.degeneracy *= degeneracyFactor
        assert round(abs(self.reaction.kinetics.A.value_si - degeneracyFactor * prefactor), 7) == 0

    def test_degeneracy_updates_kinetics_comment(self):
        """
        This method tests that a change in degeneracy will result in a modified rate constant
        """

        newDegeneracy = 8
        self.reaction.degeneracy = newDegeneracy
        assert "Multiplied by reaction path degeneracy 8.0" in self.reaction.kinetics.comment

    def test_sulfur_reaction_pairs(self):
        """
        This method tests that reaction pairs are being generated for sulfur species
        """

        self.reaction3.generate_pairs()
        assert len(self.reaction3.pairs[0]) == 2
        assert len(self.reaction3.pairs[1]) == 2

    def test_phosphorus_reaction_pairs(self):
        """
        This method tests that reaction pairs are being generated for phosphorus species
        """

        self.reaction4.generate_pairs()
        assert len(self.reaction4.pairs[0]) == 2
        assert len(self.reaction4.pairs[1]) == 2
        assert self.reaction4.pairs == self.reaction4_pairs


class TestReactionToCantera:
    """
    Contains unit tests of the Reaction class associated with forming Cantera objects.
    """

    def setup_class(self):
        """
        A method that is called prior to each unit test in this class.
        """
        # define some species:
        ch3 = Species(
            index=13,
            label="CH3",
            thermo=NASA(
                polynomials=[
                    NASAPolynomial(
                        coeffs=[
                            3.91547,
                            0.00184154,
                            3.48744e-06,
                            -3.3275e-09,
                            8.49964e-13,
                            16285.6,
                            0.351739,
                        ],
                        Tmin=(100, "K"),
                        Tmax=(1337.62, "K"),
                    ),
                    NASAPolynomial(
                        coeffs=[
                            3.54145,
                            0.00476788,
                            -1.82149e-06,
                            3.28878e-10,
                            -2.22547e-14,
                            16224,
                            1.6604,
                        ],
                        Tmin=(1337.62, "K"),
                        Tmax=(5000, "K"),
                    ),
                ],
                Tmin=(100, "K"),
                Tmax=(5000, "K"),
                comment="""
Thermo library: primaryThermoLibrary + radical(CH3)
""",
            ),
            molecule=[Molecule(smiles="[CH3]")],
        )

        ethane = Species(
            label="ethane",
            thermo=NASA(
                polynomials=[
                    NASAPolynomial(
                        coeffs=[
                            3.78033,
                            -0.00324263,
                            5.52381e-05,
                            -6.38581e-08,
                            2.28637e-11,
                            -11620.3,
                            5.21034,
                        ],
                        Tmin=(100, "K"),
                        Tmax=(954.51, "K"),
                    ),
                    NASAPolynomial(
                        coeffs=[
                            4.58983,
                            0.0141508,
                            -4.75962e-06,
                            8.60294e-10,
                            -6.21717e-14,
                            -12721.8,
                            -3.61739,
                        ],
                        Tmin=(954.51, "K"),
                        Tmax=(5000, "K"),
                    ),
                ],
                Tmin=(100, "K"),
                Tmax=(5000, "K"),
                comment="""
Thermo group additivity estimation: group(Cs-CsHHH) + gauche(Cs(CsRRR)) + other(R) + group(Cs-CsHHH) + gauche(Cs(CsRRR)) + other(R)
""",
            ),
            molecule=[Molecule(smiles="CC")],
        )

        co2 = Species(
            index=16,
            label="CO2",
            thermo=NASA(
                polynomials=[
                    NASAPolynomial(
                        coeffs=[
                            3.27861,
                            0.00274152,
                            7.16074e-06,
                            -1.08027e-08,
                            4.14282e-12,
                            -48470.3,
                            5.97937,
                        ],
                        Tmin=(100, "K"),
                        Tmax=(988.89, "K"),
                    ),
                    NASAPolynomial(
                        coeffs=[
                            4.5461,
                            0.00291913,
                            -1.15484e-06,
                            2.27654e-10,
                            -1.7091e-14,
                            -48980.4,
                            -1.43275,
                        ],
                        Tmin=(988.89, "K"),
                        Tmax=(5000, "K"),
                    ),
                ],
                Tmin=(100, "K"),
                Tmax=(5000, "K"),
                comment="""
Thermo group additivity estimation: group(Cdd-OdOd) + other(R) + group(O2d-Cd) + other(R) + group(O2d-Cd) + other(R)
""",
            ),
            molecule=[Molecule(smiles="O=C=O")],
        )

        ch4 = Species(
            index=15,
            label="CH4",
            thermo=NASA(
                polynomials=[
                    NASAPolynomial(
                        coeffs=[
                            4.20541,
                            -0.00535556,
                            2.51123e-05,
                            -2.13762e-08,
                            5.97522e-12,
                            -10161.9,
                            -0.921275,
                        ],
                        Tmin=(100, "K"),
                        Tmax=(1084.12, "K"),
                    ),
                    NASAPolynomial(
                        coeffs=[
                            0.908272,
                            0.0114541,
                            -4.57173e-06,
                            8.2919e-10,
                            -5.66314e-14,
                            -9719.98,
                            13.9931,
                        ],
                        Tmin=(1084.12, "K"),
                        Tmax=(5000, "K"),
                    ),
                ],
                Tmin=(100, "K"),
                Tmax=(5000, "K"),
                comment="""
Thermo library: primaryThermoLibrary
""",
            ),
            molecule=[Molecule(smiles="C")],
        )

        h2o = Species(
            index=27,
            label="H2O",
            thermo=NASA(
                polynomials=[
                    NASAPolynomial(
                        coeffs=[
                            4.05764,
                            -0.000787933,
                            2.90876e-06,
                            -1.47518e-09,
                            2.12838e-13,
                            -30281.6,
                            -0.311363,
                        ],
                        Tmin=(100, "K"),
                        Tmax=(1130.24, "K"),
                    ),
                    NASAPolynomial(
                        coeffs=[
                            2.84325,
                            0.00275108,
                            -7.8103e-07,
                            1.07243e-10,
                            -5.79389e-15,
                            -29958.6,
                            5.91041,
                        ],
                        Tmin=(1130.24, "K"),
                        Tmax=(5000, "K"),
                    ),
                ],
                Tmin=(100, "K"),
                Tmax=(5000, "K"),
                comment="""
Thermo library: primaryThermoLibrary
""",
            ),
            molecule=[Molecule(smiles="O")],
        )

        ar = Species(
            label="Ar",
            thermo=NASA(
                polynomials=[
                    NASAPolynomial(
                        coeffs=[2.5, 0, 0, 0, 0, -745.375, 4.37967],
                        Tmin=(200, "K"),
                        Tmax=(1000, "K"),
                    ),
                    NASAPolynomial(
                        coeffs=[2.5, 0, 0, 0, 0, -745.375, 4.37967],
                        Tmin=(1000, "K"),
                        Tmax=(6000, "K"),
                    ),
                ],
                Tmin=(200, "K"),
                Tmax=(6000, "K"),
                comment="""
Thermo library: primaryThermoLibrary
""",
            ),
            molecule=[Molecule(smiles="[Ar]")],
        )

        h2 = Species(
            index=2,
            label="H2",
            thermo=NASA(
                polynomials=[
                    NASAPolynomial(
                        coeffs=[
                            3.43536,
                            0.00021271,
                            -2.78625e-07,
                            3.40267e-10,
                            -7.76031e-14,
                            -1031.36,
                            -3.90842,
                        ],
                        Tmin=(100, "K"),
                        Tmax=(1959.08, "K"),
                    ),
                    NASAPolynomial(
                        coeffs=[
                            2.78816,
                            0.000587644,
                            1.59009e-07,
                            -5.52736e-11,
                            4.34309e-15,
                            -596.143,
                            0.112747,
                        ],
                        Tmin=(1959.08, "K"),
                        Tmax=(5000, "K"),
                    ),
                ],
                Tmin=(100, "K"),
                Tmax=(5000, "K"),
                comment="""
Thermo library: primaryThermoLibrary
""",
            ),
            molecule=[Molecule(smiles="[H][H]")],
        )

        h = Species(
            index=3,
            label="H",
            thermo=NASA(
                polynomials=[
                    NASAPolynomial(
                        coeffs=[
                            2.5,
                            -1.91243e-12,
                            2.45329e-15,
                            -1.02377e-18,
                            1.31369e-22,
                            25474.2,
                            -0.444973,
                        ],
                        Tmin=(100, "K"),
                        Tmax=(4563.27, "K"),
                    ),
                    NASAPolynomial(
                        coeffs=[
                            2.50167,
                            -1.43051e-06,
                            4.6025e-10,
                            -6.57826e-14,
                            3.52412e-18,
                            25472.7,
                            -0.455578,
                        ],
                        Tmin=(4563.27, "K"),
                        Tmax=(5000, "K"),
                    ),
                ],
                Tmin=(100, "K"),
                Tmax=(5000, "K"),
                comment="""
Thermo library: primaryThermoLibrary
""",
            ),
            molecule=[Molecule(smiles="[H]")],
        )

        oh = Species(
            index=4,
            label="OH",
            thermo=NASA(
                polynomials=[
                    NASAPolynomial(
                        coeffs=[
                            3.51457,
                            2.92773e-05,
                            -5.32163e-07,
                            1.01949e-09,
                            -3.85945e-13,
                            3414.25,
                            2.10435,
                        ],
                        Tmin=(100, "K"),
                        Tmax=(1145.75, "K"),
                    ),
                    NASAPolynomial(
                        coeffs=[
                            3.07194,
                            0.000604016,
                            -1.39783e-08,
                            -2.13446e-11,
                            2.48066e-15,
                            3579.39,
                            4.578,
                        ],
                        Tmin=(1145.75, "K"),
                        Tmax=(5000, "K"),
                    ),
                ],
                Tmin=(100, "K"),
                Tmax=(5000, "K"),
                comment="""
Thermo library: primaryThermoLibrary
""",
            ),
            molecule=[Molecule(smiles="[OH]")],
        )

        ho2 = Species(
            index=5,
            label="HO2",
            thermo=NASA(
                polynomials=[
                    NASAPolynomial(
                        coeffs=[
                            4.04594,
                            -0.00173464,
                            1.03766e-05,
                            -1.02202e-08,
                            3.34908e-12,
                            -986.754,
                            4.63581,
                        ],
                        Tmin=(100, "K"),
                        Tmax=(932.15, "K"),
                    ),
                    NASAPolynomial(
                        coeffs=[
                            3.21024,
                            0.00367942,
                            -1.27701e-06,
                            2.18045e-10,
                            -1.46338e-14,
                            -910.369,
                            8.18291,
                        ],
                        Tmin=(932.15, "K"),
                        Tmax=(5000, "K"),
                    ),
                ],
                Tmin=(100, "K"),
                Tmax=(5000, "K"),
                comment="""
Thermo group additivity estimation: group(O2s-OsH) + gauche(O2s(RR)) + other(R) + group(O2s-OsH) + gauche(O2s(RR)) + other(R) + radical(HOOJ)
""",
            ),
            molecule=[Molecule(smiles="[O]O")],
        )

        o2 = Species(
            index=6,
            label="O2",
            thermo=NASA(
                polynomials=[
                    NASAPolynomial(
                        coeffs=[
                            3.53732,
                            -0.00121572,
                            5.3162e-06,
                            -4.89446e-09,
                            1.45846e-12,
                            -1038.59,
                            4.68368,
                        ],
                        Tmin=(100, "K"),
                        Tmax=(1074.55, "K"),
                    ),
                    NASAPolynomial(
                        coeffs=[
                            3.15382,
                            0.00167804,
                            -7.69974e-07,
                            1.51275e-10,
                            -1.08782e-14,
                            -1040.82,
                            6.16756,
                        ],
                        Tmin=(1074.55, "K"),
                        Tmax=(5000, "K"),
                    ),
                ],
                Tmin=(100, "K"),
                Tmax=(5000, "K"),
                comment="""
Thermo library: primaryThermoLibrary
""",
            ),
            molecule=[Molecule(smiles="[O][O]")],
        )

        co = Species(
            index=9,
            label="CO",
            thermo=NASA(
                polynomials=[
                    NASAPolynomial(
                        coeffs=[
                            3.66965,
                            -0.00550953,
                            2.00538e-05,
                            -2.08391e-08,
                            7.43738e-12,
                            1200.77,
                            -12.4224,
                        ],
                        Tmin=(100, "K"),
                        Tmax=(884.77, "K"),
                    ),
                    NASAPolynomial(
                        coeffs=[
                            2.8813,
                            0.00231665,
                            -4.40151e-07,
                            4.75633e-11,
                            -2.78282e-15,
                            1173.45,
                            -9.65831,
                        ],
                        Tmin=(884.77, "K"),
                        Tmax=(5000, "K"),
                    ),
                ],
                Tmin=(100, "K"),
                Tmax=(5000, "K"),
                comment="""
Thermo group additivity estimation: group(Ct-CtCs) + other(R) + group(O2s-CsCs) + other(R)
""",
            ),
            molecule=[Molecule(smiles="[C-]#[O+]")],
        )

        h2o2 = Species(
            index=7,
            label="H2O2",
            thermo=NASA(
                polynomials=[
                    NASAPolynomial(
                        coeffs=[
                            3.73136,
                            0.00335071,
                            9.35033e-06,
                            -1.521e-08,
                            6.41585e-12,
                            -17721.2,
                            5.45911,
                        ],
                        Tmin=(100, "K"),
                        Tmax=(908.87, "K"),
                    ),
                    NASAPolynomial(
                        coeffs=[
                            5.41579,
                            0.00261008,
                            -4.39892e-07,
                            4.91087e-11,
                            -3.35188e-15,
                            -18303,
                            -4.02248,
                        ],
                        Tmin=(908.87, "K"),
                        Tmax=(5000, "K"),
                    ),
                ],
                Tmin=(100, "K"),
                Tmax=(5000, "K"),
                comment="""
Thermo group additivity estimation: group(O2s-OsH) + gauche(O2s(RR)) + other(R) + group(O2s-OsH) + gauche(O2s(RR)) + other(R)
""",
            ),
            molecule=[Molecule(smiles="OO")],
        )

        self.species_list = [
            ch3,
            ethane,
            co2,
            ch4,
            h2o,
            ar,
            h2,
            h,
            oh,
            ho2,
            o2,
            co,
            h2o2,
        ]

        yaml_def = """
phases:
- name: gas
  thermo: ideal-gas
  kinetics: gas
  elements: [O, H, C, Ar]
  species: [ethane, Ar, H2(2), H(3), OH(4), HO2(5), O2(6), H2O2(7), CO(9), CH3(13), CH4(15), CO2(16), H2O(27)]
  state: {T: 300.0, P: 1 atm}
species:
- name: ethane
  composition: {C: 2, H: 6}
  thermo:
    model: NASA7
    temperature-ranges: [100.0, 954.51, 5000.0]
    data:
    - [3.78033, -3.24263e-03, 5.52381e-05, -6.38581e-08, 2.28637e-11, -1.16203e+04,
      5.21034]
    - [4.58983, 0.0141508, -4.75962e-06, 8.60294e-10, -6.21717e-14, -1.27218e+04,
      -3.61739]
- name: Ar
  composition: {Ar: 1}
  thermo:
    model: NASA7
    temperature-ranges: [200.0, 6000.0]
    data:
    - [2.5, 0.0, 0.0, 0.0, 0.0, -745.375, 4.37967]
- name: H2(2)
  composition: {H: 2}
  thermo:
    model: NASA7
    temperature-ranges: [100.0, 1959.08, 5000.0]
    data:
    - [3.43536, 2.1271e-04, -2.78625e-07, 3.40267e-10, -7.76031e-14, -1031.36,
      -3.90842]
    - [2.78816, 5.87644e-04, 1.59009e-07, -5.52736e-11, 4.34309e-15, -596.143,
      0.112747]
- name: H(3)
  composition: {H: 1}
  thermo:
    model: NASA7
    temperature-ranges: [100.0, 4563.27, 5000.0]
    data:
    - [2.5, -1.91243e-12, 2.45329e-15, -1.02377e-18, 1.31369e-22, 2.54742e+04,
      -0.444973]
    - [2.50167, -1.43051e-06, 4.6025e-10, -6.57826e-14, 3.52412e-18, 2.54727e+04,
      -0.455578]
- name: OH(4)
  composition: {H: 1, O: 1}
  thermo:
    model: NASA7
    temperature-ranges: [100.0, 1145.75, 5000.0]
    data:
    - [3.51457, 2.92773e-05, -5.32163e-07, 1.01949e-09, -3.85945e-13, 3414.25,
      2.10435]
    - [3.07194, 6.04016e-04, -1.39783e-08, -2.13446e-11, 2.48066e-15, 3579.39,
      4.578]
- name: HO2(5)
  composition: {H: 1, O: 2}
  thermo:
    model: NASA7
    temperature-ranges: [100.0, 932.15, 5000.0]
    data:
    - [4.04594, -1.73464e-03, 1.03766e-05, -1.02202e-08, 3.34908e-12, -986.754,
      4.63581]
    - [3.21024, 3.67942e-03, -1.27701e-06, 2.18045e-10, -1.46338e-14, -910.369,
      8.18291]
- name: O2(6)
  composition: {O: 2}
  thermo:
    model: NASA7
    temperature-ranges: [100.0, 1074.55, 5000.0]
    data:
    - [3.53732, -1.21572e-03, 5.3162e-06, -4.89446e-09, 1.45846e-12, -1038.59,
      4.68368]
    - [3.15382, 1.67804e-03, -7.69974e-07, 1.51275e-10, -1.08782e-14, -1040.82,
      6.16756]
- name: H2O2(7)
  composition: {H: 2, O: 2}
  thermo:
    model: NASA7
    temperature-ranges: [100.0, 908.87, 5000.0]
    data:
    - [3.73136, 3.35071e-03, 9.35033e-06, -1.521e-08, 6.41585e-12, -1.77212e+04,
      5.45911]
    - [5.41579, 2.61008e-03, -4.39892e-07, 4.91087e-11, -3.35188e-15, -1.8303e+04,
      -4.02248]
- name: CO(9)
  composition: {C: 1, O: 1}
  thermo:
    model: NASA7
    temperature-ranges: [100.0, 884.77, 5000.0]
    data:
    - [3.66965, -5.50953e-03, 2.00538e-05, -2.08391e-08, 7.43738e-12, 1200.77,
      -12.4224]
    - [2.8813, 2.31665e-03, -4.40151e-07, 4.75633e-11, -2.78282e-15, 1173.45,
      -9.65831]
- name: CH3(13)
  composition: {C: 1, H: 3}
  thermo:
    model: NASA7
    temperature-ranges: [100.0, 1337.62, 5000.0]
    data:
    - [3.91547, 1.84154e-03, 3.48744e-06, -3.3275e-09, 8.49964e-13, 1.62856e+04,
      0.351739]
    - [3.54145, 4.76788e-03, -1.82149e-06, 3.28878e-10, -2.22547e-14, 1.6224e+04,
      1.6604]
- name: CH4(15)
  composition: {C: 1, H: 4}
  thermo:
    model: NASA7
    temperature-ranges: [100.0, 1084.12, 5000.0]
    data:
    - [4.20541, -5.35556e-03, 2.51123e-05, -2.13762e-08, 5.97522e-12, -1.01619e+04,
      -0.921275]
    - [0.908272, 0.0114541, -4.57173e-06, 8.2919e-10, -5.66314e-14, -9719.98,
      13.9931]
- name: CO2(16)
  composition: {C: 1, O: 2}
  thermo:
    model: NASA7
    temperature-ranges: [100.0, 988.89, 5000.0]
    data:
    - [3.27861, 2.74152e-03, 7.16074e-06, -1.08027e-08, 4.14282e-12, -4.84703e+04,
      5.97937]
    - [4.5461, 2.91913e-03, -1.15484e-06, 2.27654e-10, -1.7091e-14, -4.89804e+04,
      -1.43275]
- name: H2O(27)
  composition: {H: 2, O: 1}
  thermo:
    model: NASA7
    temperature-ranges: [100.0, 1130.24, 5000.0]
    data:
    - [4.05764, -7.87933e-04, 2.90876e-06, -1.47518e-09, 2.12838e-13, -3.02816e+04,
      -0.311363]
    - [2.84325, 2.75108e-03, -7.8103e-07, 1.07243e-10, -5.79389e-15, -2.99586e+04,
      5.91041]
reactions: 
- equation: H(3) + H(3) + M <=> H2(2) + M
  type: three-body
  rate-constant:
    A: 1.000000e+18
    b: -1.0
    Ea: 0.0
  efficiencies:
    CO2(16): 0.0 
    CH4(15): 2.0 
    ethane: 3.0 
    H2O(27): 0.0 
    H2(2): 0.0 
    Ar: 0.63  
        """

        gas = ct.Solution(yaml=yaml_def)

        self.troe = Reaction(
            index=1,
            reactants=[ch3, ch3],
            products=[ethane],
            kinetics=Troe(
                arrheniusHigh=Arrhenius(
                    A=(6.77e16, "cm^3/(mol*s)"),
                    n=-1.18,
                    Ea=(0.654, "kcal/mol"),
                    T0=(1, "K"),
                ),
                arrheniusLow=Arrhenius(
                    A=(3.4e41, "cm^6/(mol^2*s)"),
                    n=-7.03,
                    Ea=(2.762, "kcal/mol"),
                    T0=(1, "K"),
                ),
                alpha=0.619,
                T3=(73.2, "K"),
                T1=(1180, "K"),
                T2=(10000, "K"),
                efficiencies={
                    Molecule(smiles="O=C=O"): 2.0,
                    Molecule(smiles="[H][H]"): 2.0,
                    Molecule(smiles="O"): 6.0,
                    Molecule(smiles="[Ar]"): 0.7,
                    Molecule(smiles="C"): 2.0,
                    Molecule(smiles="CC"): 3.0,
                },
            ),
        )

        self.ct_troe = ct.Reaction.from_yaml(
            """
            equation: CH3(13) + CH3(13) (+M) <=> ethane (+M)  # Reaction 1
            type: falloff
            low-P-rate-constant: {A: 3.4e+41 cm^6/mol^2/s, b: -7.03, Ea: 2.762 kcal/mol}
            high-P-rate-constant: {A: 6.77e+16 cm^3/mol/s, b: -1.18, Ea: 0.654 kcal/mol}
            Troe: {A: 0.619, T3: 73.2, T1: 1180.0, T2: 1.0e+04}
            efficiencies: {CH4(15): 2.0, Ar: 0.7, ethane: 3.0, H2(2): 2.0, CO2(16): 2.0,
            H2O(27): 6.0}""",
            gas,
        )

        self.arrheniusBi = Reaction(
            index=2,
            reactants=[h, ch4],
            products=[h2, ch3],
            kinetics=Arrhenius(A=(6.6e08, "cm^3/(mol*s)"), n=1.62, Ea=(10.84, "kcal/mol"), T0=(1, "K")),
        )

        self.ct_arrheniusBi = ct.Reaction.from_yaml(
            """ 
            equation: H(3) + CH4(15) <=> H2(2) + CH3(13)
            rate-constant: {A: 6.6e+08 cm^3/mol/s, b: 1.62, Ea: 10.84 kcal/mol}
        """,
            gas,
        )

        self.arrheniusBi_irreversible = Reaction(
            index=10,
            reactants=[h, ch4],
            products=[h2, ch3],
            kinetics=Arrhenius(A=(6.6e08, "cm^3/(mol*s)"), n=1.62, Ea=(10.84, "kcal/mol"), T0=(1, "K")),
            reversible=False,
        )

        self.ct_arrheniusBi_irreversible = ct.Reaction.from_yaml(
            """
            equation: H(3) + CH4(15) => H2(2) + CH3(13)  # Reaction 3
            rate-constant: {A: 6.6e+08 cm^3/mol/s, b: 1.62, Ea: 10.84 kcal/mol}
        """,
            gas,
        )

        self.arrheniusMono = Reaction(
            index=15,
            reactants=[h2o2],
            products=[h2, o2],
            kinetics=Arrhenius(A=(6.6e03, "1/s"), n=1.62, Ea=(10.84, "kcal/mol"), T0=(1, "K")),
        )

        self.ct_arrheniusMono = ct.Reaction.from_yaml(
            """
            equation: H2O2(7) <=> H2(2) + O2(6)  # Reaction 4
            rate-constant: {A: 6600.0 1/s, b: 1.62, Ea: 10.84 kcal/mol}
        """,
            gas,
        )

        self.arrheniusTri = Reaction(
            index=20,
            reactants=[h, h, o2],
            products=[h2o2],
            kinetics=Arrhenius(
                A=(6.6e08, "cm^6/(mol^2*s)"),
                n=1.62,
                Ea=(10.84, "kcal/mol"),
                T0=(1, "K"),
            ),
        )

        self.ct_arrheniusTri = ct.Reaction.from_yaml(
            """
            equation: H(3) + H(3) + O2(6) <=> H2O2(7)  # Reaction 5
            rate-constant: {A: 6.6e+08 cm^6/mol^2/s, b: 1.62, Ea: 10.84 kcal/mol}
        """,
            gas,
        )

        self.multiArrhenius = Reaction(
            index=3,
            reactants=[oh, ho2],
            products=[h2o, o2],
            kinetics=MultiArrhenius(
                arrhenius=[
                    Arrhenius(
                        A=(1.45e13, "cm^3/(mol*s)"),
                        n=0,
                        Ea=(-0.5, "kcal/mol"),
                        T0=(1, "K"),
                    ),
                    Arrhenius(
                        A=(5e15, "cm^3/(mol*s)"),
                        n=0,
                        Ea=(17.33, "kcal/mol"),
                        T0=(1, "K"),
                    ),
                ]
            ),
        )
        self.ct_multiArrhenius = [
            ct.Reaction.from_yaml(
                """
            equation: OH(4) + HO2(5) <=> H2O(27) + O2(6)  # Reaction 6
            duplicate: true
            rate-constant: {A: 1.45e+13 cm^3/mol/s, b: 0.0, Ea: -0.5 kcal/mol}
            """,
                gas,
            ),
            ct.Reaction.from_yaml(
                """
            equation: OH(4) + HO2(5) <=> H2O(27) + O2(6)  # Reaction 7
            duplicate: true
            rate-constant: {A: 5.0e+15 cm^3/mol/s, b: 0.0, Ea: 17.33 kcal/mol}
            """,
                gas,
            ),
        ]

        self.pdepArrhenius = Reaction(
            index=4,
            reactants=[ho2, ho2],
            products=[o2, h2o2],
            kinetics=PDepArrhenius(
                pressures=([0.1, 1, 10], "atm"),
                arrhenius=[
                    Arrhenius(
                        A=(8.8e16, "cm^3/(mol*s)"),
                        n=-1.05,
                        Ea=(6461, "cal/mol"),
                        T0=(1, "K"),
                    ),
                    Arrhenius(
                        A=(8e21, "cm^3/(mol*s)"),
                        n=-2.39,
                        Ea=(11180, "cal/mol"),
                        T0=(1, "K"),
                    ),
                    Arrhenius(
                        A=(3.3e24, "cm^3/(mol*s)"),
                        n=-3.04,
                        Ea=(15610, "cal/mol"),
                        T0=(1, "K"),
                    ),
                ],
            ),
        )

        self.ct_pdepArrhenius = ct.Reaction.from_yaml(
            """
            equation: HO2(5) + HO2(5) <=> O2(6) + H2O2(7)  # Reaction 8
            type: pressure-dependent-Arrhenius
            rate-constants:
            - {P: 0.1 atm, A: 8.8e+16 cm^3/mol/s, b: -1.05, Ea: 6.461 kcal/mol}
            - {P: 1.0 atm, A: 8.0e+21 cm^3/mol/s, b: -2.39, Ea: 11.18 kcal/mol}
            - {P: 10.0 atm, A: 3.3e+24 cm^3/mol/s, b: -3.04, Ea: 15.61 kcal/mol}
        """,
            gas,
        )

        self.multiPdepArrhenius = Reaction(
            index=5,
            reactants=[ho2, ch3],
            products=[o2, ch4],
            kinetics=MultiPDepArrhenius(
                arrhenius=[
                    PDepArrhenius(
                        pressures=([0.001, 1, 3], "atm"),
                        arrhenius=[
                            Arrhenius(
                                A=(9.3e10, "cm^3/(mol*s)"),
                                n=0,
                                Ea=(0, "cal/mol"),
                                T0=(1, "K"),
                            ),
                            Arrhenius(
                                A=(8e10, "cm^3/(mol*s)"),
                                n=0,
                                Ea=(0, "cal/mol"),
                                T0=(1, "K"),
                            ),
                            Arrhenius(
                                A=(7e10, "cm^3/(mol*s)"),
                                n=0,
                                Ea=(0, "cal/mol"),
                                T0=(1, "K"),
                            ),
                        ],
                    ),
                    PDepArrhenius(
                        pressures=([0.001, 1, 3], "atm"),
                        arrhenius=[
                            Arrhenius(
                                A=(710000, "cm^3/(mol*s)"),
                                n=1.8,
                                Ea=(1133, "cal/mol"),
                                T0=(1, "K"),
                            ),
                            Arrhenius(
                                A=(880000, "cm^3/(mol*s)"),
                                n=1.77,
                                Ea=(954, "cal/mol"),
                                T0=(1, "K"),
                            ),
                            Arrhenius(
                                A=(290000, "cm^3/(mol*s)"),
                                n=1.9,
                                Ea=(397, "cal/mol"),
                                T0=(1, "K"),
                            ),
                        ],
                    ),
                ],
            ),
        )
        self.ct_multiPdepArrhenius = [
            ct.Reaction.from_yaml(
                """
            equation: HO2(5) + CH3(13) <=> O2(6) + CH4(15)  # Reaction 9
            duplicate: true
            type: pressure-dependent-Arrhenius
            rate-constants:
            - {P: 0.001 atm, A: 9.3e+10 cm^3/mol/s, b: 0.0, Ea: 0.0 kcal/mol}
            - {P: 1.0 atm, A: 8.0e+10 cm^3/mol/s, b: 0.0, Ea: 0.0 kcal/mol}
            - {P: 3.0 atm, A: 7.0e+10 cm^3/mol/s, b: 0.0, Ea: 0.0 kcal/mol}
            """,
                gas,
            ),
            ct.Reaction.from_yaml(
                """
            equation: HO2(5) + CH3(13) <=> O2(6) + CH4(15)  # Reaction 10
            duplicate: true
            type: pressure-dependent-Arrhenius
            rate-constants:
            - {P: 0.001 atm, A: 7.1e+05 cm^3/mol/s, b: 1.8, Ea: 1.133 kcal/mol}
            - {P: 1.0 atm, A: 8.8e+05 cm^3/mol/s, b: 1.77, Ea: 0.954 kcal/mol}
            - {P: 3.0 atm, A: 2.9e+05 cm^3/mol/s, b: 1.9, Ea: 0.397 kcal/mol}
            """,
                gas,
            ),
        ]

        self.chebyshev = Reaction(
            index=6,
            reactants=[h, ch3],
            products=[ch4],
            kinetics=Chebyshev(
                coeffs=[
                    [12.68, 0.3961, -0.05481, -0.003606],
                    [-0.7128, 0.731, -0.0941, -0.008587],
                    [-0.5806, 0.57, -0.05539, -0.01115],
                    [-0.4074, 0.3653, -0.0118, -0.01171],
                    [-0.2403, 0.1779, 0.01946, -0.008505],
                    [-0.1133, 0.0485, 0.03121, -0.002955],
                ],
                kunits="cm^3/(mol*s)",
                Tmin=(300, "K"),
                Tmax=(3000, "K"),
                Pmin=(0.001, "atm"),
                Pmax=(98.692, "atm"),
            ),
        )

        self.ct_chebyshev = ct.Reaction.from_yaml(
            """
            equation: H(3) + CH3(13) <=> CH4(15)  # Reaction 11
            type: Chebyshev
            temperature-range: [300.0, 3000.0]
            pressure-range: [0.001 atm, 98.692 atm]
            data:
            - [9.680, 0.3961, -0.05481, -3.606e-03]
            - [-0.7128, 0.731, -0.0941, -8.587e-03]
            - [-0.5806, 0.57, -0.05539, -0.01115]
            - [-0.4074, 0.3653, -0.0118, -0.01171]
            - [-0.2403, 0.1779, 0.01946, -8.505e-03]
            - [-0.1133, 0.0485, 0.03121, -2.955e-03]
        """,
            gas,
        )

        self.thirdBody = Reaction(
            index=7,
            reactants=[h, h],
            products=[h2],
            kinetics=ThirdBody(
                arrheniusLow=Arrhenius(A=(1e18, "cm^6/(mol^2*s)"), n=-1, Ea=(0, "kcal/mol"), T0=(1, "K")),
                efficiencies={
                    Molecule(smiles="O=C=O"): 0.0,
                    Molecule(smiles="[H][H]"): 0.0,
                    Molecule(smiles="O"): 0.0,
                    Molecule(smiles="[Ar]"): 0.63,
                    Molecule(smiles="C"): 2.0,
                    Molecule(smiles="CC"): 3.0,
                },
            ),
        )

        self.ct_thirdBody = ct.Reaction.from_yaml(
            """
            equation: H(3) + H(3) + M <=> H2(2) + M  # Reaction 12
            type: three-body
            rate-constant: {A: 1.0e+18 cm^6/mol^2/s, b: -1.0, Ea: 0.0 kcal/mol}
            efficiencies: {H2(2): 0.0, CH4(15): 2.0, CO2(16): 0.0, Ar: 0.63, H2O(27): 0.0,
                ethane: 3.0}
        """,
            gas,
        )

        self.lindemann = Reaction(
            index=8,
            reactants=[h, o2],
            products=[ho2],
            kinetics=Lindemann(
                arrheniusHigh=Arrhenius(A=(1.8e10, "cm^3/(mol*s)"), n=0, Ea=(2.385, "kcal/mol"), T0=(1, "K")),
                arrheniusLow=Arrhenius(A=(6.02e14, "cm^6/(mol^2*s)"), n=0, Ea=(3, "kcal/mol"), T0=(1, "K")),
                efficiencies={
                    Molecule(smiles="O=C=O"): 3.5,
                    Molecule(smiles="[H][H]"): 2.0,
                    Molecule(smiles="O"): 6.0,
                    Molecule(smiles="[Ar]"): 0.5,
                    Molecule(smiles="C"): 2.0,
                    Molecule(smiles="CC"): 3.0,
                    Molecule(smiles="[O][O]"): 6.0,
                },
            ),
        )

        self.ct_lindemann = ct.Reaction.from_yaml(
            """
            equation: H(3) + O2(6) (+M) <=> HO2(5) (+M)  # Reaction 13
            type: falloff
            low-P-rate-constant: {A: 6.02e+14 cm^6/mol^2/s, b: 0.0, Ea: 3.0 kcal/mol}
            high-P-rate-constant: {A: 1.8e+10 cm^3/mol/s, b: 0.0, Ea: 2.385 kcal/mol}
            efficiencies: {H2O(27): 6.0, Ar: 0.5, CH4(15): 2.0, O2(6): 6.0, ethane: 3.0,
                H2(2): 2.0, CO2(16): 3.5}
        """,
            gas,
        )

    def test_arrhenius(self):
        """
        Tests formation of cantera reactions with Arrhenius or kinetics.
        """

        rmg_objects = [
            self.arrheniusBi,
            self.arrheniusBi_irreversible,
            self.arrheniusMono,
            self.arrheniusTri,
        ]

        ct_objects = [
            self.ct_arrheniusBi,
            self.ct_arrheniusBi_irreversible,
            self.ct_arrheniusMono,
            self.ct_arrheniusTri,
        ]
        converted_ct_objects = [obj.to_cantera(self.species_list, use_chemkin_identifier=True) for obj in rmg_objects]

        for converted_obj, ct_obj in zip(converted_ct_objects, ct_objects):
            # Check that the reaction class is the same
            assert type(converted_obj) == type(ct_obj)
            # Check that the reaction string is the same
            assert repr(converted_obj) == repr(ct_obj)
            # Check that the rate is the same. arrhenius string is not going to be identical
            assert converted_obj.rate.input_data == ct_obj.rate.input_data

    def test_multi_arrhenius(self):
        """
        Tests formation of cantera reactions with MultiArrhenius kinetics.
        """
        rmg_objects = [self.multiArrhenius]
        ct_objects = [self.ct_multiArrhenius]
        converted_ct_objects = [obj.to_cantera(self.species_list, use_chemkin_identifier=True) for obj in rmg_objects]

        for converted_obj, ct_obj in zip(converted_ct_objects, ct_objects):
            # Check that the same number of reactions are produced
            assert len(converted_obj) == len(ct_obj)

            for converted_rxn, ct_rxn in zip(converted_obj, ct_obj):
                # Check that the reaction has the same type
                assert type(converted_rxn) == type(ct_rxn)
                # Check that the reaction string is the same
                assert repr(converted_rxn) == repr(ct_rxn)
                # Check that the Arrhenius rates are identical
                assert round(abs(converted_rxn.rate.pre_exponential_factor - ct_rxn.rate.pre_exponential_factor), 3) == 0
                assert round(abs(converted_rxn.rate.temperature_exponent - ct_rxn.rate.temperature_exponent), 7) == 0
                assert round(abs(converted_rxn.rate.activation_energy - ct_rxn.rate.activation_energy), 7) == 0

    def test_pdep_arrhenius(self):
        """
        Tests formation of cantera reactions with PDepArrhenius kinetics.
        """
        rmg_objects = [self.pdepArrhenius]
        ct_objects = [self.ct_pdepArrhenius]
        converted_ct_objects = [obj.to_cantera(self.species_list, use_chemkin_identifier=True) for obj in rmg_objects]

        for converted_obj, ct_obj in zip(converted_ct_objects, ct_objects):
            # Check that the reaction class is the same
            assert type(converted_obj) == type(ct_obj)
            # Check that the reaction string is the same
            assert repr(converted_obj) == repr(ct_obj)
            # Check that the Arrhenius rates are identical
            assert str(converted_obj.rates) == str(ct_obj.rates)

    def test_multi_pdep_arrhenius(self):
        """
        Tests formation of cantera reactions with MultiPDepArrhenius kinetics.
        """

        rmg_objects = [self.multiPdepArrhenius]
        ct_objects = [self.ct_multiPdepArrhenius]
        converted_ct_objects = [obj.to_cantera(self.species_list, use_chemkin_identifier=True) for obj in rmg_objects]

        for converted_obj, ct_obj in zip(converted_ct_objects, ct_objects):
            # Check that the same number of reactions are produced
            assert len(converted_obj) == len(ct_obj)

            for converted_rxn, ct_rxn in zip(converted_obj, ct_obj):
                # Check that the reaction has the same type
                assert type(converted_rxn) == type(ct_rxn)
                # Check that the reaction string is the same
                assert repr(converted_rxn) == repr(ct_rxn)
                # Check that the Arrhenius rates are identical
                assert str(converted_rxn.rates) == str(ct_rxn.rates)

    def test_chebyshev(self):
        """
        Tests formation of cantera reactions with Chebyshev kinetics.
        """
        ct_chebyshev = self.chebyshev.to_cantera(self.species_list, use_chemkin_identifier=True)
        assert type(ct_chebyshev.rate) == type(self.ct_chebyshev.rate)
        assert ct_chebyshev.rate.temperature_range == self.ct_chebyshev.rate.temperature_range
        assert ct_chebyshev.rate.pressure_range == self.ct_chebyshev.rate.pressure_range
        assert (ct_chebyshev.rate.data == self.ct_chebyshev.rate.data).all()

    def test_falloff(self):
        """
        Tests formation of cantera reactions with Falloff kinetics.
        """
        ct_troe = self.troe.to_cantera(self.species_list, use_chemkin_identifier=True)
        assert type(ct_troe.rate) == type(self.ct_troe.rate)
        assert round(abs(ct_troe.rate.low_rate.pre_exponential_factor - self.ct_troe.rate.low_rate.pre_exponential_factor), 3) == 0
        assert ct_troe.rate.low_rate.temperature_exponent == self.ct_troe.rate.low_rate.temperature_exponent
        assert ct_troe.rate.low_rate.activation_energy == self.ct_troe.rate.low_rate.activation_energy
        assert ct_troe.efficiencies == self.ct_troe.efficiencies

        ct_third_body = self.thirdBody.to_cantera(self.species_list, use_chemkin_identifier=True)
        assert type(ct_third_body.rate) == type(self.ct_thirdBody.rate)
        assert round(abs(ct_third_body.rate.pre_exponential_factor - self.ct_thirdBody.rate.pre_exponential_factor), 3) == 0
        assert ct_third_body.rate.temperature_exponent == self.ct_thirdBody.rate.temperature_exponent
        assert ct_third_body.rate.activation_energy == self.ct_thirdBody.rate.activation_energy
        assert ct_third_body.efficiencies == self.ct_thirdBody.efficiencies

        ct_lindemann = self.lindemann.to_cantera(self.species_list, use_chemkin_identifier=True)
        assert type(ct_lindemann.rate) == type(self.ct_lindemann.rate)
        assert ct_lindemann.efficiencies == self.ct_lindemann.efficiencies
        assert str(ct_lindemann.rate.low_rate) == str(self.ct_lindemann.rate.low_rate)
        assert str(ct_lindemann.rate.high_rate) == str(self.ct_lindemann.rate.high_rate)
