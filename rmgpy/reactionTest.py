#!/usr/bin/env python
# encoding: utf-8 -*-

"""
This module contains unit tests of the rmgpy.reaction module.
"""

import numpy
import unittest
from external.wip import work_in_progress

from rmgpy.species import Species, TransitionState
from rmgpy.reaction import Reaction
from rmgpy.statmech.translation import Translation, IdealGasTranslation
from rmgpy.statmech.rotation import Rotation, LinearRotor, NonlinearRotor, KRotor, SphericalTopRotor
from rmgpy.statmech.vibration import Vibration, HarmonicOscillator
from rmgpy.statmech.torsion import Torsion, HinderedRotor
from rmgpy.statmech.conformer import Conformer
from rmgpy.kinetics import Arrhenius
from rmgpy.thermo import Wilhoit
import rmgpy.constants as constants

################################################################################

class PseudoSpecies:
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
    def isIsomorphic(self, other):
        return self.label.lower() == other.label.lower()

class TestReactionIsomorphism(unittest.TestCase):
    """
    Contains unit tests of the isomorphism testing of the  Reaction class.
    """

    def makeReaction(self,reaction_string):
        """"
        Make a Reaction (containing PseudoSpecies) of from a string like 'Ab=CD'
        """
        reactants, products = reaction_string.split('=')
        reactants = [PseudoSpecies(i) for i in reactants]
        products = [PseudoSpecies(i) for i in products]
        return Reaction(reactants=reactants, products=products)
    def test1to1(self):
        r1 = self.makeReaction('A=B')
        self.assertTrue(r1.isIsomorphic(self.makeReaction('a=B')))
        self.assertTrue(r1.isIsomorphic(self.makeReaction('b=A')))
        self.assertFalse(r1.isIsomorphic(self.makeReaction('B=a'),eitherDirection=False))
        self.assertFalse(r1.isIsomorphic(self.makeReaction('A=C')))
        self.assertFalse(r1.isIsomorphic(self.makeReaction('A=BB')))
    def test1to2(self):
        r1 = self.makeReaction('A=BC')
        self.assertTrue(r1.isIsomorphic(self.makeReaction('a=Bc')))
        self.assertTrue(r1.isIsomorphic(self.makeReaction('cb=a')))
        self.assertTrue(r1.isIsomorphic(self.makeReaction('a=cb'),eitherDirection=False))
        self.assertFalse(r1.isIsomorphic(self.makeReaction('bc=a'),eitherDirection=False))
        self.assertFalse(r1.isIsomorphic(self.makeReaction('a=c')))
        self.assertFalse(r1.isIsomorphic(self.makeReaction('ab=c')))
    def test2to2(self):
        r1 = self.makeReaction('AB=CD')
        self.assertTrue(r1.isIsomorphic(self.makeReaction('ab=cd')))
        self.assertTrue(r1.isIsomorphic(self.makeReaction('ab=dc'),eitherDirection=False))
        self.assertTrue(r1.isIsomorphic(self.makeReaction('dc=ba')))
        self.assertFalse(r1.isIsomorphic(self.makeReaction('cd=ab'),eitherDirection=False))
        self.assertFalse(r1.isIsomorphic(self.makeReaction('ab=ab')))
        self.assertFalse(r1.isIsomorphic(self.makeReaction('ab=cde')))
    def test2to3(self):
        r1 = self.makeReaction('AB=CDE')
        self.assertTrue(r1.isIsomorphic(self.makeReaction('ab=cde')))
        self.assertTrue(r1.isIsomorphic(self.makeReaction('ba=edc'),eitherDirection=False))
        self.assertTrue(r1.isIsomorphic(self.makeReaction('dec=ba')))
        self.assertFalse(r1.isIsomorphic(self.makeReaction('cde=ab'),eitherDirection=False))
        self.assertFalse(r1.isIsomorphic(self.makeReaction('ab=abc')))
        self.assertFalse(r1.isIsomorphic(self.makeReaction('abe=cde')))


class TestReaction(unittest.TestCase):
    """
    Contains unit tests of the Reaction class.
    """
    
    def setUp(self):
        """
        A method that is called prior to each unit test in this class.
        """
        ethylene = Species(
            label = 'C2H4',
            conformer = Conformer(
                E0 = (44.7127, 'kJ/mol'),
                modes = [
                    IdealGasTranslation(
                        mass = (28.0313, 'amu'),
                    ),
                    NonlinearRotor(
                        inertia = (
                            [3.41526, 16.6498, 20.065],
                            'amu*angstrom^2',
                        ),
                        symmetry = 4,
                    ),
                    HarmonicOscillator(
                        frequencies = (
                            [828.397, 970.652, 977.223, 1052.93, 1233.55, 1367.56, 1465.09, 1672.25, 3098.46, 3111.7, 3165.79, 3193.54],
                            'cm^-1',
                        ),
                    ),
                ],
                spinMultiplicity = 1,
                opticalIsomers = 1,
            ),
        )
        
        hydrogen = Species(          
            label = 'H',
            conformer = Conformer(
                E0 = (211.794, 'kJ/mol'),
                modes = [
                    IdealGasTranslation(
                        mass = (1.00783, 'amu'),
                    ),
                ],
                spinMultiplicity = 2,
                opticalIsomers = 1,
            ),
        )
        
        ethyl = Species(
            label = 'C2H5',
            conformer = Conformer(
                E0 = (111.603, 'kJ/mol'),
                modes = [
                    IdealGasTranslation(
                        mass = (29.0391, 'amu'),
                    ),
                    NonlinearRotor(
                        inertia = (
                            [4.8709, 22.2353, 23.9925],
                            'amu*angstrom^2',
                        ),
                        symmetry = 1,
                    ),
                    HarmonicOscillator(
                        frequencies = (
                            [482.224, 791.876, 974.355, 1051.48, 1183.21, 1361.36, 1448.65, 1455.07, 1465.48, 2688.22, 2954.51, 3033.39, 3101.54, 3204.73],
                            'cm^-1',
                        ),
                    ),
                    HinderedRotor(
                        inertia = (1.11481, 'amu*angstrom^2'),
                        symmetry = 6,
                        barrier = (0.244029, 'kJ/mol'),
                        semiclassical = None,
                    ),
                ],
                spinMultiplicity = 2,
                opticalIsomers = 1,
            ),
        )
        
        TS = TransitionState(
            label = 'TS',
            conformer = Conformer(
                E0 = (266.694, 'kJ/mol'),
                modes = [
                    IdealGasTranslation(
                        mass = (29.0391, 'amu'),
                    ),
                    NonlinearRotor(
                        inertia = (
                            [6.78512, 22.1437, 22.2114],
                            'amu*angstrom^2',
                        ),
                        symmetry = 1,
                    ),
                    HarmonicOscillator(
                        frequencies = (
                            [412.75, 415.206, 821.495, 924.44, 982.714, 1024.16, 1224.21, 1326.36, 1455.06, 1600.35, 3101.46, 3110.55, 3175.34, 3201.88],
                            'cm^-1',
                        ),
                    ),
                ],
                spinMultiplicity = 2,
                opticalIsomers = 1,
            ),
            frequency = (-750.232, 'cm^-1'),
        )
        
        self.reaction = Reaction(
            reactants = [hydrogen, ethylene],
            products = [ethyl], 
            kinetics = Arrhenius(
                A = (501366000.0, 'cm^3/(mol*s)'),
                n = 1.637,
                Ea = (4.32508, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (2500, 'K'),
            ),
            transitionState = TS,
        )
    
        # CC(=O)O[O]
        acetylperoxy = Species(
            label='acetylperoxy',
            thermo=Wilhoit(Cp0=(4.0*constants.R,"J/(mol*K)"), CpInf=(21.0*constants.R,"J/(mol*K)"), a0=-3.95, a1=9.26, a2=-15.6, a3=8.55, B=(500.0,"K"), H0=(-6.151e+04,"J/mol"), S0=(-790.2,"J/(mol*K)")),
        )

        # C[C]=O
        acetyl = Species(
            label='acetyl',
            thermo=Wilhoit(Cp0=(4.0*constants.R,"J/(mol*K)"), CpInf=(15.5*constants.R,"J/(mol*K)"), a0=0.2541, a1=-0.4712, a2=-4.434, a3=2.25, B=(500.0,"K"), H0=(-1.439e+05,"J/mol"), S0=(-524.6,"J/(mol*K)")),
        )

        # [O][O]
        oxygen = Species(
            label='oxygen',
            thermo=Wilhoit(Cp0=(3.5*constants.R,"J/(mol*K)"), CpInf=(4.5*constants.R,"J/(mol*K)"), a0=-0.9324, a1=26.18, a2=-70.47, a3=44.12, B=(500.0,"K"), H0=(1.453e+04,"J/mol"), S0=(-12.19,"J/(mol*K)")),
        )
        
        self.reaction2 = Reaction(
            reactants=[acetyl, oxygen], 
            products=[acetylperoxy], 
            kinetics = Arrhenius(
                A = (2.65e12, 'cm^3/(mol*s)'),
                n = 0.0,
                Ea = (0.0, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (2000, 'K'),
            ),
        )
        
    def testIsIsomerization(self):
        """
        Test the Reaction.isIsomerization() method.
        """
        isomerization = Reaction(reactants=[Species()], products=[Species()])
        association = Reaction(reactants=[Species(),Species()], products=[Species()])
        dissociation = Reaction(reactants=[Species()], products=[Species(),Species()])
        bimolecular = Reaction(reactants=[Species(),Species()], products=[Species(),Species()])
        self.assertTrue(isomerization.isIsomerization())
        self.assertFalse(association.isIsomerization())
        self.assertFalse(dissociation.isIsomerization())
        self.assertFalse(bimolecular.isIsomerization())
        
    def testIsAssociation(self):
        """
        Test the Reaction.isAssociation() method.
        """
        isomerization = Reaction(reactants=[Species()], products=[Species()])
        association = Reaction(reactants=[Species(),Species()], products=[Species()])
        dissociation = Reaction(reactants=[Species()], products=[Species(),Species()])
        bimolecular = Reaction(reactants=[Species(),Species()], products=[Species(),Species()])
        self.assertFalse(isomerization.isAssociation())
        self.assertTrue(association.isAssociation())
        self.assertFalse(dissociation.isAssociation())
        self.assertFalse(bimolecular.isAssociation())

    def testIsDissociation(self):
        """
        Test the Reaction.isDissociation() method.
        """
        isomerization = Reaction(reactants=[Species()], products=[Species()])
        association = Reaction(reactants=[Species(),Species()], products=[Species()])
        dissociation = Reaction(reactants=[Species()], products=[Species(),Species()])
        bimolecular = Reaction(reactants=[Species(),Species()], products=[Species(),Species()])
        self.assertFalse(isomerization.isDissociation())
        self.assertFalse(association.isDissociation())
        self.assertTrue(dissociation.isDissociation())
        self.assertFalse(bimolecular.isDissociation())

    def testHasTemplate(self):
        """
        Test the Reaction.hasTemplate() method.
        """
        reactants = self.reaction.reactants[:]
        products = self.reaction.products[:]
        self.assertTrue(self.reaction.hasTemplate(reactants, products))
        self.assertTrue(self.reaction.hasTemplate(products, reactants))
        self.assertFalse(self.reaction2.hasTemplate(reactants, products))
        self.assertFalse(self.reaction2.hasTemplate(products, reactants))
        
        reactants.reverse()
        products.reverse()
        self.assertTrue(self.reaction.hasTemplate(reactants, products))
        self.assertTrue(self.reaction.hasTemplate(products, reactants))
        self.assertFalse(self.reaction2.hasTemplate(reactants, products))
        self.assertFalse(self.reaction2.hasTemplate(products, reactants))
        
        reactants = self.reaction2.reactants[:]
        products = self.reaction2.products[:]
        self.assertFalse(self.reaction.hasTemplate(reactants, products))
        self.assertFalse(self.reaction.hasTemplate(products, reactants))
        self.assertTrue(self.reaction2.hasTemplate(reactants, products))
        self.assertTrue(self.reaction2.hasTemplate(products, reactants))
        
        reactants.reverse()
        products.reverse()
        self.assertFalse(self.reaction.hasTemplate(reactants, products))
        self.assertFalse(self.reaction.hasTemplate(products, reactants))
        self.assertTrue(self.reaction2.hasTemplate(reactants, products))
        self.assertTrue(self.reaction2.hasTemplate(products, reactants))

    def testEnthalpyOfReaction(self):
        """
        Test the Reaction.getEnthalpyOfReaction() method.
        """
        Tlist = numpy.arange(200.0, 2001.0, 200.0, numpy.float64)
        Hlist0 = [float(v) for v in ['-146007', '-145886', '-144195', '-141973', '-139633', '-137341', '-135155', '-133093', '-131150', '-129316']]
        Hlist = self.reaction2.getEnthalpiesOfReaction(Tlist)
        for i in range(len(Tlist)):
            self.assertAlmostEqual(Hlist[i] / 1000., Hlist0[i] / 1000., 2)

    def testEntropyOfReaction(self):
        """
        Test the Reaction.getEntropyOfReaction() method.
        """
        Tlist = numpy.arange(200.0, 2001.0, 200.0, numpy.float64)
        Slist0 = [float(v) for v in ['-156.793', '-156.872', '-153.504', '-150.317', '-147.707', '-145.616', '-143.93', '-142.552', '-141.407', '-140.441']]
        Slist = self.reaction2.getEntropiesOfReaction(Tlist)
        for i in range(len(Tlist)):
            self.assertAlmostEqual(Slist[i], Slist0[i], 2)

    def testFreeEnergyOfReaction(self):
        """
        Test the Reaction.getFreeEnergyOfReaction() method.
        """
        Tlist = numpy.arange(200.0, 2001.0, 200.0, numpy.float64)
        Glist0 = [float(v) for v in ['-114648', '-83137.2', '-52092.4', '-21719.3', '8073.53', '37398.1', '66346.8', '94990.6', '123383', '151565']]
        Glist = self.reaction2.getFreeEnergiesOfReaction(Tlist)
        for i in range(len(Tlist)):
            self.assertAlmostEqual(Glist[i] / 1000., Glist0[i] / 1000., 2)

    def testEquilibriumConstantKa(self):
        """
        Test the Reaction.getEquilibriumConstant() method.
        """
        Tlist = numpy.arange(200.0, 2001.0, 200.0, numpy.float64)
        Kalist0 = [float(v) for v in ['8.75951e+29', '7.1843e+10', '34272.7', '26.1877', '0.378696', '0.0235579', '0.00334673', '0.000792389', '0.000262777', '0.000110053']]
        Kalist = self.reaction2.getEquilibriumConstants(Tlist, type='Ka')
        for i in range(len(Tlist)):
            self.assertAlmostEqual(Kalist[i] / Kalist0[i], 1.0, 4)

    def testEquilibriumConstantKc(self):
        """
        Test the Reaction.getEquilibriumConstant() method.
        """
        Tlist = numpy.arange(200.0, 2001.0, 200.0, numpy.float64)
        Kclist0 = [float(v) for v in ['1.45661e+28', '2.38935e+09', '1709.76', '1.74189', '0.0314866', '0.00235045', '0.000389568', '0.000105413', '3.93273e-05', '1.83006e-05']]
        Kclist = self.reaction2.getEquilibriumConstants(Tlist, type='Kc')
        for i in range(len(Tlist)):
            self.assertAlmostEqual(Kclist[i] / Kclist0[i], 1.0, 4)

    def testEquilibriumConstantKp(self):
        """
        Test the Reaction.getEquilibriumConstant() method.
        """
        Tlist = numpy.arange(200.0, 2001.0, 200.0, numpy.float64)
        Kplist0 = [float(v) for v in ['8.75951e+24', '718430', '0.342727', '0.000261877', '3.78696e-06', '2.35579e-07', '3.34673e-08', '7.92389e-09', '2.62777e-09', '1.10053e-09']]
        Kplist = self.reaction2.getEquilibriumConstants(Tlist, type='Kp')
        for i in range(len(Tlist)):
            self.assertAlmostEqual(Kplist[i] / Kplist0[i], 1.0, 4)

    def testStoichiometricCoefficient(self):
        """
        Test the Reaction.getStoichiometricCoefficient() method.
        """
        for reactant in self.reaction.reactants:
            self.assertEqual(self.reaction.getStoichiometricCoefficient(reactant), -1)
        for product in self.reaction.products:
            self.assertEqual(self.reaction.getStoichiometricCoefficient(product), 1)
        for reactant in self.reaction2.reactants:
            self.assertEqual(self.reaction.getStoichiometricCoefficient(reactant), 0)
        for product in self.reaction2.products:
            self.assertEqual(self.reaction.getStoichiometricCoefficient(product), 0)

    def testRateCoefficient(self):
        """
        Test the Reaction.getRateCoefficient() method.
        """
        Tlist = numpy.arange(200.0, 2001.0, 200.0, numpy.float64)
        P = 1e5
        for T in Tlist:
            self.assertAlmostEqual(self.reaction.getRateCoefficient(T, P) / self.reaction.kinetics.getRateCoefficient(T), 1.0, 6)
    
    def testGenerateReverseRateCoefficient(self):
        """
        Test the Reaction.generateReverseRateCoefficient() method.
        """
        Tlist = numpy.arange(200.0, 2001.0, 200.0, numpy.float64)
        P = 1e5
        reverseKinetics = self.reaction2.generateReverseRateCoefficient()
        for T in Tlist:
            kr0 = self.reaction2.getRateCoefficient(T, P) / self.reaction2.getEquilibriumConstant(T)
            kr = reverseKinetics.getRateCoefficient(T)
            self.assertAlmostEqual(kr0 / kr, 1.0, 0)

    def testGenerateReverseRateCoefficientArrhenius(self):
        """
        Test the Reaction.generateReverseRateCoefficient() method works for the Arrhenius format.
        """
        original_kinetics = Arrhenius(
                    A = (2.65e12, 'cm^3/(mol*s)'),
                    n = 0.0,
                    Ea = (0.0, 'kJ/mol'),
                    T0 = (1, 'K'),
                    Tmin = (300, 'K'),
                    Tmax = (2000, 'K'),
                )
        self.reaction2.kinetics = original_kinetics

        reverseKinetics = self.reaction2.generateReverseRateCoefficient()

        self.reaction2.kinetics = reverseKinetics
        # reverse reactants, products to ensure Keq is correctly computed
        self.reaction2.reactants, self.reaction2.products = self.reaction2.products, self.reaction2.reactants
        reversereverseKinetics = self.reaction2.generateReverseRateCoefficient()

        # check that reverting the reverse yields the original
        Tlist = numpy.arange(original_kinetics.Tmin.value_si, original_kinetics.Tmax.value_si, 200.0, numpy.float64)
        P = 1e5
        for T in Tlist:
            korig = original_kinetics.getRateCoefficient(T, P)
            krevrev = reversereverseKinetics.getRateCoefficient(T, P)
            self.assertAlmostEqual(korig / krevrev, 1.0, 0)

    @work_in_progress
    def testGenerateReverseRateCoefficientArrheniusEP(self):
        """
        Test the Reaction.generateReverseRateCoefficient() method works for the ArrheniusEP format.
        """
        from rmgpy.kinetics import ArrheniusEP

        original_kinetics = ArrheniusEP(
                    A = (2.65e12, 'cm^3/(mol*s)'),
                    n = 0.0,
                    alpha = 0.5,
                    E0 = (41.84, 'kJ/mol'),
                    Tmin = (300, 'K'),
                    Tmax = (2000, 'K'),
                )
        self.reaction2.kinetics = original_kinetics

        reverseKinetics = self.reaction2.generateReverseRateCoefficient()

        self.reaction2.kinetics = reverseKinetics
        # reverse reactants, products to ensure Keq is correctly computed
        self.reaction2.reactants, self.reaction2.products = self.reaction2.products, self.reaction2.reactants
        reversereverseKinetics = self.reaction2.generateReverseRateCoefficient()

        # check that reverting the reverse yields the original
        Tlist = numpy.arange(original_kinetics.Tmin, original_kinetics.Tmax, 200.0, numpy.float64)
        P = 1e5
        for T in Tlist:
            korig = original_kinetics.getRateCoefficient(T, P)
            krevrev = reversereverseKinetics.getRateCoefficient(T, P)
            self.assertAlmostEqual(korig / krevrev, 1.0, 0)

    def testGenerateReverseRateCoefficientPDepArrhenius(self):
        """
        Test the Reaction.generateReverseRateCoefficient() method works for the PDepArrhenius format.
        """
        from rmgpy.kinetics import PDepArrhenius

        arrhenius0 = Arrhenius(
            A = (1.0e6,"s^-1"),
            n = 1.0, 
            Ea = (10.0,"kJ/mol"), 
            T0 = (300.0,"K"), 
            Tmin = (300.0,"K"), 
            Tmax = (2000.0,"K"), 
            comment = """This data is completely made up""",
        )

        arrhenius1 = Arrhenius(
            A = (1.0e12,"s^-1"), 
            n = 1.0, 
            Ea = (20.0,"kJ/mol"), 
            T0 = (300.0,"K"), 
            Tmin = (300.0,"K"), 
            Tmax = (2000.0,"K"), 
            comment = """This data is completely made up""",
        )

        pressures = numpy.array([0.1, 10.0])
        arrhenius = [arrhenius0, arrhenius1]
        Tmin = 300.0
        Tmax = 2000.0
        Pmin = 0.1
        Pmax = 10.0
        comment = """This data is completely made up"""

        original_kinetics = PDepArrhenius(
            pressures = (pressures,"bar"),
            arrhenius = arrhenius,
            Tmin = (Tmin,"K"), 
            Tmax = (Tmax,"K"), 
            Pmin = (Pmin,"bar"), 
            Pmax = (Pmax,"bar"),
            comment = comment,
        )

        self.reaction2.kinetics = original_kinetics

        reverseKinetics = self.reaction2.generateReverseRateCoefficient()

        self.reaction2.kinetics = reverseKinetics
        # reverse reactants, products to ensure Keq is correctly computed
        self.reaction2.reactants, self.reaction2.products = self.reaction2.products, self.reaction2.reactants
        reversereverseKinetics = self.reaction2.generateReverseRateCoefficient()

        # check that reverting the reverse yields the original
        Tlist = numpy.arange(Tmin, Tmax, 200.0, numpy.float64)
        P = 1e5
        for T in Tlist:
            korig = original_kinetics.getRateCoefficient(T, P)
            krevrev = reversereverseKinetics.getRateCoefficient(T, P)
            self.assertAlmostEqual(korig / krevrev, 1.0, 0)


    def testGenerateReverseRateCoefficientMultiArrhenius(self):
        """
        Test the Reaction.generateReverseRateCoefficient() method works for the MultiArrhenius format.
        """
        from rmgpy.kinetics import MultiArrhenius

        pressures = numpy.array([0.1, 10.0])
        Tmin = 300.0
        Tmax = 2000.0
        Pmin = 0.1
        Pmax = 10.0
        comment = """This data is completely made up"""

        arrhenius = [
            Arrhenius(
                A = (9.3e-14,"cm^3/(molecule*s)"),
                n = 0.0,
                Ea = (4740*constants.R*0.001,"kJ/mol"),
                T0 = (1,"K"),
                Tmin = (Tmin,"K"),
                Tmax = (Tmax,"K"),
                comment = comment,
            ),
            Arrhenius(
                A = (1.4e-9,"cm^3/(molecule*s)"),
                n = 0.0,
                Ea = (11200*constants.R*0.001,"kJ/mol"),
                T0 = (1,"K"),
                Tmin = (Tmin,"K"),
                Tmax = (Tmax,"K"),
                comment = comment,
            ),
        ]

        original_kinetics = MultiArrhenius(
            arrhenius = arrhenius,
            Tmin = (Tmin,"K"),
            Tmax = (Tmax,"K"),
            comment = comment,
        )

        self.reaction2.kinetics = original_kinetics

        reverseKinetics = self.reaction2.generateReverseRateCoefficient()

        self.reaction2.kinetics = reverseKinetics
        # reverse reactants, products to ensure Keq is correctly computed
        self.reaction2.reactants, self.reaction2.products = self.reaction2.products, self.reaction2.reactants
        reversereverseKinetics = self.reaction2.generateReverseRateCoefficient()

        # check that reverting the reverse yields the original
        Tlist = numpy.arange(Tmin, Tmax, 200.0, numpy.float64)
        P = 1e5
        for T in Tlist:
            korig = original_kinetics.getRateCoefficient(T, P)
            krevrev = reversereverseKinetics.getRateCoefficient(T, P)
            self.assertAlmostEqual(korig / krevrev, 1.0, 0)
    def testTSTCalculation(self):
        """
        A test of the transition state theory k(T) calculation function,
        using the reaction H + C2H4 -> C2H5.
        """
        Tlist = 1000.0/numpy.arange(0.4, 3.35, 0.01)
        klist = numpy.array([self.reaction.calculateTSTRateCoefficient(T) for T in Tlist])
        arrhenius = Arrhenius().fitToData(Tlist, klist, kunits='m^3/(mol*s)')
        klist2 = numpy.array([arrhenius.getRateCoefficient(T) for T in Tlist])
        
        # Check that the correct Arrhenius parameters are returned
        self.assertAlmostEqual(arrhenius.A.value_si, 2265.2488, delta=1e-2)
        self.assertAlmostEqual(arrhenius.n.value_si, 1.45419, delta=1e-4)
        self.assertAlmostEqual(arrhenius.Ea.value_si, 6645.24, delta=1e-2)
        # Check that the fit is satisfactory (defined here as always within 5%)
        for i in range(len(Tlist)):
            self.assertAlmostEqual(klist[i], klist2[i], delta=5e-2 * klist[i])
        
    def testPickle(self):
        """
        Test that a Reaction object can be successfully pickled and
        unpickled with no loss of information.
        """
        import cPickle
        reaction = cPickle.loads(cPickle.dumps(self.reaction,-1))

        self.assertEqual(len(self.reaction.reactants), len(reaction.reactants))
        self.assertEqual(len(self.reaction.products), len(reaction.products))
        for reactant0, reactant in zip(self.reaction.reactants, reaction.reactants):
            self.assertAlmostEqual(reactant0.conformer.E0.value_si / 1e6, reactant.conformer.E0.value_si / 1e6, 2)
            self.assertEqual(reactant0.conformer.E0.units, reactant.conformer.E0.units)
        for product0, product in zip(self.reaction.products, reaction.products):
            self.assertAlmostEqual(product0.conformer.E0.value_si / 1e6, product.conformer.E0.value_si / 1e6, 2)
            self.assertEqual(product0.conformer.E0.units, product.conformer.E0.units)
        self.assertAlmostEqual(self.reaction.transitionState.conformer.E0.value_si / 1e6, reaction.transitionState.conformer.E0.value_si / 1e6, 2)
        self.assertEqual(self.reaction.transitionState.conformer.E0.units, reaction.transitionState.conformer.E0.units)
        self.assertAlmostEqual(self.reaction.transitionState.frequency.value_si, reaction.transitionState.frequency.value_si, 2)
        self.assertEqual(self.reaction.transitionState.frequency.units, reaction.transitionState.frequency.units)
        
        self.assertAlmostEqual(self.reaction.kinetics.A.value_si, reaction.kinetics.A.value_si, delta=1e-6)
        self.assertAlmostEqual(self.reaction.kinetics.n.value_si, reaction.kinetics.n.value_si, delta=1e-6)
        self.assertAlmostEqual(self.reaction.kinetics.T0.value_si, reaction.kinetics.T0.value_si, delta=1e-6)
        self.assertAlmostEqual(self.reaction.kinetics.Ea.value_si, reaction.kinetics.Ea.value_si, delta=1e-6)
        self.assertEqual(self.reaction.kinetics.comment, reaction.kinetics.comment)

        self.assertEqual(self.reaction.duplicate, reaction.duplicate)
        self.assertEqual(self.reaction.degeneracy, reaction.degeneracy)        

    def testOutput(self):
        """
        Test that a Reaction object can be successfully reconstructed
        from its repr() output with no loss of information.
        """
        exec('reaction = %r' % (self.reaction))

        self.assertEqual(len(self.reaction.reactants), len(reaction.reactants))
        self.assertEqual(len(self.reaction.products), len(reaction.products))
        for reactant0, reactant in zip(self.reaction.reactants, reaction.reactants):
            self.assertAlmostEqual(reactant0.conformer.E0.value_si / 1e6, reactant.conformer.E0.value_si / 1e6, 2)
            self.assertEqual(reactant0.conformer.E0.units, reactant.conformer.E0.units)
        for product0, product in zip(self.reaction.products, reaction.products):
            self.assertAlmostEqual(product0.conformer.E0.value_si / 1e6, product.conformer.E0.value_si / 1e6, 2)
            self.assertEqual(product0.conformer.E0.units, product.conformer.E0.units)
        self.assertAlmostEqual(self.reaction.transitionState.conformer.E0.value_si / 1e6, reaction.transitionState.conformer.E0.value_si / 1e6, 2)
        self.assertEqual(self.reaction.transitionState.conformer.E0.units, reaction.transitionState.conformer.E0.units)
        self.assertAlmostEqual(self.reaction.transitionState.frequency.value_si, reaction.transitionState.frequency.value_si, 2)
        self.assertEqual(self.reaction.transitionState.frequency.units, reaction.transitionState.frequency.units)
        
        self.assertAlmostEqual(self.reaction.kinetics.A.value_si, reaction.kinetics.A.value_si, delta=1e-6)
        self.assertAlmostEqual(self.reaction.kinetics.n.value_si, reaction.kinetics.n.value_si, delta=1e-6)
        self.assertAlmostEqual(self.reaction.kinetics.T0.value_si, reaction.kinetics.T0.value_si, delta=1e-6)
        self.assertAlmostEqual(self.reaction.kinetics.Ea.value_si, reaction.kinetics.Ea.value_si, delta=1e-6)
        self.assertEqual(self.reaction.kinetics.comment, reaction.kinetics.comment)
        
        self.assertEqual(self.reaction.duplicate, reaction.duplicate)
        self.assertEqual(self.reaction.degeneracy, reaction.degeneracy)   

################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
