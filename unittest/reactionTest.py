#!/usr/bin/env python
# encoding: utf-8 -*-

"""
This module contains unit tests of the rmgpy.reaction module.
"""

import numpy
import unittest

from rmgpy.species import Species, TransitionState
from rmgpy.reaction import *
from rmgpy.statmech import *
from rmgpy.kinetics import Arrhenius
from rmgpy.thermo import Wilhoit
from rmgpy.quantity import constants

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
        states = StatesModel(
            modes = [Translation(mass=0.0280313), RigidRotor(linear=False, inertia=[5.69516e-47, 2.77584e-46, 3.34536e-46], symmetry=4), HarmonicOscillator(frequencies=[834.499, 973.312, 975.369, 1067.13, 1238.46, 1379.46, 1472.29, 1691.34, 3121.57, 3136.7, 3192.46, 3220.98])],
            spinMultiplicity=1,
        )
        ethylene = Species(states=states, E0=-205882860.949)
        
        states = StatesModel(
            modes = [Translation(mass=0.00100783), HarmonicOscillator(frequencies=[])],
            spinMultiplicity=2,
        )
        hydrogen = Species(states=states, E0=-1318675.56138)
        
        states = StatesModel(
            modes = [Translation(mass=0.0290391), RigidRotor(linear=False, inertia=[8.07491e-47, 3.69475e-46, 3.9885e-46], symmetry=1), HarmonicOscillator(frequencies=[466.816, 815.399, 974.674, 1061.98, 1190.71, 1402.03, 1467, 1472.46, 1490.98, 2972.34, 2994.88, 3089.96, 3141.01, 3241.96])],
            spinMultiplicity=2,
        )
        ethyl = Species(states=states, E0=-207340036.867)
        
        states = StatesModel(
            modes = [Translation(mass=0.0290391), RigidRotor(linear=False, inertia=[1.2553e-46, 3.68827e-46, 3.80416e-46], symmetry=2), HarmonicOscillator(frequencies=[241.47, 272.706, 833.984, 961.614, 974.994, 1052.32, 1238.23, 1364.42, 1471.38, 1655.51, 3128.29, 3140.3, 3201.94, 3229.51])],
            spinMultiplicity=2,
        )
        TS = TransitionState(states=states, E0=-207188826.467, frequency=-309.3437)
        
        kinetics = Arrhenius(A=1.0e6, n=1.0, Ea=10000.0, T0=298.15, comment='These parameters are completely made up')

        self.reaction = Reaction(reactants=[hydrogen, ethylene], products=[ethyl], kinetics=kinetics, transitionState=TS)
    
        # CC(=O)O[O]
        acetylperoxy = Species(
            label='acetylperoxy',
            thermo=Wilhoit(cp0=4.0*constants.R, cpInf=21.0*constants.R, a0=-3.95, a1=9.26, a2=-15.6, a3=8.55, B=500.0, H0=-6.151e+04, S0=-790.2),
        )

        # C[C]=O
        acetyl = Species(
            label='acetyl',
            thermo=Wilhoit(cp0=4.0*constants.R, cpInf=15.5*constants.R, a0=0.2541, a1=-0.4712, a2=-4.434, a3=2.25, B=500.0, H0=-1.439e+05, S0=-524.6),
        )

        # [O][O]
        oxygen = Species(
            label='oxygen',
            thermo=Wilhoit(cp0=3.5*constants.R, cpInf=4.5*constants.R, a0=-0.9324, a1=26.18, a2=-70.47, a3=44.12, B=500.0, H0=1.453e+04, S0=-12.19),
        )
        
        self.reaction2 = Reaction(reactants=[acetyl, oxygen], products=[acetylperoxy], kinetics=Arrhenius(A=2.65e6, n=0.0, Ea=0.0*4184))
        
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

    def testTSTCalculation(self):
        """
        A test of the transition state theory k(T) calculation function,
        using the reaction H + C2H4 -> C2H5.
        """
        Tlist = 1000.0/numpy.arange(0.4, 3.35, 0.05)
        klist = self.reaction.calculateTSTRateCoefficients(Tlist, tunneling='')
        arrhenius = Arrhenius().fitToData(Tlist, klist, kunits='')
        klist2 = arrhenius.getRateCoefficients(Tlist)
        
        # Check that the correct Arrhenius parameters are returned
        self.assertAlmostEqual(arrhenius.A.value/1.07506e+07, 1.0, 3)
        self.assertAlmostEqual(arrhenius.n.value/1.47803, 1.0, 3)
        self.assertAlmostEqual(arrhenius.Ea.value/10194., 1.0, 3)
        # Check that the fit is satisfactory
        for i in range(len(Tlist)):
            self.assertTrue(abs(1 - klist2[i] / klist[i]) < 0.01)
        
    def testPickle(self):
        """
        Test that a Reaction object can be successfully pickled and
        unpickled with no loss of information.
        """
        import cPickle
        reaction = cPickle.loads(cPickle.dumps(self.reaction))

        self.assertEqual(len(self.reaction.reactants), len(reaction.reactants))
        self.assertEqual(len(self.reaction.products), len(reaction.products))
        for reactant0, reactant in zip(self.reaction.reactants, reaction.reactants):
            self.assertAlmostEqual(reactant0.E0.value / 1e6, reactant.E0.value / 1e6, 2)
            self.assertEqual(reactant0.E0.units, reactant.E0.units)
        for product0, product in zip(self.reaction.products, reaction.products):
            self.assertAlmostEqual(product0.E0.value / 1e6, product.E0.value / 1e6, 2)
            self.assertEqual(product0.E0.units, product.E0.units)
        self.assertAlmostEqual(self.reaction.transitionState.E0.value / 1e6, reaction.transitionState.E0.value / 1e6, 2)
        self.assertEqual(self.reaction.transitionState.E0.units, reaction.transitionState.E0.units)
        self.assertAlmostEqual(self.reaction.transitionState.frequency.value, reaction.transitionState.frequency.value, 2)
        self.assertEqual(self.reaction.transitionState.frequency.units, reaction.transitionState.frequency.units)
        
        self.assertEqual(self.reaction.kinetics.A.value, reaction.kinetics.A.value)
        self.assertEqual(self.reaction.kinetics.n.value, reaction.kinetics.n.value)
        self.assertEqual(self.reaction.kinetics.T0.value, reaction.kinetics.T0.value)
        self.assertEqual(self.reaction.kinetics.Ea.value, reaction.kinetics.Ea.value)
        self.assertEqual(self.reaction.kinetics.comment, reaction.kinetics.comment)

        self.assertEqual(self.reaction.thirdBody, reaction.thirdBody)
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
            self.assertAlmostEqual(reactant0.E0.value / 1e6, reactant.E0.value / 1e6, 2)
            self.assertEqual(reactant0.E0.units, reactant.E0.units)
        for product0, product in zip(self.reaction.products, reaction.products):
            self.assertAlmostEqual(product0.E0.value / 1e6, product.E0.value / 1e6, 2)
            self.assertEqual(product0.E0.units, product.E0.units)
        self.assertAlmostEqual(self.reaction.transitionState.E0.value / 1e6, reaction.transitionState.E0.value / 1e6, 2)
        self.assertEqual(self.reaction.transitionState.E0.units, reaction.transitionState.E0.units)
        self.assertAlmostEqual(self.reaction.transitionState.frequency.value, reaction.transitionState.frequency.value, 2)
        self.assertEqual(self.reaction.transitionState.frequency.units, reaction.transitionState.frequency.units)
        
        self.assertEqual(self.reaction.kinetics.A.value, reaction.kinetics.A.value)
        self.assertEqual(self.reaction.kinetics.n.value, reaction.kinetics.n.value)
        self.assertEqual(self.reaction.kinetics.T0.value, reaction.kinetics.T0.value)
        self.assertEqual(self.reaction.kinetics.Ea.value, reaction.kinetics.Ea.value)
        self.assertEqual(self.reaction.kinetics.comment, reaction.kinetics.comment)
        
        self.assertEqual(self.reaction.thirdBody, reaction.thirdBody)
        self.assertEqual(self.reaction.duplicate, reaction.duplicate)
        self.assertEqual(self.reaction.degeneracy, reaction.degeneracy)   

################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
