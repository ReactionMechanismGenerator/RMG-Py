#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2018 Prof. William H. Green (whgreen@mit.edu),           #
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

import numpy
import unittest
from external.wip import work_in_progress

from rmgpy.species import Species, TransitionState
from rmgpy.reaction import Reaction
from rmgpy.quantity import Quantity
from rmgpy.statmech.translation import Translation, IdealGasTranslation
from rmgpy.statmech.rotation import Rotation, LinearRotor, NonlinearRotor, KRotor, SphericalTopRotor
from rmgpy.statmech.vibration import Vibration, HarmonicOscillator
from rmgpy.statmech.torsion import Torsion, HinderedRotor
from rmgpy.statmech.conformer import Conformer
from rmgpy.kinetics import Arrhenius, ArrheniusEP
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
    def test2to3_usingCheckOnlyLabel(self):
        r1 = self.makeReaction('AB=CDE')
        self.assertTrue(r1.isIsomorphic(self.makeReaction('AB=CDE'),checkOnlyLabel=True))
        self.assertTrue(r1.isIsomorphic(self.makeReaction('BA=EDC'),eitherDirection=False,checkOnlyLabel=True))
        self.assertFalse(r1.isIsomorphic(self.makeReaction('Ab=CDE'),checkOnlyLabel=True))
        self.assertFalse(r1.isIsomorphic(self.makeReaction('BA=EDd'),eitherDirection=False,checkOnlyLabel=True))


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
            degeneracy = 2,
        )
        self.reaction.kinetics.comment = '''
        Multiplied by reaction path degeneracy 2.0
        '''

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

        oxygen_atom = Species().fromSMILES('[O]')
        SO2 = Species().fromSMILES('O=S=O')
        SO3 = Species().fromSMILES('O=S(=O)=O')

        self.reaction3 = Reaction(
            reactants=[oxygen_atom, SO2],
            products=[SO3],
            kinetics = Arrhenius(A=(3.7e+11, 'cm^3/(mol*s)'), n=0, Ea=(1689, 'cal/mol'), T0=(1, 'K')))
        
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
    
    def testFixBarrierHeight(self):
        """
        Test that fixBarrierHeight:
            1) raises Ea to match endothermicity of reaction
            2) forces Ea to be positive if forcePositive=True
            3) Evans-Polanyi kinetics are handled so that negative Ea if Ea<E0 are set to min(0,E0)
        """
        
        #setup
        rxn = self.reaction2.copy()
        revRxn = rxn.copy()
        revRxn.reactants = rxn.products
        revRxn.products = rxn.reactants
        
        #test that endothermicity is matched 
        rxn.fixBarrierHeight()
        Ea = rxn.kinetics.Ea.value_si
        self.assertTrue(Ea==0.0)
        
        revRxn.fixBarrierHeight()
        Ea = revRxn.kinetics.Ea.value_si
        H0 = sum([spec.getThermoData().E0.value_si for spec in rxn.products]) \
            - sum([spec.getThermoData().E0.value_si for spec in rxn.reactants])
        self.assertAlmostEqual(Ea,-H0,3)
        
        #test that Ea is forced to be positive if forcePositive is set to True
        Ea = Quantity((-10000.0,'J/mol'))
        rxn.kinetics.Ea = Ea
        rxn.fixBarrierHeight()
        self.assertTrue(rxn.kinetics.Ea.value_si==Ea.value_si)
        
        rxn.fixBarrierHeight(forcePositive=True)
        self.assertTrue(rxn.kinetics.Ea.value_si==0.0)
        
        #Test for ArrheniusEP handling
        #if calculated Ea < 0 and Ea < E0, Ea is set to min(0,E0)
        H298 = rxn.getEnthalpyOfReaction(298)
        E0s = [-1000000.0,-10.0,0.0,10.0,1000000.0]
        
        for i,E0 in enumerate(E0s):
            kinetics = ArrheniusEP(
                A = (1.0, rxn.kinetics.A.units),
                n = (0, rxn.kinetics.n.units),
                alpha = 1.0,
                E0 = (E0, 'J/mol'),
            )
            rxn.kinetics = kinetics
            rxn.fixBarrierHeight()
            Ea = rxn.kinetics.Ea.value_si
            if i < 2:
                self.assertTrue(Ea==E0)
            elif i < 4:
                self.assertTrue(Ea==0.0)
            else:
                self.assertTrue(Ea==E0+H298)
            
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

    def testGenerateReverseRateCoefficientMultiPDepArrhenius(self):
        """
        Test the Reaction.generateReverseRateCoefficient() method works for the MultiPDepArrhenius format.
        """
        from rmgpy.kinetics import PDepArrhenius, MultiPDepArrhenius

        Tmin = 350.
        Tmax = 1500.
        Pmin = 1e-1
        Pmax = 1e1
        pressures = numpy.array([1e-1,1e1])
        comment = 'CH3 + C2H6 <=> CH4 + C2H5 (Baulch 2005)'
        arrhenius = [
            PDepArrhenius(
                pressures = (pressures,"bar"),
                arrhenius = [
                    Arrhenius(
                        A = (9.3e-16,"cm^3/(molecule*s)"),
                        n = 0.0,
                        Ea = (4740*constants.R*0.001,"kJ/mol"),
                        T0 = (1,"K"),
                        Tmin = (Tmin,"K"),
                        Tmax = (Tmax,"K"),
                        comment = comment,
                    ),
                    Arrhenius(
                        A = (9.3e-14,"cm^3/(molecule*s)"),
                        n = 0.0,
                        Ea = (4740*constants.R*0.001,"kJ/mol"),
                        T0 = (1,"K"),
                        Tmin = (Tmin,"K"),
                        Tmax = (Tmax,"K"),
                        comment = comment,
                    ),
                ],
                Tmin = (Tmin,"K"), 
                Tmax = (Tmax,"K"), 
                Pmin = (Pmin,"bar"), 
                Pmax = (Pmax,"bar"),
                comment = comment,
            ),
            PDepArrhenius(
                pressures = (pressures,"bar"),
                arrhenius = [
                    Arrhenius(
                        A = (1.4e-11,"cm^3/(molecule*s)"),
                        n = 0.0,
                        Ea = (11200*constants.R*0.001,"kJ/mol"),
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
                ],
                Tmin = (Tmin,"K"), 
                Tmax = (Tmax,"K"), 
                Pmin = (Pmin,"bar"), 
                Pmax = (Pmax,"bar"),
                comment = comment,
            ),
        ]  

        original_kinetics = MultiPDepArrhenius(
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

    def testGenerateReverseRateCoefficientThirdBody(self):
        """
        Test the Reaction.generateReverseRateCoefficient() method works for the ThirdBody format.
        """

        from rmgpy.kinetics import ThirdBody

        arrheniusLow = Arrhenius(
            A = (2.62e+33,"cm^6/(mol^2*s)"), 
            n = -4.76, 
            Ea = (10.21,"kJ/mol"), 
            T0 = (1,"K"),
        )
        efficiencies = {"C": 3, "C(=O)=O": 2, "CC": 3, "O": 6, "[Ar]": 0.7, "[C]=O": 1.5, "[H][H]": 2}
        Tmin = 300.
        Tmax = 2000.
        Pmin = 0.01
        Pmax = 100.
        comment = """H + CH3 -> CH4"""
        thirdBody = ThirdBody(
            arrheniusLow = arrheniusLow,
            Tmin = (Tmin,"K"),
            Tmax = (Tmax,"K"),
            Pmin = (Pmin,"bar"),
            Pmax = (Pmax,"bar"),
            efficiencies = efficiencies,
            comment = comment,
        )
         
        original_kinetics = thirdBody

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

    def testGenerateReverseRateCoefficientLindemann(self):
        """
        Test the Reaction.generateReverseRateCoefficient() method works for the Lindemann format.
        """

        from rmgpy.kinetics import Lindemann

        arrheniusHigh = Arrhenius(
            A = (1.39e+16,"cm^3/(mol*s)"), 
            n = -0.534, 
            Ea = (2.243,"kJ/mol"), 
            T0 = (1,"K"),
        )
        arrheniusLow = Arrhenius(
            A = (2.62e+33,"cm^6/(mol^2*s)"), 
            n = -4.76, 
            Ea = (10.21,"kJ/mol"), 
            T0 = (1,"K"),
        )
        efficiencies = {"C": 3, "C(=O)=O": 2, "CC": 3, "O": 6, "[Ar]": 0.7, "[C]=O": 1.5, "[H][H]": 2}
        Tmin = 300.
        Tmax = 2000.
        Pmin = 0.01
        Pmax = 100.
        comment = """H + CH3 -> CH4"""
        lindemann = Lindemann(
            arrheniusHigh = arrheniusHigh,
            arrheniusLow = arrheniusLow,
            Tmin = (Tmin,"K"),
            Tmax = (Tmax,"K"),
            Pmin = (Pmin,"bar"),
            Pmax = (Pmax,"bar"),
            efficiencies = efficiencies,
            comment = comment,
        )
         
        original_kinetics = lindemann
        
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


    def testGenerateReverseRateCoefficientTroe(self):
        """
        Test the Reaction.generateReverseRateCoefficient() method works for the Troe format.
        """

        from rmgpy.kinetics import Troe

        arrheniusHigh = Arrhenius(
            A = (1.39e+16,"cm^3/(mol*s)"), 
            n = -0.534, 
            Ea = (2.243,"kJ/mol"), 
            T0 = (1,"K"),
        )
        arrheniusLow = Arrhenius(
            A = (2.62e+33,"cm^6/(mol^2*s)"), 
            n = -4.76, 
            Ea = (10.21,"kJ/mol"), 
            T0 = (1,"K"),
        )
        alpha = 0.783
        T3 = 74
        T1 = 2941
        T2 = 6964
        efficiencies = {"C": 3, "C(=O)=O": 2, "CC": 3, "O": 6, "[Ar]": 0.7, "[C]=O": 1.5, "[H][H]": 2}
        Tmin = 300.
        Tmax = 2000.
        Pmin = 0.01
        Pmax = 100.
        comment = """H + CH3 -> CH4"""
        troe = Troe(
            arrheniusHigh = arrheniusHigh,
            arrheniusLow = arrheniusLow,
            alpha = alpha,
            T3 = (T3,"K"),
            T1 = (T1,"K"),
            T2 = (T2,"K"),
            Tmin = (Tmin,"K"),
            Tmax = (Tmax,"K"),
            Pmin = (Pmin,"bar"),
            Pmax = (Pmax,"bar"),
            efficiencies = efficiencies,
            comment = comment,
        )
         
        original_kinetics = troe
        
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

    def testDegeneracyUpdatesRate(self):
        """
        This method tests that a change in degeneracy will result in a modified rate constant
        """

        prefactor = self.reaction.kinetics.A.value_si
        degeneracyFactor = 2
        self.reaction.degeneracy *= degeneracyFactor
        self.assertAlmostEqual(self.reaction.kinetics.A.value_si, degeneracyFactor * prefactor)

    def testDegeneracyUpdatesKineticsComment(self):
        """
        This method tests that a change in degeneracy will result in a modified rate constant
        """

        newDegeneracy = 8
        self.reaction.degeneracy = newDegeneracy
        self.assertIn('Multiplied by reaction path degeneracy 8.0', self.reaction.kinetics.comment)

    def testSulfurReactionPairs(self):
        """
        This method tests that reaction pairs are being generated for sulfur species
        """

        self.reaction3.generatePairs()
        self.assertEqual(len(self.reaction3.pairs[0]), 2)
        self.assertEqual(len(self.reaction3.pairs[1]), 2)


class TestReactionToCantera(unittest.TestCase):
    """
    Contains unit tests of the Reaction class associated with forming Cantera objects.
    """

    def setUp(self):
        """
        A method that is called prior to each unit test in this class.
        """
        from rmgpy.kinetics import Arrhenius, MultiArrhenius, PDepArrhenius, MultiPDepArrhenius, ThirdBody, Troe, Lindemann, Chebyshev
        from rmgpy.molecule import Molecule
        from rmgpy.thermo import NASA, NASAPolynomial
        import cantera as ct
        
        # define some species:
        ch3 = Species(index=13, label="CH3", thermo=NASA(polynomials=[NASAPolynomial(coeffs=[3.91547,0.00184154,3.48744e-06,-3.3275e-09,8.49964e-13,16285.6,0.351739], Tmin=(100,'K'), Tmax=(1337.62,'K')), NASAPolynomial(coeffs=[3.54145,0.00476788,-1.82149e-06,3.28878e-10,-2.22547e-14,16224,1.6604], Tmin=(1337.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), comment="""
Thermo library: primaryThermoLibrary + radical(CH3)
"""), molecule=[Molecule(SMILES="[CH3]")])

        ethane = Species(label="ethane", thermo=NASA(polynomials=[NASAPolynomial(coeffs=[3.78033,-0.00324263,5.52381e-05,-6.38581e-08,2.28637e-11,-11620.3,5.21034], Tmin=(100,'K'), Tmax=(954.51,'K')), NASAPolynomial(coeffs=[4.58983,0.0141508,-4.75962e-06,8.60294e-10,-6.21717e-14,-12721.8,-3.61739], Tmin=(954.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), comment="""
Thermo group additivity estimation: group(Cs-CsHHH) + gauche(Cs(CsRRR)) + other(R) + group(Cs-CsHHH) + gauche(Cs(CsRRR)) + other(R)
"""), molecule=[Molecule(SMILES="CC")])
                      
        co2 = Species(index=16, label="CO2", thermo=NASA(polynomials=[NASAPolynomial(coeffs=[3.27861,0.00274152,7.16074e-06,-1.08027e-08,4.14282e-12,-48470.3,5.97937], Tmin=(100,'K'), Tmax=(988.89,'K')), NASAPolynomial(coeffs=[4.5461,0.00291913,-1.15484e-06,2.27654e-10,-1.7091e-14,-48980.4,-1.43275], Tmin=(988.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), comment="""
Thermo group additivity estimation: group(Cdd-OdOd) + other(R) + group(O2d-Cd) + other(R) + group(O2d-Cd) + other(R)
"""), molecule=[Molecule(SMILES="O=C=O")])

        ch4 = Species(index=15, label="CH4", thermo=NASA(polynomials=[NASAPolynomial(coeffs=[4.20541,-0.00535556,2.51123e-05,-2.13762e-08,5.97522e-12,-10161.9,-0.921275], Tmin=(100,'K'), Tmax=(1084.12,'K')), NASAPolynomial(coeffs=[0.908272,0.0114541,-4.57173e-06,8.2919e-10,-5.66314e-14,-9719.98,13.9931], Tmin=(1084.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), comment="""
Thermo library: primaryThermoLibrary
"""), molecule=[Molecule(SMILES="C")])
        
        h2o = Species(index=27, label="H2O", thermo=NASA(polynomials=[NASAPolynomial(coeffs=[4.05764,-0.000787933,2.90876e-06,-1.47518e-09,2.12838e-13,-30281.6,-0.311363], Tmin=(100,'K'), Tmax=(1130.24,'K')), NASAPolynomial(coeffs=[2.84325,0.00275108,-7.8103e-07,1.07243e-10,-5.79389e-15,-29958.6,5.91041], Tmin=(1130.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), comment="""
Thermo library: primaryThermoLibrary
"""), molecule=[Molecule(SMILES="O")])
        
        ar = Species(label="Ar", thermo=NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,4.37967], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,4.37967], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), comment="""
Thermo library: primaryThermoLibrary
"""), molecule=[Molecule(SMILES="[Ar]")])
        
        h2 = Species(index=2, label="H2", thermo=NASA(polynomials=[NASAPolynomial(coeffs=[3.43536,0.00021271,-2.78625e-07,3.40267e-10,-7.76031e-14,-1031.36,-3.90842], Tmin=(100,'K'), Tmax=(1959.08,'K')), NASAPolynomial(coeffs=[2.78816,0.000587644,1.59009e-07,-5.52736e-11,4.34309e-15,-596.143,0.112747], Tmin=(1959.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), comment="""
Thermo library: primaryThermoLibrary
"""), molecule=[Molecule(SMILES="[H][H]")])
        
        h = Species(index=3, label="H", thermo=NASA(polynomials=[NASAPolynomial(coeffs=[2.5,-1.91243e-12,2.45329e-15,-1.02377e-18,1.31369e-22,25474.2,-0.444973], Tmin=(100,'K'), Tmax=(4563.27,'K')), NASAPolynomial(coeffs=[2.50167,-1.43051e-06,4.6025e-10,-6.57826e-14,3.52412e-18,25472.7,-0.455578], Tmin=(4563.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), comment="""
Thermo library: primaryThermoLibrary
"""), molecule=[Molecule(SMILES="[H]")])
        
        oh = Species(index=4, label="OH", thermo=NASA(polynomials=[NASAPolynomial(coeffs=[3.51457,2.92773e-05,-5.32163e-07,1.01949e-09,-3.85945e-13,3414.25,2.10435], Tmin=(100,'K'), Tmax=(1145.75,'K')), NASAPolynomial(coeffs=[3.07194,0.000604016,-1.39783e-08,-2.13446e-11,2.48066e-15,3579.39,4.578], Tmin=(1145.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), comment="""
Thermo library: primaryThermoLibrary
"""), molecule=[Molecule(SMILES="[OH]")])
        
        ho2 = Species(index=5, label="HO2", thermo=NASA(polynomials=[NASAPolynomial(coeffs=[4.04594,-0.00173464,1.03766e-05,-1.02202e-08,3.34908e-12,-986.754,4.63581], Tmin=(100,'K'), Tmax=(932.15,'K')), NASAPolynomial(coeffs=[3.21024,0.00367942,-1.27701e-06,2.18045e-10,-1.46338e-14,-910.369,8.18291], Tmin=(932.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), comment="""
Thermo group additivity estimation: group(O2s-OsH) + gauche(O2s(RR)) + other(R) + group(O2s-OsH) + gauche(O2s(RR)) + other(R) + radical(HOOJ)
"""), molecule=[Molecule(SMILES="[O]O")])
        
        o2 = Species(index=6, label="O2", thermo=NASA(polynomials=[NASAPolynomial(coeffs=[3.53732,-0.00121572,5.3162e-06,-4.89446e-09,1.45846e-12,-1038.59,4.68368], Tmin=(100,'K'), Tmax=(1074.55,'K')), NASAPolynomial(coeffs=[3.15382,0.00167804,-7.69974e-07,1.51275e-10,-1.08782e-14,-1040.82,6.16756], Tmin=(1074.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), comment="""
Thermo library: primaryThermoLibrary
"""), molecule=[Molecule(SMILES="[O][O]")])
        
        co = Species(index=9, label="CO", thermo=NASA(polynomials=[NASAPolynomial(coeffs=[3.66965,-0.00550953,2.00538e-05,-2.08391e-08,7.43738e-12,1200.77,-12.4224], Tmin=(100,'K'), Tmax=(884.77,'K')), NASAPolynomial(coeffs=[2.8813,0.00231665,-4.40151e-07,4.75633e-11,-2.78282e-15,1173.45,-9.65831], Tmin=(884.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), comment="""
Thermo group additivity estimation: group(Ct-CtCs) + other(R) + group(O2s-CsCs) + other(R)
"""), molecule=[Molecule(SMILES="[C-]#[O+]")])
        
        h2o2 = Species(index=7, label="H2O2", thermo=NASA(polynomials=[NASAPolynomial(coeffs=[3.73136,0.00335071,9.35033e-06,-1.521e-08,6.41585e-12,-17721.2,5.45911], Tmin=(100,'K'), Tmax=(908.87,'K')), NASAPolynomial(coeffs=[5.41579,0.00261008,-4.39892e-07,4.91087e-11,-3.35188e-15,-18303,-4.02248], Tmin=(908.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), comment="""
Thermo group additivity estimation: group(O2s-OsH) + gauche(O2s(RR)) + other(R) + group(O2s-OsH) + gauche(O2s(RR)) + other(R)
"""), molecule=[Molecule(SMILES="OO")])
        
        self.speciesList = [ch3,ethane,co2,ch4,h2o,ar,h2, h, oh, ho2, o2, co, h2o2]
                
        self.troe = Reaction(index=1, reactants=[ch3,ch3], products=[ethane], 
                        kinetics=Troe(arrheniusHigh=Arrhenius(A=(6.77e+16,'cm^3/(mol*s)'), n=-1.18, Ea=(0.654,'kcal/mol'), T0=(1,'K')), arrheniusLow=Arrhenius(A=(3.4e+41,'cm^6/(mol^2*s)'), n=-7.03, Ea=(2.762,'kcal/mol'), T0=(1,'K')), alpha=0.619, T3=(73.2,'K'), T1=(1180,'K'), T2=(10000,'K'), 
                                      efficiencies={Molecule(SMILES="O=C=O"): 2.0, Molecule(SMILES="[H][H]"): 2.0, Molecule(SMILES="O"): 6.0, Molecule(SMILES="[Ar]"): 0.7, Molecule(SMILES="C"): 2.0, Molecule(SMILES="CC"): 3.0}))
        
        self.ct_troe = ct.Reaction.fromCti('''falloff_reaction('CH3(13) + CH3(13) (+ M) <=> ethane (+ M)',
                 kf=[(6.770000e+16,'cm3/mol/s'), -1.18, (0.654,'kcal/mol')],
                 kf0=[(3.400000e+41,'cm6/mol2/s'), -7.03, (2.762,'kcal/mol')],
                 efficiencies='ethane:3.0 CO2(16):2.0 CH4(15):2.0 Ar:0.7 H2O(27):6.0 H2(2):2.0',
                 falloff=Troe(A=0.619, T3=73.2, T1=1180.0, T2=10000.0))''')
        
        self.arrheniusBi = Reaction(index=2, reactants=[h,ch4], products=[h2,ch3],kinetics=Arrhenius(A=(6.6e+08,'cm^3/(mol*s)'), n=1.62, Ea=(10.84,'kcal/mol'), T0=(1,'K')))
        
        self.ct_arrheniusBi = ct.Reaction.fromCti('''reaction('H(3) + CH4(15) <=> H2(2) + CH3(13)', [(6.600000e+08,'cm3/mol/s'), 1.62, (10.84,'kcal/mol')])''')
        
        self.arrheniusBi_irreversible = Reaction(index=10, reactants=[h,ch4], products=[h2,ch3],kinetics=Arrhenius(A=(6.6e+08,'cm^3/(mol*s)'), n=1.62, Ea=(10.84,'kcal/mol'), T0=(1,'K')),reversible=False)
        
        self.ct_arrheniusBi_irreversible = ct.Reaction.fromCti('''reaction('H(3) + CH4(15) => H2(2) + CH3(13)', [(6.600000e+08,'cm3/mol/s'), 1.62, (10.84,'kcal/mol')])''')
        
        self.arrheniusMono = Reaction(index=15, reactants=[h2o2], products=[h2,o2],kinetics=Arrhenius(A=(6.6e+03,'1/s'), n=1.62, Ea=(10.84,'kcal/mol'), T0=(1,'K')))
        
        self.ct_arrheniusMono = ct.Reaction.fromCti('''reaction('H2O2(7) <=> H2(2) + O2(6)', [(6.600000e+03,'1/s'), 1.62, (10.84,'kcal/mol')])''')
        
        self.arrheniusTri = Reaction(index=20, reactants=[h,h,o2], products=[h2o2],kinetics=Arrhenius(A=(6.6e+08,'cm^6/(mol^2*s)'), n=1.62, Ea=(10.84,'kcal/mol'), T0=(1,'K')))
        self.ct_arrheniusTri = ct.Reaction.fromCti('''reaction('H(3) + H(3) + O2(6) <=> H2O2(7)', [(6.6e+08, 'cm6/mol2/s'), 1.62, (10.84,'kcal/mol')])''')
        
        self.multiArrhenius = Reaction(index=3, reactants=[oh,ho2], products=[h2o,o2],
                                  kinetics=MultiArrhenius(arrhenius=[Arrhenius(A=(1.45e+13,'cm^3/(mol*s)'), n=0, Ea=(-0.5,'kcal/mol'), T0=(1,'K')), Arrhenius(A=(5e+15,'cm^3/(mol*s)'), n=0, Ea=(17.33,'kcal/mol'), T0=(1,'K'))]))
        
        self.ct_multiArrhenius = [ ct.Reaction.fromCti('''reaction('OH(4) + HO2(5) <=> H2O(27) + O2(6)', [(1.450000e+13,'cm3/mol/s'), 0.0, (-0.5,'kcal/mol')],
         options='duplicate')'''), ct.Reaction.fromCti('''reaction('OH(4) + HO2(5) <=> H2O(27) + O2(6)', [(5.000000e+15,'cm3/mol/s'), 0.0, (17.33,'kcal/mol')],
         options='duplicate')''')]
        
        self.pdepArrhenius = Reaction(index=4, reactants=[ho2,ho2], products=[o2,h2o2],
                                 kinetics = PDepArrhenius(
                                    pressures = ([0.1, 1, 10], 'atm'),
                                    arrhenius = [
                                        Arrhenius(
                                            A = (8.8e+16, 'cm^3/(mol*s)'),
                                            n = -1.05,
                                            Ea = (6461, 'cal/mol'),
                                            T0 = (1, 'K'),
                                        ),
                                        Arrhenius(
                                            A = (8e+21, 'cm^3/(mol*s)'),
                                            n = -2.39,
                                            Ea = (11180, 'cal/mol'),
                                            T0 = (1, 'K'),
                                        ),
                                        Arrhenius(
                                            A = (3.3e+24, 'cm^3/(mol*s)'),
                                            n = -3.04,
                                            Ea = (15610, 'cal/mol'),
                                            T0 = (1, 'K'),
                                        ),
                                    ],
                                ),
                                )
        
        self.ct_pdepArrhenius = ct.Reaction.fromCti('''pdep_arrhenius('HO2(5) + HO2(5) <=> O2(6) + H2O2(7)',
               [(0.1, 'atm'), (8.800000e+16, 'cm3/mol/s'), -1.05, (6.461,'kcal/mol')],
               [(1.0, 'atm'), (8.000000e+21,'cm3/mol/s'), -2.39, (11.18,'kcal/mol')],
               [(10.0, 'atm'), (3.300000e+24,'cm3/mol/s'), -3.04, (15.61,'kcal/mol')])''')
        
        self.multiPdepArrhenius = Reaction(index=5, reactants=[ho2,ch3], products=[o2,ch4],
                                      kinetics = MultiPDepArrhenius(
                                                arrhenius = [
                                                    PDepArrhenius(
                                                        pressures = ([0.001, 1, 3], 'atm'),
                                                        arrhenius = [
                                                            Arrhenius(A=(9.3e+10, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
                                                            Arrhenius(A=(8e+10, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
                                                            Arrhenius(A=(7e+10, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
                                                        ],
                                                    ),
                                                    PDepArrhenius(
                                                        pressures = ([0.001, 1, 3], 'atm'),
                                                        arrhenius = [
                                                            Arrhenius(A=(710000, 'cm^3/(mol*s)'), n=1.8, Ea=(1133, 'cal/mol'), T0=(1, 'K')),
                                                            Arrhenius(A=(880000, 'cm^3/(mol*s)'), n=1.77, Ea=(954, 'cal/mol'), T0=(1, 'K')),
                                                            Arrhenius(A=(290000, 'cm^3/(mol*s)'), n=1.9, Ea=(397, 'cal/mol'), T0=(1, 'K')),
                                                        ],
                                                    ),
                                                ],
                                            ),
                                        )
        
        self.ct_multiPdepArrhenius = [ct.Reaction.fromCti('''pdep_arrhenius('HO2(5) + CH3(13) <=> O2(6) + CH4(15)',
               [(0.001, 'atm'), (9.300000e+10, 'cm3/mol/s'), 0.0, (0.0,'kcal/mol')],
               [(1.0, 'atm'), (8.000000e+10, 'cm3/mol/s'), 0.0, (0.0,'kcal/mol')],
               [(3.0, 'atm'), (7.000000e+10, 'cm3/mol/s'), 0.0, (0.0,'kcal/mol')],
               options='duplicate')'''),
               ct.Reaction.fromCti('''pdep_arrhenius('HO2(5) + CH3(13) <=> O2(6) + CH4(15)',
               [(0.001, 'atm'), (7.100000e+05, 'cm3/mol/s'), 1.8, (1.133,'kcal/mol')],
               [(1.0, 'atm'), (8.800000e+05, 'cm3/mol/s'), 1.77, (0.954,'kcal/mol')],
               [(3.0, 'atm'), (2.900000e+05, 'cm3/mol/s'), 1.9, (0.397,'kcal/mol')],
               options='duplicate')''')]
        
        self.chebyshev = Reaction(index=6, reactants=[h,ch3], products=[ch4], kinetics=Chebyshev(coeffs=[[12.68,0.3961,-0.05481,-0.003606],[-0.7128,0.731,-0.0941,-0.008587],[-0.5806,0.57,-0.05539,-0.01115],[-0.4074,0.3653,-0.0118,-0.01171],[-0.2403,0.1779,0.01946,-0.008505],[-0.1133,0.0485,0.03121,-0.002955]], kunits='cm^3/(mol*s)', Tmin=(300,'K'), Tmax=(3000,'K'), Pmin=(0.001,'atm'), Pmax=(98.692,'atm')))
        
        self.ct_chebyshev = ct.Reaction.fromCti('''chebyshev_reaction('H(3) + CH3(13) (+ M) <=> CH4(15) (+ M)',
                   Tmin=300.0, Tmax=3000.0,
                   Pmin=(0.001, 'atm'), Pmax=(98.692, 'atm'),
                   coeffs=[[ 9.68000e+00,  3.96100e-01, -5.48100e-02, -3.60600e-03],
                           [-7.12800e-01,  7.31000e-01, -9.41000e-02, -8.58700e-03],
                           [-5.80600e-01,  5.70000e-01, -5.53900e-02, -1.11500e-02],
                           [-4.07400e-01,  3.65300e-01, -1.18000e-02, -1.17100e-02],
                           [-2.40300e-01,  1.77900e-01,  1.94600e-02, -8.50500e-03],
                           [-1.13300e-01,  4.85000e-02,  3.12100e-02, -2.95500e-03]])''')
        
        
        self.thirdBody = Reaction(index=7, reactants=[h,h], products=[h2],
                             kinetics=ThirdBody(arrheniusLow=Arrhenius(A=(1e+18,'cm^6/(mol^2*s)'), n=-1, Ea=(0,'kcal/mol'), T0=(1,'K')), 
                                                efficiencies={Molecule(SMILES="O=C=O"): 0.0, Molecule(SMILES="[H][H]"): 0.0, Molecule(SMILES="O"): 0.0, 
                                                              Molecule(SMILES="[Ar]"): 0.63, Molecule(SMILES="C"): 2.0, Molecule(SMILES="CC"): 3.0}))
        
        
        self.ct_thirdBody = ct.Reaction.fromCti('''three_body_reaction('H(3) + H(3) + M <=> H2(2) + M', [(1.000000e+18,'cm6/mol2/s'), -1.0, (0.0,'kcal/mol')],
                    efficiencies='CO2(16):0.0 CH4(15):2.0 ethane:3.0 H2O(27):0.0 H2(2):0.0 Ar:0.63')''')
        
        self.lindemann = Reaction(index=8, reactants=[h,o2], products=[ho2], 
                             kinetics=Lindemann(arrheniusHigh=Arrhenius(A=(1.8e+10,'cm^3/(mol*s)'), n=0, Ea=(2.385,'kcal/mol'), T0=(1,'K')), arrheniusLow=Arrhenius(A=(6.02e+14,'cm^6/(mol^2*s)'), n=0, Ea=(3,'kcal/mol'), T0=(1,'K')), efficiencies={Molecule(SMILES="O=C=O"): 3.5, Molecule(SMILES="[H][H]"): 2.0, Molecule(SMILES="O"): 6.0, Molecule(SMILES="[Ar]"): 0.5, Molecule(SMILES="C"): 2.0, Molecule(SMILES="CC"): 3.0, Molecule(SMILES="[O][O]"): 6.0}))
    
        self.ct_lindemann =ct.Reaction.fromCti('''falloff_reaction('H(3) + O2(6) (+ M) <=> HO2(5) (+ M)',
                 kf=[(1.800000e+10,'cm3/mol/s'), 0.0, (2.385,'kcal/mol')],
                 kf0=[(6.020000e+14,'cm6/mol2/s'), 0.0, (3.0,'kcal/mol')],
                 efficiencies='CO2(16):3.5 CH4(15):2.0 ethane:3.0 H2O(27):6.0 O2(6):6.0 H2(2):2.0 Ar:0.5')''')
        
    def testArrhenius(self):
        """
        Tests formation of cantera reactions with Arrhenius or kinetics.
        """
        
        rmgObjects = [self.arrheniusBi, self.arrheniusBi_irreversible, self.arrheniusMono, self.arrheniusTri]
        
        ctObjects = [self.ct_arrheniusBi, self.ct_arrheniusBi_irreversible, self.ct_arrheniusMono, self.ct_arrheniusTri]
        converted_ctObjects = [obj.toCantera(self.speciesList, useChemkinIdentifier = True) for obj in rmgObjects]
        
        for converted_obj, ct_obj in zip(converted_ctObjects, ctObjects):
            # Check that the reaction class is the same
            self.assertEqual(type(converted_obj), type(ct_obj))  
            # Check that the reaction string is the same
            self.assertEqual(repr(converted_obj), repr(ct_obj))
            # Check that the Arrhenius string is identical
            self.assertEqual(str(converted_obj.rate), str(ct_obj.rate))

    def testMultiArrhenius(self):
        """
        Tests formation of cantera reactions with MultiArrhenius kinetics.
        """
        rmgObjects = [self.multiArrhenius]
        ctObjects = [self.ct_multiArrhenius]
        converted_ctObjects = [obj.toCantera(self.speciesList, useChemkinIdentifier = True) for obj in rmgObjects]
                
        for converted_obj, ct_obj in zip(converted_ctObjects, ctObjects):
            # Check that the same number of reactions are produced
            self.assertEqual(len(converted_obj), len(ct_obj))  
            
            for converted_rxn, ct_rxn in zip(converted_obj, ct_obj):
                # Check that the reaction has the same type
                self.assertEqual(type(converted_rxn), type(ct_rxn))
                # Check that the reaction string is the same
                self.assertEqual(repr(converted_rxn), repr(ct_rxn))
                # Check that the Arrhenius rates are identical
                self.assertEqual(str(converted_rxn.rate), str(ct_rxn.rate))
        
    def testPDepArrhenius(self):
        """
        Tests formation of cantera reactions with PDepArrhenius kinetics.
        """
        rmgObjects = [self.pdepArrhenius]
        ctObjects = [self.ct_pdepArrhenius]
        converted_ctObjects = [obj.toCantera(self.speciesList, useChemkinIdentifier = True) for obj in rmgObjects]
        
        for converted_obj, ct_obj in zip(converted_ctObjects, ctObjects):
            # Check that the reaction class is the same
            self.assertEqual(type(converted_obj), type(ct_obj))  
            # Check that the reaction string is the same
            self.assertEqual(repr(converted_obj), repr(ct_obj))
            # Check that the Arrhenius rates are identical
            self.assertEqual(str(converted_obj.rates), str(ct_obj.rates))
    
    def testMultiPdepArrhenius(self):
        """
        Tests formation of cantera reactions with MultiPDepArrhenius kinetics.
        """
        
        rmgObjects = [self.multiPdepArrhenius]
        ctObjects = [self.ct_multiPdepArrhenius]
        converted_ctObjects = [obj.toCantera(self.speciesList, useChemkinIdentifier = True) for obj in rmgObjects]
                
        for converted_obj, ct_obj in zip(converted_ctObjects, ctObjects):
            # Check that the same number of reactions are produced
            self.assertEqual(len(converted_obj), len(ct_obj))  
            
            for converted_rxn, ct_rxn in zip(converted_obj, ct_obj):
                # Check that the reaction has the same type
                self.assertEqual(type(converted_rxn), type(ct_rxn))
                # Check that the reaction string is the same
                self.assertEqual(repr(converted_rxn), repr(ct_rxn))
                # Check that the Arrhenius rates are identical
                self.assertEqual(str(converted_rxn.rates), str(ct_rxn.rates))
        
     
    def testChebyshev(self):
        """
        Tests formation of cantera reactions with Chebyshev kinetics.
        """
        ct_chebyshev = self.chebyshev.toCantera(self.speciesList, useChemkinIdentifier = True)
        self.assertEqual(type(ct_chebyshev),type(self.ct_chebyshev))
        self.assertEqual(repr(ct_chebyshev),repr(self.ct_chebyshev))
        
        self.assertEqual(ct_chebyshev.Tmax, self.ct_chebyshev.Tmax)
        self.assertEqual(ct_chebyshev.Tmin, self.ct_chebyshev.Tmin)
        self.assertEqual(ct_chebyshev.Pmax, self.ct_chebyshev.Pmax)
        self.assertEqual(ct_chebyshev.Pmin, self.ct_chebyshev.Pmin)
        self.assertTrue((ct_chebyshev.coeffs == self.ct_chebyshev.coeffs).all())
        
        
    def testFalloff(self):
        """
        Tests formation of cantera reactions with Falloff kinetics.
        """
        ct_thirdBody = self.thirdBody.toCantera(self.speciesList, useChemkinIdentifier = True)
        self.assertEqual(type(ct_thirdBody),type(self.ct_thirdBody))
        self.assertEqual(repr(ct_thirdBody),repr(self.ct_thirdBody))
        self.assertEqual(str(ct_thirdBody.rate), str(self.ct_thirdBody.rate))
        self.assertEqual(ct_thirdBody.efficiencies, self.ct_thirdBody.efficiencies)
        
        ct_lindemann = self.lindemann.toCantera(self.speciesList, useChemkinIdentifier = True)
        self.assertEqual(type(ct_lindemann),type(self.ct_lindemann))
        self.assertEqual(repr(ct_lindemann), repr(self.ct_lindemann))
        self.assertEqual(ct_lindemann.efficiencies, self.ct_lindemann.efficiencies)
        self.assertEqual(str(ct_lindemann.low_rate), str(self.ct_lindemann.low_rate))
        self.assertEqual(str(ct_lindemann.high_rate), str(self.ct_lindemann.high_rate))
        self.assertEqual(str(ct_lindemann.falloff), str(self.ct_lindemann.falloff))
        
        
        ct_troe = self.troe.toCantera(self.speciesList, useChemkinIdentifier = True)
        self.assertEqual(type(ct_troe),type(self.ct_troe))
        self.assertEqual(repr(ct_troe), repr(self.ct_troe))
        self.assertEqual(ct_troe.efficiencies, self.ct_troe.efficiencies)
        
        self.assertEqual(str(ct_troe.low_rate), str(self.ct_troe.low_rate))
        self.assertEqual(str(ct_troe.high_rate), str(self.ct_troe.high_rate))
        self.assertEqual(str(ct_troe.falloff), str(self.ct_troe.falloff))
        
################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
