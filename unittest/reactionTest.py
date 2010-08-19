#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
import unittest
import sys
sys.path.append('.')

from chempy.species import Species, TransitionState
from chempy.reaction import *
from chempy.states import *
from chempy.kinetics import ArrheniusModel
from chempy.thermo import WilhoitModel

################################################################################

class ReactionTest(unittest.TestCase):
    """
    Contains unit tests for the chempy.reaction module, used for working with
    chemical reaction objects.
    """
    
    def testReactionThermo(self):
        """
        Tests the reaction thermodynamics functions using the reaction
        acetyl + oxygen -> acetylperoxy.
        """
        
        # CC(=O)O[O]
        acetylperoxy = Species(
            label='acetylperoxy',
            thermo=WilhoitModel(cp0=4.0*constants.R, cpInf=21.0*constants.R, a0=-3.95, a1=9.26, a2=-15.6, a3=8.55, B=500.0, H0=-6.151e+04, S0=-790.2),
        )

        # C[C]=O
        acetyl = Species(
            label='acetyl',
            thermo=WilhoitModel(cp0=4.0*constants.R, cpInf=15.5*constants.R, a0=0.2541, a1=-0.4712, a2=-4.434, a3=2.25, B=500.0, H0=-1.439e+05, S0=-524.6),
        )

        # [O][O]
        oxygen = Species(
            label='oxygen',
            thermo=WilhoitModel(cp0=3.5*constants.R, cpInf=4.5*constants.R, a0=-0.9324, a1=26.18, a2=-70.47, a3=44.12, B=500.0, H0=1.453e+04, S0=-12.19),
        )
        
        reaction = Reaction(
            reactants=[acetyl, oxygen],
            products=[acetylperoxy],
            kinetics=ArrheniusModel(A=2.65e6, n=0.0, Ea=0.0*4184),
        )

        Tlist = numpy.arange(200.0, 2001.0, 200.0, numpy.float64)

        Hlist0 = [float(v) for v in ['-146007', '-145886', '-144195', '-141973', '-139633', '-137341', '-135155', '-133093', '-131150', '-129316']]
        Slist0 = [float(v) for v in ['-156.793', '-156.872', '-153.504', '-150.317', '-147.707', '-145.616', '-143.93', '-142.552', '-141.407', '-140.441']]
        Glist0 = [float(v) for v in ['-114648', '-83137.2', '-52092.4', '-21719.3', '8073.53', '37398.1', '66346.8', '94990.6', '123383', '151565']]
        Kalist0 = [float(v) for v in ['8.75951e+29', '7.1843e+10', '34272.7', '26.1877', '0.378696', '0.0235579', '0.00334673', '0.000792389', '0.000262777', '0.000110053']]
        Kclist0 = [float(v) for v in ['1.45661e+28', '2.38935e+09', '1709.76', '1.74189', '0.0314866', '0.00235045', '0.000389568', '0.000105413', '3.93273e-05', '1.83006e-05']]
        Kplist0 = [float(v) for v in ['8.75951e+24', '718430', '0.342727', '0.000261877', '3.78696e-06', '2.35579e-07', '3.34673e-08', '7.92389e-09', '2.62777e-09', '1.10053e-09']]

        Hlist = reaction.getEnthalpiesOfReaction(Tlist)
        Slist = reaction.getEntropiesOfReaction(Tlist)
        Glist = reaction.getFreeEnergiesOfReaction(Tlist)
        Kalist = reaction.getEquilibriumConstants(Tlist, type='Ka')
        Kclist = reaction.getEquilibriumConstants(Tlist, type='Kc')
        Kplist = reaction.getEquilibriumConstants(Tlist, type='Kp')

        for i in range(len(Tlist)):
            self.assertAlmostEqual( Hlist[i] /  Hlist0[i], 1.0, 4)
            self.assertAlmostEqual( Slist[i] /  Slist0[i], 1.0, 4)
            self.assertAlmostEqual( Glist[i] /  Glist0[i], 1.0, 4)
            self.assertAlmostEqual(Kalist[i] / Kalist0[i], 1.0, 4)
            self.assertAlmostEqual(Kclist[i] / Kclist0[i], 1.0, 4)
            self.assertAlmostEqual(Kplist[i] / Kplist0[i], 1.0, 4)

    def testTSTCalculation(self):
        """
        A test of the transition state theory k(T) calculation function,
        using the reaction H + C2H4 -> C2H5.
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
        
        reaction = Reaction(reactants=[hydrogen, ethylene], products=[ethyl])
    
        import numpy
        Tlist = 1000.0/numpy.arange(0.4, 3.35, 0.05)
        klist = reaction.calculateTSTRateCoefficient(Tlist, TS, tunneling='')
        arrhenius = ArrheniusModel().fitToData(Tlist, klist)
        klist2 = arrhenius.getRateCoefficients(Tlist)

        # Check that the correct Arrhenius parameters are returned
        self.assertAlmostEqual(arrhenius.A/458.87, 1.0, 2)
        self.assertAlmostEqual(arrhenius.n/0.978, 1.0, 2)
        self.assertAlmostEqual(arrhenius.Ea/10194, 1.0, 2)
        # Check that the fit is satisfactory
        for i in range(len(Tlist)):
            self.assertTrue(abs(1 - klist2[i] / klist[i]) < 0.01)
        
        
if __name__ == '__main__':
    unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )
