#!/usr/bin/env python
# encoding: utf-8

"""
This module contains unit tests of the rmgpy.quantity module.
"""

import unittest
import math
import quantities as pq

from rmgpy.quantity import *

################################################################################

class ConstantsTest(unittest.TestCase):
    """
    Contains unit tests for the rmgpy.chem.constants module.
    """
    
    def testQuantity(self):
        """
        Unit tests for the Quantity class.
        """

        q = Quantity(10.0)
        self.assertEqual(q.value, 10.0)
        self.assertEqual(q.uncertainty, 0.0)
        self.assertEqual(q.values, None)
        self.assertEqual(q.uncertainties, None)
        self.assertEqual(q.units, "")
        self.assertEqual(q.uncertaintyType, "")

        q = Quantity((10.0,"cm"))
        self.assertEqual(q.value, 0.1)
        self.assertEqual(q.uncertainty, 0.0)
        self.assertEqual(q.values, None)
        self.assertEqual(q.uncertainties, None)
        self.assertEqual(q.units, "cm")
        self.assertEqual(q.uncertaintyType, "")

        q = Quantity((10.0,"cm","+|-",2.0))
        self.assertEqual(q.value, 0.1)
        self.assertEqual(q.uncertainty, 0.02)
        self.assertEqual(q.values, None)
        self.assertEqual(q.uncertainties, None)
        self.assertEqual(q.units, "cm")
        self.assertEqual(q.uncertaintyType, "+|-")
        
        q = Quantity((10.0,"cm","*|/",2.0))
        self.assertEqual(q.value, 0.1)
        self.assertEqual(q.uncertainty, 2.0)
        self.assertEqual(q.values, None)
        self.assertEqual(q.uncertainties, None)
        self.assertEqual(q.units, "cm")
        self.assertEqual(q.uncertaintyType, "*|/")
        
        q = Quantity(([30.0,20.0,10.0],"cm","+|-",2.0))
        self.assertEqual(q.value, 0.0)
        self.assertEqual(q.uncertainty, 0.02)
        self.assertEqual(q.values[0], 0.3)
        self.assertEqual(q.values[1], 0.2)
        self.assertEqual(q.values[2], 0.1)
        self.assertEqual(q.uncertainties, None)
        self.assertEqual(q.units, "cm")
        self.assertEqual(q.uncertaintyType, "+|-")
        
        q = Quantity(([10.0,20.0,30.0],"cm","+|-",[1.0,2.0,3.0]))
        self.assertEqual(q.value, 0.0)
        self.assertEqual(q.uncertainty, 0.0)
        self.assertEqual(q.values[0], 0.1)
        self.assertEqual(q.values[1], 0.2)
        self.assertEqual(q.values[2], 0.3)
        self.assertEqual(q.uncertainties[0], 0.01)
        self.assertEqual(q.uncertainties[1], 0.02)
        self.assertEqual(q.uncertainties[2], 0.03)
        self.assertEqual(q.units, "cm")
        self.assertEqual(q.uncertaintyType, "+|-")
################################################################################

class TestCustomUnits(unittest.TestCase):
    """
    Contains unit tests of the custom units added to the quantities package.
    These tests ensure that, when these units are used in a Quantity object,
    the resulting value has the correct conversion factor (or at least that the
    conversion factor does not change!)
    """

    def testMolecules(self):
        q = Quantity(1.0,"molecules")
        self.assertAlmostEqual(q.value / (1.0/constants.Na), 1.0, 6)
    
    def testMolecule(self):
        q = Quantity(1.0,"molecule")
        self.assertAlmostEqual(q.value / (1.0/constants.Na), 1.0, 6)
        
    def testKcalPerMol(self):
        q = Quantity(1.0,"kcal/mol")
        self.assertAlmostEqual(q.value / (1.0*4184.), 1.0, 6)
    
    def testKJPerMol(self):
        q = Quantity(1.0,"kJ/mol")
        self.assertAlmostEqual(q.value / (1.0*1000.), 1.0, 6)
    
    def testJPerMol(self):
        q = Quantity(1.0,"J/kmol")
        self.assertAlmostEqual(q.value / (1.0/1000.), 1.0, 6)
    
################################################################################

class TestConstants(unittest.TestCase):
    """
    Contains unit tests of the Constants class. These tests ensure that the
    values of the constants are correct (or at least that they do not change!)
    """
    
    def testAvogadroConstant(self):
        self.assertAlmostEqual(constants.Na / pq.constants.N_A.simplified, 1.0, 6)

    def testBoltzmannConstant(self):
        self.assertAlmostEqual(constants.kB / pq.constants.Boltzmann_constant.simplified, 1.0, 6)
    
    def testGasLawConstant(self):
        self.assertAlmostEqual(constants.R / pq.constants.R.simplified, 1.0, 6)
    
    def testPlanckConstant(self):
        self.assertAlmostEqual(constants.h / pq.constants.Planck_constant.simplified, 1.0, 6)
    
    def testSpeedOfLight(self):
        self.assertAlmostEqual(constants.c / 299792458, 1.0, 6)
    
    def testPi(self):
        self.assertAlmostEqual(constants.pi / math.pi, 1.0, 6)

################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
