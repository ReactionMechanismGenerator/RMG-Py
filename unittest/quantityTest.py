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

class TestQuantity(unittest.TestCase):
    """
    Contains unit tests of the Quantity class. Most of the unit tests focus on
    the __init__() method, which is where most of the magic happens regarding
    Quantity objects. Tests of constructing Quantity objects via other means
    (especially via repr() and pickling/unpickling) are also performed.
    """
    
    def testValue(self):
        """
        Test that the correct Quantity object is created when only a number
        is provided.
        """
        q = Quantity(10.0)
        self.assertEqual(q.value, 10.0)
        self.assertEqual(q.uncertainty, 0.0)
        self.assertEqual(q.values, None)
        self.assertEqual(q.uncertainties, None)
        self.assertEqual(q.units, "")
        self.assertEqual(q.uncertaintyType, "")
        self.assertFalse(q.isArray())
        self.assertFalse(q.isUncertaintyAdditive())
        self.assertFalse(q.isUncertaintyMultiplicative())
        self.assertEqual(q.getConversionFactorToSI(), 1.0)
        self.assertEqual(q.getConversionFactorFromSI(), 1.0)
    
    def testValueWithUnits(self):
        """
        Test that the correct Quantity object is created when a number with
        units is provided.
        """
        q = Quantity(10.0,"cm")
        self.assertEqual(q.value, 0.1)
        self.assertEqual(q.uncertainty, 0.0)
        self.assertEqual(q.values, None)
        self.assertEqual(q.uncertainties, None)
        self.assertEqual(q.units, "cm")
        self.assertEqual(q.uncertaintyType, "")
        self.assertFalse(q.isArray())
        self.assertFalse(q.isUncertaintyAdditive())
        self.assertFalse(q.isUncertaintyMultiplicative())
        self.assertEqual(q.getConversionFactorToSI(), 0.01)
        self.assertEqual(q.getConversionFactorFromSI(), 100.0)
    
    def testValueWithUncertainty(self):
        """
        Test that the correct Quantity object is created when a number with
        units and an uncertainty is provided.
        """
        q = Quantity(10.0,"cm",2.0)
        self.assertEqual(q.value, 0.1)
        self.assertEqual(q.uncertainty, 0.02)
        self.assertEqual(q.values, None)
        self.assertEqual(q.uncertainties, None)
        self.assertEqual(q.units, "cm")
        self.assertEqual(q.uncertaintyType, "+|-")
        self.assertFalse(q.isArray())
        self.assertTrue(q.isUncertaintyAdditive())
        self.assertFalse(q.isUncertaintyMultiplicative())
        self.assertEqual(q.getConversionFactorToSI(), 0.01)
        self.assertEqual(q.getConversionFactorFromSI(), 100.0)
    
    def testValueWithAdditiveUncertainty(self):
        """
        Test that the correct Quantity object is created when a number with
        units and an additive uncertainty is provided.
        """
        q = Quantity(10.0,"cm","+|-",2.0)
        self.assertEqual(q.value, 0.1)
        self.assertEqual(q.uncertainty, 0.02)
        self.assertEqual(q.values, None)
        self.assertEqual(q.uncertainties, None)
        self.assertEqual(q.units, "cm")
        self.assertEqual(q.uncertaintyType, "+|-")
        self.assertFalse(q.isArray())
        self.assertTrue(q.isUncertaintyAdditive())
        self.assertFalse(q.isUncertaintyMultiplicative())
        self.assertEqual(q.getConversionFactorToSI(), 0.01)
        self.assertEqual(q.getConversionFactorFromSI(), 100.0)
    
    def testValueWithMultiplicativeUncertainty(self):
        """
        Test that the correct Quantity object is created when a number with
        units and a multiplicative uncertainty is provided.
        """
        q = Quantity(10.0,"cm","*|/",2.0)
        self.assertEqual(q.value, 0.1)
        self.assertEqual(q.uncertainty, 2.0)
        self.assertEqual(q.values, None)
        self.assertEqual(q.uncertainties, None)
        self.assertEqual(q.units, "cm")
        self.assertEqual(q.uncertaintyType, "*|/")
        self.assertFalse(q.isArray())
        self.assertFalse(q.isUncertaintyAdditive())
        self.assertTrue(q.isUncertaintyMultiplicative())
        self.assertEqual(q.getConversionFactorToSI(), 0.01)
        self.assertEqual(q.getConversionFactorFromSI(), 100.0)
    
    def testValues(self):
        """
        Test that the correct Quantity object is created when a set of numbers
        is provided.
        """
        q = Quantity([30.0,20.0,10.0])
        self.assertEqual(q.value, 0.0)
        self.assertEqual(q.uncertainty, 0.0)
        self.assertNotEqual(q.values, None)
        self.assertEqual(q.values[0], 30.0)
        self.assertEqual(q.values[1], 20.0)
        self.assertEqual(q.values[2], 10.0)
        self.assertEqual(q.uncertainties, None)
        self.assertEqual(q.units, "")
        self.assertEqual(q.uncertaintyType, "")
        self.assertTrue(q.isArray())
        self.assertFalse(q.isUncertaintyAdditive())
        self.assertFalse(q.isUncertaintyMultiplicative())
        self.assertEqual(q.getConversionFactorToSI(), 1.0)
        self.assertEqual(q.getConversionFactorFromSI(), 1.0)
    
    def testValuesWithUnits(self):
        """
        Test that the correct Quantity object is created when a set of numbers
        with units is provided.
        """
        q = Quantity([30.0,20.0,10.0],"cm")
        self.assertEqual(q.value, 0.0)
        self.assertEqual(q.uncertainty, 0.0)
        self.assertNotEqual(q.values, None)
        self.assertEqual(q.values[0], 0.3)
        self.assertEqual(q.values[1], 0.2)
        self.assertEqual(q.values[2], 0.1)
        self.assertEqual(q.uncertainties, None)
        self.assertEqual(q.units, "cm")
        self.assertEqual(q.uncertaintyType, "")
        self.assertTrue(q.isArray())
        self.assertFalse(q.isUncertaintyAdditive())
        self.assertFalse(q.isUncertaintyMultiplicative())
        self.assertEqual(q.getConversionFactorToSI(), 0.01)
        self.assertEqual(q.getConversionFactorFromSI(), 100.0)
    
    def testValuesWithUncertainty(self):
        """
        Test that the correct Quantity object is created when a set of numbers
        with a single uncertainty is provided.
        """
        q = Quantity([30.0,20.0,10.0],"cm","+|-",2.0)
        self.assertEqual(q.value, 0.0)
        self.assertEqual(q.uncertainty, 0.02)
        self.assertNotEqual(q.values, None)
        self.assertEqual(len(q.values), 3)
        self.assertEqual(q.values[0], 0.3)
        self.assertEqual(q.values[1], 0.2)
        self.assertEqual(q.values[2], 0.1)
        self.assertEqual(q.uncertainties, None)
        self.assertEqual(q.units, "cm")
        self.assertEqual(q.uncertaintyType, "+|-")
        self.assertTrue(q.isArray())
        self.assertTrue(q.isUncertaintyAdditive())
        self.assertFalse(q.isUncertaintyMultiplicative())
        self.assertEqual(q.getConversionFactorToSI(), 0.01)
        self.assertEqual(q.getConversionFactorFromSI(), 100.0)
    
    def testValuesWithUncertainties(self):
        """
        Test that the correct Quantity object is created when a set of numbers
        with a corresponding set of uncertainties is provided.
        """
        q = Quantity([10.0,20.0,30.0],"cm","+|-",[1.0,2.0,3.0])
        self.assertEqual(q.value, 0.0)
        self.assertEqual(q.uncertainty, 0.0)
        self.assertNotEqual(q.values, None)
        self.assertEqual(len(q.values), 3)
        self.assertEqual(q.values[0], 0.1)
        self.assertEqual(q.values[1], 0.2)
        self.assertEqual(q.values[2], 0.3)
        self.assertNotEqual(q.uncertainties, None)
        self.assertEqual(len(q.uncertainties), 3)
        self.assertEqual(q.uncertainties[0], 0.01)
        self.assertEqual(q.uncertainties[1], 0.02)
        self.assertEqual(q.uncertainties[2], 0.03)
        self.assertEqual(q.units, "cm")
        self.assertEqual(q.uncertaintyType, "+|-")
        self.assertTrue(q.isArray())
        self.assertTrue(q.isUncertaintyAdditive())
        self.assertFalse(q.isUncertaintyMultiplicative())
        self.assertEqual(q.getConversionFactorToSI(), 0.01)
        self.assertEqual(q.getConversionFactorFromSI(), 100.0)
    
    def testValueWithUncertainties(self):
        """
        Test that the correct Quantity object is created when a set of numbers
        with an incorrect set of uncertainties is provided.
        """
        try:
            Quantity(10.0,"cm","+|-",[1.0,2.0])
            self.fail("Expected a QuantityError to be raised as a result of providing multiple uncertainties for a single value.")
        except QuantityError:
            pass
    
    def testValuesWithInvalidUncertainties(self):
        """
        Test that the correct Quantity object is created when a set of numbers
        with an incorrect set of uncertainties is provided.
        """
        try:
            Quantity([10.0,20.0,30.0],"cm","+|-",[1.0,2.0])
            self.fail("Expected a QuantityError to be raised as a result of providing multiple uncertainties, but not one for each of the values.")
        except QuantityError:
            pass
    
    def testInvalidUncertaintyType(self):
        """
        Test that the correct Quantity object is created when a set of numbers
        with an incorrect set of uncertainties is provided.
        """
        try:
            Quantity([10.0,20.0,30.0],"cm","+-*/",1.0)
            self.fail("Expected a QuantityError to be raised as a result of providing an invalid uncertainty type.")
        except QuantityError:
            pass
    
    def testFromTuple(self):
        """
        Test that the correct Quantity object is created when the parameters
        are passed within a tuple. This is useful e.g. when parsing a 
        Python-formatted input file via the exec() command.
        """
        q = Quantity(([10.0,20.0,30.0],"cm","+|-",[1.0,2.0,3.0]))
        self.assertEqual(q.value, 0.0)
        self.assertEqual(q.uncertainty, 0.0)
        self.assertNotEqual(q.values, None)
        self.assertEqual(len(q.values), 3)
        self.assertEqual(q.values[0], 0.1)
        self.assertEqual(q.values[1], 0.2)
        self.assertEqual(q.values[2], 0.3)
        self.assertNotEqual(q.uncertainties, None)
        self.assertEqual(len(q.uncertainties), 3)
        self.assertEqual(q.uncertainties[0], 0.01)
        self.assertEqual(q.uncertainties[1], 0.02)
        self.assertEqual(q.uncertainties[2], 0.03)
        self.assertEqual(q.units, "cm")
        self.assertEqual(q.uncertaintyType, "+|-")
        self.assertTrue(q.isArray())
        self.assertTrue(q.isUncertaintyAdditive())
        self.assertFalse(q.isUncertaintyMultiplicative())
        self.assertEqual(q.getConversionFactorToSI(), 0.01)
        self.assertEqual(q.getConversionFactorFromSI(), 100.0)
    
    def testFromQuantity(self):
        """
        Test that the correct Quantity object is created when the parameters
        are passed via another Quantity object. This is useful e.g. when
        pickling and unpickling.
        """
        q0 = Quantity(([10.0,20.0,30.0],"cm","+|-",[1.0,2.0,3.0]))
        q = Quantity(q0)
        self.assertEqual(q.value, 0.0)
        self.assertEqual(q.uncertainty, 0.0)
        self.assertNotEqual(q.values, None)
        self.assertEqual(len(q.values), 3)
        self.assertEqual(q.values[0], 0.1)
        self.assertEqual(q.values[1], 0.2)
        self.assertEqual(q.values[2], 0.3)
        self.assertNotEqual(q.uncertainties, None)
        self.assertEqual(len(q.uncertainties), 3)
        self.assertEqual(q.uncertainties[0], 0.01)
        self.assertEqual(q.uncertainties[1], 0.02)
        self.assertEqual(q.uncertainties[2], 0.03)
        self.assertEqual(q.units, "cm")
        self.assertEqual(q.uncertaintyType, "+|-")
        self.assertTrue(q.isArray())
        self.assertTrue(q.isUncertaintyAdditive())
        self.assertFalse(q.isUncertaintyMultiplicative())
        self.assertEqual(q.getConversionFactorToSI(), 0.01)
        self.assertEqual(q.getConversionFactorFromSI(), 100.0)
    
    def testPickle(self):
        """
        Test that a Quantity object can be successfully pickled and
        unpickled with no loss of information.
        """
        import cPickle
        q0 = Quantity([10.0,20.0,30.0],"cm","+|-",[1.0,2.0,3.0])
        q = cPickle.loads(cPickle.dumps(q0))
        self.assertEqual(q.value, 0.0)
        self.assertEqual(q.uncertainty, 0.0)
        self.assertNotEqual(q.values, None)
        self.assertEqual(len(q.values), 3)
        self.assertEqual(q.values[0], 0.1)
        self.assertEqual(q.values[1], 0.2)
        self.assertEqual(q.values[2], 0.3)
        self.assertNotEqual(q.uncertainties, None)
        self.assertEqual(len(q.uncertainties), 3)
        self.assertEqual(q.uncertainties[0], 0.01)
        self.assertEqual(q.uncertainties[1], 0.02)
        self.assertEqual(q.uncertainties[2], 0.03)
        self.assertEqual(q.units, "cm")
        self.assertEqual(q.uncertaintyType, "+|-")
        self.assertTrue(q.isArray())
        self.assertTrue(q.isUncertaintyAdditive())
        self.assertFalse(q.isUncertaintyMultiplicative())
        self.assertEqual(q.getConversionFactorToSI(), 0.01)
        self.assertEqual(q.getConversionFactorFromSI(), 100.0)
    
    def testOutput(self):
        """
        Test that we can reconstruct a Quantity object from its repr()
        output with no loss of information.
        """
        q0 = Quantity([10.0,20.0,30.0],"cm","+|-",[1.0,2.0,3.0])
        exec('q = Quantity{0!r}'.format(q0))
        self.assertEqual(q.value, 0.0)
        self.assertEqual(q.uncertainty, 0.0)
        self.assertNotEqual(q.values, None)
        self.assertEqual(len(q.values), 3)
        self.assertEqual(q.values[0], 0.1)
        self.assertEqual(q.values[1], 0.2)
        self.assertEqual(q.values[2], 0.3)
        self.assertNotEqual(q.uncertainties, None)
        self.assertEqual(len(q.uncertainties), 3)
        self.assertEqual(q.uncertainties[0], 0.01)
        self.assertEqual(q.uncertainties[1], 0.02)
        self.assertEqual(q.uncertainties[2], 0.03)
        self.assertEqual(q.units, "cm")
        self.assertEqual(q.uncertaintyType, "+|-")
        self.assertTrue(q.isArray())
        self.assertTrue(q.isUncertaintyAdditive())
        self.assertFalse(q.isUncertaintyMultiplicative())
        self.assertEqual(q.getConversionFactorToSI(), 0.01)
        self.assertEqual(q.getConversionFactorFromSI(), 100.0)
    
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
