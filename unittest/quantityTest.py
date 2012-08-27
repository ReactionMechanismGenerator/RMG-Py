#!/usr/bin/env python
# encoding: utf-8

"""
This module contains unit tests of the rmgpy.quantity module.
"""

import unittest
import math
import quantities as pq

import rmgpy.constants as constants
from rmgpy.quantity import *

################################################################################

class TestScalarQuantity(unittest.TestCase):
    """
    Contains unit tests of the :class:`ScalarQuantity` class.
    """
    
    def setUp(self):
        self.value = 1.0
        self.units = 'kJ/mol'
        self.uncertaintyType = '*|/'
        self.uncertainty = 1.5
    
    def test_init_simple1(self):
        """
        Test the ScalarQuantity.__init__() method for a dimensionless quantity.
        """
        quantity = Quantity(self.value)
        self.assertTrue(isinstance(quantity, ScalarQuantity))
        self.assertEqual(quantity.value, self.value)
        self.assertEqual(quantity.units, '')
        self.assertEqual(quantity.uncertaintyType, '+|-')
        self.assertEqual(quantity.uncertainty, 0.0)
        self.assertEqual(quantity.value_si, self.value)
        
    def test_init_simple2(self):
        """
        Test the ScalarQuantity.__init__() method for a quantity with units.
        """
        quantity = Quantity(self.value, self.units)
        self.assertTrue(isinstance(quantity, ScalarQuantity))
        self.assertEqual(quantity.value, self.value)
        self.assertEqual(quantity.units, self.units)
        self.assertEqual(quantity.uncertaintyType, '+|-')
        self.assertEqual(quantity.uncertainty, 0.0)
        self.assertEqual(quantity.value_si, self.value * 1000.)
        
    def test_init_simple3(self):
        """
        Test the ScalarQuantity.__init__() method for a quantity with units and
        (additive) uncertainty.
        """
        quantity = Quantity(self.value, self.units, self.uncertainty)
        self.assertTrue(isinstance(quantity, ScalarQuantity))
        self.assertEqual(quantity.value, self.value)
        self.assertEqual(quantity.units, self.units)
        self.assertEqual(quantity.uncertaintyType, '+|-')
        self.assertEqual(quantity.uncertainty, self.uncertainty)
        self.assertEqual(quantity.value_si, self.value * 1000.)
        
    def test_init_simple4(self):
        """
        Test the ScalarQuantity.__init__() method for a quantity with units and
        multiplicative uncertainty.
        """
        quantity = Quantity(self.value,  self.units, self.uncertaintyType, self.uncertainty)
        self.assertTrue(isinstance(quantity, ScalarQuantity))
        self.assertEqual(quantity.value, self.value)
        self.assertEqual(quantity.units, self.units)
        self.assertEqual(quantity.uncertaintyType, self.uncertaintyType)
        self.assertEqual(quantity.uncertainty, self.uncertainty)
        self.assertEqual(quantity.value_si, self.value * 1000.)

    def test_init_tuple1(self):
        """
        Test the ScalarQuantity.__init__() method for a dimensionless quantity, 
        specified using a tuple of arguments.
        """
        quantity = Quantity((self.value,))
        self.assertTrue(isinstance(quantity, ScalarQuantity))
        self.assertEqual(quantity.value, self.value)
        self.assertEqual(quantity.units, '')
        self.assertEqual(quantity.uncertaintyType, '+|-')
        self.assertEqual(quantity.uncertainty, 0.0)
        self.assertEqual(quantity.value_si, self.value)
        
    def test_init_tuple2(self):
        """
        Test the ScalarQuantity.__init__() method for a quantity with units, 
        specified using a tuple of arguments.
        """
        quantity = Quantity((self.value, self.units))
        self.assertTrue(isinstance(quantity, ScalarQuantity))
        self.assertEqual(quantity.value, self.value)
        self.assertEqual(quantity.units, self.units)
        self.assertEqual(quantity.uncertaintyType, '+|-')
        self.assertEqual(quantity.uncertainty, 0.0)
        self.assertEqual(quantity.value_si, self.value * 1000.)
        
    def test_init_tuple3(self):
        """
        Test the ScalarQuantity.__init__() method for a quantity with units and
        (additive) uncertainty, specified using a tuple of arguments.
        """
        quantity = Quantity((self.value, self.units, self.uncertainty))
        self.assertTrue(isinstance(quantity, ScalarQuantity))
        self.assertEqual(quantity.value, self.value)
        self.assertEqual(quantity.units, self.units)
        self.assertEqual(quantity.uncertaintyType, '+|-')
        self.assertEqual(quantity.uncertainty, self.uncertainty)
        self.assertEqual(quantity.value_si, self.value * 1000.)
        
    def test_init_tuple4(self):
        """
        Test the ScalarQuantity.__init__() method for a quantity with units and
        multiplicative uncertainty, specified using a tuple of arguments.
        """
        quantity = Quantity((self.value,  self.units, self.uncertaintyType, self.uncertainty))
        self.assertTrue(isinstance(quantity, ScalarQuantity))
        self.assertEqual(quantity.value, self.value)
        self.assertEqual(quantity.units, self.units)
        self.assertEqual(quantity.uncertaintyType, self.uncertaintyType)
        self.assertEqual(quantity.uncertainty, self.uncertainty)
        self.assertEqual(quantity.value_si, self.value * 1000.)
        
    def test_init_keyword1(self):
        """
        Test the ScalarQuantity.__init__() method for a dimensionless quantity, 
        specified using keyword arguments.
        """
        quantity = Quantity(value=self.value)
        self.assertTrue(isinstance(quantity, ScalarQuantity))
        self.assertEqual(quantity.value, self.value)
        self.assertEqual(quantity.units, '')
        self.assertEqual(quantity.uncertaintyType, '+|-')
        self.assertEqual(quantity.uncertainty, 0.0)
        self.assertEqual(quantity.value_si, self.value)
        
    def test_init_keyword2(self):
        """
        Test the ScalarQuantity.__init__() method for a quantity with units, 
        specified using keyword arguments.
        """
        quantity = Quantity(self.value, units=self.units)
        self.assertTrue(isinstance(quantity, ScalarQuantity))
        self.assertEqual(quantity.value, self.value)
        self.assertEqual(quantity.units, self.units)
        self.assertEqual(quantity.uncertaintyType, '+|-')
        self.assertEqual(quantity.uncertainty, 0.0)
        self.assertEqual(quantity.value_si, self.value * 1000.)
        
    def test_init_keyword3(self):
        """
        Test the ScalarQuantity.__init__() method for a quantity with units and
        (additive) uncertainty, specified using keyword arguments.
        """
        quantity = Quantity(self.value, self.units, uncertainty=self.uncertainty)
        self.assertTrue(isinstance(quantity, ScalarQuantity))
        self.assertEqual(quantity.value, self.value)
        self.assertEqual(quantity.units, self.units)
        self.assertEqual(quantity.uncertaintyType, '+|-')
        self.assertEqual(quantity.uncertainty, self.uncertainty)
        self.assertEqual(quantity.value_si, self.value * 1000.)
        
    def test_init_keyword4(self):
        """
        Test the ScalarQuantity.__init__() method for a quantity with units and
        multiplicative uncertainty, specified using keyword arguments.
        """
        quantity = Quantity(self.value,  self.units, self.uncertainty, uncertaintyType=self.uncertaintyType)
        self.assertTrue(isinstance(quantity, ScalarQuantity))
        self.assertEqual(quantity.value, self.value)
        self.assertEqual(quantity.units, self.units)
        self.assertEqual(quantity.uncertaintyType, self.uncertaintyType)
        self.assertEqual(quantity.uncertainty, self.uncertainty)
        self.assertEqual(quantity.value_si, self.value * 1000.)

    def test_str1(self):
        """
        Test the ScalarQuantity.__str__() method for a dimensionless quantity.
        """
        quantity = Quantity(self.value)
        self.assertEqual(str(quantity), '1')
        
    def test_str2(self):
        """
        Test the ScalarQuantity.__str__() method for a quantity with units.
        """
        quantity = Quantity(self.value, self.units)
        self.assertEqual(str(quantity), '1 kJ/mol')
        
    def test_str3(self):
        """
        Test the ScalarQuantity.__str__() method for a quantity with units and
        (additive) uncertainty.
        """
        quantity = Quantity(self.value, self.units, self.uncertainty)
        self.assertEqual(str(quantity), '1 +|- 1.5 kJ/mol')
        
    def test_str4(self):
        """
        Test the ScalarQuantity.__str__() method for a quantity with units and
        multiplicative uncertainty.
        """
        quantity = Quantity(self.value,  self.units, self.uncertaintyType, self.uncertainty)
        self.assertEqual(str(quantity), '1 *|/ 1.5 kJ/mol')

    def test_repr1(self):
        """
        Test the ScalarQuantity.__repr__() method for a dimensionless quantity.
        """
        quantity = Quantity(self.value)
        self.assertEqual(repr(quantity), '1')
        
    def test_repr2(self):
        """
        Test the ScalarQuantity.__repr__() method for a quantity with units.
        """
        quantity = Quantity(self.value, self.units)
        self.assertEqual(repr(quantity), "(1,'kJ/mol')")
        
    def test_repr3(self):
        """
        Test the ScalarQuantity.__repr__() method for a quantity with units and
        (additive) uncertainty.
        """
        quantity = Quantity(self.value, self.units, self.uncertainty)
        self.assertEqual(repr(quantity), "(1,'kJ/mol','+|-',1.5)")
        
    def test_repr4(self):
        """
        Test the ScalarQuantity.__repr__() method for a quantity with units and
        multiplicative uncertainty.
        """
        quantity = Quantity(self.value,  self.units, self.uncertaintyType, self.uncertainty)
        self.assertEqual(repr(quantity), "(1,'kJ/mol','*|/',1.5)")

################################################################################

class TestArrayQuantity(unittest.TestCase):
    """
    Contains unit tests of the :class:`ArrayQuantity` class.
    """
    
    def setUp(self):
        self.value = numpy.array([1.0,2.0,3.0])
        self.units = 'kJ/mol'
        self.uncertaintyType = '*|/'
        self.uncertainty = numpy.array([1.5,2.0,2.5])
    
    def test_init_simple1(self):
        """
        Test the ArrayQuantity.__init__() method for a dimensionless quantity.
        """
        quantity = Quantity(self.value)
        self.assertTrue(isinstance(quantity, ArrayQuantity))
        self.assertEqual(quantity.units, '')
        self.assertEqual(quantity.uncertaintyType, '+|-')
        self.assertEqual(quantity.value.shape, self.value.shape)
        for i in range(quantity.value.shape[0]):
            self.assertEqual(quantity.value[i], self.value[i])
            self.assertEqual(quantity.uncertainty[i], 0.0)
            self.assertEqual(quantity.value_si[i], self.value[i])
        
    def test_init_simple2(self):
        """
        Test the ArrayQuantity.__init__() method for a quantity with units.
        """
        quantity = Quantity(self.value, self.units)
        self.assertTrue(isinstance(quantity, ArrayQuantity))
        self.assertEqual(quantity.units, self.units)
        self.assertEqual(quantity.uncertaintyType, '+|-')
        self.assertEqual(quantity.value.shape, self.value.shape)
        for i in range(quantity.value.shape[0]):
            self.assertEqual(quantity.value[i], self.value[i])
            self.assertEqual(quantity.uncertainty[i], 0.0)
            self.assertEqual(quantity.value_si[i], self.value[i] * 1000.)
        
    def test_init_simple3(self):
        """
        Test the ArrayQuantity.__init__() method for a quantity with units and
        (additive) uncertainty.
        """
        quantity = Quantity(self.value, self.units, self.uncertainty)
        self.assertTrue(isinstance(quantity, ArrayQuantity))
        self.assertEqual(quantity.units, self.units)
        self.assertEqual(quantity.uncertaintyType, '+|-')
        self.assertEqual(quantity.value.shape, self.value.shape)
        for i in range(quantity.value.shape[0]):
            self.assertEqual(quantity.value[i], self.value[i])
            self.assertEqual(quantity.uncertainty[i], self.uncertainty[i])
            self.assertEqual(quantity.value_si[i], self.value[i] * 1000.)
        
    def test_init_simple4(self):
        """
        Test the ArrayQuantity.__init__() method for a quantity with units and
        multiplicative uncertainty.
        """
        quantity = Quantity(self.value,  self.units, self.uncertaintyType, self.uncertainty)
        self.assertTrue(isinstance(quantity, ArrayQuantity))
        self.assertEqual(quantity.units, self.units)
        self.assertEqual(quantity.uncertaintyType, self.uncertaintyType)
        self.assertEqual(quantity.value.shape, self.value.shape)
        for i in range(quantity.value.shape[0]):
            self.assertEqual(quantity.value[i], self.value[i])
            self.assertEqual(quantity.uncertainty[i], self.uncertainty[i])
            self.assertEqual(quantity.value_si[i], self.value[i] * 1000.)

    def test_init_simple5(self):
        """
        Test the ArrayQuantity.__init__() method for a quantity with units and
        (additive) uncertainty.
        """
        quantity = Quantity(self.value, self.units, 2.0)
        self.assertTrue(isinstance(quantity, ArrayQuantity))
        self.assertEqual(quantity.units, self.units)
        self.assertEqual(quantity.uncertaintyType, '+|-')
        self.assertEqual(quantity.value.shape, self.value.shape)
        for i in range(quantity.value.shape[0]):
            self.assertEqual(quantity.value[i], self.value[i])
            self.assertEqual(quantity.uncertainty[i], 2.0)
            self.assertEqual(quantity.value_si[i], self.value[i] * 1000.)

    def test_init_simple6(self):
        """
        Test the ArrayQuantity.__init__() method for a quantity with units and
        (additive) uncertainty.
        """
        try:
            quantity = Quantity(self.value, self.units, [1.5,2.0])
            self.fail("Expected a QuantityError to be raised as a result of providing multiple uncertainties, but not one for each of the values.")
        except QuantityError:
            pass

    def test_init_tuple1(self):
        """
        Test the ArrayQuantity.__init__() method for a dimensionless quantity, 
        specified using a tuple of arguments.
        """
        quantity = Quantity((self.value,))
        self.assertTrue(isinstance(quantity, ArrayQuantity))
        self.assertEqual(quantity.units, '')
        self.assertEqual(quantity.uncertaintyType, '+|-')
        self.assertEqual(quantity.value.shape, self.value.shape)
        for i in range(quantity.value.shape[0]):
            self.assertEqual(quantity.value[i], self.value[i])
            self.assertEqual(quantity.uncertainty[i], 0.0)
            self.assertEqual(quantity.value_si[i], self.value[i])
        
    def test_init_tuple2(self):
        """
        Test the ArrayQuantity.__init__() method for a quantity with units, 
        specified using a tuple of arguments.
        """
        quantity = Quantity((self.value, self.units))
        self.assertTrue(isinstance(quantity, ArrayQuantity))
        self.assertEqual(quantity.units, self.units)
        self.assertEqual(quantity.uncertaintyType, '+|-')
        self.assertEqual(quantity.value.shape, self.value.shape)
        for i in range(quantity.value.shape[0]):
            self.assertEqual(quantity.value[i], self.value[i])
            self.assertEqual(quantity.uncertainty[i], 0.0)
            self.assertEqual(quantity.value_si[i], self.value[i] * 1000.)
        
    def test_init_tuple3(self):
        """
        Test the ArrayQuantity.__init__() method for a quantity with units and
        (additive) uncertainty, specified using a tuple of arguments.
        """
        quantity = Quantity((self.value, self.units, self.uncertainty))
        self.assertTrue(isinstance(quantity, ArrayQuantity))
        self.assertEqual(quantity.units, self.units)
        self.assertEqual(quantity.uncertaintyType, '+|-')
        self.assertEqual(quantity.value.shape, self.value.shape)
        for i in range(quantity.value.shape[0]):
            self.assertEqual(quantity.value[i], self.value[i])
            self.assertEqual(quantity.uncertainty[i], self.uncertainty[i])
            self.assertEqual(quantity.value_si[i], self.value[i] * 1000.)
        
    def test_init_tuple4(self):
        """
        Test the ArrayQuantity.__init__() method for a quantity with units and
        multiplicative uncertainty, specified using a tuple of arguments.
        """
        quantity = Quantity((self.value,  self.units, self.uncertaintyType, self.uncertainty))
        self.assertTrue(isinstance(quantity, ArrayQuantity))
        self.assertEqual(quantity.units, self.units)
        self.assertEqual(quantity.uncertaintyType, self.uncertaintyType)
        self.assertEqual(quantity.value.shape, self.value.shape)
        for i in range(quantity.value.shape[0]):
            self.assertEqual(quantity.value[i], self.value[i])
            self.assertEqual(quantity.uncertainty[i], self.uncertainty[i])
            self.assertEqual(quantity.value_si[i], self.value[i] * 1000.)
        
    def test_init_keyword1(self):
        """
        Test the ArrayQuantity.__init__() method for a dimensionless quantity, 
        specified using keyword arguments.
        """
        quantity = Quantity(value=self.value)
        self.assertTrue(isinstance(quantity, ArrayQuantity))
        self.assertEqual(quantity.units, '')
        self.assertEqual(quantity.uncertaintyType, '+|-')
        self.assertEqual(quantity.value.shape, self.value.shape)
        for i in range(quantity.value.shape[0]):
            self.assertEqual(quantity.value[i], self.value[i])
            self.assertEqual(quantity.uncertainty[i], 0.0)
            self.assertEqual(quantity.value_si[i], self.value[i])
        
    def test_init_keyword2(self):
        """
        Test the ArrayQuantity.__init__() method for a quantity with units, 
        specified using keyword arguments.
        """
        quantity = Quantity(self.value, units=self.units)
        self.assertTrue(isinstance(quantity, ArrayQuantity))
        self.assertEqual(quantity.units, self.units)
        self.assertEqual(quantity.uncertaintyType, '+|-')
        self.assertEqual(quantity.value.shape, self.value.shape)
        for i in range(quantity.value.shape[0]):
            self.assertEqual(quantity.value[i], self.value[i])
            self.assertEqual(quantity.uncertainty[i], 0.0)
            self.assertEqual(quantity.value_si[i], self.value[i] * 1000.)
        
    def test_init_keyword3(self):
        """
        Test the ArrayQuantity.__init__() method for a quantity with units and
        (additive) uncertainty, specified using keyword arguments.
        """
        quantity = Quantity(self.value, self.units, uncertainty=self.uncertainty)
        self.assertTrue(isinstance(quantity, ArrayQuantity))
        self.assertEqual(quantity.units, self.units)
        self.assertEqual(quantity.uncertaintyType, '+|-')
        self.assertEqual(quantity.value.shape, self.value.shape)
        for i in range(quantity.value.shape[0]):
            self.assertEqual(quantity.value[i], self.value[i])
            self.assertEqual(quantity.uncertainty[i], self.uncertainty[i])
            self.assertEqual(quantity.value_si[i], self.value[i] * 1000.)
        
    def test_init_keyword4(self):
        """
        Test the ArrayQuantity.__init__() method for a quantity with units and
        multiplicative uncertainty, specified using keyword arguments.
        """
        quantity = Quantity(self.value,  self.units, self.uncertainty, uncertaintyType=self.uncertaintyType)
        self.assertTrue(isinstance(quantity, ArrayQuantity))
        self.assertEqual(quantity.units, self.units)
        self.assertEqual(quantity.uncertaintyType, self.uncertaintyType)
        self.assertEqual(quantity.value.shape, self.value.shape)
        for i in range(quantity.value.shape[0]):
            self.assertEqual(quantity.value[i], self.value[i])
            self.assertEqual(quantity.uncertainty[i], self.uncertainty[i])
            self.assertEqual(quantity.value_si[i], self.value[i] * 1000.)

    def test_str1(self):
        """
        Test the ArrayQuantity.__str__() method for a dimensionless quantity.
        """
        quantity = Quantity(self.value)
        self.assertEqual(str(quantity), '[1,2,3]')
        
    def test_str2(self):
        """
        Test the ArrayQuantity.__str__() method for a quantity with units.
        """
        quantity = Quantity(self.value, self.units)
        self.assertEqual(str(quantity), '[1,2,3] kJ/mol')
        
    def test_str3(self):
        """
        Test the ArrayQuantity.__str__() method for a quantity with units and
        (additive) uncertainty.
        """
        quantity = Quantity(self.value, self.units, self.uncertainty)
        self.assertEqual(str(quantity), '[1,2,3] +|- [1.5,2,2.5] kJ/mol')
        
    def test_str4(self):
        """
        Test the ArrayQuantity.__str__() method for a quantity with units and
        multiplicative uncertainty.
        """
        quantity = Quantity(self.value,  self.units, self.uncertaintyType, self.uncertainty)
        self.assertEqual(str(quantity), '[1,2,3] *|/ [1.5,2,2.5] kJ/mol')

    def test_repr1(self):
        """
        Test the ArrayQuantity.__repr__() method for a dimensionless quantity.
        """
        quantity = Quantity(self.value)
        self.assertEqual(repr(quantity), '[1,2,3]')
        
    def test_repr2(self):
        """
        Test the ArrayQuantity.__repr__() method for a quantity with units.
        """
        quantity = Quantity(self.value, self.units)
        self.assertEqual(repr(quantity), "([1,2,3],'kJ/mol')")
        
    def test_repr3(self):
        """
        Test the ArrayQuantity.__repr__() method for a quantity with units and
        (additive) uncertainty.
        """
        quantity = Quantity(self.value, self.units, self.uncertainty)
        self.assertEqual(repr(quantity), "([1,2,3],'kJ/mol','+|-',[1.5,2,2.5])")
        
    def test_repr4(self):
        """
        Test the ArrayQuantity.__repr__() method for a quantity with units and
        multiplicative uncertainty.
        """
        quantity = Quantity(self.value,  self.units, self.uncertaintyType, self.uncertainty)
        self.assertEqual(repr(quantity), "([1,2,3],'kJ/mol','*|/',[1.5,2,2.5])")

    def test_equals(self):
        """
        Tests the ``equals`` function within quantity to make sure that it rejects
        nonidentical quantity objects.
        """
        q1 = Quantity((10,"cm","+|-",1.0))
        q1b = Quantity((10,"cm","+|-",1.0))
        q2 = Quantity((10,"kJ","+|-",1.0))
        q2b = Quantity((10,"kJ","+|-",1.0))
        q3 = Quantity(([10.0,20.0,30.0],"cm","+|-",[1.0,2.0,3.0]))
        q3b = Quantity(([10.0,20.0,30.0],"cm","+|-",[1.0,2.0,3.0]))
        q4 = Quantity(([10.0,20.0,30.0],"cm","+|-",1.0))
        q4b = Quantity(([10.0,20.0,30.0],"cm","+|-",1.0))

        self.assertTrue(q1.equals(q1b))
        self.assertTrue(q2.equals(q2b))
        self.assertTrue(q3.equals(q3b))
        self.assertTrue(q4.equals(q4b))
        self.assertFalse(q1.equals(q2))
        self.assertFalse(q1.equals(q3))
        self.assertFalse(q1.equals(q4))
        self.assertFalse(q3.equals(q4))
    
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
        self.assertAlmostEqual(q.value_si / (1.0/constants.Na), 1.0, 6)
    
    def testMolecule(self):
        q = Quantity(1.0,"molecule")
        self.assertAlmostEqual(q.value_si / (1.0/constants.Na), 1.0, 6)
        
    def testKcalPerMol(self):
        q = Quantity(1.0,"kcal/mol")
        self.assertAlmostEqual(q.value_si / (1.0*4184.), 1.0, 6)
    
    def testKJPerMol(self):
        q = Quantity(1.0,"kJ/mol")
        self.assertAlmostEqual(q.value_si / (1.0*1000.), 1.0, 6)
    
    def testJPerMol(self):
        q = Quantity(1.0,"J/kmol")
        self.assertAlmostEqual(q.value_si / (1.0/1000.), 1.0, 6)

################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
