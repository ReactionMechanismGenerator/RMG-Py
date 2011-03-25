#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
import unittest

from rmgpy.chem.constants import *

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

        q = Quantity(10.0,"cm")
        self.assertEqual(q.value, 0.1)
        self.assertEqual(q.uncertainty, 0.0)
        self.assertEqual(q.values, None)
        self.assertEqual(q.uncertainties, None)
        self.assertEqual(q.units, "cm")
        self.assertEqual(q.uncertaintyType, "")

        q = Quantity(10.0,"cm","+|-",2.0)
        self.assertEqual(q.value, 0.1)
        self.assertEqual(q.uncertainty, 0.02)
        self.assertEqual(q.values, None)
        self.assertEqual(q.uncertainties, None)
        self.assertEqual(q.units, "cm")
        self.assertEqual(q.uncertaintyType, "+|-")
        
        q = Quantity(10.0,"cm","*|/",2.0)
        self.assertEqual(q.value, 0.1)
        self.assertEqual(q.uncertainty, 2.0)
        self.assertEqual(q.values, None)
        self.assertEqual(q.uncertainties, None)
        self.assertEqual(q.units, "cm")
        self.assertEqual(q.uncertaintyType, "*|/")
        
        q = Quantity([30.0,20.0,10.0],"cm","+|-",2.0)
        self.assertEqual(q.value, 0.0)
        self.assertEqual(q.uncertainty, 0.02)
        self.assertEqual(q.values[0], 0.3)
        self.assertEqual(q.values[1], 0.2)
        self.assertEqual(q.values[2], 0.1)
        self.assertEqual(q.uncertainties, None)
        self.assertEqual(q.units, "cm")
        self.assertEqual(q.uncertaintyType, "+|-")
        
        q = Quantity([10.0,20.0,30.0],"cm","+|-",[1.0,2.0,3.0])
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

if __name__ == '__main__':
    unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )
