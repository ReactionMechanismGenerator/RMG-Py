#!/usr/bin/env python
# encoding: utf-8

"""
This module contains unit tests of the rmgpy.element module.
"""

import unittest

from rmgpy.molecule.element import Element
import rmgpy.molecule.element

################################################################################

class TestElement(unittest.TestCase):
    """
    Contains unit tests of the Element class.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.element = rmgpy.molecule.element.C
        
    def testPickle(self):
        """
        Test that an Element object can be successfully pickled and
        unpickled with no loss of information.
        """
        import cPickle
        element = cPickle.loads(cPickle.dumps(self.element))
        self.assertEqual(self.element.number, element.number)
        self.assertEqual(self.element.symbol, element.symbol)
        self.assertEqual(self.element.name, element.name)
        self.assertEqual(self.element.mass, element.mass)
    
    def testOutput(self):
        """
        Test that we can reconstruct an Element object from its repr()
        output with no loss of information.
        """
        exec('element = {0!r}'.format(self.element))
        self.assertEqual(self.element.number, element.number)
        self.assertEqual(self.element.symbol, element.symbol)
        self.assertEqual(self.element.name, element.name)
        self.assertEqual(self.element.mass, element.mass)
    
    def testGetElement(self):
        """
        Test the rmgpy.elements.getElement() method.
        """
        self.assertTrue(rmgpy.molecule.element.getElement(6) is self.element)
        self.assertTrue(rmgpy.molecule.element.getElement('C') is self.element)
        
################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
