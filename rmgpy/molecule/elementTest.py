#!/usr/bin/env python
# encoding: utf-8

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2017 Prof. William H. Green (whgreen@mit.edu), 
#   Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

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

    def testGetElementIsotope(self):
        """
        Test that the rmgpy.elements.getElement() method works for isotopes.
        """
        self.assertTrue(isinstance(rmgpy.molecule.element.getElement('C', isotope=13), Element))
        self.assertTrue(isinstance(rmgpy.molecule.element.getElement(6, isotope=13), Element))

    def testChemkinName(self):
        """
        Test that retrieving the chemkin name of an element works.
        """
        d = rmgpy.molecule.element.getElement('H', isotope=2)
        self.assertEqual(d.chemkinName, 'D')

        c13 = rmgpy.molecule.element.getElement('C', isotope=13)
        self.assertEqual(c13.chemkinName, 'CI')

        o18 = rmgpy.molecule.element.getElement('O', isotope=18)
        self.assertEqual(o18.chemkinName, 'OI')
        
################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
