#!/usr/bin/env python
# encoding: utf-8

"""
This module contains unit tests of the rmgpy.atomtype module.
"""

import unittest

import rmgpy.molecule.atomtype
from rmgpy.molecule.atomtype import AtomType, getAtomType

################################################################################

class TestAtomType(unittest.TestCase):
    """
    Contains unit tests of the AtomType class.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.atomType = rmgpy.molecule.atomtype.atomTypes['Cd']
        
    def testPickle(self):
        """
        Test that an AtomType object can be successfully pickled and
        unpickled with no loss of information.
        """
        import cPickle
        atomType = cPickle.loads(cPickle.dumps(self.atomType))
        self.assertEqual(self.atomType.label, atomType.label)
        self.assertEqual(len(self.atomType.generic), len(atomType.generic))
        for item1, item2 in zip(self.atomType.generic, atomType.generic):
            self.assertEqual(item1.label, item2.label)
        self.assertEqual(len(self.atomType.specific), len(atomType.specific))
        for item1, item2 in zip(self.atomType.specific, atomType.specific):
            self.assertEqual(item1.label, item2.label)
        self.assertEqual(len(self.atomType.incrementBond), len(atomType.incrementBond))
        for item1, item2 in zip(self.atomType.incrementBond, atomType.incrementBond):
            self.assertEqual(item1.label, item2.label)
        self.assertEqual(len(self.atomType.decrementBond), len(atomType.decrementBond))
        for item1, item2 in zip(self.atomType.decrementBond, atomType.decrementBond):
            self.assertEqual(item1.label, item2.label)
        self.assertEqual(len(self.atomType.formBond), len(atomType.formBond))
        for item1, item2 in zip(self.atomType.formBond, atomType.formBond):
            self.assertEqual(item1.label, item2.label)
        self.assertEqual(len(self.atomType.breakBond), len(atomType.breakBond))
        for item1, item2 in zip(self.atomType.breakBond, atomType.breakBond):
            self.assertEqual(item1.label, item2.label)
        self.assertEqual(len(self.atomType.incrementRadical), len(atomType.incrementRadical))
        for item1, item2 in zip(self.atomType.incrementRadical, atomType.incrementRadical):
            self.assertEqual(item1.label, item2.label)
        self.assertEqual(len(self.atomType.decrementRadical), len(atomType.decrementRadical))
        for item1, item2 in zip(self.atomType.decrementRadical, atomType.decrementRadical):
            self.assertEqual(item1.label, item2.label)
        
################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
