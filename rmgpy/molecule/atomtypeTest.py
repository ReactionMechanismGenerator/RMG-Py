#!/usr/bin/env python
# encoding: utf-8

"""
This module contains unit tests of the rmgpy.molecule.atomtype module.
"""

import unittest

import rmgpy.molecule
from rmgpy.molecule import atomtype, Molecule
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
    
    def testOutput(self):
        """
        Test that we can reconstruct an AtomType object from its repr()
        with no loss of information.
        """
        exec('atomType = rmgpy.molecule.atomtype.atomTypes[{0!r}]'.format(
                                    self.atomType.__repr__().split('"')[1]))
        return self.atomType.equivalent(atomType)
    
    def testEquivalent(self):
        """
        Test the AtomType.equivalent() method.
        """
        return self.atomType.equivalent(rmgpy.molecule.atomtype.atomTypes['Cd'])
    
    def testIsSpecficCaseOf(self):
        """
        Test the AtomType.isSpecificCaseOf() method.
        """
        return self.atomType.isSpecificCaseOf(rmgpy.molecule.atomtype.atomTypes['C'])
    
    def testSetActions(self):
        """
        Test the AtomType.setActions() method.
        """
        other = rmgpy.molecule.atomtype.atomTypes['R']
        other.setActions(self.atomType.incrementBond,
                               self.atomType.decrementBond,
                               self.atomType.formBond,
                               self.atomType.breakBond,
                               self.atomType.incrementRadical,
                               self.atomType.decrementRadical,
                               self.atomType.incrementLonePair,
                               self.atomType.decrementLonePair)
        self.assertEqual(self.atomType.incrementBond, other.incrementBond)
        self.assertEqual(self.atomType.decrementBond, other.decrementBond)
        self.assertEqual(self.atomType.formBond, other.formBond)
        self.assertEqual(self.atomType.breakBond, other.breakBond)
        self.assertEqual(self.atomType.incrementRadical, other.incrementRadical)
        self.assertEqual(self.atomType.decrementRadical, other.decrementRadical)

################################################################################

class TestGetAtomType(unittest.TestCase):
    """
    Contains unit tests of the getAtomType() method.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.mol1 = Molecule().fromSMILES('COO=CC=C=CC#C')
        self.mol2 = Molecule().fromSMILES('c1ccccc1')
        self.mol3 = Molecule().fromSMILES('[H]')
        self.mol4 = Molecule().fromSMILES(
                                'O=[Si][Si][Si]=[Si]=[Si][Si]#[Si]SS=S')
        self.mol5 = Molecule().fromSMILES('[N]')
        self.mol6 = Molecule().fromSMILES('[Ar]')
        self.mol7 = Molecule().fromSMILES('[He]')
        self.mol8 = Molecule().fromSMILES('[Ne]')
    
    def testCarbonTypes(self):
        """
        Test that getAtomType() returns appropriate carbon atom types.
        """
        self.assertEqual(getAtomType(self.mol1.atoms[0],
            self.mol1.getBonds(self.mol1.atoms[0])).label, 'Cs')
        self.assertEqual(getAtomType(self.mol1.atoms[4],
            self.mol1.getBonds(self.mol1.atoms[4])).label, 'Cd')
        self.assertEqual(getAtomType(self.mol1.atoms[5],
            self.mol1.getBonds(self.mol1.atoms[5])).label, 'Cdd')
        self.assertEqual(getAtomType(self.mol1.atoms[7],
            self.mol1.getBonds(self.mol1.atoms[7])).label, 'Ct')
        self.assertEqual(getAtomType(self.mol1.atoms[3],
            self.mol1.getBonds(self.mol1.atoms[3])).label, 'CO')
        self.assertEqual(getAtomType(self.mol2.atoms[0],
            self.mol2.getBonds(self.mol2.atoms[0])).label, 'Cb')
    
    def testHydrogenType(self):
        """
        Test that getAtomType() returns the hydrogen atom type.
        """
        self.assertEqual(getAtomType(self.mol3.atoms[0],
            self.mol3.getBonds(self.mol3.atoms[0])).label, 'H')
    
    def testOxygenTypes(self):
        """
        Test that getAtomType() returns appropriate oxygen atom types.
        """
        self.assertEqual(getAtomType(self.mol1.atoms[1],
            self.mol1.getBonds(self.mol1.atoms[1])).label, 'Os')
        self.assertEqual(getAtomType(self.mol1.atoms[2],
            self.mol1.getBonds(self.mol1.atoms[2])).label, 'Od')
    
    def testSiliconTypes(self):
        """
        Test that getAtomType() returns appropriate silicon atom types.
        """
        self.assertEqual(getAtomType(self.mol4.atoms[2],
            self.mol4.getBonds(self.mol4.atoms[2])).label, 'Sis')
        self.assertEqual(getAtomType(self.mol4.atoms[3],
            self.mol4.getBonds(self.mol4.atoms[3])).label, 'Sid')
        self.assertEqual(getAtomType(self.mol4.atoms[4],
            self.mol4.getBonds(self.mol4.atoms[4])).label, 'Sidd')
        self.assertEqual(getAtomType(self.mol4.atoms[6],
            self.mol4.getBonds(self.mol4.atoms[6])).label, 'Sit')
        self.assertEqual(getAtomType(self.mol4.atoms[1],
            self.mol4.getBonds(self.mol4.atoms[1])).label, 'SiO')
    
    def testSulfurTypes(self):
        """
        Test that getAtomType() returns appropriate sulfur atom types.
        """
        self.assertEqual(getAtomType(self.mol4.atoms[8],
            self.mol4.getBonds(self.mol4.atoms[8])).label, 'Ss')
        self.assertEqual(getAtomType(self.mol4.atoms[9],
            self.mol4.getBonds(self.mol4.atoms[9])).label, 'Sd')
    
    def testNoneTypes(self):
        """
        Test that getAtomType() returns appropriate NoneTypes.
        """
        self.assertIsNone(getAtomType(self.mol5.atoms[0],
            self.mol5.getBonds(self.mol5.atoms[0])))
        self.assertIsNone(getAtomType(self.mol6.atoms[0],
            self.mol6.getBonds(self.mol6.atoms[0])))
        self.assertIsNone(getAtomType(self.mol7.atoms[0],
            self.mol7.getBonds(self.mol7.atoms[0])))
        self.assertIsNone(getAtomType(self.mol8.atoms[0],
            self.mol8.getBonds(self.mol8.atoms[0])))

################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
