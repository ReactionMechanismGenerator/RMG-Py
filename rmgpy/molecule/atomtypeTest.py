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
        other = rmgpy.molecule.atomtype.AtomType('Test', generic=['R'], specific=[])
        other.setActions(self.atomType.incrementBond,
                               self.atomType.decrementBond,
                               self.atomType.formBond,
                               self.atomType.breakBond,
                               self.atomType.incrementRadical,
                               self.atomType.decrementRadical)
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
        self.mol1 = Molecule().fromSMILES('COC(=O)CC=C=CC#C')
        # self.mol2 = Molecule().fromSMILES('c1ccccc1')
        ## the fromSMILES method currently Kekulizes, so to test Benzene we use fromAdjacencyList
        self.mol2 = Molecule().fromAdjacencyList('''1 C 0 {2,B} {6,B}
                                                    2 C 0 {1,B} {3,B}
                                                    3 C 0 {2,B} {4,B}
                                                    4 C 0 {3,B} {5,B}
                                                    5 C 0 {4,B} {6,B}
                                                    6 C 0 {1,B} {5,B}''')
        self.mol3 = Molecule().fromSMILES('[H]')
        self.mol4 = Molecule().fromSMILES(
                                'O=[Si][Si][Si]=[Si]=[Si][Si]#[Si]SS=S')
        self.mol5 = Molecule().fromSMILES('[N]')
        self.mol6 = Molecule().fromSMILES('[Ar]')
        self.mol7 = Molecule().fromSMILES('[He]')
        self.mol8 = Molecule().fromSMILES('[Ne]')
    
    def atomType(self, mol, atomID):
        atom = mol.atoms[atomID]
        type = getAtomType(atom, mol.getBonds(atom))
        if type is None:
            return type
        else:
            return type.label

    def testCarbonTypes(self):
        """
        Test that getAtomType() returns appropriate carbon atom types.
        """
        self.assertEqual(self.atomType(self.mol1, 0), 'Cs')
        self.assertEqual(self.atomType(self.mol1, 5), 'Cd')
        self.assertEqual(self.atomType(self.mol1, 6), 'Cdd')
        self.assertEqual(self.atomType(self.mol1, 8), 'Ct')
        self.assertEqual(self.atomType(self.mol1, 2), 'CO')
        self.assertEqual(self.atomType(self.mol2, 0), 'Cb')
    
    def testHydrogenType(self):
        """
        Test that getAtomType() returns the hydrogen atom type.
        """
        self.assertEqual(self.atomType(self.mol3, 0), 'H')
    
    def testOxygenTypes(self):
        """
        Test that getAtomType() returns appropriate oxygen atom types.
        """
        self.assertEqual(self.atomType(self.mol1, 1), 'Os')
        self.assertEqual(self.atomType(self.mol1, 3), 'Od')
    
    def testSiliconTypes(self):
        """
        Test that getAtomType() returns appropriate silicon atom types.
        """
        self.assertEqual(self.atomType(self.mol4, 2), 'Sis')
        self.assertEqual(self.atomType(self.mol4, 3), 'Sid')
        self.assertEqual(self.atomType(self.mol4, 4), 'Sidd')
        self.assertEqual(self.atomType(self.mol4, 6), 'Sit')
        self.assertEqual(self.atomType(self.mol4, 1), 'SiO')
    
    def testSulfurTypes(self):
        """
        Test that getAtomType() returns appropriate sulfur atom types.
        """
        self.assertEqual(self.atomType(self.mol4, 8), 'Ss')
        self.assertEqual(self.atomType(self.mol4, 9), 'Sd')
    
    def testNoneTypes(self):
        """
        Test that getAtomType() returns appropriate NoneTypes.
        """
        self.assertIsNone(self.atomType(self.mol5, 0))
        self.assertIsNone(self.atomType(self.mol6, 0))
        self.assertIsNone(self.atomType(self.mol7, 0))
        self.assertIsNone(self.atomType(self.mol8, 0))

################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
