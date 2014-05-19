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
        self.mol1 = Molecule().fromSMILES('COC(=O)CC=C=CC#C')
        # self.mol2 = Molecule().fromSMILES('c1ccccc1')
        ## the fromSMILES method currently Kekulizes, so to test Benzene we use fromAdjacencyList
        self.mol2 = Molecule().fromAdjacencyList('''1  C U0 L0 {2,B} {6,B} {7,S}
                                                    2  C U0 L0 {1,B} {3,B} {8,S}
                                                    3  C U0 L0 {2,B} {4,B} {9,S}
                                                    4  C U0 L0 {3,B} {5,B} {10,S}
                                                    5  C U0 L0 {4,B} {6,B} {11,S}
                                                    6  C U0 L0 {1,B} {5,B} {12,S}
                                                    7  H U0 L0 {1,S}
                                                    8  H U0 L0 {2,S}
                                                    9  H U0 L0 {3,S}
                                                    10 H U0 L0 {4,S}
                                                    11 H U0 L0 {5,S}
                                                    12 H U0 L0 {6,S}''')
        self.mol3 = Molecule().fromSMILES('[H]')
        self.mol4 = Molecule().fromSMILES(
                                'O=[Si][Si][Si]=[Si]=[Si][Si]#[Si]SS=S')
        self.mol5 = Molecule().fromAdjacencyList('''1 H U0 L0 {3,S}
                                                    2 H U0 L0 {3,S}
                                                    3 N U0 L0 {1,S} {2,S} {4,D}
                                                    4 N U0 L2 {3,D}''')
        self.mol6 = Molecule().fromSMILES('[Ar]')
        self.mol7 = Molecule().fromSMILES('[He]')
        self.mol8 = Molecule().fromSMILES('[Ne]')
        self.mol9 = Molecule().fromAdjacencyList('''1 N U0 L1 {2,S} {3,S} {4,S}
                                                    2 H U0 L0 {1,S}
                                                    3 H U0 L0 {1,S}
                                                    4 H U0 L0 {1,S}''')
        
        self.mol10 = Molecule().fromAdjacencyList('''1 N U1 L1 {2,S} {3,S}
                                                     2 H U0 L0 {1,S}
                                                     3 H U0 L0 {1,S}''')
        
        self.mol11 = Molecule().fromAdjacencyList('''1 N U2 L1 {2,S}
                                                     2 H U0 L0 {1,S}''')
        
        self.mol12 = Molecule().fromAdjacencyList('''1 N U0 L1 {2,T}
                                                     2 C U1 L0 {1,T}''')
        
        self.mol13 = Molecule().fromAdjacencyList('''1 N U0 L0 {2,S} {3,S} {4,S} {5,S}
                                                     2 H U0 L0 {1,S}
                                                     3 H U0 L0 {1,S}
                                                     4 H U0 L0 {1,S}
                                                     5 O U0 L3 {1,S}''')
        
        self.mol14 = Molecule().fromAdjacencyList('''1 N U0 L2 {2,D}
                                                     2 N U0 L0 {1,D} {3,D}
                                                     3 O U0 L2 {2,D}''')
        
        self.mol15 = Molecule().fromAdjacencyList('''1 N U0 L1 {2,T}
                                                     2 N U0 L0 {1,T} {3,S}
                                                     3 O U0 L3 {2,S}''')
        
        self.mol16 = Molecule().fromAdjacencyList('''1 N U0 L1 {2,D} {3,S}
                                                     2 O U0 L2 {1,D}
                                                     3 O U1 L2 {1,S}''')
        
        self.mol17 = Molecule().fromAdjacencyList('''1 N U1 L1 {2,D}
                                                     2 O U0 L2 {1,D}''')
        
        self.mol18 = Molecule().fromAdjacencyList('''1  N U0 L0 {2,B} {6,B} {7,S}
                                                     2  C U0 L0 {1,B} {3,B} {8,S}
                                                     3  C U0 L0 {2,B} {4,B} {9,S}
                                                     4  C U0 L0 {3,B} {5,B} {10,S}
                                                     5  C U0 L0 {4,B} {6,B} {11,S}
                                                     6  N U0 L1 {1,B} {5,B}
                                                     7  O U0 L3 {1,S}
                                                     8  H U0 L0 {2,S}
                                                     9  H U0 L0 {3,S}
                                                     10 H U0 L0 {4,S}
                                                     11 H U0 L0 {5,S}''')
        
    
    def atomType(self, mol, atomID):
        atom = mol.atoms[atomID]
        type = getAtomType(atom, mol.getBonds(atom))
        if type is None:
            return type
        else:
            return type.label

    def testHydrogenType(self):
        """
        Test that getAtomType() returns the hydrogen atom type.
        """
        self.assertEqual(self.atomType(self.mol3, 0), 'H')
        
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
    
    def testNitrogenTypes(self):
        """
        Test that getAtomType() returns appropriate nitrogen atom types.
        """
        self.assertEqual(self.atomType(self.mol5, 2), 'N5d')
        self.assertEqual(self.atomType(self.mol5, 3), 'N1d')
        self.assertEqual(self.atomType(self.mol9, 0), 'N3s')
        self.assertEqual(self.atomType(self.mol10, 0), 'N3s')
        self.assertEqual(self.atomType(self.mol11, 0), 'N3s')
        self.assertEqual(self.atomType(self.mol16, 0), 'N3d')
        self.assertEqual(self.atomType(self.mol17, 0), 'N3d')
        self.assertEqual(self.atomType(self.mol12, 0), 'N3t')
        self.assertEqual(self.atomType(self.mol13, 0), 'N5s')
        self.assertEqual(self.atomType(self.mol14, 1), 'N5dd')
        self.assertEqual(self.atomType(self.mol15, 1), 'N5t')
        self.assertEqual(self.atomType(self.mol18, 5), 'N3b')
        self.assertEqual(self.atomType(self.mol18, 0), 'N5b')
        
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
    
    def testOtherTypes(self):
        """
        Test that getAtomType() returns appropriate types for other misc inerts.
        """
        self.assertEqual(self.atomType(self.mol6, 0), 'Ar')
        self.assertEqual(self.atomType(self.mol7, 0), 'He')
        self.assertEqual(self.atomType(self.mol8, 0), 'Ne')

################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
