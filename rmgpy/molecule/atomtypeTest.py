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
        self.mol2 = Molecule().fromAdjacencyList('''1  C u0 p0 {2,B} {6,B} {7,S}
                                                    2  C u0 p0 {1,B} {3,B} {8,S}
                                                    3  C u0 p0 {2,B} {4,B} {9,S}
                                                    4  C u0 p0 {3,B} {5,B} {10,S}
                                                    5  C u0 p0 {4,B} {6,B} {11,S}
                                                    6  C u0 p0 {1,B} {5,B} {12,S}
                                                    7  H u0 p0 {1,S}
                                                    8  H u0 p0 {2,S}
                                                    9  H u0 p0 {3,S}
                                                    10 H u0 p0 {4,S}
                                                    11 H u0 p0 {5,S}
                                                    12 H u0 p0 {6,S}''')
        self.mol3 = Molecule().fromSMILES('[H]')
        self.mol4 = Molecule().fromSMILES(
                                'O=[Si][Si][Si]=[Si]=[Si][Si]#[Si]SS=S')
        self.mol5 = Molecule().fromAdjacencyList('''1 H u0 p0 {3,S}
                                                    2 H u0 p0 {3,S}
                                                    3 N u0 p0 c+1 {1,S} {2,S} {4,D}
                                                    4 N u0 p2 c-1 {3,D}''')
        self.mol6 = Molecule().fromSMILES('[Ar]')
        self.mol7 = Molecule().fromSMILES('[He]')
        self.mol8 = Molecule().fromSMILES('[Ne]')
        self.mol9 = Molecule().fromAdjacencyList('''1 N u0 p1 {2,S} {3,S} {4,S}
                                                    2 H u0 p0 {1,S}
                                                    3 H u0 p0 {1,S}
                                                    4 H u0 p0 {1,S}''')
        
        self.mol10 = Molecule().fromAdjacencyList('''1 N u1 p1 {2,S} {3,S}
                                                     2 H u0 p0 {1,S}
                                                     3 H u0 p0 {1,S}''')
        
        self.mol11 = Molecule().fromAdjacencyList('''1 N u2 p1 {2,S}
                                                     2 H u0 p0 {1,S}''')
        
        self.mol12 = Molecule().fromAdjacencyList('''1 N u0 p1 {2,T}
                                                     2 C u1 p0 {1,T}''')
        
        self.mol14 = Molecule().fromAdjacencyList('''1 N u0 p2 c-1 {2,D}
                                                     2 N u0 p0 c+1 {1,D} {3,D}
                                                     3 O u0 p2 {2,D}''')
        
        self.mol15 = Molecule().fromAdjacencyList('''1 N u0 p1 c0 {2,T}
                                                     2 N u0 p0 c+1 {1,T} {3,S}
                                                     3 O u0 p3 c-1 {2,S}''')
        
        self.mol16 = Molecule().fromAdjacencyList('''1 N u0 p1 {2,D} {3,S}
                                                     2 O u0 p2 {1,D}
                                                     3 O u1 p2 {1,S}''')
        
        self.mol17 = Molecule().fromAdjacencyList('''1 N u1 p1 {2,D}
                                                     2 O u0 p2 {1,D}''')
        
        self.mol18 = Molecule().fromAdjacencyList('''1  N u0 p0 c+1 {2,B} {6,B} {7,S}
                                                     2  C u0 p0 {1,B} {3,B} {8,S}
                                                     3  C u0 p0 {2,B} {4,B} {9,S}
                                                     4  C u0 p0 {3,B} {5,B} {10,S}
                                                     5  C u0 p0 {4,B} {6,B} {11,S}
                                                     6  N u0 p1 {1,B} {5,B}
                                                     7  O u0 p3 c-1 {1,S}
                                                     8  H u0 p0 {2,S}
                                                     9  H u0 p0 {3,S}
                                                     10 H u0 p0 {4,S}
                                                     11 H u0 p0 {5,S}''')

        self.mol19 = Molecule().fromSMILES('C=S')

        self.mol20 = Molecule().fromSMILES('[C-]#[O+]')

        self.mol21 = Molecule().fromAdjacencyList('''1 S u0 p3 c-1 {2,S}
                                                     2 S u0 p2 c+1 {1,S}''')

        self.mol22 = Molecule().fromAdjacencyList('''1 S u0 p3 c0''')

        self.mol23 = Molecule().fromAdjacencyList('''1 S u0 p2 c0 {2,S} {5,S}
                                                     2 S u0 p1 c+1 {1,S} {3,S} {4,S}
                                                     3 C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
                                                     4 O u0 p3 c-1 {2,S}
                                                     5 H u0 p0 c0 {1,S}
                                                     6 H u0 p0 c0 {3,S}
                                                     7 H u0 p0 c0 {3,S}
                                                     8 H u0 p0 c0 {3,S}''')

        self.mol24 = Molecule().fromAdjacencyList('''1 C u0 p0 c0 {2,D} {4,S} {5,S}
                                                     2 S u0 p2 c-1 {1,D} {3,S}
                                                     3 O u0 p2 c+1 {2,S}
                                                     4 H u0 p0 c0 {1,S}
                                                     5 H u0 p0 c0 {1,S}''')

        self.mol25 = Molecule().fromAdjacencyList('''1 S u0 p1 c0 {2,S} {5,S} {7,S} {8,S}
                                                     2 O u0 p2 c0 {1,S} {3,S}
                                                     3 S u0 p1 c0 {2,S} {4,S} {9,D}
                                                     4 O u0 p2 c0 {3,S} {6,S}
                                                     5 H u0 p0 c0 {1,S}
                                                     6 H u0 p0 c0 {4,S}
                                                     7 H u0 p0 c0 {1,S}
                                                     8 H u0 p0 c0 {1,S}
                                                     9 O u0 p2 c0 {3,D}''')

        self.mol26 = Molecule().fromAdjacencyList('''1 O u0 p3 c-1 {2,S}
                                                     2 S u0 p1 c+1 {1,S} {3,D}
                                                     3 O u0 p2 c0 {2,D}''')

        #self.mol27 = Molecule().fromAdjacencyList('''1 S u0 p1 c0 {2,B} {5,B}
        #                                             2 C u0 p0 c0 {1,B} {3,B} {6,S}
        #                                             3 C u0 p0 c0 {2,B} {4,B} {7,S}
        #                                             4 C u0 p0 c0 {3,B} {5,B} {8,S}
        #                                             5 C u0 p0 c0 {1,B} {4,B} {9,S}
        #                                             6 H u0 p0 c0 {2,S}
        #                                             7 H u0 p0 c0 {3,S}
        #                                             8 H u0 p0 c0 {4,S}
        #                                             9 H u0 p0 c0 {5,S}''')

        self.mol28 = Molecule().fromAdjacencyList('''1  O u0 p2 c0 {2,D}
                                                     2  S u0 p1 c0 {1,D} {3,D}
                                                     3  C u0 p0 c0 {2,D} {4,S} {7,S}
                                                     4  C u0 p0 c0 {3,S} {5,T}
                                                     5  S u0 p1 c0 {4,T} {6,S}
                                                     6  S u0 p0 c0 {5,S} {8,S} {9,S} {10,S} {11,S} {12,S}
                                                     7  H u0 p0 c0 {3,S}
                                                     8  H u0 p0 c0 {6,S}
                                                     9  H u0 p0 c0 {6,S}
                                                     10 H u0 p0 c0 {6,S}
                                                     11 H u0 p0 c0 {6,S}
                                                     12 H u0 p0 c0 {6,S}''')

        self.mol29 = Molecule().fromAdjacencyList('''1 C u0 p1 c-1 {2,T}
                                                     2 S u0 p1 c+1 {1,T}''')

        self.mol30 = Molecule().fromAdjacencyList('''1 S u0 p0 c0 {2,D} {3,S} {4,S} {5,S} {6,S}
                                                     2 O u0 p2 c0 {1,D}
                                                     3 H u0 p0 c0 {1,S}
                                                     4 H u0 p0 c0 {1,S}
                                                     5 H u0 p0 c0 {1,S}
                                                     6 H u0 p0 c0 {1,S}''')

        self.mol31 = Molecule().fromAdjacencyList('''1 S u0 p0 c+1 {2,S} {3,D} {4,D}
                                                     2 O u0 p3 c-1 {1,S}
                                                     3 O u0 p2 c0 {1,D}
                                                     4 O u0 p2 c0 {1,D}''')

        self.mol32 = Molecule().fromAdjacencyList('''1 O u0 p2 c0 {2,D}
                                                     2 S u0 p0 c0 {1,D} {3,D} {4,S} {5,S}
                                                     3 O u0 p2 c0 {2,D}
                                                     4 O u0 p2 c0 {2,S} {6,S}
                                                     5 O u0 p2 c0 {2,S} {7,S}
                                                     6 H u0 p0 c0 {4,S}
                                                     7 H u0 p0 c0 {5,S}''')

        self.mol33 = Molecule().fromAdjacencyList('''1 O u0 p3 c-1 {2,S}
                                                     2 S u0 p0 c+1 {1,S} {3,D} {4,D}
                                                     3 O u0 p2 c0 {2,D}
                                                     4 O u0 p2 c0 {2,D}''')

        self.mol34 = Molecule().fromAdjacencyList('''1 O u0 p2 c0 {2,D}
                                                     2 S u0 p0 c0 {1,D} {3,D} {4,D}
                                                     3 O u0 p2 c0 {2,D}
                                                     4 O u0 p2 c0 {2,D}''')

        self.mol35 = Molecule().fromAdjacencyList('''1 S u0 p0 c0 {2,T} {3,S} {4,S} {5,S}
                                                     2 N u0 p1 c0 {1,T}
                                                     3 H u0 p0 c0 {1,S}
                                                     4 H u0 p0 c0 {1,S}
                                                     5 H u0 p0 c0 {1,S}''')

        self.mol36 = Molecule().fromAdjacencyList('''1 S u0 p0 c0 {2,T} {3,D} {4,S}
                                                     2 N u0 p1 c0 {1,T}
                                                     3 O u0 p2 c0 {1,D}
                                                     4 H u0 p0 c0 {1,S}''')

        self.mol37 = Molecule().fromAdjacencyList('''1 N u0 p1 c0 {2,T}
                                                     2 S u0 p0 c0 {1,T} {3,T}
                                                     3 N u0 p1 c0 {2,T}''')

        self.mol38 = Molecule().fromSMILES('O=S=O')

        self.mol39 = Molecule().fromAdjacencyList('''1 N u0 p2 c-1 {2,S} {3,S}
                                                     2 H u0 p0 c0 {1,S}
                                                     3 N u0 p0 c+1 {1,S} {4,T}
                                                     4 C u0 p0 c0 {3,T} {5,S}
                                                     5 H u0 p0 c0 {4,S}''')

        self.mol40 = Molecule().fromAdjacencyList('''1 N u0 p0 c+1 {2,S} {3,T}
                                                     2 H u0 p0 c0 {1,S}
                                                     3 N u0 p0 c+1 {1,T} {4,S}
                                                     4 N u0 p3 c-2 {3,S}''')


        self.mol41 = Molecule().fromAdjacencyList('''1 N u0 p2 c0 {2,S}
                                                     2 H u0 p0 c0 {1,S}''')


        self.mol42 = Molecule().fromAdjacencyList('''1 N u0 p1 c0 {2,T}
                                                     2 N u0 p0 c+1 {1,T} {3,S}
                                                     3 S u0 p2 c-1 {2,S} {4,S} {5,S}
                                                     4 O u1 p2 c0 {3,S}
                                                     5 O u1 p2 c0 {3,S}''')

        self.mol43 = Molecule().fromAdjacencyList('''1 C u0 p1 c-1 {2,D} {3,S}
                                                     2 S u1 p0 c+1 {1,D} {4,S} {5,S}
                                                     3 H u0 p0 c0 {1,S}
                                                     4 H u0 p0 c0 {2,S}
                                                     5 H u0 p0 c0 {2,S}''')

        self.mol44 = Molecule().fromAdjacencyList('''1 O u0 p3 c0''')

        self.mol45 = Molecule().fromAdjacencyList('''1 O u0 p2 c0 {2,S} {5,S}
                                                     2 N u0 p0 c+1 {1,S} {3,S} {4,D}
                                                     3 O u0 p3 c-1 {2,S}
                                                     4 O u0 p2 c0 {2,D}
                                                     5 H u0 p0 c0 {1,S}''')

        self.mol49 = Molecule().fromAdjacencyList('''1 O u0 p3 c-1 {2,S}
                                                     2 O u0 p1 c+1 {1,S} {3,S} {4,S}
                                                     3 H u0 p0 c0 {2,S}
                                                     4 S u0 p2 c0 {2,S} {5,S}
                                                     5 H u0 p0 c0 {4,S}''')

        self.mol50 = Molecule().fromAdjacencyList('''1 O u0 p3 c-1 {2,S}
                                                     2 O u0 p1 c+1 {1,S} {3,D}
                                                     3 C u0 p0 c0 {2,D} {4,S} {5,S}
                                                     4 H u0 p0 c0 {3,S}
                                                     5 H u0 p0 c0 {3,S}''')

        self.mol51 = Molecule().fromAdjacencyList('''1 O u0 p2 c0 {2,S} {7,S}
                                                     2 S u0 p0 c+1 {1,S} {3,S} {4,S} {5,S} {6,S}
                                                     3 H u0 p0 c0 {2,S}
                                                     4 H u0 p0 c0 {2,S}
                                                     5 H u0 p0 c0 {2,S}
                                                     6 O u0 p3 c-1 {2,S}
                                                     7 H u0 p0 c0 {1,S}''')

        self.mol52 = Molecule().fromAdjacencyList('''1  C u0 p0 c0 {2,D} {6,S} {8,S}
                                                     2  C u0 p0 c0 {1,D} {3,S} {9,S}
                                                     3  C u0 p0 c0 {2,S} {4,S} {10,S} {11,S}
                                                     4  C u0 p0 c0 {3,S} {5,S} {6,S} {12,S}
                                                     5  O u0 p3 c-1 {4,S}
                                                     6  C u0 p0 c+1 {1,S} {4,S} {7,S}
                                                     7  H u0 p0 c0 {6,S}
                                                     8  H u0 p0 c0 {1,S}
                                                     9  H u0 p0 c0 {2,S}
                                                     10 H u0 p0 c0 {3,S}
                                                     11 H u0 p0 c0 {3,S}
                                                     12 H u0 p0 c0 {4,S}''')

        self.mol53 = Molecule().fromAdjacencyList('''1 N u0 p0 c-1 {2,D} {3,D} {4,D}
                                                     2 C u0 p0 c0 {1,D} {5,S} {6,S}
                                                     3 C u0 p0 c0 {1,D} {7,S} {8,S}
                                                     4 N u0 p0 c+1 {1,D} {9,S} {10,S}
                                                     5 H u0 p0 c0 {2,S}
                                                     6 H u0 p0 c0 {2,S}
                                                     7 H u0 p0 c0 {3,S}
                                                     8 H u0 p0 c0 {3,S}
                                                     9 H u0 p0 c0 {4,S}
                                                     10 H u0 p0 c0 {4,S}''')

        self.mol54 = Molecule().fromAdjacencyList('''1 C u0 p0 c+1 {2,S} {3,D}
                                                     2 H u0 p0 c0 {1,S}
                                                     3 C u0 p0 c0 {1,D} {4,D}
                                                     4 C u0 p1 c-1 {3,D} {5,S}
                                                     5 H u0 p0 c0 {4,S}''')

        self.mol55 = Molecule().fromAdjacencyList('''1  C u0 p0 c0 {2,B} {10,B} {11,S}
                                                     2  C u0 p0 c0 {1,B} {3,B} {12,S}
                                                     3  C u0 p0 c0 {2,B} {4,B} {13,S}
                                                     4  C u0 p0 c0 {3,B} {5,B} {9,B}
                                                     5  C u0 p0 c0 {4,B} {6,B} {14,S}
                                                     6  C u0 p0 c0 {5,B} {7,B} {15,S}
                                                     7  C u0 p0 c0 {6,B} {8,B} {16,S}
                                                     8  C u0 p0 c0 {7,B} {9,B} {17,S}
                                                     9  C u0 p0 c0 {4,B} {8,B} {10,B}
                                                     10 C u0 p0 c0 {1,B} {9,B} {18,S}
                                                     11 H u0 p0 c0 {1,S}
                                                     12 H u0 p0 c0 {2,S}
                                                     13 H u0 p0 c0 {3,S}
                                                     14 H u0 p0 c0 {5,S}
                                                     15 H u0 p0 c0 {6,S}
                                                     16 H u0 p0 c0 {7,S}
                                                     17 H u0 p0 c0 {8,S}
                                                     18 H u0 p0 c0 {10,S}''')

        self.mol56 = Molecule().fromAdjacencyList('''1 C u0 p1 c0 {2,S} {3,S}
                                                     2 H u0 p0 c0 {1,S}
                                                     3 H u0 p0 c0 {1,S}''')

        self.mol57 = Molecule().fromAdjacencyList('''1 C u0 p1 c-1 {2,S} {3,S} {4,S}
                                                     2 H u0 p0 c0 {1,S}
                                                     3 H u0 p0 c0 {1,S}
                                                     4 N u0 p0 c+1 {1,S} {5,T}
                                                     5 N u0 p1 c0 {4,T}''')

        self.mol58 = Molecule().fromAdjacencyList('''1 C u0 p1 c0 {2,D}
                                                     2 C u0 p0 c0 {1,D} {3,S} {4,S}
                                                     3 H u0 p0 c0 {2,S}
                                                     4 H u0 p0 c0 {2,S}''')

        self.mol59 = Molecule().fromAdjacencyList('''1 C u0 p1 c-1 {2,S} {3,D}
                                                     2 H u0 p0 c0 {1,S}
                                                     3 N u0 p0 c+1 {1,D} {4,D}
                                                     4 O u0 p2 c0 {3,D}''')

        self.mol60 = Molecule().fromAdjacencyList('''1 C u0 p0 c0 {2,D} {3,D}
                                                     2 C u0 p0 c+1 {1,D} {4,S}
                                                     3 C u0 p1 c-1 {1,D} {5,S}
                                                     4 H u0 p0 c0 {2,S}
                                                     5 H u0 p0 c0 {3,S}''')

        self.mol64 = Molecule().fromAdjacencyList('''1 N u0 p1 c0 {2,D} {4,S}
                                                     2 N u0 p0 c+1 {1,D} {3,D}
                                                     3 N u0 p2 c-1 {2,D}
                                                     4 H u0 p0 c0 {1,S}''')

        self.mol65 = Molecule().fromAdjacencyList('''1 N u0 p0 c0 {2,T} {3,S} {4,S}
                                                     2 C u0 p0 c0 {1,T} {5,S}
                                                     3 H u0 p0 c0 {1,S}
                                                     4 H u0 p0 c0 {1,S}
                                                     5 H u0 p0 c0 {2,S}''')

        self.mol69 = Molecule().fromAdjacencyList('''1 N u0 p0 c+1 {2,T} {3,S}
                                                     2 S u0 p2 c-1 {1,T}
                                                     3 H u0 p0 c0 {1,S}''')

        self.mol70 = Molecule().fromAdjacencyList('''1 S u0 p0 c+1 {2,D} {3,T}
                                                     2 N u0 p2 c-1 {1,D}
                                                     3 N u0 p1 c0 {1,T}''')

        #self.mol71 = Molecule().fromAdjacencyList('''1 O u0 p1 c0 {2,B} {5,B}
        #                                             2 C u0 p0 c0 {1,B} {3,B} {6,S}
        #                                             3 C u0 p0 c0 {2,B} {4,B} {7,S}
        #                                             4 C u0 p0 c0 {3,B} {5,B} {8,S}
        #                                             5 C u0 p0 c0 {1,B} {4,B} {9,S}
        #                                             6 H u0 p0 c0 {2,S}
        #                                             7 H u0 p0 c0 {3,S}
        #                                             8 H u0 p0 c0 {4,S}
        #                                             9 H u0 p0 c0 {5,S}''')

        #self.mol72 = Molecule().fromAdjacencyList('''1  N u0 p0 c0 {2,B} {6,B} {7,D}
        #                                             2  C u0 p0 {1,B} {3,B} {8,S}
        #                                             3  C u0 p0 {2,B} {4,B} {9,S}
        #                                             4  C u0 p0 {3,B} {5,B} {10,S}
        #                                             5  C u0 p0 {4,B} {6,B} {11,S}
        #                                             6  N u0 p1 {1,B} {5,B}
        #                                             7  O u0 p2 c0 {1,D}
        #                                             8  H u0 p0 {2,S}
        #                                             9  H u0 p0 {3,S}
        #                                             10 H u0 p0 {4,S}
        #                                             11 H u0 p0 {5,S}''')

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
        self.assertEqual(self.atomType(self.mol52, 5), 'Csc')
        self.assertEqual(self.atomType(self.mol1, 5), 'Cd')
        self.assertEqual(self.atomType(self.mol60, 1), 'Cdc')
        self.assertEqual(self.atomType(self.mol1, 2), 'CO')
        self.assertEqual(self.atomType(self.mol19, 0), 'CS')
        self.assertEqual(self.atomType(self.mol1, 6), 'Cdd')
        self.assertEqual(self.atomType(self.mol1, 9), 'Ct')
        self.assertEqual(self.atomType(self.mol2, 0), 'Cb')
        self.assertEqual(self.atomType(self.mol55, 3), 'Cbf')
        self.assertEqual(self.atomType(self.mol56, 0), 'C2s')
        self.assertEqual(self.atomType(self.mol57, 0), 'C2sc')
        self.assertEqual(self.atomType(self.mol58, 0), 'C2d')
        self.assertEqual(self.atomType(self.mol59, 0), 'C2dc')
        self.assertEqual(self.atomType(self.mol60, 2), 'C2dc')
        self.assertEqual(self.atomType(self.mol20, 0), 'C2tc')
        self.assertEqual(self.atomType(self.mol29, 0), 'C2tc')
    
    def testNitrogenTypes(self):
        """
        Test that getAtomType() returns appropriate nitrogen atom types.
        """
        self.assertEqual(self.atomType(self.mol40, 3), 'N0sc')
        self.assertEqual(self.atomType(self.mol41, 0), 'N1s')
        self.assertEqual(self.atomType(self.mol39, 0), 'N1sc')
        self.assertEqual(self.atomType(self.mol5, 3), 'N1dc')
        self.assertEqual(self.atomType(self.mol9, 0), 'N3s')
        self.assertEqual(self.atomType(self.mol10, 0), 'N3s')
        self.assertEqual(self.atomType(self.mol11, 0), 'N3s')
        self.assertEqual(self.atomType(self.mol16, 0), 'N3d')
        self.assertEqual(self.atomType(self.mol17, 0), 'N3d')
        self.assertEqual(self.atomType(self.mol12, 0), 'N3t')
        self.assertEqual(self.atomType(self.mol18, 5), 'N3b')
        self.assertEqual(self.atomType(self.mol5, 2), 'N5dc')
        self.assertEqual(self.atomType(self.mol64, 1), 'N5ddc')
        self.assertEqual(self.atomType(self.mol53, 0), 'N5dddc')
        self.assertEqual(self.atomType(self.mol65, 0), 'N5t')
        self.assertEqual(self.atomType(self.mol15, 1), 'N5tc')
        self.assertEqual(self.atomType(self.mol39, 2), 'N5tc')
        self.assertEqual(self.atomType(self.mol18, 0), 'N5b')
        # self.assertEqual(self.atomType(self.mol72, 0), 'N5bd')  # aromatic nitrogen currently doesn't work well in RMG. See RMG-Py #982
        
    def testOxygenTypes(self):
        """
        Test that getAtomType() returns appropriate oxygen atom types.
        """
        self.assertEqual(self.atomType(self.mol44, 0), 'Oa')
        self.assertEqual(self.atomType(self.mol45, 2), 'O0sc')
        self.assertEqual(self.atomType(self.mol49, 0), 'O0sc')
        self.assertEqual(self.atomType(self.mol1, 1),  'O2s')
        self.assertEqual(self.atomType(self.mol24, 2), 'O2sc')
        self.assertEqual(self.atomType(self.mol1, 3),  'O2d')
        self.assertEqual(self.atomType(self.mol49, 1), 'O4sc')
        self.assertEqual(self.atomType(self.mol50, 1), 'O4dc')
        self.assertEqual(self.atomType(self.mol20, 1), 'O4tc')
        # self.assertEqual(self.atomType(self.mol71, 0), 'O4b')  # aromatic oxygen currently doesn't work well in RMG. See RMG-Py #982
    
    def testSiliconTypes(self):
        """
        Test that getAtomType() returns appropriate silicon atom types.
        """
        self.assertEqual(self.atomType(self.mol4, 2), 'Sis')
        self.assertEqual(self.atomType(self.mol4, 1), 'SiO')
        self.assertEqual(self.atomType(self.mol4, 5), 'Sid')
        self.assertEqual(self.atomType(self.mol4, 4), 'Sidd')
        self.assertEqual(self.atomType(self.mol4, 7), 'Sit')
    
    def testSulfurTypes(self):
        """
        Test that getAtomType() returns appropriate sulfur atom types.
        """
        self.assertEqual(self.atomType(self.mol22, 0), 'Sa')
        self.assertEqual(self.atomType(self.mol21, 0), 'S0sc')
        self.assertEqual(self.atomType(self.mol23, 0), 'S2s')
        self.assertEqual(self.atomType(self.mol21, 1), 'S2sc')
        self.assertEqual(self.atomType(self.mol42, 2), 'S2sc')
        self.assertEqual(self.atomType(self.mol19, 1), 'S2d')
        self.assertEqual(self.atomType(self.mol24, 1), 'S2dc')
        self.assertEqual(self.atomType(self.mol69, 1), 'S2tc')
        self.assertEqual(self.atomType(self.mol25, 0), 'S4s')
        self.assertEqual(self.atomType(self.mol23, 1), 'S4sc')
        self.assertEqual(self.atomType(self.mol25, 2), 'S4d')
        self.assertEqual(self.atomType(self.mol28, 1), 'S4dd')
        self.assertEqual(self.atomType(self.mol38, 1), 'S4dd')
        self.assertEqual(self.atomType(self.mol26, 1), 'S4dc')
        # self.assertEqual(self.atomType(self.mol27, 0), 'S4b')  # aromatic sulfur currently doesn't work well in RMG. See RMG-Py #982
        self.assertEqual(self.atomType(self.mol28, 4), 'S4t')
        self.assertEqual(self.atomType(self.mol29, 1), 'S4tdc')
        self.assertEqual(self.atomType(self.mol28, 5), 'S6s')
        self.assertEqual(self.atomType(self.mol51, 1), 'S6sc')
        self.assertEqual(self.atomType(self.mol30, 0), 'S6d')
        self.assertEqual(self.atomType(self.mol32, 1), 'S6dd')
        self.assertEqual(self.atomType(self.mol34, 1), 'S6ddd')
        self.assertEqual(self.atomType(self.mol43, 1), 'S6dc')
        self.assertEqual(self.atomType(self.mol31, 0), 'S6dc')
        self.assertEqual(self.atomType(self.mol33, 1), 'S6dc')
        self.assertEqual(self.atomType(self.mol35, 0), 'S6t')
        self.assertEqual(self.atomType(self.mol36, 0), 'S6td')
        self.assertEqual(self.atomType(self.mol37, 1), 'S6tt')
        self.assertEqual(self.atomType(self.mol70, 0), 'S6tdc')
    
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
