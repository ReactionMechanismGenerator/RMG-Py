#!/usr/bin/python
# -*- coding: utf-8 -*-

import unittest
from external.wip import work_in_progress
from rmgpy.molecule.molecule import *
from rmgpy.molecule.group import Group
from rmgpy.molecule.element import getElement, elementList

import logging
logging.basicConfig(level=logging.DEBUG)

################################################################################

class TestGroupAdjLists(unittest.TestCase):
    """
    Contains adjacency list unit tests of the Graph class.
    """
    def setUp(self):
        pass

    def testFromOldAdjacencyList1(self):
        """
        adjlist: Test the Group.fromAdjacencyList() method on an old style adjacency list.
        """
        adjlist = """
1 *2 {Cs,Cd} 0 {2,{S,D}} {3,S}
2 *1 {Os,Od}  0   {1,{S,D}}
3    R!H     {0,1} {1,S}
            """
        group = Group().fromAdjacencyList(adjlist)

        atom1, atom2, atom3 = group.atoms
        self.assertTrue(group.hasBond(atom1, atom2))
        self.assertTrue(group.hasBond(atom1, atom3))
        self.assertFalse(group.hasBond(atom2, atom3))
        bond12 = atom1.bonds[atom2]
        bond13 = atom1.bonds[atom3]

        self.assertTrue(atom1.label == '*2')
        self.assertTrue(atom1.atomType[0].label in ['Cs', 'Cd'])
        self.assertTrue(atom1.atomType[1].label in ['Cs', 'Cd'])
        self.assertTrue(atom1.radicalElectrons == [0])

        self.assertTrue(atom2.label == '*1')
        self.assertTrue(atom2.atomType[0].label in ['Os', 'Od'])
        self.assertTrue(atom2.atomType[1].label in ['Os', 'Od'])
        self.assertTrue(atom2.radicalElectrons == [0])

        self.assertTrue(atom3.label == '')
        self.assertTrue(atom3.atomType[0].label == 'R!H')
        self.assertTrue(atom3.radicalElectrons == [0, 1])

        self.assertTrue(bond12.order == ['S', 'D'])
        self.assertTrue(bond13.order == ['S'])


    def testFromAdjacencyList(self):
        """
        adjlist: Test the Group.fromAdjacencyList() method.
        """
        adjlist = """
1 *2 {Cs,Cd} u0 {2,{S,D}} {3,S}
2 *1 {Os,Od} u0 {1,{S,D}}
3    R!H     u0 {1,S}
            """
        group = Group().fromAdjacencyList(adjlist)
        
        atom1, atom2, atom3 = group.atoms
        self.assertTrue(group.hasBond(atom1, atom2))
        self.assertTrue(group.hasBond(atom1, atom3))
        self.assertFalse(group.hasBond(atom2, atom3))
        bond12 = atom1.bonds[atom2]
        bond13 = atom1.bonds[atom3]
           
        self.assertTrue(atom1.label == '*2')
        self.assertTrue(atom1.atomType[0].label in ['Cs', 'Cd'])
        self.assertTrue(atom1.atomType[1].label in ['Cs', 'Cd'])
        self.assertTrue(atom1.radicalElectrons == [0])

        self.assertTrue(atom2.label == '*1')
        self.assertTrue(atom2.atomType[0].label in ['Os', 'Od'])
        self.assertTrue(atom2.atomType[1].label in ['Os', 'Od'])
        self.assertTrue(atom2.radicalElectrons == [0])
        
        self.assertTrue(atom3.label == '')
        self.assertTrue(atom3.atomType[0].label == 'R!H')
        self.assertTrue(atom3.radicalElectrons == [0])

        self.assertTrue(bond12.order == ['S', 'D'])
        self.assertTrue(bond13.order == ['S'])

    def testToAdjacencyList(self):
        """
        adjlist: Test the Group.toAdjacencyList() method.
        """
        adjlist = """
1 *2 {Cs,Cd} u0 {2,{S,D}} {3,S}
2 *1 {Os,Od} u0 {1,{S,D}}
3    R!H     u0 {1,S}
            """
        group = Group().fromAdjacencyList(adjlist)
        adjlist2 = group.toAdjacencyList()

        self.assertEqual(adjlist.strip(), adjlist2.strip())


class TestMoleculeAdjLists(unittest.TestCase):
    """
    adjlist: Contains adjacency list unit tests of the Molecule class.
    """
    
    def setUp(self):
        pass

    def testFromAdjacencyList1(self):
        """
        adjlist: Test the Molecule.fromAdjacencyList() method 1.
        """
        # molecule 1
        adjlist = """
1 *1 C u1 p0 c0  {2,S} {3,S} {4,S}
2    H u0 p0 c0  {1,S}
3    H u0 p0 c0  {1,S}
4 *2 N u0 p0 c+1 {1,S} {5,S} {6,D}
5    O u0 p3 c-1 {4,S}
6    O u0 p2 c0  {4,D}
            """
        molecule = Molecule().fromAdjacencyList(adjlist)
        
        self.assertTrue(molecule.multiplicity == 2)
        
        atom1 = molecule.atoms[0]
        atom2 = molecule.atoms[3]
        atom3 = molecule.atoms[4]
        atom4 = molecule.atoms[5]
        self.assertTrue(molecule.hasBond(atom2, atom1))
        self.assertTrue(molecule.hasBond(atom2, atom3))
        self.assertTrue(molecule.hasBond(atom2, atom4))
        self.assertFalse(molecule.hasBond(atom1, atom3))
        self.assertFalse(molecule.hasBond(atom1, atom4))
        bond21 = atom2.bonds[atom1]
        bond23 = atom2.bonds[atom3]
        bond24 = atom2.bonds[atom4]
           
        self.assertTrue(atom1.label == '*1')
        self.assertTrue(atom1.element.symbol == 'C')
        self.assertTrue(atom1.radicalElectrons == 1)
        self.assertTrue(atom1.charge == 0)
        
        self.assertTrue(atom2.label == '*2')
        self.assertTrue(atom2.element.symbol == 'N')
        self.assertTrue(atom2.radicalElectrons == 0)
        self.assertTrue(atom2.charge == 1)
        
        self.assertTrue(atom3.label == '')
        self.assertTrue(atom3.element.symbol == 'O')
        self.assertTrue(atom3.radicalElectrons == 0)
        self.assertTrue(atom3.charge == -1)
        
        self.assertTrue(atom4.label == '')
        self.assertTrue(atom4.element.symbol == 'O')
        self.assertTrue(atom4.radicalElectrons == 0)
        self.assertTrue(atom4.charge == 0)

        self.assertTrue(bond21.isSingle())
        self.assertTrue(bond23.isSingle())
        self.assertTrue(bond24.isDouble())
        
        
    def testFromAdjacencyList2(self):
        """
        adjlist: Test the Molecule.fromAdjacencyList() method 2.
        """
        # molecule 2
        adjlist = """
1 *1 C u1 p0 {2,S} {3,S} {4,S}
2    H u0 p0 {1,S}
3    H u0 p0 {1,S}
4 *2 N u0 p0 {1,S} {5,S} {6,D}
5    O u0 p3 {4,S}
6    O u0 p2 {4,D}
            """
        molecule = Molecule().fromAdjacencyList(adjlist)
        
        self.assertTrue(molecule.multiplicity == 2)
        
        atom1 = molecule.atoms[0]
        atom2 = molecule.atoms[3]
        atom3 = molecule.atoms[4]
        atom4 = molecule.atoms[5]
        self.assertTrue(molecule.hasBond(atom2, atom1))
        self.assertTrue(molecule.hasBond(atom2, atom3))
        self.assertTrue(molecule.hasBond(atom2, atom4))
        self.assertFalse(molecule.hasBond(atom1, atom3))
        self.assertFalse(molecule.hasBond(atom1, atom4))
        bond21 = atom2.bonds[atom1]
        bond23 = atom2.bonds[atom3]
        bond24 = atom2.bonds[atom4]
           
        self.assertTrue(atom1.label == '*1')
        self.assertTrue(atom1.element.symbol == 'C')
        self.assertTrue(atom1.radicalElectrons == 1)
        self.assertTrue(atom1.charge == 0)
        
        self.assertTrue(atom2.label == '*2')
        self.assertTrue(atom2.element.symbol == 'N')
        self.assertTrue(atom2.radicalElectrons == 0)
        self.assertTrue(atom2.charge == 1)
        
        self.assertTrue(atom3.label == '')
        self.assertTrue(atom3.element.symbol == 'O')
        self.assertTrue(atom3.radicalElectrons == 0)
        self.assertTrue(atom3.charge == -1)
        
        self.assertTrue(atom4.label == '')
        self.assertTrue(atom4.element.symbol == 'O')
        self.assertTrue(atom4.radicalElectrons == 0)
        self.assertTrue(atom4.charge == 0)

        self.assertTrue(bond21.isSingle())
        self.assertTrue(bond23.isSingle())
        self.assertTrue(bond24.isDouble())
        
    def testFromAdjacencyList3(self):
        """
        adjlist: Test the Molecule.fromAdjacencyList() method 3.
        """
        # molecule 3
        adjlist = """
1 *1 C u1 {2,S} {3,S} {4,S}
2    H u0 {1,S}
3    H u0 {1,S}
4 *2 N u0 {1,S} {5,S} {6,D}
5    O u0 {4,S}
6    O u0 {4,D}
            """
        molecule = Molecule().fromAdjacencyList(adjlist)
        
        self.assertTrue(molecule.multiplicity == 2)
        
        atom1 = molecule.atoms[0]
        atom2 = molecule.atoms[3]
        atom3 = molecule.atoms[4]
        atom4 = molecule.atoms[5]
        self.assertTrue(molecule.hasBond(atom2, atom1))
        self.assertTrue(molecule.hasBond(atom2, atom3))
        self.assertTrue(molecule.hasBond(atom2, atom4))
        self.assertFalse(molecule.hasBond(atom1, atom3))
        self.assertFalse(molecule.hasBond(atom1, atom4))
        bond21 = atom2.bonds[atom1]
        bond23 = atom2.bonds[atom3]
        bond24 = atom2.bonds[atom4]
           
        self.assertTrue(atom1.label == '*1')
        self.assertTrue(atom1.element.symbol == 'C')
        self.assertTrue(atom1.radicalElectrons == 1)
        self.assertTrue(atom1.charge == 0)
        
        self.assertTrue(atom2.label == '*2')
        self.assertTrue(atom2.element.symbol == 'N')
        self.assertTrue(atom2.radicalElectrons == 0)
        self.assertTrue(atom2.charge == 1)
        
        self.assertTrue(atom3.label == '')
        self.assertTrue(atom3.element.symbol == 'O')
        self.assertTrue(atom3.radicalElectrons == 0)
        self.assertTrue(atom3.charge == -1)
        
        self.assertTrue(atom4.label == '')
        self.assertTrue(atom4.element.symbol == 'O')
        self.assertTrue(atom4.radicalElectrons == 0)
        self.assertTrue(atom4.charge == 0)

        self.assertTrue(bond21.isSingle())
        self.assertTrue(bond23.isSingle())
        self.assertTrue(bond24.isDouble())
        
        
    def testFromAdjacencyList4(self):
        """
        adjlist: Test the Molecule.fromAdjacencyList() method 4.
        """
        # molecule 4
        adjlist = """
1 *1 C u1 p0 {2,S}
2 *2 N u0 p0 {1,S} {3,S} {4,D}
3    O u0 p3 {2,S}
4    O u0 p2 {2,D}
            """
        molecule = Molecule().fromAdjacencyList(adjlist, saturateH=True)
        
        self.assertTrue(molecule.multiplicity == 2)
        
        atom1 = molecule.atoms[0]
        atom2 = molecule.atoms[1]
        atom3 = molecule.atoms[2]
        atom4 = molecule.atoms[3]
        self.assertTrue(molecule.hasBond(atom2, atom1))
        self.assertTrue(molecule.hasBond(atom2, atom3))
        self.assertTrue(molecule.hasBond(atom2, atom4))
        self.assertFalse(molecule.hasBond(atom1, atom3))
        self.assertFalse(molecule.hasBond(atom1, atom4))
        bond21 = atom2.bonds[atom1]
        bond23 = atom2.bonds[atom3]
        bond24 = atom2.bonds[atom4]
           
        self.assertTrue(atom1.label == '*1')
        self.assertTrue(atom1.element.symbol == 'C')
        self.assertTrue(atom1.radicalElectrons == 1)
        self.assertTrue(atom1.charge == 0)
        
        self.assertTrue(atom2.label == '*2')
        self.assertTrue(atom2.element.symbol == 'N')
        self.assertTrue(atom2.radicalElectrons == 0)
        self.assertTrue(atom2.charge == 1)
        
        self.assertTrue(atom3.label == '')
        self.assertTrue(atom3.element.symbol == 'O')
        self.assertTrue(atom3.radicalElectrons == 0)
        self.assertTrue(atom3.charge == -1)
        
        self.assertTrue(atom4.label == '')
        self.assertTrue(atom4.element.symbol == 'O')
        self.assertTrue(atom4.radicalElectrons == 0)
        self.assertTrue(atom4.charge == 0)

        self.assertTrue(bond21.isSingle())
        self.assertTrue(bond23.isSingle())
        self.assertTrue(bond24.isDouble())


    def testToAdjacencyList(self):
        """
        adjlist: Test the Molecule.toAdjacencyList() method.
        """
        adjlist = """
1 *1 C u1 p0 c0  {2,S} {3,S} {4,S}
2    H u0 p0 c0  {1,S}
3    H u0 p0 c0  {1,S}
4 *2 N u0 p0 c+1 {1,S} {5,S} {6,D}
5    O u0 p3 c-1 {4,S}
6    O u0 p2 c0  {4,D}
            """
        molecule = Molecule().fromAdjacencyList(adjlist)
        adjlist_1 = molecule.toAdjacencyList(removeH=False)
        newMolecule = Molecule().fromAdjacencyList(adjlist_1)
        self.assertTrue(molecule.isIsomorphic(newMolecule))
        
        #self.assertEqual(adjlist_1.strip(), adjlist.strip())
        
    def testFromIntermediateAdjacencyList1(self):
        """
        Test we can read an intermediate style adjacency list with implicit hydrogens 1
        """
        adjList = """
        1 O 0 0
        """  # should be Water
        molecule = Molecule().fromAdjacencyList(adjList, saturateH=True)  # only works with saturateH=True
        self.assertEqual(molecule.getFormula(), 'H2O')

    def testFromOldAdjacencyList1(self):
        """
        Test we can read an old style adjacency list with implicit hydrogens 1
        """
        adjList = """
        1 O 0 
        """  # should be Water
        molecule = Molecule().fromAdjacencyList(adjList)  # only works with saturateH=True
        self.assertEqual(molecule.getFormula(), 'H2O')


    def testAdjacencyList(self):
        """
        adjlist: Check the adjacency list read/write functions for a full molecule.
        """
        molecule1 = Molecule().fromAdjacencyList("""
        1  C u0 {2,D} {7,S} {8,S}
        2  C u0 {1,D} {3,S} {9,S}
        3  C u0 {2,S} {4,D} {10,S}
        4  C u0 {3,D} {5,S} {11,S}
        5  C u1 {4,S} {6,S} {12,S}
        6  C u0 {5,S} {13,S} {14,S} {15,S}
        7  H u0 {1,S}
        8  H u0 {1,S}
        9  H u0 {2,S}
        10 H u0 {3,S}
        11 H u0 {4,S}
        12 H u0 {5,S}
        13 H u0 {6,S}
        14 H u0 {6,S}
        15 H u0 {6,S}
        """)
        molecule2 = Molecule().fromSMILES('C=CC=C[CH]C')
        self.assertTrue(molecule1.isIsomorphic(molecule2))
        self.assertTrue(molecule2.isIsomorphic(molecule1))

################################################################################

if __name__ == '__main__':

    unittest.main(testRunner=unittest.TextTestRunner(verbosity=3))
