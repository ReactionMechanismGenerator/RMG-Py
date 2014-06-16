#!/usr/bin/python
# -*- coding: utf-8 -*-

import unittest
from external.wip import work_in_progress
from rmgpy.molecule.molecule import *
from rmgpy.molecule.group import Group
from rmgpy.molecule.element import getElement, elementList

################################################################################

class TestGroup(unittest.TestCase):
    """
    Contains adjacency list unit tests of the Graph class.
    """

    def setUp(self):
        self.adjlist = """
1 *2 {Cs,Cd} U0 {2,{S,D}} {3,S}
2 *1 {Os,Od} U0 {1,{S,D}}
3    R!H     U0 {1,S}
            """
        self.group = Group().fromAdjacencyList(self.adjlist)


    def testFromAdjacencyList(self):
        """
        adjlist: Test the Group.fromAdjacencyList() method.
        """
        atom1, atom2, atom3 = self.group.atoms
        self.assertTrue(self.group.hasBond(atom1, atom2))
        self.assertTrue(self.group.hasBond(atom1, atom3))
        self.assertFalse(self.group.hasBond(atom2, atom3))
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
        adjlist = self.group.toAdjacencyList()
        print 'adjlist', adjlist
        self.assertEqual(adjlist.strip(), self.adjlist.strip())



class TestMolecule(unittest.TestCase):
    """
    adjlist: Contains adjacency list unit tests of the Molecule class.
    """
    
    def setUp(self):
        self.adjlist_1 = """
1 *1 C U1 L0 E0  {2,S} {3,S} {4,S}
2    H U0 L0 E0  {1,S}
3    H U0 L0 E0  {1,S}
4 *2 N U0 L0 E+1 {1,S} {5,S} {6,D}
5    O U0 L3 E-1 {4,S}
6    O U0 L2 E0  {4,D}
            """
        self.molecule = [Molecule().fromAdjacencyList(self.adjlist_1)]
        
        self.adjlist_2 = """
1 *1 C U1 L0 {2,S} {3,S} {4,S}
2    H U0 L0 {1,S}
3    H U0 L0 {1,S}
4 *2 N U0 L0 {1,S} {5,S} {6,D}
5    O U0 L3 {4,S}
6    O U0 L2 {4,D}
            """
        self.molecule.append(Molecule().fromAdjacencyList(self.adjlist_2))
        
        self.adjlist_3 = """
1 *1 C U1 {2,S} {3,S} {4,S}
2    H U0 {1,S}
3    H U0 {1,S}
4 *2 N U0 {1,S} {5,S} {6,D}
5    O U0 {4,S}
6    O U0 {4,D}
            """
        self.molecule.append(Molecule().fromAdjacencyList(self.adjlist_3))
        
        self.adjlist_4 = """
1 *1 C U1 L0 {2,S}
2 *2 N U0 L0 {1,S} {3,S} {4,D}
3    O U0 L3 {2,S}
4    O U0 L2 {2,D}
            """
        self.molecule.append(Molecule().fromAdjacencyList(self.adjlist_4,saturateH=True))
        


    def testFromAdjacencyList(self):
        """
        adjlist: Test the Molecule.fromAdjacencyList() method.
        """
        
        # molecule 1
        
        self.assertTrue(self.molecule[0].multiplicity == 2)
        
        atom1 = self.molecule[0].atoms[0]
        atom2 = self.molecule[0].atoms[3]
        atom3 = self.molecule[0].atoms[4]
        atom4 = self.molecule[0].atoms[5]
        self.assertTrue(self.molecule[0].hasBond(atom2,atom1))
        self.assertTrue(self.molecule[0].hasBond(atom2,atom3))
        self.assertTrue(self.molecule[0].hasBond(atom2,atom4))
        self.assertFalse(self.molecule[0].hasBond(atom1,atom3))
        self.assertFalse(self.molecule[0].hasBond(atom1,atom4))
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
        
        
        # molecule 2
        
        self.assertTrue(self.molecule[1].multiplicity == 2)
        
        atom1 = self.molecule[1].atoms[0]
        atom2 = self.molecule[1].atoms[3]
        atom3 = self.molecule[1].atoms[4]
        atom4 = self.molecule[1].atoms[5]
        self.assertTrue(self.molecule[1].hasBond(atom2,atom1))
        self.assertTrue(self.molecule[1].hasBond(atom2,atom3))
        self.assertTrue(self.molecule[1].hasBond(atom2,atom4))
        self.assertFalse(self.molecule[1].hasBond(atom1,atom3))
        self.assertFalse(self.molecule[1].hasBond(atom1,atom4))
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
        
        
        # molecule 3
        
        self.assertTrue(self.molecule[2].multiplicity == 2)
        
        atom1 = self.molecule[2].atoms[0]
        atom2 = self.molecule[2].atoms[3]
        atom3 = self.molecule[2].atoms[4]
        atom4 = self.molecule[2].atoms[5]
        self.assertTrue(self.molecule[2].hasBond(atom2,atom1))
        self.assertTrue(self.molecule[2].hasBond(atom2,atom3))
        self.assertTrue(self.molecule[2].hasBond(atom2,atom4))
        self.assertFalse(self.molecule[2].hasBond(atom1,atom3))
        self.assertFalse(self.molecule[2].hasBond(atom1,atom4))
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
        
        
        # molecule 4
        
        self.assertTrue(self.molecule[3].multiplicity == 2)
        
        atom1 = self.molecule[3].atoms[0]
        atom2 = self.molecule[3].atoms[1]
        atom3 = self.molecule[3].atoms[2]
        atom4 = self.molecule[3].atoms[3]
        self.assertTrue(self.molecule[3].hasBond(atom2,atom1))
        self.assertTrue(self.molecule[3].hasBond(atom2,atom3))
        self.assertTrue(self.molecule[3].hasBond(atom2,atom4))
        self.assertFalse(self.molecule[3].hasBond(atom1,atom3))
        self.assertFalse(self.molecule[3].hasBond(atom1,atom4))
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
        adjlist_1 = self.molecule[0].toAdjacencyList(removeH=False)
        newMolecule = Molecule().fromAdjacencyList(adjlist_1)
        self.assertTrue(self.molecule[0].isIsomorphic(newMolecule))
        
        #self.assertEqual(adjlist_1.strip(), self.adjlist_1.strip())
        
#    def testFromOldAdjacencyList(self):
#        """
#        Test we can read things with implicit hydrogens.
#        """
#        adjList = """
#        1 O 0 
#        """ # should be Water
#        molecule = Molecule().fromAdjacencyList(adjList, saturateH=True) # only works with saturateH=True
#        self.assertEqual(molecule.getFormula(),'H2O')



    def testAdjacencyList(self):
        """
        adjlist: Check the adjacency list read/write functions for a full molecule.
        """
        molecule1 = Molecule().fromAdjacencyList("""
        1  C U0 {2,D} {7,S} {8,S}
        2  C U0 {1,D} {3,S} {9,S}
        3  C U0 {2,S} {4,D} {10,S}
        4  C U0 {3,D} {5,S} {11,S}
        5  C U1 {4,S} {6,S} {12,S}
        6  C U0 {5,S} {13,S} {14,S} {15,S}
        7  H U0 {1,S}
        8  H U0 {1,S}
        9  H U0 {2,S}
        10 H U0 {3,S}
        11 H U0 {4,S}
        12 H U0 {5,S}
        13 H U0 {6,S}
        14 H U0 {6,S}
        15 H U0 {6,S}
        """)
        molecule2 = Molecule().fromSMILES('C=CC=C[CH]C')
        self.assertTrue(molecule1.isIsomorphic(molecule2))
        self.assertTrue(molecule2.isIsomorphic(molecule1))

################################################################################

if __name__ == '__main__':
    unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )
