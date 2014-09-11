#!/usr/bin/python
# -*- coding: utf-8 -*-

import unittest
from external.wip import work_in_progress
from rmgpy.molecule.adjlist import *
from rmgpy.molecule.adjlist import InvalidAdjacencyListError
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
1 *2 [Cs,Cd] u0 {2,[S,D]} {3,S}
2 *1 [Os,Od] u0 {1,[S,D]}
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

    def testFromAdjacencyList_multiplicity(self):
        gp = Group().fromAdjacencyList(
        """
        multiplicity [1]
        1 C u0 p0 c0
        """
        )
        self.assertEqual(len(gp.multiplicity), 1)
        self.assertEqual(gp.multiplicity[0], 1)
    
    def testFromAdjacencyList_multiplicity_list(self):
        gp = Group().fromAdjacencyList(
        """
        multiplicity [ 1, 3, 5 ]
        1 C u0 p0 c0
        """
        )
        self.assertEqual(len(gp.multiplicity), 3)
        self.assertEqual(gp.multiplicity[0], 1)
        self.assertEqual(gp.multiplicity[1], 3)
        self.assertEqual(gp.multiplicity[2], 5)
        
    def testToAdjacencyList(self):
        """
        adjlist: Test the Group.toAdjacencyList() method.
        """
        adjlist = """
1 *2 [Cs,Cd] u0 {2,[S,D]} {3,S}
2 *1 [Os,Od] u0 {1,[S,D]}
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
1 *1 C u1 {2,S} {3,S} {4,S}
2    H u0 {1,S}
3    H u0 {1,S}
4 *2 N u0 p0 c+1 {1,S} {5,S} {6,D}
5    O u0 p3 c-1 {4,S}
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
4 *2 N u0 p0 c+1 {1,S} {5,S} {6,D}
5    O u0 p3 c-1 {4,S}
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
        
        
    def testFromAdjacencyList4(self):
        """
        adjlist: Test the Molecule.fromAdjacencyList() method 4.
        """
        # molecule 4
        adjlist = """
1 *1 C u1 {2,S}
2 *2 N u0 p0 c+1 {1,S} {3,S} {4,D}
3    O u0 p3 c-1 {2,S}
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
    
    def testVariousSpinAdjlists(self):
        """
        adjlist: Test that molecules with old or intermediate adjacency list formats containing unusual 
        spin states can get converted to the proper new adjlist format.
        """
        
        adjlist_2S = """
1 C 2S 0 {2,S} {3,S}
2 H 0  0 {1,S}
3 H 0  0 {1,S}
"""
        adjlist_2S_new ="""
1 C u0 p1 c0  {2,S} {3,S}
2 H u0 p0 c0  {1,S}
3 H u0 p0 c0  {1,S}
"""
        mol_2S = Molecule().fromAdjacencyList(adjlist_2S)
        mol_2S_new = Molecule().fromAdjacencyList(adjlist_2S_new)
        self.assertTrue(mol_2S.isIsomorphic(mol_2S_new))
        
        adjlist_2T = """
1 C 2T 0 {2,S} {3,S}
2 H 0  0 {1,S}
3 H 0  0 {1,S}
"""
        adjlist_2T_new ="""
1 C u2 p0 c0  {2,S} {3,S}
2 H u0 p0 c0  {1,S}
3 H u0 p0 c0  {1,S}
"""
        mol_2T = Molecule().fromAdjacencyList(adjlist_2T)
        mol_2T_new = Molecule().fromAdjacencyList(adjlist_2T_new)
        self.assertTrue(mol_2T.isIsomorphic(mol_2T_new))

        adjlist_3D = """
1 C 3D 0 {2,S}
2 H 0  0 {1,S}
"""
        adjlist_3D_new = """
1 C u1 p1 c0  {2,S}
2 H u0 p0 c0  {1,S}
"""
        mol_3D = Molecule().fromAdjacencyList(adjlist_3D)
        mol_3D_new = Molecule().fromAdjacencyList(adjlist_3D_new)
        self.assertTrue(mol_3D.isIsomorphic(mol_3D_new))

        adjlist_3Q = """
1 N 3Q 1
"""
        adjlist_3Q_new = """
1 N u3 p1 c0 
"""
        mol_3Q = Molecule().fromAdjacencyList(adjlist_3Q)
        mol_3Q_new = Molecule().fromAdjacencyList(adjlist_3Q_new)
        self.assertTrue(mol_3Q.isIsomorphic(mol_3Q_new))
        
        adjlist_4S = """
1 C 4S 0
        """
        adjlist_4S_new = """
1 C u0 p2 c0
"""
        mol_4S = Molecule().fromAdjacencyList(adjlist_4S)
        mol_4S_new = Molecule().fromAdjacencyList(adjlist_4S_new)
        self.assertTrue(mol_4S.isIsomorphic(mol_4S_new))
        
        adjlist_4T = """
1 C 4T 0
"""
        adjlist_4T_new = """
1 C u2 p1 c0
"""
        mol_4T = Molecule().fromAdjacencyList(adjlist_4T)
        mol_4T_new = Molecule().fromAdjacencyList(adjlist_4T_new)
        self.assertTrue(mol_4T.isIsomorphic(mol_4T_new))
        
        adjlist_4V = """
1 C 4V 0
"""
        adjlist_4V_new ="""
1 C u4 p0 c0        
"""
        mol_4V = Molecule().fromAdjacencyList(adjlist_4V)
        mol_4V_new = Molecule().fromAdjacencyList(adjlist_4V_new)
        self.assertTrue(mol_4V.isIsomorphic(mol_4V_new))
        
    def testWildcardAdjlists(self):
        """
        adjlist: Test that molecule adjlists containing wildcards raise an InvalidAdjacencyListError.
        """
        # A molecule with a wildcard assignment
        wildcardAdjlist1 = "1 C u1 px c0"
        wildcardAdjlist2 = "1 C ux p2 c0"
        wildcardAdjlist3 = "1 C u1 p2 cx"
        wildcardAdjlist4 = "1 [C,N] u1 p2 c0"

        with self.assertRaises(InvalidAdjacencyListError):
            Molecule().fromAdjacencyList(wildcardAdjlist1)
        with self.assertRaises(InvalidAdjacencyListError):
            Molecule().fromAdjacencyList(wildcardAdjlist2)
        with self.assertRaises(InvalidAdjacencyListError):
            Molecule().fromAdjacencyList(wildcardAdjlist3)
        with self.assertRaises(InvalidAdjacencyListError):
            Molecule().fromAdjacencyList(wildcardAdjlist4)
            
    def testIncorrectAdjlists(self):
        """
        adjlist: Test that improperly formed adjlists raise an InvalidAdjacencyListError.
        """
        # Carbon with 1 radical and 3 lone pairs = 7 total electrons.  Should have -3 charge but doesn't
        adjlist1 = "1 C u1 p3 c0"
        
        with self.assertRaises(InvalidAdjacencyListError):
            Molecule().fromAdjacencyList(adjlist1)
        
    def testHelium(self):
        """
        adjlist: Test that the adjlist reading and writing works with Helium.
        """
        smiles = '[He]'
        inchi = 'InChI=1S/He'
        adjlist = '1 He u0 p1 c0'
        adjlist_old = '1 He 0'
        adjlist_intermediate = '1 He 0 1'
        
        mol_smiles = Molecule().fromSMILES(smiles)
        mol_inchi = Molecule().fromInChI(inchi)
        mol = Molecule().fromAdjacencyList(adjlist)
        mol_old = Molecule().fromAdjacencyList(adjlist_old)
        mol_intermediate = Molecule().fromAdjacencyList(adjlist_intermediate)
        
        # Isomorphic check
        self.assertTrue(mol_smiles.isIsomorphic(mol))
        self.assertTrue(mol_smiles.isIsomorphic(mol_inchi))
        self.assertTrue(mol_smiles.isIsomorphic(mol_old))
        self.assertTrue(mol_smiles.isIsomorphic(mol_intermediate))
        
        # Adjlist check
        self.assertEqual(mol_smiles.toAdjacencyList().strip(), adjlist)
        self.assertEqual(mol_inchi.toAdjacencyList().strip(), adjlist)
        self.assertEqual(mol.toAdjacencyList().strip(), adjlist)
        self.assertEqual(mol_old.toAdjacencyList().strip(), adjlist)
        self.assertEqual(mol_intermediate.toAdjacencyList().strip(), adjlist)
        
        self.assertEqual(mol.toSMILES(),smiles)
        self.assertEqual(mol.toInChI(),'InChI=1S/He')
        
    def testToAdjacencyList(self):
        """
        adjlist: Test the Molecule.toAdjacencyList() method.
        """
        inter_adjlist = """
1 *1 C 1 0  {2,S} {3,S} {4,S}
2    H 0 0  {1,S}
3    H 0 0  {1,S}
4 *2 N 0 0 {1,S} {5,S} {6,D}
5    O 0 3 {4,S}
6    O 0 2 {4,D}
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
        molecule2 = Molecule().fromAdjacencyList(inter_adjlist)
        adjlist_1 = molecule.toAdjacencyList(removeH=False)
        self.assertEqual(adjlist_1,molecule2.toAdjacencyList())
        newMolecule = Molecule().fromAdjacencyList(adjlist_1)
        self.assertTrue(molecule.isIsomorphic(newMolecule))
        
    def testFromIntermediateAdjacencyList1(self):
        """
        Test we can read an intermediate style adjacency list with implicit hydrogens 1
        """
        adjList = """
        1 O 0 2
        """  # should be Water
        molecule = Molecule().fromAdjacencyList(adjList, saturateH=True)  
        self.assertEqual(molecule.getFormula(), 'H2O')

    def testFromOldAdjacencyList1(self):
        """
        Test we can read an old style adjacency list with implicit hydrogens 1
        """
        adjList = """
        1 O 0 
        """  # should be Water
        molecule = Molecule().fromAdjacencyList(adjList)  
        self.assertEqual(molecule.getFormula(), 'H2O')

    def testFromOldAdjacencyList2(self):
        """
        Test we can read an old style adjacency list with implicit hydrogens 2
        """
        adjlist = """
        1 C 2S
        """  
        adjlist_new = """
        1 C u0 p1 c0 {2,S} {3,S}
        2 H u0 p0 c0 {1,S}
        3 H u0 p0 c0 {1,S}
        """
        molecule = Molecule().fromAdjacencyList(adjlist)  
        molecule_new = Molecule().fromAdjacencyList(adjlist_new)
        self.assertTrue(molecule.isIsomorphic(molecule_new))
        
    def testFromOldAdjacencyList3(self):
        """
        Test we can read an old style adjacency list with implicit hydrogens 3
        """
        adjlist = """
        1 C 0
        """  
        adjlist_new = """
        1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
        2 H u0 p0 c0 {1,S}
        3 H u0 p0 c0 {1,S}
        4 H u0 p0 c0 {1,S}
        5 H u0 p0 c0 {1,S}
        """
        molecule = Molecule().fromAdjacencyList(adjlist)  
        molecule_new = Molecule().fromAdjacencyList(adjlist_new)
        self.assertTrue(molecule.isIsomorphic(molecule_new))
        
    def testFromOldAdjacencyList4(self):
        """
        Test we can read an old style adjacency list with implicit hydrogens 4
        """
        adjlist = """
        1 O 2S
        """  
        adjlist_new = """
        1 O u0 p3 c0
        """
        molecule = Molecule().fromAdjacencyList(adjlist)  
        molecule_new = Molecule().fromAdjacencyList(adjlist_new)
        self.assertTrue(molecule.isIsomorphic(molecule_new))
    
    @work_in_progress
    def testFromOldAdjacencyList5(self):
        """
        Test we can read an old style adjacency list with implicit hydrogens 5
        """
        adjlist = """
        1 C 2S {2,T}
        2 O 2S {1,T}
        """  
        adjlist_new = """
        1 C u0 p1 c-1 {2,T}
        2 O u0 p1 c+1 {1,T}
        """
        molecule = Molecule().fromAdjacencyList(adjlist)  
        molecule_new = Molecule().fromAdjacencyList(adjlist_new)
        self.assertTrue(molecule.isIsomorphic(molecule_new))    
        # Currently the fromOldAdjacencyList cannot correctly interpret CO written in this old form
        # (I don't think any adjlists are actually formed this way.)  
        # Currently 'adjlist' will fail when the Molecule is determined to be non-neurtral in net charge.
        
    def testFromOldAdjacencyList6(self):
        """
        Test we can read an old style adjacency list with implicit hydrogens 1
        """
        adjlist = """
        1 C 4T
        """  
        adjlist_new = """
        1 C u2 p1 c0 
        """
        molecule = Molecule().fromAdjacencyList(adjlist)
        molecule_new = Molecule().fromAdjacencyList(adjlist_new)
        self.assertTrue(molecule.isIsomorphic(molecule_new))
        
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
class TestConsistencyChecker(unittest.TestCase):
    def test_check_hund_rule_fail(self):
        with self.assertRaises(InvalidAdjacencyListError):
            Molecule().fromAdjacencyList("""
            multiplicity 1
            1 C u2 p0 c0
            """, saturateH=True)
    def test_check_hund_rule_success(self):
        try:
            Molecule().fromAdjacencyList("""
            multiplicity 3
            1 C u2 p0 c0
            """, saturateH=True)
        except InvalidAdjacencyListError:
            self.fail('InvalidAdjacencyListError thrown unexpectedly!')
    
if __name__ == '__main__':

    unittest.main(testRunner=unittest.TextTestRunner(verbosity=3))
