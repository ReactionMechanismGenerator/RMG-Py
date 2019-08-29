import os
import unittest

from rmgpy.species import Species
from rmgpy.molecule import resonance
from rmgpy.molecule.atomtype import atomTypes
from rmgpy.molecule.element import getElement
from rmgpy.molecule.molecule import Atom, Bond, Molecule

import afm.fragment

class TestCuttingLabel(unittest.TestCase):

    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.cutting_label_R = afm.fragment.CuttingLabel('R')

    def test_symbol(self):

        self.assertEqual('R', self.cutting_label_R.symbol)

    def test_copy(self):

        cutting_label_R_copy = self.cutting_label_R.copy()

        self.assertEqual('R', cutting_label_R_copy.name)
        self.assertEqual(self.cutting_label_R.label, 
                         cutting_label_R_copy.label)
        self.assertEqual(self.cutting_label_R.charge, 
                         cutting_label_R_copy.charge)
        self.assertEqual(self.cutting_label_R.radicalElectrons, 
                         cutting_label_R_copy.radicalElectrons)
        self.assertEqual(self.cutting_label_R.lonePairs, 
                         cutting_label_R_copy.lonePairs)
        self.assertEqual(self.cutting_label_R.isotope, 
                         cutting_label_R_copy.isotope)

class TestFragment(unittest.TestCase):

    def setUp(self):
        """
        A function run before each unit test in this class.
        """

        # construct the first fragment
        atom_C1 = Atom(element=getElement('C'), 
                    radicalElectrons=0, 
                    charge=0, 
                    lonePairs=0)

        cutting_label_R1 = afm.fragment.CuttingLabel('R')
        cutting_label_L1 = afm.fragment.CuttingLabel('L')
        
        vertices = [
            atom_C1,
            cutting_label_R1,
            cutting_label_L1
        ]

        bonds = [
            Bond(atom_C1, cutting_label_R1),
            Bond(atom_C1, cutting_label_L1)
        ]
        
        self.fragment1 = afm.fragment.Fragment()
        for vertex in vertices: self.fragment1.addVertex(vertex)
        for bond in bonds: self.fragment1.addEdge(bond)

        # construct the second fragment
        atom_C2 = Atom(element=getElement('C'), 
                    radicalElectrons=0, 
                    charge=0, 
                    lonePairs=0)

        cutting_label_R2 = afm.fragment.CuttingLabel('R')
        cutting_label_L2 = afm.fragment.CuttingLabel('L')
        
        vertices = [
            atom_C2,
            cutting_label_R2,
            cutting_label_L2
        ]

        bonds = [
            Bond(atom_C2, cutting_label_R2),
            Bond(atom_C2, cutting_label_L2)
        ]
        
        self.fragment2 = afm.fragment.Fragment()
        for vertex in vertices: self.fragment2.addVertex(vertex)
        for bond in bonds: self.fragment2.addEdge(bond)

    def test_fragment_isomorphism(self):

        self.assertTrue(self.fragment1.isIsomorphic(self.fragment2))

    def test_from_SMILES_like_string1(self):

        # generate fragment from SMILES like string
        # the atom type is also calculated
        smiles_like = 'C'
        fragment = afm.fragment.Fragment().from_SMILES_like_string(smiles_like)
        
        # construct fragment manually
        atom_C = Atom(element=getElement('C'), 
                    radicalElectrons=0, 
                    charge=0, 
                    lonePairs=0)

        atom_H1 = Atom(element=getElement('H'), 
                    radicalElectrons=0, 
                    charge=0, 
                    lonePairs=0)

        atom_H2 = Atom(element=getElement('H'), 
                    radicalElectrons=0, 
                    charge=0, 
                    lonePairs=0)

        atom_H3 = Atom(element=getElement('H'), 
                    radicalElectrons=0, 
                    charge=0, 
                    lonePairs=0)

        atom_H4 = Atom(element=getElement('H'), 
                    radicalElectrons=0, 
                    charge=0, 
                    lonePairs=0)

        atom_C.atomType=atomTypes['Cs']
        atom_H1.atomType=atomTypes['H']
        atom_H2.atomType=atomTypes['H']
        atom_H3.atomType=atomTypes['H']
        atom_H4.atomType=atomTypes['H']

        vertices = [
            atom_C,
            atom_H1,
            atom_H2,
            atom_H3,
            atom_H4
        ]

        bonds = [
            Bond(atom_C, atom_H1, 1),
            Bond(atom_C, atom_H2, 1),
            Bond(atom_C, atom_H3, 1),
            Bond(atom_C, atom_H4, 1)
        ]
        
        expected_fragment = afm.fragment.Fragment()
        for vertex in vertices: expected_fragment.addVertex(vertex)
        for bond in bonds: expected_fragment.addEdge(bond)
	expected_fragment.update()

        self.assertTrue(expected_fragment.isIsomorphic(fragment))

    def test_from_SMILES_like_string2(self):

        # generate fragment from SMILES like string
        # the atom type is also calculated
        smiles_like = 'RCR'
        fragment = afm.fragment.Fragment().from_SMILES_like_string(smiles_like)

        atom_C = Atom(element=getElement('C'), 
                    radicalElectrons=0, 
                    charge=0, 
                    lonePairs=0)

        atom_H1 = Atom(element=getElement('H'), 
                    radicalElectrons=0, 
                    charge=0, 
                    lonePairs=0)

        atom_H2 = Atom(element=getElement('H'), 
                    radicalElectrons=0, 
                    charge=0, 
                    lonePairs=0)

        # construct fragment manually
        atom_C.atomType=atomTypes['Cs']
        atom_H1.atomType=atomTypes['H']
        atom_H2.atomType=atomTypes['H']

        cutting_label_R1 = afm.fragment.CuttingLabel('R')
        cutting_label_R2 = afm.fragment.CuttingLabel('R')
        
        vertices = [
            atom_C,
            cutting_label_R1,
            cutting_label_R2,
            atom_H1,
            atom_H2
        ]

        bonds = [
            Bond(atom_C, cutting_label_R1, 1),
            Bond(atom_C, cutting_label_R2, 1),
            Bond(atom_C, atom_H1, 1),
            Bond(atom_C, atom_H2, 1)
        ]
        
        expected_fragment = afm.fragment.Fragment()
        for vertex in vertices: expected_fragment.addVertex(vertex)
        for bond in bonds: expected_fragment.addEdge(bond)
	expected_fragment.update()

        self.assertTrue(expected_fragment.isIsomorphic(fragment))

    def test_from_SMILES_like_string3(self):

        # generate fragment from SMILES like string
        # the atom type is also calculated
        smiles_like = 'RCL'
        fragment = afm.fragment.Fragment().from_SMILES_like_string(smiles_like)

        atom_C = Atom(element=getElement('C'), 
                    radicalElectrons=0, 
                    charge=0, 
                    lonePairs=0)

        atom_H1 = Atom(element=getElement('H'), 
                    radicalElectrons=0, 
                    charge=0, 
                    lonePairs=0)

        atom_H2 = Atom(element=getElement('H'), 
                    radicalElectrons=0, 
                    charge=0, 
                    lonePairs=0)

        # construct fragment manually
        atom_C.atomType=atomTypes['Cs']
        atom_H1.atomType=atomTypes['H']
        atom_H2.atomType=atomTypes['H']

        cutting_label_R = afm.fragment.CuttingLabel('R')
        cutting_label_L = afm.fragment.CuttingLabel('L')
        
        vertices = [
            atom_C,
            cutting_label_R,
            cutting_label_L,
            atom_H1,
            atom_H2
        ]

        bonds = [
            Bond(atom_C, cutting_label_R, 1),
            Bond(atom_C, cutting_label_L, 1),
            Bond(atom_C, atom_H1, 1),
            Bond(atom_C, atom_H2, 1)
        ]
        
        expected_fragment = afm.fragment.Fragment()
        for vertex in vertices: expected_fragment.addVertex(vertex)
        for bond in bonds: expected_fragment.addEdge(bond)
	expected_fragment.update()

        self.assertTrue(expected_fragment.isIsomorphic(fragment))

    def test_isSubgraphIsomorphic1(self):

        from rmgpy.molecule.group import Group

        smiles_like = '[CH2]CR'
        fragment = afm.fragment.Fragment().from_SMILES_like_string(smiles_like)

        adj = """
1 * R u1
"""
        other = Group().fromAdjacencyList(adj)

        self.assertTrue(fragment.isSubgraphIsomorphic(other))

    def test_isSubgraphIsomorphic2(self):

        from rmgpy.molecule.group import Group

        smiles_like = '[CH2]CR'
        fragment = afm.fragment.Fragment().from_SMILES_like_string(smiles_like)

        adj = """
1 * Ct  u1 {2,T}
2   R!H u0 {1,T}
"""
        other = Group().fromAdjacencyList(adj)

        self.assertFalse(fragment.isSubgraphIsomorphic(other))

    def test_isSubgraphIsomorphic3(self):

        from rmgpy.molecule.group import Group

        smiles_like = '[CH2]CR'
        fragment = afm.fragment.Fragment().from_SMILES_like_string(smiles_like)

        adj = """
1 * R u1
"""
        other = Group().fromAdjacencyList(adj)

        # create initial map
        frag_atom_star = None
        for frag_atom in fragment.vertices:
            if isinstance(frag_atom, Atom) and frag_atom.radicalElectrons == 1:
                frag_atom_star = frag_atom
                break

        group_atom_star = other.vertices[0]
        initialMap = {frag_atom_star: group_atom_star}

        self.assertTrue(fragment.isSubgraphIsomorphic(other, initialMap=initialMap))

    def test_isSubgraphIsomorphic4(self):

        from rmgpy.molecule.group import Group

        smiles_like = '[CH2]CR'
        fragment = afm.fragment.Fragment().from_SMILES_like_string(smiles_like)

        fragment.assign_representative_molecule()

        adj = """
1 * Cs u1 {2,S}
2   N u0 {1,S}
"""
        other = Group().fromAdjacencyList(adj)

        # create initial map
        frag_atom_star = None
        for frag_atom in fragment.vertices:
            if isinstance(frag_atom, Atom) and frag_atom.radicalElectrons == 1:
                frag_atom_star = frag_atom
                break

        group_atom_star = other.vertices[0]
        initialMap = {frag_atom_star: group_atom_star}

        self.assertFalse(fragment.isSubgraphIsomorphic(other, initialMap=initialMap))

    def test_assign_representative_species(self):

        smiles_like = 'RCR'
        fragment = afm.fragment.Fragment().from_SMILES_like_string(smiles_like)

        fragment.assign_representative_species()

        expected_repr_spec = Species().fromSMILES('CCCCC')

        self.assertTrue(expected_repr_spec.isIsomorphic(fragment.species_repr))

    def test_assign_representative_molecule(self):

        smiles_like = 'RCR'
        fragment = afm.fragment.Fragment().from_SMILES_like_string(smiles_like)

        fragment.assign_representative_molecule()

        expected_repr_mol = Molecule().fromSMILES('CCCCC')

        self.assertTrue(expected_repr_mol.isIsomorphic(fragment.mol_repr))

    def test_getMolecularWeight1(self):

        fragmental_weight = self.fragment1.getMolecularWeight()
        self.assertAlmostEqual(fragmental_weight*1000, 12.01, 2)

    def test_getMolecularWeight2(self):

        smiles_like = 'RCR'
        fragment = afm.fragment.Fragment().from_SMILES_like_string(smiles_like)
        fragmental_weight = fragment.getMolecularWeight()
        self.assertAlmostEqual(fragmental_weight*1000, 14.03, 2)

    def test_updateAtomTypes(self):

        atom_C = Atom(element=getElement('C'), 
                    radicalElectrons=0, 
                    charge=0, 
                    lonePairs=0)

        atom_H1 = Atom(element=getElement('H'), 
                    radicalElectrons=0, 
                    charge=0, 
                    lonePairs=0)

        atom_H2 = Atom(element=getElement('H'), 
                    radicalElectrons=0, 
                    charge=0, 
                    lonePairs=0)

        cutting_label_R1 = afm.fragment.CuttingLabel('R')
        cutting_label_R2 = afm.fragment.CuttingLabel('R')
        
        vertices = [
            atom_C,
            cutting_label_R1,
            cutting_label_R2,
            atom_H1,
            atom_H2
        ]

        bonds = [
            Bond(atom_C, cutting_label_R1, 1),
            Bond(atom_C, cutting_label_R2, 1),
            Bond(atom_C, atom_H1, 1),
            Bond(atom_C, atom_H2, 1)
        ]
        
        fragment = afm.fragment.Fragment()
        for vertex in vertices: fragment.addVertex(vertex)
        for bond in bonds: fragment.addEdge(bond)

        fragment.updateAtomTypes()

        for v in fragment.vertices:
            if isinstance(v, Atom) and v.isCarbon():
                break

        self.assertTrue(v.atomType == atomTypes['Cs'])

    def test_update(self):

        atom_C = Atom(element=getElement('C'), 
                    radicalElectrons=0, 
                    charge=0, 
                    lonePairs=0)

        atom_H1 = Atom(element=getElement('H'), 
                    radicalElectrons=0, 
                    charge=0, 
                    lonePairs=0)

        atom_H2 = Atom(element=getElement('H'), 
                    radicalElectrons=0, 
                    charge=0, 
                    lonePairs=0)

        cutting_label_R1 = afm.fragment.CuttingLabel('R')
        cutting_label_R2 = afm.fragment.CuttingLabel('R')
        
        vertices = [
            atom_C,
            cutting_label_R1,
            cutting_label_R2,
            atom_H1,
            atom_H2
        ]

        bonds = [
            Bond(atom_C, cutting_label_R1, 1),
            Bond(atom_C, cutting_label_R2, 1),
            Bond(atom_C, atom_H1, 1),
            Bond(atom_C, atom_H2, 1)
        ]
        
        fragment = afm.fragment.Fragment()
        for vertex in vertices: fragment.addVertex(vertex)
        for bond in bonds: fragment.addEdge(bond)

        fragment.update()

        for v in fragment.vertices:
            if isinstance(v, Atom) and v.isCarbon():
                break

        self.assertTrue(v.atomType == atomTypes['Cs'])
        self.assertTrue(fragment.getNetCharge() == 0)
        self.assertTrue(fragment.multiplicity == 1)

    def test_toAdjacencyList1(self):

        # generate fragment from SMILES like string
        # the atom type is also calculated
        smiles_like = 'C'
        fragment = afm.fragment.Fragment().from_SMILES_like_string(smiles_like)

        fragment.update()
        adj = fragment.toAdjacencyList(removeH=True)
        expected_adj = """1 C u0 p0 c0
"""
        self.assertEqual(adj, expected_adj)

    def test_toAdjacencyList2(self):

        # generate fragment from SMILES like string
        # the atom type is also calculated
        smiles_like = 'CR'
        fragment = afm.fragment.Fragment().from_SMILES_like_string(smiles_like)

        fragment.update()
        adj = fragment.toAdjacencyList(removeH=True)
        expected_adj = """1 C u0 p0 c0 {2,S}
2 R u0 p0 c0 {1,S}
"""
        self.assertEqual(adj, expected_adj)

    def test_toAdjacencyList3(self):

        # generate fragment from SMILES like string
        # the atom type is also calculated
        smiles_like = 'CR'
        fragment = afm.fragment.Fragment().from_SMILES_like_string(smiles_like)

        fragment.update()
        adj = fragment.toAdjacencyList()
        expected_adj = """1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 R u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
"""
        self.assertEqual(adj, expected_adj)

    def test_toAdjacencyList4(self):

        # generate fragment from SMILES like string
        # the atom type is also calculated
        smiles_like = '[CH2]R'
        fragment = afm.fragment.Fragment().from_SMILES_like_string(smiles_like)

        fragment.update()
        adj = fragment.toAdjacencyList()
        expected_adj = """multiplicity 2
1 C u1 p0 c0 {2,S} {3,S} {4,S}
2 R u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
"""
        self.assertEqual(adj, expected_adj)

    def test_fromAdjacencyList1(self):

        adj = """multiplicity 2
1 C u1 p0 c0 {2,S} {3,S} {4,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 R u0 p0 c0 {1,S}
"""
        fragment = afm.fragment.Fragment().fromAdjacencyList(adj)

        # create expected fragment
        smiles_like = '[CH2]R'
        expected_fragment = afm.fragment.Fragment().from_SMILES_like_string(smiles_like)
        expected_fragment.update()

        self.assertTrue(isinstance(fragment.multiplicity, int))
        self.assertTrue(fragment.multiplicity == 2)
        self.assertTrue(fragment.isIsomorphic(expected_fragment))

    def test_fromAdjacencyList2(self):

        adj = """1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 R u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
"""
        fragment = afm.fragment.Fragment().fromAdjacencyList(adj)

        # create expected fragment
        smiles_like = 'CR'
        expected_fragment = afm.fragment.Fragment().from_SMILES_like_string(smiles_like)
        expected_fragment.update()

        self.assertTrue(isinstance(fragment.multiplicity, int))
        self.assertTrue(fragment.multiplicity == 1)
        self.assertTrue(fragment.isIsomorphic(expected_fragment))

    def test_fromAdjacencyList3(self):

        adj = """1 C u0 p0 c0 {2,S}
2 R u0 p0 c0 {1,S}
"""
        fragment = afm.fragment.Fragment().fromAdjacencyList(adj, saturateH=True)

        # create expected fragment
        smiles_like = 'CR'
        expected_fragment = afm.fragment.Fragment().from_SMILES_like_string(smiles_like)
        expected_fragment.update()

        self.assertTrue(isinstance(fragment.multiplicity, int))
        self.assertTrue(fragment.multiplicity == 1)
        self.assertTrue(fragment.isIsomorphic(expected_fragment))

    def test_getAromaticRings(self):

        adj = """1  C u0 p0 c0 {2,D} {6,S} {8,S}
2  C u0 p0 c0 {1,D} {3,S} {9,S}
3  C u0 p0 c0 {2,S} {4,D} {10,S}
4  C u0 p0 c0 {3,D} {5,S} {11,S}
5  C u0 p0 c0 {4,S} {6,D} {12,S}
6  C u0 p0 c0 {1,S} {5,D} {7,S}
7  R u0 p0 c0 {6,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
"""
        fragment = afm.fragment.Fragment().fromAdjacencyList(adj)

        # create expected fragment
        aromaticRings, aromaticBonds = fragment.getAromaticRings()
        self.assertEqual(len(aromaticRings), 1)
        self.assertEqual(len(aromaticRings[0]), 6)
        self.assertEqual(len(aromaticBonds), 1)
        self.assertEqual(len(aromaticBonds[0]), 6)

        aromaticRing_atomSet = set(aromaticRings[0])
        aromaticBonds_set = set(aromaticBonds[0])

        expected_aromaticRing_atomSet = set()
        for v in fragment.vertices:
            if v.isCarbon():
                expected_aromaticRing_atomSet.add(v)

        expected_aromaticBonds_set = set()
        expected_aromaticRing_atomList = list(expected_aromaticRing_atomSet)
        for i, atom1 in enumerate(expected_aromaticRing_atomList):
            for atom2 in expected_aromaticRing_atomList[i+1:]:
                try:
                    bond = fragment.getBond(atom1, atom2)
                    expected_aromaticBonds_set.add(bond)
                except ValueError:
                    pass

        self.assertEqual(aromaticRing_atomSet, expected_aromaticRing_atomSet)
        self.assertEqual(aromaticBonds_set, expected_aromaticBonds_set)

    def test_generate_resonance_structures1(self):

        adj = """1  C u0 p0 c0 {2,S} {3,S} {11,S} {12,S}
2  C u0 p0 c0 {1,S} {4,S} {13,S} {14,S}
3  C u0 p0 c0 {1,S} {5,S} {15,S} {16,S}
4  C u0 p0 c0 {2,S} {17,S} {18,S} {19,S}
5  C u0 p0 c0 {3,S} {6,S} {7,D}
6  C u0 p0 c0 {5,S} {8,D} {20,S}
7  C u1 p0 c0 {5,D} {10,S}
8  C u0 p0 c0 {6,D} {9,S} {21,S}
9  C u0 p0 c0 {8,S} {10,D} {22,S}
10 C u0 p0 c0 {7,S} {9,D} {23,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 R u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {10,S}
"""
        fragment = afm.fragment.Fragment().fromAdjacencyList(adj)

        frag_res = resonance.generate_resonance_structures(fragment, 
                                                           clar_structures=False)

        self.assertEqual(len(frag_res), 2)

    def test_generate_resonance_structures2(self):

        adj = """1  C u0 p0 c0 {2,D} {10,S} {11,S}
2  C u0 p0 c0 {1,D} {3,S} {12,S}
3  C u0 p0 c0 {2,S} {4,D} {13,S}
4  C u0 p0 c0 {3,D} {5,S} {9,S}
5  C u0 p0 c0 {4,S} {6,D} {14,S}
6  C u0 p0 c0 {5,D} {7,S} {15,S}
7  C u0 p0 c0 {6,S} {8,D} {16,S}
8  C u0 p0 c0 {7,D} {9,S} {17,S}
9  C u0 p0 c0 {4,S} {8,S} {10,D}
10 C u0 p0 c0 {1,S} {9,D} {18,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {10,S}
"""
        fragment = afm.fragment.Fragment().fromAdjacencyList(adj)

        frag_res = resonance.generate_resonance_structures(fragment, 
                                                           clar_structures=True)

        self.assertEqual(len(frag_res), 3)
