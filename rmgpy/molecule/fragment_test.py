import os
import unittest

from rmgpy.species import Species
from rmgpy.molecule import resonance
from rmgpy.molecule.atomtype import ATOMTYPES
from rmgpy.molecule.element import get_element
from rmgpy.molecule.molecule import Atom, Bond, Molecule

import rmgpy.molecule.fragment

class TestCuttingLabel(unittest.TestCase):

    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.cutting_label_R = rmgpy.molecule.fragment.CuttingLabel('R')

    def test_symbol(self):

        self.assertEqual('R', self.cutting_label_R.symbol)

    def test_copy(self):

        cutting_label_R_copy = self.cutting_label_R.copy()

        self.assertEqual('R', cutting_label_R_copy.name)
        self.assertEqual(self.cutting_label_R.label, 
                         cutting_label_R_copy.label)
        self.assertEqual(self.cutting_label_R.charge, 
                         cutting_label_R_copy.charge)
        self.assertEqual(self.cutting_label_R.radical_electrons, 
                         cutting_label_R_copy.radical_electrons)
        self.assertEqual(self.cutting_label_R.lone_pairs, 
                         cutting_label_R_copy.lone_pairs)
        self.assertEqual(self.cutting_label_R.isotope, 
                         cutting_label_R_copy.isotope)

class TestFragment(unittest.TestCase):
    def setUp(self):
        """
        A function run before each unit test in this class.
        """

        # construct the first fragment
        atom_C1 = Atom(element=get_element('C'), 
                    radical_electrons=0, 
                    charge=0, 
                    lone_pairs=0)

        cutting_label_R1 = rmgpy.molecule.fragment.CuttingLabel('R')
        cutting_label_L1 = rmgpy.molecule.fragment.CuttingLabel('L')

        vertices = [
            atom_C1,
            cutting_label_R1,
            cutting_label_L1
        ]

        bonds = [
            Bond(atom_C1, cutting_label_R1),
            Bond(atom_C1, cutting_label_L1)
        ]

        self.fragment1 = rmgpy.molecule.fragment.Fragment()
        for vertex in vertices: self.fragment1.add_vertex(vertex)
        for bond in bonds: self.fragment1.add_edge(bond)

        # construct the second fragment
        atom_C2 = Atom(element=get_element('C'), 
                    radical_electrons=0, 
                    charge=0, 
                    lone_pairs=0)

        cutting_label_R2 = rmgpy.molecule.fragment.CuttingLabel('R')
        cutting_label_L2 = rmgpy.molecule.fragment.CuttingLabel('L')

        vertices = [
            atom_C2,
            cutting_label_R2,
            cutting_label_L2
        ]

        bonds = [
            Bond(atom_C2, cutting_label_R2),
            Bond(atom_C2, cutting_label_L2)
        ]

        self.fragment2 = rmgpy.molecule.fragment.Fragment()
        for vertex in vertices: self.fragment2.add_vertex(vertex)
        for bond in bonds: self.fragment2.add_edge(bond)

    def test_fragment_isomorphism(self):

        self.assertTrue(self.fragment1.is_isomorphic(self.fragment2))

    def test_from_smiles_like_string1(self):

        # generate fragment from SMILES like string
        # the atom type is also calculated
        smiles_like = 'C'
        fragment = rmgpy.molecule.fragment.Fragment().from_smiles_like_string(smiles_like)
        
        # construct fragment manually
        atom_C = Atom(element=get_element('C'),
                    radical_electrons=0, 
                    charge=0, 
                    lone_pairs=0)

        atom_H1 = Atom(element=get_element('H'), 
                    radical_electrons=0, 
                    charge=0, 
                    lone_pairs=0)

        atom_H2 = Atom(element=get_element('H'), 
                    radical_electrons=0, 
                    charge=0, 
                    lone_pairs=0)

        atom_H3 = Atom(element=get_element('H'), 
                    radical_electrons=0, 
                    charge=0, 
                    lone_pairs=0)

        atom_H4 = Atom(element=get_element('H'), 
                    radical_electrons=0, 
                    charge=0, 
                    lone_pairs=0)

        atom_C.atomtype=ATOMTYPES['Cs']
        atom_H1.atomtype=ATOMTYPES['H']
        atom_H2.atomtype=ATOMTYPES['H']
        atom_H3.atomtype=ATOMTYPES['H']
        atom_H4.atomtype=ATOMTYPES['H']

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
        
        expected_fragment = rmgpy.molecule.fragment.Fragment()
        for vertex in vertices: expected_fragment.add_vertex(vertex)
        for bond in bonds: expected_fragment.add_edge(bond)
        expected_fragment.update()

        self.assertTrue(expected_fragment.is_isomorphic(fragment))

    def test_from_SMILES_like_string2(self):

        # generate fragment from SMILES like string
        # the atom type is also calculated
        smiles_like = 'RCR'
        fragment = rmgpy.molecule.fragment.Fragment().from_smiles_like_string(smiles_like)

        atom_C = Atom(element=get_element('C'),
                    radical_electrons=0, 
                    charge=0, 
                    lone_pairs=0)

        atom_H1 = Atom(element=get_element('H'), 
                    radical_electrons=0, 
                    charge=0, 
                    lone_pairs=0)

        atom_H2 = Atom(element=get_element('H'), 
                    radical_electrons=0, 
                    charge=0, 
                    lone_pairs=0)

        # construct fragment manually
        atom_C.atomtype=ATOMTYPES['Cs']
        atom_H1.atomtype=ATOMTYPES['H']
        atom_H2.atomtype=ATOMTYPES['H']

        cutting_label_R1 = rmgpy.molecule.fragment.CuttingLabel('R')
        cutting_label_R2 = rmgpy.molecule.fragment.CuttingLabel('R')
        
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
        
        expected_fragment = rmgpy.molecule.fragment.Fragment()
        for vertex in vertices: expected_fragment.add_vertex(vertex)
        for bond in bonds: expected_fragment.add_edge(bond)
        expected_fragment.update()

        self.assertTrue(expected_fragment.is_isomorphic(fragment))

    def test_from_SMILES_like_string3(self):

        # generate fragment from SMILES like string
        # the atom type is also calculated
        smiles_like = 'RCL'
        fragment = rmgpy.molecule.fragment.Fragment().from_smiles_like_string(smiles_like)

        atom_C = Atom(element=get_element('C'),
                    radical_electrons=0, 
                    charge=0, 
                    lone_pairs=0)

        atom_H1 = Atom(element=get_element('H'), 
                    radical_electrons=0, 
                    charge=0, 
                    lone_pairs=0)

        atom_H2 = Atom(element=get_element('H'), 
                    radical_electrons=0, 
                    charge=0, 
                    lone_pairs=0)

        # construct fragment manually
        atom_C.atomtype=ATOMTYPES['Cs']
        atom_H1.atomtype=ATOMTYPES['H']
        atom_H2.atomtype=ATOMTYPES['H']

        cutting_label_R = rmgpy.molecule.fragment.CuttingLabel('R')
        cutting_label_L = rmgpy.molecule.fragment.CuttingLabel('L')
        
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
        
        expected_fragment = rmgpy.molecule.fragment.Fragment()
        for vertex in vertices: expected_fragment.add_vertex(vertex)
        for bond in bonds: expected_fragment.add_edge(bond)
        expected_fragment.update()

        self.assertTrue(expected_fragment.is_isomorphic(fragment))

    def test_is_subgraph_isomorphic1(self):

        from rmgpy.molecule.group import Group

        smiles_like = '[CH2]CR'
        fragment = rmgpy.molecule.fragment.Fragment().from_smiles_like_string(smiles_like)

        adj =     """
                1 * R u1
                  """
        other = Group().from_adjacency_list(adj)

        self.assertTrue(fragment.is_subgraph_isomorphic(other))

    def test_is_subgraph_isomorphic2(self):

        from rmgpy.molecule.group import Group

        smiles_like = '[CH2]CR'
        fragment = rmgpy.molecule.fragment.Fragment().from_smiles_like_string(smiles_like)

        adj =     """
                1 * Ct  u1 {2,T}
                2   R!H u0 {1,T}
                  """
        other = Group().from_adjacency_list(adj)

        self.assertFalse(fragment.is_subgraph_isomorphic(other))

    def test_is_subgraph_isomorphic3(self):

        from rmgpy.molecule.group import Group

        smiles_like = '[CH2]CR'
        fragment = rmgpy.molecule.fragment.Fragment().from_smiles_like_string(smiles_like)

        adj =     """
                1 * R u1
                  """
        other = Group().from_adjacency_list(adj)

        # create initial map
        frag_atom_star = None
        for frag_atom in fragment.vertices:
            if isinstance(frag_atom, Atom) and frag_atom.radical_electrons == 1:
                frag_atom_star = frag_atom
                break

        group_atom_star = other.vertices[0]
        initial_map = {frag_atom_star: group_atom_star}

        self.assertTrue(fragment.is_subgraph_isomorphic(other, initial_map=initial_map))

    def test_is_subgraph_isomorphic4(self):

        from rmgpy.molecule.group import Group

        smiles_like = '[CH2]CR'
        fragment = rmgpy.molecule.fragment.Fragment().from_smiles_like_string(smiles_like)

        fragment.assign_representative_molecule()

        adj =     """
                1 * Cs u1 {2,S}
                2   N u0 {1,S}
                  """
        other = Group().from_adjacency_list(adj)

        # create initial map
        frag_atom_star = None
        for frag_atom in fragment.vertices:
            if isinstance(frag_atom, Atom) and frag_atom.radical_electrons == 1:
                frag_atom_star = frag_atom
                break

        group_atom_star = other.vertices[0]
        initial_map = {frag_atom_star: group_atom_star}

        self.assertFalse(fragment.is_subgraph_isomorphic(other, initial_map=initial_map))

    def test_assign_representative_species(self):

        smiles_like = 'RCR'
        fragment = rmgpy.molecule.fragment.Fragment().from_smiles_like_string(smiles_like)

        fragment.assign_representative_species()

        expected_repr_spec = Species().from_smiles('C=CC(C)(CCCCCCCC(C)(C=C)C(C)C(C)C=CC)C(C)C(C)C=CC')

        self.assertTrue(expected_repr_spec.is_isomorphic(fragment.species_repr))

    def test_assign_representative_molecule(self):

        smiles_like = 'RCR'
        fragment = rmgpy.molecule.fragment.Fragment().from_smiles_like_string(smiles_like)

        fragment.assign_representative_molecule()

        expected_repr_mol = Molecule().from_smiles('C=CC(C)(CCCCCCCC(C)(C=C)C(C)C(C)C=CC)C(C)C(C)C=CC')

        self.assertTrue(expected_repr_mol.is_isomorphic(fragment.mol_repr))

    def test_get_molecular_weight1(self):

        fragmental_weight = self.fragment1.get_molecular_weight()
        self.assertAlmostEqual(fragmental_weight*1000, 12.01, 2)

    def test_get_molecular_weight2(self):

        smiles_like = 'RCR'
        fragment = rmgpy.molecule.fragment.Fragment().from_smiles_like_string(smiles_like)
        fragmental_weight = fragment.get_molecular_weight()
        self.assertAlmostEqual(fragmental_weight*1000, 14.03, 2)

    def test_update_atomtypes(self):

        atom_C = Atom(element=get_element('C'),
                    radical_electrons=0, 
                    charge=0, 
                    lone_pairs=0)

        atom_H1 = Atom(element=get_element('H'), 
                    radical_electrons=0, 
                    charge=0, 
                    lone_pairs=0)

        atom_H2 = Atom(element=get_element('H'), 
                    radical_electrons=0, 
                    charge=0, 
                    lone_pairs=0)

        cutting_label_R1 = rmgpy.molecule.fragment.CuttingLabel('R')
        cutting_label_R2 = rmgpy.molecule.fragment.CuttingLabel('R')
        
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
        
        fragment = rmgpy.molecule.fragment.Fragment()
        for vertex in vertices: fragment.add_vertex(vertex)
        for bond in bonds: fragment.add_edge(bond)

        fragment.update_atomtypes()

        for v in fragment.vertices:
            if isinstance(v, Atom) and v.is_carbon():
                break

        self.assertTrue(v.atomtype == ATOMTYPES['Cs'])

    def test_update(self):

        atom_C = Atom(element=get_element('C'),
                    radical_electrons=0, 
                    charge=0, 
                    lone_pairs=0)

        atom_H1 = Atom(element=get_element('H'), 
                    radical_electrons=0, 
                    charge=0, 
                    lone_pairs=0)

        atom_H2 = Atom(element=get_element('H'), 
                    radical_electrons=0, 
                    charge=0, 
                    lone_pairs=0)

        cutting_label_R1 = rmgpy.molecule.fragment.CuttingLabel('R')
        cutting_label_R2 = rmgpy.molecule.fragment.CuttingLabel('R')
        
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
        
        fragment = rmgpy.molecule.fragment.Fragment()
        for vertex in vertices: fragment.add_vertex(vertex)
        for bond in bonds: fragment.add_edge(bond)

        fragment.update()

        for v in fragment.vertices:
            if isinstance(v, Atom) and v.is_carbon():
                break

        self.assertTrue(v.atomtype == ATOMTYPES['Cs'])
        self.assertTrue(fragment.get_net_charge() == 0)
        self.assertTrue(fragment.multiplicity == 1)

    def test_to_adjacency_list1(self):

        # generate fragment from smiles like string
        # the atom type is also calculated
        smiles_like = 'C'
        fragment = rmgpy.molecule.fragment.Fragment().from_smiles_like_string(smiles_like)

        fragment.update()
        adj = fragment.to_adjacency_list(remove_h=True)
        expected_adj = """1 C u0 p0 c0\n"""
        self.assertEqual(adj, expected_adj)

    def test_to_adjacency_list2(self):

        # generate fragment from smiles like string
        # removed H
        smiles_like = 'CR'
        fragment = rmgpy.molecule.fragment.Fragment().from_smiles_like_string(smiles_like)

        fragment.update()
        adj = fragment.to_adjacency_list(remove_h=True)
        expected_adj = """1 C u0 p0 c0 {2,S}\n2 R u0 p0 c0 {1,S}\n"""
        self.assertEqual(adj, expected_adj)

    def test_to_adjacency_list3(self):

        # generate fragment from smiles like string
        # with H
        smiles_like = 'CR'
        fragment = rmgpy.molecule.fragment.Fragment().from_smiles_like_string(smiles_like)

        fragment.update()
        adj = fragment.to_adjacency_list()
        expected_adj = """1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 R u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
"""
        self.assertEqual(adj, expected_adj)

    def test_to_adjacency_list4(self):

        # generate fragment from SMILES like string
        # radical species
        smiles_like = '[CH2]R'
        fragment = rmgpy.molecule.fragment.Fragment().from_smiles_like_string(smiles_like)

        fragment.update()
        adj = fragment.to_adjacency_list()
        expected_adj = """multiplicity 2
1 C u1 p0 c0 {2,S} {3,S} {4,S}
2 R u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
"""
        self.assertEqual(adj, expected_adj)

    def test_from_adjacency_list1(self):

        adj = """multiplicity 2
1 C u1 p0 c0 {2,S} {3,S} {4,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 R u0 p0 c0 {1,S}
"""
        fragment = rmgpy.molecule.fragment.Fragment().from_adjacency_list(adj)

        # create expected fragment
        smiles_like = '[CH2]R'
        expected_fragment = rmgpy.molecule.fragment.Fragment().from_smiles_like_string(smiles_like)
        expected_fragment.update()

        self.assertTrue(isinstance(fragment.multiplicity, int))
        self.assertTrue(fragment.multiplicity == 2)
        self.assertTrue(fragment.is_isomorphic(expected_fragment))

    def test_from_adjacency_list2(self):

        adj = """1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 R u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
"""
        fragment = rmgpy.molecule.fragment.Fragment().from_adjacency_list(adj)

        # create expected fragment
        smiles_like = 'CR'
        expected_fragment = rmgpy.molecule.fragment.Fragment().from_smiles_like_string(smiles_like)
        expected_fragment.update()

        self.assertTrue(fragment.multiplicity == 1)
        self.assertTrue(fragment.is_isomorphic(expected_fragment))

    def test_from_adjacency_list3(self):

        adj = """1 C u0 p0 c0 {2,S}
2 R u0 p0 c0 {1,S}
"""
        fragment = rmgpy.molecule.fragment.Fragment().from_adjacency_list(adj, saturate_h=True)

        # create expected fragment
        smiles_like = 'CR'
        expected_fragment = rmgpy.molecule.fragment.Fragment().from_smiles_like_string(smiles_like)
        expected_fragment.update()

        self.assertTrue(fragment.multiplicity == 1)
        self.assertTrue(fragment.is_isomorphic(expected_fragment))

    def test_get_aromatic_rings(self):

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
        fragment = rmgpy.molecule.fragment.Fragment().from_adjacency_list(adj)

        # create expected fragment
        aromatic_rings, aromatic_bonds = fragment.get_aromatic_rings()
        self.assertEqual(len(aromatic_rings), 1)
        self.assertEqual(len(aromatic_rings[0]), 6)
        self.assertEqual(len(aromatic_bonds), 1)
        self.assertEqual(len(aromatic_bonds[0]), 6)

        aromatic_ring_atomset = set(aromatic_rings[0])
        aromatic_bonds_set = set(aromatic_bonds[0])

        expected_aromatic_ring_atomset = set()
        for v in fragment.vertices:
            if v.is_carbon():
                expected_aromatic_ring_atomset.add(v)

        expected_aromatic_bonds_set = set()
        expected_aromatic_ring_atomlist = list(expected_aromatic_ring_atomset)
        for i, atom1 in enumerate(expected_aromatic_ring_atomlist):
            for atom2 in expected_aromatic_ring_atomlist[i+1:]:
                try:
                    bond = fragment.get_bond(atom1, atom2)
                    expected_aromatic_bonds_set.add(bond)
                except ValueError:
                    pass

        self.assertEqual(aromatic_ring_atomset, expected_aromatic_ring_atomset)
        self.assertEqual(aromatic_bonds_set, expected_aromatic_bonds_set)

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
        fragment = rmgpy.molecule.fragment.Fragment().from_adjacency_list(adj)

        frag_res = resonance.generate_resonance_structures(fragment, 
                                                           clar_structures=False)

        self.assertEqual(len(frag_res), 2)

    def test_generate_resonance_structures2(self):

        ### use ethylbenzly radical to test if it can generate aromatic fragment or not

        fragment = rmgpy.molecule.fragment.Fragment().from_smiles_like_string('c1ccccc1[CH]CR')

        frag_res = resonance.generate_resonance_structures(fragment, 
                                                           clar_structures=True)

        self.assertTrue(frag_res[0].is_aromatic())    
        self.assertEqual(len(frag_res), 4)

    def test_fragment_is_identical(self):

        self.assertTrue(self.fragment1.is_identical(self.fragment2))


    def test_fragment_get_formula(self):

        fragment = rmgpy.molecule.fragment.Fragment().from_smiles_like_string('CCR')

        self.assertTrue(fragment.get_formula()=='C2H5R')

    def test_fragment_is_linear(self):

        fragment = rmgpy.molecule.fragment.Fragment().from_smiles_like_string('C#CR')

        self.assertTrue(fragment.is_linear())

    def test_fragment_get_element_count(self):

        fragment = rmgpy.molecule.fragment.Fragment().from_smiles_like_string('[CH2]CR')

        self.assertEqual(fragment.get_element_count()['C'], 2)
        self.assertEqual(fragment.get_element_count()['H'], 4)

    def test_fragment_get_num_atoms(self):

        fragment = rmgpy.molecule.fragment.Fragment().from_smiles_like_string('[CH2]CR')

        self.assertEqual(fragment.get_num_atoms(element='C'), 2)
        self.assertEqual(fragment.get_num_atoms(element='H'), 4)

    def test_fragment_to_smiles(self):

        adj = """1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 R u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
"""

        fragment = rmgpy.molecule.fragment.Fragment().from_adjacency_list(adj)
        fragment.update()
        smiles = fragment.to_smiles()
        expected_smiles = 'CR'

        self.assertTrue(smiles, expected_smiles )

    def test_sliceitup_arom1(self):

        # test avoid cutting aromatic species at ring position
        fragment = rmgpy.molecule.fragment.Fragment().from_smiles_like_string('c1ccccc1CCCCCCCC')

        frags = rmgpy.molecule.fragment.Fragment().sliceitup_arom(fragment)

        # check element balance
        expected_element = fragment.get_element_count()
        # after cutting
        total_element = {}
        for frag in frags:
            elements = frag.get_element_count()
            for ele, num in elements.items():
                if ele in total_element:
                    total_element[ele] += num
                else:
                    total_element[ele] = num

        self.assertTrue(frags[0].is_aromatic() or frags[1].is_aromatic())
        self.assertEqual(len(frags), 2)
        self.assertEqual(expected_element, total_element)

    def test_sliceitup_arom2(self):

        # do not cut when input is aliphatic species

        fragment = rmgpy.molecule.fragment.Fragment().from_smiles_like_string('CCCCCCCCCC')
        frags = rmgpy.molecule.fragment.Fragment().sliceitup_arom(fragment)

        self.assertEqual(len(frags), 1)
        self.assertTrue(fragment.is_isomorphic(frags[0]))

    def test_sliceitup_aliph(self):

        # test avoid cutting species at ring position

        import re
        fragment = rmgpy.molecule.fragment.Fragment().from_smiles_like_string('C1=CCC=CC1CCCCCCCC')
        frags = rmgpy.molecule.fragment.Fragment().sliceitup_aliph(fragment)

        # check element balance
        expected_element = fragment.get_element_count()
        # after cutting
        total_element = {}
        for frag in frags:
            elements = frag.get_element_count()
            for ele, num in elements.items():
                if ele in total_element:
                    total_element[ele] += num
                else:
                    total_element[ele] = num

        # check whether ring still exists
        f0 = re.findall(r'\d', frags[0].to_smiles())
        f1 = re.findall(r'\d', frags[1].to_smiles())

        self.assertTrue(f0 != [] or f1 != [])
        self.assertEqual(len(frags), 2)
        self.assertTrue(expected_element == total_element)

    def test_cut_molecule1(self):

        # test output string
        fragment = rmgpy.molecule.fragment.Fragment().from_smiles_like_string('CCCCCCCCCC')
        new_frags = fragment.cut_molecule(output_smiles=True)

        self.assertEqual(len(new_frags), 2)
        self.assertTrue(isinstance(new_frags[0], str))

    def test_cut_molecule2(self):

        # test input Fragment
        fragment = rmgpy.molecule.fragment.Fragment().from_smiles_like_string('CCCCCCCCCCR')
        new_frags = fragment.cut_molecule()

        self.assertEqual(len(new_frags), 2)

    def test_calculate_symmetry_number1(self):

        #  for fragment with 1 CuttingLabel
        fragment = rmgpy.molecule.fragment.Fragment().from_smiles_like_string('CCR')
        fragment.calculate_symmetry_number()
        fsn = fragment.symmetry_number

        self.assertEqual(fsn, 3.0)

    def test_calculate_symmetry_number2(self):

        #  for fragment with 2 CuttingLabel
        fragment = rmgpy.molecule.fragment.Fragment().from_smiles_like_string('LCCR')
        fragment.calculate_symmetry_number()
        fsn = fragment.symmetry_number

        self.assertEqual(fsn, 1.0)

    def test_get_symmetry_number1(self):

        # fragment symmetry number == -1
        fragment = rmgpy.molecule.fragment.Fragment().from_smiles_like_string('CCR')
        fsn = fragment.symmetry_number

        self.assertEqual(fsn, -1)

    def test_get_symmetry_number2(self):

        # fragment symmetry number != -1
        fragment = rmgpy.molecule.fragment.Fragment().from_smiles_like_string('CCR')
        fragment.get_symmetry_number()
        fsn = fragment.symmetry_number

        # fragment.symmetry_number == -1
        self.assertTrue(fsn != -1)
        self.assertEqual(fsn, 3.0)

    def test_is_radical(self):

        fragment = rmgpy.molecule.fragment.Fragment().from_smiles_like_string('[CH2]CR')

        self.assertTrue(fragment.is_radical())

    def test_is_aromatic(self):

        fragment = rmgpy.molecule.fragment.Fragment().from_smiles_like_string('c1ccccc1CR')
        frag = fragment.generate_resonance_structures()[0]

        self.assertTrue(frag.is_aromatic())

    def test_get_representative_molecule(self):

        fragment = rmgpy.molecule.fragment.Fragment().from_smiles_like_string('CCR')
        mol_repr,_ = fragment.get_representative_molecule()
        ethane = Molecule().from_smiles('CC')

        self.assertTrue(mol_repr.is_isomorphic(ethane))

    def test_assign_representative_species(self):

        fragment = rmgpy.molecule.fragment.Fragment().from_smiles_like_string('CCR')
        fragment.assign_representative_species()

        self.assertEqual(fragment.species_repr.symmetry_number, 3.0)
        self.assertEqual(fragment.species_repr.smiles.count('C'), 14+2)

    def test_to_rdkit_mol(self):

        fragment = rmgpy.molecule.fragment.Fragment().from_smiles_like_string('CCR')
        rdmol,_ = fragment.to_rdkit_mol()

        self.assertEqual(rdmol.GetNumAtoms(), 8)
