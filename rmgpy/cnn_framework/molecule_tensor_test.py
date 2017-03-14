
import unittest
from molecule_tensor import *
from rmgpy.molecule.molecule import Molecule
import numpy as np

class Test_Molecule_Tensor(unittest.TestCase):

	def test_get_atom_attributes(self):

		mol_test = Molecule().fromSMILES('C')
		non_H_atoms = [atom0 for atom0 in mol_test.atoms if not atom0.isHydrogen()]
		
		atom_attributes_dict = get_atom_attributes(mol_test, non_H_atoms)

		self.assertEqual(len(atom_attributes_dict), 1)

		expected_attributes = [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
								1, 0, 0, 0, 0, 0,
								0, 0, 0, 0, 1,
								0,
								0,
								0]

		expected_attributes = np.array(expected_attributes, dtype=np.float32)
		for test_entry, expected_entry in zip(atom_attributes_dict[non_H_atoms[0]], expected_attributes):
			self.assertAlmostEqual(test_entry, expected_entry, 4)

	def test_is_atom_in_rings(self):

		mol_test = Molecule().fromAdjacencyList(
"""1  C u0 p0 c0 {2,S} {11,S} {12,S} {13,S}
2  C u0 p0 c0 {1,S} {3,S} {14,S} {15,S}
3  C u0 p0 c0 {2,S} {4,S} {10,S} {16,S}
4  C u0 p0 c0 {3,S} {5,S} {9,S} {17,S}
5  C u0 p0 c0 {4,S} {6,S} {18,S} {19,S}
6  S u0 p2 c0 {5,S} {7,S}
7  C u0 p0 c0 {6,S} {8,S} {10,S} {20,S}
8  C u0 p0 c0 {7,S} {9,S} {21,S} {22,S}
9  C u0 p0 c0 {4,S} {8,S} {23,S} {24,S}
10 C u0 p0 c0 {3,S} {7,S} {11,S} {25,S}
11 C u0 p0 c0 {1,S} {10,S} {26,S} {27,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {10,S}
26 H u0 p0 c0 {11,S}
27 H u0 p0 c0 {11,S}
"""
)
		atom5 = mol_test.atoms[5]
		atom5_in_rings = is_atom_in_ring(mol_test, atom5)
		expected_atom5_in_rings = [0, 0, 0, 2, 0, 0]

		self.assertEqual(atom5_in_rings, expected_atom5_in_rings)

	def test_is_bond_conjugated_1(self):

		mol_test = Molecule().fromAdjacencyList(
"""1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
"""
)
		atom0 = mol_test.atoms[0]
		atom1 = mol_test.atoms[1]
		bond = atom0.bonds[atom1]

		self.assertFalse(is_bond_conjugated(bond))

	def test_is_bond_conjugated_2(self):

		mol_test = Molecule().fromAdjacencyList(
"""1  C u0 p0 c0 {2,D} {5,S} {6,S}
2  C u0 p0 c0 {1,D} {3,S} {7,S}
3  C u0 p0 c0 {2,S} {4,D} {8,S}
4  C u0 p0 c0 {3,D} {9,S} {10,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
"""
)
		atom0 = mol_test.atoms[0]
		atom1 = mol_test.atoms[1]
		atom2 = mol_test.atoms[2]
		atom3 = mol_test.atoms[3]
		atom4 = mol_test.atoms[4]
		bond01 = atom0.bonds[atom1]
		bond12 = atom1.bonds[atom2]
		bond23 = atom2.bonds[atom3]
		bond04 = atom0.bonds[atom4]

		self.assertTrue(is_bond_conjugated(bond01))
		self.assertTrue(is_bond_conjugated(bond12))
		self.assertTrue(is_bond_conjugated(bond23))
		self.assertFalse(is_bond_conjugated(bond04))

	def test_is_bond_conjugated_3(self):

		mol_test = Molecule().fromAdjacencyList(
"""1  C u0 p0 c0 {2,D} {6,S} {7,S}
2  C u0 p0 c0 {1,D} {3,S} {8,S}
3  C u0 p0 c0 {2,S} {4,D} {9,S}
4  C u0 p0 c0 {3,D} {5,S} {10,S}
5  C u0 p0 c0 {4,S} {6,D} {11,S}
6  C u0 p0 c0 {1,S} {5,D} {12,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
"""
)
		atom0 = mol_test.atoms[0]
		atom1 = mol_test.atoms[1]
		atom6 = mol_test.atoms[6]
		bond01 = atom0.bonds[atom1]
		bond06 = atom0.bonds[atom6]

		self.assertTrue(is_bond_conjugated(bond01))
		self.assertFalse(is_bond_conjugated(bond06))

	def test_is_bond_conjugated_4(self):

		mol_test = Molecule().fromAdjacencyList(
"""1  C u0 p0 c0 {2,D} {8,S} {9,S}
2  C u0 p0 c0 {1,D} {3,S} {10,S}
3  C u0 p0 c0 {2,S} {4,S} {11,S} {12,S}
4  C u0 p0 c0 {3,S} {5,D} {13,S}
5  C u0 p0 c0 {4,D} {6,S} {14,S}
6  C u0 p0 c0 {5,S} {7,D} {15,S}
7  C u0 p0 c0 {6,D} {16,S} {17,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
"""
)
		atom0 = mol_test.atoms[0]
		atom1 = mol_test.atoms[1]
		atom2 = mol_test.atoms[2]
		atom3 = mol_test.atoms[3]
		atom4 = mol_test.atoms[4]
		atom5 = mol_test.atoms[5]
		atom6 = mol_test.atoms[6]
		bond01 = atom0.bonds[atom1]
		bond12 = atom1.bonds[atom2]
		bond45 = atom4.bonds[atom5]
		bond56 = atom5.bonds[atom6]

		self.assertFalse(is_bond_conjugated(bond01))
		self.assertFalse(is_bond_conjugated(bond12))
		self.assertTrue(is_bond_conjugated(bond45))
		self.assertTrue(is_bond_conjugated(bond56))

	def test_is_bond_in_ring(self):

		mol_test = Molecule().fromAdjacencyList(
"""1  C u0 p0 c0 {2,S} {11,S} {12,S} {13,S}
2  C u0 p0 c0 {1,S} {3,S} {14,S} {15,S}
3  C u0 p0 c0 {2,S} {4,S} {10,S} {16,S}
4  C u0 p0 c0 {3,S} {5,S} {9,S} {17,S}
5  C u0 p0 c0 {4,S} {6,S} {18,S} {19,S}
6  S u0 p2 c0 {5,S} {7,S}
7  C u0 p0 c0 {6,S} {8,S} {10,S} {20,S}
8  C u0 p0 c0 {7,S} {9,S} {21,S} {22,S}
9  C u0 p0 c0 {4,S} {8,S} {23,S} {24,S}
10 C u0 p0 c0 {3,S} {7,S} {11,S} {25,S}
11 C u0 p0 c0 {1,S} {10,S} {26,S} {27,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {10,S}
26 H u0 p0 c0 {11,S}
27 H u0 p0 c0 {11,S}
"""
)
		atom4 = mol_test.atoms[4]
		atom5 = mol_test.atoms[5]
		bond45 = atom4.bonds[atom5]
		bond45_in_rings = is_bond_in_ring(mol_test, bond45)
		expected_bond45_in_rings = [0, 0, 0, 2, 0, 0]

		self.assertEqual(bond45_in_rings, expected_bond45_in_rings)

	def test_get_bond_attributes_1(self):

		mol_test = Molecule().fromSMILES('CC')
		non_H_atoms = [atom0 for atom0 in mol_test.atoms if not atom0.isHydrogen()]
		
		bond_attributes_dict = get_bond_attributes(mol_test, non_H_atoms)

		self.assertEqual(len(bond_attributes_dict), 1)

		expected_attributes = [1, 0, 0, 0, 
								0,
								0,
								0,
								1]
		expected_attributes = np.array(expected_attributes, dtype=np.float32)
		for test_entry, expected_entry in zip(bond_attributes_dict[non_H_atoms[0].bonds[non_H_atoms[1]]], expected_attributes):
			self.assertAlmostEqual(test_entry, expected_entry, 4)

	def test_get_bond_attributes_2(self):

		mol_test = Molecule().fromAdjacencyList(
"""1  C u0 p0 c0 {2,B} {6,B} {7,S}
2  C u0 p0 c0 {1,B} {3,B} {8,S}
3  C u0 p0 c0 {2,B} {4,B} {9,S}
4  C u0 p0 c0 {3,B} {5,B} {10,S}
5  C u0 p0 c0 {4,B} {6,B} {11,S}
6  C u0 p0 c0 {1,B} {5,B} {12,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
"""
)
		non_H_atoms = [atom0 for atom0 in mol_test.atoms if not atom0.isHydrogen()]
		
		bond_attributes_dict = get_bond_attributes(mol_test, non_H_atoms)

		self.assertEqual(len(bond_attributes_dict), 6)

		expected_attributes = [0, 1, 0, 0, 
								1,
								1,
								1,
								1]
		expected_attributes = np.array(expected_attributes, dtype=np.float32)
		for test_entry, expected_entry in zip(bond_attributes_dict[non_H_atoms[0].bonds[non_H_atoms[1]]], expected_attributes):
			self.assertAlmostEqual(test_entry, expected_entry, 4)
	
	def test_get_molecule_tensor_1(self):

		mol_test = Molecule().fromSMILES('C')

		mol_tensor_test = get_molecule_tensor(mol_test)

		self.assertEqual(len(mol_tensor_test), 1)
		self.assertEqual(len(mol_tensor_test[0]), 1)
		self.assertEqual(len(mol_tensor_test[0][0]), get_attribute_vector_size())

	def test_get_molecule_tensor_2(self):

		mol_test = Molecule().fromSMILES('CC')

		mol_tensor_test = get_molecule_tensor(mol_test)

		self.assertEqual(len(mol_tensor_test), 2)
		self.assertEqual(len(mol_tensor_test[0]), 2)
		self.assertEqual(len(mol_tensor_test[0][0]), get_attribute_vector_size())


