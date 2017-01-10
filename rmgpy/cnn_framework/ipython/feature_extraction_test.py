import numpy as np
from feature_extraction import *
import unittest
from rmgpy.molecule.molecule import Molecule
from rmgpy.molecule.resonance import generateAromaticResonanceIsomers


class Test_MolecularFeatureExtractor(unittest.TestCase):

	def setUp(self):

		self.feature_extractor = MolecularFeatureExtractor()

	def test_get_i_member_ring_counts(self):

		mol0 = Molecule().fromSMILES('C1[C@H]2C[C@@H]1C2')
		i_min = 3
		i_max = 8
		i_member_ring_counts0 = self.feature_extractor.get_i_member_ring_counts(mol0, i_min, i_max)

		self.assertEqual(i_member_ring_counts0, [0,2,0,0,0,0])

		mol1 = Molecule().fromSMILES('C1C2C3CC4CC2C13C4')
		i_min = 3
		i_max = 8
		i_member_ring_counts1 = self.feature_extractor.get_i_member_ring_counts(mol1, i_min, i_max)

		self.assertEqual(i_member_ring_counts1, [0,2,2,0,0,0])

	def test_get_ij_bicyclic_counts(self):

		mol0 = Molecule().fromSMILES('C1[C@H]2C[C@@H]1C2')
		i_min = 3
		i_max = 8

		ij_bicyclic_counts0 = self.feature_extractor.get_ij_bicyclic_counts(mol0, i_min, i_max)
		expected_counts0 = [0, 0, 0, 0, 0, 0,
							1, 0, 0, 0, 0,
							0, 0, 0, 0,
							0, 0, 0,
							0, 0,
							0]

		self.assertEqual(list(ij_bicyclic_counts0), expected_counts0)

		mol1 = Molecule().fromSMILES('CCCC1C2C3CC3C12')
		i_min = 3
		i_max = 8

		ij_bicyclic_counts1 = self.feature_extractor.get_ij_bicyclic_counts(mol1, i_min, i_max)
		expected_counts1 = [0, 2, 0, 0, 0, 0,
							0, 0, 0, 0, 0,
							0, 0, 0, 0,
							0, 0, 0,
							0, 0,
							0]

		self.assertEqual(list(ij_bicyclic_counts1), expected_counts1)

	def test_get_atoms_shared_by_tri_tetra_cyclics(self):

		mol0 = Molecule().fromSMILES('C1[C@H]2C[C@@H]1C2')
		[count_in_3rings0, count_in_4_more_rings0] = self.feature_extractor.get_atoms_shared_by_tri_tetra_cyclics(mol0)
		expected_counts0 = [0,0]

		self.assertEqual([count_in_3rings0, count_in_4_more_rings0], expected_counts0)

		mol1 = Molecule().fromSMILES('C1C2CC3CC4CC1C234')
		[count_in_3rings1, count_in_4_more_rings1] = self.feature_extractor.get_atoms_shared_by_tri_tetra_cyclics(mol1)
		expected_counts1 = [0,1]

		self.assertEqual([count_in_3rings1, count_in_4_more_rings1], expected_counts1)

	def test_get_bond_order_sharedby_bicyclic_counts(self):

		mol0 = Molecule().fromSMILES('C1[C@H]2C[C@@H]1C2')
		orders = ['S', 'D', 'T', 'B']
		bond_order_sharedby_bicyclic_counts0 = self.feature_extractor.get_bond_order_sharedby_bicyclic_counts(mol0, orders)

		self.assertEqual(list(bond_order_sharedby_bicyclic_counts0), [2,0,0,0])

		mol1 = Molecule().fromSMILES('C1C2CC3CC4CC1C234')
		orders = ['S', 'D', 'T', 'B']
		bond_order_sharedby_bicyclic_counts1 = self.feature_extractor.get_bond_order_sharedby_bicyclic_counts(mol1, orders)

		self.assertEqual(list(bond_order_sharedby_bicyclic_counts1), [4,0,0,0])

	def test_get_bond_order_in_i_member_ring_counts(self):

		mol0 = Molecule().fromSMILES('C1[C@H]2C[C@@H]1C2')
		orders = ['S', 'D', 'T', 'B']
		i_min = 3
		i_max = 8
		bond_order_in_i_member_ring_counts0 = self.feature_extractor.get_bond_order_in_i_member_ring_counts(mol0, i_min, i_max, orders)

		expected_counts0 = [0, 0, 0, 0,
							8, 0, 0, 0,
							0, 0, 0, 0,
							0, 0, 0, 0,
							0, 0, 0, 0,
							0, 0, 0, 0]

		self.assertEqual(list(bond_order_in_i_member_ring_counts0), expected_counts0)

		mol1= Molecule().fromSMILES('CCCC1CC11CCC1')
		orders = ['S', 'D', 'T', 'B']
		i_min = 3
		i_max = 8
		bond_order_in_i_member_ring_counts1 = self.feature_extractor.get_bond_order_in_i_member_ring_counts(mol1, i_min, i_max, orders)

		expected_counts1 = [3, 0, 0, 0,
							4, 0, 0, 0,
							0, 0, 0, 0,
							0, 0, 0, 0,
							0, 0, 0, 0,
							0, 0, 0, 0]

		self.assertEqual(list(bond_order_in_i_member_ring_counts1), expected_counts1)

	def test_get_benzene_ring_count(self):

		mol0 = Molecule().fromSMILES('C1[C@H]2C[C@@H]1C2')
		benzene_count0 = self.feature_extractor.get_benzene_ring_count(mol0)

		self.assertEqual(benzene_count0, 0)

		mol1 = Molecule().fromSMILES('c1ccccc1CC')
		mol_aro = generateAromaticResonanceIsomers(mol1)[0]
		benzene_count1 = self.feature_extractor.get_benzene_ring_count(mol_aro)

		self.assertEqual(benzene_count1, 1)




