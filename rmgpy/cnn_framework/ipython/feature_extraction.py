import numpy as np
from rmgpy.data.thermo import bicyclicDecompositionForPolyring

class MolecularFeatureExtractor():

	def get_i_member_ring_counts(self, mol, i_min, i_max):

		sssr = mol.getSmallestSetOfSmallestRings()
		sssr_sizes = [len(ring) for ring in sssr]
		i_member_ring_counts = [sssr_sizes.count(i) for i in range(i_min, i_max+1)]

		return i_member_ring_counts


	def get_ij_bicyclic_counts(self, mol, i_min, i_max):

		ij_bicyclic_counts = np.zeros((i_max-i_min+1)*(i_max-i_min+2)/2)
		_, polyrings = mol.getDisparateRings()
		for polyring in polyrings:
			bicyclicsMergedFromRingPair, _ = bicyclicDecompositionForPolyring(polyring)
			for bicyclic in bicyclicsMergedFromRingPair:
				sssr = bicyclic.getSmallestSetOfSmallestRings()
				sssr_sizes = sorted([len(ring) for ring in sssr])
				assert len(sssr_sizes)==2
				i = sssr_sizes[0]
				j = sssr_sizes[1]
				
				idx = (i-i_min)*(i_max-i_min+1 + i_max-i_min+1 - (i-i_min-1))/2 + (j-(i_min+(i-i_min)))
				ij_bicyclic_counts[idx] += 1

		return ij_bicyclic_counts

	def get_atoms_shared_by_tri_tetra_cyclics(self, mol):

		sssr = mol.getSmallestSetOfSmallestRings()
		atoms = set([])
		for ring in sssr:
			atoms = atoms.union(set(ring))
		atoms = list(atoms)
		count_in_3rings = 0
		count_in_4_more_rings = 0
		for atom in atoms:
			count_in_ring = 0
			for ring in sssr:
				if atom in ring:
					count_in_ring += 1
			if count_in_ring == 3:
				count_in_3rings += 1
			if count_in_ring >= 4:
				count_in_4_more_rings += 1

		return [count_in_3rings, count_in_4_more_rings]

	def get_bond_order_sharedby_bicyclic_counts(self, mol, orders):

		bond_order_sharedby_bicyclic_counts = np.zeros(len(orders))

		_, polyrings = mol.getDisparateRings()
		for polyring in polyrings:
			bicyclicsMergedFromRingPair, _ = bicyclicDecompositionForPolyring(polyring)
			for bicyclic in bicyclicsMergedFromRingPair:
				sssr = bicyclic.getSmallestSetOfSmallestRings()
				assert len(sssr)==2
				# get the atoms shared by both rings
				shared_atoms = []
				ring1 = sssr[0]
				ring2 = sssr[1]
				for atom in ring1:
					if atom in ring2 and atom not in shared_atoms:
						shared_atoms.append(atom)
				for atom in ring2:
					if atom in ring1 and atom not in shared_atoms:
						shared_atoms.append(atom)
				
				# check double bond existance between those shared atoms
				shared_bonds = []
				for idx_i in range(len(shared_atoms)):
					atom_i = shared_atoms[idx_i]
					for idx_j in range(idx_i+1, len(shared_atoms)):
						atom_j = shared_atoms[idx_j]
						if atom_j in atom_i.bonds:
							if atom_i.bonds[atom_j] not in shared_bonds:
								shared_bonds.append(atom_i.bonds[atom_j])
				
				bond_orders = [bond.getOrderStr() for bond in shared_bonds]
				bond_orders_counts = np.array([bond_orders.count(order) for order in orders])
				bond_order_sharedby_bicyclic_counts += bond_orders_counts

		return bond_order_sharedby_bicyclic_counts

	def get_bond_order_in_i_member_ring_counts(self, mol, i_min, i_max, orders):

		bond_order_in_i_member_ring_counts = np.zeros(len(orders)*(i_max-i_min+1))
		sssr = mol.getSmallestSetOfSmallestRings()
		for ring in sssr:
			bonds_in_unicyclic = []
			for idx_i in range(len(ring)):
				atom_i = ring[idx_i]
				for idx_j in range(idx_i+1, len(ring)):
					atom_j = ring[idx_j]
					if atom_j in atom_i.bonds:
						if atom_i.bonds[atom_j] not in bonds_in_unicyclic:
							bonds_in_unicyclic.append(atom_i.bonds[atom_j])
		
			bond_orders = [bond.getOrderStr() for bond in bonds_in_unicyclic]
			bond_orders_counts = np.array([bond_orders.count(order) for order in orders])
			ring_size = len(ring)
			bond_order_in_i_member_ring_counts[((ring_size-i_min)*len(orders)):((ring_size-i_min+1)*len(orders))] += bond_orders_counts

		return bond_order_in_i_member_ring_counts

	def get_benzene_ring_count(self, mol):

		sssr = mol.getSmallestSetOfSmallestRings()
		benzene_count = 0
		for ring in sssr:
			is_benzene = True
			for ring_atom in ring:
				for bonded_atom, bond in ring_atom.bonds.iteritems():
					if bonded_atom in ring:
						if not bond.isBenzene():
							is_benzene = False 
							break
			if is_benzene:
				benzene_count += 1
		return benzene_count