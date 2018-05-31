
import numpy as np

def pad_molecule_tensor(molecule_tensor, final_size):
	"""
	this method takes a numpy tensor with dimension:
	non_H_atom_num * non_H_atom_num * attribute_num

	and pad it to final_size * final_size * attribute_num
	using zeros
	"""
	padded_molecule_tensor = np.zeros((final_size, final_size, molecule_tensor.shape[2]))

	padded_molecule_tensor[:molecule_tensor.shape[0], :molecule_tensor.shape[1], :] = molecule_tensor
	return padded_molecule_tensor



def get_molecule_tensor(molecule, 
						add_extra_atom_attribute=True, 
						add_extra_bond_attribute=True,
						differentiate_atom_type=True,
						differentiate_bond_type=True):

	"""
	this method takes RMG `Molecule` object and vectorize 
	it into a numpy molecule tensor with dimension: 
	non_H_atom_num * non_H_atom_num * attribute_num
	"""
	non_H_atoms = [atom0 for atom0 in molecule.atoms if not atom0.isHydrogen()]
	non_H_atom_num = len(non_H_atoms)
	
	atom_attributes_dict = get_atom_attributes(molecule, 
											   non_H_atoms, 
											   add_extra_atom_attribute,
											   differentiate_atom_type)
	bond_attributes_dict = get_bond_attributes(molecule, 
											   non_H_atoms, 
											   add_extra_bond_attribute,
											   differentiate_bond_type)

	atom_attribute_num = len(atom_attributes_dict.values()[0])
	bond_attribute_num = len(bond_attributes_dict.values()[0])
	attribute_num = atom_attribute_num + bond_attribute_num
	molecule_tensor = np.zeros((non_H_atom_num, non_H_atom_num, attribute_num))

	
	for i, atom in enumerate(non_H_atoms):
		molecule_tensor[i, i, :] = np.concatenate((atom_attributes_dict[atom], np.zeros_like(bond_attributes_dict.values()[0])))
		for bonded_atom, bond in atom.bonds.iteritems():
			if bonded_atom in non_H_atoms:
				j = non_H_atoms.index(bonded_atom)
				molecule_tensor[i, j, :] = np.concatenate((atom_attributes_dict[bonded_atom], bond_attributes_dict[bond]))

	return molecule_tensor

def get_attribute_vector_size(add_extra_atom_attribute=True, 
							  add_extra_bond_attribute=True,
							  differentiate_atom_type=True,
							  differentiate_bond_type=True):

	"""
	this method examines current feature engineering setup
	returns attibute vector size on the fly using an example molecule
	"""
	from rmgpy.molecule.molecule import Molecule

	molecule = Molecule().fromSMILES('CC')
	non_H_atoms = [atom0 for atom0 in molecule.atoms if not atom0.isHydrogen()]
	
	atom_attributes_dict = get_atom_attributes(molecule, 
											   non_H_atoms, 
											   add_extra_atom_attribute,
											   differentiate_atom_type)
	bond_attributes_dict = get_bond_attributes(molecule, 
											   non_H_atoms, 
											   add_extra_bond_attribute,
											   differentiate_bond_type)

	atom_attribute_num = len(atom_attributes_dict.values()[0])
	bond_attribute_num = len(bond_attributes_dict.values()[0])
	attribute_num = atom_attribute_num + bond_attribute_num

	return attribute_num

def get_atom_attributes(molecule, non_H_atoms, 
						add_extra_attribute=True,
						differentiate_atom_type=True):

	"""
	this method takes a molecule with hydrogen pre-removed and returns a dict
	with non-H atom as key, atom attributes as value
	"""
	
	atom_attributes_dict = {}
	for atom in non_H_atoms:

		attributes = oneHotVector(
			atom.element.number, 
			[5, 6, 7, 8, 9, 15, 16, 17, 35, 53, 999]
		)

		non_H_neighbors = [atom0 for atom0 in atom.bonds if not atom0.isHydrogen()]
		attributes += oneHotVector(
			len(non_H_neighbors),
			[0, 1, 2, 3, 4, 5]
		)
		# Add hydrogen count
		H_neighbors = [atom0 for atom0 in atom.bonds if atom0.isHydrogen()]
		attributes += oneHotVector(
			len(H_neighbors),
			[0, 1, 2, 3, 4]
		)
		# charge
		atom.updateCharge()
		attributes.append(atom.charge)
		
		#in ring
		attributes.append(molecule.isVertexInCycle(atom))
		
		if differentiate_atom_type:
			# Add boolean if aromatic atom
			is_aromatic = False
			for bonded_atom, bond in atom.bonds.iteritems():
				if bond.isBenzene():
					is_aromatic = True
					break
			attributes.append(is_aromatic)

		# add atom in i-member rings
		if add_extra_attribute:
			attributes += is_atom_in_ring(molecule, atom)

		atom_attributes_dict[atom] = np.array(attributes, dtype=np.float32)
	
	return atom_attributes_dict

def is_atom_in_ring(molecule, atom):
	"""
	Given an atom and the molecule it belongs to
	return a list of counts, each count shows
	the occurances of the atom is involved in
	3, 4, 5, 6, 7, 8-member rings

	For instance, atom is in one 3-member ring and two 
	5-member rings, the returned list would be
	[1, 0, 2, 0, 0, 0]
	"""
	SSSR = molecule.getDeterministicSmallestSetOfSmallestRings()
	atom_in_rings = [0]*6
	for ring in SSSR:
		if atom in ring:
			idx = min(len(ring) - 3, 5)
			atom_in_rings[idx] = atom_in_rings[idx] + 1

	return atom_in_rings

def get_bond_attributes(molecule, 
						non_H_atoms, 
						add_extra_attribute=True, 
						differentiate_bond_type=True):
	"""
	this method takes a molecule with hydrogen pre-removed and returns a dict
	with bond as key, bond attributes as value
	"""
	bond_attributes_dict = {}
	for atom in non_H_atoms:
		for bonded_atom, bond in atom.bonds.iteritems():
			if not bonded_atom.isHydrogen() and (bond not in bond_attributes_dict):
				
				attributes = [1]
				if differentiate_bond_type:
					attributes.extend(oneHotVector(bond.getOrderStr(),
												   ['S', 'B', 'D', 'T']))
					# Add if is aromatic
					attributes.append(bond.isBenzene())

					attributes.append(is_bond_conjugated(bond))
				
				attributes.append(molecule.__isChainInCycle([bond.atom1, bond.atom2]))

				if add_extra_attribute:
					attributes.extend(is_bond_in_ring(molecule, bond))

				# the 1 here doesn't contribute to
				# fingerprint calculation; it's just
				# a helpful indicator in matrix operation
				# to show there's a bond in between some 
				# row_idx and col_idx
				attributes.append(1)

				bond_attributes_dict[bond] = np.array(attributes, dtype=np.float32)

	if not bond_attributes_dict:
		if differentiate_bond_type:
			bond_attributes_dict['no_bond'] = np.array([0]*15, dtype=np.float32)
		else:
			bond_attributes_dict['no_bond'] = np.array([0]*9, dtype=np.float32)
	
	return bond_attributes_dict

def is_bond_conjugated(bond):
	
	atom1 = bond.atom1
	atom2 = bond.atom2
	if bond.getOrderStr() == 'S':
		atom1_has_non_single_bond = False
		atom2_has_non_single_bond = False
		for bonded_atom, bond in atom1.bonds.iteritems():
			if bond.getOrderStr() != "S":
				atom1_has_non_single_bond = True
				break
		for bonded_atom, bond in atom2.bonds.iteritems():
			if bond.getOrderStr() != "S":
				atom2_has_non_single_bond = True
				break

		return (atom1_has_non_single_bond and atom2_has_non_single_bond)

	elif bond.getOrderStr() == 'B':
		return True
	else:
		for atom3, bond23 in atom2.bonds.iteritems():
			if bond23.getOrderStr() == "S":
				for atom4, bond34 in atom3.bonds.iteritems():
					if bond34.getOrderStr() != "S":
						return True

		for atom3, bond13 in atom1.bonds.iteritems():
			if bond13.getOrderStr() == "S":
				for atom4, bond34 in atom3.bonds.iteritems():
					if bond34.getOrderStr() != "S":
						return True

		return False

def is_bond_in_ring(molecule, bond):
	"""
	Given a bond and the molecule it belongs to
	return a list of counts, each count shows
	the occurances where the bond is involved in
	3, 4, 5, 6, 7, 8-member rings

	For instance, a bond is in one 3-member ring and two 
	5-member rings, the returned list would be
	[1, 0, 2, 0, 0, 0]
	"""
	SSSR = molecule.getDeterministicSmallestSetOfSmallestRings()
	bond_in_rings = [0]*6
	atom1 = bond.atom1
	atom2 = bond.atom2
	for ring in SSSR:
		if atom1 in ring and atom1 in ring :
			idx = min(len(ring) - 3, 5)
			bond_in_rings[idx] = bond_in_rings[idx] + 1

	return bond_in_rings
	
def oneHotVector(val, lst):
	'''Converts a value to a one-hot vector based on options in lst'''
	if val not in lst:
		val = lst[-1]
	return map(lambda x: x == val, lst)