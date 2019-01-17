from .molecule cimport Atom, Molecule

cpdef dict _known_smiles_molecules
cpdef dict _known_smiles_radicals

cpdef str toInChI(Molecule mol)

cpdef str create_U_layer(Molecule mol, str auxinfo)

cpdef str toAugmentedInChI(Molecule mol)

cpdef str toInChIKey(Molecule mol)

cpdef str toAugmentedInChIKey(Molecule mol)

cpdef str toSMARTS(Molecule mol)

cpdef str toSMILES(Molecule mol)

cpdef toOBMol(Molecule mol, bint returnMapping=*)

cpdef toRDKitMol(Molecule mol, bint removeHs=*, bint returnMapping=*, bint sanitize=*)

cpdef bint is_valid_combo(list combo, Molecule mol, list distances)

cpdef list find_lowest_u_layer(Molecule mol, list u_layer, list equivalent_atoms)

cpdef Molecule generate_minimum_resonance_isomer(Molecule mol)

cpdef list get_unpaired_electrons(Molecule mol)

cpdef list compute_agglomerate_distance(list agglomerates, Molecule mol)

cpdef str create_P_layer(Molecule mol, str auxinfo)