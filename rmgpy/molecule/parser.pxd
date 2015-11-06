# global imports

cimport element as elements
cimport inchi as inchiutil

# no .pxd files for these:
#from .util cimport retrieveElementCount, VALENCES, ORDERS
#from .inchi cimport AugmentedInChI, compose_aug_inchi_key, compose_aug_inchi, INCHI_PREFIX, MULT_PREFIX, U_LAYER_PREFIX

from .molecule cimport Atom, Bond, Molecule

cpdef list BACKENDS
cpdef dict INSTALLED_BACKENDS
cpdef dict INCHI_LOOKUPS
cpdef dict SMILES_LOOKUPS


#  from <identifier> functions:

cdef Molecule __fromSMILES(Molecule mol, str smilesstr, str backend)

cdef Molecule __fromInChI(Molecule mol, str inchistr, str backend)

cdef Molecule __parse(Molecule mol, str identifier, str type_identifier, str backend)

cpdef Molecule parse_openbabel(Molecule mol, str identifier, str type_identifier)

cpdef Molecule fromInChI(Molecule mol, str inchistr, backend=*)

cpdef Molecule fromSMILES(Molecule mol, str smilesstr, str backend=*)

cpdef Molecule fromSMARTS(Molecule mol, str smartsstr)

cpdef Molecule fromAugmentedInChI(Molecule mol, aug_inchi)
    
cpdef Molecule fromRDKitMol(Molecule mol, object rdkitmol)

cpdef Molecule fromOBMol(Molecule mol, object obmol)

cdef Molecule __lookup(Molecule mol, str identifier, str type_identifier)

# parser helper functions: 

cpdef reset_lone_pairs(Molecule mol, list p_indices)

cdef Molecule fix_unsaturated_bond_to_biradical(Molecule mol, str inchi, list u_indices)

cpdef bint isUnsaturated(Molecule mol)

cpdef isCorrectlyParsed(Molecule mol, str identifier)
   
cpdef check(Molecule mol, aug_inchi)

cpdef fix_oxygen_unsaturated_bond(Molecule mol, list u_indices)

cpdef find_lowest_u_layer(Molecule mol, list u_layer, list equivalent_atoms)

cpdef fixCharge(Molecule mol, list u_indices)

cpdef fix_triplet_to_singlet(Molecule mol, list p_indices)
