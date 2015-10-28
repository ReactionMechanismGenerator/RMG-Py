# global imports

cimport element as elements

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

cdef __parse(Molecule mol, str identifier, str type_identifier, str backend)

cpdef parse_openbabel(Molecule mol, str identifier, str type_identifier)

cpdef fromInChI(Molecule mol, str inchistr, backend=*)

cpdef fromSMILES(Molecule mol, str smilesstr, str backend=*)

cpdef fromSMARTS(Molecule mol, str smartsstr)

cpdef Molecule fromAugmentedInChI(Molecule mol, aug_inchi)
    
cpdef Molecule fromRDKitMol(Molecule mol, object rdkitmol)

cpdef fromOBMol(Molecule mol, object obmol)

cdef __lookup(Molecule mol, str identifier, str type_identifier)

# parser helper functions: 

cpdef reset_lone_pairs_to_default(Molecule mol)

cdef Molecule convert_unsaturated_bond_to_biradical(Molecule mol, str inchi, list u_indices)

cpdef bint isUnsaturated(Molecule mol)
    
cpdef bint check_number_unpaired_electrons(Molecule mol)

cpdef isCorrectlyParsed(Molecule mol, str identifier)
   
cpdef check(Molecule mol, aug_inchi)

cpdef correct_O_unsaturated_bond(Molecule mol, list u_indices)

cpdef find_lowest_u_layer(Molecule mol, list u_layer, list equivalent_atoms)

cpdef fixCharge(Molecule mol, list u_indices)
