###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2020 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

cimport numpy as np

cimport rmgpy.constants as constants
cimport rmgpy.molecule.group as gr
from rmgpy.molecule.atomtype cimport AtomType
from rmgpy.molecule.element cimport Element
from rmgpy.molecule.graph cimport Vertex, Edge, Graph

################################################################################
cdef dict bond_orders 

cdef class Atom(Vertex):

    cdef public Element element
    cdef public short radical_electrons
    cdef public short charge
    cdef public str label
    cdef public AtomType atomtype
    cdef public np.ndarray coords
    cdef public short lone_pairs
    cdef public int id
    cdef public dict props
    
    cpdef bint equivalent(self, Vertex other, bint strict=?) except -2

    cpdef bint is_specific_case_of(self, Vertex other) except -2

    cpdef Vertex copy(self)

    cpdef bint is_hydrogen(self)

    cpdef bint is_non_hydrogen(self)

    cpdef bint is_carbon(self)

    cpdef bint is_oxygen(self)

    cpdef bint is_fluorine(self)

    cpdef bint is_silicon(self)

    cpdef bint is_phosphorus(self)

    cpdef bint is_sulfur(self)

    cpdef bint is_chlorine(self)

    cpdef bint is_iodine(self)

    cpdef bint is_nos(self)
    
    cpdef bint is_surface_site(self)
    
    cpdef increment_radical(self)

    cpdef decrement_radical(self)
    
    cpdef set_lone_pairs(self, int lone_pairs)
    
    cpdef increment_lone_pairs(self)
    
    cpdef decrement_lone_pairs(self)
    
    cpdef update_charge(self)

    cpdef get_total_bond_order(self)
    
################################################################################
    
cdef class Bond(Edge):

    cdef public float order

    cpdef bint equivalent(self, Edge other) except -2

    cpdef bint is_specific_case_of(self, Edge other) except -2

    cpdef str get_order_str(self)
    
    cpdef set_order_str(self, str new_order)
    
    cpdef float get_order_num(self)
    
    cpdef set_order_num(self, float new_order)

    cpdef Edge copy(self)
    
    cpdef bint is_order(self, float other_order)

    cpdef bint is_van_der_waals(self) except -2

    cpdef bint is_single(self) except -2

    cpdef bint is_double(self) except -2

    cpdef bint is_triple(self) except -2
    
    cpdef bint is_quadruple(self) except -2
    
    cpdef bint is_benzene(self) except -2

    cpdef increment_order(self)

    cpdef decrement_order(self)

    cpdef str get_bond_string(self)

################################################################################

cdef class Molecule(Graph):

    cdef public float symmetry_number
    cdef public int multiplicity
    cdef public bint reactive
    cdef public dict props
    cdef str _fingerprint
    cdef str _inchi
    cdef str _smiles

    cpdef add_atom(self, Atom atom)

    cpdef add_bond(self, Bond bond)

    cpdef dict get_bonds(self, Atom atom)

    cpdef Bond get_bond(self, Atom atom1, Atom atom2)

    cpdef bint has_atom(self, Atom atom)

    cpdef bint has_bond(self, Atom atom1, Atom atom2)

    cpdef bint contains_surface_site(self)
    
    cpdef bint is_surface_site(self)

    cpdef remove_atom(self, Atom atom)

    cpdef remove_bond(self, Bond bond)

    cpdef remove_van_der_waals_bonds(self)

    cpdef sort_atoms(self)
    
    cpdef str get_formula(self)

    cpdef short get_radical_count(self)

    cpdef short get_singlet_carbene_count(self)

    cpdef double get_molecular_weight(self)

    cpdef int get_num_atoms(self, str element=?)

    cpdef Graph copy(self, bint deep=?)

    cpdef delete_hydrogens(self)

    cpdef clear_labeled_atoms(self)

    cpdef bint contains_labeled_atom(self, str label) except -2

    cpdef list get_labeled_atoms(self, str label)

    cpdef dict get_all_labeled_atoms(self)

    cpdef dict get_element_count(self)

    cpdef bint is_isomorphic(self, Graph other, dict initial_map=?, bint generate_initial_map=?, bint save_order=?, bint strict=?) except -2

    cpdef list find_isomorphism(self, Graph other, dict initial_map=?, bint save_order=?, bint strict=?)

    cpdef bint is_subgraph_isomorphic(self, Graph other, dict initial_map=?, bint generate_initial_map=?, bint save_order=?) except -2

    cpdef list find_subgraph_isomorphisms(self, Graph other, dict initial_map=?, bint save_order=?)

    cpdef bint is_atom_in_cycle(self, Atom atom) except -2

    cpdef bint is_bond_in_cycle(self, Bond bond) except -2

    cpdef draw(self, str path)

    cpdef from_inchi(self, str inchistr, backend=?, bint raise_atomtype_exception=?)

    cpdef from_smiles(self, str smilesstr, backend=?, bint raise_atomtype_exception=?)

    cpdef from_adjacency_list(self, str adjlist, bint saturate_h=?, bint raise_atomtype_exception=?,
                              bint raise_charge_exception=?)

    cpdef from_xyz(self, np.ndarray atomic_nums, np.ndarray coordinates, float critical_distance_factor=?, bint raise_atomtype_exception=?)
    
    cpdef str to_inchi(self)

    cpdef str to_augmented_inchi(self)

    cpdef str to_inchi_key(self)

    cpdef str to_augmented_inchi_key(self)

    cpdef str to_smiles(self)

    cpdef to_adjacency_list(self, str label=?, bint remove_h=?, bint remove_lone_pairs=?, bint old_style=?)

    cpdef bint is_linear(self) except -2

    cpdef bint is_heterocyclic(self) except -2

    cpdef int count_internal_rotors(self) except -2

    cpdef double calculate_cp0(self) except -1

    cpdef double calculate_cpinf(self) except -1
    
    cpdef update_atomtypes(self, bint log_species=?, bint raise_exception=?)
    
    cpdef bint is_radical(self) except -2

    cpdef bint has_lone_pairs(self) except -2

    cpdef bint is_aryl_radical(self, list aromatic_rings=?) except -2

    cpdef float calculate_symmetry_number(self) except -1

    cpdef list generate_resonance_structures(self, bint keep_isomorphic=?, bint filter_structures=?)

    cpdef identify_ring_membership(self)

    cpdef int count_aromatic_rings(self)

    cpdef tuple get_aromatic_rings(self, list rings=?)

    cpdef list get_deterministic_sssr(self)

    cpdef kekulize(self)

    cpdef assign_atom_ids(self)

    cpdef bint atom_ids_valid(self)

    cpdef bint is_identical(self, Molecule other, bint strict=?) except -2

    cpdef dict enumerate_bonds(self)

cdef atom_id_counter
