###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2021 Prof. William H. Green (whgreen@mit.edu),           #
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

from .graph cimport Vertex, Edge, Graph
from .atomtype cimport AtomType
cimport rmgpy.molecule.molecule as mol
from cpython cimport bool
################################################################################

cdef class GroupAtom(Vertex):

    cdef public list atomtype
    cdef public list radical_electrons
    cdef public list charge
    cdef public str label
    cdef public list lone_pairs

    cdef public dict props

    cdef public list reg_dim_atm
    cdef public list reg_dim_u
    cdef public list reg_dim_r

    cpdef Vertex copy(self)

    cpdef _change_bond(self, short order)

    cpdef _form_bond(self, float order)

    cpdef _break_bond(self, float order)

    cpdef _gain_radical(self, short radical)

    cpdef _lose_radical(self, short radical)
    
    cpdef _gain_pair(self, short radical)

    cpdef _lose_pair(self, short radical)

    cpdef apply_action(self, list action)

    cpdef bint equivalent(self, Vertex other, bint strict=?) except -2

    cpdef bint is_specific_case_of(self, Vertex other) except -2

    cpdef bint is_surface_site(self) except -2

    cpdef bint is_bonded_to_surface(self) except -2

    cpdef bint is_oxygen(self)

    cpdef bint is_sulfur(self)

    cpdef bint is_nitrogen(self)

    cpdef bint is_carbon(self)

    cpdef list count_bonds(self, wildcards = ?)

    cpdef bint has_wildcards(self)

    cpdef mol.Atom make_sample_atom(self)

################################################################################

cdef class GroupBond(Edge):

    cdef public list order
    cdef public list reg_dim

    cpdef Edge copy(self)

    cpdef list get_order_str(self)
    
    cpdef set_order_str(self, list new_order)
    
    cpdef list get_order_num(self)
    
    cpdef set_order_num(self, list new_order)

    cpdef _change_bond(self, short order)

    cpdef bint is_single(self, bint wildcards = ?) except -2

    cpdef bint is_double(self, bint wildcards = ?) except -2

    cpdef bint is_triple(self, bint wildcards = ?) except -2

    cpdef bint is_benzene(self, bint wildcards = ?) except -2

    cpdef apply_action(self, list action)

    cpdef bint equivalent(self, Edge other) except -2

    cpdef bint is_specific_case_of(self, Edge other) except -2

    cpdef make_bond(self, mol.Molecule molecule, mol.Atom atom1, mol.Atom atom2)

################################################################################

cdef class Group(Graph):

    cdef public dict props
    cdef public list multiplicity

    # These read-only attribues act as a "fingerprint" for accelerating
    # subgraph isomorphism checks
    cdef public dict elementCount
    cdef public short radicalCount

    cpdef add_atom(self, GroupAtom atom)

    cpdef add_bond(self, GroupBond bond)

    cpdef dict get_bonds(self, GroupAtom atom)

    cpdef GroupBond get_bond(self, GroupAtom atom1, GroupAtom atom2)

    cpdef bint has_atom(self, GroupAtom atom)

    cpdef bint has_bond(self, GroupAtom atom1, GroupAtom atom2)

    cpdef remove_atom(self, GroupAtom atom)

    cpdef remove_bond(self, GroupBond bond)

    cpdef remove_van_der_waals_bonds(self)

    cpdef sort_atoms(self)

    cpdef list sort_by_connectivity(self, list atom_list)

    cpdef Graph copy(self, bint deep=?)

    cpdef clear_labeled_atoms(self)

    cpdef bint contains_labeled_atom(self, str label)

    cpdef list get_labeled_atoms(self, str label)

    cpdef dict get_all_labeled_atoms(self)

    cpdef dict get_element_count(self)

    cpdef from_adjacency_list(self, str adjlist)

    cpdef to_adjacency_list(self, str label=?)
    
    cpdef update_fingerprint(self)

    cpdef update_charge(self)

    cpdef bint is_isomorphic(self, Graph other, dict initial_map=?, bint generate_initial_map=?, bint save_order=?, bint strict=?) except -2

    cpdef list find_isomorphism(self, Graph other, dict initial_map=?, bint save_order=?, bint strict=?)

    cpdef bint is_subgraph_isomorphic(self, Graph other, dict initial_map=?, bint generate_initial_map=?, bint save_order=?) except -2

    cpdef list find_subgraph_isomorphisms(self, Graph other, dict initial_map=?, bint save_order=?)
    
    cpdef bint is_identical(self, Graph other, bint save_order=?)

    cpdef bint is_surface_site(self) except -2

    cpdef bint contains_surface_site(self) except -2

    cpdef list get_surface_sites(self)

    cpdef bint is_aromatic_ring(self)

    cpdef bint standardize_atomtype(self)

    cpdef bint add_explicit_ligands(self)

    cpdef GroupAtom create_and_connect_atom(self, list atomtype, GroupAtom connecting_atom, list bond_orders)

    cpdef bint standardize_group(self)

    cpdef Group add_implicit_atoms_from_atomtype(self)

    cpdef pick_wildcards(self)

    cpdef mol.Molecule make_sample_molecule(self)

    cpdef tuple classify_benzene_carbons(self, dict partners=?)

    cpdef Group add_implicit_benzene(self)

    cpdef bint is_benzene_explicit(self)

    cpdef Group merge_groups(self, Group other, bint keep_identical_labels=?)

    cpdef reset_ring_membership(self)
