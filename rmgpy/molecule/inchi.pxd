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

from .molecule cimport Atom, Bond, Molecule

cpdef tuple decompose_aug_inchi(str string)

cpdef str remove_inchi_prefix(str string)

cpdef str compose_aug_inchi(str inchi, str ulayer=*, str player=*)

cpdef str compose_aug_inchi_key(str inchi_key, str ulayer=*, str player=*)

cpdef list _parse_h_layer(str inchi)

cpdef list _parse_e_layer(str auxinfo)

cpdef list _parse_n_layer(str auxinfo)

cpdef bint _has_unexpected_lone_pairs(Molecule mol)

cpdef list _get_unpaired_electrons(Molecule mol)

cpdef Molecule _generate_minimum_resonance_isomer(Molecule mol)

cpdef list _compute_agglomerate_distance(list agglomerates, Molecule mol)

cpdef bint _is_valid_combo(list combo, Molecule mol, list distances)

cpdef list _find_lowest_u_layer(Molecule mol, list u_layer, list equivalent_atoms)

cpdef str _create_u_layer(Molecule mol, str auxinfo)

cpdef list _find_lowest_p_layer(Molecule minmol, list p_layer, list equivalent_atoms)

cpdef str _create_p_layer(Molecule mol, str auxinfo)

cpdef tuple create_augmented_layers(Molecule mol)


cpdef _fix_triplet_to_singlet(Molecule mol, list p_indices)

cpdef _convert_charge_to_unpaired_electron(Molecule mol, list u_indices)

cpdef _convert_4_atom_3_bond_path(Atom start)

cpdef _convert_3_atom_2_bond_path(Atom start, Molecule mol)

cpdef _convert_delocalized_charge_to_unpaired_electron(Molecule mol, list u_indices)

cpdef _fix_adjacent_charges(Molecule mol)

cpdef _fix_charge(Molecule mol, list u_indices)

cpdef _reset_lone_pairs(Molecule mol, list p_indices)

cpdef _fix_oxygen_unsaturated_bond(Molecule mol, list u_indices)

cpdef bint _is_unsaturated(Molecule mol)

cpdef bint _convert_unsaturated_bond_to_triplet(Bond bond)

cpdef bint _fix_mobile_h(Molecule mol, str inchi, int u1, int u2)

cpdef bint _fix_butadiene_path(Atom start, Atom end)

cpdef Molecule _fix_unsaturated_bond_to_biradical(Molecule mol, str inchi, list u_indices)

cpdef _fix_unsaturated_bond(Molecule mol, list u_indices, aug_inchi)

cpdef _check_molecule(Molecule mol, aug_inchi)

cpdef fix_molecule(Molecule mol, aug_inchi)
