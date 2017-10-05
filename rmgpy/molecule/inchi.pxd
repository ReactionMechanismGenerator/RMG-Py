###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2018 Prof. William H. Green (whgreen@mit.edu),           #
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

cpdef tuple decompose(string)

cpdef str ignore_prefix(str string)

cpdef str compose_aug_inchi(str inchi, str ulayer=*, str player=*)

cpdef str compose_aug_inchi_key(str inchi_key, str ulayer=*, str player=*)

cpdef list parse_H_layer(str inchi)

cpdef list parse_E_layer(str auxinfo)

cpdef list parse_N_layer(str auxinfo)

cpdef str create_U_layer(Molecule mol, str auxinfo)

cpdef bint is_valid_combo(list combo, Molecule mol, list distances)

cpdef list find_lowest_u_layer(Molecule mol, list u_layer, list equivalent_atoms)

cpdef Molecule generate_minimum_resonance_isomer(Molecule mol)

cpdef list get_unpaired_electrons(Molecule mol)

cpdef list compute_agglomerate_distance(list agglomerates, Molecule mol)

cpdef str create_P_layer(Molecule mol, str auxinfo)

cpdef reset_lone_pairs(Molecule mol, list p_indices)

cdef Molecule fix_unsaturated_bond_to_biradical(Molecule mol, str inchi, list u_indices)

cpdef bint isUnsaturated(Molecule mol)

cpdef check(Molecule mol, aug_inchi)

cpdef fix_oxygen_unsaturated_bond(Molecule mol, list u_indices)

cpdef fixCharge(Molecule mol, list u_indices)

cpdef fix_triplet_to_singlet(Molecule mol, list p_indices)
