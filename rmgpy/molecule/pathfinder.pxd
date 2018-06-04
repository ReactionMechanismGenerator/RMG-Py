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
from .graph cimport Vertex, Edge

cpdef list find_butadiene(Atom start, Atom end)

cpdef list find_butadiene_end_with_charge(Atom start)

cpdef list find_allyl_end_with_charge(Atom start)

cpdef list find_shortest_path(Vertex start, Vertex end, list path=*)

cpdef list add_unsaturated_bonds(list path)

cpdef list add_allyls(list path)

cpdef list add_inverse_allyls(list path)

cpdef dict compute_atom_distance(list atom_indices, Molecule mol)

cpdef list find_allyl_delocalization_paths(Atom atom1)

cpdef list find_lone_pair_radical_delocalization_paths(Atom atom1)

cpdef list find_lone_pair_multiple_bond_delocalization_paths(Atom atom1)

cpdef list find_lone_pair_radical_multiple_bond_delocalization_paths(Atom atom1)

cpdef list find_N5ddc_N5tc_delocalization_paths(Atom atom1)

cpdef list find_N5dc_delocalization_paths(Atom atom1)

cpdef bint is_NOS_able_to_gain_lone_pair(Atom atom)

cpdef bint is_NOS_able_to_lose_lone_pair(Atom atom)
