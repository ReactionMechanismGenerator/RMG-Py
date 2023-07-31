###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2023 Prof. William H. Green (whgreen@mit.edu),           #
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

from rmgpy.molecule.graph cimport Vertex, Edge, Graph
from rmgpy.molecule.molecule cimport Atom, Bond, Molecule

cpdef list populate_resonance_algorithms(dict features=?)

cpdef dict analyze_molecule(Graph mol, bint save_order=?)

cpdef list generate_resonance_structures(Graph mol, bint clar_structures=?, bint keep_isomorphic=?, bint filter_structures=?, bint save_order=?)

cpdef list _generate_resonance_structures(list mol_list, list method_list, bint keep_isomorphic=?, bint copy=?, bint save_order=?)

cpdef list generate_allyl_delocalization_resonance_structures(Graph mol)

cpdef list generate_lone_pair_multiple_bond_resonance_structures(Graph mol)

cpdef list generate_adj_lone_pair_radical_resonance_structures(Graph mol)

cpdef list generate_adj_lone_pair_multiple_bond_resonance_structures(Graph mol)

cpdef list generate_adj_lone_pair_radical_multiple_bond_resonance_structures(Graph mol)

cpdef list generate_N5dc_radical_resonance_structures(Graph mol)

cpdef list generate_isomorphic_resonance_structures(Graph mol, bint saturate_h=?)

cpdef list generate_optimal_aromatic_resonance_structures(Graph mol, dict features=?, bint save_order=?)

cpdef list generate_aromatic_resonance_structure(Graph mol, list aromatic_bonds=?, bint copy=?, bint save_order=?)

cpdef list generate_aryne_resonance_structures(Graph mol)

cpdef list generate_kekule_structure(Graph mol)

cpdef list generate_clar_structures(Graph mol, bint save_order=?)

cpdef list _clar_optimization(Graph mol, list constraints=?, max_num=?, save_order=?)

cpdef list _clar_transformation(Graph mol, list aromatic_ring)

cpdef list generate_adsorbate_shift_down_resonance_structures(Graph mol)

cpdef list generate_adsorbate_shift_up_resonance_structures(Graph mol)

cpdef list generate_adsorbate_double_shift_up_resonance_structures(Graph mol)

cpdef list generate_adsorbate_double_shift_down_resonance_structures(Graph mol)

cpdef list generate_adsorbate_conjugate_resonance_structures(Graph mol)