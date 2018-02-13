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