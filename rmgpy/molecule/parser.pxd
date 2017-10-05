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

cdef Molecule __fromSMARTS(Molecule mol, str identifier, str backend)

cdef Molecule __parse(Molecule mol, str identifier, str type_identifier, str backend)

cpdef Molecule parse_openbabel(Molecule mol, str identifier, str type_identifier)

cpdef Molecule fromInChI(Molecule mol, str inchistr, backend=*)

cpdef Molecule fromSMILES(Molecule mol, str smilesstr, str backend=*)

cpdef Molecule fromSMARTS(Molecule mol, str smartsstr, str backend=*)

cpdef Molecule fromAugmentedInChI(Molecule mol, aug_inchi)
    
cdef Molecule __lookup(Molecule mol, str identifier, str type_identifier)

# parser helper functions: 

cpdef reset_lone_pairs(Molecule mol, list p_indices)

cdef Molecule fix_unsaturated_bond_to_biradical(Molecule mol, str inchi, list u_indices)

cpdef bint isUnsaturated(Molecule mol)

cpdef isCorrectlyParsed(Molecule mol, str identifier)
   
cpdef check(Molecule mol, aug_inchi)

cpdef fix_oxygen_unsaturated_bond(Molecule mol, list u_indices)

cpdef fixCharge(Molecule mol, list u_indices)

cpdef fix_triplet_to_singlet(Molecule mol, list p_indices)
