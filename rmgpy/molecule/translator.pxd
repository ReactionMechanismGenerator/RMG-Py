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
cimport element as elements
cimport inchi as inchiutil

cpdef list BACKENDS
cpdef dict INCHI_LOOKUPS
cpdef dict SMILES_LOOKUPS

cpdef dict _known_smiles_molecules
cpdef dict _known_smiles_radicals

cpdef str toInChI(Molecule mol)

cpdef str toAugmentedInChI(Molecule mol)

cpdef str toInChIKey(Molecule mol)

cpdef str toAugmentedInChIKey(Molecule mol)

cpdef str toSMARTS(Molecule mol)

cpdef str toSMILES(Molecule mol)

cdef Molecule __fromSMILES(Molecule mol, str smilesstr, str backend)

cdef Molecule __fromInChI(Molecule mol, str inchistr, str backend)

cdef Molecule __fromSMARTS(Molecule mol, str identifier, str backend)

cdef Molecule __parse(Molecule mol, str identifier, str type_identifier, str backend)

cpdef Molecule parse_openbabel(Molecule mol, str identifier, str type_identifier)

cpdef isCorrectlyParsed(Molecule mol, str identifier)

cpdef Molecule fromInChI(Molecule mol, str inchistr, backend=*)

cpdef Molecule fromSMILES(Molecule mol, str smilesstr, str backend=*)

cpdef Molecule fromSMARTS(Molecule mol, str smartsstr, str backend=*)

cpdef Molecule fromAugmentedInChI(Molecule mol, aug_inchi)

cdef Molecule __lookup(Molecule mol, str identifier, str type_identifier)

