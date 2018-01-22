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

cpdef dict MOLECULE_LOOKUPS
cpdef dict RADICAL_LOOKUPS

cpdef str toInChI(Molecule mol, str backend=?, int aug_level=?)

cpdef str toInChIKey(Molecule mol, str backend=?, int aug_level=?)

cpdef str toSMARTS(Molecule mol, backend=?)

cpdef str toSMILES(Molecule mol, backend=?)

cpdef Molecule fromInChI(Molecule mol, str inchistr, backend=?)

cpdef Molecule fromSMILES(Molecule mol, str smilesstr, str backend=?)

cpdef Molecule fromSMARTS(Molecule mol, str smartsstr, str backend=?)

cpdef Molecule fromAugmentedInChI(Molecule mol, aug_inchi)

cpdef object _rdkit_translator(object input_object, str identifier_type, Molecule mol=?)

cpdef object _openbabel_translator(object input_object, str identifier_type, Molecule mol=?)

cdef Molecule _lookup(Molecule mol, str identifier, str identifier_type)

cpdef _is_correctly_parsed(Molecule mol, str identifier)

cdef Molecule _read(Molecule mol, str identifier, str identifier_type, str backend)

cdef str _write(Molecule mol, str identifier_type, str backend)
