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

from .graph cimport Vertex, Edge, Graph
from .atomtype cimport AtomType
cimport rmgpy.molecule.molecule as mol

################################################################################

cdef class GroupAtom(Vertex):

    cdef public list atomType
    cdef public list radicalElectrons
    cdef public list charge
    cdef public str label
    cdef public list lonePairs
    cdef public dict props

    cpdef Vertex copy(self)

    cpdef __changeBond(self, short order)

    cpdef __formBond(self, float order)

    cpdef __breakBond(self, float order)

    cpdef __gainRadical(self, short radical)

    cpdef __loseRadical(self, short radical)
    
    cpdef __gainPair(self, short radical)

    cpdef __losePair(self, short radical)

    cpdef applyAction(self, list action)

    cpdef bint equivalent(self, Vertex other) except -2

    cpdef bint isSpecificCaseOf(self, Vertex other) except -2

    cpdef bint isOxygen(self)

    cpdef bint isSulfur(self)

    cpdef list countBonds(self, wildcards = ?)

    cpdef bint hasWildcards(self)

    cpdef mol.Atom makeSampleAtom(self)

################################################################################

cdef class GroupBond(Edge):

    cdef public list order

    cpdef Edge copy(self)

    cpdef list getOrderStr(self)
    
    cpdef setOrderStr(self, list newOrder)
    
    cpdef list getOrderNum(self)
    
    cpdef setOrderNum(self, list newOrder)

    cpdef __changeBond(self, short order)

    cpdef bint isSingle(self, bint wildcards = ?) except -2

    cpdef bint isDouble(self, bint wildcards = ?) except -2

    cpdef bint isTriple(self, bint wildcards = ?) except -2

    cpdef bint isBenzene(self, bint wildcards = ?) except -2

    cpdef applyAction(self, list action)

    cpdef bint equivalent(self, Edge other) except -2

    cpdef bint isSpecificCaseOf(self, Edge other) except -2

    cpdef makeBond(self, mol.Molecule molecule, mol.Atom atom1, mol.Atom atom2)

################################################################################

cdef class Group(Graph):

    cdef public dict props
    cdef public list multiplicity

    # These read-only attribues act as a "fingerprint" for accelerating
    # subgraph isomorphism checks
    cdef public dict elementCount
    cdef public short radicalCount

    cpdef addAtom(self, GroupAtom atom)

    cpdef addBond(self, GroupBond bond)

    cpdef dict getBonds(self, GroupAtom atom)

    cpdef GroupBond getBond(self, GroupAtom atom1, GroupAtom atom2)

    cpdef bint hasAtom(self, GroupAtom atom)

    cpdef bint hasBond(self, GroupAtom atom1, GroupAtom atom2)

    cpdef removeAtom(self, GroupAtom atom)

    cpdef removeBond(self, GroupBond bond)

    cpdef sortAtoms(self)

    cpdef list sortByConnectivity(self, list atomList)

    cpdef Graph copy(self, bint deep=?)

    cpdef clearLabeledAtoms(self)

    cpdef bint containsLabeledAtom(self, str label)

    cpdef GroupAtom getLabeledAtom(self, str label)

    cpdef dict getLabeledAtoms(self)

    cpdef dict get_element_count(self)

    cpdef fromAdjacencyList(self, str adjlist)

    cpdef toAdjacencyList(self, str label=?)
    
    cpdef updateFingerprint(self)

    cpdef update_charge(self)

    cpdef bint isIsomorphic(self, Graph other, dict initialMap=?) except -2

    cpdef list findIsomorphism(self, Graph other, dict initialMap=?)

    cpdef bint isSubgraphIsomorphic(self, Graph other, dict initialMap=?) except -2

    cpdef list findSubgraphIsomorphisms(self, Graph other, dict initialMap=?)
    
    cpdef bint isIdentical(self, Graph other)

    cpdef bint isAromaticRing(self)

    cpdef bint standardizeAtomType(self)

    cpdef bint addExplicitLigands(self)

    cpdef GroupAtom createAndConnectAtom(self, list atomtype, GroupAtom connectingAtom, list bondOrders)

    cpdef bint standardizeGroup(self)

    cpdef Group addImplicitAtomsFromAtomType(self)

    cpdef pickWildcards(self)

    cpdef mol.Molecule makeSampleMolecule(self)

    cpdef tuple classifyBenzeneCarbons(self, dict partners=?)

    cpdef Group addImplicitBenzene(self)

    cpdef bint isBenzeneExplicit(self)

    cpdef Group mergeGroups(self, Group other)

    cpdef resetRingMembership(self)
