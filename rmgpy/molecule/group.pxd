################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2009-2011 by the RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

from .graph cimport Vertex, Edge, Graph
from .atomtype cimport AtomType

################################################################################

cdef class GroupAtom(Vertex):

    cdef public list atomType
    cdef public list radicalElectrons
    cdef public list charge
    cdef public str label
    cdef public list lonePairs

    cpdef Vertex copy(self)

    cpdef __changeBond(self, short order)

    cpdef __formBond(self, str order)

    cpdef __breakBond(self, str order)

    cpdef __gainRadical(self, short radical)

    cpdef __loseRadical(self, short radical)
    
    cpdef __gainPair(self, short radical)

    cpdef __losePair(self, short radical)

    cpdef applyAction(self, list action)

    cpdef bint equivalent(self, Vertex other)

    cpdef bint isSpecificCaseOf(self, Vertex other)

################################################################################

cdef class GroupBond(Edge):

    cdef public list order

    cpdef Edge copy(self)

    cpdef __changeBond(self, short order)

    cpdef applyAction(self, list action)

    cpdef bint equivalent(self, Edge other)

    cpdef bint isSpecificCaseOf(self, Edge other)

################################################################################

cdef class Group(Graph):

    cdef public list multiplicity

    # These read-only attribues act as a "fingerprint" for accelerating
    # subgraph isomorphism checks
    cdef public short carbonCount
    cdef public short nitrogenCount
    cdef public short oxygenCount
    cdef public short sulfurCount
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

    cpdef Graph copy(self, bint deep=?)

    cpdef clearLabeledAtoms(self)

    cpdef bint containsLabeledAtom(self, str label)

    cpdef GroupAtom getLabeledAtom(self, str label)

    cpdef dict getLabeledAtoms(self)

    cpdef fromAdjacencyList(self, str adjlist)

    cpdef toAdjacencyList(self, str label=?)
    
    cpdef updateFingerprint(self)

    cpdef bint isIsomorphic(self, Graph other, dict initialMap=?) except -2

    cpdef list findIsomorphism(self, Graph other, dict initialMap=?)

    cpdef bint isSubgraphIsomorphic(self, Graph other, dict initialMap=?) except -2

    cpdef list findSubgraphIsomorphisms(self, Graph other, dict initialMap=?)
    
    cpdef bint isIdentical(self, Graph other)
