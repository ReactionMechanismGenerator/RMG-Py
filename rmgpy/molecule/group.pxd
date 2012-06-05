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

################################################################################

cdef class GroupAtom(Vertex):

    cdef public list atomType
    cdef public list radicalElectrons
    cdef public list spinMultiplicity
    cdef public list charge
    cdef public str label

    cpdef copy(self)

    cpdef __changeBond(self, short order)

    cpdef __formBond(self, str order)

    cpdef __breakBond(self, str order)

    cpdef __gainRadical(self, short radical)

    cpdef __loseRadical(self, short radical)

    cpdef applyAction(self, list action)

    cpdef bint equivalent(self, Vertex other)

    cpdef bint isSpecificCaseOf(self, Vertex other)

################################################################################

cdef class GroupBond(Edge):

    cdef public list order

    cpdef copy(self)

    cpdef __changeBond(self, short order)

    cpdef applyAction(self, list action)

    cpdef bint equivalent(self, Edge other)

    cpdef bint isSpecificCaseOf(self, Edge other)

################################################################################

cdef class Group(Graph):

    cpdef addAtom(self, GroupAtom atom)

    cpdef addBond(self, GroupAtom atom1, GroupAtom atom2, GroupBond bond)

    cpdef dict getBonds(self, GroupAtom atom)

    cpdef GroupBond getBond(self, GroupAtom atom1, GroupAtom atom2)

    cpdef bint hasAtom(self, GroupAtom atom)

    cpdef bint hasBond(self, GroupAtom atom1, GroupAtom atom2)

    cpdef removeAtom(self, GroupAtom atom)

    cpdef removeBond(self, GroupAtom atom1, GroupAtom GroupAtom2)

    cpdef sortAtoms(self)

    cpdef Graph copy(self, bint deep=?)

    cpdef clearLabeledAtoms(self)

    cpdef bint containsLabeledAtom(self, str label)

    cpdef GroupAtom getLabeledAtom(self, str label)

    cpdef dict getLabeledAtoms(self)

    cpdef fromAdjacencyList(self, str adjlist)

    cpdef toAdjacencyList(self, str label=?)

    cpdef bint isIsomorphic(self, Graph other, dict initialMap=?)

    cpdef tuple findIsomorphism(self, Graph other, dict initialMap=?)

    cpdef bint isSubgraphIsomorphic(self, Graph other, dict initialMap=?)

    cpdef tuple findSubgraphIsomorphisms(self, Graph other, dict initialMap=?)
