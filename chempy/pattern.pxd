################################################################################
#
#   ChemPy - A chemistry toolkit for Python
#
#   Copyright (c) 2010 by Joshua W. Allen (jwallen@mit.edu)
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

from graph cimport Vertex, Edge, Graph

################################################################################

cpdef bint atomTypesEquivalent(str atomType1, str atomType2)

cpdef bint atomTypesSpecificCaseOf(str atomType1, str atomType2)

cpdef str getAtomType(atom, dict bonds)

################################################################################

cdef class AtomPattern(Vertex):

    cdef public list atomType
    cdef public list radicalElectrons
    cdef public list spinMultiplicity
    cdef public list implicitHydrogens
    cdef public list charge
    cdef public str label

    cpdef copy(self)

    cpdef __changeBond(self, short order)

    cpdef __formBond(self, short order)

    cpdef __breakBond(self, short order)

    cpdef __gainRadical(self, short radical)

    cpdef __loseRadical(self, short radical)

    cpdef applyAction(self, tuple action)

    cpdef bint equivalent(self, Vertex other)

    cpdef bint isSpecificCaseOf(self, Vertex other)

################################################################################

cdef class BondPattern(Edge):

    cdef public list order

    cpdef copy(self)

    cpdef __changeBond(self, short order)

    cpdef applyAction(self, tuple action)

    cpdef bint equivalent(self, Edge other)

    cpdef bint isSpecificCaseOf(self, Edge other)

################################################################################

cdef class MoleculePattern(Graph):

    cpdef addAtom(self, Atom atom)

    cpdef addBond(self, Atom atom1, Atom atom2, Bond bond)

    cpdef dict getBonds(self, Atom atom)

    cpdef Bond getBond(self, Atom atom1, Atom atom2)

    cpdef bint hasAtom(self, Atom atom)

    cpdef bint hasBond(self, Atom atom1, Atom atom2)

    cpdef removeAtom(self, Atom atom)

    cpdef removeBond(self, Atom atom1, Atom atom2)

    cpdef sortAtoms(self)

    cpdef copy(self, bint deep=?)

    cpdef clearLabeledAtoms(self)

    cpdef bint containsLabeledAtom(self, str label)

    cpdef Atom getLabeledAtom(self, str label)

	cpdef dict getLabeledAtoms(self)

    cpdef fromAdjacencyList(self, str adjlist, bint withLabel=?)

    cpdef toAdjacencyList(self)

    cpdef bint isIsomorphic(self, Graph other, dict initialMap=?)

    cpdef tuple findIsomorphism(self, Graph other, dict initialMap=?)

    cpdef bint isSubgraphIsomorphic(self, Graph other, dict initialMap=?)

    cpdef tuple findSubgraphIsomorphisms(self, Graph other, dict initialMap=?)

################################################################################

cpdef fromAdjacencyList(str adjlist, bint pattern=?, bint addH=?, bint withLabel=?)

cpdef toAdjacencyList(Graph molecule, str label=?, bint pattern=?, bint removeH=?)
