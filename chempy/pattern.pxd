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

    cpdef bint __atomTypesEquivalent(self, str atomType1, str atomType2)

    cpdef bint __atomTypesSpecificCaseOf(self, str atomType1, str atomType2)

    cpdef bint equivalent(self, other)

    cpdef bint isSpecificCaseOf(self, other)

################################################################################

cdef class BondPattern(Edge):

    cdef public list order

################################################################################

cdef class MoleculePattern(Graph):

    cdef public bint implicitHydrogens
