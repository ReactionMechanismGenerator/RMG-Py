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
from .group cimport GroupAtom, GroupBond, Group
from .element cimport Element
cimport rmgpy.constants as constants

################################################################################

cdef class Atom(Vertex):

    cdef public Element element
    cdef public short radicalElectrons
    cdef public short spinMultiplicity
    cdef public short charge
    cdef public str label
    cdef public AtomType atomType
    cdef public list coords

    cpdef bint equivalent(self, Vertex other) except -2

    cpdef bint isSpecificCaseOf(self, Vertex other) except -2

    cpdef Vertex copy(self)

    cpdef bint isHydrogen(self)

    cpdef bint isNonHydrogen(self)

    cpdef bint isCarbon(self)

    cpdef bint isOxygen(self)
    
    cpdef incrementRadical(self)

    cpdef decrementRadical(self)
    
################################################################################

cpdef object SMILEwriter
    
cdef class Bond(Edge):

    cdef public str order

    cpdef bint equivalent(self, Edge other) except -2

    cpdef bint isSpecificCaseOf(self, Edge other) except -2

    cpdef Edge copy(self)

    cpdef bint isSingle(self) except -2

    cpdef bint isDouble(self) except -2

    cpdef bint isTriple(self) except -2

    cpdef incrementOrder(self)

    cpdef decrementOrder(self)

################################################################################

cdef class Molecule(Graph):

    cdef public bint implicitHydrogens
    cdef public int symmetryNumber
    cdef public object rdMol
    cdef public int rdMolConfId
    cdef str _fingerprint
    
    cpdef str getFingerprint(self)
    
    cpdef addAtom(self, Atom atom)

    cpdef addBond(self, Bond bond)

    cpdef dict getBonds(self, Atom atom)

    cpdef Bond getBond(self, Atom atom1, Atom atom2)

    cpdef bint hasAtom(self, Atom atom)

    cpdef bint hasBond(self, Atom atom1, Atom atom2)

    cpdef removeAtom(self, Atom atom)

    cpdef removeBond(self, Bond bond)

    cpdef sortAtoms(self)

    cpdef str getFormula(self)

    cpdef short getRadicalCount(self)

    cpdef double getMolecularWeight(self)

    cpdef int getNumAtoms(self, str element=?)

    cpdef int getNumberOfRadicalElectrons(self)

    cpdef Graph copy(self, bint deep=?)

    cpdef deleteHydrogens(self)

    cpdef clearLabeledAtoms(self)

    cpdef bint containsLabeledAtom(self, str label) except -2

    cpdef Atom getLabeledAtom(self, str label)

    cpdef dict getLabeledAtoms(self)

    cpdef bint isIsomorphic(self, Graph other, dict initialMap=?) except -2

    cpdef list findIsomorphism(self, Graph other, dict initialMap=?)

    cpdef bint isSubgraphIsomorphic(self, Graph other, dict initialMap=?) except -2

    cpdef list findSubgraphIsomorphisms(self, Graph other, dict initialMap=?)

    cpdef bint isAtomInCycle(self, Atom atom) except -2

    cpdef bint isBondInCycle(self, Bond bond) except -2

    cpdef draw(self, str path)

    cpdef fromCML(self, str cmlstr)

    cpdef fromInChI(self, str inchistr)

    cpdef fromSMILES(self, str smilesstr)

    cpdef fromOBMol(self, obmol)

    cpdef fromAdjacencyList(self, str adjlist)

    cpdef str toCML(self)

    cpdef str toInChI(self)

    cpdef str toAugmentedInChI(self)

    cpdef str toInChIKey(self)

    cpdef str toAugmentedInChIKey(self)

    cpdef str toSMILES(self)

    cpdef toOBMol(self)

    cpdef toAdjacencyList(self, str label=?, bint removeH=?)

    cpdef bint isLinear(self) except -2

    cpdef int countInternalRotors(self) except -2

    cpdef double calculateCp0(self) except -1

    cpdef double calculateCpInf(self) except -1
    
    cpdef updateAtomTypes(self)
    
    cpdef bint isRadical(self) except -2
    
    cpdef list generateResonanceIsomers(self)
    
    cpdef list getAdjacentResonanceIsomers(self)

    cpdef findAllDelocalizationPaths(self, Atom atom1)

    cpdef int calculateSymmetryNumber(self) except -1
