################################################################################
#
#	RMG - Reaction Mechanism Generator
#
#	Copyright (c) 2002-2009 Prof. William H. Green (whgreen@mit.edu) and the
#	RMG Team (rmg_dev@mit.edu)
#
#	Permission is hereby granted, free of charge, to any person obtaining a
#	copy of this software and associated documentation files (the 'Software'),
#	to deal in the Software without restriction, including without limitation
#	the rights to use, copy, modify, merge, publish, distribute, sublicense,
#	and/or sell copies of the Software, and to permit persons to whom the
#	Software is furnished to do so, subject to the following conditions:
#
#	The above copyright notice and this permission notice shall be included in
#	all copies or substantial portions of the Software.
#
#	THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#	FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#	DEALINGS IN THE SOFTWARE.
#
################################################################################

cimport graph

cdef class Element:

	cdef public int number
	cdef public str name
	cdef public str symbol
	cdef public float mass
	cdef public list valence

################################################################################

cdef class AtomType:

	cdef public str label
	cdef public Element element
	cdef public str description
	cdef public object doubleBonds
	cdef public object tripleBonds
	cdef public object benzeneBonds
	cdef public object formBond
	cdef public object breakBond
	cdef public object incrementBond
	cdef public object decrementBond

	cpdef bint equivalent(AtomType self, AtomType other)
	cpdef bint isSpecificCaseOf(AtomType self, AtomType other)

################################################################################

cdef class ElectronState:

	cdef public str label
	cdef public int order
	cdef public list spin
	cdef public ElectronState increment
	cdef public ElectronState decrement

	cpdef bint equivalent(ElectronState self, ElectronState other)
	cpdef bint isSpecificCaseOf(ElectronState self, ElectronState other)

################################################################################

cdef class BondType:

	cdef public str label
	cdef public str name
	cdef public float order
	cdef public int piElectrons
	cdef public str location

	cpdef bint equivalent(BondType self, BondType other)
	cpdef bint isSpecificCaseOf(BondType self, BondType other)

################################################################################

cdef class Atom(graph.Vertex):
	
	cdef public list _atomType
	cdef public list _electronState
	cdef public int charge
	cdef public str label

	cpdef bint equivalent(Atom self, graph.Vertex other)
	cpdef bint isSpecificCaseOf(Atom self, graph.Vertex other)

################################################################################

cdef class Bond(graph.Edge):
	
	cdef public list atoms
	cdef public list _bondType

	cpdef bint equivalent(Bond self, graph.Edge other)
	cpdef bint isSpecificCaseOf(Bond self, graph.Edge other)

################################################################################

