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

cimport chem

cdef extern from "dictobject.h":
	ctypedef class __builtin__.dict [object PyDictObject]:
		pass

################################################################################

cdef class Vertex:

	cdef public short int connectivity_value_1
	cdef public short int connectivity_value_2
	cdef public short int connectivity_value_3

	cpdef bint equivalent(Vertex self, Vertex other)

################################################################################

cdef class Edge:

	cpdef bint equivalent(Edge self, Edge other)

################################################################################

cdef class Graph(dict):

	cpdef list vertices(Graph self)

	cpdef list edges(Graph self)

	cpdef addVertex(Graph self, vertex)

	cpdef addEdge(Graph self, vertices, edge)

	cpdef dict getEdges(Graph self, vertex)

	cpdef getEdge(Graph self, tuple vertices)

	cpdef bint hasEdge(self, tuple vertices)

	cpdef removeVertex(Graph self, vertex1)

	cpdef removeEdge(Graph self, vertices)

	cpdef isIsomorphic(Graph self, Graph other, dict map12_0, dict map21_0)

	cpdef isSubgraphIsomorphic(Graph self, Graph other, dict map12_0, dict map21_0)

	cpdef findSubgraphIsomorphisms(Graph self, Graph other, dict map12_0, dict map21_0)

	cpdef Graph copy(Graph self)

	cpdef Graph merge(Graph self, Graph other)

	cpdef list split(Graph self)

	cpdef list getSmallestSetOfSmallestRings(Graph self)

	cpdef bint isVertexInCycle(Graph self, chem.Atom vertex)

	cpdef bint __isChainInCycle(Graph self, list chain)

	cpdef list getAllCycles(Graph self, chem.Atom startingVertex)

	cpdef list __exploreCyclesRecursively(Graph self, list chain, list cycleList)

	cpdef set_connectivity_values(Graph self)

	cpdef sort_and_label_vertices(Graph self)

################################################################################

cpdef VF2_isomorphism(Graph graph1, Graph graph2, dict map12, dict map21,
	bint subgraph, bint findAll)

cpdef bint __VF2_feasible(Graph graph1, Graph graph2, chem.Atom vertex1,
	chem.Atom vertex2, dict map21, dict map12, dict terminals1, dict terminals2,
	bint subgraph)

cpdef bint __VF2_match(Graph graph1, Graph graph2, dict map21, dict map12,
	dict terminals1, dict terminals2, bint subgraph, bint findAll,
	list map21List, list map12List, int call_depth)

cpdef int __get_sort_label(chem.Atom vertex)

cpdef int __global_atom_sort_value(chem.Atom atom)

cpdef list __VF2_pairs(Graph graph1, Graph graph2, dict terminals1,
	dict terminals2, dict map21, dict map12)

cpdef dict __VF2_terminals(Graph graph, dict mapping)

cpdef dict __VF2_new_terminals(Graph graph, dict mapping, dict old_terminals,
	new_vertex)
