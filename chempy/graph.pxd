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

cdef class Vertex(object):

    cdef public short connectivity1
    cdef public short connectivity2
    cdef public short connectivity3
    cdef public short sortingLabel

    cpdef bint equivalent(self, Vertex other)

    cpdef bint isSpecificCaseOf(self, Vertex other)

    cpdef resetConnectivityValues(self)

cpdef short getVertexConnectivityValue(Vertex vertex) except 1 # all values should be negative

cpdef short getVertexSortingLabel(Vertex vertex) except -1 # all values should be nonnegative

################################################################################

cdef class Edge(object):

    cpdef bint equivalent(Edge self, Edge other)

    cpdef bint isSpecificCaseOf(self, Edge other)

################################################################################

cdef class Graph:

    cdef public list vertices
    cdef public dict edges

    cpdef Vertex addVertex(self, Vertex vertex)

    cpdef Edge addEdge(self, Vertex vertex1, Vertex vertex2, Edge edge)

    cpdef dict getEdges(self, Vertex vertex)

    cpdef Edge getEdge(self, Vertex vertex1, Vertex vertex2)

    cpdef bint hasVertex(self, Vertex vertex)

    cpdef bint hasEdge(self, Vertex vertex1, Vertex vertex2)

    cpdef removeVertex(self, Vertex vertex)

    cpdef removeEdge(self, Vertex vertex1, Vertex vertex2)

    cpdef Graph copy(self, bint deep=?)

    cpdef Graph merge(self, other)

    cpdef list split(self)

    cpdef resetConnectivityValues(self)

    cpdef updateConnectivityValues(self)

    cpdef sortVertices(self)

    cpdef bint isIsomorphic(self, Graph other, dict initialMap=?)

    cpdef tuple findIsomorphism(self, Graph other, dict initialMap=?)

    cpdef bint isSubgraphIsomorphic(self, Graph other, dict initialMap=?)

    cpdef tuple findSubgraphIsomorphisms(self, Graph other, dict initialMap=?)

    cpdef bint isCyclic(self)

    cpdef bint isVertexInCycle(self, Vertex vertex)

    cpdef bint isEdgeInCycle(self, Vertex vertex1, Vertex vertex2)

    cpdef bint __isChainInCycle(self, list chain)

    cpdef getAllCycles(self, Vertex startingVertex)

    cpdef __exploreCyclesRecursively(self, list chain, list cycleList)

    cpdef getSmallestSetOfSmallestRings(self)

################################################################################

cpdef VF2_isomorphism(Graph graph1, Graph graph2, bint subgraph=?, 
    bint findAll=?, dict initialMap=?)

cpdef bint __VF2_feasible(Graph graph1, Graph graph2, Vertex vertex1,
    Vertex vertex2, dict map21, dict map12, list terminals1, list terminals2,
    bint subgraph) except -2 # bint should be 0 or 1

cpdef bint __VF2_match(Graph graph1, Graph graph2, dict map21, dict map12,
    list terminals1, list terminals2, bint subgraph, bint findAll,
    list map21List, list map12List, int call_depth) except -2 # bint should be 0 or 1

cpdef list __VF2_terminals(Graph graph, dict mapping)

cpdef list __VF2_updateTerminals(Graph graph, dict mapping, list old_terminals,
    Vertex new_vertex)
