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

cdef class Vertex(object):

    cdef public dict edges

    # These attributes are used in the VF2 graph isomorphism algorithm
    cdef public long connectivity
    cdef public short sortingLabel
    cdef public bint terminal
    cdef public Vertex mapping
    cdef public bint ignore

    cpdef Vertex copy(self)

    cpdef bint equivalent(self, Vertex other) except -2

    cpdef bint isSpecificCaseOf(self, Vertex other) except -2

    cpdef resetConnectivityValues(self)

cpdef long getVertexConnectivityValue(Vertex vertex) except 1 # all values should be negative

cpdef short getVertexSortingLabel(Vertex vertex) except -1 # all values should be nonnegative

################################################################################

cdef class Edge(object):

    cdef public Vertex vertex1, vertex2
    
    cpdef Edge copy(self)

    cpdef bint equivalent(Edge self, Edge other) except -2

    cpdef bint isSpecificCaseOf(self, Edge other) except -2

    cpdef Vertex getOtherVertex(self, Vertex vertex)

################################################################################

cdef class Graph:

    cdef public list vertices

    cpdef Vertex addVertex(self, Vertex vertex)

    cpdef Edge addEdge(self, Edge edge)

    cpdef dict getEdges(self, Vertex vertex)

    cpdef Edge getEdge(self, Vertex vertex1, Vertex vertex2)

    cpdef bint hasVertex(self, Vertex vertex) except -2

    cpdef bint hasEdge(self, Vertex vertex1, Vertex vertex2) except -2

    cpdef removeVertex(self, Vertex vertex)

    cpdef removeEdge(self, Edge edge)

    cpdef updateConnectivityValues(self)
    
    cpdef Graph copy(self, bint deep=?)

    cpdef Graph merge(self, Graph other)

    cpdef list split(self)

    cpdef resetConnectivityValues(self)

    cpdef list __update(self, old_values, trial_number)

    cpdef sortVertices(self)

    cpdef bint isIsomorphic(self, Graph other, dict initialMap=?) except -2

    cpdef list findIsomorphism(self, Graph other, dict initialMap=?)

    cpdef bint isSubgraphIsomorphic(self, Graph other, dict initialMap=?) except -2

    cpdef list findSubgraphIsomorphisms(self, Graph other, dict initialMap=?)

    cpdef bint isCyclic(self) except -2

    cpdef bint isVertexInCycle(self, Vertex vertex) except -2

    cpdef bint isEdgeInCycle(self, Edge edge) except -2

    cpdef bint __isChainInCycle(self, list chain) except -2

    cpdef list getAllCyclicVertices(self)
    
    cpdef list getAllPolycyclicVertices(self)

    cpdef list getAllCycles(self, Vertex startingVertex)

    cpdef list __exploreCyclesRecursively(self, list chain, list cycles)

    cpdef list getSmallestSetOfSmallestRings(self)
    
    cpdef bint isMappingValid(self, Graph other, dict mapping) except -2
