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

cdef class Vertex(object):

    cdef public dict edges

    # These attributes are used in the VF2 graph isomorphism algorithm
    cdef public short connectivity1
    cdef public short connectivity2
    cdef public short connectivity3
    cdef public short sortingLabel
    cdef public bint terminal
    cdef public Vertex mapping
    cdef public bint ignore
    
    cpdef Vertex copy(self)

    cpdef bint equivalent(self, Vertex other) except -2

    cpdef bint isSpecificCaseOf(self, Vertex other) except -2

    cpdef resetConnectivityValues(self)

cpdef short getVertexConnectivityValue(Vertex vertex) except 1 # all values should be negative

cpdef short getVertexSortingLabel(Vertex vertex) except -1 # all values should be nonnegative

################################################################################

cdef class Edge(object):

    cdef public Vertex vertex1, vertex2
    
    cpdef Edge copy(self)

    cpdef bint equivalent(Edge self, Edge other) except -2

    cpdef bint isSpecificCaseOf(self, Edge other) except -2

    cpdef Vertex getOtherVertex(self, Vertex vertex)

################################################################################

cdef Vertex _getEdgeVertex1(Edge edge)

cdef Vertex _getEdgeVertex2(Edge edge)

cdef class Graph:

    cdef public list vertices
    
    cdef public list ordered_vertices

    cpdef Vertex addVertex(self, Vertex vertex)

    cpdef Edge addEdge(self, Edge edge)

    cpdef list getAllEdges(self)

    cpdef dict getEdges(self, Vertex vertex)

    cpdef Edge getEdge(self, Vertex vertex1, Vertex vertex2)

    cpdef bint hasVertex(self, Vertex vertex) except -2

    cpdef bint hasEdge(self, Vertex vertex1, Vertex vertex2) except -2

    cpdef removeVertex(self, Vertex vertex)

    cpdef removeEdge(self, Edge edge)

    cpdef updateConnectivityValues(self)
    
    cpdef Graph copy(self, bint deep=?)

    cpdef dict copyAndMap(self)

    cpdef Graph merge(self, Graph other)

    cpdef list split(self)

    cpdef resetConnectivityValues(self)

    cpdef sortVertices(self, bint saveOrder=?)
    
    cpdef restore_vertex_order(self)

    cpdef bint isIsomorphic(self, Graph other, dict initialMap=?, bint saveOrder=?) except -2

    cpdef list findIsomorphism(self, Graph other, dict initialMap=?, bint saveOrder=?)

    cpdef bint isSubgraphIsomorphic(self, Graph other, dict initialMap=?, bint saveOrder=?) except -2

    cpdef list findSubgraphIsomorphisms(self, Graph other, dict initialMap=?, bint saveOrder=?)

    cpdef bint isCyclic(self) except -2

    cpdef bint isVertexInCycle(self, Vertex vertex) except -2

    cpdef bint isEdgeInCycle(self, Edge edge) except -2

    cpdef bint __isChainInCycle(self, list chain) except -2

    cpdef list getAllCyclicVertices(self)
    
    cpdef list getAllPolycyclicVertices(self)
    
    cpdef list getPolycyclicRings(self)
    
    cpdef list getMonocyclicRings(self)
    
    cpdef tuple getDisparateRings(self)

    cpdef tuple _merge_cycles(self, list cycle_sets)

    cpdef list getAllCycles(self, Vertex startingVertex)

    cpdef list getAllCyclesOfSize(self, int size)

    cpdef list getAllSimpleCyclesOfSize(self, int size)

    cpdef list __exploreCyclesRecursively(self, list chain, list cycles)

    cpdef list getSmallestSetOfSmallestRings(self)

    cpdef list getRelevantCycles(self)

    cpdef list _sortCyclicVertices(self, list vertices)
    
    cpdef list getLargestRing(self, Vertex vertex)
    
    cpdef bint isMappingValid(self, Graph other, dict mapping) except -2
