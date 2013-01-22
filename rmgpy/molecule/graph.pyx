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

"""
This module contains an implementation of a graph data structure (the 
:class:`Graph` class) and functions for manipulating that graph, including 
efficient isomorphism functions. This module also contains base classes for
the vertices and edges (:class:`Vertex` and :class:`Edge`, respectively) that
are the components of a graph.
"""

import logging

################################################################################

cdef class Vertex(object):
    """
    A base class for vertices in a graph. Contains several connectivity values
    useful for accelerating isomorphism searches, as proposed by
    `Morgan (1965) <http://dx.doi.org/10.1021/c160017a018>`_.

    =================== =============== ========================================
    Attribute           Type            Description
    =================== =============== ========================================
    `connectivity1`     ``int``         The number of nearest neighbors
    `connectivity2`     ``int``         The sum of the neighbors' `connectivity1` values
    `connectivity3`     ``int``         The sum of the neighbors' `connectivity2` values
    `sortingLabel`      ``int``         An integer label used to sort the vertices
    =================== =============== ========================================
    
    """

    def __init__(self):
        self.edges = {}
        self.resetConnectivityValues()

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        d = {
            'edges': self.edges,
            'connectivity1': self.connectivity1,
            'connectivity2': self.connectivity2,
            'connectivity3': self.connectivity3,
            'sortingLabel': self.sortingLabel,
            'terminal': self.terminal,
            'mapping': self.mapping,
        }
        return (Vertex, (), d)

    def __setstate__(self, d):
        self.edges = d['edges']
        self.connectivity1 = d['connectivity1']
        self.connectivity2 = d['connectivity2']
        self.connectivity3 = d['connectivity3']
        self.sortingLabel = d['sortingLabel']
        self.terminal = d['terminal']
        self.mapping = d['mapping']

    cpdef Vertex copy(self):
        """
        Return a copy of the vertex. The default implementation assumes that no
        semantic information is associated with each vertex, and therefore
        simply returns a new :class:`Vertex` object.
        """
        new = Vertex()
        return new

    cpdef bint equivalent(self, Vertex other) except -2:
        """
        Return :data:`True` if two vertices `self` and `other` are semantically
        equivalent, or :data:`False` if not. You should reimplement this
        function in a derived class if your vertices have semantic information.
        """
        return True

    cpdef bint isSpecificCaseOf(self, Vertex other) except -2:
        """
        Return ``True`` if `self` is semantically more specific than `other`,
        or ``False`` if not. You should reimplement this function in a derived
        class if your edges have semantic information.
        """
        return True

    cpdef resetConnectivityValues(self):
        """
        Reset the cached structure information for this vertex.
        """
        self.connectivity1 = -1
        self.connectivity2 = -1
        self.connectivity3 = -1
        self.sortingLabel = -1
        self.terminal = False
        self.mapping = None

cpdef short getVertexConnectivityValue(Vertex vertex) except 1:
    """
    Return a value used to sort vertices prior to poposing candidate pairs in
    :meth:`__VF2_pairs`. The value returned is based on the vertex's
    connectivity values (and assumes that they are set properly).
    """
    return ( -256*vertex.connectivity1 - 16*vertex.connectivity2 - vertex.connectivity3 )

cpdef short getVertexSortingLabel(Vertex vertex) except -1:
    """
    Return a value used to sort vertices prior to poposing candidate pairs in
    :meth:`__VF2_pairs`. The value returned is based on the vertex's
    connectivity values (and assumes that they are set properly).
    """
    return vertex.sortingLabel

################################################################################

cdef class Edge(object):
    """
    A base class for edges in a graph. This class does *not* store the vertex
    pair that comprises the edge; that functionality would need to be included
    in the derived class.
    """

    def __init__(self, vertex1, vertex2):
        self.vertex1 = vertex1
        self.vertex2 = vertex2

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (Edge, (self.vertex1, self.vertex2))

    cpdef Edge copy(self):
        """
        Return a copy of the edge. The default implementation assumes that no
        semantic information is associated with each edge, and therefore
        simply returns a new :class:`Edge` object. Note that the vertices are
        not copied in this implementation.
        """
        new = Edge(self.vertex1, self.vertex2)
        return new

    cpdef bint equivalent(self, Edge other) except -2:
        """
        Return ``True`` if two edges `self` and `other` are semantically
        equivalent, or ``False`` if not. You should reimplement this
        function in a derived class if your edges have semantic information.
        """
        return True

    cpdef bint isSpecificCaseOf(self, Edge other) except -2:
        """
        Return ``True`` if `self` is semantically more specific than `other`,
        or ``False`` if not. You should reimplement this function in a derived
        class if your edges have semantic information.
        """
        return True

    cpdef Vertex getOtherVertex(self, Vertex vertex):
        """
        Given a vertex that makes up part of the edge, return the other vertex.
        Raise a :class:`ValueError` if the given vertex is not part of the
        edge.
        """
        if self.vertex1 is vertex:
            return self.vertex2
        elif self.vertex2 is vertex:
            return self.vertex1
        else:
            raise ValueError('The given vertex is not one of the vertices of this edge.')

################################################################################

cdef class Graph:
    """
    A graph data type. The vertices of the graph are stored in a list
    `vertices`; this provides a consistent traversal order. The edges of the
    graph are stored in a dictionary of dictionaries `edges`. A single edge can
    be accessed using ``graph.edges[vertex1][vertex2]`` or the :meth:`getEdge`
    method; in either case, an exception will be raised if the edge does not
    exist. All edges of a vertex can be accessed using ``graph.edges[vertex]``
    or the :meth:`getEdges` method.
    """

    def __init__(self, vertices=None):
        self.vertices = vertices or []
        
    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (Graph, (self.vertices,))

    cpdef Vertex addVertex(self, Vertex vertex):
        """
        Add a `vertex` to the graph. The vertex is initialized with no edges.
        """
        self.vertices.append(vertex)
        vertex.edges = dict()
        return vertex

    cpdef Edge addEdge(self, Edge edge):
        """
        Add an `edge` to the graph. The two vertices in the edge must already
        exist in the graph, or a :class:`ValueError` is raised.
        """
        if edge.vertex1 not in self.vertices or edge.vertex2 not in self.vertices:
            raise ValueError('Attempted to add edge between vertices not in the graph.')
        edge.vertex1.edges[edge.vertex2] = edge
        edge.vertex2.edges[edge.vertex1] = edge
        return edge

    cpdef dict getEdges(self, Vertex vertex):
        """
        Return a list of the edges involving the specified `vertex`.
        """
        return vertex.edges

    cpdef Edge getEdge(self, Vertex vertex1, Vertex vertex2):
        """
        Returns the edge connecting vertices `vertex1` and `vertex2`.
        """
        try:
            return vertex1.edges[vertex2]
        except KeyError:
            raise ValueError('The specified vertices are not connected by an edge in this graph.')

    cpdef bint hasVertex(self, Vertex vertex) except -2:
        """
        Returns ``True`` if `vertex` is a vertex in the graph, or ``False`` if
        not.
        """
        return vertex in self.vertices

    cpdef bint hasEdge(self, Vertex vertex1, Vertex vertex2) except -2:
        """
        Returns ``True`` if vertices `vertex1` and `vertex2` are connected
        by an edge, or ``False`` if not.
        """
        return vertex1 in self.vertices and vertex2 in vertex1.edges

    cpdef removeVertex(self, Vertex vertex):
        """
        Remove `vertex` and all edges associated with it from the graph. Does
        not remove vertices that no longer have any edges as a result of this
        removal.
        """
        cdef Vertex vertex2
        for vertex2 in vertex.edges:
            del vertex2.edges[vertex]
        vertex.edges = dict()
        self.vertices.remove(vertex)

    cpdef removeEdge(self, Edge edge):
        """
        Remove the specified `edge` from the graph.
        Does not remove vertices that no longer have any edges as a result of
        this removal.
        """
        del edge.vertex1.edges[edge.vertex2]
        del edge.vertex2.edges[edge.vertex1]

    cpdef Graph copy(self, bint deep=False):
        """
        Create a copy of the current graph. If `deep` is ``True``, a deep copy
        is made: copies of the vertices and edges are used in the new graph.
        If `deep` is ``False`` or not specified, a shallow copy is made: the
        original vertices and edges are used in the new graph.
        """
        cdef Graph other
        cdef Vertex vertex, vertex1, vertex2
        cdef Edge edge
        cdef dict edges, mapping
        cdef list vertices
        cdef int index1, index2
        
        other = Graph()
        vertices = self.vertices
        mapping = {}
        for vertex in vertices:
            if deep:
                vertex2 = other.addVertex(vertex.copy())
                mapping[vertex] = vertex2
            else:
                edges = vertex.edges
                other.addVertex(vertex)
                vertex.edges = edges
        if deep:
            for vertex1 in vertices:
                for vertex2 in vertex1.edges:
                    edge = vertex1.edges[vertex2]
                    edge = edge.copy()
                    edge.vertex1 = mapping[vertex1]
                    edge.vertex2 = mapping[vertex2]
                    other.addEdge(edge)
        return other

    cpdef Graph merge(self, Graph other):
        """
        Merge two graphs so as to store them in a single Graph object.
        """
        cdef Graph new
        cdef Vertex vertex, vertex1, vertex2
        
        # Create output graph
        new = Graph()

        # Add vertices to output graph
        for vertex in self.vertices:
            edges = vertex.edges
            new.addVertex(vertex)
            vertex.edges = edges
        for vertex in other.vertices:
            edges = vertex.edges
            new.addVertex(vertex)
            vertex.edges = edges

        return new

    cpdef list split(self):
        """
        Convert a single Graph object containing two or more unconnected graphs
        into separate graphs.
        """
        cdef Graph new1, new2
        cdef Vertex vertex, vertex1, vertex2
        cdef list verticesToMove
        cdef int index
        
        # Create potential output graphs
        new1 = self.copy()
        new2 = Graph()

        if len(self.vertices) == 0:
            return [new1]

        # Arbitrarily choose last atom as starting point
        verticesToMove = [ self.vertices[-1] ]

        # Iterate until there are no more atoms to move
        index = 0
        while index < len(verticesToMove):
            for v2 in verticesToMove[index].edges:
                if v2 not in verticesToMove:
                    verticesToMove.append(v2)
            index += 1
        
        # If all atoms are to be moved, simply return new1
        if len(new1.vertices) == len(verticesToMove):
            return [new1]

        # Copy to new graph and remove from old graph
        for vertex in verticesToMove:
            new2.vertices.append(vertex)
            new1.vertices.remove(vertex)
        
        new = [new2]
        new.extend(new1.split())
        return new

    cpdef resetConnectivityValues(self):
        """
        Reset any cached connectivity information. Call this method when you
        have modified the graph.
        """
        cdef Vertex vertex
        for vertex in self.vertices: vertex.resetConnectivityValues()
        
    cpdef updateConnectivityValues(self):
        """
        Update the connectivity values for each vertex in the graph. These are
        used to accelerate the isomorphism checking.
        """
        cdef Vertex vertex1, vertex2
        cdef short count
        
        for vertex1 in self.vertices:
            count = len(vertex1.edges)
            vertex1.connectivity1 = count
        for vertex1 in self.vertices:
            count = 0
            for vertex2 in vertex1.edges: count += vertex2.connectivity1
            vertex1.connectivity2 = count
        for vertex1 in self.vertices:
            count = 0
            for vertex2 in vertex1.edges: count += vertex2.connectivity2
            vertex1.connectivity3 = count
        
    cpdef sortVertices(self):
        """
        Sort the vertices in the graph. This can make certain operations, e.g.
        the isomorphism functions, much more efficient.
        """
        cdef Vertex vertex
        cdef int index
        # Only need to conduct sort if there is an invalid sorting label on any vertex
        for vertex in self.vertices:
            if vertex.sortingLabel < 0: break
        else:
            return
        # If we need to sort then let's also update the connecitivities so
        # we're sure they are right, since the sorting labels depend on them
        self.updateConnectivityValues()
        self.vertices.sort(key=getVertexConnectivityValue)
        for index, vertex in enumerate(self.vertices):
            vertex.sortingLabel = index

    cpdef bint isIsomorphic(self, Graph other, dict initialMap=None) except -2:
        """
        Returns :data:`True` if two graphs are isomorphic and :data:`False`
        otherwise. Uses the VF2 algorithm of Vento and Foggia.
        """
        return VF2_isomorphism(self, other, False, False, initialMap)

    cpdef list findIsomorphism(self, Graph other, dict initialMap=None):
        """
        Returns :data:`True` if `other` is subgraph isomorphic and :data:`False`
        otherwise, and the matching mapping.
        Uses the VF2 algorithm of Vento and Foggia.
        """
        return VF2_isomorphism(self, other, False, True, initialMap)

    cpdef bint isSubgraphIsomorphic(self, Graph other, dict initialMap=None) except -2:
        """
        Returns :data:`True` if `other` is subgraph isomorphic and :data:`False`
        otherwise. Uses the VF2 algorithm of Vento and Foggia.
        """
        return VF2_isomorphism(self, other, True, False, initialMap)

    cpdef list findSubgraphIsomorphisms(self, Graph other, dict initialMap=None):
        """
        Returns :data:`True` if `other` is subgraph isomorphic and :data:`False`
        otherwise. Also returns the lists all of valid mappings.

        Uses the VF2 algorithm of Vento and Foggia.
        """
        return VF2_isomorphism(self, other, True, True, initialMap)

    cpdef bint isCyclic(self) except -2:
        """
        Return ``True`` if one or more cycles are present in the graph or
        ``False`` otherwise.
        """
        cdef Vertex vertex
        for vertex in self.vertices:
            if self.isVertexInCycle(vertex):
                return True
        return False

    cpdef bint isVertexInCycle(self, Vertex vertex) except -2:
        """
        Return ``True`` if the given `vertex` is contained in one or more
        cycles in the graph, or ``False`` if not.
        """
        return self.__isChainInCycle([vertex])

    cpdef bint isEdgeInCycle(self, Edge edge) except -2:
        """
        Return :data:`True` if the edge between vertices `vertex1` and `vertex2`
        is in one or more cycles in the graph, or :data:`False` if not.
        """
        cdef list cycles
        cycles = self.getAllCycles(edge.vertex1)
        for cycle in cycles:
            if edge.vertex2 in cycle:
                return True
        return False

    cpdef bint __isChainInCycle(self, list chain) except -2:
        """
        Return ``True`` if the given `chain` of vertices is contained in one
        or more cycles or ``False`` otherwise. This function recursively calls
        itself.
        """
        cdef Vertex vertex1, vertex2
        cdef Edge edge

        vertex1 = chain[-1]
        for vertex2 in vertex1.edges:
            if vertex2 is chain[0] and len(chain) > 2:
                return True
            elif vertex2 not in chain:
                # Make the chain a little longer and explore again
                chain.append(vertex2)
                if self.__isChainInCycle(chain):
                    # We found a cycle, so the return value must be True
                    return True
                else:
                    # We did not find a cycle down this path, so remove the vertex from the chain
                    chain.remove(vertex2)
        # If we reach this point then we did not find any cycles involving this chain
        return False
    
    cpdef list getAllCyclicVertices(self):
        """ 
        Returns all vertices belonging to one or more cycles.        
        """
        cdef list cyclicVertices
        # Loop through all vertices and check whether they are cyclic
        cyclicVertices = []
        for vertex in self.vertices:
            if self.isVertexInCycle(vertex):
                cyclicVertices.append(vertex)                
        return cyclicVertices
    
    cpdef list getAllPolycyclicVertices(self):
        """
        Return all vertices belonging to two or more cycles, fused or spirocyclic.
        """
        cdef list SSSR, vertices, polycyclicVertices
        SSSR = self.getSmallestSetOfSmallestRings()
        polycyclicVertices = []
        if SSSR:            
            vertices = []
            for cycle in SSSR:
                for vertex in cycle:
                    if vertex not in vertices:
                        vertices.append(vertex)
                    else:
                        if vertex not in polycyclicVertices:
                            polycyclicVertices.append(vertex)     
        return polycyclicVertices                                    

    cpdef list getAllCycles(self, Vertex startingVertex):
        """
        Given a starting vertex, returns a list of all the cycles containing
        that vertex.
        """
        return self.__exploreCyclesRecursively([startingVertex], [])

    cpdef list __exploreCyclesRecursively(self, list chain, list cycles):
        """
        Search the graph for cycles by recursive spidering. Given a `chain`
        (list) of connected atoms and a list of `cycles` found so far, find any
        cycles involving the chain of atoms and append them to the list of
        cycles. This function recursively calls itself.
        """
        cdef Vertex vertex1, vertex2
        
        vertex1 = chain[-1]
        # Loop over each of the atoms neighboring the last atom in the chain
        for vertex2 in vertex1.edges:
            if vertex2 is chain[0] and len(chain) > 2:
                # It is the first atom in the chain, so the chain is a cycle!
                cycles.append(chain[:])
            elif vertex2 not in chain:
                # Make the chain a little longer and explore again
                chain.append(vertex2)
                cycles = self.__exploreCyclesRecursively(chain, cycles)
                # Any cycles down this path have now been found, so remove vertex2 from the chain
                chain.pop(-1)
        # At this point we should have discovered all of the cycles involving the current chain
        return cycles

    cpdef list getSmallestSetOfSmallestRings(self):
        """
        Return a list of the smallest set of smallest rings in the graph. The
        algorithm implements was adapted from a description by Fan, Panaye,
        Doucet, and Barbu (doi: 10.1021/ci00015a002)

        B. T. Fan, A. Panaye, J. P. Doucet, and A. Barbu. "Ring Perception: A
        New Algorithm for Directly Finding the Smallest Set of Smallest Rings
        from a Connection Table." *J. Chem. Inf. Comput. Sci.* **33**,
        p. 657-662 (1993).
        """
        cdef Graph graph
        cdef bint done, found
        cdef list cycleList, cycles, cycle, graphs, neighbors, verticesToRemove, vertices
        cdef Vertex vertex, rootVertex

        # Make a copy of the graph so we don't modify the original
        graph = self.copy(deep=True)
        vertices = graph.vertices[:]
        
        # Step 1: Remove all terminal vertices
        done = False
        while not done:
            verticesToRemove = []
            for vertex in graph.vertices:
                if len(vertex.edges) == 1: verticesToRemove.append(vertex)
            done = len(verticesToRemove) == 0
            # Remove identified vertices from graph
            for vertex in verticesToRemove:
                graph.removeVertex(vertex)

        # Step 2: Remove all other vertices that are not part of cycles
        verticesToRemove = []
        for vertex in graph.vertices:
            found = graph.isVertexInCycle(vertex)
            if not found:
                verticesToRemove.append(vertex)
        # Remove identified vertices from graph
        for vertex in verticesToRemove:
            graph.removeVertex(vertex)

        # Step 3: Split graph into remaining subgraphs
        graphs = graph.split()

        # Step 4: Find ring sets in each subgraph
        cycleList = []
        for graph in graphs:

            while len(graph.vertices) > 0:

                # Choose root vertex as vertex with smallest number of edges
                rootVertex = None
                for vertex in graph.vertices:
                    if rootVertex is None:
                        rootVertex = vertex
                    elif len(vertex.edges) < len(rootVertex.edges):
                        rootVertex = vertex

                # Get all cycles involving the root vertex
                cycles = graph.getAllCycles(rootVertex)
                if len(cycles) == 0:
                    # This vertex is no longer in a ring, so remove it
                    graph.removeVertex(rootVertex)
                    continue

                # Keep the smallest of the cycles found above
                cycle = cycles[0]
                for c in cycles[1:]:
                    if len(c) < len(cycle):
                        cycle = c
                cycleList.append(cycle)

                # Remove from the graph all vertices in the cycle that have only two edges
                verticesToRemove = []
                for vertex in cycle:
                    if len(vertex.edges) <= 2:
                        verticesToRemove.append(vertex)
                if len(verticesToRemove) == 0:
                    # there are no vertices in this cycle that with only two edges
                    # Remove edge between root vertex and any one vertex it is connected to
                    vertex = rootVertex.edges.keys()[0]
                    graph.removeEdge(rootVertex.edges[vertex])
                else:
                    for vertex in verticesToRemove:
                        graph.removeVertex(vertex)

        # Map atoms in cycles back to atoms in original graph
        for i in range(len(cycleList)):
            cycleList[i] = [self.vertices[vertices.index(v)] for v in cycleList[i]]

        return cycleList

    cpdef bint isMappingValid(self, Graph other, dict mapping) except -2:
        """
        Check that a proposed `mapping` of vertices from `self` to `other`
        is valid by checking that the vertices and edges involved in the
        mapping are mutually equivalent.
        """
        cdef Vertex vertex1, vertex2
        cdef list vertices1, vert
        cdef bint selfHasEdge, otherHasEdge
        cdef int i, j
        
        # Check that the mapped pairs of vertices are equivalent
        for vertex1, vertex2 in mapping.items():
            if not vertex1.equivalent(vertex2):
                return False
        
        # Check that any edges connected mapped vertices are equivalent
        vertices1 = mapping.keys()
        vertices2 = mapping.values()
        for i in range(len(vertices1)):
            for j in range(i+1, len(vertices1)):
                selfHasEdge = self.hasEdge(vertices1[i], vertices1[j])
                otherHasEdge = other.hasEdge(vertices2[i], vertices2[j])
                if selfHasEdge and otherHasEdge:
                    # Both graphs have the edge, so we must check it for equivalence
                    edge1 = self.getEdge(vertices1[i], vertices1[j])
                    edge2 = other.getEdge(vertices2[i], vertices2[j])
                    if not edge1.equivalent(edge2):
                        return False
                elif selfHasEdge or otherHasEdge:
                    # Only one of the graphs has the edge, so the mapping must be invalid
                    return False
        
        # If we're here then the vertices and edges are equivalent, so the
        # mapping is valid
        return True

################################################################################

class VF2Error(Exception):
    """
    An exception raised if an error occurs within the VF2 graph isomorphism
    algorithm. Pass a string describing the error.
    """
    pass

cpdef VF2_isomorphism(Graph graph1, Graph graph2, bint subgraph=False, bint findAll=False, dict initialMapping=None):
    """
    Use the VF2 algorithm of Vento and Foggia to evaluate the isomorphism of
    the graphs `graph1` and `graph2`. A number of options affect how the
    isomorphism check is performed:
    
    * If `subgraph` is ``True``, the function will determine if `graph2` is a
      subgraph of `graph1` instead of a full graph.
    
    * If `findAll` is ``True``, this function returns a list of valid mappings
      from `graph1` to `graph2`; each mapping is a ``dict`` with vertices from
      `graph1` as the keys and vertices from `graph2` as the values. If
      `findAll` is ``False``, this function simply returns ``True`` if a
      valid mapping was found, or ``False`` if not.
    
    * The `initialMapping` parameter is used to specify a mapping of vertices
      from `graph1` to those in `subgraph` that are fixed in the isomorphism
      check; this mapping will appear in every returned mapping. Note that no
      validation of this initial mapping is performed in this function.
    """
    cdef Vertex vertex1, vertex2
    cdef list mappingList
    cdef int callDepth
    
    # Some quick isomorphism checks based on graph sizes
    if not subgraph and len(graph2.vertices) != len(graph1.vertices):
        # The two graphs don't have the same number of vertices, so they
        # cannot be isomorphic
        return list() if findAll else False
    elif not subgraph and len(graph2.vertices) == len(graph1.vertices) == 0:
        # The two graphs don't have any vertices; this means they are
        # trivially isomorphic
        return list() if findAll else True
    elif subgraph and len(graph2.vertices) > len(graph1.vertices):
        # The second graph has more vertices than the first, so it cannot be
        # a subgraph of the first
        return list() if findAll else False

    # Initialize callDepth with the size of the smallest graph
    # Each recursive call to VF2_match will decrease it by one;
    # when the whole graph has been explored, it should reach 0
    # It should never go below zero!
    callDepth = len(graph2.vertices)

    # Sort the vertices in each graph to make the isomorphism more efficient
    graph1.sortVertices()
    graph2.sortVertices()

    # Initialize mapping by clearing any previous mapping information
    for vertex1 in graph1.vertices:
        vertex1.mapping = None
        vertex1.terminal = False
    for vertex2 in graph2.vertices:
        vertex2.mapping = None
        vertex2.terminal = False
    # Set the initial mapping if provided
    if initialMapping is not None:
        for vertex1, vertex2 in initialMapping.items():
            VF2_addToMapping(vertex1, vertex2)
        callDepth -= len(initialMapping)

    mappingList = []
    isMatch = VF2_match(graph1, graph2, subgraph, findAll, mappingList, callDepth)

    return mappingList if findAll else isMatch

cdef bint VF2_feasible(Graph graph1, Graph graph2, Vertex vertex1, Vertex vertex2, bint subgraph) except -2:
    """
    Return ``True`` if two vertices `vertex1` and `vertex2` from graphs
    `graph1` and `graph2`, respectively, are feasible matches. The `subgraph` 
    flag indicates whether or not to treat `graph2` as a subgraph of `graph1`.

    The feasibility is assessed through a series of semantic and structural
    checks. Only the combination of the semantic checks and the level 0
    structural check are both necessary and sufficient to ensure feasibility.
    (This does *not* mean that `vertex1` and `vertex2` are always a match,
    although the level 1 and level 2 checks preemptively eliminate a number of
    false positives.)
    """
    cdef Vertex vert1, vert2
    cdef Edge edge1, edge2
    cdef int term1Count, term2Count, neither1Count, neither2Count

    if not subgraph:
        # To be feasible the connectivity values must be an exact match
        if vertex1.connectivity1 != vertex2.connectivity1: return False
        if vertex1.connectivity2 != vertex2.connectivity2: return False
        if vertex1.connectivity3 != vertex2.connectivity3: return False
    
    # Semantic check #1: vertex1 and vertex2 must be equivalent
    if subgraph:
        if not vertex1.isSpecificCaseOf(vertex2): return False
    else:
        if not vertex1.equivalent(vertex2): return False
    
    # Semantic check #2: adjacent vertices to vertex1 and vertex2 that are
    # already mapped should be connected by equivalent edges
    for vert2 in vertex2.edges:
        if vert2.mapping is not None:
            vert1 = vert2.mapping
            if vert1 not in vertex1.edges:
                # The vertices are joined in graph2, but not in graph1
                return False
            edge1 = vertex1.edges[vert1]
            edge2 = vertex2.edges[vert2]
            if subgraph:
                if not edge1.isSpecificCaseOf(edge2): return False
            else:
                if not edge1.equivalent(edge2): return False

    # There could still be edges in graph1 that aren't in graph2; this is okay
    # for subgraph matching, but not for exact matching
    if not subgraph:
        for vert1 in vertex1.edges:
            if vert1.mapping is not None:
                if vert1.mapping not in vertex2.edges: 
                    # The vertices are joined in graph1, but not in graph2
                    return False

    # Count number of terminals adjacent to vertex1 and vertex2
    term1Count = 0; term2Count = 0; neither1Count = 0; neither2Count = 0
    for vert1 in vertex1.edges:
        if vert1.terminal: term1Count += 1
        elif vert1.mapping is not None: neither1Count += 1
    for vert2 in vertex2.edges:
        if vert2.terminal: term2Count += 1
        elif vert2.mapping is not None: neither2Count += 1

    # Level 2 look-ahead: the number of adjacent vertices of vertex1 and
    # vertex2 that are non-terminals must be equal
    if subgraph:
        if neither1Count < neither2Count: return False
    else:
        if neither1Count != neither2Count: return False

    # Level 1 look-ahead: the number of adjacent vertices of vertex1 and
    # vertex2 that are terminals must be equal
    if subgraph:
        if term1Count < term2Count: return False
    else:
        if term1Count != term2Count: return False

    # Level 0 look-ahead: all adjacent vertices of vertex2 already in the
    # mapping must map to adjacent vertices of vertex1
    for vert2 in vertex2.edges:
        if vert2.mapping is not None:
            if vert2.mapping not in vertex1.edges: return False
    # Also, all adjacent vertices of vertex1 already in the mapping must map to
    # adjacent vertices of vertex2, unless we are subgraph matching
    if not subgraph:
        for vert1 in vertex1.edges:
            if vert1.mapping is not None:
                if vert1.mapping not in vertex2.edges: return False

    # All of our tests have been passed, so the two vertices are a feasible pair
    return True

cdef bint VF2_match(Graph graph1, Graph graph2, bint subgraph, bint findAll, list mappingList, int callDepth) except -2:
    """
    Recursively explore two graphs `graph1` and `graph2` in search of one or
    more isomorphism relationships by attempting to map vertices to one
    another. The `subgraph` flag indicates whether or not to treat `graph2` as
    a subgraph of `graph1`. The `findAll` flag indicates whether to find all
    isomorphisms or only the first. The `mappingList` parameter stores the
    current list of found mappings. The `callDepth` parameter keeps track of
    how many matching pairs of vertices have been identified, and is used to
    know when an isomorphism is found. Returns ``True`` if at least one
    isomorphism was found or ``False`` if none were found.
    """
    cdef Vertex vertex1, vertex2
    cdef dict mapping
    
    # The call depth should never be negative!
    if callDepth < 0:
        raise VF2Error('Negative call depth encountered in VF2_match().')

    # Done if we have mapped to all vertices in graph
    if callDepth == 0:
        if findAll:
            mapping = {}
            for vertex2 in graph2.vertices:
                assert vertex2.mapping is not None
                assert vertex2.mapping.mapping is vertex2
                mapping[vertex2.mapping] = vertex2
            mappingList.append(mapping)
        return True

    # Create list of pairs of candidates for inclusion in mapping
    vertices1 = []
    for vertex2 in graph2.vertices:
        if vertex2.terminal:
            # graph2 has terminals, so graph1 also must have terminals
            for vertex1 in graph1.vertices:
                if vertex1.terminal:
                    vertices1.append(vertex1)
            break
    else:
        # graph2 does not have terminals, so graph1 also must not have terminals
        vertex2 = graph2.vertices[0]
        vertices1 = graph1.vertices[:]
    
    for vertex1 in vertices1:
        # Propose a pairing
        if VF2_feasible(graph1, graph2, vertex1, vertex2, subgraph):
            # Add proposed match to mapping
            VF2_addToMapping(vertex1, vertex2)
            # Recurse
            isMatch = VF2_match(graph1, graph2, subgraph, findAll, mappingList, callDepth-1)
            if isMatch:
                if not findAll:
                    return True
            # Undo proposed match
            VF2_removeFromMapping(vertex1, vertex2)

    # None of the proposed matches led to a complete isomorphism, so return False
    return False

cdef void VF2_addToMapping(Vertex vertex1, Vertex vertex2):
    """
    Add a pair of vertices `vertex1` and `vertex2` to the current mapping,
    and update the terminals information for the neighbors of each vertex.
    """
    cdef Vertex v
    
    # Map the vertices to one another
    vertex1.mapping = vertex2
    vertex2.mapping = vertex1
    
    # Remove these vertices from the set of terminals
    vertex1.terminal = False
    vertex2.terminal = False
    
    # Add any neighboring vertices not already in mapping to terminals
    for v in vertex1.edges:
        v.terminal = v.mapping is None
    for v in vertex2.edges:
        v.terminal = v.mapping is None

cdef void VF2_removeFromMapping(Vertex vertex1, Vertex vertex2):
    """
    Remove a pair of vertices `vertex1` and `vertex2` from the current mapping,
    and update the terminals information for the neighbors of each vertex.
    """
    cdef Vertex v, v2
    
    # Unmap the vertices from one another
    vertex1.mapping = None
    vertex2.mapping = None
    
    # Restore these vertices to the set of terminals
    for v in vertex1.edges:
        if v.mapping is not None:
            vertex1.terminal = True
            break
        else:
            vertex1.terminal = False
    for v in vertex2.edges:
        if v.mapping is not None:
            vertex2.terminal = True
            break
        else:
            vertex2.terminal = False
    
    # Recompute the terminal status of any neighboring atoms
    for v in vertex1.edges:
        if v.mapping is not None: continue
        for v2 in v.edges:
            if v2.mapping is not None:
                v.terminal = True
                break
        else:
            v.terminal = False
    for v in vertex2.edges:
        if v.mapping is not None: continue
        for v2 in v.edges:
            if v2.mapping is not None:
                v.terminal = True
                break
        else:
            v.terminal = False
