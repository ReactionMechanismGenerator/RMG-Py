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
This module contains graph ismorphism functions that implement the VF2
algorithm of Vento and Foggia.
"""

cimport cython

################################################################################

class VF2Error(Exception):
    """
    An exception raised if an error occurs within the VF2 graph isomorphism
    algorithm. Pass a string describing the error.
    """
    pass

cdef class VF2:
    """
    An implementation of the second version of the Vento-Foggia (VF2) algorithm
    for graph and subgraph isomorphism.
    """
    
    def __init__(self):
        self.graph1 = None
        self.graph2 = None

    cpdef bint isIsomorphic(self, Graph graph1, Graph graph2, dict initialMapping) except -2:
        """
        Return ``True`` if graph `graph1` is isomorphic to graph `graph2` with
        the optional initial mapping `initialMapping`, or ``False`` otherwise.
        """
        self.isomorphism(graph1, graph2, initialMapping, False, False)
        return self.isMatch
        
    cpdef list findIsomorphism(self, Graph graph1, Graph graph2, dict initialMapping):
        """
        Return a list of dicts of all valid isomorphism mappings from graph
        `graph1` to graph `graph2` with the optional initial mapping 
        `initialMapping`. If no valid isomorphisms are found, an empty list is
        returned.
        """
        self.isomorphism(graph1, graph2, initialMapping, False, True)
        return self.mappingList

    cpdef bint isSubgraphIsomorphic(self, Graph graph1, Graph graph2, dict initialMapping) except -2:
        """
        Return ``True`` if graph `graph1` is subgraph isomorphic to subgraph
        `graph2` with the optional initial mapping `initialMapping`, or
        ``False`` otherwise.
        """
        self.isomorphism(graph1, graph2, initialMapping, True, False)
        return self.isMatch

    cpdef list findSubgraphIsomorphisms(self, Graph graph1, Graph graph2, dict initialMapping):
        """
        Return a list of dicts of all valid subgraph isomorphism mappings from
        graph `graph1` to subgraph `graph2` with the optional initial mapping 
        `initialMapping`. If no valid subgraph isomorphisms are found, an empty
        list is returned.
        """
        self.isomorphism(graph1, graph2, initialMapping, True, True)
        return self.mappingList
        
    cdef isomorphism(self, Graph graph1, Graph graph2, dict initialMapping, bint subgraph, bint findAll):
        """
        Evaluate the isomorphism relationship between graphs `graph1` and
        `graph2` with optional initial mapping `initialMapping`. If `subgraph`
        is ``True``, `graph2` is treated as a possible subgraph of `graph1`.
        If `findAll` is ``True``, all isomorphisms are found; otherwise only
        the first is found.
        """
        cdef int callDepth, index1, index2
        
        if self.graph1 is not graph1:
            self.graph1 = graph1
            graph1.sortVertices()
            
        if self.graph2 is not graph2:
            self.graph2 = graph2
            graph2.sortVertices()
        
        self.initialMapping = initialMapping
        self.subgraph = subgraph
        self.findAll = findAll
    
        # Clear previous result
        self.isMatch = False
        self.mappingList = []
        
        # Some quick isomorphism checks based on graph sizes
        if not self.subgraph and len(graph2.vertices) != len(graph1.vertices):
            # The two graphs don't have the same number of vertices, so they
            # cannot be isomorphic
            return
        elif not self.subgraph and len(graph2.vertices) == len(graph1.vertices) == 0:
            # The two graphs don't have any vertices; this means they are
            # trivially isomorphic
            self.isMatch = True
            return
        elif self.subgraph and len(graph2.vertices) > len(graph1.vertices):
            # The second graph has more vertices than the first, so it cannot be
            # a subgraph of the first
            return

        # Initialize callDepth with the size of the smallest graph
        # Each recursive call to match() will decrease it by one;
        # when the whole graph has been explored, it should reach 0
        # It should never go below zero!
        callDepth = len(graph2.vertices)

        # Initialize mapping by clearing any previous mapping information
        for vertex1 in graph1.vertices:
            vertex1.mapping = None
            vertex1.terminal = False
        for vertex2 in graph2.vertices:
            vertex2.mapping = None
            vertex2.terminal = False
        # Set the initial mapping if provided
        if self.initialMapping is not None:
            for vertex1, vertex2 in self.initialMapping.items():
                self.addToMapping(vertex1, vertex2)
            callDepth -= len(self.initialMapping)
            
        self.match(callDepth)

    cdef bint match(self, int callDepth) except -2:
        """
        Recursively search for pairs of vertices to match, until all vertices
        are matched or the viable set of matches is exhausted. The `callDepth`
        parameter helps ensure we never enter an infinite loop.
        """
        cdef Vertex vertex1, vertex2
        cdef dict mapping
        cdef bint hasTerminals
        
        # The call depth should never be negative!
        if callDepth < 0:
            raise VF2Error('Negative call depth encountered in VF2_match().')

        # Done if we have mapped to all vertices in graph
        if callDepth == 0:
            if self.findAll:
                mapping = {}
                for vertex2 in self.graph2.vertices:
                    assert vertex2.mapping is not None
                    assert vertex2.mapping.mapping is vertex2
                    mapping[vertex2.mapping] = vertex2
                self.mappingList.append(mapping)
            self.isMatch = True
            return True

        # Create list of pairs of candidates for inclusion in mapping
        hasTerminals = False
        for vertex2 in self.graph2.vertices:
            if vertex2.terminal:
                # graph2 has terminals, so graph1 also must have terminals
                hasTerminals = True
                break
        else:
            vertex2 = self.graph2.vertices[0]
            
        for vertex1 in self.graph1.vertices:
            # If terminals are available, then skip vertices in the first
            # graph that are not terminals
            if hasTerminals and not vertex1.terminal: continue
            # Propose a pairing
            if self.feasible(vertex1, vertex2):
                # Add proposed match to mapping
                self.addToMapping(vertex1, vertex2)
                # Recurse
                isMatch = self.match(callDepth-1)
                if isMatch and not self.findAll:
                    return True
                # Undo proposed match
                self.removeFromMapping(vertex1, vertex2)
                
        # None of the proposed matches led to a complete isomorphism, so return False
        return False     
        
    cdef bint feasible(self, Vertex vertex1, Vertex vertex2) except -2:
        """
        Return ``True`` if vertex `vertex1` from the first graph is a feasible
        match for vertex `vertex2` from the second graph, or ``False`` if not.
        The semantic and structural relationship of the vertices is evaluated,
        including several structural "look-aheads" that cheaply eliminate many
        otherwise feasible pairs.
        """
        cdef Vertex vert1, vert2
        cdef Edge edge1, edge2
        cdef int term1Count, term2Count, neither1Count, neither2Count
        
        if not self.subgraph:
            # To be feasible the connectivity values must be an exact match
            if vertex1.connectivity != vertex2.connectivity: return False
        
        # Semantic check #1: vertex1 and vertex2 must be equivalent
        if self.subgraph:
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
                if self.subgraph:
                    if not edge1.isSpecificCaseOf(edge2): return False
                else:
                    if not edge1.equivalent(edge2): return False

        # There could still be edges in graph1 that aren't in graph2; this is okay
        # for subgraph matching, but not for exact matching
        if not self.subgraph:
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
        if self.subgraph:
            if neither1Count < neither2Count: return False
        else:
            if neither1Count != neither2Count: return False

        # Level 1 look-ahead: the number of adjacent vertices of vertex1 and
        # vertex2 that are terminals must be equal
        if self.subgraph:
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
        if not self.subgraph:
            for vert1 in vertex1.edges:
                if vert1.mapping is not None:
                    if vert1.mapping not in vertex2.edges: return False

        # All of our tests have been passed, so the two vertices are a feasible pair
        return True
    
    cdef addToMapping(self, Vertex vertex1, Vertex vertex2):
        """
        Add as valid a mapping of vertex `vertex1` from the first graph to
        vertex `vertex2` from the second graph, and update the terminals
        status accordingly.        
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
    
    cdef removeFromMapping(self, Vertex vertex1, Vertex vertex2):
        """
        Remove as valid a mapping of vertex `vertex1` from the first graph to
        vertex `vertex2` from the second graph, and update the terminals
        status accordingly.        
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
            