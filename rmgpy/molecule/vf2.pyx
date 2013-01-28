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
