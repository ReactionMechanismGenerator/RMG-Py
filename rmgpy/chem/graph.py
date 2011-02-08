#!/usr/bin/python
# -*- coding: utf-8 -*-

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

"""
This module contains an implementation of a graph data structure (the 
:class:`Graph` class) and functions for manipulating that graph, including 
efficient isomorphism functions.
"""

import cython
import logging

################################################################################

class Vertex(object):
    """
    A base class for vertices in a graph. Contains several connectivity values
    useful for accelerating isomorphism searches, as proposed by
    `Morgan (1965) <http://dx.doi.org/10.1021/c160017a018>`_.

    ==================  ========================================================
    Attribute           Description
    ==================  ========================================================
    `connectivity1`     The number of nearest neighbors
    `connectivity2`     The sum of the neighbors' `connectivity1` values
    `connectivity3`     The sum of the neighbors' `connectivity2` values
    `sortingLabel`      An integer used to sort the vertices
    ==================  ========================================================

    """

    def __init__(self):
        self.resetConnectivityValues()

    def equivalent(self, other):
        """
        Return :data:`True` if two vertices `self` and `other` are semantically
        equivalent, or :data:`False` if not. You should reimplement this
        function in a derived class if your vertices have semantic information.
        """
        return True

    def isSpecificCaseOf(self, other):
        """
        Return ``True`` if `self` is semantically more specific than `other`,
        or ``False`` if not. You should reimplement this function in a derived
        class if your edges have semantic information.
        """
        return True

    def resetConnectivityValues(self):
        """
        Reset the cached structure information for this vertex.
        """
        self.connectivity1 = -1
        self.connectivity2 = -1
        self.connectivity3 = -1
        self.sortingLabel = -1

def getVertexConnectivityValue(vertex):
    """
    Return a value used to sort vertices prior to poposing candidate pairs in
    :meth:`__VF2_pairs`. The value returned is based on the vertex's
    connectivity values (and assumes that they are set properly).
    """
    return ( -256*vertex.connectivity1 - 16*vertex.connectivity2 - vertex.connectivity3 )

def getVertexSortingLabel(vertex):
    """
    Return a value used to sort vertices prior to poposing candidate pairs in
    :meth:`__VF2_pairs`. The value returned is based on the vertex's
    connectivity values (and assumes that they are set properly).
    """
    return vertex.sortingLabel

################################################################################

class Edge(object):
    """
    A base class for edges in a graph. This class does *not* store the vertex
    pair that comprises the edge; that functionality would need to be included
    in the derived class.
    """

    def __init__(self):
        pass

    def equivalent(self, other):
        """
        Return ``True`` if two edges `self` and `other` are semantically
        equivalent, or ``False`` if not. You should reimplement this
        function in a derived class if your edges have semantic information.
        """
        return True

    def isSpecificCaseOf(self, other):
        """
        Return ``True`` if `self` is semantically more specific than `other`,
        or ``False`` if not. You should reimplement this function in a derived
        class if your edges have semantic information.
        """
        return True

################################################################################

class Graph:
    """
    A graph data type. The vertices of the graph are stored in a list
    `vertices`; this provides a consistent traversal order. The edges of the
    graph are stored in a dictionary of dictionaries `edges`. A single edge can
    be accessed using ``graph.edges[vertex1][vertex2]`` or the :meth:`getEdge`
    method; in either case, an exception will be raised if the edge does not
    exist. All edges of a vertex can be accessed using ``graph.edges[vertex]``
    or the :meth:`getEdges` method.
    """

    def __init__(self, vertices=None, edges=None):
        self.vertices = vertices or []
        self.edges = edges or {}
        
    def addVertex(self, vertex):
        """
        Add a `vertex` to the graph. The vertex is initialized with no edges.
        """
        self.vertices.append(vertex)
        self.edges[vertex] = dict()
        return vertex

    def addEdge(self, vertex1, vertex2, edge):
        """
        Add an `edge` to the graph as an edge connecting the two vertices
        `vertex1` and `vertex2`.
        """
        self.edges[vertex1][vertex2] = edge
        self.edges[vertex2][vertex1] = edge
        return edge

    def getEdges(self, vertex):
        """
        Return a list of the edges involving the specified `vertex`.
        """
        return self.edges[vertex]

    def getEdge(self, vertex1, vertex2):
        """
        Returns the edge connecting vertices `vertex1` and `vertex2`.
        """
        return self.edges[vertex1][vertex2]

    def hasVertex(self, vertex):
        """
        Returns ``True`` if `vertex` is a vertex in the graph, or ``False`` if
        not.
        """
        return vertex in self.vertices

    def hasEdge(self, vertex1, vertex2):
        """
        Returns ``True`` if vertices `vertex1` and `vertex2` are connected
        by an edge, or ``False`` if not.
        """
        return vertex2 in self.edges[vertex1] if vertex1 in self.edges else False

    def removeVertex(self, vertex):
        """
        Remove `vertex` and all edges associated with it from the graph. Does
        not remove vertices that no longer have any edges as a result of this
        removal.
        """
        for vertex2 in self.vertices:
            if vertex2 is not vertex:
                if vertex in self.edges[vertex2]:
                    del self.edges[vertex2][vertex]
        del self.edges[vertex]
        self.vertices.remove(vertex)

    def removeEdge(self, vertex1, vertex2):
        """
        Remove the edge having vertices `vertex1` and `vertex2` from the graph.
        Does not remove vertices that no longer have any edges as a result of
        this removal.
        """
        del self.edges[vertex1][vertex2]
        del self.edges[vertex2][vertex1]

    def copy(self, deep=False):
        """
        Create a copy of the current graph. If `deep` is ``True``, a deep copy
        is made: copies of the vertices and edges are used in the new graph.
        If `deep` is ``False`` or not specified, a shallow copy is made: the
        original vertices and edges are used in the new graph.
        """
        other = cython.declare(Graph)
        other = Graph()
        for vertex in self.vertices:
            other.addVertex(vertex.copy() if deep else vertex)
        for vertex1 in self.vertices:
            for vertex2 in self.edges[vertex1]:
                if deep:
                    index1 = self.vertices.index(vertex1)
                    index2 = self.vertices.index(vertex2)
                    other.addEdge(other.vertices[index1], other.vertices[index2],
                    self.edges[vertex1][vertex2].copy())
                else:
                    other.addEdge(vertex1, vertex2, self.edges[vertex1][vertex2])
        return other

    def merge(self, other):
        """
        Merge two graphs so as to store them in a single Graph object.
        """

        # Create output graph
        new = cython.declare(Graph)
        new = Graph()

        # Add vertices to output graph
        for vertex in self.vertices:
            new.addVertex(vertex)
        for vertex in other.vertices:
            new.addVertex(vertex)

        # Add edges to output graph
        for v1 in self.vertices:
            for v2 in self.edges[v1]:
                new.edges[v1][v2] = self.edges[v1][v2]
        for v1 in other.vertices:
            for v2 in other.edges[v1]:
                new.edges[v1][v2] = other.edges[v1][v2]

        return new

    def split(self):
        """
        Convert a single Graph object containing two or more unconnected graphs
        into separate graphs.
        """

        # Create potential output graphs
        new1 = cython.declare(Graph)
        new2 = cython.declare(Graph)
        verticesToMove = cython.declare(list)
        index = cython.declare(cython.int)

        new1 = self.copy()
        new2 = Graph()

        if len(self.vertices) == 0:
            return [new1]

        # Arbitrarily choose last atom as starting point
        verticesToMove = [ self.vertices[-1] ]

        # Iterate until there are no more atoms to move
        index = 0
        while index < len(verticesToMove):
            for v2 in self.edges[verticesToMove[index]]:
                if v2 not in verticesToMove:
                    verticesToMove.append(v2)
            index += 1

        # If all atoms are to be moved, simply return new1
        if len(new1.vertices) == len(verticesToMove):
            return [new1]

        # Copy to new graph
        for vertex in verticesToMove:
            new2.addVertex(vertex)
        for v1 in verticesToMove:
            for v2, edge in new1.edges[v1].iteritems():
                new2.edges[v1][v2] = edge

        # Remove from old graph
        for v1 in new2.vertices:
            for v2 in new2.edges[v1]:
                if v1 in verticesToMove and v2 in verticesToMove:
                    del new1.edges[v1][v2]
        for vertex in verticesToMove:
            new1.removeVertex(vertex)

        new = [new2]
        new.extend(new1.split())
        return new

    def resetConnectivityValues(self):
        """
        Reset any cached connectivity information. Call this method when you
        have modified the graph.
        """
        vertex = cython.declare(Vertex)
        for vertex in self.vertices: vertex.resetConnectivityValues()
        
    def updateConnectivityValues(self):
        """
        Update the connectivity values for each vertex in the graph. These are
        used to accelerate the isomorphism checking.
        """

        cython.declare(count=cython.short, edges=dict)
        cython.declare(vertex1=Vertex, vertex2=Vertex)

        assert str(self.__class__) != 'chempy.molecule.Molecule' or not self.implicitHydrogens, "%s has implicit hydrogens" % self

        for vertex1 in self.vertices:
            count = len(self.edges[vertex1])
            vertex1.connectivity1 = count
        for vertex1 in self.vertices:
            count = 0
            edges = self.edges[vertex1]
            for vertex2 in edges: count += vertex2.connectivity1
            vertex1.connectivity2 = count
        for vertex1 in self.vertices:
            count = 0
            edges = self.edges[vertex1]
            for vertex2 in edges: count += vertex2.connectivity2
            vertex1.connectivity3 = count
        
    def sortVertices(self):
        """
        Sort the vertices in the graph. This can make certain operations, e.g.
        the isomorphism functions, much more efficient.
        """
        cython.declare(index=cython.int, vertex=Vertex)
        # Only need to conduct sort if there is an invalid sorting label on any vertex
        for vertex in self.vertices:
            if vertex.sortingLabel < 0: break
        else:
            return
        self.vertices.sort(key=getVertexConnectivityValue)
        for index, vertex in enumerate(self.vertices):
            vertex.sortingLabel = index

    def isIsomorphic(self, other, initialMap=None):
        """
        Returns :data:`True` if two graphs are isomorphic and :data:`False`
        otherwise. Uses the VF2 algorithm of Vento and Foggia.
        """
        return VF2_isomorphism(self, other, False, False, initialMap)[0]

    def findIsomorphism(self, other, initialMap=None):
        """
        Returns :data:`True` if `other` is subgraph isomorphic and :data:`False`
        otherwise, and the matching mapping.
        Uses the VF2 algorithm of Vento and Foggia.
        """
        return VF2_isomorphism(self, other, False, True, initialMap)

    def isSubgraphIsomorphic(self, other, initialMap=None):
        """
        Returns :data:`True` if `other` is subgraph isomorphic and :data:`False`
        otherwise. Uses the VF2 algorithm of Vento and Foggia.
        """
        return VF2_isomorphism(self, other, True, False, initialMap)[0]

    def findSubgraphIsomorphisms(self, other, initialMap=None):
        """
        Returns :data:`True` if `other` is subgraph isomorphic and :data:`False`
        otherwise. Also returns the lists all of valid mappings.

        Uses the VF2 algorithm of Vento and Foggia.
        """
        return VF2_isomorphism(self, other, True, True, initialMap)

    def isCyclic(self):
        """
        Return :data:`True` if one or more cycles are present in the structure
        and :data:`False` otherwise.
        """
        for vertex in self.vertices:
            if self.isVertexInCycle(vertex):
                return True
        return False

    def isVertexInCycle(self, vertex):
        """
        Return :data:`True` if `vertex` is in one or more cycles in the graph,
        or :data:`False` if not.
        """
        chain = cython.declare(list)
        chain = [vertex]
        return self.__isChainInCycle(chain)

    def isEdgeInCycle(self, vertex1, vertex2):
        """
        Return :data:`True` if the edge between vertices `vertex1` and `vertex2`
        is in one or more cycles in the graph, or :data:`False` if not.
        """
        cycle_list = self.getAllCycles(vertex1)
        for cycle in cycle_list:
            if vertex2 in cycle:
                return True
        return False

    def __isChainInCycle(self, chain):
        """
        Is the `chain` in a cycle?
        Returns True/False.
        Recursively calls itself
        """
        # Note that this function no longer returns the cycle; just True/False
        vertex2 = cython.declare(Vertex)
        edge = cython.declare(Edge)
        found = cython.declare(cython.bint)

        for vertex2, edge in self.edges[chain[-1]].iteritems():
            if vertex2 is chain[0] and len(chain) > 2:
                return True
            elif vertex2 not in chain:
                # make the chain a little longer and explore again
                chain.append(vertex2)
                found = self.__isChainInCycle(chain)
                if found: return True
                # didn't find a cycle down this path (-vertex2),
                # so remove the vertex from the chain
                chain.remove(vertex2)
        return False

    def getAllCycles(self, startingVertex):
        """
        Given a starting vertex, returns a list of all the cycles containing
        that vertex.
        """
        chain = cython.declare(list)
        cycleList = cython.declare(list)

        cycleList=list()
        chain = [startingVertex]

        #chainLabels=range(len(self.keys()))
        #print "Starting at %s in graph: %s"%(self.keys().index(startingVertex),chainLabels)

        cycleList = self.__exploreCyclesRecursively(chain, cycleList)
        return cycleList

    def __exploreCyclesRecursively(self, chain, cycleList):
        """
        Finds cycles by spidering through a graph.
        Give it a chain of atoms that are connected, `chain`,
        and a list of cycles found so far `cycleList`.
        If `chain` is a cycle, it is appended to `cycleList`.
        Then chain is expanded by one atom (in each available direction)
        and the function is called again. This recursively spiders outwards
        from the starting chain, finding all the cycles.
        """
        vertex2 = cython.declare(Vertex)
        edge = cython.declare(Edge)

        # chainLabels = cython.declare(list)
        # chainLabels=[self.keys().index(v) for v in chain]
        # print "found %d so far. Chain=%s"%(len(cycleList),chainLabels)

        for vertex2, edge in self.edges[chain[-1]].iteritems():
            # vertex2 will loop through each of the atoms
            # that are bonded to the last atom in the chain.
            if vertex2 is chain[0] and len(chain) > 2:
                # it is the first atom in the chain - so the chain IS a cycle!
                cycleList.append(chain[:])
            elif vertex2 not in chain:
                # make the chain a little longer and explore again
                chain.append(vertex2)
                cycleList = self.__exploreCyclesRecursively(chain, cycleList)
                # any cycles down this path (-vertex2) have now been found,
                # so remove the vertex from the chain
                chain.pop(-1)
        return cycleList

    def getSmallestSetOfSmallestRings(self):
        """
        Return a list of the smallest set of smallest rings in the graph. The
        algorithm implements was adapted from a description by Fan, Panaye,
        Doucet, and Barbu (doi: 10.1021/ci00015a002)

        B. T. Fan, A. Panaye, J. P. Doucet, and A. Barbu. "Ring Perception: A
        New Algorithm for Directly Finding the Smallest Set of Smallest Rings
        from a Connection Table." *J. Chem. Inf. Comput. Sci.* **33**,
        p. 657-662 (1993).
        """

        graph = cython.declare(Graph)
        done = cython.declare(cython.bint)
        verticesToRemove = cython.declare(list)
        cycleList = cython.declare(list)
        cycles = cython.declare(list)
        vertex = cython.declare(Vertex)
        rootVertex = cython.declare(Vertex)
        found = cython.declare(cython.bint)
        cycle = cython.declare(list)
        graphs = cython.declare(list)

        # Make a copy of the graph so we don't modify the original
        graph = self.copy()
        
        # Step 1: Remove all terminal vertices
        done = False
        while not done:
            verticesToRemove = []
            for vertex1, value in graph.edges.iteritems():
                if len(value) == 1: verticesToRemove.append(vertex1)
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

        ### also need to remove EDGES that are not in ring

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
                    elif len(graph.edges[vertex]) < len(graph.edges[rootVertex]):
                        rootVertex = vertex

                # Get all cycles involving the root vertex
                cycles = graph.getAllCycles(rootVertex)
                if len(cycles) == 0:
                    # this vertex is no longer in a ring.
                    # remove all its edges
                    neighbours = graph.edges[rootVertex].keys()[:]
                    for vertex2 in neighbours:
                        graph.removeEdge(rootVertex, vertex2)
                    # then remove it
                    graph.removeVertex(rootVertex)
                    #print("Removed vertex that's no longer in ring")
                    continue # (pick a new root Vertex)
#                   raise Exception('Did not find expected cycle!')

                # Keep the smallest of the cycles found above
                cycle = cycles[0]
                for c in cycles[1:]:
                    if len(c) < len(cycle):
                        cycle = c
                cycleList.append(cycle)

                # Remove from the graph all vertices in the cycle that have only two edges
                verticesToRemove = []
                for vertex in cycle:
                    if len(graph.edges[vertex]) <= 2:
                        verticesToRemove.append(vertex)
                if len(verticesToRemove) == 0:
                    # there are no vertices in this cycle that with only two edges

                    # Remove edge between root vertex and any one vertex it is connected to
                    graph.removeEdge(rootVertex, graph[rootVertex].keys()[0])
                else:
                    for vertex in verticesToRemove:
                        graph.removeVertex(vertex)

        return cycleList

################################################################################

def VF2_isomorphism(graph1, graph2, subgraph=False, findAll=False, initialMap=None):
    """
    Determines if two :class:`Graph` objects `graph1` and `graph2` are
    isomorphic. A number of options affect how the isomorphism check is
    performed:

    * If `subgraph` is ``True``, the isomorphism function will treat `graph2`
      as a subgraph of `graph1`. In this instance a   subgraph can either mean a
      smaller graph (i.e. fewer vertices and/or edges) or a less specific graph.

    * If `findAll` is ``True``, all valid isomorphisms will be found and
      returned; otherwise only the first valid isomorphism will be returned.

    * The `initialMap` parameter can be used to pass a previously-established
      mapping. This mapping will be preserved in all returned valid
      isomorphisms.

    The isomorphism algorithm used is the VF2 algorithm of Vento and Foggia.
    The function returns a boolean `isMatch` indicating whether or not one or
    more valid isomorphisms have been found, and a list `mapList` of the valid
    isomorphisms, each consisting of a dictionary mapping from vertices of
    `graph1` to corresponding vertices of `graph2`.
    """

    cython.declare(isMatch=cython.bint, map12List=list, map21List=list)
    cython.declare(terminals1=list, terminals2=list, callDepth=cython.int)
    cython.declare(vert=Vertex)

    map21List = list()

    # Some quick initial checks to avoid using the full algorithm if the
    # graphs are obviously not isomorphic (based on graph size)
    if not subgraph:
        if len(graph2.vertices) != len(graph1.vertices):
            # The two graphs don't have the same number of vertices, so they
            # cannot be isomorphic
            return False, map21List
        elif len(graph1.vertices) == len(graph2.vertices) == 0:
            logging.warning("Tried matching empty graphs (returning True)")
            # The two graphs don't have any vertices; this means they are
            # trivially isomorphic
            return True, map21List
    else:
        if len(graph2.vertices) > len(graph1.vertices):
            # The second graph has more vertices than the first, so it cannot be
            # a subgraph of the first
            return False, map21List

    if initialMap is None: initialMap = {}
    map12List = list()
    
    # Initialize callDepth with the size of the largest graph
    # Each recursive call to __VF2_match will decrease it by one;
    # when the whole graph has been explored, it should reach 0
    # It should never go below zero!
    callDepth = min(len(graph1.vertices), len(graph2.vertices)) - len(initialMap)

    # Sort the vertices in each graph to make the isomorphism more efficient
    graph1.sortVertices()
    graph2.sortVertices()

    # Generate initial mapping pairs
    #   map21 = map to 2 from 1
    #   map12 = map to 1 from 2
    map21 = initialMap
    map12 = dict([(v,k) for k,v in initialMap.iteritems()])
    
    # Generate an initial set of terminals
    terminals1 = __VF2_terminals(graph1, map21)
    terminals2 = __VF2_terminals(graph2, map12)

    isMatch = __VF2_match(graph1, graph2, map21, map12, \
        terminals1, terminals2, subgraph, findAll, map21List, map12List, callDepth)

    if findAll:
        return len(map21List) > 0, map21List
    else:
        return isMatch, map21

def __VF2_feasible(graph1, graph2, vertex1, vertex2, map21, map12, terminals1,
    terminals2, subgraph):
    """
    Returns :data:`True` if two vertices `vertex1` and `vertex2` from graphs
    `graph1` and `graph2`, respectively, are feasible matches. `mapping21` and
    `mapping12` are the current state of the mapping from `graph1` to `graph2`
    and vice versa, respectively. `terminals1` and `terminals2` are lists of
    the vertices that are directly connected to the already-mapped vertices.
    `subgraph` is :data:`True` if graph2 is to be treated as a potential
    subgraph of graph1. i.e. graph1 is a specific case of graph2.

    Uses the VF2 algorithm of Vento and Foggia. The feasibility is assessed
    through a series of semantic and structural checks. Only the combination
    of the semantic checks and the level 0 structural check are both
    necessary and sufficient to ensure feasibility. (This does *not* mean that
    vertex1 and vertex2 are always a match, although the level 1 and level 2
    checks preemptively eliminate a number of false positives.)
    """

    cython.declare(vert1=Vertex, vert2=Vertex, edge1=Edge, edge2=Edge, edges1=dict, edges2=dict)
    cython.declare(i=cython.int)
    cython.declare(term1Count=cython.int, term2Count=cython.int, neither1Count=cython.int, neither2Count=cython.int)

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

    # Get edges adjacent to each vertex
    edges1 = graph1.edges[vertex1]
    edges2 = graph2.edges[vertex2]

    # Semantic check #2: adjacent vertices to vertex1 and vertex2 that are
    # already mapped should be connected by equivalent edges
    for vert2 in edges2:
        if vert2 in map12:
            vert1 = map12[vert2]
            if not vert1 in edges1: # atoms not joined in graph1
                return False
            edge1 = edges1[vert1]
            edge2 = edges2[vert2]
            if subgraph:
                if not edge1.isSpecificCaseOf(edge2): return False
            else: # exact match required
                if not edge1.equivalent(edge2): return False

    # there could still be edges in graph1 that aren't in graph2.
    # this is ok for subgraph matching, but not for exact matching
    if not subgraph:
        for vert1 in edges1:
            if vert1 in map21:
                vert2 = map21[vert1]
                if not vert2 in edges2: return False

    # Count number of terminals adjacent to vertex1 and vertex2
    term1Count = 0; term2Count = 0; neither1Count = 0; neither2Count = 0

    for vert1 in edges1:
        if vert1 in terminals1: term1Count += 1
        elif vert1 not in map21: neither1Count += 1
    for vert2 in edges2:
        if vert2 in terminals2: term2Count += 1
        elif vert2 not in map12: neither2Count += 1

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
    for vert2 in edges2:
        if vert2 in map12:
            vert1 = map12[vert2]
            if vert1 not in edges1: return False
    # Also, all adjacent vertices of vertex1 already in the mapping must map to
    # adjacent vertices of vertex2, unless we are subgraph matching
    if not subgraph:
        for vert1 in edges1:
            if vert1 in map21:
                vert2 = map21[vert1]
                if vert2 not in edges2: return False

    # All of our tests have been passed, so the two vertices are a feasible
    # pair
    return True

def __VF2_match(graph1, graph2, map21, map12, terminals1, terminals2, subgraph,
    findAll, map21List, map12List, callDepth):
    """
    A recursive function used to explore two graphs `graph1` and `graph2` for
    isomorphism by attempting to map them to one another. `mapping21` and
    `mapping12` are the current state of the mapping from `graph1` to `graph2`
    and vice versa, respectively. `terminals1` and `terminals2` are lists of
    the vertices that are directly connected to the already-mapped vertices.
    `subgraph` is :data:`True` if graph2 is to be treated as a potential
    subgraph of graph1. i.e. graph1 is a specific case of graph2.

    If findAll=True then it adds valid mappings to map21List and
    map12List, but returns False when done (or True if the initial mapping is complete)

    Uses the VF2 algorithm of Vento and Foggia, which is O(N) in spatial complexity
    and O(N**2) (best-case) to O(N! * N) (worst-case) in temporal complexity.
    """

    cython.declare(vertices1=list, new_terminals1=list, new_terminals2=list)
    cython.declare(vertex1=Vertex, vertex2=Vertex)
    cython.declare(ismatch=cython.bint)

    # Make sure we don't get cause in an infinite recursive loop
    if callDepth < 0:
        logging.error("Recursing too deep. Now %d" % callDepth)
        if callDepth < -100:
            raise Exception("Recursing infinitely deep!")

    # Done if we have mapped to all vertices in graph
    if callDepth == 0:
        if not subgraph:
            assert len(map21) == len(graph1.vertices), \
                "Calldepth mismatch: callDepth = %g, len(map21) = %g, len(map12) = %g, len(graph1.vertices) = %g, len(graph2.vertices) = %g" % (callDepth, len(map21), len(map12), len(graph1.vertices), len(graph2.vertices))
            if findAll:
                map21List.append(map21.copy())
                map12List.append(map12.copy())
            return True
        else:
            assert len(map12) == len(graph2.vertices), \
                "Calldepth mismatch: callDepth = %g, len(map21) = %g, len(map12) = %g, len(graph1.vertices) = %g, len(graph2.vertices) = %g" % (callDepth, len(map21), len(map12), len(graph1.vertices), len(graph2.vertices))
            if findAll:
                map21List.append(map21.copy())
                map12List.append(map12.copy())
            return True

    # Create list of pairs of candidates for inclusion in mapping
    # Note that the extra Python overhead is not worth making this a standalone
    # method, so we simply put it inline here
    # If we have terminals for both graphs, then use those as a basis for the
    # pairs
    if len(terminals1) > 0 and len(terminals2) > 0:
        vertices1 = terminals1
        vertex2 = terminals2[0]
    # Otherwise construct list from all *remaining* vertices (not matched)
    else:
        # vertex2 is the lowest-labelled un-mapped vertex from graph2
        # Note that this assumes that graph2.vertices is properly sorted
        vertices1 = []
        for vertex1 in graph1.vertices:
            if vertex1 not in map21:
                vertices1.append(vertex1)
        for vertex2 in graph2.vertices:
            if vertex2 not in map12:
                break
        else:
            raise Exception("Could not find a pair to propose!")
    
    for vertex1 in vertices1:
        # propose a pairing
        if __VF2_feasible(graph1, graph2, vertex1, vertex2, map21, map12, \
                terminals1, terminals2, subgraph):
            # Update mapping accordingly
            map21[vertex1] = vertex2
            map12[vertex2] = vertex1

            # update terminals
            new_terminals1 = __VF2_updateTerminals(graph1, map21, terminals1, vertex1)
            new_terminals2 = __VF2_updateTerminals(graph2, map12, terminals2, vertex2)

            # Recurse
            ismatch = __VF2_match(graph1, graph2, \
                map21, map12, new_terminals1, new_terminals2, subgraph, findAll, \
                map21List, map12List, callDepth-1)
            if ismatch:
                if not findAll:
                    return True
            # Undo proposed match
            del map21[vertex1]
            del map12[vertex2]
            # changes to 'new_terminals' will be discarded and 'terminals' is unchanged

    return False

def __VF2_terminals(graph, mapping):
    """
    For a given graph `graph` and associated partial mapping `mapping`,
    generate a list of terminals, vertices that are directly connected to
    vertices that have already been mapped.

    List is sorted (using key=__getSortLabel) before returning.
    """

    cython.declare(terminals=list)
    terminals = list()
    for vertex2 in graph.vertices:
        if vertex2 not in mapping:
            for vertex1 in mapping:
                if vertex2 in graph.edges[vertex1]:
                    terminals.append(vertex2)
                    break
    return terminals

def __VF2_updateTerminals(graph, mapping, old_terminals, new_vertex):
    """
    For a given graph `graph` and associated partial mapping `mapping`,
    *updates* a list of terminals, vertices that are directly connected to
    vertices that have already been mapped. You have to pass it the previous
    list of terminals `old_terminals` and the vertex `vertex` that has been
    added to the mapping. Returns a new *copy* of the terminals.
    """

    cython.declare(terminals=list, vertex1=Vertex, vertex2=Vertex, edges=dict)
    cython.declare(i=cython.int, sorting_label=cython.short, sorting_label2=cython.short)

    # Copy the old terminals, leaving out the new_vertex
    terminals = old_terminals[:]
    if new_vertex in terminals: terminals.remove(new_vertex)

    # Add the terminals of new_vertex
    edges = graph.edges[new_vertex]
    for vertex1 in edges:
        if vertex1 not in mapping: # only add if not already mapped
            # find spot in the sorted terminals list where we should put this vertex
            sorting_label = vertex1.sortingLabel
            i=0; sorting_label2=-1 # in case terminals list empty
            for i in range(len(terminals)):
                vertex2 = terminals[i]
                sorting_label2 = vertex2.sortingLabel
                if sorting_label2 >= sorting_label:
                    break
                # else continue going through the list of terminals
            else: # got to end of list without breaking,
                # so add one to index to make sure vertex goes at end
                i+=1
            if sorting_label2 == sorting_label: # this vertex already in terminals.
                continue # try next vertex in graph[new_vertex]

            # insert vertex in right spot in terminals
            terminals.insert(i,vertex1)

    return terminals

################################################################################
