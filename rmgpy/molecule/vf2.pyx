###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2019 Prof. William H. Green (whgreen@mit.edu),           #
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

"""
This module contains graph ismorphism functions that implement the VF2
algorithm of Vento and Foggia.  http://dx.doi.org/10.1109/TPAMI.2004.75
"""

from rmgpy.exceptions import VF2Error
from rmgpy.molecule.graph cimport Graph

################################################################################

cdef class VF2:
    """
    An implementation of the second version of the Vento-Foggia (VF2) algorithm
    for graph and subgraph isomorphism.
    """
    def __init__(self, graphA=None, graphB=None):
        self.graph1 = graphA
        self.graph2 = graphB

    @property
    def graphA(self):
        return self.graph1

    @graphA.setter
    def graphA(self, value):
        self.graph1 = value
        self.graph1.sort_vertices()

    @property
    def graphB(self):
        return self.graph2

    @graphB.setter
    def graphB(self, value):
        self.graph2 = value
        self.graph2.sort_vertices()

    cpdef bint is_isomorphic(self, Graph graph1, Graph graph2, dict initial_mapping, bint save_order=False,
                             bint strict=True) except -2:
        """
        Return ``True`` if graph `graph1` is isomorphic to graph `graph2` with
        the optional initial mapping `initial_mapping`, or ``False`` otherwise.
        """
        self.isomorphism(graph1, graph2, initial_mapping, False, False, save_order=save_order, strict=strict)
        return self.is_match

    cpdef list find_isomorphism(self, Graph graph1, Graph graph2, dict initial_mapping, bint save_order=False,
                                bint strict=True):
        """
        Return a list of dicts of all valid isomorphism mappings from graph
        `graph1` to graph `graph2` with the optional initial mapping 
        `initial_mapping`. If no valid isomorphisms are found, an empty list is
        returned.
        """
        self.isomorphism(graph1, graph2, initial_mapping, False, True, save_order=save_order, strict=strict)
        return self.mapping_list

    cpdef bint is_subgraph_isomorphic(self, Graph graph1, Graph graph2, dict initial_mapping,
                                      bint save_order=False) except -2:
        """
        Return ``True`` if graph `graph1` is subgraph isomorphic to subgraph
        `graph2` with the optional initial mapping `initial_mapping`, or
        ``False`` otherwise.
        """
        self.isomorphism(graph1, graph2, initial_mapping, True, False, save_order)
        return self.is_match

    cpdef list find_subgraph_isomorphisms(self, Graph graph1, Graph graph2, dict initial_mapping, bint save_order=False):
        """
        Return a list of dicts of all valid subgraph isomorphism mappings from
        graph `graph1` to subgraph `graph2` with the optional initial mapping 
        `initial_mapping`. If no valid subgraph isomorphisms are found, an empty
        list is returned.
        """
        self.isomorphism(graph1, graph2, initial_mapping, True, True, save_order)
        return self.mapping_list

    cdef isomorphism(self, Graph graph1, Graph graph2, dict initial_mapping, bint subgraph, bint find_all,
                     bint save_order=False, bint strict=True):
        """
        Evaluate the isomorphism relationship between graphs `graph1` and
        `graph2` with optional initial mapping `initial_mapping`. If `subgraph`
        is ``True``, `graph2` is treated as a possible subgraph of `graph1`.
        If `find_all` is ``True``, all isomorphisms are found; otherwise only
        the first is found.
        """
        cdef int call_depth, index1, index2

        if self.graph1 is not graph1:
            self.graph1 = graph1
            graph1.sort_vertices(save_order)

        if self.graph2 is not graph2:
            self.graph2 = graph2
            graph2.sort_vertices(save_order)

        self.initial_mapping = initial_mapping
        self.subgraph = subgraph
        self.find_all = find_all
        self.strict = strict

        # Clear previous result
        self.is_match = False
        self.mapping_list = []

        # Some quick isomorphism checks based on graph sizes
        if not self.subgraph and len(graph2.vertices) != len(graph1.vertices):
            # The two graphs don't have the same number of vertices, so they
            # cannot be isomorphic
            return
        elif not self.subgraph and len(graph2.vertices) == len(graph1.vertices) == 0:
            # The two graphs don't have any vertices; this means they are
            # trivially isomorphic
            self.is_match = True
            return
        elif self.subgraph and len(graph2.vertices) > len(graph1.vertices):
            # The second graph has more vertices than the first, so it cannot be
            # a subgraph of the first
            return

        # Initialize call_depth with the size of the smallest graph
        # Each recursive call to match() will decrease it by one;
        # when the whole graph has been explored, it should reach 0
        # It should never go below zero!
        call_depth = len(graph2.vertices)

        # Initialize mapping by clearing any previous mapping information
        for vertex1 in graph1.vertices:
            vertex1.mapping = None
            vertex1.terminal = False
        for vertex2 in graph2.vertices:
            vertex2.mapping = None
            vertex2.terminal = False
        # Set the initial mapping if provided
        if self.initial_mapping is not None:
            for vertex1, vertex2 in self.initial_mapping.items():
                self.add_to_mapping(vertex1, vertex2)
            call_depth -= len(self.initial_mapping)

        self.match(call_depth)

        if save_order:
            graph1.restore_vertex_order()
            graph2.restore_vertex_order()

        # We're done, so clear the mappings to prevent downstream effects
        for vertex1 in graph1.vertices:
            vertex1.mapping = None
            vertex1.terminal = False
        for vertex2 in graph2.vertices:
            vertex2.mapping = None
            vertex2.terminal = False

    cdef bint match(self, int call_depth) except -2:
        """
        Recursively search for pairs of vertices to match, until all vertices
        are matched or the viable set of matches is exhausted. The `call_depth`
        parameter helps ensure we never enter an infinite loop.
        """
        cdef Vertex vertex1, vertex2
        cdef dict mapping
        cdef bint has_terminals

        # The call depth should never be negative!
        if call_depth < 0:
            raise VF2Error('Negative call depth encountered in VF2_match().')

        # Done if we have mapped to all vertices in graph
        if call_depth == 0:
            if self.find_all:
                mapping = {}
                for vertex2 in self.graph2.vertices:
                    if vertex2.ignore:
                        continue
                    assert vertex2.mapping is not None
                    assert vertex2.mapping.mapping is vertex2
                    mapping[vertex2.mapping] = vertex2
                self.mapping_list.append(mapping)
            self.is_match = True
            return True

        # Create list of pairs of candidates for inclusion in mapping
        """
        10.1109/TPAMI.2004.75 says:
        "The set P(s) will be made of all the node pairs (n,m),
        with n belonging to T1out(s) and m to T2out(s),
        unless one of these two sets is empty. In this case,
        the set P(s) is likewise obtained by considering
        T1in(s) and T2in(s), respectively."

        But: for us, bonds are not directional, so ignore Tin(s)
        and just use Tout(s) which is what we call "terminals".
        """
        has_terminals = False
        for vertex2 in self.graph2.vertices:
            if vertex2.ignore:
                continue
            if vertex2.terminal:
                # graph2 has terminals, so graph1 also must have terminals
                has_terminals = True
                break
        else:
            """
            "In presence of not connected graphs, for some state s,
            all of the above sets may be empty. In this case,
            the set of candidate pairs making up P(s) will be
            the set Pd(s) of all the pairs of nodes not contained
            neither in G1(s) nor in G2(s)."

            So: use nodes not yet mapped.
            """
            # Take first unmapped vertex
            for vertex2 in self.graph2.vertices:
                if vertex2.mapping is None:
                    break
            else:
                raise VF2Error("Still seeking candidate pairs but all nodes in graph2 are already mapped.")

        for vertex1 in self.graph1.vertices:
            if vertex1.ignore:
                continue
            # If terminals are available, then skip vertices in the first
            # graph that are not terminals
            if has_terminals and not vertex1.terminal:
                continue
            # Otherwise take any node that is not already matched
            if vertex1.mapping is not None:
                continue
            # Propose a pairing
            if self.feasible(vertex1, vertex2):
                # Add proposed match to mapping
                self.add_to_mapping(vertex1, vertex2)
                # Recurse
                is_match = self.match(call_depth - 1)
                if is_match and not self.find_all:
                    return True
                # Undo proposed match
                self.remove_from_mapping(vertex1, vertex2)

        # None of the proposed matches led to a complete isomorphism, so return False
        return False

    cpdef bint feasible(self, Vertex vertex1, Vertex vertex2) except -2:
        """
        Return ``True`` if vertex `vertex1` from the first graph is a feasible
        match for vertex `vertex2` from the second graph, or ``False`` if not.
        The semantic and structural relationship of the vertices is evaluated,
        including several structural "look-aheads" that cheaply eliminate many
        otherwise feasible pairs.
        """
        cdef Vertex vert1, vert2
        cdef Edge edge1, edge2
        cdef int term1_count, term2_count, neither1_count, neither2_count

        if not self.subgraph:
            # To be feasible the connectivity values must be an exact match
            if vertex1.connectivity1 != vertex2.connectivity1: return False
            if vertex1.connectivity2 != vertex2.connectivity2: return False
            if vertex1.connectivity3 != vertex2.connectivity3: return False

        # Semantic check #1: vertex1 and vertex2 must be equivalent
        if self.subgraph:
            if not vertex1.is_specific_case_of(vertex2): return False
        else:
            if not vertex1.equivalent(vertex2, strict=self.strict): return False

        # Semantic check #2: adjacent vertices to vertex1 and vertex2 that are
        # already mapped should be connected by equivalent edges
        for vert2 in vertex2.edges:
            if vert2.mapping is not None:
                vert1 = vert2.mapping
                if vert1 not in vertex1.edges:
                    # The vertices are joined in graph2, but not in graph1
                    return False
                if self.strict:
                    # Check that the edges are equivalent
                    # If self.strict=False, we only care that the edge exists
                    edge1 = vertex1.edges[vert1]
                    edge2 = vertex2.edges[vert2]
                    if self.subgraph:
                        if not edge1.is_specific_case_of(edge2): return False
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
        term1_count = term2_count = neither1_count = neither2_count = 0  # note that 0 is immutable
        for vert1 in vertex1.edges:
            if vert1.terminal:
                term1_count += 1
            elif vert1.mapping is not None:
                neither1_count += 1
        for vert2 in vertex2.edges:
            if vert2.terminal:
                term2_count += 1
            elif vert2.mapping is not None:
                neither2_count += 1

        # Level 2 look-ahead: the number of adjacent vertices of vertex1 and
        # vertex2 that are non-terminals must be equal
        if self.subgraph:
            if neither1_count < neither2_count: return False
        else:
            if neither1_count != neither2_count: return False

        # Level 1 look-ahead: the number of adjacent vertices of vertex1 and
        # vertex2 that are terminals must be equal
        if self.subgraph:
            if term1_count < term2_count: return False
        else:
            if term1_count != term2_count: return False

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

    cdef add_to_mapping(self, Vertex vertex1, Vertex vertex2):
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

    cdef remove_from_mapping(self, Vertex vertex1, Vertex vertex2):
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
