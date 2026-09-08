###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2026 Prof. William H. Green (whgreen@mit.edu),           #
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
algorithm of Vento and Foggia.  https://doi.org/10.1109/TPAMI.2004.75
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
                             bint strict=True, bint check_labels=False) except -2:
        """
        Return ``True`` if graph `graph1` is isomorphic to graph `graph2` with
        the optional initial mapping `initial_mapping`, or ``False`` otherwise.
        """
        self.isomorphism(graph1, graph2, initial_mapping, False, False, False, False, save_order=save_order, strict=strict,
                          check_labels=check_labels)
        return self.is_match

    cpdef list find_isomorphism(self, Graph graph1, Graph graph2, dict initial_mapping, bint save_order=False,
                                bint strict=True, bint check_labels=False):
        """
        Return a list of dicts of all valid isomorphism mappings from graph
        `graph1` to graph `graph2` with the optional initial mapping
        `initial_mapping`. If no valid isomorphisms are found, an empty list is
        returned.
        """
        self.isomorphism(graph1, graph2, initial_mapping, False, True, False, False, save_order=save_order, strict=strict,
                          check_labels=check_labels)
        return self.mapping_list

    cpdef bint is_subgraph_isomorphic(self, Graph graph1, Graph graph2, dict initial_mapping,
                                      bint save_order=False, bint check_labels=False) except -2:
        """
        Return ``True`` if graph `graph1` is subgraph isomorphic to subgraph
        `graph2` with the optional initial mapping `initial_mapping`, or
        ``False`` otherwise.
        """
        self.isomorphism(graph1, graph2, initial_mapping, True, False, False, False, save_order, strict=True, check_labels=check_labels)
        return self.is_match

    cpdef list find_subgraph_isomorphisms(self, Graph graph1, Graph graph2, dict initial_mapping, bint save_order=False,
                                          bint check_labels=False):
        """
        Return a list of dicts of all valid subgraph isomorphism mappings from
        graph `graph1` to subgraph `graph2` with the optional initial mapping
        `initial_mapping`. If no valid subgraph isomorphisms are found, an empty
        list is returned.
        """
        self.isomorphism(graph1, graph2, initial_mapping, True, True, False, False, save_order, strict=True, check_labels=check_labels)
        return self.mapping_list

    cpdef bint is_intersection_isomorphic(self, Graph graph1, Graph graph2, dict initial_mapping,
                                      bint save_order=False, bint check_labels=False) except -2:
        """
        Return ``True`` if subgraph `graph1` is intersection isomorphic to subgraph
        `graph2` with the optional initial mapping `initial_mapping`, or
        ``False`` otherwise.
        """
        self.isomorphism(graph1, graph2, initial_mapping, False, False, True, False, save_order=save_order, strict=True, check_labels=check_labels)
        return self.is_match

    cpdef list find_intersection_isomorphisms(self, Graph graph1, Graph graph2, dict initial_mapping, bint save_order=False, bint check_labels=False):
        """
        Return a list of dicts of all valid intersection isomorphism mappings from
        subgraph `graph1` to subgraph `graph2` with the optional initial mapping
        `initial_mapping`. If no valid intersection isomorphisms are found, an empty
        list is returned.
        """
        self.isomorphism(graph1, graph2, initial_mapping, False, True, True, False, save_order=save_order, strict=True, check_labels=check_labels)
        return self.mapping_list

    cpdef list find_largest_incomplete_isomorphisms(self, Graph graph1, Graph graph2, dict initial_mapping,
                                                     bint save_order=False, bint check_labels=False, bint find_all=False):
        """
        Find the largest common (non-induced) subgraph between `graph1` and `graph2`: the largest
        partial mapping from a subset of `graph2`'s vertices into `graph1` such that every edge of
        `graph2` between two mapped vertices has a corresponding edge in `graph1`, using the same
        vertex/edge compatibility rules as subgraph isomorphism (`graph1` may have extra vertices
        and extra edges beyond what `graph2` requires; unlike subgraph isomorphism, `graph2` itself
        need not be fully covered).

        Returns a list of dicts mapping `graph1` vertices to `graph2` vertices, containing only the
        matched subset (an unmatched `graph2` vertex is simply absent from every mapping's values).
        If no vertices could be matched at all, a list containing a single empty dict is returned.

        By default (`find_all=False`) only one mapping achieving the largest size is returned. If
        `find_all` is ``True``, every mapping achieving that size is returned instead -- note this
        can be very large for graphs with substantial symmetry (e.g. a periodic crystal lattice),
        since many differently-labeled mappings can all achieve the same maximum coverage.
        """
        self.isomorphism(graph1, graph2, initial_mapping, False, find_all, False, incomplete=True,
                         save_order=save_order, strict=True, check_labels=check_labels)
        return self.mapping_list

    cdef isomorphism(self, Graph graph1, Graph graph2, dict initial_mapping, bint subgraph, bint find_all, bint intersection,
                     bint incomplete, bint save_order=False, bint strict=True, bint check_labels=False):
        """
        Evaluate the isomorphism relationship between graphs `graph1` and
        `graph2` with optional initial mapping `initial_mapping`. If `subgraph`
        is ``True``, `graph2` is treated as a possible subgraph of `graph1`.
        If `find_all` is ``True``, all isomorphisms are found; otherwise only
        the first is found. If `check_labels` is ``True``, vertices only match
        if their `label` attributes also match.
        """
        cdef int call_depth, index1, index2

        if self.graph1 is not graph1:
            self.graph1 = graph1
            graph1.sort_vertices(save_order)

        if self.graph2 is not graph2:
            self.graph2 = graph2
            graph2.sort_vertices(save_order)

        self.initial_mapping = initial_mapping
        # Incomplete (largest common subgraph) search reuses subgraph mode's feasibility rules
        # (graph1/host may have extra vertices and extra edges; graph2/pattern's edges among
        # matched vertices must be realized in graph1) -- it just no longer requires graph2 to be
        # fully covered.
        self.subgraph = subgraph or incomplete
        self.find_all = find_all
        self.intersection = intersection
        self.incomplete = incomplete
        self.strict = strict
        self.check_labels = check_labels

        # Clear previous result
        self.is_match = False
        self.mapping_list = []
        self.best_size = -1
        self.max_possible_size = 0

        if self.incomplete:
            # Neither graph need be fully covered, so the usual size fast-rejects don't apply.
            # The best achievable coverage can't exceed the number of (non-ignored) vertices on
            # either side.
            self.max_possible_size = min(
                sum([1 for vertex1 in graph1.vertices if not vertex1.ignore]),
                sum([1 for vertex2 in graph2.vertices if not vertex2.ignore]),
            )
            if self.max_possible_size == 0:
                return
        elif not self.intersection and not self.subgraph and len(graph2.vertices) != len(graph1.vertices):
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
            vertex2.excluded = False
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
            vertex2.excluded = False

    cdef bint match(self, int call_depth) except -2:
        """
        Recursively search for pairs of vertices to match, until all vertices
        are matched or the viable set of matches is exhausted. The `call_depth`
        parameter helps ensure we never enter an infinite loop.
        """
        cdef Vertex vertex1, vertex2
        cdef dict mapping
        cdef bint has_terminals
        cdef int matched_so_far

        # The call depth should never be negative!
        if call_depth < 0:
            raise VF2Error('Negative call depth encountered in VF2_match().')

        if self.incomplete:
            # Branch-and-bound: this branch can decide the fate of at most `call_depth` more
            # graph2 vertices, so it can never do better than matched_so_far + call_depth. If that
            # can't beat (or, when not enumerating every tied mapping, even tie) the best found so
            # far, there's no point exploring it further.
            matched_so_far = self.count_mapped(self.graph2)
            if matched_so_far + call_depth < self.best_size or \
               (matched_so_far + call_depth == self.best_size and not self.find_all):
                return False

        # Done if we have mapped to all vertices in graph
        if call_depth == 0:
            if self.incomplete:
                self.record_incomplete_leaf(self.graph2)
                if not self.find_all and self.best_size == self.max_possible_size:
                    # Best possible coverage reached; no other branch can do better, so stop.
                    return True
                return False
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
            if self.incomplete and vertex2.excluded:
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
            # Take first unmapped (and, for an incomplete search, not-yet-excluded) vertex
            for vertex2 in self.graph2.vertices:
                if vertex2.mapping is None and not (self.incomplete and vertex2.excluded):
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

        if self.incomplete:
            # Also try leaving vertex2 out of the match entirely, so a single unmatchable vertex
            # doesn't fail the whole search -- tried last so that real matches (which are more
            # likely to extend the best mapping found so far) are explored, and can start pruning
            # other branches, before this one.
            vertex2.excluded = True
            is_match = self.match(call_depth - 1)
            vertex2.excluded = False
            if is_match and not self.find_all:
                return True

        # None of the proposed matches led to a complete isomorphism, so return False
        return False

    cdef int count_mapped(self, Graph graph2) except -1:
        """
        Count the vertices of `graph2` (excluding any flagged `ignore`) that are currently mapped.
        Used by the incomplete/largest-common-subgraph search to score and bound partial matches.
        """
        cdef Vertex vertex2
        cdef int count

        count = 0
        for vertex2 in graph2.vertices:
            if vertex2.ignore:
                continue
            if vertex2.mapping is not None:
                count += 1
        return count

    cdef record_incomplete_leaf(self, Graph graph2):
        """
        Called once every vertex of `graph2` has either been mapped or explicitly excluded from
        the match. If this leaf's mapping is larger than the best one found so far, it becomes the
        new best (replacing any previously recorded mappings); if it ties the best one found so far
        and `find_all` was requested, it's recorded alongside it.
        """
        cdef Vertex vertex2
        cdef dict mapping
        cdef int size

        size = self.count_mapped(graph2)
        if size <= self.best_size and not (size == self.best_size and self.find_all):
            return

        mapping = {}
        for vertex2 in graph2.vertices:
            if vertex2.ignore or vertex2.mapping is None:
                continue
            mapping[vertex2.mapping] = vertex2

        if size > self.best_size:
            self.best_size = size
            self.mapping_list = [mapping]
        else:
            self.mapping_list.append(mapping)
        self.is_match = True

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

        if not self.subgraph and not self.intersection:
            # To be feasible the connectivity values must be an exact match
            if vertex1.connectivity1 != vertex2.connectivity1: return False
            if vertex1.connectivity2 != vertex2.connectivity2: return False
            if vertex1.connectivity3 != vertex2.connectivity3: return False

        # Semantic check #1: vertex1 and vertex2 must be equivalent
        if self.subgraph:
            if not vertex1.is_specific_case_of(vertex2, check_labels=self.check_labels): return False
        elif self.intersection:
            if not vertex1.has_intersection_with(vertex2, check_labels=self.check_labels): return False
        else:
            if not vertex1.equivalent(vertex2, strict=self.strict, check_labels=self.check_labels): return False

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
                    elif self.intersection:
                        if not edge1.has_intersection_with(edge2): return False
                    else:
                        if not edge1.equivalent(edge2): return False

        # There could still be edges in graph1 that aren't in graph2; this is okay
        # for subgraph matching, but not for exact matching
        if not self.subgraph and not self.intersection:
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
        # vertex2 that are non-terminals must be equal. Skipped for an incomplete
        # (largest common subgraph) search: it assumes every one of vertex2's terminal
        # neighbors will eventually need a home in vertex1's neighborhood, which isn't true when
        # graph2 need not be fully covered -- some of them may end up excluded from the match
        # instead, so vertex1 doesn't need the capacity for all of them.
        if self.incomplete:
            pass
        elif self.subgraph:
            if neither1_count < neither2_count: return False
        elif not self.intersection:
            if neither1_count != neither2_count: return False

        # Level 1 look-ahead: the number of adjacent vertices of vertex1 and
        # vertex2 that are terminals must be equal. Skipped for the same reason as above.
        if self.incomplete:
            pass
        elif self.subgraph:
            if term1_count < term2_count: return False
        elif not self.intersection:
            if term1_count != term2_count: return False

        # Level 0 look-ahead: all adjacent vertices of vertex2 already in the
        # mapping must map to adjacent vertices of vertex1
        if not self.intersection:
            for vert2 in vertex2.edges:
                if vert2.mapping is not None:
                    if vert2.mapping not in vertex1.edges: return False
        # Also, all adjacent vertices of vertex1 already in the mapping must map to
        # adjacent vertices of vertex2, unless we are subgraph matching
        if not self.subgraph and not self.intersection:
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
