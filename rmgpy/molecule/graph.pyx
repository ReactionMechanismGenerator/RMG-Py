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
This module contains an implementation of a graph data structure (the 
:class:`Graph` class) and functions for manipulating that graph, including 
efficient isomorphism functions. This module also contains base classes for
the vertices and edges (:class:`Vertex` and :class:`Edge`, respectively) that
are the components of a graph.
"""

import itertools

from rmgpy.molecule.vf2 cimport VF2

################################################################################

cdef class Vertex(object):
    """
    A base class for vertices in a graph. Contains several connectivity values
    useful for accelerating isomorphism searches, as proposed by
    `Morgan (1965) <https://doi.org/10.1021/c160017a018>`_.

    =================== =============== ========================================
    Attribute           Type            Description
    =================== =============== ========================================
    `connectivity1`     ``int``         The number of nearest neighbors
    `connectivity2`     ``int``         The sum of the neighbors' `connectivity1` values
    `connectivity3`     ``int``         The sum of the neighbors' `connectivity2` values
    `edges`             ``dict``        Dictionary of edges with keys being neighboring vertices
    `sorting_label`      ``int``         An integer label used to sort the vertices
    =================== =============== ========================================
    
    """

    def __init__(self):
        self.edges = {}
        self.reset_connectivity_values()
        self.ignore = False
        self.excluded = False

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        d = {
            'edges': self.edges,
            'connectivity1': self.connectivity1,
            'connectivity2': self.connectivity2,
            'connectivity3': self.connectivity3,
            'sorting_label': self.sorting_label,
            'terminal': self.terminal,
            'mapping': self.mapping,
        }
        return (Vertex, (), d)

    def __setstate__(self, d):
        self.edges = d['edges']
        self.connectivity1 = d['connectivity1']
        self.connectivity2 = d['connectivity2']
        self.connectivity3 = d['connectivity3']
        self.sorting_label = d['sorting_label']
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

    cpdef bint equivalent(self, Vertex other, bint strict=True, bint check_labels=False) except -2:
        """
        Return :data:`True` if two vertices `self` and `other` are semantically
        equivalent, or :data:`False` if not. You should reimplement this
        function in a derived class if your vertices have semantic information.
        If `check_labels` is ``True``, subclasses with a `label` attribute
        should also require that the labels match.
        """
        return True

    cpdef bint is_specific_case_of(self, Vertex other, bint check_labels=False) except -2:
        """
        Return ``True`` if `self` is semantically more specific than `other`,
        or ``False`` if not. You should reimplement this function in a derived
        class if your edges have semantic information. If `check_labels` is
        ``True``, subclasses with a `label` attribute should also require that
        the labels match.
        """
        return True

    cpdef reset_connectivity_values(self):
        """
        Reset the cached structure information for this vertex.
        """
        self.connectivity1 = -1
        self.connectivity2 = -1
        self.connectivity3 = -1
        self.sorting_label = -1
        self.terminal = False
        self.mapping = None

cpdef short get_vertex_connectivity_value(Vertex vertex) except 1:
    """
    Return a value used to sort vertices prior to poposing candidate pairs in
    :meth:`__VF2_pairs`. The value returned is based on the vertex's
    connectivity values (and assumes that they are set properly).
    """
    return ( -256*vertex.connectivity1 - 16*vertex.connectivity2 - vertex.connectivity3 )

cpdef short get_vertex_sorting_label(Vertex vertex) except -1:
    """
    Return a value used to sort vertices prior to poposing candidate pairs in
    :meth:`__VF2_pairs`. The value returned is based on the vertex's
    connectivity values (and assumes that they are set properly).
    """
    return vertex.sorting_label

################################################################################

cdef class Edge(object):
    """
    A base class for edges in a graph. The vertices which comprise the edge can be
    accessed using the `vertex1` and `vertex2` attributes.
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

    cpdef bint is_specific_case_of(self, Edge other) except -2:
        """
        Return ``True`` if `self` is semantically more specific than `other`,
        or ``False`` if not. You should reimplement this function in a derived
        class if your edges have semantic information.
        """
        return True

    cpdef Vertex get_other_vertex(self, Vertex vertex):
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

cdef VF2 vf2 = VF2()

cdef Vertex _get_edge_vertex1(Edge edge):
    return edge.vertex1

cdef Vertex _get_edge_vertex2(Edge edge):
    return edge.vertex2

cdef class Graph(object):
    """
    A graph data type. The vertices of the graph are stored in a list
    `vertices`; this provides a consistent traversal order. A single edge can
    be accessed using the :meth:`get_edge` method or by accessing specific
    vertices using ``vertex1.edges[vertex2]``; in either case, an exception
    will be raised if the edge does not exist. All edges of a vertex can be
    accessed using the :meth:`get_edges` method or ``vertex.edges``.
    """

    def __init__(self, vertices=None):
        self.vertices = vertices or []
        self._relevant_cycles = None

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (Graph, (self.vertices,))

    cpdef Vertex add_vertex(self, Vertex vertex):
        """
        Add a `vertex` to the graph. The vertex is initialized with no edges.
        """
        self.vertices.append(vertex)
        vertex.edges = dict()
        self._relevant_cycles = None
        return vertex

    cpdef Edge add_edge(self, Edge edge):
        """
        Add an `edge` to the graph. The two vertices in the edge must already
        exist in the graph, or a :class:`ValueError` is raised.
        """
        if edge.vertex1 not in self.vertices or edge.vertex2 not in self.vertices:
            raise ValueError('Attempted to add edge between vertices not in the graph.')
        edge.vertex1.edges[edge.vertex2] = edge
        edge.vertex2.edges[edge.vertex1] = edge
        self._relevant_cycles = None
        return edge

    cpdef list get_all_edges(self):
        """
        Returns a list of all edges in the graph.
        """
        cdef set edge_set
        cdef Vertex vertex
        cdef Edge edge

        edge_set = set()
        for vertex in self.vertices:
            for edge in vertex.edges.values():
                edge_set.add(edge)

        return list(edge_set)

    cpdef dict get_edges(self, Vertex vertex):
        """
        Return a dictionary of the edges involving the specified `vertex`.
        """
        return vertex.edges

    cpdef Edge get_edge(self, Vertex vertex1, Vertex vertex2):
        """
        Returns the edge connecting vertices `vertex1` and `vertex2`.
        """
        try:
            return vertex1.edges[vertex2]
        except KeyError:
            raise ValueError('The specified vertices are not connected by an edge in this graph.')

    cpdef bint has_vertex(self, Vertex vertex) except -2:
        """
        Returns ``True`` if `vertex` is a vertex in the graph, or ``False`` if
        not.
        """
        return vertex in self.vertices

    cpdef bint has_edge(self, Vertex vertex1, Vertex vertex2) except -2:
        """
        Returns ``True`` if vertices `vertex1` and `vertex2` are connected
        by an edge, or ``False`` if not.
        """
        return vertex1 in self.vertices and vertex2 in vertex1.edges

    cpdef remove_vertex(self, Vertex vertex):
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
        self._relevant_cycles = None

    cpdef remove_edge(self, Edge edge):
        """
        Remove the specified `edge` from the graph.
        Does not remove vertices that no longer have any edges as a result of
        this removal.
        """
        del edge.vertex1.edges[edge.vertex2]
        del edge.vertex2.edges[edge.vertex1]
        self._relevant_cycles = None

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
                vertex2 = other.add_vertex(vertex.copy())
                mapping[vertex] = vertex2
            else:
                edges = vertex.edges
                other.add_vertex(vertex)
                vertex.edges = edges
        if deep:
            for vertex1 in vertices:
                for vertex2 in vertex1.edges:
                    edge = vertex1.edges[vertex2]
                    edge = edge.copy()
                    edge.vertex1 = mapping[vertex1]
                    edge.vertex2 = mapping[vertex2]
                    other.add_edge(edge)
        return other

    cpdef dict copy_and_map(self):
        """
        Create a deep copy of the current graph, and return the dict
        'mapping'. Method was modified from Graph.copy() method
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
            vertex2 = other.add_vertex(vertex.copy())
            mapping[vertex] = vertex2

        for vertex1 in vertices:
            for vertex2 in vertex1.edges:
                edge = vertex1.edges[vertex2]
                edge = edge.copy()
                edge.vertex1 = mapping[vertex1]
                edge.vertex2 = mapping[vertex2]
                other.add_edge(edge)
        return mapping

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
            new.add_vertex(vertex)
            vertex.edges = edges

        if self is other:
            other = other.copy(deep=True)

        for vertex in other.vertices:
            edges = vertex.edges
            new.add_vertex(vertex)
            vertex.edges = edges

        return new

    cpdef list split(self):
        """
        Convert a single Graph object containing two or more unconnected graphs
        into separate graphs.
        """
        cdef Graph new1, new2
        cdef Vertex vertex, vertex1, vertex2
        cdef list vertices_to_move
        cdef int index

        # Create potential output graphs
        new1 = self.copy()
        new2 = Graph()

        if len(self.vertices) == 0:
            return [new1]

        # Arbitrarily choose last atom as starting point
        vertices_to_move = [self.vertices[-1]]

        # Iterate until there are no more atoms to move
        index = 0
        while index < len(vertices_to_move):
            for v2 in vertices_to_move[index].edges:
                if v2 not in vertices_to_move:
                    vertices_to_move.append(v2)
            index += 1

        # If all atoms are to be moved, simply return new1
        if len(new1.vertices) == len(vertices_to_move):
            return [new1]

        # Copy to new graph and remove from old graph
        for vertex in vertices_to_move:
            new2.vertices.append(vertex)
            new1.vertices.remove(vertex)

        new = [new2]
        new.extend(new1.split())
        return new

    cpdef reset_connectivity_values(self):
        """
        Reset any cached connectivity information. Call this method when you
        have modified the graph.
        """
        cdef Vertex vertex
        for vertex in self.vertices: vertex.reset_connectivity_values()

    cpdef update_connectivity_values(self):
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

    cpdef sort_vertices(self, bint save_order=False):
        """
        Sort the vertices in the graph. This can make certain operations, e.g.
        the isomorphism functions, much more efficient.
        """
        cdef Vertex vertex
        cdef int index

        if save_order:
            self.ordered_vertices = self.vertices[:]

        # Only need to conduct sort if there is an invalid sorting label on any vertex
        for vertex in self.vertices:
            if vertex.sorting_label < 0: break
        else:
            return
        # If we need to sort then let's also update the connecitivities so
        # we're sure they are right, since the sorting labels depend on them
        self.update_connectivity_values()
        self.vertices.sort(key=get_vertex_connectivity_value)
        for index, vertex in enumerate(self.vertices):
            vertex.sorting_label = index

    cpdef restore_vertex_order(self):
        """
        reorder the vertices to what they were before sorting
        if you saved the order
        """
        if not self.ordered_vertices or len(self.vertices) != len(self.ordered_vertices):
            raise ValueError('Number of vertices has changed cannot restore original vertex order')
        else:
            self.vertices = self.ordered_vertices

    cpdef bint is_isomorphic(self, Graph other, dict initial_map=None, bint generate_initial_map=False, bint save_order=False, bint strict=True, bint check_labels=False) except -2:
        """
        Returns :data:`True` if two graphs are isomorphic and :data:`False`
        otherwise. Uses the VF2 algorithm of Vento and Foggia.

        Args:
            initial_map (dict, optional): initial atom mapping to use
            generate_initial_map (bool, optional): if ``True``, initialize map by pairing atoms with same labels
            save_order (bool, optional):  if ``True``, reset atom order after performing atom isomorphism
            strict (bool, optional):     if ``False``, perform isomorphism ignoring electrons
            check_labels (bool, optional): if ``True``, atoms only match if their `label` attributes match
        """
        if generate_initial_map:
            initial_map = dict()
            for atom in self.vertices:
                if atom.label and atom.label != '':
                    for a in other.vertices:
                        if a.label == atom.label:
                            initial_map[atom] = a
                            break
                    else:
                        return False
            if not self.is_mapping_valid(other, initial_map, equivalent=True, strict=True, check_labels=check_labels):
                return False

        return vf2.is_isomorphic(self, other, initial_map, save_order=save_order, strict=strict, check_labels=check_labels)

    cpdef list find_isomorphism(self, Graph other, dict initial_map=None, bint save_order=False, bint strict=True, bint check_labels=False):
        """
        Returns :data:`True` if `other` is subgraph isomorphic and :data:`False`
        otherwise, and the matching mapping.
        Uses the VF2 algorithm of Vento and Foggia.

        Args:
            initial_map (dict, optional): initial atom mapping to use
            save_order (bool, optional):  if ``True``, reset atom order after performing atom isomorphism
            strict (bool, optional):     if ``False``, perform isomorphism ignoring electrons
            check_labels (bool, optional): if ``True``, atoms only match if their `label` attributes match
        """
        return vf2.find_isomorphism(self, other, initial_map, save_order=save_order, strict=strict, check_labels=check_labels)

    cpdef bint is_subgraph_isomorphic(self, Graph other, dict initial_map=None, bint save_order=False, bint check_labels=False) except -2:
        """
        Returns :data:`True` if `other` is subgraph isomorphic and :data:`False`
        otherwise. Uses the VF2 algorithm of Vento and Foggia.
        """
        return vf2.is_subgraph_isomorphic(self, other, initial_map, save_order=save_order, check_labels=check_labels)

    cpdef list find_subgraph_isomorphisms(self, Graph other, dict initial_map=None, bint save_order=False, bint check_labels=False):
        """
        Returns :data:`True` if `other` is subgraph isomorphic and :data:`False`
        otherwise. Also returns the lists all of valid mappings.

        Uses the VF2 algorithm of Vento and Foggia.
        """
        return vf2.find_subgraph_isomorphisms(self, other, initial_map, save_order=save_order, check_labels=check_labels)

    cpdef bint is_intersection_isomorphic(self, Graph other, dict initial_map=None, bint save_order=False, bint check_labels=False) except -2:
        """
        Returns :data:`True` if `other` is intersection isomorphic and :data:`False`
        otherwise. Uses the VF2 algorithm of Vento and Foggia.
        """
        return vf2.is_intersection_isomorphic(self, other, initial_map, save_order=save_order, check_labels=check_labels)

    cpdef list find_intersection_isomorphisms(self, Graph other, dict initial_map=None, bint save_order=False, bint check_labels=False):
        """
        Returns :data:`True` if `other` is intersection isomorphic and :data:`False`
        otherwise. Also returns the lists all of valid mappings.

        Uses the VF2 algorithm of Vento and Foggia.
        """
        return vf2.find_intersection_isomorphisms(self, other, initial_map, save_order=save_order, check_labels=check_labels)

    cpdef list find_largest_incomplete_isomorphisms(self, Graph other, dict initial_map=None, bint save_order=False, bint check_labels=False, bint find_all=False):
        """
        Find the largest common (non-induced) subgraph between `self` and `other`: the largest
        partial mapping from a subset of `other`'s vertices into `self` such that every edge of
        `other` between two mapped vertices has a corresponding edge in `self`. `self` may have
        extra vertices and extra edges beyond what `other` requires; `other` need not be fully
        covered. See :meth:`VF2.find_largest_incomplete_isomorphisms` for details.

        Uses the VF2 algorithm of Vento and Foggia.
        """
        return vf2.find_largest_incomplete_isomorphisms(self, other, initial_map, save_order=save_order, check_labels=check_labels, find_all=find_all)

    cpdef bint is_cyclic(self) except -2:
        """
        Return ``True`` if one or more cycles are present in the graph or
        ``False`` otherwise.
        """
        return len(self._get_relevant_cycles()) > 0

    cpdef bint is_vertex_in_cycle(self, Vertex vertex) except -2:
        """
        Return ``True`` if the given `vertex` is contained in one or more
        cycles in the graph, or ``False`` if not.
        """
        cdef list cycle
        for cycle in self._get_relevant_cycles():
            if vertex in cycle:
                return True
        return False

    cpdef bint is_edge_in_cycle(self, Edge edge) except -2:
        """
        Return :data:`True` if the edge between vertices `vertex1` and `vertex2`
        is in one or more cycles in the graph, or :data:`False` if not.
        """
        cdef list cycle
        cdef int i, n

        for cycle in self._get_relevant_cycles():
            n = len(cycle)
            for i in range(n):
                if ((cycle[i] is edge.vertex1 and cycle[i - 1] is edge.vertex2) or
                        (cycle[i] is edge.vertex2 and cycle[i - 1] is edge.vertex1)):
                    return True
        return False

    cdef list _get_relevant_cycles(self):
        """
        Return the graph's *relevant cycles* -- the union of the edge sets of every minimum-weight
        cycle basis -- computed and cached on first use (invalidated by add_vertex/add_edge/
        remove_vertex/remove_edge, matching how the cache below is populated).

        This is what every other cycle/ring-membership method on this class is now built on: a
        single polynomial-time computation shared across all of them, rather than each doing its
        own unbounded recursive search. See the module-level get_relevant_cycles() function below
        for the algorithm.

        Each returned cycle is a list of Vertex objects in ring-traversal order (consecutive
        entries -- including the wraparound from the last back to the first -- are bonded).
        """
        if self._relevant_cycles is None:
            self._relevant_cycles = get_relevant_cycles(self)
        return self._relevant_cycles

    cpdef list get_all_cyclic_vertices(self):
        """
        Returns all vertices belonging to one or more cycles.
        """
        cdef list cyclic_vertices, cycle
        cdef set seen
        cdef Vertex vertex

        seen = set()
        for cycle in self._get_relevant_cycles():
            seen.update(cycle)
        # Preserve self.vertices order, matching the original vertex-by-vertex-scan behavior
        cyclic_vertices = [vertex for vertex in self.vertices if vertex in seen]
        return cyclic_vertices

    cpdef list get_all_cycles(self, Vertex starting_vertex):
        """
        Given a starting vertex, returns a list of all the cycles containing
        that vertex.

        This function returns a duplicate of each cycle to preserve prior behavior
        where [0,1,2,3]
        is counted as separate from [0,3,2,1]
        """
        cdef list result, cycle
        result = []
        for cycle in self._get_relevant_cycles():
            if starting_vertex in cycle:
                result.append(cycle[:])
                result.append(list(reversed(cycle)))
        return result

    cpdef list get_all_cycles_of_size(self, int size):
        """
        Return a list of the all non-duplicate relevant rings with length 'size'.
        """
        cdef list cycle
        return [cycle[:] for cycle in self._get_relevant_cycles() if len(cycle) == size]

    cpdef list get_all_simple_cycles_of_size(self, int size):
        """
        Return a list of all non-duplicate monocyclic rings with length 'size'.

        Naive approach by eliminating polycyclic rings that are returned by
        ``getAllCyclicsOfSize``.
        """
        cdef list cycle_list
        cdef int i, internal_connectivity
        cdef Vertex vertex

        cycle_list = self.get_all_cycles_of_size(size)

        i = 0
        while i < len(cycle_list):
            for vertex in cycle_list[i]:
                internal_connectivity = sum([1 if vertex2 in cycle_list[i] else 0 for vertex2 in vertex.edges.keys()])
                if internal_connectivity > 2:
                    del cycle_list[i]
                    break
            else:
                i += 1

        return cycle_list

    cpdef list sort_cyclic_vertices(self, list vertices):
        """
        Given a list of vertices comprising a cycle, sort them such that adjacent
        entries in the list are connected to each other.
        Warning: Assumes that the cycle is elementary, ie. no bridges.
        """
        cdef list ordered
        cdef Vertex vertex

        ordered = [vertices.pop()]
        while vertices:
            for vertex in vertices:
                if vertex in ordered[-1].edges:
                    ordered.append(vertex)
                    vertices.remove(vertex)
                    break
            else:
                # No connected vertex was found
                raise RuntimeError('Could not sort cyclic vertices because '
                                   'not all vertices are connected to two '
                                   'other vertices in the input list.')

        if not self.has_edge(ordered[0], ordered[-1]):
            raise RuntimeError('Input vertices do not comprise a single cycle.')

        return ordered


    cpdef list get_largest_ring(self, Vertex vertex):
        """
        returns the largest ring containing vertex. This is typically
        useful for finding the longest path in a polycyclic ring, since
        the polycyclic rings returned from get_polycycles are not necessarily
        in order in the ring structure.
        """
        all_cycles = self.get_all_cycles(vertex)
        longest_cycle = []
        for cycle in all_cycles:
            if len(cycle) > len(longest_cycle):
                longest_cycle = cycle
        return longest_cycle

    cpdef bint is_mapping_valid(self, Graph other, dict mapping, bint equivalent=True, bint strict=True, bint check_labels=False) except -2:
        """
        Check that a proposed `mapping` of vertices from `self` to `other`
        is valid by checking that the vertices and edges involved in the
        mapping are mutually equivalent.  If equivalent is ``True`` it checks
        if atoms and edges are equivalent, if ``False`` it checks if they
        are specific cases of each other. If strict is ``True``, electrons
        and bond orders are considered, and ignored if ``False``. If
        check_labels is ``True``, atoms only match if their `label`
        attributes also match.
        """
        cdef Vertex vertex1, vertex2
        cdef list vertices1, vertices2
        cdef bint self_has_edge, other_has_edge
        cdef int i, j

        # Check that the mapped pairs of vertices compare True
        for vertex1, vertex2 in mapping.items():
            if equivalent:
                if not vertex1.equivalent(vertex2, strict=strict, check_labels=check_labels):
                    return False
            else:
                if not vertex1.is_specific_case_of(vertex2, check_labels=check_labels):
                    return False

        # Check that any edges connected mapped vertices are equivalent
        vertices1 = list(mapping.keys())
        vertices2 = list(mapping.values())
        for i in range(len(vertices1)):
            for j in range(i + 1, len(vertices1)):
                self_has_edge = self.has_edge(vertices1[i], vertices1[j])
                other_has_edge = other.has_edge(vertices2[i], vertices2[j])
                if self_has_edge and other_has_edge:
                    # Both graphs have the edge, so we must check it for equivalence
                    if strict:
                        edge1 = self.get_edge(vertices1[i], vertices1[j])
                        edge2 = other.get_edge(vertices2[i], vertices2[j])
                        if equivalent:
                            if not edge1.equivalent(edge2):
                                return False
                        else:
                            if not edge1.is_specific_case_of(edge2):
                                return False
                elif not equivalent and self_has_edge and not other_has_edge:
                    #in the subgraph case self can have edges other doesn't have
                    continue
                elif self_has_edge or other_has_edge:
                    # Only one of the graphs has the edge, so the mapping must be invalid
                    return False

        # If we're here then the vertices and edges compare True, so the
        # mapping is valid
        return True

    cpdef bint has_same_labels(self, Graph other, list ignore_labels=None) except -2:
        """
        Returns ``True`` if `self` and `other` have the same labels on
        their vertices, with the same number of vertices bearing each
        label (i.e. the multisets of vertex labels are equal). Vertices
        without a label (an empty or unset `label` attribute) are ignored.
        Returns ``False`` otherwise.

        If `ignore_labels` is given, vertices bearing any of those labels
        are excluded from the comparison entirely.
        """
        cdef dict labels1, labels2
        cdef set skip

        skip = set(ignore_labels) if ignore_labels else set()

        labels1 = {}
        for vertex in self.vertices:
            label = vertex.label
            if label and label not in skip:
                labels1[label] = labels1.get(label, 0) + 1

        labels2 = {}
        for vertex in other.vertices:
            label = vertex.label
            if label and label not in skip:
                labels2[label] = labels2.get(label, 0) + 1

        return labels1 == labels2

    cpdef list get_edges_in_cycle(self, list vertices, bint sort=False):
        """
        For a given list of atoms comprising a ring, return the set of bonds
        connecting them, in order around the ring.

        If `sort=True`, then sort the vertices to match their connectivity.
        Otherwise, assumes that they are already sorted, which is true for
        cycles returned by get_relevant_cycles or get_smallest_set_of_smallest_rings.
        """
        cdef list edges
        cdef int i, j

        if sort:
            self.sort_cyclic_vertices(vertices)

        edges = []
        for i, j in zip(range(len(vertices)), range(-1, len(vertices) - 1)):
            try:
                edges.append(self.get_edge(vertices[i], vertices[j]))
            except ValueError:
                raise ValueError('Edge does not exist between vertices in ring. '
                                 'Check that the vertices are properly ordered '
                                 'such that consecutive vertices are connected.')

        return edges

def get_relevant_cycles(graph):
    """
    Compute the relevant cycles of `graph`. A cycle never spans more than one biconnected
    component, so each component's cycles are found independently and concatenated.

    Returns a list of cycles, each a list of Vertex objects in ring-traversal order (consecutive
    entries, including the wraparound from the last entry back to the first, are bonded).
    """
    vertex_index = {vertex: i for i, vertex in enumerate(graph.vertices)}
    all_cycles = []
    for component_vertices, component_edges in _find_biconnected_components(graph, vertex_index):
        if len(component_edges) < len(component_vertices):
            # A tree (or a single bridge edge): contains no cycles at all.
            continue
        all_cycles.extend(_find_relevant_cycles_in_component(component_vertices, component_edges, vertex_index))
    return all_cycles

def _find_biconnected_components(graph, vertex_index):
    """
    Partition `graph`'s edges into biconnected components (maximal subgraphs where removing any
    single vertex leaves the rest connected) using Tarjan's algorithm, in its standard recursive
    form -- recursion depth here is bounded by the graph's DFS depth, i.e. O(V), a fundamentally
    different (and harmless) risk profile from the exponential branching this module's cycle
    detection used to do.

    The algorithm is a single depth-first traversal that tracks, for every vertex `v`:
      - `discovery_time[v]`: the order `v` was first reached in, 0, 1, 2, ...
      - `low_link[v]`: the *earliest* discovery_time reachable from anywhere in `v`'s DFS subtree
        by following at most one edge that jumps back to an already-visited ancestor (a "back
        edge"). This starts out equal to `discovery_time[v]` and only ever decreases, as `v`
        inherits the lowest low_link found among its children and among its own back edges.

    Every edge is pushed onto `edge_stack` as DFS descends: tree edges when moving to an unvisited
    neighbor, back edges when a neighbor turns out to already be an ancestor (the mirror case --
    a neighbor already visited as a *descendant* -- is skipped, since that's this same edge seen
    from its other endpoint). A back edge by itself only says the vertices between it and its
    ancestor lie on *some* cycle; it says nothing about where that component ends, since several
    back edges from different (possibly overlapping) subtrees can resolve into the same component
    or into separate ones, depending on how far back each one reaches.

    That's what `low_link` resolves. When DFS finishes exploring a child `w` of `v`, compare
    `low_link[w]` against `discovery_time[v]`:
      - If `low_link[w] >= discovery_time[v]`, nothing in `w`'s whole subtree ever found a back
        edge reaching further up than `v` -- so `v` is a cut vertex separating that subtree from
        the rest of the graph discovered so far, and everything on `edge_stack` back down to (and
        including) the edge `(v, w)` is popped off as one finished biconnected component right
        then (this also correctly closes out the last component when `v` is the root of its DFS
        tree, since `discovery_time[root] == 0` and `low_link` is never negative).
      - Otherwise, `w`'s subtree reached back past `v`, so this branch is still part of a larger,
        still-open component: `v` inherits `w`'s low_link, nothing is closed, and DFS moves on to
        `v`'s next child.

    `vertex_index` maps each vertex to its position in `graph.vertices`; it exists only to make
    neighbor traversal order -- and so the whole computation's output -- deterministic.

    Returns a list of (vertices, edges) tuples: `vertices` is a list of the component's Vertex
    objects (canonically ordered), `edges` is a set of frozenset({v1, v2}) pairs identifying its
    edges.
    """
    discovery_time = {}
    low_link = {}
    edge_stack = []
    components = []
    next_discovery_time = [0]

    def neighbors(v):
        return sorted(v.edges.keys(), key=lambda w: vertex_index[w])

    def close_component(boundary_edge):
        # Pop the edge stack down to (and including) `boundary_edge`: that's exactly the set of
        # edges discovered since this component was entered, so it's exactly one finished
        # biconnected component.
        comp_edges = set()
        comp_vertices = set()
        while edge_stack:
            popped = edge_stack.pop()
            comp_edges.add(popped)
            comp_vertices.update(popped)
            if popped == boundary_edge:
                break
        components.append((sorted(comp_vertices, key=lambda w: vertex_index[w]), comp_edges))

    def dfs(v, parent_edge):
        discovery_time[v] = low_link[v] = next_discovery_time[0]
        next_discovery_time[0] += 1
        for w in neighbors(v):
            edge = frozenset((v, w))
            if edge == parent_edge:
                continue
            if w not in discovery_time:
                edge_stack.append(edge)
                dfs(w, edge)
                if low_link[w] < low_link[v]:
                    low_link[v] = low_link[w]
                if low_link[w] >= discovery_time[v]:
                    # Nothing in w's subtree reaches back past v, so v is a cut vertex here: close
                    # off everything accumulated since entering this edge as one finished
                    # biconnected component (this also correctly closes out the component when v
                    # is the root, since discovery_time[root] == 0 and low_link[w] is always >= 0).
                    close_component(edge)
            elif discovery_time[w] < discovery_time[v]:
                # A back edge to an ancestor -- note the mirror case (w already visited, but as a
                # *descendant*, discovery_time[w] > discovery_time[v]) is deliberately not handled
                # here: that's this same edge encountered from its other endpoint, already pushed
                # once above.
                edge_stack.append(edge)
                if discovery_time[w] < low_link[v]:
                    low_link[v] = discovery_time[w]

    for start in graph.vertices:
        if start not in discovery_time:
            dfs(start, None)

    return components

def _find_relevant_cycles_in_component(component_vertices, component_edges, vertex_index):
    """
    Find the relevant cycles within a single biconnected component via Vismara's algorithm: root a
    BFS tree at every vertex of the component in turn. For each edge (y,z) whose endpoints are
    equidistant from the root (an *odd*-length cycle candidate), or each pair of vertices (y,z)
    sharing a common neighbor x with y,z equidistant from the root and one step closer than x (an
    *even*-length candidate), check whether the shortest paths from the root to y and to z share
    only the root -- if so, that pair yields a relevant cycle through the root.
    """
    neighbors_in_component = {v: [] for v in component_vertices}
    for edge in component_edges:
        v1, v2 = tuple(edge)
        neighbors_in_component[v1].append(v2)
        neighbors_in_component[v2].append(v1)
    for v in neighbors_in_component:
        neighbors_in_component[v].sort(key=lambda w: vertex_index[w])

    def edge_key(cycle):
        n = len(cycle)
        return frozenset(frozenset((cycle[i], cycle[i - 1])) for i in range(n))

    seen_edge_keys = set()
    cycles = []

    for root in component_vertices:
        dist = {root: 0}
        pred = {root: None}
        order = [root]
        i = 0
        while i < len(order):
            u = order[i]
            i += 1
            for w in neighbors_in_component[u]:
                if w not in dist:
                    dist[w] = dist[u] + 1
                    pred[w] = u
                    order.append(w)

        def path_to_root(v):
            path = []
            while v is not None:
                path.append(v)
                v = pred[v]
            return path  # v, ..., root

        def shares_only_root(y, z):
            return (set(path_to_root(y)) & set(path_to_root(z))) == {root}

        def record_candidate(y, z, middle):
            # root -> ... -> y [-> middle] -> z -> ... -> (back to root, implicitly, since this is
            # a ring) -- drop the trailing root from the z-side path since it's already the first
            # entry (from the reversed y-side path), or this cycle would list root twice
            cycle = list(reversed(path_to_root(y))) + middle + path_to_root(z)[:-1]
            key = edge_key(cycle)
            if key not in seen_edge_keys:
                seen_edge_keys.add(key)
                cycles.append(cycle)

        # Odd-length candidates: a direct edge between two vertices equidistant from the root.
        for edge in component_edges:
            y, z = tuple(edge)
            if y not in dist or z not in dist or dist[y] != dist[z] or dist[y] == 0:
                continue
            if shares_only_root(y, z):
                record_candidate(y, z, [])

        # Even-length candidates: two vertices sharing a common neighbor x, both one step closer
        # to the root than x.
        for x in component_vertices:
            if x not in dist:
                continue
            candidates = [w for w in neighbors_in_component[x] if w in dist and dist[w] == dist[x] - 1]
            for a in range(len(candidates)):
                for b in range(a + 1, len(candidates)):
                    y, z = candidates[a], candidates[b]
                    if shares_only_root(y, z):
                        record_candidate(y, z, [x])

    cycles.sort(key=lambda cycle: (len(cycle), tuple(sorted(vertex_index[v] for v in cycle))))
    return cycles
