###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2020 Prof. William H. Green (whgreen@mit.edu),           #
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

import py_rdl

from rmgpy.molecule.vf2 import VF2

################################################################################

class Vertex(object):
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
    `edges`             ``dict``        Dictionary of edges with keys being neighboring vertices
    `sorting_label`      ``int``         An integer label used to sort the vertices
    =================== =============== ========================================
    
    """

    def __init__(self):
        self.edges = {}
        self.reset_connectivity_values()
        self.ignore = False

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

    def copy(self):
        """
        Return a copy of the vertex. The default implementation assumes that no
        semantic information is associated with each vertex, and therefore
        simply returns a new :class:`Vertex` object.
        """
        new = Vertex()
        return new

    def equivalent(self, other, strict=True):
        """
        Return :data:`True` if two vertices `self` and `other` are semantically
        equivalent, or :data:`False` if not. You should reimplement this
        function in a derived class if your vertices have semantic information.
        """
        return True

    def is_specific_case_of(self, other):
        """
        Return ``True`` if `self` is semantically more specific than `other`,
        or ``False`` if not. You should reimplement this function in a derived
        class if your edges have semantic information.
        """
        return True

    def reset_connectivity_values(self):
        """
        Reset the cached structure information for this vertex.
        """
        self.connectivity1 = -1
        self.connectivity2 = -1
        self.connectivity3 = -1
        self.sorting_label = -1
        self.terminal = False
        self.mapping = None

def get_vertex_connectivity_value(vertex):
    """
    Return a value used to sort vertices prior to poposing candidate pairs in
    :meth:`__VF2_pairs`. The value returned is based on the vertex's
    connectivity values (and assumes that they are set properly).
    """
    return ( -256*vertex.connectivity1 - 16*vertex.connectivity2 - vertex.connectivity3 )

def get_vertex_sorting_label(vertex):
    """
    Return a value used to sort vertices prior to poposing candidate pairs in
    :meth:`__VF2_pairs`. The value returned is based on the vertex's
    connectivity values (and assumes that they are set properly).
    """
    return vertex.sorting_label

################################################################################

class Edge(object):
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

    def copy(self):
        """
        Return a copy of the edge. The default implementation assumes that no
        semantic information is associated with each edge, and therefore
        simply returns a new :class:`Edge` object. Note that the vertices are
        not copied in this implementation.
        """
        new = Edge(self.vertex1, self.vertex2)
        return new

    def equivalent(self, other):
        """
        Return ``True`` if two edges `self` and `other` are semantically
        equivalent, or ``False`` if not. You should reimplement this
        function in a derived class if your edges have semantic information.
        """
        return True

    def is_specific_case_of(self, other):
        """
        Return ``True`` if `self` is semantically more specific than `other`,
        or ``False`` if not. You should reimplement this function in a derived
        class if your edges have semantic information.
        """
        return True

    def get_other_vertex(self, vertex):
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

vf2 = VF2()

def _get_edge_vertex1(edge):
    return edge.vertex1

def _get_edge_vertex2(edge):
    return edge.vertex2

class Graph(object):
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

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (Graph, (self.vertices,))

    def add_vertex(self, vertex):
        """
        Add a `vertex` to the graph. The vertex is initialized with no edges.
        """
        self.vertices.append(vertex)
        vertex.edges = dict()
        return vertex

    def add_edge(self, edge):
        """
        Add an `edge` to the graph. The two vertices in the edge must already
        exist in the graph, or a :class:`ValueError` is raised.
        """
        if edge.vertex1 not in self.vertices or edge.vertex2 not in self.vertices:
            raise ValueError('Attempted to add edge between vertices not in the graph.')
        edge.vertex1.edges[edge.vertex2] = edge
        edge.vertex2.edges[edge.vertex1] = edge
        return edge

    def get_all_edges(self):
        """
        Returns a list of all edges in the graph.
        """

        edge_set = set()
        for vertex in self.vertices:
            for edge in vertex.edges.values():
                edge_set.add(edge)

        return list(edge_set)

    def get_edges(self, vertex):
        """
        Return a dictionary of the edges involving the specified `vertex`.
        """
        return vertex.edges

    def get_edge(self, vertex1, vertex2):
        """
        Returns the edge connecting vertices `vertex1` and `vertex2`.
        """
        try:
            return vertex1.edges[vertex2]
        except KeyError:
            raise ValueError('The specified vertices are not connected by an edge in this graph.')

    def has_vertex(self, vertex):
        """
        Returns ``True`` if `vertex` is a vertex in the graph, or ``False`` if
        not.
        """
        return vertex in self.vertices

    def has_edge(self, vertex1, vertex2):
        """
        Returns ``True`` if vertices `vertex1` and `vertex2` are connected
        by an edge, or ``False`` if not.
        """
        return vertex1 in self.vertices and vertex2 in vertex1.edges

    def remove_vertex(self, vertex):
        """
        Remove `vertex` and all edges associated with it from the graph. Does
        not remove vertices that no longer have any edges as a result of this
        removal.
        """
        for vertex2 in vertex.edges:
            del vertex2.edges[vertex]
        vertex.edges = dict()
        self.vertices.remove(vertex)

    def remove_edge(self, edge):
        """
        Remove the specified `edge` from the graph.
        Does not remove vertices that no longer have any edges as a result of
        this removal.
        """
        del edge.vertex1.edges[edge.vertex2]
        del edge.vertex2.edges[edge.vertex1]

    def copy(self, deep=False):
        """
        Create a copy of the current graph. If `deep` is ``True``, a deep copy
        is made: copies of the vertices and edges are used in the new graph.
        If `deep` is ``False`` or not specified, a shallow copy is made: the
        original vertices and edges are used in the new graph.
        """

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

    def copy_and_map(self):
        """
        Create a deep copy of the current graph, and return the dict
        'mapping'. Method was modified from Graph.copy() method
        """

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

    def merge(self, other):
        """
        Merge two graphs so as to store them in a single Graph object.
        """

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

    def split(self):
        """
        Convert a single Graph object containing two or more unconnected graphs
        into separate graphs.
        """

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

    def reset_connectivity_values(self):
        """
        Reset any cached connectivity information. Call this method when you
        have modified the graph.
        """
        for vertex in self.vertices: vertex.reset_connectivity_values()

    def update_connectivity_values(self):
        """
        Update the connectivity values for each vertex in the graph. These are
        used to accelerate the isomorphism checking.
        """

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

    def sort_vertices(self, save_order=False):
        """
        Sort the vertices in the graph. This can make certain operations, e.g.
        the isomorphism functions, much more efficient.
        """

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

    def restore_vertex_order(self):
        """
        reorder the vertices to what they were before sorting
        if you saved the order
        """
        if not self.ordered_vertices or len(self.vertices) != len(self.ordered_vertices):
            raise ValueError('Number of vertices has changed cannot restore original vertex order')
        else:
            self.vertices = self.ordered_vertices

    def is_isomorphic(self, other, initial_map=None, save_order=False, strict=True):
        """
        Returns :data:`True` if two graphs are isomorphic and :data:`False`
        otherwise. Uses the VF2 algorithm of Vento and Foggia.

        Args:
            initial_map (dict, optional): initial atom mapping to use
            save_order (bool, optional):  if ``True``, reset atom order after performing atom isomorphism
            strict (bool, optional):     if ``False``, perform isomorphism ignoring electrons
        """
        return vf2.is_isomorphic(self, other, initial_map, save_order=save_order, strict=strict)

    def find_isomorphism(self, other, initial_map=None, save_order=False, strict=True):
        """
        Returns :data:`True` if `other` is subgraph isomorphic and :data:`False`
        otherwise, and the matching mapping.
        Uses the VF2 algorithm of Vento and Foggia.

        Args:
            initial_map (dict, optional): initial atom mapping to use
            save_order (bool, optional):  if ``True``, reset atom order after performing atom isomorphism
            strict (bool, optional):     if ``False``, perform isomorphism ignoring electrons
        """
        return vf2.find_isomorphism(self, other, initial_map, save_order=save_order, strict=strict)

    def is_subgraph_isomorphic(self, other, initial_map=None, save_order=False):
        """
        Returns :data:`True` if `other` is subgraph isomorphic and :data:`False`
        otherwise. Uses the VF2 algorithm of Vento and Foggia.
        """
        return vf2.is_subgraph_isomorphic(self, other, initial_map, save_order=save_order)

    def find_subgraph_isomorphisms(self, other, initial_map=None, save_order=False):
        """
        Returns :data:`True` if `other` is subgraph isomorphic and :data:`False`
        otherwise. Also returns the lists all of valid mappings.

        Uses the VF2 algorithm of Vento and Foggia.
        """
        return vf2.find_subgraph_isomorphisms(self, other, initial_map, save_order=save_order)

    def is_cyclic(self):
        """
        Return ``True`` if one or more cycles are present in the graph or
        ``False`` otherwise.
        """
        for vertex in self.vertices:
            if self.is_vertex_in_cycle(vertex):
                return True
        return False

    def is_vertex_in_cycle(self, vertex):
        """
        Return ``True`` if the given `vertex` is contained in one or more
        cycles in the graph, or ``False`` if not.
        """
        return self._is_chain_in_cycle([vertex])

    def is_edge_in_cycle(self, edge):
        """
        Return :data:`True` if the edge between vertices `vertex1` and `vertex2`
        is in one or more cycles in the graph, or :data:`False` if not.
        """
        cycles = self.get_all_cycles(edge.vertex1)
        for cycle in cycles:
            if edge.vertex2 in cycle:
                return True
        return False

    def _is_chain_in_cycle(self, chain):
        """
        Return ``True`` if the given `chain` of vertices is contained in one
        or more cycles or ``False`` otherwise. This function recursively calls
        itself.
        """

        vertex1 = chain[-1]
        for vertex2 in vertex1.edges:
            if vertex2 is chain[0] and len(chain) > 2:
                return True
            elif vertex2 not in chain:
                # Make the chain a little longer and explore again
                chain.append(vertex2)
                if self._is_chain_in_cycle(chain):
                    # We found a cycle, so the return value must be True
                    return True
                else:
                    # We did not find a cycle down this path, so remove the vertex from the chain
                    chain.remove(vertex2)
        # If we reach this point then we did not find any cycles involving this chain
        return False

    def get_all_cyclic_vertices(self):
        """ 
        Returns all vertices belonging to one or more cycles.        
        """
        # Loop through all vertices and check whether they are cyclic
        cyclic_vertices = []
        for vertex in self.vertices:
            if self.is_vertex_in_cycle(vertex):
                cyclic_vertices.append(vertex)
        return cyclic_vertices

    def get_all_polycyclic_vertices(self):
        """
        Return all vertices belonging to two or more cycles, fused or spirocyclic.
        """
        sssr = self.get_smallest_set_of_smallest_rings()
        polycyclic_vertices = []
        if sssr:
            vertices = []
            for cycle in sssr:
                for vertex in cycle:
                    if vertex not in vertices:
                        vertices.append(vertex)
                    else:
                        if vertex not in polycyclic_vertices:
                            polycyclic_vertices.append(vertex)
        return polycyclic_vertices

    def get_polycycles(self):
        """
        Return a list of cycles that are polycyclic.
        In other words, merge the cycles which are fused or spirocyclic into 
        a single polycyclic cycle, and return only those cycles. 
        Cycles which are not polycyclic are not returned.
        """

        sssr = self.get_smallest_set_of_smallest_rings()
        if not sssr:
            return []

        polycyclic_vertices = self.get_all_polycyclic_vertices()

        if not polycyclic_vertices:
            # no polycyclic vertices detected
            return []
        else:
            # polycyclic vertices found, merge cycles together
            # that have common polycyclic vertices            
            continuous_cycles = []
            for vertex in polycyclic_vertices:
                # First check if it is in any existing continuous cycles
                for cycle in continuous_cycles:
                    if vertex in cycle:
                        polycyclic_cycle = cycle
                        break
                else:
                    # Otherwise create a new cycle
                    polycyclic_cycle = set()
                    continuous_cycles.append(polycyclic_cycle)

                for cycle in sssr:
                    if vertex in cycle:
                        polycyclic_cycle.update(cycle)

            # convert each set to a list
            continuous_cycles = [list(cycle) for cycle in continuous_cycles]
            return continuous_cycles

    def get_monocycles(self):
        """
        Return a list of cycles that are monocyclic.
        """

        sssr = self.get_smallest_set_of_smallest_rings()
        if not sssr:
            return []

        polycyclic_vertices = self.get_all_polycyclic_vertices()

        if not polycyclic_vertices:
            # No polycyclic_vertices detected, all the rings from get_smallest_set_of_smallest_rings
            # are monocyclic
            return sssr

        polycyclic_sssr = []
        for vertex in polycyclic_vertices:
            for cycle in sssr:
                if vertex in cycle:
                    if cycle not in polycyclic_sssr:
                        polycyclic_sssr.append(cycle)

        # remove the polycyclic cycles from the list of SSSR, leaving behind just the monocyclics
        monocyclic_cycles = sssr
        for cycle in polycyclic_sssr:
            monocyclic_cycles.remove(cycle)
        return monocyclic_cycles

    def get_disparate_cycles(self):
        """
        Get all disjoint monocyclic and polycyclic cycle clusters in the molecule.
        Takes the RC and recursively merges all cycles which share vertices.
        
        Returns: monocyclic_cycles, polycyclic_cycles
        """
        rc = self.get_relevant_cycles()

        if not rc:
            return [], []

        # Convert cycles to sets
        cycle_sets = [set(cycle_list) for cycle_list in rc]

        # Merge connected cycles
        monocyclic_cycles, polycyclic_cycles = self._merge_cycles(cycle_sets)

        # Convert cycles back to lists
        monocyclic_cycles = [list(cycle_set) for cycle_set in monocyclic_cycles]
        polycyclic_cycles = [list(cycle_set) for cycle_set in polycyclic_cycles]

        return monocyclic_cycles, polycyclic_cycles

    def _merge_cycles(self, cycle_sets):
        """
        Recursively merges cycles that share common atoms.
        
        Returns one list with unmerged cycles and one list with merged cycles.
        """
        unmerged_cycles = []
        merged_cycles = []

        # Loop through each cycle
        for cycle in cycle_sets:
            merged = False
            new = False

            # Check if it's attached to an existing merged cycle
            for m_cycle in merged_cycles:
                if not m_cycle.isdisjoint(cycle):
                    m_cycle.update(cycle)
                    merged = True
                    # It should only match one merged cycle, so we can break here
                    break
            else:
                # If it doesn't match any existing merged cycles, initiate a new one
                m_cycle = cycle.copy()
                new = True

            # Check if the new merged cycle is attached to any of the unmerged cycles
            matched = []
            for i, u_cycle in enumerate(unmerged_cycles):
                if not m_cycle.isdisjoint(u_cycle):
                    m_cycle.update(u_cycle)
                    matched.append(i)
                    merged = True
            # Remove matched cycles from list of unmerged cycles
            for i in reversed(matched):
                del unmerged_cycles[i]

            if merged and new:
                merged_cycles.append(m_cycle)
            elif not merged:
                unmerged_cycles.append(cycle)

        # If any rings were successfully merged, try to merge further
        if len(merged_cycles) > 1:
            u, m = self._merge_cycles(merged_cycles)
            merged_cycles = u + m

        return unmerged_cycles, merged_cycles

    def get_all_cycles(self, starting_vertex):
        """
        Given a starting vertex, returns a list of all the cycles containing
        that vertex.

        This function returns a duplicate of each cycle because [0,1,2,3]
        is counted as separate from [0,3,2,1]
        """
        return self._explore_cycles_recursively([starting_vertex], [])

    def get_all_cycles_of_size(self, size):
        """
        Return a list of the all non-duplicate rings with length 'size'. The
        algorithm implements was adapted from a description by Fan, Panaye,
        Doucet, and Barbu (doi: 10.1021/ci00015a002)

        B. T. Fan, A. Panaye, J. P. Doucet, and A. Barbu. "Ring Perception: A
        New Algorithm for Directly Finding the Smallest Set of Smallest Rings
        from a Connection Table." *J. Chem. Inf. Comput. Sci.* **33**,
        p. 657-662 (1993).
        """

        # Make a copy of the graph so we don't modify the original
        graph = self.copy(deep=True)
        vertices = graph.vertices[:]

        # Step 1: Remove all terminal vertices
        done = False
        while not done:
            vertices_to_remove = []
            for vertex in graph.vertices:
                if len(vertex.edges) == 1: vertices_to_remove.append(vertex)
            done = len(vertices_to_remove) == 0
            # Remove identified vertices from graph
            for vertex in vertices_to_remove:
                graph.remove_vertex(vertex)

        # Step 2: Remove all other vertices that are not part of cycles
        vertices_to_remove = []
        for vertex in graph.vertices:
            found = graph.is_vertex_in_cycle(vertex)
            if not found:
                vertices_to_remove.append(vertex)
        # Remove identified vertices from graph
        for vertex in vertices_to_remove:
            graph.remove_vertex(vertex)

        # Step 3: Split graph into remaining subgraphs
        graphs = graph.split()

        # Step 4: Find ring sets in each subgraph
        cycle_list = []

        for graph in graphs:

            while len(graph.vertices) > 0:

                # Choose root vertex as vertex with smallest number of edges
                root_vertex = None
                graph.update_connectivity_values()
                for vertex in graph.vertices:
                    if root_vertex is None:
                        root_vertex = vertex
                    elif get_vertex_connectivity_value(vertex) > get_vertex_connectivity_value(root_vertex):
                        root_vertex = vertex

                # Get all cycles involving the root vertex
                cycles = graph.get_all_cycles(root_vertex)
                if len(cycles) == 0:
                    # This vertex is no longer in a ring, so remove it
                    graph.remove_vertex(root_vertex)
                    continue

                # Keep the smallest of the cycles found above
                cycle = cycles[0]
                for c in cycles:
                    if len(c) == size: cycle_list.append(c)

                # Remove the root vertex to create single edges, note this will not
                # function properly if there is no vertex with 2 edges (i.e. cubane)
                graph.remove_vertex(root_vertex)

                # Remove from the graph all vertices in the cycle that have only one edge
                lone_carbon = True
                while lone_carbon:
                    lone_carbon = False
                    vertices_to_remove = []

                    for vertex in cycle:
                        if len(vertex.edges) == 1:
                            lone_carbon = True
                            vertices_to_remove.append(vertex)
                    else:
                        for vertex in vertices_to_remove:
                            graph.remove_vertex(vertex)

        # Map atoms in cycles back to atoms in original graph
        for i in range(len(cycle_list)):
            cycle_list[i] = [self.vertices[vertices.index(v)] for v in cycle_list[i]]

        #remove duplicates if there are more than 2 cycles:
        if len(cycle_list) < 2: return cycle_list
        cycle_set_list = [set(cycle_list[0])]
        for cycle1 in cycle_list[1:]:
            set1 = set(cycle1)
            for set2 in cycle_set_list:
                if set1 == set2:
                    break
            #not a duplicate so add it to cycle_set_list
            else:
                cycle_set_list.append(set1)

        #transform back to list of lists:
        cycle_set_list = [list(set1) for set1 in cycle_set_list]

        return cycle_set_list

    def get_all_simple_cycles_of_size(self, size):
        """
        Return a list of all non-duplicate monocyclic rings with length 'size'.

        Naive approach by eliminating polycyclic rings that are returned by
        ``getAllCyclicsOfSize``.
        """
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

    def _explore_cycles_recursively(self, chain, cycles):
        """
        Search the graph for cycles by recursive spidering. Given a `chain`
        (list) of connected atoms and a list of `cycles` found so far, find any
        cycles involving the chain of atoms and append them to the list of
        cycles. This function recursively calls itself.

        This function returns a duplicate of each cycle because [0,1,2,3]
        is counted as separate from [0,3,2,1]
        """

        vertex1 = chain[-1]
        # Loop over each of the atoms neighboring the last atom in the chain
        for vertex2 in vertex1.edges:
            if vertex2 is chain[0] and len(chain) > 2:
                # It is the first atom in the chain, so the chain is a cycle!
                cycles.append(chain[:])
            elif vertex2 not in chain:
                # Make the chain a little longer and explore again
                chain.append(vertex2)
                cycles = self._explore_cycles_recursively(chain, cycles)
                # Any cycles down this path have now been found, so remove vertex2 from the chain
                chain.pop(-1)
        # At this point we should have discovered all of the cycles involving the current chain
        return cycles

    def get_smallest_set_of_smallest_rings(self):
        """
        Returns the smallest set of smallest rings as a list of lists.
        Uses RingDecomposerLib for ring perception.

        Kolodzik, A.; Urbaczek, S.; Rarey, M.
        Unique Ring Families: A Chemically Meaningful Description
        of Molecular Ring Topologies.
        J. Chem. Inf. Model., 2012, 52 (8), pp 2013-2021

        Flachsenberg, F.; Andresen, N.; Rarey, M.
        RingDecomposerLib: An Open-Source Implementation of
        Unique Ring Families and Other Cycle Bases.
        J. Chem. Inf. Model., 2017, 57 (2), pp 122-126
        """

        graph = py_rdl.Graph.from_edges(
            self.get_all_edges(),
            _get_edge_vertex1,
            _get_edge_vertex2,
        )

        data = py_rdl.wrapper.DataInternal(graph.get_nof_nodes(), graph.get_edges().keys())
        data.calculate()

        sssr = []
        for cycle in data.get_sssr():
            sssr.append(self.sort_cyclic_vertices([graph.get_node_for_index(i) for i in cycle.nodes]))

        return sssr

    def get_relevant_cycles(self):
        """
        Returns the set of relevant cycles as a list of lists.
        Uses RingDecomposerLib for ring perception.

        Kolodzik, A.; Urbaczek, S.; Rarey, M.
        Unique Ring Families: A Chemically Meaningful Description
        of Molecular Ring Topologies.
        J. Chem. Inf. Model., 2012, 52 (8), pp 2013-2021

        Flachsenberg, F.; Andresen, N.; Rarey, M.
        RingDecomposerLib: An Open-Source Implementation of
        Unique Ring Families and Other Cycle Bases.
        J. Chem. Inf. Model., 2017, 57 (2), pp 122-126
        """

        graph = py_rdl.Graph.from_edges(
            self.get_all_edges(),
            _get_edge_vertex1,
            _get_edge_vertex2,
        )

        data = py_rdl.wrapper.DataInternal(graph.get_nof_nodes(), graph.get_edges().keys())
        data.calculate()

        rc = []
        for cycle in data.get_rcs():
            rc.append(self.sort_cyclic_vertices([graph.get_node_for_index(i) for i in cycle.nodes]))

        return rc

    def sort_cyclic_vertices(self, vertices):
        """
        Given a list of vertices comprising a cycle, sort them such that adjacent
        entries in the list are connected to each other.
        Warning: Assumes that the cycle is elementary, ie. no bridges.
        """

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

    def get_max_cycle_overlap(self):
        """
        Return the maximum number of vertices that are shared between
        any two cycles in the graph. For example, if there are only
        disparate monocycles or no cycles, the maximum overlap is zero;
        if there are "spiro" cycles, it is one; if there are "fused"
        cycles, it is two; and if there are "bridged" cycles, it is
        three.
        """

        cycles = self.get_smallest_set_of_smallest_rings()
        max_overlap = 0
        for i, j in itertools.combinations(range(len(cycles)), 2):
            overlap = len(set(cycles[i]) & set(cycles[j]))
            max_overlap = max(overlap, max_overlap)
        return max_overlap

    def get_largest_ring(self, vertex):
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

    def is_mapping_valid(self, other, mapping, equivalent=True, strict=True):
        """
        Check that a proposed `mapping` of vertices from `self` to `other`
        is valid by checking that the vertices and edges involved in the
        mapping are mutually equivalent.  If equivalent is ``True`` it checks
        if atoms and edges are equivalent, if ``False`` it checks if they
        are specific cases of each other. If strict is ``True``, electrons
        and bond orders are considered, and ignored if ``False``.
        """

        # Check that the mapped pairs of vertices compare True
        for vertex1, vertex2 in mapping.items():
            if equivalent:
                if not vertex1.equivalent(vertex2, strict=strict):
                    return False
            else:
                if not vertex1.is_specific_case_of(vertex2):
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

    def get_edges_in_cycle(self, vertices, sort=False):
        """
        For a given list of atoms comprising a ring, return the set of bonds
        connecting them, in order around the ring.

        If `sort=True`, then sort the vertices to match their connectivity.
        Otherwise, assumes that they are already sorted, which is true for
        cycles returned by get_relevant_cycles or get_smallest_set_of_smallest_rings.
        """

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
