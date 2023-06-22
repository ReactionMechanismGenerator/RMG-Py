###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2021 Prof. William H. Green (whgreen@mit.edu),           #
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
    cdef public short sorting_label
    cdef public bint terminal
    cdef public Vertex mapping
    cdef public bint ignore
    
    cpdef Vertex copy(self)

    cpdef bint equivalent(self, Vertex other, bint strict=?) except -2

    cpdef bint is_specific_case_of(self, Vertex other) except -2

    cpdef reset_connectivity_values(self)

cpdef short get_vertex_connectivity_value(Vertex vertex) except 1 # all values should be negative

cpdef short get_vertex_sorting_label(Vertex vertex) except -1 # all values should be nonnegative

################################################################################

cdef class Edge(object):

    cdef public Vertex vertex1, vertex2
    
    cpdef Edge copy(self)

    cpdef bint equivalent(Edge self, Edge other) except -2

    cpdef bint is_specific_case_of(self, Edge other) except -2

    cpdef Vertex get_other_vertex(self, Vertex vertex)

################################################################################

cdef Vertex _get_edge_vertex1(Edge edge)

cdef Vertex _get_edge_vertex2(Edge edge)

cdef class Graph(object):

    cdef public list vertices
    
    cdef public list ordered_vertices

    cpdef Vertex add_vertex(self, Vertex vertex)

    cpdef Edge add_edge(self, Edge edge)

    cpdef list get_all_edges(self)

    cpdef dict get_edges(self, Vertex vertex)

    cpdef Edge get_edge(self, Vertex vertex1, Vertex vertex2)

    cpdef bint has_vertex(self, Vertex vertex) except -2

    cpdef bint has_edge(self, Vertex vertex1, Vertex vertex2) except -2

    cpdef remove_vertex(self, Vertex vertex)

    cpdef remove_edge(self, Edge edge)

    cpdef update_connectivity_values(self)
    
    cpdef Graph copy(self, bint deep=?)

    cpdef dict copy_and_map(self)

    cpdef Graph merge(self, Graph other)

    cpdef list split(self)

    cpdef reset_connectivity_values(self)

    cpdef sort_vertices(self, bint save_order=?)
    
    cpdef restore_vertex_order(self)

    cpdef bint is_isomorphic(self, Graph other, dict initial_map=?, bint generate_initial_map=?, bint save_order=?, bint strict=?) except -2

    cpdef list find_isomorphism(self, Graph other, dict initial_map=?, bint save_order=?, bint strict=?)

    cpdef bint is_subgraph_isomorphic(self, Graph other, dict initial_map=?, bint save_order=?) except -2

    cpdef list find_subgraph_isomorphisms(self, Graph other, dict initial_map=?, bint save_order=?)

    cpdef bint is_cyclic(self) except -2

    cpdef bint is_vertex_in_cycle(self, Vertex vertex) except -2

    cpdef bint is_edge_in_cycle(self, Edge edge) except -2

    cpdef bint _is_chain_in_cycle(self, list chain) except -2

    cpdef list get_all_cyclic_vertices(self)
    
    cpdef list get_all_polycyclic_vertices(self)
    
    cpdef list get_polycycles(self)
    
    cpdef list get_monocycles(self)
    
    cpdef tuple get_disparate_cycles(self)

    cpdef tuple _merge_cycles(self, list cycle_sets)

    cpdef list get_all_cycles(self, Vertex starting_vertex)

    cpdef list get_all_cycles_of_size(self, int size)

    cpdef list get_all_simple_cycles_of_size(self, int size)

    cpdef list _explore_cycles_recursively(self, list chain, list cycles)

    cpdef list get_smallest_set_of_smallest_rings(self)

    cpdef list get_relevant_cycles(self)

    cpdef list sort_cyclic_vertices(self, list vertices)

    cpdef int get_max_cycle_overlap(self)
    
    cpdef list get_largest_ring(self, Vertex vertex)
    
    cpdef bint is_mapping_valid(self, Graph other, dict mapping, bint equivalent=?, bint strict=?) except -2

    cpdef list get_edges_in_cycle(self, list vertices, bint sort=?)
