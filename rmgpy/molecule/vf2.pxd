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

from rmgpy.molecule.graph cimport Vertex, Edge, Graph

cdef class VF2:

    cdef Graph graph1, graph2

    cpdef Graph graphA, graphB
    
    cdef dict initial_mapping
    cdef bint subgraph
    cdef bint find_all
    cdef bint strict
    
    cdef bint is_match
    cdef list mapping_list
    
    cpdef bint is_isomorphic(self, Graph graph1, Graph graph2, dict initial_mapping, bint save_order=?, bint strict=?) except -2
        
    cpdef list find_isomorphism(self, Graph graph1, Graph graph2, dict initial_mapping, bint save_order=?, bint strict=?)

    cpdef bint is_subgraph_isomorphic(self, Graph graph1, Graph graph2, dict initial_mapping, bint save_order=?) except -2

    cpdef list find_subgraph_isomorphisms(self, Graph graph1, Graph graph2, dict initial_mapping, bint save_order=?)
    
    cdef isomorphism(self, Graph graph1, Graph graph2, dict initial_mapping, bint subgraph, bint find_all, bint save_order=?, bint strict=?)

    cdef bint match(self, int call_depth) except -2
        
    cpdef bint feasible(self, Vertex vertex1, Vertex vertex2) except -2
    
    cdef add_to_mapping(self, Vertex vertex1, Vertex vertex2)
        
    cdef remove_from_mapping(self, Vertex vertex1, Vertex vertex2)
