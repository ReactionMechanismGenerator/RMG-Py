###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2018 Prof. William H. Green (whgreen@mit.edu),           #
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

from graph cimport Vertex, Edge, Graph

cdef class VF2:

    cdef Graph graph1, graph2

    cpdef Graph graphA, graphB
    
    cdef dict initialMapping
    cdef bint subgraph
    cdef bint findAll
    
    cdef bint isMatch
    cdef list mappingList
    
    cpdef bint isIsomorphic(self, Graph graph1, Graph graph2, dict initialMapping, bint saveOrder=?) except -2
        
    cpdef list findIsomorphism(self, Graph graph1, Graph graph2, dict initialMapping, bint saveOrder=?)

    cpdef bint isSubgraphIsomorphic(self, Graph graph1, Graph graph2, dict initialMapping, bint saveOrder=?) except -2

    cpdef list findSubgraphIsomorphisms(self, Graph graph1, Graph graph2, dict initialMapping, bint saveOrder=?)
    
    cdef isomorphism(self, Graph graph1, Graph graph2, dict initialMapping, bint subgraph, bint findAll, bint saveOrder=?)

    cdef bint match(self, int callDepth) except -2
        
    cpdef bint feasible(self, Vertex vertex1, Vertex vertex2) except -2
    
    cdef addToMapping(self, Vertex vertex1, Vertex vertex2)
        
    cdef removeFromMapping(self, Vertex vertex1, Vertex vertex2)
