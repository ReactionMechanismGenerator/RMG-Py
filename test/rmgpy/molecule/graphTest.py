#!/usr/bin/env python3

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


from rmgpy.molecule.graph import Edge, Graph, Vertex
import pytest


class TestGraph:
    """
    Contains unit tests of the Vertex, Edge, and Graph classes. Most of the
    functionality of Vertex and Edge is only meaningful when part of a graph,
    so we test them all together instead of having separate unit test classes
    for each.
    """

    def setup_class(self):
        """
        A function run before each unit test in this class.
        """
        vertices = [Vertex() for _ in range(6)]
        edges = [
            Edge(vertices[0], vertices[1]),
            Edge(vertices[1], vertices[2]),
            Edge(vertices[2], vertices[3]),
            Edge(vertices[3], vertices[4]),
            Edge(vertices[4], vertices[5]),
        ]

        self.graph = Graph(vertices)
        for edge in edges:
            self.graph.add_edge(edge)

    def test_vertices(self):
        """
        Test that the vertices attribute can be accessed.
        """
        vertices = self.graph.vertices
        assert isinstance(vertices, list)
        assert len(vertices) == 6

    def test_add_vertex(self):
        """
        Test the Graph.add_vertex() method.
        """
        vertex = Vertex()
        self.graph.add_vertex(vertex)
        assert vertex in self.graph.vertices
        assert vertex.edges == {}

    def test_add_edge(self):
        """
        Test the Graph.add_edge() method.
        """
        vertex1 = Vertex()
        vertex2 = Vertex()
        edge = Edge(vertex1, vertex2)
        try:
            self.graph.add_edge(edge)
            self.fail("Added edge between vertices not in graph to graph.")
        except ValueError:
            pass
        self.graph.add_vertex(vertex1)
        self.graph.add_vertex(vertex2)
        self.graph.add_edge(edge)
        assert vertex1 in self.graph.vertices
        assert vertex1 in vertex2.edges
        assert vertex2 in self.graph.vertices
        assert vertex2 in vertex1.edges
        assert vertex1.edges[vertex2] is edge
        assert vertex2.edges[vertex1] is edge

    def test_get_edge(self):
        """
        Test the Graph.get_edge() method.
        """
        vertex1 = self.graph.vertices[2]
        vertex2 = self.graph.vertices[4]
        try:
            self.graph.get_edge(vertex1, vertex2)
            self.fail("Returned an edge between vertices that should not be connected in graph.")
        except ValueError:
            pass
        vertex1 = self.graph.vertices[2]
        vertex2 = self.graph.vertices[3]
        edge = self.graph.get_edge(vertex1, vertex2)
        assert edge != None
        assert isinstance(edge, Edge)
        assert vertex1.edges[vertex2] is edge
        assert vertex2.edges[vertex1] is edge

    def test_get_edges(self):
        """
        Test the Graph.get_edges() method.
        """
        vertex1 = self.graph.vertices[2]
        edges = self.graph.get_edges(vertex1)
        assert isinstance(edges, dict)
        assert len(edges) == 2
        assert self.graph.vertices[1] in edges
        assert self.graph.vertices[3] in edges

    def test_get_all_edges(self):
        """
        Test the Graph.get_all_edges() method.
        """
        edges = self.graph.get_all_edges()
        assert isinstance(edges, list)
        assert len(edges) == 5

    def test_has_vertex(self):
        """
        Test the Graph.has_vertex() method.
        """
        vertex = Vertex()
        assert not self.graph.has_vertex(vertex)
        for v in self.graph.vertices:
            assert self.graph.has_vertex(v)

    def test_has_edge(self):
        """
        Test the Graph.has_edge() method.
        """
        vertex1 = self.graph.vertices[2]
        vertex2 = self.graph.vertices[4]
        assert not self.graph.has_edge(vertex1, vertex2)
        vertex1 = self.graph.vertices[2]
        vertex2 = self.graph.vertices[3]
        assert self.graph.has_edge(vertex1, vertex2)

    def test_remove_vertex(self):
        """
        Test the Graph.remove_vertex() method.
        """
        vertex = self.graph.vertices[2]
        assert self.graph.has_vertex(vertex)
        self.graph.remove_vertex(vertex)
        assert not self.graph.has_vertex(vertex)
        for v in self.graph.vertices:
            assert not (vertex in v.edges)

    def test_remove_edge(self):
        """
        Test the Graph.remove_edge() method.
        """
        vertex1 = self.graph.vertices[2]
        vertex2 = self.graph.vertices[3]
        assert self.graph.has_edge(vertex1, vertex2)
        edge = self.graph.get_edge(vertex1, vertex2)
        self.graph.remove_edge(edge)
        assert not (vertex1 in vertex2.edges)
        assert not (vertex2 in vertex1.edges)

    def test_copy(self):
        """
        Test the graph copy function to ensure a complete copy of the graph is
        made while preserving vertices and edges.
        """

        vertices = [Vertex() for _ in range(6)]
        edges = [
            Edge(vertices[0], vertices[1]),
            Edge(vertices[1], vertices[2]),
            Edge(vertices[2], vertices[3]),
            Edge(vertices[3], vertices[4]),
            Edge(vertices[4], vertices[5]),
        ]

        graph = Graph()
        for vertex in vertices:
            graph.add_vertex(vertex)
        for edge in edges:
            graph.add_edge(edge)

        graph2 = graph.copy()
        for vertex in graph.vertices:
            assert graph2.has_vertex(vertex)
        for v1 in graph.vertices:
            for v2 in v1.edges:
                assert graph2.has_edge(v1, v2)
                assert graph2.has_edge(v2, v1)
        assert graph2.is_isomorphic(graph)
        assert graph.is_isomorphic(graph2)

    def test_copy_and_map(self):
        """
        Test the returned dictionary points toward equivaalent vertices and edges
        """

        vertices = [Vertex() for _ in range(6)]
        edges = [
            Edge(vertices[0], vertices[1]),
            Edge(vertices[1], vertices[2]),
            Edge(vertices[2], vertices[3]),
            Edge(vertices[3], vertices[4]),
            Edge(vertices[4], vertices[5]),
        ]

        graph = Graph()
        for vertex in vertices:
            graph.add_vertex(vertex)
        for edge in edges:
            graph.add_edge(edge)

        graph_dict = graph.copy_and_map()
        graph2 = Graph(vertices=list(graph_dict.values()))

        for vertex in graph.vertices:
            assert graph2.has_vertex(graph_dict[vertex])
        for v1 in graph.vertices:
            for v2 in v1.edges:
                assert graph2.has_edge(graph_dict[v1], graph_dict[v2])
                assert graph2.has_edge(graph_dict[v2], graph_dict[v1])
        assert graph2.is_isomorphic(graph)
        assert graph.is_isomorphic(graph2)

    def test_split(self):
        """
        Test the graph split function to ensure a proper splitting of the graph
        is being done.
        """

        vertices = [Vertex() for _ in range(6)]
        edges = [
            Edge(vertices[0], vertices[1]),
            Edge(vertices[1], vertices[2]),
            Edge(vertices[2], vertices[3]),
            Edge(vertices[4], vertices[5]),
        ]

        graph = Graph()
        for vertex in vertices:
            graph.add_vertex(vertex)
        for edge in edges:
            graph.add_edge(edge)

        graphs = graph.split()

        assert len(graphs) == 2
        assert len(graphs[0].vertices) == 4 or len(graphs[0].vertices) == 2
        assert len(graphs[0].vertices) + len(graphs[1].vertices) == len(graph.vertices)

    def test_merge(self):
        """
        Test the graph merge function to ensure a proper merging of the graph
        is being done.
        """

        vertices1 = [Vertex() for _ in range(4)]
        edges1 = [
            Edge(vertices1[0], vertices1[1]),
            Edge(vertices1[1], vertices1[2]),
            Edge(vertices1[2], vertices1[3]),
        ]

        vertices2 = [Vertex() for _ in range(3)]
        edges2 = [
            Edge(vertices2[0], vertices2[1]),
            Edge(vertices2[1], vertices2[2]),
        ]

        graph1 = Graph()
        for vertex in vertices1:
            graph1.add_vertex(vertex)
        for edge in edges1:
            graph1.add_edge(edge)

        graph2 = Graph()
        for vertex in vertices2:
            graph2.add_vertex(vertex)
        for edge in edges2:
            graph2.add_edge(edge)

        graph = graph1.merge(graph2)

        assert len(graph1.vertices) + len(graph2.vertices) == len(graph.vertices)
        for vertex1 in vertices1:
            assert vertex1 in graph.vertices
            for vertex2 in vertex1.edges:
                assert vertex2 in graph.vertices
        for vertex2 in vertices2:
            assert vertex2 in graph.vertices
            for vertex1 in vertex2.edges:
                assert vertex1 in vertex2.edges

    def test_reset_connectivity_values(self):
        """
        Test the Graph.reset_connectivity_values() method.
        """
        self.graph.reset_connectivity_values()
        for vertex in self.graph.vertices:
            assert vertex.connectivity1 == -1
            assert vertex.connectivity2 == -1
            assert vertex.connectivity3 == -1
            assert vertex.sorting_label == -1

    def test_update_connectivity_values(self):
        """
        Test the Graph.update_connectivity_values() method.
        """
        self.graph.update_connectivity_values()
        assert self.graph.vertices[0].connectivity1 == 1
        assert self.graph.vertices[0].connectivity2 == 2
        assert self.graph.vertices[0].connectivity3 == 3
        assert self.graph.vertices[0].sorting_label == -1
        assert self.graph.vertices[1].connectivity1 == 2
        assert self.graph.vertices[1].connectivity2 == 3
        assert self.graph.vertices[1].connectivity3 == 6
        assert self.graph.vertices[1].sorting_label == -1
        assert self.graph.vertices[2].connectivity1 == 2
        assert self.graph.vertices[2].connectivity2 == 4
        assert self.graph.vertices[2].connectivity3 == 7
        assert self.graph.vertices[2].sorting_label == -1
        assert self.graph.vertices[3].connectivity1 == 2
        assert self.graph.vertices[3].connectivity2 == 4
        assert self.graph.vertices[3].connectivity3 == 7
        assert self.graph.vertices[3].sorting_label == -1
        assert self.graph.vertices[4].connectivity1 == 2
        assert self.graph.vertices[4].connectivity2 == 3
        assert self.graph.vertices[4].connectivity3 == 6
        assert self.graph.vertices[4].sorting_label == -1
        assert self.graph.vertices[5].connectivity1 == 1
        assert self.graph.vertices[5].connectivity2 == 2
        assert self.graph.vertices[5].connectivity3 == 3
        assert self.graph.vertices[5].sorting_label == -1

    def test_sort_vertices(self):
        """
        Test the Graph.sort_vertices() method.
        """
        self.graph.update_connectivity_values()
        self.graph.sort_vertices()
        for vertex1, vertex2 in zip(self.graph.vertices[:-1], self.graph.vertices[1:]):
            assert vertex1.sorting_label < vertex2.sorting_label
            assert vertex1.connectivity3 >= vertex2.connectivity3
            assert vertex1.connectivity2 >= vertex2.connectivity2
            assert vertex1.connectivity1 >= vertex2.connectivity1

    def test_vertex_connectivity_values(self):
        """
        Tests the vertex connectivity values as introduced by Morgan (1965).

        First CV1 is the number of neighbours
        CV2 is the sum of neighbouring CV1 values
        CV3 is the sum of neighbouring CV2 values

        Graph:     Expected (and tested) values:

        0-1-2-3-4            1-3-2-2-1   3-4-5-3-2    4-11-7-7-3
        |                    |           |             |
        5                    1           3             4

        """
        vertices = [Vertex() for _ in range(6)]
        edges = [
            Edge(vertices[0], vertices[1]),
            Edge(vertices[1], vertices[2]),
            Edge(vertices[2], vertices[3]),
            Edge(vertices[3], vertices[4]),
            Edge(vertices[1], vertices[5]),
        ]

        graph = Graph()
        for vertex in vertices:
            graph.add_vertex(vertex)
        for edge in edges:
            graph.add_edge(edge)

        graph.update_connectivity_values()

        for i, cv_ in enumerate([1, 3, 2, 2, 1, 1]):
            cv = vertices[i].connectivity1
            assert cv == cv_, "On vertex {0:d} got connectivity[0]={1:d} but expected {2:d}".format(i, cv, cv_)
        for i, cv_ in enumerate([3, 4, 5, 3, 2, 3]):
            cv = vertices[i].connectivity2
            assert cv == cv_, "On vertex {0:d} got connectivity[0]={1:d} but expected {2:d}".format(i, cv, cv_)
        for i, cv_ in enumerate([4, 11, 7, 7, 3, 4]):
            cv = vertices[i].connectivity3
            assert cv == cv_, "On vertex {0:d} got connectivity[0]={1:d} but expected {2:d}".format(i, cv, cv_)

    def test_isomorphism(self):
        """
        Check the graph isomorphism functions.
        """

        vertices1 = [Vertex() for _ in range(6)]
        edges1 = [
            Edge(vertices1[0], vertices1[1]),
            Edge(vertices1[1], vertices1[2]),
            Edge(vertices1[2], vertices1[3]),
            Edge(vertices1[3], vertices1[4]),
            Edge(vertices1[4], vertices1[5]),
        ]

        vertices2 = [Vertex() for _ in range(6)]
        edges2 = [
            Edge(vertices2[0], vertices2[1]),
            Edge(vertices2[1], vertices2[2]),
            Edge(vertices2[2], vertices2[3]),
            Edge(vertices2[3], vertices2[4]),
            Edge(vertices2[4], vertices2[5]),
        ]

        graph1 = Graph()
        for vertex in vertices1:
            graph1.add_vertex(vertex)
        for edge in edges1:
            graph1.add_edge(edge)

        graph2 = Graph()
        for vertex in vertices2:
            graph2.add_vertex(vertex)
        for edge in edges2:
            graph2.add_edge(edge)

        assert graph1.is_isomorphic(graph2)
        assert graph1.is_subgraph_isomorphic(graph2)
        assert graph2.is_isomorphic(graph1)
        assert graph2.is_subgraph_isomorphic(graph1)

    def test_isomorphism_disconnected(self):
        """
        Check the graph isomorphism for broken graphs.

        This tries to match graphs with a missing bond,
        eg. [ 0-1-2-3-4  5 ] should match [ 0-1-2-3-4  5 ]
        """

        vertices1 = [Vertex() for _ in range(6)]
        edges1 = [
            Edge(vertices1[0], vertices1[1]),
            Edge(vertices1[1], vertices1[2]),
            Edge(vertices1[2], vertices1[3]),
            Edge(vertices1[3], vertices1[4]),
            # Edge(vertices1[4], vertices1[5]),
        ]

        vertices2 = [Vertex() for _ in range(6)]
        edges2 = [
            Edge(vertices2[0], vertices2[1]),
            Edge(vertices2[1], vertices2[2]),
            Edge(vertices2[2], vertices2[3]),
            Edge(vertices2[3], vertices2[4]),
            # Edge(vertices2[4], vertices2[5]),
        ]

        graph1 = Graph()
        for vertex in vertices1:
            graph1.add_vertex(vertex)
        for edge in edges1:
            graph1.add_edge(edge)

        graph2 = Graph()
        for vertex in vertices2:
            graph2.add_vertex(vertex)
        for edge in edges2:
            graph2.add_edge(edge)

        assert graph1.is_isomorphic(graph2)
        assert graph1.is_subgraph_isomorphic(graph2)
        assert graph2.is_isomorphic(graph1)
        assert graph2.is_subgraph_isomorphic(graph1)
        assert len(graph1.find_subgraph_isomorphisms(graph2)) > 0

    def test_subgraph_isomorphism(self):
        """
        Check the subgraph isomorphism functions.
        """

        vertices1 = [Vertex() for _ in range(6)]
        edges1 = [
            Edge(vertices1[0], vertices1[1]),
            Edge(vertices1[1], vertices1[2]),
            Edge(vertices1[2], vertices1[3]),
            Edge(vertices1[3], vertices1[4]),
            Edge(vertices1[4], vertices1[5]),
        ]
        vertices2 = [Vertex() for _ in range(2)]
        edges2 = [
            Edge(vertices2[0], vertices2[1]),
        ]

        graph1 = Graph()
        for vertex in vertices1:
            graph1.add_vertex(vertex)
        for edge in edges1:
            graph1.add_edge(edge)

        graph2 = Graph()
        for vertex in vertices2:
            graph2.add_vertex(vertex)
        for edge in edges2:
            graph2.add_edge(edge)

        assert not graph1.is_isomorphic(graph2)
        assert not graph2.is_isomorphic(graph1)
        assert graph1.is_subgraph_isomorphic(graph2)

        map_list = graph1.find_subgraph_isomorphisms(graph2)
        assert len(map_list) == 10

        for mapping in map_list:
            assert graph1.is_mapping_valid(graph2, mapping)
            assert graph1.is_mapping_valid(graph2, mapping)

    def test_pickle(self):
        """
        Test that a Graph object can be successfully pickled and unpickled
        with no loss of information.
        """

        vertices = [Vertex() for _ in range(6)]
        edges = [
            Edge(vertices[0], vertices[1]),
            Edge(vertices[1], vertices[2]),
            Edge(vertices[2], vertices[3]),
            Edge(vertices[3], vertices[4]),
            Edge(vertices[4], vertices[5]),
        ]

        graph0 = Graph()
        for vertex in vertices:
            graph0.add_vertex(vertex)
        for edge in edges:
            graph0.add_edge(edge)
        graph0.update_connectivity_values()

        import pickle

        graph = pickle.loads(pickle.dumps(graph0))

        assert len(graph0.vertices) == len(graph.vertices)
        for v1, v2 in zip(graph0.vertices, graph.vertices):
            assert v1.connectivity1 == v2.connectivity1
            assert v1.connectivity2 == v2.connectivity2
            assert v1.connectivity3 == v2.connectivity3
            assert v1.sorting_label == v2.sorting_label
            assert len(v1.edges) == len(v2.edges)
        assert graph0.is_isomorphic(graph)
        assert graph.is_isomorphic(graph0)

    def test_is_cyclic(self):
        """
        Test the Graph.is_cyclic() method.
        """
        assert not self.graph.is_cyclic()
        edge = Edge(self.graph.vertices[0], self.graph.vertices[3])
        self.graph.add_edge(edge)  # To create a cycle
        assert self.graph.is_cyclic()

    def test_is_vertex_in_cycle(self):
        """
        Test the Graph.is_vertex_in_cycle() method.
        """
        for vertex in self.graph.vertices:
            assert not self.graph.is_vertex_in_cycle(vertex)
        edge = Edge(self.graph.vertices[0], self.graph.vertices[3])
        self.graph.add_edge(edge)  # To create a cycle
        for vertex in self.graph.vertices[0:4]:
            assert self.graph.is_vertex_in_cycle(vertex)
        for vertex in self.graph.vertices[4:]:
            assert not self.graph.is_vertex_in_cycle(vertex)

    def test_is_edge_in_cycle(self):
        """
        Test the Graph.is_edge_in_cycle() method.
        """
        for vertex1 in self.graph.vertices:
            for vertex2, edge in vertex1.edges.items():
                assert not self.graph.is_edge_in_cycle(edge)
        edge = Edge(self.graph.vertices[0], self.graph.vertices[3])
        self.graph.add_edge(edge)  # To create a cycle
        for vertex1 in self.graph.vertices:
            for vertex2, edge in vertex1.edges.items():
                if self.graph.vertices.index(vertex1) < 4 and self.graph.vertices.index(vertex2) < 4:
                    assert self.graph.is_edge_in_cycle(edge)
                else:
                    assert not self.graph.is_edge_in_cycle(edge)

    def test_get_all_cyclic_vertices(self):
        assert self.graph.get_all_cyclic_vertices() == []
        edge = Edge(self.graph.vertices[0], self.graph.vertices[3])
        self.graph.add_edge(edge)  # To create a cycle
        assert len(self.graph.get_all_cyclic_vertices()) == 4

    def test_get_all_polycylic_vertices(self):
        edge = Edge(self.graph.vertices[0], self.graph.vertices[3])
        self.graph.add_edge(edge)  # To create a cycle
        assert self.graph.get_all_polycyclic_vertices() == []
        edge2 = Edge(self.graph.vertices[0], self.graph.vertices[5])
        self.graph.add_edge(edge2)  # Create another cycle to generate two fused cycles
        assert len(self.graph.get_all_polycyclic_vertices()) == 2
        # Add new vertices and edges to generate a spirocyclic cycle
        vertices = [Vertex() for _ in range(2)]
        for vertex in vertices:
            self.graph.add_vertex(vertex)
        edges = [
            Edge(self.graph.vertices[5], self.graph.vertices[6]),
            Edge(self.graph.vertices[6], self.graph.vertices[7]),
            Edge(self.graph.vertices[5], self.graph.vertices[7]),
        ]
        for edge in edges:
            self.graph.add_edge(edge)
        assert len(self.graph.get_all_polycyclic_vertices()) == 3

    def test_get_all_cycles(self):
        """
        Test the Graph.get_all_cycles() method.
        """
        cycle_list = self.graph.get_all_cycles(self.graph.vertices[0])
        assert len(cycle_list) == 0
        edge = Edge(self.graph.vertices[0], self.graph.vertices[3])
        self.graph.add_edge(edge)  # To create a cycle
        cycle_list = self.graph.get_all_cycles(self.graph.vertices[0])
        assert len(cycle_list) == 2
        assert len(cycle_list[0]) == 4
        assert len(cycle_list[1]) == 4

    def test_get_all_cycles_of_size(self):
        """
        Test the Graph.getRingsOfSize() method
        """
        cycle_list = self.graph.get_all_cycles_of_size(6)
        assert len(cycle_list) == 0
        edge = Edge(self.graph.vertices[0], self.graph.vertices[3])
        self.graph.add_edge(edge)  # To create a cycle of length 4
        edge = Edge(self.graph.vertices[0], self.graph.vertices[5])
        self.graph.add_edge(edge)  # To create a cycle of length 6 and another cycle of length 4
        cycle_list = self.graph.get_all_cycles_of_size(4)
        assert len(cycle_list) == 2
        assert len(cycle_list[0]) == 4
        assert len(cycle_list[1]) == 4

    def test_get_all_simple_cycles_of_size(self):
        """
        Test the Graph.get_all_simple_cycles_of_size() method.
        """
        cycle_list = self.graph.get_all_cycles_of_size(6)
        assert len(cycle_list) == 0
        edge = Edge(self.graph.vertices[0], self.graph.vertices[3])
        self.graph.add_edge(edge)  # To create a cycle of length 4
        edge = Edge(self.graph.vertices[0], self.graph.vertices[5])
        self.graph.add_edge(edge)  # To create a cycle of length 6 and another cycle of length 4
        cycle_list = self.graph.get_all_simple_cycles_of_size(4)
        assert len(cycle_list) == 2
        assert len(cycle_list[0]) == 4
        assert len(cycle_list[1]) == 4
        cycle_list = self.graph.get_all_simple_cycles_of_size(6)
        assert len(cycle_list) == 0

    def test_get_smallest_set_of_smallest_rings(self):
        """
        Test the Graph.get_smallest_set_of_smallest_rings() method.
        """
        cycle_list = self.graph.get_smallest_set_of_smallest_rings()
        assert len(cycle_list) == 0
        edge = Edge(self.graph.vertices[0], self.graph.vertices[3])
        self.graph.add_edge(edge)  # To create a cycle
        cycle_list = self.graph.get_smallest_set_of_smallest_rings()
        assert len(cycle_list) == 1
        assert len(cycle_list[0]) == 4

    def test_get_relevant_cycles(self):
        """
        Test the Graph.get_relevant_cycles() method.
        """
        cycle_list = self.graph.get_relevant_cycles()
        assert len(cycle_list) == 0
        # Create a cycle of length 4
        edge = Edge(self.graph.vertices[0], self.graph.vertices[3])
        self.graph.add_edge(edge)
        # Create a second cycle of length 4
        edge = Edge(self.graph.vertices[0], self.graph.vertices[5])
        self.graph.add_edge(edge)
        # Create a bridge forming multiple cycles of length 4
        edge = Edge(self.graph.vertices[1], self.graph.vertices[4])
        self.graph.add_edge(edge)

        # SSSR should be 3 cycles of length 4
        cycle_list = self.graph.get_smallest_set_of_smallest_rings()
        assert len(cycle_list) == 3
        size_list = sorted([len(cycle) for cycle in cycle_list])
        assert size_list == [4, 4, 4]

        # RC should be 5 cycles of length 4
        cycle_list = self.graph.get_relevant_cycles()
        assert len(cycle_list) == 5
        size_list = sorted([len(cycle) for cycle in cycle_list])
        assert size_list == [4, 4, 4, 4, 4]

    def test_cycle_list_order_sssr(self):
        """
        Test that get_smallest_set_of_smallest_rings return vertices in the proper order.

        There are methods such as symmetry and molecule drawing which rely
        on the fact that subsequent list entries are connected.
        """
        # Create a cycle of length 5
        edge = Edge(self.graph.vertices[0], self.graph.vertices[4])
        self.graph.add_edge(edge)
        # Test SSSR
        sssr = self.graph.get_smallest_set_of_smallest_rings()
        assert len(sssr) == 1
        assert len(sssr[0]) == 5
        for i in range(5):
            assert self.graph.has_edge(sssr[0][i], sssr[0][i - 1])

    def test_cycle_list_order_relevant_cycles(self):
        """
        Test that get_relevant_cycles return vertices in the proper order.

        There are methods such as symmetry and molecule drawing which rely
        on the fact that subsequent list entries are connected.
        """
        # Create a cycle of length 5
        edge = Edge(self.graph.vertices[0], self.graph.vertices[4])
        self.graph.add_edge(edge)
        # Test RC
        rc = self.graph.get_relevant_cycles()
        assert len(rc) == 1
        assert len(rc[0]) == 5
        for i in range(5):
            assert self.graph.has_edge(rc[0][i], rc[0][i - 1])

    def test_get_polycyclic_rings(self):
        """
        Test that the Graph.get_polycycles() method returns only polycyclic rings.
        """
        vertices = [Vertex() for _ in range(27)]
        bonds = [
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 6),
            (6, 7),
            (7, 8),
            (8, 9),
            (9, 10),
            (10, 11),
            (11, 12),
            (12, 13),
            (13, 14),
            (14, 15),
            (14, 12),
            (12, 16),
            (16, 10),
            (10, 17),
            (17, 18),
            (18, 19),
            (9, 20),
            (20, 21),
            (21, 7),
            (6, 22),
            (22, 23),
            (22, 4),
            (23, 3),
            (23, 24),
            (24, 25),
            (25, 1),
        ]
        edges = []
        for bond in bonds:
            edges.append(Edge(vertices[bond[0]], vertices[bond[1]]))

        graph = Graph()
        for vertex in vertices:
            graph.add_vertex(vertex)
        for edge in edges:
            graph.add_edge(edge)
        graph.update_connectivity_values()

        sssr = graph.get_smallest_set_of_smallest_rings()
        assert len(sssr) == 6
        polycyclic_vertices = set(graph.get_all_polycyclic_vertices())
        expected_polycyclic_vertices = set([vertices[index] for index in [3, 23, 4, 22, 12]])

        assert polycyclic_vertices == expected_polycyclic_vertices

        continuous_rings = graph.get_polycycles()
        expected_continuous_rings = [
            [vertices[index] for index in [1, 2, 3, 4, 5, 6, 22, 23, 24, 25]],
            # [vertices[index] for index in [7,8,9,21,20]], # This is a nonpolycyclic ring
            [vertices[index] for index in [10, 11, 12, 13, 14, 16]],
        ]

        # Convert to sets for comparison purposes
        continuous_rings = [set(ring) for ring in continuous_rings]
        expected_continuous_rings = [set(ring) for ring in expected_continuous_rings]
        for ring in expected_continuous_rings:
            assert ring in continuous_rings

    def test_get_max_cycle_overlap(self):
        """
        Test that get_max_cycle_overlap returns the correct overlap numbers
        for different graphs.
        """

        def make_graph(edge_inds):
            nvert = max(max(inds) for inds in edge_inds) + 1
            vertices = [Vertex() for _ in range(nvert)]
            graph = Graph(vertices)
            for idx1, idx2 in edge_inds:
                graph.add_edge(Edge(vertices[idx1], vertices[idx2]))
            return graph

        linear = make_graph([(0, 1), (1, 2)])
        mono = make_graph([(0, 1), (0, 2), (1, 2), (2, 3), (3, 4), (3, 5), (4, 5)])
        spiro = make_graph([(0, 1), (0, 2), (1, 2), (2, 3), (2, 4), (3, 4)])
        fused = make_graph([(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)])
        bridged = make_graph([(0, 1), (0, 2), (1, 3), (1, 4), (2, 3), (2, 5), (4, 5)])
        cube = make_graph(
            [
                (0, 1),
                (0, 2),
                (0, 4),
                (1, 3),
                (1, 5),
                (2, 3),
                (2, 6),
                (3, 7),
                (4, 5),
                (4, 6),
                (5, 7),
                (6, 7),
            ]
        )

        assert linear.get_max_cycle_overlap() == 0
        assert mono.get_max_cycle_overlap() == 0
        assert spiro.get_max_cycle_overlap() == 1
        assert fused.get_max_cycle_overlap() == 2
        assert bridged.get_max_cycle_overlap() == 3
        # With the current algorithm for maximum overlap determination, a cube
        # only has an overlap of 2, because the set of relevant cycles
        # contains the six four-membered faces. This could be changed in the
        # future.
        assert cube.get_max_cycle_overlap() == 2

    def test_get_largest_ring(self):
        """
        Test that the Graph.get_polycycles() method returns only polycyclic rings.
        """
        vertices = [Vertex() for _ in range(27)]
        bonds = [
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 6),
            (6, 7),
            (9, 10),
            (10, 11),
            (11, 12),
            (12, 13),
            (13, 14),
            (14, 15),
            (12, 16),
            (10, 17),
            (17, 18),
            (18, 19),
            (9, 20),
            (20, 21),
            (6, 22),
            (22, 23),
            (22, 8),
            (8, 4),
            (23, 3),
            (23, 24),
            (24, 25),
            (25, 1),
        ]
        edges = []
        for bond in bonds:
            edges.append(Edge(vertices[bond[0]], vertices[bond[1]]))

        graph = Graph()
        for vertex in vertices:
            graph.add_vertex(vertex)
        for edge in edges:
            graph.add_edge(edge)
        graph.update_connectivity_values()

        rings = graph.get_polycycles()
        assert len(rings) == 1

        # ensure the last ring doesn't include vertex 8, since it isn't in the
        # longest ring. Try two different items since one might contain the vertex 8
        long_ring = graph.get_largest_ring(rings[0][0])
        long_ring2 = graph.get_largest_ring(rings[0][1])

        if len(long_ring) > len(long_ring2):
            longest_ring = long_ring
        else:
            longest_ring = long_ring2

        assert len(longest_ring) == len(rings[0]) - 1

    def test_sort_cyclic_vertices(self):
        """Test that sort_cyclic_vertices works properly for a valid input."""
        edge = Edge(self.graph.vertices[0], self.graph.vertices[5])
        self.graph.add_edge(edge)  # To create a cycle

        # Sort the vertices
        original = list(self.graph.vertices)
        ordered = self.graph.sort_cyclic_vertices(original)

        # Check that we didn't lose any vertices
        assert len(self.graph.vertices) == len(ordered), "Sorting changed the number of vertices."

        # Check that the order is different
        assert self.graph.vertices != ordered, "Sorting did not change the order of vertices."

        # Check that subsequent vertices are connected
        for i in range(5):
            assert self.graph.has_edge(ordered[i], ordered[i - 1])

    def test_sort_cyclic_vertices_invalid(self):
        """Test that sort_cyclic_vertices raises an error for an invalid input."""
        edge = Edge(self.graph.vertices[0], self.graph.vertices[4])
        self.graph.add_edge(edge)  # To create a cycle

        original = list(self.graph.vertices)

        with pytest.raises(RuntimeError, match="do not comprise a single cycle"):
            self.graph.sort_cyclic_vertices(original)

    def test_sort_cyclic_vertices_noncyclic(self):
        """Test that sort_cyclic_vertices raises an error for a noncyclic input."""
        original = list(self.graph.vertices)
        with pytest.raises(RuntimeError, match="do not comprise a single cycle"):
            self.graph.sort_cyclic_vertices(original)

    def test_sort_cyclic_vertices_unconnected(self):
        """Test that sort_cyclic_vertices raises an error for an unconnected input."""
        self.graph.add_vertex(Vertex())
        original = list(self.graph.vertices)
        with pytest.raises(RuntimeError, match="not all vertices are connected"):
            self.graph.sort_cyclic_vertices(original)
