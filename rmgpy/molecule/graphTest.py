#!/usr/bin/env python
# encoding: utf-8

import unittest

from rmgpy.molecule.graph import Edge, Graph, Vertex 

################################################################################

class TestGraph(unittest.TestCase):
    """
    Contains unit tests of the Vertex, Edge, and Graph classes. Most of the
    functionality of Vertex and Edge is only meaningful when part of a graph,
    so we test them all together instead of having separate unit test classes
    for each.
    """

    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        vertices = [Vertex() for i in range(6)]
        edges = [
            Edge(vertices[0], vertices[1]),
            Edge(vertices[1], vertices[2]),
            Edge(vertices[2], vertices[3]),
            Edge(vertices[3], vertices[4]),
            Edge(vertices[4], vertices[5]),
        ]
        
        self.graph = Graph()
        for vertex in vertices: self.graph.addVertex(vertex)
        for edge in edges: self.graph.addEdge(edge)
        
    def test_addVertex(self):
        """
        Test the Graph.addVertex() method.
        """
        vertex = Vertex()
        self.graph.addVertex(vertex)
        self.assertTrue(vertex in self.graph.vertices)
        self.assertTrue(vertex.edges == {})
        
    def test_addEdge(self):
        """
        Test the Graph.addEdge() method.
        """
        vertex1 = Vertex(); vertex2 = Vertex(); edge = Edge(vertex1, vertex2)
        try:
            self.graph.addEdge(edge)
            self.fail('Added edge between vertices not in graph to graph.')
        except ValueError:
            pass
        self.graph.addVertex(vertex1)
        self.graph.addVertex(vertex2)
        self.graph.addEdge(edge)
        self.assertTrue(vertex1 in self.graph.vertices)
        self.assertTrue(vertex1 in vertex2.edges)
        self.assertTrue(vertex2 in self.graph.vertices)
        self.assertTrue(vertex2 in vertex1.edges)
        self.assertTrue(vertex1.edges[vertex2] is edge)
        self.assertTrue(vertex2.edges[vertex1] is edge)
        
    def test_getEdge(self):
        """
        Test the Graph.getEdge() method.
        """
        vertex1 = self.graph.vertices[2]
        vertex2 = self.graph.vertices[4]
        try:
            edge = self.graph.getEdge(vertex1, vertex2)
            self.fail('Returned an edge between vertices that should not be connected in graph.')
        except ValueError:
            pass
        vertex1 = self.graph.vertices[2]
        vertex2 = self.graph.vertices[3]
        edge = self.graph.getEdge(vertex1, vertex2)
        self.assertNotEqual(edge, None)
        self.assertTrue(isinstance(edge, Edge))
        self.assertTrue(vertex1.edges[vertex2] is edge)
        self.assertTrue(vertex2.edges[vertex1] is edge)

    def test_getEdges(self):
        """
        Test the Graph.getEdges() method.
        """
        vertex1 = self.graph.vertices[2]
        edges = self.graph.getEdges(vertex1)
        self.assertTrue(isinstance(edges, dict))
        self.assertEqual(len(edges), 2)
        self.assertTrue(self.graph.vertices[1] in edges)
        self.assertTrue(self.graph.vertices[3] in edges)

    def test_hasVertex(self):
        """
        Test the Graph.hasVertex() method.
        """
        vertex = Vertex()
        self.assertFalse(self.graph.hasVertex(vertex))
        for v in self.graph.vertices:
            self.assertTrue(self.graph.hasVertex(v))

    def test_hasEdge(self):
        """
        Test the Graph.hasEdge() method.
        """
        vertex1 = self.graph.vertices[2]
        vertex2 = self.graph.vertices[4]
        self.assertFalse(self.graph.hasEdge(vertex1, vertex2))
        vertex1 = self.graph.vertices[2]
        vertex2 = self.graph.vertices[3]
        self.assertTrue(self.graph.hasEdge(vertex1, vertex2))

    def test_removeVertex(self):
        """
        Test the Graph.removeVertex() method.
        """
        vertex = self.graph.vertices[2]
        self.assertTrue(self.graph.hasVertex(vertex))
        self.graph.removeVertex(vertex)
        self.assertFalse(self.graph.hasVertex(vertex))
        for v in self.graph.vertices:
            self.assertFalse(vertex in v.edges)

    def test_removeEdge(self):
        """
        Test the Graph.removeEdge() method.
        """
        vertex1 = self.graph.vertices[2]
        vertex2 = self.graph.vertices[3]
        self.assertTrue(self.graph.hasEdge(vertex1, vertex2))
        edge = self.graph.getEdge(vertex1, vertex2)
        self.graph.removeEdge(edge)
        self.assertFalse(vertex1 in vertex2.edges)
        self.assertFalse(vertex2 in vertex1.edges)
    
    def test_copy(self):
        """
        Test the graph copy function to ensure a complete copy of the graph is
        made while preserving vertices and edges.
        """
        
        vertices = [Vertex() for i in range(6)]
        edges = [
            Edge(vertices[0], vertices[1]),
            Edge(vertices[1], vertices[2]),
            Edge(vertices[2], vertices[3]),
            Edge(vertices[3], vertices[4]),
            Edge(vertices[4], vertices[5]),
        ]

        graph = Graph()
        for vertex in vertices: graph.addVertex(vertex)
        for edge in edges: graph.addEdge(edge)
        
        graph2 = graph.copy()
        for vertex in graph.vertices:
            self.assertTrue(graph2.hasVertex(vertex))
        for v1 in graph.vertices:
            for v2 in v1.edges:
                self.assertTrue(graph2.hasEdge(v1, v2))
                self.assertTrue(graph2.hasEdge(v2, v1))
        self.assertTrue(graph2.isIsomorphic(graph))
        self.assertTrue(graph.isIsomorphic(graph2))

    def test_split(self):
        """
        Test the graph split function to ensure a proper splitting of the graph
        is being done.
        """

        vertices = [Vertex() for i in range(6)]
        edges = [
            Edge(vertices[0], vertices[1]),
            Edge(vertices[1], vertices[2]),
            Edge(vertices[2], vertices[3]),
            Edge(vertices[4], vertices[5]),
        ]

        graph = Graph()
        for vertex in vertices: graph.addVertex(vertex)
        for edge in edges: graph.addEdge(edge)
        
        graphs = graph.split()
        
        self.assertTrue(len(graphs) == 2)
        self.assertTrue(len(graphs[0].vertices) == 4 or len(graphs[0].vertices) == 2)
        self.assertTrue(len(graphs[0].vertices) + len(graphs[1].vertices) == len(graph.vertices))
    
    def test_merge(self):
        """
        Test the graph merge function to ensure a proper merging of the graph
        is being done.
        """

        vertices1 = [Vertex() for i in range(4)]
        edges1 = [
            Edge(vertices1[0], vertices1[1]),
            Edge(vertices1[1], vertices1[2]),
            Edge(vertices1[2], vertices1[3]),
        ]

        vertices2 = [Vertex() for i in range(3)]
        edges2 = [
            Edge(vertices2[0], vertices2[1]),
            Edge(vertices2[1], vertices2[2]),
        ]

        graph1 = Graph()
        for vertex in vertices1: graph1.addVertex(vertex)
        for edge in edges1: graph1.addEdge(edge)

        graph2 = Graph()
        for vertex in vertices2: graph2.addVertex(vertex)
        for edge in edges2: graph2.addEdge(edge)

        graph = graph1.merge(graph2)

        self.assertTrue(len(graph1.vertices) + len(graph2.vertices) == len(graph.vertices))
        for vertex1 in vertices1:
            self.assertTrue(vertex1 in graph.vertices)
            for vertex2 in vertex1.edges:
                self.assertTrue(vertex2 in graph.vertices)
        for vertex2 in vertices2:
            self.assertTrue(vertex2 in graph.vertices)
            for vertex1 in vertex2.edges:
                self.assertTrue(vertex1 in vertex2.edges)
    
    def test_resetConnectivityValues(self):
        """
        Test the Graph.resetConnectivityValues() method.
        """
        self.graph.resetConnectivityValues()
        for vertex in self.graph.vertices:
            self.assertEqual(vertex.connectivity1, -1)
            self.assertEqual(vertex.connectivity2, -1)
            self.assertEqual(vertex.connectivity3, -1)
            self.assertEqual(vertex.sortingLabel, -1)
    
    def test_updateConnectivityValues(self):
        """
        Test the Graph.updateConnectivityValues() method.
        """
        self.graph.updateConnectivityValues()
        self.assertEqual(self.graph.vertices[0].connectivity1, 1)
        self.assertEqual(self.graph.vertices[0].connectivity2, 2)
        self.assertEqual(self.graph.vertices[0].connectivity3, 3)
        self.assertEqual(self.graph.vertices[0].sortingLabel, -1)
        self.assertEqual(self.graph.vertices[1].connectivity1, 2)
        self.assertEqual(self.graph.vertices[1].connectivity2, 3)
        self.assertEqual(self.graph.vertices[1].connectivity3, 6)
        self.assertEqual(self.graph.vertices[1].sortingLabel, -1)
        self.assertEqual(self.graph.vertices[2].connectivity1, 2)
        self.assertEqual(self.graph.vertices[2].connectivity2, 4)
        self.assertEqual(self.graph.vertices[2].connectivity3, 7)
        self.assertEqual(self.graph.vertices[2].sortingLabel, -1)
        self.assertEqual(self.graph.vertices[3].connectivity1, 2)
        self.assertEqual(self.graph.vertices[3].connectivity2, 4)
        self.assertEqual(self.graph.vertices[3].connectivity3, 7)
        self.assertEqual(self.graph.vertices[3].sortingLabel, -1)
        self.assertEqual(self.graph.vertices[4].connectivity1, 2)
        self.assertEqual(self.graph.vertices[4].connectivity2, 3)
        self.assertEqual(self.graph.vertices[4].connectivity3, 6)
        self.assertEqual(self.graph.vertices[4].sortingLabel, -1)
        self.assertEqual(self.graph.vertices[5].connectivity1, 1)
        self.assertEqual(self.graph.vertices[5].connectivity2, 2)
        self.assertEqual(self.graph.vertices[5].connectivity3, 3)
        self.assertEqual(self.graph.vertices[5].sortingLabel, -1)
    
    def test_sortVertices(self):
        """
        Test the Graph.sortVertices() method.
        """
        self.graph.updateConnectivityValues()
        self.graph.sortVertices()
        for vertex1, vertex2 in zip(self.graph.vertices[:-1], self.graph.vertices[1:]):
            self.assertTrue(vertex1.sortingLabel < vertex2.sortingLabel)
            self.assertTrue(vertex1.connectivity3 >= vertex2.connectivity3)
            self.assertTrue(vertex1.connectivity2 >= vertex2.connectivity2)
            self.assertTrue(vertex1.connectivity1 >= vertex2.connectivity1)
    
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
        vertices = [Vertex() for i in range(6)]
        edges = [
            Edge(vertices[0], vertices[1]),
            Edge(vertices[1], vertices[2]),
            Edge(vertices[2], vertices[3]),
            Edge(vertices[3], vertices[4]),
            Edge(vertices[1], vertices[5]),
        ]
        
        graph = Graph()
        for vertex in vertices: graph.addVertex(vertex)
        for edge in edges: graph.addEdge(edge)
    
        graph.updateConnectivityValues()

        for i,cv_ in enumerate([1,3,2,2,1,1]):
            cv = vertices[i].connectivity1
            self.assertEqual(cv, cv_, "On vertex {0:d} got connectivity[0]={1:d} but expected {2:d}".format(i,cv,cv_))
        for i,cv_ in enumerate([3,4,5,3,2,3]):
            cv = vertices[i].connectivity2
            self.assertEqual(cv, cv_, "On vertex {0:d} got connectivity[0]={1:d} but expected {2:d}".format(i,cv,cv_))
        for i,cv_ in enumerate([4,11,7,7,3,4]):
            cv = vertices[i].connectivity3
            self.assertEqual(cv, cv_, "On vertex {0:d} got connectivity[0]={1:d} but expected {2:d}".format(i,cv,cv_))

    def test_isomorphism(self):
        """
        Check the graph isomorphism functions.
        """

        vertices1 = [Vertex() for i in range(6)]
        edges1 = [
            Edge(vertices1[0], vertices1[1]),
            Edge(vertices1[1], vertices1[2]),
            Edge(vertices1[2], vertices1[3]),
            Edge(vertices1[3], vertices1[4]),
            Edge(vertices1[4], vertices1[5]),
        ]

        vertices2 = [Vertex() for i in range(6)]
        edges2 = [
            Edge(vertices2[0], vertices2[1]),
            Edge(vertices2[1], vertices2[2]),
            Edge(vertices2[2], vertices2[3]),
            Edge(vertices2[3], vertices2[4]),
            Edge(vertices2[4], vertices2[5]),
        ]

        graph1 = Graph()
        for vertex in vertices1: graph1.addVertex(vertex)
        for edge in edges1: graph1.addEdge(edge)

        graph2 = Graph()
        for vertex in vertices2: graph2.addVertex(vertex)
        for edge in edges2: graph2.addEdge(edge)

        self.assertTrue(graph1.isIsomorphic(graph2))
        self.assertTrue(graph1.isSubgraphIsomorphic(graph2))
        self.assertTrue(graph2.isIsomorphic(graph1))
        self.assertTrue(graph2.isSubgraphIsomorphic(graph1))

    def test_subgraphIsomorphism(self):
        """
        Check the subgraph isomorphism functions.
        """

        vertices1 = [Vertex() for i in range(6)]
        edges1 = [
            Edge(vertices1[0], vertices1[1]),
            Edge(vertices1[1], vertices1[2]),
            Edge(vertices1[2], vertices1[3]),
            Edge(vertices1[3], vertices1[4]),
            Edge(vertices1[4], vertices1[5]),
        ]
        vertices2 = [Vertex() for i in range(2)]
        edges2 = [
            Edge(vertices2[0], vertices2[1]),
        ]

        graph1 = Graph()
        for vertex in vertices1: graph1.addVertex(vertex)
        for edge in edges1: graph1.addEdge(edge)

        graph2 = Graph()
        for vertex in vertices2: graph2.addVertex(vertex)
        for edge in edges2: graph2.addEdge(edge)

        self.assertFalse(graph1.isIsomorphic(graph2))
        self.assertFalse(graph2.isIsomorphic(graph1))
        self.assertTrue(graph1.isSubgraphIsomorphic(graph2))

        mapList = graph1.findSubgraphIsomorphisms(graph2)
        self.assertTrue(len(mapList) == 10)
        
        for mapping in mapList:
            self.assertTrue( graph1.isMappingValid(graph2,mapping) )
            self.assertTrue( graph1.isMappingValid(graph2,mapping) )
    
    def test_pickle(self):
        """
        Test that a Graph object can be successfully pickled and unpickled
        with no loss of information.
        """

        vertices = [Vertex() for i in range(6)]
        edges = [
            Edge(vertices[0], vertices[1]),
            Edge(vertices[1], vertices[2]),
            Edge(vertices[2], vertices[3]),
            Edge(vertices[3], vertices[4]),
            Edge(vertices[4], vertices[5]),
        ]

        graph0 = Graph()
        for vertex in vertices: graph0.addVertex(vertex)
        for edge in edges: graph0.addEdge(edge)
        graph0.updateConnectivityValues()
        
        import cPickle
        graph = cPickle.loads(cPickle.dumps(graph0))

        self.assertEqual(len(graph0.vertices), len(graph.vertices))
        for v1, v2 in zip(graph0.vertices, graph.vertices):
            self.assertEqual(v1.connectivity1, v2.connectivity1)
            self.assertEqual(v1.connectivity2, v2.connectivity2)
            self.assertEqual(v1.connectivity3, v2.connectivity3)
            self.assertEqual(v1.sortingLabel, v2.sortingLabel)
            self.assertEqual(len(v1.edges), len(v2.edges))
        self.assertTrue(graph0.isIsomorphic(graph))
        self.assertTrue(graph.isIsomorphic(graph0))

    def test_isCyclic(self):
        """
        Test the Graph.isCyclic() method.
        """
        self.assertFalse(self.graph.isCyclic())
        edge = Edge(self.graph.vertices[0], self.graph.vertices[3])
        self.graph.addEdge(edge) # To create a cycle
        self.assertTrue(self.graph.isCyclic())
        
    def test_isVertexInCycle(self):
        """
        Test the Graph.isVertexInCycle() method.
        """
        for vertex in self.graph.vertices:
            self.assertFalse(self.graph.isVertexInCycle(vertex))
        edge = Edge(self.graph.vertices[0], self.graph.vertices[3])
        self.graph.addEdge(edge) # To create a cycle
        for vertex in self.graph.vertices[0:4]:
            self.assertTrue(self.graph.isVertexInCycle(vertex))
        for vertex in self.graph.vertices[4:]:
            self.assertFalse(self.graph.isVertexInCycle(vertex))
        
    def test_isEdgeInCycle(self):
        """
        Test the Graph.isEdgeInCycle() method.
        """
        for vertex1 in self.graph.vertices:
            for vertex2, edge in vertex1.edges.items():
                self.assertFalse(self.graph.isEdgeInCycle(edge))
        edge = Edge(self.graph.vertices[0], self.graph.vertices[3])
        self.graph.addEdge(edge) # To create a cycle
        for vertex1 in self.graph.vertices:
            for vertex2, edge in vertex1.edges.items():
                if self.graph.vertices.index(vertex1) < 4 and self.graph.vertices.index(vertex2) < 4:
                    self.assertTrue(self.graph.isEdgeInCycle(edge))
                else:
                    self.assertFalse(self.graph.isEdgeInCycle(edge))
                    
    def test_getAllCyclicVertices(self):
        self.assertListEqual(self.graph.getAllCyclicVertices(), [])
        edge = Edge(self.graph.vertices[0], self.graph.vertices[3])
        self.graph.addEdge(edge) # To create a cycle
        self.assertEqual(len(self.graph.getAllCyclicVertices()), 4)
        
    def test_getAllPolycylicVertices(self):        
        edge = Edge(self.graph.vertices[0], self.graph.vertices[3])
        self.graph.addEdge(edge) # To create a cycle        
        self.assertListEqual(self.graph.getAllPolycyclicVertices(), []) 
        edge2 = Edge(self.graph.vertices[0], self.graph.vertices[5])
        self.graph.addEdge(edge2) # Create another cycle to generate two fused cycles 
        self.assertEqual(len(self.graph.getAllPolycyclicVertices()), 2)
        # Add new vertices and edges to generate a spirocyclic cycle      
        vertices = [Vertex() for i in range(2)]        
        for vertex in vertices: self.graph.addVertex(vertex)
        edges = [
                 Edge(self.graph.vertices[5], self.graph.vertices[6]),
                 Edge(self.graph.vertices[6], self.graph.vertices[7]),
                 Edge(self.graph.vertices[5], self.graph.vertices[7]),
                 ]
        for edge in edges: self.graph.addEdge(edge)        
        self.assertEqual(len(self.graph.getAllPolycyclicVertices()), 3)
                      
    def test_getAllCycles(self):
        """
        Test the Graph.getAllCycles() method.
        """
        cycleList = self.graph.getAllCycles(self.graph.vertices[0])
        self.assertEqual(len(cycleList), 0)
        edge = Edge(self.graph.vertices[0], self.graph.vertices[3])
        self.graph.addEdge(edge) # To create a cycle
        cycleList = self.graph.getAllCycles(self.graph.vertices[0])
        self.assertEqual(len(cycleList), 2)
        self.assertEqual(len(cycleList[0]), 4)
        self.assertEqual(len(cycleList[1]), 4)
        
    def test_getSmallestSetOfSmallestRings(self):
        """
        Test the Graph.getSmallestSetOfSmallestRings() method.
        """
        cycleList = self.graph.getSmallestSetOfSmallestRings()
        self.assertEqual(len(cycleList), 0)
        edge = Edge(self.graph.vertices[0], self.graph.vertices[3])
        self.graph.addEdge(edge) # To create a cycle
        cycleList = self.graph.getSmallestSetOfSmallestRings()
        self.assertEqual(len(cycleList), 1)
        self.assertEqual(len(cycleList[0]), 4)

################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
