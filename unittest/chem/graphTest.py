#!/usr/bin/python
# -*- coding: utf-8 -*-

import unittest

import sys
sys.path.append('.')

from rmgpy.chem.graph import *

################################################################################

class GraphCheck(unittest.TestCase):

    def testCopy(self):
        """
        Test the graph copy function to ensure a complete copy of the graph is
        made while preserving vertices and edges.
        """
        
        vertices = [Vertex() for i in range(6)]
        edges = [Edge() for i in range(5)]

        graph = Graph()
        for vertex in vertices: graph.addVertex(vertex)
        graph.addEdge(vertices[0], vertices[1], edges[0])
        graph.addEdge(vertices[1], vertices[2], edges[1])
        graph.addEdge(vertices[2], vertices[3], edges[2])
        graph.addEdge(vertices[3], vertices[4], edges[3])
        graph.addEdge(vertices[4], vertices[5], edges[4])
        
        graph2 = graph.copy()
        for vertex in graph.vertices:
            self.assertTrue(vertex in graph2.edges)
            self.assertTrue(graph2.hasVertex(vertex))
        for v1 in graph.vertices:
            for v2 in graph.edges[v1]:
                self.assertTrue(graph2.hasEdge(v1, v2))
                self.assertTrue(graph2.hasEdge(v2, v1))

    def testConnectivityValues(self):
        """
        Tests the Connectivity Values 
        as introduced by Morgan (1965)
        http://dx.doi.org/10.1021/c160017a018
        
        First CV1 is the number of neighbours
        CV2 is the sum of neighbouring CV1 values
        CV3 is the sum of neighbouring CV2 values
        
        Graph:     Expected (and tested) values:
        
        0-1-2-3-4            1-3-2-2-1   3-4-5-3-2    4-11-7-7-3
        |                    |           |             |
        5                    1           3             4
        
        """
        vertices = [Vertex() for i in range(6)]
        edges = [Edge() for i in range(5)]

        graph = Graph()
        for vertex in vertices: graph.addVertex(vertex)
        graph.addEdge(vertices[0], vertices[1], edges[0])
        graph.addEdge(vertices[1], vertices[2], edges[1])
        graph.addEdge(vertices[2], vertices[3], edges[2])
        graph.addEdge(vertices[3], vertices[4], edges[3])
        graph.addEdge(vertices[1], vertices[5], edges[4])
    
        graph.updateConnectivityValues()

        for i,cv_ in enumerate([1,3,2,2,1,1]):
            cv = vertices[i].connectivity1
            self.assertEqual(cv, cv_, "On vertex %d got connectivity[0]=%d but expected %d"%(i,cv,cv_))
        for i,cv_ in enumerate([3,4,5,3,2,3]):
            cv = vertices[i].connectivity2
            self.assertEqual(cv, cv_, "On vertex %d got connectivity[1]=%d but expected %d"%(i,cv,cv_))
        for i,cv_ in enumerate([4,11,7,7,3,4]):
            cv = vertices[i].connectivity3
            self.assertEqual(cv, cv_, "On vertex %d got connectivity[2]=%d but expected %d"%(i,cv,cv_))

    def testSplit(self):
        """
        Test the graph split function to ensure a proper splitting of the graph
        is being done.
        """

        vertices = [Vertex() for i in range(6)]
        edges = [Edge() for i in range(4)]

        graph = Graph()
        for vertex in vertices: graph.addVertex(vertex)
        graph.addEdge(vertices[0], vertices[1], edges[0])
        graph.addEdge(vertices[1], vertices[2], edges[1])
        graph.addEdge(vertices[2], vertices[3], edges[2])
        graph.addEdge(vertices[4], vertices[5], edges[3])

        graphs = graph.split()

        self.assertTrue(len(graphs) == 2)
        self.assertTrue(len(graphs[0].vertices) == 4 or len(graphs[0].vertices) == 2)
        self.assertTrue(len(graphs[0].vertices) + len(graphs[1].vertices) == len(graph.vertices))
    
    def testMerge(self):
        """
        Test the graph merge function to ensure a proper merging of the graph
        is being done.
        """

        vertices1 = [Vertex() for i in range(4)]
        edges1 = [Edge() for i in range(3)]

        vertices2 = [Vertex() for i in range(3)]
        edges2 = [Edge() for i in range(2)]

        graph1 = Graph()
        for vertex in vertices1: graph1.addVertex(vertex)
        graph1.addEdge(vertices1[0], vertices1[1], edges1[0])
        graph1.addEdge(vertices1[1], vertices1[2], edges1[1])
        graph1.addEdge(vertices1[2], vertices1[3], edges1[2])

        graph2 = Graph()
        for vertex in vertices2: graph2.addVertex(vertex)
        graph2.addEdge(vertices2[0], vertices2[1], edges2[0])
        graph2.addEdge(vertices2[1], vertices2[2], edges2[1])

        graph = graph1.merge(graph2)

        self.assertTrue(len(graph1.vertices) + len(graph2.vertices) == len(graph.vertices))
        
    def testIsomorphism(self):
        """
        Check the graph isomorphism functions.
        """

        vertices1 = [Vertex() for i in range(6)]
        edges1 = [Edge() for i in range(5)]
        vertices2 = [Vertex() for i in range(6)]
        edges2 = [Edge() for i in range(5)]

        graph1 = Graph()
        for vertex in vertices1: graph1.addVertex(vertex)
        graph1.edges[vertices1[0]] = {                          vertices1[1]: edges1[0] }
        graph1.edges[vertices1[1]] = { vertices1[0]: edges1[0], vertices1[2]: edges1[1] }
        graph1.edges[vertices1[2]] = { vertices1[1]: edges1[1], vertices1[3]: edges1[2] }
        graph1.edges[vertices1[3]] = { vertices1[2]: edges1[2], vertices1[4]: edges1[3] }
        graph1.edges[vertices1[4]] = { vertices1[3]: edges1[3], vertices1[5]: edges1[4] }
        graph1.edges[vertices1[5]] = { vertices1[4]: edges1[4] }

        graph2 = Graph()
        for vertex in vertices2: graph2.addVertex(vertex)
        graph2.edges[vertices2[0]] = {                          vertices2[1]: edges2[4] }
        graph2.edges[vertices2[1]] = { vertices2[0]: edges2[4], vertices2[2]: edges2[3] }
        graph2.edges[vertices2[2]] = { vertices2[1]: edges2[3], vertices2[3]: edges2[2] }
        graph2.edges[vertices2[3]] = { vertices2[2]: edges2[2], vertices2[4]: edges2[1] }
        graph2.edges[vertices2[4]] = { vertices2[3]: edges2[1], vertices2[5]: edges2[0] }
        graph2.edges[vertices2[5]] = { vertices2[4]: edges2[0] }

        self.assertTrue(graph1.isIsomorphic(graph2))
        self.assertTrue(graph1.isSubgraphIsomorphic(graph2))
        self.assertTrue(graph2.isIsomorphic(graph1))
        self.assertTrue(graph2.isSubgraphIsomorphic(graph1))

    def testSubgraphIsomorphism(self):
        """
        Check the subgraph isomorphism functions.
        """

        vertices1 = [Vertex() for i in range(6)]
        edges1 = [Edge() for i in range(5)]
        vertices2 = [Vertex() for i in range(2)]
        edges2 = [Edge() for i in range(1)]

        graph1 = Graph()
        for vertex in vertices1: graph1.addVertex(vertex)
        graph1.edges[vertices1[0]] = {                          vertices1[1]: edges1[0] }
        graph1.edges[vertices1[1]] = { vertices1[0]: edges1[0], vertices1[2]: edges1[1] }
        graph1.edges[vertices1[2]] = { vertices1[1]: edges1[1], vertices1[3]: edges1[2] }
        graph1.edges[vertices1[3]] = { vertices1[2]: edges1[2], vertices1[4]: edges1[3] }
        graph1.edges[vertices1[4]] = { vertices1[3]: edges1[3], vertices1[5]: edges1[4] }
        graph1.edges[vertices1[5]] = { vertices1[4]: edges1[4] }

        graph2 = Graph()
        for vertex in vertices2: graph2.addVertex(vertex)
        graph2.edges[vertices2[0]] = { vertices2[1]: edges2[0] }
        graph2.edges[vertices2[1]] = { vertices2[0]: edges2[0] }


        self.assertFalse(graph1.isIsomorphic(graph2))
        self.assertFalse(graph2.isIsomorphic(graph1))
        self.assertTrue(graph1.isSubgraphIsomorphic(graph2))

        ismatch, mapList = graph1.findSubgraphIsomorphisms(graph2)
        self.assertTrue(ismatch)
        self.assertTrue(len(mapList) == 10)

################################################################################

if __name__ == '__main__':
    unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )
