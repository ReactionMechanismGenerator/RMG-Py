#!/usr/bin/python
# -*- coding: utf-8 -*-

import unittest

import sys
sys.path.append('.')

from rmgpy.graph import *
from rmgpy.molecule import Atom, Bond

# Horrible hack because Cythonised graph class now demands vertices and edges are 
# instances of chem.Atom and chem.Bond 
Vertex = Atom
class Edge(Bond):
	def __init__(self):
		self.bondType = 'S'

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
		graph.addEdge((vertices[0], vertices[1]), edges[0])
		graph.addEdge((vertices[1], vertices[2]), edges[1])
		graph.addEdge((vertices[2], vertices[3]), edges[2])
		graph.addEdge((vertices[3], vertices[4]), edges[3])
		graph.addEdge((vertices[4], vertices[5]), edges[4])
		
		graph2 = graph.copy()
		for vertex in graph:
			self.assertTrue(vertex in graph2)
		for v1 in graph:
			for v2 in graph[v1]:
				self.assertTrue(graph2.hasEdge((v1, v2)))
				self.assertTrue(graph2.hasEdge((v2, v1)))

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
		#vertices = [Vertex() for i in range(6)]
		vertices = [Atom() for i in range(6)]
		edges = [Edge() for i in range(5)]

		graph = Graph()
		for vertex in vertices: graph.addVertex(vertex)
		graph.addEdge((vertices[0], vertices[1]), edges[0])
		graph.addEdge((vertices[1], vertices[2]), edges[1])
		graph.addEdge((vertices[2], vertices[3]), edges[2])
		graph.addEdge((vertices[3], vertices[4]), edges[3])
		graph.addEdge((vertices[1], vertices[5]), edges[4])
	
		graph.setConnectivityValues()

		for i,cv_ in enumerate([1,3,2,2,1,1]):
			cv = vertices[i].connectivity1
			#print "On vertex %d got connectivity[0] = %d and expected %d"%(i,cv,cv_)
			self.assertEqual(cv, cv_, "On vertex %d got connectivity[0]=%d but expected %d"%(i,cv,cv_))
		for i,cv_ in enumerate([3,4,5,3,2,3]):
			cv = vertices[i].connectivity2
			#print "On vertex %d got connectivity[1] = %d and expected %d"%(i,cv,cv_)
			self.assertEqual(cv, cv_, "On vertex %d got connectivity[1]=%d but expected %d"%(i,cv,cv_))
		for i,cv_ in enumerate([4,11,7,7,3,4]):
			cv = vertices[i].connectivity3
			#print "On vertex %d got connectivity[2] = %d and expected %d"%(i,cv,cv_)
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
		graph.addEdge((vertices[0], vertices[1]), edges[0])
		graph.addEdge((vertices[1], vertices[2]), edges[1])
		graph.addEdge((vertices[2], vertices[3]), edges[2])
		graph.addEdge((vertices[4], vertices[5]), edges[3])

		graphs = graph.split()

		self.assertTrue(len(graphs) == 2)
		self.assertTrue(len(graphs[0]) == 4 or len(graphs[0]) == 2)
		self.assertTrue(len(graphs[0]) + len(graphs[1]) == len(graph))
		self.assertTrue(len(graphs[0].edges()) == 3 or len(graphs[0].edges()) == 1)
		self.assertTrue(len(graphs[0].edges()) + len(graphs[1].edges()) == len(graph.edges()))

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
		graph1.addEdge((vertices1[0], vertices1[1]), edges1[0])
		graph1.addEdge((vertices1[1], vertices1[2]), edges1[1])
		graph1.addEdge((vertices1[2], vertices1[3]), edges1[2])

		graph2 = Graph()
		for vertex in vertices2: graph2.addVertex(vertex)
		graph2.addEdge((vertices2[0], vertices2[1]), edges2[0])
		graph2.addEdge((vertices2[1], vertices2[2]), edges2[1])

		graph = graph1.merge(graph2)

		self.assertTrue(len(graph1) + len(graph2) == len(graph))
		self.assertTrue(len(graph1.edges()) + len(graph2.edges()) == len(graph.edges()))

	def testIsomorphism(self):
		"""
		Check the graph isomorphism functions.
		"""

		vertices1 = [Vertex() for i in range(6)]
		edges1 = [Edge() for i in range(5)]
		vertices2 = [Vertex() for i in range(6)]
		edges2 = [Edge() for i in range(5)]

		graph1 = Graph()
		graph1[vertices1[0]] = {                          vertices1[1]: edges1[0] }
		graph1[vertices1[1]] = { vertices1[0]: edges1[0], vertices1[2]: edges1[1] }
		graph1[vertices1[2]] = { vertices1[1]: edges1[1], vertices1[3]: edges1[2] }
		graph1[vertices1[3]] = { vertices1[2]: edges1[2], vertices1[4]: edges1[3] }
		graph1[vertices1[4]] = { vertices1[3]: edges1[3], vertices1[5]: edges1[4] }
		graph1[vertices1[5]] = { vertices1[4]: edges1[4] }

		graph2 = Graph()
		graph2[vertices2[0]] = {                          vertices2[1]: edges2[4] }
		graph2[vertices2[1]] = { vertices2[0]: edges2[4], vertices2[2]: edges2[3] }
		graph2[vertices2[2]] = { vertices2[1]: edges2[3], vertices2[3]: edges2[2] }
		graph2[vertices2[3]] = { vertices2[2]: edges2[2], vertices2[4]: edges2[1] }
		graph2[vertices2[4]] = { vertices2[3]: edges2[1], vertices2[5]: edges2[0] }
		graph2[vertices2[5]] = { vertices2[4]: edges2[0] }

		self.assertTrue(graph1.isIsomorphic(graph2, {}, {}))
		self.assertTrue(graph1.isSubgraphIsomorphic(graph2, {}, {}))
		self.assertTrue(graph2.isIsomorphic(graph1, {}, {}))
		self.assertTrue(graph2.isSubgraphIsomorphic(graph1, {}, {}))

	def testSubgraphIsomorphism(self):
		"""
		Check the subgraph isomorphism functions.
		"""

		vertices1 = [Vertex() for i in range(6)]
		edges1 = [Edge() for i in range(5)]
		vertices2 = [Vertex() for i in range(2)]
		edges2 = [Edge() for i in range(1)]

		graph1 = Graph()
		graph1[vertices1[0]] = {                          vertices1[1]: edges1[0] }
		graph1[vertices1[1]] = { vertices1[0]: edges1[0], vertices1[2]: edges1[1] }
		graph1[vertices1[2]] = { vertices1[1]: edges1[1], vertices1[3]: edges1[2] }
		graph1[vertices1[3]] = { vertices1[2]: edges1[2], vertices1[4]: edges1[3] }
		graph1[vertices1[4]] = { vertices1[3]: edges1[3], vertices1[5]: edges1[4] }
		graph1[vertices1[5]] = { vertices1[4]: edges1[4] }

		graph2 = Graph()
		graph2[vertices2[0]] = { vertices2[1]: edges2[0] }
		graph2[vertices2[1]] = { vertices2[0]: edges2[0] }


		self.assertFalse(graph1.isIsomorphic(graph2, {}, {}))
		self.assertFalse(graph2.isIsomorphic(graph1, {}, {}))
		self.assertTrue(graph1.isSubgraphIsomorphic(graph2, {}, {}))

		ismatch, map21, map12 = graph1.findSubgraphIsomorphisms(graph2, {}, {})
		self.assertTrue(ismatch)
		self.assertTrue(len(map21) == len(map12) == 10)

	def testRingSearch1(self):
		"""
		Check the graph cycle identification functions.
		Tests a square with some substituents.
		
		 _|_
		| |
		 - --
		"""

		vertices = [Vertex() for i in range(8)]
		edges = [Edge() for i in range(8)]

		graph = Graph()
		for vertex in vertices: graph.addVertex(vertex)
		graph.addEdge((vertices[0], vertices[1]), edges[0])
		graph.addEdge((vertices[1], vertices[2]), edges[1])
		graph.addEdge((vertices[2], vertices[3]), edges[2])
		graph.addEdge((vertices[3], vertices[0]), edges[3])
		graph.addEdge((vertices[3], vertices[4]), edges[4])
		graph.addEdge((vertices[4], vertices[5]), edges[5])
		graph.addEdge((vertices[2], vertices[6]), edges[6])
		graph.addEdge((vertices[2], vertices[7]), edges[7])

		rings = graph.getSmallestSetOfSmallestRings()

		self.assertTrue(graph.isVertexInCycle(vertices[0]))
		self.assertTrue(graph.isVertexInCycle(vertices[1]))
		self.assertTrue(graph.isVertexInCycle(vertices[2]))
		self.assertTrue(graph.isVertexInCycle(vertices[3]))
		self.assertFalse(graph.isVertexInCycle(vertices[4]))
		self.assertFalse(graph.isVertexInCycle(vertices[5]))
		self.assertFalse(graph.isVertexInCycle(vertices[6]))
		self.assertFalse(graph.isVertexInCycle(vertices[7]))

#		self.assertTrue(graph.isEdgeInCycle(edges[0]))
#		self.assertTrue(graph.isEdgeInCycle(edges[1]))
#		self.assertTrue(graph.isEdgeInCycle(edges[2]))
#		self.assertTrue(graph.isEdgeInCycle(edges[3]))
#		self.assertFalse(graph.isEdgeInCycle(edges[4]))
#		self.assertFalse(graph.isEdgeInCycle(edges[5]))
#		self.assertFalse(graph.isEdgeInCycle(edges[6]))

		self.assertTrue(len(rings) == 1)
		self.assertTrue(len(rings[0]) == 4)
		
	def testRingSearch2(self):
		"""
		Check the graph cycle identification functions. This tests a fused ring
		system of two squares that share an edge.
		"""

		vertices = [Vertex() for i in range(6)]
		edges = [Edge() for i in range(7)]

		graph = Graph()
		for vertex in vertices: graph.addVertex(vertex)
		graph.addEdge((vertices[0], vertices[1]), edges[0])
		graph.addEdge((vertices[1], vertices[2]), edges[1])
		graph.addEdge((vertices[2], vertices[3]), edges[2])
		graph.addEdge((vertices[3], vertices[0]), edges[3])
		graph.addEdge((vertices[3], vertices[4]), edges[4])
		graph.addEdge((vertices[4], vertices[5]), edges[5])
		graph.addEdge((vertices[5], vertices[2]), edges[6])

		rings = graph.getSmallestSetOfSmallestRings()

		self.assertTrue(graph.isVertexInCycle(vertices[0]))
		self.assertTrue(graph.isVertexInCycle(vertices[1]))
		self.assertTrue(graph.isVertexInCycle(vertices[2]))
		self.assertTrue(graph.isVertexInCycle(vertices[3]))
		self.assertTrue(graph.isVertexInCycle(vertices[4]))
		self.assertTrue(graph.isVertexInCycle(vertices[5]))
		
#		self.assertTrue(graph.isEdgeInCycle(edges[0]))
#		self.assertTrue(graph.isEdgeInCycle(edges[1]))
#		self.assertTrue(graph.isEdgeInCycle(edges[2]))
#		self.assertTrue(graph.isEdgeInCycle(edges[3]))
#		self.assertTrue(graph.isEdgeInCycle(edges[4]))
#		self.assertTrue(graph.isEdgeInCycle(edges[5]))
#		self.assertTrue(graph.isEdgeInCycle(edges[6]))

		self.assertTrue(len(rings) == 2)
		self.assertTrue(len(rings[0]) == 4)
		self.assertTrue(len(rings[1]) == 4)

		self.assertFalse(rings[0][0] in rings[1] and rings[0][1] in rings[1] and rings[0][2] in rings[1] and rings[0][3] in rings[1])

	def testRingSearch3(self):
		"""
		Check the graph cycle identification functions. This tests a system
		of a triangle and a square connected by a two-vertex bridge.
		"""

		vertices = [Vertex() for i in range(9)]
		edges = [Edge() for i in range(10)]

		graph = Graph()
		for vertex in vertices: graph.addVertex(vertex)
		graph.addEdge((vertices[0], vertices[1]), edges[0])
		graph.addEdge((vertices[1], vertices[2]), edges[1])
		graph.addEdge((vertices[2], vertices[3]), edges[2])
		graph.addEdge((vertices[3], vertices[0]), edges[3])
		graph.addEdge((vertices[3], vertices[4]), edges[4])
		graph.addEdge((vertices[4], vertices[5]), edges[5])
		graph.addEdge((vertices[5], vertices[6]), edges[6])
		graph.addEdge((vertices[6], vertices[7]), edges[7])
		graph.addEdge((vertices[7], vertices[8]), edges[8])
		graph.addEdge((vertices[8], vertices[6]), edges[9])

		rings = graph.getSmallestSetOfSmallestRings()

		self.assertTrue(graph.isVertexInCycle(vertices[0]))
		self.assertTrue(graph.isVertexInCycle(vertices[1]))
		self.assertTrue(graph.isVertexInCycle(vertices[2]))
		self.assertTrue(graph.isVertexInCycle(vertices[3]))
		self.assertFalse(graph.isVertexInCycle(vertices[4]))
		self.assertFalse(graph.isVertexInCycle(vertices[5]))
		self.assertTrue(graph.isVertexInCycle(vertices[6]))
		self.assertTrue(graph.isVertexInCycle(vertices[7]))
		self.assertTrue(graph.isVertexInCycle(vertices[8]))

#		self.assertTrue(graph.isEdgeInCycle(edges[0]))
#		self.assertTrue(graph.isEdgeInCycle(edges[1]))
#		self.assertTrue(graph.isEdgeInCycle(edges[2]))
#		self.assertTrue(graph.isEdgeInCycle(edges[3]))
#		self.assertFalse(graph.isEdgeInCycle(edges[4]))
#		self.assertFalse(graph.isEdgeInCycle(edges[5]))
#		self.assertFalse(graph.isEdgeInCycle(edges[6]))
#		self.assertTrue(graph.isEdgeInCycle(edges[7]))
#		self.assertTrue(graph.isEdgeInCycle(edges[8]))
#		self.assertTrue(graph.isEdgeInCycle(edges[9]))

		self.assertTrue(len(rings) == 2)
		self.assertTrue(len(rings[0]) == 3 or len(rings[0]) == 4)
		self.assertTrue(len(rings[0]) + len(rings[1]) == 7)
		
	def testRingSearch4(self):
		"""
		Check the graph cycle identification functions. The test graph is a
		cube of vertices.
		"""

		vertices = [Vertex() for i in range(8)]
		edges = [Edge() for i in range(12)]

		graph = Graph()
		for vertex in vertices: graph.addVertex(vertex)
		graph.addEdge((vertices[0], vertices[1]), edges[0])
		graph.addEdge((vertices[1], vertices[2]), edges[1])
		graph.addEdge((vertices[2], vertices[3]), edges[2])
		graph.addEdge((vertices[3], vertices[0]), edges[3])
		graph.addEdge((vertices[4], vertices[5]), edges[4])
		graph.addEdge((vertices[5], vertices[6]), edges[5])
		graph.addEdge((vertices[6], vertices[7]), edges[6])
		graph.addEdge((vertices[7], vertices[4]), edges[7])
		graph.addEdge((vertices[4], vertices[0]), edges[8])
		graph.addEdge((vertices[5], vertices[1]), edges[9])
		graph.addEdge((vertices[6], vertices[2]), edges[10])
		graph.addEdge((vertices[7], vertices[3]), edges[11])

		rings = graph.getSmallestSetOfSmallestRings()

		self.assertTrue(graph.isVertexInCycle(vertices[0]))
		self.assertTrue(graph.isVertexInCycle(vertices[1]))
		self.assertTrue(graph.isVertexInCycle(vertices[2]))
		self.assertTrue(graph.isVertexInCycle(vertices[3]))
		self.assertTrue(graph.isVertexInCycle(vertices[4]))
		self.assertTrue(graph.isVertexInCycle(vertices[5]))
		self.assertTrue(graph.isVertexInCycle(vertices[6]))
		self.assertTrue(graph.isVertexInCycle(vertices[7]))

#		self.assertTrue(graph.isEdgeInCycle(edges[0]))
#		self.assertTrue(graph.isEdgeInCycle(edges[1]))
#		self.assertTrue(graph.isEdgeInCycle(edges[2]))
#		self.assertTrue(graph.isEdgeInCycle(edges[3]))
#		self.assertTrue(graph.isEdgeInCycle(edges[4]))
#		self.assertTrue(graph.isEdgeInCycle(edges[5]))
#		self.assertTrue(graph.isEdgeInCycle(edges[6]))
#		self.assertTrue(graph.isEdgeInCycle(edges[7]))
#		self.assertTrue(graph.isEdgeInCycle(edges[8]))
#		self.assertTrue(graph.isEdgeInCycle(edges[9]))
#		self.assertTrue(graph.isEdgeInCycle(edges[10]))
#		self.assertTrue(graph.isEdgeInCycle(edges[11]))

		# Not exactly sure how many SSSR rings ought to be found for this case
		self.assertTrue(len(rings) == 5 or len(rings) == 6)



################################################################################

if __name__ == '__main__':
	unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )

