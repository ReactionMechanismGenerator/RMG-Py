#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
#   ChemPy - A chemistry toolkit for Python
#
#   Copyright (c) 2010 by Joshua W. Allen (jwallen@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

"""
This module contains an implementation of a graph data structure (the 
:class:`Graph` class) and functions for manipulating that graph, including 
efficient isomorphism functions.
"""

import cython
import log as logging

################################################################################

class Vertex(object):
	"""
	A base class for vertices in a graph. Contains several attributes useful
	for accelerating isomorphism searches, as proposed by Morgan (1965); see
	http://dx.doi.org/10.1021/c160017a018 for more information.

	==================  ========================================================
	Attribute           Description
	==================  ========================================================
	`connectivity1`     The number of nearest neighbors
	`connectivity2`     The sum of the neighbors' `connectivity1` values
	`connectivity3`     The sum of the neighbors' `connectivity2` values
	`sorting_label`     An integer label useful for sorting vertices in a
	                    desired manner
	==================  ========================================================

	"""

	def __init__(self):
		self.resetCachedStructureInfo()
	
	def equivalent(self, other):
		"""
		Return :data:`True` if two vertices `self` and `other` are semantically
		equivalent, or :data:`False` if not. You should reimplement this
		function in a derived class if your vertices have semantic information.
		"""
		return True

	def resetCachedStructureInfo(self):
		"""
		Reset the cached structure information for this vertex.
		"""
		self.connectivity1 = -1
		self.connectivity2 = -1
		self.connectivity3 = -1
		self.sorting_label = -1

def __getSortLabel(vertex):
	"""
	Used to sort vertices prior to poposing candidate pairs in :method:`__VF2_pairs`

	This returns the `sorting_label` that is stored on the vertex. It should have been
	put there recently by a call to :method:`Graph.sortAndLabelVertices()`
	"""
	return vertex.sorting_label

def globalAtomSortValue(atom):
	"""
	Used to sort atoms prior to proposing candidate pairs in :method:`__VF2_pairs`
	The lowest (or most negative) values will be first in the list when you sort,
	so should be given to the atoms you want to explore first.
	For now, that is (roughly speaking) the most connected atoms. This definitely helps for large graphs
	but bizarrely the opposite ordering seems to help small graphs. Not sure about subggraphs...

	Assumes that atom.connictivityN (N=1..3) are all up to date.
	"""
	#return hash(atom)  # apparently random?
	#return (atom.connectivity[0] ) # not unique enough
	return ( -256*atom.connectivity1 - 16*atom.connectivity2 - atom.connectivity3 )

################################################################################

class Edge(object):
	"""
	A base class for edges in a graph. This class does *not* store the vertex
	pair that comprises the edge; that functionality would need to be included
	in the derived class.
	"""

	def __init__(self):
		pass

	def equivalent(self, other):
		"""
		Return :data:`True` if two edges `self` and `other` are semantically
		equivalent, or :data:`False` if not. You should reimplement this
		function in a derived class if your edges have semantic information.
		"""
		return True

################################################################################

class Graph(dict):
	"""
	A representation of a graph using a dictionary of dictionaries. The keys
	of the outer dictionary are the vertices, while edges are accessed via
	self[vertex1][vertex2].
	"""
	
	def __init__(self, vertices=None, edges=None):
		self.clear()
		if vertices is not None:
			for v in vertices: self.addVertex(v)
		if edges is not None:
			for e in edges: self.addEdge(e.atoms, e)
		self.resetCachedStructureInfo()

	def __reduce__(self):
		return (Graph, (), None, None, self.iteritems())
		
	def resetCachedStructureInfo(self):
		"""Reset any cached structural information.
		
		Call this method when you have modified the graph or structure,
		so that any information (eg. connectivity values, ring locations) that
		we are cacheing, is reset."""
		vert = cython.declare(Vertex)
		for vert in self: vert.resetCachedStructureInfo()
		
	def vertices(self):
		"""
		Return a list of the vertices in the graph.
		"""
		return self.keys()
	
	def edges(self):
		"""
		Return a list of the edges in the graph.
		"""
		edgelist = cython.declare(list)
		pairslist = cython.declare(list)
		
		edgelist = list()
		pairslist = list()
		for v1 in self:
			for v2 in self[v1]:
				if (v1, v2) not in pairslist:
					edgelist.append(self[v1][v2])
					pairslist.append((v1,v2))
					pairslist.append((v2,v1))
		return edgelist
	
	def addVertex(self, vertex):
		"""
		Add a `vertex` to the graph. The vertex is initialized with no edges.
		"""
		self[vertex] = dict()
		return vertex

	def addEdge(self, vertices, edge):
		"""
		Add an `edge` to the graph as an edge connecting the two vertices
		specified in the 2-tuple `vertices`.
		"""
		v1, v2 = vertices
		self[v1][v2] = edge
		self[v2][v1] = edge
		return edge

	def getEdges(self, vertex):
		"""
		Return a list of the edges involving the specified `vertex`.
		"""
		return self[vertex]

	def getEdge(self, vertices):
		"""
		Returns the edge connecting vertices in 2-tuple `vertices`, or None if
		no edge exists.
		"""
		v1, v2 = vertices
		return self[v1][v2] if self.hasEdge(vertices) else None

	def hasEdge(self, vertices):
		"""
		Returns :data:`True` if vertices in the 2-tuple `vertices` are connected
		by an edge, and :data:`False` if not.
		"""
		v1, v2 = vertices
		if v1 in self:
			return v2 in self[v1]
		return False

	def removeVertex(self, vertex1):
		"""
		Remove `vertex1` and all edges associated with it from the graph. Does
		not remove vertices that no longer have any edges as a result of this
		removal.
		"""
		for vertex2 in self:
			if vertex2 is not vertex1:
				if vertex1 in self[vertex2]:
					del self[vertex2][vertex1]
		del self[vertex1]

	def removeEdge(self, vertices):
		"""
		Remove the edge having vertices as specified in the 2-tuple `vertices`
		from the graph. Does not remove vertices that no longer have any edges
		as a result of this removal.
		"""
		v1, v2 = vertices
		del self[v1][v2]
		del self[v2][v1]

	def __isomorphism(self, other, subgraph, findAll, map12_0=None, map21_0=None):
		if map12_0 == None: map12_0 = dict()
		if map21_0 == None: map21_0 = dict()
		return VF2_isomorphism(self, other, map21_0, map12_0,
			subgraph, findAll)
		
	def isIsomorphic(self, other, map12_0=None, map21_0=None):
		"""
		Returns :data:`True` if two graphs are isomorphic and :data:`False`
		otherwise. Uses the VF2 algorithm of Vento and Foggia.
		"""
		if len(self) != len(other): return False
		ismatch, map21, map12 = self.__isomorphism(other, False, False, map12_0, map21_0)
		return ismatch

	def findIsomorphism(self, other, map12_0=None, map21_0=None):
		"""
		Returns :data:`True` if `other` is subgraph isomorphic and :data:`False`
		otherwise, and the matching mapping.
		Uses the VF2 algorithm of Vento and Foggia.
		"""
		return self.__isomorphism(other, False, False, map12_0, map21_0)

	def isSubgraphIsomorphic(self, other, map12_0=None, map21_0=None):
		"""
		Returns :data:`True` if `other` is subgraph isomorphic and :data:`False`
		otherwise. Uses the VF2 algorithm of Vento and Foggia.
		"""
		ismatch, map21, map12 = self.__isomorphism(other, True, False, map12_0, map21_0)
		return ismatch

	def findSubgraphIsomorphisms(self, other, map12_0=None, map21_0=None):
		"""
		Returns :data:`True` if `other` is subgraph isomorphic and :data:`False`
		otherwise. Also returns the lists all of valid mappings.
		
		Uses the VF2 algorithm of Vento and Foggia.
		"""
		return self.__isomorphism(other, True, True, map12_0, map21_0)

	def copy(self):
		"""
		Create a copy of the current graph. The vertices used in both graphs
		are the SAME vertices, not copies of each other, so it should not be used
		for copying structures (instead do structure.copy()) 
		"""
		other = cython.declare(Graph)
		other = Graph()
		for vertex in self:
			other.addVertex(vertex)
		for v1 in self:
			for v2 in self[v1]:
				other[v1][v2] = self[v1][v2]
		
		return other

	def merge(self, other):
		"""
		Merge two graphs so as to store them in a single Graph object.
		"""

		# Create output graph
		new = cython.declare(Graph)
		new = Graph()
		
		# Add vertices to output graph
		for vertex in self:
			new.addVertex(vertex)
		for vertex in other:
			new.addVertex(vertex)

		# Add edges to output graph
		for v1 in self:
			for v2 in self[v1]:
				new[v1][v2] = self[v1][v2]
		for v1 in other:
			for v2 in other[v1]:
				new[v1][v2] = other[v1][v2]

		return new

	def split(self):
		"""
		Convert a single Graph object containing two or more unconnected graphs
		into separate graphs.
		"""
		
		# Create potential output graphs
		new1 = cython.declare(Graph)
		new2 = cython.declare(Graph)
		verticesToMove = cython.declare(list)
		index = cython.declare(cython.int)
		
		new1 = self.copy()
		new2 = Graph()
		
		if len(self.vertices()) == 0:
			return [new1]

		# Arbitrarily choose last atom as starting point
		verticesToMove = [ self.vertices()[-1] ]

		# Iterate until there are no more atoms to move
		index = 0
		while index < len(verticesToMove):
			for v2 in self.getEdges(verticesToMove[index]):
				if v2 not in verticesToMove:
					verticesToMove.append(v2)
			index += 1

		# If all atoms are to be moved, simply return new1
		if len(new1.vertices()) == len(verticesToMove):
			return [new1]

		# Copy to new graph
		for vertex in verticesToMove:
			new2.addVertex(vertex)
		for v1 in verticesToMove:
			for v2, edge in new1[v1].iteritems():
				new2[v1][v2] = edge

		# Remove from old graph
		for v1 in new2:
			for v2 in new2[v1]:
				if v1 in verticesToMove and v2 in verticesToMove:
					del new1[v1][v2]
		for vertex in verticesToMove:
			new1.removeVertex(vertex)

		new = [new2]
		new.extend(new1.split())
		return new

	def getSmallestSetOfSmallestRings(self):
		"""
		Return a list of the smallest set of smallest rings in the graph. The
		algorithm implements was adapted from a description by Fan, Panaye,
		Doucet, and Barbu (doi: 10.1021/ci00015a002)

		B. T. Fan, A. Panaye, J. P. Doucet, and A. Barbu. "Ring Perception: A
		New Algorithm for Directly Finding the Smallest Set of Smallest Rings
		from a Connection Table." *J. Chem. Inf. Comput. Sci.* **33**, 
		p. 657-662 (1993).
		"""

		graph = cython.declare(Graph)
		done = cython.declare(cython.bint)
		verticesToRemove = cython.declare(list)
		cycleList = cython.declare(list)
		cycles = cython.declare(list)
		vertex = cython.declare(Vertex)
		rootVertex = cython.declare(Vertex)
		found = cython.declare(cython.bint)
		cycle = cython.declare(list)
		graphs = cython.declare(list)

		# Make a copy of the graph so we don't modify the original
		graph = self.copy()
		
		# Step 1: Remove all terminal vertices
		done = False
		while not done:
			verticesToRemove = []
			for vertex1, value in graph.iteritems():
				if len(value) == 1: verticesToRemove.append(vertex1)
			done = len(verticesToRemove) == 0
			# Remove identified vertices from graph
			for vertex in verticesToRemove:
				graph.removeVertex(vertex)
		
		# Step 2: Remove all other vertices that are not part of cycles
		verticesToRemove = []
		for vertex in graph:
			found = graph.isVertexInCycle(vertex)
			if not found:
				verticesToRemove.append(vertex)
		# Remove identified vertices from graph
		for vertex in verticesToRemove:
			graph.removeVertex(vertex)
			
		### also need to remove EDGES that are not in ring
		
		# Step 3: Split graph into remaining subgraphs
		graphs = graph.split()
		
		# Step 4: Find ring sets in each subgraph
		cycleList = []
		for graph in graphs:
			
			while len(graph) > 0:
				
				# Choose root vertex as vertex with smallest number of edges
				rootVertex = None
				for vertex in graph:
					if rootVertex is None:
						rootVertex = vertex
					elif len(graph[vertex]) < len(graph[rootVertex]):
						rootVertex = vertex

				# Get all cycles involving the root vertex
				cycles = graph.getAllCycles(rootVertex)
				if len(cycles) == 0:
					# this vertex is no longer in a ring.
					# remove all its edges
					neighbours = graph[rootVertex].keys()[:]
					for vertex2 in neighbours:
						graph.removeEdge((rootVertex, vertex2))
					# then remove it
					graph.removeVertex(rootVertex)
					#print("Removed vertex that's no longer in ring")
					continue # (pick a new root Vertex)
#					raise Exception('Did not find expected cycle!')

				# Keep the smallest of the cycles found above
				cycle = cycles[0]
				for c in cycles[1:]:
					if len(c) < len(cycle):
						cycle = c
				cycleList.append(cycle)

				# Remove from the graph all vertices in the cycle that have only two edges
				verticesToRemove = []
				for vertex in cycle:
					if len(graph[vertex]) <= 2:
						verticesToRemove.append(vertex)
				if len(verticesToRemove) == 0:
					# there are no vertices in this cycle that with only two edges
					
					# Remove edge between root vertex and any one vertex it is connected to
					graph.removeEdge((rootVertex, graph[rootVertex].keys()[0]))
				else:
					for vertex in verticesToRemove:
						graph.removeVertex(vertex)
						
		return cycleList

	def isVertexInCycle(self, vertex):
		""" 
		Is `vertex` in a cycle?
		Returns :data:`True` if it is in a cycle, else :data:`False`.
		"""
		chain = cython.declare(list)
		chain = [vertex]
		return self.__isChainInCycle(chain)

	def __isChainInCycle(self, chain):
		""" 
		Is the `chain` in a cycle? 
		Returns True/False.
		Recursively calls itself
		"""
		# Note that this function no longer returns the cycle; just True/False
		vertex2 = cython.declare(Vertex)
		edge = cython.declare(Edge)
		found = cython.declare(cython.bint)
		
		for vertex2, edge in self[chain[-1]].iteritems():
			if vertex2 is chain[0] and len(chain) > 2:
				return True
			elif vertex2 not in chain:
				# make the chain a little longer and explore again
				chain.append(vertex2)
				found = self.__isChainInCycle(chain)
				if found: return True
				# didn't find a cycle down this path (-vertex2),
				# so remove the vertex from the chain
				chain.remove(vertex2)
		return False

	def getAllCycles(self, startingVertex):
		"""
		Given a starting vertex, returns a list of all the cycles containing 
		that vertex.
		"""
		chain = cython.declare(list)
		cycleList = cython.declare(list)
		
		cycleList=list()
		chain = [startingVertex]
		
		#chainLabels=range(len(self.keys()))
		#print "Starting at %s in graph: %s"%(self.keys().index(startingVertex),chainLabels)
		
		cycleList = self.__exploreCyclesRecursively(chain, cycleList)
		return cycleList
		
	def __exploreCyclesRecursively(self, chain, cycleList):
		"""
		Finds cycles by spidering through a graph.
		Give it a chain of atoms that are connected, `chain`,
		and a list of cycles found so far `cycleList`.
		If `chain` is a cycle, it is appended to `cycleList`.
		Then chain is expanded by one atom (in each available direction)
		and the function is called again. This recursively spiders outwards
		from the starting chain, finding all the cycles.
		"""
		vertex2 = cython.declare(Vertex)
		edge = cython.declare(Edge)
		
		# chainLabels = cython.declare(list)
		# chainLabels=[self.keys().index(v) for v in chain] 
		# print "found %d so far. Chain=%s"%(len(cycleList),chainLabels)
		
		for vertex2, edge in self[chain[-1]].iteritems():
			# vertex2 will loop through each of the atoms 
			# that are bonded to the last atom in the chain.
			if vertex2 is chain[0] and len(chain) > 2:
				# it is the first atom in the chain - so the chain IS a cycle!
				cycleList.append(chain[:])
			elif vertex2 not in chain:
				# make the chain a little longer and explore again
				chain.append(vertex2)
				cycleList = self.__exploreCyclesRecursively(chain, cycleList)
				# any cycles down this path (-vertex2) have now been found,
				# so remove the vertex from the chain
				chain.pop(-1)
		return cycleList

	def setConnectivityValues(self):
		"""
		Set the connectivity values for each vertex in the graph. These are
		used to accelerate the isomorphism checking.
		"""
		
		count = cython.declare(cython.short)
		vert1 = cython.declare(Vertex)
		vert2 = cython.declare(Vertex)
		
		for vert1 in self:
			count = 0
			count = len(self[vert1])
			vert1.connectivity1 = count	
		for vert1 in self:
			count = 0
			for vert2 in self[vert1]: count += vert2.connectivity1
			vert1.connectivity2 = count
		for vert1 in self:
			count = 0
			for vert2 in self[vert1]: count += vert2.connectivity2
			vert1.connectivity3 = count
		
		
	def sortAndLabelVertices(self):
		"""
		Sort the vertices according to `globalAtomSortValue(atom)`,
		and record the sorting index on the vertices so they know what order 
		they were in. These are stored in `vertex.sorting_label`.
		"""
		# get the vertices into a list (so order is preserved)
		ordered_vertices = cython.declare(list)
		ordered_vertices = self.vertices()
		## sort the list according to something wise
		ordered_vertices.sort(key=globalAtomSortValue)
		# vertices with the same globalAtomSortValue can be in
		# an arbitary order, as long as it remains constant
		# so we record the ordering index ON the vertices
		i = cython.declare(cython.int)
		vertex = cython.declare(Vertex)
		for i in range(len(ordered_vertices)):
			vertex=ordered_vertices[i]
			vertex.sorting_label = i
			


################################################################################

def VF2_isomorphism(graph1, graph2, map12, map21, subgraph=False, findAll=False):
	"""
	Returns :data:`True` if two :class:`Graph` objects are isomorphic and
	:data:`False` otherwise. Uses the VF2 algorithm of Vento and Foggia. 

	If `subgraph` is :data:`True` then graph2 is checked for being a potential
	subgraph of graph1. i.e. graph1 is a specific case of graph2. 

	`findAll` is used to specify whether all isomorphisms should be returned, or
	only the first.

 	Returns tuple (is_match, map12, map21)  or if `findAll` is :data:`True` 
	then it returns tuple (is_match, list_of_map12s, list_of_map21s)"""

	map12List = cython.declare(list)
	map21List = cython.declare(list)
	ismatch = cython.declare(cython.bint)
	terminals1 = cython.declare(list)
	terminals2 = cython.declare(list)
	call_depth = cython.declare(cython.int)

	map12List = list()
	map21List = list()
	
	if not subgraph:
		if len(graph2) != len(graph1):
			#logging.debug("Tried matching graphs of different sizes!")
			return False, None, None # is_match, map12, map21
		elif len(graph2) == 0 == len(graph1):
			logging.warning("Tried matching empty graphs (returning True)")
			# This occurs (at least) when testing bond symmetry of H2
			return True, None, None
	else:
		if len(graph2)>len(graph1):
			#logging.debug("Tried matching small graph to larger subgraph")
			return False, None, None
	
	# start call_depth off as the size of the largest graph.
	# each recursive call to __VF2_match will decrease it by one,
	# until, when the whole graph has been explored, it should reach 0
	# it should never go below zero!
	call_depth = len(graph1)

	# update the connectivity values (before sorting by them)
	# sort the vertices according to something wise (based on connectivity value), and 
	# record the sorting order on each vertex (as vertex.sorting_label)
	cython.declare(vert=Vertex)
	vert = graph1.iterkeys().next() # just check the first atom in the graph
	if vert.sorting_label < 0:
		graph1.setConnectivityValues()
		graph1.sortAndLabelVertices()
	vert = graph2.iterkeys().next() # just check the first atom in the graph
	if vert.sorting_label < 0:
		graph2.setConnectivityValues()
		graph2.sortAndLabelVertices()
	
	# graph1.setConnectivityValues() # could probably run this less often elsewhere
	# graph2.setConnectivityValues() # as the values don't change often
	
	# graph1.sortAndLabelVertices()
	# graph2.sortAndLabelVertices()
	
	terminals1 = __VF2_terminals(graph1, map21)
	terminals2 = __VF2_terminals(graph2, map12)
	
	ismatch = __VF2_match(graph1, graph2, map21, map12, \
		 terminals1, terminals2, subgraph, findAll, map21List, map12List, call_depth)

	if findAll:
		return len(map21List) > 0, map21List, map12List
	else:
		return ismatch, map12, map21

def __VF2_feasible(graph1, graph2, vertex1, vertex2, map21, map12, terminals1,
	terminals2, subgraph):
	"""
	Returns :data:`True` if two vertices `vertex1` and `vertex2` from graphs
	`graph1` and `graph2`, respectively, are feasible matches. `mapping21` and
	`mapping12` are the current state of the mapping from `graph1` to `graph2`
	and vice versa, respectively. `terminals1` and `terminals2` are lists of
	the vertices that are directly connected to the already-mapped vertices.
	`subgraph` is :data:`True` if graph2 is to be treated as a potential
	subgraph of graph1. i.e. graph1 is a specific case of graph2.

	Uses the VF2 algorithm of Vento and Foggia. The feasibility is assessed
	through a series of semantic and structural checks. Only the combination
	of the semantic checks and the level 0 structural check are both
	necessary and sufficient to ensure feasibility. (This does *not* mean that
	vertex1 and vertex2 are always a match, although the level 1 and level 2
	checks preemptively eliminate a number of false positives.)
	"""

	vert1 = cython.declare(Vertex)
	vert2 = cython.declare(Vertex)
	edge1 = cython.declare(Edge)
	edge2 = cython.declare(Edge)
	edges1 = cython.declare(dict)
	edges2 = cython.declare(dict)
	i = cython.declare(cython.int)

	# Richard's Connectivity Value check
	# not sure where this is best done. Is it more specific or more general?
	if not subgraph: # then exact match required.
		if vertex1.connectivity1 != vertex2.connectivity1: return False
		if vertex1.connectivity2 != vertex2.connectivity2: return False
		if vertex1.connectivity3 != vertex2.connectivity3: return False
		
	# Semantic check #1: vertex1 and vertex2 must be equivalent
	if subgraph:
		if not vertex1.isSpecificCaseOf(vertex2):
			return False		
	else: # exact match required
		if not vertex1.equivalent(vertex2): 
		# Warning - I think current atom.equivalent returns True too often
			return False
		
	# Semantic check #2: adjacent vertices to vertex1 and vertex2 that are
	# already mapped should be connected by equivalent edges
	edges1 = graph1[vertex1]
	edges2 = graph2[vertex2]
		
	for vert2 in edges2:  # nb. 1 is specific case of 2
		if vert2 in map12:
			vert1 = map12[vert2]
			if not vert1 in edges1: # atoms not joined in graph1
				return False 
			edge1 = edges1[vert1]
			edge2 = edges2[vert2]
			if subgraph:
				if not edge1.isSpecificCaseOf(edge2):
					return False
			else: # exact match required
				if not edge1.equivalent(edge2):
				# Warning - I think bond.equivalent() returns True too often
					return False
	
	# there could still be edges in graph1 that aren't in graph2.
	# this is ok for subgraph matching, but not for exact matching
	if not subgraph:
		for vert1 in edges1:
			if vert1 in map21:
				vert2 = map21[vert1]
				if not vert2 in edges2: # atoms not joined in graph1
					return False 
				# if they exist in graph2, then 
				# we've already checked they're the same type as in graph1
	
	# Count number of terminals adjacent to vertex1 and vertex2
	term1Count = cython.declare(cython.int)
	term2Count = cython.declare(cython.int)
	neither1Count = cython.declare(cython.int)
	neither2Count = cython.declare(cython.int)

	term1Count = 0; term2Count = 0; neither1Count = 0; neither2Count = 0

	for vert1 in edges1:
		if vert1 in terminals1:
			term1Count += 1
		elif vert1 not in map21:
			neither1Count += 1
	for vert2 in edges2:
		if vert2 in terminals2:
			term2Count += 1
		elif vert2 not in map12:
			neither2Count += 1

	# Level 2 look-ahead: the number of adjacent vertices of vertex1 and
	# vertex2 that are non-terminals must be equal
	if subgraph:
		if neither1Count < neither2Count:
			return False
	else:
		if neither1Count != neither2Count:
			return False

	# Level 1 look-ahead: the number of adjacent vertices of vertex1 and
	# vertex2 that are terminals must be equal
	if subgraph:
		if term1Count < term2Count:
			return False
	else:
		if term1Count != term2Count:
			return False

	# Level 0 look-ahead: all adjacent vertices of vertex2 already in the
	# mapping must map to adjacent vertices of vertex1
	for vert2 in edges2:
		if vert2 in map12:
			vert1 = map12[vert2]
			if vert1 not in edges1:
				return False
	# ...AND all adjacent vertices of vertex1 already in the
	# mapping must map to adjacent vertices of vertex2, unless we are subgraph matching.
	if not subgraph:
		for vert1 in edges1:
			if vert1 in map21:
				vert2 = map21[vert1]
				if vert2 not in edges2:
					return False
	return True

def __VF2_match(graph1, graph2, map21, map12, terminals1, terminals2, subgraph,
	findAll, map21List, map12List, call_depth):
	"""
	A recursive function used to explore two graphs `graph1` and `graph2` for
	isomorphism by attempting to map them to one another. `mapping21` and
	`mapping12` are the current state of the mapping from `graph1` to `graph2`
	and vice versa, respectively. `terminals1` and `terminals2` are lists of
	the vertices that are directly connected to the already-mapped vertices.
	`subgraph` is :data:`True` if graph2 is to be treated as a potential
	subgraph of graph1. i.e. graph1 is a specific case of graph2.
	
	If findAll=True then it adds valid mappings to map21List and 
	map12List, but returns False when done (or True if the initial mapping is complete)

	Uses the VF2 algorithm of Vento and Foggia, which is O(N) in spatial complexity
	and O(N**2) (best-case) to O(N! * N) (worst-case) in temporal complexity.
	"""

	new_terminals1 = cython.declare(list)
	new_terminals2 = cython.declare(list)
	vertex1 = cython.declare(Vertex)
	vertex2 = cython.declare(Vertex)
	ismatch = cython.declare(cython.bint)
	pairs = cython.declare(list)

	if call_depth < 0:
		logging.error("Recursing too deep. Now %d"%call_depth)
		if call_depth < -100:
			raise Exception("Recursing infinitely deep!")
			logging.error('******************************************')
			logging.error("** RETURNING TRUE FOR ISOMORPHISM MATCH **")
			logging.error("**   JUST TO GET OUT OF INFINITE LOOP   **")
			logging.error('******************************************')
			return True
	
	# Done if we have mapped to all vertices in graph
	if len(map21) >= len(graph1):
		if findAll:
			map21List.append(map21.copy())
			map12List.append(map12.copy())
			#logging.verbose("Adding valid mapping to mapList")
		return True
	if len(map12) >= len(graph2) and subgraph:
		if findAll:
			map21List.append(map21.copy())
			map12List.append(map12.copy())
			#logging.verbose("Adding valid mapping to mapList")
		return True
	
	# Create list of pairs of candidates for inclusion in mapping
	pairs = __VF2_pairs(graph1, graph2, terminals1, terminals2, map21, map12)
	
	for vertex1, vertex2 in pairs:
		# propose a pairing
		if __VF2_feasible(graph1, graph2, vertex1, vertex2, map21, map12, \
				terminals1, terminals2, subgraph):
			# Update mapping accordingly
			map21[vertex1] = vertex2
			map12[vertex2] = vertex1
			
			# update terminals
			new_terminals1 = __VF2_new_terminals(graph1, map21, terminals1, vertex1)
			new_terminals2 = __VF2_new_terminals(graph2, map12, terminals2, vertex2)

			# Recurse
			ismatch = __VF2_match(graph1, graph2, \
				map21, map12, new_terminals1, new_terminals2, subgraph, findAll, \
				map21List, map12List, call_depth-1)
			if ismatch:
				if not findAll:
					return True
			# Undo proposed match
			del map21[vertex1]
			del map12[vertex2]
			# changes to 'new_terminals' will be discarded and 'terminals' is unchanged
			
	return False
	
def __VF2_pairs(graph1, graph2, terminals1, terminals2, map21, map12):
	"""
	Create a list of pairs of candidates for inclusion in the VF2 mapping. If
	there are a nonzero number of terminals in each graph, the candidates are
	selected to be one terminal from the first graph and all terminals from the
	second graph. If there are no terminals, the candidates are	selected to be
	one vertex from the first graph and all vertices from the second graph.
	"""
	pairs = cython.declare(list)
	vertex1 = cython.declare(Vertex)
	vertex2 = cython.declare(Vertex)
	terminal1 = cython.declare(Vertex)
	terminal2 = cython.declare(Vertex)
	list_to_sort = cython.declare(list)
	lowest_label = cython.declare(cython.short)
	this_label = cython.declare(cython.short)
	
	pairs = list()
	
	# Construct list from terminals if possible
	if len(terminals1) > 0 and len(terminals2) > 0:
		#list_to_sort = terminals2
		#list_to_sort.sort(key=__getSortLabel)
		#terminal2 = list_to_sort[0]
		terminal2 = terminals2[0]
		#list_to_sort = terminals1
		#list_to_sort.sort(key=__getSortLabel)
		
		for terminal1 in terminals1:
			pairs.append([terminal1, terminal2])
	# Otherwise construct list from all *remaining* vertices (not matched)
	else:
		# vertex2 is the lowest-labelled un-mapped vertex from graph2
		list_to_sort = graph2.keys()
		lowest_label = 32766 # hopefully we don't have more unmapped atoms than this!
		for vertex1 in list_to_sort: # just using vertex1 as a temporary variable
			this_label = vertex1.sorting_label
			if this_label < lowest_label:
				if not vertex1 in map12:
					lowest_label = this_label
					vertex2 = vertex1
			
		# pair with all vertex1s
		list_to_sort = graph1.keys() 
		list_to_sort.sort(key=__getSortLabel)
		for vertex1 in list_to_sort:
			if vertex1 not in map21: # exclude already mapped vertices
				pairs.append([vertex1, vertex2])
	
	return pairs

def __VF2_terminals(graph, mapping):
	"""
	For a given graph `graph` and associated partial mapping `mapping`,
	generate a list of terminals, vertices that are directly connected to
	vertices that have already been mapped.
	
	List is sorted (using key=__getSortLabel) before returning.
	"""
	
	terminals = cython.declare(list)
	vertex = cython.declare(Vertex)
	vert = cython.declare(Vertex)
	
	terminals = []# list()
	
	for vertex in mapping:
		for vert in graph[vertex]:
			if vert not in mapping:
				if vert not in terminals:
					terminals.append(vert)
	terminals.sort(key=__getSortLabel)
	return terminals

def __VF2_new_terminals(graph, mapping, old_terminals, new_vertex):
	"""
	For a given graph `graph` and associated partial mapping `mapping`,
	UPDATES a list of terminals, vertices that are directly connected to
	vertices that have already been mapped. You have to pass it the previous 
	list of terminals `old_terminals` and the vertex `vertex` that has been added 
	to the mapping. Returns a new COPY of the terminals.
	"""
	
	vertex = cython.declare(Vertex)
	vertex2 = cython.declare(Vertex)
	sorting_label = cython.declare(cython.short)
	sorting_label2 = cython.declare(cython.short)
	terminals = cython.declare(list)
	i = cython.declare(int)
	
	# copy the old terminals, leaving out the new_vertex
	#terminals = [v for v in old_terminals if not v is new_vertex]
	terminals = old_terminals[:]
	if new_vertex in terminals:
		terminals.remove(new_vertex)
	
	# Add the terminals of new_vertex
	for vertex in graph[new_vertex]:
		if vertex not in mapping: # only add if not already mapped
	
	## the next block of code is equivalent to these two lines:
	#		if vertex not in terminals: terminals.append(vertex)
	#terminals.sort(key=__getSortLabel)
			
			# find spot in the sorted terminals list where we should put this vertex
			sorting_label = vertex.sorting_label
			i=0; sorting_label2=-1 # in case terminals list empty
			for i in range(len(terminals)):
				vertex2 = terminals[i]
				sorting_label2 = vertex2.sorting_label
				if sorting_label2 >= sorting_label:
					break  
				# else continue going through the list of terminals
			else: # got to end of list without breaking, 
				# so add one to index to make sure vertex goes at end
				i+=1
			if sorting_label2 == sorting_label: # this vertex already in terminals.
				continue # try next vertex in graph[new_vertex]
			
			# insert vertex in right spot in terminals
			terminals.insert(i,vertex)
	
	return terminals

################################################################################
