#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
#	RMG - Reaction Mechanism Generator
#
#	Copyright (c) 2002-2009 Prof. William H. Green (whgreen@mit.edu) and the
#	RMG Team (rmg_dev@mit.edu)
#
#	Permission is hereby granted, free of charge, to any person obtaining a
#	copy of this software and associated documentation files (the 'Software'),
#	to deal in the Software without restriction, including without limitation
#	the rights to use, copy, modify, merge, publish, distribute, sublicense,
#	and/or sell copies of the Software, and to permit persons to whom the
#	Software is furnished to do so, subject to the following conditions:
#
#	The above copyright notice and this permission notice shall be included in
#	all copies or substantial portions of the Software.
#
#	THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#	FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#	DEALINGS IN THE SOFTWARE.
#
################################################################################

"""
Contains an implementation of a graph data structure (the :class:`Graph` class)
and functions for manipulating that graph, including isomorphism functions.
"""

import chem

import cython
import logging

################################################################################

class Vertex:

	def __init__(self):
		self.connectivity = [0, 0, 0]
	
	def equivalent(self, other):
		return True

################################################################################

class Edge:

	def __init__(self):
		pass

	def equivalent(self, other):
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

	def __reduce__(self):
		return (Graph, (), None, None, self.iteritems())
		
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

	def isIsomorphic(self, other, map12_0, map21_0):
		"""
		Returns :data:`True` if two graphs are isomorphic and :data:`False`
		otherwise. Uses the VF2 algorithm of Vento and Foggia.
		"""
		if len(self) != len(other): return False
		ismatch, map21, map12 = VF2_isomorphism(self, other, map21_0, map12_0, False, False)
		return ismatch

	def isSubgraphIsomorphic(self, other, map12_0, map21_0):
		"""
		Returns :data:`True` if `other` is subgraph isomorphic and :data:`False`
		otherwise. Uses the VF2 algorithm of Vento and Foggia.
		"""
		ismatch, map21, map12 = VF2_isomorphism(self, other, map21_0, map12_0, True, False)
		return ismatch

	def findSubgraphIsomorphisms(self, other, map12_0, map21_0):
		"""
		Returns :data:`True` if `other` is subgraph isomorphic and :data:`False`
		otherwise. Uses the VF2 algorithm of Vento and Foggia.
		"""
		return VF2_isomorphism(self, other, map21_0, map12_0, True, True)

	def copy(self):
		"""
		Create a copy of the current graph.
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
		vertex = cython.declare(chem.Atom)
		rootVertex = cython.declare(chem.Atom)
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
		vertex2 = cython.declare(chem.Atom)
		edge = cython.declare(chem.Bond)
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
		vertex2 = cython.declare(chem.Atom)
		edge = cython.declare(chem.Bond)
		chainLabels = cython.declare(list)

		chainLabels=[self.keys().index(v) for v in chain]
		#print "found %d so far. Chain=%s"%(len(cycleList),chainLabels)
		
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
		Sets the Extended Connectivity values as introduced by Morgan (1965)
		http://dx.doi.org/10.1021/c160017a018
		
		First CV1 is the number of neighbours  (Morgan proposed non-Hydrogen neighbours)
		CV2 is the sum of neighbouring CV1 values
		CV3 is the sum of neighbouring CV2 values
		"""
		
		count = cython.declare(cython.int)
		vert1 = cython.declare(chem.Atom)
		vert2 = cython.declare(chem.Atom)
		
		for i in range(3):
			for vert1 in self:
				count = 0
				if i == 0:
					count = len(self[vert1])
				else:
					for vert2 in self[vert1]: count += vert1.connectivity[i-1]
				vert1.connectivity[i] = count
			
			
	def sort_and_label_vertices(self):
		"""
		Sort the vertices according to something wise,
		and record the sorting index on the vertices so they know what order 
		they were in. These are stored in `vertex.sorting_label`.
		"""
		# get the vertices into a list (so order is preserved)
		ordered_vertices = cython.declare(list)
		ordered_vertices = self.vertices()
		## sort the list according to something wise
		ordered_vertices.sort(key=global_atom_sort_value)
		# vertices with the same __global_atom_sort_value can be in 
		# an arbitary order, as long as it remains constant
		# so we record the ordering index ON the vertices
		for i in range(len(ordered_vertices)):
			(ordered_vertices[i]).sorting_label = i
    	
################################################################################

def VF2_isomorphism(graph1, graph2, map12, map21, subgraph=False, findAll=False):
	"""
	Returns :data:`True` if two :class:`Graph` objects are isomorphic and
	:data:`False` otherwise. Uses the VF2 algorithm of Vento and Foggia. If
	`subgraph` is :data:`True` then graph2 is checked for being a potential
	subgraph of graph1. `findAll` is used to specify whether all isomorphisms
	should be returned,  or only the first.
	
	Returns tuple (is_match, map12, map21)
	"""

	map12List = cython.declare(list)
	map21List = cython.declare(list)
	ismatch = cython.declare(cython.bint)
	terminals1 = cython.declare(dict)
	terminals2 = cython.declare(dict)
	call_depth = cython.declare(cython.int)

	map12List = list()
	map21List = list()
	
	if not subgraph:
		if len(graph2) != len(graph1):
			logging.debug("Tried matching graphs of different sizes!")
			return False, None, None # is_match, map12, map21
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
	graph1.setConnectivityValues() # could probably run this less often elsewhere
	graph2.setConnectivityValues() # as the values don't change often
	
	# sort the vertices according to something wise (based on connectivity value), and 
	# record the sorting order on each vertex (as vertex.sorting_label)
	graph1.sort_and_label_vertices()
	graph2.sort_and_label_vertices()
	
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
	subgraph of graph1.

	Uses the VF2 algorithm of Vento and Foggia. The feasibility is assessed
	through a series of semantic and structural checks. Only the combination
	of the semantic checks and the level 0 structural check are both
	necessary and sufficient to ensure feasibility. (This does *not* mean that
	vertex1 and vertex2 are always a match, although the level 1 and level 2
	checks preemptively eliminate a number of false positives.)
	"""

	vert1 = cython.declare(chem.Atom)
	vert2 = cython.declare(chem.Atom)
	edge1 = cython.declare(chem.Bond)
	edge2 = cython.declare(chem.Bond)
	edges1 = cython.declare(dict)
	edges2 = cython.declare(dict)

	# Richard's Connectivity Value check
	# not sure where this is best done. Is it more specific or more general?
	if not subgraph:
		for i in range(3):
			if vertex1.connectivity[i] != vertex2.connectivity[i]: return False
		
	# Semantic check #1: vertex1 and vertex2 must be equivalent
	if not vertex1.equivalent(vertex2):
		return False
	
	# Semantic check #2: adjacent vertices to vertex1 and vertex2 that are
	# already mapped should be connected by equivalent edges
	edges1 = graph1[vertex1]
	edges2 = graph2[vertex2]
		
	for vert1 in edges1:
	# for vert1, edge1 in edges1.iteritems(): # if you uncomment this..**
		if vert1 in map21:
			vert2 = map21[vert1]
			if not vert2 in edges2:
				return False
			edge1 = edges1[vert1] # **..then remove this
			edge2 = edges2[vert2]
			if not edge1.equivalent(edge2):
				return False
	
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

	# Level 0 look-ahead: all adjacent vertices of vertex1 already in the
	# mapping must map to adjacent vertices of vertex2
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
	subgraph of graph1.

	Uses the VF2 algorithm of Vento and Foggia, which is O(N) in spatial complexity
	and O(N**2) (best-case) to O(N! * N) (worst-case) in temporal complexity.
	"""

	new_terminals1 = cython.declare(dict)
	new_terminals2 = cython.declare(dict)
	vertex1 = cython.declare(chem.Atom)
	vertex2 = cython.declare(chem.Atom)
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
	
	# Done if we have mapped to all vertices in graph2
	if len(map12) >= len(graph2) or len(map21) >= len(graph1):
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
				if findAll:
					map21List.append(map21.copy())
					map12List.append(map12.copy())
				else:
					return True
			# Undo proposed match
			del map21[vertex1]
			del map12[vertex2]
			# changes to 'new_terminals' will be discarded and 'terminals' is unchanged

	return False

def __get_sort_label(vertex):
	"""
	Used to sort vertices prior to poposing candidate pairs in :method:`__VF2_pairs`

	This returns the `sorting_label` that is stored on the vertex. It should have been
	put there recently by a call to :method:`Graph.sort_and_label_vertices()`
	"""
	return vertex.sorting_label
		
def global_atom_sort_value(atom):
	"""
	Used to sort atoms prior to poposing candidate pairs in :method:`__VF2_pairs` 
	The lowest (or most negative) values will be first in the list when you sort, 
	so should be given to the atoms you want to explore first. 
	For now, that is (roughly speaking) the most connected atoms. This definitely helps for large graphs
	but bizarrely the opposite ordering seems to help small graphs. Not sure about subggraphs...
	"""
	#return hash(atom)  # apparently random?
	#return (atom.connectivity[0] ) # not unique enough
	return ( - (atom.connectivity[0]<<8)  # left shift 8 = multiply by 2**8=256
			 - (atom.connectivity[1]<<4)  # left shift 4 = multiply by 2**4=16
			 - (atom.connectivity[2])
			)
	
def __VF2_pairs(graph1, graph2, terminals1, terminals2, map21, map12):
	"""
	Create a list of pairs of candidates for inclusion in the VF2 mapping. If
	there are a nonzero number of terminals in each graph, the candidates are
	selected to be one terminal from the first graph and all terminals from the
	second graph. If there are no terminals, the candidates are	selected to be
	one vertex from the first graph and all vertices from the second graph.
	"""

	pairs = cython.declare(list)
	vertex1 = cython.declare(chem.Atom)
	vertex2 = cython.declare(chem.Atom)
	terminal1 = cython.declare(chem.Atom)
	terminal2 = cython.declare(chem.Atom)
	list_to_sort = cython.declare(list)

	pairs = list()

	# Construct list from terminals if possible
	if len(terminals1) > 0 and len(terminals2) > 0:
		list_to_sort = terminals2.keys()
		list_to_sort.sort(key=__get_sort_label)
		terminal2 = list_to_sort[0]
		list_to_sort = terminals1.keys()
		list_to_sort.sort(key=__get_sort_label)
		
		for terminal1 in list_to_sort:
			pairs.append([terminal1, terminal2])
	# Otherwise construct list from all *remaining* vertices (not matched)
	else:
		list_to_sort = graph2.keys()
		# remove already mapped vertices
		for vertex2 in map12:
			list_to_sort.remove(vertex2)
		list_to_sort.sort(key=__get_sort_label)
		vertex2 = list_to_sort[0]  # take first vertex2
		# pair with all vertex1s
		list_to_sort = graph1.keys() 
		list_to_sort.sort(key=__get_sort_label)
		for vertex1 in list_to_sort:
			if vertex1 not in map21: # exclude already mapped vertices
				pairs.append([vertex1, vertex2])
	
	return pairs

def __VF2_terminals(graph, mapping):
	"""
	For a given graph `graph` and associated partial mapping `mapping`,
	generate a list of terminals, vertices that are directly connected to
	vertices that have already been mapped.
	"""

	terminals = cython.declare(dict)
	vertex = cython.declare(chem.Atom)
	vert = cython.declare(chem.Atom)

	terminals = dict()

	for vertex in mapping:
		for vert in graph[vertex]:
			if vert not in mapping:
				terminals[vert] = True
	return terminals

def __VF2_new_terminals(graph, mapping, old_terminals, new_vertex):
	"""
	For a given graph `graph` and associated partial mapping `mapping`,
	UPDATES a list of terminals, vertices that are directly connected to
	vertices that have already been mapped. You have to pass it the previous 
	list of terminals `old_terminals` and the vertex `vertex` that has been added 
	to the mapping. Returns a new copy of the terminals.
	"""
	
	terminals = cython.declare(dict)
	terminals = dict()

	# copy the old terminals, leaving out the new_vertex
	for vertex in old_terminals:
		if not vertex is new_vertex: 
			terminals[vertex] = True
	
	# add the terminals of new_vertex
	for vertex in graph[new_vertex]:
		if vertex not in mapping:
			terminals[vertex] = True
			
	return terminals

################################################################################

