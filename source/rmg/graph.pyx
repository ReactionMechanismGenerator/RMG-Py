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

cdef extern from "dictobject.h":
	ctypedef class __builtin__.dict [object PyDictObject]:
		pass

cdef class Vertex:

	def __init__(self):
		pass

	cpdef bint equivalent(Vertex self, Vertex other):
		return True

cdef class Edge:

	def __init__(self):
		pass

	cpdef bint equivalent(Edge self, Edge other):
		return True

cdef class Graph(dict):
	"""
	A representation of a graph using a dictionary of dictionaries. The keys
	of the outer dictionary are the vertices, while edges are accessed via
	self[vertex1][vertex2].
	"""

	def __init__(self, vertices=None, edges=None):
		self.clear()

	cpdef list vertices(Graph self):
		"""
		Return a list of the vertices in the graph.
		"""
		return self.keys()
	
	cpdef list edges(Graph self):
		"""
		Return a list of the edges in the graph.
		"""
		cdef list edgelist = list()
		cdef list pairslist = list()
		for v1 in self:
			for v2 in self[v1]:
				if (v1, v2) not in pairslist:
					edgelist.append(self[v1][v2])
					pairslist.append((v1,v2))
					pairslist.append((v2,v1))
		return edgelist
	
	cpdef addVertex(Graph self, vertex):
		"""
		Add a `vertex` to the graph. The vertex is initialized with no edges.
		"""
		self[vertex] = dict()
		return vertex

	cpdef addEdge(Graph self, vertices, edge):
		"""
		Add an `edge` to the graph as an edge connecting the two vertices
		specified in the 2-tuple `vertices`.
		"""
		v1, v2 = vertices
		self[v1][v2] = edge
		self[v2][v1] = edge
		return edge

	cpdef dict getEdges(Graph self, vertex):
		"""
		Return a list of the edges involving the specified `vertex`.
		"""
		return self[vertex]

	cpdef getEdge(Graph self, tuple vertices):
		"""
		Returns the edge connecting vertices in 2-tuple `vertices`, or None if
		no edge exists.
		"""
		v1, v2 = vertices
		return self[v1][v2] if self.hasEdge(vertices) else None

	cpdef bint hasEdge(self, tuple vertices):
		"""
		Returns :data:`True` if vertices in the 2-tuple `vertices` are connected
		by an edge, and :data:`False` if not.
		"""
		v1, v2 = vertices
		if v1 in self:
			return v2 in self[v1]
		return False

	cpdef removeVertex(Graph self, vertex1):
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

	cpdef removeEdge(Graph self, vertices):
		"""
		Remove the edge having vertices as specified in the 2-tuple `vertices
		from the graph. Does not remove vertices that no longer have any edges
		as a result of this removal.
		"""
		v1, v2 = vertices
		del self[v1][v2]
		del self[v2][v1]

	cpdef isIsomorphic(Graph self, Graph other, dict map12_0, dict map21_0):
		"""
		Returns :data:`True` if two graphs are isomorphic and :data:`False`
		otherwise. Uses the VF2 algorithm of Vento and Foggia.
		"""
		if len(self) != len(other): return False
		ismatch, map21, map12 = VF2_isomorphism(self, other, map21_0, map12_0, False, False)
		return ismatch

	cpdef isSubgraphIsomorphic(Graph self, Graph other, dict map12_0, dict map21_0):
		"""
		Returns :data:`True` if `other` is subgraph isomorphic and :data:`False`
		otherwise. Uses the VF2 algorithm of Vento and Foggia.
		"""
		ismatch, map21, map12 = VF2_isomorphism(self, other, map21_0, map12_0, True, False)
		return ismatch

	cpdef findSubgraphIsomorphisms(Graph self, Graph other, dict map12_0, dict map21_0):
		"""
		Returns :data:`True` if `other` is subgraph isomorphic and :data:`False`
		otherwise. Uses the VF2 algorithm of Vento and Foggia.
		"""
		return VF2_isomorphism(self, other, map21_0, map12_0, True, True)

	cpdef Graph copy(Graph self):
		"""
		Create a copy of the current graph.
		"""
		cdef Graph other = Graph()

		for vertex in self:
			other.addVertex(vertex)
		for v1 in self:
			for v2 in self[v1]:
				other[v1][v2] = self[v1][v2]
		
		return other

	cpdef Graph merge(Graph self, Graph other):
		"""
		Merge two graphs so as to store them in a single Graph object.
		"""

		# Create output graph
		cdef Graph new = Graph()

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

	cpdef list split(Graph self):
		"""
		Convert a single Graph object containing two or more unconnected graphs
		into separate graphs.
		"""
		
		# Create potential output graphs
		cdef Graph new1 = self.copy()
		cdef Graph new2 = Graph()
		cdef list verticesToMove
		cdef int index

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

	cpdef list getSmallestSetOfSmallestRings(Graph self):
		"""
		Return a list of the smallest set of smallest rings in the graph. The
		algorithm implements was adapted from a description by Fan, Panaye,
		Doucet, and Barbu (doi: 10.1021/ci00015a002)

		B. T. Fan, A. Panaye, J. P. Doucet, and A. Barbu. "Ring Perception: A
		New Algorithm for Directly Finding the Smallest Set of Smallest Rings
		from a Connection Table." *J. Chem. Inf. Comput. Sci.* **33**, 
		p. 657-662 (1993).
		"""

		# Make a copy of the graph so we don't modify the original
		cdef Graph graph = self.copy()
		cdef bint done
		cdef list verticesToMove
		cdef list cycleList
		cdef list cycles

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
			found, cycle = graph.isVertexInCycle(vertex)
			if not found:
				verticesToRemove.append(vertex)
		# Remove identified vertices from graph
		for vertex in verticesToRemove:
			graph.removeVertex(vertex)

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
				cycles = list()
				graph.getAllCycles(rootVertex, cycles)
				if len(cycles) == 0:
					raise Exception('Did not find expected cycle!')

				# Keep the smallest of the cycles found above
				cycle = cycles[0]
				for c in cycles[1:]:
					if len(c) < len(cycle):
						cycle = c
				cycleList.append(cycle)

				# Remove all vertices in the cycle from the graph that have only two edges
				verticesToRemove = []
				for vertex in cycle:
					if len(graph[vertex]) <= 2:
						verticesToRemove.append(vertex)
				if len(verticesToRemove) == 0:
					# Remove edge between root vertex and any vertex it is connected to
					graph.removeEdge((rootVertex, graph[rootVertex].keys()[0]))
				else:
					for vertex in verticesToRemove:
						graph.removeVertex(vertex)


		return cycleList

	cpdef isVertexInCycle(Graph self, cycle):

		if not isinstance(cycle, list): cycle = [cycle]

		for vertex2, edge in self[cycle[-1]].iteritems():
			if vertex2 is cycle[0] and len(cycle) > 2:
				return True, cycle
			elif vertex2 not in cycle:
				cycle.append(vertex2)
				found, c = self.isVertexInCycle(cycle)
				if found: return True, cycle
				cycle.remove(vertex2)

		return False, []


	cpdef getAllCycles(Graph self, cycle, cycleList):

		if not isinstance(cycle, list): cycle = [cycle]

		for vertex2, edge in self[cycle[-1]].iteritems():
			if vertex2 is cycle[0] and len(cycle) > 2:
				cycleList.append(cycle[:])
			elif vertex2 not in cycle:
				cycle.append(vertex2)
				self.getAllCycles(cycle, cycleList)
				cycle.pop(-1)

################################################################################

cpdef VF2_isomorphism(Graph graph1, Graph graph2, dict map12, dict map21, \
	bint subgraph=False, bint findAll=False):
	"""
	Returns :data:`True` if two graphs are isomorphic and :data:`False`
	otherwise. Uses the VF2 algorithm of Vento and Foggia. `subgraph` is
	:data:`True` if graph2 is a potential subgraph of graph1. `findAll` is
	used to specify whether all isomorphisms should be returned, or only the
	first.
	"""

	cdef list map12List = list(), map21List = list()
	cdef dict terminals1, terminals2

	terminals1 = __VF2_terminals(graph1, map21)
	terminals2 = __VF2_terminals(graph2, map12)

	ismatch = __VF2_match(graph1, graph2, map21, map12, \
		terminals1, terminals2, subgraph, findAll, map21List, map12List)

	if findAll:
		return len(map21List) > 0, map21List, map12List
	else:
		return ismatch, map12, map21

cdef __VF2_feasible(Graph graph1, Graph graph2, object vertex1, object vertex2, \
	dict map21, dict map12, dict terminals1, dict terminals2, bint subgraph):
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
	cdef dict edges1 = graph1[vertex1]
	cdef dict edges2 = graph2[vertex2]

	# Semantic check #1: vertex1 and vertex2 must be equivalent
	if not vertex1.equivalent(vertex2):
		return False

	# Semantic check #2: adjacent vertices to vertex1 and vertex2 that are
	# already mapped should be connected by equivalent edges
	for vert1, edge1 in edges1.iteritems():
		if vert1 in map21:
			vert2 = map21[vert1]
			if not vert2 in edges2:
				return False
			edge2 = edges2[vert2]
			if not edge1.equivalent(edge2):
				return False

	# Count number of terminals adjacent to vertex1 and vertex2
	cdef int term1Count = 0, term2Count = 0, neither1Count = 0, neither2Count = 0
	
	for vert1, edge1 in edges1.iteritems():
		if vert1 in terminals1:
			term1Count += 1
		elif vert1 not in map21:
			neither1Count += 1
	for vert2, edge2 in edges2.iteritems():
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
	for vert1, edge1 in edges1.iteritems():
		if vert1 in map21:
			vert2 = map21[vert1]
			if vert2 not in edges2:
				return False

	return True

cdef bint __VF2_match(Graph graph1, Graph graph2, dict map21, dict map12, \
	dict terminals1, dict terminals2, bint subgraph, bint findAll, \
	list map21List, list map12List):
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

	# Done if we have mapped to all vertices in graph2
	if len(map12) >= len(graph2) or len(map21) >= len(graph1):
		return True

	# Create list of pairs of candidates for inclusion in mapping
	cdef list pairs = __VF2_pairs(graph1, graph2, terminals1, terminals2)

	for vertex1, vertex2 in pairs:
		if __VF2_feasible(graph1, graph2, vertex1, vertex2, map21, map12, \
				terminals1, terminals2, subgraph):
			# Update mapping and terminals accordingly
			map21[vertex1] = vertex2
			map12[vertex2] = vertex1
			terminals1 = __VF2_terminals(graph1, map21)
			terminals2 = __VF2_terminals(graph2, map12)
			# Recurse
			ismatch = __VF2_match(graph1, graph2, \
				map21, map12, terminals1, terminals2, subgraph, findAll, \
				map21List, map12List)
			if ismatch:
				if findAll:
					map21List.append(map21.copy())
					map12List.append(map12.copy())
				else:
					return True
			# Undo proposed match
			del map21[vertex1]
			del map12[vertex2]
			terminals1 = __VF2_terminals(graph1, map21)
			terminals2 = __VF2_terminals(graph2, map12)

	return False

cdef list __VF2_pairs(Graph graph1, Graph graph2, dict terminals1, dict terminals2):
	"""
	Create a list of pairs of candidates for inclusion in the VF2 mapping. If
	there are a nonzero number of terminals in each graph, the candidates are
	selected to be one terminal from the first graph and all terminals from the
	second graph. If there are no terminals, the candidates are	selected to be
	one vertex from the first graph and all vertices from the second graph.
	"""

	cdef list pairs = list()

	# Construct list from terminals if possible
	if len(terminals1) > 0 and len(terminals2) > 0:
		terminal2 = terminals2.keys()[0]
		for terminal1 in terminals1:
			pairs.append([terminal1, terminal2])
	# Otherwise construct list from all remaining vertices
	else:
		vertex2 = graph2.keys()[0]
		for vertex1 in graph1:
			pairs.append([vertex1, vertex2])

	return pairs

cdef dict __VF2_terminals(Graph graph, dict mapping):
	"""
	For a given graph `graph` and associated partial mapping `mapping`,
	generate a list of terminals, vertices that are directly connected to
	vertices that have already been mapped.
	"""

	cdef dict terminals = dict()
	for vertex in mapping:
		for vert, edge in graph[vertex].iteritems():
			if vert not in mapping:
				terminals[vert] = True
	return terminals

################################################################################

