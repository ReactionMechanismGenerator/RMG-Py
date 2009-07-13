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

def VF2_isomorphism(graph1, graph2, map21, map12, subgraph=False, findAll=False):
	"""
	Returns :data:`True` if two graphs are isomorphic and :data:`False`
	otherwise. Uses the VF2 algorithm of Vento and Foggia. `subgraph` is 
	:data:`True` if graph2 is a potential subgraph of graph1. `findAll` is
	used to specify whether all isomorphisms should be returned, or only the
	first.
	"""
	if map12 is None: map12 = {}
	if map21 is None: map21 = {}

	map12List = []; map21List = []
	terminals1 = __VF2_terminals(graph1, map21)
	terminals2 = __VF2_terminals(graph2, map12)
	
	ismatch = __VF2_match(graph1, graph2, map21, map12, \
		terminals1, terminals2, subgraph, findAll, map21List, map12List)
	
	if findAll:
		return len(map21List) > 0, map21List, map12List
	else:
		return ismatch, map21, map12

def __VF2_feasible(graph1, graph2, vertex1, vertex2, \
	map21, map12, terminals1, terminals2, subgraph):
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
	edges1 = graph1[vertex1]
	edges2 = graph2[vertex2]
	
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
	term1Count = 0; term2Count = 0
	neither1Count = 0; neither2Count = 0
	
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

def __VF2_match(graph1, graph2, map21, map12, \
	terminals1, terminals2, subgraph, findAll, \
	map21List, map12List):
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
	pairs = __VF2_pairs(graph1, graph2, terminals1, terminals2)
	
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

def __VF2_pairs(graph1, graph2, terminals1, terminals2):
	"""
	Create a list of pairs of candidates for inclusion in the VF2 mapping. If
	there are a nonzero number of terminals in each graph, the candidates are
	selected to be one terminal from the first graph and all terminals from the
	second graph. If there are no terminals, the candidates are	selected to be
	one vertex from the first graph and all vertices from the second graph.
	"""

	pairs = []
	
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

def __VF2_terminals(graph, mapping):
	"""
	For a given graph `graph` and associated partial mapping `mapping`, 
	generate a list of terminals, vertices that are directly connected to
	vertices that have already been mapped.
	"""

	terminals = {}
	for vertex in mapping:
		for vert, edge in graph[vertex].iteritems():
			if vert not in mapping:
				terminals[vert] = True
	return terminals

	