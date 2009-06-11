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
Contains functions for working with the graph data type, in particular functions
for graph and subgraph isomorphism comparisons. To use the functions in this
module, the graphs must be represented as dictionaries of dictionaries, where
both the outer dictionary and inner dictionaries use vertices as keys. Semantic
checks are implemented by calls to :meth:`equivalent()` methods of vertices
and edges.

The primary isomorphism algorithm is the VF2 algorithm of Vento and Foggia, 
which is O(N) in spatial complexity	and O(N**2) (best-case) to O(N! * N) 
(worst-case) in temporal complexity. For more information see the following
references:	

L. P. Cordella, P. Foggia, C. Sansone, and M. Vento. "Performance Evaluation
of the VF Graph Matching Algorithm." Proc. of the 10th ICIAP, IEEE Computer
Society Press, p. 1172-1177 (1999).

L. P. Cordella, P. Foggia, C. Sansone, and M. Vento. "An Improved Algorithm for
Matching Large Graphs." 3rd IAPR-TC15 Workshop on Graph-based Representations 
in Pattern Recognition, Cuen, p. 149-156 (2007).
"""

class Vertex:
	"""
	A generic vertex (node) of a graph. The vertex is required to have a
	:meth:`equivalent()` function for comparing semantic properties between
	pairs of vertices.
	"""
	def __init__(self, color):
		self.color = color
		
	def equivalent(self, other):
		"""
		Used to determine if two vertices have equivalent semantic properties.
		In this	case, two vertices are equivalent if they have the same color.
		"""
		return self.color == other.color
		
################################################################################

class Edge:
	"""
	A generic edge of a graph. The edge is required to have a
	:meth:`equivalent()` function for comparing semantic properties between
	pairs of edges.
	"""
	
	def __init__(self, color):
		self.color = color
		
	def equivalent(self, other):
		"""
		Used to determine if two edges have equivalent semantic properties.
		In this	case, two edges are equivalent if they have the same color.
		"""
		return self.color == other.color
		
################################################################################

class Graph:
	"""
	A representation of a graph data structure. Internally the graph is 
	represented as a dictionary of dictionaries. If a vertex is in the graph
	it will be in the outer dictionary. If two vertices in the graph are
	connected by an edge, each edge will be in the inner dictionary.
	"""

	def __init__(self):
		self.graph = {}
		
	def addVertex(self, vertex):
		"""
		Add `vertex` to the graph as a vertex. The vertex is initialized with
		no edges.
		"""
		self.graph[vertex] = {}
		return vertex
		
	def addEdge(self, vertex1, vertex2, edge):
		"""
		Add `edge` to the graph as an edge connecting vertices `vertex1` and
		`vertex2`, which must already be present in the graph.
		"""
		self.graph[vertex1][vertex2] = edge
		self.graph[vertex2][vertex1] = edge
		return edge
		
	def hasEdge(self, vertex1, vertex2):
		"""
		Returns true if vertices `vertex1` and `vertex2`, are in the graph and
		are connected by an edge.
		"""
		if vertex1 in self.graph.keys():
			if vertex2 in self.graph[vertex1].keys():
				return True
		return False
	
	def isIsomorphic(self, other):
		"""
		Returns :data:`True` if two graphs are isomorphic and :data:`False`
		otherwise. Uses the VF2 algorithm of Vento and Foggia.
		"""
		return VF2_isomorphic(self.graph, other.graph, False)
		
	def isSubgraphIsomorphic(self, other, map12={}, map21={}):
		"""
		Returns :data:`True` if `other` is subgraph isomorphic and :data:`False`
		otherwise. Uses the VF2 algorithm of Vento and Foggia.
		"""
		return VF2_isomorphic(self.graph, other.graph, True, False)

	def findSubgraphIsomorphisms(self, other):
		"""
		Returns :data:`True` if `other` is subgraph isomorphic and :data:`False`
		otherwise. Uses the VF2 algorithm of Vento and Foggia.
		"""
		return VF2_isomorphic(self.graph, other.graph, True, True)


################################################################################

class ChemGraph:
	"""
	A representation of a chemical species or fragment using a graph data 
	structure. The vertices represent atoms, while the edges represent bonds.
	This class can be used to represent a resonance form of a chemical species
	or a functional group. Atom iteration is possible via the `atoms` method,
	while bond iteration is possible via the `bonds` method.
	
	Internally the graph is represented as a dictionary of dictionaries. If a 
	vertex is in the graph it will be in the outer dictionary. If two vertices 
	in the graph are connected by an edge, each edge will be in the inner 
	dictionary.
	"""

	def __init__(self, atoms=[], bonds=[]):
		self.initialize(atoms, bonds)
		
	def atoms(self):
		"""
		Return a list of the atoms in the structure.
		"""
		return self.graph.keys()
		
	def bonds(self):
		"""
		Return a list of the bonds in the structure.
		"""
		bondlist = []
		for atom1 in self.graph:
			for atom2 in self.graph[atom1]:
				bond = self.graph[atom1][atom2]
				if bond not in bondlist:
					bondlist.append(bond)
		return bondlist
	
	def addAtom(self, atom):
		"""
		Add `atom` to the graph as a vertex. The atom is initialized with
		no edges.
		"""
		self.graph[atom] = {}
		
	def addBond(self, bond):
		"""
		Add `bond` to the graph as an edge connecting atoms `atom1` and
		`atom2`, which must already be present in the graph.
		"""
		atom1 = bond.atoms[0]; atom2 = bond.atoms[1]
		self.graph[atom1][atom2] = bond
		self.graph[atom2][atom1] = bond
		
	def hasBond(self, atom1, atom2):
		"""
		Returns true if atoms `atom1` and `atom2`, are in the graph and
		are connected by a bond.
		"""
		if atom1 in self.graph.keys():
			if atom2 in self.graph[atom1].keys():
				return True
		return False
	
	def removeAtom(self, atom):
		"""
		Remove `atom` from the graph as a vertex. Also removes all bonds
		associated with `atom`. Does not remove atoms that no longer have any
		bonds as a result of this removal.
		"""
		if atom not in self.graph: return
		for atom2 in self.graph:
			if atom2 is not atom:
				if atom in self.graph[atom2]:
					del self.graph[atom2][atom]
		del self.graph[atom]

	def removeBond(self, bond):
		"""
		Remove `bond` from the graph. Does not remove atoms that no longer have
		any bonds as a result of this removal.
		"""
		atom1 = bond.atoms[0]; atom2 = bond.atoms[1]
		del self.graph[atom1][atom2]
		del self.graph[atom2][atom1]

	def isIsomorphic(self, other):
		"""
		Returns :data:`True` if two graphs are isomorphic and :data:`False`
		otherwise. Uses the VF2 algorithm of Vento and Foggia.
		"""
		return VF2_isomorphic(self.graph, other.graph, False, False)
		
	def isSubgraphIsomorphic(self, other, map12={}, map21={}):
		"""
		Returns :data:`True` if `other` is subgraph isomorphic and :data:`False`
		otherwise. Uses the VF2 algorithm of Vento and Foggia.
		"""
		return VF2_isomorphic(self.graph, other.graph, True, False, map12, map21)
	
	def findSubgraphIsomorphisms(self, other):
		"""
		Returns :data:`True` if `other` is subgraph isomorphic and :data:`False`
		otherwise. Uses the VF2 algorithm of Vento and Foggia.
		"""
		return VF2_isomorphic(self.graph, other.graph, True, True)

	def initialize(self, atoms, bonds):
		"""
		Rebuild the `graph` data member based on the lists of atoms and bonds
		provided in `atoms` and `bonds`, respectively.
		"""
		self.graph = {}
		for atom in atoms:
			self.graph[atom] = {}
		for bond in bonds:
			self.graph[bond.atoms[0]][bond.atoms[1]] = bond
			self.graph[bond.atoms[1]][bond.atoms[0]] = bond
	
################################################################################

def VF2_isomorphic(graph1, graph2, subgraph=False, findAll=False, mapping12={}, mapping21={}):	
	"""
	Returns :data:`True` if two graphs are isomorphic and :data:`False`
	otherwise. Uses the VF2 algorithm of Vento and Foggia. `subgraph` is 
	:data:`True` if graph2 is a potential subgraph of graph1. `findAll` is
	used to specify whether all isomorphisms should be returned, or only the
	first.
	"""
	mappings12 = []; mappings21 = []
	terminals1 = __VF2_terminals(graph1, mapping21)
	terminals2 = __VF2_terminals(graph2, mapping12)
	ismatch, map12, map21 = __VF2_match(graph1, graph2, mapping21, mapping12, terminals1, terminals2, subgraph, findAll, mappings21, mappings12)
	if findAll:
		return (mappings12 > 0), mappings12, mappings21
	else:
		return ismatch, map12, map21

def __VF2_feasible(graph1, graph2, vertex1, vertex2, mapping21, mapping12, \
                 terminals1, terminals2, subgraph):
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
		if vert1 in mapping21.keys():
			vert2 = mapping21[vert1]
			if not vert2 in edges2.keys():
				return False
			edge2 = edges2[vert2]
			if not edge1.equivalent(edge2):
				return False
	
	# Count number of terminals adjacent to vertex1 and vertex2
	term1Count = 0; term2Count = 0
	neither1Count = 0; neither2Count = 0
	for vert1, edge1 in edges1.iteritems():
		if vert1 in terminals1.keys():
			term1Count += 1
		elif vert1 not in mapping21.keys():
			neither1Count += 1
	for vert2, edge2 in edges2.iteritems():
		if vert2 in terminals2.keys():
			term2Count += 1
		elif vert2 not in mapping12.keys():
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
		if vert1 in mapping21.keys():
			vert2 = mapping21[vert1]
			if vert2 not in edges2.keys():
				return False
	
	return True

def __VF2_match(graph1, graph2, mapping21, mapping12, terminals1, terminals2, subgraph, findAll, mappingList21=[], mappingList12=[]):
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
	if len(mapping12) == len(graph2) or len(mapping21) == len(graph1):
		return True, mapping12, mapping21
		
	# Create list of pairs of candidates for inclusion in mapping
	pairs = __VF2_pairs(graph1, graph2, terminals1, terminals2)
	
	for vertex1, vertex2 in pairs:
		if __VF2_feasible(graph1, graph2, vertex1, vertex2, mapping21, mapping12, terminals1, terminals2, subgraph):
			# Update mapping and terminals accordingly
			mapping21[vertex1] = vertex2
			mapping12[vertex2] = vertex1
			terminals1 = __VF2_terminals(graph1, mapping21)
			terminals2 = __VF2_terminals(graph2, mapping12)
			# Recurse
			ismatch, mapping12_0, mapping21_0 = __VF2_match(graph1, graph2, \
					mapping21, mapping12, terminals1, terminals2, subgraph, findAll, mappingList21, mappingList12)
			if ismatch:
				if findAll:
					mappingList21.append(mapping21_0.copy())
					mappingList12.append(mapping12_0.copy())
				else:
					return True, mapping12_0, mapping21_0
			# Undo proposed match
			del mapping21[vertex1]
			del mapping12[vertex2]
			terminals1 = __VF2_terminals(graph1, mapping21)
			terminals2 = __VF2_terminals(graph2, mapping12)


	return False, [], []

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
		for terminal2 in terminals2:
			if len(pairs) == 0:
				for terminal1 in terminals1:
					pairs.append([terminal1, terminal2])
	# Otherwise construct list from all remaining vertices
	else:
		for vertex2 in graph2:
			if vertex2 in graph2.keys() and len(pairs) == 0:
				for vertex1 in graph1:
					if vertex1 in graph1.keys():
						pairs.append([vertex1, vertex2])

	return pairs

def __VF2_terminals(graph, mapping):
	"""
	For a given graph `graph` and associated partial mapping `mapping`, 
	generate a list of terminals, vertices that are directly connected to
	vertices that have already been mapped.
	"""

	terminals = {}
	for vertex in graph:
		if vertex in mapping:
			for vert, edge in graph[vertex].iteritems():
				if vert not in mapping:
					terminals[vert] = True
	return terminals
		
################################################################################

if __name__ == '__main__':

	#graph1 = Graph()
	#vertex1 = graph1.addVertex(Vertex('red'))
	#vertex2 = graph1.addVertex(Vertex('blue'))
	#vertex3 = graph1.addVertex(Vertex('red'))
	#vertex4 = graph1.addVertex(Vertex('blue'))
	#vertex5 = graph1.addVertex(Vertex('red'))
	#vertex6 = graph1.addVertex(Vertex('blue'))
	#edge1 = graph1.addEdge(vertex1, vertex2, Edge('black'))
	#edge2 = graph1.addEdge(vertex2, vertex3, Edge('white'))
	#edge3 = graph1.addEdge(vertex3, vertex4, Edge('black'))
	#edge4 = graph1.addEdge(vertex4, vertex5, Edge('white'))
	#edge5 = graph1.addEdge(vertex5, vertex6, Edge('black'))
	#edge6 = graph1.addEdge(vertex6, vertex1, Edge('white'))
	
	#graph2 = Graph()
	#vertex4 = graph2.addVertex(Vertex('blue'))
	#vertex5 = graph2.addVertex(Vertex('blue'))
	#vertex6 = graph2.addVertex(Vertex('blue'))
	#vertex1 = graph2.addVertex(Vertex('red'))
	#vertex2 = graph2.addVertex(Vertex('red'))
	#vertex3 = graph2.addVertex(Vertex('red'))
	#edge3 = graph2.addEdge(vertex5, vertex2, Edge('black'))
	#edge1 = graph2.addEdge(vertex4, vertex1, Edge('black'))
	#edge5 = graph2.addEdge(vertex6, vertex3, Edge('black'))
	#edge4 = graph2.addEdge(vertex3, vertex5, Edge('white'))
	#edge2 = graph2.addEdge(vertex2, vertex4, Edge('white'))
	#edge6 = graph2.addEdge(vertex6, vertex1, Edge('white'))
	
	graph1 = Graph()
	vertex1 = graph1.addVertex(Vertex('red'))
	vertex2 = graph1.addVertex(Vertex('blue'))
	vertex3 = graph1.addVertex(Vertex('red'))
	vertex4 = graph1.addVertex(Vertex('blue'))
	vertex5 = graph1.addVertex(Vertex('red'))
	vertex6 = graph1.addVertex(Vertex('blue'))
	edge1 = graph1.addEdge(vertex1, vertex2, Edge('black'))
	edge2 = graph1.addEdge(vertex2, vertex3, Edge('white'))
	edge3 = graph1.addEdge(vertex3, vertex4, Edge('black'))
	edge4 = graph1.addEdge(vertex4, vertex1, Edge('white'))
	edge5 = graph1.addEdge(vertex4, vertex5, Edge('black'))
	edge6 = graph1.addEdge(vertex4, vertex6, Edge('white'))
	
	graph2 = Graph()
	vertex4 = graph2.addVertex(Vertex('blue'))
	vertex5 = graph2.addVertex(Vertex('blue'))
	#vertex6 = graph2.addVertex(Vertex('blue'))
	vertex1 = graph2.addVertex(Vertex('red'))
	vertex2 = graph2.addVertex(Vertex('red'))
	#vertex3 = graph2.addVertex(Vertex('red'))
	edge1 = graph2.addEdge(vertex1, vertex4, Edge('black'))
	edge2 = graph2.addEdge(vertex4, vertex2, Edge('white'))
	edge3 = graph2.addEdge(vertex2, vertex5, Edge('black'))
	edge4 = graph2.addEdge(vertex5, vertex1, Edge('white'))
	#edge5 = graph2.addEdge(vertex5, vertex3, Edge('black'))
	#edge6 = graph2.addEdge(vertex5, vertex6, Edge('white'))
	
	
	#print graph1.isSubgraphIsomorphic(graph2)
	print graph1.findSubgraphIsomorphisms(graph2)
