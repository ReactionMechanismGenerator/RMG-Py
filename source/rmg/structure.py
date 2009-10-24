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
Contains classes describing chemical structures, which are made up of atoms and
bonds.
"""

import logging
import pybel
import openbabel

import chem
import graph

################################################################################

class InvalidAdjacencyListException(Exception):
	"""
	An exception used when parsing an adjacency list to indicate that it is
	invalid. The label of the adjacency list is stored in the `label` attribute.
	"""

	def __init__(self, label):
		self.label = label

	def __str__(self):
		return 'Invalid adjacency list: ' + self.label


################################################################################

class Structure:
	"""
	A representation of a chemical species using a graph data structure. The
	vertices represent atoms, while the edges represent bonds. The attributes
	are:

	================  ==========================================================
	Attributes        Description
	================  ==========================================================
	`graph`           A graph representation of the atoms and bonds of the
	                  structure
	`symmetryNumber`  The combined external and internal symmetry number of the
	                  structure
	================  ==========================================================

	"""
	
	def __init__(self, atoms=None, bonds=None, SMILES=None):
		self.initialize(atoms or [], bonds or [])
		if SMILES:
			self.fromSMILES(SMILES)

	def __repr__(self):
		message = "Structure(SMILES='%s')"%(self.toSMILES())
		return message

	def atoms(self):
		"""
		Return a list of the atoms in the structure.
		"""
		return self.graph.vertices()

	def bonds(self):
		"""
		Return a list of the bonds in the structure.
		"""
		return self.graph.edges()

	def addAtom(self, atom):
		"""
		Add `atom` to the graph as a vertex. The atom is initialized with
		no edges.
		"""
		return self.graph.addVertex(atom)

	def addBond(self, bond):
		"""
		Add `bond` to the graph as an edge connecting atoms `atom1` and
		`atom2`, which must already be present in the graph.
		"""
		atom1, atom2 = bond.atoms
		return self.graph.addEdge((atom1, atom2), bond)

	def getBonds(self, atom):
		"""
		Return a dictionary of neighbouring atoms and their bonds to the specified `atom`.
		
		The keys in the dictionary are the neighbouring atoms.
		The values are the bonds themselves.
		`getBonds(atom).values()` will return just the bonds
		
		"""
		return self.graph.getEdges(atom)

	def getBond(self, atom1, atom2):
		"""
		Returns the bond connecting atoms `atom1` and `atom2` if it exists, or
		:data:`None` if not.
		"""
		return self.graph.getEdge((atom1, atom2))

	def hasBond(self, atom1, atom2):
		"""
		Returns true if atoms `atom1` and `atom2`, are in the graph and
		are connected by a bond.
		"""
		return self.graph.hasEdge((atom1, atom2))

	def removeAtom(self, atom):
		"""
		Remove `atom` from the graph as a vertex. Also removes all bonds
		associated with `atom`. Does not remove atoms that no longer have any
		bonds as a result of this removal.
		"""
		self.graph.removeVertex(atom)

	def removeBond(self, bond):
		"""
		Remove `bond` from the graph. Does not remove atoms that no longer have
		any bonds as a result of this removal.
		"""
		atom1, atom2 = bond.atoms
		return self.graph.removeEdge((atom1, atom2))

	def isIsomorphic(self, other, map12=None, map21=None):
		"""
		Returns :data:`True` if two graphs are isomorphic and :data:`False`
		otherwise. Uses the VF2 algorithm of Vento and Foggia.
		"""
		if map12 is None: map12 = dict()
		if map21 is None: map21 = dict()
		return self.graph.isIsomorphic(other.graph, map12, map21)

	def isSubgraphIsomorphic(self, other, map12=None, map21=None):
		"""
		Returns :data:`True` if `other` is subgraph isomorphic and :data:`False`
		otherwise. Uses the VF2 algorithm of Vento and Foggia.
		"""
		if map12 is None: map12 = dict()
		if map21 is None: map21 = dict()
		return self.graph.isSubgraphIsomorphic(other.graph, map12, map21)

	def findSubgraphIsomorphisms(self, other, map12=None, map21=None):
		"""
		Returns :data:`True` if `other` is subgraph isomorphic and :data:`False`
		otherwise. Uses the VF2 algorithm of Vento and Foggia.
		"""
		if map12 is None: map12 = dict()
		if map21 is None: map21 = dict()
		return self.graph.findSubgraphIsomorphisms(other.graph, map12, map21)

	def initialize(self, atoms, bonds):
		"""
		Rebuild the `graph` data member based on the lists of atoms and bonds
		provided in `atoms` and `bonds`, respectively.
		"""
		self.graph = graph.Graph()

		if atoms is None or bonds is None:
			return

		for atom in atoms:
			self.addAtom(atom)
		for bond in bonds:
			self.addBond(bond)

	def copy(self):
		"""
		Create a copy of the current Structure.
		"""
		atoms = []; bonds = []
		for atom in self.atoms():
			atoms.append(atom.copy())
		for bond in self.bonds():
			newBond = bond.copy()
			bonds.append(newBond)
			index1 = self.atoms().index(bond.atoms[0])
			index2 = self.atoms().index(bond.atoms[1])
			newBond.atoms = [atoms[index1], atoms[index2]]

		return Structure(atoms, bonds)

	def merge(self, other):
		"""
		Merge two graphs so as to store them in a single Graph object.
		"""
		structure = Structure()
		structure.graph = self.graph.merge(other.graph)
		return structure

	def split(self):
		"""
		Convert a single Graph object containing two or more unconnected graphs
		into separate graphs.
		"""
		graphs = self.graph.split()
		structures = []
		for g in graphs:
			structure = Structure()
			structure.graph = g
			structures.append(structure)
		return structures

	def getSmallestSetOfSmallestRings(self):
		"""
		Return the smallest set of smallest rings for the structure.
		"""
		return self.graph.getSmallestSetOfSmallestRings()

	def getFormula(self):
		"""
		Return the molecular formula for the structure.
		"""
		mol = pybel.Molecule(self.toOBMol())
		return mol.formula

	def fromAdjacencyList(self, adjlist):
		"""
		Convert a string adjacency list `adjlist` into a structure object.
		"""

		try:

			atoms = []; bonds = []; atomdict = {}; bonddict = {}

			lines = adjlist.splitlines()

			label = lines[0]

			for line in lines[1:]:

				data = line.split()

				# Skip if blank line
				if len(data) == 0:
					continue

				# First item is index for atom
				# Sometimes these have a trailing period (as if in a numbered list),
				# so remove it just in case
				aid = int(data[0].strip('.'))

				# If second item is '*', the atom is the center atom
				center = ''
				index = 1
				if data[1][0] == '*':
					center = data[1]
					index = 2

				# Next is the element or functional group element
				# A list can be specified with the {,} syntax
				elements = data[index]
				if elements[0] == '{':
					elements = elements[1:-1].split(',')
				else:
					elements = [elements]

				# Next is the electron state
				elecState = data[index+1].upper()
				if elecState[0] == '{':
					elecState = elecState[1:-1].split(',')
				else:
					elecState = [elecState]

				# Create a new atom based on the above information
				atom = chem.Atom(elements, elecState, 0, center)

				# Add the atom to the list
				atoms.append(atom)
				atomdict[aid] = atom

				bonddict[aid] = {}

				# Process list of bonds
				for datum in data[index+2:]:

					# Sometimes commas are used to delimit bonds in the bond list,
					# so strip them just in case
					datum = datum.strip(',')

					aid2, comma, btype = datum[1:-1].partition(',')
					aid2 = int(aid2)

					if btype[0] == '{':
						btype = btype[1:-1].split(',')
					else:
						btype = [btype]

					if aid2 in atomdict:
						bond = chem.Bond([atomdict[aid], atomdict[aid2]], btype)
						bonds.append(bond)

					bonddict[aid][aid2] = btype

		except Exception, e:
			raise InvalidAdjacencyListException(label)

		# Check consistency using bonddict
		for atom1 in bonddict:
			for atom2 in bonddict[atom1]:
				if atom2 not in bonddict:
					raise InvalidAdjacencyListException(label)
				elif atom1 not in bonddict[atom2]:
					raise InvalidAdjacencyListException(label)
				elif bonddict[atom1][atom2] != bonddict[atom2][atom1]:
					raise InvalidAdjacencyListException(label)


		# Create and return functional group or species
		self.initialize(atoms, bonds)



	def fromCML(self, cmlstr):
		"""
		Convert a string of CML `cmlstr` to a Structure object.
		"""
		cmlstr = cmlstr.replace('\t', '')
		mol = pybel.readstring('cml', cmlstr)
		self.fromOBMol(mol.OBMol)

	def fromInChI(self, inchistr):
		"""
		Convert an InChI string `inchistr` to a Structure object.
		"""
		mol = pybel.readstring('inchi', inchistr)
		self.fromOBMol(mol.OBMol)

	def fromSMILES(self, smilesstr):
		"""
		Convert a SMILES string `smilesstr` to a Structure object.
		"""
		mol = pybel.readstring('smiles', smilesstr)
		self.fromOBMol(mol.OBMol)

	def fromOBMol(self, obmol):
		"""
		Convert an OpenBabel OBMol object `obmol` to a Structure object.
		"""

		atoms = []; bonds = []

		# Add hydrogen atoms to complete molecule if needed
		obmol.AddHydrogens()

		# Iterate through atoms in obmol
		for i in range(0, obmol.NumAtoms()):
			obatom = obmol.GetAtom(i + 1)

			# Use atomic number as key for element
			number = obatom.GetAtomicNum()

			# Process spin multiplicity
			electron = obatom.GetSpinMultiplicity()
			if electron == 0: electron = '0'
			elif electron == 1:	electron = '2S'
			elif electron == 2:	electron = '1'
			elif electron == 3:	electron = '2T'

			atom = chem.Atom(chem.elements[number].symbol, chem.electronStates[electron])
			atoms.append(atom)

			# Add bonds by iterating again through atoms
			for j in range(0, i):
				obatom2 = obmol.GetAtom(j + 1)
				obbond = obatom.GetBond(obatom2)
				if obbond is not None:
					order = ''

					# Process bond type
					if obbond.IsSingle(): order = 'S'
					elif obbond.IsDouble(): order = 'D'
					elif obbond.IsTriple(): order = 'T'
					elif obbond.IsAromatic(): order = 'B'

					bond = chem.Bond([atoms[i], atoms[j]], chem.bondTypes[order])
					bonds.append(bond)

		# Create the graph from the atom and bond lists
		self.initialize(atoms, bonds)

	def toAdjacencyList(self, label=''):
		"""
		Convert the structure object to an adjacency list. The `label` parameter
		is an optional string to put as the first line of the adjacency list;
		if set to the empty string, this line will be omitted.
		"""

		adjlist = ''

		if label != '': adjlist += label + '\n'

		atoms = self.atoms()

		for i, atom in enumerate(atoms):

			# Atom number
			adjlist += str(i+1) + ' '

			# Atom label
			if atom.label != '':
				adjlist += atom.label + ' '

			# Atom type(s)
			if atom.atomType.__class__ == list:
				adjlist += '{' + atom.atomType[0].label
				for atomType in atom.atomType[1:]:
					adjlist += ',' + atomType.label
				adjlist += '} '
			else:
				adjlist += atom.atomType.label + ' '

			# Electron state(s)
			if atom.electronState.__class__ == list:
				adjlist += '{' + atom.electronState[0].label
				for electronState in atom.electronState[1:]:
					adjlist += ',' + electronState.label
				adjlist += '} '
			else:
				adjlist += atom.electronState.label + ' '

			# Bonds list
			for atom2, bond in self.getBonds(atom).iteritems():
				adjlist += '{' + str(atoms.index(atom2)+1) + ','

				# Bond type(s)
				if bond.bondType.__class__ == list:
					adjlist += '{' + bond.bondType[0].label
					for bondType in bond.bondType[1:]:
						adjlist += ',' + bondType.label
					adjlist += '}'
				else:
					adjlist += bond.bondType.label

				adjlist += '} '

			# Each atom begins on a new list
			adjlist += '\n'

		return adjlist


	def toOBMol(self):
		"""
		Convert a Structure object to an OpenBabel OBMol object.
		"""
		atoms = self.atoms(); bonds = self.bonds()

		obmol = openbabel.OBMol()
		for atom in atoms:
			a = obmol.NewAtom()
			a.SetAtomicNum(atom.atomType.element.number)
		for bond in bonds:
			index1 = atoms.index(bond.atoms[0])
			index2 = atoms.index(bond.atoms[1])
			order = bond.bondType.order
			if order == 1.5: order = 5
			obmol.AddBond(index1+1, index2+1, int(order))

		obmol.AssignSpinMultiplicity(True)

		return obmol

	def toCML(self):
		"""
		Convert a Structure object to CML.
		"""
		mol = pybel.Molecule(self.toOBMol())
		cml = mol.write('cml').strip()
		return '\n'.join([l for l in cml.split('\n') if l.strip()])

	def toInChI(self):
		"""
		Convert a Structure object to an InChI string.
		"""
		# This version does not write a warning to stderr if stereochemistry is undefined
		obmol = self.toOBMol()
		obConversion = openbabel.OBConversion()
		obConversion.SetOutFormat('inchi')
		obConversion.SetOptions('w', openbabel.OBConversion.OUTOPTIONS)
		return obConversion.WriteString(obmol).strip()
		# This version writes a warning to stderr if stereochemistry is undefined
		#mol = pybel.Molecule(self.toOBMol())
		#return mol.write('inchi').strip()

	def toSMILES(self):
		"""
		Convert a Structure object to an SMILES string.
		"""
		mol = pybel.Molecule(self.toOBMol())
		return mol.write('smiles').strip()

	def toXML(self, dom, root):
		"""
		Convert a Structure object to an XML DOM tree.
		"""
		import xml.dom.minidom
		cml = dom.createElement('cml')
		root.appendChild(cml)

		dom2 = xml.dom.minidom.parseString(self.toCML())
		cml.appendChild(dom2.documentElement)

	def toDOT(self, label=''):
		"""
		Convert a Structure object to a graph image via DOT. This is useful
		for visualizing the functional groups that make up the databases. The
		output is a string containing a graph in DOT format, which can be
		passed to graphviz to produce an image; neato is recommended.

		Atoms are visualized as vertices in the outputted graph. Vertices are
		labeled with the atom type(s) of each corresponding	atom. Labeled atoms
		('\*', '\*1', etc.) are color-coded, with a unique color for each label.
		Bonds are indicated with edges; multiple bonds are represented by
		multiple edges between the same pair of vertices. The edge line style is
		used to denote further semantic information: dashed lines indicate
		optional higher-order bonds, while dotted lines indicate benzene bonds.
		"""

		import pydot

		graph = pydot.Dot(size='5,4', rankdir='LR',
				graph_type='graph', simplify=True, fontsize='8',
				overlap='true', dpi='85',center="True")

		# List of atoms (vertices)
		for i, atom in enumerate(self.atoms()):
			# Generate vertex label from atom type labels
			label = ','.join([atomType.label for atomType in atom._atomType])
			# Labeled atoms are color coded
			color = 'black'
			if atom.label != '':
				colors = {'*': 'red', '*1': 'red', '*2': 'green', '*3': 'blue', '*4': 'yellow', '*5': 'purple', '*6': 'orange', '*7': 'magenta', '*8': 'cyan'}
				color = colors[atom.label]
			# Create and add vertex to graph
			node = pydot.Node(str(i+1), label=label, color=color, fontcolor=color)
			graph.add_node(node)
			
		# List of bonds (edges)
		for i, bond in enumerate(self.bonds()):
			index1 = self.atoms().index(bond.atoms[0])
			index2 = self.atoms().index(bond.atoms[1])
			
			single = False; double = False; triple = False; benzene = False
			for type in bond._bondType:
				if type.order == 1: single = True
				if type.order == 2: double = True
				if type.order == 3: triple = True
				if type.order == 1.5: benzene = True

			label = []
			if single: label.append('S')
			if double: label.append('D')
			if triple: label.append('T')
			if benzene: label.append('B')
			label = ','.join(label)

			# Create and add edge to graph
			edge = pydot.Edge(str(index1+1), str(index2+1), label=label, len='2')
			graph.add_edge(edge)

		return graph
	
	def simplifyAtomTypes(self):
		"""
		Iterate through the atoms in the structure, setting them to be equal
		to their element.
		"""
		for atom1 in self.atoms():
			# Only works for single atom types, not lists
			if atom1.atomType.__class__ == list:
				continue
			# Skip generic atom types
			if atom1.atomType.label == 'R' or atom1.atomType.label == 'R!H':
				continue
			# Reset atom type to that of element
			atom1.atomType = chem.atomTypes[atom1.atomType.element.symbol]

	def updateAtomTypes(self):
		"""
		Iterate through the atoms in the structure, checking their atom types
		to ensure they are correct (i.e. accurately describe their local bond
		environment) and complete (i.e. are as detailed as possible).
		"""

		# NOTE: Does not yet process CO atom type!

		for atom1 in self.atoms():
			# Only works for single atom types, not lists
			if atom1.atomType.__class__ == list:
				continue
			# Skip generic atom types
			if atom1.atomType.label == 'R' or atom1.atomType.label == 'R!H':
				continue
			# Count numbers of each bond type
			single = 0; double = 0; triple = 0; benzene = 0; carbonyl = False
			for atom2, bond12 in self.getBonds(atom1).iteritems():
				if bond12.isSingle(): single += 1
				elif bond12.isDouble(): double += 1
				elif bond12.isTriple(): triple += 1
				elif bond12.isBenzene(): benzene += 1
				if atom1.isCarbon() and atom2.isOxygen() and bond12.isDouble():
					carbonyl = True

			# Use counts to determine proper atom type
			atomType = atom1.atomType.element.symbol
			if atomType == 'C':
				if triple == 1: atomType = 'Ct'
				elif single == 3 or single == 4: atomType = 'Cs'
				elif double == 2: atomType = 'Cdd'
				elif carbonyl: atomType = 'CO'
				elif double == 1: atomType = 'Cd'
				elif benzene == 1 or benzene == 2: atomType = 'Cb'
				elif benzene == 3: atomType = 'Cbf'
			elif atomType == 'O':
				if single == 1 or single == 2: atomType = 'Os'
				elif double == 1: atomType = 'Od'

			# Do nothing if suggested and specified atom types are identical
			if atom1.atomType.label == atomType:
				pass
			# Do nothing if specified atom type is 'Cbf' and suggested is 'Cb'
			elif atom1.atomType.label == 'Cbf' and atomType == 'Cb':
				pass
			# Do nothing if specified atom type is 'Cdd' and suggested is 'CO'
			elif atom1.atomType.label == 'Cdd' and atomType == 'CO':
				pass
			# Make change if specified atom type is element
			elif atom1.atomType.label == atom1.atomType.element.symbol:
				#logging.warning('Changed "' + atom1.atomType.label + '" to "' + atomType + '".')
				atom1.atomType = chem.atomTypes[atomType]
			# Else print warning
			else:
				logging.warning('Suggested atom type "' + atomType + '" does not match specified atom type "' + atom1.atomType.label + '".')

	def getRadicalCount(self):
		"""
		Get the number of radicals in the structure.
		"""
		radical = 0
		for atom in self.atoms():
			radical += atom.electronState.order
		return radical

	def getAdjacentResonanceIsomers(self):
		"""
		Generate all of the resonance isomers formed by one allyl radical shift.
		"""

		isomers = []

		# Radicals
		if self.getRadicalCount() > 0:
			# Iterate over radicals in structure
			for atom in self.atoms():
				paths = self.findAllDelocalizationPaths(atom)
				for path in paths:

					atom1, atom2, atom3, bond12, bond23 = path

					# Adjust to (potentially) new resonance isomer
					atom1.decreaseFreeElectron()
					atom3.increaseFreeElectron()
					bond12.increaseOrder()
					bond23.decreaseOrder()

					# Make a copy of isomer
					isomer = self.copy()

					# Restore current isomer
					atom1.increaseFreeElectron()
					atom3.decreaseFreeElectron()
					bond12.decreaseOrder()
					bond23.increaseOrder()

					# Append to isomer list if unique
					isomers.append(isomer)

		return isomers

	def findAllDelocalizationPaths(self, atom1):
		"""
		Find all the delocalization paths allyl to the radical center indicated
		by `atom1`. Used to generate resonance isomers.
		"""

		# No paths if atom1 is not a radical
		if atom1.electronState.order <= 0:
			return []

		# Find all delocalization paths
		paths = []
		for atom2, bond12 in self.getBonds(atom1).iteritems():
			# Vinyl bond must be capable of gaining an order
			if bond12.canIncreaseOrder():
				for atom3, bond23 in self.getBonds(atom2).iteritems():
					# Allyl bond must be capable of losing an order without
					# breaking
					if atom1 is not atom3 and bond23.canDecreaseOrder():
						paths.append([atom1, atom2, atom3, bond12, bond23])
		return paths

	def clearLabeledAtoms(self):
		"""
		Remove the labels from all atoms in the structure.
		"""
		for atom in self.atoms():
			atom.label = ''

	def containsLabeledAtom(self, label):
		"""
		Return :data:`True` if the structure contains an atom with the label
		`label` and :data:`False` otherwise.
		"""
		for atom in self.atoms():
			if atom.label == label: return True
		return False

	def getLabeledAtom(self, label):
		"""
		Return the atoms in functional group structure that are labeled, i.e.
		the center atoms in the structure.
		"""
		for atom in self.atoms():
			if atom.label == label: return atom
		return None

	def getLabeledAtoms(self):
		"""
		Return the atoms in functional group structure that are labeled, i.e.
		the center atoms in the structure.
		"""
		atoms = {}
		for atom in self.atoms():
			if atom.isCenter(): atoms[atom.label] = atom
		return atoms

	def isAtomInCycle(self, atom):
		"""
		Return :data:`True` if `atom` is in one or more cycles in the structure,
		and :data:`False` if not.
		"""
		return self.graph.isVertexInCycle(atom)

	def isBondInCycle(self, bond):
		"""
		Return :data:`True` if `atom` is in one or more cycles in the structure,
		and :data:`False` if not.
		"""
		atom1 = bond.atoms[0]
		atom2 = bond.atoms[1]
		cycle_list = self.graph.getAllCycles(atom1)
		for cycle in cycle_list:
			if atom2 in cycle:
				return True
		return False
		
		# You could have both atoms of the bond be in separate cycles,
		# eg. C6H11-C6H11
		# /\_/\
		# || ||
		# \/ \/
		# so this doesn't work: 
		## return self.isAtomInCycle(bond.atoms[0]) and self.isAtomInCycle(bond.atoms[1])

	def getAllCycles(self, atom):
		"""
		Given a starting `atom`, return a list of all cycles this atom is found
		in.
		"""
		return self.graph.getAllCycles(atom)

	def isCyclic(self):
		"""
		Return :data:`True` if one or more cycles are present in the structure
		and :data:`False` otherwise.
		"""
		for atom in self.atoms():
			if self.isAtomInCycle(atom):
				return True
		return False

	def calculateNumberOfRotors(self):
		"""
		Return the number of rotors in the molecule.
		"""
		# # this is tempting, but returns the number of rotors that lead to 
		# # different conformers, not the total number of internal rotors
		# obmol = self.toOBMol()
		# rotors = obmol.NumRotors()
		# return rotors
		
		if self.isLinear(): return 0
		rotors = 0
		for bond in self.bonds():
			if not bond.isSingle(): continue # only count single bonds
			if self.isBondInCycle(bond): continue # ring bonds no good either
			
			for atom in bond.atoms:
				# an atom has only one bond, so is terminal
				if len(self.graph[atom])==1: break
			else: # didn't break
				rotors+=1
		return rotors
		
		
	def isLinear(self):
		"""
		Return :data:`True` if the structure is linear and :data:`False` otherwise.
		"""
		atoms = self.atoms()
		atomCount = len(atoms)
		
		# True if two atoms
		if atomCount == 2:
			return True
		# False if one atom
		if atomCount == 1:
			return False
		# False if cyclic
		if self.isCyclic():
			return False
			
		# True if all bonds are double bonds
		for bond in self.bonds(): 
			if not bond.isDouble(): break
		else: # didn't break
			return True
		
		# True if alternating single-triple bonds
		for atom in atoms:
			#print 'atom',atom.atomType.label
			bonds = self.getBonds(atom).values()
			#print 'bonds',[b.bondType.label for b in bonds]
			if len(bonds)==1:
				continue # ok, next atom
			if len(bonds)>2:
				break # fail!
			if bonds[0].isSingle() and bonds[1].isTriple():
				continue # ok, next atom
			if bonds[1].isSingle() and bonds[0].isTriple():
				continue # ok, next atom
			break # fail if we haven't continued
		else:
			# didn't fail
			return True
			
		# not returned yet? must be nonlinear
		return False

	def calculateAtomSymmetryNumber(self, atom):
		"""
		Return the symmetry number centered at `atom` in the structure. The
		`atom` of interest must not be in a cycle.
		"""
		symmetryNumber = 1

		single = 0; double = 0; triple = 0; benzene = 0
		numNeighbors = 0
		for atom2, bond in self.getBonds(atom).iteritems():
			if bond.isSingle(): single += 1
			elif bond.isDouble(): double += 1
			elif bond.isTriple(): triple += 1
			elif bond.isBenzene(): benzene += 1

			numNeighbors += 1
		
		# If atom has zero or one neighbors, the symmetry number is 1
		if numNeighbors < 2: return symmetryNumber

		# Create temporary structures for each functional group attached to atom
		structure = Structure(self.atoms(), self.bonds())
		bondsToRemove = [bond for atom2, bond in structure.graph[atom].iteritems()]
		for bond in bondsToRemove: structure.removeBond(bond)
		structure.removeAtom(atom)
		groups = structure.split()

		# Determine equivalence of functional groups around atom
		groupIsomorphism = dict([(group, dict()) for group in groups])
		for group1 in groups:
			for group2 in groups:
				if group1 is not group2 and group2 not in groupIsomorphism[group1]:
					groupIsomorphism[group1][group2] = group1.isIsomorphic(group2)
					groupIsomorphism[group2][group1] = groupIsomorphism[group1][group2]
				elif group1 is group2:
					groupIsomorphism[group1][group1] = True
		count = [sum([int(groupIsomorphism[group1][group2]) for group2 in groups]) for group1 in groups]
		for i in range(count.count(2) / 2):
			count.remove(2)
		for i in range(count.count(3) / 3):
			count.remove(3); count.remove(3)
		for i in range(count.count(4) / 4):
			count.remove(4); count.remove(4); count.remove(4)
		count.sort(); count.reverse()
		
		if atom.getFreeElectronCount() == 0:
			if single == 4:
				# Four single bonds
				if count == [4]: symmetryNumber *= 12
				elif count == [3, 1]: symmetryNumber *= 3
				elif count == [2, 2]: symmetryNumber *= 2
				elif count == [2, 1, 1]: symmetryNumber *= 1
				elif count == [1, 1, 1, 1]: symmetryNumber *= 1
			elif single == 2:
				# Two single bonds
				if count == [2]: symmetryNumber *= 2
			elif double == 2:
				# Two double bonds
				if count == [2]: symmetryNumber *= 2
		elif atom.getFreeElectronCount() == 1:
			if single == 3:
				# Three single bonds
				if count == [3]: symmetryNumber *= 6
				elif count == [2, 1]: symmetryNumber *= 2
				elif count == [1, 1, 1]: symmetryNumber *= 1
		elif atom.getFreeElectronCount() == 2:
			if single == 2:
				# Two single bonds
				if count == [2]: symmetryNumber *= 2

		return symmetryNumber

	def calculateBondSymmetryNumber(self, bond):
		"""
		Return the symmetry number centered at `bond` in the structure.
		"""
		symmetryNumber = 1
		if bond.isSingle() or bond.isDouble() or bond.isTriple():
			if bond.atoms[0].equivalent(bond.atoms[1]):
				# An O-O bond is considered to be an "optical isomer" and so no
				# symmetry correction will be applied
				if bond.atoms[0].atomType.label == 'Os' and bond.atoms[1].atomType.label == 'Os' and \
					bond.atoms[0].getFreeElectronCount() == 0 and bond.atoms[1].getFreeElectronCount() == 0:
					pass
				else:
					structure = Structure(self.atoms(), self.bonds())
					structure.removeBond(bond)
					structure1, structure2 = structure.split()


					if bond.atoms[0] in structure1.atoms(): structure1.removeAtom(bond.atoms[0])
					if bond.atoms[1] in structure1.atoms(): structure1.removeAtom(bond.atoms[1])
					if bond.atoms[0] in structure2.atoms(): structure2.removeAtom(bond.atoms[0])
					if bond.atoms[1] in structure2.atoms(): structure2.removeAtom(bond.atoms[1])
					groups1 = structure1.split()
					groups2 = structure2.split()

					# Test functional groups for symmetry
					if len(groups1) == len(groups2) == 1:
						if groups1[0].isIsomorphic(groups2[0]): symmetryNumber *= 2
					elif len(groups1) == len(groups2) == 2:
						if groups1[0].isIsomorphic(groups2[0]) and groups1[1].isIsomorphic(groups2[1]): symmetryNumber *= 2
						elif groups1[1].isIsomorphic(groups2[0]) and groups1[0].isIsomorphic(groups2[1]): symmetryNumber *= 2
					elif len(groups1) == len(groups2) == 3:
						if groups1[0].isIsomorphic(groups2[0]) and groups1[1].isIsomorphic(groups2[1]) and groups1[2].isIsomorphic(groups2[2]): symmetryNumber *= 2
						elif groups1[0].isIsomorphic(groups2[0]) and groups1[1].isIsomorphic(groups2[2]) and groups1[2].isIsomorphic(groups2[1]): symmetryNumber *= 2
						elif groups1[0].isIsomorphic(groups2[1]) and groups1[1].isIsomorphic(groups2[2]) and groups1[2].isIsomorphic(groups2[0]): symmetryNumber *= 2
						elif groups1[0].isIsomorphic(groups2[1]) and groups1[1].isIsomorphic(groups2[0]) and groups1[2].isIsomorphic(groups2[2]): symmetryNumber *= 2
						elif groups1[0].isIsomorphic(groups2[2]) and groups1[1].isIsomorphic(groups2[0]) and groups1[2].isIsomorphic(groups2[1]): symmetryNumber *= 2
						elif groups1[0].isIsomorphic(groups2[2]) and groups1[1].isIsomorphic(groups2[1]) and groups1[2].isIsomorphic(groups2[0]): symmetryNumber *= 2

		return symmetryNumber

	def calculateAxisSymmetryNumber(self):
		"""
		Get the axis symmetry number correction. The "axis" refers to a series
		of two or more cumulated double bonds (e.g. C=C=C, etc.). Corrections
		for single C=C bonds are handled in getBondSymmetryNumber().
		
		Each axis (C=C=C) has the potential to double the symmetry number.
		If an end has 0 or 1 groups (eg. =C=CJJ or =C=C-R) then it cannot 
		alter the axis symmetry and is disregarded::
		
			A=C=C=C..        A-C=C=C=C-A
			
			  s=1                s=1
		
		If an end has 2 groups that are different then it breaks the symmetry 
		and the symmetry for that axis is 1, no matter what's at the other end::
		
			A\               A\         /A
			  T=C=C=C=C-A      T=C=C=C=T
			B/               A/         \B
			      s=1             s=1
		
		If you have one or more ends with 2 groups, and neither end breaks the 
		symmetry, then you have an axis symmetry number of 2::
		
			A\         /B      A\         
			  C=C=C=C=C          C=C=C=C-B
			A/         \B      A/         
			      s=2                s=2
		"""

		symmetryNumber = 1

		# List all double bonds in the structure
		doubleBonds = []
		for bond in self.bonds():
			if bond.isDouble(): doubleBonds.append(bond)

		# Search for adjacent double bonds
		cumulatedBonds = []
		for i, bond1 in enumerate(doubleBonds):
			for bond2 in doubleBonds[i+1:]:
				if bond1.atoms[0] in bond2.atoms or bond1.atoms[1] in bond2.atoms:
					listToAddTo = None
					for cumBonds in cumulatedBonds:
						if bond1 in cumBonds or bond2 in cumBonds:
							listToAddTo = cumBonds
					if listToAddTo:
						if bond1 not in listToAddTo: listToAddTo.append(bond1)
						if bond2 not in listToAddTo: listToAddTo.append(bond2)
					else:
						cumulatedBonds.append([bond1, bond2])

		# For each set of adjacent double bonds, check for axis symmetry
		for bonds in cumulatedBonds:
			
			# Do nothing if less than two cumulated bonds
			if len(bonds) < 2: continue

			# Do nothing if axis is in cycle
			#found = False
			#for bond in bonds:
			#	if self.isBondInCycle(bond): found = True
			#if found: continue

			# Find terminal atoms in axis
			# Terminal atoms labelled T:  T=C=C=C=T
			terminalAtoms = []
			for bond12 in bonds:
				atom1, atom2 = bond12.atoms
				# if atom is also in one of the OTHER bonds in the axis, 
				# then it's not a terminal
				for bond in bonds:
					if bond is not bond12:
						if atom1 in bond.atoms: atom1 = None
						if atom2 in bond.atoms: atom2 = None
				if atom1 is not None: terminalAtoms.append(atom1)	
				if atom2 is not None: terminalAtoms.append(atom2)
			if len(terminalAtoms) != 2: continue # when does this occur??
			
			# Remove axis from (copy of) structure
			structure = Structure(self.atoms(), self.bonds())
			for bond in bonds:
				structure.removeBond(bond)
			atomsToRemove = []
			for atom in structure.atoms():
				if len(structure.graph[atom]) == 0: # it's not bonded to anything
					atomsToRemove.append(atom)
			for atom in atomsToRemove: structure.removeAtom(atom)

			# Split remaining fragments of structure
			end_fragments = structure.split()
			# you may have only one end fragment,
			# eg. if you started with H2C=C=C..
			
			# 
			# there can be two groups at each end     A\         /B
			#                                           T=C=C=C=T
			#                                         A/         \B
			
			# to start with nothing has broken symmetry about the axis
			symmetry_broken=False 
			for fragment in end_fragments: # a fragment is one end of the axis
				
				# remove the atom that was at the end of the axis and split what's left into groups
				for atom in terminalAtoms:
					if atom in fragment.atoms(): fragment.removeAtom(atom)
				groups = fragment.split()
				
				# If end has only one group then it can't contribute to (nor break) axial symmetry
				#   Eg. this has no axis symmetry:   A-T=C=C=C=T-A
				# so we remove this end from the list of interesting end fragments
				if len(groups)==1:
					end_fragments.remove(fragment)
					continue # next end fragment
				if len(groups)==2:
					if not groups[0].isIsomorphic(groups[1]):
						# this end has broken the symmetry of the axis
						symmetry_broken = True
						
			# If there are end fragments left that can contribute to symmetry,
			# and none of them broke it, then double the symmetry number
			if end_fragments and not symmetry_broken:
				symmetryNumber *= 2
					
		return symmetryNumber

	def calculateCyclicSymmetryNumber(self):
		"""
		Get the symmetry number correction for cyclic regions of a molecule.
		For complicated fused rings the smallest set of smallest rings is used.
		"""

		symmetryNumber = 1

		# Get symmetry number for each ring in structure
		rings = self.getSmallestSetOfSmallestRings()
		for ring in rings:

			# Make copy of structure
			structure = Structure(self.atoms(), self.bonds())

			# Remove bonds of ring from structure
			for i, atom1 in enumerate(ring):
				for atom2 in ring[i+1:]:
					if structure.hasBond(atom1, atom2):
						bond = self.getBond(atom1, atom2)
						structure.removeBond(bond)

			structures = structure.split()
			groups = []
			for struct in structures:
				for atom in ring:
					if atom in struct.atoms(): struct.removeAtom(atom)
				groups.append(struct.split())

			# Find equivalent functional groups on ring
			equivalentGroups = []
			for group in groups:
				found = False
				for eqGroup in equivalentGroups:
					if not found:
						if group.isIsomorphic(eqGroup[0]):
							eqGroup.append(group)
							found = True
				if not found:
					equivalentGroups.append([group])

			# Find equivalent bonds on ring
			equivalentBonds = []
			for i, atom1 in enumerate(ring):
				for atom2 in ring[i+1:]:
					if self.hasBond(atom1, atom2):
						bond = self.getBond(atom1, atom2)
						found = False
						for eqBond in equivalentBonds:
							if not found:
								if bond.equivalent(eqBond[0]):
									eqBond.append(group)
									found = True
						if not found:
							equivalentBonds.append([bond])

			# Find maximum number of equivalent groups and bonds
			maxEquivalentGroups = 0
			for groups in equivalentGroups:
				if len(groups) > maxEquivalentGroups:
					maxEquivalentGroups = len(groups)
			maxEquivalentBonds = 0
			for bonds in equivalentBonds:
				if len(bonds) > maxEquivalentBonds:
					maxEquivalentBonds = len(bonds)

			if maxEquivalentGroups == maxEquivalentBonds == len(ring):
				symmetryNumber *= len(ring)
			else:
				symmetryNumber *= max(maxEquivalentGroups, maxEquivalentBonds)

			print len(ring), maxEquivalentGroups, maxEquivalentBonds, symmetryNumber


		return symmetryNumber

	def calculateSymmetryNumber(self):
		"""
		Return the symmetry number for the structure. The symmetry number
		includes both external and internal modes.
		"""
		symmetryNumber = 1

		for atom in self.atoms():
			if not self.isAtomInCycle(atom):
				symmetryNumber *= self.calculateAtomSymmetryNumber(atom)

		for bond in self.bonds():
			if not self.isBondInCycle(bond):
				symmetryNumber *= self.calculateBondSymmetryNumber(bond)

		symmetryNumber *= self.calculateAxisSymmetryNumber()

		#if self.isCyclic():
		#	symmetryNumber *= self.calculateCyclicSymmetryNumber()

		self.symmetryNumber = symmetryNumber

################################################################################

if __name__ == '__main__':

	structure1 = Structure()
	atom1 = structure1.addAtom(chem.Atom('C', '0'))
	atom2 = structure1.addAtom(chem.Atom('C', '0'))
	atom3 = structure1.addAtom(chem.Atom('C', '0'))
	atom4 = structure1.addAtom(chem.Atom('C', '0'))
	atom5 = structure1.addAtom(chem.Atom('C', '0'))
	atom6 = structure1.addAtom(chem.Atom('C', '0'))
	bond1 = structure1.addBond(chem.Bond([atom1, atom2], 'S'))
	bond2 = structure1.addBond(chem.Bond([atom2, atom3], 'S'))
	bond3 = structure1.addBond(chem.Bond([atom3, atom4], 'S'))
	bond4 = structure1.addBond(chem.Bond([atom4, atom5], 'S'))
	bond5 = structure1.addBond(chem.Bond([atom5, atom6], 'S'))

	structure2 = Structure()
	atom1 = structure2.addAtom(chem.Atom('C', '0'))
	atom2 = structure2.addAtom(chem.Atom('C', '0'))
	atom3 = structure2.addAtom(chem.Atom('C', '0'))
	atom4 = structure2.addAtom(chem.Atom('C', '0'))
	bond1 = structure2.addBond(chem.Bond([atom1, atom2], 'S'))
	bond2 = structure2.addBond(chem.Bond([atom2, atom3], 'S'))
	bond3 = structure2.addBond(chem.Bond([atom3, atom4], 'S'))

	for i in range(10000):
		match, map21List, map12List = structure1.findSubgraphIsomorphisms(structure2)
