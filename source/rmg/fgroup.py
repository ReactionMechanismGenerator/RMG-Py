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
Contains classes describing entities relating to functional groups.
"""

import quantities as pq
from chem import ElectronState, electronStates, BondType, bondTypes, ChemGraph

################################################################################

class Element:
    """
	Represent a single functional group element. Each element has a unique
	string `label`, a string `type` representing the symbol of the underlying
	chemical element, and a string `description` of the element.
	
	This class is specifically for properties that all atoms of the same 
	element share. Ideally there is only one instance of this class for each 
	element.
	"""

    def __init__(self, label, type, description, valence):
        """Initialize a functional group element.

        Parameters:
        label -- The unique label for the element
        type -- The symbol of the underlying chemical element
        description -- A description of the element
        """
        self.label = label
        self.type = type
        self.description = description
        self.valence = valence
        
################################################################################

def loadElements():
	"""
	Loads entries into a dictionary of elements. The dictionary created by
	this function is always available at :data:`rmg.fgroup.elements`.
	"""
	elements = {}
	elements['R'] = Element('R', 'R', 'generic functional group', 0)
	elements['R!H'] = Element('R!H', 'R!H', 'generic non-hydrogen functional group', 0)
	elements['H'] = Element('H', 'H', 'hydrogen atom', 1)
	elements['C'] = Element('C', 'C', 'carbon atom with undefined bonds', 4)
	elements['Ct'] = Element('Ct', 'C', 'carbon atom with one triple bond and one single bond', 4)
	elements['Cs'] = Element('Cs', 'C', 'carbon atom with four single bonds', 4)
	elements['Cd'] = Element('Cd', 'C', 'carbon atom with one double bond and the rest undefined', 4)
	elements['Cdd'] = Element('Cdd', 'C', 'carbon atom with two double bonds', 4)
	elements['Cds'] = Element('Cds', 'C', 'carbon atom with one double bond and two single bonds', 4)
	elements['Cb'] = Element('Cb', 'C', 'carbon atom belonging to a benzene ring', 4)
	elements['Cbf'] = Element('Cbf', 'C', 'carbon atom belonging to a fused benzene ring', 4)
	elements['O'] = Element('O', 'O', 'oxygen atom with undefined bonds', 2)
	elements['Os'] = Element('Os', 'O', 'oxygen atom with two single bonds', 2)
	elements['Od'] = Element('Od', 'O', 'oxygen atom with one double bond', 2)
	elements['CO'] = Element('CO', 'C', 'non-central carbon atom bonded with a double bond to a non-central O', 2)
	return elements
	
# The dictionary of elements, accessible by atomic number, symbol, or name.
elements = loadElements()

################################################################################

class Atom:
	"""
	Represent an atom in a functional group. Each atom is defined by an 
	`element` (stored internally as an :class:`Element` object), a `electronState`
	(stored internally as an :class:`ElectronState` object), a numeric `charge`,
	and a boolean `center` which is :data:`True` if the atom is the center atom
	for the functional group and :data:`False` if not.
	"""	

	def __init__(self, element=None, electronState=None, charge=0, center=False):
		"""
		Initialize an atom object.
		"""
		self.setElement(element)
		self.setElectronState(electronState)
		self.charge = charge
		self.center = center
		
	def setElement(self, element):
		"""
		Set the element that this atom represents. The *element* parameter can
		be an :class:`Element` object or a string containing the element label
		In all cases *element* will be converted to and stored as an 
		:class:`Element` object.
		"""
		if element.__class__ == str or element.__class__ == unicode:
			element = elements[element]
		elif element.__class__ == list:
			labels = element
			element = []
			for label in labels:
				element.append(elements[label])
			if len(element) == 1:
				element = element[0]
		self.element = element

	def setElectronState(self, electronState):
		"""
		Set the electron state that this atom represents. The *electronState* 
		parameter can be an :class:`ElectronState` object or a string 
		representing the label of the desired electron state. In all cases 
		*electronState*	will be converted to and stored as an 
		:class:`ElectronState` object.
		"""
		if electronState.__class__ == str or electronState.__class__ == unicode:
			electronState = electronStates[electronState]
		elif electronState.__class__ == list:
			labels = electronState
			electronState = []
			for label in labels:
				electronState.append(electronStates[label])
			if len(electronState) == 1:
				electronState = electronState[0]
		self.electronState = electronState

################################################################################

class Bond:
	"""
	Represent a bond in a functional group. Each bond has a list `atoms` of 
	length two containing the two atoms in the bond and a `bondType` object,
	stored internally as a :class:`BondType` object.
	"""	

	def __init__(self, atoms, bondType=''):
		self.setBondType(bondType)
		self.atoms = atoms

	def setBondType(self, bondType):
		"""
		Set the bond type that this bond represents. The *bondType* 
		parameter can be a :class:`BondType` object, a number representing
		the bond order, or a string representing the label of the desired bond 
		type. In all cases *bondType* will be converted to and stored as a
		:class:`BondType` object.
		"""
		if bondType.__class__ == str or bondType.__class__ == unicode:
			bondType = bondTypes[bondType]
		elif bondType.__class__ == list:
			labels = bondType
			bondType = []
			for label in labels:
				bondType.append(bondTypes[label])
			if len(bondType) == 1:
				bondType = bondType[0]
		self.bondType = bondType

################################################################################

class Structure(ChemGraph):
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
		
	def fromAdjacencyList(self, adjlist):
		"""
		Convert a string adjacency list `adjlist` into a functional group 
		object.
		"""
		
		atoms = []; bonds = []; atomdict = {}
				
		lines = adjlist.splitlines()
		
		for line in lines[1:]:
			
			data = line.split()
			
			# First item is index for atom
			aid = int(data[0])
			
			# If second item is '*', the atom is the center atom
			center = False
			index = 1
			if data[1] == '*':
				center = True
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
			atom = Atom(elements, elecState, 0, center)
				
			# Add the atom to the list
			atoms.append(atom)
			atomdict[aid] = atom
			
			# Process list of bonds
			for datum in data[index+2:]:
				aid2, comma, btype = datum[1:-1].partition(',')
				aid2 = int(aid2)
				
				if btype[0] == '{':
					btype = btype[1:-1].split(',')
				else:
					btype = [btype]
				
				if aid2 in atomdict:
					bond = Bond([atomdict[aid], atomdict[aid2]], btype)
					bonds.append(bond)
			
		# Create and return functional group or species
		self.initialize(atoms, bonds)
	
################################################################################

class FunctionalGroup:
	"""
	Represent a functional group. Each functional group has a unique string 
	`label`. The `structure` variable contains a :class:`Structure` object 
	representing the functional group connectivity.
	"""	

	def __init__(self, label='', structure=None):
		"""
		Initialize a functional group.
		"""
		self.label = label
		if structure is None:
			self.structure = Structure()
		else:
			self.structure = structure
		
################################################################################
