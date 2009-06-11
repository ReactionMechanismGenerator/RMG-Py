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
Contains classes describing chemical entities: elements, atoms, bonds, species, etc.
"""

import logging
import quantities as pq
import pybel
import openbabel

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

class Element:
    """
	Represent a single chemical element. Each element has an atomic
	`number`, a `name`, a `symbol`, an atomic `mass`, and a `valence`, a list
	of the possible numbers of bonds allowed.

	This class is specifically for properties that all atoms of the same element
	share. Ideally there is only one instance of this class for each element.
	"""

    def __init__(self, number, name, symbol, mass, valence):
        """Initialize a chemical element.

        Parameters:
        number -- The atomic number of the element
        name -- The name of the element
        symbol -- The symbol for the element
        mass -- The atomic mass of the element
        valence -- The valency of the element
        """
        self.number = number
        self.name = name
        self.symbol = symbol
        self.mass = mass
        self.valence = valence

################################################################################

def loadElements():
	"""
	Loads entries into a dictionary of elements.
	"""

	# Chemical elements
	elements = {}
	elements[1] = elements['H'] = elements['hydrogen'] = Element(1, 'hydrogen', 'H', pq.Quantity(1.00794, 'g/mol'), 1)
	elements[2] = elements['He'] = elements['helium'] = Element(2, 'helium', 'He', pq.Quantity(4.002602, 'g/mol'), 0)
	elements[6] = elements['C'] = elements['carbon'] = Element(6, 'carbon', 'C', pq.Quantity(12.0107, 'g/mol'), 4)
	elements[7] = elements['N'] = elements['nitrogen'] = Element(7, 'nitrogen', 'N', pq.Quantity(14.00674, 'g/mol'), [3,5])
	elements[8] = elements['O'] = elements['oxygen'] = Element(8, 'oxygen', 'O', pq.Quantity(15.9994, 'g/mol'), 2)
	elements[9] = elements['F'] = elements['fluorine'] = Element(9, 'fluorine', 'F', pq.Quantity(18.998403, 'g/mol'), 1)
	elements[10] = elements['Ne'] = elements['neon'] = Element(10, 'neon', 'Ne', pq.Quantity(20.1797, 'g/mol'), 0)
	elements[14] = elements['Si'] = elements['silicon'] = Element(14, 'silicon', 'Si', pq.Quantity(28.0855, 'g/mol'), 4)
	elements[15] = elements['P'] = elements['phosphorus'] = Element(15, 'phosphorus', 'P', pq.Quantity(30.973761, 'g/mol'), [3,5])
	elements[16] = elements['S'] = elements['sulfur'] = Element(16, 'sulfur', 'S', pq.Quantity(32.065, 'g/mol'), [2,6])
	elements[17] = elements['Cl'] = elements['chlorine'] = Element(17, 'chlorine', 'Cl', pq.Quantity(35.453, 'g/mol'), 1)
	elements[18] = elements['Ar'] = elements['argon'] = Element(18, 'argon', 'Ar', pq.Quantity(39.348, 'g/mol'), 0)
	elements[35] = elements['Br'] = elements['bromine'] = Element(35, 'bromine', 'Br', pq.Quantity(79.904, 'g/mol'), 1)
	elements[53] = elements['I'] = elements['iodine'] = Element(53, 'iodine', 'I', pq.Quantity(126.90447, 'g/mol'), 1)
	elements[0] = elements['R'] = Element(0, '', 'R', pq.Quantity(0.0, 'g/mol'), 0)

	return elements

# The dictionary of elements, accessible by atomic number, symbol, or name.
elements = loadElements()

################################################################################

class AtomType:
	"""
	Represent a single atom type by its chemical element and, optionally, some
	information about the local bond structure around that element. Each
	element has a unique string `label`, the underlying chemical `element`, and
	a string `description` of the element.

	This class is specifically for properties that all atoms of the same
	element share. Ideally there is only one instance of this class for each
	element.
	"""

	def __init__(self, label, element, description):
		"""
		Initialize an atom type.
		"""
		self.label = label
		self.element = element
		self.description = description

################################################################################

def loadAtomTypes():
	"""
	Loads entries into a dictionary of atom types.
	"""

	# Functional group atom types
	atomTypes = {}
	atomTypes['H']		= AtomType('H', 	elements['H'], 	'hydrogen')
	atomTypes['C']		= AtomType('C', 	elements['C'], 	'carbon')
	atomTypes['N']		= AtomType('N', 	elements['N'], 	'nitrogen')
	atomTypes['O']		= AtomType('O', 	elements['O'], 	'oxygen')
	atomTypes['F']		= AtomType('F', 	elements['F'], 	'fluorine')
	atomTypes['Ne'] 	= AtomType('Ne', 	elements['Ne'], 'neon')
	atomTypes['Si'] 	= AtomType('Si', 	elements['Si'], 'silicon')
	atomTypes['P']		= AtomType('P', 	elements['P'],	'phosphorus')
	atomTypes['S']		= AtomType('S', 	elements['S'],	'sulfur')
	atomTypes['Cl'] 	= AtomType('Cl', 	elements['Cl'], 'chlorine')
	atomTypes['Ar'] 	= AtomType('Ar', 	elements['Ar'], 'argon')
	atomTypes['Br'] 	= AtomType('Br', 	elements['Br'], 'bromine')
	atomTypes['I']		= AtomType('I', 	elements['I'],	'iodine')
	atomTypes['R']		= AtomType('R', 	None, 			'generic functional group')
	atomTypes['R!H'] 	= AtomType('R!H',	None, 			'generic non-hydrogen functional group')
	atomTypes['Ct'] 	= AtomType('Ct', 	elements['C'], 	'carbon with one triple bond and one single bond')
	atomTypes['Cs'] 	= AtomType('Cs', 	elements['C'], 	'carbon with four single bonds')
	atomTypes['Cd'] 	= AtomType('Cd', 	elements['C'], 	'carbon with one double bond and the rest undefined')
	atomTypes['Cdd'] 	= AtomType('Cdd', 	elements['C'], 	'carbon with two double bonds')
	atomTypes['Cds'] 	= AtomType('Cds', 	elements['C'], 	'carbon with one double bond and two single bonds')
	atomTypes['Cb'] 	= AtomType('Cb', 	elements['C'], 	'carbon belonging to a benzene ring')
	atomTypes['Cbf'] 	= AtomType('Cbf', 	elements['C'], 	'carbon belonging to a fused benzene ring')
	atomTypes['Os'] 	= AtomType('Os', 	elements['O'], 	'oxygen with two single bonds')
	atomTypes['Od'] 	= AtomType('Od', 	elements['O'], 	'oxygen with one double bond')
	atomTypes['CO'] 	= AtomType('CO', 	elements['C'], 	'non-central carbon bonded with a double bond to a non-central O')

	return atomTypes

# The dictionary of elements, accessible by atomic number, symbol, or name.
atomTypes = loadAtomTypes()

################################################################################

class ElectronState:
	"""
	Represent a single free electron state (none, radical, etc.) Each state is
	defined by a unique string `label`; the `order`, or number of
	free electrons; and a `spin` multiplicity.

	This class is specifically for properties that all free electron states
	share. Ideally there is only one instance of this class for each free
	electron state.
	"""

	def __init__(self, label, order, spin=''):
		"""
		Initialize a free electron state.

		Parameters:
		label -- A string unique to this free electron state
		order -- The number of free electrons
		spin -- The spin state for polyradicals (singlet = '1', triplet = '3')
		"""
		self.label = label
		self.order = order
		self.spin = spin

############################################################################

def loadElectronStates():
	"""
	Loads entries into a dictionary of free electron states. The dictionary
	created by this function is always available at
	:data:`rmg.chem.electronStates`.
	"""
	electronStates = {}
	electronStates['0'] = ElectronState('0', 0, 1)
	electronStates['1'] = ElectronState('1', 1, 2)
	electronStates['2'] = ElectronState('2', 2, [1,3])
	electronStates['2S'] = ElectronState('2S', 2, 1)
	electronStates['2T'] = ElectronState('2T', 2, 3)
	electronStates['3'] = ElectronState('3', 3, [2,4])
	electronStates['4'] = ElectronState('4', 4, [1,3,5])
	return electronStates

# The dictionary of electron states, accessible by label.
electronStates = loadElectronStates()

################################################################################

class BondType:
	"""
	Represent a type of chemical bond. Each bond type has a unique string
	`label`; a unique string `name`; a numeric bond `order`; an integral
	`piElectrons`, the number of pi electrons; and a string `location` with
	bond geometry information (i.e. 'cis' or 'trans').

	This class is specifically for properties that all bonds of the same type
	share. Ideally there is only one instance of this class for each bond type.
	"""

	def __init__(self, label, name, order, piElectrons, location=''):
		"""Initialize a chemical bond.

		Parameters:
		label -- A short string unique to this bond type
		name -- A longer string unique to this bond type
		order -- The bond order (1 = single, 2 = double, 3 = triple, etc.)
		piElectrons -- The number of pi electrons in the bond
		location -- Bond location information ('cis' or 'trans')
		"""
		self.label = label
		self.name = name
		self.order = order
		self.piElectrons = piElectrons
		self.location = location

	def __repr__(self):
		"""x.__repr__() <==> repr(x)"""
		return self.__module__+".BondType('%s','%s',%g,%g,location='%s')"%(self.label,self.name,self.order,self.piElectrons,self.location)

################################################################################

def loadBondTypes():
	"""
	Loads entries into a dictionary of bond types. The dictionary created by
	this function is always available at :data:`rmg.chem.bondTypes`.
	"""

	bondTypes = {}
	bondTypes[1] = bondTypes['S'] = bondTypes['single'] = BondType('S', 'single', 1, 0, '')
	bondTypes[2] = bondTypes['D'] = bondTypes['double'] = BondType('D', 'double', 2, 2, '')
	bondTypes['Dcis'] = BondType('Dcis', 'double_cis', 2, 2, 'cis')
	bondTypes['Dtrans'] = BondType('Dtrans', 'double_trans', 2, 2, 'trans')
	bondTypes[3] = bondTypes['T'] = bondTypes['triple'] = BondType('T', 'triple', 3, 4, '')
	bondTypes[1.5] = bondTypes['B'] = bondTypes['benzene'] = BondType('B', 'benzene', 1.5, 1, '')
	return bondTypes

# The dictionary of bond types, accessible by label, order, or name.
bondTypes = loadBondTypes()

################################################################################

class Atom(object):
	"""
	Represent an atom in a chemical species or functional group.
	"""

	def __init__(self, atomType=None, electronState=None, charge=0, label=''):
		"""
		Initialize an atom object.
		"""
		self.atomType = atomType
		self.electronState = electronState
		self.charge = charge
		self.label = label
		
	def getAtomType(self):
		if len(self._atomType) == 1: return self._atomType[0]
		else: return self._atomType

	def setAtomType(self, atomType):
		"""
		Set the element that this atom represents. The `electronState`
		parameter is an :class:`FGElement` object or a string representing the
		label of the desired element. If the parameter is a string, it will be
		converted to an	:class:`FGElement` object.
		"""
		self._atomType = []
		if atomType.__class__ != list:
			atomType = [atomType]
		for atom in atomType:
			if (atom.__class__ == str or atom.__class__ == unicode) and atom in atomTypes:
				atom = atomTypes[atom]
			if atom.__class__ == AtomType:
				self._atomType.append(atom)
			else:
				raise Exception('Invalid parameter "atomType".')
		
	atomType = property(getAtomType, setAtomType)

	def getElectronState(self):
		if len(self._electronState) == 1: return self._electronState[0]
		else: return self._electronState

	def setElectronState(self, electronState):
		"""
		Set the electron state that this atom represents. The *electronState*
		parameter can be an :class:`ElectronState` object or a string
		representing the label of the desired electron state. In all cases
		*electronState*	will be converted to and stored as an
		:class:`ElectronState` object.
		"""
		if electronState.__class__ != list:
			electronState = [electronState]
		self._electronState = []
		for state in electronState:
			if (state.__class__ == str or state.__class__ == unicode) and state in electronStates:
				state = electronStates[state]
			if state.__class__ == ElectronState:
				self._electronState.append(state)
			else:
				raise Exception('Invalid parameter "electronState".')

	electronState = property(getElectronState, setElectronState)

	def equivalent(self, other):
		"""
		Return :data:`True` if `self` and `other` are equivalent atoms, i.e.
		they have the same element and electronic state.
		"""

		atomTypesMatch = False; electronStatesMatch = False

		for atomType1 in self._atomType:
			for atomType2 in other._atomType:

				# If either is a generic atom type, then always return True
				if atomType1.label == 'R' or atomType2.label == 'R':
					atomTypesMatch = True
				# If either is a generic non-hydrogen atom type, then return
				# True if any atom type in the remaining one is non-hydrogen
				elif atomType1.label == 'R!H':
					if atomType2.label != 'H': atomTypesMatch = True
				elif atomType2.label == 'R!H':
					if atomType1.label != 'H': atomTypesMatch = True
				# If either represents an element without surrounding bond info,
				# match remaining to any with the same element
				elif atomType1.label == atomType1.element.symbol == \
						atomType2.element.symbol:
					atomTypesMatch = True
				elif atomType2.label == atomType2.element.symbol == \
						atomType1.element.symbol:
					atomTypesMatch = True
				# Special case: 'Cd' matches any of 'Cd', 'Cdd', or 'Cds'
				elif atomType1.label == 'Cd' and (atomType2.label == 'Cd' or \
						atomType2.label == 'Cdd' or atomType2.label == 'Cds'):
					atomTypesMatch = True
				elif atomType2.label == 'Cd' and (atomType1.label == 'Cd' or \
						atomType1.label == 'Cdd' or atomType1.label == 'Cds'):
					atomTypesMatch = True
				# Otherwise labels must match exactly
				elif atomType1.label == atomType2.label:
					atomTypesMatch = True

		for elecState1 in self._electronState:
			for elecState2 in other._electronState:

				if elecState1.label == '2' and (elecState2.label == '2' or elecState2.label == '2S' or elecState2.label == '2T'):
					electronStatesMatch = True
				elif (elecState1.label == '2' or elecState1.label == '2S' or elecState1.label == '2T') and elecState2.label == '2':
					electronStatesMatch = True
				elif elecState1.label == elecState2.label:
					electronStatesMatch = True

		return (atomTypesMatch and electronStatesMatch)
	
	def copy(self):
		"""
		Generate a copy of the current atom.
		"""
		return Atom(self.atomType, self.electronState, self.charge, self.label)

	def isCenter(self):
		"""
		Return :data:`True` if the atom is a center atom and :data:`False`
		otherwise.
		"""
		return len(self.label) > 0

	def isElement(self, symbol):
		"""
		Return :data:`True` if the atom has the element with the symbol
		`symbol` and :data:`False` otherwise.
		"""
		if len(self._atomType) == 0 or len(self._atomType) > 1: return False
		elif self.atomType.element is None: return False
		else: return self.atomType.element.symbol == symbol

	def isHydrogen(self):
		"""
		Return :data:`True` if the atom is a hydrogen atom and :data:`False`
		otherwise.
		"""
		return self.isElement('H')

	def isCarbon(self):
		"""
		Return :data:`True` if the atom is a carbon atom and :data:`False`
		otherwise.
		"""
		return self.isElement('C')

	def isOxygen(self):
		"""
		Return :data:`True` if the atom is an oxygen atom and :data:`False`
		otherwise.
		"""
		return self.isElement('O')

	def isNonHydrogen(self):
		"""
		Return :data:`True` if the atom is not a hydrogen atom and :data:`False`
		otherwise.
		"""
		return self.atomType.element.symbol != 'H'

	def hasFreeElectron(self):
		"""
		Return :data:`True` if the atom has one or more unpaired electrons and
		:data:`False` otherwise.
		"""
		return self.electronState.order > 0

	def freeElectronCount(self):
		"""
		Return the number of unpaired electrons.
		"""
		return self.electronState.order

	def increaseFreeElectron(self):
		"""
		Increase the number of unpaired electrons on this atom by one.
		"""
		if self.electronState.label == '0': self.electronState = electronStates['1']
		elif self.electronState.label == '1': self.electronState = electronStates['2']
		elif self.electronState.label == '2': self.electronState = electronStates['3']
		elif self.electronState.label == '2S': self.electronState = electronStates['3']
		elif self.electronState.label == '2T': self.electronState = electronStates['3']
		else:
			logging.exception('Cannot increase the radical number of this atom.')

	def decreaseFreeElectron(self):
		"""
		Decrease the number of unpaird electrons on this atom by one.
		"""
		if self.electronState.label == '1': self.electronState = electronStates['0']
		elif self.electronState.label == '2': self.electronState = electronStates['1']
		elif self.electronState.label == '2S': self.electronState = electronStates['1']
		elif self.electronState.label == '2T': self.electronState = electronStates['1']
		elif self.electronState.label == '3': self.electronState = electronStates['2']
		else:
			logging.exception('Cannot decrease the radical number of this atom.')

################################################################################

class Bond:
	"""
	Represent a bond in a chemical species. Each bond has a list `atoms` of
	length two containing the two atoms in the bond and a `bondType` object,
	stored internally as a :class:`BondType` object.
	"""

	def __init__(self, atoms, bondType=''):
		self.setBondType(bondType)
		self.atoms = atoms

	def __repr__(self):
		"""x.__repr__() <==> repr(x)"""
		return self.__module__+".Bond(%s,bondType='%s')"%(repr(self.atoms),self.bondType.label)

	def getBondType(self):
		if len(self._bondType) == 1: return self._bondType[0]
		else: return self._bondType

	def setBondType(self, bondType):
		"""
		Set the bond type that this bond represents. The `bondType`
		parameter can be a :class:`BondType` object, a number representing
		the bond order, or a string representing the label of the desired bond
		type. In all cases `bondType` will be converted to and stored as a
		:class:`BondType` object.
		"""
		if bondType.__class__ != list:
			bondType = [bondType]
		self._bondType = []
		for bond in bondType:
			if (bond.__class__ == str or bond.__class__ == unicode) and bond in bondTypes:
				bond = bondTypes[bond]
			if bond.__class__ == BondType:
				self._bondType.append(bond)
			else:
				raise Exception('Invalid parameter "bondType".')

	bondType = property(getBondType, setBondType)

	def equivalent(self, other):
		"""
		Return :data:`True` if `self` and `other` are equivalent bonds, i.e.
		they have the same bond type.
		"""
		for bondType1 in self._bondType:
			for bondType2 in other._bondType:
				if bondType1 == bondType2:
					return True

		return False

	def copy(self):
		"""
		Generate a copy of the current bond.
		"""
		return Bond(self.atoms, self.bondType)

	def isSingle(self):
		"""
		Return :data:`True` if the bond represents a single bond and
		:data:`False` otherwise.
		"""
		return self.bondType.label == 'S'

	def isDouble(self):
		"""
		Return :data:`True` if the bond represents a double bond and
		:data:`False` otherwise.
		"""
		return self.bondType.label == 'D'

	def isTriple(self):
		"""
		Return :data:`True` if the bond represents a triple bond and
		:data:`False` otherwise.
		"""
		return self.bondType.label == 'T'

	def isBenzene(self):
		"""
		Return :data:`True` if the bond represents a benzene (aromatic) bond and
		:data:`False` otherwise.
		"""
		return self.bondType.label == 'B'

	def canIncreaseOrder(self):
		"""
		Return :data:`True` if the bond order can be increased by one.
		"""
		return self.bondType.label == 'S' or self.bondType.label == 'D'

	def canDecreaseOrder(self):
		"""
		Return :data:`True` if the bond order can be decreased by one without
		breaking.
		"""
		return self.bondType.label == 'D' or self.bondType.label == 'T'

	def increaseOrder(self):
		"""
		Increase the bond order by one.
		"""
		if self.bondType.label == 'S': self.bondType = bondTypes['D']
		elif self.bondType.label == 'D': self.bondType = bondTypes['T']
		else:
			logging.exception('Cannot increase the bond order of this bond.')

	def decreaseOrder(self):
		"""
		Decrease the bond order by one.
		"""
		if self.bondType.label == 'D': self.bondType = bondTypes['S']
		elif self.bondType.label == 'T': self.bondType = bondTypes['D']
		else:
			logging.exception('Cannot decrease the bond order of this bond.')


################################################################################

class Structure(graph.ChemGraph):
	"""
	A representation of a chemical species using a graph data structure. The
	vertices represent atoms, while the edges represent bonds.
	"""

	def __init__(self, atoms=[], bonds=[]):
		graph.ChemGraph.initialize(self, atoms, bonds)

	def fromAdjacencyList(self, adjlist):
		"""
		Convert a string adjacency list `adjlist` into a structure object.
		"""

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
			atom = Atom(elements, elecState, 0, center)

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
					bond = Bond([atomdict[aid], atomdict[aid2]], btype)
					bonds.append(bond)

				bonddict[aid][aid2] = btype

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
		self.updateAtomTypes()


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

			atom = Atom(elements[number].symbol, electronStates[electron])
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

					bond = Bond([atoms[i], atoms[j]], bondTypes[order])
					bonds.append(bond)

		# Create the graph from the atom and bond lists
		self.initialize(atoms, bonds)

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
			obmol.AddBond(index1+1, index2+1, order)

		return obmol

	def toCML(self):
		"""
		Convert a Structure object to CML.
		"""
		mol = pybel.Molecule(self.toOBMol())
		return mol.write('cml').strip()

	def toInChI(self):
		"""
		Convert a Structure object to an InChI string.
		"""
		mol = pybel.Molecule(self.toOBMol())
		return mol.write('inchi').strip()

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
		cml = dom.createElement('cml')
		root.appendChild(cml)

		dom2 = xml.dom.minidom.parseString(self.toCML())
		cml.appendChild(dom2.documentElement)

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
			for atom2, bond12 in self.graph[atom1].iteritems():
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
				elif double == 1 and (single == 1 or single == 2): atomType = 'Cds'
				elif double == 1: atomType = 'Cd'
				elif benzene == 1 or benzene == 2: atomType = 'Cb'
				elif benzene == 3: atomType = 'Cbf'
			elif atomType == 'O':
				if single == 1 or single == 2: atomType = 'Os'
				elif double == 1: atomType = 'Od'

			# Do nothing if suggested and specified atom types are identical
			if atom1.atomType.label == atomType:
				pass
			# Do nothing if suggested atom type is element
			elif atomType == atom1.atomType.element.symbol or atomType == 'Cd':
				pass
			# Do nothing if specified atom type is 'Cds' or 'Cdd' and suggested is 'Cd'
			elif (atom1.atomType.label == 'Cds' or atom1.atomType.label == 'Cdd') and atomType == 'Cd':
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
				atom1.atomType = atomTypes[atomType]
			# Make change if specified atom type is 'Cd' and suggested is 'Cds' or 'Cdd'
			elif atom1.atomType.label == 'Cd' and (atomType == 'Cds' or atomType == 'Cdd' or atomType == 'CO'):
				#logging.warning('Changed "' + atom1.atomType.label + '" to "' + atomType + '".')
				atom1.atomType = atomTypes[atomType]
			# Else print warning
			else:
				logging.warning('Suggested atom type "' + atomType + '" does not match specified atom type "' + atom1.atomType.label + '".')
			
	def radicalCount(self):
		"""
		Get the number of radicals in the structure.
		"""
		radical = 0
		for atom in self.atoms():
			radical += atom.electronState.order
		return radical

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
		for atom2, bond12 in self.graph[atom1].iteritems():
			# Vinyl bond must be capable of gaining an order
			if bond12.canIncreaseOrder():
				for atom3, bond23 in self.graph[atom2].iteritems():
					# Allyl bond must be capable of losing an order without
					# breaking
					if atom1 is not atom3 and bond23.canDecreaseOrder():
						paths.append([atom1, atom2, atom3, bond12, bond23])
		return paths

	def getCenter(self):
		"""
		Return the center atom of the functional group structure.
		"""
		for atom in self.atoms():
			if atom.isCenter(): return atom
		return None

################################################################################

if __name__ == '__main__':
	
	print ''
	
	print 'Atom types available:'
	for key, atomType in atomTypes.iteritems():
		print '\t' + str(key)
	print ''

	print 'Free electron states available:'
	for key, electronState in electronStates.iteritems():
		print '\t' + str(key)
	print ''
	
	print 'Bond types available:'
	for key, bondType in bondTypes.iteritems():
		print '\t' + str(key)
	print ''
	
