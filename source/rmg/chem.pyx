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
Contains classes describing simple chemical entities: elements, atoms, bonds, etc.
"""

import logging

cimport graph
#import graph

################################################################################

cdef class Element:
	"""
	Represent a single chemical element. Each element has an atomic
	`number`, a `name`, a `symbol`, an atomic `mass`, and a `valence`, a list
	of the possible numbers of bonds allowed.

	This class is specifically for properties that all atoms of the same element
	share. Ideally there is only one instance of this class for each element.
	"""

	cdef public int number
	cdef public str name
	cdef public str symbol
	cdef public float mass
	cdef public list valence

	def __init__(self, number, name, symbol, mass, valence):
		"""
		Initialize a chemical element.
		"""
		self.number = number
		self.name = name
		self.symbol = symbol
		self.mass = mass
		self.valence = valence

	def __reduce__(self):
		"""
		Used for pickling.
		"""
		return (Element, (self.number, self.name, self.symbol, self.mass, self.valence))

################################################################################

def loadElements():
	"""
	Loads entries into a dictionary of elements.
	"""

	# Chemical elements
	elements = {}
	elements[0] = elements['R'] = Element(0, '', 'R', 0.0, [0])
	elements[1] = elements['H'] = elements['hydrogen'] = Element(1, 'hydrogen', 'H', 0.00100794, [1])
	elements[2] = elements['He'] = elements['helium'] = Element(2, 'helium', 'He', 0.004002602, [0])
	elements[6] = elements['C'] = elements['carbon'] = Element(6, 'carbon', 'C', 0.0120107, [4])
	elements[7] = elements['N'] = elements['nitrogen'] = Element(7, 'nitrogen', 'N', 0.01400674, [3,5])
	elements[8] = elements['O'] = elements['oxygen'] = Element(8, 'oxygen', 'O', 0.0159994, [2])
	elements[9] = elements['F'] = elements['fluorine'] = Element(9, 'fluorine', 'F', 0.018998403, [1])
	elements[10] = elements['Ne'] = elements['neon'] = Element(10, 'neon', 'Ne', 0.0201797, [0])
	elements[14] = elements['Si'] = elements['silicon'] = Element(14, 'silicon', 'Si', 0.0280855, [4])
	elements[15] = elements['P'] = elements['phosphorus'] = Element(15, 'phosphorus', 'P', 0.030973761, [3,5])
	elements[16] = elements['S'] = elements['sulfur'] = Element(16, 'sulfur', 'S', 0.032065, [2,6])
	elements[17] = elements['Cl'] = elements['chlorine'] = Element(17, 'chlorine', 'Cl', 0.035453, [1])
	elements[18] = elements['Ar'] = elements['argon'] = Element(18, 'argon', 'Ar', 0.039348, [0])
	elements[35] = elements['Br'] = elements['bromine'] = Element(35, 'bromine', 'Br', 0.079904, [1])
	elements[53] = elements['I'] = elements['iodine'] = Element(53, 'iodine', 'I', 0.12690447, [1])
	
	return elements

# The dictionary of elements, accessible by atomic number, symbol, or name.
elements = loadElements()

################################################################################

cdef class AtomType:
	"""A type of atom such as Cd, a carbon atom with all single bonds.
	
	Represent a single atom type by its chemical element and, optionally, some
	information about the local bond structure around that element. Each
	element has a unique string `label`, the underlying chemical `element`, and
	a string `description` of the element.

	This class is specifically for properties that all atoms of the same
	element share. Ideally there is only one instance of this class for each
	element.
	"""

	cdef public str label
	cdef public object element
	cdef public str description

	def __init__(self, label, element, description):
		"""
		Initialize an atom type.
		"""
		self.label = label
		self.element = element
		self.description = description

	def __reduce__(self):
		"""
		Used for pickling.
		"""
		return (AtomType, (self.label, self.element, self.description))

	cpdef bint equivalent(AtomType self, AtomType other):
		"""
		Returns :data:`True` if two atom types are equivalent or :data:`False`
		otherwise. Respects wildcards, e.g. returns True for {R!H}=={C}
		"""

		# If either is a generic atom type, then always return True
		if self.label == 'R' or other.label == 'R':
			return True
		# If either is a generic non-hydrogen atom type, then return
		# True if any atom type in the remaining one is non-hydrogen
		elif self.label == 'R!H':
			if other.label != 'H':
				#logging.debug('I think %s == %s'%(self.label, other.label))
				return True
			else: 
				return False
		elif other.label == 'R!H':
			if self.label != 'H':
				#logging.debug('I think %s == %s'%(self.label, other.label))
				return True
			else: 
				return False
		# If either represents an element without surrounding bond info,
		# match remaining to any with the same element
		elif self.label == self.element.symbol == \
				other.element.symbol:
			return True
		elif other.label == other.element.symbol == \
				self.element.symbol:
			return True
		# Special case: 'Cd' matches any of 'Cd', 'Cdd', 'Cds', or 'CO'
		#elif self.label == 'Cd' and (other.label == 'Cd' or \
		#		other.label == 'Cdd' or other.label == 'Cds' or other.label == 'CO'):
		#	return True
		#elif other.label == 'Cd' and (self.label == 'Cd' or \
		#		self.label == 'Cdd' or self.label == 'Cds' or self.label == 'CO'):
		#	return True
		# Special case: 'Sid' matches any of 'Sid', 'Sidd', 'Sids', or 'SiO'
		#elif self.label == 'Sid' and (other.label == 'Sid' or \
		#		other.label == 'Sidd' or other.label == 'Sids' or other.label == 'SiO'):
		#	return True
		#elif other.label == 'Sid' and (self.label == 'Sid' or \
		#		self.label == 'Sidd' or self.label == 'Sids' or self.label == 'SiO'):
		#	return True
		# Otherwise labels must match exactly
		elif self.label == other.label:
			return True

		return False

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
	atomTypes['Cd'] 	= AtomType('Cd', 	elements['C'], 	'carbon with one double bond and two single bonds')
	atomTypes['Cdd'] 	= AtomType('Cdd', 	elements['C'], 	'carbon with two double bonds')
	atomTypes['Cb'] 	= AtomType('Cb', 	elements['C'], 	'carbon belonging to a benzene ring')
	atomTypes['Cbf'] 	= AtomType('Cbf', 	elements['C'], 	'carbon belonging to a fused benzene ring')
	atomTypes['CO'] 	= AtomType('CO', 	elements['C'], 	'non-central carbon bonded with a double bond to a non-central oxygen')
	atomTypes['Os'] 	= AtomType('Os', 	elements['O'], 	'oxygen with two single bonds')
	atomTypes['Od'] 	= AtomType('Od', 	elements['O'], 	'oxygen with one double bond')

	atomTypes['Sit'] 	= AtomType('Sit', 	elements['Si'], 'silicon with one triple bond and one single bond')
	atomTypes['Sis'] 	= AtomType('Sis', 	elements['Si'], 'silicon with four single bonds')
	atomTypes['Sid'] 	= AtomType('Sids', 	elements['Si'], 'silicon with one double bond and two single bonds')
	atomTypes['Sidd'] 	= AtomType('Sidd', 	elements['Si'], 'silicon with two double bonds')
	atomTypes['Sib'] 	= AtomType('Sib', 	elements['Si'], 'silicon belonging to a benzene ring')
	atomTypes['Sibf'] 	= AtomType('Sibf', 	elements['Si'], 'silicon belonging to a fused benzene ring')
	atomTypes['SiO'] 	= AtomType('SiO', 	elements['Si'], 'non-central silicon bonded with a double bond to a non-central oxygen')

	return atomTypes

# The dictionary of elements, accessible by atomic number, symbol, or name.
atomTypes = loadAtomTypes()

################################################################################

cdef class ElectronState:
	"""
	Represent a single free electron state (none, radical, etc.) Each state is
	defined by a unique string `label`; the `order`, or number of
	free electrons; and a `spin` multiplicity.

	This class is specifically for properties that all free electron states
	share. Ideally there is only one instance of this class for each free
	electron state.
	"""

	cdef public str label
	cdef public int order
	cdef public list spin

	def __init__(self, label, order, spin):
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

	def __reduce__(self):
		"""
		Used for pickling.
		"""
		return (ElectronState, (self.label, self.order, self.spin))

	cpdef bint equivalent(ElectronState self, ElectronState other):
		"""
		Returns :data:`True` if two electron states are equivalent or
		:data:`False` otherwise.
		"""

		if self.label == '2' and (other.label == '2' or other.label == '2S' or other.label == '2T'):
			return True
		elif (self.label == '2' or self.label == '2S' or self.label == '2T') and other.label == '2':
			return True
		elif self.label == other.label:
			return True
		
		return False

############################################################################

def loadElectronStates():
	"""
	Loads entries into a dictionary of free electron states. The dictionary
	created by this function is always available at
	:data:`rmg.chem.electronStates`.
	"""
	electronStates = {}
	electronStates['0'] = ElectronState('0', 0, [1])
	electronStates['1'] = ElectronState('1', 1, [2])
	electronStates['2'] = ElectronState('2', 2, [1,3])
	electronStates['2S'] = ElectronState('2S', 2, [1])
	electronStates['2T'] = ElectronState('2T', 2, [3])
	electronStates['3'] = ElectronState('3', 3, [2,4])
	electronStates['4'] = ElectronState('4', 4, [1,3,5])
	return electronStates

# The dictionary of electron states, accessible by label.
electronStates = loadElectronStates()

################################################################################

cdef class BondType:
	"""
	Represent a type of chemical bond. Each bond type has a unique string
	`label`; a unique string `name`; a numeric bond `order`; an integral
	`piElectrons`, the number of pi electrons; and a string `location` with
	bond geometry information (i.e. 'cis' or 'trans').

	This class is specifically for properties that all bonds of the same type
	share. Ideally there is only one instance of this class for each bond type.
	"""

	cdef public str label
	cdef public str name
	cdef public float order
	cdef public int piElectrons
	cdef public str location

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

	def __reduce__(self):
		"""
		Used for pickling.
		"""
		return (BondType, (self.label, self.name, self.order, self.piElectrons, self.location))

	cpdef bint equivalent(BondType self, BondType other):
		"""
		Returns :data:`True` if two bond types are equivalent or :data:`False`
		otherwise. The method is strictly a comparison of the bond type via the
		`label` attribute. Unlike with atom types, no fuzzy matching of bond
		types is currently allowed (i.e. 'D' is different from 'Dcis' and
		'Dtrans').
		"""
		return self.label == other.label

	def __repr__(self):
		"""x.__repr__() <==> repr(x)"""
		return +"%s.BondType('%s','%s',%g,%g,location='%s')"%(self.__module__,self.label,self.name,self.order,self.piElectrons,self.location)

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

cdef class Atom:
	"""
	Represent an atom in a chemical species or functional group. The `atomType`
	and `electronState` attributes contain lists of the allowed atom types and
	electron states, respectively, for this atom. The `charge` attribute stores
	the resulting formal charge. The `label` attribute can be used to tag
	individual atoms, e.g. center atoms or reactive sites in functional groups.
	"""

## moved to chem.pyxd
#	cdef public list _atomType
#	cdef public list _electronState
#	cdef public int charge
#	cdef public str label

	def __init__(self, atomType='R', electronState='0', charge=0, label=''):
		"""
		Initialize an atom object.
		"""
		self.atomType = atomType
		self.electronState = electronState
		self.charge = charge
		self.label = label
		
		# for Extended Connectivity; as introduced by Morgan (1965)
		# http://dx.doi.org/10.1021/c160017a018
		self.connectivity_value_1 = 0
		self.connectivity_value_2 = 0
		self.connectivity_value_3 = 0

	def __reduce__(self):
		"""
		Used for pickling.
		"""
		return (Atom, (self.atomType, self.electronState, self.charge, self.label))

	#def __repr__(self):
	#	"""x.__repr__() <==> repr(x)"""
	#	return "Atom(atomType=%s,electronState=%s,charge=%s,label='%s')" % (self.atomType, self.electronState, self.charge, self.label)

	def getAtomType(self):
		"""
		Get the list of allowed atom types. If the list is of length 1, the
		lone item in the list is returned instead.
		"""
		if len(self._atomType) == 1: return self._atomType[0]
		else: return self._atomType

	def setAtomType(self, atomType):
		"""
		Set the atom type that this atom represents. The `atomType`
		parameter is any of:
		
		* A string containing the label of a single atom type
		
		* An :class:`AtomType` object representing the atom type
		
		* A list containing one or more of each of the above

		In all cases, the data will be stored internally as a list of
		:class:`AtomType` objects.
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
				raise Exception('Invalid atom type "' + str(atomType) + '".')
		
	atomType = property(getAtomType, setAtomType)

	def getElectronState(self):
		"""
		Get the list of allowed electronic states. If the list is of length 1,
		the lone item in the list is returned instead.
		"""
		if len(self._electronState) == 1: return self._electronState[0]
		else: return self._electronState

	def setElectronState(self, electronState):
		"""
		Set the electron state that this atom represents. The `electronState`
		parameter is any of:

		* A string containing the label of a single electron state

		* An :class:`ElectronState` object representing the electron state

		* A list containing one or more of each of the above

		In all cases, the data will be stored internally as a list of
		:class:`ElectronState` objects.
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

	cpdef bint equivalent(Atom self, Atom other):
		"""
		Return :data:`True` if `self` and `other` are equivalent atoms, i.e.
		they have the same element and electronic state.
		"""

		cdef bint atomTypesMatch = False
		cdef bint electronStatesMatch = False
		cdef AtomType atomType1, atomType2
		cdef ElectronState elecState1, elecState2

		for atomType1 in self._atomType:
			for atomType2 in other._atomType:
				if atomType1.equivalent(atomType2): atomTypesMatch = True

		for elecState1 in self._electronState:
			for elecState2 in other._electronState:
				if elecState1.equivalent(elecState2): electronStatesMatch = True

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

	def getElement(self):
		"""
		Return the element that this atom represents. If there are multiple
		atom types, then this method returns :data:`None`.
		"""
		if len(self._atomType) == 0 or len(self._atomType) > 1: return None
		return self.atomType.element

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
		return not self.isHydrogen()

	def hasFreeElectron(self):
		"""
		Return :data:`True` if the atom has one or more unpaired electrons and
		:data:`False` otherwise.
		"""
		return self.electronState.order > 0

	def getFreeElectronCount(self):
		"""
		Return the number of unpaired electrons.
		"""
		return self.electronState.order

	def canIncreaseFreeElectron(self):
		"""
		Return :data:`True` if the number of unpaired electrons on this atom can
		be increased by one. This method restricts the radical order to three or
		fewer.
		"""
		return self.getFreeElectronCount() < 4

	def canDecreaseFreeElectron(self):
		"""
		Return :data:`True` if the number of unpaired electrons on this atom can
		be decreased by one. This method requires that there be at least one
		free electron on this atom.
		"""
		return self.getFreeElectronCount() > 0

	def increaseFreeElectron(self):
		"""
		Increase the number of unpaired electrons on this atom by one.
		"""
		if self.electronState.label == '0': self.electronState = electronStates['1']
		elif self.electronState.label == '1': self.electronState = electronStates['2']
		elif self.electronState.label == '2': self.electronState = electronStates['3']
		elif self.electronState.label == '2S': self.electronState = electronStates['3']
		elif self.electronState.label == '2T': self.electronState = electronStates['3']
		elif self.electronState.label == '3': self.electronState = electronStates['4']
		else:
			raise Exception('Cannot increase the radical number of this atom.')

	def decreaseFreeElectron(self):
		"""
		Decrease the number of unpaird electrons on this atom by one.
		"""
		if self.electronState.label == '1': self.electronState = electronStates['0']
		elif self.electronState.label == '2': self.electronState = electronStates['1']
		elif self.electronState.label == '2S': self.electronState = electronStates['1']
		elif self.electronState.label == '2T': self.electronState = electronStates['1']
		elif self.electronState.label == '3': self.electronState = electronStates['2']
		elif self.electronState.label == '4': self.electronState = electronStates['3']
		else:
			raise Exception('Cannot decrease the radical number of this atom.')

################################################################################

cdef class Bond:
	"""
	Represent a bond in a chemical species. Each bond has a list `atoms` of
	length two containing the two atoms in the bond and a `bondType` object,
	stored internally as a :class:`BondType` object.
	"""
# moved to pxd file:
#	cdef public list atoms
#	cdef public list _bondType
#	cdef public list atoms
#	cdef public list _bondType

	def __init__(self, atoms, bondType='S'):
		self.bondType = bondType
		self.atoms = atoms

	def __repr__(self):
		"""x.__repr__() <==> repr(x)"""
		#return self.__module__+".Bond(%s,bondType='%s')"%(repr(self.atoms),self.bondType.label)
		return "Bond(%s,bondType='%s')"%(repr(self.atoms),self.bondType.label)

	def __reduce__(self):
		"""
		Used for pickling.
		"""
		return (Bond, (self.atoms, self.bondType))

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

	cpdef bint equivalent(Bond self, Bond other):
		"""
		Return :data:`True` if `self` and `other` are equivalent bonds, i.e.
		they have the same bond type.
		"""
		cdef BondType bondType1, bondType2

		for bondType1 in self._bondType:
			for bondType2 in other._bondType:
				if bondType1.equivalent(bondType2): return True
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
		return self.bondType.order == 1

	def isDouble(self):
		"""
		Return :data:`True` if the bond represents a double bond and
		:data:`False` otherwise.
		"""
		return self.bondType.order == 2

	def isTriple(self):
		"""
		Return :data:`True` if the bond represents a triple bond and
		:data:`False` otherwise.
		"""
		return self.bondType.order == 3

	def isBenzene(self):
		"""
		Return :data:`True` if the bond represents a benzene (aromatic) bond and
		:data:`False` otherwise.
		"""
		return self.bondType.order == 1.5

	def canIncreaseOrder(self):
		"""
		Return :data:`True` if the bond order can be increased by one.
		"""
		return self.isSingle() or self.isDouble()

	def canDecreaseOrder(self):
		"""
		Return :data:`True` if the bond order can be decreased by one without
		breaking.
		"""
		return self.isDouble() or self.isTriple()

	def increaseOrder(self):
		"""
		Increase the bond order by one.
		"""
		if self.isSingle():		self.bondType = bondTypes['D']
		elif self.isDouble():	self.bondType = bondTypes['T']
		else:
			logging.exception('Cannot increase the bond order of this bond.')

	def decreaseOrder(self):
		"""
		Decrease the bond order by one.
		"""
		if self.isDouble():		self.bondType = bondTypes['S']
		elif self.isTriple():	self.bondType = bondTypes['D']
		else:
			logging.exception('Cannot decrease the bond order of this bond.')


################################################################################
