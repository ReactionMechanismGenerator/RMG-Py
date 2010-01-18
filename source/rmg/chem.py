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
Contains classes describing simple chemical entities:

* :class:`Element` - A chemical element

* :class:`AtomType` - A type of atom, combining element information with local bonding information

* :class:`ElectronState` - A free electron state for an atom, describing the number of free electrons and the spin multiplicity

* :class:`BondType` - A type of bond, reflecting the bond order

* :class:`Atom` - A chemical atom, combining one (or more) atom types with one (or more) electron states

* :class:`Bond` - A chemical bond, combining one (or more) bond types

This module can be compiled using Cython to a shared library, which provides a
significant speed boost to running in pure Python mode.

"""

import cython

################################################################################

class Element:
	"""
	A single chemical element. The attributes are:

	=========  =================================================================
	Attribute  Description
	=========  =================================================================
	`number`   The atomic number of the element
	`name`     The name of the element
	`symbol`   The symbol used for the element
	`mass`     The mass of the element in kg/mol
	`valence`  A list of the allowed numbers of bonds to this element
	=========  =================================================================
	
	This class is specifically for properties that all atoms of the same element
	share. Ideally there is only one instance of this class for each element.
	"""

	def __init__(self, number, name, symbol, mass, valence):
		"""
		Initialize a chemical element.
		"""
		self.number = number
		self.name = name
		self.symbol = intern(symbol)
		self.mass = mass
		self.valence = valence

	def __repr__(self):
		"""
		Return a representation that can be used to reconstruct the object.
		"""
		return "Element(%s, '%s', '%s', %s, %s)" % (self.number, self.name, self.symbol, self.mass, self.valence)
	
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
cython.declare(AtomType_R=str)
AtomType_R = intern('R')
cython.declare(AtomType_RnotH=str)
AtomType_RnotH = intern('R!H')
cython.declare(AtomType_H=str)
AtomType_H = intern('H')

class AtomType:
	"""
	A type of atom which combines the element information with information
	about the local bonding structure. The attributes are:

	=============== ============================================================
	Attribute       Description
	=============== ============================================================
	`label`         A unique string identifier for the atom type
	`element`       The :class:`Element` object this atom type represents, if
	                any
	`description`   A one-line description of the atom type
	`doubleBonds`   The number of double bonds to this atom type
	`tripleBonds`   The number of triple bonds to this atom type
	`benzeneBonds`  The number of benzene bonds to this atom type
	`formBond`      The atom type(s) resulting from an adjacent bond formation
	`breakBond`     The atom type(s) resulting from an adjacent bond breaking
	`incrementBond` The atom type(s) resulting from an adjacent bond order
	                increment
	`decrementBond` The atom type(s) resulting from an adjacent bond order
	                decrement
	=============== ============================================================

	This class is specifically for properties that all atoms of the same
	type share. Ideally there is only one instance of this class for each
	atom type.
	"""

	def __init__(self, label='', element=None, doubleBonds=None, tripleBonds=None, benzeneBonds=None, description=''):
		"""
		Initialize an atom type.
		"""
		self.label = intern(label)
		self.element = element
		self.doubleBonds = doubleBonds
		self.tripleBonds = tripleBonds
		self.benzeneBonds = benzeneBonds
		self.description = description
		self.formBond = None
		self.breakBond = None
		self.incrementBond = None
		self.decrementBond = None

	def __repr__(self):
		"""
		Return a representation that can be used to reconstruct the object.
		"""
		return "AtomType('%s', %s, %s, %s, %s, '%s')" % (self.label, self.element, self.doubleBonds, self.tripleBonds, self.benzeneBonds, self.description)
	
	def __reduce__(self):
		"""
		Used for pickling.
		"""
		return (AtomType, (self.label, self.element, self.doubleBonds, self.tripleBonds, self.benzeneBonds, self.description))

	def setActions(self, formBond=None, breakBond=None, incrementBond=None, decrementBond=None):
		"""
		Set the atom types that result from the application of the
		relevant chemical actions. For atom types, these actions are:

		* FORM_BOND - form a single bond between this atom and another atom

		* BREAK_BOND - break a single bond between this atom and another atom

		* INCREMENT_BOND - increment the order of a bond involving this atom by one

		* DECREMENT_BOND - decrement the order of a bond involving this atom by one

		If an action cannot be applied to an atom type, it should be set
		to :data:`None`. If multiple atom types can result from an action, a
		list should be used; RMG will eliminate atom types based on the
		resulting structure, if possible.
		"""
		self.formBond = formBond
		self.breakBond = breakBond
		self.incrementBond = incrementBond
		self.decrementBond = decrementBond

	def equivalent(self, other):
		"""
		Returns :data:`True` if two atom types are equivalent or :data:`False`
		otherwise. Respects wildcards, e.g. returns True for {R!H}=={C}
		"""
		
		# If either is a generic atom type, then always return True
		if self.label is AtomType_R or other.label is AtomType_R: 
			return True
		# If either is a generic non-hydrogen atom type, then return
		# True if any atom type in the remaining one is non-hydrogen
		elif self.label is AtomType_RnotH:
			if not other.label is AtomType_H:
				#logging.debug('I think %s == %s'%(self.label, other.label))
				return True
			else: 
				return False
		elif other.label is AtomType_RnotH:
			if not self.label is AtomType_H:
				#logging.debug('I think %s == %s'%(self.label, other.label))
				return True
			else: 
				return False
		# If either represents an element without surrounding bond info,
		# match remaining to any with the same element
		elif self.label is self.element.symbol is other.element.symbol:
			return True
		elif other.label is other.element.symbol is self.element.symbol:
			return True
		# Otherwise labels must match exactly
		elif self.label is other.label:
			return True
		return False

################################################################################

def loadAtomTypes():
	"""
	Loads entries into a dictionary of atom types.
	"""

	# Functional group atom types
	atomTypes = {}

	atomTypes['H']		= AtomType('H', 	elements['H'], 	description='hydrogen')
	atomTypes['C']		= AtomType('C', 	elements['C'], 	description='carbon')
	atomTypes['N']		= AtomType('N', 	elements['N'], 	description='nitrogen')
	atomTypes['O']		= AtomType('O', 	elements['O'], 	description='oxygen')
	atomTypes['F']		= AtomType('F', 	elements['F'], 	description='fluorine')
	atomTypes['Ne'] 	= AtomType('Ne', 	elements['Ne'], description='neon')
	atomTypes['Si'] 	= AtomType('Si', 	elements['Si'], description='silicon')
	atomTypes['P']		= AtomType('P', 	elements['P'],	description='phosphorus')
	atomTypes['S']		= AtomType('S', 	elements['S'],	description='sulfur')
	atomTypes['Cl'] 	= AtomType('Cl', 	elements['Cl'], description='chlorine')
	atomTypes['Ar'] 	= AtomType('Ar', 	elements['Ar'], description='argon')
	atomTypes['Br'] 	= AtomType('Br', 	elements['Br'], description='bromine')
	atomTypes['I']		= AtomType('I', 	elements['I'],	description='iodine')

	atomTypes['R']		= AtomType('R', 	None, 			description='generic functional group')
	atomTypes['R!H'] 	= AtomType('R!H',	None, 			description='generic non-hydrogen functional group')

	atomTypes['Cs'] 	= AtomType('Cs', 	elements['C'], 	doubleBonds=0, tripleBonds=0, benzeneBonds=0, description='carbon with four single bonds')
	atomTypes['Cd'] 	= AtomType('Cd', 	elements['C'], 	doubleBonds=1, tripleBonds=0, benzeneBonds=0, description='carbon with one double bond and two single bonds')
	atomTypes['Cdd'] 	= AtomType('Cdd', 	elements['C'], 	doubleBonds=2, tripleBonds=0, benzeneBonds=0, description='carbon with two double bonds')
	atomTypes['Ct'] 	= AtomType('Ct', 	elements['C'], 	doubleBonds=0, tripleBonds=1, benzeneBonds=0, description='carbon with one triple bond and one single bond')
	atomTypes['Cb'] 	= AtomType('Cb', 	elements['C'], 	doubleBonds=0, tripleBonds=0, benzeneBonds=2, description='carbon belonging to a benzene ring')
	atomTypes['Cbf'] 	= AtomType('Cbf', 	elements['C'], 	doubleBonds=0, tripleBonds=0, benzeneBonds=3, description='carbon belonging to a fused benzene ring')

	atomTypes['Os'] 	= AtomType('Os', 	elements['O'], 	doubleBonds=0, tripleBonds=0, benzeneBonds=0, description='oxygen with two single bonds')
	atomTypes['Od'] 	= AtomType('Od', 	elements['O'], 	doubleBonds=1, tripleBonds=0, benzeneBonds=0, description='oxygen with one double bond')

	atomTypes['Sis'] 	= AtomType('Sis', 	elements['Si'], doubleBonds=0, tripleBonds=0, benzeneBonds=0, description='silicon with four single bonds')
	atomTypes['Sid'] 	= AtomType('Sid', 	elements['Si'], doubleBonds=1, tripleBonds=0, benzeneBonds=0, description='silicon with one double bond and two single bonds')
	atomTypes['Sidd'] 	= AtomType('Sidd', 	elements['Si'], doubleBonds=2, tripleBonds=0, benzeneBonds=0, description='silicon with two double bonds')
	atomTypes['Sit'] 	= AtomType('Sit', 	elements['Si'], doubleBonds=0, tripleBonds=1, benzeneBonds=0, description='silicon with one triple bond and one single bond')
	atomTypes['Sib'] 	= AtomType('Sib', 	elements['Si'], doubleBonds=0, tripleBonds=0, benzeneBonds=2, description='silicon belonging to a benzene ring')
	atomTypes['Sibf'] 	= AtomType('Sibf', 	elements['Si'], doubleBonds=0, tripleBonds=0, benzeneBonds=3, description='silicon belonging to a fused benzene ring')

	atomTypes['H'].setActions(formBond=atomTypes['H'], breakBond=atomTypes['H'], incrementBond=None, decrementBond=None)
	atomTypes['C'].setActions(formBond=atomTypes['C'], breakBond=atomTypes['C'], incrementBond=atomTypes['C'], decrementBond=atomTypes['C'])
	atomTypes['O'].setActions(formBond=atomTypes['O'], breakBond=atomTypes['O'], incrementBond=atomTypes['O'], decrementBond=atomTypes['O'])

	atomTypes['R'].setActions(formBond=atomTypes['R'], breakBond=atomTypes['R'], incrementBond=atomTypes['R'], decrementBond=atomTypes['R'])
	atomTypes['R!H'].setActions(formBond=atomTypes['R!H'], breakBond=atomTypes['R!H'], incrementBond=atomTypes['R!H'], decrementBond=atomTypes['R!H'])

	atomTypes['Cs'].setActions(formBond=atomTypes['Cs'], breakBond=atomTypes['Cs'], incrementBond=atomTypes['Cd'], decrementBond=None)
	atomTypes['Cd'].setActions(formBond=atomTypes['Cd'], breakBond=atomTypes['Cd'], incrementBond=[atomTypes['Cdd'],atomTypes['Ct']], decrementBond=atomTypes['Cs'])
	atomTypes['Cdd'].setActions(formBond=atomTypes['Cdd'], breakBond=atomTypes['Cdd'], incrementBond=None, decrementBond=atomTypes['Cd'])
	atomTypes['Ct'].setActions(formBond=atomTypes['Ct'], breakBond=atomTypes['Ct'], incrementBond=None, decrementBond=atomTypes['Cd'])
	atomTypes['Cb'].setActions(formBond=atomTypes['Cb'], breakBond=atomTypes['Cb'], incrementBond=None, decrementBond=None)
	atomTypes['Cbf'].setActions(formBond=atomTypes['Cbf'], breakBond=atomTypes['Cbf'], incrementBond=None, decrementBond=None)

	atomTypes['Os'].setActions(formBond=atomTypes['Os'], breakBond=atomTypes['Os'], incrementBond=atomTypes['Od'], decrementBond=None)
	atomTypes['Od'].setActions(formBond=atomTypes['Od'], breakBond=atomTypes['Od'], incrementBond=None, decrementBond=atomTypes['Os'])

	return atomTypes

# The dictionary of elements, accessible by atomic number, symbol, or name.
atomTypes = loadAtomTypes()

################################################################################

class ElectronState:
	"""
	A single free electron state for an atom. The attributes are:

	===========  ===============================================================
	Attribute    Description
	===========  ===============================================================
	`label`      A unique string identifier for the electron state
	`order`      The number of free electrons
	`spin`       A list of the allowed spin states for polyradicals (singlet =
	             1, doublet = 2, etc.)
	`increment`  The electron state that results from incrementing the radical
	             count, if any
	`decrement`  The electron state that results from decrementing the radical
	             count, if any
	===========  ===============================================================

	This class is specifically for properties that all free electron states
	share. Ideally there is only one instance of this class for each free
	electron state.
	"""

	def __init__(self, label='', order=0, spin=None, increment=None, decrement=None):
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
		self.increment = increment
		self.decrement = decrement

	def __repr__(self):
		"""
		Return a representation that can be used to reconstruct the object.
		"""
		return "ElectronState('%s', %s, %s)" % (self.label, self.order, self.spin)
	
	def __reduce__(self):
		"""
		Used for pickling.
		"""
		return (ElectronState, (self.label, self.order, self.spin, None, None))

	def setActions(self, increment=None, decrement=None):
		"""
		Set the electron states that result from the application of the
		relevant chemical actions. For electron states, these actions are:

		* INCREMENT_RADICAL - increment the free electron count by one

		* DECREMENT_RADICAL - decrement the free electron count by one

		If an action cannot be applied to an electron state, it should be set
		to :data:`None`.
		"""
		self.increment = increment
		self.decrement = decrement

	def equivalent(self, other):
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

	# Set increment and decrement attributes
	electronStates['0'].setActions(increment=electronStates['1'])
	electronStates['1'].setActions(increment=electronStates['2'], decrement=electronStates['0'])
	electronStates['2'].setActions(increment=electronStates['3'], decrement=electronStates['1'])
	electronStates['2S'].setActions(increment=electronStates['3'], decrement=electronStates['1'])
	electronStates['2T'].setActions(increment=electronStates['3'], decrement=electronStates['1'])
	electronStates['3'].setActions(increment=electronStates['4'], decrement=electronStates['2'])
	electronStates['4'].setActions(decrement=electronStates['3'])

	return electronStates

# The dictionary of electron states, accessible by label.
electronStates = loadElectronStates()

################################################################################

class BondType:
	"""
	A type of chemical bond. The attributes are:

	=============  =============================================================
	Attribute      Description
	=============  =============================================================
	`label`        A unique string identifier for the bond type
	`name`         A longer descriptor of the bond type
	`order`        The bond order
	`piElectrons`  The number of pi electrons that contribute to the bond
	`location`     Geometric description of bond if needed (i.e. 'cis' or 'trans')
	=============  =============================================================

	This class is specifically for properties that all bonds of the same type
	share. Ideally there is only one instance of this class for each bond type.
	"""

	def __init__(self, label='', name='', order=0, piElectrons=0, location=''):
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
		"""
		Return a representation that can be used to reconstruct the object.
		"""
		return "BondType('%s', '%s', %g, %g, location='%s')" % (self.label, self.name, self.order, self.piElectrons, self.location)

	def __reduce__(self):
		"""
		Used for pickling.
		"""
		return (BondType, (self.label, self.name, self.order, self.piElectrons, self.location))

	def equivalent(self, other):
		"""
		Returns :data:`True` if two bond types are equivalent or :data:`False`
		otherwise. The method is strictly a comparison of the bond type via the
		`label` attribute. Unlike with atom types, no fuzzy matching of bond
		types is currently allowed (i.e. 'D' is different from 'Dcis' and
		'Dtrans').
		"""
		return self.label == other.label

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

class InvalidChemicalActionException(Exception):

	"""
	An exception used attempting to carry out a chemical action (bond formation
	or breaking, increasing or decreasing free radical count on an atom, etc.)
	that is not permitted.
	"""

	def __init__(self, msg):
		self.msg = msg

	def __str__(self):
		return self.msg

################################################################################

class Atom(object):
	"""
	Represent an atom in a chemical species or functional group. The attributes
	are:

	===============  ===========================================================
	Attribute        Description
	===============  ===========================================================
	`atomType`       The atom's type or a list of the allowed types
	`electronState`  The atom's free electron state or a list of the allowed
	                 states
	`charge`         The formal charge of the atom (always zero)
	`label`          A string used to tag individual atoms, e.g. center atoms or
	                 reactive sites
	===============  ===========================================================

	"""

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
		self.connectivity1 = -1
		self.connectivity2 = -1
		self.connectivity3 = -1
		self.sorting_label = -1
	
	def __repr__(self):
		"""
		Return a representation that can be used to reconstruct the object.
		"""
		return "Atom(atomType=%s, electronState=%s, charge=%s, label='%s')" % (self.atomType, self.electronState, self.charge, self.label)

	def __reduce__(self):
		"""
		Used for pickling.
		"""
		return (Atom, ([a.label for a in self._atomType], [e.label for e in self._electronState], self.charge, self.label))

	def equals(self, other):
		"""
		Equality comparison.
		"""
		return (self.atomType.label == other.atomType.label and
			self.electronState.label == other.electronState.label and
			self.charge == other.charge and
			self.label == other.label)

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

	def equivalent(self, other):
		"""
		Return :data:`True` if `self` and `other` are equivalent atoms, i.e.
		they have the same atom type and electronic state. If one or both of
		the atoms being tested for equivalence have multiple atom types,
		:data:`True` is returned if *any* atom type in one atom matches *any*
		atom type in the other. This is also done for electron states. Both an
		atom types match and an electron state match must be found in order to
		return :data:`True`.
		"""

		atomTypesMatch = cython.declare(cython.bint, False)
		electronStatesMatch = cython.declare(cython.bint, False)
		
		atomType1 = cython.declare(AtomType)
		atomType2 = cython.declare(AtomType)
		elecState1 = cython.declare(ElectronState)
		elecState2 = cython.declare(ElectronState)
		
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
		otherwise. An atom is a center atom if it has any form of label.
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
		`symbol` and :data:`False` otherwise. Returns False if atom is a wildcard
		(e.g. in a functional group)
		"""
		if len(self._atomType) == 0 or len(self._atomType) > 1: return False
		elif self.atomType.element is None: return False
		else: return self.atomType.element is elements[symbol]

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
		for electronState in self._electronState:
			if electronState.increment is None: return False
		return True
	
	def canDecreaseFreeElectron(self):
		"""
		Return :data:`True` if the number of unpaired electrons on this atom can
		be decreased by one. This method requires that there be at least one
		free electron on this atom.
		"""
		for electronState in self._electronState:
			if electronState.decrement is None: return False
		return True
	
	def increaseFreeElectron(self):
		"""
		Increase the number of unpaired electrons on this atom by one.
		"""
		for i in range(len(self._electronState)):
			if self._electronState[i].increment is None:
				raise InvalidChemicalActionException('Cannot increase the radical number of this atom.')
			else:
				self._electronState[i] = self._electronState[i].increment

	def decreaseFreeElectron(self):
		"""
		Decrease the number of unpaired electrons on this atom by one.
		"""
		for i in range(len(self._electronState)):
			if self._electronState[i].decrement is None:
				raise InvalidChemicalActionException('Cannot decrease the radical number of this atom.')
			else:
				self._electronState[i] = self._electronState[i].decrement

	def __applyAction(self, action, bonds, unique):

		newAtomTypes = []

		# Collect possible new atom types
		for atomType in self._atomType:
			if action == 'FORM_BOND':			newAtomType = atomType.formBond
			elif action == 'BREAK_BOND':		newAtomType = atomType.breakBond
			elif action == 'INCREMENT_BOND':	newAtomType = atomType.incrementBond
			elif action == 'DECREMENT_BOND':	newAtomType = atomType.decrementBond
			if newAtomType is None:
				pass
			elif isinstance(newAtomType, list):
				newAtomTypes.extend(newAtomType)
			else:
				newAtomTypes.append(newAtomType)

		# Eliminate duplicates
		newAtomTypes = list(set(newAtomTypes))

		# Count numbers of each higher-order bond type
		wildcardBond = False
		doubleBonds = 0; tripleBonds = 0; benzeneBonds = 0
		for atom2, bond12 in bonds.iteritems():
			if isinstance(bond12.bondType, list): wildcardBond = True
			elif bond12.isDouble(): doubleBonds += 1
			elif bond12.isTriple(): tripleBonds += 1
			elif bond12.isBenzene(): benzeneBonds += 1

		# Eliminate impossible atom types
		if not wildcardBond:
			atomTypesToRemove = []
			for atomType in newAtomTypes:
				if (atomType.doubleBonds != doubleBonds and atomType.doubleBonds is not None) or \
					(atomType.tripleBonds != tripleBonds and atomType.tripleBonds is not None) or \
					(atomType.benzeneBonds != 0 and benzeneBonds == 0 and atomType.benzeneBonds is not None) or \
					(atomType.benzeneBonds == 0 and benzeneBonds != 0 and atomType.benzeneBonds is not None):
					atomTypesToRemove.append(atomType)
			for atomType in atomTypesToRemove:
				newAtomTypes.remove(atomType)
			
		# Raise exception if we don't have at least one atom type remaining
		if len(newAtomTypes) == 0:
			raise InvalidChemicalActionException('No atom types remaining after %s action with original atom type(s) %s.' % (action, self._atomType))
		# Raise exception if we have more than one atom type remaining and the
		# unique flag is set
		elif len(newAtomTypes) > 1 and unique:
			raise InvalidChemicalActionException('Unable to determine a unique atom type after %s action with original atom type(s) %s.' % (action, self._atomType))

		# Set atom type(s) to atom
		self._atomType = newAtomTypes


	def formBond(self, bonds, unique):
		"""
		Modify the atom type(s) of this atom in response to a FORM_BOND 
		chemical action.
		"""
		self.__applyAction('FORM_BOND', bonds, unique)

	def breakBond(self, bonds, unique):
		"""
		Modify the atom type(s) of this atom in response to a BREAK_BOND
		chemical action.
		"""
		self.__applyAction('BREAK_BOND', bonds, unique)

	def incrementBond(self, bonds, unique):
		"""
		Modify the atom type(s) of this atom in response to a INCREMENT_BOND
		chemical action.
		"""
		self.__applyAction('INCREMENT_BOND', bonds, unique)

	def decrementBond(self, bonds, unique):
		"""
		Modify the atom type(s) of this atom in response to a DECREMENT_BOND
		chemical action.
		"""
		self.__applyAction('DECREMENT_BOND', bonds, unique)


################################################################################

class Bond(object):
	"""
	A chemical bond between atoms. The attributes are:

	===========  ===========================================================
	Attribute    Description
	===========  ===========================================================
	`atoms`      A list of the two atoms linked by the bond
	`bondType`   The bond's type or a list of the allowed types
	===========  ===========================================================

	"""

	def __init__(self, atoms=[None, None], bondType='S'):
		self.bondType = bondType
		self.atoms = atoms

	def __repr__(self):
		"""
		Return a representation that can be used to reconstruct the object.
		"""
		return "Bond(%s, bondType='%s')" % (repr(self.atoms), self.bondType.label)

	def __reduce__(self):
		"""
		Used for pickling.
		"""
		return (Bond, (self.atoms, [b.label for b in self._bondType]))

	def equals(self, other):
		"""
		Equality comparison.
		"""
		return (((self.atoms[0].equals(other.atoms[0]) and self.atoms[1].equals(other.atoms[1])) or
			(self.atoms[0].equals(other.atoms[1]) or self.atoms[1].equals(other.atoms[0]))) and
			self.bondType.label == other.bondType.label)

	def getBondType(self):
		"""
		Get the list of allowed bond types. If the list is of length 1, the
		lone item in the list is returned instead.
		"""
		if len(self._bondType) == 1: return self._bondType[0]
		else: return self._bondType

	def setBondType(self, bondType):
		"""
		Set the bond type that this atom represents. The `bondType`
		parameter is any of:

		* A number representing the bond order

		* A string containing the label of a single bond type

		* An :class:`BondType` object representing the bond type

		* A list containing one or more of each of the above

		In all cases, the data will be stored internally as a list of
		:class:`BondType` objects.
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
		they have the same bond type. If one or both of the bonds being tested
		for equivalence have multiple bond types, :data:`True` is returned if
		*any* bond type in one bond matches *any* bond type in the other.
		"""
		bondType1 = cython.declare(BondType)
		bondType2 = cython.declare(BondType)
		
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
		for bondType in self._bondType:
			if bondType.order != 1 and bondType.order != 2: return False
		return True

	def canDecreaseOrder(self):
		"""
		Return :data:`True` if the bond order can be decreased by one without
		breaking.
		"""
		for bondType in self._bondType:
			if bondType.order != 2 and bondType.order != 3: return False
		return True
	
	def increaseOrder(self):
		"""
		Increase the bond order by one.
		"""
		for i in range(len(self._bondType)):
			if self._bondType[i].order == 1:		self._bondType[i] = bondTypes['D']
			elif self._bondType[i].order == 2:		self._bondType[i] = bondTypes['T']
			else:
				raise InvalidChemicalActionException('Cannot increase the bond order of this bond.')

	def decreaseOrder(self):
		"""
		Decrease the bond order by one.
		"""
		for i in range(len(self._bondType)):
			if self._bondType[i].order == 2:		self._bondType[i] = bondTypes['S']
			elif self._bondType[i].order == 3:		self._bondType[i] = bondTypes['D']
			else:
				raise InvalidChemicalActionException('Cannot decrease the bond order of this bond.')


################################################################################
