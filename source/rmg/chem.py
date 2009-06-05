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

import quantities as pq
import pybel
import openbabel

import graph

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
		return self.atomType == other.atomType and self.electronState == other.electronState

	def copy(self):
		"""
		Generate a copy of the current atom.
		"""
		return Atom(self.element, self.electronState, self.charge, self.label)

	def isCenter(self):
		"""
		Return :data:`True` if the atom is a center atom and :data:`False`
		otherwise.
		"""
		return len(self.label) > 0

	def isHydrogen(self):
		"""
		Return :data:`True` if the atom is a hydrogen atom and :data:`False`
		otherwise.
		"""
		return self.atomType.element.symbol == 'H'

	def isNonHydrogen(self):
		"""
		Return :data:`True` if the atom is not a hydrogen atom and :data:`False`
		otherwise.
		"""
		return self.atomType.element.symbol != 'H'

	def increaseRadical(self):
		"""
		Increase the number of radical electrons on this atom by one.
		"""
		if self.electronState.label == '0': self.electronState = electronStates['1']
		elif self.electronState.label == '1': self.electronState = electronStates['2']
		elif self.electronState.label == '2': self.electronState = electronStates['3']
		elif self.electronState.label == '2S': self.electronState = electronStates['3']
		elif self.electronState.label == '2T': self.electronState = electronStates['3']
		else:
			logging.exception('Cannot increase the radical number of this atom.')

	def decreaseRadical(self):
		"""
		Decrease the number of radical electrons on this atom by one.
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
		return self.bondType == other.bondType

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

		atoms = []; bonds = []; atomdict = {}

		lines = adjlist.splitlines()

		for line in lines[1:]:

			data = line.split()

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

	def radicalCount(self):
		"""
		Get the number of radicals in the structure.
		"""
		radical = 0
		for atom in self.atoms():
			radical += atom.electronState.order
		return radical

	def generateResonanceIsomers(self):
		"""
		Generate a list of all of the resonance isomers of this structure.
		"""

		isomers = [self]

		# Radicals
		if self.radicalCount() > 0:
			# Iterate over resonance isomers
			index = 0
			while index < len(isomers):
				isomer = isomers[index]
				# Iterate over radicals in structure; use indices because the
				# pointer to isomer (but not its contents) may be changing
				for i in range(0, len(isomer.atoms())):
					atom = isomer.atoms()[i]
					paths = isomer.findAllDelocalizationPaths(atom)
					for path in paths:

						# Make a copy of isomer
						oldIsomer = isomer.copy()
						isomers[index] = oldIsomer
						newIsomer = isomer
						isomer = oldIsomer

						# Adjust to (potentially) new resonance isomer
						atom1, atom2, atom3, bond12, bond23 = path
						atom1.decreaseRadical()
						atom3.increaseRadical()
						bond12.increaseOrder()
						bond23.decreaseOrder()

						# Append to isomer list if unique
						found = False
						for isom in isomers:
							if isom.isIsomorphic(newIsomer): found = True
						if not found:
							isomers.append(newIsomer)

				# Move to next resonance isomer
				index += 1

		return isomers

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

	def getThermoData(self):
		"""
		Calculate the thermodynamic parameters for this structure using the
		thermo database to look up group additivity values for each heavy
		(i.e. non-hydrogen) atom.
		"""

		# Initialize thermo data
		self.thermo = ThermoGAData()

		# Iterate over heavy (non-hydrogen) atoms
		for atom in self.atoms():
			if atom.isNonHydrogen():
				thermoData = Structure.thermoDatabase.getThermoData(self, atom)
				if thermoData is not None:
					self.thermo += thermoData

################################################################################

class Species:
	"""
	Represent a chemical species (including all of its resonance forms). Each
	species has a unique integer `id` assigned automatically by RMG and a
	not-necessarily unique string `label`. The *structure* variable contains a
	list of :class:`Structure` objects representing each resonance form. The
	`reactive` flag is :data:`True` if the species can react and :data:`False`
	if it is inert.
	"""

	# A static counter for the number of species created since the RMG job began.
	numSpecies = 0

	def __init__(self, label='', structure=None, reactive=True):
		"""
		Initialize a Species object.
		"""
		Species.numSpecies += 1
		self.id = Species.numSpecies
		self.label = label
		self.structure = [structure]
		self.reactive = reactive

		self.thermo = None
		self.lennardJones = None
		self.spectralData = None

		if structure is not None:
			self.__generateData()

	def __generateData(self):
		"""
		Generate supplemental parameters and information about the species:
		resonance isomers, thermodynamic data, etc.
		"""
		self.structure = self.structure[0].generateResonanceIsomers()
		self.getThermoData()

	def getThermoData(self):
		"""
		Generate thermodynamic data for the species by use of the thermo
		database.
		"""
		for structure in self.structure:
			structure.getThermoData()

	def __str__(self):
		"""
		Return a string representation of the species, in the form 'label(id)'.
		"""
		return self.label + '(' + str(self.id) + ')'

################################################################################

class ThermoGAData:
	"""
	A set of thermodynamic parameters as determined from Benson's group
	additivity data. The attributes are:

	- `H298` = the standard enthalpy of formation at 298 K in J/mol

	- `S298` = the standard entropy of formation at 298 K in J/mol*K

	- `Cp300` = the standard heat capacity at 300 K in J/mol*K

	- `Cp400` = the standard heat capacity at 400 K in J/mol*K

	- `Cp500` = the standard heat capacity at 500 K in J/mol*K

	- `Cp600` = the standard heat capacity at 600 K in J/mol*K

	- `Cp800` = the standard heat capacity at 800 K in J/mol*K

	- `Cp1000` = the standard heat capacity at 1000 K in J/mol*K

	- `Cp1500` = the standard heat capacity at 1500 K in J/mol*K

	- `comment` = a string describing the source of the data
	"""

	CpTlist = pq.Quantity([300.0, 400.0, 500.0, 600.0, 800.0, 1000.0, 1500.0], 'K')

	def __init__(self, H298=0.0, S298=0.0, Cp300=0.0, Cp400=0.0, Cp500=0.0, \
	             Cp600=0.0, Cp800=0.0, Cp1000=0.0, Cp1500=0.0, comment=''):
		"""Initialize a set of group additivity thermodynamic data."""

		self.H298 = H298
		self.S298 = S298
		self.Cp = [Cp300, Cp400, Cp500, Cp600, Cp800, Cp1000, Cp1500]
		self.comment = comment

	def fromDatabase(self, data, comment):
		"""
		Process a list of numbers `data` and associated description `comment`
		generated while reading from a thermodynamic database.
		"""

		if len(data) != 12:
			raise Exception('Invalid list of thermo data; should be a list of numbers of length 12.')

		H298, S298, Cp300, Cp400, Cp500, Cp600, Cp800, Cp1000, Cp1500, \
		  dH, dS, dCp = data
		self.H298 = pq.UncertainQuantity(H298, 'kcal/mol', dH)
		#self.H298.units = 'J/mol'
		self.S298 = pq.UncertainQuantity(S298, 'cal/(mol*K)', dS)
		#self.S298.units = 'J/(mol*K)'
		self.Cp = pq.UncertainQuantity([Cp300, Cp400, Cp500, Cp600, Cp800, Cp1000, Cp1500], 'kcal/(mol*K)', [dCp, dCp, dCp, dCp, dCp, dCp, dCp])
		#self.Cp.units = 'J/(mol*K)'

		#self.Cp[0] = pq.UncertainQuantity(Cp300, 'kcal/(mol*K)', dCp)
		#self.Cp[0].units = 'J/(mol*K)'
		#self.Cp[1] = pq.UncertainQuantity(Cp400, 'kcal/(mol*K)', dCp)
		#self.Cp[1].units = 'J/(mol*K)'
		#self.Cp[2] = pq.UncertainQuantity(Cp500, 'kcal/(mol*K)', dCp)
		#self.Cp[2].units = 'J/(mol*K)'
		#self.Cp[3] = pq.UncertainQuantity(Cp600, 'kcal/(mol*K)', dCp)
		#self.Cp[3].units = 'J/(mol*K)'
		#self.Cp[4] = pq.UncertainQuantity(Cp800, 'kcal/(mol*K)', dCp)
		#self.Cp[4].units = 'J/(mol*K)'
		#self.Cp[5] = pq.UncertainQuantity(Cp1000, 'kcal/(mol*K)', dCp)
		#self.Cp[5].units = 'J/(mol*K)'
		#self.Cp[6] = pq.UncertainQuantity(Cp1500, 'kcal/(mol*K)', dCp)
		#self.Cp[6].units = 'J/(mol*K)'
		self.comment = comment

	def __add__(self, other):
		"""
		Add two sets of thermodynamic data together. All parameters are
		additive.
		"""
		new = ThermoGAData()
		new.H298 = self.H298 + other.H298
		new.S298 = self.S298 + other.S298
		new.Cp = []
		for i in range(len(self.Cp)):
			new.Cp.append(self.Cp + other.Cp)
		new.comment = self.comment + '; ' + other.comment
		return new

	def heatCapacity(self, T):
		"""
		Return the heat capacity at the specified temperature `T`. This is done
		via linear interpolation between the provided values.
		"""
		T.units = 'K'; T = float(T)
		if T < 300.0:
			raise TemperatureOutOfRangeException('No thermodynamic data available for T < 300 K.')
		# Use Cp(1500 K) if T > 1500 K
		elif T > 1500.0: T = 1500.0

		Cpfun = scipy.interpolate.interp1d(ThermoGAData.CpTlist, self.Cp)
		return pq.Quantity(Cpfun(T), 'J/(mol*K)')

	def enthalpy(self, T):
		"""
		Return the enthalpy of formation at the specified temperature `T`.
		"""
		T.units = 'K'; T = float(T)
		if T < 300.0:
			raise TemperatureOutOfRangeException('No thermodynamic data available for T < 300 K.')

	def getCpLinearization(self):
		slope = []; intercept = []
		for i in range(0, len(self.Cp)-1):
			slope.append((self.Cp[i+1] - self.Cp[i]) / (ThermoGAData.CpTlist[i+1] - ThermoGAData.CpTlist[i]))
			print slope, ThermoGAData.CpTlist[i], slope * ThermoGAData.CpTlist[i]
			quit()
			intercept.append(self.Cp[i] - slope[i] * ThermoGAData.CpTlist[i])
		return slope, intercept

	def toXML(self, dom, root):

		thermo = dom.createElement('thermo')
		root.appendChild(thermo)

		self.valueToXML(dom, thermo, 'enthalpyOfFormation', self.H298, '298 K')
		self.valueToXML(dom, thermo, 'entropyOfFormation', self.S298, '298 K')
		self.valueToXML(dom, thermo, 'heatCapacity', self.Cp[0], '300 K')
		self.valueToXML(dom, thermo, 'heatCapacity', self.Cp[1], '400 K')
		self.valueToXML(dom, thermo, 'heatCapacity', self.Cp[2], '500 K')
		self.valueToXML(dom, thermo, 'heatCapacity', self.Cp[3], '600 K')
		self.valueToXML(dom, thermo, 'heatCapacity', self.Cp[4], '800 K')
		self.valueToXML(dom, thermo, 'heatCapacity', self.Cp[5], '1000 K')
		self.valueToXML(dom, thermo, 'heatCapacity', self.Cp[6], '1500 K')

		element = dom.createElement('comment')
		thermo.appendChild(element)
		comment = dom.createTextNode(self.comment)
		element.appendChild(comment)


	def valueToXML(self, dom, root, name, value, temp):
		element = dom.createElement(name)
		root.appendChild(element)

		units = str(value.units).split()[1]

		element.setAttribute('temperature', temp)
		element.setAttribute('units', units)
		element.setAttribute('uncertainty', str(value.uncertainty))

		valueNode = dom.createTextNode(str(float(value)))
		element.appendChild(valueNode)

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
	
