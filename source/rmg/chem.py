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
Contains classes describing chemical items.
"""

import quantities as pq
import pybel
import openbabel

################################################################################

class Element:
    """
	Represent a single chemical element.

    This class is specifically for properties that all atoms of the same
    element share. Ideally there is only one instance of this class for each
    atom variant (isotope, valence, etc.).
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
	return elements
	
elements = loadElements()

################################################################################

class ElectronState:
    """Represent a single free electron state (none, radical, etc.)

    This class is specifically for properties that all free electron states
    share. Ideally there is only one instance of this class for each
    free electron state.
    """

    def __init__(self, label, order, spin=''):
        """
		Initialize a free electron state.

        Parameters:
        label -- A string unique to this free electron state
        order -- The number of free electrons
        spin -- The spin state for polyradicals (singlet = 'S', triplet = 'T')
        """
        self.label = label
        self.order = order
        self.spin = spin

################################################################################

def loadElectronStates():
	"""
	Loads entries into a dictionary of free electron states.
	"""
	electronStates = {}
	electronStates['0'] = ElectronState('0', 0, '')
	electronStates['1'] = ElectronState('1', 1, '')
	electronStates['2S'] = ElectronState('2S', 2, 'S')
	electronStates['2T'] = ElectronState('2T', 2, 'T')
	electronStates['3'] = ElectronState('3', 3, '')
	electronStates['4'] = ElectronState('4', 4, '')
	return electronStates

electronStates = loadElectronStates()

################################################################################
		
class Atom:
	"""
	Represent an atom in a chemical species.
	"""	

	def __init__(self, element=None, electronState=None, charge=0):
		"""
		Initialize an atom object.
		"""
		self.setElement(element)
		self.setElectronState(electronState)
		self.charge = charge

	def setElement(self, element):
		"""
		Set the element that this atom represents.
		"""
		if element.__class__ == str or element.__class__ == unicode:
			element = elements[element]
		if element.__class__ != Element:
			raise Exception('Invalid parameter "element".')
		self.element = element

	def setElectronState(self, electronState):
		"""
		Set the electron state that this atom represents.
		"""
		if electronState.__class__ == str or electronState.__class__ == unicode:
			electronState = electronStates[electronState]
		if electronState.__class__ != ElectronState:
			raise Exception('Invalid parameter "electronState".')
		self.electronState = electronState

################################################################################

class BondType:
    """
	Represent a type of chemical bond.

    This class is specifically for properties that all bonds of the same
    type share. Ideally there is only one instance of this class for each
    bond type.
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

################################################################################

def loadBondTypes():
	"""
	Loads entries into a dictionary of bond types.
	"""
	bondTypes = {}
	bondTypes[1] = bondTypes['S'] = bondTypes['single'] = BondType('S', 'single', 1, 0, '')
	bondTypes[2] = bondTypes['D'] = bondTypes['double'] = BondType('D', 'double', 2, 2, '')
	bondTypes['Dcis'] = BondType('Dcis', 'double_cis', 2, 2, 'cis')
	bondTypes['Dtrans'] = BondType('Dtrans', 'double_trans', 2, 2, 'trans')
	bondTypes[3] = bondTypes['T'] = bondTypes['triple'] = BondType('T', 'triple', 3, 4, '')
	bondTypes[1.5] = bondTypes['B'] = bondTypes['benzene'] = BondType('B', 'benzene', 1.5, 1, '')
	return bondTypes
	
bondTypes = loadBondTypes()

################################################################################

class Bond:
	"""
	Represent a bond in a chemical species.
	"""	

	def __init__(self, atom1, atom2, bondType=''):
		self.setBondType(bondType)
		self.atoms = [atom1, atom2]

	def setBondType(self, bondType):
		"""
		Set the bond type that this bond represents.
		"""
		if bondType.__class__ == str or bondType.__class__ == unicode:
			bondType = bondTypes[bondType]
		if bondType.__class__ != BondType:
			raise Exception('Invalid parameter "bondType".')
		self.bondType = bondType

################################################################################

class Structure:
	"""
	Represent the chemical structure of a single resonance form of a chemical
	species as a graph.
	"""	

	def __init__(self, atoms=[], bonds=[]):
		self.atoms = atoms
		self.bonds = bonds
		self.updateGraph()
	
	def updateGraph(self):
		self.graph = {}
		for atom in self.atoms:
			self.graph[atom] = []
		for bond in self.bonds:
			for atom in bond.atoms:
				self.graph[atom].append(bond)
	
	def fromOBMol(self, obmol):
		"""
		Convert an OpenBabel OBMol object to a Structure object.
		"""
		
		self.atoms = []; self.bonds = []; self.graph = {}
		
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
			
			atom = Atom(elements[number], electronStates[electron])
			self.atoms.append(atom)
			
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
					
					bond = Bond(self.atoms[i], self.atoms[j], bondTypes[order])
					self.bonds.append(bond)
		
		# Create the graph from the atom and bond lists
		self.updateGraph()
		
	def toOBMol(self):
		"""
		Convert a Structure object to an OpenBabel OBMol object.
		"""
		obmol = openbabel.OBMol()
		for atom in self.atoms:
			a = obmol.NewAtom()
			a.SetAtomicNum(atom.element.number)
		for bond in self.bonds:
			index1 = self.atoms.index(bond.atoms[0])
			index2 = self.atoms.index(bond.atoms[1])
			order = bond.bondType.order
			if order == 1.5: order = 5
			obmol.AddBond(index1+1, index2+1, order)
		
		return obmol
	
	def toInChI(self):
		mol = pybel.Molecule(self.toOBMol())
		return mol.write('inchi').strip()
		
################################################################################

class Species:
	"""
	Represent a chemical species (including all of its resonance forms).
	"""	

	numSpecies = 0

	def __init__(self, label='', structure=None, reactive=True):
		Species.numSpecies += 1
		self.id = Species.numSpecies
		self.label = label
		self.structure = structure
		self.reactive = reactive
		
	def toInChI(self):
		return self.structure.toInChI()
	
	def __str__(self):
		return self.label + '(' + str(self.id) + ')'
		
################################################################################

if __name__ == '__main__':
	
	print ''
	
	print 'Elements available:'
	for key, element in elements.iteritems():
		print '\t' + str(key) + ' ' + element.symbol
	print ''
	
	print 'Free electron states available:'
	for key, electronState in electronStates.iteritems():
		print '\t' + electronState.label
	print ''
	
	print 'Bond types available:'
	for key, bondType in bondTypes.iteritems():
		print '\t' + str(key) + ' ' + bondType.label
	print ''
	