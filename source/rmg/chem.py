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
	Loads entries into a dictionary of elements. The dictionary created by
	this function is always available at :data:`rmg.chem.elements`.
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
	elements[0] = elements['R'] = Element(0, '', 'R', pq.Quantity(0.0, 'g/mol'), 0)
	return elements
	
# The dictionary of elements, accessible by atomic number, symbol, or name.
elements = loadElements()

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
		
class Atom:
	"""
	Represent an atom in a chemical species. Each atom is defined by an 
	`element` (stored internally as an :class:`Element` object), a `electronState`
	(stored internally as an :class:`ElectronState` object), and a numeric `charge`.
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
		Set the element that this atom represents. The *element* parameter can
		be an :class:`Element` object, an integer representing the atomic 
		number, or a string containing the element symbol or name. In all 
		cases *element* will be converted to and stored as an :class:`Element` 
		object.
		"""
		if element.__class__ == str or element.__class__ == unicode:
			element = elements[element]
		if element.__class__ != Element:
			raise Exception('Invalid parameter "element".')
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
		if electronState.__class__ != ElectronState:
			raise Exception('Invalid parameter "electronState".')
		self.electronState = electronState
	
	def equivalent(self, other):
		"""
		Return :data:`True` if `self` and `other` are equivalent atoms, i.e. 
		they have the same element and electronic state.
		"""
		return self.element == other.element and self.electronState == other.electronState

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
		if bondType.__class__ != BondType:
			raise Exception('Invalid parameter "bondType".')
		self.bondType = bondType
	
	def equivalent(self, other):
		"""
		Return :data:`True` if `self` and `other` are equivalent bonds, i.e. 
		they have the same bond type.
		"""
		return self.bondType == other.bondType

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
		
	def addBond(self, atom1, atom2, bond):
		"""
		Add `bond` to the graph as an edge connecting atoms `atom1` and
		`atom2`, which must already be present in the graph.
		"""
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
	
	def isIsomorphic(self, other):
		"""
		Returns :data:`True` if two graphs are isomorphic and :data:`False`
		otherwise. Uses the VF2 algorithm of Vento and Foggia.
		"""
		return graph.VF2_isomorphic(self.graph, other.graph, False)
		
	def isSubgraphIsomorphic(self, other):
		"""
		Returns :data:`True` if `other` is subgraph isomorphic and :data:`False`
		otherwise. Uses the VF2 algorithm of Vento and Foggia.
		"""
		return graph.VF2_isomorphic(self.graph, other.graph, True)
	
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

class Structure(ChemGraph):
	"""
	A representation of a chemical species using a graph data structure. The 
	vertices represent atoms, while the edges represent bonds.
	"""

	def __init__(self, atoms=[], bonds=[]):
		self.initialize(atoms, bonds)
		
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
		
	def fromAdjacencyList(self, adjlist):
		"""
		Convert a string adjacency list `adjlist` into a chemical structure.
		"""
		
		atoms = []; bonds = []; atomdict = {}
				
		lines = adjlist.splitlines()
		label = lines[0]
		for line in lines[1:]:
			
			data = line.split()
			
			# First item is index for atom
			aid = int(data[0])
			
			# Next is the element or functional group element
			element = data[1]
				
			# Next is the electron state
			elecState = data[2].upper()
			
			# Create a new atom based on the above information
			atom = Atom(element, elecState, 0)
				
			# Add the atom to the list
			atoms.append(atom)
			atomdict[aid] = atom
			
			# Process list of bonds
			for datum in data[3:]:
				aid2, comma, btype = datum[1:-1].partition(',')
				aid2 = int(aid2)
				if aid2 in atomdict:
					bond = Bond([atomdict[aid], atomdict[aid2]], btype)
					bonds.append(bond)
			
		# Create species structure
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
			
			atom = Atom(elements[number], electronStates[electron])
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
			a.SetAtomicNum(atom.element.number)
		for bond in bonds:
			index1 = atoms.index(bond.atoms[0])
			index2 = atoms.index(bond.atoms[1])
			order = bond.bondType.order
			if order == 1.5: order = 5
			obmol.AddBond(index1+1, index2+1, order)
		
		return obmol
	
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
		if structure is None:
			self.structure = Structure()
		else:
			self.structure = structure
		self.reactive = reactive
		
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
	
	#print ''
	
	#print 'Elements available:'
	#for key, element in elements.iteritems():
		#print '\t' + str(key) + ' ' + element.symbol
	#print ''
	
	#print 'Free electron states available:'
	#for key, electronState in electronStates.iteritems():
		#print '\t' + electronState.label
	#print ''
	
	#print 'Bond types available:'
	#for key, bondType in bondTypes.iteritems():
		#print '\t' + str(key) + ' ' + bondType.label
	#print ''
	
	cml1 = """
<molecule>
<atomArray>
<atom id="a1" elementType="C" />
<atom id="a2" elementType="C" />
<atom id="a3" elementType="C" />
<atom id="a4" elementType="C" />
<atom id="a5" elementType="C" />
<atom id="a6" elementType="Cl" />
</atomArray>
<bondArray>
<bond atomRefs2="a1 a2" order="2" />
<bond atomRefs2="a2 a3" order="1" />
<bond atomRefs2="a3 a4" order="2" />
<bond atomRefs2="a4 a5" order="1" />
<bond atomRefs2="a5 a6" order="1" />
</bondArray>
</molecule>
	"""
	structure1 = Structure()
	structure1.fromCML(cml1)
	print structure1.toInChI(), structure1.toSMILES()
	
	cml2 = """
<molecule>
<atomArray>
<atom id="a3" elementType="C" />
<atom id="a6" elementType="C" />
<atom id="a4" elementType="C" />
<atom id="a1" elementType="F" />
<atom id="a5" elementType="C" />
<atom id="a2" elementType="C" />
</atomArray>
<bondArray>
<bond atomRefs2="a2 a3" order="1" />
<bond atomRefs2="a4 a5" order="1" />
<bond atomRefs2="a1 a2" order="1" />
<bond atomRefs2="a5 a6" order="2" />
<bond atomRefs2="a3 a4" order="2" />
</bondArray>
</molecule>
	"""
	structure2 = Structure()
	structure2.fromCML(cml2)
	print structure2.toInChI(), structure2.toSMILES()
	
	print structure1.isIsomorphic(structure2)
	
