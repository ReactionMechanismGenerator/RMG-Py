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

import unittest

import sys
sys.path.append('../source')

import cPickle

import rmg.main as main
import rmg.data as data
import rmg.species as species
import rmg.reaction as reaction
import rmg.structure as structure
import rmg.thermo as thermo

def loadThermoDatabases(databasePath):
	"""
	Create and load the thermodynamics databases.
	"""
	import os.path
	databasePath += '/'

	# Create and load thermo databases
	species.thermoDatabase = species.ThermoDatabaseSet()
	species.thermoDatabase.load(databasePath)

	# Create and load forbidden structures
	species.forbiddenStructures = data.Dictionary()
	species.forbiddenStructures.load(os.path.join(databasePath, 'forbiddenStructure.txt'))
	species.forbiddenStructures.toStructure()

def loadKineticsDatabases(databasePath, only_families=False):
	"""
	Create and load the kinetics databases (reaction families).
	If only_families is a list like ['H_Abstraction'] then only families in this
	list will be loaded.
	"""
	reaction.kineticsDatabase = reaction.ReactionFamilySet()
	reaction.kineticsDatabase.load(databasePath, only_families=only_families)

################################################################################

class RestartCheck(unittest.TestCase):
	"""
	Tests the ability of various RMG classes to be pickled and unpickled
	successfully.
	"""

	def testAtomRestart(self):
		"""
		Restart tests for the :class:`rmg.chem.Atom` class.
		"""
		import rmg.chem as chem

		# Create
		atom0 = chem.Atom(atomType='O', electronState='1', charge=0, label='*1')
		# Pickle
		f = open('test.pkl', 'wb'); cPickle.dump(atom0, f); f.close()
		# Unpickle
		f = open('test.pkl', 'rb'); atom = cPickle.load(f); f.close()
		# Compare
		self.assertTrue(atom0.equals(atom))

	def testBondRestart(self):
		"""
		Restart tests for the :class:`rmg.chem.Bond` class.
		"""
		import rmg.chem as chem

		# Create
		atom1 = chem.Atom(atomType='C', electronState='0', charge=0, label='*1')
		atom2 = chem.Atom(atomType='O', electronState='1', charge=0, label='*2')
		bond0 = chem.Bond(atoms=[atom1, atom2], bondType='D')
		# Pickle
		f = open('test.pkl', 'wb'); cPickle.dump(bond0, f); f.close()
		# Unpickle
		f = open('test.pkl', 'rb'); bond = cPickle.load(f); f.close()
		# Compare
		self.assertTrue(bond0.equals(bond))

	def testStructureRestart(self):
		"""
		Restart tests for the :class:`rmg.structure.Structure` class.
		"""
		import rmg.structure as structure

		# Create
		structure1 = structure.Structure(SMILES='C=CC=C[CH]CC=O')
		structure1.symmetryNumber = 2
		# Pickle
		f = open('test.pkl', 'wb'); cPickle.dump(structure1, f); f.close()
		# Unpickle
		f = open('test.pkl', 'rb'); structure2 = cPickle.load(f); f.close()
		# Compare
		self.assertTrue(structure1.symmetryNumber == structure2.symmetryNumber)

	def testSpeciesRestart(self):
		"""
		Restart tests for the :class:`rmg.species.Species` class.
		"""
		import rmg.structure as structure
		import rmg.species as species

		# Create
		smiles = 'C=CC=C[CH]C'
		struct0 = structure.Structure(SMILES=smiles)
		spec0 = species.makeNewSpecies(struct0, label=smiles, reactive=True)

		# Pickle
		f = open('test.pkl', 'wb'); cPickle.dump(spec0, f); f.close()
		# Unpickle
		f = open('test.pkl', 'rb'); spec = cPickle.load(f); f.close()
		# Compare
		self.assertTrue(spec0.id == spec.id == 1)
		self.assertTrue(spec0.label == spec.label == smiles)
		self.assertTrue(spec0.reactive == spec.reactive == True)
		self.assertTrue(len(spec0.structure) == len(spec.structure) == 3)
		for struct1 in spec0.structure:
			match = False
			for struct2 in spec.structure:
				if struct1.isIsomorphic(struct2): match = True
			self.assertTrue(match)
		self.assertTrue(spec0.thermoData.equals(spec.thermoData))
		#self.assertTrue(spec0.thermoSnapshot == spec.thermoSnapshot)
		self.assertTrue(spec0.spectralData == spec.spectralData)
		self.assertTrue(spec0.lennardJones == spec.lennardJones)

	def testThermoRestart(self):
		"""
		Restart tests for the :class:`rmg.thermo.ThermoGAData`,
		:class:`rmg.thermo.ThermoWilhoitData`, and 
		:class:`rmg.thermo.ThermoNASAData` classes.
		"""
		import rmg.thermo as thermo
		
		# Create
		thermoGAData0 = thermo.ThermoGAData(H298=-3.08, S298=64.27, 
			Cp=[12.28, 14.34, 16.30, 18.05, 20.92, 23.08, 26.39], comment='acetyl')
		thermoWilhoitData0 = thermo.convertGAtoWilhoit(thermoGAData0, 6, 1, False)
		thermoNASAData0 = thermo.convertWilhoitToNASA(thermoWilhoitData0)
		# Pickle
		f = open('test.pkl', 'wb');
		cPickle.dump(thermoGAData0, f)
		cPickle.dump(thermoWilhoitData0, f)
		cPickle.dump(thermoNASAData0, f)
		f.close()
		# Unpickle
		f = open('test.pkl', 'rb')
		thermoGAData = cPickle.load(f)
		thermoWilhoitData = cPickle.load(f)
		thermoNASAData = cPickle.load(f)
		f.close()
		# Compare
		self.assertTrue(thermoGAData0.equals(thermoGAData))
		self.assertTrue(thermoWilhoitData0.equals(thermoWilhoitData))
		self.assertTrue(thermoNASAData0.equals(thermoNASAData))

	def testKineticsRestart(self):
		"""
		Restart tests for the :class:`rmg.kinetics.ArrheniusKinetics` and
		:class:`rmg.kinetics.ArrheniusEPKinetics` classes.
		"""
		import rmg.kinetics as kinetics

		# Create
		aKinetics0 = kinetics.ArrheniusKinetics(A=1.29e6, Ea=24.57, n=2.0)
		aepKinetics0 = kinetics.ArrheniusEPKinetics(A=1.29e6, E0=24.57, n=2.0, alpha=1.5)
		# Pickle
		f = open('test.pkl', 'wb');
		cPickle.dump(aKinetics0, f)
		cPickle.dump(aepKinetics0, f)
		f.close()
		# Unpickle
		f = open('test.pkl', 'rb')
		aKinetics = cPickle.load(f)
		aepKinetics = cPickle.load(f)
		f.close()
		# Compare
		self.assertTrue(aKinetics0.equals(aKinetics))
		self.assertTrue(aepKinetics0.equals(aepKinetics))

	def testReactionRestart(self):
		"""
		Restart tests for the :class:`rmg.reaction.Reaction` class.
		"""
		self.fail('Not yet implemented -- but it needs to be!')

################################################################################

if __name__ == '__main__':

	# Show debug messages as databases are loading
	main.initializeLog(10)

	# Load databases
	databasePath = '../data/RMG_database'
	loadThermoDatabases(databasePath)
	#loadKineticsDatabases(databasePath)
	
	# Show info messages during tests
	main.initializeLog(20)

	# Conduct tests
	unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )