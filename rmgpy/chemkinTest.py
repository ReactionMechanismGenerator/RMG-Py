import unittest
import os
from chemkin import *
import rmgpy
###################################################

class ChemkinTest(unittest.TestCase):

	def testReadTemplateReactionFamilyForMinimalExample(self):
		"""
		This example is mainly to test if family info can be correctly 
		parsed from comments like '!Template reaction: R_Recombination'.
		"""
		folder = os.path.join(os.path.dirname(rmgpy.__file__),'test_data/chemkin/chemkin_py')
		
		chemkinPath = os.path.join(folder, 'minimal', 'chem.inp')
		dictionaryPath = os.path.join(folder,'minimal', 'species_dictionary.txt')

		# loadChemkinFile
		species, reactions = loadChemkinFile(chemkinPath, dictionaryPath) 

		reaction1 = reactions[0]
		self.assertEqual(reaction1.family, "R_Recombination")

		reaction2 = reactions[1]
		self.assertEqual(reaction2.family, "H_Abstraction")
	
	def testReadTemplateReactionFamilyForPDDExample(self):
		"""
		This example is mainly to ensure comments like 
		'! Kinetics were estimated in this direction instead 
		of the reverse because:' or '! This direction matched 
		an entry in H_Abstraction, the other was just an estimate.'
		won't interfere reaction family info retrival.
		"""
		folder = os.path.join(os.path.dirname(rmgpy.__file__),'test_data/chemkin/chemkin_py')
		
		chemkinPath = os.path.join(folder, 'pdd', 'chem.inp')
		dictionaryPath = os.path.join(folder,'pdd', 'species_dictionary.txt')

		# loadChemkinFile
		species, reactions = loadChemkinFile(chemkinPath, dictionaryPath) 

		reaction1 = reactions[0]
		self.assertEqual(reaction1.family, "H_Abstraction")

		reaction2 = reactions[1]
		self.assertEqual(reaction2.family, "H_Abstraction")
		
	def testTransportDataReadAndWrite(self):
		"""
		Test that we can write to chemkin and recreate the same transport object
		"""
		from rmgpy.species import Species
		from rmgpy.molecule import Molecule
		from rmgpy.transport import TransportData
		
		Ar = Species(label="Ar", transportData=TransportData(shapeIndex=0, epsilon=(1134.93,'J/mol'), sigma=(3.33,'angstrom'), dipoleMoment=(0,'De'), polarizability=(0,'angstrom^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""))
		
		Ar_write = Species(label="Ar")
		folder = os.path.join(os.path.dirname(rmgpy.__file__),'test_data')
		
		tempTransportPath = os.path.join(folder, 'tran_temp.dat')
		
		saveTransportFile(tempTransportPath, [Ar])
		speciesDict = {'Ar': Ar_write}
		loadTransportFile(tempTransportPath, speciesDict)
		self.assertEqual(repr(Ar),repr(Ar_write))
		
		os.remove(tempTransportPath)
		

	def testUseChemkinNames(self):
		"""
		Test that the official chemkin names are used as labels for the created Species objects.
		"""

		folder = os.path.join(os.path.dirname(rmgpy.__file__),'test_data/chemkin/chemkin_py')
		
		chemkinPath = os.path.join(folder, 'minimal', 'chem.inp')
		dictionaryPath = os.path.join(folder,'minimal', 'species_dictionary.txt')

		# loadChemkinFile
		species, reactions = loadChemkinFile(chemkinPath, dictionaryPath, useChemkinNames=True) 

		expected = [
			'Ar',
			'He',
			'Ne',
			'N2',
			'ethane',
			'CH3',
			'C2H5',
			'C'
		]

		for spc, label in zip(species, expected):
			self.assertEqual(spc.label, label)

	def testReactantN2IsReactiveAndGetsRightSpeciesIdentifier(self):
		"""
		Test that after loading chemkin files, species such as N2, which is in the default 
		inert list of RMG, should be treated as reactive species and given right species
		Identifier when it's reacting in reactions.
		"""
		folder = os.path.join(os.path.dirname(rmgpy.__file__),'test_data/chemkin/chemkin_py')
		
		chemkinPath = os.path.join(folder, 'NC', 'chem.inp')
		dictionaryPath = os.path.join(folder,'NC', 'species_dictionary.txt')

		# loadChemkinFile
		species, reactions = loadChemkinFile(chemkinPath, dictionaryPath, useChemkinNames=True) 

		for n2 in species:
		    if n2.label == 'N2':
		        break
		self.assertTrue(n2.reactive)

		self.assertEqual(getSpeciesIdentifier(n2), 'N2(35)') 

	

