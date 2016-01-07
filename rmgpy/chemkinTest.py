import unittest
import os
from chemkin import *
###################################################

class ChemkinTest(unittest.TestCase):

	def testReadTemplateReactionFamilyForMinimalExample(self):
		"""
		This example is mainly to test if family info can be correctly 
		parsed from comments like '!Template reaction: R_Recombination'.
		"""
		folder = os.path.join(os.getcwd(),'rmgpy/tools/data/chemkin/chemkin_py')
		
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
		folder = os.path.join(os.getcwd(),'rmgpy/tools/data/chemkin/chemkin_py')
		
		chemkinPath = os.path.join(folder, 'pdd', 'chem.inp')
		dictionaryPath = os.path.join(folder,'pdd', 'species_dictionary.txt')

		# loadChemkinFile
		species, reactions = loadChemkinFile(chemkinPath, dictionaryPath) 

		reaction1 = reactions[0]
		self.assertEqual(reaction1.family, "H_Abstraction")

		reaction2 = reactions[1]
		self.assertEqual(reaction2.family, "H_Abstraction")