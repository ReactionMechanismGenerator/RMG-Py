#!/usr/bin/python
# -*- coding: utf-8 -*-
 
import unittest

import sys
sys.path.append('../source')

import rmg.thermo
from rmg.structure import *
from rmg.species import *

################################################################################

class SpeciesCheck(unittest.TestCase):                          

	def testResonance(self):
		"""
		Check that the resonance form generator is working correctly.
		"""
		species = Species()
		species.fromSMILES('C=CC=CC=CC=C[CH]C')
		species.getResonanceIsomers()
		self.assertTrue(len(species.structure) == 5, "Found %d structures, expected 5"%len(species.structure) )
		for structure in species.structure:
			self.assertTrue(structure.getFormula() == 'C10H13')
		
		species = Species()
		species.fromSMILES('C=CC=CC=CC=C[CH]C=C')
		species.getResonanceIsomers()
		self.assertTrue(len(species.structure) == 3)
		for structure in species.structure:
			self.assertTrue(structure.getFormula() == 'C11H13')
	
	def testMakeNewSpecies(self):
	
		structure1 = Structure()
		structure1.fromSMILES('C=CC=C[CH]C')
		
		structure2 = Structure()
		structure2.fromSMILES('C[CH]C=CC=C')
		
		species1, isNew = makeNewSpecies(structure1)
		species2, isNew = makeNewSpecies(structure2)
		
		self.assertTrue(species1 is species2)
				
################################################################################

if __name__ == '__main__':
	unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )