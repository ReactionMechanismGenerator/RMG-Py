#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
from unittest import TestCase

from rmgpy import settings
from species import Species
from rmgpy.data.solvation import *
from rmgpy.molecule.molecule import Molecule

###################################################

class TestSoluteDatabase(TestCase):
	
	def runTest(self):
		pass

	def testSoluteGeneration(self):
		
		self.database = SoluteDatabase()
		self.database.load(os.path.join(settings['database.directory'], 'thermo'))
		
		self.testCases = [
		
		['cyclopentane',		'C1CCCC1',		0.10,	0,	0.263,	2.477,	0],
		['methylcyclopentane',	'C1(CCCC1)C',	0.10,	0,	0.225,	2.816,	0],
		['cyclohexane',			'C1CCCCC1',		0.10,	0,	0.305,	2.964,	0]
		
		]
		
		for name, smiles, S, B, E, L, A in self.testCases:
			species = Species(molecule=[Molecule(SMILES=smiles)])
			soluteData = self.database.getSoluteData(Species(molecule=[species.molecule[0]]))
			print self.assertEqual(soluteData.S, S)
			print self.assertEqual(soluteData.B, B)
			print self.assertEqual(soluteData.E, E)
			print self.assertEqual(soluteData.L, L)
			print self.assertEqual(soluteData.A, A)
	
#####################################################

if __name__ == '__main__':
	myTest = TestSoluteDatabase()
	myTest.testSoluteGeneration()