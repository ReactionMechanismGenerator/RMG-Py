#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
from unittest import TestCase

from rmgpy import settings
from rmgpy.species import Species
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
		
		# from RMG-Java test run
		['octane', 'CCCCCCCC', 0.127,	0.085, 0.04,	3.766,	0.0030	, 1.2358], 	
		['water', 'O', 0.524,	0.378,	0.309,	0.802,	0.348,	0.1673],
		
		
		# from RMG-Java (Jalan et al supplementary data):
		#['1,2-ethanediol', 'C(CO)O', 0.823, 0.685, 0.327, 2.572, 0.693, None]
		['1-decanol', 'C(CCCCCCCCC)O', 0.449, 0.385, 0.205, 5.614, 0.348, None],
		#['acetylaldehyde', 'CC=O', 0.622,	0.423,	0.171,	1.415,	0.0030,	0.4061],
		['acetic acid', 'C(C)(=O)O', 0.508, 0.411, 0.152, 1.873, 0.591, None],
		#['anthracene', 'C2=CC=CC3=CC1=CC=CC=C1C=C23', 1.181, 0.181, 1.648, 7.316, 0.003, None]
		#['dimethyl ether']
		
		# from Abraham:
		#['cyclopentane',		'C1CCCC1',		0.10,	0,	0.263,	2.477,	0, None],		#['methylcyclopentane',	'C1(CCCC1)C',	0.10,	0,	0.225,	2.816,	0],
		#['cyclohexane',			'C1CCCCC1',		0.10,	0,	0.305,	2.964,	0]
		]
		
		for name, smiles, S, B, E, L, A, V in self.testCases:
			species = Species(molecule=[Molecule(SMILES=smiles)])
			soluteData = self.database.getSoluteData(Species(molecule=[species.molecule[0]]))
			print name, soluteData
			print self.assertAlmostEqual(soluteData.S, S)
			print self.assertAlmostEqual(soluteData.B, B)
			print self.assertAlmostEqual(soluteData.E, E)
			print self.assertAlmostEqual(soluteData.L, L)
			print self.assertAlmostEqual(soluteData.A, A)
	
#####################################################

if __name__ == '__main__':
	myTest = TestSoluteDatabase()
	myTest.testSoluteGeneration()