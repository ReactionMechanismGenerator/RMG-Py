#!/usr/bin/python
# -*- coding: utf-8 -*-
 
import unittest

import sys
sys.path.append('../source')

from rmg.structure import *

################################################################################

class StructureCheck(unittest.TestCase):                          

	def testIsomorphism(self):
		"""
		Check the graph isomorphism functions.
		"""
		structure1 = Structure()
		structure1.fromSMILES('C=CC=C[CH]C')
		
		structure2 = Structure()
		structure2.fromSMILES('C[CH]C=CC=C')
		
		self.assertTrue(structure1.isIsomorphic(structure2))
		self.assertTrue(structure2.isIsomorphic(structure1))
		
	def testSubgraphIsomorphism(self):
	
		structure1 = Structure()
		structure1.fromSMILES('C=CC=C[CH]C')
		
		structure2 = Structure()
		structure2.fromAdjacencyList("""
		1 C 0 {2,D}
		2 C 0 {1,D}
		""")
		
		self.assertTrue(structure1.isSubgraphIsomorphic(structure2))
		match, map21, map12 = structure1.findSubgraphIsomorphisms(structure2)
		self.assertTrue(match)
		self.assertTrue(len(map21) == len(map12) == 4)
		for mapA, mapB in zip(map21, map12):
			self.assertTrue(len(mapA) == len(mapB) == min(len(structure1.atoms()), len(structure2.atoms())))
			for key, value in mapA.iteritems():
				self.assertTrue(value in mapB)
				self.assertTrue(key is mapB[value])
				self.assertTrue(key in structure1.atoms())
				self.assertTrue(value in structure2.atoms())
			for key, value in mapB.iteritems():
				self.assertTrue(value in mapA)
				self.assertTrue(key is mapA[value])
				self.assertTrue(key in structure2.atoms())
				self.assertTrue(value in structure1.atoms())
				
			
################################################################################

if __name__ == '__main__':
	unittest.main()