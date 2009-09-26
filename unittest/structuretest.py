#!/usr/bin/python
# -*- coding: utf-8 -*-
 
import unittest

import sys
sys.path.append('../source')

from rmg.structure import *

################################################################################

class StructureCheck(unittest.TestCase):  

	def testSSSR(self):
		"""
		Check the graph's Smallest Set of Smallest Rings function
		"""
		structure = Structure()
		structure.fromSMILES('C(CC1C(C(CCCCCCCC)C1c1ccccc1)c1ccccc1)CCCCCC')
		#http://cactus.nci.nih.gov/chemical/structure/C(CC1C(C(CCCCCCCC)C1c1ccccc1)c1ccccc1)CCCCCC/image
		sssr = structure.getSmallestSetOfSmallestRings()
		self.assertEqual( len(sssr), 3)

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
		self.assertTrue(len(map21) == len(map12) == 4, "len(map21) = %d, len(map12) = %d, should both = 4"%(len(map21),len(map12)))
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

	def testSubgraphIsomorphismAgain(self):

		structure1 = Structure()
		#structure1.fromSMILES('C=CC=CCC')
		structure1.fromAdjacencyList("""
		1 * C 0 {2,D} {7,S} {8,S}
		2 C 0 {1,D} {3,S} {9,S}
		3 C 0 {2,S} {4,D} {10,S}
		4 C 0 {3,D} {5,S} {11,S}
		5 C 0 {4,S} {6,S} {12,S} {13,S}
		6 C 0 {5,S} {14,S} {15,S} {16,S}
		7 H 0 {1,S}
		8 H 0 {1,S}
		9 H 0 {2,S}
		10 H 0 {3,S}
		11 H 0 {4,S}
		12 H 0 {5,S}
		13 H 0 {5,S}
		14 H 0 {6,S}
		15 H 0 {6,S}
		16 H 0 {6,S}
		""")

		structure2 = Structure()
		structure2.fromAdjacencyList("""
		1 * C 0 {2,D} {3,S} {4,S}
		2 C 0 {1,D}
		3 H 0 {1,S}
		4 H 0 {1,S}
		""")

		labeled1 = structure1.getLabeledAtoms().values()[0]
		labeled2 = structure2.getLabeledAtoms().values()[0]

		map21_0 = {labeled2: labeled1}; map12_0 = {labeled1: labeled2}
		self.assertTrue(structure1.isSubgraphIsomorphic(structure2, map12_0, map21_0))

		map21_0 = {labeled2: labeled1}; map12_0 = {labeled1: labeled2}
		match, map21, map12 = structure1.findSubgraphIsomorphisms(structure2, map12_0, map21_0)
		self.assertTrue(match)
		self.assertTrue(len(map21) == len(map12) == 2)
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

	def testIsInCycle(self):

		# ethane
		struct = Structure()
		struct.fromSMILES('CC')
		for atom in struct.atoms():
			self.assertFalse(struct.isAtomInCycle(atom))
		for bond in struct.bonds():
			self.assertFalse(struct.isBondInCycle(bond))

		# cyclohexane
		struct = Structure()
		struct.fromInChI('InChI=1/C6H12/c1-2-4-6-5-3-1/h1-6H2')
		for atom in struct.atoms():
			if atom.atomType.element.symbol == 'H':
				self.assertFalse(struct.isAtomInCycle(atom))
			elif atom.atomType.element.symbol == 'C':
				self.assertTrue(struct.isAtomInCycle(atom))
		for bond in struct.bonds():
			atom1, atom2 = bond.atoms
			if atom1.atomType.element.symbol == 'C' and atom2.atomType.element.symbol == 'C':
				self.assertTrue(struct.isBondInCycle(bond))
			else:
				self.assertFalse(struct.isBondInCycle(bond))
			
#		struct = Structure()
#		struct.fromInChI('InChI=1/C6H12/c1-2-4-6-5-3-1/h1-6H2')
#		symmetryNumber = struct.calculateCyclicSymmetryNumber()
#		self.assertEqual(symmetryNumber, 12)

	def testAtomSymmetryNumber(self):

		# methane
		struct = Structure()
		struct.fromSMILES('C')
		symmetryNumber = 1
		for atom in struct.atoms():
			symmetryNumber *= struct.calculateAtomSymmetryNumber(atom)
		self.assertEqual(symmetryNumber, 12)
		
		# methyl
		struct = Structure()
		struct.fromSMILES('[CH3]')
		symmetryNumber = 1
		for atom in struct.atoms():
			symmetryNumber *= struct.calculateAtomSymmetryNumber(atom)
		self.assertEqual(symmetryNumber, 6)
		
		# ethane
		struct = Structure()
		struct.fromSMILES('CC')
		symmetryNumber = 1
		for atom in struct.atoms():
			symmetryNumber *= struct.calculateAtomSymmetryNumber(atom)
		self.assertEqual(symmetryNumber, 9)

		# propane
		struct = Structure()
		struct.fromSMILES('CCC')
		symmetryNumber = 1
		for atom in struct.atoms():
			symmetryNumber *= struct.calculateAtomSymmetryNumber(atom)
		self.assertEqual(symmetryNumber, 18)
		
		# isobutane
		struct = Structure()
		struct.fromSMILES('CC(C)C')
		symmetryNumber = 1
		for atom in struct.atoms():
			symmetryNumber *= struct.calculateAtomSymmetryNumber(atom)
		self.assertEqual(symmetryNumber, 81)

	def testBondSymmetryNumber(self):

		# ethane
		struct = Structure()
		struct.fromSMILES('CC')
		symmetryNumber = 1
		for bond in struct.bonds():
			symmetryNumber *= struct.calculateBondSymmetryNumber(bond)
		self.assertEqual(symmetryNumber, 2)
		
		# propane
		struct = Structure()
		struct.fromSMILES('CCC')
		symmetryNumber = 1
		for bond in struct.bonds():
			symmetryNumber *= struct.calculateBondSymmetryNumber(bond)
		self.assertEqual(symmetryNumber, 1)

		# butane
		struct = Structure()
		struct.fromSMILES('CCCC')
		symmetryNumber = 1
		for bond in struct.bonds():
			symmetryNumber *= struct.calculateBondSymmetryNumber(bond)
		self.assertEqual(symmetryNumber, 2)

		# ethylene
		struct = Structure()
		struct.fromSMILES('C=C')
		symmetryNumber = 1
		for bond in struct.bonds():
			symmetryNumber *= struct.calculateBondSymmetryNumber(bond)
		self.assertEqual(symmetryNumber, 2)

		# acetylene
		struct = Structure()
		struct.fromSMILES('C#C')
		symmetryNumber = 1
		for bond in struct.bonds():
			symmetryNumber *= struct.calculateBondSymmetryNumber(bond)
		self.assertEqual(symmetryNumber, 2)

	def testAxisSymmetryNumber(self):

		# allene
		struct = Structure()
		struct.fromSMILES('C=C=C')
		symmetryNumber = struct.calculateAxisSymmetryNumber()
		self.assertEqual(symmetryNumber, 2)

		# cumulene
		struct = Structure()
		struct.fromSMILES('C=C=C=C')
		symmetryNumber = struct.calculateAxisSymmetryNumber()
		self.assertEqual(symmetryNumber, 2)

#	def testCyclicSymmetryNumber(self):
#
#		# cyclohexane
#		struct = Structure()
#		struct.fromInChI('InChI=1/C6H12/c1-2-4-6-5-3-1/h1-6H2')
#		symmetryNumber = struct.calculateCyclicSymmetryNumber()
#		self.assertEqual(symmetryNumber, 12)

	def testSymmetryNumber(self):

		# ethane
		struct = Structure()
		struct.fromSMILES('CC')
		struct.calculateSymmetryNumber()
		self.assertEqual(struct.symmetryNumber, 18)

		
################################################################################
from timeit import Timer



if __name__ == '__main__':
	
	startup = """gc.enable() # enable garbage collection in timeit
import sys
sys.path.append('../source')
from rmg.structure import Structure
structure1 = Structure()
structure1.fromSMILES('C=CC=C[CH]C')
structure2 = Structure()
structure2.fromSMILES('C[CH]C=CC=C')
"""
	test = "structure1.isIsomorphic(structure2)"
	print "Timing isIsomorphic:"
	t = Timer(test,startup)
	print "  took %.3f milliseconds"%min(t.repeat(repeat=10,number=1000))
	
	
	unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )