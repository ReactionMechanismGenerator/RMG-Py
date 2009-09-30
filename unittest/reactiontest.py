#!/usr/bin/python
# -*- coding: utf-8 -*-
 
import unittest

import sys
sys.path.append('../source')

from rmg.structure import Structure
from rmg.species import makeNewSpecies
from rmg.reaction import *
import rmg.reaction as reaction

# Run this whether being run as __main__ or called by other unit test suite:



	
	
################################################################################

class ReactionCheck(unittest.TestCase):                          

	def testMakeNewReaction(self):
		"""Create a new reaction and check you can identify the reverse"""
		structure1 = Structure()
		structure1.fromAdjacencyList("""
		1 C 0 {2,D} {7,S} {8,S}
		2 C 0 {1,D} {3,S} {9,S}
		3 C 0 {2,S} {4,D} {10,S}
		4 C 0 {3,D} {5,S} {11,S}
		5 *1 C 0 {4,S} {6,S} {12,S} {13,S}
		6 C 0 {5,S} {14,S} {15,S} {16,S}
		7 H 0 {1,S}
		8 H 0 {1,S}
		9 H 0 {2,S}
		10 H 0 {3,S}
		11 H 0 {4,S}
		12 *2 H 0 {5,S}
		13 H 0 {5,S}
		14 H 0 {6,S}
		15 H 0 {6,S}
		16 H 0 {6,S}
		""")
		
		structure2 = Structure()
		structure2.fromAdjacencyList("""
		1 *3 H 1
		""")
		
		structure3 = Structure()
		structure3.fromAdjacencyList("""
		1 C 0 {2,D} {7,S} {8,S}
		2 C 0 {1,D} {3,S} {9,S}
		3 C 0 {2,S} {4,D} {10,S}
		4 C 0 {3,D} {5,S} {11,S}
		5 *3 C 1 {4,S} {6,S} {12,S}
		6 C 0 {5,S} {13,S} {14,S} {15,S}
		7 H 0 {1,S}
		8 H 0 {1,S}
		9 H 0 {2,S}
		10 H 0 {3,S}
		11 H 0 {4,S}
		12 H 0 {5,S}
		13 H 0 {6,S}
		14 H 0 {6,S}
		15 H 0 {6,S}
		""")
		
		structure4 = Structure()
		structure4.fromAdjacencyList("""
		1 *1 H 0 {2,S}
		2 *2 H 0 {1,S}
		""")
		
		C6H10 = makeNewSpecies(structure1)
		H = makeNewSpecies(structure2)
		C6H9 = makeNewSpecies(structure3)
		H2 = makeNewSpecies(structure4)
		
		# wipe the reaction list
		reaction.reactionList=[]
		
		reaction1, isNew = makeNewReaction([C6H9, H2], [C6H10, H], \
			[C6H9.structure[0], H2.structure[0]], \
			[C6H10.structure[0], H.structure[0]], \
			None)
		self.assertFalse(reaction1 is None)
		self.assertTrue(isNew)
		
		reaction2, isNew = makeNewReaction([C6H10, H], [C6H9, H2], \
			[C6H10.structure[0], H.structure[0]], \
			[C6H9.structure[0], H2.structure[0]], \
			None)
		
		self.assertTrue(reaction1 is reaction2)
		self.assertFalse(isNew)
				
		
class ReactionSetCheck(unittest.TestCase): 
	
	def loadDatabase(self,only_families=[]):
		# Load database
		databasePath = '../data/RMG_database'
		if only_families: 
			print "Only loading reaction families %s"%only_families
		reaction.kineticsDatabase = reaction.ReactionFamilySet()
		reaction.kineticsDatabase.load(databasePath, only_families=only_families)

	
	def testEthaneExtrusionFromEther(self):
		"""A reaction that was giving trouble in a diesel simulation"""
		self.loadDatabase(only_families=['1,3_Insertion_ROR'])
		structure1 = Structure()
		structure1.fromSMILES("c1c([CH]C(C=CCCCCCC)OOCc2c3ccccc3ccc2)cccc1")
		species1 = makeNewSpecies(structure1)
		
		structure2 = Structure()
		structure2.fromSMILES("c1ccc(c2ccccc12)COO")
		species2 = makeNewSpecies(structure2)
		
		structure3 = Structure()
		structure3.fromSMILES("c1cc(C=C[C]=CCCCCCC)ccc1")
		species3 = makeNewSpecies(structure3)
		
		# wipe the reaction list
		reaction.reactionList=[]
		#print "FORWARD"
		rxns = reaction.kineticsDatabase.getReactions([species2,species3])
		#for rxn in rxns:
		#	print 'Reaction family:',rxn.family
		#	print 'Reaction:',rxn
		#	print 'Kinetics:',rxn.kinetics
		#	print 'bestKinetics:',rxn.bestKinetics
		#	print
		self.assertEqual(len(rxns),18, "Was expecting to make 18 reactions for %s + %s"%(species2,species3))
			
		# wipe the reaction list
		reaction.reactionList=[]
		#print "REVERSE"
		rxns = reaction.kineticsDatabase.getReactions([species1])
		#for rxn in rxns:
		#	print 'Reaction family:',rxn.family
		#	print 'Reaction:',rxn
		#	print 'Kinetics:',rxn.kinetics
		#	print 'bestKinetics:',rxn.bestKinetics
		#	print
		self.assertEqual(len(rxns),1, "Was expecting to make 1 reaction for %s"%(species1))
		rxn = rxns[0]
		self.assertEqual(len(rxn.reactants),2,"Reaction wasn't bimolecular")
		self.assertTrue(species2 in rxn.reactants, "Didn't make %s"%species2)
		self.assertTrue(species3 in rxn.reactants, "Didn't make %s"%species3)

	def testIntraHmigration(self):
		"""An intra-H migration reaction"""
		self.loadDatabase(only_families=['intra_H_migration'])
		structure1 = Structure()
		structure1.fromSMILES("[CH](CCCc1ccccc1)CCCCCC")
		species1 = makeNewSpecies(structure1)

		structure2 = Structure()
		structure2.fromSMILES("C(CCCC[CH]c1ccccc1)CCCC")
		species2 = makeNewSpecies(structure2)
			
		# wipe the reaction list
		reaction.reactionList=[]

		rxns = reaction.kineticsDatabase.getReactions([species1])
		#for rxn in rxns:
		#	print 'Reaction family:',rxn.family
		#	print 'Reaction:',rxn
		#	print 'Kinetics:',rxn.kinetics
		#	print 'bestKinetics:',rxn.bestKinetics
		#	print
		
		all_products = []
		for rxn in rxns:
			self.assertEqual(str(rxn.family),'Intra H migration',"Was trying to test 'Intra H migration' but made a reaction from family %s"%rxn.family)
			self.assertEqual(len(rxn.reactants),1,"Reaction %s wasn't unimolecular"%rxn)
			all_products.extend(rxn.products)
		self.assertTrue(species2 in all_products, "None of the reactions made %s"%(species2))

		#self.assertEqual(len(rxns),1, "Was expecting to make 1 reaction for %s"%(species1))



################################################################################

if __name__ == '__main__':
	unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )