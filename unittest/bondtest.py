#!/usr/bin/python
# -*- coding: utf-8 -*-
 
import unittest

import sys
sys.path.append('../source')

from rmg.chem import *

################################################################################

class BondTypeEquivalencyCheck(unittest.TestCase):                          

	def testSelfEquivalence(self):
		"""
		Check that all atom types are equivalent to themselves.
		"""
		for key, value in bondTypes.iteritems():
			self.assertTrue(value.equivalent(value))
			
#	def testSpecialEquivalence(self):
#		"""
#		Check that certain special cases are equivalent.
#		"""
#		self.assertTrue(bondTypes['D'].equivalent(bondTypes['Dcis']))
#		self.assertTrue(bondTypes['D'].equivalent(bondTypes['Dtrans']))
#		self.assertTrue(bondTypes['Dcis'].equivalent(bondTypes['D']))
#		self.assertTrue(bondTypes['Dtrans'].equivalent(bondTypes['D']))
						
		
################################################################################

class BondCheck(unittest.TestCase):                          

	def testEquivalency(self):
		"""
		Checks involving bond equivalency.
		"""
		bond1 = Bond([None, None], ['S', 'D', 'T'])
		bond2 = Bond([None, None], ['D'])
		self.assertTrue(bond1.equivalent(bond2))
		self.assertTrue(bond2.equivalent(bond1))
						
	def testCopy(self):
		"""
		Check that atom copy operations are successful.
		"""
		for label, bondType in bondTypes.iteritems():
			bond1 = Bond([None, None], bondType)
			bond2 = bond1.copy()
			self.assertTrue(bond2 is not None)
			self.assertTrue(bond1 is not bond2)
			self.assertTrue(bond1.bondType == bond2.bondType)
			self.assertTrue(bond1.atoms == bond2.atoms)

	def testBondTypeIdentification(self):
		"""
		Check the Bond.isSingle(), etc. functions.
		"""
		for label1, bondType in bondTypes.iteritems():
			bond = Bond([None, None], bondType)
			if bondType.order == 1:
				self.assertTrue(bond.isSingle())
				self.assertFalse(bond.isDouble())
				self.assertFalse(bond.isTriple())
				self.assertFalse(bond.isBenzene())
			elif bondType.order == 2:
				self.assertFalse(bond.isSingle())
				self.assertTrue(bond.isDouble())
				self.assertFalse(bond.isTriple())
				self.assertFalse(bond.isBenzene())
			elif bondType.order == 3:
				self.assertFalse(bond.isSingle())
				self.assertFalse(bond.isDouble())
				self.assertTrue(bond.isTriple())
				self.assertFalse(bond.isBenzene())
			elif bondType.order == 1.5:
				self.assertFalse(bond.isSingle())
				self.assertFalse(bond.isDouble())
				self.assertFalse(bond.isTriple())
				self.assertTrue(bond.isBenzene())

	def testBondTypeManipulation(self):
		"""
		Check the bond type manipulation functions.
		"""
		for label1, bondType in bondTypes.iteritems():
			bond = Bond([None, None], bondType)
	
				
			if bondType.order == 1:
				self.assertTrue(bond.canIncreaseOrder())
				self.assertFalse(bond.canDecreaseOrder())
			elif bondType.order == 2:
				self.assertTrue(bond.canIncreaseOrder())
				self.assertTrue(bond.canDecreaseOrder())
			elif bondType.order == 3:
				self.assertFalse(bond.canIncreaseOrder())
				self.assertTrue(bond.canDecreaseOrder())
			else:
				self.assertFalse(bond.canIncreaseOrder())
				self.assertFalse(bond.canDecreaseOrder())
			
			if bond.canIncreaseOrder():
				order = bondType.order
				bond.increaseOrder()
				self.assertTrue(bond.bondType.order == order + 1)
				bond.decreaseOrder()
				self.assertTrue(bond.bondType.order == order)
				
			if bond.canDecreaseOrder():
				order = bondType.order
				bond.decreaseOrder()
				self.assertTrue(bond.bondType.order == order - 1)
				bond.increaseOrder()
				self.assertTrue(bond.bondType.order == order)
			

################################################################################

if __name__ == '__main__':
	unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )
	
	