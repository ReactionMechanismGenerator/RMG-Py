#!/usr/bin/python
# -*- coding: utf-8 -*-
 
import unittest

import sys
sys.path.append('../source')

from rmg.chem import *

################################################################################

class AtomTypeEquivalencyCheck(unittest.TestCase):                          

	def testSelfEquivalence(self):
		"""
		Check that all atom types are equivalent to themselves.
		"""
		for key, value in atomTypes.iteritems():
			self.assertTrue(value.equivalent(value))
			
	def testElementEquivalence(self):
		"""
		Check that all element atom types are equivalent to all other atom
		types with the same element.
		"""
		
		for key1, value1 in atomTypes.iteritems():
			if value1.element is not None:
				if value1.label == value1.element.symbol:
					for key2, value2 in atomTypes.iteritems():
						if value2.element is not None:
							if value2.element.symbol == value1.element.symbol and value2 is not value1:
								self.assertTrue(value1.equivalent(value2))
								self.assertTrue(value2.equivalent(value1))
	
	def testGenericEquivalence(self):
		"""
		Check that the generic atom type 'R' is equivalent to all other atom
		types and that the generic atom type 'R!H' is equivalent to all other
		atom types except 'H'.
		"""
		
		for key, value in atomTypes.iteritems():
			if value is not atomTypes['R'] and value is not atomTypes['R!H']:
				self.assertTrue(value.equivalent(atomTypes['R']))
				self.assertTrue(atomTypes['R'].equivalent(value))
				if value is not atomTypes['H']:
					self.assertTrue(value.equivalent(atomTypes['R!H']))
					self.assertTrue(atomTypes['R!H'].equivalent(value))
		self.assertFalse(atomTypes['H'].equivalent(atomTypes['R!H']))
		self.assertFalse(atomTypes['R!H'].equivalent(atomTypes['H']))
		
	
	def testIsSpecificCaseOf(self):
		"""
		Check the isSpecificCaseOf() method.
		"""
		
		for name, a in atomTypes.iteritems():
			
			# everything is in R and itself
			self.assertTrue(a.isSpecificCaseOf(atomTypes['R']) )
			self.assertTrue(a.isSpecificCaseOf(a) )
			
			# R is in nothing but R
			if name is not 'R':
				self.assertFalse(atomTypes['R'].isSpecificCaseOf(a), "R should not be specific case of %s"%name)
				
			# H is not in R!H
			if name is 'H': 
				self.assertFalse(a.isSpecificCaseOf(atomTypes['R!H']))
			elif name is not 'R': # but everything else is, except R
				self.assertTrue(a.isSpecificCaseOf(atomTypes['R!H']))
				
			# R!H is in nothing but R and itself
			if name not in ['R','R!H']:
				self.assertFalse(atomTypes['R!H'].isSpecificCaseOf(a))
			
			# maybe there are others I should test
			if a.element and a.element.symbol is 'C':
				self.assertTrue(a.isSpecificCaseOf(atomTypes['C']))
				if name != 'C': # then something specific like 'Cd'
					self.assertFalse(atomTypes['C'].isSpecificCaseOf(a))
					

################################################################################

class ElectronStateEquivalencyCheck(unittest.TestCase):                          

	def testSelfEquivalence(self):
		"""
		Check that all atom types are equivalent to themselves.
		"""
		for key, value in electronStates.iteritems():
			self.assertTrue(value.equivalent(value))
			
	def testSpecialEquivalence(self):
		"""
		Check that certain special cases are equivalent.
		"""
		self.assertTrue(electronStates['2'].equivalent(electronStates['2S']))
		self.assertTrue(electronStates['2'].equivalent(electronStates['2T']))
		self.assertTrue(electronStates['2S'].equivalent(electronStates['2']))
		self.assertTrue(electronStates['2T'].equivalent(electronStates['2']))
						
		
################################################################################

class AtomCheck(unittest.TestCase):                          

	def testEquivalency(self):
		"""
		Checks involving atom equivalency.
		"""
		atom1 = Atom(['Cs', 'Cd'], '0', 0, '1*')
		atom2 = Atom(['C'], '0', 0, '2*')
		self.assertTrue(atom1.equivalent(atom2))
		self.assertTrue(atom2.equivalent(atom1))
						
	def testCopy(self):
		"""
		Check that atom copy operations are successful.
		"""
		for atomType in atomTypes:
			for electronState in electronStates:
				atom1 = Atom(atomType, electronState, -1, '1*')
				atom2 = atom1.copy()
				self.assertTrue(atom2 is not None)
				self.assertTrue(atom1 is not atom2)
				self.assertTrue(atom1.atomType == atom2.atomType)
				self.assertTrue(atom1.electronState == atom2.electronState)
				self.assertTrue(atom1.charge == atom2.charge)
				self.assertTrue(atom1.label == atom2.label)
				
	def testCenter(self):
		"""
		Check the Atom.isCenter() function.
		"""
		self.assertTrue(Atom('C', '0', 0, '1*').isCenter())
		self.assertTrue(Atom('C', '0', 0, '*').isCenter())
		self.assertFalse(Atom('C', '0', 0, '').isCenter())

	def testElement(self):
		"""
		Check the Atom.isElement() function.
		"""
		for label1, atomType in atomTypes.iteritems():
			for label2, electronState in electronStates.iteritems():
				atom = Atom(atomType, electronState, 0, '1*')
				if atomType.element is not None:
					self.assertTrue(atom.isElement(atomType.element.symbol),
						"%s.isElement('%s') returned False"%(atom,atomType.element.symbol))
					if atomType.element.symbol == 'C':
						self.assertTrue(atom.isCarbon())
					if atomType.element.symbol == 'O':
						self.assertTrue(atom.isOxygen())
					if atomType.element.symbol == 'H':
						self.assertTrue(atom.isHydrogen())
						self.assertFalse(atom.isNonHydrogen())
					else:
						self.assertFalse(atom.isHydrogen())
						self.assertTrue(atom.isNonHydrogen())
				

	def testFreeElectron(self):
		"""
		Check the Atom.isElement() function.
		"""
		for label1, atomType in atomTypes.iteritems():
			for label2, electronState in electronStates.iteritems():
				atom = Atom(atomType, electronState, 0, '1*')
				self.assertTrue(atom.getFreeElectronCount() == electronState.order)
				
				if electronState.order == 0:
					self.assertTrue(atom.canIncreaseFreeElectron())
					self.assertFalse(atom.canDecreaseFreeElectron())
				elif electronState.order >= 4:
					self.assertFalse(atom.canIncreaseFreeElectron())
					self.assertTrue(atom.canDecreaseFreeElectron())
				else:
					self.assertTrue(atom.canIncreaseFreeElectron())
					self.assertTrue(atom.canDecreaseFreeElectron())
				
				if atom.canIncreaseFreeElectron():
					order = electronState.order
					atom.increaseFreeElectron()
					self.assertTrue(atom.electronState.order == order + 1)
					atom.decreaseFreeElectron()
					self.assertTrue(atom.electronState.order == order)
					
				if atom.canDecreaseFreeElectron():
					order = electronState.order
					atom.decreaseFreeElectron()
					self.assertTrue(atom.electronState.order == order - 1)
					atom.increaseFreeElectron()
					self.assertTrue(atom.electronState.order == order)
				
################################################################################

if __name__ == '__main__':
	AtomTypeEquivalencyCheck('testIsSpecificCaseOf').debug()
	
	unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )