#!/usr/bin/python
# -*- coding: utf-8 -*-
 
import unittest

import sys
sys.path.append('.')

from rmgpy.data.base import *
from rmgpy.chem.pattern import MoleculePattern

thermoDatabase = 'database/output/RMG_Database/thermo_groups'

################################################################################

class DatabaseCheck(unittest.TestCase):                          

	def testDatabaseLoad(self):
		"""
		Check the database load functions.
		"""
		dictstr = thermoDatabase + '/Group_Dictionary.txt'
		treestr = thermoDatabase + '/Group_Tree.txt'
		libstr = thermoDatabase + '/Group_Library.txt'
		
		database = Database()
		database.load(dictstr, treestr, libstr)
		
		# All nodes in library should be in tree and dictionary
		# All nodes in tree should be in dictionary
		for node in database.library:
			self.assertTrue(node in database.tree.parent)
			self.assertTrue(node in database.tree.children)
			self.assertTrue(node in database.dictionary)
		for node in database.tree.parent:
			self.assertTrue(node in database.tree.children)
			self.assertTrue(node in database.dictionary)
		for node in database.tree.children:
			self.assertTrue(node in database.tree.parent)
			self.assertTrue(node in database.dictionary)
		
		# All parents in tree should be in tree
		for node in database.tree.parent:
			parentNode = database.tree.parent[node]
			if parentNode is not None:
				self.assertTrue(parentNode in database.tree.parent)
				self.assertTrue(parentNode in database.tree.children)
		
		# All children in tree should be in tree
		for node in database.tree.children:
			for childNode in database.tree.children[node]:
				if childNode is not None:
					self.assertTrue(childNode in database.tree.parent)
					self.assertTrue(childNode in database.tree.children)
		
		# All values in dictionary should be chemical structures
		for node in database.dictionary:
			self.assertTrue(isinstance(database.dictionary[node], MoleculePattern))

################################################################################

if __name__ == '__main__':
	unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )
	
	