#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import unittest
from external.wip import work_in_progress

from rmgpy.data.base import Entry, Database
from rmgpy.molecule import Group

################################################################################

class TestBaseDatabase(unittest.TestCase):
    """
    Contains unit tests for the base class of rmgpy.data.
    """
    
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        # Set up a dummy database
        self.database = Database()


    def testMatchNodeToStructure(self):
        """
        Test that the MatchNodeToStructure family works properly.
        """
        entry1 = Entry(
            item = Group().fromAdjacencyList(
        """
        1 *3 C  1 {2,D} {3,S}
        2    C  0 {1,D}
        3 *5 Cd 0 {1,S} {4,D}
        4    C  0 {3,D}
        """)
        )        
        
        entry2 = Entry(
            item= Group().fromAdjacencyList(
        """
        1 *3 C  1 {2,D} {3,S}
        2 *5 C  0 {1,D}
        3    Cd 0 {1,S} {4,D}
        4    C  0 {3,D}
        """)
        )
        
        entry3 = Entry(
            item = Group().fromAdjacencyList(
        """
        1 *3 C  1 {2,D} {3,S}
        2    C  0 {1,D}
        3    Cd 0 {1,S} {4,D}
        4    C  0 {3,D}
        """)
        )
        # The group should match to itself
        self.assertTrue(self.database.matchNodeToStructure(entry1,entry1.item,atoms=entry1.item.getLabeledAtoms()))  
        
        # These groups should not match each other
        self.assertFalse(self.database.matchNodeToStructure(entry1,entry2.item,atoms=entry2.item.getLabeledAtoms()))  

        # entry1 contains more labels than entry3, therefore cannot be matched by entry3
        self.assertFalse(self.database.matchNodeToStructure(entry3,entry1.item,atoms=entry1.item.getLabeledAtoms()))

        # entry3 contains fewer labels than entry1, therefore it can be matched
        self.assertTrue(self.database.matchNodeToStructure(entry1,entry3.item,atoms=entry3.item.getLabeledAtoms()))
        
    def testMatchNodeToNode(self):
        """
        Test that nodes can match other nodes.
        """
        entry1 = Entry(
            item = Group().fromAdjacencyList(
        """
        1 *1 R!H 1
        """)
        )        
        
        entry2 = Entry(
            item= Group().fromAdjacencyList(
        """
        1 *1 Cb 1
        """)
        )
        self.assertTrue(self.database.matchNodeToNode(entry1,entry1))
        self.assertFalse(self.database.matchNodeToNode(entry1,entry2))
################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))

