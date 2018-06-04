#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2017 Prof. William H. Green (whgreen@mit.edu), 
#   Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

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
        1 *3 C  u1 {2,D} {3,S}
        2    C  u0 {1,D}
        3 *5 Cd u0 {1,S} {4,D}
        4    C  u0 {3,D}
        """)
        )        
        
        entry2 = Entry(
            item= Group().fromAdjacencyList(
        """
        1 *3 C  u1 {2,D} {3,S}
        2 *5 C  u0 {1,D}
        3    Cd u0 {1,S} {4,D}
        4    C  u0 {3,D}
        """)
        )
        
        entry3 = Entry(
            item = Group().fromAdjacencyList(
        """
        1 *3 C  u1 {2,D} {3,S}
        2    C  u0 {1,D}
        3    Cd u0 {1,S} {4,D}
        4    C  u0 {3,D}
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
        1 *1 R!H u1
        """)
        )        
        
        entry2 = Entry(
            item= Group().fromAdjacencyList(
        """
        1 *1 Cb u1
        """)
        )
        self.assertTrue(self.database.matchNodeToNode(entry1,entry1))
        self.assertFalse(self.database.matchNodeToNode(entry1,entry2))



################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))

