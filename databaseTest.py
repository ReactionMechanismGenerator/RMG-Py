"""
This scripts runs tests on the database 
"""
import os.path
import logging
import unittest
from external.wip import work_in_progress
from rmgpy import settings
from rmgpy.data.rmg import RMGDatabase
from copy import copy, deepcopy

class TestDatabase(unittest.TestCase):
    """
    Contains unit tests for the database for rigorous error checking.
    """
    @classmethod
    def setUpClass(cls):
        """
        Load the database before running the tests.
        """
        databaseDirectory = settings['database.directory']
        cls.database = RMGDatabase()
        cls.database.load(databaseDirectory, kineticsFamilies='all')
        
    def test_kinetics_checkCorrectNumberofNodesInRules(self):
        """
        This test ensures that each rate rule contains the proper number of nodes according to the family it originates.
        """
        
        for family_name, family in self.database.kinetics.families.iteritems():
            expectedNumberNodes = len(family.getRootTemplate())
            for label, entries in family.rules.entries.iteritems():
                for entry in entries:
                    nodes = label.split(';')
                    self.assertEqual(len(nodes), expectedNumberNodes, "Wrong number of groups or semicolons in family {family} rule {entry}.  Should be {num_nodes}".format(family=family_name,entry=entry,num_nodes=expectedNumberNodes))
                    
    def test_kinetics_checkNodesInRulesFoundInGroups(self):
        """
        This test ensures that each rate rule contains nodes that exist in the groups
        """
        for family_name, family in self.database.kinetics.families.iteritems():
            for label, entries in family.rules.entries.iteritems():
                for entry in entries:
                    nodes = label.split(';')
                    for node in nodes:
                        self.assertTrue(node in family.groups.entries, "In {family} family, no group definition found for label {label} in rule {entry}".format(family=family_name, label=node, entry=entry))
                                        
    def test_kinetics_checkGroupsFoundInTree(self):
        """
        This test checks whether groups are found in the tree.
        """
        for family_name, family in self.database.kinetics.families.iteritems():
            for nodeName, nodeGroup in family.groups.entries.iteritems():
                ascendParent = nodeGroup
                # Check whether the node has proper parents unless it is the top reactant or product node
                while ascendParent not in family.groups.top and ascendParent not in family.forwardTemplate.products:
                    child = ascendParent
                    ascendParent = ascendParent.parent
                    self.assertTrue(ascendParent is not None, "Group {group} in {family} family was found in the tree without a proper parent.".format(group=child,family=family_name))
                    self.assertTrue(child in ascendParent.children, "Group {group} in {family} family was found in the tree without a proper parent.".format(group=nodeName,family=family_name))
                    
    def test_kinetics_checkGroupsNonidentical(self):
        """
        This test checks that the groups are non-identical.
        """
        from rmgpy.data.base import Database
        for family_name, originalFamily in self.database.kinetics.families.iteritems():
            family = Database()
            family.entries = originalFamily.groups.entries
            entriesCopy = copy(family.entries)
            for nodeName, nodeGroup in family.entries.iteritems():
                del entriesCopy[nodeName]
                for nodeNameOther, nodeGroupOther in entriesCopy.iteritems():
                    self.assertFalse(family.matchNodeToNode(nodeGroup, nodeGroupOther), "Group {group} in {family} family was found to be identical to group {groupOther}".format(group=nodeName, family=family_name, groupOther=nodeNameOther))
    
    @work_in_progress
    def test_kinetics_checkChildParentRelationships(self):
        """
        This test checks that groups' parent-child relationships are correct in the database.
        """
        from rmgpy.data.base import Database
        for family_name, originalFamily in self.database.kinetics.families.iteritems():
            family = Database()
            family.entries = originalFamily.groups.entries
            entriesCopy = copy(family.entries)
            for nodeName, nodeGroup in family.entries.iteritems():
                ascendParent = nodeGroup
                # Check whether the node has proper parents unless it is the top reactant or product node
                while ascendParent not in originalFamily.groups.top and ascendParent not in originalFamily.forwardTemplate.products:
                    child = ascendParent
                    ascendParent = ascendParent.parent
                    # The parent should be more general than the child
                    self.assertTrue(family.matchNodeToChild(ascendParent,child), 
                                    "In {family} family, group {parent} is not a proper parent of its child {child}.".format(family=family_name, parent=ascendParent, child=child))
                    
                    

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
