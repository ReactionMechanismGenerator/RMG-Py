"""
This scripts runs tests on the database 
"""
import os.path
import logging
import unittest
from external.wip import work_in_progress
from rmgpy import settings
from rmgpy.data.rmg import RMGDatabase

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
            self.assertTrue(family_name in familyNumberNodes, "{family} family does not exist in test's dictionary of families. You may need to update the unit test.".format(family=family_name))
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
                        
    @work_in_progress                    
    def test_kinetics_checkGroupsFoundInTree(self):
        """
        This test checks whether groups are found in the tree.
        """
        for family_name, family in self.database.kinetics.families.iteritems():
            for nodeName, nodeGroup in family.groups.entries.iteritems():
                ascendParent = nodeGroup
                while ascendParent not in family.groups.top and ascendParent not in family.forwardTemplate.products:
                    child = ascendParent
                    ascendParent = ascendParent.parent
                    self.assertTrue(ascendParent is not None, "Group {group} in {family} family was found in the tree without a proper parent.".format(group=child,family=family_name))
                    self.assertTrue(child in ascendParent.children, "Group {group} in {family} family was found in the tree without a proper parent.".format(group=nodeName,family=family_name))


if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
