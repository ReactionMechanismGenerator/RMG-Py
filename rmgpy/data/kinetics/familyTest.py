import unittest
from rmgpy.data.kinetics.database import KineticsDatabase
import os.path
from rmgpy.molecule.group import Group
###################################################

class TestFamily(unittest.TestCase):


    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        # Set up a dummy database
        dir_path = os.path.dirname(os.path.realpath(__file__))
        self.database = KineticsDatabase()
        self.database.loadFamilies(os.path.join(dir_path,"family_test_data"), families=['intra_H_migration'])
        self.family = self.database.families['intra_H_migration']

    def testGetBackboneRoots(self):
        """
        Test the getBackboneRoots() function
        """
        backbones = self.family.getBackboneRoots()
        self.assertEquals(backbones[0].label, "RnH")

    def testGetEndRoots(self):
        """
        Test the getEndRoots() function
        """
        ends = self.family.getEndRoots()
        self.assertEquals(len(ends), 2)
        self.assertIn(self.family.groups.entries["Y_rad_out"], ends)
        self.assertIn(self.family.groups.entries["XH_out"], ends)

    def testGetTopLevelGroups(self):
        """
        Test the getTopLevelGroups() function
        """
        topGroups = self.family.getTopLevelGroups(self.family.groups.entries["RnH"])
        self.assertEquals(len(topGroups), 2)
        self.assertIn(self.family.groups.entries["R5Hall"], topGroups)
        self.assertIn(self.family.groups.entries["R6Hall"], topGroups)