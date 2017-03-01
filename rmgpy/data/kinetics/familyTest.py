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

    def testMergeGroups(self):
        """
        Test the mergeGroups() function
        """
        #basic test of merging a backbone and end group
        backbone1 = Group().fromAdjacencyList("""
1 *1 R!H u1 {2,S}
2 *4 R!H u0 {1,S} {3,S}
3 *6 R!H u0 {2,S} {4,S}
4 *5 R!H u0 {3,S} {5,S}
5 *2 R!H u0 {4,S} {6,S}
6 *3 H   u0 {5,S}
""")

        end1 = Group().fromAdjacencyList("""
1 *2 Cs u0 {2,S} {3,S}
2 *3 H  u0 {1,S}
3    S  u0 {1,S}
""")
        desiredMerge1 = Group().fromAdjacencyList("""
1 *1 R!H u1 {2,S}
2 *4 R!H u0 {1,S} {3,S}
3 *6 R!H u0 {2,S} {4,S}
4 *5 R!H u0 {3,S} {5,S}
5 *2 Cs  u0 {4,S} {6,S} {7,S}
6 *3 H   u0 {5,S}
7    S   u0 {5,S}
""")

        mergedGroup = self.family.mergeGroups(backbone1, end1)
        self.assertTrue(mergedGroup.isIdentical(desiredMerge1))

        #test it works when there is a cyclical structure to the backbone

        backbone2 = Group().fromAdjacencyList("""
1 *1 R!H u1 {2,S} {4,S}
2 *4 R!H u0 {1,S} {3,S}
3 *2 R!H u0 {2,S} {4,S}
4 *3 R!H u0 {3,S} {1,S}
""")

        end2 = Group().fromAdjacencyList("""
1 *2 Os u0 {2,S}
2 *3 Cs u0 {1,S}
""")
        desiredMerge2 = Group().fromAdjacencyList("""
1 *1 R!H u1 {2,S} {4,S}
2 *4 R!H u0 {1,S} {3,S}
3 *2 Os u0 {2,S} {4,S}
4 *3 Cs u0 {3,S} {1,S}
""")
        mergedGroup = self.family.mergeGroups(backbone2, end2)
        self.assertTrue(mergedGroup.isIdentical(desiredMerge2))
