import unittest
from rmgpy.species import Species
import rmgpy.tools.data as dt

class ComparisonBundleTest(unittest.TestCase):
    def testCheckAndMakeConsistent(self):
        """
        Tests the checkAndMakeConsistent function works. It is called during initialization of ComparisonBundle
        """

        acetylene=Species().fromSMILES('C#C')
        x1=dt.GenericData(label='time', data=[1,2,3], species=None, reaction=None, units='s', index=None)
        y1=dt.GenericData(label='oldModel', data=[4,5,6], species=acetylene, reaction=None, units='molFrac', index=None)
        x2=dt.GenericData(label='time', data=[7,8,9], species=None, reaction=None, units=None, index=None)
        y2=dt.GenericData(label='newModel', data=[10,11,12], species=None, reaction=None, units=None, index=None)

        test1=dt.ComparisonBundle(title="test", xDataList=[x2, x1], yDataList=[y2, y1])

        #tests that the head attribute can copy slave attributes
        self.assertEqual(test1.yUnits, 'molFrac')
        self.assertEqual(test1.xUnits, 's')
        self.assertTrue(test1.species.isIsomorphic(acetylene))

        #test that head attribute gets distributed to xData and yData
        for x in test1.xDataList:
            self.assertEqual(x.units, 's')
        for y in test1.yDataList:
            self.assertEqual(y.units,'molFrac')
            self.assertTrue(y.species.isIsomorphic(acetylene))

    def testAddDataSet(self):

        #define some variables for the tests
        acetylene=Species().fromSMILES('C#C')
        methane=Species().fromSMILES('C')
        x1=dt.GenericData(label='time', data=[1,2,3], species=None, reaction=None, units='s', index=None)
        y1=dt.GenericData(label='oldModel', data=[4,5,6], species=acetylene, reaction=None, units='molFrac', index=None)
        x2=dt.GenericData(label='time', data=[7,8,9], species=None, reaction=None, units=None, index=None)
        y2=dt.GenericData(label='newModel', data=[10,11,12], species=None, reaction=None, units=None, index=None)
        x3=dt.GenericData(label='time', data=[1,2,3], species=None, reaction=None, units='s', index=None)
        y3=dt.GenericData(label='expt', data=[4,5,6], species=acetylene, reaction=None, units='molFrac', index=None)

        #tests that addSet Function works correctly
        test2=dt.ComparisonBundle(title="AddTest", xDataList=[x2, x1], yDataList=[y2, y1])
        test2.addDataSet(x3,y3)
        self.assertTrue(test2.xDataList[-1]==x3)
        self.assertTrue(test2.yDataList[-1]==y3)