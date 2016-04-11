import unittest
from rmgpy.species import Species
import rmgpy.tools.data as dt
import numpy

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

    def testRemoveDataSet(self):

        #define some variables for the tests
        acetylene=Species().fromSMILES('C#C')
        x1=dt.GenericData(label='time', data=[1,2,3], species=None, reaction=None, units='s', index=None)
        y1=dt.GenericData(label='oldModel', data=[4,5,6], species=acetylene, reaction=None, units='molFrac', index=None)
        x2=dt.GenericData(label='time', data=[7,8,9], species=None, reaction=None, units=None, index=None)
        y2=dt.GenericData(label='newModel', data=[10,11,12], species=None, reaction=None, units=None, index=None)
        x3=dt.GenericData(label='time', data=[1,2,3], species=None, reaction=None, units='s', index=None)
        y3=dt.GenericData(label='expt', data=[4,5,6], species=acetylene, reaction=None, units='molFrac', index=None)

        #Check to see that class catches the incorrect units
        test3=dt.ComparisonBundle(title="RemoveTest", xDataList=[x2, x1, x3], yDataList=[y2, y1, y3])

        #tests that the removeSet Function works correctly for yLabel
        test3.removeDataSet(yLabel="oldModel")
        self.assertEqual(test3.xDataList, [x2,x3])
        self.assertEqual(test3.yDataList, [y2,y3])

        #test the renumbering of indicies works
        for xData, yData in zip(test3.xDataList, test3.yDataList):
            self.assertEqual(xData.index, test3.xDataList.index(xData))
            self.assertEqual(yData.index, test3.yDataList.index(yData))

        #tests that the removeSet Function works correctly for index
        test3.removeDataSet(index=1)
        self.assertEquals(test3.xDataList, [x2])
        self.assertEquals(test3.yDataList, [y2])

    def testCheckConsistent(self):

        #define some variables for the tests
        acetylene=Species().fromSMILES('C#C')
        methane=Species().fromSMILES('C')
        x1=dt.GenericData(label='time', data=[1,2,3], species=None, reaction=None, units='s', index=None)
        y1=dt.GenericData(label='oldModel', data=[4,5,6], species=acetylene, reaction=None, units='molFrac', index=None)
        x2=dt.GenericData(label='time', data=[7,8,9], species=None, reaction=None, units=None, index=None)
        y2=dt.GenericData(label='newModel', data=[10,11,12], species=None, reaction=None, units=None, index=None)
        x4=dt.GenericData(label='time', data=[1,2,3], species=None, reaction=None, units='ms', index=None)
        y4=dt.GenericData(label='expt', data=[4,5,6], species=methane, reaction=None, units='molFrac', index=None)

        #Check to see that class catches the incorrect units
        test4=dt.ComparisonBundle(title="CheckTest", xDataList=[x2, x1], yDataList=[y2, y1])
        self.assertRaises(AssertionError, test4.addDataSet,x4,y4)

        #Check to see that the class catches incorrect species
        test4=dt.ComparisonBundle(title="CheckTest", xDataList=[x2, x1], yDataList=[y2, y1])
        y4.species=acetylene
        self.assertRaises(AssertionError, test4.addDataSet,x4,y4)

    def testMakeLog(self):

        acetylene=Species().fromSMILES('C#C')
        x1=dt.GenericData(label='time', data=[1,2,3], species=None, reaction=None, units='s', index=None)
        y1=dt.GenericData(label='oldModel', data=[0,1,10], species=acetylene, reaction=None, units='molFrac', index=None)
        x2=dt.GenericData(label='time', data=[7,8,9], species=None, reaction=None, units=None, index=None)
        y2=dt.GenericData(label='newModel', data=[-1,1,100], species=None, reaction=None, units='molFrac', index=None)

        test5=dt.ComparisonBundle(title="CheckLog", xDataList=[x2, x1], yDataList=[y2, y1])
        test5.makeLog(10)

        #What should be outputted (should remove first data point and print correct logs)
        x2New=numpy.array([8,9])
        y2New=numpy.array([0,1])
        x1New=numpy.array([2,3])
        y1New=numpy.array([0,2])

        #Check values
        self.assertAlmostEqual(test5.xDataList[0].data.all(), x2New.all())
        self.assertAlmostEqual(test5.yDataList[0].data.all(), y2New.all())
        self.assertAlmostEqual(test5.xDataList[1].data.all(), x1New.all())
        self.assertAlmostEqual(test5.yDataList[1].data.all(), y1New.all())

        #Check units
        self.assertEqual(test5.yUnits, 'log10 molFrac')
        self.assertEqual(test5.yDataList[0].units, 'log10 molFrac')
        self.assertEqual(test5.yDataList[1].units, 'log10 molFrac')