#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2018 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

import unittest
from external.wip import work_in_progress

from rmgpy.molecule.group import ActionError, GroupAtom, GroupBond, Group
from rmgpy.molecule.atomtype import atomTypes
from rmgpy.molecule import Molecule
import element as elements

################################################################################

class TestGroupAtom(unittest.TestCase):
    """
    Contains unit tests of the GroupAtom class.
    """

    def setUp(self):
        """
        A method called before each unit test in this class.
        """
        self.atom = GroupAtom(atomType=[atomTypes['Cd']], radicalElectrons=[1], charge=[0], label='*1', lonePairs=[0])
    
    def testApplyActionBreakBond(self):
        """
        Test the GroupAtom.applyAction() method for a BREAK_BOND action.
        """
        action = ['BREAK_BOND', '*1', 1, '*2']
        for label, atomType in atomTypes.iteritems():
            atom0 = GroupAtom(atomType=[atomType], radicalElectrons=[1], charge=[0], label='*1', lonePairs=[0])
            atom = atom0.copy()
            try:
                atom.applyAction(action)
                self.assertEqual(len(atom.atomType), len(atomType.breakBond))
                for a in atomType.breakBond:
                    self.assertTrue(a in atom.atomType)
                self.assertEqual(atom0.radicalElectrons, atom.radicalElectrons)
                self.assertEqual(atom0.charge, atom.charge)
                self.assertEqual(atom0.label, atom.label)
                self.assertEqual(atom0.label, atom.label)
                self.assertEqual(atom0.lonePairs, atom.lonePairs)
            except ActionError:
                self.assertEqual(len(atomType.breakBond), 0)
    
    def testApplyActionFormBond(self):
        """
        Test the GroupAtom.applyAction() method for a FORM_BOND action.
        """
        action = ['FORM_BOND', '*1', 1, '*2']
        for label, atomType in atomTypes.iteritems():
            atom0 = GroupAtom(atomType=[atomType], radicalElectrons=[1], charge=[0], label='*1', lonePairs=[0])
            atom = atom0.copy()
            try:
                atom.applyAction(action)
                self.assertEqual(len(atom.atomType), len(atomType.formBond))
                for a in atomType.formBond:
                    self.assertTrue(a in atom.atomType)
                self.assertEqual(atom0.radicalElectrons, atom.radicalElectrons)
                self.assertEqual(atom0.charge, atom.charge)
                self.assertEqual(atom0.label, atom.label)
                self.assertEqual(atom0.lonePairs, atom.lonePairs)
            except ActionError:
                self.assertEqual(len(atomType.formBond), 0)
    
    def testApplyActionIncrementBond(self):
        """
        Test the GroupAtom.applyAction() method for a CHANGE_BOND action.
        """
        action = ['CHANGE_BOND', '*1', 1, '*2']
        for label, atomType in atomTypes.iteritems():
            atom0 = GroupAtom(atomType=[atomType], radicalElectrons=[1], charge=[0], label='*1', lonePairs=[0])
            atom = atom0.copy()
            try:
                atom.applyAction(action)
                self.assertEqual(len(atom.atomType), len(atomType.incrementBond))
                for a in atomType.incrementBond:
                    self.assertTrue(a in atom.atomType)
                self.assertEqual(atom0.radicalElectrons, atom.radicalElectrons)
                self.assertEqual(atom0.charge, atom.charge)
                self.assertEqual(atom0.label, atom.label)
                self.assertEqual(atom0.lonePairs, atom.lonePairs)
            except ActionError:
                self.assertEqual(len(atomType.incrementBond), 0)
    
    def testApplyActionDecrementBond(self):
        """
        Test the GroupAtom.applyAction() method for a CHANGE_BOND action.
        """
        action = ['CHANGE_BOND', '*1', -1, '*2']
        for label, atomType in atomTypes.iteritems():
            atom0 = GroupAtom(atomType=[atomType], radicalElectrons=[1], charge=[0], label='*1', lonePairs=[0])
            atom = atom0.copy()
            try:
                atom.applyAction(action)
                self.assertEqual(len(atom.atomType), len(atomType.decrementBond))
                for a in atomType.decrementBond:
                    self.assertTrue(a in atom.atomType)
                self.assertEqual(atom0.radicalElectrons, atom.radicalElectrons)
                self.assertEqual(atom0.charge, atom.charge)
                self.assertEqual(atom0.label, atom.label)
                self.assertEqual(atom0.lonePairs, atom.lonePairs)
            except ActionError:
                self.assertEqual(len(atomType.decrementBond), 0)
    
    def testApplyActionGainRadical(self):
        """
        Test the GroupAtom.applyAction() method for a GAIN_RADICAL action.
        """
        action = ['GAIN_RADICAL', '*1', 1]
        for label, atomType in atomTypes.iteritems():
            atom0 = GroupAtom(atomType=[atomType], radicalElectrons=[1], charge=[0], label='*1', lonePairs=[0])
            atom = atom0.copy()
            try:
                atom.applyAction(action)
                self.assertEqual(len(atom.atomType), len(atomType.incrementRadical))
                for a in atomType.incrementRadical:
                    self.assertTrue(a in atom.atomType, "GAIN_RADICAL on {0} gave {1} not {2}".format(atomType, atom.atomType, atomType.incrementRadical))
                self.assertEqual(atom0.radicalElectrons, [r - 1 for r in atom.radicalElectrons])
                self.assertEqual(atom0.charge, atom.charge)
                self.assertEqual(atom0.label, atom.label)
                self.assertEqual(atom0.lonePairs, atom.lonePairs)
            except ActionError:
                self.assertEqual(len(atomType.incrementRadical), 0)

        #test when radicals unspecified
        group = Group().fromAdjacencyList("""
        1 R ux
        """) #ux causes a wildcare for radicals
        atom1 = group.atoms[0]
        atom1.applyAction(action)
        self.assertListEqual(atom1.radicalElectrons, [1,2,3,4])
    
    def testApplyActionLoseRadical(self):
        """
        Test the GroupAtom.applyAction() method for a LOSE_RADICAL action.
        """
        action = ['LOSE_RADICAL', '*1', 1]
        for label, atomType in atomTypes.iteritems():
            atom0 = GroupAtom(atomType=[atomType], radicalElectrons=[1], charge=[0], label='*1', lonePairs=[0])
            atom = atom0.copy()
            try:
                atom.applyAction(action)
                self.assertEqual(len(atom.atomType), len(atomType.decrementRadical))
                for a in atomType.incrementRadical:
                    self.assertTrue(a in atom.atomType, "LOSE_RADICAL on {0} gave {1} not {2}".format(atomType, atom.atomType, atomType.decrementRadical))
                self.assertEqual(atom0.radicalElectrons, [r + 1 for r in atom.radicalElectrons])
                self.assertEqual(atom0.charge, atom.charge)
                self.assertEqual(atom0.label, atom.label)
                self.assertEqual(atom0.lonePairs, atom.lonePairs)
            except ActionError:
                self.assertEqual(len(atomType.decrementRadical), 0)

        #test when radicals unspecified
        group = Group().fromAdjacencyList("""
        1 R ux
        """) #ux causes a wildcare for radicals
        atom1 = group.atoms[0]
        atom1.applyAction(action)
        self.assertListEqual(atom1.radicalElectrons, [0,1,2,3])

    def testApplyActionGainPair(self):
        """
        Test the GroupAtom.applyAction() method for a GAIN_PAIR action when lonePairs is either specified or not.
        """

        action = ['GAIN_PAIR', '*1', 1]

        #lonePairs specified:
        for label, atomType in atomTypes.iteritems():
            atom0 = GroupAtom(atomType=[atomType], radicalElectrons=[1], charge=[0], label='*1', lonePairs=[0])
            atom = atom0.copy()
            try:
                atom.applyAction(action)
                self.assertEqual(len(atom.atomType), len(atomType.incrementLonePair))
                for a in atomType.incrementLonePair:
                    self.assertTrue(a in atom.atomType,
                                    "GAIN_PAIR on {0} gave {1} not {2}".format(atomType, atom.atomType,
                                                                                  atomType.incrementLonePair))
                self.assertEqual(atom0.radicalElectrons, atom.radicalElectrons)
                self.assertEqual(atom0.charge, atom.charge)
                self.assertEqual(atom0.label, atom.label)
                self.assertEqual(atom0.lonePairs, [r - 1 for r in atom.lonePairs])
            except ActionError:
                self.assertEqual(len(atomType.incrementLonePair), 0)

        #lonePairs unspecified:
        for label, atomType in atomTypes.iteritems():
            atom0 = GroupAtom(atomType=[atomType], radicalElectrons=[1], charge=[0], label='*1')
            atom = atom0.copy()
            try:
                atom.applyAction(action)
                self.assertEqual(len(atom.atomType), len(atomType.incrementLonePair))
                for a in atomType.incrementLonePair:
                    self.assertTrue(a in atom.atomType,
                                    "GAIN_PAIR on {0} gave {1} not {2}".format(atomType, atom.atomType,
                                                                               atomType.incrementLonePair))
                self.assertEqual(atom0.radicalElectrons, atom.radicalElectrons)
                self.assertEqual(atom0.charge, atom.charge)
                self.assertEqual(atom0.label, atom.label)
                self.assertEqual([0,1,2,3], [r - 1 for r in atom.lonePairs])
            except ActionError:
                self.assertEqual(len(atomType.incrementLonePair), 0)

    def testApplyActionLosePair(self):
        """
        Test the GroupAtom.applyAction() method for a LOSE_PAIR action when lonePairs is either specified or not.
        """

        action = ['LOSE_PAIR', '*1', 1]

        # lonePairs specified:
        for label, atomType in atomTypes.iteritems():
            atom0 = GroupAtom(atomType=[atomType], radicalElectrons=[1], charge=[0], label='*1', lonePairs=[1])
            atom = atom0.copy()
            try:
                atom.applyAction(action)
                self.assertEqual(len(atom.atomType), len(atomType.decrementLonePair))
                for a in atomType.decrementLonePair:
                    self.assertTrue(a in atom.atomType,
                                    "LOSE_PAIR on {0} gave {1} not {2}".format(atomType, atom.atomType,
                                                                                  atomType.decrementLonePair))
                self.assertEqual(atom0.radicalElectrons, atom.radicalElectrons)
                self.assertEqual(atom0.charge, atom.charge)
                self.assertEqual(atom0.label, atom.label)
                self.assertEqual(atom0.lonePairs, [r + 1 for r in atom.lonePairs])
            except ActionError:
                self.assertEqual(len(atomType.decrementLonePair), 0)

        #lonePairs unspecified:
        for label, atomType in atomTypes.iteritems():
            atom0 = GroupAtom(atomType=[atomType], radicalElectrons=[1], charge=[0], label='*1')
            atom = atom0.copy()
            try:
                atom.applyAction(action)
                self.assertEqual(len(atom.atomType), len(atomType.decrementLonePair))
                for a in atomType.decrementLonePair:
                    self.assertTrue(a in atom.atomType,
                                    "LOSE_PAIR on {0} gave {1} not {2}".format(atomType, atom.atomType,
                                                                               atomType.decrementLonePair))
                self.assertEqual(atom0.radicalElectrons, atom.radicalElectrons)
                self.assertEqual(atom0.charge, atom.charge)
                self.assertEqual(atom0.label, atom.label)
                self.assertEqual([1,2,3,4], [r + 1 for r in atom.lonePairs])
            except ActionError:
                self.assertEqual(len(atomType.decrementLonePair), 0)

    def testEquivalent(self):
        """
        Test the GroupAtom.equivalent() method.
        """
        for label1, atomType1 in atomTypes.iteritems():
            for label2, atomType2 in atomTypes.iteritems():
                atom1 = GroupAtom(atomType=[atomType1], radicalElectrons=[1], charge=[0], label='*1', lonePairs=[0])
                atom2 = GroupAtom(atomType=[atomType2], radicalElectrons=[1], charge=[0], label='*1', lonePairs=[0])
                if label1 == label2 or atomType2 in atomType1.generic or atomType1 in atomType2.generic:
                    self.assertTrue(atom1.equivalent(atom2), '{0!s} is not equivalent to {1!s}'.format(atom1, atom2))
                    self.assertTrue(atom2.equivalent(atom1), '{0!s} is not equivalent to {1!s}'.format(atom2, atom1))
                else:
                    self.assertFalse(atom1.equivalent(atom2), '{0!s} is equivalent to {1!s}'.format(atom1, atom2))
                    self.assertFalse(atom2.equivalent(atom1), '{0!s} is equivalent to {1!s}'.format(atom2, atom1))
            # Now see if charge and radical count are checked properly
            for charge in range(3):
                for radicals in range(2):
                    for lonePair in range(2):
                        atom3 = GroupAtom(atomType=[atomType1], radicalElectrons=[radicals], charge=[charge], label='*1', lonePairs=[lonePair])
                        if radicals == 1 and charge == 0 and lonePair == 0:
                            self.assertTrue(atom1.equivalent(atom3), '{0!s} is not equivalent to {1!s}'.format(atom1, atom3))
                            self.assertTrue(atom1.equivalent(atom3), '{0!s} is not equivalent to {1!s}'.format(atom3, atom1))
                        else:
                            self.assertFalse(atom1.equivalent(atom3), '{0!s} is equivalent to {1!s}'.format(atom1, atom3))
                            self.assertFalse(atom1.equivalent(atom3), '{0!s} is equivalent to {1!s}'.format(atom3, atom1))

    def testIsSpecificCaseOf(self):
        """
        Test the GroupAtom.isSpecificCaseOf() method.
        """
        for label1, atomType1 in atomTypes.iteritems():
            for label2, atomType2 in atomTypes.iteritems():
                atom1 = GroupAtom(atomType=[atomType1], radicalElectrons=[1], charge=[0], label='*1', lonePairs=[0])
                atom2 = GroupAtom(atomType=[atomType2], radicalElectrons=[1], charge=[0], label='*1', lonePairs=[0])
                # And make more generic types of these two atoms
                atom1gen = GroupAtom(atomType=[atomType1], radicalElectrons=[0, 1], charge=[0, 1], label='*1', lonePairs=[0, 1])
                atom2gen = GroupAtom(atomType=[atomType2], radicalElectrons=[0, 1], charge=[0, 1], label='*1', lonePairs=[0, 1])
                if label1 == label2 or atomType2 in atomType1.generic:
                    self.assertTrue(atom1.isSpecificCaseOf(atom2), '{0!s} is not a specific case of {1!s}'.format(atom1, atom2))
                    self.assertTrue(atom1.isSpecificCaseOf(atom2gen), '{0!s} is not a specific case of {1!s}'.format(atom1, atom2gen))
                    self.assertFalse(atom1gen.isSpecificCaseOf(atom2), '{0!s} is a specific case of {1!s}'.format(atom1gen, atom2))
                else:
                    self.assertFalse(atom1.isSpecificCaseOf(atom2), '{0!s} is a specific case of {1!s}'.format(atom1, atom2))
                    self.assertFalse(atom1.isSpecificCaseOf(atom2gen), '{0!s} is a specific case of {1!s}'.format(atom1, atom2gen))
                    self.assertFalse(atom1gen.isSpecificCaseOf(atom2), '{0!s} is a specific case of {1!s}'.format(atom1gen, atom2))
    
    def testCopy(self):
        """
        Test the GroupAtom.copy() method.
        """
        atom = self.atom.copy()
        self.assertEqual(len(self.atom.atomType), len(atom.atomType))
        self.assertEqual(self.atom.atomType[0].label, atom.atomType[0].label)
        self.assertEqual(self.atom.radicalElectrons, atom.radicalElectrons)
        self.assertEqual(self.atom.charge, atom.charge)
        self.assertEqual(self.atom.label, atom.label)
        self.assertEqual(self.atom.lonePairs, atom.lonePairs)
    
    def testPickle(self):
        """
        Test that a GroupAtom object can be successfully pickled and
        unpickled with no loss of information.
        """
        import cPickle
        atom = cPickle.loads(cPickle.dumps(self.atom))
        self.assertEqual(len(self.atom.atomType), len(atom.atomType))
        self.assertEqual(self.atom.atomType[0].label, atom.atomType[0].label)
        self.assertEqual(self.atom.radicalElectrons, atom.radicalElectrons)
        self.assertEqual(self.atom.charge, atom.charge)
        self.assertEqual(self.atom.label, atom.label)
        self.assertEqual(self.atom.lonePairs, atom.lonePairs)

    def testCountBonds(self):
        """
        Tests the countBonds function
        """
        adjlist = """
1 *2 C u0     {2,[D,T]} {3,S}
2 *3 C u0     {1,[D,T]} {4,B}
3    C ux     {1,S} {5,D}
4    C u[0,1] {2,B}
5    O u0     {3,D}
"""
        test = Group().fromAdjacencyList(adjlist)
        #returns a list of [single, allDouble, rDouble, oDouble, sDouble, triple, benzene]
        self.assertListEqual([1,0,0,0,0,0,0], test.atoms[0].countBonds())
        self.assertListEqual([1,1,1,0,0,1,0], test.atoms[0].countBonds(wildcards = True))
        self.assertListEqual([0,0,0,0,0,0,1], test.atoms[3].countBonds())
        self.assertListEqual([1,1,0,1,0,0,0], test.atoms[2].countBonds())

    def testHasWildcards(self):
        """
        Tests the GroupAtom.hasWildcards() method
        """
        self.assertFalse(self.atom.hasWildcards())
        adjlist = """
1 *2 C     u0     {2,[D,T]} {3,S}
2 *3 C     u0     {1,[D,T]} {4,S}
3    C     ux     {1,S} {5,S}
4    C     u[0,1] {2,S}
5    [C,O] u0     {3,S}
"""
        group = Group().fromAdjacencyList(adjlist)
        for index, atom in enumerate(group.atoms):
            self.assertTrue(atom.hasWildcards(), 'GroupAtom with index {0} should have wildcards, but does not'.format(index))

    def testMakeSampleAtom(self):
        """
        Tests the GroupAtom.makeSampleAtom() method
        """
        newAtom = self.atom.makeSampleAtom()

        self.assertEquals(newAtom.element, elements.__dict__['C'])
        self.assertEquals(newAtom.radicalElectrons, 1)
        self.assertEquals(newAtom.charge, 0)
        self.assertEquals(newAtom.lonePairs, 0)
################################################################################

class TestGroupBond(unittest.TestCase):
    """
    Contains unit tests of the GroupBond class.
    """

    def setUp(self):
        """
        A method called before each unit test in this class.
        """
        self.bond = GroupBond(None, None, order=[2])
        self.orderList = [[1], [2], [3], [1.5], [1,2], [2,1], [2,3], [1,2,3]]
    
    def testGetOrderStr(self):
        """
        test the Bond.getOrderStr() method
        """
        bond = GroupBond(None,None,order = [1,2,3,1.5])
        self.assertEqual(bond.getOrderStr(),['S','D','T','B'])
        
    def testSetOrderStr(self):
        """
        test the Bond.setOrderStr() method
        """
        
        self.bond.setOrderStr(["B",'T'])
        self.assertEqual(set(self.bond.order), set([3,1.5]))
    
    def testGetOrderNum(self):
        """
        test the Bond.getOrderNum() method
        """
        self.assertEqual(self.bond.getOrderNum(),[2])
        
    def testSetOrderNum(self):
        """
        test the Bond.setOrderNum() method
        """
        
        self.bond.setOrderNum([3,1,2])
        self.assertEqual(self.bond.getOrderStr(),['T','S','D'])

    def testIsSingle(self):
        """
        test the Bond.isSingle() method
        """
        self.bond.setOrderNum([1])
        self.assertTrue(self.bond.isSingle())

        #test interaction with wildcards
        self.bond.setOrderNum([1, 2, 3 , 1.5])
        self.assertFalse(self.bond.isSingle(wildcards = False))
        self.assertTrue(self.bond.isSingle(wildcards = True))

    def testIsDouble(self):
        """
        test the Bond.isDouble() method
        """
        self.bond.setOrderNum([2])
        self.assertTrue(self.bond.isDouble())

        #test interaction with wildcards
        self.bond.setOrderNum([1, 2, 3 , 1.5])
        self.assertFalse(self.bond.isDouble(wildcards = False))
        self.assertTrue(self.bond.isDouble(wildcards = True))
    
    def testIsTriple(self):
        """
        test the Bond.isTriple() method
        """
        self.bond.setOrderNum([3])
        self.assertTrue(self.bond.isTriple())

        #test interaction with wildcards
        self.bond.setOrderNum([1, 2, 3 , 1.5])
        self.assertFalse(self.bond.isTriple(wildcards = False))
        self.assertTrue(self.bond.isTriple(wildcards = True))

    def testIsBenzene(self):
        """
        test the Bond.isBenzene() method
        """
        self.bond.setOrderNum([1.5])
        self.assertTrue(self.bond.isBenzene())

        #test interaction with wildcards
        self.bond.setOrderNum([1, 2, 3 , 1.5])
        self.assertFalse(self.bond.isBenzene(wildcards = False))
        self.assertTrue(self.bond.isBenzene(wildcards = True))

    def testApplyActionBreakBond(self):
        """
        Test the GroupBond.applyAction() method for a BREAK_BOND action.
        """
        action = ['BREAK_BOND', '*1', 1, '*2']
        for order0 in self.orderList:
            bond0 = GroupBond(None, None, order=order0)
            bond = bond0.copy()
            try:
                bond.applyAction(action)
                self.fail('GroupBond.applyAction() unexpectedly processed a BREAK_BOND action.')
            except ActionError:
                pass
    
    def testApplyActionFormBond(self):
        """
        Test the GroupBond.applyAction() method for a FORM_BOND action.
        """
        action = ['FORM_BOND', '*1', 1, '*2']
        for order0 in self.orderList:
            bond0 = GroupBond(None, None, order=order0)
            bond = bond0.copy()
            try:
                bond.applyAction(action)
                self.fail('GroupBond.applyAction() unexpectedly processed a FORM_BOND action.')
            except ActionError:
                pass
    
    def testApplyActionIncrementBond(self):
        """
        Test the GroupBond.applyAction() method for a CHANGE_BOND action.
        """
        action = ['CHANGE_BOND', '*1', 1, '*2']
        for order0 in self.orderList:
            bond0 = GroupBond(None, None, order=order0)
            bond = bond0.copy()
            try:
                bond.applyAction(action)
            except ActionError:
                self.assertTrue(3 in order0 or 1.5 in order0)
                
    def testApplyActionDecrementBond(self):
        """
        Test the GroupBond.applyAction() method for a CHANGE_BOND action.
        """
        action = ['CHANGE_BOND', '*1', -1, '*2']
        for order0 in self.orderList:
            bond0 = GroupBond(None, None, order=order0)
            bond = bond0.copy()
            try:
                bond.applyAction(action)
            except ActionError:
                self.assertTrue(1 in order0 or 1.5 in order0)
            
    def testApplyActionGainRadical(self):
        """
        Test the GroupBond.applyAction() method for a GAIN_RADICAL action.
        """
        action = ['GAIN_RADICAL', '*1', 1]
        for order0 in self.orderList:
            bond0 = GroupBond(None, None, order=order0)
            bond = bond0.copy()
            try:
                bond.applyAction(action)
                self.fail('GroupBond.applyAction() unexpectedly processed a GAIN_RADICAL action.')
            except ActionError:
                pass
    
    def testApplyActionLoseRadical(self):
        """
        Test the GroupBond.applyAction() method for a LOSE_RADICAL action.
        """
        action = ['LOSE_RADICAL', '*1', 1]
        for order0 in self.orderList:
            bond0 = GroupBond(None, None, order=order0)
            bond = bond0.copy()
            try:
                bond.applyAction(action)
                self.fail('GroupBond.applyAction() unexpectedly processed a LOSE_RADICAL action.')
            except ActionError:
                pass
    
    def testEquivalent(self):
        """
        Test the GroupBond.equivalent() method.
        """
        for order1 in self.orderList:
            for order2 in self.orderList:
                bond1 = GroupBond(None, None, order=order1)
                bond2 = GroupBond(None, None, order=order2)
                if order1 == order2 or (all([o in order2 for o in order1]) and all([o in order1 for o in order2])):
                    self.assertTrue(bond1.equivalent(bond2))
                    self.assertTrue(bond2.equivalent(bond1))
                else:
                    self.assertFalse(bond1.equivalent(bond2))
                    self.assertFalse(bond2.equivalent(bond1))
    
    def testIsSpecificCaseOf(self):
        """
        Test the GroupBond.isSpecificCaseOf() method.
        """
        for order1 in self.orderList:
            for order2 in self.orderList:
                bond1 = GroupBond(None, None, order=order1)
                bond2 = GroupBond(None, None, order=order2)
                if order1 == order2 or all([o in order2 for o in order1]):
                    self.assertTrue(bond1.isSpecificCaseOf(bond2))
                else:
                    self.assertFalse(bond1.isSpecificCaseOf(bond2))
                
    def testCopy(self):
        """
        Test the GroupBond.copy() method.
        """
        bond = self.bond.copy()
        self.assertEqual(len(self.bond.order), len(bond.order))
        self.assertEqual(self.bond.order, bond.order)
    
    def testPickle(self):
        """
        Test that a GroupBond object can be successfully pickled and
        unpickled with no loss of information.
        """
        import cPickle
        bond = cPickle.loads(cPickle.dumps(self.bond))
        self.assertEqual(len(self.bond.order), len(bond.order))
        self.assertEqual(self.bond.order, bond.order)

################################################################################

class TestGroup(unittest.TestCase):
    """
    Contains unit tests of the Graph class.
    """

    def setUp(self):
        self.adjlist = """
1 *2 [Cs,Cd]   u0 {2,[S,D]} {3,S}
2 *1 [O2s,O2d] u0 {1,[S,D]}
3    R!H       u0 {1,S}
            """
        self.group = Group().fromAdjacencyList(self.adjlist)
        
    def testClearLabeledAtoms(self):
        """
        Test the Group.clearLabeledAtoms() method.
        """
        self.group.clearLabeledAtoms()
        for atom in self.group.atoms:
            self.assertEqual(atom.label, '')

    def testContainsLabeledAtom(self):
        """
        Test the Group.containsLabeledAtom() method.
        """
        for atom in self.group.atoms:
            if atom.label != '':
                self.assertTrue(self.group.containsLabeledAtom(atom.label))
        self.assertFalse(self.group.containsLabeledAtom('*3'))
        self.assertFalse(self.group.containsLabeledAtom('*4'))
        self.assertFalse(self.group.containsLabeledAtom('*5'))
        self.assertFalse(self.group.containsLabeledAtom('*6'))
        
    def testGetLabeledAtom(self):
        """
        Test the Group.getLabeledAtom() method.
        """
        for atom in self.group.atoms:
            if atom.label != '':
                self.assertEqual(atom, self.group.getLabeledAtom(atom.label))
        try:
            self.group.getLabeledAtom('*3')
            self.fail('Unexpected successful return from Group.getLabeledAtom() with invalid atom label.')
        except ValueError:
            pass
            
    def testGetLabeledAtoms(self):
        """
        Test the Group.getLabeledAtoms() method.
        """
        labeled = self.group.getLabeledAtoms()
        for atom in self.group.atoms:
            if atom.label != '':
                self.assertTrue(atom.label in labeled)
                self.assertTrue(atom in labeled.values())
            else:
                self.assertFalse(atom.label in labeled)
                self.assertFalse(atom in labeled.values())

    def testFromAdjacencyList(self):
        """
        Test the Group.fromAdjacencyList() method.
        """
        atom1, atom2, atom3 = self.group.atoms
        self.assertTrue(self.group.hasBond(atom1,atom2))
        self.assertTrue(self.group.hasBond(atom1,atom3))
        self.assertFalse(self.group.hasBond(atom2,atom3))
        bond12 = atom1.bonds[atom2]
        bond13 = atom1.bonds[atom3]
           
        self.assertTrue(atom1.label == '*2')
        self.assertTrue(atom1.atomType[0].label in ['Cs','Cd'])
        self.assertTrue(atom1.atomType[1].label in ['Cs','Cd'])
        self.assertTrue(atom1.radicalElectrons == [0])
        
        self.assertTrue(atom2.label == '*1')
        self.assertTrue(atom2.atomType[0].label in ['O2s','O2d'])
        self.assertTrue(atom2.atomType[1].label in ['O2s','O2d'])
        self.assertTrue(atom2.radicalElectrons == [0])
        
        self.assertTrue(atom3.label == '')
        self.assertTrue(atom3.atomType[0].label == 'R!H')
        self.assertTrue(atom3.radicalElectrons == [0])

        self.assertTrue(bond12.order == [1,2])
        self.assertTrue(bond13.isSingle())

    def testToAdjacencyList(self):
        """
        Test the Group.toAdjacencyList() method.
        """
        adjlist = self.group.toAdjacencyList()
        self.assertEqual(adjlist.strip(), self.adjlist.strip(),adjlist)

    def testIsIsomorphic(self):
        """
        Test the Group.isIsomorphic() method.
        """
        adjlist = """
1  *1 [O2s,O2d] u0 {3,[S,D]}
2     R!H       u0 {3,S}
3  *2 [Cs,Cd]   u0 {1,[S,D]} {2,S}
            """
        group = Group().fromAdjacencyList(adjlist)
        self.assertTrue(self.group.isIsomorphic(group))
        self.assertTrue(group.isIsomorphic(self.group))
        
    def testFindIsomorphism(self):
        """
        Test the Group.findIsomorphism() method.
        """
        adjlist = """
1  *1 [O2s,O2d] u0 {3,[S,D]}
2     R!H       u0 {3,S}
3  *2 [Cs,Cd]   u0 {1,[S,D]} {2,S}
            """
        group = Group().fromAdjacencyList(adjlist)
        result = self.group.findIsomorphism(group)
        self.assertEqual(len(result), 1)
        for atom1, atom2 in result[0].items():
            self.assertTrue(atom1 in self.group.atoms)
            self.assertTrue(atom2 in group.atoms)
            self.assertTrue(atom1.equivalent(atom2))
            for atom3 in atom1.bonds:
                atom4 = result[0][atom3]
                self.assertTrue(atom4 in atom2.bonds)
                self.assertTrue(atom3.equivalent(atom4))
                bond1 = atom1.bonds[atom3]
                bond2 = atom2.bonds[atom4]
                self.assertTrue(bond1.equivalent(bond2))
        
    def testIsSubgraphIsomorphic(self):
        """
        Test the Group.isSubgraphIsomorphic() method.
        """
        adjlist = """
1  *1 [Cs,Cd] u0
            """
        group = Group().fromAdjacencyList(adjlist)
        self.assertTrue(self.group.isSubgraphIsomorphic(group))
        self.assertFalse(group.isIsomorphic(self.group))
        
    def testFindSubgraphIsomorphisms(self):
        """
        Test the Group.findSubgraphIsomorphisms() method.
        """
        adjlist = """
1  *1 [Cs,Cd] u0
            """
        group = Group().fromAdjacencyList(adjlist)
        result = self.group.findSubgraphIsomorphisms(group)
        self.assertEqual(len(result), 1)
        for atom1, atom2 in result[0].iteritems():
            self.assertTrue(atom1 in self.group.atoms)
            self.assertTrue(atom2 in group.atoms)
            self.assertTrue(atom1.equivalent(atom2))
    
    def testGenerateExtensions(self):
        """
        test that appropriate group extensions are being generated
        """
        
        testGrp = Group().fromAdjacencyList("""
1 *2 C u0 {2,[S,D]} 
2 *1 C u[0,1] {1,[S,D]} {3,S}
3    R!H     u0 {2,S}
            """)
        
        extensions = testGrp.getExtensions(R=[atomTypes[i] for i in ['C','O','H']])
        extensions = [a[0] for a in extensions]
        
        ans = ['1 *2 C   u0     {2,[S,D]} {4,[S,D,T,B]}\n2 *1 C   u[0,1] {1,[S,D]} {3,S}\n3    R!H u0     {2,S}\n4    R!H ux     {1,[S,D,T,B]}\n',
 '1 *2 C   u0     {2,S}\n2 *1 C   u[0,1] {1,S} {3,S}\n3    R!H u0     {2,S}\n',
 '1 *2 C   u0     {2,D}\n2 *1 C   u[0,1] {1,D} {3,S}\n3    R!H u0     {2,S}\n',
 '1 *2 C   u0 {2,[S,D]}\n2 *1 C   u0 {1,[S,D]} {3,S}\n3    R!H u0 {2,S}\n',
 '1 *2 C   u0 {2,[S,D]}\n2 *1 C   u1 {1,[S,D]} {3,S}\n3    R!H u0 {2,S}\n',
 '1 *2 C   u0     {2,[S,D]}\n2 *1 C   u[0,1] {1,[S,D]} {3,S} {4,[S,D,T,B]}\n3    R!H u0     {2,S}\n4    R!H ux     {2,[S,D,T,B]}\n',
 '1 *2 C   u0     {2,S}\n2 *1 C   u[0,1] {1,S} {3,S}\n3    R!H u0     {2,S}\n',
 '1 *2 C   u0     {2,D}\n2 *1 C   u[0,1] {1,D} {3,S}\n3    R!H u0     {2,S}\n',
 '1 *2 C u0     {2,[S,D]}\n2 *1 C u[0,1] {1,[S,D]} {3,S}\n3    C u0     {2,S}\n',
 '1 *2 C u0     {2,[S,D]}\n2 *1 C u[0,1] {1,[S,D]} {3,S}\n3    O u0     {2,S}\n',
 '1 *2 C   u0     {2,[S,D]}\n2 *1 C   u[0,1] {1,[S,D]} {3,S}\n3    R!H u0     {2,S} {4,[S,D,T,B]}\n4    R!H ux     {3,[S,D,T,B]}\n',
 '1 *2 C   u0     {2,[S,D]} {3,[S,D,T,B]}\n2 *1 C   u[0,1] {1,[S,D]} {3,S}\n3    R!H u0     {1,[S,D,T,B]} {2,S}\n']
        ans = [Group().fromAdjacencyList(k) for k in ans]

        for v in ans:
            boos = [ext.isIdentical(v) and ext.isSubgraphIsomorphic(v,generateInitialMap=True) for ext in extensions]
            self.assertTrue(any(boos),'generated extensions did not match expected extensions')
    
    def testGeneratedExtensionsSubgraphs(self):
        testGrp = Group().fromAdjacencyList("""
1 *2 C u0 {2,[S,D]} 
2 *1 C u[0,1] {1,[S,D]} {3,S}
3    R!H     u0 {2,S}
            """)
        extensions = testGrp.getExtensions(R=[atomTypes[i] for i in ['C','O','H']])
        extensions = [a[0] for a in extensions]
        
        for ext in extensions:
            self.assertTrue(ext.isSubgraphIsomorphic(testGrp,generateInitialMap=True))
        
    def testPickle(self):
        """
        Test that a Group object can be successfully pickled and
        unpickled with no loss of information.
        """
        import cPickle
        group = cPickle.loads(cPickle.dumps(self.group))
        
        self.assertEqual(len(self.group.atoms), len(group.atoms))
        for atom0, atom in zip(group.atoms, self.group.atoms):
            self.assertTrue(atom0.equivalent(atom))

        self.assertTrue(self.group.isIsomorphic(group))
        self.assertTrue(group.isIsomorphic(self.group))

    def testCreateAndConnectAtom(self):
        """
        Tests createAndConnectAtom method
        """
        adjlist1 = """
1  *1 C u0 {2,S}
2  *2 C u0 {1,S}
"""

        answer1 = """
1  *1 C  u0 {2,S} {3,B}
2  *2 C  u0 {1,S}
3     Cb u0 {1,B}
"""

        group1 = Group().fromAdjacencyList(adjlist1)
        answer1 = Group().fromAdjacencyList(answer1)
        atom1 = group1.getLabeledAtom("*1")
        newAtom = group1.createAndConnectAtom(atomtypes = ["Cb"], connectingAtom = atom1, bondOrders = ["B"])
        self.assertTrue(group1.isIsomorphic(answer1))

        answer2 = """
1  *1 C       u0 {2,S} {3,[S,D]}
2  *2 C       u0 {1,S}
3     [Cs,Cd] u0 {1,[S,D]}
"""

        #Test that wildcards work alright
        group2 = Group().fromAdjacencyList(adjlist1)
        answer2 = Group().fromAdjacencyList(answer2)
        atom1 = group2.getLabeledAtom("*1")
        newAtom = group2.createAndConnectAtom(atomtypes = ["Cs", "Cd"], connectingAtom = atom1, bondOrders = ["S","D"])
        self.assertTrue(group2.isIsomorphic(answer2))



    def testAddImplicitAtomsFromAtomType(self):
        """
        test Group.addImplicitAtomsFromAtomType() method
        """
        #basic test adding oDouble
        adjlist1 = """
1  *1 CO u0
            """

        adjlist2 = """
1  *1 CO u0 {2,D}
2     O  u0 {1,D}
            """

        group1 = Group().fromAdjacencyList(adjlist1)
        group2 = Group().fromAdjacencyList(adjlist2)

        newGroup = group1.addImplicitAtomsFromAtomType()
        self.assertTrue(group2.isIsomorphic(newGroup))
        #testing the allDouble match (more complicated
        adjlist3 = """
1  *1 Cdd u0
            """

        adjlist4 = """
1  *1 Cdd u0 {2,D} {3,D}
2     C   u0 {1,D}
3     C   u0 {1,D}
            """
        group3 = Group().fromAdjacencyList(adjlist3)
        group4 = Group().fromAdjacencyList(adjlist4)

        newGroup =group3.addImplicitAtomsFromAtomType()
        self.assertTrue(group4.isIsomorphic(newGroup))
        #test adding a triple bond
        adjlist5 = """
1  *1 Ct u0
            """

        adjlist6 = """
1  *1 Ct u0 {2,T}
2     C  u0 {1,T}
            """
        group5 = Group().fromAdjacencyList(adjlist5)
        group6 = Group().fromAdjacencyList(adjlist6)

        newGroup =group5.addImplicitAtomsFromAtomType()
        self.assertTrue(group6.isIsomorphic(newGroup))
        #test addition of lone pairs
        adjlist7 = """
1  *1 N1dc u0
            """

        adjlist8 = """
1  *1 N1dc u0 p2 {2,D}
2     C    u0 {1,D}
            """
        group7 = Group().fromAdjacencyList(adjlist7)
        group8 = Group().fromAdjacencyList(adjlist8)

        newGroup = group7.addImplicitAtomsFromAtomType()
        self.assertTrue(group8.isIsomorphic(newGroup))

        #test multiple implicit atoms at a time
        adjlist9 = """
1  *1 Cd u0 {2,S}
2     Ct u0 {1,S}
            """

        adjlist10 = """
1  *1 C  u0 {2,S} {3,D}
2     Ct u0 {1,S} {4,T}
3     C  u0 {1,D}
4     C  u0 {2,T}
            """
        group9 = Group().fromAdjacencyList(adjlist9)
        group10 = Group().fromAdjacencyList(adjlist10)

        newGroup =group9.addImplicitAtomsFromAtomType()
        self.assertTrue(group10.isIsomorphic(newGroup))

    def testClassifyBenzeneCarbons(self):
        """
        Tests the method classifyingBenzeneCarbons
        """

        #This tests that we classify Cb atom types correctly
        adjlist1 = """
1  *1 Cb u0 {2,B}
2  *2 Cb u0 {1,B}
"""
        group1 = Group().fromAdjacencyList(adjlist1)

        (cbAtomList, cbfAtomList, cbfAtomList1, cbfAtomList2, connectedCbfs)=group1.classifyBenzeneCarbons()
        self.assertEquals(len(cbfAtomList), 0)
        for atom in group1.atoms:
            self.assertIn(atom, cbAtomList)

        #This tests that we classify Cbf atomtypes correctly
        adjlist2 = """
1 *1 Cbf u0
"""
        group2 = Group().fromAdjacencyList(adjlist2)
        (cbAtomList, cbfAtomList, cbfAtomList1, cbfAtomList2, connectedCbfs)=group2.classifyBenzeneCarbons()
        self.assertIn(group2.atoms[0], cbfAtomList)
        self.assertIn(group2.atoms[0], cbfAtomList1)

        #This tests that we can classify Cb atoms based on bonding and not just atomtype
        adjlist3 = """
1 *1 C u0 {2,B}
2 *2 C u0 {1,B} {3,B}
3 *3 C u0 {2,B}
"""
        group3 = Group().fromAdjacencyList(adjlist3)
        (cbAtomList, cbfAtomList, cbfAtomList1, cbfAtomList2, connectedCbfs)=group3.classifyBenzeneCarbons()
        for atom in group3.atoms:
            self.assertIn(atom, cbAtomList)

        #This tests that we can classify Cbf1 atoms based on bonding and not just atomtype
        adjlist4 = """
1 *1 C u0 {2,B} {3,B} {4,B}
2 *2 C u0 {1,B}
3 *3 C u0 {1,B}
4 *4 C u0 {1,B}
"""
        group4 = Group().fromAdjacencyList(adjlist4)
        (cbAtomList, cbfAtomList, cbfAtomList1, cbfAtomList2, connectedCbfs)=group4.classifyBenzeneCarbons()
        self.assertEquals(len(cbfAtomList), 1)
        self.assertEquals(len(cbfAtomList1), 1)

        #This tests that we can classify Cbf2 atoms. In the following partial group, we should have:
        #one Cbf2 atom, two Cbf1 atoms, and 5 Cb atoms
        adjlist5 = """
1 *1 C u0 {2,B} {3,B} {4,B}
2 *2 C u0 {1,B} {5,B} {6,B}
3 *3 C u0 {1,B} {7,B} {8,B}
4 *4 C u0 {1,B}
5 *5 C u0 {2,B}
6 *6 C u0 {2,B}
7 *7 C u0 {3,B}
8 *8 C u0 {3,B}
"""
        group5 = Group().fromAdjacencyList(adjlist5)
        (cbAtomList, cbfAtomList, cbfAtomList1, cbfAtomList2, connectedCbfs)=group5.classifyBenzeneCarbons()
        self.assertEquals(len(cbfAtomList1), 2)
        self.assertEquals(len(cbfAtomList2), 1)
        self.assertEquals(len(cbAtomList), 5)

        #Tests that we can classify connected Cbfs correctly. *1 should be connected to both *2 and *3
        atom1 = group5.getLabeledAtom("*1")
        atom2 = group5.getLabeledAtom("*2")
        atom3 = group5.getLabeledAtom("*3")
        self.assertIn(atom2, connectedCbfs[atom1])
        self.assertIn(atom3, connectedCbfs[atom1])
        self.assertIn(atom1, connectedCbfs[atom2])
        self.assertIn(atom1, connectedCbfs[atom3])

    def testSortByConnectivity(self):
        """
        Tests sortByConnectivity method
        """

        #Basic test, we should get *1, *3 *2
        adjlist1 = """
1 *1 C u0 {3,B}
2 *2 C u0 {3,B}
3 *3 C u0 {1,B} {2,B}
"""
        group1 = Group().fromAdjacencyList(adjlist1)
        orderedAtoms = group1.sortByConnectivity(group1.atoms)
        self.assertEquals([x.label for x in orderedAtoms], ["*1", "*3", "*2"])

        #Check a detached case, we should get *1, *3, *4, *2, *5
        adjlist2 = """
1 *1 C u0 {3,B}
2 *2 C u0 {4,S} {5,B}
3 *3 C u0 {1,B} {4,B}
4 *4 C u0 {3,B} {2,S}
5 *5 C u0 {2,B}
"""
        group2 = Group().fromAdjacencyList(adjlist2)
        orderedAtoms = group2.sortByConnectivity(group2.atoms)
        self.assertEquals([x.label for x in orderedAtoms], ["*1", "*3", "*4", "*2", "*5"])


    def testAddImplicitBenzene(self):
        """
        Test the Group.addImplicitBenzene method
        """

        #tests it can make a benzene molecule
        adjlist1 = """
1  *1 Cb u0 {2,B}
2  *2 Cb u0 {1,B}
            """
        #tests it can make a bi-phenyl
        adjlist2 = """
1  *1 Cb u0 {2,S}
2  *2 Cb u0 {1,S}
            """
        #tests it can make a napthalene
        adjlist3 = """
1  *1 Cbf u0
            """

        #Test handling of Cbf2 atoms
        adjlist4 = """
1  *1 Cbf u0 p2 c0 {2,B}
2  *2 Cbf u0 p0 c0 {1,B} {3,B}
3  *3 Cbf u0 p0 c0 {2,B}
    """

        #test handling of heteroatoms and wildcards
        adjlist5 = """
1 *1 Cbf u0 {2,B} {3,B} {4,B}
2    R!H u0 {1,B}
3    R!H u0 {1,B}
4    R!H u0 {1,B}
    """
        adjlist6 = """
1  *1 Cbf u0 p2 c0 {2,B}
2  *2 Cb u0 p0 c0 {1,B} {3,B}
3  *3 Cb u0 p0 c0 {2,B} {4,S}
4  *4 O  u0 p0 c0 {3,S}
    """

        adjlist7="""
1 *1 Cb u0 {4,B}
2 *2 Cb u0 {3,B}
3 *3 Cb u0 {4,B} {2,B}
4 *4 Cb u0 {1,B} {3,B}
"""

        benzene ="""
1 C u0 {2,B} {6,B}
2 C u0 {1,B} {3,B}
3 C u0 {2,B} {4,B}
4 C u0 {3,B} {5,B}
5 C u0 {4,B} {6,B}
6 C u0 {5,B} {1,B}
        """

        biphenyl ="""
1  C u0 {2,B} {6,B} {7,S}
2  C u0 {1,B} {3,B}
3  C u0 {2,B} {4,B}
4  C u0 {3,B} {5,B}
5  C u0 {4,B} {6,B}
6  C u0 {5,B} {1,B}
7  C u0 {8,B} {12,B} {1,S}
8  C u0 {7,B} {9,B}
9  C u0 {8,B} {10,B}
10 C u0 {9,B} {11,B}
11 C u0 {10,B} {12,B}
12 C u0 {11,B} {7,B}
        """

        naphthalene ="""
1  C u0 {2,B} {10,B}
2  C u0 {1,B} {3,B}
3  C u0 {2,B} {4,B}
4  C u0 {3,B} {5,B} {9,B}
5  C u0 {4,B} {6,B}
6  C u0 {5,B} {7,B}
7  C u0 {6,B} {8,B}
8  C u0 {7,B} {9,B}
9  C u0 {4,B} {8,B} {10,B}
10 C u0 {1,B} {9,B}
        """

        phenanthrene = """
1  Cbf u0 p2 c0 {2,B} {7,B} {11,B}
2  Cbf u0 p0 c0 {1,B} {3,B} {5,B}
3  Cbf u0 p0 c0 {2,B} {4,B} {6,B}
4  C   u0 {3,B} {8,B} {14,B}
5  C   u0 {2,B} {9,B}
6  C   u0 {3,B} {12,B}
7  C   u0 {1,B} {8,B}
8  C   u0 {4,B} {7,B}
9  C   u0 {5,B} {10,B}
10 C   u0 {9,B} {11,B}
11 C   u0 {1,B} {10,B}
12 C   u0 {6,B} {13,B}
13 C   u0 {12,B} {14,B}
14 C   u0 {4,B} {13,B}
    """

        answer5 = """
1  *1 Cbf u0 {2,B} {3,B} {4,B}
2     R!H u0 {1,B} {5,B}
3     R!H u0 {1,B} {7,B} {10,B}
4     R!H u0 {1,B} {8,B}
5     Cb  u0 {2,B} {6,B}
6     Cb  u0 {5,B} {7,B}
7     Cb  u0 {3,B} {6,B}
8     Cb  u0 {4,B} {9,B}
9     Cb  u0 {8,B} {10,B}
10    Cb  u0 {3,B} {9,B}
"""
        answer6="""
1  *1 Cbf u0 p2 c0 {2,B} {5,B} {8,B}
2  *2 Cb  u0 p0 c0 {1,B} {3,B} {11,B}
3  *3 Cb  u0 p0 c0 {2,B} {4,S} {7,B}
4  *4 O   u0 p0 c0 {3,S}
5     Cb  u0 {1,B} {6,B}
6     Cb  u0 {5,B} {7,B}
7     Cb  u0 {3,B} {6,B}
8     Cb  u0 {1,B} {9,B}
9     Cb  u0 {8,B} {10,B}
10    Cb  u0 {9,B} {11,B}
11    Cb  u0 {2,B} {10,B}
"""

        group1 = Group().fromAdjacencyList(adjlist1)
        group2 = Group().fromAdjacencyList(adjlist2)
        group3 = Group().fromAdjacencyList(adjlist3)
        group4 = Group().fromAdjacencyList(adjlist4)
        group5 = Group().fromAdjacencyList(adjlist5)
        group6 = Group().fromAdjacencyList(adjlist6)
        group7 = Group().fromAdjacencyList(adjlist7)

        benzeneGroup = Group().fromAdjacencyList(benzene)
        biphenylGroup = Group().fromAdjacencyList(biphenyl)
        naphthaleneGroup = Group().fromAdjacencyList(naphthalene)
        phenanthreneGroup = Group().fromAdjacencyList(phenanthrene)
        answer5 = Group().fromAdjacencyList(answer5)
        answer6 = Group().fromAdjacencyList(answer6)

        group1 = group1.addImplicitBenzene()
        self.assertTrue(benzeneGroup.isIsomorphic(group1))
        group2 = group2.addImplicitBenzene()
        self.assertTrue(biphenylGroup.isIsomorphic(group2))
        group3 = group3.addImplicitBenzene()
        self.assertTrue(naphthaleneGroup.isIsomorphic(group3))
        group4 = group4.addImplicitBenzene()
        self.assertTrue(phenanthreneGroup.isIsomorphic(group4))
        group5 = group5.addImplicitBenzene()
        self.assertTrue(answer5.isIsomorphic(group5))
        group6 = group6.addImplicitBenzene()
        self.assertTrue(answer6.isIsomorphic(group6))
        group7 = group7.addImplicitBenzene()
        self.assertTrue(benzeneGroup.isIsomorphic(group7))

    def testPickWildcards(self):
        """
        Test the Group.pickWildCards function
        """
        #The following tests are for picking optimal bond orders when there are bond wilcards
        #test that Cb/Cbf atoms with [D,B] chooses [B] instead of [D] bonds
        adjlist1 = """
    1 *1 R!H       u1 {2,[D,B]}
    2 *2 [Cbf,Cdd] u0 {1,[D,B]} {3,[D,B]}
    3 *3 [Cb,Cd]   u0 {2,[D,B]} {4,S}
    4 *4 R!H       u0 {3,S} {5,S}
    5 *5 H         u0 {4,S}
"""
        group1 = Group().fromAdjacencyList(adjlist1)
        group1.pickWildcards()
        atoms = group1.atoms
        self.assertTrue(atoms[0].bonds[atoms[1]].isBenzene())
        self.assertTrue(atoms[1].bonds[atoms[2]].isBenzene())

        adjlist2 = """
    1 *1 R!H       u1 {2,[S,D]} {4,[S,D]}
    2 *2 [CO,Cdd]  u0 {1,[S,D]} {3,[S,D]}
    3 *3 [O2d,Cd]  u0 {2,[S,D]}
    4 *4 [Cdd,Cd]  u0 {1,[S,D]}
"""
        group2 = Group().fromAdjacencyList(adjlist2)
        group2.pickWildcards()
        atoms = group2.atoms
        self.assertTrue(atoms[1].bonds[atoms[2]].isDouble)
        self.assertTrue(atoms[0].bonds[atoms[3]].isDouble)

    def testMakeSampleMolecule(self):
        """
        Test the Group.makeSampleMolecule method
        """

        def performSampMoleComparison(adjlist, answer_smiles):
            """
            Creates a sample molecule from the adjlist and returns if it is isomorphic to a molecule created from
            the inputted smiles
            """
            group = Group().fromAdjacencyList(adjlist)
            result = group.makeSampleMolecule()
            return (result.isIsomorphic(Molecule().fromSMILES(answer_smiles)))
########################################################################################################################
        #tests adding implicit atoms
        adjlist = """
1  *1 Cd u0
            """
        answer_smiles = 'C=C'
        self.assertTrue(performSampMoleComparison(adjlist, answer_smiles))

        #test creating implicit benzene atoms
        adjlist2 = """
1  *1 Cbf u0 {2,B}
2     Cbf u0 {1,B}
            """

        group2 = Group().fromAdjacencyList(adjlist2)
        result2 = group2.makeSampleMolecule()
        naphthaleneMolecule = Molecule().fromSMILES('C1=CC=C2C=CC=CC2=C1')
        resonanceList2=naphthaleneMolecule.generate_resonance_structures()
        self.assertTrue(any([result2.isIsomorphic(x) for x in resonanceList2]))

        #test the creation of a positively charged species
        adjlist = """
1  *1 N5sc u0
        """
        answer_smiles = '[NH4+]'
        self.assertTrue(performSampMoleComparison(adjlist, answer_smiles))

        #test the creation of a negatively charged species
        adjlist = """
1  *1 N1sc u0
        """
        answer_smiles = '[NH2-]'
        self.assertTrue(performSampMoleComparison(adjlist, answer_smiles))

        #test creation of charged species when some single bonds present
        adjlist = """
1 *2 [N5sc,N5dc] u0 {2,S} {3,S}
2 *3 R!H         u1 {1,S}
3 *4 H           u0 {1,S}
"""
        answer_smiles = '[NH3+][CH2]'
        self.assertTrue(performSampMoleComparison(adjlist, answer_smiles))

    def testIsBenzeneExplicit(self):
        """
        Test the Group.isBenzeneExplicit method
        """
        adjlist1 = """
1  *1 Cb u0 {2,B}
2  *2 Cb u0 {1,B}
        """
        group1 = Group().fromAdjacencyList(adjlist1)
        self.assertFalse(group1.isBenzeneExplicit())

        benzene ="""
1 C u0 {2,B} {6,B}
2 C u0 {1,B} {3,B}
3 C u0 {2,B} {4,B}
4 C u0 {3,B} {5,B}
5 C u0 {4,B} {6,B}
6 C u0 {5,B} {1,B}
        """
        benzene = Group().fromAdjacencyList(benzene)
        self.assertTrue(benzene.isBenzeneExplicit())

    def test_repr_png(self):
        """Test that a png representation can be created."""
        adjlist = """
1 *1 [C,Cd,Ct,CO,CS,Cb] u1 {2,[S,D,T,B]}
2 *2 [C,Cd,Ct,CO,CS,Cb] u0 {1,[S,D,T,B]} {3,[S,D,T,B]}
3 *3 [C,Cd,Ct,CO,CS,Cb] u0 {2,[S,D,T,B]} {4,[S,D,T,B]}
4 *4 [C,Cd,Ct,CO,CS,Cb] u0 {3,[S,D,T,B]}
        """
        group = Group().fromAdjacencyList(adjlist)
        result = group._repr_png_()
        self.assertIsNotNone(result)

    def testDrawGroup(self):
        """Test that the draw method returns the expected pydot graph."""
        adjlist = """
1 *1 [C,Cd,Ct,CO,CS,Cb] u1 {2,[S,D,T,B]}
2 *2 [C,Cd,Ct,CO,CS,Cb] u0 {1,[S,D,T,B]} {3,[S,D,T,B]}
3 *3 [C,Cd,Ct,CO,CS,Cb] u0 {2,[S,D,T,B]} {4,[S,D,T,B]}
4 *4 [C,Cd,Ct,CO,CS,Cb] u0 {3,[S,D,T,B]}
        """
        # Use of tabs in the expected string is intentional
        expected = """
graph G {
	graph [dpi=52];
	node [label="\N"];
	1	 [fontname=Helvetica,
		fontsize=16,
		label="*1 C,Cd,Ct,CO,CS,Cb"];
	2	 [fontname=Helvetica,
		fontsize=16,
		label="*2 C,Cd,Ct,CO,CS,Cb"];
	1 -- 2	 [fontname=Helvetica,
		fontsize=16,
		label="S,D,T,B"];
	3	 [fontname=Helvetica,
		fontsize=16,
		label="*3 C,Cd,Ct,CO,CS,Cb"];
	2 -- 3	 [fontname=Helvetica,
		fontsize=16,
		label="S,D,T,B"];
	4	 [fontname=Helvetica,
		fontsize=16,
		label="*4 C,Cd,Ct,CO,CS,Cb"];
	3 -- 4	 [fontname=Helvetica,
		fontsize=16,
		label="S,D,T,B"];
}
        """
        group = Group().fromAdjacencyList(adjlist)
        result = group.draw('canon')
        self.assertEqual(''.join(result.split()), ''.join(expected.split()))

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

        mergedGroup = backbone1.mergeGroups(end1)
        self.assertTrue(mergedGroup.isIdentical(desiredMerge1))

        #test it works when there is a cyclical structure to the backbone

        backbone2 = Group().fromAdjacencyList("""
1 *1 R!H u1 {2,S} {4,S}
2 *4 R!H u0 {1,S} {3,S}
3 *2 R!H u0 {2,S} {4,S}
4 *3 R!H u0 {3,S} {1,S}
""")

        end2 = Group().fromAdjacencyList("""
1 *2 O2s u0 {2,S}
2 *3 Cs  u0 {1,S}
""")
        desiredMerge2 = Group().fromAdjacencyList("""
1 *1 R!H u1 {2,S} {4,S}
2 *4 R!H u0 {1,S} {3,S}
3 *2 O2s u0 {2,S} {4,S}
4 *3 Cs  u0 {3,S} {1,S}
""")
        mergedGroup = backbone2.mergeGroups(end2)
        self.assertTrue(mergedGroup.isIdentical(desiredMerge2))

    def test_get_element_count(self):
        """Test that we can count elements properly."""
        group = Group().fromAdjacencyList("""
1 R!H u0 {2,S}
2 [Cs,Cd,Ct,Cb] u0 {1,S} {3,S}
3 [Cs,Cd,Ct,Cb,O2s,S2s] u0 {2,S} {4,S}
4 N1s u0 {3,S}
""")
        expected = {'C': 1, 'N': 1}
        result = group.get_element_count()
        self.assertEqual(expected, result)

################################################################################

if __name__ == '__main__':
    unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )
