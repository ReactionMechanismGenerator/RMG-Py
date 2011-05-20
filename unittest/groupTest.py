#!/usr/bin/python
# -*- coding: utf-8 -*-

import unittest

from rmgpy.group import *
from rmgpy.atomtype import atomTypes

################################################################################

class TestGroupAtom(unittest.TestCase):
    """
    Contains unit tests of the GroupAtom class.
    """

    def setUp(self):
        """
        A method called before each unit test in this class.
        """
        self.atom = GroupAtom(atomType=[atomTypes['Cd']], radicalElectrons=[1], spinMultiplicity=[2], charge=[0], label='*1')
    
    def testApplyActionBreakBond(self):
        """
        Test the GroupAtom.applyAction() method for a BREAK_BOND action.
        """
        action = ['BREAK_BOND', '*1', 'S', '*2']
        for label, atomType in atomTypes.iteritems():
            atom0 = GroupAtom(atomType=[atomType], radicalElectrons=[1], spinMultiplicity=[2], charge=[0], label='*1')
            atom = atom0.copy()
            try:
                atom.applyAction(action)
                self.assertEqual(len(atom.atomType), len(atomType.breakBond))
                for a in atomType.breakBond:
                    self.assertTrue(a in atom.atomType)
                self.assertEqual(atom0.radicalElectrons, atom.radicalElectrons)
                self.assertEqual(atom0.spinMultiplicity, atom.spinMultiplicity)
                self.assertEqual(atom0.charge, atom.charge)
                self.assertEqual(atom0.label, atom.label)
            except ActionError:
                self.assertEqual(len(atomType.breakBond), 0)
    
    def testApplyActionFormBond(self):
        """
        Test the GroupAtom.applyAction() method for a FORM_BOND action.
        """
        action = ['FORM_BOND', '*1', 'S', '*2']
        for label, atomType in atomTypes.iteritems():
            atom0 = GroupAtom(atomType=[atomType], radicalElectrons=[1], spinMultiplicity=[2], charge=[0], label='*1')
            atom = atom0.copy()
            try:
                atom.applyAction(action)
                self.assertEqual(len(atom.atomType), len(atomType.formBond))
                for a in atomType.formBond:
                    self.assertTrue(a in atom.atomType)
                self.assertEqual(atom0.radicalElectrons, atom.radicalElectrons)
                self.assertEqual(atom0.spinMultiplicity, atom.spinMultiplicity)
                self.assertEqual(atom0.charge, atom.charge)
                self.assertEqual(atom0.label, atom.label)
            except ActionError:
                self.assertEqual(len(atomType.formBond), 0)
    
    def testApplyActionIncrementBond(self):
        """
        Test the GroupAtom.applyAction() method for a CHANGE_BOND action.
        """
        action = ['CHANGE_BOND', '*1', 1, '*2']
        for label, atomType in atomTypes.iteritems():
            atom0 = GroupAtom(atomType=[atomType], radicalElectrons=[1], spinMultiplicity=[2], charge=[0], label='*1')
            atom = atom0.copy()
            try:
                atom.applyAction(action)
                self.assertEqual(len(atom.atomType), len(atomType.incrementBond))
                for a in atomType.incrementBond:
                    self.assertTrue(a in atom.atomType)
                self.assertEqual(atom0.radicalElectrons, atom.radicalElectrons)
                self.assertEqual(atom0.spinMultiplicity, atom.spinMultiplicity)
                self.assertEqual(atom0.charge, atom.charge)
                self.assertEqual(atom0.label, atom.label)
            except ActionError:
                self.assertEqual(len(atomType.incrementBond), 0)
    
    def testApplyActionDecrementBond(self):
        """
        Test the GroupAtom.applyAction() method for a CHANGE_BOND action.
        """
        action = ['CHANGE_BOND', '*1', -1, '*2']
        for label, atomType in atomTypes.iteritems():
            atom0 = GroupAtom(atomType=[atomType], radicalElectrons=[1], spinMultiplicity=[2], charge=[0], label='*1')
            atom = atom0.copy()
            try:
                atom.applyAction(action)
                self.assertEqual(len(atom.atomType), len(atomType.decrementBond))
                for a in atomType.decrementBond:
                    self.assertTrue(a in atom.atomType)
                self.assertEqual(atom0.radicalElectrons, atom.radicalElectrons)
                self.assertEqual(atom0.spinMultiplicity, atom.spinMultiplicity)
                self.assertEqual(atom0.charge, atom.charge)
                self.assertEqual(atom0.label, atom.label)
            except ActionError:
                self.assertEqual(len(atomType.decrementBond), 0)
    
    def testApplyActionGainRadical(self):
        """
        Test the GroupAtom.applyAction() method for a GAIN_RADICAL action.
        """
        action = ['GAIN_RADICAL', '*1', 1]
        for label, atomType in atomTypes.iteritems():
            atom0 = GroupAtom(atomType=[atomType], radicalElectrons=[1], spinMultiplicity=[2], charge=[0], label='*1')
            atom = atom0.copy()
            try:
                atom.applyAction(action)
                self.assertEqual(len(atom.atomType), len(atomType.incrementRadical))
                for a in atomType.incrementRadical:
                    self.assertTrue(a in atom.atomType)
                self.assertEqual(atom0.radicalElectrons, [r - 1 for r in atom.radicalElectrons])
                self.assertEqual(atom0.spinMultiplicity, [s - 1 for s in atom.spinMultiplicity])
                self.assertEqual(atom0.charge, atom.charge)
                self.assertEqual(atom0.label, atom.label)
            except ActionError:
                self.assertEqual(len(atomType.incrementRadical), 0)
    
    def testApplyActionLoseRadical(self):
        """
        Test the GroupAtom.applyAction() method for a LOSE_RADICAL action.
        """
        action = ['LOSE_RADICAL', '*1', 1]
        for label, atomType in atomTypes.iteritems():
            atom0 = GroupAtom(atomType=[atomType], radicalElectrons=[1], spinMultiplicity=[2], charge=[0], label='*1')
            atom = atom0.copy()
            try:
                atom.applyAction(action)
                self.assertEqual(len(atom.atomType), len(atomType.incrementRadical))
                for a in atomType.incrementRadical:
                    self.assertTrue(a in atom.atomType)
                self.assertEqual(atom0.radicalElectrons, [r + 1 for r in atom.radicalElectrons])
                self.assertEqual(atom0.spinMultiplicity, [s + 1 for s in atom.spinMultiplicity])
                self.assertEqual(atom0.charge, atom.charge)
                self.assertEqual(atom0.label, atom.label)
            except ActionError:
                self.assertEqual(len(atomType.incrementRadical), 0)
    
    def testEquivalent(self):
        """
        Test the GroupAtom.equivalent() method.
        """
        for label1, atomType1 in atomTypes.iteritems():
            for label2, atomType2 in atomTypes.iteritems():
                atom1 = GroupAtom(atomType=[atomType1], radicalElectrons=[1], spinMultiplicity=[2], charge=[0], label='*1')
                atom2 = GroupAtom(atomType=[atomType2], radicalElectrons=[1], spinMultiplicity=[2], charge=[0], label='*1')
                if label1 == label2 or atomType2 in atomType1.generic or atomType1 in atomType2.generic:
                    self.assertTrue(atom1.equivalent(atom2))
                    self.assertTrue(atom2.equivalent(atom1))
                else:
                    self.assertFalse(atom1.equivalent(atom2))
                    self.assertFalse(atom2.equivalent(atom1))
    
    def testIsSpecificCaseOf(self):
        """
        Test the GroupAtom.isSpecificCaseOf() method.
        """
        for label1, atomType1 in atomTypes.iteritems():
            for label2, atomType2 in atomTypes.iteritems():
                atom1 = GroupAtom(atomType=[atomType1], radicalElectrons=[1], spinMultiplicity=[2], charge=[0], label='*1')
                atom2 = GroupAtom(atomType=[atomType2], radicalElectrons=[1], spinMultiplicity=[2], charge=[0], label='*1')
                if label1 == label2 or atomType2 in atomType1.generic:
                    self.assertTrue(atom1.isSpecificCaseOf(atom2))
                else:
                    self.assertFalse(atom1.isSpecificCaseOf(atom2))
    
    def testCopy(self):
        """
        Test the GroupAtom.copy() method.
        """
        atom = self.atom.copy()
        self.assertEqual(len(self.atom.atomType), len(atom.atomType))
        self.assertEqual(self.atom.atomType[0].label, atom.atomType[0].label)
        self.assertEqual(self.atom.radicalElectrons, atom.radicalElectrons)
        self.assertEqual(self.atom.spinMultiplicity, atom.spinMultiplicity)
        self.assertEqual(self.atom.charge, atom.charge)
        self.assertEqual(self.atom.label, atom.label)
    
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
        self.assertEqual(self.atom.spinMultiplicity, atom.spinMultiplicity)
        self.assertEqual(self.atom.charge, atom.charge)
        self.assertEqual(self.atom.label, atom.label)
        
    def testOutput(self):
        """
        Test that we can reconstruct a GroupAtom object from its repr()
        output with no loss of information.
        """
        exec('atom = {0!r}'.format(self.atom))
        self.assertEqual(len(self.atom.atomType), len(atom.atomType))
        self.assertEqual(self.atom.atomType[0].label, atom.atomType[0].label)
        self.assertEqual(self.atom.radicalElectrons, atom.radicalElectrons)
        self.assertEqual(self.atom.spinMultiplicity, atom.spinMultiplicity)
        self.assertEqual(self.atom.charge, atom.charge)
        self.assertEqual(self.atom.label, atom.label)

################################################################################

class TestGroupBond(unittest.TestCase):
    """
    Contains unit tests of the GroupBond class.
    """

    def setUp(self):
        """
        A method called before each unit test in this class.
        """
        self.bond = GroupBond(order=['D'])
        self.orderList = [['S'], ['D'], ['T'], ['B'], ['S','D'], ['D','S'], ['D','T'], ['S','D','T']]
    
    def testApplyActionBreakBond(self):
        """
        Test the GroupBond.applyAction() method for a BREAK_BOND action.
        """
        action = ['BREAK_BOND', '*1', 'S', '*2']
        for order0 in self.orderList:
            bond0 = GroupBond(order=order0)
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
        action = ['FORM_BOND', '*1', 'S', '*2']
        for order0 in self.orderList:
            bond0 = GroupBond(order=order0)
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
            bond0 = GroupBond(order=order0)
            bond = bond0.copy()
            try:
                bond.applyAction(action)
            except ActionError:
                self.assertTrue('T' in order0 or 'B' in order0)
                
    def testApplyActionDecrementBond(self):
        """
        Test the GroupBond.applyAction() method for a CHANGE_BOND action.
        """
        action = ['CHANGE_BOND', '*1', -1, '*2']
        for order0 in self.orderList:
            bond0 = GroupBond(order=order0)
            bond = bond0.copy()
            try:
                bond.applyAction(action)
            except ActionError:
                self.assertTrue('S' in order0 or 'B' in order0)
            
    def testApplyActionGainRadical(self):
        """
        Test the GroupBond.applyAction() method for a GAIN_RADICAL action.
        """
        action = ['GAIN_RADICAL', '*1', 1]
        for order0 in self.orderList:
            bond0 = GroupBond(order=order0)
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
            bond0 = GroupBond(order=order0)
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
                bond1 = GroupBond(order=order1)
                bond2 = GroupBond(order=order2)
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
                bond1 = GroupBond(order=order1)
                bond2 = GroupBond(order=order2)
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
        
    def testOutput(self):
        """
        Test that we can reconstruct a GroupBond object from its repr()
        output with no loss of information.
        """
        exec('bond = {0!r}'.format(self.bond))
        self.assertEqual(len(self.bond.order), len(bond.order))
        self.assertEqual(self.bond.order, bond.order)

################################################################################

class TestGroup(unittest.TestCase):
    """
    Contains unit tests of the Graph class.
    """

    def setUp(self):
        self.adjlist = """
1  *2 {Cs,Cd} 0 {2,{S,D}} {3,S}
2  *1 {Os,Od} 0 {1,{S,D}}
3     R!H 0 {1,S}
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
        bond12 = self.group.bonds[atom1][atom2]
        bond13 = self.group.bonds[atom1][atom3]
           
        self.assertTrue(atom1.label == '*2')
        self.assertTrue(atom1.atomType[0].label in ['Cs','Cd'])
        self.assertTrue(atom1.atomType[1].label in ['Cs','Cd'])
        self.assertTrue(atom1.radicalElectrons == [0])
        
        self.assertTrue(atom2.label == '*1')
        self.assertTrue(atom2.atomType[0].label in ['Os','Od'])
        self.assertTrue(atom2.atomType[1].label in ['Os','Od'])
        self.assertTrue(atom2.radicalElectrons == [0])
        
        self.assertTrue(atom3.label == '')
        self.assertTrue(atom3.atomType[0].label == 'R!H')
        self.assertTrue(atom3.radicalElectrons == [0])

        self.assertTrue(bond12.order == ['S','D'])
        self.assertTrue(bond13.order == ['S'])

    def testToAdjacencyList(self):
        """
        Test the Group.toAdjacencyList() method.
        """
        adjlist = self.group.toAdjacencyList()
        self.assertEqual(adjlist.strip(), self.adjlist.strip())

    def testIsIsomorphic(self):
        """
        Test the Group.isIsomorphic() method.
        """
        adjlist = """
1  *1 {Os,Od} 0 {3,{S,D}}
2     R!H 0 {3,S}
3  *2 {Cs,Cd} 0 {1,{S,D}} {2,S}
            """
        group = Group().fromAdjacencyList(adjlist)
        self.assertTrue(self.group.isIsomorphic(group))
        self.assertTrue(group.isIsomorphic(self.group))
        
    def testFindIsomorphism(self):
        """
        Test the Group.findIsomorphism() method.
        """
        adjlist = """
1  *1 {Os,Od} 0 {3,{S,D}}
2     R!H 0 {3,S}
3  *2 {Cs,Cd} 0 {1,{S,D}} {2,S}
            """
        group = Group().fromAdjacencyList(adjlist)
        result = self.group.findIsomorphism(group)
        self.assertTrue(result[0])
        self.assertEqual(len(result[1]), 1)
        for atom1, atom2 in result[1][0].iteritems():
            self.assertTrue(atom1 in self.group.atoms)
            self.assertTrue(atom2 in group.atoms)
            self.assertTrue(atom1.equivalent(atom2))
            for atom3 in self.group.bonds[atom1]:
                atom4 = result[1][0][atom3]
                self.assertTrue(atom4 in group.bonds[atom2])
                self.assertTrue(atom3.equivalent(atom4))
                bond1 = self.group.bonds[atom1][atom3]
                bond2 = group.bonds[atom2][atom4]
                self.assertTrue(bond1.equivalent(bond2))
        
    def testIsSubgraphIsomorphic(self):
        """
        Test the Group.isSubgraphIsomorphic() method.
        """
        adjlist = """
1  *1 {Cs,Cd} 0
            """
        group = Group().fromAdjacencyList(adjlist)
        self.assertTrue(self.group.isSubgraphIsomorphic(group))
        self.assertFalse(group.isIsomorphic(self.group))
        
    def testFindSubgraphIsomorphisms(self):
        """
        Test the Group.findSubgraphIsomorphisms() method.
        """
        adjlist = """
1  *1 {Cs,Cd} 0
            """
        group = Group().fromAdjacencyList(adjlist)
        result = self.group.findSubgraphIsomorphisms(group)
        self.assertTrue(result[0])
        self.assertEqual(len(result[1]), 1)
        for atom1, atom2 in result[1][0].iteritems():
            self.assertTrue(atom1 in self.group.atoms)
            self.assertTrue(atom2 in group.atoms)
            self.assertTrue(atom1.equivalent(atom2))
        
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
        
################################################################################

if __name__ == '__main__':
    unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )