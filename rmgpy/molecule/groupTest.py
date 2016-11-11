#!/usr/bin/python
# -*- coding: utf-8 -*-

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
        self.atom = GroupAtom(atomType=[atomTypes['Cd']], radicalElectrons=[1], charge=[0], label='*1')
    
    def testApplyActionBreakBond(self):
        """
        Test the GroupAtom.applyAction() method for a BREAK_BOND action.
        """
        action = ['BREAK_BOND', '*1', 'S', '*2']
        for label, atomType in atomTypes.iteritems():
            atom0 = GroupAtom(atomType=[atomType], radicalElectrons=[1], charge=[0], label='*1')
            atom = atom0.copy()
            try:
                atom.applyAction(action)
                self.assertEqual(len(atom.atomType), len(atomType.breakBond))
                for a in atomType.breakBond:
                    self.assertTrue(a in atom.atomType)
                self.assertEqual(atom0.radicalElectrons, atom.radicalElectrons)
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
            atom0 = GroupAtom(atomType=[atomType], radicalElectrons=[1], charge=[0], label='*1')
            atom = atom0.copy()
            try:
                atom.applyAction(action)
                self.assertEqual(len(atom.atomType), len(atomType.formBond))
                for a in atomType.formBond:
                    self.assertTrue(a in atom.atomType)
                self.assertEqual(atom0.radicalElectrons, atom.radicalElectrons)
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
            atom0 = GroupAtom(atomType=[atomType], radicalElectrons=[1], charge=[0], label='*1')
            atom = atom0.copy()
            try:
                atom.applyAction(action)
                self.assertEqual(len(atom.atomType), len(atomType.incrementBond))
                for a in atomType.incrementBond:
                    self.assertTrue(a in atom.atomType)
                self.assertEqual(atom0.radicalElectrons, atom.radicalElectrons)
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
            atom0 = GroupAtom(atomType=[atomType], radicalElectrons=[1], charge=[0], label='*1')
            atom = atom0.copy()
            try:
                atom.applyAction(action)
                self.assertEqual(len(atom.atomType), len(atomType.decrementBond))
                for a in atomType.decrementBond:
                    self.assertTrue(a in atom.atomType)
                self.assertEqual(atom0.radicalElectrons, atom.radicalElectrons)
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
            atom0 = GroupAtom(atomType=[atomType], radicalElectrons=[1], charge=[0], label='*1')
            atom = atom0.copy()
            try:
                atom.applyAction(action)
                self.assertEqual(len(atom.atomType), len(atomType.incrementRadical))
                for a in atomType.incrementRadical:
                    self.assertTrue(a in atom.atomType, "GAIN_RADICAL on {0} gave {1} not {2}".format(atomType, atom.atomType, atomType.incrementRadical))
                self.assertEqual(atom0.radicalElectrons, [r - 1 for r in atom.radicalElectrons])
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
            atom0 = GroupAtom(atomType=[atomType], radicalElectrons=[1], charge=[0], label='*1')
            atom = atom0.copy()
            try:
                atom.applyAction(action)
                self.assertEqual(len(atom.atomType), len(atomType.decrementRadical))
                for a in atomType.incrementRadical:
                    self.assertTrue(a in atom.atomType, "LOSE_RADICAL on {0} gave {1} not {2}".format(atomType, atom.atomType, atomType.decrementRadical))
                self.assertEqual(atom0.radicalElectrons, [r + 1 for r in atom.radicalElectrons])
                self.assertEqual(atom0.charge, atom.charge)
                self.assertEqual(atom0.label, atom.label)
            except ActionError:
                self.assertEqual(len(atomType.decrementRadical), 0)
    
    def testEquivalent(self):
        """
        Test the GroupAtom.equivalent() method.
        """
        for label1, atomType1 in atomTypes.iteritems():
            for label2, atomType2 in atomTypes.iteritems():
                atom1 = GroupAtom(atomType=[atomType1], radicalElectrons=[1], charge=[0], label='*1')
                atom2 = GroupAtom(atomType=[atomType2], radicalElectrons=[1], charge=[0], label='*1')
                if label1 == label2 or atomType2 in atomType1.generic or atomType1 in atomType2.generic:
                    self.assertTrue(atom1.equivalent(atom2), '{0!s} is not equivalent to {1!s}'.format(atom1, atom2))
                    self.assertTrue(atom2.equivalent(atom1), '{0!s} is not equivalent to {1!s}'.format(atom2, atom1))
                else:
                    self.assertFalse(atom1.equivalent(atom2), '{0!s} is equivalent to {1!s}'.format(atom1, atom2))
                    self.assertFalse(atom2.equivalent(atom1), '{0!s} is equivalent to {1!s}'.format(atom2, atom1))
            # Now see if charge and radical count are checked properly
            for charge in range(3):
                for radicals in range(2):
                    atom3 = GroupAtom(atomType=[atomType1], radicalElectrons=[radicals], charge=[charge], label='*1')
                    if radicals == 1 and charge == 0:
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
                atom1 = GroupAtom(atomType=[atomType1], radicalElectrons=[1], charge=[0], label='*1')
                atom2 = GroupAtom(atomType=[atomType2], radicalElectrons=[1], charge=[0], label='*1')
                # And make more generic types of these two atoms
                atom1gen = GroupAtom(atomType=[atomType1], radicalElectrons=[0, 1], charge=[0, 1], label='*1')
                atom2gen = GroupAtom(atomType=[atomType2], radicalElectrons=[0, 1], charge=[0, 1], label='*1')
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
################################################################################

class TestGroupBond(unittest.TestCase):
    """
    Contains unit tests of the GroupBond class.
    """

    def setUp(self):
        """
        A method called before each unit test in this class.
        """
        self.bond = GroupBond(None, None, order=['D'])
        self.orderList = [['S'], ['D'], ['T'], ['B'], ['S','D'], ['D','S'], ['D','T'], ['S','D','T']]
    
    def testApplyActionBreakBond(self):
        """
        Test the GroupBond.applyAction() method for a BREAK_BOND action.
        """
        action = ['BREAK_BOND', '*1', 'S', '*2']
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
        action = ['FORM_BOND', '*1', 'S', '*2']
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
                self.assertTrue('T' in order0 or 'B' in order0)
                
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
                self.assertTrue('S' in order0 or 'B' in order0)
            
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
1 *2 [Cs,Cd] u0 {2,[S,D]} {3,S}
2 *1 [Os,Od] u0 {1,[S,D]}
3    R!H     u0 {1,S}
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
        self.assertEqual(adjlist.strip(), self.adjlist.strip(),adjlist)

    def testIsIsomorphic(self):
        """
        Test the Group.isIsomorphic() method.
        """
        adjlist = """
1  *1 [Os,Od] u0 {3,[S,D]}
2     R!H     u0 {3,S}
3  *2 [Cs,Cd] u0 {1,[S,D]} {2,S}
            """
        group = Group().fromAdjacencyList(adjlist)
        self.assertTrue(self.group.isIsomorphic(group))
        self.assertTrue(group.isIsomorphic(self.group))
        
    def testFindIsomorphism(self):
        """
        Test the Group.findIsomorphism() method.
        """
        adjlist = """
1  *1 [Os,Od] u0 {3,[S,D]}
2     R!H     u0 {3,S}
3  *2 [Cs,Cd] u0 {1,[S,D]} {2,S}
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
2     C   u0 {1,T}
            """
        group5 = Group().fromAdjacencyList(adjlist5)
        group6 = Group().fromAdjacencyList(adjlist6)

        newGroup =group5.addImplicitAtomsFromAtomType()
        self.assertTrue(group6.isIsomorphic(newGroup))
        #test addition of lone pairs
        adjlist7 = """
1  *1 N1d u0
            """

        adjlist8 = """
1  *1 N1d u0 p2 {2,D}
2     C   u0 {1,D}
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
1  *1 C u0 {2,S} {3,D}
2     Ct u0 {1,S} {4,T}
3     C  u0 {1,D}
4     C  u0 {2,T}
            """
        group9 = Group().fromAdjacencyList(adjlist9)
        group10 = Group().fromAdjacencyList(adjlist10)

        newGroup =group9.addImplicitAtomsFromAtomType()
        self.assertTrue(group10.isIsomorphic(newGroup))

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
4  *4 O u0 p0 c0 {3,S}
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
12  C u0 {11,B} {7,B}
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

    def testMakeSampleMolecule(self):
        """
        Test the Group.makeSampleMolecule method
        """

        # result = self.group.makeSampleMolecule()
        # print result.multiplicity
        # self.assertTrue(result.isIsomorphic(Molecule().fromSMILES('OCC')))

        #tests adding implicit atoms
        adjlist1 = """
1  *1 Cd u0
            """

        group1 = Group().fromAdjacencyList(adjlist1)
        result1 = group1.makeSampleMolecule()
        self.assertTrue(result1.isIsomorphic(Molecule().fromSMILES('C=C')))

        #test creating implicit benzene atoms
        adjlist2 = """
1  *1 Cbf u0 {2,B}
2     Cbf u0 {1,B}
            """

        group2 = Group().fromAdjacencyList(adjlist2)
        result2 = group2.makeSampleMolecule()
        naphthaleneMolecule = Molecule().fromSMILES('C1=CC=C2C=CC=CC2=C1')
        resonanceList2=naphthaleneMolecule.generateResonanceIsomers()
        self.assertTrue(any([result2.isIsomorphic(x) for x in resonanceList2]))

        #test the creation of a charged species
        adjlist3 = """
1  *1 N5s u0
        """

        group3 = Group().fromAdjacencyList(adjlist3)
        result3 = group3.makeSampleMolecule()
        self.assertTrue(result3.isIsomorphic(Molecule().fromSMILES('[NH4+]')))

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


################################################################################

if __name__ == '__main__':
    unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )
