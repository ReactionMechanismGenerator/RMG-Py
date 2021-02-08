#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2020 Prof. William H. Green (whgreen@mit.edu),           #
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

import rmgpy.molecule.element as elements
from rmgpy.molecule import Molecule
from rmgpy.molecule.atomtype import ATOMTYPES
from rmgpy.molecule.group import ActionError, GroupAtom, GroupBond, Group

################################################################################


class TestGroupAtom(unittest.TestCase):
    """
    Contains unit tests of the GroupAtom class.
    """

    def setUp(self):
        """
        A method called before each unit test in this class.
        """
        self.atom = GroupAtom(atomtype=[ATOMTYPES['Cd']], radical_electrons=[1], charge=[0], label='*1', lone_pairs=[0])

    def test_apply_action_break_bond(self):
        """
        Test the GroupAtom.apply_action() method for a BREAK_BOND action.
        """
        action = ['BREAK_BOND', '*1', 1, '*2']
        for label, atomtype in ATOMTYPES.items():
            atom0 = GroupAtom(atomtype=[atomtype], radical_electrons=[1], charge=[0], label='*1', lone_pairs=[0])
            atom = atom0.copy()
            try:
                atom.apply_action(action)
                self.assertEqual(len(atom.atomtype), len(atomtype.break_bond))
                for a in atomtype.break_bond:
                    self.assertTrue(a in atom.atomtype)
                self.assertEqual(atom0.radical_electrons, atom.radical_electrons)
                self.assertEqual(atom0.charge, atom.charge)
                self.assertEqual(atom0.label, atom.label)
                self.assertEqual(atom0.label, atom.label)
                self.assertEqual(atom0.lone_pairs, atom.lone_pairs)
            except ActionError:
                self.assertEqual(len(atomtype.break_bond), 0)

    def test_apply_action_form_bond(self):
        """
        Test the GroupAtom.apply_action() method for a FORM_BOND action.
        """
        action = ['FORM_BOND', '*1', 1, '*2']
        for label, atomtype in ATOMTYPES.items():
            atom0 = GroupAtom(atomtype=[atomtype], radical_electrons=[1], charge=[0], label='*1', lone_pairs=[0])
            atom = atom0.copy()
            try:
                atom.apply_action(action)
                self.assertEqual(len(atom.atomtype), len(atomtype.form_bond))
                for a in atomtype.form_bond:
                    self.assertTrue(a in atom.atomtype)
                self.assertEqual(atom0.radical_electrons, atom.radical_electrons)
                self.assertEqual(atom0.charge, atom.charge)
                self.assertEqual(atom0.label, atom.label)
                self.assertEqual(atom0.lone_pairs, atom.lone_pairs)
            except ActionError:
                self.assertEqual(len(atomtype.form_bond), 0)

    def test_apply_action_increment_bond(self):
        """
        Test the GroupAtom.apply_action() method for a CHANGE_BOND action.
        """
        action = ['CHANGE_BOND', '*1', 1, '*2']
        for label, atomtype in ATOMTYPES.items():
            atom0 = GroupAtom(atomtype=[atomtype], radical_electrons=[1], charge=[0], label='*1', lone_pairs=[0])
            atom = atom0.copy()
            try:
                atom.apply_action(action)
                self.assertEqual(len(atom.atomtype), len(atomtype.increment_bond))
                for a in atomtype.increment_bond:
                    self.assertTrue(a in atom.atomtype)
                self.assertEqual(atom0.radical_electrons, atom.radical_electrons)
                self.assertEqual(atom0.charge, atom.charge)
                self.assertEqual(atom0.label, atom.label)
                self.assertEqual(atom0.lone_pairs, atom.lone_pairs)
            except ActionError:
                self.assertEqual(len(atomtype.increment_bond), 0)

    def test_apply_action_decrement_bond(self):
        """
        Test the GroupAtom.apply_action() method for a CHANGE_BOND action.
        """
        action = ['CHANGE_BOND', '*1', -1, '*2']
        for label, atomtype in ATOMTYPES.items():
            atom0 = GroupAtom(atomtype=[atomtype], radical_electrons=[1], charge=[0], label='*1', lone_pairs=[0])
            atom = atom0.copy()
            try:
                atom.apply_action(action)
                self.assertEqual(len(atom.atomtype), len(atomtype.decrement_bond))
                for a in atomtype.decrement_bond:
                    self.assertTrue(a in atom.atomtype)
                self.assertEqual(atom0.radical_electrons, atom.radical_electrons)
                self.assertEqual(atom0.charge, atom.charge)
                self.assertEqual(atom0.label, atom.label)
                self.assertEqual(atom0.lone_pairs, atom.lone_pairs)
            except ActionError:
                self.assertEqual(len(atomtype.decrement_bond), 0)

    def test_apply_action_gain_radical(self):
        """
        Test the GroupAtom.apply_action() method for a GAIN_RADICAL action.
        """
        action = ['GAIN_RADICAL', '*1', 1]
        for label, atomtype in ATOMTYPES.items():
            atom0 = GroupAtom(atomtype=[atomtype], radical_electrons=[1], charge=[0], label='*1', lone_pairs=[0])
            atom = atom0.copy()
            try:
                atom.apply_action(action)
                self.assertEqual(len(atom.atomtype), len(atomtype.increment_radical))
                for a in atomtype.increment_radical:
                    self.assertTrue(a in atom.atomtype,
                                    "GAIN_RADICAL on {0} gave {1} not {2}".format(atomtype, atom.atomtype,
                                                                                  atomtype.increment_radical))
                self.assertEqual(atom0.radical_electrons, [r - 1 for r in atom.radical_electrons])
                self.assertEqual(atom0.charge, atom.charge)
                self.assertEqual(atom0.label, atom.label)
                self.assertEqual(atom0.lone_pairs, atom.lone_pairs)
            except ActionError:
                self.assertEqual(len(atomtype.increment_radical), 0)

        # test when radicals unspecified
        group = Group().from_adjacency_list("""
        1 R ux
        """)  # ux causes a wildcard for radicals
        atom1 = group.atoms[0]
        atom1.apply_action(action)
        self.assertListEqual(atom1.radical_electrons, [1, 2, 3, 4])

    def test_apply_action_lose_radical(self):
        """
        Test the GroupAtom.apply_action() method for a LOSE_RADICAL action.
        """
        action = ['LOSE_RADICAL', '*1', 1]
        for label, atomtype in ATOMTYPES.items():
            atom0 = GroupAtom(atomtype=[atomtype], radical_electrons=[1], charge=[0], label='*1', lone_pairs=[0])
            atom = atom0.copy()
            try:
                atom.apply_action(action)
                self.assertEqual(len(atom.atomtype), len(atomtype.decrement_radical))
                for a in atomtype.increment_radical:
                    self.assertTrue(a in atom.atomtype,
                                    "LOSE_RADICAL on {0} gave {1} not {2}".format(atomtype, atom.atomtype,
                                                                                  atomtype.decrement_radical))
                self.assertEqual(atom0.radical_electrons, [r + 1 for r in atom.radical_electrons])
                self.assertEqual(atom0.charge, atom.charge)
                self.assertEqual(atom0.label, atom.label)
                self.assertEqual(atom0.lone_pairs, atom.lone_pairs)
            except ActionError:
                self.assertEqual(len(atomtype.decrement_radical), 0)

        # test when radicals unspecified
        group = Group().from_adjacency_list("""
        1 R ux
        """)  # ux causes a wildcard for radicals
        atom1 = group.atoms[0]
        atom1.apply_action(action)
        self.assertListEqual(atom1.radical_electrons, [0, 1, 2, 3])

    def test_apply_action_gain_pair(self):
        """
        Test the GroupAtom.apply_action() method for a GAIN_PAIR action when lone_pairs is either specified or not.
        """

        action = ['GAIN_PAIR', '*1', 1]

        # lone_pairs specified:
        for label, atomtype in ATOMTYPES.items():
            atom0 = GroupAtom(atomtype=[atomtype], radical_electrons=[1], charge=[0], label='*1', lone_pairs=[0])
            atom = atom0.copy()
            try:
                atom.apply_action(action)
                self.assertEqual(len(atom.atomtype), len(atomtype.increment_lone_pair))
                for a in atomtype.increment_lone_pair:
                    self.assertTrue(a in atom.atomtype,
                                    "GAIN_PAIR on {0} gave {1} not {2}".format(atomtype, atom.atomtype,
                                                                               atomtype.increment_lone_pair))
                self.assertEqual(atom0.radical_electrons, atom.radical_electrons)
                self.assertEqual(atom0.charge, atom.charge)
                self.assertEqual(atom0.label, atom.label)
                self.assertEqual(atom0.lone_pairs, [r - 1 for r in atom.lone_pairs])
            except ActionError:
                self.assertEqual(len(atomtype.increment_lone_pair), 0)

        # lone_pairs unspecified:
        for label, atomtype in ATOMTYPES.items():
            atom0 = GroupAtom(atomtype=[atomtype], radical_electrons=[1], charge=[0], label='*1')
            atom = atom0.copy()
            try:
                atom.apply_action(action)
                self.assertEqual(len(atom.atomtype), len(atomtype.increment_lone_pair))
                for a in atomtype.increment_lone_pair:
                    self.assertTrue(a in atom.atomtype,
                                    "GAIN_PAIR on {0} gave {1} not {2}".format(atomtype, atom.atomtype,
                                                                               atomtype.increment_lone_pair))
                self.assertEqual(atom0.radical_electrons, atom.radical_electrons)
                self.assertEqual(atom0.charge, atom.charge)
                self.assertEqual(atom0.label, atom.label)
                self.assertEqual([0, 1, 2, 3], [r - 1 for r in atom.lone_pairs])
            except ActionError:
                self.assertEqual(len(atomtype.increment_lone_pair), 0)

    def test_apply_action_lose_pair(self):
        """
        Test the GroupAtom.apply_action() method for a LOSE_PAIR action when lone_pairs is either specified or not.
        """

        action = ['LOSE_PAIR', '*1', 1]

        # lone_pairs specified:
        for label, atomtype in ATOMTYPES.items():
            atom0 = GroupAtom(atomtype=[atomtype], radical_electrons=[1], charge=[0], label='*1', lone_pairs=[1])
            atom = atom0.copy()
            try:
                atom.apply_action(action)
                self.assertEqual(len(atom.atomtype), len(atomtype.decrement_lone_pair))
                for a in atomtype.decrement_lone_pair:
                    self.assertTrue(a in atom.atomtype,
                                    "LOSE_PAIR on {0} gave {1} not {2}".format(atomtype, atom.atomtype,
                                                                               atomtype.decrement_lone_pair))
                self.assertEqual(atom0.radical_electrons, atom.radical_electrons)
                self.assertEqual(atom0.charge, atom.charge)
                self.assertEqual(atom0.label, atom.label)
                self.assertEqual(atom0.lone_pairs, [r + 1 for r in atom.lone_pairs])
            except ActionError:
                self.assertEqual(len(atomtype.decrement_lone_pair), 0)

        # lone_pairs unspecified:
        for label, atomtype in ATOMTYPES.items():
            atom0 = GroupAtom(atomtype=[atomtype], radical_electrons=[1], charge=[0], label='*1')
            atom = atom0.copy()
            try:
                atom.apply_action(action)
                self.assertEqual(len(atom.atomtype), len(atomtype.decrement_lone_pair))
                for a in atomtype.decrement_lone_pair:
                    self.assertTrue(a in atom.atomtype,
                                    "LOSE_PAIR on {0} gave {1} not {2}".format(atomtype, atom.atomtype,
                                                                               atomtype.decrement_lone_pair))
                self.assertEqual(atom0.radical_electrons, atom.radical_electrons)
                self.assertEqual(atom0.charge, atom.charge)
                self.assertEqual(atom0.label, atom.label)
                self.assertEqual([1, 2, 3, 4], [r + 1 for r in atom.lone_pairs])
            except ActionError:
                self.assertEqual(len(atomtype.decrement_lone_pair), 0)

    def test_equivalent(self):
        """
        Test the GroupAtom.equivalent() method.
        """
        for label1, atomType1 in ATOMTYPES.items():
            for label2, atomType2 in ATOMTYPES.items():
                atom1 = GroupAtom(atomtype=[atomType1], radical_electrons=[1], charge=[0], label='*1', lone_pairs=[0])
                atom2 = GroupAtom(atomtype=[atomType2], radical_electrons=[1], charge=[0], label='*1', lone_pairs=[0])
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
                        atom3 = GroupAtom(atomtype=[atomType1], radical_electrons=[radicals], charge=[charge],
                                          label='*1', lone_pairs=[lonePair])
                        if radicals == 1 and charge == 0 and lonePair == 0:
                            self.assertTrue(atom1.equivalent(atom3),
                                            '{0!s} is not equivalent to {1!s}'.format(atom1, atom3))
                            self.assertTrue(atom1.equivalent(atom3),
                                            '{0!s} is not equivalent to {1!s}'.format(atom3, atom1))
                        else:
                            self.assertFalse(atom1.equivalent(atom3),
                                             '{0!s} is equivalent to {1!s}'.format(atom1, atom3))
                            self.assertFalse(atom1.equivalent(atom3),
                                             '{0!s} is equivalent to {1!s}'.format(atom3, atom1))

    def test_is_specific_case_of(self):
        """
        Test the GroupAtom.is_specific_case_of() method.
        """
        for label1, atomType1 in ATOMTYPES.items():
            for label2, atomType2 in ATOMTYPES.items():
                atom1 = GroupAtom(atomtype=[atomType1], radical_electrons=[1], charge=[0], label='*1', lone_pairs=[0])
                atom2 = GroupAtom(atomtype=[atomType2], radical_electrons=[1], charge=[0], label='*1', lone_pairs=[0])
                # And make more generic types of these two atoms
                atom1gen = GroupAtom(atomtype=[atomType1], radical_electrons=[0, 1], charge=[0, 1], label='*1',
                                     lone_pairs=[0, 1])
                atom2gen = GroupAtom(atomtype=[atomType2], radical_electrons=[0, 1], charge=[0, 1], label='*1',
                                     lone_pairs=[0, 1])
                if label1 == label2 or atomType2 in atomType1.generic:
                    self.assertTrue(atom1.is_specific_case_of(atom2),
                                    '{0!s} is not a specific case of {1!s}'.format(atom1, atom2))
                    self.assertTrue(atom1.is_specific_case_of(atom2gen),
                                    '{0!s} is not a specific case of {1!s}'.format(atom1, atom2gen))
                    self.assertFalse(atom1gen.is_specific_case_of(atom2),
                                     '{0!s} is a specific case of {1!s}'.format(atom1gen, atom2))
                else:
                    self.assertFalse(atom1.is_specific_case_of(atom2),
                                     '{0!s} is a specific case of {1!s}'.format(atom1, atom2))
                    self.assertFalse(atom1.is_specific_case_of(atom2gen),
                                     '{0!s} is a specific case of {1!s}'.format(atom1, atom2gen))
                    self.assertFalse(atom1gen.is_specific_case_of(atom2),
                                     '{0!s} is a specific case of {1!s}'.format(atom1gen, atom2))

    def test_copy(self):
        """
        Test the GroupAtom.copy() method.
        """
        atom = self.atom.copy()
        self.assertEqual(len(self.atom.atomtype), len(atom.atomtype))
        self.assertEqual(self.atom.atomtype[0].label, atom.atomtype[0].label)
        self.assertEqual(self.atom.radical_electrons, atom.radical_electrons)
        self.assertEqual(self.atom.charge, atom.charge)
        self.assertEqual(self.atom.label, atom.label)
        self.assertEqual(self.atom.lone_pairs, atom.lone_pairs)

    def test_pickle(self):
        """
        Test that a GroupAtom object can be successfully pickled and
        unpickled with no loss of information.
        """
        import pickle
        atom = pickle.loads(pickle.dumps(self.atom))
        self.assertEqual(len(self.atom.atomtype), len(atom.atomtype))
        self.assertEqual(self.atom.atomtype[0].label, atom.atomtype[0].label)
        self.assertEqual(self.atom.radical_electrons, atom.radical_electrons)
        self.assertEqual(self.atom.charge, atom.charge)
        self.assertEqual(self.atom.label, atom.label)
        self.assertEqual(self.atom.lone_pairs, atom.lone_pairs)

    def test_count_bonds(self):
        """
        Tests the count_bonds function
        """
        adjlist = """
1 *2 C u0     {2,[D,T]} {3,S}
2 *3 C u0     {1,[D,T]} {4,B}
3    C ux     {1,S} {5,D}
4    C u[0,1] {2,B}
5    O u0     {3,D}
6    C u0     {7,Q}
7    C u0     {6,Q}
"""
        test = Group().from_adjacency_list(adjlist)
        # returns a list of [single, allDouble, rDouble, oDouble, sDouble, triple, quadruple, benzene]
        self.assertListEqual([1, 0, 0, 0, 0, 0, 0, 0], test.atoms[0].count_bonds())
        self.assertListEqual([1, 1, 1, 0, 0, 1, 0, 0], test.atoms[0].count_bonds(wildcards=True))
        self.assertListEqual([0, 0, 0, 0, 0, 0, 0, 1], test.atoms[3].count_bonds())
        self.assertListEqual([1, 1, 0, 1, 0, 0, 0, 0], test.atoms[2].count_bonds())
        self.assertListEqual([0, 0, 0, 0, 0, 0, 1, 0], test.atoms[5].count_bonds())

    def test_has_wildcards(self):
        """
        Tests the GroupAtom.has_wildcards() method
        """
        self.assertFalse(self.atom.has_wildcards())
        adjlist = """
1 *2 C     u0     {2,[D,T]} {3,S}
2 *3 C     u0     {1,[D,T]} {4,S}
3    C     ux     {1,S} {5,S}
4    C     u[0,1] {2,S}
5    [C,O] u0     {3,S}
"""
        group = Group().from_adjacency_list(adjlist)
        for index, atom in enumerate(group.atoms):
            self.assertTrue(atom.has_wildcards(),
                            'GroupAtom with index {0} should have wildcards, but does not'.format(index))

    def test_make_sample_atom(self):
        """
        Tests the GroupAtom.make_sample_atom() method
        """
        new_atom = self.atom.make_sample_atom()

        self.assertEquals(new_atom.element, elements.__dict__['C'])
        self.assertEquals(new_atom.radical_electrons, 1)
        self.assertEquals(new_atom.charge, 0)
        self.assertEquals(new_atom.lone_pairs, 0)


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
        self.orderList = [[1], [2], [3], [1.5], [1, 2], [2, 1], [2, 3], [1, 2, 3]]  # todo : unit tests for vdw

    def test_get_order_str(self):
        """
        test the Bond.get_order_str() method
        """
        bond = GroupBond(None, None, order=[1, 2, 3, 1.5])
        self.assertEqual(bond.get_order_str(), ['S', 'D', 'T', 'B'])

    def test_set_order_str(self):
        """
        test the Bond.set_order_str() method
        """

        self.bond.set_order_str(["B", 'T'])
        self.assertEqual(set(self.bond.order), {3, 1.5})

    def test_get_order_num(self):
        """
        test the Bond.get_order_num() method
        """
        self.assertEqual(self.bond.get_order_num(), [2])

    def test_set_order_num(self):
        """
        test the Bond.set_order_num() method
        """

        self.bond.set_order_num([3, 1, 2])
        self.assertEqual(self.bond.get_order_str(), ['T', 'S', 'D'])

    def test_is_single(self):
        """
        test the Bond.is_single() method
        """
        self.bond.set_order_num([1])
        self.assertTrue(self.bond.is_single())

        # test interaction with wildcards
        self.bond.set_order_num([1, 2, 3, 1.5])
        self.assertFalse(self.bond.is_single(wildcards=False))
        self.assertTrue(self.bond.is_single(wildcards=True))

    def test_is_double(self):
        """
        test the Bond.is_double() method
        """
        self.bond.set_order_num([2])
        self.assertTrue(self.bond.is_double())

        # test interaction with wildcards
        self.bond.set_order_num([1, 2, 3, 1.5])
        self.assertFalse(self.bond.is_double(wildcards=False))
        self.assertTrue(self.bond.is_double(wildcards=True))

    def test_is_triple(self):
        """
        test the Bond.is_triple() method
        """
        self.bond.set_order_num([3])
        self.assertTrue(self.bond.is_triple())

        # test interaction with wildcards
        self.bond.set_order_num([1, 2, 3, 1.5])
        self.assertFalse(self.bond.is_triple(wildcards=False))
        self.assertTrue(self.bond.is_triple(wildcards=True))

    def test_is_benzene(self):
        """
        test the Bond.is_benzene() method
        """
        self.bond.set_order_num([1.5])
        self.assertTrue(self.bond.is_benzene())

        # test interaction with wildcards
        self.bond.set_order_num([1, 2, 3, 1.5])
        self.assertFalse(self.bond.is_benzene(wildcards=False))
        self.assertTrue(self.bond.is_benzene(wildcards=True))

    def test_apply_action_break_bond(self):
        """
        Test the GroupBond.apply_action() method for a BREAK_BOND action.
        """
        action = ['BREAK_BOND', '*1', 1, '*2']
        for order0 in self.orderList:
            bond0 = GroupBond(None, None, order=order0)
            bond = bond0.copy()
            try:
                bond.apply_action(action)
                self.fail('GroupBond.apply_action() unexpectedly processed a BREAK_BOND action.')
            except ActionError:
                pass

    def test_apply_action_form_bond(self):
        """
        Test the GroupBond.apply_action() method for a FORM_BOND action.

        Tests that forming a bond between things already bonded, raises
        an ActionError
        """
        action = ['FORM_BOND', '*1', 1, '*2']
        for order0 in self.orderList:
            bond0 = GroupBond(None, None, order=order0)
            bond = bond0.copy()
            try:
                bond.apply_action(action)
                self.fail('GroupBond.apply_action() unexpectedly processed a FORM_BOND action.')
            except ActionError:
                pass

    def test_apply_action_increment_bond(self):
        """
        Test the GroupBond.apply_action() method for a CHANGE_BOND action.
        """
        action = ['CHANGE_BOND', '*1', 1, '*2']
        for order0 in self.orderList:
            bond0 = GroupBond(None, None, order=order0)
            bond = bond0.copy()
            try:
                bond.apply_action(action)
            except ActionError:
                self.assertTrue(3 in order0 or 1.5 in order0)

    def test_apply_action_decrement_bond(self):
        """
        Test the GroupBond.apply_action() method for a CHANGE_BOND action.
        """
        action = ['CHANGE_BOND', '*1', -1, '*2']
        for order0 in self.orderList:
            bond0 = GroupBond(None, None, order=order0)
            bond = bond0.copy()
            try:
                bond.apply_action(action)
            except ActionError:
                self.assertTrue(1 in order0 or 1.5 in order0)

    def test_apply_action_gain_radical(self):
        """
        Test the GroupBond.apply_action() method for a GAIN_RADICAL action.
        """
        action = ['GAIN_RADICAL', '*1', 1]
        for order0 in self.orderList:
            bond0 = GroupBond(None, None, order=order0)
            bond = bond0.copy()
            try:
                bond.apply_action(action)
                self.fail('GroupBond.apply_action() unexpectedly processed a GAIN_RADICAL action.')
            except ActionError:
                pass

    def test_apply_action_lose_radical(self):
        """
        Test the GroupBond.apply_action() method for a LOSE_RADICAL action.
        """
        action = ['LOSE_RADICAL', '*1', 1]
        for order0 in self.orderList:
            bond0 = GroupBond(None, None, order=order0)
            bond = bond0.copy()
            try:
                bond.apply_action(action)
                self.fail('GroupBond.apply_action() unexpectedly processed a LOSE_RADICAL action.')
            except ActionError:
                pass

    def test_equivalent(self):
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

    def test_is_specific_case_of(self):
        """
        Test the GroupBond.is_specific_case_of() method.
        """
        for order1 in self.orderList:
            for order2 in self.orderList:
                bond1 = GroupBond(None, None, order=order1)
                bond2 = GroupBond(None, None, order=order2)
                if order1 == order2 or all([o in order2 for o in order1]):
                    self.assertTrue(bond1.is_specific_case_of(bond2))
                else:
                    self.assertFalse(bond1.is_specific_case_of(bond2))

    def test_copy(self):
        """
        Test the GroupBond.copy() method.
        """
        bond = self.bond.copy()
        self.assertEqual(len(self.bond.order), len(bond.order))
        self.assertEqual(self.bond.order, bond.order)

    def test_pickle(self):
        """
        Test that a GroupBond object can be successfully pickled and
        unpickled with no loss of information.
        """
        import pickle
        bond = pickle.loads(pickle.dumps(self.bond))
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
        self.group = Group().from_adjacency_list(self.adjlist)

    def test_clear_labeled_atoms(self):
        """
        Test the Group.clear_labeled_atoms() method.
        """
        self.group.clear_labeled_atoms()
        for atom in self.group.atoms:
            self.assertEqual(atom.label, '')

    def test_contains_labeled_atom(self):
        """
        Test the Group.contains_labeled_atom() method.
        """
        for atom in self.group.atoms:
            if atom.label != '':
                self.assertTrue(self.group.contains_labeled_atom(atom.label))
        self.assertFalse(self.group.contains_labeled_atom('*3'))
        self.assertFalse(self.group.contains_labeled_atom('*4'))
        self.assertFalse(self.group.contains_labeled_atom('*5'))
        self.assertFalse(self.group.contains_labeled_atom('*6'))

    def test_contains_surface_site(self):
        """
        Test the Group.contains_surface_site() method.
        """
        self.assertFalse(self.group.contains_surface_site())
        surface_group = Group().from_adjacency_list("""
1 *1 X u0 {2,[S,D]}
2 *2 R u0 {1,[S,D]}
""")
        self.assertTrue(surface_group.contains_surface_site())

    def test_is_surface_site(self):
        """
        Test the Group.is_surface_site() method.
        """
        self.assertFalse(self.group.is_surface_site())
        surface_group = Group().from_adjacency_list("""
1 *1 X u0 {2,[S,D]}
2 *2 R u0 {1,[S,D]}
""")
        self.assertFalse(surface_group.is_surface_site())
        surface_site = Group().from_adjacency_list("1 *1 X u0")
        self.assertTrue(surface_site.is_surface_site())

    def test_get_labeled_atom(self):
        """
        Test the Group.get_labeled_atoms() method.
        """
        for atom in self.group.atoms:
            if atom.label != '':
                self.assertEqual(atom, self.group.get_labeled_atoms(atom.label)[0])
        try:
            self.group.get_labeled_atoms('*3')[0]
            self.fail('Unexpected successful return from Group.get_labeled_atoms() with invalid atom label.')
        except ValueError:
            pass

    def test_get_labeled_atoms(self):
        """
        Test the Group.get_all_labeled_atoms() method.
        """
        labeled = self.group.get_all_labeled_atoms()
        for atom in self.group.atoms:
            if atom.label != '':
                self.assertTrue(atom.label in labeled)
                self.assertTrue(atom in list(labeled.values()))
            else:
                self.assertFalse(atom.label in labeled)
                self.assertFalse(atom in list(labeled.values()))

    def test_from_adjacency_list(self):
        """
        Test the Group.from_adjacency_list() method.
        """
        atom1, atom2, atom3 = self.group.atoms
        self.assertTrue(self.group.has_bond(atom1, atom2))
        self.assertTrue(self.group.has_bond(atom1, atom3))
        self.assertFalse(self.group.has_bond(atom2, atom3))
        bond12 = atom1.bonds[atom2]
        bond13 = atom1.bonds[atom3]

        self.assertTrue(atom1.label == '*2')
        self.assertTrue(atom1.atomtype[0].label in ['Cs', 'Cd'])
        self.assertTrue(atom1.atomtype[1].label in ['Cs', 'Cd'])
        self.assertTrue(atom1.radical_electrons == [0])

        self.assertTrue(atom2.label == '*1')
        self.assertTrue(atom2.atomtype[0].label in ['O2s', 'O2d'])
        self.assertTrue(atom2.atomtype[1].label in ['O2s', 'O2d'])
        self.assertTrue(atom2.radical_electrons == [0])

        self.assertTrue(atom3.label == '')
        self.assertTrue(atom3.atomtype[0].label == 'R!H')
        self.assertTrue(atom3.radical_electrons == [0])

        self.assertTrue(bond12.order == [1, 2])
        self.assertTrue(bond13.is_single())

    def test_to_adjacency_list(self):
        """
        Test the Group.to_adjacency_list() method.
        """
        adjlist = self.group.to_adjacency_list()
        self.assertEqual(adjlist.strip(), self.adjlist.strip(), adjlist)

    def test_is_isomorphic(self):
        """
        Test the Group.is_isomorphic() method.
        """
        adjlist = """
1  *1 [O2s,O2d] u0 {3,[S,D]}
2     R!H       u0 {3,S}
3  *2 [Cs,Cd]   u0 {1,[S,D]} {2,S}
"""
        group = Group().from_adjacency_list(adjlist)
        self.assertTrue(self.group.is_isomorphic(group))
        self.assertTrue(group.is_isomorphic(self.group))

    def test_find_isomorphism(self):
        """
        Test the Group.find_isomorphism() method.
        """
        adjlist = """
1  *1 [O2s,O2d] u0 {3,[S,D]}
2     R!H       u0 {3,S}
3  *2 [Cs,Cd]   u0 {1,[S,D]} {2,S}
"""
        group = Group().from_adjacency_list(adjlist)
        result = self.group.find_isomorphism(group)
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

    def test_is_subgraph_isomorphic(self):
        """
        Test the Group.is_subgraph_isomorphic() method.
        """
        adjlist = """
1  *1 [Cs,Cd] u0
"""
        group = Group().from_adjacency_list(adjlist)
        self.assertTrue(self.group.is_subgraph_isomorphic(group))
        self.assertFalse(group.is_isomorphic(self.group))

    def test_find_subgraph_isomorphisms(self):
        """
        Test the Group.find_subgraph_isomorphisms() method.
        """
        adjlist = """
1  *1 [Cs,Cd] u0
            """
        group = Group().from_adjacency_list(adjlist)
        result = self.group.find_subgraph_isomorphisms(group)
        self.assertEqual(len(result), 1)
        for atom1, atom2 in result[0].items():
            self.assertTrue(atom1 in self.group.atoms)
            self.assertTrue(atom2 in group.atoms)
            self.assertTrue(atom1.equivalent(atom2))

    def test_generate_extensions(self):
        """
        test that appropriate group extensions are being generated
        """

        test_grp = Group().from_adjacency_list("""
1 *2 C u0 r0 {2,[S,D]} 
2 *1 C u[0,1] {1,[S,D]} {3,S}
3    R!H     u0 r1 {2,S}
""")

        ans = [
            '1 *2 C   u0     r0 {2,[S,D]} {4,[S,D,T,B]}\n2 *1 C   u[0,1] {1,[S,D]} {3,S}\n3    R!H u0     r1 {2,S}\n4    R!H ux     {1,[S,D,T,B]}\n',
            '1 *2 C   u0 r0 {2,[S,D]}\n2 *1 C   u0 {1,[S,D]} {3,S}\n3    R!H u0 r1 {2,S}\n',
            '1 *2 C   u0 r0 {2,[S,D]}\n2 *1 C   u1 {1,[S,D]} {3,S}\n3    R!H u0 r1 {2,S}\n',
            '1 *2 C   u0     r0 {2,[S,D]}\n2 *1 C   u[0,1] r1 {1,[S,D]} {3,S}\n3    R!H u0     r1 {2,S}\n',
            '1 *2 C   u0     r0 {2,[S,D]}\n2 *1 C   u[0,1] {1,[S,D]} {3,S} {4,[S,D,T,B]}\n3    R!H u0     r1 {2,S}\n4    R!H ux     {2,[S,D,T,B]}\n',
            '1 *2 C   u0     r0 {2,S}\n2 *1 C   u[0,1] {1,S} {3,S}\n3    R!H u0     r1 {2,S}\n',
            '1 *2 C   u0     r0 {2,D}\n2 *1 C   u[0,1] {1,D} {3,S}\n3    R!H u0     r1 {2,S}\n',
            '1 *2 C u0     r0 {2,[S,D]}\n2 *1 C u[0,1] {1,[S,D]} {3,S}\n3    C u0     r1 {2,S}\n',
            '1 *2 C u0     r0 {2,[S,D]}\n2 *1 C u[0,1] {1,[S,D]} {3,S}\n3    O u0     r1 {2,S}\n',
            '1 *2 C   u0     r0 {2,[S,D]}\n2 *1 C   u[0,1] {1,[S,D]} {3,S}\n3    R!H u0     r1 {2,S} {4,[S,D,T,B]}\n4    R!H ux     {3,[S,D,T,B]}\n',
            '1 *2 C   u0     r0 {2,[S,D]} {3,[S,D,T,B]}\n2 *1 C   u[0,1] {1,[S,D]} {3,S}\n3    R!H u0     r1 {1,[S,D,T,B]} {2,S}\n'
        ]
        ans = [Group().from_adjacency_list(k) for k in ans]

        extensions = test_grp.get_extensions(r=[ATOMTYPES[i] for i in ['C', 'O', 'H']])
        extensions = [a[0] for a in extensions]

        self.assertEqual(len(extensions), len(ans))

        for v in ans:
            boos = [ext.is_identical(v) and ext.is_subgraph_isomorphic(v, generate_initial_map=True) for ext in extensions]
            self.assertTrue(any(boos), 'generated extensions did not match expected extensions')

    def test_generated_extensions_subgraphs(self):
        test_grp = Group().from_adjacency_list("""
1 *2 C u0 {2,[S,D]} 
2 *1 C u[0,1] {1,[S,D]} {3,S}
3    R!H     u0 {2,S}
""")
        extensions = test_grp.get_extensions(r=[ATOMTYPES[i] for i in ['C', 'O', 'H']])
        extensions = [a[0] for a in extensions]

        for ext in extensions:
            self.assertTrue(ext.is_subgraph_isomorphic(test_grp, generate_initial_map=True))

    def test_pickle(self):
        """
        Test that a Group object can be successfully pickled and
        unpickled with no loss of information.
        """
        import pickle
        group = pickle.loads(pickle.dumps(self.group))

        self.assertEqual(len(self.group.atoms), len(group.atoms))
        for atom0, atom in zip(group.atoms, self.group.atoms):
            self.assertTrue(atom0.equivalent(atom))

        self.assertTrue(self.group.is_isomorphic(group))
        self.assertTrue(group.is_isomorphic(self.group))

    def test_create_and_connect_atom(self):
        """
        Tests create_and_connect_atom method
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

        group1 = Group().from_adjacency_list(adjlist1)
        answer1 = Group().from_adjacency_list(answer1)
        atom1 = group1.get_labeled_atoms("*1")[0]
        new_atom = group1.create_and_connect_atom(atomtypes=["Cb"], connecting_atom=atom1, bond_orders=["B"])
        self.assertTrue(group1.is_isomorphic(answer1))

        answer2 = """
1  *1 C       u0 {2,S} {3,[S,D]}
2  *2 C       u0 {1,S}
3     [Cs,Cd] u0 {1,[S,D]}
"""

        # Test that wildcards work alright
        group2 = Group().from_adjacency_list(adjlist1)
        answer2 = Group().from_adjacency_list(answer2)
        atom1 = group2.get_labeled_atoms("*1")[0]
        new_atom = group2.create_and_connect_atom(atomtypes=["Cs", "Cd"], connecting_atom=atom1, bond_orders=["S", "D"])
        self.assertTrue(group2.is_isomorphic(answer2))

    def test_add_implicit_atoms_from_atom_type(self):
        """
        test Group.add_implicit_atoms_from_atomtype() method
        """
        # basic test adding oDouble
        adjlist1 = """
1  *1 CO u0
"""

        adjlist2 = """
1  *1 CO u0 {2,D}
2     O  u0 {1,D}
"""

        group1 = Group().from_adjacency_list(adjlist1)
        group2 = Group().from_adjacency_list(adjlist2)

        new_group = group1.add_implicit_atoms_from_atomtype()
        self.assertTrue(group2.is_isomorphic(new_group))
        # testing the allDouble match (more complicated
        adjlist3 = """
1  *1 Cdd u0
"""

        adjlist4 = """
1  *1 Cdd u0 {2,D} {3,D}
2     C   u0 {1,D}
3     C   u0 {1,D}
"""
        group3 = Group().from_adjacency_list(adjlist3)
        group4 = Group().from_adjacency_list(adjlist4)

        new_group = group3.add_implicit_atoms_from_atomtype()
        self.assertTrue(group4.is_isomorphic(new_group))
        # test adding a triple bond
        adjlist5 = """
1  *1 Ct u0
"""

        adjlist6 = """
1  *1 Ct u0 {2,T}
2     C  u0 {1,T}
"""
        group5 = Group().from_adjacency_list(adjlist5)
        group6 = Group().from_adjacency_list(adjlist6)

        new_group = group5.add_implicit_atoms_from_atomtype()
        self.assertTrue(group6.is_isomorphic(new_group))
        # test addition of lone pairs
        adjlist7 = """
1  *1 N1dc u0
"""

        adjlist8 = """
1  *1 N1dc u0 p2 {2,D}
2     C    u0 {1,D}
"""
        group7 = Group().from_adjacency_list(adjlist7)
        group8 = Group().from_adjacency_list(adjlist8)

        new_group = group7.add_implicit_atoms_from_atomtype()
        self.assertTrue(group8.is_isomorphic(new_group))

        # test multiple implicit atoms at a time
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
        group9 = Group().from_adjacency_list(adjlist9)
        group10 = Group().from_adjacency_list(adjlist10)

        new_group = group9.add_implicit_atoms_from_atomtype()
        self.assertTrue(group10.is_isomorphic(new_group))

    def test_classify_benzene_carbons(self):
        """
        Tests the method classifyingBenzeneCarbons
        """

        # This tests that we classify Cb atom types correctly
        adjlist1 = """
1  *1 Cb u0 {2,B}
2  *2 Cb u0 {1,B}
"""
        group1 = Group().from_adjacency_list(adjlist1)

        (cbAtomList, cbfAtomList, cbfAtomList1, cbfAtomList2, connectedCbfs) = group1.classify_benzene_carbons()
        self.assertEquals(len(cbfAtomList), 0)
        for atom in group1.atoms:
            self.assertIn(atom, cbAtomList)

        # This tests that we classify Cbf atomtypes correctly
        adjlist2 = """
1 *1 Cbf u0
"""
        group2 = Group().from_adjacency_list(adjlist2)
        (cbAtomList, cbfAtomList, cbfAtomList1, cbfAtomList2, connectedCbfs) = group2.classify_benzene_carbons()
        self.assertIn(group2.atoms[0], cbfAtomList)
        self.assertIn(group2.atoms[0], cbfAtomList1)

        # This tests that we can classify Cb atoms based on bonding and not just atomtype
        adjlist3 = """
1 *1 C u0 {2,B}
2 *2 C u0 {1,B} {3,B}
3 *3 C u0 {2,B}
"""
        group3 = Group().from_adjacency_list(adjlist3)
        (cbAtomList, cbfAtomList, cbfAtomList1, cbfAtomList2, connectedCbfs) = group3.classify_benzene_carbons()
        for atom in group3.atoms:
            self.assertIn(atom, cbAtomList)

        # This tests that we can classify Cbf1 atoms based on bonding and not just atomtype
        adjlist4 = """
1 *1 C u0 {2,B} {3,B} {4,B}
2 *2 C u0 {1,B}
3 *3 C u0 {1,B}
4 *4 C u0 {1,B}
"""
        group4 = Group().from_adjacency_list(adjlist4)
        (cbAtomList, cbfAtomList, cbfAtomList1, cbfAtomList2, connectedCbfs) = group4.classify_benzene_carbons()
        self.assertEquals(len(cbfAtomList), 1)
        self.assertEquals(len(cbfAtomList1), 1)

        # This tests that we can classify Cbf2 atoms. In the following partial group, we should have:
        # one Cbf2 atom, two Cbf1 atoms, and 5 Cb atoms
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
        group5 = Group().from_adjacency_list(adjlist5)
        (cbAtomList, cbfAtomList, cbfAtomList1, cbfAtomList2, connectedCbfs) = group5.classify_benzene_carbons()
        self.assertEquals(len(cbfAtomList1), 2)
        self.assertEquals(len(cbfAtomList2), 1)
        self.assertEquals(len(cbAtomList), 5)

        # Tests that we can classify connected Cbfs correctly. *1 should be connected to both *2 and *3
        atom1 = group5.get_labeled_atoms("*1")[0]
        atom2 = group5.get_labeled_atoms("*2")[0]
        atom3 = group5.get_labeled_atoms("*3")[0]
        self.assertIn(atom2, connectedCbfs[atom1])
        self.assertIn(atom3, connectedCbfs[atom1])
        self.assertIn(atom1, connectedCbfs[atom2])
        self.assertIn(atom1, connectedCbfs[atom3])

    def test_sort_by_connectivity(self):
        """
        Tests sort_by_connectivity method
        """

        # Basic test, we should get *1, *3 *2
        adjlist1 = """
1 *1 C u0 {3,B}
2 *2 C u0 {3,B}
3 *3 C u0 {1,B} {2,B}
"""
        group1 = Group().from_adjacency_list(adjlist1)
        ordered_atoms = group1.sort_by_connectivity(group1.atoms)
        self.assertEquals([x.label for x in ordered_atoms], ["*1", "*3", "*2"])

        # Check a detached case, we should get *1, *3, *4, *2, *5
        adjlist2 = """
1 *1 C u0 {3,B}
2 *2 C u0 {4,S} {5,B}
3 *3 C u0 {1,B} {4,B}
4 *4 C u0 {3,B} {2,S}
5 *5 C u0 {2,B}
"""
        group2 = Group().from_adjacency_list(adjlist2)
        ordered_atoms = group2.sort_by_connectivity(group2.atoms)
        self.assertEquals([x.label for x in ordered_atoms], ["*1", "*3", "*4", "*2", "*5"])

    def test_add_implicit_benzene(self):
        """
        Test the Group.add_implicit_benzene method
        """

        # tests it can make a benzene molecule
        adjlist1 = """
1  *1 Cb u0 {2,B}
2  *2 Cb u0 {1,B}
"""
        # tests it can make a bi-phenyl
        adjlist2 = """
1  *1 Cb u0 {2,S}
2  *2 Cb u0 {1,S}
"""
        # tests it can make a napthalene
        adjlist3 = """
1  *1 Cbf u0
"""

        # Test handling of Cbf2 atoms
        adjlist4 = """
1  *1 Cbf u0 p2 c0 {2,B}
2  *2 Cbf u0 p0 c0 {1,B} {3,B}
3  *3 Cbf u0 p0 c0 {2,B}
"""

        # test handling of heteroatoms and wildcards
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

        adjlist7 = """
1 *1 Cb u0 {4,B}
2 *2 Cb u0 {3,B}
3 *3 Cb u0 {4,B} {2,B}
4 *4 Cb u0 {1,B} {3,B}
"""

        benzene = """
1 C u0 {2,B} {6,B}
2 C u0 {1,B} {3,B}
3 C u0 {2,B} {4,B}
4 C u0 {3,B} {5,B}
5 C u0 {4,B} {6,B}
6 C u0 {5,B} {1,B}
"""

        biphenyl = """
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

        naphthalene = """
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
        answer6 = """
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

        group1 = Group().from_adjacency_list(adjlist1)
        group2 = Group().from_adjacency_list(adjlist2)
        group3 = Group().from_adjacency_list(adjlist3)
        group4 = Group().from_adjacency_list(adjlist4)
        group5 = Group().from_adjacency_list(adjlist5)
        group6 = Group().from_adjacency_list(adjlist6)
        group7 = Group().from_adjacency_list(adjlist7)

        benzene_group = Group().from_adjacency_list(benzene)
        biphenyl_group = Group().from_adjacency_list(biphenyl)
        naphthalene_group = Group().from_adjacency_list(naphthalene)
        phenanthrene_group = Group().from_adjacency_list(phenanthrene)
        answer5 = Group().from_adjacency_list(answer5)
        answer6 = Group().from_adjacency_list(answer6)

        group1 = group1.add_implicit_benzene()
        self.assertTrue(benzene_group.is_isomorphic(group1))
        group2 = group2.add_implicit_benzene()
        self.assertTrue(biphenyl_group.is_isomorphic(group2))
        group3 = group3.add_implicit_benzene()
        self.assertTrue(naphthalene_group.is_isomorphic(group3))
        group4 = group4.add_implicit_benzene()
        self.assertTrue(phenanthrene_group.is_isomorphic(group4))
        group5 = group5.add_implicit_benzene()
        self.assertTrue(answer5.is_isomorphic(group5))
        group6 = group6.add_implicit_benzene()
        self.assertTrue(answer6.is_isomorphic(group6))
        group7 = group7.add_implicit_benzene()
        self.assertTrue(benzene_group.is_isomorphic(group7))

    def test_pick_wildcards(self):
        """
        Test the Group.pickWildCards function
        """
        # The following tests are for picking optimal bond orders when there are bond wilcards
        # test that Cb/Cbf atoms with [D,B] chooses [B] instead of [D] bonds
        adjlist1 = """
    1 *1 R!H       u1 {2,[D,B]}
    2 *2 [Cbf,Cdd] u0 {1,[D,B]} {3,[D,B]}
    3 *3 [Cb,Cd]   u0 {2,[D,B]} {4,S}
    4 *4 R!H       u0 {3,S} {5,S}
    5 *5 H         u0 {4,S}
"""
        group1 = Group().from_adjacency_list(adjlist1)
        group1.pick_wildcards()
        atoms = group1.atoms
        self.assertTrue(atoms[0].bonds[atoms[1]].is_benzene())
        self.assertTrue(atoms[1].bonds[atoms[2]].is_benzene())

        adjlist2 = """
    1 *1 R!H       u1 {2,[S,D]} {4,[S,D]}
    2 *2 [CO,Cdd]  u0 {1,[S,D]} {3,[S,D]}
    3 *3 [O2d,Cd]  u0 {2,[S,D]}
    4 *4 [Cdd,Cd]  u0 {1,[S,D]}
"""
        group2 = Group().from_adjacency_list(adjlist2)
        group2.pick_wildcards()
        atoms = group2.atoms
        self.assertTrue(atoms[1].bonds[atoms[2]].is_double)
        self.assertTrue(atoms[0].bonds[atoms[3]].is_double)

    def test_make_sample_molecule(self):
        """
        Test the Group.make_sample_molecule method
        """

        def perform_samp_mole_comparison(_adjlist, _answer_smiles):
            """
            Creates a sample molecule from the adjlist and returns if it is isomorphic to a molecule created from
            the inputted smiles
            """
            group = Group().from_adjacency_list(_adjlist)
            result = group.make_sample_molecule()
            return result.is_isomorphic(Molecule().from_smiles(_answer_smiles))

        # tests adding implicit atoms
        adjlist = """
1  *1 Cd u0
"""
        answer_smiles = 'C=C'
        self.assertTrue(perform_samp_mole_comparison(adjlist, answer_smiles))

        # test creating implicit benzene atoms
        adjlist2 = """
1  *1 Cbf u0 {2,B}
2     Cbf u0 {1,B}
"""

        group2 = Group().from_adjacency_list(adjlist2)
        result2 = group2.make_sample_molecule()
        naphthalene_molecule = Molecule().from_smiles('C1=CC=C2C=CC=CC2=C1')
        resonance_list2 = naphthalene_molecule.generate_resonance_structures()
        self.assertTrue(any([result2.is_isomorphic(x) for x in resonance_list2]))

        # test the creation of a positively charged species
        adjlist = """
1  *1 N5sc u0
"""
        answer_smiles = '[NH4+]'
        self.assertTrue(perform_samp_mole_comparison(adjlist, answer_smiles))

        adjlist = """
1  *1 P5dc u0
"""
        answer_smiles = '[PH2+]=C'
        self.assertTrue(perform_samp_mole_comparison(adjlist, answer_smiles))

        # test the creation of a negatively charged species
        adjlist = """
1  *1 N1sc u0
"""
        answer_smiles = '[NH2-]'
        self.assertTrue(perform_samp_mole_comparison(adjlist, answer_smiles))

        adjlist = """
1  *1 P5sc u0
"""
        answer_smiles = '[PH6-]'
        self.assertTrue(perform_samp_mole_comparison(adjlist, answer_smiles))

        # test creation of charged species when some single bonds present
        adjlist = """
1 *2 [N5sc,N5dc] u0 {2,S} {3,S}
2 *3 R!H         u1 {1,S}
3 *4 H           u0 {1,S}
"""
        answer_smiles = '[NH3+][CH2]'
        self.assertTrue(perform_samp_mole_comparison(adjlist, answer_smiles))

    def test_is_benzene_explicit(self):
        """
        Test the Group.is_benzene_explicit method
        """
        adjlist1 = """
1  *1 Cb u0 {2,B}
2  *2 Cb u0 {1,B}
"""
        group1 = Group().from_adjacency_list(adjlist1)
        self.assertFalse(group1.is_benzene_explicit())

        benzene = """
1 C u0 {2,B} {6,B}
2 C u0 {1,B} {3,B}
3 C u0 {2,B} {4,B}
4 C u0 {3,B} {5,B}
5 C u0 {4,B} {6,B}
6 C u0 {5,B} {1,B}
"""
        benzene = Group().from_adjacency_list(benzene)
        self.assertTrue(benzene.is_benzene_explicit())

    def test_repr_png(self):
        """Test that a png representation can be created."""
        adjlist = """
1 *1 [C,Cd,Ct,CO,CS,Cb] u1 {2,[S,D,T,B]}
2 *2 [C,Cd,Ct,CO,CS,Cb] u0 {1,[S,D,T,B]} {3,[S,D,T,B]}
3 *3 [C,Cd,Ct,CO,CS,Cb] u0 {2,[S,D,T,B]} {4,[S,D,T,B]}
4 *4 [C,Cd,Ct,CO,CS,Cb] u0 {3,[S,D,T,B]}
"""
        group = Group().from_adjacency_list(adjlist)
        result = group._repr_png_()
        self.assertIsNotNone(result)

    def test_draw_group(self):
        """Test that the draw method returns the expected pydot graph."""
        adjlist = """
1 *1 [C,Cd,Ct,CO,CS,Cb] u1 {2,[S,D,T,B]}
2 *2 [C,Cd,Ct,CO,CS,Cb] u0 {1,[S,D,T,B]} {3,[S,D,T,B]}
3 *3 [C,Cd,Ct,CO,CS,Cb] u0 {2,[S,D,T,B]} {4,[S,D,T,B]}
4 *4 [C,Cd,Ct,CO,CS,Cb] u0 {3,[S,D,T,B]}
"""
        # Use of tabs in the expected string is intentional
        expected = rb"""
graph G {
	graph [dpi=52];
	node [label="\N"];
	1	 [fontname=Helvetica,
		fontsize=16,
		label="1 *1 C,Cd,Ct,CO,CS,Cb"];
	2	 [fontname=Helvetica,
		fontsize=16,
		label="2 *2 C,Cd,Ct,CO,CS,Cb"];
	1 -- 2	 [fontname=Helvetica,
		fontsize=16,
		label="S,D,T,B"];
	3	 [fontname=Helvetica,
		fontsize=16,
		label="3 *3 C,Cd,Ct,CO,CS,Cb"];
	2 -- 3	 [fontname=Helvetica,
		fontsize=16,
		label="S,D,T,B"];
	4	 [fontname=Helvetica,
		fontsize=16,
		label="4 *4 C,Cd,Ct,CO,CS,Cb"];
	3 -- 4	 [fontname=Helvetica,
		fontsize=16,
		label="S,D,T,B"];
}
        """
        group = Group().from_adjacency_list(adjlist)
        result = group.draw('canon')
        self.assertEqual(b''.join(result.split()), b''.join(expected.split()))

    def test_merge_groups(self):
        """
        Test the merge_groups() function
        """
        # basic test of merging a backbone and end group
        backbone1 = Group().from_adjacency_list("""
1 *1 R!H u1 {2,S}
2 *4 R!H u0 {1,S} {3,S}
3 *6 R!H u0 {2,S} {4,S}
4 *5 R!H u0 {3,S} {5,S}
5 *2 R!H u0 {4,S} {6,S}
6 *3 H   u0 {5,S}
""")

        end1 = Group().from_adjacency_list("""
1 *2 Cs u0 {2,S} {3,S}
2 *3 H  u0 {1,S}
3    S  u0 {1,S}
""")
        desired_merge1 = Group().from_adjacency_list("""
1 *1 R!H u1 {2,S}
2 *4 R!H u0 {1,S} {3,S}
3 *6 R!H u0 {2,S} {4,S}
4 *5 R!H u0 {3,S} {5,S}
5 *2 Cs  u0 {4,S} {6,S} {7,S}
6 *3 H   u0 {5,S}
7    S   u0 {5,S}
""")

        merged_group = backbone1.merge_groups(end1)
        self.assertTrue(merged_group.is_identical(desired_merge1))

        # test it works when there is a cyclical structure to the backbone

        backbone2 = Group().from_adjacency_list("""
1 *1 R!H u1 {2,S} {4,S}
2 *4 R!H u0 {1,S} {3,S}
3 *2 R!H u0 {2,S} {4,S}
4 *3 R!H u0 {3,S} {1,S}
""")

        end2 = Group().from_adjacency_list("""
1 *2 O2s u0 {2,S}
2 *3 Cs  u0 {1,S}
""")
        desired_merge2 = Group().from_adjacency_list("""
1 *1 R!H u1 {2,S} {4,S}
2 *4 R!H u0 {1,S} {3,S}
3 *2 O2s u0 {2,S} {4,S}
4 *3 Cs  u0 {3,S} {1,S}
""")
        merged_group = backbone2.merge_groups(end2)
        self.assertTrue(merged_group.is_identical(desired_merge2))

    def test_get_element_count(self):
        """Test that we can count elements properly."""
        group = Group().from_adjacency_list("""
1 R!H u0 {2,S}
2 [Cs,Cd,Ct,Cb] u0 {1,S} {3,S}
3 [Cs,Cd,Ct,Cb,O2s,S2s] u0 {2,S} {4,S}
4 N3s u0 {3,S} {5,S}
5 P1s u0 {4,S} 
""")
        expected = {'C': 1, 'N': 1, 'P': 1}
        result = group.get_element_count()
        self.assertEqual(expected, result)


################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
