#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2023 Prof. William H. Green (whgreen@mit.edu),           #
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


from rmgpy.exceptions import InchiException
from rmgpy.molecule.element import get_element, element_list
from rmgpy.molecule.group import Group, ActionError
from rmgpy.molecule.molecule import Atom, Bond, Molecule
import pytest


class TestAtom:
    """
    Contains unit tests of the Atom class.
    """

    def setup_class(self):
        """
        A method called before each unit test in this class.
        """
        self.atom = Atom(
            element=get_element("C"),
            radical_electrons=1,
            charge=0,
            label="*1",
            lone_pairs=0,
        )

        self.atom1 = Atom(element=get_element("C"), radical_electrons=0, lone_pairs=0)
        self.atom2 = Atom(element=get_element("C"), radical_electrons=0, lone_pairs=0)
        self.atom3 = Atom(element=get_element("C"), radical_electrons=1, lone_pairs=0)
        self.atom4 = Atom(element=get_element("H"), radical_electrons=1, lone_pairs=0)

    def test_mass(self):
        """
        Test the Atom.mass property.
        """
        assert self.atom.mass == self.atom.element.mass

    def test_number(self):
        """
        Test the Atom.number property.
        """
        assert self.atom.number == self.atom.element.number

    def test_symbol(self):
        """
        Test the Atom.symbol property.
        """
        assert self.atom.symbol == self.atom.element.symbol

    def test_equality(self):
        """Test that we can perform equality comparison with Atom objects"""
        assert self.atom1 == self.atom1
        assert self.atom1 != self.atom2
        assert self.atom1 != self.atom3
        assert self.atom1 != self.atom4

    def test_less_than(self):
        """Test that we can perform less than comparison with Atom objects"""
        assert not (self.atom1 < self.atom2)  # Because the sorting keys should be identical
        assert self.atom2 < self.atom3
        assert self.atom4 < self.atom1

    def test_greater_than(self):
        """Test that we can perform greater than comparison with Atom objects"""
        assert not (self.atom2 > self.atom1)  # Because the sorting keys should be identical
        assert self.atom3 > self.atom1
        assert self.atom1 > self.atom4

    def test_hash(self):
        """Test behavior of Atom hashing using dictionaries and sets"""
        # Test dictionary behavior
        assert len(dict.fromkeys([self.atom1, self.atom2, self.atom3, self.atom4])) == 4

        # Test set behavior
        assert len({self.atom1, self.atom2, self.atom3, self.atom4}) == 4

    def test_is_hydrogen(self):
        """
        Test the Atom.is_hydrogen() method.
        """
        for element in element_list:
            atom = Atom(element=element, radical_electrons=1, charge=0, label="*1", lone_pairs=0)
            if element.symbol == "H":
                assert atom.is_hydrogen()
            else:
                assert not atom.is_hydrogen()

    def test_is_non_hydrogen(self):
        """
        Test the Atom.is_non_hydrogen() method.
        """
        for element in element_list:
            atom = Atom(element=element, radical_electrons=1, charge=0, label="*1", lone_pairs=0)
            if element.symbol == "H":
                assert not atom.is_non_hydrogen()
            else:
                assert atom.is_non_hydrogen(), "Atom {0!r} isn't reporting is_non_hydrogen()".format(atom)

    def test_is_halogen(self):
        """
        Test the Atom.is_halogen() method.
        """
        for element in element_list:
            atom = Atom(element=element, radical_electrons=1, charge=0, label="*1", lone_pairs=3)
            if element.symbol in ["F", "Cl", "Br", "I"]:
                assert atom.is_halogen()
            else:
                assert not atom.is_halogen(), "Atom {0!r} is reporting is_halogen(), but it shouldn't be".format(atom)

    def test_is_carbon(self):
        """
        Test the Atom.is_carbon() method.
        """
        for element in element_list:
            atom = Atom(element=element, radical_electrons=1, charge=0, label="*1", lone_pairs=0)
            if element.symbol == "C":
                assert atom.is_carbon()
            else:
                assert not atom.is_carbon()

    def test_is_silicon(self):
        """
        Test the Atom.is_silicon() method.
        """
        for element in element_list:
            atom = Atom(element=element, radical_electrons=1, charge=0, label="*1", lone_pairs=0)
            if element.symbol == "Si":
                assert atom.is_silicon()
            else:
                assert not atom.is_silicon()

    def test_is_oxygen(self):
        """
        Test the Atom.is_oxygen() method.
        """
        for element in element_list:
            atom = Atom(element=element, radical_electrons=2, charge=0, label="*1", lone_pairs=2)
            if element.symbol == "O":
                assert atom.is_oxygen()
            else:
                assert not atom.is_oxygen()

    def test_is_nitrogen(self):
        """
        Test the Atom.is_nitrogen() method.
        """
        for element in element_list:
            atom = Atom(element=element, radical_electrons=1, charge=0, label="*1", lone_pairs=1)
            if element.symbol == "N":
                assert atom.is_nitrogen()
            else:
                assert not atom.is_nitrogen()

    def test_is_phosphorus(self):
        """
        Test the Atom.is_phosphorus() method.
        """
        for element in element_list:
            atom = Atom(element=element, radical_electrons=1, charge=0, label="*1", lone_pairs=1)
            if element.symbol == "P":
                assert atom.is_phosphorus()
            else:
                assert not atom.is_phosphorus()

    def test_is_sulfur(self):
        """
        Test the Atom.is_sulfur() method.
        """
        for element in element_list:
            atom = Atom(element=element, radical_electrons=1, charge=0, label="*1", lone_pairs=2)
            if element.symbol == "S":
                assert atom.is_sulfur()
            else:
                assert not atom.is_sulfur()

    def test_is_fluorine(self):
        """
        Test the Atom.is_fluorine() method.
        """
        for element in element_list:
            atom = Atom(element=element, radical_electrons=1, charge=0, label="*1", lone_pairs=3)
            if element.symbol == "F":
                assert atom.is_fluorine()
            else:
                assert not atom.is_fluorine()

    def test_is_chlorine(self):
        """
        Test the Atom.is_chlorine() method.
        """
        for element in element_list:
            atom = Atom(element=element, radical_electrons=1, charge=0, label="*1", lone_pairs=3)
            if element.symbol == "Cl":
                assert atom.is_chlorine()
            else:
                assert not atom.is_chlorine()

    def test_is_bromine(self):
        """
        Test the Atom.is_bromine() method.
        """
        for element in element_list:
            atom = Atom(element=element, radical_electrons=1, charge=0, label="*1", lone_pairs=3)
            if element.symbol == "Br":
                assert atom.is_bromine()
            else:
                assert not atom.is_bromine()

    def test_is_iodine(self):
        """
        Test the Atom.is_iodine() method.
        """
        for element in element_list:
            atom = Atom(element=element, radical_electrons=1, charge=0, label="*1", lone_pairs=3)
            if element.symbol == "I":
                assert atom.is_iodine()
            else:
                assert not atom.is_iodine()

    def test_is_nos(self):
        """
        Test the Atom.is_nos() method.
        """
        for element in element_list:
            atom = Atom(element=element, radical_electrons=1, charge=0, label="*1", lone_pairs=2)
            if element.symbol in ["N", "O", "S"]:
                assert atom.is_nos()
            else:
                assert not atom.is_nos()

    def test_is_surface_site(self):
        """
        Test the Atom.is_surface_site() method.
        """
        for element in element_list:
            atom = Atom(element=element, radical_electrons=0, charge=0, label="*1", lone_pairs=0)
            if element.symbol == "X":
                assert atom.is_surface_site()
            else:
                assert not atom.is_surface_site()

    def test_is_bonded_to_surface(self):
        """
        Test the Atom.is_bonded_to_surface_site() method.
        """

        adsorbate = Molecule(smiles="*=O")  # X=O
        for atom in adsorbate.atoms:
            if atom.is_surface_site():
                assert not atom.is_bonded_to_surface()
            else:
                assert atom.is_bonded_to_surface()

    def test_is_bonded_to_halogen(self):
        """
        Test the Atom.is_bonded_to_halogen() method.
        """

        cf4 = Molecule(smiles="FC(F)(F)F")  # CF4
        ch4 = Molecule(smiles="C")  # CH4

        for atom in cf4.atoms:
            if atom.is_halogen():
                assert not atom.is_bonded_to_halogen()
            else:
                assert atom.is_bonded_to_halogen()

        for atom in ch4.atoms:
            assert not atom.is_bonded_to_halogen()

    def test_increment_radical(self):
        """
        Test the Atom.increment_radical() method.
        """
        radical_electrons = self.atom.radical_electrons
        self.atom.increment_radical()
        assert self.atom.radical_electrons == radical_electrons + 1

    def test_decrement_radical(self):
        """
        Test the Atom.decrement_radical() method.
        """
        radical_electrons = self.atom.radical_electrons
        self.atom.decrement_radical()
        assert self.atom.radical_electrons == radical_electrons - 1

    def test_apply_action_break_bond(self):
        """
        Test the Atom.apply_action() method for a BREAK_BOND action.
        """
        action = ["BREAK_BOND", "*1", 1, "*2"]
        for element in element_list:
            atom0 = Atom(element=element, radical_electrons=1, charge=0, label="*1", lone_pairs=0)
            atom = atom0.copy()
            atom.apply_action(action)
            assert atom0.element == atom.element
            assert atom0.radical_electrons == atom.radical_electrons
            assert atom0.charge == atom.charge
            assert atom0.label == atom.label

    def test_apply_action_form_bond(self):
        """
        Test the Atom.apply_action() method for a FORM_BOND action.
        """
        action = ["FORM_BOND", "*1", 1, "*2"]
        for element in element_list:
            atom0 = Atom(element=element, radical_electrons=1, charge=0, label="*1", lone_pairs=0)
            atom = atom0.copy()
            atom.apply_action(action)
            assert atom0.element == atom.element
            assert atom0.radical_electrons == atom.radical_electrons
            assert atom0.charge == atom.charge
            assert atom0.label == atom.label

    def test_apply_action_increment_bond(self):
        """
        Test the Atom.apply_action() method for a CHANGE_BOND action.
        """
        action = ["CHANGE_BOND", "*1", 1, "*2"]
        for element in element_list:
            atom0 = Atom(element=element, radical_electrons=1, charge=0, label="*1", lone_pairs=0)
            atom = atom0.copy()
            atom.apply_action(action)
            assert atom0.element == atom.element
            assert atom0.radical_electrons == atom.radical_electrons
            assert atom0.charge == atom.charge
            assert atom0.label == atom.label

    def test_apply_action_decrement_bond(self):
        """
        Test the Atom.apply_action() method for a CHANGE_BOND action.
        """
        action = ["CHANGE_BOND", "*1", -1, "*2"]
        for element in element_list:
            atom0 = Atom(element=element, radical_electrons=1, charge=0, label="*1", lone_pairs=0)
            atom = atom0.copy()
            atom.apply_action(action)
            assert atom0.element == atom.element
            assert atom0.radical_electrons == atom.radical_electrons
            assert atom0.charge == atom.charge
            assert atom0.label == atom.label

    def test_apply_action_gain_radical(self):
        """
        Test the Atom.apply_action() method for a GAIN_RADICAL action.
        """
        action = ["GAIN_RADICAL", "*1", 1]
        for element in element_list:
            atom0 = Atom(element=element, radical_electrons=1, charge=0, label="*1", lone_pairs=0)
            atom = atom0.copy()
            atom.apply_action(action)
            assert atom0.element == atom.element
            assert atom0.radical_electrons == atom.radical_electrons - 1
            assert atom0.charge == atom.charge
            assert atom0.label == atom.label

    def test_apply_action_lose_radical(self):
        """
        Test the Atom.apply_action() method for a LOSE_RADICAL action.
        """
        action = ["LOSE_RADICAL", "*1", 1]
        for element in element_list:
            atom0 = Atom(element=element, radical_electrons=1, charge=0, label="*1", lone_pairs=0)
            atom = atom0.copy()
            atom.apply_action(action)
            assert atom0.element == atom.element
            assert atom0.radical_electrons == atom.radical_electrons + 1
            assert atom0.charge == atom.charge
            assert atom0.label == atom.label

    def test_equivalent(self):
        """
        Test the Atom.equivalent() method.
        """
        for index1, element1 in enumerate(element_list[0:10]):
            for index2, element2 in enumerate(element_list[0:10]):
                atom1 = Atom(
                    element=element1,
                    radical_electrons=1,
                    charge=0,
                    label="*1",
                    lone_pairs=0,
                )
                atom2 = Atom(
                    element=element2,
                    radical_electrons=1,
                    charge=0,
                    label="*1",
                    lone_pairs=0,
                )
                if index1 == index2:
                    assert atom1.equivalent(atom2)
                    assert atom2.equivalent(atom1)
                else:
                    assert not atom1.equivalent(atom2)
                    assert not atom2.equivalent(atom1)

    def test_is_specific_case_of(self):
        """
        Test the Atom.is_specific_case_of() method.
        """
        for index1, element1 in enumerate(element_list[0:10]):
            for index2, element2 in enumerate(element_list[0:10]):
                atom1 = Atom(
                    element=element1,
                    radical_electrons=1,
                    charge=0,
                    label="*1",
                    lone_pairs=0,
                )
                atom2 = Atom(
                    element=element2,
                    radical_electrons=1,
                    charge=0,
                    label="*1",
                    lone_pairs=0,
                )
                if index1 == index2:
                    assert atom1.is_specific_case_of(atom2)
                else:
                    assert not atom1.is_specific_case_of(atom2)

    def test_copy(self):
        """
        Test the Atom.copy() method.
        """
        atom = self.atom.copy()
        assert self.atom.element.symbol == atom.element.symbol
        assert self.atom.atomtype == atom.atomtype
        assert self.atom.radical_electrons == atom.radical_electrons
        assert self.atom.charge == atom.charge
        assert self.atom.label == atom.label

    def test_pickle(self):
        """
        Test that a Atom object can be successfully pickled and
        unpickled with no loss of information.
        """
        import pickle

        atom = pickle.loads(pickle.dumps(self.atom))
        assert self.atom.element.symbol == atom.element.symbol
        assert self.atom.atomtype == atom.atomtype
        assert self.atom.radical_electrons == atom.radical_electrons
        assert self.atom.charge == atom.charge
        assert self.atom.label == atom.label

    def test_isotope_equivalent(self):
        """
        Test the Atom.equivalent() method for non-normal isotopes
        """

        atom1 = Atom(element=get_element("H"))
        atom2 = Atom(element=get_element("H", 2))
        atom3 = Atom(element=get_element("H"))

        assert not atom1.equivalent(atom2)
        assert atom1.equivalent(atom3)

    def test_get_bond_orders_for_atom(self):
        """
        Test Atom.get_total_bond_order for all carbons in naphthalene
        """

        m = Molecule().from_smiles("C12C(C=CC=C1)=CC=CC=2")
        isomers = m.generate_resonance_structures()
        for isomer in isomers:
            for atom in isomer.atoms:
                if atom.symbol == "C":
                    assert atom.get_total_bond_order() == 4.0


class TestBond:
    """
    Contains unit tests of the Bond class.
    """

    def setup_method(self):
        """
        A method called before each unit test in this class.
        """
        self.bond = Bond(atom1=None, atom2=None, order=2)
        self.orderList = [1, 2, 3, 4, 1.5, 0.30000000000000004]

        self.bond1 = Bond(atom1=None, atom2=None, order=1)
        self.bond2 = Bond(atom1=None, atom2=None, order=1)
        self.bond3 = Bond(atom1=None, atom2=None, order=2)
        self.bond4 = Bond(atom1=None, atom2=None, order=3)

    def test_equality(self):
        """Test that we can perform equality comparison with Bond objects"""
        assert self.bond1 == self.bond1
        assert self.bond1 != self.bond2
        assert self.bond1 != self.bond3
        assert self.bond1 != self.bond4

    def test_less_than(self):
        """Test that we can perform less than comparison with Bond objects"""
        assert not (self.bond1 < self.bond2)  # Because the sorting keys should be identical
        assert self.bond2 < self.bond3
        assert self.bond3 < self.bond4

    def test_greater_than(self):
        """Test that we can perform greater than comparison with Bond objects"""
        assert not (self.bond2 > self.bond1)  # Because the sorting keys should be identical
        assert self.bond3 > self.bond1
        assert self.bond4 > self.bond1

    def test_hash(self):
        """Test behavior of Bond hashing using dictionaries and sets"""
        # Test dictionary behavior
        assert len(dict.fromkeys([self.bond1, self.bond2, self.bond3, self.bond4])) == 4

        # Test set behavior
        assert len({self.bond1, self.bond2, self.bond3, self.bond4}) == 4

    def test_get_order_str(self):
        """
        test the Bond.get_order_str() method
        """

        assert self.bond.get_order_str() == "D"

    def test_set_order_str(self):
        """
        test the Bond.set_order_str() method
        """

        self.bond.set_order_str("B")
        assert self.bond.order == 1.5

    def test_get_order_num(self):
        """
        test the Bond.get_order_num() method
        """
        assert self.bond.get_order_num() == 2

    def test_set_order_num(self):
        """
        test the Bond.set_order_num() method
        """

        self.bond.set_order_num(3)
        assert self.bond.get_order_str() == "T"

    def test_is_order(self):
        """
        Test the Bond.is_order() method.
        """
        for order in self.orderList:
            bond = Bond(None, None, order=order)
            assert bond.is_order(round(order, 2))

    def test_is_single(self):
        """
        Test the Bond.is_single() method.
        """
        for order in self.orderList:
            bond = Bond(None, None, order=order)
            if order == 1:
                assert bond.is_single()
            else:
                assert not bond.is_single()

    def test_is_single_can_take_floating_point_addition(self):
        """
        Test the Bond.is_single() method with taking floating point addition
        roundoff errors
        """
        new_order = 0.1 + 0.3 * 3
        assert new_order != 1

        self.bond.set_order_num(new_order)
        assert self.bond.is_single()

    def test_is_double(self):
        """
        Test the Bond.is_double() method.
        """
        for order in self.orderList:
            bond = Bond(None, None, order=order)
            if order == 2:
                assert bond.is_double()
            else:
                assert not bond.is_double()

    def test_is_triple(self):
        """
        Test the Bond.is_triple() method.
        """
        for order in self.orderList:
            bond = Bond(None, None, order=order)
            if order == 3:
                assert bond.is_triple()
            else:
                assert not bond.is_triple()

    def test_is_benzene(self):
        """
        Test the Bond.is_benzene() method.
        """
        for order in self.orderList:
            bond = Bond(None, None, order=order)
            if order == 1.5:
                assert bond.is_benzene()
            else:
                assert not bond.is_benzene()

    def test_is_quadruple(self):
        """
        Test the Bond.is_quadruple() method.
        """
        for order in self.orderList:
            bond = Bond(None, None, order=order)
            if order == 4:
                assert bond.is_quadruple()
            else:
                assert not bond.is_quadruple()

    def test_increment_order(self):
        """
        Test the Bond.increment_order() method.
        """
        for order in self.orderList:
            bond = Bond(None, None, order=order)
            try:
                bond.increment_order()
                if order == 1:
                    assert bond.is_double()
                elif order == 2:
                    assert bond.is_triple()
                elif order == 3:
                    assert bond.is_quadruple()
            except ActionError:
                assert order >= 4  # or benzene??

    def test_decrement_order(self):
        """
        Test the Bond.decrement_order() method.
        """
        for order in self.orderList:
            bond = Bond(None, None, order=order)
            try:
                bond.decrement_order()
                if order == 2:
                    assert bond.is_single()
                elif order == 3:
                    assert bond.is_double()
                elif order == "Q":
                    assert bond.is_triple()
            except ActionError:
                assert order < 1

    def test_apply_action_break_bond(self):
        """
        Test the Bond.apply_action() method for a BREAK_BOND action.
        """
        action = ["BREAK_BOND", "*1", 1, "*2"]
        for order0 in self.orderList:
            bond0 = Bond(None, None, order=order0)
            bond = bond0.copy()
            try:
                bond.apply_action(action)
                assert False, "Bond.apply_action() unexpectedly processed a BREAK_BOND action " "with order {0}.".format(order0)
            except ActionError:
                pass

    def test_apply_action_form_bond(self):
        """
        Test the Bond.apply_action() method for a FORM_BOND action.
        """
        action = ["FORM_BOND", "*1", 1, "*2"]
        for order0 in self.orderList:
            bond0 = Bond(None, None, order=order0)
            bond = bond0.copy()
            try:
                bond.apply_action(action)
                assert False, "Bond.apply_action() unexpectedly processed a FORM_BOND action " "with order {0}.".format(order0)
            except ActionError:
                pass

    def test_apply_action_increment_bond(self):
        """
        Test the Bond.apply_action() method for a CHANGE_BOND action.
        """
        action = ["CHANGE_BOND", "*1", 1, "*2"]
        for order0 in self.orderList:
            bond0 = Bond(None, None, order=order0)
            bond = bond0.copy()
            try:
                bond.apply_action(action)
            except ActionError:
                assert 4 <= order0, "Test failed with order {0}".format(order0)

    def test_apply_action_decrement_bond(self):
        """
        Test the Bond.apply_action() method for a CHANGE_BOND action.
        """
        action = ["CHANGE_BOND", "*1", -1, "*2"]
        for order0 in self.orderList:
            bond0 = Bond(None, None, order=order0)
            bond = bond0.copy()
            try:
                bond.apply_action(action)
            except ActionError:
                assert order0 < 1, "Test failed with order {0}".format(order0)

    def test_apply_action_gain_radical(self):
        """
        Test the Bond.apply_action() method for a GAIN_RADICAL action.
        """
        action = ["GAIN_RADICAL", "*1", 1]
        for order0 in self.orderList:
            bond0 = Bond(None, None, order=order0)
            bond = bond0.copy()
            try:
                bond.apply_action(action)
                assert False, "Bond.apply_action() unexpectedly processed a GAIN_RADICAL action " "with order {0}.".format(order0)
            except ActionError:
                pass

    def test_apply_action_lose_radical(self):
        """
        Test the Bond.apply_action() method for a LOSE_RADICAL action.
        """
        action = ["LOSE_RADICAL", "*1", 1]
        for order0 in self.orderList:
            bond0 = Bond(None, None, order=order0)
            bond = bond0.copy()
            try:
                bond.apply_action(action)
                assert False, "Bond.apply_action() unexpectedly processed a LOSE_RADICAL action " "with order {0}.".format(order0)
            except ActionError:
                pass

    def test_equivalent(self):
        """
        Test the GroupBond.equivalent() method.
        """
        for order1 in self.orderList:
            for order2 in self.orderList:
                bond1 = Bond(None, None, order=order1)
                bond2 = Bond(None, None, order=order2)
                if order1 == order2:
                    assert bond1.equivalent(bond2)
                    assert bond2.equivalent(bond1)
                else:
                    assert not bond1.equivalent(bond2)
                    assert not bond2.equivalent(bond1)

    def test_is_specific_case_of(self):
        """
        Test the Bond.is_specific_case_of() method.
        """
        for order1 in self.orderList:
            for order2 in self.orderList:
                bond1 = Bond(None, None, order=order1)
                bond2 = Bond(None, None, order=order2)
                if order1 == order2:
                    assert bond1.is_specific_case_of(bond2)
                else:
                    assert not bond1.is_specific_case_of(bond2)

    def test_copy(self):
        """
        Test the Bond.copy() method.
        """
        bond = self.bond.copy()
        assert self.bond.order == bond.order

    def test_pickle(self):
        """
        Test that a Bond object can be successfully pickled and
        unpickled with no loss of information.
        """
        import pickle

        bond = pickle.loads(pickle.dumps(self.bond))
        assert self.bond.order == bond.order

    def test_update_lone_pairs(self):
        """
        Test that update_lone_pairs works as expected
        """
        mol_n1sc_n5t = Molecule().from_adjacency_list(
            """
            1 N u0 p0 c+1 {2,T} {4,S}
            2 N u0 p0 c+1 {1,T} {3,S}
            3 N u0 p3 c-2 {2,S}
            4 H u0 p0 c0 {1,S}"""
        )
        mol_n1s = Molecule().from_adjacency_list(
            """
            1 N u0 p2 c0 {2,S}
            2 H u0 p0 c0 {1,S}"""
        )
        mol_n3s = Molecule().from_adjacency_list(
            """
            multiplicity 3
            1 N u2 p1 c0 {2,S}
            2 H u0 p0 c0 {1,S}"""
        )
        mol_n3b = Molecule().from_adjacency_list(
            """
            1  N u0 p1 c0 {2,D} {6,S}
            2  C u0 p0 c0 {1,D} {3,S} {7,S}
            3  C u0 p0 c0 {2,S} {4,D} {8,S}
            4  C u0 p0 c0 {3,D} {5,S} {9,S}
            5  C u0 p0 c0 {4,S} {6,D} {10,S}
            6  C u0 p0 c0 {1,S} {5,D} {11,S}
            7  H u0 p0 c0 {2,S}
            8  H u0 p0 c0 {3,S}
            9  H u0 p0 c0 {4,S}
            10 H u0 p0 c0 {5,S}
            11 H u0 p0 c0 {6,S}"""
        )
        mol_n5s = Molecule().from_adjacency_list(
            """
            multiplicity 2
            1 N u1 p0 c+1 {2,S} {3,S} {4,S}
            2 H u0 p0 c0 {1,S}
            3 H u0 p0 c0 {1,S}
            4 O u0 p3 c-1 {1,S}"""
        )
        mol_n5d = Molecule().from_adjacency_list(
            """
            1 N u0 p0 c+1 {2,D} {3,S} {4,S}
            2 O u0 p2 c0 {1,D}
            3 O u0 p2 c0 {1,S} {5,S}
            4 O u0 p3 c-1 {1,S}
            5 H u0 p0 c0 {3,S}"""
        )
        mol_n5dd = Molecule().from_adjacency_list(
            """
            1 N u0 p2 c-1 {2,D}
            2 N u0 p0 c+1 {1,D} {3,D}
            3 O u0 p2 c0 {2,D}"""
        )
        mol_ch2_s = Molecule().from_adjacency_list(
            """
            1 C u0 p1 c0 {2,S} {3,S}
            2 H u0 p0 c0 {1,S}
            3 H u0 p0 c0 {1,S}"""
        )
        mol_carbonyl = Molecule().from_adjacency_list(
            """
            1 O u0 p2 c0 {2,D}
            2 C u0 p0 c0 {1,D} {3,S} {4,S}
            3 H u0 p0 c0 {2,S}
            4 H u0 p0 c0 {2,S}"""
        )

        mol_n1sc_n5t.update_lone_pairs()
        mol_n1s.update_lone_pairs()
        mol_n3s.update_lone_pairs()
        mol_n3b.update_lone_pairs()
        mol_n5s.update_lone_pairs()
        mol_n5d.update_lone_pairs()
        mol_n5dd.update_lone_pairs()
        mol_ch2_s.update_lone_pairs()
        mol_carbonyl.update_lone_pairs()

        assert mol_n1sc_n5t.atoms[0].lone_pairs == 0
        assert mol_n1sc_n5t.atoms[2].lone_pairs == 3
        assert mol_n1s.atoms[0].lone_pairs == 2
        assert mol_n3s.atoms[0].lone_pairs == 1
        assert mol_n3b.atoms[0].lone_pairs == 1
        assert mol_n5s.atoms[0].lone_pairs == 0
        assert mol_n5s.atoms[3].lone_pairs == 3
        assert mol_n5d.atoms[0].lone_pairs == 0
        assert mol_n5d.atoms[1].lone_pairs == 2
        assert mol_n5d.atoms[2].lone_pairs == 2
        assert mol_n5d.atoms[3].lone_pairs == 3
        assert mol_n5dd.atoms[0].lone_pairs == 2
        assert mol_n5dd.atoms[1].lone_pairs == 0
        assert mol_n5dd.atoms[2].lone_pairs == 2
        assert mol_ch2_s.atoms[0].lone_pairs == 1
        assert mol_carbonyl.atoms[0].lone_pairs == 2
        assert mol_carbonyl.atoms[1].lone_pairs == 0

    def test_get_bond_string(self):
        """Test that bond objects can return a bond string"""
        bond = Bond(
            atom1=Atom(element=get_element(1)),
            atom2=Atom(element=get_element(6)),
            order=1,
        )
        assert bond.get_bond_string() == "C-H"


class TestMolecule:
    """
    Contains unit tests of the Molecule class.
    """

    def setup_method(self):
        self.adjlist_1 = """
1 *1 C u1 p0 c0  {2,S} {3,S} {4,S}
2    H u0 p0 c0  {1,S}
3    H u0 p0 c0  {1,S}
4 *2 N u0 p0 c+1 {1,S} {5,S} {6,D}
5    O u0 p3 c-1 {4,S}
6    O u0 p2 c0  {4,D}
            """
        self.molecule = [Molecule().from_adjacency_list(self.adjlist_1)]

        self.adjlist_2 = """
1 *1 C u1 p0 {2,S}
2 *2 N u0 p0 c+1 {1,S} {3,S} {4,D}
3    O u0 p3 c-1 {2,S}
4    O u0 p2 {2,D}
            """
        self.molecule.append(Molecule().from_adjacency_list(self.adjlist_2, saturate_h=True))

        self.mHBonds = Molecule().from_smiles("C(NC=O)OO")

        self.mol1 = Molecule(smiles="C")
        self.mol2 = Molecule(smiles="C")
        self.mol3 = Molecule(smiles="CC")

    def test_equality(self):
        """Test that we can perform equality comparison with Molecule objects"""
        assert self.mol1 == self.mol1
        assert self.mol1 == self.mol2
        assert self.mol1 != self.mol3

    def test_less_than(self):
        """Test that we can perform less than comparison with Molecule objects"""
        assert not (self.mol1 < self.mol2)  # Because the sorting keys should be identical
        assert self.mol1 < self.mol3

    def test_greater_than(self):
        """Test that we can perform greater than comparison with Molecule objects"""
        assert not (self.mol2 > self.mol1)  # Because the sorting keys should be identical
        assert self.mol3 > self.mol1

    def test_hash(self):
        """Test behavior of Molecule hashing using dictionaries and sets"""
        # Test dictionary behavior
        assert len(dict.fromkeys([self.mol1, self.mol2, self.mol3])) == 2

        # Test set behavior
        assert len({self.mol1, self.mol2, self.mol3}) == 2

    def test_clear_labeled_atoms(self):
        """
        Test the Molecule.clear_labeled_atoms() method.
        """
        self.molecule[0].clear_labeled_atoms()
        for atom in self.molecule[0].atoms:
            assert atom.label == ""

    def test_contains_labeled_atom(self):
        """
        Test the Molecule.contains_labeled_atom() method.
        """
        for atom in self.molecule[0].atoms:
            if atom.label != "":
                assert self.molecule[0].contains_labeled_atom(atom.label)
        assert not self.molecule[0].contains_labeled_atom("*3")
        assert not self.molecule[0].contains_labeled_atom("*4")
        assert not self.molecule[0].contains_labeled_atom("*5")
        assert not self.molecule[0].contains_labeled_atom("*6")

    def test_get_labeled_atom(self):
        """
        Test the Molecule.get_labeled_atoms() method.
        """
        for atom in self.molecule[0].atoms:
            if atom.label != "":
                assert atom == self.molecule[0].get_labeled_atoms(atom.label)[0]
        try:
            self.molecule[0].get_labeled_atoms("*3")[0]
            assert False, "Unexpected successful return from Molecule.get_labeled_atoms() with invalid atom label."
        except ValueError:
            pass

    def test_get_labeled_atoms(self):
        """
        Test the Molecule.get_all_labeled_atoms() method.
        """
        labeled = self.molecule[0].get_all_labeled_atoms()
        for atom in self.molecule[0].atoms:
            if atom.label != "":
                assert atom.label in labeled
                assert atom in list(labeled.values())
            else:
                assert not (atom.label in labeled)
                assert not (atom in list(labeled.values()))

        multiple_label_molecule = Molecule().from_adjacency_list(
            """
1 * C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2 * C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3 * C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4 * C u0 p0 c0 {2,S} {12,S} {13,S} {14,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 *1 H u0 p0 c0 {2,S}
8 *1 H u0 p0 c0 {2,S}
9 H u0 p0 c0 {3,S}
10 *1 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
"""
        )
        labeled = multiple_label_molecule.get_all_labeled_atoms()
        assert "*" in labeled
        assert "*1" in labeled
        assert len(labeled["*"]) == 4
        assert len(labeled["*1"]) == 3

    def test_get_formula(self):
        """
        Test the Molecule.get_all_labeled_atoms() method.
        """
        assert self.molecule[0].get_formula() == "CH2NO2"
        assert self.molecule[1].get_formula() == "CH2NO2"

    def test_radical_count(self):
        """
        Test the Molecule.get_radical_count() method.
        """
        assert self.molecule[0].get_radical_count() == sum([atom.radical_electrons for atom in self.molecule[0].atoms])
        assert self.molecule[1].get_radical_count() == sum([atom.radical_electrons for atom in self.molecule[1].atoms])

    def test_get_molecular_weight(self):
        """
        Test the Molecule.get_molecular_weight() method.
        """
        assert round(abs(self.molecule[0].get_molecular_weight() * 1000 - 60.03), 2) == 0
        assert round(abs(self.molecule[1].get_molecular_weight() * 1000 - 60.03), 2) == 0

    def test_from_adjacency_list(self):
        """
        Test the Molecule.from_adjacency_list() method.
        """

        # molecule 1

        assert self.molecule[0].multiplicity == 2

        atom1 = self.molecule[0].atoms[0]
        atom2 = self.molecule[0].atoms[3]
        atom3 = self.molecule[0].atoms[4]
        atom4 = self.molecule[0].atoms[5]
        assert self.molecule[0].has_bond(atom2, atom1)
        assert self.molecule[0].has_bond(atom2, atom3)
        assert self.molecule[0].has_bond(atom2, atom4)
        assert not self.molecule[0].has_bond(atom1, atom3)
        assert not self.molecule[0].has_bond(atom1, atom4)
        bond21 = atom2.bonds[atom1]
        bond23 = atom2.bonds[atom3]
        bond24 = atom2.bonds[atom4]

        assert atom1.label == "*1"
        assert atom1.element.symbol == "C"
        assert atom1.radical_electrons == 1
        assert atom1.charge == 0

        assert atom2.label == "*2"
        assert atom2.element.symbol == "N"
        assert atom2.radical_electrons == 0
        assert atom2.charge == 1

        assert atom3.label == ""
        assert atom3.element.symbol == "O"
        assert atom3.radical_electrons == 0
        assert atom3.charge == -1

        assert atom4.label == ""
        assert atom4.element.symbol == "O"
        assert atom4.radical_electrons == 0
        assert atom4.charge == 0

        assert bond21.is_single()
        assert bond23.is_single()
        assert bond24.is_double()

        # molecule 2

        assert self.molecule[1].multiplicity == 2

        atom1 = self.molecule[1].atoms[0]
        atom2 = self.molecule[1].atoms[1]
        atom3 = self.molecule[1].atoms[2]
        atom4 = self.molecule[1].atoms[3]
        assert self.molecule[1].has_bond(atom2, atom1)
        assert self.molecule[1].has_bond(atom2, atom3)
        assert self.molecule[1].has_bond(atom2, atom4)
        assert not self.molecule[1].has_bond(atom1, atom3)
        assert not self.molecule[1].has_bond(atom1, atom4)
        bond21 = atom2.bonds[atom1]
        bond23 = atom2.bonds[atom3]
        bond24 = atom2.bonds[atom4]

        assert atom1.label == "*1"
        assert atom1.element.symbol == "C"
        assert atom1.radical_electrons == 1
        assert atom1.charge == 0

        assert atom2.label == "*2"
        assert atom2.element.symbol == "N"
        assert atom2.radical_electrons == 0
        assert atom2.charge == 1

        assert atom3.label == ""
        assert atom3.element.symbol == "O"
        assert atom3.radical_electrons == 0
        assert atom3.charge == -1

        assert atom4.label == ""
        assert atom4.element.symbol == "O"
        assert atom4.radical_electrons == 0
        assert atom4.charge == 0

        assert bond21.is_single()
        assert bond23.is_single()
        assert bond24.is_double()

    def test_to_adjacency_list(self):
        """
        Test the Molecule.to_adjacency_list() method.
        """
        adjlist_1 = self.molecule[0].to_adjacency_list(remove_h=False)
        new_molecule = Molecule().from_adjacency_list(adjlist_1)
        assert self.molecule[0].is_isomorphic(new_molecule)

    def test_isomorphism(self):
        """
        Check the graph isomorphism functions.
        """
        molecule1 = Molecule().from_smiles("C=CC=C[CH]C")
        molecule2 = Molecule().from_smiles("C[CH]C=CC=C")
        assert molecule1.is_isomorphic(molecule2)
        assert molecule2.is_isomorphic(molecule1)

        molecule1 = Molecule().from_adjacency_list(
            """
multiplicity 2
1  *1 C u0 p0 c0 {2,D} {8,S} {9,S}
2  C u0 p0 c0 {1,D} {3,S} {10,S}
3  C u0 p0 c0 {2,S} {4,D} {11,S}
4  C u0 p0 c0 {3,D} {5,S} {12,S}
5  C u1 p0 c0 {4,S} {6,S} {7,S}
6  H u0 p0 c0 {5,S}
7  C u0 p0 c0 {5,S} {13,S} {14,S} {15,S}
8  *2 H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}"""
        )
        molecule2 = Molecule().from_adjacency_list(
            """
multiplicity 2
1  *1 C u0 p0 c0 {2,D} {13,S} {9,S}
2  C u0 p0 c0 {1,D} {3,S} {10,S}
3  C u0 p0 c0 {2,S} {4,D} {11,S}
4  C u0 p0 c0 {3,D} {5,S} {12,S}
5  C u1 p0 c0 {4,S} {6,S} {7,S}
6  H u0 p0 c0 {5,S}
7  C u0 p0 c0 {5,S} {8,S} {14,S} {15,S}
8  H u0 p0 c0 {7,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 *2 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}"""
        )

        assert molecule1.is_isomorphic(molecule2, generate_initial_map=True)
        assert molecule2.is_isomorphic(molecule1, generate_initial_map=True)

    def test_subgraph_isomorphism(self):
        """
        Check the graph isomorphism functions.
        """
        molecule = Molecule().from_smiles("C=CC=C[CH]C")
        group = Group().from_adjacency_list(
            """
        1 Cd u0 p0 c0 {2,D}
        2 Cd u0 p0 c0 {1,D}
        """
        )

        assert molecule.is_subgraph_isomorphic(group)
        mappings = molecule.find_subgraph_isomorphisms(group)
        assert len(mappings) == 4

        for mapping in mappings:
            assert len(mapping) == min(len(molecule.atoms), len(group.atoms))
            for key, value in mapping.items():
                assert key in molecule.atoms
                assert value in group.atoms

    def test_subgraph_isomorphism_again(self):
        molecule = Molecule()
        molecule.from_adjacency_list(
            """
        1 * C u0 p0 c0 {2,D} {7,S} {8,S}
        2   C u0 p0 c0 {1,D} {3,S} {9,S}
        3   C u0 p0 c0 {2,S} {4,D} {10,S}
        4   C u0 p0 c0 {3,D} {5,S} {11,S}
        5   C u0 p0 c0 {4,S} {6,S} {12,S} {13,S}
        6   C u0 p0 c0 {5,S} {14,S} {15,S} {16,S}
        7   H u0 p0 c0 {1,S}
        8   H u0 p0 c0 {1,S}
        9   H u0 p0 c0 {2,S}
        10  H u0 p0 c0 {3,S}
        11  H u0 p0 c0 {4,S}
        12  H u0 p0 c0 {5,S}
        13  H u0 p0 c0 {5,S}
        14  H u0 p0 c0 {6,S}
        15  H u0 p0 c0 {6,S}
        16  H u0 p0 c0 {6,S}
        """
        )

        group = Group()
        group.from_adjacency_list(
            """
        1 * C u0 p0 c0 {2,D} {3,S} {4,S}
        2   C u0 p0 c0 {1,D}
        3   H u0 p0 c0 {1,S}
        4   H u0 p0 c0 {1,S}
        """
        )

        labeled1 = list(molecule.get_all_labeled_atoms().values())[0]
        labeled2 = list(group.get_all_labeled_atoms().values())[0]

        initial_map = {labeled1: labeled2}
        assert molecule.is_subgraph_isomorphic(group, initial_map)

        initial_map = {labeled1: labeled2}
        mappings = molecule.find_subgraph_isomorphisms(group, initial_map)
        assert len(mappings) == 2
        for mapping in mappings:
            assert len(mapping) == min(len(molecule.atoms), len(group.atoms))
            for key, value in mapping.items():
                assert key in molecule.atoms
                assert value in group.atoms

    def test_subgraph_isomorphism_many_labels(self):
        molecule = Molecule()  # specific case (species)
        molecule.from_adjacency_list(
            """
1 *1 C  u1 p0 c0 {2,S} {3,S} {4,S}
2    C  u0 p0 c0 {1,S} {3,S} {5,S} {6,S}
3    C  u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
4    H  u0 p0 c0 {1,S}
5    H  u0 p0 c0 {2,S}
6    H  u0 p0 c0 {2,S}
7    H  u0 p0 c0 {3,S}
8    H  u0 p0 c0 {3,S}
        """
        )

        group = Group()  # general case (functional group)
        group.from_adjacency_list(
            """
1 *1 C   u1 p0 c0 {2,S}, {3,S}
2    R!H u0 p0 c0 {1,S}
3    R!H u0 p0 c0 {1,S}
        """
        )

        labeled1 = molecule.get_all_labeled_atoms()
        labeled2 = group.get_all_labeled_atoms()
        initial_map = {}
        for label, atom1 in labeled1.items():
            initial_map[atom1] = labeled2[label]
        assert molecule.is_subgraph_isomorphic(group, initial_map)

        mappings = molecule.find_subgraph_isomorphisms(group, initial_map)
        assert len(mappings) == 2
        for mapping in mappings:
            assert len(mapping) == min(len(molecule.atoms), len(group.atoms))
            for key, value in mapping.items():
                assert key in molecule.atoms
                assert value in group.atoms

    def test_subgraph_isomorphism_rings(self):
        molecule = Molecule(smiles="C1CCCC1CCC")
        group_no_ring = Group().from_adjacency_list(
            """
1 *1 C u0 p0 c0 r0
        """
        )
        group_ring = Group().from_adjacency_list(
            """
1 *1 C u0 p0 c0 r1
        """
        )

        assert molecule.is_subgraph_isomorphic(group_no_ring)
        mapping = molecule.find_subgraph_isomorphisms(group_no_ring)
        assert len(mapping) == 3
        assert molecule.is_subgraph_isomorphic(group_ring)
        mapping = molecule.find_subgraph_isomorphisms(group_ring)
        assert len(mapping) == 5

    def test_lax_isomorphism(self):
        """Test that we can do isomorphism comparison with strict=False"""
        mol1 = Molecule().from_adjacency_list(
            """
multiplicity 2
1 O u0 p2 c0 {3,D}
2 C u1 p0 c0 {3,S} {4,S} {5,S}
3 C u0 p0 c0 {1,D} {2,S} {6,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
        """
        )

        mol2 = Molecule().from_adjacency_list(
            """
multiplicity 2
1 O u1 p2 c0 {3,S}
2 C u0 p0 c0 {3,D} {4,S} {5,S}
3 C u0 p0 c0 {1,S} {2,D} {6,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
        """
        )

        assert mol1.is_isomorphic(mol2, strict=False)

    def test_adjacency_list(self):
        """
        Check the adjacency list read/write functions for a full molecule.
        """
        molecule1 = Molecule().from_adjacency_list(
            """
        1  C u0 p0 c0 {2,D} {7,S} {8,S}
        2  C u0 p0 c0 {1,D} {3,S} {9,S}
        3  C u0 p0 c0 {2,S} {4,D} {10,S}
        4  C u0 p0 c0 {3,D} {5,S} {11,S}
        5  C u1 {4,S} {6,S} {12,S}
        6  C u0 p0 c0 {5,S} {13,S} {14,S} {15,S}
        7  H u0 p0 c0 {1,S}
        8  H u0 p0 c0 {1,S}
        9  H u0 p0 c0 {2,S}
        10 H u0 p0 c0 {3,S}
        11 H u0 p0 c0 {4,S}
        12 H u0 p0 c0 {5,S}
        13 H u0 p0 c0 {6,S}
        14 H u0 p0 c0 {6,S}
        15 H u0 p0 c0 {6,S}
        """
        )
        molecule2 = Molecule().from_smiles("C=CC=C[CH]C")
        assert molecule1.is_isomorphic(molecule2)
        assert molecule2.is_isomorphic(molecule1)

    def test_generate_h_bonded_structures(self):
        """
        Test that the correct set of Hydrogen Bonded structures are generated
        """
        correct_set = [
            """1  C u0 p0 c0 {2,S} {4,S} {6,S} {7,S}
2  N u0 p1 c0 {1,S} {3,S} {8,S}
3  C u0 p0 c0 {2,S} {9,D} {10,S}
4  O u0 p2 c0 {1,S} {5,S}
5  O u0 p2 c0 {4,S} {8,H} {11,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S} {5,H}
9  O u0 p2 c0 {3,D}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
""",
            """1  C u0 p0 c0 {2,S} {4,S} {6,S} {7,S}
2  N u0 p1 c0 {1,S} {3,S} {8,S} {11,H}
3  C u0 p0 c0 {2,S} {9,D} {10,S}
4  O u0 p2 c0 {1,S} {5,S}
5  O u0 p2 c0 {4,S} {11,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  O u0 p2 c0 {3,D}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {2,H} {5,S}
""",
            """1  C u0 p0 c0 {2,S} {4,S} {6,S} {7,S}
2  N u0 p1 c0 {1,S} {3,S} {8,S} {11,H}
3  C u0 p0 c0 {2,S} {9,D} {10,S}
4  O u0 p2 c0 {1,S} {5,S}
5  O u0 p2 c0 {4,S} {8,H} {11,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S} {5,H}
9  O u0 p2 c0 {3,D}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {2,H} {5,S}
""",
            """1  C u0 p0 c0 {2,S} {4,S} {6,S} {7,S}
2  N u0 p1 c0 {1,S} {3,S} {8,S}
3  C u0 p0 c0 {2,S} {9,D} {10,S}
4  O u0 p2 c0 {1,S} {5,S}
5  O u0 p2 c0 {4,S} {11,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  O u0 p2 c0 {3,D} {11,H}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S} {9,H}
""",
            """1  C u0 p0 c0 {2,S} {4,S} {6,S} {7,S}
2  N u0 p1 c0 {1,S} {3,S} {8,S}
3  C u0 p0 c0 {2,S} {9,D} {10,S}
4  O u0 p2 c0 {1,S} {5,S}
5  O u0 p2 c0 {4,S} {8,H} {11,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S} {5,H}
9  O u0 p2 c0 {3,D} {11,H}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S} {9,H}
""",
        ]

        mols = [Molecule().from_adjacency_list(k) for k in correct_set]

        assert set(mols) == set(self.mHBonds.generate_h_bonded_structures())

    def test_remove_h_bonds(self):
        """
        test that remove HBonds removes all hydrogen bonds from a given molecule
        """
        test_mol = self.mHBonds.generate_h_bonded_structures()[0]
        test_mol.remove_h_bonds()

        for i, atm1 in enumerate(test_mol.atoms):
            for j, atm2 in enumerate(test_mol.atoms):
                if j < i and test_mol.has_bond(atm1, atm2):
                    bd = test_mol.get_bond(atm1, atm2)
                    assert round(abs(bd.order - 0.1), 7) != 0

    def test_sssr(self):
        """
        Test the Molecule.get_smallest_set_of_smallest_rings() method with a complex
        polycyclic molecule.
        """
        molecule = Molecule()
        molecule.from_smiles("C(CC1C(C(CCCCCCCC)C1c1ccccc1)c1ccccc1)CCCCCC")
        # http://cactus.nci.nih.gov/chemical/structure/C(CC1C(C(CCCCCCCC)C1c1ccccc1)c1ccccc1)CCCCCC/image
        sssr = molecule.get_smallest_set_of_smallest_rings()
        assert len(sssr) == 3

    def test_is_in_cycle_ethane(self):
        """
        Test the Molecule is_atom_in_cycle() and is_bond_in_cycle() methods with ethane.
        """
        molecule = Molecule().from_smiles("CC")
        for atom in molecule.atoms:
            assert not molecule.is_atom_in_cycle(atom)
        for atom1 in molecule.atoms:
            for atom2, bond in atom1.bonds.items():
                assert not molecule.is_bond_in_cycle(bond)

    def test_is_in_cycle_cyclohexane(self):
        """
        Test the Molecule is_atom_in_cycle() and is_bond_in_cycle() methods with cyclohexane.
        """
        molecule = Molecule().from_inchi("InChI=1/C6H12/c1-2-4-6-5-3-1/h1-6H2")
        for atom in molecule.atoms:
            if atom.is_hydrogen():
                assert not molecule.is_atom_in_cycle(atom)
            elif atom.is_carbon():
                assert molecule.is_atom_in_cycle(atom)
        for atom1 in molecule.atoms:
            for atom2, bond in atom1.bonds.items():
                if atom1.is_carbon() and atom2.is_carbon():
                    assert molecule.is_bond_in_cycle(bond)
                else:
                    assert not molecule.is_bond_in_cycle(bond)

    def test_from_smiles_h(self):
        """
        Make sure that H radical is produced properly from its SMILES
        representation.
        """
        molecule = Molecule(smiles="[H]")
        assert len(molecule.atoms) == 1
        h = molecule.atoms[0]
        assert h.is_hydrogen()
        assert h.radical_electrons == 1

    def test_from_inchi_h(self):
        """
        Make sure that H radical is produced properly from its InChI
        representation.
        """
        molecule = Molecule().from_inchi("InChI=1/H")
        assert len(molecule.atoms) == 1
        h = molecule.atoms[0]
        assert h.is_hydrogen()
        assert h.radical_electrons == 1

    def test_pickle(self):
        """
        Test that a Molecule object can be successfully pickled and
        unpickled with no loss of information.
        """
        molecule0 = Molecule().from_smiles("C=CC=C[CH2]C")
        molecule0.update()
        import pickle

        molecule = pickle.loads(pickle.dumps(molecule0))

        assert len(molecule0.atoms) == len(molecule.atoms)
        assert molecule0.get_formula() == molecule.get_formula()
        assert molecule0.is_isomorphic(molecule)
        assert molecule.is_isomorphic(molecule0)

    def test_radical_ch(self):
        """
        Test that the species [CH] has one radical electrons and a spin multiplicity of 2.
        """
        molecule = Molecule().from_smiles("[CH]")
        assert molecule.atoms[0].radical_electrons == 1
        assert molecule.multiplicity == 2
        assert molecule.get_radical_count() == 1

    def test_radical_ch2(self):
        """
        Test that the species [CH2] has two radical electrons and a spin multiplicity of 3.
        """
        molecule = Molecule().from_smiles("[CH2]")
        assert molecule.atoms[0].radical_electrons == 2
        assert molecule.multiplicity == 3
        assert molecule.get_radical_count() == 2

    def test_radical_ch2ch2ch2(self):
        """
        Test radical count on [CH2]C[CH2]
        """
        molecule = Molecule().from_smiles("[CH2]C[CH2]")
        assert molecule.get_radical_count() == 2

    def test_singlet_carbene(self):
        """Test radical and carbene count on singlet carbene."""
        mol = Molecule().from_adjacency_list(
            """
1 C u0 p1 {2,S}
2 C u0 p1 {1,S}
""",
            saturate_h=True,
        )
        assert mol.get_radical_count() == 0
        assert mol.get_singlet_carbene_count() == 2

    def test_triplet_carbene(self):
        """Test radical and carbene count on triplet carbene."""
        mol = Molecule().from_adjacency_list(
            """
1 C u2 p0 {2,S}
2 C u0 p1 {1,S}
""",
            saturate_h=True,
        )
        assert mol.get_radical_count() == 2
        assert mol.get_singlet_carbene_count() == 1

    def test_singlet_carbon(self):
        """Test that get_singlet_carbene_count returns 1 for singlet carbon atom."""
        mol = Molecule().from_adjacency_list("1 C u0 p2")
        assert mol.get_singlet_carbene_count() == 1

    def test_smiles(self):
        """
        Test that we can generate a few SMILES strings as expected
        """
        test_strings = [
            "[C-]#[O+]",
            "[C]",
            "[CH]",
            "OO",
            "[H][H]",
            "[H]",
            "[He]",
            "[O]",
            "O",
            "[CH3]",
            "C",
            "[OH]",
            "CCC",
            "CC",
            "N#N",
            "[O]O",
            "C[CH2]",
            "[Ar]",
            "CCCC",
            "O=C=O",
            "[C]#N",
        ]
        for s in test_strings:
            molecule = Molecule(smiles=s)
            assert s == molecule.to_smiles()

    def test_kekule_to_smiles(self):
        """
        Test that we can print SMILES strings of Kekulized structures

        The first two are different Kekule forms of the same thing.
        """
        test_cases = {
            "CC1=C(O)C=CC=C1": """
1 C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
2 C u0 p0 c0 {1,S} {3,D} {4,S}
3 C u0 p0 c0 {2,D} {5,S} {8,S}
4 C u0 p0 c0 {2,S} {7,D} {12,S}
5 C u0 p0 c0 {3,S} {6,D} {13,S}
6 C u0 p0 c0 {5,D} {7,S} {14,S}
7 C u0 p0 c0 {4,D} {6,S} {15,S}
8 O u0 p2 c0 {3,S} {16,S}
9 H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
""",
            "CC1=CC=CC=C1O": """
1 C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
2 C u0 p0 c0 {1,S} {3,S} {4,D}
3 C u0 p0 c0 {2,S} {5,D} {8,S}
4 C u0 p0 c0 {2,D} {7,S} {15,S}
5 C u0 p0 c0 {3,D} {6,S} {12,S}
6 C u0 p0 c0 {5,S} {7,D} {13,S}
7 C u0 p0 c0 {4,S} {6,D} {14,S}
8 O u0 p2 c0 {3,S} {16,S}
9 H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {8,S}
""",
            "CC1=CC=CC=C1": """
1  C u0 p0 c0 {2,D} {6,S} {7,S}
2  C u0 p0 c0 {1,D} {3,S} {8,S}
3  C u0 p0 c0 {2,S} {4,D} {9,S}
4  C u0 p0 c0 {3,D} {5,S} {10,S}
5  C u0 p0 c0 {4,S} {6,D} {11,S}
6  C u0 p0 c0 {1,S} {5,D} {12,S}
7  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
""",
        }
        for smiles, adjlist in test_cases.items():
            m = Molecule().from_adjacency_list(adjlist)
            s = m.to_smiles()
            assert s == smiles, "Generated SMILES string {0} instead of {1}".format(s, smiles)

    def test_kekule_round_trip_smiles(self):
        """
        Test that we can round-trip SMILES strings of Kekulized aromatics
        """
        test_strings = [
            "CC1=CC=CC=C1O",
            "CC1=C(O)C=CC=C1",
            # 'Cc1ccccc1O',  # this will fail because it is Kekulized during from_smiles()
        ]
        for s in test_strings:
            molecule = Molecule(smiles=s)
            assert s == molecule.to_smiles(), "Started with {0} but ended with {1}".format(s, molecule.to_smiles())

    def test_inchi_key(self):
        """
        Test that InChI Key generation is working properly.
        """
        molecule = Molecule().from_inchi("InChI=1S/C7H12/c1-2-7-4-3-6(1)5-7/h6-7H,1-5H2")
        key = molecule.to_inchi_key()
        assert key == "UMRZSTCPUPJPOJ-UHFFFAOYSA-N"

    def test_augmented_inchi(self):
        """
        Test the Augmented InChI generation
        """
        mol = Molecule().from_adjacency_list(
            """
            1     C     u1 p0 c0 {2,S}
            2     C     u1 p0 c0 {1,S}
        """,
            saturate_h=True,
        )

        assert mol.to_augmented_inchi() == "InChI=1S/C2H4/c1-2/h1-2H2/u1,2"

    def test_augmented_inchi_key(self):
        """
        Test the Augmented InChI Key generation
        """
        mol = Molecule().from_adjacency_list(
            """
            1     C     u1 p0 c0 {2,S}
            2     C     u1 p0 c0 {1,S}
        """,
            saturate_h=True,
        )

        assert mol.to_augmented_inchi_key() == "VGGSQFUCUMXWEO-UHFFFAOYSA-N-u1,2"

    def test_linear_methane(self):
        """
        Test the Molecule.is_linear() method.
        """
        assert not Molecule().from_smiles("C").is_linear()

    def test_linear_ethane(self):
        """
        Test the Molecule.is_linear() method.
        """
        assert not Molecule().from_smiles("CC").is_linear()

    def test_linear_propane(self):
        """
        Test the Molecule.is_linear() method.
        """
        assert not Molecule().from_smiles("CCC").is_linear()

    def test_linear_neopentane(self):
        """
        Test the Molecule.is_linear() method.
        """
        assert not Molecule().from_smiles("CC(C)(C)C").is_linear()

    def test_linear_hydrogen(self):
        """
        Test the Molecule.is_linear() method.
        """
        assert not Molecule().from_smiles("[H]").is_linear()

    def test_linear_oxygen(self):
        """
        Test the Molecule.is_linear() method.
        """
        assert Molecule().from_smiles("O=O").is_linear()

    def test_linear_carbon_dioxide(self):
        """
        Test the Molecule.is_linear() method.
        """
        assert Molecule().from_smiles("O=C=O").is_linear()

    def test_linear_acetylene(self):
        """
        Test the Molecule.is_linear() method.
        """
        assert Molecule().from_smiles("C#C").is_linear()

    def test_linear135_hexatriyne(self):
        """
        Test the Molecule.is_linear() method.
        """
        assert Molecule().from_smiles("C#CC#CC#C").is_linear()

    def test_aromatic_benzene(self):
        """
        Test the Molecule.is_aromatic() method for Benzene.
        """
        m = Molecule().from_smiles("C1=CC=CC=C1")
        isomers = m.generate_resonance_structures()
        assert any(isomer.is_aromatic() for isomer in isomers)

    def test_aromatic_naphthalene(self):
        """
        Test the Molecule.is_aromatic() method for Naphthalene.
        """
        m = Molecule().from_smiles("C12C(C=CC=C1)=CC=CC=2")
        isomers = m.generate_resonance_structures()
        assert any(isomer.is_aromatic() for isomer in isomers)

    def test_aromatic_cyclohexane(self):
        """
        Test the Molecule.is_aromatic() method for Cyclohexane.
        """
        m = Molecule().from_smiles("C1CCCCC1")
        isomers = m.generate_resonance_structures()
        assert not any(isomer.is_aromatic() for isomer in isomers)

    def test_heterocyclic_cyclohexanol(self):
        """
        Test the Molecule.is_heterocyclic() method for Cyclohexanol.
        """
        assert not Molecule().from_smiles("OC1CCCCC1").is_heterocyclic()

    def test_heterocyclic_furan(self):
        """
        Test the Molecule.is_heterocyclic() method for Furan.
        """
        assert Molecule().from_smiles("C1C=COC=1").is_heterocyclic()

    def test_heterocyclic_pyridine(self):
        """
        Test the Molecule.is_heterocyclic() method for Pyridine.
        """
        assert Molecule().from_smiles("c1cccnc1").is_heterocyclic()

    def test_count_internal_rotors_ethane(self):
        """
        Test the Molecule.count_internal_rotors() method.
        """
        assert Molecule().from_smiles("CC").count_internal_rotors() == 1

    def test_count_internal_rotors_propane(self):
        """
        Test the Molecule.count_internal_rotors() method.
        """
        assert Molecule().from_smiles("CCC").count_internal_rotors() == 2

    def test_count_internal_rotors_neopentane(self):
        """
        Test the Molecule.count_internal_rotors() method.
        """
        assert Molecule().from_smiles("CC(C)(C)C").count_internal_rotors() == 4

    def test_count_internal_rotors_methyl_cyclohexane(self):
        """
        Test the Molecule.count_internal_rotors() method.
        """
        assert Molecule().from_smiles("C1CCCC1C").count_internal_rotors() == 1

    def test_count_internal_rotors_ethylene(self):
        """
        Test the Molecule.count_internal_rotors() method.
        """
        assert Molecule().from_smiles("C=C").count_internal_rotors() == 0

    def test_count_internal_rotors_acetylene(self):
        """
        Test the Molecule.count_internal_rotors() method.
        """
        assert Molecule().from_smiles("C#C").count_internal_rotors() == 0

    def test_carbene_identifiers(self):
        """
        Test that singlet carbene molecules, bearing an electron pair rather than unpaired electrons
        are correctly converted into rdkit molecules and identifiers.
        """

        ch2_t = """
        multiplicity 3
        1 C u2 p0 c0 {2,S} {3,S}
        2 H u0 p0 c0 {1,S}
        3 H u0 p0 c0 {1,S}
        """

        mol = Molecule().from_adjacency_list(ch2_t)

        assert mol.to_augmented_inchi() == "InChI=1S/CH2/h1H2/u1,1"
        assert mol.to_smiles() == "[CH2]"

        ch2_s = """
        multiplicity 1
        1 C u0 p1 c0 {2,S} {3,S}
        2 H u0 p0 c0 {1,S}
        3 H u0 p0 c0 {1,S}
        """

        mol = Molecule().from_adjacency_list(ch2_s)
        assert mol.to_augmented_inchi() == "InChI=1S/CH2/h1H2/lp1"
        assert mol.to_smiles() == "[CH2]"

    def test_get_symmetry_number(self):
        """
        Test that the symmetry number getter works properly
        """

        mol = Molecule().from_smiles("C")

        assert 12 == mol.get_symmetry_number()

        empty = Molecule()
        assert 1 == empty.get_symmetry_number()

    def test_molecule_props(self):
        """
        Test a key-value pair is added to the props attribute of Molecule.
        """
        self.molecule[0].props["foo"] = "bar"
        assert isinstance(self.molecule[0].props, dict)
        assert self.molecule[0].props["foo"] == "bar"

    def test_molecule_props_object_attribute(self):
        """
        Test that Molecule's props dictionaries are independent of each other.

        Create a test in which is checked whether props is an object attribute rather
        than a class attribute
        """
        spc2 = Molecule()
        self.molecule[0].props["foo"] = "bar"
        spc3 = Molecule()
        spc3.props["foo"] = "bla"
        assert self.molecule[0].props["foo"] == "bar"
        assert spc2.props == {}
        assert spc3.props == {"foo": "bla"}

    @pytest.mark.skip(reason="WIP")
    def test_count_internal_rotors_dimethyl_acetylene(self):
        """
        Test the Molecule.count_internal_rotors() method for dimethylacetylene.

        This is a "hard" test that currently fails.
        """
        assert Molecule().from_smiles("CC#CC").count_internal_rotors() == 1

    def test_saturate_aromatic_radical(self):
        """
        Test that the Molecule.saturate() method works properly for an indenyl radical
        containing Benzene bonds
        """
        indenyl = Molecule().from_adjacency_list(
            """
multiplicity 2
1  C u0 p0 c0 {2,B} {3,S} {4,B}
2  C u0 p0 c0 {1,B} {5,B} {6,S}
3  C u0 p0 c0 {1,S} {7,D} {11,S}
4  C u0 p0 c0 {1,B} {8,B} {12,S}
5  C u0 p0 c0 {2,B} {9,B} {15,S}
6  C u1 p0 c0 {2,S} {7,S} {16,S}
7  C u0 p0 c0 {3,D} {6,S} {10,S}
8  C u0 p0 c0 {4,B} {9,B} {13,S}
9  C u0 p0 c0 {5,B} {8,B} {14,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
"""
        )
        indene = Molecule().from_adjacency_list(
            """
1  C u0 p0 c0 {2,B} {3,S} {4,B}
2  C u0 p0 c0 {1,B} {5,B} {6,S}
3  C u0 p0 c0 {1,S} {7,D} {11,S}
4  C u0 p0 c0 {1,B} {8,B} {12,S}
5  C u0 p0 c0 {2,B} {9,B} {15,S}
6  C u0 p0 c0 {2,S} {7,S} {16,S} {17,S}
7  C u0 p0 c0 {3,D} {6,S} {10,S}
8  C u0 p0 c0 {4,B} {9,B} {13,S}
9  C u0 p0 c0 {5,B} {8,B} {14,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
"""
        )
        saturated_molecule = indenyl.copy(deep=True)
        saturated_molecule.saturate_radicals()
        assert saturated_molecule.is_isomorphic(indene)

    def test_replace_halogen_with_hydrogen(self):
        """
        Test that the Molecule.replace_halogen_with_hydrogen() method works properly for various halogenated molecules
        """
        testCases = [
            # halogenated molecule SMILES, hydrogenated (halogens replaced with hydrogens) molecule SMILES
            ["[F]", "[H]"],
            ["Cl", "[H][H]"],
            ["[Br][Br]", "[H][H]"],
            ["Fc1c(Cl)c(Br)c(I)cc1", "c1ccccc1"],
            ["F[CH]COC(Cl)(Cl)", "[CH2]COC"],
        ]

        for smiles1, smiles2 in testCases:
            mol_halogenated = Molecule().from_smiles(smiles1)
            mol_replaced = mol_halogenated.copy(deep=True)
            mol_replaced.replace_halogen_with_hydrogen()
            mol_replaced_check = Molecule().from_smiles(smiles2)
            assert mol_replaced.is_isomorphic(mol_replaced_check)

    def test_surface_molecules(self):
        """
        Test that we can identify surface molecules.
        """
        adsorbed = Molecule().from_adjacency_list(
            """
                                                1 H u0 p0 c0 {2,S}
                                                2 X u0 p0 c0 {1,S}
                                                """
        )
        assert adsorbed.contains_surface_site()
        gas = Molecule().from_adjacency_list(
            """
                                        1 H u0 p0 c0 {2,S}
                                        2 H u0 p0 c0 {1,S}
                                        """
        )
        assert not gas.contains_surface_site()

        surface_site = Molecule().from_adjacency_list(
            """
                                                1 X u0 p0 c0
                                                """
        )
        assert surface_site.is_surface_site()
        assert not (adsorbed.is_surface_site())
        assert not (gas.is_surface_site())

    def test_is_multidentate(self):
        """
        Test that we can identify a multidentate adsorbate
        """
        monodentate = Molecule().from_adjacency_list(
            """
                                                1 H u0 p0 c0 {2,S}
                                                2 X u0 p0 c0 {1,S}
                                                """
        )


        adjlist = """
1 X u0 p0 c0 {3,S}
2 X u0 p0 c0 {4,S}
3 C u0 p0 c0 {1,S} {4,T}
4 C u0 p0 c0 {2,S} {3,T}
            """

        bidentate = Molecule().from_adjacency_list(adjlist)
        assert not monodentate.is_multidentate()
        assert bidentate.is_multidentate()

    def test_malformed_augmented_inchi(self):
        """Test that augmented inchi without InChI layer raises Exception."""
        malform_aug_inchi = "foo"
        with pytest.raises(InchiException):
            Molecule().from_augmented_inchi(malform_aug_inchi)

    def test_malformed_augmented_inchi_wrong_inchi_layer(self):
        """Test that augmented inchi with wrong layer is caught."""
        malform_aug_inchi = "InChI=1S/CH3/h1H2"
        with pytest.raises(Exception):
            Molecule().from_augmented_inchi(malform_aug_inchi)

    def test_malformed_augmented_inchi_wrong_mult(self):
        """Test that augmented inchi with wrong layer is caught."""
        malform_aug_inchi = "InChI=1S/CH3/h1H3"
        with pytest.raises(Exception):
            Molecule().from_augmented_inchi(malform_aug_inchi)

    def test_malformed_augmented_inchi_wrong_indices(self):
        """Test that augmented inchi with wrong layer is caught."""
        malform_aug_inchi = "InChI=1S/C6H6/c1-3-5-6-4-2/h1,6H,2,5H2/u4,1"
        with pytest.raises(Exception):
            Molecule().from_augmented_inchi(malform_aug_inchi)

    def test_update_lone_pairs(self):
        adjlist = """
1 Si u0 p1 c0 {2,S} {3,S}
2 H  u0 p0 c0 {1,S}
3 H  u0 p0 c0 {1,S}
"""

        mol = Molecule().from_adjacency_list(adjlist)
        mol.update_lone_pairs()
        lp = 0
        for atom in mol.atoms:
            lp += atom.lone_pairs
        assert lp == 1

    def test_large_mol_update(self):
        adjlist = """
1  C u0 p0 c0 {7,S} {33,S} {34,S} {35,S}
2  C u0 p0 c0 {8,S} {36,S} {37,S} {38,S}
3  C u0 p0 c0 {5,S} {9,D} {39,S}
4  C u0 p0 c0 {6,S} {10,D} {40,S}
5  C u0 p0 c0 {3,S} {17,S} {41,S} {85,S}
6  C u0 p0 c0 {4,S} {18,D} {42,S}
7  C u0 p0 c0 {1,S} {11,S} {43,S} {44,S}
8  C u0 p0 c0 {2,S} {12,S} {45,S} {46,S}
9  C u0 p0 c0 {3,D} {31,S} {47,S}
10 C u0 p0 c0 {4,D} {32,S} {48,S}
11 C u0 p0 c0 {7,S} {19,S} {51,S} {52,S}
12 C u0 p0 c0 {8,S} {20,S} {53,S} {54,S}
13 C u0 p0 c0 {18,S} {32,S} {50,S} {86,S}
14 C u0 p0 c0 {17,D} {31,S} {49,S}
15 C u0 p0 c0 {17,S} {25,S} {63,S} {64,S}
16 C u0 p0 c0 {18,S} {26,S} {65,S} {66,S}
17 C u0 p0 c0 {5,S} {14,D} {15,S}
18 C u0 p0 c0 {6,D} {13,S} {16,S}
19 C u0 p0 c0 {11,S} {23,S} {55,S} {56,S}
20 C u0 p0 c0 {12,S} {24,S} {57,S} {58,S}
21 C u0 p0 c0 {25,S} {29,S} {75,S} {76,S}
22 C u0 p0 c0 {26,S} {30,S} {77,S} {78,S}
23 C u0 p0 c0 {19,S} {27,S} {71,S} {72,S}
24 C u0 p0 c0 {20,S} {28,S} {73,S} {74,S}
25 C u0 p0 c0 {15,S} {21,S} {59,S} {60,S}
26 C u0 p0 c0 {16,S} {22,S} {61,S} {62,S}
27 C u0 p0 c0 {23,S} {29,S} {79,S} {80,S}
28 C u0 p0 c0 {24,S} {30,S} {81,S} {82,S}
29 C u0 p0 c0 {21,S} {27,S} {67,S} {68,S}
30 C u0 p0 c0 {22,S} {28,S} {69,S} {70,S}
31 C u0 p0 c0 {9,S} {14,S} {32,S} {83,S}
32 C u0 p0 c0 {10,S} {13,S} {31,S} {84,S}
33 H u0 p0 c0 {1,S}
34 H u0 p0 c0 {1,S}
35 H u0 p0 c0 {1,S}
36 H u0 p0 c0 {2,S}
37 H u0 p0 c0 {2,S}
38 H u0 p0 c0 {2,S}
39 H u0 p0 c0 {3,S}
40 H u0 p0 c0 {4,S}
41 H u0 p0 c0 {5,S}
42 H u0 p0 c0 {6,S}
43 H u0 p0 c0 {7,S}
44 H u0 p0 c0 {7,S}
45 H u0 p0 c0 {8,S}
46 H u0 p0 c0 {8,S}
47 H u0 p0 c0 {9,S}
48 H u0 p0 c0 {10,S}
49 H u0 p0 c0 {14,S}
50 H u0 p0 c0 {13,S}
51 H u0 p0 c0 {11,S}
52 H u0 p0 c0 {11,S}
53 H u0 p0 c0 {12,S}
54 H u0 p0 c0 {12,S}
55 H u0 p0 c0 {19,S}
56 H u0 p0 c0 {19,S}
57 H u0 p0 c0 {20,S}
58 H u0 p0 c0 {20,S}
59 H u0 p0 c0 {25,S}
60 H u0 p0 c0 {25,S}
61 H u0 p0 c0 {26,S}
62 H u0 p0 c0 {26,S}
63 H u0 p0 c0 {15,S}
64 H u0 p0 c0 {15,S}
65 H u0 p0 c0 {16,S}
66 H u0 p0 c0 {16,S}
67 H u0 p0 c0 {29,S}
68 H u0 p0 c0 {29,S}
69 H u0 p0 c0 {30,S}
70 H u0 p0 c0 {30,S}
71 H u0 p0 c0 {23,S}
72 H u0 p0 c0 {23,S}
73 H u0 p0 c0 {24,S}
74 H u0 p0 c0 {24,S}
75 H u0 p0 c0 {21,S}
76 H u0 p0 c0 {21,S}
77 H u0 p0 c0 {22,S}
78 H u0 p0 c0 {22,S}
79 H u0 p0 c0 {27,S}
80 H u0 p0 c0 {27,S}
81 H u0 p0 c0 {28,S}
82 H u0 p0 c0 {28,S}
83 H u0 p0 c0 {31,S}
84 H u0 p0 c0 {32,S}
85 H u0 p0 c0 {5,S}
86 H u0 p0 c0 {13,S}
        """
        mol = Molecule().from_adjacency_list(adjlist)

        mol.reset_connectivity_values()

        try:
            mol.update_connectivity_values()
        except OverflowError:
            assert False, "update_connectivity_values() raised OverflowError unexpectedly!"

    def test_large_mol_creation(self):
        """
        Test molecules between C1 to C201 in 10 carbon intervals to make
        sure that overflow errors are not generated.
        """
        for i in range(1, 202, 10):
            smi = "C" * i
            try:
                Molecule(smiles=smi)
            except OverflowError:
                assert False, "Creation of C{} failed!".format(i)

    def test_get_polycyclic_rings(self):
        """
        Test that polycyclic rings within a molecule are returned properly in the function
        `Graph().get_polycycles()`
        """
        # norbornane
        m1 = Molecule(smiles="C1CC2CCC1C2")
        polyrings1 = m1.get_polycycles()
        assert len(polyrings1) == 1
        ring = polyrings1[0]
        assert len(ring) == 7  # 7 carbons in cycle

        # dibenzyl
        m2 = Molecule(smiles="C1=CC=C(C=C1)CCC1C=CC=CC=1")
        polyrings2 = m2.get_polycycles()
        assert len(polyrings2) == 0

        # spiro[2.5]octane
        m3 = Molecule(smiles="C1CCC2(CC1)CC2")
        polyrings3 = m3.get_polycycles()
        assert len(polyrings3) == 1
        ring = polyrings3[0]
        assert len(ring) == 8

        # 1-phenyl norbornane
        m4 = Molecule(smiles="C1=CC=C(C=C1)C12CCC(CC1)C2")
        polyrings4 = m4.get_polycycles()
        assert len(polyrings4) == 1
        ring = polyrings4[0]
        assert len(ring) == 7

    def test_get_monocyclic_rings(self):
        """
        Test that monocyclic rings within a molecule are returned properly in the function
        `Graph().get_monocycles()`
        """
        m1 = Molecule(smiles="C(CCCC1CCCCC1)CCCC1CCCC1")
        monorings = m1.get_monocycles()
        assert len(monorings) == 2

        m2 = Molecule(smiles="C(CCC1C2CCC1CC2)CC1CCC1")
        monorings = m2.get_monocycles()
        assert len(monorings) == 1
        assert len(monorings[0]) == 4

        m3 = Molecule(smiles="CCCCC")
        monorings = m3.get_monocycles()
        assert len(monorings) == 0

    def test_get_disparate_rings(self):
        """
        Test that monocyclic rings within a molecule are returned properly in the function
        `Graph().get_disparate_cycles()`
        """

        # norbornane
        m1 = Molecule(smiles="C1CC2CCC1C2")
        monorings, polyrings = m1.get_disparate_cycles()
        assert len(monorings) == 0
        assert len(polyrings) == 1
        assert len(polyrings[0]) == 7  # 7 carbons in cycle

        # norbornane + cyclobutane on chain
        m2 = Molecule(smiles="C(CCC1C2CCC1CC2)CC1CCC1")
        monorings, polyrings = m2.get_disparate_cycles()
        assert len(monorings) == 1
        assert len(polyrings) == 1
        assert len(monorings[0]) == 4
        assert len(polyrings[0]) == 7

        # spiro-octane + cyclobutane on chain
        m3 = Molecule(smiles="C1CCC2(CC1)CC2CCCCC1CCC1")
        monorings, polyrings = m3.get_disparate_cycles()
        assert len(polyrings) == 1
        assert len(monorings) == 1
        assert len(monorings[0]) == 4
        assert len(polyrings[0]) == 8

        # butane
        m4 = Molecule(smiles="CCCC")
        monorings, polyrings = m4.get_disparate_cycles()
        assert len(monorings) == 0
        assert len(polyrings) == 0

        # benzene + cyclopropane on chain + cyclopropane on chain
        m5 = Molecule(smiles="C1=CC=C(CCCC2CC2)C(=C1)CCCCCC1CC1")
        monorings, polyrings = m5.get_disparate_cycles()
        assert len(monorings) == 3
        assert len(polyrings) == 0

        # octacene
        m6 = Molecule(smiles="c1ccc2cc3cc4cc5cc6cc7cc8ccccc8cc7cc6cc5cc4cc3cc2c1")
        monorings, polyrings = m6.get_disparate_cycles()
        assert len(monorings) == 0
        assert len(polyrings) == 1
        assert len(polyrings[0]) == 34

        # JP-10
        m7 = Molecule(smiles="C1CC2C3CCC(C3)C2C1")
        monorings, polyrings = m7.get_disparate_cycles()
        assert len(monorings) == 0
        assert len(polyrings) == 1
        assert len(polyrings[0]) == 10

    def test_get_smallest_set_of_smallest_rings(self):
        """
        Test that SSSR within a molecule are returned properly in the function
        `Graph().get_smallest_set_of_smallest_rings()`
        """

        m1 = Molecule(smiles="C12CCC1C3CC2CC3")
        sssr1 = m1.get_smallest_set_of_smallest_rings()
        sssr1_sizes = sorted([len(ring) for ring in sssr1])
        sssr1_sizes_expected = [4, 5, 5]
        assert sssr1_sizes == sssr1_sizes_expected

        m2 = Molecule(smiles="C1(CC2)C(CC3)CC3C2C1")
        sssr2 = m2.get_smallest_set_of_smallest_rings()
        sssr2_sizes = sorted([len(ring) for ring in sssr2])
        sssr2_sizes_expected = [5, 5, 6]
        assert sssr2_sizes == sssr2_sizes_expected

        m3 = Molecule(smiles="C1(CC2)C2C(CCCC3)C3C1")
        sssr3 = m3.get_smallest_set_of_smallest_rings()
        sssr3_sizes = sorted([len(ring) for ring in sssr3])
        sssr3_sizes_expected = [4, 5, 6]
        assert sssr3_sizes == sssr3_sizes_expected

        m4 = Molecule(smiles="C12=CC=CC=C1C3=C2C=CC=C3")
        sssr4 = m4.get_smallest_set_of_smallest_rings()
        sssr4_sizes = sorted([len(ring) for ring in sssr4])
        sssr4_sizes_expected = [4, 6, 6]
        assert sssr4_sizes == sssr4_sizes_expected

        m5 = Molecule(smiles="C12=CC=CC=C1CC3=C(C=CC=C3)C2")
        sssr5 = m5.get_smallest_set_of_smallest_rings()
        sssr5_sizes = sorted([len(ring) for ring in sssr5])
        sssr5_sizes_expected = [6, 6, 6]
        assert sssr5_sizes == sssr5_sizes_expected

    def test_get_deterministic_smallest_set_of_smallest_rings_case1(self):
        """
        Test fused tricyclic can be decomposed into single rings more
        deterministically
        """
        smiles = "C1C2C3C=CCCC2C13"

        previous_num_shared_atoms_list = None
        # repeat 100 time to test non-deterministic behavior
        for _ in range(100):
            mol = Molecule().from_smiles(smiles)
            sssr_det = mol.get_deterministic_sssr()

            num_shared_atoms_list = []
            for i, ring_i in enumerate(sssr_det):
                for j in range(i + 1, len(sssr_det)):
                    ring_j = sssr_det[j]
                    num_shared_atoms = len(set(ring_i).intersection(ring_j))

                    num_shared_atoms_list.append(num_shared_atoms)

            num_shared_atoms_list = sorted(num_shared_atoms_list)

            if previous_num_shared_atoms_list is None:
                previous_num_shared_atoms_list = num_shared_atoms_list
                continue
            assert num_shared_atoms_list == previous_num_shared_atoms_list
            previous_num_shared_atoms_list = num_shared_atoms_list

    def test_get_deterministic_smallest_set_of_smallest_rings_case2(self):
        """
        Test if two possible smallest rings can join the smallest set
        the method can pick one of them deterministically using sum of
        atomic numbers along the rings.
        In this test case and with currect method setup, ring (CCSCCCCC)
        will be picked rather than ring(CCCOCC).
        """

        smiles = "C1=CC2C3CSC(CO3)C2C1"

        previous_atom_symbols_list = None
        # repeat 100 time to test non-deterministic behavior
        for _ in range(100):
            mol = Molecule().from_smiles(smiles)
            sssr_det = mol.get_deterministic_sssr()

            atom_symbols_list = []
            for ring in sssr_det:
                atom_symbols = sorted([a.element.symbol for a in ring])
                atom_symbols_list.append(atom_symbols)

            atom_symbols_list = sorted(atom_symbols_list)

            if previous_atom_symbols_list is None:
                previous_atom_symbols_list = atom_symbols_list
                continue
            assert atom_symbols_list == previous_atom_symbols_list
            previous_atom_symbols_list = atom_symbols_list

    def test_get_deterministic_smallest_set_of_smallest_rings_case3(self):
        """
        Test if two possible smallest rings can join the smallest set
        the method can pick one of them deterministically when their
        sum of atomic numbers along the rings are also equal to each other.

        To break the tie, one option we have is to consider adding contributions
        from other parts of the molecule, such as atomic number weighted connectivity
        value and differentiate bond orders when calculating connectivity values.
        """
        smiles = "C=1CC2C3CSC(O[Si]3)C2C1"

        previous_atom_symbols_list = None
        # repeat 100 time to test non-deterministic behavior
        for _ in range(100):
            mol = Molecule().from_smiles(smiles)
            sssr_det = mol.get_deterministic_sssr()

            atom_symbols_list = []
            for ring in sssr_det:
                atom_symbols = sorted([a.element.symbol for a in ring])
                atom_symbols_list.append(atom_symbols)

            atom_symbols_list = sorted(atom_symbols_list)

            if previous_atom_symbols_list is None:
                previous_atom_symbols_list = atom_symbols_list
                continue
            assert atom_symbols_list == previous_atom_symbols_list
            previous_atom_symbols_list = atom_symbols_list

    def test_to_group(self):
        """
        Test if we can convert a Molecule object into a Group object.
        """
        mol = Molecule().from_smiles("CC(C)CCCC(C)C1CCC2C3CC=C4CC(O)CCC4(C)C3CCC12C")  # cholesterol
        mol.atoms[0].label = "*1"
        mol.atoms[1].label = "*2"
        group = mol.to_group()

        assert isinstance(group, Group)

        assert len(mol.atoms) == len(group.atoms)

        molbondcount = sum([1 for atom in mol.atoms for _ in atom.edges.items()])
        groupbondcount = sum([1 for atom in group.atoms for _ in atom.edges.items()])
        assert molbondcount == groupbondcount

        for i, molAt in enumerate(mol.atoms):
            group_atom = group.atoms[i]
            assert group_atom.label == molAt.label
            atom_types = [groupAtomType.equivalent(molAt.atomtype) for groupAtomType in group_atom.atomtype]
            assert any(atom_types)

    def test_to_adjacency_list_with_isotopes(self):
        """
        Test the Molecule.to_adjacency_list() method works for atoms with unexpected isotopes.
        """

        mol = Molecule().from_smiles("CC")
        mol.atoms[0].element = get_element("C", 13)

        table = str.maketrans({"\n": None, " ": None})  # Translation table to remove whitespace

        adjlist = mol.to_adjacency_list().translate(table)
        adjlist_exp = """
        1 C u0 p0 c0 i13 {2,S} {3,S} {4,S} {5,S}
        2 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
        3 H u0 p0 c0 {1,S}
        4 H u0 p0 c0 {1,S}
        5 H u0 p0 c0 {1,S}
        6 H u0 p0 c0 {2,S}
        7 H u0 p0 c0 {2,S}
        8 H u0 p0 c0 {2,S}
        """.translate(
            table
        )

        assert adjlist == adjlist_exp

        mol = Molecule().from_smiles("CC")
        mol.atoms[2].element = get_element("H", 2)

        adjlist = mol.to_adjacency_list().translate(table)
        adjlist_exp = """
        1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
        2 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
        3 H u0 p0 c0 i2 {1,S}
        4 H u0 p0 c0 {1,S}
        5 H u0 p0 c0 {1,S}
        6 H u0 p0 c0 {2,S}
        7 H u0 p0 c0 {2,S}
        8 H u0 p0 c0 {2,S}
        """.translate(
            table
        )

        assert adjlist == adjlist_exp

        mol = Molecule().from_smiles("OC")
        mol.atoms[0].element = get_element("O", 18)

        adjlist = mol.to_adjacency_list().translate(table)
        adjlist_exp = """
        1 O u0 p2 c0 i18 {2,S} {3,S}
        2 C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
        3 H u0 p0 c0 {1,S}
        4 H u0 p0 c0 {2,S}
        5 H u0 p0 c0 {2,S}
        6 H u0 p0 c0 {2,S}
        """.translate(
            table
        )

        assert adjlist == adjlist_exp

    def test_from_adjacency_list_with_isotopes(self):
        """
        Test the Molecule.from_adjacency_list() method works for atoms with unexpected isotopes.
        """

        exp = Molecule().from_smiles("CC")
        exp.atoms[0].element = get_element("C", 13)

        adjlist_calc = """
        1 C u0 p0 c0 i13 {2,S} {3,S} {4,S} {5,S}
        2 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
        3 H u0 p0 c0 {1,S}
        4 H u0 p0 c0 {1,S}
        5 H u0 p0 c0 {1,S}
        6 H u0 p0 c0 {2,S}
        7 H u0 p0 c0 {2,S}
        8 H u0 p0 c0 {2,S}
        """
        calc = Molecule().from_adjacency_list(adjlist_calc)

        assert exp.is_isomorphic(calc)

        exp = Molecule().from_smiles("CC")
        exp.atoms[2].element = get_element("H", 2)

        adjlist_calc = """
        1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
        2 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
        3 H u0 p0 c0 i2 {1,S}
        4 H u0 p0 c0 {1,S}
        5 H u0 p0 c0 {1,S}
        6 H u0 p0 c0 {2,S}
        7 H u0 p0 c0 {2,S}
        8 H u0 p0 c0 {2,S}
        """
        calc = Molecule().from_adjacency_list(adjlist_calc)

        assert exp.is_isomorphic(calc)

        exp = Molecule().from_smiles("OC")
        exp.atoms[0].element = get_element("O", 18)

        adjlist_calc = """
        1 O u0 p2 c0 i18 {2,S} {3,S}
        2 C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
        3 H u0 p0 c0 {1,S}
        4 H u0 p0 c0 {2,S}
        5 H u0 p0 c0 {2,S}
        6 H u0 p0 c0 {2,S}
        """
        calc = Molecule().from_adjacency_list(adjlist_calc)

        assert exp.is_isomorphic(calc)

    def test_aromaticity_perception_benzene(self):
        """Test aromaticity perception via get_aromatic_rings for benzene."""
        mol = Molecule(smiles="c1ccccc1")
        aromatic_atoms, aromatic_bonds = mol.get_aromatic_rings()
        assert len(aromatic_atoms) == 1
        assert len(aromatic_bonds) == 1
        for bond in aromatic_bonds[0]:
            assert bond.atom1 in aromatic_atoms[0] and bond.atom2 in aromatic_atoms[0]

    def test_aromaticity_perception_tetralin(self):
        """Test aromaticity perception via get_aromatic_rings for tetralin."""
        mol = Molecule(smiles="c1ccc2c(c1)CCCC2")
        aromatic_atoms, aromatic_bonds = mol.get_aromatic_rings()
        assert len(aromatic_atoms) == 1
        assert len(aromatic_bonds) == 1
        for bond in aromatic_bonds[0]:
            assert bond.atom1 in aromatic_atoms[0] and bond.atom2 in aromatic_atoms[0]

    def test_aromaticity_perception_biphenyl(self):
        """Test aromaticity perception via get_aromatic_rings for biphenyl."""
        mol = Molecule(smiles="c1ccc(cc1)c2ccccc2")
        aromatic_atoms, aromatic_bonds = mol.get_aromatic_rings()
        assert len(aromatic_atoms) == 2
        assert len(aromatic_bonds) == 2
        for index in range(len(aromatic_atoms)):
            for bond in aromatic_bonds[index]:
                assert bond.atom1 in aromatic_atoms[index] and bond.atom2 in aromatic_atoms[index]

    def test_aromaticity_perception_azulene(self):
        """Test aromaticity perception via get_aromatic_rings for azulene."""
        mol = Molecule(smiles="c1cccc2cccc2c1")
        aromatic_atoms, aromatic_bonds = mol.get_aromatic_rings()
        assert len(aromatic_atoms) == 0
        assert len(aromatic_bonds) == 0

    def test_aromaticity_perception_furan(self):
        """Test aromaticity perception via get_aromatic_rings for furan."""
        mol = Molecule(smiles="c1ccoc1")
        aromatic_atoms, aromatic_bonds = mol.get_aromatic_rings()
        assert len(aromatic_atoms) == 0
        assert len(aromatic_bonds) == 0

    def test_aromaticity_perception_benzonaphthalene(self):
        """Test aromaticity perception via get_aromatic_rings for benzonaphthalene with multiple fused bonds."""
        mol = Molecule(smiles="c1cc2ccc3ccc(c1)c2c3")
        aromatic_atoms, aromatic_bonds = mol.get_aromatic_rings()
        assert len(aromatic_atoms) == 1
        assert len(aromatic_bonds) == 1

    def test_aromaticity_perception_save_order(self):
        """Test aromaticity perception via get_aromatic_rings for phenyl radical without changing atom order."""
        mol = Molecule().from_adjacency_list(
            """multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {7,D}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,D} {8,S}
4  C u0 p0 c0 {3,D} {5,S} {9,S}
5  C u0 p0 c0 {4,S} {6,D} {10,S}
6  C u0 p0 c0 {5,D} {7,S} {11,S}
7  C u0 p0 c0 {1,D} {6,S} {12,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""
        )
        aromatic_atoms, aromatic_bonds = mol.get_aromatic_rings(save_order=True)
        assert len(aromatic_atoms) == 1
        assert len(aromatic_bonds) == 1
        # A quick check for non-changed atom order is to check
        # if the first atom becomes the oxygen atom after calling `get_aromatic_rings`
        assert not mol.atoms[0].is_oxygen()

    def test_aryl_radical_true(self):
        """Test aryl radical perception for phenyl radical."""
        mol = Molecule(smiles="[c]1ccccc1")
        assert mol.is_aryl_radical()

    def test_has_halogen(self):
        """Test Molecule.has_halogen() method."""
        mol1 = Molecule(smiles="CCCCl")
        mol2 = Molecule(smiles="CCCCC")
        assert mol1.has_halogen()
        assert not mol2.has_halogen()

    def test_aryl_radical_false(self):
        """Test aryl radical perception for benzyl radical."""
        mol = Molecule(smiles="[CH2]c1ccccc1")
        assert not mol.is_aryl_radical()

    def test_aryl_radical_birad(self):
        """Test aryl radical perception for biradical species.

        This is a case that is not properly handled right now, since a single boolean cannot
        characterize multiple radicals. In such cases, the method will return false if
        any of the radicals is not an aryl radical."""
        mol = Molecule(smiles="[CH2]c1c[c]ccc1")
        assert not mol.is_aryl_radical()

    def test_identical_true(self):
        """Test that the is_identical returns True with butane"""
        mol = Molecule(smiles="CCCC")
        mol.assign_atom_ids()
        mol_copy = mol.copy(deep=True)
        assert mol.is_isomorphic(mol_copy)
        assert mol.is_identical(mol_copy)

    def test_identical_true2(self):
        """Test that is_identical with strict=False returns True with allyl"""
        mol = Molecule(smiles="C=C[CH2]")
        mol.assign_atom_ids()
        res = mol.generate_resonance_structures(keep_isomorphic=True)
        assert len(res) == 2

        mol2 = res[1]
        assert mol.is_isomorphic(mol2)
        assert not mol.is_identical(mol2)
        assert mol.is_identical(mol2, strict=False)

    def test_identical_false(self):
        """Test that the is_identical returns False with butane"""
        mol = Molecule(smiles="CCCC")
        mol.assign_atom_ids()
        mol_copy = mol.copy(deep=True)
        # Remove a hydrogen from mol
        a = mol.atoms[-1]

        mol.remove_atom(a)
        # Remove a different hydrogen from mol_copy
        b = mol_copy.atoms[-2]

        mol_copy.remove_atom(b)

        assert mol.is_isomorphic(mol_copy)
        assert not mol.is_identical(mol_copy)

    def test_identical_false2(self):
        """Test that the is_identical method returns False with ethene"""
        # Manually test addition of H radical to ethene
        reactant1 = Molecule(smiles="C=C")
        carbons = [atom for atom in reactant1.atoms if atom.symbol == "C"]
        carbons[0].label = "*1"
        carbons[1].label = "*2"
        reactant2 = Molecule(smiles="[H]")
        reactant2.atoms[0].label = "*3"
        # Merge reactants
        mol = reactant1.merge(reactant2)
        mol.assign_atom_ids()
        mol_copy = mol.copy(deep=True)
        # Manually perform R_Addition_MultipleBond of *3 to *1
        labeled_atoms = mol.get_all_labeled_atoms()
        mol.get_bond(labeled_atoms["*1"], labeled_atoms["*2"]).decrement_order()
        mol.add_bond(Bond(labeled_atoms["*1"], labeled_atoms["*3"], order="S"))
        labeled_atoms["*2"].increment_radical()
        labeled_atoms["*3"].decrement_radical()
        # Manually perform R_Addition_MultipleBond of *3 to *2
        labeled_atoms = mol_copy.get_all_labeled_atoms()
        mol_copy.get_bond(labeled_atoms["*1"], labeled_atoms["*2"]).decrement_order()
        mol_copy.add_bond(Bond(labeled_atoms["*2"], labeled_atoms["*3"], order="S"))
        labeled_atoms["*1"].increment_radical()
        labeled_atoms["*3"].decrement_radical()

        assert mol.is_isomorphic(mol_copy)
        assert not mol.is_identical(mol_copy)

    def test_atom_id_valid(self):
        """see if the atomIDVvalid method properly returns True"""
        mol = Molecule(smiles="CCCC")
        for index, atom in enumerate(mol.atoms):
            atom.id = index
        assert mol.atom_ids_valid()

    def test_atom_id_valid2(self):
        """see if the atomIDVvalid method properly returns False"""
        mol = Molecule(smiles="CCCC")
        for index, atom in enumerate(mol.atoms):
            atom.id = index
        mol.atoms[3].id = 4
        assert not mol.atom_ids_valid()

    def test_atom_id_valid3(self):
        """see if the atomIDVvalid method properly returns False"""
        mol = Molecule(smiles="CCCC")
        assert not mol.atom_ids_valid()

    def test_assign_atom_id(self):
        """see if the assignAtomID method properly labels molecule"""
        mol = Molecule(smiles="CCCC")
        mol.assign_atom_ids()
        assert mol.atom_ids_valid()

    def test_fingerprint_property(self):
        """Test that the Molecule.fingerprint property works"""
        # Test getting fingerprint
        assert self.molecule[0].fingerprint == "C01H02N01O02S00"

        # Test setting fingerprint
        self.molecule[0].fingerprint = "nitronate"
        assert self.molecule[0].fingerprint == "nitronate"

    def test_fingerprint_property_more_elements(self):
        """Test that the Molecule.fingerprint property is consistent with many elements"""
        mol1 = Molecule().from_adjacency_list(
            """
1 Cl u0 p3 c0 {2,S}
2 C  u0 p0 c0 {1,S} {3,S} {4,S} {6,S}
3 F  u0 p3 c0 {2,S}
4 O  u0 p2 c0 {2,S} {5,S}
5 O  u0 p2 c0 {4,S} {7,S}
6 H  u0 p0 c0 {2,S}
7 H  u0 p0 c0 {5,S}
"""
        )
        mol2 = Molecule().from_adjacency_list(
            """
1 Cl u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {4,S} {5,S}
4 O u0 p2 c0 {3,S} {7,S}
5 C u0 p0 c0 {1,S} {2,S} {3,S} {6,S}
6 H u0 p0 c0 {5,S}
7 H u0 p0 c0 {4,S}
"""
        )
        # Confirm that atom orders are different
        mol1_atoms = "".join([atom.symbol for atom in mol1.atoms])
        mol2_atoms = "".join([atom.symbol for atom in mol2.atoms])
        assert mol1_atoms != mol2_atoms

        # Test getting fingerprint
        expected = "C01H02N00O02S00Cl01F01"
        assert mol1.fingerprint == expected
        assert mol2.fingerprint == expected

    def test_saturate_unfilled_valence(self):
        """
        Test the saturateUnfilledValence for an aromatic and nonaromatic case
        """
        # test butane
        expected = Molecule(smiles="CCCC")
        test = expected.copy(deep=True)
        test.delete_hydrogens()

        hydrogens = 0
        for atom in test.atoms:
            if atom.is_hydrogen():
                hydrogens += 1
        assert hydrogens == 0

        test.saturate_unfilled_valence()

        hydrogens = 0
        for atom in test.atoms:
            if atom.is_hydrogen():
                hydrogens += 1
        assert hydrogens == 10

        test.update()
        assert expected.is_isomorphic(test)

        # test benzene
        expected = Molecule(smiles="c1ccccc1")
        test = expected.copy(deep=True)
        test.delete_hydrogens()
        hydrogens = 0
        for atom in test.atoms:
            if atom.is_hydrogen():
                hydrogens += 1
        assert hydrogens == 0

        test.saturate_unfilled_valence()

        hydrogens = 0
        for atom in test.atoms:
            if atom.is_hydrogen():
                hydrogens += 1
        assert hydrogens == 6

        test.update()
        assert expected.is_isomorphic(test)

    def test_get_element_count(self):
        """Test that we can count elements properly."""
        mol1 = Molecule(smiles="c1ccccc1")
        expected1 = {"C": 6, "H": 6}
        result1 = mol1.get_element_count()
        assert expected1 == result1

        mol2 = Molecule(smiles="CS(C)(=O)=O")
        expected2 = {"C": 2, "H": 6, "O": 2, "S": 1}
        result2 = mol2.get_element_count()
        assert expected2 == result2

        mol3 = Molecule(smiles="CCN")
        expected3 = {"C": 2, "H": 7, "N": 1}
        result3 = mol3.get_element_count()
        assert expected3 == result3

    def test_ring_perception(self):
        """Test that identifying ring membership of atoms works properly."""
        mol = Molecule(smiles="c12ccccc1cccc2")
        mol.identify_ring_membership()
        for atom in mol.atoms:
            if atom.element == "C":
                assert atom.props["inRing"]
            elif atom.element == "H":
                assert not atom.props["inRing"]

    def test_enumerate_bonds(self):
        """Test that generating a count of bond labels works properly."""
        adj_list = """
        1  O u0 p2 c0 {4,S} {23,S} {24,H}
        2  O u0 p2 c0 {8,S} {23,H} {24,S}
        3  C u0 p0 c0 {4,S} {6,S} {14,S} {15,S}
        4  C u0 p0 c0 {1,S} {3,S} {16,S} {17,S}
        5  C u0 p0 c0 {7,S} {12,S} {18,S} {19,S}
        6  C u0 p0 c0 {3,S} {8,B} {9,B}
        7  C u0 p0 c0 {5,S} {9,B} {10,B}
        8  C u0 p0 c0 {2,S} {6,B} {11,B}
        9  C u0 p0 c0 {6,B} {7,B} {22,S}
        10 C u0 p0 c0 {7,B} {11,B} {20,S}
        11 C u0 p0 c0 {8,B} {10,B} {21,S}
        12 C u0 p0 c0 {5,S} {13,T}
        13 C u0 p0 c0 {12,T} {25,S}
        14 H u0 p0 c0 {3,S}
        15 H u0 p0 c0 {3,S}
        16 H u0 p0 c0 {4,S}
        17 H u0 p0 c0 {4,S}
        18 H u0 p0 c0 {5,S}
        19 H u0 p0 c0 {5,S}
        20 H u0 p0 c0 {10,S}
        21 H u0 p0 c0 {11,S}
        22 H u0 p0 c0 {9,S}
        23 H u0 p0 c0 {1,S} {2,H}
        24 H u0 p0 c0 {1,H} {2,S}
        25 H u0 p0 c0 {13,S}
        """
        mol = Molecule().from_adjacency_list(adj_list)
        bonds = mol.enumerate_bonds()
        assert bonds["C#C"] == 1
        assert bonds["C-C"] == 4
        assert bonds["C-H"] == 10
        assert bonds["C-O"] == 2
        assert bonds["C:C"] == 6
        assert bonds["H-O"] == 2
        assert bonds["H~O"] == 2

    def test_count_aromatic_rings(self):
        """Test that we can count aromatic rings correctly."""
        mol = Molecule(smiles="c12ccccc1cccc2")
        out = mol.generate_resonance_structures()
        result = [m.count_aromatic_rings() for m in out]

        assert result == [2, 1, 0]

    def test_remove_van_der_waals_bonds(self):
        """Test we can remove a van-der-Waals bond"""
        adjlist = """
1 X  u0 p0 c0 {2,vdW}
2 H  u0 p0 c0 {1,vdW}, {3,S}
3 H  u0 p0 c0 {2,S}
"""
        mol = Molecule().from_adjacency_list(adjlist)
        assert len(mol.get_all_edges()) == 2
        mol.remove_van_der_waals_bonds()
        assert len(mol.get_all_edges()) == 1
