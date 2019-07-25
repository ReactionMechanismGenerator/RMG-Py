#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2019 Prof. William H. Green (whgreen@mit.edu),           #
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
from .molecule import Atom, Bond, Molecule
from .group import Group, ActionError
from .element import getElement, elementList

################################################################################
class TestAtom(unittest.TestCase):
    """
    Contains unit tests of the Atom class.
    """

    def setUp(self):
        """
        A method called before each unit test in this class.
        """
        self.atom = Atom(element=getElement('C'), radicalElectrons=1, charge=0, label='*1', lonePairs=0)
    
    def testMass(self):
        """
        Test the Atom.mass property.
        """
        self.assertTrue(self.atom.mass == self.atom.element.mass)
    
    def testNumber(self):
        """
        Test the Atom.number property.
        """
        self.assertTrue(self.atom.number == self.atom.element.number)
    
    def testSymbol(self):
        """
        Test the Atom.symbol property.
        """
        self.assertTrue(self.atom.symbol == self.atom.element.symbol)
    
    def testIsHydrogen(self):
        """
        Test the Atom.isHydrogen() method.
        """
        for element in elementList:
            atom = Atom(element=element, radicalElectrons=1, charge=0, label='*1', lonePairs=0)
            if element.symbol == 'H':
                self.assertTrue(atom.isHydrogen())
            else:
                self.assertFalse(atom.isHydrogen())
    
    def testIsNonHydrogen(self):
        """
        Test the Atom.isNonHydrogen() method.
        """
        for element in elementList:
            atom = Atom(element=element, radicalElectrons=1, charge=0, label='*1', lonePairs=0)
            if element.symbol == 'H':
                self.assertFalse(atom.isNonHydrogen())
            else:
                self.assertTrue(atom.isNonHydrogen(), "Atom {0!r} isn't reporting isNonHydrogen()".format(atom))

    def testIsCarbon(self):
        """
        Test the Atom.isCarbon() method.
        """
        for element in elementList:
            atom = Atom(element=element, radicalElectrons=1, charge=0, label='*1', lonePairs=0)
            if element.symbol == 'C':
                self.assertTrue(atom.isCarbon())
            else:
                self.assertFalse(atom.isCarbon())

    def testIsSilicon(self):
        """
        Test the Atom.isSilicon() method.
        """
        for element in elementList:
            atom = Atom(element=element, radicalElectrons=1, charge=0, label='*1', lonePairs=0)
            if element.symbol == 'Si':
                self.assertTrue(atom.isSilicon())
            else:
                self.assertFalse(atom.isSilicon())

    def testIsOxygen(self):
        """
        Test the Atom.isOxygen() method.
        """
        for element in elementList:
            atom = Atom(element=element, radicalElectrons=2, charge=0, label='*1', lonePairs=2)
            if element.symbol == 'O':
                self.assertTrue(atom.isOxygen())
            else:
                self.assertFalse(atom.isOxygen())

    def testIsNitrogen(self):
        """
        Test the Atom.isNitrogen() method.
        """
        for element in elementList:
            atom = Atom(element=element, radicalElectrons=1, charge=0, label='*1', lonePairs=1)
            if element.symbol == 'N':
                self.assertTrue(atom.isNitrogen())
            else:
                self.assertFalse(atom.isNitrogen())

    def testIsSulfur(self):
        """
        Test the Atom.isSulfur() method.
        """
        for element in elementList:
            atom = Atom(element=element, radicalElectrons=1, charge=0, label='*1', lonePairs=2)
            if element.symbol == 'S':
                self.assertTrue(atom.isSulfur())
            else:
                self.assertFalse(atom.isSulfur())

    def testIsFluorine(self):
        """
        Test the Atom.isFluorine() method.
        """
        for element in elementList:
            atom = Atom(element=element, radicalElectrons=1, charge=0, label='*1', lonePairs=3)
            if element.symbol == 'F':
                self.assertTrue(atom.isFluorine())
            else:
                self.assertFalse(atom.isFluorine())

    def testIsChlorine(self):
        """
        Test the Atom.isChlorine() method.
        """
        for element in elementList:
            atom = Atom(element=element, radicalElectrons=1, charge=0, label='*1', lonePairs=3)
            if element.symbol == 'Cl':
                self.assertTrue(atom.isChlorine())
            else:
                self.assertFalse(atom.isChlorine())

    def testIsIodine(self):
        """
        Test the Atom.isIodine() method.
        """
        for element in elementList:
            atom = Atom(element=element, radicalElectrons=1, charge=0, label='*1', lonePairs=3)
            if element.symbol == 'I':
                self.assertTrue(atom.isIodine())
            else:
                self.assertFalse(atom.isIodine())

    def testIsNOS(self):
        """
        Test the Atom.isNOS() method.
        """
        for element in elementList:
            atom = Atom(element=element, radicalElectrons=1, charge=0, label='*1', lonePairs=2)
            if element.symbol in ['N', 'O', 'S']:
                self.assertTrue(atom.isNOS())
            else:
                self.assertFalse(atom.isNOS())

    def testIsSurfaceSite(self):
        """
        Test the Atom.isSurfaceSite() method.
        """
        for element in elementList:
            atom = Atom(element=element, radicalElectrons=0, charge=0, label='*1', lonePairs=0)
            if element.symbol == 'X':
                self.assertTrue(atom.isSurfaceSite())
            else:
                self.assertFalse(atom.isSurfaceSite())

    def testIncrementRadical(self):
        """
        Test the Atom.incrementRadical() method.
        """
        radicalElectrons = self.atom.radicalElectrons
        self.atom.incrementRadical()
        self.assertEqual(self.atom.radicalElectrons, radicalElectrons + 1)
    
    def testDecrementRadical(self):
        """
        Test the Atom.decrementRadical() method.
        """
        radicalElectrons = self.atom.radicalElectrons
        self.atom.decrementRadical()
        self.assertEqual(self.atom.radicalElectrons, radicalElectrons - 1)
           
    def testApplyActionBreakBond(self):
        """
        Test the Atom.applyAction() method for a BREAK_BOND action.
        """
        action = ['BREAK_BOND', '*1', 1, '*2']
        for element in elementList:
            atom0 = Atom(element=element, radicalElectrons=1, charge=0, label='*1', lonePairs=0)
            atom = atom0.copy()
            atom.applyAction(action)
            self.assertEqual(atom0.element, atom.element)
            self.assertEqual(atom0.radicalElectrons, atom.radicalElectrons)
            self.assertEqual(atom0.charge, atom.charge)
            self.assertEqual(atom0.label, atom.label)
    
    def testApplyActionFormBond(self):
        """
        Test the Atom.applyAction() method for a FORM_BOND action.
        """
        action = ['FORM_BOND', '*1', 1, '*2']
        for element in elementList:
            atom0 = Atom(element=element, radicalElectrons=1, charge=0, label='*1', lonePairs=0)
            atom = atom0.copy()
            atom.applyAction(action)
            self.assertEqual(atom0.element, atom.element)
            self.assertEqual(atom0.radicalElectrons, atom.radicalElectrons)
            self.assertEqual(atom0.charge, atom.charge)
            self.assertEqual(atom0.label, atom.label)
    
    def testApplyActionIncrementBond(self):
        """
        Test the Atom.applyAction() method for a CHANGE_BOND action.
        """
        action = ['CHANGE_BOND', '*1', 1, '*2']
        for element in elementList:
            atom0 = Atom(element=element, radicalElectrons=1, charge=0, label='*1', lonePairs=0)
            atom = atom0.copy()
            atom.applyAction(action)
            self.assertEqual(atom0.element, atom.element)
            self.assertEqual(atom0.radicalElectrons, atom.radicalElectrons)
            self.assertEqual(atom0.charge, atom.charge)
            self.assertEqual(atom0.label, atom.label)
    
    def testApplyActionDecrementBond(self):
        """
        Test the Atom.applyAction() method for a CHANGE_BOND action.
        """
        action = ['CHANGE_BOND', '*1', -1, '*2']
        for element in elementList:
            atom0 = Atom(element=element, radicalElectrons=1, charge=0, label='*1', lonePairs=0)
            atom = atom0.copy()
            atom.applyAction(action)
            self.assertEqual(atom0.element, atom.element)
            self.assertEqual(atom0.radicalElectrons, atom.radicalElectrons)
            self.assertEqual(atom0.charge, atom.charge)
            self.assertEqual(atom0.label, atom.label)
    
    def testApplyActionGainRadical(self):
        """
        Test the Atom.applyAction() method for a GAIN_RADICAL action.
        """
        action = ['GAIN_RADICAL', '*1', 1]
        for element in elementList:
            atom0 = Atom(element=element, radicalElectrons=1, charge=0, label='*1', lonePairs=0)
            atom = atom0.copy()
            atom.applyAction(action)
            self.assertEqual(atom0.element, atom.element)
            self.assertEqual(atom0.radicalElectrons, atom.radicalElectrons - 1)
            self.assertEqual(atom0.charge, atom.charge)
            self.assertEqual(atom0.label, atom.label)
    
    def testApplyActionLoseRadical(self):
        """
        Test the Atom.applyAction() method for a LOSE_RADICAL action.
        """
        action = ['LOSE_RADICAL', '*1', 1]
        for element in elementList:
            atom0 = Atom(element=element, radicalElectrons=1, charge=0, label='*1', lonePairs=0)
            atom = atom0.copy()
            atom.applyAction(action)
            self.assertEqual(atom0.element, atom.element)
            self.assertEqual(atom0.radicalElectrons, atom.radicalElectrons + 1)
            self.assertEqual(atom0.charge, atom.charge)
            self.assertEqual(atom0.label, atom.label)
    
    def testEquivalent(self):
        """
        Test the Atom.equivalent() method.
        """
        for index1, element1 in enumerate(elementList[0:10]):
            for index2, element2 in enumerate(elementList[0:10]):
                atom1 = Atom(element=element1, radicalElectrons=1, charge=0, label='*1', lonePairs=0)
                atom2 = Atom(element=element2, radicalElectrons=1, charge=0, label='*1', lonePairs=0)
                if index1 == index2:
                    self.assertTrue(atom1.equivalent(atom2))
                    self.assertTrue(atom2.equivalent(atom1))
                else:
                    self.assertFalse(atom1.equivalent(atom2))
                    self.assertFalse(atom2.equivalent(atom1))
    
    def testIsSpecificCaseOf(self):
        """
        Test the Atom.isSpecificCaseOf() method.
        """
        for index1, element1 in enumerate(elementList[0:10]):
            for index2, element2 in enumerate(elementList[0:10]):
                atom1 = Atom(element=element1, radicalElectrons=1, charge=0, label='*1', lonePairs=0)
                atom2 = Atom(element=element2, radicalElectrons=1, charge=0, label='*1', lonePairs=0)
                if index1 == index2:
                    self.assertTrue(atom1.isSpecificCaseOf(atom2))
                else:
                    self.assertFalse(atom1.isSpecificCaseOf(atom2))
    
    def testCopy(self):
        """
        Test the Atom.copy() method.
        """
        atom = self.atom.copy()
        self.assertEqual(self.atom.element.symbol, atom.element.symbol)
        self.assertEqual(self.atom.atomType, atom.atomType)
        self.assertEqual(self.atom.radicalElectrons, atom.radicalElectrons)
        self.assertEqual(self.atom.charge, atom.charge)
        self.assertEqual(self.atom.label, atom.label)
    
    def testPickle(self):
        """
        Test that a Atom object can be successfully pickled and
        unpickled with no loss of information.
        """
        import cPickle
        atom = cPickle.loads(cPickle.dumps(self.atom))
        self.assertEqual(self.atom.element.symbol, atom.element.symbol)
        self.assertEqual(self.atom.atomType, atom.atomType)
        self.assertEqual(self.atom.radicalElectrons, atom.radicalElectrons)
        self.assertEqual(self.atom.charge, atom.charge)
        self.assertEqual(self.atom.label, atom.label)
    
    def testIsotopeEquivalent(self):
        """
        Test the Atom.equivalent() method for non-normal isotopes
        """

        atom1 = Atom(element=getElement('H'))
        atom2 = Atom(element=getElement('H', 2))
        atom3 = Atom(element=getElement('H'))

        self.assertFalse(atom1.equivalent(atom2))
        self.assertTrue(atom1.equivalent(atom3))

    def testGetBondOrdersForAtom(self):
        """
        Test Atom.getBondOrdersForAtom for all carbons in naphthalene
        """
        
        m = Molecule().fromSMILES('C12C(C=CC=C1)=CC=CC=2')
        isomers = m.generate_resonance_structures()
        for isomer in isomers:
            for atom in isomer.atoms:
                if atom.symbol == 'C':
                    self.assertEqual(atom.getBondOrdersForAtom(), 4.0)


################################################################################

class TestBond(unittest.TestCase):
    """
    Contains unit tests of the Bond class.
    """

    def setUp(self):
        """
        A method called before each unit test in this class.
        """
        self.bond = Bond(atom1=None, atom2=None, order=2)
        self.orderList = [1,2,3,4,1.5, 0.30000000000000004]
    
    def testGetOrderStr(self):
        """
        test the Bond.getOrderStr() method
        """
        
        self.assertEqual(self.bond.getOrderStr(),'D')
        
    def testSetOrderStr(self):
        """
        test the Bond.setOrderStr() method
        """
        
        self.bond.setOrderStr("B")
        self.assertEqual(self.bond.order, 1.5)
    
    def testGetOrderNum(self):
        """
        test the Bond.getOrderNum() method
        """
        self.assertEqual(self.bond.getOrderNum(),2)
        
    def testSetOrderNum(self):
        """
        test the Bond.setOrderNum() method
        """
        
        self.bond.setOrderNum(3)
        self.assertEqual(self.bond.getOrderStr(),'T')
    
        
    def testIsOrder(self):
        """
        Test the Bond.isOrder() method.
        """
        for order in self.orderList:
            bond = Bond(None, None, order=order)
            self.assertTrue(bond.isOrder(round(order,2)))
        
        
    def testIsSingle(self):
        """
        Test the Bond.isSingle() method.
        """
        for order in self.orderList:
            bond = Bond(None, None, order=order)
            if order == 1:
                self.assertTrue(bond.isSingle())
            else:
                self.assertFalse(bond.isSingle())
    
    def testIsSingleCanTakeFloatingPointAddition(self):
        """
        Test the Bond.isSingle() method with taking floating point addition
        roundoff errors
        """
        new_order = 0.1 + 0.3*3
        self.assertNotEqual(new_order, 1)
        
        self.bond.setOrderNum(new_order)
        self.assertTrue(self.bond.isSingle())
        
    def testIsDouble(self):
        """
        Test the Bond.isDouble() method.
        """
        for order in self.orderList:
            bond = Bond(None, None, order=order)
            if order == 2:
                self.assertTrue(bond.isDouble())
            else:
                self.assertFalse(bond.isDouble())

    def testIsTriple(self):
        """
        Test the Bond.isTriple() method.
        """
        for order in self.orderList:
            bond = Bond(None, None, order=order)
            if order == 3:
                self.assertTrue(bond.isTriple())
            else:
                self.assertFalse(bond.isTriple())

    def testIsBenzene(self):
        """
        Test the Bond.isBenzene() method.
        """
        for order in self.orderList:
            bond = Bond(None, None, order=order)
            if order == 1.5:
                self.assertTrue(bond.isBenzene())
            else:
                self.assertFalse(bond.isBenzene())

    def testIsQuadruple(self):
        """
        Test the Bond.isQuadruple() method.
        """
        for order in self.orderList:
            bond = Bond(None, None, order=order)
            if order == 4:
                self.assertTrue(bond.isQuadruple())
            else:
                self.assertFalse(bond.isQuadruple())
                
    def testIncrementOrder(self):
        """
        Test the Bond.incrementOrder() method.
        """
        for order in self.orderList:
            bond = Bond(None, None, order=order)
            try:
                bond.incrementOrder()
                if order == 1: 
                    self.assertTrue(bond.isDouble())
                elif order == 2: 
                    self.assertTrue(bond.isTriple())
                elif order == 3:
                    self.assertTrue(bond.isQuadruple())
            except ActionError:
                self.assertTrue(order >= 4) # or benzene??
        
    def testDecrementOrder(self):
        """
        Test the Bond.decrementOrder() method.
        """
        for order in self.orderList:
            bond = Bond(None, None, order=order)
            try:
                bond.decrementOrder()
                if order == 2: 
                    self.assertTrue(bond.isSingle())
                elif order == 3: 
                    self.assertTrue(bond.isDouble())
                elif order == 'Q':
                    self.assertTrue(bond.isTriple())
            except ActionError:
                self.assertTrue(order < 1)
                
    def testApplyActionBreakBond(self):
        """
        Test the Bond.applyAction() method for a BREAK_BOND action.
        """
        action = ['BREAK_BOND', '*1', 1, '*2']
        for order0 in self.orderList:
            bond0 = Bond(None, None, order=order0)
            bond = bond0.copy()
            try:
                bond.applyAction(action)
                self.fail('Bond.applyAction() unexpectedly processed a BREAK_BOND action with order {0}.'.format(order0))
            except ActionError:
                pass
    
    def testApplyActionFormBond(self):
        """
        Test the Bond.applyAction() method for a FORM_BOND action.
        """
        action = ['FORM_BOND', '*1', 1, '*2']
        for order0 in self.orderList:
            bond0 = Bond(None, None, order=order0)
            bond = bond0.copy()
            try:
                bond.applyAction(action)
                self.fail('Bond.applyAction() unexpectedly processed a FORM_BOND action with order {0}.'.format(order0))
            except ActionError:
                pass
    
    def testApplyActionIncrementBond(self):
        """
        Test the Bond.applyAction() method for a CHANGE_BOND action.
        """
        action = ['CHANGE_BOND', '*1', 1, '*2']
        for order0 in self.orderList:
            bond0 = Bond(None, None, order=order0)
            bond = bond0.copy()
            try:
                bond.applyAction(action)
            except ActionError:
                self.assertTrue(4 <= order0,'Test failed with order {0}'.format(order0))
                
    def testApplyActionDecrementBond(self):
        """
        Test the Bond.applyAction() method for a CHANGE_BOND action.
        """
        action = ['CHANGE_BOND', '*1', -1, '*2']
        for order0 in self.orderList:
            bond0 = Bond(None, None, order=order0)
            bond = bond0.copy()
            try:
                bond.applyAction(action)
            except ActionError:
                self.assertTrue(order0 < 1,'Test failed with order {0}'.format(order0))
            
    def testApplyActionGainRadical(self):
        """
        Test the Bond.applyAction() method for a GAIN_RADICAL action.
        """
        action = ['GAIN_RADICAL', '*1', 1]
        for order0 in self.orderList:
            bond0 = Bond(None, None, order=order0)
            bond = bond0.copy()
            try:
                bond.applyAction(action)
                self.fail('Bond.applyAction() unexpectedly processed a GAIN_RADICAL action with order {0}.'.format(order0))
            except ActionError:
                pass
    
    def testApplyActionLoseRadical(self):
        """
        Test the Bond.applyAction() method for a LOSE_RADICAL action.
        """
        action = ['LOSE_RADICAL', '*1', 1]
        for order0 in self.orderList:
            bond0 = Bond(None, None, order=order0)
            bond = bond0.copy()
            try:
                bond.applyAction(action)
                self.fail('Bond.applyAction() unexpectedly processed a LOSE_RADICAL action with order {0}.'.format(order0))
            except ActionError:
                pass
    
    def testEquivalent(self):
        """
        Test the GroupBond.equivalent() method.
        """
        for order1 in self.orderList:
            for order2 in self.orderList:
                bond1 = Bond(None, None, order=order1)
                bond2 = Bond(None, None, order=order2)
                if order1 == order2:
                    self.assertTrue(bond1.equivalent(bond2))
                    self.assertTrue(bond2.equivalent(bond1))
                else:
                    self.assertFalse(bond1.equivalent(bond2))
                    self.assertFalse(bond2.equivalent(bond1))
    
    def testIsSpecificCaseOf(self):
        """
        Test the Bond.isSpecificCaseOf() method.
        """
        for order1 in self.orderList:
            for order2 in self.orderList:
                bond1 = Bond(None, None, order=order1)
                bond2 = Bond(None, None, order=order2)
                if order1 == order2:
                    self.assertTrue(bond1.isSpecificCaseOf(bond2))
                else:
                    self.assertFalse(bond1.isSpecificCaseOf(bond2))
                
    def testCopy(self):
        """
        Test the Bond.copy() method.
        """
        bond = self.bond.copy()
        self.assertEqual(self.bond.order, bond.order)
    
    def testPickle(self):
        """
        Test that a Bond object can be successfully pickled and
        unpickled with no loss of information.
        """
        import cPickle
        bond = cPickle.loads(cPickle.dumps(self.bond))
        self.assertEqual(self.bond.order, bond.order)

    def testUpdateLonePairs(self):
        """
        Test that updateLonePairs works as expected
        """
        mol_N1sc_N5t = Molecule().fromAdjacencyList("""
            1 N u0 p0 c+1 {2,T} {4,S}
            2 N u0 p0 c+1 {1,T} {3,S}
            3 N u0 p3 c-2 {2,S}
            4 H u0 p0 c0 {1,S}""")
        mol_N1s = Molecule().fromAdjacencyList("""
            1 N u0 p2 c0 {2,S}
            2 H u0 p0 c0 {1,S}""")
        mol_N3s = Molecule().fromAdjacencyList("""
            multiplicity 3
            1 N u2 p1 c0 {2,S}
            2 H u0 p0 c0 {1,S}""")
        mol_N3b = Molecule().fromAdjacencyList("""
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
            11 H u0 p0 c0 {6,S}""")
        mol_N5s = Molecule().fromAdjacencyList("""
            multiplicity 2
            1 N u1 p0 c+1 {2,S} {3,S} {4,S}
            2 H u0 p0 c0 {1,S}
            3 H u0 p0 c0 {1,S}
            4 O u0 p3 c-1 {1,S}""")
        mol_N5d = Molecule().fromAdjacencyList("""
            1 N u0 p0 c+1 {2,D} {3,S} {4,S}
            2 O u0 p2 c0 {1,D}
            3 O u0 p2 c0 {1,S} {5,S}
            4 O u0 p3 c-1 {1,S}
            5 H u0 p0 c0 {3,S}""")
        mol_N5dd = Molecule().fromAdjacencyList("""
            1 N u0 p2 c-1 {2,D}
            2 N u0 p0 c+1 {1,D} {3,D}
            3 O u0 p2 c0 {2,D}""")
        mol_CH2_S = Molecule().fromAdjacencyList("""
            1 C u0 p1 c0 {2,S} {3,S}
            2 H u0 p0 c0 {1,S}
            3 H u0 p0 c0 {1,S}""")
        mol_carbonyl = Molecule().fromAdjacencyList("""
            1 O u0 p2 c0 {2,D}
            2 C u0 p0 c0 {1,D} {3,S} {4,S}
            3 H u0 p0 c0 {2,S}
            4 H u0 p0 c0 {2,S}""")

        mol_N1sc_N5t.updateLonePairs()
        mol_N1s.updateLonePairs()
        mol_N3s.updateLonePairs()
        mol_N3b.updateLonePairs()
        mol_N5s.updateLonePairs()
        mol_N5d.updateLonePairs()
        mol_N5dd.updateLonePairs()
        mol_CH2_S.updateLonePairs()
        mol_carbonyl.updateLonePairs()

        self.assertEqual(mol_N1sc_N5t.atoms[0].lonePairs, 0)
        self.assertEqual(mol_N1sc_N5t.atoms[2].lonePairs, 3)
        self.assertEqual(mol_N1s.atoms[0].lonePairs, 2)
        self.assertEqual(mol_N3s.atoms[0].lonePairs, 1)
        self.assertEqual(mol_N3b.atoms[0].lonePairs, 1)
        self.assertEqual(mol_N5s.atoms[0].lonePairs, 0)
        self.assertEqual(mol_N5s.atoms[3].lonePairs, 3)
        self.assertEqual(mol_N5d.atoms[0].lonePairs, 0)
        self.assertEqual(mol_N5d.atoms[1].lonePairs, 2)
        self.assertEqual(mol_N5d.atoms[2].lonePairs, 2)
        self.assertEqual(mol_N5d.atoms[3].lonePairs, 3)
        self.assertEqual(mol_N5dd.atoms[0].lonePairs, 2)
        self.assertEqual(mol_N5dd.atoms[1].lonePairs, 0)
        self.assertEqual(mol_N5dd.atoms[2].lonePairs, 2)
        self.assertEqual(mol_CH2_S.atoms[0].lonePairs, 1)
        self.assertEqual(mol_carbonyl.atoms[0].lonePairs, 2)
        self.assertEqual(mol_carbonyl.atoms[1].lonePairs, 0)

    def test_get_bond_string(self):
        """Test that bond objects can return a bond string"""
        bond = Bond(atom1=Atom(element=getElement(1)), atom2=Atom(element=getElement(6)), order=1)
        self.assertEqual(bond.get_bond_string(), 'C-H')
        
################################################################################

class TestMolecule(unittest.TestCase):
    """
    Contains unit tests of the Molecule class.
    """
    
    def setUp(self):
        self.adjlist_1 = """
1 *1 C u1 p0 c0  {2,S} {3,S} {4,S}
2    H u0 p0 c0  {1,S}
3    H u0 p0 c0  {1,S}
4 *2 N u0 p0 c+1 {1,S} {5,S} {6,D}
5    O u0 p3 c-1 {4,S}
6    O u0 p2 c0  {4,D}
            """
        self.molecule = [Molecule().fromAdjacencyList(self.adjlist_1)]
        
        self.adjlist_2 = """
1 *1 C u1 p0 {2,S}
2 *2 N u0 p0 c+1 {1,S} {3,S} {4,D}
3    O u0 p3 c-1 {2,S}
4    O u0 p2 {2,D}
            """
        self.molecule.append(Molecule().fromAdjacencyList(self.adjlist_2,saturateH=True))
        
        
        self.mHBonds = Molecule().fromSMILES('C(NC=O)OO')
        
    def testClearLabeledAtoms(self):
        """
        Test the Molecule.clearLabeledAtoms() method.
        """
        self.molecule[0].clearLabeledAtoms()
        for atom in self.molecule[0].atoms:
            self.assertEqual(atom.label, '')

    def testContainsLabeledAtom(self):
        """
        Test the Molecule.containsLabeledAtom() method.
        """
        for atom in self.molecule[0].atoms:
            if atom.label != '':
                self.assertTrue(self.molecule[0].containsLabeledAtom(atom.label))
        self.assertFalse(self.molecule[0].containsLabeledAtom('*3'))
        self.assertFalse(self.molecule[0].containsLabeledAtom('*4'))
        self.assertFalse(self.molecule[0].containsLabeledAtom('*5'))
        self.assertFalse(self.molecule[0].containsLabeledAtom('*6'))
        
    def testGetLabeledAtom(self):
        """
        Test the Molecule.getLabeledAtom() method.
        """
        for atom in self.molecule[0].atoms:
            if atom.label != '':
                self.assertEqual(atom, self.molecule[0].getLabeledAtom(atom.label))
        try:
            self.molecule[0].getLabeledAtom('*3')
            self.fail('Unexpected successful return from Molecule.getLabeledAtom() with invalid atom label.')
        except ValueError:
            pass
            
    def testGetLabeledAtoms(self):
        """
        Test the Molecule.getLabeledAtoms() method.
        """
        labeled = self.molecule[0].getLabeledAtoms()
        for atom in self.molecule[0].atoms:
            if atom.label != '':
                self.assertTrue(atom.label in labeled)
                self.assertTrue(atom in labeled.values())
            else:
                self.assertFalse(atom.label in labeled)
                self.assertFalse(atom in labeled.values())
        
        multipleLabelMolecule = Molecule().fromAdjacencyList("""
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
""")
        labeled = multipleLabelMolecule.getLabeledAtoms()
        self.assertTrue('*' in labeled)
        self.assertTrue('*1' in labeled)
        self.assertEqual(len(labeled['*']),4)
        self.assertEqual(len(labeled['*1']),3)
        
    def testGetFormula(self):
        """
        Test the Molecule.getLabeledAtoms() method.
        """
        self.assertEqual(self.molecule[0].getFormula(), 'CH2NO2')
        self.assertEqual(self.molecule[1].getFormula(), 'CH2NO2')


    def testRadicalCount(self):
        """
        Test the Molecule.getRadicalCount() method.
        """
        self.assertEqual( self.molecule[0].getRadicalCount(), sum([atom.radicalElectrons for atom in self.molecule[0].atoms]) )
        self.assertEqual( self.molecule[1].getRadicalCount(), sum([atom.radicalElectrons for atom in self.molecule[1].atoms]) )
        
    def testGetMolecularWeight(self):
        """
        Test the Molecule.getMolecularWeight() method.
        """
        self.assertAlmostEqual(self.molecule[0].getMolecularWeight() * 1000, 60.03, 2)
        self.assertAlmostEqual(self.molecule[1].getMolecularWeight() * 1000, 60.03, 2)

    def testFromAdjacencyList(self):
        """
        Test the Molecule.fromAdjacencyList() method.
        """
        
        # molecule 1
        
        self.assertTrue(self.molecule[0].multiplicity == 2)
        
        atom1 = self.molecule[0].atoms[0]
        atom2 = self.molecule[0].atoms[3]
        atom3 = self.molecule[0].atoms[4]
        atom4 = self.molecule[0].atoms[5]
        self.assertTrue(self.molecule[0].hasBond(atom2,atom1))
        self.assertTrue(self.molecule[0].hasBond(atom2,atom3))
        self.assertTrue(self.molecule[0].hasBond(atom2,atom4))
        self.assertFalse(self.molecule[0].hasBond(atom1,atom3))
        self.assertFalse(self.molecule[0].hasBond(atom1,atom4))
        bond21 = atom2.bonds[atom1]
        bond23 = atom2.bonds[atom3]
        bond24 = atom2.bonds[atom4]
           
        self.assertTrue(atom1.label == '*1')
        self.assertTrue(atom1.element.symbol == 'C')
        self.assertTrue(atom1.radicalElectrons == 1)
        self.assertTrue(atom1.charge == 0)
        
        self.assertTrue(atom2.label == '*2')
        self.assertTrue(atom2.element.symbol == 'N')
        self.assertTrue(atom2.radicalElectrons == 0)
        self.assertTrue(atom2.charge == 1)
        
        self.assertTrue(atom3.label == '')
        self.assertTrue(atom3.element.symbol == 'O')
        self.assertTrue(atom3.radicalElectrons == 0)
        self.assertTrue(atom3.charge == -1)
        
        self.assertTrue(atom4.label == '')
        self.assertTrue(atom4.element.symbol == 'O')
        self.assertTrue(atom4.radicalElectrons == 0)
        self.assertTrue(atom4.charge == 0)

        self.assertTrue(bond21.isSingle())
        self.assertTrue(bond23.isSingle())
        self.assertTrue(bond24.isDouble())
        
        # molecule 2
        
        self.assertTrue(self.molecule[1].multiplicity == 2)
        
        atom1 = self.molecule[1].atoms[0]
        atom2 = self.molecule[1].atoms[1]
        atom3 = self.molecule[1].atoms[2]
        atom4 = self.molecule[1].atoms[3]
        self.assertTrue(self.molecule[1].hasBond(atom2,atom1))
        self.assertTrue(self.molecule[1].hasBond(atom2,atom3))
        self.assertTrue(self.molecule[1].hasBond(atom2,atom4))
        self.assertFalse(self.molecule[1].hasBond(atom1,atom3))
        self.assertFalse(self.molecule[1].hasBond(atom1,atom4))
        bond21 = atom2.bonds[atom1]
        bond23 = atom2.bonds[atom3]
        bond24 = atom2.bonds[atom4]
           
        self.assertTrue(atom1.label == '*1')
        self.assertTrue(atom1.element.symbol == 'C')
        self.assertTrue(atom1.radicalElectrons == 1)
        self.assertTrue(atom1.charge == 0)
        
        self.assertTrue(atom2.label == '*2')
        self.assertTrue(atom2.element.symbol == 'N')
        self.assertTrue(atom2.radicalElectrons == 0)
        self.assertTrue(atom2.charge == 1)
        
        self.assertTrue(atom3.label == '')
        self.assertTrue(atom3.element.symbol == 'O')
        self.assertTrue(atom3.radicalElectrons == 0)
        self.assertTrue(atom3.charge == -1)
        
        self.assertTrue(atom4.label == '')
        self.assertTrue(atom4.element.symbol == 'O')
        self.assertTrue(atom4.radicalElectrons == 0)
        self.assertTrue(atom4.charge == 0)

        self.assertTrue(bond21.isSingle())
        self.assertTrue(bond23.isSingle())
        self.assertTrue(bond24.isDouble())
        
        

    def testToAdjacencyList(self):
        """
        Test the Molecule.toAdjacencyList() method.
        """
        adjlist_1 = self.molecule[0].toAdjacencyList(removeH=False)
        newMolecule = Molecule().fromAdjacencyList(adjlist_1)
        self.assertTrue(self.molecule[0].isIsomorphic(newMolecule))
        
        #self.assertEqual(adjlist_1.strip(), self.adjlist_1.strip())
        
#    def testFromOldAdjacencyList(self):
#        """
#        Test we can read things with implicit hydrogens.
#        """
#        adjList = """
#        1 O 0 
#        """ # should be Water
#        molecule = Molecule().fromAdjacencyList(adjList, saturateH=True) # only works with saturateH=True
#        self.assertEqual(molecule.getFormula(),'H2O')

    def testIsomorphism(self):
        """
        Check the graph isomorphism functions.
        """
        molecule1 = Molecule().fromSMILES('C=CC=C[CH]C')
        molecule2 = Molecule().fromSMILES('C[CH]C=CC=C')
        self.assertTrue(molecule1.isIsomorphic(molecule2))
        self.assertTrue(molecule2.isIsomorphic(molecule1))
        
        molecule1 = Molecule().fromAdjacencyList("""
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
15 H u0 p0 c0 {7,S}""")
        molecule2 = Molecule().fromAdjacencyList("""
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
15 H u0 p0 c0 {7,S}""")
        
        self.assertTrue(molecule1.isIsomorphic(molecule2,generateInitialMap=True))
        self.assertTrue(molecule2.isIsomorphic(molecule1,generateInitialMap=True))
    
        
    def testSubgraphIsomorphism(self):
        """
        Check the graph isomorphism functions.
        """
        molecule = Molecule().fromSMILES('C=CC=C[CH]C')
        group = Group().fromAdjacencyList("""
        1 Cd u0 p0 c0 {2,D}
        2 Cd u0 p0 c0 {1,D}
        """)

        self.assertTrue(molecule.isSubgraphIsomorphic(group))
        mapping = molecule.findSubgraphIsomorphisms(group)
        self.assertTrue(len(mapping) == 4, "len(mapping) = %d, should be = 4" % (len(mapping)))
        for map in mapping:
            self.assertTrue(len(map) == min(len(molecule.atoms), len(group.atoms)))
            for key, value in map.iteritems():
                self.assertTrue(key in molecule.atoms)
                self.assertTrue(value in group.atoms)

    def testSubgraphIsomorphismAgain(self):
        molecule = Molecule()
        molecule.fromAdjacencyList("""
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
        """)

        group = Group()
        group.fromAdjacencyList("""
        1 * C u0 p0 c0 {2,D} {3,S} {4,S}
        2   C u0 p0 c0 {1,D}
        3   H u0 p0 c0 {1,S}
        4   H u0 p0 c0 {1,S}
        """)

        labeled1 = molecule.getLabeledAtoms().values()[0]
        labeled2 = group.getLabeledAtoms().values()[0]

        initialMap = {labeled1: labeled2}
        self.assertTrue(molecule.isSubgraphIsomorphic(group, initialMap))

        initialMap = {labeled1: labeled2}
        mapping = molecule.findSubgraphIsomorphisms(group, initialMap)
        self.assertTrue(len(mapping) == 2,  "len(mapping) = %d, should be = 2" % (len(mapping)))
        for map in mapping:
            self.assertTrue(len(map) == min(len(molecule.atoms), len(group.atoms)))
            for key, value in map.iteritems():
                self.assertTrue(key in molecule.atoms)
                self.assertTrue(value in group.atoms)

    def testSubgraphIsomorphismManyLabels(self):
        molecule = Molecule() # specific case (species)
        molecule.fromAdjacencyList("""
1 *1 C  u1 p0 c0 {2,S} {3,S} {4,S}
2    C  u0 p0 c0 {1,S} {3,S} {5,S} {6,S}
3    C  u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
4    H  u0 p0 c0 {1,S}
5    H  u0 p0 c0 {2,S}
6    H  u0 p0 c0 {2,S}
7    H  u0 p0 c0 {3,S}
8    H  u0 p0 c0 {3,S}
        """)

        group = Group() # general case (functional group)
        group.fromAdjacencyList("""
1 *1 C   u1 p0 c0 {2,S}, {3,S}
2    R!H u0 p0 c0 {1,S}
3    R!H u0 p0 c0 {1,S}
        """)

        labeled1 = molecule.getLabeledAtoms()
        labeled2 = group.getLabeledAtoms()
        initialMap = {}
        for label,atom1 in labeled1.iteritems():
            initialMap[atom1] = labeled2[label]
        self.assertTrue(molecule.isSubgraphIsomorphic(group, initialMap))

        mapping = molecule.findSubgraphIsomorphisms(group, initialMap)
        self.assertEqual(len(mapping), 2)
        for map in mapping:
            self.assertTrue(len(map) == min(len(molecule.atoms), len(group.atoms)))
            for key, value in map.iteritems():
                self.assertTrue(key in molecule.atoms)
                self.assertTrue(value in group.atoms)

    def testSubgraphIsomorphismRings(self):
        molecule = Molecule(SMILES='C1CCCC1CCC')
        groupNoRing = Group().fromAdjacencyList("""
1 *1 C u0 p0 c0 r0
        """)
        groupRing = Group().fromAdjacencyList("""
1 *1 C u0 p0 c0 r1
        """)

        self.assertTrue(molecule.isSubgraphIsomorphic(groupNoRing))
        mapping = molecule.findSubgraphIsomorphisms(groupNoRing)
        self.assertEqual(len(mapping), 3)
        self.assertTrue(molecule.isSubgraphIsomorphic(groupRing))
        mapping = molecule.findSubgraphIsomorphisms(groupRing)
        self.assertEqual(len(mapping), 5)

    def test_lax_isomorphism(self):
        """Test that we can do isomorphism comparison with strict=False"""
        mol1 = Molecule().fromAdjacencyList("""
multiplicity 2
1 O u0 p2 c0 {3,D}
2 C u1 p0 c0 {3,S} {4,S} {5,S}
3 C u0 p0 c0 {1,D} {2,S} {6,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
        """)

        mol2 = Molecule().fromAdjacencyList("""
multiplicity 2
1 O u1 p2 c0 {3,S}
2 C u0 p0 c0 {3,D} {4,S} {5,S}
3 C u0 p0 c0 {1,S} {2,D} {6,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
        """)

        self.assertTrue(mol1.isIsomorphic(mol2, strict=False))

    def testAdjacencyList(self):
        """
        Check the adjacency list read/write functions for a full molecule.
        """
        molecule1 = Molecule().fromAdjacencyList("""
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
        """)
        molecule2 = Molecule().fromSMILES('C=CC=C[CH]C')
        self.assertTrue(molecule1.isIsomorphic(molecule2))
        self.assertTrue(molecule2.isIsomorphic(molecule1))
    
    def test_generate_H_bonded_structures(self):
        """
        Test that the correct set of Hydrogen Bonded structures are generated
        """
        correctSet = ['1  C u0 p0 c0 {2,S} {4,S} {6,S} {7,S}\n2  N u0 p1 c0 {1,S} {3,S} {8,S}\n3  C u0 p0 c0 {2,S} {9,D} {10,S}\n4  O u0 p2 c0 {1,S} {5,S}\n5  O u0 p2 c0 {4,S} {8,H} {11,S}\n6  H u0 p0 c0 {1,S}\n7  H u0 p0 c0 {1,S}\n8  H u0 p0 c0 {2,S} {5,H}\n9  O u0 p2 c0 {3,D}\n10 H u0 p0 c0 {3,S}\n11 H u0 p0 c0 {5,S}\n',
 '1  C u0 p0 c0 {2,S} {4,S} {6,S} {7,S}\n2  N u0 p1 c0 {1,S} {3,S} {8,S} {11,H}\n3  C u0 p0 c0 {2,S} {9,D} {10,S}\n4  O u0 p2 c0 {1,S} {5,S}\n5  O u0 p2 c0 {4,S} {11,S}\n6  H u0 p0 c0 {1,S}\n7  H u0 p0 c0 {1,S}\n8  H u0 p0 c0 {2,S}\n9  O u0 p2 c0 {3,D}\n10 H u0 p0 c0 {3,S}\n11 H u0 p0 c0 {2,H} {5,S}\n',
 '1  C u0 p0 c0 {2,S} {4,S} {6,S} {7,S}\n2  N u0 p1 c0 {1,S} {3,S} {8,S} {11,H}\n3  C u0 p0 c0 {2,S} {9,D} {10,S}\n4  O u0 p2 c0 {1,S} {5,S}\n5  O u0 p2 c0 {4,S} {8,H} {11,S}\n6  H u0 p0 c0 {1,S}\n7  H u0 p0 c0 {1,S}\n8  H u0 p0 c0 {2,S} {5,H}\n9  O u0 p2 c0 {3,D}\n10 H u0 p0 c0 {3,S}\n11 H u0 p0 c0 {2,H} {5,S}\n',
 '1  C u0 p0 c0 {2,S} {4,S} {6,S} {7,S}\n2  N u0 p1 c0 {1,S} {3,S} {8,S}\n3  C u0 p0 c0 {2,S} {9,D} {10,S}\n4  O u0 p2 c0 {1,S} {5,S}\n5  O u0 p2 c0 {4,S} {11,S}\n6  H u0 p0 c0 {1,S}\n7  H u0 p0 c0 {1,S}\n8  H u0 p0 c0 {2,S}\n9  O u0 p2 c0 {3,D} {11,H}\n10 H u0 p0 c0 {3,S}\n11 H u0 p0 c0 {5,S} {9,H}\n',
 '1  C u0 p0 c0 {2,S} {4,S} {6,S} {7,S}\n2  N u0 p1 c0 {1,S} {3,S} {8,S}\n3  C u0 p0 c0 {2,S} {9,D} {10,S}\n4  O u0 p2 c0 {1,S} {5,S}\n5  O u0 p2 c0 {4,S} {8,H} {11,S}\n6  H u0 p0 c0 {1,S}\n7  H u0 p0 c0 {1,S}\n8  H u0 p0 c0 {2,S} {5,H}\n9  O u0 p2 c0 {3,D} {11,H}\n10 H u0 p0 c0 {3,S}\n11 H u0 p0 c0 {5,S} {9,H}\n']
        
        mols = [Molecule().fromAdjacencyList(k) for k in correctSet]
        
        self.assertEqual(set(mols),set(self.mHBonds.generate_H_bonded_structures()))
    
    def test_remove_H_bonds(self):
        """
        test that remove HBonds removes all hydrogen bonds from a given molecule
        """
        testMol = self.mHBonds.generate_H_bonded_structures()[0]
        testMol.remove_H_bonds()
        
        for i,atm1 in enumerate(testMol.atoms):
            for j,atm2 in enumerate(testMol.atoms):
                if j<i and testMol.hasBond(atm1,atm2):
                    bd = testMol.getBond(atm1,atm2)
                    self.assertNotAlmostEqual(bd.order,0.1)
                    
    def testSSSR(self):
        """
        Test the Molecule.getSmallestSetOfSmallestRings() method with a complex
        polycyclic molecule.
        """
        molecule = Molecule()
        molecule.fromSMILES('C(CC1C(C(CCCCCCCC)C1c1ccccc1)c1ccccc1)CCCCCC')
        #http://cactus.nci.nih.gov/chemical/structure/C(CC1C(C(CCCCCCCC)C1c1ccccc1)c1ccccc1)CCCCCC/image
        sssr = molecule.getSmallestSetOfSmallestRings()
        self.assertEqual( len(sssr), 3)

    def testIsInCycleEthane(self):
        """
        Test the Molecule.isInCycle() method with ethane.
        """
        molecule = Molecule().fromSMILES('CC')
        for atom in molecule.atoms:
            self.assertFalse(molecule.isAtomInCycle(atom))
        for atom1 in molecule.atoms:
            for atom2, bond in atom1.bonds.items():
                self.assertFalse(molecule.isBondInCycle(bond))

    def testIsInCycleCyclohexane(self):
        """
        Test the Molecule.isInCycle() method with ethane.
        """
        molecule = Molecule().fromInChI('InChI=1/C6H12/c1-2-4-6-5-3-1/h1-6H2')
        for atom in molecule.atoms:
            if atom.isHydrogen():
                self.assertFalse(molecule.isAtomInCycle(atom))
            elif atom.isCarbon():
                self.assertTrue(molecule.isAtomInCycle(atom))
        for atom1 in molecule.atoms:
            for atom2, bond in atom1.bonds.items():
                if atom1.isCarbon() and atom2.isCarbon():
                    self.assertTrue(molecule.isBondInCycle(bond))
                else:
                    self.assertFalse(molecule.isBondInCycle(bond))
        
    def testFromSMILESH(self):
        """
        Make sure that H radical is produced properly from its SMILES
        representation.
        """
        molecule = Molecule(SMILES='[H]')
        self.assertEqual(len(molecule.atoms), 1)
        H = molecule.atoms[0]
        self.assertTrue(H.isHydrogen())
        self.assertEqual(H.radicalElectrons, 1)

    def testFromInChIH(self):
        """
        Make sure that H radical is produced properly from its InChI
        representation.
        """
        molecule = Molecule().fromInChI('InChI=1/H')
        self.assertEqual(len(molecule.atoms), 1)
        H = molecule.atoms[0]
        self.assertTrue(H.isHydrogen())
        self.assertEqual(H.radicalElectrons, 1)

    def testPickle(self):
        """
        Test that a Molecule object can be successfully pickled and
        unpickled with no loss of information.
        """
        molecule0 = Molecule().fromSMILES('C=CC=C[CH2]C')
        molecule0.update()
        import cPickle
        molecule = cPickle.loads(cPickle.dumps(molecule0))
        
        self.assertEqual(len(molecule0.atoms), len(molecule.atoms))
        self.assertEqual(molecule0.getFormula(), molecule.getFormula())
        self.assertTrue(molecule0.isIsomorphic(molecule))
        self.assertTrue(molecule.isIsomorphic(molecule0))

    def testRadicalCH(self):
        """
        Test that the species [CH] has one radical electrons and a spin multiplicity of 2.
        """
        molecule = Molecule().fromSMILES('[CH]')
        self.assertEqual(molecule.atoms[0].radicalElectrons, 1)
        self.assertEqual(molecule.multiplicity, 2)
        self.assertEqual(molecule.getRadicalCount(), 1)

    def testRadicalCH2(self):
        """
        Test that the species [CH2] has two radical electrons and a spin multiplicity of 3.
        """
        molecule = Molecule().fromSMILES('[CH2]')
        self.assertEqual(molecule.atoms[0].radicalElectrons, 2)
        self.assertEqual(molecule.multiplicity, 3)
        self.assertEqual(molecule.getRadicalCount(), 2)
        
    def testRadicalCH2CH2CH2(self):
        """
        Test radical count on [CH2]C[CH2]
        """
        molecule = Molecule().fromSMILES('[CH2]C[CH2]')
        self.assertEqual(molecule.getRadicalCount(), 2)

    def testSingletCarbene(self):
        """Test radical and carbene count on singlet carbene."""
        mol = Molecule().fromAdjacencyList("""
1 C u0 p1 {2,S}
2 C u0 p1 {1,S}
""", saturateH=True)
        self.assertEqual(mol.getRadicalCount(), 0)
        self.assertEqual(mol.getSingletCarbeneCount(), 2)

    def testTripletCarbene(self):
        """Test radical and carbene count on triplet carbene."""
        mol = Molecule().fromAdjacencyList("""
1 C u2 p0 {2,S}
2 C u0 p1 {1,S}
""", saturateH=True)
        self.assertEqual(mol.getRadicalCount(), 2)
        self.assertEqual(mol.getSingletCarbeneCount(), 1)

    def testSingletCarbon(self):
        """Test that getSingletCarbeneCount returns 1 for singlet carbon atom."""
        mol = Molecule().fromAdjacencyList('1 C u0 p2')
        self.assertEqual(mol.getSingletCarbeneCount(), 1)
        
    def testSMILES(self):
        """
        Test that we can generate a few SMILES strings as expected
        """
        import rmgpy.molecule
        test_strings = ['[C-]#[O+]', '[C]', '[CH]', 'OO', '[H][H]', '[H]',
                       '[He]', '[O]', 'O', '[CH3]', 'C', '[OH]', 'CCC',
                       'CC', 'N#N', '[O]O', 'C[CH2]', '[Ar]', 'CCCC',
                       'O=C=O', '[C]#N',
                       ]
        for s in test_strings:
            molecule = Molecule(SMILES=s)
            self.assertEqual(s, molecule.toSMILES())

    def testKekuleToSMILES(self):
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
        for smiles, adjlist in test_cases.iteritems():
            m = Molecule().fromAdjacencyList(adjlist)
            s = m.toSMILES()
            self.assertEqual(s, smiles, "Generated SMILES string {0} instead of {1}".format(s, smiles))

    def testKekuleRoundTripSMILES(self):
        """
        Test that we can round-trip SMILES strings of Kekulized aromatics
        """
        test_strings = [
            'CC1=CC=CC=C1O',
            'CC1=C(O)C=CC=C1',
            # 'Cc1ccccc1O',  # this will fail because it is Kekulized during fromSMILES()
        ]
        for s in test_strings:
            molecule = Molecule(SMILES=s)
            self.assertEqual(s, molecule.toSMILES(), "Started with {0} but ended with {1}".format(s, molecule.toSMILES()))

    def testInChIKey(self):
        """
        Test that InChI Key generation is working properly.
        """
        molecule = Molecule().fromInChI('InChI=1S/C7H12/c1-2-7-4-3-6(1)5-7/h6-7H,1-5H2')
        key = molecule.toInChIKey()
        self.assertEqual(key, 'UMRZSTCPUPJPOJ-UHFFFAOYSA-N')
        
    def testAugmentedInChI(self):
        """
        Test the Augmented InChI generation
        """
        mol = Molecule().fromAdjacencyList("""
            1     C     u1 p0 c0 {2,S}
            2     C     u1 p0 c0 {1,S}
        """, saturateH=True)
        
        self.assertEqual(mol.toAugmentedInChI(), 'InChI=1S/C2H4/c1-2/h1-2H2/u1,2')
        
    def testAugmentedInChIKey(self):
        """
        Test the Augmented InChI Key generation
        """
        mol = Molecule().fromAdjacencyList("""
            1     C     u1 p0 c0 {2,S}
            2     C     u1 p0 c0 {1,S}
        """, saturateH=True)
        
        self.assertEqual(mol.toAugmentedInChIKey(), 'VGGSQFUCUMXWEO-UHFFFAOYSA-N-u1,2')

    def testLinearMethane(self):
        """
        Test the Molecule.isLinear() method.
        """
        self.assertFalse(Molecule().fromSMILES('C').isLinear())
    
    def testLinearEthane(self):
        """
        Test the Molecule.isLinear() method.
        """
        self.assertFalse(Molecule().fromSMILES('CC').isLinear())
    
    def testLinearPropane(self):
        """
        Test the Molecule.isLinear() method.
        """
        self.assertFalse(Molecule().fromSMILES('CCC').isLinear())
    
    def testLinearNeopentane(self):
        """
        Test the Molecule.isLinear() method.
        """
        self.assertFalse(Molecule().fromSMILES('CC(C)(C)C').isLinear())
    
    def testLinearHydrogen(self):
        """
        Test the Molecule.isLinear() method.
        """
        self.assertFalse(Molecule().fromSMILES('[H]').isLinear())
    
    def testLinearOxygen(self):
        """
        Test the Molecule.isLinear() method.
        """
        self.assertTrue(Molecule().fromSMILES('O=O').isLinear())
    
    def testLinearCarbonDioxide(self):
        """
        Test the Molecule.isLinear() method.
        """
        self.assertTrue(Molecule().fromSMILES('O=C=O').isLinear())
    
    def testLinearAcetylene(self):
        """
        Test the Molecule.isLinear() method.
        """
        self.assertTrue(Molecule().fromSMILES('C#C').isLinear())
    
    def testLinear135Hexatriyne(self):
        """
        Test the Molecule.isLinear() method.
        """
        self.assertTrue(Molecule().fromSMILES('C#CC#CC#C').isLinear())

    def testAromaticBenzene(self):
        """
        Test the Molecule.isAromatic() method for Benzene.
        """
        m = Molecule().fromSMILES('C1=CC=CC=C1')
        isomers = m.generate_resonance_structures()
        self.assertTrue(any(isomer.isAromatic() for isomer in isomers))

    def testAromaticNaphthalene(self):
        """
        Test the Molecule.isAromatic() method for Naphthalene.
        """
        m = Molecule().fromSMILES('C12C(C=CC=C1)=CC=CC=2')
        isomers = m.generate_resonance_structures()
        self.assertTrue(any(isomer.isAromatic() for isomer in isomers))
                        
    def testAromaticCyclohexane(self):
        """
        Test the Molecule.isAromatic() method for Cyclohexane.
        """
        m = Molecule().fromSMILES('C1CCCCC1')
        isomers = m.generate_resonance_structures()
        self.assertFalse(any(isomer.isAromatic() for isomer in isomers))

    def testHeterocyclicCyclohexanol(self):
        """
        Test the Molecule.isHeterocyclic() method for Cyclohexanol.
        """
        self.assertFalse(Molecule().fromSMILES('OC1CCCCC1').isHeterocyclic())

    def testHeterocyclicFuran(self):
        """
        Test the Molecule.isHeterocyclic() method for Furan.
        """
        self.assertTrue(Molecule().fromSMILES('C1C=COC=1').isHeterocyclic())

    def testHeterocyclicPyridine(self):
        """
        Test the Molecule.isHeterocyclic() method for Pyridine.
        """
        self.assertTrue(Molecule().fromSMILES('c1cccnc1').isHeterocyclic())

    def testCountInternalRotorsEthane(self):
        """
        Test the Molecule.countInternalRotors() method.
        """
        self.assertEqual(Molecule().fromSMILES('CC').countInternalRotors(), 1)
        
    def testCountInternalRotorsPropane(self):
        """
        Test the Molecule.countInternalRotors() method.
        """
        self.assertEqual(Molecule().fromSMILES('CCC').countInternalRotors(), 2)
    
    def testCountInternalRotorsNeopentane(self):
        """
        Test the Molecule.countInternalRotors() method.
        """
        self.assertEqual(Molecule().fromSMILES('CC(C)(C)C').countInternalRotors(), 4)
    
    def testCountInternalRotorsMethylCyclohexane(self):
        """
        Test the Molecule.countInternalRotors() method.
        """
        self.assertEqual(Molecule().fromSMILES('C1CCCC1C').countInternalRotors(), 1)
    
    def testCountInternalRotorsEthylene(self):
        """
        Test the Molecule.countInternalRotors() method.
        """
        self.assertEqual(Molecule().fromSMILES('C=C').countInternalRotors(), 0)
    
    def testCountInternalRotorsAcetylene(self):
        """
        Test the Molecule.countInternalRotors() method.
        """
        self.assertEqual(Molecule().fromSMILES('C#C').countInternalRotors(), 0)
    
    def testCarbeneIdentifiers(self):
        """
        Test that singlet carbene molecules, bearing an electron pair rather than unpaired electrons
        are correctly converted into rdkit molecules and identifiers.
        """
        

        ch2_t = '''
        multiplicity 3
        1 C u2 p0 c0 {2,S} {3,S}
        2 H u0 p0 c0 {1,S}
        3 H u0 p0 c0 {1,S}
        '''
        
        mol = Molecule().fromAdjacencyList(ch2_t)
    
        self.assertEqual( mol.toAugmentedInChI(), 'InChI=1S/CH2/h1H2/u1,1')
        self.assertEqual( mol.toSMILES(), '[CH2]')
        

        ch2_s = '''
        multiplicity 1
        1 C u0 p1 c0 {2,S} {3,S}
        2 H u0 p0 c0 {1,S}
        3 H u0 p0 c0 {1,S}
        '''
        
        mol = Molecule().fromAdjacencyList(ch2_s)
        self.assertEqual( mol.toAugmentedInChI(), 'InChI=1S/CH2/h1H2/lp1')
        self.assertEqual( mol.toSMILES(), '[CH2]')
        
        
    def testGetSymmetryNumber(self):
        """
        Test that the symmetry number getter works properly
        """
        
        mol = Molecule().fromSMILES('C')
        
        self.assertEquals(12, mol.getSymmetryNumber())
        
        empty = Molecule()
        self.assertEquals(1, empty.getSymmetryNumber())
    
    def testMoleculeProps(self):
        """
        Test a key-value pair is added to the props attribute of Molecule.
        """
        self.molecule[0].props['foo'] = 'bar'
        self.assertIsInstance(self.molecule[0].props, dict)
        self.assertEquals(self.molecule[0].props['foo'], 'bar')
        
    def testMoleculeProps_object_attribute(self):
        """
        Test that Molecule's props dictionaries are independent of each other.
        
        Create a test in which is checked whether props is an object attribute rather
        than a class attribute
        """
        spc2 = Molecule()
        self.molecule[0].props['foo'] = 'bar'
        spc3 = Molecule()
        spc3.props['foo'] = 'bla'
        self.assertEquals(self.molecule[0].props['foo'], 'bar')
        self.assertDictEqual(spc2.props, {})
        self.assertDictEqual(spc3.props, {'foo': 'bla'})
        
    @work_in_progress
    def testCountInternalRotorsDimethylAcetylene(self):
        """
        Test the Molecule.countInternalRotors() method for dimethylacetylene.
        
        This is a "hard" test that currently fails.
        """
        self.assertEqual(Molecule().fromSMILES('CC#CC').countInternalRotors(), 1)
        
    def testSaturateAromaticRadical(self):
        """
        Test that the Molecule.saturate() method works properly for an indenyl radical
        containing Benzene bonds
        """
        indenyl = Molecule().fromAdjacencyList("""
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
""")
        indene = Molecule().fromAdjacencyList("""
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
""")
        saturated_molecule = indenyl.copy(deep=True)
        saturated_molecule.saturate_radicals()
        self.assertTrue(saturated_molecule.isIsomorphic(indene))

    def testSurfaceMolecules(self):
        """
        Test that we can identify surface molecules.
        """
        adsorbed = Molecule().fromAdjacencyList("""
                                                1 H u0 p0 c0 {2,S}
                                                2 X u0 p0 c0 {1,S}
                                                """)
        self.assertTrue(adsorbed.containsSurfaceSite())
        gas = Molecule().fromAdjacencyList("""
                                        1 H u0 p0 c0 {2,S}
                                        2 H u0 p0 c0 {1,S}
                                        """)
        self.assertFalse(gas.containsSurfaceSite())

        surface_site = Molecule().fromAdjacencyList("""
                                                1 X u0 p0 c0
                                                """)
        self.assertTrue((surface_site.isSurfaceSite()))
        self.assertFalse((adsorbed.isSurfaceSite()))
        self.assertFalse((gas.isSurfaceSite()))

    def testMalformedAugmentedInChI(self):
        """Test that augmented inchi without InChI layer raises Exception."""
        from .inchi import InchiException

        malform_aug_inchi = 'foo'
        with self.assertRaises(InchiException):
            mol = Molecule().fromAugmentedInChI(malform_aug_inchi)

    def testMalformedAugmentedInChI_Wrong_InChI_Layer(self):
        """Test that augmented inchi with wrong layer is caught."""
        malform_aug_inchi = 'InChI=1S/CH3/h1H2'
        with self.assertRaises(Exception):
            mol = Molecule().fromAugmentedInChI(malform_aug_inchi)

    def testMalformedAugmentedInChI_Wrong_Mult(self):
        """Test that augmented inchi with wrong layer is caught."""
        malform_aug_inchi = 'InChI=1S/CH3/h1H3'
        with self.assertRaises(Exception):
            mol = Molecule().fromAugmentedInChI(malform_aug_inchi)

    def testMalformedAugmentedInChI_Wrong_Indices(self):
        """Test that augmented inchi with wrong layer is caught."""
        malform_aug_inchi = 'InChI=1S/C6H6/c1-3-5-6-4-2/h1,6H,2,5H2/u4,1'
        with self.assertRaises(Exception):
            mol = Molecule().fromAugmentedInChI(malform_aug_inchi)

    def testUpdateLonePairs(self):
        adjlist = """
1 Si u0 p1 c0 {2,S} {3,S}
2 H  u0 p0 c0 {1,S}
3 H  u0 p0 c0 {1,S}
"""

        mol = Molecule().fromAdjacencyList(adjlist)
        mol.updateLonePairs()
        lp = 0
        for atom in mol.atoms:
            lp += atom.lonePairs
        self.assertEqual(lp, 1)
                    
    def testLargeMolUpdate(self):
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
        mol = Molecule().fromAdjacencyList(adjlist)

        mol.resetConnectivityValues()

        try:
            mol.updateConnectivityValues()
        except OverflowError:
            self.fail("updateConnectivityValues() raised OverflowError unexpectedly!")

    def testLargeMolCreation(self):
        """
        Test molecules between C1 to C201 in 10 carbon intervals to make
        sure that overflow errors are not generated.
        """
        for i in xrange(1,202,10):
            smi = 'C'*i
            try:
                m = Molecule(SMILES=smi)
            except OverflowError:
                self.fail('Creation of C{} failed!'.format(i))

    def testGetPolycyclicRings(self):
        """
        Test that polycyclic rings within a molecule are returned properly in the function
        `Graph().getPolycyclicRings()`
        """
        # norbornane
        m1 = Molecule(SMILES='C1CC2CCC1C2')
        polyrings1 = m1.getPolycyclicRings()
        self.assertEqual(len(polyrings1), 1)
        ring = polyrings1[0]
        self.assertEqual(len(ring),7)  # 7 carbons in cycle
        
        # dibenzyl
        m2 = Molecule(SMILES='C1=CC=C(C=C1)CCC1C=CC=CC=1')
        polyrings2 = m2.getPolycyclicRings()
        self.assertEqual(len(polyrings2), 0)
        
        # spiro[2.5]octane
        m3 = Molecule(SMILES='C1CCC2(CC1)CC2')
        polyrings3 = m3.getPolycyclicRings()
        self.assertEqual(len(polyrings3), 1)
        ring = polyrings3[0]
        self.assertEqual(len(ring),8)
        
        # 1-phenyl norbornane
        m4 = Molecule(SMILES='C1=CC=C(C=C1)C12CCC(CC1)C2')
        polyrings4 = m4.getPolycyclicRings()
        self.assertEqual(len(polyrings4), 1)
        ring = polyrings4[0]
        self.assertEqual(len(ring),7)
        
    def testGetMonocyclicRings(self):
        """
        Test that monocyclic rings within a molecule are returned properly in the function
        `Graph().getMonocyclicRings()`
        """
        m1 = Molecule(SMILES='C(CCCC1CCCCC1)CCCC1CCCC1')
        monorings = m1.getMonocyclicRings()
        self.assertEqual(len(monorings),2)
        
        m2 = Molecule(SMILES='C(CCC1C2CCC1CC2)CC1CCC1')
        monorings = m2.getMonocyclicRings()
        self.assertEqual(len(monorings),1)
        self.assertEqual(len(monorings[0]),4)
        
        m3 = Molecule(SMILES='CCCCC')
        monorings = m3.getMonocyclicRings()
        self.assertEqual(len(monorings),0)
        
    def testGetDisparateRings(self):
        """
        Test that monocyclic rings within a molecule are returned properly in the function
        `Graph().getDisparateRings()`
        """
        
        # norbornane
        m1 = Molecule(SMILES='C1CC2CCC1C2')
        monorings, polyrings = m1.getDisparateRings()
        self.assertEqual(len(monorings), 0)
        self.assertEqual(len(polyrings), 1)
        self.assertEqual(len(polyrings[0]), 7)  # 7 carbons in cycle

        # norbornane + cyclobutane on chain
        m2 = Molecule(SMILES='C(CCC1C2CCC1CC2)CC1CCC1')
        monorings, polyrings = m2.getDisparateRings()
        self.assertEqual(len(monorings), 1)
        self.assertEqual(len(polyrings), 1)
        self.assertEqual(len(monorings[0]), 4)
        self.assertEqual(len(polyrings[0]), 7)

        # spiro-octane + cyclobutane on chain
        m3 = Molecule(SMILES='C1CCC2(CC1)CC2CCCCC1CCC1')
        monorings, polyrings = m3.getDisparateRings()
        self.assertEqual(len(polyrings), 1)
        self.assertEqual(len(monorings), 1)
        self.assertEqual(len(monorings[0]), 4)
        self.assertEqual(len(polyrings[0]), 8)

        # butane
        m4 = Molecule(SMILES='CCCC')
        monorings, polyrings = m4.getDisparateRings()
        self.assertEqual(len(monorings), 0)
        self.assertEqual(len(polyrings), 0)

        # benzene + cyclopropane on chain + cyclopropane on chain
        m5 = Molecule(SMILES='C1=CC=C(CCCC2CC2)C(=C1)CCCCCC1CC1')
        monorings, polyrings = m5.getDisparateRings()
        self.assertEqual(len(monorings), 3)
        self.assertEqual(len(polyrings), 0)

        # octacene
        m6 = Molecule(SMILES='c1ccc2cc3cc4cc5cc6cc7cc8ccccc8cc7cc6cc5cc4cc3cc2c1')
        monorings, polyrings = m6.getDisparateRings()
        self.assertEqual(len(monorings), 0)
        self.assertEqual(len(polyrings), 1)
        self.assertEqual(len(polyrings[0]), 34)

        # JP-10
        m7 = Molecule(SMILES='C1CC2C3CCC(C3)C2C1')
        monorings, polyrings = m7.getDisparateRings()
        self.assertEqual(len(monorings), 0)
        self.assertEqual(len(polyrings), 1)
        self.assertEqual(len(polyrings[0]), 10)

    def testGetSmallestSetOfSmallestRings(self):
        """
        Test that SSSR within a molecule are returned properly in the function
        `Graph().getSmallestSetOfSmallestRings()`
        """

        m1 = Molecule(SMILES='C12CCC1C3CC2CC3')
        sssr1 = m1.getSmallestSetOfSmallestRings()
        sssr1_sizes = sorted([len(ring) for ring in sssr1])
        sssr1_sizes_expected = [4, 5, 5]
        self.assertEqual(sssr1_sizes, sssr1_sizes_expected)
        
        m2 = Molecule(SMILES='C1(CC2)C(CC3)CC3C2C1')
        sssr2 = m2.getSmallestSetOfSmallestRings()
        sssr2_sizes = sorted([len(ring) for ring in sssr2])
        sssr2_sizes_expected = [5, 5, 6]
        self.assertEqual(sssr2_sizes, sssr2_sizes_expected)
        
        
        m3 = Molecule(SMILES='C1(CC2)C2C(CCCC3)C3C1')
        sssr3 = m3.getSmallestSetOfSmallestRings()
        sssr3_sizes = sorted([len(ring) for ring in sssr3])
        sssr3_sizes_expected = [4, 5, 6]
        self.assertEqual(sssr3_sizes, sssr3_sizes_expected)
        
        m4 = Molecule(SMILES='C12=CC=CC=C1C3=C2C=CC=C3')
        sssr4 = m4.getSmallestSetOfSmallestRings()
        sssr4_sizes = sorted([len(ring) for ring in sssr4])
        sssr4_sizes_expected = [4, 6, 6]
        self.assertEqual(sssr4_sizes, sssr4_sizes_expected)
        
        m5 = Molecule(SMILES='C12=CC=CC=C1CC3=C(C=CC=C3)C2')
        sssr5 = m5.getSmallestSetOfSmallestRings()
        sssr5_sizes = sorted([len(ring) for ring in sssr5])
        sssr5_sizes_expected = [6, 6, 6]
        self.assertEqual(sssr5_sizes, sssr5_sizes_expected)
    
    def testGetDeterministicSmallestSetOfSmallestRingsCase1(self):
        """
        Test fused tricyclic can be decomposed into single rings more 
        deterministically
        """
        smiles = 'C1C2C3C=CCCC2C13'

        previous_num_shared_atoms_list = None
        # repeat 100 time to test non-deterministic behavior
        for _ in range(100):
            mol =  Molecule().fromSMILES(smiles)
            sssr_det = mol.getDeterministicSmallestSetOfSmallestRings()

            
            num_shared_atoms_list = []
            for i, ring_i in enumerate(sssr_det):
                for j in range(i+1, len(sssr_det)):
                    ring_j = sssr_det[j]
                    num_shared_atoms = len(set(ring_i).intersection(ring_j))

                    num_shared_atoms_list.append(num_shared_atoms)

            num_shared_atoms_list = sorted(num_shared_atoms_list)
            
            if previous_num_shared_atoms_list is None:
                previous_num_shared_atoms_list = num_shared_atoms_list
                continue
            self.assertEqual(num_shared_atoms_list, previous_num_shared_atoms_list)
            previous_num_shared_atoms_list = num_shared_atoms_list

    def testGetDeterministicSmallestSetOfSmallestRingsCase2(self):
        """
        Test if two possible smallest rings can join the smallest set
        the method can pick one of them deterministically using sum of 
        atomic numbers along the rings.
        In this test case and with currect method setup, ring (CCSCCCCC)
        will be picked rather than ring(CCCOCC).
        """

        smiles = 'C1=CC2C3CSC(CO3)C2C1'

        previous_atom_symbols_list = None
        # repeat 100 time to test non-deterministic behavior
        for _ in range(100):
            mol =  Molecule().fromSMILES(smiles)
            sssr_det = mol.getDeterministicSmallestSetOfSmallestRings()

            atom_symbols_list = []
            for ring in sssr_det:
                atom_symbols = sorted([a.element.symbol for a in ring])
                atom_symbols_list.append(atom_symbols)

            atom_symbols_list = sorted(atom_symbols_list)

            if previous_atom_symbols_list is None:
                previous_atom_symbols_list = atom_symbols_list
                continue
            self.assertEqual(atom_symbols_list, previous_atom_symbols_list)
            previous_atom_symbols_list = atom_symbols_list

    @work_in_progress
    def testGetDeterministicSmallestSetOfSmallestRingsCase3(self):
        """
        Test if two possible smallest rings can join the smallest set
        the method can pick one of them deterministically when their
        sum of atomic numbers along the rings are also equal to each other.
        
        To break the tie, one option we have is to consider adding contributions
        from other parts of the molecule, such as atomic number weighted connectivity
        value and differentiate bond orders when calculating connectivity values.
        """
        smiles = 'C=1CC2C3CSC(O[Si]3)C2C1'

        previous_atom_symbols_list = None
        # repeat 100 time to test non-deterministic behavior
        for _ in range(100):
            mol =  Molecule().fromSMILES(smiles)
            sssr_det = mol.getDeterministicSmallestSetOfSmallestRings()

            atom_symbols_list = []
            for ring in sssr_det:
                atom_symbols = sorted([a.element.symbol for a in ring])
                atom_symbols_list.append(atom_symbols)

            atom_symbols_list = sorted(atom_symbols_list)

            if previous_atom_symbols_list is None:
                previous_atom_symbols_list = atom_symbols_list
                continue
            self.assertEqual(atom_symbols_list, previous_atom_symbols_list)
            previous_atom_symbols_list = atom_symbols_list

    def testToGroup(self):
        """
        Test if we can convert a Molecule object into a Group object.
        """
        mol = Molecule().fromSMILES('CC(C)CCCC(C)C1CCC2C3CC=C4CC(O)CCC4(C)C3CCC12C')#cholesterol
        group = mol.toGroup()
        
        self.assertTrue(isinstance(group, Group))
        
        self.assertEquals(len(mol.atoms), len(group.atoms))

        molbondcount = sum([1 for atom in mol.atoms for bondedAtom, bond in atom.edges.iteritems()])
        groupbondcount = sum([1 for atom in group.atoms for bondedAtom, bond in atom.edges.iteritems()])
        self.assertEquals(molbondcount, groupbondcount)

        for i, molAt in enumerate(mol.atoms):
            groupAtom = group.atoms[i]
            atomTypes = [groupAtomType.equivalent(molAt.atomType) for groupAtomType in groupAtom.atomType]
            self.assertTrue(any(atomTypes))

    def testToAdjacencyListWithIsotopes(self):
        """
        Test the Molecule.toAdjacencyList() method works for atoms with unexpected isotopes.
        """

        mol = Molecule().fromSMILES('CC')
        mol.atoms[0].element = getElement('C', 13)

        adjlist = mol.toAdjacencyList().translate(None, '\n ')
        adjlistExp = """
        1 C u0 p0 c0 i13 {2,S} {3,S} {4,S} {5,S}
        2 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
        3 H u0 p0 c0 {1,S}
        4 H u0 p0 c0 {1,S}
        5 H u0 p0 c0 {1,S}
        6 H u0 p0 c0 {2,S}
        7 H u0 p0 c0 {2,S}
        8 H u0 p0 c0 {2,S}
        """.translate(None, '\n ')
        
        self.assertEquals(adjlist, adjlistExp)

        mol = Molecule().fromSMILES('CC')
        mol.atoms[2].element = getElement('H', 2)

        adjlist = mol.toAdjacencyList().translate(None, '\n ')
        adjlistExp = """
        1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
        2 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
        3 H u0 p0 c0 i2 {1,S}
        4 H u0 p0 c0 {1,S}
        5 H u0 p0 c0 {1,S}
        6 H u0 p0 c0 {2,S}
        7 H u0 p0 c0 {2,S}
        8 H u0 p0 c0 {2,S}
        """.translate(None, '\n ')
        
        self.assertEquals(adjlist, adjlistExp)


        mol = Molecule().fromSMILES('OC')
        mol.atoms[0].element = getElement('O', 18)

        adjlist = mol.toAdjacencyList().translate(None, '\n ')
        adjlistExp = """
        1 O u0 p2 c0 i18 {2,S} {3,S}
        2 C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
        3 H u0 p0 c0 {1,S}
        4 H u0 p0 c0 {2,S}
        5 H u0 p0 c0 {2,S}
        6 H u0 p0 c0 {2,S}
        """.translate(None, '\n ')
        
        self.assertEquals(adjlist, adjlistExp)

    def testFromAdjacencyListWithIsotopes(self):
        """
        Test the Molecule.fromAdjacencyList() method works for atoms with unexpected isotopes.
        """

        exp = Molecule().fromSMILES('CC')
        exp.atoms[0].element = getElement('C', 13)

        adjlistCalc = """
        1 C u0 p0 c0 i13 {2,S} {3,S} {4,S} {5,S}
        2 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
        3 H u0 p0 c0 {1,S}
        4 H u0 p0 c0 {1,S}
        5 H u0 p0 c0 {1,S}
        6 H u0 p0 c0 {2,S}
        7 H u0 p0 c0 {2,S}
        8 H u0 p0 c0 {2,S}
        """
        calc = Molecule().fromAdjacencyList(adjlistCalc)
        
        self.assertTrue(exp.isIsomorphic(calc))

        exp = Molecule().fromSMILES('CC')
        exp.atoms[2].element = getElement('H', 2)

        adjlistCalc = """
        1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
        2 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
        3 H u0 p0 c0 i2 {1,S}
        4 H u0 p0 c0 {1,S}
        5 H u0 p0 c0 {1,S}
        6 H u0 p0 c0 {2,S}
        7 H u0 p0 c0 {2,S}
        8 H u0 p0 c0 {2,S}
        """
        calc = Molecule().fromAdjacencyList(adjlistCalc)
        
        self.assertTrue(exp.isIsomorphic(calc))

        exp = Molecule().fromSMILES('OC')
        exp.atoms[0].element = getElement('O', 18)

        adjlistCalc = """
        1 O u0 p2 c0 i18 {2,S} {3,S}
        2 C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
        3 H u0 p0 c0 {1,S}
        4 H u0 p0 c0 {2,S}
        5 H u0 p0 c0 {2,S}
        6 H u0 p0 c0 {2,S}
        """
        calc = Molecule().fromAdjacencyList(adjlistCalc)
        
        self.assertTrue(exp.isIsomorphic(calc))

    def testAromaticityPerceptionBenzene(self):
        """Test aromaticity perception via getAromaticRings for benzene."""
        mol = Molecule(SMILES='c1ccccc1')
        aromaticAtoms, aromaticBonds = mol.getAromaticRings()
        self.assertEqual(len(aromaticAtoms), 1)
        self.assertEqual(len(aromaticBonds), 1)
        for bond in aromaticBonds[0]:
            self.assertTrue(bond.atom1 in aromaticAtoms[0] and bond.atom2 in aromaticAtoms[0])

    def testAromaticityPerceptionTetralin(self):
        """Test aromaticity perception via getAromaticRings for tetralin."""
        mol = Molecule(SMILES='c1ccc2c(c1)CCCC2')
        aromaticAtoms, aromaticBonds = mol.getAromaticRings()
        self.assertEqual(len(aromaticAtoms), 1)
        self.assertEqual(len(aromaticBonds), 1)
        for bond in aromaticBonds[0]:
            self.assertTrue(bond.atom1 in aromaticAtoms[0] and bond.atom2 in aromaticAtoms[0])

    def testAromaticityPerceptionBiphenyl(self):
        """Test aromaticity perception via getAromaticRings for biphenyl."""
        mol = Molecule(SMILES='c1ccc(cc1)c2ccccc2')
        aromaticAtoms, aromaticBonds = mol.getAromaticRings()
        self.assertEqual(len(aromaticAtoms), 2)
        self.assertEqual(len(aromaticBonds), 2)
        for index in range(len(aromaticAtoms)):
            for bond in aromaticBonds[index]:
                self.assertTrue(bond.atom1 in aromaticAtoms[index] and bond.atom2 in aromaticAtoms[index])

    def testAromaticityPerceptionAzulene(self):
        """Test aromaticity perception via getAromaticRings for azulene."""
        mol = Molecule(SMILES='c1cccc2cccc2c1')
        aromaticAtoms, aromaticBonds = mol.getAromaticRings()
        self.assertEqual(len(aromaticAtoms), 0)
        self.assertEqual(len(aromaticBonds), 0)

    def testAromaticityPerceptionFuran(self):
        """Test aromaticity perception via getAromaticRings for furan."""
        mol = Molecule(SMILES='c1ccoc1')
        aromaticAtoms, aromaticBonds = mol.getAromaticRings()
        self.assertEqual(len(aromaticAtoms), 0)
        self.assertEqual(len(aromaticBonds), 0)

    def testArylRadicalTrue(self):
        """Test aryl radical perception for phenyl radical."""
        mol = Molecule(SMILES='[c]1ccccc1')
        self.assertTrue(mol.isArylRadical())

    def testArylRadicalFalse(self):
        """Test aryl radical perception for benzyl radical."""
        mol = Molecule(SMILES='[CH2]c1ccccc1')
        self.assertFalse(mol.isArylRadical())

    def testArylRadicalBirad(self):
        """Test aryl radical perception for biradical species.

        This is a case that is not properly handled right now, since a single boolean cannot
        characterize multiple radicals. In such cases, the method will return false if
        any of the radicals is not an aryl radical."""
        mol = Molecule(SMILES='[CH2]c1c[c]ccc1')
        self.assertFalse(mol.isArylRadical())

    def testIdenticalTrue(self):
        """Test that the isIdentical returns True with butane"""
        mol = Molecule(SMILES='CCCC')
        mol.assignAtomIDs()
        molCopy = mol.copy(deep=True)
        self.assertTrue(mol.isIsomorphic(molCopy))
        self.assertTrue(mol.isIdentical(molCopy))

    def testIdenticalTrue2(self):
        """Test that isIdentical with strict=False returns True with allyl"""
        mol = Molecule(SMILES='C=C[CH2]')
        mol.assignAtomIDs()
        res = mol.generate_resonance_structures(keep_isomorphic=True)
        self.assertEqual(len(res), 2)

        mol2 = res[1]
        self.assertTrue(mol.isIsomorphic(mol2))
        self.assertFalse(mol.isIdentical(mol2))
        self.assertTrue(mol.isIdentical(mol2, strict=False))

    def testIdenticalFalse(self):
        """Test that the isIdentical returns False with butane"""
        mol = Molecule(SMILES='CCCC')
        mol.assignAtomIDs()
        molCopy = mol.copy(deep=True)
        # Remove a hydrogen from mol
        a = mol.atoms[-1]

        mol.removeAtom(a)
        # Remove a different hydrogen from molCopy
        b = molCopy.atoms[-2]

        molCopy.removeAtom(b)

        self.assertTrue(mol.isIsomorphic(molCopy))
        self.assertFalse(mol.isIdentical(molCopy))

    def testIdenticalFalse2(self):
        """Test that the isIdentical method returns False with ethene"""
        # Manually test addition of H radical to ethene
        reactant1 = Molecule(SMILES='C=C')
        carbons = [atom for atom in reactant1.atoms if atom.symbol == 'C']
        carbons[0].label = '*1'
        carbons[1].label = '*2'
        reactant2 = Molecule(SMILES='[H]')
        reactant2.atoms[0].label = '*3'
        # Merge reactants
        mol = reactant1.merge(reactant2)
        mol.assignAtomIDs()
        molCopy = mol.copy(deep=True)
        # Manually perform R_Addition_MultipleBond of *3 to *1
        labeledAtoms = mol.getLabeledAtoms()
        mol.getBond(labeledAtoms['*1'], labeledAtoms['*2']).decrementOrder()
        mol.addBond(Bond(labeledAtoms['*1'], labeledAtoms['*3'], order='S'))
        labeledAtoms['*2'].incrementRadical()
        labeledAtoms['*3'].decrementRadical()
        # Manually perform R_Addition_MultipleBond of *3 to *2
        labeledAtoms = molCopy.getLabeledAtoms()
        molCopy.getBond(labeledAtoms['*1'], labeledAtoms['*2']).decrementOrder()
        molCopy.addBond(Bond(labeledAtoms['*2'], labeledAtoms['*3'], order='S'))
        labeledAtoms['*1'].incrementRadical()
        labeledAtoms['*3'].decrementRadical()

        self.assertTrue(mol.isIsomorphic(molCopy))
        self.assertFalse(mol.isIdentical(molCopy))

    def testatomidvalid(self):
        """see if the atomIDVvalid method properly returns True"""
        mol = Molecule(SMILES='CCCC')
        for index, atom in enumerate(mol.atoms):
            atom.id =index
        self.assertTrue(mol.atomIDValid())

    def testatomidvalid2(self):
        """see if the atomIDVvalid method properly returns False"""
        mol = Molecule(SMILES='CCCC')
        for index, atom in enumerate(mol.atoms):
            atom.index =index
        mol.atoms[3].index = 4
        self.assertFalse(mol.atomIDValid())

    def testatomidvalid2(self):
        """see if the atomIDVvalid method properly returns False"""
        mol = Molecule(SMILES='CCCC')
        self.assertFalse(mol.atomIDValid())

    def testassignatomid(self):
        """see if the assignAtomID method properly labels molecule"""
        mol = Molecule(SMILES='CCCC')
        mol.assignAtomIDs()
        self.assertTrue(mol.atomIDValid())

    def testFingerprintProperty(self):
        """Test that the Molecule.fingerprint property works"""
        # Test getting fingerprint
        self.assertEqual(self.molecule[0].fingerprint, 'CH2NO2')

        # Test setting fingerprint
        self.molecule[0].fingerprint = 'nitronate'
        self.assertEqual(self.molecule[0].fingerprint, 'nitronate')

    def testSaturateUnfilledValence(self):
        """
        Test the saturateUnfilledValence for an aromatic and nonaromatic case
        """
        #test butane
        expected = Molecule(SMILES='CCCC')
        test = expected.copy(deep = True)
        test.deleteHydrogens()

        hydrogens = 0
        for atom in test.atoms:
            if atom.isHydrogen(): hydrogens +=1
        self.assertEquals(hydrogens, 0)

        test.saturate_unfilled_valence()

        hydrogens = 0
        for atom in test.atoms:
            if atom.isHydrogen(): hydrogens +=1
        self.assertEquals(hydrogens, 10)

        test.update()
        self.assertTrue(expected.isIsomorphic(test))

        #test benzene
        expected = Molecule(SMILES='c1ccccc1')
        test = expected.copy(deep = True)
        test.deleteHydrogens()
        hydrogens = 0
        for atom in test.atoms:
            if atom.isHydrogen(): hydrogens +=1
        self.assertEquals(hydrogens, 0)

        test.saturate_unfilled_valence()

        hydrogens = 0
        for atom in test.atoms:
            if atom.isHydrogen(): hydrogens +=1
        self.assertEquals(hydrogens, 6)

        test.update()
        self.assertTrue(expected.isIsomorphic(test))

    def test_get_element_count(self):
        """Test that we can count elements properly."""
        mol1 = Molecule(SMILES='c1ccccc1')
        expected1 = {'C': 6, 'H': 6}
        result1 = mol1.get_element_count()
        self.assertEqual(expected1, result1)

        mol2 = Molecule(SMILES='CS(C)(=O)=O')
        expected2 = {'C': 2, 'H': 6, 'O': 2, 'S': 1}
        result2 = mol2.get_element_count()
        self.assertEqual(expected2, result2)

        mol3 = Molecule(SMILES='CCN')
        expected3 = {'C': 2, 'H': 7, 'N': 1}
        result3 = mol3.get_element_count()
        self.assertEqual(expected3, result3)

    def testRingPerception(self):
        """Test that identifying ring membership of atoms works properly."""
        mol = Molecule(SMILES='c12ccccc1cccc2')
        mol.identifyRingMembership()
        for atom in mol.atoms:
            if atom.element == 'C':
                self.assertTrue(atom.props['inRing'])
            elif atom.element == 'H':
                self.assertFalse(atom.props['inRing'])

    def test_enumerate_bonds(self):
        """Test that generating a count of bond labels works properly."""
        adj_list = '''
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
        '''
        mol = Molecule().fromAdjacencyList(adj_list)
        bonds = mol.enumerate_bonds()
        self.assertEqual(bonds['C#C'], 1)
        self.assertEqual(bonds['C-C'], 4)
        self.assertEqual(bonds['C-H'], 10)
        self.assertEqual(bonds['C-O'], 2)
        self.assertEqual(bonds['C:C'], 6)
        self.assertEqual(bonds['H-O'], 2)
        self.assertEqual(bonds['H~O'], 2)

################################################################################

if __name__ == '__main__':
    unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )
