#!/usr/bin/python
# -*- coding: utf-8 -*-

import unittest

from external.wip import work_in_progress
from .molecule import Atom, Bond, Molecule
from .group import Group, ActionError
from .element import getElement, elementList
from .resonance import generateAromaticResonanceIsomers


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
                self.assertTrue(atom.isNonHydrogen())

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

    def testIsOxygen(self):
        """
        Test the Atom.isOxygen() method.
        """
        for element in elementList:
            atom = Atom(element=element, radicalElectrons=1, charge=0, label='*1', lonePairs=2)
            if element.symbol == 'O':
                self.assertTrue(atom.isOxygen())
            else:
                self.assertFalse(atom.isOxygen())

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
        action = ['BREAK_BOND', '*1', 'S', '*2']
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
        action = ['FORM_BOND', '*1', 'S', '*2']
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

################################################################################

class TestBond(unittest.TestCase):
    """
    Contains unit tests of the Bond class.
    """

    def setUp(self):
        """
        A method called before each unit test in this class.
        """
        self.bond = Bond(atom1=None, atom2=None, order='D')
        self.orderList = ['S','D','T','B']
    
    def testIsSingle(self):
        """
        Test the Bond.isSingle() method.
        """
        for order in self.orderList:
            bond = Bond(None, None, order=order)
            if order == 'S':
                self.assertTrue(bond.isSingle())
            else:
                self.assertFalse(bond.isSingle())
    
    def testIsDouble(self):
        """
        Test the Bond.isDouble() method.
        """
        for order in self.orderList:
            bond = Bond(None, None, order=order)
            if order == 'D':
                self.assertTrue(bond.isDouble())
            else:
                self.assertFalse(bond.isDouble())

    def testIsTriple(self):
        """
        Test the Bond.isTriple() method.
        """
        for order in self.orderList:
            bond = Bond(None, None, order=order)
            if order == 'T':
                self.assertTrue(bond.isTriple())
            else:
                self.assertFalse(bond.isTriple())

    def testIsBenzene(self):
        """
        Test the Bond.isBenzene() method.
        """
        for order in self.orderList:
            bond = Bond(None, None, order=order)
            if order == 'B':
                self.assertTrue(bond.isBenzene())
            else:
                self.assertFalse(bond.isBenzene())

    def testIncrementOrder(self):
        """
        Test the Bond.incrementOrder() method.
        """
        for order in self.orderList:
            bond = Bond(None, None, order=order)
            try:
                bond.incrementOrder()
                if order == 'S': 
                    self.assertTrue(bond.isDouble())
                elif order == 'D': 
                    self.assertTrue(bond.isTriple())
            except ActionError:
                self.assertTrue(order in ['T','B'])
        
    def testDecrementOrder(self):
        """
        Test the Bond.decrementOrder() method.
        """
        for order in self.orderList:
            bond = Bond(None, None, order=order)
            try:
                bond.decrementOrder()
                if order == 'D': 
                    self.assertTrue(bond.isSingle())
                elif order == 'T': 
                    self.assertTrue(bond.isDouble())
            except ActionError:
                self.assertTrue(order in ['S','B'])
                
    def testApplyActionBreakBond(self):
        """
        Test the Bond.applyAction() method for a BREAK_BOND action.
        """
        action = ['BREAK_BOND', '*1', 'S', '*2']
        for order0 in self.orderList:
            bond0 = Bond(None, None, order=order0)
            bond = bond0.copy()
            try:
                bond.applyAction(action)
                self.fail('Bond.applyAction() unexpectedly processed a BREAK_BOND action.')
            except ActionError:
                pass
    
    def testApplyActionFormBond(self):
        """
        Test the Bond.applyAction() method for a FORM_BOND action.
        """
        action = ['FORM_BOND', '*1', 'S', '*2']
        for order0 in self.orderList:
            bond0 = Bond(None, None, order=order0)
            bond = bond0.copy()
            try:
                bond.applyAction(action)
                self.fail('Bond.applyAction() unexpectedly processed a FORM_BOND action.')
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
                self.assertTrue('T' == order0 or 'B' == order0)
                
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
                self.assertTrue('S' == order0 or 'B' == order0)
            
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
                self.fail('Bond.applyAction() unexpectedly processed a GAIN_RADICAL action.')
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
                self.fail('Bond.applyAction() unexpectedly processed a LOSE_RADICAL action.')
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
        Test that the species [CH] has three radical electrons and a spin multiplicity of 4.
        """
        molecule = Molecule().fromSMILES('[CH]')
        self.assertEqual(molecule.atoms[0].radicalElectrons, 3)
        self.assertEqual(molecule.multiplicity, 4)
        self.assertEqual(molecule.getRadicalCount(), 3)

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
        
    def testSMILES(self):
        """
        Test that we can generate a few SMILES strings as expected
        """
        import rmgpy.molecule
        test_strings = ['[C-]#[O+]', '[C]', '[CH]', 'OO', '[H][H]', '[H]',
                       '[He]', '[O]', 'O', '[CH3]', 'C', '[OH]', 'CCC',
                       'CC', 'N#N', '[O]O', 'C[CH2]', '[Ar]', 'CCCC',
                       'O=C=O', 'N#[C]',
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
                    "CC1C=CC=CC=1O":"""
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
                        16 H u0 p0 c0 {8,S}""",
                    "CC1=CC=CC=C1O":"""
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
                        16 H u0 p0 c0 {8,S}""",
                    "CC1C=CC=CC=1":"""
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
                        15 H u0 p0 c0 {7,S}"""
                    }
        for smiles, adjlist in test_cases.iteritems():
            m = Molecule().fromAdjacencyList(adjlist)
            s = m.toSMILES()
            self.assertEqual(s, smiles, "Generated SMILES string {0} instead of {1}".format(s, smiles))
        

    def testKekuleRoundTripSMILES(self):
        """
        Test that we can round-trip SMILES strings of Kekulized aromatics
        """
        import rmgpy.molecule
        test_strings = [
                       'CC1=CC=CC=C1O', 'CC1C=CC=CC=1O',
                       # 'Cc1ccccc1O', # this will fail because it is Kekulized during fromSMILES()
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
        self.assertEqual(key, 'UMRZSTCPUPJPOJ-UHFFFAOYSA')
        
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
        
        self.assertEqual(mol.toAugmentedInChIKey(), 'VGGSQFUCUMXWEO-UHFFFAOYSA-u1,2')

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
        isomers = m.generateResonanceIsomers()
        self.assertTrue(any(isomer.isAromatic() for isomer in isomers))

    def testAromaticNaphthalene(self):
        """
        Test the Molecule.isAromatic() method for Naphthalene.
        """
        m = Molecule().fromSMILES('C12C(C=CC=C1)=CC=CC=2')
        isomers = m.generateResonanceIsomers()
        self.assertTrue(any(isomer.isAromatic() for isomer in isomers))
                        
    def testAromaticCyclohexane(self):
        """
        Test the Molecule.isAromatic() method for Cyclohexane.
        """
        m = Molecule().fromSMILES('C1CCCCC1')
        isomers = m.generateResonanceIsomers()
        self.assertFalse(any(isomer.isAromatic() for isomer in isomers))
         
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
        saturated_molecule.saturate()
        self.assertTrue(saturated_molecule.isIsomorphic(indene))
        
    def testFusedAromatic1(self):
        """Test we can make aromatic perylene from both adjlist and SMILES"""
        perylene = Molecule().fromAdjacencyList("""
1  C u0 p0 c0 {3,B} {6,B} {7,B}
2  C u0 p0 c0 {4,B} {5,B} {8,B}
3  C u0 p0 c0 {1,B} {4,B} {11,B}
4  C u0 p0 c0 {2,B} {3,B} {12,B}
5  C u0 p0 c0 {2,B} {6,B} {15,B}
6  C u0 p0 c0 {1,B} {5,B} {16,B}
7  C u0 p0 c0 {1,B} {9,B} {10,B}
8  C u0 p0 c0 {2,B} {13,B} {14,B}
9  C u0 p0 c0 {7,B} {17,B} {22,S}
10 C u0 p0 c0 {7,B} {18,B} {23,S}
11 C u0 p0 c0 {3,B} {18,B} {25,S}
12 C u0 p0 c0 {4,B} {19,B} {26,S}
13 C u0 p0 c0 {8,B} {19,B} {28,S}
14 C u0 p0 c0 {8,B} {20,B} {29,S}
15 C u0 p0 c0 {5,B} {20,B} {31,S}
16 C u0 p0 c0 {6,B} {17,B} {32,S}
17 C u0 p0 c0 {9,B} {16,B} {21,S}
18 C u0 p0 c0 {10,B} {11,B} {24,S}
19 C u0 p0 c0 {12,B} {13,B} {27,S}
20 C u0 p0 c0 {14,B} {15,B} {30,S}
21 H u0 p0 c0 {17,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {18,S}
25 H u0 p0 c0 {11,S}
26 H u0 p0 c0 {12,S}
27 H u0 p0 c0 {19,S}
28 H u0 p0 c0 {13,S}
29 H u0 p0 c0 {14,S}
30 H u0 p0 c0 {20,S}
31 H u0 p0 c0 {15,S}
32 H u0 p0 c0 {16,S}
""")
        perylene2 = Molecule().fromSMILES('c1cc2cccc3c4cccc5cccc(c(c1)c23)c54')
        for isomer in generateAromaticResonanceIsomers(perylene2):
            if perylene.isIsomorphic(isomer):
                break
        else:  # didn't break
            self.fail("{} isn't isomorphic with any aromatic forms of {}".format(
                            perylene.toSMILES(),
                            perylene2.toSMILES()
                        ))

    def testFusedAromatic2(self):
        """Test we can make aromatic naphthalene from both adjlist and SMILES"""
        naphthalene = Molecule().fromAdjacencyList("""
1  C u0 p0 c0 {2,B} {3,B} {4,B}
2  C u0 p0 c0 {1,B} {5,B} {6,B}
3  C u0 p0 c0 {1,B} {8,B} {13,S}
4  C u0 p0 c0 {1,B} {9,B} {14,S}
5  C u0 p0 c0 {2,B} {10,B} {17,S}
6  C u0 p0 c0 {2,B} {7,B} {18,S}
7  C u0 p0 c0 {6,B} {8,B} {11,S}
8  C u0 p0 c0 {3,B} {7,B} {12,S}
9  C u0 p0 c0 {4,B} {10,B} {15,S}
10 C u0 p0 c0 {5,B} {9,B} {16,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
""")
        naphthalene2 = Molecule().fromSMILES('C1=CC=C2C=CC=CC2=C1')
        for isomer in generateAromaticResonanceIsomers(naphthalene2):
            if naphthalene.isIsomorphic(isomer):
                break
        else:  # didn't break
            self.fail("{} isn't isomorphic with any aromatic forms of {}".format(
                            naphthalene.toSMILES(),
                            naphthalene2.toSMILES()
                        ))

    def testFusedAromatic3(self):
        """Test we can make aromatic pyrene_rad from both adjlist and SMILES"""
        pyrene_rad = Molecule().fromAdjacencyList("""
multiplicity 2
1  C u0 p0 c0 {2,B} {3,B} {5,B}
2  C u0 p0 c0 {1,B} {4,B} {6,S}
3  C u0 p0 c0 {1,B} {8,B} {9,B}
4  C u0 p0 c0 {2,B} {10,B} {11,S}
5  C u0 p0 c0 {1,B} {7,B} {15,S}
6  C u0 p0 c0 {2,S} {12,S} {16,D}
7  C u0 p0 c0 {5,B} {13,B} {17,S}
8  C u0 p0 c0 {3,B} {13,B} {19,S}
9  C u0 p0 c0 {3,B} {10,B} {20,S}
10 C u0 p0 c0 {4,B} {9,B} {21,S}
11 C u1 p0 c0 {4,S} {14,S} {22,S}
12 C u0 p0 c0 {6,S} {14,D} {24,S}
13 C u0 p0 c0 {7,B} {8,B} {18,S}
14 C u0 p0 c0 {11,S} {12,D} {23,S}
15 C u0 p0 c0 {5,S} {16,D} {25,S}
16 C u0 p0 c0 {6,D} {15,D}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {13,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {11,S}
23 H u0 p0 c0 {14,S}
24 H u0 p0 c0 {12,S}
25 H u0 p0 c0 {15,S}
""")
        pyrene_rad2 = Molecule().fromSMILES('[C]1C=C2C=CC=C3C=CC4=CC=CC=1C4=C23')
        for isomer in pyrene_rad2.generateResonanceIsomers():
            if pyrene_rad.isIsomorphic(isomer):
                break
        else:  # didn't break
            self.fail("{} isn't isomorphic with any aromatic forms of {}".format(
                            pyrene_rad.toSMILES(),
                            pyrene_rad2.toSMILES()
                        ))

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

    def testRDKitMolAtomMapping(self):
        """
        Test that the atom mapping returned by toRDKitMol contains the correct
        atom indices of the atoms of the molecule when hydrogens are removed.
        """
        from .generator import toRDKitMol

        adjlist = '''
1 H u0 p0 c0 {2,S}
2 C u0 p0 c0 {1,S} {3,S} {4,S} {5,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {2,S}
5 O u0 p2 c0 {2,S} {6,S}
6 H u0 p0 c0 {5,S}
        '''

        mol = Molecule().fromAdjacencyList(adjlist)
        rdkitmol, rdAtomIndices = toRDKitMol(mol, removeHs=True, returnMapping=True)

        heavy_atoms = [at for at in mol.atoms if at.number != 1]
        for at1 in heavy_atoms:
            for at2 in heavy_atoms:
                if mol.hasBond(at1, at2):
                    try:
                        rdkitmol.GetBondBetweenAtoms(rdAtomIndices[at1],rdAtomIndices[at2])
                    except RuntimeError:
                        self.fail("RDKit failed in finding the bond in the original atom!")
    
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
        self.assertEqual(len(polyrings[0]),7)  # 7 carbons in cycle
        
        m2 = Molecule(SMILES='C(CCC1C2CCC1CC2)CC1CCC1')
        monorings, polyrings = m2.getDisparateRings()
        self.assertEqual(len(monorings),1)
        self.assertEqual(len(polyrings),1)
        self.assertEqual(len(monorings[0]),4)
        self.assertEqual(len(polyrings[0]),7)
        
        
        m3 = Molecule(SMILES='C1CCC2(CC1)CC2CCCCC1CCC1')
        monorings, polyrings = m3.getDisparateRings()
        self.assertEqual(len(polyrings), 1)
        self.assertEqual(len(monorings),1)
        self.assertEqual(len(monorings[0]),4)
        self.assertEqual(len(polyrings[0]),8)
        
        m4 = Molecule(SMILES='CCCC')
        monorings, polyrings = m4.getDisparateRings()
        self.assertEqual(len(monorings),0)
        self.assertEqual(len(polyrings),0)
        
        m5 = Molecule(SMILES='C1=CC=C(CCCC2CC2)C(=C1)CCCCCC1CC1')
        monorings, polyrings = m5.getDisparateRings()
        self.assertEqual(len(monorings),3)
        self.assertEqual(len(polyrings),0)

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

################################################################################

if __name__ == '__main__':
    unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )
