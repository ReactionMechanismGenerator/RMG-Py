#!/usr/bin/python
# -*- coding: utf-8 -*-

import unittest

from rmgpy.molecule import *
from rmgpy.group import Group
from rmgpy.element import getElement, elementList

################################################################################

class TestAtom(unittest.TestCase):
    """
    Contains unit tests of the Atom class.
    """

    def setUp(self):
        """
        A method called before each unit test in this class.
        """
        self.atom = Atom(element=getElement('C'), radicalElectrons=1, spinMultiplicity=2, implicitHydrogens=3, charge=0, label='*1')
    
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
            atom = Atom(element=element, radicalElectrons=1, spinMultiplicity=2, implicitHydrogens=0, charge=0, label='*1')
            if element.symbol == 'H':
                self.assertTrue(atom.isHydrogen())
            else:
                self.assertFalse(atom.isHydrogen())
    
    def testIsNonHydrogen(self):
        """
        Test the Atom.isNonHydrogen() method.
        """
        for element in elementList:
            atom = Atom(element=element, radicalElectrons=1, spinMultiplicity=2, implicitHydrogens=0, charge=0, label='*1')
            if element.symbol == 'H':
                self.assertFalse(atom.isNonHydrogen())
            else:
                self.assertTrue(atom.isNonHydrogen())

    def testIsCarbon(self):
        """
        Test the Atom.isCarbon() method.
        """
        for element in elementList:
            atom = Atom(element=element, radicalElectrons=1, spinMultiplicity=2, implicitHydrogens=0, charge=0, label='*1')
            if element.symbol == 'C':
                self.assertTrue(atom.isCarbon())
            else:
                self.assertFalse(atom.isCarbon())

    def testIsOxygen(self):
        """
        Test the Atom.isOxygen() method.
        """
        for element in elementList:
            atom = Atom(element=element, radicalElectrons=1, spinMultiplicity=2, implicitHydrogens=0, charge=0, label='*1')
            if element.symbol == 'O':
                self.assertTrue(atom.isOxygen())
            else:
                self.assertFalse(atom.isOxygen())

    def testIncrementRadical(self):
        """
        Test the Atom.incrementRadical() method.
        """
        radicalElectrons = self.atom.radicalElectrons
        spinMultiplicity = self.atom.spinMultiplicity
        self.atom.incrementRadical()
        self.assertEqual(self.atom.radicalElectrons, radicalElectrons + 1)
        self.assertEqual(self.atom.spinMultiplicity, spinMultiplicity + 1)
    
    def testDecrementRadical(self):
        """
        Test the Atom.decrementRadical() method.
        """
        radicalElectrons = self.atom.radicalElectrons
        spinMultiplicity = self.atom.spinMultiplicity
        self.atom.decrementRadical()
        self.assertEqual(self.atom.radicalElectrons, radicalElectrons - 1)
        self.assertEqual(self.atom.spinMultiplicity, spinMultiplicity - 1)
           
    def testApplyActionBreakBond(self):
        """
        Test the Atom.applyAction() method for a BREAK_BOND action.
        """
        action = ['BREAK_BOND', '*1', 'S', '*2']
        for element in elementList:
            atom0 = Atom(element=element, radicalElectrons=1, spinMultiplicity=2, implicitHydrogens=3, charge=0, label='*1')
            atom = atom0.copy()
            atom.applyAction(action)
            self.assertEqual(atom0.element, atom.element)
            self.assertEqual(atom0.radicalElectrons, atom.radicalElectrons)
            self.assertEqual(atom0.spinMultiplicity, atom.spinMultiplicity)
            self.assertEqual(atom0.implicitHydrogens, atom.implicitHydrogens)
            self.assertEqual(atom0.charge, atom.charge)
            self.assertEqual(atom0.label, atom.label)
    
    def testApplyActionFormBond(self):
        """
        Test the Atom.applyAction() method for a FORM_BOND action.
        """
        action = ['FORM_BOND', '*1', 'S', '*2']
        for element in elementList:
            atom0 = Atom(element=element, radicalElectrons=1, spinMultiplicity=2, implicitHydrogens=3, charge=0, label='*1')
            atom = atom0.copy()
            atom.applyAction(action)
            self.assertEqual(atom0.element, atom.element)
            self.assertEqual(atom0.radicalElectrons, atom.radicalElectrons)
            self.assertEqual(atom0.spinMultiplicity, atom.spinMultiplicity)
            self.assertEqual(atom0.implicitHydrogens, atom.implicitHydrogens)
            self.assertEqual(atom0.charge, atom.charge)
            self.assertEqual(atom0.label, atom.label)
    
    def testApplyActionIncrementBond(self):
        """
        Test the Atom.applyAction() method for a CHANGE_BOND action.
        """
        action = ['CHANGE_BOND', '*1', 1, '*2']
        for element in elementList:
            atom0 = Atom(element=element, radicalElectrons=1, spinMultiplicity=2, implicitHydrogens=3, charge=0, label='*1')
            atom = atom0.copy()
            atom.applyAction(action)
            self.assertEqual(atom0.element, atom.element)
            self.assertEqual(atom0.radicalElectrons, atom.radicalElectrons)
            self.assertEqual(atom0.spinMultiplicity, atom.spinMultiplicity)
            self.assertEqual(atom0.implicitHydrogens, atom.implicitHydrogens)
            self.assertEqual(atom0.charge, atom.charge)
            self.assertEqual(atom0.label, atom.label)
    
    def testApplyActionDecrementBond(self):
        """
        Test the Atom.applyAction() method for a CHANGE_BOND action.
        """
        action = ['CHANGE_BOND', '*1', -1, '*2']
        for element in elementList:
            atom0 = Atom(element=element, radicalElectrons=1, spinMultiplicity=2, implicitHydrogens=3, charge=0, label='*1')
            atom = atom0.copy()
            atom.applyAction(action)
            self.assertEqual(atom0.element, atom.element)
            self.assertEqual(atom0.radicalElectrons, atom.radicalElectrons)
            self.assertEqual(atom0.spinMultiplicity, atom.spinMultiplicity)
            self.assertEqual(atom0.implicitHydrogens, atom.implicitHydrogens)
            self.assertEqual(atom0.charge, atom.charge)
            self.assertEqual(atom0.label, atom.label)
    
    def testApplyActionGainRadical(self):
        """
        Test the Atom.applyAction() method for a GAIN_RADICAL action.
        """
        action = ['GAIN_RADICAL', '*1', 1]
        for element in elementList:
            atom0 = Atom(element=element, radicalElectrons=1, spinMultiplicity=2, implicitHydrogens=3, charge=0, label='*1')
            atom = atom0.copy()
            atom.applyAction(action)
            self.assertEqual(atom0.element, atom.element)
            self.assertEqual(atom0.radicalElectrons, atom.radicalElectrons - 1)
            self.assertEqual(atom0.spinMultiplicity, atom.spinMultiplicity - 1)
            self.assertEqual(atom0.implicitHydrogens, atom.implicitHydrogens)
            self.assertEqual(atom0.charge, atom.charge)
            self.assertEqual(atom0.label, atom.label)
    
    def testApplyActionLoseRadical(self):
        """
        Test the Atom.applyAction() method for a LOSE_RADICAL action.
        """
        action = ['LOSE_RADICAL', '*1', 1]
        for element in elementList:
            atom0 = Atom(element=element, radicalElectrons=1, spinMultiplicity=2, implicitHydrogens=3, charge=0, label='*1')
            atom = atom0.copy()
            atom.applyAction(action)
            self.assertEqual(atom0.element, atom.element)
            self.assertEqual(atom0.radicalElectrons, atom.radicalElectrons + 1)
            self.assertEqual(atom0.spinMultiplicity, atom.spinMultiplicity + 1)
            self.assertEqual(atom0.implicitHydrogens, atom.implicitHydrogens)
            self.assertEqual(atom0.charge, atom.charge)
            self.assertEqual(atom0.label, atom.label)
    
    def testEquivalent(self):
        """
        Test the Atom.equivalent() method.
        """
        for index1, element1 in enumerate(elementList[0:10]):
            for index2, element2 in enumerate(elementList[0:10]):
                atom1 = Atom(element=element1, radicalElectrons=1, spinMultiplicity=2, charge=0, label='*1')
                atom2 = Atom(element=element2, radicalElectrons=1, spinMultiplicity=2, charge=0, label='*1')
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
                atom1 = Atom(element=element1, radicalElectrons=1, spinMultiplicity=2, charge=0, label='*1')
                atom2 = Atom(element=element2, radicalElectrons=1, spinMultiplicity=2, charge=0, label='*1')
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
        self.assertEqual(self.atom.spinMultiplicity, atom.spinMultiplicity)
        self.assertEqual(self.atom.implicitHydrogens, atom.implicitHydrogens)
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
        self.assertEqual(self.atom.spinMultiplicity, atom.spinMultiplicity)
        self.assertEqual(self.atom.implicitHydrogens, atom.implicitHydrogens)
        self.assertEqual(self.atom.charge, atom.charge)
        self.assertEqual(self.atom.label, atom.label)
        
    def testOutput(self):
        """
        Test that we can reconstruct a Atom object from its repr()
        output with no loss of information.
        """
        exec('atom = {0!r}'.format(self.atom))
        self.assertEqual(self.atom.element.symbol, atom.element.symbol)
        self.assertEqual(self.atom.atomType, atom.atomType)
        self.assertEqual(self.atom.radicalElectrons, atom.radicalElectrons)
        self.assertEqual(self.atom.spinMultiplicity, atom.spinMultiplicity)
        self.assertEqual(self.atom.implicitHydrogens, atom.implicitHydrogens)
        self.assertEqual(self.atom.charge, atom.charge)
        self.assertEqual(self.atom.label, atom.label)

################################################################################

class TestBond(unittest.TestCase):
    """
    Contains unit tests of the Bond class.
    """

    def setUp(self):
        """
        A method called before each unit test in this class.
        """
        self.bond = Bond(order='D')
        self.orderList = ['S','D','T','B']
    
    def testIsSingle(self):
        """
        Test the Bond.isSingle() method.
        """
        for order in self.orderList:
            bond = Bond(order=order)
            if order == 'S':
                self.assertTrue(bond.isSingle())
            else:
                self.assertFalse(bond.isSingle())
    
    def testIsDouble(self):
        """
        Test the Bond.isDouble() method.
        """
        for order in self.orderList:
            bond = Bond(order=order)
            if order == 'D':
                self.assertTrue(bond.isDouble())
            else:
                self.assertFalse(bond.isDouble())

    def testIsTriple(self):
        """
        Test the Bond.isTriple() method.
        """
        for order in self.orderList:
            bond = Bond(order=order)
            if order == 'T':
                self.assertTrue(bond.isTriple())
            else:
                self.assertFalse(bond.isTriple())

    def testIsBenzene(self):
        """
        Test the Bond.isBenzene() method.
        """
        for order in self.orderList:
            bond = Bond(order=order)
            if order == 'B':
                self.assertTrue(bond.isBenzene())
            else:
                self.assertFalse(bond.isBenzene())

    def testIncrementOrder(self):
        """
        Test the Bond.incrementOrder() method.
        """
        for order in self.orderList:
            bond = Bond(order=order)
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
            bond = Bond(order=order)
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
            bond0 = Bond(order=order0)
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
            bond0 = Bond(order=order0)
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
            bond0 = Bond(order=order0)
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
            bond0 = Bond(order=order0)
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
            bond0 = Bond(order=order0)
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
            bond0 = Bond(order=order0)
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
                bond1 = Bond(order=order1)
                bond2 = Bond(order=order2)
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
                bond1 = Bond(order=order1)
                bond2 = Bond(order=order2)
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
        
    def testOutput(self):
        """
        Test that we can reconstruct a Bond object from its repr()
        output with no loss of information.
        """
        exec('bond = {0!r}'.format(self.bond))
        self.assertEqual(self.bond.order, bond.order)

################################################################################

class TestMolecule(unittest.TestCase):
    """
    Contains unit tests of the Molecule class.
    """
    
    def setUp(self):
        self.adjlist = """
1  *2 C     1 {2,D} {3,S}
2  *1 O     0 {1,D}
3     C     0 {1,S}
            """
        self.molecule = Molecule().fromAdjacencyList(self.adjlist)
        
    def testClearLabeledAtoms(self):
        """
        Test the Molecule.clearLabeledAtoms() method.
        """
        self.molecule.clearLabeledAtoms()
        for atom in self.molecule.atoms:
            self.assertEqual(atom.label, '')

    def testContainsLabeledAtom(self):
        """
        Test the Molecule.containsLabeledAtom() method.
        """
        for atom in self.molecule.atoms:
            if atom.label != '':
                self.assertTrue(self.molecule.containsLabeledAtom(atom.label))
        self.assertFalse(self.molecule.containsLabeledAtom('*3'))
        self.assertFalse(self.molecule.containsLabeledAtom('*4'))
        self.assertFalse(self.molecule.containsLabeledAtom('*5'))
        self.assertFalse(self.molecule.containsLabeledAtom('*6'))
        
    def testGetLabeledAtom(self):
        """
        Test the Molecule.getLabeledAtom() method.
        """
        for atom in self.molecule.atoms:
            if atom.label != '':
                self.assertEqual(atom, self.molecule.getLabeledAtom(atom.label))
        try:
            self.molecule.getLabeledAtom('*3')
            self.fail('Unexpected successful return from Molecule.getLabeledAtom() with invalid atom label.')
        except ValueError:
            pass
            
    def testGetLabeledAtoms(self):
        """
        Test the Molecule.getLabeledAtoms() method.
        """
        labeled = self.molecule.getLabeledAtoms()
        for atom in self.molecule.atoms:
            if atom.label != '':
                self.assertTrue(atom.label in labeled)
                self.assertTrue(atom in labeled.values())
            else:
                self.assertFalse(atom.label in labeled)
                self.assertFalse(atom in labeled.values())

    def testGetFormula(self):
        """
        Test the Molecule.getLabeledAtoms() method.
        """
        self.assertEqual(self.molecule.getFormula(), 'C2H3O')

    def testGetMolecularWeight(self):
        """
        Test the Molecule.getMolecularWeight() method.
        """
        self.assertAlmostEqual(self.molecule.getMolecularWeight() * 1000, 43.04, 2)

    def testFromAdjacencyList(self):
        """
        Test the Molecule.fromAdjacencyList() method.
        """
        atom1, atom2, atom3 = self.molecule.atoms
        self.assertTrue(self.molecule.hasBond(atom1,atom2))
        self.assertTrue(self.molecule.hasBond(atom1,atom3))
        self.assertFalse(self.molecule.hasBond(atom2,atom3))
        bond12 = self.molecule.bonds[atom1][atom2]
        bond13 = self.molecule.bonds[atom1][atom3]
           
        self.assertTrue(atom1.label == '*2')
        self.assertTrue(atom1.element.symbol == 'C')
        self.assertTrue(atom1.radicalElectrons == 1)
        self.assertTrue(atom1.spinMultiplicity == 2)
        self.assertTrue(atom1.implicitHydrogens == 0)
        self.assertTrue(atom1.charge == 0)
        
        self.assertTrue(atom2.label == '*1')
        self.assertTrue(atom2.element.symbol == 'O')
        self.assertTrue(atom2.radicalElectrons == 0)
        self.assertTrue(atom2.spinMultiplicity == 1)
        self.assertTrue(atom2.implicitHydrogens == 0)
        self.assertTrue(atom2.charge == 0)
        
        self.assertTrue(atom3.label == '')
        self.assertTrue(atom3.element.symbol == 'C')
        self.assertTrue(atom3.radicalElectrons == 0)
        self.assertTrue(atom3.spinMultiplicity == 1)
        self.assertTrue(atom3.implicitHydrogens == 3)
        self.assertTrue(atom3.charge == 0)

        self.assertTrue(bond12.isDouble())
        self.assertTrue(bond13.isSingle())

    def testToAdjacencyList(self):
        """
        Test the Molecule.toAdjacencyList() method.
        """
        adjlist = self.molecule.toAdjacencyList(removeH=True)
        self.assertEqual(adjlist.strip(), self.adjlist.strip())

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
        1 Cd 0 {2,D}
        2 Cd 0 {1,D}
        """)

        self.assertTrue(molecule.isSubgraphIsomorphic(group))
        match, mapping = molecule.findSubgraphIsomorphisms(group)
        self.assertTrue(match)
        self.assertTrue(len(mapping) == 4, "len(mapping) = %d, should be = 4" % (len(mapping)))
        for map in mapping:
            self.assertTrue(len(map) == min(len(molecule.atoms), len(group.atoms)))
            for key, value in map.iteritems():
                self.assertTrue(key in molecule.atoms)
                self.assertTrue(value in group.atoms)

    def testSubgraphIsomorphismAgain(self):
        molecule = Molecule()
        molecule.fromAdjacencyList("""
        1 * C 0 {2,D} {7,S} {8,S}
        2 C 0 {1,D} {3,S} {9,S}
        3 C 0 {2,S} {4,D} {10,S}
        4 C 0 {3,D} {5,S} {11,S}
        5 C 0 {4,S} {6,S} {12,S} {13,S}
        6 C 0 {5,S} {14,S} {15,S} {16,S}
        7 H 0 {1,S}
        8 H 0 {1,S}
        9 H 0 {2,S}
        10 H 0 {3,S}
        11 H 0 {4,S}
        12 H 0 {5,S}
        13 H 0 {5,S}
        14 H 0 {6,S}
        15 H 0 {6,S}
        16 H 0 {6,S}
        """)

        group = Group()
        group.fromAdjacencyList("""
        1 * C 0 {2,D} {3,S} {4,S}
        2   C 0 {1,D}
        3   H 0 {1,S}
        4   H 0 {1,S}
        """)

        molecule.makeHydrogensExplicit()

        labeled1 = molecule.getLabeledAtoms().values()[0]
        labeled2 = group.getLabeledAtoms().values()[0]

        initialMap = {labeled1: labeled2}
        self.assertTrue(molecule.isSubgraphIsomorphic(group, initialMap))

        initialMap = {labeled1: labeled2}
        match, mapping = molecule.findSubgraphIsomorphisms(group, initialMap)
        self.assertTrue(match)
        self.assertTrue(len(mapping) == 2,  "len(mapping) = %d, should be = 2" % (len(mapping)))
        for map in mapping:
            self.assertTrue(len(map) == min(len(molecule.atoms), len(group.atoms)))
            for key, value in map.iteritems():
                self.assertTrue(key in molecule.atoms)
                self.assertTrue(value in group.atoms)

    def testSubgraphIsomorphismManyLabels(self):
        molecule = Molecule() # specific case (species)
        molecule.fromAdjacencyList("""
1 *1 C  1 {2,S} {3,S}
2    C  0 {1,S} {3,S}
3    C  0 {1,S} {2,S}
        """)

        group = Group() # general case (functional group)
        group.fromAdjacencyList("""
1 *1 C 1 {2,S}, {3,S}
2    R 0 {1,S}
3    R 0 {1,S}
        """)

        labeled1 = molecule.getLabeledAtoms()
        labeled2 = group.getLabeledAtoms()
        initialMap = {}
        for label,atom1 in labeled1.iteritems():
            initialMap[atom1] = labeled2[label]
        self.assertTrue(molecule.isSubgraphIsomorphic(group, initialMap))

        match, mapping = molecule.findSubgraphIsomorphisms(group, initialMap)
        self.assertTrue(match)
        self.assertTrue(len(mapping) == 1)
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
        1 C 0       {2,D}
        2 C 0 {1,D} {3,S}
        3 C 0 {2,S} {4,D}
        4 C 0 {3,D} {5,S}
        5 C 1 {4,S} {6,S}
        6 C 0 {5,S}
        """)
        molecule2 = Molecule().fromSMILES('C=CC=C[CH]C')
        
        molecule1.makeHydrogensExplicit()
        molecule2.makeHydrogensExplicit()
        self.assertTrue(molecule1.isIsomorphic(molecule2))
        self.assertTrue(molecule2.isIsomorphic(molecule1))
        
        molecule1.makeHydrogensImplicit()
        molecule2.makeHydrogensImplicit()
        self.assertTrue(molecule1.isIsomorphic(molecule2))
        self.assertTrue(molecule2.isIsomorphic(molecule1))
        
        molecule1.makeHydrogensExplicit()
        molecule2.makeHydrogensImplicit()
        self.assertTrue(molecule1.isIsomorphic(molecule2))
        self.assertTrue(molecule2.isIsomorphic(molecule1))
        
        molecule1.makeHydrogensImplicit()
        molecule2.makeHydrogensExplicit()
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
        for atom1 in molecule.bonds:
            for atom2 in molecule.bonds[atom1]:
                self.assertFalse(molecule.isBondInCycle(atom1, atom2))

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
        for atom1 in molecule.bonds:
            for atom2 in molecule.bonds[atom1]:
                if atom1.isCarbon() and atom2.isCarbon():
                    self.assertTrue(molecule.isBondInCycle(atom1, atom2))
                else:
                    self.assertFalse(molecule.isBondInCycle(atom1, atom2))

    def testImplicitHydrogens(self):
        """
        Test that a molecule can be converted to and from implicit hydrogen
        mode with no loss of data.
        """
        self.assertTrue(self.molecule.implicitHydrogens)
        self.assertEqual(len(self.molecule.atoms), 3)
        self.assertEqual(self.molecule.atoms[0].implicitHydrogens, 0)
        self.assertEqual(self.molecule.atoms[1].implicitHydrogens, 0)
        self.assertEqual(self.molecule.atoms[2].implicitHydrogens, 3)
        self.assertEqual(sum([1 for atom in self.molecule.atoms if atom.isHydrogen()]), 0)
        self.assertEqual(self.molecule.getFormula(), 'C2H3O')
        self.assertAlmostEqual(self.molecule.getMolecularWeight() * 1000, 43.04, 2)
        
        self.molecule.makeHydrogensExplicit()
        self.assertFalse(self.molecule.implicitHydrogens)
        self.assertEqual(len(self.molecule.atoms), 6)
        self.assertEqual(self.molecule.atoms[0].implicitHydrogens, 0)
        self.assertEqual(self.molecule.atoms[1].implicitHydrogens, 0)
        self.assertEqual(self.molecule.atoms[2].implicitHydrogens, 0)
        self.assertEqual(sum([1 for atom in self.molecule.atoms if atom.isHydrogen()]), 3)
        self.assertEqual(self.molecule.getFormula(), 'C2H3O')
        self.assertAlmostEqual(self.molecule.getMolecularWeight() * 1000, 43.04, 2)
        
        self.molecule.makeHydrogensImplicit()
        self.assertTrue(self.molecule.implicitHydrogens)
        self.assertEqual(len(self.molecule.atoms), 3)
        self.assertEqual(self.molecule.atoms[0].implicitHydrogens, 3)
        self.assertEqual(self.molecule.atoms[1].implicitHydrogens, 0)
        self.assertEqual(self.molecule.atoms[2].implicitHydrogens, 0)
        self.assertEqual(sum([1 for atom in self.molecule.atoms if atom.isHydrogen()]), 0)
        self.assertEqual(self.molecule.getFormula(), 'C2H3O')
        self.assertAlmostEqual(self.molecule.getMolecularWeight() * 1000, 43.04, 2)
        
    def testFromSMILESH(self):
        """
        Make sure that H radical is produced properly from its SMILES
        representation.
        """
        molecule = Molecule(SMILES='[H]')
        self.assertTrue(len(molecule.atoms) == 1)
        H = molecule.atoms[0]
        self.assertTrue(H.isHydrogen())
        self.assertTrue(H.radicalElectrons == 1)

    def testFromInChIH(self):
        """
        Make sure that H radical is produced properly from its InChI
        representation.
        """
        molecule = Molecule(InChI='InChI=1/H')
        self.assertTrue(len(molecule.atoms) == 1)
        H = molecule.atoms[0]
        self.assertTrue(H.isHydrogen())
        self.assertTrue(H.radicalElectrons == 1)

    def testPickle(self):
        """
        Test that a Molecule object can be successfully pickled and
        unpickled with no loss of information.
        """
        molecule0 = Molecule().fromSMILES('C=CC=C[CH2]C')
        molecule0.updateAtomTypes()
        molecule0.updateConnectivityValues()
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
        self.assertEqual(molecule.atoms[0].spinMultiplicity, 4)

    def testRadicalCH2(self):
        """
        Test that the species [CH2] has two radical electrons and a spin multiplicity of 3.
        """
        molecule = Molecule().fromSMILES('[CH2]')
        self.assertEqual(molecule.atoms[0].radicalElectrons, 2)
        self.assertEqual(molecule.atoms[0].spinMultiplicity, 3)
        
    def testInChIKey(self):
        """
        Test that InChI Key generation is working properly.
        """
        molecule = Molecule().fromInChI('InChI=1S/C7H12/c1-2-7-4-3-6(1)5-7/h6-7H,1-5H2')
        key = molecule.toInChIKey()
        self.assertEqual(key, 'UMRZSTCPUPJPOJ-UHFFFAOYSA')
################################################################################

class TestMoleculeSymmetry(unittest.TestCase):
    """
    Contains unit tests of the methods for computing symmetry numbers for a
    given Molecule object.
    """
    
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
    
    def testCountInternalRotorsDimethylAcetylene(self):
        """
        Test the Molecule.countInternalRotors() method for dimethylacetylene.
        This is a "hard" test that currently fails.
        """
        self.assertEqual(Molecule().fromSMILES('CC#CC').countInternalRotors(), 1)
    
    def testAtomSymmetryNumberMethane(self):
        """
        Test the Molecule.calculateAtomSymmtryNumber() method.
        """
        molecule = Molecule().fromSMILES('C')
        symmetryNumber = 1
        for atom in molecule.atoms:
            if not molecule.isAtomInCycle(atom):
                symmetryNumber *= molecule.calculateAtomSymmetryNumber(atom)
        self.assertEqual(symmetryNumber, 12)
        
    def testAtomSymmetryNumberMethyl(self):
        """
        Test the Molecule.calculateAtomSymmtryNumber() method.
        """
        molecule = Molecule().fromSMILES('[CH3]')
        symmetryNumber = 1
        for atom in molecule.atoms:
            if not molecule.isAtomInCycle(atom):
                symmetryNumber *= molecule.calculateAtomSymmetryNumber(atom)
        self.assertEqual(symmetryNumber, 6)
        
    def testAtomSymmetryNumberEthane(self):
        """
        Test the Molecule.calculateAtomSymmtryNumber() method.
        """
        molecule = Molecule().fromSMILES('CC')
        symmetryNumber = 1
        for atom in molecule.atoms:
            if not molecule.isAtomInCycle(atom):
                symmetryNumber *= molecule.calculateAtomSymmetryNumber(atom)
        self.assertEqual(symmetryNumber, 9)
    
    def testAtomSymmetryNumberPropane(self):
        """
        Test the Molecule.calculateAtomSymmtryNumber() method.
        """
        molecule = Molecule().fromSMILES('CCC')
        symmetryNumber = 1
        for atom in molecule.atoms:
            if not molecule.isAtomInCycle(atom):
                symmetryNumber *= molecule.calculateAtomSymmetryNumber(atom)
        self.assertEqual(symmetryNumber, 18)
    
    def testAtomSymmetryNumberIsobutane(self):
        """
        Test the Molecule.calculateAtomSymmtryNumber() method.
        """
        molecule = Molecule().fromSMILES('CC(C)C')
        symmetryNumber = 1
        for atom in molecule.atoms:
            if not molecule.isAtomInCycle(atom):
                symmetryNumber *= molecule.calculateAtomSymmetryNumber(atom)
        self.assertEqual(symmetryNumber, 81)

    def testBondSymmetryNumberEthane(self):
        """
        Test the Molecule.calculateBondSymmtryNumber() method.
        """
        molecule = Molecule().fromSMILES('CC')
        symmetryNumber = 1
        for atom1 in molecule.bonds:
            for atom2 in molecule.bonds[atom1]:
                if molecule.atoms.index(atom1) < molecule.atoms.index(atom2):
                    symmetryNumber *= molecule.calculateBondSymmetryNumber(atom1, atom2)
        self.assertEqual(symmetryNumber, 2)
        
    def testBondSymmetryNumberPropane(self):
        """
        Test the Molecule.calculateBondSymmtryNumber() method.
        """
        molecule = Molecule().fromSMILES('CCC')
        symmetryNumber = 1
        for atom1 in molecule.bonds:
            for atom2 in molecule.bonds[atom1]:
                if molecule.atoms.index(atom1) < molecule.atoms.index(atom2):
                    symmetryNumber *= molecule.calculateBondSymmetryNumber(atom1, atom2)
        self.assertEqual(symmetryNumber, 1)
    
    def testBondSymmetryNumberButane(self):
        """
        Test the Molecule.calculateBondSymmtryNumber() method.
        """
        molecule = Molecule().fromSMILES('CCCC')
        symmetryNumber = 1
        for atom1 in molecule.bonds:
            for atom2 in molecule.bonds[atom1]:
                if molecule.atoms.index(atom1) < molecule.atoms.index(atom2):
                    symmetryNumber *= molecule.calculateBondSymmetryNumber(atom1, atom2)
        self.assertEqual(symmetryNumber, 2)
    
    def testBondSymmetryNumberEthylene(self):
        """
        Test the Molecule.calculateBondSymmtryNumber() method.
        """
        molecule = Molecule().fromSMILES('C=C')
        symmetryNumber = 1
        for atom1 in molecule.bonds:
            for atom2 in molecule.bonds[atom1]:
                if molecule.atoms.index(atom1) < molecule.atoms.index(atom2):
                    symmetryNumber *= molecule.calculateBondSymmetryNumber(atom1, atom2)
        self.assertEqual(symmetryNumber, 2)
    
    def testBondSymmetryNumberAcetylene(self):
        """
        Test the Molecule.calculateBondSymmtryNumber() method.
        """
        molecule = Molecule().fromSMILES('C#C')
        symmetryNumber = 1
        for atom1 in molecule.bonds:
            for atom2 in molecule.bonds[atom1]:
                if molecule.atoms.index(atom1) < molecule.atoms.index(atom2):
                    symmetryNumber *= molecule.calculateBondSymmetryNumber(atom1, atom2)
        self.assertEqual(symmetryNumber, 2)
    
    def testAxisSymmetryNumberEthylene(self):
        """
        Test the Molecule.calculateAxisSymmtryNumber() method.
        """
        self.assertEqual(Molecule().fromSMILES('C=C').calculateAxisSymmetryNumber(), 2)
        
    def testAxisSymmetryNumberPropadiene(self):
        """
        Test the Molecule.calculateAxisSymmtryNumber() method.
        """
        self.assertEqual(Molecule().fromSMILES('C=C=C').calculateAxisSymmetryNumber(), 2)
        
    def testAxisSymmetryNumberButatriene(self):
        """
        Test the Molecule.calculateAxisSymmtryNumber() method.
        """
        self.assertEqual(Molecule().fromSMILES('C=C=C=C').calculateAxisSymmetryNumber(), 2)
        
    def testAxisSymmetryNumberButatrienyl(self):
        """
        Test the Molecule.calculateAxisSymmtryNumber() method.
        """
        self.assertEqual(Molecule().fromSMILES('C=C=C=[CH]').calculateAxisSymmetryNumber(), 1)
    
    def testAxisSymmetryNumberPropadienyl(self):
        """
        Test the Molecule.calculateAxisSymmtryNumber() method.
        """
        self.assertEqual(Molecule().fromSMILES('C=C=[C]').calculateAxisSymmetryNumber(), 2)
        
    def testAxisSymmetryNumber12Butadienyl(self):
        """
        Test the Molecule.calculateAxisSymmtryNumber() method.
        """
        self.assertEqual(Molecule().fromSMILES('CC=C=[C]').calculateAxisSymmetryNumber(), 1)
    
    def testAxisSymmetryNumber12Hexadienyl(self):
        """
        Test the Molecule.calculateAxisSymmtryNumber() method.
        """
        self.assertEqual(Molecule().fromSMILES('C=C=CCCC').calculateAxisSymmetryNumber(), 1)
    
    def testAxisSymmetryNumber1(self):
        """
        Test the Molecule.calculateAxisSymmtryNumber() method.
        """
        self.assertEqual(Molecule().fromSMILES('CC(C)=C=C(CC)CC').calculateAxisSymmetryNumber(), 2)
    
    def testAxisSymmetryNumber2(self):
        """
        Test the Molecule.calculateAxisSymmtryNumber() method.
        """
        self.assertEqual(Molecule().fromSMILES('C=C=C(C(C(C(C=C=C)=C=C)=C=C)=C=C)').calculateAxisSymmetryNumber(), 2)
    
    def testAxisSymmetryNumber3(self):
        """
        Test the Molecule.calculateAxisSymmtryNumber() method.
        """
        self.assertEqual(Molecule().fromSMILES('C=C=[C]C(C)(C)[C]=C=C').calculateAxisSymmetryNumber(), 4)
    
    def testAxisSymmetryNumber4(self):
        """
        Test the Molecule.calculateAxisSymmtryNumber() method.
        """
        self.assertEqual(Molecule().fromSMILES('C=C=C=O').calculateAxisSymmetryNumber(), 2)
    
    def testAxisSymmetryNumber5(self):
        """
        Test the Molecule.calculateAxisSymmtryNumber() method.
        """
        self.assertEqual(Molecule().fromSMILES('CC=C=C=O').calculateAxisSymmetryNumber(), 1)
    
    def testAxisSymmetryNumber5(self):
        """
        Test the Molecule.calculateAxisSymmtryNumber() method.
        """
        self.assertEqual(Molecule().fromSMILES('C=C=C=N').calculateAxisSymmetryNumber(), 1)
    
    def testAxisSymmetryNumber5(self):
        """
        Test the Molecule.calculateAxisSymmtryNumber() method.
        """
        self.assertEqual(Molecule().fromSMILES('C=C=C=[N]').calculateAxisSymmetryNumber(), 2)
    
#   def testCyclicSymmetryNumber(self):
#
#		# cyclohexane
#       molecule = Molecule().fromInChI('InChI=1/C6H12/c1-2-4-6-5-3-1/h1-6H2')
#       molecule.makeHydrogensExplicit()
#       symmetryNumber = molecule.calculateCyclicSymmetryNumber()
#       self.assertEqual(symmetryNumber, 12)

    def testTotalSymmetryNumberEthane(self):
        """
        Test the Molecule.calculateSymmtryNumber() method.
        """
        self.assertEqual(Molecule().fromSMILES('CC').calculateSymmetryNumber(), 18)
    
    def testTotalSymmetryNumber1(self):
        """
        Test the Molecule.calculateSymmtryNumber() method.
        """
        self.assertEqual(Molecule().fromSMILES('C=C=[C]C(C)(C)[C]=C=C').calculateSymmetryNumber(), '???')
    
    def testTotalSymmetryNumber2(self):
        """
        Test the Molecule.calculateSymmtryNumber() method.
        """
        self.assertEqual(Molecule().fromSMILES('C(=CC(c1ccccc1)C([CH]CCCCCC)C=Cc1ccccc1)[CH]CCCCCC').calculateSymmetryNumber(), 1)
    
    def testSymmetryNumberHydroxyl(self):
        """
        Test the Molecule.calculateSymmtryNumber() method.
        """
        self.assertEqual(Molecule().fromSMILES('[OH]').calculateSymmetryNumber(), 1)
       
    def testSymmetryNumberOxygen(self):
        """
        Test the Molecule.calculateSymmtryNumber() method.
        """
        self.assertEqual(Molecule().fromSMILES('O=O').calculateAxisSymmetryNumber(), 1)
        self.assertEqual(Molecule().fromSMILES('O=O').calculateSymmetryNumber(), 2)
        
    def testSymmetryNumberDicarbon(self):
        """
        Test the Molecule.calculateSymmtryNumber() method.
        """
        self.assertEqual(Molecule().fromSMILES('[C]#[C]').calculateSymmetryNumber(), 2)
    
    def testSymmetryNumberHydrogen(self):
        """
        Test the Molecule.calculateSymmtryNumber() method.
        """
        self.assertEqual(Molecule().fromSMILES('[H][H]').calculateSymmetryNumber(), 2)
    
    def testSymmetryNumberAcetylene(self):
        """
        Test the Molecule.calculateSymmtryNumber() method.
        """
        self.assertEqual(Molecule().fromSMILES('C#C').calculateSymmetryNumber(), 2)
    
    def testSymmetryNumberButadiyne(self):
        """
        Test the Molecule.calculateSymmtryNumber() method.
        """
        self.assertEqual(Molecule().fromSMILES('C#CC#C').calculateSymmetryNumber(), 2)
    
    def testSymmetryNumberMethane(self):
        """
        Test the Molecule.calculateSymmtryNumber() method.
        """
        self.assertEqual(Molecule().fromSMILES('C').calculateSymmetryNumber(), 12)
    
    def testSymmetryNumberFormaldehyde(self):
        """
        Test the Molecule.calculateSymmtryNumber() method.
        """
        self.assertEqual(Molecule().fromSMILES('C=O').calculateSymmetryNumber(), 2)
    
    def testSymmetryNumberMethyl(self):
        """
        Test the Molecule.calculateSymmtryNumber() method.
        """
        self.assertEqual(Molecule().fromSMILES('[CH3]').calculateSymmetryNumber(), 6)
    
    def testSymmetryNumberWater(self):
        """
        Test the Molecule.calculateSymmtryNumber() method.
        """
        self.assertEqual(Molecule().fromSMILES('O').calculateSymmetryNumber(), 2)
    
    def testSymmetryNumberEthylene(self):
        """
        Test the Molecule.calculateSymmtryNumber() method.
        """
        self.assertEqual(Molecule().fromSMILES('C=C').calculateSymmetryNumber(), 4)
    
    def testSymmetryNumberEthenyl(self):
        """
        Test the Molecule.calculateSymmtryNumber() method.
        """
        self.assertEqual(Molecule().fromSMILES('C=[CH]').calculateSymmetryNumber(), 1)
    
    def testSymmetryNumberCyclic(self):
        """
        Test the Molecule.calculateSymmtryNumber() method.
        """
        self.assertEqual(Molecule().fromSMILES('C1=C=C=1').calculateSymmetryNumber(), '6?')
    
################################################################################

if __name__ == '__main__':
    unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )