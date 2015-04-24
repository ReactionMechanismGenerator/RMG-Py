#!/usr/bin/python
# -*- coding: utf-8 -*-

import unittest
from external.wip import work_in_progress
from rmgpy.molecule.molecule import Atom, Bond, Molecule, ActionError
from rmgpy.molecule.group import Group
from rmgpy.molecule.element import getElement, elementList

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
        molecule = Molecule(InChI='InChI=1/H')
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
        test_strings =['[C-]#[O+]', '[C]', '[CH]', 'OO', '[H][H]', '[H]', '[He]', '[O]', 'O', '[CH3]', 'C', '[OH]', 'CCC', 'CC', 'N#N', '[O]O', 'C[CH2]', '[Ar]', 'CCCC','O=C=O','[C]#N']
        for s in test_strings:
            molecule = Molecule(SMILES=s)
            self.assertEqual(s,molecule.toSMILES())

    def testInChIKey(self):
        """
        Test that InChI Key generation is working properly.
        """
        molecule = Molecule().fromInChI('InChI=1S/C7H12/c1-2-7-4-3-6(1)5-7/h6-7H,1-5H2')
        key = molecule.toInChIKey()
        self.assertEqual(key, 'UMRZSTCPUPJPOJ-UHFFFAOYSA')
        
    def testAugmentedInChI(self):
        """
        Test that Augmented InChI generation is printing the /mult layer
        """
        mol = Molecule().fromAdjacencyList("""
            1     C     u1 p0 c0 {2,S}
            2     C     u1 p0 c0 {1,S}
        """, saturateH=True)
        
        self.assertEqual(mol.toAugmentedInChI(), 'InChI=1S/C2H4/c1-2/h1-2H2/mult3')
        
    def testAugmentedInChIKey(self):
        """
        Test that Augmented InChI Key generation is printing the mult layer
        """
        mol = Molecule().fromAdjacencyList("""
            1     C     u1 p0 c0 {2,S}
            2     C     u1 p0 c0 {1,S}
        """, saturateH=True)
        
        self.assertEqual(mol.toAugmentedInChIKey(), 'VGGSQFUCUMXWEO-UHFFFAOYSA-mult3')

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
    
        self.assertEqual( mol.toAugmentedInChI(), 'InChI=1S/CH2/h1H2/mult3')
        self.assertEqual( mol.toSMILES(), '[CH2]')
        

        ch2_s = '''
        multiplicity 1
        1 C u0 p1 c0 {2,S} {3,S}
        2 H u0 p0 c0 {1,S}
        3 H u0 p0 c0 {1,S}
        '''
        
        mol = Molecule().fromAdjacencyList(ch2_s)
        self.assertEqual( mol.toAugmentedInChI(), 'InChI=1S/CH2/h1H2/mult1')
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
        

################################################################################

if __name__ == '__main__':
    unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )