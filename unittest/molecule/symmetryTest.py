#!/usr/bin/python
# -*- coding: utf-8 -*-

import unittest

from rmgpy.molecule.molecule import *
from rmgpy.molecule.symmetry import *

################################################################################

class TestMoleculeSymmetry(unittest.TestCase):
    """
    Contains unit tests of the methods for computing symmetry numbers for a
    given Molecule object.
    """
        
    def testAtomSymmetryNumberMethane(self):
        """
        Test the Molecule.calculateAtomSymmtryNumber() method.
        """
        molecule = Molecule().fromSMILES('C')
        symmetryNumber = 1
        for atom in molecule.atoms:
            if not molecule.isAtomInCycle(atom):
                symmetryNumber *= calculateAtomSymmetryNumber(molecule, atom)
        self.assertEqual(symmetryNumber, 12)
        
    def testAtomSymmetryNumberMethyl(self):
        """
        Test the Molecule.calculateAtomSymmtryNumber() method.
        """
        molecule = Molecule().fromSMILES('[CH3]')
        symmetryNumber = 1
        for atom in molecule.atoms:
            if not molecule.isAtomInCycle(atom):
                symmetryNumber *= calculateAtomSymmetryNumber(molecule, atom)
        self.assertEqual(symmetryNumber, 6)
        
    def testAtomSymmetryNumberEthane(self):
        """
        Test the Molecule.calculateAtomSymmtryNumber() method.
        """
        molecule = Molecule().fromSMILES('CC')
        symmetryNumber = 1
        for atom in molecule.atoms:
            if not molecule.isAtomInCycle(atom):
                symmetryNumber *= calculateAtomSymmetryNumber(molecule, atom)
        self.assertEqual(symmetryNumber, 9)
    
    def testAtomSymmetryNumberPropane(self):
        """
        Test the Molecule.calculateAtomSymmtryNumber() method.
        """
        molecule = Molecule().fromSMILES('CCC')
        symmetryNumber = 1
        for atom in molecule.atoms:
            if not molecule.isAtomInCycle(atom):
                symmetryNumber *= calculateAtomSymmetryNumber(molecule, atom)
        self.assertEqual(symmetryNumber, 18)
    
    def testAtomSymmetryNumberIsobutane(self):
        """
        Test the Molecule.calculateAtomSymmtryNumber() method.
        """
        molecule = Molecule().fromSMILES('CC(C)C')
        symmetryNumber = 1
        for atom in molecule.atoms:
            if not molecule.isAtomInCycle(atom):
                symmetryNumber *= calculateAtomSymmetryNumber(molecule, atom)
        self.assertEqual(symmetryNumber, 81)

    def testBondSymmetryNumberEthane(self):
        """
        Test the Molecule.calculateBondSymmtryNumber() method.
        """
        molecule = Molecule().fromSMILES('CC')
        symmetryNumber = 1
        for atom1 in molecule.atoms:
            for atom2 in atom1.bonds:
                if molecule.atoms.index(atom1) < molecule.atoms.index(atom2):
                    symmetryNumber *= calculateBondSymmetryNumber(molecule, atom1, atom2)
        self.assertEqual(symmetryNumber, 2)
        
    def testBondSymmetryNumberPropane(self):
        """
        Test the Molecule.calculateBondSymmtryNumber() method.
        """
        molecule = Molecule().fromSMILES('CCC')
        symmetryNumber = 1
        for atom1 in molecule.atoms:
            for atom2 in atom1.bonds:
                if molecule.atoms.index(atom1) < molecule.atoms.index(atom2):
                    symmetryNumber *= calculateBondSymmetryNumber(molecule, atom1, atom2)
        self.assertEqual(symmetryNumber, 1)
    
    def testBondSymmetryNumberButane(self):
        """
        Test the Molecule.calculateBondSymmtryNumber() method.
        """
        molecule = Molecule().fromSMILES('CCCC')
        symmetryNumber = 1
        for atom1 in molecule.atoms:
            for atom2 in atom1.bonds:
                if molecule.atoms.index(atom1) < molecule.atoms.index(atom2):
                    symmetryNumber *= calculateBondSymmetryNumber(molecule, atom1, atom2)
        self.assertEqual(symmetryNumber, 2)
    
    def testBondSymmetryNumberEthylene(self):
        """
        Test the Molecule.calculateBondSymmtryNumber() method.
        """
        molecule = Molecule().fromSMILES('C=C')
        symmetryNumber = 1
        for atom1 in molecule.atoms:
            for atom2 in atom1.bonds:
                if molecule.atoms.index(atom1) < molecule.atoms.index(atom2):
                    symmetryNumber *= calculateBondSymmetryNumber(molecule, atom1, atom2)
        self.assertEqual(symmetryNumber, 2)
    
    def testBondSymmetryNumberAcetylene(self):
        """
        Test the Molecule.calculateBondSymmtryNumber() method.
        """
        molecule = Molecule().fromSMILES('C#C')
        symmetryNumber = 1
        for atom1 in molecule.atoms:
            for atom2 in atom1.bonds:
                if molecule.atoms.index(atom1) < molecule.atoms.index(atom2):
                    symmetryNumber *= calculateBondSymmetryNumber(molecule, atom1, atom2)
        self.assertEqual(symmetryNumber, 2)
    
    def testAxisSymmetryNumberEthylene(self):
        """
        Test the Molecule.calculateAxisSymmtryNumber() method.
        """
        molecule = Molecule().fromSMILES('C=C')
        self.assertEqual(calculateAxisSymmetryNumber(molecule), 2)
        
    def testAxisSymmetryNumberPropadiene(self):
        """
        Test the Molecule.calculateAxisSymmtryNumber() method.
        """
        molecule = Molecule().fromSMILES('C=C=C')
        self.assertEqual(calculateAxisSymmetryNumber(molecule), 2)
        
    def testAxisSymmetryNumberButatriene(self):
        """
        Test the Molecule.calculateAxisSymmtryNumber() method.
        """
        molecule = Molecule().fromSMILES('C=C=C=C')
        self.assertEqual(calculateAxisSymmetryNumber(molecule), 2)
        
    def testAxisSymmetryNumberButatrienyl(self):
        """
        Test the Molecule.calculateAxisSymmtryNumber() method.
        """
        molecule = Molecule().fromSMILES('C=C=C=[CH]')
        self.assertEqual(calculateAxisSymmetryNumber(molecule), 1)
    
    def testAxisSymmetryNumberPropadienyl(self):
        """
        Test the Molecule.calculateAxisSymmtryNumber() method.
        """
        molecule = Molecule().fromSMILES('C=C=[C]')
        self.assertEqual(calculateAxisSymmetryNumber(molecule), 2)
        
    def testAxisSymmetryNumber12Butadienyl(self):
        """
        Test the Molecule.calculateAxisSymmtryNumber() method.
        """
        molecule = Molecule().fromSMILES('CC=C=[C]')
        self.assertEqual(calculateAxisSymmetryNumber(molecule), 1)
    
    def testAxisSymmetryNumber12Hexadienyl(self):
        """
        Test the Molecule.calculateAxisSymmtryNumber() method.
        """
        molecule = Molecule().fromSMILES('C=C=CCCC')
        self.assertEqual(calculateAxisSymmetryNumber(molecule), 1)
    
    def testAxisSymmetryNumber1(self):
        """
        Test the Molecule.calculateAxisSymmtryNumber() method.
        """
        molecule = Molecule().fromSMILES('CC(C)=C=C(CC)CC')
        self.assertEqual(calculateAxisSymmetryNumber(molecule), 2)
    
    def testAxisSymmetryNumber2(self):
        """
        Test the Molecule.calculateAxisSymmtryNumber() method.
        """
        molecule = Molecule().fromSMILES('C=C=C(C(C(C(C=C=C)=C=C)=C=C)=C=C)')
        self.assertEqual(calculateAxisSymmetryNumber(molecule), 2)
    
    def testAxisSymmetryNumber3(self):
        """
        Test the Molecule.calculateAxisSymmtryNumber() method.
        """
        molecule = Molecule().fromSMILES('C=C=[C]C(C)(C)[C]=C=C')
        self.assertEqual(calculateAxisSymmetryNumber(molecule), 4)
    
    def testAxisSymmetryNumber4(self):
        """
        Test the Molecule.calculateAxisSymmtryNumber() method.
        """
        molecule = Molecule().fromSMILES('C=C=C=O')
        self.assertEqual(calculateAxisSymmetryNumber(molecule), 2)
    
    def testAxisSymmetryNumber5(self):
        """
        Test the Molecule.calculateAxisSymmtryNumber() method.
        """
        molecule = Molecule().fromSMILES('CC=C=C=O')
        self.assertEqual(calculateAxisSymmetryNumber(molecule), 1)
    
    def testAxisSymmetryNumber6(self):
        """
        Test the Molecule.calculateAxisSymmtryNumber() method.
        """
        molecule = Molecule().fromSMILES('C=C=C=N')
        self.assertEqual(calculateAxisSymmetryNumber(molecule), 1)
    
    def testAxisSymmetryNumber7(self):
        """
        Test the Molecule.calculateAxisSymmtryNumber() method.
        """
        molecule = Molecule().fromSMILES('C=C=C=[N]')
        self.assertEqual(calculateAxisSymmetryNumber(molecule), 2)
    
    def testAxisSymmetryOxygenSinglet(self):
        """
        Test the Molecule.calculateAxisSymmtryNumber() method.
        """
        molecule = Molecule().fromSMILES('O=O')
        self.assertEqual(calculateAxisSymmetryNumber(molecule), 1)
    
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