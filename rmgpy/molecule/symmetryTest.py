#!/usr/bin/python
# -*- coding: utf-8 -*-

import unittest
from external.wip import work_in_progress

from rmgpy.molecule.molecule import Molecule
from rmgpy.molecule.symmetry import calculateAtomSymmetryNumber, calculateAxisSymmetryNumber, calculateBondSymmetryNumber, calculateCyclicSymmetryNumber
from rmgpy.species import Species
from rmgpy.molecule.resonance import generateAromaticResonanceStructures
################################################################################

class TestMoleculeSymmetry(unittest.TestCase):
    """
    Contains unit tests of the methods for computing symmetry numbers for a
    given Molecule object.
    """
        
    def testAtomSymmetryNumberMethane(self):
        """
        Test the Molecule.calculateAtomSymmetryNumber() on CH4
        """
        molecule = Molecule().fromSMILES('C')
        symmetryNumber = 1
        for atom in molecule.atoms:
            if not molecule.isAtomInCycle(atom):
                symmetryNumber *= calculateAtomSymmetryNumber(molecule, atom)
        self.assertEqual(symmetryNumber, 12)
        
    def testAtomSymmetryNumberMethyl(self):
        """
        Test the Molecule.calculateAtomSymmetryNumber() on [CH3]
        """
        molecule = Molecule().fromSMILES('[CH3]')
        symmetryNumber = 1
        for atom in molecule.atoms:
            if not molecule.isAtomInCycle(atom):
                symmetryNumber *= calculateAtomSymmetryNumber(molecule, atom)
        self.assertEqual(symmetryNumber, 6)
        
    def testAtomSymmetryNumberEthane(self):
        """
        Test the Molecule.calculateAtomSymmetryNumber() on CC
        """
        molecule = Molecule().fromSMILES('CC')
        symmetryNumber = 1
        for atom in molecule.atoms:
            if not molecule.isAtomInCycle(atom):
                symmetryNumber *= calculateAtomSymmetryNumber(molecule, atom)
        self.assertEqual(symmetryNumber, 9)
    
    def testAtomSymmetryNumberPropane(self):
        """
        Test the Molecule.calculateAtomSymmetryNumber() on CCC
        """
        molecule = Molecule().fromSMILES('CCC')
        symmetryNumber = 1
        for atom in molecule.atoms:
            if not molecule.isAtomInCycle(atom):
                symmetryNumber *= calculateAtomSymmetryNumber(molecule, atom)
        self.assertEqual(symmetryNumber, 18)
    
    def testAtomSymmetryNumberIsobutane(self):
        """
        Test the Molecule.calculateAtomSymmetryNumber() on CC(C)C
        """
        molecule = Molecule().fromSMILES('CC(C)C')
        symmetryNumber = 1
        for atom in molecule.atoms:
            if not molecule.isAtomInCycle(atom):
                symmetryNumber *= calculateAtomSymmetryNumber(molecule, atom)
        self.assertEqual(symmetryNumber, 81)

    def testBondSymmetryNumberEthane(self):
        """
        Test the Molecule.calculateBondSymmetryNumber() on CC
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
        Test the Molecule.calculateBondSymmetryNumber() on CCC
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
        Test the Molecule.calculateBondSymmetryNumber() on CCCC
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
        Test the Molecule.calculateBondSymmetryNumber() on C=C
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
        Test the Molecule.calculateBondSymmetryNumber() on C#C
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
        Test the Molecule.calculateAxisSymmetryNumber() on C=C
        """
        molecule = Molecule().fromSMILES('C=C')
        self.assertEqual(calculateAxisSymmetryNumber(molecule), 2)
        
    def testAxisSymmetryNumberPropadiene(self):
        """
        Test the Molecule.calculateAxisSymmetryNumber() on C=C=C
        """
        molecule = Molecule().fromSMILES('C=C=C')
        self.assertEqual(calculateAxisSymmetryNumber(molecule), 2)
        
    def testAxisSymmetryNumberButatriene(self):
        """
        Test the Molecule.calculateAxisSymmetryNumber() on C=C=C=C
        """
        molecule = Molecule().fromSMILES('C=C=C=C')
        self.assertEqual(calculateAxisSymmetryNumber(molecule), 2)
        
    def testAxisSymmetryNumberButatrienyl(self):
        """
        Test the Molecule.calculateAxisSymmetryNumber() on C=C=C=[CH]
        """
        molecule = Molecule().fromSMILES('C=C=C=[CH]')
        self.assertEqual(calculateAxisSymmetryNumber(molecule), 1)
    
    def testAxisSymmetryNumberPropadienyl(self):
        """
        Test the Molecule.calculateAxisSymmetryNumber() on C=C=[C]
        """
        molecule = Molecule().fromSMILES('C=C=[C]')
        self.assertEqual(calculateAxisSymmetryNumber(molecule), 2)
        
    def testAxisSymmetryNumber12Butadienyl(self):
        """
        Test the Molecule.calculateAxisSymmetryNumber() on CC=C=[C]
        """
        molecule = Molecule().fromSMILES('CC=C=[C]')
        self.assertEqual(calculateAxisSymmetryNumber(molecule), 1)
    
    def testAxisSymmetryNumber12Hexadienyl(self):
        """
        Test the Molecule.calculateAxisSymmetryNumber() on C=C=CCCC
        """
        molecule = Molecule().fromSMILES('C=C=CCCC')
        self.assertEqual(calculateAxisSymmetryNumber(molecule), 1)
    
    def testAxisSymmetryNumber1(self):
        """
        Test the Molecule.calculateAxisSymmetryNumber() on CC(C)=C=C(CC)CC
        """
        molecule = Molecule().fromSMILES('CC(C)=C=C(CC)CC')
        self.assertEqual(calculateAxisSymmetryNumber(molecule), 2)
    
    def testAxisSymmetryNumber2(self):
        """
        Test the Molecule.calculateAxisSymmetryNumber() on C=C=C(C(C(C(C=C=C)=C=C)=C=C)=C=C)
        """
        molecule = Molecule().fromSMILES('C=C=C(C(C(C(C=C=C)=C=C)=C=C)=C=C)')
        self.assertEqual(calculateAxisSymmetryNumber(molecule), 2)
    
    @work_in_progress
    def testAxisSymmetryNumber3(self):
        """
        Test the Molecule.calculateAxisSymmetryNumber() on C=C=[C]C(C)(C)[C]=C=C
        """
        molecule = Molecule().fromSMILES('C=C=[C]C(C)(C)[C]=C=C')
        self.assertEqual(calculateAxisSymmetryNumber(molecule), 4)
    
    def testAxisSymmetryNumber4(self):
        """
        Test the Molecule.calculateAxisSymmetryNumber() on C=C=C=O
        """
        molecule = Molecule().fromSMILES('C=C=C=O')
        self.assertEqual(calculateAxisSymmetryNumber(molecule), 2)
    
    def testAxisSymmetryNumber5(self):
        """
        Test the Molecule.calculateAxisSymmetryNumber() on CC=C=C=O
        """
        molecule = Molecule().fromSMILES('CC=C=C=O')
        self.assertEqual(calculateAxisSymmetryNumber(molecule), 1)
    
    def testAxisSymmetryNumber6(self):
        """
        Test the Molecule.calculateAxisSymmetryNumber() on C=C=C=N
        """
        molecule = Molecule().fromSMILES('C=C=C=N')
        self.assertEqual(calculateAxisSymmetryNumber(molecule), 1)
    
    def testAxisSymmetryNumber7(self):
        """
        Test the Molecule.calculateAxisSymmetryNumber() on C=C=C=[N]
        """
        molecule = Molecule().fromSMILES('C=C=C=[N]')
        self.assertEqual(calculateAxisSymmetryNumber(molecule), 2)
    
    def testAxisSymmetryOxygenSinglet(self):
        """
        Test the Molecule.calculateAxisSymmetryNumber() on O=O
        """
        molecule = Molecule().fromSMILES('O=O')
        self.assertEqual(calculateAxisSymmetryNumber(molecule), 1)
    
    def testCyclicSymmetryNumberCyclohexane(self):
        """
        Test the Molecule.calculateCyclicSymmetryNumber() on C1CCCCC1
        """
        molecule = Molecule().fromSMILES('C1CCCCC1')
        symmetryNumber = calculateCyclicSymmetryNumber(molecule)
        self.assertEqual(symmetryNumber, 12)

    @work_in_progress
    def testCyclicSymmetryNumberCyclohexanone(self):
        """
        Test the Molecule.calculateCyclicSymmetryNumber() on C1CCCCC1=O
        """
        molecule = Molecule().fromSMILES('C1CCCCC1=O')
        symmetryNumber = calculateCyclicSymmetryNumber(molecule)
        self.assertEqual(symmetryNumber, 2)
        
    @work_in_progress
    def testCyclicSymmetryNumberCyclohexan_tri_one(self):
        """
        Test the Molecule.calculateCyclicSymmetryNumber() on C1CCC(=O)C(=O)C1=O
        """
        molecule = Molecule().fromSMILES('C1CCC(=O)C(=O)C1=O')
        symmetryNumber = calculateCyclicSymmetryNumber(molecule)
        self.assertEqual(symmetryNumber, 2)

    def testTotalSymmetryNumberBenzene(self):
        """
        Test the Species.getSymmetryNumber() (total symmetry) on c1ccccc1
        """
        molecule = Molecule().fromSMILES('c1ccccc1')
        species = Species(molecule=[molecule])
        symmetryNumber = species.getSymmetryNumber()
        self.assertEqual(symmetryNumber, 12)

    def testTotalSymmetryNumberbute_di_yl(self):
        """
        Test the Species.getSymmetryNumber() (total symmetry) on [CH2][CH]C=C
        """
        species = Species().fromSMILES('[CH2][CH]C=C')
        symmetryNumber = species.getSymmetryNumber()
        self.assertEqual(symmetryNumber, 2)

    def testTotalSymmetryNumberAllyl(self):
        """
        Test the Species.getSymmetryNumber() (total symmetry) on Allyl, [CH2]C=C
        """
        molecule = Molecule().fromSMILES('[CH2]C=C')
        species = Species(molecule=[molecule])
        symmetryNumber = species.getSymmetryNumber()
        self.assertEqual(symmetryNumber, 2)

    def testTotalSymmetryNumberPentenyl(self):
        """
        Test the Species.getSymmetryNumer() for C[CH]C=CC
        and ensures that it is differet than the molecule object
        """
        spc = Species(molecule=[Molecule().fromSMILES('C[CH]C=CC')])
        symmetryNumber = spc.getSymmetryNumber()
        self.assertEqual(symmetryNumber, 18, 'the symmetry number for C[CH]C=CC is 18, but RMG returned {}'.format(symmetryNumber))

    def testSpeciesSymmetryNumberIsNotMoleculeSymmetryNumber(self):
        """
        Tests that the species symmetry number can be different from the molecule symmetry number

        This molecule's resonance isomer hybrid should return more symmetry
        than the base molecule object.
        """
        molecule = Molecule().fromSMILES('C[CH]C=CC')
        species = Species(molecule=[molecule])
        self.assertEqual(molecule.getSymmetryNumber() * 2, species.getSymmetryNumber())

    def testAxisSymmetryNumberAllyl(self):
        """
        Test the Molecule.calculateAxisSymmetryNumber() on [CH2]C=C
        """
        spc = Species(molecule=[Molecule().fromSMILES('[CH2]C=C')])
        molecule = spc.getResonanceHybrid()
        self.assertEqual(calculateAxisSymmetryNumber(molecule), 1)

    def testMiddleCarbonAtomSymmetryNumberAllylUsingResonanceHybrid(self):
        """
        Test the Molecule.calculateAtomSymmetryNumber() for the middle carbon
        on [CH2]C=C using the resonance hybrid structure
        """
        molecule = Molecule().fromSMILES('[CH2]C=C')
        species = Species(molecule=[molecule])
        resonanceHybrid = species.getResonanceHybrid()
        for atom in resonanceHybrid.atoms:
            if atom.symbol == 'C':
                number_carbon_bonds = sum([1 for bond in atom.bonds if bond.symbol=='C'])
                if number_carbon_bonds == 2:
                    atom.label = 'center'
                    symmetryNumber = calculateAtomSymmetryNumber(resonanceHybrid, atom)
                    self.assertEqual(symmetryNumber, 2)
                    pass

    def testEdgeCarbonAtomSymmetryNumberAllyl(self):
        """
        Test the Molecule.calculateAtomSymmetryNumber() for the edge carbons
        on [CH2]C=C
        """
        spc = Species(molecule=[Molecule().fromSMILES('[CH2]C=C')])
        molecule = spc.getResonanceHybrid()
        for atom in molecule.atoms:
            if atom.symbol =='C':
                number_carbon_bonds = sum([1 for bond in atom.bonds if bond.symbol=='C'])
                if number_carbon_bonds == 1:
                    atom.label = 'edge'
                    symmetryNumber = calculateAtomSymmetryNumber(molecule, atom)
                    self.assertEqual(symmetryNumber, 1)
        
    def testEdgeCarbonAtomSymmetryNumberAllylUsingResonanceHybrid(self):
        """
        Test the Molecule.calculateAtomSymmetryNumber() for the edge carbons
        on [CH2]C=C
        """
        molecule = Molecule().fromSMILES('[CH2]C=C')
        species = Species(molecule=[molecule])
        resonanceHybrid = species.getResonanceHybrid()
        for atom in resonanceHybrid.atoms:
            if atom.symbol =='C':
                number_carbon_bonds = sum([1 for bond in atom.bonds if bond.symbol=='C'])
                if number_carbon_bonds == 1:
                    atom.label = 'edge'
                    symmetryNumber = calculateAtomSymmetryNumber(resonanceHybrid, atom)
                    self.assertEqual(symmetryNumber, 1)

    def testAtomSymmetryNumberAllyl(self):
        """
        Test the Molecule.calculateAtomSymmetryNumber() on [CH2]C=C
        """
        spc = Species(molecule=[Molecule().fromSMILES('[CH2]C=C')])
        molecule = spc.getResonanceHybrid()
        symmetryNumber = 1
        for atom in molecule.atoms:
            if not molecule.isAtomInCycle(atom):
                symmetryNumber *= calculateAtomSymmetryNumber(molecule, atom)
        self.assertEqual(symmetryNumber, 2)

    def testBondSymmetryNumberAllyl(self):
        """
        Test the Molecule.calculateBondSymmetryNumber() on [CH2]C=C
        """
        spc = Species(molecule=[Molecule().fromSMILES('[CH2]C=C')])
        molecule = spc.getResonanceHybrid()
        symmetryNumber = 1
        for atom1 in molecule.atoms:
            for atom2 in atom1.bonds:
                if molecule.atoms.index(atom1) < molecule.atoms.index(atom2):
                    symmetryNumber *= calculateBondSymmetryNumber(molecule, atom1, atom2)
        self.assertEqual(symmetryNumber, 1)

    @work_in_progress
    def testTotalSymmetryNumberToluene(self):
        """
        Test the Species.getSymmetryNumber() (total symmetry) on c1ccccc1C
        """
        molecule = Molecule().fromSMILES('c1ccccc1C')
        species = Species(molecule=[molecule])
        symmetryNumber = species.getSymmetryNumber()
        self.assertEqual(symmetryNumber, 6)

    @work_in_progress
    def testTotalSymmetryNumberSpecialCyclic(self):
        """
        Test the Species.getSymmetryNumber() (total symmetry) from issue # 332
        """
        molecule = Molecule().fromSMILES('C1(C(C(C(C(C1C2CCC2)C3CCC3)C4CCC4)C5CCC5)C6CCC6)C7CCC7')
        species = Species(molecule=[molecule])
        symmetryNumber = species.getSymmetryNumber()
        self.assertEqual(symmetryNumber, 12)

    @work_in_progress
    def testTotalSymmetryNumberChlorobenzene(self):
        """
        Test the Species.getSymmetryNumber() (total symmetry) on c1ccccc1Cl
        """
        molecule = Molecule().fromSMILES('c1ccccc1Cl')
        species = Species(molecule=[molecule])
        symmetryNumber = species.getSymmetryNumber()
        self.assertEqual(symmetryNumber, 2)

    @work_in_progress
    def testTotalSymmetryNumberPhenoxyKecle(self):
        """
        Test the Species.getSymmetryNumber() (total symmetry) on c1ccccc1[O]
        using kecle structure
        """
        molecule = Molecule().fromSMILES('c1ccccc1[O]')
        species = Species(molecule=[molecule])
        symmetryNumber = species.getSymmetryNumber()
        self.assertEqual(symmetryNumber, 2)

    @work_in_progress
    def testTotalSymmetryNumberPhenoxyBenzene(self):
        """
        Test symmetry on c1ccccc1[O] using phenoxy benzene structure
        """
        molecule = Molecule().fromSMILES('c1ccccc1[O]')
        species = Species(molecule=[molecule])
        aromatic_molecule = generateAromaticResonanceStructures(molecule)[0]
        symmetryNumber = aromatic_molecule.getSymmetryNumber()
        self.assertEqual(symmetryNumber, 2)

    def testTotalSymmetryNumber12Dimethylbenzene(self):
        """
        Test the Species.getSymmetryNumber() (total symmetry) on Cc1ccccc1C
        """
        molecule = Molecule().fromSMILES('Cc1ccccc1C')
        species = Species(molecule=[molecule])
        symmetryNumber = species.getSymmetryNumber()
        self.assertEqual(symmetryNumber, 18)

    @work_in_progress
    def testTotalSymmetryNumber14Dimethylbenzene(self):
        """
        Test the Species.getSymmetryNumber() (total symmetry) on Cc1ccc(C)cc1
        """
        molecule = Molecule().fromSMILES('Cc1ccc(C)cc1')
        species = Species(molecule=[molecule])
        symmetryNumber = species.getSymmetryNumber()
        self.assertEqual(symmetryNumber, 36)
        
    def testTotalSymmetryNumberEthane(self):
        """
        Test the Molecule.calculateSymmetryNumber() on CC
        """
        self.assertEqual(Species().fromSMILES('CC').getSymmetryNumber(), 18)
    
    def testTotalSymmetryNumber1(self):
        """
        Test the Molecule.calculateSymmetryNumber() on C=C=[C]C(C)(C)[C]=C=C
        """
        self.assertEqual(Species().fromSMILES('C=C=[C]C(C)(C)[C]=C=C').getSymmetryNumber(), 18)
    
    @work_in_progress
    def testTotalSymmetryNumber2(self):
        """
        Test the Molecule.calculateSymmetryNumber() on C(=CC(c1ccccc1)C([CH]CCCCCC)C=Cc1ccccc1)[CH]CCCCCC
        """
        self.assertEqual(Species().fromSMILES('C(=CC(c1ccccc1)C([CH]CCCCCC)C=Cc1ccccc1)[CH]CCCCCC').getSymmetryNumber(), '36?')
    
    def testSymmetryNumberHydroxyl(self):
        """
        Test the Molecule.calculateSymmetryNumber() on [OH]
        """
        self.assertEqual(Species().fromSMILES('[OH]').getSymmetryNumber(), 1)
       
    def testSymmetryNumberOxygen(self):
        """
        Test the Molecule.calculateSymmetryNumber() on O=O
        """
        self.assertEqual(Species().fromSMILES('O=O').getSymmetryNumber(), 2)
        
    def testSymmetryNumberDicarbon(self):
        """
        Test the Molecule.calculateSymmetryNumber() on [C]#[C]
        """
        self.assertEqual(Species().fromSMILES('[C]#[C]').getSymmetryNumber(), 2)
    
    def testSymmetryNumberHydrogen(self):
        """
        Test the Molecule.calculateSymmetryNumber() on [H][H]
        """
        self.assertEqual(Species().fromSMILES('[H][H]').getSymmetryNumber(), 2)
    
    def testSymmetryNumberAcetylene(self):
        """
        Test the Molecule.calculateSymmetryNumber() on C#C
        """
        self.assertEqual(Species().fromSMILES('C#C').getSymmetryNumber(), 2)
    
    def testSymmetryNumberButadiyne(self):
        """
        Test the Molecule.calculateSymmetryNumber() on C#CC#C
        """
        self.assertEqual(Species().fromSMILES('C#CC#C').getSymmetryNumber(), 2)
    
    def testSymmetryNumberMethane(self):
        """
        Test the Molecule.calculateSymmetryNumber() on CH4
        """
        self.assertEqual(Species().fromSMILES('C').getSymmetryNumber(), 12)
    
    def testSymmetryNumberFormaldehyde(self):
        """
        Test the Molecule.calculateSymmetryNumber() on C=O
        """
        self.assertEqual(Species().fromSMILES('C=O').getSymmetryNumber(), 2)
    
    def testSymmetryNumberMethyl(self):
        """
        Test the Molecule.calculateSymmetryNumber() on [CH3]
        """
        self.assertEqual(Species().fromSMILES('[CH3]').getSymmetryNumber(), 6)
    
    def testSymmetryNumberWater(self):
        """
        Test the Molecule.calculateSymmetryNumber() on H2O
        """
        self.assertEqual(Species().fromSMILES('O').getSymmetryNumber(), 2)
    
    def testSymmetryNumberEthylene(self):
        """
        Test the Molecule.calculateSymmetryNumber() on C=C
        """
        self.assertEqual(Species().fromSMILES('C=C').getSymmetryNumber(), 4)
    
    def testSymmetryNumberEthenyl(self):
        """
        Test the Molecule.calculateSymmetryNumber() on C=[CH].
        """
        self.assertEqual(Species().fromSMILES('C=[CH]').getSymmetryNumber(), 1)
    
    def testSymmetryNumberCyclic(self):
        """
        Test the Molecule.calculateSymmetryNumber() on C1=C=C=1
        """
        self.assertEqual(Species().fromSMILES('C1=C=C=1').getSymmetryNumber(), 6)
    
################################################################################

if __name__ == '__main__':
    unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )
