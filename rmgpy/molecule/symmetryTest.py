#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2018 Prof. William H. Green (whgreen@mit.edu),           #
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

from rmgpy.molecule.molecule import Molecule
from rmgpy.molecule.symmetry import calculateAtomSymmetryNumber, calculateAxisSymmetryNumber, calculateBondSymmetryNumber, calculateCyclicSymmetryNumber
from rmgpy.species import Species
from rmgpy.molecule.resonance import generate_aromatic_resonance_structures
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
    
    def testAtomSymmetryNumberEthanewithDeuteriumTritium(self):
        """
        Test the Molecule.calculateAtomSymmetryNumber() on CC(D)(T)

        This is meant to test whether chirality is accounted for, which
        should half the symmetry term.

        The total number is 1.5 because the methyl group contributes *3 and
        the chiral center contributes *0.5
        """
        molecule = Molecule().fromAdjacencyList(
"""
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3 D u0 p0 c0 {1,S}
4 T u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}

""")
        symmetryNumber = 1
        for atom in molecule.atoms:
            if not molecule.isAtomInCycle(atom):
                symmetryNumber *= calculateAtomSymmetryNumber(molecule, atom)
        self.assertAlmostEqual(symmetryNumber, 1.5)

    def testAtomSymmetryNumberWithTwoChiralCenters(self):
        """
        Test the Molecule.calculateAtomSymmetryNumber() on [CH2]CC([CH2])C(C)C=C

        This is meant to test whether chirality is accounted for, which
        should half the symmetry term.

        The molecule has one methyl group (*3), two CH2dot groups (*2 each), and
        two chiral centers (*0.5), leading to a total atom symmetry number of 3
        """
        molecule = Molecule().fromAdjacencyList(
"""
multiplicity 3
1  C u1 p0 c0 {2,S} {3,S} {4,S}
2  H u0 p0 c0 {1,S}
3  H u0 p0 c0 {1,S}
4  C u0 p0 c0 {1,S} {5,S} {13,S} {14,S}
5  C u0 p0 c0 {4,S} {6,S} {9,S} {15,S}
6  C u1 p0 c0 {5,S} {7,S} {8,S}
7  H u0 p0 c0 {6,S}
8  H u0 p0 c0 {6,S}
9  C u0 p0 c0 {5,S} {10,S} {11,S} {16,S}
10 C u0 p0 c0 {9,S} {17,S} {18,S} {19,S}
11 C u0 p0 c0 {9,S} {12,D} {20,S}
12 C u0 p0 c0 {11,D} {21,S} {22,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {10,S}
20 H u0 p0 c0 {11,S}
21 H u0 p0 c0 {12,S}
22 H u0 p0 c0 {12,S}
""")
        symmetryNumber = 1
        for atom in molecule.atoms:
            if not molecule.isAtomInCycle(atom):
                symmetryNumber *= calculateAtomSymmetryNumber(molecule, atom)
        self.assertAlmostEqual(symmetryNumber, 3)

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

    def testCyclicSymmetryNumberCyclohexanone(self):
        """
        Test the Molecule.calculateCyclicSymmetryNumber() on C1CCCCC1=O
        """
        molecule = Molecule().fromSMILES('C1CCCCC1=O')
        symmetryNumber = calculateCyclicSymmetryNumber(molecule)
        self.assertEqual(symmetryNumber, 2)
        
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

    def testTotalSymmetryNumberToluene(self):
        """
        Test the Species.getSymmetryNumber() (total symmetry) on c1ccccc1C
        """
        molecule = Molecule().fromSMILES('c1ccccc1C')
        species = Species(molecule=[molecule])
        symmetryNumber = species.getSymmetryNumber()
        self.assertEqual(symmetryNumber, 6)

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
        aromatic_molecule = generate_aromatic_resonance_structures(molecule)[0]
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

    def test_symmetry_number_S1_total(self):
        """
        Test the Molecule.calculateSymmetryNumber() on [CH]1CCC1CC1CC1
        """
        self.assertEqual(Species().fromSMILES('[CH]1CCC1CC1CC1').getSymmetryNumber(),1)

    def testCyclicSymmetryNumberS1(self):
        """
        Test the Molecule.calculateCyclicSymmetryNumber() on [CH]1CCC1CC1CC1
        """
        molecule = Molecule().fromSMILES('[CH]1CCC1CC1CC1')
        symmetryNumber = calculateCyclicSymmetryNumber(molecule)
        self.assertEqual(symmetryNumber, 1)
    
    def testCyclicSymmetryNumberMethylCycloPropane(self):
        """
        Test the Molecule.calculateCyclicSymmetryNumber() on CC1CC1
        """
        molecule = Molecule().fromSMILES('CC1CC1')
        symmetryNumber = calculateCyclicSymmetryNumber(molecule)
        self.assertEqual(symmetryNumber, 1)

    def testCyclicSymmetryNumberMethylCycloPropene(self):
        """
        Test the Molecule.calculateCyclicSymmetryNumber() on C=C1CC1
        """
        molecule = Molecule().fromSMILES('C=C1CC1')
        symmetryNumber = calculateCyclicSymmetryNumber(molecule)
        self.assertEqual(symmetryNumber, 2)


    def testCyclicSymmetryNumberMethylCycloButene(self):
        """
        Test the Molecule.calculateCyclicSymmetryNumber() on C=C1CCC1
        """
        molecule = Molecule().fromSMILES('C=C1CCC1')
        symmetryNumber = calculateCyclicSymmetryNumber(molecule)
        self.assertEqual(symmetryNumber, 2)

    def testCyclicSymmetryNumberMethylCycloButane(self):
        """
        Test the Molecule.calculateCyclicSymmetryNumber() on CC1CCC1
        """
        molecule = Molecule().fromSMILES('CC1CCC1')
        symmetryNumber = calculateCyclicSymmetryNumber(molecule)
        self.assertEqual(symmetryNumber, 1)

    def testCyclicSymmetryNumberDiMethylCycloButane(self):
        """
        Test the Molecule.calculateCyclicSymmetryNumber() on CC1CC(C)C1
        """
        molecule = Molecule().fromSMILES('CC1CC(C)C1')
        symmetryNumber = calculateCyclicSymmetryNumber(molecule)
        self.assertEqual(symmetryNumber, 4)


    def test_symmetry_number_dimethylcylcobutane_total(self):
        """
        Test the Molecule.calculateSymmetryNumber() on CC1CC(C)C1
        """
        self.assertEqual(Species().fromSMILES('CC1CC(C)C1').getSymmetryNumber(),36)
################################################################################

if __name__ == '__main__':
    unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )
