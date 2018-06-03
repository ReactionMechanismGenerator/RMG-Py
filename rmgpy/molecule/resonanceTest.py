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

from .molecule import Molecule

from .resonance import *
from .resonance import _clar_optimization, _clar_transformation

class ResonanceTest(unittest.TestCase):

    def testAllylShift(self):
        """Test allyl shift for hexadienyl radical"""
        molList = generate_resonance_structures(Molecule(SMILES="C=C[CH]C=CC"))
        self.assertEqual(len(molList), 3)

    def testTriradAllylShift(self):
        """Test allyl shift for a tri-rad carbon"""
        molList = generate_resonance_structures(Molecule(SMILES="[C]N=N"))
        self.assertEqual(len(molList), 3)
        self.assertTrue(any([any([atom.charge != 0 for atom in mol.vertices]) for mol in molList]))  # expecting [C]=[N+.][NH-]

    def testOxime(self):
        """Test resonance structure generation for CC=N[O] radical

        Simple case for lone pair <=> radical resonance"""
        molList = generate_resonance_structures(Molecule(SMILES="CC=N[O]"))
        self.assertEqual(len(molList), 3)
        self.assertTrue(any([any([atom.charge != 0 for atom in mol.vertices]) for mol in molList]))

    def testRingAllylShift(self):
        """Test allyl shift for a cyclic species with heteroatoms"""
        molList = generate_resonance_structures(Molecule(SMILES="[CH]1C=NC=N1"))
        self.assertEqual(len(molList), 3)

    def testCarbeneAllylShift(self):
        """Test allyl shift for a carbene species"""
        molList = generate_resonance_structures(Molecule(SMILES="[C]=C=O"))
        self.assertEqual(len(molList), 2)

    def testCH2CHCHO(self):
        """Test resonance structure generation for C=C[CH][O] bi-radical

        Test case for allyl bi-radical resonance"""
        molList = generate_resonance_structures(Molecule(SMILES="C=C[CH][O]"))
        self.assertEqual(len(molList), 3)

    def CH2NO(self):
        """Test combined resonance transitions of allyl-shift and lonePair-radical"""
        molList = generate_resonance_structures(Molecule(SMILES="[CH2]N=O"))
        self.assertEqual(len(molList), 3)

    def CH3S2O2(self):
        """Test combined resonance transitions of one_pair_radical_multiple_bond"""
        molList = generate_resonance_structures(Molecule(SMILES="CSS(=O)[O]"))
        self.assertEqual(len(molList), 3)

    def N2SO2(self):
        """Test the resonance transitions of a species with several hereroatoms and several multiple bonds"""
        molList = generate_resonance_structures(Molecule(SMILES="[N-]=[N+]=S(=O)=O"))
        self.assertEqual(len(molList), 2)

    def NSH(self):
        """Test that a resonance structure with a minimal octet deviation but higher charge span is filtered out"""
        molList = generate_resonance_structures(Molecule(SMILES="N#S"))
        self.assertEqual(len(molList), 1)
        self.assertTrue(all([atom.charge == 0 for atom in molList[0].vertices]))

    def testNCO(self):
        """Test resonance structure generation for NCO

        NCO should only have two resonance structures [N.]=C=O <=> N#C[O.], and not a third structure which has
        the same octet deviation, has a charge separation, but no ne radical site: [N+.]#C[O-]"""
        mol_list = generate_resonance_structures(Molecule(SMILES="[N]=C=O"))
        self.assertEqual(len(mol_list), 2)
        self.assertTrue(all([all([atom.charge == 0 for atom in mol.vertices]) for mol in mol_list]))  # none of the
        # structures should be charged

    def testNO2(self):
        """Test resonance structure generation for [O]N=O radical

        Test case for the lone pair <=> radical resonance transition.
        Also tests that the filtering function allows charge separation when the radical site is changed."""
        molList = generate_resonance_structures(Molecule(SMILES="[O]N=O"))
        self.assertEqual(len(molList), 2)
        self.assertTrue(any([any([atom.charge != 0 for atom in mol.vertices]) for mol in molList]))  # one of the
        # structures should be charged

    def testN2O(self):
        """Test resonance structure generation for N#[N+][O-]

        A classic N5ddc <=> N5tc resonance transition"""
        molList = generate_resonance_structures(Molecule(SMILES="N#[N+][O-]"))
        self.assertEqual(len(molList), 2)
        self.assertTrue(all([any([atom.charge != 0 for atom in mol.vertices]) for mol in molList]))  # both structures
        # should have some charged atoms

        sbonds = 0
        dbonds = 0
        tbonds = 0
        for mol in molList:
            for atom in mol.atoms:
                for bond in atom.bonds.itervalues():
                    if bond.isSingle():
                        sbonds += 1
                    elif bond.isDouble():
                        dbonds += 1
                    elif bond.isTriple():
                        tbonds += 1
        self.assertEqual(sbonds / 2, 1)  # each bond is counted twice above
        self.assertEqual(dbonds / 2, 2)
        self.assertEqual(tbonds / 2, 1)

    def testAzide(self):
        """Test resonance structure generation for ethyl azide

        Simple case for N5ddc <=> N5tc resonance
        Azides are described by three resonance structures: N=[N+]=[N-] <=> [NH-][N+]#N <=> [NH+]#[N+][N-2]
        However, since the third does not contribute to reactivity and has a higher charge span, it is filtered out"""
        molList = generate_resonance_structures(Molecule(SMILES="CCN=[N+]=[N-]"))
        self.assertEqual(len(molList), 2)
        self.assertTrue(all([any([atom.charge != 0 for atom in mol.vertices]) for mol in molList]))

    def test_ozone(self):
        """Test resonance structure generation for O3, S3 and SO2.

        Compare that these iso-electronic structures have the same number of resonance structures"""
        mol_list_1 = generate_resonance_structures(Molecule(SMILES="[O-][O+]=O"))
        self.assertEqual(len(mol_list_1), 1)
        mol_list_2 = generate_resonance_structures(Molecule(SMILES="O=S=O"))
        self.assertEqual(len(mol_list_2), 1)
        mol_list_3 = generate_resonance_structures(Molecule(SMILES="S=S=S"))
        self.assertEqual(len(mol_list_3), 1)

    def testHCOvsHCS(self):
        """Test resonance structure generation for [CH]=O and [CH]=S

        These iso-electronic structures have a different(!) number of resonance structures"""
        mol_list_1 = generate_resonance_structures(Molecule(SMILES="[CH]=O"))
        self.assertEqual(len(mol_list_1), 1)
        mol_list_2 = generate_resonance_structures(Molecule(SMILES="[CH]=S"))
        self.assertEqual(len(mol_list_2), 2)

    def testNO(self):
        """Test that an incorrect NO structure [::N][::O.] is correctly identified as [:N.]=[::O]

        The incorrect structure could be generated from HON (O[::N]) during an RMG run, and should be identified as NO.
        The original structure should be kept as unreactive (appended at the end of the molecule list)"""
        mol_list = generate_resonance_structures(Molecule().fromAdjacencyList("""multiplicity 2
                                                                                 1 N u0 p2 c0 {2,S}
                                                                                 2 O u1 p2 c0 {1,S}"""))
        self.assertEqual(len(mol_list), 2)
        self.assertTrue(mol_list[0].reactive)
        self.assertFalse(mol_list[1].reactive)
        self.assertEqual(mol_list[0].vertices[0].lonePairs + mol_list[0].vertices[1].lonePairs, 3)
        self.assertEqual(mol_list[1].vertices[0].lonePairs + mol_list[1].vertices[1].lonePairs, 4)

    def testN5dc(self):
        """Test the N5dc resonance transformation

        We should see N=[N+]([O])([O-]) <=> [NH-][N+]([O])=O
        Two isomorphic structures should be included in molList: N=[N+]([O])([O-]) <=> N=[N+]([O-])([O])"""
        mol = Molecule(SMILES="N=[N+]([O-])[O]")
        mol_list = generate_resonance_structures(mol, keep_isomorphic=True)
        self.assertEqual(len(mol_list), 6)
        isomorphic_counter = 0
        negatively_charged_nitrogen = 0
        for mol1 in mol_list:
            if mol1.isIsomorphic(mol):
                isomorphic_counter += 1
            for atom in mol1.vertices:
                if atom.isNitrogen() and atom.charge < 0:
                    negatively_charged_nitrogen += 1
        self.assertEquals(isomorphic_counter, 2)
        self.assertEquals(negatively_charged_nitrogen, 2)

    def testStyryl1(self):
        """Test resonance structure generation for styryl, with radical on branch

        In this case, the radical can be delocalized into the aromatic ring"""
        molList = generate_resonance_structures(Molecule(SMILES="c1ccccc1[C]=C"))
        self.assertEqual(len(molList), 4)

    def testStyryl2(self):
        """Test resonance structure generation for styryl, with radical on ring

        In this case, the radical can be delocalized into the aromatic ring"""
        molList = generate_resonance_structures(Molecule(SMILES="C=C=C1C=C[CH]C=C1"))
        self.assertEqual(len(molList), 4)

    def testNaphthyl(self):
        """Test resonance structure generation for naphthyl radical

        In this case, the radical is orthogonal to the pi-orbital plane and cannot delocalize"""
        molList = generate_resonance_structures(Molecule(SMILES="c12[c]cccc1cccc2"))
        self.assertEqual(len(molList), 4)

    def testMethylNapthalene(self):
        """Test resonance structure generation for methyl naphthalene

        Example of stable polycyclic aromatic species"""
        molList = generate_resonance_structures(Molecule(SMILES="CC1=CC=CC2=CC=CC=C12"))
        self.assertEqual(len(molList), 4)

    def testMethylPhenanthrene(self):
        """Test resonance structure generation for methyl phenanthrene

        Example of stable polycyclic aromatic species"""
        molList = generate_resonance_structures(Molecule(SMILES="CC1=CC=CC2C3=CC=CC=C3C=CC=21"))
        self.assertEqual(len(molList), 3)

    def testMethylPhenanthreneRadical(self):
        """Test resonance structure generation for methyl phenanthrene radical

        Example radical polycyclic aromatic species where the radical can delocalize"""
        molList = generate_resonance_structures(Molecule(SMILES="[CH2]C1=CC=CC2C3=CC=CC=C3C=CC=21"))
        self.assertEqual(len(molList), 9)

    def testAromaticWithLonePairResonance(self):
        """Test resonance structure generation for aromatic species with lone pair <=> radical resonance"""
        molList = generate_resonance_structures(Molecule(SMILES="c1ccccc1CC=N[O]"))
        self.assertEqual(len(molList), 6)

    def testAromaticWithNResonance(self):
        """Test resonance structure generation for aromatic species with N5ddc <=> N5tc resonance"""
        molList = generate_resonance_structures(Molecule(SMILES="c1ccccc1CCN=[N+]=[N-]"))
        self.assertEqual(len(molList), 4)
        # TODO: this test cannot be run because RDKit (which checks for aromaticity) cannot process hyper-valence N

    def testNoClarStructures(self):
        """Test that we can turn off Clar structure generation."""
        molList = generate_resonance_structures(Molecule(SMILES='C1=CC=CC2C3=CC=CC=C3C=CC=21'), clar_structures=False)
        self.assertEqual(len(molList), 2)

    def testC13H11Rad(self):
        """Test resonance structure generation for p-methylbenzylbenzene radical

        Has multiple resonance structures that break aromaticity of a ring"""
        molList = generate_resonance_structures(Molecule(SMILES="[CH](c1ccccc1)c1ccc(C)cc1"))
        self.assertEqual(len(molList), 6)

    def testC8H8(self):
        """Test resonance structure generation for 5,6-dimethylene-1,3-cyclohexadiene

        Example of molecule that RDKit considers aromatic, but RMG does not"""
        molList = generate_resonance_structures(Molecule(SMILES="C=C1C=CC=CC1=C"))
        self.assertEqual(len(molList), 1)

    def testC8H7J(self):
        """Test resonance structure generation for 5,6-dimethylene-1,3-cyclohexadiene radical

        Example of molecule that RDKit considers aromatic, but RMG does not"""
        molList = generate_resonance_structures(Molecule(SMILES="C=C1C=CC=CC1=[CH]"))
        self.assertEqual(len(molList), 1)

    def testC8H7J2(self):
        """Test resonance structure generation for 5,6-dimethylene-1,3-cyclohexadiene radical

        Example of molecule that RDKit considers aromatic, but RMG does not"""
        molList = generate_resonance_structures(Molecule(SMILES="C=C1C=[C]C=CC1=C"))
        self.assertEqual(len(molList), 1)

    def test_C9H9_aro(self):
        """Test cyclopropyl benzene radical, aromatic SMILES"""
        mol = Molecule(SMILES="[CH]1CC1c1ccccc1")
        molList = generate_resonance_structures(mol)
        self.assertEqual(len(molList), 2)
    
    def test_C9H9_kek(self):
        """Test cyclopropyl benzene radical, kekulized SMILES"""
        mol = Molecule(SMILES="[CH]1CC1C1C=CC=CC=1")
        molList = generate_resonance_structures(mol)
        self.assertEqual(len(molList), 2)

    def test_Benzene_aro(self):
        """Test benzene, aromatic SMILES"""
        mol = Molecule(SMILES="c1ccccc1")
        molList = generate_resonance_structures(mol)
        self.assertEqual(len(molList), 2)
    
    def test_Benzene_kek(self):
        """Test benzene, kekulized SMILES"""
        mol = Molecule(SMILES="C1C=CC=CC=1")
        molList = generate_resonance_structures(mol)
        self.assertEqual(len(molList), 2)

    def test_C9H11_aro(self):
        """Test propylbenzene radical, aromatic SMILES"""
        mol = Molecule(SMILES="[CH2]CCc1ccccc1")
        molList = generate_resonance_structures(mol)
        self.assertEqual(len(molList), 2)

    def test_C10H11_aro(self):
        """Test cyclobutylbenzene radical, aromatic SMILES"""
        mol = Molecule(SMILES="[CH]1CCC1c1ccccc1")
        molList = generate_resonance_structures(mol)
        self.assertEqual(len(molList), 2)

    def test_C9H10_aro(self):
        """Test cyclopropylbenzene, aromatic SMILES"""
        mol = Molecule(SMILES="C1CC1c1ccccc1")
        molList = generate_resonance_structures(mol)
        self.assertEqual(len(molList), 2)

    def test_C10H12_aro(self):
        """Test cyclopropylmethyl benzene, aromatic SMILES"""
        mol = Molecule(SMILES="C1CC1c1c(C)cccc1")
        molList = generate_resonance_structures(mol)
        self.assertEqual(len(molList), 3)

    def test_C9H10_aro_2(self):
        """Test cyclopropyl benzene, generate aromatic resonance isomers"""
        mol = Molecule(SMILES="C1CC1c1ccccc1")
        molList = generate_aromatic_resonance_structures(mol)
        self.assertEqual(len(molList), 1)

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
        for isomer in generate_aromatic_resonance_structures(perylene2):
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
        for isomer in generate_aromatic_resonance_structures(naphthalene2):
            if naphthalene.isIsomorphic(isomer):
                break
        else:  # didn't break
            self.fail("{} isn't isomorphic with any aromatic forms of {}".format(
                naphthalene.toSMILES(),
                naphthalene2.toSMILES()
            ))

    def testAromaticResonanceStructures(self):
        """Test that generate_aromatic_resonance_structures gives consistent output

        Check that we get the same resonance structure regardless of which structure we start with"""
        # Kekulized form, radical on methyl
        struct1 = Molecule().fromAdjacencyList("""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,D} {7,S}
2  C u0 p0 c0 {1,S} {4,S} {8,D}
3  C u0 p0 c0 {1,D} {5,S} {11,S}
4  C u0 p0 c0 {2,S} {9,D} {10,S}
5  C u0 p0 c0 {3,S} {6,D} {15,S}
6  C u0 p0 c0 {5,D} {12,S} {16,S}
7  C u0 p0 c0 {1,S} {12,D} {18,S}
8  C u0 p0 c0 {2,D} {13,S} {19,S}
9  C u0 p0 c0 {4,D} {14,S} {22,S}
10 C u0 p0 c0 {4,S} {11,D} {23,S}
11 C u0 p0 c0 {3,S} {10,D} {24,S}
12 C u0 p0 c0 {6,S} {7,D} {17,S}
13 C u0 p0 c0 {8,S} {14,D} {20,S}
14 C u0 p0 c0 {9,S} {13,D} {21,S}
15 C u1 p0 c0 {5,S} {25,S} {26,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {12,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {13,S}
21 H u0 p0 c0 {14,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {11,S}
25 H u0 p0 c0 {15,S}
26 H u0 p0 c0 {15,S}
""")
        # Kekulized form, radical on ring
        struct2 = Molecule().fromAdjacencyList("""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {7,D}
2  C u0 p0 c0 {1,S} {4,S} {8,D}
3  C u0 p0 c0 {1,S} {5,S} {11,D}
4  C u0 p0 c0 {2,S} {9,S} {10,D}
5  C u0 p0 c0 {3,S} {6,S} {15,D}
6  C u0 p0 c0 {5,S} {12,D} {16,S}
7  C u0 p0 c0 {1,D} {12,S} {17,S}
8  C u0 p0 c0 {2,D} {13,S} {18,S}
9  C u0 p0 c0 {4,S} {14,D} {19,S}
10 C u0 p0 c0 {4,D} {11,S} {20,S}
11 C u0 p0 c0 {3,D} {10,S} {21,S}
12 C u0 p0 c0 {6,D} {7,S} {22,S}
13 C u1 p0 c0 {8,S} {14,S} {23,S}
14 C u0 p0 c0 {9,D} {13,S} {24,S}
15 C u0 p0 c0 {5,D} {25,S} {26,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {10,S}
21 H u0 p0 c0 {11,S}
22 H u0 p0 c0 {12,S}
23 H u0 p0 c0 {13,S}
24 H u0 p0 c0 {14,S}
25 H u0 p0 c0 {15,S}
26 H u0 p0 c0 {15,S}
""")
        # Aromatic form
        struct3 = Molecule().fromAdjacencyList("""
multiplicity 2
1  C u0 p0 c0 {2,B} {3,B} {7,B}
2  C u0 p0 c0 {1,B} {4,B} {8,B}
3  C u0 p0 c0 {1,B} {5,B} {11,B}
4  C u0 p0 c0 {2,B} {9,B} {10,B}
5  C u0 p0 c0 {3,B} {6,B} {15,S}
6  C u0 p0 c0 {5,B} {12,B} {16,S}
7  C u0 p0 c0 {1,B} {12,B} {18,S}
8  C u0 p0 c0 {2,B} {13,B} {19,S}
9  C u0 p0 c0 {4,B} {14,B} {22,S}
10 C u0 p0 c0 {4,B} {11,B} {23,S}
11 C u0 p0 c0 {3,B} {10,B} {24,S}
12 C u0 p0 c0 {6,B} {7,B} {17,S}
13 C u0 p0 c0 {8,B} {14,B} {20,S}
14 C u0 p0 c0 {9,B} {13,B} {21,S}
15 C u1 p0 c0 {5,S} {25,S} {26,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {12,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {13,S}
21 H u0 p0 c0 {14,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {11,S}
25 H u0 p0 c0 {15,S}
26 H u0 p0 c0 {15,S}
""")
        result1 = generate_aromatic_resonance_structures(struct1)
        result2 = generate_aromatic_resonance_structures(struct2)
        result3 = generate_aromatic_resonance_structures(struct3)

        self.assertEqual(len(result1), 1)
        self.assertEqual(len(result2), 1)
        self.assertEqual(len(result3), 1)

        self.assertTrue(result1[0].isIsomorphic(result2[0]))
        self.assertTrue(result1[0].isIsomorphic(result3[0]))

    @work_in_progress
    def testBridgedAromatic(self):
        """Test that we can handle bridged aromatics.

        This is affected by how we perceive rings. Using getSmallestSetOfSmallestRings gives
        non-deterministic output, so using getAllCyclesOfSize allows this test to pass."""
        mol = Molecule(SMILES='c12c3cccc1c3ccc2')
        arom = Molecule().fromAdjacencyList("""
1  C u0 p0 c0 {2,B} {3,B} {8,B}
2  C u0 p0 c0 {1,B} {4,B} {5,B}
3  C u0 p0 c0 {1,B} {4,S} {6,B}
4  C u0 p0 c0 {2,B} {3,S} {7,B}
5  C u0 p0 c0 {2,B} {9,B} {11,S}
6  C u0 p0 c0 {3,B} {9,B} {13,S}
7  C u0 p0 c0 {4,B} {10,B} {14,S}
8  C u0 p0 c0 {1,B} {10,B} {16,S}
9  C u0 p0 c0 {5,B} {6,B} {12,S}
10 C u0 p0 c0 {7,B} {8,B} {15,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {8,S}
""")

        out = generate_resonance_structures(mol)

        self.assertEqual(len(out), 3)
        self.assertTrue(arom.isIsomorphic(out[1]))

    def testPolycyclicAromaticWithNonAromaticRing(self):
        """Test that we can make aromatic resonance structures when there is a pseudo-aromatic ring.

        This applies in cases where RDKit misidentifies one ring as aromatic, but there are other
        rings in the molecule that are actually aromatic."""
        mol = Molecule(SMILES='c1c2cccc1C(=C)C=[C]2')
        arom = Molecule().fromAdjacencyList("""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,B} {5,B}
2  C u0 p0 c0 {1,S} {8,S} {9,D}
3  C u0 p0 c0 {4,B} {6,B} {10,S}
4  C u0 p0 c0 {1,B} {3,B} {11,S}
5  C u0 p0 c0 {1,B} {7,B} {14,S}
6  C u0 p0 c0 {3,B} {7,B} {12,S}
7  C u0 p0 c0 {5,B} {6,B} {13,S}
8  C u0 p0 c0 {2,S} {10,D} {15,S}
9  C u0 p0 c0 {2,D} {16,S} {17,S}
10 C u1 p0 c0 {3,S} {8,D}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {9,S}
""")

        out = generate_resonance_structures(mol)

        self.assertEqual(len(out), 2)
        self.assertTrue(arom.isIsomorphic(out[1]))

    def testPolycyclicAromaticWithNonAromaticRing2(self):
        """Test that we can make aromatic resonance structures when there is a pseudo-aromatic ring.

        This applies in cases where RDKit misidentifies one ring as aromatic, but there are other
        rings in the molecule that are actually aromatic."""
        mol = Molecule(SMILES='C=C(C1=CC2=C(C=C1C=C3)C4=CC5=CC=CC=C5C=C4C=C2)C3=C')
        arom = Molecule().fromAdjacencyList("""
1  C u0 p0 c0 {4,S} {6,B} {11,B}
2  C u0 p0 c0 {3,B} {5,B} {12,B}
3  C u0 p0 c0 {2,B} {9,B} {13,B}
4  C u0 p0 c0 {1,S} {10,S} {23,D}
5  C u0 p0 c0 {2,B} {11,B} {20,B}
6  C u0 p0 c0 {1,B} {12,B} {15,S}
7  C u0 p0 c0 {8,B} {13,B} {17,B}
8  C u0 p0 c0 {7,B} {14,B} {18,B}
9  C u0 p0 c0 {3,B} {14,B} {19,B}
10 C u0 p0 c0 {4,S} {16,S} {24,D}
11 C u0 p0 c0 {1,B} {5,B} {25,S}
12 C u0 p0 c0 {2,B} {6,B} {26,S}
13 C u0 p0 c0 {3,B} {7,B} {29,S}
14 C u0 p0 c0 {8,B} {9,B} {34,S}
15 C u0 p0 c0 {6,S} {16,D} {27,S}
16 C u0 p0 c0 {10,S} {15,D} {28,S}
17 C u0 p0 c0 {7,B} {21,B} {30,S}
18 C u0 p0 c0 {8,B} {22,B} {33,S}
19 C u0 p0 c0 {9,B} {20,B} {35,S}
20 C u0 p0 c0 {5,B} {19,B} {36,S}
21 C u0 p0 c0 {17,B} {22,B} {31,S}
22 C u0 p0 c0 {18,B} {21,B} {32,S}
23 C u0 p0 c0 {4,D} {37,S} {38,S}
24 C u0 p0 c0 {10,D} {39,S} {40,S}
25 H u0 p0 c0 {11,S}
26 H u0 p0 c0 {12,S}
27 H u0 p0 c0 {15,S}
28 H u0 p0 c0 {16,S}
29 H u0 p0 c0 {13,S}
30 H u0 p0 c0 {17,S}
31 H u0 p0 c0 {21,S}
32 H u0 p0 c0 {22,S}
33 H u0 p0 c0 {18,S}
34 H u0 p0 c0 {14,S}
35 H u0 p0 c0 {19,S}
36 H u0 p0 c0 {20,S}
37 H u0 p0 c0 {23,S}
38 H u0 p0 c0 {23,S}
39 H u0 p0 c0 {24,S}
40 H u0 p0 c0 {24,S}
""")

        out = generate_resonance_structures(mol)

        self.assertEqual(len(out), 4)
        self.assertTrue(arom.isIsomorphic(out[1]))

    def testKekulizeBenzene(self):
        """Test that we can kekulize benzene."""
        arom = Molecule().fromAdjacencyList("""
1  C u0 p0 c0 {2,B} {6,B} {7,S}
2  C u0 p0 c0 {1,B} {3,B} {8,S}
3  C u0 p0 c0 {2,B} {4,B} {9,S}
4  C u0 p0 c0 {3,B} {5,B} {10,S}
5  C u0 p0 c0 {4,B} {6,B} {11,S}
6  C u0 p0 c0 {1,B} {5,B} {12,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
""")
        keku = Molecule().fromAdjacencyList("""
1  C u0 p0 c0 {2,D} {6,S} {7,S}
2  C u0 p0 c0 {1,D} {3,S} {8,S}
3  C u0 p0 c0 {2,S} {4,D} {9,S}
4  C u0 p0 c0 {3,D} {5,S} {10,S}
5  C u0 p0 c0 {4,S} {6,D} {11,S}
6  C u0 p0 c0 {1,S} {5,D} {12,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
""")
        out = generate_kekule_structure(arom)

        self.assertEqual(len(out), 1)
        self.assertTrue(out[0].isIsomorphic(keku))

    def testKekulizeNaphthalene(self):
        """Test that we can kekulize naphthalene."""
        arom = Molecule().fromAdjacencyList("""
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
        out = generate_kekule_structure(arom)

        self.assertEqual(len(out), 1)
        self.assertFalse(out[0].isAromatic())

        bonds = set()
        for atom in out[0].atoms:
            for bond in atom.bonds.itervalues():
                bonds.add(bond)
        dBonds = 0
        for bond in bonds:
            if bond.isDouble():
                dBonds += 1

        self.assertEqual(dBonds, 5)

    def testKekulizePhenanthrene(self):
        """Test that we can kekulize phenanthrene."""
        arom = Molecule().fromAdjacencyList("""
1  C u0 p0 c0 {2,B} {3,B} {5,B}
2  C u0 p0 c0 {1,B} {4,B} {9,B}
3  C u0 p0 c0 {1,B} {6,B} {10,B}
4  C u0 p0 c0 {2,B} {7,B} {8,B}
5  C u0 p0 c0 {1,B} {12,B} {17,S}
6  C u0 p0 c0 {3,B} {7,B} {18,S}
7  C u0 p0 c0 {4,B} {6,B} {19,S}
8  C u0 p0 c0 {4,B} {13,B} {20,S}
9  C u0 p0 c0 {2,B} {14,B} {23,S}
10 C u0 p0 c0 {3,B} {11,B} {24,S}
11 C u0 p0 c0 {10,B} {12,B} {15,S}
12 C u0 p0 c0 {5,B} {11,B} {16,S}
13 C u0 p0 c0 {8,B} {14,B} {21,S}
14 C u0 p0 c0 {9,B} {13,B} {22,S}
15 H u0 p0 c0 {11,S}
16 H u0 p0 c0 {12,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {13,S}
22 H u0 p0 c0 {14,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {10,S}
""")
        out = generate_kekule_structure(arom)

        self.assertEqual(len(out), 1)
        self.assertFalse(out[0].isAromatic())

        bonds = set()
        for atom in out[0].atoms:
            for bond in atom.bonds.itervalues():
                bonds.add(bond)
        dBonds = 0
        for bond in bonds:
            if bond.isDouble():
                dBonds += 1

        self.assertEqual(dBonds, 7)

    def testKekulizePyrene(self):
        """Test that we can kekulize pyrene."""
        arom = Molecule().fromAdjacencyList("""
1  C u0 p0 c0 {2,B} {3,B} {6,B}
2  C u0 p0 c0 {1,B} {4,B} {5,B}
3  C u0 p0 c0 {1,B} {7,B} {8,B}
4  C u0 p0 c0 {2,B} {9,B} {10,B}
5  C u0 p0 c0 {2,B} {11,B} {12,B}
6  C u0 p0 c0 {1,B} {13,B} {14,B}
7  C u0 p0 c0 {3,B} {15,B} {18,S}
8  C u0 p0 c0 {3,B} {9,B} {19,S}
9  C u0 p0 c0 {4,B} {8,B} {20,S}
10 C u0 p0 c0 {4,B} {16,B} {21,S}
11 C u0 p0 c0 {5,B} {16,B} {23,S}
12 C u0 p0 c0 {5,B} {13,B} {24,S}
13 C u0 p0 c0 {6,B} {12,B} {25,S}
14 C u0 p0 c0 {6,B} {15,B} {26,S}
15 C u0 p0 c0 {7,B} {14,B} {17,S}
16 C u0 p0 c0 {10,B} {11,B} {22,S}
17 H u0 p0 c0 {15,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {16,S}
23 H u0 p0 c0 {11,S}
24 H u0 p0 c0 {12,S}
25 H u0 p0 c0 {13,S}
26 H u0 p0 c0 {14,S}
""")
        out = generate_kekule_structure(arom)

        self.assertEqual(len(out), 1)
        self.assertFalse(out[0].isAromatic())

        bonds = set()
        for atom in out[0].atoms:
            for bond in atom.bonds.itervalues():
                bonds.add(bond)
        dBonds = 0
        for bond in bonds:
            if bond.isDouble():
                dBonds += 1

        self.assertEqual(dBonds, 8)

    def testKekulizeCorannulene(self):
        """Test that we can kekulize corannulene."""
        arom = Molecule().fromAdjacencyList("""
1  C u0 p0 c0 {2,B} {5,B} {8,B}
2  C u0 p0 c0 {1,B} {3,B} {10,B}
3  C u0 p0 c0 {2,B} {4,B} {9,B}
4  C u0 p0 c0 {3,B} {5,B} {6,B}
5  C u0 p0 c0 {1,B} {4,B} {7,B}
6  C u0 p0 c0 {4,B} {12,B} {13,B}
7  C u0 p0 c0 {5,B} {14,B} {15,B}
8  C u0 p0 c0 {1,B} {16,B} {20,B}
9  C u0 p0 c0 {3,B} {11,B} {17,B}
10 C u0 p0 c0 {2,B} {18,B} {19,B}
11 C u0 p0 c0 {9,B} {12,B} {21,S}
12 C u0 p0 c0 {6,B} {11,B} {22,S}
13 C u0 p0 c0 {6,B} {14,B} {23,S}
14 C u0 p0 c0 {7,B} {13,B} {24,S}
15 C u0 p0 c0 {7,B} {16,B} {25,S}
16 C u0 p0 c0 {8,B} {15,B} {26,S}
17 C u0 p0 c0 {9,B} {18,B} {27,S}
18 C u0 p0 c0 {10,B} {17,B} {28,S}
19 C u0 p0 c0 {10,B} {20,B} {29,S}
20 C u0 p0 c0 {8,B} {19,B} {30,S}
21 H u0 p0 c0 {11,S}
22 H u0 p0 c0 {12,S}
23 H u0 p0 c0 {13,S}
24 H u0 p0 c0 {14,S}
25 H u0 p0 c0 {15,S}
26 H u0 p0 c0 {16,S}
27 H u0 p0 c0 {17,S}
28 H u0 p0 c0 {18,S}
29 H u0 p0 c0 {19,S}
30 H u0 p0 c0 {20,S}
""")
        out = generate_kekule_structure(arom)

        self.assertEqual(len(out), 1)
        self.assertFalse(out[0].isAromatic())

        bonds = set()
        for atom in out[0].atoms:
            for bond in atom.bonds.itervalues():
                bonds.add(bond)
        dBonds = 0
        for bond in bonds:
            if bond.isDouble():
                dBonds += 1

        self.assertEqual(dBonds, 10)

    def testKekulizeCoronene(self):
        """Test that we can kekulize coronene."""
        arom = Molecule().fromAdjacencyList("""
1  C u0 p0 c0 {2,B} {6,B} {12,B}
2  C u0 p0 c0 {1,B} {3,B} {7,B}
3  C u0 p0 c0 {2,B} {4,B} {8,B}
4  C u0 p0 c0 {3,B} {5,B} {9,B}
5  C u0 p0 c0 {4,B} {6,B} {10,B}
6  C u0 p0 c0 {1,B} {5,B} {11,B}
7  C u0 p0 c0 {2,B} {14,B} {15,B}
8  C u0 p0 c0 {3,B} {16,B} {17,B}
9  C u0 p0 c0 {4,B} {18,B} {19,B}
10 C u0 p0 c0 {5,B} {20,B} {21,B}
11 C u0 p0 c0 {6,B} {22,B} {23,B}
12 C u0 p0 c0 {1,B} {13,B} {24,B}
13 C u0 p0 c0 {12,B} {14,B} {25,S}
14 C u0 p0 c0 {7,B} {13,B} {26,S}
15 C u0 p0 c0 {7,B} {16,B} {27,S}
16 C u0 p0 c0 {8,B} {15,B} {28,S}
17 C u0 p0 c0 {8,B} {18,B} {29,S}
18 C u0 p0 c0 {9,B} {17,B} {30,S}
19 C u0 p0 c0 {9,B} {20,B} {31,S}
20 C u0 p0 c0 {10,B} {19,B} {32,S}
21 C u0 p0 c0 {10,B} {22,B} {33,S}
22 C u0 p0 c0 {11,B} {21,B} {34,S}
23 C u0 p0 c0 {11,B} {24,B} {35,S}
24 C u0 p0 c0 {12,B} {23,B} {36,S}
25 H u0 p0 c0 {13,S}
26 H u0 p0 c0 {14,S}
27 H u0 p0 c0 {15,S}
28 H u0 p0 c0 {16,S}
29 H u0 p0 c0 {17,S}
30 H u0 p0 c0 {18,S}
31 H u0 p0 c0 {19,S}
32 H u0 p0 c0 {20,S}
33 H u0 p0 c0 {21,S}
34 H u0 p0 c0 {22,S}
35 H u0 p0 c0 {23,S}
36 H u0 p0 c0 {24,S}
""")
        out = generate_kekule_structure(arom)

        self.assertEqual(len(out), 1)
        self.assertFalse(out[0].isAromatic())

        bonds = set()
        for atom in out[0].atoms:
            for bond in atom.bonds.itervalues():
                bonds.add(bond)
        dBonds = 0
        for bond in bonds:
            if bond.isDouble():
                dBonds += 1

        self.assertEqual(dBonds, 12)

    def testKekulizeBridgedAromatic(self):
        """Test that we can kekulize a bridged polycyclic aromatic species."""
        arom = Molecule().fromAdjacencyList("""
1  C u0 p0 c0 {2,B} {3,S} {6,B}
2  C u0 p0 c0 {1,B} {3,B} {11,S}
3  C u0 p0 c0 {1,S} {2,B} {4,B}
4  C u0 p0 c0 {3,B} {5,B} {12,S}
5  C u0 p0 c0 {4,B} {6,B} {10,B}
6  C u0 p0 c0 {1,B} {5,B} {7,B}
7  C u0 p0 c0 {6,B} {8,B} {13,S}
8  C u0 p0 c0 {7,B} {9,B} {14,S}
9  C u0 p0 c0 {8,B} {10,B} {15,S}
10 C u0 p0 c0 {5,B} {9,B} {16,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {10,S}
""")
        out = generate_kekule_structure(arom)

        self.assertEqual(len(out), 1)
        self.assertFalse(out[0].isAromatic())

        bonds = set()
        for atom in out[0].atoms:
            for bond in atom.bonds.itervalues():
                bonds.add(bond)
        dBonds = 0
        for bond in bonds:
            if bond.isDouble():
                dBonds += 1

        self.assertEqual(dBonds, 5)

    def testKekulizeResonanceIsomer(self):
        """
        Tests that an aromatic molecule returns at least one Kekulized resonance isomer.
        
        A molecule formed using an aromatic adjacency list returns both
        the aromatic and a kekulized form as resonance isomers.
        """
        toluene = Molecule().fromAdjacencyList("""
1  H 0 {2,S}
2  C 0 {3,S} {9,S} {10,S} {1,S}
3  C 0 {4,B} {8,B} {2,S}
4  C 0 {3,B} {5,B} {11,S}
5  C 0 {4,B} {6,B} {12,S}
6  C 0 {5,B} {7,B} {13,S}
7  C 0 {6,B} {8,B} {14,S}
8  C 0 {3,B} {7,B} {15,S}
9  H 0 {2,S}
10  H 0 {2,S}
11  H 0 {4,S}
12  H 0 {5,S}
13  H 0 {6,S}
14  H 0 {7,S}
15  H 0 {8,S}""")
        
        toluene_kekulized = Molecule().fromAdjacencyList("""
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
""")
        kekulized_isomer = generate_kekule_structure(toluene)[0]
        self.assertTrue(kekulized_isomer.isIsomorphic(toluene_kekulized))

        for isomer in generate_resonance_structures(toluene):
            if isomer.isIsomorphic(toluene_kekulized):
                break
        else:  # didn't brake
            self.assertTrue(False, "Didn't find the Kekulized toulene in the result of getResonanceIsomers()")

    def testMultipleKekulizedResonanceIsomers(self):
        """Test we can make both Kekule structures of o-cresol"""

        adjlist_aromatic = """
1 C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
2 C u0 p0 c0 {1,S} {3,B} {4,B}
3 C u0 p0 c0 {2,B} {5,B} {8,S}
4 C u0 p0 c0 {2,B} {7,B} {15,S}
5 C u0 p0 c0 {3,B} {6,B} {12,S}
6 C u0 p0 c0 {5,B} {7,B} {13,S}
7 C u0 p0 c0 {4,B} {6,B} {14,S}
8 O u0 p2 c0 {3,S} {16,S}
9 H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {8,S}
"""
        molecule = Molecule().fromAdjacencyList(adjlist_aromatic)
        self.assertTrue(molecule.isAromatic(), "Starting molecule should be aromatic")
        isomers = generate_resonance_structures(molecule)
        self.assertEqual(len(isomers), 3, "Didn't generate 3 resonance isomers")
        self.assertFalse(isomers[1].isAromatic(), "Second resonance isomer shouldn't be aromatic")
        self.assertFalse(isomers[2].isAromatic(), "Third resonance isomer shouldn't be aromatic")
        self.assertFalse(isomers[1].isIsomorphic(isomers[2]), "Second and third resonance isomers should be different")

    def testMultipleKekulizedResonanceIsomersRad(self):
        """Test we can make all resonance structures of o-cresol radical"""

        adjlist_aromatic = """
1 C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
2 C u0 p0 c0 {1,S} {3,B} {4,B}
3 C u0 p0 c0 {2,B} {5,B} {8,S}
4 C u0 p0 c0 {2,B} {7,B} {15,S}
5 C u0 p0 c0 {3,B} {6,B} {12,S}
6 C u0 p0 c0 {5,B} {7,B} {13,S}
7 C u0 p0 c0 {4,B} {6,B} {14,S}
8 O u1 p2 c0 {3,S}
9 H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {4,S}
"""
        molecule = Molecule().fromAdjacencyList(adjlist_aromatic)
        self.assertTrue(molecule.isAromatic(), "Starting molecule should be aromatic")
        molList = generate_resonance_structures(molecule)
        self.assertEqual(len(molList), 6, "Expected 6 resonance structures, but generated {0}.".format(len(molList)))
        aromatic = 0
        for mol in molList:
            if mol.isAromatic():
                aromatic += 1
        self.assertEqual(aromatic, 1, "Should only have 1 aromatic resonance structure")

    @work_in_progress
    def testKekulizedResonanceIsomersFused(self):
        """Test we can make aromatic and Kekulized resonance isomers of 2-methylanthracen-1-ol
        
        This fused ring PAH will be harder"""

        kekulized1 = """multiplicity 1
1 C u0 p0 c0 {2,S} {17,S} {18,S} {19,S}
2 C u0 p0 c0 {1,S} {7,S} {10,D}
3 C u0 p0 c0 {4,S} {7,D} {9,S}
4 C u0 p0 c0 {3,S} {8,S} {11,D}
5 C u0 p0 c0 {6,S} {8,D} {12,S}
6 C u0 p0 c0 {5,S} {9,D} {13,S}
7 C u0 p0 c0 {2,S} {3,D} {16,S}
8 C u0 p0 c0 {4,S} {5,D} {22,S}
9 C u0 p0 c0 {3,S} {6,D} {27,S}
10 C u0 p0 c0 {2,D} {11,S} {20,S}
11 C u0 p0 c0 {4,D} {10,S} {21,S}
12 C u0 p0 c0 {5,S} {14,D} {23,S}
13 C u0 p0 c0 {6,S} {15,D} {26,S}
14 C u0 p0 c0 {12,D} {15,S} {24,S}
15 C u0 p0 c0 {13,D} {14,S} {25,S}
16 O u0 p2 c0 {7,S} {28,S}
17 H u0 p0 c0 {1,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {1,S}
20 H u0 p0 c0 {10,S}
21 H u0 p0 c0 {11,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {12,S}
24 H u0 p0 c0 {14,S}
25 H u0 p0 c0 {15,S}
26 H u0 p0 c0 {13,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {16,S}
"""
        kekulized2 = """multiplicity 1
1 C u0 p0 c0 {2,S} {17,S} {18,S} {19,S}
2 C u0 p0 c0 {1,S} {7,D} {10,S}
3 C u0 p0 c0 {4,S} {7,S} {9,D}
4 C u0 p0 c0 {3,S} {8,D} {11,S}
5 C u0 p0 c0 {6,S} {8,S} {12,D}
6 C u0 p0 c0 {5,S} {9,S} {13,D}
7 C u0 p0 c0 {2,D} {3,S} {16,S}
8 C u0 p0 c0 {4,D} {5,S} {22,S}
9 C u0 p0 c0 {3,D} {6,S} {27,S}
10 C u0 p0 c0 {2,S} {11,D} {20,S}
11 C u0 p0 c0 {4,S} {10,D} {21,S}
12 C u0 p0 c0 {5,D} {14,S} {23,S}
13 C u0 p0 c0 {6,D} {15,S} {26,S}
14 C u0 p0 c0 {12,S} {15,D} {24,S}
15 C u0 p0 c0 {13,S} {14,D} {25,S}
16 O u0 p2 c0 {7,S} {28,S}
17 H u0 p0 c0 {1,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {1,S}
20 H u0 p0 c0 {10,S}
21 H u0 p0 c0 {11,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {12,S}
24 H u0 p0 c0 {14,S}
25 H u0 p0 c0 {15,S}
26 H u0 p0 c0 {13,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {16,S}
"""
        kekulized3 = """multiplicity 1
1 C u0 p0 c0 {2,S} {17,S} {18,S} {19,S}
2 C u0 p0 c0 {1,S} {7,D} {10,S}
3 C u0 p0 c0 {4,S} {7,S} {9,D}
4 C u0 p0 c0 {3,S} {8,D} {11,S}
5 C u0 p0 c0 {6,D} {8,S} {12,S}
6 C u0 p0 c0 {5,D} {9,S} {13,S}
7 C u0 p0 c0 {2,D} {3,S} {16,S}
8 C u0 p0 c0 {4,D} {5,S} {20,S}
9 C u0 p0 c0 {3,D} {6,S} {21,S}
10 C u0 p0 c0 {2,S} {11,D} {22,S}
11 C u0 p0 c0 {4,S} {10,D} {23,S}
12 C u0 p0 c0 {5,S} {14,D} {24,S}
13 C u0 p0 c0 {6,S} {15,D} {25,S}
14 C u0 p0 c0 {12,D} {15,S} {26,S}
15 C u0 p0 c0 {13,D} {14,S} {27,S}
16 O u0 p2 c0 {7,S} {28,S}
17 H u0 p0 c0 {1,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {1,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {11,S}
24 H u0 p0 c0 {12,S}
25 H u0 p0 c0 {13,S}
26 H u0 p0 c0 {14,S}
27 H u0 p0 c0 {15,S}
28 H u0 p0 c0 {16,S}
"""
        kekulized4 = """multiplicity 1
1  C u0 p0 c0 {2,S} {17,S} {18,S} {19,S}
2  C u0 p0 c0 {1,S} {7,D} {10,S}
3  C u0 p0 c0 {4,D} {7,S} {9,S}
4  C u0 p0 c0 {3,D} {8,S} {11,S}
5  C u0 p0 c0 {6,S} {8,D} {12,S}
6  C u0 p0 c0 {5,S} {9,D} {13,S}
7  C u0 p0 c0 {2,D} {3,S} {16,S}
8  C u0 p0 c0 {4,S} {5,D} {20,S}
9  C u0 p0 c0 {3,S} {6,D} {21,S}
10 C u0 p0 c0 {2,S} {11,D} {22,S}
11 C u0 p0 c0 {4,S} {10,D} {23,S}
12 C u0 p0 c0 {5,S} {14,D} {24,S}
13 C u0 p0 c0 {6,S} {15,D} {25,S}
14 C u0 p0 c0 {12,D} {15,S} {26,S}
15 C u0 p0 c0 {13,D} {14,S} {27,S}
16 O u0 p2 c0 {7,S} {28,S}
17 H u0 p0 c0 {1,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {1,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {11,S}
24 H u0 p0 c0 {12,S}
25 H u0 p0 c0 {13,S}
26 H u0 p0 c0 {14,S}
27 H u0 p0 c0 {15,S}
28 H u0 p0 c0 {16,S}
"""
        m1 = Molecule().fromAdjacencyList(kekulized1)
        m2 = Molecule().fromAdjacencyList(kekulized2)
        m3 = Molecule().fromAdjacencyList(kekulized3)
        m4 = Molecule().fromAdjacencyList(kekulized4)
        resonance_forms = (m1, m2, m3, m4)

        for starting in resonance_forms:
            self.assertFalse(starting.isAromatic(), "Starting molecule should not be aromatic")

            isomers = generate_resonance_structures(starting)
            # print "starting with {0!r} I generated these:".format(starting)
            # print repr(isomers)
            for isomer in isomers:
                if isomer.isAromatic():
                    break
            else:  # didn't break
                self.fail("None of the generated resonance isomers {0!r} are aromatic".format(isomers))

            for generated in isomers:
                for expected in resonance_forms:
                    if generated.isIsomorphic(expected):
                        break
                else:  # didn't break
                    if generated.isAromatic():
                        continue  # because the aromatic isomer isn't in our resonance_forms list
                    self.fail("Generated a resonance form {0!r} that was not expected!\n{1}\nAlthough that may be a bug in the unit test (not sure I got them all)".format(generated, generated.toAdjacencyList()))

            for expected in resonance_forms:
                for generated in isomers:
                    if expected.isIsomorphic(generated):
                        break
                else:  # didn't break
                    self.fail(("Expected a resonance form {0!r} that was not generated.\n"
                              "Only generated these:\n{1}").format(expected, '\n'.join([repr(g) for g in isomers])))

    def testKeepIsomorphicStructuresFunctionsWhenTrue(self):
        """Test that keep_isomorphic works for resonance structure generation when True."""
        mol = Molecule(SMILES='C=C[CH2]')
        mol.assignAtomIDs()
        out = generate_resonance_structures(mol, keep_isomorphic=True)

        self.assertEqual(len(out), 2)
        self.assertTrue(out[0].isIsomorphic(out[1]))
        self.assertFalse(out[0].isIdentical(out[1]))

    def testKeepIsomorphicStructuresFunctionsWhenFalse(self):
        """Test that keep_isomorphic works for resonance structure generation when False."""
        mol = Molecule(SMILES='C=C[CH2]')
        mol.assignAtomIDs()
        out = generate_resonance_structures(mol, keep_isomorphic=False)

        self.assertEqual(len(out), 1)

    def testFalseNegativeAromaticityPerception(self):
        """Test that we obtain the correct aromatic structure for a monocyclic aromatic that RDKit mis-identifies."""
        mol = Molecule(SMILES='[CH2]C=C1C=CC(=C)C=C1')
        out = generate_resonance_structures(mol)

        aromatic = Molecule().fromAdjacencyList("""
multiplicity 2
1  C u0 p0 c0 {4,B} {5,B} {7,S}
2  C u0 p0 c0 {3,B} {6,B} {8,S}
3  C u0 p0 c0 {2,B} {4,B} {10,S}
4  C u0 p0 c0 {1,B} {3,B} {11,S}
5  C u0 p0 c0 {1,B} {6,B} {13,S}
6  C u0 p0 c0 {2,B} {5,B} {14,S}
7  C u0 p0 c0 {1,S} {9,D} {12,S}
8  C u1 p0 c0 {2,S} {15,S} {16,S}
9  C u0 p0 c0 {7,D} {17,S} {18,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {9,S}
""")

        self.assertEqual(len(out), 5)
        self.assertTrue(any([m.isIsomorphic(aromatic) for m in out]))

    def testFalseNegativePolycyclicAromaticityPerception(self):
        """Test that we generate proper structures for a polycyclic aromatic that RDKit mis-identifies."""
        mol = Molecule(SMILES='C=C1C=CC=C2C=C[CH]C=C12')
        out = generate_resonance_structures(mol)

        clar = Molecule().fromAdjacencyList("""
multiplicity 2
1  C u0 p0 c0 {2,B} {3,B} {7,S}
2  C u0 p0 c0 {1,B} {5,B} {6,S}
3  C u0 p0 c0 {1,B} {4,B} {11,S}
4  C u0 p0 c0 {3,B} {8,B} {13,S}
5  C u0 p0 c0 {2,B} {8,B} {15,S}
6  C u0 p0 c0 {2,S} {9,D} {16,S}
7  C u0 p0 c0 {1,S} {10,D} {18,S}
8  C u0 p0 c0 {4,B} {5,B} {14,S}
9  C u0 p0 c0 {6,D} {10,S} {17,S}
10 C u0 p0 c0 {7,D} {9,S} {12,S}
11 C u1 p0 c0 {3,S} {19,S} {20,S}
12 H u0 p0 c0 {10,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {11,S}
20 H u0 p0 c0 {11,S}
""")

        self.assertEqual(len(out), 6)
        self.assertTrue(any([m.isIsomorphic(clar) for m in out]))

    def testFalseNegativePolycylicAromaticityPerception2(self):
        """Test that we obtain the correct aromatic structure for a polycylic aromatic that RDKit mis-identifies."""
        mol = Molecule(SMILES='[CH2]C=C1C=CC(=C)C2=C1C=CC=C2')
        out = generate_resonance_structures(mol)

        aromatic = Molecule().fromAdjacencyList("""
multiplicity 2
1 C u0 p0 c0 {2,B} {4,B} {8,B}
2 C u0 p0 c0 {1,B} {3,B} {7,B}
3 C u0 p0 c0 {2,B} {5,B} {9,S}
4 C u0 p0 c0 {1,B} {6,B} {12,S}
5 C u0 p0 c0 {3,B} {6,B} {15,S}
6 C u0 p0 c0 {4,B} {5,B} {16,S}
7 C u0 p0 c0 {2,B} {10,B} {17,S}
8 C u0 p0 c0 {1,B} {11,B} {20,S}
9 C u0 p0 c0 {3,S} {13,D} {14,S}
10 C u0 p0 c0 {7,B} {11,B} {18,S}
11 C u0 p0 c0 {8,B} {10,B} {19,S}
12 C u1 p0 c0 {4,S} {21,S} {22,S}
13 C u0 p0 c0 {9,D} {23,S} {24,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {11,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {12,S}
22 H u0 p0 c0 {12,S}
23 H u0 p0 c0 {13,S}
24 H u0 p0 c0 {13,S}
""")

        self.assertEqual(len(out), 7)
        self.assertTrue(any([m.isIsomorphic(aromatic) for m in out]))


class ClarTest(unittest.TestCase):
    """
    Contains unit tests for Clar structure methods.
    """

    def testClarTransformation(self):
        """Test that clarTransformation generates an aromatic ring."""
        mol = Molecule().fromSMILES('c1ccccc1')
        sssr = mol.getSmallestSetOfSmallestRings()
        _clar_transformation(mol, sssr[0])
        mol.updateAtomTypes()

        self.assertTrue(mol.isAromatic())

    def testClarOptimization(self):
        """Test to ensure pi electrons are conserved during optimization"""
        mol = Molecule().fromSMILES('C1=CC=C2C=CC=CC2=C1')  # Naphthalene
        output = _clar_optimization(mol)

        for molecule, asssr, bonds, solution in output:

            # Count pi electrons in molecule
            pi = 0
            for bond in bonds:
                if bond.isDouble():
                    pi += 2

            # Count pi electrons in solution
            y = solution[0:len(asssr)]
            x = solution[len(asssr):]
            pi_solution = 6 * sum(y) + 2 * sum(x)

            # Check that both counts give 10 pi electrons
            self.assertEqual(pi, 10)
            self.assertEqual(pi_solution, 10)

            # Check that we only assign 1 aromatic sextet
            self.assertEqual(sum(y), 1)

    def testPhenanthrene(self):
        """Test that we generate 1 Clar structure for phenanthrene."""
        mol = Molecule().fromSMILES('C1=CC=C2C(C=CC3=CC=CC=C32)=C1')
        newmol = generate_clar_structures(mol)

        struct = Molecule().fromAdjacencyList("""1  C u0 p0 c0 {2,S} {3,B} {5,B}
2  C u0 p0 c0 {1,S} {4,B} {9,B}
3  C u0 p0 c0 {1,B} {6,S} {10,B}
4  C u0 p0 c0 {2,B} {7,S} {8,B}
5  C u0 p0 c0 {1,B} {12,B} {17,S}
6  C u0 p0 c0 {3,S} {7,D} {18,S}
7  C u0 p0 c0 {4,S} {6,D} {19,S}
8  C u0 p0 c0 {4,B} {13,B} {20,S}
9  C u0 p0 c0 {2,B} {14,B} {23,S}
10 C u0 p0 c0 {3,B} {11,B} {24,S}
11 C u0 p0 c0 {10,B} {12,B} {15,S}
12 C u0 p0 c0 {5,B} {11,B} {16,S}
13 C u0 p0 c0 {8,B} {14,B} {21,S}
14 C u0 p0 c0 {9,B} {13,B} {22,S}
15 H u0 p0 c0 {11,S}
16 H u0 p0 c0 {12,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {13,S}
22 H u0 p0 c0 {14,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {10,S}
""")

        self.assertEqual(len(newmol), 1)
        self.assertTrue(newmol[0].isIsomorphic(struct))

    def testPhenalene(self):
        """Test that we generate 2 Clar structures for phenalene.

        Case where there is one non-aromatic ring."""
        mol = Molecule().fromSMILES('C1=CC2=CC=CC3CC=CC(=C1)C=32')
        newmol = generate_clar_structures(mol)

        struct1 = Molecule().fromAdjacencyList("""1  C u0 p0 c0 {2,S} {6,S} {14,S} {15,S}
2  C u0 p0 c0 {1,S} {3,S} {7,D}
3  C u0 p0 c0 {2,S} {4,B} {5,B}
4  C u0 p0 c0 {3,B} {9,B} {10,S}
5  C u0 p0 c0 {3,B} {8,S} {11,B}
6  C u0 p0 c0 {1,S} {8,D} {16,S}
7  C u0 p0 c0 {2,D} {13,S} {21,S}
8  C u0 p0 c0 {5,S} {6,D} {22,S}
9  C u0 p0 c0 {4,B} {12,B} {18,S}
10 C u0 p0 c0 {4,S} {13,D} {19,S}
11 C u0 p0 c0 {5,B} {12,B} {23,S}
12 C u0 p0 c0 {9,B} {11,B} {17,S}
13 C u0 p0 c0 {7,S} {10,D} {20,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {1,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {12,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {10,S}
20 H u0 p0 c0 {13,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {11,S}
""")
        struct2 = Molecule().fromAdjacencyList("""1  C u0 p0 c0 {2,S} {6,S} {14,S} {15,S}
2  C u0 p0 c0 {1,S} {3,B} {7,B}
3  C u0 p0 c0 {2,B} {4,B} {5,S}
4  C u0 p0 c0 {3,B} {9,S} {10,B}
5  C u0 p0 c0 {3,S} {8,S} {11,D}
6  C u0 p0 c0 {1,S} {8,D} {16,S}
7  C u0 p0 c0 {2,B} {13,B} {21,S}
8  C u0 p0 c0 {5,S} {6,D} {22,S}
9  C u0 p0 c0 {4,S} {12,D} {18,S}
10 C u0 p0 c0 {4,B} {13,B} {19,S}
11 C u0 p0 c0 {5,D} {12,S} {23,S}
12 C u0 p0 c0 {9,D} {11,S} {17,S}
13 C u0 p0 c0 {7,B} {10,B} {20,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {1,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {12,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {10,S}
20 H u0 p0 c0 {13,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {11,S}
""")

        self.assertEqual(len(newmol), 2)
        self.assertTrue(newmol[0].isIsomorphic(struct1) or newmol[0].isIsomorphic(struct2))
        self.assertTrue(newmol[1].isIsomorphic(struct2) or newmol[1].isIsomorphic(struct1))
        self.assertFalse(newmol[0].isIsomorphic(newmol[1]))

    def testCorannulene(self):
        """Test that we generate 5 Clar structures for corannulene

        Case where linear relaxation does not give an integer solution"""
        mol = Molecule().fromSMILES('C1=CC2=CC=C3C=CC4=C5C6=C(C2=C35)C1=CC=C6C=C4')
        newmol = generate_clar_structures(mol)

        struct = Molecule().fromAdjacencyList("""1  C u0 p0 c0 {2,S} {5,B} {8,B}
2  C u0 p0 c0 {1,S} {3,B} {10,B}
3  C u0 p0 c0 {2,B} {4,S} {9,B}
4  C u0 p0 c0 {3,S} {5,S} {6,D}
5  C u0 p0 c0 {1,B} {4,S} {7,B}
6  C u0 p0 c0 {4,D} {12,S} {13,S}
7  C u0 p0 c0 {5,B} {14,S} {15,B}
8  C u0 p0 c0 {1,B} {16,B} {20,S}
9  C u0 p0 c0 {3,B} {11,S} {17,B}
10 C u0 p0 c0 {2,B} {18,B} {19,S}
11 C u0 p0 c0 {9,S} {12,D} {21,S}
12 C u0 p0 c0 {6,S} {11,D} {22,S}
13 C u0 p0 c0 {6,S} {14,D} {23,S}
14 C u0 p0 c0 {7,S} {13,D} {24,S}
15 C u0 p0 c0 {7,B} {16,B} {25,S}
16 C u0 p0 c0 {8,B} {15,B} {26,S}
17 C u0 p0 c0 {9,B} {18,B} {27,S}
18 C u0 p0 c0 {10,B} {17,B} {28,S}
19 C u0 p0 c0 {10,S} {20,D} {29,S}
20 C u0 p0 c0 {8,S} {19,D} {30,S}
21 H u0 p0 c0 {11,S}
22 H u0 p0 c0 {12,S}
23 H u0 p0 c0 {13,S}
24 H u0 p0 c0 {14,S}
25 H u0 p0 c0 {15,S}
26 H u0 p0 c0 {16,S}
27 H u0 p0 c0 {17,S}
28 H u0 p0 c0 {18,S}
29 H u0 p0 c0 {19,S}
30 H u0 p0 c0 {20,S}
""")

        self.assertEqual(len(newmol), 5)
        self.assertTrue(newmol[0].isIsomorphic(struct))
        self.assertTrue(newmol[1].isIsomorphic(struct))
        self.assertTrue(newmol[2].isIsomorphic(struct))
        self.assertTrue(newmol[3].isIsomorphic(struct))
        self.assertTrue(newmol[4].isIsomorphic(struct))

    def testExocyclicDB(self):
        """Test that Clar structure generation doesn't modify exocyclic double bonds

        Important for cases where RDKit considers rings to be aromatic by counting pi-electron contributions
        from exocyclic double bonds, while they don't actually contribute to aromaticity"""

        mol = Molecule(SMILES="C=C1C=CC=CC1=C")
        newmol = generate_clar_structures(mol)

        self.assertEquals(len(newmol), 0)
