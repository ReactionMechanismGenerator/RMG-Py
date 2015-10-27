import unittest
from external.wip import work_in_progress

from .molecule import Molecule

from .resonance import *

class ResonanceTest(unittest.TestCase):

    def test_C9H9_aro(self):
        """CyclopropylBenzene-radical, aromatic bonds"""
        mol = Molecule(SMILES="[CH]1CC1c1ccccc1")
        generateResonanceIsomers(mol)
    
    def test_C9H9_kek(self):
        """CyclopropylBenzene-radical, kekulized bonds"""
        mol = Molecule(SMILES="[CH]1CC1C1C=CC=CC=1")
        generateResonanceIsomers(mol)

    def test_Benzene_aro(self):
        mol = Molecule(SMILES="c1ccccc1")
        generateResonanceIsomers(mol)
    
    def test_Benzene_kek(self):
        mol = Molecule(SMILES="C1C=CC=CC=1")
        generateResonanceIsomers(mol)

    def test_C9H11_aro(self):
        """PropylBenzene-radical"""
        mol = Molecule(SMILES="[CH2]CCc1ccccc1")
        generateResonanceIsomers(mol)

    def test_C10H11_aro(self):
        """CyclobutylBenzene-radical"""
        mol = Molecule(SMILES="[CH]1CCC1c1ccccc1")
        generateResonanceIsomers(mol)

    def test_C9H10_aro(self):
        """CyclopropylBenzene, aromatic bonds"""
        mol = Molecule(SMILES="C1CC1c1ccccc1")
        generateResonanceIsomers(mol)

    def test_C10H12_aro(self):
        """CyclopropylMethylBenzene"""
        mol = Molecule(SMILES="C1CC1c1c(C)cccc1")
        generateResonanceIsomers(mol)

    def test_C9H10_aro(self):
        """CyclopropylBenzene, generate aro resonance isomers"""
        mol = Molecule(SMILES="C1CC1c1ccccc1")
        getAromaticResonanceIsomers(mol)

    @work_in_progress
    def testMultipleKekulizedResonanceIsomers(self):
        "Test we can make both Kekulized resonance isomers of 2-Hydroxy-1-methylbenzene"

        adjlist_aromatic = """multiplicity 1
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
        kekulized_isomer = getKekulizedResonanceIsomers(toluene)[0]
        self.assertTrue(kekulized_isomer.isIsomorphic(toluene_kekulized))

        for isomer in generateResonanceIsomers(toluene):
            if isomer.isIsomorphic(toluene_kekulized):
                break
        else:  # didn't brake
            self.assertTrue(False, "Didn't find the Kekulized toulene in the result of getResonanceIsomers()")

    @work_in_progress
    def testMultipleKekulizedResonanceIsomers(self):
        "Test we can make both Kekulized resonance isomers of 2-Hydroxy-1-methylbenzene"

        adjlist_aromatic = """multiplicity 1
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
        isomers = generateResonanceIsomers(molecule)
        self.assertEqual(len(isomers), 3, "Didn't generate 3 resonance isomers")
        self.assertFalse(isomers[1].isAromatic(), "Second resonance isomer shouldn't be aromatic")
        self.assertFalse(isomers[2].isAromatic(), "Third resonance isomer shouldn't be aromatic")
        self.assertFalse(isomers[1].isIsomorphic(isomers[2]), "Second and third resonance isomers should be different")

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

            isomers = generateResonanceIsomers(starting)
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
