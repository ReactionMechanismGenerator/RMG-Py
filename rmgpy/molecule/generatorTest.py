import re
import unittest
from external.wip import work_in_progress

from rmgpy.species import Species
from .molecule import Atom, Molecule
from .inchi import P_LAYER_PREFIX, U_LAYER_PREFIX
from .generator import *

class CreateULayerTest(unittest.TestCase):
    def testC4H6(self):
        """
        Test that 3-butene-1,2-diyl biradical is always resulting in the 
        same u-layer, regardless of the original order.
        """

        # radical positions 3 and 4
        adjlist1 = """
1  C u0 p0 c0 {2,D} {5,S} {6,S}
2  C u0 p0 c0 {1,D} {3,S} {7,S}
3  C u1 p0 c0 {2,S} {4,S} {8,S}
4  C u1 p0 c0 {3,S} {9,S} {10,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}

        """        

        # radical positions 1 and 2
        adjlist2 = """
1  C u1 p0 c0 {2,S} {5,S} {6,S}
2  C u1 p0 c0 {1,S} {3,S} {7,S}
3  C u0 p0 c0 {2,S} {4,D} {8,S}
4  C u0 p0 c0 {3,D} {9,S} {10,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
        """

        u_layers = []
        for adjlist in [adjlist1, adjlist2]:
            mol = Molecule().fromAdjacencyList(adjlist)
            u_layer = create_augmented_layers(mol)[0]
            u_layers.append(u_layer)

        self.assertEquals(u_layers[0], u_layers[1])

class InChIGenerationTest(unittest.TestCase):
    def compare(self, adjlist, aug_inchi):
        spc = Species(molecule=[Molecule().fromAdjacencyList(adjlist)])
        spc.generateResonanceIsomers()

        ignore_prefix = r"(InChI=1+)(S*)/"

        exp = re.split(ignore_prefix, aug_inchi)[-1]
        comp = re.split(ignore_prefix, spc.getAugmentedInChI())[-1]
        self.assertEquals(exp, comp)

    def test_C5H5(self):
        """
        Test that the unpaired electron of 1,3-cyclopentadienyl radical always
        ends up on the 1-carbon atom.
        """

        adjlist = """
1 C 0 {2,D} {5,S}
2 C 0 {1,D} {3,S} 
3 C 0 {2,S} {4,D} 
4 C 0 {3,D} {5,S} 
5 C 1 {4,S} {1,S}
        """

        aug_inchi = 'InChI=1S/C5H5/c1-2-4-5-3-1/h1-5H/u1'
        self.compare(adjlist, aug_inchi)


    def test_C7H8(self):
        """Looks a lot like toluene but with 1 double bond replaced by a biradical."""

        """unpaired electrons on tertiary carbon, and on carbon in para position."""
        adjlist = """
1  C u1 p0 c0 {2,S} {3,S} {4,S}
2  C u0 p0 c0 {1,S} {7,D} {10,S}
3  C u0 p0 c0 {1,S} {6,D} {11,S}
4  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
5  C u1 p0 c0 {6,S} {7,S} {15,S}
6  C u0 p0 c0 {3,D} {5,S} {8,S}
7  C u0 p0 c0 {2,D} {5,S} {9,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
        """

        aug_inchi = 'InChI=1S/C7H8/c1-7-5-3-2-4-6-7/h2-6H,1H3/u2,3'
        self.compare(adjlist, aug_inchi)
    
    def test_C8H8(self):
        """Looks a lot like cycloctene but with 1 double bond replaced by a biradical."""

        adjlist = """
1  C u0 p0 c0 {2,S} {5,D} {9,S}
2  C u0 p0 c0 {1,S} {3,D} {10,S}
3  C u0 p0 c0 {2,D} {4,S} {11,S}
4  C u1 p0 c0 {3,S} {6,S} {12,S}
5  C u0 p0 c0 {1,D} {8,S} {14,S}
6  C u1 p0 c0 {4,S} {7,S} {15,S}
7  C u0 p0 c0 {6,S} {8,D} {13,S}
8  C u0 p0 c0 {5,S} {7,D} {16,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {8,S}
        """

        aug_inchi = 'InChI=1S/C8H8/c1-2-4-6-8-7-5-3-1/h1-8H/u1,2'
        self.compare(adjlist, aug_inchi)

    def test_benzyne(self):

        adjlist = """
1  C u0 p0 c0 {2,T} {6,S}
2  C u0 p0 c0 {1,T} {3,S}
3  C u0 p0 c0 {2,S} {4,D} {7,S}
4  C u0 p0 c0 {3,D} {5,S} {8,S}
5  C u0 p0 c0 {4,S} {6,D} {9,S}
6  C u0 p0 c0 {1,S} {5,D} {10,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
        """
        benzatetraene = 'InChI=1S/C6H4/c1-2-4-6-5-3-1/h1-4H'
        aug_inchi = 'InChI=1S/C6H4/c1-2-4-6-5-3-1/h1-4H'
        self.compare(adjlist, aug_inchi)

    def test_H(self):
        adjlist = """
multiplicity 2
1 H u1 p0 c0
"""
        aug_inchi = 'InChI=1S/H/u1'
        self.compare(adjlist, aug_inchi)


    def test_C6H8(self):
        """
        Test that the 2 unpaired electrons of .CC(=C)C(C.)=C
        do not end up at the same side of the central C-C bond.
        """
        adjlist = """
1 C 0 {2,D}
2 C 0 {1,D} {3,S} {4,S}
3 C 1 {2,S}
4 C 0 {2,S} {5,S} {6,D}
5 C 1 {4,S}
6 C 0 {4,D}
        """

        aug_inchi = 'InChI=1S/C6H8/c1-5(2)6(3)4/h1-4H2/u1,3'
        self.compare(adjlist, aug_inchi)


    def test_C6H10_tetrarad(self):
        adjlist = """
1  C u1 p0 c0 {3,S} {7,S} {8,S}
2  C u1 p0 c0 {4,S} {9,S} {10,S}
3  C u1 p0 c0 {1,S} {5,S} {11,S}
4  C u1 p0 c0 {2,S} {6,S} {12,S}
5  C u0 p0 c0 {3,S} {6,S} {13,S} {14,S}
6  C u0 p0 c0 {4,S} {5,S} {15,S} {16,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
        """

        aug_inchi = 'InChI=1S/C6H10/c1-3-5-6-4-2/h3-4H,1-2,5-6H2/u1,2,3,4'
        self.compare(adjlist, aug_inchi)

    def test_Buta13diyl_triplet(self):
        """
        C=CC.C.
        """
        adjlist = """
        multiplicity 3
1  C u1 p0 c0 {2,S} {5,S} {6,S}
2  C u1 p0 c0 {1,S} {3,S} {7,S}
3  C u0 p0 c0 {2,S} {4,D} {8,S}
4  C u0 p0 c0 {3,D} {9,S} {10,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
"""

        aug_inchi = 'InChI=1S/C4H6/c1-3-4-2/h3-4H,1-2H2/u1,2'
        self.compare(adjlist, aug_inchi)

    def test_CH2O2(self):

        adjlist = """
1 C 1 {2,S} {3,S}
2 O 0 {1,S}
3 O 1 {1,S}
"""

        aug_inchi = 'InChI=1/CH2O2/c2-1-3/h1H,(H,2,3)/u1,2'
        self.compare(adjlist, aug_inchi)

    def test_C7H10(self):
        adjlist = """

        1 C 1 {2,S}
2 C 0 {1,S} {3,D} {4,S}
3 C 0 {2,D}
4 C 0 {2,S} {5,S}
5 C 1 {4,S} {6,S} {7,S}
6 C 1 {5,S}
7 C 1 {5,S}
"""

        aug_inchi = 'InChI=1S/C7H10/c1-6(2)5-7(3)4/h1-5H2/u1,2,3,6'
        self.compare(adjlist, aug_inchi)

    def test_C5H6O(self):

        adjlist = """
1 C 1 {2,S}
2 C 0 {1,S} {3,D}
3 C 0 {2,D} {4,S} {5,S}
4 O 1 {3,S}
5 C 0 {3,S} {6,D}
6 C 0 {5,D}
"""

        aug_inchi = 'InChI=1S/C5H6O/c1-3-5(6)4-2/h3-4H,1-2H2/u1,3'
        self.compare(adjlist, aug_inchi)

    def test_C7H9(self):

        adjlist = """
1 C 0 {4,D} 
2 C 0 {5,D}
3 C 1 {6,S}
4 C 0 {1,D} {7,S}
5 C 0 {2,D} {7,S}
6 C 1 {3,S} {7,S}
7 C 1 {4,S} {5,S} {6,S}
"""

        aug_inchi = 'InChI=1S/C7H9/c1-4-7(5-2)6-3/h4-6H,1-3H2/u1,4,7'
        self.compare(adjlist, aug_inchi)

    def test_C7H9(self):

        adjlist = """
1 C 0 {5,D}
2 C 1 {6,S}
3 C 1 {7,S}
4 C 0 {8,D}
5 C 0 {1,D} {9,S}
6 C 1 {2,S} {10,S}
7 C 1 {3,S} {11,S}
8 C 0 {4,D} {11,S}
9 C 0 {5,S} {11,S}
10 C 0 {6,S} {11,S}
11 C 0 {7,S} {8,S} {9,S} {10,S}
"""

        aug_inchi = 'InChI=1S/C11H16/c1-5-9-11(7-3,8-4)10-6-2/h5-8H,1-4,9-10H2/u1,3,5,7'
        self.compare(adjlist, aug_inchi)

    def test_singlet_vs_closed_shell(self):
        adjlist_singlet = """
1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 C u0 p0 c0 {1,D} {3,S} {5,S}
3 C u0 p1 c0 {1,S} {2,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
        """

        adjlist_closed_shell = """
1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 C u0 p0 c0 {1,D} {3,D}
3 C u0 p0 c0 {1,S} {2,D} {5,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {3,S}
        """

        singlet = Species(molecule=[Molecule().fromAdjacencyList(adjlist_singlet)])
        singlet.generateResonanceIsomers()
        closed_shell = Species(molecule=[Molecule().fromAdjacencyList(adjlist_closed_shell)])
        closed_shell.generateResonanceIsomers()

        singlet_aug_inchi = singlet.getAugmentedInChI()
        closed_shell_aug_inchi = closed_shell.getAugmentedInChI()
        self.assertTrue(singlet_aug_inchi != closed_shell_aug_inchi)

    def test_C6H5(self):
        """Test that the u-layer of phenyl shows atom 1."""
        adjlist = """
multiplicity 2
1  C u0 p0 c0 {2,D} {3,S} {10,S}
2  C u0 p0 c0 {1,D} {5,S} {7,S}
3  C u0 p0 c0 {1,S} {6,D} {8,S}
4  C u0 p0 c0 {5,D} {6,S} {11,S}
5  C u0 p0 c0 {2,S} {4,D} {9,S}
6  C u1 p0 c0 {3,D} {4,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {4,S}
"""

        aug_inchi = 'InChI=1S/C6H5/c1-2-4-6-5-3-1/h1-5H/u1'
        self.compare(adjlist, aug_inchi)

class ExpectedLonePairsTest(unittest.TestCase):

    def test_SingletCarbon(self):
        mol = Molecule(atoms=[Atom(element='C', lonePairs=1)])
        unexpected = has_unexpected_lone_pairs(mol)
        self.assertTrue(unexpected)

    def test_NormalCarbon(self):
        mol = Molecule(atoms=[Atom(element='C', lonePairs=0)])
        unexpected = has_unexpected_lone_pairs(mol)
        self.assertFalse(unexpected)

    def test_NormalOxygen(self):
        mol = Molecule(atoms=[Atom(element='O', lonePairs=2)])
        unexpected = has_unexpected_lone_pairs(mol)
        self.assertFalse(unexpected)

    def test_Oxygen_3LP(self):
        mol = Molecule(atoms=[Atom(element='O', lonePairs=3)])
        unexpected = has_unexpected_lone_pairs(mol)
        self.assertTrue(unexpected)

class CreateAugmentedLayersTest(unittest.TestCase):
    def test_Methane(self):
        smi = 'C'
        mol = Molecule().fromSMILES(smi)
        ulayer, player = create_augmented_layers(mol)
        self.assertTrue(not ulayer)
        self.assertTrue(not player)

    def test_SingletMethylene(self):
        adjlist = """
multiplicity 1
1 C u0 p1 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
"""
        mol = Molecule().fromAdjacencyList(adjlist)
        ulayer, player = create_augmented_layers(mol)
        self.assertTrue(not ulayer)
        self.assertEquals(P_LAYER_PREFIX + '1', player)

    def test_TripletMethylene(self):
        adjlist = """
multiplicity 3
1 C u2 p0 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
"""
        mol = Molecule().fromAdjacencyList(adjlist)
        ulayer, player = create_augmented_layers(mol)
        self.assertEquals(U_LAYER_PREFIX + '1,1', ulayer)
        self.assertTrue(not player)

    @work_in_progress
    def test_Nitrate(self):
        """
        Test that N atom in the p-layer has correct symbol.
        """
        
        adjlist = """
1 O u0 p2 c0 {4,D}
2 O u0 p3 c-1 {4,S}
3 O u0 p3 c-1 {4,S}
4 N u0 p0 c+1 {1,D} {2,S} {3,S}
"""
        mol = Molecule().fromAdjacencyList(adjlist)
        ulayer, player = create_augmented_layers(mol)
        self.assertTrue(not ulayer)
        self.assertTrue(player.contains(P_LAYER_PREFIX + '1(0)'))



if __name__ == '__main__':
    unittest.main()