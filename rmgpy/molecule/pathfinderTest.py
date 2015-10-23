import unittest

from rmgpy.molecule import Molecule
from rmgpy.molecule.pathfinder import *


class FindButadieneTest(unittest.TestCase):
    def test_13butadiene(self):

        mol = Molecule().fromSMILES("C=CC=C")#1,3-butadiene

        start, end = mol.atoms[0], mol.atoms[3]
        path = find_butadiene(start, end)
        self.assertIsNotNone(path)

    def test_acrolein(self):
        mol = Molecule().fromSMILES("C=CC=O")#Acrolein

        start, end = mol.atoms[0], mol.atoms[3]
        path = find_butadiene(start, end)
        self.assertIsNotNone(path)

        start, end = mol.atoms[0], mol.atoms[4]#wrong end
        path = find_butadiene(start, end)
        self.assertIsNone(path)

        start, end = mol.atoms[-1], mol.atoms[3]#wrong start
        path = find_butadiene(start, end)
        self.assertIsNone(path)
    
    def test_135hexatriene(self):
        mol = Molecule().fromSMILES("C=CC=CC=C")#1,3,5-hexatriene

        start, end = mol.atoms[0], mol.atoms[5]
        path = find_butadiene(start, end)
        self.assertIsNotNone(path)        

    def test_13cyclohexadiene(self):

        adjlist = """
1  C u0 p0 c0 {2,D} {6,S} {7,S}
2  C u0 p0 c0 {1,D} {3,S} {8,S}
3  C u0 p0 c0 {2,S} {4,D} {9,S}
4  C u0 p0 c0 {3,D} {5,S} {10,S}
5  C u0 p0 c0 {4,S} {6,S} {11,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {13,S} {14,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
        """
        mol = Molecule().fromAdjacencyList(adjlist)#1,3-cyclohexadiene

        start, end = mol.atoms[0], mol.atoms[3]
        path = find_butadiene(start, end)
        self.assertIsNotNone(path)     


    def test_14cyclohexadiene(self):

        adjlist = """
1  C u0 p0 c0 {2,D} {6,S} {7,S}
2  C u0 p0 c0 {1,D} {3,S} {8,S}
3  C u0 p0 c0 {2,S} {4,S} {9,S} {10,S}
4  C u0 p0 c0 {3,S} {5,D} {11,S}
5  C u0 p0 c0 {4,D} {6,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {13,S} {14,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
        """

        mol = Molecule().fromAdjacencyList(adjlist)#1,4-cyclohexadiene

        start, end = mol.atoms[0], mol.atoms[3]
        path = find_butadiene(start, end)
        self.assertIsNone(path)   

    def test_Benzene(self):
        mol = Molecule().fromSMILES("C1=CC=CC=C1")#benzene

        start, end = mol.atoms[0], mol.atoms[5]
        path = find_butadiene(start, end)
        self.assertIsNotNone(path)    
    
    def test_C4H4(self):
        mol = Molecule().fromSMILES("C=C=C=C")#C4H4

        start, end = mol.atoms[0], mol.atoms[3]
        path = find_butadiene(start, end)
        self.assertIsNotNone(path)    
