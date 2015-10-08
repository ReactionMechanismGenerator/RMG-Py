import unittest

from rmgpy.molecule import Molecule
from rmgpy.molecule.parser import fromAugmentedInChI

class InChITest(unittest.TestCase):

    def compare(self, inchi, mult, u_indices = None):
        u_layer = ','.join([str(i) for i in u_indices]) if u_indices else None

        if u_layer:
            aug_inchi = 'InChI=1/' + inchi  + '/mult' + str(mult) + '/u' + u_layer
        else: 
            aug_inchi = 'InChI=1/' + inchi  + '/mult' + str(mult)

        mol = Molecule()
        mol = fromAugmentedInChI(mol, aug_inchi)
        self.assertEqual(mol.getNumberOfRadicalElectrons(), mult - 1)
        return mol

    def test_Ethane_parsing(self):
        inchi = 'C2H6/c1-2/h1-2H3'
        mult = 1

        aug_inchi = 'InChI=1/' + inchi  + '/mult' + str(mult)

        # assert aug_inchi == '', aug_inchi
        self.compare(inchi, mult)
        
    def test_Ethyl_parsing(self):
        inchi = 'C2H5/c1-2/h1H2,2H3'
        mult = 2
        self.compare(inchi, mult)

    def test_CH3_parsing(self):
        inchi = 'CH3/h1H3'
        mult = 2
        self.compare(inchi, mult)

    def test_H2_parsing(self):
        inchi = 'H2/h1H'
        mult = 1
        self.compare(inchi, mult)

    def test_C2H4_biradical_parsing(self):
        inchi = 'C2H4/c1-2/h1-2H2'
        mult = 3
        u_indices = [1,2]
        self.compare(inchi, mult, u_indices)

    def test_C2H3_triradical_parsing(self):
        inchi = 'C2H3/c1-2/h1H,2H2'
        mult = 4
        u_indices = [1,1,2]
        self.compare(inchi, mult, u_indices)

    def test_C3H6_biradical_parsing(self):
        inchi = 'C3H6/c1-3-2/h1-3H2'
        mult = 3
        u_indices = [1,2]
        self.compare(inchi, mult, u_indices)

    def testC2H3O3(self):
        adjlist = '''
        1 C u0 p0 c0 {2,D} {6,S} {7,S}
        2 C u0 p0 c0 {1,D} {3,S} {5,S}
        3 O u1 p2 c0 {2,S}
        4 O u0 p2 c0 {5,S} {8,S}
        5 O u0 p2 c0 {2,S} {4,S}
        6 H u0 p0 c0 {1,S}
        7 H u0 p0 c0 {1,S}
        8 H u0 p0 c0 {4,S}
        '''
        inchi = 'C2H3O3/c1-2(3)5-4/h4H,1H2'
        mult = 2
        self.compare(inchi, mult)

    def testC2H2(self):
        inchi = 'C2H2/c1-2/h1-2H'
        mult = 3
        u_indices = [1,2]
        mol = self.compare(inchi, mult,  u_indices)

    def testO2(self):
        inchi = 'O2/c1-2'
        mult = 3
        u_indices = [1,2]
        self.compare(inchi, mult, u_indices)

    def testTriRadicalZwitterMult4(self):
        inchi = 'C6H11/c1-3-5-6-4-2/h5H,1-4,6H2'
        mult = 4
        u_indices = [1,3,2]
        self.compare(inchi, mult, u_indices)

    def testTriRadicalDoubleBondMult4(self):
        inchi = 'C4H7/c1-3-4-2/h3H,1-2,4H2'
        mult = 4
        u_indices = [1,2,3]
        self.compare(inchi, mult, u_indices)

    def testTriRadical2DoubleBondMult4(self):
        inchi = 'C6H9/c1-4-6(3)5-2/h1,4-6H,2H2,3H3'
        mult = 4
        u_indices = [1, 5, 2]
        self.compare(inchi, mult, u_indices)

    def testQuadriRadicalDoubleBondZwitterMult5(self):
        inchi = 'C8H14/c1-4-6-7-8(3)5-2/h5-6,8H,1-2,4,7H2,3H3'
        mult = 5
        u_indices = [1, 6, 2, 5]
        mol = self.compare(inchi, mult, u_indices)

    def testQuadri2DoubleBondMult5(self):
        inchi = 'C8H14/c1-5-7(3)8(4)6-2/h5-8H,1-2H2,3-4H3'
        mult = 5
        u_indices = [1, 5, 6, 2]
        self.compare(inchi, mult, u_indices)

    def testC2H3O3(self):
        adjlist = """
        1 C u0 p0 c0 {2,D} {3,S} {5,S}
        2 C u0 p0 c0 {1,D} {6,S} {7,S}
        3 O u0 p2 c0 {1,S} {4,S}
        4 O u0 p2 c0 {3,S} {8,S}
        5 O u1 p2 c0 {1,S}
        6 H u0 p0 c0 {2,S}
        7 H u0 p0 c0 {2,S}
        8 H u0 p0 c0 {4,S}
        """

        mol = Molecule().fromAdjacencyList(adjlist)

        inchi = 'C2H3O3/c1-2(3)5-4/h4H,1H2'
        mult = 2
        self.compare(inchi, mult)

    def testC5H6O(self):
        inchi = 'C5H6O/c6-5-3-1-2-4-5/h1-3,5H,4H2'
        mult = 3
        u_indices = [3, 6]
        self.compare(inchi, mult, u_indices)

    def testC5H6O_2(self):
        inchi = 'C5H6O/c1-5-3-2-4-6-5/h2-5H,1H2'
        mult = 3
        u_indices = [1,3]
        self.compare(inchi, mult, u_indices)

    def testC5H6O_3(self):
        inchi = 'C5H6O/c1-5-3-2-4-6-5/h2-5H,1H2'
        mult = 5
        u_indices = [1,2,3,4]
        self.compare(inchi, mult, u_indices)

    def testCO(self):
        inchi = 'CO/c1-2'
        mult = 1
        mol = self.compare(inchi, mult)

        assert mol.atoms[1].lonePairs == 1 # Oxygen

        assert mol.atoms[0].charge == -1
        assert mol.atoms[1].charge == +1

    def testMethylene(self):
        inchi = 'CH2/h1H2'

        mult = 1
        self.compare(inchi, mult)        

        mult = 3
        self.compare(inchi, mult)
    

    def testC4H6O(self):
        inchi = 'C4H6O/c1-2-3-4-5/h2H,3H2,1H3'
        mult = 3
        u_indices = [2,5]
        mol = self.compare(inchi, mult, u_indices)
        for at in mol.atoms:
            if at.isOxygen():
                assert at.lonePairs == 2
    
    def testC6H6(self):
        inchi = 'C6H6/c1-3-5-6-4-2/h1,6H,2,5H2'
        mult = 3
        u_indices = [3, 1]
        mol = self.compare(inchi, mult, u_indices)

    def testC4H6O_2(self):
        inchi = 'C4H6O/c1-2-3-4-5/h2,4H,1,3H2'
        mult = 3
        u_indices = [4, 5]
        mol = self.compare(inchi, mult, u_indices)

    def test_CO_triplet(self):
        from rmgpy.species import Species

        adjlist = """
        multiplicity 3
        1 C u2 p0 c0 {2,D}
        2 O u0 p2 c0 {1,D}

        """
        spc = Species(molecule=[Molecule().fromAdjacencyList(adjlist)])
        aug_inchi = spc.getAugmentedInChI()

        self.assertEqual(Species(molecule=[Molecule().fromAugmentedInChI(aug_inchi)]).isIsomorphic(spc), True)
        
    def test_CCCO_triplet(self):
        from rmgpy.species import Species

        adjlist = """
        multiplicity 3
1 C u0 p0 c0 {2,D} {5,S} {6,S}
2 C u0 p0 c0 {1,D} {3,S} {7,S}
3 C u1 p0 c0 {2,S} {4,S} {8,S}
4 O u1 p2 c0 {3,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {3,S}
        """
        mol = Molecule().fromAdjacencyList(adjlist)
        
        spc = Species(molecule=[mol])
        spc.generateResonanceIsomers()
        aug_inchi = spc.getAugmentedInChI()

        self.assertEqual(Species(molecule=[Molecule().fromAugmentedInChI(aug_inchi)]).isIsomorphic(spc), True)


    def testCreateULayer(self):
        from rmgpy.molecule.parser import createULayer

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

        mol1 = Molecule().fromAdjacencyList(adjlist1)

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

        mol2 = Molecule().fromAdjacencyList(adjlist2)
        self.assertEqual(createULayer(mol1), createULayer(mol2))

    def test_find_4_atom_3_bond_path(self):
        from rmgpy.molecule.parser import find_4_atom_3_bond_path

        mol = Molecule().fromSMILES("C=CC=C")#1,3-butadiene

        start, end = mol.atoms[0], mol.atoms[3]
        path = find_4_atom_3_bond_path(start, end)
        self.assertIsNotNone(path)

        mol = Molecule().fromSMILES("C=CC=O")#Acrolein

        start, end = mol.atoms[0], mol.atoms[3]
        path = find_4_atom_3_bond_path(start, end)
        self.assertIsNotNone(path)

        start, end = mol.atoms[0], mol.atoms[4]#wrong end
        path = find_4_atom_3_bond_path(start, end)
        self.assertIsNone(path)

        start, end = mol.atoms[-1], mol.atoms[3]#wrong start
        path = find_4_atom_3_bond_path(start, end)
        self.assertIsNone(path)

        mol = Molecule().fromSMILES("C=CC=CC=C")#1,3,5-hexatriene

        start, end = mol.atoms[0], mol.atoms[5]
        path = find_4_atom_3_bond_path(start, end)
        self.assertIsNotNone(path)        

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
        path = find_4_atom_3_bond_path(start, end)
        self.assertIsNotNone(path)     

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
        path = find_4_atom_3_bond_path(start, end)
        self.assertIsNone(path)   

        mol = Molecule().fromSMILES("C1=CC=CC=C1")#benzene

        start, end = mol.atoms[0], mol.atoms[5]
        path = find_4_atom_3_bond_path(start, end)
        self.assertIsNotNone(path)    

        mol = Molecule().fromSMILES("C=C=C=C")#C4H4

        start, end = mol.atoms[0], mol.atoms[3]
        path = find_4_atom_3_bond_path(start, end)
        self.assertIsNotNone(path)    


    def testC3H4(self):
        inchi = 'C3H4/c1-3-2/h1,3H,2H2/'
        mult = 3
        u_indices = [1, 1]
        mol = self.compare(inchi, mult, u_indices)

    def test_Buta13diyl_triplet(self):
        from rmgpy.species import Species
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
        mol = Molecule().fromAdjacencyList(adjlist)
        spc = Species(molecule=[mol])

        aug_inchi = spc.getAugmentedInChI()
        self.assertEqual(Species(molecule=[Molecule().fromAugmentedInChI(aug_inchi)]).isIsomorphic(spc), True)


    def test_C6H8O2(self):
        inchi = 'C6H8O2/c1-3-5(7)6(8)4-2/h3-6H,1-2H2'
        mult = 3
        u_indices = [7,8]
        self.compare(inchi, mult, u_indices)

    def test_C3H3O3(self):
        inchi = 'C3H3O3/c1-2-5-3-6-4/h1-3H'
        mult = 4
        u_indices = [1,3,4]
        self.compare(inchi, mult, u_indices)

    def test_CH2O2(self):
        inchi = 'CH2O2/c2-1-3/h1H,(H,2,3)'
        mult = 3
        u_indices = [1,3]
        self.compare(inchi, mult, u_indices)

    def test_C2H2O3(self):
        inchi = 'C2H2O3/c1-5-2(3)4/h1H2'
        mult = 3
        u_indices = [1,3]
        self.compare(inchi, mult, u_indices)

    def test_C3H4O4(self):
        inchi = 'C3H4O4/c4-3(5)1-2-7-6/h1-3,6H'
        mult = 3
        u_indices = [4,5]
        self.compare(inchi, mult, u_indices)

    def test_C6H6O4(self):
        inchi = 'InChI=1S/C6H6O4/c1-2-4-9-6(7)3-5-10-8/h2-3H,1,5H2'
        mult = 5
        u_indices = [1,3,4,8]
        self.compare(inchi, mult, u_indices)

if __name__ == '__main__':
    unittest.main()