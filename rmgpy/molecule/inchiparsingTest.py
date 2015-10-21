import unittest

from rmgpy.molecule import Molecule
from rmgpy.molecule.parser import fromAugmentedInChI
from rmgpy.species import Species

class InChITest(unittest.TestCase):

    def compare(self, inchi, mult, u_indices = None):
        from rmgpy.molecule.util import retrieveElementCount, VALENCES, ORDERS

        u_layer = ','.join([str(i) for i in u_indices]) if u_indices else None

        if u_layer:
            aug_inchi = 'InChI=1/' + inchi  + '/mult' + str(mult) + '/u' + u_layer
        else: 
            aug_inchi = 'InChI=1/' + inchi  + '/mult' + str(mult)

        mol = Molecule()
        mol = fromAugmentedInChI(mol, aug_inchi)
        self.assertEqual(mol.getNumberOfRadicalElectrons(), mult - 1)

        for at in mol.atoms:
            order = 0
            bonds = at.edges.values()
            for bond in bonds:
                order += ORDERS[bond.order]

            self.assertTrue((order + at.radicalElectrons + 2*at.lonePairs + at.charge) == VALENCES[at.symbol])

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
        u_indices = [2,4]
        mol = self.compare(inchi, mult, u_indices)
        for at in mol.atoms:
            if at.isOxygen():
                self.assertTrue(at.lonePairs == 2)
    
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

        adjlist = """
        multiplicity 3
        1 C u2 p0 c0 {2,D}
        2 O u0 p2 c0 {1,D}

        """
        spc = Species(molecule=[Molecule().fromAdjacencyList(adjlist)])
        aug_inchi = spc.getAugmentedInChI()

        self.assertEqual(Species(molecule=[Molecule().fromAugmentedInChI(aug_inchi)]).isIsomorphic(spc), True)
        
    def test_CCCO_triplet(self):

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
        from rmgpy.molecule.parser import createULayer, toRDKitMol
        from rdkit import Chem

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

        m1 = toRDKitMol(mol1)
        inchi1 , auxinfo1 = Chem.MolToInchiAndAuxInfo(m1, options='-SNon')

        mol2 = Molecule().fromAdjacencyList(adjlist2)
        m2 = toRDKitMol(mol2)
        inchi2 , auxinfo2 = Chem.MolToInchiAndAuxInfo(m2, options='-SNon')

        print '\n'.join([auxinfo1, auxinfo2])
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
        spc = Species(molecule=[Molecule().fromAugmentedInChI('InChI=1/'+inchi+'mult3/u1,1')])
        spc.generateResonanceIsomers()

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

    def test_C3H2O3(self):

        inchi = 'InChI=1S/C3H2O3/c1-2-3(4)6-5/h1H2'
        mult = 3
        u_indices = [2,5]

        aug_inchi = inchi+'/mult'+str(mult) + '/u2,5'
        spc = Species(molecule=[Molecule().fromAugmentedInChI(aug_inchi)])

        self.compare(inchi, mult, u_indices)

    def test_C6H6O6(self):
        inchi = 'C6H6O6/c7-6(2-5-12-9)10-3-1-4-11-8/h1,7H,4-5H2'
        mult = 5
        u_indices = [2,3,8,9]
        self.compare(inchi, mult, u_indices)

    def test_C3H2(self):
        inchi = 'C3H2/c1-3-2/h1-2H'
        mult = 3
        u_indices = [1,1]
        self.compare(inchi, mult, u_indices)

    def test_Isomorphic_Different_InChIs(self):
        from rmgpy.molecule.parser import createULayer, toRDKitMol
        from rdkit import Chem

        adjlist1 = """
multiplicity 3
1  C u1 p0 c0 {2,S} {9,S} {10,S}
2  C u1 p0 c0 {1,S} {3,S} {11,S}
3  C u0 p0 c0 {2,S} {4,S} {5,S} {12,S}
4  C u0 p0 c0 {3,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {6,S} {7,S} {16,S}
6  C u0 p0 c0 {5,S} {17,S} {18,S} {19,S}
7  C u0 p0 c0 {5,S} {8,D} {20,S}
8  C u0 p0 c0 {7,D} {21,S} {22,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
        """

        adjlist2 = """
multiplicity 3
1  C u0 p0 c0 {2,D} {9,S} {10,S}
2  C u0 p0 c0 {1,D} {3,S} {11,S}
3  C u0 p0 c0 {2,S} {4,S} {5,S} {12,S}
4  C u0 p0 c0 {3,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {6,S} {7,S} {16,S}
6  C u0 p0 c0 {5,S} {17,S} {18,S} {19,S}
7  C u1 p0 c0 {5,S} {8,S} {20,S}
8  C u1 p0 c0 {7,S} {21,S} {22,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
        """

        spc1 = Species(molecule=[Molecule().fromAdjacencyList(adjlist1)])
        spc2 = Species(molecule=[Molecule().fromAdjacencyList(adjlist2)])

        inchi1 = spc1.getAugmentedInChI()
        inchi2 = spc2.getAugmentedInChI()

        print inchi1, inchi2

        m1 = toRDKitMol(spc1.molecule[0])
        inchi1 , auxinfo1 = Chem.MolToInchiAndAuxInfo(m1, options='-SNon')
        
        m2 = toRDKitMol(spc2.molecule[0])
        inchi2 , auxinfo2 = Chem.MolToInchiAndAuxInfo(m2, options='-SNon')

        print '\n'.join([auxinfo1, auxinfo2])


    def test_group_adjacent_unpaired_electrons(self):
        from rmgpy.molecule.parser import parse_E_layer, toRDKitMol, group_adjacent_unpaired_electrons
        from rdkit import Chem


        adjlist = """

1 C 0 {4,D} 
2 C 0 {5,D}
3 C 1 {6,S}
4 C 0 {1,D} {7,S}
5 C 0 {2,D} {7,S}
6 C 1 {3,S} {7,S}
7 C 1 {4,S} {5,S} {6,S}
        """

        mol = Molecule().fromAdjacencyList(adjlist)
        u_layer = [7, 6, 3]

        m = toRDKitMol(mol)
        inchi , auxinfo = Chem.MolToInchiAndAuxInfo(m, options='-SNon')

        equivalent_atoms = parse_E_layer(auxinfo)
        pairs = group_adjacent_unpaired_electrons(mol, u_layer, equivalent_atoms)
        print pairs
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

        mol = Molecule().fromAdjacencyList(adjlist)
        u_layer = [7, 3, 6, 2]

        m = toRDKitMol(mol)
        inchi , auxinfo = Chem.MolToInchiAndAuxInfo(m, options='-SNon')

        equivalent_atoms = parse_E_layer(auxinfo)
        pairs = group_adjacent_unpaired_electrons(mol, u_layer, equivalent_atoms)
        print pairs

    def test_generate_combos(self):
        from rmgpy.molecule.parser import generate_combos
        
        group = [2, 6]
        equivalent_atoms = [[1,2],[3,4],[5,6], [7,8], [9,10]]
        expected =[[1,5], [2,6], [1,6], [2,5]]
        combos = generate_combos(group, equivalent_atoms)
        self.assertTrue(len(combos) == len(expected) and sorted(combos) == sorted(expected))

        group = [2, 6]
        equivalent_atoms = [[2, 4, 5, 6]]

        combos = generate_combos(group, equivalent_atoms)
        expected = [[2,4], [2,5], [2,6], [4,5], [4,6],[5,6]]        
        self.assertTrue(len(combos) == len(expected) and sorted(combos) == sorted(expected))
    
    def test_find_lowest_u_layer(self):
        from rmgpy.molecule.parser import parse_E_layer, toRDKitMol, group_adjacent_unpaired_electrons, find_lowest_u_layer
        from rdkit import Chem


        adjlist = """

1 C 0 {4,D} 
2 C 0 {5,D}
3 C 1 {6,S}
4 C 0 {1,D} {7,S}
5 C 0 {2,D} {7,S}
6 C 1 {3,S} {7,S}
7 C 1 {4,S} {5,S} {6,S}
        """

        mol = Molecule().fromAdjacencyList(adjlist)
        u_layer = [7, 6, 3]

        m = toRDKitMol(mol)
        inchi , auxinfo = Chem.MolToInchiAndAuxInfo(m, options='-SNon')
        equivalent_atoms = parse_E_layer(auxinfo)
        new_u_layer = find_lowest_u_layer(mol, u_layer, equivalent_atoms)
        self.assertEquals([1, 4, 7], new_u_layer)

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

        mol = Molecule().fromAdjacencyList(adjlist)
        u_layer = [7, 3, 6, 2]

        m = toRDKitMol(mol)
        inchi , auxinfo = Chem.MolToInchiAndAuxInfo(m, options='-SNon')
        equivalent_atoms = parse_E_layer(auxinfo)
        new_u_layer = find_lowest_u_layer(mol, u_layer, equivalent_atoms)
        self.assertEquals([1, 3, 5, 7], new_u_layer)

    def test_C3H4(self):
        inchi = 'InChI=1S/C3H4/c1-3-2/h1,3H,2H2'
        mult = 3
        u_indices = [1,1]
        self.compare(inchi, mult, u_indices)
    
    def test_C6H8(self):
        inchi = 'InChI=1S/C6H8/c1-3-5-6-4-2/h1,4H,2,5-6H2'
        mult = 5
        u_indices = [1,1,3,3]
        self.compare(inchi, mult, u_indices)

    def test_C6H10(self):
        inchi = 'InChI=1S/C6H10/c1-3-5-6-4-2/h3-4H,1-2,5-6H2'
        mult = 3
        u_indices = [1,3]
        self.compare(inchi, mult, u_indices)


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

        mol = Molecule().fromAdjacencyList(adjlist)

        aug_inchi = mol.toAugmentedInChI()
        self.assertTrue(aug_inchi.endswith('/u1,2,3,4'))
    

class ParseELayerTest(unittest.TestCase):
    def test_no_equivalence_layer(self):
        """Test that the absence of an E-layer results in an empty list."""
        from rmgpy.molecule.parser import parse_E_layer

        auxinfo = "AuxInfo=1/0/N:1/rA:1C/rB:/rC:;"
        e_layer = parse_E_layer(auxinfo)
        self.assertFalse(e_layer)

    def test_C8H22(self):
        from rmgpy.molecule.parser import parse_E_layer

        auxinfo = "AuxInfo=1/0/N:1,8,4,6,2,7,3,5/E:(1,2)(3,4)(5,6)(7,8)/rA:8C.2C.2CCCCCC/rB:s1;s2;s3;s3;s5;s5;d7;/rC:;;;;;;;;"
        e_layer = parse_E_layer(auxinfo)
        expected = [[1, 2], [3, 4], [5, 6], [7, 8]]
        self.assertTrue(len(e_layer) == len(expected) and sorted(e_layer) == sorted(expected))

    def test_C7H17(self):
        from rmgpy.molecule.parser import parse_E_layer

        auxinfo = "AuxInfo=1/0/N:3,5,7,2,4,6,1/E:(1,2,3)(4,5,6)/rA:7CCCCCCC/rB:s1;d2;s1;d4;s1;d6;/rC:;;;;;;;"
        e_layer = parse_E_layer(auxinfo)
        expected = [[1, 2, 3], [4, 5, 6]]
        self.assertTrue(len(e_layer) == len(expected) and sorted(e_layer) == sorted(expected))


class ParseNLayerTest(unittest.TestCase):
    def test_OCCC(self):
       from rmgpy.molecule.parser import parse_N_layer
       auxinfo = "AuxInfo=1/0/N:4,3,2,1/rA:4OCCC/rB:s1;s2;s3;/rC:;;;;"
       n_layer = parse_N_layer(auxinfo)
       expected = [4,3,2,1]
       self.assertTrue(len(n_layer) == len(expected) and sorted(n_layer) == sorted(expected))


if __name__ == '__main__':
    unittest.main()