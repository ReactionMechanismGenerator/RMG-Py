#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2020 Prof. William H. Green (whgreen@mit.edu),           #
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

from rmgpy.molecule import Molecule
from rmgpy.molecule.pathfinder import compute_atom_distance, find_adj_lone_pair_multiple_bond_delocalization_paths, \
    find_adj_lone_pair_radical_delocalization_paths, find_adj_lone_pair_radical_multiple_bond_delocalization_paths, \
    find_allyl_delocalization_paths, find_allyl_end_with_charge, find_butadiene, find_butadiene_end_with_charge, \
    find_lone_pair_multiple_bond_paths, find_N5dc_radical_delocalization_paths, find_shortest_path


class FindButadieneTest(unittest.TestCase):
    def test_13butadiene(self):
        mol = Molecule().from_smiles("C=CC=C")  # 1,3-butadiene

        start, end = mol.atoms[0], mol.atoms[3]
        path = find_butadiene(start, end)
        self.assertIsNotNone(path)

    def test_acrolein(self):
        mol = Molecule().from_smiles("C=CC=O")  # Acrolein

        start, end = mol.atoms[0], mol.atoms[3]
        path = find_butadiene(start, end)
        self.assertIsNotNone(path)

        start, end = mol.atoms[0], mol.atoms[4]  # wrong end
        path = find_butadiene(start, end)
        self.assertIsNone(path)

        start, end = mol.atoms[-1], mol.atoms[3]  # wrong start
        path = find_butadiene(start, end)
        self.assertIsNone(path)

    def test_135hexatriene(self):
        mol = Molecule().from_smiles("C=CC=CC=C")  # 1,3,5-hexatriene

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
        mol = Molecule().from_adjacency_list(adjlist)  # 1,3-cyclohexadiene

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

        mol = Molecule().from_adjacency_list(adjlist)  # 1,4-cyclohexadiene

        start, end = mol.atoms[0], mol.atoms[3]
        path = find_butadiene(start, end)
        self.assertIsNone(path)

    def test_benzene(self):
        mol = Molecule().from_smiles("C1=CC=CC=C1")  # benzene

        start, end = mol.atoms[0], mol.atoms[5]
        path = find_butadiene(start, end)
        self.assertIsNotNone(path)

    def test_c4h4(self):
        mol = Molecule().from_smiles("C=C=C=C")  # C4H4

        start, end = mol.atoms[0], mol.atoms[3]
        path = find_butadiene(start, end)
        self.assertIsNotNone(path)


class FindAllylEndWithChargeTest(unittest.TestCase):
    def test_c2h2o3(self):
        adjlist = """
1 C u0 p0 c0 {5,D} {6,S} {7,S}
2 C u0 p0 c0 {3,D} {4,S} {5,S}
3 O u0 p2 c0 {2,D}
4 O u0 p3 c-1 {2,S}
5 O u0 p1 c+1 {1,D} {2,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {1,S}
        """

        mol = Molecule().from_adjacency_list(adjlist)
        start = mol.atoms[2]
        paths = find_allyl_end_with_charge(start)
        idx_path = sorted([[mol.atoms.index(atom) + 1 for atom in path[0::2]] for path in paths])

        expected_idx_path = [[3, 2, 4], [3, 2, 5]]
        self.assertEquals(idx_path, expected_idx_path)

    def test_c3h2(self):
        inchi = "InChI=1S/C3H2/c1-3-2/h1-2H"
        mol = Molecule().from_inchi(inchi)
        start = mol.atoms[0]
        path = find_allyl_end_with_charge(start)[0]
        idx_path = [mol.atoms.index(atom) + 1 for atom in path[0::2]]

        expected_idx_path = [1, 3, 2]
        self.assertEquals(idx_path, expected_idx_path)

    def test_c3h4(self):
        inchi = "InChI=1S/C3H4/c1-3-2/h1,3H,2H2"
        mol = Molecule().from_inchi(inchi)
        start = mol.atoms[0]
        path = find_allyl_end_with_charge(start)[0]
        idx_path = [mol.atoms.index(atom) + 1 for atom in path[0::2]]

        expected_idx_path = [1, 3, 2]
        self.assertEquals(idx_path, expected_idx_path)

    def test_c3h2o3(self):
        adjlist = """
1 C u0 p0 c0 {2,D} {7,S} {8,S}
2 C u0 p0 c0 {1,D} {3,D}
3 C u0 p0 c0 {2,D} {4,S} {6,S}
4 O u0 p3 c-1 {3,S}
5 O u0 p2 c0 {6,D}
6 O u0 p1 c+1 {3,S} {5,D}
7 H u0 p0 c0 {1,S}
8 H u0 p0 c0 {1,S}
        """

        mol = Molecule().from_adjacency_list(adjlist)
        start = mol.atoms[1]
        paths = find_allyl_end_with_charge(start)
        idx_paths = sorted([[mol.atoms.index(atom) + 1 for atom in path[0::2]] for path in paths])
        idx_paths = sorted(idx_paths)

        expected_idx_paths = [[2, 3, 4], [2, 3, 6]]
        self.assertEquals(idx_paths, expected_idx_paths)

    def test_c3h4o4(self):
        inchi = "InChI=1S/C3H4O4/c4-3(5)1-2-7-6/h1-3,6H"
        mol = Molecule().from_inchi(inchi)
        start = mol.atoms[6]
        path = find_allyl_end_with_charge(start)[0]
        idx_path = [mol.atoms.index(atom) + 1 for atom in path[0::2]]

        expected_idx_path = [7, 2, 1]
        self.assertEquals(idx_path, expected_idx_path)

    def test_c5h6o(self):
        inchi = "InChI=1S/C5H6O/c6-5-3-1-2-4-5/h1-3,5H,4H2"
        mol = Molecule().from_inchi(inchi)
        start = mol.atoms[1]
        path = find_allyl_end_with_charge(start)[0]
        idx_path = [mol.atoms.index(atom) + 1 for atom in path[0::2]]

        expected_idx_path = [2, 1, 3]
        self.assertEquals(idx_path, expected_idx_path)


class FindButadieneEndWithChargeTest(unittest.TestCase):
    def test_co(self):
        adjlist = """
1 C u0 p1 c-1 {2,T}
2 O u0 p1 c+1 {1,T}
        """

        mol = Molecule().from_adjacency_list(adjlist)
        start = mol.atoms[0]
        path = find_butadiene_end_with_charge(start)
        idx_path = [mol.atoms.index(atom) + 1 for atom in path[0::2]]

        expected_idx_path = [1, 2]
        self.assertEquals(idx_path, expected_idx_path)

    def test_c2h2o3(self):
        adjlist = """
1 C u0 p0 c0 {5,D} {6,S} {7,S}
2 C u0 p0 c0 {3,D} {4,S} {5,S}
3 O u0 p2 c0 {2,D}
4 O u0 p3 c-1 {2,S}
5 O u0 p1 c+1 {1,D} {2,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {1,S}
        """

        mol = Molecule().from_adjacency_list(adjlist)
        start = mol.atoms[0]
        path = find_butadiene_end_with_charge(start)
        idx_path = [mol.atoms.index(atom) + 1 for atom in path[0::2]]

        expected_idx_path = [1, 5]
        self.assertEquals(idx_path, expected_idx_path)

    def test_c3h2o3(self):
        adjlist = """
1 C u0 p0 c0 {2,D} {7,S} {8,S}
2 C u0 p0 c0 {1,D} {3,D}
3 C u0 p0 c0 {2,D} {4,S} {6,S}
4 O u0 p3 c-1 {3,S}
5 O u0 p2 c0 {6,D}
6 O u0 p1 c+1 {3,S} {5,D}
7 H u0 p0 c0 {1,S}
8 H u0 p0 c0 {1,S}
        """

        mol = Molecule().from_adjacency_list(adjlist)
        start = mol.atoms[4]
        path = find_butadiene_end_with_charge(start)
        idx_path = [mol.atoms.index(atom) + 1 for atom in path[0::2]]

        expected_idx_path = [5, 6]
        self.assertEquals(idx_path, expected_idx_path)

    def test_c4h6o(self):
        adjlist = """
1  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
2  C u0 p1 c-1 {1,S} {3,S} {9,S}
3  C u0 p0 c0 {2,S} {4,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {5,T}
5  O u0 p1 c+1 {4,T}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
        """

        mol = Molecule().from_adjacency_list(adjlist)
        start = mol.atoms[3]
        path = find_butadiene_end_with_charge(start)
        idx_path = [mol.atoms.index(atom) + 1 for atom in path[0::2]]

        expected_idx_path = [4, 5]
        self.assertEquals(idx_path, expected_idx_path)

    def test_c5h6o2(self):
        adjlist = """
1  C u0 p1 c-1 {5,S} {7,S} {8,S}
2  C u0 p0 c0 {3,D} {4,S} {9,S}
3  C u0 p0 c0 {2,D} {5,S} {10,S}
4  C u0 p0 c0 {2,S} {6,D} {11,S}
5  C u0 p0 c0 {1,S} {3,S} {6,S} {12,S}
6  O u0 p1 c+1 {4,D} {5,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
        """

        mol = Molecule().from_adjacency_list(adjlist)
        start = mol.atoms[2]
        path = find_butadiene_end_with_charge(start)
        idx_path = [mol.atoms.index(atom) + 1 for atom in path[0::2]]

        expected_idx_path = [3, 2, 4, 6]
        self.assertEquals(idx_path, expected_idx_path)

    def test_c6h6o4(self):
        inchi = "InChI=1S/C6H6O4/c1-2-4-9-6(7)3-5-10-8/h2-3H,1,5H2"
        mol = Molecule().from_inchi(inchi)
        start = mol.atoms[0]
        path = find_butadiene_end_with_charge(start)
        idx_path = [mol.atoms.index(atom) + 1 for atom in path[0::2]]

        expected_idx_path = [1, 2, 4, 9]
        self.assertEquals(idx_path, expected_idx_path)

    def test_c6h6o6(self):
        inchi = "InChI=1S/C6H6O6/c7-6(2-5-12-9)10-3-1-4-11-8/h1,7H,4-5H2"
        mol = Molecule().from_inchi(inchi)
        start = mol.atoms[2]
        path = find_butadiene_end_with_charge(start)
        idx_path = [mol.atoms.index(atom) + 1 for atom in path[0::2]]

        expected_idx_path = [3, 10]
        self.assertEquals(idx_path, expected_idx_path)


class ShortestPathTest(unittest.TestCase):

    def test_ccc(self):
        smi = 'CCC'
        mol = Molecule().from_smiles(smi)
        start = mol.atoms[0]
        end = mol.atoms[2]

        path = find_shortest_path(start, end)
        self.assertEquals(len(path), 3)

    def test_cyclohexane(self):
        smi = 'C1CCCCC1'
        mol = Molecule().from_smiles(smi)
        start = mol.atoms[0]
        end = mol.atoms[2]

        path = find_shortest_path(start, end)
        self.assertEquals(len(path), 3)

    def test_bicyclo420octane(self):
        smi = 'C12CCC1CCCC2'
        mol = Molecule().from_smiles(smi)
        start = mol.atoms[0]
        end = mol.atoms[4]

        path = find_shortest_path(start, end)
        self.assertEquals(len(path), 3)


class DistanceComputingTest(unittest.TestCase):

    def test_2_atoms(self):
        smi = 'CCC'
        mol = Molecule().from_smiles(smi)
        atom_indices = [1, 2]
        distances = compute_atom_distance(atom_indices, mol)

        expected = {(1, 2): 1}
        self.assertEquals(distances, expected)

    def test_3_atoms(self):
        smi = 'CCC'
        mol = Molecule().from_smiles(smi)
        atom_indices = [1, 2, 3]
        distances = compute_atom_distance(atom_indices, mol)

        expected = {
            (1, 2): 1,
            (1, 3): 2,
            (2, 3): 1,
        }
        self.assertEquals(distances, expected)


class FindAllylDelocalizationPathsTest(unittest.TestCase):
    """
    test the find_allyl_delocalization_paths method
    """

    def test_allyl_radical(self):
        smiles = "[CH2]C=C"
        mol = Molecule().from_smiles(smiles)
        paths = find_allyl_delocalization_paths(mol.atoms[0])
        self.assertTrue(paths)

    def test_nitrogenated_birad(self):
        smiles = '[N]C=[CH]'
        mol = Molecule().from_smiles(smiles)
        paths = find_allyl_delocalization_paths(mol.atoms[0])
        self.assertTrue(paths)


class FindLonePairMultipleBondPathsTest(unittest.TestCase):
    """
    test the find_lone_pair_multiple_bond_paths method
    """

    def test_azide(self):
        smiles = "[N-]=[N+]=N"
        mol = Molecule().from_smiles(smiles)
        paths = find_lone_pair_multiple_bond_paths(mol.atoms[2])
        self.assertTrue(paths)

    def test_nh2cho(self):
        smiles = 'NC=O'
        mol = Molecule().from_smiles(smiles)
        paths = find_lone_pair_multiple_bond_paths(mol.atoms[0])
        self.assertTrue(paths)

    def test_n2oa(self):
        smiles = "[N-]=[N+]=O"
        mol = Molecule().from_smiles(smiles)
        paths = find_lone_pair_multiple_bond_paths(mol.atoms[0])
        self.assertTrue(paths)

    def test_n2ob(self):
        smiles = "N#[N+][O-]"
        mol = Molecule().from_smiles(smiles)
        paths = find_lone_pair_multiple_bond_paths(mol.atoms[2])
        self.assertTrue(paths)

    def test_hn3(self):
        smiles = "[NH-][N+]#N"
        mol = Molecule().from_smiles(smiles)
        paths = find_lone_pair_multiple_bond_paths(mol.atoms[0])
        self.assertTrue(paths)

    def test_sn2(self):
        smiles = "OS(O)=[N+]=[N-]"
        mol = Molecule().from_smiles(smiles)
        paths = find_lone_pair_multiple_bond_paths(mol.atoms[2])
        self.assertTrue(paths)

    def test_h2nnoo(self):
        smiles = "N[N+]([O-])=O"
        mol = Molecule().from_smiles(smiles)
        paths = find_lone_pair_multiple_bond_paths(mol.atoms[0])
        self.assertTrue(paths)


class FindAdjLonePairRadicalDelocalizationPaths(unittest.TestCase):
    """
    test the find_lone_pair_radical_delocalization_paths method
    """

    def test_no2a(self):
        smiles = "[O]N=O"
        mol = Molecule().from_smiles(smiles)
        paths = find_adj_lone_pair_radical_delocalization_paths(mol.atoms[0])
        self.assertTrue(paths)

    def test_no2b(self):
        smiles = "[O-][N+]=O"
        mol = Molecule().from_smiles(smiles)
        paths = find_adj_lone_pair_radical_delocalization_paths(mol.atoms[1])
        self.assertTrue(paths)

    def test_hoso(self):
        smiles = "[O]SO"
        mol = Molecule().from_smiles(smiles)
        paths = find_adj_lone_pair_radical_delocalization_paths(mol.atoms[0])
        self.assertTrue(paths)

    def test_double_bond(self):
        adj = """multiplicity 2
                 1 O u1 p1 c+1 {2,D}
                 2 N u0 p2 c-1 {1,D}"""
        mol = Molecule().from_adjacency_list(adj)
        paths = find_adj_lone_pair_radical_delocalization_paths(mol.atoms[0])
        self.assertTrue(paths)


class FindAdjLonePairMultipleBondDelocalizationPaths(unittest.TestCase):
    """
    test the find_lone_pair_multiple_bond_delocalization_paths method
    """

    def test_sho3(self):
        smiles = "O=[SH](=O)[O]"
        mol = Molecule().from_smiles(smiles)
        paths = find_adj_lone_pair_multiple_bond_delocalization_paths(mol.atoms[0])
        self.assertTrue(paths)


class FindAdjLonePairRadicalMultipleBondDelocalizationPaths(unittest.TestCase):
    """
    test the find_lone_pair_radical_multiple_bond_delocalization_paths method
    """

    def test_ns(self):
        smiles = "N#[S]"
        mol = Molecule().from_smiles(smiles)
        paths = find_adj_lone_pair_radical_multiple_bond_delocalization_paths(mol.atoms[1])
        self.assertTrue(paths)

    def test_hso3(self):
        smiles = "O[S](=O)=O"
        mol = Molecule().from_smiles(smiles)
        paths = find_adj_lone_pair_radical_multiple_bond_delocalization_paths(mol.atoms[1])
        self.assertTrue(paths)


class FindN5dcRadicalDelocalizationPaths(unittest.TestCase):
    """
    test the find_N5dc_radical_delocalization_paths method
    """

    def test_hnnoo(self):
        smiles = "N=[N+]([O])([O-])"
        mol = Molecule().from_smiles(smiles)
        paths = find_N5dc_radical_delocalization_paths(mol.atoms[1])
        self.assertTrue(paths)
