#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2023 Prof. William H. Green (whgreen@mit.edu),           #
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


from rmgpy.molecule.molecule import Molecule
from rmgpy.molecule.resonance import (
    _clar_optimization,
    _clar_transformation,
    generate_clar_structures,
    generate_kekule_structure,
    generate_optimal_aromatic_resonance_structures,
    generate_resonance_structures,
)


class ResonanceTest:
    def test_allyl_shift(self):
        """Test allyl shift for hexadienyl radical"""
        mol_list = generate_resonance_structures(Molecule(smiles="C=C[CH]C=CC"))
        assert len(mol_list) == 3

    def test_trirad_allyl_shift(self):
        """Test allyl shift for a tri-rad carbon"""
        mol_list = generate_resonance_structures(Molecule(smiles="[C]N=N"))
        assert len(mol_list) == 3
        assert any([any([atom.charge != 0 for atom in mol.vertices]) for mol in mol_list])  # expecting [C]=[N+.][NH-]

    def test_oxime(self):
        """Test resonance structure generation for CC=N[O] radical

        Simple case for lone pair <=> radical resonance"""
        mol_list = generate_resonance_structures(Molecule(smiles="CC=N[O]"))
        assert len(mol_list) == 3
        assert any([any([atom.charge != 0 for atom in mol.vertices]) for mol in mol_list])

    def test_ring_allyl_shift(self):
        """Test allyl shift for a cyclic species with heteroatoms"""
        mol_list = generate_resonance_structures(Molecule(smiles="[CH]1C=NC=N1"))
        assert len(mol_list) == 5

    def test_carbene_allyl_shift(self):
        """Test allyl shift for a carbene species"""
        mol_list = generate_resonance_structures(Molecule(smiles="[C]=C=O"))
        assert len(mol_list) == 2

    def test_ch2chcho(self):
        """Test resonance structure generation for C=C[CH][O] bi-radical

        Test case for allyl bi-radical resonance"""
        mol_list = generate_resonance_structures(Molecule(smiles="C=C[CH][O]"))
        assert len(mol_list) == 3

    def ch2no(self):
        """Test combined resonance transitions of allyl-shift and lonePair-radical"""
        mol_list = generate_resonance_structures(Molecule(smiles="[CH2]N=O"))
        assert len(mol_list) == 3

    def ch3s2o2(self):
        """Test combined resonance transitions of one_pair_radical_multiple_bond"""
        mol_list = generate_resonance_structures(Molecule(smiles="CSS(=O)[O]"))
        assert len(mol_list) == 3

    def n2so2(self):
        """Test the resonance transitions of a species with several hereroatoms and several multiple bonds"""
        mol_list = generate_resonance_structures(Molecule(smiles="[N-]=[N+]=S(=O)=O"))
        assert len(mol_list) == 2

    def nsh(self):
        """Test that a resonance structure with a minimal octet deviation but higher charge span is filtered out"""
        mol_list = generate_resonance_structures(Molecule(smiles="N#S"))
        assert len(mol_list) == 1
        assert all([atom.charge == 0 for atom in mol_list[0].vertices])

    def test_nco(self):
        """Test resonance structure generation for NCO

        NCO should only have two resonance structures [N.]=C=O <=> N#C[O.], and not a third structure which has
        the same octet deviation, has a charge separation, but no ne radical site: [N+.]#C[O-]
        """
        mol_list = generate_resonance_structures(Molecule(smiles="[N]=C=O"))
        assert len(mol_list) == 2
        assert all([all([atom.charge == 0 for atom in mol.vertices]) for mol in mol_list])  # none of the
        # structures should be charged

    def test_no2(self):
        """Test resonance structure generation for [O]N=O radical

        Test case for the lone pair <=> radical resonance transition.
        Also tests that the filtering function allows charge separation when the radical site is changed.
        """
        mol_list = generate_resonance_structures(Molecule(smiles="[O]N=O"))
        assert len(mol_list) == 2
        assert any([any([atom.charge != 0 for atom in mol.vertices]) for mol in mol_list])  # one of the
        # structures should be charged

    def test_n2o(self):
        """Test resonance structure generation for N#[N+][O-]

        A classic N5ddc <=> N5tc resonance transition"""
        mol_list = generate_resonance_structures(Molecule(smiles="N#[N+][O-]"))
        assert len(mol_list) == 2
        assert all([any([atom.charge != 0 for atom in mol.vertices]) for mol in mol_list])  # both structures
        # should have some charged atoms

        sbonds = 0
        dbonds = 0
        tbonds = 0
        for mol in mol_list:
            for atom in mol.atoms:
                for bond in atom.bonds.values():
                    if bond.is_single():
                        sbonds += 1
                    elif bond.is_double():
                        dbonds += 1
                    elif bond.is_triple():
                        tbonds += 1
        assert sbonds / 2 == 1  # each bond is counted twice above
        assert dbonds / 2 == 2
        assert tbonds / 2 == 1

    def test_azide(self):
        """Test resonance structure generation for ethyl azide

        Simple case for N5ddc <=> N5tc resonance
        Azides are described by three resonance structures: N=[N+]=[N-] <=> [NH-][N+]#N <=> [NH+]#[N+][N-2]
        """
        mol_list = generate_resonance_structures(Molecule(smiles="CCN=[N+]=[N-]"))
        assert len(mol_list) == 2
        assert all([any([atom.charge != 0 for atom in mol.vertices]) for mol in mol_list])

    def test_ozone(self):
        """Test resonance structure generation for O3, S3 and SO2.

        Compare that these iso-electronic structures have the same number of resonance structures
        """
        mol_list_1 = generate_resonance_structures(Molecule(smiles="[O-][O+]=O"))
        assert len(mol_list_1) == 1
        mol_list_2 = generate_resonance_structures(Molecule(smiles="O=S=O"))
        assert len(mol_list_2) == 1
        mol_list_3 = generate_resonance_structures(Molecule(smiles="S=S=S"))
        assert len(mol_list_3) == 1

    def test_hco_vs_hcs(self):
        """Test resonance structure generation for [CH]=O and [CH]=S

        These iso-electronic structures have a different(!) number of resonance structures
        """
        mol_list_1 = generate_resonance_structures(Molecule(smiles="[CH]=O"))
        assert len(mol_list_1) == 1
        mol_list_2 = generate_resonance_structures(Molecule(smiles="[CH]=S"))
        assert len(mol_list_2) == 2

    def test_no(self):
        """Test that an incorrect NO structure [::N][::O.] is correctly identified as [:N.]=[::O]

        The incorrect structure could be generated from HON (O[::N]) during an RMG run, and should be identified as NO.
        The original structure should be kept as unreactive (appended at the end of the molecule list)
        """
        mol_list = generate_resonance_structures(
            Molecule().from_adjacency_list(
                """multiplicity 2
                                                                                 1 N u0 p2 c0 {2,S}
                                                                                 2 O u1 p2 c0 {1,S}"""
            )
        )
        assert len(mol_list) == 2
        assert mol_list[0].reactive
        assert not mol_list[1].reactive
        assert mol_list[0].vertices[0].lone_pairs + mol_list[0].vertices[1].lone_pairs == 3
        assert mol_list[1].vertices[0].lone_pairs + mol_list[1].vertices[1].lone_pairs == 4

    def test_n5dc_radical(self):
        """Test the N5dc radical resonance transformation

        We should see N=[N+]([O])([O-]) <=> [NH-][N+]([O])=O
        Two isomorphic structures should be included in mol_list: N=[N+]([O])([O-]) <=> N=[N+]([O-])([O])
        """
        mol = Molecule(smiles="N=[N+]([O-])[O]")
        mol_list = generate_resonance_structures(mol, keep_isomorphic=True)
        assert len(mol_list) == 6
        isomorphic_counter = 0
        negatively_charged_nitrogen = 0
        for mol1 in mol_list:
            if mol1.is_isomorphic(mol):
                isomorphic_counter += 1
            for atom in mol1.vertices:
                if atom.is_nitrogen() and atom.charge < 0:
                    negatively_charged_nitrogen += 1
        assert isomorphic_counter == 2
        assert negatively_charged_nitrogen == 2

    def test_n5dc(self):
        """Test the N5dc resonance transformation

        We should see N[N+]([O-])=O <=> N[N+](=O)[O-], which are isomorphic"""
        mol = Molecule(smiles="N[N+]([O-])=O")
        mol_list = generate_resonance_structures(mol, keep_isomorphic=True)
        assert len(mol_list) == 2
        assert mol_list[0].is_isomorphic(mol_list[1])

    def test_styryl1(self):
        """Test resonance structure generation for styryl, with radical on branch

        In this case, the radical can be delocalized into the aromatic ring"""
        mol_list = generate_resonance_structures(Molecule(smiles="c1ccccc1[C]=C"))
        assert len(mol_list) == 4

    def test_styryl2(self):
        """Test resonance structure generation for styryl, with radical on ring

        In this case, the radical can be delocalized into the aromatic ring"""
        mol_list = generate_resonance_structures(Molecule(smiles="C=C=C1C=C[CH]C=C1"))
        assert len(mol_list) == 3

    def test_naphthyl(self):
        """Test resonance structure generation for naphthyl radical

        In this case, the radical is orthogonal to the pi-orbital plane and cannot delocalize
        """
        mol_list = generate_resonance_structures(Molecule(smiles="c12[c]cccc1cccc2"))
        assert len(mol_list) == 4

    def test_methyl_napthalene(self):
        """Test resonance structure generation for methyl naphthalene

        Example of stable polycyclic aromatic species"""
        mol_list = generate_resonance_structures(Molecule(smiles="CC1=CC=CC2=CC=CC=C12"))
        assert len(mol_list) == 4

    def test_methyl_phenanthrene(self):
        """Test resonance structure generation for methyl phenanthrene

        Example of stable polycyclic aromatic species"""
        mol_list = generate_resonance_structures(Molecule(smiles="CC1=CC=CC2C3=CC=CC=C3C=CC=21"))
        assert len(mol_list) == 3

    def test_methyl_phenanthrene_radical(self):
        """Test resonance structure generation for methyl phenanthrene radical

        Example radical polycyclic aromatic species where the radical can delocalize"""
        mol_list = generate_resonance_structures(Molecule(smiles="[CH2]C1=CC=CC2C3=CC=CC=C3C=CC=21"))
        assert len(mol_list) == 9

    def test_aromatic_with_lone_pair_resonance(self):
        """Test resonance structure generation for aromatic species with lone pair <=> radical resonance"""
        mol_list = generate_resonance_structures(Molecule(smiles="c1ccccc1CC=N[O]"))
        assert len(mol_list) == 4

    def test_aromatic_with_n_resonance(self):
        """Test resonance structure generation for aromatic species with lone pair resonance"""
        mol_list = generate_resonance_structures(Molecule(smiles="c1ccccc1CCN=[N+]=[N-]"))
        assert len(mol_list) == 2

    def test_aromatic_with_o_n_resonance(self):
        """Test resonance structure generation for aromatic species with heteroatoms

        This test was specifically designed to recreate RMG-Py issue #1598.
        Key conditions: having heteroatoms, starting with aromatic structure, and keep_isomorphic=True
        """
        mol_list = generate_resonance_structures(
            Molecule().from_adjacency_list(
                """
multiplicity 2
1 O u0 p2 c0 {4,S} {16,S}
2 O u1 p2 c0 {4,S}
3 O u0 p3 c-1 {4,S}
4 N u0 p0 c+1 {1,S} {2,S} {3,S} {5,S}
5 C u0 p0 c0 {4,S} {6,B} {7,B}
6 C u0 p0 c0 {5,B} {8,B} {11,S}
7 C u0 p0 c0 {5,B} {10,B} {15,S}
8 C u0 p0 c0 {6,B} {9,B} {12,S}
9 C u0 p0 c0 {8,B} {10,B} {13,S}
10 C u0 p0 c0 {7,B} {9,B} {14,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {1,S}
"""
            ),
            keep_isomorphic=True,
        )
        assert len(mol_list) == 1

    def test_no_clar_structures(self):
        """Test that we can turn off Clar structure generation."""
        mol_list = generate_resonance_structures(Molecule(smiles="C1=CC=CC2C3=CC=CC=C3C=CC=21"), clar_structures=False)
        assert len(mol_list) == 2

    def test_c13h11_rad(self):
        """Test resonance structure generation for p-methylbenzylbenzene radical

        Has multiple resonance structures that break aromaticity of a ring"""
        mol_list = generate_resonance_structures(Molecule(smiles="[CH](c1ccccc1)c1ccc(C)cc1"))
        assert len(mol_list) == 6

    def test_c8h8(self):
        """Test resonance structure generation for 5,6-dimethylene-1,3-cyclohexadiene

        Example of molecule that RDKit considers aromatic, but RMG does not"""
        mol_list = generate_resonance_structures(Molecule(smiles="C=C1C=CC=CC1=C"))
        assert len(mol_list) == 1

    def test_c8h7_j(self):
        """Test resonance structure generation for 5,6-dimethylene-1,3-cyclohexadiene radical

        Example of molecule that RDKit considers aromatic, but RMG does not"""
        mol_list = generate_resonance_structures(Molecule(smiles="C=C1C=CC=CC1=[CH]"))
        assert len(mol_list) == 1

    def test_c8h7_j2(self):
        """Test resonance structure generation for 5,6-dimethylene-1,3-cyclohexadiene radical

        Example of molecule that RDKit considers aromatic, but RMG does not"""
        mol_list = generate_resonance_structures(Molecule(smiles="C=C1C=[C]C=CC1=C"))
        assert len(mol_list) == 1

    def test_c9h9_aro(self):
        """Test cyclopropyl benzene radical, aromatic SMILES"""
        mol = Molecule(smiles="[CH]1CC1c1ccccc1")
        mol_list = generate_resonance_structures(mol)
        assert len(mol_list) == 2

    def test_c9h9_kek(self):
        """Test cyclopropyl benzene radical, kekulized SMILES"""
        mol = Molecule(smiles="[CH]1CC1C1C=CC=CC=1")
        mol_list = generate_resonance_structures(mol)
        assert len(mol_list) == 2

    def test_benzene_aro(self):
        """Test benzene, aromatic SMILES"""
        mol = Molecule(smiles="c1ccccc1")
        mol_list = generate_resonance_structures(mol)
        assert len(mol_list) == 2

    def test_benzene_kek(self):
        """Test benzene, kekulized SMILES"""
        mol = Molecule(smiles="C1C=CC=CC=1")
        mol_list = generate_resonance_structures(mol)
        assert len(mol_list) == 2

    def test_c9h11_aro(self):
        """Test propylbenzene radical, aromatic SMILES"""
        mol = Molecule(smiles="[CH2]CCc1ccccc1")
        mol_list = generate_resonance_structures(mol)
        assert len(mol_list) == 2

    def test_c10h11_aro(self):
        """Test cyclobutylbenzene radical, aromatic SMILES"""
        mol = Molecule(smiles="[CH]1CCC1c1ccccc1")
        mol_list = generate_resonance_structures(mol)
        assert len(mol_list) == 2

    def test_c9h10_aro(self):
        """Test cyclopropylbenzene, aromatic SMILES"""
        mol = Molecule(smiles="C1CC1c1ccccc1")
        mol_list = generate_resonance_structures(mol)
        assert len(mol_list) == 2

    def test_c10h12_aro(self):
        """Test cyclopropylmethyl benzene, aromatic SMILES"""
        mol = Molecule(smiles="C1CC1c1c(C)cccc1")
        mol_list = generate_resonance_structures(mol)
        assert len(mol_list) == 2

    def test_c9h10_aro_2(self):
        """Test cyclopropyl benzene, generate aromatic resonance isomers"""
        mol = Molecule(smiles="C1CC1c1ccccc1")
        mol_list = generate_optimal_aromatic_resonance_structures(mol)
        assert len(mol_list) == 1

    def test_aryne_1_ring(self):
        """Test aryne resonance for benzyne"""
        mol1 = Molecule(smiles="C1=CC=C=C=C1")
        mol2 = Molecule(smiles="C1C#CC=CC=1")

        mol_list1 = generate_resonance_structures(mol1)
        assert len(mol_list1) == 2

        mol_list2 = generate_resonance_structures(mol2)
        assert len(mol_list2) == 2

        assert mol_list1[1].is_isomorphic(mol2)
        assert mol_list2[1].is_isomorphic(mol1)

    def test_aryne_2_rings(self):
        """Test aryne resonance in naphthyne"""
        mol1 = Molecule(smiles="C12=CC=C=C=C1C=CC=C2")
        mol2 = Molecule(smiles="C12C#CC=CC=1C=CC=C2")

        mol_list1 = generate_resonance_structures(mol1)
        assert len(mol_list1) == 2
        assert sum(1 for mol in mol_list1 if mol.reactive) == 2

        mol_list2 = generate_resonance_structures(mol2)
        assert len(mol_list2) == 3
        assert sum(1 for mol in mol_list2 if mol.reactive) == 2

        # Check that they both have an aromatic resonance form
        assert mol_list1[1].is_isomorphic(mol_list2[0])

    def test_aryne_3_rings(self):
        """Test aryne resonance in phenanthryne"""
        mol = Molecule(smiles="C12C#CC=CC=1C=CC3=C2C=CC=C3")

        mol_list = generate_resonance_structures(mol)
        assert len(mol_list) == 5
        assert sum(1 for mol in mol_list if mol.reactive) == 4

    def test_fused_aromatic1(self):
        """Test we can make aromatic perylene from both adjlist and SMILES"""
        perylene = Molecule().from_adjacency_list(
            """
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
"""
        )
        perylene2 = Molecule().from_smiles("c1cc2cccc3c4cccc5cccc(c(c1)c23)c54")
        for isomer in generate_optimal_aromatic_resonance_structures(perylene2):
            if perylene.is_isomorphic(isomer):
                break
        else:  # didn't break
            assert False, "{} isn't isomorphic with any aromatic forms of {}".format(perylene.to_smiles(), perylene2.to_smiles())

    def test_fused_aromatic2(self):
        """Test we can make aromatic naphthalene from both adjlist and SMILES"""
        naphthalene = Molecule().from_adjacency_list(
            """
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
"""
        )
        naphthalene2 = Molecule().from_smiles("C1=CC=C2C=CC=CC2=C1")
        for isomer in generate_optimal_aromatic_resonance_structures(naphthalene2):
            if naphthalene.is_isomorphic(isomer):
                break
        else:  # didn't break
            assert False, "{} isn't isomorphic with any aromatic forms of {}".format(naphthalene.to_smiles(), naphthalene2.to_smiles())

    def test_aromatic_resonance_structures(self):
        """Test that generate_optimal_aromatic_resonance_structures gives consistent output

        Check that we get the same resonance structure regardless of which structure we start with
        """
        # Kekulized form, radical on methyl
        struct1 = Molecule().from_adjacency_list(
            """
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
"""
        )
        # Kekulized form, radical on ring
        struct2 = Molecule().from_adjacency_list(
            """
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
"""
        )
        # Aromatic form
        struct3 = Molecule().from_adjacency_list(
            """
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
"""
        )
        result1 = generate_optimal_aromatic_resonance_structures(struct1)
        result2 = generate_optimal_aromatic_resonance_structures(struct2)
        result3 = generate_optimal_aromatic_resonance_structures(struct3)

        assert len(result1) == 1
        assert len(result2) == 1
        assert len(result3) == 1

        assert result1[0].is_isomorphic(result2[0])
        assert result1[0].is_isomorphic(result3[0])

    def test_bridged_aromatic(self):
        """Test that we can handle bridged aromatics.

        This is affected by how we perceive rings. Using get_smallest_set_of_smallest_rings gives
        non-deterministic output, so using get_all_cycles_of_size allows this test to pass.

        Update: Highly-strained fused rings are no longer considered aromatic."""
        mol = Molecule(smiles="c12c3cccc1c3ccc2")
        arom = Molecule().from_adjacency_list(
            """
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
"""
        )

        out = generate_resonance_structures(mol)

        assert len(out) == 1
        assert not arom.is_isomorphic(out[0])

    def test_polycyclic_aromatic_with_non_aromatic_ring(self):
        """Test that we can make aromatic resonance structures when there is a pseudo-aromatic ring.

        This applies in cases where RDKit misidentifies one ring as aromatic, but there are other
        rings in the molecule that are actually aromatic.

        Update: Highly-strained fused rings are no longer considered aromatic."""
        mol = Molecule(smiles="c1c2cccc1C(=C)C=[C]2")
        arom = Molecule().from_adjacency_list(
            """
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
"""
        )

        out = generate_resonance_structures(mol)

        assert len(out) == 5
        assert not any(arom.is_isomorphic(res) for res in out)

    def test_polycyclic_aromatic_with_non_aromatic_ring2(self):
        """Test that we can make aromatic resonance structures when there is a pseudo-aromatic ring.

        This applies in cases where RDKit misidentifies one ring as aromatic, but there are other
        rings in the molecule that are actually aromatic."""
        mol = Molecule(smiles="C=C(C1=CC2=C(C=C1C=C3)C4=CC5=CC=CC=C5C=C4C=C2)C3=C")
        arom = Molecule().from_adjacency_list(
            """
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
"""
        )

        out = generate_resonance_structures(mol)

        assert len(out) == 4
        assert arom.is_isomorphic(out[0])

    def test_kekulize_benzene(self):
        """Test that we can kekulize benzene."""
        arom = Molecule().from_adjacency_list(
            """
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
"""
        )
        keku = Molecule().from_adjacency_list(
            """
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
"""
        )
        out = generate_kekule_structure(arom)

        assert len(out) == 1
        assert out[0].is_isomorphic(keku)

    def test_kekulize_naphthalene(self):
        """Test that we can kekulize naphthalene."""
        arom = Molecule().from_adjacency_list(
            """
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
"""
        )
        out = generate_kekule_structure(arom)

        assert len(out) == 1
        assert not out[0].is_aromatic()

        bonds = set()
        for atom in out[0].atoms:
            for bond in atom.bonds.values():
                bonds.add(bond)
        d_bonds = 0
        for bond in bonds:
            if bond.is_double():
                d_bonds += 1

        assert d_bonds == 5

    def test_kekulize_phenanthrene(self):
        """Test that we can kekulize phenanthrene."""
        arom = Molecule().from_adjacency_list(
            """
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
"""
        )
        out = generate_kekule_structure(arom)

        assert len(out) == 1
        assert not out[0].is_aromatic()

        bonds = set()
        for atom in out[0].atoms:
            for bond in atom.bonds.values():
                bonds.add(bond)
        d_bonds = 0
        for bond in bonds:
            if bond.is_double():
                d_bonds += 1

        assert d_bonds == 7

    def test_kekulize_pyrene(self):
        """Test that we can kekulize pyrene."""
        arom = Molecule().from_adjacency_list(
            """
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
"""
        )
        out = generate_kekule_structure(arom)

        assert len(out) == 1
        assert not out[0].is_aromatic()

        bonds = set()
        for atom in out[0].atoms:
            for bond in atom.bonds.values():
                bonds.add(bond)
        d_bonds = 0
        for bond in bonds:
            if bond.is_double():
                d_bonds += 1

        assert d_bonds == 8

    def test_kekulize_corannulene(self):
        """Test that we can kekulize corannulene."""
        arom = Molecule().from_adjacency_list(
            """
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
"""
        )
        out = generate_kekule_structure(arom)

        assert len(out) == 1
        assert not out[0].is_aromatic()

        bonds = set()
        for atom in out[0].atoms:
            for bond in atom.bonds.values():
                bonds.add(bond)
        d_bonds = 0
        for bond in bonds:
            if bond.is_double():
                d_bonds += 1

        assert d_bonds == 10

    def test_kekulize_coronene(self):
        """Test that we can kekulize coronene."""
        arom = Molecule().from_adjacency_list(
            """
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
"""
        )
        out = generate_kekule_structure(arom)

        assert len(out) == 1
        assert not out[0].is_aromatic()

        bonds = set()
        for atom in out[0].atoms:
            for bond in atom.bonds.values():
                bonds.add(bond)
        d_bonds = 0
        for bond in bonds:
            if bond.is_double():
                d_bonds += 1

        assert d_bonds == 12

    def test_kekulize_bridged_aromatic(self):
        """Test that we can kekulize a bridged polycyclic aromatic species."""
        arom = Molecule().from_adjacency_list(
            """
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
"""
        )
        out = generate_kekule_structure(arom)

        assert len(out) == 1
        assert not out[0].is_aromatic()

        bonds = set()
        for atom in out[0].atoms:
            for bond in atom.bonds.values():
                bonds.add(bond)
        d_bonds = 0
        for bond in bonds:
            if bond.is_double():
                d_bonds += 1

        assert d_bonds == 5

    def test_multiple_kekulized_resonance_isomers_rad(self):
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
        molecule = Molecule().from_adjacency_list(adjlist_aromatic)
        assert molecule.is_aromatic(), "Starting molecule should be aromatic"
        mol_list = generate_resonance_structures(molecule)
        assert len(mol_list) == 4, "Expected 4 resonance structures, but generated {0}.".format(len(mol_list))
        aromatic = 0
        for mol in mol_list:
            if mol.is_aromatic():
                aromatic += 1
        assert aromatic == 1, "Should only have 1 aromatic resonance structure"

    def test_keep_isomorphic_structures_functions_when_true(self):
        """Test that keep_isomorphic works for resonance structure generation when True."""
        mol = Molecule(smiles="C=C[CH2]")
        mol.assign_atom_ids()
        out = generate_resonance_structures(mol, keep_isomorphic=True)

        assert len(out) == 2
        assert out[0].is_isomorphic(out[1])
        assert not out[0].is_identical(out[1])

    def test_keep_isomorphic_structures_functions_when_false(self):
        """Test that keep_isomorphic works for resonance structure generation when False."""
        mol = Molecule(smiles="C=C[CH2]")
        mol.assign_atom_ids()
        out = generate_resonance_structures(mol, keep_isomorphic=False)

        assert len(out) == 1

    def test_false_negative_aromaticity_perception(self):
        """Test that we obtain the correct aromatic structure for a monocyclic aromatic that RDKit mis-identifies."""
        mol = Molecule(smiles="[CH2]C=C1C=CC(=C)C=C1")
        out = generate_resonance_structures(mol)

        aromatic = Molecule().from_adjacency_list(
            """
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
"""
        )

        assert len(out) == 4
        assert any([m.is_isomorphic(aromatic) for m in out])

    def test_false_negative_polycyclic_aromaticity_perception(self):
        """Test that we generate proper structures for a polycyclic aromatic that RDKit mis-identifies."""
        mol = Molecule(smiles="C=C1C=CC=C2C=C[CH]C=C12")
        out = generate_resonance_structures(mol)

        clar = Molecule().from_adjacency_list(
            """
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
"""
        )

        assert len(out) == 6
        assert any([m.is_isomorphic(clar) for m in out])

    def test_false_negative_polycylic_aromaticity_perception2(self):
        """Test that we obtain the correct aromatic structure for a polycylic aromatic that RDKit mis-identifies."""
        mol = Molecule(smiles="[CH2]C=C1C=CC(=C)C2=C1C=CC=C2")
        out = generate_resonance_structures(mol)

        aromatic = Molecule().from_adjacency_list(
            """
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
"""
        )

        assert len(out) == 7
        assert any([m.is_isomorphic(aromatic) for m in out])

    def test_inconsistent_aromatic_structure_generation(self):
        """Test an unusual case of inconsistent aromaticity perception.

        Update: Highly-strained fused rings are no longer considered aromatic.
        That prevents the inconsistent aromatic structure for this molecule."""
        mol1 = Molecule().from_adjacency_list(
            """
multiplicity 2
1  C u0 p0 c0 {2,S} {6,S} {11,S} {12,S}
2  C u0 p0 c0 {1,S} {3,D} {4,S}
3  C u0 p0 c0 {2,D} {5,S} {7,S}
4  C u0 p0 c0 {2,S} {7,D} {10,S}
5  C u1 p0 c0 {3,S} {8,S} {9,S}
6  C u0 p0 c0 {1,S} {8,D} {13,S}
7  C u0 p0 c0 {3,S} {4,D} {17,S}
8  C u0 p0 c0 {5,S} {6,D} {14,S}
9  C u0 p0 c0 {5,S} {10,D} {15,S}
10 C u0 p0 c0 {4,S} {9,D} {16,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {7,S}
"""
        )

        mol2 = Molecule().from_adjacency_list(
            """
multiplicity 2
1  C u0 p0 c0 {2,S} {6,S} {11,S} {12,S}
2  C u0 p0 c0 {1,S} {3,D} {4,S}
3  C u0 p0 c0 {2,D} {5,S} {7,S}
4  C u0 p0 c0 {2,S} {7,D} {9,S}
5  C u1 p0 c0 {3,S} {8,S} {10,S}
6  C u0 p0 c0 {1,S} {8,D} {13,S}
7  C u0 p0 c0 {3,S} {4,D} {16,S}
8  C u0 p0 c0 {5,S} {6,D} {15,S}
9  C u0 p0 c0 {4,S} {10,D} {17,S}
10 C u0 p0 c0 {5,S} {9,D} {14,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
"""
        )

        # These two slightly different adjlists should be the same structure
        assert mol1.is_isomorphic(mol2)

        # However, they give different resonance structures
        res1 = generate_resonance_structures(mol1)
        res2 = generate_resonance_structures(mol2)
        assert res1 == res2

    def test_resonance_without_changing_atom_order1(self):
        """Test generating resonance structures without changing the atom order"""
        mol = Molecule().from_adjacency_list(
            """multiplicity 2
1  C u1 p0 c0 {2,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {3,S} {4,D}
3  C u0 p0 c0 {2,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {2,D} {5,S} {6,S}
5  C u0 p0 c0 {4,S} {13,S} {14,S} {15,S}
6  C u0 p0 c0 {4,S} {7,S} {16,S} {17,S}
7  O u0 p2 c0 {6,S} {18,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}"""
        )

        # Note: if save_order = False, atoms will be sorted
        # and the O atom will be reindexed to 1, which should be
        # true regardless of RMG's version and environment.
        # A copy is used to avoid the original mol's atoms are sorted.
        res_mols = mol.copy(deep=True).generate_resonance_structures(save_order=True)

        # Assign atom ids
        for molecule in [mol] + res_mols:
            for idx, atom in enumerate(molecule.atoms):
                atom.id = idx

        # Comparing atom symbol as its nearest neighbors
        for res_mol in res_mols:
            for atom1, atom2 in zip(mol.atoms, res_mol.atoms):
                assert atom1.element.symbol == atom2.element.symbol
                atom1_nb = {nb.id for nb in list(atom1.bonds.keys())}
                atom2_nb = {nb.id for nb in list(atom2.bonds.keys())}
                assert atom1_nb == atom2_nb

    def test_resonance_without_changing_atom_order2(self):
        """Test generating resonance structures for aromatic molecules without changing the atom order"""
        mol = Molecule().from_adjacency_list(
            """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {7,D}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,D} {8,S}
4  C u0 p0 c0 {3,D} {5,S} {9,S}
5  C u0 p0 c0 {4,S} {6,D} {10,S}
6  C u0 p0 c0 {5,D} {7,S} {11,S}
7  C u0 p0 c0 {1,D} {6,S} {12,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""
        )

        # Note: if save_order = False, atoms will be sorted
        # and the O atom will be reindexed to 1, which should be
        # true regardless of RMG's version and environment.
        # However, at commit 819779 and earlier versions, the atom order will be sorted
        # regardless the value of `save_order`
        res_mols = mol.copy(deep=True).generate_resonance_structures(save_order=True)

        # Assign atom ids
        for molecule in [mol] + res_mols:
            for idx, atom in enumerate(molecule.atoms):
                atom.id = idx

        # Comparing atom symbol as its nearest neighbors
        for res_mol in res_mols:
            for atom1, atom2 in zip(mol.atoms, res_mol.atoms):
                assert atom1.element.symbol == atom2.element.symbol
                atom1_nb = {nb.id for nb in list(atom1.bonds.keys())}
                atom2_nb = {nb.id for nb in list(atom2.bonds.keys())}
                assert atom1_nb == atom2_nb

    def test_adsorbate_resonance_cc1(self):
        """Test if all three resonance structures for X#CC#X are generated"""
        adjlist = """
1 X u0 p0 c0 {3,T}
2 X u0 p0 c0 {4,T}
3 C u0 p0 c0 {1,T} {4,S}
4 C u0 p0 c0 {2,T} {3,S}
        """

        mol = Molecule().from_adjacency_list(adjlist)

        mol_list = generate_resonance_structures(mol)
        assert len(mol_list) == 3

    def test_adsorbate_resonance_cc2(self):
        """Test if all three resonance structures for X=C=C=X are generated"""
        adjlist = """
1 X u0 p0 c0 {3,D}
2 X u0 p0 c0 {4,D}
3 C u0 p0 c0 {1,D} {4,D}
4 C u0 p0 c0 {2,D} {3,D}
        """

        mol = Molecule().from_adjacency_list(adjlist)

        mol_list = generate_resonance_structures(mol)
        assert len(mol_list) == 3

    def test_adsorbate_resonance_cc3(self):
        """Test if all three resonance structures for XC#CX are generated"""
        adjlist = """
1 X u0 p0 c0 {3,S}
2 X u0 p0 c0 {4,S}
3 C u0 p0 c0 {1,S} {4,T}
4 C u0 p0 c0 {2,S} {3,T}
        """

        mol = Molecule().from_adjacency_list(adjlist)

        mol_list = generate_resonance_structures(mol)
        assert len(mol_list) == 3

    def test_adsorbate_resonance_chchc(self):
        """Test if both resonance structures for XCHCHXC are generated"""
        adjlist = """
1 X u0 p0 c0 {3,T}
2 X u0 p0 c0 {5,S}
3 C u0 p0 c0 {1,T} {4,S}
4 C u0 p0 c0 {3,S} {5,D} {6,S}
5 C u0 p0 c0 {2,S} {4,D} {7,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
        """

        mol = Molecule().from_adjacency_list(adjlist)

        mol_list = generate_resonance_structures(mol)
        assert len(mol_list) == 2

class ClarTest:
    """
    Contains unit tests for Clar structure methods.
    """

    def test_clar_transformation(self):
        """Test that clarTransformation generates an aromatic ring."""
        mol = Molecule().from_smiles("c1ccccc1")
        sssr = mol.get_smallest_set_of_smallest_rings()
        _clar_transformation(mol, sssr[0])
        mol.update_atomtypes()

        assert mol.is_aromatic()

    def test_clar_optimization(self):
        """Test to ensure pi electrons are conserved during optimization"""
        mol = Molecule().from_smiles("C1=CC=C2C=CC=CC2=C1")  # Naphthalene
        output = _clar_optimization(mol)

        for molecule, asssr, bonds, solution in output:
            # Count pi electrons in molecule
            pi = 0
            for bond in bonds:
                if bond.is_double():
                    pi += 2

            # Count pi electrons in solution
            y = solution[0 : len(asssr)]
            x = solution[len(asssr) :]
            pi_solution = 6 * sum(y) + 2 * sum(x)

            # Check that both counts give 10 pi electrons
            assert pi == 10
            assert pi_solution == 10

            # Check that we only assign 1 aromatic sextet
            assert sum(y) == 1

    def test_phenanthrene(self):
        """Test that we generate 1 Clar structure for phenanthrene."""
        mol = Molecule().from_smiles("C1=CC=C2C(C=CC3=CC=CC=C32)=C1")
        newmol = generate_clar_structures(mol)

        struct = Molecule().from_adjacency_list(
            """1  C u0 p0 c0 {2,S} {3,B} {5,B}
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
"""
        )

        assert len(newmol) == 1
        assert newmol[0].is_isomorphic(struct)

    def test_phenalene(self):
        """Test that we generate 2 Clar structures for phenalene.

        Case where there is one non-aromatic ring."""
        mol = Molecule().from_smiles("C1=CC2=CC=CC3CC=CC(=C1)C=32")
        newmol = generate_clar_structures(mol)

        struct1 = Molecule().from_adjacency_list(
            """1  C u0 p0 c0 {2,S} {6,S} {14,S} {15,S}
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
"""
        )
        struct2 = Molecule().from_adjacency_list(
            """1  C u0 p0 c0 {2,S} {6,S} {14,S} {15,S}
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
"""
        )

        assert len(newmol) == 2
        assert newmol[0].is_isomorphic(struct1) or newmol[0].is_isomorphic(struct2)
        assert newmol[1].is_isomorphic(struct2) or newmol[1].is_isomorphic(struct1)
        assert not newmol[0].is_isomorphic(newmol[1])

    def test_corannulene(self):
        """Test that we generate 5 Clar structures for corannulene

        Case where linear relaxation does not give an integer solution"""
        mol = Molecule().from_smiles("C1=CC2=CC=C3C=CC4=C5C6=C(C2=C35)C1=CC=C6C=C4")
        newmol = generate_clar_structures(mol)

        struct = Molecule().from_adjacency_list(
            """1  C u0 p0 c0 {2,S} {5,B} {8,B}
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
"""
        )

        assert len(newmol) == 5
        assert newmol[0].is_isomorphic(struct)
        assert newmol[1].is_isomorphic(struct)
        assert newmol[2].is_isomorphic(struct)
        assert newmol[3].is_isomorphic(struct)
        assert newmol[4].is_isomorphic(struct)

    def test_exocyclic_db(self):
        """Test that Clar structure generation doesn't modify exocyclic double bonds

        Important for cases where RDKit considers rings to be aromatic by counting pi-electron contributions
        from exocyclic double bonds, while they don't actually contribute to aromaticity
        """

        mol = Molecule(smiles="C=C1C=CC=CC1=C")
        newmol = generate_clar_structures(mol)

        assert len(newmol) == 0

    def test_surface_o(self):
        """Test resonance structure generation for surface adsorbed O=X

        Should not crash."""
        # See https://github.com/cfgoldsmith/RMG-Py/issues/43
        mol_list = generate_resonance_structures(
            Molecule().from_adjacency_list(
                """OX
1 X u0 p0 c0 {2,D}
2 O u0 p2 c0 {1,D}"""
            )
        )
        assert len(mol_list) == 1

    def test_surface_radical(self):
        """Test resonance structure generation for radical on surface

        Should not put radical on X atom"""
        mol_list = generate_resonance_structures(Molecule(smiles="*=C[CH]C"))
        assert len(mol_list) == 1

    def test_sulfur_triple_bond(self):
        """
        Test the prevention of S#S formation through the find_lone_pair_multiplebond_paths and
        find_adj_lone_pair_multiple_bond_delocalization_paths
        """
        mol_list = generate_resonance_structures(Molecule(smiles="S1SSS1"), filter_structures=False)
        assert len(mol_list) == 10
