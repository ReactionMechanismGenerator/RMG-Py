#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2021 Prof. William H. Green (whgreen@mit.edu),           #
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
from rmgpy.molecule.resonance import generate_optimal_aromatic_resonance_structures
from rmgpy.molecule.symmetry import (
    calculate_atom_symmetry_number,
    calculate_axis_symmetry_number,
    calculate_bond_symmetry_number,
    calculate_cyclic_symmetry_number,
    _indistinguishable,
)
from rmgpy.species import Species


class TestMoleculeSymmetry:
    """
    Contains unit tests of the methods for computing symmetry numbers for a
    given Molecule object.
    """

    def test_atom_symmetry_number_methane(self):
        """
        Test the Molecule.calculate_atom_symmetry_number() on CH4
        """
        molecule = Molecule().from_smiles("C")
        symmetry_number = 1
        for atom in molecule.atoms:
            if not molecule.is_atom_in_cycle(atom):
                symmetry_number *= calculate_atom_symmetry_number(molecule, atom)
        assert symmetry_number == 12

    def test_atom_symmetry_number_methyl(self):
        """
        Test the Molecule.calculate_atom_symmetry_number() on [CH3]
        """
        molecule = Molecule().from_smiles("[CH3]")
        symmetry_number = 1
        for atom in molecule.atoms:
            if not molecule.is_atom_in_cycle(atom):
                symmetry_number *= calculate_atom_symmetry_number(molecule, atom)
        assert symmetry_number == 6

    def test_atom_symmetry_number_ethane(self):
        """
        Test the Molecule.calculate_atom_symmetry_number() on CC
        """
        molecule = Molecule().from_smiles("CC")
        symmetry_number = 1
        for atom in molecule.atoms:
            if not molecule.is_atom_in_cycle(atom):
                symmetry_number *= calculate_atom_symmetry_number(molecule, atom)
        assert symmetry_number == 9

    def test_atom_symmetry_number_ethane_with_deuterium_tritium(self):
        """
        Test the Molecule.calculate_atom_symmetry_number() on CC(D)(T)

        This is meant to test whether chirality is accounted for, which
        should half the symmetry term.

        The total number is 1.5 because the methyl group contributes *3 and
        the chiral center contributes *0.5
        """
        molecule = Molecule().from_adjacency_list(
            """
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3 D u0 p0 c0 {1,S}
4 T u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
"""
        )
        symmetry_number = 1
        for atom in molecule.atoms:
            if not molecule.is_atom_in_cycle(atom):
                symmetry_number *= calculate_atom_symmetry_number(molecule, atom)
        assert round(abs(symmetry_number - 1.5), 7) == 0

    def test_atom_symmetry_number_with_two_chiral_centers(self):
        """
        Test the Molecule.calculate_atom_symmetry_number() on [CH2]CC([CH2])C(C)C=C

        This is meant to test whether chirality is accounted for, which
        should half the symmetry term.

        The molecule has one methyl group (*3), two CH2dot groups (*2 each), and
        two chiral centers (*0.5), leading to a total atom symmetry number of 3
        """
        molecule = Molecule().from_adjacency_list(
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
"""
        )
        symmetry_number = 1
        for atom in molecule.atoms:
            if not molecule.is_atom_in_cycle(atom):
                symmetry_number *= calculate_atom_symmetry_number(molecule, atom)
        assert round(abs(symmetry_number - 3), 7) == 0

    def test_atom_symmetry_number_propane(self):
        """
        Test the Molecule.calculate_atom_symmetry_number() on CCC
        """
        molecule = Molecule().from_smiles("CCC")
        symmetry_number = 1
        for atom in molecule.atoms:
            if not molecule.is_atom_in_cycle(atom):
                symmetry_number *= calculate_atom_symmetry_number(molecule, atom)
        assert symmetry_number == 18

    def test_atom_symmetry_number_isobutane(self):
        """
        Test the Molecule.calculate_atom_symmetry_number() on CC(C)C
        """
        molecule = Molecule().from_smiles("CC(C)C")
        symmetry_number = 1
        for atom in molecule.atoms:
            if not molecule.is_atom_in_cycle(atom):
                symmetry_number *= calculate_atom_symmetry_number(molecule, atom)
        assert symmetry_number == 81

    def test_bond_symmetry_number_ethane(self):
        """
        Test the Molecule.calculate_bond_symmetry_number() on CC
        """
        molecule = Molecule().from_smiles("CC")
        symmetry_number = 1
        for atom1 in molecule.atoms:
            for atom2 in list(atom1.bonds):
                if molecule.atoms.index(atom1) < molecule.atoms.index(atom2):
                    symmetry_number *= calculate_bond_symmetry_number(molecule, atom1, atom2)
        assert symmetry_number == 2

    def test_bond_symmetry_number_propane(self):
        """
        Test the Molecule.calculate_bond_symmetry_number() on CCC
        """
        molecule = Molecule().from_smiles("CCC")
        symmetry_number = 1
        for atom1 in molecule.atoms:
            for atom2 in list(atom1.bonds):
                if molecule.atoms.index(atom1) < molecule.atoms.index(atom2):
                    symmetry_number *= calculate_bond_symmetry_number(molecule, atom1, atom2)
        assert symmetry_number == 1

    def test_bond_symmetry_number_butane(self):
        """
        Test the Molecule.calculate_bond_symmetry_number() on CCCC
        """
        molecule = Molecule().from_smiles("CCCC")
        symmetry_number = 1
        for atom1 in molecule.atoms:
            for atom2 in list(atom1.bonds):
                if molecule.atoms.index(atom1) < molecule.atoms.index(atom2):
                    symmetry_number *= calculate_bond_symmetry_number(molecule, atom1, atom2)
        assert symmetry_number == 2

    def test_bond_symmetry_number_ethylene(self):
        """
        Test the Molecule.calculate_bond_symmetry_number() on C=C
        """
        molecule = Molecule().from_smiles("C=C")
        symmetry_number = 1
        for atom1 in molecule.atoms:
            for atom2 in list(atom1.bonds):
                if molecule.atoms.index(atom1) < molecule.atoms.index(atom2):
                    symmetry_number *= calculate_bond_symmetry_number(molecule, atom1, atom2)
        assert symmetry_number == 2

    def test_bond_symmetry_number_acetylene(self):
        """
        Test the Molecule.calculate_bond_symmetry_number() on C#C
        """
        molecule = Molecule().from_smiles("C#C")
        symmetry_number = 1
        for atom1 in molecule.atoms:
            for atom2 in list(atom1.bonds):
                if molecule.atoms.index(atom1) < molecule.atoms.index(atom2):
                    symmetry_number *= calculate_bond_symmetry_number(molecule, atom1, atom2)
        assert symmetry_number == 2

    def test_axis_symmetry_number_ethylene(self):
        """
        Test the Molecule.calculate_axis_symmetry_number() on C=C
        """
        molecule = Molecule().from_smiles("C=C")
        assert calculate_axis_symmetry_number(molecule) == 2

    def test_axis_symmetry_number_propadiene(self):
        """
        Test the Molecule.calculate_axis_symmetry_number() on C=C=C
        """
        molecule = Molecule().from_smiles("C=C=C")
        assert calculate_axis_symmetry_number(molecule) == 2

    def test_axis_symmetry_number_butatriene(self):
        """
        Test the Molecule.calculate_axis_symmetry_number() on C=C=C=C
        """
        molecule = Molecule().from_smiles("C=C=C=C")
        assert calculate_axis_symmetry_number(molecule) == 2

    def test_axis_symmetry_number_butatrienyl(self):
        """
        Test the Molecule.calculate_axis_symmetry_number() on C=C=C=[CH]
        """
        molecule = Molecule().from_smiles("C=C=C=[CH]")
        assert calculate_axis_symmetry_number(molecule) == 1

    def test_axis_symmetry_number_propadienyl(self):
        """
        Test the Molecule.calculate_axis_symmetry_number() on C=C=[C]
        """
        molecule = Molecule().from_smiles("C=C=[C]")
        assert calculate_axis_symmetry_number(molecule) == 2

    def test_axis_symmetry_number_12_butadienyl(self):
        """
        Test the Molecule.calculate_axis_symmetry_number() on CC=C=[C]
        """
        molecule = Molecule().from_smiles("CC=C=[C]")
        assert calculate_axis_symmetry_number(molecule) == 1

    def test_axis_symmetry_number_12_hexadienyl(self):
        """
        Test the Molecule.calculate_axis_symmetry_number() on C=C=CCCC
        """
        molecule = Molecule().from_smiles("C=C=CCCC")
        assert calculate_axis_symmetry_number(molecule) == 1

    def test_axis_symmetry_number1(self):
        """
        Test the Molecule.calculate_axis_symmetry_number() on CC(C)=C=C(CC)CC
        """
        molecule = Molecule().from_smiles("CC(C)=C=C(CC)CC")
        assert calculate_axis_symmetry_number(molecule) == 2

    def test_axis_symmetry_number2(self):
        """
        Test the Molecule.calculate_axis_symmetry_number() on C=C=C(C(C(C(C=C=C)=C=C)=C=C)=C=C)
        """
        molecule = Molecule().from_smiles("C=C=C(C(C(C(C=C=C)=C=C)=C=C)=C=C)")
        assert calculate_axis_symmetry_number(molecule) == 2

    @pytest.mark.skip(reason="WIP")
    def test_axis_symmetry_number3(self):
        """
        Test the Molecule.calculate_axis_symmetry_number() on C=C=[C]C(C)(C)[C]=C=C
        """
        molecule = Molecule().from_smiles("C=C=[C]C(C)(C)[C]=C=C")
        assert calculate_axis_symmetry_number(molecule) == 4

    def test_axis_symmetry_number4(self):
        """
        Test the Molecule.calculate_axis_symmetry_number() on C=C=C=O
        """
        molecule = Molecule().from_smiles("C=C=C=O")
        assert calculate_axis_symmetry_number(molecule) == 2

    def test_axis_symmetry_number5(self):
        """
        Test the Molecule.calculate_axis_symmetry_number() on CC=C=C=O
        """
        molecule = Molecule().from_smiles("CC=C=C=O")
        assert calculate_axis_symmetry_number(molecule) == 1

    def test_axis_symmetry_number6(self):
        """
        Test the Molecule.calculate_axis_symmetry_number() on C=C=C=N
        """
        molecule = Molecule().from_smiles("C=C=C=N")
        assert calculate_axis_symmetry_number(molecule) == 1

    def test_axis_symmetry_number7(self):
        """
        Test the Molecule.calculate_axis_symmetry_number() on C=C=C=[N]
        """
        molecule = Molecule().from_smiles("C=C=C=[N]")
        assert calculate_axis_symmetry_number(molecule) == 2

    def test_axis_symmetry_oxygen_singlet(self):
        """
        Test the Molecule.calculate_axis_symmetry_number() on O=O
        """
        molecule = Molecule().from_smiles("O=O")
        assert calculate_axis_symmetry_number(molecule) == 1

    def test_cyclic_symmetry_number_cyclohexane(self):
        """
        Test the Molecule.calculate_cyclic_symmetry_number() on C1CCCCC1
        """
        molecule = Molecule().from_smiles("C1CCCCC1")
        symmetry_number = calculate_cyclic_symmetry_number(molecule)
        assert symmetry_number == 12

    def test_cyclic_symmetry_number_cyclohexanone(self):
        """
        Test the Molecule.calculate_cyclic_symmetry_number() on C1CCCCC1=O
        """
        molecule = Molecule().from_smiles("C1CCCCC1=O")
        symmetry_number = calculate_cyclic_symmetry_number(molecule)
        assert symmetry_number == 2

    def test_cyclic_symmetry_number_cyclohexan_tri_one(self):
        """
        Test the Molecule.calculate_cyclic_symmetry_number() on C1CCC(=O)C(=O)C1=O
        """
        molecule = Molecule().from_smiles("C1CCC(=O)C(=O)C1=O")
        symmetry_number = calculate_cyclic_symmetry_number(molecule)
        assert symmetry_number == 2

    def test_total_symmetry_number_benzene(self):
        """
        Test the Species.get_symmetry_number() (total symmetry) on c1ccccc1
        """
        molecule = Molecule().from_smiles("c1ccccc1")
        species = Species(molecule=[molecule])
        symmetry_number = species.get_symmetry_number()
        assert symmetry_number == 12

    def test_total_symmetry_numberbute_di_yl(self):
        """
        Test the Species.get_symmetry_number() (total symmetry) on [CH2][CH]C=C
        """
        species = Species().from_smiles("[CH2][CH]C=C")
        symmetry_number = species.get_symmetry_number()
        assert symmetry_number == 2

    def test_total_symmetry_number_allyl(self):
        """
        Test the Species.get_symmetry_number() (total symmetry) on Allyl, [CH2]C=C
        """
        molecule = Molecule().from_smiles("[CH2]C=C")
        species = Species(molecule=[molecule])
        symmetry_number = species.get_symmetry_number()
        assert symmetry_number == 2

    def test_total_symmetry_number_pentenyl(self):
        """
        Test the Species.getSymmetryNumer() for C[CH]C=CC
        and ensures that it is differet than the molecule object
        """
        spc = Species(molecule=[Molecule().from_smiles("C[CH]C=CC")])
        symmetry_number = spc.get_symmetry_number()
        assert symmetry_number == 18, "the symmetry number for C[CH]C=CC is 18, but RMG returned {}".format(symmetry_number)

    def test_species_symmetry_number_is_not_molecule_symmetry_number(self):
        """
        Tests that the species symmetry number can be different from the molecule symmetry number

        This molecule's resonance isomer hybrid should return more symmetry
        than the base molecule object.
        """
        molecule = Molecule().from_smiles("C[CH]C=CC")
        species = Species(molecule=[molecule])
        assert molecule.get_symmetry_number() * 2 == species.get_symmetry_number()

    def test_axis_symmetry_number_allyl(self):
        """
        Test the Molecule.calculate_axis_symmetry_number() on [CH2]C=C
        """
        spc = Species(molecule=[Molecule().from_smiles("[CH2]C=C")])
        molecule = spc.get_resonance_hybrid()
        assert calculate_axis_symmetry_number(molecule) == 1

    def test_middle_carbon_atom_symmetry_number_allyl_using_resonance_hybrid(self):
        """
        Test the Molecule.calculate_atom_symmetry_number() for the middle carbon
        on [CH2]C=C using the resonance hybrid structure
        """
        molecule = Molecule().from_smiles("[CH2]C=C")
        species = Species(molecule=[molecule])
        resonance_hybrid = species.get_resonance_hybrid()
        for atom in resonance_hybrid.atoms:
            if atom.symbol == "C":
                number_carbon_bonds = sum([1 for bond in atom.bonds if bond.symbol == "C"])
                if number_carbon_bonds == 2:
                    atom.label = "center"
                    symmetry_number = calculate_atom_symmetry_number(resonance_hybrid, atom)
                    assert symmetry_number == 2
                    pass

    def test_edge_carbon_atom_symmetry_number_allyl(self):
        """
        Test the Molecule.calculate_atom_symmetry_number() for the edge carbons
        on [CH2]C=C
        """
        spc = Species(molecule=[Molecule().from_smiles("[CH2]C=C")])
        molecule = spc.get_resonance_hybrid()
        for atom in molecule.atoms:
            if atom.symbol == "C":
                number_carbon_bonds = sum([1 for bond in atom.bonds if bond.symbol == "C"])
                if number_carbon_bonds == 1:
                    atom.label = "edge"
                    symmetry_number = calculate_atom_symmetry_number(molecule, atom)
                    assert symmetry_number == 1

    def test_edge_carbon_atom_symmetry_number_allyl_using_resonance_hybrid(self):
        """
        Test the Molecule.calculate_atom_symmetry_number() for the edge carbons
        on [CH2]C=C
        """
        molecule = Molecule().from_smiles("[CH2]C=C")
        species = Species(molecule=[molecule])
        resonance_hybrid = species.get_resonance_hybrid()
        for atom in resonance_hybrid.atoms:
            if atom.symbol == "C":
                number_carbon_bonds = sum([1 for bond in atom.bonds if bond.symbol == "C"])
                if number_carbon_bonds == 1:
                    atom.label = "edge"
                    symmetry_number = calculate_atom_symmetry_number(resonance_hybrid, atom)
                    assert symmetry_number == 1

    def test_atom_symmetry_number_allyl(self):
        """
        Test the Molecule.calculate_atom_symmetry_number() on [CH2]C=C
        """
        spc = Species(molecule=[Molecule().from_smiles("[CH2]C=C")])
        molecule = spc.get_resonance_hybrid()
        symmetry_number = 1
        for atom in molecule.atoms:
            if not molecule.is_atom_in_cycle(atom):
                symmetry_number *= calculate_atom_symmetry_number(molecule, atom)
        assert symmetry_number == 2

    def test_bond_symmetry_number_allyl(self):
        """
        Test the Molecule.calculate_bond_symmetry_number() on [CH2]C=C
        """
        spc = Species(molecule=[Molecule().from_smiles("[CH2]C=C")])
        molecule = spc.get_resonance_hybrid()
        symmetry_number = 1
        for atom1 in molecule.atoms:
            for atom2 in atom1.bonds:
                if molecule.atoms.index(atom1) < molecule.atoms.index(atom2):
                    symmetry_number *= calculate_bond_symmetry_number(molecule, atom1, atom2)
        assert symmetry_number == 1

    def test_total_symmetry_number_toluene(self):
        """
        Test the Species.get_symmetry_number() (total symmetry) on c1ccccc1C
        """
        molecule = Molecule().from_smiles("c1ccccc1C")
        species = Species(molecule=[molecule])
        symmetry_number = species.get_symmetry_number()
        assert symmetry_number == 6

    def test_total_symmetry_number_special_cyclic(self):
        """
        Test the Species.get_symmetry_number() (total symmetry) from issue # 332
        """
        molecule = Molecule().from_smiles("C1(C(C(C(C(C1C2CCC2)C3CCC3)C4CCC4)C5CCC5)C6CCC6)C7CCC7")
        species = Species(molecule=[molecule])
        symmetry_number = species.get_symmetry_number()
        assert symmetry_number == 12

    def test_total_symmetry_number_chlorobenzene(self):
        """
        Test the Species.get_symmetry_number() (total symmetry) on c1ccccc1Cl
        """
        molecule = Molecule().from_smiles("c1ccccc1Cl")
        species = Species(molecule=[molecule])
        symmetry_number = species.get_symmetry_number()
        assert symmetry_number == 2

    def test_total_symmetry_number_phenoxy_kecle(self):
        """
        Test the Species.get_symmetry_number() (total symmetry) on c1ccccc1[O]
        using kecle structure
        """
        molecule = Molecule().from_smiles("c1ccccc1[O]")
        species = Species(molecule=[molecule])
        symmetry_number = species.get_symmetry_number()
        assert symmetry_number == 2

    def test_total_symmetry_number_phenoxy_benzene(self):
        """
        Test symmetry on c1ccccc1[O] using phenoxy benzene structure
        """
        molecule = Molecule().from_smiles("c1ccccc1[O]")
        aromatic_molecule = generate_optimal_aromatic_resonance_structures(molecule)[0]
        symmetry_number = aromatic_molecule.get_symmetry_number()
        assert symmetry_number == 2

    def test_total_symmetry_number12_dimethylbenzene(self):
        """
        Test the Species.get_symmetry_number() (total symmetry) on Cc1ccccc1C
        """
        molecule = Molecule().from_smiles("Cc1ccccc1C")
        species = Species(molecule=[molecule])
        symmetry_number = species.get_symmetry_number()
        assert symmetry_number == 18

    def test_total_symmetry_number14_dimethylbenzene(self):
        """
        Test the Species.get_symmetry_number() (total symmetry) on Cc1ccc(C)cc1
        """
        molecule = Molecule().from_smiles("Cc1ccc(C)cc1")
        species = Species(molecule=[molecule])
        symmetry_number = species.get_symmetry_number()
        assert symmetry_number == 36

    def test_total_symmetry_number_ethane(self):
        """
        Test the Molecule.calculate_symmetry_number() on CC
        """
        assert Species().from_smiles("CC").get_symmetry_number() == 18

    def test_total_symmetry_number1(self):
        """
        Test the Molecule.calculate_symmetry_number() on C=C=[C]C(C)(C)[C]=C=C
        """
        assert Species().from_smiles("C=C=[C]C(C)(C)[C]=C=C").get_symmetry_number() == 18

    @pytest.mark.skip(reason="WIP")
    def test_total_symmetry_number2(self):
        """
        Test the Molecule.calculate_symmetry_number() on C(=CC(c1ccccc1)C([CH]CCCCCC)C=Cc1ccccc1)[CH]CCCCCC
        """
        assert Species().from_smiles("C(=CC(c1ccccc1)C([CH]CCCCCC)C=Cc1ccccc1)[CH]CCCCCC").get_symmetry_number() == "36?"

    def test_symmetry_number_hydroxyl(self):
        """
        Test the Molecule.calculate_symmetry_number() on [OH]
        """
        assert Species().from_smiles("[OH]").get_symmetry_number() == 1

    def test_symmetry_number_oxygen(self):
        """
        Test the Molecule.calculate_symmetry_number() on O=O
        """
        assert Species().from_smiles("O=O").get_symmetry_number() == 2

    def test_symmetry_number_dicarbon(self):
        """
        Test the Molecule.calculate_symmetry_number() on [C]#[C]
        """
        assert Species().from_smiles("[C]#[C]").get_symmetry_number() == 2

    def test_symmetry_number_hydrogen(self):
        """
        Test the Molecule.calculate_symmetry_number() on [H][H]
        """
        assert Species().from_smiles("[H][H]").get_symmetry_number() == 2

    def test_symmetry_number_acetylene(self):
        """
        Test the Molecule.calculate_symmetry_number() on C#C
        """
        assert Species().from_smiles("C#C").get_symmetry_number() == 2

    def test_symmetry_number_butadiyne(self):
        """
        Test the Molecule.calculate_symmetry_number() on C#CC#C
        """
        assert Species().from_smiles("C#CC#C").get_symmetry_number() == 2

    def test_symmetry_number_methane(self):
        """
        Test the Molecule.calculate_symmetry_number() on CH4
        """
        assert Species().from_smiles("C").get_symmetry_number() == 12

    def test_symmetry_number_formaldehyde(self):
        """
        Test the Molecule.calculate_symmetry_number() on C=O
        """
        assert Species().from_smiles("C=O").get_symmetry_number() == 2

    def test_symmetry_number_methyl(self):
        """
        Test the Molecule.calculate_symmetry_number() on [CH3]
        """
        assert Species().from_smiles("[CH3]").get_symmetry_number() == 6

    def test_symmetry_number_water(self):
        """
        Test the Molecule.calculate_symmetry_number() on H2O
        """
        assert Species().from_smiles("O").get_symmetry_number() == 2

    def test_symmetry_number_ethylene(self):
        """
        Test the Molecule.calculate_symmetry_number() on C=C
        """
        assert Species().from_smiles("C=C").get_symmetry_number() == 4

    def test_symmetry_number_ethenyl(self):
        """
        Test the Molecule.calculate_symmetry_number() on C=[CH].
        """
        assert Species().from_smiles("C=[CH]").get_symmetry_number() == 1

    def test_symmetry_number_cyclic(self):
        """
        Test the Molecule.calculate_symmetry_number() on C1=C=C=1
        """
        assert Species().from_smiles("C1=C=C=1").get_symmetry_number() == 6

    def test_symmetry_number__s1_total(self):
        """
        Test the Molecule.calculate_symmetry_number() on [CH]1CCC1CC1CC1
        """
        assert Species().from_smiles("[CH]1CCC1CC1CC1").get_symmetry_number() == 1

    def test_cyclic_symmetry_number_s1(self):
        """
        Test the Molecule.calculate_cyclic_symmetry_number() on [CH]1CCC1CC1CC1
        """
        molecule = Molecule().from_smiles("[CH]1CCC1CC1CC1")
        symmetry_number = calculate_cyclic_symmetry_number(molecule)
        assert symmetry_number == 1

    def test_cyclic_symmetry_number_methyl_cyclo_propane(self):
        """
        Test the Molecule.calculate_cyclic_symmetry_number() on CC1CC1
        """
        molecule = Molecule().from_smiles("CC1CC1")
        symmetry_number = calculate_cyclic_symmetry_number(molecule)
        assert symmetry_number == 1

    def test_cyclic_symmetry_number_methyl_cyclo_propene(self):
        """
        Test the Molecule.calculate_cyclic_symmetry_number() on C=C1CC1
        """
        molecule = Molecule().from_smiles("C=C1CC1")
        symmetry_number = calculate_cyclic_symmetry_number(molecule)
        assert symmetry_number == 2

    def test_cyclic_symmetry_number_methyl_cyclo_butene(self):
        """
        Test the Molecule.calculate_cyclic_symmetry_number() on C=C1CCC1
        """
        molecule = Molecule().from_smiles("C=C1CCC1")
        symmetry_number = calculate_cyclic_symmetry_number(molecule)
        assert symmetry_number == 2

    def test_cyclic_symmetry_number_methyl_cyclo_butane(self):
        """
        Test the Molecule.calculate_cyclic_symmetry_number() on CC1CCC1
        """
        molecule = Molecule().from_smiles("CC1CCC1")
        symmetry_number = calculate_cyclic_symmetry_number(molecule)
        assert symmetry_number == 1

    def test_cyclic_symmetry_number_di_methyl_cyclo_butane(self):
        """
        Test the Molecule.calculate_cyclic_symmetry_number() on CC1CC(C)C1
        """
        molecule = Molecule().from_smiles("CC1CC(C)C1")
        symmetry_number = calculate_cyclic_symmetry_number(molecule)
        assert symmetry_number == 4

    def test_symmetry_number_dimethylcylcobutane_total(self):
        """
        Test the Molecule.calculate_symmetry_number() on CC1CC(C)C1
        """
        assert Species().from_smiles("CC1CC(C)C1").get_symmetry_number() == 36

    def test_indistinguishable(self):
        """
        Test that the _indistinguishable function works properly
        """
        mol = Molecule().from_smiles("c1ccccc1")
        assert _indistinguishable(mol.atoms[0], mol.atoms[1])
        assert _indistinguishable(mol.atoms[0], mol.atoms[2])
        assert _indistinguishable(mol.atoms[0], mol.atoms[3])

    def test_indistinguishable_2(self):
        """
        Test that the _indistinguishable function works properly
        """
        mol = Molecule().from_smiles("c1ccccc1[O]")
        # Carbon with O different from other carbons (with H)
        assert not _indistinguishable(mol.atoms[5], mol.atoms[4])
        assert not _indistinguishable(mol.atoms[5], mol.atoms[3])
        assert not _indistinguishable(mol.atoms[5], mol.atoms[2])
        # Ortho carbons are the same
        assert _indistinguishable(mol.atoms[0], mol.atoms[4])
        # Meta carbons are the same
        assert _indistinguishable(mol.atoms[1], mol.atoms[3])
        # O is different from H
        assert not _indistinguishable(mol.atoms[6], mol.atoms[7])
