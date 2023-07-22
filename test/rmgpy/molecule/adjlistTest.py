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

import logging


from rmgpy.molecule.adjlist import InvalidAdjacencyListError
from rmgpy.molecule.group import Group
from rmgpy.molecule.molecule import Molecule
import pytest

logging.basicConfig(level=logging.DEBUG)


class TestGroupAdjLists:
    """
    Contains adjacency list unit tests of the Graph class.
    """

    def setup_class(self):
        pass

    def test_from_old_adjacency_list1(self):
        """
        adjlist: Test the Group.from_adjacency_list() method on an old style adjacency list.
        """
        adjlist = """
1 *2 {Cs,Cd} 0 {2,{S,D}} {3,S}
2 *1 {O2s,O2d}  0   {1,{S,D}}
3    R!H     {0,1} {1,S}
            """
        group = Group().from_adjacency_list(adjlist)

        atom1, atom2, atom3 = group.atoms
        assert group.has_bond(atom1, atom2)
        assert group.has_bond(atom1, atom3)
        assert not group.has_bond(atom2, atom3)
        bond12 = atom1.bonds[atom2]
        bond13 = atom1.bonds[atom3]

        assert atom1.label == "*2"
        assert atom1.atomtype[0].label in ["Cs", "Cd"]
        assert atom1.atomtype[1].label in ["Cs", "Cd"]
        assert atom1.radical_electrons == [0]

        assert atom2.label == "*1"
        assert atom2.atomtype[0].label in ["O2s", "O2d"]
        assert atom2.atomtype[1].label in ["O2s", "O2d"]
        assert atom2.radical_electrons == [0]

        assert atom3.label == ""
        assert atom3.atomtype[0].label == "R!H"
        assert atom3.radical_electrons == [0, 1]

        assert bond12.order == [1, 2]
        assert bond13.is_single()

    def test_from_adjacency_list(self):
        """
        adjlist: Test the Group.from_adjacency_list() method.
        """
        adjlist = """
1 *2 [Cs,Cd]   u0 {2,[S,D]} {3,S}
2 *1 [O2s,O2d] u0 {1,[S,D]}
3    R!H       u0 {1,S}
            """
        group = Group().from_adjacency_list(adjlist)

        atom1, atom2, atom3 = group.atoms
        assert group.has_bond(atom1, atom2)
        assert group.has_bond(atom1, atom3)
        assert not group.has_bond(atom2, atom3)
        bond12 = atom1.bonds[atom2]
        bond13 = atom1.bonds[atom3]

        assert atom1.label == "*2"
        assert atom1.atomtype[0].label in ["Cs", "Cd"]
        assert atom1.atomtype[1].label in ["Cs", "Cd"]
        assert atom1.radical_electrons == [0]

        assert atom2.label == "*1"
        assert atom2.atomtype[0].label in ["O2s", "O2d"]
        assert atom2.atomtype[1].label in ["O2s", "O2d"]
        assert atom2.radical_electrons == [0]

        assert atom3.label == ""
        assert atom3.atomtype[0].label == "R!H"
        assert atom3.radical_electrons == [0]

        assert bond12.order == [1, 2]
        assert bond13.is_single()

    def test_from_adjacency_list_multiplicity(self):
        gp = Group().from_adjacency_list(
            """
            multiplicity [1]
            1 C u0 p0 c0
            """
        )
        assert len(gp.multiplicity) == 1
        assert gp.multiplicity[0] == 1

    def test_from_adjacency_list_multiplicity_list(self):
        gp = Group().from_adjacency_list(
            """
            multiplicity [ 1, 3, 5 ]
            1 C u0 p0 c0
            """
        )
        assert len(gp.multiplicity) == 3
        assert gp.multiplicity[0] == 1
        assert gp.multiplicity[1] == 3
        assert gp.multiplicity[2] == 5

    def test_to_adjacency_list(self):
        """
        adjlist: Test the Group.to_adjacency_list() method.
        """
        adjlist = """
1 *2 [Cs,Cd]   u0 {2,[S,D]} {3,S}
2 *1 [O2s,O2d] u0 {1,[S,D]}
3    R!H       u0 {1,S}
            """
        group = Group().from_adjacency_list(adjlist)
        adjlist2 = group.to_adjacency_list()

        assert adjlist.strip() == adjlist2.strip()

    def test_atom_props(self):
        """Test that the atom props attribute can be properly read and written."""
        adjlist = """
1 *1 R!H u1 r0 {2,S}
2 *4 R!H u0 r0 {1,S} {3,S}
3 *2 Cb  u0 r1 {2,S} {4,B}
4 *3 Cb  u0 r1 {3,B}
        """
        group = Group().from_adjacency_list(adjlist)
        for atom in group.atoms:
            if atom.atomtype[0].label == "R!H":
                assert not atom.props["inRing"]
            elif atom.atomtype[0].label == "Cb":
                assert atom.props["inRing"]
        adjlist2 = group.to_adjacency_list()

        assert adjlist.strip() == adjlist2.strip()

    def test_metal_facet_site_morphology(self):
        adjlist1 = """metal Cu
facet 111
1 *3 X u0 p0 c0 s"atop" m"terrace"  {2,S} {4,S}
2 *1 O u0 p2 c0 {1,S} {3,R}
3 *2 H u0 p0 c0 {2,R} {4,R}
4 *4 X u0 p0 c0 s"hcp" m"terrace" {3,R} {1,S}"""

        adjlist2 = """multiplicity [1]
metal [Cu, Fe, CuO2 ]
facet [111, 211, 1101, 110, ]
1 *3 X u0 p0 c0 s["atop","fcc"] m"terrace"  {2,S} {4,S}
2 *1 O u0 p2 c0 {1,S} {3,R}
3 *2 H u0 p0 c0 {2,R} {4,R}
4 *4 X u0 p0 c0 s"hcp" m["terrace","sc"] {3,R} {1,S}"""

        adjlist1test = """metal Cu
facet 111
1 *3 X u0 p0 c0 s"atop" m"terrace" {2,S} {4,S}
2 *1 O u0 p2 c0 {1,S} {3,R}
3 *2 H u0 p0 c0 {2,R} {4,R}
4 *4 X u0 p0 c0 s"hcp" m"terrace" {1,S} {3,R}"""

        adjlist2test = """multiplicity [1]
metal [Cu,Fe,CuO2]
facet [111,211,1101,110]
1 *3 X u0 p0 c0 s["atop","fcc"] m"terrace" {2,S} {4,S}
2 *1 O u0 p2 c0 {1,S} {3,R}
3 *2 H u0 p0 c0 {2,R} {4,R}
4 *4 X u0 p0 c0 s"hcp" m["terrace","sc"] {1,S} {3,R}"""
        mol = Molecule().from_adjacency_list(adjlist1, check_consistency=False)
        group = Group().from_adjacency_list(adjlist2)

        assert mol.metal == "Cu"
        assert mol.facet == "111"
        assert group.metal == ["Cu", "Fe", "CuO2"]
        assert group.facet == ["111", "211", "1101", "110"]

        assert mol.atoms[0].site == "atop"
        assert mol.atoms[3].site == "hcp"
        assert mol.atoms[0].morphology == "terrace"
        assert mol.atoms[3].morphology == "terrace"

        assert group.atoms[0].site == ["atop", "fcc"]
        assert group.atoms[3].site == ["hcp"]
        assert group.atoms[0].morphology == ["terrace"]
        assert group.atoms[3].morphology == ["terrace", "sc"]

        assert mol.to_adjacency_list().strip() == adjlist1test.strip()
        assert group.to_adjacency_list().strip() == adjlist2test.strip()


class TestMoleculeAdjLists:
    """
    adjlist: Contains adjacency list unit tests of the Molecule class.
    """

    def setup_class(self):
        pass

    def test_from_adjacency_list1(self):
        """
        adjlist: Test the Molecule.from_adjacency_list() method 1.
        """
        # molecule 1
        adjlist = """
1 *1 C u1 p0 c0  {2,S} {3,S} {4,S}
2    H u0 p0 c0  {1,S}
3    H u0 p0 c0  {1,S}
4 *2 N u0 p0 c+1 {1,S} {5,S} {6,D}
5    O u0 p3 c-1 {4,S}
6    O u0 p2 c0  {4,D}
            """
        molecule = Molecule().from_adjacency_list(adjlist)

        assert molecule.multiplicity == 2

        atom1 = molecule.atoms[0]
        atom2 = molecule.atoms[3]
        atom3 = molecule.atoms[4]
        atom4 = molecule.atoms[5]
        assert molecule.has_bond(atom2, atom1)
        assert molecule.has_bond(atom2, atom3)
        assert molecule.has_bond(atom2, atom4)
        assert not molecule.has_bond(atom1, atom3)
        assert not molecule.has_bond(atom1, atom4)
        bond21 = atom2.bonds[atom1]
        bond23 = atom2.bonds[atom3]
        bond24 = atom2.bonds[atom4]

        assert atom1.label == "*1"
        assert atom1.element.symbol == "C"
        assert atom1.radical_electrons == 1
        assert atom1.charge == 0

        assert atom2.label == "*2"
        assert atom2.element.symbol == "N"
        assert atom2.radical_electrons == 0
        assert atom2.charge == 1

        assert atom3.label == ""
        assert atom3.element.symbol == "O"
        assert atom3.radical_electrons == 0
        assert atom3.charge == -1

        assert atom4.label == ""
        assert atom4.element.symbol == "O"
        assert atom4.radical_electrons == 0
        assert atom4.charge == 0

        assert bond21.is_single()
        assert bond23.is_single()
        assert bond24.is_double()

    def test_from_adjacency_list2(self):
        """
        adjlist: Test the Molecule.from_adjacency_list() method 2.
        """
        # molecule 2
        adjlist = """
1 *1 C u1 {2,S} {3,S} {4,S}
2    H u0 {1,S}
3    H u0 {1,S}
4 *2 N u0 p0 c+1 {1,S} {5,S} {6,D}
5    O u0 p3 c-1 {4,S}
6    O u0 p2 {4,D}
            """
        molecule = Molecule().from_adjacency_list(adjlist)

        assert molecule.multiplicity == 2

        atom1 = molecule.atoms[0]
        atom2 = molecule.atoms[3]
        atom3 = molecule.atoms[4]
        atom4 = molecule.atoms[5]
        assert molecule.has_bond(atom2, atom1)
        assert molecule.has_bond(atom2, atom3)
        assert molecule.has_bond(atom2, atom4)
        assert not molecule.has_bond(atom1, atom3)
        assert not molecule.has_bond(atom1, atom4)
        bond21 = atom2.bonds[atom1]
        bond23 = atom2.bonds[atom3]
        bond24 = atom2.bonds[atom4]

        assert atom1.label == "*1"
        assert atom1.element.symbol == "C"
        assert atom1.radical_electrons == 1
        assert atom1.charge == 0

        assert atom2.label == "*2"
        assert atom2.element.symbol == "N"
        assert atom2.radical_electrons == 0
        assert atom2.charge == 1

        assert atom3.label == ""
        assert atom3.element.symbol == "O"
        assert atom3.radical_electrons == 0
        assert atom3.charge == -1

        assert atom4.label == ""
        assert atom4.element.symbol == "O"
        assert atom4.radical_electrons == 0
        assert atom4.charge == 0

        assert bond21.is_single()
        assert bond23.is_single()
        assert bond24.is_double()

    def test_from_adjacency_list3(self):
        """
        adjlist: Test the Molecule.from_adjacency_list() method 3.
        """
        # molecule 3
        adjlist = """
1 *1 C u1 {2,S} {3,S} {4,S}
2    H u0 {1,S}
3    H u0 {1,S}
4 *2 N u0 p0 c+1 {1,S} {5,S} {6,D}
5    O u0 p3 c-1 {4,S}
6    O u0 p2 {4,D}
            """
        molecule = Molecule().from_adjacency_list(adjlist)

        assert molecule.multiplicity == 2

        atom1 = molecule.atoms[0]
        atom2 = molecule.atoms[3]
        atom3 = molecule.atoms[4]
        atom4 = molecule.atoms[5]
        assert molecule.has_bond(atom2, atom1)
        assert molecule.has_bond(atom2, atom3)
        assert molecule.has_bond(atom2, atom4)
        assert not molecule.has_bond(atom1, atom3)
        assert not molecule.has_bond(atom1, atom4)
        bond21 = atom2.bonds[atom1]
        bond23 = atom2.bonds[atom3]
        bond24 = atom2.bonds[atom4]

        assert atom1.label == "*1"
        assert atom1.element.symbol == "C"
        assert atom1.radical_electrons == 1
        assert atom1.charge == 0

        assert atom2.label == "*2"
        assert atom2.element.symbol == "N"
        assert atom2.radical_electrons == 0
        assert atom2.charge == 1

        assert atom3.label == ""
        assert atom3.element.symbol == "O"
        assert atom3.radical_electrons == 0
        assert atom3.charge == -1

        assert atom4.label == ""
        assert atom4.element.symbol == "O"
        assert atom4.radical_electrons == 0
        assert atom4.charge == 0

        assert bond21.is_single()
        assert bond23.is_single()
        assert bond24.is_double()

    def test_from_adjacency_list4(self):
        """
        adjlist: Test the Molecule.from_adjacency_list() method 4.
        """
        # molecule 4
        adjlist = """
1 *1 C u1 {2,S}
2 *2 N u0 p0 c+1 {1,S} {3,S} {4,D}
3    O u0 p3 c-1 {2,S}
4    O u0 p2 {2,D}
            """
        molecule = Molecule().from_adjacency_list(adjlist, saturate_h=True)

        assert molecule.multiplicity == 2

        atom1 = molecule.atoms[0]
        atom2 = molecule.atoms[1]
        atom3 = molecule.atoms[2]
        atom4 = molecule.atoms[3]
        assert molecule.has_bond(atom2, atom1)
        assert molecule.has_bond(atom2, atom3)
        assert molecule.has_bond(atom2, atom4)
        assert not molecule.has_bond(atom1, atom3)
        assert not molecule.has_bond(atom1, atom4)
        bond21 = atom2.bonds[atom1]
        bond23 = atom2.bonds[atom3]
        bond24 = atom2.bonds[atom4]

        assert atom1.label == "*1"
        assert atom1.element.symbol == "C"
        assert atom1.radical_electrons == 1
        assert atom1.charge == 0

        assert atom2.label == "*2"
        assert atom2.element.symbol == "N"
        assert atom2.radical_electrons == 0
        assert atom2.charge == 1

        assert atom3.label == ""
        assert atom3.element.symbol == "O"
        assert atom3.radical_electrons == 0
        assert atom3.charge == -1

        assert atom4.label == ""
        assert atom4.element.symbol == "O"
        assert atom4.radical_electrons == 0
        assert atom4.charge == 0

        assert bond21.is_single()
        assert bond23.is_single()
        assert bond24.is_double()

    def test_from_adjacency_list5(self):
        """
        adjlist: Test if from_adjacency_list works when saturateH is turned on
        and test molecule is fused aromatics.
        """
        # molecule 5
        adjlist = """
1  * C u0 p0 c0 {2,B} {3,B} {4,B}
2    C u0 p0 c0 {1,B} {5,B} {6,B}
3    C u0 p0 c0 {1,B} {8,B} {13,S}
4    C u0 p0 c0 {1,B} {9,B}
5    C u0 p0 c0 {2,B} {10,B}
6    C u0 p0 c0 {2,B} {7,B}
7    C u0 p0 c0 {6,B} {8,B} {11,S}
8    C u0 p0 c0 {3,B} {7,B} {12,S}
9    C u0 p0 c0 {4,B} {10,B}
10   C u0 p0 c0 {5,B} {9,B}
11   H u0 p0 c0 {7,S}
12   H u0 p0 c0 {8,S}
13   H u0 p0 c0 {3,S}
            """
        molecule = Molecule().from_adjacency_list(adjlist, saturate_h=True)

        assert molecule.multiplicity == 1

        atom1 = molecule.atoms[0]
        atom2 = molecule.atoms[1]
        atom3 = molecule.atoms[2]
        atom7 = molecule.atoms[6]
        atom11 = molecule.atoms[10]
        bond21 = atom2.bonds[atom1]
        bond13 = atom1.bonds[atom3]
        bond7_11 = atom7.bonds[atom11]

        assert atom1.label == "*"
        assert atom1.element.symbol == "C"
        assert atom1.radical_electrons == 0
        assert atom1.charge == 0

        assert atom2.label == ""
        assert atom2.element.symbol == "C"
        assert atom2.radical_electrons == 0
        assert atom2.charge == 0

        assert bond21.is_benzene()
        assert bond13.is_benzene()
        assert bond7_11.is_single()

    def test_various_spin_adjlists(self):
        """
        adjlist: Test that molecules with old or intermediate adjacency list formats containing unusual
        spin states can get converted to the proper new adjlist format.
        """

        adjlist_2s = """
1 C 2S 0 {2,S} {3,S}
2 H 0  0 {1,S}
3 H 0  0 {1,S}
"""
        adjlist_2s_new = """
1 C u0 p1 c0  {2,S} {3,S}
2 H u0 p0 c0  {1,S}
3 H u0 p0 c0  {1,S}
"""
        mol_2s = Molecule().from_adjacency_list(adjlist_2s)
        mol_2s_new = Molecule().from_adjacency_list(adjlist_2s_new)
        assert mol_2s.is_isomorphic(mol_2s_new)

        adjlist_2t = """
1 C 2T 0 {2,S} {3,S}
2 H 0  0 {1,S}
3 H 0  0 {1,S}
"""
        adjlist_2t_new = """
1 C u2 p0 c0  {2,S} {3,S}
2 H u0 p0 c0  {1,S}
3 H u0 p0 c0  {1,S}
"""
        mol_2t = Molecule().from_adjacency_list(adjlist_2t)
        mol_2t_new = Molecule().from_adjacency_list(adjlist_2t_new)
        assert mol_2t.is_isomorphic(mol_2t_new)

        adjlist_3d = """
1 C 3D 0 {2,S}
2 H 0  0 {1,S}
"""
        adjlist_3d_new = """
1 C u1 p1 c0  {2,S}
2 H u0 p0 c0  {1,S}
"""
        mol_3d = Molecule().from_adjacency_list(adjlist_3d)
        mol_3d_new = Molecule().from_adjacency_list(adjlist_3d_new)
        assert mol_3d.is_isomorphic(mol_3d_new)

        adjlist_3q = """
1 N 3Q 1
"""
        adjlist_3q_new = """
1 N u3 p1 c0 
"""
        mol_3q = Molecule().from_adjacency_list(adjlist_3q)
        mol_3q_new = Molecule().from_adjacency_list(adjlist_3q_new)
        assert mol_3q.is_isomorphic(mol_3q_new)

        adjlist_4s = """
1 C 4S 0
        """
        adjlist_4s_new = """
1 C u0 p2 c0
"""
        mol_4s = Molecule().from_adjacency_list(adjlist_4s)
        mol_4s_new = Molecule().from_adjacency_list(adjlist_4s_new)
        assert mol_4s.is_isomorphic(mol_4s_new)

        adjlist_4t = """
1 C 4T 0
"""
        adjlist_4t_new = """
1 C u2 p1 c0
"""
        mol_4t = Molecule().from_adjacency_list(adjlist_4t)
        mol_4t_new = Molecule().from_adjacency_list(adjlist_4t_new)
        assert mol_4t.is_isomorphic(mol_4t_new)

        adjlist_4v = """
1 C 4V 0
"""
        adjlist_4v_new = """
1 C u4 p0 c0        
"""
        mol_4v = Molecule().from_adjacency_list(adjlist_4v)
        mol_4v_new = Molecule().from_adjacency_list(adjlist_4v_new)
        assert mol_4v.is_isomorphic(mol_4v_new)

    def test_wildcard_adjlists(self):
        """
        adjlist: Test that molecule adjlists containing wildcards raise an InvalidAdjacencyListError.
        """
        # A molecule with a wildcard assignment
        wildcard_adjlist1 = "1 C u1 px c0"
        wildcard_adjlist2 = "1 C ux p2 c0"
        wildcard_adjlist3 = "1 C u1 p2 cx"
        wildcard_adjlist4 = "1 [C,N] u1 p2 c0"

        with pytest.raises(InvalidAdjacencyListError):
            Molecule().from_adjacency_list(wildcard_adjlist1)
        with pytest.raises(InvalidAdjacencyListError):
            Molecule().from_adjacency_list(wildcard_adjlist2)
        with pytest.raises(InvalidAdjacencyListError):
            Molecule().from_adjacency_list(wildcard_adjlist3)
        with pytest.raises(InvalidAdjacencyListError):
            Molecule().from_adjacency_list(wildcard_adjlist4)

    def test_incorrect_adjlists(self):
        """
        adjlist: Test that improperly formed adjlists raise an InvalidAdjacencyListError.
        """
        # Carbon with 1 radical and 2 lone pairs = 5 total electrons.  Should have -1 charge but doesn't
        adjlist1 = "1 C u1 p2 c0"

        with pytest.raises(InvalidAdjacencyListError):
            Molecule().from_adjacency_list(adjlist1)

    def test_helium(self):
        """
        adjlist: Test that the adjlist reading and writing works with Helium.
        """
        smiles = "[He]"
        inchi = "InChI=1S/He"
        adjlist = "1 He u0 p1 c0"
        adjlist_old = "1 He 0"
        adjlist_intermediate = "1 He 0 1"

        mol_smiles = Molecule().from_smiles(smiles)
        mol_inchi = Molecule().from_inchi(inchi)
        mol = Molecule().from_adjacency_list(adjlist)
        mol_old = Molecule().from_adjacency_list(adjlist_old)
        mol_intermediate = Molecule().from_adjacency_list(adjlist_intermediate)

        # Isomorphic check
        assert mol_smiles.is_isomorphic(mol)
        assert mol_smiles.is_isomorphic(mol_inchi)
        assert mol_smiles.is_isomorphic(mol_old)
        assert mol_smiles.is_isomorphic(mol_intermediate)

        # Adjlist check
        assert mol_smiles.to_adjacency_list().strip() == adjlist
        assert mol_inchi.to_adjacency_list().strip() == adjlist
        assert mol.to_adjacency_list().strip() == adjlist
        assert mol_old.to_adjacency_list().strip() == adjlist
        assert mol_intermediate.to_adjacency_list().strip() == adjlist

        assert mol.to_smiles() == smiles
        assert mol.to_inchi() == "InChI=1S/He"

    def test_to_adjacency_list(self):
        """
        adjlist: Test the Molecule.to_adjacency_list() method.
        """
        inter_adjlist = """
1 *1 C 1 0  {2,S} {3,S} {4,S}
2    H 0 0  {1,S}
3    H 0 0  {1,S}
4 *2 N 0 0 {1,S} {5,S} {6,D}
5    O 0 3 {4,S}
6    O 0 2 {4,D}
            """

        adjlist = """
1 *1 C u1 p0 c0  {2,S} {3,S} {4,S}
2    H u0 p0 c0  {1,S}
3    H u0 p0 c0  {1,S}
4 *2 N u0 p0 c+1 {1,S} {5,S} {6,D}
5    O u0 p3 c-1 {4,S}
6    O u0 p2 c0  {4,D}
            """
        molecule = Molecule().from_adjacency_list(adjlist)
        molecule2 = Molecule().from_adjacency_list(inter_adjlist)
        adjlist_1 = molecule.to_adjacency_list(remove_h=False)
        assert adjlist_1 == molecule2.to_adjacency_list()
        new_molecule = Molecule().from_adjacency_list(adjlist_1)
        assert molecule.is_isomorphic(new_molecule)

    def test_to_adjacency_list_for_non_integer_bonds(self):
        """
        Test the adjacency list can be created for molecules with bond orders
        that don't fit into single, double, triple, or benzene
        """
        from rmgpy.molecule.molecule import Atom, Bond, Molecule

        atom1 = Atom(element="H", lone_pairs=0)
        atom2 = Atom(element="H", lone_pairs=0)
        bond = Bond(atom1, atom2, 0.5)
        mol = Molecule(multiplicity=1)
        mol.add_atom(atom1)
        mol.add_atom(atom2)
        mol.add_bond(bond)
        adjlist = mol.to_adjacency_list()
        assert "H" in adjlist
        assert "{1,0.5}" in adjlist

    @pytest.mark.skip(reason="WIP")
    def test_from_adjacency_list_for_non_integer_bonds(self):
        """
        Test molecule can be created from the adjacency list for molecules with bond orders
        that don't fit into single, double, triple, or benzene.

        This test is a work in progress since currently reading one of these
        objects thows an `InvalidAdjacencyListError`. Since the number radical
        electrons is an integer, having fractional bonds leads to this error.
        Fixing it would require switching radical electrons to floats.
        """
        from rmgpy.molecule.molecule import Molecule

        adjlist = """
        1 H u1 p2 c0 {2,0.5}
        2 H u1 p2 c0 {1,0.5}
        """
        mol = Molecule().from_adjacency_list(adjlist)
        atom0 = mol.atoms[0]
        assert len(atom0.bonds) == 1
        bond = next(iter(atom0.bonds.values()))
        assert round(abs(bond[0].get_order_num() - 0.5), 7) == 0

    def test_from_intermediate_adjacency_list1(self):
        """
        Test we can read an intermediate style adjacency list with implicit hydrogens 1
        """
        adjlist = """
        1 O 0 2
        """  # should be Water
        molecule = Molecule().from_adjacency_list(adjlist, saturate_h=True)
        assert molecule.get_formula() == "H2O"

    def test_from_old_adjacency_list1(self):
        """
        Test we can read an old style adjacency list with implicit hydrogens 1
        """
        adjlist = """
        1 O 0 
        """  # should be Water
        molecule = Molecule().from_adjacency_list(adjlist)
        assert molecule.get_formula() == "H2O"

    def test_from_old_adjacency_list2(self):
        """
        Test we can read an old style adjacency list with implicit hydrogens 2
        """
        adjlist = """
        1 C 2S
        """
        adjlist_new = """
        1 C u0 p1 c0 {2,S} {3,S}
        2 H u0 p0 c0 {1,S}
        3 H u0 p0 c0 {1,S}
        """
        molecule = Molecule().from_adjacency_list(adjlist)
        molecule_new = Molecule().from_adjacency_list(adjlist_new)
        assert molecule.is_isomorphic(molecule_new)

    def test_from_old_adjacency_list3(self):
        """
        Test we can read an old style adjacency list with implicit hydrogens 3
        """
        adjlist = """
        1 C 0
        """
        adjlist_new = """
        1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
        2 H u0 p0 c0 {1,S}
        3 H u0 p0 c0 {1,S}
        4 H u0 p0 c0 {1,S}
        5 H u0 p0 c0 {1,S}
        """
        molecule = Molecule().from_adjacency_list(adjlist)
        molecule_new = Molecule().from_adjacency_list(adjlist_new)
        assert molecule.is_isomorphic(molecule_new)

    def test_from_old_adjacency_list4(self):
        """
        Test we can read an old style adjacency list with implicit hydrogens 4
        """
        adjlist = """
        1 O 2S
        """
        adjlist_new = """
        1 O u0 p3 c0
        """
        molecule = Molecule().from_adjacency_list(adjlist)
        molecule_new = Molecule().from_adjacency_list(adjlist_new)
        assert molecule.is_isomorphic(molecule_new)

    @pytest.mark.skip(reason="WIP")
    def test_from_old_adjacency_list5(self):
        """
        Test we can read an old style adjacency list with implicit hydrogens 5
        """
        adjlist = """
        1 C 2S {2,T}
        2 O 2S {1,T}
        """
        adjlist_new = """
        1 C u0 p1 c-1 {2,T}
        2 O u0 p1 c+1 {1,T}
        """
        molecule = Molecule().from_adjacency_list(adjlist)
        molecule_new = Molecule().from_adjacency_list(adjlist_new)
        assert molecule.is_isomorphic(molecule_new)
        # Currently the from_old_adjacency_list cannot correctly interpret CO written in this old form
        # (I don't think any adjlists are actually formed this way.)
        # Currently 'adjlist' will fail when the Molecule is determined to be non-neurtral in net charge.

    def test_from_old_adjacency_list6(self):
        """
        Test we can read an old style adjacency list with implicit hydrogens 1
        """
        adjlist = """
        1 C 4T
        """
        adjlist_new = """
        1 C u2 p1 c0 
        """
        molecule = Molecule().from_adjacency_list(adjlist)
        molecule_new = Molecule().from_adjacency_list(adjlist_new)
        assert molecule.is_isomorphic(molecule_new)

    def test_adjacency_list(self):
        """
        adjlist: Check the adjacency list read/write functions for a full molecule.
        """
        molecule1 = Molecule().from_adjacency_list(
            """
        1  C u0 {2,D} {7,S} {8,S}
        2  C u0 {1,D} {3,S} {9,S}
        3  C u0 {2,S} {4,D} {10,S}
        4  C u0 {3,D} {5,S} {11,S}
        5  C u1 {4,S} {6,S} {12,S}
        6  C u0 {5,S} {13,S} {14,S} {15,S}
        7  H u0 {1,S}
        8  H u0 {1,S}
        9  H u0 {2,S}
        10 H u0 {3,S}
        11 H u0 {4,S}
        12 H u0 {5,S}
        13 H u0 {6,S}
        14 H u0 {6,S}
        15 H u0 {6,S}
        """
        )
        molecule2 = Molecule().from_smiles("C=CC=C[CH]C")
        assert molecule1.is_isomorphic(molecule2)
        assert molecule2.is_isomorphic(molecule1)

        # Test that charges are correctly stored and written with adjacency lists
        adjlist3 = """
1 C u0 p1 c-1 {2,T}
2 O u0 p1 c+1 {1,T}
"""
        molecule3 = Molecule().from_adjacency_list(adjlist3)
        assert molecule3.atoms[0].charge == -1
        assert molecule3.atoms[1].charge == 1
        adjlist4 = molecule3.to_adjacency_list()
        assert adjlist3.strip() == adjlist4.strip()

    def test_group_adjacency_list(self):
        """
        adjlist: Check the adjacency list read/write functions for a full molecule.
        """
        adjlist = """1 C u0 {2,D}
2 O u1 p1 c[-1,0,+1] {1,D}
"""
        group = Group().from_adjacency_list(
            """
        1 C u0 {2,D} 
        2 O u1 p1 c[-1,0,+1] {1,D}
        """
        )
        assert adjlist == group.to_adjacency_list()

    def test_to_old_ajacency_list(self):
        """
        adjlist: Check that we can convert back to old style adjacency list
        """
        molecule2 = Molecule().from_smiles("C=CC=C[CH]C")
        string = """1 C 0 {2,D}
2 C 0 {1,D} {3,S}
3 C 0 {2,S} {4,D}
4 C 0 {3,D} {5,S}
5 C 1 {4,S} {6,S}
6 C 0 {5,S}"""
        assert molecule2.to_adjacency_list(remove_h=True, old_style=True).strip() == string.strip()


class TestConsistencyChecker:
    def test_check_hund_rule_fail(self):
        with pytest.raises(InvalidAdjacencyListError):
            Molecule().from_adjacency_list(
                """
            multiplicity 1
            1 C u2 p0 c0
            """,
                saturate_h=True,
            )

    def test_check_hund_rule_success(self):
        try:
            Molecule().from_adjacency_list(
                """
            multiplicity 3
            1 C u2 p0 c0
            """,
                saturate_h=True,
            )
        except InvalidAdjacencyListError:
            self.fail("InvalidAdjacencyListError thrown unexpectedly!")

    def test_check_multiplicity(self):
        """
        adjlist: Check that RMG allows different electron spins in the same molecule with multiplicity = 2s + 1
        """
        # [N] radical:
        try:
            Molecule().from_adjacency_list(
                """multiplicity 4
                                            1 N u3 p1 c0"""
            )
        except InvalidAdjacencyListError:
            self.fail("InvalidAdjacencyListError thrown unexpectedly for N tri-rad!")

        # A general molecule with 4 radicals, multiplicity 5:
        try:
            Molecule().from_adjacency_list(
                """multiplicity 5
                                            1 O u1 p2 c0 {2,S}
                                            2 C u1 p0 c0 {1,S} {3,S} {4,S}
                                            3 H u0 p0 c0 {2,S}
                                            4 N u1 p1 c0 {2,S} {5,S}
                                            5 O u1 p2 c0 {4,S}"""
            )
        except InvalidAdjacencyListError:
            self.fail("InvalidAdjacencyListError thrown unexpectedly for a molecule with 4 radicals, multiplicity 5")

        # A general molecule with 4 radicals, multiplicity 3:
        try:
            Molecule().from_adjacency_list(
                """multiplicity 3
                                            1 O u1 p2 c0 {2,S}
                                            2 C u1 p0 c0 {1,S} {3,S} {4,S}
                                            3 H u0 p0 c0 {2,S}
                                            4 N u1 p1 c0 {2,S} {5,S}
                                            5 O u1 p2 c0 {4,S}"""
            )
        except InvalidAdjacencyListError:
            self.fail("InvalidAdjacencyListError thrown unexpectedly for a molecule with 4 radicals, multiplicity 3")

        # [N]=C=[N] singlet:
        try:
            Molecule().from_adjacency_list(
                """multiplicity 1
                                            1 N u1 p1 c0 {2,D}
                                            2 C u0 p0 c0 {1,D} {3,D}
                                            3 N u1 p1 c0 {2,D}"""
            )
        except InvalidAdjacencyListError:
            self.fail("InvalidAdjacencyListError thrown unexpectedly for singlet [N]=C=[N]!")
