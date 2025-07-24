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

"""
This module contains unit tests of the rmgpy.molecule.atomtype module.
"""

import os
import os.path
import itertools

from rmgpy.molecule import Molecule, Atom, Bond
from rmgpy.molecule.draw import MoleculeDrawer
from rmgpy.species import Species


class TestMoleculeDrawer:
    """
    Contains unit tests of the MoleculeDrawer class.
    """

    def setup_class(self):
        self.drawer = MoleculeDrawer()
        self.molecule = Molecule(smiles="CC(=O)CC")

    def test_draw_png(self):
        """
        Test we can create PNG files from molecules.
        """
        try:
            from cairocffi import ImageSurface
        except ImportError:
            from cairo import ImageSurface
        path = "test_molecule.png"
        if os.path.exists(path):
            os.unlink(path)
        surface, _cr, (_xoff, _yoff, width, height) = self.drawer.draw(self.molecule, file_format="png", target=path)
        assert os.path.exists(path), "File doesn't exist"
        assert width > height
        os.unlink(path)
        assert isinstance(surface, ImageSurface)

    def test_draw_pdf(self):
        """
        Test we can create PDF files from molecules.
        """
        try:
            from cairocffi import PDFSurface
        except ImportError:
            from cairo import PDFSurface
        path = "test_molecule.pdf"
        if os.path.exists(path):
            os.unlink(path)
        surface, _cr, (_xoff, _yoff, width, height) = self.drawer.draw(self.molecule, file_format="pdf", target=path)
        assert isinstance(surface, PDFSurface)
        assert width > height
        os.unlink(path)

    def test_draw_polycycle(self):
        """
        Test we can draw a polycyclic molecule
        """
        try:
            from cairocffi import PDFSurface
        except ImportError:
            from cairo import PDFSurface
        path = "test_molecule.pdf"
        if os.path.exists(path):
            os.unlink(path)
        polycycle = Molecule(smiles="C123CC4CC1COCC2CCC34")
        surface, _cr, (_xoff, _yoff, width, height) = self.drawer.draw(polycycle, file_format="pdf", target=path)
        assert isinstance(surface, PDFSurface)
        assert width > height
        os.unlink(path)

    def test_draw_pdf_without_file(self):
        """
        Test we can create PDF surface without a temporary file (newer versions of PyCairo?)
        """
        try:
            from cairocffi import PDFSurface
        except ImportError:
            from cairo import PDFSurface
        surface, _cr, (_xoff, _yoff, width, height) = self.drawer.draw(self.molecule, file_format="pdf")
        assert isinstance(surface, PDFSurface)
        assert width > height

    def test_draw_non_standard_bonds(self):
        spec = Species().from_smiles("[CH2]C=C[CH2]")
        hybrid = spec.get_resonance_hybrid()
        try:
            from cairocffi import PDFSurface
        except ImportError:
            from cairo import PDFSurface
        surface, _cr, (_xoff, _yoff, width, height) = self.drawer.draw(hybrid, file_format="pdf")
        assert width > height
        assert isinstance(surface, PDFSurface)

    def test_draw_hydrogen_bond_adsorbate(self):
        molecule = Molecule().from_adjacency_list(
            """
1  O u0 p3 c-1 {2,S} {10,H}
2  N u0 p0 c+1 {1,S} {3,D} {4,S}
3  O u0 p2 c0 {2,D}
4  O u0 p2 c0 {2,S} {7,S}
5  N u0 p1 c0 {6,S} {8,S} {9,S} {7,H}
6  O u0 p2 c0 {5,S} {10,S}
7  H u0 p0 c0 {4,S} {5,H}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S} {1,H}
11 X u0 p0 c0
        """
        )
        try:
            from cairocffi import PDFSurface
        except ImportError:
            from cairo import PDFSurface
        surface, _cr, (_xoff, _yoff, _width, _height) = self.drawer.draw(molecule, file_format="pdf")
        assert isinstance(surface, PDFSurface)

    def test_draw_bidentate_adsorbate(self):
        try:
            from cairocffi import ImageSurface
        except ImportError:
            from cairo import ImageSurface

        test_molecules = [
            Molecule().from_adjacency_list(
            """
1  C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
4  X u0 p0 c0 {1,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  X u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
        """),
            Molecule().from_adjacency_list(
"""
1  C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
4  X u0 p0 c0 {1,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 X u0 p0 c0 {3,S}
"""),
        ]
        for molecule in [Molecule(smiles="CC(=O)CCO"), Molecule(smiles="C1CC=CC1CC"), Molecule(smiles="C=CCC(O)=O")]:
            bondable = [a for a in molecule.atoms if a.is_non_hydrogen() and any(b.is_hydrogen() for b in a.bonds)]
            for b1, b2 in itertools.combinations(bondable, 2):
                # find a hydrogen atom bonded to each of the two atoms
                for h1 in b1.bonds:
                    if h1.is_hydrogen():
                        break
                for h2 in b2.bonds:
                    if h2.is_hydrogen():
                        break
                molecule.remove_atom(h1)
                molecule.remove_atom(h2)
                x1 = Atom(element='X', radical_electrons=0, charge=0, label='', lone_pairs=0)
                x2 = Atom(element='X', radical_electrons=0, charge=0, label='', lone_pairs=0)
                molecule.add_atom(x1)
                molecule.add_atom(x2)
                molecule.add_bond(Bond(b1, x1, order=1))
                molecule.add_bond(Bond(b2, x2, order=1))
                test_molecules.append(molecule.copy(deep=True))
                molecule.remove_atom(x1)
                molecule.remove_atom(x2)
                molecule.add_atom(h1)
                molecule.add_atom(h2)
                molecule.add_bond(Bond(b1, h1, order=1))
                molecule.add_bond(Bond(b2, h2, order=1))

            for b1, b2, b3 in itertools.combinations(bondable, 3):
                # find a hydrogen atom bonded to each of the two atoms
                for h1 in b1.bonds:
                    if h1.is_hydrogen():
                        break
                for h2 in b2.bonds:
                    if h2.is_hydrogen():
                        break
                for h3 in b3.bonds:
                    if h3.is_hydrogen():
                        break
                molecule.remove_atom(h1)
                molecule.remove_atom(h2)
                molecule.remove_atom(h3)
                x1 = Atom(element='X', radical_electrons=0, charge=0, label='', lone_pairs=0)
                x2 = Atom(element='X', radical_electrons=0, charge=0, label='', lone_pairs=0)
                x3 = Atom(element='X', radical_electrons=0, charge=0, label='', lone_pairs=0)
                molecule.add_atom(x1)
                molecule.add_atom(x2)
                molecule.add_atom(x3)
                molecule.add_bond(Bond(b1, x1, order=1))
                molecule.add_bond(Bond(b2, x2, order=1))
                molecule.add_bond(Bond(b3, x3, order=1))
                test_molecules.append(molecule.copy(deep=True))
                molecule.remove_atom(x1)
                molecule.remove_atom(x2)
                molecule.remove_atom(x3)
                molecule.add_atom(h1)
                molecule.add_atom(h2)
                molecule.add_atom(h3)
                molecule.add_bond(Bond(b1, h1, order=1))
                molecule.add_bond(Bond(b2, h2, order=1))
                molecule.add_bond(Bond(b3, h3, order=1))
            
        test_molecules.append(Molecule().from_adjacency_list(
"""
1  C  u0 p0 c0 {2,S} {3,S} {9,S} {10,S}
2  X u0 p0 c0 {1,S}
3  C  u0 p0 c0 {1,S} {4,S} {5,S} {11,S}
4  X u0 p0 c0 {3,S}
5  C  u0 p0 c0 {3,S} {6,S} {7,S} {12,S}
6  X u0 p0 c0 {5,S}
7  C  u0 p0 c0 {5,S} {8,S} {13,S} {14,S}
8  X u0 p0 c0 {7,S}
9  H  u0 p0 c0 {1,S}
10 H  u0 p0 c0 {1,S}
11 H  u0 p0 c0 {3,S}
12 H  u0 p0 c0 {5,S}
13 H  u0 p0 c0 {7,S}
14 H  u0 p0 c0 {7,S}
"""))
        test_molecules.append(Molecule().from_adjacency_list(
"""
1 O u0 p2 c0 {4,S} {5,S}
2 O u0 p2 c0 {5,S} {6,S}
3 O u0 p2 c0 {6,S} {26,S}
4 C u0 p0 c0 {1,S} {7,S} {8,S} {19,S}
5 C u0 p0 c0 {1,S} {2,S} {10,S} {21,S}
6 C u0 p0 c0 {2,S} {3,S} {9,S} {20,S}
7 C u0 p0 c0 {4,S} {15,S} {16,S} {24,S}
8 C u0 p0 c0 {4,S} {17,S} {18,S} {25,S}
9 C u0 p0 c0 {6,S} {11,S} {12,S} {22,S}
10 C u0 p0 c0 {5,S} {13,S} {14,S} {23,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {8,S}
19 X u0 p0 c0 {4,S}
20 X u0 p0 c0 {6,S}
21 X u0 p0 c0 {5,S}
22 X u0 p0 c0 {9,S}
23 X u0 p0 c0 {10,S}
24 X u0 p0 c0 {7,S}
25 X u0 p0 c0 {8,S}
26 X u0 p0 c0 {3,S}
"""))
        test_molecules.append(Molecule(smiles="*CC(*)(C*)OCC#*"))
        test_molecules.append(Molecule(smiles="*CC(*)(C*)C*"))
        for number, molecule in enumerate(test_molecules, 1):
            path = f"test_polydentate_{number}.png"
            if os.path.exists(path):
                os.unlink(path)
            self.drawer.clear()
            surface, _cr, (_xoff, _yoff, width, height) = self.drawer.draw(molecule, file_format="png", target=path)
            assert os.path.exists(path), "File doesn't exist"
            os.unlink(path)
            assert isinstance(surface, ImageSurface)

    def test_draw_bidentate_with_charge_separation(self):
        molecule = Molecule().from_adjacency_list(
            """
1 X u0 p0 c0 {3,S}
2 X u0 p0 c0 {4,D}
3 O u0 p2 c0 {1,S} {4,S}
4 N u0 p0 c+1 {3,S} {2,D} {5,S}
5 O u0 p3 c-1 {4,S}
        """
        )
        try:
            from cairocffi import PDFSurface
        except ImportError:
            from cairo import PDFSurface
        surface, _cr, (_xoff, _yoff, _width, _height) = self.drawer.draw(molecule, file_format="pdf")
        assert isinstance(surface, PDFSurface)
