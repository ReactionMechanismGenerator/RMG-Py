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

import os
import shutil
import unittest

from rdkit import Chem
from rmgpy.molecule import Molecule
from rmgpy.qm.main import QMSettings
from rmgpy.qm.molecule import Geometry

scratch_dir = os.path.join(os.path.abspath(os.path.dirname(__file__)), "scratch")


class TestQMMolecule(unittest.TestCase):
    """
    Contains unit tests for the Geometry class in the qm module.
    """

    def test_rd_embed_common_species(self):
        """
        Test that rd_embed() works for various commmon molecules.
        """
        for smi in [
            "O",
            "C",
            "[CH2]",
            "[CH3]",
            "[OH]",
            "O=C=O",
            "[CH3+]",
            "[OH-]",
            "N",
            "N#N",
            "[NH2]",
            "CO[O]",
        ]:
            mol = Molecule().from_smiles(smi)
            geom = Geometry(
                settings=QMSettings(
                    scratchDirectory=scratch_dir,
                ),
                unique_id="0",
                molecule=mol,
            )
            rdmol, _ = geom.rd_build()
            rdmol, _ = geom.rd_embed(rdmol, num_conf_attempts=20)
            self.assertEqual(rdmol.GetNumConformers(), 20)
        shutil.rmtree(scratch_dir)

    def test_rd_embed_uncommon_species(self):
        """
        Test that rd_embed() works for species that ETKDG fails to embed conformers for.
        """
        for smi in [
            "[N:1](=[C-:2][C@@:3]1([H:8])[O:4][C@@:5]2([H:9])[C:6]([H:10])([H:11])[C+:7]12)[H:12]",
            "[N+:1]1#[C:2][C@@:3]2([H:8])[O:4][C-:5]([H:9])[C:6]([H:10])([H:11])[C@@:7]12[H:12]",
        ]:
            # RMG has difficulty converting these charged species into RDKit Mol.
            # Directly, generating rdmol using RDKit here.
            # So far, this shouldn't be a concern since RMG doesn't work for charged species
            # And these species are only to test the `rd_embed`` method
            geom = Geometry(
                settings=QMSettings(
                    scratchDirectory=scratch_dir,
                ),
                unique_id="0",
                molecule=None,
            )
            rdmol = Chem.rdmolops.AddHs(Chem.MolFromSmiles(smi))
            rdmol, _ = geom.rd_embed(rdmol, num_conf_attempts=20)
            self.assertEqual(rdmol.GetNumConformers(), 20)
        shutil.rmtree(scratch_dir)
