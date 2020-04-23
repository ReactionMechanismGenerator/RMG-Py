#!/usr/bin/env python
# -*- coding: utf-8 -*-

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

"""
This script contains unit tests of the :mod:`arkane.reference` module.
"""

import os
import unittest

from arkane.encorr.reference import ReferenceSpecies
from rmgpy.species import Species


################################################################################
FILE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(os.path.dirname(FILE_DIR), 'data')
SPECIES_DIR = os.path.join(DATA_DIR, 'species')


class TestReferenceSpecies(unittest.TestCase):
    """
    A class for testing that the ReferenceSpecies class functions properly
    """

    @classmethod
    def setUpClass(cls):
        cls.methane = Species(smiles='C')
        cls.ethane = Species(smiles='CC')
        cls.propane = Species(smiles='CCC')

        cls.thermo_data = ThermoData(H298=(100.0, 'kJ/mol'), S298=(100.0, 'J/(mol*K)'))
        cls.xyz_dict = {'symbols': ('H', 'H'), 'isotopes': (1, 1), 'coords': ((0.0, 0.0, 0.0), (0.708, 0.0, 0.0))}
        cls.unscaled_freqs = ArrayQuantity([3000.0], 'cm^-1')
        cls.electronic_energy = ScalarQuantity(10.0, 'J/mol')
        cls.t1_diagnostic = 0.101

    def test_instantiate_reference_species(self):
        """
        Test that a ReferenceSpecies object can be instantiated with the minimal acceptable input, and throws an error
        if the minimal acceptable input is not given.
        """
        ref_spcs = ReferenceSpecies(species=self.ethane)
        self.assertEqual(ref_spcs.smiles, 'CC')
        self.assertEqual(ref_spcs.inchi_key, self.ethane.molecule[0].to_inchi_key())

        ref_from_smiles = ReferenceSpecies(smiles='CCC')
        self.assertEqual(ref_from_smiles.smiles, 'CCC')
        self.assertEqual(ref_from_smiles.inchi_key, self.propane.molecule[0].to_inchi_key())

        ref_from_inchi = ReferenceSpecies(inchi='InChI=1S/C2H6/c1-2/h1-2H3')
        self.assertEqual(ref_from_inchi.smiles, 'CC')
        self.assertEqual(ref_from_inchi.inchi_key, self.ethane.molecule[0].to_inchi_key())

        ref_from_adj = ReferenceSpecies(adjacency_list='1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}\n2 H u0 p0 c0 {1,S}\n'
                                                       '3 H u0 p0 c0 {1,S}\n4 H u0 p0 c0 {1,S}\n5 H u0 p0 c0 {1,S}\n')
        self.assertEqual(ref_from_adj.smiles, 'C')
        self.assertEqual(ref_from_adj.inchi_key, self.methane.molecule[0].to_inchi_key())

        with self.assertRaises(ValueError):
            ReferenceSpecies()

    def load_ref_from_yaml(self):
        """
        Test that the example ReferenceSpecies YAML file can be loaded
        """
        ref_spcs = ReferenceSpecies.__new__(ReferenceSpecies)
        ref_spcs.load_yaml(os.path.join(SPECIES_DIR, 'reference_species_example.yml'))

        self.assertEqual(ref_spcs.smiles, 'C#C[CH2]')
        self.assertEqual(ref_spcs.label, 'example_reference_species')
        self.assertIsInstance(ref_spcs.calculated_data, CalculatedDataEntry)

    def test_save_ref_to_yaml(self):
        """
        Test that a ReferenceSpecies object can be saved to a YAML file successfully
        """
        label = 'test_reference_species'
        ref_spcs = ReferenceSpecies(species=self.ethane, label=label)
        self.assertEqual(ref_spcs.label, label)
        ref_spcs.save_yaml(path=SPECIES_DIR)

        loaded_ref = ReferenceSpecies.__new__(ReferenceSpecies)
        load_path = os.path.join(SPECIES_DIR, f'{label}.yml')
        loaded_ref.load_yaml(path=load_path)

        self.assertEqual(loaded_ref.smiles, 'CC')

        # Finally, delete this newly created file
        os.remove(load_path)


if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
