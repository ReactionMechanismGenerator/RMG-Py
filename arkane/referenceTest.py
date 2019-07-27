#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2019 Prof. William H. Green (whgreen@mit.edu),           #
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
import shutil

from arkane.reference import ReferenceSpecies, ReferenceDataEntry, CalculatedDataEntry, ReferenceDatabase
from rmgpy.species import Species
from rmgpy.statmech import Conformer
from rmgpy.thermo import ThermoData
from rmgpy.thermo.nasa import NASA
from rmgpy.thermo.wilhoit import Wilhoit


################################################################################
FILE_DIR = os.path.dirname(os.path.abspath(__file__))


class TestReferenceSpecies(unittest.TestCase):
    """
    A class for testing that the ReferenceSpecies class functions properly
    """

    def setUp(self):
        self.methane = Species(SMILES='C')
        self.ethane = Species(SMILES='CC')
        self.propane = Species(SMILES='CCC')

        self.thermo_data = ThermoData(H298=(100.0, 'kJ/mol'), S298=(100.0, 'J/(mol*K)'))

    def test_instantiate_reference_species(self):
        """
        Test that a ReferenceSpecies object can be instantiated with the minimal acceptable input, and throws an error
        if the minimal acceptable input is not given.
        """
        ref_spcs = ReferenceSpecies(species=self.ethane)
        self.assertEqual(ref_spcs.smiles, 'CC')
        self.assertEqual(ref_spcs.inchi_key, self.ethane.molecule[0].toInChIKey())

        ref_from_smiles = ReferenceSpecies(smiles='CCC')
        self.assertEqual(ref_from_smiles.smiles, 'CCC')
        self.assertEqual(ref_from_smiles.inchi_key, self.propane.molecule[0].toInChIKey())

        ref_from_inchi = ReferenceSpecies(inchi='InChI=1S/C2H6/c1-2/h1-2H3')
        self.assertEqual(ref_from_inchi.smiles, 'CC')
        self.assertEqual(ref_from_inchi.inchi_key, self.ethane.molecule[0].toInChIKey())

        ref_from_adj = ReferenceSpecies(adjacency_list='1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}\n2 H u0 p0 c0 {1,S}\n'
                                                       '3 H u0 p0 c0 {1,S}\n4 H u0 p0 c0 {1,S}\n5 H u0 p0 c0 {1,S}\n')
        self.assertEqual(ref_from_adj.smiles, 'C')
        self.assertEqual(ref_from_adj.inchi_key, self.methane.molecule[0].toInChIKey())

        with self.assertRaises(ValueError):
            _ = ReferenceSpecies()

    def load_ref_from_yaml(self):
        """
        Test that the example ReferenceSpecies YAML file can be loaded
        """
        ref_spcs = ReferenceSpecies.__new__(ReferenceSpecies)
        ref_spcs.load_yaml(os.path.join(FILE_DIR, 'data', 'species', 'reference_species_example.yml'))

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
        ref_spcs.save_yaml(path=os.path.join(FILE_DIR, 'data'))

        loaded_ref = ReferenceSpecies.__new__(ReferenceSpecies)
        load_path = os.path.join(FILE_DIR, 'data', 'species', '{0}.yml'.format(label))
        loaded_ref.load_yaml(path=load_path)

        self.assertEqual(loaded_ref.smiles, 'CC')

        # Finally, delete this newly created file
        os.remove(load_path)

    def test_reference_data_entry(self):
        """
        Test that the ReferenceDataEntry class functions properly and enforces the standard for storing data
        """
        data_entry = ReferenceDataEntry(self.thermo_data)
        self.assertIsInstance(data_entry.thermo_data, ThermoData)
        self.assertEqual(data_entry.thermo_data.H298.value_si, 100000.0)

        with self.assertRaises(ValueError):
            _ = ReferenceDataEntry({'H298': (100.0, 'kJ/mol')})

    def test_calculated_data_entry(self):
        """
        Test that the CalculatedDataEntry class functions properly and enforces the standard for storing data
        """
        data_entry = CalculatedDataEntry(Conformer(), NASA(), self.thermo_data)
        self.assertEqual(data_entry.thermo_data.H298.value_si, 100000.0)

        data_entry = CalculatedDataEntry(Conformer(), Wilhoit(), self.thermo_data)
        self.assertEqual(data_entry.thermo_data.H298.value_si, 100000.0)

        with self.assertRaises(ValueError):
            _ = CalculatedDataEntry({'xyz': '0 0 0'}, NASA(), self.thermo_data)

        with self.assertRaises(ValueError):
            _ = CalculatedDataEntry(Conformer(), {'coeffs': []}, self.thermo_data)

        with self.assertRaises(ValueError):
            _ = CalculatedDataEntry(Conformer(), NASA(), {'H298': (100.0, 'kJ/mol')})


class TestReferenceDatabase(unittest.TestCase):
    """
    Test that the ReferenceDatabase class functions properly
    """

    def test_load_main_reference_set(self):
        """
        Test that the main reference set can be loaded properly
        """
        database = ReferenceDatabase()
        database.load()
        self.assertIn('main', database.reference_sets)
        self.assertIsInstance(database.reference_sets['main'][0], ReferenceSpecies)

        # Also test that calling load again appends a new set in the database
        data_dir = os.path.join(FILE_DIR, 'data')
        testing_dir = os.path.join(data_dir, 'testing_set')
        example_ref_file = os.path.join(data_dir, 'species', 'reference_species_example.yml')
        spcs_dir = os.path.join(testing_dir, '0')
        spcs_file = os.path.join(spcs_dir, '0.yml')
        if os.path.exists(testing_dir):  # Delete the testing directory if it existed previously
            shutil.rmtree(testing_dir)
        os.mkdir(testing_dir)
        os.mkdir(spcs_dir)
        shutil.copyfile(example_ref_file, spcs_file)
        database.load(paths=[testing_dir])
        self.assertIn('main', database.reference_sets)
        self.assertIn('testing_set', database.reference_sets)

        # Finally, remove the testing directory
        shutil.rmtree(testing_dir)


if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
