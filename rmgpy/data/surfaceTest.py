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
from unittest import TestCase, TestLoader, TextTestRunner

from rmgpy import settings
from rmgpy.data.surface import MetalDatabase
from rmgpy.data.base import Entry
from rmgpy.exceptions import DatabaseError
from rmgpy.quantity import Energy, SurfaceConcentration

###################################################

class TestMetalDatabase(TestCase):

    def setUp(self):
        self.database = MetalDatabase()
        self.database.load(os.path.join(settings['database.directory'], 'surface'))

    def tearDown(self):
        """
        Reset the database & parameters
        """
        import rmgpy.data.rmg
        rmgpy.data.rmg.database = None

    def test_load_metal_library(self):
        """Test we can obtain metal parameters from a library"""

        test_entry = Entry(
            index=1,
            label="Pt111",
            binding_energies={
                'H': Energy(-2.75367887E+00, 'eV/molecule'),
                'C': Energy(-7.02515507E+00, 'eV/molecule'),
                'N': Energy(-4.63224568E+00, 'eV/molecule'),
                'O': Energy(-3.81153179E+00, 'eV/molecule'),
            },
            surface_site_density=SurfaceConcentration(2.483E-09, 'mol/cm^2'),
            facet="111",
            metal="Pt",
            short_desc=u"fcc",
            long_desc=u"""
        Calculated by Katrin Blondal and Bjarne Kreitz at Brown University
            """,
        )

        self.assertEqual(repr(self.database.get_binding_energies(test_entry.label)), repr(test_entry.binding_energies))
        self.assertEqual(repr(self.database.get_surface_site_density(test_entry.label)), repr(test_entry.surface_site_density))

    def test_write_entry_to_database(self):
        """Test we can write an entry to the database"""

        test_entry = Entry(
            index=100,
            label="Me111",
            binding_energies={
                'H': Energy(0., 'eV/molecule'),
                'C': Energy(0., 'eV/molecule'),
                'N': Energy(0., 'eV/molecule'),
                'O': Energy(0., 'eV/molecule'),
            },
            surface_site_density=SurfaceConcentration(0., 'mol/cm^2'),
            facet="111",
            metal="Me",
            short_desc=u"fcc",
            long_desc=u"""
        Test
            """,
        )

        MetalLib = self.database.libraries['surface']
        self.database.add_entry(test_entry)

        # test to see if the entry was added
        self.assertEqual(repr(self.database.get_binding_energies(test_entry.label)), repr(test_entry.binding_energies))
        self.assertEqual(repr(self.database.get_surface_site_density(test_entry.label)), repr(test_entry.surface_site_density))

        # write the new entry
        self.database.save(os.path.join(settings['database.directory'], 'surface'))
        # MetalLib.save_entry(os.path.join(settings['database.directory'], 'surface/libraries/metal.py'), test_entry)

        # test to see if entry was written
        with open(os.path.join(settings['database.directory'], 'surface/libraries/metal.py'), "r") as f:
            if "Me111" in f.read():
                self.database.remove_entry(test_entry)
                self.database.save(os.path.join(settings['database.directory'], 'surface'))
            else:
                raise DatabaseError("Unable to write entry to database.")


    def test_load_from_label(self):
        """Test we can obtain metal parameters from a string"""

        test_pt111 = "Pt111"
        self.assertIsNotNone(self.database.get_binding_energies(test_pt111))

        test_notexistent = "Pt000"
        with self.assertRaises(DatabaseError):
            self.database.get_binding_energies(test_notexistent)

    def test_load_all_entries_on_one_metal(self):
        """Test we can load all entries from the database on one metal"""

        self.assertGreaterEqual(len(self.database.get_all_entries_on_metal("Pt")), 2)
        self.assertGreaterEqual(len(self.database.get_all_entries_on_metal("Ni")), 2)
        self.assertGreaterEqual(len(self.database.get_all_entries_on_metal("Co")), 2)

#####################################################


if __name__ == '__main__':
    suite = TestLoader().loadTestsFromTestCase(TestMetalDatabase)
    TextTestRunner(verbosity=2).run(suite)
