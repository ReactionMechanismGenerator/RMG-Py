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
import unittest

import numpy as np

import rmgpy
from rmgpy.data.rmg import RMGDatabase
from rmgpy.tools.uncertainty import Uncertainty


class TestUncertainty(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """This method is run once before all tests in this class."""
        test_dir = rmgpy.settings['test_data.directory']

        data_dir = os.path.join(test_dir, 'testing_database')
        chem_dir = os.path.join(test_dir, 'parsing_data')
        chemkin_file = os.path.join(chem_dir, 'chem_annotated.inp')
        spc_dict = os.path.join(chem_dir, 'species_dictionary.txt')

        cls.uncertainty = Uncertainty(output_directory='chemDir')
        cls.uncertainty.load_model(chemkin_file, spc_dict)

        # load database properly
        cls.uncertainty.database = RMGDatabase()
        cls.uncertainty.database.load(
            data_dir,
            kinetics_families=['1,2_shiftC', '6_membered_central_C-C_shift', 'Disproportionation', 'H_Abstraction',
                              'Intra_ene_reaction', 'intra_H_migration', 'Intra_R_Add_Exo_scission',
                              'intra_substitutionS_isomerization', 'R_Addition_MultipleBond', 'R_Recombination'],
            kinetics_depositories=['training'],
            thermo_libraries=['primaryThermoLibrary'],
            reaction_libraries=['GRI-Mech3.0'],
        )

        # Prepare the database by loading training reactions and averaging the rate rules verbosely
        for family in cls.uncertainty.database.kinetics.families.values():
            if not family.auto_generated:
                family.add_rules_from_training(thermo_database=cls.uncertainty.database.thermo)
                family.fill_rules_by_averaging_up(verbose=True)

    @classmethod
    def tearDownClass(cls):
        """This method is run once after all tests in this class."""
        # Reset module level database
        import rmgpy.data.rmg
        rmgpy.data.rmg.database = None

    def test_uncertainty_assignment(self):
        """
        Test that the thermo and kinetic parameter uncertainties can be properly assigned.
        """
        # Step 1: parse comments for sources
        self.uncertainty.extract_sources_from_model()
        self.assertEqual(len(self.uncertainty.species_sources_dict), len(self.uncertainty.species_list))
        self.assertEqual(len(self.uncertainty.reaction_sources_dict), len(self.uncertainty.reaction_list))

        # Step 2: compile sources to obtain overall list
        self.uncertainty.compile_all_sources()

        # Check thermo sources
        grp_expected = {
            'O2s-CsH', 'O2s-OsH', 'O2d-Cd', 'Cds-OdHH', 'Cds-(Cdd-O2d)HH', 'Ct-CtH', 'Cds-CdsHH', 'Cs-OsHHH',
            'Cs-CsHHH', 'Cdd-CdsOd'
        }
        rad_expected = {'Acetyl', 'HOOJ', 'Cds_P', 'CCJ', 'CsJOH', 'CJ3'}
        other_expected = {'ketene', 'R'}
        self.assertEqual(set(self.uncertainty.all_thermo_sources), {'GAV', 'Library', 'QM'})
        self.assertEqual(set(self.uncertainty.all_thermo_sources['GAV']), {'group', 'radical', 'other'})
        grp = set([e.label for e in self.uncertainty.all_thermo_sources['GAV']['group']])
        rad = set([e.label for e in self.uncertainty.all_thermo_sources['GAV']['radical']])
        other = set([e.label for e in self.uncertainty.all_thermo_sources['GAV']['other']])
        self.assertEqual(grp, grp_expected)
        self.assertEqual(rad, rad_expected)
        self.assertEqual(other, other_expected)
        self.assertEqual(sorted(self.uncertainty.all_thermo_sources['Library']), [0, 1, 5, 13, 14])
        self.assertFalse(self.uncertainty.all_thermo_sources['QM'])

        # Check kinetics sources
        rr_expected = {
            'O_rad/NonDeO;O_Csrad', 'Ct_rad/Ct;O_Csrad', 'O_atom_triplet;O_Csrad', 'C_rad/H2/Cd;O_Csrad',
            'CH2_triplet;O_Csrad', 'H_rad;O_Csrad', 'O_rad/NonDeC;O_Csrad', 'C_rad/Cs3;O_Csrad', 'Cd_pri_rad;O_Csrad',
            'O2b;O_Csrad', 'O_pri_rad;Cmethyl_Csrad', 'C_rad/H/NonDeC;O_Csrad', 'O_pri_rad;O_Csrad',
            'C_methyl;O_Csrad', 'C_rad/H2/Cs;O_Csrad', 'C_rad/H2/O;O_Csrad', 'CO_pri_rad;O_Csrad'
        }
        self.assertEqual(set(self.uncertainty.all_kinetic_sources), {'Rate Rules', 'Training', 'Library', 'PDep'})
        self.assertEqual(list(self.uncertainty.all_kinetic_sources['Rate Rules'].keys()), ['Disproportionation'])
        rr = set([e.label for e in self.uncertainty.all_kinetic_sources['Rate Rules']['Disproportionation']])
        self.assertEqual(rr, rr_expected)
        self.assertEqual(list(self.uncertainty.all_kinetic_sources['Training'].keys()), ['Disproportionation'])
        self.assertEqual(self.uncertainty.all_kinetic_sources['Library'], [0])
        self.assertEqual(self.uncertainty.all_kinetic_sources['PDep'], [4])

        # Step 3: assign and propagate uncertainties
        self.uncertainty.assign_parameter_uncertainties()

        thermo_unc = self.uncertainty.thermo_input_uncertainties
        kinetic_unc = self.uncertainty.kinetic_input_uncertainties

        np.testing.assert_allclose(thermo_unc, [1.5, 1.5, 2.0, 1.9, 3.1, 1.5, 1.9, 2.0, 2.0, 1.9, 2.2, 1.9, 2.0, 1.5],
                                   rtol=1e-4)
        np.testing.assert_allclose(kinetic_unc, [0.5, 1.5, 5.806571, 0.5, 2.0], rtol=1e-4)
