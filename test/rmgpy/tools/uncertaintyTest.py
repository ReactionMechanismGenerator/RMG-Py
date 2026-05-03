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

import os


import numpy as np

import rmgpy
from rmgpy.data.rmg import RMGDatabase
from rmgpy.tools.uncertainty import Uncertainty


class TestUncertainty:
    @classmethod
    def setup_class(cls):
        """This method is run once before all tests in this class."""
        test_dir = rmgpy.settings["test_data.directory"]

        data_dir = os.path.join(test_dir, "testing_database")
        chem_dir = os.path.join(test_dir, "parsing_data")
        chemkin_file = os.path.join(chem_dir, "chem_annotated.inp")
        spc_dict = os.path.join(chem_dir, "species_dictionary.txt")

        cls.uncertainty = Uncertainty(output_directory="chemDir")
        cls.uncertainty.load_model(chemkin_file, spc_dict)

        # load database properly
        cls.uncertainty.database = RMGDatabase()
        cls.uncertainty.database.load(
            data_dir,
            kinetics_families=[
                "1,2_shiftC",
                "6_membered_central_C-C_shift",
                "Disproportionation",
                "H_Abstraction",
                "Intra_ene_reaction",
                "intra_H_migration",
                "Intra_R_Add_Exo_scission",
                "intra_substitutionS_isomerization",
                "R_Addition_MultipleBond",
                "R_Recombination",
            ],
            kinetics_depositories=["training"],
            thermo_libraries=["primaryThermoLibrary"],
            reaction_libraries=["GRI-Mech3.0"],
        )

        # Prepare the database by loading training reactions and averaging the rate rules verbosely
        for family in cls.uncertainty.database.kinetics.families.values():
            if not family.auto_generated:
                family.add_rules_from_training(thermo_database=cls.uncertainty.database.thermo)
                family.fill_rules_by_averaging_up(verbose=True)

    @classmethod
    def teardown_class(cls):
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
        assert len(self.uncertainty.species_sources_dict) == len(self.uncertainty.species_list)
        assert len(self.uncertainty.reaction_sources_dict) == len(self.uncertainty.reaction_list)

        # Step 2: compile sources to obtain overall list
        self.uncertainty.compile_all_sources()

        # Check thermo sources
        grp_expected = {
            "O2s-CsH",
            "O2s-OsH",
            "O2d-Cd",
            "Cds-OdHH",
            "Cds-(Cdd-O2d)HH",
            "Ct-CtH",
            "Cds-CdsHH",
            "Cs-OsHHH",
            "Cs-CsHHH",
            "Cdd-CdsOd",
            "Cds-CdsCsH",
            "Cs-(Cds-Cds)HHH",
            "Cs-CsCsHH",
        }
        rad_expected = {
            "Acetyl",
            "HOOJ",
            "Cds_P",
            "CCJ",
            "CsJOH",
            "CJ3",
            "Allyl_P",
            "CCJC",
            "CH3",
            "RCCJ",
        }
        other_expected = {"ketene", "R"}
        assert set(self.uncertainty.all_thermo_sources) == {"GAV", "Library", "QM", "ADS", "Surface_Library"}
        assert set(self.uncertainty.all_thermo_sources["GAV"]) == {"group", "radical", "other"}
        grp = set([e.label for e in self.uncertainty.all_thermo_sources["GAV"]["group"]])
        rad = set([e.label for e in self.uncertainty.all_thermo_sources["GAV"]["radical"]])
        other = set([e.label for e in self.uncertainty.all_thermo_sources["GAV"]["other"]])
        assert grp == grp_expected
        assert rad == rad_expected
        assert other == other_expected
        assert sorted(self.uncertainty.all_thermo_sources["Library"]) == [0, 1, 5, 13, 16]
        assert not self.uncertainty.all_thermo_sources["QM"]

        # Check kinetics sources
        Disproportionation_rr_expected = {
            'Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_N-1R!H-inRing_N-Sp-6C-4C',
            'Root_Ext-2R!H-R_2R!H->C_4R->C',
        }
        H_Abstraction_rr_expected = {
            'C/H3/Cs;C_methyl',
            'C/H3/Cs\\H2\\Cs|O;Cd_Cd\\H2_rad/Cs',
            'C/H3/Cs\\H2\\O;C_methyl',
            'C/H3/Cs\\H3;C_rad/H2/Cs\\H3',
            'C/H3/Cs\\H3;C_rad/H2/Cs\\H\\Cs\\Cs|O',
            'C/H3/Cs\\H3;Cd_Cd\\H2_pri_rad',
        }
        assert set(self.uncertainty.all_kinetic_sources) == {"Rate Rules", "Training", "Library", "PDep", "Surface_Library"}
        assert set(self.uncertainty.all_kinetic_sources["Rate Rules"].keys()) == {"Disproportionation", "H_Abstraction"}
        rr = set([e.label for e in self.uncertainty.all_kinetic_sources["Rate Rules"]["Disproportionation"]])
        assert rr == Disproportionation_rr_expected
        rr = set([e.label for e in self.uncertainty.all_kinetic_sources["Rate Rules"]["H_Abstraction"]])
        assert rr == H_Abstraction_rr_expected
        assert set(self.uncertainty.all_kinetic_sources["Training"].keys()) == {"Disproportionation", "H_Abstraction"}
        assert self.uncertainty.all_kinetic_sources["Library"] == [0]
        assert self.uncertainty.all_kinetic_sources["PDep"] == [6]

        # Step 3: assign and propagate uncertainties
        self.uncertainty.assign_parameter_uncertainties()

        thermo_unc = self.uncertainty.thermo_input_uncertainties
        kinetic_unc = self.uncertainty.kinetic_input_uncertainties

        np.testing.assert_allclose(
            thermo_unc,
            [1.5, 1.5, 2.61966, 2.51994, 2.23886, 1.5, 2.30761, 2.41611, 2.61966, 2.51994, 2.61966, 2.51994, 2.61966, 1.5, 2.23886, 2.30761, 1.5, 2.61966, 2.07366, 2.19376, 2.19376, 2.30761, 1.94616, 2.07366, 2.07366],
            rtol=1e-4,
        )
        np.testing.assert_allclose(
            kinetic_unc,
            [0.5, 1.118, 1.9783, 1.9783, 1.5363, 0.5, 2.0, 5.9369, 5.9369, 0.5],
            rtol=1e-4
        )

    def test_specific_species_uncertainties(self):
        """
        Test uncertainties for a few specific examples
        """

        expected_results = {  # order is (total_uncertainty, [group_names], [group_counts])
            'CCCC': (2.5199409675625337, ['Cs-CsCsHH', 'Cs-CsHHH'], [2, 2]),
            'CCCCCCCCCC': (6.091048438487417, ['Cs-CsCsHH', 'Cs-CsHHH'], [8, 2]),
            'CC(OO)CC': (2.5199409675625337, ['O2s-OsCs', 'O2s-OsH', 'Cs-CsCsOsH', 'Cs-CsCsHH', 'Cs-CsHHH'], [1, 1, 1, 1, 2]),
            'C=NCC': (2.07365649035707, ['N3d-CdCs', 'Cs-(N3dCd)CsHH', 'Cs-CsHHH', 'Cd-N3dHH'], [1, 1, 1, 1]),
            'C=C': (2.07365649035707, ['Cds-CdsHH'], [2]),
            'C*': (7.271261019245562, ['CH3'], [1]),  # Gas library + radical + adsorption correction
            'O=[CH]*': (7.150786643440007, ['Cds-OdHH', 'HCdsJO'], [1, 1]),  # GAV + radical + adsorption correction
        }

        uncertainty = rmgpy.tools.uncertainty.Uncertainty()
        uncertainty.database = self.uncertainty.database  # use the same database as the main test
        new_species_list = [rmgpy.species.Species(smiles=spc) for spc in expected_results.keys()]
        for spc in new_species_list:
            spc.thermo = self.uncertainty.database.thermo.get_thermo_data(spc)
            if not isinstance(spc.thermo, rmgpy.thermo.NASA):
                spc.thermo = spc.thermo.to_nasa(Tmin=298, Tmax=3000, Tint=1000)
        uncertainty.species_list = new_species_list
        uncertainty.reaction_list = []  # no need to test kinetics here

        uncertainty.extract_sources_from_model()  # this will populate the sources dict with the new species
        uncertainty.assign_parameter_uncertainties()  # this will assign the uncertainties based on the sources and the database

        for i, sp in enumerate(uncertainty.species_list):
            source = uncertainty.species_sources_dict[sp]

            actual_group_data = []
            for group_type, group_entries in source['GAV'].items():
                for group, count in group_entries:
                    actual_group_data.append((group.label, count))
            expected_result = expected_results[sp.smiles]
            expected_group_data = sorted(zip(expected_result[1], expected_result[2]))
            actual_group_data = sorted(actual_group_data)
            assert np.isclose(uncertainty.thermo_input_uncertainties[i], expected_result[0], rtol=1e-4)
            assert actual_group_data == expected_group_data
