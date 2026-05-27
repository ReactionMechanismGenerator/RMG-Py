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
import copy

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

        cls.uncertainty = Uncertainty(output_directory=os.path.abspath(os.path.join(os.path.dirname(__file__), "chemDir")))
        cls.uncertainty.load_model(chemkin_file, spc_dict)
        for i in range(len(cls.uncertainty.species_list)):
            cls.uncertainty.species_list[i].index = i  # local analysis depends on species being indexed

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

        expected_uncorrelated_thermo_uncertainties = np.array([1.5, 1.5, 2.61966, 2.51994, 2.23886, 1.5, 2.30761, 2.41611, 2.61966, 2.51994, 2.61966, 2.51994, 2.61966, 1.5, 2.23886, 2.30761, 1.5, 2.61966, 2.07366, 2.19376, 2.19376, 2.30761, 1.94616, 2.07366, 2.07366])
        expected_uncorrelated_kinetic_uncertainties = np.array([0.5, 1.118, 1.9783, 1.9783, 1.5363, 0.5, 2.0, 1.5363, 1.5363, 0.5])
        np.testing.assert_allclose(thermo_unc, expected_uncorrelated_thermo_uncertainties, rtol=1e-4)
        np.testing.assert_allclose(kinetic_unc, expected_uncorrelated_kinetic_uncertainties, rtol=1e-4)

        # ---------------------------- Now repeat for assign_intermediate_uncertainties -----------------------------
        # uncorrelated
        self.uncertainty.assign_intermediate_uncertainties(correlated=False)
        intermediate_thermo_unc = self.uncertainty.thermo_intermediate_uncertainties
        intermediate_kinetic_unc = self.uncertainty.kinetic_intermediate_uncertainties
        np.testing.assert_allclose(intermediate_thermo_unc, expected_uncorrelated_thermo_uncertainties, rtol=1e-4)
        np.testing.assert_allclose(intermediate_kinetic_unc, expected_uncorrelated_kinetic_uncertainties, rtol=1e-4)

        # correlated
        self.uncertainty.assign_intermediate_uncertainties(correlated=True)
        
        # do a spot check on some of the intermediates (dG/dq) these are derivatives, not uncertainties
        # Thermo library example
        assert self.uncertainty.thermo_intermediate_uncertainties[0].keys() == {'Library O(0)'}
        assert self.uncertainty.thermo_intermediate_uncertainties[0]['Library O(0)'] == 1

        # Thermo GAV example
        assert tuple(sorted(self.uncertainty.thermo_intermediate_uncertainties[2].keys())) == ('Estimation HO2(2)', 'Group(group) O2s-OsH', 'Group(other) R', 'Group(radical) HOOJ')
        assert self.uncertainty.thermo_intermediate_uncertainties[2]['Estimation HO2(2)'] == 1
        assert self.uncertainty.thermo_intermediate_uncertainties[2]['Group(group) O2s-OsH'] == 2
        assert self.uncertainty.thermo_intermediate_uncertainties[2]['Group(other) R'] == 2
        assert self.uncertainty.thermo_intermediate_uncertainties[2]['Group(radical) HOOJ'] == 1

        # Thermo library + GAV
        assert tuple(sorted(self.uncertainty.thermo_intermediate_uncertainties[14].keys())) == ('Estimation CH3(14)', 'Group(radical) CH3', 'Library CH4(16)')
        assert self.uncertainty.thermo_intermediate_uncertainties[14]['Estimation CH3(14)'] == 1
        assert self.uncertainty.thermo_intermediate_uncertainties[14]['Group(radical) CH3'] == 1
        assert self.uncertainty.thermo_intermediate_uncertainties[14]['Library CH4(16)'] == 1

        # Kinetics library
        assert self.uncertainty.kinetic_intermediate_uncertainties[0].keys() == {'Library O(0)+H2O2(3)<=>OH(1)+HO2(2)'}

        # Rate rule (exact)
        assert tuple(sorted(self.uncertainty.kinetic_intermediate_uncertainties[1].keys())) == ('Estimation Family CH3(14)+PC3H7(15)<=>CH4(16)+CH2CH2CH2(17)', 'Rate Rule H_Abstraction C/H3/Cs;C_methyl')
        assert self.uncertainty.kinetic_intermediate_uncertainties[1]['Estimation Family CH3(14)+PC3H7(15)<=>CH4(16)+CH2CH2CH2(17)'] == 1
        assert self.uncertainty.kinetic_intermediate_uncertainties[1]['Rate Rule H_Abstraction C/H3/Cs;C_methyl'] == 1
        
        # Rate rule (non-exact, multiple rule weights)
        assert tuple(sorted(self.uncertainty.kinetic_intermediate_uncertainties[3].keys())) == (
            'Estimation Family C2H3(20)+C3H8(19)<=>C2H4(11)+PC3H7(15)',
            'Estimation Nonexact C2H3(20)+C3H8(19)<=>C2H4(11)+PC3H7(15)',
            'Rate Rule H_Abstraction C/H3/Cs\\H2\\Cs|O;Cd_Cd\\H2_rad/Cs',
            'Rate Rule H_Abstraction C/H3/Cs\\H3;Cd_Cd\\H2_pri_rad',
        )
        assert self.uncertainty.kinetic_intermediate_uncertainties[3]['Estimation Family C2H3(20)+C3H8(19)<=>C2H4(11)+PC3H7(15)'] == 1
        assert np.isclose(self.uncertainty.kinetic_intermediate_uncertainties[3]['Estimation Nonexact C2H3(20)+C3H8(19)<=>C2H4(11)+PC3H7(15)'], 0.4771212547, rtol=1e-4)
        assert self.uncertainty.kinetic_intermediate_uncertainties[3]['Rate Rule H_Abstraction C/H3/Cs\\H2\\Cs|O;Cd_Cd\\H2_rad/Cs'] == 0.5
        assert self.uncertainty.kinetic_intermediate_uncertainties[3]['Rate Rule H_Abstraction C/H3/Cs\\H3;Cd_Cd\\H2_pri_rad'] == 0.5

        # Training reaction
        assert self.uncertainty.kinetic_intermediate_uncertainties[5].keys() == {'Training H_Abstraction CH3(14)+C2H6(18)<=>CH4(16)+C2H5(12)'}
        assert self.uncertainty.kinetic_intermediate_uncertainties[5]['Training H_Abstraction CH3(14)+C2H6(18)<=>CH4(16)+C2H5(12)'] == 1

        # PDEP
        assert self.uncertainty.kinetic_intermediate_uncertainties[6].keys() == {'PDep HCCO(10)(+M)<=>O(0)+C2H(8)(+M)'}
        assert self.uncertainty.kinetic_intermediate_uncertainties[6]['PDep HCCO(10)(+M)<=>O(0)+C2H(8)(+M)'] == 1

        # correlated uncertainties should match uncorrelated, so check diagonal of covariance matrix
        thermo_covariance = np.sqrt(self.uncertainty.get_thermo_covariance_matrix().diagonal())
        kinetic_covariance = np.sqrt(self.uncertainty.get_kinetic_covariance_matrix().diagonal())
        assert np.isclose(thermo_covariance, expected_uncorrelated_thermo_uncertainties, rtol=1e-4).all()
        assert np.isclose(kinetic_covariance, expected_uncorrelated_kinetic_uncertainties, rtol=1e-4).all()

    def test_source_correlations(self):
        # Check some examples of different species containing the same sources

        # ------------------------------------------------------------------------------
        # Make sure CH3 (Library + Radical) has a library index/value in common with CH4
        i_CH4 = rmgpy.tools.uncertainty.get_i_thing(rmgpy.species.Species(smiles='C'), self.uncertainty.species_list)
        assert i_CH4 >= 0

        i_CH3 = rmgpy.tools.uncertainty.get_i_thing(rmgpy.species.Species(smiles='[CH3]'), self.uncertainty.species_list)
        assert i_CH3 >= 0

        self.uncertainty.extract_sources_from_model()
        self.uncertainty.assign_parameter_uncertainties(correlated=True)
        
        src1 = self.uncertainty.species_sources_dict[self.uncertainty.species_list[i_CH4]]  # CH4
        src2 = self.uncertainty.species_sources_dict[self.uncertainty.species_list[i_CH3]]  # CH3

        assert 'Library' in src1
        assert 'Library' in src2
        assert 'GAV' in src2
        assert src1['Library'] == src2['Library']  # make sure they refer to the same library source

        # -----------------------------------------------------------------------------
        # Make sure CH3X (Library + GAV + Adsorption Correction) has a library index/value in common with CH4 (Library)
        i_CH3X = rmgpy.tools.uncertainty.get_i_thing(rmgpy.species.Species(smiles='C*'), self.uncertainty.species_list)
        assert i_CH3X == -1
        # This is not in the model, so add it to the species list
        CH3X = rmgpy.species.Species(smiles='C*')
        CH3X.thermo = self.uncertainty.database.thermo.get_thermo_data(CH3X)
        self.uncertainty.species_list.append(CH3X)
        i_CH3X = rmgpy.tools.uncertainty.get_i_thing(CH3X, self.uncertainty.species_list)
        assert i_CH3X >= 0

        self.uncertainty.extract_sources_from_model()
        self.uncertainty.assign_parameter_uncertainties(correlated=True)

        src1 = self.uncertainty.species_sources_dict[self.uncertainty.species_list[i_CH4]]  # CH4
        src2 = self.uncertainty.species_sources_dict[self.uncertainty.species_list[i_CH3X]]  # CH3X

        assert 'Library' in src1
        assert 'Library' in src2
        assert 'ADS' in src2
        assert 'GAV' in src2
        assert src1['Library'] == src2['Library']  # make sure they refer to the same library source
        self.uncertainty.species_list.pop()  # remove the extra species so it doesn't affect other tests

    def test_local_analysis(self):
        """
        Test to run uncorrelated and then correlated local_analysis and make sure the results are expected
        """
        # variances are listed in decreasing order
        # names are listed in order of decreasing variance contribution
        expected_uncorrelated_total_variance = 1.8329056941266446
        expected_uncorrelated_thermo_variances = np.array([0.17092419, 0.09781627, 0.06186124, 0.04856985, 0.00391013, 0.00306632, 0.00041446, 9.953e-05])
        expected_uncorrelated_kinetics_variances = np.array([1.1311145, 0.15888459, 0.085189, 0.02022449, 0.01687337, 0.01605351, 0.01588366, 0.0010829, 0.00080811, 0.00012957])
        expected_correlated_total_variance = 1.7732795017083922
        expected_correlated_thermo_variances = np.array([0.09145902, 0.07672388, 0.04856985, 0.04573167, 0.03236887, 0.01747643, 0.01098087, 0.00143231, 0.0013764, 0.00031352, 0.00028968, 0.00014685, 0.0001338, 3.263e-05])
        expected_correlated_kinetics_variances = np.array([0.53202843, 0.47926886, 0.11981721, 0.11321232, 0.06070094, 0.04059757, 0.02176716, 0.01687337, 0.0161796, 0.01605351, 0.007471, 0.00673013, 0.0040449, 0.00253735, 0.00253735, 0.00168253, 0.00136045, 0.00136045, 0.00080811, 0.00050935, 0.00045884, 0.00012957, 0.00011471])
        expected_uncorrelated_thermo_labels = [
            'dln[C2H6(18)]/dG[CH(4)]',
            'dln[C2H6(18)]/dG[C2H3(20)]',
            'dln[C2H6(18)]/dG[C2H6(18)]',
            'dln[C2H6(18)]/dG[CH2(5)]',
            'dln[C2H6(18)]/dG[CH4(16)]',
            'dln[C2H6(18)]/dG[CH3(14)]',
            'dln[C2H6(18)]/dG[C2H4(11)]',
            'dln[C2H6(18)]/dG[C2H5(12)]',
        ]
        expected_uncorrelated_kinetics_labels = [
            'k8: C2H5(12)+CH3CHCH3(21)<=>C2H6(18)+C3H6(22)',
            'k3: C2H6(18)+PC3H7(15)<=>C2H5(12)+C3H8(19)',
            'k4: C2H3(20)+C3H8(19)<=>C2H4(11)+PC3H7(15)',
            'k2: CH3(14)+PC3H7(15)<=>CH4(16)+CH2CH2CH2(17)',
            'k1: O(0)+H2O2(3)<=>OH(1)+HO2(2)',
            'k7: HCCO(10)(+M)<=>O(0)+C2H(8)(+M)',
            'k5: CH3(14)+C3H8(19)<=>CH4(16)+PC3H7(15)',
            'k9: C3H5(24)+CH2CH2CH2(17)<=>C3H5(23)+C3H6(22)',
            'k10: CH3(14)+C2H5(12)<=>CH4(16)+C2H4(11)',
            'k6: CH3(14)+C2H6(18)<=>CH4(16)+C2H5(12)',
        ]
        expected_correlated_thermo_labels = [
            'Library CH4(16)',
            'Estimation CH(4)',
            'Library CH2(5)',
            'Estimation C2H3(20)',
            'Estimation C2H6(18)',
            'Group(radical) CJ3',
            'Group(radical) CCJ',
            'Group(group) Cs-CsHHH',
            'Estimation CH3(14)',
            'Group(radical) CH3',
            'Group(other) R',
            'Estimation C2H4(11)',
            'Group(group) Cds-CdsHH',
            'Estimation C2H5(12)',
        ]
        expected_correlated_kinetics_labels = [
            'Estimation Nonexact C2H5(12)+CH3CHCH3(21)<=>C2H6(18)+C3H6(22)',
            'Estimation Family C2H5(12)+CH3CHCH3(21)<=>C2H6(18)+C3H6(22)',
            'Rate Rule Disproportionation Root_Ext-2R!H-R_2R!H->C_4R->C',
            'Estimation Nonexact C2H6(18)+PC3H7(15)<=>C2H5(12)+C3H8(19)',
            'Estimation Nonexact C2H3(20)+C3H8(19)<=>C2H4(11)+PC3H7(15)',
            'Estimation Family C2H6(18)+PC3H7(15)<=>C2H5(12)+C3H8(19)',
            'Estimation Family C2H3(20)+C3H8(19)<=>C2H4(11)+PC3H7(15)',
            'Library O(0)+H2O2(3)<=>OH(1)+HO2(2)',
            'Estimation Family CH3(14)+PC3H7(15)<=>CH4(16)+CH2CH2CH2(17)',
            'PDep HCCO(10)(+M)<=>O(0)+C2H(8)(+M)',
            'Estimation Nonexact CH3(14)+C3H8(19)<=>CH4(16)+PC3H7(15)',
            'Estimation Family CH3(14)+C3H8(19)<=>CH4(16)+PC3H7(15)',
            'Rate Rule H_Abstraction C/H3/Cs;C_methyl',
            'Rate Rule H_Abstraction C/H3/Cs\\H3;C_rad/H2/Cs\\H\\Cs\\Cs|O',
            'Rate Rule H_Abstraction C/H3/Cs\\H3;C_rad/H2/Cs\\H3',
            'Rate Rule H_Abstraction C/H3/Cs\\H2\\O;C_methyl',
            'Rate Rule H_Abstraction C/H3/Cs\\H3;Cd_Cd\\H2_pri_rad',
            'Rate Rule H_Abstraction C/H3/Cs\\H2\\Cs|O;Cd_Cd\\H2_rad/Cs',
            'Training Disproportionation CH3(14)+C2H5(12)<=>CH4(16)+C2H4(11)',
            'Estimation Nonexact C3H5(24)+CH2CH2CH2(17)<=>C3H5(23)+C3H6(22)',
            'Estimation Family C3H5(24)+CH2CH2CH2(17)<=>C3H5(23)+C3H6(22)',
            'Training H_Abstraction CH3(14)+C2H6(18)<=>CH4(16)+C2H5(12)',
            'Rate Rule Disproportionation Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_N-1R!H-inRing_N-Sp-6C-4C',
        ]

        sensitive_species = [self.uncertainty.species_list[18]]

        # uncorrelated analysis first
        self.uncertainty.assign_parameter_uncertainties()
        output = self.uncertainty.local_analysis(sensitive_species=sensitive_species)
        total_variance, kinetic_uncertainty, thermo_uncertainty = output[sensitive_species[0]]
        assert np.isclose(total_variance, expected_uncorrelated_total_variance)

        # order of kinetic or thermo uncertainty is not guaranteed, this sorts by contribution
        kinetic_variances = [r[2] for r in kinetic_uncertainty]
        kinetics_names = [r[0] for r in kinetic_uncertainty]
        sorted_kinetics_names = [x for _, x in sorted(zip(kinetic_variances, kinetics_names))][::-1]
        sorted_kinetic_variances = sorted(kinetic_variances, reverse=True)
        assert np.isclose(sorted_kinetic_variances, expected_uncorrelated_kinetics_variances).all()
        assert sorted_kinetics_names == expected_uncorrelated_kinetics_labels

        thermo_variances = [s[2] for s in thermo_uncertainty]
        thermo_names = [s[0] for s in thermo_uncertainty]
        sorted_thermo_names = [x for _, x in sorted(zip(thermo_variances, thermo_names))][::-1]
        sorted_thermo_variances = sorted(thermo_variances, reverse=True)
        assert np.isclose(sorted_thermo_variances, expected_uncorrelated_thermo_variances).all()
        assert sorted_thermo_names == expected_uncorrelated_thermo_labels

        # now repeat for correlated analysis
        self.uncertainty.assign_parameter_uncertainties(correlated=True)
        output = self.uncertainty.local_analysis(sensitive_species=sensitive_species, correlated=True)
        total_variance, kinetic_uncertainty, thermo_uncertainty = output[sensitive_species[0]]
        assert np.isclose(total_variance, expected_correlated_total_variance)

        # order of kinetic or thermo uncertainty is not guaranteed, this sorts by contribution
        kinetic_variances = [r[2] for r in kinetic_uncertainty]
        kinetics_names = [r[0] for r in kinetic_uncertainty]
        sorted_kinetic_variances = sorted(kinetic_variances, reverse=True)
        sorted_kinetics_names = [x for _, x in sorted(zip(kinetic_variances, kinetics_names))][::-1]
        assert np.isclose(sorted_kinetic_variances, expected_correlated_kinetics_variances).all()
        assert sorted_kinetics_names == expected_correlated_kinetics_labels

        thermo_variances = [s[2] for s in thermo_uncertainty]
        thermo_names = [s[0] for s in thermo_uncertainty]
        sorted_thermo_variances = sorted(thermo_variances, reverse=True)
        sorted_thermo_names = [x for _, x in sorted(zip(thermo_variances, thermo_names))][::-1]
        assert np.isclose(sorted_thermo_variances, expected_correlated_thermo_variances).all()
        assert sorted_thermo_names == expected_correlated_thermo_labels

        # -------------------- repeat the exact same test for new formulation --------------------------
        # uncorrelated analysis first
        self.uncertainty.assign_intermediate_uncertainties()
        output = self.uncertainty.local_analysis_intermediate(sensitive_species=sensitive_species)
        total_variance, kinetic_uncertainty, thermo_uncertainty = output[sensitive_species[0]]
        assert np.isclose(total_variance, expected_uncorrelated_total_variance)

        # order of kinetic or thermo uncertainty is not guaranteed, this sorts by contribution
        kinetic_variances = [r[2] for r in kinetic_uncertainty]
        kinetics_names = [r[0] for r in kinetic_uncertainty]
        sorted_kinetics_names = [x for _, x in sorted(zip(kinetic_variances, kinetics_names))][::-1]
        sorted_kinetic_variances = sorted(kinetic_variances, reverse=True)
        assert np.isclose(sorted_kinetic_variances, expected_uncorrelated_kinetics_variances).all()
        assert sorted_kinetics_names == expected_uncorrelated_kinetics_labels

        thermo_variances = [s[2] for s in thermo_uncertainty]
        thermo_names = [s[0] for s in thermo_uncertainty]
        sorted_thermo_names = [x for _, x in sorted(zip(thermo_variances, thermo_names))][::-1]
        sorted_thermo_variances = sorted(thermo_variances, reverse=True)
        assert np.isclose(sorted_thermo_variances, expected_uncorrelated_thermo_variances).all()
        assert sorted_thermo_names == expected_uncorrelated_thermo_labels

        # now repeat for correlated analysis
        self.uncertainty.assign_intermediate_uncertainties(correlated=True)
        output = self.uncertainty.local_analysis_intermediate(sensitive_species=sensitive_species, correlated=True)
        total_variance, kinetic_uncertainty, thermo_uncertainty = output[sensitive_species[0]]
        assert np.isclose(total_variance, expected_correlated_total_variance)

        # order of kinetic or thermo uncertainty is not guaranteed, this sorts by contribution
        kinetic_variances = [r[2] for r in kinetic_uncertainty]
        kinetics_names = [r[0] for r in kinetic_uncertainty]
        sorted_kinetic_variances = sorted(kinetic_variances, reverse=True)
        sorted_kinetics_names = [x for _, x in sorted(zip(kinetic_variances, kinetics_names))][::-1]
        assert np.isclose(sorted_kinetic_variances, expected_correlated_kinetics_variances).all()
        assert sorted_kinetics_names == expected_correlated_kinetics_labels

        thermo_variances = [s[2] for s in thermo_uncertainty]
        thermo_names = [s[0] for s in thermo_uncertainty]
        sorted_thermo_variances = sorted(thermo_variances, reverse=True)
        sorted_thermo_names = [x for _, x in sorted(zip(thermo_variances, thermo_names))][::-1]
        assert np.isclose(sorted_thermo_variances, expected_correlated_thermo_variances).all()
        assert sorted_thermo_names == expected_correlated_thermo_labels

    def test_covariance_matrices(self):
        """
        Test that the covariance matrices are being constructed correctly, and that the correlated uncertainties are different from the uncorrelated ones
        """

        # have to add an extra reaction to see any kinetic correlations
        # copy reaction 4 and change the index so it is a new reaction, but with the same source (rate rule) as the original reaction
        extra_reaction = copy.deepcopy(self.uncertainty.reaction_list[4])
        self.uncertainty.reaction_list.append(extra_reaction)
        try:  # this will still error out if there's a problem, but will reset the reaction list so it doesn't affect other tests
            self.uncertainty.extract_sources_from_model()  # this will assign the same source to the new reaction as the original reaction

            self.uncertainty.assign_parameter_uncertainties(correlated=False)
            uncorrelated_thermo_inputs = np.array(self.uncertainty.thermo_input_uncertainties)
            uncorrelated_kinetic_inputs = np.array(self.uncertainty.kinetic_input_uncertainties)

            self.uncertainty.assign_intermediate_uncertainties(correlated=False)
            uncorrelated_thermo_covariance = self.uncertainty.get_thermo_covariance_matrix()
            uncorrelated_kinetic_covariance = self.uncertainty.get_kinetic_covariance_matrix()
            
            self.uncertainty.assign_intermediate_uncertainties(correlated=True)
            correlated_thermo_covariance = self.uncertainty.get_thermo_covariance_matrix()
            correlated_kinetic_covariance = self.uncertainty.get_kinetic_covariance_matrix()
            Sigma_ww_thermo = self.uncertainty._get_intermediate_thermo_covariance_matrix()
            Sigma_ww_kinetics = self.uncertainty._get_intermediate_kinetics_covariance_matrix()
        finally:
            self.uncertainty.reaction_list.pop()  # remove the extra reaction so it doesn't affect other tests
        
        # check that the diagonal elements of the correlated and uncorrelated covariance matrices are the same and equal to the squares of the input uncertainties
        np.testing.assert_allclose(np.diag(uncorrelated_thermo_covariance), np.float_power(uncorrelated_thermo_inputs, 2.0), rtol=1e-4)
        np.testing.assert_allclose(np.diag(correlated_thermo_covariance), np.float_power(uncorrelated_thermo_inputs, 2.0), rtol=1e-4)
        np.testing.assert_allclose(np.diag(uncorrelated_kinetic_covariance), np.float_power(uncorrelated_kinetic_inputs, 2.0), rtol=1e-4)
        np.testing.assert_allclose(np.diag(correlated_kinetic_covariance), np.float_power(uncorrelated_kinetic_inputs, 2.0), rtol=1e-4)

        # check that the off-diagonal elements of the uncorrelated covariance matrix are zero
        off_diagonal_kinetic_uncorrelated = uncorrelated_kinetic_covariance - np.diag(np.diag(uncorrelated_kinetic_covariance))
        assert np.allclose(off_diagonal_kinetic_uncorrelated, 0, atol=1e-8)
        off_diagonal_thermo_uncorrelated = uncorrelated_thermo_covariance - np.diag(np.diag(uncorrelated_thermo_covariance))
        assert np.allclose(off_diagonal_thermo_uncorrelated, 0, atol=1e-8)

        # check that the off-diagonal elements of the correlated covariance matrix are not all zero
        off_diagonal_kinetic_correlated = correlated_kinetic_covariance - np.diag(np.diag(correlated_kinetic_covariance))
        assert not np.allclose(off_diagonal_kinetic_correlated, 0, atol=1e-8)
        off_diagonal_thermo_correlated = correlated_thermo_covariance - np.diag(np.diag(correlated_thermo_covariance))
        assert not np.allclose(off_diagonal_thermo_correlated, 0, atol=1e-8)

        # check that the correlated covariance matrices are symmetric
        assert np.allclose(correlated_kinetic_covariance, correlated_kinetic_covariance.T, atol=1e-8)
        assert np.allclose(correlated_thermo_covariance, correlated_thermo_covariance.T, atol=1e-8)
        assert np.allclose(Sigma_ww_kinetics, Sigma_ww_kinetics.T, atol=1e-8)
        assert np.allclose(Sigma_ww_thermo, Sigma_ww_thermo.T, atol=1e-8)

        # check that the matrix is positive semi-definite by confirming that all eigenvalues are non-negative
        kinetic_eigenvalues = np.linalg.eigvals(correlated_kinetic_covariance)
        assert np.all(kinetic_eigenvalues >= -1e-8)  # allow for small numerical errors
        thermo_eigenvalues = np.linalg.eigvals(correlated_thermo_covariance)
        assert np.all(thermo_eigenvalues >= -1e-8)  # allow for small numerical errors
        intermediate_kinetic_eigenvalues = np.linalg.eigvals(Sigma_ww_kinetics)
        assert np.all(intermediate_kinetic_eigenvalues >= -1e-8)
        intermediate_thermo_eigenvalues = np.linalg.eigvals(Sigma_ww_thermo)
        assert np.all(intermediate_thermo_eigenvalues >= -1e-8)

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
