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

"""
This script contains unit tests of the :mod:`arkane.reference` module.
"""

import os

import shutil

from arkane.encorr.isodesmic import ErrorCancelingSpecies
from arkane.modelchem import LevelOfTheory
from arkane.encorr.reference import (
    ReferenceSpecies,
    ReferenceDataEntry,
    CalculatedDataEntry,
    ReferenceDatabase,
)
from rmgpy.species import Species
from rmgpy.thermo import ThermoData
import pytest


FILE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(os.path.dirname(FILE_DIR), "..", "..", "arkane", "data")
SPECIES_DIR = os.path.join(DATA_DIR, "species")


class TestReferenceSpecies:
    """
    A class for testing that the ReferenceSpecies class functions properly
    """

    @classmethod
    def setup_class(cls):
        cls.methane = Species(smiles="C")
        cls.ethane = Species(smiles="CC")
        cls.propane = Species(smiles="CCC")

        cls.thermo_data = ThermoData(H298=(100.0, "kJ/mol"), S298=(100.0, "J/(mol*K)"))
        cls.xyz_dict = {
            "symbols": ("H", "H"),
            "isotopes": (1, 1),
            "coords": ((0.0, 0.0, 0.0), (0.708, 0.0, 0.0)),
        }
        cls.t1_diagnostic = 0.101

    def test_instantiate_reference_species(self):
        """
        Test that a ReferenceSpecies object can be instantiated with the minimal acceptable input, and throws an error
        if the minimal acceptable input is not given.
        """
        ref_spcs = ReferenceSpecies(species=self.ethane)
        assert ref_spcs.smiles == "CC"
        assert ref_spcs.inchi_key == self.ethane.molecule[0].to_inchi_key()

        ref_from_smiles = ReferenceSpecies(smiles="CCC")
        assert ref_from_smiles.smiles == "CCC"
        assert ref_from_smiles.inchi_key == self.propane.molecule[0].to_inchi_key()

        ref_from_inchi = ReferenceSpecies(inchi="InChI=1S/C2H6/c1-2/h1-2H3")
        assert ref_from_inchi.smiles == "CC"
        assert ref_from_inchi.inchi_key == self.ethane.molecule[0].to_inchi_key()

        ref_from_adj = ReferenceSpecies(
            adjacency_list="1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}\n2 H u0 p0 c0 {1,S}\n" "3 H u0 p0 c0 {1,S}\n4 H u0 p0 c0 {1,S}\n5 H u0 p0 c0 {1,S}\n"
        )
        assert ref_from_adj.smiles == "C"
        assert ref_from_adj.inchi_key == self.methane.molecule[0].to_inchi_key()

        with pytest.raises(ValueError):
            ReferenceSpecies()

    def load_ref_from_yaml(self):
        """
        Test that the example ReferenceSpecies YAML file can be loaded
        """
        ref_spcs = ReferenceSpecies.__new__(ReferenceSpecies)
        ref_spcs.load_yaml(os.path.join(SPECIES_DIR, "reference_species_example.yml"))

        assert ref_spcs.smiles == "C#C[CH2]"
        assert ref_spcs.label == "example_reference_species"
        assert isinstance(ref_spcs.calculated_data, CalculatedDataEntry)

    def test_save_ref_to_yaml(self):
        """
        Test that a ReferenceSpecies object can be saved to a YAML file successfully
        """
        label = "test_reference_species"
        ref_spcs = ReferenceSpecies(species=self.ethane, label=label)
        assert ref_spcs.label == label
        ref_spcs.save_yaml(path=SPECIES_DIR)

        loaded_ref = ReferenceSpecies.__new__(ReferenceSpecies)
        load_path = os.path.join(SPECIES_DIR, f"{label}.yml")
        loaded_ref.load_yaml(path=load_path)

        assert loaded_ref.smiles == "CC"

        # Finally, delete this newly created file
        os.remove(load_path)

    def test_reference_data_entry(self):
        """
        Test that the ReferenceDataEntry class functions properly and enforces the standard for storing data
        """
        data_entry = ReferenceDataEntry(self.thermo_data)
        assert isinstance(data_entry.thermo_data, ThermoData)
        assert data_entry.thermo_data.H298.value_si == 100000.0

        with pytest.raises(ValueError):
            ReferenceDataEntry({"H298": (100.0, "kJ/mol")})

    def test_calculated_data_entry(self):
        """
        Test that the CalculatedDataEntry class functions properly and enforces the standard for storing data
        """
        data_entry = CalculatedDataEntry(self.thermo_data, xyz_dict=self.xyz_dict, t1_diagnostic=self.t1_diagnostic)
        assert data_entry.thermo_data.H298.value_si == 100000.0
        assert isinstance(data_entry.xyz_dict, dict)

        data_entry_minimal = CalculatedDataEntry(self.thermo_data)
        assert isinstance(data_entry_minimal.thermo_data, ThermoData)


class TestReferenceDatabase:
    """
    Test that the ReferenceDatabase class functions properly
    """

    @classmethod
    def setup_class(cls):
        cls.database = ReferenceDatabase()
        cls.database.load()

    def test_load_main_reference_set(self):
        """
        Test that the main reference set can be loaded properly
        """
        assert "main" in self.database.reference_sets
        assert isinstance(self.database.reference_sets["main"][0], ReferenceSpecies)

        # Also test that calling load again appends a new set in the database
        data_dir = os.path.join(DATA_DIR)
        testing_dir = os.path.join(data_dir, "testing_set")
        example_ref_file = os.path.join(data_dir, "species", "reference_species_example.yml")
        spcs_file = os.path.join(testing_dir, "0.yml")
        if os.path.exists(testing_dir):  # Delete the testing directory if it existed previously
            shutil.rmtree(testing_dir)
        os.mkdir(testing_dir)
        shutil.copyfile(example_ref_file, spcs_file)
        self.database.load(paths=[testing_dir])
        assert "main" in self.database.reference_sets
        assert "testing_set" in self.database.reference_sets

        # Finally, remove the testing directory
        shutil.rmtree(testing_dir)

    def test_extract_level_of_theory(self):
        """
        Test that a given level of theory can be extracted from the reference set database
        """
        # Create a quick example database
        ref_data_1 = ReferenceDataEntry(ThermoData(H298=(100, "kJ/mol", "+|-", 2)))
        ref_data_2 = ReferenceDataEntry(ThermoData(H298=(25, "kcal/mol", "+|-", 1)))

        calc_data_1 = CalculatedDataEntry(ThermoData(H298=(110, "kJ/mol")))
        calc_data_2 = CalculatedDataEntry(ThermoData(H298=(120, "kJ/mol")))

        ethane = ReferenceSpecies(
            smiles="CC",
            reference_data={"precise": ref_data_1, "less_precise": ref_data_2},
            calculated_data={
                LevelOfTheory("good_chem"): calc_data_1,
                LevelOfTheory("bad_chem"): calc_data_2,
            },
            preferred_reference="less_precise",
        )

        propane = ReferenceSpecies(
            smiles="CCC",
            reference_data={"precise": ref_data_1, "less_precise": ref_data_2},
            calculated_data={
                LevelOfTheory("good_chem"): calc_data_1,
                LevelOfTheory("bad_chem"): calc_data_2,
            },
        )

        butane = ReferenceSpecies(
            smiles="CCCC",
            reference_data={"precise": ref_data_1, "less_precise": ref_data_2},
            calculated_data={LevelOfTheory("bad_chem"): calc_data_2},
        )

        database = ReferenceDatabase()
        database.reference_sets = {
            "testing_1": [ethane, butane],
            "testing_2": [propane],
        }

        model_chem_list = database.extract_level_of_theory(LevelOfTheory("good_chem"))
        assert len(model_chem_list) == 2
        assert isinstance(model_chem_list[0], ErrorCancelingSpecies)

        for spcs in model_chem_list:
            smiles = spcs.molecule.to_smiles()
            assert smiles not in ["CCCC"]
            assert smiles in ["CC", "CCC"]

            if smiles == "CC":  # Test that `less_precise` is the source since it was set manually as preferred
                assert round(abs(spcs.high_level_hf298.value_si - 25.0 * 4184.0), 7) == 0

            if smiles == "CCC":  # Test that `precise` is the source since it has the lowest uncertainty
                assert round(abs(spcs.high_level_hf298.value_si - 100.0 * 1000.0), 7) == 0

    def test_list_available_chemistry(self):
        """
        Test that a set of available levels of theory can be return for the reference database
        """
        level_of_theory_list = self.database.list_available_chemistry()
        assert LevelOfTheory(method="wb97m-v", basis="def2-tzvpd", software="qchem") in level_of_theory_list

    def test_get_species_from_index(self):
        """
        Test that we can retrieve a list of species with specific indices
        """
        test_indices = [5, 309, 105]
        retrieved_species = self.database.get_species_from_index(test_indices)
        assert len(retrieved_species) == 3
        assert retrieved_species[1].index == 309

    def test_get_species_from_label(self):
        """
        Test that we can retrieve a list of species with specific labels
        """
        test_labels = ["1-Butene", "Acetic acid", "Ethanol"]
        retrieved_species = self.database.get_species_from_label(test_labels)
        assert len(retrieved_species) == 3
        assert retrieved_species[0].label == "1-Butene"
