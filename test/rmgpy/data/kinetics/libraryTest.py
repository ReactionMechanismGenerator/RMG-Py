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
import os.path
import shutil


import rmgpy
from rmgpy import settings
from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.data.kinetics.library import LibraryReaction
from rmgpy.kinetics import Arrhenius, Troe, PDepArrhenius
from rmgpy.kinetics.model import PDepKineticsModel


class TestLibrary:
    @classmethod
    def setup_class(cls):
        """
        A function run ONCE before all unit tests in this class.
        """
        # Set up a dummy database
        cls.database = KineticsDatabase()
        cls.database.load_libraries(
            os.path.join(
                settings["test_data.directory"],
                "testing_database",
                "kinetics",
                "libraries",
            ),
            libraries=None,
        )  # this loads all of them
        cls.libraries = cls.database.libraries

    def test_get_library_reactions(self):
        """
        test that get_library_reactions loads reactions correctly
        """
        lib_rxns = self.libraries["GRI-Mech3.0"].get_library_reactions()
        for rxn in lib_rxns:
            assert isinstance(rxn, LibraryReaction)

        lib_rxns = self.libraries["ethane-oxidation"].get_library_reactions()  # should have no direct library reactions
        for rxn in lib_rxns:
            if isinstance(rxn.kinetics, PDepKineticsModel):
                assert isinstance(rxn, LibraryReaction)  # can load pdep as networks yet so load as libraries
            else:
                assert isinstance(rxn, TemplateReaction)  # all reactions are template based

    def test_save_library(self):
        """
        This tests the library.save method by writing a new temporary file and
        loading it and comparing the original and copied reactions
        """
        for library_name in ["ethane-oxidation", "surface-example"]:
            copy_path = os.path.join(
                settings["test_data.directory"],
                "testing_database",
                "kinetics",
                "libraries",
                library_name + "-copy",
            )
            if os.path.exists(copy_path):
                logging.warning(f"Removing existing directory {copy_path}.")
                shutil.rmtree(copy_path)
            os.makedirs(copy_path)
            try:
                self.libraries[library_name].save(os.path.join(copy_path, "reactions.py"))
                self.database.load_libraries(
                    os.path.join(
                        settings["test_data.directory"],
                        "testing_database",
                        "kinetics",
                        "libraries",
                    ),
                    libraries=None,
                )  # this loads all of them, including the new copy
                original_rxns = self.database.libraries[library_name].get_library_reactions()
                copy_rxns = self.database.libraries[library_name + "-copy"].get_library_reactions()
                for i in range(len(original_rxns)):
                    if repr(original_rxns[i]).strip() != repr(copy_rxns[i]).strip():  # if it's representation looks different...
                        assert isinstance(copy_rxns[i], TemplateReaction)  # ...it must be a TemplateReaction
                    # Todo: we need to check more carefully for equivalence,
                    # eg. is the metal attribute saved for surface reactions?
            finally:
                shutil.rmtree(copy_path)

    def test_loading_external_kinetic_library(self):
        """This tests loading a kinetic library which is not in the RMG-database repo"""
        kinetic_lib_in_db_path = os.path.join(settings["database.directory"], "kinetics", "libraries", "NOx2018")
        kinetic_lib_in_test_dir_path = os.path.join(os.path.dirname(__file__), "..", "..", "test_data", "copied_kinetic_lib")
        os.makedirs(kinetic_lib_in_test_dir_path)
        for file_name in ["reactions.py", "dictionary.txt"]:
            shutil.copyfile(
                src=os.path.join(kinetic_lib_in_db_path, file_name),
                dst=os.path.join(kinetic_lib_in_test_dir_path, file_name),
            )
        self.database.load_libraries(path="", libraries=[kinetic_lib_in_test_dir_path])
        assert "copied_kinetic_lib" in list(self.database.libraries.keys())
        shutil.rmtree(kinetic_lib_in_test_dir_path, ignore_errors=True)

    def test_generate_high_p_limit_kinetics(self):
        """
        Test that a :class:Arrhenius kinetics object representing the high pressure limit rate
        is returned from Troe/Lindmann/PDepArrhenius/Chebyshev kinetic classes
        """
        lib_rxns = self.libraries["lib_net"].get_library_reactions()
        for rxn in lib_rxns:
            assert rxn.network_kinetics is None
            logging.debug("Processing reaction {0}".format(rxn))
            success = rxn.generate_high_p_limit_kinetics()
            if (isinstance(rxn.kinetics, PDepArrhenius) and rxn.kinetics.pressures.value_si[-1] < 9000000) or not rxn.is_unimolecular():
                # generate_high_p_limit_kinetics() should return `False` if the reaction is not unimolecular
                # or if it is a PDepArrhenius or Chebyshev with Pmax < 90 bar
                assert not success
            else:
                assert success
                if isinstance(rxn.kinetics, Arrhenius):
                    # If the library reaction is already an Arrhenius expression, network_kinetics isn't generated
                    assert rxn.network_kinetics is None
                else:
                    assert isinstance(rxn.network_kinetics, Arrhenius)
                    if isinstance(rxn.kinetics, Troe):
                        # This block quantitative tests the "H + CH2 <=> CH3" reaction
                        # from the test library test_data/testing_database/kinetics/libraries/lib_net/reactions.py
                        # 1. Check that the T exponent in the modified Arrhenius (the "n") equals to 0
                        assert round(abs(rxn.network_kinetics.n.value_si - 0), 7) == 0
                        # 2. Check that the pre-exponential factor equals to 6e+8 m^3/(mol*s)
                        assert round(abs(int(rxn.network_kinetics.A.value_si) - 6e8), 7) == 0
