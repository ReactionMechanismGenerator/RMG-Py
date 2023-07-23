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
This module contains unit tests of the :mod:`arkane.main` module.
"""

import logging
import os

import pytest

from arkane import Arkane
from arkane.common import clean_dir

from warnings import warn


@pytest.mark.functional
class TestArkaneExamples:
    """
    Run all of Arkane's examples, and report which one failed
    """

    @classmethod
    def setup_class(cls):
        """A function that is run ONCE before all unit tests in this class."""
        cls.base_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "..", "examples", "arkane")
        cls.test_base_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "arkane", "data")
        cls.failed = []
        cls.example_types = ["species", "reactions", "explorer", "networks", "bac"]

    def test_arkane_examples(self):
        for example_type in self.example_types:
            example_type_path = os.path.join(self.base_path, example_type)
            for example in sorted(os.listdir(example_type_path)):
                if example == "thermo_demo":
                    warn("Skipping thermo_demo test - no longer working.", RuntimeWarning)
                    continue
                path = os.path.join(example_type_path, example)
                arkane = Arkane(input_file=os.path.join(path, "input.py"), output_directory=path)
                arkane.plot = example_type != "bac"  # Don't plot BAC examples because they require a lot of memory
                logging.info("running {}".format(example))
                arkane.execute()
                with open(os.path.join(path, "arkane.log"), "r") as f:
                    log = f.readlines()
                for line in log[::-1]:
                    if "execution terminated" in line:
                        break
                else:
                    self.failed.append([example_type, example])
        error_message = "Arkane example(s) failed: "
        for example_type, example_name in self.failed:
            error_message += "{1} in {0}; ".format(example_name, example_type)
        assert not self.failed, error_message

    def test_arkane_two_parameter_arrhenius_fit(self):
        test_path = os.path.join(self.test_base_path, "two_parameter_arrhenius_fit")
        file_to_remove = ["output.py", "chem.inp", "supporting_information.csv"]
        for file in file_to_remove:
            if os.path.exists(os.path.join(test_path, file)):
                os.remove(os.path.join(test_path, file))
        arkane = Arkane(input_file=os.path.join(test_path, "input.py"), output_directory=test_path)
        arkane.plot = False
        arkane.save_rmg_libraries = False
        arkane.execute()
        with open(os.path.join(test_path, "output.py"), "r") as f:
            output = f.readlines()
        reverse_output = output[::-1]
        for i, line in enumerate(reverse_output):
            if "kinetics fitted using" in line:
                msg_output = line.rstrip()
                n_output = reverse_output[i - 5].split("=")[1].strip().replace(",", "")
                break
        msg_expected = "# kinetics fitted using the two-parameter Arrhenius equation k = A * exp(-Ea/RT)"
        assert msg_output == msg_expected
        n_expected = "0"
        assert n_output == n_expected

    @classmethod
    def teardown_class(cls):
        """A function that is run ONCE after all unit tests in this class."""
        cls.extensions_to_delete = ["pdf", "csv", "txt", "inp"]
        cls.files_to_delete = ["arkane.log", "output.py", "supporting_information.csv"]
        cls.files_to_keep = ["README.txt"]  # files to keep that have extensions marked for deletion
        cls.base_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "..", "examples", "arkane")
        for example_type in cls.example_types:
            example_type_path = os.path.join(cls.base_path, example_type)
            for example in os.listdir(example_type_path):
                # clean working folder from all previous test output
                example_path = os.path.join(example_type_path, example)
                clean_dir(
                    base_dir_path=example_path,
                    files_to_delete=cls.files_to_delete,
                    file_extensions_to_delete=cls.extensions_to_delete,
                    files_to_keep=cls.files_to_keep,
                    sub_dir_to_keep=["r0"],
                )
        cls.test_base_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "arkane", "data")
        tests = ["two_parameter_arrhenius_fit"]
        for test in tests:
            test_path = os.path.join(cls.test_base_path, test)
            clean_dir(
                base_dir_path=test_path,
                files_to_delete=cls.files_to_delete,
                file_extensions_to_delete=cls.extensions_to_delete,
            )
