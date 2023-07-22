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
import os
import os.path
import shutil


import rmgpy
from rmgpy.tools.simulate import run_simulation


class SimulateTest:
    def setup_class(self):
        """This method is run once before each unit test"""
        # Disable logging
        logging.disable(logging.WARNING)

    def test_minimal(self):
        """Test that we can simlulate a SimpleReactor with sensitivity"""
        folder = os.path.join(os.path.dirname(rmgpy.__file__), "tools", "data", "sim", "simple")

        input_file = os.path.join(folder, "input.py")
        chemkin_file = os.path.join(folder, "chem.inp")
        dict_file = os.path.join(folder, "species_dictionary.txt")

        run_simulation(input_file, chemkin_file, dict_file)

        simfile = os.path.join(folder, "solver", "simulation_1_13.csv")
        sensfile = os.path.join(folder, "solver", "sensitivity_1_SPC_1.csv")

        assert os.path.isfile(simfile)
        assert os.path.isfile(sensfile)

        shutil.rmtree(os.path.join(folder, "solver"))
        os.remove(os.path.join(folder, "simulate.log"))

    def test_liquid(self):
        """Test that we can simulate a LiquidReactor with sensitivity"""
        folder = os.path.join(os.path.dirname(rmgpy.__file__), "tools", "data", "sim", "liquid")

        input_file = os.path.join(folder, "input.py")
        chemkin_file = os.path.join(folder, "chem.inp")
        dict_file = os.path.join(folder, "species_dictionary.txt")

        run_simulation(input_file, chemkin_file, dict_file, diffusion_limited=False)

        simfile = os.path.join(folder, "solver", "simulation_1_28.csv")
        sensfile = os.path.join(folder, "solver", "sensitivity_1_SPC_1.csv")

        assert os.path.isfile(simfile)
        assert os.path.isfile(sensfile)

        shutil.rmtree(os.path.join(folder, "solver"))
        os.remove(os.path.join(folder, "simulate.log"))

    def test_mb_sampled(self):
        """Test that we can simulate an MBSampledReactor"""
        folder = os.path.join(os.path.dirname(rmgpy.__file__), "tools", "data", "sim", "mbSampled")

        input_file = os.path.join(folder, "input.py")
        chemkin_file = os.path.join(folder, "chem.inp")
        dict_file = os.path.join(folder, "species_dictionary.txt")

        run_simulation(input_file, chemkin_file, dict_file)

        simfile = os.path.join(folder, "solver", "simulation_1_30.csv")

        assert os.path.isfile(simfile)

        shutil.rmtree(os.path.join(folder, "solver"))
        os.remove(os.path.join(folder, "simulate.log"))

    def teardown_class(self):
        import rmgpy.data.rmg

        rmgpy.data.rmg.database = None

        # Reset logging
        logging.disable(logging.NOTSET)
