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
import os.path
import shutil

import pytest

import rmgpy
from rmgpy.tools.fluxdiagram import create_flux_diagram


@pytest.mark.functional
class FluxDiagramTest:
    def test_avi_simple(self):
        folder = os.path.join(os.path.dirname(rmgpy.__file__), "tools", "data", "flux")

        input_file = os.path.join(folder, "input_simple.py")
        chemkin_file = os.path.join(folder, "chemkin", "chem.inp")
        dict_file = os.path.join(folder, "chemkin", "species_dictionary.txt")
        settings = {
            "max_node_count": 50,
            "max_edge_count": 50,
            "concentration_tol": 1e-6,
            "species_rate_tol": 1e-6,
            "max_node_pen_width": 10.0,
            "max_edge_pen_width": 10.0,
            "radius": 2,
            "central_reaction_count": 2,
            "time_step": 10**0.1,
        }
        create_flux_diagram(
            input_file,
            chemkin_file,
            dict_file,
            central_species_list=[1],
            superimpose=True,
            settings=settings,
        )

        outputdir = os.path.join(folder, "flux")
        simfile = os.path.join(outputdir, "1", "flux_diagram.avi")

        speciesdir = os.path.join(folder, "species")

        assert os.path.isfile(simfile)

        shutil.rmtree(outputdir)
        shutil.rmtree(speciesdir)

    def test_avi_liquid(self):
        folder = os.path.join(os.path.dirname(rmgpy.__file__), "tools", "data", "flux")

        input_file = os.path.join(folder, "input_liquid.py")
        chemkin_file = os.path.join(folder, "chemkin", "chem.inp")
        dict_file = os.path.join(folder, "chemkin", "species_dictionary.txt")
        create_flux_diagram(input_file, chemkin_file, dict_file, diffusion_limited=False)

        outputdir = os.path.join(folder, "flux")
        simfile = os.path.join(outputdir, "1", "flux_diagram.avi")

        speciesdir = os.path.join(folder, "species")

        assert os.path.isfile(simfile)

        shutil.rmtree(outputdir)
        shutil.rmtree(speciesdir)

    def teardown_class(self):
        import rmgpy.data.rmg

        rmgpy.data.rmg.database = None
