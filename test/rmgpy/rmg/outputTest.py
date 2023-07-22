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
import shutil


from rmgpy.rmg.model import CoreEdgeReactionModel, ReactionModel
from rmgpy.rmg.output import save_output_html
from rmgpy.chemkin import load_chemkin_file


class TestOutput:
    def test_save_output_html(self):
        """
        This example is to test if an HTML file can be generated
        for the provided chemkin model.
        """
        folder = os.path.join(os.path.dirname(__file__), "..", "test_data", "saveOutputHTML")

        chemkin_path = os.path.join(folder, "eg6", "chem_annotated.inp")
        dictionary_path = os.path.join(folder, "eg6", "species_dictionary.txt")

        # load_chemkin_file
        species, reactions = load_chemkin_file(chemkin_path, dictionary_path)

        # convert it into a reaction model:
        core = ReactionModel(species, reactions)
        cerm = CoreEdgeReactionModel(core)

        out = os.path.join(folder, "output.html")
        save_output_html(out, cerm)

        assert os.path.isfile(out)
        os.remove(out)
        shutil.rmtree(os.path.join(folder, "species"))
