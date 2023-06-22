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

import os.path
import shutil
import unittest

from rmgpy.tools.diffmodels import execute


class DiffModelsTest(unittest.TestCase):

    def test_identical_models(self):
        folder = os.path.join(os.getcwd(), 'rmgpy/tools/data/diffmodels')

        chemkin1 = os.path.join(folder, 'chem1.inp')
        species_dict1 = os.path.join(folder, 'species_dictionary1.txt')

        chemkin2 = os.path.join(folder, 'chem2.inp')
        species_dict2 = os.path.join(folder, 'species_dictionary2.txt')

        kwargs = {
            'wd': folder,
        }

        execute(chemkin1, species_dict1, None, chemkin2, species_dict2, None, **kwargs)

        shutil.rmtree(os.path.join(folder, 'species1'))
        shutil.rmtree(os.path.join(folder, 'species2'))
        os.remove(os.path.join(folder, 'diff.html'))

    def test_surface_models(self):
        folder = os.path.join(os.getcwd(), 'rmgpy/tools/data/diffmodels/surf_model')

        chemkin_gas1 = os.path.join(folder, 'chem_gas1.inp')
        chemkin_surf1 = os.path.join(folder, 'chem_surface1.inp')
        species_dict1 = os.path.join(folder, 'species_dictionary1.txt')

        chemkin_gas2 = os.path.join(folder, 'chem_gas2.inp')
        chemkin_surf2 = os.path.join(folder, 'chem_surface2.inp')
        species_dict2 = os.path.join(folder, 'species_dictionary2.txt')

        kwargs = {
            'wd': folder,
            'surface_path1': chemkin_surf1,
            'surface_path2': chemkin_surf2,
        }

        execute(chemkin_gas1, species_dict1, None,
                chemkin_gas2, species_dict2, None, **kwargs)

        shutil.rmtree(os.path.join(folder, 'species1'))
        shutil.rmtree(os.path.join(folder, 'species2'))
        os.remove(os.path.join(folder, 'diff.html'))
