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
from rmgpy.chemkin import load_chemkin_file


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

        results = execute(chemkin1, species_dict1, None, chemkin2, species_dict2, None, **kwargs)
        assert len(results) == 6, 'Wrong number of return items'

        common_species = results[0]
        unique_species1 = results[1]
        unique_species2 = results[2]
        common_reactions = results[3]
        unique_reactions1 = results[4]
        unique_reactions2 = results[5]

        assert len(common_species) == 13, 'different number of expected common species'
        assert len(unique_species1) == 0, 'model 1 has unique species'
        assert len(unique_species2) == 0, 'model 2 has unique species'

        assert len(common_reactions) == 20, 'different number of expected common reactions'
        assert len(unique_reactions1) == 0, 'model 1 has unique reactions'
        assert len(unique_reactions2) == 0, 'model 2 has unique reactions'


        # check that it made the diff.html
        assert os.path.exists(os.path.join(folder, 'diff.html'))

        shutil.rmtree(os.path.join(folder, 'species1'))
        shutil.rmtree(os.path.join(folder, 'species2'))
        os.remove(os.path.join(folder, 'diff.html'))

    
    def test_different_models(self):
        folder = os.path.join(os.getcwd(), 'rmgpy/tools/data/diffmodels')

        chemkin1 = os.path.join(folder, 'chem1.inp')
        species_dict1 = os.path.join(folder, 'species_dictionary1.txt')

        chemkin3 = os.path.join(folder, 'chem3.inp')
        species_dict3 = os.path.join(folder, 'species_dictionary3.txt')

        kwargs = {
            'wd': folder,
        }

        results = execute(chemkin1, species_dict1, None, chemkin3, species_dict3, None, **kwargs)
        assert len(results) == 6, 'Wrong number of return items'

        common_species = results[0]
        unique_species1 = results[1]
        unique_species2 = results[2]
        common_reactions = results[3]
        unique_reactions1 = results[4]
        unique_reactions2 = results[5]

        assert len(common_species) != 13, 'different number of expected common species'
        assert len(unique_species1) > 0, 'model 1 has no unique species'
        assert len(unique_species2) > 0, 'model 2 has no unique species'

        assert len(common_reactions) != 20, 'different number of expected common reactions'
        assert len(unique_reactions1) > 0, 'model 1 has no unique reactions'
        # assert len(unique_reactions2) > 0, 'model 2 has no unique reactions'


        # Assert that this gives different results
        shutil.rmtree(os.path.join(folder, 'species1'))
        shutil.rmtree(os.path.join(folder, 'species2'))
        os.remove(os.path.join(folder, 'diff.html'))

    def test_identical_surface_models(self):
        folder = os.path.join(os.getcwd(), 'rmgpy/tools/data/diffmodels')

        chemkin1 = os.path.join(folder, 'chem_annotated-gas1.inp')
        surface1 = os.path.join(folder, 'chem_annotated-surface1.inp')
        species_dict1 = os.path.join(folder, 'surf_species_dictionary1.txt')

        chemkin2 = os.path.join(folder, 'chem_annotated-gas1.inp')
        surface2 = os.path.join(folder, 'chem_annotated-surface1.inp')
        species_dict2 = os.path.join(folder, 'surf_species_dictionary1.txt')
        
        gas_species_list, gas_reaction_list = load_chemkin_file(chemkin1, dictionary_path=species_dict1)
        surface_species_list, surface_reaction_list = load_chemkin_file(surface1, dictionary_path=species_dict1)
        n_species = len(gas_species_list) + len(surface_species_list)
        n_reactions = len(gas_reaction_list) + len(surface_reaction_list)


        kwargs = {
            'wd': folder,
            'surface1': surface1,
            'surface2': surface2
        }

        results = execute(chemkin1, species_dict1, None, chemkin2, species_dict2, None, **kwargs)
        assert len(results) == 6, 'Wrong number of return items'

        common_species = results[0]
        unique_species1 = results[1]
        unique_species2 = results[2]
        common_reactions = results[3]
        unique_reactions1 = results[4]
        unique_reactions2 = results[5]

        assert len(common_species) == n_species, 'different number of expected common species'
        assert len(unique_species1) == 0, 'model 1 has unique species'
        assert len(unique_species2) == 0, 'model 2 has unique species'

        assert len(common_reactions) == n_reactions, 'different number of expected common reactions'
        assert len(unique_reactions1) == 0, 'model 1 has unique reactions'
        assert len(unique_reactions2) == 0, 'model 2 has unique reactions'


        # check that it made the diff.html
        assert os.path.exists(os.path.join(folder, 'diff.html'))

        shutil.rmtree(os.path.join(folder, 'species1'))
        shutil.rmtree(os.path.join(folder, 'species2'))
        os.remove(os.path.join(folder, 'diff.html'))
