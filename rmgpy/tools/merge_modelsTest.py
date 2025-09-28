#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2019 Prof. William H. Green (whgreen@mit.edu),           #
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

import unittest
import os
import os.path

from rmgpy.tools.merge_models import get_models_to_merge, combine_models


class MergeModelsTest(unittest.TestCase):

    def test_merge_different_models(self):
        folder = os.path.join(os.getcwd(), 'rmgpy/tools/data/diffmodels')

        chemkin3 = os.path.join(folder, 'chem3.inp')
        speciesDict3 = os.path.join(folder, 'species_dictionary3.txt')

        chemkin2 = os.path.join(folder, 'chem2.inp')
        speciesDict2 = os.path.join(folder, 'species_dictionary2.txt')

        models = get_models_to_merge(((chemkin3, speciesDict3, None), (chemkin2, speciesDict2, None)))
        final_model = combine_models(models)
        species = final_model.species
        reactions = final_model.reactions

        # make sure all species are included
        self.assertEqual(len(species), 15)

        # make sure indexes are not unnecessarily redone
        for s in species:
            if s.label == 'CH2O':
                self.assertEqual(s.index, 150)
            elif s.label == 'CH3':
                self.assertEqual(s.index, -1)
            elif s.label == 'C3H7':
                self.assertEqual(s.index, 14)

        # make sure indexes are redone when there is a conflict
        H_index = False
        for s in species:
            if s.label == 'H':
                if isinstance(H_index, bool):
                    H_index = s.index
                else:
                    # found second matching label, make sure index different
                    self.assertNotEqual(s.index, H_index)
                    break
        else:
            raise Exception("Could not find two species identical labels")

        # make sure reaction rates come from first model
        for r in reactions:
            if len(r.reactants) == 2 and r.reactants[0].label == 'CH3' and\
                                         r.reactants[1].label == 'CH3':
                self.assertAlmostEqual(r.kinetics.A.value_si, 8.260e+9, places=0,
                                       msg="Kinetics did not match from first input model")
