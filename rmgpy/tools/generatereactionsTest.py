#!/usr/bin/env python3

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

import os.path
import shutil
import unittest

from nose.plugins.attrib import attr

import rmgpy
from rmgpy.rmg.main import RMG
from rmgpy.tools.generatereactions import execute


@attr('functional')
class GenerateReactionsTest(unittest.TestCase):

    def test(self):
        folder = os.path.join(os.path.dirname(rmgpy.__file__), 'tools/data/generate')

        input_file = os.path.join(folder, 'input.py')

        rmg = RMG(input_file=input_file, output_directory=folder)
        rmg = execute(rmg)

        self.assertIsNotNone(rmg)
        self.assertIsNotNone(rmg.reaction_model.output_species_list)
        self.assertIsNotNone(rmg.reaction_model.output_reaction_list)

        shutil.rmtree(os.path.join(folder, 'pdep'))
        os.remove(os.path.join(folder, 'restart_from_seed.py'))

    def test_duplicate_reaction(self):
        """
        Test that the radical addition reaction

        HCJ=O + CH2O = [CH2]OC=O

        present in the reaction library "Methylformate",
        only appears once in the model.

        """

        from rmgpy.reaction import Reaction
        from rmgpy.molecule import Molecule
        folder = os.path.join(os.path.dirname(rmgpy.__file__), 'tools/data/generate/duplicates')

        input_file = os.path.join(folder, 'input.py')

        rmg = RMG(input_file=input_file, output_directory=folder)
        rmg = execute(rmg)

        self.assertIsNotNone(rmg)

        rxn_flagged = Reaction(reactants=[Molecule(smiles='[CH]=O'), Molecule(smiles='C=O')],
                               products=[Molecule(smiles='[CH2]OC=O')])

        count = 0
        for reaction in rmg.reaction_model.core.reactions:
            if reaction.is_isomorphic(rxn_flagged):
                count += 1

        self.assertEquals(count, 1)

        shutil.rmtree(os.path.join(folder, 'pdep'))
        os.remove(os.path.join(folder, 'restart_from_seed.py'))

    def test_library_reaction_enters_core(self):
        """
        Test that a reaction from a Reaction Library enters the core
        right after the initialization step if all the input species are 
        present in that reaction.
        
        The following reaction from the Methylformate library
        
        HCjO + CH2O <=> Fmoml
        
        should appear in the model if HCjO, CH2O and Fmoml are all used as input species
        """
        from rmgpy.reaction import Reaction
        from rmgpy.molecule import Molecule
        folder = os.path.join(os.path.dirname(rmgpy.__file__), 'tools/data/generate/libraryReaction')

        input_file = os.path.join(folder, 'input.py')

        rmg = RMG(input_file=input_file, output_directory=folder)
        rmg = execute(rmg)

        self.assertIsNotNone(rmg)

        # Assert that the flagged reaction occurs
        rxn_flagged = Reaction(reactants=[Molecule(smiles='[CH]=O'), Molecule(smiles='C=O')],
                               products=[Molecule(smiles='[CH2]OC=O')])

        count = 0
        for reaction in rmg.reaction_model.core.reactions:
            if reaction.is_isomorphic(rxn_flagged):
                count += 1

        self.assertEquals(count, 1)

        # Assert that the core only has 1 reaction
        self.assertEquals(len(rmg.reaction_model.core.reactions), 1)
        shutil.rmtree(os.path.join(folder, 'pdep'))
        os.remove(os.path.join(folder, 'restart_from_seed.py'))

    def setUp(self):
        import rmgpy.data.rmg
        rmgpy.data.rmg.database = None

    def tearDown(self):
        """
        Reset the loaded database
        """
        import rmgpy.data.rmg
        rmgpy.data.rmg.database = None
