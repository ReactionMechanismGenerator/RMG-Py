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

import itertools
import os
import unittest 
import numpy as np

from rmgpy import settings
from rmgpy.data.kinetics import TemplateReaction
from rmgpy.data.rmg import RMGDatabase
from rmgpy.species import Species

from rmgpy.rmg.main import RMG
from rmgpy.rmg.react import react, react_all

###################################################

TESTFAMILIES = ['H_Abstraction', 'R_Recombination', 'Disproportionation', 'R_Addition_MultipleBond']


class TestReact(unittest.TestCase):

    def setUp(self):
        """
        A method that is run before each unit test in this class.
        """
        # set-up RMG object
        self.rmg = RMG()

        # load kinetic database and forbidden structures
        self.rmg.database = RMGDatabase()
        path = os.path.join(settings['database.directory'])

        # forbidden structure loading
        self.rmg.database.loadForbiddenStructures(os.path.join(path, 'forbiddenStructures.py'))
        # kinetics family loading
        self.rmg.database.loadKinetics(os.path.join(path, 'kinetics'),
                                       kineticsFamilies=TESTFAMILIES,
                                       reactionLibraries=[]
                                       )

    def testReact(self):
        """
        Test that the ``react`` function works in serial
        """
        procnum = 1

        spc_a = Species().fromSMILES('[OH]')
        spcs = [Species().fromSMILES('CC'), Species().fromSMILES('[CH3]')]
        spc_tuples = [((spc_a, spc), ['H_Abstraction']) for spc in spcs]

        reaction_list = list(itertools.chain.from_iterable(react(spc_tuples, procnum)))
        self.assertIsNotNone(reaction_list)
        self.assertEqual(len(reaction_list), 3)
        self.assertTrue(all([isinstance(rxn, TemplateReaction) for rxn in reaction_list]))

    def testReactParallel(self):
        """
        Test that the ``react`` function works in parallel using Python multiprocessing
        """
        import rmgpy.rmg.main
        rmgpy.rmg.main.maxproc = 2
        procnum = 2

        spc_a = Species().fromSMILES('[OH]')
        spcs = [Species().fromSMILES('CC'), Species().fromSMILES('[CH3]')]
        spc_tuples = [((spc_a, spc), ['H_Abstraction']) for spc in spcs]

        reaction_list = list(itertools.chain.from_iterable(react(spc_tuples, procnum)))
        self.assertIsNotNone(reaction_list)
        self.assertEqual(len(reaction_list), 3)
        self.assertTrue(all([isinstance(rxn, TemplateReaction) for rxn in reaction_list]))

        # Reset module level maxproc back to default
        rmgpy.rmg.main.maxproc = 1

    def testReactAll(self):
        """
        Test that the ``react_all`` function works in serial
        """
        procnum = 1

        spcs = [
                Species().fromSMILES('C=C'),
                Species().fromSMILES('[CH3]'),
                Species().fromSMILES('[OH]'),
                Species().fromSMILES('CCCCCCCCCCC')
                ]

        n = len(spcs)
        reaction_list, spc_tuples = react_all(spcs, n, np.ones(n), np.ones([n, n]), np.ones([n, n, n]), procnum)
        self.assertIsNotNone(reaction_list)
        self.assertEqual(len(reaction_list), 34)
        self.assertEqual(len(spc_tuples), 34)

        flat_rxn_list = list(itertools.chain.from_iterable(reaction_list))
        self.assertEqual(len(flat_rxn_list), 44)
        self.assertTrue(all([isinstance(rxn, TemplateReaction) for rxn in flat_rxn_list]))

    def testReactAllParallel(self):
        """
        Test that the ``react_all`` function works in parallel using Python multiprocessing
        """
        import rmgpy.rmg.main
        rmgpy.rmg.main.maxproc = 2
        procnum = 2

        spcs = [
                Species().fromSMILES('C=C'),
                Species().fromSMILES('[CH3]'),
                Species().fromSMILES('[OH]'),
                Species().fromSMILES('CCCCCCCCCCC')
                ]

        n = len(spcs)
        reaction_list, spc_tuples = react_all(spcs, n, np.ones(n), np.ones([n, n]), np.ones([n, n, n]), procnum)
        self.assertIsNotNone(reaction_list)
        self.assertEqual(len(reaction_list), 94)
        self.assertEqual(len(spc_tuples), 94)

        flat_rxn_list = list(itertools.chain.from_iterable(reaction_list))
        self.assertEqual(len(flat_rxn_list), 44)
        self.assertTrue(all([isinstance(rxn, TemplateReaction) for rxn in flat_rxn_list]))

        # Reset module level maxproc back to default
        rmgpy.rmg.main.maxproc = 1

    def tearDown(self):
        """
        Reset the loaded database
        """
        import rmgpy.data.rmg
        rmgpy.data.rmg.database = None


if __name__ == '__main__':
    unittest.main()
