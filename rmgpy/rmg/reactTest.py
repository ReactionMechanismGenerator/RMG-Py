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

import os
import unittest 
import numpy as np

from rmgpy import settings
from rmgpy.data.kinetics import TemplateReaction
from rmgpy.data.rmg import RMGDatabase
from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.species import Species

from rmgpy.rmg.main import RMG
from rmgpy.rmg.react import react, reactAll

###################################################

TESTFAMILY = 'H_Abstraction'

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
                                       kineticsFamilies=[TESTFAMILY],
                                       reactionLibraries=[]
                                       )

    def testReact(self):
        """
        Test that reaction generation from the available families works.
        """
        spcA = Species().fromSMILES('[OH]')
        spcs = [Species().fromSMILES('CC'), Species().fromSMILES('[CH3]')]
        spcTuples = [(spcA, spc, ['H_Abstraction']) for spc in spcs]

        reactionList = list(react(*spcTuples))
        self.assertIsNotNone(reactionList)
        self.assertTrue(all([isinstance(rxn, TemplateReaction) for rxn in reactionList]))

    def testReactMultiproc(self):
        """
        Test that reaction generation from the available families works with python multiprocessing.
        """
        import rmgpy.rmg.main
        rmgpy.rmg.main.maxproc = 2

        spcA = Species().fromSMILES('[OH]')
        spcs = [Species().fromSMILES('CC'), Species().fromSMILES('[CH3]')]
        spcTuples = [(spcA, spc, ['H_Abstraction']) for spc in spcs]

        reactionList = list(react(*spcTuples))
        self.assertIsNotNone(reactionList)
        self.assertTrue(all([isinstance(rxn, TemplateReaction) for rxn in reactionList]))

    def testReactAll(self):
        """
        Test that the reactAll function works.
        """

        spcs = [
                Species().fromSMILES('CC'),
                Species().fromSMILES('[CH3]'),
                Species().fromSMILES('[OH]')
                ]

        N = len(spcs)
        rxns = reactAll(spcs, N, np.ones(N), np.ones([N,N]))
        self.assertIsNotNone(rxns)
        self.assertTrue(all([isinstance(rxn, TemplateReaction) for rxn in rxns]))


    def tearDown(self):
        """
        Reset the loaded database
        """
        import rmgpy.data.rmg
        rmgpy.data.rmg.database = None

if __name__ == '__main__':
    unittest.main()
