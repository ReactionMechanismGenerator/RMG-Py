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
from rmgpy.species import Species

from rmgpy.rmg.main import RMG
from rmgpy.rmg.react import react, react_all

###################################################

TESTFAMILY = ['H_Abstraction','R_Recombination','Intra_Disproportionation','Intra_RH_Add_Endocyclic',
        'Singlet_Carbene_Intra_Disproportionation','Intra_ene_reaction','Disproportionation',
        '1,4_Linear_birad_scission','R_Addition_MultipleBond','2+2_cycloaddition_Cd','Diels_alder_addition',
        'Intra_RH_Add_Exocyclic','Intra_Retro_Diels_alder_bicyclic','Intra_2+2_cycloaddition_Cd',
        'Birad_recombination','Intra_Diels_alder_monocyclic','1,4_Cyclic_birad_scission',
        '1,2_Insertion_carbene','1,2_Insertion_CO']

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
                                       kineticsFamilies=TESTFAMILY,
                                       reactionLibraries=[]
                                       )

#    def testReact(self):
#        """
#        Test that reaction generation from the available families works.
#        """
#        spcA = Species().fromSMILES('[OH]')
#        spcs = [Species().fromSMILES('CC'), Species().fromSMILES('[CH3]')]
#        spcTuples = [(spcA, spc, ['H_Abstraction']) for spc in spcs]
#
#        reactionList = list(react(*spcTuples))
#        self.assertIsNotNone(reactionList)
#        self.assertTrue(all([isinstance(rxn, TemplateReaction) for rxn in reactionList]))

    def testReactMultiproc(self):
        """
        Test that reaction generation from the available families works with python multiprocessing.
        """
        import rmgpy.rmg.main
        rmgpy.rmg.main.maxproc = 2

        spcA = Species().fromSMILES('[OH]')
        spcs = [Species().fromSMILES('CC'), Species().fromSMILES('[CH3]')]
        spcTuples = [((spcA, spc), ['H_Abstraction']) for spc in spcs]

        reactionList = list(react(*spcTuples))
        self.assertIsNotNone(reactionList)
        self.assertTrue(all([isinstance(rxn, TemplateReaction) for rxn in reactionList]))

    def testReactAll(self):
        """
        Test that the reactAll function works.
        """
        import rmgpy.rmg.main
        rmgpy.rmg.main.maxproc = 2

        spcs = [
                Species().fromSMILES('CC'),
                Species().fromSMILES('[CH3]'),
                Species().fromSMILES('[OH]'),
                Species().fromSMILES('CCCCCCCCCCC')
                ]

        N = len(spcs)
        rxns = react_all(spcs, N, np.ones(N), np.ones([N, N]), np.ones([N, N, N]))
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
