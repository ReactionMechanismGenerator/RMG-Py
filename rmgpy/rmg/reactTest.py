#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2017 Prof. William H. Green (whgreen@mit.edu), 
#   Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

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
from rmgpy.rmg.react import react, reactAll, deflate, deflateReaction

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
        spcTuples = [(spcA, spc) for spc in spcs]

        reactionList = list(react(*spcTuples))
        self.assertIsNotNone(reactionList)
        self.assertTrue(all([isinstance(rxn, TemplateReaction) for rxn in reactionList]))

    def testDeflate(self):
        """
        Test that reaction deflate function works.
        """
        molA = Species().fromSMILES('[OH]')
        molB = Species().fromSMILES('CC')
        molC = Species().fromSMILES('[CH3]')

        reactants = [molA, molB]

        # both reactants were already part of the core:
        reactantIndices = [1, 2]

        rxn = Reaction(reactants=[molA, molB], products=[molC],
        pairs=[(molA, molC), (molB, molC)])

        deflate([rxn], reactants, reactantIndices)

        for spc, t in zip(rxn.reactants, [int, int]):
            self.assertTrue(isinstance(spc, t))
        self.assertEquals(rxn.reactants, reactantIndices)
        for spc in rxn.products:
            self.assertTrue(isinstance(spc, Species))

        # one of the reactants was not yet part of the core:
        reactantIndices = [-1, 2]

        rxn = Reaction(reactants=[molA, molB], products=[molC],
                pairs=[(molA, molC), (molB, molC)])

        deflate([rxn], reactants, reactantIndices)

        for spc, t in zip(rxn.reactants, [Species, int]):
            self.assertTrue(isinstance(spc, t))
        for spc in rxn.products:
            self.assertTrue(isinstance(spc, Species))

    def testReactStoreIndices(self):
        """
        Test that reaction generation keeps track of the original species indices.
        """

        indices = {'[OH]':1, 'CC':2, '[CH3]':3}

        # make it bidirectional so that we can look-up indices as well:
        revd=dict([reversed(i) for i in indices.items()])
        indices.update(revd)

        spcA = Species(index=indices['[OH]']).fromSMILES('[OH]')
        spcs = [Species(index=indices['CC']).fromSMILES('CC'),
                Species(index=indices['[CH3]']).fromSMILES('[CH3]')]

        spcTuples = [(spcA, spc) for spc in spcs]

        reactionList = list(react(*spcTuples))
        self.assertIsNotNone(reactionList)
        self.assertEquals(len(reactionList), 3)
        for rxn in reactionList:
            for i, reactant in enumerate(rxn.reactants):
                rxn.reactants[i] = Molecule().fromSMILES(indices[reactant])
            self.assertTrue(rxn.isBalanced())

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

    def testDeflateReaction(self):
        """
        Test if the deflateReaction function works.
        """

        molA = Species().fromSMILES('[OH]')
        molB = Species().fromSMILES('CC')
        molC = Species().fromSMILES('[CH3]')

        # both reactants were already part of the core:
        reactantIndices = [1, 2]
        molDict = {molA.molecule[0]: 1, molB.molecule[0]: 2}

        rxn = Reaction(reactants=[molA, molB], products=[molC],
        pairs=[(molA, molC), (molB, molC)])

        deflateReaction(rxn, molDict)

        for spc, t in zip(rxn.reactants, [int, int]):
            self.assertTrue(isinstance(spc, t))
        self.assertEquals(rxn.reactants, reactantIndices)
        for spc in rxn.products:
            self.assertTrue(isinstance(spc, Species))

        # one of the reactants was not yet part of the core:
        reactantIndices = [-1, 2]
        molDict = {molA.molecule[0]: molA, molB.molecule[0]: 2}

        rxn = Reaction(reactants=[molA, molB], products=[molC],
                pairs=[(molA, molC), (molB, molC)])

        deflateReaction(rxn, molDict)

        for spc, t in zip(rxn.reactants, [Species, int]):
            self.assertTrue(isinstance(spc, t), 'Species {} is not of type {}'.format(spc,t))
        for spc in rxn.products:
            self.assertTrue(isinstance(spc, Species))


    def tearDown(self):
        """
        Reset the loaded database
        """
        import rmgpy.data.rmg
        rmgpy.data.rmg.database = None

if __name__ == '__main__':
    unittest.main()
