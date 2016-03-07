#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2010 Prof. William H. Green (whgreen@mit.edu) and the
#   RMG Team (rmg_dev@mit.edu)
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

from rmgpy import settings
from rmgpy.data.rmg import RMGDatabase, database
from rmgpy.rmg.main import RMG
from rmgpy.reaction import Reaction
from rmgpy.rmg.react import react
from rmgpy.rmg.model import *

###################################################


class TestSpecies(unittest.TestCase):
    """
    Contains unit tests of the Species class.
    """

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
        self.rmg.database.loadThermo(os.path.join(path, 'thermo'))
        

    def testGetThermoData(self):
        """
        Test that getThermoData method of Species works.
        """
        spc = Species().fromSMILES('CCC')

        self.assertFalse(spc.thermo)
        spc.getThermoData()
        self.assertTrue(spc.thermo)
        thermo = spc.thermo
        spc.getThermoData()

        self.assertEquals(id(thermo), id(spc.thermo))
        
        spc.thermo = None
        spc.getThermoData()
        self.assertNotEquals(id(thermo), id(spc.thermo))

    def tearDown(self):
        """
        Reset the loaded database
        """
        import rmgpy.data.rmg
        rmgpy.data.rmg.database = None

class TestCoreEdgeReactionModel(unittest.TestCase):
    """
    Contains unit tests of the CoreEdgeReactionModel class.
    """

    def setUp(self):
        """
        A method that is run before each unit test in this class.
        """
        TESTFAMILY = 'H_Abstraction'

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

    def testMakeNewSpecies(self):
        """
        Test that CoreEdgeReactionModel.makeNewSpecies method correctly stores the unique species.
        """

        # adding 3 unique species:
        cerm = CoreEdgeReactionModel()

        spcs = [Species().fromSMILES('[OH]'), 
                Species().fromSMILES('CC'),
                Species().fromSMILES('[CH3]')]

        for spc in spcs:
            cerm.makeNewSpecies(spc)

        self.assertEquals(len(cerm.speciesDict), len(spcs))    
        self.assertEquals(len(cerm.indexSpeciesDict), len(spcs))

        # adding 3 unique, and 1 already existing species:
        cerm = CoreEdgeReactionModel()

        spcs = [Species().fromSMILES('[OH]'), 
                Species().fromSMILES('CC'),
                Species().fromSMILES('[CH3]'),
                Species().fromSMILES('CC')]#duplicate species

        for spc in spcs:
            cerm.makeNewSpecies(spc)

        self.assertEquals(len(cerm.speciesDict), len(spcs) - 1)    
        self.assertEquals(len(cerm.indexSpeciesDict), len(spcs) - 1)

    def testMakeNewReaction(self):
        """
        Test that CoreEdgeReactionModel.makeNewReaction method correctly works.
        """

        spcA = Species().fromSMILES('[OH]')
        spcs = [Species().fromSMILES('CC'), Species().fromSMILES('[CH3]')]
        spcTuples = [(spcA, spc) for spc in spcs]

        rxns = list(react(*spcTuples))

        cerm = CoreEdgeReactionModel()

        for rxn in rxns:
            cerm.makeNewReaction(rxn)

        """
        3 expected H-abstraction reactions:
            OH + CC = H2O + C[CH2]
            OH + [CH3] = H2O + [CH2]
            OH + [CH3] = [O] + C
        """

        # count no. of entries in reactionDict:
        counter = 0
        for fam, v1 in cerm.reactionDict.iteritems():
            for key2, v2 in v1.iteritems():
                for key3, rxnList in v2.iteritems():
                    counter += len(rxnList)

        self.assertEquals(counter, 3)

    def testInflate(self):
        """
        Test that CoreEdgeReactionModel.inflate method correctly works.
        """
        spcA = Species().fromSMILES('[OH]')
        spcs = [Species().fromSMILES('CC'), Species().fromSMILES('[CH3]')]
        spcTuples = [(spcA, spc) for spc in spcs]

        rxns = list(react(*spcTuples))

        cerm = CoreEdgeReactionModel()

        for rxn in rxns:
            cerm.makeNewReaction(rxn)

        """
        3 expected H-abstraction reactions:
            OH + CC = H2O + C[CH2]
            OH + [CH3] = H2O + [CH2]
            OH + [CH3] = [O] + C
        """
        for i, rxn in enumerate(rxns):
            rxns[i] = cerm.inflate(rxn)

        for rxn in rxns:
            self.assertTrue(rxn.isBalanced())


    def tearDown(self):
        """
        Reset the loaded database
        """
        import rmgpy.data.rmg
        rmgpy.data.rmg.database = None


if __name__ == '__main__':
    unittest.main()
