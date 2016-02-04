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

from .main import RMG
from .model import Species
from rmgpy import settings
from rmgpy.data.kinetics import TemplateReaction
from rmgpy.data.rmg import RMGDatabase, database
from rmgpy.molecule import Molecule

from .model import *

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

    def testReactMolecules(self):
        """
        Test that reaction generation for Molecule objects works.
        """
        
        molecules = [Molecule(SMILES='CC'), Molecule(SMILES='[CH3]')]

        reactionList = reactMolecules(molecules, TESTFAMILY)
        
        self.assertIsNotNone(reactionList)
        self.assertTrue(all([isinstance(rxn, TemplateReaction) for rxn in reactionList]))

    def testReactSpecies(self):
        """
        Test that reaction generation for Species objects works.
        """
        spcs = [Species().fromSMILES('CC'), Species().fromSMILES('[CH3]')]

        reactionList = reactSpecies(spcs, TESTFAMILY)
        self.assertIsNotNone(reactionList)
        self.assertTrue(all([isinstance(rxn, TemplateReaction) for rxn in reactionList]))

    def testReactFamily(self):
        """
        Test that reaction generation from a family works.
        """
        spcA = Species().fromSMILES('[OH]')
        spcs = [Species().fromSMILES('CC'), Species().fromSMILES('[CH3]')]

        reactionList = reactFamily(TESTFAMILY, spcA, spcs)
        self.assertIsNotNone(reactionList)
        self.assertTrue(all([isinstance(rxn, TemplateReaction) for rxn in reactionList]))

    def tearDown(self):
        """
        Reset the loaded database
        """
        import rmgpy.data.rmg
        rmgpy.data.rmg.database = None
