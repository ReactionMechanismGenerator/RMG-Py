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

from rmgpy import settings
from rmgpy.data.rmg import RMGDatabase, database
from rmgpy.rmg.main import RMG
from rmgpy.reaction import Reaction
from rmgpy.rmg.react import react
from rmgpy.rmg.model import *
from rmgpy.data.base import ForbiddenStructures

###################################################


class TestSpecies(unittest.TestCase):
    """
    Contains unit tests of the Species class.
    """
    @classmethod
    def setUpClass(cls):
        """
        A method that is run before each unit test in this class.
        """
        # set-up RMG object
        cls.rmg = RMG()

        # load kinetic database and forbidden structures
        cls.rmg.database = RMGDatabase()
        path = os.path.join(settings['database.directory'])

        # forbidden structure loading
        cls.rmg.database.loadThermo(os.path.join(path, 'thermo'))
        

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

    @classmethod
    def tearDownClass(cls):
        """
        Reset the loaded database
        """
        import rmgpy.data.rmg
        rmgpy.data.rmg.database = None

class TestCoreEdgeReactionModel(unittest.TestCase):
    """
    Contains unit tests of the CoreEdgeReactionModel class.
    """

    @classmethod
    def setUpClass(cls):
        """
        A method that is run before each unit test in this class.
        """
        TESTFAMILY = 'H_Abstraction'

        # set-up RMG object
        rmg = RMG()

        # load kinetic database and forbidden structures
        rmg.database = RMGDatabase()
        path=os.path.join(settings['test_data.directory'], 'testing_database')


        # kinetics family loading
        rmg.database.loadKinetics(os.path.join(path, 'kinetics'),
                                       kineticsFamilies=[TESTFAMILY],
                                       reactionLibraries=[]
                                       )
        #load empty forbidden structures to avoid any dependence on forbidden structures
        #for these tests
        for family in rmg.database.kinetics.families.values():
            family.forbidden = ForbiddenStructures()
        rmg.database.forbiddenStructures = ForbiddenStructures()

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

    def test_checkForExistingReaction_elminates_duplicate(self):
        """
        Test that checkForExistingReaction catches duplicate reactions
        """
        from rmgpy.data.kinetics.family import TemplateReaction
        cerm = CoreEdgeReactionModel()

        # make species' objects
        spcA = Species().fromSMILES('[H]')
        spcB = Species().fromSMILES('C=C[CH2]C')
        spcC = Species().fromSMILES('C=C=CC')
        spcD = Species().fromSMILES('[H][H]')
        spcA.label = '[H]'
        spcB.label = 'C=C[CH2]C'
        spcC.label = 'C=C=CC'
        spcD.label = '[H][H]'
        spcB.generateResonanceIsomers()

        cerm.addSpeciesToCore(spcA)
        cerm.addSpeciesToCore(spcB)
        cerm.addSpeciesToCore(spcC)
        cerm.addSpeciesToCore(spcD)

        reaction_in_model = TemplateReaction(reactants=[spcA,spcB],
                                             products = [spcC, spcD],
                                                family = 'H_Abstraction',
                                                template = ['Csd','H'])
        reaction_in_model.reactants.sort()
        reaction_in_model.products.sort()

        reaction_to_add = TemplateReaction(reactants=[spcA,spcB],
                                             products = [spcC, spcD],
                                                family = 'H_Abstraction',
                                                template = ['Csd','H'])
        cerm.addReactionToCore(reaction_in_model)
        cerm.registerReaction(reaction_in_model)

        found, rxn = cerm.checkForExistingReaction(reaction_to_add)

        self.assertTrue(found, 'checkForExistingReaction failed to identify existing reaction')

    def test_checkForExistingReaction_keeps_different_template_reactions(self):
        """
        Test that checkForExistingReaction keeps reactions with different templates
        """
        from rmgpy.data.kinetics.family import TemplateReaction
        cerm = CoreEdgeReactionModel()

        # make species' objects
        spcA = Species().fromSMILES('[H]')
        spcB = Species().fromSMILES('C=C[CH2]C')
        spcC = Species().fromSMILES('C=C=CC')
        spcD = Species().fromSMILES('[H][H]')
        spcA.label = '[H]'
        spcB.label = 'C=C[CH2]C'
        spcC.label = 'C=C=CC'
        spcD.label = '[H][H]'
        spcB.generateResonanceIsomers()

        cerm.addSpeciesToCore(spcA)
        cerm.addSpeciesToCore(spcB)
        cerm.addSpeciesToCore(spcC)
        cerm.addSpeciesToCore(spcD)
        
        reaction_in_model = TemplateReaction(reactants=[spcA,spcB],
                                             products = [spcC, spcD],
                                                family = 'H_Abstraction',
                                                template = ['Csd','H'])
        reaction_in_model.reactants.sort()
        reaction_in_model.products.sort()

        reaction_to_add = TemplateReaction(reactants=[spcA,spcB],
                                             products = [spcC, spcD],
                                                family = 'H_Abstraction',
                                                template = ['Cs12345','H'])
        cerm.addReactionToCore(reaction_in_model)
        cerm.registerReaction(reaction_in_model)

        found, rxn = cerm.checkForExistingReaction(reaction_to_add)

        self.assertFalse(found, 'checkForExistingReaction failed to identify reactions with different templates')

    def test_checkForExistingReaction_removes_duplicates_in_opposite_directions(self):
        """
        Test that checkForExistingReaction will remove duplicates if reaction.reverse
        would be a duplicate
        """
        from rmgpy.data.kinetics.family import TemplateReaction
        cerm = CoreEdgeReactionModel()

        # make species' objects
        spcA = Species().fromSMILES('[H]')
        spcB = Species().fromSMILES('C=C[CH2]C')
        spcC = Species().fromSMILES('C=C=CC')
        spcD = Species().fromSMILES('[H][H]')
        spcA.label = '[H]'
        spcB.label = 'C=C[CH2]C'
        spcC.label = 'C=C=CC'
        spcD.label = '[H][H]'
        spcB.generateResonanceIsomers()

        cerm.addSpeciesToCore(spcA)
        cerm.addSpeciesToCore(spcB)
        cerm.addSpeciesToCore(spcC)
        cerm.addSpeciesToCore(spcD)

        reaction_in_model = TemplateReaction(reactants=[spcA,spcB],
                                             products = [spcC, spcD],
                                                family = 'H_Abstraction',
                                                template = ['Csd','H'])
        reaction_in_model.reactants.sort()
        reaction_in_model.products.sort()

        reaction_to_add = TemplateReaction(reactants=[spcA,spcB],
                                             products = [spcC, spcD],
                                                family = 'H_Abstraction',
                                                template = ['Csd','H'])
        cerm.addReactionToCore(reaction_in_model)
        cerm.registerReaction(reaction_in_model)

        found, rxn = cerm.checkForExistingReaction(reaction_to_add)

        self.assertTrue(found, 'checkForExistingReaction failed to identify existing reaction')

    def test_checkForExistingReaction_keeps_different_template_reactions(self):
        """
        Test that checkForExistingReaction keeps reactions with different templates
        """
        from rmgpy.data.kinetics.family import TemplateReaction
        cerm = CoreEdgeReactionModel()

        # make species' objects
        spcA = Species().fromSMILES('[H]')
        spcB = Species().fromSMILES('C=C[CH2]C')
        spcC = Species().fromSMILES('C=C=CC')
        spcD = Species().fromSMILES('[H][H]')
        spcA.label = '[H]'
        spcB.label = 'C=C[CH2]C'
        spcC.label = 'C=C=CC'
        spcD.label = '[H][H]'
        spcB.generateResonanceIsomers()

        cerm.addSpeciesToCore(spcA)
        cerm.addSpeciesToCore(spcB)
        cerm.addSpeciesToCore(spcC)
        cerm.addSpeciesToCore(spcD)
        
        reaction_in_model = TemplateReaction(reactants=[spcA,spcB],
                                             products = [spcC, spcD],
                                                family = 'H_Abstraction',
                                                template = ['Csd','H'])
        reaction_in_model.reactants.sort()
        reaction_in_model.products.sort()

        reaction_to_add = TemplateReaction(reactants=[spcA,spcB],
                                             products = [spcC, spcD],
                                                family = 'H_Abstraction',
                                                template = ['Cs12345','H'])
        cerm.addReactionToCore(reaction_in_model)
        cerm.registerReaction(reaction_in_model)

        found, rxn = cerm.checkForExistingReaction(reaction_to_add)

        self.assertFalse(found, 'checkForExistingReaction failed to identify reactions with different templates')

    def test_checkForExistingReaction_removes_duplicates_in_opposite_directions(self):
        """
        Test that checkForExistingReaction will remove duplicates if reaction.reverse
        would be a duplicate
        """
        from rmgpy.data.kinetics.family import TemplateReaction
        cerm = CoreEdgeReactionModel()

        # make species' objects
        s1 = Species().fromSMILES("[H]")
        s2 = Species().fromSMILES("CC")
        s3 = Species().fromSMILES("[H][H]")
        s4 = Species().fromSMILES("C[CH2]")
        s1.label = 'H'
        s2.label = 'CC'
        s3.label = 'HH'
        s4.label = 'C[CH2]'

        rxn_f = TemplateReaction(reactants = [s1, s2],
                                    products = [s3, s4],
                                    template = ['C/H3/Cs\H3','H_rad'],
                                    degeneracy = 6,
                                    family = 'H_Abstraction',
                                    reverse = TemplateReaction(reactants = [s3, s4],
                                                            products = [s1, s2],
                                                            template = ['H2', 'C_rad/H2/Cs\\H3'],
                                                            degeneracy = 2,
                                                            family = 'H_Abstraction',
                                                            )
                                    )

        rxn_r = TemplateReaction(reactants = [s3, s4],
                                    products = [s1, s2],
                                    template = ['H2', 'C_rad/H2/Cs\\H3'],
                                    degeneracy = 2,
                                    family = 'H_Abstraction',
                                    reverse = TemplateReaction(reactants = [s1, s2],
                                                            products = [s3, s4],
                                                            template = ['C/H3/Cs\H3','H_rad'],
                                                            degeneracy = 6,
                                                            family = 'H_Abstraction',
                                                            )
                                    )

        rxn_f.reactants.sort()
        rxn_f.products.sort()

        cerm.addReactionToCore(rxn_f)
        cerm.registerReaction(rxn_f)

        reactions = cerm.searchRetrieveReactions(rxn_r)
        self.assertEqual(1, len(reactions),'cerm.searchRetrieveReactions could not identify reverse reaction')

        found, rxn = cerm.checkForExistingReaction(rxn_r)

        self.assertTrue(found, 'checkForExistingReaction failed to identify existing reaction when it is in the reverse direction')
        self.assertEqual(rxn, rxn_f)

    @classmethod
    def tearDownClass(cls):
        """
        Reset the loaded database
        """
        import rmgpy.data.rmg
        rmgpy.data.rmg.database = None


if __name__ == '__main__':
    unittest.main()
