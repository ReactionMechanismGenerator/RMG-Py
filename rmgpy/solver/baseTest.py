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

import unittest
import pickle
import os.path
import numpy
from rmgpy.tools.loader import loadRMGPyJob
import rmgpy
from rmgpy.rmg.settings import ModelSettings, SimulatorSettings
from rmgpy.solver.base import *

class ConcentrationPrinter:
    def __init__(self):
        self.species_names = []
        self.data = []

    def update(self, subject):
        self.data.append((subject.t , subject.coreSpeciesConcentrations))

class ReactionSystemTest(unittest.TestCase):

    def setUp(self):
        self.listener = ConcentrationPrinter()

        folder = os.path.join(os.path.dirname(rmgpy.__file__),'solver/files/listener/')
        inputFile = os.path.join(folder, 'input.py')
        chemkinFile = os.path.join(folder, 'chemkin/chem.inp')
        spc_dict = os.path.join(folder, 'chemkin/species_dictionary.txt')

        self.rmg = loadRMGPyJob(inputFile, chemkinFile, spc_dict, generateImages=False)

    def testSurfaceInitialization(self):
        """
        test that initialize_surface is correctly removing species and reactions when
        they are no longer consistent with the surface (due to other species/reactions moving to the 
        bulk core)
        """
        reactionSystem = self.rmg.reactionSystems[0]
        reactionSystem.attach(self.listener)
        reactionModel = self.rmg.reactionModel
        
        coreSpecies = reactionModel.core.species
        coreReactions = reactionModel.core.reactions
        surfaceSpecies = [coreSpecies[7],coreSpecies[6]] 
        surfaceReactions = [coreReactions[0],coreReactions[2],coreReactions[3]]
        
        reactionSystem.initializeModel(coreSpecies,coreReactions,
                                       reactionModel.edge.species,reactionModel.edge.reactions,surfaceSpecies,surfaceReactions)
        
        self.assertEquals(len(surfaceSpecies),1) #only H should be left
        self.assertEquals(len(surfaceReactions),2) #all the reactions with H should stay
        
    
    def testSurfaceLayeringConstraint(self):
        """
        test that the correct maximum under the surface layering constraint is being
        found
        """
        reactionSystem = self.rmg.reactionSystems[0]
        reactionSystem.attach(self.listener)
        reactionModel = self.rmg.reactionModel
        coreSpecies = reactionModel.core.species
        coreReactions = reactionModel.core.reactions
        
        edgeSpecies = [coreSpecies[6],coreSpecies[7]]
        edgeReactions = coreReactions[1:]
        surfaceSpecies = [coreSpecies[5]] 
        surfaceReactions = [coreReactions[0]]
        coreSpecies = coreSpecies[0:6]+[coreSpecies[8]]
        coreReactions = surfaceReactions[:]
        reactionSystem.numCoreReactions = 1
        reactionSystem.numCoreSpecies = 7
        
        reactionSystem.initializeModel(coreSpecies,coreReactions,
                                       edgeSpecies,edgeReactions,surfaceSpecies,surfaceReactions)
        
        self.assertEquals(len(reactionSystem.surfaceSpeciesIndices),1) #surfaceSpeciesIndices calculated correctly
        self.assertEquals(reactionSystem.surfaceSpeciesIndices[0],5) #surfaceSpeciesIndices calculated correctly

        inds = reactionSystem.getLayeringIndices()

        self.assertEquals(inds[0],1) #worked correctly
        self.assertEquals(inds[1],2)
    
    def testAddReactionsToSurface(self):
        """
        tests that addReactionsToSurface gives the correct surfaceSpecies and surfaceReactions lists after being called
        """
        reactionSystem = self.rmg.reactionSystems[0]
        reactionSystem.attach(self.listener)
        reactionModel = self.rmg.reactionModel
        species = reactionModel.core.species
        reactions = reactionModel.core.reactions
        
        coreSpecies = species[0:6]
        coreReactions = [reactions[0]]
        surfaceSpecies = []
        surfaceReactions = []
        edgeSpecies = species[6:]
        edgeReactions = reactions[1:]
        
        reactionSystem.initializeModel(coreSpecies,coreReactions,
                                       edgeSpecies,edgeReactions,surfaceSpecies,surfaceReactions)
        
        newSurfaceReactions = edgeReactions
        newSurfaceReactionInds = [edgeReactions.index(i) for i in newSurfaceReactions]
        
        surfaceSpecies,surfaceReactions=reactionSystem.addReactionsToSurface(newSurfaceReactions,newSurfaceReactionInds,surfaceSpecies,surfaceReactions,edgeSpecies)
        
        self.assertEqual(set(surfaceSpecies),set(edgeSpecies)) #all edge species should now be in the surface
        self.assertEqual(set(surfaceReactions),set(edgeReactions)) #all edge reactions should now be in the surface
    
    def testGenerateWeights(self):
        """
        test that proper weight vectors are generated
        """
        reactionSystem = self.rmg.reactionSystems[0]
        reactionSystem.attach(self.listener)
        reactionModel = self.rmg.reactionModel
        species = reactionModel.core.species
        reactions = reactionModel.core.reactions
        
        coreSpecies = species
        coreReactions = reactions
        surfaceSpecies = []
        surfaceReactions = []
        edgeSpecies = []
        edgeReactions = []
        
        reactionSystem.modelSettings = ModelSettings()
        
        reactionSystem.modelSettings.speciesWeights = {species[0].label:2.0}
        reactionSystem.modelSettings.radicalChangeWeightBase = 2.0
        
        reactionSystem.initializeModel(coreSpecies,coreReactions,
                                       edgeSpecies,edgeReactions,surfaceSpecies,surfaceReactions)
        
        spcWts = numpy.ones(reactionSystem.numCoreSpecies)
        spcWts[0] = 2.0
        
        for ind,rxn in enumerate(coreReactions):
            rdif = 0
            for reactant in rxn.reactants:
                rdif += reactant.getRadicalCount()
            for prod in rxn.products:
                rdif -= prod.getRadicalCount()
            rdif = abs(rdif)
            self.assertEqual(2.0**rdif,reactionSystem.reactionWeightVec[ind])
        
        self.assertSequenceEqual(spcWts.tolist(),reactionSystem.speciesWeightVec.tolist())
        
        
    def testEdgeDynamicsNumberCalculations(self):
        """
        test an edge dynamics number calculation
        """
        reactionSystem = self.rmg.reactionSystems[0]
        reactionSystem.attach(self.listener)
        reactionModel = self.rmg.reactionModel
        species = reactionModel.core.species
        reactions = reactionModel.core.reactions
        
        coreSpecies = species[0:6]
        coreReactions = [reactions[0]]
        surfaceSpecies = []
        surfaceReactions = []
        edgeSpecies = species[6:]
        edgeReactions = reactions[1:]
        
        reactionSystem.initializeModel(coreSpecies,coreReactions,
                                       edgeSpecies,edgeReactions,surfaceSpecies,surfaceReactions)
        
        reactionSystem.modelSettings = ModelSettings()
        reactionSystem.simulatorSettings = SimulatorSettings()
        
        reactionSystem.step(1.0e-15)
        
        logging.info('core spcs was {0}'.format(len(reactionSystem.coreSpeciesProductionRates)))
        logging.info(reactionSystem.reactantIndices)
        logging.info(reactionSystem.productIndices)
        
        dynNums = reactionSystem.calculateEdgeDynamicsNumbers()
        
        #the selected reaction (first edge reaction) has 4 and 5 from the core as reactants
        #and 6 and 8 for products
        
        rate = reactionSystem.edgeReactionRates[0]
        
        prod5 = reactionSystem.coreSpeciesProductionRates[5]
        prod4 = reactionSystem.coreSpeciesProductionRates[4]
        loss5 = reactionSystem.coreSpeciesConsumptionRates[5]
        loss4 = reactionSystem.coreSpeciesConsumptionRates[4]
        logging.info((prod5,prod4,loss5,loss4,rate))
        
        lnDyn5 = numpy.log(prod5/loss5)
        lnDyn5r = numpy.log(prod5/(loss5+rate))
        
        lnDyn4 = numpy.log(prod4/loss4)
        lnDyn4r = numpy.log(prod4/(rate+loss4))
        logging.info((lnDyn5,lnDyn5r,lnDyn4,lnDyn4r))
        
        Dyn = abs(lnDyn5-lnDyn5r) + abs(lnDyn4-lnDyn4r)
        
        self.assertAlmostEqual(Dyn,dynNums[0],5)
    
    def testSurfaceDynamicsNumberCalculations(self):
        """
        test a surface dynamics number calculation
        """
        reactionSystem = self.rmg.reactionSystems[0]
        reactionSystem.attach(self.listener)
        reactionModel = self.rmg.reactionModel
        species = reactionModel.core.species
        reactions = reactionModel.core.reactions
        
        coreSpecies = species
        coreReactions = reactions
        surfaceSpecies = species[5:]
        surfaceReactions = reactions
        edgeSpecies = []
        edgeReactions = []
        
        reactionSystem.initializeModel(coreSpecies,coreReactions,
                                       edgeSpecies,edgeReactions,surfaceSpecies,surfaceReactions)
        
        reactionSystem.modelSettings = ModelSettings()
        reactionSystem.simulatorSettings = SimulatorSettings()
        
        reactionSystem.step(1.0e-15)
        
        dynNums = reactionSystem.calculateSurfaceDynamicsNumbers()
        
        #the selected reaction (first core reaction) has 6 + 7 from the core as reactants
        #and 4 as a product, the rate is negative
        #6 and 7 are in the surface so they are ignored
        
        rate = reactionSystem.coreReactionRates[2]
        prod4 = reactionSystem.coreSpeciesProductionRates[4]
        loss4 = reactionSystem.coreSpeciesConsumptionRates[4]
   
        lnDyn4 = numpy.log(prod4) - numpy.log(loss4)
        lnDyn4r = numpy.log(prod4) - numpy.log(loss4-rate)
        
        Dyn = abs(lnDyn4-lnDyn4r)
        
        self.assertAlmostEqual(Dyn,dynNums[2],5)
        
    def testAttachDetach(self):
        """
        Test that a ReactionSystem listener can be attached/detached.
        """
        #create observable

        reactionSystem = self.rmg.reactionSystems[0]
        reactionSystem.attach(self.listener)
        self.assertNotEqual(reactionSystem._observers, [])
        
        reactionSystem.detach(self.listener)    
        self.assertEquals(reactionSystem._observers, [])

    def testListen(self):
        """
        Test that data can be retrieved from an attached ReactionSystem listener.
        """
        #create observable
        reactionSystem = self.rmg.reactionSystems[0]
        reactionSystem.attach(self.listener)

        reactionModel = self.rmg.reactionModel

        self.assertEqual(self.listener.data, [])
        
        modelSettings = ModelSettings(toleranceMoveToCore=1,toleranceKeepInEdge=0,toleranceInterruptSimulation=1)
        simulatorSettings = SimulatorSettings()
        
        # run simulation:
        terminated, obj,sspcs,srxns = reactionSystem.simulate(
            coreSpecies = reactionModel.core.species,
            coreReactions = reactionModel.core.reactions,
            edgeSpecies = reactionModel.edge.species,
            edgeReactions = reactionModel.edge.reactions,
            surfaceSpecies = [],
            surfaceReactions = [],
            modelSettings = modelSettings,
            simulatorSettings = simulatorSettings,
        ) 

        self.assertNotEqual(self.listener.data, [])

    def testPickle(self):
        """
        Test that a ReactionSystem object can be un/pickled.
        """
        rxnSys1 = self.rmg.reactionSystems[0]
        rxnSys = pickle.loads(pickle.dumps(rxnSys1))

        self.assertIsNotNone(rxnSys)
        self.assertTrue(isinstance(rxnSys, rmgpy.solver.simple.SimpleReactor))
        self.assertEqual(rxnSys.T.value_si, rxnSys1.T.value_si)
        self.assertEqual(rxnSys.P.value_si, rxnSys1.P.value_si)
        self.assertEqual(rxnSys.termination[0].conversion, rxnSys1.termination[0].conversion)
        self.assertEqual(rxnSys.termination[1].time.value_si, rxnSys1.termination[1].time.value_si)
        
        
if __name__ == '__main__':
    unittest.main()
