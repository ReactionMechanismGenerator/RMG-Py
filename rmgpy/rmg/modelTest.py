#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2018 Prof. William H. Green (whgreen@mit.edu),           #
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

from rmgpy import settings
from rmgpy.data.rmg import RMGDatabase, database
from rmgpy.rmg.main import RMG
from rmgpy.reaction import Reaction
from rmgpy.rmg.react import react
from rmgpy.rmg.model import *
from rmgpy.data.base import ForbiddenStructures
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.data.thermo import *
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
    
    def testAddNewSurfaceObjects(self):
        """
        basic test that surface movement object management works properly
        """
        #create object with ReactionSystem behavior
        class rsys:
            pass
        class item:
            pass
        T = item()
        P = item()
        T.value_si = 1000.0
        P.value_si = 101000.0
        rsys.T = T
        rsys.P = P
        
        cerm = CoreEdgeReactionModel()
        
        spcA = Species().fromSMILES('[OH]')
        spcs = [Species().fromSMILES('CC'), Species().fromSMILES('[CH3]')]
        spcTuples = [(spcA, spc) for spc in spcs]
        
        rxns = list(react(*spcTuples))
        rxns += list(react(*[(spcs[0],spcs[1])]))
        
        for rxn in rxns:
            cerm.makeNewReaction(rxn)
        
        cerm.core.species = [spcA]+spcs
        
        corerxns = []
        edgerxns = []
        edgespcs = set()
        for rxn in rxns:
            if set(rxn.reactants+rxn.products) <= set(cerm.core.species):
                corerxns.append(rxn)
            else:
                edgespcs |= set(cerm.core.species)-set(rxn.reactants+rxn.products)
                edgerxns.append(rxn)
                
        cerm.edge.species += list(edgespcs)
        
        cerm.core.reactions = corerxns
        cerm.edge.reactions = edgerxns
        
        cerm.surface.species = []
        cerm.surface.reactions = []
        
        newSurfaceReactions = [cerm.edge.reactions[0]]
        newSurfaceSpecies = []
        obj = newSurfaceReactions
        
        cerm.addNewSurfaceObjects(obj,newSurfaceSpecies,newSurfaceReactions,rsys)
        
        empty = set()
        
        self.assertEqual(cerm.newSurfaceSpcsAdd,empty)
        self.assertEqual(cerm.newSurfaceSpcsLoss,empty)
        self.assertEqual(cerm.newSurfaceRxnsLoss,empty)
        self.assertEqual(cerm.newSurfaceRxnsAdd,set([cerm.edge.reactions[0]]))
        
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
    
    def testThermoFilterSpecies(self):
        """
        test that thermoFilterSpecies leaves species alone if if toleranceThermoKeepInEdge
        is high and removes them if if toleranceThermoKeepInEdge is low
        """
        
        cerm = CoreEdgeReactionModel()

        spcs = [Species().fromSMILES('[OH]'), 
                Species().fromSMILES('C'),
                Species().fromSMILES('[CH3]'),
                Species().fromSMILES('[CH2]'),
                Species().fromSMILES('O')]
        
        for spc in spcs:
            cerm.makeNewSpecies(spc,label=spc.molecule[0].toSMILES())
            spc.label = spc.molecule[0].toSMILES()
        
        thermoDict = {'[OH]':NASA(polynomials=[NASAPolynomial(coeffs=[3.51457,2.92787e-05,-5.32168e-07,1.0195e-09,-3.85947e-13,3414.25,2.10435], Tmin=(100,'K'), Tmax=(1145.75,'K')), NASAPolynomial(coeffs=[3.07194,0.000604014,-1.39775e-08,-2.13448e-11,2.48067e-15,3579.39,4.578], Tmin=(1145.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.3945,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""OH(D)""", comment="""Thermo library: primaryThermoLibrary"""),
        'C' : NASA(polynomials=[NASAPolynomial(coeffs=[4.20541,-0.00535551,2.51121e-05,-2.1376e-08,5.97513e-12,-10161.9,-0.921259], Tmin=(100,'K'), Tmax=(1084.13,'K')), NASAPolynomial(coeffs=[0.908298,0.011454,-4.57171e-06,8.29185e-10,-5.66309e-14,-9719.99,13.9929], Tmin=(1084.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-84.435,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""CH4""", comment="""Thermo library: primaryThermoLibrary"""),
        '[CH3]' : NASA(polynomials=[NASAPolynomial(coeffs=[3.67359,0.00201095,5.73022e-06,-6.87117e-09,2.54386e-12,16445,1.60456], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.28572,0.0072399,-2.98714e-06,5.95685e-10,-4.67154e-14,16775.6,8.48007], Tmin=(1000,'K'), Tmax=(3500,'K'))], Tmin=(200,'K'), Tmax=(3500,'K'), E0=(136.42,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""CH3""", comment="""Thermo library: GRI-Mech3.0"""),
        '[CH2]': NASA(polynomials=[NASAPolynomial(coeffs=[4.01192,-0.000154978,3.26298e-06,-2.40422e-09,5.69497e-13,45867.7,0.533201], Tmin=(100,'K'), Tmax=(1104.62,'K')), NASAPolynomial(coeffs=[3.14983,0.00296674,-9.76056e-07,1.54115e-10,-9.50338e-15,46058.1,4.77808], Tmin=(1104.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(T)""", comment="""Thermo library: primaryThermoLibrary"""),
        'O' : NASA(polynomials=[NASAPolynomial(coeffs=[4.05764,-0.000787936,2.90877e-06,-1.47519e-09,2.12842e-13,-30281.6,-0.311364], Tmin=(100,'K'), Tmax=(1130.24,'K')), NASAPolynomial(coeffs=[2.84325,0.00275109,-7.81031e-07,1.07244e-10,-5.79392e-15,-29958.6,5.91042], Tmin=(1130.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-251.755,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""H2O""", comment="""Thermo library: primaryThermoLibrary"""),
        }
        
        for spc in spcs[:3]:
            cerm.addSpeciesToCore(spc)
            
        reaction = TemplateReaction(
            reactants = [spcs[0],spcs[2]],
            products = [spcs[-1],spcs[-2]],
            degeneracy = 1,
            reversible = True,
            family = 'H_Abstraction',
        )
        
        cerm.processNewReactions(newReactions=[reaction],newSpecies=[]) #adds CH2 and O to edge
        
        for spc in cerm.core.species+cerm.edge.species:
            spc.thermo = thermoDict[spc.molecule[0].toSMILES()] #assign thermo
            
        cerm.setThermodynamicFilteringParameters(Tmax=300.0, 
                                                 toleranceThermoKeepSpeciesInEdge=1000.0,
                                                 minCoreSizeForPrune=0,
                                                 maximumEdgeSpecies=1,
                                                 reactionSystems=[])
        
        cerm.thermoFilterSpecies(cerm.edge.species) #should not do anythinb because toleranceThermoKeepSpeciesInEdge is high
        
        
        difset = set([x.molecule[0].toSMILES() for x in cerm.edge.species])-set([x.molecule[0].toSMILES() for x in cerm.core.species])
        
        self.assertEquals(len(difset),2) #no change in edge
        
        cerm.setThermodynamicFilteringParameters(Tmax=300.0, 
                                                 toleranceThermoKeepSpeciesInEdge=0.0,
                                                 minCoreSizeForPrune=0,
                                                 maximumEdgeSpecies=1,
                                                 reactionSystems=[])
        
        cerm.thermoFilterSpecies(cerm.edge.species) #should remove stuff since CH2 and O have high thermo
        
        difset = set([x.molecule[0].toSMILES() for x in cerm.edge.species])-set([x.molecule[0].toSMILES() for x in cerm.core.species])
        
        self.assertLess(len(difset),2) #edge is smaller
        
    def testThermoFilterDown(self):
        """
        test that thermoFilterDown with maximumEdgeSpecies = 1 reduces
        the edge to one species
        """
        cerm = CoreEdgeReactionModel()

        spcs = [Species().fromSMILES('[OH]'), 
                Species().fromSMILES('C'),
                Species().fromSMILES('[CH3]'),
                Species().fromSMILES('[CH2]'),
                Species().fromSMILES('O')]
        
        for spc in spcs:
            cerm.makeNewSpecies(spc,label=spc.molecule[0].toSMILES())
            spc.label = spc.molecule[0].toSMILES()
        
        thermoDict = {'[OH]':NASA(polynomials=[NASAPolynomial(coeffs=[3.51457,2.92787e-05,-5.32168e-07,1.0195e-09,-3.85947e-13,3414.25,2.10435], Tmin=(100,'K'), Tmax=(1145.75,'K')), NASAPolynomial(coeffs=[3.07194,0.000604014,-1.39775e-08,-2.13448e-11,2.48067e-15,3579.39,4.578], Tmin=(1145.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.3945,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""OH(D)""", comment="""Thermo library: primaryThermoLibrary"""),
        'C' : NASA(polynomials=[NASAPolynomial(coeffs=[4.20541,-0.00535551,2.51121e-05,-2.1376e-08,5.97513e-12,-10161.9,-0.921259], Tmin=(100,'K'), Tmax=(1084.13,'K')), NASAPolynomial(coeffs=[0.908298,0.011454,-4.57171e-06,8.29185e-10,-5.66309e-14,-9719.99,13.9929], Tmin=(1084.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-84.435,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""CH4""", comment="""Thermo library: primaryThermoLibrary"""),
        '[CH3]' : NASA(polynomials=[NASAPolynomial(coeffs=[3.67359,0.00201095,5.73022e-06,-6.87117e-09,2.54386e-12,16445,1.60456], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.28572,0.0072399,-2.98714e-06,5.95685e-10,-4.67154e-14,16775.6,8.48007], Tmin=(1000,'K'), Tmax=(3500,'K'))], Tmin=(200,'K'), Tmax=(3500,'K'), E0=(136.42,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""CH3""", comment="""Thermo library: GRI-Mech3.0"""),
        '[CH2]': NASA(polynomials=[NASAPolynomial(coeffs=[4.01192,-0.000154978,3.26298e-06,-2.40422e-09,5.69497e-13,45867.7,0.533201], Tmin=(100,'K'), Tmax=(1104.62,'K')), NASAPolynomial(coeffs=[3.14983,0.00296674,-9.76056e-07,1.54115e-10,-9.50338e-15,46058.1,4.77808], Tmin=(1104.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(T)""", comment="""Thermo library: primaryThermoLibrary"""),
        'O' : NASA(polynomials=[NASAPolynomial(coeffs=[4.05764,-0.000787936,2.90877e-06,-1.47519e-09,2.12842e-13,-30281.6,-0.311364], Tmin=(100,'K'), Tmax=(1130.24,'K')), NASAPolynomial(coeffs=[2.84325,0.00275109,-7.81031e-07,1.07244e-10,-5.79392e-15,-29958.6,5.91042], Tmin=(1130.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-251.755,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""H2O""", comment="""Thermo library: primaryThermoLibrary"""),
        }
        
        for spc in spcs[:3]:
            cerm.addSpeciesToCore(spc)
            
        reaction = TemplateReaction(
            reactants = [spcs[0],spcs[2]],
            products = [spcs[-1],spcs[-2]],
            degeneracy = 1,
            reversible = True,
            family = 'H_Abstraction',
        )
        
        cerm.processNewReactions(newReactions=[reaction],newSpecies=[]) #add CH2 and O to edge
        
        for spc in cerm.core.species+cerm.edge.species:
            spc.thermo = thermoDict[spc.molecule[0].toSMILES()] #assign thermo
        
        cerm.setThermodynamicFilteringParameters(Tmax=300.0, 
                                                 toleranceThermoKeepSpeciesInEdge=1000.0,
                                                 minCoreSizeForPrune=0,
                                                 maximumEdgeSpecies=1,
                                                 reactionSystems=[])
           
        difset = set([x.molecule[0].toSMILES() for x in cerm.edge.species])-set([x.molecule[0].toSMILES() for x in cerm.core.species])
        
        self.assertEquals(len(difset),2) #no change because toleranceThermoKeepSpeciesInEdge is high
        
        cerm.thermoFilterDown(maximumEdgeSpecies=1)
        
        difset = set([x.molecule[0].toSMILES() for x in cerm.edge.species])-set([x.molecule[0].toSMILES() for x in cerm.core.species])
        
        self.assertEquals(len(difset),1) #should be one because we thermo filtered down to one edge species
        
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

    def test_checkForExistingReaction_eliminates_identical_reactions(self):
        """
        Test that checkForExistingReaction catches identical reactions.
        """
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
        spcB.generate_resonance_structures()

        cerm.addSpeciesToCore(spcA)
        cerm.addSpeciesToCore(spcB)
        cerm.addSpeciesToCore(spcC)
        cerm.addSpeciesToCore(spcD)

        reaction_in_model = TemplateReaction(reactants=[spcA, spcB],
                                             products=[spcC, spcD],
                                             family='H_Abstraction',
                                             template=['Csd', 'H'])
        reaction_in_model.reactants.sort()
        reaction_in_model.products.sort()

        reaction_to_add = TemplateReaction(reactants=[spcA, spcB],
                                           products=[spcC, spcD],
                                           family='H_Abstraction',
                                           template=['Csd', 'H'])
        cerm.addReactionToCore(reaction_in_model)
        cerm.registerReaction(reaction_in_model)

        found, rxn = cerm.checkForExistingReaction(reaction_to_add)

        self.assertTrue(found, 'checkForExistingReaction failed to identify existing reaction')

    def test_checkForExistingReaction_keeps_identical_reactions_with_duplicate_flag(self):
        """
        Test that checkForExistingReaction keeps reactions with different templates and duplicate=True.
        """
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
        spcB.generate_resonance_structures()

        cerm.addSpeciesToCore(spcA)
        cerm.addSpeciesToCore(spcB)
        cerm.addSpeciesToCore(spcC)
        cerm.addSpeciesToCore(spcD)
        
        reaction_in_model = TemplateReaction(reactants=[spcA, spcB],
                                             products=[spcC, spcD],
                                             family='H_Abstraction',
                                             template=['Csd', 'H'],
                                             duplicate=True)
        reaction_in_model.reactants.sort()
        reaction_in_model.products.sort()

        reaction_to_add = TemplateReaction(reactants=[spcA, spcB],
                                           products=[spcC, spcD],
                                           family='H_Abstraction',
                                           template=['Cs12345', 'H'],
                                           duplicate=True)
        cerm.addReactionToCore(reaction_in_model)
        cerm.registerReaction(reaction_in_model)

        found, rxn = cerm.checkForExistingReaction(reaction_to_add)

        self.assertFalse(found, 'checkForExistingReaction failed to identify duplicate template reactions')

    def test_checkForExistingReaction_eliminates_identical_reactions_without_duplicate_flag(self):
        """
        Test that checkForExistingReaction eliminates reactions with different templates and duplicate=false
        """
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
        spcB.generate_resonance_structures()

        cerm.addSpeciesToCore(spcA)
        cerm.addSpeciesToCore(spcB)
        cerm.addSpeciesToCore(spcC)
        cerm.addSpeciesToCore(spcD)

        reaction_in_model = TemplateReaction(reactants=[spcA, spcB],
                                             products=[spcC, spcD],
                                             family='H_Abstraction',
                                             template=['Csd', 'H'],
                                             duplicate=False)
        reaction_in_model.reactants.sort()
        reaction_in_model.products.sort()

        reaction_to_add = TemplateReaction(reactants=[spcA, spcB],
                                           products=[spcC, spcD],
                                           family='H_Abstraction',
                                           template=['Cs12345', 'H'],
                                           duplicate=False)
        cerm.addReactionToCore(reaction_in_model)
        cerm.registerReaction(reaction_in_model)

        found, rxn = cerm.checkForExistingReaction(reaction_to_add)

        self.assertTrue(found, 'checkForExistingReaction failed to eliminate reactions without duplicate tag')

    def test_checkForExistingReaction_removes_duplicates_in_opposite_directions(self):
        """
        Test that checkForExistingReaction removes duplicate reverse reactions
        """
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

        rxn_f = TemplateReaction(reactants=[s1, s2],
                                 products=[s3, s4],
                                 template=['C/H3/Cs/H3', 'H_rad'],
                                 degeneracy=6,
                                 family='H_Abstraction',
                                 reverse=TemplateReaction(reactants=[s3, s4],
                                                          products=[s1, s2],
                                                          template=['H2', 'C_rad/H2/Cs/H3'],
                                                          degeneracy=2,
                                                          family='H_Abstraction')
                                 )

        rxn_r = TemplateReaction(reactants=[s3, s4],
                                 products=[s1, s2],
                                 template=['H2', 'C_rad/H2/Cs/H3'],
                                 degeneracy=2,
                                 family='H_Abstraction',
                                 reverse=TemplateReaction(reactants=[s1, s2],
                                                          products=[s3, s4],
                                                          template=['C/H3/Cs/H3', 'H_rad'],
                                                          degeneracy=6,
                                                          family='H_Abstraction')
                                 )

        rxn_f.reactants.sort()
        rxn_f.products.sort()

        cerm.addReactionToCore(rxn_f)
        cerm.registerReaction(rxn_f)

        reactions = cerm.searchRetrieveReactions(rxn_r)
        self.assertEqual(1, len(reactions), 'cerm.searchRetrieveReactions could not identify reverse reaction')

        found, rxn = cerm.checkForExistingReaction(rxn_r)

        self.assertTrue(found, 'checkForExistingReaction failed to identify existing reaction in the reverse direction')
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
