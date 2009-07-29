#!/usr/bin/python
# -*- coding: utf-8 -*-
 

import unittest

import pylab
import numpy
import math
import sys

sys.path.append('../source')
import rmg
from rmg.species import *
from rmg.reaction import *
from rmg.model import *



import os
import shutil




################################################################################

def initializeSimulation1():
    
    speciesA = Species(1, 'A')
    speciesB = Species(2, 'B')
    
    reactionAB = Reaction([speciesA], [speciesB])
    reactionAB.kinetics = [ArrheniusEPKinetics(1.0, 0.0, 0.0, 0.0)]
    
    model = CoreEdgeReactionModel()
    model.addSpeciesToCore(speciesA)
    model.addSpeciesToCore(speciesB)
    model.addReactionToCore(reactionAB)
    model.termination.append(TerminationTime(10.0))
    model.absoluteTolerance = 1e-24
    model.relativeTolerance = 1e-12
    
    system = BatchReactor()
    system.equationOfState = IdealGas()
    system.temperatureModel = TemperatureModel()
    system.temperatureModel.setIsothermal(1000.0)
    system.pressureModel = PressureModel()
    system.pressureModel.setIsobaric(1.0E5)
    
    system.initialConcentration[speciesA] = 1.0
    
    return speciesA, speciesB, reactionAB, model, system

################################################################################

def initializeSimulation2():
    
    speciesA = Species(1, 'A')
    speciesB = Species(2, 'B')
    speciesC = Species(3, 'C')
    
    reactionABC = Reaction([speciesA, speciesB], [speciesC])
    reactionABC.kinetics = [ArrheniusEPKinetics(1.0, 0.0, 0.0, 0.0)]
    
    model = CoreEdgeReactionModel()
    model.addSpeciesToCore(speciesA)
    model.addSpeciesToCore(speciesB)
    model.addSpeciesToCore(speciesC)
    model.addReactionToCore(reactionABC)
    model.termination.append(TerminationTime(10.0))
    model.absoluteTolerance = 1e-24
    model.relativeTolerance = 1e-12
    
    system = BatchReactor()
    system.equationOfState = IdealGas()
    system.temperatureModel = TemperatureModel()
    system.temperatureModel.setIsothermal(1000.0)
    system.pressureModel = PressureModel()
    system.pressureModel.setIsobaric(1.0)
    
    return speciesA, speciesB, speciesC, reactionABC, model, system

################################################################################

def simulate(model, system):

    t, y, valid, species = system.simulate(model)
        
    # Reshape y into a matrix rather than a list of lists
    y0 = numpy.zeros((len(t), len(y[0])), float)
    for i, u in enumerate(y):
        for j, v in enumerate(u):
            y0[i,j] = v
            
    return t, y0
    
def postprocess(t, y, model=None):
    """Make concentration plot and show. 
    
    If passed a model, will use the core species to make a legend"""
    pylab.plot(t[1:], y[1:,3:])
    if model:
        pylab.legend([s.label for s in model.core.species])
    pylab.xlabel('Time (s)')
    pylab.ylabel('Concentration (mol/m^3)')
    pylab.show()
    
    #print t, y

################################################################################

    #reload(cti)
    
class SimulationCheck(unittest.TestCase):
    testfolder = 'canteraloadertest'
    
    def setUp(self):
        """setUp gets called before each test"""
        from rmg import constants
        import rmg.cantera_loader
        import ctml_writer as cti
        
        constants.scratchDir = os.path.join(self.testfolder,'temp')
        if os.path.isdir(constants.scratchDir): shutil.rmtree(constants.scratchDir)
        os.makedirs(constants.scratchDir)

        if cti._species:
            reload(cti)
        if rmg.cantera_loader._species:
            reload(rmg.cantera_loader)
            
        rmg.initializeLog(verbose=20)
        
    def tearDown(self):
        """tearDown gets called after each test"""
        pass

    def testLoadCantera(self):
        import sys, os
        import rmg.cantera_loader 
        filename = 'canteraHXD13.cti'
        filepath = os.path.join(self.testfolder,filename)
        model = rmg.cantera_loader.loadCanteraFile(filepath)
        self.assertTrue(len(model.core.species) == 21)
        self.assertTrue(len(model.core.reactions) == 33)
        
    def testLoadChemkin(self):
        import sys, os
        import rmg.cantera_loader 
        
        filename = 'chemkinHXD13.inp'
        filepath = os.path.join(self.testfolder,filename)
        model = rmg.cantera_loader.loadChemkinFile(filepath)
        self.assertTrue(len(model.core.species) == 21) 
        self.assertTrue(len(model.core.reactions) == 33)
        
    def apply_default_settings(self,model,system):
        """Applies a bunch of default settings to the model and system"""
        model.termination.append(TerminationTime(10.0))
        model.absoluteTolerance = 1e-24
        model.relativeTolerance = 1e-12
        system.equationOfState = IdealGas()
        system.temperatureModel = TemperatureModel()
        system.temperatureModel.setIsothermal(1000.0)
        system.pressureModel = PressureModel()
        system.pressureModel.setIsobaric(1.0E5)

    def testCanteraAtoB(self):
        """
        Simulation one is a simple isomerization reaction A --> B, with the
        thermodynamics designed for an equilibrium of all B. This occurs in an
        isothermal, isobaric, homogeneous batch reactor.
        """
        import sys, os
        import rmg.cantera_loader 
        
        filename = 'canteraA>B.cti'
        filepath = os.path.join(self.testfolder,filename)
        model = rmg.cantera_loader.loadCanteraFile(filepath)   
        system = BatchReactor()
        self.apply_default_settings(model,system)
        speciesA = rmg.cantera_loader._speciesByName['A'].getRmgSpecies()
        system.initialConcentration[speciesA] = 1.0
        
        t, y = simulate(model, system)
        postprocess(t, y, model)
        
        # Check equilibrium
        self.assertTrue(y[-1,3] < 0.0001 * y[0,3])
        self.assertTrue(y[-1,4] > 0.9999 * y[0,3])
        
        # Check kinetics
        for i in range(len(t)):
            if abs(t[i] - 1.0) < 0.0001: # find the timestep(s) close to 1 s
                self.assertAlmostEqual(y[i,3] / y[0,3], math.exp(-1.0), 4)
            # check it all the way along
            self.assertAlmostEqual(y[i,3] / y[0,3], math.exp(-1.0*t[i]), 4)
        
    def testChemkinAtoB(self):
        """
        Simulation one is a simple isomerization reaction A --> B, with the
        thermodynamics designed for an equilibrium of all B. This occurs in an
        isothermal, isobaric, homogeneous batch reactor.
        """
        import sys, os
        import rmg.cantera_loader 
        
        filename = 'chemkinA>B.inp'
        filepath = os.path.join(self.testfolder,filename)
        model = rmg.cantera_loader.loadChemkinFile(filepath)   
        system = BatchReactor()
        self.apply_default_settings(model,system)
        speciesA = rmg.cantera_loader._speciesByName['A'].getRmgSpecies()
        system.initialConcentration[speciesA] = 1.0
        
        t, y = simulate(model, system)
        postprocess(t, y, model)
        
        # Check equilibrium
        self.assertTrue(y[-1,3] < 0.0001 * y[0,3])
        self.assertTrue(y[-1,4] > 0.9999 * y[0,3])
        
        # Check kinetics
        for i in range(len(t)):
            if abs(t[i] - 1.0) < 0.0001: # find the timestep(s) close to 1 s
                self.assertAlmostEqual(y[i,3] / y[0,3], math.exp(-1.0), 4)
            # check it all the way along
            self.assertAlmostEqual(y[i,3] / y[0,3], math.exp(-1.0*t[i]), 4)

    def testChemkinAtoB_rev(self):
        """
        Simulation two is a simple isomerization reaction A <-> B, with the
        thermodynamics designed for an equimolar equilibrium. This occurs in an
        isothermal, isobaric, homogeneous batch reactor.
        """
        import sys, os
        import rmg.cantera_loader 
        filename = 'chemkinA=B.inp'
        filepath = os.path.join(self.testfolder,filename)
        model = rmg.cantera_loader.loadChemkinFile(filepath)           
        system = BatchReactor()
        self.apply_default_settings(model,system)
        
        speciesA = rmg.cantera_loader._speciesByName['A'].getRmgSpecies()
        system.initialConcentration[speciesA] = 1.0
        
        t, y = simulate(model, system)
        postprocess(t, y, model)
        # Check equilibrium
        self.assertAlmostEqual(y[-1,3], y[-1,4], 4)
        
    def testChemkinAtoB_rev_2B(self):
        """
        Simulation three is a simple isomerization reaction A <-> B, with the
        thermodynamics designed for an equilibrium ratio of 1:2 A:B. This occurs
        in an isothermal, isobaric, homogeneous batch reactor.
        """
        import sys, os
        import rmg.cantera_loader 
        filename = 'chemkinA=B_2B.inp'
        filepath = os.path.join(self.testfolder,filename)
        model = rmg.cantera_loader.loadChemkinFile(filepath)           
        system = BatchReactor()
        self.apply_default_settings(model,system)
        
        speciesA = rmg.cantera_loader._speciesByName['A'].getRmgSpecies()
        system.initialConcentration[speciesA] = 1.0
        
        t, y = simulate(model, system)
        postprocess(t, y, model)
        
        # Check equilibrium
        self.assertAlmostEqual(2.0 * y[-1,3], y[-1,4], 4)


    def testCanteraAtoB_againstCantera(self):
        """
        Simulation one is a simple isomerization reaction A --> B, with the
        thermodynamics designed for an equilibrium of all B. This occurs in an
        isothermal, isobaric, homogeneous batch reactor. 
        
        Tests by comparing with Cantera
        """
        
        import sys, os
        import rmg.cantera_loader 
        
        # run it in RMG
        filename = 'canteraA>B.cti'
        filepath = os.path.join(self.testfolder,filename)
        model = rmg.cantera_loader.loadCanteraFile(filepath)   
        system = BatchReactor()
        self.apply_default_settings(model,system)
        
        speciesA = rmg.cantera_loader._speciesByName['A'].getRmgSpecies()
        system.initialConcentration[speciesA] = 1.0
        
        t, y = simulate(model, system)
        postprocess(t, y, model)
        
        # Check equilibrium
        self.assertTrue(y[-1,3] < 0.0001 * y[0,3])
        self.assertTrue(y[-1,4] > 0.9999 * y[0,3])
        
        # Check kinetics
        for i in range(len(t)):
            if abs(t[i] - 1.0) < 0.0001: # find the timestep(s) close to 1 s
                self.assertAlmostEqual(y[i,3] / y[0,3], math.exp(-1.0), 4)
            # check it all the way along
            self.assertAlmostEqual(y[i,3] / y[0,3], math.exp(-1.0*t[i]), 4)
            
        # Now test it against Cantera
            
class DontRunTheseTestsYet:
    def testHexadiene(self):
        """1,3-Hexadiene cantera file in an isobaric isothermal batch reactor"""
        import sys, os
        import rmg.cantera_loader 
        filename = 'canteraHXD13.cti'
        filepath = os.path.join(testfolder,filename)
        model = rmg.cantera_loader.loadCanteraFile(filepath)
        
        model.termination.append(rmg.model.TerminationTime(10.0))
        model.absoluteTolerance = 1e-24
        model.relativeTolerance = 1e-12
    
        Temp = 1350.0 # K
        Pres = 101325.0 # Pa
        system = rmg.model.BatchReactor()
        system.equationOfState = rmg.model.IdealGas()
        system.temperatureModel = rmg.model.TemperatureModel()
        system.temperatureModel.setIsothermal(Temp)
        system.pressureModel = rmg.model.PressureModel()
        system.pressureModel.setIsobaric(Pres)
    
        hxd = rmg.cantera_loader._speciesByName['HXD13(1)'].getRmgSpecies()
        model.termination.append(rmg.model.TerminationConversion(hxd,0.9))
        
        Vol = system.equationOfState.getVolume(T=Temp,P=Pres,N=1.0)
        system.initialConcentration[hxd] = 1.0 / Vol
    
        t, y = simulate(model, system)
        postprocess(t, y)
        
        concentrationDict=dict()
        concentrationList=y[-1,3:]/y[-1,1]
        for i,spec in enumerate(model.core.species):
            concentrationDict[spec]=concentrationList[i]
            
        rxn = model.core.reactions[0]
        print rxn.getRate(Temp,Pres,concentrationDict)
    


        
    def testSimulation2A(self):
        """
        Simulation one is an association reaction A + B --> C, with the
        thermodynamics designed for an equilibrium of all C. This occurs in an
        isothermal, isobaric, homogeneous batch reactor.
        """
    
        speciesA, speciesB, speciesC, reactionABC, model, system = initializeSimulation2()
        
        speciesD = Species(4, 'D')
        model.addSpeciesToCore(speciesD)
        system.initialConcentration[speciesD] = 0.0
        
        system.initialConcentration[speciesA] = 1.0
        system.initialConcentration[speciesB] = 1.0
    
        speciesA.thermoData = ThermoGAData(0.0, 0.0, [0, 0, 0, 0, 0, 0, 0])
        speciesB.thermoData = ThermoGAData(0.0, 0.0, [0, 0, 0, 0, 0, 0, 0])
        speciesC.thermoData = ThermoGAData(0.0, 0.0, [0, 0, 0, 0, 0, 0, 0])
        
        t, y = simulate(model, system)
        postprocess(t, y)
        
        # Check equilibrium
        #self.assertTrue(y[-1,3] < 0.0001 * y[0,3])
        #self.assertTrue(y[-1,4] > 0.9999 * y[0,3])
        
        # Check kinetics
        #for i in range(len(t)):
        #   if abs(t[i] - 1.0) < 0.0001:
        #       self.assertAlmostEqual(y[i,3] / y[0,3], math.exp(-1.0), 4)
        
    
        
################################################################################

if __name__ == '__main__':
 #  unittest.main()
    suite = unittest.TestLoader().loadTestsFromTestCase(SimulationCheck)
    unittest.TextTestRunner(verbosity=2).run(suite)
    
    