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
    
    defaults = { 'time': 10.0,
                 'T': 1000.0,
                 'P': 1.0E5
               }
    
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
        
        pylab.figure(1)
    
    def tearDown(self):
        """tearDown gets called after each test"""
        pass
    def apply_default_settings(self,model,system):
        """Applies a bunch of default settings to the model and system"""
        model.termination.append(TerminationTime(self.defaults['time']))
        model.absoluteTolerance = 1e-24
        model.relativeTolerance = 1e-12
        system.equationOfState = IdealGas()
        system.temperatureModel = TemperatureModel()
        system.temperatureModel.setIsothermal(self.defaults['T'])
        system.pressureModel = PressureModel()
        system.pressureModel.setIsobaric(self.defaults['P'])
    def test1LoadCantera(self):
        import sys, os
        import rmg.cantera_loader 
        filename = 'canteraHXD13.cti'
        filepath = os.path.join(self.testfolder,filename)
        model = rmg.cantera_loader.loadCanteraFile(filepath)
        self.assertTrue(len(model.core.species) == 21)
        self.assertTrue(len(model.core.reactions) == 33)
    
    def test2LoadChemkin(self):
        import sys, os
        import rmg.cantera_loader 
        
        filename = 'chemkinHXD13.inp'
        filepath = os.path.join(self.testfolder,filename)
        model = rmg.cantera_loader.loadChemkinFile(filepath)
        self.assertTrue(len(model.core.species) == 21) 
        self.assertTrue(len(model.core.reactions) == 33)
    
    def test3CanteraAtoB(self):
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
    
    def test4ChemkinAtoB(self):
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
    
    def test5ChemkinAtoB_rev(self):
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
    
    def test6ChemkinAtoB_rev_2B(self):
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
    
    def test7CanteraAtoB_againstCantera(self):
        """
        Simulation one is a simple isomerization reaction A <-> B, with the
        thermodynamics designed for an equilibrium of equal A and B. This occurs in an
        isothermal, isobaric, homogeneous batch reactor. 
        
        Tests by comparing with Cantera
        """
        pylab.figure(2)
        import sys, os
        import rmg.cantera_loader 
        
        # run it in RMG
        filename = 'canteraA=B_2B.cti'
        filepath = os.path.join(self.testfolder,filename)
        model = rmg.cantera_loader.loadCanteraFile(filepath)   
        system = BatchReactor()
        self.apply_default_settings(model,system)
        
        speciesA = rmg.cantera_loader._speciesByName['A'].getRmgSpecies()
        system.initialConcentration[speciesA] = 1.0 # what units? doesn't match P/V
        
        rmg_t, rmg_y = simulate(model, system)
        postprocess(rmg_t, rmg_y, model)
        
        # Now test it against Cantera
        import Cantera
        import Cantera.Reactor 
        
        # load the mechanism into gas
        gas = Cantera.importPhase(filepath,'chem')
        iA = gas.speciesIndex('A') # identify species A
        
        # set the inital gas conditions
        gas.set(T=self.defaults['T'], P=self.defaults['P'])
        gas.setMoleFractions("A:1")
        
        # create the environment
        gasAir = Cantera.Air()
        gasAir.set(T=self.defaults['T'], P=self.defaults['P'])
        # create a reactor for the batch reactor
        # and a reservoir for the environment
        reactor = Cantera.Reactor.Reactor(gas, volume = 1.0)
        environment = Cantera.Reactor.Reservoir(gasAir)
        # Define a wall between the reactor and the environment, and
        # make it flexible, so that the pressure in the reactor is held
        # at the environment pressure.
        wall = Cantera.Reactor.Wall(reactor,environment)
        wall.set(K = 1.0e12)   # set expansion parameter. dV/dt = KA(P_1 - P_2)
        wall.set(A = 1.0)
        wall.setHeatTransferCoeff(1.0e15) # W/m2/K
        
        # put the reactor into a reactor network
        sim = Cantera.Reactor.ReactorNet([reactor]) 
        sim.setInitialTime(0.0)
        maxtime = self.defaults['time']
        
        step=0
        cantera_t=list()
        cantera_y=list()
#        while sim.time()<maxtime:
#            step+=1
        sim.step(maxtime) # take one free step to clear the max step size in the ode solver
        for time in rmg_t:
            if time: # don't integrate if time==0
                sim.advance(time) # now try to hit the rmg_t spot on!
            cantera_t.append(sim.time())
            output=[gas.pressure(), reactor.volume(), gas.temperature()]
            #output.extend(gas.massFractions())
            output.extend([gas.moleFraction(i) for i in range(gas.nSpecies())])
            cantera_y.append(output)
        cantera_t = numpy.array(cantera_t)
        cantera_y = numpy.array(cantera_y)
        postprocess(cantera_t, cantera_y)
        
        #check the results are the same shape
        self.assert_( cantera_y.shape == rmg_y.shape )
        
        #compare pressures
        for i in rmg_y[:,0]/cantera_y[:,0]:
            self.assertAlmostEqual(i,1.0,3)
        #compare volumes
        for i in rmg_y[:,1]/cantera_y[:,1]:
            self.assertAlmostEqual(i,1.0,3)
        #compare temperatures
        for i in rmg_y[:,2]/cantera_y[:,2]:
            self.assertAlmostEqual(i,1.0,3)      
        
        #compare concentration profiles
        ratio = rmg_y[1:,3:] / cantera_y[1:,3:]
        for i in ratio.reshape(ratio.size,1):
            self.assertAlmostEqual(i,1.0,3)
 
    def test8Chemkin2AtoB_2BagainstCantera(self):
        """
        Simulation one is a simple isomerization reaction A <-> B, with the
        thermodynamics designed for an equilibrium of equal A and B. This occurs in an
        isothermal, isobaric, homogeneous batch reactor. 
        
        Tests by comparing with Cantera
        """
        pylab.figure(3)
        import sys, os
        import rmg.cantera_loader 
        
        # run it in RMG
        #filename = 'cantera2A=B_2B.cti'
        filename = 'chemkin2A=B_2B.inp'
        filepath = os.path.join(self.testfolder,filename)
        model = rmg.cantera_loader.loadChemkinFile(filepath)   
        system = BatchReactor()
        self.apply_default_settings(model,system)
        
        speciesA = rmg.cantera_loader._speciesByName['A'].getRmgSpecies()
        system.initialConcentration[speciesA] = 1.0 # what units? doesn't match P/V
        
        rmg_t, rmg_y = simulate(model, system)
        postprocess(rmg_t, rmg_y, model)
        
        # Now test it against Cantera
        import Cantera
        import Cantera.Reactor 
        
        canterafilename = os.path.splitext(filename)[0] + '.cti'
        canterafilepath = os.path.join(self.testfolder,'temp',canterafilename)
        # load the mechanism into gas
        gas = Cantera.importPhase(canterafilepath,'chem')
        iA = gas.speciesIndex('A') # identify species A
        
        # set the inital gas conditions
        gas.set(T=self.defaults['T'], P=self.defaults['P'])
        gas.setMoleFractions("A:1")
        
        # CANTERA USES THESE SI UNITS:
        # Property SI unit 
        # Temperature K 
        # Pressure Pa 
        # Density kg/m3 
        # Quantity kmol 
        # Concentration kmol/m3 
        # Viscosity Pa-s 
        # Thermal Conductivity W/m-K 
        
        #check some units
        import quantities as pq
        rmg_k = model.core.reactions[0].kinetics[0].getRateConstant(
                                            T=self.defaults['T'] ) # m3/mol/s
        rmg_k *= float(pq.quantity.Quantity(1.0,'m^3/mol/s').simplified)
        cantera_k = gas.fwdRateConstants()[0] # m3/kmol/s
        cantera_k *= float(pq.quantity.Quantity(1.0,'m^3/kmol/s').simplified)
        
        self.assertAlmostEqual( rmg_k, cantera_k, 5 )
        
        # create the environment
        gasAir = Cantera.Air()
        gasAir.set(T=self.defaults['T'], P=self.defaults['P'])
        # create a reactor for the batch reactor
        # and a reservoir for the environment
        reactor = Cantera.Reactor.Reactor(gas, volume = 1.0)
        environment = Cantera.Reactor.Reservoir(gasAir)
        # Define a wall between the reactor and the environment, and
        # make it flexible, so that the pressure in the reactor is held
        # at the environment pressure.
        wall = Cantera.Reactor.Wall(reactor,environment)
        wall.set(K = 1.0e12)   # set expansion parameter. dV/dt = KA(P_1 - P_2)
        wall.set(A = 1.0)
        wall.setHeatTransferCoeff(1.0e15) # W/m2/K
        
        # put the reactor into a reactor network
        sim = Cantera.Reactor.ReactorNet([reactor]) 
        sim.setInitialTime(0.0)
        maxtime = self.defaults['time']
        
        step=0
        cantera_t=list()
        cantera_y=list()
#        while sim.time()<maxtime:
#            step+=1
        sim.step(maxtime) # take one free step to clear the max step size in the ode solver
        for time in rmg_t:
            if time: # don't integrate if time==0
                sim.advance(time) # now try to hit the rmg_t spot on!
            cantera_t.append(sim.time())
            output=[gas.pressure(), reactor.volume(), gas.temperature()]
            #output.extend(gas.massFractions())
            output.extend([gas.moleFraction(i) for i in range(gas.nSpecies())])
            cantera_y.append(output)
        cantera_t = numpy.array(cantera_t)
        cantera_y = numpy.array(cantera_y)
        postprocess(cantera_t, cantera_y)
        
        #check the results are the same shape
        self.assert_( cantera_y.shape == rmg_y.shape )
        
        #compare pressures
        for i in rmg_y[:,0]/cantera_y[:,0]:
            self.assertAlmostEqual(i,1.0,3)
        #compare volumes
        for i in rmg_y[:,1]/cantera_y[:,1]:
            self.assertAlmostEqual(i,1.0,3)
        #compare temperatures
        for i in rmg_y[:,2]/cantera_y[:,2]:
            self.assertAlmostEqual(i,1.0,3)      
        
        #compare concentration profiles
        ratio = rmg_y[1:,3:] / cantera_y[1:,3:]
        for i in ratio.reshape(ratio.size,1):
            self.assertAlmostEqual(i,1.0,3)

            
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
  # this is handy to manually interrupt a test so you can 
  # play with the PDB debugger
    SimulationCheck('test8Chemkin2AtoB_2BagainstCantera').debug()