#!/usr/bin/python
# -*- coding: utf-8 -*-

import unittest
import numpy
import os

import rmgpy.quantity

from rmgpy.molecule import Molecule
from rmgpy.species import Species
from rmgpy.reaction import Reaction
from rmgpy.kinetics import Arrhenius
from rmgpy.thermo import ThermoData
from rmgpy.solver.simple import SimpleReactor
from rmgpy.solver.base import TerminationTime, TerminationConversion
import rmgpy.constants as constants
from rmgpy.chemkin import loadChemkinFile

################################################################################

class SimpleReactorCheck(unittest.TestCase):

    def testSolve(self):
        """
        Test the simple batch reactor with a simple kinetic model. Here we
        choose a kinetic model consisting of the hydrogen abstraction reaction
        CH4 + C2H5 <=> CH3 + C2H6.
        """
        CH4 = Species(
            molecule=[Molecule().fromSMILES("C")],
            thermo=ThermoData(Tdata=([300,400,500,600,800,1000,1500],"K"), Cpdata=([ 8.615, 9.687,10.963,12.301,14.841,16.976,20.528],"cal/(mol*K)"), H298=(-17.714,"kcal/mol"), S298=(44.472,"cal/(mol*K)"))
            )
        CH3 = Species(
            molecule=[Molecule().fromSMILES("[CH3]")],
            thermo=ThermoData(Tdata=([300,400,500,600,800,1000,1500],"K"), Cpdata=([ 9.397,10.123,10.856,11.571,12.899,14.055,16.195],"cal/(mol*K)"), H298=(  9.357,"kcal/mol"), S298=(45.174,"cal/(mol*K)"))
            )
        C2H6 = Species(
            molecule=[Molecule().fromSMILES("CC")],
            thermo=ThermoData(Tdata=([300,400,500,600,800,1000,1500],"K"), Cpdata=([12.684,15.506,18.326,20.971,25.500,29.016,34.595],"cal/(mol*K)"), H298=(-19.521,"kcal/mol"), S298=(54.799,"cal/(mol*K)"))
            )
        C2H5 = Species(
            molecule=[Molecule().fromSMILES("C[CH2]")],
            thermo=ThermoData(Tdata=([300,400,500,600,800,1000,1500],"K"), Cpdata=([11.635,13.744,16.085,18.246,21.885,24.676,29.107],"cal/(mol*K)"), H298=( 29.496,"kcal/mol"), S298=(56.687,"cal/(mol*K)"))
            )
        

        rxn1 = Reaction(reactants=[C2H6,CH3], products=[C2H5,CH4], kinetics=Arrhenius(A=(686.375*6,'m^3/(mol*s)'), n=4.40721, Ea=(7.82799,'kcal/mol'), T0=(298.15,'K')))

        coreSpecies = [CH4,CH3,C2H6,C2H5]
        edgeSpecies = []
        coreReactions = [rxn1]
        edgeReactions = []

        T = 1000; P = 1.0e5
        rxnSystem = SimpleReactor(T, P, initialMoleFractions={C2H5: 0.1, CH3: 0.1, CH4: 0.4, C2H6: 0.4}, termination=[])

        rxnSystem.initializeModel(coreSpecies, coreReactions, edgeSpecies, edgeReactions)

        tlist = numpy.array([10**(i/10.0) for i in range(-130, -49)], numpy.float64)

        # Integrate to get the solution at each time point
        t = []; y = []; reactionRates = []; speciesRates = []
        for t1 in tlist:
            rxnSystem.advance(t1)
            t.append(rxnSystem.t)
            # You must make a copy of y because it is overwritten by DASSL at
            # each call to advance()
            y.append(rxnSystem.y.copy())
            reactionRates.append(rxnSystem.coreReactionRates.copy())
            speciesRates.append(rxnSystem.coreSpeciesRates.copy())

        # Convert the solution vectors to numpy arrays
        t = numpy.array(t, numpy.float64)
        y = numpy.array(y, numpy.float64)
        reactionRates = numpy.array(reactionRates, numpy.float64)
        speciesRates = numpy.array(speciesRates, numpy.float64)
        V = constants.R * rxnSystem.T.value_si * numpy.sum(y) / rxnSystem.P.value_si       

        # Check that we're computing the species fluxes correctly
        for i in range(t.shape[0]):
            self.assertAlmostEqual(reactionRates[i,0], speciesRates[i,0], delta=1e-6*reactionRates[i,0])
            self.assertAlmostEqual(reactionRates[i,0], -speciesRates[i,1], delta=1e-6*reactionRates[i,0])
            self.assertAlmostEqual(reactionRates[i,0], -speciesRates[i,2], delta=1e-6*reactionRates[i,0])
            self.assertAlmostEqual(reactionRates[i,0], speciesRates[i,3], delta=1e-6*reactionRates[i,0])
        
        # Check that we've reached equilibrium 
        self.assertAlmostEqual(reactionRates[-1,0], 0.0, delta=1e-2)
        
        #######        
        # Unit test for the jacobian function:
        # Solve a reaction system and check if the analytical jacobian matches the finite difference jacobian
        
        H2 = Species(
            molecule=[Molecule().fromSMILES("[H][H]")],
            thermo=ThermoData(Tdata=([300,400,500,600,800,1000,1500],"K"), Cpdata=([6.89,6.97,6.99,7.01,7.08,7.22,7.72],"cal/(mol*K)"), H298=( 0,"kcal/mol"), S298=(31.23,"cal/(mol*K)"))
            )
        
        rxnList = []
        rxnList.append(Reaction(reactants=[C2H6], products=[CH3,CH3], kinetics=Arrhenius(A=(686.375*6,'1/s'), n=4.40721, Ea=(7.82799,'kcal/mol'), T0=(298.15,'K'))))
        rxnList.append(Reaction(reactants=[CH3,CH3], products=[C2H6], kinetics=Arrhenius(A=(686.375*6,'m^3/(mol*s)'), n=4.40721, Ea=(7.82799,'kcal/mol'), T0=(298.15,'K'))))
        
        rxnList.append(Reaction(reactants=[C2H6,CH3], products=[C2H5,CH4], kinetics=Arrhenius(A=(46.375*6,'m^3/(mol*s)'), n=3.40721, Ea=(6.82799,'kcal/mol'), T0=(298.15,'K'))))        
        rxnList.append(Reaction(reactants=[C2H5,CH4], products=[C2H6,CH3], kinetics=Arrhenius(A=(46.375*6,'m^3/(mol*s)'), n=3.40721, Ea=(6.82799,'kcal/mol'), T0=(298.15,'K'))))        
        
        rxnList.append(Reaction(reactants=[C2H5,CH4], products=[CH3,CH3,CH3], kinetics=Arrhenius(A=(246.375*6,'m^3/(mol*s)'), n=1.40721, Ea=(3.82799,'kcal/mol'), T0=(298.15,'K'))))       
        rxnList.append(Reaction(reactants=[CH3,CH3,CH3], products=[C2H5,CH4], kinetics=Arrhenius(A=(246.375*6,'m^6/(mol^2*s)'), n=1.40721, Ea=(3.82799,'kcal/mol'), T0=(298.15,'K'))))#        
        
        rxnList.append(Reaction(reactants=[C2H6,CH3,CH3], products=[C2H5,C2H5,H2], kinetics=Arrhenius(A=(146.375*6,'m^6/(mol^2*s)'), n=2.40721, Ea=(8.82799,'kcal/mol'), T0=(298.15,'K'))))
        rxnList.append(Reaction(reactants=[C2H5,C2H5,H2], products=[C2H6,CH3,CH3], kinetics=Arrhenius(A=(146.375*6,'m^6/(mol^2*s)'), n=2.40721, Ea=(8.82799,'kcal/mol'), T0=(298.15,'K'))))
        
        rxnList.append(Reaction(reactants=[C2H6,C2H6], products=[CH3,CH4,C2H5], kinetics=Arrhenius(A=(1246.375*6,'m^3/(mol*s)'), n=0.40721, Ea=(8.82799,'kcal/mol'), T0=(298.15,'K'))))
        rxnList.append(Reaction(reactants=[CH3,CH4,C2H5], products=[C2H6,C2H6], kinetics=Arrhenius(A=(46.375*6,'m^6/(mol^2*s)'), n=0.10721, Ea=(8.82799,'kcal/mol'), T0=(298.15,'K'))))
        

        for rxn in rxnList:
            coreSpecies = [CH4,CH3,C2H6,C2H5,H2]
            edgeSpecies = []
            coreReactions = [rxn]
            
            rxnSystem0 = SimpleReactor(T,P,initialMoleFractions={CH4:0.2,CH3:0.1,C2H6:0.35,C2H5:0.15, H2:0.2},termination=[])
            rxnSystem0.initializeModel(coreSpecies, coreReactions, edgeSpecies, edgeReactions)
            dydt0 = rxnSystem0.residual(0.0, rxnSystem0.y, numpy.zeros(rxnSystem0.y.shape))[0]
            numCoreSpecies = len(coreSpecies)
            dN = .000001*sum(rxnSystem0.y)
            dN_array = dN*numpy.eye(numCoreSpecies)
            
            dydt = []
            for i in range(numCoreSpecies):
                rxnSystem0.y[i] += dN 
                dydt.append(rxnSystem0.residual(0.0, rxnSystem0.y, numpy.zeros(rxnSystem0.y.shape))[0])
                rxnSystem0.y[i] -= dN  # reset y to original y0
            
            # Let the solver compute the jacobian       
            solverJacobian = rxnSystem0.jacobian(0.0, rxnSystem0.y, dydt0, 0.0)     
            # Compute the jacobian using finite differences
            jacobian = numpy.zeros((numCoreSpecies, numCoreSpecies))
            for i in range(numCoreSpecies):
                for j in range(numCoreSpecies):
                    jacobian[i,j] = (dydt[j][i]-dydt0[i])/dN
                    self.assertAlmostEqual(jacobian[i,j], solverJacobian[i,j], delta=abs(1e-4*jacobian[i,j]))
        
        #print 'Solver jacobian'
        #print solverJacobian
        #print 'Numerical jacobian'
        #print jacobian
        
        ###
        # Unit test for the compute rate derivative
        rxnList = []
        rxnList.append(Reaction(reactants=[C2H6], products=[CH3,CH3], kinetics=Arrhenius(A=(686.375e6,'1/s'), n=4.40721, Ea=(7.82799,'kcal/mol'), T0=(298.15,'K')))) 
        rxnList.append(Reaction(reactants=[C2H6,CH3], products=[C2H5,CH4], kinetics=Arrhenius(A=(46.375*6,'m^3/(mol*s)'), n=3.40721, Ea=(6.82799,'kcal/mol'), T0=(298.15,'K'))))        
        rxnList.append(Reaction(reactants=[C2H6,CH3,CH3], products=[C2H5,C2H5,H2], kinetics=Arrhenius(A=(146.375*6,'m^6/(mol^2*s)'), n=2.40721, Ea=(8.82799,'kcal/mol'), T0=(298.15,'K'))))
        
        
        coreSpecies = [CH4,CH3,C2H6,C2H5,H2]
        edgeSpecies = []
        coreReactions = rxnList
        
        rxnSystem0 = SimpleReactor(T,P,initialMoleFractions={CH4:0.2,CH3:0.1,C2H6:0.35,C2H5:0.15, H2:0.2},termination=[])
        rxnSystem0.initializeModel(coreSpecies, coreReactions, edgeSpecies, edgeReactions)
        dfdt0 = rxnSystem0.residual(0.0, rxnSystem0.y, numpy.zeros(rxnSystem0.y.shape))[0]
        solver_dfdk = rxnSystem0.computeRateDerivative()
        #print 'Solver d(dy/dt)/dk'
        #print solver_dfdk
        
        integrationTime = 1e-8
        rxnSystem0.termination.append(TerminationTime((integrationTime,'s')))
        rxnSystem0.simulate(coreSpecies, coreReactions, [], [], 0, 1, 0)

        y0 = rxnSystem0.y
        
        dfdk = numpy.zeros((numCoreSpecies,len(rxnList)))   # d(dy/dt)/dk
        
        for i in range(len(rxnList)):
            k0 = rxnList[i].getRateCoefficient(T,P)
            rxnList[i].kinetics.A.value_si = rxnList[i].kinetics.A.value_si*(1+1e-3)               
            dk = rxnList[i].getRateCoefficient(T,P) - k0

            rxnSystem = SimpleReactor(T,P,initialMoleFractions={CH4:0.2,CH3:0.1,C2H6:0.35,C2H5:0.15, H2:0.2},termination=[])
            rxnSystem.initializeModel(coreSpecies, coreReactions, edgeSpecies, edgeReactions)

            dfdt = rxnSystem.residual(0.0, rxnSystem.y, numpy.zeros(rxnSystem.y.shape))[0]  
            dfdk[:,i]=(dfdt-dfdt0)/dk          
            
            
            rxnSystem.termination.append(TerminationTime((integrationTime,'s')))
            rxnSystem.simulate(coreSpecies, coreReactions, [], [], 0, 1, 0)
            
            rxnList[i].kinetics.A.value_si = rxnList[i].kinetics.A.value_si/(1+1e-3)  # reset A factor
            
        for i in range(numCoreSpecies):
            for j in range(len(rxnList)):
                self.assertAlmostEqual(dfdk[i,j], solver_dfdk[i,j], delta=abs(1e-3*dfdk[i,j]))
            
        #print 'Numerical d(dy/dt)/dk'    
        #print dfdk
        
        
        
#        # Visualize the simulation results
#        import pylab
#        fig = pylab.figure(figsize=(6,6))
#        pylab.subplot(2,1,1)
#        pylab.semilogx(t, y)
#        pylab.ylabel('Concentration (mol/m$^\\mathdefault{3}$)')
#        pylab.legend(['CH4', 'CH3', 'C2H6', 'C2H5'], loc=4)
#        pylab.subplot(2,1,2)
#        pylab.semilogx(t, speciesRates)
#        pylab.legend(['CH4', 'CH3', 'C2H6', 'C2H5'], loc=4)
#        pylab.xlabel('Time (s)')
#        pylab.ylabel('Rate (mol/m$^\\mathdefault{3}$*s)')
#        fig.subplots_adjust(left=0.12, bottom=0.10, right=0.95, top=0.95, wspace=0.20, hspace=0.35)
#        pylab.show()


    def testColliderModel(self):
        """
        Test the solver's ability to simulate a model with collision efficiencies.
        """
        chemFile = os.path.join(os.path.dirname(__file__),'files','collider_model','chem.inp')
        dictionaryFile = os.path.join(os.path.dirname(__file__),'files','collider_model','species_dictionary.txt')
        speciesList, reactionList = loadChemkinFile(chemFile, dictionaryFile)

        
        smilesDict = {'H':'[H]','HO2':'[O]O','O2':'[O][O]','Ar':'[Ar]','N2':'N#N','CO2':'O=C=O','CH3':'[CH3]','CH4':'C'}
        speciesDict = {}
        for name, smiles in smilesDict.iteritems():
            mol = Molecule(SMILES=smiles)
            for species in speciesList:
                if species.isIsomorphic(mol):
                    speciesDict[name] = species
                    break
                
                
        T = 1000 # K
        P = 10 # Pa            
        initialMoleFractions = {speciesDict['O2']:0.5,
                        speciesDict['H']:0.5,
                        speciesDict['CO2']:1.0, 
                        speciesDict['Ar']:4.0}
        
        # Initialize the model
        rxnSystem = SimpleReactor(T,P,initialMoleFractions=initialMoleFractions,termination=None)
        rxnSystem.initializeModel(speciesList, reactionList, [], [])
        
        # Advance to time = 0.1 s
        rxnSystem.advance(0.1)
        # Compare simulated mole fractions with expected mole fractions from CHEMKIN
        simulatedMoleFracs = rxnSystem.y/numpy.sum(rxnSystem.y)
        expectedMoleFracs = numpy.array([0.6666667,0,0,0,0.1666667,0, 0.08333333,0.08333333,2.466066000000000E-10,0,0,0,0,0])
        for i in range(len(simulatedMoleFracs)):
            self.assertAlmostEqual(simulatedMoleFracs[i],expectedMoleFracs[i])
            
        # Advance to time = 5 s
        rxnSystem.advance(5)
        # Compare simulated mole fractions with expected mole fractions from CHEMKIN
        simulatedMoleFracs = rxnSystem.y/numpy.sum(rxnSystem.y)
        expectedMoleFracs = numpy.array([0.6666667,0,0,0, 0.1666667,0,0.08333332,0.08333332,1.233033000000000E-08,0,0,0,0,0])
        for i in range(len(simulatedMoleFracs)):
            self.assertAlmostEqual(simulatedMoleFracs[i],expectedMoleFracs[i])
        
        # Try a new set of conditions
        
        T = 850 # K
        P = 200 # Pa        
        initialMoleFractions = {speciesDict['O2']:0.5,
                        speciesDict['H']:1,
                        speciesDict['CO2']:1, 
                        speciesDict['N2']:4,
                        speciesDict['CH3']:1}
        
        # Initialize the model
        rxnSystem = SimpleReactor(T,P,initialMoleFractions=initialMoleFractions,termination=None)
        rxnSystem.initializeModel(speciesList, reactionList, [], [])
        
        # Advance to time = 5 s
        rxnSystem.advance(5)  

        # Compare simulated mole fractions with expected mole fractions from CHEMKIN
        simulatedMoleFracs = rxnSystem.y/numpy.sum(rxnSystem.y)
        expectedMoleFracs = numpy.array([0,0,0,0.5487241, 0.137181,0, 0.1083234, 0.0685777, 1.280687000000000E-05,  0,0,0,   0.1083362, 0.02884481])
        for i in range(len(simulatedMoleFracs)):
            self.assertAlmostEqual(simulatedMoleFracs[i],expectedMoleFracs[i])
            