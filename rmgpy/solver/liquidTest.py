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
from rmgpy.solver.liquid import LiquidReactor
from rmgpy.solver.base import TerminationTime, TerminationConversion
import rmgpy.constants as constants
from rmgpy.chemkin import loadChemkinFile

################################################################################

class LiquidReactorCheck(unittest.TestCase):

    def setUp(self):
        """
        Here we choose a kinetic model consisting of the hydrogen abstraction reaction
        CH4 + C2H5 <=> CH3 + C2H6.
        """

        Tlist = [300,400,500,600,800,1000,1500]
        self.CH4 = Species(
            molecule=[Molecule().fromSMILES("C")],
            thermo=ThermoData(Tdata=(Tlist,"K"), Cpdata=([ 8.615, 9.687,10.963,12.301,14.841,16.976,20.528],"cal/(mol*K)"), H298=(-17.714,"kcal/mol"), S298=(44.472,"cal/(mol*K)"))
            )
        self.CH3 = Species(
            molecule=[Molecule().fromSMILES("[CH3]")],
            thermo=ThermoData(Tdata=(Tlist,"K"), Cpdata=([ 9.397,10.123,10.856,11.571,12.899,14.055,16.195],"cal/(mol*K)"), H298=(  9.357,"kcal/mol"), S298=(45.174,"cal/(mol*K)"))
            )
        self.C2H6 = Species(
            molecule=[Molecule().fromSMILES("CC")],
            thermo=ThermoData(Tdata=(Tlist,"K"), Cpdata=([12.684,15.506,18.326,20.971,25.500,29.016,34.595],"cal/(mol*K)"), H298=(-19.521,"kcal/mol"), S298=(54.799,"cal/(mol*K)"))
            )
        self.C2H5 = Species(
            molecule=[Molecule().fromSMILES("C[CH2]")],
            thermo=ThermoData(Tdata=(Tlist,"K"), Cpdata=([11.635,13.744,16.085,18.246,21.885,24.676,29.107],"cal/(mol*K)"), H298=( 29.496,"kcal/mol"), S298=(56.687,"cal/(mol*K)"))
            )

        self.H2 = Species(
            molecule=[Molecule().fromSMILES("[H][H]")],
            thermo=ThermoData(Tdata=(Tlist,"K"), Cpdata=([6.89,6.97,6.99,7.01,7.08,7.22,7.72],"cal/(mol*K)"), H298=( 0,"kcal/mol"), S298=(31.23,"cal/(mol*K)"))
            )
        
        self.T = 1000
        
    def testComputeFlux(self):
        """
        Test the liquid batch reactor with a simple kinetic model. 
        """
        
        rxn1 = Reaction(reactants=[self.C2H6,self.CH3], products=[self.C2H5,self.CH4], kinetics=Arrhenius(A=(686.375*6,'m^3/(mol*s)'), n=4.40721, Ea=(7.82799,'kcal/mol'), T0=(298.15,'K')))

        coreSpecies = [self.CH4,self.CH3,self.C2H6,self.C2H5]
        edgeSpecies = []
        coreReactions = [rxn1]
        edgeReactions = []

        
        c0={self.C2H5: 0.1, self.CH3: 0.1, self.CH4: 0.4, self.C2H6: 0.4}


        rxnSystem = LiquidReactor(self.T, c0, termination=[])

        rxnSystem.initializeModel(coreSpecies, coreReactions, edgeSpecies, edgeReactions)

        tlist = numpy.array([10**(i/10.0) for i in xrange(-130, -49)], numpy.float64)

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
        for i in xrange(t.shape[0]):
            self.assertAlmostEqual(reactionRates[i,0], speciesRates[i,0], delta=1e-6*reactionRates[i,0])
            self.assertAlmostEqual(reactionRates[i,0], -speciesRates[i,1], delta=1e-6*reactionRates[i,0])
            self.assertAlmostEqual(reactionRates[i,0], -speciesRates[i,2], delta=1e-6*reactionRates[i,0])
            self.assertAlmostEqual(reactionRates[i,0], speciesRates[i,3], delta=1e-6*reactionRates[i,0])
        
        # Check that we've reached equilibrium 
        self.assertAlmostEqual(reactionRates[-1,0], 0.0, delta=1e-2)
        

    def test_jacobian(self):
        """
        Unit test for the jacobian function:
        Solve a reaction system and check if the analytical jacobian matches the finite difference jacobian

        """

        coreSpecies = [self.CH4,self.CH3,self.C2H6,self.C2H5]
        edgeSpecies = []

        rxn1 = Reaction(reactants=[self.C2H6,self.CH3], products=[self.C2H5,self.CH4], kinetics=Arrhenius(A=(686.375*6,'m^3/(mol*s)'), n=4.40721, Ea=(7.82799,'kcal/mol'), T0=(298.15,'K')))
        coreReactions = [rxn1]
        edgeReactions = []
        numCoreSpecies = len(coreSpecies)

        rxnList = []
        rxnList.append(Reaction(reactants=[self.C2H6], products=[self.CH3,self.CH3], kinetics=Arrhenius(A=(686.375*6,'1/s'), n=4.40721, Ea=(7.82799,'kcal/mol'), T0=(298.15,'K'))))
        rxnList.append(Reaction(reactants=[self.CH3,self.CH3], products=[self.C2H6], kinetics=Arrhenius(A=(686.375*6,'m^3/(mol*s)'), n=4.40721, Ea=(7.82799,'kcal/mol'), T0=(298.15,'K'))))
        
        rxnList.append(Reaction(reactants=[self.C2H6,self.CH3], products=[self.C2H5,self.CH4], kinetics=Arrhenius(A=(46.375*6,'m^3/(mol*s)'), n=3.40721, Ea=(6.82799,'kcal/mol'), T0=(298.15,'K'))))        
        rxnList.append(Reaction(reactants=[self.C2H5,self.CH4], products=[self.C2H6,self.CH3], kinetics=Arrhenius(A=(46.375*6,'m^3/(mol*s)'), n=3.40721, Ea=(6.82799,'kcal/mol'), T0=(298.15,'K'))))        
        
        rxnList.append(Reaction(reactants=[self.C2H5,self.CH4], products=[self.CH3,self.CH3,self.CH3], kinetics=Arrhenius(A=(246.375*6,'m^3/(mol*s)'), n=1.40721, Ea=(3.82799,'kcal/mol'), T0=(298.15,'K'))))       
        rxnList.append(Reaction(reactants=[self.CH3,self.CH3,self.CH3], products=[self.C2H5,self.CH4], kinetics=Arrhenius(A=(246.375*6,'m^6/(mol^2*s)'), n=1.40721, Ea=(3.82799,'kcal/mol'), T0=(298.15,'K'))))#        
        
        rxnList.append(Reaction(reactants=[self.C2H6,self.CH3,self.CH3], products=[self.C2H5,self.C2H5,self.H2], kinetics=Arrhenius(A=(146.375*6,'m^6/(mol^2*s)'), n=2.40721, Ea=(8.82799,'kcal/mol'), T0=(298.15,'K'))))
        rxnList.append(Reaction(reactants=[self.C2H5,self.C2H5,self.H2], products=[self.C2H6,self.CH3,self.CH3], kinetics=Arrhenius(A=(146.375*6,'m^6/(mol^2*s)'), n=2.40721, Ea=(8.82799,'kcal/mol'), T0=(298.15,'K'))))
        
        rxnList.append(Reaction(reactants=[self.C2H6,self.C2H6], products=[self.CH3,self.CH4,self.C2H5], kinetics=Arrhenius(A=(1246.375*6,'m^3/(mol*s)'), n=0.40721, Ea=(8.82799,'kcal/mol'), T0=(298.15,'K'))))
        rxnList.append(Reaction(reactants=[self.CH3,self.CH4,self.C2H5], products=[self.C2H6,self.C2H6], kinetics=Arrhenius(A=(46.375*6,'m^6/(mol^2*s)'), n=0.10721, Ea=(8.82799,'kcal/mol'), T0=(298.15,'K'))))
        

        for rxn in rxnList:
            coreSpecies = [self.CH4,self.CH3,self.C2H6,self.C2H5,self.H2]
            edgeSpecies = []
            coreReactions = [rxn]
            
            c0={self.CH4:0.2,self.CH3:0.1,self.C2H6:0.35,self.C2H5:0.15, self.H2:0.2}
            rxnSystem0 = LiquidReactor(self.T, c0,termination=[])
            rxnSystem0.initializeModel(coreSpecies, coreReactions, edgeSpecies, edgeReactions)
            dydt0 = rxnSystem0.residual(0.0, rxnSystem0.y, numpy.zeros(rxnSystem0.y.shape))[0]
            
            dN = .000001*sum(rxnSystem0.y)
            dN_array = dN*numpy.eye(numCoreSpecies)
            
            dydt = []
            for i in xrange(numCoreSpecies):
                rxnSystem0.y[i] += dN 
                dydt.append(rxnSystem0.residual(0.0, rxnSystem0.y, numpy.zeros(rxnSystem0.y.shape))[0])
                rxnSystem0.y[i] -= dN  # reset y to original y0
            
            # Let the solver compute the jacobian       
            solverJacobian = rxnSystem0.jacobian(0.0, rxnSystem0.y, dydt0, 0.0)     
            # Compute the jacobian using finite differences
            jacobian = numpy.zeros((numCoreSpecies, numCoreSpecies))
            for i in xrange(numCoreSpecies):
                for j in xrange(numCoreSpecies):
                    jacobian[i,j] = (dydt[j][i]-dydt0[i])/dN
                    self.assertAlmostEqual(jacobian[i,j], solverJacobian[i,j], delta=abs(1e-4*jacobian[i,j]))
        
        #print 'Solver jacobian'
        #print solverJacobian
        #print 'Numerical jacobian'
        #print jacobian

     
    def test_compute_derivative(self):

        rxnList = []
        rxnList.append(Reaction(reactants=[self.C2H6], products=[self.CH3,self.CH3], kinetics=Arrhenius(A=(686.375e6,'1/s'), n=4.40721, Ea=(7.82799,'kcal/mol'), T0=(298.15,'K')))) 
        rxnList.append(Reaction(reactants=[self.C2H6,self.CH3], products=[self.C2H5,self.CH4], kinetics=Arrhenius(A=(46.375*6,'m^3/(mol*s)'), n=3.40721, Ea=(6.82799,'kcal/mol'), T0=(298.15,'K'))))        
        rxnList.append(Reaction(reactants=[self.C2H6,self.CH3,self.CH3], products=[self.C2H5,self.C2H5,self.H2], kinetics=Arrhenius(A=(146.375*6,'m^6/(mol^2*s)'), n=2.40721, Ea=(8.82799,'kcal/mol'), T0=(298.15,'K'))))
        
        
        coreSpecies = [self.CH4,self.CH3,self.C2H6,self.C2H5, self.H2]
        edgeSpecies = []
        coreReactions = rxnList
        edgeReactions = []
        numCoreSpecies = len(coreSpecies)
        
        c0={self.CH4:0.2,self.CH3:0.1,self.C2H6:0.35,self.C2H5:0.15, self.H2:0.2}

        rxnSystem0 = LiquidReactor(self.T, c0,termination=[])
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
        
        c0={self.CH4:0.2,self.CH3:0.1,self.C2H6:0.35,self.C2H5:0.15, self.H2:0.2}

        for i in xrange(len(rxnList)):
            k0 = rxnList[i].getRateCoefficient(self.T)
            rxnList[i].kinetics.A.value_si = rxnList[i].kinetics.A.value_si*(1+1e-3)               
            dk = rxnList[i].getRateCoefficient(self.T) - k0

            rxnSystem = LiquidReactor(self.T, c0,termination=[])
            rxnSystem.initializeModel(coreSpecies, coreReactions, edgeSpecies, edgeReactions)

            dfdt = rxnSystem.residual(0.0, rxnSystem.y, numpy.zeros(rxnSystem.y.shape))[0]  
            dfdk[:,i]=(dfdt-dfdt0)/dk          
            
            
            rxnSystem.termination.append(TerminationTime((integrationTime,'s')))
            rxnSystem.simulate(coreSpecies, coreReactions, [], [], 0, 1, 0)
            
            rxnList[i].kinetics.A.value_si = rxnList[i].kinetics.A.value_si/(1+1e-3)  # reset A factor
            
        for i in xrange(numCoreSpecies):
            for j in xrange(len(rxnList)):
                self.assertAlmostEqual(dfdk[i,j], solver_dfdk[i,j], delta=abs(1e-3*dfdk[i,j]))