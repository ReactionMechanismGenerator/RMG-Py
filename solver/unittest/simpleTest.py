#!/usr/bin/python
# -*- coding: utf-8 -*-

import unittest
import numpy

import quantities
quantities.set_default_units('si')
quantities.UnitQuantity('kilocalorie', 1000.0*quantities.cal, symbol='kcal')
quantities.UnitQuantity('kilojoule', 1000.0*quantities.J, symbol='kJ')

from chempy.molecule import Molecule
from chempy.species import Species
from chempy.reaction import Reaction
from chempy.kinetics import ArrheniusModel
from chempy.thermo import ThermoGAModel
from rmgsolver.simple import SimpleReactor

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
            thermo=ThermoGAModel(Tdata=([300,400,500,600,800,1000,1500],"K"), Cpdata=([ 8.615, 9.687,10.963,12.301,14.841,16.976,20.528],"cal/mol*K"), H298=(-17.714,"kcal/mol"), S298=(44.472,"cal/mol*K"))
            )
        CH3 = Species(
            molecule=[Molecule().fromSMILES("[CH3]")],
            thermo=ThermoGAModel(Tdata=([300,400,500,600,800,1000,1500],"K"), Cpdata=([ 9.397,10.123,10.856,11.571,12.899,14.055,16.195],"cal/mol*K"), H298=(  9.357,"kcal/mol"), S298=(45.174,"cal/mol*K"))
            )
        C2H6 = Species(
            molecule=[Molecule().fromSMILES("CC")],
            thermo=ThermoGAModel(Tdata=([300,400,500,600,800,1000,1500],"K"), Cpdata=([12.684,15.506,18.326,20.971,25.500,29.016,34.595],"cal/mol*K"), H298=(-19.521,"kcal/mol"), S298=(54.799,"cal/mol*K"))
            )
        C2H5 = Species(
            molecule=[Molecule().fromSMILES("C[CH2]")],
            thermo=ThermoGAModel(Tdata=([300,400,500,600,800,1000,1500],"K"), Cpdata=([11.635,13.744,16.085,18.246,21.885,24.676,29.107],"cal/mol*K"), H298=( 29.496,"kcal/mol"), S298=(56.687,"cal/mol*K"))
            )

        rxn1 = Reaction(reactants=[C2H6,CH3], products=[C2H5,CH4], kinetics=ArrheniusModel(A=686.375*6, n=4.40721, Ea=7.82799*4184., T0=298.15))

        coreSpecies = [CH4,CH3,C2H6,C2H5]
        edgeSpecies = []
        coreReactions = [rxn1]
        edgeReactions = []

        T = 1000; P = 1.0e5; C = P / 8.314472 / T
        rxnSystem = SimpleReactor(T, P, C0={C2H5: 0.1*C, CH3: 0.1*C, CH4: 0.4*C, C2H6: 0.4*C})
        rxnSystem.initialize(coreSpecies, coreReactions, edgeSpecies, edgeReactions)

        tlist = numpy.array([10**(i/10.0) for i in range(-130, -49)], numpy.float64)

        # Integrate to get the solution at each time point
        t = []; y = []; reactionRates = []; speciesRates = []
        for t1 in tlist:
            rxnSystem.advance(t1)
            t.append(rxnSystem.t)
            # You must make a copy of y because it is overwritten by DASSL at
            # each call to advance()
            y.append(rxnSystem.y.copy())
            reactionRates.append(rxnSystem.reactionRates.copy())
            speciesRates.append(rxnSystem.speciesRates.copy())

        # Convert the solution vectors to numpy arrays
        t = numpy.array(t, numpy.float64)
        y = numpy.array(y, numpy.float64)
        reactionRates = numpy.array(reactionRates, numpy.float64)
        speciesRates = numpy.array(speciesRates, numpy.float64)

        import pylab
        fig = pylab.figure(figsize=(6,6))
        pylab.subplot(2,1,1)
        pylab.semilogx(t, y)
        pylab.ylabel('Concentration (mol/m$^\\mathdefault{3}$)')
        pylab.legend(['CH4', 'CH3', 'C2H6', 'C2H5'], loc=4)
        pylab.subplot(2,1,2)
        pylab.semilogx(t, speciesRates)
        pylab.legend(['CH4', 'CH3', 'C2H6', 'C2H5'], loc=4)
        pylab.xlabel('Time (s)')
        pylab.ylabel('Rate (mol/m$^\\mathdefault{3}$*s)')
        fig.subplots_adjust(left=0.12, bottom=0.10, right=0.95, top=0.95, wspace=0.20, hspace=0.35)
        pylab.show()

################################################################################

if __name__ == '__main__':
    unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )