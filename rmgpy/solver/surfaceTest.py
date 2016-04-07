#!/usr/bin/python
# -*- coding: utf-8 -*-

import unittest
import numpy

import rmgpy.quantity

from rmgpy.molecule import Molecule
from rmgpy.species import Species
from rmgpy.reaction import Reaction
from rmgpy.kinetics import SurfaceArrhenius
from rmgpy.thermo import ThermoData
from rmgpy.solver.surface import SurfaceReactor
from rmgpy.solver.base import TerminationTime, TerminationConversion
import rmgpy.constants as constants

################################################################################


class SurfaceReactorCheck(unittest.TestCase):
    def testSolveH2(self):
        """
        Test the surface batch reactor with a dissociative adsorption of H2
        
        Here we choose a kinetic model consisting of the dissociative adsorption reaction
        H2 + 2X <=> 2 HX
        """
        H2 = Species(
            molecule=[Molecule().fromSMILES("[H][H]")],
            thermo=ThermoData(Tdata=([300, 400, 500, 600, 800, 1000, 1500],
                                     "K"),
                              Cpdata=([6.955, 6.955, 6.956, 6.961, 7.003,
                                       7.103, 7.502], "cal/(mol*K)"),
                              H298=(0, "kcal/mol"),
                              S298=(31.129  , "cal/(mol*K)")))
        X = Species(
            molecule=[Molecule().fromAdjacencyList("1 X u0 p0")],
            thermo=ThermoData(Tdata=([300, 400, 500, 600, 800, 1000, 1500],
                                     "K"),
                              Cpdata=([0., 0., 0., 0., 0., 0., 0.], "cal/(mol*K)"),
                              H298=(0.0, "kcal/mol"),
                              S298=(0.0, "cal/(mol*K)")))
        HX = Species(
            molecule=[Molecule().fromAdjacencyList("1 H u0 p0 {2,S} \n 2 X u0 p0 {1,S}")],
            thermo=ThermoData(Tdata=([300, 400, 500, 600, 800, 1000, 1500],
                                     "K"),
                              Cpdata=([1.50, 2.58, 3.40, 4.00, 4.73, 5.13, 5.57], "cal/(mol*K)"),
                              H298=(-11.26, "kcal/mol"),
                              S298=(0.44, "cal/(mol*K)")))

        rxn1 = Reaction(reactants=[H2, X, X],
                        products=[HX, HX],
                        kinetics=SurfaceArrhenius(A=(9.05e18, 'cm^5/(mol^2*s)'),
                                           n=0.5,
                                           Ea=(5.0, 'kJ/mol'),
                                           T0=(1.0, 'K')))

        coreSpecies = [H2, X, HX]
        edgeSpecies = []
        coreReactions = [rxn1]
        edgeReactions = []

        T = 1000
        initialP = 1.0e5
        rxnSystem = SurfaceReactor(
            T, initialP,
            initialGasMoleFractions={H2: 1.0},
            initialSurfaceCoverages={X: 1.0},
            surfaceVolumeRatio=(1e1, 'm^-1'),
            surfaceSiteDensity=(2.72e-9, 'mol/cm^2'),
            termination=[])

        rxnSystem.initializeModel(coreSpecies, coreReactions, edgeSpecies,
                                  edgeReactions)

        tlist = numpy.logspace(-13, -5, 81, dtype=numpy.float64)

        # Integrate to get the solution at each time point
        t = []
        y = []
        reactionRates = []
        speciesRates = []
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
        V = constants.R * rxnSystem.T.value_si * numpy.sum(y) / rxnSystem.initialP.value_si

        # Check that we're computing the species fluxes correctly
        for i in range(t.shape[0]):
            self.assertAlmostEqual(reactionRates[i, 0], -1.0 * speciesRates[i, 0],
                                   delta=1e-6 * reactionRates[i, 0])
            self.assertAlmostEqual(reactionRates[i, 0], -0.5 * speciesRates[i, 1],
                                   delta=1e-6 * reactionRates[i, 0])
            self.assertAlmostEqual(reactionRates[i, 0], 0.5 * speciesRates[i, 2],
                                   delta=1e-6 * reactionRates[i, 0])

        # Check that we've reached equilibrium
        self.assertAlmostEqual(reactionRates[-1, 0], 0.0, delta=1e-2)

        # Visualize the simulation results
        import pylab
        fig = pylab.figure(figsize=(6, 6))
        pylab.subplot(2, 1, 1)
        pylab.semilogx(t, y[:, 2])
        pylab.ylabel('Concentration (mol/m$^\\mathdefault{3 or 2}$)')
        pylab.legend(['HX'], loc=4)
        pylab.subplot(2, 1, 2)
        pylab.semilogx(t, speciesRates)
        pylab.legend(['H2', 'X', 'HX'], loc=4)
        pylab.xlabel('Time (s)')
        pylab.ylabel('Rate (mol/m$^\\mathdefault{3 or 2}$*s)')
        #fig.subplots_adjust(left=0.21, bottom=0.10, right=0.95, top=0.95, wspace=0.20, hspace=0.35)
        pylab.tight_layout()
        #pylab.show()
        pylab.savefig('surfaceTest.pdf')

        return

