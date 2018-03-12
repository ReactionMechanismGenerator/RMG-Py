#!/usr/bin/python
# -*- coding: utf-8 -*-

import unittest
import numpy

import rmgpy.quantity

from rmgpy.molecule import Molecule
from rmgpy.species import Species
from rmgpy.reaction import Reaction
from rmgpy.kinetics import SurfaceArrhenius, StickingCoefficient
from rmgpy.thermo import ThermoData, NASA, NASAPolynomial
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
        We use a SurfaceArrhenius for the rate expression.
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
        surfaceSpecies = []
        surfaceReactions = []

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
                                  edgeReactions, surfaceSpecies, surfaceReactions)

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
        pylab.savefig('surfaceTestH2.pdf')

        return


    def testSolveCH3(self):
        """
        Test the surface batch reactor with a nondissociative adsorption of CH3
        
        Here we choose a kinetic model consisting of the  adsorption reaction
        CH3 + X <=>  CH3X
        We use a sticking coefficient for the rate expression.
        """

        CH3 = Species(
            molecule=[Molecule().fromSMILES("[CH3]")],
            thermo=NASA(polynomials=[NASAPolynomial(coeffs=[3.91547, 0.00184155, 3.48741e-06, -3.32746e-09, 8.49953e-13, 16285.6, 0.351743], Tmin=(100, 'K'), Tmax=(1337.63, 'K')),
                                     NASAPolynomial(coeffs=[3.54146, 0.00476786, -1.82148e-06, 3.28876e-10, -2.22545e-14, 16224, 1.66032], Tmin=(1337.63, 'K'), Tmax=(5000, 'K'))],
                       Tmin=(100, 'K'), Tmax=(5000, 'K'), E0=(135.382, 'kJ/mol'),
                       comment="""Thermo library: primaryThermoLibrary + radical(CH3)"""
                       ),
            molecularWeight=(15.0345, 'amu'),
                    )

        X = Species(
            molecule=[Molecule().fromAdjacencyList("1 X u0 p0")],
            thermo=NASA(polynomials=[NASAPolynomial(coeffs=[0, 0, 0, 0, 0, 0, 0], Tmin=(298, 'K'), Tmax=(1000, 'K')),
                                    NASAPolynomial(coeffs=[0, 0, 0, 0, 0, 0, 0], Tmin=(1000, 'K'), Tmax=(2000, 'K'))],
                        Tmin=(298, 'K'), Tmax=(2000, 'K'), E0=(-6.19426, 'kJ/mol'),
                        comment="""Thermo library: surfaceThermo""")
                    )

        CH3X = Species(
            molecule=[Molecule().fromAdjacencyList("1 H u0 p0 {2,S} \n 2 X u0 p0 {1,S}")],
            thermo=NASA(polynomials=[NASAPolynomial(coeffs=[-0.552219, 0.026442, -3.55617e-05, 2.60044e-08, -7.52707e-12, -4433.47, 0.692144], Tmin=(298, 'K'), Tmax=(1000, 'K')),
                                     NASAPolynomial(coeffs=[3.62557, 0.00739512, -2.43797e-06, 1.86159e-10, 3.6485e-14, -5187.22, -18.9668], Tmin=(1000, 'K'), Tmax=(2000, 'K'))],
                        Tmin=(298, 'K'), Tmax=(2000, 'K'), E0=(-39.1285, 'kJ/mol'),
                        comment="""Thermo library: surfaceThermo""")
                    )

        rxn1 = Reaction(reactants=[CH3, X],
                        products=[CH3X],
                        kinetics=StickingCoefficient(A=0.1, n=0, Ea=(0, 'kcal/mol'),
                                                     T0=(1, 'K'),
                                                     Tmin=(200, 'K'), Tmax=(3000, 'K'),
                                                     comment="""Exact match found for rate rule (Adsorbate;VacantSite)"""
                                                     )
#                        kinetics=SurfaceArrhenius(A=(2.7e10, 'cm^3/(mol*s)'),
#                                           n=0.5,
#                                           Ea=(5.0, 'kJ/mol'),
#                                           T0=(1.0, 'K'))
                        )
        coreSpecies = [CH3, X, CH3X]
        edgeSpecies = []
        coreReactions = [rxn1]
        edgeReactions = []
        surfaceSpecies = []
        surfaceReactions = []

        T = 800.
        initialP = 1.0e5
        rxnSystem = SurfaceReactor(
            T, initialP,
            initialGasMoleFractions={CH3: 1.0},
            initialSurfaceCoverages={X: 1.0},
            surfaceVolumeRatio=(1., 'm^-1'),
            surfaceSiteDensity=(2.72e-9, 'mol/cm^2'),
            termination=[])
        # in chemkin, the sites are mostly occupied in about 1e-8 seconds.

        rxnSystem.initializeModel(coreSpecies, coreReactions, edgeSpecies,
                                  edgeReactions, surfaceSpecies, surfaceReactions)

        tlist = numpy.logspace(-13, -5, 81, dtype=numpy.float64)

        print "Surface site density:", rxnSystem.surfaceSiteDensity.value_si

        print "rxn1 rate coefficient", rxn1.getSurfaceRateCoefficient(rxnSystem.T.value_si,
                                            rxnSystem.surfaceSiteDensity.value_si
                                            )

        # Integrate to get the solution at each time point
        t = []
        y = []
        reactionRates = []
        speciesRates = []
        t.append(rxnSystem.t)
        # You must make a copy of y because it is overwritten by DASSL at
        # each call to advance()
        y.append(rxnSystem.y.copy())
        reactionRates.append(rxnSystem.coreReactionRates.copy())
        speciesRates.append(rxnSystem.coreSpeciesRates.copy())
        print "time: ", t
        print "moles:", y
        print "reaction rates:", reactionRates
        print "species rates:", speciesRates
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
            self.assertAlmostEqual(reactionRates[i, 0], -speciesRates[i, 0],
                                   delta=1e-6 * reactionRates[i, 0])
            self.assertAlmostEqual(reactionRates[i, 0], -speciesRates[i, 1],
                                   delta=1e-6 * reactionRates[i, 0])
            self.assertAlmostEqual(reactionRates[i, 0], speciesRates[i, 2],
                                   delta=1e-6 * reactionRates[i, 0])

        # Check that we've reached equilibrium by the end
        self.assertAlmostEqual(reactionRates[-1, 0], 0.0, delta=1e-2)

        # Visualize the simulation results
        import pylab
        fig = pylab.figure(figsize=(6, 6))
        pylab.subplot(2, 1, 1)
        pylab.semilogx(t, y[:, 2])
        pylab.ylabel('Concentration (mol/m$^\\mathdefault{3 or 2}$)')
        pylab.legend(['CH3X'], loc=4)
        pylab.subplot(2, 1, 2)
        pylab.semilogx(t, speciesRates)
        pylab.legend(['CH3', 'X', 'CH3X'], loc=4)
        pylab.xlabel('Time (s)')
        pylab.ylabel('Rate (mol/m$^\\mathdefault{3 or 2}$*s)')
        #fig.subplots_adjust(left=0.21, bottom=0.10, right=0.95, top=0.95, wspace=0.20, hspace=0.35)
        pylab.tight_layout()
        #pylab.show()
        pylab.savefig('surfaceTestCH3.pdf')

        return

        # Dump of core reactions from a test simulation:
        '''
        [ Reaction(
        index = 1,
        reactants = [
            Species(
                index = 1,
                label = 'methyl',
                thermo = NASA(
                    polynomials = [
                        NASAPolynomial(
                            coeffs = [3.91547, 0.00184155, 3.48741e-06, -3.32746e-09, 8.49953e-13, 16285.6, 0.351743],
                            Tmin = (100, 'K'),
                            Tmax = (1337.63, 'K'),
                        ),
                        NASAPolynomial(
                            coeffs = [3.54146, 0.00476786, -1.82148e-06, 3.28876e-10, -2.22545e-14, 16224, 1.66032],
                            Tmin = (1337.63, 'K'),
                            Tmax = (5000, 'K'),
                        ),
                    ],
                    Tmin = (100, 'K'),
                    Tmax = (5000, 'K'),
                    E0 = (135.382, 'kJ/mol'),
                    comment = 'Thermo library: primaryThermoLibrary + radical(CH3)',
                ),
                conformer = Conformer(E0=(135.382, 'kJ/mol'), modes=[]),
                molecule = [
                    Molecule(SMILES='[CH3]'),
                ],
                transportData = TransportData(
                    shapeIndex = 2,
                    epsilon = (1197.29, 'J/mol'),
                    sigma = (3.8, 'angstroms'),
                    dipoleMoment = (0, 'C*m'),
                    polarizability = (0, 'angstroms^3'),
                    rotrelaxcollnum = 0,
                    comment = 'GRI-Mech',
                ),
                molecularWeight = (15.0345, 'amu'),
                energyTransferModel = SingleExponentialDown(alpha0=(3.5886, 'kJ/mol'), T0=(300, 'K'), n=0.85),
            ),
            Species(
                index = 2,
                label = 'site',
                thermo = NASA(
                    polynomials = [
                        NASAPolynomial(coeffs=[0, 0, 0, 0, 0, 0, 0], Tmin=(298, 'K'), Tmax=(1000, 'K')),
                        NASAPolynomial(coeffs=[0, 0, 0, 0, 0, 0, 0], Tmin=(1000, 'K'), Tmax=(2000, 'K')),
                    ],
                    Tmin = (298, 'K'),
                    Tmax = (2000, 'K'),
                    E0 = (-6.19426, 'kJ/mol'),
                    comment = 'Thermo library: surfaceThermo',
                ),
                conformer = Conformer(E0=(-6.19426, 'kJ/mol'), modes=[]),
                molecule = [
                    Molecule(SMILES='[Ni]'),
                ],
                transportData = TransportData(
                    shapeIndex = 2,
                    epsilon = (1235.53, 'J/mol'),
                    sigma = (3.758e-10, 'm'),
                    dipoleMoment = (0, 'C*m'),
                    polarizability = (0, 'angstroms^3'),
                    rotrelaxcollnum = 0,
                    comment = 'Epsilon & sigma estimated with fixed Lennard Jones Parameters. This is the fallback method! Try improving transport databases!',
                ),
                molecularWeight = (0, 'amu'),
                energyTransferModel = SingleExponentialDown(alpha0=(3.5886, 'kJ/mol'), T0=(300, 'K'), n=0.85),
            ),
        ],
        products = [
            Species(
                index = 3,
                label = 'C[Ni]',
                thermo = NASA(
                    polynomials = [
                        NASAPolynomial(
                            coeffs = [-0.552219, 0.026442, -3.55617e-05, 2.60044e-08, -7.52707e-12, -4433.47, 0.692144],
                            Tmin = (298, 'K'),
                            Tmax = (1000, 'K'),
                        ),
                        NASAPolynomial(
                            coeffs = [3.62557, 0.00739512, -2.43797e-06, 1.86159e-10, 3.6485e-14, -5187.22, -18.9668],
                            Tmin = (1000, 'K'),
                            Tmax = (2000, 'K'),
                        ),
                    ],
                    Tmin = (298, 'K'),
                    Tmax = (2000, 'K'),
                    E0 = (-39.1285, 'kJ/mol'),
                    comment = 'Thermo library: surfaceThermo',
                ),
                conformer = Conformer(E0=(-39.1285, 'kJ/mol'), modes=[]),
                molecule = [
                    Molecule(SMILES='C[Ni]'),
                ],
                transportData = TransportData(
                    shapeIndex = 2,
                    epsilon = (920.412, 'J/mol'),
                    sigma = (4.443e-10, 'm'),
                    dipoleMoment = (0, 'C*m'),
                    polarizability = (0, 'angstroms^3'),
                    rotrelaxcollnum = 0,
                    comment = 'Epsilon & sigma estimated with fixed Lennard Jones Parameters. This is the fallback method! Try improving transport databases!',
                ),
                molecularWeight = (15.0345, 'amu'),
                energyTransferModel = SingleExponentialDown(alpha0=(3.5886, 'kJ/mol'), T0=(300, 'K'), n=0.85),
            ),
        ],
        kinetics = StickingCoefficient(
            A = 0.1,
            n = 0,
            Ea = (0, 'kcal/mol'),
            T0 = (1, 'K'),
            Tmin = (200, 'K'),
            Tmax = (3000, 'K'),
            comment = 'Exact match found for rate rule (Adsorbate;VacantSite)',
        ),
        pairs = [
            (
                Species(
                    index = 1,
                    label = 'methyl',
                    thermo = NASA(
                        polynomials = [
                            NASAPolynomial(
                                coeffs = [3.91547, 0.00184155, 3.48741e-06, -3.32746e-09, 8.49953e-13, 16285.6, 0.351743],
                                Tmin = (100, 'K'),
                                Tmax = (1337.63, 'K'),
                            ),
                            NASAPolynomial(
                                coeffs = [3.54146, 0.00476786, -1.82148e-06, 3.28876e-10, -2.22545e-14, 16224, 1.66032],
                                Tmin = (1337.63, 'K'),
                                Tmax = (5000, 'K'),
                            ),
                        ],
                        Tmin = (100, 'K'),
                        Tmax = (5000, 'K'),
                        E0 = (135.382, 'kJ/mol'),
                        comment = 'Thermo library: primaryThermoLibrary + radical(CH3)',
                    ),
                    conformer = Conformer(E0=(135.382, 'kJ/mol'), modes=[]),
                    molecule = [
                        Molecule(SMILES='[CH3]'),
                    ],
                    transportData = TransportData(
                        shapeIndex = 2,
                        epsilon = (1197.29, 'J/mol'),
                        sigma = (3.8, 'angstroms'),
                        dipoleMoment = (0, 'C*m'),
                        polarizability = (0, 'angstroms^3'),
                        rotrelaxcollnum = 0,
                        comment = 'GRI-Mech',
                    ),
                    molecularWeight = (15.0345, 'amu'),
                    energyTransferModel = SingleExponentialDown(alpha0=(3.5886, 'kJ/mol'), T0=(300, 'K'), n=0.85),
                ),
                Species(
                    index = 3,
                    label = 'C[Ni]',
                    thermo = NASA(
                        polynomials = [
                            NASAPolynomial(
                                coeffs = [-0.552219, 0.026442, -3.55617e-05, 2.60044e-08, -7.52707e-12, -4433.47, 0.692144],
                                Tmin = (298, 'K'),
                                Tmax = (1000, 'K'),
                            ),
                            NASAPolynomial(
                                coeffs = [3.62557, 0.00739512, -2.43797e-06, 1.86159e-10, 3.6485e-14, -5187.22, -18.9668],
                                Tmin = (1000, 'K'),
                                Tmax = (2000, 'K'),
                            ),
                        ],
                        Tmin = (298, 'K'),
                        Tmax = (2000, 'K'),
                        E0 = (-39.1285, 'kJ/mol'),
                        comment = 'Thermo library: surfaceThermo',
                    ),
                    conformer = Conformer(E0=(-39.1285, 'kJ/mol'), modes=[]),
                    molecule = [
                        Molecule(SMILES='C[Ni]'),
                    ],
                    transportData = TransportData(
                        shapeIndex = 2,
                        epsilon = (920.412, 'J/mol'),
                        sigma = (4.443e-10, 'm'),
                        dipoleMoment = (0, 'C*m'),
                        polarizability = (0, 'angstroms^3'),
                        rotrelaxcollnum = 0,
                        comment = 'Epsilon & sigma estimated with fixed Lennard Jones Parameters. This is the fallback method! Try improving transport databases!',
                    ),
                    molecularWeight = (15.0345, 'amu'),
                    energyTransferModel = SingleExponentialDown(alpha0=(3.5886, 'kJ/mol'), T0=(300, 'K'), n=0.85),
                ),
            ),
            (
                Species(
                    index = 2,
                    label = 'site',
                    thermo = NASA(
                        polynomials = [
                            NASAPolynomial(coeffs=[0, 0, 0, 0, 0, 0, 0], Tmin=(298, 'K'), Tmax=(1000, 'K')),
                            NASAPolynomial(coeffs=[0, 0, 0, 0, 0, 0, 0], Tmin=(1000, 'K'), Tmax=(2000, 'K')),
                        ],
                        Tmin = (298, 'K'),
                        Tmax = (2000, 'K'),
                        E0 = (-6.19426, 'kJ/mol'),
                        comment = 'Thermo library: surfaceThermo',
                    ),
                    conformer = Conformer(E0=(-6.19426, 'kJ/mol'), modes=[]),
                    molecule = [
                        Molecule(SMILES='[Ni]'),
                    ],
                    transportData = TransportData(
                        shapeIndex = 2,
                        epsilon = (1235.53, 'J/mol'),
                        sigma = (3.758e-10, 'm'),
                        dipoleMoment = (0, 'C*m'),
                        polarizability = (0, 'angstroms^3'),
                        rotrelaxcollnum = 0,
                        comment = 'Epsilon & sigma estimated with fixed Lennard Jones Parameters. This is the fallback method! Try improving transport databases!',
                    ),
                    molecularWeight = (0, 'amu'),
                    energyTransferModel = SingleExponentialDown(alpha0=(3.5886, 'kJ/mol'), T0=(300, 'K'), n=0.85),
                ),
                Species(
                    index = 3,
                    label = 'C[Ni]',
                    thermo = NASA(
                        polynomials = [
                            NASAPolynomial(
                                coeffs = [-0.552219, 0.026442, -3.55617e-05, 2.60044e-08, -7.52707e-12, -4433.47, 0.692144],
                                Tmin = (298, 'K'),
                                Tmax = (1000, 'K'),
                            ),
                            NASAPolynomial(
                                coeffs = [3.62557, 0.00739512, -2.43797e-06, 1.86159e-10, 3.6485e-14, -5187.22, -18.9668],
                                Tmin = (1000, 'K'),
                                Tmax = (2000, 'K'),
                            ),
                        ],
                        Tmin = (298, 'K'),
                        Tmax = (2000, 'K'),
                        E0 = (-39.1285, 'kJ/mol'),
                        comment = 'Thermo library: surfaceThermo',
                    ),
                    conformer = Conformer(E0=(-39.1285, 'kJ/mol'), modes=[]),
                    molecule = [
                        Molecule(SMILES='C[Ni]'),
                    ],
                    transportData = TransportData(
                        shapeIndex = 2,
                        epsilon = (920.412, 'J/mol'),
                        sigma = (4.443e-10, 'm'),
                        dipoleMoment = (0, 'C*m'),
                        polarizability = (0, 'angstroms^3'),
                        rotrelaxcollnum = 0,
                        comment = 'Epsilon & sigma estimated with fixed Lennard Jones Parameters. This is the fallback method! Try improving transport databases!',
                    ),
                    molecularWeight = (15.0345, 'amu'),
                    energyTransferModel = SingleExponentialDown(alpha0=(3.5886, 'kJ/mol'), T0=(300, 'K'), n=0.85),
                ),
            ),
        ],
    )]
            '''


