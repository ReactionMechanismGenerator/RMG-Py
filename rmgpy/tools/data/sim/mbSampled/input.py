# Data sources
database(
    thermoLibraries=['BurkeH2O2', 'primaryThermoLibrary', 'DFT_QCI_thermo', 'CBS_QB3_1dHR', 'FFCM1(-)'],
    reactionLibraries=[('BurkeH2O2inArHe', False), ('FFCM1(-)', False)],
    seedMechanisms=[],
    kineticsDepositories=['training'],
    kineticsFamilies='default',
    kineticsEstimator='rate rules',
)

generatedSpeciesConstraints(
    allowed=['seed mechanisms', 'input species', 'reaction libraries'],
    maximumOxygenAtoms=4,
    maximumCarbonAtoms=18,
    maximumRadicalElectrons=1,
    allowSingletO2=False,
)

# List of species
species(
    label='C6H5',
    reactive=True,
    structure=adjacencyList(
        """
        multiplicity 2
        1 C u0 p0 c0 {2,B} {3,B} {8,S}
        2 C u0 p0 c0 {1,B} {4,B} {7,S}
        3 C u0 p0 c0 {1,B} {5,B} {9,S}
        4 C u0 p0 c0 {2,B} {6,B} {10,S}
        5 C u0 p0 c0 {3,B} {6,B} {11,S}
        6 C u1 p0 c0 {4,B} {5,B}
        7 H u0 p0 c0 {2,S}
        8 H u0 p0 c0 {1,S}
        9 H u0 p0 c0 {3,S}
        10 H u0 p0 c0 {4,S}
        11 H u0 p0 c0 {5,S}
        """),
)

species(
    label='C2H4',
    reactive=True,
    structure=SMILES('C=C'),
)

species(
    label='He',
    reactive=False,
    structure=SMILES("[He]"),
)

mbsampledReactor(
    temperature=(796, 'K'),
    pressure=(10, 'torr'),
    initialMoleFractions={
        "C6H5": 0.7e+12,
        "C2H4": 2.68e+16,
        "He": 9.38e+16,
    },
    mbsamplingRate=3500,
    terminationTime=(0.01, 's'),
    constantSpecies=['C2H4', 'He']
)

simulator(
    atol=1e-10,
    rtol=1e-6,
    sens_atol=1e-6,
    sens_rtol=1e-4,
)

model(
    toleranceKeepInEdge=0.01,
    toleranceMoveToCore=0.1,
    toleranceInterruptSimulation=1E8,
    maximumEdgeSpecies=100000,
    filterReactions=True,
)

options(
    units='si',
    saveSimulationProfiles=True,
    generateOutputHTML=False,
    generatePlots=False,
    saveEdgeSpecies=False,
)
