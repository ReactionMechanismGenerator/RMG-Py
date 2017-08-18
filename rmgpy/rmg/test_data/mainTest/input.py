database(
    thermoLibraries = ['primaryThermoLibrary'],
    reactionLibraries = [],
    seedMechanisms = [],
    kineticsDepositories = ['training'],
    kineticsFamilies = ['R_Recombination'],
    kineticsEstimator = 'rate rules',
)

species(
    label='ethane',
    reactive=True,
    structure=SMILES("CC"),
)

simpleReactor(
    temperature=(1350,'K'),
    pressure=(1.0,'bar'),
    initialMoleFractions={
        "ethane": 1.0,
    },
    terminationConversion={
        'ethane': 0.000000000001,
    },
    terminationTime=(1e6,'s'),
)

simulator(
    atol=1e-16,
    rtol=1e-8
)

model(
    toleranceKeepInEdge=0.0,
    toleranceMoveToCore=0.2,
    toleranceInterruptSimulation=0.2,
)
options(
    name='testSeed',
    units='si',
    generateSeedEachIteration=True,
    saveSeedToDatabase=True,
    saveRestartPeriod=None,
    generateOutputHTML=False,
    generatePlots=False,
    saveEdgeSpecies=False,
    saveSimulationProfiles=False,
)
