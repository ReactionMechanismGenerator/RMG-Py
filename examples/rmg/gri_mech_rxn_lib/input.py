# Data sources
database(
    thermoLibraries = ['primaryThermoLibrary'],
    reactionLibraries = [('GRI-Mech3.0',False)],#False = unused reactions will not be  added to the chemkin file
    seedMechanisms = [],
    kineticsDepositories = ['training'],
    kineticsFamilies = 'default',
    kineticsEstimator = 'rate rules',
)

# List of species
species(
    label='ethane',
    reactive=True,
    structure=SMILES("CC"),
)

# Reaction systems
simpleReactor(
    temperature=(1350,'K'),
    pressure=(1.0,'bar'),
    initialMoleFractions={
        "ethane": 1.0,
    },
    terminationConversion={
        'ethane': 0.1,
    },
    terminationTime=(1e-1,'s'),
)

simulator(
    atol=1e-16,
    rtol=1e-8,
)

model(
    toleranceKeepInEdge=0.0,
    toleranceMoveToCore=0.1,
    toleranceInterruptSimulation=0.1,
    maximumEdgeSpecies=100000
)

options(
    units='si',
    saveRestartPeriod=None,
    generateOutputHTML=False,
    generatePlots=False,
    saveEdgeSpecies=True,
    saveSimulationProfiles=True,
)
