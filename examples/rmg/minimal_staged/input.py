# Data sources
database(
    thermoLibraries = ['primaryThermoLibrary'],
    reactionLibraries = [],
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
        'ethane': 0.9,
    },
    terminationTime=(1e6,'s'),
)

simulator( #first stage
    atol=1e-16,
    rtol=1e-8,
)

simulator( #second stage
    atol=1e-16,
    rtol=1e-8,
)

model(#first stage
    toleranceKeepInEdge=0.0,
    toleranceMoveToCore=1.0,
    toleranceInterruptSimulation=1.0,
    maximumEdgeSpecies=100000,
    maxNumObjsPerIter=1,#number of objects (species, reactions, PDepNetworks) that can be taken from one simulation
    maxNumSpecies = 10, #after this number of species is hit will move to next stage
)

model( #second stage
    toleranceKeepInEdge=0.0,
    toleranceMoveToCore=0.1,
    toleranceInterruptSimulation=1.0,
    maximumEdgeSpecies=100000,
    maxNumObjsPerIter=3,
    maxNumSpecies = 100,
)


options(
    units='si',
    saveRestartPeriod=None,
    generateOutputHTML=True,
    generatePlots=False,
    saveEdgeSpecies=True,
    saveSimulationProfiles=True,
)
