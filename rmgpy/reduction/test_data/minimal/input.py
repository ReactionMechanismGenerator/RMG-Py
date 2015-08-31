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
    # terminationConversion={
    #     'ethane': 0.9,
    # },
    terminationTime=(1e-3,'s'),
)

# simpleReactor(
#     temperature=(1750,'K'),
#     pressure=(10.0,'bar'),
#     initialMoleFractions={
#         "ethane": 1.0,
#     },
#     # terminationConversion={
#     #     'ethane': 0.9,
#     # },
#     terminationTime=(1e-2,'s'),
# )

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
    generatePlots=False,
    saveEdgeSpecies=True,
    saveSimulationProfiles=True,
)
