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
    label='n-decylbz',
    reactive=True,
    structure=SMILES("CCCCCCCCCCc1ccccc1"),
)
species(
    label='n-C11',
    reactive=True,
    structure=SMILES("CCCCCCCCCCC"),
)
species(
    label='n-C13',
    reactive=True,
    structure=SMILES("CCCCCCCCCCCCC"),
)
species(
    label='n-C16',
    reactive=True,
    structure=SMILES("CCCCCCCCCCCCCCCC"),
)
species(
    label='n-C19',
    reactive=True,
    structure=SMILES("CCCCCCCCCCCCCCCCCCC"),
)
species(
    label='n-C21',
    reactive=True,
    structure=SMILES("CCCCCCCCCCCCCCCCCCCCC"),
)
species(
    label='1M-napthalene',
    reactive=True,
    structure=SMILES("Cc1cccc2ccccc12"),
)
species(
    label='O2',
    reactive=True,
    structure=SMILES("[O][O]"),
)

# Reaction systems
simpleReactor(
    temperature=(500,'K'),
    pressure=(201.0,'bar'),
    initialMoleFractions={
        "n-C11": 0.05,
        "n-C13": 0.19,
        "n-C16": 0.25,
        "n-C19": 0.18,
        "n-C21": 0.10,
        "n-decylbz": 0.12,
        "1M-napthalene": 0.1,
        "O2": 0.05,
    },
    terminationConversion={
        'O2': 0.5,
    },
    terminationTime=(3600,'s'),
)

simulator(
    atol=1e-16,
    rtol=1e-8,
)

model(
    toleranceKeepInEdge=0.0001,
    toleranceMoveToCore=0.01,
    toleranceInterruptSimulation=0.5,
    maximumEdgeSpecies=1000
)

options(
    units='si',
    saveRestartPeriod=None,
    drawMolecules=False,
    generatePlots=False,
)
