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

species(
    label='O2',
    reactive=True,
    structure=SMILES('[O][O]')
)

species(
   label='N2',
   reactive=False,
   structure=SMILES('N#N'),
)

# Reaction systems
simpleReactor(
    temperature=[(1000,'K'),(1500,'K')],
    pressure=[(1.0,'bar'),(10.0,'bar')],
    nSimsTerm=12, 
    initialMoleFractions={
        "ethane": [0.05,0.15],
        "O2": 0.1,
        "N2": 0.9,
    },
    terminationConversion={
        'ethane': 0.1,
    },
    terminationTime=(1e1,'s'),
    balanceSpecies = "N2",
)

simulator(
    atol=1e-16,
    rtol=1e-8,
)

model(
    toleranceKeepInEdge=0.0,
    toleranceMoveToCore=0.1,
    toleranceInterruptSimulation=0.1,
    maximumEdgeSpecies=100000,
    filterReactions=True,
)

options(
    units='si',
    saveRestartPeriod=None,
    generateOutputHTML=False,
    generatePlots=False,
    saveEdgeSpecies=False,
    saveSimulationProfiles=False,
)
