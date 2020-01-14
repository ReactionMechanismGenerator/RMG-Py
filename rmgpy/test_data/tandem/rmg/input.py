# Data sources
database(
    thermoLibraries = ['primaryThermoLibrary'],
    reactionLibraries = [],
    seedMechanisms = [],
    kineticsDepositories = ['training'],
    kineticsFamilies = ['H_Abstraction','Disproportionation','R_Recombination'],
    kineticsEstimator = 'rate rules',
)

# List of species
species(
    label='H2',
    reactive=True,
    structure=SMILES("[H][H]"),
)
species(
    label='H',
    reactive=True,
    structure=SMILES("[H]"),
)
species(
   label='O2',
   reactive=True,
   structure=SMILES("[O][O]"),
)

# Reaction systems
simpleReactor(
    temperature=(1000, 'K'),
    pressure=(1.0, 'bar'),
    initialMoleFractions={'H2':.67, 'O2':.33, 'H': 0.001},
    terminationConversion={'H2': 0.1},
    terminationTime=(1,'s'),
)

simulator(atol=1e-16, rtol=1e-8)

model(
    toleranceKeepInEdge=0.0,
    toleranceMoveToCore=0.01,
    toleranceInterruptSimulation=0.01,
    maximumEdgeSpecies=100000,
)

options(
    units='si',
    generateOutputHTML=False,
    generatePlots=False,
    saveEdgeSpecies=False,
    saveSimulationProfiles=False,
)
