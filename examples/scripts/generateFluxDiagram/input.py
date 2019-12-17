# Data sources
database(
    thermoLibraries = ['primaryThermoLibrary'],
    reactionLibraries = [],
    seedMechanisms = [],
    kineticsDepositories = ['training'],
    kineticsFamilies = ['H_Abstraction','Disproportionation','R_Recombination',
                        'Birad_recombination', 'Birad_R_Recombination'],
    kineticsEstimator = 'rate rules',
)

# List of species
species(
    label='H2',
    reactive=True,
    structure=SMILES("[H][H]"),
)
species(
   label='O2',
   reactive=True,
   structure=SMILES("[O][O]"),
)

# Reaction systems
simpleReactor(
    temperature=(1000,'K'),
    pressure=(1.0,'bar'),
    initialMoleFractions={
        'H2':.67, 'O2':.33,
    },
    terminationConversion={
        'H2': 0.9,
    },
    terminationTime=(1e6,'s'),
)

simulator(
    atol=1e-16,
    rtol=1e-8,
)

model(
    toleranceKeepInEdge=0.0,
    toleranceMoveToCore=0.001,
    toleranceInterruptSimulation=0.001,
    maximumEdgeSpecies=100000,
)

options(
    units='si',
    generateOutputHTML=True,
    generatePlots=False,
    saveEdgeSpecies=True,
    saveSimulationProfiles=True,
)
