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
constantVIdealGasReactor(
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
    toleranceMoveToCore=0.01,
    toleranceInterruptSimulation=0.01,
)

options(
    units='si',
    saveEdgeSpecies=True,
    generateRMSEachIter=False,
)