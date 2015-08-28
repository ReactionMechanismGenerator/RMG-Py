# Data sources
database(
    thermoLibraries = ['primaryThermoLibrary'],
    reactionLibraries = [],
    seedMechanisms = [],
    kineticsDepositories = ['training'], 
    kineticsFamilies = 'default',
    kineticsEstimator = 'rate rules',
)

# Constraints on generated species
generatedSpeciesConstraints(
    maximumCarbonAtoms = 7,
)

# List of species
species(
    label='n-heptane',
    reactive=True,
    structure=SMILES("CCCCCCC"),
)

species(
    label='Ar',
    reactive=False,
    structure=SMILES("[Ar]"),
)

simpleReactor(
    temperature=(1000,'K'),
    pressure=(400,'Pa'),
    initialMoleFractions={
        "n-heptane": 0.02,
        "Ar": 0.98,
    },
    terminationConversion={
        'n-heptane': 0.99,
    },
    terminationTime=(1e6,'s'),
)

simpleReactor(
    temperature=(1500,'K'),
    pressure=(400,'Pa'),
    initialMoleFractions={
        "n-heptane": 0.02,
        "Ar": 0.98,
    },
    terminationConversion={
        'n-heptane': 0.99,
    },
    terminationTime=(1e6,'s'),
)

simpleReactor(
    temperature=(2000,'K'),
    pressure=(400,'Pa'),
    initialMoleFractions={
        "n-heptane": 0.02,
        "Ar": 0.98,
    },
    terminationConversion={
        'n-heptane': 0.99,
    },
    terminationTime=(1e6,'s'),
)

simulator(
    atol=1e-16,
    rtol=1e-8,
)

model(
    toleranceKeepInEdge=0.0,
    toleranceMoveToCore=0.01,
    toleranceInterruptSimulation=0.1,
    maximumEdgeSpecies=100000
)

pressureDependence(
    method='modified strong collision',
    maximumGrainSize=(0.5,'kcal/mol'),
    minimumNumberOfGrains=250,
    temperatures=(300,3000,'K',8),
    pressures=(0.001,100,'bar',5),
    interpolation=('Chebyshev', 6, 4),
)

options(
    units='si',
    saveRestartPeriod=(1,'hour'),
    drawMolecules=False,
    generatePlots=False,
)
