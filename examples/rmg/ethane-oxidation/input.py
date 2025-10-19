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
    structure=SMILES("[O][O]"),
)


species(
    label='Ar',
    reactive=False,
    structure=SMILES("[Ar]"),
)

# Reaction systems
simpleReactor(
    temperature=(1500,'K'),
    pressure=(2.0,'bar'),
    initialMoleFractions={
        "ethane": 1,
        "O2": 3.5,
        "Ar": 1,
    },
    terminationConversion={
        'O2': 0.001,
    },
    terminationTime=(1,'s'),
)

simulator(
    atol=1e-16,
    rtol=1e-8,
)

model(
    toleranceKeepInEdge=0.0001,
    toleranceMoveToCore=0.5,
    maximumEdgeSpecies=1000
)

pressureDependence(
    method='modified strong collision',
    maximumGrainSize=(0.5,'kcal/mol'),
    minimumNumberOfGrains=250,
    temperatures=(300,3000,'K',8),
    pressures=(0.001,100,'bar',5),
    interpolation=('Chebyshev', 6, 4),
    completedNetworks=['C2H6'],
)

options(
    units='si',
    generateOutputHTML=False,
    generatePlots=False,
    generatePESDiagrams=True,
)
