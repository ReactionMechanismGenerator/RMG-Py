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
species(
   label='N2',
   reactive=False,
   structure=SMILES("N#N"),
)

# Reaction systems
simpleReactor(
    temperature=(1000,'K'),
    pressure=(1.0,'bar'),
    initialMoleFractions={
        'H2':.57, 'O2':.33,
        'N2':.1,
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
    saveRestartPeriod=None,
    generateOutputHTML=True,
    generatePlots=True,
    saveEdgeSpecies=True,
    saveSimulationProfiles=True,
    generateSeedEachIteration=True,
    saveSeedToDatabase=True,
    verboseComments=True,
)

pressureDependence(
        method='modified strong collision',
        maximumGrainSize=(0.5,'kcal/mol'),
        minimumNumberOfGrains=250,
        temperatures=(300,2000,'K',8),
        pressures=(0.01,100,'bar',5),
        interpolation=('Chebyshev', 6, 4),
        maximumAtoms=16,
)

quantumMechanics(
        software='mopac',
        method='pm3',
        fileStore='QMfiles',
        scratchDirectory = None,
        onlyCyclics = False,
        maxRadicalNumber = 0,
        )

