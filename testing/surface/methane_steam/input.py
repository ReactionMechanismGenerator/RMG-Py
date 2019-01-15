# Data sources
database(
    thermoLibraries=['surfaceThermo', 'primaryThermoLibrary'],
    reactionLibraries = [],
    seedMechanisms = [],
    kineticsDepositories = ['training'],
    kineticsFamilies = 'default',
    kineticsEstimator = 'rate rules',
)

# List of species
species(
    label='methyl',
    reactive=True,
    structure=SMILES("[CH3]"),
)

species(
    label='site',
    reactive=True,
    structure=adjacencyList("1 X u0"),
)

# Reaction systems
surfaceReactor(
    temperature=(1350,'K'),
    initialPressure=(1.0, 'bar'),
    initialGasMoleFractions={
        "methyl": 1.0,
    },
    initialSurfaceCoverages={
        "site": 1.0,
    },
    surfaceVolumeRatio = (10., 'm^-1'),
    surfaceSiteDensity=(2.9e-9, 'mol/cm^2'),
    terminationTime=(1e-3,'s'),
)

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
    generateOutputHTML=True,
    generatePlots=False,
    saveEdgeSpecies=True,
    saveSimulationProfiles=True,
    verboseComments=True,
)
