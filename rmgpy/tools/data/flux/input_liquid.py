# Data sources
database(
    thermoLibraries = ['primaryThermoLibrary'],
    reactionLibraries = [],
    seedMechanisms = [],
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
    label='heptane',
    reactive=True,
    structure=SMILES("CCCCCCC"),
)

# Reaction systems
liquidReactor(
    temperature=(700,'K'),
    initialConcentrations={
        "ethane": (1.0e-3, 'mol/cm^3'),
    },
    terminationTime=(1e6,'s'),
    constantSpecies=['ethane']
)

solvation(
    solvent='heptane'
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
    generatePlots=False,
    saveSimulationProfiles=True,
)
