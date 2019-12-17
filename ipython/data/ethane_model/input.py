# Data sources
database(
    thermoLibraries = ['DFT_QCI_thermo', 'primaryThermoLibrary'],
    reactionLibraries = [],
    seedMechanisms = [],
    kineticsDepositories = ['training'],
    kineticsFamilies = 'default',
    kineticsEstimator = 'rate rules',
)

generatedSpeciesConstraints(
    allowed=['input species','seed mechanisms','reaction libraries'],
    maximumRadicalElectrons = 2,
    maximumCarbonAtoms = 10,
)

# List of species
species(
    label='ethane',
    reactive=True,
    structure=SMILES("CC"),
)

species(
    label='methane',
    reactive=True,
    structure=SMILES("C"),
)

# Reaction systems
simpleReactor(
    temperature=(1300,'K'),
    pressure=(1.0,'bar'),
    initialMoleFractions={
        "ethane": 1.0,
    },
    terminationTime=(0.5,'ms'),
    sensitivity=['ethane','methane']
)

simulator(
    atol=1e-16,
    rtol=1e-8,
)

model(
    toleranceMoveToCore=0.01,
    filterReactions=True,
)

options(
    units='si',
    generateOutputHTML=False,
    generatePlots=False,
    saveEdgeSpecies=True,
    saveSimulationProfiles=True,
)

