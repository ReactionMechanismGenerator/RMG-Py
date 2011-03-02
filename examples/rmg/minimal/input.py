# Data sources
database(
    thermo_groups='output/RMG_database/thermo_groups',
    thermo_libraries=[
        'output/RMG_database/thermo_libraries/primaryThermoLibrary',
        #'output/RMG_database/thermo_libraries/GRI-Mech3.0',
    ],
    kinetics_groups='output/RMG_database/kinetics_groups',
    #seed_mechanisms=[
    #    'output/RMG_database/kinetics_libraries/GRI-Mech3.0',
    #],
    frequencies_groups='output/RMG_database/frequencies_groups',
)

# List of species
species(
    label='ethane',
    reactive=True,
    structure=SMILES("CC"),
)

# Reaction systems
simpleReactor(
    temperature=(1350,'K'),
    pressure=(1.0,'bar'),
    initialMoleFractions={
        "ethane": 1.0,
    },
)

termination(
    conversion={
        'ethane': 0.9,
    },
    time=(1e6,'s'),
)

simulator(
    atol=1e-16,
    rtol=1e-8,
)

model(
    toleranceKeepInEdge=0.0,
    toleranceMoveToCore=0.1,
    toleranceInterruptSimulation=1.0,
    maximumEdgeSpecies=100000
)

options(
    units='si',
    saveRestart=True,
    drawMolecules=False,
    generatePlots=False,
)
