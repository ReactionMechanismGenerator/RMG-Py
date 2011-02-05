# Data sources
database(
    thermo_groups='output/RMG_Database/thermo_groups',
    thermo_libraries=[
        'output/RMG_Database/thermo_libraries/primaryThermoLibrary',
        'output/RMG_Database/thermo_libraries/GRI-Mech3.0',
    ],
    kinetics_groups='output/RMG_Database/kinetics_groups',
    seed_mechanisms=['output/RMG_Database/kinetics_libraries/GRI-Mech3.0'],
)

# List of species
species(
    label='O2',     # oxygen
    reactive=True,
    structure=SMILES("O=O"),
)
species(
    label='C8H18i', # isooctane
    reactive=True,
    structure=SMILES("CC(C)CC(C)(C)C"),
)
species(
    label='C2H6On', # ethanol
    reactive=True,
    structure=SMILES("CCO"),
)
species(
    label='C7H8t',  # toluene
    reactive=True,
    structure=SMILES("Cc1ccccc1"),
)
species(
    label='C6H12n', # hex-1-ene
    reactive=True,
    structure=SMILES("CCCCC=C"),
)
species(
    label='Ar',    # argon
    reactive=False,
    structure=SMILES("[Ar]"),
)

# Reaction systems
batchReactor(
    volume=(1.0,'m^3'),
    area=(1.0,'m^2'),
    physicalPropertyModel="idealGas",
    temperatureModel='isothermal',
    pressureModel='isobaric',
    initialConditions={
        "T": (900,'K'),
        "P": (10.0,'bar'),
        "O2": 0.00707000960551,
        "C8H18i": 6.89999461306e-05,
        "C2H6On": 0.0018620014036,
        "C7H8t": 4.80001013417e-05,
        "C6H12n": 2.1000044337e-05,
        "Ar": 0.990929988899,
    },
    reservoirConditions={
        "T": (900,'K'),
        "P": (10.0,'bar'),
        "O2": 0.21,
        "N2": 0.79,
    },
)



termination(
    conversion={
        'C2H6On': 0.9,
    },
    time=(2e2,'s'),
)
   
simulator(
    atol=1e-16,
    rtol=1e-8,
)

model(
    toleranceKeepInEdge=0.0,
    toleranceMoveToCore=0.1,
    toleranceInterruptSimulation=1.0,
    maximumEdgeSpecies=100000,
)

#pressureDependence(
#    method='modified strong collision',
#    grainSize=(2.0,'kcal/mol'),
#    numberOfGrains=200,
#    temperatures=(300,'K',2000,'K',8),
#    pressures=(0.01,'bar',100,'bar',5),
#    interpolation=('Chebyshev', 4, 4),
#)

options(
    units='si',
    saveRestart='daily',
    drawMolecules=False,
    generatePlots=False,
)
