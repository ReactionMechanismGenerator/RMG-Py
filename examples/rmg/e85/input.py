# Data sources
database(
    thermoLibraries = ['primaryThermoLibrary', 'GRI-Mech3.0'],
    reactionLibraries = [],
    seedMechanisms = ['GRI-Mech3.0'],
    kineticsDepositories = ['training'],
    kineticsFamilies = ['!Intra_Disproportionation','!Substitution_O'],
    kineticsEstimator = 'rate rules',
)

# Constraints on generated species
generatedSpeciesConstraints(
    allowed=['seed mechanisms'],
    maximumRadicalElectrons = 2,
)

# List of species
species(
    label='O2',     # oxygen
    reactive=True,
    structure=SMILES("[O][O]"),
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
simpleReactor(
    temperature=(1000,'K'),
    pressure=(1.0,'bar'),
    initialMoleFractions={
        "O2": 0.00707000960551,
        "C8H18i": 6.89999461306e-05,
        "C2H6On": 0.0018620014036,
        "C7H8t": 4.80001013417e-05,
        "C6H12n": 2.1000044337e-05,
        "Ar": 0.990929988899,
    },
    terminationConversion = {'C2H6On': 0.9},
    terminationTime = (2e2,'s'),
)
   
simulator(
    atol=1e-16,
    rtol=1e-8,
)

model(
    toleranceKeepInEdge=0.0,
    toleranceMoveToCore=0.01,
    toleranceInterruptSimulation=1.0,
    maximumEdgeSpecies=1000000,
)

#pressureDependence(
#    method='reservoir state',
#    grainSize=(2.0,'kcal/mol'),
#    numberOfGrains=200,
#    temperatures=(300,2000,'K',8),
#    pressures=(0.01,100,'bar',5),
#    interpolation=('Chebyshev', 6, 4),
#)

options(
    units='si',
    saveRestartPeriod=(1,'day'),
    drawMolecules=False,
    generatePlots=False,
)

