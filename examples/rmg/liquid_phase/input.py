# Data sources
database(
    thermoLibraries = ['primaryThermoLibrary'],
    reactionLibraries = [],
    seedMechanisms = [],
    kineticsDepositories = ['training'],
    kineticsFamilies = ['!Intra_Disproportionation'],
    kineticsEstimator = 'rate rules',
)

# List of species
species(
    label='octane',
    reactive=True,
    structure=SMILES("C(CCCCC)CC"),
)

species(
    label='oxygen',
    reactive=True,
    structure=SMILES("[O][O]"),
)

# Reaction systems
liquidReactor(
    temperature=(500,'K'),
    initialConcentrations={
        "octane": (6.154e-3,'mol/cm^3'),
        "oxygen": (4.953e-6,'mol/cm^3')
    },
    terminationConversion={
        'octane': 0.9,
    },
    terminationTime=(1e6,'s'),
)

solvation(
	solvent='octane',
	viscosity = (1e-3, "cP"),
	temperature = (298, "K"), # T at which diffusion limits on reaction rates evaluated for chemkin output
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
    drawMolecules=False,
    generatePlots=False,
    saveConcentrationProfiles=True,
)
