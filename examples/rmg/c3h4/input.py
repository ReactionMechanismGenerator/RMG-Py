# Data sources
database(
	thermoLibraries = ['primaryThermoLibrary', 'GRI-Mech3.0-N'],
	reactionLibraries = [],
	seedMechanisms = ['GRI-Mech3.0-N'],
	kineticsDepositories = ['training'], 
	kineticsFamilies = 'default',
	kineticsEstimator = 'rate rules',
)

generatedSpeciesConstraints(
    allowed = ['seed mechanisms'],
    maximumRadicalElectrons = 4,
)

# List of species
species(
	label='CH2',
	reactive=True,
	structure=SMILES("[CH2]"),
)
species(
	label='C2H2',
	reactive=True,
	structure=SMILES("C#C"),
)
species(
	label='N2',
	reactive=False,
	structure=InChI("InChI=1/N2/c1-2"),
)

# Reaction systems
simpleReactor(
	temperature=(1350,'K'),
	pressure=(1.0,'bar'),
	initialMoleFractions={
		"CH2": 0.001,
		"C2H2": 0.099,
		"N2": 0.9,
	},
	terminationConversion={
		'C2H2': 0.9,
	},
	terminationTime=(1e12,'s'),
)

simulator(
	atol=1e-16,
	rtol=1e-8,
)

model(
	toleranceKeepInEdge=1e-5,
	toleranceMoveToCore=0.01,
	toleranceInterruptSimulation=0.1,
	maximumEdgeSpecies=9999999
)

options(
	units='si',
	saveRestartPeriod=None,
	generateOutputHTML=False,
	generatePlots=False,
)
