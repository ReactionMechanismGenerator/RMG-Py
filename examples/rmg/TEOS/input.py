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
	label='TEOS',
	reactive=True,
	structure=InChI("InChI=1/C8H20O4Si/c1-5-9-13(10-6-2,11-7-3)12-8-4/h5-8H2,1-4H3"),
)
species(
	label='EtOH',
	reactive=True,
	structure=InChI("InChI=1/C2H6O/c1-2-3/h3H,2H2,1H3"),
)
species(
	label='Ar',
	reactive=False,
	structure=InChI("InChI=1/Ar"),
)

# Reaction systems
simpleReactor(
	temperature=(1500,'K'),
	pressure=(1.0,'bar'),
	initialMoleFractions={
		"TEOS": 0.1,
		"EtOH": 0.3,
		"Ar": 0.6,
	},
	terminationConversion={
		'TEOS': 0.5,
	},
	terminationTime=(1e-3,'s'),
)

simulator(
	atol=1e-16,
	rtol=1e-8,
)

model(
	toleranceKeepInEdge=1e-4,
	toleranceMoveToCore=0.1,
	toleranceInterruptSimulation=1.0,
	maximumEdgeSpecies=9999999
)

options(
	units='si',
	saveRestartPeriod=None,
	drawMolecules=False,
	generatePlots=False,
)
