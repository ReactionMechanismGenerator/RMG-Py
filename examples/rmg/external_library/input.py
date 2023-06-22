# Sometimes it is useful to load an RMG thermo or kinetic library that is not within the RMG-database repository.
# One use case could be, e.g., running RMG on a server where it is installed in a shared folder to which users have a read-only access.
# Here we show an example of running RMG with a kinetic library external to its database.
# THe same approach can be applied for kinetic seed mechanisms and for thermodynamic libraries.

# Data sources
database(
	thermoLibraries = ['primaryThermoLibrary'],
	reactionLibraries = ['custom_library'],  # Specify a FULL PATH if needed. In this case, the library is located in the run folder, and a relative path (just the folder name) is enough.
	seedMechanisms = [],
	kineticsDepositories = ['training'], 
	kineticsFamilies = 'default',
	kineticsEstimator = 'rate rules',
)

# List of species
species(
	label='CH2',
	structure=SMILES("[CH2]"),
)
species(
	label='C2H2',
	structure=SMILES("C#C"),
)
species(
	label='N2',
	reactive=False,
	structure=SMILES("N#N"),
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
	terminationConversion={'C2H2': 0.5},
	terminationTime=(0.5,'s'),
)
model(
	toleranceKeepInEdge=0,
	toleranceMoveToCore=0.05,
	toleranceInterruptSimulation=0.05,
	maximumEdgeSpecies=100000
)

simulator(
	atol=1e-16,
	rtol=1e-8,
)

options(
	units='si',
	generateOutputHTML=False,
	generatePlots=False,
)

generatedSpeciesConstraints(
    allowed = ['seed mechanisms', 'reaction libraries'],
    maximumRadicalElectrons = 2,
)

