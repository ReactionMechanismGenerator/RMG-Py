options(
    title = 'Oxidation',
    tolerance = 0.5,
)

# observables
observable(
	label = 'OH',
    structure=SMILES('[OH]'),
)


# species definition used in the reactor setup specification
species(
	label = 'OH',
	structure=SMILES('[OH]'),
)

species(
	label = 'N2',
    structure=SMILES("N#N"),
)

species(
	label = 'O2',
	structure=SMILES('[O][O]'),
)

species(
	label = 'propane',
	structure=SMILES('CCC'),
)

# reactor setups
reactorSetups(
    reactorTypes=['IdealGasReactor'],
    terminationTimes=([100.0],'s'),
    initialMoleFractionsList=[{
        "propane": 2.0/7.0,
        "O2": 1.0,
        "N2":4.0,
    }],
    temperatures=([725.0],'K'),
    pressures=([10.0],'bar'),
)
