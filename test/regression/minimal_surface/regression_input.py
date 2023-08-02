options(
    title = 'minimal_surface',
    tolerance = 0.5,
)

# observables
observable(
	label = 'CH4',
    structure=SMILES('C'),
)

observable(
	label = 'O2',
    structure=SMILES('[O][O]'),
)

observable(
	label = 'X',
    structure=adjacencyList(
        """
1 X u0 p0 c0
"""),
)


# List of species
species(
    label='CH4',
    structure=SMILES("[CH4]"),
)
species(
   label='O2',
   structure=adjacencyList(
       """
1 O u1 p2 c0 {2,S}
2 O u1 p2 c0 {1,S}
"""),
)
species(
    label='N2',
    structure=SMILES("N#N"),
)
species(
	label='X',
    structure=adjacencyList(
        """
1 X u0 p0 c0
"""),
)


# reactor setups
reactorSetups(
    reactorTypes=['IdealGasReactor'],
    terminationTimes=([1.0], 's'),
    initialMoleFractionsList=[{
        'CH4': 0.15,
        'O2': 0.15,
        'N2': 0.7,
    }],
    initialSurfaceMoleFractionsList=[{
        'X': 1.0,
    }],

    temperatures=([1000.0], 'K'),
    pressures=([1.0], 'bar'),
)
