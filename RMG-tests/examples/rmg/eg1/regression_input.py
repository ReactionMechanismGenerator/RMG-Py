
options(
    title = 'EG1',
    tolerance = 0.05
)

# observables
observable(
	label = 'ethane',
    structure=adjacencyList(
        """
        1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
		2 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
		3 H u0 p0 c0 {1,S}
		4 H u0 p0 c0 {1,S}
		5 H u0 p0 c0 {1,S}
		6 H u0 p0 c0 {2,S}
		7 H u0 p0 c0 {2,S}
		8 H u0 p0 c0 {2,S}
        """),
)

observable(
	label = 'methyl',
    structure=SMILES("[CH3]"),
)

# species definition used in the reactor setup specification
species(
	label = 'argon',
    structure=SMILES("[Ar]"),
)

# reactor setups
reactorSetups(
    reactorTypes=['IdealGasReactor'],
    terminationTimes=([1e3],'s'),
    initialMoleFractionsList=[{
        "ethane": .95,
        "argon": 0.05,
    }],
    temperatures=([1350],'K'),
    pressures=([1.0],'bar'),
)