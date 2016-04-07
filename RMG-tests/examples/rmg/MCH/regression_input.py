
options(
    title = 'MCH',
    tolerance = 0.05
)

observable(
    label='MCH',
    structure=SMILES("CC1CCCCC1"),
)

species(
    label='O2',
    structure=adjacencyList(
        """
        multiplicity 3
        1 O u1 p2 c0 {2,S}
        2 O u1 p2 c0 {1,S}
        """),
)

species(
    label='N2',
    structure=SMILES("N#N"),
)

# reactor setups
reactorSetups(
    reactorTypes=['IdealGasReactor'],
    terminationTimes=([100],'s'),
    initialMoleFractionsList=[{
        "MCH": 0.01,
        "O2": 0.105,
        "N2": 0.885,
    }],
    temperatures=([2000],'K'),
    pressures=([20.],'bar'),
)