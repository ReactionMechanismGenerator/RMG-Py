
options(
    title = 'NC',
    tolerance = 0.05
)

observable(
    label='NC',
    structure=SMILES("NC"),
)

species(
    label='O2',
    structure=SMILES("[O][O]"),
)

species(
    label='Ar',
    structure=adjacencyList(
                """
                1 Ar u0 p4 c0
                """),
)


# reactor setups
reactorSetups(
    reactorTypes=['IdealGasReactor'],
    terminationTimes=([0.002],'s'),
    initialMoleFractionsList=[{
        "NC": 0.0005,
        "O2": 0.002,
        "Ar": 0.9975,
    }],
    temperatures=([1500],'K'),
    pressures=([2.],'atm'),
)