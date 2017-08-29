
options(
    title = 'EG7',
    tolerance = 0.05
)

observable(
    label='ethane',
    structure=SMILES("CC"),
)

# reactor setups
reactorSetups(
    reactorTypes=['IdealGasReactor'],
    terminationTimes=([1e-1],'s'),
    initialMoleFractionsList=[{
        "ethane": 1.0,
    }],
    temperatures=([1350],'K'),
    pressures=([1.],'bar'),
)