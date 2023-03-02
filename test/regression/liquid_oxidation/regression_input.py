options(
    title='liquid_oxidation',
    tolerance=0.1
)

observable(
    label='pentane',
    structure=SMILES('CCCCC')
)

species(
    label='oxygen',
    structure=SMILES('[O][O]'),
)

reactorSetups(
    reactorTypes=['IdealGasReactor'],
    terminationTimes=([1e3], 's'),
    initialMoleFractionsList=[{
        'pentane': 0.9,
        'oxygen': 0.1,
    }],
    temperatures=([600], 'K'),
    pressures=([1.0], 'bar'),
)
