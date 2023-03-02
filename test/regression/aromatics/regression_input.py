options(
    title = 'Aromatics',
    tolerance = 0.5,
)

# observables
observable(
    label = 'benzene',
    structure=SMILES('c1ccccc1'),
)


# species definition used in the reactor setup specification
species(
    label = 'benzene',
    structure=SMILES('c1ccccc1'),
)

species(
    label = 'ethyne',
    structure=SMILES('C#C'),
)

# reactor setups
reactorSetups(
    reactorTypes=['IdealGasReactor'],
    terminationTimes=([100],'s'),
    initialMoleFractionsList=[{
        "benzene": 0.2,
        "ethyne": 0.8,
    }],
    temperatures=([1500.0],'K'),
    pressures=([1.0],'bar'),
)
