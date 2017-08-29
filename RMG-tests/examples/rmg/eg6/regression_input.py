
options(
    title = 'EG6',
    tolerance = 0.05
)

observable(
    label='ethane',
    structure=SMILES("CC"),
)

species(
    label='O2',
    structure=SMILES("[O][O]"),
)

# species definition used in the reactor setup specification
species(
	label = 'Ar',
    structure=SMILES("[Ar]"),
)

# reactor setups
reactorSetups(
    reactorTypes=['IdealGasReactor'],
    terminationTimes=([1e0],'s'),
    initialMoleFractionsList=[{
        "ethane": 1,
        "O2": 3.5,
        "Ar": 1,
    }],
    temperatures=([1500],'K'),
    pressures=([2.],'bar'),
)