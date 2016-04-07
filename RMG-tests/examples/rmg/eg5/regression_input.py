
options(
    title = 'EG5',
    tolerance = 0.05
)

observable(
    label='n-heptane',
    structure=SMILES("CCCCCCC"),
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
        "n-heptane": 0.02,
        "Ar": 0.98,
    }],
    temperatures=([1500],'K'),
    pressures=([400.],'Pa'),
)