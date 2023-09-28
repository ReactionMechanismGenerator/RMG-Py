
options(
    title='SO2',
    tolerance=0.1
)

observable(
    label='SO2',
    structure=SMILES("O=S=O"),
)
species(
    label='H2S',
    structure=SMILES("S"),
)

species(
    label='O2',
    structure=SMILES("[O][O]"),
)

species(
    label='N2',
    structure=SMILES("N#N"),
)

reactorSetups(
    reactorTypes=['IdealGasReactor'],
    terminationTimes=([0.01],'s'),
    initialMoleFractionsList=[{
        "H2S": 0.000756,
        "O2": 0.001290,
        "N2": 0.997954}],
    temperatures=([900],'K'),
    pressures=([30.],'bar'),
)
