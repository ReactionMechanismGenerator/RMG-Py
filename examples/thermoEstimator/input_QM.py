database(
    thermoLibraries = ['primaryThermoLibrary', 'GRI-Mech3.0']   
)

species(
    label='Cineole',
    structure=SMILES('CC12CCC(CC1)C(C)(C)O2'),
)

species(
    label='JP10',
    structure=SMILES('C1CC2C3CCC(C3)C2C1'),
)

species(
    label='chx',
    structure=SMILES('C1CCCCC1'),
)

species(
    label='chpent',
    structure=SMILES('C1CCCC1'),
)

species(
    label='chb',
    structure=SMILES('C1CCC1'),
)

species(
    label='chprop',
    structure=SMILES('C1CC1'),
)

species(
    label='chhept',
    structure=SMILES('C1CCCCCC1'),
)


quantumMechanics(
    software='mopac',
    method='pm3',
    onlyCyclics = True,
    maxRadicalNumber = 0,
)
