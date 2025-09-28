database(
    thermoLibraries = ['primaryThermoLibrary', 'GRI-Mech3.0']   
)

species(
    label='DIPK',
    structure=SMILES("CC(C)C(=O)C(C)C"),
)
species(
    label='O2',
    structure=SMILES("[O][O]"),
)
species(
    label='R_tert',
    structure=SMILES("CC(C)C(=O)[C](C)C"),
)
species(
    label='R_pri',
    structure=SMILES("CC(C)C(=O)C(C)[CH2]"),
)
species(
    label='Cineole',
    structure=SMILES('CC12CCC(CC1)C(C)(C)O2'),
)

quantumMechanics(
    software='mopac',
    method='pm3',
    onlyCyclics = True,
    maxRadicalNumber = 0,
)
