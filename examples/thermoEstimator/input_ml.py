database(
    thermoLibraries=['primaryThermoLibrary'],
)

mlEstimator(
    thermo=True,
    # Name of folder containing ML architecture and parameters in database
    name='main',
    # Limits on atom numbers
    minHeavyAtoms=1,
    maxHeavyAtoms=None,
    minCarbonAtoms=0,
    maxCarbonAtoms=None,
    minOxygenAtoms=0,
    maxOxygenAtoms=None,
    minNitrogenAtoms=0,
    maxNitrogenAtoms=None,
    # Limits on cycles
    onlyCyclics=False,
    onlyHeterocyclics=False,    # If onlyHeterocyclics is True, the machine learning estimator is restricted to only
                                # heterocyclics species regardless of onlyCyclics setting.
                                # But onlyCyclics should also be True if onlyHeterocyclics is True.
    minCycleOverlap=0,  # specifies the minimum number of atoms that must be shared between any two cycles
                        # If minCycleOverlap is greater than zero, the machine learning estimator is restricted to
                        # only cyclic species with the specified minimum cyclic overlap regardless of onlyCyclics
                        # setting.
    # If the estimated uncertainty of the thermo prediction is greater than
    # any of these values, then don't use the ML estimate
    H298UncertaintyCutoff=(3.0, 'kcal/mol'),
    S298UncertaintyCutoff=(2.0, 'cal/(mol*K)'),
    CpUncertaintyCutoff=(2.0, 'cal/(mol*K)')
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
    label='spirooctane',
    structure=SMILES('C1CCC2(C1)CCC2'),
)
species(
    label='bicyclooctane',
    structure=SMILES('C1CC2CCC1CC2'),
)
species(
    label='benzene',
    structure=SMILES('c1ccccc1'),
)
species(
    label='naphthalene',
    structure=SMILES('c12ccccc1cccc2'),
)
