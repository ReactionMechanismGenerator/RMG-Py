# Data sources
database(
    thermoLibraries = ['primaryThermoLibrary','DFT_QCI_thermo','GRI-Mech3.0'],
    reactionLibraries = [('Methylformate',False),('Glarborg/highP',False)],
    seedMechanisms = ['Glarborg/C2'],
    kineticsDepositories = ['training'],
    kineticsFamilies = ['!Intra_Disproportionation'],
    kineticsEstimator = 'rate rules',
)

# List of species
species(
    label='Mfmt',
    reactive=True,
    structure=SMILES("COC=O"),
)
species(
    label='O2',
    reactive=True,
    structure=SMILES("[O][O]"),
)
species(
    label='C2H',
    reactive=True,
    structure=SMILES("C#[C]"),
)
species(
    label='CH',
    reactive=True,
    structure=adjacencyList(
        """
        1     C     3 {2,S}
        2     H     0 {1,S}
        """),
)
species(
    label='H2O',
    reactive=True,
    structure=SMILES("O"),
)
species(
    label='H2',
    reactive=True,
    structure=SMILES("[H][H]"),
)
species(
    label='CO',
    reactive=True,
    structure=SMILES("[C]=O"),
)
species(
    label='CO2',
    reactive=True,
    structure=SMILES("C(=O)=O"),
)
species(
    label='CH4',
    reactive=True,
    structure=SMILES("C"),
)
species(
    label='CH3',
    reactive=True,
    structure=SMILES("[CH3]"),
)
species(
    label='CH3OH',
    reactive=True,
    structure=SMILES("CO"),
)
species(
    label='C2H4',
    reactive=True,
    structure=SMILES("C=C"),
)
species(
    label='C2H2',
    reactive=True,
    structure=SMILES("C#C"),
)
species(
    label='CH2O',
    reactive=True,
    structure=SMILES("C=O"),
)
species(
    label='CH3CHO',
    reactive=True,
    structure=SMILES("CC=O"),
)


# Bath gas
species(
    label='Ar',
    reactive=False,
    structure=InChI("InChI=1S/Ar"),
)

# Reaction systems
simpleReactor(
    temperature=(650,'K'),
    pressure=(1.0,'bar'),
    initialMoleFractions={
        "Mfmt": 0.01,
        "O2": 0.02,
        "Ar": 0.08,
    },
    terminationTime=(0.5,'s'),
)
simpleReactor(
    temperature=(1350,'K'),
    pressure=(3.0,'bar'),
    initialMoleFractions={
        "Mfmt": 0.01,
        "O2": 0.02,
        "Ar": 0.97,
    },
    terminationTime=(0.5,'s'),
)
simpleReactor(
    temperature=(1950,'K'),
    pressure=(10.0,'bar'),
    initialMoleFractions={
        "Mfmt": 0.01,
        "O2": 0.02,
        "Ar": 0.97,
    },
    terminationTime=(0.5,'s'),
)

simulator(
    atol=1e-22,
    rtol=1e-8,
)

model(
    toleranceKeepInEdge=0.0,
    toleranceMoveToCore=0.0005,
    toleranceInterruptSimulation=1.0,
    maximumEdgeSpecies=100000
)

pressureDependence(
    method='modified strong collision', # 'reservoir state'
    maximumGrainSize=(1.0,'kcal/mol'),
    minimumNumberOfGrains=200,
    temperatures=(290,3500,'K',8),
    pressures=(0.02,100,'bar',5),
    interpolation=('Chebyshev', 6, 4),
    maximumAtoms=16,
)

options(
    units='si',
    saveRestartPeriod=(1,'hour'),
    drawMolecules=False,
    generatePlots=False,
    saveConcentrationProfiles=False,
)
