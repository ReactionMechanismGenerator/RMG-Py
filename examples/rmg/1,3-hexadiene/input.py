# Data sources
database(
    thermoLibraries = ['primaryThermoLibrary', 'GRI-Mech3.0'],
    reactionLibraries = [],
    seedMechanisms = [],
    kineticsDepositories = ['training'], #  'all', 'default'==['training'], [], 
    kineticsFamilies = ['!Intra_Disproportionation'],
    kineticsEstimator = 'rate rules',
)

# List of species
species(
    label='HXD13',
    reactive=True,
    structure=CML(
        """
        <molecule>
            <atomArray>
                <atom id="a1" elementType="C" />
                <atom id="a2" elementType="C" />
                <atom id="a3" elementType="C" />
                <atom id="a4" elementType="C" />
                <atom id="a5" elementType="C" />
                <atom id="a6" elementType="C" />
            </atomArray>
            <bondArray>
                <bond atomRefs2="a1 a2" order="D" />
                <bond atomRefs2="a2 a3" order="S" />
                <bond atomRefs2="a3 a4" order="D" />
                <bond atomRefs2="a4 a5" order="S" />
                <bond atomRefs2="a5 a6" order="S" />
            </bondArray>
        </molecule>
        """),
)
species(
    label='CH4',
    reactive=True,
    structure=SMILES("C"),
)
species(
    label='H2',
    reactive=True,
    structure=adjacencyList(
        """
        1 H 0 {2,S}
        2 H 0 {1,S}
        """),
)
species(
    label='N2',
    reactive=False,
    structure=InChI("InChI=1/N2/c1-2"),
)

# Reaction systems
simpleReactor(
    temperature=(1350,'K'),
    pressure=(1.0,'bar'),
    initialMoleFractions={
        "HXD13": 6.829e-4,
        "CH4": 0.104,
        "H2": 0.0156,
        "N2": 0.8797,
    },
    terminationConversion={
        'HXD13': 0.9,
    },
    terminationTime=(1e0,'s'),
)

simulator(
    atol=1e-16,
    rtol=1e-8,
)

model(
    toleranceKeepInEdge=0.0,
    toleranceMoveToCore=0.5,
    toleranceInterruptSimulation=0.5,
    maximumEdgeSpecies=100000
)

quantumMechanics(
    software='mopac',
    fileStore='QMfiles', # relative to where you run it? defaults to inside the output folder.
    scratchDirectory = None, # not currently used
    onlyCyclics = True,
    maxRadicalNumber = 0,
    )

pressureDependence(
    method='modified strong collision',
    maximumGrainSize=(0.5,'kcal/mol'),
    minimumNumberOfGrains=250,
    temperatures=(300,2000,'K',8),
    pressures=(0.01,100,'bar',5),
    interpolation=('Chebyshev', 6, 4),
)

options(
    units='si',
    saveRestartPeriod=(1,'hour'),
    drawMolecules=False,
    generatePlots=False,
)
