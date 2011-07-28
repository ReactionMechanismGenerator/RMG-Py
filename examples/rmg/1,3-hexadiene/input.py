# Data sources
database(
    '/$RMGpy/RMG-database/input',
    thermoLibraries = ['primaryThermoLibrary', 'GRI-Mech3.0'],
    reactionLibraries = ['GRI-Mech3.0'],
    seedMechanisms = ['Glarborg/C1'],
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
)

termination(
    conversion={
        'HXD13': 0.9,
    },
    time=(1e0,'s'),
)

simulator(
    atol=1e-16,
    rtol=1e-8,
)

model(
    toleranceKeepInEdge=0.0,
    toleranceMoveToCore=0.1,
    toleranceInterruptSimulation=1.0,
    maximumEdgeSpecies=100000
)

pressureDependence(
    method='modified strong collision',
    maximumGrainSize=(2.0,'kcal/mol'),
    minimumNumberOfGrains=200,
    temperatures=(300,2000,'K',8),
    pressures=(0.01,100,'bar',5),
    interpolation=('Chebyshev', 6, 4),
)

options(
    units='si',
    saveRestart=False,
    drawMolecules=False,
    generatePlots=False,
)
