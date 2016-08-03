# Data sources
database(
    thermoLibraries = ['primaryThermoLibrary'],
    reactionLibraries = [('Methylformate', False)],
    kineticsFamilies = ['R_Addition_MultipleBond'],
    kineticsEstimator = 'rate rules',
)

# List of species

species(
    label='HCjO',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 C u1 p0 c0 {2,D} {3,S}
2 O u0 p2 c0 {1,D}
3 H u0 p0 c0 {1,S}
        """),
)

species(
    label='CH2O',
    reactive=True,
    structure=adjacencyList(
        """
1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 O u0 p2 c0 {1,D}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
        """),
)

species(
    label='Fmoml',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 C u1 p0 c0 {3,S} {5,S} {6,S}
2 C u0 p0 c0 {3,S} {4,D} {7,S}
3 O u0 p2 c0 {1,S} {2,S}
4 O u0 p2 c0 {2,D}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {2,S}
        """),
)


simpleReactor(
    temperature=(1399,'K'),
    pressure=(1.93,'atm'),
    initialMoleFractions={
        "HCjO": 1,
    },
    terminationConversion={
        'HCjO': 0.999,
    },
    terminationTime=(1e-3,'s'),
)

