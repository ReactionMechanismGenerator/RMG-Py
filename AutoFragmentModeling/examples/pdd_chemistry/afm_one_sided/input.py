database(
    thermoLibraries = ['primaryThermoLibrary', 'SulfurLibrary','DFT_QCI_thermo'],
    reactionLibraries = [],
    seedMechanisms = [],
    kineticsDepositories = 'default',
    kineticsFamilies = ['H_Abstraction','R_Addition_MultipleBond','R_Recombination',
    'Disproportionation'],
    kineticsEstimator = 'rate rules',
)

generatedSpeciesConstraints(
    maximumRadicalElectrons = 1,
)

species(
    label='arccccr',
    reactive=True,
    structure=fragment_adj("""1  C u0 p0 c0 {2,S} {3,S} {11,S} {12,S}
2  C u0 p0 c0 {1,S} {4,S} {13,S} {14,S}
3  C u0 p0 c0 {1,S} {5,S} {15,S} {16,S}
4  C u0 p0 c0 {2,S} {17,S} {18,S} {19,S}
5  C u0 p0 c0 {3,S} {6,S} {7,D}
6  C u0 p0 c0 {5,S} {8,D} {20,S}
7  C u0 p0 c0 {5,D} {10,S} {24,S}
8  C u0 p0 c0 {6,D} {9,S} {21,S}
9  C u0 p0 c0 {8,S} {10,D} {22,S}
10 C u0 p0 c0 {7,S} {9,D} {23,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 R u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {7,S}
"""),
)

species(
    label='rccccr',
    reactive=True,
    structure=fragment_adj("""1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {10,S} {13,S} {14,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  R u0 p0 c0 {3,S}
10 R u0 p0 c0 {4,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
"""),
)

species(
    label='rcccc',
    reactive=True,
    structure=fragment_adj("""1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {9,S} {13,S} {14,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  R u0 p0 c0 {4,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
"""),
)


simpleReactor(
    temperature = (623,"K"),  # 350 C
    pressure = (1050,"bar"),
    initialMoleFractions={
        "arccccr": 1,
        "rccccr": 1,
        "rcccc": 1
    },
    terminationTime = (1000,"h"),
    terminationConversion = {
        "arccccr": 0.6,
    },
)

simulator(
    atol = 1e-16,
    rtol = 1e-08,
    sens_atol = 1e-06,
    sens_rtol = 0.0001,
)

model(
    toleranceKeepInEdge = 0.05,
    toleranceMoveToCore = 0.2,
    toleranceInterruptSimulation = 1000000000,
    filterReactions=True,
)

options(
    units = "si",
    saveRestartPeriod = None,
    generateOutputHTML = True,
    generatePlots = False,
    saveSimulationProfiles = False,
    saveEdgeSpecies = True,
    verboseComments = False,
)
