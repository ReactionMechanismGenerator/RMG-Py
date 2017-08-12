database(
    thermoLibraries = ['primaryThermoLibrary'],
    reactionLibraries = [],
    seedMechanisms = [],
    kineticsDepositories = 'default',
    kineticsFamilies = 'default',
    kineticsEstimator = 'rate rules',
)

species(
    label = "oxygen",
    reactive = True,
    structure = adjacencyList(
"""
multiplicity 3
1 O u1 p2 c0 {2,S}
2 O u1 p2 c0 {1,S}
"""),
)

species(
    label = "hexane",
    reactive = True,
    structure = adjacencyList(
"""
1  C u0 p0 c0 {2,S} {7,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {3,S} {10,S} {11,S}
3  C u0 p0 c0 {2,S} {4,S} {12,S} {13,S}
4  C u0 p0 c0 {3,S} {5,S} {14,S} {15,S}
5  C u0 p0 c0 {4,S} {6,S} {16,S} {17,S}
6  C u0 p0 c0 {5,S} {18,S} {19,S} {20,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
"""),
)

liquidReactor(
    temperature = (400,"K"),
    initialConcentrations={
        "hexane": (7.598e-3,'mol/cm^3'),
        "oxygen": (4.953e-6,'mol/cm^3'),
    },
    terminationTime = (10,"s"),
)

solvation(
    solvent='hexane'
)

simulator(
    atol = 1e-16,
    rtol = 1e-08,
)

model(
    toleranceMoveToCore = 0.001,
    toleranceKeepInEdge = 1e-9,
    toleranceInterruptSimulation = 0.1,
    maximumEdgeSpecies = 100000,
    minCoreSizeForPrune = 50,
    minSpeciesExistIterationsForPrune = 2,
    filterReactions = 0,
)

options(
    units = "si",
    saveRestartPeriod = None,
    generateOutputHTML = False,
    generatePlots = False,
    saveSimulationProfiles = False,
    saveEdgeSpecies = True,
    verboseComments = False,
)