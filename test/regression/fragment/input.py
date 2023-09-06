database(
    thermoLibraries=["primaryThermoLibrary", "SulfurLibrary", "DFT_QCI_thermo"],
    reactionLibraries=[],
    seedMechanisms=[],
    kineticsDepositories=["training"],
    kineticsFamilies=[
        "H_Abstraction",
        "R_Addition_MultipleBond",
        "R_Recombination",
        "Disproportionation",
    ],
    kineticsEstimator="rate rules",
)

generatedSpeciesConstraints(
    maximumRadicalElectrons=1,
)

species(
    label="LCCCC",
    reactive=True,
    structure=fragment_adj(
        """1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {12,S} {13,S} {14,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  L u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
"""
    ),
)

species(
    label="RCCCC",
    reactive=True,
    structure=fragment_adj(
        """1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {12,S} {13,S} {14,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  R u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
"""
    ),
)


simpleReactor(
    temperature=(1600, "K"),  # Temp referance to rmg example input file
    pressure=(40, "bar"),
    initialMoleFractions={"LCCCC": 1, "RCCCC": 1},
    terminationTime=(1e6, "s"),
    terminationConversion={
        "RCCCC": 0.99,
    },
)

simulator(
    atol=1e-16,
    rtol=1e-08,
)

model(
    #    toleranceKeepInEdge = 0.05,
    toleranceMoveToCore=0.01,
    toleranceInterruptSimulation=0.01,
    maxNumSpecies=10,
    filterReactions=True,
)

options(
    units="si",
    saveRestartPeriod=None,
    generateOutputHTML=False,
    generatePlots=False,
    saveSimulationProfiles=False,
    saveEdgeSpecies=False,
    verboseComments=False,
)
