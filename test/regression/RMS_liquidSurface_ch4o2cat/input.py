# Data sources
database(
    thermoLibraries=[
        "surfaceThermoPt111",
        "primaryThermoLibrary",
        "thermo_DFT_CCSDTF12_BAC",
        "DFT_QCI_thermo",
    ],  # 'surfaceThermoPt' is the default. Thermo data is derived using bindingEnergies for other metals
    reactionLibraries=[
        ("Surface/CPOX_Pt/Deutschmann2006_adjusted", False)
    ],  # when Ni is used change the library to Surface/Deutschmann_Ni
    seedMechanisms=[],
    kineticsDepositories=["training"],
    kineticsFamilies=["surface", "default"],
    kineticsEstimator="rate rules",
)

catalystProperties(metal="Pt111")

species(
    label="CH4",
    reactive=True,
    structure=SMILES("[CH4]"),
)

species(
    label="O2",
    reactive=True,
    structure=adjacencyList(
        """
1 O u1 p2 c0 {2,S}
2 O u1 p2 c0 {1,S}
"""
    ),
)

species(
    label="pentane",
    reactive=False,
    structure=SMILES("CCCCC"),
)

species(
    label="vacantX",
    reactive=True,
    structure=adjacencyList("1 X u0"),
)

# Reaction systems
liquidSurfaceReactor(
    temperature=(600, "K"),
    initialConcentrations={
        "CH4": (1.0e-4, "mol/cm^3"),
        "O2": (1.0e-4, "mol/cm^3"),
        "pentane": (3.0e-3, "mol/cm^3"),
    },
    initialSurfaceCoverages={
        "vacantX": 1.0,
    },
    surfaceVolumeRatio=(1.0e5, "m^-1"),
    terminationConversion={
        "CH4": 0.99,
    },
    terminationTime=(0.1, "s"),
)

solvation(solvent="pentane")

simulator(
    atol=1e-18,
    rtol=1e-12,
)

model(
    toleranceMoveToCore=0.01,
    toleranceKeepInEdge=0.001,
    toleranceInterruptSimulation=1e8,
    maximumEdgeSpecies=100,
    minCoreSizeForPrune=10,
    minSpeciesExistIterationsForPrune=2,
    maxNumObjsPerIter=3,
    maxNumSpecies=10,
)

options(
    units="si",
    saveEdgeSpecies=True,
)
