# Data sources
database(
    thermoLibraries=['surfaceThermoPt111', 'primaryThermoLibrary', 'thermo_DFT_CCSDTF12_BAC','DFT_QCI_thermo', 'electrocatThermo'],
    reactionLibraries = [('Surface/CPOX_Pt/Deutschmann2006_adjusted', False)],
    seedMechanisms = [],
    kineticsDepositories = ['training'],
    kineticsFamilies = ['surface','electrochem'],
    kineticsEstimator = 'rate rules',

)

catalystProperties(
    metal = 'Ag111'
)

# List of species
species(
    label='CO2',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {2,D}
2 C u0 p0 c0 {1,D} {3,D}
3 O u0 p2 c0 {2,D}
"""),
)


species(
   label='proton',
   reactive=True,
   structure=adjacencyList(
       """
1 H u0 p0 c+1
"""),
)

species(
    label='vacantX',
    reactive=True,
    structure=adjacencyList("1 X u0"),
)

species(
   label='H',
   reactive=True,
   structure=adjacencyList(
       """
1 H u1 p0 c0
"""),
)

# species(
#    label='water',
#    reactive=True,
#    structure=adjacencyList(
#        """
# 1 H u0 p0 c0 {3,S}
# 2 H u0 p0 c0 {3,S}
# 3 O u0 p2 c0 {1,S} {2,S}
# """),
# )

# species(
#     label='CO2X2',
#     reactive=True,
#     structure=adjacencyList("""
# 1 O u0 p2 c0 {3,S} {5,S}
# 2 O u0 p2 c0 {3,D}
# 3 C u0 p0 c0 {1,S} {2,D} {4,S}
# 4 X u0 p0 c0 {3,S}
# 5 X u0 p0 c0 {1,S}
# """),
# )

species(
    label='CO2X',
    reactive=True,
    structure=adjacencyList("""
1 O u0 p2 c0 {3,D}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,D} {2,D}
4 X u0 p0 c0
"""),
)

species(
    label='CHO2X',
    reactive=True,
    structure=adjacencyList("""
1 O u0 p2 c0 {3,S} {5,S}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,S} {2,D} {4,S}
4 H u0 p0 c0 {3,S}
5 X u0 p0 c0 {1,S}
"""),
)

# species(
#     label='HX',
#     reactive=True,
#     structure=adjacencyList(
#         """
# 1 X u0 {2,S}
# 2 H u0 p0 c0 {1,S}
# """),
# )

liquidSurfaceReactor(
    temperature=(300,'K'),
    liqPotential=(0,'V'),
    surfPotential=(-1.2,'V'),
    initialConcentrations={
        # "CO2": (1e-3,'mol/cm^3'),
        "proton": (1e-3,'mol/m^3'),
        # "water": (0.055, 'mol/cm^3'),
    },
	initialSurfaceCoverages={
        # "HX": 0.5,
        # # "CXO2": 0.0,
        "CHO2X": 1.0,
        # "vacantX": 1.0,
    },
    surfaceVolumeRatio=(1.0e5, 'm^-1'),
    terminationTime=(1.0,'sec'),
    # terminationConversion={'CO2': 0.90},
    # constantSpecies=["proton"],
 )

solvation(
	solvent='water'
)

simulator(
    atol=1e-16,
    rtol=1e-8,
)

model(
    toleranceKeepInEdge=1E-16,
    toleranceMoveToCore=1E-8,
    toleranceRadMoveToCore=1E-10,
    toleranceInterruptSimulation=1E10,
    filterReactions=False,
    maximumEdgeSpecies=10000,
    toleranceBranchReactionToCore=1E-8,
    branchingIndex=0.5,
    branchingRatioMax=1.0,
)

options(
    units='si',
    generateOutputHTML=True,
    generatePlots=True,
    saveEdgeSpecies=True,
    saveSimulationProfiles=False,
)

generatedSpeciesConstraints(
    allowed=['input species','reaction libraries'],
    maximumSurfaceSites=2,
    maximumCarbonAtoms=4,
    maximumOxygenAtoms=4,
    maximumRadicalElectrons=1,
)

