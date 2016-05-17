# Data sources
database(
    thermoLibraries=['surfaceThermo', 'primaryThermoLibrary'],
    reactionLibraries = [('Surface', True)],
    seedMechanisms = [],
    kineticsDepositories = ['training'],
    kineticsFamilies = 'default',
    kineticsEstimator = 'rate rules',
)

# List of species
#species(
#    label='methyl',
#    reactive=True,
#    structure=SMILES("[CH3]"),
#)

# List of species
#species(
#    label='methylene',
#    reactive=True,
#    structure=SMILES("[CH2]"),
#)

species(
    label='methane',
    reactive=True,
    structure=SMILES("[CH4]"),
)

#species(
#   label='c2h4',
#   reactive=True,
#   structure=adjacencyList(
#       """
#1 C u0 p0 c0 {2,D} {3,S} {4,S}
#2 C u0 p0 c0 {1,D} {5,S} {6,S}
#3 H u0 p0 c0 {1,S}
#4 H u0 p0 c0 {1,S}
#5 H u0 p0 c0 {2,S}
#6 H u0 p0 c0 {2,S}
#"""),
#)

species(
   label='o2',
   reactive=True,
   structure=adjacencyList(
       """
1 O u1 p2 c0 {2,S}
2 O u1 p2 c0 {1,S}
"""),
)


species(
    label='site',
    reactive=True,
    structure=adjacencyList("1 X u0"),
)

# Reaction systems
surfaceReactor(
    temperature=(1350,'K'),
    initialPressure=(1.0, 'bar'),
    initialGasMoleFractions={
#        "methyl": 1.0,
        "methane": 1.0,
	"o2": 1.0,
    },
    initialSurfaceCoverages={
        "site": 2.0,
    },
    surfaceVolumeRatio=(1.e5, 'm^-1'),
    surfaceSiteDensity=(2.9e-9, 'mol/cm^2'),
    terminationTime=(1e-1, 's'),
)

simulator(
    atol=1e-20,
    rtol=1e-12,
)

model(
    toleranceKeepInEdge=0.0,
    toleranceMoveToCore=1e-6,
    toleranceInterruptSimulation=0.1,
    maximumEdgeSpecies=100000
)

options(
    units='si',
    saveRestartPeriod=None,
    generateOutputHTML=True,
    generatePlots=False,
    saveEdgeSpecies=True,
    saveSimulationProfiles=True,
)
