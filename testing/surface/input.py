# Data sources
database(
    thermoLibraries=['surfaceThermo', 'primaryThermoLibrary'],
    reactionLibraries = [('Deutschmann_Ni', True)],
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



species(
    label='CH4',
    reactive=True,
    structure=SMILES("[CH4]"),
)

#species(
#    label='water',
#    reactive=True,
#    structure=adjacencyList(
#       """
#1 O u0 p2 {2,S} {3,S} {4,vdW}
#2 H u0 p0 {1,S}
#3 H u0 p0 {1,S}
#4 X u0 p0 {1,vdW}
#"""),
#)

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
   label='O2',
   reactive=True,
   structure=adjacencyList(
       """
1 O u1 p2 c0 {2,S}
2 O u1 p2 c0 {1,S}
"""),
)

species(
    label='CO2',
    reactive=True,
    structure=SMILES("O=C=O"),
)

species(
    label='H2O',
    reactive=True,
    structure=SMILES("O"),
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
        "CH4": 1.0,
        "O2": 0.0,
        "CO2": 1.0,
        "H2O": 0.0,
    },
    initialSurfaceCoverages={
        "site": 2.0,
    },
    surfaceVolumeRatio=(1.e5, 'm^-1'),
    surfaceSiteDensity=(2.9e-9, 'mol/cm^2'),
    terminationConversion = { "CH4":0.5,},
    terminationTime=(1e-2, 's'),
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
