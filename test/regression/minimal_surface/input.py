# Data sources
database(
    thermoLibraries=['surfaceThermoPt111', 'primaryThermoLibrary', 'thermo_DFT_CCSDTF12_BAC','DFT_QCI_thermo'], # 'surfaceThermoPt' is the default. Thermo data is derived using bindingEnergies for other metals 
    reactionLibraries = [('Surface/CPOX_Pt/Deutschmann2006_adjusted', False)], # when Ni is used change the library to Surface/Deutschmann_Ni 
    seedMechanisms = [],
    kineticsDepositories = ['training'],
    kineticsFamilies = ['surface', 'default'],
    kineticsEstimator = 'rate rules',
)
catalystProperties(
    bindingEnergies = {
        'H': (-2.75368,'eV/molecule'),
        'C': (-7.02516,'eV/molecule'),
        'N': (-4.63225,'eV/molecule'),
        'O': (-3.81153,'eV/molecule'),
    },
    surfaceSiteDensity=(2.483e-9, 'mol/cm^2'),
    coverageDependence=False
)
species(
    label='CH4',
    reactive=True,
    structure=SMILES("C"),
)
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
    label='N2',
    reactive=False,
    structure=SMILES("N#N"),
)
species(
    label='X',
    reactive=True,
    structure=adjacencyList("1 X u0"),
)
# added for training
species(
    label='CO2X',
    reactive=True,
    structure=adjacencyList(
    """
    1 O u0 p2 c0 {3,D}
    2 O u0 p2 c0 {3,D}
    3 C u0 p0 c0 {1,D} {2,D}
    4 X u0 p0 c0
    """),
)
species(
    label='COX',
    reactive=True,
    structure=adjacencyList(
    """
    1 O u0 p2 c0 {2,D}
    2 C u0 p0 c0 {1,D} {3,D}
    3 X u0 p0 c0 {2,D}
    """),
)
species(
    label='OX',
    reactive=True,
    structure=adjacencyList(
    """
    1 O u0 p2 c0 {2,D}
    2 X u0 p0 c0 {1,D}
    """),
)
# If you would like to forbid the bidentate form of absorbed CO2 from your model,
# use the following `CO2_bidentate` forbidden structure
# forbidden(
#     label='CO2_bidentate',
#     structure=SMILES("O=C(*)O*"),
# )
#----------
# Reaction systems
surfaceReactor(
    temperature=(1000, 'K'),
    initialPressure=(1.0, 'bar'),
    initialGasMoleFractions={
        "CH4": 0.15,
        "O2": 0.15,
        "N2": 0.7,
    },
    initialSurfaceCoverages={
        "X": 1.0,
    },
    surfaceVolumeRatio=(1.e5, 'm^-1'),
    terminationConversion={"CH4": 0.95,},
    terminationTime=(0.1, 's'),
    terminationRateRatio=0.01,
)
simulator(
    atol=1e-18,
    rtol=1e-12,
)
model(
    toleranceKeepInEdge=0.0,
    toleranceMoveToCore=1e-1,
    toleranceInterruptSimulation=0.1,
    maximumEdgeSpecies=100000,
)
options(
    units='si',
    generateOutputHTML=False,
    generatePlots=False, # Enable to make plots of core and edge size etc. But takes a lot of the total runtime!
    saveEdgeSpecies=True,
    saveSimulationProfiles=False,
)
generatedSpeciesConstraints(
    allowed=['input species','reaction libraries'],
    maximumCarbonAtoms=2, 
    maximumOxygenAtoms=2,
    maximumSurfaceSites=2,
)


