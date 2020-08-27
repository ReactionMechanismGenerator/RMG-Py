# Data sources
database(
    thermoLibraries=['surfaceThermoPt', 'primaryThermoLibrary', 'thermo_DFT_CCSDTF12_BAC'], 
    reactionLibraries = [('Surface/Deutschmann_Ni', True)], # when Pt is used change the library to Surface/CPOX_Pt/Deutschmann2006
    seedMechanisms = [],
    kineticsDepositories = ['training'],
    kineticsFamilies = ['surface','default'],
    kineticsEstimator = 'rate rules',
)

catalystProperties(
    bindingEnergies = {  # values for Ni(111)
                        'H': (-2.892, 'eV/molecule'),
                        'O': (-4.989, 'eV/molecule'),
                        'C': (-6.798, 'eV/molecule'),
                        'N': (-5.164, 'eV/molecule'), 
                      },
    surfaceSiteDensity=(3.148e-9, 'mol/cm^2'), # values for Ni(111)
)

# List of species

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
    label='H2',
    reactive=True,
    structure=SMILES("[H][H]"),
)

species(
    label='CO',
    reactive=True,
    structure=SMILES("[C-]#[O+]"),
)

species(
    label='C2H6',
    reactive=True,
    structure=SMILES("CC"),
)

species(
    label='CH2O',
    reactive=True,
    structure=SMILES("C=O"),
)

species(
    label='CH3',
    reactive=True,
    structure=SMILES("[CH3]"),
)

species(
    label='C3H8',
    reactive=True,
    structure=SMILES("CCC"),
)

species(
    label='H',
    reactive=True,
    structure=SMILES("[H]"),
)

species(
    label='C2H5',
    reactive=True,
    structure=SMILES("C[CH2]"),
)

species(
    label='CH3OH',
    reactive=True,
    structure=SMILES("CO"),
)

species(
    label='HCO',
    reactive=True,
    structure=SMILES("[CH]=O"),
)

species(
    label='CH3CHO',
    reactive=True,
    structure=SMILES("CC=O"),
)

species(
    label='OH',
    reactive=True,
    structure=SMILES("[OH]"),
)

species(
    label='C2H4',
    reactive=True,
    structure=SMILES("C=C"),
)

#-------
species(
    label='site',
    reactive=True,
    structure=adjacencyList("1 X u0"),
)
#----------
# Reaction systems
surfaceReactor(
    temperature=(1000,'K'),
    initialPressure=(1.0, 'bar'),
    initialGasMoleFractions={
        "CH4": 1.0,
        "O2": 0.0,
        "CO2": 1.2,
        "H2O": 1.2,
        "H2": 0.0,
        "CH3OH": 0.0,
        "C2H4": 0.0,
    },
    initialSurfaceCoverages={
        "site": 1.0,
    },
    surfaceVolumeRatio=(1.e5, 'm^-1'),
    terminationConversion = { "CH4":0.9,},
    terminationTime=(0.01, 's'),
)

simulator(
    atol=1e-18,
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
    generateOutputHTML=True,
    generatePlots=False, # Enable to make plots of core and edge size etc.. But takes a lot of the total runtime!
    saveEdgeSpecies=True,
    saveSimulationProfiles=True,
    verboseComments=True,
)
