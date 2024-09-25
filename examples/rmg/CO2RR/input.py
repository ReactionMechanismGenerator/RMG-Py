# Data sources
database(
    thermoLibraries=['surfaceThermoPt111', 'primaryThermoLibrary', 'thermo_DFT_CCSDTF12_BAC','DFT_QCI_thermo', 'electrocatThermo', 
    # 'CO2RR_Adsorbates_Ag111'
    ],
    reactionLibraries = [('Surface/CPOX_Pt/Deutschmann2006_adjusted', False)],
    seedMechanisms = [],
    kineticsDepositories = ['training'],
    kineticsFamilies = ['electrochem', 
                        # 'surface',
                        'Surface_Abstraction',
                        'Surface_Abstraction_vdW',
                        'Surface_Abstraction_Single_vdW',
                        'Surface_Abstraction_Beta_double_vdW',
                        'Surface_Adsorption_Dissociative',
                        'Surface_Adsorption_Dissociative_Double',
                        'Surface_Adsorption_vdW',
                        'Surface_Dissociation',
                        'Surface_Dissociation_Double_vdW',
                        'Surface_Dissociation_vdW',
                        'Surface_EleyRideal_Addition_Multiple_Bond',
                        'Surface_Migration',
                        ],
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

species(
    label='CO2HX',
    reactive=True,
    structure=adjacencyList("""
1 O u0 p2 c0 {2,S} {4,S}
2 C u0 p0 c0 {1,S} {3,D} {5,S}
3 O u0 p2 c0 {2,D}
4 H u0 p0 c0 {1,S}
5 X u0 p0 c0 {2,S}

"""),
)

species(
    label='OCX',
    reactive=True,
    structure=adjacencyList("""
1 O u0 p2 c0 {2,D} 
2 C u0 p0 c0 {1,D} {3,D}
3 X u0 p0 c0 {2,D}
"""),
)

species(
    label='OX',
    reactive=True,
    structure=adjacencyList("""
1 O u0 p2 c0 {2,D} 
2 X u0 p0 c0 {1,D}
"""),
)

species(
    label='CH2O2X',
    reactive=True,
    structure=adjacencyList("""
1 O u0 p2 c0 {3,S} {5,S}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,S} {2,D} {4,S}
4 H u0 p0 c0 {3,S}
5 H u0 p0 c0 {1,S}
6 X u0 p0 c0
"""),
)

species(
    label='CHOX',
    reactive=True,
    structure=adjacencyList("""
1 O u0 p2 c0 {2,D}
2 C u0 p0 c0 {1,D} {3,S} {4,S}
3 H u0 p0 c0 {2,S}
4 X u0 p0 c0 {2,S}
"""),
)

species(
    label='CH2OX',
    reactive=True,
    structure=adjacencyList("""
1 O u0 p2 c0 {2,D}
2 C u0 p0 c0 {1,D} {3,S} {4,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {2,S}
5 X u0 p0 c0
"""),
)


forbidden(
        label='CO2-bidentate',
        structure=adjacencyList(
                """
                1 O u0 p2 c0 {2,D}
                2 C u0 p0 c0 {1,D} {3,S} {4,S}
                3 X u0 p0 c0 {2,S}
                4 O u0 p2 c0 {2,S} {5,S}
                5 X u0 p0 c0 {4,S}
                """
        )
)

liquidSurfaceReactor(
    temperature=(300,'K'),
    liqPotential=(0,'V'),
    surfPotential=(-2.0,'V'),
    initialConcentrations={
        "CO2": (1e-3,'mol/cm^3'),
        "proton": (1e-4,'mol/m^3'),
    },
	initialSurfaceCoverages={
        # "HX": 0.5,
        # # "CXO2": 0.0,
        "CHO2X": 0.1,
        "CO2HX": 0.1,
        "vacantX": 0.1,
        "CO2X": 0.4,
        'OX': 0.1,
        'OCX': 0.1,
        'CH2O2X': 0.05,
        'CHOX': 0.04,
        'CH2OX': 0.01
    },
    surfaceVolumeRatio=(1.0e5, 'm^-1'),
    terminationTime=(1.0e3,'sec'),
    # terminationConversion={'CO2': 0.90},
    # constantSpecies=["proton"],
 )

liquidSurfaceReactor(
    temperature=(300,'K'),
    liqPotential=(0,'V'),
    surfPotential=(-1.5,'V'),
    initialConcentrations={
        "CO2": (1e-3,'mol/cm^3'),
        "proton": (1e-4,'mol/m^3'),
    },
	initialSurfaceCoverages={
        # "HX": 0.5,
        # # "CXO2": 0.0,
        "CHO2X": 0.1,
        "CO2HX": 0.1,
        "vacantX": 0.1,
        "CO2X": 0.4,
        'OX': 0.1,
        'OCX': 0.1,
        'CH2O2X': 0.05,
        'CHOX': 0.04,
        'CH2OX': 0.01
    },
    surfaceVolumeRatio=(1.0e5, 'm^-1'),
    terminationTime=(1.0e3,'sec'),
    # terminationConversion={'CO2': 0.90},
    # constantSpecies=["proton"],
 )

liquidSurfaceReactor(
    temperature=(300,'K'),
    liqPotential=(0,'V'),
    surfPotential=(-1.0,'V'),
    initialConcentrations={
        "CO2": (1e-3,'mol/cm^3'),
        "proton": (1e-4,'mol/m^3'),
    },
	initialSurfaceCoverages={
        # "HX": 0.5,
        # # "CXO2": 0.0,
        "CHO2X": 0.1,
        "CO2HX": 0.1,
        "vacantX": 0.1,
        "CO2X": 0.4,
        'OX': 0.1,
        'OCX': 0.1,
        'CH2O2X': 0.05,
        'CHOX': 0.04,
        'CH2OX': 0.01
    },
    surfaceVolumeRatio=(1.0e5, 'm^-1'),
    terminationTime=(1.0e3,'sec'),
    # terminationConversion={'CO2': 0.90},
    # constantSpecies=["proton"],
 )

# liquidSurfaceReactor(
#     temperature=(300,'K'),
#     liqPotential=(0,'V'),
#     surfPotential=(-0.5,'V'),
#     initialConcentrations={
#         "CO2": (1e-3,'mol/cm^3'),
#         "proton": (1e-4,'mol/m^3'),
#     },
# 	initialSurfaceCoverages={
#         # "HX": 0.5,
#         # # "CXO2": 0.0,
#         "CHO2X": 0.1,
#         "CO2HX": 0.1,
#         "vacantX": 0.1,
#         "CO2X": 0.4,
#         'OX': 0.1,
#         'OCX': 0.1,
#         'CH2O2X': 0.05,
#         'CHOX': 0.04,
#         'CH2OX': 0.01
#     },
#     surfaceVolumeRatio=(1.0e5, 'm^-1'),
#     terminationTime=(1.0e3,'sec'),
#     # terminationConversion={'CO2': 0.90},
#     # constantSpecies=["proton"],
#  )

solvation(
	solvent='water'
)

simulator(
    atol=1e-16,
    rtol=1e-8,
)

model(
    toleranceKeepInEdge=1E-16,
    toleranceMoveToCore=1E-3,
    toleranceRadMoveToCore=1E-6,
    toleranceInterruptSimulation=1E6,
    filterReactions=False,
    maximumEdgeSpecies=5000,
    toleranceBranchReactionToCore=1E-6,
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
    maximumCarbonAtoms=3,
    maximumOxygenAtoms=2,
    maximumRadicalElectrons=1,
)
