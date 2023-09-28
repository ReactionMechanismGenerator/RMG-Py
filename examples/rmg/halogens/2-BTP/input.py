# Example input file for 2-BTP/CF3Br in methane flames

database(
    thermoLibraries = ['primaryThermoLibrary', 'DFT_QCI_thermo', 'thermo_DFT_CCSDTF12_BAC',
                    'halogens', 'CHOBr_G4','CHOF_G4','CHOFBr_G4', '2-BTP', '2-BTP_G4', 'Fluorine',
                    'CBS_QB3_1dHR','NISTThermoLibrary','JetSurF2.0'],
    #reactionLibraries = ['2-BTP/seed'], #optional: use NISTs 2-BTP kinetics library as seed or reaction library
    reactionLibraries = [],
    seedMechanisms = ['FFCM1(-)'], #optional: use FFCM1(-) kinetics library as seed for methane combustion
    kineticsDepositories = ['training'],
    kineticsFamilies = ['default','halogens'], # load both the `default` and `halogens` recommended kinetics families
    kineticsEstimator = 'rate rules',
)

species(
    label = "2-BTP",
    reactive = True,
    structure = SMILES("FC(F)(F)C(Br)=C")
)

species(
    label = "CF3Br",
    reactive = True,
    structure = SMILES("FC(F)(F)Br")
)

species(
    label = 'CH4',
    reactive = True,
    structure = SMILES('C')
)

species(
    label = "O2",
    reactive = True,
    structure = adjacencyList(
    """
    multiplicity 3
    1 O u1 p2 c0 {2,S}
    2 O u1 p2 c0 {1,S}
""" 
))

species(
    label = "N2",
    reactive = False,
    structure = adjacencyList(
    """
    1 N u0 p1 c0 {2,T}
    2 N u0 p1 c0 {1,T}
    """
))

# reactor for 2-BTP in methane flame
simpleReactor(
        temperature=[(1000,'K'),(2200,'K')],
        pressure=(1.0,'bar'),
        nSims=24,
        initialMoleFractions={
        "2-BTP": [0.0,0.04],
        "CH4" : [0.05,0.15],
        "O2": 0.21,
        "N2": 0.78,
        },
        terminationConversion={
        'CH4': 0.99,
        },
        terminationTime=(1,'s'),
        )

# pressue dependence is not yet implemented for halogen chemistry
# pressureDependence(
#     method='modified strong collision',
#     maximumGrainSize=(0.5,'kcal/mol'),
#     minimumNumberOfGrains=250,
#     temperatures=(300,2200,'K',8),
#     pressures=(0.01,100,'bar',5),
#     interpolation=('Chebyshev', 6, 4),
# )


simulator(
    atol = 1e-16,
    rtol = 1e-08,
    sens_atol = 1e-06,
    sens_rtol = 0.0001,
)

model(
    toleranceMoveToCore = 0.25,
    toleranceInterruptSimulation = 1,
    maximumEdgeSpecies = 5e5,
    filterReactions = True,
    filterThreshold = 5e8,
    minCoreSizeForPrune=400,
    minSpeciesExistIterationsForPrune=4,
)


# model(
#     toleranceMoveToCore=1,
#     toleranceInterruptSimulation=1e8,
#     toleranceKeepInEdge=0.05,
#     maximumEdgeSpecies=200000,
#     minCoreSizeForPrune=50,
#     minSpeciesExistIterationsForPrune=2,
#     filterReactions = True,
#     filterThreshold=5e8
#     )

generatedSpeciesConstraints(
    allowed=['input species','seed mechanisms','reaction libraries'],
    maximumCarbonAtoms=12,
    maximumOxygenAtoms=6,
    # maximumHeavyAtoms=24,
    maximumRadicalElectrons=2,
    maximumSingletCarbenes=1,
    maximumCarbeneRadicals=0,
    allowSingletO2 = False,
)

options(
    units = "si",
    saveRestartPeriod = None,
    generateOutputHTML = True,
    generatePlots = False,
    saveSimulationProfiles = True,
    saveEdgeSpecies = False,
    keepIrreversible = True,
    verboseComments = False,
)
