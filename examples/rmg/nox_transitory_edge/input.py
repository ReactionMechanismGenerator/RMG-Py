database(
    thermoLibraries = ['BurkeH2O2','FFCM1(-)','primaryThermoLibrary','thermo_DFT_CCSDTF12_BAC','CBS_QB3_1dHR','DFT_QCI_thermo','primaryNS', 'NitrogenCurran'],
    reactionLibraries = ['primaryNitrogenLibrary'],
    seedMechanisms = ['BurkeH2O2inN2','FFCM1(-)'],
    kineticsDepositories = 'default',
    kineticsFamilies = 'default',
    kineticsEstimator = 'rate rules',
)

species(
    label='ethane',
    reactive=True,	
    structure=SMILES("CC"),
)
species(
    label='O2',
    structure=SMILES("[O][O]"),
)
species(
    label='N2',
    reactive=False,
    structure=adjacencyList("""
    1 N u0 p1 c0 {2,T}
	2 N u0 p1 c0 {1,T}
	"""),
)

species(
   label='NO',
   reactive=True,
   structure=SMILES("[N]=O"),
)

constantTPIdealGasReactor(
    temperature=(900,'K'),
    pressure=(10.0,'bar'),
    initialMoleFractions={
        "N2": 4,
        "O2": 1,
        "ethane": 2.0/7.0,
    },

    terminationConversion={
        'ethane': .90,
    },
   terminationTime=(40,'s'),
)

simulator(
    atol=1e-16,
    rtol=1e-8,
)

model(
    toleranceMoveToCore=0.3,
    toleranceInterruptSimulation=0.3,
    maxNumObjsPerIter=1,
    terminateAtMaxObjects=True,
    filterReactions=False,
    #toleranceBranchReactionToCore=0.001,
    #branchingIndex=0.5,
    #branchingRatioMax=1.0,
    toleranceTransitoryDict={'NO':0.2},
)

options(
    units='si',
    generateSeedEachIteration=False,
    generateOutputHTML=False,
    generatePlots=False,
    saveSimulationProfiles=False,
    verboseComments=False,
    saveEdgeSpecies=True,
    keepIrreversible=False,
)

generatedSpeciesConstraints(
    allowed=['input species','seed mechanisms','reaction libraries'],
    #maximumCarbonAtoms=5,
    #maximumOxygenAtoms=8,
    #maximumNitrogenAtoms=0,
    maximumSiliconAtoms=0,
    maximumSulfurAtoms=0,
    maximumRadicalElectrons=2,
)
