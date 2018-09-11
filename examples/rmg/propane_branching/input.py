database(
    thermoLibraries = ['BurkeH2O2','primaryThermoLibrary','thermo_DFT_CCSDTF12_BAC','CBS_QB3_1dHR','DFT_QCI_thermo'],
    reactionLibraries = [],
    seedMechanisms = ['BurkeH2O2inN2','ERC-FoundationFuelv0.9'],
    kineticsDepositories = 'default', 
    kineticsFamilies = 'default',
    kineticsEstimator = 'rate rules',
)


species(
    label='propane',
    reactive=True,	
    structure=SMILES("CCC"),
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


simpleReactor(
    temperature=(700,'K'),
    pressure=(10.0,'bar'),
    initialMoleFractions={
        "N2": 4,
        "O2": 1,
        "propane": 1./5.0,
    },

    terminationConversion={
        'propane': .90,
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
    filterReactions=True,
    toleranceBranchReactionToCore=0.001,
    branchingIndex=0.5,
    branchingRatioMax=1.0,
)

options(
    units='si',
    generateSeedEachIteration=False,
    generateOutputHTML=False,
    generatePlots=False,
    saveSimulationProfiles=False,
    verboseComments=False,
    saveEdgeSpecies=False,
    keepIrreversible=False,
)

generatedSpeciesConstraints(
    allowed=['input species','seed mechanisms','reaction libraries'],
    maximumCarbonAtoms=5,
    maximumOxygenAtoms=8,
    maximumNitrogenAtoms=0,
    maximumSiliconAtoms=0,
    maximumSulfurAtoms=0,
    maximumRadicalElectrons=2,
)
