database(
    thermoLibraries = ['primaryThermoLibrary', 'DFT_QCI_thermo'],
    reactionLibraries = [],
    seedMechanisms = [],
    kineticsDepositories = 'default',
    kineticsFamilies = ['H_Abstraction','R_Addition_MultipleBond','R_Recombination',
    'Disproportionation', 'Intra_R_Add_Exocyclic', 'Intra_R_Add_Endocyclic'],
    kineticsEstimator = 'rate rules',
)

generatedSpeciesConstraints(
    allowed=['input species','seed mechanisms','reaction libraries'],
    maximumRadicalElectrons = 1,
    maximumCarbonAtoms = 10,
)

species(
    label = "benzene",
    structure = SMILES("c1ccccc1")
)

species(
    label = 'ethyne',
    structure = SMILES("C#C")
)

simpleReactor(
    temperature = (1500,"K"),
    pressure = (1,"bar"),
    initialMoleFractions={
        "benzene": 0.1,
        "ethyne": 0.9,
    },
    terminationTime = (100,"s"),
    terminationConversion = {
        "benzene": 0.5,
    },
)

simulator(
    atol = 1e-16,
    rtol = 1e-08,
)

model(
    toleranceMoveToCore = 0.1,
    toleranceInterruptSimulation = 0.1,
    filterReactions = True,
    maxNumSpecies = 15,
)

options(
    units = "si",
    saveRestartPeriod = None,
    generateOutputHTML = False,
    generatePlots = False,
    saveSimulationProfiles = False,
    saveEdgeSpecies = True,
    verboseComments = False,
)

