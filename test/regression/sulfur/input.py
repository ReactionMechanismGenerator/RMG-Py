# Tests sulfur chemistry and the following functionalities:
# Filter reactions, Adding multiple species, Multiple reactors, Species constraints, Libraries and Seed

database(
    thermoLibraries=['BurkeH2O2','thermo_DFT_CCSDTF12_BAC','primaryNS','primaryThermoLibrary'],
    reactionLibraries=['primarySulfurLibrary'],
    seedMechanisms=['BurkeH2O2inN2'],
    kineticsDepositories=['training'],
    kineticsFamilies='default',
    kineticsEstimator='rate rules',
)

species(
    label='H2S',
    reactive=True,
    structure=SMILES('S'),
)

species(
    label='O2',
    reactive=True,
    structure=SMILES('[O][O]'),
)

species(
    label='N2',
    reactive=False,
    structure=SMILES('N#N'),
)

simpleReactor(
    temperature=(500,'K'),
    pressure=(30,'bar'),
    initialMoleFractions={
        "H2S": 0.000756,
        "O2": 0.001290,
        "N2": 0.997954,
    },
    terminationTime=(0.1,'s'),
)
simpleReactor(
    temperature=(900,'K'),
    pressure=(30,'bar'),
    initialMoleFractions={
        "H2S": 0.000756,
        "O2": 0.001290,
        "N2": 0.997954,
    },
    terminationTime=(0.1,'s'),
)

simulator(
    atol=1e-16,
    rtol=1e-8,
)

model(
    toleranceKeepInEdge=0,
    toleranceMoveToCore=0.1,
    toleranceInterruptSimulation=0.1,
    maximumEdgeSpecies=300000,
    filterReactions=True,
)

options(
    generateOutputHTML=False,
    generatePlots=False,
    saveEdgeSpecies=True,
    saveSimulationProfiles=False,
    saveRestartPeriod=None,
)

generatedSpeciesConstraints(
    allowed=['input species','seed mechanisms','reaction libraries'],
    maximumCarbonAtoms=0,
    maximumOxygenAtoms=4,
    maximumNitrogenAtoms=2,
    maximumSiliconAtoms=0,
    maximumSulfurAtoms=3,
    maximumHeavyAtoms=5,
    maximumRadicalElectrons=1,
    allowSingletO2=False,
)
