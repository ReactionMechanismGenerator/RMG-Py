# Tests nitrogen chemistry and the following functionalities:
# Filter reactions, Adding multiple species, Thermo pruning, Species constraints,
# Libraries and Seed, Rate ratio termination

database(
    thermoLibraries=['BurkeH2O2','primaryThermoLibrary','thermo_DFT_CCSDTF12_BAC','DFT_QCI_thermo','FFCM1(-)','NitrogenCurran'],
    reactionLibraries=[('FFCM1(-)',False),'Nitrogen_Dean_and_Bozzelli'],
    seedMechanisms=['BurkeH2O2inArHe'],
    kineticsDepositories=['training'],
    kineticsFamilies='default',
    kineticsEstimator='rate rules',
)

species(
    label='NC',
    reactive=True,
    structure=SMILES("NC"),
)

species(
    label='O2',
    reactive=True,
    structure=SMILES("[O][O]"),
)

species(
    label='Ar',
    reactive=False,
    structure=adjacencyList('1 Ar u0 p4 c0'),
)

simpleReactor(
    temperature=(1399,'K'),
    pressure=(1.93,'atm'),
    initialMoleFractions={'NC': 0.0005,
                          'O2': 0.002,
                          'Ar': 0.9975},
    terminationRateRatio=0.01,
)

simulator(
    atol=1e-16,
    rtol=1e-8,
)

model(
    toleranceKeepInEdge=0,
    toleranceMoveToCore=0.1,
    toleranceInterruptSimulation=0.1,
    toleranceThermoKeepSpeciesInEdge=0.5,
    maximumEdgeSpecies=10000,
    maxNumObjsPerIter=3,
    terminateAtMaxObjects=True,
    filterReactions=True,
)

options(
    units='si',
    saveEdgeSpecies=True,
)

generatedSpeciesConstraints(
    allowed=['input species','seed mechanisms','reaction libraries'],
    maximumCarbonAtoms=2,
    maximumOxygenAtoms=2,
    maximumNitrogenAtoms=2,
    maximumSiliconAtoms=0,
    maximumSulfurAtoms=0,
    maximumHeavyAtoms=3,
    maximumRadicalElectrons=1,
    allowSingletO2=False,
)
