#Example for testing oxidation chemistry (propane), adding multiple species at a time, ranged reactors, staging, 
#species constraints, max species termination, and rate ratio termination

database(
    thermoLibraries = ['primaryThermoLibrary','DFT_QCI_thermo'],
    reactionLibraries = [],
    seedMechanisms = ['BurkeH2O2inN2','ERC-FoundationFuelv0.9'],
    kineticsDepositories = ['training'],
    kineticsFamilies = ['H_Abstraction','R_Addition_MultipleBond','intra_H_migration','R_Recombination','Disproportionation'],
    kineticsEstimator = 'rate rules',
)

species(
    label='propane',
    reactive=True,
    structure=SMILES('CCC'),
)

species(
    label='O2',
    reactive=True,
    structure=SMILES('[O][O]')
)

species(
   label='N2',
   reactive=False,
   structure=SMILES('N#N'),
)

simpleReactor(
    temperature=[(700,'K'),(750,'K')],
    pressure=(10.0,'bar'),
    nSims=2,
    initialMoleFractions={
        "propane": 2.0/7.0,
        "O2": 1.0,
        "N2": 4.0,
    },
    terminationRateRatio=0.5,
)

simulator(
    atol=1e-16,
    rtol=1e-8,
)

simulator(
    atol=1e-16,
    rtol=1e-8,
)

model(
    toleranceMoveToCore=0.01,
    toleranceInterruptSimulation=0.01,
    filterReactions=True,
    maxNumObjsPerIter=1,
    maxNumSpecies=10,
)

model(
    toleranceKeepInEdge=0.0,
    toleranceMoveToCore=0.25,
    toleranceInterruptSimulation=0.25,
    filterReactions=True,
    maxNumObjsPerIter=5,
    terminateAtMaxObjects=True,
)

options(
    units='si',
    saveRestartPeriod=None,
    generateOutputHTML=False,
    generatePlots=False,
    saveEdgeSpecies=True,
    saveSimulationProfiles=False,
)

generatedSpeciesConstraints(
    allowed=['input species', 'seed mechanisms', 'reaction libraries'],
    maximumCarbonAtoms=3,
    maximumOxygenAtoms=4,
    maximumNitrogenAtoms=0,
    maximumSiliconAtoms=0,
    maximumSulfurAtoms=0,
    maximumRadicalElectrons=1,
    maximumSingletCarbenes=0,
    maximumCarbeneRadicals=0,
)
