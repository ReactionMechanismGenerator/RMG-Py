"""
This example is for demonstrating the use of RMG for molecular growth
We use mass weighted fluxes to bias the algorithm toward adding larger species
"""
database(
    thermoLibraries = ['Klippenstein_Glarborg2016','BurkeH2O2','primaryThermoLibrary','thermo_DFT_CCSDTF12_BAC','CBS_QB3_1dHR','DFT_QCI_thermo'],
    reactionLibraries = ['Klippenstein_Glarborg2016','BurkeH2O2inN2'],
    seedMechanisms = [],
    kineticsDepositories = 'default', 
    kineticsFamilies = 'default',
    kineticsEstimator = 'rate rules',
)


species(
    label='butadiene',
    reactive=True,	
    structure=SMILES("C=CC=C"),
)


simpleReactor(
    temperature=(700,'K'),
    pressure=(10.0,'bar'),
    initialMoleFractions={
        "butadiene": 1,
    },

    terminationTime=(1e5,'s'),
    terminationAvgMW=150.0, #terminate simulation when the average molecule weight reaches terminationAvgMW amu
)

simulator(
    atol=1e-16,
    rtol=1e-8,
)

model(
    toleranceMoveToCore=1.0,
    toleranceInterruptSimulation=1.0,
    maxNumObjsPerIter=1,
    fluxBasis='mass', #use mass weighting
    massIndex=2.0,    #weight fluxes by mass to the power of massIndex
    terminateAtMaxObjects=True,
    filterReactions=True,
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
    maximumRadicalElectrons=2,
    maximumHeavyAtoms=15,
)
