database(
    thermoLibraries = ['primaryThermoLibrary'],
    reactionLibraries = [],
    seedMechanisms = ['GRI-Mech3.0'],
    kineticsDepositories = ['training'],
    kineticsFamilies = ['R_Recombination'],
    kineticsEstimator = 'rate rules',
)

species(
    label='ethane',
    reactive=True,
    structure=SMILES("CC"),
)

species(
        label='N2',
        reactive=False,
        structure=SMILES("N#N"))

simpleReactor(
    temperature=(1350,'K'),
    pressure=[(1.0,'bar'),(10.0,'bar')],
    initialMoleFractions={
        "ethane": [0.5,1.0],
        "N2":1.0,
    },
    terminationConversion={
        'ethane': 0.000000000001,
    },
    terminationTime=(1e6,'s'),
    terminationRateRatio=0.01,
    balanceSpecies='N2',
)

liquidReactor(
    temperature=(1350,'K'),
    initialConcentrations={
        "ethane": [(0.5,'mol/m^3'),(1.0,'mol/m^3')],
        "N2":(1.0,'mol/m^3'),
    },
    terminationConversion={
        'ethane': 0.000000000001,
    },
    terminationTime=(1e6,'s'),
    terminationRateRatio=0.01,
)

liquidReactor(    
    temperature=[(1000,'K'),(1200,'K')],
    initialConcentrations={
        "ethane": (1.0,'mol/m^3'),
        "N2":(1.0,'mol/m^3'),
    },
    terminationConversion={
    	'ethane': 0.000000000001,
    },
    terminationTime=(1e6,'s'),
    terminationRateRatio=0.01,
)

simulator(
    atol=1e-16,
    rtol=1e-8
)

model(
    toleranceKeepInEdge=0.0,
    toleranceMoveToCore=0.2,
    toleranceInterruptSimulation=0.2,
)
options(
    name='testSeed',
    units='si',
    generateSeedEachIteration=True,
    saveSeedToDatabase=True,
    saveRestartPeriod=None,
    generateOutputHTML=False,
    generatePlots=False,
    saveEdgeSpecies=False,
    saveSimulationProfiles=False,
)

generatedSpeciesConstraints(allowed=['seed mechanisms','reaction libraries'],
maximumRadicalElectrons=3,maximumCarbeneRadicals=3,maximumSingletCarbenes=3)
