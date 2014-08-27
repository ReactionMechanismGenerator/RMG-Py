# Data sources
database(
    thermoLibraries = ['KlippensteinH2O2', 'primaryThermoLibrary','DFT_QCI_thermo','CH','CHN','CHO','CHON','CN','NISTThermoLibrary','thermo_DFT_CCSDTF12_BAC','GRI-Mech3.0-N'],
    reactionLibraries = [('Nitrogen_Dean_and_Bozzelli',False)], 
    seedMechanisms = ['ERC-FoundationFuelv0.9'],
    kineticsDepositories = ['training'],
    kineticsFamilies = 'default',
    kineticsEstimator = 'rate rules',
)

# Constraints on generated species
generatedSpeciesConstraints(
    allowed = ['seed mechanisms', 'reaction libraries'],
    #maximumCarbonAtoms = 7,
    #maximumHydrogenAtoms = 8,
    #maximumOxygenAtoms = 5,
    maximumNitrogenAtoms = 2,
    #maximumSiliconAtoms = 0,
    #maximumSulfurAtoms = 0,
    #maximumHeavyAtoms = 3,
    maximumRadicalElectrons = 2,
)

# List of species
species(
    label='CH3NO2',
    reactive=True,
        structure=adjacencyList(
        """
        1 C u0 p0     {2,S} {3,S} {4,S} {5,S}
        2 H u0 p0     {1,S}
        3 H u0 p0     {1,S}
        4 H u0 p0     {1,S}
        5 N u0 p0 c+1 {1,S} {6,D} {7,S}
        6 O u0 p2     {5,D}
        7 O u0 p3 c-1 {5,S}
        """),
)

species(
    label='O2',
    reactive=True,
        structure=adjacencyList(
        """
        1 O u1 p2 {2,S}
        2 O u1 p2 {1,S}
        """),
)

species(
    label='N2',
    reactive=True,
        structure=adjacencyList(
        """
        1 N u0 p1 {2,T}
        2 N u0 p1 {1,T}
        """),
)

# Reaction systems

simpleReactor(
    temperature=(1500,'K'),
    pressure=(10.0,'bar'),
    initialMoleFractions={
        "CH3NO2": 0.1,
        "O2": 0.21,
        "N2": 0.69,
    },
    terminationConversion={
        'CH3NO2': 0.1,
    },
)

model(
    toleranceKeepInEdge=1e-5,
    toleranceMoveToCore=0.1,
    toleranceInterruptSimulation=0.1,
    maximumEdgeSpecies=10000
)

options(
    units='si',
    saveRestartPeriod=None,
    drawMolecules=False,
    generatePlots=False,
)
