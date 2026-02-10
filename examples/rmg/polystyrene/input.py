# 1. Database
database(
    thermoLibraries=['primaryThermoLibrary', 'BurkeH2O2', 'thermo_DFT_CCSDTF12_BAC', 'DFT_QCI_thermo', 'CBS_QB3_1dHR'],
    reactionLibraries=['primaryH2O2', 'NOx2018'],
    transportLibraries=['OneDMinN2', 'PrimaryTransportLibrary', 'NOx2018', 'GRI-Mech'],
    seedMechanisms=[],
    kineticsDepositories='default',
    kineticsFamilies='default',
    kineticsEstimator='rate rules',
)

# 2. Species Definitions
species(
    label='N2',
    reactive=False,
    structure=SMILES("N#N")
)

# 3. Polymer Definition
polymer(
    label='PS',
    monomer='[CH2][CH]c1ccccc1',
    end_groups=['[CH3]', '[H]'], 
    cutoff=3,
    Mn=5000.0,
    Mw=6000.0,
    initial_mass=1.0,
)

# 4. Polymer Phase Definition
pp = polymer_phase(
    label='Polystyrene_Melt',
    species=['PS', 'N2'], 
    solvent='N2',
    density=(1050.0, 'kg/m^3'),
)

# 5. Hybrid Polymer Reactor
hybridPolymerReactor(
    temperature=(1000.0, 'K'),
    pressure=(1.0, 'bar'),
    initialMoles={
        'N2': 0.99, 
        'PS': 0.01, 
    },
    polymerPhase=pp,
    terminationTime=(0.1, 's'),
    sensitivity=None,
    constant_gas_volume=False,
)

# 6. Model and Simulator Settings
model(
    toleranceMoveToCore=0.1,
    toleranceInterruptSimulation=0.1,
    filterReactions=False,
    filterThreshold=100000000.0,
    maxNumObjsPerIter=1,
    terminateAtMaxObjects=False,
)

simulator(atol=1e-16, rtol=1e-08, sens_atol=1e-06, sens_rtol=0.0001)

# 7. PDep
#pressureDependence(
#    method='modified strong collision',
#    maximumGrainSize=(2.0, 'kJ/mol'),
#    minimumNumberOfGrains=250,
#    temperatures=(300, 2500, 'K', 10),
#    pressures=(0.1, 100, 'bar', 10),
#    interpolation=('Chebyshev', 6, 4),
#    maximumAtoms=16,
#)

options(
    name='Seed',
    generateSeedEachIteration=True,
    saveSeedToDatabase=False,
    units='si',
    generateOutputHTML=True,
    generatePlots=False,
    saveSimulationProfiles=True,
    verboseComments=False,
    saveEdgeSpecies=True,
    keepIrreversible=False,
    trimolecularProductReversible=True,
    wallTime='00:00:00:00',
    saveSeedModulus=-1,
)

generatedSpeciesConstraints(
    allowed=['input species', 'seed mechanisms', 'reaction libraries'],
    maximumCarbonAtoms=8,
    maximumOxygenAtoms=4,
    maximumNitrogenAtoms=0,
    maximumSiliconAtoms=0,
    maximumSulfurAtoms=0,
    maximumHeavyAtoms=5,
    maximumRadicalElectrons=2,
    maximumSingletCarbenes=1,
    maximumCarbeneRadicals=0,
    allowSingletO2=True,
)
