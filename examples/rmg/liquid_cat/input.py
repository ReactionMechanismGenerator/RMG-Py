# Data sources
database(
    thermoLibraries=['surfaceThermoPt111', 'primaryThermoLibrary', 'thermo_DFT_CCSDTF12_BAC','DFT_QCI_thermo'], # 'surfaceThermoPt' is the default. Thermo data is derived using bindingEnergies for other metals
<<<<<<< HEAD
    reactionLibraries = [],
=======
    reactionLibraries = [('Surface/CPOX_Pt/Deutschmann2006', False)], # when Ni is used change the library to Surface/Deutschmann_Ni
>>>>>>> 2899d7045 (add test example for liquidSurfaceReactor)
    seedMechanisms = [],
    kineticsDepositories = ['training'],
    kineticsFamilies = ['surface','default'],
    kineticsEstimator = 'rate rules',

)

catalystProperties(
    metal = 'Pt111'
)

# List of species
species(
    label='ethylene-glycol',
    reactive=True,
    structure=SMILES("OCCO"),
)


species(
   label='O2',
   reactive=True,
   structure=adjacencyList(
       """
1 O u1 p2 c0 {2,S}
2 O u1 p2 c0 {1,S}
"""),
)

species(
    label='vacantX',
    reactive=True,
    structure=adjacencyList("1 X u0"),
)


liquidSurfaceReactor(
    temperature=(450,'K'),
    initialConcentrations={
        "ethylene-glycol": (1.7935e-2,'mol/cm^3'),
        "O2": (4.953e-6,'mol/cm^3'),
    },
	initialSurfaceCoverages={
        "vacantX": 1.0,
    },
    surfaceVolumeRatio=(1.0e5, 'm^-1'),
    terminationTime=(10.0,'sec'),
    terminationConversion={'ethylene-glycol': 0.90},
    constantSpecies=['O2'],
)

solvation(
	solvent='ethane-1,2-diol'
)

simulator(
    atol=1e-16,
    rtol=1e-8,
)

model(
    toleranceKeepInEdge=1E-20,
    toleranceMoveToCore=0.1,
    toleranceInterruptSimulation=0.1,
    maximumEdgeSpecies=100000
)

options(
    units='si',
)
