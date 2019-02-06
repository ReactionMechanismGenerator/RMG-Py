title = 'Addition of CH3 across a double bond in CH2O'

description = \
"""
This example illustrates how more complex explorer jobs work.  In this case the source channel involves two reactants
and since CH3 can add across the double bond two ways this results in two pressure dependent networks.  
"""
database(
    thermoLibraries = ['primaryThermoLibrary','Klippenstein_Glarborg2016','thermo_DFT_CCSDTF12_BAC','CBS_QB3_1dHR','DFT_QCI_thermo'],
    reactionLibraries = [],
    kineticsDepositories = ['training'],
    kineticsFamilies = 'default',
    kineticsEstimator = 'rate rules',
)



species(
label='CH3',
structure=SMILES('[CH3]'),
)

species(
label='CH2O',
structure=SMILES('O=C'),
)

species(
    label = 'N2',
    structure = SMILES('N#N'),
    molecularWeight = (28.04,"g/mol"),
    collisionModel = TransportData(sigma=(3.70,'angstrom'), epsilon=(94.9,'K')),
    reactive = False
)

pressureDependence(
    label = 'CH2O+CH3',
    Tmin = (300.0,'K'), Tmax = (1200,'K'), Tcount = 7,
    Pmin = (1.0,'atm'), Pmax = (10,'atm'), Pcount = 7,
    maximumGrainSize = (0.5,'kcal/mol'), 
    minimumGrainCount = 500, 
    method = 'modified strong collision',
    interpolationModel = ('pdeparrhenius'), 
    activeKRotor = True,
    rmgmode = False, 
)

explorer(
    source=['CH3','CH2O'],
    explore_tol=1.0e-1,
    bathGas={'N2':1.0},
    maximumRadicalElectrons=1,
)
