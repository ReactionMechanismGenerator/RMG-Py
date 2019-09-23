title = 'methoxy decomposition to H + CH2O'

description = \
"""
This example illustrates how to manually set up an Arkane input file for a small P-dep reaction system (using only the
RRHO assumption, and without tunneling, although this can be easily implemented). Such a calculation is desirable if
the user wishes to supply experimentally determined frequencies, for example. Although some commented notes below may
be useful, see http://reactionmechanismgenerator.github.io/RMG-Py/users/arkane/index.html for more documented information about
Arkane and creating input files.
(information pertaining this file is adopted by Dames and Golden, 2013, JPCA 117 (33) 7686-96.)
"""
database(
    thermoLibraries = ['primaryThermoLibrary'],
    reactionLibraries = [],
    kineticsDepositories = ['training'],
    kineticsFamilies = 'default',
    kineticsEstimator = 'rate rules',
)

transitionState(
    label = 'TS3',
    E0 = (34.1,'kcal/mol'),  # this INCLUDES the ZPE. Note that other energy units are also possible (e.g., kJ/mol)
    spinMultiplicity = 2,
    opticalIsomers = 1, 
    frequency = (-967,'cm^-1'),
    modes = [   # these modes are used to compute the partition functions
        HarmonicOscillator(frequencies=([466,581,1169,1242,1499,1659,2933,3000],'cm^-1')),
        NonlinearRotor(rotationalConstant=([0.970, 1.029, 3.717],"cm^-1"),symmetry=1, quantum=False),
        IdealGasTranslation(mass=(31.01843,"g/mol")) #this must be included for every species/ts
    ],

)

transitionState(
    label = 'TS2',
    E0 = (38.9,'kcal/mol'),
    spinMultiplicity = 2,
    opticalIsomers = 1, 
    frequency = (-1934,'cm^-1'),
    modes = [
        HarmonicOscillator(frequencies=([792, 987 ,1136, 1142, 1482 ,2441 ,3096, 3183],'cm^-1')),
        NonlinearRotor(rotationalConstant=([0.928,0.962,5.807],"cm^-1"),symmetry=1, quantum=False),
        IdealGasTranslation(mass=(31.01843,"g/mol")) 
    ],

)
transitionState(
    label = 'TS1',
    E0 = (39.95,'kcal/mol'), 
    spinMultiplicity = 2,
    opticalIsomers = 1, 
    frequency = (-1756,'cm^-1'),
    modes = [
        HarmonicOscillator(frequencies=([186 ,626 ,1068, 1234, 1474, 1617, 2994 ,3087],'cm^-1')),
        NonlinearRotor(rotationalConstant=([0.966,0.986,5.253],"cm^-1"),symmetry=1, quantum=False),
        IdealGasTranslation(mass=(31.01843,"g/mol")) 
    ],

)

species(
    label = 'methoxy',
    structure = SMILES('C[O]'),
    E0 = (9.44,'kcal/mol'),
    modes = [
        HarmonicOscillator(frequencies=([758,960,1106 ,1393,1403,1518,2940,3019,3065],'cm^-1')),
        NonlinearRotor(rotationalConstant=([0.916, 0.921, 5.251],"cm^-1"),symmetry=3, quantum=False),
        IdealGasTranslation(mass=(31.01843,"g/mol")),
    ],
    spinMultiplicity = 3.88, # 3+exp(-89/T)
    opticalIsomers = 1,
    molecularWeight = (31.01843,'amu'),
    collisionModel = TransportData(sigma=(3.69e-10,'m'), epsilon=(4.0,'kJ/mol')),
    energyTransferModel = SingleExponentialDown(alpha0=(0.956,'kJ/mol'), T0=(300,'K'), n=0.95),
)

species(
    label = 'CH2O',
    structure = SMILES('C=O'),
    E0 = (28.69,'kcal/mol'),
    molecularWeight = (30.0106,"g/mol"),
    collisionModel = TransportData(sigma=(3.69e-10,'m'), epsilon=(4.0,'kJ/mol')),
    energyTransferModel = SingleExponentialDown(alpha0=(0.956,'kJ/mol'), T0=(300,'K'), n=0.95),
    spinMultiplicity = 1, 
    opticalIsomers = 1,
    modes = [
        HarmonicOscillator(frequencies=([1180,1261,1529,1764,2931,2999],'cm^-1')), 
        NonlinearRotor(rotationalConstant=([1.15498821005263, 1.3156969584727, 9.45570474524524],"cm^-1"),symmetry=2, quantum=False),
        IdealGasTranslation(mass=(30.0106,"g/mol")),
    ],
)

species(
    label = 'H',
    structure = SMILES('[H]'),
    E0 = (0.000,'kcal/mol'),
    molecularWeight = (1.00783,"g/mol"),
    collisionModel = TransportData(sigma=(3.69e-10,'m'), epsilon=(4.0,'kJ/mol')),
    energyTransferModel = SingleExponentialDown(alpha0=(0.956,'kJ/mol'), T0=(300,'K'), n=0.95),
    modes = [
        IdealGasTranslation(mass=(1.00783,"g/mol")),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,

)

species(
    label = 'CH2Ob',  #this is a special system with two chemically equivalent product channels. Thus, different labels are used.
    structure = SMILES('C=O'),
    E0 = (28.69,'kcal/mol'),
    molecularWeight = (30.0106,"g/mol"),
    collisionModel = TransportData(sigma=(3.69e-10,'m'), epsilon=(4.0,'kJ/mol')),
    energyTransferModel = SingleExponentialDown(alpha0=(0.956,'kJ/mol'), T0=(300,'K'), n=0.95),
    spinMultiplicity = 1, 
    opticalIsomers = 1,
    modes = [
        HarmonicOscillator(frequencies=([1180,1261,1529,1764,2931,2999],'cm^-1')), 
        NonlinearRotor(rotationalConstant=([1.15498821005263, 1.3156969584727, 9.45570474524524],"cm^-1"),symmetry=2, quantum=False),
        IdealGasTranslation(mass=(30.0106,"g/mol")),
    ],
)

species(
    label = 'Hb',
    structure = SMILES('[H]'),
    E0 = (0.0001,'kcal/mol'),
    molecularWeight = (1.00783,"g/mol"),
    collisionModel = TransportData(sigma=(3.69e-10,'m'), epsilon=(4.0,'kJ/mol')),
    energyTransferModel = SingleExponentialDown(alpha0=(0.956,'kJ/mol'), T0=(300,'K'), n=0.95),
    modes = [
        IdealGasTranslation(mass=(1.00783,"g/mol")),
    ],
    spinMultiplicity = 2, 
    opticalIsomers = 1,

)
species(
    label = 'CH2OH',
    structure = SMILES('[CH2]O'),
    E0 = (0.00,'kcal/mol'),
    molecularWeight = (31.01843,"g/mol"),	
    modes = [
        HarmonicOscillator(frequencies=([418,595, 1055, 1198, 1368, 1488, 3138, 3279, 3840],'cm^-1')),
        # below is an example of how to include hindered rotors
        #HinderedRotor(inertia=(5.75522e-47,'kg*m^2'), symmetry=1, barrier=(22427.8,'J/mol'), semiclassical=False),
        NonlinearRotor(rotationalConstant=([0.868,0.993,6.419],"cm^-1"),symmetry=1, quantum=False),
        IdealGasTranslation(mass=(31.01843,"g/mol")),
    ],
    spinMultiplicity = 2, 
    opticalIsomers = 2,
    collisionModel = TransportData(sigma=(3.69e-10,'m'), epsilon=(4.0,'kJ/mol')),
    energyTransferModel = SingleExponentialDown(alpha0=(0.956,'kJ/mol'), T0=(300,'K'), n=0.95),
)

species(
    label = 'He',
#    freqScaleFactor = 1, # TypeError: species() got an unexpected keyword argument 'freqScaleFactor'.
    structure = SMILES('[He]'),
    reactive=False,
    molecularWeight = (4.003,'amu'),
    collisionModel = TransportData(sigma=(2.55e-10,'m'), epsilon=(0.0831,'kJ/mol')),
    energyTransferModel = SingleExponentialDown(alpha0=(0.956,'kJ/mol'), T0=(300,'K'), n=0.95),
)

reaction(
    label = 'CH2O+H=Methoxy',
    reactants = ['CH2O','H'],
    products = ['methoxy'],
    transitionState = 'TS3',
)

reaction(
   label = 'CH2OH = CH2Ob+Hb',
    reactants = ['CH2OH'],
   products = ['CH2Ob', 'Hb'],
    transitionState = 'TS1',
)

reaction(
    label = 'CH2OH = Methoxy',
    products = ['methoxy'],
    reactants = ['CH2OH'],
    transitionState = 'TS2',
)

kinetics('CH2O+H=Methoxy')
kinetics('CH2OH = Methoxy')
kinetics('CH2OH = CH2Ob+Hb' )

pressureDependence(
    label = 'methoxy',
    Tmin = (450,'K'), Tmax = (1200,'K'), Tcount = 4, 
    Tlist = ([450,500,678,700],'K'),
    Pmin = (0.01,'atm'), Pmax = (1000,'atm'), Pcount = 7,
    Plist = ([0.01,0.1,1,3,10,100,1000],'atm'), 
    maximumGrainSize = (0.5,'kcal/mol'), 
    minimumGrainCount = 500, 
    method = 'modified strong collision',
    #Other methods include: 'reservoir state', 'chemically-significant eigenvalues', 
    interpolationModel = ('pdeparrhenius'), 
    activeKRotor = True,
#    active_j_rotor = False,  # causes Arkane to crash
    rmgmode = False, 
)

explorer(
    source=['methoxy'],
    explore_tol=0.01,
    energy_tol=8e1,
    flux_tol=1e-6,
     bathGas={'He':1.0},
     maximumRadicalElectrons=2,
)
