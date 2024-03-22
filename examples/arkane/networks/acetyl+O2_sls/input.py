################################################################################
#
#   Arkane input file for acetyl + O2 pressure-dependent reaction network
#
################################################################################

title = 'acetyl + oxygen'

description = \
"""
The chemically-activated reaction of acetyl with oxygen. This system is of
interest in atmospheric chemistry as a step in the conversion of acetaldehyde
to the secondary pollutant peroxyacetylnitrate (PAN); it is also potentially
important in the ignition chemistry of ethanol.
"""

species(
    label = 'acetylperoxy',
    structure = SMILES('CC(=O)O[O]'),
    E0 = (-34.6,'kcal/mol'),
    modes = [
        IdealGasTranslation(mass=(75.04,"g/mol")),
        NonlinearRotor(inertia=([54.2977,104.836,156.05],"amu*angstrom^2"), symmetry=1),
        HarmonicOscillator(frequencies=([319.695,500.474,536.674,543.894,727.156,973.365,1037.77,1119.72,1181.55,1391.11,1449.53,1454.72,1870.51,3037.12,3096.93,3136.39],"cm^-1")),
        HinderedRotor(inertia=(7.38359,"amu*angstrom^2"), symmetry=1, fourier=([[-1.95191,-11.8215,0.740041,-0.049118,-0.464522],[0.000227764,0.00410782,-0.000805364,-0.000548218,-0.000266277]],"kJ/mol")),
        HinderedRotor(inertia=(2.94723,"amu*angstrom^2"), symmetry=3, fourier=([[0.130647,0.0401507,-2.54582,-0.0436065,-0.120982],[-0.000701659,-0.000989654,0.00783349,-0.00140978,-0.00145843]],"kJ/mol")),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (75.04,"g/mol"),
    collisionModel = TransportData(sigma=(5.09,'angstrom'), epsilon=(473,'K')),
    energyTransferModel = SingleExponentialDown(
        alpha0 = (0.5718,'kcal/mol'),
        T0 = (300,'K'),
        n = 0.85,
    ),
)

species(
    label = 'hydroperoxylvinoxy',
    structure = SMILES('[CH2]C(=O)OO'),
    E0 = (-32.4,'kcal/mol'),
    modes = [
        IdealGasTranslation(mass=(75.04,"g/mol")),
        NonlinearRotor(inertia=([44.8034,110.225,155.029],"u*angstrom**2"), symmetry=1),
        HarmonicOscillator(frequencies=([318.758,420.907,666.223,675.962,752.824,864.66,998.471,1019.54,1236.21,1437.91,1485.74,1687.9,3145.44,3262.88,3434.34],"cm^-1")),
        HinderedRotor(inertia=(1.68464,"u*angstrom**2"), symmetry=2, fourier=([[0.359649,-16.1155,-0.593311,1.72918,0.256314],[-7.42981e-06,-0.000238057,3.29276e-05,-6.62608e-05,8.8443e-05]],"kJ/mol")),
        HinderedRotor(inertia=(8.50433,"u*angstrom**2"), symmetry=1, fourier=([[-7.53504,-23.4471,-3.3974,-0.940593,-0.313674],[-4.58248,-2.47177,0.550012,1.03771,0.844226]],"kJ/mol")),
        HinderedRotor(inertia=(0.803309,"u*angstrom**2"), symmetry=1, fourier=([[-8.65946,-3.97888,-1.13469,-0.402197,-0.145101],[4.41884e-05,4.83249e-05,1.30275e-05,-1.31353e-05,-6.66878e-06]],"kJ/mol")),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (75.04,"g/mol"),
    collisionModel = TransportData(sigma=(5.09,'angstrom'), epsilon=(473,'K')),
    energyTransferModel = SingleExponentialDown(
        alpha0 = (0.5718,'kcal/mol'),
        T0 = (300,'K'),
        n = 0.85,
    ),
)

species(
    label = 'acetyl',
    structure = SMILES('C[C]=O'),
    E0 = (0.0,'kcal/mol'),  #(-20.5205,"kJ/mol")
    modes = [
        IdealGasTranslation(mass=(43.05,"g/mol")),
        NonlinearRotor(inertia=([5.94518,50.8166,53.6436],"u*angstrom**2"), symmetry=1),
        HarmonicOscillator(frequencies=([464.313,845.126,1010.54,1038.43,1343.54,1434.69,1442.25,1906.18,2985.46,3076.57,3079.46],"cm^-1")),
        HinderedRotor(inertia=(1.61752,"u*angstrom**2"), symmetry=3, barrier=(2.00242,"kJ/mol")),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
)

species(
    label = 'oxygen',
    structure = SMILES('[O][O]'),
    E0 = (0.0,'kcal/mol'),  #(-5.74557,"kJ/mol")
    modes = [
        IdealGasTranslation(mass=(32.00,"g/mol")),
        LinearRotor(inertia=(11.6056,"u*angstrom**2"), symmetry=2),
        HarmonicOscillator(frequencies=([1621.54],"cm^-1")),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
)

species(
    label = 'ketene',
    structure = SMILES('C=C=O'),
    E0 = (-6.6,'kcal/mol'),
    modes = [
        IdealGasTranslation(mass=(42.0106,"g/mol")),
        NonlinearRotor(inertia=([1.76922,48.8411,50.6103],"u*angstrom**2"), symmetry=2),
        HarmonicOscillator(frequencies=([441.622,548.317,592.155,981.379,1159.66,1399.86,2192.1,3150.02,3240.58],"cm^-1")),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

species(
    label = 'lactone',
    structure = SMILES('C1OC1(=O)'),
    E0 = (-30.8,'kcal/mol'),
    modes = [
        IdealGasTranslation(mass=(58.0055,"g/mol")),
        NonlinearRotor(inertia=([20.1707,62.4381,79.1616],"u*angstrom**2"), symmetry=1),
        HarmonicOscillator(frequencies=([484.387,527.771,705.332,933.372,985.563,1050.26,1107.93,1176.43,1466.54,1985.09,3086.23,3186.46],"cm^-1")),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
)



species(
    label = 'hydroxyl',
    structure = SMILES('[OH]'),
    E0 = (0.0,'kcal/mol'),
    modes = [
        IdealGasTranslation(mass=(17.0027,"g/mol")),
        LinearRotor(inertia=(0.899988,"u*angstrom**2"), symmetry=1),
        HarmonicOscillator(frequencies=([3676.39],"cm^-1")),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
)

species(
    label = 'hydroperoxyl',
    structure = SMILES('O[O]'),
    E0 = (0.0,'kcal/mol'),
    modes = [
        IdealGasTranslation(mass=(32.9977,"g/mol")),
        NonlinearRotor(inertia=([0.811564,14.9434,15.755],"u*angstrom**2"), symmetry=1),
        HarmonicOscillator(frequencies=([1156.77,1424.28,3571.81],"cm^-1")),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
)

species(
    label = 'nitrogen',
    structure = SMILES('N#N'),
    molecularWeight = (28.04,"g/mol"),
    collisionModel = TransportData(sigma=(3.70,'angstrom'), epsilon=(94.9,'K')),
    reactive = False
)

################################################################################

transitionState(
    label = 'entrance1',
    E0 = (0.0,'kcal/mol'),
)

transitionState(
    label = 'isom1',
    E0 = (-5.8,'kcal/mol'),
    modes = [
        IdealGasTranslation(mass=(75.04,"g/mol")),
        NonlinearRotor(inertia=([49.3418,103.697,149.682],"u*angstrom**2"), symmetry=1, quantum=False),
        HarmonicOscillator(frequencies=([148.551,306.791,484.573,536.709,599.366,675.538,832.594,918.413,1022.28,1031.45,1101.01,1130.05,1401.51,1701.26,1844.17,3078.6,3163.07],"cm^-1"), quantum=True),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    frequency = (-1679.04,'cm^-1'),
)

transitionState(
    label = 'exit1',
    E0 = (0.6,'kcal/mol'),
    modes = [
        IdealGasTranslation(mass=(75.04,"g/mol")),
        NonlinearRotor(inertia=([55.4256, 136.1886, 188.2442],"u*angstrom**2"), symmetry=1),
        HarmonicOscillator(frequencies=([59.9256,204.218,352.811,466.297,479.997,542.345,653.897,886.657,1017.91,1079.17,1250.02,1309.14,1370.36,1678.52,2162.41,3061.53,3135.55],"cm^-1")),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    frequency=(-1053.25,'cm^-1'),
)

transitionState(
    label = 'exit2',
    E0 = (-4.6,'kcal/mol'),
    modes = [
        IdealGasTranslation(mass=(75.0082,"g/mol"), quantum=False),
        NonlinearRotor(inertia=([51.7432,119.373,171.117],"u*angstrom**2"), symmetry=1, quantum=False),
        HarmonicOscillator(frequencies=([250.311,383.83,544.382,578.988,595.324,705.422,964.712,1103.25,1146.91,1415.98,1483.33,1983.79,3128,3143.84,3255.29],"cm^-1")),
        HinderedRotor(inertia=(9.35921,"u*angstrom**2"), symmetry=1, fourier=([[-11.2387,-12.5928,-3.87844,1.13314,-0.358812],[-1.59863,-8.03329,-5.05235,3.13723,2.45989]],"kJ/mol")),
        HinderedRotor(inertia=(0.754698,"u*angstrom**2"), symmetry=1, barrier=(47.7645,"kJ/mol")),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    frequency=(-404.271,'cm^-1'),
)

transitionState(
    label = 'exit3',
    E0 = (-7.2,'kcal/mol'),
    modes = [
        IdealGasTranslation(mass=(75.04,"g/mol")),
        NonlinearRotor(inertia=([53.2821, 120.4050, 170.1570],"u*angstrom**2"), symmetry=1),
        HarmonicOscillator(frequencies=([152.155,181.909,311.746,348.646,608.487,624.378,805.347,948.875,995.256,996.982,1169.1,1412.6,1834.83,3124.43,3245.2,3634.45],"cm^-1")),
        HinderedRotor(inertia=(0.813269,"u*angstrom**2"), symmetry=1, fourier=([[-1.15338,-2.18664,-0.669531,-0.11502,-0.0512599],[0.00245222,0.0107485,0.00334564,-0.00165288,-0.0028674]],"kJ/mol")),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    frequency=(-618.234,'cm^-1'),
)

################################################################################

reaction(
    label = 'entrance1',
    reactants = ['acetyl', 'oxygen'],
    products = ['acetylperoxy'],
    transitionState = 'entrance1',
    kinetics = Arrhenius(A=(2.65e6,'m^3/(mol*s)'), n=0.0, Ea=(0.0,'kcal/mol'), T0=(1,"K")),
)

reaction(
    label = 'isom1',
    reactants = ['acetylperoxy'],
    products = ['hydroperoxylvinoxy'],
    transitionState = 'isom1',
    tunneling = 'Eckart',
)

reaction(
    label = 'exit1',
    reactants = ['acetylperoxy'],
    products = ['ketene', 'hydroperoxyl'],
    transitionState = 'exit1',
    tunneling = 'Eckart',
)

reaction(
    label = 'exit2',
    reactants = ['hydroperoxylvinoxy'],
    products = ['ketene', 'hydroperoxyl'],
    transitionState = 'exit2',
    tunneling = 'Eckart',
)

reaction(
    label = 'exit3',
    reactants = ['hydroperoxylvinoxy'],
    products = ['lactone', 'hydroxyl'],
    transitionState = 'exit3',
    tunneling = 'Eckart',
)

#################################################################################

network(
    label = 'acetyl + O2',
    isomers = [
        'acetylperoxy',
        'hydroperoxylvinoxy',
    ],
    reactants = [
        ('acetyl', 'oxygen'),
        ('ketene', 'hydroperoxyl'),
	('lactone', 'hydroxyl')
    ],
    bathGas = {
        'nitrogen': 1.0,
    }
)

#################################################################################

pressureDependence(
    'acetyl + O2',
    Tmin=(300.0,'K'), Tmax=(2000.0,'K'), Tcount=8,
    Pmin=(0.01,'bar'), Pmax=(100.0,'bar'), Pcount=5,
    maximumGrainSize = (1.0,'kcal/mol'),
    minimumGrainCount = 250,
    method = 'simulation least squares',
    interpolationModel = ('chebyshev', 6, 4),
    activeJRotor = True,
)
