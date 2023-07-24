species(
    label = 'CH2OH(4)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {2,S} {5,S}
2 C u1 p0 c0 {1,S} {3,S} {4,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {1,S}
"""),
    E0 = (-33.4407,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1082.68,3535.97],'cm^-1')),
        HinderedRotor(inertia=(0.00270437,'amu*angstrom^2'), symmetry=1, barrier=(24.1813,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (31.034,'amu'),
    collisionModel = TransportData(epsilon=(4,'kJ/mol'), sigma=(3.69e-10,'m')),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.71178,0.0019307,2.12356e-05,-3.03177e-08,1.24887e-11,-4007.46,7.2919], Tmin=(100,'K'), Tmax=(895,'K')), NASAPolynomial(coeffs=[6.05623,0.00302187,1.71339e-08,-6.96092e-11,5.18063e-15,-4890.48,-6.34725], Tmin=(895,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-33.4407,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-OsHHH) + radical(CsJOH)"""),
)

species(
    label = 'methoxy(1)',
    structure = adjacencyList("""multiplicity 2
1 O u1 p2 c0 {2,S}
2 C u0 p0 c0 {1,S} {3,S} {4,S} {5,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
"""),
    E0 = (-0.0825145,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (31.034,'amu'),
    collisionModel = TransportData(epsilon=(4,'kJ/mol'), sigma=(3.69e-10,'m')),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.00135,-0.00415675,3.26351e-05,-3.71113e-08,1.35707e-11,-6.15226,6.81374], Tmin=(100,'K'), Tmax=(916.887,'K')), NASAPolynomial(coeffs=[4.01624,0.0062681,-1.58066e-06,2.44601e-10,-1.70333e-14,-449.814,4.33868], Tmin=(916.887,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-0.0825145,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-OsHHH) + radical(H3COJ)"""),
)

species(
    label = 'H(2)',
    structure = adjacencyList("""multiplicity 2
1 H u1 p0 c0
"""),
    E0 = (211.805,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (1.00797,'amu'),
    collisionModel = TransportData(epsilon=(4,'kJ/mol'), sigma=(3.69e-10,'m')),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,4.30732e-14,-5.28913e-17,2.10457e-20,-2.5714e-24,25474.2,-0.444973], Tmin=(100,'K'), Tmax=(4996.86,'K')), NASAPolynomial(coeffs=[-48.6862,0.0190728,8.49269e-07,-9.90486e-10,9.34419e-14,103971,356.208], Tmin=(4996.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(211.805,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""H""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'CH2O(3)',
    structure = adjacencyList("""1 O u0 p2 c0 {2,D}
2 C u0 p0 c0 {1,D} {3,S} {4,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (-118.609,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (30.026,'amu'),
    collisionModel = TransportData(epsilon=(4,'kJ/mol'), sigma=(3.69e-10,'m')),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.32289,-0.00506325,2.15155e-05,-1.76521e-08,4.31813e-12,-14279,2.39243], Tmin=(100,'K'), Tmax=(1402.29,'K')), NASAPolynomial(coeffs=[3.18,0.00955592,-6.27298e-06,1.33554e-09,-9.68405e-14,-15075.2,4.31051], Tmin=(1402.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-118.609,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-OdHH)"""),
)

species(
    label = 'He',
    structure = adjacencyList("""1 He u0 p1 c0
"""),
    E0 = (-6.19738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (4.0026,'amu'),
    collisionModel = TransportData(epsilon=(0.0831,'kJ/mol'), sigma=(2.55e-10,'m')),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,0.928724], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,0.928724], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-6.19738,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""He""", comment="""Thermo library: primaryThermoLibrary"""),
)

transitionState(
    label = 'TS1',
    E0 = (101.089,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (119.197,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (124.212,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['H(2)', 'CH2O(3)'],
    products = ['methoxy(1)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(1.5752e+09,'cm^3/(mol*s)'), n=1.36838, Ea=(18.6413,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), Tmax=(2500,'K'), comment="""Fitted to 50 data points; dA = *|/ 1.06822, dn = +|- 0.00866998, dEa = +|- 0.0471817 kJ/molReaction library: 'kineticsjobs'"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CH2OH(4)'],
    products = ['methoxy(1)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(5.69456e+11,'s^-1'), n=0.318893, Ea=(163.386,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), Tmax=(2500,'K'), comment="""Fitted to 50 data points; dA = *|/ 1.03007, dn = +|- 0.00389241, dEa = +|- 0.0211824 kJ/molReaction library: 'kineticsjobs'"""),
)

reaction(
    label = 'reaction1',
    reactants = ['CH2OH(4)'],
    products = ['H(2)', 'CH2O(3)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(5.47273e+10,'s^-1'), n=0.86939, Ea=(168.401,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), Tmax=(2500,'K'), comment="""Fitted to 50 data points; dA = *|/ 1.0571, dn = +|- 0.00729546, dEa = +|- 0.0397016 kJ/molReaction library: 'kineticsjobs'"""),
)

network(
    label = 'PDepNetwork #1',
    isomers = [
        'CH2OH(4)',
        'methoxy(1)',
    ],
    reactants = [
        ('H(2)', 'CH2O(3)'),
    ],
    bathGas = {
        'He': 1,
    },
)

pressureDependence(
    label = 'PDepNetwork #1',
    Tmin = (450,'K'),
    Tmax = (700,'K'),
    Tcount = 4,
    Tlist = ([450,500,678,700],'K'),
    Pmin = (0.0101325,'bar'),
    Pmax = (1013.25,'bar'),
    Pcount = 7,
    Plist = ([0.01,0.1,1,3,10,100,1000],'atm'),
    maximumGrainSize = (0.5,'kcal/mol'),
    minimumGrainCount = 500,
    method = 'modified strong collision',
    interpolationModel = ('pdeparrhenius',),
    activeKRotor = True,
    activeJRotor = True,
    rmgmode = True,
)

