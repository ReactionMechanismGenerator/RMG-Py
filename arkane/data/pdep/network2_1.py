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
    label = '[OH](8)',
    structure = adjacencyList("""multiplicity 2
1 O u1 p2 c0 {2,S}
2 H u0 p0 c0 {1,S}
"""),
    E0 = (28.3945,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3668.68],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (17.0074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.51457,2.92774e-05,-5.32163e-07,1.01949e-09,-3.85945e-13,3414.25,2.10435], Tmin=(100,'K'), Tmax=(1145.75,'K')), NASAPolynomial(coeffs=[3.07194,0.000604016,-1.39783e-08,-2.13446e-11,2.48066e-15,3579.39,4.578], Tmin=(1145.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.3945,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""OH(D)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH2](9)',
    structure = adjacencyList("""multiplicity 3
1 C u2 p0 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
"""),
    E0 = (381.37,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1066.91,2790.98,3622.37],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.01192,-0.000154978,3.26298e-06,-2.40422e-09,5.69497e-13,45867.7,0.533201], Tmin=(100,'K'), Tmax=(1104.62,'K')), NASAPolynomial(coeffs=[3.14983,0.00296674,-9.76056e-07,1.54115e-10,-9.50338e-15,46058.1,4.77808], Tmin=(1104.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH]O(10)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {2,S} {3,S}
2 C u2 p0 c0 {1,S} {4,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (205.907,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,402.687,3356.16],'cm^-1')),
        HinderedRotor(inertia=(0.0105041,'amu*angstrom^2'), symmetry=1, barrier=(23.1304,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (30.026,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.76003,0.00295745,8.86361e-06,-1.33922e-08,5.33443e-12,24775.7,6.76103], Tmin=(100,'K'), Tmax=(943.113,'K')), NASAPolynomial(coeffs=[5.07487,0.00326008,-9.68494e-07,1.67782e-10,-1.21781e-14,24266.2,-0.891511], Tmin=(943.113,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(205.907,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(78.9875,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-OsHHH) + radical(OsCsJ2H_triplet)"""),
)

species(
    label = '[CH2][O](7)',
    structure = adjacencyList("""multiplicity 3
1 O u1 p2 c0 {2,S}
2 C u1 p0 c0 {1,S} {3,S} {4,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (185.5,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (30.026,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.8129,-0.00139725,2.72241e-05,-3.809e-08,1.59451e-11,22322.5,7.75997], Tmin=(100,'K'), Tmax=(884.263,'K')), NASAPolynomial(coeffs=[6.9814,-0.00114411,2.05208e-06,-4.58163e-10,3.18292e-14,21191.9,-10.3609], Tmin=(884.263,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(185.5,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-OsHHH) + radical(H3COJ) + radical(CsJOH)"""),
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
    E0 = (-0.717018,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-5.75102,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (279.245,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (282.026,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (261.618,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['CH2OH(4)'],
    products = ['H(2)', 'CH2O(3)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5.51244e+10,'s^-1'), n=0.868564, Ea=(168.41,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2500,'K'), comment="""Fitted to 59 data points; dA = *|/ 1.05152, dn = +|- 0.00659302, dEa = +|- 0.0362682 kJ/molReaction library: 'kineticsjobs'"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CH2OH(4)'],
    products = ['methoxy(1)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(5.63501e+11,'s^-1'), n=0.320211, Ea=(163.376,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2500,'K'), comment="""Fitted to 59 data points; dA = *|/ 1.02731, dn = +|- 0.00353557, dEa = +|- 0.0194492 kJ/molReaction library: 'kineticsjobs'"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[OH](8)', '[CH2](9)'],
    products = ['CH2OH(4)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(15.4803,'m^3/(mol*s)'), n=1.88017, Ea=(5.1666,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [O_rad;Birad] for rate rule [O_pri_rad;Birad]
Euclidian distance = 1.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(2)', '[CH]O(10)'],
    products = ['CH2OH(4)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(2)', '[CH2][O](7)'],
    products = ['CH2OH(4)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.80515e+06,'m^3/(mol*s)'), n=0.314888, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O"""),
)

network(
    label = 'PDepNetwork #2',
    isomers = [
        'CH2OH(4)',
    ],
    reactants = [
        ('H(2)', 'CH2O(3)'),
    ],
    bathGas = {
        'He': 1,
    },
)

pressureDependence(
    label = 'PDepNetwork #2',
    Tmin = (450,'K'),
    Tmax = (1200,'K'),
    Tcount = 4,
    Tlist = ([450,800,1000,1200],'K'),
    Pmin = (0.0101325,'bar'),
    Pmax = (1013.25,'bar'),
    Pcount = 3,
    Plist = ([0.01,1,1000],'atm'),
    maximumGrainSize = (0.5,'kcal/mol'),
    minimumGrainCount = 500,
    method = 'modified strong collision',
    interpolationModel = ('pdeparrhenius',),
    activeKRotor = True,
    activeJRotor = True,
    rmgmode = True,
)

