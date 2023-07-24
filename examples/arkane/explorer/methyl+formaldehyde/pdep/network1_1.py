species(
    label = '[CH2]OC(5)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {2,S} {3,S}
2 C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
3 C u1 p0 c0 {1,S} {7,S} {8,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {3,S}
"""),
    E0 = (-10.4064,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,1747.94],'cm^-1')),
        HinderedRotor(inertia=(0.112552,'amu*angstrom^2'), symmetry=1, barrier=(5.6503,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112551,'amu*angstrom^2'), symmetry=1, barrier=(5.6503,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (45.0605,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2741.03,'J/mol'), sigma=(4.97958,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=428.14 K, Pc=50.37 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.19193,0.015635,6.95077e-09,-5.41305e-09,1.83143e-12,-1220.84,11.0224], Tmin=(100,'K'), Tmax=(1282.76,'K')), NASAPolynomial(coeffs=[4.8276,0.0158408,-6.43852e-06,1.16148e-09,-7.83638e-14,-2077.05,1.02219], Tmin=(1282.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-10.4064,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-OsHHH) + group(Cs-OsHHH) + radical(CsJOCH3)"""),
)

species(
    label = 'CH2O(2)',
    structure = adjacencyList("""1 O u0 p2 c0 {2,D}
2 C u0 p0 c0 {1,D} {3,S} {4,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (-118.547,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (30.026,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = ThermoData(Tdata=([300,400,500,600,800,1000,1500],'K'), Cpdata=([35.4385,39.2459,43.7646,48.1997,55.9401,61.965,61.965],'J/(mol*K)'), H298=(-108.575,'kJ/mol'), S298=(218.834,'J/(mol*K)'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), E0=(-118.547,'kJ/mol'), comment="""Thermo group additivity estimation: group(Cds-OdHH)"""),
)

species(
    label = 'CH3(1)',
    structure = adjacencyList("""multiplicity 2
1 C u1 p0 c0 {2,S} {3,S} {4,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
"""),
    E0 = (135.366,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([546.865,1439.42,1439.62,4000,4000,4000],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (15.0346,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = ThermoData(Tdata=([300,400,500,600,800,1000,1500],'K'), Cpdata=([9.331,9.988,10.611,11.274,12.457,13.427,15.105],'cal/(mol*K)'), H298=(34.893,'kcal/mol'), S298=(46.3704,'cal/(mol*K)'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), E0=(135.366,'kJ/mol'), comment="""Thermo library: primaryThermoLibrary + radical(CH3)"""),
)

species(
    label = '[CH2](7)',
    structure = adjacencyList("""1 C u0 p1 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
"""),
    E0 = (419.091,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1369.93,2896,2896.02],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.10264,-0.00144066,5.45062e-06,-3.57996e-09,7.56173e-13,50400.6,-0.411757], Tmin=(100,'K'), Tmax=(1442.38,'K')), NASAPolynomial(coeffs=[2.62652,0.00394756,-1.49921e-06,2.54531e-10,-1.6295e-14,50691.7,6.78352], Tmin=(1442.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(419.091,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(S)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH2]O(8)',
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.71178,0.0019307,2.12356e-05,-3.03177e-08,1.24887e-11,-4007.46,7.2919], Tmin=(100,'K'), Tmax=(895,'K')), NASAPolynomial(coeffs=[6.05623,0.00302187,1.71339e-08,-6.96092e-11,5.18063e-15,-4890.48,-6.34725], Tmin=(895,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-33.4407,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-OsHHH) + radical(CsJOH)"""),
)

species(
    label = 'N2',
    structure = adjacencyList("""1 N u0 p1 c0 {2,T}
2 N u0 p1 c0 {1,T}
"""),
    E0 = (-8.64289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0137,'amu'),
    collisionModel = TransportData(epsilon=(789.043,'J/mol'), sigma=(3.7,'angstrom')),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53101,-0.000123661,-5.02999e-07,2.43531e-09,-1.40881e-12,-1046.98,2.96747], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.95258,0.0013969,-4.92632e-07,7.8601e-11,-4.60755e-15,-923.949,5.87189], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-8.64289,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""N2""", comment="""Thermo library: primaryThermoLibrary"""),
)

transitionState(
    label = 'TS1',
    E0 = (-41.7742,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (297.369,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['CH2O(2)', 'CH3(1)'],
    products = ['[CH2]OC(5)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(0.00491391,'m^3/(mol*s)'), n=2.86227, Ea=(29.6879,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;CsJ-HHH] + [Od_CO-HH;YJ] for rate rule [Od_CO-HH;CsJ-HHH]
Euclidian distance = 3.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2](7)', '[CH2]O(8)'],
    products = ['[CH2]OC(5)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(71881.9,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

network(
    label = 'PDepNetwork #1',
    isomers = [
        '[CH2]OC(5)',
    ],
    reactants = [
        ('CH2O(2)', 'CH3(1)'),
    ],
    bathGas = {
        'N2': 1,
    },
)

pressureDependence(
    label = 'PDepNetwork #1',
    Tmin = (300,'K'),
    Tmax = (1200,'K'),
    Tcount = 7,
    Tlist = ([1200,800,600,480,400,342.857,300],'K'),
    Pmin = (1,'atm'),
    Pmax = (10,'atm'),
    Pcount = 7,
    Plist = ([1.01325,1.48725,2.18298,3.20418,4.70309,6.90319,10.1325],'bar'),
    maximumGrainSize = (0.5,'kcal/mol'),
    minimumGrainCount = 500,
    method = 'modified strong collision',
    interpolationModel = ('pdeparrhenius',),
    activeKRotor = True,
    activeJRotor = True,
    rmgmode = True,
)

