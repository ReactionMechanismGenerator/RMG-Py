species(
    label = 'O=C=C(O)C[C-]=[OH+](1267)',
    structure = adjacencyList("""1  O u0 p2 c0 {5,S} {10,S}
2  O u0 p1 c+1 {6,D} {11,S}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {5,S} {6,S} {8,S} {9,S}
5  C u0 p0 c0 {1,S} {4,S} {7,D}
6  C u0 p1 c-1 {2,D} {4,S}
7  C u0 p0 c0 {3,D} {5,D}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
"""),
    E0 = (247.8,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2120,512.5,787.5,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.530354,0.0785477,-0.000110257,7.88773e-08,-2.20479e-11,29926.4,33.1627], Tmin=(100,'K'), Tmax=(882.459,'K')), NASAPolynomial(coeffs=[13.8424,0.0182043,-7.68135e-06,1.38145e-09,-9.24765e-14,27577,-29.3952], Tmin=(882.459,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(247.8,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O4dc-C2dcH) + missing(O2d-Cdd) + group(Cs-CsCsHH) + group(Cds-(Cdd-O2d)CsOs) + group(CsJ2_singlet-CsH) + missing(Cdd-CdO2d)"""),
)

species(
    label = 'CO(15)',
    structure = adjacencyList("""1 O u0 p1 c+1 {2,T}
2 C u0 p1 c-1 {1,T}
"""),
    E0 = (-119.219,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2084.51],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.01,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(762.44,'J/mol'), sigma=(3.69,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.76,'angstroms^3'), rotrelaxcollnum=4.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.5971,-0.00102424,2.83336e-06,-1.75825e-09,3.42586e-13,-14343.2,3.45822], Tmin=(100,'K'), Tmax=(1669.94,'K')), NASAPolynomial(coeffs=[2.92797,0.00181929,-8.35302e-07,1.51268e-10,-9.88862e-15,-14292.7,6.51151], Tmin=(1669.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-119.219,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""CO""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = 'C=C(O)C=O(1218)',
    structure = adjacencyList("""1 O u0 p2 c0 {3,S} {9,S}
2 O u0 p2 c0 {5,D}
3 C u0 p0 c0 {1,S} {4,D} {5,S}
4 C u0 p0 c0 {3,D} {6,S} {7,S}
5 C u0 p0 c0 {2,D} {3,S} {8,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {1,S}
"""),
    E0 = (-287.896,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(0.96167,'amu*angstrom^2'), symmetry=1, barrier=(22.1107,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.958778,'amu*angstrom^2'), symmetry=1, barrier=(22.0442,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (72.0626,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3730.82,'J/mol'), sigma=(5.21133,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=582.75 K, Pc=59.81 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.79427,0.0406221,-2.75025e-05,3.77778e-10,4.12279e-12,-34539.5,15.0215], Tmin=(100,'K'), Tmax=(1004,'K')), NASAPolynomial(coeffs=[14.2809,0.00753011,-2.94618e-06,5.955e-10,-4.58424e-14,-37886.2,-49.4496], Tmin=(1004,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-287.896,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H)"""),
)

species(
    label = 'C=C(O)C(=O)[C-]=[OH+](1306)',
    structure = adjacencyList("""1  O u0 p2 c0 {4,S} {10,S}
2  O u0 p1 c+1 {7,D} {11,S}
3  O u0 p2 c0 {5,D}
4  C u0 p0 c0 {1,S} {5,S} {6,D}
5  C u0 p0 c0 {3,D} {4,S} {7,S}
6  C u0 p0 c0 {4,D} {8,S} {9,S}
7  C u0 p1 c-1 {2,D} {5,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
"""),
    E0 = (171.113,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.893271,0.0705198,-8.90388e-05,5.8384e-08,-1.52017e-11,20690.2,32.5997], Tmin=(100,'K'), Tmax=(938.772,'K')), NASAPolynomial(coeffs=[12.6041,0.0206207,-9.30726e-06,1.76204e-09,-1.22763e-13,18491.5,-23.158], Tmin=(938.772,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(171.113,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O4dc-C2dcH) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-OdCsCs) + group(Cds-CdsHH) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = 'C=[C-][OH+]C(O)=C=O(1307)',
    structure = adjacencyList("""1  O u0 p1 c+1 {4,S} {6,S} {8,S}
2  O u0 p2 c0 {4,S} {11,S}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {1,S} {2,S} {7,D}
5  C u0 p0 c0 {6,D} {9,S} {10,S}
6  C u0 p1 c-1 {1,S} {5,D}
7  C u0 p0 c0 {3,D} {4,D}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {2,S}
"""),
    E0 = (171.511,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73714,0.0556172,-9.31145e-05,8.32827e-08,-2.88013e-11,20703.8,16.1621], Tmin=(100,'K'), Tmax=(844.788,'K')), NASAPolynomial(coeffs=[7.04767,0.0194684,-9.39052e-06,1.79287e-09,-1.22966e-13,20199.2,-6.23847], Tmin=(844.788,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(171.511,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + missing(O2d-Cdd) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(CsJ2_singlet-CsH) + missing(Cdd-CdO2d)"""),
)

species(
    label = 'O=CC(=O)C[C-]=[OH+](1266)',
    structure = adjacencyList("""1  O u0 p1 c+1 {7,D} {11,S}
2  O u0 p2 c0 {5,D}
3  O u0 p2 c0 {6,D}
4  C u0 p0 c0 {5,S} {7,S} {8,S} {9,S}
5  C u0 p0 c0 {2,D} {4,S} {6,S}
6  C u0 p0 c0 {3,D} {5,S} {10,S}
7  C u0 p1 c-1 {1,D} {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {1,S}
"""),
    E0 = (197.431,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33214,0.065426,-9.6598e-05,8.80632e-08,-3.2853e-11,23835,33.5529], Tmin=(100,'K'), Tmax=(780.091,'K')), NASAPolynomial(coeffs=[5.26,0.0348407,-1.77035e-05,3.47644e-09,-2.44522e-13,23540,17.6157], Tmin=(780.091,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(197.431,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O4dc-C2dcH) + group(Cs-CsCsHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)H) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = 'O=C=C(O)C=[C-][OH2+](1308)',
    structure = adjacencyList("""1  O u0 p1 c+1 {6,S} {9,S} {10,S}
2  O u0 p2 c0 {4,S} {11,S}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {2,S} {5,S} {7,D}
5  C u0 p0 c0 {4,S} {6,D} {8,S}
6  C u0 p1 c-1 {1,S} {5,D}
7  C u0 p0 c0 {3,D} {4,D}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
"""),
    E0 = (147.915,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,350,440,435,1725,3010,987.5,1337.5,450,1655,2120,512.5,787.5,180,180,180,180,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.218966,'amu*angstrom^2'), symmetry=1, barrier=(5.03447,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.215219,'amu*angstrom^2'), symmetry=1, barrier=(4.94831,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.220301,'amu*angstrom^2'), symmetry=1, barrier=(5.06514,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.722975,0.0790556,-0.000136384,1.16887e-07,-3.82154e-11,17901.4,23.1492], Tmin=(100,'K'), Tmax=(873.012,'K')), NASAPolynomial(coeffs=[11.0031,0.0183746,-8.79115e-06,1.6349e-09,-1.09232e-13,16624,-22.0866], Tmin=(873.012,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(147.915,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + missing(O2d-Cdd) + group(Cds-(Cdd-O2d)CsOs) + group(Cds-CdsCsH) + group(CsJ2_singlet-CsH) + missing(Cdd-CdO2d)"""),
)

species(
    label = 'O=C=C=O(1310)',
    structure = adjacencyList("""1 O u0 p2 c0 {3,D}
2 O u0 p2 c0 {4,D}
3 C u0 p0 c0 {1,D} {4,D}
4 C u0 p0 c0 {2,D} {3,D}
"""),
    E0 = (63.731,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (56.0201,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5355,0.0240331,-4.2614e-05,3.87327e-08,-1.34307e-11,7697.1,25.5769], Tmin=(100,'K'), Tmax=(857.393,'K')), NASAPolynomial(coeffs=[4.78325,0.00773531,-3.93445e-06,7.52107e-10,-5.12482e-14,7525.26,16.3243], Tmin=(857.393,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(63.731,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(87.302,'J/(mol*K)'), label="""OCCO(S)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'C=[C-][OH2+](1309)',
    structure = adjacencyList("""1 O u0 p1 c+1 {3,S} {4,S} {5,S}
2 C u0 p0 c0 {3,D} {6,S} {7,S}
3 C u0 p1 c-1 {1,S} {2,D}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
"""),
    E0 = (322.572,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.34922,0.0191947,-3.34721e-05,3.95254e-08,-1.69764e-11,38815.2,9.15218], Tmin=(100,'K'), Tmax=(831.136,'K')), NASAPolynomial(coeffs=[0.881082,0.0184291,-9.27107e-06,1.80958e-09,-1.26146e-13,39662.2,23.2302], Tmin=(831.136,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(322.572,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cds-CdsHH) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = 'N2',
    structure = adjacencyList("""1 N u0 p1 c0 {2,T}
2 N u0 p1 c0 {1,T}
"""),
    E0 = (-8.69489,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0137,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(810.913,'J/mol'), sigma=(3.621,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.76,'angstroms^3'), rotrelaxcollnum=4.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.61263,-0.00100893,2.49899e-06,-1.43376e-09,2.58637e-13,-1051.1,2.6527], Tmin=(100,'K'), Tmax=(1817.03,'K')), NASAPolynomial(coeffs=[2.97588,0.00164144,-7.19736e-07,1.25381e-10,-7.91548e-15,-1025.83,5.5377], Tmin=(1817.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.69489,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""N2""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = 'Ne',
    structure = adjacencyList("""1 Ne u0 p4 c0
"""),
    E0 = (-6.19738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (20.1801,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1235.53,'J/mol'), sigma=(3.758e-10,'m'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with fixed Lennard Jones Parameters. This is the fallback method! Try improving transport databases!"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,3.35532], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,3.35532], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-6.19738,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""Ne""", comment="""Thermo library: primaryThermoLibrary"""),
)

transitionState(
    label = 'TS1',
    E0 = (100.023,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (192.957,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (191.609,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (210.828,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (191.339,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (234.935,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O=C=C(O)C[C-]=[OH+](1267)'],
    products = ['CO(15)', 'C=C(O)C=O(1218)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(1.69762e+27,'s^-1'), n=-4.4977, Ea=(5.68115,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.14486523882846947, var=132.63092240559882, Tref=1000.0, N=66, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root in family Retroene."""),
)

reaction(
    label = 'reaction2',
    reactants = ['O=C=C(O)C[C-]=[OH+](1267)'],
    products = ['C=C(O)C(=O)[C-]=[OH+](1306)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(7.50839e+11,'s^-1'), n=0.136197, Ea=(98.6152,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.01703312323488117, var=1.3519057502930525, Tref=1000.0, N=4, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R in family 1,3_sigmatropic_rearrangement."""),
)

reaction(
    label = 'reaction3',
    reactants = ['O=C=C(O)C[C-]=[OH+](1267)'],
    products = ['C=[C-][OH+]C(O)=C=O(1307)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(7.50839e+11,'s^-1'), n=0.136197, Ea=(97.2673,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.01703312323488117, var=1.3519057502930525, Tref=1000.0, N=4, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R in family 1,3_sigmatropic_rearrangement."""),
)

reaction(
    label = 'reaction4',
    reactants = ['O=C=C(O)C[C-]=[OH+](1267)'],
    products = ['O=CC(=O)C[C-]=[OH+](1266)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4.7171e-39,'s^-1'), n=14.8922, Ea=(116.486,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.01101064541230846, var=9.816149634195241, Tref=1000.0, N=17, data_mean=0.0, correlation='Root_3R!H->C_N-3C-inRing_N-1R!H->N_Ext-3C-R',), comment="""Estimated from node Root_3R!H->C_N-3C-inRing_N-1R!H->N_Ext-3C-R in family Ketoenol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['O=C=C(O)C=[C-][OH2+](1308)'],
    products = ['O=C=C(O)C[C-]=[OH+](1267)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.90624e-39,'s^-1'), n=15.0595, Ea=(196.882,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.0038146995811367415, var=4.810563762111573, Tref=1000.0, N=6, data_mean=0.0, correlation='Root_3R!H->C_N-3C-inRing_N-1R!H->N_Ext-3C-R_Ext-5R!H-R_N-6R!H->O',), comment="""Estimated from node Root_3R!H->C_N-3C-inRing_N-1R!H->N_Ext-3C-R_Ext-5R!H-R_N-6R!H->O in family Ketoenol.
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction6',
    reactants = ['O=C=C(O)C[C-]=[OH+](1267)'],
    products = ['O=C=C=O(1310)', 'C=[C-][OH2+](1309)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.46545e+14,'s^-1'), n=0.20628, Ea=(140.593,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(3000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R!H->C',), comment="""Estimated from node Root_N-4R!H->C in family Retroene."""),
)

network(
    label = 'PDepNetwork #883',
    isomers = [
        'O=C=C(O)C[C-]=[OH+](1267)',
    ],
    reactants = [
        ('CO(15)', 'C=C(O)C=O(1218)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #883',
    Tmin = (300,'K'),
    Tmax = (2200,'K'),
    Tcount = 8,
    Tlist = ([302.51,323.546,371.247,459.823,619.914,913.864,1434.46,2073.82],'K'),
    Pmin = (0.01,'bar'),
    Pmax = (100,'bar'),
    Pcount = 5,
    Plist = ([0.0125282,0.0667467,1,14.982,79.8202],'bar'),
    maximumGrainSize = (0.5,'kcal/mol'),
    minimumGrainCount = 250,
    method = 'modified strong collision',
    interpolationModel = ('Chebyshev', 6, 4),
    activeKRotor = True,
    activeJRotor = True,
    rmgmode = True,
)

