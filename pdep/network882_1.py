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
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,2782.5,750,1395,475,1775,1000,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33214,0.065426,-9.6598e-05,8.80632e-08,-3.2853e-11,23835,33.5529], Tmin=(100,'K'), Tmax=(780.091,'K')), NASAPolynomial(coeffs=[5.26,0.0348407,-1.77035e-05,3.47644e-09,-2.44522e-13,23540,17.6157], Tmin=(780.091,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(197.431,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O4dc-C2dcH) + group(Cs-CsCsHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)H) + group(CsJ2_singlet-CsH)"""),
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
    label = 'O=CC[C-]=[OH+](1048)',
    structure = adjacencyList("""1 O u0 p1 c+1 {5,D} {9,S}
2 O u0 p2 c0 {4,D}
3 C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
4 C u0 p0 c0 {2,D} {3,S} {8,S}
5 C u0 p1 c-1 {1,D} {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {1,S}
"""),
    E0 = (305.423,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,180,180,180,1118.95,1119.04,1121.12],'cm^-1')),
        HinderedRotor(inertia=(0.00320101,'amu*angstrom^2'), symmetry=1, barrier=(2.83691,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.49418,'amu*angstrom^2'), symmetry=1, barrier=(57.3461,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (72.0626,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.43187,0.0380089,-4.41375e-05,3.87369e-08,-1.5179e-11,36787,27.1749], Tmin=(100,'K'), Tmax=(745.21,'K')), NASAPolynomial(coeffs=[3.51782,0.0274308,-1.32858e-05,2.58503e-09,-1.82039e-13,36757,23.1399], Tmin=(745.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(305.423,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O4dc-C2dcH) + group(Cs-CsCsHH) + group(Cds-OdCsH) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = 'C=C(C=O)[O+]=[C-]O(1311)',
    structure = adjacencyList("""1  O u0 p1 c+1 {4,S} {7,D}
2  O u0 p2 c0 {7,S} {11,S}
3  O u0 p2 c0 {6,D}
4  C u0 p0 c0 {1,S} {5,D} {6,S}
5  C u0 p0 c0 {4,D} {8,S} {9,S}
6  C u0 p0 c0 {3,D} {4,S} {10,S}
7  C u0 p1 c-1 {1,D} {2,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {2,S}
"""),
    E0 = (60.8196,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.999014,0.0823556,-0.000163167,1.68999e-07,-6.47619e-11,7407.17,34.6127], Tmin=(100,'K'), Tmax=(847.798,'K')), NASAPolynomial(coeffs=[1.38248,0.0417223,-2.2584e-05,4.43612e-09,-3.07274e-13,8737.41,41.0548], Tmin=(847.798,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(60.8196,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = 'O=C=COC[C-]=[OH+](1312)',
    structure = adjacencyList("""1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p1 c+1 {6,D} {11,S}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
5  C u0 p0 c0 {1,S} {7,D} {10,S}
6  C u0 p1 c-1 {2,D} {4,S}
7  C u0 p0 c0 {3,D} {5,D}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {2,S}
"""),
    E0 = (294.132,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,2120,512.5,787.5,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.607078,0.072607,-8.66195e-05,5.14444e-08,-1.18985e-11,35499.9,32.8037], Tmin=(100,'K'), Tmax=(1063.54,'K')), NASAPolynomial(coeffs=[15.7277,0.0157378,-6.41197e-06,1.16722e-09,-8.01772e-14,32283.6,-41.0761], Tmin=(1063.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(294.132,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O4dc-C2dcH) + missing(O2d-Cdd) + group(Cs-CsOsHH) + group(Cds-(Cdd-O2d)OsH) + group(CsJ2_singlet-CsH) + missing(Cdd-CdO2d)"""),
)

species(
    label = 'C=[C-][OH+]C(=O)C=O(1313)',
    structure = adjacencyList("""1  O u0 p1 c+1 {4,S} {7,S} {8,S}
2  O u0 p2 c0 {4,D}
3  O u0 p2 c0 {5,D}
4  C u0 p0 c0 {1,S} {2,D} {5,S}
5  C u0 p0 c0 {3,D} {4,S} {9,S}
6  C u0 p0 c0 {7,D} {10,S} {11,S}
7  C u0 p1 c-1 {1,S} {6,D}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (253.939,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2782.5,750,1395,475,1775,1000,2950,3100,1380,975,1025,1650,180,180,180,180,969.817,1705.86,2561.32,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0963874,'amu*angstrom^2'), symmetry=1, barrier=(2.21614,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0963874,'amu*angstrom^2'), symmetry=1, barrier=(2.21614,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0963874,'amu*angstrom^2'), symmetry=1, barrier=(2.21614,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.74915,0.0602065,-0.000109873,1.12311e-07,-4.32424e-11,30612.5,19.3547], Tmin=(100,'K'), Tmax=(838.579,'K')), NASAPolynomial(coeffs=[2.18191,0.0332966,-1.72959e-05,3.38056e-09,-2.34698e-13,31413.5,22.5518], Tmin=(838.579,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(253.939,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cds-CdsCsCs) + group(Cds-O2d(Cds-O2d)H) + group(Cds-CdsHH) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = 'O=CC(O)=C[C-]=[OH+](1314)',
    structure = adjacencyList("""1  O u0 p2 c0 {4,S} {10,S}
2  O u0 p1 c+1 {7,D} {11,S}
3  O u0 p2 c0 {6,D}
4  C u0 p0 c0 {1,S} {5,D} {6,S}
5  C u0 p0 c0 {4,D} {7,S} {8,S}
6  C u0 p0 c0 {3,D} {4,S} {9,S}
7  C u0 p1 c-1 {2,D} {5,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
"""),
    E0 = (183.034,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,2782.5,750,1395,475,1775,1000,375.879,375.88,375.88,375.88,375.88,375.88],'cm^-1')),
        HinderedRotor(inertia=(0.239655,'amu*angstrom^2'), symmetry=1, barrier=(24.0276,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0884013,'amu*angstrom^2'), symmetry=1, barrier=(8.86313,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.361358,'amu*angstrom^2'), symmetry=1, barrier=(36.2295,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.849025,0.0748716,-0.000102206,7.34563e-08,-2.12723e-11,22122.5,31.5125], Tmin=(100,'K'), Tmax=(840.251,'K')), NASAPolynomial(coeffs=[11.3615,0.0248282,-1.2872e-05,2.57939e-09,-1.84871e-13,20355.8,-17.3749], Tmin=(840.251,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(183.034,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O4dc-C2dcH) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)H) + group(CsJ2_singlet-CsH)"""),
)

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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3319.97,'J/mol'), sigma=(5.949e-10,'m'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with fixed Lennard Jones Parameters. This is the fallback method! Try improving transport databases!"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.530354,0.0785477,-0.000110257,7.88773e-08,-2.20479e-11,29926.4,33.1627], Tmin=(100,'K'), Tmax=(882.459,'K')), NASAPolynomial(coeffs=[13.8424,0.0182043,-7.68135e-06,1.38145e-09,-9.24765e-14,27577,-29.3952], Tmin=(882.459,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(247.8,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O4dc-C2dcH) + missing(O2d-Cdd) + group(Cs-CsCsHH) + group(Cds-(Cdd-O2d)CsOs) + group(CsJ2_singlet-CsH) + missing(Cdd-CdO2d)"""),
)

species(
    label = 'O=CC(=O)C=[C-][OH2+](1315)',
    structure = adjacencyList("""1  O u0 p1 c+1 {7,S} {10,S} {11,S}
2  O u0 p2 c0 {4,D}
3  O u0 p2 c0 {6,D}
4  C u0 p0 c0 {2,D} {5,S} {6,S}
5  C u0 p0 c0 {4,S} {7,D} {8,S}
6  C u0 p0 c0 {3,D} {4,S} {9,S}
7  C u0 p1 c-1 {1,S} {5,D}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
"""),
    E0 = (97.546,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,375,552.5,462.5,1710,3010,987.5,1337.5,450,1655,2782.5,750,1395,475,1775,1000,180,180,180,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.11563,'amu*angstrom^2'), symmetry=1, barrier=(2.65856,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114876,'amu*angstrom^2'), symmetry=1, barrier=(2.64123,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114844,'amu*angstrom^2'), symmetry=1, barrier=(2.64048,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52868,0.0658628,-0.000122308,1.25151e-07,-4.83736e-11,11809.9,23.5267], Tmin=(100,'K'), Tmax=(829.379,'K')), NASAPolynomial(coeffs=[2.51602,0.034842,-1.87128e-05,3.7056e-09,-2.59228e-13,12549.2,24.3927], Tmin=(829.379,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(97.546,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-O2d)H) + group(CsJ2_singlet-CsH)"""),
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
    E0 = (61.2981,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (345.639,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (127.985,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (484.952,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (216.464,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (189.609,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (218.415,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (147.45,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (332.838,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O=CC(=O)C[C-]=[OH+](1266)'],
    products = ['CO(15)', 'C=C(O)C=O(1218)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(1.69762e+27,'s^-1'), n=-4.4977, Ea=(9.73853,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.14486523882846947, var=132.63092240559882, Tref=1000.0, N=66, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root in family Retroene."""),
)

reaction(
    label = 'reaction2',
    reactants = ['CO(15)', 'O=CC[C-]=[OH+](1048)'],
    products = ['O=CC(=O)C[C-]=[OH+](1266)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(0.0591985,'m^3/(mol*s)'), n=2.3675, Ea=(305.306,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO;R_H]
Euclidian distance = 0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction3',
    reactants = ['O=CC(=O)C[C-]=[OH+](1266)'],
    products = ['C=C(C=O)[O+]=[C-]O(1311)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(7.50839e+11,'s^-1'), n=0.136197, Ea=(76.4259,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.01703312323488117, var=1.3519057502930525, Tref=1000.0, N=4, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R in family 1,3_sigmatropic_rearrangement."""),
)

reaction(
    label = 'reaction4',
    reactants = ['O=C=COC[C-]=[OH+](1312)'],
    products = ['O=CC(=O)C[C-]=[OH+](1266)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(336.692,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R!H-inRing_N-2R!H->C',), comment="""Estimated from node Root_N-1R!H-inRing_N-2R!H->C in family 1,3_sigmatropic_rearrangement."""),
)

reaction(
    label = 'reaction5',
    reactants = ['C=[C-][OH+]C(=O)C=O(1313)'],
    products = ['O=CC(=O)C[C-]=[OH+](1266)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(7.50839e+11,'s^-1'), n=0.136197, Ea=(108.396,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.01703312323488117, var=1.3519057502930525, Tref=1000.0, N=4, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R in family 1,3_sigmatropic_rearrangement."""),
)

reaction(
    label = 'reaction6',
    reactants = ['O=CC(O)=C[C-]=[OH+](1314)'],
    products = ['O=CC(=O)C[C-]=[OH+](1266)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(5.51726e-18,'s^-1'), n=8.84109, Ea=(152.446,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R!H->C_N-3C-inRing_N-1R!H->N_Ext-3C-R_Ext-5R!H-R_Ext-1C-R_N-Sp-6R!H#5R!H_5R!H->C',), comment="""Estimated from node Root_3R!H->C_N-3C-inRing_N-1R!H->N_Ext-3C-R_Ext-5R!H-R_Ext-1C-R_N-Sp-6R!H#5R!H_5R!H->C in family Ketoenol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['O=C=C(O)C[C-]=[OH+](1267)'],
    products = ['O=CC(=O)C[C-]=[OH+](1266)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(4.7171e-39,'s^-1'), n=14.8922, Ea=(116.486,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.01101064541230846, var=9.816149634195241, Tref=1000.0, N=17, data_mean=0.0, correlation='Root_3R!H->C_N-3C-inRing_N-1R!H->N_Ext-3C-R',), comment="""Estimated from node Root_3R!H->C_N-3C-inRing_N-1R!H->N_Ext-3C-R in family Ketoenol."""),
)

reaction(
    label = 'reaction8',
    reactants = ['O=CC(=O)C=[C-][OH2+](1315)'],
    products = ['O=CC(=O)C[C-]=[OH+](1266)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.94404e-30,'s^-1'), n=12.4085, Ea=(195.776,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.0010171218693333563, var=8.259939515435356, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_3R!H->C_N-3C-inRing_N-1R!H->N_Ext-3C-R_Ext-5R!H-R_6R!H->O',), comment="""Estimated from node Root_3R!H->C_N-3C-inRing_N-1R!H->N_Ext-3C-R_Ext-5R!H-R_6R!H->O in family Ketoenol.
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction9',
    reactants = ['O=CC(=O)C[C-]=[OH+](1266)'],
    products = ['O=C=C=O(1310)', 'C=[C-][OH2+](1309)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.94505e+09,'s^-1'), n=0.601596, Ea=(281.278,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.07591793873856303, var=4.052420941329839, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_4R!H->C_N-1R!H->O_2R!H->C_5R!H->O',), comment="""Estimated from node Root_4R!H->C_N-1R!H->O_2R!H->C_5R!H->O in family Retroene."""),
)

network(
    label = 'PDepNetwork #882',
    isomers = [
        'O=CC(=O)C[C-]=[OH+](1266)',
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
    label = 'PDepNetwork #882',
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

