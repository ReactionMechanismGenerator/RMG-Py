species(
    label = '[O]OC(=CCO)OC(=O)O(1705)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {9,S} {10,S}
2  O u0 p2 c0 {7,S} {14,S}
3  O u0 p2 c0 {6,S} {9,S}
4  O u0 p2 c0 {10,S} {15,S}
5  O u0 p2 c0 {10,D}
6  O u1 p2 c0 {3,S}
7  C u0 p0 c0 {2,S} {8,S} {11,S} {12,S}
8  C u0 p0 c0 {7,S} {9,D} {13,S}
9  C u0 p0 c0 {1,S} {3,S} {8,D}
10 C u0 p0 c0 {1,S} {4,S} {5,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {4,S}
"""),
    E0 = (-601.381,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,350,440,435,1725,180,180,180,180,1596.65,1600,2914.95,3200],'cm^-1')),
        HinderedRotor(inertia=(0.150333,'amu*angstrom^2'), symmetry=1, barrier=(3.45646,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150333,'amu*angstrom^2'), symmetry=1, barrier=(3.45646,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150333,'amu*angstrom^2'), symmetry=1, barrier=(3.45646,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150333,'amu*angstrom^2'), symmetry=1, barrier=(3.45646,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150333,'amu*angstrom^2'), symmetry=1, barrier=(3.45646,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150333,'amu*angstrom^2'), symmetry=1, barrier=(3.45646,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (149.079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.234193,0.103281,-0.000146095,1.01442e-07,-2.36454e-11,-72186.5,35.4467], Tmin=(100,'K'), Tmax=(621.784,'K')), NASAPolynomial(coeffs=[12.826,0.0346426,-1.76096e-05,3.46052e-09,-2.4362e-13,-74107.9,-23.7461], Tmin=(621.784,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-601.381,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsCs) + group(Cds-OdOsOs) + radical(ROOJ)"""),
)

species(
    label = 'CO2(4)',
    structure = adjacencyList("""1 O u0 p2 c0 {3,D}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,D} {2,D}
"""),
    E0 = (-403.131,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([459.169,1086.67,1086.67,2300.05],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0094,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(1622.99,'J/mol'), sigma=(3.941,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.27791,0.00275772,7.12823e-06,-1.07859e-08,4.14249e-12,-48475.6,5.97853], Tmin=(100,'K'), Tmax=(988.175,'K')), NASAPolynomial(coeffs=[4.55068,0.00290734,-1.14646e-06,2.25805e-10,-1.69533e-14,-48986,-1.45644], Tmin=(988.175,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-403.131,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(62.3585,'J/(mol*K)'), label="""CO2""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = '[O]OC(=O)CCO(1660)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {6,S} {12,S}
2  O u0 p2 c0 {4,S} {7,S}
3  O u0 p2 c0 {7,D}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {6,S} {7,S} {8,S} {9,S}
6  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
7  C u0 p0 c0 {2,S} {3,D} {5,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {1,S}
"""),
    E0 = (-358.314,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,492.5,1135,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,180,180,1069.85,1204.02,1600,3200],'cm^-1')),
        HinderedRotor(inertia=(0.152636,'amu*angstrom^2'), symmetry=1, barrier=(3.50941,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152636,'amu*angstrom^2'), symmetry=1, barrier=(3.50941,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152636,'amu*angstrom^2'), symmetry=1, barrier=(3.50941,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152636,'amu*angstrom^2'), symmetry=1, barrier=(3.50941,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.069,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4330.82,'J/mol'), sigma=(5.54449,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=676.46 K, Pc=57.65 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.73685,0.0795945,-0.000132574,1.1403e-07,-3.71033e-11,-42985,25.4208], Tmin=(100,'K'), Tmax=(913.763,'K')), NASAPolynomial(coeffs=[9.03895,0.0241137,-1.00819e-05,1.74099e-09,-1.10552e-13,-43703.3,-9.51139], Tmin=(913.763,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-358.314,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-O2s(Cds-O2d)) + group(O2s-(Os-CdOd)H) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsOs) + radical(C(=O)OOJ)"""),
)

species(
    label = 'C=CC(O)(O[O])OC(=O)O(1725)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {7,S} {10,S}
2  O u0 p2 c0 {7,S} {14,S}
3  O u0 p2 c0 {6,S} {7,S}
4  O u0 p2 c0 {10,S} {15,S}
5  O u0 p2 c0 {10,D}
6  O u1 p2 c0 {3,S}
7  C u0 p0 c0 {1,S} {2,S} {3,S} {8,S}
8  C u0 p0 c0 {7,S} {9,D} {11,S}
9  C u0 p0 c0 {8,D} {12,S} {13,S}
10 C u0 p0 c0 {1,S} {4,S} {5,D}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {4,S}
"""),
    E0 = (-699.233,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (149.079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0165382,0.0916451,-0.000108775,6.78271e-08,-1.70454e-11,-83958.2,35.3511], Tmin=(100,'K'), Tmax=(963.735,'K')), NASAPolynomial(coeffs=[14.5417,0.0313581,-1.49415e-05,2.9169e-09,-2.07183e-13,-86757.9,-34.1876], Tmin=(963.735,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-699.233,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(O2s-OsH) + group(Cs-CsOsOsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-OdOsOs) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(=O)C(CO)C(=O)O(1699)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {8,S} {14,S}
2  O u0 p2 c0 {9,S} {15,S}
3  O u0 p2 c0 {6,S} {10,S}
4  O u0 p2 c0 {9,D}
5  O u0 p2 c0 {10,D}
6  O u1 p2 c0 {3,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {11,S}
8  C u0 p0 c0 {1,S} {7,S} {12,S} {13,S}
9  C u0 p0 c0 {2,S} {4,D} {7,S}
10 C u0 p0 c0 {3,S} {5,D} {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {2,S}
"""),
    E0 = (-736.041,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (149.079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.218789,0.0970234,-0.00013374,9.87672e-08,-2.90028e-11,-88377.2,35.9184], Tmin=(100,'K'), Tmax=(835.538,'K')), NASAPolynomial(coeffs=[13.8415,0.0297158,-1.29125e-05,2.36553e-09,-1.60169e-13,-90726.9,-29.3888], Tmin=(835.538,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-736.041,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-O2s(Cds-O2d)) + group(O2s-(Cds-O2d)H) + group(O2s-(Os-CdOd)H) + group(Cs-CsCsCsH) + group(Cs-CsOsHH) + group(Cds-OdCsOs) + group(Cds-OdCsOs) + radical(C(=O)OOJ)"""),
)

species(
    label = 'O(10)',
    structure = adjacencyList("""multiplicity 3
1 O u2 p2 c0
"""),
    E0 = (243.005,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (15.9994,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1235.53,'J/mol'), sigma=(3.758e-10,'m'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with fixed Lennard Jones Parameters. This is the fallback method! Try improving transport databases!"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,1.49298e-14,-2.05327e-17,9.23927e-21,-1.28064e-24,29226.7,5.11107], Tmin=(100,'K'), Tmax=(3959.16,'K')), NASAPolynomial(coeffs=[2.5,1.19009e-10,-4.23987e-14,6.68965e-18,-3.94352e-22,29226.7,5.11107], Tmin=(3959.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.005,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""O""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = 'O=C(O)OC(=O)[CH]CO(1726)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {8,S} {9,S}
2  O u0 p2 c0 {6,S} {13,S}
3  O u0 p2 c0 {9,S} {14,S}
4  O u0 p2 c0 {8,D}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {2,S} {7,S} {10,S} {11,S}
7  C u1 p0 c0 {6,S} {8,S} {12,S}
8  C u0 p0 c0 {1,S} {4,D} {7,S}
9  C u0 p0 c0 {1,S} {3,S} {5,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {3,S}
"""),
    E0 = (-763.283,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (133.079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.424213,0.0829613,-0.000103214,6.73472e-08,-1.7654e-11,-91676.6,35.7444], Tmin=(100,'K'), Tmax=(926.715,'K')), NASAPolynomial(coeffs=[13.3084,0.0273485,-1.31971e-05,2.58996e-09,-1.84318e-13,-94064.6,-25.4337], Tmin=(926.715,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-763.283,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-O2d)) + group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsOs) + group(Cds-OdOsOs) + radical(CCJCO)"""),
)

species(
    label = 'O=C(O)O[C]1OOC1CO(1727)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {3,S} {7,S}
2  O u0 p2 c0 {9,S} {10,S}
3  O u0 p2 c0 {1,S} {9,S}
4  O u0 p2 c0 {8,S} {14,S}
5  O u0 p2 c0 {10,S} {15,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {11,S}
8  C u0 p0 c0 {4,S} {7,S} {12,S} {13,S}
9  C u1 p0 c0 {2,S} {3,S} {7,S}
10 C u0 p0 c0 {2,S} {5,S} {6,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
"""),
    E0 = (-595.377,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (149.079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.157677,0.0833785,-7.94837e-05,3.76501e-08,-7.16132e-12,-71468.1,31.9555], Tmin=(100,'K'), Tmax=(1252.65,'K')), NASAPolynomial(coeffs=[17.3047,0.0286246,-1.39187e-05,2.75645e-09,-1.97422e-13,-75764,-54.6317], Tmin=(1252.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-595.377,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(Cs-CsCsOsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cds-OdOsOs) + ring(12dioxetane) + radical(Cs_P)"""),
)

species(
    label = 'OCC=C1OOO[C](O)O1(1728)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {9,S} {10,S}
2  O u0 p2 c0 {5,S} {9,S}
3  O u0 p2 c0 {5,S} {10,S}
4  O u0 p2 c0 {7,S} {14,S}
5  O u0 p2 c0 {2,S} {3,S}
6  O u0 p2 c0 {10,S} {15,S}
7  C u0 p0 c0 {4,S} {8,S} {11,S} {12,S}
8  C u0 p0 c0 {7,S} {9,D} {13,S}
9  C u0 p0 c0 {1,S} {2,S} {8,D}
10 C u1 p0 c0 {1,S} {3,S} {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (-748.128,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (149.079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.108607,0.0771194,-4.91488e-05,9.89879e-09,5.31107e-13,-89820.3,34.3593], Tmin=(100,'K'), Tmax=(1291.89,'K')), NASAPolynomial(coeffs=[21.0965,0.0294746,-1.47416e-05,2.93503e-09,-2.09675e-13,-96802.2,-79.1908], Tmin=(1291.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-748.128,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsOs) + group(Cs-OsOsOsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsCs) + ring(Cyclohexanone) + radical(CsOOOring)"""),
)

species(
    label = 'O=C(O)OC1([CH]CO)OO1(1729)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {7,S} {10,S}
2  O u0 p2 c0 {3,S} {7,S}
3  O u0 p2 c0 {2,S} {7,S}
4  O u0 p2 c0 {8,S} {14,S}
5  O u0 p2 c0 {10,S} {15,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {2,S} {3,S} {9,S}
8  C u0 p0 c0 {4,S} {9,S} {11,S} {12,S}
9  C u1 p0 c0 {7,S} {8,S} {13,S}
10 C u0 p0 c0 {1,S} {5,S} {6,D}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
"""),
    E0 = (-614.592,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (149.079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.773144,0.0897917,-9.2336e-05,4.65768e-08,-8.98765e-12,-73733.5,39.4056], Tmin=(100,'K'), Tmax=(1364.56,'K')), NASAPolynomial(coeffs=[22.9695,0.0148464,-4.0739e-06,5.83842e-10,-3.51674e-14,-79715.3,-80.6946], Tmin=(1364.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-614.592,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(Cs-CsOsOsOs) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cds-OdOsOs) + ring(dioxirane) + radical(CCJCOOH)"""),
)

species(
    label = '[O]C1(O)OOC(=CCO)O1(1730)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {7,S} {10,S}
2  O u0 p2 c0 {3,S} {7,S}
3  O u0 p2 c0 {2,S} {10,S}
4  O u0 p2 c0 {8,S} {14,S}
5  O u0 p2 c0 {7,S} {15,S}
6  O u1 p2 c0 {7,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
8  C u0 p0 c0 {4,S} {9,S} {11,S} {12,S}
9  C u0 p0 c0 {8,S} {10,D} {13,S}
10 C u0 p0 c0 {1,S} {3,S} {9,D}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
"""),
    E0 = (-453.194,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (149.079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60596,0.0635602,-2.46899e-05,-7.1598e-08,8.3674e-11,-54431.7,28.8455], Tmin=(100,'K'), Tmax=(460.682,'K')), NASAPolynomial(coeffs=[5.51294,0.04633,-2.29411e-05,4.52788e-09,-3.22222e-13,-54968.8,11.1022], Tmin=(460.682,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-453.194,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsOsOs) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsCs) + ring(Cyclopentane) + radical(OCOJ)"""),
)

species(
    label = 'O2(2)',
    structure = adjacencyList("""multiplicity 3
1 O u1 p2 c0 {2,S}
2 O u1 p2 c0 {1,S}
"""),
    E0 = (-8.62178,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1483.7],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (31.9988,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(887.157,'J/mol'), sigma=(3.467,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53763,-0.00122825,5.36751e-06,-4.93119e-09,1.45951e-12,-1037.99,4.67181], Tmin=(100,'K'), Tmax=(1087.72,'K')), NASAPolynomial(coeffs=[3.16428,0.00169452,-8.00324e-07,1.59027e-10,-1.14889e-14,-1048.45,6.08296], Tmin=(1087.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.62178,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""O2""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = 'O=C(O)O[C]=CCO(1731)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {7,S} {8,S}
2  O u0 p2 c0 {5,S} {12,S}
3  O u0 p2 c0 {7,S} {13,S}
4  O u0 p2 c0 {7,D}
5  C u0 p0 c0 {2,S} {6,S} {9,S} {10,S}
6  C u0 p0 c0 {5,S} {8,D} {11,S}
7  C u0 p0 c0 {1,S} {3,S} {4,D}
8  C u1 p0 c0 {1,S} {6,D}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
"""),
    E0 = (-481,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,1685,370,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.624849,0.0715127,-7.3914e-05,3.80459e-08,-7.73887e-12,-57727.1,30.8569], Tmin=(100,'K'), Tmax=(1191.89,'K')), NASAPolynomial(coeffs=[15.9936,0.019935,-9.00307e-06,1.73891e-09,-1.23438e-13,-61390.6,-45.9862], Tmin=(1191.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-481,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdOsOs) + radical(C=CJO)"""),
)

species(
    label = 'O=C(O)O(913)',
    structure = adjacencyList("""1 O u0 p2 c0 {4,S} {5,S}
2 O u0 p2 c0 {4,S} {6,S}
3 O u0 p2 c0 {4,D}
4 C u0 p0 c0 {1,S} {2,S} {3,D}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
"""),
    E0 = (-625.106,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (62.0248,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.98121,0.0114384,2.77965e-05,-4.94724e-08,2.0801e-11,-75136.3,10.9456], Tmin=(100,'K'), Tmax=(956.744,'K')), NASAPolynomial(coeffs=[11.8952,0.00171628,-1.48077e-07,9.26264e-11,-1.38642e-14,-78102.7,-38.2537], Tmin=(956.744,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-625.106,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(124.717,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(O2s-(Cds-O2d)H) + group(Cds-OdOsOs)"""),
)

species(
    label = '[O]OC#CCO(1732)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {4,S} {9,S}
2 O u0 p2 c0 {3,S} {6,S}
3 O u1 p2 c0 {2,S}
4 C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
5 C u0 p0 c0 {4,S} {6,T}
6 C u0 p0 c0 {2,S} {5,T}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {1,S}
"""),
    E0 = (146.056,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (87.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.86962,0.049814,-7.19766e-05,5.71342e-08,-1.76802e-11,17640.4,22.4511], Tmin=(100,'K'), Tmax=(915.064,'K')), NASAPolynomial(coeffs=[7.7323,0.0171027,-6.74295e-06,1.14831e-09,-7.32575e-14,16864,-3.69192], Tmin=(915.064,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(146.056,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-OsCt) + group(O2s-OsH) + group(Cs-CtOsHH) + group(Ct-CtCs) + group(Ct-CtOs) + radical(ROOJ)"""),
)

species(
    label = 'O=C(O)OC(=[C]CO)OO(1733)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {8,S} {9,S}
2  O u0 p2 c0 {5,S} {8,S}
3  O u0 p2 c0 {7,S} {13,S}
4  O u0 p2 c0 {9,S} {14,S}
5  O u0 p2 c0 {2,S} {15,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {3,S} {10,S} {11,S} {12,S}
8  C u0 p0 c0 {1,S} {2,S} {10,D}
9  C u0 p0 c0 {1,S} {4,S} {6,D}
10 C u1 p0 c0 {7,S} {8,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
"""),
    E0 = (-515.544,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,3580,3650,1210,1345,900,1100,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,1685,370,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (149.079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.853916,0.116136,-0.00018763,1.59654e-07,-5.39077e-11,-61839.8,37.3372], Tmin=(100,'K'), Tmax=(781.379,'K')), NASAPolynomial(coeffs=[13.571,0.0351579,-1.84814e-05,3.65201e-09,-2.56639e-13,-63876.3,-27.3025], Tmin=(781.379,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-515.544,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsCs) + group(Cds-OdOsOs) + radical(Cds_S)"""),
)

species(
    label = 'OH(9)',
    structure = adjacencyList("""multiplicity 2
1 O u1 p2 c0 {2,S}
2 H u0 p0 c0 {1,S}
"""),
    E0 = (28.372,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3287.46],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (17.0074,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(1235.53,'J/mol'), sigma=(3.758e-10,'m'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with fixed Lennard Jones Parameters. This is the fallback method! Try improving transport databases!"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.4858,0.00133395,-4.70038e-06,5.64372e-09,-2.06315e-12,3411.96,1.99787], Tmin=(100,'K'), Tmax=(1005.25,'K')), NASAPolynomial(coeffs=[2.88225,0.0010387,-2.35657e-07,1.40242e-11,6.34477e-16,3669.56,5.59056], Tmin=(1005.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.372,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""OH""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = 'O=C(O)OC(=O)C=CO(1734)',
    structure = adjacencyList("""1  O u0 p2 c0 {7,S} {9,S}
2  O u0 p2 c0 {8,S} {12,S}
3  O u0 p2 c0 {9,S} {13,S}
4  O u0 p2 c0 {7,D}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {7,S} {8,D} {10,S}
7  C u0 p0 c0 {1,S} {4,D} {6,S}
8  C u0 p0 c0 {2,S} {6,D} {11,S}
9  C u0 p0 c0 {1,S} {3,S} {5,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
"""),
    E0 = (-830.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.204049,0.0766054,-9.4966e-05,5.16562e-08,-9.9114e-12,-99770.9,30.648], Tmin=(100,'K'), Tmax=(985.47,'K')), NASAPolynomial(coeffs=[21.6122,0.00195007,3.73975e-08,-9.7691e-12,-1.81794e-15,-104585,-75.3361], Tmin=(985.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-830.73,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-O2d)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-O2d)H) + missing(Cd-COCdH0) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-OdOsOs)"""),
)

species(
    label = '[O]C(=O)OC(=CCO)OO(1735)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {9,S} {10,S}
2  O u0 p2 c0 {4,S} {9,S}
3  O u0 p2 c0 {7,S} {14,S}
4  O u0 p2 c0 {2,S} {15,S}
5  O u1 p2 c0 {10,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {3,S} {8,S} {11,S} {12,S}
8  C u0 p0 c0 {7,S} {9,D} {13,S}
9  C u0 p0 c0 {1,S} {2,S} {8,D}
10 C u0 p0 c0 {1,S} {5,S} {6,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
"""),
    E0 = (-509.638,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,350,440,435,1725,180,180,180,390.961,763.349,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155131,'amu*angstrom^2'), symmetry=1, barrier=(3.56678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155131,'amu*angstrom^2'), symmetry=1, barrier=(3.56678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155131,'amu*angstrom^2'), symmetry=1, barrier=(3.56678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155131,'amu*angstrom^2'), symmetry=1, barrier=(3.56678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155131,'amu*angstrom^2'), symmetry=1, barrier=(3.56678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155131,'amu*angstrom^2'), symmetry=1, barrier=(3.56678,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (149.079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.547924,0.109564,-0.000176704,1.5308e-07,-5.25463e-11,-61140.6,36.2402], Tmin=(100,'K'), Tmax=(799.622,'K')), NASAPolynomial(coeffs=[11.6865,0.0369307,-1.90067e-05,3.72308e-09,-2.60204e-13,-62731.7,-17.7627], Tmin=(799.622,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-509.638,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsCs) + group(Cds-OdOsOs) + radical(OC=OOJ)"""),
)

species(
    label = '[O]CC=C(OO)OC(=O)O(1736)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {9,S} {10,S}
2  O u0 p2 c0 {4,S} {9,S}
3  O u0 p2 c0 {10,S} {14,S}
4  O u0 p2 c0 {2,S} {15,S}
5  O u1 p2 c0 {7,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {5,S} {8,S} {11,S} {12,S}
8  C u0 p0 c0 {7,S} {9,D} {13,S}
9  C u0 p0 c0 {1,S} {2,S} {8,D}
10 C u0 p0 c0 {1,S} {3,S} {6,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
"""),
    E0 = (-527.681,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,350,440,435,1725,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (149.079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.478663,0.109076,-0.000176444,1.57772e-07,-5.63963e-11,-63314.3,36.3742], Tmin=(100,'K'), Tmax=(780.582,'K')), NASAPolynomial(coeffs=[9.97014,0.0422461,-2.24886e-05,4.47824e-09,-3.1644e-13,-64540.7,-8.85414], Tmin=(780.582,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-527.681,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsCs) + group(Cds-OdOsOs) + radical(CCOJ)"""),
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
    E0 = (-128.36,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (72.1047,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-52.3647,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-47.2571,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-100.078,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-120.581,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-99.0094,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (26.3443,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-16.6006,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (39.5652,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (1.78541,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-63.6337,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (32.9261,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (17.9409,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]OC(=CCO)OC(=O)O(1705)'],
    products = ['CO2(4)', '[O]OC(=O)CCO(1660)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(2.46545e+14,'s^-1'), n=0.20628, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(3000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R!H->C',), comment="""Estimated from node Root_N-4R!H->C in family Retroene."""),
)

reaction(
    label = 'reaction2',
    reactants = ['[O]OC(=CCO)OC(=O)O(1705)'],
    products = ['C=CC(O)(O[O])OC(=O)O(1725)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.30791e+26,'s^-1'), n=-3.67984, Ea=(200.465,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.4163866947181062, var=128.1927359939506, Tref=1000.0, N=20, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root in family 1,3_sigmatropic_rearrangement."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[O]OC(=CCO)OC(=O)O(1705)'],
    products = ['[O]OC(=O)C(CO)C(=O)O(1699)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(7.50839e+11,'s^-1'), n=0.136197, Ea=(75.9955,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.01703312323488117, var=1.3519057502930525, Tref=1000.0, N=4, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R in family 1,3_sigmatropic_rearrangement."""),
)

reaction(
    label = 'reaction4',
    reactants = ['O(10)', 'O=C(O)OC(=O)[CH]CO(1726)'],
    products = ['[O]OC(=CCO)OC(=O)O(1705)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(15399.6,'m^3/(mol*s)'), n=0.932885, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [O_sec_rad;O_birad] for rate rule [O_rad/OneDe;O_birad]
Euclidian distance = 1.0
family: Birad_R_Recombination
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]OC(=CCO)OC(=O)O(1705)'],
    products = ['O=C(O)O[C]1OOC1CO(1727)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.391e+11,'s^-1'), n=0.287, Ea=(28.2817,'kJ/mol'), T0=(1,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Backbone1_N-2R!H-inRing_Ext-4R!H-R_N-2R!H->C',), comment="""Estimated from node Backbone1_N-2R!H-inRing_Ext-4R!H-R_N-2R!H->C in family Intra_R_Add_Endocyclic."""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]OC(=CCO)OC(=O)O(1705)'],
    products = ['OCC=C1OOO[C](O)O1(1728)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(132.911,'s^-1'), n=2.22721, Ea=(7.77918,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-1.425017045770312, var=4.564120841035978, Tref=1000.0, N=3, data_mean=0.0, correlation='Backbone3_N-Sp-4R!H=1R!H_Sp-2R!H-1R!H_Ext-3R!H-R',), comment="""Estimated from node Backbone3_N-Sp-4R!H=1R!H_Sp-2R!H-1R!H_Ext-3R!H-R in family Intra_R_Add_Endocyclic."""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]OC(=CCO)OC(=O)O(1705)'],
    products = ['O=C(O)OC1([CH]CO)OO1(1729)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(3.36539e+09,'s^-1'), n=0.84129, Ea=(29.3508,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.031430614527560546, var=6.099077214950256, Tref=1000.0, N=84, data_mean=0.0, correlation='Backbone1',), comment="""Estimated from node Backbone1 in family Intra_R_Add_Exocyclic."""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]OC(=CCO)OC(=O)O(1705)'],
    products = ['[O]C1(O)OOC(=CCO)O1(1730)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.234039,'s^-1'), n=3.38729, Ea=(154.704,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.23174178600410858, var=45.068923389250436, Tref=1000.0, N=140, data_mean=0.0, correlation='Backbone3',), comment="""Estimated from node Backbone3 in family Intra_R_Add_Exocyclic."""),
)

reaction(
    label = 'reaction9',
    reactants = ['O2(2)', 'O=C(O)O[C]=CCO(1731)'],
    products = ['[O]OC(=CCO)OC(=O)O(1705)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.21695e+08,'m^3/(mol*s)'), n=-0.428622, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing',), comment="""Estimated from node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing in family R_Recombination.
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]OC(=CCO)OC(=O)O(1705)'],
    products = ['O=C(O)O(913)', '[O]OC#CCO(1732)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(6.63512e+11,'s^-1'), n=0, Ea=(167.925,'kJ/mol'), T0=(1,'K'), Tmin=(500,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_4R!H->C_1R!H->O_Ext-2R!H-R_7R!H->O',), comment="""Estimated from node Root_4R!H->C_1R!H->O_Ext-2R!H-R_7R!H->O in family Retroene."""),
)

reaction(
    label = 'reaction11',
    reactants = ['O=C(O)OC(=[C]CO)OO(1733)'],
    products = ['[O]OC(=CCO)OC(=O)O(1705)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;XH_out] for rate rule [R4H_DSS;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 2.23606797749979
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]OC(=CCO)OC(=O)O(1705)'],
    products = ['OH(9)', 'O=C(O)OC(=O)C=CO(1734)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(23000,'s^-1'), n=2.11, Ea=(64.7265,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5H;O_rad_out;Cs_H_out_H/NonDeO] for rate rule [R5H_SSMS;O_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]C(=O)OC(=CCO)OO(1735)'],
    products = ['[O]OC(=CCO)OC(=O)O(1705)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(247806,'s^-1'), n=1.46258, Ea=(69.5427,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_SSSSS;O_rad_out;XH_out] for rate rule [R6H_SSSSS;O_rad_out;O_H_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O]CC=C(OO)OC(=O)O(1736)'],
    products = ['[O]OC(=CCO)OC(=O)O(1705)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.33753e+06,'s^-1'), n=1.02312, Ea=(72.6006,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H;O_rad_out;XH_out] for rate rule [R6H_RSMSR;O_rad_out;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #1181',
    isomers = [
        '[O]OC(=CCO)OC(=O)O(1705)',
    ],
    reactants = [
        ('CO2(4)', '[O]OC(=O)CCO(1660)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1181',
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

