species(
    label = '[O]OCC(=O)C(=O)C(=O)O(1618)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {6,S} {7,S}
2  O u0 p2 c0 {10,S} {13,S}
3  O u0 p2 c0 {8,D}
4  O u0 p2 c0 {9,D}
5  O u0 p2 c0 {10,D}
6  O u1 p2 c0 {1,S}
7  C u0 p0 c0 {1,S} {8,S} {11,S} {12,S}
8  C u0 p0 c0 {3,D} {7,S} {9,S}
9  C u0 p0 c0 {4,D} {8,S} {10,S}
10 C u0 p0 c0 {2,S} {5,D} {9,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {2,S}
"""),
    E0 = (-561.206,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,365,385,505,600,445,480,1700,1720,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (147.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.425322,0.10837,-0.000187625,1.67021e-07,-5.74593e-11,-67348.8,33.824], Tmin=(100,'K'), Tmax=(836.562,'K')), NASAPolynomial(coeffs=[11.3232,0.0327332,-1.71078e-05,3.32478e-09,-2.29551e-13,-68633.4,-16.6891], Tmin=(836.562,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-561.206,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-O2d)H) + group(O2s-OsH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)(Cds-O2d)) + group(Cds-O2d(Cds-O2d)O2s) + radical(ROOJ)"""),
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
    label = '[O]OCC(O)=C=O(1445)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {6,S} {10,S}
3  O u1 p2 c0 {1,S}
4  O u0 p2 c0 {7,D}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
6  C u0 p0 c0 {2,S} {5,S} {7,D}
7  C u0 p0 c0 {4,D} {6,D}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {2,S}
"""),
    E0 = (-174.621,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2120,512.5,787.5,180],'cm^-1')),
        HinderedRotor(inertia=(0.802099,'amu*angstrom^2'), symmetry=1, barrier=(18.4418,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.800815,'amu*angstrom^2'), symmetry=1, barrier=(18.4123,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.798681,'amu*angstrom^2'), symmetry=1, barrier=(18.3633,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.053,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4036.53,'J/mol'), sigma=(4.8121,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=630.50 K, Pc=82.2 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.395357,0.0813886,-0.000128624,9.68091e-08,-2.74456e-11,-20874.2,23.1948], Tmin=(100,'K'), Tmax=(971.254,'K')), NASAPolynomial(coeffs=[15.5342,0.00949717,-2.85599e-06,3.6488e-10,-1.68126e-14,-23364.8,-47.0823], Tmin=(971.254,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-174.621,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)OsHH) + group(Cds-(Cdd-O2d)CsOs) + missing(Cdd-CdO2d) + radical(ROOJ)"""),
)

species(
    label = 'C=C(OO[O])C(=O)C(=O)O(1645)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {3,S} {7,S}
2  O u0 p2 c0 {9,S} {13,S}
3  O u0 p2 c0 {1,S} {6,S}
4  O u0 p2 c0 {8,D}
5  O u0 p2 c0 {9,D}
6  O u1 p2 c0 {3,S}
7  C u0 p0 c0 {1,S} {8,S} {10,D}
8  C u0 p0 c0 {4,D} {7,S} {9,S}
9  C u0 p0 c0 {2,S} {5,D} {8,S}
10 C u0 p0 c0 {7,D} {11,S} {12,S}
11 H u0 p0 c0 {10,S}
12 H u0 p0 c0 {10,S}
13 H u0 p0 c0 {2,S}
"""),
    E0 = (-335.707,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,375,552.5,462.5,1710,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (147.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.334391,0.0862945,-0.000112263,7.66621e-08,-2.10876e-11,-40249.2,33.9736], Tmin=(100,'K'), Tmax=(883.26,'K')), NASAPolynomial(coeffs=[13.0429,0.0287416,-1.45229e-05,2.88947e-09,-2.0671e-13,-42494.1,-25.7598], Tmin=(883.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-335.707,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-O2d)H) + group(O2s-OsOs) + group(O2s-OsH) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)O2s) + group(Cds-CdsHH) + radical(ROOJ)"""),
)

species(
    label = '[O]OCC(=C=O)OC(=O)O(1646)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {8,S} {9,S}
2  O u0 p2 c0 {5,S} {7,S}
3  O u0 p2 c0 {9,S} {13,S}
4  O u0 p2 c0 {9,D}
5  O u1 p2 c0 {2,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {2,S} {8,S} {11,S} {12,S}
8  C u0 p0 c0 {1,S} {7,S} {10,D}
9  C u0 p0 c0 {1,S} {3,S} {4,D}
10 C u0 p0 c0 {6,D} {8,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {3,S}
"""),
    E0 = (-539.561,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2120,512.5,787.5,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (147.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.828602,0.111435,-0.000169808,1.26904e-07,-3.68109e-11,-64725,32.4306], Tmin=(100,'K'), Tmax=(851.752,'K')), NASAPolynomial(coeffs=[18.4746,0.020783,-1.01632e-05,1.94958e-09,-1.35103e-13,-68013.3,-57.5988], Tmin=(851.752,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-539.561,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-OsCs) + group(O2s-(Cds-O2d)H) + group(O2s-OsH) + missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)OsHH) + group(Cds-(Cdd-O2d)CsOs) + group(Cds-OdOsOs) + missing(Cdd-CdO2d) + radical(ROOJ)"""),
)

species(
    label = '[O]OCOC(=C=O)C(=O)O(1647)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {7,S} {8,S}
2  O u0 p2 c0 {5,S} {7,S}
3  O u0 p2 c0 {9,S} {13,S}
4  O u0 p2 c0 {9,D}
5  O u1 p2 c0 {2,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {2,S} {11,S} {12,S}
8  C u0 p0 c0 {1,S} {9,S} {10,D}
9  C u0 p0 c0 {3,S} {4,D} {8,S}
10 C u0 p0 c0 {6,D} {8,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {3,S}
"""),
    E0 = (-445.424,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2120,512.5,787.5,180,180,180,180,1600,1679.39,2815.05,3200],'cm^-1')),
        HinderedRotor(inertia=(0.152446,'amu*angstrom^2'), symmetry=1, barrier=(3.50504,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152446,'amu*angstrom^2'), symmetry=1, barrier=(3.50504,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152446,'amu*angstrom^2'), symmetry=1, barrier=(3.50504,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152446,'amu*angstrom^2'), symmetry=1, barrier=(3.50504,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152446,'amu*angstrom^2'), symmetry=1, barrier=(3.50504,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (147.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.991322,0.101207,-0.000122903,6.31955e-08,-1.03558e-11,-53384.5,35.1068], Tmin=(100,'K'), Tmax=(933.474,'K')), NASAPolynomial(coeffs=[26.8951,0.00315603,1.95774e-07,-1.0859e-10,7.06783e-15,-59525,-102.514], Tmin=(933.474,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-445.424,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-(Cds-O2d)H) + group(O2s-OsH) + missing(O2d-Cdd) + group(Cs-OsOsHH) + group(Cds-(Cdd-O2d)CsOs) + group(Cds-O2d(Cds-Cds)O2s) + missing(Cdd-CdO2d) + radical(ROOJ)"""),
)

species(
    label = '[O]OCC(=O)C(=C=O)OO(1648)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {3,S} {9,S}
2  O u0 p2 c0 {5,S} {7,S}
3  O u0 p2 c0 {1,S} {13,S}
4  O u0 p2 c0 {8,D}
5  O u1 p2 c0 {2,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {2,S} {8,S} {11,S} {12,S}
8  C u0 p0 c0 {4,D} {7,S} {9,S}
9  C u0 p0 c0 {1,S} {8,S} {10,D}
10 C u0 p0 c0 {6,D} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {3,S}
"""),
    E0 = (-165.559,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,350,440,435,1725,2120,512.5,787.5,180,181.741],'cm^-1')),
        HinderedRotor(inertia=(0.373906,'amu*angstrom^2'), symmetry=1, barrier=(8.59684,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.373481,'amu*angstrom^2'), symmetry=1, barrier=(8.60079,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00513219,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.374163,'amu*angstrom^2'), symmetry=1, barrier=(8.6165,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.41948,'amu*angstrom^2'), symmetry=1, barrier=(32.6889,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (147.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.464996,0.109486,-0.000191095,1.70741e-07,-5.89357e-11,-19762.2,35.1175], Tmin=(100,'K'), Tmax=(833.289,'K')), NASAPolynomial(coeffs=[11.4652,0.0325693,-1.72685e-05,3.37566e-09,-2.33824e-13,-21068.3,-16.1698], Tmin=(833.289,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-165.559,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(O2s-OsH) + missing(O2d-Cdd) + group(Cs-(Cds-O2d)OsHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-(Cdd-O2d)CsOs) + missing(Cdd-CdO2d) + radical(ROOJ)"""),
)

species(
    label = '[O]OCC(=O)OC(O)=C=O(1609)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {8,S} {9,S}
2  O u0 p2 c0 {5,S} {7,S}
3  O u0 p2 c0 {9,S} {13,S}
4  O u0 p2 c0 {8,D}
5  O u1 p2 c0 {2,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {2,S} {8,S} {11,S} {12,S}
8  C u0 p0 c0 {1,S} {4,D} {7,S}
9  C u0 p0 c0 {1,S} {3,S} {10,D}
10 C u0 p0 c0 {6,D} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {3,S}
"""),
    E0 = (-488.678,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2120,512.5,787.5,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (147.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.16639,0.109254,-0.00013117,6.96341e-08,-1.34783e-11,-58529.7,36.6216], Tmin=(100,'K'), Tmax=(1467.84,'K')), NASAPolynomial(coeffs=[32.7391,-0.0049538,5.04519e-06,-1.09168e-09,7.64443e-14,-66720.6,-138.169], Tmin=(1467.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-488.678,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + missing(O2d-Cdd) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsOs) + group(Cds-(Cdd-O2d)OsOs) + missing(Cdd-CdO2d) + radical(ROOJ)"""),
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
    label = '[O]CC(=O)C(=O)C(=O)O(1649)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {9,S} {12,S}
2  O u1 p2 c0 {6,S}
3  O u0 p2 c0 {7,D}
4  O u0 p2 c0 {8,D}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {2,S} {7,S} {10,S} {11,S}
7  C u0 p0 c0 {3,D} {6,S} {8,S}
8  C u0 p0 c0 {4,D} {7,S} {9,S}
9  C u0 p0 c0 {1,S} {5,D} {8,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {1,S}
"""),
    E0 = (-540.983,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,365,385,505,600,445,480,1700,1720,180,180,180,202.486,1566.33,1600,2933.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.146309,'amu*angstrom^2'), symmetry=1, barrier=(3.36393,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.146309,'amu*angstrom^2'), symmetry=1, barrier=(3.36393,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.146309,'amu*angstrom^2'), symmetry=1, barrier=(3.36393,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.146309,'amu*angstrom^2'), symmetry=1, barrier=(3.36393,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (131.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.218182,0.0922116,-0.000155888,1.37802e-07,-4.73544e-11,-64937.7,30.7742], Tmin=(100,'K'), Tmax=(832.216,'K')), NASAPolynomial(coeffs=[10.1073,0.0293854,-1.5081e-05,2.92175e-09,-2.01806e-13,-66054,-11.9365], Tmin=(832.216,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-540.983,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)(Cds-O2d)) + group(Cds-O2d(Cds-O2d)O2s) + radical(C=OCOJ)"""),
)

species(
    label = 'O=C(O)C(=O)[C]1COOO1(1650)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {3,S} {7,S}
2  O u0 p2 c0 {3,S} {8,S}
3  O u0 p2 c0 {1,S} {2,S}
4  O u0 p2 c0 {10,S} {13,S}
5  O u0 p2 c0 {9,D}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {8,S} {11,S} {12,S}
8  C u1 p0 c0 {2,S} {7,S} {9,S}
9  C u0 p0 c0 {5,D} {8,S} {10,S}
10 C u0 p0 c0 {4,S} {6,D} {9,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {4,S}
"""),
    E0 = (-346.411,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (147.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.00570539,0.0769065,-7.5923e-05,3.81158e-08,-7.45267e-12,-41510.5,31.9763], Tmin=(100,'K'), Tmax=(1296.72,'K')), NASAPolynomial(coeffs=[18.4199,0.0181767,-5.75691e-06,8.95921e-10,-5.59104e-14,-46124.1,-61.0214], Tmin=(1296.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-346.411,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)O2s) + ring(123trioxolane) + radical(C2CsJOO)"""),
)

species(
    label = 'O=C(O)[C]1OOOCC1=O(1651)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {3,S} {7,S}
2  O u0 p2 c0 {3,S} {8,S}
3  O u0 p2 c0 {1,S} {2,S}
4  O u0 p2 c0 {10,S} {13,S}
5  O u0 p2 c0 {9,D}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {9,S} {11,S} {12,S}
8  C u1 p0 c0 {2,S} {9,S} {10,S}
9  C u0 p0 c0 {5,D} {7,S} {8,S}
10 C u0 p0 c0 {4,S} {6,D} {8,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {4,S}
"""),
    E0 = (-401.973,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (147.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.956629,0.0485486,1.23573e-05,-4.98004e-08,2.06171e-11,-48220.9,29.2567], Tmin=(100,'K'), Tmax=(1074.24,'K')), NASAPolynomial(coeffs=[19.0042,0.0239671,-1.28305e-05,2.76376e-09,-2.10886e-13,-54557.5,-70.5505], Tmin=(1074.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-401.973,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(O2s-(Cds-O2d)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsCs) + group(Cds-OdCsOs) + ring(Cyclohexanone) + radical(C2CsJOO)"""),
)

species(
    label = 'O=C1COOO[C](O)C1=O(1652)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {3,S} {7,S}
2  O u0 p2 c0 {3,S} {10,S}
3  O u0 p2 c0 {1,S} {2,S}
4  O u0 p2 c0 {10,S} {13,S}
5  O u0 p2 c0 {8,D}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {8,S} {11,S} {12,S}
8  C u0 p0 c0 {5,D} {7,S} {9,S}
9  C u0 p0 c0 {6,D} {8,S} {10,S}
10 C u1 p0 c0 {2,S} {4,S} {9,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {4,S}
"""),
    E0 = (-291.555,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (147.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.530521,0.0242385,0.000156086,-2.67282e-07,1.17657e-10,-34892.3,34.2851], Tmin=(100,'K'), Tmax=(910.726,'K')), NASAPolynomial(coeffs=[44.0165,-0.0303323,2.12706e-05,-4.11486e-09,2.64926e-13,-48470.7,-202.504], Tmin=(910.726,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-291.555,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsOs) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)Cs) + ring(Cycloheptane) + radical(Cs_P)"""),
)

species(
    label = '[O]C1(C(=O)C(=O)O)COO1(1653)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {2,S} {7,S}
2  O u0 p2 c0 {1,S} {8,S}
3  O u0 p2 c0 {10,S} {13,S}
4  O u1 p2 c0 {7,S}
5  O u0 p2 c0 {9,D}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
8  C u0 p0 c0 {2,S} {7,S} {11,S} {12,S}
9  C u0 p0 c0 {5,D} {7,S} {10,S}
10 C u0 p0 c0 {3,S} {6,D} {9,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {3,S}
"""),
    E0 = (-468.636,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (147.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.53273,0.0885861,-9.52518e-05,4.84105e-08,-9.43293e-12,-56191.3,30.0468], Tmin=(100,'K'), Tmax=(1266.76,'K')), NASAPolynomial(coeffs=[23.5822,0.0124394,-5.08505e-06,9.57887e-10,-6.80013e-14,-62300.9,-91.9962], Tmin=(1266.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-468.636,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(Cs-CsCsOsOs) + group(Cs-CsOsHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)O2s) + ring(12dioxetane) + radical(C=OCOJ)"""),
)

species(
    label = '[O]C1(C(=O)O)OOCC1=O(1654)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {2,S} {7,S}
2  O u0 p2 c0 {1,S} {8,S}
3  O u0 p2 c0 {10,S} {13,S}
4  O u1 p2 c0 {7,S}
5  O u0 p2 c0 {9,D}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {9,S} {11,S} {12,S}
9  C u0 p0 c0 {5,D} {7,S} {8,S}
10 C u0 p0 c0 {3,S} {6,D} {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {3,S}
"""),
    E0 = (-571.577,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (147.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.461008,0.0503373,3.61076e-05,-1.02394e-07,4.75812e-11,-68592.3,31.1309], Tmin=(100,'K'), Tmax=(944.663,'K')), NASAPolynomial(coeffs=[27.9499,0.00315537,1.12262e-06,-1.43696e-10,-4.99493e-15,-76874.1,-116.268], Tmin=(944.663,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-571.577,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(Cs-CsCsOsOs) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsCs) + group(Cds-OdCsOs) + ring(Cyclopentane) + radical(C=OCOJ)"""),
)

species(
    label = '[O]C1(O)OOCC(=O)C1=O(1655)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {2,S} {7,S}
2  O u0 p2 c0 {1,S} {8,S}
3  O u0 p2 c0 {7,S} {13,S}
4  O u1 p2 c0 {7,S}
5  O u0 p2 c0 {10,D}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {3,S} {4,S} {9,S}
8  C u0 p0 c0 {2,S} {10,S} {11,S} {12,S}
9  C u0 p0 c0 {6,D} {7,S} {10,S}
10 C u0 p0 c0 {5,D} {8,S} {9,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {3,S}
"""),
    E0 = (-470.842,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (147.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.12921,0.163194,-0.000343155,3.63977e-07,-1.46847e-10,-56436.4,20.7077], Tmin=(100,'K'), Tmax=(775.859,'K')), NASAPolynomial(coeffs=[4.32533,0.0699971,-4.71265e-05,1.00673e-08,-7.3412e-13,-55634.5,2.82879], Tmin=(775.859,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-470.842,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsOs) + group(Cs-(Cds-O2d)OsHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)Cs) + ring(six-sidedoubles) + radical(C=OCOJ)"""),
)

species(
    label = '[O]OC=C(O)C(=O)C(=O)O(1619)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {7,S} {13,S}
2  O u0 p2 c0 {10,S} {12,S}
3  O u0 p2 c0 {6,S} {9,S}
4  O u0 p2 c0 {8,D}
5  O u0 p2 c0 {10,D}
6  O u1 p2 c0 {3,S}
7  C u0 p0 c0 {1,S} {8,S} {9,D}
8  C u0 p0 c0 {4,D} {7,S} {10,S}
9  C u0 p0 c0 {3,S} {7,D} {11,S}
10 C u0 p0 c0 {2,S} {5,D} {8,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {1,S}
"""),
    E0 = (-580.477,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,492.5,1135,1000,350,440,435,1725,375,552.5,462.5,1710,3010,987.5,1337.5,450,1655,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (147.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.509184,0.100631,-0.000136041,8.9272e-08,-2.26873e-11,-69654,32.7356], Tmin=(100,'K'), Tmax=(970.738,'K')), NASAPolynomial(coeffs=[19.536,0.0180341,-8.41188e-06,1.62175e-09,-1.1443e-13,-73545.8,-63.3759], Tmin=(970.738,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-580.477,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-O2d)H) + group(O2s-OsH) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-O2d)O2s) + radical(ROOJ)"""),
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
    label = '[CH2]C(=O)C(=O)C(=O)O(1468)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {7,S} {11,S}
2  O u0 p2 c0 {6,D}
3  O u0 p2 c0 {5,D}
4  O u0 p2 c0 {7,D}
5  C u0 p0 c0 {3,D} {6,S} {7,S}
6  C u0 p0 c0 {2,D} {5,S} {8,S}
7  C u0 p0 c0 {1,S} {4,D} {5,S}
8  C u1 p0 c0 {6,S} {9,S} {10,S}
9  H u0 p0 c0 {8,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {1,S}
"""),
    E0 = (-427.873,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,365,385,505,600,445,480,1700,1720,3000,3100,440,815,1455,1000,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4849.38,'J/mol'), sigma=(5.95975,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=757.46 K, Pc=51.98 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.689781,0.0800539,-0.000131579,1.14811e-07,-3.93808e-11,-51349,27.6622], Tmin=(100,'K'), Tmax=(814.501,'K')), NASAPolynomial(coeffs=[9.55244,0.02623,-1.34885e-05,2.62941e-09,-1.82845e-13,-52451.1,-11.1793], Tmin=(814.501,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-427.873,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)HHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)(Cds-O2d)) + group(Cds-O2d(Cds-O2d)O2s) + radical(CJCC=O)"""),
)

species(
    label = '[O]OC=C=O(111)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {2,S} {4,S}
2 O u1 p2 c0 {1,S}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {1,S} {5,D} {6,S}
5 C u0 p0 c0 {3,D} {4,D}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (63.7216,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (73.0274,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.46746,0.0359239,-5.52506e-05,4.42272e-08,-1.39769e-11,7717.02,15.4679], Tmin=(100,'K'), Tmax=(798.729,'K')), NASAPolynomial(coeffs=[7.55928,0.00970831,-4.67347e-06,8.90219e-10,-6.13029e-14,6926.46,-7.80988], Tmin=(798.729,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(63.7216,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""ketenylperoxy""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'O=C=C(O)O(1297)',
    structure = adjacencyList("""1 O u0 p2 c0 {4,S} {6,S}
2 O u0 p2 c0 {4,S} {7,S}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {1,S} {2,S} {5,D}
5 C u0 p0 c0 {3,D} {4,D}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {2,S}
"""),
    E0 = (-383.509,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,350,440,435,1725,2120,512.5,787.5],'cm^-1')),
        HinderedRotor(inertia=(1.3496,'amu*angstrom^2'), symmetry=1, barrier=(31.03,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.35264,'amu*angstrom^2'), symmetry=1, barrier=(31.0999,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (74.0354,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3234,0.0319077,4.00991e-05,-1.14681e-07,5.8394e-11,-46003.4,13.6979], Tmin=(100,'K'), Tmax=(885.499,'K')), NASAPolynomial(coeffs=[31.5632,-0.0301777,1.90444e-05,-3.79828e-09,2.59031e-13,-54280.2,-145.01], Tmin=(885.499,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-383.509,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + missing(O2d-Cdd) + group(Cds-(Cdd-O2d)OsOs) + missing(Cdd-CdO2d)"""),
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
    label = 'O=CC(=O)C(=O)C(=O)O(1656)',
    structure = adjacencyList("""1  O u0 p2 c0 {8,S} {11,S}
2  O u0 p2 c0 {7,D}
3  O u0 p2 c0 {6,D}
4  O u0 p2 c0 {8,D}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {3,D} {7,S} {8,S}
7  C u0 p0 c0 {2,D} {6,S} {9,S}
8  C u0 p0 c0 {1,S} {4,D} {6,S}
9  C u0 p0 c0 {5,D} {7,S} {10,S}
10 H u0 p0 c0 {9,S}
11 H u0 p0 c0 {1,S}
"""),
    E0 = (-670.803,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.055,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.271027,0.0937719,-0.000173196,1.59127e-07,-5.52656e-11,-80555.9,29.1296], Tmin=(100,'K'), Tmax=(861.904,'K')), NASAPolynomial(coeffs=[9.2009,0.0274246,-1.43873e-05,2.76654e-09,-1.8839e-13,-81170.2,-7.25836], Tmin=(861.904,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-670.803,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cds-O2d(Cds-O2d)(Cds-O2d)) + group(Cds-O2d(Cds-O2d)(Cds-O2d)) + group(Cds-O2d(Cds-O2d)O2s) + group(Cds-O2d(Cds-O2d)H)"""),
)

species(
    label = '[O]C(=O)C(=O)C(=O)COO(1657)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {2,S} {7,S}
2  O u0 p2 c0 {1,S} {13,S}
3  O u0 p2 c0 {8,D}
4  O u0 p2 c0 {9,D}
5  O u1 p2 c0 {10,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {8,S} {11,S} {12,S}
8  C u0 p0 c0 {3,D} {7,S} {9,S}
9  C u0 p0 c0 {4,D} {8,S} {10,S}
10 C u0 p0 c0 {5,S} {6,D} {9,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {2,S}
"""),
    E0 = (-451.622,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2850,1437.5,1250,1305,750,350,365,385,505,600,445,480,1700,1720,180,180,180,180,1588.21,1600,2924.45,3200],'cm^-1')),
        HinderedRotor(inertia=(0.148848,'amu*angstrom^2'), symmetry=1, barrier=(3.42231,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148848,'amu*angstrom^2'), symmetry=1, barrier=(3.42231,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148848,'amu*angstrom^2'), symmetry=1, barrier=(3.42231,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148848,'amu*angstrom^2'), symmetry=1, barrier=(3.42231,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148848,'amu*angstrom^2'), symmetry=1, barrier=(3.42231,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (147.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.696492,0.114854,-0.000202861,1.80716e-07,-6.21592e-11,-54159.6,35.1222], Tmin=(100,'K'), Tmax=(828.12,'K')), NASAPolynomial(coeffs=[12.6524,0.0313275,-1.7064e-05,3.36695e-09,-2.34298e-13,-55717.4,-22.8175], Tmin=(828.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-451.622,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-O2d)H) + group(O2s-OsH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)(Cds-O2d)) + group(Cds-O2d(Cds-O2d)O2s) + radical(C=OC=OOJ)"""),
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
    E0 = (-146.178,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (205.238,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-25.821,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (28.7945,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (320.076,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (2.86512,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (92.6136,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (68.4154,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-11.3814,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (99.0367,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-78.0446,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-120.496,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-29.7005,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-40.387,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-45.9033,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (159.557,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-11.2046,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (11.183,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]OCC(=O)C(=O)C(=O)O(1618)'],
    products = ['CO2(4)', '[O]OCC(O)=C=O(1445)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(2.46545e+14,'s^-1'), n=0.20628, Ea=(24.4368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(3000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R!H->C',), comment="""Estimated from node Root_N-4R!H->C in family Retroene."""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=C(OO[O])C(=O)C(=O)O(1645)'],
    products = ['[O]OCC(=O)C(=O)C(=O)O(1618)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.30791e+26,'s^-1'), n=-3.67984, Ea=(150.354,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.4163866947181062, var=128.1927359939506, Tref=1000.0, N=20, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root in family 1,3_sigmatropic_rearrangement."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[O]OCC(=C=O)OC(=O)O(1646)'],
    products = ['[O]OCC(=O)C(=O)C(=O)O(1618)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(7.50839e+11,'s^-1'), n=0.136197, Ea=(123.148,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.01703312323488117, var=1.3519057502930525, Tref=1000.0, N=4, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R in family 1,3_sigmatropic_rearrangement."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O]OCOC(=C=O)C(=O)O(1647)'],
    products = ['[O]OCC(=O)C(=O)C(=O)O(1618)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(7.50839e+11,'s^-1'), n=0.136197, Ea=(83.6266,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.01703312323488117, var=1.3519057502930525, Tref=1000.0, N=4, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R in family 1,3_sigmatropic_rearrangement."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]OCC(=O)C(=C=O)OO(1648)'],
    products = ['[O]OCC(=O)C(=O)C(=O)O(1618)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.30791e+26,'s^-1'), n=-3.67984, Ea=(95.044,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.4163866947181062, var=128.1927359939506, Tref=1000.0, N=20, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root in family 1,3_sigmatropic_rearrangement."""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]OCC(=O)OC(O)=C=O(1609)'],
    products = ['[O]OCC(=O)C(=O)C(=O)O(1618)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(7.50839e+11,'s^-1'), n=0.136197, Ea=(100.952,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.01703312323488117, var=1.3519057502930525, Tref=1000.0, N=4, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R in family 1,3_sigmatropic_rearrangement."""),
)

reaction(
    label = 'reaction7',
    reactants = ['O(10)', '[O]CC(=O)C(=O)C(=O)O(1649)'],
    products = ['[O]OCC(=O)C(=O)C(=O)O(1618)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(15399.6,'m^3/(mol*s)'), n=0.932885, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]OCC(=O)C(=O)C(=O)O(1618)'],
    products = ['O=C(O)C(=O)[C]1COOO1(1650)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(6.35e+09,'s^-1'), n=0.19, Ea=(239.03,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Backbone2_N-Sp-3R!H=1R!H_N-1R!H-inRing_Sp-4R!H-3R!H_Ext-2R!H-R_Ext-6R!H-R',), comment="""Estimated from node Backbone2_N-Sp-3R!H=1R!H_N-1R!H-inRing_Sp-4R!H-3R!H_Ext-2R!H-R_Ext-6R!H-R in family Intra_R_Add_Endocyclic."""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]OCC(=O)C(=O)C(=O)O(1618)'],
    products = ['O=C(O)[C]1OOOCC1=O(1651)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4e+08,'s^-1'), n=0.19, Ea=(159.233,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Backbone3_N-Sp-4R!H=1R!H_Sp-2R!H-1R!H_Ext-3R!H-R_Ext-7R!H-R',), comment="""Estimated from node Backbone3_N-Sp-4R!H=1R!H_Sp-2R!H-1R!H_Ext-3R!H-R_Ext-7R!H-R in family Intra_R_Add_Endocyclic.
Ea raised from 151.8 to 159.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]OCC(=O)C(=O)C(=O)O(1618)'],
    products = ['O=C1COOO[C](O)C1=O(1652)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(5.649e+12,'s^-1'), n=0.287, Ea=(269.651,'kJ/mol'), T0=(1,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Backbone4_N-1R!H-inRing_N-7R!H->C',), comment="""Estimated from node Backbone4_N-1R!H-inRing_N-7R!H->C in family Intra_R_Add_Endocyclic.
Ea raised from 264.5 to 269.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]OCC(=O)C(=O)C(=O)O(1618)'],
    products = ['[O]C1(C(=O)C(=O)O)COO1(1653)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.29992e-11,'s^-1'), n=6.41391, Ea=(92.57,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.060263897582951434, var=9.993740724871635, Tref=1000.0, N=74, data_mean=0.0, correlation='Backbone2',), comment="""Estimated from node Backbone2 in family Intra_R_Add_Exocyclic.
Ea raised from 90.2 to 92.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]OCC(=O)C(=O)C(=O)O(1618)'],
    products = ['[O]C1(C(=O)O)OOCC1=O(1654)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(0.234039,'s^-1'), n=3.38729, Ea=(50.1186,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.23174178600410858, var=45.068923389250436, Tref=1000.0, N=140, data_mean=0.0, correlation='Backbone3',), comment="""Estimated from node Backbone3 in family Intra_R_Add_Exocyclic."""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]OCC(=O)C(=O)C(=O)O(1618)'],
    products = ['[O]C1(O)OOCC(=O)C1=O(1655)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.87342e+07,'s^-1'), n=1.00417, Ea=(140.914,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0654062900291529, var=6.730759255238982, Tref=1000.0, N=75, data_mean=0.0, correlation='Backbone4',), comment="""Estimated from node Backbone4 in family Intra_R_Add_Exocyclic."""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O]OC=C(O)C(=O)C(=O)O(1619)'],
    products = ['[O]OCC(=O)C(=O)C(=O)O(1618)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(3.39196e-35,'s^-1'), n=13.7658, Ea=(149.499,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.005573865076554089, var=7.000467719691818, Tref=1000.0, N=12, data_mean=0.0, correlation='Root_3R!H->C_N-3C-inRing_N-1R!H->N_Ext-3C-R_Ext-5R!H-R',), comment="""Estimated from node Root_3R!H->C_N-3C-inRing_N-1R!H->N_Ext-3C-R_Ext-5R!H-R in family Ketoenol."""),
)

reaction(
    label = 'reaction15',
    reactants = ['O2(2)', '[CH2]C(=O)C(=O)C(=O)O(1468)'],
    products = ['[O]OCC(=O)C(=O)C(=O)O(1618)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.96276e+06,'m^3/(mol*s)'), n=-0.119415, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R_Sp-5R!H-4R!H',), comment="""Estimated from node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R_Sp-5R!H-4R!H in family R_Recombination.
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O]OCC(=O)C(=O)C(=O)O(1618)'],
    products = ['[O]OC=C=O(111)', 'O=C=C(O)O(1297)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(4.02168e+12,'s^-1'), n=0.0883205, Ea=(330.171,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.021492990850914453, var=5.303057773824753, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_4R!H->C_N-1R!H->O_2R!H->C_Ext-1C-R_N-7R!H->N_Ext-4C-R_7BrCClFILiOPSSi->O',), comment="""Estimated from node Root_4R!H->C_N-1R!H->O_2R!H->C_Ext-1C-R_N-7R!H->N_Ext-4C-R_7BrCClFILiOPSSi->O in family Retroene.
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O]OCC(=O)C(=O)C(=O)O(1618)'],
    products = ['OH(9)', 'O=CC(=O)C(=O)C(=O)O(1656)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(3.32e+07,'s^-1'), n=1.69, Ea=(159.41,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_SS_O;O_rad_out;Cs_H_out_H/OneDe] for rate rule [R3H_SS_O;O_rad_out;Cs_H_out_H/CO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[O]C(=O)C(=O)C(=O)COO(1657)'],
    products = ['[O]OCC(=O)C(=O)C(=O)O(1618)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.25419e+06,'s^-1'), n=1.03067, Ea=(72.2138,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7H;O_rad_out;XH_out] for rate rule [R7H;O_rad_out;O_H_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #1110',
    isomers = [
        '[O]OCC(=O)C(=O)C(=O)O(1618)',
    ],
    reactants = [
        ('CO2(4)', '[O]OCC(O)=C=O(1445)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1110',
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

