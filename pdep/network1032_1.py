species(
    label = '[CH2]C(=O)C(=O)CO(1485)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {4,S} {12,S}
2  O u0 p2 c0 {5,D}
3  O u0 p2 c0 {6,D}
4  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
5  C u0 p0 c0 {2,D} {4,S} {6,S}
6  C u0 p0 c0 {3,D} {5,S} {7,S}
7  C u1 p0 c0 {6,S} {10,S} {11,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {1,S}
"""),
    E0 = (-282.719,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,365,385,505,600,445,480,1700,1720,3000,3100,440,815,1455,1000,343.384,1240.35],'cm^-1')),
        HinderedRotor(inertia=(0.106567,'amu*angstrom^2'), symmetry=1, barrier=(8.79485,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0154236,'amu*angstrom^2'), symmetry=1, barrier=(1.27953,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.109606,'amu*angstrom^2'), symmetry=1, barrier=(8.83007,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.301786,'amu*angstrom^2'), symmetry=1, barrier=(24.3226,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25463,0.0638513,-7.30318e-05,4.52939e-08,-1.1503e-11,-33907.3,25.4825], Tmin=(100,'K'), Tmax=(947.19,'K')), NASAPolynomial(coeffs=[10.3794,0.025317,-1.20071e-05,2.34213e-09,-1.66271e-13,-35635.9,-18.044], Tmin=(947.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-282.719,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)Cs) + radical(CJCC=O)"""),
)

species(
    label = 'CH2O(20)',
    structure = adjacencyList("""1 O u0 p2 c0 {2,D}
2 C u0 p0 c0 {1,D} {3,S} {4,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (-119.277,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (30.026,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2167.5,'J/mol'), sigma=(4.1028,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=338.56 K, Pc=71.21 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.13878,-0.00469509,2.25728e-05,-2.09848e-08,6.36114e-12,-14349.3,3.23829], Tmin=(100,'K'), Tmax=(1041.97,'K')), NASAPolynomial(coeffs=[2.36098,0.00766801,-3.19768e-06,6.0472e-10,-4.27514e-14,-14279.5,10.4456], Tmin=(1041.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-119.277,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""CH2O""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C=C(O)[C]=O(1233)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {3,S} {8,S}
2 O u0 p2 c0 {5,D}
3 C u0 p0 c0 {1,S} {4,D} {5,S}
4 C u0 p0 c0 {3,D} {6,S} {7,S}
5 C u1 p0 c0 {2,D} {3,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {1,S}
"""),
    E0 = (-127.274,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(1.03719,'amu*angstrom^2'), symmetry=1, barrier=(23.847,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0374,'amu*angstrom^2'), symmetry=1, barrier=(23.8518,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0546,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3730.82,'J/mol'), sigma=(5.21133,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=582.75 K, Pc=59.81 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.43271,0.0523177,-6.55172e-05,3.78376e-08,-8.27375e-12,-15211.5,15.8439], Tmin=(100,'K'), Tmax=(1137.99,'K')), NASAPolynomial(coeffs=[15.2759,0.00365943,-1.38033e-06,2.64519e-10,-1.95358e-14,-18362.2,-52.7311], Tmin=(1137.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-127.274,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(C=CCJ=O)"""),
)

species(
    label = 'O=[C]CC(=O)CO(1515)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {4,S} {12,S}
2  O u0 p2 c0 {6,D}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
5  C u0 p0 c0 {6,S} {7,S} {10,S} {11,S}
6  C u0 p0 c0 {2,D} {4,S} {5,S}
7  C u1 p0 c0 {3,D} {5,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {1,S}
"""),
    E0 = (-312.953,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.823447,0.0753463,-0.000115262,9.55108e-08,-3.12615e-11,-37530.4,25.9028], Tmin=(100,'K'), Tmax=(835.314,'K')), NASAPolynomial(coeffs=[9.77749,0.0250542,-1.16368e-05,2.18059e-09,-1.48514e-13,-38767.6,-14.1357], Tmin=(835.314,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-312.953,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)(Cds-O2d)HH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + radical(C=OCCJ=O)"""),
)

species(
    label = 'C=C([O])C(=C)OO(1516)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {12,S}
3  O u1 p2 c0 {5,S}
4  C u0 p0 c0 {1,S} {5,S} {6,D}
5  C u0 p0 c0 {3,S} {4,S} {7,D}
6  C u0 p0 c0 {4,D} {10,S} {11,S}
7  C u0 p0 c0 {5,D} {8,S} {9,S}
8  H u0 p0 c0 {7,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {2,S}
"""),
    E0 = (-71.565,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,325,375,415,465,420,450,1700,1750,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.684355,'amu*angstrom^2'), symmetry=1, barrier=(15.7347,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.684219,'amu*angstrom^2'), symmetry=1, barrier=(15.7315,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.684326,'amu*angstrom^2'), symmetry=1, barrier=(15.734,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.736475,0.0729126,-9.24337e-05,6.05715e-08,-1.56305e-11,-8490.7,23.7442], Tmin=(100,'K'), Tmax=(951.219,'K')), NASAPolynomial(coeffs=[13.4501,0.0194504,-8.12791e-06,1.48544e-09,-1.01462e-13,-10909.4,-36.9558], Tmin=(951.219,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-71.565,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C([C]=O)OCO(1517)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {4,S} {12,S}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
5  C u0 p0 c0 {1,S} {6,D} {7,S}
6  C u0 p0 c0 {5,D} {10,S} {11,S}
7  C u1 p0 c0 {3,D} {5,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {2,S}
"""),
    E0 = (-300.547,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.44304,0.0879325,-9.95378e-05,4.92992e-08,-8.79531e-12,-35924,29.7383], Tmin=(100,'K'), Tmax=(1641.96,'K')), NASAPolynomial(coeffs=[28.0607,-0.00521652,4.99304e-06,-1.03336e-09,6.96487e-14,-42745,-118.498], Tmin=(1641.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-300.547,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-OsOsHH) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(C=CCJ=O)"""),
)

species(
    label = 'CH2(19)',
    structure = adjacencyList("""multiplicity 3
1 C u2 p0 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
"""),
    E0 = (381.37,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1066.91,2790.99,3622.37],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1235.53,'J/mol'), sigma=(3.758e-10,'m'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with fixed Lennard Jones Parameters. This is the fallback method! Try improving transport databases!"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.01192,-0.000154979,3.26298e-06,-2.40422e-09,5.69497e-13,45867.7,0.5332], Tmin=(100,'K'), Tmax=(1104.57,'K')), NASAPolynomial(coeffs=[3.14983,0.00296674,-9.76055e-07,1.54115e-10,-9.50337e-15,46058.1,4.77807], Tmin=(1104.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'O=[C]C(=O)CO(1518)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {4,S} {9,S}
2 O u0 p2 c0 {5,D}
3 O u0 p2 c0 {6,D}
4 C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
5 C u0 p0 c0 {2,D} {4,S} {6,S}
6 C u1 p0 c0 {3,D} {5,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {1,S}
"""),
    E0 = (-271.191,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,1855,455,950,241.657],'cm^-1')),
        HinderedRotor(inertia=(0.361842,'amu*angstrom^2'), symmetry=1, barrier=(14.8005,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00288628,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.357007,'amu*angstrom^2'), symmetry=1, barrier=(14.8001,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (87.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.90459,0.0480492,-5.78665e-05,3.66854e-08,-9.34293e-12,-32542.9,20.8079], Tmin=(100,'K'), Tmax=(953.387,'K')), NASAPolynomial(coeffs=[9.53906,0.0160181,-7.47057e-06,1.44526e-09,-1.02127e-13,-33998.7,-15.6596], Tmin=(953.387,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-271.191,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)H) + radical(CCCJ=O)"""),
)

species(
    label = '[O]C(CO)=C1CO1(1519)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {4,S} {6,S}
2  O u0 p2 c0 {5,S} {12,S}
3  O u1 p2 c0 {7,S}
4  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
5  C u0 p0 c0 {2,S} {7,S} {10,S} {11,S}
6  C u0 p0 c0 {1,S} {4,S} {7,D}
7  C u0 p0 c0 {3,S} {5,S} {6,D}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {2,S}
"""),
    E0 = (-241.63,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11559,0.0513572,-2.02179e-05,-2.41504e-08,1.71956e-11,-28946.3,22.9451], Tmin=(100,'K'), Tmax=(912.997,'K')), NASAPolynomial(coeffs=[18.6487,0.00682676,-9.96608e-08,-1.09102e-10,6.92028e-15,-33493.4,-67.4153], Tmin=(912.997,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-241.63,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + ring(methyleneoxirane) + radical(C=C(C)OJ)"""),
)

species(
    label = '[O]C1=C(CO)OC1(1520)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {4,S} {6,S}
2  O u0 p2 c0 {5,S} {12,S}
3  O u1 p2 c0 {7,S}
4  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
5  C u0 p0 c0 {2,S} {6,S} {10,S} {11,S}
6  C u0 p0 c0 {1,S} {5,S} {7,D}
7  C u0 p0 c0 {3,S} {4,S} {6,D}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {2,S}
"""),
    E0 = (-266.862,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28827,0.0451612,-1.00301e-06,-4.29573e-08,2.32112e-11,-31985,23.0959], Tmin=(100,'K'), Tmax=(928.211,'K')), NASAPolynomial(coeffs=[18.5898,0.00781115,-7.74393e-07,6.52702e-11,-8.05099e-15,-36799.8,-67.7195], Tmin=(928.211,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-266.862,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + ring(Cyclobutene) + radical(C=C(C)OJ)"""),
)

species(
    label = '[O]C1(CO)CC1=O(1521)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {6,S} {12,S}
2  O u1 p2 c0 {4,S}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
5  C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
6  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
7  C u0 p0 c0 {3,D} {4,S} {5,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {1,S}
"""),
    E0 = (-119.836,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06615,0.0642544,-7.57193e-05,4.78798e-08,-1.2011e-11,-14307.1,24.2319], Tmin=(100,'K'), Tmax=(976.682,'K')), NASAPolynomial(coeffs=[11.9174,0.0198126,-7.46396e-06,1.28914e-09,-8.51053e-14,-16426.7,-27.8629], Tmin=(976.682,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-119.836,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsCs) + ring(cyclopropanone) + radical(CC(C)(C=O)OJ)"""),
)

species(
    label = 'C=C([O])C(O)=CO(1522)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {4,S} {12,S}
2  O u0 p2 c0 {6,S} {11,S}
3  O u1 p2 c0 {5,S}
4  C u0 p0 c0 {1,S} {5,S} {6,D}
5  C u0 p0 c0 {3,S} {4,S} {7,D}
6  C u0 p0 c0 {2,S} {4,D} {8,S}
7  C u0 p0 c0 {5,D} {9,S} {10,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {1,S}
"""),
    E0 = (-404.852,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,325,375,415,465,420,450,1700,1750,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.04113,'amu*angstrom^2'), symmetry=1, barrier=(23.9376,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03813,'amu*angstrom^2'), symmetry=1, barrier=(23.8685,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04376,'amu*angstrom^2'), symmetry=1, barrier=(23.998,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.50757,0.0910907,-0.000110308,5.90203e-08,-1.12556e-11,-48468.3,27.6091], Tmin=(100,'K'), Tmax=(1565.85,'K')), NASAPolynomial(coeffs=[26.3383,-0.00592887,7.43008e-06,-1.66482e-09,1.18932e-13,-54015.2,-109.084], Tmin=(1565.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-404.852,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
)

species(
    label = 'CH2CO(27)',
    structure = adjacencyList("""1 O u0 p2 c0 {3,D}
2 C u0 p0 c0 {3,D} {4,S} {5,S}
3 C u0 p0 c0 {1,D} {2,D}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
"""),
    E0 = (-60.9858,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2120,512.5,787.5],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0366,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2433.35,'J/mol'), sigma=(4.2737,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=380.08 K, Pc=70.74 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.483,0.00833006,5.93357e-06,-1.33931e-08,5.73096e-12,-7313.53,6.51903], Tmin=(100,'K'), Tmax=(954.874,'K')), NASAPolynomial(coeffs=[5.8825,0.00580165,-1.91261e-06,3.35906e-10,-2.37389e-14,-8114.75,-6.74227], Tmin=(954.874,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-60.9858,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""ketene""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'O=[C]CO(129)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {3,S} {7,S}
2 O u0 p2 c0 {4,D}
3 C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
4 C u1 p0 c0 {2,D} {3,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {1,S}
"""),
    E0 = (-177.837,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(1.20761,'amu*angstrom^2'), symmetry=1, barrier=(27.7654,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.21002,'amu*angstrom^2'), symmetry=1, barrier=(27.8206,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (59.044,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.8007,0.0149395,2.57533e-05,-4.93228e-08,2.10358e-11,-21335.3,13.9997], Tmin=(100,'K'), Tmax=(958.142,'K')), NASAPolynomial(coeffs=[12.515,0.00304265,-4.8669e-07,1.51278e-10,-1.80956e-14,-24512.3,-39.3154], Tmin=(958.142,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-177.837,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CsCJ=O)"""),
)

species(
    label = 'CC(=O)C(=O)[CH]O(1523)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {7,S} {12,S}
2  O u0 p2 c0 {5,D}
3  O u0 p2 c0 {6,D}
4  C u0 p0 c0 {5,S} {8,S} {9,S} {10,S}
5  C u0 p0 c0 {2,D} {4,S} {6,S}
6  C u0 p0 c0 {3,D} {5,S} {7,S}
7  C u1 p0 c0 {1,S} {6,S} {11,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {1,S}
"""),
    E0 = (-355.823,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4145,0.0612323,-6.80599e-05,4.31738e-08,-1.15277e-11,-42706.3,23.3419], Tmin=(100,'K'), Tmax=(892.235,'K')), NASAPolynomial(coeffs=[8.5977,0.0290292,-1.39211e-05,2.72207e-09,-1.933e-13,-43988.1,-10.4938], Tmin=(892.235,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-355.823,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)Cs) + radical(OCJC=O)"""),
)

species(
    label = 'CC(=O)C(=O)C[O](1524)',
    structure = adjacencyList("""multiplicity 2
1  O u1 p2 c0 {4,S}
2  O u0 p2 c0 {6,D}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
5  C u0 p0 c0 {7,S} {10,S} {11,S} {12,S}
6  C u0 p0 c0 {2,D} {4,S} {7,S}
7  C u0 p0 c0 {3,D} {5,S} {6,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (-250.575,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,365,385,505,600,445,480,1700,1720,403.596,407.659,2325.34],'cm^-1')),
        HinderedRotor(inertia=(0.0524641,'amu*angstrom^2'), symmetry=1, barrier=(6.09332,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0509284,'amu*angstrom^2'), symmetry=1, barrier=(6.04469,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152092,'amu*angstrom^2'), symmetry=1, barrier=(17.5398,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.65713,0.0555992,-5.64733e-05,3.33184e-08,-8.43831e-12,-30056.3,23.4933], Tmin=(100,'K'), Tmax=(928.173,'K')), NASAPolynomial(coeffs=[7.77774,0.0292224,-1.38465e-05,2.70154e-09,-1.91828e-13,-31192.5,-5.57887], Tmin=(928.173,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-250.575,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)Cs) + radical(C=OCOJ)"""),
)

species(
    label = '[CH2]C(=O)O(130)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {3,S} {7,S}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,S} {2,D} {4,S}
4 C u1 p0 c0 {3,S} {5,S} {6,S}
5 H u0 p0 c0 {4,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {1,S}
"""),
    E0 = (-247.753,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3000,3100,440,815,1455,1000,412.99,413,413.009,413.011],'cm^-1')),
        HinderedRotor(inertia=(0.000988382,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000988378,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (59.044,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.92689,0.0179294,4.49451e-06,-1.67266e-08,6.86775e-12,-29754.3,12.035], Tmin=(100,'K'), Tmax=(1057.99,'K')), NASAPolynomial(coeffs=[7.9836,0.0118142,-5.27089e-06,1.04332e-09,-7.61682e-14,-31552,-16.0851], Tmin=(1057.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-247.753,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""CH2COOH""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'O=[C]C(=O)CCO(1217)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {5,S} {12,S}
2  O u0 p2 c0 {6,D}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {5,S} {6,S} {8,S} {9,S}
5  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
6  C u0 p0 c0 {2,D} {4,S} {7,S}
7  C u1 p0 c0 {3,D} {6,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {1,S}
"""),
    E0 = (-311.182,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,1855,455,950,244.928,245.168],'cm^-1')),
        HinderedRotor(inertia=(0.00280568,'amu*angstrom^2'), symmetry=1, barrier=(0.120457,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16282,'amu*angstrom^2'), symmetry=1, barrier=(6.92241,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157761,'amu*angstrom^2'), symmetry=1, barrier=(6.92303,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.241927,'amu*angstrom^2'), symmetry=1, barrier=(10.2819,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00477,0.0729253,-0.000114268,1.00507e-07,-3.48562e-11,-37325.5,26.097], Tmin=(100,'K'), Tmax=(834.43,'K')), NASAPolynomial(coeffs=[7.39545,0.0298063,-1.43136e-05,2.71927e-09,-1.86544e-13,-37957.4,-0.973347], Tmin=(834.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-311.182,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)H) + radical(CCCJ=O)"""),
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
    label = 'C=C=C([O])CO(1525)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {3,S} {11,S}
2  O u1 p2 c0 {4,S}
3  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
4  C u0 p0 c0 {2,S} {3,S} {6,D}
5  C u0 p0 c0 {6,D} {9,S} {10,S}
6  C u0 p0 c0 {4,D} {5,D}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {1,S}
"""),
    E0 = (-77.7518,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2950,3100,1380,975,1025,1650,540,610,2055,180,4000],'cm^-1')),
        HinderedRotor(inertia=(0.292221,'amu*angstrom^2'), symmetry=1, barrier=(6.71874,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.291977,'amu*angstrom^2'), symmetry=1, barrier=(6.71313,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0812,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.67008,0.0485362,-4.93828e-05,2.64069e-08,-5.59981e-12,-9265.18,11.4876], Tmin=(100,'K'), Tmax=(1148.89,'K')), NASAPolynomial(coeffs=[11.312,0.0149666,-5.55391e-06,9.74253e-10,-6.55975e-14,-11480.7,-36.3673], Tmin=(1148.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-77.7518,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cd-(Cd)CddHH) + group(Cdd-CdsCds) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C1=C(CO)OO1(1526)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {6,S}
3  O u0 p2 c0 {4,S} {12,S}
4  C u0 p0 c0 {3,S} {5,S} {8,S} {9,S}
5  C u0 p0 c0 {1,S} {4,S} {6,D}
6  C u0 p0 c0 {2,S} {5,D} {7,S}
7  C u1 p0 c0 {6,S} {10,S} {11,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {3,S}
"""),
    E0 = (66.1108,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.42258,0.0510852,-4.08437e-05,1.65452e-08,-2.69329e-12,8048.63,25.9987], Tmin=(100,'K'), Tmax=(1455.51,'K')), NASAPolynomial(coeffs=[13.0712,0.0190726,-7.85245e-06,1.43416e-09,-9.78027e-14,4657.7,-34.5716], Tmin=(1455.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(66.1108,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + ring(Cyclobutene) + radical(C=C(O)CJ)"""),
)

species(
    label = 'C=C1OC1([O])CO(1527)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {4,S} {6,S}
2  O u0 p2 c0 {5,S} {12,S}
3  O u1 p2 c0 {4,S}
4  C u0 p0 c0 {1,S} {3,S} {5,S} {6,S}
5  C u0 p0 c0 {2,S} {4,S} {8,S} {9,S}
6  C u0 p0 c0 {1,S} {4,S} {7,D}
7  C u0 p0 c0 {6,D} {10,S} {11,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {2,S}
"""),
    E0 = (-157.778,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.25106,0.0796515,-8.7456e-05,4.31484e-08,-7.58979e-12,-18756.2,27.2256], Tmin=(100,'K'), Tmax=(1719.85,'K')), NASAPolynomial(coeffs=[22.751,-0.00123616,4.95184e-06,-1.14536e-09,8.05721e-14,-23305.3,-90.8088], Tmin=(1719.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-157.778,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(methyleneoxirane) + radical(C=CC(C)(O)OJ)"""),
)

species(
    label = 'CC(=O)[C]=O(257)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {4,D}
2 O u0 p2 c0 {5,D}
3 C u0 p0 c0 {4,S} {6,S} {7,S} {8,S}
4 C u0 p0 c0 {1,D} {3,S} {5,S}
5 C u1 p0 c0 {2,D} {4,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {3,S}
"""),
    E0 = (-125.938,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (71.0546,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.68578,0.0288798,-2.18988e-05,8.65015e-09,-1.4294e-12,-15099.6,16.029], Tmin=(100,'K'), Tmax=(1378.83,'K')), NASAPolynomial(coeffs=[7.54034,0.0147966,-6.578e-06,1.24252e-09,-8.6297e-14,-16438.3,-8.95089], Tmin=(1378.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-125.938,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)HHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)H) + radical(CCCJ=O)"""),
)

species(
    label = '[CH]=C(O)C(=O)CO(1528)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {4,S} {10,S}
2  O u0 p2 c0 {6,S} {11,S}
3  O u0 p2 c0 {5,D}
4  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
5  C u0 p0 c0 {3,D} {4,S} {6,S}
6  C u0 p0 c0 {2,S} {5,S} {7,D}
7  C u1 p0 c0 {6,D} {12,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-231.208,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,350,440,435,1725,3120,650,792.5,1650,287.118],'cm^-1')),
        HinderedRotor(inertia=(0.240111,'amu*angstrom^2'), symmetry=1, barrier=(14.0462,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.240111,'amu*angstrom^2'), symmetry=1, barrier=(14.0462,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.240111,'amu*angstrom^2'), symmetry=1, barrier=(14.0462,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.240111,'amu*angstrom^2'), symmetry=1, barrier=(14.0462,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.578707,0.072446,-8.5142e-05,4.90141e-08,-1.09406e-11,-27682.3,24.728], Tmin=(100,'K'), Tmax=(1102.54,'K')), NASAPolynomial(coeffs=[16.6312,0.0142081,-5.91026e-06,1.10581e-09,-7.75633e-14,-31222,-54.2831], Tmin=(1102.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-231.208,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'C=C(O)C([O])=CO(1529)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {4,S} {12,S}
2  O u0 p2 c0 {6,S} {11,S}
3  O u1 p2 c0 {5,S}
4  C u0 p0 c0 {1,S} {5,S} {7,D}
5  C u0 p0 c0 {3,S} {4,S} {6,D}
6  C u0 p0 c0 {2,S} {5,D} {8,S}
7  C u0 p0 c0 {4,D} {9,S} {10,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {1,S}
"""),
    E0 = (-404.852,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.50757,0.0910907,-0.000110308,5.90203e-08,-1.12556e-11,-48468.3,27.6091], Tmin=(100,'K'), Tmax=(1565.85,'K')), NASAPolynomial(coeffs=[26.3383,-0.00592887,7.43008e-06,-1.66482e-09,1.18932e-13,-54015.2,-109.084], Tmin=(1565.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-404.852,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C(O)C(=O)C[O](1479)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {5,S} {12,S}
2  O u1 p2 c0 {4,S}
3  O u0 p2 c0 {6,D}
4  C u0 p0 c0 {2,S} {6,S} {8,S} {9,S}
5  C u0 p0 c0 {1,S} {6,S} {7,D}
6  C u0 p0 c0 {3,D} {4,S} {5,S}
7  C u0 p0 c0 {5,D} {10,S} {11,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {1,S}
"""),
    E0 = (-234.571,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,375,552.5,462.5,1710,2950,3100,1380,975,1025,1650,283.42,285.041,289.512],'cm^-1')),
        HinderedRotor(inertia=(0.00209889,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.211563,'amu*angstrom^2'), symmetry=1, barrier=(12.1212,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.211408,'amu*angstrom^2'), symmetry=1, barrier=(12.1205,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.890155,0.0660378,-7.24238e-05,3.99394e-08,-8.65169e-12,-28098.4,23.9874], Tmin=(100,'K'), Tmax=(1127.78,'K')), NASAPolynomial(coeffs=[14.7349,0.0169328,-7.11115e-06,1.33068e-09,-9.30343e-14,-31221.1,-44.4705], Tmin=(1127.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-234.571,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsHH) + radical(C=OCOJ)"""),
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
    E0 = (-48.2245,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (36.6034,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (262.181,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (22.7103,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (289.467,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (27.3394,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (29.7727,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (59.4517,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-16.8073,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-17.6503,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-41.0475,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-36.8484,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (199.83,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (17.2467,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (344.541,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (261.369,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (22.9631,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-49.4908,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (140.209,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (5.93906,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-17.5763,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(=O)C(=O)CO(1485)'],
    products = ['CH2O(20)', 'C=C(O)[C]=O(1233)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(2.46545e+14,'s^-1'), n=0.20628, Ea=(55.2064,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(3000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R!H->C',), comment="""Estimated from node Root_N-4R!H->C in family Retroene."""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(=O)C(=O)CO(1485)'],
    products = ['O=[C]CC(=O)CO(1515)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.08731e+10,'s^-1'), n=0.796, Ea=(140.034,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCCJ;CsJ;CO] + [cCCJ;CsJ-HH;C] for rate rule [cCCJ;CsJ-HH;CO]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=C([O])C(=C)OO(1516)'],
    products = ['[CH2]C(=O)C(=O)CO(1485)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.30791e+26,'s^-1'), n=-3.67984, Ea=(154.458,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.4163866947181062, var=128.1927359939506, Tref=1000.0, N=20, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root in family 1,3_sigmatropic_rearrangement."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]C(=O)C(=O)CO(1485)'],
    products = ['C=C([C]=O)OCO(1517)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(7.50839e+11,'s^-1'), n=0.136197, Ea=(126.141,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.01703312323488117, var=1.3519057502930525, Tref=1000.0, N=4, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R in family 1,3_sigmatropic_rearrangement."""),
)

reaction(
    label = 'reaction5',
    reactants = ['CH2(19)', 'O=[C]C(=O)CO(1518)'],
    products = ['[CH2]C(=O)C(=O)CO(1485)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.73526e+06,'m^3/(mol*s)'), n=0.377732, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_rad/OneDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]C(=O)C(=O)CO(1485)'],
    products = ['[O]C(CO)=C1CO1(1519)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(6.95187e+18,'s^-1'), n=-1.62813, Ea=(130.77,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.08275800860939324, var=25.73102631834909, Tref=1000.0, N=5, data_mean=0.0, correlation='Backbone0_N-2R!H-inRing_N-1R!H-inRing_Sp-2R!H-1R!H',), comment="""Estimated from node Backbone0_N-2R!H-inRing_N-1R!H-inRing_Sp-2R!H-1R!H in family Intra_R_Add_Endocyclic."""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]C(=O)C(=O)CO(1485)'],
    products = ['[O]C1=C(CO)OC1(1520)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.57562e+06,'s^-1'), n=1.6342, Ea=(133.204,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.025159707188840464, var=2.1106144811101477, Tref=1000.0, N=3, data_mean=0.0, correlation='Backbone1_N-2R!H-inRing_Ext-2R!H-R',), comment="""Estimated from node Backbone1_N-2R!H-inRing_Ext-2R!H-R in family Intra_R_Add_Endocyclic."""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]C(=O)C(=O)CO(1485)'],
    products = ['[O]C1(CO)CC1=O(1521)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(3.36539e+09,'s^-1'), n=0.84129, Ea=(162.883,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.031430614527560546, var=6.099077214950256, Tref=1000.0, N=84, data_mean=0.0, correlation='Backbone1',), comment="""Estimated from node Backbone1 in family Intra_R_Add_Exocyclic.
Ea raised from 162.5 to 162.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=C([O])C(O)=CO(1522)'],
    products = ['[CH2]C(=O)C(=O)CO(1485)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4.7171e-39,'s^-1'), n=14.8922, Ea=(208.757,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.01101064541230846, var=9.816149634195241, Tref=1000.0, N=17, data_mean=0.0, correlation='Root_3R!H->C_N-3C-inRing_N-1R!H->N_Ext-3C-R',), comment="""Estimated from node Root_3R!H->C_N-3C-inRing_N-1R!H->N_Ext-3C-R in family Ketoenol."""),
)

reaction(
    label = 'reaction10',
    reactants = ['CH2CO(27)', 'O=[C]CO(129)'],
    products = ['[CH2]C(=O)C(=O)CO(1485)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.7315,'m^3/(mol*s)'), n=1.87352, Ea=(41.8843,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Ck_Cds-HH;CJ] for rate rule [Ck_Cds-HH;CO_rad/NonDe]
Euclidian distance = 3.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C(=O)C(=O)CO(1485)'],
    products = ['CC(=O)C(=O)[CH]O(1523)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(304,'s^-1'), n=2.77, Ea=(62.3834,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 324 used for R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeO
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeO]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['CC(=O)C(=O)C[O](1524)'],
    products = ['[CH2]C(=O)C(=O)CO(1485)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(63.58,'s^-1','*|/',3.66), n=2.81162, Ea=(34.4385,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 3 used for R5H_CC(O2d)CC;O_rad_out;Cs_H_out_2H
Exact match found for rate rule [R5H_CC(O2d)CC;O_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['CH2CO(27)', '[CH2]C(=O)O(130)'],
    products = ['[CH2]C(=O)C(=O)CO(1485)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.79e-11,'m^3/(mol*s)'), n=3.97, Ea=(329.281,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [doublebond;R_OH] for rate rule [Cdd_Cd_2H;R_OH]
Euclidian distance = 2.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C(=O)C(=O)CO(1485)'],
    products = ['O=[C]C(=O)CCO(1217)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(7.50839e+11,'s^-1'), n=0.136197, Ea=(120.678,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.01703312323488117, var=1.3519057502930525, Tref=1000.0, N=4, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R in family 1,3_sigmatropic_rearrangement."""),
)

reaction(
    label = 'reaction15',
    reactants = ['O(10)', 'C=C=C([O])CO(1525)'],
    products = ['[CH2]C(=O)C(=O)CO(1485)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/OneDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C(=O)C(=O)CO(1485)'],
    products = ['[CH2]C1=C(CO)OO1(1526)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.57562e+06,'s^-1'), n=1.6342, Ea=(364.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.025159707188840464, var=2.1106144811101477, Tref=1000.0, N=3, data_mean=0.0, correlation='Backbone1_N-2R!H-inRing_Ext-2R!H-R',), comment="""Estimated from node Backbone1_N-2R!H-inRing_Ext-2R!H-R in family Intra_R_Add_Endocyclic."""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C(=O)C(=O)CO(1485)'],
    products = ['C=C1OC1([O])CO(1527)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(3.36539e+09,'s^-1'), n=0.84129, Ea=(126.394,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.031430614527560546, var=6.099077214950256, Tref=1000.0, N=84, data_mean=0.0, correlation='Backbone1',), comment="""Estimated from node Backbone1 in family Intra_R_Add_Exocyclic."""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C(=O)C(=O)CO(1485)'],
    products = ['CH2O(20)', 'CC(=O)[C]=O(257)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.46545e+14,'s^-1'), n=0.20628, Ea=(53.9401,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(3000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R!H->C',), comment="""Estimated from node Root_N-4R!H->C in family Retroene."""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=C(O)C(=O)CO(1528)'],
    products = ['[CH2]C(=O)C(=O)CO(1485)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C(=O)C(=O)CO(1485)'],
    products = ['C=C(O)C([O])=CO(1529)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(3.22e+08,'s^-1'), n=1.09, Ea=(109.37,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4H_SSS;O_rad_out;Cs_H_out_H/NonDeO] for rate rule [R4H_S(Cd)SS;O_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C=C(O)C(=O)C[O](1479)'],
    products = ['[CH2]C(=O)C(=O)CO(1485)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(22219.5,'s^-1'), n=1.84103, Ea=(37.7069,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5H_CCC;O_rad_out;XH_out] + [R5H_CCC(Cd);Y_rad_out;XH_out] for rate rule [R5H_CCC(Cd);O_rad_out;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #1032',
    isomers = [
        '[CH2]C(=O)C(=O)CO(1485)',
    ],
    reactants = [
        ('CH2O(20)', 'C=C(O)[C]=O(1233)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1032',
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

