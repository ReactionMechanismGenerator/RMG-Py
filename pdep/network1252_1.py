species(
    label = 'O=CC=C(OO)OCO(1849)',
    structure = adjacencyList("""1  O u0 p2 c0 {6,S} {7,S}
2  O u0 p2 c0 {4,S} {7,S}
3  O u0 p2 c0 {6,S} {14,S}
4  O u0 p2 c0 {2,S} {15,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {1,S} {3,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {8,D}
8  C u0 p0 c0 {7,D} {9,S} {12,S}
9  C u0 p0 c0 {5,D} {8,S} {13,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
"""),
    E0 = (-496.046,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,2782.5,750,1395,475,1775,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.087,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.02759,0.0994492,-0.00011155,5.85887e-08,-1.17479e-11,-59470,33.7205], Tmin=(100,'K'), Tmax=(1234.07,'K')), NASAPolynomial(coeffs=[26.0613,0.0116447,-4.82364e-06,9.32269e-10,-6.75672e-14,-66155.9,-102.665], Tmin=(1234.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-496.046,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-OsOsHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H)"""),
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
    label = 'O=CCC(=O)OO(1774)',
    structure = adjacencyList("""1  O u0 p2 c0 {2,S} {6,S}
2  O u0 p2 c0 {1,S} {11,S}
3  O u0 p2 c0 {6,D}
4  O u0 p2 c0 {7,D}
5  C u0 p0 c0 {6,S} {7,S} {8,S} {9,S}
6  C u0 p0 c0 {1,S} {3,D} {5,S}
7  C u0 p0 c0 {4,D} {5,S} {10,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {2,S}
"""),
    E0 = (-462.407,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4128.08,'J/mol'), sigma=(5.50619,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=644.80 K, Pc=56.11 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09842,0.0661832,-8.97924e-05,6.45239e-08,-1.83523e-11,-55512.3,23.6145], Tmin=(100,'K'), Tmax=(863.539,'K')), NASAPolynomial(coeffs=[11.2501,0.0191618,-8.11789e-06,1.47254e-09,-9.9282e-14,-57265.7,-23.8723], Tmin=(863.539,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-462.407,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-O2d)) + group(O2s-(Os-CdOd)H) + group(Cs-(Cds-O2d)(Cds-O2d)HH) + group(Cds-OdCsOs) + group(Cds-OdCsH)"""),
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
    label = 'C=C(OO)OCO(1865)',
    structure = adjacencyList("""1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {5,S} {12,S}
4  O u0 p2 c0 {2,S} {13,S}
5  C u0 p0 c0 {1,S} {3,S} {8,S} {9,S}
6  C u0 p0 c0 {1,S} {2,S} {7,D}
7  C u0 p0 c0 {6,D} {10,S} {11,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
"""),
    E0 = (-372.558,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2950,3100,1380,975,1025,1650,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.931631,0.0828224,-8.91814e-05,4.38427e-08,-7.88997e-12,-44608.7,30.3138], Tmin=(100,'K'), Tmax=(1599.53,'K')), NASAPolynomial(coeffs=[24.9544,0.00182873,2.02044e-06,-5.24389e-10,3.7729e-14,-50809.8,-100.228], Tmin=(1599.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-372.558,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-OsOsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH)"""),
)

species(
    label = 'O=CC(CO)C(=O)OO(1866)',
    structure = adjacencyList("""1  O u0 p2 c0 {3,S} {8,S}
2  O u0 p2 c0 {7,S} {14,S}
3  O u0 p2 c0 {1,S} {15,S}
4  O u0 p2 c0 {8,D}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {7,S} {8,S} {9,S} {10,S}
7  C u0 p0 c0 {2,S} {6,S} {11,S} {12,S}
8  C u0 p0 c0 {1,S} {4,D} {6,S}
9  C u0 p0 c0 {5,D} {6,S} {13,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {3,S}
"""),
    E0 = (-664.743,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.087,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0428553,0.0899472,-0.000114082,7.71519e-08,-2.08057e-11,-79810.1,32.9155], Tmin=(100,'K'), Tmax=(906.644,'K')), NASAPolynomial(coeffs=[13.9969,0.0283841,-1.22298e-05,2.25938e-09,-1.54882e-13,-82340.4,-33.0373], Tmin=(906.644,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-664.743,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-O2d)) + group(O2s-CsH) + group(O2s-(Os-CdOd)H) + group(Cs-CsCsCsH) + group(Cs-CsOsHH) + group(Cds-OdCsOs) + group(Cds-OdCsH)"""),
)

species(
    label = 'O=CC(O)C(=O)OCO(1867)',
    structure = adjacencyList("""1  O u0 p2 c0 {7,S} {8,S}
2  O u0 p2 c0 {6,S} {15,S}
3  O u0 p2 c0 {7,S} {14,S}
4  O u0 p2 c0 {8,D}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {2,S} {8,S} {9,S} {10,S}
7  C u0 p0 c0 {1,S} {3,S} {11,S} {12,S}
8  C u0 p0 c0 {1,S} {4,D} {6,S}
9  C u0 p0 c0 {5,D} {6,S} {13,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {2,S}
"""),
    E0 = (-895.002,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.087,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.184482,0.0679769,-2.85296e-05,-2.35666e-08,1.69483e-11,-107492,35.0895], Tmin=(100,'K'), Tmax=(961.644,'K')), NASAPolynomial(coeffs=[22.1259,0.0146565,-4.54771e-06,8.4091e-10,-6.43916e-14,-113467,-79.0296], Tmin=(961.644,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-895.002,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-OsOsHH) + group(Cds-OdCsOs) + group(Cds-OdCsH)"""),
)

species(
    label = 'OC=C=C(OO)OCO(1868)',
    structure = adjacencyList("""1  O u0 p2 c0 {6,S} {7,S}
2  O u0 p2 c0 {5,S} {7,S}
3  O u0 p2 c0 {6,S} {13,S}
4  O u0 p2 c0 {8,S} {14,S}
5  O u0 p2 c0 {2,S} {15,S}
6  C u0 p0 c0 {1,S} {3,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {9,D}
8  C u0 p0 c0 {4,S} {9,D} {12,S}
9  C u0 p0 c0 {7,D} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
"""),
    E0 = (-440.772,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,3580,3650,1210,1345,900,1100,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,540,610,2055,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.087,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.98237,0.118769,-0.000141491,7.43434e-08,-1.40907e-11,-52731.2,38.3811], Tmin=(100,'K'), Tmax=(1537.41,'K')), NASAPolynomial(coeffs=[34.9477,-0.00685828,7.36556e-06,-1.60385e-09,1.12701e-13,-61210,-150.568], Tmin=(1537.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-440.772,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-OsOsHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'CH2OH(32)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {2,S} {5,S}
2 C u1 p0 c0 {1,S} {3,S} {4,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {1,S}
"""),
    E0 = (-26.6187,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,915.949,3793.45],'cm^-1')),
        HinderedRotor(inertia=(0.0105829,'amu*angstrom^2'), symmetry=1, barrier=(6.29861,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (31.034,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3045.92,'J/mol'), sigma=(4.69067,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=475.77 K, Pc=66.97 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.62325,0.00544082,9.71812e-06,-1.58618e-08,6.40702e-12,-3185.33,6.71472], Tmin=(100,'K'), Tmax=(942.605,'K')), NASAPolynomial(coeffs=[5.31889,0.0054271,-1.68868e-06,2.88831e-10,-2.02657e-14,-3824.05,-3.05791], Tmin=(942.605,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-26.6187,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""CH2OH""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[O]C=CC(=O)OO(1813)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {2,S} {6,S}
2  O u0 p2 c0 {1,S} {10,S}
3  O u0 p2 c0 {6,D}
4  O u1 p2 c0 {7,S}
5  C u0 p0 c0 {6,S} {7,D} {8,S}
6  C u0 p0 c0 {1,S} {3,D} {5,S}
7  C u0 p0 c0 {4,S} {5,D} {9,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {2,S}
"""),
    E0 = (-278.785,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,4000,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.280532,'amu*angstrom^2'), symmetry=1, barrier=(6.44999,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.280229,'amu*angstrom^2'), symmetry=1, barrier=(6.44301,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.281616,'amu*angstrom^2'), symmetry=1, barrier=(6.47492,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30371,0.0590968,-8.67305e-05,5.65448e-08,-1.30326e-11,-33432.6,20.6679], Tmin=(100,'K'), Tmax=(825.773,'K')), NASAPolynomial(coeffs=[15.032,0.00187594,3.56095e-07,-1.55673e-10,1.3902e-14,-36016.3,-44.8507], Tmin=(825.773,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-278.785,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-O2d)) + group(O2s-(Cds-Cd)H) + group(O2s-(Os-CdOd)H) + missing(Cd-COCdH0) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsOsH) + radical(C=COJ)"""),
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
    label = '[O]C=CC(=O)OCO(1870)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {5,S} {7,S}
2  O u0 p2 c0 {5,S} {13,S}
3  O u0 p2 c0 {7,D}
4  O u1 p2 c0 {8,S}
5  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
6  C u0 p0 c0 {7,S} {8,D} {11,S}
7  C u0 p0 c0 {1,S} {3,D} {6,S}
8  C u0 p0 c0 {4,S} {6,D} {12,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {2,S}
"""),
    E0 = (-529.223,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,1291.11,2795.47,3839.61,4000,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.438897,'amu*angstrom^2'), symmetry=1, barrier=(10.0911,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.438897,'amu*angstrom^2'), symmetry=1, barrier=(10.0911,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.438897,'amu*angstrom^2'), symmetry=1, barrier=(10.0911,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.438897,'amu*angstrom^2'), symmetry=1, barrier=(10.0911,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.973656,0.0496301,-8.93358e-06,-4.2906e-08,2.52497e-11,-63526.2,26.8498], Tmin=(100,'K'), Tmax=(926.142,'K')), NASAPolynomial(coeffs=[23.0793,-0.00195039,3.51591e-06,-6.93681e-10,4.13235e-14,-69503.3,-88.2643], Tmin=(926.142,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-529.223,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-OsOsHH) + missing(Cd-COCdH0) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsOsH) + radical(C=COJ)"""),
)

species(
    label = 'HO2(12)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {2,S} {3,S}
2 O u1 p2 c0 {1,S}
3 H u0 p0 c0 {1,S}
"""),
    E0 = (2.67648,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1112.81,1388.53,3298.45],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (33.0068,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3486.33,'J/mol'), sigma=(4.38953,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=544.56 K, Pc=93.53 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.02956,-0.00263991,1.52233e-05,-1.71675e-08,6.26753e-12,322.677,4.84426], Tmin=(100,'K'), Tmax=(923.908,'K')), NASAPolynomial(coeffs=[4.15132,0.00191149,-4.11289e-07,6.34993e-11,-4.86415e-15,83.4268,3.09349], Tmin=(923.908,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(2.67648,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HO2""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = 'O=CC=[C]OCO(1871)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {4,S} {7,S}
2  O u0 p2 c0 {4,S} {12,S}
3  O u0 p2 c0 {6,D}
4  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
5  C u0 p0 c0 {6,S} {7,D} {10,S}
6  C u0 p0 c0 {3,D} {5,S} {11,S}
7  C u1 p0 c0 {1,S} {5,D}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {2,S}
"""),
    E0 = (-223.66,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,2782.5,750,1395,475,1775,1000,1685,370,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.862175,'amu*angstrom^2'), symmetry=1, barrier=(19.8231,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.861543,'amu*angstrom^2'), symmetry=1, barrier=(19.8086,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.861073,'amu*angstrom^2'), symmetry=1, barrier=(19.7978,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.861899,'amu*angstrom^2'), symmetry=1, barrier=(19.8168,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.55771,0.0579454,-1.86073e-05,-3.64201e-08,2.30555e-11,-26759.8,27.3245], Tmin=(100,'K'), Tmax=(940.073,'K')), NASAPolynomial(coeffs=[24.549,0.000526376,1.74551e-06,-3.14064e-10,1.32939e-14,-33244,-97.4338], Tmin=(940.073,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-223.66,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-OsOsHH) + group(Cd-Cd(CO)H) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + radical(C=CJO)"""),
)

species(
    label = '[CH2]OC(=CC=O)OO(1872)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {5,S} {8,S}
2  O u0 p2 c0 {3,S} {5,S}
3  O u0 p2 c0 {2,S} {13,S}
4  O u0 p2 c0 {7,D}
5  C u0 p0 c0 {1,S} {2,S} {6,D}
6  C u0 p0 c0 {5,D} {7,S} {9,S}
7  C u0 p0 c0 {4,D} {6,S} {10,S}
8  C u1 p0 c0 {1,S} {11,S} {12,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {3,S}
"""),
    E0 = (-112.082,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,2782.5,750,1395,475,1775,1000,3000,3100,440,815,1455,1000,180,1142.47],'cm^-1')),
        HinderedRotor(inertia=(0.231757,'amu*angstrom^2'), symmetry=1, barrier=(5.32855,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.75077,'amu*angstrom^2'), symmetry=1, barrier=(17.2617,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.228599,'amu*angstrom^2'), symmetry=1, barrier=(5.25593,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.753455,'amu*angstrom^2'), symmetry=1, barrier=(17.3234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.75468,'amu*angstrom^2'), symmetry=1, barrier=(17.3516,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.5311,0.0838708,-0.000123892,1.02959e-07,-3.52196e-11,-13362.5,28.7051], Tmin=(100,'K'), Tmax=(710.105,'K')), NASAPolynomial(coeffs=[9.41999,0.0338007,-1.81281e-05,3.66638e-09,-2.63137e-13,-14624.9,-11.1357], Tmin=(710.105,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-112.082,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-OsHHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=COCJ)"""),
)

species(
    label = 'H(5)',
    structure = adjacencyList("""multiplicity 2
1 H u1 p0 c0
"""),
    E0 = (211.792,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (1.00797,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(2172.56,'J/mol'), sigma=(4.12675,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=339.35 K, Pc=70.14 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,1.49298e-14,-2.05327e-17,9.23927e-21,-1.28064e-24,25472.7,-0.459566], Tmin=(100,'K'), Tmax=(3959.16,'K')), NASAPolynomial(coeffs=[2.5,1.19009e-10,-4.23987e-14,6.68965e-18,-3.94352e-22,25472.7,-0.459565], Tmin=(3959.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(211.792,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""H""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = '[O]COC(=CC=O)OO(1873)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {6,S} {7,S}
2  O u0 p2 c0 {3,S} {7,S}
3  O u0 p2 c0 {2,S} {14,S}
4  O u1 p2 c0 {6,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {8,D}
8  C u0 p0 c0 {7,D} {9,S} {12,S}
9  C u0 p0 c0 {5,D} {8,S} {13,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {3,S}
"""),
    E0 = (-268.648,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,2782.5,750,1395,475,1775,1000,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.376968,0.0863187,-9.88257e-05,5.7979e-08,-1.38782e-11,-32185.9,29.1723], Tmin=(100,'K'), Tmax=(999.011,'K')), NASAPolynomial(coeffs=[13.941,0.0320092,-1.72812e-05,3.56234e-09,-2.60628e-13,-34896,-36.2528], Tmin=(999.011,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-268.648,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-OsOsHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(OCOJ)"""),
)

species(
    label = '[O]OC(=CC=O)OCO(1874)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {6,S} {7,S}
2  O u0 p2 c0 {6,S} {14,S}
3  O u0 p2 c0 {5,S} {7,S}
4  O u0 p2 c0 {9,D}
5  O u1 p2 c0 {3,S}
6  C u0 p0 c0 {1,S} {2,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {3,S} {8,D}
8  C u0 p0 c0 {7,D} {9,S} {12,S}
9  C u0 p0 c0 {4,D} {8,S} {13,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {2,S}
"""),
    E0 = (-344.041,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,2782.5,750,1395,475,1775,1000,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.626438,0.0942976,-0.000111108,6.17184e-08,-1.31058e-11,-41205.6,33.0364], Tmin=(100,'K'), Tmax=(1165.73,'K')), NASAPolynomial(coeffs=[23.4552,0.0116661,-4.78295e-06,9.12723e-10,-6.55813e-14,-46820.1,-86.8365], Tmin=(1165.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-344.041,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-OsOsHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(ROOJ)"""),
)

species(
    label = 'O=CC=C(OO)O[CH]O(1875)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {6,S} {8,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {8,S} {13,S}
4  O u0 p2 c0 {2,S} {14,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {1,S} {2,S} {7,D}
7  C u0 p0 c0 {6,D} {9,S} {10,S}
8  C u1 p0 c0 {1,S} {3,S} {11,S}
9  C u0 p0 c0 {5,D} {7,S} {12,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
"""),
    E0 = (-307.086,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,3615,1277.5,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.669746,0.110962,-0.000175176,1.40891e-07,-4.47717e-11,-36773.4,32.2313], Tmin=(100,'K'), Tmax=(773.454,'K')), NASAPolynomial(coeffs=[15.2726,0.0285106,-1.52656e-05,3.05271e-09,-2.1657e-13,-39239.4,-40.585], Tmin=(773.454,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-307.086,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-OsOsHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(OCJO)"""),
)

species(
    label = 'HCO(16)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {2,D}
2 C u1 p0 c0 {1,D} {3,S}
3 H u0 p0 c0 {2,S}
"""),
    E0 = (32.5964,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1116.11,1837.47,2716.03],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (29.018,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2167.5,'J/mol'), sigma=(4.1028,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=338.56 K, Pc=71.21 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.06159,-0.0018195,9.37895e-06,-8.24742e-09,2.33851e-12,3918.59,3.98389], Tmin=(100,'K'), Tmax=(1105.12,'K')), NASAPolynomial(coeffs=[3.05335,0.00410464,-1.74962e-06,3.28532e-10,-2.28985e-14,4002.52,8.32038], Tmin=(1105.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(32.5964,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HCO""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH]=C(OO)OCO(1876)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {5,S} {10,S}
4  O u0 p2 c0 {2,S} {12,S}
5  C u0 p0 c0 {1,S} {3,S} {8,S} {9,S}
6  C u0 p0 c0 {1,S} {2,S} {7,D}
7  C u1 p0 c0 {6,D} {11,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {4,S}
"""),
    E0 = (-125.462,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.828401,'amu*angstrom^2'), symmetry=1, barrier=(19.0466,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.828393,'amu*angstrom^2'), symmetry=1, barrier=(19.0464,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.828407,'amu*angstrom^2'), symmetry=1, barrier=(19.0467,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.829024,'amu*angstrom^2'), symmetry=1, barrier=(19.0609,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.829229,'amu*angstrom^2'), symmetry=1, barrier=(19.0656,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.744365,0.0834681,-9.66383e-05,5.03918e-08,-9.61326e-12,-14900.8,30.0308], Tmin=(100,'K'), Tmax=(1491.2,'K')), NASAPolynomial(coeffs=[25.3126,-0.000994754,2.97659e-06,-6.93852e-10,4.95565e-14,-21052.3,-100.66], Tmin=(1491.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-125.462,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-OsOsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = '[O]C=C=C(OO)OCO(1877)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {6,S} {7,S}
2  O u0 p2 c0 {4,S} {7,S}
3  O u0 p2 c0 {6,S} {13,S}
4  O u0 p2 c0 {2,S} {14,S}
5  O u1 p2 c0 {8,S}
6  C u0 p0 c0 {1,S} {3,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {9,D}
8  C u0 p0 c0 {5,S} {9,D} {12,S}
9  C u0 p0 c0 {7,D} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
"""),
    E0 = (-299.309,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,540,610,2055,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.31556,0.108748,-0.000127195,6.60522e-08,-1.24599e-11,-35745.3,37.4006], Tmin=(100,'K'), Tmax=(1525.73,'K')), NASAPolynomial(coeffs=[32.5183,-0.00428673,5.27873e-06,-1.15871e-09,8.11841e-14,-43847.6,-137.088], Tmin=(1525.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-299.309,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-OsOsHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ)"""),
)

species(
    label = 'O=C[C]C(OO)OCO(1879)',
    structure = adjacencyList("""1  O u0 p2 c0 {6,S} {7,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {7,S} {14,S}
4  O u0 p2 c0 {2,S} {15,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
7  C u0 p0 c0 {1,S} {3,S} {11,S} {12,S}
8  C u0 p0 c0 {5,D} {9,S} {13,S}
9  C u0 p1 c0 {6,S} {8,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
"""),
    E0 = (-232.633,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.087,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.661808,0.097825,-0.000108843,5.92252e-08,-1.25626e-11,-27807.4,46.1753], Tmin=(100,'K'), Tmax=(1154.62,'K')), NASAPolynomial(coeffs=[21.4399,0.0212568,-9.37068e-06,1.79079e-09,-1.26754e-13,-32911.2,-63.6303], Tmin=(1154.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-232.633,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-OsOsHH) + group(Cds-OdCsH) + group(CsJ2_singlet-CsH)"""),
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
    E0 = (-308.045,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-1.71202,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-246.794,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-218.222,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-143.068,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-120.645,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-306.89,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-36.2249,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (101.049,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (127.903,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (52.5093,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (89.4648,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (91.8931,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (97.2414,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-24.6379,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O=CC=C(OO)OCO(1849)'],
    products = ['CH2O(20)', 'O=CCC(=O)OO(1774)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(2.46545e+14,'s^-1'), n=0.20628, Ea=(3.24294,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(3000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R!H->C',), comment="""Estimated from node Root_N-4R!H->C in family Retroene."""),
)

reaction(
    label = 'reaction2',
    reactants = ['CO(15)', 'C=C(OO)OCO(1865)'],
    products = ['O=CC=C(OO)OCO(1849)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(0.118397,'m^3/(mol*s)'), n=2.3675, Ea=(305.306,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO;R_H] for rate rule [CO;Cd_pri]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction3',
    reactants = ['O=CC=C(OO)OCO(1849)'],
    products = ['O=CC(CO)C(=O)OO(1866)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(7.50839e+11,'s^-1'), n=0.136197, Ea=(64.4931,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.01703312323488117, var=1.3519057502930525, Tref=1000.0, N=4, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R in family 1,3_sigmatropic_rearrangement."""),
)

reaction(
    label = 'reaction4',
    reactants = ['O=CC=C(OO)OCO(1849)'],
    products = ['O=CC(O)C(=O)OCO(1867)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.30791e+26,'s^-1'), n=-3.67984, Ea=(93.065,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.4163866947181062, var=128.1927359939506, Tref=1000.0, N=20, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root in family 1,3_sigmatropic_rearrangement."""),
)

reaction(
    label = 'reaction5',
    reactants = ['OC=C=C(OO)OCO(1868)'],
    products = ['O=CC=C(OO)OCO(1849)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(3.39196e-35,'s^-1'), n=13.7658, Ea=(112.945,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.005573865076554089, var=7.000467719691818, Tref=1000.0, N=12, data_mean=0.0, correlation='Root_3R!H->C_N-3C-inRing_N-1R!H->N_Ext-3C-R_Ext-5R!H-R',), comment="""Estimated from node Root_3R!H->C_N-3C-inRing_N-1R!H->N_Ext-3C-R_Ext-5R!H-R in family Ketoenol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH2OH(32)', '[O]C=CC(=O)OO(1813)'],
    products = ['O=CC=C(OO)OCO(1849)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(106477,'m^3/(mol*s)'), n=0.348287, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.0108230153501, var=2.70964383578, Tref=1000.0, N=19, data_mean=0.0, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R',), comment="""Estimated from node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R in family R_Recombination."""),
)

reaction(
    label = 'reaction7',
    reactants = ['OH(9)', '[O]C=CC(=O)OCO(1870)'],
    products = ['O=CC=C(OO)OCO(1849)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2e+07,'m^3/(mol*s)'), n=1.78837e-07, Ea=(9.20349,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_N-2R->C_N-3R!H->N_2OS->O',), comment="""Estimated from node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_N-2R->C_N-3R!H->N_2OS->O in family R_Recombination."""),
)

reaction(
    label = 'reaction8',
    reactants = ['HO2(12)', 'O=CC=[C]OCO(1871)'],
    products = ['O=CC=C(OO)OCO(1849)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(6.08474e+07,'m^3/(mol*s)'), n=-0.428622, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing',), comment="""Estimated from node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing in family R_Recombination."""),
)

reaction(
    label = 'reaction9',
    reactants = ['OH(9)', '[CH2]OC(=CC=O)OO(1872)'],
    products = ['O=CC=C(OO)OCO(1849)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(7.7e+07,'m^3/(mol*s)'), n=4.95181e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_2R->C_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1CNOS->N_1COS->O_2R->C_Ext-2C-R in family R_Recombination."""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(5)', '[O]COC(=CC=O)OO(1873)'],
    products = ['O=CC=C(OO)OCO(1849)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.80515e+06,'m^3/(mol*s)'), n=0.314888, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O in family R_Recombination."""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(5)', '[O]OC(=CC=O)OCO(1874)'],
    products = ['O=CC=C(OO)OCO(1849)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(5250.69,'m^3/(mol*s)'), n=1.27262, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O in family R_Recombination."""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(5)', 'O=CC=C(OO)O[CH]O(1875)'],
    products = ['O=CC=C(OO)OCO(1849)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.88213e+06,'m^3/(mol*s)'), n=0.314663, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.379272271586, var=0.88677526262, Tref=1000.0, N=15, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O in family R_Recombination."""),
)

reaction(
    label = 'reaction13',
    reactants = ['HCO(16)', '[CH]=C(OO)OCO(1876)'],
    products = ['O=CC=C(OO)OCO(1849)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.81e+07,'m^3/(mol*s)'), n=-1.42812e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_N-Sp-3R!H-2R_3R!H->O_Ext-1C-R_4R!H->C',), comment="""Estimated from node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_N-Sp-3R!H-2R_3R!H->O_Ext-1C-R_4R!H->C in family R_Recombination."""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(5)', '[O]C=C=C(OO)OCO(1877)'],
    products = ['O=CC=C(OO)OCO(1849)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.52804e+49,'m^3/(mol*s)'), n=-12.7885, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R_Ext-3CS-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R_Ext-3CS-R in family R_Recombination."""),
)

reaction(
    label = 'reaction15',
    reactants = ['O=C[C]C(OO)OCO(1879)'],
    products = ['O=CC=C(OO)OCO(1849)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(3.95159e+18,'s^-1'), n=-1.97082, Ea=(23.237,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.8305765209480382, var=17.041866923497537, Tref=1000.0, N=4, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root in family Singlet_Carbene_Intra_Disproportionation."""),
)

network(
    label = 'PDepNetwork #1252',
    isomers = [
        'O=CC=C(OO)OCO(1849)',
    ],
    reactants = [
        ('CH2O(20)', 'O=CCC(=O)OO(1774)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1252',
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

