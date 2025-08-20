species(
    label = '[O]OC(=CCO)OO[OH+][O-](1694)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p1 c+1 {3,S} {6,S} {14,S}
2  O u0 p2 c0 {3,S} {10,S}
3  O u0 p2 c0 {1,S} {2,S}
4  O u0 p2 c0 {8,S} {15,S}
5  O u0 p2 c0 {7,S} {10,S}
6  O u0 p3 c-1 {1,S}
7  O u1 p2 c0 {5,S}
8  C u0 p0 c0 {4,S} {9,S} {11,S} {12,S}
9  C u0 p0 c0 {8,S} {10,D} {13,S}
10 C u0 p0 c0 {2,S} {5,S} {9,D}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {4,S}
"""),
    E0 = (-226.103,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,350,440,435,1725,412.075,2644.08,4000,4000,4000,4000,4000,4000,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.513312,'amu*angstrom^2'), symmetry=1, barrier=(15.1281,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.513312,'amu*angstrom^2'), symmetry=1, barrier=(15.1281,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.513312,'amu*angstrom^2'), symmetry=1, barrier=(15.1281,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.513312,'amu*angstrom^2'), symmetry=1, barrier=(15.1281,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.513312,'amu*angstrom^2'), symmetry=1, barrier=(15.1281,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.513312,'amu*angstrom^2'), symmetry=1, barrier=(15.1281,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (153.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.239334,0.105376,-0.000219503,2.26947e-07,-8.5338e-11,-27080.2,26.7474], Tmin=(100,'K'), Tmax=(867.377,'K')), NASAPolynomial(coeffs=[0.779885,0.0480373,-2.54979e-05,4.93419e-09,-3.36803e-13,-25110.9,36.1094], Tmin=(867.377,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-226.103,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-CsH) + group(O2s-OsH) + missing(O0sc-O4sc) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsCs) + radical(ROOJ)"""),
)

species(
    label = 'O3(14)',
    structure = adjacencyList("""1 O u0 p1 c+1 {2,S} {3,D}
2 O u0 p3 c-1 {1,S}
3 O u0 p2 c0 {1,D}
"""),
    E0 = (132.542,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([923.818,923.85,924.126],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (47.9982,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1971.36,'J/mol'), sigma=(5.118e-10,'m'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with fixed Lennard Jones Parameters. This is the fallback method! Try improving transport databases!"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.46261,0.00958278,-7.08736e-06,1.36337e-09,2.96965e-13,16061.5,12.1419], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.42937,0.00182038,-7.70561e-07,1.49929e-10,-1.07556e-14,15235.3,-3.26639], Tmin=(1000,'K'), Tmax=(5000,'K'))], Tmin=(200,'K'), Tmax=(5000,'K'), E0=(132.542,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""O3""", comment="""Thermo library: primaryThermoLibrary"""),
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
    label = 'C=CC(O)(O[O])OO[OH+][O-](1747)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p1 c+1 {5,S} {6,S} {14,S}
2  O u0 p2 c0 {5,S} {8,S}
3  O u0 p2 c0 {8,S} {15,S}
4  O u0 p2 c0 {7,S} {8,S}
5  O u0 p2 c0 {1,S} {2,S}
6  O u0 p3 c-1 {1,S}
7  O u1 p2 c0 {4,S}
8  C u0 p0 c0 {2,S} {3,S} {4,S} {9,S}
9  C u0 p0 c0 {8,S} {10,D} {11,S}
10 C u0 p0 c0 {9,D} {12,S} {13,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {10,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {3,S}
"""),
    E0 = (-353.454,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (153.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0468749,0.10352,-0.000195623,1.90313e-07,-6.92981e-11,-42384.3,25.1964], Tmin=(100,'K'), Tmax=(859.286,'K')), NASAPolynomial(coeffs=[5.29201,0.0414191,-2.14357e-05,4.13476e-09,-2.82956e-13,-41894.4,8.78245], Tmin=(859.286,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-353.454,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-OsH) + missing(O0sc-O4sc) + group(Cs-CsOsOsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(=O)C(CO)O[OH+][O-](1748)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p1 c+1 {2,S} {6,S} {14,S}
2  O u0 p2 c0 {1,S} {8,S}
3  O u0 p2 c0 {9,S} {15,S}
4  O u0 p2 c0 {7,S} {10,S}
5  O u0 p2 c0 {10,D}
6  O u0 p3 c-1 {1,S}
7  O u1 p2 c0 {4,S}
8  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
9  C u0 p0 c0 {3,S} {8,S} {12,S} {13,S}
10 C u0 p0 c0 {4,S} {5,D} {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {3,S}
"""),
    E0 = (-545.575,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (153.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0568412,0.096336,-0.000155522,1.29527e-07,-4.18453e-11,-65478.2,25.779], Tmin=(100,'K'), Tmax=(855.33,'K')), NASAPolynomial(coeffs=[12.6647,0.0256246,-1.18405e-05,2.20332e-09,-1.48657e-13,-67244,-31.2081], Tmin=(855.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-545.575,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-CsH) + group(O2s-O2s(Cds-O2d)) + group(O2s-(Os-CdOd)H) + missing(O0sc-O4sc) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsHH) + group(Cds-OdCsOs) + radical(C(=O)OOJ)"""),
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
    label = 'O=C([CH]CO)OO[OH+][O-](1749)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p1 c+1 {3,S} {6,S} {13,S}
2  O u0 p2 c0 {3,S} {9,S}
3  O u0 p2 c0 {1,S} {2,S}
4  O u0 p2 c0 {7,S} {14,S}
5  O u0 p2 c0 {9,D}
6  O u0 p3 c-1 {1,S}
7  C u0 p0 c0 {4,S} {8,S} {10,S} {11,S}
8  C u1 p0 c0 {7,S} {9,S} {12,S}
9  C u0 p0 c0 {2,S} {5,D} {8,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {4,S}
"""),
    E0 = (-532.7,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,273.277,1059.98,1059.98,1059.98,1059.98,1059.98,1059.98,1059.98,1059.98,1059.98,1059.98,1059.98,1059.98,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00510022,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00510022,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00510022,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00510022,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00510022,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00510022,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (137.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33821,0.0644763,-7.10098e-05,4.47561e-08,-1.20323e-11,-63978.2,21.5484], Tmin=(100,'K'), Tmax=(879.208,'K')), NASAPolynomial(coeffs=[8.38739,0.032406,-1.62956e-05,3.26903e-09,-2.35675e-13,-65217.8,-11.5524], Tmin=(879.208,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-532.7,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-O2s(Cds-O2d)) + group(O2s-CsH) + group(O2s-(Cds-Cd)(Cds-Cd)) + missing(O0sc-O4sc) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsOs) + radical(CCJCO)"""),
)

species(
    label = '[O][OH+][O-](990)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p1 c+1 {2,S} {3,S} {4,S}
2 O u1 p2 c0 {1,S}
3 O u0 p3 c-1 {1,S}
4 H u0 p0 c0 {1,S}
"""),
    E0 = (27.8905,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([4000,4000,4000,4000,4000,4000],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (49.0062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.6884,-0.00949619,1.54366e-05,-9.39464e-09,1.89056e-12,3324.75,-15.5454], Tmin=(100,'K'), Tmax=(1618.8,'K')), NASAPolynomial(coeffs=[2.24103,0.00323,-2.54436e-06,5.59114e-10,-4.02653e-14,3242.02,-5.26235], Tmin=(1618.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(27.8905,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)(Cds-Cd)) + missing(O0sc-O4sc) + radical(CsOJ)"""),
)

species(
    label = 'OCC=C1OOO1(1750)',
    structure = adjacencyList("""1  O u0 p2 c0 {4,S} {7,S}
2  O u0 p2 c0 {4,S} {7,S}
3  O u0 p2 c0 {5,S} {11,S}
4  O u0 p2 c0 {1,S} {2,S}
5  C u0 p0 c0 {3,S} {6,S} {8,S} {9,S}
6  C u0 p0 c0 {5,S} {7,D} {10,S}
7  C u0 p0 c0 {1,S} {2,S} {6,D}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {3,S}
"""),
    E0 = (28.4617,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.12111,0.0323533,1.90854e-06,-2.2841e-08,9.91982e-12,3498.48,25.2138], Tmin=(100,'K'), Tmax=(1056.43,'K')), NASAPolynomial(coeffs=[10.7468,0.0196966,-8.52243e-06,1.66457e-09,-1.20698e-13,559.762,-22.1567], Tmin=(1056.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.4617,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-OsOs) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsCs) + ring(Cyclobutane)"""),
)

species(
    label = '[O-][OH+]OO[C]1OOC1CO(1751)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p1 c+1 {6,S} {7,S} {14,S}
2  O u0 p2 c0 {3,S} {8,S}
3  O u0 p2 c0 {2,S} {10,S}
4  O u0 p2 c0 {6,S} {10,S}
5  O u0 p2 c0 {9,S} {15,S}
6  O u0 p2 c0 {1,S} {4,S}
7  O u0 p3 c-1 {1,S}
8  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
9  C u0 p0 c0 {5,S} {8,S} {12,S} {13,S}
10 C u1 p0 c0 {3,S} {4,S} {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {5,S}
"""),
    E0 = (-249.598,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (153.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.274788,0.094128,-0.000162241,1.55303e-07,-5.78842e-11,-29897.4,21.4995], Tmin=(100,'K'), Tmax=(811.15,'K')), NASAPolynomial(coeffs=[6.01307,0.0421792,-2.24386e-05,4.45497e-09,-3.13208e-13,-30050.2,-0.18703], Tmin=(811.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-249.598,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)(Cds-Cd)) + missing(O0sc-O4sc) + group(Cs-CsCsOsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + ring(12dioxetane) + radical(Cs_P)"""),
)

species(
    label = '[O-][OH+]OOC1([CH]CO)OO1(1752)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p1 c+1 {5,S} {7,S} {14,S}
2  O u0 p2 c0 {3,S} {8,S}
3  O u0 p2 c0 {2,S} {8,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {1,S} {4,S}
6  O u0 p2 c0 {9,S} {15,S}
7  O u0 p3 c-1 {1,S}
8  C u0 p0 c0 {2,S} {3,S} {4,S} {10,S}
9  C u0 p0 c0 {6,S} {10,S} {11,S} {12,S}
10 C u1 p0 c0 {8,S} {9,S} {13,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (-268.813,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (153.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.192822,0.09044,-0.000138964,1.15373e-07,-3.75599e-11,-32200,25.9069], Tmin=(100,'K'), Tmax=(855.613,'K')), NASAPolynomial(coeffs=[10.7903,0.0298603,-1.34117e-05,2.47114e-09,-1.66319e-13,-33609.5,-21.2065], Tmin=(855.613,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-268.813,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-CsH) + missing(O0sc-O4sc) + group(Cs-CsOsOsOs) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(dioxirane) + radical(CCJCOOH)"""),
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
    label = '[O-][OH+]OO[C]=CCO(1753)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p1 c+1 {3,S} {5,S} {12,S}
2  O u0 p2 c0 {6,S} {13,S}
3  O u0 p2 c0 {1,S} {4,S}
4  O u0 p2 c0 {3,S} {8,S}
5  O u0 p3 c-1 {1,S}
6  C u0 p0 c0 {2,S} {7,S} {9,S} {10,S}
7  C u0 p0 c0 {6,S} {8,D} {11,S}
8  C u1 p0 c0 {4,S} {7,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {2,S}
"""),
    E0 = (-105.722,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,1685,370,180,1290.59,3911.82,3911.82,3911.82,3911.82,3911.82,3911.82,3911.82,3911.82,3911.82],'cm^-1')),
        HinderedRotor(inertia=(0.189389,'amu*angstrom^2'), symmetry=1, barrier=(18.0387,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.189389,'amu*angstrom^2'), symmetry=1, barrier=(18.0387,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.189389,'amu*angstrom^2'), symmetry=1, barrier=(18.0387,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.189389,'amu*angstrom^2'), symmetry=1, barrier=(18.0387,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.189389,'amu*angstrom^2'), symmetry=1, barrier=(18.0387,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (121.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40572,0.0691803,-0.000127231,1.28499e-07,-4.87192e-11,-12633.7,21.102], Tmin=(100,'K'), Tmax=(846.465,'K')), NASAPolynomial(coeffs=[2.52285,0.0358224,-1.83608e-05,3.56555e-09,-2.46279e-13,-11816.9,21.8406], Tmin=(846.465,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-105.722,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)(Cds-Cd)) + missing(O0sc-O4sc) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=CJO)"""),
)

species(
    label = '[O-][OH+]OOC(=[C]CO)OO(1754)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p1 c+1 {3,S} {7,S} {13,S}
2  O u0 p2 c0 {3,S} {9,S}
3  O u0 p2 c0 {1,S} {2,S}
4  O u0 p2 c0 {6,S} {9,S}
5  O u0 p2 c0 {8,S} {14,S}
6  O u0 p2 c0 {4,S} {15,S}
7  O u0 p3 c-1 {1,S}
8  C u0 p0 c0 {5,S} {10,S} {11,S} {12,S}
9  C u0 p0 c0 {2,S} {4,S} {10,D}
10 C u1 p0 c0 {8,S} {9,D}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (-140.266,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,1685,370,180,180,180,180,180,4000,4000,4000,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00223502,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00223502,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00223502,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00223502,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00223502,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00223502,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00223502,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (153.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.065556,0.113673,-0.000240037,2.47667e-07,-9.29776e-11,-16747,27.5559], Tmin=(100,'K'), Tmax=(865.044,'K')), NASAPolynomial(coeffs=[1.29789,0.0489895,-2.6644e-05,5.19453e-09,-3.55793e-13,-14798.6,33.8008], Tmin=(865.044,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-140.266,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-CsH) + group(O2s-OsH) + missing(O0sc-O4sc) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsCs) + radical(Cds_S)"""),
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
    label = 'O=C(C=CO)OO[OH+][O-](1755)',
    structure = adjacencyList("""1  O u0 p1 c+1 {3,S} {6,S} {12,S}
2  O u0 p2 c0 {3,S} {8,S}
3  O u0 p2 c0 {1,S} {2,S}
4  O u0 p2 c0 {9,S} {13,S}
5  O u0 p2 c0 {8,D}
6  O u0 p3 c-1 {1,S}
7  C u0 p0 c0 {8,S} {9,D} {10,S}
8  C u0 p0 c0 {2,S} {5,D} {7,S}
9  C u0 p0 c0 {4,S} {7,D} {11,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {4,S}
"""),
    E0 = (-600.147,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.98625,0.0597032,-6.83696e-05,3.63009e-08,-7.30581e-12,-72066.9,16.9223], Tmin=(100,'K'), Tmax=(1233.5,'K')), NASAPolynomial(coeffs=[17.8247,0.0050987,-1.96677e-06,4.11766e-10,-3.18682e-14,-76220.9,-67.8467], Tmin=(1233.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-600.147,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-O2s(Cds-O2d)) + group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + missing(O0sc-O4sc) + missing(Cd-COCdH0) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsOsH)"""),
)

species(
    label = '[O-][O+]OOC(=CCO)OO(1756)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {4,S} {10,S}
2  O u0 p2 c0 {5,S} {10,S}
3  O u0 p2 c0 {8,S} {14,S}
4  O u0 p2 c0 {1,S} {6,S}
5  O u0 p2 c0 {2,S} {15,S}
6  O u1 p1 c+1 {4,S} {7,S}
7  O u0 p3 c-1 {6,S}
8  C u0 p0 c0 {3,S} {9,S} {11,S} {12,S}
9  C u0 p0 c0 {8,S} {10,D} {13,S}
10 C u0 p0 c0 {1,S} {2,S} {9,D}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {5,S}
"""),
    E0 = (-160.184,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,350,440,435,1725,180,180,180,358.876,3169.49,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00520152,'amu*angstrom^2'), symmetry=1, barrier=(0.119639,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00520152,'amu*angstrom^2'), symmetry=1, barrier=(0.119639,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00520152,'amu*angstrom^2'), symmetry=1, barrier=(0.119639,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00520152,'amu*angstrom^2'), symmetry=1, barrier=(0.119639,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00520152,'amu*angstrom^2'), symmetry=1, barrier=(0.119639,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00520152,'amu*angstrom^2'), symmetry=1, barrier=(0.119639,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00520152,'amu*angstrom^2'), symmetry=1, barrier=(0.119639,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (153.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0742497,0.106807,-0.000214276,2.1824e-07,-8.21021e-11,-19144,26.2542], Tmin=(100,'K'), Tmax=(853.617,'K')), NASAPolynomial(coeffs=[2.45816,0.0475789,-2.57496e-05,5.04916e-09,-3.4861e-13,-17800.1,25.3862], Tmin=(853.617,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-160.184,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-OsH) + missing(O0sc-O4sc) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsCs) + radical(CsOJ)"""),
)

species(
    label = '[O]CC=C(OO)OO[OH+][O-](1757)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p1 c+1 {4,S} {7,S} {14,S}
2  O u0 p2 c0 {4,S} {10,S}
3  O u0 p2 c0 {5,S} {10,S}
4  O u0 p2 c0 {1,S} {2,S}
5  O u0 p2 c0 {3,S} {15,S}
6  O u1 p2 c0 {8,S}
7  O u0 p3 c-1 {1,S}
8  C u0 p0 c0 {6,S} {9,S} {11,S} {12,S}
9  C u0 p0 c0 {8,S} {10,D} {13,S}
10 C u0 p0 c0 {2,S} {3,S} {9,D}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {5,S}
"""),
    E0 = (-152.403,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,350,440,435,1725,378.799,4000,4000,4000,4000,4000,4000,4000,4000,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.193725,'amu*angstrom^2'), symmetry=1, barrier=(29.4381,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.193725,'amu*angstrom^2'), symmetry=1, barrier=(29.4381,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.193725,'amu*angstrom^2'), symmetry=1, barrier=(29.4381,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.193725,'amu*angstrom^2'), symmetry=1, barrier=(29.4381,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.193725,'amu*angstrom^2'), symmetry=1, barrier=(29.4381,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.193725,'amu*angstrom^2'), symmetry=1, barrier=(29.4381,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (153.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.315491,0.106518,-0.000228346,2.44784e-07,-9.48251e-11,-18221.7,26.5736], Tmin=(100,'K'), Tmax=(857.425,'K')), NASAPolynomial(coeffs=[-2.25271,0.0559906,-3.06003e-05,6.00862e-09,-4.14579e-13,-15483.5,51.9676], Tmin=(857.425,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-152.403,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-CsH) + group(O2s-OsH) + missing(O0sc-O4sc) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsCs) + radical(CCOJ)"""),
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
    E0 = (-15.0006,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (139.852,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (68.7471,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-48.5643,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (233.891,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-32.7248,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-30.0353,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (63.1953,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (81.5813,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (16.1622,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (86.8977,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (97.7367,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]OC(=CCO)OO[OH+][O-](1694)'],
    products = ['O3(14)', '[O]OC(=O)CCO(1660)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(2.46545e+14,'s^-1'), n=0.20628, Ea=(33.5637,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(3000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R!H->C',), comment="""Estimated from node Root_N-4R!H->C in family Retroene."""),
)

reaction(
    label = 'reaction2',
    reactants = ['[O]OC(=CCO)OO[OH+][O-](1694)'],
    products = ['C=CC(O)(O[O])OO[OH+][O-](1747)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.30791e+26,'s^-1'), n=-3.67984, Ea=(188.416,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.4163866947181062, var=128.1927359939506, Tref=1000.0, N=20, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root in family 1,3_sigmatropic_rearrangement."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[O]OC(=CCO)OO[OH+][O-](1694)'],
    products = ['[O]OC(=O)C(CO)O[OH+][O-](1748)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.30791e+26,'s^-1'), n=-3.67984, Ea=(117.311,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.4163866947181062, var=128.1927359939506, Tref=1000.0, N=20, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root in family 1,3_sigmatropic_rearrangement."""),
)

reaction(
    label = 'reaction4',
    reactants = ['O(10)', 'O=C([CH]CO)OO[OH+][O-](1749)'],
    products = ['[O]OC(=CCO)OO[OH+][O-](1694)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(15399.6,'m^3/(mol*s)'), n=0.932885, Ea=(63.592,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [O_sec_rad;O_birad] for rate rule [O_rad/OneDe;O_birad]
Euclidian distance = 1.0
family: Birad_R_Recombination
Ea raised from -2.0 to 63.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]OC(=CCO)OO[OH+][O-](1694)'],
    products = ['[O][OH+][O-](990)', 'OCC=C1OOO1(1750)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(3.11355e+11,'s^-1'), n=0, Ea=(282.455,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3OO;Y_rad_intra;OOR]
Euclidian distance = 0
family: Cyclic_Ether_Formation
Ea raised from 281.2 to 282.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]OC(=CCO)OO[OH+][O-](1694)'],
    products = ['[O-][OH+]OO[C]1OOC1CO(1751)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.391e+11,'s^-1'), n=0.287, Ea=(15.8396,'kJ/mol'), T0=(1,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Backbone1_N-2R!H-inRing_Ext-4R!H-R_N-2R!H->C',), comment="""Estimated from node Backbone1_N-2R!H-inRing_Ext-4R!H-R_N-2R!H->C in family Intra_R_Add_Endocyclic."""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]OC(=CCO)OO[OH+][O-](1694)'],
    products = ['[O-][OH+]OOC1([CH]CO)OO1(1752)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(3.36539e+09,'s^-1'), n=0.84129, Ea=(18.529,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.031430614527560546, var=6.099077214950256, Tref=1000.0, N=84, data_mean=0.0, correlation='Backbone1',), comment="""Estimated from node Backbone1 in family Intra_R_Add_Exocyclic."""),
)

reaction(
    label = 'reaction8',
    reactants = ['O2(2)', '[O-][OH+]OO[C]=CCO(1753)'],
    products = ['[O]OC(=CCO)OO[OH+][O-](1694)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.21695e+08,'m^3/(mol*s)'), n=-0.428622, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing',), comment="""Estimated from node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing in family R_Recombination.
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O-][OH+]OOC(=[C]CO)OO(1754)'],
    products = ['[O]OC(=CCO)OO[OH+][O-](1694)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;XH_out] for rate rule [R4H_DSS;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 2.23606797749979
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]OC(=CCO)OO[OH+][O-](1694)'],
    products = ['OH(9)', 'O=C(C=CO)OO[OH+][O-](1755)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(23000,'s^-1'), n=2.11, Ea=(64.7265,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5H;O_rad_out;Cs_H_out_H/NonDeO] for rate rule [R5H_SSMS;O_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O-][O+]OOC(=CCO)OO(1756)'],
    products = ['[O]OC(=CCO)OO[OH+][O-](1694)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(123903,'s^-1'), n=1.46258, Ea=(69.5427,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_SSSSS;O_rad_out;XH_out] for rate rule [R6H_SSSSS;O_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]CC=C(OO)OO[OH+][O-](1757)'],
    products = ['[O]OC(=CCO)OO[OH+][O-](1694)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.33753e+06,'s^-1'), n=1.02312, Ea=(72.6006,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H;O_rad_out;XH_out] for rate rule [R6H_RSMSR;O_rad_out;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #1170',
    isomers = [
        '[O]OC(=CCO)OO[OH+][O-](1694)',
    ],
    reactants = [
        ('O3(14)', '[O]OC(=O)CCO(1660)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1170',
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

