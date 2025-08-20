species(
    label = '[O]OC=C(O)CC(=O)O(820)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {7,S} {14,S}
2  O u0 p2 c0 {8,S} {13,S}
3  O u0 p2 c0 {5,S} {9,S}
4  O u0 p2 c0 {8,D}
5  O u1 p2 c0 {3,S}
6  C u0 p0 c0 {7,S} {8,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {6,S} {9,D}
8  C u0 p0 c0 {2,S} {4,D} {6,S}
9  C u0 p0 c0 {3,S} {7,D} {12,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {1,S}
"""),
    E0 = (-505.677,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,200,800,1000,1200,1400,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.265204,0.0807559,-7.8361e-05,3.59229e-08,-6.34439e-12,-60654.3,33.0158], Tmin=(100,'K'), Tmax=(1388.57,'K')), NASAPolynomial(coeffs=[23.0231,0.01367,-5.89152e-06,1.12954e-09,-8.01296e-14,-67121.8,-86.9818], Tmin=(1388.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-505.677,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-O2d)H) + group(O2s-OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-OdCsOs) + group(Cds-CdsOsH) + radical(ROOJ)"""),
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
    label = 'C=C(O)CO[O](645)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {5,S} {11,S}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
5  C u0 p0 c0 {2,S} {4,S} {6,D}
6  C u0 p0 c0 {5,D} {9,S} {10,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {2,S}
"""),
    E0 = (-138.061,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2950,3100,1380,975,1025,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.532508,'amu*angstrom^2'), symmetry=1, barrier=(12.2434,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.532647,'amu*angstrom^2'), symmetry=1, barrier=(12.2466,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.533286,'amu*angstrom^2'), symmetry=1, barrier=(12.2613,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (89.07,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4075.04,'J/mol'), sigma=(5.36414,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=636.51 K, Pc=59.91 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23404,0.0562499,-6.01604e-05,3.28655e-08,-7.0062e-12,-16501.4,22.4257], Tmin=(100,'K'), Tmax=(1153.59,'K')), NASAPolynomial(coeffs=[13.5536,0.0135331,-4.61689e-06,7.66992e-10,-5.00784e-14,-19343.8,-38.7696], Tmin=(1153.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-138.061,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(ROOJ)"""),
)

species(
    label = 'CC(O)=CO[O](644)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {5,S} {11,S}
2  O u0 p2 c0 {3,S} {6,S}
3  O u1 p2 c0 {2,S}
4  C u0 p0 c0 {5,S} {7,S} {8,S} {9,S}
5  C u0 p0 c0 {1,S} {4,S} {6,D}
6  C u0 p0 c0 {2,S} {5,D} {10,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {1,S}
"""),
    E0 = (-140.728,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.632468,'amu*angstrom^2'), symmetry=1, barrier=(14.5417,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.634013,'amu*angstrom^2'), symmetry=1, barrier=(14.5772,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.634144,'amu*angstrom^2'), symmetry=1, barrier=(14.5802,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (89.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.830754,0.0593933,-6.44122e-05,3.40449e-08,-6.78263e-12,-16802.7,23.2584], Tmin=(100,'K'), Tmax=(1373.82,'K')), NASAPolynomial(coeffs=[16.4657,0.00766352,-1.15401e-06,5.91453e-11,3.72472e-16,-20512.9,-55.005], Tmin=(1373.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-140.728,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(ROOJ)"""),
)

species(
    label = 'H2O(11)',
    structure = adjacencyList("""1 O u0 p2 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
"""),
    E0 = (-251.721,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1591.57,3573.8,3573.81],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (18.0153,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(6727.26,'J/mol'), sigma=(2.641,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.76,'angstroms^3'), rotrelaxcollnum=4.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.99882,-0.000554834,2.76775e-06,-1.55666e-09,3.02332e-13,-30274.6,-0.0308953], Tmin=(100,'K'), Tmax=(1281.43,'K')), NASAPolynomial(coeffs=[3.19561,0.0019524,-1.67116e-07,-2.97944e-11,4.45144e-15,-30068.7,4.04333], Tmin=(1281.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-251.721,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""H2O""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = '[O]OC=C(O)C=C=O(937)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {5,S} {11,S}
2  O u0 p2 c0 {3,S} {7,S}
3  O u1 p2 c0 {2,S}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {1,S} {6,S} {7,D}
6  C u0 p0 c0 {5,S} {8,D} {9,S}
7  C u0 p0 c0 {2,S} {5,D} {10,S}
8  C u0 p0 c0 {4,D} {6,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {1,S}
"""),
    E0 = (-125.624,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,492.5,1135,1000,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,2120,512.5,787.5,180],'cm^-1')),
        HinderedRotor(inertia=(0.756629,'amu*angstrom^2'), symmetry=1, barrier=(17.3964,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.760992,'amu*angstrom^2'), symmetry=1, barrier=(17.4967,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.756366,'amu*angstrom^2'), symmetry=1, barrier=(17.3903,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.313785,0.0792491,-0.000106408,6.70108e-08,-1.55858e-11,-14974.6,26.9709], Tmin=(100,'K'), Tmax=(886.619,'K')), NASAPolynomial(coeffs=[17.9107,0.00885938,-2.54602e-06,3.6344e-10,-2.12944e-14,-18448.7,-57.8014], Tmin=(886.619,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-125.624,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + missing(O2d-Cdd) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(ROOJ)"""),
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
    label = '[O]OC=C(O)O(938)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {5,S} {8,S}
2 O u0 p2 c0 {5,S} {9,S}
3 O u0 p2 c0 {4,S} {6,S}
4 O u1 p2 c0 {3,S}
5 C u0 p0 c0 {1,S} {2,S} {6,D}
6 C u0 p0 c0 {3,S} {5,D} {7,S}
7 H u0 p0 c0 {6,S}
8 H u0 p0 c0 {1,S}
9 H u0 p0 c0 {2,S}
"""),
    E0 = (-256.073,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,492.5,1135,1000,350,440,435,1725,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.792583,'amu*angstrom^2'), symmetry=1, barrier=(18.223,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.789555,'amu*angstrom^2'), symmetry=1, barrier=(18.1534,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.78976,'amu*angstrom^2'), symmetry=1, barrier=(18.1581,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (91.0428,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.561023,0.0666093,-8.99863e-05,5.44487e-08,-1.20098e-11,-30666.8,22], Tmin=(100,'K'), Tmax=(1287.14,'K')), NASAPolynomial(coeffs=[19.0416,-0.00220892,3.48245e-06,-8.35886e-10,6.29717e-14,-34481,-68.159], Tmin=(1287.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-256.073,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(ROOJ)"""),
)

species(
    label = 'C=C(O)C(O[O])C(=O)O(813)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {7,S} {14,S}
3  O u0 p2 c0 {8,S} {13,S}
4  O u0 p2 c0 {8,D}
5  O u1 p2 c0 {1,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {10,S}
7  C u0 p0 c0 {2,S} {6,S} {9,D}
8  C u0 p0 c0 {3,S} {4,D} {6,S}
9  C u0 p0 c0 {7,D} {11,S} {12,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {2,S}
"""),
    E0 = (-527.073,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (133.079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.213631,0.0829353,-8.62432e-05,4.36394e-08,-8.53838e-12,-63232.3,34.2697], Tmin=(100,'K'), Tmax=(1257.14,'K')), NASAPolynomial(coeffs=[20.9269,0.0156697,-5.98252e-06,1.07656e-09,-7.41279e-14,-68547.6,-72.5587], Tmin=(1257.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-527.073,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-O2d)H) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsOs) + group(Cds-OdCsOs) + group(Cds-CdsHH) + radical(ROOJ)"""),
)

species(
    label = 'C=C(O)OC(O)=CO[O](939)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {6,S} {7,S}
2  O u0 p2 c0 {6,S} {14,S}
3  O u0 p2 c0 {7,S} {13,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u1 p2 c0 {4,S}
6  C u0 p0 c0 {1,S} {2,S} {8,D}
7  C u0 p0 c0 {1,S} {3,S} {9,D}
8  C u0 p0 c0 {4,S} {6,D} {10,S}
9  C u0 p0 c0 {7,D} {11,S} {12,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {2,S}
"""),
    E0 = (-276.753,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,492.5,1135,1000,325,375,415,465,420,450,1700,1750,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.766128,0.107733,-0.000157887,1.11159e-07,-2.94259e-11,-33116.6,32.3904], Tmin=(100,'K'), Tmax=(799.776,'K')), NASAPolynomial(coeffs=[19.5019,0.0166893,-6.49804e-06,1.10851e-09,-7.13403e-14,-36688.8,-62.9272], Tmin=(799.776,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-276.753,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ)"""),
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
    label = '[O]C=C(O)CC(=O)O(940)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {6,S} {12,S}
2  O u0 p2 c0 {7,S} {13,S}
3  O u0 p2 c0 {7,D}
4  O u1 p2 c0 {8,S}
5  C u0 p0 c0 {6,S} {7,S} {9,S} {10,S}
6  C u0 p0 c0 {1,S} {5,S} {8,D}
7  C u0 p0 c0 {2,S} {3,D} {5,S}
8  C u0 p0 c0 {4,S} {6,D} {11,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {2,S}
"""),
    E0 = (-640.714,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.459143,0.0589234,-1.12977e-05,-4.05626e-08,2.24556e-11,-76915.6,27.0094], Tmin=(100,'K'), Tmax=(978.132,'K')), NASAPolynomial(coeffs=[23.6561,0.00885066,-3.19525e-06,7.29431e-10,-6.34114e-14,-83596.1,-95.3419], Tmin=(978.132,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-640.714,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-O2d)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-OdCsOs) + group(Cds-CdsOsH) + radical(C=COJ)"""),
)

species(
    label = 'O=C(O)CC1(O)[CH]OO1(941)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {3,S} {6,S}
2  O u0 p2 c0 {6,S} {13,S}
3  O u0 p2 c0 {1,S} {8,S}
4  O u0 p2 c0 {9,S} {14,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7  C u0 p0 c0 {6,S} {9,S} {10,S} {11,S}
8  C u1 p0 c0 {3,S} {6,S} {12,S}
9  C u0 p0 c0 {4,S} {5,D} {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {4,S}
"""),
    E0 = (-449.12,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (133.079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.465114,0.0859554,-8.84989e-05,4.4178e-08,-8.4955e-12,-53845.5,28.6477], Tmin=(100,'K'), Tmax=(1282.5,'K')), NASAPolynomial(coeffs=[22.3513,0.0147928,-5.26742e-06,9.12609e-10,-6.16517e-14,-59697.9,-87.1053], Tmin=(1282.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-449.12,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(Cs-CsCsOsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsOs) + ring(12dioxetane) + radical(CCsJOO)"""),
)

species(
    label = 'O[C]1CC(O)=COOO1(942)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {5,S} {8,S}
2  O u0 p2 c0 {5,S} {9,S}
3  O u0 p2 c0 {7,S} {14,S}
4  O u0 p2 c0 {8,S} {13,S}
5  O u0 p2 c0 {1,S} {2,S}
6  C u0 p0 c0 {7,S} {8,S} {10,S} {11,S}
7  C u0 p0 c0 {3,S} {6,S} {9,D}
8  C u1 p0 c0 {1,S} {4,S} {6,S}
9  C u0 p0 c0 {2,S} {7,D} {12,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {3,S}
"""),
    E0 = (-184.283,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (133.079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.465896,0.0143521,0.000210493,-3.45853e-07,1.52345e-10,-21977.1,33.5434], Tmin=(100,'K'), Tmax=(899.766,'K')), NASAPolynomial(coeffs=[51.861,-0.0463418,3.19554e-05,-6.31458e-09,4.19174e-13,-38017.8,-246.723], Tmin=(899.766,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-184.283,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(Cycloheptane) + radical(Cs_P)"""),
)

species(
    label = 'O=C(O)C[C](O)C1OO1(943)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {2,S} {6,S}
2  O u0 p2 c0 {1,S} {6,S}
3  O u0 p2 c0 {8,S} {14,S}
4  O u0 p2 c0 {9,S} {13,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {1,S} {2,S} {8,S} {12,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {11,S}
8  C u1 p0 c0 {3,S} {6,S} {7,S}
9  C u0 p0 c0 {4,S} {5,D} {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {3,S}
"""),
    E0 = (-454.48,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (133.079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.173577,0.083109,-9.5144e-05,5.67033e-08,-1.33741e-11,-54522.5,31.8751], Tmin=(100,'K'), Tmax=(1036.44,'K')), NASAPolynomial(coeffs=[15.4705,0.0240741,-9.70735e-06,1.74953e-09,-1.19018e-13,-57693.5,-42.4717], Tmin=(1036.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-454.48,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(Cs-CsCsOsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsOs) + ring(dioxirane) + radical(C2CsJOH)"""),
)

species(
    label = '[O]C1(O)CC(O)=COO1(944)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {2,S} {7,S}
2  O u0 p2 c0 {1,S} {9,S}
3  O u0 p2 c0 {7,S} {13,S}
4  O u0 p2 c0 {8,S} {14,S}
5  O u1 p2 c0 {7,S}
6  C u0 p0 c0 {7,S} {8,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {3,S} {5,S} {6,S}
8  C u0 p0 c0 {4,S} {6,S} {9,D}
9  C u0 p0 c0 {2,S} {8,D} {12,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
"""),
    E0 = (-407.243,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (133.079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.23895,0.076735,-7.5858e-05,3.76222e-08,-7.07388e-12,-48813.3,29.5008], Tmin=(100,'K'), Tmax=(1477.47,'K')), NASAPolynomial(coeffs=[19.2972,0.0142838,-2.74822e-06,2.53799e-10,-9.76204e-15,-53542.7,-68.8441], Tmin=(1477.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-407.243,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsOsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(34dihydro12dioxin) + radical(CCOJ)"""),
)

species(
    label = '[O]OCC(=O)CC(=O)O(658)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {5,S} {7,S}
2  O u0 p2 c0 {9,S} {14,S}
3  O u0 p2 c0 {8,D}
4  O u0 p2 c0 {9,D}
5  O u1 p2 c0 {1,S}
6  C u0 p0 c0 {8,S} {9,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {8,S} {12,S} {13,S}
8  C u0 p0 c0 {3,D} {6,S} {7,S}
9  C u0 p0 c0 {2,S} {4,D} {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {2,S}
"""),
    E0 = (-515.156,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (133.079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.511338,0.082969,-0.000108001,8.02897e-08,-2.46998e-11,-61839,31.2643], Tmin=(100,'K'), Tmax=(786.006,'K')), NASAPolynomial(coeffs=[9.84589,0.0354634,-1.73382e-05,3.38916e-09,-2.39473e-13,-63306.4,-11.5215], Tmin=(786.006,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-515.156,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-O2d)H) + group(O2s-OsH) + group(Cs-(Cds-O2d)(Cds-O2d)HH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsCs) + group(Cds-OdCsOs) + radical(ROOJ)"""),
)

species(
    label = '[O]OC=C(O)C=C(O)O(945)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {6,S} {14,S}
2  O u0 p2 c0 {8,S} {13,S}
3  O u0 p2 c0 {8,S} {12,S}
4  O u0 p2 c0 {5,S} {9,S}
5  O u1 p2 c0 {4,S}
6  C u0 p0 c0 {1,S} {7,S} {9,D}
7  C u0 p0 c0 {6,S} {8,D} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {7,D}
9  C u0 p0 c0 {4,S} {6,D} {11,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {1,S}
"""),
    E0 = (-418.194,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3615,3650,1210,1277.5,1345,900,1000,1100,492.5,1135,1000,325,375,415,465,420,450,1700,1750,2995,3025,975,1000,1300,1375,400,500,1630,1680,180],'cm^-1')),
        HinderedRotor(inertia=(0.865763,'amu*angstrom^2'), symmetry=1, barrier=(19.9056,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.861815,'amu*angstrom^2'), symmetry=1, barrier=(19.8148,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.85991,'amu*angstrom^2'), symmetry=1, barrier=(19.771,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.861944,'amu*angstrom^2'), symmetry=1, barrier=(19.8178,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.860021,'amu*angstrom^2'), symmetry=1, barrier=(19.7736,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (133.079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.03761,0.116968,-0.000157344,9.46619e-08,-2.08082e-11,-50065.9,33.2003], Tmin=(100,'K'), Tmax=(1282.81,'K')), NASAPolynomial(coeffs=[30.7965,-0.00395975,5.7438e-06,-1.36412e-09,1.02205e-13,-56964,-127.435], Tmin=(1282.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-418.194,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(ROOJ)"""),
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
    label = '[CH]=C(O)CC(=O)O(946)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {5,S} {10,S}
2  O u0 p2 c0 {6,S} {12,S}
3  O u0 p2 c0 {6,D}
4  C u0 p0 c0 {5,S} {6,S} {8,S} {9,S}
5  C u0 p0 c0 {1,S} {4,S} {7,D}
6  C u0 p0 c0 {2,S} {3,D} {4,S}
7  C u1 p0 c0 {5,D} {11,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {2,S}
"""),
    E0 = (-326.288,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3120,650,792.5,1650,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.996038,0.0537704,-2.69227e-05,-8.71587e-09,7.79555e-12,-39124.5,24.2195], Tmin=(100,'K'), Tmax=(1046.75,'K')), NASAPolynomial(coeffs=[17.8484,0.0145344,-6.75564e-06,1.4053e-09,-1.06689e-13,-44031.1,-64.4379], Tmin=(1046.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-326.288,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-OdCsOs) + group(Cds-CdsHH) + radical(Cds_P)"""),
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
    label = 'C=C(O)O(947)',
    structure = adjacencyList("""1 O u0 p2 c0 {3,S} {7,S}
2 O u0 p2 c0 {3,S} {8,S}
3 C u0 p0 c0 {1,S} {2,S} {4,D}
4 C u0 p0 c0 {3,D} {5,S} {6,S}
5 H u0 p0 c0 {4,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {1,S}
8 H u0 p0 c0 {2,S}
"""),
    E0 = (-323.78,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (60.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13794,0.0445303,-4.90871e-05,2.46123e-08,-4.41263e-12,-38822.8,14.5979], Tmin=(100,'K'), Tmax=(1681.79,'K')), NASAPolynomial(coeffs=[14.232,-3.33456e-05,2.62933e-06,-6.33079e-10,4.54539e-14,-41329.1,-49.7376], Tmin=(1681.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-323.78,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-CdsCsCs) + group(Cds-CdsHH)"""),
)

species(
    label = '[O]OC#CO(948)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {4,S} {6,S}
2 O u0 p2 c0 {3,S} {5,S}
3 O u1 p2 c0 {2,S}
4 C u0 p0 c0 {1,S} {5,T}
5 C u0 p0 c0 {2,S} {4,T}
6 H u0 p0 c0 {1,S}
"""),
    E0 = (199.907,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (73.0274,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.49787,0.0359027,-5.74249e-05,4.8708e-08,-1.63977e-11,24094.6,16.3582], Tmin=(100,'K'), Tmax=(790.068,'K')), NASAPolynomial(coeffs=[6.88575,0.0111936,-5.77806e-06,1.13252e-09,-7.92127e-14,23479.1,-3.28424], Tmin=(790.068,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(199.907,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(124.717,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCt) + group(O2s-CtH) + group(O2s-OsH) + group(Ct-CtOs) + group(Ct-CtOs) + radical(ROOJ)"""),
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
    label = 'O=C=C(O)CC(=O)O(949)',
    structure = adjacencyList("""1  O u0 p2 c0 {6,S} {12,S}
2  O u0 p2 c0 {7,S} {11,S}
3  O u0 p2 c0 {7,D}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {6,S} {7,S} {9,S} {10,S}
6  C u0 p0 c0 {1,S} {5,S} {8,D}
7  C u0 p0 c0 {2,S} {3,D} {5,S}
8  C u0 p0 c0 {4,D} {6,D}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {1,S}
"""),
    E0 = (-608.984,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.323458,0.071665,-7.65387e-05,3.98537e-08,-7.98205e-12,-73103.6,27.0204], Tmin=(100,'K'), Tmax=(1234.38,'K')), NASAPolynomial(coeffs=[18.7516,0.0119467,-3.96737e-06,6.57842e-10,-4.3403e-14,-77652.9,-65.7643], Tmin=(1234.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-608.984,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-O2d)H) + missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-(Cdd-O2d)CsOs) + group(Cds-OdCsOs) + missing(Cdd-CdO2d)"""),
)

species(
    label = 'O=CC(O)=CC(=O)O(950)',
    structure = adjacencyList("""1  O u0 p2 c0 {5,S} {12,S}
2  O u0 p2 c0 {7,S} {11,S}
3  O u0 p2 c0 {7,D}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {1,S} {6,D} {8,S}
6  C u0 p0 c0 {5,D} {7,S} {9,S}
7  C u0 p0 c0 {2,S} {3,D} {6,S}
8  C u0 p0 c0 {4,D} {5,S} {10,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {1,S}
"""),
    E0 = (-617.635,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3204,0.0636719,-8.45529e-05,5.60319e-08,-1.47518e-11,-74191.8,21.548], Tmin=(100,'K'), Tmax=(924.693,'K')), NASAPolynomial(coeffs=[11.9981,0.0174842,-9.63106e-06,2.01766e-09,-1.48914e-13,-76166.5,-29.1298], Tmin=(924.693,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-617.635,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-O2d)H) + group(Cds-Cds(Cds-O2d)O2s) + missing(Cd-COCdH0) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-O2d(Cds-Cds)H)"""),
)

species(
    label = '[O]C(=O)CC(O)=COO(951)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {3,S} {8,S}
2  O u0 p2 c0 {7,S} {13,S}
3  O u0 p2 c0 {1,S} {14,S}
4  O u1 p2 c0 {9,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {7,S} {9,S} {10,S} {11,S}
7  C u0 p0 c0 {2,S} {6,S} {8,D}
8  C u0 p0 c0 {1,S} {7,D} {12,S}
9  C u0 p0 c0 {4,S} {5,D} {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {3,S}
"""),
    E0 = (-431.976,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,200,800,960,1120,1280,1440,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.271595,0.0765206,-6.84137e-05,2.89368e-08,-4.82963e-12,-51816.3,30.4899], Tmin=(100,'K'), Tmax=(1429.13,'K')), NASAPolynomial(coeffs=[19.7418,0.0220246,-1.12146e-05,2.25405e-09,-1.61923e-13,-57381.3,-70.3947], Tmin=(1429.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-431.976,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-O2d)H) + group(O2s-OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-OdCsOs) + group(Cds-CdsOsH) + radical(CCOJ)"""),
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
    E0 = (-176.021,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (93.1094,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (104.579,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (325.113,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-68.8206,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (82.1839,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-84.8185,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-33.2064,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (128.608,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-124.854,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-52.5273,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-57.5529,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-1.12994,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-22.0189,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (60.9609,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (228.456,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-16.919,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-130.309,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-46.8721,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]OC=C(O)CC(=O)O(820)'],
    products = ['CO2(4)', 'C=C(O)CO[O](645)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(2.46545e+14,'s^-1'), n=0.20628, Ea=(16.765,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(3000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R!H->C',), comment="""Estimated from node Root_N-4R!H->C in family Retroene."""),
)

reaction(
    label = 'reaction2',
    reactants = ['CO2(4)', 'CC(O)=CO[O](644)'],
    products = ['[O]OC=C(O)CC(=O)O(820)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(0.01228,'m^3/(mol*s)'), n=2.86279, Ea=(324.078,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=5.3453611489250746e-05, var=2.530939774000289, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-5CbCdCsHN->Cs_N-4CbCdCsHN->H_Ext-4CsN-R',), comment="""Estimated from node Root_N-5CbCdCsHN->Cs_N-4CbCdCsHN->H_Ext-4CsN-R in family 1,3_Insertion_CO2.
Multiplied by reaction path degeneracy 6.0"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H2O(11)', '[O]OC=C(O)C=C=O(937)'],
    products = ['[O]OC=C(O)CC(=O)O(820)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(0.000454348,'m^3/(mol*s)'), n=2.96333, Ea=(169.034,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_Cdd;H_OH] for rate rule [cco_HDe;H_OH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CH2CO(27)', '[O]OC=C(O)O(938)'],
    products = ['[O]OC=C(O)CC(=O)O(820)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.58e-11,'m^3/(mol*s)'), n=3.97, Ea=(329.281,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [doublebond;R_OH] for rate rule [cco_2H;Cd_sec_OH]
Euclidian distance = 2.8284271247461903
Multiplied by reaction path degeneracy 2.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]OC=C(O)CC(=O)O(820)'],
    products = ['C=C(O)C(O[O])C(=O)O(813)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(7.50839e+11,'s^-1'), n=0.136197, Ea=(123.966,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.01703312323488117, var=1.3519057502930525, Tref=1000.0, N=4, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R in family 1,3_sigmatropic_rearrangement."""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=C(O)OC(O)=CO[O](939)'],
    products = ['[O]OC=C(O)CC(=O)O(820)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(7.50839e+11,'s^-1'), n=0.136197, Ea=(46.0462,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.01703312323488117, var=1.3519057502930525, Tref=1000.0, N=4, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R in family 1,3_sigmatropic_rearrangement."""),
)

reaction(
    label = 'reaction7',
    reactants = ['O(10)', '[O]C=C(O)CC(=O)O(940)'],
    products = ['[O]OC=C(O)CC(=O)O(820)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(15399.6,'m^3/(mol*s)'), n=0.932885, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [O_sec_rad;O_birad] for rate rule [O_rad/OneDe;O_birad]
Euclidian distance = 1.0
family: Birad_R_Recombination
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]OC=C(O)CC(=O)O(820)'],
    products = ['O=C(O)CC1(O)[CH]OO1(941)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(6.12241e+07,'s^-1'), n=0.973219, Ea=(159.58,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.004397713362318132, var=0.35833417126470263, Tref=1000.0, N=6, data_mean=0.0, correlation='Backbone1_N-2R!H-inRing_Ext-4R!H-R_Ext-4R!H-R_Ext-5R!H-R',), comment="""Estimated from node Backbone1_N-2R!H-inRing_Ext-4R!H-R_Ext-4R!H-R_Ext-5R!H-R in family Intra_R_Add_Endocyclic."""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]OC=C(O)CC(=O)O(820)'],
    products = ['O[C]1CC(O)=COOO1(942)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(5.649e+12,'s^-1'), n=0.287, Ea=(321.394,'kJ/mol'), T0=(1,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Backbone4_N-1R!H-inRing_N-7R!H->C',), comment="""Estimated from node Backbone4_N-1R!H-inRing_N-7R!H->C in family Intra_R_Add_Endocyclic.
Ea raised from 314.4 to 321.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]OC=C(O)CC(=O)O(820)'],
    products = ['O=C(O)C[C](O)C1OO1(943)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.36539e+09,'s^-1'), n=0.84129, Ea=(67.9324,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.031430614527560546, var=6.099077214950256, Tref=1000.0, N=84, data_mean=0.0, correlation='Backbone1',), comment="""Estimated from node Backbone1 in family Intra_R_Add_Exocyclic."""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]OC=C(O)CC(=O)O(820)'],
    products = ['[O]C1(O)CC(O)=COO1(944)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.87342e+07,'s^-1'), n=1.00417, Ea=(140.259,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0654062900291529, var=6.730759255238982, Tref=1000.0, N=75, data_mean=0.0, correlation='Backbone4',), comment="""Estimated from node Backbone4 in family Intra_R_Add_Exocyclic."""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]OC=C(O)CC(=O)O(820)'],
    products = ['[O]OCC(=O)CC(=O)O(658)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3.39196e-35,'s^-1'), n=13.7658, Ea=(135.233,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.005573865076554089, var=7.000467719691818, Tref=1000.0, N=12, data_mean=0.0, correlation='Root_3R!H->C_N-3C-inRing_N-1R!H->N_Ext-3C-R_Ext-5R!H-R',), comment="""Estimated from node Root_3R!H->C_N-3C-inRing_N-1R!H->N_Ext-3C-R_Ext-5R!H-R in family Ketoenol."""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]OC=C(O)C=C(O)O(945)'],
    products = ['[O]OC=C(O)CC(=O)O(820)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.10345e-17,'s^-1'), n=8.84109, Ea=(104.174,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R!H->C_N-3C-inRing_N-1R!H->N_Ext-3C-R_Ext-5R!H-R_Ext-1C-R_N-Sp-6R!H#5R!H_5R!H->C',), comment="""Estimated from node Root_3R!H->C_N-3C-inRing_N-1R!H->N_Ext-3C-R_Ext-5R!H-R_Ext-1C-R_N-Sp-6R!H#5R!H_5R!H->C in family Ketoenol.
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction14',
    reactants = ['O2(2)', '[CH]=C(O)CC(=O)O(946)'],
    products = ['[O]OC=C(O)CC(=O)O(820)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(6.90212e-07,'m^3/(mol*s)'), n=3.72998, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_N-Sp-4R!H-2R',), comment="""Estimated from node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_N-Sp-4R!H-2R in family R_Recombination.
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[O]OC=C(O)CC(=O)O(820)'],
    products = ['[O]OC=C=O(111)', 'C=C(O)O(947)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.46545e+14,'s^-1'), n=0.20628, Ea=(253.747,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(3000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R!H->C',), comment="""Estimated from node Root_N-4R!H->C in family Retroene."""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O]OC=C(O)CC(=O)O(820)'],
    products = ['[O]OC#CO(948)', 'C=C(O)O(947)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.78363e+11,'s^-1'), n=0.169307, Ea=(421.243,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_4R!H->C_N-1R!H->O_2R!H->C_5R!H->O_Ext-2C-R_N-7R!H->C',), comment="""Estimated from node Root_4R!H->C_N-1R!H->O_2R!H->C_5R!H->O_Ext-2C-R_N-7R!H->C in family Retroene."""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O]OC=C(O)CC(=O)O(820)'],
    products = ['OH(9)', 'O=C=C(O)CC(=O)O(949)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(3.36833e+07,'s^-1'), n=1.41, Ea=(175.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Y_rad_out;Cd_H_out_doubleC] for rate rule [R3H_SS_O;O_rad_out;Cd_H_out_doubleC]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[O]OC=C(O)CC(=O)O(820)'],
    products = ['OH(9)', 'O=CC(O)=CC(=O)O(950)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(59826.2,'s^-1'), n=1.86494, Ea=(62.4776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R5H_SSMS;O_rad_out;Cs_H_out_H/CO]
Euclidian distance = 2.449489742783178
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[O]C(=O)CC(O)=COO(951)'],
    products = ['[O]OC=C(O)CC(=O)O(820)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.25419e+06,'s^-1'), n=1.03067, Ea=(72.2138,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7H;O_rad_out;XH_out] for rate rule [R7H;O_rad_out;O_H_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #748',
    isomers = [
        '[O]OC=C(O)CC(=O)O(820)',
    ],
    reactants = [
        ('CO2(4)', 'C=C(O)CO[O](645)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #748',
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

