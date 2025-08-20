species(
    label = 'CC(=O)CCOO(522)',
    structure = adjacencyList("""1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {15,S}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {5,S} {7,S} {8,S} {9,S}
5  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
6  C u0 p0 c0 {7,S} {12,S} {13,S} {14,S}
7  C u0 p0 c0 {3,D} {4,S} {6,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {2,S}
"""),
    E0 = (-346.516,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,375,552.5,462.5,1710,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.41809,'amu*angstrom^2'), symmetry=1, barrier=(9.6127,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.041899,'amu*angstrom^2'), symmetry=1, barrier=(0.963341,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00523072,'amu*angstrom^2'), symmetry=1, barrier=(9.62255,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.418482,'amu*angstrom^2'), symmetry=1, barrier=(9.62174,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.972,'amu*angstrom^2'), symmetry=1, barrier=(45.3401,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.105,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4053.87,'J/mol'), sigma=(5.82012,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=633.21 K, Pc=46.66 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45003,0.0654138,-4.72459e-05,-1.77575e-08,3.87688e-11,-41593.8,23.7243], Tmin=(100,'K'), Tmax=(508,'K')), NASAPolynomial(coeffs=[6.42069,0.0411772,-1.96848e-05,3.82034e-09,-2.69079e-13,-42291.1,1.21761], Tmin=(508,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-346.516,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs)"""),
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
    label = 'CC(=O)CC[O](532)',
    structure = adjacencyList("""multiplicity 2
1  O u1 p2 c0 {4,S}
2  O u0 p2 c0 {6,D}
3  C u0 p0 c0 {4,S} {6,S} {7,S} {8,S}
4  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
5  C u0 p0 c0 {6,S} {11,S} {12,S} {13,S}
6  C u0 p0 c0 {2,D} {3,S} {5,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (-192.316,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,375,552.5,462.5,1710,180,180,2734.44],'cm^-1')),
        HinderedRotor(inertia=(0.00503921,'amu*angstrom^2'), symmetry=1, barrier=(7.15181,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0462977,'amu*angstrom^2'), symmetry=1, barrier=(1.06448,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.311993,'amu*angstrom^2'), symmetry=1, barrier=(7.17333,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (87.0972,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3912.27,'J/mol'), sigma=(5.70694,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=611.09 K, Pc=47.76 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.7085,0.0356906,-1.52108e-05,2.18329e-09,-6.02754e-14,-23150.4,15.4207], Tmin=(100,'K'), Tmax=(2776.22,'K')), NASAPolynomial(coeffs=[38.6935,-0.00349493,-1.01554e-07,1.09365e-11,4.24126e-15,-46899.9,-196.873], Tmin=(2776.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-192.316,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + radical(CCOJ)"""),
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
    label = 'CCCOO(528)',
    structure = adjacencyList("""1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {13,S}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
4  C u0 p0 c0 {1,S} {3,S} {8,S} {9,S}
5  C u0 p0 c0 {3,S} {10,S} {11,S} {12,S}
6  H u0 p0 c0 {3,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {2,S}
"""),
    E0 = (-200.879,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,330.429],'cm^-1')),
        HinderedRotor(inertia=(0.0015423,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111517,'amu*angstrom^2'), symmetry=1, barrier=(8.66269,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111544,'amu*angstrom^2'), symmetry=1, barrier=(8.66196,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111724,'amu*angstrom^2'), symmetry=1, barrier=(8.66459,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (76.0945,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72478,0.0452898,-2.97046e-05,1.0046e-08,-1.39809e-12,-24074.5,20.2775], Tmin=(100,'K'), Tmax=(1639.51,'K')), NASAPolynomial(coeffs=[11.0854,0.0224522,-8.81028e-06,1.54981e-09,-1.02549e-13,-27143.8,-29.5099], Tmin=(1639.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-200.879,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), label="""npropylooh""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C=C(C)OCOO(529)',
    structure = adjacencyList("""1  O u0 p2 c0 {4,S} {6,S}
2  O u0 p2 c0 {3,S} {4,S}
3  O u0 p2 c0 {2,S} {15,S}
4  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
5  C u0 p0 c0 {6,S} {10,S} {11,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {7,D}
7  C u0 p0 c0 {6,D} {13,S} {14,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {3,S}
"""),
    E0 = (-310.202,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.105,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0463997,0.0647489,-3.85868e-06,-6.23555e-08,3.43581e-11,-37145.9,25.7207], Tmin=(100,'K'), Tmax=(930.284,'K')), NASAPolynomial(coeffs=[27.26,0.00464042,1.30883e-06,-3.06671e-10,1.3513e-14,-44671.4,-116.837], Tmin=(930.284,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-310.202,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-OsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH)"""),
)

species(
    label = 'CC(O)=CCOO(530)',
    structure = adjacencyList("""1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {6,S} {14,S}
3  O u0 p2 c0 {1,S} {15,S}
4  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
5  C u0 p0 c0 {6,S} {10,S} {11,S} {12,S}
6  C u0 p0 c0 {2,S} {5,S} {7,D}
7  C u0 p0 c0 {4,S} {6,D} {13,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {3,S}
"""),
    E0 = (-326.091,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,466.714],'cm^-1')),
        HinderedRotor(inertia=(0.402434,'amu*angstrom^2'), symmetry=1, barrier=(9.25275,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.745635,'amu*angstrom^2'), symmetry=1, barrier=(17.1436,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.745233,'amu*angstrom^2'), symmetry=1, barrier=(17.1344,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.744882,'amu*angstrom^2'), symmetry=1, barrier=(17.1263,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.744862,'amu*angstrom^2'), symmetry=1, barrier=(17.1258,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.105,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.42923,0.07414,-7.13633e-05,3.55995e-08,-7.08802e-12,-39087.5,26.3634], Tmin=(100,'K'), Tmax=(1213.52,'K')), NASAPolynomial(coeffs=[15.5338,0.0243523,-9.82196e-06,1.79074e-09,-1.22985e-13,-42753.4,-49.4305], Tmin=(1213.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-326.091,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH)"""),
)

species(
    label = 'C=C(O)CCOO(531)',
    structure = adjacencyList("""1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {6,S} {14,S}
3  O u0 p2 c0 {1,S} {15,S}
4  C u0 p0 c0 {5,S} {6,S} {8,S} {9,S}
5  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
6  C u0 p0 c0 {2,S} {4,S} {7,D}
7  C u0 p0 c0 {6,D} {12,S} {13,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {3,S}
"""),
    E0 = (-319.403,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,3615,1277.5,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,350,440,435,1725,2950,3100,1380,975,1025,1650,298.243,880.509],'cm^-1')),
        HinderedRotor(inertia=(0.769308,'amu*angstrom^2'), symmetry=1, barrier=(17.6879,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.769292,'amu*angstrom^2'), symmetry=1, barrier=(17.6875,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.769326,'amu*angstrom^2'), symmetry=1, barrier=(17.6883,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.769293,'amu*angstrom^2'), symmetry=1, barrier=(17.6876,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.769275,'amu*angstrom^2'), symmetry=1, barrier=(17.6871,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.105,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0368769,0.0781008,-7.63264e-05,3.75394e-08,-7.22509e-12,-38265,27.829], Tmin=(100,'K'), Tmax=(1270.34,'K')), NASAPolynomial(coeffs=[18.6329,0.0195464,-7.18643e-06,1.25522e-09,-8.4477e-14,-42989.6,-66.3359], Tmin=(1270.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-319.403,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH)"""),
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
    label = '[CH2]CC(C)=O(450)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {4,D}
2  C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
3  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
4  C u0 p0 c0 {1,D} {2,S} {3,S}
5  C u1 p0 c0 {2,S} {11,S} {12,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (-46.6048,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,375,552.5,462.5,1710,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.086264,'amu*angstrom^2'), symmetry=1, barrier=(1.98338,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0858642,'amu*angstrom^2'), symmetry=1, barrier=(1.97419,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0859664,'amu*angstrom^2'), symmetry=1, barrier=(1.97654,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0978,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28203,0.0410003,-2.92978e-05,1.25297e-08,-2.50954e-12,-5546.14,17.7899], Tmin=(100,'K'), Tmax=(1074.54,'K')), NASAPolynomial(coeffs=[5.39524,0.0294113,-1.31201e-05,2.4927e-09,-1.74337e-13,-6215.19,2.54662], Tmin=(1074.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-46.6048,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + radical(CJCC=O)"""),
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
    label = 'CC(=O)CCO[O](533)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {7,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {5,S} {7,S} {8,S} {9,S}
5  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
6  C u0 p0 c0 {7,S} {12,S} {13,S} {14,S}
7  C u0 p0 c0 {2,D} {4,S} {6,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (-194.512,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,375,552.5,462.5,1710,180,2337.58],'cm^-1')),
        HinderedRotor(inertia=(0.332839,'amu*angstrom^2'), symmetry=1, barrier=(7.65262,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.33276,'amu*angstrom^2'), symmetry=1, barrier=(7.65081,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.332825,'amu*angstrom^2'), symmetry=1, barrier=(7.6523,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.332704,'amu*angstrom^2'), symmetry=1, barrier=(7.64951,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.097,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05225,0.0719006,-0.000102043,9.047e-08,-3.25306e-11,-23295,25.7871], Tmin=(100,'K'), Tmax=(818.417,'K')), NASAPolynomial(coeffs=[5.29226,0.038656,-1.81625e-05,3.44782e-09,-2.37614e-13,-23569.7,8.74309], Tmin=(818.417,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-194.512,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + radical(ROOJ)"""),
)

species(
    label = 'CH3CO(35)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {3,D}
2 C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
3 C u1 p0 c0 {1,D} {2,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
"""),
    E0 = (-21.2938,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,300.038,2569.12],'cm^-1')),
        HinderedRotor(inertia=(0.00185183,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (43.0446,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2980.68,'J/mol'), sigma=(4.94225,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=465.58 K, Pc=56.03 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.63678,0.00450844,1.80295e-05,-2.2597e-08,7.88625e-12,-2544.95,9.39483], Tmin=(100,'K'), Tmax=(1011.57,'K')), NASAPolynomial(coeffs=[4.394,0.0104824,-4.12727e-06,7.69396e-10,-5.45058e-14,-3157,3.46495], Tmin=(1011.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-21.2938,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""acetyl""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]COO(159)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {2,S} {3,S}
2 O u0 p2 c0 {1,S} {9,S}
3 C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
4 C u1 p0 c0 {3,S} {7,S} {8,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {2,S}
"""),
    E0 = (34.3631,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.258052,'amu*angstrom^2'), symmetry=1, barrier=(5.93312,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.258018,'amu*angstrom^2'), symmetry=1, barrier=(5.93233,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0141803,'amu*angstrom^2'), symmetry=1, barrier=(16.6287,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (61.0599,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.44847,0.0313149,-2.37472e-05,9.5034e-09,-1.55883e-12,4191.03,17.4641], Tmin=(100,'K'), Tmax=(1426.96,'K')), NASAPolynomial(coeffs=[8.63689,0.0139676,-5.51185e-06,9.8389e-10,-6.62176e-14,2424.92,-14.5917], Tmin=(1426.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(34.3631,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), label="""CH2CH2OOH""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'CC([O])=CCOO(534)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {14,S}
3  O u1 p2 c0 {7,S}
4  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
5  C u0 p0 c0 {7,S} {10,S} {11,S} {12,S}
6  C u0 p0 c0 {4,S} {7,D} {13,S}
7  C u0 p0 c0 {3,S} {5,S} {6,D}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {2,S}
"""),
    E0 = (-188.286,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,350,440,435,1725,341.873,342.742],'cm^-1')),
        HinderedRotor(inertia=(0.102654,'amu*angstrom^2'), symmetry=1, barrier=(8.5361,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102447,'amu*angstrom^2'), symmetry=1, barrier=(8.54121,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102759,'amu*angstrom^2'), symmetry=1, barrier=(8.54694,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.247438,'amu*angstrom^2'), symmetry=1, barrier=(20.6499,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.097,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03884,0.067214,-6.82583e-05,3.86535e-08,-9.07396e-12,-22540.7,25.6834], Tmin=(100,'K'), Tmax=(1016.08,'K')), NASAPolynomial(coeffs=[10.528,0.0298594,-1.3115e-05,2.47436e-09,-1.72608e-13,-24469.1,-20.248], Tmin=(1016.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-188.286,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ)"""),
)

species(
    label = 'CH3(21)',
    structure = adjacencyList("""multiplicity 2
1 C u1 p0 c0 {2,S} {3,S} {4,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
"""),
    E0 = (136.55,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([532.915,1389.47,1392.77,2779.27,3448.38,3448.47],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (15.0346,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1235.53,'J/mol'), sigma=(3.758e-10,'m'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with fixed Lennard Jones Parameters. This is the fallback method! Try improving transport databases!"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.94801,0.000827589,8.34936e-06,-9.82641e-09,3.80107e-12,16425.4,0.336652], Tmin=(100,'K'), Tmax=(660.437,'K')), NASAPolynomial(coeffs=[3.2217,0.00522646,-1.64125e-06,2.58224e-10,-1.62579e-14,16521.3,3.53938], Tmin=(660.437,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(136.55,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""CH3""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'O=[C]CCOO(536)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {11,S}
3  O u0 p2 c0 {6,D}
4  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
5  C u0 p0 c0 {4,S} {6,S} {9,S} {10,S}
6  C u1 p0 c0 {3,D} {5,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {2,S}
"""),
    E0 = (-121.386,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(0.3921,'amu*angstrom^2'), symmetry=1, barrier=(9.01516,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.391917,'amu*angstrom^2'), symmetry=1, barrier=(9.01094,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.39205,'amu*angstrom^2'), symmetry=1, barrier=(9.014,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.39203,'amu*angstrom^2'), symmetry=1, barrier=(9.01354,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (89.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05407,0.0752736,-0.000135797,1.26935e-07,-4.46351e-11,-14503.3,22.6872], Tmin=(100,'K'), Tmax=(874.567,'K')), NASAPolynomial(coeffs=[6.2427,0.0277304,-1.34137e-05,2.51242e-09,-1.68808e-13,-14500.2,3.55676], Tmin=(874.567,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-121.386,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCCJ=O)"""),
)

species(
    label = 'C=C([O])CCOO(537)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {14,S}
3  O u1 p2 c0 {6,S}
4  C u0 p0 c0 {5,S} {6,S} {8,S} {9,S}
5  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
6  C u0 p0 c0 {3,S} {4,S} {7,D}
7  C u0 p0 c0 {6,D} {12,S} {13,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {2,S}
"""),
    E0 = (-181.598,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,350,440,435,1725,2950,3100,1380,975,1025,1650,329.544,329.63,329.641],'cm^-1')),
        HinderedRotor(inertia=(0.00155138,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166461,'amu*angstrom^2'), symmetry=1, barrier=(12.8348,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166457,'amu*angstrom^2'), symmetry=1, barrier=(12.8349,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.320218,'amu*angstrom^2'), symmetry=1, barrier=(24.6925,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.097,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.752388,0.0700473,-6.97956e-05,3.67837e-08,-7.8179e-12,-21723.1,26.7605], Tmin=(100,'K'), Tmax=(1133.12,'K')), NASAPolynomial(coeffs=[13.4148,0.0253467,-1.06203e-05,1.96724e-09,-1.3615e-13,-24592.6,-35.9106], Tmin=(1133.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-181.598,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
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
    E0 = (-84.0617,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (197.013,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-113.274,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-110.14,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-115.407,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (35.9543,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (97.1631,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (92.9519,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (103.388,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (95.0466,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (110.077,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction5',
    reactants = ['OH(9)', 'CC(=O)CC[O](532)'],
    products = ['CC(=O)CCOO(522)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(2e+07,'m^3/(mol*s)'), n=1.78837e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_N-2R->C_N-3R!H->N_2OS->O',), comment="""Estimated from node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_N-2R->C_N-3R!H->N_2OS->O in family R_Recombination."""),
)

reaction(
    label = 'reaction1',
    reactants = ['CO(15)', 'CCCOO(528)'],
    products = ['CC(=O)CCOO(522)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(538,'cm^3/(mol*s)'), n=3.29, Ea=(437.228,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [CO;Cs_Cs] for rate rule [CO;C_methyl_C_pri]
Euclidian distance = 1.0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=C(C)OCOO(529)'],
    products = ['CC(=O)CCOO(522)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(7.50839e+11,'s^-1'), n=0.136197, Ea=(117.046,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.01703312323488117, var=1.3519057502930525, Tref=1000.0, N=4, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R in family 1,3_sigmatropic_rearrangement."""),
)

reaction(
    label = 'reaction3',
    reactants = ['CC(O)=CCOO(530)'],
    products = ['CC(=O)CCOO(522)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(5.51726e-18,'s^-1'), n=8.84109, Ea=(136.068,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R!H->C_N-3C-inRing_N-1R!H->N_Ext-3C-R_Ext-5R!H-R_Ext-1C-R_N-Sp-6R!H#5R!H_5R!H->C',), comment="""Estimated from node Root_3R!H->C_N-3C-inRing_N-1R!H->N_Ext-3C-R_Ext-5R!H-R_Ext-1C-R_N-Sp-6R!H#5R!H_5R!H->C in family Ketoenol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['C=C(O)CCOO(531)'],
    products = ['CC(=O)CCOO(522)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.34803e-38,'s^-1'), n=14.7488, Ea=(124.114,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R!H->C_N-3C-inRing_N-1R!H->N_Ext-1C-R_5R!H->C_N-5C-inRing_Ext-5C-R_N-6R!H->O_N-6BrCClFINPSSi->N_Sp-6C-5C_Ext-6C-R',), comment="""Estimated from node Root_3R!H->C_N-3C-inRing_N-1R!H->N_Ext-1C-R_5R!H->C_N-5C-inRing_Ext-5C-R_N-6R!H->O_N-6BrCClFINPSSi->N_Sp-6C-5C_Ext-6C-R in family Ketoenol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['HO2(12)', '[CH2]CC(C)=O(450)'],
    products = ['CC(=O)CCOO(522)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(308.407,'m^3/(mol*s)'), n=0.967216, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.0049236047651, var=0.0302625193734, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R',), comment="""Estimated from node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R in family R_Recombination."""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(5)', 'CC(=O)CCO[O](533)'],
    products = ['CC(=O)CCOO(522)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(5250.69,'m^3/(mol*s)'), n=1.27262, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O in family R_Recombination."""),
)

reaction(
    label = 'reaction8',
    reactants = ['CH3CO(35)', '[CH2]COO(159)'],
    products = ['CC(=O)CCOO(522)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(6.30526e+07,'m^3/(mol*s)'), n=-0.283333, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-6.50053689359e-11, var=0.305422193575, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_3R!H->O_Ext-1C-R_5R!H->C',), comment="""Estimated from node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_3R!H->O_Ext-1C-R_5R!H->C in family R_Recombination."""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(5)', 'CC([O])=CCOO(534)'],
    products = ['CC(=O)CCOO(522)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R in family R_Recombination."""),
)

reaction(
    label = 'reaction10',
    reactants = ['CH3(21)', 'O=[C]CCOO(536)'],
    products = ['CC(=O)CCOO(522)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4.52905e+07,'m^3/(mol*s)'), n=-0.17, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-8.93790333559e-10, var=0.843006791659, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_3R!H->O',), comment="""Estimated from node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_3R!H->O in family R_Recombination."""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(5)', 'C=C([O])CCOO(537)'],
    products = ['CC(=O)CCOO(522)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(9.17499e+07,'m^3/(mol*s)'), n=0.115342, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_3R!H->C_Sp-3C-2CN_Ext-3C-R_N-Sp-4R!H=3C_Sp-4R!H-3C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_3R!H->C_Sp-3C-2CN_Ext-3C-R_N-Sp-4R!H=3C_Sp-4R!H-3C in family R_Recombination."""),
)

network(
    label = 'PDepNetwork #559',
    isomers = [
        'CC(=O)CCOO(522)',
    ],
    reactants = [
        ('OH(9)', 'CC(=O)CC[O](532)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #559',
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

