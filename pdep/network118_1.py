species(
    label = 'C#COO[OH+][O-](172)',
    structure = adjacencyList("""1 O u0 p1 c+1 {2,S} {4,S} {7,S}
2 O u0 p2 c0 {1,S} {3,S}
3 O u0 p2 c0 {2,S} {5,S}
4 O u0 p3 c-1 {1,S}
5 C u0 p0 c0 {3,S} {6,T}
6 C u0 p0 c0 {5,T} {8,S}
7 H u0 p0 c0 {1,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (75.212,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,750,770,3400,2100,253.126,255.475,264.678,325.932,4000,4000,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00470194,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00470194,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00470194,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (90.0348,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.3404,0.0441091,-8.64887e-05,8.75665e-08,-3.32736e-11,9098.29,6.66639], Tmin=(100,'K'), Tmax=(830.154,'K')), NASAPolynomial(coeffs=[4.01192,0.0191193,-1.07336e-05,2.15561e-09,-1.51541e-13,9404.34,2.42828], Tmin=(830.154,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(75.212,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-OsCt) + missing(O0sc-O4sc) + group(Ct-CtOs) + group(Ct-CtH)"""),
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
    label = 'C#CO[O](106)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {2,S} {3,S}
2 O u1 p2 c0 {1,S}
3 C u0 p0 c0 {1,S} {4,T}
4 C u0 p0 c0 {3,T} {5,S}
5 H u0 p0 c0 {4,S}
"""),
    E0 = (341.373,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,750,770,3400,2100,207.446,207.514],'cm^-1')),
        HinderedRotor(inertia=(1.42756,'amu*angstrom^2'), symmetry=1, barrier=(43.6062,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (57.0281,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.84092,0.0270254,-4.16766e-05,3.32106e-08,-1.02953e-11,41097.9,11.7897], Tmin=(100,'K'), Tmax=(872.95,'K')), NASAPolynomial(coeffs=[6.69656,0.00695041,-3.04407e-06,5.47518e-10,-3.61754e-14,40516.5,-5.76216], Tmin=(872.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(341.373,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCt) + group(O2s-OsH) + group(Ct-CtOs) + group(Ct-CtH) + radical(ROOJ)"""),
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
    label = 'C#COO[O+][O-](989)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {2,S} {5,S}
2 O u0 p2 c0 {1,S} {3,S}
3 O u1 p1 c+1 {2,S} {4,S}
4 O u0 p3 c-1 {3,S}
5 C u0 p0 c0 {1,S} {6,T}
6 C u0 p0 c0 {5,T} {7,S}
7 H u0 p0 c0 {6,S}
"""),
    E0 = (293.136,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,750,770,3400,2100,458.689,1375.09,4000,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.771145,'amu*angstrom^2'), symmetry=1, barrier=(17.7301,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.771871,'amu*angstrom^2'), symmetry=1, barrier=(17.7468,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.765449,'amu*angstrom^2'), symmetry=1, barrier=(17.5992,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (89.0268,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.43797,0.0419494,-8.58606e-05,8.77812e-08,-3.3467e-11,35305,5.99005], Tmin=(100,'K'), Tmax=(829.076,'K')), NASAPolynomial(coeffs=[4.06833,0.0170598,-1.00295e-05,2.03836e-09,-1.43938e-13,35619.7,1.95858], Tmin=(829.076,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(293.136,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(145.503,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-OsCt) + missing(O0sc-O4sc) + group(Ct-CtOs) + group(Ct-CtH) + radical(CsOJ)"""),
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
    label = 'HCCO(34)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {3,D}
2 C u1 p0 c0 {3,D} {4,S}
3 C u0 p0 c0 {1,D} {2,D}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (165.191,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,231.114,1089.61,3388.99],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (41.0287,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2433.35,'J/mol'), sigma=(4.2737,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=380.08 K, Pc=70.74 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.47174,0.0113766,-1.1308e-05,6.13951e-09,-1.35088e-12,19887.1,6.87336], Tmin=(100,'K'), Tmax=(1096.2,'K')), NASAPolynomial(coeffs=[5.39217,0.00436893,-1.71893e-06,3.07751e-10,-2.08706e-14,19466,-2.56797], Tmin=(1096.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(165.191,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""HCCO""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[O]O[OH+][O-](991)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p1 c+1 {2,S} {3,S} {5,S}
2 O u0 p2 c0 {1,S} {4,S}
3 O u0 p3 c-1 {1,S}
4 O u1 p2 c0 {2,S}
5 H u0 p0 c0 {1,S}
"""),
    E0 = (-109.594,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([4000,4000,4000,4000,4000,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.51379,'amu*angstrom^2'), symmetry=1, barrier=(11.813,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (65.0056,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.7725,0.0130465,-4.38867e-05,5.59007e-08,-2.29167e-11,-13180.5,-0.0882193], Tmin=(100,'K'), Tmax=(888.654,'K')), NASAPolynomial(coeffs=[-0.567748,0.0116732,-6.27454e-06,1.20637e-09,-8.10388e-14,-11583.5,24.9839], Tmin=(888.654,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-109.594,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-OsH) + missing(O0sc-O4sc) + radical(ROOJ)"""),
)

species(
    label = 'C2H(22)',
    structure = adjacencyList("""multiplicity 2
1 C u0 p0 c0 {2,T} {3,S}
2 C u1 p0 c0 {1,T}
3 H u0 p0 c0 {1,S}
"""),
    E0 = (556.809,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (25.0293,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(2353.39,'J/mol'), sigma=(4.42547,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=367.59 K, Pc=61.61 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.01451,0.0139906,-3.08142e-05,3.10833e-08,-1.10944e-11,66983,5.75944], Tmin=(100,'K'), Tmax=(918.727,'K')), NASAPolynomial(coeffs=[3.14384,0.0049849,-2.32626e-06,4.0876e-10,-2.5579e-14,67315.6,7.08562], Tmin=(918.727,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(556.809,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(62.3585,'J/(mol*K)'), label="""HC2""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[C]#COO[OH+][O-](992)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p1 c+1 {2,S} {4,S} {7,S}
2 O u0 p2 c0 {1,S} {3,S}
3 O u0 p2 c0 {2,S} {5,S}
4 O u0 p3 c-1 {1,S}
5 C u0 p0 c0 {3,S} {6,T}
6 C u1 p0 c0 {5,T}
7 H u0 p0 c0 {1,S}
"""),
    E0 = (412.356,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,183.688,184.093,184.498,972.84,3575.21,4000,4000,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.314637,'amu*angstrom^2'), symmetry=1, barrier=(7.33713,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.314637,'amu*angstrom^2'), symmetry=1, barrier=(7.33713,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.314637,'amu*angstrom^2'), symmetry=1, barrier=(7.33713,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (89.0268,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.35659,0.0476948,-0.000108977,1.13839e-07,-4.24947e-11,49643.1,7.35235], Tmin=(100,'K'), Tmax=(877.246,'K')), NASAPolynomial(coeffs=[2.61064,0.0182289,-1.01898e-05,1.98179e-09,-1.34408e-13,50687.7,12.3682], Tmin=(877.246,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(412.356,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(145.503,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-OsCt) + missing(O0sc-O4sc) + group(Ct-CtOs) + group(Ct-CtH) + radical(Acetyl)"""),
)

species(
    label = '[C]=COO[OH+][O-](993)',
    structure = adjacencyList("""1 O u0 p1 c+1 {2,S} {4,S} {7,S}
2 O u0 p2 c0 {1,S} {3,S}
3 O u0 p2 c0 {2,S} {5,S}
4 O u0 p3 c-1 {1,S}
5 C u0 p0 c0 {3,S} {6,D} {8,S}
6 C u0 p1 c0 {5,D}
7 H u0 p0 c0 {1,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (229.652,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,186.615,194.266,449.7,1784.66,4000,4000,4000,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.221888,'amu*angstrom^2'), symmetry=1, barrier=(5.35704,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.221888,'amu*angstrom^2'), symmetry=1, barrier=(5.35704,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.221888,'amu*angstrom^2'), symmetry=1, barrier=(5.35704,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (90.0348,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.25901,0.0548577,-0.000135266,1.49824e-07,-5.82022e-11,27667.6,7.89908], Tmin=(100,'K'), Tmax=(869.449,'K')), NASAPolynomial(coeffs=[-0.6964,0.0269396,-1.54778e-05,3.05637e-09,-2.09749e-13,29750.6,30.7675], Tmin=(869.449,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(229.652,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + missing(O0sc-O4sc) + group(Cds-CdOsH) + group(CdJ2_singlet-Cds)"""),
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
    E0 = (-82.9694,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (155.275,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (316.154,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (4.30694,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (258.441,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (435.374,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (90.7321,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C#COO[OH+][O-](172)'],
    products = ['O3(14)', 'CH2CO(27)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(2.46545e+14,'s^-1'), n=0.20628, Ea=(30.5927,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(3000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-4R!H->C',), comment="""Estimated from node Root_N-4R!H->C in family Retroene."""),
)

reaction(
    label = 'reaction2',
    reactants = ['HO2(12)', 'C#CO[O](106)'],
    products = ['C#COO[OH+][O-](172)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(106477,'m^3/(mol*s)'), n=0.348287, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.0108230153501, var=2.70964383578, Tref=1000.0, N=19, data_mean=0.0, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R',), comment="""Estimated from node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R in family R_Recombination."""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(5)', 'C#COO[O+][O-](989)'],
    products = ['C#COO[OH+][O-](172)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(5250.69,'m^3/(mol*s)'), n=1.27262, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O in family R_Recombination."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O][OH+][O-](990)', 'HCCO(34)'],
    products = ['C#COO[OH+][O-](172)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.81e+06,'m^3/(mol*s)'), n=-5.80997e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R_N-2R->C',), comment="""Estimated from node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R_N-2R->C in family R_Recombination."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]O[OH+][O-](991)', 'C2H(22)'],
    products = ['C#COO[OH+][O-](172)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(106477,'m^3/(mol*s)'), n=0.348287, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.0108230153501, var=2.70964383578, Tref=1000.0, N=19, data_mean=0.0, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R',), comment="""Estimated from node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R in family R_Recombination."""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(5)', '[C]#COO[OH+][O-](992)'],
    products = ['C#COO[OH+][O-](172)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.81e+08,'m^3/(mol*s)'), n=-8.9089e-10, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_3R!H->C_N-Sp-3C-2CN',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_3R!H->C_N-Sp-3C-2CN in family R_Recombination."""),
)

reaction(
    label = 'reaction7',
    reactants = ['[C]=COO[OH+][O-](993)'],
    products = ['C#COO[OH+][O-](172)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(3.95159e+18,'s^-1'), n=-1.97082, Ea=(49.8542,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.8305765209480382, var=17.041866923497537, Tref=1000.0, N=4, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root in family Singlet_Carbene_Intra_Disproportionation."""),
)

network(
    label = 'PDepNetwork #118',
    isomers = [
        'C#COO[OH+][O-](172)',
    ],
    reactants = [
        ('O3(14)', 'CH2CO(27)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #118',
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

