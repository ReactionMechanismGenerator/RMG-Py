species(
    label = 'ethane(1)',
    structure = SMILES('CC'),
    E0 = (-105.416,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,180,899.833,899.917,1718.13,2306.62,2310.39,2310.6,2310.8],'cm^-1')),
        HinderedRotor(inertia=(0.00591986,'amu*angstrom^2'), symmetry=1, barrier=(3.39703,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (30.069,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2097.75,'J/mol'), sigma=(4.302,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.05996,-0.00393683,5.60245e-05,-6.65913e-08,2.57665e-11,-12679.2,3.54102], Tmin=(10,'K'), Tmax=(730.382,'K')), NASAPolynomial(coeffs=[-1.00362,0.0270461,-1.4284e-05,3.67956e-09,-3.72725e-13,-12026.3,25.7852], Tmin=(730.382,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-105.407,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), label="""C2H6_0""", comment="""Thermo library: t3_thermo"""),
)

species(
    label = 'CH3_0(5)',
    structure = SMILES('[CH3]'),
    E0 = (132.783,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([480.666,1347.66,1347.71,2724.95,3382.65,3382.8],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (15.0345,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.98397,0.000880651,9.08622e-06,-1.19151e-08,5.2772e-12,15970.3,0.235907], Tmin=(10,'K'), Tmax=(562.242,'K')), NASAPolynomial(coeffs=[3.46132,0.00459902,-8.34018e-07,-1.52309e-10,4.68974e-14,16029.1,2.45645], Tmin=(562.242,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(132.778,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""CH3_0""", comment="""Thermo library: t3_thermo"""),
)

species(
    label = 'CH2(S)(3)',
    structure = SMILES('[CH2]'),
    E0 = (419.091,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1369.93,2896.01,2896.02],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.10264,-0.00144068,5.45069e-06,-3.58002e-09,7.56192e-13,50400.6,-0.411765], Tmin=(100,'K'), Tmax=(1442.36,'K')), NASAPolynomial(coeffs=[2.62648,0.00394763,-1.49924e-06,2.54539e-10,-1.62956e-14,50691.8,6.78378], Tmin=(1442.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(419.091,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(S)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'CH4(4)',
    structure = SMILES('C'),
    E0 = (-84.435,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1092.7,1092.7,1649.06,2067.52,2067.53,2067.53,2067.54,3478.74,4000],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (16.0425,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.20542,-0.00535557,2.51123e-05,-2.13763e-08,5.97524e-12,-10161.9,-0.921279], Tmin=(100,'K'), Tmax=(1084.12,'K')), NASAPolynomial(coeffs=[0.908266,0.0114541,-4.57174e-06,8.29192e-10,-5.66315e-14,-9719.97,13.9931], Tmin=(1084.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-84.435,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""CH4""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'H(7)',
    structure = SMILES('[H]'),
    E0 = (211.805,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (1.00794,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,9.24385e-15,-1.3678e-17,6.66185e-21,-1.00107e-24,25474.2,-0.444973], Tmin=(100,'K'), Tmax=(3459.6,'K')), NASAPolynomial(coeffs=[2.5,9.20456e-12,-3.58608e-15,6.15199e-19,-3.92042e-23,25474.2,-0.444973], Tmin=(3459.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(211.805,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""H""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'C[CH2](6)',
    structure = SMILES('C[CH2]'),
    E0 = (99.8726,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,1681.65,1683.61,1683.66,1683.82,1684.41],'cm^-1')),
        HinderedRotor(inertia=(0.184397,'amu*angstrom^2'), symmetry=1, barrier=(9.11786,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (29.0611,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.81667,-0.000564623,3.72036e-05,-4.21198e-08,1.48964e-11,12022.7,6.97543], Tmin=(100,'K'), Tmax=(939.904,'K')), NASAPolynomial(coeffs=[3.20633,0.0137961,-4.48777e-06,7.6705e-10,-5.26061e-14,11617.8,7.11806], Tmin=(939.904,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(99.8726,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo library: t3_thermo + radical(CCJ)"""),
)

species(
    label = 'N2',
    structure = SMILES('N#N'),
    E0 = (-8.64289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0135,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(2684.29,'J/mol'), sigma=(3.461,'angstroms'), dipoleMoment=(1.781,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""OneDMinN2"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53101,-0.000123661,-5.02999e-07,2.43531e-09,-1.40881e-12,-1046.98,2.96747], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.95258,0.0013969,-4.92632e-07,7.8601e-11,-4.60755e-15,-923.949,5.87189], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-8.64289,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""N2""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'Ar',
    structure = SMILES('[Ar]'),
    E0 = (-6.19738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (39.348,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1134.93,'J/mol'), sigma=(3.33,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,4.37967], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,4.37967], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-6.19738,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""Ar""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'He',
    structure = SMILES('[He]'),
    E0 = (-6.19738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (4.0026,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(84.8076,'J/mol'), sigma=(2.576,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,0.928724], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,0.928724], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-6.19738,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""He""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'Ne',
    structure = SMILES('[Ne]'),
    E0 = (-6.19738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (20.1797,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1235.53,'J/mol'), sigma=(3.758e-10,'m'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with fixed Lennard Jones Parameters. This is the fallback method! Try improving transport databases!"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,3.35532], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,3.35532], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-6.19738,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""Ne""", comment="""Thermo library: primaryThermoLibrary"""),
)

transitionState(
    label = 'TS1',
    E0 = (266.131,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (334.656,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (311.677,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['CH3_0(5)', 'CH3_0(5)'],
    products = ['ethane(1)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(9.45e+14,'cm^3/(mol*s)'), n=-0.538, Ea=(135.1,'cal/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2000,'K'), comment="""Matched reaction 9 CH3 + CH3 <=> C2H6 in R_Recombination/training
This reaction matched rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing]
family: R_Recombination"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(S)(3)', 'CH4(4)'],
    products = ['ethane(1)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.74695e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;Cs_H] for rate rule [carbene;C_methane]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(7)', 'C[CH2](6)'],
    products = ['ethane(1)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1e+14,'cm^3/(mol*s)','+|-',1e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""Matched reaction 58 H + C2H5 <=> C2H6-2 in R_Recombination/training
This reaction matched rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_3R!H->C_Sp-3C-2CN]
family: R_Recombination"""),
)

network(
    label = 'PDepNetwork #3',
    isomers = [
        'ethane(1)',
    ],
    reactants = [
        ('CH3_0(5)', 'CH3_0(5)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ar': 0.25,
        'He': 0.25,
        'Ne': 0.25,
    },
)

pressureDependence(
    label = 'PDepNetwork #3',
    Tmin = (300,'K'),
    Tmax = (2200,'K'),
    Tcount = 10,
    Tlist = ([301.603,314.817,343.437,392.555,471.896,599.244,806.147,1141.38,1635.51,2117.45],'K'),
    Pmin = (0.01,'bar'),
    Pmax = (100,'bar'),
    Pcount = 10,
    Plist = ([0.0105834,0.0165191,0.0385289,0.1236,0.486554,2.05527,8.09061,25.9546,60.5359,94.488],'bar'),
    maximumGrainSize = (0.5,'kcal/mol'),
    minimumGrainCount = 250,
    method = 'modified strong collision',
    interpolationModel = ('Chebyshev', 6, 4),
    activeKRotor = True,
    activeJRotor = True,
    rmgmode = True,
)

