#!/usr/bin/env python
# encoding: utf-8

name = "/home/alongd/Code/RMG-Py/importer/NOx2018"
shortDesc = u"/home/alongd/Code/RMG-Py/importer/NOx2018/kinetics.txt"
longDesc = u"""
Unknown source
"""
entry(
    index = 1,
    label = "H + O2 <=> O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(15286, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H + O2 <=> O + OH""",
)

entry(
    index = 2,
    label = "O + H2 <=> OH + H",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(3.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(7948, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(8.8e+14, 'cm^3/(mol*s)'), n=0, Ea=(19175, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is O + H2 <=> OH + H""",
)

entry(
    index = 3,
    label = "OH + H2 <=> H + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.2e+08, 'cm^3/(mol*s)'),
        n = 1.51,
        Ea = (3430, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is OH + H2 <=> H + H2O""",
)

entry(
    index = 4,
    label = "OH + OH <=> O + H2O",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(2e+07, 'cm^3/(mol*s)'), n=1.651, Ea=(631, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (2.6e+11, 'cm^3/(mol*s)'),
                n = -0.057,
                Ea = (-827, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is OH + OH <=> O + H2O""",
)

entry(
    index = 5,
    label = "H2 <=> H + H",
    degeneracy = 1,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(
            A = (4.6e+19, 'cm^3/(mol*s)'),
            n = -1.4,
            Ea = (104380, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        efficiencies = {'[C-]#[O+]': 1.9, '[H][H]': 2.5, 'O=C=O': 3.8, 'O': 12, '[Ar]': 0},
    ),
    shortDesc = u"""The chemkin file reaction is H2 <=> H + H""",
)

entry(
    index = 6,
    label = "H2 + AR <=> H + H + AR",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.8e+18, 'cm^3/(mol*s)'),
        n = -1.1,
        Ea = (104380, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is H2 + AR <=> H + H + AR""",
)

entry(
    index = 7,
    label = "H + O <=> OH",
    degeneracy = 1,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(A=(4.7e+18, 'cm^6/(mol^2*s)'), n=-1, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        efficiencies = {'[C-]#[O+]': 1.9, '[H][H]': 2.5, 'O=C=O': 3.8, 'O': 12, '[Ar]': 0.75},
    ),
    shortDesc = u"""The chemkin file reaction is H + O <=> OH""",
)

entry(
    index = 8,
    label = "O + O <=> O2",
    degeneracy = 1,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(
            A = (1.9e+13, 'cm^6/(mol^2*s)'),
            n = 0,
            Ea = (-1788, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        efficiencies = {'[C-]#[O+]': 1.9, '[H][H]': 2.5, 'O=C=O': 3.8, 'O': 12, '[Ar]': 0},
    ),
    shortDesc = u"""The chemkin file reaction is O + O <=> O2""",
)

entry(
    index = 9,
    label = "H2O <=> H + OH",
    degeneracy = 1,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(
            A = (6.1e+27, 'cm^3/(mol*s)'),
            n = -3.322,
            Ea = (120790, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        efficiencies = {'O=C=O': 3.8, 'O': 0, '[H][H]': 3, '[O][O]': 1.5, 'N#N': 2, '[C-]#[O+]': 1.9},
    ),
    shortDesc = u"""The chemkin file reaction is H2O <=> H + OH""",
)

entry(
    index = 10,
    label = "H2O + H2O <=> H + OH + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1e+26, 'cm^3/(mol*s)'),
        n = -2.44,
        Ea = (120180, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is H2O + H2O <=> H + OH + H2O""",
)

entry(
    index = 11,
    label = "H + O2 <=> HO2",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(4.7e+12, 'cm^3/(mol*s)'), n=0.44, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (6.366e+20, 'cm^6/(mol^2*s)'),
            n = -1.72,
            Ea = (524.8, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.5,
        T3 = (1e-30, 'K'),
        T1 = (1e+30, 'K'),
        efficiencies = {'O=C=O': 3.8, 'O': 14, '[H][H]': 2, '[O][O]': 0.78, '[C-]#[O+]': 1.9, '[Ar]': 0.67},
    ),
    shortDesc = u"""The chemkin file reaction is H + O2 <=> HO2""",
)

entry(
    index = 12,
    label = "HO2 + H <=> H2 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.8e+06, 'cm^3/(mol*s)'),
        n = 2.09,
        Ea = (-1451, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is HO2 + H <=> H2 + O2""",
)

entry(
    index = 13,
    label = "HO2 + H <=> OH + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.1e+13, 'cm^3/(mol*s)'), n=0, Ea=(295, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HO2 + H <=> OH + OH""",
)

entry(
    index = 14,
    label = "HO2 + H <=> H2O + O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.4e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HO2 + H <=> H2O + O""",
)

entry(
    index = 15,
    label = "HO2 + O <=> O2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.9e+10, 'cm^3/(mol*s)'), n=1, Ea=(-724, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HO2 + O <=> O2 + OH""",
)

entry(
    index = 16,
    label = "HO2 + OH <=> H2O + O2",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(
                A = (1.9e+20, 'cm^3/(mol*s)'),
                n = -2.49,
                Ea = (584, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.2e+09, 'cm^3/(mol*s)'),
                n = 1.24,
                Ea = (-1310, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is HO2 + OH <=> H2O + O2""",
)

entry(
    index = 17,
    label = "HO2 + HO2 <=> H2O2 + O2",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(
                A = (1.2e+09, 'cm^3/(mol*s)'),
                n = 0.7712,
                Ea = (-1825, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.3e+12, 'cm^3/(mol*s)'),
                n = 0.295,
                Ea = (7397, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is HO2 + HO2 <=> H2O2 + O2""",
)

entry(
    index = 18,
    label = "H2O2 <=> OH + OH",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2e+12, 's^-1'), n=0.9, Ea=(48749, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (2.5e+24, 'cm^3/(mol*s)'),
            n = -2.3,
            Ea = (48749, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.43,
        T3 = (1e-30, 'K'),
        T1 = (1e+30, 'K'),
        efficiencies = {'OO': 7.7, 'O=C=O': 1.6, 'O': 7.5, '[H][H]': 3.7, '[O][O]': 1.2, 'N#N': 1.5, '[C-]#[O+]': 2.8, '[Ar]': 1},
    ),
    shortDesc = u"""The chemkin file reaction is H2O2 <=> OH + OH""",
)

entry(
    index = 19,
    label = "H2O2 + H <=> H2O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(3970, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2O2 + H <=> H2O + OH""",
)

entry(
    index = 20,
    label = "H2O2 + H <=> HO2 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(7950, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2O2 + H <=> HO2 + H2""",
)

entry(
    index = 21,
    label = "H2O2 + O <=> HO2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.6e+06, 'cm^3/(mol*s)'), n=2, Ea=(3970, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2O2 + O <=> HO2 + OH""",
)

entry(
    index = 22,
    label = "H2O2 + OH <=> HO2 + H2O",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(1.7e+12, 'cm^3/(mol*s)'), n=0, Ea=(318, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(7.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(7270, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is H2O2 + OH <=> HO2 + H2O""",
)

entry(
    index = 23,
    label = "CO + O <=> CO2",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.8e+10, 'cm^3/(mol*s)'), n=0, Ea=(2384, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.4e+24, 'cm^6/(mol^2*s)'),
            n = -2.79,
            Ea = (4191, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 1,
        T3 = (1e-30, 'K'),
        T1 = (1e+30, 'K'),
        T2 = (1e+30, 'K'),
        efficiencies = {'[C-]#[O+]': 1.9, '[H][H]': 2.5, 'O=C=O': 3.8, 'O': 12},
    ),
    shortDesc = u"""The chemkin file reaction is CO + O <=> CO2""",
)

entry(
    index = 24,
    label = "CO + OH <=> CO2 + H",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01315, 0.1315, 1.315, 13.158, 131.58], 'atm'),
        arrhenius = [
            Arrhenius(A=(210000, 'cm^3/(mol*s)'), n=1.9, Ea=(-1064, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (250000, 'cm^3/(mol*s)'),
                n = 1.88,
                Ea = (-1043, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(A=(870000, 'cm^3/(mol*s)'), n=1.73, Ea=(-685, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(6.8e+06, 'cm^3/(mol*s)'), n=1.48, Ea=(48, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2.3e+07, 'cm^3/(mol*s)'), n=1.35, Ea=(974, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CO + OH <=> CO2 + H""",
)

entry(
    index = 25,
    label = "CO + OH <=> HOCO",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.013158, 0.13158, 1.3158, 13.158, 131.58], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (1.7e+15, 'cm^3/(mol*s)'),
                n = -2.68,
                Ea = (859, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (5.9e+18, 'cm^3/(mol*s)'),
                n = -3.35,
                Ea = (887, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.6e+20, 'cm^3/(mol*s)'),
                n = -3.5,
                Ea = (1309, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (7.1e+20, 'cm^3/(mol*s)'),
                n = -3.32,
                Ea = (1763, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.1e+20, 'cm^3/(mol*s)'),
                n = -2.78,
                Ea = (2056, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CO + OH <=> HOCO""",
)

entry(
    index = 26,
    label = "CO + HO2 <=> CO2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (160000, 'cm^3/(mol*s)'),
        n = 2.18,
        Ea = (17943, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CO + HO2 <=> CO2 + OH""",
)

entry(
    index = 27,
    label = "CO + O2 <=> CO2 + O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.7e+12, 'cm^3/(mol*s)'), n=0, Ea=(60500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CO + O2 <=> CO2 + O""",
)

entry(
    index = 28,
    label = "CO + H2O2 <=> HOCO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(36000, 'cm^3/(mol*s)'), n=2.5, Ea=(28660, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CO + H2O2 <=> HOCO + OH""",
)

entry(
    index = 29,
    label = "HOCO <=> CO2 + H",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(8.2e+11, 's^-1'), n=0.413, Ea=(35335, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (6e+26, 'cm^3/(mol*s)'),
            n = -3.148,
            Ea = (37116, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.39,
        T3 = (1e-30, 'K'),
        T1 = (1e+30, 'K'),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is HOCO <=> CO2 + H""",
)

entry(
    index = 30,
    label = "HOCO + H <=> CO2 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.1e+17, 'cm^3/(mol*s)'),
        n = -1.3475,
        Ea = (555, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is HOCO + H <=> CO2 + H2""",
)

entry(
    index = 31,
    label = "HOCO + H <=> CO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (6e+15, 'cm^3/(mol*s)'),
        n = -0.525,
        Ea = (2125, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is HOCO + H <=> CO + H2O""",
)

entry(
    index = 32,
    label = "HOCO + O <=> CO2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HOCO + O <=> CO2 + OH""",
)

entry(
    index = 33,
    label = "HOCO + OH <=> CO2 + H2O",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(4.6e+12, 'cm^3/(mol*s)'), n=0, Ea=(-89, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(9.5e+06, 'cm^3/(mol*s)'), n=2, Ea=(-89, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is HOCO + OH <=> CO2 + H2O""",
)

entry(
    index = 34,
    label = "HOCO + HO2 <=> CO2 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HOCO + HO2 <=> CO2 + H2O2""",
)

entry(
    index = 35,
    label = "HOCO + O2 <=> CO2 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+09, 'cm^3/(mol*s)'), n=1, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HOCO + O2 <=> CO2 + HO2""",
)

entry(
    index = 36,
    label = "CH2O <=> HCO + H",
    degeneracy = 1,
    kinetics = Lindemann(
        arrheniusHigh = Arrhenius(A=(8e+15, 's^-1'), n=0, Ea=(87726, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.3e+36, 'cm^3/(mol*s)'),
            n = -5.5,
            Ea = (93932, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is CH2O <=> HCO + H""",
)

entry(
    index = 37,
    label = "CH2O <=> CO + H2",
    degeneracy = 1,
    kinetics = Lindemann(
        arrheniusHigh = Arrhenius(A=(3.7e+13, 's^-1'), n=0, Ea=(71969, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4.4e+38, 'cm^3/(mol*s)'),
            n = -6.1,
            Ea = (93932, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is CH2O <=> CO + H2""",
)

entry(
    index = 38,
    label = "CH2O + H <=> HCO + H2",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiPDepArrhenius(
        arrhenius = [
            PDepArrhenius(
                pressures = ([0.04, 1, 10], 'atm'),
                arrhenius = [
                    Arrhenius(
                        A = (7.4e+23, 'cm^3/(mol*s)'),
                        n = -2.732,
                        Ea = (16379, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.4e+23, 'cm^3/(mol*s)'),
                        n = -2.355,
                        Ea = (17519, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (7.3e+23, 'cm^3/(mol*s)'),
                        n = -2.665,
                        Ea = (17634, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            PDepArrhenius(
                pressures = ([0.04, 1, 10], 'atm'),
                arrhenius = [
                    Arrhenius(
                        A = (2.1e+10, 'cm^3/(mol*s)'),
                        n = 1.057,
                        Ea = (3720, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.6e+15, 'cm^3/(mol*s)'),
                        n = -0.444,
                        Ea = (5682, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (4.2e+09, 'cm^3/(mol*s)'),
                        n = 1.294,
                        Ea = (3591, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2O + H <=> HCO + H2""",
)

entry(
    index = 39,
    label = "CH2O + H <=> H + CO + H2",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.04, 1, 10], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (7.2e+08, 'cm^3/(mol*s)'),
                n = 1.903,
                Ea = (11733, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (5.1e+07, 'cm^3/(mol*s)'),
                n = 2.182,
                Ea = (11524, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.1e+09, 'cm^3/(mol*s)'),
                n = 1.812,
                Ea = (13163, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2O + H <=> H + CO + H2""",
)

entry(
    index = 40,
    label = "CH2O + O <=> HCO + OH",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(
                A = (5.6e+31, 'cm^3/(mol*s)'),
                n = -5.189,
                Ea = (19968, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.4e+15, 'cm^3/(mol*s)'),
                n = -0.53,
                Ea = (4011, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2O + O <=> HCO + OH""",
)

entry(
    index = 41,
    label = "CH2O + O <=> H + CO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.5e+21, 'cm^3/(mol*s)'),
        n = -1.903,
        Ea = (22674, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2O + O <=> H + CO + OH""",
)

entry(
    index = 42,
    label = "CH2O + OH <=> HCO + H2O",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.04, 1, 10], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (3.6e+09, 'cm^3/(mol*s)'),
                n = 1.167,
                Ea = (-206, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.9e+09, 'cm^3/(mol*s)'),
                n = 1.256,
                Ea = (-302, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.1e+09, 'cm^3/(mol*s)'),
                n = 1.33,
                Ea = (-392, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2O + OH <=> HCO + H2O""",
)

entry(
    index = 43,
    label = "CH2O + OH <=> H + CO + H2O",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.04, 1, 10], 'atm'),
        arrhenius = [
            Arrhenius(A=(7e+10, 'cm^3/(mol*s)'), n=0.911, Ea=(8646, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (7.2e+10, 'cm^3/(mol*s)'),
                n = 0.892,
                Ea = (9310, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (8.4e+10, 'cm^3/(mol*s)'),
                n = 0.879,
                Ea = (9843, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2O + OH <=> H + CO + H2O""",
)

entry(
    index = 44,
    label = "CH2O + HO2 <=> HCO + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.4e+08, 'cm^3/(mol*s)'),
        n = 1.298,
        Ea = (12129, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2O + HO2 <=> HCO + H2O2""",
)

entry(
    index = 45,
    label = "CH2O + HO2 <=> H + CO + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.5e+14, 'cm^3/(mol*s)'),
        n = 0.027,
        Ea = (30120, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2O + HO2 <=> H + CO + H2O2""",
)

entry(
    index = 46,
    label = "CH2O + O2 <=> HCO + HO2",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(
                A = (1.8e+16, 'cm^3/(mol*s)'),
                n = -0.639,
                Ea = (45400, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (6.6e+08, 'cm^3/(mol*s)'),
                n = 1.36,
                Ea = (37324, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2O + O2 <=> HCO + HO2""",
)

entry(
    index = 47,
    label = "CH2O + O2 <=> H + CO + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.4e+15, 'cm^3/(mol*s)'),
        n = 0.027,
        Ea = (56388, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2O + O2 <=> H + CO + HO2""",
)

entry(
    index = 48,
    label = "HCO <=> H + CO",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(4.9e+16, 's^-1'), n=-0.93, Ea=(19724, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (7.4e+21, 'cm^3/(mol*s)'),
            n = -2.36,
            Ea = (19383, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.103,
        T3 = (139, 'K'),
        T1 = (10900, 'K'),
        T2 = (4550, 'K'),
        efficiencies = {'C': 5, 'O=C=O': 3, 'O': 15, '[H][H]': 2, '[He]': 1.3, '[O][O]': 1.5, 'N#N': 1.5, '[C-]#[O+]': 1.5},
    ),
    shortDesc = u"""The chemkin file reaction is HCO <=> H + CO""",
)

entry(
    index = 49,
    label = "HCO + H <=> CO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCO + H <=> CO + H2""",
)

entry(
    index = 50,
    label = "HCO + O <=> CO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCO + O <=> CO + OH""",
)

entry(
    index = 51,
    label = "HCO + O <=> CO2 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCO + O <=> CO2 + H""",
)

entry(
    index = 52,
    label = "HCO + OH <=> CO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCO + OH <=> CO + H2O""",
)

entry(
    index = 53,
    label = "HCO + O2 <=> CO + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (6.9e+06, 'cm^3/(mol*s)'),
        n = 1.9,
        Ea = (-1369, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is HCO + O2 <=> CO + HO2""",
)

entry(
    index = 54,
    label = "HCO + HO2 <=> CO2 + OH + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCO + HO2 <=> CO2 + OH + H""",
)

entry(
    index = 55,
    label = "HCO + HCO <=> CO + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.7e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCO + HCO <=> CO + CH2O""",
)

entry(
    index = 56,
    label = "CH3 + H <=> CH4",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(
            A = (2.3e+14, 'cm^3/(mol*s)'),
            n = 0.032,
            Ea = (144, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        arrheniusLow = Arrhenius(
            A = (2.7e+35, 'cm^6/(mol^2*s)'),
            n = -5.345,
            Ea = (3380, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.395,
        T3 = (163.5, 'K'),
        T1 = (4250.3, 'K'),
        T2 = (1.25368e+06, 'K'),
        efficiencies = {'C': 3.85, 'O=C=O': 4, 'CC': 4.5, 'O': 10, '[H][H]': 4, '[He]': 2, 'N#N': 1.4, '[C-]#[O+]': 1.4, '[Ar]': 0.61},
    ),
    shortDesc = u"""The chemkin file reaction is CH3 + H <=> CH4""",
)

entry(
    index = 57,
    label = "CH4 + H <=> CH3 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4100, 'cm^3/(mol*s)'), n=3.156, Ea=(8755, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH4 + H <=> CH3 + H2""",
)

entry(
    index = 58,
    label = "CH4 + O <=> CH3 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(440000, 'cm^3/(mol*s)'), n=2.5, Ea=(6577, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH4 + O <=> CH3 + OH""",
)

entry(
    index = 59,
    label = "CH4 + OH <=> CH3 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+06, 'cm^3/(mol*s)'), n=2.182, Ea=(2506, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH4 + OH <=> CH3 + H2O""",
)

entry(
    index = 60,
    label = "CH4 + HO2 <=> CH3 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(47000, 'cm^3/(mol*s)'), n=2.5, Ea=(21000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH4 + HO2 <=> CH3 + H2O2""",
)

entry(
    index = 61,
    label = "CH4 + O2 <=> CH3 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (200000, 'cm^3/(mol*s)'),
        n = 2.745,
        Ea = (51714, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH4 + O2 <=> CH3 + HO2""",
)

entry(
    index = 62,
    label = "CH4 + CH2 <=> CH3 + CH3",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(4.3e+12, 'cm^3/(mol*s)'), n=0, Ea=(10030, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(4.3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH4 + CH2 <=> CH3 + CH3""",
)

entry(
    index = 63,
    label = "CH4 + CH <=> C2H4 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(-400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH4 + CH <=> C2H4 + H""",
)

entry(
    index = 64,
    label = "CH3 <=> CH + H2",
    degeneracy = 1,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(A=(3.1e+15, 'cm^3/(mol*s)'), n=0, Ea=(80871, 'cal/mol'), T0=(1, 'K')),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is CH3 <=> CH + H2""",
)

entry(
    index = 65,
    label = "CH3 <=> CH2 + H",
    degeneracy = 1,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(A=(2.2e+15, 'cm^3/(mol*s)'), n=0, Ea=(82659, 'cal/mol'), T0=(1, 'K')),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is CH3 <=> CH2 + H""",
)

entry(
    index = 66,
    label = "CH3 + H <=> CH2 + H2",
    degeneracy = 1,
    duplicate = True,
    kinetics = Arrhenius(
        A = (1.2e+06, 'cm^3/(mol*s)'),
        n = 2.43,
        Ea = (11941, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3 + H <=> CH2 + H2""",
)

entry(
    index = 67,
    label = "CH2(S) + H2 <=> CH3 + H",
    degeneracy = 1,
    duplicate = True,
    kinetics = Arrhenius(A=(7.2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2(S) + H2 <=> CH3 + H""",
)

entry(
    index = 68,
    label = "CH3 + O <=> CH2O + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.9e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + O <=> CH2O + H""",
)

entry(
    index = 69,
    label = "CH3 + O <=> H2 + CO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + O <=> H2 + CO + H""",
)

entry(
    index = 70,
    label = "CH3 + OH <=> CH2 + H2O",
    degeneracy = 1,
    duplicate = True,
    kinetics = Arrhenius(A=(43000, 'cm^3/(mol*s)'), n=2.568, Ea=(3997, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + OH <=> CH2 + H2O""",
)

entry(
    index = 71,
    label = "CH3 + OH <=> CH3OH",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.001316, 0.013158, 0.131579, 1.31579, 13.1579, 131.579], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (6.9e+30, 'cm^3/(mol*s)'),
                n = -6.63794,
                Ea = (2829, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.1e+32, 'cm^3/(mol*s)'),
                n = -6.63695,
                Ea = (3364, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.5e+32, 'cm^3/(mol*s)'),
                n = -6.36057,
                Ea = (3954, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (5.6e+30, 'cm^3/(mol*s)'),
                n = -5.64842,
                Ea = (4214, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.4e+27, 'cm^3/(mol*s)'),
                n = -4.33275,
                Ea = (3685, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.3e+22, 'cm^3/(mol*s)'),
                n = -2.66369,
                Ea = (2451, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH3 + OH <=> CH3OH""",
)

entry(
    index = 72,
    label = "CH3 + OH <=> CH2(S) + H2O",
    degeneracy = 1,
    duplicate = True,
    kinetics = PDepArrhenius(
        pressures = ([0.001316, 0.013158, 0.131579, 1.31579, 13.1579, 131.579], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (1.1e+14, 'cm^3/(mol*s)'),
                n = -0.45845,
                Ea = (-496, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.2e+14, 'cm^3/(mol*s)'),
                n = -0.53832,
                Ea = (-220, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.2e+15, 'cm^3/(mol*s)'),
                n = -0.72747,
                Ea = (600, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (4.3e+15, 'cm^3/(mol*s)'),
                n = -0.85972,
                Ea = (1888, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (4.4e+14, 'cm^3/(mol*s)'),
                n = -0.53864,
                Ea = (2932, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (6.1e+10, 'cm^3/(mol*s)'),
                n = 0.5956,
                Ea = (2923, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH3 + OH <=> CH2(S) + H2O""",
)

entry(
    index = 73,
    label = "CH3 + OH <=> H2 + CH2O",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.001316, 0.013158, 0.131579, 1.31579, 13.1579, 131.579], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (3.9e+09, 'cm^3/(mol*s)'),
                n = 0.25392,
                Ea = (-1221, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2e+10, 'cm^3/(mol*s)'),
                n = 0.06025,
                Ea = (-624, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.8e+11, 'cm^3/(mol*s)'),
                n = -0.24957,
                Ea = (498, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (3.6e+12, 'cm^3/(mol*s)'),
                n = -0.53245,
                Ea = (2042, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.2e+12, 'cm^3/(mol*s)'),
                n = -0.43166,
                Ea = (3415, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.4e+09, 'cm^3/(mol*s)'),
                n = 0.45344,
                Ea = (3791, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH3 + OH <=> H2 + CH2O""",
)

entry(
    index = 74,
    label = "CH3 + OH <=> H + CH2OH",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.001316, 0.013158, 0.131579, 1.31579, 13.1579, 131.579], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (8.4e+09, 'cm^3/(mol*s)'),
                n = 0.96279,
                Ea = (3230, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (8.4e+09, 'cm^3/(mol*s)'),
                n = 0.96279,
                Ea = (3230, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1e+10, 'cm^3/(mol*s)'),
                n = 0.94201,
                Ea = (3295, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (5.6e+10, 'cm^3/(mol*s)'),
                n = 0.73966,
                Ea = (3971, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (5.5e+11, 'cm^3/(mol*s)'),
                n = 0.4862,
                Ea = (5443, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.5e+10, 'cm^3/(mol*s)'),
                n = 0.9092,
                Ea = (6402, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH3 + OH <=> H + CH2OH""",
)

entry(
    index = 75,
    label = "CH3 + OH <=> H + CH3O",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.001316, 0.013158, 0.131579, 1.31579, 13.1579, 131.579], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (7.9e+08, 'cm^3/(mol*s)'),
                n = 1.06509,
                Ea = (11859, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (7.9e+08, 'cm^3/(mol*s)'),
                n = 1.06509,
                Ea = (11859, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (7.9e+08, 'cm^3/(mol*s)'),
                n = 1.06509,
                Ea = (11859, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (7.9e+08, 'cm^3/(mol*s)'),
                n = 1.06457,
                Ea = (11859, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1e+09, 'cm^3/(mol*s)'),
                n = 1.03413,
                Ea = (11970, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (3.1e+09, 'cm^3/(mol*s)'),
                n = 0.92189,
                Ea = (12981, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH3 + OH <=> H + CH3O""",
)

entry(
    index = 76,
    label = "CH3 + OH <=> H2 + HCOH",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.001316, 0.013158, 0.131579, 1.31579, 13.1579, 131.579], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (1.2e+09, 'cm^3/(mol*s)'),
                n = 0.83024,
                Ea = (-2323, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (6.4e+09, 'cm^3/(mol*s)'),
                n = 0.63305,
                Ea = (-1701, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (8e+10, 'cm^3/(mol*s)'),
                n = 0.33964,
                Ea = (-565, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (6.5e+11, 'cm^3/(mol*s)'),
                n = 0.11155,
                Ea = (932, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.1e+11, 'cm^3/(mol*s)'),
                n = 0.29509,
                Ea = (2200, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (9.4e+07, 'cm^3/(mol*s)'),
                n = 1.28631,
                Ea = (2424, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH3 + OH <=> H2 + HCOH""",
)

entry(
    index = 77,
    label = "CH3 + HO2 <=> CH3O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1e+12, 'cm^3/(mol*s)'),
        n = 0.2688,
        Ea = (-688, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3 + HO2 <=> CH3O + OH""",
)

entry(
    index = 78,
    label = "CH3 + O2 <=> CH3O + O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(28297, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + O2 <=> CH3O + O""",
)

entry(
    index = 79,
    label = "CH3 + O2 <=> CH2O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.9e+11, 'cm^3/(mol*s)'), n=0, Ea=(9842, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + O2 <=> CH2O + OH""",
)

entry(
    index = 80,
    label = "CH3 + O2 <=> CH3OO",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 1, 10, 20, 50, 100], 'atm'),
        arrhenius = [
            Arrhenius(A=(6.8e+24, 'cm^3/(mol*s)'), n=-3, Ea=(0, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(5e+22, 'cm^3/(mol*s)'), n=-3.85, Ea=(2000, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (3.4e+21, 'cm^3/(mol*s)'),
                n = -3.2,
                Ea = (2300, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (4.1e+20, 'cm^3/(mol*s)'),
                        n = -2.94,
                        Ea = (1900, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (3.3e+29, 'cm^3/(mol*s)'),
                        n = -5.6,
                        Ea = (6850, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (2.8e+18, 'cm^3/(mol*s)'),
                        n = -2.2,
                        Ea = (1400, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (5.6e+28, 'cm^3/(mol*s)'),
                        n = -5.25,
                        Ea = (6850, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (1.1e+19, 'cm^3/(mol*s)'),
                        n = -2.3,
                        Ea = (1800, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (4.1e+30, 'cm^3/(mol*s)'),
                        n = -5.7,
                        Ea = (8750, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH3 + O2 <=> CH3OO""",
)

entry(
    index = 81,
    label = "CH2O + CH3 <=> HCO + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (350000, 'cm^3/(mol*s)'),
        n = 2.157,
        Ea = (6234, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2O + CH3 <=> HCO + CH4""",
)

entry(
    index = 82,
    label = "CH2O + CH3 <=> H + CO + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.9e+11, 'cm^3/(mol*s)'),
        n = 0.887,
        Ea = (24224, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2O + CH3 <=> H + CO + CH4""",
)

entry(
    index = 83,
    label = "CH3 + HCO <=> CH4 + CO",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.05, 0.1, 1, 10], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (2.9e+18, 'cm^3/(mol*s)'),
                n = -1.84,
                Ea = (2134, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (8.7e+18, 'cm^3/(mol*s)'),
                n = -1.97,
                Ea = (2684, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.8e+20, 'cm^3/(mol*s)'),
                n = -2.3,
                Ea = (4781, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.1e+21, 'cm^3/(mol*s)'),
                n = -2.45,
                Ea = (7417, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH3 + HCO <=> CH4 + CO""",
)

entry(
    index = 84,
    label = "CH3 + HCO <=> CH2CO + H2",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.05, 0.1, 1, 10], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (6.1e+06, 'cm^3/(mol*s)'),
                n = 1.24,
                Ea = (-1733, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.1e+07, 'cm^3/(mol*s)'),
                n = 1.18,
                Ea = (-1303, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(A=(4.9e+08, 'cm^3/(mol*s)'), n=0.75, Ea=(842, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (1.6e+11, 'cm^3/(mol*s)'),
                n = 0.109,
                Ea = (4387, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH3 + HCO <=> CH2CO + H2""",
)

entry(
    index = 85,
    label = "CH3 + CH3 <=> C2H5 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(16055, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + CH3 <=> C2H5 + H""",
)

entry(
    index = 86,
    label = "CH3 + CH2 <=> C2H4 + H",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(
                A = (1.2e+15, 'cm^3/(mol*s)'),
                n = -0.3432,
                Ea = (153, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH3 + CH2 <=> C2H4 + H""",
)

entry(
    index = 87,
    label = "CH3 + CH <=> C2H3 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + CH <=> C2H3 + H""",
)

entry(
    index = 88,
    label = "CH3 + C <=> C2H2 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + C <=> C2H2 + H""",
)

entry(
    index = 89,
    label = "CH2 <=> CH + H",
    degeneracy = 1,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(A=(5.6e+15, 'cm^3/(mol*s)'), n=0, Ea=(89000, 'cal/mol'), T0=(1, 'K')),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is CH2 <=> CH + H""",
)

entry(
    index = 90,
    label = "CH2 <=> C + H2",
    degeneracy = 1,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(
            A = (5.8e+12, 'cm^3/(mol*s)'),
            n = 0.5,
            Ea = (68500, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is CH2 <=> C + H2""",
)

entry(
    index = 91,
    label = "CH2 + H <=> CH + H2",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(1.6e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2 + H <=> CH + H2""",
)

entry(
    index = 92,
    label = "CH2 + O <=> CO + H + H",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(1.2e+14, 'cm^3/(mol*s)'), n=0, Ea=(536, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2 + O <=> CO + H + H""",
)

entry(
    index = 93,
    label = "CH2 + O <=> CO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+13, 'cm^3/(mol*s)'), n=0, Ea=(536, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2 + O <=> CO + H2""",
)

entry(
    index = 94,
    label = "CH2 + OH <=> CH2O + H",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(
                A = (2.8e+13, 'cm^3/(mol*s)'),
                n = 0.1228,
                Ea = (-161, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2 + OH <=> CH2O + H""",
)

entry(
    index = 95,
    label = "CH2 + OH <=> CH + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (860000, 'cm^3/(mol*s)'),
        n = 2.019,
        Ea = (6776, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2 + OH <=> CH + H2O""",
)

entry(
    index = 96,
    label = "CH2 + O2 <=> CO2 + H + H",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.1e+09, 'cm^3/(mol*s)'),
        n = 0.9929,
        Ea = (-269, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2 + O2 <=> CO2 + H + H""",
)

entry(
    index = 97,
    label = "CH2 + O2 <=> CH2O + O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.2e+09, 'cm^3/(mol*s)'),
        n = 1.08,
        Ea = (1196, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2 + O2 <=> CH2O + O""",
)

entry(
    index = 98,
    label = "CH2 + CO2 <=> CO + CH2O",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(1000, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2 + CO2 <=> CO + CH2O""",
)

entry(
    index = 99,
    label = "CH2 + CH2 <=> C2H2 + H + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+13, 'cm^3/(mol*s)'), n=0.0022, Ea=(8, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2 + CH2 <=> C2H2 + H + H""",
)

entry(
    index = 100,
    label = "CH2 + CH2 <=> C2H2 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+13, 'cm^3/(mol*s)'), n=0.0022, Ea=(8, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2 + CH2 <=> C2H2 + H2""",
)

entry(
    index = 101,
    label = "CH2 + CH <=> C2H2 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2 + CH <=> C2H2 + H""",
)

entry(
    index = 102,
    label = "CH2 + C <=> C2H + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2 + C <=> C2H + H""",
)

entry(
    index = 103,
    label = "CH2(S) + N2 <=> CH2 + N2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+13, 'cm^3/(mol*s)'), n=0, Ea=(430, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2(S) + N2 <=> CH2 + N2""",
)

entry(
    index = 104,
    label = "CH2(S) + AR <=> CH2 + AR",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.5e+13, 'cm^3/(mol*s)'), n=0, Ea=(884, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2(S) + AR <=> CH2 + AR""",
)

entry(
    index = 105,
    label = "CH2(S) + H <=> CH2 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2(S) + H <=> CH2 + H""",
)

entry(
    index = 106,
    label = "CH2(S) + O2 <=> CH2 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2(S) + O2 <=> CH2 + O2""",
)

entry(
    index = 107,
    label = "CH2(S) + H2O <=> CH2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2(S) + H2O <=> CH2 + H2O""",
)

entry(
    index = 108,
    label = "CH + H <=> C + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH + H <=> C + H2""",
)

entry(
    index = 109,
    label = "CH + O <=> CO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.7e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH + O <=> CO + H""",
)

entry(
    index = 110,
    label = "CH + OH <=> HCO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH + OH <=> HCO + H""",
)

entry(
    index = 111,
    label = "CH + OH <=> H + CO + H",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.8e+23, 'cm^3/(mol*s)'),
        n = -2.473,
        Ea = (19927, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH + OH <=> H + CO + H""",
)

entry(
    index = 112,
    label = "CH + OH <=> C + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+07, 'cm^3/(mol*s)'), n=2, Ea=(3000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH + OH <=> C + H2O""",
)

entry(
    index = 113,
    label = "CH + O2 <=> HCO + O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.7e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH + O2 <=> HCO + O""",
)

entry(
    index = 114,
    label = "CH + O2 <=> H + CO + O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2e+23, 'cm^3/(mol*s)'),
        n = -2.473,
        Ea = (19927, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH + O2 <=> H + CO + O""",
)

entry(
    index = 115,
    label = "CH + H2O <=> CH2O + H",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.8e+07, 'cm^3/(mol*s)'),
        n = 1.59,
        Ea = (-2610, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH + H2O <=> CH2O + H""",
)

entry(
    index = 116,
    label = "CH + CO2 <=> HCO + CO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8.8e+06, 'cm^3/(mol*s)'),
        n = 1.75,
        Ea = (-1040, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH + CO2 <=> HCO + CO""",
)

entry(
    index = 117,
    label = "CH + CO2 <=> H + CO + CO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.3e+16, 'cm^3/(mol*s)'),
        n = -0.723,
        Ea = (18887, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH + CO2 <=> H + CO + CO""",
)

entry(
    index = 118,
    label = "CH + CH2O <=> CH2CO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.5e+13, 'cm^3/(mol*s)'), n=0, Ea=(-517, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH + CH2O <=> CH2CO + H""",
)

entry(
    index = 119,
    label = "C + OH <=> CO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C + OH <=> CO + H""",
)

entry(
    index = 120,
    label = "C + O2 <=> CO + O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C + O2 <=> CO + O""",
)

entry(
    index = 121,
    label = "CH3OH <=> CH2(S) + H2O",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.1e+18, 's^-1'), n=-1.017, Ea=(91712, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.4e+23, 'cm^3/(mol*s)'),
            n = -8.3446,
            Ea = (99596, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.9922,
        T3 = (943, 'K'),
        T1 = (47310, 'K'),
        T2 = (47110, 'K'),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is CH3OH <=> CH2(S) + H2O""",
)

entry(
    index = 122,
    label = "CH2OH + H <=> CH3OH",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(4.3e+15, 'cm^3/(mol*s)'), n=-0.79, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.844e+37, 'cm^6/(mol^2*s)'),
            n = -6.21,
            Ea = (1333, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.25,
        T3 = (210, 'K'),
        T1 = (1434, 'K'),
        T2 = (1e+30, 'K'),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is CH2OH + H <=> CH3OH""",
)

entry(
    index = 123,
    label = "CH3O + H <=> CH3OH",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2.4e+12, 'cm^3/(mol*s)'), n=0.515, Ea=(50, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4.66e+41, 'cm^6/(mol^2*s)'),
            n = -7.44,
            Ea = (14080, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.7,
        T3 = (100, 'K'),
        T1 = (90000, 'K'),
        T2 = (10000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, 'N#N': 1, '[C-]#[O+]': 1.5},
    ),
    shortDesc = u"""The chemkin file reaction is CH3O + H <=> CH3OH""",
)

entry(
    index = 124,
    label = "CH3OH + H <=> CH2OH + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(66000, 'cm^3/(mol*s)'), n=2.728, Ea=(4449, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OH + H <=> CH2OH + H2""",
)

entry(
    index = 125,
    label = "CH3OH + H <=> CH3O + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(41000, 'cm^3/(mol*s)'), n=2.658, Ea=(9221, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OH + H <=> CH3O + H2""",
)

entry(
    index = 126,
    label = "CH3OH + O <=> CH2OH + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.1e+13, 'cm^3/(mol*s)'), n=0, Ea=(5305, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OH + O <=> CH2OH + OH""",
)

entry(
    index = 127,
    label = "CH3OH + O <=> CH3O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.7e+12, 'cm^3/(mol*s)'), n=0, Ea=(5305, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OH + O <=> CH3O + OH""",
)

entry(
    index = 128,
    label = "CH3OH + OH <=> CH2OH + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.5e+08, 'cm^3/(mol*s)'),
        n = 1.4434,
        Ea = (113, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3OH + OH <=> CH2OH + H2O""",
)

entry(
    index = 129,
    label = "CH3OH + OH <=> CH3O + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.7e+07, 'cm^3/(mol*s)'),
        n = 1.4434,
        Ea = (113, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3OH + OH <=> CH3O + H2O""",
)

entry(
    index = 130,
    label = "CH3OH + HO2 <=> CH2OH + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (0.00035, 'cm^3/(mol*s)'),
        n = 4.85,
        Ea = (10346, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3OH + HO2 <=> CH2OH + H2O2""",
)

entry(
    index = 131,
    label = "CH3OH + HO2 <=> CH3O + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (0.0015, 'cm^3/(mol*s)'),
        n = 4.61,
        Ea = (15928, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3OH + HO2 <=> CH3O + H2O2""",
)

entry(
    index = 132,
    label = "CH3OH + O2 <=> CH2OH + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (360000, 'cm^3/(mol*s)'),
        n = 2.27,
        Ea = (42760, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3OH + O2 <=> CH2OH + HO2""",
)

entry(
    index = 133,
    label = "CH3O + HO2 <=> CH3OH + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.4e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O + HO2 <=> CH3OH + O2""",
)

entry(
    index = 134,
    label = "CH3OH + CH3 <=> CH3O + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.013, 'cm^3/(mol*s)'), n=4.161, Ea=(6002, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OH + CH3 <=> CH3O + CH4""",
)

entry(
    index = 135,
    label = "CH3OH + CH3 <=> CH2OH + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.74, 'cm^3/(mol*s)'), n=3.781, Ea=(7183, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OH + CH3 <=> CH2OH + CH4""",
)

entry(
    index = 136,
    label = "CH2OH <=> CH2O + H",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(7.4e+10, 's^-1'), n=0.811, Ea=(39559, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.5e+21, 'cm^3/(mol*s)'),
            n = -1.99,
            Ea = (23983, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.844,
        T3 = (900, 'K'),
        T1 = (1, 'K'),
        T2 = (3315, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'O': 6, '[H][H]': 2, '[He]': 0.67, '[O][O]': 1, '[C-]#[O+]': 1.5, '[Ar]': 0.85},
    ),
    shortDesc = u"""The chemkin file reaction is CH2OH <=> CH2O + H""",
)

entry(
    index = 137,
    label = "CH2OH + H <=> CH2O + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+06, 'cm^3/(mol*s)'), n=1.86, Ea=(147, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2OH + H <=> CH2O + H2""",
)

entry(
    index = 138,
    label = "CH2OH + O <=> CH2O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(-693, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2OH + O <=> CH2O + OH""",
)

entry(
    index = 139,
    label = "CH2OH + OH <=> CH2O + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2OH + OH <=> CH2O + H2O""",
)

entry(
    index = 140,
    label = "CH2OH + HO2 <=> CH2O + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2OH + HO2 <=> CH2O + H2O2""",
)

entry(
    index = 141,
    label = "CH2OH + O2 <=> CH2O + HO2",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(7.2e+13, 'cm^3/(mol*s)'), n=0, Ea=(3736, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2.9e+16, 'cm^3/(mol*s)'), n=-1.5, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2OH + O2 <=> CH2O + HO2""",
)

entry(
    index = 142,
    label = "CH2OH + CH3 <=> C2H4 + H2O",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.001, 0.01, 0.1, 1, 10, 100], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (6.3e+24, 'cm^3/(mol*s)'),
                n = -3.7134,
                Ea = (2798, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.2e+25, 'cm^3/(mol*s)'),
                n = -3.7867,
                Ea = (3001, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (3.2e+27, 'cm^3/(mol*s)'),
                n = -4.45,
                Ea = (5345, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (7.2e+29, 'cm^3/(mol*s)'),
                n = -5.0344,
                Ea = (9245, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.7e+27, 'cm^3/(mol*s)'),
                n = -4.1839,
                Ea = (11152, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (3.9e+17, 'cm^3/(mol*s)'),
                n = -1.3688,
                Ea = (8978, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2OH + CH3 <=> C2H4 + H2O""",
)

entry(
    index = 143,
    label = "CH2OH + HCO <=> CH3OH + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2OH + HCO <=> CH3OH + CO""",
)

entry(
    index = 144,
    label = "CH2OH + HCO <=> CH2O + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2OH + HCO <=> CH2O + CH2O""",
)

entry(
    index = 145,
    label = "CH2OH + CH2O <=> CH3OH + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5500, 'cm^3/(mol*s)'), n=2.81, Ea=(5862, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2OH + CH2O <=> CH3OH + HCO""",
)

entry(
    index = 146,
    label = "CH2OH + CH2O <=> CH3OH + H + CO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.3e+13, 'cm^3/(mol*s)'),
        n = 0.337,
        Ea = (25789, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2OH + CH2O <=> CH3OH + H + CO""",
)

entry(
    index = 147,
    label = "CH2OH + CH2OH <=> CH3OH + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2OH + CH2OH <=> CH3OH + CH2O""",
)

entry(
    index = 148,
    label = "CH2OH + CH3O <=> CH3OH + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2OH + CH3O <=> CH3OH + CH2O""",
)

entry(
    index = 149,
    label = "CH3O <=> CH2O + H",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.1e+10, 's^-1'), n=1.21, Ea=(24069, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (6e+16, 'cm^3/(mol*s)'),
            n = -0.547,
            Ea = (18012, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.341,
        T3 = (28, 'K'),
        T1 = (1000, 'K'),
        T2 = (2339, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'O': 6, '[H][H]': 2, '[He]': 0.67, '[O][O]': 1, '[C-]#[O+]': 1.5, '[Ar]': 0.85},
    ),
    shortDesc = u"""The chemkin file reaction is CH3O <=> CH2O + H""",
)

entry(
    index = 150,
    label = "CH3O + H <=> CH2O + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.6e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(-519, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O + H <=> CH2O + H2""",
)

entry(
    index = 151,
    label = "CH3O + O <=> CH2O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O + O <=> CH2O + OH""",
)

entry(
    index = 152,
    label = "CH3O + OH <=> CH2O + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O + OH <=> CH2O + H2O""",
)

entry(
    index = 153,
    label = "CH3O + HO2 <=> CH2O + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O + HO2 <=> CH2O + H2O2""",
)

entry(
    index = 154,
    label = "CH3O + O2 <=> CH2O + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.48, 'cm^3/(mol*s)'), n=3.567, Ea=(-1055, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O + O2 <=> CH2O + HO2""",
)

entry(
    index = 155,
    label = "CH3O + CO <=> CH3 + CO2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (9.5e+25, 'cm^3/(mol*s)'),
        n = -4.93,
        Ea = (9080, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3O + CO <=> CH3 + CO2""",
)

entry(
    index = 156,
    label = "CH3O + CH3 <=> CH2O + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O + CH3 <=> CH2O + CH4""",
)

entry(
    index = 157,
    label = "CH3O + CH2O <=> CH3OH + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.6e+11, 'cm^3/(mol*s)'), n=0, Ea=(2285, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O + CH2O <=> CH3OH + HCO""",
)

entry(
    index = 158,
    label = "CH3O + CH3O <=> CH3OH + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O + CH3O <=> CH3OH + CH2O""",
)

entry(
    index = 159,
    label = "CH3OOH <=> CH3O + OH",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(4.1e+19, 's^-1'), n=-1.153, Ea=(44226, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.9e+42, 'cm^3/(mol*s)'),
            n = -7.502,
            Ea = (46730, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.8375,
        T3 = (36562, 'K'),
        T1 = (498.8, 'K'),
        T2 = (9990, 'K'),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is CH3OOH <=> CH3O + OH""",
)

entry(
    index = 160,
    label = "CH3OOH + H <=> CH2OOH + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.4e+10, 'cm^3/(mol*s)'), n=0, Ea=(1860, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OOH + H <=> CH2OOH + H2""",
)

entry(
    index = 161,
    label = "CH3OOH + H <=> CH3OO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.4e+10, 'cm^3/(mol*s)'), n=0, Ea=(1860, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OOH + H <=> CH3OO + H2""",
)

entry(
    index = 162,
    label = "CH3OOH + H <=> CH3O + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+10, 'cm^3/(mol*s)'), n=0, Ea=(1860, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OOH + H <=> CH3O + H2O""",
)

entry(
    index = 163,
    label = "CH3OOH + O <=> CH2OOH + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(4750, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OOH + O <=> CH2OOH + OH""",
)

entry(
    index = 164,
    label = "CH3OOH + O <=> CH3OO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.7e+12, 'cm^3/(mol*s)'), n=0, Ea=(4750, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OOH + O <=> CH3OO + OH""",
)

entry(
    index = 165,
    label = "CH3OOH + OH <=> CH3OO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.1e+12, 'cm^3/(mol*s)'), n=0, Ea=(-437, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OOH + OH <=> CH3OO + H2O""",
)

entry(
    index = 166,
    label = "CH3OOH + OH <=> CH2OOH + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.2e+11, 'cm^3/(mol*s)'), n=0, Ea=(-258, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OOH + OH <=> CH2OOH + H2O""",
)

entry(
    index = 167,
    label = "CH3OOH + HO2 <=> CH3OO + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(41000, 'cm^3/(mol*s)'), n=2.5, Ea=(10206, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OOH + HO2 <=> CH3OO + H2O2""",
)

entry(
    index = 168,
    label = "CH3OO + H <=> CH3O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OO + H <=> CH3O + OH""",
)

entry(
    index = 169,
    label = "CH3OO + O <=> CH3O + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.9e+10, 'cm^3/(mol*s)'), n=1, Ea=(-724, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OO + O <=> CH3O + O2""",
)

entry(
    index = 170,
    label = "CH3OO + OH <=> CH3OH + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OO + OH <=> CH3OH + O2""",
)

entry(
    index = 171,
    label = "CH3OO + HO2 <=> CH3OOH + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.5e+11, 'cm^3/(mol*s)'), n=0, Ea=(-1490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OO + HO2 <=> CH3OOH + O2""",
)

entry(
    index = 172,
    label = "CH3OO + CH3 <=> CH3O + CH3O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.1e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1411, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OO + CH3 <=> CH3O + CH3O""",
)

entry(
    index = 173,
    label = "CH3OOH + CH3 <=> CH3OO + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (0.0016, 'cm^3/(mol*s)'),
        n = 4.322,
        Ea = (-235, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3OOH + CH3 <=> CH3OO + CH4""",
)

entry(
    index = 174,
    label = "CH3OO + CH2OH <=> CH2O + CH3OOH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OO + CH2OH <=> CH2O + CH3OOH""",
)

entry(
    index = 175,
    label = "CH3OO + HCO <=> CH3O + H + CO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OO + HCO <=> CH3O + H + CO2""",
)

entry(
    index = 176,
    label = "CH3OO + CO <=> CH3O + CO2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (160000, 'cm^3/(mol*s)'),
        n = 2.18,
        Ea = (17940, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3OO + CO <=> CH3O + CO2""",
)

entry(
    index = 177,
    label = "CH3OO + CH2O <=> CH3OOH + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(41000, 'cm^3/(mol*s)'), n=2.5, Ea=(10206, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OO + CH2O <=> CH3OOH + HCO""",
)

entry(
    index = 178,
    label = "CH3OO + CH2O <=> CH3OOH + H + CO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.5e+14, 'cm^3/(mol*s)'),
        n = 0.027,
        Ea = (30133, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3OO + CH2O <=> CH3OOH + H + CO""",
)

entry(
    index = 179,
    label = "CH3OO + CH3O <=> CH2O + CH3OOH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OO + CH3O <=> CH2O + CH3OOH""",
)

entry(
    index = 180,
    label = "CH3OO + CH3OH <=> CH3OOH + CH2OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(19400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OO + CH3OH <=> CH3OOH + CH2OH""",
)

entry(
    index = 181,
    label = "CH3OO + CH3OO <=> CH3O + CH3O + O2",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(
                A = (1.1e+18, 'cm^3/(mol*s)'),
                n = -2.4,
                Ea = (1800, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(A=(7e+10, 'cm^3/(mol*s)'), n=0, Ea=(800, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH3OO + CH3OO <=> CH3O + CH3O + O2""",
)

entry(
    index = 182,
    label = "CH3OO + CH3OO <=> CH3OH + CH2O + O2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2e+11, 'cm^3/(mol*s)'),
        n = -0.55,
        Ea = (-1600, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3OO + CH3OO <=> CH3OH + CH2O + O2""",
)

entry(
    index = 183,
    label = "CH2OOH <=> CH2O + OH",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.04, 1, 10, 100], 'atm'),
        arrhenius = [
            Arrhenius(A=(9.6e+10, 's^-1'), n=-0.925, Ea=(1567, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2.4e+12, 's^-1'), n=-0.925, Ea=(1567, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2.5e+13, 's^-1'), n=-0.927, Ea=(1579, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(7e+14, 's^-1'), n=-1.064, Ea=(1744, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2OOH <=> CH2O + OH""",
)

entry(
    index = 184,
    label = "HCOH <=> CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.2e+13, 's^-1'), n=0, Ea=(32109, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCOH <=> CH2O""",
)

entry(
    index = 185,
    label = "HCOH <=> CO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.6e+13, 's^-1'), n=0, Ea=(49465, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCOH <=> CO + H2""",
)

entry(
    index = 186,
    label = "HCOH + H <=> CH2O + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCOH + H <=> CH2O + H""",
)

entry(
    index = 187,
    label = "HCOH + H <=> HCO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCOH + H <=> HCO + H2""",
)

entry(
    index = 188,
    label = "HCOH + H <=> H + CO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCOH + H <=> H + CO + H2""",
)

entry(
    index = 189,
    label = "HCOH + O <=> CO2 + H + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCOH + O <=> CO2 + H + H""",
)

entry(
    index = 190,
    label = "HCOH + OH <=> HCO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCOH + OH <=> HCO + H2O""",
)

entry(
    index = 191,
    label = "HCOH + OH <=> H + CO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCOH + OH <=> H + CO + H2O""",
)

entry(
    index = 192,
    label = "CH2O + OH <=> HOCH2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (6.3e+06, 'cm^3/(mol*s)'),
        n = 1.63,
        Ea = (4282, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2O + OH <=> HOCH2O""",
)

entry(
    index = 193,
    label = "HOCH2O <=> HOCHO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 's^-1'), n=0, Ea=(14900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HOCH2O <=> HOCHO + H""",
)

entry(
    index = 194,
    label = "CH3 + CH3 <=> C2H6",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(
            A = (9.5e+14, 'cm^3/(mol*s)'),
            n = -0.538,
            Ea = (179, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        arrheniusLow = Arrhenius(
            A = (1.269e+41, 'cm^6/(mol^2*s)'),
            n = -7,
            Ea = (2762, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.62,
        T3 = (73, 'K'),
        T1 = (1180, 'K'),
        T2 = (1e+30, 'K'),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is CH3 + CH3 <=> C2H6""",
)

entry(
    index = 195,
    label = "C2H6 + H <=> C2H5 + H2",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(7400, 'cm^3/(mol*s)'), n=3.1, Ea=(5340, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(3.3e+14, 'cm^3/(mol*s)'), n=0, Ea=(13667, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H6 + H <=> C2H5 + H2""",
)

entry(
    index = 196,
    label = "C2H6 + O <=> C2H5 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(180000, 'cm^3/(mol*s)'), n=2.8, Ea=(5800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H6 + O <=> C2H5 + OH""",
)

entry(
    index = 197,
    label = "C2H6 + OH <=> C2H5 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.6e+06, 'cm^3/(mol*s)'),
        n = 2.224,
        Ea = (741, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H6 + OH <=> C2H5 + H2O""",
)

entry(
    index = 198,
    label = "C2H6 + HO2 <=> C2H5 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(87000, 'cm^3/(mol*s)'), n=2.65, Ea=(18900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H6 + HO2 <=> C2H5 + H2O2""",
)

entry(
    index = 199,
    label = "C2H6 + O2 <=> C2H5 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.9e+07, 'cm^3/(mol*s)'),
        n = 1.9,
        Ea = (49548, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H6 + O2 <=> C2H5 + HO2""",
)

entry(
    index = 200,
    label = "C2H6 + CH3 <=> C2H5 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(35, 'cm^3/(mol*s)'), n=3.44, Ea=(10384, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H6 + CH3 <=> C2H5 + CH4""",
)

entry(
    index = 201,
    label = "C2H6 + CH2(S) <=> C2H5 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H6 + CH2(S) <=> C2H5 + CH3""",
)

entry(
    index = 202,
    label = "C2H6 + CH3OO <=> C2H5 + CH3OOH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(19, 'cm^3/(mol*s)'), n=3.64, Ea=(17100, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H6 + CH3OO <=> C2H5 + CH3OOH""",
)

entry(
    index = 203,
    label = "C2H4 + H <=> C2H5",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(
            A = (1.4e+09, 'cm^3/(mol*s)'),
            n = 1.463,
            Ea = (1355, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        arrheniusLow = Arrhenius(
            A = (2e+39, 'cm^6/(mol^2*s)'),
            n = -6.642,
            Ea = (5769, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -0.569,
        T3 = (299, 'K'),
        T1 = (9147, 'K'),
        T2 = (152.4, 'K'),
        efficiencies = {'[C-]#[O+]': 1.5, '[H][H]': 2, 'O=C=O': 3, 'O': 10, 'N#N': 1.2},
    ),
    shortDesc = u"""The chemkin file reaction is C2H4 + H <=> C2H5""",
)

entry(
    index = 204,
    label = "C2H5 + H <=> C2H6",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(
            A = (5.2e+17, 'cm^3/(mol*s)'),
            n = -0.99,
            Ea = (1580, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        arrheniusLow = Arrhenius(
            A = (1.99e+41, 'cm^6/(mol^2*s)'),
            n = -7.08,
            Ea = (6685, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.8422,
        T3 = (125, 'K'),
        T1 = (2219, 'K'),
        T2 = (6882, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, 'N#N': 1, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C2H5 + H <=> C2H6""",
)

entry(
    index = 205,
    label = "C2H5 + O <=> CH3 + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5 + O <=> CH3 + CH2O""",
)

entry(
    index = 206,
    label = "C2H5 + O <=> CH3CHO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5 + O <=> CH3CHO + H""",
)

entry(
    index = 207,
    label = "C2H5 + O <=> C2H4 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5 + O <=> C2H4 + OH""",
)

entry(
    index = 208,
    label = "C2H5 + OH <=> C2H4 + H2O",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.001, 0.01, 0.1, 1, 10, 100], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (1.3e+19, 'cm^3/(mol*s)'),
                n = -1.96,
                Ea = (273, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.2e+19, 'cm^3/(mol*s)'),
                n = -1.9533,
                Ea = (239, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (4.1e+19, 'cm^3/(mol*s)'),
                n = -2.1007,
                Ea = (625, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (7.9e+22, 'cm^3/(mol*s)'),
                n = -2.9892,
                Ea = (3863, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.8e+24, 'cm^3/(mol*s)'),
                n = -3.3287,
                Ea = (7749, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (4.7e+18, 'cm^3/(mol*s)'),
                n = -1.5805,
                Ea = (7999, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H5 + OH <=> C2H4 + H2O""",
)

entry(
    index = 209,
    label = "C2H5 + OH <=> CH3 + CH2OH",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.001, 0.01, 0.1, 1, 10, 100], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (9.2e+17, 'cm^3/(mol*s)'),
                n = -1.2994,
                Ea = (2505, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.1e+18, 'cm^3/(mol*s)'),
                n = -1.3206,
                Ea = (2569, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (5.7e+18, 'cm^3/(mol*s)'),
                n = -1.5182,
                Ea = (3185, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (6.5e+21, 'cm^3/(mol*s)'),
                n = -2.3515,
                Ea = (6023, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.9e+25, 'cm^3/(mol*s)'),
                n = -3.2495,
                Ea = (10576, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (6.5e+22, 'cm^3/(mol*s)'),
                n = -2.4427,
                Ea = (12647, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H5 + OH <=> CH3 + CH2OH""",
)

entry(
    index = 210,
    label = "C2H5 + HO2 <=> CH3CH2O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5 + HO2 <=> CH3CH2O + OH""",
)

entry(
    index = 211,
    label = "C2H5 + O2 <=> CH3CH2OO",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (1.3e+42, 'cm^3/(mol*s)'),
                n = -11.12,
                Ea = (5137, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(A=(2e+43, 'cm^3/(mol*s)'), n=-11.3, Ea=(5485, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (1.2e+44, 'cm^3/(mol*s)'),
                n = -11.36,
                Ea = (5850, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (3.2e+44, 'cm^3/(mol*s)'),
                n = -11.32,
                Ea = (6198, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.3e+45, 'cm^3/(mol*s)'),
                n = -11.33,
                Ea = (6761, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.2e+45, 'cm^3/(mol*s)'),
                n = -11.15,
                Ea = (7163, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (4.3e+44, 'cm^3/(mol*s)'),
                n = -10.83,
                Ea = (7564, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (4.9e+43, 'cm^3/(mol*s)'),
                n = -10.37,
                Ea = (7810, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(A=(4e+42, 'cm^3/(mol*s)'), n=-9.86, Ea=(8124, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (1.2e+40, 'cm^3/(mol*s)'),
                n = -8.95,
                Ea = (7857, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.2e+37, 'cm^3/(mol*s)'),
                n = -7.95,
                Ea = (7525, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.6e+34, 'cm^3/(mol*s)'),
                n = -6.88,
                Ea = (6913, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.6e+30, 'cm^3/(mol*s)'),
                n = -5.56,
                Ea = (5909, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H5 + O2 <=> CH3CH2OO""",
)

entry(
    index = 212,
    label = "C2H5 + O2 <=> C2H4 + HO2",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100], 'atm'),
        arrhenius = [
            Arrhenius(A=(9.1, 'cm^3/(mol*s)'), n=2.87, Ea=(-5099, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(12, 'cm^3/(mol*s)'), n=2.84, Ea=(-5029, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(28, 'cm^3/(mol*s)'), n=2.73, Ea=(-4780, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(110, 'cm^3/(mol*s)'), n=2.56, Ea=(-4380, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(960, 'cm^3/(mol*s)'), n=2.3, Ea=(-3735, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(13000, 'cm^3/(mol*s)'), n=1.98, Ea=(-2933, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (490000, 'cm^3/(mol*s)'),
                n = 1.54,
                Ea = (-1790, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.4e+07, 'cm^3/(mol*s)'),
                n = 1.07,
                Ea = (-497.7, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.5e+09, 'cm^3/(mol*s)'),
                n = 0.51,
                Ea = (1157, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.4e+11, 'cm^3/(mol*s)'),
                n = 0.04,
                Ea = (2789, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (3.1e+12, 'cm^3/(mol*s)'),
                n = -0.31,
                Ea = (4501, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (5.3e+12, 'cm^3/(mol*s)'),
                n = -0.33,
                Ea = (5728, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.9e+11, 'cm^3/(mol*s)'),
                n = 0.14,
                Ea = (6373, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H5 + O2 <=> C2H4 + HO2""",
)

entry(
    index = 213,
    label = "C2H5 + O2 <=> CH2CH2OOH",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (3.2e+21, 'cm^3/(mol*s)'),
                n = -5.53,
                Ea = (-83.5, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (4.6e+23, 'cm^3/(mol*s)'),
                n = -6.12,
                Ea = (586.3, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.9e+25, 'cm^3/(mol*s)'),
                n = -6.6,
                Ea = (1279, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.9e+24, 'cm^3/(mol*s)'),
                n = -6.19,
                Ea = (1229, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (5.9e+25, 'cm^3/(mol*s)'),
                n = -6.49,
                Ea = (2026, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.9e+25, 'cm^3/(mol*s)'),
                n = -6.26,
                Ea = (2449, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (4.6e+26, 'cm^3/(mol*s)'),
                n = -6.47,
                Ea = (3598, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (4.7e+26, 'cm^3/(mol*s)'),
                n = -6.29,
                Ea = (4518, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (3.6e+26, 'cm^3/(mol*s)'),
                n = -6.03,
                Ea = (5715, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (6.3e+25, 'cm^3/(mol*s)'),
                n = -5.58,
                Ea = (6793, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (6.6e+23, 'cm^3/(mol*s)'),
                n = -4.74,
                Ea = (7756, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (7.7e+20, 'cm^3/(mol*s)'),
                n = -3.63,
                Ea = (8319, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.9e+15, 'cm^3/(mol*s)'),
                n = -1.72,
                Ea = (8034, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H5 + O2 <=> CH2CH2OOH""",
)

entry(
    index = 214,
    label = "C2H5 + O2 <=> cC2H4O + OH",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (1.4e-05, 'cm^3/(mol*s)'),
                n = 4.2,
                Ea = (-5618, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.9e-05, 'cm^3/(mol*s)'),
                n = 4.16,
                Ea = (-5537, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(A=(5e-05, 'cm^3/(mol*s)'), n=4.04, Ea=(-5260, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (0.00023, 'cm^3/(mol*s)'),
                n = 3.85,
                Ea = (-4825, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (0.0026, 'cm^3/(mol*s)'),
                n = 3.55,
                Ea = (-4121, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(A=(0.052, 'cm^3/(mol*s)'), n=3.18, Ea=(-3238, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(4.1, 'cm^3/(mol*s)'), n=2.65, Ea=(-1928, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(610, 'cm^3/(mol*s)'), n=2.04, Ea=(-370.9, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(490000, 'cm^3/(mol*s)'), n=1.22, Ea=(1802, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(5e+08, 'cm^3/(mol*s)'), n=0.39, Ea=(4218, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (8.8e+11, 'cm^3/(mol*s)'),
                n = -0.49,
                Ea = (7190, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.7e+14, 'cm^3/(mol*s)'),
                n = -1.09,
                Ea = (9936, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (9.7e+14, 'cm^3/(mol*s)'),
                n = -1.22,
                Ea = (12500, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H5 + O2 <=> cC2H4O + OH""",
)

entry(
    index = 215,
    label = "C2H6 + HCO <=> C2H5 + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.3, 'cm^3/(mol*s)'), n=3.74, Ea=(16933, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H6 + HCO <=> C2H5 + CH2O""",
)

entry(
    index = 216,
    label = "C2H5 + HCO <=> C2H6 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5 + HCO <=> C2H6 + CO""",
)

entry(
    index = 217,
    label = "C2H5 + HCO <=> CH3 + CH2CHO",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.001, 0.01, 0.1, 1, 10, 100], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (9.2e+17, 'cm^3/(mol*s)'),
                n = -1.2994,
                Ea = (2505, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.1e+18, 'cm^3/(mol*s)'),
                n = -1.3206,
                Ea = (2569, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (5.7e+18, 'cm^3/(mol*s)'),
                n = -1.5182,
                Ea = (3185, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (6.5e+21, 'cm^3/(mol*s)'),
                n = -2.3515,
                Ea = (6023, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.9e+25, 'cm^3/(mol*s)'),
                n = -3.2495,
                Ea = (10576, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (6.5e+22, 'cm^3/(mol*s)'),
                n = -2.4427,
                Ea = (12647, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H5 + HCO <=> CH3 + CH2CHO""",
)

entry(
    index = 218,
    label = "C2H5 + CH3 <=> C2H4 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5 + CH3 <=> C2H4 + CH4""",
)

entry(
    index = 219,
    label = "C2H5 + CH3OO <=> CH3CH2O + CH3O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.1e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1410, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5 + CH3OO <=> CH3CH2O + CH3O""",
)

entry(
    index = 220,
    label = "C2H5 + C2H5 <=> C2H6 + C2H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5 + C2H5 <=> C2H6 + C2H4""",
)

entry(
    index = 221,
    label = "C2H3 + H <=> C2H4",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.9e+13, 'cm^3/(mol*s)'), n=0.2, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(A=(2.1e+24, 'cm^6/(mol^2*s)'), n=-1.3, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        alpha = 0.5,
        T3 = (1e-30, 'K'),
        T1 = (1e+30, 'K'),
        T2 = (1e+30, 'K'),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is C2H3 + H <=> C2H4""",
)

entry(
    index = 222,
    label = "C2H4 <=> H2CC + H2",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(8e+12, 's^-1'), n=0.44, Ea=(88800, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (7e+50, 'cm^3/(mol*s)'),
            n = -9.31,
            Ea = (99900, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.735,
        T3 = (180, 'K'),
        T1 = (1035, 'K'),
        T2 = (5417, 'K'),
        efficiencies = {'O': 6, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C2H4 <=> H2CC + H2""",
)

entry(
    index = 223,
    label = "C2H4 + H <=> C2H3 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(240, 'cm^3/(mol*s)'), n=3.62, Ea=(11266, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4 + H <=> C2H3 + H2""",
)

entry(
    index = 224,
    label = "C2H4 + O <=> CH3 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.9e+17, 'cm^3/(mol*s)'),
        n = -1.717,
        Ea = (2891, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H4 + O <=> CH3 + HCO""",
)

entry(
    index = 225,
    label = "C2H4 + O <=> CH3 + H + CO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.5e+27, 'cm^3/(mol*s)'),
        n = -4.19,
        Ea = (22819, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H4 + O <=> CH3 + H + CO""",
)

entry(
    index = 226,
    label = "C2H4 + O <=> CH3CO + H",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8.7e+12, 'cm^3/(mol*s)'),
        n = -0.4843,
        Ea = (1957, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H4 + O <=> CH3CO + H""",
)

entry(
    index = 227,
    label = "C2H4 + O <=> CH2CHO + H",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (9.2e+09, 'cm^3/(mol*s)'),
        n = 0.9475,
        Ea = (1723, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H4 + O <=> CH2CHO + H""",
)

entry(
    index = 228,
    label = "C2H4 + O <=> CH2 + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.8e+06, 'cm^3/(mol*s)'),
        n = 1.991,
        Ea = (2858, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H4 + O <=> CH2 + CH2O""",
)

entry(
    index = 229,
    label = "C2H4 + O <=> CH2CO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.1e+17, 'cm^3/(mol*s)'),
        n = -1.831,
        Ea = (3177, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H4 + O <=> CH2CO + H2""",
)

entry(
    index = 230,
    label = "C2H4 + OH <=> C2H3 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.13, 'cm^3/(mol*s)'), n=4.2, Ea=(-860, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4 + OH <=> C2H3 + H2O""",
)

entry(
    index = 231,
    label = "C2H4 + OH <=> CH3 + CH2O",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.025, 0.1, 1, 10, 100], 'atm'),
        arrhenius = [
            Arrhenius(A=(5.4, 'cm^3/(mol*s)'), n=2.92, Ea=(-1733, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(32, 'cm^3/(mol*s)'), n=2.71, Ea=(-1172, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(560, 'cm^3/(mol*s)'), n=2.36, Ea=(-181, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(180000, 'cm^3/(mol*s)'), n=1.68, Ea=(2061, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (2.4e+09, 'cm^3/(mol*s)'),
                n = 0.56,
                Ea = (6007, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.8e+13, 'cm^3/(mol*s)'),
                n = -0.5,
                Ea = (11455, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H4 + OH <=> CH3 + CH2O""",
)

entry(
    index = 232,
    label = "C2H4 + OH <=> CH3CHO + H",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.025, 0.1, 1, 10, 100], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (2.4e-07, 'cm^3/(mol*s)'),
                n = 5.3,
                Ea = (-2051, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (8.7e-05, 'cm^3/(mol*s)'),
                n = 4.57,
                Ea = (-618, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(A=(0.4, 'cm^3/(mol*s)'), n=3.54, Ea=(1882, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(0.024, 'cm^3/(mol*s)'), n=3.91, Ea=(1723, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (8.3e+08, 'cm^3/(mol*s)'),
                n = 1.01,
                Ea = (10507, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (6.8e+09, 'cm^3/(mol*s)'),
                n = 0.81,
                Ea = (13867, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H4 + OH <=> CH3CHO + H""",
)

entry(
    index = 233,
    label = "C2H4 + OH <=> CH2CHOH + H",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.025, 0.1, 1, 10, 100], 'atm'),
        arrhenius = [
            Arrhenius(A=(10000, 'cm^3/(mol*s)'), n=2.6, Ea=(4121, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(11000, 'cm^3/(mol*s)'), n=2.6, Ea=(4129, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(15000, 'cm^3/(mol*s)'), n=2.56, Ea=(4238, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(320000, 'cm^3/(mol*s)'), n=2.19, Ea=(5256, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (1.9e+08, 'cm^3/(mol*s)'),
                n = 1.43,
                Ea = (7829, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (8.6e+10, 'cm^3/(mol*s)'),
                n = 0.75,
                Ea = (11491, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H4 + OH <=> CH2CHOH + H""",
)

entry(
    index = 234,
    label = "C2H4 + OH <=> CH2CH2OH",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiPDepArrhenius(
        arrhenius = [
            PDepArrhenius(
                pressures = ([0.01, 0.025, 0.1, 1, 10, 100], 'atm'),
                arrhenius = [
                    Arrhenius(
                        A = (1.4e+47, 'cm^3/(mol*s)'),
                        n = -11.64,
                        Ea = (11099, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(6e+37, 'cm^3/(mol*s)'), n=-9.76, Ea=(1995, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(6e+37, 'cm^3/(mol*s)'), n=-9.65, Ea=(2363, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(6e+37, 'cm^3/(mol*s)'), n=-8.14, Ea=(8043, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (6e+37, 'cm^3/(mol*s)'),
                        n = -7.77,
                        Ea = (10736, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (6e+37, 'cm^3/(mol*s)'),
                        n = -7.44,
                        Ea = (14269, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            PDepArrhenius(
                pressures = ([0.01, 0.025, 0.1, 1, 10, 100], 'atm'),
                arrhenius = [
                    Arrhenius(
                        A = (1.4e+47, 'cm^3/(mol*s)'),
                        n = -11.64,
                        Ea = (11099, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(5e+37, 'cm^3/(mol*s)'), n=-8.68, Ea=(5355, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (2.6e+35, 'cm^3/(mol*s)'),
                        n = -7.79,
                        Ea = (5017, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (7.3e+31, 'cm^3/(mol*s)'),
                        n = -6.91,
                        Ea = (2855, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(3e+26, 'cm^3/(mol*s)'), n=-4.87, Ea=(2297, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (2.8e+19, 'cm^3/(mol*s)'),
                        n = -2.41,
                        Ea = (1011, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H4 + OH <=> CH2CH2OH""",
)

entry(
    index = 235,
    label = "C2H4 + HO2 <=> CH2CH2OOH",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (3.1e+20, 'cm^3/(mol*s)'),
                n = -5.24,
                Ea = (11030, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.7e+20, 'cm^3/(mol*s)'),
                n = -5.14,
                Ea = (11430, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.6e+21, 'cm^3/(mol*s)'),
                n = -5.24,
                Ea = (12190, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.2e+23, 'cm^3/(mol*s)'),
                n = -5.53,
                Ea = (12990, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (6.7e+25, 'cm^3/(mol*s)'),
                n = -6.13,
                Ea = (14140, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.2e+27, 'cm^3/(mol*s)'),
                n = -6.37,
                Ea = (14980, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (8.7e+28, 'cm^3/(mol*s)'),
                n = -6.6,
                Ea = (16010, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (3.4e+29, 'cm^3/(mol*s)'),
                n = -6.55,
                Ea = (16800, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.7e+29, 'cm^3/(mol*s)'),
                n = -6.27,
                Ea = (17530, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.2e+28, 'cm^3/(mol*s)'),
                n = -5.71,
                Ea = (17940, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.5e+26, 'cm^3/(mol*s)'),
                n = -4.82,
                Ea = (18070, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.1e+23, 'cm^3/(mol*s)'),
                n = -3.77,
                Ea = (17820, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (4.4e+18, 'cm^3/(mol*s)'),
                n = -2.17,
                Ea = (16840, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H4 + HO2 <=> CH2CH2OOH""",
)

entry(
    index = 236,
    label = "C2H4 + HO2 <=> cC2H4O + OH",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100], 'atm'),
        arrhenius = [
            Arrhenius(A=(10000, 'cm^3/(mol*s)'), n=2.42, Ea=(12050, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(11000, 'cm^3/(mol*s)'), n=2.41, Ea=(12060, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(11000, 'cm^3/(mol*s)'), n=2.41, Ea=(12070, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(12000, 'cm^3/(mol*s)'), n=2.41, Ea=(12080, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(13000, 'cm^3/(mol*s)'), n=2.39, Ea=(12120, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(20000, 'cm^3/(mol*s)'), n=2.34, Ea=(12230, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(60000, 'cm^3/(mol*s)'), n=2.2, Ea=(12530, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (420000, 'cm^3/(mol*s)'),
                n = 1.96,
                Ea = (13090, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.3e+07, 'cm^3/(mol*s)'),
                n = 1.54,
                Ea = (14120, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (9.3e+08, 'cm^3/(mol*s)'),
                n = 1.02,
                Ea = (15470, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.1e+11, 'cm^3/(mol*s)'),
                n = 0.45,
                Ea = (17220, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.1e+12, 'cm^3/(mol*s)'),
                n = 0.11,
                Ea = (18750, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.2e+12, 'cm^3/(mol*s)'),
                n = 0.16,
                Ea = (19980, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H4 + HO2 <=> cC2H4O + OH""",
)

entry(
    index = 237,
    label = "C2H4 + O2 <=> C2H3 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.1e+13, 'cm^3/(mol*s)'), n=0, Ea=(60010, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4 + O2 <=> C2H3 + HO2""",
)

entry(
    index = 238,
    label = "C2H4 + CH3 <=> C2H3 + CH4",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(980, 'cm^3/(mol*s)'), n=2.947, Ea=(15148, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (8.1e-05, 'cm^3/(mol*s)'),
                n = 4.417,
                Ea = (8836, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H4 + CH3 <=> C2H3 + CH4""",
)

entry(
    index = 239,
    label = "C2H4 + CH2(S) <=> C2H3 + CH3",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiPDepArrhenius(
        arrhenius = [
            PDepArrhenius(
                pressures = ([0.01, 0.1, 1, 10, 100], 'atm'),
                arrhenius = [
                    Arrhenius(
                        A = (1.8e+19, 'cm^3/(mol*s)'),
                        n = -1.95,
                        Ea = (6787, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.7e+19, 'cm^3/(mol*s)'),
                        n = -1.8,
                        Ea = (4310, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (4.2e+24, 'cm^3/(mol*s)'),
                        n = -3.19,
                        Ea = (9759, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (7.9e+24, 'cm^3/(mol*s)'),
                        n = -3.08,
                        Ea = (13894, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (7.4e+29, 'cm^3/(mol*s)'),
                        n = -4.28,
                        Ea = (23849, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            PDepArrhenius(
                pressures = ([0.01, 0.1, 1, 10, 100], 'atm'),
                arrhenius = [
                    Arrhenius(
                        A = (4.3e+12, 'cm^3/(mol*s)'),
                        n = 0.19,
                        Ea = (-110, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(2.3e+11, 'cm^3/(mol*s)'), n=0.54, Ea=(48, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(4.9e+09, 'cm^3/(mol*s)'), n=1.02, Ea=(600, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (1.5e+08, 'cm^3/(mol*s)'),
                        n = 1.33,
                        Ea = (1228, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (8.1e+10, 'cm^3/(mol*s)'),
                        n = 0.55,
                        Ea = (5507, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H4 + CH2(S) <=> C2H3 + CH3""",
)

entry(
    index = 240,
    label = "C2H2 + H <=> C2H3",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(
            A = (1.7e+10, 'cm^3/(mol*s)'),
            n = 1.266,
            Ea = (2709, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        arrheniusLow = Arrhenius(
            A = (6.3e+31, 'cm^6/(mol^2*s)'),
            n = -4.664,
            Ea = (3780, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.7878,
        T3 = (-10212, 'K'),
        T1 = (1e-30, 'K'),
        efficiencies = {'[C-]#[O+]': 2, '[H][H]': 2, 'O=C=O': 3, 'O': 5},
    ),
    shortDesc = u"""The chemkin file reaction is C2H2 + H <=> C2H3""",
)

entry(
    index = 241,
    label = "C2H3 + H <=> C2H2 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3 + H <=> C2H2 + H2""",
)

entry(
    index = 242,
    label = "C2H3 + O <=> CH2CO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3 + O <=> CH2CO + H""",
)

entry(
    index = 243,
    label = "C2H3 + OH <=> C2H2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3 + OH <=> C2H2 + H2O""",
)

entry(
    index = 244,
    label = "C2H3 + HO2 <=> CH2CHO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3 + HO2 <=> CH2CHO + OH""",
)

entry(
    index = 245,
    label = "C2H3 + O2 <=> CH2CHOO",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiPDepArrhenius(
        arrhenius = [
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(
                        A = (1.6e+24, 'cm^3/(mol*s)'),
                        n = -5.45,
                        Ea = (9662, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (3.5e+56, 'cm^3/(mol*s)'),
                        n = -15.01,
                        Ea = (19160, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.3e+64, 'cm^3/(mol*s)'),
                        n = -16.97,
                        Ea = (21290, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (3.3e+61, 'cm^3/(mol*s)'),
                        n = -15.79,
                        Ea = (20150, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (7.3e+53, 'cm^3/(mol*s)'),
                        n = -13.11,
                        Ea = (17300, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (4.2e+48, 'cm^3/(mol*s)'),
                        n = -11.21,
                        Ea = (16000, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2.3e+43, 'cm^3/(mol*s)'),
                        n = -9.38,
                        Ea = (14810, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (3.4e+39, 'cm^3/(mol*s)'),
                        n = -8.04,
                        Ea = (14360, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(
                        A = (1.8e-09, 'cm^3/(mol*s)'),
                        n = 4.15,
                        Ea = (-4707, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2.4e+22, 'cm^3/(mol*s)'),
                        n = -4.52,
                        Ea = (2839, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(2e+26, 'cm^3/(mol*s)'), n=-5.43, Ea=(2725, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (6.1e+28, 'cm^3/(mol*s)'),
                        n = -5.89,
                        Ea = (3154, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2.1e+29, 'cm^3/(mol*s)'),
                        n = -5.8,
                        Ea = (3520, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (3.5e+28, 'cm^3/(mol*s)'),
                        n = -5.37,
                        Ea = (3636, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (3.3e+27, 'cm^3/(mol*s)'),
                        n = -4.95,
                        Ea = (3610, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(1e+27, 'cm^3/(mol*s)'), n=-4.72, Ea=(3680, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H3 + O2 <=> CH2CHOO""",
)

entry(
    index = 246,
    label = "C2H3 + O2 <=> CHCHO + OH",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiPDepArrhenius(
        arrhenius = [
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(
                        A = (3.9e+11, 'cm^3/(mol*s)'),
                        n = -0.11,
                        Ea = (2131, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(1.1e+09, 'cm^3/(mol*s)'), n=0.55, Ea=(46, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(8.5e+08, 'cm^3/(mol*s)'), n=0.56, Ea=(1, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(2.8e+14, 'cm^3/(mol*s)'), n=-1.83, Ea=(5, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (2.6e+20, 'cm^3/(mol*s)'),
                        n = -2.84,
                        Ea = (7530, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(9.2e+14, 'cm^3/(mol*s)'), n=-2.26, Ea=(0, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (6.1e+25, 'cm^3/(mol*s)'),
                        n = -4.21,
                        Ea = (13050, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.7e+30, 'cm^3/(mol*s)'),
                        n = -5.35,
                        Ea = (18430, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(A=(9.9e+11, 'cm^3/(mol*s)'), n=-0.66, Ea=(-1, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (6.9e+14, 'cm^3/(mol*s)'),
                        n = -1.16,
                        Ea = (4542, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2.8e+13, 'cm^3/(mol*s)'),
                        n = -0.72,
                        Ea = (3479, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(5e+11, 'cm^3/(mol*s)'), n=-0.14, Ea=(1995, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (2.4e+10, 'cm^3/(mol*s)'),
                        n = 0.23,
                        Ea = (1573, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.7e+14, 'cm^3/(mol*s)'),
                        n = -0.82,
                        Ea = (4450, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.4e+11, 'cm^3/(mol*s)'),
                        n = 0.05,
                        Ea = (3774, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (3.2e+11, 'cm^3/(mol*s)'),
                        n = -0.02,
                        Ea = (5338, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H3 + O2 <=> CHCHO + OH""",
)

entry(
    index = 247,
    label = "C2H3 + O2 <=> CH2CO + OH",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiPDepArrhenius(
        arrhenius = [
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(A=(870, 'cm^3/(mol*s)'), n=2.41, Ea=(6061, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(890, 'cm^3/(mol*s)'), n=2.41, Ea=(6078, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(940, 'cm^3/(mol*s)'), n=2.4, Ea=(6112, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1100, 'cm^3/(mol*s)'), n=2.39, Ea=(6180, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1100, 'cm^3/(mol*s)'), n=2.38, Ea=(6179, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1400, 'cm^3/(mol*s)'), n=2.36, Ea=(6074, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (2.5e+06, 'cm^3/(mol*s)'),
                        n = 1.42,
                        Ea = (8480, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.7e+10, 'cm^3/(mol*s)'),
                        n = 0.36,
                        Ea = (12010, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(A=(0.18, 'cm^3/(mol*s)'), n=3.12, Ea=(1331, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(0.21, 'cm^3/(mol*s)'), n=3.11, Ea=(1383, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(0.27, 'cm^3/(mol*s)'), n=3.08, Ea=(1496, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(0.53, 'cm^3/(mol*s)'), n=3.01, Ea=(1777, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.4, 'cm^3/(mol*s)'), n=2.9, Ea=(2225, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(0.42, 'cm^3/(mol*s)'), n=2.93, Ea=(2052, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (0.00012, 'cm^3/(mol*s)'),
                        n = 4.21,
                        Ea = (2043, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(0.0013, 'cm^3/(mol*s)'), n=3.97, Ea=(3414, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H3 + O2 <=> CH2CO + OH""",
)

entry(
    index = 248,
    label = "C2H3 + O2 <=> CH2CHO + O",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiPDepArrhenius(
        arrhenius = [
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(
                        A = (7.2e+20, 'cm^3/(mol*s)'),
                        n = -2.67,
                        Ea = (6742, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(7e+20, 'cm^3/(mol*s)'), n=-2.67, Ea=(6713, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(9e+20, 'cm^3/(mol*s)'), n=-2.7, Ea=(6724, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (6.5e+20, 'cm^3/(mol*s)'),
                        n = -2.65,
                        Ea = (6489, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (4.1e+20, 'cm^3/(mol*s)'),
                        n = -2.53,
                        Ea = (6406, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.6e+23, 'cm^3/(mol*s)'),
                        n = -3.22,
                        Ea = (8697, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2.9e+25, 'cm^3/(mol*s)'),
                        n = -3.77,
                        Ea = (11530, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (9.3e+25, 'cm^3/(mol*s)'),
                        n = -3.8,
                        Ea = (13910, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(
                        A = (1.2e+10, 'cm^3/(mol*s)'),
                        n = 0.62,
                        Ea = (-278, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.3e+10, 'cm^3/(mol*s)'),
                        n = 0.62,
                        Ea = (-248, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(1.5e+10, 'cm^3/(mol*s)'), n=0.6, Ea=(-163, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.8e+10, 'cm^3/(mol*s)'), n=0.58, Ea=(38, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(8.9e+09, 'cm^3/(mol*s)'), n=0.67, Ea=(248, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(6.7e+09, 'cm^3/(mol*s)'), n=0.72, Ea=(778, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (1.4e+09, 'cm^3/(mol*s)'),
                        n = 0.92,
                        Ea = (1219, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (7.1e+07, 'cm^3/(mol*s)'),
                        n = 1.28,
                        Ea = (1401, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H3 + O2 <=> CH2CHO + O""",
)

entry(
    index = 249,
    label = "C2H3 + O2 <=> C2H2 + HO2",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiPDepArrhenius(
        arrhenius = [
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(
                        A = (1.1e+07, 'cm^3/(mol*s)'),
                        n = 1.28,
                        Ea = (3322, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (7.8e+06, 'cm^3/(mol*s)'),
                        n = 1.33,
                        Ea = (3216, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.2e+07, 'cm^3/(mol*s)'),
                        n = 1.27,
                        Ea = (3311, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2.2e+07, 'cm^3/(mol*s)'),
                        n = 1.19,
                        Ea = (3367, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(1.1e+08, 'cm^3/(mol*s)'), n=1, Ea=(3695, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (1.3e+11, 'cm^3/(mol*s)'),
                        n = 0.12,
                        Ea = (5872, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.2e+09, 'cm^3/(mol*s)'),
                        n = 0.82,
                        Ea = (5617, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.1e+17, 'cm^3/(mol*s)'),
                        n = -1.45,
                        Ea = (12230, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(A=(48, 'cm^3/(mol*s)'), n=2.75, Ea=(-796, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(52, 'cm^3/(mol*s)'), n=2.73, Ea=(-768, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(56, 'cm^3/(mol*s)'), n=2.73, Ea=(-659, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(46, 'cm^3/(mol*s)'), n=2.76, Ea=(-493, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(3.8, 'cm^3/(mol*s)'), n=3.07, Ea=(-601, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(5.5, 'cm^3/(mol*s)'), n=3.07, Ea=(86, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(4.5e+08, 'cm^3/(mol*s)'), n=0, Ea=(955, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(20, 'cm^3/(mol*s)'), n=2.94, Ea=(1847, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H3 + O2 <=> C2H2 + HO2""",
)

entry(
    index = 250,
    label = "C2H3 + O2 <=> OCHCHO + H",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiPDepArrhenius(
        arrhenius = [
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(
                        A = (4.8e+14, 'cm^3/(mol*s)'),
                        n = -1.03,
                        Ea = (912, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(5e+14, 'cm^3/(mol*s)'), n=-1.04, Ea=(923, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (6.4e+14, 'cm^3/(mol*s)'),
                        n = -1.07,
                        Ea = (983, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (3.7e+15, 'cm^3/(mol*s)'),
                        n = -1.29,
                        Ea = (1441, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2.4e+18, 'cm^3/(mol*s)'),
                        n = -2.13,
                        Ea = (3234, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.3e+15, 'cm^3/(mol*s)'),
                        n = -1.09,
                        Ea = (2393, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (3.6e+33, 'cm^3/(mol*s)'),
                        n = -6.5,
                        Ea = (14910, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (3.3e+31, 'cm^3/(mol*s)'),
                        n = -5.76,
                        Ea = (16250, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(
                        A = (0.00028, 'cm^3/(mol*s)'),
                        n = 4.04,
                        Ea = (-7019, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (0.00035, 'cm^3/(mol*s)'),
                        n = 4.01,
                        Ea = (-6978, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (0.00097, 'cm^3/(mol*s)'),
                        n = 3.89,
                        Ea = (-6768, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(0.5, 'cm^3/(mol*s)'), n=3.15, Ea=(-5496, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (130000, 'cm^3/(mol*s)'),
                        n = 1.67,
                        Ea = (-2931, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (4.5e+15, 'cm^3/(mol*s)'),
                        n = -3.08,
                        Ea = (-4836, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(3.8e+10, 'cm^3/(mol*s)'), n=0.22, Ea=(941, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(2.8e+08, 'cm^3/(mol*s)'), n=0.83, Ea=(858, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H3 + O2 <=> OCHCHO + H""",
)

entry(
    index = 251,
    label = "C2H3 + O2 <=> CH2O + HCO",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiPDepArrhenius(
        arrhenius = [
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(
                        A = (1.4e+36, 'cm^3/(mol*s)'),
                        n = -7.6,
                        Ea = (12610, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2.6e+15, 'cm^3/(mol*s)'),
                        n = -1.28,
                        Ea = (513, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2.2e+36, 'cm^3/(mol*s)'),
                        n = -7.57,
                        Ea = (12490, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (3e+35, 'cm^3/(mol*s)'),
                        n = -7.32,
                        Ea = (11820, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.6e+36, 'cm^3/(mol*s)'),
                        n = -7.47,
                        Ea = (12460, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (5.8e+35, 'cm^3/(mol*s)'),
                        n = -7.2,
                        Ea = (13430, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (3.5e+20, 'cm^3/(mol*s)'),
                        n = -2.57,
                        Ea = (5578, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (3e+33, 'cm^3/(mol*s)'),
                        n = -6.28,
                        Ea = (16000, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(
                        A = (1.4e+36, 'cm^3/(mol*s)'),
                        n = -7.6,
                        Ea = (12610, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2.5e+15, 'cm^3/(mol*s)'),
                        n = -1.28,
                        Ea = (513, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (5.3e+15, 'cm^3/(mol*s)'),
                        n = -1.29,
                        Ea = (521, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (6.8e+15, 'cm^3/(mol*s)'),
                        n = -1.31,
                        Ea = (646, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.1e+16, 'cm^3/(mol*s)'),
                        n = -1.36,
                        Ea = (1066, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2.8e+15, 'cm^3/(mol*s)'),
                        n = -1.18,
                        Ea = (1429, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.1e+69, 'cm^3/(mol*s)'),
                        n = -19.23,
                        Ea = (14760, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(4.7e+10, 'cm^3/(mol*s)'), n=0.19, Ea=(831, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H3 + O2 <=> CH2O + HCO""",
)

entry(
    index = 252,
    label = "C2H3 + O2 <=> CH2O + H + CO",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiPDepArrhenius(
        arrhenius = [
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(
                        A = (6.5e+36, 'cm^3/(mol*s)'),
                        n = -7.6,
                        Ea = (12640, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (6.3e+36, 'cm^3/(mol*s)'),
                        n = -7.6,
                        Ea = (12610, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (5.1e+36, 'cm^3/(mol*s)'),
                        n = -7.57,
                        Ea = (12490, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (7.1e+35, 'cm^3/(mol*s)'),
                        n = -7.32,
                        Ea = (11820, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (3.7e+36, 'cm^3/(mol*s)'),
                        n = -7.47,
                        Ea = (12460, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.3e+36, 'cm^3/(mol*s)'),
                        n = -7.2,
                        Ea = (13430, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (8.3e+20, 'cm^3/(mol*s)'),
                        n = -2.57,
                        Ea = (5578, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (7.1e+33, 'cm^3/(mol*s)'),
                        n = -6.28,
                        Ea = (16000, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(
                        A = (1.2e+16, 'cm^3/(mol*s)'),
                        n = -1.28,
                        Ea = (515, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.2e+16, 'cm^3/(mol*s)'),
                        n = -1.28,
                        Ea = (513, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.3e+16, 'cm^3/(mol*s)'),
                        n = -1.29,
                        Ea = (521, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.6e+16, 'cm^3/(mol*s)'),
                        n = -1.31,
                        Ea = (646, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2.4e+16, 'cm^3/(mol*s)'),
                        n = -1.36,
                        Ea = (1066, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (6.6e+15, 'cm^3/(mol*s)'),
                        n = -1.18,
                        Ea = (1429, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2.7e+69, 'cm^3/(mol*s)'),
                        n = -19.23,
                        Ea = (14760, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(1.1e+11, 'cm^3/(mol*s)'), n=0.19, Ea=(831, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H3 + O2 <=> CH2O + H + CO""",
)

entry(
    index = 253,
    label = "C2H3 + O2 <=> CO + CH3O",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiPDepArrhenius(
        arrhenius = [
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(
                        A = (8.2e+18, 'cm^3/(mol*s)'),
                        n = -2.66,
                        Ea = (3201, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (4.1e+14, 'cm^3/(mol*s)'),
                        n = -1.32,
                        Ea = (886, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (4.3e+14, 'cm^3/(mol*s)'),
                        n = -1.33,
                        Ea = (901, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=-0.33, Ea=(-748, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.9e+12, 'cm^3/(mol*s)'), n=-3, Ea=(-8995, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.9e+24, 'cm^3/(mol*s)'), n=-5.63, Ea=(2, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (1.1e+18, 'cm^3/(mol*s)'),
                        n = -2.22,
                        Ea = (5178, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (5.8e+32, 'cm^3/(mol*s)'),
                        n = -6.45,
                        Ea = (16810, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(
                        A = (1.3e+09, 'cm^3/(mol*s)'),
                        n = 0.18,
                        Ea = (-1717, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (6e+11, 'cm^3/(mol*s)'),
                        n = -2.93,
                        Ea = (-9564, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2.9e+11, 'cm^3/(mol*s)'),
                        n = -2.93,
                        Ea = (-10120, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (5.8e+21, 'cm^3/(mol*s)'),
                        n = -3.54,
                        Ea = (4772, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(5e+15, 'cm^3/(mol*s)'), n=-1.62, Ea=(1849, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (9.3e+16, 'cm^3/(mol*s)'),
                        n = -1.96,
                        Ea = (3324, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1e+72, 'cm^3/(mol*s)'),
                        n = -20.69,
                        Ea = (15860, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.1e+09, 'cm^3/(mol*s)'),
                        n = 0.31,
                        Ea = (1024, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H3 + O2 <=> CO + CH3O""",
)

entry(
    index = 254,
    label = "C2H3 + O2 <=> CO2 + CH3",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiPDepArrhenius(
        arrhenius = [
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(
                        A = (2.4e+35, 'cm^3/(mol*s)'),
                        n = -7.76,
                        Ea = (12630, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.7e+35, 'cm^3/(mol*s)'),
                        n = -7.72,
                        Ea = (12520, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (4.5e+34, 'cm^3/(mol*s)'),
                        n = -7.55,
                        Ea = (12140, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (7.3e+31, 'cm^3/(mol*s)'),
                        n = -6.7,
                        Ea = (10440, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (3.6e+35, 'cm^3/(mol*s)'),
                        n = -7.75,
                        Ea = (12830, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2.1e+35, 'cm^3/(mol*s)'),
                        n = -7.53,
                        Ea = (14050, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (3.8e+18, 'cm^3/(mol*s)'),
                        n = -2.44,
                        Ea = (5408, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.2e+32, 'cm^3/(mol*s)'),
                        n = -6.32,
                        Ea = (16190, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(
                        A = (6.3e+13, 'cm^3/(mol*s)'),
                        n = -1.16,
                        Ea = (406, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (6.2e+13, 'cm^3/(mol*s)'),
                        n = -1.16,
                        Ea = (401, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (6.1e+13, 'cm^3/(mol*s)'),
                        n = -1.16,
                        Ea = (397, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (5.3e+13, 'cm^3/(mol*s)'),
                        n = -1.14,
                        Ea = (447, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.5e+14, 'cm^3/(mol*s)'),
                        n = -1.26,
                        Ea = (988, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=-1.11, Ea=(1409, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (1.4e+70, 'cm^3/(mol*s)'),
                        n = -20.11,
                        Ea = (15430, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(9.2e+08, 'cm^3/(mol*s)'), n=0.25, Ea=(855, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H3 + O2 <=> CO2 + CH3""",
)

entry(
    index = 255,
    label = "C2H3 + CH2O <=> C2H4 + HCO",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiPDepArrhenius(
        arrhenius = [
            PDepArrhenius(
                pressures = ([0.001, 0.01, 0.1, 1, 10, 100, 1000], 'atm'),
                arrhenius = [
                    Arrhenius(
                        A = (1.1e+07, 'cm^3/(mol*s)'),
                        n = 1.09,
                        Ea = (1807, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2.5e+07, 'cm^3/(mol*s)'),
                        n = 0.993,
                        Ea = (1995, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2.5e+08, 'cm^3/(mol*s)'),
                        n = 0.704,
                        Ea = (2596, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.4e+10, 'cm^3/(mol*s)'),
                        n = 0.209,
                        Ea = (3934, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (3.5e+13, 'cm^3/(mol*s)'),
                        n = -0.726,
                        Ea = (6944, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (3.3e+14, 'cm^3/(mol*s)'),
                        n = -0.866,
                        Ea = (10966, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(17, 'cm^3/(mol*s)'), n=3.17, Ea=(9400, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            PDepArrhenius(
                pressures = ([0.001, 0.01, 0.1, 1, 10, 100, 1000], 'atm'),
                arrhenius = [
                    Arrhenius(
                        A = (-2.3e+16, 'cm^3/(mol*s)'),
                        n = -1.269,
                        Ea = (20617, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (-5.2e+16, 'cm^3/(mol*s)'),
                        n = -1.366,
                        Ea = (20805, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (-1.5e+18, 'cm^3/(mol*s)'),
                        n = -1.769,
                        Ea = (22524, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (-8.5e+19, 'cm^3/(mol*s)'),
                        n = -2.264,
                        Ea = (23862, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (-4.4e+23, 'cm^3/(mol*s)'),
                        n = -3.278,
                        Ea = (27795, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (-4.2e+24, 'cm^3/(mol*s)'),
                        n = -3.418,
                        Ea = (31817, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (-2.1e+11, 'cm^3/(mol*s)'),
                        n = 0.618,
                        Ea = (30251, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H3 + CH2O <=> C2H4 + HCO""",
)

entry(
    index = 256,
    label = "C2H3 + CH2O <=> C2H4 + H + CO",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.001, 0.01, 0.1, 1, 10, 100, 1000], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (2.3e+16, 'cm^3/(mol*s)'),
                n = -1.269,
                Ea = (20617, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (5.2e+16, 'cm^3/(mol*s)'),
                n = -1.366,
                Ea = (20805, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.5e+18, 'cm^3/(mol*s)'),
                n = -1.769,
                Ea = (22524, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (8.5e+19, 'cm^3/(mol*s)'),
                n = -2.264,
                Ea = (23862, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (4.4e+23, 'cm^3/(mol*s)'),
                n = -3.278,
                Ea = (27795, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (4.2e+24, 'cm^3/(mol*s)'),
                n = -3.418,
                Ea = (31817, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.1e+11, 'cm^3/(mol*s)'),
                n = 0.618,
                Ea = (30251, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H3 + CH2O <=> C2H4 + H + CO""",
)

entry(
    index = 257,
    label = "C2H3 + HCO <=> C2H4 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3 + HCO <=> C2H4 + CO""",
)

entry(
    index = 258,
    label = "C2H3 + CH3 <=> C2H2 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+12, 'cm^3/(mol*s)'), n=0, Ea=(-765, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3 + CH3 <=> C2H2 + CH4""",
)

entry(
    index = 259,
    label = "C2H3 + CH <=> CH2 + C2H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3 + CH <=> CH2 + C2H2""",
)

entry(
    index = 260,
    label = "C2H3 + CH3OH <=> C2H4 + CH3O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+06, 'cm^3/(mol*s)'), n=1.51, Ea=(26630, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3 + CH3OH <=> C2H4 + CH3O""",
)

entry(
    index = 261,
    label = "C2H3 + CH3OH <=> C2H4 + CH2OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.018, 'cm^3/(mol*s)'), n=4.02, Ea=(23370, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3 + CH3OH <=> C2H4 + CH2OH""",
)

entry(
    index = 262,
    label = "C2H3 + C2H3 <=> C2H4 + C2H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3 + C2H3 <=> C2H4 + C2H2""",
)

entry(
    index = 263,
    label = "C2H3 + C2H <=> C2H2 + C2H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3 + C2H <=> C2H2 + C2H2""",
)

entry(
    index = 264,
    label = "C2H2 <=> C2H + H",
    degeneracy = 1,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(
            A = (9.1e+30, 'cm^3/(mol*s)'),
            n = -3.7,
            Ea = (127138, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        efficiencies = {'[C-]#[O+]': 2, '[H][H]': 2, 'O=C=O': 3, 'O': 5},
    ),
    shortDesc = u"""The chemkin file reaction is C2H2 <=> C2H + H""",
)

entry(
    index = 265,
    label = "C2H2 <=> H2CC",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(18000, 's^-1'), n=3.51, Ea=(43300, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (2.5e+15, 'cm^3/(mol*s)'),
            n = -0.64,
            Ea = (49700, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.5,
        T3 = (1e-30, 'K'),
        T1 = (1e+30, 'K'),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is C2H2 <=> H2CC""",
)

entry(
    index = 266,
    label = "C2H + H2 <=> C2H2 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(410000, 'cm^3/(mol*s)'), n=2.39, Ea=(864, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H + H2 <=> C2H2 + H""",
)

entry(
    index = 267,
    label = "C2H2 + O <=> HCCO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.4e+07, 'cm^3/(mol*s)'), n=2, Ea=(1900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H2 + O <=> HCCO + H""",
)

entry(
    index = 268,
    label = "C2H2 + O <=> CH2 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.1e+06, 'cm^3/(mol*s)'), n=2, Ea=(1900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H2 + O <=> CH2 + CO""",
)

entry(
    index = 269,
    label = "C2H2 + OH <=> CH3 + CO",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.025, 0.1, 1, 10, 100], 'atm'),
        arrhenius = [
            Arrhenius(A=(480000, 'cm^3/(mol*s)'), n=1.68, Ea=(-330, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(4.4e+06, 'cm^3/(mol*s)'), n=1.4, Ea=(227, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (7.7e+07, 'cm^3/(mol*s)'),
                n = 1.05,
                Ea = (1115, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.3e+09, 'cm^3/(mol*s)'),
                n = 0.73,
                Ea = (2579, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (4.3e+08, 'cm^3/(mol*s)'),
                n = 0.92,
                Ea = (3736, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(A=(830000, 'cm^3/(mol*s)'), n=1.77, Ea=(4697, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H2 + OH <=> CH3 + CO""",
)

entry(
    index = 270,
    label = "C2H2 + OH <=> HCCOH + H",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.025, 0.1, 1, 10, 100], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (280000, 'cm^3/(mol*s)'),
                n = 2.28,
                Ea = (12419, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (750000, 'cm^3/(mol*s)'),
                n = 2.16,
                Ea = (12547, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.8e+06, 'cm^3/(mol*s)'),
                n = 2.04,
                Ea = (12669, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(A=(2.4e+06, 'cm^3/(mol*s)'), n=2, Ea=(12713, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (3.2e+06, 'cm^3/(mol*s)'),
                n = 1.97,
                Ea = (12810, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (7.4e+06, 'cm^3/(mol*s)'),
                n = 1.89,
                Ea = (13603, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H2 + OH <=> HCCOH + H""",
)

entry(
    index = 271,
    label = "C2H2 + OH <=> CHCHOH",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiPDepArrhenius(
        arrhenius = [
            PDepArrhenius(
                pressures = ([0.01, 0.025, 0.1, 1, 10, 100], 'atm'),
                arrhenius = [
                    Arrhenius(
                        A = (2.9e+64, 'cm^3/(mol*s)'),
                        n = -18.57,
                        Ea = (10009, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (4.7e+59, 'cm^3/(mol*s)'),
                        n = -16.87,
                        Ea = (9087, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.2e+28, 'cm^3/(mol*s)'),
                        n = -5.56,
                        Ea = (3724, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.9e+44, 'cm^3/(mol*s)'),
                        n = -11.38,
                        Ea = (6299, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.5e+24, 'cm^3/(mol*s)'),
                        n = -4.06,
                        Ea = (3261, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (6.2e+20, 'cm^3/(mol*s)'),
                        n = -2.8,
                        Ea = (2831, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            PDepArrhenius(
                pressures = ([0.01, 0.025, 0.1, 1, 10, 100], 'atm'),
                arrhenius = [
                    Arrhenius(
                        A = (2.6e+33, 'cm^3/(mol*s)'),
                        n = -7.36,
                        Ea = (6392, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (4.4e+32, 'cm^3/(mol*s)'),
                        n = -7.02,
                        Ea = (5933, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (6.4e+42, 'cm^3/(mol*s)'),
                        n = -9.96,
                        Ea = (11737, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (3.5e+31, 'cm^3/(mol*s)'),
                        n = -6.2,
                        Ea = (6635, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (4.5e+31, 'cm^3/(mol*s)'),
                        n = -5.92,
                        Ea = (8761, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.6e+29, 'cm^3/(mol*s)'),
                        n = -4.91,
                        Ea = (9734, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H2 + OH <=> CHCHOH""",
)

entry(
    index = 272,
    label = "C2H2 + OH <=> CH2CO + H",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.025, 0.1, 1, 10, 100], 'atm'),
        arrhenius = [
            Arrhenius(A=(1700, 'cm^3/(mol*s)'), n=2.56, Ea=(-844, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(15000, 'cm^3/(mol*s)'), n=2.28, Ea=(-292, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(300000, 'cm^3/(mol*s)'), n=1.92, Ea=(598, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (7.5e+06, 'cm^3/(mol*s)'),
                n = 1.55,
                Ea = (2106, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (5.1e+06, 'cm^3/(mol*s)'),
                n = 1.65,
                Ea = (3400, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(A=(15000, 'cm^3/(mol*s)'), n=2.45, Ea=(4477, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H2 + OH <=> CH2CO + H""",
)

entry(
    index = 273,
    label = "C2H2 + OH <=> C2H + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.6e+06, 'cm^3/(mol*s)'),
        n = 2.14,
        Ea = (17060, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H2 + OH <=> C2H + H2O""",
)

entry(
    index = 274,
    label = "C2H2 + HO2 <=> CH2CHOO",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiPDepArrhenius(
        arrhenius = [
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(A=(5e+06, 'cm^3/(mol*s)'), n=-1.02, Ea=(9152, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (6e+17, 'cm^3/(mol*s)'),
                        n = -3.82,
                        Ea = (10790, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2.5e+48, 'cm^3/(mol*s)'),
                        n = -12.82,
                        Ea = (25220, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (4.1e+50, 'cm^3/(mol*s)'),
                        n = -13.07,
                        Ea = (27220, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (9.1e+46, 'cm^3/(mol*s)'),
                        n = -11.57,
                        Ea = (26880, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (4.6e+43, 'cm^3/(mol*s)'),
                        n = -10.24,
                        Ea = (26930, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (5.6e+38, 'cm^3/(mol*s)'),
                        n = -8.49,
                        Ea = (26210, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2.5e+35, 'cm^3/(mol*s)'),
                        n = -7.26,
                        Ea = (26390, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(
                        A = (1.9e+26, 'cm^3/(mol*s)'),
                        n = -8.34,
                        Ea = (9249, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (5.3e+129, 'cm^3/(mol*s)'),
                        n = -41.74,
                        Ea = (35930, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2e+18, 'cm^3/(mol*s)'),
                        n = -3.67,
                        Ea = (10480, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (4.9e+21, 'cm^3/(mol*s)'),
                        n = -4.37,
                        Ea = (12220, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.9e+22, 'cm^3/(mol*s)'),
                        n = -4.28,
                        Ea = (13080, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2.1e+21, 'cm^3/(mol*s)'),
                        n = -3.78,
                        Ea = (13380, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.4e+20, 'cm^3/(mol*s)'),
                        n = -3.3,
                        Ea = (13410, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.4e+19, 'cm^3/(mol*s)'),
                        n = -2.91,
                        Ea = (13420, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H2 + HO2 <=> CH2CHOO""",
)

entry(
    index = 275,
    label = "C2H2 + HO2 <=> CHCHO + OH",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiPDepArrhenius(
        arrhenius = [
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(
                        A = (5.5e+09, 'cm^3/(mol*s)'),
                        n = 0.91,
                        Ea = (18500, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (5.9e+09, 'cm^3/(mol*s)'),
                        n = 0.9,
                        Ea = (18550, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (6.8e+09, 'cm^3/(mol*s)'),
                        n = 0.88,
                        Ea = (18640, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.6e+10, 'cm^3/(mol*s)'),
                        n = 0.77,
                        Ea = (19040, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (3.5e+09, 'cm^3/(mol*s)'),
                        n = 0.99,
                        Ea = (18810, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (5.4e+10, 'cm^3/(mol*s)'),
                        n = 0.61,
                        Ea = (20740, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (3.7e+08, 'cm^3/(mol*s)'),
                        n = 1.23,
                        Ea = (15960, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.5e+11, 'cm^3/(mol*s)'),
                        n = 0.48,
                        Ea = (17730, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(
                        A = (2.4e+07, 'cm^3/(mol*s)'),
                        n = 1.54,
                        Ea = (14690, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2.5e+07, 'cm^3/(mol*s)'),
                        n = 1.54,
                        Ea = (14700, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2.6e+07, 'cm^3/(mol*s)'),
                        n = 1.54,
                        Ea = (14730, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2.5e+07, 'cm^3/(mol*s)'),
                        n = 1.56,
                        Ea = (14790, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.5e+08, 'cm^3/(mol*s)'),
                        n = 1.32,
                        Ea = (15090, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.6e+08, 'cm^3/(mol*s)'),
                        n = 1.36,
                        Ea = (15420, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.7e+07, 'cm^3/(mol*s)'),
                        n = 1.59,
                        Ea = (15910, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (7.2e+06, 'cm^3/(mol*s)'),
                        n = 1.73,
                        Ea = (16020, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H2 + HO2 <=> CHCHO + OH""",
)

entry(
    index = 276,
    label = "C2H2 + HO2 <=> CH2CO + OH",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiPDepArrhenius(
        arrhenius = [
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(
                        A = (6.3e-07, 'cm^3/(mol*s)'),
                        n = 4.75,
                        Ea = (15530, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (6.7e-07, 'cm^3/(mol*s)'),
                        n = 4.74,
                        Ea = (15550, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (4.2e-07, 'cm^3/(mol*s)'),
                        n = 4.81,
                        Ea = (15410, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (5.3e-07, 'cm^3/(mol*s)'),
                        n = 4.78,
                        Ea = (15460, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(1e-06, 'cm^3/(mol*s)'), n=4.69, Ea=(15640, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (4.7e-05, 'cm^3/(mol*s)'),
                        n = 4.22,
                        Ea = (16780, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(0.9, 'cm^3/(mol*s)'), n=2.97, Ea=(19730, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(3600, 'cm^3/(mol*s)'), n=1.97, Ea=(23010, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(
                        A = (1.3e-14, 'cm^3/(mol*s)'),
                        n = 6.58,
                        Ea = (10270, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.3e-14, 'cm^3/(mol*s)'),
                        n = 6.59,
                        Ea = (10330, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(4e-14, 'cm^3/(mol*s)'), n=6.36, Ea=(10270, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (3.3e-15, 'cm^3/(mol*s)'),
                        n = 6.7,
                        Ea = (10090, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(8.7e-21, 'cm^3/(mol*s)'), n=8.3, Ea=(8107, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (8.4e-22, 'cm^3/(mol*s)'),
                        n = 8.76,
                        Ea = (8804, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (6.9e-14, 'cm^3/(mol*s)'),
                        n = 6.67,
                        Ea = (13130, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (6.6e-12, 'cm^3/(mol*s)'),
                        n = 6.15,
                        Ea = (14730, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H2 + HO2 <=> CH2CO + OH""",
)

entry(
    index = 277,
    label = "C2H2 + HO2 <=> CH2CHO + O",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiPDepArrhenius(
        arrhenius = [
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(A=(5.5, 'cm^3/(mol*s)'), n=1.19, Ea=(12880, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (1.1e+08, 'cm^3/(mol*s)'),
                        n = 0.77,
                        Ea = (13600, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.2e+07, 'cm^3/(mol*s)'),
                        n = 1.09,
                        Ea = (13050, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(3e+07, 'cm^3/(mol*s)'), n=0.98, Ea=(13310, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (2e+74, 'cm^3/(mol*s)'),
                        n = -16.33,
                        Ea = (109200, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (7.5e+14, 'cm^3/(mol*s)'),
                        n = -1.17,
                        Ea = (18350, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (8.6e+18, 'cm^3/(mol*s)'),
                        n = -2.27,
                        Ea = (22230, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (5.8e+18, 'cm^3/(mol*s)'),
                        n = -2.09,
                        Ea = (24350, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(
                        A = (0.00029, 'cm^3/(mol*s)'),
                        n = 4.16,
                        Ea = (7736, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(6100, 'cm^3/(mol*s)'), n=3.81, Ea=(8394, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (0.00054, 'cm^3/(mol*s)'),
                        n = 4.09,
                        Ea = (8044, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (0.00025, 'cm^3/(mol*s)'),
                        n = 4.19,
                        Ea = (8203, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(66000, 'cm^3/(mol*s)'), n=1.85, Ea=(12360, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(0.29, 'cm^3/(mol*s)'), n=3.38, Ea=(10590, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(2, 'cm^3/(mol*s)'), n=3.17, Ea=(11740, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(0.11, 'cm^3/(mol*s)'), n=3.52, Ea=(11980, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H2 + HO2 <=> CH2CHO + O""",
)

entry(
    index = 278,
    label = "C2H2 + HO2 <=> OCHCHO + H",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiPDepArrhenius(
        arrhenius = [
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(
                        A = (8.5e+07, 'cm^3/(mol*s)'),
                        n = 0.48,
                        Ea = (11720, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (7.4e+07, 'cm^3/(mol*s)'),
                        n = 0.5,
                        Ea = (11690, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (7.9e+07, 'cm^3/(mol*s)'),
                        n = 0.49,
                        Ea = (11700, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2.2e+09, 'cm^3/(mol*s)'),
                        n = 0.06,
                        Ea = (12470, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (7e+49, 'cm^3/(mol*s)'),
                        n = -10.18,
                        Ea = (77110, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (4.1e+16, 'cm^3/(mol*s)'),
                        n = -2.03,
                        Ea = (17630, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (9.4e+16, 'cm^3/(mol*s)'),
                        n = -2.03,
                        Ea = (19590, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (5.9e+21, 'cm^3/(mol*s)'),
                        n = -3.32,
                        Ea = (25030, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(
                        A = (2.4e-06, 'cm^3/(mol*s)'),
                        n = 4.43,
                        Ea = (5578, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(2e-06, 'cm^3/(mol*s)'), n=4.45, Ea=(5564, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (1.8e-06, 'cm^3/(mol*s)'),
                        n = 4.46,
                        Ea = (5654, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2.2e-05, 'cm^3/(mol*s)'),
                        n = 4.17,
                        Ea = (6416, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (770000, 'cm^3/(mol*s)'),
                        n = 1.18,
                        Ea = (11340, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(0.02, 'cm^3/(mol*s)'), n=3.38, Ea=(8696, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(0.0061, 'cm^3/(mol*s)'), n=3.53, Ea=(9217, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(0.068, 'cm^3/(mol*s)'), n=3.27, Ea=(10760, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H2 + HO2 <=> OCHCHO + H""",
)

entry(
    index = 279,
    label = "C2H2 + HO2 <=> CH2O + HCO",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiPDepArrhenius(
        arrhenius = [
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(
                        A = (3.9e+13, 'cm^3/(mol*s)'),
                        n = -1.17,
                        Ea = (13750, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(4.3, 'cm^3/(mol*s)'), n=2.64, Ea=(7253, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (2.6e-06, 'cm^3/(mol*s)'),
                        n = 4.34,
                        Ea = (4525, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (3.3e+102, 'cm^3/(mol*s)'),
                        n = -24.18,
                        Ea = (138600, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (5.2e+15, 'cm^3/(mol*s)'),
                        n = -1.75,
                        Ea = (15180, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (7.3e+35, 'cm^3/(mol*s)'),
                        n = -7.77,
                        Ea = (26970, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.8e+28, 'cm^3/(mol*s)'),
                        n = -5.3,
                        Ea = (25130, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2.5e+16, 'cm^3/(mol*s)'),
                        n = -1.7,
                        Ea = (20030, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(A=(8.4, 'cm^3/(mol*s)'), n=2.56, Ea=(7382, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (1.6e+13, 'cm^3/(mol*s)'),
                        n = -1.05,
                        Ea = (13520, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(6.9e+09, 'cm^3/(mol*s)'), n=0, Ea=(11720, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (8.1e+07, 'cm^3/(mol*s)'),
                        n = 0.6,
                        Ea = (10850, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(3.5, 'cm^3/(mol*s)'), n=2.69, Ea=(8025, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (9.8e+06, 'cm^3/(mol*s)'),
                        n = 0.91,
                        Ea = (11710, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(18000, 'cm^3/(mol*s)'), n=1.7, Ea=(11250, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (4.3e-06, 'cm^3/(mol*s)'),
                        n = 4.31,
                        Ea = (6829, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H2 + HO2 <=> CH2O + HCO""",
)

entry(
    index = 280,
    label = "C2H2 + HO2 <=> CH2O + H + CO",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiPDepArrhenius(
        arrhenius = [
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(
                        A = (9.1e+13, 'cm^3/(mol*s)'),
                        n = -1.17,
                        Ea = (13750, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(9.9, 'cm^3/(mol*s)'), n=2.64, Ea=(7253, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (6.1e-06, 'cm^3/(mol*s)'),
                        n = 4.34,
                        Ea = (4525, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (7.8e+102, 'cm^3/(mol*s)'),
                        n = -24.18,
                        Ea = (138600, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.2e+16, 'cm^3/(mol*s)'),
                        n = -1.75,
                        Ea = (15180, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.7e+36, 'cm^3/(mol*s)'),
                        n = -7.77,
                        Ea = (26970, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (4.1e+28, 'cm^3/(mol*s)'),
                        n = -5.3,
                        Ea = (25130, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (5.8e+16, 'cm^3/(mol*s)'),
                        n = -1.7,
                        Ea = (20030, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(A=(20, 'cm^3/(mol*s)'), n=2.56, Ea=(7382, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (3.6e+13, 'cm^3/(mol*s)'),
                        n = -1.05,
                        Ea = (13520, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(1.6e+10, 'cm^3/(mol*s)'), n=0, Ea=(11720, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (1.9e+08, 'cm^3/(mol*s)'),
                        n = 0.6,
                        Ea = (10850, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(8.3, 'cm^3/(mol*s)'), n=2.69, Ea=(8025, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (2.3e+07, 'cm^3/(mol*s)'),
                        n = 0.91,
                        Ea = (11710, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(42000, 'cm^3/(mol*s)'), n=1.7, Ea=(11250, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1e-05, 'cm^3/(mol*s)'), n=4.31, Ea=(6829, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H2 + HO2 <=> CH2O + H + CO""",
)

entry(
    index = 281,
    label = "C2H2 + HO2 <=> CO + CH3O",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiPDepArrhenius(
        arrhenius = [
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(A=(3.5e+11, 'cm^3/(mol*s)'), n=0, Ea=(49510, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (2.8e+08, 'cm^3/(mol*s)'),
                        n = 0.01,
                        Ea = (11920, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (8.1e+07, 'cm^3/(mol*s)'),
                        n = 0.18,
                        Ea = (11650, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (8.9e+69, 'cm^3/(mol*s)'),
                        n = -15.85,
                        Ea = (102500, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (5.7e+12, 'cm^3/(mol*s)'),
                        n = -1.25,
                        Ea = (14570, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (3.3e+23, 'cm^3/(mol*s)'),
                        n = -4.45,
                        Ea = (21210, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2.4e+22, 'cm^3/(mol*s)'),
                        n = -3.96,
                        Ea = (22650, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.2e+18, 'cm^3/(mol*s)'),
                        n = -2.57,
                        Ea = (22360, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(A=(29000, 'cm^3/(mol*s)'), n=1.23, Ea=(9903, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (9.7e-07, 'cm^3/(mol*s)'),
                        n = 4.15,
                        Ea = (5173, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.8e-08, 'cm^3/(mol*s)'),
                        n = 4.62,
                        Ea = (4517, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (540000, 'cm^3/(mol*s)'),
                        n = 0.86,
                        Ea = (10700, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (0.00054, 'cm^3/(mol*s)'),
                        n = 3.42,
                        Ea = (7218, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(290, 'cm^3/(mol*s)'), n=1.84, Ea=(10460, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(8.1, 'cm^3/(mol*s)'), n=2.3, Ea=(10560, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (0.00069, 'cm^3/(mol*s)'),
                        n = 3.42,
                        Ea = (9329, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H2 + HO2 <=> CO + CH3O""",
)

entry(
    index = 282,
    label = "C2H2 + HO2 <=> CO2 + CH3",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiPDepArrhenius(
        arrhenius = [
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(
                        A = (1.2e-07, 'cm^3/(mol*s)'),
                        n = 4.31,
                        Ea = (4614, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.1e-07, 'cm^3/(mol*s)'),
                        n = 4.32,
                        Ea = (4622, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.8e+142, 'cm^3/(mol*s)'),
                        n = -35.04,
                        Ea = (188700, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (4e+84, 'cm^3/(mol*s)'),
                        n = -19.8,
                        Ea = (119800, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=-1.6, Ea=(14980, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (8.6e+28, 'cm^3/(mol*s)'),
                        n = -6.15,
                        Ea = (24030, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.3e+27, 'cm^3/(mol*s)'),
                        n = -5.42,
                        Ea = (25380, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.7e+15, 'cm^3/(mol*s)'),
                        n = -1.8,
                        Ea = (20370, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(A=(2e+08, 'cm^3/(mol*s)'), n=0, Ea=(11790, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(2e+08, 'cm^3/(mol*s)'), n=0, Ea=(11780, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (160000, 'cm^3/(mol*s)'),
                        n = 0.95,
                        Ea = (10200, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.4e+06, 'cm^3/(mol*s)'),
                        n = 0.68,
                        Ea = (10810, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(0.0093, 'cm^3/(mol*s)'), n=3, Ea=(7659, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(19000, 'cm^3/(mol*s)'), n=1.26, Ea=(11230, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(290, 'cm^3/(mol*s)'), n=1.79, Ea=(11240, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (3.9e-07, 'cm^3/(mol*s)'),
                        n = 4.21,
                        Ea = (7314, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H2 + HO2 <=> CO2 + CH3""",
)

entry(
    index = 283,
    label = "C2H2 + O2 <=> HCO + HCO",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(6.1e+12, 'cm^3/(mol*s)'), n=0, Ea=(53250, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (1.7e+07, 'cm^3/(mol*s)'),
                n = 1.67,
                Ea = (70960, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H2 + O2 <=> HCO + HCO""",
)

entry(
    index = 284,
    label = "C2H2 + O2 <=> HCO + H + CO",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(
                A = (6.7e+33, 'cm^3/(mol*s)'),
                n = -5.633,
                Ea = (82336, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.1e+26, 'cm^3/(mol*s)'),
                n = -3.525,
                Ea = (73959, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H2 + O2 <=> HCO + H + CO""",
)

entry(
    index = 285,
    label = "C2H2 + O2 <=> H + CO + H + CO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.218e+32, 'cm^3/(mol*s)'),
        n = -4.869,
        Ea = (93010.6, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H2 + O2 <=> H + CO + H + CO""",
)

entry(
    index = 286,
    label = "C2H2 + CH2(S) <=> C2H2 + CH2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H2 + CH2(S) <=> C2H2 + CH2""",
)

entry(
    index = 287,
    label = "H2CC + H <=> C2H2 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2CC + H <=> C2H2 + H""",
)

entry(
    index = 288,
    label = "H2CC + OH <=> CH2CO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2CC + OH <=> CH2CO + H""",
)

entry(
    index = 289,
    label = "H2CC + O2 <=> CH2 + CO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2CC + O2 <=> CH2 + CO2""",
)

entry(
    index = 290,
    label = "C2 + H2 <=> C2H + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(400000, 'cm^3/(mol*s)'), n=2.4, Ea=(1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2 + H2 <=> C2H + H""",
)

entry(
    index = 291,
    label = "C2H + O <=> CH + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H + O <=> CH + CO""",
)

entry(
    index = 292,
    label = "C2H + OH <=> HCCO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H + OH <=> HCCO + H""",
)

entry(
    index = 293,
    label = "C2H + OH <=> C2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+07, 'cm^3/(mol*s)'), n=2, Ea=(8000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H + OH <=> C2 + H2O""",
)

entry(
    index = 294,
    label = "C2H + O2 <=> CO + CO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.7e+13, 'cm^3/(mol*s)'), n=-0.16, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H + O2 <=> CO + CO + H""",
)

entry(
    index = 295,
    label = "C2H + CH4 <=> CH3 + C2H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.2e+12, 'cm^3/(mol*s)'), n=0, Ea=(976, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H + CH4 <=> CH3 + C2H2""",
)

entry(
    index = 296,
    label = "C2 <=> C + C",
    degeneracy = 1,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(A=(1.5e+16, 'cm^3/(mol*s)'), n=0, Ea=(142300, 'cal/mol'), T0=(1, 'K')),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is C2 <=> C + C""",
)

entry(
    index = 297,
    label = "C2 + O <=> C + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2 + O <=> C + CO""",
)

entry(
    index = 298,
    label = "C2 + OH <=> C2O + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2 + OH <=> C2O + H""",
)

entry(
    index = 299,
    label = "C2 + O2 <=> CO + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+12, 'cm^3/(mol*s)'), n=0, Ea=(980, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2 + O2 <=> CO + CO""",
)

entry(
    index = 300,
    label = "CH3CH2OH <=> CH2OH + CH3",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.001, 0.01, 0.1, 1, 10, 100], 'atm'),
        arrhenius = [
            Arrhenius(A=(1.3e+51, 's^-1'), n=-10.59, Ea=(100869, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(5.2e+59, 's^-1'), n=-13.98, Ea=(99850, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.6e+66, 's^-1'), n=-15.3, Ea=(105331, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(5.6e+64, 's^-1'), n=-14.47, Ea=(107039, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.6e+58, 's^-1'), n=-12.29, Ea=(105708, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.8e+47, 's^-1'), n=-8.96, Ea=(101002, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH3CH2OH <=> CH2OH + CH3""",
)

entry(
    index = 301,
    label = "CH3CH2OH <=> C2H5 + OH",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.001, 0.01, 0.1, 1, 10, 100], 'atm'),
        arrhenius = [
            Arrhenius(A=(8.1e+46, 's^-1'), n=-11.33, Ea=(110991, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.8e+56, 's^-1'), n=-13.49, Ea=(107178, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(4.7e+63, 's^-1'), n=-14.99, Ea=(109561, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.5e+65, 's^-1'), n=-14.89, Ea=(112282, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2.8e+61, 's^-1'), n=-13.4, Ea=(113016, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(6.2e+51, 's^-1'), n=-10.34, Ea=(109879, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH3CH2OH <=> C2H5 + OH""",
)

entry(
    index = 302,
    label = "CH3CH2OH <=> C2H4 + H2O",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.001, 0.01, 0.1, 1, 10, 100], 'atm'),
        arrhenius = [
            Arrhenius(A=(3.4e+59, 's^-1'), n=-14.22, Ea=(83625, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2.6e+57, 's^-1'), n=-13.29, Ea=(85214, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.7e+52, 's^-1'), n=-11.52, Ea=(84698, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(5.2e+43, 's^-1'), n=-8.9, Ea=(81461, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(4.6e+32, 's^-1'), n=-5.6, Ea=(76019, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(3.8e+20, 's^-1'), n=-2.06, Ea=(69426, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH3CH2OH <=> C2H4 + H2O""",
)

entry(
    index = 303,
    label = "CH3CHOH + H <=> CH3CH2OH",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.1, 1, 10, 100], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (9.9e+42, 'cm^3/(mol*s)'),
                n = -10.77,
                Ea = (8942, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.9e+55, 'cm^3/(mol*s)'),
                n = -13.56,
                Ea = (14306, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.2e+53, 'cm^3/(mol*s)'),
                n = -12.33,
                Ea = (14505, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.7e+50, 'cm^3/(mol*s)'),
                n = -11.04,
                Ea = (15896, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.6e+40, 'cm^3/(mol*s)'),
                n = -7.82,
                Ea = (12916, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH3CHOH + H <=> CH3CH2OH""",
)

entry(
    index = 304,
    label = "CH3CH2OH + H <=> CH3CHOH + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8800, 'cm^3/(mol*s)'), n=2.68, Ea=(2913, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OH + H <=> CH3CHOH + H2""",
)

entry(
    index = 305,
    label = "CH3CH2OH + H <=> CH2CH2OH + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5300, 'cm^3/(mol*s)'), n=2.81, Ea=(7491, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OH + H <=> CH2CH2OH + H2""",
)

entry(
    index = 306,
    label = "CH3CH2OH + H <=> CH3CH2O + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(950, 'cm^3/(mol*s)'), n=3.14, Ea=(8696, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OH + H <=> CH3CH2O + H2""",
)

entry(
    index = 307,
    label = "CH3CH2OH + O <=> CH2CH2OH + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(970, 'cm^3/(mol*s)'), n=3.23, Ea=(4660, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OH + O <=> CH2CH2OH + OH""",
)

entry(
    index = 308,
    label = "CH3CH2OH + O <=> CH3CHOH + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(150000, 'cm^3/(mol*s)'), n=2.47, Ea=(876, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OH + O <=> CH3CHOH + OH""",
)

entry(
    index = 309,
    label = "CH3CH2OH + O <=> CH3CH2O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0015, 'cm^3/(mol*s)'), n=4.7, Ea=(1730, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OH + O <=> CH3CH2O + OH""",
)

entry(
    index = 310,
    label = "CH3CH2OH + OH <=> CH3CHOH + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(450, 'cm^3/(mol*s)'), n=3.11, Ea=(-2666, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OH + OH <=> CH3CHOH + H2O""",
)

entry(
    index = 311,
    label = "CH3CH2OH + OH <=> CH2CH2OH + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9400, 'cm^3/(mol*s)'), n=2.67, Ea=(-1004, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OH + OH <=> CH2CH2OH + H2O""",
)

entry(
    index = 312,
    label = "CH3CH2OH + OH <=> CH3CH2O + H2O",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(2400, 'cm^3/(mol*s)'), n=2.82, Ea=(-691, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (7.9e+07, 'cm^3/(mol*s)'),
                n = 1.18,
                Ea = (-303, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH3CH2OH + OH <=> CH3CH2O + H2O""",
)

entry(
    index = 313,
    label = "CH3CH2OH + HO2 <=> CH3CHOH + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8200, 'cm^3/(mol*s)'), n=2.55, Ea=(10750, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OH + HO2 <=> CH3CHOH + H2O2""",
)

entry(
    index = 314,
    label = "CH3CH2OH + HO2 <=> CH2CH2OH + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(12000, 'cm^3/(mol*s)'), n=2.55, Ea=(15750, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OH + HO2 <=> CH2CH2OH + H2O2""",
)

entry(
    index = 315,
    label = "CH3CH2OH + HO2 <=> CH3CH2O + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(24000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OH + HO2 <=> CH3CH2O + H2O2""",
)

entry(
    index = 316,
    label = "CH3CH2OH + CH3 <=> CH3CHOH + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(20, 'cm^3/(mol*s)'), n=3.37, Ea=(7630, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OH + CH3 <=> CH3CHOH + CH4""",
)

entry(
    index = 317,
    label = "CH3CH2OH + CH3 <=> CH2CH2OH + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2, 'cm^3/(mol*s)'), n=3.57, Ea=(7717, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OH + CH3 <=> CH2CH2OH + CH4""",
)

entry(
    index = 318,
    label = "CH3CH2OH + CH3 <=> CH3CH2O + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(330, 'cm^3/(mol*s)'), n=3.3, Ea=(12283, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OH + CH3 <=> CH3CH2O + CH4""",
)

entry(
    index = 319,
    label = "CH3CHOH <=> CH3CHO + H",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(6.2e+09, 's^-1'), n=1.31, Ea=(33778, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(A=(1.8e+16, 'cm^3/(mol*s)'), n=0, Ea=(20782, 'cal/mol'), T0=(1, 'K')),
        alpha = 0.187,
        T3 = (65.2, 'K'),
        T1 = (2568, 'K'),
        T2 = (41226, 'K'),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is CH3CHOH <=> CH3CHO + H""",
)

entry(
    index = 320,
    label = "CH3CHOH <=> CH2CHOH + H",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(6.4e+09, 's^-1'), n=1.33, Ea=(35974, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(A=(8.2e+14, 'cm^3/(mol*s)'), n=0, Ea=(21517, 'cal/mol'), T0=(1, 'K')),
        alpha = 0.473,
        T3 = (10, 'K'),
        T1 = (2218, 'K'),
        T2 = (2615, 'K'),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is CH3CHOH <=> CH2CHOH + H""",
)

entry(
    index = 321,
    label = "CH3CHOH <=> CH3 + CH2O",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2.2e+09, 's^-1'), n=1.18, Ea=(33987, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(A=(5.9e+15, 'cm^3/(mol*s)'), n=0, Ea=(21333, 'cal/mol'), T0=(1, 'K')),
        alpha = 0.124,
        T3 = (1, 'K'),
        T1 = (1729, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is CH3CHOH <=> CH3 + CH2O""",
)

entry(
    index = 322,
    label = "CH3CHOH + H <=> CH2CHOH + H2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.1e+12, 'cm^3/(mol*s)'),
        n = 0.2728,
        Ea = (-334, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3CHOH + H <=> CH2CHOH + H2""",
)

entry(
    index = 323,
    label = "CH3CHOH + H <=> C2H4 + H2O",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.001, 0.01, 0.1, 1, 10, 100], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (1.2e+17, 'cm^3/(mol*s)'),
                n = -1.166,
                Ea = (284, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.2e+17, 'cm^3/(mol*s)'),
                n = -1.162,
                Ea = (266, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.8e+17, 'cm^3/(mol*s)'),
                n = -1.216,
                Ea = (386, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.6e+20, 'cm^3/(mol*s)'),
                n = -2.079,
                Ea = (3148, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (9.3e+23, 'cm^3/(mol*s)'),
                n = -2.996,
                Ea = (7954, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.6e+20, 'cm^3/(mol*s)'),
                n = -1.812,
                Ea = (9448, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH3CHOH + H <=> C2H4 + H2O""",
)

entry(
    index = 324,
    label = "CH3CHOH + H <=> CH3 + CH2OH",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.001, 0.01, 0.1, 1, 10, 100], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (1.4e+17, 'cm^3/(mol*s)'),
                n = -0.912,
                Ea = (3081, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.5e+17, 'cm^3/(mol*s)'),
                n = -0.923,
                Ea = (3116, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (4.5e+17, 'cm^3/(mol*s)'),
                n = -1.052,
                Ea = (3509, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.3e+20, 'cm^3/(mol*s)'),
                n = -1.795,
                Ea = (5893, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (5e+24, 'cm^3/(mol*s)'),
                n = -2.949,
                Ea = (10754, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (4e+23, 'cm^3/(mol*s)'),
                n = -2.527,
                Ea = (13637, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH3CHOH + H <=> CH3 + CH2OH""",
)

entry(
    index = 325,
    label = "CH3CHOH + H <=> C2H5 + OH",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.001, 0.01, 0.1, 1, 10, 100], 'atm'),
        arrhenius = [
            Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0.021, Ea=(4442, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (4.4e+13, 'cm^3/(mol*s)'),
                n = 0.01,
                Ea = (4476, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.1e+14, 'cm^3/(mol*s)'),
                n = -0.095,
                Ea = (4790, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.6e+16, 'cm^3/(mol*s)'),
                n = -0.697,
                Ea = (6677, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (6.8e+20, 'cm^3/(mol*s)'),
                n = -1.943,
                Ea = (11331, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (6.3e+21, 'cm^3/(mol*s)'),
                n = -2.106,
                Ea = (15269, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH3CHOH + H <=> C2H5 + OH""",
)

entry(
    index = 326,
    label = "CH3CHOH + O <=> CH3CHO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHOH + O <=> CH3CHO + OH""",
)

entry(
    index = 327,
    label = "CH3CHOH + OH <=> CH3CHO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHOH + OH <=> CH3CHO + H2O""",
)

entry(
    index = 328,
    label = "CH3CHOH + HO2 <=> CH3CHO + OH + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHOH + HO2 <=> CH3CHO + OH + OH""",
)

entry(
    index = 329,
    label = "CH3CHOH + O2 <=> CH3CHO + HO2",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.001, 0.01, 0.1, 1, 10, 100], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (5.3e+17, 'cm^3/(mol*s)'),
                n = -1.637,
                Ea = (838, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (5.3e+17, 'cm^3/(mol*s)'),
                n = -1.637,
                Ea = (838, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (5.3e+17, 'cm^3/(mol*s)'),
                n = -1.637,
                Ea = (838, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (5.3e+17, 'cm^3/(mol*s)'),
                n = -1.638,
                Ea = (839, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.5e+18, 'cm^3/(mol*s)'),
                n = -1.771,
                Ea = (1120, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (3.8e+20, 'cm^3/(mol*s)'),
                n = -2.429,
                Ea = (3090, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH3CHOH + O2 <=> CH3CHO + HO2""",
)

entry(
    index = 330,
    label = "CH3CHOH + O2 <=> CH2CHOH + HO2",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.001, 0.01, 0.1, 1, 10, 100], 'atm'),
        arrhenius = [
            Arrhenius(A=(510, 'cm^3/(mol*s)'), n=2.495, Ea=(-414, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(510, 'cm^3/(mol*s)'), n=2.496, Ea=(-414, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(530, 'cm^3/(mol*s)'), n=2.49, Ea=(-402, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(760, 'cm^3/(mol*s)'), n=2.45, Ea=(-296, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(8900, 'cm^3/(mol*s)'), n=2.146, Ea=(470, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (440000, 'cm^3/(mol*s)'),
                n = 1.699,
                Ea = (2330, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH3CHOH + O2 <=> CH2CHOH + HO2""",
)

entry(
    index = 331,
    label = "CH2CH2OH <=> CH2CHOH + H",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.0013, 1, 100], 'atm'),
        arrhenius = [
            Arrhenius(A=(2.7e+15, 's^-1'), n=-1.92, Ea=(29383, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(3.3e+28, 's^-1'), n=-5.26, Ea=(35583, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2.7e+27, 's^-1'), n=-4.44, Ea=(37205, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CH2OH <=> CH2CHOH + H""",
)

entry(
    index = 332,
    label = "CH2CH2OH + H <=> C2H4 + H2O",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.001, 0.01, 0.1, 1, 10, 100], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (1.7e+17, 'cm^3/(mol*s)'),
                n = -1.184,
                Ea = (335, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.6e+17, 'cm^3/(mol*s)'),
                n = -1.176,
                Ea = (299, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.7e+18, 'cm^3/(mol*s)'),
                n = -1.461,
                Ea = (1107, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.6e+22, 'cm^3/(mol*s)'),
                n = -2.599,
                Ea = (5235, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (6.5e+23, 'cm^3/(mol*s)'),
                n = -2.883,
                Ea = (9307, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (3.6e+16, 'cm^3/(mol*s)'),
                n = -0.716,
                Ea = (8767, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CH2OH + H <=> C2H4 + H2O""",
)

entry(
    index = 333,
    label = "CH2CH2OH + H <=> CH3 + CH2OH",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.001, 0.01, 0.1, 1, 10, 100], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (1.5e+17, 'cm^3/(mol*s)'),
                n = -0.903,
                Ea = (3024, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.9e+17, 'cm^3/(mol*s)'),
                n = -0.935,
                Ea = (3120, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.5e+18, 'cm^3/(mol*s)'),
                n = -1.243,
                Ea = (4062, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.9e+22, 'cm^3/(mol*s)'),
                n = -2.3,
                Ea = (7693, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.8e+25, 'cm^3/(mol*s)'),
                n = -3.1,
                Ea = (12454, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (7.5e+20, 'cm^3/(mol*s)'),
                n = -1.693,
                Ea = (13429, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CH2OH + H <=> CH3 + CH2OH""",
)

entry(
    index = 334,
    label = "CH2CH2OH + H <=> C2H5 + OH",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.001, 0.01, 0.1, 1, 10, 100], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (3.6e+13, 'cm^3/(mol*s)'),
                n = 0.05139,
                Ea = (4302, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (4.6e+13, 'cm^3/(mol*s)'),
                n = 0.02101,
                Ea = (4392, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (3.4e+14, 'cm^3/(mol*s)'),
                n = -0.21686,
                Ea = (5113, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (9.2e+17, 'cm^3/(mol*s)'),
                n = -1.15762,
                Ea = (8193, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.6e+22, 'cm^3/(mol*s)'),
                n = -2.27331,
                Ea = (13261, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (8.1e+19, 'cm^3/(mol*s)'),
                n = -1.50969,
                Ea = (15534, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CH2OH + H <=> C2H5 + OH""",
)

entry(
    index = 335,
    label = "CH2CH2OH + H <=> CH3CH2OH",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(
            A = (5.2e+17, 'cm^3/(mol*s)'),
            n = -0.99,
            Ea = (1580, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        arrheniusLow = Arrhenius(
            A = (2e+41, 'cm^6/(mol^2*s)'),
            n = -7.08,
            Ea = (6685, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.8422,
        T3 = (125, 'K'),
        T1 = (2219, 'K'),
        T2 = (6882, 'K'),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is CH2CH2OH + H <=> CH3CH2OH""",
)

entry(
    index = 336,
    label = "CH2CH2OH + O <=> CH2O + CH2OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CH2OH + O <=> CH2O + CH2OH""",
)

entry(
    index = 337,
    label = "CH2CH2OH + OH <=> CH2CHOH + H2O",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.001, 0.01, 0.1, 1, 10, 100], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (1.3e+19, 'cm^3/(mol*s)'),
                n = -1.96,
                Ea = (273, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.2e+19, 'cm^3/(mol*s)'),
                n = -1.9533,
                Ea = (239, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (4.1e+19, 'cm^3/(mol*s)'),
                n = -2.1007,
                Ea = (625, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (7.9e+22, 'cm^3/(mol*s)'),
                n = -2.9892,
                Ea = (3863, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.8e+24, 'cm^3/(mol*s)'),
                n = -3.3287,
                Ea = (7749, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (4.7e+18, 'cm^3/(mol*s)'),
                n = -1.5805,
                Ea = (7999, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CH2OH + OH <=> CH2CHOH + H2O""",
)

entry(
    index = 338,
    label = "CH2CH2OH + HO2 <=> CH3CH2OH + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CH2OH + HO2 <=> CH3CH2OH + O2""",
)

entry(
    index = 339,
    label = "CH2CH2OH + HO2 => CH2OH + CH2O + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CH2OH + HO2 => CH2OH + CH2O + OH""",
)

entry(
    index = 340,
    label = "CH2CH2OH + O2 <=> HOCH2CH2OO",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.013, 0.1, 1, 10, 100], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (1.5e+44, 'cm^3/(mol*s)'),
                n = -11.15,
                Ea = (5523, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (4.9e+42, 'cm^3/(mol*s)'),
                n = -10.34,
                Ea = (5913, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (6.4e+38, 'cm^3/(mol*s)'),
                n = -8.77,
                Ea = (5859, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (5.6e+32, 'cm^3/(mol*s)'),
                n = -6.58,
                Ea = (5046, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (4.2e+26, 'cm^3/(mol*s)'),
                n = -4.46,
                Ea = (3940, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CH2OH + O2 <=> HOCH2CH2OO""",
)

entry(
    index = 341,
    label = "CH2CH2OH + O2 <=> CH2CHOH + HO2",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiPDepArrhenius(
        arrhenius = [
            PDepArrhenius(
                pressures = ([0.013, 0.1, 1, 10, 100], 'atm'),
                arrhenius = [
                    Arrhenius(
                        A = (1.3e+53, 'cm^3/(mol*s)'),
                        n = -11.88,
                        Ea = (35927, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=-0.79, Ea=(877, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (3.6e+13, 'cm^3/(mol*s)'),
                        n = -0.88,
                        Ea = (3074, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (4.4e+20, 'cm^3/(mol*s)'),
                        n = -2.85,
                        Ea = (8516, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.9e+30, 'cm^3/(mol*s)'),
                        n = -5.51,
                        Ea = (16616, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            PDepArrhenius(
                pressures = ([0.013, 0.1, 1, 10, 100], 'atm'),
                arrhenius = [
                    Arrhenius(
                        A = (2.3e+10, 'cm^3/(mol*s)'),
                        n = -0.15,
                        Ea = (-791, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2.8e+61, 'cm^3/(mol*s)'),
                        n = -14.17,
                        Ea = (43492, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(6000, 'cm^3/(mol*s)'), n=-10, Ea=(199, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(6000, 'cm^3/(mol*s)'), n=-10, Ea=(199, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(6000, 'cm^3/(mol*s)'), n=-10, Ea=(199, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CH2OH + O2 <=> CH2CHOH + HO2""",
)

entry(
    index = 342,
    label = "CH2CH2OH + O2 <=> CH2O + CH2O + OH",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.013, 0.1, 1, 10, 100], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (5.6e+22, 'cm^3/(mol*s)'),
                n = -3.95,
                Ea = (1210, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.4e+24, 'cm^3/(mol*s)'),
                n = -4.31,
                Ea = (2664, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (4.4e+24, 'cm^3/(mol*s)'),
                n = -4.36,
                Ea = (4396, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(A=(3e+25, 'cm^3/(mol*s)'), n=-4.5, Ea=(6763, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (1.2e+29, 'cm^3/(mol*s)'),
                n = -5.44,
                Ea = (11323, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CH2OH + O2 <=> CH2O + CH2O + OH""",
)

entry(
    index = 343,
    label = "CH3CH2O <=> CH3 + CH2O",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(6.3e+10, 's^-1'), n=0.93, Ea=(17098, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4.7e+25, 'cm^3/(mol*s)'),
            n = 0.93,
            Ea = (16532, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.426,
        T3 = (0.3, 'K'),
        T1 = (2278, 'K'),
        T2 = (100000, 'K'),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is CH3CH2O <=> CH3 + CH2O""",
)

entry(
    index = 344,
    label = "CH3CHO + H <=> CH3CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.6e+07, 'cm^3/(mol*s)'),
        n = 1.71,
        Ea = (7090, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3CHO + H <=> CH3CH2O""",
)

entry(
    index = 345,
    label = "CH3CH2O + H <=> CH2OH + CH3",
    degeneracy = 1,
    kinetics = Lindemann(
        arrheniusHigh = Arrhenius(
            A = (2.6e+18, 'cm^3/(mol*s)'),
            n = -1.05,
            Ea = (5128, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        arrheniusLow = Arrhenius(A=(3e+11, 'cm^6/(mol^2*s)'), n=0.893, Ea=(17, 'cal/mol'), T0=(1, 'K')),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is CH3CH2O + H <=> CH2OH + CH3""",
)

entry(
    index = 346,
    label = "CH3CH2O + H <=> CH3CHO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.5e+09, 'cm^3/(mol*s)'), n=1.15, Ea=(673, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2O + H <=> CH3CHO + H2""",
)

entry(
    index = 347,
    label = "CH3CH2O + OH <=> CH3CHO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2O + OH <=> CH3CHO + H2O""",
)

entry(
    index = 348,
    label = "CH3CH2O + O2 <=> CH3CHO + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.5e+10, 'cm^3/(mol*s)'), n=0, Ea=(645, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2O + O2 <=> CH3CHO + HO2""",
)

entry(
    index = 349,
    label = "CH3CH2O + CO <=> C2H5 + CO2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (9.5e+25, 'cm^3/(mol*s)'),
        n = -4.93,
        Ea = (9080, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3CH2O + CO <=> C2H5 + CO2""",
)

entry(
    index = 350,
    label = "CH3CHO <=> CH3 + HCO",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2.7e+22, 's^-1'), n=-1.74, Ea=(86355, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.1e+59, 'cm^3/(mol*s)'),
            n = -11.3,
            Ea = (95912, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.18,
        T3 = (462, 'K'),
        T1 = (167730, 'K'),
        T2 = (1.58e+06, 'K'),
        efficiencies = {'C': 4.23, 'O=C=O': 2.86, 'CC': 4.23, 'O': 8.57, '[H][H]': 2.86, 'N#N': 1.43, '[C-]#[O+]': 2.14},
    ),
    shortDesc = u"""The chemkin file reaction is CH3CHO <=> CH3 + HCO""",
)

entry(
    index = 351,
    label = "CH3CHO <=> CH2CO + H2",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.05, 0.1, 1, 10], 'atm'),
        arrhenius = [
            Arrhenius(A=(4e+44, 's^-1'), n=-10.07, Ea=(87428, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(7.4e+44, 's^-1'), n=-10.05, Ea=(88422, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(8.5e+44, 's^-1'), n=-9.77, Ea=(90905, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2.2e+45, 's^-1'), n=-9.55, Ea=(94879, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH3CHO <=> CH2CO + H2""",
)

entry(
    index = 352,
    label = "CH3CHO <=> CH4 + CO",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.05, 0.1, 1, 10], 'atm'),
        arrhenius = [
            Arrhenius(A=(5.1e+45, 's^-1'), n=-9.85, Ea=(89018, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.4e+45, 's^-1'), n=-9.65, Ea=(87925, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.9e+45, 's^-1'), n=-9.43, Ea=(89415, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.6e+45, 's^-1'), n=-9.1, Ea=(92793, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH3CHO <=> CH4 + CO""",
)

entry(
    index = 353,
    label = "CH3CHO <=> CH2CHOH",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.05, 0.1, 1, 10], 'atm'),
        arrhenius = [
            Arrhenius(A=(7.3e+45, 's^-1'), n=-10.04, Ea=(78785, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2.9e+45, 's^-1'), n=-9.86, Ea=(78884, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.1e+46, 's^-1'), n=-9.76, Ea=(81964, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2.8e+45, 's^-1'), n=-9.35, Ea=(84645, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH3CHO <=> CH2CHOH""",
)

entry(
    index = 354,
    label = "CH3CHO + H <=> CH3CO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(130000, 'cm^3/(mol*s)'), n=2.58, Ea=(1219, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHO + H <=> CH3CO + H2""",
)

entry(
    index = 355,
    label = "CH3CHO + H <=> CH2CHO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2700, 'cm^3/(mol*s)'), n=3.1, Ea=(5203, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHO + H <=> CH2CHO + H2""",
)

entry(
    index = 356,
    label = "CH3CHO + O <=> CH3CO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(1808, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHO + O <=> CH3CO + OH""",
)

entry(
    index = 357,
    label = "CH3CHO + O <=> CH2CHO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(90000, 'cm^3/(mol*s)'), n=2.8, Ea=(5800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHO + O <=> CH2CHO + OH""",
)

entry(
    index = 358,
    label = "CH3CHO + OH <=> CH3CO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(-709, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHO + OH <=> CH3CO + H2O""",
)

entry(
    index = 359,
    label = "CH3CHO + OH <=> CH2CHO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.5e+13, 'cm^3/(mol*s)'), n=0, Ea=(5313, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHO + OH <=> CH2CHO + H2O""",
)

entry(
    index = 360,
    label = "CH3CHO + HO2 <=> CH3CO + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+13, 'cm^3/(mol*s)'), n=0, Ea=(16293, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHO + HO2 <=> CH3CO + H2O2""",
)

entry(
    index = 361,
    label = "CH3CHO + HO2 <=> CH2CHO + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.1e+13, 'cm^3/(mol*s)'), n=0, Ea=(23248, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHO + HO2 <=> CH2CHO + H2O2""",
)

entry(
    index = 362,
    label = "CH3CHO + O2 <=> CH3CO + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(120000, 'cm^3/(mol*s)'), n=2.5, Ea=(37554, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHO + O2 <=> CH3CO + HO2""",
)

entry(
    index = 363,
    label = "CH3CHO + O2 <=> CH2CHO + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.5e+07, 'cm^3/(mol*s)'),
        n = 1.9,
        Ea = (49548, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3CHO + O2 <=> CH2CHO + HO2""",
)

entry(
    index = 364,
    label = "CH3CHO + CH3 <=> CH3CO + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.5e-08, 'cm^3/(mol*s)'),
        n = 6.21,
        Ea = (1629, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3CHO + CH3 <=> CH3CO + CH4""",
)

entry(
    index = 365,
    label = "CH3CHO + CH3 <=> CH2CHO + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.18, 'cm^3/(mol*s)'), n=3.44, Ea=(10384, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHO + CH3 <=> CH2CHO + CH4""",
)

entry(
    index = 366,
    label = "CH3CHO + CH3O <=> CH3CO + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(2981, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHO + CH3O <=> CH3CO + CH3OH""",
)

entry(
    index = 367,
    label = "CH3CHO + CH3O <=> CH2CHO + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+11, 'cm^3/(mol*s)'), n=0, Ea=(7100, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHO + CH3O <=> CH2CHO + CH3OH""",
)

entry(
    index = 368,
    label = "CH3CHO + CH3OO <=> CH3CO + CH3OOH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+13, 'cm^3/(mol*s)'), n=0, Ea=(16293, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHO + CH3OO <=> CH3CO + CH3OOH""",
)

entry(
    index = 369,
    label = "CH3CHO + CH3OO <=> CH2CHO + CH3OOH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.1e+13, 'cm^3/(mol*s)'), n=0, Ea=(23248, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHO + CH3OO <=> CH2CHO + CH3OOH""",
)

entry(
    index = 370,
    label = "CH3CHO + C2H3 <=> CH2CHO + C2H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.2, 'cm^3/(mol*s)'), n=3.96, Ea=(25990, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHO + C2H3 <=> CH2CHO + C2H4""",
)

entry(
    index = 371,
    label = "CH3CHO + C2H3 <=> CH3CO + C2H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.16, 'cm^3/(mol*s)'), n=3.62, Ea=(16810, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHO + C2H3 <=> CH3CO + C2H4""",
)

entry(
    index = 372,
    label = "cC2H4O <=> CH2CHO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+13, 's^-1'), n=0.2, Ea=(71780, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC2H4O <=> CH2CHO + H""",
)

entry(
    index = 373,
    label = "cC2H4O <=> CH3 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.6e+13, 's^-1'), n=0.4, Ea=(61880, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC2H4O <=> CH3 + HCO""",
)

entry(
    index = 374,
    label = "cC2H4O <=> CH3CO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 's^-1'), n=0.25, Ea=(65310, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC2H4O <=> CH3CO + H""",
)

entry(
    index = 375,
    label = "cC2H4O <=> CH2CO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.6e+12, 's^-1'), n=-0.2, Ea=(63030, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC2H4O <=> CH2CO + H2""",
)

entry(
    index = 376,
    label = "cC2H4O <=> CH3CHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.2e+12, 's^-1'), n=-0.75, Ea=(46424, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC2H4O <=> CH3CHO""",
)

entry(
    index = 377,
    label = "cC2H4O <=> C2H2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.6e+12, 's^-1'), n=0.06, Ea=(69530, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC2H4O <=> C2H2 + H2O""",
)

entry(
    index = 378,
    label = "cC2H4O + H <=> CH3CHO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(10950, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC2H4O + H <=> CH3CHO + H""",
)

entry(
    index = 379,
    label = "cC2H4O + H <=> cC2H3O + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(8310, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC2H4O + H <=> cC2H3O + H2""",
)

entry(
    index = 380,
    label = "cC2H4O + H <=> C2H3 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+09, 'cm^3/(mol*s)'), n=0, Ea=(5000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC2H4O + H <=> C2H3 + H2O""",
)

entry(
    index = 381,
    label = "cC2H4O + H <=> C2H4 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.5e+10, 'cm^3/(mol*s)'), n=0, Ea=(5000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC2H4O + H <=> C2H4 + OH""",
)

entry(
    index = 382,
    label = "cC2H4O + O <=> cC2H3O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.9e+12, 'cm^3/(mol*s)'), n=0, Ea=(5250, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC2H4O + O <=> cC2H3O + OH""",
)

entry(
    index = 383,
    label = "cC2H4O + OH <=> cC2H3O + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(3610, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC2H4O + OH <=> cC2H3O + H2O""",
)

entry(
    index = 384,
    label = "cC2H4O + HO2 <=> cC2H3O + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+12, 'cm^3/(mol*s)'), n=0, Ea=(17000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC2H4O + HO2 <=> cC2H3O + H2O2""",
)

entry(
    index = 385,
    label = "cC2H4O + O2 <=> cC2H3O + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(61500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC2H4O + O2 <=> cC2H3O + HO2""",
)

entry(
    index = 386,
    label = "cC2H4O + CH3 <=> cC2H3O + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.1e+12, 'cm^3/(mol*s)'), n=0, Ea=(11830, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC2H4O + CH3 <=> cC2H3O + CH4""",
)

entry(
    index = 387,
    label = "CH2CHOH + H <=> CH2CHO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1500, 'cm^3/(mol*s)'), n=3.077, Ea=(7230, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHOH + H <=> CH2CHO + H2""",
)

entry(
    index = 388,
    label = "CH2CHOH + H <=> CHCHOH + H2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.5e+07, 'cm^3/(mol*s)'),
        n = 2.03,
        Ea = (15180, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHOH + H <=> CHCHOH + H2""",
)

entry(
    index = 389,
    label = "CH2CHOH + O <=> CH2OH + HCO",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(3.9e+12, 'cm^3/(mol*s)'), n=0, Ea=(1494, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(6.2e+13, 'cm^3/(mol*s)'), n=0, Ea=(6855, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHOH + O <=> CH2OH + HCO""",
)

entry(
    index = 390,
    label = "CH2CHOH + O <=> CH2OH + H + CO",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(
                A = (3.7e+23, 'cm^3/(mol*s)'),
                n = -2.473,
                Ea = (26782, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.3e+22, 'cm^3/(mol*s)'),
                n = -2.473,
                Ea = (21421, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHOH + O <=> CH2OH + H + CO""",
)

entry(
    index = 391,
    label = "CH2CHOH + O <=> CH2CHO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.6e+07, 'cm^3/(mol*s)'), n=2, Ea=(4400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHOH + O <=> CH2CHO + OH""",
)

entry(
    index = 392,
    label = "CH2CHOH + OH <=> CHCHOH + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.13, 'cm^3/(mol*s)'), n=4.2, Ea=(-860, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHOH + OH <=> CHCHOH + H2O""",
)

entry(
    index = 393,
    label = "CH2CHOH + OH <=> CH2CHO + H2O",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(2400, 'cm^3/(mol*s)'), n=2.82, Ea=(-691, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (7.9e+07, 'cm^3/(mol*s)'),
                n = 1.18,
                Ea = (-303, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHOH + OH <=> CH2CHO + H2O""",
)

entry(
    index = 394,
    label = "CH2CHOH + HO2 <=> CH2CHO + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.6e+12, 'cm^3/(mol*s)'), n=0, Ea=(16293, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHOH + HO2 <=> CH2CHO + H2O2""",
)

entry(
    index = 395,
    label = "CH2CHOH + HO2 <=> CH3CHO + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(150000, 'cm^3/(mol*s)'), n=1.67, Ea=(6810, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHOH + HO2 <=> CH3CHO + HO2""",
)

entry(
    index = 396,
    label = "CH2CHOH + O2 => CH2O + HCO + OH",
    degeneracy = 1,
    duplicate = True,
    reversible = False,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(
                A = (3.5e+07, 'cm^3/(mol*s)'),
                n = 1.8,
                Ea = (39000, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (-2.1e+17, 'cm^3/(mol*s)'),
                n = -0.673,
                Ea = (58927, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHOH + O2 => CH2O + HCO + OH""",
)

entry(
    index = 397,
    label = "CH2CHOH + O2 => CH2O + H + CO + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (2.1e+17, 'cm^3/(mol*s)'),
        n = -0.673,
        Ea = (58927, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHOH + O2 => CH2O + H + CO + OH""",
)

entry(
    index = 398,
    label = "CHCHOH <=> HCCOH + H",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.04, 1, 10, 100], 'atm'),
        arrhenius = [
            Arrhenius(A=(4.4e+29, 's^-1'), n=-6.153, Ea=(51383, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.1e+31, 's^-1'), n=-6.153, Ea=(51383, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.5e+32, 's^-1'), n=-6.168, Ea=(52239, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(5.5e+29, 's^-1'), n=-5.057, Ea=(52377, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CHCHOH <=> HCCOH + H""",
)

entry(
    index = 399,
    label = "CHCHOH + H <=> CH2CHO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CHCHOH + H <=> CH2CHO + H""",
)

entry(
    index = 400,
    label = "CHCHOH + H <=> HCCOH + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CHCHOH + H <=> HCCOH + H2""",
)

entry(
    index = 401,
    label = "CHCHOH + O <=> OCHCHO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CHCHOH + O <=> OCHCHO + H""",
)

entry(
    index = 402,
    label = "CHCHOH + OH <=> HCCOH + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CHCHOH + OH <=> HCCOH + H2O""",
)

entry(
    index = 403,
    label = "CHCHOH + O2 <=> OCHCHO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(-187, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CHCHOH + O2 <=> OCHCHO + OH""",
)

entry(
    index = 404,
    label = "CHCHOH + O2 <=> HOCHO + HCO",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(3.3e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(-1.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(-187, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CHCHOH + O2 <=> HOCHO + HCO""",
)

entry(
    index = 405,
    label = "CHCHOH + O2 <=> HOCHO + H + CO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.1e+22, 'cm^3/(mol*s)'),
        n = -2.498,
        Ea = (20266, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CHCHOH + O2 <=> HOCHO + H + CO""",
)

entry(
    index = 406,
    label = "CHCHOH + CH2O <=> CH2CHOH + HCO",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(5400, 'cm^3/(mol*s)'), n=2.81, Ea=(5860, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (-3.2e+13, 'cm^3/(mol*s)'),
                n = 0.337,
                Ea = (25787, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CHCHOH + CH2O <=> CH2CHOH + HCO""",
)

entry(
    index = 407,
    label = "CHCHOH + CH2O <=> CH2CHOH + H + CO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.2e+13, 'cm^3/(mol*s)'),
        n = 0.337,
        Ea = (25787, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CHCHOH + CH2O <=> CH2CHOH + H + CO""",
)

entry(
    index = 408,
    label = "CHCHOH + HCO <=> CH2CHOH + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CHCHOH + HCO <=> CH2CHOH + CO""",
)

entry(
    index = 409,
    label = "CHCHOH + CH3 <=> HCCOH + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CHCHOH + CH3 <=> HCCOH + CH4""",
)

entry(
    index = 410,
    label = "cC2H3O <=> CH2CHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.7e+31, 's^-1'), n=-6.9, Ea=(14994, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC2H3O <=> CH2CHO""",
)

entry(
    index = 411,
    label = "cC2H3O <=> CH2CO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 's^-1'), n=0, Ea=(14863, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC2H3O <=> CH2CO + H""",
)

entry(
    index = 412,
    label = "cC2H3O <=> CH3 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.1e+12, 's^-1'), n=0, Ea=(14280, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC2H3O <=> CH3 + CO""",
)

entry(
    index = 413,
    label = "CH3CO <=> CH3 + CO",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.025, 0.1, 1, 10, 100], 'atm'),
        arrhenius = [
            Arrhenius(A=(6.9e+14, 's^-1'), n=-1.97, Ea=(14585, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2.4e+15, 's^-1'), n=-2, Ea=(14805, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2e+16, 's^-1'), n=-2.09, Ea=(15197, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(6.5e+18, 's^-1'), n=-2.52, Ea=(16436, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(8.2e+19, 's^-1'), n=-2.55, Ea=(17263, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.3e+20, 's^-1'), n=-2.32, Ea=(18012, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH3CO <=> CH3 + CO""",
)

entry(
    index = 414,
    label = "CH2CO + H <=> CH3CO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.3e+08, 'cm^3/(mol*s)'),
        n = 1.61,
        Ea = (2627, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2CO + H <=> CH3CO""",
)

entry(
    index = 415,
    label = "CH3CO + H <=> CH3 + HCO",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(2.1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (-1.2e+23, 'cm^3/(mol*s)'),
                n = -2.473,
                Ea = (19927, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH3CO + H <=> CH3 + HCO""",
)

entry(
    index = 416,
    label = "CH3CO + H <=> CH3 + H + CO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.2e+23, 'cm^3/(mol*s)'),
        n = -2.473,
        Ea = (19927, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3CO + H <=> CH3 + H + CO""",
)

entry(
    index = 417,
    label = "CH3CO + H <=> CH2CO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CO + H <=> CH2CO + H2""",
)

entry(
    index = 418,
    label = "CH3CO + O <=> CH3 + CO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.6e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CO + O <=> CH3 + CO2""",
)

entry(
    index = 419,
    label = "CH3CO + O <=> CH2CO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CO + O <=> CH2CO + OH""",
)

entry(
    index = 420,
    label = "CH3CO + OH <=> CH2CO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CO + OH <=> CH2CO + H2O""",
)

entry(
    index = 421,
    label = "CH3CO + O2 <=> CH3C(O)OO",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.1, 1, 10], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (3.6e+31, 'cm^3/(mol*s)'),
                n = -4.769,
                Ea = (2188, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.8e+34, 'cm^3/(mol*s)'),
                n = -7.21,
                Ea = (6060, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (4.9e+31, 'cm^3/(mol*s)'),
                n = -6.087,
                Ea = (6541, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH3CO + O2 <=> CH3C(O)OO""",
)

entry(
    index = 422,
    label = "CH3CO + O2 <=> CH2CO + HO2",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.1, 1, 10], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (1.3e+08, 'cm^3/(mol*s)'),
                n = 1.986,
                Ea = (228, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (3.6e+10, 'cm^3/(mol*s)'),
                n = 0.544,
                Ea = (3721, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (7.7e+13, 'cm^3/(mol*s)'),
                n = -0.335,
                Ea = (7510, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH3CO + O2 <=> CH2CO + HO2""",
)

entry(
    index = 423,
    label = "CH3CO + O2 <=> CH2O + CO + OH",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.1, 1, 10], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (5.1e+22, 'cm^3/(mol*s)'),
                n = -3.524,
                Ea = (3255, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (4.9e+23, 'cm^3/(mol*s)'),
                n = -3.712,
                Ea = (5895, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (4.8e+22, 'cm^3/(mol*s)'),
                n = -3.303,
                Ea = (8598, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH3CO + O2 <=> CH2O + CO + OH""",
)

entry(
    index = 424,
    label = "CH3CO + CH3 <=> C2H6 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CO + CH3 <=> C2H6 + CO""",
)

entry(
    index = 425,
    label = "CH3CO + CH3 <=> CH2CO + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CO + CH3 <=> CH2CO + CH4""",
)

entry(
    index = 426,
    label = "CH3CO + CH3OO <=> CH3 + CO2 + CH3O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CO + CH3OO <=> CH3 + CO2 + CH3O""",
)

entry(
    index = 427,
    label = "CH2CHO <=> CH2CO + H",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.1, 1, 10, 100], 'atm'),
        arrhenius = [
            Arrhenius(A=(2.4e+25, 's^-1'), n=-4.8, Ea=(43424, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2.4e+30, 's^-1'), n=-5.86, Ea=(46114, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.3e+34, 's^-1'), n=-6.57, Ea=(49454, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(3.5e+36, 's^-1'), n=-6.92, Ea=(52979, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.2e+36, 's^-1'), n=-6.48, Ea=(55191, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHO <=> CH2CO + H""",
)

entry(
    index = 428,
    label = "CH2CHO <=> CH3 + CO",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.025, 0.1, 1, 10, 100], 'atm'),
        arrhenius = [
            Arrhenius(A=(1.2e+30, 's^-1'), n=-6.07, Ea=(41332, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.5e+31, 's^-1'), n=-6.27, Ea=(42478, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(6.4e+32, 's^-1'), n=-6.57, Ea=(44282, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(6.5e+34, 's^-1'), n=-6.87, Ea=(47191, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2.2e+35, 's^-1'), n=-6.76, Ea=(49548, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2.2e+33, 's^-1'), n=-5.97, Ea=(50448, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHO <=> CH3 + CO""",
)

entry(
    index = 429,
    label = "CH2CHO + H <=> CH3 + HCO",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiPDepArrhenius(
        arrhenius = [
            PDepArrhenius(
                pressures = ([0.01, 0.1, 1, 10, 100], 'atm'),
                arrhenius = [
                    Arrhenius(
                        A = (5e+20, 'cm^3/(mol*s)'),
                        n = -2.063,
                        Ea = (3994, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (6.5e+21, 'cm^3/(mol*s)'),
                        n = -2.371,
                        Ea = (4936, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1e+33, 'cm^3/(mol*s)'),
                        n = -5.363,
                        Ea = (17064, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.8e+40, 'cm^3/(mol*s)'),
                        n = -7.368,
                        Ea = (24518, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.2e+32, 'cm^3/(mol*s)'),
                        n = -4.93,
                        Ea = (22682, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            PDepArrhenius(
                pressures = ([0.01, 0.1, 1, 10, 100], 'atm'),
                arrhenius = [
                    Arrhenius(
                        A = (3e+28, 'cm^3/(mol*s)'),
                        n = -4.169,
                        Ea = (12362, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (4e+29, 'cm^3/(mol*s)'),
                        n = -4.477,
                        Ea = (13304, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.2e+27, 'cm^3/(mol*s)'),
                        n = -3.851,
                        Ea = (9085, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2.7e+28, 'cm^3/(mol*s)'),
                        n = -4.075,
                        Ea = (13264, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.97e+24, 'cm^3/(mol*s)'),
                        n = -2.822,
                        Ea = (14304.9, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHO + H <=> CH3 + HCO""",
)

entry(
    index = 430,
    label = "CH2CHO + H <=> CH3 + H + CO",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.1, 1, 10, 100], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (1.1e+27, 'cm^3/(mol*s)'),
                n = -3.408,
                Ea = (23047, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.5e+28, 'cm^3/(mol*s)'),
                n = -3.716,
                Ea = (23989, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.1e+32, 'cm^3/(mol*s)'),
                n = -4.773,
                Ea = (27620, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.7e+35, 'cm^3/(mol*s)'),
                n = -5.573,
                Ea = (32381, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (4.5e+30, 'cm^3/(mol*s)'),
                n = -4.166,
                Ea = (33356, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHO + H <=> CH3 + H + CO""",
)

entry(
    index = 431,
    label = "CH2CHO + H <=> CH3CO + H",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.001, 0.01, 0.1, 1, 10, 100], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (3.6e+13, 'cm^3/(mol*s)'),
                n = 0.05139,
                Ea = (4301, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (4.6e+13, 'cm^3/(mol*s)'),
                n = 0.02101,
                Ea = (4392, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (3.4e+14, 'cm^3/(mol*s)'),
                n = -0.21686,
                Ea = (5113, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (9.2e+17, 'cm^3/(mol*s)'),
                n = -1.15762,
                Ea = (8193, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.6e+22, 'cm^3/(mol*s)'),
                n = -2.27331,
                Ea = (13261, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (8.1e+19, 'cm^3/(mol*s)'),
                n = -1.50969,
                Ea = (15534, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHO + H <=> CH3CO + H""",
)

entry(
    index = 432,
    label = "CH2CHO + O <=> CH2O + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHO + O <=> CH2O + HCO""",
)

entry(
    index = 433,
    label = "CH2CHO + O <=> CH2O + H + CO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3e+23, 'cm^3/(mol*s)'),
        n = -2.473,
        Ea = (19927, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHO + O <=> CH2O + H + CO""",
)

entry(
    index = 434,
    label = "CH2CHO + OH <=> CH2CO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHO + OH <=> CH2CO + H2O""",
)

entry(
    index = 435,
    label = "CH2CHO + OH <=> CH2OH + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHO + OH <=> CH2OH + HCO""",
)

entry(
    index = 436,
    label = "CH2CHO + OH <=> CH2OH + H + CO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (6e+22, 'cm^3/(mol*s)'),
        n = -2.473,
        Ea = (19927, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHO + OH <=> CH2OH + H + CO""",
)

entry(
    index = 437,
    label = "CH2CHO + HO2 <=> CH2O + HCO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHO + HO2 <=> CH2O + HCO + OH""",
)

entry(
    index = 438,
    label = "CH2CHO + O2 <=> CH2O + CO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+10, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHO + O2 <=> CH2O + CO + OH""",
)

entry(
    index = 439,
    label = "CH2CHO + CH2 <=> C2H4 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHO + CH2 <=> C2H4 + HCO""",
)

entry(
    index = 440,
    label = "CH2CHO + CH2 <=> C2H4 + H + CO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3e+23, 'cm^3/(mol*s)'),
        n = -2.473,
        Ea = (19927, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHO + CH2 <=> C2H4 + H + CO""",
)

entry(
    index = 441,
    label = "CH2CHO + CH <=> C2H3 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHO + CH <=> C2H3 + HCO""",
)

entry(
    index = 442,
    label = "CH2CHO + CH <=> C2H3 + H + CO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (6e+23, 'cm^3/(mol*s)'),
        n = -2.473,
        Ea = (19927, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHO + CH <=> C2H3 + H + CO""",
)

entry(
    index = 443,
    label = "CHCHO + H <=> CH2CO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CHCHO + H <=> CH2CO + H""",
)

entry(
    index = 444,
    label = "CHCHO + O2 <=> CO2 + H + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.1e+09, 'cm^3/(mol*s)'),
        n = 0.9929,
        Ea = (-269, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CHCHO + O2 <=> CO2 + H + HCO""",
)

entry(
    index = 445,
    label = "CHCHO + O2 <=> OCHCHO + O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.3e+06, 'cm^3/(mol*s)'),
        n = 2.4202,
        Ea = (1604, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CHCHO + O2 <=> OCHCHO + O""",
)

entry(
    index = 446,
    label = "CH2 + CO <=> CH2CO",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(8.1e+11, 'cm^3/(mol*s)'), n=0.5, Ea=(4510, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (2.69e+33, 'cm^6/(mol^2*s)'),
            n = -5.11,
            Ea = (7095, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.5907,
        T3 = (275, 'K'),
        T1 = (1226, 'K'),
        T2 = (5185, 'K'),
        efficiencies = {'O': 6, 'N#N': 1, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is CH2 + CO <=> CH2CO""",
)

entry(
    index = 447,
    label = "CH2CO + H <=> CH3 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (7.8e+08, 'cm^3/(mol*s)'),
        n = 1.45,
        Ea = (2780, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2CO + H <=> CH3 + CO""",
)

entry(
    index = 448,
    label = "CH2CO + H <=> HCCO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+07, 'cm^3/(mol*s)'), n=2, Ea=(10000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CO + H <=> HCCO + H2""",
)

entry(
    index = 449,
    label = "CH2CO + O <=> CO2 + CH2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(1350, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CO + O <=> CO2 + CH2""",
)

entry(
    index = 450,
    label = "CH2CO + OH <=> CH2OH + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1013, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CO + OH <=> CH2OH + CO""",
)

entry(
    index = 451,
    label = "CH2CO + OH <=> CH3 + CO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.7e+11, 'cm^3/(mol*s)'), n=0, Ea=(-1013, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CO + OH <=> CH3 + CO2""",
)

entry(
    index = 452,
    label = "CH2CO + OH <=> HCCO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+07, 'cm^3/(mol*s)'), n=2, Ea=(3000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CO + OH <=> HCCO + H2O""",
)

entry(
    index = 453,
    label = "CH2CO + CH2(S) <=> C2H4 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.6e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CO + CH2(S) <=> C2H4 + CO""",
)

entry(
    index = 454,
    label = "HCCOH + H <=> HCCO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+07, 'cm^3/(mol*s)'), n=2, Ea=(1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCCOH + H <=> HCCO + H2""",
)

entry(
    index = 455,
    label = "HCCOH + OH <=> HCCO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+07, 'cm^3/(mol*s)'), n=2, Ea=(1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCCOH + OH <=> HCCO + H2O""",
)

entry(
    index = 456,
    label = "CH + CO <=> HCCO",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (2.7e+28, 'cm^6/(mol^2*s)'),
            n = -3.74,
            Ea = (1936, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.5757,
        T3 = (237, 'K'),
        T1 = (1652, 'K'),
        T2 = (5069, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, 'N#N': 1, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is CH + CO <=> HCCO""",
)

entry(
    index = 457,
    label = "HCCO + H <=> CH2(S) + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.5e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCCO + H <=> CH2(S) + CO""",
)

entry(
    index = 458,
    label = "HCCO + O <=> CO + CO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCCO + O <=> CO + CO + H""",
)

entry(
    index = 459,
    label = "HCCO + OH <=> C2O + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(14000, 'cm^3/(mol*s)'), n=2.65, Ea=(1472, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCCO + OH <=> C2O + H2O""",
)

entry(
    index = 460,
    label = "HCCO + OH <=> CH2CO + O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.3e+06, 'cm^3/(mol*s)'),
        n = 1.99,
        Ea = (11280, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is HCCO + OH <=> CH2CO + O""",
)

entry(
    index = 461,
    label = "HCCO + OH <=> HCOH + CO",
    degeneracy = 1,
    duplicate = True,
    kinetics = PDepArrhenius(
        pressures = ([1, 10, 100], 'atm'),
        arrhenius = [
            Arrhenius(A=(3e+16, 'cm^3/(mol*s)'), n=-0.935, Ea=(659, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (1.1e+18, 'cm^3/(mol*s)'),
                n = -1.392,
                Ea = (1395, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (3.2e+18, 'cm^3/(mol*s)'),
                n = -1.523,
                Ea = (1627, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is HCCO + OH <=> HCOH + CO""",
)

entry(
    index = 462,
    label = "HCCO + OH <=> HCOH + CO",
    degeneracy = 1,
    duplicate = True,
    kinetics = PDepArrhenius(
        pressures = ([1, 10, 100], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (8.7e+19, 'cm^3/(mol*s)'),
                n = -1.792,
                Ea = (5994, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (3.5e+22, 'cm^3/(mol*s)'),
                n = -2.475,
                Ea = (9163, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.3e+24, 'cm^3/(mol*s)'),
                n = -2.902,
                Ea = (10522, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is HCCO + OH <=> HCOH + CO""",
)

entry(
    index = 463,
    label = "HCCO + OH <=> HCOH + CO",
    degeneracy = 1,
    duplicate = True,
    kinetics = Arrhenius(A=(2.9e+12, 'cm^3/(mol*s)'), n=0.37, Ea=(-24, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCCO + OH <=> HCOH + CO""",
)

entry(
    index = 464,
    label = "HCCO + OH <=> CH2O + CO",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([1, 100], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (1.2e+21, 'cm^3/(mol*s)'),
                n = -2.459,
                Ea = (2528, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(A=(1.1e+08, 'cm^3/(mol*s)'), n=0.11, Ea=(52, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is HCCO + OH <=> CH2O + CO""",
)

entry(
    index = 465,
    label = "HCCO + OH <=> OCHCO + H",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.1, 1, 10, 100], 'atm'),
        arrhenius = [
            Arrhenius(A=(2.6e+08, 'cm^3/(mol*s)'), n=1.41, Ea=(845, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2.6e+08, 'cm^3/(mol*s)'), n=1.41, Ea=(845, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2.6e+08, 'cm^3/(mol*s)'), n=1.41, Ea=(849, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(3e+08, 'cm^3/(mol*s)'), n=1.4, Ea=(917, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (8.2e+08, 'cm^3/(mol*s)'),
                n = 1.28,
                Ea = (1531, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is HCCO + OH <=> OCHCO + H""",
)

entry(
    index = 466,
    label = "HCCO + OH <=> CO2 + CH2",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiPDepArrhenius(
        arrhenius = [
            PDepArrhenius(
                pressures = ([0.1, 1, 10], 'atm'),
                arrhenius = [
                    Arrhenius(
                        A = (1.7e+15, 'cm^3/(mol*s)'),
                        n = -1.19,
                        Ea = (-521, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (7.2e+27, 'cm^3/(mol*s)'),
                        n = -5.023,
                        Ea = (2468, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1e+91, 'cm^3/(mol*s)'),
                        n = -20.137,
                        Ea = (114841, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            PDepArrhenius(
                pressures = ([0.1, 1, 10], 'atm'),
                arrhenius = [
                    Arrhenius(
                        A = (-7.4e+17, 'cm^3/(mol*s)'),
                        n = -1.92,
                        Ea = (1686, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.1e+21, 'cm^3/(mol*s)'),
                        n = -2.28,
                        Ea = (16960, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2.3e+65, 'cm^3/(mol*s)'),
                        n = -16.078,
                        Ea = (19592, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            PDepArrhenius(
                pressures = ([0.01, 0.1, 1, 10], 'atm'),
                arrhenius = [
                    Arrhenius(A=(1e+19, 'cm^3/(mol*s)'), n=-2.08, Ea=(44, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.4e+19, 'cm^3/(mol*s)'), n=-2.12, Ea=(88, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(7.1e+19, 'cm^3/(mol*s)'), n=-2.3, Ea=(824, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (1.8e+20, 'cm^3/(mol*s)'),
                        n = -2.34,
                        Ea = (2421, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is HCCO + OH <=> CO2 + CH2""",
)

entry(
    index = 467,
    label = "HCCO + O2 <=> CO2 + CO + H",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.9e+12, 'cm^3/(mol*s)'),
        n = -0.142,
        Ea = (1150, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is HCCO + O2 <=> CO2 + CO + H""",
)

entry(
    index = 468,
    label = "HCCO + O2 <=> CO + CO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.6e+11, 'cm^3/(mol*s)'),
        n = -0.02,
        Ea = (1020, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is HCCO + O2 <=> CO + CO + OH""",
)

entry(
    index = 469,
    label = "HCCO + O2 <=> HCO + CO + O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(220, 'cm^3/(mol*s)'), n=2.69, Ea=(3540, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCCO + O2 <=> HCO + CO + O""",
)

entry(
    index = 470,
    label = "HCCO + O2 <=> H + CO + CO + O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.3e+12, 'cm^3/(mol*s)'),
        n = 0.217,
        Ea = (23467, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is HCCO + O2 <=> H + CO + CO + O""",
)

entry(
    index = 471,
    label = "HCCO + CH2 <=> C2H3 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCCO + CH2 <=> C2H3 + CO""",
)

entry(
    index = 472,
    label = "HCCO + CH <=> C2H2 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCCO + CH <=> C2H2 + CO""",
)

entry(
    index = 473,
    label = "HCCO + HCCO <=> C2H2 + CO + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCCO + HCCO <=> C2H2 + CO + CO""",
)

entry(
    index = 474,
    label = "C2O <=> C + CO",
    degeneracy = 1,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(A=(2e+15, 'cm^3/(mol*s)'), n=0, Ea=(44200, 'cal/mol'), T0=(1, 'K')),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is C2O <=> C + CO""",
)

entry(
    index = 475,
    label = "C2O + H <=> CH + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2O + H <=> CH + CO""",
)

entry(
    index = 476,
    label = "C2O + O <=> CO + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2O + O <=> CO + CO""",
)

entry(
    index = 477,
    label = "C2O + OH <=> CO + CO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2O + OH <=> CO + CO + H""",
)

entry(
    index = 478,
    label = "C2O + O2 <=> CO + CO + O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(2600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2O + O2 <=> CO + CO + O""",
)

entry(
    index = 479,
    label = "C2O + O2 <=> CO + CO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(2600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2O + O2 <=> CO + CO2""",
)

entry(
    index = 480,
    label = "C2O + C <=> CO + C2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2O + C <=> CO + C2""",
)

entry(
    index = 481,
    label = "CH3CH2OOH <=> CH3CH2O + OH",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.1, 1, 100], 'atm'),
        arrhenius = [
            Arrhenius(A=(6.1e+58, 's^-1'), n=-14.05, Ea=(54131, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(9.3e+52, 's^-1'), n=-11.91, Ea=(53378, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.4e+33, 's^-1'), n=-5.27, Ea=(48696, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH3CH2OOH <=> CH3CH2O + OH""",
)

entry(
    index = 482,
    label = "CH3CH2OOH + H <=> CH3CHOOH + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.5e+10, 'cm^3/(mol*s)'), n=0, Ea=(1860, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OOH + H <=> CH3CHOOH + H2""",
)

entry(
    index = 483,
    label = "CH3CH2OOH + H <=> CH3CH2OO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.3e+10, 'cm^3/(mol*s)'), n=0, Ea=(1860, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OOH + H <=> CH3CH2OO + H2""",
)

entry(
    index = 484,
    label = "CH3CH2OOH + H <=> CH3CH2O + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+10, 'cm^3/(mol*s)'), n=0, Ea=(1860, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OOH + H <=> CH3CH2O + H2O""",
)

entry(
    index = 485,
    label = "CH3CH2OOH + O <=> CH3CHOOH + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(4750, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OOH + O <=> CH3CHOOH + OH""",
)

entry(
    index = 486,
    label = "CH3CH2OOH + O <=> CH3CH2OO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.7e+12, 'cm^3/(mol*s)'), n=0, Ea=(4750, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OOH + O <=> CH3CH2OO + OH""",
)

entry(
    index = 487,
    label = "CH3CH2OOH + OH <=> CH3CHOOH + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.2e+11, 'cm^3/(mol*s)'), n=0, Ea=(-258, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OOH + OH <=> CH3CHOOH + H2O""",
)

entry(
    index = 488,
    label = "CH3CH2OOH + OH <=> CH3CH2OO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.1e+12, 'cm^3/(mol*s)'), n=0, Ea=(-437, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OOH + OH <=> CH3CH2OO + H2O""",
)

entry(
    index = 489,
    label = "CH3CH2OOH + HO2 <=> CH3CH2OO + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(41000, 'cm^3/(mol*s)'), n=2.5, Ea=(10206, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OOH + HO2 <=> CH3CH2OO + H2O2""",
)

entry(
    index = 490,
    label = "CH3CHOOH <=> CH3CHO + OH",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([1, 10, 100], 'atm'),
        arrhenius = [
            Arrhenius(A=(3.5e+12, 's^-1'), n=-0.947, Ea=(979, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(3.5e+13, 's^-1'), n=-0.947, Ea=(980, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(5.8e+14, 's^-1'), n=-1.012, Ea=(1068, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH3CHOOH <=> CH3CHO + OH""",
)

entry(
    index = 491,
    label = "CH3CH2OO <=> CH2CH2OOH",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100], 'atm'),
        arrhenius = [
            Arrhenius(A=(3.2e+31, 's^-1'), n=-8.25, Ea=(29360, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(3.5e+30, 's^-1'), n=-7.88, Ea=(29330, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.5e+29, 's^-1'), n=-7.37, Ea=(29210, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(3.5e+27, 's^-1'), n=-6.77, Ea=(29000, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(3.6e+25, 's^-1'), n=-6.04, Ea=(28780, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.6e+24, 's^-1'), n=-5.51, Ea=(28800, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.4e+21, 's^-1'), n=-4.4, Ea=(28410, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2.9e+19, 's^-1'), n=-3.73, Ea=(28490, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.3e+17, 's^-1'), n=-2.81, Ea=(28500, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(5.3e+14, 's^-1'), n=-1.9, Ea=(28470, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(4.7e+13, 's^-1'), n=-1.4, Ea=(28970, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(4.2e+12, 's^-1'), n=-0.92, Ea=(29380, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.9e+08, 's^-1'), n=0.57, Ea=(28590, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH3CH2OO <=> CH2CH2OOH""",
)

entry(
    index = 492,
    label = "CH3CH2OO <=> C2H4 + HO2",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100], 'atm'),
        arrhenius = [
            Arrhenius(A=(1.9e+46, 's^-1'), n=-11.85, Ea=(36440, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(4.2e+46, 's^-1'), n=-11.88, Ea=(36820, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(3.6e+46, 's^-1'), n=-11.77, Ea=(37100, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.7e+46, 's^-1'), n=-11.58, Ea=(37330, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(4.4e+45, 's^-1'), n=-11.28, Ea=(37570, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(8.1e+44, 's^-1'), n=-10.94, Ea=(37780, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(4.6e+43, 's^-1'), n=-10.43, Ea=(37910, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(8.7e+41, 's^-1'), n=-9.77, Ea=(37860, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(8.7e+39, 's^-1'), n=-9.01, Ea=(37780, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(7.2e+36, 's^-1'), n=-7.95, Ea=(37240, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(4.3e+33, 's^-1'), n=-6.84, Ea=(36660, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.9e+30, 's^-1'), n=-5.71, Ea=(35910, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.6e+26, 's^-1'), n=-4.37, Ea=(34840, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH3CH2OO <=> C2H4 + HO2""",
)

entry(
    index = 493,
    label = "CH3CH2OO <=> cC2H4O + OH",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100], 'atm'),
        arrhenius = [
            Arrhenius(A=(2e+49, 's^-1'), n=-13.32, Ea=(38820, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.7e+50, 's^-1'), n=-13.52, Ea=(39510, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(6.6e+50, 's^-1'), n=-13.62, Ea=(40180, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(3.9e+48, 's^-1'), n=-12.85, Ea=(39830, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(7.6e+48, 's^-1'), n=-12.82, Ea=(40620, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(9.7e+46, 's^-1'), n=-12.11, Ea=(40640, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(9.9e+46, 's^-1'), n=-11.94, Ea=(41670, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.3e+45, 's^-1'), n=-11.2, Ea=(42020, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.5e+44, 's^-1'), n=-10.71, Ea=(43040, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1e+42, 's^-1'), n=-9.86, Ea=(43640, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2.5e+39, 's^-1'), n=-8.87, Ea=(44290, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.9e+36, 's^-1'), n=-7.75, Ea=(44660, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2.5e+31, 's^-1'), n=-6.1, Ea=(44560, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH3CH2OO <=> cC2H4O + OH""",
)

entry(
    index = 494,
    label = "CH3CH2OO + H <=> CH3CH2O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OO + H <=> CH3CH2O + OH""",
)

entry(
    index = 495,
    label = "CH3CH2OO + O <=> CH3CH2O + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.9e+10, 'cm^3/(mol*s)'), n=1, Ea=(-724, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OO + O <=> CH3CH2O + O2""",
)

entry(
    index = 496,
    label = "CH3CH2OO + OH <=> CH3CH2OH + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.1e+15, 'cm^3/(mol*s)'), n=-0.81, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OO + OH <=> CH3CH2OH + O2""",
)

entry(
    index = 497,
    label = "CH3CH2OO + HO2 <=> CH3CH2OOH + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.6e+11, 'cm^3/(mol*s)'), n=0, Ea=(-1267, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OO + HO2 <=> CH3CH2OOH + O2""",
)

entry(
    index = 498,
    label = "CH3CH2OO + CO <=> CH3CH2O + CO2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (160000, 'cm^3/(mol*s)'),
        n = 2.18,
        Ea = (17940, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3CH2OO + CO <=> CH3CH2O + CO2""",
)

entry(
    index = 499,
    label = "CH3CH2OO + CH3 <=> CH3CH2O + CH3O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.1e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1411, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OO + CH3 <=> CH3CH2O + CH3O""",
)

entry(
    index = 500,
    label = "CH3CH2OO + CH4 <=> CH3CH2OOH + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(12, 'cm^3/(mol*s)'), n=3.69, Ea=(21300, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OO + CH4 <=> CH3CH2OOH + CH3""",
)

entry(
    index = 501,
    label = "CH3CH2OO + CH3OH <=> CH3CH2OOH + CH2OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(19400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OO + CH3OH <=> CH3CH2OOH + CH2OH""",
)

entry(
    index = 502,
    label = "CH3CH2OO + CH2O <=> CH3CH2OOH + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(41000, 'cm^3/(mol*s)'), n=2.5, Ea=(10206, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OO + CH2O <=> CH3CH2OOH + HCO""",
)

entry(
    index = 503,
    label = "CH3CH2OO + CH2O <=> CH3CH2OOH + H + CO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.5e+14, 'cm^3/(mol*s)'),
        n = 0.027,
        Ea = (30133, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3CH2OO + CH2O <=> CH3CH2OOH + H + CO""",
)

entry(
    index = 504,
    label = "CH3CH2OO + C2H6 <=> CH3CH2OOH + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.6, 'cm^3/(mol*s)'), n=3.76, Ea=(17200, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OO + C2H6 <=> CH3CH2OOH + C2H5""",
)

entry(
    index = 505,
    label = "CH3CH2OO + C2H5 <=> CH3CH2O + CH3CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.1e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1411, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OO + C2H5 <=> CH3CH2O + CH3CH2O""",
)

entry(
    index = 506,
    label = "CH3CH2OO + CH3CHO <=> CH3CH2OOH + CH3CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+13, 'cm^3/(mol*s)'), n=0, Ea=(16293, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OO + CH3CHO <=> CH3CH2OOH + CH3CO""",
)

entry(
    index = 507,
    label = "CH3CH2OO + CH3CHO <=> CH3CH2OOH + CH2CHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.1e+13, 'cm^3/(mol*s)'), n=0, Ea=(23248, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OO + CH3CHO <=> CH3CH2OOH + CH2CHO""",
)

entry(
    index = 508,
    label = "CH3CH2OO + CH3CH2OO <=> CH3CH2O + CH3CH2O + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.2e+10, 'cm^3/(mol*s)'), n=0, Ea=(46, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OO + CH3CH2OO <=> CH3CH2O + CH3CH2O + O2""",
)

entry(
    index = 509,
    label = "CH3CH2OO + CH3CH2OO <=> CH3CHO + CH3CH2OH + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.6e+10, 'cm^3/(mol*s)'), n=0, Ea=(46, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OO + CH3CH2OO <=> CH3CHO + CH3CH2OH + O2""",
)

entry(
    index = 510,
    label = "CH2CH2OOH <=> cC2H4O + OH",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100], 'atm'),
        arrhenius = [
            Arrhenius(A=(2.2e+24, 's^-1'), n=-5.76, Ea=(12410, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(7.3e+26, 's^-1'), n=-6.39, Ea=(13340, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.3e+29, 's^-1'), n=-6.91, Ea=(14240, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.7e+28, 's^-1'), n=-6.45, Ea=(14230, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2.2e+30, 's^-1'), n=-6.94, Ea=(15220, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.4e+30, 's^-1'), n=-6.7, Ea=(15540, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1e+32, 's^-1'), n=-7.1, Ea=(16610, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(6.3e+31, 's^-1'), n=-6.87, Ea=(17080, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2e+31, 's^-1'), n=-6.53, Ea=(17550, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.2e+30, 's^-1'), n=-6, Ea=(17750, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(4e+27, 's^-1'), n=-5.08, Ea=(17550, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(7.8e+24, 's^-1'), n=-4.12, Ea=(17130, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(3.2e+21, 's^-1'), n=-2.97, Ea=(16400, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CH2OOH <=> cC2H4O + OH""",
)

entry(
    index = 511,
    label = "CH2CHOOH <=> CH2CHO + OH",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([1, 10, 50, 100], 'atm'),
        arrhenius = [
            Arrhenius(A=(2e+35, 's^-1'), n=-6.7, Ea=(47450, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.1e+28, 's^-1'), n=-4.15, Ea=(46190, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2.8e+26, 's^-1'), n=-3.5, Ea=(46340, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2.2e+17, 's^-1'), n=-0.42, Ea=(44622, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHOOH <=> CH2CHO + OH""",
)

entry(
    index = 512,
    label = "CH2CHOOH + H <=> CH2CHOO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.3e+10, 'cm^3/(mol*s)'), n=0, Ea=(1860, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHOOH + H <=> CH2CHOO + H2""",
)

entry(
    index = 513,
    label = "CH2CHOOH + H <=> CH2CHO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+10, 'cm^3/(mol*s)'), n=0, Ea=(1860, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHOOH + H <=> CH2CHO + H2O""",
)

entry(
    index = 514,
    label = "CH2CHOOH + O <=> CH2CHOO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.7e+12, 'cm^3/(mol*s)'), n=0, Ea=(4750, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHOOH + O <=> CH2CHOO + OH""",
)

entry(
    index = 515,
    label = "CH2CHOOH + OH <=> CH2CHOO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.1e+12, 'cm^3/(mol*s)'), n=0, Ea=(-437, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHOOH + OH <=> CH2CHOO + H2O""",
)

entry(
    index = 516,
    label = "CH2CHOOH + HO2 <=> CH2CHOO + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(41000, 'cm^3/(mol*s)'), n=2.5, Ea=(10206, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHOOH + HO2 <=> CH2CHOO + H2O2""",
)

entry(
    index = 517,
    label = "CH2CHOO <=> CHCHO + OH",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiPDepArrhenius(
        arrhenius = [
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(A=(3.6e+49, 's^-1'), n=-12.13, Ea=(67420, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.4e+36, 's^-1'), n=-9.92, Ea=(41220, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(4.2e+40, 's^-1'), n=-10.53, Ea=(43670, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(3.8e+46, 's^-1'), n=-10.72, Ea=(51900, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.6e+49, 's^-1'), n=-11.24, Ea=(54150, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(2.4e+51, 's^-1'), n=-11.64, Ea=(56980, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(2e+54, 's^-1'), n=-12.22, Ea=(61840, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(9.5e+195, 's^-1'), n=-52.27, Ea=(163500, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(A=(1.2e+56, 's^-1'), n=-14.81, Ea=(60700, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(2.3e+40, 's^-1'), n=-9.39, Ea=(50420, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.6e+43, 's^-1'), n=-9.99, Ea=(50290, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(2.3e+124, 's^-1'), n=-36.77, Ea=(70100, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.9e+103, 's^-1'), n=-29.49, Ea=(65410, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(6e+86, 's^-1'), n=-23.81, Ea=(62170, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.5e+57, 's^-1'), n=-13.94, Ea=(55390, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.8e+34, 's^-1'), n=-6.4, Ea=(50000, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHOO <=> CHCHO + OH""",
)

entry(
    index = 518,
    label = "CH2CHOO <=> CH2CHO + O",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiPDepArrhenius(
        arrhenius = [
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(A=(2.7e+180, 's^-1'), n=-48.19, Ea=(169300, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(3.9e+38, 's^-1'), n=-8.69, Ea=(42770, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(4.6e+47, 's^-1'), n=-11.21, Ea=(47050, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(7.6e+81, 's^-1'), n=-21.28, Ea=(65080, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.9e+68, 's^-1'), n=-16.83, Ea=(60680, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(2e+55, 's^-1'), n=-12.69, Ea=(55840, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.1e+53, 's^-1'), n=-11.79, Ea=(56690, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(4.3e+48, 's^-1'), n=-10.31, Ea=(56090, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(A=(1.5e+30, 's^-1'), n=-6.64, Ea=(41110, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(9.7e-12, 's^-1'), n=5.96, Ea=(22890, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(4e+22, 's^-1'), n=-3.71, Ea=(36270, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(2.4e+33, 's^-1'), n=-6.62, Ea=(41280, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(6.4e+31, 's^-1'), n=-5.96, Ea=(41260, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(2.1e+29, 's^-1'), n=-5.1, Ea=(40710, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(4.7e+27, 's^-1'), n=-4.5, Ea=(40530, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(6e+25, 's^-1'), n=-3.85, Ea=(40120, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHOO <=> CH2CHO + O""",
)

entry(
    index = 519,
    label = "CH2CHOO <=> OCHCHO + H",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiPDepArrhenius(
        arrhenius = [
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(A=(6.4e+80, 's^-1'), n=-22.2, Ea=(51750, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(3.3e+65, 's^-1'), n=-17.01, Ea=(48090, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(6e+51, 's^-1'), n=-12.62, Ea=(43000, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.5e+44, 's^-1'), n=-10.12, Ea=(40790, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.3e+59, 's^-1'), n=-14.33, Ea=(51390, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(4.9e+26, 's^-1'), n=-4.67, Ea=(34320, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(2.1e+33, 's^-1'), n=-6.38, Ea=(39520, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.3e+32, 's^-1'), n=-5.92, Ea=(40660, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(A=(1.2e+28, 's^-1'), n=-6.01, Ea=(28740, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.4e+25, 's^-1'), n=-4.8, Ea=(28940, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(2.9e+20, 's^-1'), n=-3.29, Ea=(27550, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.6e+19, 's^-1'), n=-2.82, Ea=(27620, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.9e+22, 's^-1'), n=-3.54, Ea=(29980, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(7.5e+29, 's^-1'), n=-5.75, Ea=(34490, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(7.1e+61, 's^-1'), n=-16.16, Ea=(43280, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.1e+19, 's^-1'), n=-2.56, Ea=(29670, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHOO <=> OCHCHO + H""",
)

entry(
    index = 520,
    label = "CH2CHOO <=> CH2CO + OH",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiPDepArrhenius(
        arrhenius = [
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(A=(1.2e+47, 's^-1'), n=-12.28, Ea=(75330, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(8.4e+09, 's^-1'), n=-2.06, Ea=(33720, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(61000, 's^-1'), n=0.17, Ea=(34220, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.5e+19, 's^-1'), n=-3.61, Ea=(43060, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(2.1e+33, 's^-1'), n=-7.39, Ea=(51610, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(4.4e+36, 's^-1'), n=-7.99, Ea=(54680, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.2e+37, 's^-1'), n=-7.8, Ea=(56460, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(9.1e+35, 's^-1'), n=-7.21, Ea=(57550, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(A=(230, 's^-1'), n=-0.73, Ea=(25710, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.8e-23, 's^-1'), n=7.84, Ea=(20190, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(3.8e+63, 's^-1'), n=-20.44, Ea=(43420, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(3.2e+27, 's^-1'), n=-7.76, Ea=(37230, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(2.3e-05, 's^-1'), n=3.47, Ea=(31560, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(0.11, 's^-1'), n=2.64, Ea=(34160, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(560, 's^-1'), n=1.7, Ea=(36450, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.1e+07, 's^-1'), n=0.52, Ea=(38670, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHOO <=> CH2CO + OH""",
)

entry(
    index = 521,
    label = "CH2CHOO <=> CH2O + HCO",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiPDepArrhenius(
        arrhenius = [
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(A=(1.7e+174, 's^-1'), n=-55.52, Ea=(60320, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(9e+66, 's^-1'), n=-17.25, Ea=(48120, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.8e+43, 's^-1'), n=-9.87, Ea=(37960, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(8.6e+33, 's^-1'), n=-6.88, Ea=(34370, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(7.3e+171, 's^-1'), n=-43.53, Ea=(191900, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1e+32, 's^-1'), n=-6.06, Ea=(35500, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.9e+34, 's^-1'), n=-6.57, Ea=(38510, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(5.7e+29, 's^-1'), n=-5.19, Ea=(36800, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(A=(2.3e+35, 's^-1'), n=-7.97, Ea=(31280, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(2.1e+26, 's^-1'), n=-4.96, Ea=(28780, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.5e+20, 's^-1'), n=-3.08, Ea=(26630, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.1e+130, 's^-1'), n=-39.38, Ea=(54700, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(2.4e+34, 's^-1'), n=-6.87, Ea=(35700, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(2.2e+175, 's^-1'), n=-53.78, Ea=(68500, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.1e+185, 's^-1'), n=-54.22, Ea=(88990, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(470, 's^-1'), n=1.81, Ea=(18100, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHOO <=> CH2O + HCO""",
)

entry(
    index = 522,
    label = "CH2CHOO <=> CH2O + H + CO",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiPDepArrhenius(
        arrhenius = [
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(A=(3.9e+174, 's^-1'), n=-55.52, Ea=(60320, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(2.1e+67, 's^-1'), n=-17.25, Ea=(48120, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(4.3e+43, 's^-1'), n=-9.87, Ea=(37960, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(2e+34, 's^-1'), n=-6.88, Ea=(34370, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.7e+172, 's^-1'), n=-43.53, Ea=(191900, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(2.4e+32, 's^-1'), n=-6.06, Ea=(35500, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(4.3e+34, 's^-1'), n=-6.57, Ea=(38510, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.3e+30, 's^-1'), n=-5.19, Ea=(36800, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(A=(5.3e+35, 's^-1'), n=-7.97, Ea=(31280, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(4.9e+26, 's^-1'), n=-4.96, Ea=(28780, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(3.4e+20, 's^-1'), n=-3.08, Ea=(26630, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(2.5e+130, 's^-1'), n=-39.38, Ea=(54700, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(5.5e+34, 's^-1'), n=-6.87, Ea=(35700, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(5.1e+175, 's^-1'), n=-53.78, Ea=(68500, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(2.5e+185, 's^-1'), n=-54.22, Ea=(88990, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1100, 's^-1'), n=1.81, Ea=(18100, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHOO <=> CH2O + H + CO""",
)

entry(
    index = 523,
    label = "CH2CHOO <=> CO + CH3O",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiPDepArrhenius(
        arrhenius = [
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(A=(5.2e+33, 's^-1'), n=-7.92, Ea=(31320, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.3e+98, 's^-1'), n=-27.09, Ea=(64060, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.8e+33, 's^-1'), n=-7.27, Ea=(33760, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(3.8e+33, 's^-1'), n=-7.2, Ea=(35100, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.3e+79, 's^-1'), n=-19.61, Ea=(74870, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(4.1e+32, 's^-1'), n=-6.62, Ea=(37210, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(6.9e+44, 's^-1'), n=-10.04, Ea=(47030, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(29, 's^-1'), n=2.492, Ea=(21710, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(A=(2.3e+129, 's^-1'), n=-41.86, Ea=(45850, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(2.4e+28, 's^-1'), n=-5.99, Ea=(30540, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(8.7e-50, 's^-1'), n=16.63, Ea=(-3900, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.2e-39, 's^-1'), n=13.61, Ea=(-1317, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(8.8e+86, 's^-1'), n=-23.08, Ea=(61010, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1300, 's^-1'), n=1.44, Ea=(18660, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(2e+17, 's^-1'), n=-2.23, Ea=(28590, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(29, 's^-1'), n=2.492, Ea=(21710, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHOO <=> CO + CH3O""",
)

entry(
    index = 524,
    label = "CH2CHOO <=> CO2 + CH3",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiPDepArrhenius(
        arrhenius = [
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(A=(5.1e+33, 's^-1'), n=-7.95, Ea=(31290, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.2e+118, 's^-1'), n=-33.13, Ea=(73790, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(8.6e+32, 's^-1'), n=-7.21, Ea=(33550, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(3.3e+33, 's^-1'), n=-7.22, Ea=(34990, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(3.5e-79, 's^-1'), n=25.01, Ea=(-21020, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(8.2e+32, 's^-1'), n=-6.76, Ea=(37270, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(7e+37, 's^-1'), n=-8.06, Ea=(42200, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(6.7e+48, 's^-1'), n=-11.657, Ea=(44080, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            PDepArrhenius(
                pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
                arrhenius = [
                    Arrhenius(A=(4.2e+122, 's^-1'), n=-39.75, Ea=(43640, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(2e+29, 's^-1'), n=-6.29, Ea=(30920, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(5.1e-66, 's^-1'), n=21.37, Ea=(-11110, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.8e-47, 's^-1'), n=15.85, Ea=(-5283, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(3.8e+32, 's^-1'), n=-6.8, Ea=(35690, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(4.6, 's^-1'), n=2.1, Ea=(17170, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(3.5e+14, 's^-1'), n=-1.58, Ea=(26470, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1900, 's^-1'), n=2.081, Ea=(25118, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHOO <=> CO2 + CH3""",
)

entry(
    index = 525,
    label = "CH2CHOO + H <=> CH2CHO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHOO + H <=> CH2CHO + OH""",
)

entry(
    index = 526,
    label = "CH2CHOO + O <=> CH2CHO + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(-145, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHOO + O <=> CH2CHO + O2""",
)

entry(
    index = 527,
    label = "CH2CHOO + OH <=> CH2CHOH + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+15, 'cm^3/(mol*s)'), n=-0.6, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHOO + OH <=> CH2CHOH + O2""",
)

entry(
    index = 528,
    label = "CH2CHOO + OH <=> CH2CHO + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+11, 'cm^3/(mol*s)'), n=0.6, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHOO + OH <=> CH2CHO + HO2""",
)

entry(
    index = 529,
    label = "CH2CHOO + HO2 <=> CH2CHOOH + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.5e+11, 'cm^3/(mol*s)'), n=0, Ea=(-1391, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHOO + HO2 <=> CH2CHOOH + O2""",
)

entry(
    index = 530,
    label = "CH2CHOO + CO <=> CH2CHO + CO2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (160000, 'cm^3/(mol*s)'),
        n = 2.18,
        Ea = (17940, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHOO + CO <=> CH2CHO + CO2""",
)

entry(
    index = 531,
    label = "CH2CHOO + CH3 <=> CH2CHO + CH3O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.1e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1411, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHOO + CH3 <=> CH2CHO + CH3O""",
)

entry(
    index = 532,
    label = "CH2CHOO + CH4 <=> CH2CHOOH + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(47000, 'cm^3/(mol*s)'), n=2.5, Ea=(21000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHOO + CH4 <=> CH2CHOOH + CH3""",
)

entry(
    index = 533,
    label = "CH2CHOO + CH3OH <=> CH2CHOOH + CH2OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(19400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHOO + CH3OH <=> CH2CHOOH + CH2OH""",
)

entry(
    index = 534,
    label = "CH2CHOO + CH2O <=> CH2CHOOH + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(41000, 'cm^3/(mol*s)'), n=2.5, Ea=(10206, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHOO + CH2O <=> CH2CHOOH + HCO""",
)

entry(
    index = 535,
    label = "CH2CHOO + C2H6 <=> CH2CHOOH + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.6, 'cm^3/(mol*s)'), n=3.76, Ea=(17200, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHOO + C2H6 <=> CH2CHOOH + C2H5""",
)

entry(
    index = 536,
    label = "OCHCHO <=> CH2O + CO",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.009869, 0.04935, 0.09869, 0.4935, 0.9869, 4.935, 9.869], 'atm'),
        arrhenius = [
            Arrhenius(A=(4.2e+53, 's^-1'), n=-12.5, Ea=(70845, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(5.1e+54, 's^-1'), n=-12.6, Ea=(73012, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1e+55, 's^-1'), n=-12.6, Ea=(73877, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(4.5e+55, 's^-1'), n=-12.6, Ea=(75869, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(8e+55, 's^-1'), n=-12.6, Ea=(76713, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.1e+56, 's^-1'), n=-12.2, Ea=(77643, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(5.5e+56, 's^-1'), n=-12.6, Ea=(79964, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is OCHCHO <=> CH2O + CO""",
)

entry(
    index = 537,
    label = "OCHCHO <=> HCOH + CO",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.009869, 0.04935, 0.09869, 0.4935, 0.9869, 4.935, 9.869], 'atm'),
        arrhenius = [
            Arrhenius(A=(8.4e+52, 's^-1'), n=-12.6, Ea=(72393, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(8.3e+54, 's^-1'), n=-12.9, Ea=(75113, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(4.4e+55, 's^-1'), n=-13, Ea=(76257, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.3e+57, 's^-1'), n=-13.2, Ea=(78851, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2.6e+57, 's^-1'), n=-13.2, Ea=(79754, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1e+57, 's^-1'), n=-12.9, Ea=(81161, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(5.7e+59, 's^-1'), n=-13.3, Ea=(83539, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is OCHCHO <=> HCOH + CO""",
)

entry(
    index = 538,
    label = "OCHCHO <=> CO + CO + H2",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.009869, 0.04935, 0.09869, 0.4935, 0.9869, 4.935, 9.869], 'atm'),
        arrhenius = [
            Arrhenius(A=(6e+51, 's^-1'), n=-12.1, Ea=(71854, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.4e+54, 's^-1'), n=-12.5, Ea=(74751, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.8e+55, 's^-1'), n=-12.7, Ea=(76137, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.3e+57, 's^-1'), n=-13, Ea=(78972, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(6.1e+57, 's^-1'), n=-13.1, Ea=(80147, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(5.8e+57, 's^-1'), n=-12.9, Ea=(81871, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(3.4e+59, 's^-1'), n=-13.3, Ea=(84294, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is OCHCHO <=> CO + CO + H2""",
)

entry(
    index = 539,
    label = "OCHCHO <=> HCO + HCO",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.009869, 0.04935, 0.09869, 0.4935, 0.9869, 4.935, 9.869], 'atm'),
        arrhenius = [
            Arrhenius(A=(1e+42, 's^-1'), n=-9.7, Ea=(73534, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(6e+48, 's^-1'), n=-11.1, Ea=(77462, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.7e+51, 's^-1'), n=-11.6, Ea=(79111, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(5.3e+55, 's^-1'), n=-12.5, Ea=(82774, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.9e+57, 's^-1'), n=-12.8, Ea=(84321, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2.2e+59, 's^-1'), n=-13.1, Ea=(87258, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(3e+60, 's^-1'), n=-13.3, Ea=(88993, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is OCHCHO <=> HCO + HCO""",
)

entry(
    index = 540,
    label = "OCHCHO + H <=> CH2O + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(4300, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is OCHCHO + H <=> CH2O + HCO""",
)

entry(
    index = 541,
    label = "OCHCHO + H <=> CH2O + H + CO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.2e+23, 'cm^3/(mol*s)'),
        n = -2.473,
        Ea = (24227, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is OCHCHO + H <=> CH2O + H + CO""",
)

entry(
    index = 542,
    label = "OCHCHO + O <=> OCHCO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.2e+11, 'cm^3/(mol*s)'),
        n = 0.57,
        Ea = (2760, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is OCHCHO + O <=> OCHCO + OH""",
)

entry(
    index = 543,
    label = "OCHCHO + OH <=> OCHCO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+06, 'cm^3/(mol*s)'), n=2, Ea=(-1630, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is OCHCHO + OH <=> OCHCO + H2O""",
)

entry(
    index = 544,
    label = "OCHCHO + HO2 <=> HOCHO + CO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (0.00033, 'cm^3/(mol*s)'),
        n = 3.995,
        Ea = (300, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is OCHCHO + HO2 <=> HOCHO + CO + OH""",
)

entry(
    index = 545,
    label = "OCHCHO + HO2 <=> OCHCO + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(41000, 'cm^3/(mol*s)'), n=2.5, Ea=(10206, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is OCHCHO + HO2 <=> OCHCO + H2O2""",
)

entry(
    index = 546,
    label = "OCHCHO + O2 <=> OCHCO + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(240000, 'cm^3/(mol*s)'), n=2.5, Ea=(36461, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is OCHCHO + O2 <=> OCHCO + HO2""",
)

entry(
    index = 547,
    label = "OCHCO <=> HCO + CO",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.1, 1, 100], 'atm'),
        arrhenius = [
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(3.8e+12, 's^-1'), n=0, Ea=(8610, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(-8e+21, 's^-1'), n=-2.359, Ea=(27420, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            Arrhenius(A=(3.8e+13, 's^-1'), n=0, Ea=(8665, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(4.1e+14, 's^-1'), n=0, Ea=(8765, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.1e+14, 's^-1'), n=0.133, Ea=(10140, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is OCHCO <=> HCO + CO""",
)

entry(
    index = 548,
    label = "OCHCO <=> H + CO + CO",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.1, 1, 100], 'atm'),
        arrhenius = [
            Arrhenius(A=(8e+21, 's^-1'), n=-2.359, Ea=(27420, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2.3e+23, 's^-1'), n=-2.473, Ea=(28592, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2.5e+24, 's^-1'), n=-2.473, Ea=(28692, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.4e+24, 's^-1'), n=-2.419, Ea=(30991, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is OCHCO <=> H + CO + CO""",
)

entry(
    index = 549,
    label = "OCHCO + O2 <=> CO + CO2 + OH",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.1, 1], 'atm'),
        arrhenius = [
            Arrhenius(A=(1.6e+14, 'cm^3/(mol*s)'), n=0, Ea=(1540, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.1e+14, 'cm^3/(mol*s)'), n=0, Ea=(1300, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(3.3e+14, 'cm^3/(mol*s)'), n=0, Ea=(2075, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is OCHCO + O2 <=> CO + CO2 + OH""",
)

entry(
    index = 550,
    label = "HOCHO <=> CO + H2O",
    degeneracy = 1,
    kinetics = Lindemann(
        arrheniusHigh = Arrhenius(A=(7.5e+14, 's^-1'), n=0, Ea=(68710, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(A=(4.1e+15, 'cm^3/(mol*s)'), n=0, Ea=(52980, 'cal/mol'), T0=(1, 'K')),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is HOCHO <=> CO + H2O""",
)

entry(
    index = 551,
    label = "HOCHO <=> CO2 + H2",
    degeneracy = 1,
    kinetics = Lindemann(
        arrheniusHigh = Arrhenius(A=(4.5e+13, 's^-1'), n=0, Ea=(68240, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(A=(1.7e+15, 'cm^3/(mol*s)'), n=0, Ea=(51110, 'cal/mol'), T0=(1, 'K')),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is HOCHO <=> CO2 + H2""",
)

entry(
    index = 552,
    label = "HOCHO + H <=> HOCO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(230, 'cm^3/(mol*s)'), n=3.272, Ea=(4858, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HOCHO + H <=> HOCO + H2""",
)

entry(
    index = 553,
    label = "HOCHO + H <=> OCHO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (420000, 'cm^3/(mol*s)'),
        n = 2.255,
        Ea = (14091, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is HOCHO + H <=> OCHO + H2""",
)

entry(
    index = 554,
    label = "HOCHO + O <=> HOCO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(51, 'cm^3/(mol*s)'), n=3.422, Ea=(4216, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HOCHO + O <=> HOCO + OH""",
)

entry(
    index = 555,
    label = "HOCHO + O <=> OCHO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (170000, 'cm^3/(mol*s)'),
        n = 2.103,
        Ea = (9880, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is HOCHO + O <=> OCHO + OH""",
)

entry(
    index = 556,
    label = "HOCHO + OH <=> HOCO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (7.8e-06, 'cm^3/(mol*s)'),
        n = 5.57,
        Ea = (-2365, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is HOCHO + OH <=> HOCO + H2O""",
)

entry(
    index = 557,
    label = "HOCHO + OH <=> OCHO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.9e-05, 'cm^3/(mol*s)'),
        n = 4.91,
        Ea = (-5067, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is HOCHO + OH <=> OCHO + H2O""",
)

entry(
    index = 558,
    label = "HOCHO + HO2 <=> HOCO + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.47, 'cm^3/(mol*s)'), n=3.975, Ea=(16787, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HOCHO + HO2 <=> HOCO + H2O2""",
)

entry(
    index = 559,
    label = "HOCHO + HO2 <=> OCHO + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(39, 'cm^3/(mol*s)'), n=3.08, Ea=(25206, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HOCHO + HO2 <=> OCHO + H2O2""",
)

entry(
    index = 560,
    label = "HOCO + HO2 <=> HOCHO + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HOCO + HO2 <=> HOCHO + O2""",
)

entry(
    index = 561,
    label = "HOCHO + O2 <=> OCHO + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(63000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HOCHO + O2 <=> OCHO + HO2""",
)

entry(
    index = 562,
    label = "OCHO <=> CO2 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+10, 's^-1'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is OCHO <=> CO2 + H""",
)

entry(
    index = 563,
    label = "OCHO + O2 <=> CO2 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is OCHO + O2 <=> CO2 + HO2""",
)

entry(
    index = 564,
    label = "HOCH2CH2OO <=> CH2O + CH2O + OH",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.013, 0.1, 1, 10, 100], 'atm'),
        arrhenius = [
            Arrhenius(A=(5.6e+13, 's^-1'), n=-1.9, Ea=(14338, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2.2e+14, 's^-1'), n=-1.92, Ea=(14870, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2e+15, 's^-1'), n=-2.03, Ea=(15913, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(5e+16, 's^-1'), n=-2.26, Ea=(17552, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(3.8e+18, 's^-1'), n=-2.6, Ea=(19972, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is HOCH2CH2OO <=> CH2O + CH2O + OH""",
)

entry(
    index = 565,
    label = "HOCH2CH2OO <=> CH2CHOH + HO2",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.013, 0.1, 1, 10, 100], 'atm'),
        arrhenius = [
            Arrhenius(A=(2e+09, 's^-1'), n=-1.01, Ea=(13160, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2.1e+09, 's^-1'), n=-0.81, Ea=(13598, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.4e+10, 's^-1'), n=-0.78, Ea=(14836, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(8.5e+11, 's^-1'), n=-1.01, Ea=(17045, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(5.8e+14, 's^-1'), n=-1.51, Ea=(20561, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is HOCH2CH2OO <=> CH2CHOH + HO2""",
)

entry(
    index = 566,
    label = "HOCH2CH2OO + HO2 => CH2OOH + CH2OH + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+11, 'cm^3/(mol*s)'), n=0, Ea=(-1490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HOCH2CH2OO + HO2 => CH2OOH + CH2OH + O2""",
)

entry(
    index = 567,
    label = "HOCH2CH2OO + CH2O => CH2OOH + CH2OH + HCO",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(41000, 'cm^3/(mol*s)'), n=2.5, Ea=(10206, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HOCH2CH2OO + CH2O => CH2OOH + CH2OH + HCO""",
)

entry(
    index = 568,
    label = "HOCH2CH2OO + C2H4 => CH2O + CH2OH + CH3CHO",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.2e+12, 'cm^3/(mol*s)'), n=0, Ea=(17200, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HOCH2CH2OO + C2H4 => CH2O + CH2OH + CH3CHO""",
)

entry(
    index = 569,
    label = "NH2 + H <=> NH3",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.6e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.6e+22, 'cm^6/(mol^2*s)'),
            n = -1.76,
            Ea = (0, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.5,
        T3 = (1e-30, 'K'),
        T1 = (1e+30, 'K'),
        T2 = (1e+30, 'K'),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is NH2 + H <=> NH3""",
)

entry(
    index = 570,
    label = "NH3 + H <=> NH2 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (640000, 'cm^3/(mol*s)'),
        n = 2.39,
        Ea = (10171, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NH3 + H <=> NH2 + H2""",
)

entry(
    index = 571,
    label = "NH3 + O <=> NH2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (9.4e+06, 'cm^3/(mol*s)'),
        n = 1.94,
        Ea = (6460, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NH3 + O <=> NH2 + OH""",
)

entry(
    index = 572,
    label = "NH3 + OH <=> NH2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+06, 'cm^3/(mol*s)'), n=2.04, Ea=(566, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NH3 + OH <=> NH2 + H2O""",
)

entry(
    index = 573,
    label = "NH3 + HO2 <=> NH2 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+11, 'cm^3/(mol*s)'), n=0, Ea=(22000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NH3 + HO2 <=> NH2 + H2O2""",
)

entry(
    index = 574,
    label = "NH + H2 <=> NH2 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.1e+13, 'cm^3/(mol*s)'), n=0, Ea=(15417, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NH + H2 <=> NH2 + H""",
)

entry(
    index = 575,
    label = "NH2 + O <=> HNO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.6e+14, 'cm^3/(mol*s)'), n=-0.5, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NH2 + O <=> HNO + H""",
)

entry(
    index = 576,
    label = "NH2 + O <=> NH + OH",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(0.86, 'cm^3/(mol*s)'), n=4.01, Ea=(1673, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is NH2 + O <=> NH + OH""",
)

entry(
    index = 577,
    label = "NH2 + OH <=> NH + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.3e+06, 'cm^3/(mol*s)'),
        n = 1.949,
        Ea = (-217, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NH2 + OH <=> NH + H2O""",
)

entry(
    index = 578,
    label = "NH2 + HO2 <=> NH3 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(17000, 'cm^3/(mol*s)'), n=1.55, Ea=(2027, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NH2 + HO2 <=> NH3 + O2""",
)

entry(
    index = 579,
    label = "NH2 + HO2 <=> H2NO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.5e+17, 'cm^3/(mol*s)'),
        n = -1.28,
        Ea = (1166, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NH2 + HO2 <=> H2NO + OH""",
)

entry(
    index = 580,
    label = "NH2 + HO2 <=> HNO + H2O",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(1.6e+07, 'cm^3/(mol*s)'), n=0.55, Ea=(525, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (5.7e+15, 'cm^3/(mol*s)'),
                n = -1.12,
                Ea = (707, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is NH2 + HO2 <=> HNO + H2O""",
)

entry(
    index = 581,
    label = "NH2 + HO2 <=> HON + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.1e+07, 'cm^3/(mol*s)'), n=0.64, Ea=(811, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NH2 + HO2 <=> HON + H2O""",
)

entry(
    index = 582,
    label = "NH2 + O2 <=> H2NO + O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.6e+11, 'cm^3/(mol*s)'),
        n = 0.4872,
        Ea = (29050, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NH2 + O2 <=> H2NO + O""",
)

entry(
    index = 583,
    label = "NH2 + O2 <=> HNO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (0.029, 'cm^3/(mol*s)'),
        n = 3.764,
        Ea = (18185, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NH2 + O2 <=> HNO + OH""",
)

entry(
    index = 584,
    label = "NH2 + NH2 <=> NH3 + NH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.6, 'cm^3/(mol*s)'), n=3.53, Ea=(552, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NH2 + NH2 <=> NH3 + NH""",
)

entry(
    index = 585,
    label = "NH2 + NH <=> NH3 + N",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9600, 'cm^3/(mol*s)'), n=2.46, Ea=(107, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NH2 + NH <=> NH3 + N""",
)

entry(
    index = 586,
    label = "NH2 + N <=> N2 + H + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NH2 + N <=> N2 + H + H""",
)

entry(
    index = 587,
    label = "NH2 + HNO <=> NH3 + NO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(590, 'cm^3/(mol*s)'), n=2.95, Ea=(-3469, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NH2 + HNO <=> NH3 + NO""",
)

entry(
    index = 588,
    label = "NH2 + NO <=> N2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.6e+19, 'cm^3/(mol*s)'),
        n = -2.369,
        Ea = (870, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NH2 + NO <=> N2 + H2O""",
)

entry(
    index = 589,
    label = "NH2 + NO <=> NNH + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.3e+10, 'cm^3/(mol*s)'),
        n = 0.294,
        Ea = (-866, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NH2 + NO <=> NNH + OH""",
)

entry(
    index = 590,
    label = "NH2 + HONO <=> NH3 + NO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(71, 'cm^3/(mol*s)'), n=3.02, Ea=(-4940, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NH2 + HONO <=> NH3 + NO2""",
)

entry(
    index = 591,
    label = "NH2 + NO2 <=> H2NO + NO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8.6e+11, 'cm^3/(mol*s)'),
        n = 0.11,
        Ea = (-1186, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NH2 + NO2 <=> H2NO + NO""",
)

entry(
    index = 592,
    label = "NH2 + NO2 <=> N2O + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.2e+11, 'cm^3/(mol*s)'),
        n = 0.11,
        Ea = (-1186, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NH2 + NO2 <=> N2O + H2O""",
)

entry(
    index = 593,
    label = "NH + H <=> N + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NH + H <=> N + H2""",
)

entry(
    index = 594,
    label = "NH + O <=> NO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NH + O <=> NO + H""",
)

entry(
    index = 595,
    label = "NH + OH <=> HNO + H",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.2e+14, 'cm^3/(mol*s)'),
        n = -0.376,
        Ea = (-46, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NH + OH <=> HNO + H""",
)

entry(
    index = 596,
    label = "NH + OH <=> N + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.6e+07, 'cm^3/(mol*s)'),
        n = 1.733,
        Ea = (-576, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NH + OH <=> N + H2O""",
)

entry(
    index = 597,
    label = "NH + O2 <=> HNO + O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(13850, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NH + O2 <=> HNO + O""",
)

entry(
    index = 598,
    label = "NH + O2 <=> NO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.9e+10, 'cm^3/(mol*s)'), n=0, Ea=(1530, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NH + O2 <=> NO + OH""",
)

entry(
    index = 599,
    label = "NH + NH <=> NH2 + N",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.57, 'cm^3/(mol*s)'), n=3.88, Ea=(342, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NH + NH <=> NH2 + N""",
)

entry(
    index = 600,
    label = "NH + N <=> N2 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NH + N <=> N2 + H""",
)

entry(
    index = 601,
    label = "NH + NO <=> N2O + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.7e+15, 'cm^3/(mol*s)'), n=-0.78, Ea=(20, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NH + NO <=> N2O + H""",
)

entry(
    index = 602,
    label = "NH + NO <=> N2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.8e+14, 'cm^3/(mol*s)'), n=-0.78, Ea=(20, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NH + NO <=> N2 + OH""",
)

entry(
    index = 603,
    label = "NH + HONO <=> NH2 + NO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NH + HONO <=> NH2 + NO2""",
)

entry(
    index = 604,
    label = "NH + NO2 <=> N2O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.1e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NH + NO2 <=> N2O + OH""",
)

entry(
    index = 605,
    label = "NH + NO2 <=> HNO + NO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.9e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NH + NO2 <=> HNO + NO""",
)

entry(
    index = 606,
    label = "N + OH <=> NO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is N + OH <=> NO + H""",
)

entry(
    index = 607,
    label = "N + O2 <=> NO + O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.4e+09, 'cm^3/(mol*s)'), n=1, Ea=(6280, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is N + O2 <=> NO + O""",
)

entry(
    index = 608,
    label = "N + NO <=> N2 + O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.4e+12, 'cm^3/(mol*s)'), n=0.14, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is N + NO <=> N2 + O""",
)

entry(
    index = 609,
    label = "NNH <=> N2 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+09, 's^-1'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NNH <=> N2 + H""",
)

entry(
    index = 610,
    label = "NNH + H <=> N2 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NNH + H <=> N2 + H2""",
)

entry(
    index = 611,
    label = "NNH + O <=> N2O + H",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.9e+14, 'cm^3/(mol*s)'),
        n = -0.274,
        Ea = (-22, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NNH + O <=> N2O + H""",
)

entry(
    index = 612,
    label = "NNH + O <=> N2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.2e+13, 'cm^3/(mol*s)'),
        n = 0.145,
        Ea = (-217, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NNH + O <=> N2 + OH""",
)

entry(
    index = 613,
    label = "NNH + O <=> NH + NO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.2e+11, 'cm^3/(mol*s)'),
        n = 0.381,
        Ea = (-409, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NNH + O <=> NH + NO""",
)

entry(
    index = 614,
    label = "NNH + OH <=> N2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NNH + OH <=> N2 + H2O""",
)

entry(
    index = 615,
    label = "NNH + O2 <=> N2 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.6e+14, 'cm^3/(mol*s)'),
        n = -0.385,
        Ea = (-13, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NNH + O2 <=> N2 + HO2""",
)

entry(
    index = 616,
    label = "NNH + NH <=> N2 + NH2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NNH + NH <=> N2 + NH2""",
)

entry(
    index = 617,
    label = "NNH + NH2 <=> N2 + NH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NNH + NH2 <=> N2 + NH3""",
)

entry(
    index = 618,
    label = "NNH + NO <=> N2 + HNO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NNH + NO <=> N2 + HNO""",
)

entry(
    index = 619,
    label = "NH2OH <=> NH2 + OH",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.4e+20, 's^-1'), n=-1.31, Ea=(64080, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.4e+37, 'cm^3/(mol*s)'),
            n = -5.96,
            Ea = (66783, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.31,
        T3 = (1e-30, 'K'),
        T1 = (1e+30, 'K'),
        T2 = (1e+30, 'K'),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is NH2OH <=> NH2 + OH""",
)

entry(
    index = 620,
    label = "NH2OH + H <=> HNOH + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.8e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(6249, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NH2OH + H <=> HNOH + H2""",
)

entry(
    index = 621,
    label = "NH2OH + H <=> H2NO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(5067, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NH2OH + H <=> H2NO + H2""",
)

entry(
    index = 622,
    label = "NH2OH + O <=> HNOH + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.3e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(3865, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NH2OH + O <=> HNOH + OH""",
)

entry(
    index = 623,
    label = "NH2OH + O <=> H2NO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(3010, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NH2OH + O <=> H2NO + OH""",
)

entry(
    index = 624,
    label = "NH2OH + OH <=> HNOH + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(15000, 'cm^3/(mol*s)'), n=2.61, Ea=(-3537, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NH2OH + OH <=> HNOH + H2O""",
)

entry(
    index = 625,
    label = "NH2OH + OH <=> H2NO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (150000, 'cm^3/(mol*s)'),
        n = 2.28,
        Ea = (-1296, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NH2OH + OH <=> H2NO + H2O""",
)

entry(
    index = 626,
    label = "NH2OH + NH2 <=> HNOH + NH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.11, 'cm^3/(mol*s)'), n=4, Ea=(-97, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NH2OH + NH2 <=> HNOH + NH3""",
)

entry(
    index = 627,
    label = "NH2OH + NH2 <=> H2NO + NH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.5, 'cm^3/(mol*s)'), n=3.42, Ea=(-1013, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NH2OH + NH2 <=> H2NO + NH3""",
)

entry(
    index = 628,
    label = "NH2OH + NH <=> HNOH + NH2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0029, 'cm^3/(mol*s)'), n=4.4, Ea=(1564, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NH2OH + NH <=> HNOH + NH2""",
)

entry(
    index = 629,
    label = "NH2OH + NH <=> H2NO + NH2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0015, 'cm^3/(mol*s)'), n=4.6, Ea=(2424, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NH2OH + NH <=> H2NO + NH2""",
)

entry(
    index = 630,
    label = "NH2OH + HO2 <=> HNOH + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(29000, 'cm^3/(mol*s)'), n=2.69, Ea=(9557, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NH2OH + HO2 <=> HNOH + H2O2""",
)

entry(
    index = 631,
    label = "NH2OH + HO2 <=> H2NO + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(14000, 'cm^3/(mol*s)'), n=2.69, Ea=(6418, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NH2OH + HO2 <=> H2NO + H2O2""",
)

entry(
    index = 632,
    label = "H2NO <=> HNO + H",
    degeneracy = 1,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(
            A = (2.8e+24, 'cm^3/(mol*s)'),
            n = -2.83,
            Ea = (64915, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        efficiencies = {'O': 10},
    ),
    shortDesc = u"""The chemkin file reaction is H2NO <=> HNO + H""",
)

entry(
    index = 633,
    label = "H2NO <=> HNOH",
    degeneracy = 1,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(A=(1.1e+29, 'cm^3/(mol*s)'), n=-4, Ea=(44000, 'cal/mol'), T0=(1, 'K')),
        efficiencies = {'O': 10},
    ),
    shortDesc = u"""The chemkin file reaction is H2NO <=> HNOH""",
)

entry(
    index = 634,
    label = "H2NO + H <=> HNO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+07, 'cm^3/(mol*s)'), n=2, Ea=(2000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2NO + H <=> HNO + H2""",
)

entry(
    index = 635,
    label = "H2NO + H <=> NH2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2NO + H <=> NH2 + OH""",
)

entry(
    index = 636,
    label = "H2NO + O <=> HNO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+07, 'cm^3/(mol*s)'), n=2, Ea=(2000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2NO + O <=> HNO + OH""",
)

entry(
    index = 637,
    label = "H2NO + OH <=> HNO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2NO + OH <=> HNO + H2O""",
)

entry(
    index = 638,
    label = "H2NO + HO2 <=> HNO + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(29000, 'cm^3/(mol*s)'), n=2.69, Ea=(-1600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2NO + HO2 <=> HNO + H2O2""",
)

entry(
    index = 639,
    label = "H2NO + O2 <=> HNO + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(230, 'cm^3/(mol*s)'), n=2.994, Ea=(16500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2NO + O2 <=> HNO + HO2""",
)

entry(
    index = 640,
    label = "H2NO + NH2 <=> HNO + NH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+12, 'cm^3/(mol*s)'), n=0, Ea=(1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2NO + NH2 <=> HNO + NH3""",
)

entry(
    index = 641,
    label = "H2NO + NO <=> HNO + HNO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(20000, 'cm^3/(mol*s)'), n=2, Ea=(13000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2NO + NO <=> HNO + HNO""",
)

entry(
    index = 642,
    label = "H2NO + NO2 <=> HONO + HNO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(44000, 'cm^3/(mol*s)'), n=2.64, Ea=(4040, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2NO + NO2 <=> HONO + HNO""",
)

entry(
    index = 643,
    label = "HNOH <=> HNO + H",
    degeneracy = 1,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(
            A = (2e+24, 'cm^3/(mol*s)'),
            n = -2.84,
            Ea = (58934, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        efficiencies = {'O': 10},
    ),
    shortDesc = u"""The chemkin file reaction is HNOH <=> HNO + H""",
)

entry(
    index = 644,
    label = "HNOH + H <=> NH2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HNOH + H <=> NH2 + OH""",
)

entry(
    index = 645,
    label = "HNOH + H <=> HNO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.8e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(378, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HNOH + H <=> HNO + H2""",
)

entry(
    index = 646,
    label = "HNOH + O <=> HNO + OH",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(7e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(3.3e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(-358, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is HNOH + O <=> HNO + OH""",
)

entry(
    index = 647,
    label = "HNOH + OH <=> HNO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+06, 'cm^3/(mol*s)'), n=2, Ea=(-1192, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HNOH + OH <=> HNO + H2O""",
)

entry(
    index = 648,
    label = "HNOH + HO2 <=> HNO + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(29000, 'cm^3/(mol*s)'), n=2.69, Ea=(-1600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HNOH + HO2 <=> HNO + H2O2""",
)

entry(
    index = 649,
    label = "HNOH + HO2 <=> NH2OH + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(29000, 'cm^3/(mol*s)'), n=2.69, Ea=(-1600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HNOH + HO2 <=> NH2OH + O2""",
)

entry(
    index = 650,
    label = "HNOH + O2 <=> HNO + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+12, 'cm^3/(mol*s)'), n=0, Ea=(25000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HNOH + O2 <=> HNO + HO2""",
)

entry(
    index = 651,
    label = "HNOH + NH2 <=> NH3 + HNO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.8e+06, 'cm^3/(mol*s)'),
        n = 1.94,
        Ea = (-1152, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is HNOH + NH2 <=> NH3 + HNO""",
)

entry(
    index = 652,
    label = "HNOH + NO2 <=> HONO + HNO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+11, 'cm^3/(mol*s)'), n=0, Ea=(2000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HNOH + NO2 <=> HONO + HNO""",
)

entry(
    index = 653,
    label = "NO + H <=> HNO",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.5e+15, 'cm^3/(mol*s)'), n=-0.41, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (2.4e+14, 'cm^6/(mol^2*s)'),
            n = 0.206,
            Ea = (-1550, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.82,
        T3 = (1e-30, 'K'),
        T1 = (1e+30, 'K'),
        T2 = (1e+30, 'K'),
        efficiencies = {'N#N': 1.6},
    ),
    shortDesc = u"""The chemkin file reaction is NO + H <=> HNO""",
)

entry(
    index = 654,
    label = "HNO + H <=> NO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.6e+10, 'cm^3/(mol*s)'), n=0.94, Ea=(495, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HNO + H <=> NO + H2""",
)

entry(
    index = 655,
    label = "HNO + O <=> NO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HNO + O <=> NO + OH""",
)

entry(
    index = 656,
    label = "HNO + OH <=> NO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.2e+09, 'cm^3/(mol*s)'),
        n = 1.189,
        Ea = (334, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is HNO + OH <=> NO + H2O""",
)

entry(
    index = 657,
    label = "HNO + HO2 <=> HNO2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2000, 'cm^3/(mol*s)'), n=2.36, Ea=(8980, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HNO + HO2 <=> HNO2 + OH""",
)

entry(
    index = 658,
    label = "HNO + O2 <=> HO2 + NO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(16000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HNO + O2 <=> HO2 + NO""",
)

entry(
    index = 659,
    label = "HNO + HNO <=> N2O + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+08, 'cm^3/(mol*s)'), n=0, Ea=(3100, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HNO + HNO <=> N2O + H2O""",
)

entry(
    index = 660,
    label = "HNO + NO2 <=> HONO + NO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(44000, 'cm^3/(mol*s)'), n=2.64, Ea=(4040, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HNO + NO2 <=> HONO + NO""",
)

entry(
    index = 661,
    label = "NO + HO2 <=> NO2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.1e+12, 'cm^3/(mol*s)'), n=0, Ea=(-497, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NO + HO2 <=> NO2 + OH""",
)

entry(
    index = 662,
    label = "NO + O <=> NO2",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.3e+15, 'cm^3/(mol*s)'), n=-0.75, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4.72e+24, 'cm^6/(mol^2*s)'),
            n = -2.87,
            Ea = (1550, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.75,
        T3 = (1000, 'K'),
        T1 = (100000, 'K'),
        T2 = (1e+30, 'K'),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is NO + O <=> NO2""",
)

entry(
    index = 663,
    label = "NO2 + H <=> NO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+14, 'cm^3/(mol*s)'), n=0, Ea=(362, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NO2 + H <=> NO + OH""",
)

entry(
    index = 664,
    label = "NO2 + O <=> NO + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.1e+14, 'cm^3/(mol*s)'), n=-0.52, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NO2 + O <=> NO + O2""",
)

entry(
    index = 665,
    label = "NO2 + HO2 <=> HONO + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.9, 'cm^3/(mol*s)'), n=3.32, Ea=(3044, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NO2 + HO2 <=> HONO + O2""",
)

entry(
    index = 666,
    label = "NO2 + HO2 <=> HNO2 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(19, 'cm^3/(mol*s)'), n=3.26, Ea=(4983, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NO2 + HO2 <=> HNO2 + O2""",
)

entry(
    index = 667,
    label = "NO2 + NO2 <=> NO + NO + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(27599, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NO2 + NO2 <=> NO + NO + O2""",
)

entry(
    index = 668,
    label = "NO2 + NO2 <=> NO3 + NO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (9.6e+09, 'cm^3/(mol*s)'),
        n = 0.73,
        Ea = (20900, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NO2 + NO2 <=> NO3 + NO""",
)

entry(
    index = 669,
    label = "NO + OH <=> HONO",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.1e+14, 'cm^3/(mol*s)'), n=-0.3, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.392e+23, 'cm^6/(mol^2*s)'),
            n = -2.5,
            Ea = (0, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.75,
        T3 = (1e-30, 'K'),
        T1 = (1e+30, 'K'),
        T2 = (1e+30, 'K'),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is NO + OH <=> HONO""",
)

entry(
    index = 670,
    label = "NO2 + H2 <=> HONO + H",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(18, 'cm^3/(mol*s)'), n=3.51, Ea=(26300, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(24, 'cm^3/(mol*s)'), n=3.62, Ea=(35800, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is NO2 + H2 <=> HONO + H""",
)

entry(
    index = 671,
    label = "HONO + H <=> HNO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.6e+10, 'cm^3/(mol*s)'),
        n = 0.86,
        Ea = (5000, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is HONO + H <=> HNO + OH""",
)

entry(
    index = 672,
    label = "HONO + H <=> NO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8.1e+06, 'cm^3/(mol*s)'),
        n = 1.89,
        Ea = (3850, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is HONO + H <=> NO + H2O""",
)

entry(
    index = 673,
    label = "HONO + O <=> NO2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+13, 'cm^3/(mol*s)'), n=0, Ea=(5960, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HONO + O <=> NO2 + OH""",
)

entry(
    index = 674,
    label = "HONO + OH <=> NO2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-520, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HONO + OH <=> NO2 + H2O""",
)

entry(
    index = 675,
    label = "HONO + NO2 <=> HONO2 + NO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+11, 'cm^3/(mol*s)'), n=0, Ea=(32700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HONO + NO2 <=> HONO2 + NO""",
)

entry(
    index = 676,
    label = "HONO + HONO <=> NO + NO2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.35, 'cm^3/(mol*s)'), n=3.64, Ea=(12140, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HONO + HONO <=> NO + NO2 + H2O""",
)

entry(
    index = 677,
    label = "HNO2 <=> HONO",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2.5e+14, 's^-1'), n=0, Ea=(32300, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(A=(3.1e+18, 'cm^3/(mol*s)'), n=0, Ea=(31500, 'cal/mol'), T0=(1, 'K')),
        alpha = 1.149,
        T3 = (1e-30, 'K'),
        T1 = (3125, 'K'),
        T2 = (1e+30, 'K'),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is HNO2 <=> HONO""",
)

entry(
    index = 678,
    label = "NO2 + H2 <=> HNO2 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(240, 'cm^3/(mol*s)'), n=3.15, Ea=(31100, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NO2 + H2 <=> HNO2 + H""",
)

entry(
    index = 679,
    label = "HNO2 + O <=> NO2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(2000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HNO2 + O <=> NO2 + OH""",
)

entry(
    index = 680,
    label = "HNO2 + OH <=> NO2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HNO2 + OH <=> NO2 + H2O""",
)

entry(
    index = 681,
    label = "NO2 + O <=> NO3",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.5e+12, 'cm^3/(mol*s)'), n=0.24, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(A=(2.5e+20, 'cm^6/(mol^2*s)'), n=-1.5, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        alpha = 0.71,
        T3 = (1e-30, 'K'),
        T1 = (1700, 'K'),
        T2 = (1e+30, 'K'),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is NO2 + O <=> NO3""",
)

entry(
    index = 682,
    label = "NO3 + H <=> NO2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NO3 + H <=> NO2 + OH""",
)

entry(
    index = 683,
    label = "NO3 + O <=> NO2 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NO3 + O <=> NO2 + O2""",
)

entry(
    index = 684,
    label = "NO3 + OH <=> NO2 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NO3 + OH <=> NO2 + HO2""",
)

entry(
    index = 685,
    label = "NO3 + HO2 <=> NO2 + O2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NO3 + HO2 <=> NO2 + O2 + OH""",
)

entry(
    index = 686,
    label = "NO3 + NO2 <=> NO + NO2 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+10, 'cm^3/(mol*s)'), n=0, Ea=(2940, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NO3 + NO2 <=> NO + NO2 + O2""",
)

entry(
    index = 687,
    label = "NO2 + OH <=> HONO2",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(A=(2.938e+25, 'cm^6/(mol^2*s)'), n=-3, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        alpha = 0.4,
        T3 = (1e-30, 'K'),
        T1 = (1e+30, 'K'),
        T2 = (1e+30, 'K'),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is NO2 + OH <=> HONO2""",
)

entry(
    index = 688,
    label = "HONO2 + H <=> H2 + NO3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.6e+08, 'cm^3/(mol*s)'),
        n = 1.5,
        Ea = (16400, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is HONO2 + H <=> H2 + NO3""",
)

entry(
    index = 689,
    label = "HONO2 + H <=> H2O + NO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(61, 'cm^3/(mol*s)'), n=3.3, Ea=(6285, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HONO2 + H <=> H2O + NO2""",
)

entry(
    index = 690,
    label = "HONO2 + H <=> OH + HONO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(380000, 'cm^3/(mol*s)'), n=2.3, Ea=(6976, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HONO2 + H <=> OH + HONO""",
)

entry(
    index = 691,
    label = "HONO2 + OH <=> H2O + NO3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+10, 'cm^3/(mol*s)'), n=0, Ea=(-1240, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HONO2 + OH <=> H2O + NO3""",
)

entry(
    index = 692,
    label = "N2O <=> N2 + O",
    degeneracy = 1,
    kinetics = Lindemann(
        arrheniusHigh = Arrhenius(A=(9.9e+10, 's^-1'), n=0, Ea=(57901, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(A=(6e+14, 'cm^3/(mol*s)'), n=0, Ea=(57444, 'cal/mol'), T0=(1, 'K')),
        efficiencies = {'[O][O]': 1.4, 'O': 12, 'N#N': 1.7},
    ),
    shortDesc = u"""The chemkin file reaction is N2O <=> N2 + O""",
)

entry(
    index = 693,
    label = "N2O + H <=> N2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (6.4e+07, 'cm^3/(mol*s)'),
        n = 1.835,
        Ea = (13492, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is N2O + H <=> N2 + OH""",
)

entry(
    index = 694,
    label = "N2O + O <=> NO + NO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.2e+13, 'cm^3/(mol*s)'), n=0, Ea=(27679, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is N2O + O <=> NO + NO""",
)

entry(
    index = 695,
    label = "N2O + O <=> N2 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.2e+13, 'cm^3/(mol*s)'), n=0, Ea=(27679, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is N2O + O <=> N2 + O2""",
)

entry(
    index = 696,
    label = "N2O + OH <=> N2 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.013, 'cm^3/(mol*s)'), n=4.72, Ea=(36560, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is N2O + OH <=> N2 + HO2""",
)

entry(
    index = 697,
    label = "N2O + OH <=> HNO + NO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (0.00012, 'cm^3/(mol*s)'),
        n = 4.33,
        Ea = (25080, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is N2O + OH <=> HNO + NO""",
)

entry(
    index = 698,
    label = "N2O + NO <=> NO2 + N2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (530000, 'cm^3/(mol*s)'),
        n = 2.23,
        Ea = (46280, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is N2O + NO <=> NO2 + N2""",
)

entry(
    index = 699,
    label = "NH2 + NH2 <=> N2H2 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.7e+08, 'cm^3/(mol*s)'),
        n = 1.62,
        Ea = (11783, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NH2 + NH2 <=> N2H2 + H2""",
)

entry(
    index = 700,
    label = "NH2 + NH2 <=> H2NN + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(72000, 'cm^3/(mol*s)'), n=1.88, Ea=(8802, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NH2 + NH2 <=> H2NN + H2""",
)

entry(
    index = 701,
    label = "NH2 + NH <=> N2H2 + H",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.3e+14, 'cm^3/(mol*s)'),
        n = -0.272,
        Ea = (-77, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NH2 + NH <=> N2H2 + H""",
)

entry(
    index = 702,
    label = "HNOH + NH2 <=> N2H3 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(10, 'cm^3/(mol*s)'), n=3.46, Ea=(-467, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HNOH + NH2 <=> N2H3 + OH""",
)

entry(
    index = 703,
    label = "HNOH + NH2 <=> H2NN + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8.8e+16, 'cm^3/(mol*s)'),
        n = -1.08,
        Ea = (1113, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is HNOH + NH2 <=> H2NN + H2O""",
)

entry(
    index = 704,
    label = "NH2 + NH2 <=> N2H4",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(
            A = (5.6e+14, 'cm^3/(mol*s)'),
            n = -0.414,
            Ea = (66, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        arrheniusLow = Arrhenius(
            A = (1.6e+34, 'cm^6/(mol^2*s)'),
            n = -5.49,
            Ea = (1987, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.31,
        T3 = (1e-30, 'K'),
        T1 = (1e+30, 'K'),
        T2 = (1e+30, 'K'),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is NH2 + NH2 <=> N2H4""",
)

entry(
    index = 705,
    label = "N2H4 <=> H2NN + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.4e+14, 's^-1'), n=0, Ea=(74911, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is N2H4 <=> H2NN + H2""",
)

entry(
    index = 706,
    label = "N2H4 + H <=> N2H3 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(2500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is N2H4 + H <=> N2H3 + H2""",
)

entry(
    index = 707,
    label = "N2H4 + H <=> NH3 + NH2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(230000, 'cm^3/(mol*s)'), n=1.42, Ea=(8202, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is N2H4 + H <=> NH3 + NH2""",
)

entry(
    index = 708,
    label = "N2H4 + O <=> N2H3 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.5e+11, 'cm^3/(mol*s)'), n=0, Ea=(-1270, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is N2H4 + O <=> N2H3 + OH""",
)

entry(
    index = 709,
    label = "N2H4 + O <=> N2H2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.9e+11, 'cm^3/(mol*s)'), n=0, Ea=(-1270, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is N2H4 + O <=> N2H2 + H2O""",
)

entry(
    index = 710,
    label = "N2H4 + OH <=> N2H3 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+13, 'cm^3/(mol*s)'), n=0, Ea=(-318, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is N2H4 + OH <=> N2H3 + H2O""",
)

entry(
    index = 711,
    label = "N2H4 + NH2 <=> N2H3 + NH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.76, 'cm^3/(mol*s)'), n=4, Ea=(4048, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is N2H4 + NH2 <=> N2H3 + NH3""",
)

entry(
    index = 712,
    label = "N2H4 + NO <=> N2H3 + HNO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(60, 'cm^3/(mol*s)'), n=3.16, Ea=(30845, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is N2H4 + NO <=> N2H3 + HNO""",
)

entry(
    index = 713,
    label = "N2H4 + NO2 <=> N2H3 + HONO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(82, 'cm^3/(mol*s)'), n=3.13, Ea=(8860, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is N2H4 + NO2 <=> N2H3 + HONO""",
)

entry(
    index = 714,
    label = "N2H4 + NO2 <=> N2H3 + HNO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.024, 'cm^3/(mol*s)'), n=4.14, Ea=(7946, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is N2H4 + NO2 <=> N2H3 + HNO2""",
)

entry(
    index = 715,
    label = "N2H2 + H <=> N2H3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+14, 'cm^3/(mol*s)'), n=0, Ea=(3871, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is N2H2 + H <=> N2H3""",
)

entry(
    index = 716,
    label = "N2H3 + H <=> N2H2 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(-10, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is N2H3 + H <=> N2H2 + H2""",
)

entry(
    index = 717,
    label = "N2H3 + O <=> N2H2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(-646, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is N2H3 + O <=> N2H2 + OH""",
)

entry(
    index = 718,
    label = "N2H3 + O <=> NH2 + HNO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is N2H3 + O <=> NH2 + HNO""",
)

entry(
    index = 719,
    label = "N2H3 + O => NH2 + NO + H",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is N2H3 + O => NH2 + NO + H""",
)

entry(
    index = 720,
    label = "N2H3 + OH <=> N2H2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+06, 'cm^3/(mol*s)'), n=2, Ea=(-1192, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is N2H3 + OH <=> N2H2 + H2O""",
)

entry(
    index = 721,
    label = "N2H3 + OH <=> H2NN + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is N2H3 + OH <=> H2NN + H2O""",
)

entry(
    index = 722,
    label = "N2H3 + OH <=> NH3 + HNO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(15000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is N2H3 + OH <=> NH3 + HNO""",
)

entry(
    index = 723,
    label = "N2H3 + HO2 <=> N2H2 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(14000, 'cm^3/(mol*s)'), n=2.69, Ea=(-1600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is N2H3 + HO2 <=> N2H2 + H2O2""",
)

entry(
    index = 724,
    label = "N2H3 + HO2 <=> N2H4 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(920000, 'cm^3/(mol*s)'), n=1.94, Ea=(2126, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is N2H3 + HO2 <=> N2H4 + O2""",
)

entry(
    index = 725,
    label = "N2H3 + NH2 <=> N2H2 + NH3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (920000, 'cm^3/(mol*s)'),
        n = 1.94,
        Ea = (-1152, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is N2H3 + NH2 <=> N2H2 + NH3""",
)

entry(
    index = 726,
    label = "N2H3 + NH2 <=> H2NN + NH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is N2H3 + NH2 <=> H2NN + NH3""",
)

entry(
    index = 727,
    label = "N2H3 + NH <=> N2H2 + NH2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is N2H3 + NH <=> N2H2 + NH2""",
)

entry(
    index = 728,
    label = "N2H2 <=> NNH + H",
    degeneracy = 1,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(
            A = (1.9e+27, 'cm^3/(mol*s)'),
            n = -3.05,
            Ea = (66107, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        efficiencies = {'O': 7},
    ),
    shortDesc = u"""The chemkin file reaction is N2H2 <=> NNH + H""",
)

entry(
    index = 729,
    label = "N2H2 + H <=> NNH + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.1e+14, 'cm^3/(mol*s)'), n=0, Ea=(3128, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is N2H2 + H <=> NNH + H2""",
)

entry(
    index = 730,
    label = "N2H2 + O <=> NNH + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.3e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(497, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is N2H2 + O <=> NNH + OH""",
)

entry(
    index = 731,
    label = "N2H2 + O <=> NH2 + NO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is N2H2 + O <=> NH2 + NO""",
)

entry(
    index = 732,
    label = "N2H2 + OH <=> NNH + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(59, 'cm^3/(mol*s)'), n=3.4, Ea=(1360, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is N2H2 + OH <=> NNH + H2O""",
)

entry(
    index = 733,
    label = "N2H2 + NH2 <=> NNH + NH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.088, 'cm^3/(mol*s)'), n=4.05, Ea=(1610, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is N2H2 + NH2 <=> NNH + NH3""",
)

entry(
    index = 734,
    label = "N2H2 + NH <=> NNH + NH2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+06, 'cm^3/(mol*s)'), n=2, Ea=(-1192, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is N2H2 + NH <=> NNH + NH2""",
)

entry(
    index = 735,
    label = "N2H2 + NO <=> N2O + NH2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+12, 'cm^3/(mol*s)'), n=0, Ea=(11922, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is N2H2 + NO <=> N2O + NH2""",
)

entry(
    index = 736,
    label = "H2NN <=> NNH + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.4e+26, 's^-1'), n=-4.83, Ea=(46228, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2NN <=> NNH + H""",
)

entry(
    index = 737,
    label = "H2NN <=> N2 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.5e+14, 's^-1'), n=0, Ea=(52785, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2NN <=> N2 + H2""",
)

entry(
    index = 738,
    label = "H2NN <=> N2H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+14, 's^-1'), n=0, Ea=(46931, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2NN <=> N2H2""",
)

entry(
    index = 739,
    label = "H2NN + H <=> NNH + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.8e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(-894, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2NN + H <=> NNH + H2""",
)

entry(
    index = 740,
    label = "H2NN + H <=> N2H2 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2NN + H <=> N2H2 + H""",
)

entry(
    index = 741,
    label = "H2NN + O <=> NNH + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.3e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(-894, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2NN + O <=> NNH + OH""",
)

entry(
    index = 742,
    label = "H2NN + O <=> NH2 + NO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2NN + O <=> NH2 + NO""",
)

entry(
    index = 743,
    label = "H2NN + OH <=> NNH + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+06, 'cm^3/(mol*s)'), n=2, Ea=(-1192, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2NN + OH <=> NNH + H2O""",
)

entry(
    index = 744,
    label = "H2NN + OH => NH2 + NO + H",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2NN + OH => NH2 + NO + H""",
)

entry(
    index = 745,
    label = "H2NN + HO2 => NH2 + NO + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(9e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2NN + HO2 => NH2 + NO + OH""",
)

entry(
    index = 746,
    label = "H2NN + HO2 <=> NNH + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(29000, 'cm^3/(mol*s)'), n=2.69, Ea=(-1600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2NN + HO2 <=> NNH + H2O2""",
)

entry(
    index = 747,
    label = "H2NN + O2 <=> NH2 + NO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(5961, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2NN + O2 <=> NH2 + NO2""",
)

entry(
    index = 748,
    label = "H2NN + NH2 <=> NNH + NH3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.8e+06, 'cm^3/(mol*s)'),
        n = 1.94,
        Ea = (-1152, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is H2NN + NH2 <=> NNH + NH3""",
)

entry(
    index = 749,
    label = "HCN <=> H + CN",
    degeneracy = 1,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(
            A = (3.4e+35, 'cm^3/(mol*s)'),
            n = -5.13,
            Ea = (133000, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        efficiencies = {'[O][O]': 1.5, 'O': 10, 'N#N': 0},
    ),
    shortDesc = u"""The chemkin file reaction is HCN <=> H + CN""",
)

entry(
    index = 750,
    label = "HCN + N2 <=> H + CN + N2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.6e+26, 'cm^3/(mol*s)'),
        n = -2.6,
        Ea = (124890, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is HCN + N2 <=> H + CN + N2""",
)

entry(
    index = 751,
    label = "HCN <=> HNC",
    degeneracy = 1,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(
            A = (1.6e+26, 'cm^3/(mol*s)'),
            n = -3.23,
            Ea = (54900, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        efficiencies = {'O=C=O': 2, 'O': 7, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is HCN <=> HNC""",
)

entry(
    index = 752,
    label = "CN + H2 <=> HCN + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(110000, 'cm^3/(mol*s)'), n=2.6, Ea=(1908, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CN + H2 <=> HCN + H""",
)

entry(
    index = 753,
    label = "HCN + O <=> NCO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(14000, 'cm^3/(mol*s)'), n=2.64, Ea=(4980, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCN + O <=> NCO + H""",
)

entry(
    index = 754,
    label = "HCN + O <=> NH + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3500, 'cm^3/(mol*s)'), n=2.64, Ea=(4980, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCN + O <=> NH + CO""",
)

entry(
    index = 755,
    label = "HCN + O <=> CN + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.2e+10, 'cm^3/(mol*s)'),
        n = 0.4,
        Ea = (20665, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is HCN + O <=> CN + OH""",
)

entry(
    index = 756,
    label = "HCN + OH <=> CN + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.9e+06, 'cm^3/(mol*s)'),
        n = 1.83,
        Ea = (10300, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is HCN + OH <=> CN + H2O""",
)

entry(
    index = 757,
    label = "HCN + OH <=> HOCN + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(59000, 'cm^3/(mol*s)'), n=2.4, Ea=(12500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCN + OH <=> HOCN + H""",
)

entry(
    index = 758,
    label = "HCN + OH <=> HNCO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.002, 'cm^3/(mol*s)'), n=4, Ea=(1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCN + OH <=> HNCO + H""",
)

entry(
    index = 759,
    label = "HCN + OH <=> NH2 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.00078, 'cm^3/(mol*s)'), n=4, Ea=(4000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCN + OH <=> NH2 + CO""",
)

entry(
    index = 760,
    label = "HCN + O2 <=> CN + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(75100, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCN + O2 <=> CN + HO2""",
)

entry(
    index = 761,
    label = "HCN + CN <=> NCCN + H",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.5e+07, 'cm^3/(mol*s)'),
        n = 1.71,
        Ea = (1530, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is HCN + CN <=> NCCN + H""",
)

entry(
    index = 762,
    label = "HNC + H <=> HCN + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(3600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HNC + H <=> HCN + H""",
)

entry(
    index = 763,
    label = "HNC + O <=> NH + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.6e+12, 'cm^3/(mol*s)'), n=0, Ea=(2200, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HNC + O <=> NH + CO""",
)

entry(
    index = 764,
    label = "HNC + OH <=> HNCO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.6e+12, 'cm^3/(mol*s)'), n=0, Ea=(-479, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HNC + OH <=> HNCO + H""",
)

entry(
    index = 765,
    label = "HNC + OH <=> CN + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(300, 'cm^3/(mol*s)'), n=3.16, Ea=(10600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HNC + OH <=> CN + H2O""",
)

entry(
    index = 766,
    label = "CN + O <=> CO + N",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.9e+12, 'cm^3/(mol*s)'), n=0.46, Ea=(723, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CN + O <=> CO + N""",
)

entry(
    index = 767,
    label = "CN + OH <=> NCO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+15, 'cm^3/(mol*s)'), n=-0.437, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CN + OH <=> NCO + H""",
)

entry(
    index = 768,
    label = "CN + O2 <=> NCO + O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(-417, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CN + O2 <=> NCO + O""",
)

entry(
    index = 769,
    label = "CN + O2 <=> NO + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.4e+12, 'cm^3/(mol*s)'), n=0, Ea=(-417, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CN + O2 <=> NO + CO""",
)

entry(
    index = 770,
    label = "CN + NO <=> NCO + N",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(42100, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CN + NO <=> NCO + N""",
)

entry(
    index = 771,
    label = "CN + NO2 <=> NCO + NO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.3e+15, 'cm^3/(mol*s)'),
        n = -0.752,
        Ea = (344, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CN + NO2 <=> NCO + NO""",
)

entry(
    index = 772,
    label = "CN + NO2 <=> CO + N2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.9e+14, 'cm^3/(mol*s)'),
        n = -0.752,
        Ea = (344, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CN + NO2 <=> CO + N2O""",
)

entry(
    index = 773,
    label = "CN + NO2 <=> N2 + CO2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.7e+14, 'cm^3/(mol*s)'),
        n = -0.752,
        Ea = (344, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CN + NO2 <=> N2 + CO2""",
)

entry(
    index = 774,
    label = "CN + HNO <=> HCN + NO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CN + HNO <=> HCN + NO""",
)

entry(
    index = 775,
    label = "CN + HONO <=> HCN + NO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CN + HONO <=> HCN + NO2""",
)

entry(
    index = 776,
    label = "CN + HNCO <=> HCN + NCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.2e+13, 'cm^3/(mol*s)'), n=0, Ea=(9670, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CN + HNCO <=> HCN + NCO""",
)

entry(
    index = 777,
    label = "CN + HNCO <=> HNCN + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.7e+13, 'cm^3/(mol*s)'), n=0, Ea=(9430, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CN + HNCO <=> HNCN + CO""",
)

entry(
    index = 778,
    label = "CN + NCO <=> NCN + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+13, 'cm^3/(mol*s)'), n=0, Ea=(-308, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CN + NCO <=> NCN + CO""",
)

entry(
    index = 779,
    label = "HNCO <=> CO + NH",
    degeneracy = 1,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(A=(1.1e+16, 'cm^3/(mol*s)'), n=0, Ea=(86000, 'cal/mol'), T0=(1, 'K')),
        efficiencies = {'N#N': 1.5},
    ),
    shortDesc = u"""The chemkin file reaction is HNCO <=> CO + NH""",
)

entry(
    index = 780,
    label = "HNCO + H <=> NH2 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(36000, 'cm^3/(mol*s)'), n=2.49, Ea=(2345, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HNCO + H <=> NH2 + CO""",
)

entry(
    index = 781,
    label = "HNCO + H <=> NCO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+07, 'cm^3/(mol*s)'), n=1.66, Ea=(13900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HNCO + H <=> NCO + H2""",
)

entry(
    index = 782,
    label = "HNCO + O <=> NCO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.2e+06, 'cm^3/(mol*s)'),
        n = 2.11,
        Ea = (11430, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is HNCO + O <=> NCO + OH""",
)

entry(
    index = 783,
    label = "HNCO + O <=> NH + CO2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (9.6e+07, 'cm^3/(mol*s)'),
        n = 1.41,
        Ea = (8520, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is HNCO + O <=> NH + CO2""",
)

entry(
    index = 784,
    label = "HNCO + O <=> HNO + CO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.5e+08, 'cm^3/(mol*s)'),
        n = 1.57,
        Ea = (44012, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is HNCO + O <=> HNO + CO""",
)

entry(
    index = 785,
    label = "HNCO + OH <=> NCO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.5e+07, 'cm^3/(mol*s)'), n=1.5, Ea=(3600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HNCO + OH <=> NCO + H2O""",
)

entry(
    index = 786,
    label = "HNCO + HO2 <=> NCO + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+11, 'cm^3/(mol*s)'), n=0, Ea=(22000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HNCO + HO2 <=> NCO + H2O2""",
)

entry(
    index = 787,
    label = "HNCO + O2 <=> HNO + CO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(58900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HNCO + O2 <=> HNO + CO2""",
)

entry(
    index = 788,
    label = "HNCO + H2O <=> NH3 + CO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(48500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HNCO + H2O <=> NH3 + CO2""",
)

entry(
    index = 789,
    label = "HNCO + NH <=> NH2 + NCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(23700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HNCO + NH <=> NH2 + NCO""",
)

entry(
    index = 790,
    label = "HNCO + HNCO <=> HNCNH + CO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.9e+11, 'cm^3/(mol*s)'), n=0, Ea=(42100, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HNCO + HNCO <=> HNCNH + CO2""",
)

entry(
    index = 791,
    label = "HOCN + H <=> HNCO + H",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.1e+08, 'cm^3/(mol*s)'),
        n = 0.84,
        Ea = (1917, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is HOCN + H <=> HNCO + H""",
)

entry(
    index = 792,
    label = "HOCN + H <=> NH2 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.2e+08, 'cm^3/(mol*s)'),
        n = 0.61,
        Ea = (2076, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is HOCN + H <=> NH2 + CO""",
)

entry(
    index = 793,
    label = "HOCN + H <=> H2 + NCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(6617, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HOCN + H <=> H2 + NCO""",
)

entry(
    index = 794,
    label = "HOCN + O <=> OH + NCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(4133, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HOCN + O <=> OH + NCO""",
)

entry(
    index = 795,
    label = "HOCN + OH <=> H2O + NCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+06, 'cm^3/(mol*s)'), n=2, Ea=(-248, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HOCN + OH <=> H2O + NCO""",
)

entry(
    index = 796,
    label = "HOCN + NH2 <=> NCO + NH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(920000, 'cm^3/(mol*s)'), n=1.94, Ea=(3646, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HOCN + NH2 <=> NCO + NH3""",
)

entry(
    index = 797,
    label = "HCNO <=> HCN + O",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.1, 1, 10], 'atm'),
        arrhenius = [
            Arrhenius(A=(2e+30, 's^-1'), n=-6.03, Ea=(60733, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(4.2e+31, 's^-1'), n=-6.12, Ea=(61210, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(5.9e+31, 's^-1'), n=-5.85, Ea=(61935, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is HCNO <=> HCN + O""",
)

entry(
    index = 798,
    label = "HCNO + H <=> HCN + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (7.2e+10, 'cm^3/(mol*s)'),
        n = 0.841,
        Ea = (8612, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is HCNO + H <=> HCN + OH""",
)

entry(
    index = 799,
    label = "HCNO + O <=> HCO + NO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCNO + O <=> HCO + NO""",
)

entry(
    index = 800,
    label = "HCNO + O <=> NCO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCNO + O <=> NCO + OH""",
)

entry(
    index = 801,
    label = "HCNO + OH <=> CO + H2NO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(-1490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCNO + OH <=> CO + H2NO""",
)

entry(
    index = 802,
    label = "HCNO + OH <=> HCO + HNO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCNO + OH <=> HCO + HNO""",
)

entry(
    index = 803,
    label = "HCNO + CN <=> HCN + NCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCNO + CN <=> HCN + NCO""",
)

entry(
    index = 804,
    label = "HCNO + NCO <=> HCN + NO + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.6e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCNO + NCO <=> HCN + NO + CO""",
)

entry(
    index = 805,
    label = "NCO <=> N + CO",
    degeneracy = 1,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(A=(2.2e+14, 'cm^3/(mol*s)'), n=0, Ea=(54050, 'cal/mol'), T0=(1, 'K')),
        efficiencies = {'N#N': 1.5},
    ),
    shortDesc = u"""The chemkin file reaction is NCO <=> N + CO""",
)

entry(
    index = 806,
    label = "NCO + H <=> CO + NH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.2e+13, 'cm^3/(mol*s)'), n=0, Ea=(1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NCO + H <=> CO + NH""",
)

entry(
    index = 807,
    label = "NCO + O <=> NO + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+15, 'cm^3/(mol*s)'), n=-0.5, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NCO + O <=> NO + CO""",
)

entry(
    index = 808,
    label = "NCO + OH <=> HON + CO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.3e+12, 'cm^3/(mol*s)'),
        n = -0.07,
        Ea = (5126, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NCO + OH <=> HON + CO""",
)

entry(
    index = 809,
    label = "NCO + OH <=> H + CO + NO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8.3e+12, 'cm^3/(mol*s)'),
        n = -0.05,
        Ea = (18042, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NCO + OH <=> H + CO + NO""",
)

entry(
    index = 810,
    label = "NCO + HO2 <=> HNCO + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NCO + HO2 <=> HNCO + O2""",
)

entry(
    index = 811,
    label = "NCO + O2 <=> NO + CO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+12, 'cm^3/(mol*s)'), n=0, Ea=(20000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NCO + O2 <=> NO + CO2""",
)

entry(
    index = 812,
    label = "NCO + NO <=> N2O + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+19, 'cm^3/(mol*s)'), n=-2.19, Ea=(1743, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NCO + NO <=> N2O + CO""",
)

entry(
    index = 813,
    label = "NCO + NO <=> N2 + CO2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.5e+21, 'cm^3/(mol*s)'),
        n = -2.74,
        Ea = (1824, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NCO + NO <=> N2 + CO2""",
)

entry(
    index = 814,
    label = "NCO + NO2 <=> CO + NO + NO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.5e+11, 'cm^3/(mol*s)'), n=0, Ea=(-707, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NCO + NO2 <=> CO + NO + NO""",
)

entry(
    index = 815,
    label = "NCO + NO2 <=> CO2 + N2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+12, 'cm^3/(mol*s)'), n=0, Ea=(-707, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NCO + NO2 <=> CO2 + N2O""",
)

entry(
    index = 816,
    label = "NCO + HNO <=> HNCO + NO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NCO + HNO <=> HNCO + NO""",
)

entry(
    index = 817,
    label = "NCO + HONO <=> HNCO + NO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.6e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NCO + HONO <=> HNCO + NO2""",
)

entry(
    index = 818,
    label = "NCO + NH3 <=> HNCO + NH2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(28000, 'cm^3/(mol*s)'), n=2.48, Ea=(980, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NCO + NH3 <=> HNCO + NH2""",
)

entry(
    index = 819,
    label = "NCO + N <=> N2 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NCO + N <=> N2 + CO""",
)

entry(
    index = 820,
    label = "NCO + NCO <=> CO + CO + N2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NCO + NCO <=> CO + CO + N2""",
)

entry(
    index = 821,
    label = "CH3CN <=> CH2CN + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.9e+14, 's^-1'), n=0, Ea=(94940, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CN <=> CH2CN + H""",
)

entry(
    index = 822,
    label = "CH3CN + H <=> HCN + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.4e+10, 'cm^3/(mol*s)'), n=0.8, Ea=(6800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CN + H <=> HCN + CH3""",
)

entry(
    index = 823,
    label = "CH3CN + H <=> HNC + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.8e+15, 'cm^3/(mol*s)'),
        n = -0.32,
        Ea = (20030, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3CN + H <=> HNC + CH3""",
)

entry(
    index = 824,
    label = "CH3CN + H <=> CH2CN + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(60000, 'cm^3/(mol*s)'), n=3.01, Ea=(8522, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CN + H <=> CH2CN + H2""",
)

entry(
    index = 825,
    label = "CH3CN + O <=> CH3 + NCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+09, 'cm^3/(mol*s)'), n=1.8, Ea=(8130, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CN + O <=> CH3 + NCO""",
)

entry(
    index = 826,
    label = "CH3CN + O <=> CH2CN + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.7e+08, 'cm^3/(mol*s)'),
        n = 1.18,
        Ea = (14360, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3CN + O <=> CH2CN + OH""",
)

entry(
    index = 827,
    label = "CH3CN + OH <=> CH2CN + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+07, 'cm^3/(mol*s)'), n=2, Ea=(2000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CN + OH <=> CH2CN + H2O""",
)

entry(
    index = 828,
    label = "CH3CN + CH3 <=> CH2CN + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(7000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CN + CH3 <=> CH2CN + CH4""",
)

entry(
    index = 829,
    label = "CH3CN + CN <=> CH2CN + HCN",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(2000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CN + CN <=> CH2CN + HCN""",
)

entry(
    index = 830,
    label = "CH2CN + O <=> CH2O + CN",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+12, 'cm^3/(mol*s)'), n=0.64, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CN + O <=> CH2O + CN""",
)

entry(
    index = 831,
    label = "CH2OH + CN <=> CH2CN + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2OH + CN <=> CH2CN + OH""",
)

entry(
    index = 832,
    label = "CO + NO2 <=> NO + CO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+13, 'cm^3/(mol*s)'), n=0, Ea=(33800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CO + NO2 <=> NO + CO2""",
)

entry(
    index = 833,
    label = "CO + N2O <=> N2 + CO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.7e+11, 'cm^3/(mol*s)'), n=0, Ea=(20237, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CO + N2O <=> N2 + CO2""",
)

entry(
    index = 834,
    label = "CO2 + N <=> NO + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+12, 'cm^3/(mol*s)'), n=0, Ea=(20000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CO2 + N <=> NO + CO""",
)

entry(
    index = 835,
    label = "CO2 + CN <=> NCO + CO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.7e+06, 'cm^3/(mol*s)'),
        n = 2.16,
        Ea = (26900, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CO2 + CN <=> NCO + CO""",
)

entry(
    index = 836,
    label = "HOCO + NO <=> CO + HONO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HOCO + NO <=> CO + HONO""",
)

entry(
    index = 837,
    label = "CH2O + NO2 <=> HONO + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.4e-07, 'cm^3/(mol*s)'),
        n = 5.64,
        Ea = (9220, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2O + NO2 <=> HONO + HCO""",
)

entry(
    index = 838,
    label = "CH2O + NO2 <=> HNO2 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.11, 'cm^3/(mol*s)'), n=4.22, Ea=(19850, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2O + NO2 <=> HNO2 + HCO""",
)

entry(
    index = 839,
    label = "CH2O + CN <=> HCO + HCN",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1700, 'cm^3/(mol*s)'), n=2.72, Ea=(-1427, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2O + CN <=> HCO + HCN""",
)

entry(
    index = 840,
    label = "CH2O + NCO <=> HNCO + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2O + NCO <=> HNCO + HCO""",
)

entry(
    index = 841,
    label = "HCO + NO <=> HNO + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.9e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCO + NO <=> HNO + CO""",
)

entry(
    index = 842,
    label = "HCO + HNO <=> NO + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.58, 'cm^3/(mol*s)'), n=3.84, Ea=(115, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCO + HNO <=> NO + CH2O""",
)

entry(
    index = 843,
    label = "HCO + NO2 <=> NO + CO2 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCO + NO2 <=> NO + CO2 + H""",
)

entry(
    index = 844,
    label = "HCO + NO2 <=> HONO + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCO + NO2 <=> HONO + CO""",
)

entry(
    index = 845,
    label = "HCO + NO2 <=> NO + CO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCO + NO2 <=> NO + CO + OH""",
)

entry(
    index = 846,
    label = "HCO + NCO <=> HNCO + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCO + NCO <=> HNCO + CO""",
)

entry(
    index = 847,
    label = "CH4 + NH2 <=> CH3 + NH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1500, 'cm^3/(mol*s)'), n=3.01, Ea=(9940, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH4 + NH2 <=> CH3 + NH3""",
)

entry(
    index = 848,
    label = "CH4 + NO2 <=> CH3 + HONO",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(0.11, 'cm^3/(mol*s)'), n=4.28, Ea=(26300, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(74, 'cm^3/(mol*s)'), n=3.42, Ea=(33100, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH4 + NO2 <=> CH3 + HONO""",
)

entry(
    index = 849,
    label = "CH4 + NO2 <=> CH3 + HNO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.4, 'cm^3/(mol*s)'), n=4.18, Ea=(31200, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH4 + NO2 <=> CH3 + HNO2""",
)

entry(
    index = 850,
    label = "CH4 + CN <=> CH3 + HCN",
    degeneracy = 1,
    kinetics = Arrhenius(A=(860000, 'cm^3/(mol*s)'), n=2.3, Ea=(-32, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH4 + CN <=> CH3 + HCN""",
)

entry(
    index = 851,
    label = "CH4 + NCO <=> CH3 + HNCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(8120, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH4 + NCO <=> CH3 + HNCO""",
)

entry(
    index = 852,
    label = "CH3 + NH2 <=> CH4 + NH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.8e+06, 'cm^3/(mol*s)'),
        n = 1.94,
        Ea = (9210, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3 + NH2 <=> CH4 + NH""",
)

entry(
    index = 853,
    label = "CH3 + NH2 <=> CH2 + NH3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.6e+06, 'cm^3/(mol*s)'),
        n = 1.87,
        Ea = (7570, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3 + NH2 <=> CH2 + NH3""",
)

entry(
    index = 854,
    label = "CH3 + NH <=> CH2NH + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + NH <=> CH2NH + H""",
)

entry(
    index = 855,
    label = "CH3 + NH <=> N + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(820000, 'cm^3/(mol*s)'), n=1.87, Ea=(5852, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + NH <=> N + CH4""",
)

entry(
    index = 856,
    label = "CH3 + N <=> H2CN + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + N <=> H2CN + H""",
)

entry(
    index = 857,
    label = "CH3 + H2NO <=> CH3O + NH2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + H2NO <=> CH3O + NH2""",
)

entry(
    index = 858,
    label = "CH3 + H2NO <=> CH4 + HNO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.6e+06, 'cm^3/(mol*s)'),
        n = 1.87,
        Ea = (2961, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3 + H2NO <=> CH4 + HNO""",
)

entry(
    index = 859,
    label = "CH3 + HNO <=> NO + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.5e+11, 'cm^3/(mol*s)'), n=0.76, Ea=(348, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + HNO <=> NO + CH4""",
)

entry(
    index = 860,
    label = "CH3 + HNO <=> CH3NO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8100, 'cm^3/(mol*s)'), n=2.4, Ea=(6160, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + HNO <=> CH3NO + H""",
)

entry(
    index = 861,
    label = "CH3 + NO <=> HCN + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.15, 'cm^3/(mol*s)'), n=3.52, Ea=(3950, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + NO <=> HCN + H2O""",
)

entry(
    index = 862,
    label = "CH3 + NO <=> H2CN + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.15, 'cm^3/(mol*s)'), n=3.52, Ea=(3950, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + NO <=> H2CN + OH""",
)

entry(
    index = 863,
    label = "CH3 + NO2 <=> CH3O + NO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + NO2 <=> CH3O + NO""",
)

entry(
    index = 864,
    label = "CH3 + CN <=> CH2CN + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + CN <=> CH2CN + H""",
)

entry(
    index = 865,
    label = "CH3 + HOCN <=> CH3CN + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(2000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + HOCN <=> CH3CN + OH""",
)

entry(
    index = 866,
    label = "CH2 + N2 <=> HCN + NH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(74000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2 + N2 <=> HCN + NH""",
)

entry(
    index = 867,
    label = "CH2 + N <=> HCN + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2 + N <=> HCN + H""",
)

entry(
    index = 868,
    label = "CH2 + NO <=> HCNO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.1e+12, 'cm^3/(mol*s)'), n=0, Ea=(-378, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2 + NO <=> HCNO + H""",
)

entry(
    index = 869,
    label = "CH2 + NO <=> HCN + OH",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(3.9e+11, 'cm^3/(mol*s)'), n=0, Ea=(-378, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2 + NO <=> HCN + OH""",
)

entry(
    index = 870,
    label = "CH2 + NO2 <=> CH2O + NO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.9e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2 + NO2 <=> CH2O + NO""",
)

entry(
    index = 871,
    label = "CH2(S) + NH3 <=> CH2NH2 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2(S) + NH3 <=> CH2NH2 + H""",
)

entry(
    index = 872,
    label = "CH2(S) + NO <=> CH2 + NO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2(S) + NO <=> CH2 + NO""",
)

entry(
    index = 873,
    label = "CH2(S) + N2O <=> CH2O + N2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2(S) + N2O <=> CH2O + N2""",
)

entry(
    index = 874,
    label = "CH + N2 <=> H + NCN",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (5.9e+08, 'cm^3/(mol*s)'),
                n = 1.06,
                Ea = (15960, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (5.9e+08, 'cm^3/(mol*s)'),
                n = 1.06,
                Ea = (15950, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (6.2e+08, 'cm^3/(mol*s)'),
                n = 1.05,
                Ea = (15960, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(A=(9e+08, 'cm^3/(mol*s)'), n=1.01, Ea=(16120, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (2.5e+09, 'cm^3/(mol*s)'),
                n = 0.89,
                Ea = (16620, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (9.2e+09, 'cm^3/(mol*s)'),
                n = 0.75,
                Ea = (17410, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(A=(3e+10, 'cm^3/(mol*s)'), n=0.62, Ea=(18480, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (3.8e+10, 'cm^3/(mol*s)'),
                n = 0.62,
                Ea = (19460, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (9.5e+09, 'cm^3/(mol*s)'),
                n = 0.81,
                Ea = (20340, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH + N2 <=> H + NCN""",
)

entry(
    index = 875,
    label = "CH + N2 <=> HNCN",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (6e+23, 'cm^3/(mol*s)'),
                n = -4.41,
                Ea = (14410, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (7.3e+23, 'cm^3/(mol*s)'),
                n = -4.3,
                Ea = (14760, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (7.3e+23, 'cm^3/(mol*s)'),
                n = -4.17,
                Ea = (15200, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(A=(4.8e+23, 'cm^3/(mol*s)'), n=-4, Ea=(15570, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (1.4e+23, 'cm^3/(mol*s)'),
                n = -3.74,
                Ea = (15820, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.7e+22, 'cm^3/(mol*s)'),
                n = -3.38,
                Ea = (15840, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (6.8e+20, 'cm^3/(mol*s)'),
                n = -2.9,
                Ea = (15690, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.8e+19, 'cm^3/(mol*s)'),
                n = -2.37,
                Ea = (15430, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (3.1e+17, 'cm^3/(mol*s)'),
                n = -1.78,
                Ea = (15240, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH + N2 <=> HNCN""",
)

entry(
    index = 876,
    label = "CH + NH3 <=> H2CN + H + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(-630, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH + NH3 <=> H2CN + H + H""",
)

entry(
    index = 877,
    label = "CH + N <=> CN + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH + N <=> CN + H""",
)

entry(
    index = 878,
    label = "CH + NO <=> CO + NH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.1e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH + NO <=> CO + NH""",
)

entry(
    index = 879,
    label = "CH + NO <=> NCO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH + NO <=> NCO + H""",
)

entry(
    index = 880,
    label = "CH + NO <=> HCN + O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.9e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH + NO <=> HCN + O""",
)

entry(
    index = 881,
    label = "CH + NO <=> CN + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.1e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH + NO <=> CN + OH""",
)

entry(
    index = 882,
    label = "CH + NO <=> HCO + N",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH + NO <=> HCO + N""",
)

entry(
    index = 883,
    label = "CH + NO2 <=> HCO + NO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH + NO2 <=> HCO + NO""",
)

entry(
    index = 884,
    label = "CN + N <=> C + N2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.9e+14, 'cm^3/(mol*s)'), n=-0.4, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CN + N <=> C + N2""",
)

entry(
    index = 885,
    label = "CH + N2O <=> HCN + NO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.9e+13, 'cm^3/(mol*s)'), n=0, Ea=(-511, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH + N2O <=> HCN + NO""",
)

entry(
    index = 886,
    label = "C + NO <=> CN + O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C + NO <=> CN + O""",
)

entry(
    index = 887,
    label = "C + NO <=> CO + N",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C + NO <=> CO + N""",
)

entry(
    index = 888,
    label = "C + N2O <=> CN + NO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C + N2O <=> CN + NO""",
)

entry(
    index = 889,
    label = "CH3OH + NO2 <=> HONO + CH2OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(150, 'cm^3/(mol*s)'), n=3.32, Ea=(20035, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OH + NO2 <=> HONO + CH2OH""",
)

entry(
    index = 890,
    label = "CH3OH + NO2 <=> HNO2 + CH2OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2400, 'cm^3/(mol*s)'), n=2.9, Ea=(27470, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OH + NO2 <=> HNO2 + CH2OH""",
)

entry(
    index = 891,
    label = "CH3O + HNO <=> NO + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O + HNO <=> NO + CH3OH""",
)

entry(
    index = 892,
    label = "CH3O + NO <=> HNO + CH2O",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(7.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(2017, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2.5e+18, 'cm^3/(mol*s)'), n=-2.56, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH3O + NO <=> HNO + CH2O""",
)

entry(
    index = 893,
    label = "CH3O + NO2 <=> HONO + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+12, 'cm^3/(mol*s)'), n=0, Ea=(2285, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O + NO2 <=> HONO + CH2O""",
)

entry(
    index = 894,
    label = "CH2OH + HNO <=> NO + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2OH + HNO <=> NO + CH3OH""",
)

entry(
    index = 895,
    label = "CH2OH + NO <=> HNCO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2OH + NO <=> HNCO + H2O""",
)

entry(
    index = 896,
    label = "CH2OH + NO2 <=> HONO + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2OH + NO2 <=> HONO + CH2O""",
)

entry(
    index = 897,
    label = "CH3OO + NO <=> CH3O + NO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.4e+12, 'cm^3/(mol*s)'), n=0, Ea=(-715, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OO + NO <=> CH3O + NO2""",
)

entry(
    index = 898,
    label = "C2H6 + NH2 <=> C2H5 + NH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(45, 'cm^3/(mol*s)'), n=3.46, Ea=(5600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H6 + NH2 <=> C2H5 + NH3""",
)

entry(
    index = 899,
    label = "C2H6 + NO2 <=> C2H5 + HONO",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(3.3, 'cm^3/(mol*s)'), n=3.84, Ea=(23900, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(85, 'cm^3/(mol*s)'), n=3.45, Ea=(32000, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H6 + NO2 <=> C2H5 + HONO""",
)

entry(
    index = 900,
    label = "C2H6 + NO2 <=> C2H5 + HNO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(320, 'cm^3/(mol*s)'), n=3.19, Ea=(26500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H6 + NO2 <=> C2H5 + HNO2""",
)

entry(
    index = 901,
    label = "C2H6 + CN <=> C2H5 + HCN",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+08, 'cm^3/(mol*s)'), n=1.8, Ea=(-994, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H6 + CN <=> C2H5 + HCN""",
)

entry(
    index = 902,
    label = "C2H6 + NCO <=> C2H5 + HNCO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.5e-09, 'cm^3/(mol*s)'),
        n = 6.89,
        Ea = (-2910, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H6 + NCO <=> C2H5 + HNCO""",
)

entry(
    index = 903,
    label = "C2H5 + N <=> C2H4 + NH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5 + N <=> C2H4 + NH""",
)

entry(
    index = 904,
    label = "C2H5 + N <=> CH3 + H2CN",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5 + N <=> CH3 + H2CN""",
)

entry(
    index = 905,
    label = "C2H5 + NO2 <=> CH3CH2O + NO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=-0.2, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5 + NO2 <=> CH3CH2O + NO""",
)

entry(
    index = 906,
    label = "C2H4 + NH2 <=> C2H3 + NH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.3e+12, 'cm^3/(mol*s)'), n=0, Ea=(10274, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4 + NH2 <=> C2H3 + NH3""",
)

entry(
    index = 907,
    label = "C2H3 + NO <=> HCN + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (7e+21, 'cm^3/(mol*s)'),
        n = -3.382,
        Ea = (1025, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H3 + NO <=> HCN + CH2O""",
)

entry(
    index = 908,
    label = "C2H3 + HONO <=> C2H4 + NO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(810000, 'cm^3/(mol*s)'), n=1.87, Ea=(5504, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3 + HONO <=> C2H4 + NO2""",
)

entry(
    index = 909,
    label = "C2H3 + HNO2 <=> C2H4 + NO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(810000, 'cm^3/(mol*s)'), n=1.87, Ea=(4838, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3 + HNO2 <=> C2H4 + NO2""",
)

entry(
    index = 910,
    label = "C2H3 + NO2 <=> CH2CHO + NO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.7e+14, 'cm^3/(mol*s)'), n=-0.6, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3 + NO2 <=> CH2CHO + NO""",
)

entry(
    index = 911,
    label = "C2H2 + NCO <=> HCCO + HCN",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.4e+12, 'cm^3/(mol*s)'), n=0, Ea=(1815, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H2 + NCO <=> HCCO + HCN""",
)

entry(
    index = 912,
    label = "C2H + NH3 <=> C2H2 + NH2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.2e+12, 'cm^3/(mol*s)'), n=0, Ea=(-735, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H + NH3 <=> C2H2 + NH2""",
)

entry(
    index = 913,
    label = "C2H + NO <=> HCN + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(570, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H + NO <=> HCN + CO""",
)

entry(
    index = 914,
    label = "C2H + NO <=> CN + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(570, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H + NO <=> CN + HCO""",
)

entry(
    index = 915,
    label = "C2H + NO2 <=> HCCO + NO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(-258, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H + NO2 <=> HCCO + NO""",
)

entry(
    index = 916,
    label = "C2 + N2 <=> CN + CN",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.5e+13, 'cm^3/(mol*s)'), n=0, Ea=(41730, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2 + N2 <=> CN + CN""",
)

entry(
    index = 917,
    label = "C2 + NO <=> C2O + N",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.3e+13, 'cm^3/(mol*s)'), n=0, Ea=(8640, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2 + NO <=> C2O + N""",
)

entry(
    index = 918,
    label = "CH3CH2O + NO <=> CH3CHO + HNO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.6e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2O + NO <=> CH3CHO + HNO""",
)

entry(
    index = 919,
    label = "CH3CH2O + NO2 <=> CH3CHO + HONO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2O + NO2 <=> CH3CHO + HONO""",
)

entry(
    index = 920,
    label = "CH2CH2OH + NO2 => CH2O + CH2OH + NO",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CH2OH + NO2 => CH2O + CH2OH + NO""",
)

entry(
    index = 921,
    label = "CH2CHO + NO <=> HCN + HOCHO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (7e+21, 'cm^3/(mol*s)'),
        n = -3.382,
        Ea = (1025, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHO + NO <=> HCN + HOCHO""",
)

entry(
    index = 922,
    label = "CH2CHO + NO2 <=> CH2CO + HONO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+15, 'cm^3/(mol*s)'), n=-0.68, Ea=(1430, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHO + NO2 <=> CH2CO + HONO""",
)

entry(
    index = 923,
    label = "CH3CO + NO2 => CH3 + CO2 + NO",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CO + NO2 => CH3 + CO2 + NO""",
)

entry(
    index = 924,
    label = "HCCO + N <=> HCN + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCCO + N <=> HCN + CO""",
)

entry(
    index = 925,
    label = "HCCO + NO <=> HCNO + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(-676, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCCO + NO <=> HCNO + CO""",
)

entry(
    index = 926,
    label = "HCCO + NO <=> HCN + CO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.1e+12, 'cm^3/(mol*s)'), n=0, Ea=(-676, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCCO + NO <=> HCN + CO2""",
)

entry(
    index = 927,
    label = "HCCO + NO2 <=> HCNO + CO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCCO + NO2 <=> HCNO + CO2""",
)

entry(
    index = 928,
    label = "C2O + NO <=> CO + NCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(670, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2O + NO <=> CO + NCO""",
)

entry(
    index = 929,
    label = "C2O + NO2 <=> CO2 + NCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.1e+13, 'cm^3/(mol*s)'), n=0, Ea=(125, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2O + NO2 <=> CO2 + NCO""",
)

entry(
    index = 930,
    label = "NCN <=> C + N2",
    degeneracy = 1,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(A=(8.9e+14, 'cm^3/(mol*s)'), n=0, Ea=(62100, 'cal/mol'), T0=(1, 'K')),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is NCN <=> C + N2""",
)

entry(
    index = 931,
    label = "H + NCN <=> HNCN",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (3.9e+23, 'cm^3/(mol*s)'),
                n = -4.34,
                Ea = (5347, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (3.3e+25, 'cm^3/(mol*s)'),
                n = -4.71,
                Ea = (4102, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (5.6e+27, 'cm^3/(mol*s)'),
                n = -5.13,
                Ea = (3741, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.7e+29, 'cm^3/(mol*s)'),
                n = -5.36,
                Ea = (3947, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.5e+30, 'cm^3/(mol*s)'),
                n = -5.43,
                Ea = (4415, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.6e+30, 'cm^3/(mol*s)'),
                n = -5.34,
                Ea = (4870, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.2e+30, 'cm^3/(mol*s)'),
                n = -5.09,
                Ea = (5275, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.9e+29, 'cm^3/(mol*s)'),
                n = -4.72,
                Ea = (5476, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (5.1e+27, 'cm^3/(mol*s)'),
                n = -4.15,
                Ea = (5370, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is H + NCN <=> HNCN""",
)

entry(
    index = 932,
    label = "NCN + H <=> HCN + N",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (2.2e+11, 'cm^3/(mol*s)'),
                n = 0.71,
                Ea = (5321, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.2e+11, 'cm^3/(mol*s)'),
                n = 0.71,
                Ea = (5321, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.2e+11, 'cm^3/(mol*s)'),
                n = 0.71,
                Ea = (5321, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.2e+11, 'cm^3/(mol*s)'),
                n = 0.71,
                Ea = (5321, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.2e+11, 'cm^3/(mol*s)'),
                n = 0.71,
                Ea = (5321, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.2e+11, 'cm^3/(mol*s)'),
                n = 0.71,
                Ea = (5321, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.2e+11, 'cm^3/(mol*s)'),
                n = 0.71,
                Ea = (5322, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(A=(2.3e+11, 'cm^3/(mol*s)'), n=0.7, Ea=(5327, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (2.5e+11, 'cm^3/(mol*s)'),
                n = 0.69,
                Ea = (5371, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is NCN + H <=> HCN + N""",
)

entry(
    index = 933,
    label = "NCN + H <=> HNC + N",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100], 'atm'),
        arrhenius = [
            Arrhenius(A=(0.00039, 'cm^3/(mol*s)'), n=4.7, Ea=(2440, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(0.00039, 'cm^3/(mol*s)'), n=4.7, Ea=(2440, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(0.0004, 'cm^3/(mol*s)'), n=4.7, Ea=(2440, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(0.0004, 'cm^3/(mol*s)'), n=4.7, Ea=(2438, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (0.00043, 'cm^3/(mol*s)'),
                n = 4.69,
                Ea = (2434, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (0.00049, 'cm^3/(mol*s)'),
                n = 4.67,
                Ea = (2423, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (0.00071, 'cm^3/(mol*s)'),
                n = 4.62,
                Ea = (2408, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(A=(0.0017, 'cm^3/(mol*s)'), n=4.52, Ea=(2622, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(0.0096, 'cm^3/(mol*s)'), n=4.32, Ea=(3641, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is NCN + H <=> HNC + N""",
)

entry(
    index = 934,
    label = "NCN + O <=> CN + NO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.5e+13, 'cm^3/(mol*s)'), n=0.17, Ea=(-34, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NCN + O <=> CN + NO""",
)

entry(
    index = 935,
    label = "NCN + OH <=> NCNOH",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.03, 0.05, 0.1, 0.3, 1, 3, 10, 30, 100], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (1.6e+31, 'cm^3/(mol*s)'),
                n = -6.65,
                Ea = (2718, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.9e+31, 'cm^3/(mol*s)'),
                n = -6.59,
                Ea = (2940, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (3.8e+31, 'cm^3/(mol*s)'),
                n = -6.55,
                Ea = (3042, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.2e+32, 'cm^3/(mol*s)'),
                n = -6.7,
                Ea = (3421, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.6e+32, 'cm^3/(mol*s)'),
                n = -6.51,
                Ea = (3578, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.8e+32, 'cm^3/(mol*s)'),
                n = -6.37,
                Ea = (3924, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (5.4e+31, 'cm^3/(mol*s)'),
                n = -6.08,
                Ea = (4106, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (6.5e+30, 'cm^3/(mol*s)'),
                n = -5.67,
                Ea = (4217, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.2e+29, 'cm^3/(mol*s)'),
                n = -5.11,
                Ea = (4086, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.5e+27, 'cm^3/(mol*s)'),
                n = -4.35,
                Ea = (3691, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is NCN + OH <=> NCNOH""",
)

entry(
    index = 936,
    label = "NCN + OH <=> NCO + NH",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.03, 0.05, 0.1, 0.3, 1, 3, 10, 30, 100], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (8.6e+14, 'cm^3/(mol*s)'),
                n = -0.95,
                Ea = (734, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.6e+15, 'cm^3/(mol*s)'),
                n = -1.08,
                Ea = (1128, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (5.4e+15, 'cm^3/(mol*s)'),
                n = -1.17,
                Ea = (1391, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.7e+16, 'cm^3/(mol*s)'),
                n = -1.3,
                Ea = (1843, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.4e+17, 'cm^3/(mol*s)'),
                n = -1.55,
                Ea = (2791, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.7e+18, 'cm^3/(mol*s)'),
                n = -1.83,
                Ea = (4143, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(A=(1e+19, 'cm^3/(mol*s)'), n=-2.03, Ea=(5607, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (2.2e+19, 'cm^3/(mol*s)'),
                n = -2.08,
                Ea = (7339, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (6.4e+18, 'cm^3/(mol*s)'),
                n = -1.88,
                Ea = (8866, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (6.3e+16, 'cm^3/(mol*s)'),
                n = -1.25,
                Ea = (10220, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is NCN + OH <=> NCO + NH""",
)

entry(
    index = 937,
    label = "NCN + OH <=> HCN + NO",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.03, 0.05, 0.1, 0.3, 1, 3, 10, 30, 100], 'atm'),
        arrhenius = [
            Arrhenius(A=(290000, 'cm^3/(mol*s)'), n=2.04, Ea=(1505, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(690000, 'cm^3/(mol*s)'), n=1.94, Ea=(1748, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (1.2e+06, 'cm^3/(mol*s)'),
                n = 1.87,
                Ea = (1902, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.9e+06, 'cm^3/(mol*s)'),
                n = 1.76,
                Ea = (2163, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.9e+07, 'cm^3/(mol*s)'),
                n = 1.54,
                Ea = (2727, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.6e+08, 'cm^3/(mol*s)'),
                n = 1.22,
                Ea = (3593, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (4.1e+09, 'cm^3/(mol*s)'),
                n = 0.89,
                Ea = (4624, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (7.1e+10, 'cm^3/(mol*s)'),
                n = 0.56,
                Ea = (5985, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (3.9e+11, 'cm^3/(mol*s)'),
                n = 0.38,
                Ea = (7329, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.6e+11, 'cm^3/(mol*s)'),
                n = 0.48,
                Ea = (8655, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is NCN + OH <=> HCN + NO""",
)

entry(
    index = 938,
    label = "NCN + O2 <=> NCO + NO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+12, 'cm^3/(mol*s)'), n=0, Ea=(23167, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NCN + O2 <=> NCO + NO""",
)

entry(
    index = 939,
    label = "NCN + NO <=> CN + N2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.9e+12, 'cm^3/(mol*s)'), n=0, Ea=(6280, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NCN + NO <=> CN + N2O""",
)

entry(
    index = 940,
    label = "NCN + NCN <=> CN + CN + N2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.7e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NCN + NCN <=> CN + CN + N2""",
)

entry(
    index = 941,
    label = "NCN + H2 <=> HNCN + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.1e+13, 'cm^3/(mol*s)'), n=0, Ea=(24100, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NCN + H2 <=> HNCN + H""",
)

entry(
    index = 942,
    label = "HNCN + O <=> HNC + NO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+14, 'cm^3/(mol*s)'), n=-0.05, Ea=(72, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HNCN + O <=> HNC + NO""",
)

entry(
    index = 943,
    label = "HNCN + O <=> NH + NCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.6e+13, 'cm^3/(mol*s)'), n=-0.05, Ea=(72, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HNCN + O <=> NH + NCO""",
)

entry(
    index = 944,
    label = "HNCN + O <=> CN + HNO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.4e+12, 'cm^3/(mol*s)'), n=-0.05, Ea=(72, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HNCN + O <=> CN + HNO""",
)

entry(
    index = 945,
    label = "HNCN + OH <=> NCN + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (100000, 'cm^3/(mol*s)'),
        n = 2.48,
        Ea = (-1887, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is HNCN + OH <=> NCN + H2O""",
)

entry(
    index = 946,
    label = "HNCN + O2 <=> NCN + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.3e+08, 'cm^3/(mol*s)'),
        n = 1.28,
        Ea = (24200, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is HNCN + O2 <=> NCN + HO2""",
)

entry(
    index = 947,
    label = "NCNOH <=> NCO + NH",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.03, 0.05, 0.1, 0.3, 1, 3, 10, 30, 100], 'atm'),
        arrhenius = [
            Arrhenius(A=(2.1e+35, 's^-1'), n=-7.73, Ea=(56420, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(4e+35, 's^-1'), n=-7.67, Ea=(56870, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(4e+35, 's^-1'), n=-7.61, Ea=(57080, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(8.5e+36, 's^-1'), n=-7.93, Ea=(57930, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(8.6e+34, 's^-1'), n=-7.2, Ea=(57900, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(5.8e+36, 's^-1'), n=-7.62, Ea=(59640, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1e+36, 's^-1'), n=-7.27, Ea=(60540, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(7.1e+34, 's^-1'), n=-6.81, Ea=(61640, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(6.2e+32, 's^-1'), n=-6.1, Ea=(62400, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.9e+29, 's^-1'), n=-4.97, Ea=(62850, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is NCNOH <=> NCO + NH""",
)

entry(
    index = 948,
    label = "NCNOH => NCNO + H",
    degeneracy = 1,
    reversible = False,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.03, 0.05, 0.1, 0.3, 1, 3, 10, 30, 100], 'atm'),
        arrhenius = [
            Arrhenius(A=(9.9e-28, 's^-1'), n=8.75, Ea=(50680, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1500, 's^-1'), n=0.49, Ea=(71770, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.6e-16, 's^-1'), n=5.85, Ea=(56960, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(0.00016, 's^-1'), n=2.58, Ea=(64160, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(0.93, 's^-1'), n=1.66, Ea=(64190, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2.1e+10, 's^-1'), n=-1.12, Ea=(66840, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.5e+15, 's^-1'), n=-2.31, Ea=(67060, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(5.9e+20, 's^-1'), n=-3.63, Ea=(68410, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.7e+24, 's^-1'), n=-4.35, Ea=(69140, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(9.2e+26, 's^-1'), n=-4.81, Ea=(69960, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is NCNOH => NCNO + H""",
)

entry(
    index = 949,
    label = "NCNOH <=> HCN + NO",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.03, 0.05, 0.1, 0.3, 1, 3, 10, 30, 100], 'atm'),
        arrhenius = [
            Arrhenius(A=(1.5e+23, 's^-1'), n=-4.81, Ea=(52570, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.2e+24, 's^-1'), n=-4.88, Ea=(53810, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(5.4e+24, 's^-1'), n=-4.97, Ea=(54490, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.2e+27, 's^-1'), n=-5.53, Ea=(55960, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.7e+27, 's^-1'), n=-5.35, Ea=(56960, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(4.2e+30, 's^-1'), n=-6.14, Ea=(59260, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2.1e+31, 's^-1'), n=-6.14, Ea=(60460, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(3.4e+31, 's^-1'), n=-5.99, Ea=(61650, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(5.3e+30, 's^-1'), n=-5.57, Ea=(62460, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(4.3e+28, 's^-1'), n=-4.78, Ea=(62950, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is NCNOH <=> HCN + NO""",
)

entry(
    index = 950,
    label = "NCNOH + H <=> HNCN + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NCNOH + H <=> HNCN + OH""",
)

entry(
    index = 951,
    label = "NCNOH + O => NCNO + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.2e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NCNOH + O => NCNO + OH""",
)

entry(
    index = 952,
    label = "NCNOH + O <=> CN + HONO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NCNOH + O <=> CN + HONO""",
)

entry(
    index = 953,
    label = "NCNOH + OH => NCNO + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NCNOH + OH => NCNO + H2O""",
)

entry(
    index = 954,
    label = "NCNO + H => HNC + NO",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NCNO + H => HNC + NO""",
)

entry(
    index = 955,
    label = "NCNO + O => NCO + NO",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NCNO + O => NCO + NO""",
)

entry(
    index = 956,
    label = "HNCNH + H <=> HNCN + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.8e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(7322, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HNCNH + H <=> HNCN + H2""",
)

entry(
    index = 957,
    label = "HNCNH + O <=> HNCN + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.4e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(4630, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HNCNH + O <=> HNCN + OH""",
)

entry(
    index = 958,
    label = "HNCNH + OH <=> HNCN + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+06, 'cm^3/(mol*s)'), n=2, Ea=(-89, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HNCNH + OH <=> HNCN + H2O""",
)

entry(
    index = 959,
    label = "CH3 + NO <=> CH3NO",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(9e+12, 'cm^3/(mol*s)'), n=0, Ea=(192, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (2.5e+16, 'cm^6/(mol^2*s)'),
            n = 0,
            Ea = (-2841, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 5,
        T3 = (1e-30, 'K'),
        T1 = (120, 'K'),
        T2 = (1e+30, 'K'),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is CH3 + NO <=> CH3NO""",
)

entry(
    index = 960,
    label = "CH3NO + H <=> CH2NO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.4e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(378, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NO + H <=> CH2NO + H2""",
)

entry(
    index = 961,
    label = "CH3NO + O <=> CH2NO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.3e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(3616, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NO + O <=> CH2NO + OH""",
)

entry(
    index = 962,
    label = "CH3NO + OH <=> CH2NO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.6e+06, 'cm^3/(mol*s)'), n=2, Ea=(-1192, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NO + OH <=> CH2NO + H2O""",
)

entry(
    index = 963,
    label = "CH3NO + CH3 <=> CH2NO + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(790000, 'cm^3/(mol*s)'), n=1.87, Ea=(5415, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NO + CH3 <=> CH2NO + CH4""",
)

entry(
    index = 964,
    label = "CH3NO + NH2 <=> CH2NO + NH3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.8e+06, 'cm^3/(mol*s)'),
        n = 1.94,
        Ea = (1073, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3NO + NH2 <=> CH2NO + NH3""",
)

entry(
    index = 965,
    label = "CH3NO + O <=> CH3 + NO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+06, 'cm^3/(mol*s)'), n=2.08, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NO + O <=> CH3 + NO2""",
)

entry(
    index = 966,
    label = "CH3NO + OH <=> CH3 + HONO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(994, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NO + OH <=> CH3 + HONO""",
)

entry(
    index = 967,
    label = "CH2NO <=> HNCO + H",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.1, 1, 10], 'atm'),
        arrhenius = [
            Arrhenius(A=(6.9e+41, 's^-1'), n=-9.3, Ea=(51702, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2.3e+42, 's^-1'), n=-9.11, Ea=(53838, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.7e+38, 's^-1'), n=-7.64, Ea=(53579, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2NO <=> HNCO + H""",
)

entry(
    index = 968,
    label = "CH2NO + H <=> CH3 + NO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2NO + H <=> CH3 + NO""",
)

entry(
    index = 969,
    label = "CH2NO + H <=> HCNO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.8e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(-894, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2NO + H <=> HCNO + H2""",
)

entry(
    index = 970,
    label = "CH2NO + O <=> CH2O + NO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2NO + O <=> CH2O + NO""",
)

entry(
    index = 971,
    label = "CH2NO + O <=> HCNO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.3e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(-894, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2NO + O <=> HCNO + OH""",
)

entry(
    index = 972,
    label = "CH2NO + OH <=> CH2OH + NO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2NO + OH <=> CH2OH + NO""",
)

entry(
    index = 973,
    label = "CH2NO + OH <=> HCNO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+06, 'cm^3/(mol*s)'), n=2, Ea=(-1192, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2NO + OH <=> HCNO + H2O""",
)

entry(
    index = 974,
    label = "CH2NO + O2 <=> CH2O + NO2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.1e+23, 'cm^3/(mol*s)'),
        n = -3.29,
        Ea = (3895, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2NO + O2 <=> CH2O + NO2""",
)

entry(
    index = 975,
    label = "CH2NO + CH3 <=> C2H5 + NO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2NO + CH3 <=> C2H5 + NO""",
)

entry(
    index = 976,
    label = "CH2NO + CH3 <=> HCNO + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.6e+06, 'cm^3/(mol*s)'),
        n = 1.87,
        Ea = (-1113, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2NO + CH3 <=> HCNO + CH4""",
)

entry(
    index = 977,
    label = "CH2NO + NH2 <=> CH2NH2 + NO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2NO + NH2 <=> CH2NH2 + NO""",
)

entry(
    index = 978,
    label = "CH2NO + NH2 <=> HCNO + NH3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.8e+06, 'cm^3/(mol*s)'),
        n = 1.94,
        Ea = (-1152, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2NO + NH2 <=> HCNO + NH3""",
)

entry(
    index = 979,
    label = "C2H5 + NO <=> C2H5NO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.4e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5 + NO <=> C2H5NO""",
)

entry(
    index = 980,
    label = "C2H5NO + H => C2H4 + NO + H2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(9.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(9220, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5NO + H => C2H4 + NO + H2""",
)

entry(
    index = 981,
    label = "C2H5NO + H <=> CH3CHNO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.4e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(378, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5NO + H <=> CH3CHNO + H2""",
)

entry(
    index = 982,
    label = "C2H5NO + O => C2H4 + NO + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.1e-07, 'cm^3/(mol*s)'), n=6.5, Ea=(274, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5NO + O => C2H4 + NO + OH""",
)

entry(
    index = 983,
    label = "C2H5NO + O <=> CH3CHNO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.3e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(3616, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5NO + O <=> CH3CHNO + OH""",
)

entry(
    index = 984,
    label = "C2H5NO + OH => C2H4 + NO + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(9.2e+06, 'cm^3/(mol*s)'), n=2, Ea=(990, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5NO + OH => C2H4 + NO + H2O""",
)

entry(
    index = 985,
    label = "C2H5NO + OH <=> CH3CHNO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.6e+06, 'cm^3/(mol*s)'), n=2, Ea=(-1192, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5NO + OH <=> CH3CHNO + H2O""",
)

entry(
    index = 986,
    label = "C2H5NO + HO2 => C2H4 + NO + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(110000, 'cm^3/(mol*s)'), n=2.5, Ea=(16850, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5NO + HO2 => C2H4 + NO + H2O2""",
)

entry(
    index = 987,
    label = "C2H5NO + O2 => C2H4 + NO + HO2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(730000, 'cm^3/(mol*s)'), n=2.5, Ea=(49160, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5NO + O2 => C2H4 + NO + HO2""",
)

entry(
    index = 988,
    label = "CH3CHNO <=> C2H4 + NO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 's^-1'), n=0, Ea=(22000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHNO <=> C2H4 + NO""",
)

entry(
    index = 989,
    label = "CH2CHNO + H <=> CHCHNO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.4e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(378, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHNO + H <=> CHCHNO + H2""",
)

entry(
    index = 990,
    label = "CH2CHNO + H <=> C2H3 + HNO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(2782, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHNO + H <=> C2H3 + HNO""",
)

entry(
    index = 991,
    label = "CH2CHNO + O <=> CHCHNO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.3e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(3616, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHNO + O <=> CHCHNO + OH""",
)

entry(
    index = 992,
    label = "CH2CHNO + O <=> C2H3 + NO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+08, 'cm^3/(mol*s)'), n=2.08, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHNO + O <=> C2H3 + NO2""",
)

entry(
    index = 993,
    label = "CH2CHNO + OH <=> CHCHNO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.6e+06, 'cm^3/(mol*s)'), n=2, Ea=(-1192, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHNO + OH <=> CHCHNO + H2O""",
)

entry(
    index = 994,
    label = "CH2CHNO + OH <=> C2H3 + HONO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(994, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHNO + OH <=> C2H3 + HONO""",
)

entry(
    index = 995,
    label = "CHCHNO <=> C2H2 + NO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 's^-1'), n=0, Ea=(890, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CHCHNO <=> C2H2 + NO""",
)

entry(
    index = 996,
    label = "CH3NO2 <=> CH3 + NO2",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.8e+16, 's^-1'), n=0, Ea=(58500, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.259e+17, 'cm^3/(mol*s)'),
            n = 0,
            Ea = (42000, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.183,
        T3 = (1e-30, 'K'),
        T1 = (1e+30, 'K'),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is CH3NO2 <=> CH3 + NO2""",
)

entry(
    index = 997,
    label = "CH3NO2 + H <=> CH3 + HNO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.3e+12, 'cm^3/(mol*s)'), n=0, Ea=(3730, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NO2 + H <=> CH3 + HNO2""",
)

entry(
    index = 998,
    label = "CH3NO2 + H <=> CH3NO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.4e+12, 'cm^3/(mol*s)'), n=0, Ea=(3730, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NO2 + H <=> CH3NO + OH""",
)

entry(
    index = 999,
    label = "CH3NO2 + H <=> CH2NO2 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.9e+13, 'cm^3/(mol*s)'), n=0, Ea=(9220, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NO2 + H <=> CH2NO2 + H2""",
)

entry(
    index = 1000,
    label = "CH3NO2 + O <=> CH2NO2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.5e+13, 'cm^3/(mol*s)'), n=0, Ea=(5350, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NO2 + O <=> CH2NO2 + OH""",
)

entry(
    index = 1001,
    label = "CH3NO2 + OH <=> CH2NO2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(500000, 'cm^3/(mol*s)'), n=2, Ea=(1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NO2 + OH <=> CH2NO2 + H2O""",
)

entry(
    index = 1002,
    label = "CH3NO2 + HO2 <=> CH2NO2 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+12, 'cm^3/(mol*s)'), n=0, Ea=(23000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NO2 + HO2 <=> CH2NO2 + H2O2""",
)

entry(
    index = 1003,
    label = "CH3NO2 + O2 <=> CH2NO2 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(57000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NO2 + O2 <=> CH2NO2 + HO2""",
)

entry(
    index = 1004,
    label = "CH3NO2 + CH3 <=> CH2NO2 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.55, 'cm^3/(mol*s)'), n=4, Ea=(8300, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NO2 + CH3 <=> CH2NO2 + CH4""",
)

entry(
    index = 1005,
    label = "CH3NO2 + CH3O <=> CH2NO2 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+11, 'cm^3/(mol*s)'), n=0, Ea=(7000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NO2 + CH3O <=> CH2NO2 + CH3OH""",
)

entry(
    index = 1006,
    label = "CH3NO2 + NO2 <=> CH2NO2 + HONO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+11, 'cm^3/(mol*s)'), n=0, Ea=(32000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NO2 + NO2 <=> CH2NO2 + HONO""",
)

entry(
    index = 1007,
    label = "CH2NO2 <=> CH2O + NO",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.05, 1], 'atm'),
        arrhenius = [
            Arrhenius(A=(5e+11, 's^-1'), n=0, Ea=(36000, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1e+13, 's^-1'), n=0, Ea=(36000, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2NO2 <=> CH2O + NO""",
)

entry(
    index = 1008,
    label = "CH2NO2 + H <=> CH3 + NO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2NO2 + H <=> CH3 + NO2""",
)

entry(
    index = 1009,
    label = "CH2NO2 + O <=> CH2O + NO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2NO2 + O <=> CH2O + NO2""",
)

entry(
    index = 1010,
    label = "CH2NO2 + OH <=> CH2OH + NO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2NO2 + OH <=> CH2OH + NO2""",
)

entry(
    index = 1011,
    label = "CH2NO2 + OH <=> CH2O + HONO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2NO2 + OH <=> CH2O + HONO""",
)

entry(
    index = 1012,
    label = "CH2NO2 + CH3 <=> C2H5 + NO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2NO2 + CH3 <=> C2H5 + NO2""",
)

entry(
    index = 1013,
    label = "C2H5NO2 <=> C2H5 + NO2",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.046, 0.1, 1, 10], 'atm'),
        arrhenius = [
            Arrhenius(A=(3e+62, 's^-1'), n=-15.03, Ea=(71312, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(8.9e+64, 's^-1'), n=-15.52, Ea=(73513, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(5.6e+65, 's^-1'), n=-15.64, Ea=(74502, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2.2e+66, 's^-1'), n=-15.49, Ea=(76756, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(8.6e+63, 's^-1'), n=-14.48, Ea=(77543, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H5NO2 <=> C2H5 + NO2""",
)

entry(
    index = 1014,
    label = "C2H5NO2 <=> C2H4 + HNO2",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.046, 0.1, 1, 10], 'atm'),
        arrhenius = [
            Arrhenius(A=(2.8e+80, 's^-1'), n=-20.74, Ea=(74131, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.1e+77, 's^-1'), n=-19.55, Ea=(73632, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(9.7e+74, 's^-1'), n=-18.87, Ea=(73231, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(6.9e+67, 's^-1'), n=-16.52, Ea=(71461, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1e+59, 's^-1'), n=-13.71, Ea=(68719, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H5NO2 <=> C2H4 + HNO2""",
)

entry(
    index = 1015,
    label = "C2H5NO2 <=> CH2CHNO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.9e+51, 's^-1'), n=-20.019, Ea=(92377, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5NO2 <=> CH2CHNO + H2O""",
)

entry(
    index = 1016,
    label = "C2H5NO2 <=> CH3CH2ONO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.1e+10, 's^-1'), n=1, Ea=(60660, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5NO2 <=> CH3CH2ONO""",
)

entry(
    index = 1017,
    label = "C2H5NO2 + H <=> CH2CH2NO2 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(9220, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5NO2 + H <=> CH2CH2NO2 + H2""",
)

entry(
    index = 1018,
    label = "C2H5NO2 + H <=> CH3CHNO2 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.6e+07, 'cm^3/(mol*s)'),
        n = 1.65,
        Ea = (2827, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H5NO2 + H <=> CH3CHNO2 + H2""",
)

entry(
    index = 1019,
    label = "C2H5NO2 + H <=> C2H5NO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.4e+12, 'cm^3/(mol*s)'), n=0, Ea=(3730, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5NO2 + H <=> C2H5NO + OH""",
)

entry(
    index = 1020,
    label = "C2H5NO2 + O <=> CH2CH2NO2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.1e-07, 'cm^3/(mol*s)'), n=6.5, Ea=(274, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5NO2 + O <=> CH2CH2NO2 + OH""",
)

entry(
    index = 1021,
    label = "C2H5NO2 + O <=> CH3CHNO2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.9e+07, 'cm^3/(mol*s)'),
        n = 1.85,
        Ea = (1824, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H5NO2 + O <=> CH3CHNO2 + OH""",
)

entry(
    index = 1022,
    label = "C2H5NO2 + OH <=> CH2CH2NO2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.2e+06, 'cm^3/(mol*s)'), n=2, Ea=(990, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5NO2 + OH <=> CH2CH2NO2 + H2O""",
)

entry(
    index = 1023,
    label = "C2H5NO2 + OH <=> CH3CHNO2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.6e+11, 'cm^3/(mol*s)'), n=0.15, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5NO2 + OH <=> CH3CHNO2 + H2O""",
)

entry(
    index = 1024,
    label = "C2H5NO2 + OH <=> CH3CH2OH + NO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+10, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5NO2 + OH <=> CH3CH2OH + NO2""",
)

entry(
    index = 1025,
    label = "C2H5NO2 + HO2 <=> CH2CH2NO2 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(110000, 'cm^3/(mol*s)'), n=2.5, Ea=(16850, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5NO2 + HO2 <=> CH2CH2NO2 + H2O2""",
)

entry(
    index = 1026,
    label = "C2H5NO2 + HO2 <=> CH3CHNO2 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8200, 'cm^3/(mol*s)'), n=2.55, Ea=(10750, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5NO2 + HO2 <=> CH3CHNO2 + H2O2""",
)

entry(
    index = 1027,
    label = "C2H5NO2 + O2 <=> CH2CH2NO2 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(730000, 'cm^3/(mol*s)'), n=2.5, Ea=(49160, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5NO2 + O2 <=> CH2CH2NO2 + HO2""",
)

entry(
    index = 1028,
    label = "C2H5NO2 + CH3 <=> CH2CH2NO2 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(220, 'cm^3/(mol*s)'), n=3.18, Ea=(9622, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5NO2 + CH3 <=> CH2CH2NO2 + CH4""",
)

entry(
    index = 1029,
    label = "C2H5NO2 + CH3 <=> CH3CHNO2 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(730, 'cm^3/(mol*s)'), n=2.99, Ea=(7948, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5NO2 + CH3 <=> CH3CHNO2 + CH4""",
)

entry(
    index = 1030,
    label = "C2H4 + NO2 <=> CH2CH2NO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.4e+11, 'cm^3/(mol*s)'), n=0, Ea=(14000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4 + NO2 <=> CH2CH2NO2""",
)

entry(
    index = 1031,
    label = "CH3CHNO2 <=> C2H4 + NO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 's^-1'), n=0, Ea=(22000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHNO2 <=> C2H4 + NO2""",
)

entry(
    index = 1032,
    label = "CH3O + NO <=> CH3ONO",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(6e+14, 'cm^3/(mol*s)'), n=-0.6, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (8.14e+25, 'cm^6/(mol^2*s)'),
            n = -2.8,
            Ea = (0, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 1,
        T3 = (1e-30, 'K'),
        T1 = (900, 'K'),
        T2 = (1e+30, 'K'),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is CH3O + NO <=> CH3ONO""",
)

entry(
    index = 1033,
    label = "CH3ONO + H <=> CH3OH + NO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+11, 'cm^3/(mol*s)'), n=0, Ea=(1900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3ONO + H <=> CH3OH + NO""",
)

entry(
    index = 1034,
    label = "CH3ONO + H <=> CH2O + H2 + NO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.4e+11, 'cm^3/(mol*s)'), n=0, Ea=(1900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3ONO + H <=> CH2O + H2 + NO""",
)

entry(
    index = 1035,
    label = "CH3ONO + O <=> CH3O + NO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(5210, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3ONO + O <=> CH3O + NO2""",
)

entry(
    index = 1036,
    label = "CH3ONO + OH <=> CH3OH + NO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+13, 'cm^3/(mol*s)'), n=0, Ea=(3505, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3ONO + OH <=> CH3OH + NO2""",
)

entry(
    index = 1037,
    label = "CH3CH2O + NO <=> CH3CH2ONO",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.2e+13, 'cm^3/(mol*s)'), n=0, Ea=(-143, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(A=(9.43e+19, 'cm^6/(mol^2*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        alpha = 0.6,
        T3 = (1e-30, 'K'),
        T1 = (1e+30, 'K'),
        T2 = (1e+30, 'K'),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is CH3CH2O + NO <=> CH3CH2ONO""",
)

entry(
    index = 1038,
    label = "CH3CH2ONO + OH <=> CH3CH2OH + NO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.6e+14, 'cm^3/(mol*s)'), n=0, Ea=(3505, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2ONO + OH <=> CH3CH2OH + NO2""",
)

entry(
    index = 1039,
    label = "CH3O + NO2 <=> CH3ONO2",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2.2e+15, 'cm^3/(mol*s)'), n=-0.88, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (2.911e+23, 'cm^6/(mol^2*s)'),
            n = -1.74,
            Ea = (0, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.6,
        T3 = (1e-30, 'K'),
        T1 = (1e+30, 'K'),
        T2 = (1e+30, 'K'),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is CH3O + NO2 <=> CH3ONO2""",
)

entry(
    index = 1040,
    label = "CH3ONO2 + H <=> CH3O + HONO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3ONO2 + H <=> CH3O + HONO""",
)

entry(
    index = 1041,
    label = "CH3ONO2 + O <=> CH3O + NO3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.5e+13, 'cm^3/(mol*s)'), n=0, Ea=(5260, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3ONO2 + O <=> CH3O + NO3""",
)

entry(
    index = 1042,
    label = "CH3ONO2 + OH <=> CH3O + HONO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.9e+11, 'cm^3/(mol*s)'), n=0, Ea=(2027, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3ONO2 + OH <=> CH3O + HONO2""",
)

entry(
    index = 1043,
    label = "CH3CH2O + NO2 <=> CH3CH2ONO2",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(5.1e+15, 'cm^3/(mol*s)'), n=-1, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(A=(5.88e+30, 'cm^6/(mol^2*s)'), n=-4, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        alpha = 0.6,
        T3 = (1e-30, 'K'),
        T1 = (1e+30, 'K'),
        T2 = (1e+30, 'K'),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is CH3CH2O + NO2 <=> CH3CH2ONO2""",
)

entry(
    index = 1044,
    label = "CH3CH2ONO2 + OH <=> CH3CH2O + HONO2",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(2.2e+12, 'cm^3/(mol*s)'), n=0, Ea=(2140, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(3.2e+10, 'cm^3/(mol*s)'), n=0, Ea=(-250, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH3CH2ONO2 + OH <=> CH3CH2O + HONO2""",
)

entry(
    index = 1045,
    label = "CH3 + NH2 <=> CH3NH2",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.1, 1, 10], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (1.3e+54, 'cm^3/(mol*s)'),
                n = -12.72,
                Ea = (15608, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (5.1e+52, 'cm^3/(mol*s)'),
                n = -11.99,
                Ea = (16790, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.6e+47, 'cm^3/(mol*s)'),
                n = -10.15,
                Ea = (15687, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH3 + NH2 <=> CH3NH2""",
)

entry(
    index = 1046,
    label = "CH3NH2 <=> CH2NH + H2",
    degeneracy = 1,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(107260, 'cal/mol'), T0=(1, 'K')),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is CH3NH2 <=> CH2NH + H2""",
)

entry(
    index = 1047,
    label = "CH3NH2 + H <=> CH2NH2 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.6e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(5464, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NH2 + H <=> CH2NH2 + H2""",
)

entry(
    index = 1048,
    label = "CH3NH2 + H <=> CH3NH + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.8e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(9706, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NH2 + H <=> CH3NH + H2""",
)

entry(
    index = 1049,
    label = "CH3NH2 + O <=> CH2NH2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(5196, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NH2 + O <=> CH2NH2 + OH""",
)

entry(
    index = 1050,
    label = "CH3NH2 + O <=> CH3NH + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.3e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(6348, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NH2 + O <=> CH3NH + OH""",
)

entry(
    index = 1051,
    label = "CH3NH2 + OH <=> CH2NH2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NH2 + OH <=> CH2NH2 + H2O""",
)

entry(
    index = 1052,
    label = "CH3NH2 + OH <=> CH3NH + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+06, 'cm^3/(mol*s)'), n=2, Ea=(447, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NH2 + OH <=> CH3NH + H2O""",
)

entry(
    index = 1053,
    label = "CH3NH2 + CH3 <=> CH2NH2 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.5e+06, 'cm^3/(mol*s)'),
        n = 1.87,
        Ea = (9170, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3NH2 + CH3 <=> CH2NH2 + CH4""",
)

entry(
    index = 1054,
    label = "CH3NH2 + CH3 <=> CH3NH + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.6e+06, 'cm^3/(mol*s)'),
        n = 1.87,
        Ea = (8842, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3NH2 + CH3 <=> CH3NH + CH4""",
)

entry(
    index = 1055,
    label = "CH3NH2 + NH2 <=> CH2NH2 + NH3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.8e+06, 'cm^3/(mol*s)'),
        n = 1.94,
        Ea = (5494, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3NH2 + NH2 <=> CH2NH2 + NH3""",
)

entry(
    index = 1056,
    label = "CH3NH2 + NH2 <=> CH3NH + NH3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.8e+06, 'cm^3/(mol*s)'),
        n = 1.94,
        Ea = (7143, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3NH2 + NH2 <=> CH3NH + NH3""",
)

entry(
    index = 1057,
    label = "CH3 + NH2 <=> CH2NH2 + H",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.1, 1, 10], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (1.1e+13, 'cm^3/(mol*s)'),
                n = -0.13,
                Ea = (9905, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.4e+14, 'cm^3/(mol*s)'),
                n = -0.43,
                Ea = (11107, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(A=(7.4e+12, 'cm^3/(mol*s)'), n=0, Ea=(12071, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH3 + NH2 <=> CH2NH2 + H""",
)

entry(
    index = 1058,
    label = "CH3 + NH2 <=> CH3NH + H",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.1, 1, 10], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (1.2e+13, 'cm^3/(mol*s)'),
                n = -0.15,
                Ea = (16144, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (4.4e+13, 'cm^3/(mol*s)'),
                n = -0.31,
                Ea = (16641, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.4e+14, 'cm^3/(mol*s)'),
                n = -0.42,
                Ea = (17863, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH3 + NH2 <=> CH3NH + H""",
)

entry(
    index = 1059,
    label = "CH3 + NH2 <=> CH2NH + H2",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.1, 1, 10], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (2.1e+11, 'cm^3/(mol*s)'),
                n = -0.1,
                Ea = (19095, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (4.8e+11, 'cm^3/(mol*s)'),
                n = -0.2,
                Ea = (19403, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.9e+12, 'cm^3/(mol*s)'),
                n = -0.4,
                Ea = (20506, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH3 + NH2 <=> CH2NH + H2""",
)

entry(
    index = 1060,
    label = "CH2NH2 <=> CH2NH + H",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.1, 1, 10], 'atm'),
        arrhenius = [
            Arrhenius(A=(1.1e+45, 's^-1'), n=-10.24, Ea=(47817, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2.4e+48, 's^-1'), n=-10.82, Ea=(52040, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(3.2e+46, 's^-1'), n=-9.95, Ea=(53530, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2NH2 <=> CH2NH + H""",
)

entry(
    index = 1061,
    label = "CH2NH2 + H <=> CH2NH + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.8e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(-894, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2NH2 + H <=> CH2NH + H2""",
)

entry(
    index = 1062,
    label = "CH2NH2 + O <=> CH2O + NH2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2NH2 + O <=> CH2O + NH2""",
)

entry(
    index = 1063,
    label = "CH2NH2 + O <=> CH2NH + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.3e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(-894, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2NH2 + O <=> CH2NH + OH""",
)

entry(
    index = 1064,
    label = "CH2NH2 + OH <=> CH2OH + NH2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2NH2 + OH <=> CH2OH + NH2""",
)

entry(
    index = 1065,
    label = "CH2NH2 + OH <=> CH2NH + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+06, 'cm^3/(mol*s)'), n=2, Ea=(-1192, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2NH2 + OH <=> CH2NH + H2O""",
)

entry(
    index = 1066,
    label = "CH2NH2 + O2 <=> CH2NH + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+22, 'cm^3/(mol*s)'), n=-3.09, Ea=(6756, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2NH2 + O2 <=> CH2NH + HO2""",
)

entry(
    index = 1067,
    label = "CH2NH2 + CH3 <=> C2H5 + NH2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(2702, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2NH2 + CH3 <=> C2H5 + NH2""",
)

entry(
    index = 1068,
    label = "CH2NH2 + CH3 <=> CH2NH + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.6e+06, 'cm^3/(mol*s)'),
        n = 1.87,
        Ea = (-626, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2NH2 + CH3 <=> CH2NH + CH4""",
)

entry(
    index = 1069,
    label = "CH3NH <=> CH2NH + H",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.1, 1, 10], 'atm'),
        arrhenius = [
            Arrhenius(A=(1.6e+36, 's^-1'), n=-7.92, Ea=(36342, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.3e+42, 's^-1'), n=-9.24, Ea=(41340, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2.3e+44, 's^-1'), n=-9.51, Ea=(45244, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH3NH <=> CH2NH + H""",
)

entry(
    index = 1070,
    label = "CH3NH + H <=> CH2NH + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.2e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(-894, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NH + H <=> CH2NH + H2""",
)

entry(
    index = 1071,
    label = "CH3NH + O <=> CH2NH + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(-894, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NH + O <=> CH2NH + OH""",
)

entry(
    index = 1072,
    label = "CH3NH + OH <=> CH2NH + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.6e+06, 'cm^3/(mol*s)'), n=2, Ea=(-1192, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NH + OH <=> CH2NH + H2O""",
)

entry(
    index = 1073,
    label = "CH3NH + CH3 <=> CH2NH + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.4e+06, 'cm^3/(mol*s)'),
        n = 1.87,
        Ea = (-1113, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3NH + CH3 <=> CH2NH + CH4""",
)

entry(
    index = 1074,
    label = "CH2NH + H <=> H2CN + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(7322, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2NH + H <=> H2CN + H2""",
)

entry(
    index = 1075,
    label = "CH2NH + H <=> HCNH + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(6130, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2NH + H <=> HCNH + H2""",
)

entry(
    index = 1076,
    label = "CH2NH + O <=> H2CN + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(4630, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2NH + O <=> H2CN + OH""",
)

entry(
    index = 1077,
    label = "CH2NH + O <=> HCNH + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.2e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(5404, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2NH + O <=> HCNH + OH""",
)

entry(
    index = 1078,
    label = "CH2NH + O <=> CH2O + NH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+06, 'cm^3/(mol*s)'), n=2.08, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2NH + O <=> CH2O + NH""",
)

entry(
    index = 1079,
    label = "CH2NH + OH <=> H2CN + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+06, 'cm^3/(mol*s)'), n=2, Ea=(-89, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2NH + OH <=> H2CN + H2O""",
)

entry(
    index = 1080,
    label = "CH2NH + OH <=> HCNH + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+06, 'cm^3/(mol*s)'), n=2, Ea=(457, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2NH + OH <=> HCNH + H2O""",
)

entry(
    index = 1081,
    label = "CH2NH + CH3 <=> H2CN + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(820000, 'cm^3/(mol*s)'), n=1.87, Ea=(7123, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2NH + CH3 <=> H2CN + CH4""",
)

entry(
    index = 1082,
    label = "CH2NH + CH3 <=> HCNH + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(530000, 'cm^3/(mol*s)'), n=1.87, Ea=(9687, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2NH + CH3 <=> HCNH + CH4""",
)

entry(
    index = 1083,
    label = "CH2NH + NH2 <=> H2CN + NH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(920000, 'cm^3/(mol*s)'), n=1.94, Ea=(4441, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2NH + NH2 <=> H2CN + NH3""",
)

entry(
    index = 1084,
    label = "CH2NH + NH2 <=> HCNH + NH3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.8e+06, 'cm^3/(mol*s)'),
        n = 1.94,
        Ea = (6090, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2NH + NH2 <=> HCNH + NH3""",
)

entry(
    index = 1085,
    label = "H2CN <=> HCN + H",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.1, 1, 10], 'atm'),
        arrhenius = [
            Arrhenius(A=(1.3e+29, 's^-1'), n=-6.03, Ea=(29894, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(6e+31, 's^-1'), n=-6.46, Ea=(32110, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(3.5e+29, 's^-1'), n=-5.46, Ea=(32547, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is H2CN <=> HCN + H""",
)

entry(
    index = 1086,
    label = "H2CN + H <=> HCN + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(-894, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2CN + H <=> HCN + H2""",
)

entry(
    index = 1087,
    label = "H2CN + O <=> HCN + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(-894, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2CN + O <=> HCN + OH""",
)

entry(
    index = 1088,
    label = "H2CN + OH <=> HCN + H2O",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiPDepArrhenius(
        arrhenius = [
            PDepArrhenius(
                pressures = ([0.1, 1, 10], 'atm'),
                arrhenius = [
                    Arrhenius(
                        A = (2.1e+17, 'cm^3/(mol*s)'),
                        n = -1.68,
                        Ea = (318, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.5e+19, 'cm^3/(mol*s)'),
                        n = -2.18,
                        Ea = (2166, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (9.5e+21, 'cm^3/(mol*s)'),
                        n = -2.91,
                        Ea = (5633, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            PDepArrhenius(
                pressures = ([0.1, 1, 10], 'atm'),
                arrhenius = [
                    Arrhenius(A=(1.2e+06, 'cm^3/(mol*s)'), n=2, Ea=(-1192, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.2e+06, 'cm^3/(mol*s)'), n=2, Ea=(-1192, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.2e+06, 'cm^3/(mol*s)'), n=2, Ea=(-1192, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is H2CN + OH <=> HCN + H2O""",
)

entry(
    index = 1089,
    label = "H2CN + O2 <=> CH2O + NO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+12, 'cm^3/(mol*s)'), n=0, Ea=(5961, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2CN + O2 <=> CH2O + NO""",
)

entry(
    index = 1090,
    label = "H2CN + NH2 <=> HCN + NH3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (920000, 'cm^3/(mol*s)'),
        n = 1.94,
        Ea = (-1152, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is H2CN + NH2 <=> HCN + NH3""",
)

entry(
    index = 1091,
    label = "H2CN + NH <=> HCN + NH2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(-894, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2CN + NH <=> HCN + NH2""",
)

entry(
    index = 1092,
    label = "H2CN + N <=> CH2 + N2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2CN + N <=> CH2 + N2""",
)

entry(
    index = 1093,
    label = "HCNH <=> HCN + H",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.1, 1, 10], 'atm'),
        arrhenius = [
            Arrhenius(A=(7.7e+25, 's^-1'), n=-5.2, Ea=(21986, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(6.1e+28, 's^-1'), n=-5.69, Ea=(24271, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(6.2e+26, 's^-1'), n=-4.77, Ea=(24818, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is HCNH <=> HCN + H""",
)

entry(
    index = 1094,
    label = "HCNH + H <=> H2CN + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCNH + H <=> H2CN + H""",
)

entry(
    index = 1095,
    label = "HCNH + H <=> HCN + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(-894, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCNH + H <=> HCN + H2""",
)

entry(
    index = 1096,
    label = "HCNH + O <=> HNCO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCNH + O <=> HNCO + H""",
)

entry(
    index = 1097,
    label = "HCNH + O <=> HCN + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(-894, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCNH + O <=> HCN + OH""",
)

entry(
    index = 1098,
    label = "HCNH + OH <=> HCN + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+06, 'cm^3/(mol*s)'), n=2, Ea=(-1192, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCNH + OH <=> HCN + H2O""",
)

entry(
    index = 1099,
    label = "HCNH + CH3 <=> HCN + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (820000, 'cm^3/(mol*s)'),
        n = 1.87,
        Ea = (-1113, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is HCNH + CH3 <=> HCN + CH4""",
)

entry(
    index = 1100,
    label = "CH3CH2NH2 <=> C2H4 + NH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.2e+67, 's^-1'), n=-15.944, Ea=(99348, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2NH2 <=> C2H4 + NH3""",
)

entry(
    index = 1101,
    label = "C2H5 + NH2 <=> CH3CH2NH2",
    degeneracy = 1,
    kinetics = Lindemann(
        arrheniusHigh = Arrhenius(A=(7.2e+12, 'cm^3/(mol*s)'), n=0.42, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (2.2e+30, 'cm^6/(mol^2*s)'),
            n = -3.85,
            Ea = (0, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is C2H5 + NH2 <=> CH3CH2NH2""",
)

entry(
    index = 1102,
    label = "CH3CHNH2 + H <=> CH3CH2NH2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+13, 'cm^3/(mol*s)'), n=0.22, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHNH2 + H <=> CH3CH2NH2""",
)

entry(
    index = 1103,
    label = "CH2CH2NH2 + H <=> CH3CH2NH2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.4e+13, 'cm^3/(mol*s)'), n=0.16, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CH2NH2 + H <=> CH3CH2NH2""",
)

entry(
    index = 1104,
    label = "CH3CH2NH2 + H <=> CH2CH2NH2 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+07, 'cm^3/(mol*s)'), n=1.8, Ea=(5100, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2NH2 + H <=> CH2CH2NH2 + H2""",
)

entry(
    index = 1105,
    label = "CH3CH2NH2 + H <=> CH3CHNH2 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.6e+07, 'cm^3/(mol*s)'),
        n = 1.65,
        Ea = (2830, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3CH2NH2 + H <=> CH3CHNH2 + H2""",
)

entry(
    index = 1106,
    label = "CH3CH2NH2 + H <=> CH3CH2NH + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.8e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(9700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2NH2 + H <=> CH3CH2NH + H2""",
)

entry(
    index = 1107,
    label = "CH3CH2NH2 + O <=> CH2CH2NH2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.4e+07, 'cm^3/(mol*s)'), n=1.7, Ea=(5460, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2NH2 + O <=> CH2CH2NH2 + OH""",
)

entry(
    index = 1108,
    label = "CH3CH2NH2 + O <=> CH3CHNH2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(1275, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2NH2 + O <=> CH3CHNH2 + OH""",
)

entry(
    index = 1109,
    label = "CH3CH2NH2 + O <=> CH3CH2NH + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.3e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(6348, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2NH2 + O <=> CH3CH2NH + OH""",
)

entry(
    index = 1110,
    label = "CH3CH2NH2 + OH <=> CH2CH2NH2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.6e+12, 'cm^3/(mol*s)'), n=0, Ea=(1300, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2NH2 + OH <=> CH2CH2NH2 + H2O""",
)

entry(
    index = 1111,
    label = "CH3CH2NH2 + OH <=> CH3CHNH2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2NH2 + OH <=> CH3CHNH2 + H2O""",
)

entry(
    index = 1112,
    label = "CH3CH2NH2 + OH <=> CH3CH2NH + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+06, 'cm^3/(mol*s)'), n=2, Ea=(447, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2NH2 + OH <=> CH3CH2NH + H2O""",
)

entry(
    index = 1113,
    label = "CH3CH2NH2 + HO2 <=> CH2CH2NH2 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(12000, 'cm^3/(mol*s)'), n=2.55, Ea=(15750, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2NH2 + HO2 <=> CH2CH2NH2 + H2O2""",
)

entry(
    index = 1114,
    label = "CH3CH2NH2 + HO2 <=> CH3CHNH2 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8200, 'cm^3/(mol*s)'), n=2.55, Ea=(10750, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2NH2 + HO2 <=> CH3CHNH2 + H2O2""",
)

entry(
    index = 1115,
    label = "CH3CH2NH2 + CH3 <=> CH2CH2NH2 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(220, 'cm^3/(mol*s)'), n=3.18, Ea=(9620, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2NH2 + CH3 <=> CH2CH2NH2 + CH4""",
)

entry(
    index = 1116,
    label = "CH3CH2NH2 + CH3 <=> CH3CHNH2 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(730, 'cm^3/(mol*s)'), n=2.99, Ea=(7950, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2NH2 + CH3 <=> CH3CHNH2 + CH4""",
)

entry(
    index = 1117,
    label = "CH3CH2NH2 + CH3 <=> CH3CH2NH + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.6e+06, 'cm^3/(mol*s)'),
        n = 1.87,
        Ea = (8842, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3CH2NH2 + CH3 <=> CH3CH2NH + CH4""",
)

entry(
    index = 1118,
    label = "CH3CH2NH2 + NH2 <=> CH2CH2NH2 + NH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(220, 'cm^3/(mol*s)'), n=3.18, Ea=(9620, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2NH2 + NH2 <=> CH2CH2NH2 + NH3""",
)

entry(
    index = 1119,
    label = "CH3CH2NH2 + NH2 <=> CH3CHNH2 + NH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(730, 'cm^3/(mol*s)'), n=2.99, Ea=(7950, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2NH2 + NH2 <=> CH3CHNH2 + NH3""",
)

entry(
    index = 1120,
    label = "CH3CH2NH2 + NH2 <=> CH3CH2NH + NH3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.8e+06, 'cm^3/(mol*s)'),
        n = 1.94,
        Ea = (7140, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3CH2NH2 + NH2 <=> CH3CH2NH + NH3""",
)

entry(
    index = 1121,
    label = "C2H4 + NH2 <=> CH2CH2NH2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+11, 'cm^3/(mol*s)'), n=0, Ea=(3955, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4 + NH2 <=> CH2CH2NH2""",
)

entry(
    index = 1122,
    label = "CH2CH2NH2 + H <=> CH2CHNH2 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CH2NH2 + H <=> CH2CHNH2 + H2""",
)

entry(
    index = 1123,
    label = "CH2CH2NH2 + O <=> CH2O + CH2NH2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CH2NH2 + O <=> CH2O + CH2NH2""",
)

entry(
    index = 1124,
    label = "CH2CH2NH2 + OH <=> CH2CHNH2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CH2NH2 + OH <=> CH2CHNH2 + H2O""",
)

entry(
    index = 1125,
    label = "CH2CH2NH2 + HO2 => CH2O + OH + CH2NH2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CH2NH2 + HO2 => CH2O + OH + CH2NH2""",
)

entry(
    index = 1126,
    label = "CH2CH2NH2 + O2 <=> CH2CHNH2 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.7e+16, 'cm^3/(mol*s)'),
        n = -1.63,
        Ea = (3418, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2CH2NH2 + O2 <=> CH2CHNH2 + HO2""",
)

entry(
    index = 1127,
    label = "CH2CH2NH2 + HCO <=> CH3CH2NH2 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CH2NH2 + HCO <=> CH3CH2NH2 + CO""",
)

entry(
    index = 1128,
    label = "CH2CH2NH2 + CH3 <=> CH2CHNH2 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+13, 'cm^3/(mol*s)'), n=-0.32, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CH2NH2 + CH3 <=> CH2CHNH2 + CH4""",
)

entry(
    index = 1129,
    label = "CH2CHNH2 + H <=> CH3CHNH2",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(
            A = (1.4e+09, 'cm^3/(mol*s)'),
            n = 1.463,
            Ea = (1355, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        arrheniusLow = Arrhenius(
            A = (2e+39, 'cm^6/(mol^2*s)'),
            n = -6.642,
            Ea = (5769, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -0.569,
        T3 = (299, 'K'),
        T1 = (9147, 'K'),
        T2 = (152.4, 'K'),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHNH2 + H <=> CH3CHNH2""",
)

entry(
    index = 1130,
    label = "CH3CHNH2 <=> CH3CHNH + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.1e+45, 's^-1'), n=-10.24, Ea=(47817, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHNH2 <=> CH3CHNH + H""",
)

entry(
    index = 1131,
    label = "CH3CHNH2 + H <=> CH2CHNH2 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.9e+08, 'cm^3/(mol*s)'), n=1.7, Ea=(588, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHNH2 + H <=> CH2CHNH2 + H2""",
)

entry(
    index = 1132,
    label = "CH3CHNH2 + H <=> CH3 + CH2NH2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8.4e+16, 'cm^3/(mol*s)'),
        n = -0.891,
        Ea = (2903, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3CHNH2 + H <=> CH3 + CH2NH2""",
)

entry(
    index = 1133,
    label = "CH3CHNH2 + H <=> C2H4 + NH3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.7e+21, 'cm^3/(mol*s)'),
        n = -3.02,
        Ea = (2845, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3CHNH2 + H <=> C2H4 + NH3""",
)

entry(
    index = 1134,
    label = "CH3CHNH2 + H <=> C2H5 + NH2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHNH2 + H <=> C2H5 + NH2""",
)

entry(
    index = 1135,
    label = "CH3CHNH2 + O <=> CH3 + H2NCHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHNH2 + O <=> CH3 + H2NCHO""",
)

entry(
    index = 1136,
    label = "CH3CHNH2 + O <=> CH2CHNH2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHNH2 + O <=> CH2CHNH2 + OH""",
)

entry(
    index = 1137,
    label = "CH3CHNH2 + OH <=> CH2CHNH2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHNH2 + OH <=> CH2CHNH2 + H2O""",
)

entry(
    index = 1138,
    label = "CH3CHNH2 + HO2 => CH3 + OH + H2NCHO",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHNH2 + HO2 => CH3 + OH + H2NCHO""",
)

entry(
    index = 1139,
    label = "CH3CHNH2 + O2 <=> CH2CHNH2 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (6.7e+20, 'cm^3/(mol*s)'),
        n = -3.02,
        Ea = (2504, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3CHNH2 + O2 <=> CH2CHNH2 + HO2""",
)

entry(
    index = 1140,
    label = "CH3CHNH2 + HCO <=> CH3CH2NH2 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHNH2 + HCO <=> CH3CH2NH2 + CO""",
)

entry(
    index = 1141,
    label = "CH3CHNH2 + CH3 <=> CH2CHNH2 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(-769, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHNH2 + CH3 <=> CH2CHNH2 + CH4""",
)

entry(
    index = 1142,
    label = "CH3CH2NH <=> CH2NH + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.9e+10, 's^-1'), n=0, Ea=(23500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2NH <=> CH2NH + CH3""",
)

entry(
    index = 1143,
    label = "CH3CH2NH <=> CH3CHNH + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.6e+36, 's^-1'), n=-7.92, Ea=(36342, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2NH <=> CH3CHNH + H""",
)

entry(
    index = 1144,
    label = "CH3CH2NH + H <=> CH3 + CH2NH2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.4e+12, 'cm^3/(mol*s)'),
        n = 0.701,
        Ea = (346, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3CH2NH + H <=> CH3 + CH2NH2""",
)

entry(
    index = 1145,
    label = "CH3CH2NH + H <=> CH3CHNH + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.2e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(-894, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2NH + H <=> CH3CHNH + H2""",
)

entry(
    index = 1146,
    label = "CH3CH2NH + O <=> CH3CHNH + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(-894, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2NH + O <=> CH3CHNH + OH""",
)

entry(
    index = 1147,
    label = "CH3CH2NH + OH <=> CH3CHNH + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.6e+06, 'cm^3/(mol*s)'), n=2, Ea=(-1192, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2NH + OH <=> CH3CHNH + H2O""",
)

entry(
    index = 1148,
    label = "CH3CH2NH + CH3 <=> CH3CHNH + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.4e+06, 'cm^3/(mol*s)'),
        n = 1.87,
        Ea = (-1113, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3CH2NH + CH3 <=> CH3CHNH + CH4""",
)

entry(
    index = 1149,
    label = "CHCHNH2 + H <=> CH2CHNH2",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.9e+13, 'cm^3/(mol*s)'), n=0.2, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(A=(2.1e+24, 'cm^6/(mol^2*s)'), n=-1.3, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        alpha = 0.5,
        T3 = (1e-30, 'K'),
        T1 = (1e+30, 'K'),
        T2 = (1e+30, 'K'),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is CHCHNH2 + H <=> CH2CHNH2""",
)

entry(
    index = 1150,
    label = "CH2CNH2 + H <=> CH2CHNH2",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.9e+13, 'cm^3/(mol*s)'), n=0.2, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(A=(2.1e+24, 'cm^6/(mol^2*s)'), n=-1.3, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        alpha = 0.5,
        T3 = (1e-30, 'K'),
        T1 = (1e+30, 'K'),
        T2 = (1e+30, 'K'),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is CH2CNH2 + H <=> CH2CHNH2""",
)

entry(
    index = 1151,
    label = "CH2CHNH2 + H <=> CHCHNH2 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(240, 'cm^3/(mol*s)'), n=3.63, Ea=(11266, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHNH2 + H <=> CHCHNH2 + H2""",
)

entry(
    index = 1152,
    label = "CH2CHNH2 + H <=> CH2CNH2 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(240, 'cm^3/(mol*s)'), n=3.63, Ea=(11266, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHNH2 + H <=> CH2CNH2 + H2""",
)

entry(
    index = 1153,
    label = "CH2CHNH2 + H <=> CH2CHNH + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.8e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(9700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHNH2 + H <=> CH2CHNH + H2""",
)

entry(
    index = 1154,
    label = "CH2CHNH2 + O <=> CH3 + H2NCO",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(3.9e+12, 'cm^3/(mol*s)'), n=0, Ea=(1494, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(6.2e+13, 'cm^3/(mol*s)'), n=0, Ea=(6855, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHNH2 + O <=> CH3 + H2NCO""",
)

entry(
    index = 1155,
    label = "CH2CHNH2 + O <=> CH2CHNH + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.3e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(6348, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHNH2 + O <=> CH2CHNH + OH""",
)

entry(
    index = 1156,
    label = "CH2CHNH2 + OH <=> CHCHNH2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.13, 'cm^3/(mol*s)'), n=4.2, Ea=(-860, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHNH2 + OH <=> CHCHNH2 + H2O""",
)

entry(
    index = 1157,
    label = "CH2CHNH2 + OH <=> CH2CNH2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.13, 'cm^3/(mol*s)'), n=4.2, Ea=(-860, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHNH2 + OH <=> CH2CNH2 + H2O""",
)

entry(
    index = 1158,
    label = "CH2CHNH2 + OH <=> CH2CHNH + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+06, 'cm^3/(mol*s)'), n=2, Ea=(447, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHNH2 + OH <=> CH2CHNH + H2O""",
)

entry(
    index = 1159,
    label = "CH2CHNH2 + CH3 <=> CHCHNH2 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+07, 'cm^3/(mol*s)'), n=1.56, Ea=(16630, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHNH2 + CH3 <=> CHCHNH2 + CH4""",
)

entry(
    index = 1160,
    label = "CH2CHNH2 + CH3 <=> CH2CNH2 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+07, 'cm^3/(mol*s)'), n=1.56, Ea=(16630, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHNH2 + CH3 <=> CH2CNH2 + CH4""",
)

entry(
    index = 1161,
    label = "CH2CHNH2 + CH3 <=> CH2CHNH + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.6e+06, 'cm^3/(mol*s)'),
        n = 1.87,
        Ea = (8842, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHNH2 + CH3 <=> CH2CHNH + CH4""",
)

entry(
    index = 1162,
    label = "CH2CHNH2 + NH2 <=> CHCHNH2 + NH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.3e+12, 'cm^3/(mol*s)'), n=0, Ea=(10274, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHNH2 + NH2 <=> CHCHNH2 + NH3""",
)

entry(
    index = 1163,
    label = "CH2CHNH2 + NH2 <=> CH2CNH2 + NH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.3e+12, 'cm^3/(mol*s)'), n=0, Ea=(10274, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHNH2 + NH2 <=> CH2CNH2 + NH3""",
)

entry(
    index = 1164,
    label = "CH2CHNH2 + NH2 <=> CH2CHNH + NH3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.8e+06, 'cm^3/(mol*s)'),
        n = 1.94,
        Ea = (7143, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHNH2 + NH2 <=> CH2CHNH + NH3""",
)

entry(
    index = 1165,
    label = "CH3CHNH <=> CH2CHNH2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+18, 's^-1'), n=-2.4965, Ea=(67995, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHNH <=> CH2CHNH2""",
)

entry(
    index = 1166,
    label = "CH2CHNH + H <=> CH3CHNH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.8e+13, 'cm^3/(mol*s)'),
        n = 0.18,
        Ea = (-125, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHNH + H <=> CH3CHNH""",
)

entry(
    index = 1167,
    label = "CH3 + HCNH <=> CH3CHNH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + HCNH <=> CH3CHNH""",
)

entry(
    index = 1168,
    label = "CH3CHNH + H <=> CH2CHNH2 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHNH + H <=> CH2CHNH2 + H""",
)

entry(
    index = 1169,
    label = "CH3CHNH + H <=> CH3CNH + H2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.7e+13, 'cm^3/(mol*s)'),
        n = -0.35,
        Ea = (3000, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3CHNH + H <=> CH3CNH + H2""",
)

entry(
    index = 1170,
    label = "CH3CHNH + H <=> CH2CHNH + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.9e+12, 'cm^3/(mol*s)'), n=0.4, Ea=(5359, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHNH + H <=> CH2CHNH + H2""",
)

entry(
    index = 1171,
    label = "CH3CHNH + H <=> CH3CHN + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(7322, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHNH + H <=> CH3CHN + H2""",
)

entry(
    index = 1172,
    label = "CH3CHNH + O <=> CH3CNH + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.8e+18, 'cm^3/(mol*s)'),
        n = -1.9,
        Ea = (2975, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3CHNH + O <=> CH3CNH + OH""",
)

entry(
    index = 1173,
    label = "CH3CHNH + O <=> CH2CHNH + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.7e+13, 'cm^3/(mol*s)'),
        n = -0.2,
        Ea = (3556, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3CHNH + O <=> CH2CHNH + OH""",
)

entry(
    index = 1174,
    label = "CH3CHNH + O <=> CH3CHN + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(4630, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHNH + O <=> CH3CHN + OH""",
)

entry(
    index = 1175,
    label = "CH3CHNH + OH <=> CH3CNH + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.4e+11, 'cm^3/(mol*s)'),
        n = 0.3,
        Ea = (-1000, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3CHNH + OH <=> CH3CNH + H2O""",
)

entry(
    index = 1176,
    label = "CH3CHNH + OH <=> CH2CHNH + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=-0.6, Ea=(800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHNH + OH <=> CH2CHNH + H2O""",
)

entry(
    index = 1177,
    label = "CH3CHNH + OH <=> CH3CHN + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+06, 'cm^3/(mol*s)'), n=2, Ea=(-89, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHNH + OH <=> CH3CHN + H2O""",
)

entry(
    index = 1178,
    label = "CH3CHNH + CH3 <=> CH3CNH + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.9e-07, 'cm^3/(mol*s)'), n=5.8, Ea=(2200, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHNH + CH3 <=> CH3CNH + CH4""",
)

entry(
    index = 1179,
    label = "CH3CHNH + CH3 <=> CH2CHNH + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(25, 'cm^3/(mol*s)'), n=3.15, Ea=(5727, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHNH + CH3 <=> CH2CHNH + CH4""",
)

entry(
    index = 1180,
    label = "CH3CHNH + CH3 <=> CH3CHN + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(820000, 'cm^3/(mol*s)'), n=1.87, Ea=(7123, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHNH + CH3 <=> CH3CHN + CH4""",
)

entry(
    index = 1181,
    label = "CH3CHNH + NH2 <=> CH3CHN + NH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(920000, 'cm^3/(mol*s)'), n=1.94, Ea=(4441, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHNH + NH2 <=> CH3CHN + NH3""",
)

entry(
    index = 1182,
    label = "NH2 + C2H2 <=> CHCHNH2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (7.8e-18, 'cm^3/(mol*s)'),
        n = 8.31,
        Ea = (7430, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NH2 + C2H2 <=> CHCHNH2""",
)

entry(
    index = 1183,
    label = "CHCNH2 + H <=> CHCHNH2",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(
            A = (1.7e+10, 'cm^3/(mol*s)'),
            n = 1.266,
            Ea = (2709, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        arrheniusLow = Arrhenius(
            A = (6.3e+31, 'cm^6/(mol^2*s)'),
            n = -4.664,
            Ea = (3780, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.7878,
        T3 = (-10212, 'K'),
        T1 = (1e+30, 'K'),
        efficiencies = {'[C-]#[O+]': 2, '[H][H]': 2, 'O=C=O': 3, 'O': 5},
    ),
    shortDesc = u"""The chemkin file reaction is CHCNH2 + H <=> CHCHNH2""",
)

entry(
    index = 1184,
    label = "CHCHNH2 + H <=> CHCNH2 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CHCHNH2 + H <=> CHCNH2 + H2""",
)

entry(
    index = 1185,
    label = "CHCHNH2 + OH <=> CHCNH2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CHCHNH2 + OH <=> CHCNH2 + H2O""",
)

entry(
    index = 1186,
    label = "CHCHNH2 + O2 <=> OCHCHO + NH2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CHCHNH2 + O2 <=> OCHCHO + NH2""",
)

entry(
    index = 1187,
    label = "CHCHNH2 + CH3 <=> CHCNH2 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CHCHNH2 + CH3 <=> CHCNH2 + CH4""",
)

entry(
    index = 1188,
    label = "CHCNH2 + H <=> CH2CNH2",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(
            A = (1.7e+10, 'cm^3/(mol*s)'),
            n = 1.266,
            Ea = (2709, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        arrheniusLow = Arrhenius(
            A = (6.3e+31, 'cm^6/(mol^2*s)'),
            n = -4.664,
            Ea = (3780, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.7878,
        T3 = (-10212, 'K'),
        T1 = (1e+30, 'K'),
        efficiencies = {'[C-]#[O+]': 2, '[H][H]': 2, 'O=C=O': 3, 'O': 5},
    ),
    shortDesc = u"""The chemkin file reaction is CHCNH2 + H <=> CH2CNH2""",
)

entry(
    index = 1189,
    label = "CH2CNH2 + H <=> CHCNH2 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CNH2 + H <=> CHCNH2 + H2""",
)

entry(
    index = 1190,
    label = "CH2CNH2 + O <=> CH2CO + NH2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CNH2 + O <=> CH2CO + NH2""",
)

entry(
    index = 1191,
    label = "CH2CNH2 + OH <=> CHCNH2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CNH2 + OH <=> CHCNH2 + H2O""",
)

entry(
    index = 1192,
    label = "CH2CNH2 + O2 <=> OCHCHO + NH2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CNH2 + O2 <=> OCHCHO + NH2""",
)

entry(
    index = 1193,
    label = "CH2CNH2 + CH3 <=> CHCNH2 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CNH2 + CH3 <=> CHCNH2 + CH4""",
)

entry(
    index = 1194,
    label = "CH2CHNH + H <=> CH3 + HCNH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHNH + H <=> CH3 + HCNH""",
)

entry(
    index = 1195,
    label = "CH2CHNH + H <=> CH3CNH + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHNH + H <=> CH3CNH + H""",
)

entry(
    index = 1196,
    label = "CH2CHNH + H <=> CH2CNH + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHNH + H <=> CH2CNH + H2""",
)

entry(
    index = 1197,
    label = "CH2CHNH + O <=> CH2CNH + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHNH + O <=> CH2CNH + OH""",
)

entry(
    index = 1198,
    label = "CH2CHNH + OH <=> CH2CNH + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHNH + OH <=> CH2CNH + H2O""",
)

entry(
    index = 1199,
    label = "CH2CHNH + OH <=> CH2OH + HCNH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHNH + OH <=> CH2OH + HCNH""",
)

entry(
    index = 1200,
    label = "CH2CHNH + O2 <=> CH2O + CO + NH2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.7e+17, 'cm^3/(mol*s)'),
        n = -1.757,
        Ea = (11067, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHNH + O2 <=> CH2O + CO + NH2""",
)

entry(
    index = 1201,
    label = "CH3CNH <=> CH3 + HNC",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.5e+18, 's^-1'), n=-2.52, Ea=(33000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CNH <=> CH3 + HNC""",
)

entry(
    index = 1202,
    label = "CH3CNH <=> CH3CN + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.7e+25, 's^-1'), n=-5.2, Ea=(24000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CNH <=> CH3CN + H""",
)

entry(
    index = 1203,
    label = "CH3CNH + H <=> CH3 + HCNH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CNH + H <=> CH3 + HCNH""",
)

entry(
    index = 1204,
    label = "CH3CNH + H <=> CH2CNH + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CNH + H <=> CH2CNH + H2""",
)

entry(
    index = 1205,
    label = "CH3CNH + H <=> CH3CN + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(-894, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CNH + H <=> CH3CN + H2""",
)

entry(
    index = 1206,
    label = "CH3CNH + O <=> CH3 + HNCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.6e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CNH + O <=> CH3 + HNCO""",
)

entry(
    index = 1207,
    label = "CH3CNH + O <=> CH2CNH + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CNH + O <=> CH2CNH + OH""",
)

entry(
    index = 1208,
    label = "CH3CNH + O <=> CH3CN + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(-894, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CNH + O <=> CH3CN + OH""",
)

entry(
    index = 1209,
    label = "CH3CNH + OH <=> CH2CNH + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CNH + OH <=> CH2CNH + H2O""",
)

entry(
    index = 1210,
    label = "CH3CNH + OH <=> CH3CN + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+06, 'cm^3/(mol*s)'), n=2, Ea=(-1192, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CNH + OH <=> CH3CN + H2O""",
)

entry(
    index = 1211,
    label = "CH3CNH + O2 <=> CH2O + CO + NH2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.9e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CNH + O2 <=> CH2O + CO + NH2""",
)

entry(
    index = 1212,
    label = "CH3CNH + CH3 <=> CH2CNH + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CNH + CH3 <=> CH2CNH + CH4""",
)

entry(
    index = 1213,
    label = "CH3CNH + CH3 <=> CH3CN + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (820000, 'cm^3/(mol*s)'),
        n = 1.87,
        Ea = (-1113, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3CNH + CH3 <=> CH3CN + CH4""",
)

entry(
    index = 1214,
    label = "CH3 + HCN <=> CH3CHN",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(9900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + HCN <=> CH3CHN""",
)

entry(
    index = 1215,
    label = "CH3CHN + H <=> CH3CN + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(-894, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHN + H <=> CH3CN + H2""",
)

entry(
    index = 1216,
    label = "CH3CHN + H <=> CH2CHN + H2",
    degeneracy = 1,
    duplicate = True,
    kinetics = Arrhenius(A=(9e+13, 'cm^3/(mol*s)'), n=0, Ea=(15100, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHN + H <=> CH2CHN + H2""",
)

entry(
    index = 1217,
    label = "CH2CHN(S) + H2 <=> CH3CHN + H",
    degeneracy = 1,
    duplicate = True,
    kinetics = Arrhenius(A=(7.2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHN(S) + H2 <=> CH3CHN + H""",
)

entry(
    index = 1218,
    label = "CH3CHN + O <=> CH3CN + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(-894, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHN + O <=> CH3CN + OH""",
)

entry(
    index = 1219,
    label = "CH3CHN + OH <=> CH3CN + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+06, 'cm^3/(mol*s)'), n=2, Ea=(-1192, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHN + OH <=> CH3CN + H2O""",
)

entry(
    index = 1220,
    label = "CH3CHN + OH <=> CH2CHN + H2O",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(1100, 'cm^3/(mol*s)'), n=3, Ea=(2780, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (4.4e+13, 'cm^3/(mol*s)'),
                n = -0.3485,
                Ea = (-727, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH3CHN + OH <=> CH2CHN + H2O""",
)

entry(
    index = 1221,
    label = "CH3CHN + NH2 <=> CH3CN + NH3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (920000, 'cm^3/(mol*s)'),
        n = 1.94,
        Ea = (-1152, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3CHN + NH2 <=> CH3CN + NH3""",
)

entry(
    index = 1222,
    label = "CHCNH2 + H <=> CHCNH + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.8e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(9706, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CHCNH2 + H <=> CHCNH + H2""",
)

entry(
    index = 1223,
    label = "CHCNH2 + O <=> CHCNH + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.3e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(6348, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CHCNH2 + O <=> CHCNH + OH""",
)

entry(
    index = 1224,
    label = "CHCNH2 + O <=> HCCO + NH2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.4e+07, 'cm^3/(mol*s)'), n=2, Ea=(1900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CHCNH2 + O <=> HCCO + NH2""",
)

entry(
    index = 1225,
    label = "CHCNH2 + OH <=> CHCNH + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CHCNH2 + OH <=> CHCNH + H2O""",
)

entry(
    index = 1226,
    label = "CHCNH2 + CH3 <=> CHCNH + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.6e+06, 'cm^3/(mol*s)'),
        n = 1.87,
        Ea = (8842, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CHCNH2 + CH3 <=> CHCNH + CH4""",
)

entry(
    index = 1227,
    label = "CHCNH2 + NH2 <=> CHCNH + NH3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.8e+06, 'cm^3/(mol*s)'),
        n = 1.94,
        Ea = (7143, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CHCNH2 + NH2 <=> CHCNH + NH3""",
)

entry(
    index = 1228,
    label = "CH2CNH <=> CH3CN",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.5e+13, 's^-1'), n=0, Ea=(70300, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CNH <=> CH3CN""",
)

entry(
    index = 1229,
    label = "CH2CNH + H <=> CH3CN + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CNH + H <=> CH3CN + H""",
)

entry(
    index = 1230,
    label = "CH2CNH + H <=> CH3 + HNC",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.3e+10, 'cm^3/(mol*s)'),
        n = 0.851,
        Ea = (2840, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2CNH + H <=> CH3 + HNC""",
)

entry(
    index = 1231,
    label = "CH2CNH + H <=> CHCNH + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+07, 'cm^3/(mol*s)'), n=2, Ea=(10000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CNH + H <=> CHCNH + H2""",
)

entry(
    index = 1232,
    label = "CH2CNH + H <=> CH2CN + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(7322, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CNH + H <=> CH2CN + H2""",
)

entry(
    index = 1233,
    label = "CH2CNH + O <=> CH2 + HNCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(1350, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CNH + O <=> CH2 + HNCO""",
)

entry(
    index = 1234,
    label = "CH2CNH + O <=> CHCNH + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+07, 'cm^3/(mol*s)'), n=2, Ea=(10000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CNH + O <=> CHCNH + OH""",
)

entry(
    index = 1235,
    label = "CH2CNH + O <=> CH2CN + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(4630, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CNH + O <=> CH2CN + OH""",
)

entry(
    index = 1236,
    label = "CH2CNH + OH <=> CH2OH + HNC",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1013, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CNH + OH <=> CH2OH + HNC""",
)

entry(
    index = 1237,
    label = "CH2CNH + OH <=> CH3 + HNCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.7e+11, 'cm^3/(mol*s)'), n=0, Ea=(-1013, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CNH + OH <=> CH3 + HNCO""",
)

entry(
    index = 1238,
    label = "CH2CNH + OH <=> CHCNH + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+07, 'cm^3/(mol*s)'), n=2, Ea=(3000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CNH + OH <=> CHCNH + H2O""",
)

entry(
    index = 1239,
    label = "CH2CNH + OH <=> CH2CN + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+06, 'cm^3/(mol*s)'), n=2, Ea=(-89, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CNH + OH <=> CH2CN + H2O""",
)

entry(
    index = 1240,
    label = "CH2CNH + CH3 <=> CH2CN + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(820000, 'cm^3/(mol*s)'), n=1.87, Ea=(7123, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CNH + CH3 <=> CH2CN + CH4""",
)

entry(
    index = 1241,
    label = "CH2CNH + NH2 <=> CH2CN + NH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(920000, 'cm^3/(mol*s)'), n=1.94, Ea=(4441, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CNH + NH2 <=> CH2CN + NH3""",
)

entry(
    index = 1242,
    label = "CH2CHN + H <=> CH3 + HCN",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHN + H <=> CH3 + HCN""",
)

entry(
    index = 1243,
    label = "CH2CHN + O <=> CH2O + HCN",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHN + O <=> CH2O + HCN""",
)

entry(
    index = 1244,
    label = "CH2CHN + O2 <=> CH2O + HNCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHN + O2 <=> CH2O + HNCO""",
)

entry(
    index = 1245,
    label = "CH2CHN(S) <=> CH2CHN",
    degeneracy = 1,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        efficiencies = {'[H]': 0},
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHN(S) <=> CH2CHN""",
)

entry(
    index = 1246,
    label = "CH2CHN(S) + H <=> CH2CHN + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHN(S) + H <=> CH2CHN + H""",
)

entry(
    index = 1247,
    label = "CH2CHN(S) <=> c-C2H3N",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 's^-1'), n=0, Ea=(4000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHN(S) <=> c-C2H3N""",
)

entry(
    index = 1248,
    label = "CH2CHN(S) <=> CH3CN",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 's^-1'), n=0, Ea=(8000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHN(S) <=> CH3CN""",
)

entry(
    index = 1249,
    label = "CH2CHN(S) + O => HCO + HCN + H",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHN(S) + O => HCO + HCN + H""",
)

entry(
    index = 1250,
    label = "CH2CHN(S) + OH => CH2O + HCN + H",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHN(S) + OH => CH2O + HCN + H""",
)

entry(
    index = 1251,
    label = "c-C2H3N <=> CH3CN",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.7e+13, 's^-1'), n=0, Ea=(41500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is c-C2H3N <=> CH3CN""",
)

entry(
    index = 1252,
    label = "c-C2H3N + H <=> CH2NCH2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (9.8e+09, 'cm^3/(mol*s)'),
        n = 1.212,
        Ea = (1969, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is c-C2H3N + H <=> CH2NCH2""",
)

entry(
    index = 1253,
    label = "c-C2H3N + H <=> CH2CHNH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.1e+10, 'cm^3/(mol*s)'),
        n = 1.229,
        Ea = (2422, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is c-C2H3N + H <=> CH2CHNH""",
)

entry(
    index = 1254,
    label = "c-C2H3N + O => H2CN + HCO",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is c-C2H3N + O => H2CN + HCO""",
)

entry(
    index = 1255,
    label = "c-C2H3N + O => C2H3 + NO",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is c-C2H3N + O => C2H3 + NO""",
)

entry(
    index = 1256,
    label = "c-C2H3N + OH => H2CN + CH2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is c-C2H3N + OH => H2CN + CH2O""",
)

entry(
    index = 1257,
    label = "CHCNH + H <=> CH2 + HNC",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.5e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CHCNH + H <=> CH2 + HNC""",
)

entry(
    index = 1258,
    label = "CHCNH + O <=> H + CO + HNC",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CHCNH + O <=> H + CO + HNC""",
)

entry(
    index = 1259,
    label = "CHCNH + OH <=> HCO + HCNH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CHCNH + OH <=> HCO + HCNH""",
)

entry(
    index = 1260,
    label = "CHCNH + O2 <=> HNCO + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.9e+12, 'cm^3/(mol*s)'),
        n = -0.142,
        Ea = (1150, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CHCNH + O2 <=> HNCO + HCO""",
)

entry(
    index = 1261,
    label = "CHCNH + O2 <=> HNC + CO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.6e+11, 'cm^3/(mol*s)'),
        n = -0.02,
        Ea = (1020, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CHCNH + O2 <=> HNC + CO + OH""",
)

entry(
    index = 1262,
    label = "CHCNH + O2 <=> HNC + HCO + O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(220, 'cm^3/(mol*s)'), n=2.69, Ea=(3540, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CHCNH + O2 <=> HNC + HCO + O""",
)

entry(
    index = 1263,
    label = "H2NCHO <=> CO + NH3",
    degeneracy = 1,
    kinetics = Lindemann(
        arrheniusHigh = Arrhenius(A=(1e+14, 's^-1'), n=0, Ea=(75514, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(A=(8.3e+14, 'cm^3/(mol*s)'), n=0, Ea=(49084, 'cal/mol'), T0=(1, 'K')),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is H2NCHO <=> CO + NH3""",
)

entry(
    index = 1264,
    label = "H2NCHO <=> HCO + NH2",
    degeneracy = 1,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(A=(1.4e+16, 'cm^3/(mol*s)'), n=0, Ea=(72900, 'cal/mol'), T0=(1, 'K')),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is H2NCHO <=> HCO + NH2""",
)

entry(
    index = 1265,
    label = "H2NCHO <=> H2NCO + H",
    degeneracy = 1,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(A=(4.6e+15, 'cm^3/(mol*s)'), n=0, Ea=(64200, 'cal/mol'), T0=(1, 'K')),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is H2NCHO <=> H2NCO + H""",
)

entry(
    index = 1266,
    label = "H2NCHO + H <=> H2NCO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+13, 'cm^3/(mol*s)'), n=0, Ea=(6955, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2NCHO + H <=> H2NCO + H2""",
)

entry(
    index = 1267,
    label = "H2NCHO + H <=> HCO + NH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(19100, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2NCHO + H <=> HCO + NH3""",
)

entry(
    index = 1268,
    label = "H2NCHO + O <=> H2NCO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(5196, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2NCHO + O <=> H2NCO + OH""",
)

entry(
    index = 1269,
    label = "H2NCHO + OH <=> H2NCO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2NCHO + OH <=> H2NCO + H2O""",
)

entry(
    index = 1270,
    label = "H2NCHO + CH3 <=> H2NCO + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(700000, 'cm^3/(mol*s)'), n=2, Ea=(9000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2NCHO + CH3 <=> H2NCO + CH4""",
)

entry(
    index = 1271,
    label = "H2NCHO + NH2 <=> H2NCO + NH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+06, 'cm^3/(mol*s)'), n=2, Ea=(5000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2NCHO + NH2 <=> H2NCO + NH3""",
)

entry(
    index = 1272,
    label = "H2NCO <=> CO + NH2",
    degeneracy = 1,
    kinetics = Lindemann(
        arrheniusHigh = Arrhenius(A=(5.9e+12, 's^-1'), n=0, Ea=(25000, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(21700, 'cal/mol'), T0=(1, 'K')),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is H2NCO <=> CO + NH2""",
)

entry(
    index = 1273,
    label = "H2NCO + H <=> HNCO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2NCO + H <=> HNCO + H2""",
)

entry(
    index = 1274,
    label = "H2NCO + O <=> HNCO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2NCO + O <=> HNCO + OH""",
)

entry(
    index = 1275,
    label = "H2NCO + OH <=> HNCO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2NCO + OH <=> HNCO + H2O""",
)

entry(
    index = 1276,
    label = "CH3NHCH2 + H <=> CH3NHCH3",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(
            A = (5.2e+17, 'cm^3/(mol*s)'),
            n = -0.99,
            Ea = (1580, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        arrheniusLow = Arrhenius(
            A = (1.99e+41, 'cm^6/(mol^2*s)'),
            n = -7.08,
            Ea = (6685, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.8422,
        T3 = (125, 'K'),
        T1 = (2219, 'K'),
        T2 = (6882, 'K'),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is CH3NHCH2 + H <=> CH3NHCH3""",
)

entry(
    index = 1277,
    label = "CH3NCH3 + H <=> CH3NHCH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NCH3 + H <=> CH3NHCH3""",
)

entry(
    index = 1278,
    label = "CH3NHCH3 + H <=> CH3NHCH2 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.6e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(5464, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NHCH3 + H <=> CH3NHCH2 + H2""",
)

entry(
    index = 1279,
    label = "CH3NHCH3 + H <=> CH3NCH3 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.8e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(9706, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NHCH3 + H <=> CH3NCH3 + H2""",
)

entry(
    index = 1280,
    label = "CH3NHCH3 + O <=> CH3NHCH2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.1e+12, 'cm^3/(mol*s)'), n=0, Ea=(556, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NHCH3 + O <=> CH3NHCH2 + OH""",
)

entry(
    index = 1281,
    label = "CH3NHCH3 + O <=> CH3NCH3 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+12, 'cm^3/(mol*s)'), n=0, Ea=(556, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NHCH3 + O <=> CH3NCH3 + OH""",
)

entry(
    index = 1282,
    label = "CH3NHCH3 + OH <=> CH3NHCH2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NHCH3 + OH <=> CH3NHCH2 + H2O""",
)

entry(
    index = 1283,
    label = "CH3NHCH3 + OH <=> CH3NCH3 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.9e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NHCH3 + OH <=> CH3NCH3 + H2O""",
)

entry(
    index = 1284,
    label = "CH3NHCH3 + CH3 <=> CH3NHCH2 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.5e+06, 'cm^3/(mol*s)'),
        n = 1.87,
        Ea = (9170, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3NHCH3 + CH3 <=> CH3NHCH2 + CH4""",
)

entry(
    index = 1285,
    label = "CH3NHCH3 + CH3 <=> CH3NCH3 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.6e+06, 'cm^3/(mol*s)'),
        n = 1.87,
        Ea = (8842, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3NHCH3 + CH3 <=> CH3NCH3 + CH4""",
)

entry(
    index = 1286,
    label = "CH3NHCH3 + NH2 <=> CH3NHCH2 + NH3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.8e+06, 'cm^3/(mol*s)'),
        n = 1.94,
        Ea = (5494, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3NHCH3 + NH2 <=> CH3NHCH2 + NH3""",
)

entry(
    index = 1287,
    label = "CH3NHCH3 + NH2 <=> CH3NCH3 + NH3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.8e+06, 'cm^3/(mol*s)'),
        n = 1.94,
        Ea = (7143, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3NHCH3 + NH2 <=> CH3NCH3 + NH3""",
)

entry(
    index = 1288,
    label = "CH3NHCH2 <=> CH3 + CH2NH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.8e+43, 's^-1'), n=-10.302, Ea=(37459, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NHCH2 <=> CH3 + CH2NH""",
)

entry(
    index = 1289,
    label = "CH3NHCH2 <=> CH3NCH2 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.9e+44, 's^-1'), n=-10.314, Ea=(46803, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NHCH2 <=> CH3NCH2 + H""",
)

entry(
    index = 1290,
    label = "CH3NHCH2 + H <=> CH3NCH2 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.8e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(-894, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NHCH2 + H <=> CH3NCH2 + H2""",
)

entry(
    index = 1291,
    label = "CH3NHCH2 + O <=> CH2O + CH3NH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NHCH2 + O <=> CH2O + CH3NH""",
)

entry(
    index = 1292,
    label = "CH3NHCH2 + O <=> CH3NCH2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.3e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(-894, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NHCH2 + O <=> CH3NCH2 + OH""",
)

entry(
    index = 1293,
    label = "CH3NHCH2 + OH <=> CH2OH + CH3NH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NHCH2 + OH <=> CH2OH + CH3NH""",
)

entry(
    index = 1294,
    label = "CH3NHCH2 + OH <=> CH3NCH2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+06, 'cm^3/(mol*s)'), n=2, Ea=(-1192, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NHCH2 + OH <=> CH3NCH2 + H2O""",
)

entry(
    index = 1295,
    label = "CH3NHCH2 + CH3 <=> C2H5 + CH3NH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(2702, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NHCH2 + CH3 <=> C2H5 + CH3NH""",
)

entry(
    index = 1296,
    label = "CH3NHCH2 + CH3 <=> CH3NCH2 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.6e+06, 'cm^3/(mol*s)'),
        n = 1.87,
        Ea = (-626, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3NHCH2 + CH3 <=> CH3NCH2 + CH4""",
)

entry(
    index = 1297,
    label = "CH3NCH3 <=> CH3NCH2 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.6e+35, 's^-1'), n=-7.544, Ea=(38425, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NCH3 <=> CH3NCH2 + H""",
)

entry(
    index = 1298,
    label = "CH3NCH3 + H <=> CH3NCH2 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.2e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NCH3 + H <=> CH3NCH2 + H2""",
)

entry(
    index = 1299,
    label = "CH3NCH3 + O <=> CH3NO + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NCH3 + O <=> CH3NO + CH3""",
)

entry(
    index = 1300,
    label = "CH3NCH3 + OH <=> CH3NCH2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NCH3 + OH <=> CH3NCH2 + H2O""",
)

entry(
    index = 1301,
    label = "CH3NCH3 + O2 <=> CH3NO + CH3O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+09, 'cm^3/(mol*s)'), n=1, Ea=(6000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NCH3 + O2 <=> CH3NO + CH3O""",
)

entry(
    index = 1302,
    label = "CH3NCH3 + CH3 <=> CH3NCH2 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NCH3 + CH3 <=> CH3NCH2 + CH4""",
)

entry(
    index = 1303,
    label = "CH2NCH2 + H <=> CH3NCH2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.8e+13, 'cm^3/(mol*s)'),
        n = 0.18,
        Ea = (-125, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2NCH2 + H <=> CH3NCH2""",
)

entry(
    index = 1304,
    label = "CH3NCH2 + H <=> CH2NCH2 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.6e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(5464, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NCH2 + H <=> CH2NCH2 + H2""",
)

entry(
    index = 1305,
    label = "CH3NCH2 + H <=> CH3NCH + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(6130, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NCH2 + H <=> CH3NCH + H2""",
)

entry(
    index = 1306,
    label = "CH3NCH2 + O <=> CH2NCH2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(5196, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NCH2 + O <=> CH2NCH2 + OH""",
)

entry(
    index = 1307,
    label = "CH3NCH2 + O <=> CH3NCH + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.2e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(5404, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NCH2 + O <=> CH3NCH + OH""",
)

entry(
    index = 1308,
    label = "CH3NCH2 + OH <=> CH2NCH2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NCH2 + OH <=> CH2NCH2 + H2O""",
)

entry(
    index = 1309,
    label = "CH3NCH2 + OH <=> CH3NCH + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+06, 'cm^3/(mol*s)'), n=2, Ea=(457, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NCH2 + OH <=> CH3NCH + H2O""",
)

entry(
    index = 1310,
    label = "CH3NCH2 + CH3 <=> CH2NCH2 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.5e+06, 'cm^3/(mol*s)'),
        n = 1.87,
        Ea = (9170, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3NCH2 + CH3 <=> CH2NCH2 + CH4""",
)

entry(
    index = 1311,
    label = "CH3NCH2 + CH3 <=> CH3NCH + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(530000, 'cm^3/(mol*s)'), n=1.87, Ea=(9687, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NCH2 + CH3 <=> CH3NCH + CH4""",
)

entry(
    index = 1312,
    label = "CH3NCH2 + NH2 <=> CH2NCH2 + NH3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.8e+06, 'cm^3/(mol*s)'),
        n = 1.94,
        Ea = (5494, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3NCH2 + NH2 <=> CH2NCH2 + NH3""",
)

entry(
    index = 1313,
    label = "CH3NCH2 + NH2 <=> CH3NCH + NH3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.8e+06, 'cm^3/(mol*s)'),
        n = 1.94,
        Ea = (6090, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3NCH2 + NH2 <=> CH3NCH + NH3""",
)

entry(
    index = 1314,
    label = "CH2NCH2 <=> CH3NCH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+45, 's^-1'), n=-10.068, Ea=(66111, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2NCH2 <=> CH3NCH""",
)

entry(
    index = 1315,
    label = "CH2NCH2 + H <=> CH3 + H2CN",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2NCH2 + H <=> CH3 + H2CN""",
)

entry(
    index = 1316,
    label = "CH2NCH2 + O <=> CH2O + H2CN",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2NCH2 + O <=> CH2O + H2CN""",
)

entry(
    index = 1317,
    label = "CH2NCH2 + OH <=> CH2OH + H2CN",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2NCH2 + OH <=> CH2OH + H2CN""",
)

entry(
    index = 1318,
    label = "CH3NCH <=> CH3 + HCN",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.1e+15, 's^-1'), n=-2.375, Ea=(14942, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NCH <=> CH3 + HCN""",
)

entry(
    index = 1319,
    label = "CH3NCH + H <=> CH2NCH2 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NCH + H <=> CH2NCH2 + H""",
)

entry(
    index = 1320,
    label = "CH3NCH + O => CH3 + NCO + H",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(7e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3NCH + O => CH3 + NCO + H""",
)

