#!/usr/bin/env python
# encoding: utf-8

name = "/home/alongd/Code/RMG-Py/importer/Hashemi"
shortDesc = u"/home/alongd/Code/RMG-Py/importer/Hashemi/mech.inp"
longDesc = u"""
Unknown source
"""
entry(
    index = 1,
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
    index = 2,
    label = "H + O2 <=> O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(15286, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H + O2 <=> O + OH""",
)

entry(
    index = 3,
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
    index = 4,
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
    index = 5,
    label = "OH + OH <=> O + H2O",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(
                A = (2.015e+07, 'cm^3/(mol*s)'),
                n = 1.651,
                Ea = (630.7, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.595e+11, 'cm^3/(mol*s)'),
                n = -0.057,
                Ea = (-827, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is OH + OH <=> O + H2O""",
)

entry(
    index = 6,
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
    index = 7,
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
    index = 8,
    label = "H + O <=> OH",
    degeneracy = 1,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(A=(4.7e+18, 'cm^6/(mol^2*s)'), n=-1, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        efficiencies = {'[C-]#[O+]': 1.9, '[H][H]': 2.5, 'O=C=O': 3.8, 'O': 12, '[Ar]': 0.75},
    ),
    shortDesc = u"""The chemkin file reaction is H + O <=> OH""",
)

entry(
    index = 9,
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
    index = 10,
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
    index = 11,
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
                A = (5.11e+07, 'cm^3/(mol*s)'),
                n = 2.182,
                Ea = (11524, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.16e+09, 'cm^3/(mol*s)'),
                n = 1.812,
                Ea = (13163, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2O + H <=> H + CO + H2""",
)

entry(
    index = 39,
    label = "CH2O + H <=> HCO + H2",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.04, 1, 10], 'atm'),
        arrhenius = [
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (7.435e+23, 'cm^3/(mol*s)'),
                        n = -2.732,
                        Ea = (16379.1, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2.098e+10, 'cm^3/(mol*s)'),
                        n = 1.057,
                        Ea = (3719.6, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (1.357e+23, 'cm^3/(mol*s)'),
                        n = -2.355,
                        Ea = (17518.7, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.639e+15, 'cm^3/(mol*s)'),
                        n = -0.444,
                        Ea = (5682.4, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (7.314e+23, 'cm^3/(mol*s)'),
                        n = -2.665,
                        Ea = (17634.6, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (4.189e+09, 'cm^3/(mol*s)'),
                        n = 1.294,
                        Ea = (3590.9, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2O + H <=> HCO + H2""",
)

entry(
    index = 40,
    label = "CH2O + O <=> H + CO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.5141e+21, 'cm^3/(mol*s)'),
        n = -1.903,
        Ea = (22674, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2O + O <=> H + CO + OH""",
)

entry(
    index = 41,
    label = "CH2O + O <=> HCO + OH",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(
                A = (5.634e+31, 'cm^3/(mol*s)'),
                n = -5.189,
                Ea = (19968.2, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.358e+15, 'cm^3/(mol*s)'),
                n = -0.53,
                Ea = (4010.7, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2O + O <=> HCO + OH""",
)

entry(
    index = 42,
    label = "CH2O + O2 <=> HCO + HO2",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(240000, 'cm^3/(mol*s)'), n=2.5, Ea=(36461, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (-1.437e+15, 'cm^3/(mol*s)'),
                n = 0.027,
                Ea = (56388.3, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2O + O2 <=> HCO + HO2""",
)

entry(
    index = 43,
    label = "CH2O + O2 <=> H + CO + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.437e+15, 'cm^3/(mol*s)'),
        n = 0.027,
        Ea = (56388.3, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2O + O2 <=> H + CO + HO2""",
)

entry(
    index = 44,
    label = "CH2O + OH <=> H + CO + H2O",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.04, 1, 10], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (7.03e+10, 'cm^3/(mol*s)'),
                n = 0.911,
                Ea = (8646, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (7.24e+10, 'cm^3/(mol*s)'),
                n = 0.892,
                Ea = (9310, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (8.37e+10, 'cm^3/(mol*s)'),
                n = 0.879,
                Ea = (9843, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2O + OH <=> H + CO + H2O""",
)

entry(
    index = 45,
    label = "CH2O + OH <=> HCO + H2O",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.04, 1, 10], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (3.559e+09, 'cm^3/(mol*s)'),
                n = 1.167,
                Ea = (-205.7, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.854e+09, 'cm^3/(mol*s)'),
                n = 1.256,
                Ea = (-302.4, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.07e+09, 'cm^3/(mol*s)'),
                n = 1.33,
                Ea = (-391.7, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2O + OH <=> HCO + H2O""",
)

entry(
    index = 46,
    label = "CH2O + HO2 <=> HCO + H2O2",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(41000, 'cm^3/(mol*s)'), n=2.5, Ea=(10206, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (-2.454e+14, 'cm^3/(mol*s)'),
                n = 0.027,
                Ea = (30133.3, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2O + HO2 <=> HCO + H2O2""",
)

entry(
    index = 47,
    label = "CH2O + HO2 <=> H + CO + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.4543e+14, 'cm^3/(mol*s)'),
        n = 0.027,
        Ea = (30120, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2O + HO2 <=> H + CO + H2O2""",
)

entry(
    index = 48,
    label = "CH2O + CH3 <=> HCO + CH4",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(32, 'cm^3/(mol*s)'), n=3.36, Ea=(4310, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (-1.916e+11, 'cm^3/(mol*s)'),
                n = 0.887,
                Ea = (24237.3, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2O + CH3 <=> HCO + CH4""",
)

entry(
    index = 49,
    label = "CH2O + CH3 <=> H + CO + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.9155e+11, 'cm^3/(mol*s)'),
        n = 0.887,
        Ea = (24224, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2O + CH3 <=> H + CO + CH4""",
)

entry(
    index = 50,
    label = "HCO <=> H + CO",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(4.93e+16, 's^-1'), n=-0.93, Ea=(19724, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (7.43e+21, 'cm^3/(mol*s)'),
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
    index = 51,
    label = "HCO + H <=> CO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCO + H <=> CO + H2""",
)

entry(
    index = 52,
    label = "HCO + O <=> CO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCO + O <=> CO + OH""",
)

entry(
    index = 53,
    label = "HCO + O <=> CO2 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCO + O <=> CO2 + H""",
)

entry(
    index = 54,
    label = "HCO + OH <=> CO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCO + OH <=> CO + H2O""",
)

entry(
    index = 55,
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
    index = 56,
    label = "HCO + HO2 <=> CO2 + OH + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCO + HO2 <=> CO2 + OH + H""",
)

entry(
    index = 57,
    label = "HCO + HCO <=> CO + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.7e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCO + HCO <=> CO + CH2O""",
)

entry(
    index = 58,
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
    index = 59,
    label = "CH4 + H <=> CH3 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4100, 'cm^3/(mol*s)'), n=3.156, Ea=(8755, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH4 + H <=> CH3 + H2""",
)

entry(
    index = 60,
    label = "CH4 + O <=> CH3 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(440000, 'cm^3/(mol*s)'), n=2.5, Ea=(6577, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH4 + O <=> CH3 + OH""",
)

entry(
    index = 61,
    label = "CH4 + OH <=> CH3 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+06, 'cm^3/(mol*s)'), n=2.182, Ea=(2506, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH4 + OH <=> CH3 + H2O""",
)

entry(
    index = 62,
    label = "CH4 + HO2 <=> CH3 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(47000, 'cm^3/(mol*s)'), n=2.5, Ea=(21000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH4 + HO2 <=> CH3 + H2O2""",
)

entry(
    index = 63,
    label = "CH4 + O2 <=> CH3 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (203000, 'cm^3/(mol*s)'),
        n = 2.745,
        Ea = (51714, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH4 + O2 <=> CH3 + HO2""",
)

entry(
    index = 64,
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
    index = 65,
    label = "CH4 + CH <=> C2H4 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(-400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH4 + CH <=> C2H4 + H""",
)

entry(
    index = 66,
    label = "CH3 <=> CH + H2",
    degeneracy = 1,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(A=(3.1e+15, 'cm^3/(mol*s)'), n=0, Ea=(80871, 'cal/mol'), T0=(1, 'K')),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is CH3 <=> CH + H2""",
)

entry(
    index = 67,
    label = "CH3 <=> CH2 + H",
    degeneracy = 1,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(A=(2.2e+15, 'cm^3/(mol*s)'), n=0, Ea=(82659, 'cal/mol'), T0=(1, 'K')),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is CH3 <=> CH2 + H""",
)

entry(
    index = 68,
    label = "CH2 + H2 <=> CH3 + H",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(440000, 'cm^3/(mol*s)'), n=2.39, Ea=(7350, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(7.2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2 + H2 <=> CH3 + H""",
)

entry(
    index = 69,
    label = "CH3 + O <=> CH2O + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.9e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + O <=> CH2O + H""",
)

entry(
    index = 70,
    label = "CH3 + O <=> H2 + CO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + O <=> H2 + CO + H""",
)

entry(
    index = 71,
    label = "CH3 + OH <=> CH2 + H2O",
    degeneracy = 1,
    duplicate = True,
    kinetics = Arrhenius(A=(43000, 'cm^3/(mol*s)'), n=2.568, Ea=(3997, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + OH <=> CH2 + H2O""",
)

entry(
    index = 72,
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
    index = 73,
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
    index = 74,
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
    index = 75,
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
    index = 76,
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
    index = 77,
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
    index = 78,
    label = "CH3 + HO2 <=> CH3O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0.2688, Ea=(688, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + HO2 <=> CH3O + OH""",
)

entry(
    index = 79,
    label = "CH3 + O2 <=> CH3O + O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(28297, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + O2 <=> CH3O + O""",
)

entry(
    index = 80,
    label = "CH3 + O2 <=> CH2O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.9e+11, 'cm^3/(mol*s)'), n=0, Ea=(9842, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + O2 <=> CH2O + OH""",
)

entry(
    index = 81,
    label = "CH3 + O2 <=> CH3OO",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiPDepArrhenius(
        arrhenius = [
            PDepArrhenius(
                pressures = ([1, 10, 20, 50, 100], 'atm'),
                arrhenius = [
                    Arrhenius(A=(5e+22, 'cm^3/(mol*s)'), n=-3.85, Ea=(2000, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (3.35e+21, 'cm^3/(mol*s)'),
                        n = -3.2,
                        Ea = (2300, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (3.25e+29, 'cm^3/(mol*s)'),
                        n = -5.6,
                        Ea = (6850, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2.83e+18, 'cm^3/(mol*s)'),
                        n = -2.2,
                        Ea = (1400, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.05e+19, 'cm^3/(mol*s)'),
                        n = -2.3,
                        Ea = (1800, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            PDepArrhenius(
                pressures = ([20, 50, 100], 'atm'),
                arrhenius = [
                    Arrhenius(
                        A = (4.1e+20, 'cm^3/(mol*s)'),
                        n = -2.94,
                        Ea = (1900, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (5.6e+28, 'cm^3/(mol*s)'),
                        n = -5.25,
                        Ea = (6850, 'cal/mol'),
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
    index = 82,
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
    index = 83,
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
    index = 84,
    label = "CH3 + CH3 <=> C2H5 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(16055, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + CH3 <=> C2H5 + H""",
)

entry(
    index = 85,
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
    index = 86,
    label = "CH3 + CH <=> C2H3 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + CH <=> C2H3 + H""",
)

entry(
    index = 87,
    label = "CH3 + C <=> C2H2 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + C <=> C2H2 + H""",
)

entry(
    index = 88,
    label = "CH2 <=> CH + H",
    degeneracy = 1,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(A=(5.6e+15, 'cm^3/(mol*s)'), n=0, Ea=(89000, 'cal/mol'), T0=(1, 'K')),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is CH2 <=> CH + H""",
)

entry(
    index = 89,
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
    index = 90,
    label = "CH2 + H <=> CH + H2",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(1.2e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2 + H <=> CH + H2""",
)

entry(
    index = 91,
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
    index = 92,
    label = "CH2 + O <=> CO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+13, 'cm^3/(mol*s)'), n=0, Ea=(536, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2 + O <=> CO + H2""",
)

entry(
    index = 93,
    label = "CH2 + OH <=> CH2O + H",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(
                A = (2.8e+13, 'cm^3/(mol*s)'),
                n = 0.1228,
                Ea = (161, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2 + OH <=> CH2O + H""",
)

entry(
    index = 94,
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
    index = 95,
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
    index = 96,
    label = "CH2 + O2 <=> CH2O + O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.3e+06, 'cm^3/(mol*s)'),
        n = 2.4202,
        Ea = (1604, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2 + O2 <=> CH2O + O""",
)

entry(
    index = 97,
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
    index = 98,
    label = "CH2 + CH2 <=> C2H2 + H + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+13, 'cm^3/(mol*s)'), n=0.0022, Ea=(8, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2 + CH2 <=> C2H2 + H + H""",
)

entry(
    index = 99,
    label = "CH2 + CH2 <=> C2H2 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+13, 'cm^3/(mol*s)'), n=0.0022, Ea=(8, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2 + CH2 <=> C2H2 + H2""",
)

entry(
    index = 100,
    label = "CH2 + CH <=> C2H2 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2 + CH <=> C2H2 + H""",
)

entry(
    index = 101,
    label = "CH2 + C <=> C2H + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2 + C <=> C2H + H""",
)

entry(
    index = 102,
    label = "CH2(S) + N2 <=> CH2 + N2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+13, 'cm^3/(mol*s)'), n=0, Ea=(430, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2(S) + N2 <=> CH2 + N2""",
)

entry(
    index = 103,
    label = "CH2(S) + AR <=> CH2 + AR",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.5e+13, 'cm^3/(mol*s)'), n=0, Ea=(884, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2(S) + AR <=> CH2 + AR""",
)

entry(
    index = 104,
    label = "CH2(S) + H <=> CH2 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2(S) + H <=> CH2 + H""",
)

entry(
    index = 105,
    label = "CH2(S) + O2 <=> CH2 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2(S) + O2 <=> CH2 + O2""",
)

entry(
    index = 106,
    label = "CH2(S) + H2O <=> CH2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2(S) + H2O <=> CH2 + H2O""",
)

entry(
    index = 107,
    label = "CH + H <=> C + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH + H <=> C + H2""",
)

entry(
    index = 108,
    label = "CH + O <=> CO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.7e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH + O <=> CO + H""",
)

entry(
    index = 109,
    label = "CH + OH <=> HCO + H",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (-1.796e+23, 'cm^3/(mol*s)'),
                n = -2.473,
                Ea = (19927.3, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH + OH <=> HCO + H""",
)

entry(
    index = 110,
    label = "CH + OH <=> H + CO + H",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.796e+23, 'cm^3/(mol*s)'),
        n = -2.473,
        Ea = (19927.3, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH + OH <=> H + CO + H""",
)

entry(
    index = 111,
    label = "CH + OH <=> C + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+07, 'cm^3/(mol*s)'), n=2, Ea=(3000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH + OH <=> C + H2O""",
)

entry(
    index = 112,
    label = "CH + O2 <=> HCO + O",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(3.3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (-1.975e+23, 'cm^3/(mol*s)'),
                n = -2.473,
                Ea = (19927.3, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH + O2 <=> HCO + O""",
)

entry(
    index = 113,
    label = "CH + O2 <=> H + CO + O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.975e+23, 'cm^3/(mol*s)'),
        n = -2.473,
        Ea = (19927.3, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH + O2 <=> H + CO + O""",
)

entry(
    index = 114,
    label = "CH + H2O <=> CH2O + H",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8.5e+08, 'cm^3/(mol*s)'),
        n = 1.144,
        Ea = (-2051, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH + H2O <=> CH2O + H""",
)

entry(
    index = 115,
    label = "CH + CO2 <=> HCO + CO",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(
                A = (8.8e+06, 'cm^3/(mol*s)'),
                n = 1.75,
                Ea = (-1040, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (-5.268e+16, 'cm^3/(mol*s)'),
                n = -0.723,
                Ea = (18887.3, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH + CO2 <=> HCO + CO""",
)

entry(
    index = 116,
    label = "CH + CO2 <=> H + CO + CO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.268e+16, 'cm^3/(mol*s)'),
        n = -0.723,
        Ea = (18887.3, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH + CO2 <=> H + CO + CO""",
)

entry(
    index = 117,
    label = "CH + CH2O <=> CH2CO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.5e+13, 'cm^3/(mol*s)'), n=0, Ea=(-517, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH + CH2O <=> CH2CO + H""",
)

entry(
    index = 118,
    label = "C + OH <=> CO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C + OH <=> CO + H""",
)

entry(
    index = 119,
    label = "C + O2 <=> CO + O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C + O2 <=> CO + O""",
)

entry(
    index = 120,
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
    index = 121,
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
    index = 122,
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
    index = 123,
    label = "CH3OH + H <=> CH2OH + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(66000, 'cm^3/(mol*s)'), n=2.728, Ea=(4449, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OH + H <=> CH2OH + H2""",
)

entry(
    index = 124,
    label = "CH3OH + H <=> CH3O + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(41000, 'cm^3/(mol*s)'), n=2.658, Ea=(9221, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OH + H <=> CH3O + H2""",
)

entry(
    index = 125,
    label = "CH3OH + O <=> CH2OH + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.1e+13, 'cm^3/(mol*s)'), n=0, Ea=(5305, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OH + O <=> CH2OH + OH""",
)

entry(
    index = 126,
    label = "CH3OH + O <=> CH3O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.7e+12, 'cm^3/(mol*s)'), n=0, Ea=(5305, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OH + O <=> CH3O + OH""",
)

entry(
    index = 127,
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
    index = 128,
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
    index = 129,
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
    index = 130,
    label = "CH3OH + HO2 <=> CH3O + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (0.0015, 'cm^3/(mol*s)'),
        n = 4.61,
        Ea = (15828, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3OH + HO2 <=> CH3O + H2O2""",
)

entry(
    index = 131,
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
    index = 132,
    label = "CH3O + HO2 <=> CH3OH + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.4e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O + HO2 <=> CH3OH + O2""",
)

entry(
    index = 133,
    label = "CH2OH <=> CH2O + H",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(7.37e+10, 's^-1'), n=0.811, Ea=(39558.7, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.5e+21, 'cm^3/(mol*s)'),
            n = -1.99,
            Ea = (23983.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.844,
        T3 = (900, 'K'),
        T1 = (1, 'K'),
        T2 = (3315, 'K'),
        efficiencies = {'O=C=O': 2, 'O': 6, '[H][H]': 2, '[He]': 0.67, '[C-]#[O+]': 1.5, '[Ar]': 0.85},
    ),
    shortDesc = u"""The chemkin file reaction is CH2OH <=> CH2O + H""",
)

entry(
    index = 134,
    label = "CH2OH + H <=> CH2O + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+06, 'cm^3/(mol*s)'), n=1.86, Ea=(147, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2OH + H <=> CH2O + H2""",
)

entry(
    index = 135,
    label = "CH2OH + O <=> CH2O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(-693, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2OH + O <=> CH2O + OH""",
)

entry(
    index = 136,
    label = "CH2OH + OH <=> CH2O + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2OH + OH <=> CH2O + H2O""",
)

entry(
    index = 137,
    label = "CH2OH + HO2 <=> CH2O + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2OH + HO2 <=> CH2O + H2O2""",
)

entry(
    index = 138,
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
    index = 139,
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
    index = 140,
    label = "CH2OH + HCO <=> CH3OH + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2OH + HCO <=> CH3OH + CO""",
)

entry(
    index = 141,
    label = "CH2OH + HCO <=> CH2O + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2OH + HCO <=> CH2O + CH2O""",
)

entry(
    index = 142,
    label = "CH3OH + HCO <=> CH2OH + CH2O",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(
                A = (4.91e-14, 'cm^3/(mol*s)'),
                n = 7.37,
                Ea = (9730.8, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (-0.0002939, 'cm^3/(mol*s)'),
                n = 4.897,
                Ea = (29658.1, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH3OH + HCO <=> CH2OH + CH2O""",
)

entry(
    index = 143,
    label = "CH3OH + H + CO <=> CH2OH + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8.316e-05, 'cm^6/(mol^2*s)'),
        n = 5.095,
        Ea = (14283.2, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3OH + H + CO <=> CH2OH + CH2O""",
)

entry(
    index = 144,
    label = "CH2OH + CH2OH <=> CH3OH + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2OH + CH2OH <=> CH3OH + CH2O""",
)

entry(
    index = 145,
    label = "CH2OH + CH3O <=> CH3OH + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2OH + CH3O <=> CH3OH + CH2O""",
)

entry(
    index = 146,
    label = "CH3O <=> CH2O + H",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.13e+10, 's^-1'), n=1.21, Ea=(24068.8, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (6.02e+16, 'cm^3/(mol*s)'),
            n = -0.547,
            Ea = (18011.9, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.341,
        T3 = (28, 'K'),
        T1 = (1000, 'K'),
        T2 = (2339, 'K'),
        efficiencies = {'O=C=O': 2, 'O': 6, '[H][H]': 2, '[He]': 0.67, '[O][O]': 1, '[C-]#[O+]': 1.5, '[Ar]': 0.85},
    ),
    shortDesc = u"""The chemkin file reaction is CH3O <=> CH2O + H""",
)

entry(
    index = 147,
    label = "CH3O + H <=> CH2O + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.6e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(-519, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O + H <=> CH2O + H2""",
)

entry(
    index = 148,
    label = "CH3O + O <=> CH2O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O + O <=> CH2O + OH""",
)

entry(
    index = 149,
    label = "CH3O + OH <=> CH2O + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O + OH <=> CH2O + H2O""",
)

entry(
    index = 150,
    label = "CH3O + HO2 <=> CH2O + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O + HO2 <=> CH2O + H2O2""",
)

entry(
    index = 151,
    label = "CH3O + O2 <=> CH2O + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.48, 'cm^3/(mol*s)'), n=3.567, Ea=(-1055, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O + O2 <=> CH2O + HO2""",
)

entry(
    index = 152,
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
    index = 153,
    label = "CH3O + CH3 <=> CH2O + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O + CH3 <=> CH2O + CH4""",
)

entry(
    index = 154,
    label = "CH3OH + CH3 <=> CH3O + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (0.0125, 'cm^3/(mol*s)'),
        n = 4.161,
        Ea = (6002, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3OH + CH3 <=> CH3O + CH4""",
)

entry(
    index = 155,
    label = "CH3OH + CH3 <=> CH2OH + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.19e-07, 'cm^3/(mol*s)'),
        n = 5.58,
        Ea = (3896.3, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3OH + CH3 <=> CH2OH + CH4""",
)

entry(
    index = 156,
    label = "CH3OH + CH2(S) <=> CH2OH + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (0.000354, 'cm^3/(mol*s)'),
        n = 4.19,
        Ea = (3600.4, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3OH + CH2(S) <=> CH2OH + CH3""",
)

entry(
    index = 157,
    label = "CH3OH + CH2(S) <=> CH3O + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (0.000136, 'cm^3/(mol*s)'),
        n = 4.26,
        Ea = (5626, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3OH + CH2(S) <=> CH3O + CH3""",
)

entry(
    index = 158,
    label = "CH3OH + HCO <=> CH3O + CH2O",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(
                A = (0.00121, 'cm^3/(mol*s)'),
                n = 4.34,
                Ea = (23046.1, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (-7.243e+06, 'cm^3/(mol*s)'),
                n = 1.867,
                Ea = (42973.4, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH3OH + HCO <=> CH3O + CH2O""",
)

entry(
    index = 159,
    label = "CH3OH + H + CO <=> CH3O + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.049e+06, 'cm^6/(mol^2*s)'),
        n = 2.065,
        Ea = (27598.5, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3OH + H + CO <=> CH3O + CH2O""",
)

entry(
    index = 160,
    label = "CH3O + CH3O <=> CH3OH + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O + CH3O <=> CH3OH + CH2O""",
)

entry(
    index = 161,
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
    index = 162,
    label = "CH3OOH + H <=> CH2OOH + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.4e+10, 'cm^3/(mol*s)'), n=0, Ea=(1860, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OOH + H <=> CH2OOH + H2""",
)

entry(
    index = 163,
    label = "CH3OOH + H <=> CH3OO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.4e+10, 'cm^3/(mol*s)'), n=0, Ea=(1860, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OOH + H <=> CH3OO + H2""",
)

entry(
    index = 164,
    label = "CH3OOH + H <=> CH3O + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+10, 'cm^3/(mol*s)'), n=0, Ea=(1860, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OOH + H <=> CH3O + H2O""",
)

entry(
    index = 165,
    label = "CH3OOH + O <=> CH2OOH + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(4750, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OOH + O <=> CH2OOH + OH""",
)

entry(
    index = 166,
    label = "CH3OOH + O <=> CH3OO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.7e+12, 'cm^3/(mol*s)'), n=0, Ea=(4750, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OOH + O <=> CH3OO + OH""",
)

entry(
    index = 167,
    label = "CH3OOH + OH <=> CH3OO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.1e+12, 'cm^3/(mol*s)'), n=0, Ea=(-437, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OOH + OH <=> CH3OO + H2O""",
)

entry(
    index = 168,
    label = "CH3OOH + OH <=> CH2OOH + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.2e+11, 'cm^3/(mol*s)'), n=0, Ea=(-258, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OOH + OH <=> CH2OOH + H2O""",
)

entry(
    index = 169,
    label = "CH3OOH + HO2 <=> CH3OO + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(41000, 'cm^3/(mol*s)'), n=2.5, Ea=(10206, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OOH + HO2 <=> CH3OO + H2O2""",
)

entry(
    index = 170,
    label = "CH3OO + H <=> CH3O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OO + H <=> CH3O + OH""",
)

entry(
    index = 171,
    label = "CH3OO + O <=> CH3O + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.85e+10, 'cm^3/(mol*s)'), n=1, Ea=(-724, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OO + O <=> CH3O + O2""",
)

entry(
    index = 172,
    label = "CH3OO + OH <=> CH3OH + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OO + OH <=> CH3OH + O2""",
)

entry(
    index = 173,
    label = "CH3OO + HO2 <=> CH3OOH + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.5e+11, 'cm^3/(mol*s)'), n=0, Ea=(-1490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OO + HO2 <=> CH3OOH + O2""",
)

entry(
    index = 174,
    label = "CH3OO + CH3 <=> CH3O + CH3O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.1e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1411, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OO + CH3 <=> CH3O + CH3O""",
)

entry(
    index = 175,
    label = "CH3OO + CH4 <=> CH3OOH + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (0.00445, 'cm^3/(mol*s)'),
        n = 4.691,
        Ea = (19868, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3OO + CH4 <=> CH3OOH + CH3""",
)

entry(
    index = 176,
    label = "CH3OO + CH2OH <=> CH2O + CH3OOH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OO + CH2OH <=> CH2O + CH3OOH""",
)

entry(
    index = 177,
    label = "CH3OO + HCO <=> CH3O + H + CO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OO + HCO <=> CH3O + H + CO2""",
)

entry(
    index = 178,
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
    index = 179,
    label = "CH3OO + CH2O <=> CH3OOH + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.981e+09, 'cm^3/(mol*s)'),
        n = 1.111,
        Ea = (12499.5, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3OO + CH2O <=> CH3OOH + HCO""",
)

entry(
    index = 180,
    label = "CH3OO + CH2O <=> CH3OOH + H + CO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.454e+14, 'cm^3/(mol*s)'),
        n = 0.027,
        Ea = (30133.3, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3OO + CH2O <=> CH3OOH + H + CO""",
)

entry(
    index = 181,
    label = "CH3OO + CH3O <=> CH2O + CH3OOH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OO + CH3O <=> CH2O + CH3OOH""",
)

entry(
    index = 182,
    label = "CH3OO + CH3OH <=> CH3OOH + CH2OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(19400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OO + CH3OH <=> CH3OOH + CH2OH""",
)

entry(
    index = 183,
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
    index = 184,
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
    index = 185,
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
    index = 186,
    label = "HCOH <=> CO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.56e+13, 's^-1'), n=0, Ea=(49465, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCOH <=> CO + H2""",
)

entry(
    index = 187,
    label = "HCOH <=> CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.15e+13, 's^-1'), n=0, Ea=(32109, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCOH <=> CH2O""",
)

entry(
    index = 188,
    label = "HCOH + H <=> CH2O + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCOH + H <=> CH2O + H""",
)

entry(
    index = 189,
    label = "HCOH + H <=> HCO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCOH + H <=> HCO + H2""",
)

entry(
    index = 190,
    label = "HCOH + H <=> H + CO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCOH + H <=> H + CO + H2""",
)

entry(
    index = 191,
    label = "HCOH + O <=> CO2 + H + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCOH + O <=> CO2 + H + H""",
)

entry(
    index = 192,
    label = "HCOH + OH <=> HCO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCOH + OH <=> HCO + H2O""",
)

entry(
    index = 193,
    label = "HCOH + OH <=> H + CO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCOH + OH <=> H + CO + H2O""",
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
            Arrhenius(A=(7350, 'cm^3/(mol*s)'), n=3.1, Ea=(5340.02, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (3.26e+14, 'cm^3/(mol*s)'),
                n = 0,
                Ea = (13666.8, 'cal/mol'),
                T0 = (1, 'K'),
            ),
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
        A = (1.61e+06, 'cm^3/(mol*s)'),
        n = 2.224,
        Ea = (740.73, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H6 + OH <=> C2H5 + H2O""",
)

entry(
    index = 198,
    label = "C2H6 + HO2 <=> C2H5 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(26, 'cm^3/(mol*s)'), n=3.37, Ea=(15900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H6 + HO2 <=> C2H5 + H2O2""",
)

entry(
    index = 199,
    label = "C2H6 + O2 <=> C2H5 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.92e+07, 'cm^3/(mol*s)'),
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
    label = "C2H6 + CH3OO <=> CH3OOH + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(19, 'cm^3/(mol*s)'), n=3.64, Ea=(17100, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H6 + CH3OO <=> CH3OOH + C2H5""",
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
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.5e+15, 'cm^3/(mol*s)'), n=-1, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(A=(3.3e+31, 'cm^6/(mol^2*s)'), n=-4.9, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        alpha = 1e-30,
        T3 = (540, 'K'),
        T1 = (1e+30, 'K'),
        T2 = (1e-30, 'K'),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is C2H5 + O2 <=> CH3CH2OO""",
)

entry(
    index = 212,
    label = "C2H6 + HCO <=> C2H5 + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.3, 'cm^3/(mol*s)'), n=3.74, Ea=(16933, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H6 + HCO <=> C2H5 + CH2O""",
)

entry(
    index = 213,
    label = "C2H5 + HCO <=> CH2CHO + CH3",
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
    shortDesc = u"""The chemkin file reaction is C2H5 + HCO <=> CH2CHO + CH3""",
)

entry(
    index = 214,
    label = "C2H5 + CH3 <=> C2H4 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5 + CH3 <=> C2H4 + CH4""",
)

entry(
    index = 215,
    label = "C2H5 + HCO <=> C2H6 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5 + HCO <=> C2H6 + CO""",
)

entry(
    index = 216,
    label = "C2H5 + CH3OO <=> CH3O + CH3CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.1e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1410, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5 + CH3OO <=> CH3O + CH3CH2O""",
)

entry(
    index = 217,
    label = "C2H5 + C2H5 <=> C2H6 + C2H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5 + C2H5 <=> C2H6 + C2H4""",
)

entry(
    index = 218,
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
    index = 219,
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
    index = 220,
    label = "C2H4 + H <=> C2H3 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(240, 'cm^3/(mol*s)'), n=3.62, Ea=(11266, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4 + H <=> C2H3 + H2""",
)

entry(
    index = 221,
    label = "C2H4 + O <=> CH3 + HCO",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(3.9e+12, 'cm^3/(mol*s)'), n=0, Ea=(1494, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (-2.335e+22, 'cm^3/(mol*s)'),
                n = -2.473,
                Ea = (21421.3, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(A=(6.2e+13, 'cm^3/(mol*s)'), n=0, Ea=(6855, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (-3.711e+23, 'cm^3/(mol*s)'),
                n = -2.473,
                Ea = (26782.3, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H4 + O <=> CH3 + HCO""",
)

entry(
    index = 222,
    label = "C2H4 + O <=> CH3 + H + CO",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(
                A = (2.335e+22, 'cm^3/(mol*s)'),
                n = -2.473,
                Ea = (21421.3, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (3.711e+23, 'cm^3/(mol*s)'),
                n = -2.473,
                Ea = (26782.3, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H4 + O <=> CH3 + H + CO""",
)

entry(
    index = 223,
    label = "C2H4 + O <=> CH2CHO + H",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(1.7e+12, 'cm^3/(mol*s)'), n=0, Ea=(1494, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(6855, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H4 + O <=> CH2CHO + H""",
)

entry(
    index = 224,
    label = "C2H4 + OH <=> C2H3 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.13, 'cm^3/(mol*s)'), n=4.2, Ea=(-860, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4 + OH <=> C2H3 + H2O""",
)

entry(
    index = 225,
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
    index = 226,
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
    index = 227,
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
    index = 228,
    label = "C2H4 + OH <=> CH2CH2OH",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiPDepArrhenius(
        arrhenius = [
            PDepArrhenius(
                pressures = ([0.01, 0.025, 0.1, 1, 10, 100], 'atm'),
                arrhenius = [
                    Arrhenius(
                        A = (2.76e+47, 'cm^3/(mol*s)'),
                        n = -11.64,
                        Ea = (11099, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (6.02e+37, 'cm^3/(mol*s)'),
                        n = -9.76,
                        Ea = (1995, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (6.02e+37, 'cm^3/(mol*s)'),
                        n = -9.65,
                        Ea = (2363, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (6.02e+37, 'cm^3/(mol*s)'),
                        n = -8.14,
                        Ea = (8043, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (6.02e+37, 'cm^3/(mol*s)'),
                        n = -7.77,
                        Ea = (10736, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (6.02e+37, 'cm^3/(mol*s)'),
                        n = -7.44,
                        Ea = (14269, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            PDepArrhenius(
                pressures = ([0.025, 0.1, 1, 10, 100], 'atm'),
                arrhenius = [
                    Arrhenius(
                        A = (4.96e+37, 'cm^3/(mol*s)'),
                        n = -8.68,
                        Ea = (5355, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2.56e+35, 'cm^3/(mol*s)'),
                        n = -7.79,
                        Ea = (5017, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (7.29e+31, 'cm^3/(mol*s)'),
                        n = -6.91,
                        Ea = (2855, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (3.02e+26, 'cm^3/(mol*s)'),
                        n = -4.87,
                        Ea = (2297, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2.79e+19, 'cm^3/(mol*s)'),
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
    index = 229,
    label = "C2H4 + HO2 <=> C2H5 + O2",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01316, 0.98692, 9.86923, 98.6923], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (2.8e+12, 'cm^3/(mol*s)'),
                n = -0.447,
                Ea = (12658.5, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (3.407e+08, 'cm^3/(mol*s)'),
                n = 0.818,
                Ea = (12089.9, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (9.504e+14, 'cm^3/(mol*s)'),
                n = -1.012,
                Ea = (16582, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.553e+22, 'cm^3/(mol*s)'),
                n = -3.074,
                Ea = (22635.2, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H4 + HO2 <=> C2H5 + O2""",
)

entry(
    index = 230,
    label = "C2H4 + HO2 <=> cC2H4O + OH",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01316, 0.98692, 9.86923, 98.6923], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (33420, 'cm^3/(mol*s)'),
                n = 2.311,
                Ea = (13735.5, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.796e+06, 'cm^3/(mol*s)'),
                n = 1.809,
                Ea = (14760.6, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.181e+09, 'cm^3/(mol*s)'),
                n = 1.006,
                Ea = (16673.5, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.892e+15, 'cm^3/(mol*s)'),
                n = -0.736,
                Ea = (21411.3, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H4 + HO2 <=> cC2H4O + OH""",
)

entry(
    index = 231,
    label = "C2H4 + HO2 <=> CH3CH2OO",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01316, 0.98692, 9.86923, 98.6923, 1000], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (2.603e+47, 'cm^3/(mol*s)'),
                n = -12.45,
                Ea = (18714.7, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (6.596e+47, 'cm^3/(mol*s)'),
                n = -12.159,
                Ea = (20887.1, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.602e+45, 'cm^3/(mol*s)'),
                n = -11.101,
                Ea = (21253.5, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.005e+41, 'cm^3/(mol*s)'),
                n = -9.543,
                Ea = (21399.2, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (11.24, 'cm^3/(mol*s)'),
                n = 2.95,
                Ea = (7093.14, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H4 + HO2 <=> CH3CH2OO""",
)

entry(
    index = 232,
    label = "C2H4 + HO2 <=> CH2CH2OOH",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01316, 0.98692, 9.86923, 98.6923, 1000], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (1.323e+07, 'cm^3/(mol*s)'),
                n = -0.142,
                Ea = (11075.9, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (143800, 'cm^3/(mol*s)'),
                n = 0.878,
                Ea = (7815.54, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (3.377e+16, 'cm^3/(mol*s)'),
                n = -2.209,
                Ea = (13289.2, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (6.496e+26, 'cm^3/(mol*s)'),
                n = -4.893,
                Ea = (19947.9, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (28180, 'cm^3/(mol*s)'),
                n = 2.487,
                Ea = (14734, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H4 + HO2 <=> CH2CH2OOH""",
)

entry(
    index = 233,
    label = "C2H4 + O2 <=> C2H3 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.1e+13, 'cm^3/(mol*s)'), n=0, Ea=(60010, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4 + O2 <=> C2H3 + HO2""",
)

entry(
    index = 234,
    label = "C2H4 + CH3 <=> C2H3 + CH4",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(976, 'cm^3/(mol*s)'), n=2.947, Ea=(15148, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (8.13e-05, 'cm^3/(mol*s)'),
                n = 4.417,
                Ea = (8835.8, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H4 + CH3 <=> C2H3 + CH4""",
)

entry(
    index = 235,
    label = "C2H4 + CH2(S) <=> C2H3 + CH3",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiPDepArrhenius(
        arrhenius = [
            PDepArrhenius(
                pressures = ([0.01, 0.1, 1, 10, 100], 'atm'),
                arrhenius = [
                    Arrhenius(
                        A = (1.77e+19, 'cm^3/(mol*s)'),
                        n = -1.95,
                        Ea = (6787, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.68e+19, 'cm^3/(mol*s)'),
                        n = -1.8,
                        Ea = (4310, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (4.16e+24, 'cm^3/(mol*s)'),
                        n = -3.19,
                        Ea = (9759, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (7.89e+24, 'cm^3/(mol*s)'),
                        n = -3.08,
                        Ea = (13894, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (7.36e+29, 'cm^3/(mol*s)'),
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
                    Arrhenius(A=(2.26e+11, 'cm^3/(mol*s)'), n=0.54, Ea=(48, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (4.92e+09, 'cm^3/(mol*s)'),
                        n = 1.02,
                        Ea = (600, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.47e+08, 'cm^3/(mol*s)'),
                        n = 1.33,
                        Ea = (1228, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (8.11e+10, 'cm^3/(mol*s)'),
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
    index = 236,
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
    index = 237,
    label = "C2H3 + H <=> C2H2 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3 + H <=> C2H2 + H2""",
)

entry(
    index = 238,
    label = "C2H3 + O <=> CH2CO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3 + O <=> CH2CO + H""",
)

entry(
    index = 239,
    label = "C2H3 + OH <=> C2H2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3 + OH <=> C2H2 + H2O""",
)

entry(
    index = 240,
    label = "C2H3 + HO2 <=> CH2CHO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3 + HO2 <=> CH2CHO + OH""",
)

entry(
    index = 241,
    label = "C2H3 + O2 <=> CH2CHOO",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
        arrhenius = [
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (1.6e+24, 'cm^3/(mol*s)'),
                        n = -5.45,
                        Ea = (9662, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.8e-09, 'cm^3/(mol*s)'),
                        n = 4.15,
                        Ea = (-4707, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (3.5e+56, 'cm^3/(mol*s)'),
                        n = -15.01,
                        Ea = (19160, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2.4e+22, 'cm^3/(mol*s)'),
                        n = -4.52,
                        Ea = (2839, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (1.3e+64, 'cm^3/(mol*s)'),
                        n = -16.97,
                        Ea = (21290, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(2e+26, 'cm^3/(mol*s)'), n=-5.43, Ea=(2725, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (3.3e+61, 'cm^3/(mol*s)'),
                        n = -15.79,
                        Ea = (20150, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (6.1e+28, 'cm^3/(mol*s)'),
                        n = -5.89,
                        Ea = (3154, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (7.3e+53, 'cm^3/(mol*s)'),
                        n = -13.11,
                        Ea = (17300, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2.1e+29, 'cm^3/(mol*s)'),
                        n = -5.8,
                        Ea = (3520, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (4.2e+48, 'cm^3/(mol*s)'),
                        n = -11.21,
                        Ea = (16000, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (3.5e+28, 'cm^3/(mol*s)'),
                        n = -5.37,
                        Ea = (3636, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (2.3e+43, 'cm^3/(mol*s)'),
                        n = -9.38,
                        Ea = (14810, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (3.3e+27, 'cm^3/(mol*s)'),
                        n = -4.95,
                        Ea = (3610, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (3.4e+39, 'cm^3/(mol*s)'),
                        n = -8.04,
                        Ea = (14360, 'cal/mol'),
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
    index = 242,
    label = "C2H3 + O2 <=> CHCHO + OH",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
        arrhenius = [
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (3.9e+11, 'cm^3/(mol*s)'),
                        n = -0.11,
                        Ea = (2131, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(9.9e+11, 'cm^3/(mol*s)'), n=-0.66, Ea=(-1, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(1.1e+09, 'cm^3/(mol*s)'), n=0.55, Ea=(46, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (6.9e+14, 'cm^3/(mol*s)'),
                        n = -1.16,
                        Ea = (4542, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(8.5e+08, 'cm^3/(mol*s)'), n=0.56, Ea=(1, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (2.8e+13, 'cm^3/(mol*s)'),
                        n = -0.72,
                        Ea = (3479, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(2.8e+14, 'cm^3/(mol*s)'), n=-1.83, Ea=(5, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(5e+11, 'cm^3/(mol*s)'), n=-0.14, Ea=(1995, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (2.6e+20, 'cm^3/(mol*s)'),
                        n = -2.84,
                        Ea = (7530, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2.4e+10, 'cm^3/(mol*s)'),
                        n = 0.23,
                        Ea = (1573, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(9.2e+14, 'cm^3/(mol*s)'), n=-2.26, Ea=(0, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (1.7e+14, 'cm^3/(mol*s)'),
                        n = -0.82,
                        Ea = (4450, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (6.1e+25, 'cm^3/(mol*s)'),
                        n = -4.21,
                        Ea = (13050, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.4e+11, 'cm^3/(mol*s)'),
                        n = 0.05,
                        Ea = (3774, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (1.7e+30, 'cm^3/(mol*s)'),
                        n = -5.35,
                        Ea = (18430, 'cal/mol'),
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
    index = 243,
    label = "C2H3 + O2 <=> CH2CO + OH",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
        arrhenius = [
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(870, 'cm^3/(mol*s)'), n=2.41, Ea=(6061, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(0.18, 'cm^3/(mol*s)'), n=3.12, Ea=(1331, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(890, 'cm^3/(mol*s)'), n=2.41, Ea=(6078, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(0.21, 'cm^3/(mol*s)'), n=3.11, Ea=(1383, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(940, 'cm^3/(mol*s)'), n=2.4, Ea=(6112, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(0.27, 'cm^3/(mol*s)'), n=3.08, Ea=(1496, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(1100, 'cm^3/(mol*s)'), n=2.39, Ea=(6180, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(0.53, 'cm^3/(mol*s)'), n=3.01, Ea=(1777, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(1100, 'cm^3/(mol*s)'), n=2.38, Ea=(6179, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.4, 'cm^3/(mol*s)'), n=2.9, Ea=(2225, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(1400, 'cm^3/(mol*s)'), n=2.36, Ea=(6074, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(0.42, 'cm^3/(mol*s)'), n=2.93, Ea=(2052, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (2.5e+06, 'cm^3/(mol*s)'),
                        n = 1.42,
                        Ea = (8480, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (0.00012, 'cm^3/(mol*s)'),
                        n = 4.21,
                        Ea = (2043, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (1.7e+10, 'cm^3/(mol*s)'),
                        n = 0.36,
                        Ea = (12010, 'cal/mol'),
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
    index = 244,
    label = "C2H3 + O2 <=> CH2CHO + O",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
        arrhenius = [
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (7.2e+20, 'cm^3/(mol*s)'),
                        n = -2.67,
                        Ea = (6742, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.2e+10, 'cm^3/(mol*s)'),
                        n = 0.62,
                        Ea = (-278, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(7e+20, 'cm^3/(mol*s)'), n=-2.67, Ea=(6713, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (1.3e+10, 'cm^3/(mol*s)'),
                        n = 0.62,
                        Ea = (-248, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(9e+20, 'cm^3/(mol*s)'), n=-2.7, Ea=(6724, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.5e+10, 'cm^3/(mol*s)'), n=0.6, Ea=(-163, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (6.5e+20, 'cm^3/(mol*s)'),
                        n = -2.65,
                        Ea = (6489, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(1.8e+10, 'cm^3/(mol*s)'), n=0.58, Ea=(38, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (4.1e+20, 'cm^3/(mol*s)'),
                        n = -2.53,
                        Ea = (6406, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(8.9e+09, 'cm^3/(mol*s)'), n=0.67, Ea=(248, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (1.6e+23, 'cm^3/(mol*s)'),
                        n = -3.22,
                        Ea = (8697, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(6.7e+09, 'cm^3/(mol*s)'), n=0.72, Ea=(778, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (2.9e+25, 'cm^3/(mol*s)'),
                        n = -3.77,
                        Ea = (11530, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.4e+09, 'cm^3/(mol*s)'),
                        n = 0.92,
                        Ea = (1219, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (9.3e+25, 'cm^3/(mol*s)'),
                        n = -3.8,
                        Ea = (13910, 'cal/mol'),
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
    index = 245,
    label = "C2H3 + O2 <=> C2H2 + HO2",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
        arrhenius = [
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (1.1e+07, 'cm^3/(mol*s)'),
                        n = 1.28,
                        Ea = (3322, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(48, 'cm^3/(mol*s)'), n=2.75, Ea=(-796, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (7.8e+06, 'cm^3/(mol*s)'),
                        n = 1.33,
                        Ea = (3216, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(52, 'cm^3/(mol*s)'), n=2.73, Ea=(-768, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (1.2e+07, 'cm^3/(mol*s)'),
                        n = 1.27,
                        Ea = (3311, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(56, 'cm^3/(mol*s)'), n=2.73, Ea=(-659, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (2.2e+07, 'cm^3/(mol*s)'),
                        n = 1.19,
                        Ea = (3367, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(46, 'cm^3/(mol*s)'), n=2.76, Ea=(-493, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(1.1e+08, 'cm^3/(mol*s)'), n=1, Ea=(3695, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(3.8, 'cm^3/(mol*s)'), n=3.07, Ea=(-601, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (1.3e+11, 'cm^3/(mol*s)'),
                        n = 0.12,
                        Ea = (5872, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(5.5, 'cm^3/(mol*s)'), n=3.07, Ea=(86, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (1.2e+09, 'cm^3/(mol*s)'),
                        n = 0.82,
                        Ea = (5617, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(4.5e+08, 'cm^3/(mol*s)'), n=0, Ea=(955, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (1.1e+17, 'cm^3/(mol*s)'),
                        n = -1.45,
                        Ea = (12230, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(20, 'cm^3/(mol*s)'), n=2.94, Ea=(1847, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H3 + O2 <=> C2H2 + HO2""",
)

entry(
    index = 246,
    label = "C2H3 + O2 <=> OCHCHO + H",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
        arrhenius = [
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (4.8e+14, 'cm^3/(mol*s)'),
                        n = -1.03,
                        Ea = (912, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (0.00028, 'cm^3/(mol*s)'),
                        n = 4.04,
                        Ea = (-7019, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(5e+14, 'cm^3/(mol*s)'), n=-1.04, Ea=(923, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (0.00035, 'cm^3/(mol*s)'),
                        n = 4.01,
                        Ea = (-6978, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (6.4e+14, 'cm^3/(mol*s)'),
                        n = -1.07,
                        Ea = (983, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (0.00097, 'cm^3/(mol*s)'),
                        n = 3.89,
                        Ea = (-6768, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (3.7e+15, 'cm^3/(mol*s)'),
                        n = -1.29,
                        Ea = (1441, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(0.5, 'cm^3/(mol*s)'), n=3.15, Ea=(-5496, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (2.4e+18, 'cm^3/(mol*s)'),
                        n = -2.13,
                        Ea = (3234, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (130000, 'cm^3/(mol*s)'),
                        n = 1.67,
                        Ea = (-2931, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (1.3e+15, 'cm^3/(mol*s)'),
                        n = -1.09,
                        Ea = (2393, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (4.5e+15, 'cm^3/(mol*s)'),
                        n = -3.08,
                        Ea = (-4836, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (3.6e+33, 'cm^3/(mol*s)'),
                        n = -6.5,
                        Ea = (14910, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(3.8e+10, 'cm^3/(mol*s)'), n=0.22, Ea=(941, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (3.3e+31, 'cm^3/(mol*s)'),
                        n = -5.76,
                        Ea = (16250, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(2.8e+08, 'cm^3/(mol*s)'), n=0.83, Ea=(858, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H3 + O2 <=> OCHCHO + H""",
)

entry(
    index = 247,
    label = "C2H3 + O2 <=> CH2O + HCO",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (2.8e+36, 'cm^3/(mol*s)'),
                n = -7.6,
                Ea = (12610, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (5.1e+15, 'cm^3/(mol*s)'),
                n = -1.28,
                Ea = (513, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (2.2e+36, 'cm^3/(mol*s)'),
                        n = -7.57,
                        Ea = (12490, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (5.3e+15, 'cm^3/(mol*s)'),
                        n = -1.29,
                        Ea = (521, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (3e+35, 'cm^3/(mol*s)'),
                        n = -7.32,
                        Ea = (11820, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (6.8e+15, 'cm^3/(mol*s)'),
                        n = -1.31,
                        Ea = (646, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (1.6e+36, 'cm^3/(mol*s)'),
                        n = -7.47,
                        Ea = (12460, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.1e+16, 'cm^3/(mol*s)'),
                        n = -1.36,
                        Ea = (1066, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (5.8e+35, 'cm^3/(mol*s)'),
                        n = -7.2,
                        Ea = (13430, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2.8e+15, 'cm^3/(mol*s)'),
                        n = -1.18,
                        Ea = (1429, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (3.5e+20, 'cm^3/(mol*s)'),
                        n = -2.57,
                        Ea = (5578, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.1e+69, 'cm^3/(mol*s)'),
                        n = -19.23,
                        Ea = (14760, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (3e+33, 'cm^3/(mol*s)'),
                        n = -6.28,
                        Ea = (16000, 'cal/mol'),
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
    index = 248,
    label = "C2H3 + O2 <=> CH2O + H + CO",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
        arrhenius = [
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (6.5e+36, 'cm^3/(mol*s)'),
                        n = -7.6,
                        Ea = (12640, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.2e+16, 'cm^3/(mol*s)'),
                        n = -1.28,
                        Ea = (515, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (6.3e+36, 'cm^3/(mol*s)'),
                        n = -7.6,
                        Ea = (12610, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.2e+16, 'cm^3/(mol*s)'),
                        n = -1.28,
                        Ea = (513, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (5.1e+36, 'cm^3/(mol*s)'),
                        n = -7.57,
                        Ea = (12490, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.3e+16, 'cm^3/(mol*s)'),
                        n = -1.29,
                        Ea = (521, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (7.1e+35, 'cm^3/(mol*s)'),
                        n = -7.32,
                        Ea = (11820, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.6e+16, 'cm^3/(mol*s)'),
                        n = -1.31,
                        Ea = (646, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (3.7e+36, 'cm^3/(mol*s)'),
                        n = -7.47,
                        Ea = (12460, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2.4e+16, 'cm^3/(mol*s)'),
                        n = -1.36,
                        Ea = (1066, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (1.3e+36, 'cm^3/(mol*s)'),
                        n = -7.2,
                        Ea = (13430, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (6.6e+15, 'cm^3/(mol*s)'),
                        n = -1.18,
                        Ea = (1429, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (8.3e+20, 'cm^3/(mol*s)'),
                        n = -2.57,
                        Ea = (5578, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2.7e+69, 'cm^3/(mol*s)'),
                        n = -19.23,
                        Ea = (14760, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (7.1e+33, 'cm^3/(mol*s)'),
                        n = -6.28,
                        Ea = (16000, 'cal/mol'),
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
    index = 249,
    label = "C2H3 + O2 <=> CH3O + CO",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
        arrhenius = [
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (8.2e+18, 'cm^3/(mol*s)'),
                        n = -2.66,
                        Ea = (3201, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.3e+09, 'cm^3/(mol*s)'),
                        n = 0.18,
                        Ea = (-1717, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (4.1e+14, 'cm^3/(mol*s)'),
                        n = -1.32,
                        Ea = (886, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (6e+11, 'cm^3/(mol*s)'),
                        n = -2.93,
                        Ea = (-9564, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (4.3e+14, 'cm^3/(mol*s)'),
                        n = -1.33,
                        Ea = (901, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2.9e+11, 'cm^3/(mol*s)'),
                        n = -2.93,
                        Ea = (-10120, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=-0.33, Ea=(-748, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (5.8e+21, 'cm^3/(mol*s)'),
                        n = -3.54,
                        Ea = (4772, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(1.9e+12, 'cm^3/(mol*s)'), n=-3, Ea=(-8995, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(5e+15, 'cm^3/(mol*s)'), n=-1.62, Ea=(1849, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(1.9e+24, 'cm^3/(mol*s)'), n=-5.63, Ea=(2, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (9.3e+16, 'cm^3/(mol*s)'),
                        n = -1.96,
                        Ea = (3324, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (1.1e+18, 'cm^3/(mol*s)'),
                        n = -2.22,
                        Ea = (5178, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1e+72, 'cm^3/(mol*s)'),
                        n = -20.69,
                        Ea = (15860, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (5.8e+32, 'cm^3/(mol*s)'),
                        n = -6.45,
                        Ea = (16810, 'cal/mol'),
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
    shortDesc = u"""The chemkin file reaction is C2H3 + O2 <=> CH3O + CO""",
)

entry(
    index = 250,
    label = "C2H3 + O2 <=> CO2 + CH3",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
        arrhenius = [
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (2.4e+35, 'cm^3/(mol*s)'),
                        n = -7.76,
                        Ea = (12630, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (6.3e+13, 'cm^3/(mol*s)'),
                        n = -1.16,
                        Ea = (406, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (1.7e+35, 'cm^3/(mol*s)'),
                        n = -7.72,
                        Ea = (12520, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (6.2e+13, 'cm^3/(mol*s)'),
                        n = -1.16,
                        Ea = (401, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (4.5e+34, 'cm^3/(mol*s)'),
                        n = -7.55,
                        Ea = (12140, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (6.1e+13, 'cm^3/(mol*s)'),
                        n = -1.16,
                        Ea = (397, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (7.3e+31, 'cm^3/(mol*s)'),
                        n = -6.7,
                        Ea = (10440, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (5.3e+13, 'cm^3/(mol*s)'),
                        n = -1.14,
                        Ea = (447, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (3.6e+35, 'cm^3/(mol*s)'),
                        n = -7.75,
                        Ea = (12830, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.5e+14, 'cm^3/(mol*s)'),
                        n = -1.26,
                        Ea = (988, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (2.1e+35, 'cm^3/(mol*s)'),
                        n = -7.53,
                        Ea = (14050, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=-1.11, Ea=(1409, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (3.8e+18, 'cm^3/(mol*s)'),
                        n = -2.44,
                        Ea = (5408, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.4e+70, 'cm^3/(mol*s)'),
                        n = -20.11,
                        Ea = (15430, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (1.2e+32, 'cm^3/(mol*s)'),
                        n = -6.32,
                        Ea = (16190, 'cal/mol'),
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
    index = 251,
    label = "C2H3 + CH2O <=> C2H4 + HCO",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.001, 0.01, 0.1, 1, 10, 100, 1000], 'atm'),
        arrhenius = [
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (1.11e+07, 'cm^3/(mol*s)'),
                        n = 1.09,
                        Ea = (1807.2, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (-2.333e+16, 'cm^3/(mol*s)'),
                        n = -1.269,
                        Ea = (20616.8, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (2.47e+07, 'cm^3/(mol*s)'),
                        n = 0.993,
                        Ea = (1994.9, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (-5.192e+16, 'cm^3/(mol*s)'),
                        n = -1.366,
                        Ea = (20804.5, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (2.47e+08, 'cm^3/(mol*s)'),
                        n = 0.704,
                        Ea = (2596.2, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (-1.479e+18, 'cm^3/(mol*s)'),
                        n = -1.769,
                        Ea = (22523.5, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (1.42e+10, 'cm^3/(mol*s)'),
                        n = 0.209,
                        Ea = (3934.2, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (-8.5e+19, 'cm^3/(mol*s)'),
                        n = -2.264,
                        Ea = (23861.5, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (3.45e+13, 'cm^3/(mol*s)'),
                        n = -0.726,
                        Ea = (6944.3, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (-4.364e+23, 'cm^3/(mol*s)'),
                        n = -3.278,
                        Ea = (27795.2, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (3.31e+14, 'cm^3/(mol*s)'),
                        n = -0.866,
                        Ea = (10965.7, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (-4.187e+24, 'cm^3/(mol*s)'),
                        n = -3.418,
                        Ea = (31816.6, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(16.5, 'cm^3/(mol*s)'), n=3.17, Ea=(9399.8, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (-2.087e+11, 'cm^3/(mol*s)'),
                        n = 0.618,
                        Ea = (30250.7, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H3 + CH2O <=> C2H4 + HCO""",
)

entry(
    index = 252,
    label = "C2H3 + CH2O <=> C2H4 + H + CO",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.001, 0.01, 0.1, 1, 10, 100, 1000], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (2.333e+16, 'cm^3/(mol*s)'),
                n = -1.269,
                Ea = (20616.8, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (5.192e+16, 'cm^3/(mol*s)'),
                n = -1.366,
                Ea = (20804.5, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.479e+18, 'cm^3/(mol*s)'),
                n = -1.769,
                Ea = (22523.5, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (8.5e+19, 'cm^3/(mol*s)'),
                n = -2.264,
                Ea = (23861.5, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (4.364e+23, 'cm^3/(mol*s)'),
                n = -3.278,
                Ea = (27795.2, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (4.187e+24, 'cm^3/(mol*s)'),
                n = -3.418,
                Ea = (31816.6, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.087e+11, 'cm^3/(mol*s)'),
                n = 0.618,
                Ea = (30250.7, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H3 + CH2O <=> C2H4 + H + CO""",
)

entry(
    index = 253,
    label = "C2H3 + HCO <=> C2H4 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3 + HCO <=> C2H4 + CO""",
)

entry(
    index = 254,
    label = "C2H3 + CH3 <=> C2H2 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+12, 'cm^3/(mol*s)'), n=0, Ea=(-765, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3 + CH3 <=> C2H2 + CH4""",
)

entry(
    index = 255,
    label = "C2H3 + CH <=> CH2 + C2H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3 + CH <=> CH2 + C2H2""",
)

entry(
    index = 256,
    label = "C2H3 + C2H3 <=> C2H4 + C2H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3 + C2H3 <=> C2H4 + C2H2""",
)

entry(
    index = 257,
    label = "C2H3 + C2H <=> C2H2 + C2H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3 + C2H <=> C2H2 + C2H2""",
)

entry(
    index = 258,
    label = "C2H3 + CH3OH <=> C2H4 + CH3O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.03e+06, 'cm^3/(mol*s)'),
        n = 1.51,
        Ea = (26630, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H3 + CH3OH <=> C2H4 + CH3O""",
)

entry(
    index = 259,
    label = "C2H3 + CH3OH <=> C2H4 + CH2OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (0.0175, 'cm^3/(mol*s)'),
        n = 4.02,
        Ea = (23370, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H3 + CH3OH <=> C2H4 + CH2OH""",
)

entry(
    index = 260,
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
    index = 261,
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
    index = 262,
    label = "C2H + H2 <=> C2H2 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(410000, 'cm^3/(mol*s)'), n=2.39, Ea=(864, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H + H2 <=> C2H2 + H""",
)

entry(
    index = 263,
    label = "C2H2 + O <=> HCCO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.4e+07, 'cm^3/(mol*s)'), n=2, Ea=(1900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H2 + O <=> HCCO + H""",
)

entry(
    index = 264,
    label = "C2H2 + O <=> CH2 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.1e+06, 'cm^3/(mol*s)'), n=2, Ea=(1900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H2 + O <=> CH2 + CO""",
)

entry(
    index = 265,
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
    index = 266,
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
    index = 267,
    label = "C2H2 + OH <=> CHCHOH",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.025, 0.1, 1, 10, 100], 'atm'),
        arrhenius = [
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (2.9e+64, 'cm^3/(mol*s)'),
                        n = -18.57,
                        Ea = (10009, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2.6e+33, 'cm^3/(mol*s)'),
                        n = -7.36,
                        Ea = (6392, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (4.7e+59, 'cm^3/(mol*s)'),
                        n = -16.87,
                        Ea = (9087, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (4.4e+32, 'cm^3/(mol*s)'),
                        n = -7.02,
                        Ea = (5933, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (1.2e+28, 'cm^3/(mol*s)'),
                        n = -5.56,
                        Ea = (3724, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (6.4e+42, 'cm^3/(mol*s)'),
                        n = -9.96,
                        Ea = (11737, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (1.9e+44, 'cm^3/(mol*s)'),
                        n = -11.38,
                        Ea = (6299, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (3.5e+31, 'cm^3/(mol*s)'),
                        n = -6.2,
                        Ea = (6635, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (1.5e+24, 'cm^3/(mol*s)'),
                        n = -4.06,
                        Ea = (3261, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (4.5e+31, 'cm^3/(mol*s)'),
                        n = -5.92,
                        Ea = (8761, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (6.2e+20, 'cm^3/(mol*s)'),
                        n = -2.8,
                        Ea = (2831, 'cal/mol'),
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
    index = 268,
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
    index = 269,
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
    index = 270,
    label = "C2H2 + HO2 <=> CH2CHOO",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
        arrhenius = [
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (4.99e+06, 'cm^3/(mol*s)'),
                        n = -1.02,
                        Ea = (9152, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.88e+26, 'cm^3/(mol*s)'),
                        n = -8.34,
                        Ea = (9249, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (6.02e+17, 'cm^3/(mol*s)'),
                        n = -3.82,
                        Ea = (10790, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (5.26e+129, 'cm^3/(mol*s)'),
                        n = -41.74,
                        Ea = (35930, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (2.47e+48, 'cm^3/(mol*s)'),
                        n = -12.82,
                        Ea = (25220, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.96e+18, 'cm^3/(mol*s)'),
                        n = -3.67,
                        Ea = (10480, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (4.06e+50, 'cm^3/(mol*s)'),
                        n = -13.07,
                        Ea = (27220, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (4.93e+21, 'cm^3/(mol*s)'),
                        n = -4.37,
                        Ea = (12220, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (9.08e+46, 'cm^3/(mol*s)'),
                        n = -11.57,
                        Ea = (26880, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.92e+22, 'cm^3/(mol*s)'),
                        n = -4.28,
                        Ea = (13080, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (4.6e+43, 'cm^3/(mol*s)'),
                        n = -10.24,
                        Ea = (26930, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2.11e+21, 'cm^3/(mol*s)'),
                        n = -3.78,
                        Ea = (13380, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (5.61e+38, 'cm^3/(mol*s)'),
                        n = -8.49,
                        Ea = (26210, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.39e+20, 'cm^3/(mol*s)'),
                        n = -3.3,
                        Ea = (13410, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (2.53e+35, 'cm^3/(mol*s)'),
                        n = -7.26,
                        Ea = (26390, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.42e+19, 'cm^3/(mol*s)'),
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
    index = 271,
    label = "C2H2 + HO2 <=> CHCHO + OH",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
        arrhenius = [
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (5.49e+09, 'cm^3/(mol*s)'),
                        n = 0.91,
                        Ea = (18500, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2.41e+07, 'cm^3/(mol*s)'),
                        n = 1.54,
                        Ea = (14690, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (5.93e+09, 'cm^3/(mol*s)'),
                        n = 0.9,
                        Ea = (18550, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2.48e+07, 'cm^3/(mol*s)'),
                        n = 1.54,
                        Ea = (14700, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (6.8e+09, 'cm^3/(mol*s)'),
                        n = 0.88,
                        Ea = (18640, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2.63e+07, 'cm^3/(mol*s)'),
                        n = 1.54,
                        Ea = (14730, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (1.56e+10, 'cm^3/(mol*s)'),
                        n = 0.77,
                        Ea = (19040, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2.5e+07, 'cm^3/(mol*s)'),
                        n = 1.56,
                        Ea = (14790, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (3.48e+09, 'cm^3/(mol*s)'),
                        n = 0.99,
                        Ea = (18810, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.47e+08, 'cm^3/(mol*s)'),
                        n = 1.32,
                        Ea = (15090, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (5.39e+10, 'cm^3/(mol*s)'),
                        n = 0.61,
                        Ea = (20740, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.61e+08, 'cm^3/(mol*s)'),
                        n = 1.36,
                        Ea = (15420, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (3.7e+08, 'cm^3/(mol*s)'),
                        n = 1.23,
                        Ea = (15960, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.67e+07, 'cm^3/(mol*s)'),
                        n = 1.59,
                        Ea = (15910, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (1.45e+11, 'cm^3/(mol*s)'),
                        n = 0.48,
                        Ea = (17730, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (7.21e+06, 'cm^3/(mol*s)'),
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
    index = 272,
    label = "C2H2 + HO2 <=> CH2CO + OH",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
        arrhenius = [
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (6.25e-07, 'cm^3/(mol*s)'),
                        n = 4.75,
                        Ea = (15530, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.31e-14, 'cm^3/(mol*s)'),
                        n = 6.58,
                        Ea = (10270, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (6.7e-07, 'cm^3/(mol*s)'),
                        n = 4.74,
                        Ea = (15550, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.29e-14, 'cm^3/(mol*s)'),
                        n = 6.59,
                        Ea = (10330, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (4.18e-07, 'cm^3/(mol*s)'),
                        n = 4.81,
                        Ea = (15410, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (3.99e-14, 'cm^3/(mol*s)'),
                        n = 6.36,
                        Ea = (10270, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (5.28e-07, 'cm^3/(mol*s)'),
                        n = 4.78,
                        Ea = (15460, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (3.28e-15, 'cm^3/(mol*s)'),
                        n = 6.7,
                        Ea = (10090, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (1.04e-06, 'cm^3/(mol*s)'),
                        n = 4.69,
                        Ea = (15640, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (8.71e-21, 'cm^3/(mol*s)'),
                        n = 8.3,
                        Ea = (8107, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (4.68e-05, 'cm^3/(mol*s)'),
                        n = 4.22,
                        Ea = (16780, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (8.36e-22, 'cm^3/(mol*s)'),
                        n = 8.76,
                        Ea = (8804, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(0.899, 'cm^3/(mol*s)'), n=2.97, Ea=(19730, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (6.87e-14, 'cm^3/(mol*s)'),
                        n = 6.67,
                        Ea = (13130, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(3580, 'cm^3/(mol*s)'), n=1.97, Ea=(23010, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (6.63e-12, 'cm^3/(mol*s)'),
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
    index = 273,
    label = "C2H2 + HO2 <=> CH2CHO + O",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
        arrhenius = [
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (5.5e+06, 'cm^3/(mol*s)'),
                        n = 1.19,
                        Ea = (12880, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (0.000294, 'cm^3/(mol*s)'),
                        n = 4.16,
                        Ea = (7736, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (1.16e+08, 'cm^3/(mol*s)'),
                        n = 0.77,
                        Ea = (13600, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (0.00614, 'cm^3/(mol*s)'),
                        n = 3.81,
                        Ea = (8394, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (1.2e+07, 'cm^3/(mol*s)'),
                        n = 1.09,
                        Ea = (13050, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (0.000544, 'cm^3/(mol*s)'),
                        n = 4.09,
                        Ea = (8044, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (3.02e+07, 'cm^3/(mol*s)'),
                        n = 0.98,
                        Ea = (13310, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (0.000248, 'cm^3/(mol*s)'),
                        n = 4.19,
                        Ea = (8203, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (1.98e+74, 'cm^3/(mol*s)'),
                        n = -16.33,
                        Ea = (109200, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(65700, 'cm^3/(mol*s)'), n=1.85, Ea=(12360, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (7.5e+14, 'cm^3/(mol*s)'),
                        n = -1.17,
                        Ea = (18350, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(0.292, 'cm^3/(mol*s)'), n=3.38, Ea=(10590, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (8.63e+18, 'cm^3/(mol*s)'),
                        n = -2.27,
                        Ea = (22230, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(1.95, 'cm^3/(mol*s)'), n=3.17, Ea=(11740, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (5.78e+18, 'cm^3/(mol*s)'),
                        n = -2.09,
                        Ea = (24350, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(0.11, 'cm^3/(mol*s)'), n=3.52, Ea=(11980, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H2 + HO2 <=> CH2CHO + O""",
)

entry(
    index = 274,
    label = "C2H2 + HO2 <=> OCHCHO + H",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
        arrhenius = [
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (8.51e+07, 'cm^3/(mol*s)'),
                        n = 0.48,
                        Ea = (11720, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2.43e-06, 'cm^3/(mol*s)'),
                        n = 4.43,
                        Ea = (5578, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (7.43e+07, 'cm^3/(mol*s)'),
                        n = 0.5,
                        Ea = (11690, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(2e-06, 'cm^3/(mol*s)'), n=4.45, Ea=(5564, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (7.91e+07, 'cm^3/(mol*s)'),
                        n = 0.49,
                        Ea = (11700, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.81e-06, 'cm^3/(mol*s)'),
                        n = 4.46,
                        Ea = (5654, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (2.18e+09, 'cm^3/(mol*s)'),
                        n = 0.06,
                        Ea = (12470, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2.24e-05, 'cm^3/(mol*s)'),
                        n = 4.17,
                        Ea = (6416, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (7e+49, 'cm^3/(mol*s)'),
                        n = -10.18,
                        Ea = (77110, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (765000, 'cm^3/(mol*s)'),
                        n = 1.18,
                        Ea = (11340, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (4.06e+16, 'cm^3/(mol*s)'),
                        n = -2.03,
                        Ea = (17630, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(0.0201, 'cm^3/(mol*s)'), n=3.38, Ea=(8696, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (9.38e+16, 'cm^3/(mol*s)'),
                        n = -2.03,
                        Ea = (19590, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (0.00606, 'cm^3/(mol*s)'),
                        n = 3.53,
                        Ea = (9217, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (5.91e+21, 'cm^3/(mol*s)'),
                        n = -3.32,
                        Ea = (25030, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (0.0676, 'cm^3/(mol*s)'),
                        n = 3.27,
                        Ea = (10760, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H2 + HO2 <=> OCHCHO + H""",
)

entry(
    index = 275,
    label = "C2H2 + HO2 <=> CH2O + HCO",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
        arrhenius = [
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (3.9e+13, 'cm^3/(mol*s)'),
                        n = -1.17,
                        Ea = (13750, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(8.43, 'cm^3/(mol*s)'), n=2.56, Ea=(7382, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(4.26, 'cm^3/(mol*s)'), n=2.64, Ea=(7253, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (1.56e+13, 'cm^3/(mol*s)'),
                        n = -1.05,
                        Ea = (13520, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (2.59e-06, 'cm^3/(mol*s)'),
                        n = 4.34,
                        Ea = (4525, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(6.9e+09, 'cm^3/(mol*s)'), n=0, Ea=(11720, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (3.33e+102, 'cm^3/(mol*s)'),
                        n = -24.18,
                        Ea = (138600, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (8.07e+07, 'cm^3/(mol*s)'),
                        n = 0.6,
                        Ea = (10850, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (5.22e+15, 'cm^3/(mol*s)'),
                        n = -1.75,
                        Ea = (15180, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(3.54, 'cm^3/(mol*s)'), n=2.69, Ea=(8025, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (7.32e+35, 'cm^3/(mol*s)'),
                        n = -7.77,
                        Ea = (26970, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (9.84e+06, 'cm^3/(mol*s)'),
                        n = 0.91,
                        Ea = (11710, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (1.78e+28, 'cm^3/(mol*s)'),
                        n = -5.3,
                        Ea = (25130, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(17900, 'cm^3/(mol*s)'), n=1.7, Ea=(11250, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (2.47e+16, 'cm^3/(mol*s)'),
                        n = -1.7,
                        Ea = (20030, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (4.32e-06, 'cm^3/(mol*s)'),
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
    index = 276,
    label = "C2H2 + HO2 <=> CH2O + H + CO",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
        arrhenius = [
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (9.1e+13, 'cm^3/(mol*s)'),
                        n = -1.17,
                        Ea = (13750, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(19.7, 'cm^3/(mol*s)'), n=2.56, Ea=(7382, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(9.94, 'cm^3/(mol*s)'), n=2.64, Ea=(7253, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (3.63e+13, 'cm^3/(mol*s)'),
                        n = -1.05,
                        Ea = (13520, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (6.05e-06, 'cm^3/(mol*s)'),
                        n = 4.34,
                        Ea = (4525, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(1.61e+10, 'cm^3/(mol*s)'), n=0, Ea=(11720, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (7.77e+102, 'cm^3/(mol*s)'),
                        n = -24.18,
                        Ea = (138600, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.88e+08, 'cm^3/(mol*s)'),
                        n = 0.6,
                        Ea = (10850, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (1.22e+16, 'cm^3/(mol*s)'),
                        n = -1.75,
                        Ea = (15180, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(8.26, 'cm^3/(mol*s)'), n=2.69, Ea=(8025, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (1.71e+36, 'cm^3/(mol*s)'),
                        n = -7.77,
                        Ea = (26970, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2.3e+07, 'cm^3/(mol*s)'),
                        n = 0.91,
                        Ea = (11710, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (4.14e+28, 'cm^3/(mol*s)'),
                        n = -5.3,
                        Ea = (25130, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(41900, 'cm^3/(mol*s)'), n=1.7, Ea=(11250, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (5.77e+16, 'cm^3/(mol*s)'),
                        n = -1.7,
                        Ea = (20030, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.01e-05, 'cm^3/(mol*s)'),
                        n = 4.31,
                        Ea = (6829, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H2 + HO2 <=> CH2O + H + CO""",
)

entry(
    index = 277,
    label = "C2H2 + HO2 <=> CO + CH3O",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
        arrhenius = [
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(3.54e+11, 'cm^3/(mol*s)'), n=0, Ea=(49510, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(28900, 'cm^3/(mol*s)'), n=1.23, Ea=(9903, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (2.78e+08, 'cm^3/(mol*s)'),
                        n = 0.01,
                        Ea = (11920, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (9.67e-07, 'cm^3/(mol*s)'),
                        n = 4.15,
                        Ea = (5173, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (8.06e+07, 'cm^3/(mol*s)'),
                        n = 0.18,
                        Ea = (11650, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.84e-08, 'cm^3/(mol*s)'),
                        n = 4.62,
                        Ea = (4517, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (8.94e+69, 'cm^3/(mol*s)'),
                        n = -15.85,
                        Ea = (102500, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (538000, 'cm^3/(mol*s)'),
                        n = 0.86,
                        Ea = (10700, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (5.66e+12, 'cm^3/(mol*s)'),
                        n = -1.25,
                        Ea = (14570, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (0.000537, 'cm^3/(mol*s)'),
                        n = 3.42,
                        Ea = (7218, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (3.3e+23, 'cm^3/(mol*s)'),
                        n = -4.45,
                        Ea = (21210, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(286, 'cm^3/(mol*s)'), n=1.84, Ea=(10460, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (2.43e+22, 'cm^3/(mol*s)'),
                        n = -3.96,
                        Ea = (22650, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(8.11, 'cm^3/(mol*s)'), n=2.3, Ea=(10560, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (1.17e+18, 'cm^3/(mol*s)'),
                        n = -2.57,
                        Ea = (22360, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (0.000686, 'cm^3/(mol*s)'),
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
    index = 278,
    label = "C2H2 + HO2 <=> CO2 + CH3",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
        arrhenius = [
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (1.15e-07, 'cm^3/(mol*s)'),
                        n = 4.31,
                        Ea = (4614, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(2.01e+08, 'cm^3/(mol*s)'), n=0, Ea=(11790, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (1.1e-07, 'cm^3/(mol*s)'),
                        n = 4.32,
                        Ea = (4622, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(2.01e+08, 'cm^3/(mol*s)'), n=0, Ea=(11780, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (1.75e+142, 'cm^3/(mol*s)'),
                        n = -35.04,
                        Ea = (188700, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (155000, 'cm^3/(mol*s)'),
                        n = 0.95,
                        Ea = (10200, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (3.96e+84, 'cm^3/(mol*s)'),
                        n = -19.8,
                        Ea = (119800, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.38e+06, 'cm^3/(mol*s)'),
                        n = 0.68,
                        Ea = (10810, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (5.02e+13, 'cm^3/(mol*s)'),
                        n = -1.6,
                        Ea = (14980, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(0.00929, 'cm^3/(mol*s)'), n=3, Ea=(7659, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (8.56e+28, 'cm^3/(mol*s)'),
                        n = -6.15,
                        Ea = (24030, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(18600, 'cm^3/(mol*s)'), n=1.26, Ea=(11230, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (1.28e+27, 'cm^3/(mol*s)'),
                        n = -5.42,
                        Ea = (25380, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(289, 'cm^3/(mol*s)'), n=1.79, Ea=(11240, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (1.71e+15, 'cm^3/(mol*s)'),
                        n = -1.8,
                        Ea = (20370, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
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
    index = 279,
    label = "C2H2 + O2 <=> HCO + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.008e+18, 'cm^3/(mol*s)'),
        n = -1.601,
        Ea = (55209.9, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H2 + O2 <=> HCO + HCO""",
)

entry(
    index = 280,
    label = "C2H2 + O2 <=> HCO + H + CO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.651e+25, 'cm^3/(mol*s)'),
        n = -3.259,
        Ea = (74127.2, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H2 + O2 <=> HCO + H + CO""",
)

entry(
    index = 281,
    label = "C2H2 + O2 <=> H + CO + H + CO",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(
                A = (6.094e+26, 'cm^3/(mol*s)'),
                n = -3.276,
                Ea = (110815, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.185e+32, 'cm^3/(mol*s)'),
                n = -4.946,
                Ea = (93104.6, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H2 + O2 <=> H + CO + H + CO""",
)

entry(
    index = 282,
    label = "C2H2 + CH2(S) <=> C2H2 + CH2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H2 + CH2(S) <=> C2H2 + CH2""",
)

entry(
    index = 283,
    label = "H2CC + H <=> C2H2 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2CC + H <=> C2H2 + H""",
)

entry(
    index = 284,
    label = "H2CC + OH <=> CH2CO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2CC + OH <=> CH2CO + H""",
)

entry(
    index = 285,
    label = "H2CC + O2 <=> CH2 + CO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2CC + O2 <=> CH2 + CO2""",
)

entry(
    index = 286,
    label = "C2 + H2 <=> C2H + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(400000, 'cm^3/(mol*s)'), n=2.4, Ea=(1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2 + H2 <=> C2H + H""",
)

entry(
    index = 287,
    label = "C2H + O <=> CH + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H + O <=> CH + CO""",
)

entry(
    index = 288,
    label = "C2H + OH <=> HCCO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H + OH <=> HCCO + H""",
)

entry(
    index = 289,
    label = "C2H + OH <=> C2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+07, 'cm^3/(mol*s)'), n=2, Ea=(8000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H + OH <=> C2 + H2O""",
)

entry(
    index = 290,
    label = "C2H + O2 <=> CO + CO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.7e+13, 'cm^3/(mol*s)'), n=-0.16, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H + O2 <=> CO + CO + H""",
)

entry(
    index = 291,
    label = "C2H + CH4 <=> CH3 + C2H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.2e+12, 'cm^3/(mol*s)'), n=0, Ea=(976, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H + CH4 <=> CH3 + C2H2""",
)

entry(
    index = 292,
    label = "C2 <=> C + C",
    degeneracy = 1,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(A=(1.5e+16, 'cm^3/(mol*s)'), n=0, Ea=(142300, 'cal/mol'), T0=(1, 'K')),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is C2 <=> C + C""",
)

entry(
    index = 293,
    label = "C2 + O <=> C + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2 + O <=> C + CO""",
)

entry(
    index = 294,
    label = "C2 + OH <=> C2O + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2 + OH <=> C2O + H""",
)

entry(
    index = 295,
    label = "C2 + O2 <=> CO + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+12, 'cm^3/(mol*s)'), n=0, Ea=(980, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2 + O2 <=> CO + CO""",
)

entry(
    index = 296,
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
    index = 297,
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
    index = 298,
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
    index = 299,
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
    index = 300,
    label = "CH3CH2O + H <=> CH3CH2OH",
    degeneracy = 1,
    kinetics = Lindemann(
        arrheniusHigh = Arrhenius(A=(3.1e+11, 'cm^3/(mol*s)'), n=0.894, Ea=(13, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.8e+51, 'cm^6/(mol^2*s)'),
            n = -15.55,
            Ea = (11101, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is CH3CH2O + H <=> CH3CH2OH""",
)

entry(
    index = 301,
    label = "CH3CH2OH + H <=> CH3CHOH + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8800, 'cm^3/(mol*s)'), n=2.68, Ea=(2913, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OH + H <=> CH3CHOH + H2""",
)

entry(
    index = 302,
    label = "CH3CH2OH + H <=> CH2CH2OH + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5300, 'cm^3/(mol*s)'), n=2.81, Ea=(7491, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OH + H <=> CH2CH2OH + H2""",
)

entry(
    index = 303,
    label = "CH3CH2OH + H <=> CH3CH2O + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(950, 'cm^3/(mol*s)'), n=3.14, Ea=(8696, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OH + H <=> CH3CH2O + H2""",
)

entry(
    index = 304,
    label = "CH3CH2OH + O <=> CH2CH2OH + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(970, 'cm^3/(mol*s)'), n=3.23, Ea=(4660, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OH + O <=> CH2CH2OH + OH""",
)

entry(
    index = 305,
    label = "CH3CH2OH + O <=> CH3CHOH + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(150000, 'cm^3/(mol*s)'), n=2.47, Ea=(876, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OH + O <=> CH3CHOH + OH""",
)

entry(
    index = 306,
    label = "CH3CH2OH + O <=> CH3CH2O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0015, 'cm^3/(mol*s)'), n=4.7, Ea=(1730, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OH + O <=> CH3CH2O + OH""",
)

entry(
    index = 307,
    label = "CH3CH2OH + OH <=> CH3CHOH + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(450, 'cm^3/(mol*s)'), n=3.11, Ea=(-2666, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OH + OH <=> CH3CHOH + H2O""",
)

entry(
    index = 308,
    label = "CH3CH2OH + OH <=> CH2CH2OH + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9400, 'cm^3/(mol*s)'), n=2.67, Ea=(-1004, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OH + OH <=> CH2CH2OH + H2O""",
)

entry(
    index = 309,
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
    index = 310,
    label = "CH3CH2OH + HO2 <=> CH3CHOH + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8200, 'cm^3/(mol*s)'), n=2.55, Ea=(10750, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OH + HO2 <=> CH3CHOH + H2O2""",
)

entry(
    index = 311,
    label = "CH3CH2OH + HO2 <=> CH2CH2OH + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(12000, 'cm^3/(mol*s)'), n=2.55, Ea=(15750, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OH + HO2 <=> CH2CH2OH + H2O2""",
)

entry(
    index = 312,
    label = "CH3CH2OH + HO2 <=> CH3CH2O + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(24000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OH + HO2 <=> CH3CH2O + H2O2""",
)

entry(
    index = 313,
    label = "CH3CH2OH + CH3 <=> CH3CHOH + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(20, 'cm^3/(mol*s)'), n=3.37, Ea=(7630, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OH + CH3 <=> CH3CHOH + CH4""",
)

entry(
    index = 314,
    label = "CH3CH2OH + CH3 <=> CH2CH2OH + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2, 'cm^3/(mol*s)'), n=3.57, Ea=(7717, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OH + CH3 <=> CH2CH2OH + CH4""",
)

entry(
    index = 315,
    label = "CH3CH2OH + CH3 <=> CH3CH2O + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(330, 'cm^3/(mol*s)'), n=3.3, Ea=(12283, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OH + CH3 <=> CH3CH2O + CH4""",
)

entry(
    index = 316,
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
    index = 317,
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
    index = 318,
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
    index = 319,
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
    index = 320,
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
    index = 321,
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
    index = 322,
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
    index = 323,
    label = "CH3CHOH + O <=> CH3CHO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHOH + O <=> CH3CHO + OH""",
)

entry(
    index = 324,
    label = "CH3CHOH + OH <=> CH3CHO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHOH + OH <=> CH3CHO + H2O""",
)

entry(
    index = 325,
    label = "CH3CHOH + HO2 <=> CH3CHO + OH + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHOH + HO2 <=> CH3CHO + OH + OH""",
)

entry(
    index = 326,
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
    index = 327,
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
    index = 328,
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
    index = 329,
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
    index = 330,
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
    index = 331,
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
    index = 332,
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
    index = 333,
    label = "CH2CH2OH + O <=> CH2O + CH2OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CH2OH + O <=> CH2O + CH2OH""",
)

entry(
    index = 334,
    label = "CH2CH2OH + OH <=> CH2CHOH + H2O",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.001, 0.01, 0.1, 1, 10, 100], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (1.2926e+19, 'cm^3/(mol*s)'),
                n = -1.96,
                Ea = (272.7, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.2184e+19, 'cm^3/(mol*s)'),
                n = -1.9533,
                Ea = (238.8, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (4.1052e+19, 'cm^3/(mol*s)'),
                n = -2.1007,
                Ea = (625.4, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (7.9406e+22, 'cm^3/(mol*s)'),
                n = -2.9892,
                Ea = (3862.6, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.7926e+24, 'cm^3/(mol*s)'),
                n = -3.3287,
                Ea = (7748.8, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (4.6961e+18, 'cm^3/(mol*s)'),
                n = -1.5805,
                Ea = (7999.2, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CH2OH + OH <=> CH2CHOH + H2O""",
)

entry(
    index = 335,
    label = "CH2CH2OH + HO2 <=> CH3CH2OH + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CH2OH + HO2 <=> CH3CH2OH + O2""",
)

entry(
    index = 336,
    label = "CH2CH2OH + HO2 => CH2OH + CH2O + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CH2OH + HO2 => CH2OH + CH2O + OH""",
)

entry(
    index = 337,
    label = "CH2CH2OH + O2 <=> CH2CHOH + HO2",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.013, 0.1, 1, 10, 100], 'atm'),
        arrhenius = [
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (1.3e+53, 'cm^3/(mol*s)'),
                        n = -11.88,
                        Ea = (35927, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (2.3e+10, 'cm^3/(mol*s)'),
                        n = -0.15,
                        Ea = (-791, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=-0.79, Ea=(877, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (2.8e+61, 'cm^3/(mol*s)'),
                        n = -14.17,
                        Ea = (43492, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (3.6e+13, 'cm^3/(mol*s)'),
                        n = -0.88,
                        Ea = (3074, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(6000, 'cm^3/(mol*s)'), n=-10, Ea=(199, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (4.4e+20, 'cm^3/(mol*s)'),
                        n = -2.85,
                        Ea = (8516, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(6000, 'cm^3/(mol*s)'), n=-10, Ea=(199, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (1.9e+30, 'cm^3/(mol*s)'),
                        n = -5.51,
                        Ea = (16616, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(A=(6000, 'cm^3/(mol*s)'), n=-10, Ea=(199, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CH2OH + O2 <=> CH2CHOH + HO2""",
)

entry(
    index = 338,
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
    index = 339,
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
    index = 340,
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
    index = 341,
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
    index = 342,
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
    index = 343,
    label = "CH3CH2O + H <=> CH3CHO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.5e+09, 'cm^3/(mol*s)'), n=1.15, Ea=(673, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2O + H <=> CH3CHO + H2""",
)

entry(
    index = 344,
    label = "CH3CH2O + OH <=> CH3CHO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2O + OH <=> CH3CHO + H2O""",
)

entry(
    index = 345,
    label = "CH3CH2O + O2 <=> CH3CHO + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.5e+10, 'cm^3/(mol*s)'), n=0, Ea=(645, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2O + O2 <=> CH3CHO + HO2""",
)

entry(
    index = 346,
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
    index = 347,
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
    index = 348,
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
    index = 349,
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
    index = 350,
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
    index = 351,
    label = "CH3CHO + H <=> CH3CO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(130000, 'cm^3/(mol*s)'), n=2.58, Ea=(1219, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHO + H <=> CH3CO + H2""",
)

entry(
    index = 352,
    label = "CH3CHO + H <=> CH2CHO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2700, 'cm^3/(mol*s)'), n=3.1, Ea=(5203, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHO + H <=> CH2CHO + H2""",
)

entry(
    index = 353,
    label = "CH3CHO + O <=> CH3CO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(1808, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHO + O <=> CH3CO + OH""",
)

entry(
    index = 354,
    label = "CH3CHO + O <=> CH2CHO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(90000, 'cm^3/(mol*s)'), n=2.8, Ea=(5800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHO + O <=> CH2CHO + OH""",
)

entry(
    index = 355,
    label = "CH3CHO + OH <=> CH3CO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(-709, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHO + OH <=> CH3CO + H2O""",
)

entry(
    index = 356,
    label = "CH3CHO + OH <=> CH2CHO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.5e+13, 'cm^3/(mol*s)'), n=0, Ea=(5313, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHO + OH <=> CH2CHO + H2O""",
)

entry(
    index = 357,
    label = "CH3CHO + HO2 <=> CH3CO + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+13, 'cm^3/(mol*s)'), n=0, Ea=(16293, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHO + HO2 <=> CH3CO + H2O2""",
)

entry(
    index = 358,
    label = "CH3CHO + HO2 <=> CH2CHO + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.1e+13, 'cm^3/(mol*s)'), n=0, Ea=(23248, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHO + HO2 <=> CH2CHO + H2O2""",
)

entry(
    index = 359,
    label = "CH3CHO + O2 <=> CH3CO + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(120000, 'cm^3/(mol*s)'), n=2.5, Ea=(37554, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHO + O2 <=> CH3CO + HO2""",
)

entry(
    index = 360,
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
    index = 361,
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
    index = 362,
    label = "CH3CHO + CH3 <=> CH2CHO + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.18, 'cm^3/(mol*s)'), n=3.44, Ea=(10384, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHO + CH3 <=> CH2CHO + CH4""",
)

entry(
    index = 363,
    label = "CH3CHO + CH3O <=> CH3CO + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(2981, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHO + CH3O <=> CH3CO + CH3OH""",
)

entry(
    index = 364,
    label = "CH3CHO + CH3O <=> CH2CHO + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+11, 'cm^3/(mol*s)'), n=0, Ea=(7100, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHO + CH3O <=> CH2CHO + CH3OH""",
)

entry(
    index = 365,
    label = "CH3CHO + CH3OO <=> CH3CO + CH3OOH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+13, 'cm^3/(mol*s)'), n=0, Ea=(16293, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHO + CH3OO <=> CH3CO + CH3OOH""",
)

entry(
    index = 366,
    label = "CH3CHO + CH3OO <=> CH2CHO + CH3OOH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.1e+13, 'cm^3/(mol*s)'), n=0, Ea=(23248, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHO + CH3OO <=> CH2CHO + CH3OOH""",
)

entry(
    index = 367,
    label = "cC2H4O <=> CH2CHO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+13, 's^-1'), n=0.2, Ea=(71780, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC2H4O <=> CH2CHO + H""",
)

entry(
    index = 368,
    label = "cC2H4O <=> CH3 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.6e+13, 's^-1'), n=0.4, Ea=(61880, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC2H4O <=> CH3 + HCO""",
)

entry(
    index = 369,
    label = "cC2H4O <=> CH3CO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 's^-1'), n=0.25, Ea=(65310, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC2H4O <=> CH3CO + H""",
)

entry(
    index = 370,
    label = "cC2H4O <=> CH2CO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.6e+12, 's^-1'), n=-0.2, Ea=(63030, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC2H4O <=> CH2CO + H2""",
)

entry(
    index = 371,
    label = "cC2H4O <=> CH3CHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.2e+12, 's^-1'), n=-0.75, Ea=(46424, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC2H4O <=> CH3CHO""",
)

entry(
    index = 372,
    label = "cC2H4O <=> C2H2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.6e+12, 's^-1'), n=0.06, Ea=(69530, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC2H4O <=> C2H2 + H2O""",
)

entry(
    index = 373,
    label = "cC2H4O + H <=> CH3CHO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(10950, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC2H4O + H <=> CH3CHO + H""",
)

entry(
    index = 374,
    label = "cC2H4O + H <=> cC2H3O + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(8310, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC2H4O + H <=> cC2H3O + H2""",
)

entry(
    index = 375,
    label = "cC2H4O + H <=> C2H3 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+09, 'cm^3/(mol*s)'), n=0, Ea=(5000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC2H4O + H <=> C2H3 + H2O""",
)

entry(
    index = 376,
    label = "cC2H4O + H <=> C2H4 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.5e+10, 'cm^3/(mol*s)'), n=0, Ea=(5000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC2H4O + H <=> C2H4 + OH""",
)

entry(
    index = 377,
    label = "cC2H4O + O <=> cC2H3O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.9e+12, 'cm^3/(mol*s)'), n=0, Ea=(5250, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC2H4O + O <=> cC2H3O + OH""",
)

entry(
    index = 378,
    label = "cC2H4O + OH <=> cC2H3O + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(3610, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC2H4O + OH <=> cC2H3O + H2O""",
)

entry(
    index = 379,
    label = "cC2H4O + HO2 <=> cC2H3O + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+12, 'cm^3/(mol*s)'), n=0, Ea=(17000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC2H4O + HO2 <=> cC2H3O + H2O2""",
)

entry(
    index = 380,
    label = "cC2H4O + O2 <=> cC2H3O + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(61500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC2H4O + O2 <=> cC2H3O + HO2""",
)

entry(
    index = 381,
    label = "cC2H4O + CH3 <=> cC2H3O + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.1e+12, 'cm^3/(mol*s)'), n=0, Ea=(11830, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC2H4O + CH3 <=> cC2H3O + CH4""",
)

entry(
    index = 382,
    label = "CH2CHOH + H <=> CH2CHO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1500, 'cm^3/(mol*s)'), n=3.077, Ea=(7230, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHOH + H <=> CH2CHO + H2""",
)

entry(
    index = 383,
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
    index = 384,
    label = "CH2CHOH + O <=> CH2OH + HCO",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(3.9e+12, 'cm^3/(mol*s)'), n=0, Ea=(1494, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (-2.335e+22, 'cm^3/(mol*s)'),
                n = -2.473,
                Ea = (21421.3, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(A=(6.2e+13, 'cm^3/(mol*s)'), n=0, Ea=(6855, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (-3.711e+23, 'cm^3/(mol*s)'),
                n = -2.473,
                Ea = (26782.3, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHOH + O <=> CH2OH + HCO""",
)

entry(
    index = 385,
    label = "CH2CHOH + O <=> CH2OH + H + CO",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(
                A = (3.711e+23, 'cm^3/(mol*s)'),
                n = -2.473,
                Ea = (26782.3, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.335e+22, 'cm^3/(mol*s)'),
                n = -2.473,
                Ea = (21421.3, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHOH + O <=> CH2OH + H + CO""",
)

entry(
    index = 386,
    label = "CH2CHOH + O <=> CH2CHO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.6e+07, 'cm^3/(mol*s)'), n=2, Ea=(4400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHOH + O <=> CH2CHO + OH""",
)

entry(
    index = 387,
    label = "CH2CHOH + OH <=> CHCHOH + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.13, 'cm^3/(mol*s)'), n=4.2, Ea=(-860, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHOH + OH <=> CHCHOH + H2O""",
)

entry(
    index = 388,
    label = "CH2CHOH + OH <=> CH2CHO + H2O",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(
                A = (2372, 'cm^3/(mol*s)'),
                n = 2.82,
                Ea = (-691.48, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (7.905e+07, 'cm^3/(mol*s)'),
                n = 1.18,
                Ea = (-303.17, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHOH + OH <=> CH2CHO + H2O""",
)

entry(
    index = 389,
    label = "CH2CHOH + HO2 <=> CH2CHO + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.6e+12, 'cm^3/(mol*s)'), n=0, Ea=(16293, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHOH + HO2 <=> CH2CHO + H2O2""",
)

entry(
    index = 390,
    label = "CH2CHOH + HO2 <=> CH3CHO + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(150000, 'cm^3/(mol*s)'), n=1.67, Ea=(6810, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHOH + HO2 <=> CH3CHO + HO2""",
)

entry(
    index = 391,
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
                A = (-2.095e+17, 'cm^3/(mol*s)'),
                n = -0.673,
                Ea = (58927.3, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHOH + O2 => CH2O + HCO + OH""",
)

entry(
    index = 392,
    label = "CH2CHOH + O2 => CH2O + H + CO + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (2.095e+17, 'cm^3/(mol*s)'),
        n = -0.673,
        Ea = (58927.3, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHOH + O2 => CH2O + H + CO + OH""",
)

entry(
    index = 393,
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
    index = 394,
    label = "CHCHOH + H <=> CH2CHO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CHCHOH + H <=> CH2CHO + H""",
)

entry(
    index = 395,
    label = "CHCHOH + H <=> HCCOH + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CHCHOH + H <=> HCCOH + H2""",
)

entry(
    index = 396,
    label = "CHCHOH + O <=> OCHCHO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CHCHOH + O <=> OCHCHO + H""",
)

entry(
    index = 397,
    label = "CHCHOH + OH <=> HCCOH + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CHCHOH + OH <=> HCCOH + H2O""",
)

entry(
    index = 398,
    label = "CHCHOH + O2 <=> OCHCHO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(-187, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CHCHOH + O2 <=> OCHCHO + OH""",
)

entry(
    index = 399,
    label = "CHCHOH + O2 <=> HOCHO + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.848e+14, 'cm^3/(mol*s)'),
        n = -0.586,
        Ea = (1237.3, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CHCHOH + O2 <=> HOCHO + HCO""",
)

entry(
    index = 400,
    label = "CHCHOH + O2 <=> HOCHO + H + CO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.112e+22, 'cm^3/(mol*s)'),
        n = -2.498,
        Ea = (20265.6, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CHCHOH + O2 <=> HOCHO + H + CO""",
)

entry(
    index = 401,
    label = "CHCHOH + CH2O <=> CH2CHOH + HCO",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(
                A = (1.549e+11, 'cm^3/(mol*s)'),
                n = 0.413,
                Ea = (8163.9, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (5.117e+13, 'cm^3/(mol*s)'),
                n = -0.005,
                Ea = (14641.6, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CHCHOH + CH2O <=> CH2CHOH + HCO""",
)

entry(
    index = 402,
    label = "CHCHOH + CH2O <=> CH2CHOH + H + CO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.232e+13, 'cm^3/(mol*s)'),
        n = 0.337,
        Ea = (25787.3, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CHCHOH + CH2O <=> CH2CHOH + H + CO""",
)

entry(
    index = 403,
    label = "CHCHOH + HCO <=> CH2CHOH + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CHCHOH + HCO <=> CH2CHOH + CO""",
)

entry(
    index = 404,
    label = "CHCHOH + CH3 <=> HCCOH + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CHCHOH + CH3 <=> HCCOH + CH4""",
)

entry(
    index = 405,
    label = "cC2H3O <=> CH2CHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.7e+31, 's^-1'), n=-6.9, Ea=(14994, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC2H3O <=> CH2CHO""",
)

entry(
    index = 406,
    label = "cC2H3O <=> CH2CO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 's^-1'), n=0, Ea=(14863, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC2H3O <=> CH2CO + H""",
)

entry(
    index = 407,
    label = "cC2H3O <=> CH3 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.1e+12, 's^-1'), n=0, Ea=(14280, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC2H3O <=> CH3 + CO""",
)

entry(
    index = 408,
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
    index = 409,
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
    index = 410,
    label = "CH3CO + H <=> CH3 + HCO",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(2.1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (-1.257e+23, 'cm^3/(mol*s)'),
                n = -2.473,
                Ea = (19927.3, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH3CO + H <=> CH3 + HCO""",
)

entry(
    index = 411,
    label = "CH3CO + H <=> CH3 + H + CO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.257e+23, 'cm^3/(mol*s)'),
        n = -2.473,
        Ea = (19927.3, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3CO + H <=> CH3 + H + CO""",
)

entry(
    index = 412,
    label = "CH3CO + H <=> CH2CO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CO + H <=> CH2CO + H2""",
)

entry(
    index = 413,
    label = "CH3CO + O <=> CH3 + CO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.6e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CO + O <=> CH3 + CO2""",
)

entry(
    index = 414,
    label = "CH3CO + O <=> CH2CO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CO + O <=> CH2CO + OH""",
)

entry(
    index = 415,
    label = "CH3CO + OH <=> CH2CO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CO + OH <=> CH2CO + H2O""",
)

entry(
    index = 416,
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
    index = 417,
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
    index = 418,
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
    index = 419,
    label = "CH3CO + CH3 <=> C2H6 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CO + CH3 <=> C2H6 + CO""",
)

entry(
    index = 420,
    label = "CH3CO + CH3 <=> CH2CO + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CO + CH3 <=> CH2CO + CH4""",
)

entry(
    index = 421,
    label = "CH3CO + CH3OO <=> CH3 + CO2 + CH3O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CO + CH3OO <=> CH3 + CO2 + CH3O""",
)

entry(
    index = 422,
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
    index = 423,
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
    index = 424,
    label = "CH2CHO + H <=> CH3 + HCO",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.1, 1, 10, 100], 'atm'),
        arrhenius = [
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (1.9e+17, 'cm^3/(mol*s)'),
                        n = -0.935,
                        Ea = (3120, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (-1.137e+27, 'cm^3/(mol*s)'),
                        n = -3.408,
                        Ea = (23047.3, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (2.5e+18, 'cm^3/(mol*s)'),
                        n = -1.243,
                        Ea = (4062, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (-1.496e+28, 'cm^3/(mol*s)'),
                        n = -3.716,
                        Ea = (23989.3, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (1.9e+22, 'cm^3/(mol*s)'),
                        n = -2.3,
                        Ea = (7693, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (-1.137e+32, 'cm^3/(mol*s)'),
                        n = -4.773,
                        Ea = (27620.3, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (2.8e+25, 'cm^3/(mol*s)'),
                        n = -3.1,
                        Ea = (12454, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (-1.676e+35, 'cm^3/(mol*s)'),
                        n = -5.573,
                        Ea = (32381.3, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (7.5e+20, 'cm^3/(mol*s)'),
                        n = -1.693,
                        Ea = (13429, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (-4.489e+30, 'cm^3/(mol*s)'),
                        n = -4.166,
                        Ea = (33356.3, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHO + H <=> CH3 + HCO""",
)

entry(
    index = 425,
    label = "CH2CHO + H <=> CH3 + H + CO",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.1, 1, 10, 100], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (1.137e+27, 'cm^3/(mol*s)'),
                n = -3.408,
                Ea = (23047.3, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.496e+28, 'cm^3/(mol*s)'),
                n = -3.716,
                Ea = (23989.3, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.137e+32, 'cm^3/(mol*s)'),
                n = -4.773,
                Ea = (27620.3, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.676e+35, 'cm^3/(mol*s)'),
                n = -5.573,
                Ea = (32381.3, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (4.489e+30, 'cm^3/(mol*s)'),
                n = -4.166,
                Ea = (33356.3, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHO + H <=> CH3 + H + CO""",
)

entry(
    index = 426,
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
    index = 427,
    label = "CH2CHO + O <=> CH2O + HCO",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (-2.993e+23, 'cm^3/(mol*s)'),
                n = -2.473,
                Ea = (19927.3, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHO + O <=> CH2O + HCO""",
)

entry(
    index = 428,
    label = "CH2CHO + O <=> CH2O + H + CO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.993e+23, 'cm^3/(mol*s)'),
        n = -2.473,
        Ea = (19927.3, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHO + O <=> CH2O + H + CO""",
)

entry(
    index = 429,
    label = "CH2CHO + OH <=> CH2CO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHO + OH <=> CH2CO + H2O""",
)

entry(
    index = 430,
    label = "CH2CHO + OH <=> CH2OH + HCO",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (-5.986e+22, 'cm^3/(mol*s)'),
                n = -2.473,
                Ea = (19927.3, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHO + OH <=> CH2OH + HCO""",
)

entry(
    index = 431,
    label = "CH2CHO + OH <=> CH2OH + H + CO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.986e+22, 'cm^3/(mol*s)'),
        n = -2.473,
        Ea = (19927.3, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHO + OH <=> CH2OH + H + CO""",
)

entry(
    index = 432,
    label = "CH2CHO + HO2 <=> CH2O + HCO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHO + HO2 <=> CH2O + HCO + OH""",
)

entry(
    index = 433,
    label = "CH2CHO + O2 <=> CH2O + CO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+10, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHO + O2 <=> CH2O + CO + OH""",
)

entry(
    index = 434,
    label = "CH2CHO + CH2 <=> C2H4 + HCO",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (-2.993e+23, 'cm^3/(mol*s)'),
                n = -2.473,
                Ea = (19927.3, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHO + CH2 <=> C2H4 + HCO""",
)

entry(
    index = 435,
    label = "CH2CHO + CH2 <=> C2H4 + H + CO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.993e+23, 'cm^3/(mol*s)'),
        n = -2.473,
        Ea = (19927.3, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHO + CH2 <=> C2H4 + H + CO""",
)

entry(
    index = 436,
    label = "CH2CHO + CH <=> C2H3 + HCO",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (-5.986e+23, 'cm^3/(mol*s)'),
                n = -2.473,
                Ea = (19927.3, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHO + CH <=> C2H3 + HCO""",
)

entry(
    index = 437,
    label = "CH2CHO + CH <=> C2H3 + H + CO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.986e+23, 'cm^3/(mol*s)'),
        n = -2.473,
        Ea = (19927.3, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHO + CH <=> C2H3 + H + CO""",
)

entry(
    index = 438,
    label = "CHCHO + H <=> CH2CO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CHCHO + H <=> CH2CO + H""",
)

entry(
    index = 439,
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
    index = 440,
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
    index = 441,
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
    index = 442,
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
    index = 443,
    label = "CH2CO + H <=> HCCO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+07, 'cm^3/(mol*s)'), n=2, Ea=(10000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CO + H <=> HCCO + H2""",
)

entry(
    index = 444,
    label = "CH2CO + O <=> CO2 + CH2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(1350, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CO + O <=> CO2 + CH2""",
)

entry(
    index = 445,
    label = "CH2CO + OH <=> CH2OH + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1013, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CO + OH <=> CH2OH + CO""",
)

entry(
    index = 446,
    label = "CH2CO + OH <=> CH3 + CO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.7e+11, 'cm^3/(mol*s)'), n=0, Ea=(-1013, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CO + OH <=> CH3 + CO2""",
)

entry(
    index = 447,
    label = "CH2CO + OH <=> HCCO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+07, 'cm^3/(mol*s)'), n=2, Ea=(3000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CO + OH <=> HCCO + H2O""",
)

entry(
    index = 448,
    label = "CH2CO + CH2(S) <=> C2H4 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.6e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CO + CH2(S) <=> C2H4 + CO""",
)

entry(
    index = 449,
    label = "HCCOH + H <=> HCCO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+07, 'cm^3/(mol*s)'), n=2, Ea=(1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCCOH + H <=> HCCO + H2""",
)

entry(
    index = 450,
    label = "HCCOH + OH <=> HCCO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+07, 'cm^3/(mol*s)'), n=2, Ea=(1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCCOH + OH <=> HCCO + H2O""",
)

entry(
    index = 451,
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
    index = 452,
    label = "HCCO + H <=> CH2(S) + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.5e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCCO + H <=> CH2(S) + CO""",
)

entry(
    index = 453,
    label = "HCCO + O <=> CO + CO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCCO + O <=> CO + CO + H""",
)

entry(
    index = 454,
    label = "HCCO + OH <=> C2O + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(14390, 'cm^3/(mol*s)'), n=2.65, Ea=(1472, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCCO + OH <=> C2O + H2O""",
)

entry(
    index = 455,
    label = "HCCO + OH <=> CH2CO + O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.265e+06, 'cm^3/(mol*s)'),
        n = 1.99,
        Ea = (11280, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is HCCO + OH <=> CH2CO + O""",
)

entry(
    index = 456,
    label = "HCCO + OH <=> HCOH + CO",
    degeneracy = 1,
    duplicate = True,
    kinetics = PDepArrhenius(
        pressures = ([1, 10, 100], 'atm'),
        arrhenius = [
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (3.032e+16, 'cm^3/(mol*s)'),
                        n = -0.935,
                        Ea = (659.4, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (8.694e+19, 'cm^3/(mol*s)'),
                        n = -1.792,
                        Ea = (5994.3, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (1.143e+18, 'cm^3/(mol*s)'),
                        n = -1.392,
                        Ea = (1395.1, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (3.49e+22, 'cm^3/(mol*s)'),
                        n = -2.475,
                        Ea = (9162.6, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(
                        A = (3.224e+18, 'cm^3/(mol*s)'),
                        n = -1.523,
                        Ea = (1626.7, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.261e+24, 'cm^3/(mol*s)'),
                        n = -2.902,
                        Ea = (10522.1, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is HCCO + OH <=> HCOH + CO""",
)

entry(
    index = 457,
    label = "HCCO + OH <=> HCOH + CO",
    degeneracy = 1,
    duplicate = True,
    kinetics = Arrhenius(
        A = (2.873e+12, 'cm^3/(mol*s)'),
        n = 0.37,
        Ea = (-24, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is HCCO + OH <=> HCOH + CO""",
)

entry(
    index = 458,
    label = "HCCO + OH <=> CH2O + CO",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([1, 100], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (1.187e+21, 'cm^3/(mol*s)'),
                n = -2.459,
                Ea = (2527.6, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(A=(1.1e+08, 'cm^3/(mol*s)'), n=0.11, Ea=(52, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is HCCO + OH <=> CH2O + CO""",
)

entry(
    index = 459,
    label = "HCCO + OH <=> OCHCO + H",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.1, 1, 10, 100], 'atm'),
        arrhenius = [
            Arrhenius(
                A = (2.632e+08, 'cm^3/(mol*s)'),
                n = 1.41,
                Ea = (845, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.632e+08, 'cm^3/(mol*s)'),
                n = 1.41,
                Ea = (845, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.638e+08, 'cm^3/(mol*s)'),
                n = 1.41,
                Ea = (849, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.963e+08, 'cm^3/(mol*s)'),
                n = 1.4,
                Ea = (917, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (8.19e+08, 'cm^3/(mol*s)'),
                n = 1.28,
                Ea = (1531, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is HCCO + OH <=> OCHCO + H""",
)

entry(
    index = 460,
    label = "HCCO + OH <=> CO2 + CH2",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiPDepArrhenius(
        arrhenius = [
            PDepArrhenius(
                pressures = ([0.1, 1, 10], 'atm'),
                arrhenius = [
                    MultiArrhenius(
                        arrhenius = [
                            Arrhenius(
                                A = (1.698e+15, 'cm^3/(mol*s)'),
                                n = -1.19,
                                Ea = (-521, 'cal/mol'),
                                T0 = (1, 'K'),
                            ),
                            Arrhenius(
                                A = (-7.407e+17, 'cm^3/(mol*s)'),
                                n = -1.92,
                                Ea = (1686, 'cal/mol'),
                                T0 = (1, 'K'),
                            ),
                        ],
                    ),
                    MultiArrhenius(
                        arrhenius = [
                            Arrhenius(
                                A = (7.292e+27, 'cm^3/(mol*s)'),
                                n = -5.023,
                                Ea = (2468, 'cal/mol'),
                                T0 = (1, 'K'),
                            ),
                            Arrhenius(
                                A = (1.116e+21, 'cm^3/(mol*s)'),
                                n = -2.28,
                                Ea = (16960.4, 'cal/mol'),
                                T0 = (1, 'K'),
                            ),
                        ],
                    ),
                    MultiArrhenius(
                        arrhenius = [
                            Arrhenius(
                                A = (5.974e+15, 'cm^3/(mol*s)'),
                                n = -0.64,
                                Ea = (363, 'cal/mol'),
                                T0 = (1, 'K'),
                            ),
                            Arrhenius(
                                A = (-2.577e+19, 'cm^3/(mol*s)'),
                                n = -1.64,
                                Ea = (3539, 'cal/mol'),
                                T0 = (1, 'K'),
                            ),
                        ],
                    ),
                ],
            ),
            PDepArrhenius(
                pressures = ([0.01, 0.1, 1, 10], 'atm'),
                arrhenius = [
                    Arrhenius(
                        A = (1.018e+19, 'cm^3/(mol*s)'),
                        n = -2.08,
                        Ea = (44, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.397e+19, 'cm^3/(mol*s)'),
                        n = -2.12,
                        Ea = (88, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (7.106e+19, 'cm^3/(mol*s)'),
                        n = -2.3,
                        Ea = (824, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                    Arrhenius(
                        A = (1.789e+20, 'cm^3/(mol*s)'),
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
    index = 461,
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
    index = 462,
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
    index = 463,
    label = "HCCO + O2 <=> HCO + CO + O",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(220, 'cm^3/(mol*s)'), n=2.69, Ea=(3540, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (-1.317e+12, 'cm^3/(mol*s)'),
                n = 0.217,
                Ea = (23467.3, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is HCCO + O2 <=> HCO + CO + O""",
)

entry(
    index = 464,
    label = "HCCO + O2 <=> H + CO + CO + O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.317e+12, 'cm^3/(mol*s)'),
        n = 0.217,
        Ea = (23467.3, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is HCCO + O2 <=> H + CO + CO + O""",
)

entry(
    index = 465,
    label = "HCCO + CH2 <=> C2H3 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCCO + CH2 <=> C2H3 + CO""",
)

entry(
    index = 466,
    label = "HCCO + CH <=> C2H2 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCCO + CH <=> C2H2 + CO""",
)

entry(
    index = 467,
    label = "HCCO + HCCO <=> C2H2 + CO + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCCO + HCCO <=> C2H2 + CO + CO""",
)

entry(
    index = 468,
    label = "C2O <=> C + CO",
    degeneracy = 1,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(A=(2e+15, 'cm^3/(mol*s)'), n=0, Ea=(44200, 'cal/mol'), T0=(1, 'K')),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is C2O <=> C + CO""",
)

entry(
    index = 469,
    label = "C2O + H <=> CH + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2O + H <=> CH + CO""",
)

entry(
    index = 470,
    label = "C2O + O <=> CO + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2O + O <=> CO + CO""",
)

entry(
    index = 471,
    label = "C2O + OH <=> CO + CO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2O + OH <=> CO + CO + H""",
)

entry(
    index = 472,
    label = "C2O + O2 <=> CO + CO + O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(2600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2O + O2 <=> CO + CO + O""",
)

entry(
    index = 473,
    label = "C2O + O2 <=> CO + CO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(2600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2O + O2 <=> CO + CO2""",
)

entry(
    index = 474,
    label = "C2O + C <=> CO + C2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2O + C <=> CO + C2""",
)

entry(
    index = 475,
    label = "CH3CH2OOH <=> CH3CH2O + OH",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.1, 1, 100], 'atm'),
        arrhenius = [
            Arrhenius(A=(6.05e+58, 's^-1'), n=-14.05, Ea=(54131, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(9.26e+52, 's^-1'), n=-11.91, Ea=(53378, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.38e+33, 's^-1'), n=-5.27, Ea=(48696, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH3CH2OOH <=> CH3CH2O + OH""",
)

entry(
    index = 476,
    label = "CH3CH2OOH + H <=> CH3CHOOH + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.5e+10, 'cm^3/(mol*s)'), n=0, Ea=(1860, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OOH + H <=> CH3CHOOH + H2""",
)

entry(
    index = 477,
    label = "CH3CH2OOH + H <=> CH3CH2OO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.3e+10, 'cm^3/(mol*s)'), n=0, Ea=(1860, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OOH + H <=> CH3CH2OO + H2""",
)

entry(
    index = 478,
    label = "CH3CH2OOH + H <=> CH3CH2O + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+10, 'cm^3/(mol*s)'), n=0, Ea=(1860, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OOH + H <=> CH3CH2O + H2O""",
)

entry(
    index = 479,
    label = "CH3CH2OOH + O <=> CH3CHOOH + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(4750, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OOH + O <=> CH3CHOOH + OH""",
)

entry(
    index = 480,
    label = "CH3CH2OOH + O <=> CH3CH2OO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.7e+12, 'cm^3/(mol*s)'), n=0, Ea=(4750, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OOH + O <=> CH3CH2OO + OH""",
)

entry(
    index = 481,
    label = "CH3CH2OOH + OH <=> CH3CHOOH + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.2e+11, 'cm^3/(mol*s)'), n=0, Ea=(-258, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OOH + OH <=> CH3CHOOH + H2O""",
)

entry(
    index = 482,
    label = "CH3CH2OOH + OH <=> CH3CH2OO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.1e+12, 'cm^3/(mol*s)'), n=0, Ea=(-437, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OOH + OH <=> CH3CH2OO + H2O""",
)

entry(
    index = 483,
    label = "CH3CH2OOH + HO2 <=> CH3CH2OO + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(41000, 'cm^3/(mol*s)'), n=2.5, Ea=(10206, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OOH + HO2 <=> CH3CH2OO + H2O2""",
)

entry(
    index = 484,
    label = "CH3CHOOH <=> CH3CHO + OH",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([1, 10, 100], 'atm'),
        arrhenius = [
            Arrhenius(A=(3.5e+12, 's^-1'), n=-0.947, Ea=(979, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(3.5e+13, 's^-1'), n=-0.947, Ea=(980, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(5.75e+14, 's^-1'), n=-1.012, Ea=(1068, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH3CHOOH <=> CH3CHO + OH""",
)

entry(
    index = 485,
    label = "CH3CH2OO <=> CH2CH2OOH",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100], 'atm'),
        arrhenius = [
            Arrhenius(A=(4.29e+08, 's^-1'), n=-0.97, Ea=(22110, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2.77e+07, 's^-1'), n=-0.52, Ea=(21990, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(351000, 's^-1'), n=0.16, Ea=(21700, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1580, 's^-1'), n=0.99, Ea=(21290, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(0.651, 's^-1'), n=2.17, Ea=(20650, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(430, 's^-1'), n=1.31, Ea=(21600, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(0.00732, 's^-1'), n=2.96, Ea=(20660, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(3.79e-08, 's^-1'), n=4.77, Ea=(19590, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(7.38e-08, 's^-1'), n=4.78, Ea=(20030, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(8.17e-08, 's^-1'), n=4.87, Ea=(20420, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2.27e-14, 's^-1'), n=7.09, Ea=(18990, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH3CH2OO <=> CH2CH2OOH""",
)

entry(
    index = 486,
    label = "CH3CH2OO + H <=> CH3CH2O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OO + H <=> CH3CH2O + OH""",
)

entry(
    index = 487,
    label = "CH3CH2OO + O <=> CH3CH2O + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.85e+10, 'cm^3/(mol*s)'), n=1, Ea=(-724, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OO + O <=> CH3CH2O + O2""",
)

entry(
    index = 488,
    label = "CH3CH2OO + OH <=> CH3CH2OH + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OO + OH <=> CH3CH2OH + O2""",
)

entry(
    index = 489,
    label = "CH3CH2OO + HO2 <=> CH3CH2OOH + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.5e+11, 'cm^3/(mol*s)'), n=0, Ea=(-1391, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OO + HO2 <=> CH3CH2OOH + O2""",
)

entry(
    index = 490,
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
    index = 491,
    label = "CH3CH2OO + CH3 <=> CH3CH2O + CH3O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.1e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1411, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OO + CH3 <=> CH3CH2O + CH3O""",
)

entry(
    index = 492,
    label = "CH3CH2OO + CH4 <=> CH3CH2OOH + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(47000, 'cm^3/(mol*s)'), n=2.5, Ea=(21000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OO + CH4 <=> CH3CH2OOH + CH3""",
)

entry(
    index = 493,
    label = "CH3CH2OO + CH3OH <=> CH3CH2OOH + CH2OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(19400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OO + CH3OH <=> CH3CH2OOH + CH2OH""",
)

entry(
    index = 494,
    label = "CH3CH2OO + CH2O <=> CH3CH2OOH + HCO",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(41000, 'cm^3/(mol*s)'), n=2.5, Ea=(10206, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (-2.454e+14, 'cm^3/(mol*s)'),
                n = 0.027,
                Ea = (30133.3, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH3CH2OO + CH2O <=> CH3CH2OOH + HCO""",
)

entry(
    index = 495,
    label = "CH3CH2OO + CH2O <=> CH3CH2OOH + H + CO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.454e+14, 'cm^3/(mol*s)'),
        n = 0.027,
        Ea = (30133.3, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3CH2OO + CH2O <=> CH3CH2OOH + H + CO""",
)

entry(
    index = 496,
    label = "CH3CH2OO + C2H5 <=> CH3CH2O + CH3CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.1e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1411, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OO + C2H5 <=> CH3CH2O + CH3CH2O""",
)

entry(
    index = 497,
    label = "CH3CH2OO + C2H6 <=> CH3CH2OOH + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.6, 'cm^3/(mol*s)'), n=3.76, Ea=(17200, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OO + C2H6 <=> CH3CH2OOH + C2H5""",
)

entry(
    index = 498,
    label = "CH3CH2OO + CH3CHO <=> CH3CH2OOH + CH3CO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.4e+19, 'cm^3/(mol*s)'),
        n = -2.2,
        Ea = (14030, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3CH2OO + CH3CHO <=> CH3CH2OOH + CH3CO""",
)

entry(
    index = 499,
    label = "CH3CH2OO + CH3CHO <=> CH3CH2OOH + CH2CHO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.3e+11, 'cm^3/(mol*s)'),
        n = 0.4,
        Ea = (14864, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3CH2OO + CH3CHO <=> CH3CH2OOH + CH2CHO""",
)

entry(
    index = 500,
    label = "CH3CH2OO + CH3CH2OO <=> CH3CH2O + CH3CH2O + O2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.9e+11, 'cm^3/(mol*s)'),
        n = -0.27,
        Ea = (408, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3CH2OO + CH3CH2OO <=> CH3CH2O + CH3CH2O + O2""",
)

entry(
    index = 501,
    label = "CH3CH2OO + CH3CH2OO <=> CH3CHO + CH3CH2OH + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.3e+09, 'cm^3/(mol*s)'), n=0, Ea=(-850, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CH2OO + CH3CH2OO <=> CH3CHO + CH3CH2OH + O2""",
)

entry(
    index = 502,
    label = "CH2CH2OOH <=> cC2H4O + OH",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01316, 0.98692, 9.86923, 98.6923, 1000], 'atm'),
        arrhenius = [
            Arrhenius(A=(3.342e+27, 's^-1'), n=-6.117, Ea=(15373.9, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(4.361e+28, 's^-1'), n=-5.83, Ea=(17202.1, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2.689e+32, 's^-1'), n=-6.633, Ea=(20310.7, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(7.544e+35, 's^-1'), n=-7.331, Ea=(23906.8, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(6.672e+10, 's^-1'), n=0.637, Ea=(15974.2, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CH2OOH <=> cC2H4O + OH""",
)

entry(
    index = 503,
    label = "CH2CHOOH <=> CH2CHO + OH",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([1, 10, 50, 100], 'atm'),
        arrhenius = [
            Arrhenius(A=(2e+35, 's^-1'), n=-6.7, Ea=(47450, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.12e+28, 's^-1'), n=-4.15, Ea=(46190, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2.8e+26, 's^-1'), n=-3.5, Ea=(46340, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2.22e+17, 's^-1'), n=-0.42, Ea=(44622, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHOOH <=> CH2CHO + OH""",
)

entry(
    index = 504,
    label = "CH2CHOOH + H <=> CH2CHOO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.3e+10, 'cm^3/(mol*s)'), n=0, Ea=(1860, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHOOH + H <=> CH2CHOO + H2""",
)

entry(
    index = 505,
    label = "CH2CHOOH + H <=> CH2CHO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+10, 'cm^3/(mol*s)'), n=0, Ea=(1860, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHOOH + H <=> CH2CHO + H2O""",
)

entry(
    index = 506,
    label = "CH2CHOOH + O <=> CH2CHOO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.7e+12, 'cm^3/(mol*s)'), n=0, Ea=(4750, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHOOH + O <=> CH2CHOO + OH""",
)

entry(
    index = 507,
    label = "CH2CHOOH + OH <=> CH2CHOO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.1e+12, 'cm^3/(mol*s)'), n=0, Ea=(-437, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHOOH + OH <=> CH2CHOO + H2O""",
)

entry(
    index = 508,
    label = "CH2CHOOH + HO2 <=> CH2CHOO + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(41000, 'cm^3/(mol*s)'), n=2.5, Ea=(10206, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHOOH + HO2 <=> CH2CHOO + H2O2""",
)

entry(
    index = 509,
    label = "CH2CHOO <=> CHCHO + OH",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
        arrhenius = [
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(3.64e+49, 's^-1'), n=-12.13, Ea=(67420, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.17e+56, 's^-1'), n=-14.81, Ea=(60700, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(1.44e+36, 's^-1'), n=-9.92, Ea=(41220, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(2.32e+40, 's^-1'), n=-9.39, Ea=(50420, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(4.18e+40, 's^-1'), n=-10.53, Ea=(43670, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.61e+43, 's^-1'), n=-9.99, Ea=(50290, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(3.79e+46, 's^-1'), n=-10.72, Ea=(51900, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(2.33e+124, 's^-1'), n=-36.77, Ea=(70100, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(1.6e+49, 's^-1'), n=-11.24, Ea=(54150, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.88e+103, 's^-1'), n=-29.49, Ea=(65410, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(2.38e+51, 's^-1'), n=-11.64, Ea=(56980, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(5.96e+86, 's^-1'), n=-23.81, Ea=(62170, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(2e+54, 's^-1'), n=-12.22, Ea=(61840, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.51e+57, 's^-1'), n=-13.94, Ea=(55390, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(9.54e+195, 's^-1'), n=-52.27, Ea=(163500, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.79e+34, 's^-1'), n=-6.4, Ea=(50000, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHOO <=> CHCHO + OH""",
)

entry(
    index = 510,
    label = "CH2CHOO <=> CH2CHO + O",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
        arrhenius = [
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(2.7e+180, 's^-1'), n=-48.19, Ea=(169300, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.47e+30, 's^-1'), n=-6.64, Ea=(41110, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(3.9e+38, 's^-1'), n=-8.69, Ea=(42770, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(9.65e-12, 's^-1'), n=5.96, Ea=(22890, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(4.57e+47, 's^-1'), n=-11.21, Ea=(47050, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(3.95e+22, 's^-1'), n=-3.71, Ea=(36270, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(7.62e+81, 's^-1'), n=-21.28, Ea=(65080, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(2.39e+33, 's^-1'), n=-6.62, Ea=(41280, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(1.86e+68, 's^-1'), n=-16.83, Ea=(60680, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(6.37e+31, 's^-1'), n=-5.96, Ea=(41260, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(2.02e+55, 's^-1'), n=-12.69, Ea=(55840, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(2.13e+29, 's^-1'), n=-5.1, Ea=(40710, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(1.11e+53, 's^-1'), n=-11.79, Ea=(56690, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(4.66e+27, 's^-1'), n=-4.5, Ea=(40530, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(4.3e+48, 's^-1'), n=-10.31, Ea=(56090, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(5.99e+25, 's^-1'), n=-3.85, Ea=(40120, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHOO <=> CH2CHO + O""",
)

entry(
    index = 511,
    label = "CH2CHOO <=> OCHCHO + H",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
        arrhenius = [
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(6.41e+80, 's^-1'), n=-22.2, Ea=(51750, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.19e+28, 's^-1'), n=-6.01, Ea=(28740, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(3.31e+65, 's^-1'), n=-17.01, Ea=(48090, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.4e+25, 's^-1'), n=-4.8, Ea=(28940, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(5.98e+51, 's^-1'), n=-12.62, Ea=(43000, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(2.91e+20, 's^-1'), n=-3.29, Ea=(27550, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(1.48e+44, 's^-1'), n=-10.12, Ea=(40790, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.58e+19, 's^-1'), n=-2.82, Ea=(27620, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(1.26e+59, 's^-1'), n=-14.33, Ea=(51390, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.93e+22, 's^-1'), n=-3.54, Ea=(29980, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(4.93e+26, 's^-1'), n=-4.67, Ea=(34320, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(7.51e+29, 's^-1'), n=-5.75, Ea=(34490, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(2.06e+33, 's^-1'), n=-6.38, Ea=(39520, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(7.14e+61, 's^-1'), n=-16.16, Ea=(43280, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(1.3e+32, 's^-1'), n=-5.92, Ea=(40660, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.14e+19, 's^-1'), n=-2.56, Ea=(29670, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHOO <=> OCHCHO + H""",
)

entry(
    index = 512,
    label = "CH2CHOO <=> CH2CO + OH",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
        arrhenius = [
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(1.15e+47, 's^-1'), n=-12.28, Ea=(75330, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(231, 's^-1'), n=-0.73, Ea=(25710, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(8.43e+09, 's^-1'), n=-2.06, Ea=(33720, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.83e-23, 's^-1'), n=7.84, Ea=(20190, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(60600, 's^-1'), n=0.17, Ea=(34220, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(3.82e+63, 's^-1'), n=-20.44, Ea=(43420, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(1.51e+19, 's^-1'), n=-3.61, Ea=(43060, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(3.18e+27, 's^-1'), n=-7.76, Ea=(37230, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(2.13e+33, 's^-1'), n=-7.39, Ea=(51610, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(2.32e-05, 's^-1'), n=3.47, Ea=(31560, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(4.44e+36, 's^-1'), n=-7.99, Ea=(54680, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(0.106, 's^-1'), n=2.64, Ea=(34160, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(1.19e+37, 's^-1'), n=-7.8, Ea=(56460, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(562, 's^-1'), n=1.7, Ea=(36450, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(9.08e+35, 's^-1'), n=-7.21, Ea=(57550, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.11e+07, 's^-1'), n=0.52, Ea=(38670, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHOO <=> CH2CO + OH""",
)

entry(
    index = 513,
    label = "CH2CHOO <=> CH2O + HCO",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
        arrhenius = [
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(1.66e+174, 's^-1'), n=-55.52, Ea=(60320, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(2.27e+35, 's^-1'), n=-7.97, Ea=(31280, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(9.03e+66, 's^-1'), n=-17.25, Ea=(48120, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(2.08e+26, 's^-1'), n=-4.96, Ea=(28780, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(1.82e+43, 's^-1'), n=-9.87, Ea=(37960, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.45e+20, 's^-1'), n=-3.08, Ea=(26630, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(8.64e+33, 's^-1'), n=-6.88, Ea=(34370, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.06e+130, 's^-1'), n=-39.38, Ea=(54700, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(7.29e+171, 's^-1'), n=-43.53, Ea=(191900, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(2.35e+34, 's^-1'), n=-6.87, Ea=(35700, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(1.03e+32, 's^-1'), n=-6.06, Ea=(35500, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(2.18e+175, 's^-1'), n=-53.78, Ea=(68500, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(1.85e+34, 's^-1'), n=-6.57, Ea=(38510, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.07e+185, 's^-1'), n=-54.22, Ea=(88990, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(5.7e+29, 's^-1'), n=-5.19, Ea=(36800, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(468, 's^-1'), n=1.81, Ea=(18100, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHOO <=> CH2O + HCO""",
)

entry(
    index = 514,
    label = "CH2CHOO <=> CH2O + H + CO",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
        arrhenius = [
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(3.88e+174, 's^-1'), n=-55.52, Ea=(60320, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(5.29e+35, 's^-1'), n=-7.97, Ea=(31280, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(2.11e+67, 's^-1'), n=-17.25, Ea=(48120, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(4.85e+26, 's^-1'), n=-4.96, Ea=(28780, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(4.26e+43, 's^-1'), n=-9.87, Ea=(37960, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(3.37e+20, 's^-1'), n=-3.08, Ea=(26630, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(2.02e+34, 's^-1'), n=-6.88, Ea=(34370, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(2.46e+130, 's^-1'), n=-39.38, Ea=(54700, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(1.7e+172, 's^-1'), n=-43.53, Ea=(191900, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(5.49e+34, 's^-1'), n=-6.87, Ea=(35700, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(2.4e+32, 's^-1'), n=-6.06, Ea=(35500, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(5.09e+175, 's^-1'), n=-53.78, Ea=(68500, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(4.32e+34, 's^-1'), n=-6.57, Ea=(38510, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(2.49e+185, 's^-1'), n=-54.22, Ea=(88990, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(1.33e+30, 's^-1'), n=-5.19, Ea=(36800, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1090, 's^-1'), n=1.81, Ea=(18100, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHOO <=> CH2O + H + CO""",
)

entry(
    index = 515,
    label = "CH2CHOO <=> CO + CH3O",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
        arrhenius = [
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(5.2e+33, 's^-1'), n=-7.92, Ea=(31320, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(2.31e+129, 's^-1'), n=-41.86, Ea=(45850, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(1.26e+98, 's^-1'), n=-27.09, Ea=(64060, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(2.42e+28, 's^-1'), n=-5.99, Ea=(30540, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(1.8e+33, 's^-1'), n=-7.27, Ea=(33760, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(8.69e-50, 's^-1'), n=16.63, Ea=(-3900, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(3.83e+33, 's^-1'), n=-7.2, Ea=(35100, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.19e-39, 's^-1'), n=13.61, Ea=(-1317, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(1.28e+79, 's^-1'), n=-19.61, Ea=(74870, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(8.8e+86, 's^-1'), n=-23.08, Ea=(61010, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(4.07e+32, 's^-1'), n=-6.62, Ea=(37210, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1270, 's^-1'), n=1.44, Ea=(18660, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(6.86e+44, 's^-1'), n=-10.04, Ea=(47030, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.97e+17, 's^-1'), n=-2.23, Ea=(28590, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(-10700, 's^-1'), n=1.33, Ea=(15620, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.16e-07, 's^-1'), n=4.81, Ea=(12010, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHOO <=> CO + CH3O""",
)

entry(
    index = 516,
    label = "CH2CHOO <=> CO2 + CH3",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.1, 0.316, 1, 3.16, 10, 31.6, 100], 'atm'),
        arrhenius = [
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(5.09e+33, 's^-1'), n=-7.95, Ea=(31290, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(4.2e+122, 's^-1'), n=-39.75, Ea=(43640, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(1.21e+118, 's^-1'), n=-33.13, Ea=(73790, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.96e+29, 's^-1'), n=-6.29, Ea=(30920, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(8.56e+32, 's^-1'), n=-7.21, Ea=(33550, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(5.1e-66, 's^-1'), n=21.37, Ea=(-11110, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(3.27e+33, 's^-1'), n=-7.22, Ea=(34990, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(1.76e-47, 's^-1'), n=15.85, Ea=(-5283, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(3.49e-79, 's^-1'), n=25.01, Ea=(-21020, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(3.82e+32, 's^-1'), n=-6.8, Ea=(35690, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(8.16e+32, 's^-1'), n=-6.76, Ea=(37270, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(4.62, 's^-1'), n=2.1, Ea=(17170, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(7.01e+37, 's^-1'), n=-8.06, Ea=(42200, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(3.49e+14, 's^-1'), n=-1.58, Ea=(26470, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(-2510, 's^-1'), n=1.41, Ea=(14420, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(A=(4.05e-09, 's^-1'), n=5.14, Ea=(10480, 'cal/mol'), T0=(1, 'K')),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHOO <=> CO2 + CH3""",
)

entry(
    index = 517,
    label = "CH2CHOO + H <=> CH2CHO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHOO + H <=> CH2CHO + OH""",
)

entry(
    index = 518,
    label = "CH2CHOO + O <=> CH2CHO + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(-145, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHOO + O <=> CH2CHO + O2""",
)

entry(
    index = 519,
    label = "CH2CHOO + OH <=> CH2CHOH + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+15, 'cm^3/(mol*s)'), n=-0.6, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHOO + OH <=> CH2CHOH + O2""",
)

entry(
    index = 520,
    label = "CH2CHOO + OH <=> CH2CHO + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+11, 'cm^3/(mol*s)'), n=0.6, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHOO + OH <=> CH2CHO + HO2""",
)

entry(
    index = 521,
    label = "CH2CHOO + HO2 <=> CH2CHOOH + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.5e+11, 'cm^3/(mol*s)'), n=0, Ea=(-1391, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHOO + HO2 <=> CH2CHOOH + O2""",
)

entry(
    index = 522,
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
    index = 523,
    label = "CH2CHOO + CH3 <=> CH2CHO + CH3O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.1e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1411, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHOO + CH3 <=> CH2CHO + CH3O""",
)

entry(
    index = 524,
    label = "CH2CHOO + CH4 <=> CH2CHOOH + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(47000, 'cm^3/(mol*s)'), n=2.5, Ea=(21000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHOO + CH4 <=> CH2CHOOH + CH3""",
)

entry(
    index = 525,
    label = "CH2CHOO + CH3OH <=> CH2CHOOH + CH2OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(19400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHOO + CH3OH <=> CH2CHOOH + CH2OH""",
)

entry(
    index = 526,
    label = "CH2CHOO + CH2O <=> CH2CHOOH + HCO",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(41000, 'cm^3/(mol*s)'), n=2.5, Ea=(10206, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (-2.454e+14, 'cm^3/(mol*s)'),
                n = 0.027,
                Ea = (30133.3, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHOO + CH2O <=> CH2CHOOH + HCO""",
)

entry(
    index = 527,
    label = "CH2CHOO + CH2O <=> CH2CHOOH + H + CO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.454e+14, 'cm^3/(mol*s)'),
        n = 0.027,
        Ea = (30133.3, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHOO + CH2O <=> CH2CHOOH + H + CO""",
)

entry(
    index = 528,
    label = "CH2CHOO + C2H6 <=> CH2CHOOH + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.6, 'cm^3/(mol*s)'), n=3.76, Ea=(17200, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHOO + C2H6 <=> CH2CHOOH + C2H5""",
)

entry(
    index = 529,
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
    index = 530,
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
    index = 531,
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
    index = 532,
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
    index = 533,
    label = "OCHCHO + H <=> CH2O + HCO",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(5.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(4300, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (-3.232e+23, 'cm^3/(mol*s)'),
                n = -2.473,
                Ea = (24227.3, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is OCHCHO + H <=> CH2O + HCO""",
)

entry(
    index = 534,
    label = "OCHCHO + H <=> CH2O + H + CO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.232e+23, 'cm^3/(mol*s)'),
        n = -2.473,
        Ea = (24227.3, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is OCHCHO + H <=> CH2O + H + CO""",
)

entry(
    index = 535,
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
    index = 536,
    label = "OCHCHO + OH <=> OCHCO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+06, 'cm^3/(mol*s)'), n=2, Ea=(-1630, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is OCHCHO + OH <=> OCHCO + H2O""",
)

entry(
    index = 537,
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
    index = 538,
    label = "OCHCHO + HO2 <=> OCHCO + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(41000, 'cm^3/(mol*s)'), n=2.5, Ea=(10206, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is OCHCHO + HO2 <=> OCHCO + H2O2""",
)

entry(
    index = 539,
    label = "OCHCHO + O2 <=> OCHCO + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(240000, 'cm^3/(mol*s)'), n=2.5, Ea=(36461, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is OCHCHO + O2 <=> OCHCO + HO2""",
)

entry(
    index = 540,
    label = "OCHCO <=> HCO + CO",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.1, 1, 100], 'atm'),
        arrhenius = [
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(3.8e+12, 's^-1'), n=0, Ea=(8610, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (-7.988e+21, 's^-1'),
                        n = -2.359,
                        Ea = (27419.6, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(3.8e+13, 's^-1'), n=0, Ea=(8665, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (-2.275e+23, 's^-1'),
                        n = -2.473,
                        Ea = (28592.3, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(4.1e+14, 's^-1'), n=0, Ea=(8765, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (-2.454e+24, 's^-1'),
                        n = -2.473,
                        Ea = (28692.3, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
            MultiArrhenius(
                arrhenius = [
                    Arrhenius(A=(1.1e+14, 's^-1'), n=0.133, Ea=(10140, 'cal/mol'), T0=(1, 'K')),
                    Arrhenius(
                        A = (-1.391e+24, 's^-1'),
                        n = -2.419,
                        Ea = (30990.9, 'cal/mol'),
                        T0 = (1, 'K'),
                    ),
                ],
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is OCHCO <=> HCO + CO""",
)

entry(
    index = 541,
    label = "OCHCO <=> H + CO + CO",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01, 0.1, 1, 100], 'atm'),
        arrhenius = [
            Arrhenius(A=(7.988e+21, 's^-1'), n=-2.359, Ea=(27419.6, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2.275e+23, 's^-1'), n=-2.473, Ea=(28592.3, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2.454e+24, 's^-1'), n=-2.473, Ea=(28692.3, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.391e+24, 's^-1'), n=-2.419, Ea=(30990.9, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is OCHCO <=> H + CO + CO""",
)

entry(
    index = 542,
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
    index = 543,
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
    index = 544,
    label = "HOCH2O <=> HOCHO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 's^-1'), n=0, Ea=(14900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HOCH2O <=> HOCHO + H""",
)

entry(
    index = 545,
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
    index = 546,
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
    index = 547,
    label = "HOCHO + H <=> HOCO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(230, 'cm^3/(mol*s)'), n=3.272, Ea=(4858, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HOCHO + H <=> HOCO + H2""",
)

entry(
    index = 548,
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
    index = 549,
    label = "HOCHO + O <=> HOCO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(51, 'cm^3/(mol*s)'), n=3.422, Ea=(4216, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HOCHO + O <=> HOCO + OH""",
)

entry(
    index = 550,
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
    index = 551,
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
    index = 552,
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
    index = 553,
    label = "HOCHO + HO2 <=> HOCO + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.47, 'cm^3/(mol*s)'), n=3.975, Ea=(16787, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HOCHO + HO2 <=> HOCO + H2O2""",
)

entry(
    index = 554,
    label = "HOCHO + HO2 <=> OCHO + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(39, 'cm^3/(mol*s)'), n=3.08, Ea=(25206, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HOCHO + HO2 <=> OCHO + H2O2""",
)

entry(
    index = 555,
    label = "HOCO + HO2 <=> HOCHO + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HOCO + HO2 <=> HOCHO + O2""",
)

entry(
    index = 556,
    label = "HOCHO + O2 <=> OCHO + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(63000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HOCHO + O2 <=> OCHO + HO2""",
)

entry(
    index = 557,
    label = "OCHO <=> CO2 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+10, 's^-1'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is OCHO <=> CO2 + H""",
)

entry(
    index = 558,
    label = "OCHO + O2 <=> CO2 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is OCHO + O2 <=> CO2 + HO2""",
)

entry(
    index = 559,
    label = "CH3C(O)OOH <=> CH3C(O)O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+14, 's^-1'), n=0, Ea=(40142, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3C(O)OOH <=> CH3C(O)O + OH""",
)

entry(
    index = 560,
    label = "CH3C(O)OOH + H <=> CH3C(O)OO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.4e+10, 'cm^3/(mol*s)'), n=0, Ea=(1860, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3C(O)OOH + H <=> CH3C(O)OO + H2""",
)

entry(
    index = 561,
    label = "CH3C(O)OOH + O <=> CH3C(O)OO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.7e+12, 'cm^3/(mol*s)'), n=0, Ea=(4750, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3C(O)OOH + O <=> CH3C(O)OO + OH""",
)

entry(
    index = 562,
    label = "CH3C(O)OOH + OH <=> CH3C(O)OO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.1e+12, 'cm^3/(mol*s)'), n=0, Ea=(-437, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3C(O)OOH + OH <=> CH3C(O)OO + H2O""",
)

entry(
    index = 563,
    label = "CH3C(O)OOH + HO2 <=> CH3C(O)OO + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(41000, 'cm^3/(mol*s)'), n=2.5, Ea=(10206, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3C(O)OOH + HO2 <=> CH3C(O)OO + H2O2""",
)

entry(
    index = 564,
    label = "CH3C(O)OO + H <=> CH3C(O)O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3C(O)OO + H <=> CH3C(O)O + OH""",
)

entry(
    index = 565,
    label = "CH3C(O)OO + O <=> CH3 + CO2 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3C(O)OO + O <=> CH3 + CO2 + O2""",
)

entry(
    index = 566,
    label = "CH3C(O)OO + O <=> CH3O + CO + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3C(O)OO + O <=> CH3O + CO + O2""",
)

entry(
    index = 567,
    label = "CH3C(O)OO + OH <=> CH3C(O)O + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+11, 'cm^3/(mol*s)'), n=0.6, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3C(O)OO + OH <=> CH3C(O)O + HO2""",
)

entry(
    index = 568,
    label = "CH3C(O)OO + HO2 <=> CH3C(O)O + OH + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.9e+11, 'cm^3/(mol*s)'), n=0, Ea=(-1950, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3C(O)OO + HO2 <=> CH3C(O)O + OH + O2""",
)

entry(
    index = 569,
    label = "CH3C(O)OO + HO2 <=> CH3C(O)OOH + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+11, 'cm^3/(mol*s)'), n=0, Ea=(-1950, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3C(O)OO + HO2 <=> CH3C(O)OOH + O2""",
)

entry(
    index = 570,
    label = "CH3C(O)OO + CH3OO <=> CH3C(O)O + CH3O + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.1e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3C(O)OO + CH3OO <=> CH3C(O)O + CH3O + O2""",
)

entry(
    index = 571,
    label = "CH3CHO + CH3C(O)OO <=> CH3CO + CH3C(O)OOH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+13, 'cm^3/(mol*s)'), n=0, Ea=(16293, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHO + CH3C(O)OO <=> CH3CO + CH3C(O)OOH""",
)

entry(
    index = 572,
    label = "CH3CHO + CH3C(O)OO <=> CH2CHO + CH3C(O)OOH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.1e+13, 'cm^3/(mol*s)'), n=0, Ea=(23248, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHO + CH3C(O)OO <=> CH2CHO + CH3C(O)OOH""",
)

entry(
    index = 573,
    label = "CH3C(O)O <=> CH3 + CO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.2e+11, 's^-1'), n=0.29, Ea=(4579, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3C(O)O <=> CH3 + CO2""",
)

entry(
    index = 574,
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
    index = 575,
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
    index = 576,
    label = "HOCH2CH2OO + HO2 => CH2OOH + CH2OH + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+11, 'cm^3/(mol*s)'), n=0, Ea=(-1490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HOCH2CH2OO + HO2 => CH2OOH + CH2OH + O2""",
)

entry(
    index = 577,
    label = "HOCH2CH2OO + CH2O => CH2OOH + CH2OH + HCO",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(41000, 'cm^3/(mol*s)'), n=2.5, Ea=(10206, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HOCH2CH2OO + CH2O => CH2OOH + CH2OH + HCO""",
)

entry(
    index = 578,
    label = "HOCH2CH2OO + C2H4 => CH2O + CH2OH + CH3CHO",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.2e+12, 'cm^3/(mol*s)'), n=0, Ea=(17200, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HOCH2CH2OO + C2H4 => CH2O + CH2OH + CH3CHO""",
)

