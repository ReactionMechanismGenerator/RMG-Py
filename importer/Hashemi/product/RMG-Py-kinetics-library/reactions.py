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
    index = 63,
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
    index = 64,
    label = "CH4 + H <=> CH3 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4100, 'cm^3/(mol*s)'), n=3.156, Ea=(8755, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH4 + H <=> CH3 + H2""",
)

entry(
    index = 65,
    label = "CH4 + O <=> CH3 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(440000, 'cm^3/(mol*s)'), n=2.5, Ea=(6577, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH4 + O <=> CH3 + OH""",
)

entry(
    index = 66,
    label = "CH4 + OH <=> CH3 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+06, 'cm^3/(mol*s)'), n=2.182, Ea=(2506, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH4 + OH <=> CH3 + H2O""",
)

entry(
    index = 67,
    label = "CH4 + HO2 <=> CH3 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(47000, 'cm^3/(mol*s)'), n=2.5, Ea=(21000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH4 + HO2 <=> CH3 + H2O2""",
)

entry(
    index = 68,
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
    index = 71,
    label = "CH4 + CH <=> C2H4 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(-400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH4 + CH <=> C2H4 + H""",
)

entry(
    index = 72,
    label = "CH3 <=> CH + H2",
    degeneracy = 1,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(A=(3.1e+15, 'cm^3/(mol*s)'), n=0, Ea=(80871, 'cal/mol'), T0=(1, 'K')),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is CH3 <=> CH + H2""",
)

entry(
    index = 77,
    label = "CH3 + O <=> H2 + CO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + O <=> H2 + CO + H""",
)

entry(
    index = 79,
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
    index = 92,
    label = "CH3 + CH3 <=> C2H5 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(16055, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + CH3 <=> C2H5 + H""",
)

entry(
    index = 120,
    label = "CH + H <=> C + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH + H <=> C + H2""",
)

entry(
    index = 121,
    label = "CH + O <=> CO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.7e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH + O <=> CO + H""",
)

entry(
    index = 124,
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
    index = 125,
    label = "CH + OH <=> C + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+07, 'cm^3/(mol*s)'), n=2, Ea=(3000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH + OH <=> C + H2O""",
)

entry(
    index = 128,
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
    index = 132,
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
    index = 134,
    label = "C + OH <=> CO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C + OH <=> CO + H""",
)

entry(
    index = 135,
    label = "C + O2 <=> CO + O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C + O2 <=> CO + O""",
)

entry(
    index = 214,
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
    index = 215,
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
    index = 216,
    label = "C2H6 + O <=> C2H5 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(180000, 'cm^3/(mol*s)'), n=2.8, Ea=(5800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H6 + O <=> C2H5 + OH""",
)

entry(
    index = 217,
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
    index = 218,
    label = "C2H6 + HO2 <=> C2H5 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(26, 'cm^3/(mol*s)'), n=3.37, Ea=(15900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H6 + HO2 <=> C2H5 + H2O2""",
)

entry(
    index = 219,
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
    index = 220,
    label = "C2H6 + CH3 <=> C2H5 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(35, 'cm^3/(mol*s)'), n=3.44, Ea=(10384, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H6 + CH3 <=> C2H5 + CH4""",
)

entry(
    index = 223,
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
    index = 224,
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
    index = 227,
    label = "C2H5 + O <=> C2H4 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5 + O <=> C2H4 + OH""",
)

entry(
    index = 228,
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
    index = 234,
    label = "C2H5 + CH3 <=> C2H4 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5 + CH3 <=> C2H4 + CH4""",
)

entry(
    index = 237,
    label = "C2H5 + C2H5 <=> C2H6 + C2H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5 + C2H5 <=> C2H6 + C2H4""",
)

entry(
    index = 245,
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
    index = 254,
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
    index = 315,
    label = "C2H + O <=> CH + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H + O <=> CH + CO""",
)

entry(
    index = 318,
    label = "C2H + O2 <=> CO + CO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.7e+13, 'cm^3/(mol*s)'), n=-0.16, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H + O2 <=> CO + CO + H""",
)

