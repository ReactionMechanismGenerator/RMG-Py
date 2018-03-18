#!/usr/bin/env python
# encoding: utf-8

name = "/home/alongd/Code/RMG-Py/importer/JetSurF2.0"
shortDesc = u"/home/alongd/Code/RMG-Py/importer/JetSurF2.0/Mech_JetSurF2.0.txt"
longDesc = u"""
Unknown source
"""
entry(
    index = 1,
    label = "H + O2 <=> O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.644e+16, 'cm^3/(mol*s)'),
        n = -0.6707,
        Ea = (17041, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is H + O2 <=> O + OH""",
)

entry(
    index = 2,
    label = "O + H2 <=> H + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(45890, 'cm^3/(mol*s)'), n=2.7, Ea=(6260, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is O + H2 <=> H + OH""",
)

entry(
    index = 3,
    label = "OH + H2 <=> H + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.734e+08, 'cm^3/(mol*s)'),
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
    kinetics = Arrhenius(A=(39730, 'cm^3/(mol*s)'), n=2.4, Ea=(-2110, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is OH + OH <=> O + H2O""",
)

entry(
    index = 5,
    label = "H + H <=> H2",
    degeneracy = 1,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(A=(1.78e+18, 'cm^6/(mol^2*s)'), n=-1, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        efficiencies = {'[H][H]': 0, 'O=C=O': 0, 'O': 0, '[Ar]': 0.63},
    ),
    shortDesc = u"""The chemkin file reaction is H + H <=> H2""",
)

entry(
    index = 6,
    label = "H + H + H2 <=> H2 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+16, 'cm^6/(mol^2*s)'), n=-0.6, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H + H + H2 <=> H2 + H2""",
)

entry(
    index = 7,
    label = "H + H + H2O <=> H2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.624e+19, 'cm^6/(mol^2*s)'),
        n = -1.25,
        Ea = (0, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is H + H + H2O <=> H2 + H2O""",
)

entry(
    index = 8,
    label = "H + H + CO2 <=> H2 + CO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.5e+20, 'cm^6/(mol^2*s)'), n=-2, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H + H + CO2 <=> H2 + CO2""",
)

entry(
    index = 9,
    label = "H + OH <=> H2O",
    degeneracy = 1,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(A=(4.4e+22, 'cm^6/(mol^2*s)'), n=-2, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        efficiencies = {'[C-]#[O+]': 1.75, '[H][H]': 2, 'O=C=O': 3.6, 'O': 6.3, '[Ar]': 0.38},
    ),
    shortDesc = u"""The chemkin file reaction is H + OH <=> H2O""",
)

entry(
    index = 10,
    label = "O + H <=> OH",
    degeneracy = 1,
    duplicate = True,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(A=(9.428e+18, 'cm^6/(mol^2*s)'), n=-1, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        efficiencies = {'[C-]#[O+]': 1.75, '[H][H]': 2, 'O=C=O': 3.6, 'O': 12, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is O + H <=> OH""",
)

entry(
    index = 11,
    label = "O + O <=> O2",
    degeneracy = 1,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(A=(1.2e+17, 'cm^6/(mol^2*s)'), n=-1, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        efficiencies = {'[C-]#[O+]': 1.75, '[H][H]': 2.4, 'O=C=O': 3.6, 'O': 15.4, '[Ar]': 0.83},
    ),
    shortDesc = u"""The chemkin file reaction is O + O <=> O2""",
)

entry(
    index = 12,
    label = "H + O2 <=> HO2",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(5.116e+12, 'cm^3/(mol*s)'), n=0.44, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (6.328e+19, 'cm^6/(mol^2*s)'),
            n = -1.4,
            Ea = (0, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.5,
        T3 = (1e-30, 'K'),
        T1 = (1e+30, 'K'),
        efficiencies = {'[C-]#[O+]': 1.09, '[O][O]': 0.85, 'O=C=O': 2.18, 'O': 11.89, '[Ar]': 0.4},
    ),
    shortDesc = u"""The chemkin file reaction is H + O2 <=> HO2""",
)

entry(
    index = 13,
    label = "H2 + O2 <=> HO2 + H",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (591600, 'cm^3/(mol*s)'),
        n = 2.433,
        Ea = (53502, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is H2 + O2 <=> HO2 + H""",
)

entry(
    index = 14,
    label = "OH + OH <=> H2O2",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.11e+14, 'cm^3/(mol*s)'), n=-0.37, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (2.01e+17, 'cm^6/(mol^2*s)'),
            n = -0.584,
            Ea = (-2293, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.7346,
        T3 = (94, 'K'),
        T1 = (1756, 'K'),
        T2 = (5182, 'K'),
        efficiencies = {'[C-]#[O+]': 1.75, '[H][H]': 2, 'O=C=O': 3.6, 'O': 6, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is OH + OH <=> H2O2""",
)

entry(
    index = 15,
    label = "HO2 + H <=> O + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.97e+12, 'cm^3/(mol*s)'), n=0, Ea=(671, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HO2 + H <=> O + H2O""",
)

entry(
    index = 16,
    label = "HO2 + H <=> OH + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.485e+13, 'cm^3/(mol*s)'), n=0, Ea=(295, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HO2 + H <=> OH + OH""",
)

entry(
    index = 17,
    label = "HO2 + O <=> OH + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HO2 + O <=> OH + O2""",
)

entry(
    index = 18,
    label = "HO2 + HO2 <=> O2 + H2O2",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(1.3e+11, 'cm^3/(mol*s)'), n=0, Ea=(-1630, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (3.658e+14, 'cm^3/(mol*s)'),
                n = 0,
                Ea = (12000, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is HO2 + HO2 <=> O2 + H2O2""",
)

entry(
    index = 19,
    label = "HO2 + OH <=> H2O + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.89e+13, 'cm^3/(mol*s)'), n=0, Ea=(-500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HO2 + OH <=> H2O + O2""",
)

entry(
    index = 20,
    label = "H2O2 + H <=> HO2 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.05e+06, 'cm^3/(mol*s)'), n=2, Ea=(5200, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2O2 + H <=> HO2 + H2""",
)

entry(
    index = 21,
    label = "H2O2 + H <=> OH + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.41e+13, 'cm^3/(mol*s)'), n=0, Ea=(3970, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2O2 + H <=> OH + H2O""",
)

entry(
    index = 22,
    label = "H2O2 + O <=> OH + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.63e+06, 'cm^3/(mol*s)'), n=2, Ea=(3970, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2O2 + O <=> OH + HO2""",
)

entry(
    index = 23,
    label = "H2O2 + OH <=> HO2 + H2O",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(2e+12, 'cm^3/(mol*s)'), n=0, Ea=(427, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (2.67e+41, 'cm^3/(mol*s)'),
                n = -7,
                Ea = (37600, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is H2O2 + OH <=> HO2 + H2O""",
)

entry(
    index = 24,
    label = "CO + O <=> CO2",
    degeneracy = 1,
    kinetics = Lindemann(
        arrheniusHigh = Arrhenius(A=(1.362e+10, 'cm^3/(mol*s)'), n=0, Ea=(2384, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.173e+24, 'cm^6/(mol^2*s)'),
            n = -2.79,
            Ea = (4191, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        efficiencies = {'[C-]#[O+]': 1.75, '[H][H]': 2, 'O=C=O': 3.6, 'O': 12, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is CO + O <=> CO2""",
)

entry(
    index = 25,
    label = "CO + OH <=> CO2 + H",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(
                A = (70460, 'cm^3/(mol*s)'),
                n = 2.053,
                Ea = (-355.67, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (5.757e+12, 'cm^3/(mol*s)'),
                n = -0.664,
                Ea = (331.83, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CO + OH <=> CO2 + H""",
)

entry(
    index = 26,
    label = "CO + O2 <=> CO2 + O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.119e+12, 'cm^3/(mol*s)'),
        n = 0,
        Ea = (47700, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CO + O2 <=> CO2 + O""",
)

entry(
    index = 27,
    label = "CO + HO2 <=> CO2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (157000, 'cm^3/(mol*s)'),
        n = 2.18,
        Ea = (17942.6, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CO + HO2 <=> CO2 + OH""",
)

entry(
    index = 28,
    label = "HCO + H <=> CO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCO + H <=> CO + H2""",
)

entry(
    index = 29,
    label = "HCO + O <=> CO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCO + O <=> CO + OH""",
)

entry(
    index = 30,
    label = "HCO + O <=> CO2 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCO + O <=> CO2 + H""",
)

entry(
    index = 31,
    label = "HCO + OH <=> CO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.02e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCO + OH <=> CO + H2O""",
)

entry(
    index = 32,
    label = "HCO <=> CO + H",
    degeneracy = 1,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(
            A = (1.87e+17, 'cm^3/(mol*s)'),
            n = -1,
            Ea = (17000, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        efficiencies = {'[C-]#[O+]': 1.75, '[H][H]': 2, 'O=C=O': 3.6, 'O': 0},
    ),
    shortDesc = u"""The chemkin file reaction is HCO <=> CO + H""",
)

entry(
    index = 33,
    label = "HCO + H2O <=> CO + H + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.244e+18, 'cm^3/(mol*s)'),
        n = -1,
        Ea = (17000, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is HCO + H2O <=> CO + H + H2O""",
)

entry(
    index = 34,
    label = "HCO + O2 <=> CO + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.204e+10, 'cm^3/(mol*s)'),
        n = 0.807,
        Ea = (-727, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is HCO + O2 <=> CO + HO2""",
)

entry(
    index = 35,
    label = "CO + H2 <=> CH2O",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(
            A = (4.3e+07, 'cm^3/(mol*s)'),
            n = 1.5,
            Ea = (79600, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        arrheniusLow = Arrhenius(
            A = (5.07e+27, 'cm^6/(mol^2*s)'),
            n = -3.42,
            Ea = (84350, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.932,
        T3 = (197, 'K'),
        T1 = (1540, 'K'),
        T2 = (10300, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is CO + H2 <=> CH2O""",
)

entry(
    index = 36,
    label = "C + OH <=> CO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C + OH <=> CO + H""",
)

entry(
    index = 37,
    label = "C + O2 <=> CO + O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(576, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C + O2 <=> CO + O""",
)

entry(
    index = 38,
    label = "CH + H <=> C + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH + H <=> C + H2""",
)

entry(
    index = 39,
    label = "CH + O <=> CO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.7e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH + O <=> CO + H""",
)

entry(
    index = 40,
    label = "CH + OH <=> HCO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH + OH <=> HCO + H""",
)

entry(
    index = 41,
    label = "CH + H2 <=> CH2 + H",
    degeneracy = 1,
    duplicate = True,
    kinetics = Arrhenius(
        A = (1.107e+08, 'cm^3/(mol*s)'),
        n = 1.79,
        Ea = (1670, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH + H2 <=> CH2 + H""",
)

entry(
    index = 42,
    label = "CH + H2O <=> CH2O + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.71e+12, 'cm^3/(mol*s)'), n=0, Ea=(-755, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH + H2O <=> CH2O + H""",
)

entry(
    index = 43,
    label = "CH + O2 <=> HCO + O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH + O2 <=> HCO + O""",
)

entry(
    index = 44,
    label = "CH + CO <=> HCCO",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (2.69e+28, 'cm^6/(mol^2*s)'),
            n = -3.74,
            Ea = (1936, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.5757,
        T3 = (237, 'K'),
        T1 = (1652, 'K'),
        T2 = (5069, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is CH + CO <=> HCCO""",
)

entry(
    index = 45,
    label = "CH + CO2 <=> HCO + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.4e+12, 'cm^3/(mol*s)'), n=0, Ea=(690, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH + CO2 <=> HCO + CO""",
)

entry(
    index = 46,
    label = "HCO + H <=> CH2O",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(
            A = (1.09e+12, 'cm^3/(mol*s)'),
            n = 0.48,
            Ea = (-260, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        arrheniusLow = Arrhenius(
            A = (1.35e+24, 'cm^6/(mol^2*s)'),
            n = -2.57,
            Ea = (1425, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.7824,
        T3 = (271, 'K'),
        T1 = (2755, 'K'),
        T2 = (6570, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is HCO + H <=> CH2O""",
)

entry(
    index = 47,
    label = "CH2 + H <=> CH3",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2.5e+16, 'cm^3/(mol*s)'), n=-0.8, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.2e+27, 'cm^6/(mol^2*s)'),
            n = -3.14,
            Ea = (1230, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.68,
        T3 = (78, 'K'),
        T1 = (1995, 'K'),
        T2 = (5590, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is CH2 + H <=> CH3""",
)

entry(
    index = 48,
    label = "CH2 + O <=> HCO + H",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(8e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2 + O <=> HCO + H""",
)

entry(
    index = 49,
    label = "CH2 + OH <=> CH2O + H",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2 + OH <=> CH2O + H""",
)

entry(
    index = 50,
    label = "CH2 + OH <=> CH + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.13e+07, 'cm^3/(mol*s)'), n=2, Ea=(3000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2 + OH <=> CH + H2O""",
)

entry(
    index = 51,
    label = "CH2 + H2 <=> H + CH3",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(500000, 'cm^3/(mol*s)'), n=2, Ea=(7230, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(7e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2 + H2 <=> H + CH3""",
)

entry(
    index = 52,
    label = "CH2 + O2 <=> HCO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.06e+13, 'cm^3/(mol*s)'), n=0, Ea=(1500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2 + O2 <=> HCO + OH""",
)

entry(
    index = 53,
    label = "CH2 + O2 <=> CO2 + H + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.64e+12, 'cm^3/(mol*s)'), n=0, Ea=(1500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2 + O2 <=> CO2 + H + H""",
)

entry(
    index = 54,
    label = "CH2 + HO2 <=> CH2O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2 + HO2 <=> CH2O + OH""",
)

entry(
    index = 55,
    label = "CH2 + C <=> C2H + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2 + C <=> C2H + H""",
)

entry(
    index = 56,
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
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is CH2 + CO <=> CH2CO""",
)

entry(
    index = 57,
    label = "CH2 + CH <=> C2H2 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2 + CH <=> C2H2 + H""",
)

entry(
    index = 58,
    label = "CH2 + CH2 <=> C2H2 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2 + CH2 <=> C2H2 + H2""",
)

entry(
    index = 59,
    label = "CH2* + N2 <=> CH2 + N2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.5e+13, 'cm^3/(mol*s)'), n=0, Ea=(600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2* + N2 <=> CH2 + N2""",
)

entry(
    index = 60,
    label = "CH2* + AR <=> CH2 + AR",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+12, 'cm^3/(mol*s)'), n=0, Ea=(600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2* + AR <=> CH2 + AR""",
)

entry(
    index = 61,
    label = "CH2* + H <=> CH + H2",
    degeneracy = 1,
    duplicate = True,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2* + H <=> CH + H2""",
)

entry(
    index = 62,
    label = "CH2* + O <=> CO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2* + O <=> CO + H2""",
)

entry(
    index = 63,
    label = "CH2* + O2 <=> H + OH + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2* + O2 <=> H + OH + CO""",
)

entry(
    index = 64,
    label = "CH2* + O2 <=> CO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2* + O2 <=> CO + H2O""",
)

entry(
    index = 65,
    label = "CH2* + H2O <=> CH3OH",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (2.7e+38, 'cm^6/(mol^2*s)'),
            n = -6.3,
            Ea = (3100, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.1507,
        T3 = (134, 'K'),
        T1 = (2383, 'K'),
        T2 = (7265, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5},
    ),
    shortDesc = u"""The chemkin file reaction is CH2* + H2O <=> CH3OH""",
)

entry(
    index = 66,
    label = "CH2* + H2O <=> CH2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2* + H2O <=> CH2 + H2O""",
)

entry(
    index = 67,
    label = "CH2* + CO <=> CH2 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2* + CO <=> CH2 + CO""",
)

entry(
    index = 68,
    label = "CH2* + CO2 <=> CH2 + CO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2* + CO2 <=> CH2 + CO2""",
)

entry(
    index = 69,
    label = "CH2* + CO2 <=> CH2O + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2* + CO2 <=> CH2O + CO""",
)

entry(
    index = 70,
    label = "CH2O + H <=> CH2OH",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(
            A = (5.4e+11, 'cm^3/(mol*s)'),
            n = 0.454,
            Ea = (3600, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        arrheniusLow = Arrhenius(
            A = (1.27e+32, 'cm^6/(mol^2*s)'),
            n = -4.82,
            Ea = (6530, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.7187,
        T3 = (103, 'K'),
        T1 = (1291, 'K'),
        T2 = (4160, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5},
    ),
    shortDesc = u"""The chemkin file reaction is CH2O + H <=> CH2OH""",
)

entry(
    index = 71,
    label = "CH2O + H <=> CH3O",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(
            A = (5.4e+11, 'cm^3/(mol*s)'),
            n = 0.454,
            Ea = (2600, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        arrheniusLow = Arrhenius(
            A = (2.2e+30, 'cm^6/(mol^2*s)'),
            n = -4.8,
            Ea = (5560, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.758,
        T3 = (94, 'K'),
        T1 = (1555, 'K'),
        T2 = (4200, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5},
    ),
    shortDesc = u"""The chemkin file reaction is CH2O + H <=> CH3O""",
)

entry(
    index = 72,
    label = "CH2O + H <=> HCO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.3e+10, 'cm^3/(mol*s)'),
        n = 1.05,
        Ea = (3275, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2O + H <=> HCO + H2""",
)

entry(
    index = 73,
    label = "CH2O + O <=> HCO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.9e+13, 'cm^3/(mol*s)'), n=0, Ea=(3540, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2O + O <=> HCO + OH""",
)

entry(
    index = 74,
    label = "CH2O + OH <=> HCO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.43e+09, 'cm^3/(mol*s)'),
        n = 1.18,
        Ea = (-447, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2O + OH <=> HCO + H2O""",
)

entry(
    index = 75,
    label = "CH2O + O2 <=> HCO + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(40000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2O + O2 <=> HCO + HO2""",
)

entry(
    index = 76,
    label = "CH2O + HO2 <=> HCO + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(8000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2O + HO2 <=> HCO + H2O2""",
)

entry(
    index = 77,
    label = "CH2O + CH <=> CH2CO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.46e+13, 'cm^3/(mol*s)'), n=0, Ea=(-515, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2O + CH <=> CH2CO + H""",
)

entry(
    index = 78,
    label = "CH3 + H <=> CH4",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(
            A = (1.27e+16, 'cm^3/(mol*s)'),
            n = -0.63,
            Ea = (383, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        arrheniusLow = Arrhenius(
            A = (2.477e+33, 'cm^6/(mol^2*s)'),
            n = -4.76,
            Ea = (2440, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.783,
        T3 = (74, 'K'),
        T1 = (2941, 'K'),
        T2 = (6964, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is CH3 + H <=> CH4""",
)

entry(
    index = 79,
    label = "CH3 + O <=> CH2O + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.43e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + O <=> CH2O + H""",
)

entry(
    index = 80,
    label = "CH3 + OH <=> CH3OH",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(6.3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (2.7e+38, 'cm^6/(mol^2*s)'),
            n = -6.3,
            Ea = (3100, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.2105,
        T3 = (83.5, 'K'),
        T1 = (5398, 'K'),
        T2 = (8370, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5},
    ),
    shortDesc = u"""The chemkin file reaction is CH3 + OH <=> CH3OH""",
)

entry(
    index = 81,
    label = "CH3 + OH <=> CH2 + H2O",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(5.6e+07, 'cm^3/(mol*s)'), n=1.6, Ea=(5420, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2.501e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH3 + OH <=> CH2 + H2O""",
)

entry(
    index = 82,
    label = "CH3 + O2 <=> O + CH3O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.083e+13, 'cm^3/(mol*s)'),
        n = 0,
        Ea = (28800, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3 + O2 <=> O + CH3O""",
)

entry(
    index = 83,
    label = "CH3 + O2 <=> OH + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.6e+10, 'cm^3/(mol*s)'), n=0, Ea=(8940, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + O2 <=> OH + CH2O""",
)

entry(
    index = 84,
    label = "CH3 + HO2 <=> CH4 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + HO2 <=> CH4 + O2""",
)

entry(
    index = 85,
    label = "CH3 + HO2 <=> CH3O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.34e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + HO2 <=> CH3O + OH""",
)

entry(
    index = 86,
    label = "CH3 + H2O2 <=> CH4 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(24500, 'cm^3/(mol*s)'), n=2.47, Ea=(5180, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + H2O2 <=> CH4 + HO2""",
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
    label = "CH3 + CH <=> C2H3 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + CH <=> C2H3 + H""",
)

entry(
    index = 89,
    label = "CH3 + HCO <=> CH4 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.48e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + HCO <=> CH4 + CO""",
)

entry(
    index = 90,
    label = "CH3 + CH2O <=> CH4 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3320, 'cm^3/(mol*s)'), n=2.81, Ea=(5860, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + CH2O <=> CH4 + HCO""",
)

entry(
    index = 91,
    label = "CH3 + CH2 <=> C2H4 + H",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.2e+13, 'cm^3/(mol*s)'), n=0, Ea=(-570, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH3 + CH2 <=> C2H4 + H""",
)

entry(
    index = 92,
    label = "CH3 + CH3 <=> C2H6",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(
            A = (2.12e+16, 'cm^3/(mol*s)'),
            n = -0.97,
            Ea = (620, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        arrheniusLow = Arrhenius(
            A = (1.77e+50, 'cm^6/(mol^2*s)'),
            n = -9.67,
            Ea = (6220, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.5325,
        T3 = (151, 'K'),
        T1 = (1038, 'K'),
        T2 = (4970, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is CH3 + CH3 <=> C2H6""",
)

entry(
    index = 93,
    label = "CH3 + CH3 <=> H + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.99e+12, 'cm^3/(mol*s)'),
        n = 0.1,
        Ea = (10600, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3 + CH3 <=> H + C2H5""",
)

entry(
    index = 94,
    label = "CH3 + HCCO <=> C2H4 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + HCCO <=> C2H4 + CO""",
)

entry(
    index = 95,
    label = "CH3 + C2H <=> C3H3 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.41e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + C2H <=> C3H3 + H""",
)

entry(
    index = 96,
    label = "CH3O + H <=> CH3OH",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (8.6e+28, 'cm^6/(mol^2*s)'),
            n = -4,
            Ea = (3025, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.8902,
        T3 = (144, 'K'),
        T1 = (2838, 'K'),
        T2 = (45569, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5},
    ),
    shortDesc = u"""The chemkin file reaction is CH3O + H <=> CH3OH""",
)

entry(
    index = 97,
    label = "CH3O + H <=> CH2OH + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.4e+06, 'cm^3/(mol*s)'), n=1.6, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O + H <=> CH2OH + H""",
)

entry(
    index = 98,
    label = "CH3O + H <=> CH2O + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O + H <=> CH2O + H2""",
)

entry(
    index = 99,
    label = "CH3O + H <=> CH3 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O + H <=> CH3 + OH""",
)

entry(
    index = 100,
    label = "CH3O + H <=> CH2* + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O + H <=> CH2* + H2O""",
)

entry(
    index = 101,
    label = "CH3O + O <=> CH2O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O + O <=> CH2O + OH""",
)

entry(
    index = 102,
    label = "CH3O + OH <=> CH2O + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O + OH <=> CH2O + H2O""",
)

entry(
    index = 103,
    label = "CH3O + O2 <=> CH2O + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.28e-13, 'cm^3/(mol*s)'),
        n = 7.6,
        Ea = (-3530, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3O + O2 <=> CH2O + HO2""",
)

entry(
    index = 104,
    label = "CH2OH + H <=> CH3OH",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3e+31, 'cm^6/(mol^2*s)'),
            n = -4.8,
            Ea = (3300, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.7679,
        T3 = (338, 'K'),
        T1 = (1812, 'K'),
        T2 = (5081, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5},
    ),
    shortDesc = u"""The chemkin file reaction is CH2OH + H <=> CH3OH""",
)

entry(
    index = 105,
    label = "CH2OH + H <=> CH2O + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2OH + H <=> CH2O + H2""",
)

entry(
    index = 106,
    label = "CH2OH + H <=> CH3 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2OH + H <=> CH3 + OH""",
)

entry(
    index = 107,
    label = "CH2OH + H <=> CH2* + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2OH + H <=> CH2* + H2O""",
)

entry(
    index = 108,
    label = "CH2OH + O <=> CH2O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2OH + O <=> CH2O + OH""",
)

entry(
    index = 109,
    label = "CH2OH + OH <=> CH2O + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2OH + OH <=> CH2O + H2O""",
)

entry(
    index = 110,
    label = "CH2OH + O2 <=> CH2O + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2OH + O2 <=> CH2O + HO2""",
)

entry(
    index = 111,
    label = "CH4 + H <=> CH3 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (6.6e+08, 'cm^3/(mol*s)'),
        n = 1.62,
        Ea = (10840, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH4 + H <=> CH3 + H2""",
)

entry(
    index = 112,
    label = "CH4 + O <=> CH3 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.02e+09, 'cm^3/(mol*s)'),
        n = 1.5,
        Ea = (8600, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH4 + O <=> CH3 + OH""",
)

entry(
    index = 113,
    label = "CH4 + OH <=> CH3 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+08, 'cm^3/(mol*s)'), n=1.6, Ea=(3120, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH4 + OH <=> CH3 + H2O""",
)

entry(
    index = 114,
    label = "CH4 + CH <=> C2H4 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH4 + CH <=> C2H4 + H""",
)

entry(
    index = 115,
    label = "CH4 + CH2 <=> CH3 + CH3",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(2.46e+06, 'cm^3/(mol*s)'), n=2, Ea=(8270, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(-570, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH4 + CH2 <=> CH3 + CH3""",
)

entry(
    index = 116,
    label = "CH4 + C2H <=> C2H2 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.81e+12, 'cm^3/(mol*s)'), n=0, Ea=(500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH4 + C2H <=> C2H2 + CH3""",
)

entry(
    index = 117,
    label = "CH3OH + H <=> CH2OH + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+07, 'cm^3/(mol*s)'), n=2.1, Ea=(4870, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OH + H <=> CH2OH + H2""",
)

entry(
    index = 118,
    label = "CH3OH + H <=> CH3O + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.2e+06, 'cm^3/(mol*s)'), n=2.1, Ea=(4870, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OH + H <=> CH3O + H2""",
)

entry(
    index = 119,
    label = "CH3OH + O <=> CH2OH + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(388000, 'cm^3/(mol*s)'), n=2.5, Ea=(3100, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OH + O <=> CH2OH + OH""",
)

entry(
    index = 120,
    label = "CH3OH + O <=> CH3O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(130000, 'cm^3/(mol*s)'), n=2.5, Ea=(5000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OH + O <=> CH3O + OH""",
)

entry(
    index = 121,
    label = "CH3OH + OH <=> CH2OH + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.44e+06, 'cm^3/(mol*s)'), n=2, Ea=(-840, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OH + OH <=> CH2OH + H2O""",
)

entry(
    index = 122,
    label = "CH3OH + OH <=> CH3O + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.3e+06, 'cm^3/(mol*s)'), n=2, Ea=(1500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OH + OH <=> CH3O + H2O""",
)

entry(
    index = 123,
    label = "CH3OH + CH3 <=> CH2OH + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+07, 'cm^3/(mol*s)'), n=1.5, Ea=(9940, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OH + CH3 <=> CH2OH + CH4""",
)

entry(
    index = 124,
    label = "CH3OH + CH3 <=> CH3O + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+07, 'cm^3/(mol*s)'), n=1.5, Ea=(9940, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OH + CH3 <=> CH3O + CH4""",
)

entry(
    index = 125,
    label = "C2H + H <=> C2H2",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1e+17, 'cm^3/(mol*s)'), n=-1, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.75e+33, 'cm^6/(mol^2*s)'),
            n = -4.8,
            Ea = (1900, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.6464,
        T3 = (132, 'K'),
        T1 = (1315, 'K'),
        T2 = (5566, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C2H + H <=> C2H2""",
)

entry(
    index = 126,
    label = "C2H + O <=> CH + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H + O <=> CH + CO""",
)

entry(
    index = 127,
    label = "C2H + OH <=> H + HCCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H + OH <=> H + HCCO""",
)

entry(
    index = 128,
    label = "C2H + O2 <=> HCO + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(1500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H + O2 <=> HCO + CO""",
)

entry(
    index = 129,
    label = "C2H + H2 <=> H + C2H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(490000, 'cm^3/(mol*s)'), n=2.5, Ea=(560, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H + H2 <=> H + C2H2""",
)

entry(
    index = 130,
    label = "C2O + H <=> CH + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2O + H <=> CH + CO""",
)

entry(
    index = 131,
    label = "C2O + O <=> CO + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2O + O <=> CO + CO""",
)

entry(
    index = 132,
    label = "C2O + OH <=> CO + CO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2O + OH <=> CO + CO + H""",
)

entry(
    index = 133,
    label = "C2O + O2 <=> CO + CO + O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2O + O2 <=> CO + CO + O""",
)

entry(
    index = 134,
    label = "HCCO + H <=> CH2* + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCCO + H <=> CH2* + CO""",
)

entry(
    index = 135,
    label = "HCCO + O <=> H + CO + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCCO + O <=> H + CO + CO""",
)

entry(
    index = 136,
    label = "HCCO + O2 <=> OH + CO + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.6e+12, 'cm^3/(mol*s)'), n=0, Ea=(854, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCCO + O2 <=> OH + CO + CO""",
)

entry(
    index = 137,
    label = "HCCO + CH <=> C2H2 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCCO + CH <=> C2H2 + CO""",
)

entry(
    index = 138,
    label = "HCCO + CH2 <=> C2H3 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCCO + CH2 <=> C2H3 + CO""",
)

entry(
    index = 139,
    label = "HCCO + HCCO <=> C2H2 + CO + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCCO + HCCO <=> C2H2 + CO + CO""",
)

entry(
    index = 140,
    label = "HCCO + OH <=> C2O + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCCO + OH <=> C2O + H2O""",
)

entry(
    index = 141,
    label = "C2H2 <=> H2CC",
    degeneracy = 1,
    kinetics = Lindemann(
        arrheniusHigh = Arrhenius(A=(8e+14, 's^-1'), n=-0.52, Ea=(50750, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (2.45e+15, 'cm^3/(mol*s)'),
            n = -0.64,
            Ea = (49700, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, 'C=C': 2.5, 'C#C': 2.5},
    ),
    shortDesc = u"""The chemkin file reaction is C2H2 <=> H2CC""",
)

entry(
    index = 142,
    label = "C2H3 <=> C2H2 + H",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.86e+08, 's^-1'), n=1.62, Ea=(37048.2, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (2.565e+27, 'cm^3/(mol*s)'),
            n = -3.4,
            Ea = (35798.7, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 1.9816,
        T3 = (5383.7, 'K'),
        T1 = (4.2932, 'K'),
        T2 = (-0.0795, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, 'C=C': 3, 'C#C': 3, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C2H3 <=> C2H2 + H""",
)

entry(
    index = 143,
    label = "C2H2 + O <=> C2H + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.6e+19, 'cm^3/(mol*s)'),
        n = -1.41,
        Ea = (28950, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H2 + O <=> C2H + OH""",
)

entry(
    index = 144,
    label = "C2H2 + O <=> CH2 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.08e+06, 'cm^3/(mol*s)'), n=2, Ea=(1900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H2 + O <=> CH2 + CO""",
)

entry(
    index = 145,
    label = "C2H2 + O <=> HCCO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.632e+07, 'cm^3/(mol*s)'), n=2, Ea=(1900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H2 + O <=> HCCO + H""",
)

entry(
    index = 146,
    label = "C2H2 + OH <=> CH2CO + H",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (0.000218, 'cm^3/(mol*s)'),
        n = 4.5,
        Ea = (-1000, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H2 + OH <=> CH2CO + H""",
)

entry(
    index = 147,
    label = "C2H2 + OH <=> HCCOH + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(504000, 'cm^3/(mol*s)'), n=2.3, Ea=(13500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H2 + OH <=> HCCOH + H""",
)

entry(
    index = 148,
    label = "C2H2 + OH <=> C2H + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.37e+07, 'cm^3/(mol*s)'), n=2, Ea=(14000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H2 + OH <=> C2H + H2O""",
)

entry(
    index = 149,
    label = "C2H2 + OH <=> CH3 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.000483, 'cm^3/(mol*s)'), n=4, Ea=(-2000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H2 + OH <=> CH3 + CO""",
)

entry(
    index = 150,
    label = "C2H2 + HCO <=> C2H3 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+07, 'cm^3/(mol*s)'), n=2, Ea=(6000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H2 + HCO <=> C2H3 + CO""",
)

entry(
    index = 151,
    label = "C2H2 + CH2 <=> C3H3 + H",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(1.2e+13, 'cm^3/(mol*s)'), n=0, Ea=(6620, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H2 + CH2 <=> C3H3 + H""",
)

entry(
    index = 152,
    label = "C2H2 + C2H <=> C4H2 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H2 + C2H <=> C4H2 + H""",
)

entry(
    index = 153,
    label = "C2H2 + C2H <=> nC4H3",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(
            A = (8.3e+10, 'cm^3/(mol*s)'),
            n = 0.899,
            Ea = (-363, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        arrheniusLow = Arrhenius(
            A = (1.24e+31, 'cm^6/(mol^2*s)'),
            n = -4.718,
            Ea = (1871, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 1,
        T3 = (100, 'K'),
        T1 = (5613, 'K'),
        T2 = (13387, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, 'C=C': 2.5, 'C#C': 2.5},
    ),
    shortDesc = u"""The chemkin file reaction is C2H2 + C2H <=> nC4H3""",
)

entry(
    index = 154,
    label = "C2H2 + C2H <=> iC4H3",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(
            A = (8.3e+10, 'cm^3/(mol*s)'),
            n = 0.899,
            Ea = (-363, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        arrheniusLow = Arrhenius(
            A = (1.24e+31, 'cm^6/(mol^2*s)'),
            n = -4.718,
            Ea = (1871, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 1,
        T3 = (100, 'K'),
        T1 = (5613, 'K'),
        T2 = (13387, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, 'C=C': 2.5, 'C#C': 2.5},
    ),
    shortDesc = u"""The chemkin file reaction is C2H2 + C2H <=> iC4H3""",
)

entry(
    index = 155,
    label = "C2H2 + HCCO <=> C3H3 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(3000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H2 + HCCO <=> C3H3 + CO""",
)

entry(
    index = 156,
    label = "C2H2 + CH3 <=> pC3H4 + H",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.56e+09, 'cm^3/(mol*s)'),
        n = 1.1,
        Ea = (13644, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H2 + CH3 <=> pC3H4 + H""",
)

entry(
    index = 157,
    label = "C2H2 + CH3 <=> aC3H4 + H",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.14e+09, 'cm^3/(mol*s)'),
        n = 0.86,
        Ea = (22153, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H2 + CH3 <=> aC3H4 + H""",
)

entry(
    index = 158,
    label = "C2H2 + CH3 <=> CH3CCH2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.99e+22, 'cm^3/(mol*s)'),
        n = -4.39,
        Ea = (18850, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H2 + CH3 <=> CH3CCH2""",
)

entry(
    index = 159,
    label = "C2H2 + CH3 <=> CH3CHCH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.2e+35, 'cm^3/(mol*s)'),
        n = -7.76,
        Ea = (13300, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H2 + CH3 <=> CH3CHCH""",
)

entry(
    index = 160,
    label = "C2H2 + CH3 <=> aC3H5",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.68e+53, 'cm^3/(mol*s)'),
        n = -12.82,
        Ea = (35730, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H2 + CH3 <=> aC3H5""",
)

entry(
    index = 161,
    label = "H2CC + H <=> C2H2 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2CC + H <=> C2H2 + H""",
)

entry(
    index = 162,
    label = "H2CC + OH <=> CH2CO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2CC + OH <=> CH2CO + H""",
)

entry(
    index = 163,
    label = "H2CC + O2 <=> HCO + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2CC + O2 <=> HCO + HCO""",
)

entry(
    index = 164,
    label = "H2CC + C2H2 <=> C4H4",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(
            A = (350000, 'cm^3/(mol*s)'),
            n = 2.055,
            Ea = (-2400, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        arrheniusLow = Arrhenius(
            A = (1.4e+60, 'cm^6/(mol^2*s)'),
            n = -12.599,
            Ea = (7417, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.98,
        T3 = (56, 'K'),
        T1 = (580, 'K'),
        T2 = (4164, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, 'C=C': 3, 'C#C': 3},
    ),
    shortDesc = u"""The chemkin file reaction is H2CC + C2H2 <=> C4H4""",
)

entry(
    index = 165,
    label = "H2CC + C2H4 <=> C4H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2CC + C2H4 <=> C4H6""",
)

entry(
    index = 166,
    label = "CH2CO + H <=> CH2CHO",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(
            A = (3.3e+14, 'cm^3/(mol*s)'),
            n = -0.06,
            Ea = (8500, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        arrheniusLow = Arrhenius(
            A = (3.8e+41, 'cm^6/(mol^2*s)'),
            n = -7.64,
            Ea = (11900, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.337,
        T3 = (1707, 'K'),
        T1 = (3200, 'K'),
        T2 = (4131, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, 'C=C': 3, 'C#C': 3, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is CH2CO + H <=> CH2CHO""",
)

entry(
    index = 167,
    label = "CH2CO + H <=> HCCO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(8000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CO + H <=> HCCO + H2""",
)

entry(
    index = 168,
    label = "CH2CO + H <=> CH3 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.5e+09, 'cm^3/(mol*s)'),
        n = 1.43,
        Ea = (2690, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2CO + H <=> CH3 + CO""",
)

entry(
    index = 169,
    label = "CH2CO + O <=> HCCO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(8000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CO + O <=> HCCO + OH""",
)

entry(
    index = 170,
    label = "CH2CO + O <=> CH2 + CO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.75e+12, 'cm^3/(mol*s)'), n=0, Ea=(1350, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CO + O <=> CH2 + CO2""",
)

entry(
    index = 171,
    label = "CH2CO + OH <=> HCCO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(2000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CO + OH <=> HCCO + H2O""",
)

entry(
    index = 172,
    label = "HCCOH + H <=> CH2CO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCCOH + H <=> CH2CO + H""",
)

entry(
    index = 173,
    label = "C2H3 + H <=> C2H4",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(
            A = (6.08e+12, 'cm^3/(mol*s)'),
            n = 0.27,
            Ea = (280, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        arrheniusLow = Arrhenius(
            A = (1.4e+30, 'cm^6/(mol^2*s)'),
            n = -3.86,
            Ea = (3320, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.782,
        T3 = (207.5, 'K'),
        T1 = (2663, 'K'),
        T2 = (6095, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, 'C=C': 3, 'C#C': 3, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C2H3 + H <=> C2H4""",
)

entry(
    index = 174,
    label = "C2H3 + H <=> C2H2 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3 + H <=> C2H2 + H2""",
)

entry(
    index = 175,
    label = "C2H3 + H <=> H2CC + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3 + H <=> H2CC + H2""",
)

entry(
    index = 176,
    label = "C2H3 + O <=> CH2CO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3 + O <=> CH2CO + H""",
)

entry(
    index = 177,
    label = "C2H3 + O <=> CH3 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3 + O <=> CH3 + CO""",
)

entry(
    index = 178,
    label = "C2H3 + OH <=> C2H2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.011e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3 + OH <=> C2H2 + H2O""",
)

entry(
    index = 179,
    label = "C2H3 + O2 <=> C2H2 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.34e+06, 'cm^3/(mol*s)'),
        n = 1.61,
        Ea = (-383.4, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H3 + O2 <=> C2H2 + HO2""",
)

entry(
    index = 180,
    label = "C2H3 + O2 <=> CH2CHO + O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+11, 'cm^3/(mol*s)'), n=0.29, Ea=(11, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3 + O2 <=> CH2CHO + O""",
)

entry(
    index = 181,
    label = "C2H3 + O2 <=> HCO + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.6e+16, 'cm^3/(mol*s)'),
        n = -1.39,
        Ea = (1010, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H3 + O2 <=> HCO + CH2O""",
)

entry(
    index = 182,
    label = "C2H3 + HO2 <=> CH2CHO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3 + HO2 <=> CH2CHO + OH""",
)

entry(
    index = 183,
    label = "C2H3 + H2O2 <=> C2H4 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.21e+10, 'cm^3/(mol*s)'), n=0, Ea=(-596, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3 + H2O2 <=> C2H4 + HO2""",
)

entry(
    index = 184,
    label = "C2H3 + HCO <=> C2H4 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.033e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3 + HCO <=> C2H4 + CO""",
)

entry(
    index = 185,
    label = "C2H3 + HCO <=> C2H3CHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3 + HCO <=> C2H3CHO""",
)

entry(
    index = 186,
    label = "C2H3 + CH3 <=> C2H2 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.92e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3 + CH3 <=> C2H2 + CH4""",
)

entry(
    index = 187,
    label = "C2H3 + CH3 <=> C3H6",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2.5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4.27e+58, 'cm^6/(mol^2*s)'),
            n = -11.94,
            Ea = (9769.8, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.175,
        T3 = (1340.6, 'K'),
        T1 = (60000, 'K'),
        T2 = (10139.8, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, 'C=C': 3, 'C#C': 3, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C2H3 + CH3 <=> C3H6""",
)

entry(
    index = 188,
    label = "C2H3 + CH3 <=> aC3H5 + H",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.5e+24, 'cm^3/(mol*s)'),
        n = -2.83,
        Ea = (18618, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H3 + CH3 <=> aC3H5 + H""",
)

entry(
    index = 189,
    label = "C2H3 + C2H2 <=> C4H4 + H",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2e+18, 'cm^3/(mol*s)'),
        n = -1.68,
        Ea = (10600, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H3 + C2H2 <=> C4H4 + H""",
)

entry(
    index = 190,
    label = "C2H3 + C2H2 <=> nC4H5",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (9.3e+38, 'cm^3/(mol*s)'),
        n = -8.76,
        Ea = (12000, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H3 + C2H2 <=> nC4H5""",
)

entry(
    index = 191,
    label = "C2H3 + C2H2 <=> iC4H5",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.6e+46, 'cm^3/(mol*s)'),
        n = -10.98,
        Ea = (18600, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H3 + C2H2 <=> iC4H5""",
)

entry(
    index = 192,
    label = "C2H3 + C2H3 <=> C4H6",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.5e+42, 'cm^3/(mol*s)'),
        n = -8.84,
        Ea = (12483, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H3 + C2H3 <=> C4H6""",
)

entry(
    index = 193,
    label = "C2H3 + C2H3 <=> iC4H5 + H",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.2e+22, 'cm^3/(mol*s)'),
        n = -2.44,
        Ea = (13654, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H3 + C2H3 <=> iC4H5 + H""",
)

entry(
    index = 194,
    label = "C2H3 + C2H3 <=> nC4H5 + H",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.4e+20, 'cm^3/(mol*s)'),
        n = -2.04,
        Ea = (15361, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H3 + C2H3 <=> nC4H5 + H""",
)

entry(
    index = 195,
    label = "C2H3 + C2H3 <=> C2H2 + C2H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.6e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3 + C2H3 <=> C2H2 + C2H4""",
)

entry(
    index = 196,
    label = "CH2CHO <=> CH3 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.8e+41, 's^-1'), n=-9.147, Ea=(46900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHO <=> CH3 + CO""",
)

entry(
    index = 197,
    label = "CH2CHO + H <=> CH3CHO",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.2e+39, 'cm^6/(mol^2*s)'),
            n = -7.297,
            Ea = (4700, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.55,
        T3 = (8900, 'K'),
        T1 = (4350, 'K'),
        T2 = (7244, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, 'C=C': 3, 'C#C': 3},
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHO + H <=> CH3CHO""",
)

entry(
    index = 198,
    label = "CH2CHO + H <=> CH3CO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHO + H <=> CH3CO + H""",
)

entry(
    index = 199,
    label = "CH2CHO + H <=> CH3 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHO + H <=> CH3 + HCO""",
)

entry(
    index = 200,
    label = "CH2CHO + H <=> CH2CO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(4000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHO + H <=> CH2CO + H2""",
)

entry(
    index = 201,
    label = "CH2CHO + O <=> CH2CO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(4000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHO + O <=> CH2CO + OH""",
)

entry(
    index = 202,
    label = "CH2CHO + OH <=> CH2CO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(2000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHO + OH <=> CH2CO + H2O""",
)

entry(
    index = 203,
    label = "CH2CHO + O2 <=> CH2CO + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.4e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHO + O2 <=> CH2CO + HO2""",
)

entry(
    index = 204,
    label = "CH2CHO + O2 <=> CH2O + CO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+10, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHO + O2 <=> CH2O + CO + OH""",
)

entry(
    index = 205,
    label = "CH3 + CO <=> CH3CO",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(
            A = (4.85e+07, 'cm^3/(mol*s)'),
            n = 1.65,
            Ea = (6150, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        arrheniusLow = Arrhenius(
            A = (7.8e+30, 'cm^6/(mol^2*s)'),
            n = -5.395,
            Ea = (8600, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.258,
        T3 = (598, 'K'),
        T1 = (21002, 'K'),
        T2 = (1773, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, 'C=C': 3, 'C#C': 3, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is CH3 + CO <=> CH3CO""",
)

entry(
    index = 206,
    label = "CH3CO + H <=> CH3CHO",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(9.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.85e+44, 'cm^6/(mol^2*s)'),
            n = -8.569,
            Ea = (5500, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 1,
        T3 = (2900, 'K'),
        T1 = (2900, 'K'),
        T2 = (5132, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, 'C=C': 3, 'C#C': 3},
    ),
    shortDesc = u"""The chemkin file reaction is CH3CO + H <=> CH3CHO""",
)

entry(
    index = 207,
    label = "CH3CO + H <=> CH3 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CO + H <=> CH3 + HCO""",
)

entry(
    index = 208,
    label = "CH3CO + O <=> CH2CO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.9e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CO + O <=> CH2CO + OH""",
)

entry(
    index = 209,
    label = "CH3CO + O <=> CH3 + CO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.5e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CO + O <=> CH3 + CO2""",
)

entry(
    index = 210,
    label = "CH3CO + OH <=> CH2CO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CO + OH <=> CH2CO + H2O""",
)

entry(
    index = 211,
    label = "CH3CO + OH <=> CH3 + CO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CO + OH <=> CH3 + CO + OH""",
)

entry(
    index = 212,
    label = "CH3CO + HO2 <=> CH3 + CO2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CO + HO2 <=> CH3 + CO2 + OH""",
)

entry(
    index = 213,
    label = "CH3CO + H2O2 <=> CH3CHO + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+11, 'cm^3/(mol*s)'), n=0, Ea=(8226, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CO + H2O2 <=> CH3CHO + HO2""",
)

entry(
    index = 214,
    label = "CH3 + HCO <=> CH3CHO",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (2.2e+48, 'cm^6/(mol^2*s)'),
            n = -9.588,
            Ea = (5100, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.6173,
        T3 = (13.076, 'K'),
        T1 = (2078, 'K'),
        T2 = (5093, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, 'C=C': 3, 'C#C': 3},
    ),
    shortDesc = u"""The chemkin file reaction is CH3 + HCO <=> CH3CHO""",
)

entry(
    index = 215,
    label = "CH3CHO + H <=> CH3CO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.1e+09, 'cm^3/(mol*s)'),
        n = 1.16,
        Ea = (2400, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3CHO + H <=> CH3CO + H2""",
)

entry(
    index = 216,
    label = "CH3CHO + H <=> CH4 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+10, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHO + H <=> CH4 + HCO""",
)

entry(
    index = 217,
    label = "CH3CHO + O <=> CH3CO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(1800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHO + O <=> CH3CO + OH""",
)

entry(
    index = 218,
    label = "CH3CHO + OH <=> CH3CO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.35e+10, 'cm^3/(mol*s)'),
        n = 0.73,
        Ea = (-1110, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3CHO + OH <=> CH3CO + H2O""",
)

entry(
    index = 219,
    label = "CH3CHO + CH3 <=> CH3CO + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e-06, 'cm^3/(mol*s)'), n=5.6, Ea=(2460, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHO + CH3 <=> CH3CO + CH4""",
)

entry(
    index = 220,
    label = "CH3CHO + HCO <=> CO + HCO + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+12, 'cm^3/(mol*s)'), n=0, Ea=(10400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHO + HCO <=> CO + HCO + CH4""",
)

entry(
    index = 221,
    label = "CH3CHO + O2 <=> CH3CO + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(39100, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHO + O2 <=> CH3CO + HO2""",
)

entry(
    index = 222,
    label = "CH2OCH2 <=> CH3 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.63e+13, 's^-1'), n=0, Ea=(57200, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2OCH2 <=> CH3 + HCO""",
)

entry(
    index = 223,
    label = "CH2OCH2 <=> CH3CHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.26e+13, 's^-1'), n=0, Ea=(57200, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2OCH2 <=> CH3CHO""",
)

entry(
    index = 224,
    label = "CH2OCH2 <=> CH4 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.21e+13, 's^-1'), n=0, Ea=(57200, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2OCH2 <=> CH4 + CO""",
)

entry(
    index = 225,
    label = "CH2OCH2 + H <=> CH2OCH + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(8300, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2OCH2 + H <=> CH2OCH + H2""",
)

entry(
    index = 226,
    label = "CH2OCH2 + H <=> C2H3 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+09, 'cm^3/(mol*s)'), n=0, Ea=(5000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2OCH2 + H <=> C2H3 + H2O""",
)

entry(
    index = 227,
    label = "CH2OCH2 + H <=> C2H4 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.51e+10, 'cm^3/(mol*s)'), n=0, Ea=(5000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2OCH2 + H <=> C2H4 + OH""",
)

entry(
    index = 228,
    label = "CH2OCH2 + O <=> CH2OCH + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.91e+12, 'cm^3/(mol*s)'), n=0, Ea=(5250, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2OCH2 + O <=> CH2OCH + OH""",
)

entry(
    index = 229,
    label = "CH2OCH2 + OH <=> CH2OCH + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.78e+13, 'cm^3/(mol*s)'), n=0, Ea=(3610, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2OCH2 + OH <=> CH2OCH + H2O""",
)

entry(
    index = 230,
    label = "CH2OCH2 + CH3 <=> CH2OCH + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.07e+12, 'cm^3/(mol*s)'), n=0, Ea=(11830, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2OCH2 + CH3 <=> CH2OCH + CH4""",
)

entry(
    index = 231,
    label = "CH2OCH <=> CH3 + CO",
    degeneracy = 1,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(A=(3.16e+14, 'cm^3/(mol*s)'), n=0, Ea=(12000, 'cal/mol'), T0=(1, 'K')),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is CH2OCH <=> CH3 + CO""",
)

entry(
    index = 232,
    label = "CH2OCH <=> CH2CHO",
    degeneracy = 1,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(A=(5e+09, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is CH2OCH <=> CH2CHO""",
)

entry(
    index = 233,
    label = "CH2OCH <=> CH2CO + H",
    degeneracy = 1,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(8000, 'cal/mol'), T0=(1, 'K')),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is CH2OCH <=> CH2CO + H""",
)

entry(
    index = 234,
    label = "C2H4 <=> H2 + H2CC",
    degeneracy = 1,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(A=(7.86e+14, 'cm^3/(mol*s)'), n=0, Ea=(54245, 'cal/mol'), T0=(1, 'K')),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C2H4 <=> H2 + H2CC""",
)

entry(
    index = 235,
    label = "C2H4 + H <=> C2H5",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(
            A = (1.367e+09, 'cm^3/(mol*s)'),
            n = 1.463,
            Ea = (1355, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        arrheniusLow = Arrhenius(
            A = (2.027e+39, 'cm^6/(mol^2*s)'),
            n = -6.642,
            Ea = (5769, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -0.569,
        T3 = (299, 'K'),
        T1 = (9147, 'K'),
        T2 = (-152.4, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C2H4 + H <=> C2H5""",
)

entry(
    index = 236,
    label = "C2H4 + H <=> C2H3 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.07e+07, 'cm^3/(mol*s)'),
        n = 1.9,
        Ea = (12950, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H4 + H <=> C2H3 + H2""",
)

entry(
    index = 237,
    label = "C2H4 + O <=> C2H3 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.51e+07, 'cm^3/(mol*s)'),
        n = 1.9,
        Ea = (3740, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H4 + O <=> C2H3 + OH""",
)

entry(
    index = 238,
    label = "C2H4 + O <=> CH3 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.92e+07, 'cm^3/(mol*s)'),
        n = 1.83,
        Ea = (220, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H4 + O <=> CH3 + HCO""",
)

entry(
    index = 239,
    label = "C2H4 + O <=> CH2 + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(384000, 'cm^3/(mol*s)'), n=1.83, Ea=(220, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4 + O <=> CH2 + CH2O""",
)

entry(
    index = 240,
    label = "C2H4 + OH <=> C2H3 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.6e+06, 'cm^3/(mol*s)'), n=2, Ea=(2500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4 + OH <=> C2H3 + H2O""",
)

entry(
    index = 241,
    label = "C2H4 + HCO <=> C2H5 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+07, 'cm^3/(mol*s)'), n=2, Ea=(8000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4 + HCO <=> C2H5 + CO""",
)

entry(
    index = 242,
    label = "C2H4 + CH <=> aC3H4 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4 + CH <=> aC3H4 + H""",
)

entry(
    index = 243,
    label = "C2H4 + CH <=> pC3H4 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4 + CH <=> pC3H4 + H""",
)

entry(
    index = 244,
    label = "C2H4 + CH2 <=> aC3H5 + H",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(6000, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H4 + CH2 <=> aC3H5 + H""",
)

entry(
    index = 245,
    label = "C2H4 + CH2* <=> H2CC + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4 + CH2* <=> H2CC + CH4""",
)

entry(
    index = 246,
    label = "C2H4 + CH3 <=> C2H3 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(227000, 'cm^3/(mol*s)'), n=2, Ea=(9200, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4 + CH3 <=> C2H3 + CH4""",
)

entry(
    index = 247,
    label = "C2H4 + CH3 <=> nC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.3e+11, 'cm^3/(mol*s)'), n=0, Ea=(7700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4 + CH3 <=> nC3H7""",
)

entry(
    index = 248,
    label = "C2H4 + C2H <=> C4H4 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4 + C2H <=> C4H4 + H""",
)

entry(
    index = 249,
    label = "C2H4 + O2 <=> C2H3 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.22e+13, 'cm^3/(mol*s)'), n=0, Ea=(60800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4 + O2 <=> C2H3 + HO2""",
)

entry(
    index = 250,
    label = "C2H4 + C2H3 <=> C4H7",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (7.93e+38, 'cm^3/(mol*s)'),
        n = -8.47,
        Ea = (14220, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H4 + C2H3 <=> C4H7""",
)

entry(
    index = 251,
    label = "C2H4 + HO2 <=> CH2OCH2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.82e+12, 'cm^3/(mol*s)'), n=0, Ea=(17100, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4 + HO2 <=> CH2OCH2 + OH""",
)

entry(
    index = 252,
    label = "C2H5 + H <=> C2H6",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(
            A = (5.21e+17, 'cm^3/(mol*s)'),
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
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C2H5 + H <=> C2H6""",
)

entry(
    index = 253,
    label = "C2H5 + H <=> C2H4 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5 + H <=> C2H4 + H2""",
)

entry(
    index = 254,
    label = "C2H5 + O <=> CH3 + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.604e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5 + O <=> CH3 + CH2O""",
)

entry(
    index = 255,
    label = "C2H5 + O <=> CH3CHO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.02e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5 + O <=> CH3CHO + H""",
)

entry(
    index = 256,
    label = "C2H5 + O2 <=> C2H4 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+10, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5 + O2 <=> C2H4 + HO2""",
)

entry(
    index = 257,
    label = "C2H5 + HO2 <=> C2H6 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5 + HO2 <=> C2H6 + O2""",
)

entry(
    index = 258,
    label = "C2H5 + HO2 <=> C2H4 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5 + HO2 <=> C2H4 + H2O2""",
)

entry(
    index = 259,
    label = "C2H5 + HO2 <=> CH3 + CH2O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5 + HO2 <=> CH3 + CH2O + OH""",
)

entry(
    index = 260,
    label = "C2H5 + H2O2 <=> C2H6 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.7e+09, 'cm^3/(mol*s)'), n=0, Ea=(974, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5 + H2O2 <=> C2H6 + HO2""",
)

entry(
    index = 261,
    label = "C2H5 + CH3 <=> C3H8",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(4.9e+14, 'cm^3/(mol*s)'), n=-0.5, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (6.8e+61, 'cm^6/(mol^2*s)'),
            n = -13.42,
            Ea = (6000, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 1,
        T3 = (1000, 'K'),
        T1 = (1433.9, 'K'),
        T2 = (5328.8, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C2H5 + CH3 <=> C3H8""",
)

entry(
    index = 262,
    label = "C2H5 + C2H3 <=> C4H81",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.55e+56, 'cm^6/(mol^2*s)'),
            n = -11.79,
            Ea = (8984.5, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.198,
        T3 = (2277.9, 'K'),
        T1 = (60000, 'K'),
        T2 = (5723.2, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C2H5 + C2H3 <=> C4H81""",
)

entry(
    index = 263,
    label = "C2H5 + C2H3 <=> aC3H5 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.9e+32, 'cm^3/(mol*s)'),
        n = -5.22,
        Ea = (19747, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H5 + C2H3 <=> aC3H5 + CH3""",
)

entry(
    index = 264,
    label = "C2H6 + H <=> C2H5 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.15e+08, 'cm^3/(mol*s)'),
        n = 1.9,
        Ea = (7530, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H6 + H <=> C2H5 + H2""",
)

entry(
    index = 265,
    label = "C2H6 + O <=> C2H5 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8.98e+07, 'cm^3/(mol*s)'),
        n = 1.92,
        Ea = (5690, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H6 + O <=> C2H5 + OH""",
)

entry(
    index = 266,
    label = "C2H6 + OH <=> C2H5 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.54e+06, 'cm^3/(mol*s)'),
        n = 2.12,
        Ea = (870, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H6 + OH <=> C2H5 + H2O""",
)

entry(
    index = 267,
    label = "C2H6 + CH2* <=> C2H5 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(-550, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H6 + CH2* <=> C2H5 + CH3""",
)

entry(
    index = 268,
    label = "C2H6 + CH3 <=> C2H5 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (6.14e+06, 'cm^3/(mol*s)'),
        n = 1.74,
        Ea = (10450, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H6 + CH3 <=> C2H5 + CH4""",
)

entry(
    index = 269,
    label = "C3H3 + H <=> pC3H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H3 + H <=> pC3H4""",
)

entry(
    index = 270,
    label = "C3H3 + H <=> aC3H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H3 + H <=> aC3H4""",
)

entry(
    index = 271,
    label = "C3H3 + O <=> CH2O + C2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H3 + O <=> CH2O + C2H""",
)

entry(
    index = 272,
    label = "C3H3 + O2 <=> CH2CO + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+10, 'cm^3/(mol*s)'), n=0, Ea=(2868, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H3 + O2 <=> CH2CO + HCO""",
)

entry(
    index = 273,
    label = "C3H3 + HO2 <=> OH + CO + C2H3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H3 + HO2 <=> OH + CO + C2H3""",
)

entry(
    index = 274,
    label = "C3H3 + HO2 <=> aC3H4 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H3 + HO2 <=> aC3H4 + O2""",
)

entry(
    index = 275,
    label = "C3H3 + HO2 <=> pC3H4 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H3 + HO2 <=> pC3H4 + O2""",
)

entry(
    index = 276,
    label = "C3H3 + HCO <=> aC3H4 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H3 + HCO <=> aC3H4 + CO""",
)

entry(
    index = 277,
    label = "C3H3 + HCO <=> pC3H4 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H3 + HCO <=> pC3H4 + CO""",
)

entry(
    index = 278,
    label = "C3H3 + HCCO <=> C4H4 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H3 + HCCO <=> C4H4 + CO""",
)

entry(
    index = 279,
    label = "C3H3 + CH <=> iC4H3 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H3 + CH <=> iC4H3 + H""",
)

entry(
    index = 280,
    label = "C3H3 + CH2 <=> C4H4 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H3 + CH2 <=> C4H4 + H""",
)

entry(
    index = 281,
    label = "C3H3 + CH3 <=> C4H612",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (2.6e+57, 'cm^6/(mol^2*s)'),
            n = -11.94,
            Ea = (9770, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.175,
        T3 = (1340.6, 'K'),
        T1 = (60000, 'K'),
        T2 = (9769.8, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C3H3 + CH3 <=> C4H612""",
)

entry(
    index = 282,
    label = "C3H3 + C2H2 <=> C5H5",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (6.87e+55, 'cm^3/(mol*s)'),
        n = -12.5,
        Ea = (42025, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C3H3 + C2H2 <=> C5H5""",
)

entry(
    index = 283,
    label = "C3H3 + C3H3 => C6H5 + H",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H3 + C3H3 => C6H5 + H""",
)

entry(
    index = 284,
    label = "C3H3 + C3H3 => C6H6",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H3 + C3H3 => C6H6""",
)

entry(
    index = 285,
    label = "C3H3 + C4H4 <=> C6H5CH2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (653000, 'cm^3/(mol*s)'),
        n = 1.28,
        Ea = (-4611, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C3H3 + C4H4 <=> C6H5CH2""",
)

entry(
    index = 286,
    label = "C3H3 + C4H6 <=> C6H5CH3 + H",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (653000, 'cm^3/(mol*s)'),
        n = 1.28,
        Ea = (-4611, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C3H3 + C4H6 <=> C6H5CH3 + H""",
)

entry(
    index = 287,
    label = "aC3H4 + H <=> C3H3 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+06, 'cm^3/(mol*s)'), n=2, Ea=(5500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is aC3H4 + H <=> C3H3 + H2""",
)

entry(
    index = 288,
    label = "aC3H4 + H <=> CH3CHCH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.4e+29, 'cm^3/(mol*s)'),
        n = -6.09,
        Ea = (16300, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is aC3H4 + H <=> CH3CHCH""",
)

entry(
    index = 289,
    label = "aC3H4 + H <=> CH3CCH2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (9.46e+42, 'cm^3/(mol*s)'),
        n = -9.43,
        Ea = (11190, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is aC3H4 + H <=> CH3CCH2""",
)

entry(
    index = 290,
    label = "aC3H4 + H <=> aC3H5",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.52e+59, 'cm^3/(mol*s)'),
        n = -13.54,
        Ea = (26949, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is aC3H4 + H <=> aC3H5""",
)

entry(
    index = 291,
    label = "aC3H4 + O <=> C2H4 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+07, 'cm^3/(mol*s)'), n=1.8, Ea=(1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is aC3H4 + O <=> C2H4 + CO""",
)

entry(
    index = 292,
    label = "aC3H4 + OH <=> C3H3 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.3e+06, 'cm^3/(mol*s)'), n=2, Ea=(2000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is aC3H4 + OH <=> C3H3 + H2O""",
)

entry(
    index = 293,
    label = "aC3H4 + CH3 <=> C3H3 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+12, 'cm^3/(mol*s)'), n=0, Ea=(7700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is aC3H4 + CH3 <=> C3H3 + CH4""",
)

entry(
    index = 294,
    label = "aC3H4 + CH3 <=> iC4H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+11, 'cm^3/(mol*s)'), n=0, Ea=(7500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is aC3H4 + CH3 <=> iC4H7""",
)

entry(
    index = 295,
    label = "aC3H4 + C2H <=> C2H2 + C3H3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is aC3H4 + C2H <=> C2H2 + C3H3""",
)

entry(
    index = 296,
    label = "pC3H4 <=> cC3H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+44, 's^-1'), n=-9.92, Ea=(69250, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is pC3H4 <=> cC3H4""",
)

entry(
    index = 297,
    label = "pC3H4 <=> aC3H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.15e+60, 's^-1'), n=-13.93, Ea=(91117, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is pC3H4 <=> aC3H4""",
)

entry(
    index = 298,
    label = "pC3H4 + H <=> aC3H4 + H",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (6.27e+17, 'cm^3/(mol*s)'),
        n = -0.91,
        Ea = (10079, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is pC3H4 + H <=> aC3H4 + H""",
)

entry(
    index = 299,
    label = "pC3H4 + H <=> CH3CCH2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.66e+47, 'cm^3/(mol*s)'),
        n = -10.58,
        Ea = (13690, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is pC3H4 + H <=> CH3CCH2""",
)

entry(
    index = 300,
    label = "pC3H4 + H <=> CH3CHCH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.5e+28, 'cm^3/(mol*s)'),
        n = -5.74,
        Ea = (4300, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is pC3H4 + H <=> CH3CHCH""",
)

entry(
    index = 301,
    label = "pC3H4 + H <=> aC3H5",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.91e+60, 'cm^3/(mol*s)'),
        n = -14.37,
        Ea = (31644, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is pC3H4 + H <=> aC3H5""",
)

entry(
    index = 302,
    label = "pC3H4 + H <=> C3H3 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+06, 'cm^3/(mol*s)'), n=2, Ea=(5500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is pC3H4 + H <=> C3H3 + H2""",
)

entry(
    index = 303,
    label = "pC3H4 + C3H3 <=> aC3H4 + C3H3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (6.14e+06, 'cm^3/(mol*s)'),
        n = 1.74,
        Ea = (10450, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is pC3H4 + C3H3 <=> aC3H4 + C3H3""",
)

entry(
    index = 304,
    label = "pC3H4 + O <=> HCCO + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.3e+12, 'cm^3/(mol*s)'), n=0, Ea=(2250, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is pC3H4 + O <=> HCCO + CH3""",
)

entry(
    index = 305,
    label = "pC3H4 + O <=> C2H4 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(2250, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is pC3H4 + O <=> C2H4 + CO""",
)

entry(
    index = 306,
    label = "pC3H4 + OH <=> C3H3 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+06, 'cm^3/(mol*s)'), n=2, Ea=(100, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is pC3H4 + OH <=> C3H3 + H2O""",
)

entry(
    index = 307,
    label = "pC3H4 + C2H <=> C2H2 + C3H3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is pC3H4 + C2H <=> C2H2 + C3H3""",
)

entry(
    index = 308,
    label = "pC3H4 + CH3 <=> C3H3 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(7700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is pC3H4 + CH3 <=> C3H3 + CH4""",
)

entry(
    index = 309,
    label = "cC3H4 <=> aC3H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.89e+41, 's^-1'), n=-9.17, Ea=(49594, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC3H4 <=> aC3H4""",
)

entry(
    index = 310,
    label = "aC3H5 + H <=> C3H6",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.33e+60, 'cm^6/(mol^2*s)'),
            n = -12,
            Ea = (5967.8, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.02,
        T3 = (1096.6, 'K'),
        T1 = (1096.6, 'K'),
        T2 = (6859.5, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is aC3H5 + H <=> C3H6""",
)

entry(
    index = 311,
    label = "aC3H5 + H <=> aC3H4 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is aC3H5 + H <=> aC3H4 + H2""",
)

entry(
    index = 312,
    label = "aC3H5 + O <=> C2H3CHO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is aC3H5 + O <=> C2H3CHO + H""",
)

entry(
    index = 313,
    label = "aC3H5 + OH <=> C2H3CHO + H + H",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.2e+32, 'cm^3/(mol*s)'),
        n = -5.16,
        Ea = (30126, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is aC3H5 + OH <=> C2H3CHO + H + H""",
)

entry(
    index = 314,
    label = "aC3H5 + OH <=> aC3H4 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is aC3H5 + OH <=> aC3H4 + H2O""",
)

entry(
    index = 315,
    label = "aC3H5 + O2 <=> aC3H4 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.99e+15, 'cm^3/(mol*s)'),
        n = -1.4,
        Ea = (22428, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is aC3H5 + O2 <=> aC3H4 + HO2""",
)

entry(
    index = 316,
    label = "aC3H5 + O2 <=> CH3CO + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.19e+15, 'cm^3/(mol*s)'),
        n = -1.01,
        Ea = (20128, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is aC3H5 + O2 <=> CH3CO + CH2O""",
)

entry(
    index = 317,
    label = "aC3H5 + O2 <=> C2H3CHO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.82e+13, 'cm^3/(mol*s)'),
        n = -0.41,
        Ea = (22859, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is aC3H5 + O2 <=> C2H3CHO + OH""",
)

entry(
    index = 318,
    label = "aC3H5 + HO2 <=> C3H6 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.66e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is aC3H5 + HO2 <=> C3H6 + O2""",
)

entry(
    index = 319,
    label = "aC3H5 + HO2 <=> OH + C2H3 + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.6e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is aC3H5 + HO2 <=> OH + C2H3 + CH2O""",
)

entry(
    index = 320,
    label = "aC3H5 + HCO <=> C3H6 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is aC3H5 + HCO <=> C3H6 + CO""",
)

entry(
    index = 321,
    label = "aC3H5 + CH3 <=> C4H81",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(
            A = (1e+14, 'cm^3/(mol*s)'),
            n = -0.32,
            Ea = (-262.3, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        arrheniusLow = Arrhenius(
            A = (3.91e+60, 'cm^6/(mol^2*s)'),
            n = -12.81,
            Ea = (6250, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.104,
        T3 = (1606, 'K'),
        T1 = (60000, 'K'),
        T2 = (6118.4, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is aC3H5 + CH3 <=> C4H81""",
)

entry(
    index = 322,
    label = "aC3H5 + CH3 <=> aC3H4 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+12, 'cm^3/(mol*s)'), n=-0.32, Ea=(-131, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is aC3H5 + CH3 <=> aC3H4 + CH4""",
)

entry(
    index = 323,
    label = "aC3H5 <=> CH3CCH2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.06e+56, 's^-1'), n=-14.08, Ea=(75868, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is aC3H5 <=> CH3CCH2""",
)

entry(
    index = 324,
    label = "aC3H5 <=> CH3CHCH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+51, 's^-1'), n=-13.02, Ea=(73300, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is aC3H5 <=> CH3CHCH""",
)

entry(
    index = 325,
    label = "aC3H5 + C2H2 <=> lC5H7",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8.38e+30, 'cm^3/(mol*s)'),
        n = -6.242,
        Ea = (12824, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is aC3H5 + C2H2 <=> lC5H7""",
)

entry(
    index = 326,
    label = "CH3CCH2 <=> CH3CHCH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.5e+48, 's^-1'), n=-12.71, Ea=(53900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CCH2 <=> CH3CHCH""",
)

entry(
    index = 327,
    label = "CH3CCH2 + H <=> pC3H4 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.34e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CCH2 + H <=> pC3H4 + H2""",
)

entry(
    index = 328,
    label = "CH3CCH2 + O <=> CH3 + CH2CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CCH2 + O <=> CH3 + CH2CO""",
)

entry(
    index = 329,
    label = "CH3CCH2 + OH <=> CH3 + CH2CO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CCH2 + OH <=> CH3 + CH2CO + H""",
)

entry(
    index = 330,
    label = "CH3CCH2 + O2 <=> CH3CO + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CCH2 + O2 <=> CH3CO + CH2O""",
)

entry(
    index = 331,
    label = "CH3CCH2 + HO2 <=> CH3 + CH2CO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CCH2 + HO2 <=> CH3 + CH2CO + OH""",
)

entry(
    index = 332,
    label = "CH3CCH2 + HCO <=> C3H6 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CCH2 + HCO <=> C3H6 + CO""",
)

entry(
    index = 333,
    label = "CH3CCH2 + CH3 <=> pC3H4 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CCH2 + CH3 <=> pC3H4 + CH4""",
)

entry(
    index = 334,
    label = "CH3CCH2 + CH3 <=> iC4H8",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CCH2 + CH3 <=> iC4H8""",
)

entry(
    index = 335,
    label = "CH3CHCH + H <=> pC3H4 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.34e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHCH + H <=> pC3H4 + H2""",
)

entry(
    index = 336,
    label = "CH3CHCH + O <=> C2H4 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHCH + O <=> C2H4 + HCO""",
)

entry(
    index = 337,
    label = "CH3CHCH + OH <=> C2H4 + HCO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHCH + OH <=> C2H4 + HCO + H""",
)

entry(
    index = 338,
    label = "CH3CHCH + O2 <=> CH3CHO + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHCH + O2 <=> CH3CHO + HCO""",
)

entry(
    index = 339,
    label = "CH3CHCH + HO2 <=> C2H4 + HCO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHCH + HO2 <=> C2H4 + HCO + OH""",
)

entry(
    index = 340,
    label = "CH3CHCH + HCO <=> C3H6 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHCH + HCO <=> C3H6 + CO""",
)

entry(
    index = 341,
    label = "CH3CHCH + CH3 <=> pC3H4 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHCH + CH3 <=> pC3H4 + CH4""",
)

entry(
    index = 342,
    label = "C3H6 + H <=> nC3H7",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(
            A = (1.33e+13, 'cm^3/(mol*s)'),
            n = 0,
            Ea = (3260.7, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        arrheniusLow = Arrhenius(
            A = (6.26e+38, 'cm^6/(mol^2*s)'),
            n = -6.66,
            Ea = (7000, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 1,
        T3 = (1000, 'K'),
        T1 = (1310, 'K'),
        T2 = (48097, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C3H6 + H <=> nC3H7""",
)

entry(
    index = 343,
    label = "C3H6 + H <=> iC3H7",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(
            A = (1.33e+13, 'cm^3/(mol*s)'),
            n = 0,
            Ea = (1559.8, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        arrheniusLow = Arrhenius(
            A = (8.7e+42, 'cm^6/(mol^2*s)'),
            n = -7.5,
            Ea = (4721.8, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 1,
        T3 = (1000, 'K'),
        T1 = (645.4, 'K'),
        T2 = (6844.3, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C3H6 + H <=> iC3H7""",
)

entry(
    index = 344,
    label = "C3H6 + H <=> C2H4 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8e+21, 'cm^3/(mol*s)'),
        n = -2.39,
        Ea = (11180, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C3H6 + H <=> C2H4 + CH3""",
)

entry(
    index = 345,
    label = "C3H6 + H <=> aC3H5 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(173000, 'cm^3/(mol*s)'), n=2.5, Ea=(2490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6 + H <=> aC3H5 + H2""",
)

entry(
    index = 346,
    label = "C3H6 + H <=> CH3CCH2 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(400000, 'cm^3/(mol*s)'), n=2.5, Ea=(9790, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6 + H <=> CH3CCH2 + H2""",
)

entry(
    index = 347,
    label = "C3H6 + H <=> CH3CHCH + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(804000, 'cm^3/(mol*s)'), n=2.5, Ea=(12283, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6 + H <=> CH3CHCH + H2""",
)

entry(
    index = 348,
    label = "C3H6 + O <=> CH2CO + CH3 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+07, 'cm^3/(mol*s)'), n=1.65, Ea=(327, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6 + O <=> CH2CO + CH3 + H""",
)

entry(
    index = 349,
    label = "C3H6 + O <=> C2H3CHO + H + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+07, 'cm^3/(mol*s)'), n=1.65, Ea=(327, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6 + O <=> C2H3CHO + H + H""",
)

entry(
    index = 350,
    label = "C3H6 + O <=> C2H5 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.5e+07, 'cm^3/(mol*s)'),
        n = 1.65,
        Ea = (-972, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C3H6 + O <=> C2H5 + HCO""",
)

entry(
    index = 351,
    label = "C3H6 + O <=> aC3H5 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+11, 'cm^3/(mol*s)'), n=0.7, Ea=(5880, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6 + O <=> aC3H5 + OH""",
)

entry(
    index = 352,
    label = "C3H6 + O <=> CH3CCH2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+10, 'cm^3/(mol*s)'), n=0.7, Ea=(7630, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6 + O <=> CH3CCH2 + OH""",
)

entry(
    index = 353,
    label = "C3H6 + O <=> CH3CHCH + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.21e+11, 'cm^3/(mol*s)'),
        n = 0.7,
        Ea = (8960, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C3H6 + O <=> CH3CHCH + OH""",
)

entry(
    index = 354,
    label = "C3H6 + OH <=> aC3H5 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.1e+06, 'cm^3/(mol*s)'), n=2, Ea=(-298, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6 + OH <=> aC3H5 + H2O""",
)

entry(
    index = 355,
    label = "C3H6 + OH <=> CH3CCH2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.1e+06, 'cm^3/(mol*s)'), n=2, Ea=(1450, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6 + OH <=> CH3CCH2 + H2O""",
)

entry(
    index = 356,
    label = "C3H6 + OH <=> CH3CHCH + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.14e+06, 'cm^3/(mol*s)'), n=2, Ea=(2778, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6 + OH <=> CH3CHCH + H2O""",
)

entry(
    index = 357,
    label = "C3H6 + HO2 <=> aC3H5 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7130, 'cm^3/(mol*s)'), n=2.77, Ea=(14913, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6 + HO2 <=> aC3H5 + H2O2""",
)

entry(
    index = 358,
    label = "C3H6 + CH3 <=> aC3H5 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.2, 'cm^3/(mol*s)'), n=3.5, Ea=(5675, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6 + CH3 <=> aC3H5 + CH4""",
)

entry(
    index = 359,
    label = "C3H6 + CH3 <=> CH3CCH2 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.84, 'cm^3/(mol*s)'), n=3.5, Ea=(11660, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6 + CH3 <=> CH3CCH2 + CH4""",
)

entry(
    index = 360,
    label = "C3H6 + CH3 <=> CH3CHCH + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.35, 'cm^3/(mol*s)'), n=3.5, Ea=(12848, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6 + CH3 <=> CH3CHCH + CH4""",
)

entry(
    index = 361,
    label = "C3H6 + C2H3 <=> C4H6 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.23e+11, 'cm^3/(mol*s)'), n=0, Ea=(5000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6 + C2H3 <=> C4H6 + CH3""",
)

entry(
    index = 362,
    label = "C3H6 + HO2 <=> CH3CHOCH2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.09e+12, 'cm^3/(mol*s)'), n=0, Ea=(14200, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6 + HO2 <=> CH3CHOCH2 + OH""",
)

entry(
    index = 363,
    label = "C2H3CHO + H <=> C2H4 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.08e+11, 'cm^3/(mol*s)'),
        n = 0.454,
        Ea = (5820, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H3CHO + H <=> C2H4 + HCO""",
)

entry(
    index = 364,
    label = "C2H3CHO + O <=> C2H3 + OH + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(3540, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3CHO + O <=> C2H3 + OH + CO""",
)

entry(
    index = 365,
    label = "C2H3CHO + O <=> CH2O + CH2CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.9e+07, 'cm^3/(mol*s)'), n=1.8, Ea=(220, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3CHO + O <=> CH2O + CH2CO""",
)

entry(
    index = 366,
    label = "C2H3CHO + OH <=> C2H3 + H2O + CO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.43e+09, 'cm^3/(mol*s)'),
        n = 1.18,
        Ea = (-447, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H3CHO + OH <=> C2H3 + H2O + CO""",
)

entry(
    index = 367,
    label = "C2H3CHO + CH3 <=> CH2CHCO + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(11000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3CHO + CH3 <=> CH2CHCO + CH4""",
)

entry(
    index = 368,
    label = "C2H3CHO + C2H3 <=> C4H6 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.8e+21, 'cm^3/(mol*s)'),
        n = -2.44,
        Ea = (14720, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H3CHO + C2H3 <=> C4H6 + HCO""",
)

entry(
    index = 369,
    label = "CH2CHCO <=> C2H3 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 's^-1'), n=0, Ea=(27000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHCO <=> C2H3 + CO""",
)

entry(
    index = 370,
    label = "CH2CHCO + H <=> C2H3CHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHCO + H <=> C2H3CHO""",
)

entry(
    index = 371,
    label = "CH3CHOCH2 <=> CH3CH2CHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.84e+14, 's^-1'), n=0, Ea=(58500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHOCH2 <=> CH3CH2CHO""",
)

entry(
    index = 372,
    label = "CH3CHOCH2 <=> C2H5 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.45e+13, 's^-1'), n=0, Ea=(58500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHOCH2 <=> C2H5 + HCO""",
)

entry(
    index = 373,
    label = "CH3CHOCH2 <=> CH3 + CH2CHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.45e+13, 's^-1'), n=0, Ea=(58800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHOCH2 <=> CH3 + CH2CHO""",
)

entry(
    index = 374,
    label = "CH3CHOCH2 <=> CH3COCH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.01e+14, 's^-1'), n=0, Ea=(59900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHOCH2 <=> CH3COCH3""",
)

entry(
    index = 375,
    label = "CH3CHOCH2 <=> CH3 + CH3CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.54e+13, 's^-1'), n=0, Ea=(59900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHOCH2 <=> CH3 + CH3CO""",
)

entry(
    index = 376,
    label = "iC3H7 + H <=> C3H8",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.7e+58, 'cm^6/(mol^2*s)'),
            n = -12.08,
            Ea = (11263.7, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.649,
        T3 = (1213.1, 'K'),
        T1 = (1213.1, 'K'),
        T2 = (13369.7, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is iC3H7 + H <=> C3H8""",
)

entry(
    index = 377,
    label = "iC3H7 + H <=> CH3 + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.4e+28, 'cm^3/(mol*s)'),
        n = -3.94,
        Ea = (15916, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is iC3H7 + H <=> CH3 + C2H5""",
)

entry(
    index = 378,
    label = "iC3H7 + H <=> C3H6 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.2e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is iC3H7 + H <=> C3H6 + H2""",
)

entry(
    index = 379,
    label = "iC3H7 + O <=> CH3CHO + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is iC3H7 + O <=> CH3CHO + CH3""",
)

entry(
    index = 380,
    label = "iC3H7 + OH <=> C3H6 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is iC3H7 + OH <=> C3H6 + H2O""",
)

entry(
    index = 381,
    label = "iC3H7 + O2 <=> C3H6 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is iC3H7 + O2 <=> C3H6 + HO2""",
)

entry(
    index = 382,
    label = "iC3H7 + HO2 <=> CH3CHO + CH3 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is iC3H7 + HO2 <=> CH3CHO + CH3 + OH""",
)

entry(
    index = 383,
    label = "iC3H7 + HCO <=> C3H8 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is iC3H7 + HCO <=> C3H8 + CO""",
)

entry(
    index = 384,
    label = "iC3H7 + CH3 <=> CH4 + C3H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.2e+14, 'cm^3/(mol*s)'), n=-0.68, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is iC3H7 + CH3 <=> CH4 + C3H6""",
)

entry(
    index = 385,
    label = "nC3H7 + H <=> C3H8",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.01e+48, 'cm^6/(mol^2*s)'),
            n = -9.32,
            Ea = (5833.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.498,
        T3 = (1314, 'K'),
        T1 = (1314, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is nC3H7 + H <=> C3H8""",
)

entry(
    index = 386,
    label = "nC3H7 + H <=> C2H5 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.7e+24, 'cm^3/(mol*s)'),
        n = -2.92,
        Ea = (12505, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is nC3H7 + H <=> C2H5 + CH3""",
)

entry(
    index = 387,
    label = "nC3H7 + H <=> C3H6 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is nC3H7 + H <=> C3H6 + H2""",
)

entry(
    index = 388,
    label = "nC3H7 + O <=> C2H5 + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is nC3H7 + O <=> C2H5 + CH2O""",
)

entry(
    index = 389,
    label = "nC3H7 + OH <=> C3H6 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is nC3H7 + OH <=> C3H6 + H2O""",
)

entry(
    index = 390,
    label = "nC3H7 + O2 <=> C3H6 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+10, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is nC3H7 + O2 <=> C3H6 + HO2""",
)

entry(
    index = 391,
    label = "nC3H7 + HO2 <=> C2H5 + OH + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is nC3H7 + HO2 <=> C2H5 + OH + CH2O""",
)

entry(
    index = 392,
    label = "nC3H7 + HCO <=> C3H8 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is nC3H7 + HCO <=> C3H8 + CO""",
)

entry(
    index = 393,
    label = "nC3H7 + CH3 <=> CH4 + C3H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is nC3H7 + CH3 <=> CH4 + C3H6""",
)

entry(
    index = 394,
    label = "C3H8 + H <=> H2 + nC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.3e+06, 'cm^3/(mol*s)'),
        n = 2.54,
        Ea = (6756, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C3H8 + H <=> H2 + nC3H7""",
)

entry(
    index = 395,
    label = "C3H8 + H <=> H2 + iC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+06, 'cm^3/(mol*s)'), n=2.4, Ea=(4471, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H8 + H <=> H2 + iC3H7""",
)

entry(
    index = 396,
    label = "C3H8 + O <=> nC3H7 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(190000, 'cm^3/(mol*s)'), n=2.68, Ea=(3716, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H8 + O <=> nC3H7 + OH""",
)

entry(
    index = 397,
    label = "C3H8 + O <=> iC3H7 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(47600, 'cm^3/(mol*s)'), n=2.71, Ea=(2106, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H8 + O <=> iC3H7 + OH""",
)

entry(
    index = 398,
    label = "C3H8 + OH <=> nC3H7 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1400, 'cm^3/(mol*s)'), n=2.66, Ea=(527, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H8 + OH <=> nC3H7 + H2O""",
)

entry(
    index = 399,
    label = "C3H8 + OH <=> iC3H7 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(27000, 'cm^3/(mol*s)'), n=2.39, Ea=(393, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H8 + OH <=> iC3H7 + H2O""",
)

entry(
    index = 400,
    label = "C3H8 + O2 <=> nC3H7 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(50930, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H8 + O2 <=> nC3H7 + HO2""",
)

entry(
    index = 401,
    label = "C3H8 + O2 <=> iC3H7 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(47590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H8 + O2 <=> iC3H7 + HO2""",
)

entry(
    index = 402,
    label = "C3H8 + HO2 <=> nC3H7 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(61100, 'cm^3/(mol*s)'), n=2.65, Ea=(17496, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H8 + HO2 <=> nC3H7 + H2O2""",
)

entry(
    index = 403,
    label = "C3H8 + HO2 <=> iC3H7 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7130, 'cm^3/(mol*s)'), n=2.77, Ea=(14913, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H8 + HO2 <=> iC3H7 + H2O2""",
)

entry(
    index = 404,
    label = "C3H8 + CH3 <=> CH4 + nC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.903, 'cm^3/(mol*s)'), n=3.65, Ea=(7153, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H8 + CH3 <=> CH4 + nC3H7""",
)

entry(
    index = 405,
    label = "C3H8 + CH3 <=> CH4 + iC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.51, 'cm^3/(mol*s)'), n=3.46, Ea=(5480, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H8 + CH3 <=> CH4 + iC3H7""",
)

entry(
    index = 406,
    label = "C4H2 + H <=> nC4H3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.1e+42, 'cm^3/(mol*s)'),
        n = -8.72,
        Ea = (15300, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C4H2 + H <=> nC4H3""",
)

entry(
    index = 407,
    label = "C4H2 + H <=> iC4H3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.1e+30, 'cm^3/(mol*s)'),
        n = -4.92,
        Ea = (10800, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C4H2 + H <=> iC4H3""",
)

entry(
    index = 408,
    label = "C4H2 + OH <=> H2C4O + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.6e+12, 'cm^3/(mol*s)'), n=0, Ea=(-410, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H2 + OH <=> H2C4O + H""",
)

entry(
    index = 409,
    label = "C4H2 + C2H <=> C6H2 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H2 + C2H <=> C6H2 + H""",
)

entry(
    index = 410,
    label = "C4H2 + C2H <=> C6H3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.5e+37, 'cm^3/(mol*s)'),
        n = -7.68,
        Ea = (7100, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C4H2 + C2H <=> C6H3""",
)

entry(
    index = 411,
    label = "H2C4O + H <=> C2H2 + HCCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(3000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2C4O + H <=> C2H2 + HCCO""",
)

entry(
    index = 412,
    label = "H2C4O + OH <=> CH2CO + HCCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+07, 'cm^3/(mol*s)'), n=2, Ea=(2000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2C4O + OH <=> CH2CO + HCCO""",
)

entry(
    index = 413,
    label = "nC4H3 <=> iC4H3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.1e+43, 's^-1'), n=-9.49, Ea=(53000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is nC4H3 <=> iC4H3""",
)

entry(
    index = 414,
    label = "nC4H3 + H <=> iC4H3 + H",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.5e+20, 'cm^3/(mol*s)'),
        n = -1.67,
        Ea = (10800, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is nC4H3 + H <=> iC4H3 + H""",
)

entry(
    index = 415,
    label = "nC4H3 + H <=> C2H2 + H2CC",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (6.3e+25, 'cm^3/(mol*s)'),
        n = -3.34,
        Ea = (10014, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is nC4H3 + H <=> C2H2 + H2CC""",
)

entry(
    index = 416,
    label = "nC4H3 + H <=> C4H4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2e+47, 'cm^3/(mol*s)'),
        n = -10.26,
        Ea = (13070, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is nC4H3 + H <=> C4H4""",
)

entry(
    index = 417,
    label = "nC4H3 + H <=> C4H2 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is nC4H3 + H <=> C4H2 + H2""",
)

entry(
    index = 418,
    label = "nC4H3 + OH <=> C4H2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is nC4H3 + OH <=> C4H2 + H2O""",
)

entry(
    index = 419,
    label = "nC4H3 + C2H2 <=> l-C6H4 + H",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.5e+14, 'cm^3/(mol*s)'),
        n = -0.56,
        Ea = (10600, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is nC4H3 + C2H2 <=> l-C6H4 + H""",
)

entry(
    index = 420,
    label = "nC4H3 + C2H2 <=> C6H5",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (9.6e+70, 'cm^3/(mol*s)'),
        n = -17.77,
        Ea = (31300, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is nC4H3 + C2H2 <=> C6H5""",
)

entry(
    index = 421,
    label = "nC4H3 + C2H2 <=> o-C6H4 + H",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (6.9e+46, 'cm^3/(mol*s)'),
        n = -10.01,
        Ea = (30100, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is nC4H3 + C2H2 <=> o-C6H4 + H""",
)

entry(
    index = 422,
    label = "iC4H3 + H <=> C2H2 + H2CC",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.8e+23, 'cm^3/(mol*s)'),
        n = -2.55,
        Ea = (10780, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is iC4H3 + H <=> C2H2 + H2CC""",
)

entry(
    index = 423,
    label = "iC4H3 + H <=> C4H4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.4e+43, 'cm^3/(mol*s)'),
        n = -9.01,
        Ea = (12120, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is iC4H3 + H <=> C4H4""",
)

entry(
    index = 424,
    label = "iC4H3 + H <=> C4H2 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is iC4H3 + H <=> C4H2 + H2""",
)

entry(
    index = 425,
    label = "iC4H3 + OH <=> C4H2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is iC4H3 + OH <=> C4H2 + H2O""",
)

entry(
    index = 426,
    label = "iC4H3 + O2 <=> HCCO + CH2CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.86e+16, 'cm^3/(mol*s)'), n=-1.8, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is iC4H3 + O2 <=> HCCO + CH2CO""",
)

entry(
    index = 427,
    label = "C4H4 + H <=> nC4H5",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.3e+51, 'cm^3/(mol*s)'),
        n = -11.92,
        Ea = (16500, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C4H4 + H <=> nC4H5""",
)

entry(
    index = 428,
    label = "C4H4 + H <=> iC4H5",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.9e+51, 'cm^3/(mol*s)'),
        n = -11.92,
        Ea = (17700, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C4H4 + H <=> iC4H5""",
)

entry(
    index = 429,
    label = "C4H4 + H <=> nC4H3 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (665000, 'cm^3/(mol*s)'),
        n = 2.53,
        Ea = (12240, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C4H4 + H <=> nC4H3 + H2""",
)

entry(
    index = 430,
    label = "C4H4 + H <=> iC4H3 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(333000, 'cm^3/(mol*s)'), n=2.53, Ea=(9240, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H4 + H <=> iC4H3 + H2""",
)

entry(
    index = 431,
    label = "C4H4 + OH <=> nC4H3 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.1e+07, 'cm^3/(mol*s)'), n=2, Ea=(3430, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H4 + OH <=> nC4H3 + H2O""",
)

entry(
    index = 432,
    label = "C4H4 + OH <=> iC4H3 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.55e+07, 'cm^3/(mol*s)'), n=2, Ea=(430, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H4 + OH <=> iC4H3 + H2O""",
)

entry(
    index = 433,
    label = "C4H4 + O <=> C3H3 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+08, 'cm^3/(mol*s)'), n=1.45, Ea=(-860, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H4 + O <=> C3H3 + HCO""",
)

entry(
    index = 434,
    label = "C4H4 + C2H <=> l-C6H4 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H4 + C2H <=> l-C6H4 + H""",
)

entry(
    index = 435,
    label = "nC4H5 <=> iC4H5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.5e+67, 's^-1'), n=-16.89, Ea=(59100, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is nC4H5 <=> iC4H5""",
)

entry(
    index = 436,
    label = "nC4H5 + H <=> iC4H5 + H",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.1e+26, 'cm^3/(mol*s)'),
        n = -3.35,
        Ea = (17423, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is nC4H5 + H <=> iC4H5 + H""",
)

entry(
    index = 437,
    label = "nC4H5 + H <=> C4H4 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is nC4H5 + H <=> C4H4 + H2""",
)

entry(
    index = 438,
    label = "nC4H5 + OH <=> C4H4 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is nC4H5 + OH <=> C4H4 + H2O""",
)

entry(
    index = 439,
    label = "nC4H5 + HCO <=> C4H6 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is nC4H5 + HCO <=> C4H6 + CO""",
)

entry(
    index = 440,
    label = "nC4H5 + HO2 <=> C2H3 + CH2CO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.6e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is nC4H5 + HO2 <=> C2H3 + CH2CO + OH""",
)

entry(
    index = 441,
    label = "nC4H5 + H2O2 <=> C4H6 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.21e+10, 'cm^3/(mol*s)'), n=0, Ea=(-596, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is nC4H5 + H2O2 <=> C4H6 + HO2""",
)

entry(
    index = 442,
    label = "nC4H5 + HO2 <=> C4H6 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is nC4H5 + HO2 <=> C4H6 + O2""",
)

entry(
    index = 443,
    label = "nC4H5 + O2 <=> CH2CHCHCHO + O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+11, 'cm^3/(mol*s)'), n=0.29, Ea=(11, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is nC4H5 + O2 <=> CH2CHCHCHO + O""",
)

entry(
    index = 444,
    label = "nC4H5 + O2 <=> HCO + C2H3CHO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (9.2e+16, 'cm^3/(mol*s)'),
        n = -1.39,
        Ea = (1010, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is nC4H5 + O2 <=> HCO + C2H3CHO""",
)

entry(
    index = 445,
    label = "nC4H5 + C2H2 <=> C6H6 + H",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.6e+16, 'cm^3/(mol*s)'),
        n = -1.33,
        Ea = (5400, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is nC4H5 + C2H2 <=> C6H6 + H""",
)

entry(
    index = 446,
    label = "nC4H5 + C2H3 <=> C6H6 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.84e-13, 'cm^3/(mol*s)'),
        n = 7.07,
        Ea = (-3611, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is nC4H5 + C2H3 <=> C6H6 + H2""",
)

entry(
    index = 447,
    label = "iC4H5 + H <=> C4H4 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is iC4H5 + H <=> C4H4 + H2""",
)

entry(
    index = 448,
    label = "iC4H5 + H <=> C3H3 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(2000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is iC4H5 + H <=> C3H3 + CH3""",
)

entry(
    index = 449,
    label = "iC4H5 + OH <=> C4H4 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is iC4H5 + OH <=> C4H4 + H2O""",
)

entry(
    index = 450,
    label = "iC4H5 + HCO <=> C4H6 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is iC4H5 + HCO <=> C4H6 + CO""",
)

entry(
    index = 451,
    label = "iC4H5 + HO2 <=> C4H6 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is iC4H5 + HO2 <=> C4H6 + O2""",
)

entry(
    index = 452,
    label = "iC4H5 + HO2 <=> C2H3 + CH2CO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.6e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is iC4H5 + HO2 <=> C2H3 + CH2CO + OH""",
)

entry(
    index = 453,
    label = "iC4H5 + H2O2 <=> C4H6 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.21e+10, 'cm^3/(mol*s)'), n=0, Ea=(-596, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is iC4H5 + H2O2 <=> C4H6 + HO2""",
)

entry(
    index = 454,
    label = "iC4H5 + O2 <=> CH2CO + CH2CHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.16e+10, 'cm^3/(mol*s)'), n=0, Ea=(2500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is iC4H5 + O2 <=> CH2CO + CH2CHO""",
)

entry(
    index = 455,
    label = "C4H5-2 <=> iC4H5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.5e+67, 's^-1'), n=-16.89, Ea=(59100, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H5-2 <=> iC4H5""",
)

entry(
    index = 456,
    label = "iC4H5 + H <=> C4H5-2 + H",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.1e+26, 'cm^3/(mol*s)'),
        n = -3.35,
        Ea = (17423, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is iC4H5 + H <=> C4H5-2 + H""",
)

entry(
    index = 457,
    label = "C4H5-2 + HO2 <=> OH + C2H2 + CH3CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H5-2 + HO2 <=> OH + C2H2 + CH3CO""",
)

entry(
    index = 458,
    label = "C4H5-2 + O2 <=> CH3CO + CH2CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.16e+10, 'cm^3/(mol*s)'), n=0, Ea=(2500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H5-2 + O2 <=> CH3CO + CH2CO""",
)

entry(
    index = 459,
    label = "C4H5-2 + C2H2 <=> C6H6 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+14, 'cm^3/(mol*s)'), n=0, Ea=(25000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H5-2 + C2H2 <=> C6H6 + H""",
)

entry(
    index = 460,
    label = "C4H5-2 + C2H4 <=> C5H6 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+14, 'cm^3/(mol*s)'), n=0, Ea=(25000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H5-2 + C2H4 <=> C5H6 + CH3""",
)

entry(
    index = 461,
    label = "C4H6 <=> iC4H5 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.7e+36, 's^-1'), n=-6.27, Ea=(112353, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H6 <=> iC4H5 + H""",
)

entry(
    index = 462,
    label = "C4H6 <=> nC4H5 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.3e+44, 's^-1'), n=-8.62, Ea=(123608, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H6 <=> nC4H5 + H""",
)

entry(
    index = 463,
    label = "C4H6 <=> C4H4 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.5e+15, 's^-1'), n=0, Ea=(94700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H6 <=> C4H4 + H2""",
)

entry(
    index = 464,
    label = "C4H6 + H <=> nC4H5 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.33e+06, 'cm^3/(mol*s)'),
        n = 2.53,
        Ea = (12240, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C4H6 + H <=> nC4H5 + H2""",
)

entry(
    index = 465,
    label = "C4H6 + H <=> iC4H5 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(665000, 'cm^3/(mol*s)'), n=2.53, Ea=(9240, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H6 + H <=> iC4H5 + H2""",
)

entry(
    index = 466,
    label = "C4H6 + H <=> C2H4 + C2H3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.46e+30, 'cm^3/(mol*s)'),
        n = -4.34,
        Ea = (21647, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C4H6 + H <=> C2H4 + C2H3""",
)

entry(
    index = 467,
    label = "C4H6 + H <=> pC3H4 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+12, 'cm^3/(mol*s)'), n=0, Ea=(7000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H6 + H <=> pC3H4 + CH3""",
)

entry(
    index = 468,
    label = "C4H6 + H <=> aC3H4 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+12, 'cm^3/(mol*s)'), n=0, Ea=(7000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H6 + H <=> aC3H4 + CH3""",
)

entry(
    index = 469,
    label = "C4H6 + O <=> nC4H5 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.5e+06, 'cm^3/(mol*s)'), n=1.9, Ea=(3740, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H6 + O <=> nC4H5 + OH""",
)

entry(
    index = 470,
    label = "C4H6 + O <=> iC4H5 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.5e+06, 'cm^3/(mol*s)'), n=1.9, Ea=(3740, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H6 + O <=> iC4H5 + OH""",
)

entry(
    index = 471,
    label = "C4H6 + O <=> CH3CHCHCO + H",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.5e+08, 'cm^3/(mol*s)'),
        n = 1.45,
        Ea = (-860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C4H6 + O <=> CH3CHCHCO + H""",
)

entry(
    index = 472,
    label = "C4H6 + O <=> CH2CHCHCHO + H",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.5e+08, 'cm^3/(mol*s)'),
        n = 1.45,
        Ea = (-860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C4H6 + O <=> CH2CHCHCHO + H""",
)

entry(
    index = 473,
    label = "C4H6 + OH <=> nC4H5 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.2e+06, 'cm^3/(mol*s)'), n=2, Ea=(3430, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H6 + OH <=> nC4H5 + H2O""",
)

entry(
    index = 474,
    label = "C4H6 + OH <=> iC4H5 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.1e+06, 'cm^3/(mol*s)'), n=2, Ea=(430, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H6 + OH <=> iC4H5 + H2O""",
)

entry(
    index = 475,
    label = "C4H6 + HO2 <=> C4H6O25 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+12, 'cm^3/(mol*s)'), n=0, Ea=(14000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H6 + HO2 <=> C4H6O25 + OH""",
)

entry(
    index = 476,
    label = "C4H6 + HO2 <=> C2H3CHOCH2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(14000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H6 + HO2 <=> C2H3CHOCH2 + OH""",
)

entry(
    index = 477,
    label = "C4H6 + CH3 <=> nC4H5 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+14, 'cm^3/(mol*s)'), n=0, Ea=(22800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H6 + CH3 <=> nC4H5 + CH4""",
)

entry(
    index = 478,
    label = "C4H6 + CH3 <=> iC4H5 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(19800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H6 + CH3 <=> iC4H5 + CH4""",
)

entry(
    index = 479,
    label = "C4H6 + C2H3 <=> nC4H5 + C2H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(22800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H6 + C2H3 <=> nC4H5 + C2H4""",
)

entry(
    index = 480,
    label = "C4H6 + C2H3 <=> iC4H5 + C2H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.5e+13, 'cm^3/(mol*s)'), n=0, Ea=(19800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H6 + C2H3 <=> iC4H5 + C2H4""",
)

entry(
    index = 481,
    label = "C4H6 + C3H3 <=> nC4H5 + aC3H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(22500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H6 + C3H3 <=> nC4H5 + aC3H4""",
)

entry(
    index = 482,
    label = "C4H6 + C3H3 <=> iC4H5 + aC3H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(19500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H6 + C3H3 <=> iC4H5 + aC3H4""",
)

entry(
    index = 483,
    label = "C4H6 + aC3H5 <=> nC4H5 + C3H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(22500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H6 + aC3H5 <=> nC4H5 + C3H6""",
)

entry(
    index = 484,
    label = "C4H6 + aC3H5 <=> iC4H5 + C3H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(19500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H6 + aC3H5 <=> iC4H5 + C3H6""",
)

entry(
    index = 485,
    label = "C4H6 + C2H3 <=> C6H6 + H2 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.62e+11, 'cm^3/(mol*s)'), n=0, Ea=(3240, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H6 + C2H3 <=> C6H6 + H2 + H""",
)

entry(
    index = 486,
    label = "C4H612 <=> iC4H5 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.2e+15, 's^-1'), n=0, Ea=(92600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H612 <=> iC4H5 + H""",
)

entry(
    index = 487,
    label = "C4H612 + H <=> C4H6 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(4000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H612 + H <=> C4H6 + H""",
)

entry(
    index = 488,
    label = "C4H612 + H <=> iC4H5 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(170000, 'cm^3/(mol*s)'), n=2.5, Ea=(2490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H612 + H <=> iC4H5 + H2""",
)

entry(
    index = 489,
    label = "C4H612 + H <=> aC3H4 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(2000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H612 + H <=> aC3H4 + CH3""",
)

entry(
    index = 490,
    label = "C4H612 + H <=> pC3H4 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(2000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H612 + H <=> pC3H4 + CH3""",
)

entry(
    index = 491,
    label = "C4H612 + CH3 <=> iC4H5 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+13, 'cm^3/(mol*s)'), n=0, Ea=(18500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H612 + CH3 <=> iC4H5 + CH4""",
)

entry(
    index = 492,
    label = "C4H612 + O <=> CH2CO + C2H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+08, 'cm^3/(mol*s)'), n=1.65, Ea=(327, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H612 + O <=> CH2CO + C2H4""",
)

entry(
    index = 493,
    label = "C4H612 + O <=> iC4H5 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+11, 'cm^3/(mol*s)'), n=0.7, Ea=(5880, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H612 + O <=> iC4H5 + OH""",
)

entry(
    index = 494,
    label = "C4H612 + OH <=> iC4H5 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.1e+06, 'cm^3/(mol*s)'), n=2, Ea=(-298, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H612 + OH <=> iC4H5 + H2O""",
)

entry(
    index = 495,
    label = "C4H612 <=> C4H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 's^-1'), n=0, Ea=(65000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H612 <=> C4H6""",
)

entry(
    index = 496,
    label = "C4H6-2 <=> C4H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 's^-1'), n=0, Ea=(65000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H6-2 <=> C4H6""",
)

entry(
    index = 497,
    label = "C4H6-2 <=> C4H612",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 's^-1'), n=0, Ea=(67000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H6-2 <=> C4H612""",
)

entry(
    index = 498,
    label = "C4H6-2 + H <=> C4H612 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(4000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H6-2 + H <=> C4H612 + H""",
)

entry(
    index = 499,
    label = "C4H6-2 + H <=> C4H5-2 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(340000, 'cm^3/(mol*s)'), n=2.5, Ea=(2490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H6-2 + H <=> C4H5-2 + H2""",
)

entry(
    index = 500,
    label = "C4H6-2 + H <=> CH3 + pC3H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(260000, 'cm^3/(mol*s)'), n=2.5, Ea=(1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H6-2 + H <=> CH3 + pC3H4""",
)

entry(
    index = 501,
    label = "C4H6-2 <=> H + C4H5-2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+15, 's^-1'), n=0, Ea=(87300, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H6-2 <=> H + C4H5-2""",
)

entry(
    index = 502,
    label = "C4H6-2 + CH3 <=> C4H5-2 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.4e+14, 'cm^3/(mol*s)'), n=0, Ea=(18500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H6-2 + CH3 <=> C4H5-2 + CH4""",
)

entry(
    index = 503,
    label = "C2H3CHOCH2 <=> C4H6O23",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+14, 's^-1'), n=0, Ea=(50600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3CHOCH2 <=> C4H6O23""",
)

entry(
    index = 504,
    label = "C4H6O23 <=> CH3CHCHCHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.95e+13, 's^-1'), n=0, Ea=(49400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H6O23 <=> CH3CHCHCHO""",
)

entry(
    index = 505,
    label = "C4H6O23 <=> C2H4 + CH2CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.75e+15, 's^-1'), n=0, Ea=(69300, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H6O23 <=> C2H4 + CH2CO""",
)

entry(
    index = 506,
    label = "C4H6O23 <=> C2H2 + CH2OCH2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(75800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H6O23 <=> C2H2 + CH2OCH2""",
)

entry(
    index = 507,
    label = "C4H6O25 <=> C4H4O + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.3e+12, 's^-1'), n=0, Ea=(48500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H6O25 <=> C4H4O + H2""",
)

entry(
    index = 508,
    label = "C4H4O <=> CO + pC3H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.78e+15, 's^-1'), n=0, Ea=(77500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H4O <=> CO + pC3H4""",
)

entry(
    index = 509,
    label = "C4H4O <=> C2H2 + CH2CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.01e+14, 's^-1'), n=0, Ea=(77500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H4O <=> C2H2 + CH2CO""",
)

entry(
    index = 510,
    label = "CH3CHCHCHO <=> C3H6 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.9e+14, 's^-1'), n=0, Ea=(69000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHCHCHO <=> C3H6 + CO""",
)

entry(
    index = 511,
    label = "CH3CHCHCHO + H <=> CH2CHCHCHO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(170000, 'cm^3/(mol*s)'), n=2.5, Ea=(2490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHCHCHO + H <=> CH2CHCHCHO + H2""",
)

entry(
    index = 512,
    label = "CH3CHCHCHO + H <=> CH3CHCHCO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(100000, 'cm^3/(mol*s)'), n=2.5, Ea=(2490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHCHCHO + H <=> CH3CHCHCO + H2""",
)

entry(
    index = 513,
    label = "CH3CHCHCHO + H <=> CH3 + C2H3CHO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4e+21, 'cm^3/(mol*s)'),
        n = -2.39,
        Ea = (11180, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3CHCHCHO + H <=> CH3 + C2H3CHO""",
)

entry(
    index = 514,
    label = "CH3CHCHCHO + H <=> C3H6 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4e+21, 'cm^3/(mol*s)'),
        n = -2.39,
        Ea = (11180, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3CHCHCHO + H <=> C3H6 + HCO""",
)

entry(
    index = 515,
    label = "CH3CHCHCHO + CH3 <=> CH2CHCHCHO + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.1, 'cm^3/(mol*s)'), n=3.5, Ea=(5675, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHCHCHO + CH3 <=> CH2CHCHCHO + CH4""",
)

entry(
    index = 516,
    label = "CH3CHCHCHO + CH3 <=> CH3CHCHCO + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.1, 'cm^3/(mol*s)'), n=3.5, Ea=(5675, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHCHCHO + CH3 <=> CH3CHCHCO + CH4""",
)

entry(
    index = 517,
    label = "CH3CHCHCHO + C2H3 <=> CH2CHCHCHO + C2H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.21, 'cm^3/(mol*s)'), n=3.5, Ea=(4682, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHCHCHO + C2H3 <=> CH2CHCHCHO + C2H4""",
)

entry(
    index = 518,
    label = "CH3CHCHCHO + C2H3 <=> CH3CHCHCO + C2H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.11, 'cm^3/(mol*s)'), n=3.5, Ea=(4682, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHCHCHO + C2H3 <=> CH3CHCHCO + C2H4""",
)

entry(
    index = 519,
    label = "CH3CHCHCO <=> CH3CHCH + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 's^-1'), n=0, Ea=(30000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHCHCO <=> CH3CHCH + CO""",
)

entry(
    index = 520,
    label = "CH3CHCHCO + H <=> CH3CHCHCHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHCHCO + H <=> CH3CHCHCHO""",
)

entry(
    index = 521,
    label = "CH2CHCHCHO <=> aC3H5 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 's^-1'), n=0, Ea=(25000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHCHCHO <=> aC3H5 + CO""",
)

entry(
    index = 522,
    label = "CH2CHCHCHO + H <=> CH3CHCHCHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHCHCHO + H <=> CH3CHCHCHO""",
)

entry(
    index = 523,
    label = "C4H7 <=> C4H6 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.48e+53, 's^-1'), n=-12.3, Ea=(52000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H7 <=> C4H6 + H""",
)

entry(
    index = 524,
    label = "C4H7 + H <=> C4H81",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.01e+48, 'cm^6/(mol^2*s)'),
            n = -9.32,
            Ea = (5833.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.498,
        T3 = (1314, 'K'),
        T1 = (1314, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C4H7 + H <=> C4H81""",
)

entry(
    index = 525,
    label = "C4H7 + H <=> CH3 + aC3H5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+21, 'cm^3/(mol*s)'), n=-2, Ea=(11000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H7 + H <=> CH3 + aC3H5""",
)

entry(
    index = 526,
    label = "C4H7 + H <=> C4H6 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H7 + H <=> C4H6 + H2""",
)

entry(
    index = 527,
    label = "C4H7 + O2 <=> C4H6 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H7 + O2 <=> C4H6 + HO2""",
)

entry(
    index = 528,
    label = "C4H7 + HO2 <=> CH2O + OH + aC3H5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H7 + HO2 <=> CH2O + OH + aC3H5""",
)

entry(
    index = 529,
    label = "C4H7 + HCO <=> C4H81 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H7 + HCO <=> C4H81 + CO""",
)

entry(
    index = 530,
    label = "C4H7 + CH3 <=> C4H6 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H7 + CH3 <=> C4H6 + CH4""",
)

entry(
    index = 531,
    label = "iC4H7 + H <=> iC4H8",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.33e+60, 'cm^6/(mol^2*s)'),
            n = -12,
            Ea = (5967.8, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.02,
        T3 = (1096.6, 'K'),
        T1 = (1096.6, 'K'),
        T2 = (6859.5, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is iC4H7 + H <=> iC4H8""",
)

entry(
    index = 532,
    label = "iC4H7 + H <=> CH3CCH2 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.6e+45, 'cm^3/(mol*s)'),
        n = -8.19,
        Ea = (37890, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is iC4H7 + H <=> CH3CCH2 + CH3""",
)

entry(
    index = 533,
    label = "iC4H7 + O <=> CH2O + CH3CCH2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is iC4H7 + O <=> CH2O + CH3CCH2""",
)

entry(
    index = 534,
    label = "iC4H7 + HO2 <=> CH3CCH2 + CH2O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is iC4H7 + HO2 <=> CH3CCH2 + CH2O + OH""",
)

entry(
    index = 535,
    label = "C4H81 + H <=> pC4H9",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(
            A = (1.33e+13, 'cm^3/(mol*s)'),
            n = 0,
            Ea = (3260.7, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        arrheniusLow = Arrhenius(
            A = (6.26e+38, 'cm^6/(mol^2*s)'),
            n = -6.66,
            Ea = (7000, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 1,
        T3 = (1000, 'K'),
        T1 = (1310, 'K'),
        T2 = (48097, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C4H81 + H <=> pC4H9""",
)

entry(
    index = 536,
    label = "C4H81 + H <=> sC4H9",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(
            A = (1.33e+13, 'cm^3/(mol*s)'),
            n = 0,
            Ea = (1559.8, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        arrheniusLow = Arrhenius(
            A = (8.7e+42, 'cm^6/(mol^2*s)'),
            n = -7.5,
            Ea = (4721.8, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 1,
        T3 = (1000, 'K'),
        T1 = (645.4, 'K'),
        T2 = (6844.3, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C4H81 + H <=> sC4H9""",
)

entry(
    index = 537,
    label = "C4H81 + H <=> C2H4 + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.6e+22, 'cm^3/(mol*s)'),
        n = -2.39,
        Ea = (11180, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C4H81 + H <=> C2H4 + C2H5""",
)

entry(
    index = 538,
    label = "C4H81 + H <=> C3H6 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.2e+22, 'cm^3/(mol*s)'),
        n = -2.39,
        Ea = (11180, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C4H81 + H <=> C3H6 + CH3""",
)

entry(
    index = 539,
    label = "C4H81 + H <=> C4H7 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(650000, 'cm^3/(mol*s)'), n=2.54, Ea=(6756, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H81 + H <=> C4H7 + H2""",
)

entry(
    index = 540,
    label = "C4H81 + O <=> nC3H7 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.3e+08, 'cm^3/(mol*s)'),
        n = 1.45,
        Ea = (-402, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C4H81 + O <=> nC3H7 + HCO""",
)

entry(
    index = 541,
    label = "C4H81 + O <=> C4H7 + OH",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(1.5e+13, 'cm^3/(mol*s)'), n=0, Ea=(5760, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(4470, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C4H81 + O <=> C4H7 + OH""",
)

entry(
    index = 542,
    label = "C4H81 + OH <=> C4H7 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(700, 'cm^3/(mol*s)'), n=2.66, Ea=(527, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H81 + OH <=> C4H7 + H2O""",
)

entry(
    index = 543,
    label = "C4H81 + O2 <=> C4H7 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(50930, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H81 + O2 <=> C4H7 + HO2""",
)

entry(
    index = 544,
    label = "C4H81 + HO2 <=> C4H7 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(14340, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H81 + HO2 <=> C4H7 + H2O2""",
)

entry(
    index = 545,
    label = "C4H81 + CH3 <=> C4H7 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.45, 'cm^3/(mol*s)'), n=3.65, Ea=(7153, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H81 + CH3 <=> C4H7 + CH4""",
)

entry(
    index = 546,
    label = "C4H82 + H <=> sC4H9",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(
            A = (1.33e+13, 'cm^3/(mol*s)'),
            n = 0,
            Ea = (1559.8, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        arrheniusLow = Arrhenius(
            A = (8.7e+42, 'cm^6/(mol^2*s)'),
            n = -7.5,
            Ea = (4721.8, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 1,
        T3 = (1000, 'K'),
        T1 = (645.4, 'K'),
        T2 = (6844.3, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C4H82 + H <=> sC4H9""",
)

entry(
    index = 547,
    label = "C4H82 + H <=> C4H7 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(340000, 'cm^3/(mol*s)'), n=2.5, Ea=(2490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H82 + H <=> C4H7 + H2""",
)

entry(
    index = 548,
    label = "C4H82 + O <=> C2H4 + CH3CHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+08, 'cm^3/(mol*s)'), n=1.65, Ea=(327, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H82 + O <=> C2H4 + CH3CHO""",
)

entry(
    index = 549,
    label = "C4H82 + OH <=> C4H7 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.2e+06, 'cm^3/(mol*s)'), n=2, Ea=(-298, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H82 + OH <=> C4H7 + H2O""",
)

entry(
    index = 550,
    label = "C4H82 + O2 <=> C4H7 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(53300, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H82 + O2 <=> C4H7 + HO2""",
)

entry(
    index = 551,
    label = "C4H82 + HO2 <=> C4H7 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(14200, 'cm^3/(mol*s)'), n=2.77, Ea=(14913, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H82 + HO2 <=> C4H7 + H2O2""",
)

entry(
    index = 552,
    label = "C4H82 + CH3 <=> C4H7 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.4, 'cm^3/(mol*s)'), n=3.5, Ea=(5675, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H82 + CH3 <=> C4H7 + CH4""",
)

entry(
    index = 553,
    label = "iC4H8 + H <=> iC4H9",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(
            A = (1.33e+13, 'cm^3/(mol*s)'),
            n = 0,
            Ea = (3260.7, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        arrheniusLow = Arrhenius(
            A = (6.26e+38, 'cm^6/(mol^2*s)'),
            n = -6.66,
            Ea = (7000, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 1,
        T3 = (1000, 'K'),
        T1 = (1310, 'K'),
        T2 = (48097, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is iC4H8 + H <=> iC4H9""",
)

entry(
    index = 554,
    label = "iC4H8 + H <=> iC4H7 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.2e+06, 'cm^3/(mol*s)'),
        n = 2.54,
        Ea = (6760, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is iC4H8 + H <=> iC4H7 + H2""",
)

entry(
    index = 555,
    label = "iC4H8 + H <=> C3H6 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8e+21, 'cm^3/(mol*s)'),
        n = -2.39,
        Ea = (11180, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is iC4H8 + H <=> C3H6 + CH3""",
)

entry(
    index = 556,
    label = "iC4H8 + O <=> CH3 + CH3 + CH2CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+08, 'cm^3/(mol*s)'), n=1.65, Ea=(327, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is iC4H8 + O <=> CH3 + CH3 + CH2CO""",
)

entry(
    index = 557,
    label = "iC4H8 + O <=> iC3H7 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.5e+07, 'cm^3/(mol*s)'),
        n = 1.65,
        Ea = (-972, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is iC4H8 + O <=> iC3H7 + HCO""",
)

entry(
    index = 558,
    label = "iC4H8 + O <=> iC4H7 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(290000, 'cm^3/(mol*s)'), n=2.5, Ea=(3640, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is iC4H8 + O <=> iC4H7 + OH""",
)

entry(
    index = 559,
    label = "iC4H8 + OH <=> iC4H7 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.5e+08, 'cm^3/(mol*s)'), n=1.53, Ea=(775, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is iC4H8 + OH <=> iC4H7 + H2O""",
)

entry(
    index = 560,
    label = "iC4H8 + HO2 <=> iC4H7 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(20000, 'cm^3/(mol*s)'), n=2.55, Ea=(15500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is iC4H8 + HO2 <=> iC4H7 + H2O2""",
)

entry(
    index = 561,
    label = "iC4H8 + O2 <=> iC4H7 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.7e+13, 'cm^3/(mol*s)'), n=0, Ea=(50900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is iC4H8 + O2 <=> iC4H7 + HO2""",
)

entry(
    index = 562,
    label = "iC4H8 + CH3 <=> iC4H7 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.91, 'cm^3/(mol*s)'), n=3.65, Ea=(7150, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is iC4H8 + CH3 <=> iC4H7 + CH4""",
)

entry(
    index = 563,
    label = "C2H4 + C2H5 <=> pC4H9",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.5e+11, 'cm^3/(mol*s)'), n=0, Ea=(7300, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4 + C2H5 <=> pC4H9""",
)

entry(
    index = 564,
    label = "pC4H9 + H <=> C4H10",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.01e+48, 'cm^6/(mol^2*s)'),
            n = -9.32,
            Ea = (5833.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.498,
        T3 = (1314, 'K'),
        T1 = (1314, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is pC4H9 + H <=> C4H10""",
)

entry(
    index = 565,
    label = "pC4H9 + H <=> C2H5 + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.7e+24, 'cm^3/(mol*s)'),
        n = -2.92,
        Ea = (12505, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is pC4H9 + H <=> C2H5 + C2H5""",
)

entry(
    index = 566,
    label = "pC4H9 + H <=> C4H81 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is pC4H9 + H <=> C4H81 + H2""",
)

entry(
    index = 567,
    label = "pC4H9 + O <=> nC3H7 + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is pC4H9 + O <=> nC3H7 + CH2O""",
)

entry(
    index = 568,
    label = "pC4H9 + OH <=> C4H81 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is pC4H9 + OH <=> C4H81 + H2O""",
)

entry(
    index = 569,
    label = "pC4H9 + O2 <=> C4H81 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.7e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is pC4H9 + O2 <=> C4H81 + HO2""",
)

entry(
    index = 570,
    label = "pC4H9 + HO2 <=> nC3H7 + OH + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is pC4H9 + HO2 <=> nC3H7 + OH + CH2O""",
)

entry(
    index = 571,
    label = "pC4H9 + HCO <=> C4H10 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is pC4H9 + HCO <=> C4H10 + CO""",
)

entry(
    index = 572,
    label = "pC4H9 + CH3 <=> C4H81 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is pC4H9 + CH3 <=> C4H81 + CH4""",
)

entry(
    index = 573,
    label = "C3H6 + CH3 <=> sC4H9",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.7e+11, 'cm^3/(mol*s)'), n=0, Ea=(7403.6, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (2.31e+28, 'cm^6/(mol^2*s)'),
            n = -4.27,
            Ea = (1831, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.565,
        T3 = (60000, 'K'),
        T1 = (534.2, 'K'),
        T2 = (3007.2, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C3H6 + CH3 <=> sC4H9""",
)

entry(
    index = 574,
    label = "sC4H9 + H <=> C4H10",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.7e+58, 'cm^6/(mol^2*s)'),
            n = -12.08,
            Ea = (11263.7, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.649,
        T3 = (1213.1, 'K'),
        T1 = (1213.1, 'K'),
        T2 = (13369.7, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is sC4H9 + H <=> C4H10""",
)

entry(
    index = 575,
    label = "sC4H9 + H <=> C2H5 + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.4e+28, 'cm^3/(mol*s)'),
        n = -3.94,
        Ea = (15916, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is sC4H9 + H <=> C2H5 + C2H5""",
)

entry(
    index = 576,
    label = "sC4H9 + H <=> C4H81 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.2e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is sC4H9 + H <=> C4H81 + H2""",
)

entry(
    index = 577,
    label = "sC4H9 + H <=> C4H82 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.1e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is sC4H9 + H <=> C4H82 + H2""",
)

entry(
    index = 578,
    label = "sC4H9 + O <=> CH3CHO + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is sC4H9 + O <=> CH3CHO + C2H5""",
)

entry(
    index = 579,
    label = "sC4H9 + OH <=> C4H81 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is sC4H9 + OH <=> C4H81 + H2O""",
)

entry(
    index = 580,
    label = "sC4H9 + OH <=> C4H82 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is sC4H9 + OH <=> C4H82 + H2O""",
)

entry(
    index = 581,
    label = "sC4H9 + O2 <=> C4H81 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.1e+10, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is sC4H9 + O2 <=> C4H81 + HO2""",
)

entry(
    index = 582,
    label = "sC4H9 + O2 <=> C4H82 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is sC4H9 + O2 <=> C4H82 + HO2""",
)

entry(
    index = 583,
    label = "sC4H9 + HO2 <=> CH3CHO + C2H5 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is sC4H9 + HO2 <=> CH3CHO + C2H5 + OH""",
)

entry(
    index = 584,
    label = "sC4H9 + HCO <=> C4H10 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is sC4H9 + HCO <=> C4H10 + CO""",
)

entry(
    index = 585,
    label = "sC4H9 + CH3 <=> CH4 + C4H81",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.2e+14, 'cm^3/(mol*s)'), n=-0.68, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is sC4H9 + CH3 <=> CH4 + C4H81""",
)

entry(
    index = 586,
    label = "sC4H9 + CH3 <=> CH4 + C4H82",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.5e+14, 'cm^3/(mol*s)'), n=-0.68, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is sC4H9 + CH3 <=> CH4 + C4H82""",
)

entry(
    index = 587,
    label = "C3H6 + CH3 <=> iC4H9",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(9.6e+10, 'cm^3/(mol*s)'), n=0, Ea=(8003.6, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.3e+28, 'cm^6/(mol^2*s)'),
            n = -4.27,
            Ea = (2431.1, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.565,
        T3 = (60000, 'K'),
        T1 = (534.2, 'K'),
        T2 = (3007.2, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C3H6 + CH3 <=> iC4H9""",
)

entry(
    index = 588,
    label = "iC4H9 + H <=> iC4H10",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.27e+56, 'cm^6/(mol^2*s)'),
            n = -11.74,
            Ea = (6430.8, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.506,
        T3 = (1266.6, 'K'),
        T1 = (1266.6, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is iC4H9 + H <=> iC4H10""",
)

entry(
    index = 589,
    label = "iC4H9 + H <=> iC3H7 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.9e+35, 'cm^3/(mol*s)'),
        n = -5.83,
        Ea = (22470, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is iC4H9 + H <=> iC3H7 + CH3""",
)

entry(
    index = 590,
    label = "iC4H9 + H <=> iC4H8 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is iC4H9 + H <=> iC4H8 + H2""",
)

entry(
    index = 591,
    label = "iC4H9 + O <=> iC3H7 + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is iC4H9 + O <=> iC3H7 + CH2O""",
)

entry(
    index = 592,
    label = "iC4H9 + OH <=> iC4H8 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is iC4H9 + OH <=> iC4H8 + H2O""",
)

entry(
    index = 593,
    label = "iC4H9 + O2 <=> iC4H8 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+10, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is iC4H9 + O2 <=> iC4H8 + HO2""",
)

entry(
    index = 594,
    label = "iC4H9 + HO2 <=> iC3H7 + CH2O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.41e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is iC4H9 + HO2 <=> iC3H7 + CH2O + OH""",
)

entry(
    index = 595,
    label = "iC4H9 + HCO <=> iC4H10 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is iC4H9 + HCO <=> iC4H10 + CO""",
)

entry(
    index = 596,
    label = "iC4H9 + CH3 <=> iC4H8 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+12, 'cm^3/(mol*s)'), n=-0.32, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is iC4H9 + CH3 <=> iC4H8 + CH4""",
)

entry(
    index = 597,
    label = "tC4H9 <=> iC4H8 + H",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(8.3e+13, 's^-1'), n=0, Ea=(38150.4, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.9e+41, 'cm^3/(mol*s)'),
            n = -7.36,
            Ea = (36631.7, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.293,
        T3 = (649, 'K'),
        T1 = (60000, 'K'),
        T2 = (3425.9, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is tC4H9 <=> iC4H8 + H""",
)

entry(
    index = 598,
    label = "tC4H9 + H <=> iC4H10",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.47e+61, 'cm^6/(mol^2*s)'),
            n = -12.94,
            Ea = (8000, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        T3 = (1456.4, 'K'),
        T1 = (1000, 'K'),
        T2 = (10000.5, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is tC4H9 + H <=> iC4H10""",
)

entry(
    index = 599,
    label = "tC4H9 + H <=> iC4H8 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.42e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is tC4H9 + H <=> iC4H8 + H2""",
)

entry(
    index = 600,
    label = "tC4H9 + O <=> iC4H8 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is tC4H9 + O <=> iC4H8 + OH""",
)

entry(
    index = 601,
    label = "tC4H9 + O <=> CH3COCH3 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is tC4H9 + O <=> CH3COCH3 + CH3""",
)

entry(
    index = 602,
    label = "tC4H9 + OH <=> iC4H8 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is tC4H9 + OH <=> iC4H8 + H2O""",
)

entry(
    index = 603,
    label = "tC4H9 + O2 <=> iC4H8 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.8e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is tC4H9 + O2 <=> iC4H8 + HO2""",
)

entry(
    index = 604,
    label = "tC4H9 + HO2 <=> CH3 + CH3COCH3 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is tC4H9 + HO2 <=> CH3 + CH3COCH3 + OH""",
)

entry(
    index = 605,
    label = "tC4H9 + HCO <=> iC4H10 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is tC4H9 + HCO <=> iC4H10 + CO""",
)

entry(
    index = 606,
    label = "tC4H9 + CH3 <=> iC4H8 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.8e+15, 'cm^3/(mol*s)'), n=-1, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is tC4H9 + CH3 <=> iC4H8 + CH4""",
)

entry(
    index = 607,
    label = "CH3COCH3 + H <=> H2 + CH2CO + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.3e+06, 'cm^3/(mol*s)'),
        n = 2.54,
        Ea = (6756, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3COCH3 + H <=> H2 + CH2CO + CH3""",
)

entry(
    index = 608,
    label = "CH3COCH3 + O <=> OH + CH2CO + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(190000, 'cm^3/(mol*s)'), n=2.68, Ea=(3716, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3COCH3 + O <=> OH + CH2CO + CH3""",
)

entry(
    index = 609,
    label = "CH3COCH3 + OH <=> H2O + CH2CO + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.2e+07, 'cm^3/(mol*s)'), n=1.8, Ea=(934, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3COCH3 + OH <=> H2O + CH2CO + CH3""",
)

entry(
    index = 610,
    label = "CH3 + CH3CO <=> CH3COCH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+15, 'cm^3/(mol*s)'), n=-0.8, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + CH3CO <=> CH3COCH3""",
)

entry(
    index = 611,
    label = "nC3H7 + CH3 <=> C4H10",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.93e+14, 'cm^3/(mol*s)'), n=-0.32, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (2.68e+61, 'cm^6/(mol^2*s)'),
            n = -13.24,
            Ea = (6000, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 1,
        T3 = (1000, 'K'),
        T1 = (1433.9, 'K'),
        T2 = (5328.8, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is nC3H7 + CH3 <=> C4H10""",
)

entry(
    index = 612,
    label = "C2H5 + C2H5 <=> C4H10",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.88e+14, 'cm^3/(mol*s)'), n=-0.5, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (2.61e+61, 'cm^6/(mol^2*s)'),
            n = -13.42,
            Ea = (6000, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 1,
        T3 = (1000, 'K'),
        T1 = (1433.9, 'K'),
        T2 = (5328.8, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C2H5 + C2H5 <=> C4H10""",
)

entry(
    index = 613,
    label = "C4H10 + H <=> pC4H9 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(920000, 'cm^3/(mol*s)'), n=2.54, Ea=(6756, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H10 + H <=> pC4H9 + H2""",
)

entry(
    index = 614,
    label = "C4H10 + H <=> sC4H9 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+06, 'cm^3/(mol*s)'), n=2.4, Ea=(4471, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H10 + H <=> sC4H9 + H2""",
)

entry(
    index = 615,
    label = "C4H10 + O <=> pC4H9 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.9e+06, 'cm^3/(mol*s)'), n=2.4, Ea=(5500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H10 + O <=> pC4H9 + OH""",
)

entry(
    index = 616,
    label = "C4H10 + O <=> sC4H9 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(430000, 'cm^3/(mol*s)'), n=2.6, Ea=(2580, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H10 + O <=> sC4H9 + OH""",
)

entry(
    index = 617,
    label = "C4H10 + OH <=> pC4H9 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.3e+07, 'cm^3/(mol*s)'), n=1.8, Ea=(954, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H10 + OH <=> pC4H9 + H2O""",
)

entry(
    index = 618,
    label = "C4H10 + OH <=> sC4H9 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.4e+06, 'cm^3/(mol*s)'), n=2, Ea=(-596, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H10 + OH <=> sC4H9 + H2O""",
)

entry(
    index = 619,
    label = "C4H10 + O2 <=> pC4H9 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(50930, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H10 + O2 <=> pC4H9 + HO2""",
)

entry(
    index = 620,
    label = "C4H10 + O2 <=> sC4H9 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+13, 'cm^3/(mol*s)'), n=0, Ea=(47590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H10 + O2 <=> sC4H9 + HO2""",
)

entry(
    index = 621,
    label = "C4H10 + HO2 <=> pC4H9 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(61100, 'cm^3/(mol*s)'), n=2.65, Ea=(17496, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H10 + HO2 <=> pC4H9 + H2O2""",
)

entry(
    index = 622,
    label = "C4H10 + HO2 <=> sC4H9 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(14200, 'cm^3/(mol*s)'), n=2.77, Ea=(14913, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H10 + HO2 <=> sC4H9 + H2O2""",
)

entry(
    index = 623,
    label = "C4H10 + CH3 <=> pC4H9 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.903, 'cm^3/(mol*s)'), n=3.65, Ea=(7153, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H10 + CH3 <=> pC4H9 + CH4""",
)

entry(
    index = 624,
    label = "C4H10 + CH3 <=> sC4H9 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3, 'cm^3/(mol*s)'), n=3.46, Ea=(5480, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H10 + CH3 <=> sC4H9 + CH4""",
)

entry(
    index = 625,
    label = "iC3H7 + CH3 <=> iC4H10",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.4e+15, 'cm^3/(mol*s)'), n=-0.68, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4.16e+61, 'cm^6/(mol^2*s)'),
            n = -13.33,
            Ea = (3903.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.931,
        T3 = (60000, 'K'),
        T1 = (1265.3, 'K'),
        T2 = (5469.8, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is iC3H7 + CH3 <=> iC4H10""",
)

entry(
    index = 626,
    label = "iC4H10 + H <=> iC4H9 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.8e+06, 'cm^3/(mol*s)'),
        n = 2.54,
        Ea = (6760, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is iC4H10 + H <=> iC4H9 + H2""",
)

entry(
    index = 627,
    label = "iC4H10 + H <=> tC4H9 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(600000, 'cm^3/(mol*s)'), n=2.4, Ea=(2580, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is iC4H10 + H <=> tC4H9 + H2""",
)

entry(
    index = 628,
    label = "iC4H10 + O <=> iC4H9 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(430000, 'cm^3/(mol*s)'), n=2.5, Ea=(3640, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is iC4H10 + O <=> iC4H9 + OH""",
)

entry(
    index = 629,
    label = "iC4H10 + O <=> tC4H9 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(157000, 'cm^3/(mol*s)'), n=2.5, Ea=(1110, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is iC4H10 + O <=> tC4H9 + OH""",
)

entry(
    index = 630,
    label = "iC4H10 + OH <=> iC4H9 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.3e+08, 'cm^3/(mol*s)'), n=1.53, Ea=(775, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is iC4H10 + OH <=> iC4H9 + H2O""",
)

entry(
    index = 631,
    label = "iC4H10 + OH <=> tC4H9 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.73e+10, 'cm^3/(mol*s)'), n=0.51, Ea=(64, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is iC4H10 + OH <=> tC4H9 + H2O""",
)

entry(
    index = 632,
    label = "iC4H10 + HO2 <=> iC4H9 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(30000, 'cm^3/(mol*s)'), n=2.55, Ea=(15500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is iC4H10 + HO2 <=> iC4H9 + H2O2""",
)

entry(
    index = 633,
    label = "iC4H10 + HO2 <=> tC4H9 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1240, 'cm^3/(mol*s)'), n=2.77, Ea=(10500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is iC4H10 + HO2 <=> tC4H9 + H2O2""",
)

entry(
    index = 634,
    label = "iC4H10 + O2 <=> iC4H9 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(50900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is iC4H10 + O2 <=> iC4H9 + HO2""",
)

entry(
    index = 635,
    label = "iC4H10 + O2 <=> tC4H9 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(44000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is iC4H10 + O2 <=> tC4H9 + HO2""",
)

entry(
    index = 636,
    label = "iC4H10 + CH3 <=> iC4H9 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.36, 'cm^3/(mol*s)'), n=3.65, Ea=(7150, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is iC4H10 + CH3 <=> iC4H9 + CH4""",
)

entry(
    index = 637,
    label = "iC4H10 + CH3 <=> tC4H9 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.9, 'cm^3/(mol*s)'), n=3.46, Ea=(4600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is iC4H10 + CH3 <=> tC4H9 + CH4""",
)

entry(
    index = 638,
    label = "C6H2 + H <=> C6H3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.1e+30, 'cm^3/(mol*s)'),
        n = -4.92,
        Ea = (10800, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C6H2 + H <=> C6H3""",
)

entry(
    index = 639,
    label = "C6H3 + H <=> C4H2 + C2H2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.8e+23, 'cm^3/(mol*s)'),
        n = -2.55,
        Ea = (10780, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C6H3 + H <=> C4H2 + C2H2""",
)

entry(
    index = 640,
    label = "C6H3 + H <=> l-C6H4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.4e+43, 'cm^3/(mol*s)'),
        n = -9.01,
        Ea = (12120, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C6H3 + H <=> l-C6H4""",
)

entry(
    index = 641,
    label = "C6H3 + H <=> C6H2 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H3 + H <=> C6H2 + H2""",
)

entry(
    index = 642,
    label = "C6H3 + OH <=> C6H2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H3 + OH <=> C6H2 + H2O""",
)

entry(
    index = 643,
    label = "l-C6H4 + H <=> C6H5",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.7e+78, 'cm^3/(mol*s)'),
        n = -19.72,
        Ea = (31400, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is l-C6H4 + H <=> C6H5""",
)

entry(
    index = 644,
    label = "l-C6H4 + H <=> o-C6H4 + H",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.4e+54, 'cm^3/(mol*s)'),
        n = -11.7,
        Ea = (34500, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is l-C6H4 + H <=> o-C6H4 + H""",
)

entry(
    index = 645,
    label = "l-C6H4 + H <=> C6H3 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.33e+06, 'cm^3/(mol*s)'),
        n = 2.53,
        Ea = (9240, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is l-C6H4 + H <=> C6H3 + H2""",
)

entry(
    index = 646,
    label = "l-C6H4 + OH <=> C6H3 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.1e+06, 'cm^3/(mol*s)'), n=2, Ea=(430, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is l-C6H4 + OH <=> C6H3 + H2O""",
)

entry(
    index = 647,
    label = "C4H2 + C2H2 <=> o-C6H4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5e+78, 'cm^3/(mol*s)'),
        n = -19.31,
        Ea = (67920, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C4H2 + C2H2 <=> o-C6H4""",
)

entry(
    index = 648,
    label = "o-C6H4 + OH <=> CO + C5H5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is o-C6H4 + OH <=> CO + C5H5""",
)

entry(
    index = 649,
    label = "C6H5 + CH3 <=> C6H5CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.38e+13, 'cm^3/(mol*s)'), n=0, Ea=(46, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H5 + CH3 <=> C6H5CH3""",
)

entry(
    index = 650,
    label = "C6H5CH3 + O2 <=> C6H5CH2 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+14, 'cm^3/(mol*s)'), n=0, Ea=(42992, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H5CH3 + O2 <=> C6H5CH2 + HO2""",
)

entry(
    index = 651,
    label = "C6H5CH3 + OH <=> C6H5CH2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.62e+13, 'cm^3/(mol*s)'), n=0, Ea=(2770, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H5CH3 + OH <=> C6H5CH2 + H2O""",
)

entry(
    index = 652,
    label = "C6H5CH3 + OH <=> C6H4CH3 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.333e+08, 'cm^3/(mol*s)'),
        n = 1.42,
        Ea = (1450, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C6H5CH3 + OH <=> C6H4CH3 + H2O""",
)

entry(
    index = 653,
    label = "C6H5CH3 + H <=> C6H5CH2 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.259e+14, 'cm^3/(mol*s)'), n=0, Ea=(8359, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H5CH3 + H <=> C6H5CH2 + H2""",
)

entry(
    index = 654,
    label = "C6H5CH3 + H <=> C6H6 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.93e+06, 'cm^3/(mol*s)'),
        n = 2.17,
        Ea = (4163, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C6H5CH3 + H <=> C6H6 + CH3""",
)

entry(
    index = 655,
    label = "C6H5CH3 + O <=> OC6H4CH3 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(3795, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H5CH3 + O <=> OC6H4CH3 + H""",
)

entry(
    index = 656,
    label = "C6H5CH3 + CH3 <=> C6H5CH2 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.16e+11, 'cm^3/(mol*s)'), n=0, Ea=(9500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H5CH3 + CH3 <=> C6H5CH2 + CH4""",
)

entry(
    index = 657,
    label = "C6H5CH3 + C6H5 <=> C6H5CH2 + C6H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.103e+12, 'cm^3/(mol*s)'), n=0, Ea=(4400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H5CH3 + C6H5 <=> C6H5CH2 + C6H6""",
)

entry(
    index = 658,
    label = "C6H5CH3 + HO2 <=> C6H5CH2 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.975e+11, 'cm^3/(mol*s)'),
        n = 0,
        Ea = (14069, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C6H5CH3 + HO2 <=> C6H5CH2 + H2O2""",
)

entry(
    index = 659,
    label = "C6H5CH3 + HO2 <=> C6H4CH3 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.42e+12, 'cm^3/(mol*s)'), n=0, Ea=(28810, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H5CH3 + HO2 <=> C6H4CH3 + H2O2""",
)

entry(
    index = 660,
    label = "C6H5CH2 + H <=> C6H5CH3",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.1e+103, 'cm^6/(mol^2*s)'),
            n = -24.63,
            Ea = (14590, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.431,
        T3 = (383, 'K'),
        T1 = (152, 'K'),
        T2 = (4730, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5},
    ),
    shortDesc = u"""The chemkin file reaction is C6H5CH2 + H <=> C6H5CH3""",
)

entry(
    index = 661,
    label = "C6H5CH2 + H <=> C6H5 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.5e+66, 'cm^3/(mol*s)'),
        n = -13.94,
        Ea = (64580, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C6H5CH2 + H <=> C6H5 + CH3""",
)

entry(
    index = 662,
    label = "C6H5CH2 + O <=> C6H5CHO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H5CH2 + O <=> C6H5CHO + H""",
)

entry(
    index = 663,
    label = "C6H5CH2 + OH <=> C6H5CH2OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H5CH2 + OH <=> C6H5CH2OH""",
)

entry(
    index = 664,
    label = "C6H5CH2 + HO2 <=> C6H5CHO + H + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H5CH2 + HO2 <=> C6H5CHO + H + OH""",
)

entry(
    index = 665,
    label = "C6H5CH2 + C6H5OH <=> C6H5CH3 + C6H5O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.05e+11, 'cm^3/(mol*s)'), n=0, Ea=(9500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H5CH2 + C6H5OH <=> C6H5CH3 + C6H5O""",
)

entry(
    index = 666,
    label = "C6H5CH2 + HOC6H4CH3 <=> C6H5CH3 + OC6H4CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.05e+11, 'cm^3/(mol*s)'), n=0, Ea=(9500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H5CH2 + HOC6H4CH3 <=> C6H5CH3 + OC6H4CH3""",
)

entry(
    index = 667,
    label = "C6H5CH2OH + OH <=> C6H5CHO + H2O + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H5CH2OH + OH <=> C6H5CHO + H2O + H""",
)

entry(
    index = 668,
    label = "C6H5CH2OH + H <=> C6H5CHO + H2 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+13, 'cm^3/(mol*s)'), n=0, Ea=(8235, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H5CH2OH + H <=> C6H5CHO + H2 + H""",
)

entry(
    index = 669,
    label = "C6H5CH2OH + H <=> C6H6 + CH2OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+13, 'cm^3/(mol*s)'), n=0, Ea=(5148, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H5CH2OH + H <=> C6H6 + CH2OH""",
)

entry(
    index = 670,
    label = "C6H5CH2OH + C6H5 <=> C6H5CHO + C6H6 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.4e+12, 'cm^3/(mol*s)'), n=0, Ea=(4400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H5CH2OH + C6H5 <=> C6H5CHO + C6H6 + H""",
)

entry(
    index = 671,
    label = "C6H5 + HCO <=> C6H5CHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H5 + HCO <=> C6H5CHO""",
)

entry(
    index = 672,
    label = "C6H5CHO <=> C6H5CO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.98e+15, 's^-1'), n=0, Ea=(86900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H5CHO <=> C6H5CO + H""",
)

entry(
    index = 673,
    label = "C6H5CHO + O2 <=> C6H5CO + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.02e+13, 'cm^3/(mol*s)'), n=0, Ea=(38950, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H5CHO + O2 <=> C6H5CO + HO2""",
)

entry(
    index = 674,
    label = "C6H5CHO + OH <=> C6H5CO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.35e+10, 'cm^3/(mol*s)'),
        n = 0.73,
        Ea = (-1110, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C6H5CHO + OH <=> C6H5CO + H2O""",
)

entry(
    index = 675,
    label = "C6H5CHO + H <=> C6H5CO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.1e+09, 'cm^3/(mol*s)'),
        n = 1.16,
        Ea = (2400, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C6H5CHO + H <=> C6H5CO + H2""",
)

entry(
    index = 676,
    label = "C6H5CHO + H <=> C6H6 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.93e+06, 'cm^3/(mol*s)'),
        n = 2.17,
        Ea = (4163, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C6H5CHO + H <=> C6H6 + HCO""",
)

entry(
    index = 677,
    label = "C6H5CHO + O <=> C6H5CO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(1800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H5CHO + O <=> C6H5CO + OH""",
)

entry(
    index = 678,
    label = "C6H5CHO + C6H5CH2 <=> C6H5CO + C6H5CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e-06, 'cm^3/(mol*s)'), n=5.6, Ea=(2460, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H5CHO + C6H5CH2 <=> C6H5CO + C6H5CH3""",
)

entry(
    index = 679,
    label = "C6H5CHO + CH3 <=> C6H5CO + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e-06, 'cm^3/(mol*s)'), n=5.6, Ea=(2460, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H5CHO + CH3 <=> C6H5CO + CH4""",
)

entry(
    index = 680,
    label = "C6H5CHO + C6H5 <=> C6H5CO + C6H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.103e+12, 'cm^3/(mol*s)'), n=0, Ea=(4400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H5CHO + C6H5 <=> C6H5CO + C6H6""",
)

entry(
    index = 681,
    label = "C6H5CO + H2O2 <=> C6H5CHO + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+11, 'cm^3/(mol*s)'), n=0, Ea=(8226, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H5CO + H2O2 <=> C6H5CHO + HO2""",
)

entry(
    index = 682,
    label = "OC6H4CH3 + H <=> HOC6H4CH3",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4e+93, 'cm^6/(mol^2*s)'),
            n = -21.84,
            Ea = (13880, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.043,
        T3 = (304.2, 'K'),
        T1 = (60000, 'K'),
        T2 = (5896.4, 'K'),
        efficiencies = {'[C-]#[O+]': 1.5, 'C': 2, '[H][H]': 2, 'O=C=O': 2, 'O': 6},
    ),
    shortDesc = u"""The chemkin file reaction is OC6H4CH3 + H <=> HOC6H4CH3""",
)

entry(
    index = 683,
    label = "OC6H4CH3 + H <=> C6H5O + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.93e+06, 'cm^3/(mol*s)'),
        n = 2.17,
        Ea = (4163, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is OC6H4CH3 + H <=> C6H5O + CH3""",
)

entry(
    index = 684,
    label = "OC6H4CH3 + O <=> C6H4O2 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is OC6H4CH3 + O <=> C6H4O2 + CH3""",
)

entry(
    index = 685,
    label = "HOC6H4CH3 + OH <=> OC6H4CH3 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HOC6H4CH3 + OH <=> OC6H4CH3 + H2O""",
)

entry(
    index = 686,
    label = "HOC6H4CH3 + H <=> OC6H4CH3 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.15e+14, 'cm^3/(mol*s)'), n=0, Ea=(12400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HOC6H4CH3 + H <=> OC6H4CH3 + H2""",
)

entry(
    index = 687,
    label = "HOC6H4CH3 + H <=> C6H5CH3 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.21e+13, 'cm^3/(mol*s)'), n=0, Ea=(7910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HOC6H4CH3 + H <=> C6H5CH3 + OH""",
)

entry(
    index = 688,
    label = "HOC6H4CH3 + H <=> C6H5OH + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+13, 'cm^3/(mol*s)'), n=0, Ea=(5148, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HOC6H4CH3 + H <=> C6H5OH + CH3""",
)

entry(
    index = 689,
    label = "C6H5CO <=> C6H5 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.27e+14, 's^-1'), n=0, Ea=(29013, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H5CO <=> C6H5 + CO""",
)

entry(
    index = 690,
    label = "C6H5 + H <=> C6H6",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (6.6e+75, 'cm^6/(mol^2*s)'),
            n = -16.3,
            Ea = (7000, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 1,
        T3 = (0.1, 'K'),
        T1 = (584.9, 'K'),
        T2 = (6113, 'K'),
        efficiencies = {'[C-]#[O+]': 1.5, 'C': 2, '[H][H]': 2, 'O=C=O': 2, 'O': 6},
    ),
    shortDesc = u"""The chemkin file reaction is C6H5 + H <=> C6H6""",
)

entry(
    index = 691,
    label = "C6H6 + OH <=> C6H5 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (398500, 'cm^3/(mol*s)'),
        n = 2.286,
        Ea = (1058, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C6H6 + OH <=> C6H5 + H2O""",
)

entry(
    index = 692,
    label = "C6H6 + OH <=> C6H5OH + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+13, 'cm^3/(mol*s)'), n=0, Ea=(10600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H6 + OH <=> C6H5OH + H""",
)

entry(
    index = 693,
    label = "C6H6 + O <=> C6H5O + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.39e+13, 'cm^3/(mol*s)'), n=0, Ea=(4910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H6 + O <=> C6H5O + H""",
)

entry(
    index = 694,
    label = "C6H6 + O <=> C5H5 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.39e+13, 'cm^3/(mol*s)'), n=0, Ea=(4530, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H6 + O <=> C5H5 + HCO""",
)

entry(
    index = 695,
    label = "C6H5 + H2 <=> C6H6 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(57070, 'cm^3/(mol*s)'), n=2.43, Ea=(6273, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H5 + H2 <=> C6H6 + H""",
)

entry(
    index = 696,
    label = "C6H5 <=> o-C6H4 + H",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(4.3e+12, 's^-1'), n=0.616, Ea=(77313, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1e+84, 'cm^3/(mol*s)'),
            n = -18.866,
            Ea = (90064, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.902,
        T3 = (696, 'K'),
        T1 = (358, 'K'),
        T2 = (3856, 'K'),
        efficiencies = {'[C-]#[O+]': 1.5, 'C': 2, '[H][H]': 2, 'O=C=O': 2, 'O': 6},
    ),
    shortDesc = u"""The chemkin file reaction is C6H5 <=> o-C6H4 + H""",
)

entry(
    index = 697,
    label = "C6H5 + H <=> o-C6H4 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+11, 'cm^3/(mol*s)'), n=1.1, Ea=(24500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H5 + H <=> o-C6H4 + H2""",
)

entry(
    index = 698,
    label = "C6H5 + O2 <=> C6H5O + O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(6120, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H5 + O2 <=> C6H5O + O""",
)

entry(
    index = 699,
    label = "C6H5 + O2 <=> C6H4O2 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(8980, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H5 + O2 <=> C6H4O2 + H""",
)

entry(
    index = 700,
    label = "C6H5 + O <=> C5H5 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H5 + O <=> C5H5 + CO""",
)

entry(
    index = 701,
    label = "C6H5 + OH <=> C6H5O + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H5 + OH <=> C6H5O + H""",
)

entry(
    index = 702,
    label = "C6H5 + HO2 <=> C6H5O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H5 + HO2 <=> C6H5O + OH""",
)

entry(
    index = 703,
    label = "C6H5 + HO2 <=> C6H6 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H5 + HO2 <=> C6H6 + O2""",
)

entry(
    index = 704,
    label = "C6H5 + CH4 <=> C6H6 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (0.00389, 'cm^3/(mol*s)'),
        n = 4.57,
        Ea = (5256, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C6H5 + CH4 <=> C6H6 + CH3""",
)

entry(
    index = 705,
    label = "C6H5 + C2H6 <=> C6H6 + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.1e+11, 'cm^3/(mol*s)'), n=0, Ea=(4443, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H5 + C2H6 <=> C6H6 + C2H5""",
)

entry(
    index = 706,
    label = "C6H5 + CH2O <=> C6H6 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(85500, 'cm^3/(mol*s)'), n=2.19, Ea=(38, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H5 + CH2O <=> C6H6 + HCO""",
)

entry(
    index = 707,
    label = "C6H4O2 <=> C5H4O + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.4e+11, 's^-1'), n=0, Ea=(59000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H4O2 <=> C5H4O + CO""",
)

entry(
    index = 708,
    label = "C6H4O2 + H <=> CO + C5H5O(1,3)",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.3e+09, 'cm^3/(mol*s)'),
        n = 1.45,
        Ea = (3900, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C6H4O2 + H <=> CO + C5H5O(1,3)""",
)

entry(
    index = 709,
    label = "C6H4O2 + O <=> CO + CO + C2H2 + CH2CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(5000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H4O2 + O <=> CO + CO + C2H2 + CH2CO""",
)

entry(
    index = 710,
    label = "C6H5O + H <=> C5H5 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(12000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H5O + H <=> C5H5 + HCO""",
)

entry(
    index = 711,
    label = "C6H5O + H <=> C5H6 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H5O + H <=> C5H6 + CO""",
)

entry(
    index = 712,
    label = "C6H5O <=> CO + C5H5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.76e+54, 's^-1'), n=-12.06, Ea=(72800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H5O <=> CO + C5H5""",
)

entry(
    index = 713,
    label = "C6H5O + O <=> C6H4O2 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.6e+10, 'cm^3/(mol*s)'), n=0.47, Ea=(795, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H5O + O <=> C6H4O2 + H""",
)

entry(
    index = 714,
    label = "C6H5OH <=> C5H6 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 's^-1'), n=0, Ea=(60808, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H5OH <=> C5H6 + CO""",
)

entry(
    index = 715,
    label = "C6H5OH + OH <=> C6H5O + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.95e+06, 'cm^3/(mol*s)'), n=2, Ea=(-1312, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H5OH + OH <=> C6H5O + H2O""",
)

entry(
    index = 716,
    label = "C6H5OH + H <=> C6H5O + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.15e+14, 'cm^3/(mol*s)'), n=0, Ea=(12398, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H5OH + H <=> C6H5O + H2""",
)

entry(
    index = 717,
    label = "C6H5OH + O <=> C6H5O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.81e+13, 'cm^3/(mol*s)'), n=0, Ea=(7352, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H5OH + O <=> C6H5O + OH""",
)

entry(
    index = 718,
    label = "C6H5OH + C2H3 <=> C6H5O + C2H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H5OH + C2H3 <=> C6H5O + C2H4""",
)

entry(
    index = 719,
    label = "C6H5OH + nC4H5 <=> C6H5O + C4H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H5OH + nC4H5 <=> C6H5O + C4H6""",
)

entry(
    index = 720,
    label = "C6H5OH + C6H5 <=> C6H5O + C6H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.91e+12, 'cm^3/(mol*s)'), n=0, Ea=(4400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H5OH + C6H5 <=> C6H5O + C6H6""",
)

entry(
    index = 721,
    label = "C5H6 + H <=> C2H2 + aC3H5",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (7.74e+36, 'cm^3/(mol*s)'),
        n = -6.18,
        Ea = (32890, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C5H6 + H <=> C2H2 + aC3H5""",
)

entry(
    index = 722,
    label = "C5H6 + H <=> lC5H7",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8.27e+126, 'cm^3/(mol*s)'),
        n = -32.3,
        Ea = (82348, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C5H6 + H <=> lC5H7""",
)

entry(
    index = 723,
    label = "C5H6 + H <=> C5H5 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.03e+08, 'cm^3/(mol*s)'),
        n = 1.71,
        Ea = (5590, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C5H6 + H <=> C5H5 + H2""",
)

entry(
    index = 724,
    label = "C5H6 + O <=> C5H5 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(47700, 'cm^3/(mol*s)'), n=2.71, Ea=(1106, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H6 + O <=> C5H5 + OH""",
)

entry(
    index = 725,
    label = "C5H6 + O <=> C5H5O(1,3) + H",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(
                A = (8.91e+12, 'cm^3/(mol*s)'),
                n = -0.15,
                Ea = (590, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (5.6e+12, 'cm^3/(mol*s)'),
                n = -0.06,
                Ea = (200, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C5H6 + O <=> C5H5O(1,3) + H""",
)

entry(
    index = 726,
    label = "C5H6 + O <=> nC4H5 + CO + H",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8.7e+51, 'cm^3/(mol*s)'),
        n = -11.09,
        Ea = (33240, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C5H6 + O <=> nC4H5 + CO + H""",
)

entry(
    index = 727,
    label = "C5H6 + OH <=> C5H5 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.08e+06, 'cm^3/(mol*s)'), n=2, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H6 + OH <=> C5H5 + H2O""",
)

entry(
    index = 728,
    label = "C5H6 + HO2 <=> C5H5 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(11000, 'cm^3/(mol*s)'), n=2.6, Ea=(12900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H6 + HO2 <=> C5H5 + H2O2""",
)

entry(
    index = 729,
    label = "C5H6 + O2 <=> C5H5 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(37150, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H6 + O2 <=> C5H5 + HO2""",
)

entry(
    index = 730,
    label = "C5H6 + HCO <=> C5H5 + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.08e+08, 'cm^3/(mol*s)'),
        n = 1.9,
        Ea = (16000, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C5H6 + HCO <=> C5H5 + CH2O""",
)

entry(
    index = 731,
    label = "C5H6 + CH3 <=> C5H5 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.18, 'cm^3/(mol*s)'), n=4, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H6 + CH3 <=> C5H5 + CH4""",
)

entry(
    index = 732,
    label = "C5H5 + H <=> C5H6",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4.4e+80, 'cm^6/(mol^2*s)'),
            n = -18.28,
            Ea = (12994, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.068,
        T3 = (400.7, 'K'),
        T1 = (4135.8, 'K'),
        T2 = (5501.9, 'K'),
        efficiencies = {'[C-]#[O+]': 1.5, 'C': 2, '[H][H]': 2, 'O=C=O': 2, 'O': 6},
    ),
    shortDesc = u"""The chemkin file reaction is C5H5 + H <=> C5H6""",
)

entry(
    index = 733,
    label = "C5H5 + O2 <=> C5H5O(2,4) + O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (7.78e+15, 'cm^3/(mol*s)'),
        n = -0.73,
        Ea = (48740, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C5H5 + O2 <=> C5H5O(2,4) + O""",
)

entry(
    index = 734,
    label = "C5H5 + O <=> C5H5O(2,4)",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.12e-12, 'cm^3/(mol*s)'),
        n = 5.87,
        Ea = (-17310, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C5H5 + O <=> C5H5O(2,4)""",
)

entry(
    index = 735,
    label = "C5H5 + O <=> C5H4O + H",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.81e+13, 'cm^3/(mol*s)'),
        n = -0.02,
        Ea = (20, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C5H5 + O <=> C5H4O + H""",
)

entry(
    index = 736,
    label = "C5H5 + O <=> nC4H5 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.2e+13, 'cm^3/(mol*s)'),
        n = -0.17,
        Ea = (440, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C5H5 + O <=> nC4H5 + CO""",
)

entry(
    index = 737,
    label = "C5H5 + OH <=> C5H4OH + H",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.51e+57, 'cm^3/(mol*s)'),
        n = -12.18,
        Ea = (48350, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C5H5 + OH <=> C5H4OH + H""",
)

entry(
    index = 738,
    label = "C5H5 + OH <=> C5H5O(2,4) + H",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.36e+51, 'cm^3/(mol*s)'),
        n = -10.46,
        Ea = (57100, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C5H5 + OH <=> C5H5O(2,4) + H""",
)

entry(
    index = 739,
    label = "C5H5 + HO2 <=> C5H5O(2,4) + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (6.27e+29, 'cm^3/(mol*s)'),
        n = -4.69,
        Ea = (11650, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C5H5 + HO2 <=> C5H5O(2,4) + OH""",
)

entry(
    index = 740,
    label = "C5H5 + OH <=> C5H5OH",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(
                A = (6.49e+14, 'cm^3/(mol*s)'),
                n = -0.85,
                Ea = (-2730, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.15e+43, 'cm^3/(mol*s)'),
                n = -8.76,
                Ea = (18730, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (1.06e+59, 'cm^3/(mol*s)'),
                n = -13.08,
                Ea = (33450, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C5H5 + OH <=> C5H5OH""",
)

entry(
    index = 741,
    label = "C5H5 + O2 <=> C5H4O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.8e+12, 'cm^3/(mol*s)'),
        n = 0.08,
        Ea = (18000, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C5H5 + O2 <=> C5H4O + OH""",
)

entry(
    index = 742,
    label = "C5H5OH + H <=> C5H5O(2,4) + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.15e+14, 'cm^3/(mol*s)'), n=0, Ea=(15400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H5OH + H <=> C5H5O(2,4) + H2""",
)

entry(
    index = 743,
    label = "C5H5OH + H <=> C5H4OH + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(120000, 'cm^3/(mol*s)'), n=2.5, Ea=(1492, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H5OH + H <=> C5H4OH + H2""",
)

entry(
    index = 744,
    label = "C5H5OH + OH <=> C5H5O(2,4) + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H5OH + OH <=> C5H5O(2,4) + H2O""",
)

entry(
    index = 745,
    label = "C5H5OH + OH <=> C5H4OH + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.08e+06, 'cm^3/(mol*s)'), n=2, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H5OH + OH <=> C5H4OH + H2O""",
)

entry(
    index = 746,
    label = "C5H5O(2,4) + H <=> C5H5OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H5O(2,4) + H <=> C5H5OH""",
)

entry(
    index = 747,
    label = "C5H5O(2,4) <=> C5H4O + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 's^-1'), n=0, Ea=(30000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H5O(2,4) <=> C5H4O + H""",
)

entry(
    index = 748,
    label = "C5H5O(2,4) + O2 <=> C5H4O + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H5O(2,4) + O2 <=> C5H4O + HO2""",
)

entry(
    index = 749,
    label = "C5H4O + H <=> C5H5O(1,3)",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(2000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H4O + H <=> C5H5O(1,3)""",
)

entry(
    index = 750,
    label = "C5H5O(1,3) <=> c-C4H5 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 's^-1'), n=0, Ea=(36000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H5O(1,3) <=> c-C4H5 + CO""",
)

entry(
    index = 751,
    label = "C5H5O(1,3) + O2 <=> C5H4O + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H5O(1,3) + O2 <=> C5H4O + HO2""",
)

entry(
    index = 752,
    label = "C5H4OH <=> C5H4O + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.1e+13, 's^-1'), n=0, Ea=(48000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H4OH <=> C5H4O + H""",
)

entry(
    index = 753,
    label = "C5H4O <=> C2H2 + C2H2 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.2e+41, 's^-1'), n=-7.87, Ea=(98700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H4O <=> C2H2 + C2H2 + CO""",
)

entry(
    index = 754,
    label = "C5H4O + H <=> CO + c-C4H5",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.3e+09, 'cm^3/(mol*s)'),
        n = 1.45,
        Ea = (3900, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C5H4O + H <=> CO + c-C4H5""",
)

entry(
    index = 755,
    label = "C5H4O + O <=> CO + HCO + C3H3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (6.2e+08, 'cm^3/(mol*s)'),
        n = 1.45,
        Ea = (-858, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C5H4O + O <=> CO + HCO + C3H3""",
)

entry(
    index = 756,
    label = "c-C4H5 + H <=> C4H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is c-C4H5 + H <=> C4H6""",
)

entry(
    index = 757,
    label = "c-C4H5 + H <=> C2H4 + C2H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is c-C4H5 + H <=> C2H4 + C2H2""",
)

entry(
    index = 758,
    label = "c-C4H5 + O <=> CH2CHO + C2H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is c-C4H5 + O <=> CH2CHO + C2H2""",
)

entry(
    index = 759,
    label = "c-C4H5 + O2 <=> CH2CHO + CH2CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.8e+11, 'cm^3/(mol*s)'), n=0, Ea=(19000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is c-C4H5 + O2 <=> CH2CHO + CH2CO""",
)

entry(
    index = 760,
    label = "c-C4H5 <=> C4H4 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+12, 's^-1'), n=0, Ea=(52000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is c-C4H5 <=> C4H4 + H""",
)

entry(
    index = 761,
    label = "c-C4H5 <=> C2H3 + C2H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+12, 's^-1'), n=0, Ea=(58000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is c-C4H5 <=> C2H3 + C2H2""",
)

entry(
    index = 762,
    label = "aC3H5 + C2H3 <=> lC5H7 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is aC3H5 + C2H3 <=> lC5H7 + H""",
)

entry(
    index = 763,
    label = "lC5H7 + O <=> C2H3CHO + C2H3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is lC5H7 + O <=> C2H3CHO + C2H3""",
)

entry(
    index = 764,
    label = "lC5H7 + OH <=> C2H3CHO + C2H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is lC5H7 + OH <=> C2H3CHO + C2H4""",
)

entry(
    index = 765,
    label = "PXC5H9 <=> C5H8-13 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.48e+53, 's^-1'), n=-12.3, Ea=(52000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC5H9 <=> C5H8-13 + H""",
)

entry(
    index = 766,
    label = "PXC5H9 + H <=> C5H10",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.01e+48, 'cm^6/(mol^2*s)'),
            n = -9.32,
            Ea = (5833.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.498,
        T3 = (1314, 'K'),
        T1 = (1314, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC5H9 + H <=> C5H10""",
)

entry(
    index = 767,
    label = "PXC5H9 + H <=> CH3 + C4H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+21, 'cm^3/(mol*s)'), n=-2, Ea=(11000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC5H9 + H <=> CH3 + C4H7""",
)

entry(
    index = 768,
    label = "PXC5H9 + HO2 <=> CH2O + OH + C4H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC5H9 + HO2 <=> CH2O + OH + C4H7""",
)

entry(
    index = 769,
    label = "PXC5H9 + HCO <=> C5H10 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC5H9 + HCO <=> C5H10 + CO""",
)

entry(
    index = 770,
    label = "PXC5H9 <=> C2H4 + aC3H5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.04e+12, 's^-1'), n=-0.37, Ea=(25124, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC5H9 <=> C2H4 + aC3H5""",
)

entry(
    index = 771,
    label = "cC5H9 <=> PXC5H9",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.15e+12, 's^-1'), n=-0.3, Ea=(33721, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC5H9 <=> PXC5H9""",
)

entry(
    index = 772,
    label = "cC5H9 <=> cC5H8 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.53e+11, 's^-1'), n=-0.55, Ea=(33140, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC5H9 <=> cC5H8 + H""",
)

entry(
    index = 773,
    label = "SXC5H9 <=> C5H8-13 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.5e+08, 's^-1'), n=-1.35, Ea=(32487, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC5H9 <=> C5H8-13 + H""",
)

entry(
    index = 774,
    label = "SXC5H9 <=> C3H6 + C2H3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.41e+11, 's^-1'), n=0.56, Ea=(37213, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC5H9 <=> C3H6 + C2H3""",
)

entry(
    index = 775,
    label = "SXC5H9 <=> C5H8-14 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.69e-09, 's^-1'), n=-1.17, Ea=(37097, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC5H9 <=> C5H8-14 + H""",
)

entry(
    index = 776,
    label = "SXC5H9 <=> PXCH2-3-1C4H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.52e+08, 's^-1'), n=-1.42, Ea=(14609, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC5H9 <=> PXCH2-3-1C4H7""",
)

entry(
    index = 777,
    label = "PXCH2-3-1C4H7 <=> C4H6 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.6e+11, 's^-1'), n=0.41, Ea=(31254, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXCH2-3-1C4H7 <=> C4H6 + CH3""",
)

entry(
    index = 778,
    label = "SAXC5H9 <=> CH3 + C4H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.7e+12, 's^-1'), n=-0.1, Ea=(35891, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SAXC5H9 <=> CH3 + C4H6""",
)

entry(
    index = 779,
    label = "C5H10 <=> C2H5 + aC3H5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.24e+22, 's^-1'), n=-1.94, Ea=(75470, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10 <=> C2H5 + aC3H5""",
)

entry(
    index = 780,
    label = "C5H10 <=> C3H6 + C2H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.62e+06, 's^-1'), n=1.81, Ea=(53454, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10 <=> C3H6 + C2H4""",
)

entry(
    index = 781,
    label = "C5H10 + H <=> C2H4 + nC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8e+21, 'cm^3/(mol*s)'),
        n = -2.39,
        Ea = (11180, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C5H10 + H <=> C2H4 + nC3H7""",
)

entry(
    index = 782,
    label = "C5H10 + H <=> C3H6 + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.6e+22, 'cm^3/(mol*s)'),
        n = -2.39,
        Ea = (11180, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C5H10 + H <=> C3H6 + C2H5""",
)

entry(
    index = 783,
    label = "C5H10 + H <=> SAXC5H9 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(54000, 'cm^3/(mol*s)'), n=2.5, Ea=(-1900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10 + H <=> SAXC5H9 + H2""",
)

entry(
    index = 784,
    label = "C5H10 + H <=> PXC5H9 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0323, 'cm^3/(mol*s)'), n=4.7, Ea=(3679, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10 + H <=> PXC5H9 + H2""",
)

entry(
    index = 785,
    label = "C5H10 + H <=> SXC5H9 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0317, 'cm^3/(mol*s)'), n=4.65, Ea=(1340, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10 + H <=> SXC5H9 + H2""",
)

entry(
    index = 786,
    label = "C5H10 + O <=> pC4H9 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.3e+08, 'cm^3/(mol*s)'),
        n = 1.45,
        Ea = (-402, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C5H10 + O <=> pC4H9 + HCO""",
)

entry(
    index = 787,
    label = "C5H10 + O <=> PXC5H9 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.745, 'cm^3/(mol*s)'), n=4.17, Ea=(2766, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10 + O <=> PXC5H9 + OH""",
)

entry(
    index = 788,
    label = "C5H10 + OH <=> PXC5H9 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(700, 'cm^3/(mol*s)'), n=2.66, Ea=(527, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10 + OH <=> PXC5H9 + H2O""",
)

entry(
    index = 789,
    label = "C5H10 + O2 <=> PXC5H9 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(50930, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10 + O2 <=> PXC5H9 + HO2""",
)

entry(
    index = 790,
    label = "C5H10 + HO2 <=> PXC5H9 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(14340, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10 + HO2 <=> PXC5H9 + H2O2""",
)

entry(
    index = 791,
    label = "C5H10 + CH3 <=> PXC5H9 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.58e-08, 'cm^3/(mol*s)'),
        n = 6.08,
        Ea = (6223, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C5H10 + CH3 <=> PXC5H9 + CH4""",
)

entry(
    index = 792,
    label = "PXC5H11 <=> C2H4 + nC3H7",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1e+13, 's^-1'), n=0, Ea=(28366.4, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (7.1e-35, 'cm^3/(mol*s)'),
            n = 15.411,
            Ea = (-600, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -5.91,
        T3 = (333, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC5H11 <=> C2H4 + nC3H7""",
)

entry(
    index = 793,
    label = "SXC5H11 <=> C3H6 + C2H5",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(8e+12, 's^-1'), n=0, Ea=(27392.8, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.7e-33, 'cm^3/(mol*s)'),
            n = 14.91,
            Ea = (-600, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -6.53,
        T3 = (333, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SXC5H11 <=> C3H6 + C2H5""",
)

entry(
    index = 794,
    label = "S2XC5H11 <=> C4H81 + CH3",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.8e+13, 's^-1'), n=0, Ea=(29348, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4e-39, 'cm^3/(mol*s)'),
            n = 16.782,
            Ea = (-600.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -7.03,
        T3 = (314, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S2XC5H11 <=> C4H81 + CH3""",
)

entry(
    index = 795,
    label = "SXC5H11 + H <=> NC5H12",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.7e+58, 'cm^6/(mol^2*s)'),
            n = -12.08,
            Ea = (11263.7, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.649,
        T3 = (1213.1, 'K'),
        T1 = (1213.1, 'K'),
        T2 = (13369.7, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SXC5H11 + H <=> NC5H12""",
)

entry(
    index = 796,
    label = "SXC5H11 + H <=> nC3H7 + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.4e+28, 'cm^3/(mol*s)'),
        n = -3.94,
        Ea = (15916, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is SXC5H11 + H <=> nC3H7 + C2H5""",
)

entry(
    index = 797,
    label = "SXC5H11 + H <=> C5H10 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.2e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC5H11 + H <=> C5H10 + H2""",
)

entry(
    index = 798,
    label = "SXC5H11 + O <=> CH3CHO + nC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC5H11 + O <=> CH3CHO + nC3H7""",
)

entry(
    index = 799,
    label = "SXC5H11 + OH <=> C5H10 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC5H11 + OH <=> C5H10 + H2O""",
)

entry(
    index = 800,
    label = "SXC5H11 + O2 <=> C5H10 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC5H11 + O2 <=> C5H10 + HO2""",
)

entry(
    index = 801,
    label = "SXC5H11 + HO2 <=> CH3CHO + nC3H7 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC5H11 + HO2 <=> CH3CHO + nC3H7 + OH""",
)

entry(
    index = 802,
    label = "SXC5H11 + HCO <=> NC5H12 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC5H11 + HCO <=> NC5H12 + CO""",
)

entry(
    index = 803,
    label = "SXC5H11 + CH3 <=> CH4 + C5H10",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.2e+14, 'cm^3/(mol*s)'), n=-0.68, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC5H11 + CH3 <=> CH4 + C5H10""",
)

entry(
    index = 804,
    label = "PXC5H11 + H <=> NC5H12",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.01e+48, 'cm^6/(mol^2*s)'),
            n = -9.32,
            Ea = (5833.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.498,
        T3 = (1314, 'K'),
        T1 = (1314, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC5H11 + H <=> NC5H12""",
)

entry(
    index = 805,
    label = "PXC5H11 + H <=> nC3H7 + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.7e+24, 'cm^3/(mol*s)'),
        n = -2.92,
        Ea = (12505, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is PXC5H11 + H <=> nC3H7 + C2H5""",
)

entry(
    index = 806,
    label = "PXC5H11 + H <=> C5H10 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC5H11 + H <=> C5H10 + H2""",
)

entry(
    index = 807,
    label = "PXC5H11 + O <=> pC4H9 + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC5H11 + O <=> pC4H9 + CH2O""",
)

entry(
    index = 808,
    label = "PXC5H11 + OH <=> C5H10 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC5H11 + OH <=> C5H10 + H2O""",
)

entry(
    index = 809,
    label = "PXC5H11 + O2 <=> C5H10 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+10, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC5H11 + O2 <=> C5H10 + HO2""",
)

entry(
    index = 810,
    label = "PXC5H11 + HCO <=> NC5H12 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC5H11 + HCO <=> NC5H12 + CO""",
)

entry(
    index = 811,
    label = "PXC5H11 + CH3 <=> C5H10 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC5H11 + CH3 <=> C5H10 + CH4""",
)

entry(
    index = 812,
    label = "pC4H9 + CH3 <=> NC5H12",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.93e+14, 'cm^3/(mol*s)'), n=-0.32, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is pC4H9 + CH3 <=> NC5H12""",
)

entry(
    index = 813,
    label = "nC3H7 + C2H5 <=> NC5H12",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.88e+14, 'cm^3/(mol*s)'), n=-0.5, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is nC3H7 + C2H5 <=> NC5H12""",
)

entry(
    index = 814,
    label = "NC5H12 + OH <=> PXC5H11 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.73e+07, 'cm^3/(mol*s)'),
        n = 1.81,
        Ea = (868.3, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC5H12 + OH <=> PXC5H11 + H2O""",
)

entry(
    index = 815,
    label = "NC5H12 + OH <=> SXC5H11 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.41e+10, 'cm^3/(mol*s)'),
        n = 0.94,
        Ea = (504.7, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC5H12 + OH <=> SXC5H11 + H2O""",
)

entry(
    index = 816,
    label = "NC5H12 + OH <=> S2XC5H11 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.72e+06, 'cm^3/(mol*s)'),
        n = 1.81,
        Ea = (-1015.4, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC5H12 + OH <=> S2XC5H11 + H2O""",
)

entry(
    index = 817,
    label = "NC5H12 + O2 <=> PXC5H11 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(50930, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H12 + O2 <=> PXC5H11 + HO2""",
)

entry(
    index = 818,
    label = "NC5H12 + O2 <=> SXC5H11 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+13, 'cm^3/(mol*s)'), n=0, Ea=(47590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H12 + O2 <=> SXC5H11 + HO2""",
)

entry(
    index = 819,
    label = "NC5H12 + O2 <=> S2XC5H11 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(47590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H12 + O2 <=> S2XC5H11 + HO2""",
)

entry(
    index = 820,
    label = "NC5H12 + HO2 <=> PXC5H11 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(61100, 'cm^3/(mol*s)'), n=2.65, Ea=(17496, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H12 + HO2 <=> PXC5H11 + H2O2""",
)

entry(
    index = 821,
    label = "NC5H12 + HO2 <=> SXC5H11 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(14200, 'cm^3/(mol*s)'), n=2.77, Ea=(14913, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H12 + HO2 <=> SXC5H11 + H2O2""",
)

entry(
    index = 822,
    label = "NC5H12 + HO2 <=> S2XC5H11 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7130, 'cm^3/(mol*s)'), n=2.77, Ea=(14913, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H12 + HO2 <=> S2XC5H11 + H2O2""",
)

entry(
    index = 823,
    label = "NC5H12 + H <=> PXC5H11 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0645, 'cm^3/(mol*s)'), n=4.7, Ea=(3679, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H12 + H <=> PXC5H11 + H2""",
)

entry(
    index = 824,
    label = "NC5H12 + H <=> SXC5H11 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0634, 'cm^3/(mol*s)'), n=4.65, Ea=(1340, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H12 + H <=> SXC5H11 + H2""",
)

entry(
    index = 825,
    label = "NC5H12 + H <=> S2XC5H11 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0317, 'cm^3/(mol*s)'), n=4.65, Ea=(1340, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H12 + H <=> S2XC5H11 + H2""",
)

entry(
    index = 826,
    label = "NC5H12 + O <=> PXC5H11 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.49, 'cm^3/(mol*s)'), n=4.17, Ea=(2766, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H12 + O <=> PXC5H11 + OH""",
)

entry(
    index = 827,
    label = "NC5H12 + O <=> SXC5H11 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(36.5, 'cm^3/(mol*s)'), n=3.75, Ea=(825, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H12 + O <=> SXC5H11 + OH""",
)

entry(
    index = 828,
    label = "NC5H12 + O <=> S2XC5H11 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(18.2, 'cm^3/(mol*s)'), n=3.75, Ea=(825, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H12 + O <=> S2XC5H11 + OH""",
)

entry(
    index = 829,
    label = "NC5H12 + CH3 <=> PXC5H11 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.16e-08, 'cm^3/(mol*s)'),
        n = 6.08,
        Ea = (6223, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC5H12 + CH3 <=> PXC5H11 + CH4""",
)

entry(
    index = 830,
    label = "NC5H12 + CH3 <=> SXC5H11 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.82e-07, 'cm^3/(mol*s)'),
        n = 5.89,
        Ea = (4768, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC5H12 + CH3 <=> SXC5H11 + CH4""",
)

entry(
    index = 831,
    label = "NC5H12 + CH3 <=> S2XC5H11 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.41e-07, 'cm^3/(mol*s)'),
        n = 5.89,
        Ea = (4768, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC5H12 + CH3 <=> S2XC5H11 + CH4""",
)

entry(
    index = 832,
    label = "PXC6H11 <=> C6H10-13 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.48e+53, 's^-1'), n=-12.3, Ea=(52000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC6H11 <=> C6H10-13 + H""",
)

entry(
    index = 833,
    label = "PXC6H11 + H <=> C6H12",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.01e+48, 'cm^6/(mol^2*s)'),
            n = -9.32,
            Ea = (5833.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.498,
        T3 = (1314, 'K'),
        T1 = (1314, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC6H11 + H <=> C6H12""",
)

entry(
    index = 834,
    label = "PXC6H11 + H <=> CH3 + PXC5H9",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+21, 'cm^3/(mol*s)'), n=-2, Ea=(11000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC6H11 + H <=> CH3 + PXC5H9""",
)

entry(
    index = 835,
    label = "PXC6H11 + H <=> C6H10-13 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC6H11 + H <=> C6H10-13 + H2""",
)

entry(
    index = 836,
    label = "PXC6H11 + O2 <=> C6H10-13 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC6H11 + O2 <=> C6H10-13 + HO2""",
)

entry(
    index = 837,
    label = "PXC6H11 + HO2 <=> CH2O + OH + PXC5H9",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC6H11 + HO2 <=> CH2O + OH + PXC5H9""",
)

entry(
    index = 838,
    label = "PXC6H11 + HCO <=> C6H12 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC6H11 + HCO <=> C6H12 + CO""",
)

entry(
    index = 839,
    label = "PXC6H11 + CH3 <=> C6H10-13 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC6H11 + CH3 <=> C6H10-13 + CH4""",
)

entry(
    index = 840,
    label = "C6H12 <=> aC3H5 + nC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.07e+23, 's^-1'), n=-2.03, Ea=(74958, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12 <=> aC3H5 + nC3H7""",
)

entry(
    index = 841,
    label = "C6H12 <=> C3H6 + C3H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.08e+06, 's^-1'), n=1.65, Ea=(53752, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12 <=> C3H6 + C3H6""",
)

entry(
    index = 842,
    label = "C6H12 + H <=> C2H4 + pC4H9",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8e+21, 'cm^3/(mol*s)'),
        n = -2.39,
        Ea = (11180, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C6H12 + H <=> C2H4 + pC4H9""",
)

entry(
    index = 843,
    label = "C6H12 + H <=> C3H6 + nC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.6e+22, 'cm^3/(mol*s)'),
        n = -2.39,
        Ea = (11180, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C6H12 + H <=> C3H6 + nC3H7""",
)

entry(
    index = 844,
    label = "C6H12 + H <=> PXC6H11 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0323, 'cm^3/(mol*s)'), n=4.7, Ea=(3679, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12 + H <=> PXC6H11 + H2""",
)

entry(
    index = 845,
    label = "C6H12 + H <=> SXC6H11 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0317, 'cm^3/(mol*s)'), n=4.65, Ea=(1340, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12 + H <=> SXC6H11 + H2""",
)

entry(
    index = 846,
    label = "C6H12 + H <=> S2XC6H11 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0317, 'cm^3/(mol*s)'), n=4.65, Ea=(1340, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12 + H <=> S2XC6H11 + H2""",
)

entry(
    index = 847,
    label = "C6H12 + H <=> SAXC6H11 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(54000, 'cm^3/(mol*s)'), n=2.5, Ea=(-1900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12 + H <=> SAXC6H11 + H2""",
)

entry(
    index = 848,
    label = "C6H12 + O <=> PXC5H11 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.3e+08, 'cm^3/(mol*s)'),
        n = 1.45,
        Ea = (-402, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C6H12 + O <=> PXC5H11 + HCO""",
)

entry(
    index = 849,
    label = "C6H12 + O <=> PXC6H11 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.745, 'cm^3/(mol*s)'), n=4.17, Ea=(2766, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12 + O <=> PXC6H11 + OH""",
)

entry(
    index = 850,
    label = "C6H12 + OH <=> PXC6H11 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(700, 'cm^3/(mol*s)'), n=2.66, Ea=(527, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12 + OH <=> PXC6H11 + H2O""",
)

entry(
    index = 851,
    label = "C6H12 + O2 <=> PXC6H11 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(50930, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12 + O2 <=> PXC6H11 + HO2""",
)

entry(
    index = 852,
    label = "C6H12 + HO2 <=> PXC6H11 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(14340, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12 + HO2 <=> PXC6H11 + H2O2""",
)

entry(
    index = 853,
    label = "C6H12 + CH3 <=> PXC6H11 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.58e-08, 'cm^3/(mol*s)'),
        n = 6.08,
        Ea = (6223, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C6H12 + CH3 <=> PXC6H11 + CH4""",
)

entry(
    index = 854,
    label = "PXC6H13 <=> C2H4 + pC4H9",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.02e+12, 's^-1'), n=0.3, Ea=(27273.6, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (7.1e-35, 'cm^3/(mol*s)'),
            n = 15.411,
            Ea = (-600, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -5.91,
        T3 = (333, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC6H13 <=> C2H4 + pC4H9""",
)

entry(
    index = 855,
    label = "SXC6H13 <=> C3H6 + nC3H7",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(4.47e+11, 's^-1'), n=0.57, Ea=(28044.5, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.7e-33, 'cm^3/(mol*s)'),
            n = 14.91,
            Ea = (-600, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -6.53,
        T3 = (333, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SXC6H13 <=> C3H6 + nC3H7""",
)

entry(
    index = 856,
    label = "S2XC6H13 <=> C2H5 + C4H81",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.55e+12, 's^-1'), n=0.29, Ea=(28296.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4.5e-26, 'cm^3/(mol*s)'),
            n = 13.09,
            Ea = (-600.5, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -0.74,
        T3 = (308, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S2XC6H13 <=> C2H5 + C4H81""",
)

entry(
    index = 857,
    label = "S2XC6H13 <=> C5H10 + CH3",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(8.13e+10, 's^-1'), n=0.78, Ea=(29648, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4e-39, 'cm^3/(mol*s)'),
            n = 16.782,
            Ea = (-600.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -7.03,
        T3 = (314, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S2XC6H13 <=> C5H10 + CH3""",
)

entry(
    index = 858,
    label = "PXC6H13 + H <=> pC4H9 + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.7e+24, 'cm^3/(mol*s)'),
        n = -2.92,
        Ea = (12505, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is PXC6H13 + H <=> pC4H9 + C2H5""",
)

entry(
    index = 859,
    label = "PXC6H13 + H <=> C6H12 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC6H13 + H <=> C6H12 + H2""",
)

entry(
    index = 860,
    label = "PXC6H13 + O <=> PXC5H11 + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC6H13 + O <=> PXC5H11 + CH2O""",
)

entry(
    index = 861,
    label = "PXC6H13 + OH <=> C6H12 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC6H13 + OH <=> C6H12 + H2O""",
)

entry(
    index = 862,
    label = "PXC6H13 + O2 <=> C6H12 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+10, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC6H13 + O2 <=> C6H12 + HO2""",
)

entry(
    index = 863,
    label = "PXC6H13 + HO2 <=> PXC5H11 + OH + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC6H13 + HO2 <=> PXC5H11 + OH + CH2O""",
)

entry(
    index = 864,
    label = "PXC6H13 + HCO <=> NC6H14 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC6H13 + HCO <=> NC6H14 + CO""",
)

entry(
    index = 865,
    label = "PXC6H13 + CH3 <=> C6H12 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC6H13 + CH3 <=> C6H12 + CH4""",
)

entry(
    index = 866,
    label = "SXC6H13 + H <=> NC6H14",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.7e+58, 'cm^6/(mol^2*s)'),
            n = -12.08,
            Ea = (11263.7, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.649,
        T3 = (1213.1, 'K'),
        T1 = (1213.1, 'K'),
        T2 = (13369.7, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SXC6H13 + H <=> NC6H14""",
)

entry(
    index = 867,
    label = "SXC6H13 + H <=> pC4H9 + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.4e+28, 'cm^3/(mol*s)'),
        n = -3.94,
        Ea = (15916, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is SXC6H13 + H <=> pC4H9 + C2H5""",
)

entry(
    index = 868,
    label = "SXC6H13 + H <=> C6H12 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.2e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC6H13 + H <=> C6H12 + H2""",
)

entry(
    index = 869,
    label = "SXC6H13 + O <=> CH3CHO + pC4H9",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC6H13 + O <=> CH3CHO + pC4H9""",
)

entry(
    index = 870,
    label = "SXC6H13 + OH <=> C6H12 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC6H13 + OH <=> C6H12 + H2O""",
)

entry(
    index = 871,
    label = "SXC6H13 + O2 <=> C6H12 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC6H13 + O2 <=> C6H12 + HO2""",
)

entry(
    index = 872,
    label = "SXC6H13 + HO2 <=> CH3CHO + pC4H9 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC6H13 + HO2 <=> CH3CHO + pC4H9 + OH""",
)

entry(
    index = 873,
    label = "SXC6H13 + HCO <=> NC6H14 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC6H13 + HCO <=> NC6H14 + CO""",
)

entry(
    index = 874,
    label = "SXC6H13 + CH3 <=> CH4 + C6H12",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.2e+14, 'cm^3/(mol*s)'), n=-0.68, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC6H13 + CH3 <=> CH4 + C6H12""",
)

entry(
    index = 875,
    label = "S2XC6H13 + O2 <=> C6H12 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is S2XC6H13 + O2 <=> C6H12 + HO2""",
)

entry(
    index = 876,
    label = "PXC5H11 + CH3 <=> NC6H14",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.93e+14, 'cm^3/(mol*s)'), n=-0.32, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC5H11 + CH3 <=> NC6H14""",
)

entry(
    index = 877,
    label = "pC4H9 + C2H5 <=> NC6H14",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.88e+14, 'cm^3/(mol*s)'), n=-0.5, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is pC4H9 + C2H5 <=> NC6H14""",
)

entry(
    index = 878,
    label = "nC3H7 + nC3H7 <=> NC6H14",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.88e+14, 'cm^3/(mol*s)'), n=-0.5, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is nC3H7 + nC3H7 <=> NC6H14""",
)

entry(
    index = 879,
    label = "NC6H14 + OH <=> PXC6H13 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.73e+07, 'cm^3/(mol*s)'),
        n = 1.81,
        Ea = (868.3, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC6H14 + OH <=> PXC6H13 + H2O""",
)

entry(
    index = 880,
    label = "NC6H14 + OH <=> SXC6H13 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.41e+10, 'cm^3/(mol*s)'),
        n = 0.94,
        Ea = (504.7, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC6H14 + OH <=> SXC6H13 + H2O""",
)

entry(
    index = 881,
    label = "NC6H14 + OH <=> S2XC6H13 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.14e+07, 'cm^3/(mol*s)'),
        n = 1.81,
        Ea = (-1015.4, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC6H14 + OH <=> S2XC6H13 + H2O""",
)

entry(
    index = 882,
    label = "NC6H14 + O2 <=> PXC6H13 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(50930, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6H14 + O2 <=> PXC6H13 + HO2""",
)

entry(
    index = 883,
    label = "NC6H14 + O2 <=> SXC6H13 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+13, 'cm^3/(mol*s)'), n=0, Ea=(47590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6H14 + O2 <=> SXC6H13 + HO2""",
)

entry(
    index = 884,
    label = "NC6H14 + O2 <=> S2XC6H13 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+13, 'cm^3/(mol*s)'), n=0, Ea=(47590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6H14 + O2 <=> S2XC6H13 + HO2""",
)

entry(
    index = 885,
    label = "NC6H14 + HO2 <=> PXC6H13 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(61100, 'cm^3/(mol*s)'), n=2.65, Ea=(17496, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6H14 + HO2 <=> PXC6H13 + H2O2""",
)

entry(
    index = 886,
    label = "NC6H14 + HO2 <=> SXC6H13 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(14200, 'cm^3/(mol*s)'), n=2.77, Ea=(14913, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6H14 + HO2 <=> SXC6H13 + H2O2""",
)

entry(
    index = 887,
    label = "NC6H14 + HO2 <=> S2XC6H13 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(14200, 'cm^3/(mol*s)'), n=2.77, Ea=(14913, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6H14 + HO2 <=> S2XC6H13 + H2O2""",
)

entry(
    index = 888,
    label = "NC6H14 + H <=> PXC6H13 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0645, 'cm^3/(mol*s)'), n=4.7, Ea=(3679, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6H14 + H <=> PXC6H13 + H2""",
)

entry(
    index = 889,
    label = "NC6H14 + H <=> SXC6H13 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0634, 'cm^3/(mol*s)'), n=4.65, Ea=(1340, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6H14 + H <=> SXC6H13 + H2""",
)

entry(
    index = 890,
    label = "NC6H14 + H <=> S2XC6H13 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0634, 'cm^3/(mol*s)'), n=4.65, Ea=(1340, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6H14 + H <=> S2XC6H13 + H2""",
)

entry(
    index = 891,
    label = "NC6H14 + O <=> PXC6H13 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.49, 'cm^3/(mol*s)'), n=4.17, Ea=(2766, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6H14 + O <=> PXC6H13 + OH""",
)

entry(
    index = 892,
    label = "NC6H14 + O <=> SXC6H13 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(36.5, 'cm^3/(mol*s)'), n=3.75, Ea=(825, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6H14 + O <=> SXC6H13 + OH""",
)

entry(
    index = 893,
    label = "NC6H14 + O <=> S2XC6H13 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(36.5, 'cm^3/(mol*s)'), n=3.75, Ea=(825, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6H14 + O <=> S2XC6H13 + OH""",
)

entry(
    index = 894,
    label = "NC6H14 + CH3 <=> PXC6H13 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.16e-08, 'cm^3/(mol*s)'),
        n = 6.08,
        Ea = (6223, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC6H14 + CH3 <=> PXC6H13 + CH4""",
)

entry(
    index = 895,
    label = "NC6H14 + CH3 <=> SXC6H13 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.82e-07, 'cm^3/(mol*s)'),
        n = 5.89,
        Ea = (4768, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC6H14 + CH3 <=> SXC6H13 + CH4""",
)

entry(
    index = 896,
    label = "NC6H14 + CH3 <=> S2XC6H13 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.82e-07, 'cm^3/(mol*s)'),
        n = 5.89,
        Ea = (4768, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC6H14 + CH3 <=> S2XC6H13 + CH4""",
)

entry(
    index = 897,
    label = "PXC7H13 + H <=> C7H14",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.01e+48, 'cm^6/(mol^2*s)'),
            n = -9.32,
            Ea = (5833.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.498,
        T3 = (1314, 'K'),
        T1 = (1314, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC7H13 + H <=> C7H14""",
)

entry(
    index = 898,
    label = "PXC7H13 + H <=> CH3 + PXC6H11",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+21, 'cm^3/(mol*s)'), n=-2, Ea=(11000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC7H13 + H <=> CH3 + PXC6H11""",
)

entry(
    index = 899,
    label = "PXC7H13 + HO2 <=> CH2O + OH + PXC6H11",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC7H13 + HO2 <=> CH2O + OH + PXC6H11""",
)

entry(
    index = 900,
    label = "PXC7H13 + HCO <=> C7H14 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC7H13 + HCO <=> C7H14 + CO""",
)

entry(
    index = 901,
    label = "C2H4 + PXC5H9 <=> PXC7H13",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+11, 'cm^3/(mol*s)'), n=0, Ea=(7300, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4 + PXC5H9 <=> PXC7H13""",
)

entry(
    index = 902,
    label = "C7H14 <=> pC4H9 + aC3H5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.07e+23, 's^-1'), n=-2.03, Ea=(74958, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14 <=> pC4H9 + aC3H5""",
)

entry(
    index = 903,
    label = "C7H14 <=> C4H81 + C3H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.08e+06, 's^-1'), n=1.65, Ea=(53752, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14 <=> C4H81 + C3H6""",
)

entry(
    index = 904,
    label = "C7H14 + H <=> C2H4 + PXC5H11",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8e+21, 'cm^3/(mol*s)'),
        n = -2.39,
        Ea = (11180, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C7H14 + H <=> C2H4 + PXC5H11""",
)

entry(
    index = 905,
    label = "C7H14 + H <=> C3H6 + pC4H9",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.6e+22, 'cm^3/(mol*s)'),
        n = -2.39,
        Ea = (11180, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C7H14 + H <=> C3H6 + pC4H9""",
)

entry(
    index = 906,
    label = "C7H14 + H <=> PXC7H13 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0323, 'cm^3/(mol*s)'), n=4.7, Ea=(3679, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14 + H <=> PXC7H13 + H2""",
)

entry(
    index = 907,
    label = "C7H14 + O <=> PXC6H13 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.3e+08, 'cm^3/(mol*s)'),
        n = 1.45,
        Ea = (-402, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C7H14 + O <=> PXC6H13 + HCO""",
)

entry(
    index = 908,
    label = "C7H14 + O <=> PXC7H13 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.745, 'cm^3/(mol*s)'), n=4.17, Ea=(2766, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14 + O <=> PXC7H13 + OH""",
)

entry(
    index = 909,
    label = "C7H14 + OH <=> PXC7H13 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(700, 'cm^3/(mol*s)'), n=2.66, Ea=(527, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14 + OH <=> PXC7H13 + H2O""",
)

entry(
    index = 910,
    label = "C7H14 + O2 <=> PXC7H13 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(50930, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14 + O2 <=> PXC7H13 + HO2""",
)

entry(
    index = 911,
    label = "C7H14 + HO2 <=> PXC7H13 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(14340, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14 + HO2 <=> PXC7H13 + H2O2""",
)

entry(
    index = 912,
    label = "C7H14 + CH3 <=> PXC7H13 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.58e-08, 'cm^3/(mol*s)'),
        n = 6.08,
        Ea = (6223, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C7H14 + CH3 <=> PXC7H13 + CH4""",
)

entry(
    index = 913,
    label = "PXC7H15 <=> C2H4 + PXC5H11",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(7.94e+11, 's^-1'), n=0.33, Ea=(27210, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (2.8e-44, 'cm^3/(mol*s)'),
            n = 18.729,
            Ea = (-602.5, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -14.66,
        T3 = (219, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC7H15 <=> C2H4 + PXC5H11""",
)

entry(
    index = 914,
    label = "SXC7H15 <=> pC4H9 + C3H6",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(5.01e+11, 's^-1'), n=0.56, Ea=(28092.2, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (8.9e-39, 'cm^3/(mol*s)'),
            n = 16.934,
            Ea = (-602.5, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -25.27,
        T3 = (223, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SXC7H15 <=> pC4H9 + C3H6""",
)

entry(
    index = 915,
    label = "S2XC7H15 <=> nC3H7 + C4H81",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2.95e+12, 's^-1'), n=0.31, Ea=(28257.1, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (2e-38, 'cm^3/(mol*s)'),
            n = 16.814,
            Ea = (-602.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -20.96,
        T3 = (221, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S2XC7H15 <=> nC3H7 + C4H81""",
)

entry(
    index = 916,
    label = "S2XC7H15 <=> C6H12 + CH3",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.1e+11, 's^-1'), n=0.75, Ea=(29401.6, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.06e-42, 'cm^3/(mol*s)'),
            n = 18.004,
            Ea = (-602.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -20.94,
        T3 = (217, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S2XC7H15 <=> C6H12 + CH3""",
)

entry(
    index = 917,
    label = "S3XC7H15 <=> C2H5 + C5H10",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(5.89e+12, 's^-1'), n=0.31, Ea=(28257.1, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (2.1e-38, 'cm^3/(mol*s)'),
            n = 16.897,
            Ea = (-602.5, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -27.54,
        T3 = (224, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S3XC7H15 <=> C2H5 + C5H10""",
)

entry(
    index = 918,
    label = "PXC7H15 + H <=> NC7H16",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.01e+48, 'cm^6/(mol^2*s)'),
            n = -9.32,
            Ea = (5833.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.498,
        T3 = (1314, 'K'),
        T1 = (1314, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC7H15 + H <=> NC7H16""",
)

entry(
    index = 919,
    label = "PXC7H15 + H <=> PXC5H11 + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.7e+24, 'cm^3/(mol*s)'),
        n = -2.92,
        Ea = (12505, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is PXC7H15 + H <=> PXC5H11 + C2H5""",
)

entry(
    index = 920,
    label = "PXC7H15 + H <=> C7H14 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC7H15 + H <=> C7H14 + H2""",
)

entry(
    index = 921,
    label = "PXC7H15 + O <=> PXC6H13 + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC7H15 + O <=> PXC6H13 + CH2O""",
)

entry(
    index = 922,
    label = "PXC7H15 + OH <=> C7H14 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC7H15 + OH <=> C7H14 + H2O""",
)

entry(
    index = 923,
    label = "PXC7H15 + O2 <=> C7H14 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+10, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC7H15 + O2 <=> C7H14 + HO2""",
)

entry(
    index = 924,
    label = "PXC7H15 + HO2 <=> PXC6H13 + OH + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC7H15 + HO2 <=> PXC6H13 + OH + CH2O""",
)

entry(
    index = 925,
    label = "PXC7H15 + HCO <=> NC7H16 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC7H15 + HCO <=> NC7H16 + CO""",
)

entry(
    index = 926,
    label = "PXC7H15 + CH3 <=> C7H14 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC7H15 + CH3 <=> C7H14 + CH4""",
)

entry(
    index = 927,
    label = "SXC7H15 + H <=> PXC5H11 + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.4e+28, 'cm^3/(mol*s)'),
        n = -3.94,
        Ea = (15916, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is SXC7H15 + H <=> PXC5H11 + C2H5""",
)

entry(
    index = 928,
    label = "SXC7H15 + H <=> C7H14 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.2e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC7H15 + H <=> C7H14 + H2""",
)

entry(
    index = 929,
    label = "SXC7H15 + O <=> CH3CHO + PXC5H11",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC7H15 + O <=> CH3CHO + PXC5H11""",
)

entry(
    index = 930,
    label = "SXC7H15 + OH <=> C7H14 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC7H15 + OH <=> C7H14 + H2O""",
)

entry(
    index = 931,
    label = "SXC7H15 + O2 <=> C7H14 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC7H15 + O2 <=> C7H14 + HO2""",
)

entry(
    index = 932,
    label = "SXC7H15 + HO2 <=> CH3CHO + PXC5H11 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC7H15 + HO2 <=> CH3CHO + PXC5H11 + OH""",
)

entry(
    index = 933,
    label = "SXC7H15 + HCO <=> NC7H16 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC7H15 + HCO <=> NC7H16 + CO""",
)

entry(
    index = 934,
    label = "SXC7H15 + CH3 <=> CH4 + C7H14",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.2e+14, 'cm^3/(mol*s)'), n=-0.68, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC7H15 + CH3 <=> CH4 + C7H14""",
)

entry(
    index = 935,
    label = "S2XC7H15 + O2 <=> C7H14 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is S2XC7H15 + O2 <=> C7H14 + HO2""",
)

entry(
    index = 936,
    label = "S3XC7H15 + O2 <=> C7H14 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is S3XC7H15 + O2 <=> C7H14 + HO2""",
)

entry(
    index = 937,
    label = "PXC6H13 + CH3 <=> NC7H16",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.93e+14, 'cm^3/(mol*s)'), n=-0.32, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC6H13 + CH3 <=> NC7H16""",
)

entry(
    index = 938,
    label = "PXC5H11 + C2H5 <=> NC7H16",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.88e+14, 'cm^3/(mol*s)'), n=-0.5, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC5H11 + C2H5 <=> NC7H16""",
)

entry(
    index = 939,
    label = "pC4H9 + nC3H7 <=> NC7H16",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.88e+14, 'cm^3/(mol*s)'), n=-0.5, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is pC4H9 + nC3H7 <=> NC7H16""",
)

entry(
    index = 940,
    label = "NC7H16 + OH <=> PXC7H15 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.73e+07, 'cm^3/(mol*s)'),
        n = 1.81,
        Ea = (868.3, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC7H16 + OH <=> PXC7H15 + H2O""",
)

entry(
    index = 941,
    label = "NC7H16 + OH <=> SXC7H15 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.41e+10, 'cm^3/(mol*s)'),
        n = 0.94,
        Ea = (504.7, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC7H16 + OH <=> SXC7H15 + H2O""",
)

entry(
    index = 942,
    label = "NC7H16 + OH <=> S2XC7H15 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.14e+07, 'cm^3/(mol*s)'),
        n = 1.81,
        Ea = (-1015.4, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC7H16 + OH <=> S2XC7H15 + H2O""",
)

entry(
    index = 943,
    label = "NC7H16 + OH <=> S3XC7H15 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.62e+11, 'cm^3/(mol*s)'),
        n = 0.32,
        Ea = (846.5, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC7H16 + OH <=> S3XC7H15 + H2O""",
)

entry(
    index = 944,
    label = "NC7H16 + O2 <=> PXC7H15 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(50930, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + O2 <=> PXC7H15 + HO2""",
)

entry(
    index = 945,
    label = "NC7H16 + O2 <=> SXC7H15 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+13, 'cm^3/(mol*s)'), n=0, Ea=(47590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + O2 <=> SXC7H15 + HO2""",
)

entry(
    index = 946,
    label = "NC7H16 + O2 <=> S2XC7H15 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+13, 'cm^3/(mol*s)'), n=0, Ea=(47590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + O2 <=> S2XC7H15 + HO2""",
)

entry(
    index = 947,
    label = "NC7H16 + O2 <=> S3XC7H15 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(47590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + O2 <=> S3XC7H15 + HO2""",
)

entry(
    index = 948,
    label = "NC7H16 + HO2 <=> PXC7H15 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(61100, 'cm^3/(mol*s)'), n=2.65, Ea=(17496, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + HO2 <=> PXC7H15 + H2O2""",
)

entry(
    index = 949,
    label = "NC7H16 + HO2 <=> SXC7H15 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(14200, 'cm^3/(mol*s)'), n=2.77, Ea=(14913, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + HO2 <=> SXC7H15 + H2O2""",
)

entry(
    index = 950,
    label = "NC7H16 + HO2 <=> S2XC7H15 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(14200, 'cm^3/(mol*s)'), n=2.77, Ea=(14913, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + HO2 <=> S2XC7H15 + H2O2""",
)

entry(
    index = 951,
    label = "NC7H16 + HO2 <=> S3XC7H15 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7130, 'cm^3/(mol*s)'), n=2.77, Ea=(14913, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + HO2 <=> S3XC7H15 + H2O2""",
)

entry(
    index = 952,
    label = "NC7H16 + H <=> PXC7H15 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0645, 'cm^3/(mol*s)'), n=4.7, Ea=(3679, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + H <=> PXC7H15 + H2""",
)

entry(
    index = 953,
    label = "NC7H16 + H <=> SXC7H15 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0634, 'cm^3/(mol*s)'), n=4.65, Ea=(1340, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + H <=> SXC7H15 + H2""",
)

entry(
    index = 954,
    label = "NC7H16 + H <=> S2XC7H15 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0634, 'cm^3/(mol*s)'), n=4.65, Ea=(1340, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + H <=> S2XC7H15 + H2""",
)

entry(
    index = 955,
    label = "NC7H16 + H <=> S3XC7H15 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0317, 'cm^3/(mol*s)'), n=4.65, Ea=(1340, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + H <=> S3XC7H15 + H2""",
)

entry(
    index = 956,
    label = "NC7H16 + O <=> PXC7H15 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.49, 'cm^3/(mol*s)'), n=4.17, Ea=(2766, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + O <=> PXC7H15 + OH""",
)

entry(
    index = 957,
    label = "NC7H16 + O <=> SXC7H15 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(36.5, 'cm^3/(mol*s)'), n=3.75, Ea=(825, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + O <=> SXC7H15 + OH""",
)

entry(
    index = 958,
    label = "NC7H16 + O <=> S2XC7H15 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(36.5, 'cm^3/(mol*s)'), n=3.75, Ea=(825, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + O <=> S2XC7H15 + OH""",
)

entry(
    index = 959,
    label = "NC7H16 + O <=> S3XC7H15 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(18.2, 'cm^3/(mol*s)'), n=3.75, Ea=(825, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + O <=> S3XC7H15 + OH""",
)

entry(
    index = 960,
    label = "NC7H16 + CH3 <=> PXC7H15 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.16e-08, 'cm^3/(mol*s)'),
        n = 6.08,
        Ea = (6223, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC7H16 + CH3 <=> PXC7H15 + CH4""",
)

entry(
    index = 961,
    label = "NC7H16 + CH3 <=> SXC7H15 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.82e-07, 'cm^3/(mol*s)'),
        n = 5.89,
        Ea = (4768, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC7H16 + CH3 <=> SXC7H15 + CH4""",
)

entry(
    index = 962,
    label = "NC7H16 + CH3 <=> S2XC7H15 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.82e-07, 'cm^3/(mol*s)'),
        n = 5.89,
        Ea = (4768, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC7H16 + CH3 <=> S2XC7H15 + CH4""",
)

entry(
    index = 963,
    label = "NC7H16 + CH3 <=> S3XC7H15 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.41e-07, 'cm^3/(mol*s)'),
        n = 5.89,
        Ea = (4768, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC7H16 + CH3 <=> S3XC7H15 + CH4""",
)

entry(
    index = 964,
    label = "PXC8H15 + H <=> C8H16",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.01e+48, 'cm^6/(mol^2*s)'),
            n = -9.32,
            Ea = (5833.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.498,
        T3 = (1314, 'K'),
        T1 = (1314, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC8H15 + H <=> C8H16""",
)

entry(
    index = 965,
    label = "PXC8H15 + H <=> CH3 + PXC7H13",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+21, 'cm^3/(mol*s)'), n=-2, Ea=(11000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC8H15 + H <=> CH3 + PXC7H13""",
)

entry(
    index = 966,
    label = "PXC8H15 + HO2 <=> CH2O + OH + PXC7H13",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC8H15 + HO2 <=> CH2O + OH + PXC7H13""",
)

entry(
    index = 967,
    label = "PXC8H15 + HCO <=> C8H16 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC8H15 + HCO <=> C8H16 + CO""",
)

entry(
    index = 968,
    label = "C2H4 + PXC6H11 <=> PXC8H15",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+11, 'cm^3/(mol*s)'), n=0, Ea=(7300, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4 + PXC6H11 <=> PXC8H15""",
)

entry(
    index = 969,
    label = "C8H16 <=> PXC5H11 + aC3H5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.07e+23, 's^-1'), n=-2.03, Ea=(74958, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C8H16 <=> PXC5H11 + aC3H5""",
)

entry(
    index = 970,
    label = "C8H16 <=> C5H10 + C3H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.08e+06, 's^-1'), n=1.65, Ea=(53752, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C8H16 <=> C5H10 + C3H6""",
)

entry(
    index = 971,
    label = "C8H16 + H <=> C2H4 + PXC6H13",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8e+21, 'cm^3/(mol*s)'),
        n = -2.39,
        Ea = (11180, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C8H16 + H <=> C2H4 + PXC6H13""",
)

entry(
    index = 972,
    label = "C8H16 + H <=> C3H6 + PXC5H11",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.6e+22, 'cm^3/(mol*s)'),
        n = -2.39,
        Ea = (11180, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C8H16 + H <=> C3H6 + PXC5H11""",
)

entry(
    index = 973,
    label = "C8H16 + H <=> PXC8H15 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0323, 'cm^3/(mol*s)'), n=4.7, Ea=(3679, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C8H16 + H <=> PXC8H15 + H2""",
)

entry(
    index = 974,
    label = "C8H16 + O <=> PXC7H15 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.3e+08, 'cm^3/(mol*s)'),
        n = 1.45,
        Ea = (-402, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C8H16 + O <=> PXC7H15 + HCO""",
)

entry(
    index = 975,
    label = "C8H16 + O <=> PXC8H15 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.745, 'cm^3/(mol*s)'), n=4.17, Ea=(2766, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C8H16 + O <=> PXC8H15 + OH""",
)

entry(
    index = 976,
    label = "C8H16 + OH <=> PXC8H15 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(700, 'cm^3/(mol*s)'), n=2.66, Ea=(527, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C8H16 + OH <=> PXC8H15 + H2O""",
)

entry(
    index = 977,
    label = "C8H16 + O2 <=> PXC8H15 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(50930, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C8H16 + O2 <=> PXC8H15 + HO2""",
)

entry(
    index = 978,
    label = "C8H16 + HO2 <=> PXC8H15 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(14340, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C8H16 + HO2 <=> PXC8H15 + H2O2""",
)

entry(
    index = 979,
    label = "C8H16 + CH3 <=> PXC8H15 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.58e-08, 'cm^3/(mol*s)'),
        n = 6.08,
        Ea = (6223, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C8H16 + CH3 <=> PXC8H15 + CH4""",
)

entry(
    index = 980,
    label = "PXC8H17 <=> C2H4 + PXC6H13",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(9.12e+11, 's^-1'), n=0.31, Ea=(27237.8, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.8e-57, 'cm^3/(mol*s)'),
            n = 23.463,
            Ea = (-602.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -2.46,
        T3 = (206, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC8H17 <=> C2H4 + PXC6H13""",
)

entry(
    index = 981,
    label = "SXC8H17 <=> PXC5H11 + C3H6",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(6.03e+10, 's^-1'), n=0.84, Ea=(27820, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1e-43, 'cm^3/(mol*s)'),
            n = 18.591,
            Ea = (-602.5, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -43.32,
        T3 = (200, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SXC8H17 <=> PXC5H11 + C3H6""",
)

entry(
    index = 982,
    label = "S2XC8H17 <=> pC4H9 + C4H81",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2.04e+13, 's^-1'), n=0.04, Ea=(28493.6, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3e-43, 'cm^3/(mol*s)'),
            n = 18.43,
            Ea = (-602.8, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -34.47,
        T3 = (208, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S2XC8H17 <=> pC4H9 + C4H81""",
)

entry(
    index = 983,
    label = "S2XC8H17 <=> C7H14 + CH3",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(9.55e+09, 's^-1'), n=1.08, Ea=(29387.7, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.3e-46, 'cm^3/(mol*s)'),
            n = 19.133,
            Ea = (-602.7, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -34.36,
        T3 = (210, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S2XC8H17 <=> C7H14 + CH3""",
)

entry(
    index = 984,
    label = "S3XC8H17 <=> nC3H7 + C5H10",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(5.5e+11, 's^-1'), n=0.55, Ea=(28084.3, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.1e-43, 'cm^3/(mol*s)'),
            n = 18.418,
            Ea = (-602.9, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -32.13,
        T3 = (207, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S3XC8H17 <=> nC3H7 + C5H10""",
)

entry(
    index = 985,
    label = "S3XC8H17 <=> C6H12 + C2H5",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(6.76e+09, 's^-1'), n=1.11, Ea=(27023.2, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (8.2e-43, 'cm^3/(mol*s)'),
            n = 18.276,
            Ea = (-602.9, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -30.04,
        T3 = (210, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S3XC8H17 <=> C6H12 + C2H5""",
)

entry(
    index = 986,
    label = "PXC8H17 + H <=> NC8H18",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.01e+48, 'cm^6/(mol^2*s)'),
            n = -9.32,
            Ea = (5833.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.498,
        T3 = (1314, 'K'),
        T1 = (1314, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC8H17 + H <=> NC8H18""",
)

entry(
    index = 987,
    label = "PXC8H17 + H <=> PXC6H13 + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.7e+24, 'cm^3/(mol*s)'),
        n = -2.92,
        Ea = (12505, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is PXC8H17 + H <=> PXC6H13 + C2H5""",
)

entry(
    index = 988,
    label = "PXC8H17 + H <=> C8H16 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC8H17 + H <=> C8H16 + H2""",
)

entry(
    index = 989,
    label = "PXC8H17 + O <=> PXC7H15 + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC8H17 + O <=> PXC7H15 + CH2O""",
)

entry(
    index = 990,
    label = "PXC8H17 + OH <=> C8H16 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC8H17 + OH <=> C8H16 + H2O""",
)

entry(
    index = 991,
    label = "PXC8H17 + O2 <=> C8H16 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+10, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC8H17 + O2 <=> C8H16 + HO2""",
)

entry(
    index = 992,
    label = "PXC8H17 + HO2 <=> PXC7H15 + OH + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC8H17 + HO2 <=> PXC7H15 + OH + CH2O""",
)

entry(
    index = 993,
    label = "PXC8H17 + HCO <=> NC8H18 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC8H17 + HCO <=> NC8H18 + CO""",
)

entry(
    index = 994,
    label = "PXC8H17 + CH3 <=> C8H16 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC8H17 + CH3 <=> C8H16 + CH4""",
)

entry(
    index = 995,
    label = "SXC8H17 + H <=> NC8H18",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.7e+58, 'cm^6/(mol^2*s)'),
            n = -12.08,
            Ea = (11263.7, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.649,
        T3 = (1213.1, 'K'),
        T1 = (1213.1, 'K'),
        T2 = (13369.7, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SXC8H17 + H <=> NC8H18""",
)

entry(
    index = 996,
    label = "SXC8H17 + H <=> PXC6H13 + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.4e+28, 'cm^3/(mol*s)'),
        n = -3.94,
        Ea = (15916, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is SXC8H17 + H <=> PXC6H13 + C2H5""",
)

entry(
    index = 997,
    label = "SXC8H17 + H <=> C8H16 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.2e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC8H17 + H <=> C8H16 + H2""",
)

entry(
    index = 998,
    label = "SXC8H17 + O <=> CH3CHO + PXC6H13",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC8H17 + O <=> CH3CHO + PXC6H13""",
)

entry(
    index = 999,
    label = "SXC8H17 + OH <=> C8H16 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC8H17 + OH <=> C8H16 + H2O""",
)

entry(
    index = 1000,
    label = "SXC8H17 + O2 <=> C8H16 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC8H17 + O2 <=> C8H16 + HO2""",
)

entry(
    index = 1001,
    label = "SXC8H17 + HO2 <=> CH3CHO + PXC6H13 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC8H17 + HO2 <=> CH3CHO + PXC6H13 + OH""",
)

entry(
    index = 1002,
    label = "SXC8H17 + HCO <=> NC8H18 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC8H17 + HCO <=> NC8H18 + CO""",
)

entry(
    index = 1003,
    label = "SXC8H17 + CH3 <=> CH4 + C8H16",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.2e+14, 'cm^3/(mol*s)'), n=-0.68, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC8H17 + CH3 <=> CH4 + C8H16""",
)

entry(
    index = 1004,
    label = "S2XC8H17 + O2 <=> C8H16 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is S2XC8H17 + O2 <=> C8H16 + HO2""",
)

entry(
    index = 1005,
    label = "S3XC8H17 + O2 <=> C8H16 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is S3XC8H17 + O2 <=> C8H16 + HO2""",
)

entry(
    index = 1006,
    label = "PXC7H15 + CH3 <=> NC8H18",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.93e+14, 'cm^3/(mol*s)'), n=-0.32, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC7H15 + CH3 <=> NC8H18""",
)

entry(
    index = 1007,
    label = "PXC6H13 + C2H5 <=> NC8H18",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.88e+14, 'cm^3/(mol*s)'), n=-0.5, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC6H13 + C2H5 <=> NC8H18""",
)

entry(
    index = 1008,
    label = "PXC5H11 + nC3H7 <=> NC8H18",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.88e+14, 'cm^3/(mol*s)'), n=-0.5, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC5H11 + nC3H7 <=> NC8H18""",
)

entry(
    index = 1009,
    label = "pC4H9 + pC4H9 <=> NC8H18",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.88e+14, 'cm^3/(mol*s)'), n=-0.5, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is pC4H9 + pC4H9 <=> NC8H18""",
)

entry(
    index = 1010,
    label = "NC8H18 + OH <=> PXC8H17 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.73e+07, 'cm^3/(mol*s)'),
        n = 1.81,
        Ea = (868.3, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC8H18 + OH <=> PXC8H17 + H2O""",
)

entry(
    index = 1011,
    label = "NC8H18 + OH <=> SXC8H17 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.41e+10, 'cm^3/(mol*s)'),
        n = 0.94,
        Ea = (504.7, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC8H18 + OH <=> SXC8H17 + H2O""",
)

entry(
    index = 1012,
    label = "NC8H18 + OH <=> S2XC8H17 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.14e+07, 'cm^3/(mol*s)'),
        n = 1.81,
        Ea = (-1015.4, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC8H18 + OH <=> S2XC8H17 + H2O""",
)

entry(
    index = 1013,
    label = "NC8H18 + OH <=> S3XC8H17 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.12e+12, 'cm^3/(mol*s)'),
        n = 0.32,
        Ea = (846.5, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC8H18 + OH <=> S3XC8H17 + H2O""",
)

entry(
    index = 1014,
    label = "NC8H18 + O2 <=> PXC8H17 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(50930, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC8H18 + O2 <=> PXC8H17 + HO2""",
)

entry(
    index = 1015,
    label = "NC8H18 + O2 <=> SXC8H17 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+13, 'cm^3/(mol*s)'), n=0, Ea=(47590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC8H18 + O2 <=> SXC8H17 + HO2""",
)

entry(
    index = 1016,
    label = "NC8H18 + O2 <=> S2XC8H17 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+13, 'cm^3/(mol*s)'), n=0, Ea=(47590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC8H18 + O2 <=> S2XC8H17 + HO2""",
)

entry(
    index = 1017,
    label = "NC8H18 + O2 <=> S3XC8H17 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+13, 'cm^3/(mol*s)'), n=0, Ea=(47590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC8H18 + O2 <=> S3XC8H17 + HO2""",
)

entry(
    index = 1018,
    label = "NC8H18 + HO2 <=> PXC8H17 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(61100, 'cm^3/(mol*s)'), n=2.65, Ea=(17496, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC8H18 + HO2 <=> PXC8H17 + H2O2""",
)

entry(
    index = 1019,
    label = "NC8H18 + HO2 <=> SXC8H17 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(14200, 'cm^3/(mol*s)'), n=2.77, Ea=(14913, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC8H18 + HO2 <=> SXC8H17 + H2O2""",
)

entry(
    index = 1020,
    label = "NC8H18 + HO2 <=> S2XC8H17 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(14200, 'cm^3/(mol*s)'), n=2.77, Ea=(14913, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC8H18 + HO2 <=> S2XC8H17 + H2O2""",
)

entry(
    index = 1021,
    label = "NC8H18 + HO2 <=> S3XC8H17 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(14200, 'cm^3/(mol*s)'), n=2.77, Ea=(14913, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC8H18 + HO2 <=> S3XC8H17 + H2O2""",
)

entry(
    index = 1022,
    label = "NC8H18 + H <=> PXC8H17 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0645, 'cm^3/(mol*s)'), n=4.7, Ea=(3679, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC8H18 + H <=> PXC8H17 + H2""",
)

entry(
    index = 1023,
    label = "NC8H18 + H <=> SXC8H17 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0634, 'cm^3/(mol*s)'), n=4.65, Ea=(1340, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC8H18 + H <=> SXC8H17 + H2""",
)

entry(
    index = 1024,
    label = "NC8H18 + H <=> S2XC8H17 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0634, 'cm^3/(mol*s)'), n=4.65, Ea=(1340, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC8H18 + H <=> S2XC8H17 + H2""",
)

entry(
    index = 1025,
    label = "NC8H18 + H <=> S3XC8H17 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0634, 'cm^3/(mol*s)'), n=4.65, Ea=(1340, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC8H18 + H <=> S3XC8H17 + H2""",
)

entry(
    index = 1026,
    label = "NC8H18 + O <=> PXC8H17 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.49, 'cm^3/(mol*s)'), n=4.17, Ea=(2766, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC8H18 + O <=> PXC8H17 + OH""",
)

entry(
    index = 1027,
    label = "NC8H18 + O <=> SXC8H17 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(36.5, 'cm^3/(mol*s)'), n=3.75, Ea=(825, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC8H18 + O <=> SXC8H17 + OH""",
)

entry(
    index = 1028,
    label = "NC8H18 + O <=> S2XC8H17 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(36.5, 'cm^3/(mol*s)'), n=3.75, Ea=(825, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC8H18 + O <=> S2XC8H17 + OH""",
)

entry(
    index = 1029,
    label = "NC8H18 + O <=> S3XC8H17 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(36.5, 'cm^3/(mol*s)'), n=3.75, Ea=(825, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC8H18 + O <=> S3XC8H17 + OH""",
)

entry(
    index = 1030,
    label = "NC8H18 + CH3 <=> PXC8H17 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.16e-08, 'cm^3/(mol*s)'),
        n = 6.08,
        Ea = (6223, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC8H18 + CH3 <=> PXC8H17 + CH4""",
)

entry(
    index = 1031,
    label = "NC8H18 + CH3 <=> SXC8H17 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.82e-07, 'cm^3/(mol*s)'),
        n = 5.89,
        Ea = (4768, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC8H18 + CH3 <=> SXC8H17 + CH4""",
)

entry(
    index = 1032,
    label = "NC8H18 + CH3 <=> S2XC8H17 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.82e-07, 'cm^3/(mol*s)'),
        n = 5.89,
        Ea = (4768, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC8H18 + CH3 <=> S2XC8H17 + CH4""",
)

entry(
    index = 1033,
    label = "NC8H18 + CH3 <=> S3XC8H17 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.82e-07, 'cm^3/(mol*s)'),
        n = 5.89,
        Ea = (4768, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC8H18 + CH3 <=> S3XC8H17 + CH4""",
)

entry(
    index = 1034,
    label = "PXC9H17 + H <=> C9H18",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.01e+48, 'cm^6/(mol^2*s)'),
            n = -9.32,
            Ea = (5833.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.498,
        T3 = (1314, 'K'),
        T1 = (1314, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC9H17 + H <=> C9H18""",
)

entry(
    index = 1035,
    label = "PXC9H17 + H <=> CH3 + PXC8H15",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+21, 'cm^3/(mol*s)'), n=-2, Ea=(11000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC9H17 + H <=> CH3 + PXC8H15""",
)

entry(
    index = 1036,
    label = "PXC9H17 + HO2 <=> CH2O + OH + PXC8H15",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC9H17 + HO2 <=> CH2O + OH + PXC8H15""",
)

entry(
    index = 1037,
    label = "PXC9H17 + HCO <=> C9H18 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC9H17 + HCO <=> C9H18 + CO""",
)

entry(
    index = 1038,
    label = "C2H4 + PXC7H13 <=> PXC9H17",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+11, 'cm^3/(mol*s)'), n=0, Ea=(7300, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4 + PXC7H13 <=> PXC9H17""",
)

entry(
    index = 1039,
    label = "C9H18 <=> PXC6H13 + aC3H5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.07e+23, 's^-1'), n=-2.03, Ea=(74958, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C9H18 <=> PXC6H13 + aC3H5""",
)

entry(
    index = 1040,
    label = "C9H18 <=> C6H12 + C3H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.08e+06, 's^-1'), n=1.65, Ea=(53752, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C9H18 <=> C6H12 + C3H6""",
)

entry(
    index = 1041,
    label = "C9H18 + H <=> C2H4 + PXC7H15",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8e+21, 'cm^3/(mol*s)'),
        n = -2.39,
        Ea = (11180, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C9H18 + H <=> C2H4 + PXC7H15""",
)

entry(
    index = 1042,
    label = "C9H18 + H <=> C3H6 + PXC6H13",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.6e+22, 'cm^3/(mol*s)'),
        n = -2.39,
        Ea = (11180, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C9H18 + H <=> C3H6 + PXC6H13""",
)

entry(
    index = 1043,
    label = "C9H18 + H <=> PXC9H17 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0323, 'cm^3/(mol*s)'), n=4.7, Ea=(3679, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C9H18 + H <=> PXC9H17 + H2""",
)

entry(
    index = 1044,
    label = "C9H18 + O <=> PXC8H17 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.3e+08, 'cm^3/(mol*s)'),
        n = 1.45,
        Ea = (-402, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C9H18 + O <=> PXC8H17 + HCO""",
)

entry(
    index = 1045,
    label = "C9H18 + O <=> PXC9H17 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.745, 'cm^3/(mol*s)'), n=4.17, Ea=(2766, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C9H18 + O <=> PXC9H17 + OH""",
)

entry(
    index = 1046,
    label = "C9H18 + OH <=> PXC9H17 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(700, 'cm^3/(mol*s)'), n=2.66, Ea=(527, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C9H18 + OH <=> PXC9H17 + H2O""",
)

entry(
    index = 1047,
    label = "C9H18 + O2 <=> PXC9H17 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(50930, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C9H18 + O2 <=> PXC9H17 + HO2""",
)

entry(
    index = 1048,
    label = "C9H18 + HO2 <=> PXC9H17 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(14340, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C9H18 + HO2 <=> PXC9H17 + H2O2""",
)

entry(
    index = 1049,
    label = "C9H18 + CH3 <=> PXC9H17 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.58e-08, 'cm^3/(mol*s)'),
        n = 6.08,
        Ea = (6223, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C9H18 + CH3 <=> PXC9H17 + CH4""",
)

entry(
    index = 1050,
    label = "PXC9H19 <=> C2H4 + PXC7H15",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(9.12e+11, 's^-1'), n=0.31, Ea=(27237.8, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.8e-57, 'cm^3/(mol*s)'),
            n = 23.463,
            Ea = (-602.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -2.46,
        T3 = (206, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC9H19 <=> C2H4 + PXC7H15""",
)

entry(
    index = 1051,
    label = "SXC9H19 <=> C3H6 + PXC6H13",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(6.03e+10, 's^-1'), n=0.84, Ea=(27820, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1e-43, 'cm^3/(mol*s)'),
            n = 18.591,
            Ea = (-602.5, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -43.32,
        T3 = (200, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SXC9H19 <=> C3H6 + PXC6H13""",
)

entry(
    index = 1052,
    label = "S2XC9H19 <=> PXC5H11 + C4H81",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2.04e+13, 's^-1'), n=0.04, Ea=(28493.6, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3e-43, 'cm^3/(mol*s)'),
            n = 18.43,
            Ea = (-602.8, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -34.47,
        T3 = (208, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S2XC9H19 <=> PXC5H11 + C4H81""",
)

entry(
    index = 1053,
    label = "S2XC9H19 <=> C8H16 + CH3",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(9.55e+09, 's^-1'), n=1.08, Ea=(29387.7, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.3e-46, 'cm^3/(mol*s)'),
            n = 19.133,
            Ea = (-602.7, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -34.36,
        T3 = (210, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S2XC9H19 <=> C8H16 + CH3""",
)

entry(
    index = 1054,
    label = "S3XC9H19 <=> pC4H9 + C5H10",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(5.5e+11, 's^-1'), n=0.55, Ea=(28084.3, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.1e-43, 'cm^3/(mol*s)'),
            n = 18.418,
            Ea = (-602.9, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -32.13,
        T3 = (207, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S3XC9H19 <=> pC4H9 + C5H10""",
)

entry(
    index = 1055,
    label = "S3XC9H19 <=> C7H14 + C2H5",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(6.76e+09, 's^-1'), n=1.11, Ea=(27023.2, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (8.2e-43, 'cm^3/(mol*s)'),
            n = 18.276,
            Ea = (-602.9, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -30.04,
        T3 = (210, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S3XC9H19 <=> C7H14 + C2H5""",
)

entry(
    index = 1056,
    label = "S4XC9H19 <=> nC3H7 + C6H12",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.1e+12, 's^-1'), n=0.55, Ea=(28084.3, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (6.2e-43, 'cm^3/(mol*s)'),
            n = 18.418,
            Ea = (-602.9, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -32.13,
        T3 = (207, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S4XC9H19 <=> nC3H7 + C6H12""",
)

entry(
    index = 1057,
    label = "PXC9H19 + H <=> NC9H20",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.01e+48, 'cm^6/(mol^2*s)'),
            n = -9.32,
            Ea = (5833.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.498,
        T3 = (1314, 'K'),
        T1 = (1314, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC9H19 + H <=> NC9H20""",
)

entry(
    index = 1058,
    label = "PXC9H19 + H <=> PXC7H15 + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.7e+24, 'cm^3/(mol*s)'),
        n = -2.92,
        Ea = (12505, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is PXC9H19 + H <=> PXC7H15 + C2H5""",
)

entry(
    index = 1059,
    label = "PXC9H19 + H <=> C9H18 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC9H19 + H <=> C9H18 + H2""",
)

entry(
    index = 1060,
    label = "PXC9H19 + O <=> PXC8H17 + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC9H19 + O <=> PXC8H17 + CH2O""",
)

entry(
    index = 1061,
    label = "PXC9H19 + OH <=> C9H18 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC9H19 + OH <=> C9H18 + H2O""",
)

entry(
    index = 1062,
    label = "PXC9H19 + O2 <=> C9H18 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+10, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC9H19 + O2 <=> C9H18 + HO2""",
)

entry(
    index = 1063,
    label = "PXC9H19 + HO2 <=> PXC8H17 + OH + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC9H19 + HO2 <=> PXC8H17 + OH + CH2O""",
)

entry(
    index = 1064,
    label = "PXC9H19 + HCO <=> NC9H20 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC9H19 + HCO <=> NC9H20 + CO""",
)

entry(
    index = 1065,
    label = "PXC9H19 + CH3 <=> C9H18 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC9H19 + CH3 <=> C9H18 + CH4""",
)

entry(
    index = 1066,
    label = "SXC9H19 + H <=> NC9H20",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.7e+58, 'cm^6/(mol^2*s)'),
            n = -12.08,
            Ea = (11263.7, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.649,
        T3 = (1213.1, 'K'),
        T1 = (1213.1, 'K'),
        T2 = (13369.7, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SXC9H19 + H <=> NC9H20""",
)

entry(
    index = 1067,
    label = "SXC9H19 + H <=> PXC7H15 + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.4e+28, 'cm^3/(mol*s)'),
        n = -3.94,
        Ea = (15916, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is SXC9H19 + H <=> PXC7H15 + C2H5""",
)

entry(
    index = 1068,
    label = "SXC9H19 + H <=> C9H18 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.2e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC9H19 + H <=> C9H18 + H2""",
)

entry(
    index = 1069,
    label = "SXC9H19 + O <=> CH3CHO + PXC7H15",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC9H19 + O <=> CH3CHO + PXC7H15""",
)

entry(
    index = 1070,
    label = "SXC9H19 + OH <=> C9H18 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC9H19 + OH <=> C9H18 + H2O""",
)

entry(
    index = 1071,
    label = "SXC9H19 + O2 <=> C9H18 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC9H19 + O2 <=> C9H18 + HO2""",
)

entry(
    index = 1072,
    label = "SXC9H19 + HO2 <=> CH3CHO + PXC7H15 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC9H19 + HO2 <=> CH3CHO + PXC7H15 + OH""",
)

entry(
    index = 1073,
    label = "SXC9H19 + HCO <=> NC9H20 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC9H19 + HCO <=> NC9H20 + CO""",
)

entry(
    index = 1074,
    label = "SXC9H19 + CH3 <=> CH4 + C9H18",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.2e+14, 'cm^3/(mol*s)'), n=-0.68, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC9H19 + CH3 <=> CH4 + C9H18""",
)

entry(
    index = 1075,
    label = "S2XC9H19 + O2 <=> C9H18 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is S2XC9H19 + O2 <=> C9H18 + HO2""",
)

entry(
    index = 1076,
    label = "S3XC9H19 + O2 <=> C9H18 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is S3XC9H19 + O2 <=> C9H18 + HO2""",
)

entry(
    index = 1077,
    label = "S4XC9H19 + O2 <=> C9H18 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is S4XC9H19 + O2 <=> C9H18 + HO2""",
)

entry(
    index = 1078,
    label = "PXC8H17 + CH3 <=> NC9H20",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.93e+14, 'cm^3/(mol*s)'), n=-0.32, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC8H17 + CH3 <=> NC9H20""",
)

entry(
    index = 1079,
    label = "PXC7H15 + C2H5 <=> NC9H20",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.88e+14, 'cm^3/(mol*s)'), n=-0.5, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC7H15 + C2H5 <=> NC9H20""",
)

entry(
    index = 1080,
    label = "PXC6H13 + nC3H7 <=> NC9H20",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.88e+14, 'cm^3/(mol*s)'), n=-0.5, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC6H13 + nC3H7 <=> NC9H20""",
)

entry(
    index = 1081,
    label = "PXC5H11 + pC4H9 <=> NC9H20",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.88e+14, 'cm^3/(mol*s)'), n=-0.5, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC5H11 + pC4H9 <=> NC9H20""",
)

entry(
    index = 1082,
    label = "NC9H20 + OH <=> PXC9H19 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.73e+07, 'cm^3/(mol*s)'),
        n = 1.81,
        Ea = (868.3, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC9H20 + OH <=> PXC9H19 + H2O""",
)

entry(
    index = 1083,
    label = "NC9H20 + OH <=> SXC9H19 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.41e+10, 'cm^3/(mol*s)'),
        n = 0.94,
        Ea = (504.7, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC9H20 + OH <=> SXC9H19 + H2O""",
)

entry(
    index = 1084,
    label = "NC9H20 + OH <=> S2XC9H19 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.14e+07, 'cm^3/(mol*s)'),
        n = 1.81,
        Ea = (-1015.4, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC9H20 + OH <=> S2XC9H19 + H2O""",
)

entry(
    index = 1085,
    label = "NC9H20 + OH <=> S3XC9H19 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.12e+12, 'cm^3/(mol*s)'),
        n = 0.32,
        Ea = (846.5, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC9H20 + OH <=> S3XC9H19 + H2O""",
)

entry(
    index = 1086,
    label = "NC9H20 + OH <=> S4XC9H19 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.62e+11, 'cm^3/(mol*s)'),
        n = 0.32,
        Ea = (846.5, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC9H20 + OH <=> S4XC9H19 + H2O""",
)

entry(
    index = 1087,
    label = "NC9H20 + O2 <=> PXC9H19 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(50930, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC9H20 + O2 <=> PXC9H19 + HO2""",
)

entry(
    index = 1088,
    label = "NC9H20 + O2 <=> SXC9H19 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+13, 'cm^3/(mol*s)'), n=0, Ea=(47590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC9H20 + O2 <=> SXC9H19 + HO2""",
)

entry(
    index = 1089,
    label = "NC9H20 + O2 <=> S2XC9H19 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+13, 'cm^3/(mol*s)'), n=0, Ea=(47590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC9H20 + O2 <=> S2XC9H19 + HO2""",
)

entry(
    index = 1090,
    label = "NC9H20 + O2 <=> S3XC9H19 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+13, 'cm^3/(mol*s)'), n=0, Ea=(47590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC9H20 + O2 <=> S3XC9H19 + HO2""",
)

entry(
    index = 1091,
    label = "NC9H20 + O2 <=> S4XC9H19 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(47590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC9H20 + O2 <=> S4XC9H19 + HO2""",
)

entry(
    index = 1092,
    label = "NC9H20 + HO2 <=> PXC9H19 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(61100, 'cm^3/(mol*s)'), n=2.65, Ea=(17496, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC9H20 + HO2 <=> PXC9H19 + H2O2""",
)

entry(
    index = 1093,
    label = "NC9H20 + HO2 <=> SXC9H19 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(14200, 'cm^3/(mol*s)'), n=2.77, Ea=(14913, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC9H20 + HO2 <=> SXC9H19 + H2O2""",
)

entry(
    index = 1094,
    label = "NC9H20 + HO2 <=> S2XC9H19 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(14200, 'cm^3/(mol*s)'), n=2.77, Ea=(14913, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC9H20 + HO2 <=> S2XC9H19 + H2O2""",
)

entry(
    index = 1095,
    label = "NC9H20 + HO2 <=> S3XC9H19 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(14200, 'cm^3/(mol*s)'), n=2.77, Ea=(14913, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC9H20 + HO2 <=> S3XC9H19 + H2O2""",
)

entry(
    index = 1096,
    label = "NC9H20 + HO2 <=> S4XC9H19 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7130, 'cm^3/(mol*s)'), n=2.77, Ea=(14913, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC9H20 + HO2 <=> S4XC9H19 + H2O2""",
)

entry(
    index = 1097,
    label = "NC9H20 + H <=> PXC9H19 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0645, 'cm^3/(mol*s)'), n=4.7, Ea=(3679, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC9H20 + H <=> PXC9H19 + H2""",
)

entry(
    index = 1098,
    label = "NC9H20 + H <=> SXC9H19 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0634, 'cm^3/(mol*s)'), n=4.65, Ea=(1340, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC9H20 + H <=> SXC9H19 + H2""",
)

entry(
    index = 1099,
    label = "NC9H20 + H <=> S2XC9H19 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0634, 'cm^3/(mol*s)'), n=4.65, Ea=(1340, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC9H20 + H <=> S2XC9H19 + H2""",
)

entry(
    index = 1100,
    label = "NC9H20 + H <=> S3XC9H19 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0634, 'cm^3/(mol*s)'), n=4.65, Ea=(1340, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC9H20 + H <=> S3XC9H19 + H2""",
)

entry(
    index = 1101,
    label = "NC9H20 + H <=> S4XC9H19 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0317, 'cm^3/(mol*s)'), n=4.65, Ea=(1340, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC9H20 + H <=> S4XC9H19 + H2""",
)

entry(
    index = 1102,
    label = "NC9H20 + O <=> PXC9H19 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.49, 'cm^3/(mol*s)'), n=4.17, Ea=(2766, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC9H20 + O <=> PXC9H19 + OH""",
)

entry(
    index = 1103,
    label = "NC9H20 + O <=> SXC9H19 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(36.5, 'cm^3/(mol*s)'), n=3.75, Ea=(825, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC9H20 + O <=> SXC9H19 + OH""",
)

entry(
    index = 1104,
    label = "NC9H20 + O <=> S2XC9H19 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(36.5, 'cm^3/(mol*s)'), n=3.75, Ea=(825, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC9H20 + O <=> S2XC9H19 + OH""",
)

entry(
    index = 1105,
    label = "NC9H20 + O <=> S3XC9H19 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(36.5, 'cm^3/(mol*s)'), n=3.75, Ea=(825, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC9H20 + O <=> S3XC9H19 + OH""",
)

entry(
    index = 1106,
    label = "NC9H20 + O <=> S4XC9H19 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(18.2, 'cm^3/(mol*s)'), n=3.75, Ea=(825, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC9H20 + O <=> S4XC9H19 + OH""",
)

entry(
    index = 1107,
    label = "NC9H20 + CH3 <=> PXC9H19 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.16e-08, 'cm^3/(mol*s)'),
        n = 6.08,
        Ea = (6223, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC9H20 + CH3 <=> PXC9H19 + CH4""",
)

entry(
    index = 1108,
    label = "NC9H20 + CH3 <=> SXC9H19 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.82e-07, 'cm^3/(mol*s)'),
        n = 5.89,
        Ea = (4768, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC9H20 + CH3 <=> SXC9H19 + CH4""",
)

entry(
    index = 1109,
    label = "NC9H20 + CH3 <=> S2XC9H19 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.82e-07, 'cm^3/(mol*s)'),
        n = 5.89,
        Ea = (4768, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC9H20 + CH3 <=> S2XC9H19 + CH4""",
)

entry(
    index = 1110,
    label = "NC9H20 + CH3 <=> S3XC9H19 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.82e-07, 'cm^3/(mol*s)'),
        n = 5.89,
        Ea = (4768, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC9H20 + CH3 <=> S3XC9H19 + CH4""",
)

entry(
    index = 1111,
    label = "NC9H20 + CH3 <=> S4XC9H19 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.41e-07, 'cm^3/(mol*s)'),
        n = 5.89,
        Ea = (4768, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC9H20 + CH3 <=> S4XC9H19 + CH4""",
)

entry(
    index = 1112,
    label = "PXC10H19 + H <=> C10H20",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.01e+48, 'cm^6/(mol^2*s)'),
            n = -9.32,
            Ea = (5833.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.498,
        T3 = (1314, 'K'),
        T1 = (1314, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC10H19 + H <=> C10H20""",
)

entry(
    index = 1113,
    label = "PXC10H19 + H <=> CH3 + PXC9H17",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+21, 'cm^3/(mol*s)'), n=-2, Ea=(11000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC10H19 + H <=> CH3 + PXC9H17""",
)

entry(
    index = 1114,
    label = "PXC10H19 + HO2 <=> CH2O + OH + PXC9H17",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC10H19 + HO2 <=> CH2O + OH + PXC9H17""",
)

entry(
    index = 1115,
    label = "PXC10H19 + HCO <=> C10H20 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC10H19 + HCO <=> C10H20 + CO""",
)

entry(
    index = 1116,
    label = "C2H4 + PXC8H15 <=> PXC10H19",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+11, 'cm^3/(mol*s)'), n=0, Ea=(7300, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4 + PXC8H15 <=> PXC10H19""",
)

entry(
    index = 1117,
    label = "C10H20 <=> PXC7H15 + aC3H5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.07e+23, 's^-1'), n=-2.03, Ea=(74958, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C10H20 <=> PXC7H15 + aC3H5""",
)

entry(
    index = 1118,
    label = "C10H20 <=> C7H14 + C3H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.08e+06, 's^-1'), n=1.65, Ea=(53752, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C10H20 <=> C7H14 + C3H6""",
)

entry(
    index = 1119,
    label = "C10H20 + H <=> C2H4 + PXC8H17",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8e+21, 'cm^3/(mol*s)'),
        n = -2.39,
        Ea = (11180, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C10H20 + H <=> C2H4 + PXC8H17""",
)

entry(
    index = 1120,
    label = "C10H20 + H <=> C3H6 + PXC7H15",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.6e+22, 'cm^3/(mol*s)'),
        n = -2.39,
        Ea = (11180, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C10H20 + H <=> C3H6 + PXC7H15""",
)

entry(
    index = 1121,
    label = "C10H20 + H <=> PXC10H19 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0323, 'cm^3/(mol*s)'), n=4.7, Ea=(3679, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C10H20 + H <=> PXC10H19 + H2""",
)

entry(
    index = 1122,
    label = "C10H20 + O <=> PXC9H19 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.3e+08, 'cm^3/(mol*s)'),
        n = 1.45,
        Ea = (-402, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C10H20 + O <=> PXC9H19 + HCO""",
)

entry(
    index = 1123,
    label = "C10H20 + O <=> PXC10H19 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.745, 'cm^3/(mol*s)'), n=4.17, Ea=(2766, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C10H20 + O <=> PXC10H19 + OH""",
)

entry(
    index = 1124,
    label = "C10H20 + OH <=> PXC10H19 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(700, 'cm^3/(mol*s)'), n=2.66, Ea=(527, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C10H20 + OH <=> PXC10H19 + H2O""",
)

entry(
    index = 1125,
    label = "C10H20 + O2 <=> PXC10H19 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(50930, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C10H20 + O2 <=> PXC10H19 + HO2""",
)

entry(
    index = 1126,
    label = "C10H20 + HO2 <=> PXC10H19 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(14340, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C10H20 + HO2 <=> PXC10H19 + H2O2""",
)

entry(
    index = 1127,
    label = "C10H20 + CH3 <=> PXC10H19 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.58e-08, 'cm^3/(mol*s)'),
        n = 6.08,
        Ea = (6223, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C10H20 + CH3 <=> PXC10H19 + CH4""",
)

entry(
    index = 1128,
    label = "PXC10H21 <=> C2H4 + PXC8H17",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(9.12e+11, 's^-1'), n=0.31, Ea=(27237.8, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.8e-57, 'cm^3/(mol*s)'),
            n = 23.463,
            Ea = (-602.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -2.46,
        T3 = (206, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC10H21 <=> C2H4 + PXC8H17""",
)

entry(
    index = 1129,
    label = "SXC10H21 <=> C3H6 + PXC7H15",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(6.03e+10, 's^-1'), n=0.84, Ea=(27820, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1e-43, 'cm^3/(mol*s)'),
            n = 18.591,
            Ea = (-602.5, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -43.32,
        T3 = (200, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SXC10H21 <=> C3H6 + PXC7H15""",
)

entry(
    index = 1130,
    label = "S2XC10H21 <=> PXC6H13 + C4H81",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2.04e+13, 's^-1'), n=0.04, Ea=(28493.6, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3e-43, 'cm^3/(mol*s)'),
            n = 18.43,
            Ea = (-602.8, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -34.47,
        T3 = (208, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S2XC10H21 <=> PXC6H13 + C4H81""",
)

entry(
    index = 1131,
    label = "S2XC10H21 <=> C9H18 + CH3",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(9.55e+09, 's^-1'), n=1.08, Ea=(29387.7, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.3e-46, 'cm^3/(mol*s)'),
            n = 19.133,
            Ea = (-602.7, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -34.36,
        T3 = (210, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S2XC10H21 <=> C9H18 + CH3""",
)

entry(
    index = 1132,
    label = "S3XC10H21 <=> PXC5H11 + C5H10",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(5.5e+11, 's^-1'), n=0.55, Ea=(28084.3, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.1e-43, 'cm^3/(mol*s)'),
            n = 18.418,
            Ea = (-602.9, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -32.13,
        T3 = (207, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S3XC10H21 <=> PXC5H11 + C5H10""",
)

entry(
    index = 1133,
    label = "S3XC10H21 <=> C8H16 + C2H5",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(6.76e+09, 's^-1'), n=1.11, Ea=(27023.2, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (8.2e-43, 'cm^3/(mol*s)'),
            n = 18.276,
            Ea = (-602.9, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -30.04,
        T3 = (210, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S3XC10H21 <=> C8H16 + C2H5""",
)

entry(
    index = 1134,
    label = "S4XC10H21 <=> pC4H9 + C6H12",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(5.5e+11, 's^-1'), n=0.55, Ea=(28084.3, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.1e-43, 'cm^3/(mol*s)'),
            n = 18.418,
            Ea = (-602.9, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -32.13,
        T3 = (207, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S4XC10H21 <=> pC4H9 + C6H12""",
)

entry(
    index = 1135,
    label = "S4XC10H21 <=> C7H14 + nC3H7",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(5.5e+11, 's^-1'), n=0.55, Ea=(28084.3, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.1e-43, 'cm^3/(mol*s)'),
            n = 18.418,
            Ea = (-602.9, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -32.13,
        T3 = (207, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S4XC10H21 <=> C7H14 + nC3H7""",
)

entry(
    index = 1136,
    label = "PXC10H21 + H <=> NC10H22",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.01e+48, 'cm^6/(mol^2*s)'),
            n = -9.32,
            Ea = (5833.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.498,
        T3 = (1314, 'K'),
        T1 = (1314, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC10H21 + H <=> NC10H22""",
)

entry(
    index = 1137,
    label = "PXC10H21 + H <=> PXC8H17 + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.7e+24, 'cm^3/(mol*s)'),
        n = -2.92,
        Ea = (12505, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is PXC10H21 + H <=> PXC8H17 + C2H5""",
)

entry(
    index = 1138,
    label = "PXC10H21 + H <=> C10H20 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC10H21 + H <=> C10H20 + H2""",
)

entry(
    index = 1139,
    label = "PXC10H21 + O <=> PXC9H19 + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC10H21 + O <=> PXC9H19 + CH2O""",
)

entry(
    index = 1140,
    label = "PXC10H21 + OH <=> C10H20 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC10H21 + OH <=> C10H20 + H2O""",
)

entry(
    index = 1141,
    label = "PXC10H21 + O2 <=> C10H20 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+10, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC10H21 + O2 <=> C10H20 + HO2""",
)

entry(
    index = 1142,
    label = "PXC10H21 + HO2 <=> PXC9H19 + OH + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC10H21 + HO2 <=> PXC9H19 + OH + CH2O""",
)

entry(
    index = 1143,
    label = "PXC10H21 + HCO <=> NC10H22 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC10H21 + HCO <=> NC10H22 + CO""",
)

entry(
    index = 1144,
    label = "PXC10H21 + CH3 <=> C10H20 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC10H21 + CH3 <=> C10H20 + CH4""",
)

entry(
    index = 1145,
    label = "SXC10H21 + H <=> NC10H22",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.7e+58, 'cm^6/(mol^2*s)'),
            n = -12.08,
            Ea = (11263.7, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.649,
        T3 = (1213.1, 'K'),
        T1 = (1213.1, 'K'),
        T2 = (13369.7, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SXC10H21 + H <=> NC10H22""",
)

entry(
    index = 1146,
    label = "SXC10H21 + H <=> PXC8H17 + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.4e+28, 'cm^3/(mol*s)'),
        n = -3.94,
        Ea = (15916, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is SXC10H21 + H <=> PXC8H17 + C2H5""",
)

entry(
    index = 1147,
    label = "SXC10H21 + H <=> C10H20 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.2e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC10H21 + H <=> C10H20 + H2""",
)

entry(
    index = 1148,
    label = "SXC10H21 + O <=> CH3CHO + PXC8H17",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC10H21 + O <=> CH3CHO + PXC8H17""",
)

entry(
    index = 1149,
    label = "SXC10H21 + OH <=> C10H20 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC10H21 + OH <=> C10H20 + H2O""",
)

entry(
    index = 1150,
    label = "SXC10H21 + O2 <=> C10H20 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC10H21 + O2 <=> C10H20 + HO2""",
)

entry(
    index = 1151,
    label = "SXC10H21 + HO2 <=> CH3CHO + PXC8H17 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC10H21 + HO2 <=> CH3CHO + PXC8H17 + OH""",
)

entry(
    index = 1152,
    label = "SXC10H21 + HCO <=> NC10H22 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC10H21 + HCO <=> NC10H22 + CO""",
)

entry(
    index = 1153,
    label = "SXC10H21 + CH3 <=> CH4 + C10H20",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.2e+14, 'cm^3/(mol*s)'), n=-0.68, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC10H21 + CH3 <=> CH4 + C10H20""",
)

entry(
    index = 1154,
    label = "S2XC10H21 + O2 <=> C10H20 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is S2XC10H21 + O2 <=> C10H20 + HO2""",
)

entry(
    index = 1155,
    label = "S3XC10H21 + O2 <=> C10H20 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is S3XC10H21 + O2 <=> C10H20 + HO2""",
)

entry(
    index = 1156,
    label = "S4XC10H21 + O2 <=> C10H20 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is S4XC10H21 + O2 <=> C10H20 + HO2""",
)

entry(
    index = 1157,
    label = "PXC9H19 + CH3 <=> NC10H22",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.93e+14, 'cm^3/(mol*s)'), n=-0.32, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC9H19 + CH3 <=> NC10H22""",
)

entry(
    index = 1158,
    label = "PXC8H17 + C2H5 <=> NC10H22",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.88e+14, 'cm^3/(mol*s)'), n=-0.5, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC8H17 + C2H5 <=> NC10H22""",
)

entry(
    index = 1159,
    label = "PXC7H15 + nC3H7 <=> NC10H22",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.88e+14, 'cm^3/(mol*s)'), n=-0.5, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC7H15 + nC3H7 <=> NC10H22""",
)

entry(
    index = 1160,
    label = "PXC6H13 + pC4H9 <=> NC10H22",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.88e+14, 'cm^3/(mol*s)'), n=-0.5, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC6H13 + pC4H9 <=> NC10H22""",
)

entry(
    index = 1161,
    label = "PXC5H11 + PXC5H11 <=> NC10H22",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.88e+14, 'cm^3/(mol*s)'), n=-0.5, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC5H11 + PXC5H11 <=> NC10H22""",
)

entry(
    index = 1162,
    label = "NC10H22 + OH <=> PXC10H21 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.73e+07, 'cm^3/(mol*s)'),
        n = 1.81,
        Ea = (868.3, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC10H22 + OH <=> PXC10H21 + H2O""",
)

entry(
    index = 1163,
    label = "NC10H22 + OH <=> SXC10H21 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.41e+10, 'cm^3/(mol*s)'),
        n = 0.94,
        Ea = (504.7, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC10H22 + OH <=> SXC10H21 + H2O""",
)

entry(
    index = 1164,
    label = "NC10H22 + OH <=> S2XC10H21 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.14e+07, 'cm^3/(mol*s)'),
        n = 1.81,
        Ea = (-1015.4, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC10H22 + OH <=> S2XC10H21 + H2O""",
)

entry(
    index = 1165,
    label = "NC10H22 + OH <=> S3XC10H21 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.12e+12, 'cm^3/(mol*s)'),
        n = 0.32,
        Ea = (846.5, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC10H22 + OH <=> S3XC10H21 + H2O""",
)

entry(
    index = 1166,
    label = "NC10H22 + OH <=> S4XC10H21 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.12e+12, 'cm^3/(mol*s)'),
        n = 0.32,
        Ea = (846.5, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC10H22 + OH <=> S4XC10H21 + H2O""",
)

entry(
    index = 1167,
    label = "NC10H22 + O2 <=> PXC10H21 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(50930, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC10H22 + O2 <=> PXC10H21 + HO2""",
)

entry(
    index = 1168,
    label = "NC10H22 + O2 <=> SXC10H21 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+13, 'cm^3/(mol*s)'), n=0, Ea=(47590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC10H22 + O2 <=> SXC10H21 + HO2""",
)

entry(
    index = 1169,
    label = "NC10H22 + O2 <=> S2XC10H21 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+13, 'cm^3/(mol*s)'), n=0, Ea=(47590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC10H22 + O2 <=> S2XC10H21 + HO2""",
)

entry(
    index = 1170,
    label = "NC10H22 + O2 <=> S3XC10H21 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+13, 'cm^3/(mol*s)'), n=0, Ea=(47590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC10H22 + O2 <=> S3XC10H21 + HO2""",
)

entry(
    index = 1171,
    label = "NC10H22 + O2 <=> S4XC10H21 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+13, 'cm^3/(mol*s)'), n=0, Ea=(47590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC10H22 + O2 <=> S4XC10H21 + HO2""",
)

entry(
    index = 1172,
    label = "NC10H22 + HO2 <=> PXC10H21 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(61100, 'cm^3/(mol*s)'), n=2.65, Ea=(17496, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC10H22 + HO2 <=> PXC10H21 + H2O2""",
)

entry(
    index = 1173,
    label = "NC10H22 + HO2 <=> SXC10H21 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(14200, 'cm^3/(mol*s)'), n=2.77, Ea=(14913, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC10H22 + HO2 <=> SXC10H21 + H2O2""",
)

entry(
    index = 1174,
    label = "NC10H22 + HO2 <=> S2XC10H21 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(14200, 'cm^3/(mol*s)'), n=2.77, Ea=(14913, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC10H22 + HO2 <=> S2XC10H21 + H2O2""",
)

entry(
    index = 1175,
    label = "NC10H22 + HO2 <=> S3XC10H21 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(14200, 'cm^3/(mol*s)'), n=2.77, Ea=(14913, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC10H22 + HO2 <=> S3XC10H21 + H2O2""",
)

entry(
    index = 1176,
    label = "NC10H22 + HO2 <=> S4XC10H21 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(14200, 'cm^3/(mol*s)'), n=2.77, Ea=(14913, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC10H22 + HO2 <=> S4XC10H21 + H2O2""",
)

entry(
    index = 1177,
    label = "NC10H22 + H <=> PXC10H21 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0645, 'cm^3/(mol*s)'), n=4.7, Ea=(3679, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC10H22 + H <=> PXC10H21 + H2""",
)

entry(
    index = 1178,
    label = "NC10H22 + H <=> SXC10H21 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0634, 'cm^3/(mol*s)'), n=4.65, Ea=(1340, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC10H22 + H <=> SXC10H21 + H2""",
)

entry(
    index = 1179,
    label = "NC10H22 + H <=> S2XC10H21 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0634, 'cm^3/(mol*s)'), n=4.65, Ea=(1340, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC10H22 + H <=> S2XC10H21 + H2""",
)

entry(
    index = 1180,
    label = "NC10H22 + H <=> S3XC10H21 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0634, 'cm^3/(mol*s)'), n=4.65, Ea=(1340, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC10H22 + H <=> S3XC10H21 + H2""",
)

entry(
    index = 1181,
    label = "NC10H22 + H <=> S4XC10H21 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0634, 'cm^3/(mol*s)'), n=4.65, Ea=(1340, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC10H22 + H <=> S4XC10H21 + H2""",
)

entry(
    index = 1182,
    label = "NC10H22 + O <=> PXC10H21 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.49, 'cm^3/(mol*s)'), n=4.17, Ea=(2766, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC10H22 + O <=> PXC10H21 + OH""",
)

entry(
    index = 1183,
    label = "NC10H22 + O <=> SXC10H21 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(36.5, 'cm^3/(mol*s)'), n=3.75, Ea=(825, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC10H22 + O <=> SXC10H21 + OH""",
)

entry(
    index = 1184,
    label = "NC10H22 + O <=> S2XC10H21 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(36.5, 'cm^3/(mol*s)'), n=3.75, Ea=(825, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC10H22 + O <=> S2XC10H21 + OH""",
)

entry(
    index = 1185,
    label = "NC10H22 + O <=> S3XC10H21 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(36.5, 'cm^3/(mol*s)'), n=3.75, Ea=(825, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC10H22 + O <=> S3XC10H21 + OH""",
)

entry(
    index = 1186,
    label = "NC10H22 + O <=> S4XC10H21 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(36.5, 'cm^3/(mol*s)'), n=3.75, Ea=(825, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC10H22 + O <=> S4XC10H21 + OH""",
)

entry(
    index = 1187,
    label = "NC10H22 + CH3 <=> PXC10H21 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.16e-08, 'cm^3/(mol*s)'),
        n = 6.08,
        Ea = (6223, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC10H22 + CH3 <=> PXC10H21 + CH4""",
)

entry(
    index = 1188,
    label = "NC10H22 + CH3 <=> SXC10H21 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.82e-07, 'cm^3/(mol*s)'),
        n = 5.89,
        Ea = (4768, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC10H22 + CH3 <=> SXC10H21 + CH4""",
)

entry(
    index = 1189,
    label = "NC10H22 + CH3 <=> S2XC10H21 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.82e-07, 'cm^3/(mol*s)'),
        n = 5.89,
        Ea = (4768, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC10H22 + CH3 <=> S2XC10H21 + CH4""",
)

entry(
    index = 1190,
    label = "NC10H22 + CH3 <=> S3XC10H21 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.82e-07, 'cm^3/(mol*s)'),
        n = 5.89,
        Ea = (4768, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC10H22 + CH3 <=> S3XC10H21 + CH4""",
)

entry(
    index = 1191,
    label = "NC10H22 + CH3 <=> S4XC10H21 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.82e-07, 'cm^3/(mol*s)'),
        n = 5.89,
        Ea = (4768, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC10H22 + CH3 <=> S4XC10H21 + CH4""",
)

entry(
    index = 1192,
    label = "PXC11H21 + H <=> C11H22",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.01e+48, 'cm^6/(mol^2*s)'),
            n = -9.32,
            Ea = (5833.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.498,
        T3 = (1314, 'K'),
        T1 = (1314, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC11H21 + H <=> C11H22""",
)

entry(
    index = 1193,
    label = "PXC11H21 + H <=> CH3 + PXC10H19",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+21, 'cm^3/(mol*s)'), n=-2, Ea=(11000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC11H21 + H <=> CH3 + PXC10H19""",
)

entry(
    index = 1194,
    label = "PXC11H21 + HO2 <=> CH2O + OH + PXC10H19",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC11H21 + HO2 <=> CH2O + OH + PXC10H19""",
)

entry(
    index = 1195,
    label = "PXC11H21 + HCO <=> C11H22 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC11H21 + HCO <=> C11H22 + CO""",
)

entry(
    index = 1196,
    label = "C2H4 + PXC9H17 <=> PXC11H21",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+11, 'cm^3/(mol*s)'), n=0, Ea=(7300, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4 + PXC9H17 <=> PXC11H21""",
)

entry(
    index = 1197,
    label = "C11H22 <=> PXC8H17 + aC3H5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.07e+23, 's^-1'), n=-2.03, Ea=(74958, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C11H22 <=> PXC8H17 + aC3H5""",
)

entry(
    index = 1198,
    label = "C11H22 <=> C8H16 + C3H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.08e+06, 's^-1'), n=1.65, Ea=(53752, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C11H22 <=> C8H16 + C3H6""",
)

entry(
    index = 1199,
    label = "C11H22 + H <=> C2H4 + PXC9H19",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8e+21, 'cm^3/(mol*s)'),
        n = -2.39,
        Ea = (11180, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C11H22 + H <=> C2H4 + PXC9H19""",
)

entry(
    index = 1200,
    label = "C11H22 + H <=> C3H6 + PXC8H17",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.6e+22, 'cm^3/(mol*s)'),
        n = -2.39,
        Ea = (11180, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C11H22 + H <=> C3H6 + PXC8H17""",
)

entry(
    index = 1201,
    label = "C11H22 + H <=> PXC11H21 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0323, 'cm^3/(mol*s)'), n=4.7, Ea=(3679, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C11H22 + H <=> PXC11H21 + H2""",
)

entry(
    index = 1202,
    label = "C11H22 + O <=> PXC10H21 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.3e+08, 'cm^3/(mol*s)'),
        n = 1.45,
        Ea = (-402, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C11H22 + O <=> PXC10H21 + HCO""",
)

entry(
    index = 1203,
    label = "C11H22 + O <=> PXC11H21 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.745, 'cm^3/(mol*s)'), n=4.17, Ea=(2766, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C11H22 + O <=> PXC11H21 + OH""",
)

entry(
    index = 1204,
    label = "C11H22 + OH <=> PXC11H21 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(700, 'cm^3/(mol*s)'), n=2.66, Ea=(527, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C11H22 + OH <=> PXC11H21 + H2O""",
)

entry(
    index = 1205,
    label = "C11H22 + O2 <=> PXC11H21 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(50930, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C11H22 + O2 <=> PXC11H21 + HO2""",
)

entry(
    index = 1206,
    label = "C11H22 + HO2 <=> PXC11H21 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(14340, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C11H22 + HO2 <=> PXC11H21 + H2O2""",
)

entry(
    index = 1207,
    label = "C11H22 + CH3 <=> PXC11H21 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.58e-08, 'cm^3/(mol*s)'),
        n = 6.08,
        Ea = (6223, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C11H22 + CH3 <=> PXC11H21 + CH4""",
)

entry(
    index = 1208,
    label = "PXC11H23 <=> C2H4 + PXC9H19",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(9.12e+11, 's^-1'), n=0.31, Ea=(27237.8, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.8e-57, 'cm^3/(mol*s)'),
            n = 23.463,
            Ea = (-602.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -2.46,
        T3 = (206, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC11H23 <=> C2H4 + PXC9H19""",
)

entry(
    index = 1209,
    label = "SXC11H23 <=> PXC8H17 + C3H6",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(6.03e+10, 's^-1'), n=0.84, Ea=(27820, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1e-43, 'cm^3/(mol*s)'),
            n = 18.591,
            Ea = (-602.5, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -43.32,
        T3 = (200, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SXC11H23 <=> PXC8H17 + C3H6""",
)

entry(
    index = 1210,
    label = "S2XC11H23 <=> PXC7H15 + C4H81",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2.04e+13, 's^-1'), n=0.04, Ea=(28493.6, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3e-43, 'cm^3/(mol*s)'),
            n = 18.43,
            Ea = (-602.8, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -34.47,
        T3 = (208, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S2XC11H23 <=> PXC7H15 + C4H81""",
)

entry(
    index = 1211,
    label = "S2XC11H23 <=> C10H20 + CH3",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(9.55e+09, 's^-1'), n=1.08, Ea=(29387.7, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.3e-46, 'cm^3/(mol*s)'),
            n = 19.133,
            Ea = (-602.7, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -34.36,
        T3 = (210, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S2XC11H23 <=> C10H20 + CH3""",
)

entry(
    index = 1212,
    label = "S3XC11H23 <=> PXC6H13 + C5H10",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(5.5e+11, 's^-1'), n=0.55, Ea=(28084.3, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.1e-43, 'cm^3/(mol*s)'),
            n = 18.418,
            Ea = (-602.9, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -32.13,
        T3 = (207, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S3XC11H23 <=> PXC6H13 + C5H10""",
)

entry(
    index = 1213,
    label = "S3XC11H23 <=> C9H18 + C2H5",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(6.76e+09, 's^-1'), n=1.11, Ea=(27023.2, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (8.2e-43, 'cm^3/(mol*s)'),
            n = 18.276,
            Ea = (-602.9, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -30.04,
        T3 = (210, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S3XC11H23 <=> C9H18 + C2H5""",
)

entry(
    index = 1214,
    label = "S4XC11H23 <=> PXC5H11 + C6H12",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(5.5e+11, 's^-1'), n=0.55, Ea=(28084.3, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.1e-43, 'cm^3/(mol*s)'),
            n = 18.418,
            Ea = (-602.9, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -32.13,
        T3 = (207, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S4XC11H23 <=> PXC5H11 + C6H12""",
)

entry(
    index = 1215,
    label = "S4XC11H23 <=> C8H16 + nC3H7",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(5.5e+11, 's^-1'), n=0.55, Ea=(28084.3, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.1e-43, 'cm^3/(mol*s)'),
            n = 18.418,
            Ea = (-602.9, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -32.13,
        T3 = (207, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S4XC11H23 <=> C8H16 + nC3H7""",
)

entry(
    index = 1216,
    label = "S5XC11H23 <=> pC4H9 + C7H14",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.1e+12, 's^-1'), n=0.55, Ea=(28084.3, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (6.2e-43, 'cm^3/(mol*s)'),
            n = 18.418,
            Ea = (-602.9, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -32.13,
        T3 = (207, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S5XC11H23 <=> pC4H9 + C7H14""",
)

entry(
    index = 1217,
    label = "PXC11H23 + H <=> NC11H24",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.01e+48, 'cm^6/(mol^2*s)'),
            n = -9.32,
            Ea = (5833.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.498,
        T3 = (1314, 'K'),
        T1 = (1314, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC11H23 + H <=> NC11H24""",
)

entry(
    index = 1218,
    label = "PXC11H23 + H <=> PXC9H19 + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.7e+24, 'cm^3/(mol*s)'),
        n = -2.92,
        Ea = (12505, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is PXC11H23 + H <=> PXC9H19 + C2H5""",
)

entry(
    index = 1219,
    label = "PXC11H23 + H <=> C11H22 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC11H23 + H <=> C11H22 + H2""",
)

entry(
    index = 1220,
    label = "PXC11H23 + O <=> PXC10H21 + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC11H23 + O <=> PXC10H21 + CH2O""",
)

entry(
    index = 1221,
    label = "PXC11H23 + OH <=> C11H22 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC11H23 + OH <=> C11H22 + H2O""",
)

entry(
    index = 1222,
    label = "PXC11H23 + O2 <=> C11H22 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+10, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC11H23 + O2 <=> C11H22 + HO2""",
)

entry(
    index = 1223,
    label = "PXC11H23 + HO2 <=> PXC10H21 + OH + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC11H23 + HO2 <=> PXC10H21 + OH + CH2O""",
)

entry(
    index = 1224,
    label = "PXC11H23 + HCO <=> NC11H24 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC11H23 + HCO <=> NC11H24 + CO""",
)

entry(
    index = 1225,
    label = "PXC11H23 + CH3 <=> C11H22 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC11H23 + CH3 <=> C11H22 + CH4""",
)

entry(
    index = 1226,
    label = "SXC11H23 + H <=> NC11H24",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.7e+58, 'cm^6/(mol^2*s)'),
            n = -12.08,
            Ea = (11263.7, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.649,
        T3 = (1213.1, 'K'),
        T1 = (1213.1, 'K'),
        T2 = (13369.7, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SXC11H23 + H <=> NC11H24""",
)

entry(
    index = 1227,
    label = "SXC11H23 + H <=> PXC9H19 + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.4e+28, 'cm^3/(mol*s)'),
        n = -3.94,
        Ea = (15916, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is SXC11H23 + H <=> PXC9H19 + C2H5""",
)

entry(
    index = 1228,
    label = "SXC11H23 + H <=> C11H22 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.2e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC11H23 + H <=> C11H22 + H2""",
)

entry(
    index = 1229,
    label = "SXC11H23 + O <=> CH3CHO + PXC9H19",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC11H23 + O <=> CH3CHO + PXC9H19""",
)

entry(
    index = 1230,
    label = "SXC11H23 + OH <=> C11H22 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC11H23 + OH <=> C11H22 + H2O""",
)

entry(
    index = 1231,
    label = "SXC11H23 + O2 <=> C11H22 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC11H23 + O2 <=> C11H22 + HO2""",
)

entry(
    index = 1232,
    label = "SXC11H23 + HO2 <=> CH3CHO + PXC9H19 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC11H23 + HO2 <=> CH3CHO + PXC9H19 + OH""",
)

entry(
    index = 1233,
    label = "SXC11H23 + HCO <=> NC11H24 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC11H23 + HCO <=> NC11H24 + CO""",
)

entry(
    index = 1234,
    label = "SXC11H23 + CH3 <=> CH4 + C11H22",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.2e+14, 'cm^3/(mol*s)'), n=-0.68, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC11H23 + CH3 <=> CH4 + C11H22""",
)

entry(
    index = 1235,
    label = "S2XC11H23 + O2 <=> C11H22 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is S2XC11H23 + O2 <=> C11H22 + HO2""",
)

entry(
    index = 1236,
    label = "S3XC11H23 + O2 <=> C11H22 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is S3XC11H23 + O2 <=> C11H22 + HO2""",
)

entry(
    index = 1237,
    label = "S4XC11H23 + O2 <=> C11H22 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is S4XC11H23 + O2 <=> C11H22 + HO2""",
)

entry(
    index = 1238,
    label = "S5XC11H23 + O2 <=> C11H22 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is S5XC11H23 + O2 <=> C11H22 + HO2""",
)

entry(
    index = 1239,
    label = "PXC10H21 + CH3 <=> NC11H24",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.93e+14, 'cm^3/(mol*s)'), n=-0.32, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC10H21 + CH3 <=> NC11H24""",
)

entry(
    index = 1240,
    label = "PXC9H19 + C2H5 <=> NC11H24",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.88e+14, 'cm^3/(mol*s)'), n=-0.5, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC9H19 + C2H5 <=> NC11H24""",
)

entry(
    index = 1241,
    label = "PXC8H17 + nC3H7 <=> NC11H24",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.88e+14, 'cm^3/(mol*s)'), n=-0.5, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC8H17 + nC3H7 <=> NC11H24""",
)

entry(
    index = 1242,
    label = "PXC7H15 + pC4H9 <=> NC11H24",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.88e+14, 'cm^3/(mol*s)'), n=-0.5, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC7H15 + pC4H9 <=> NC11H24""",
)

entry(
    index = 1243,
    label = "PXC6H13 + PXC5H11 <=> NC11H24",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.88e+14, 'cm^3/(mol*s)'), n=-0.5, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC6H13 + PXC5H11 <=> NC11H24""",
)

entry(
    index = 1244,
    label = "NC11H24 + OH <=> PXC11H23 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.73e+07, 'cm^3/(mol*s)'),
        n = 1.81,
        Ea = (868.3, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC11H24 + OH <=> PXC11H23 + H2O""",
)

entry(
    index = 1245,
    label = "NC11H24 + OH <=> SXC11H23 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.41e+10, 'cm^3/(mol*s)'),
        n = 0.94,
        Ea = (504.7, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC11H24 + OH <=> SXC11H23 + H2O""",
)

entry(
    index = 1246,
    label = "NC11H24 + OH <=> S2XC11H23 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.14e+07, 'cm^3/(mol*s)'),
        n = 1.81,
        Ea = (-1015.4, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC11H24 + OH <=> S2XC11H23 + H2O""",
)

entry(
    index = 1247,
    label = "NC11H24 + OH <=> S3XC11H23 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.12e+12, 'cm^3/(mol*s)'),
        n = 0.32,
        Ea = (846.5, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC11H24 + OH <=> S3XC11H23 + H2O""",
)

entry(
    index = 1248,
    label = "NC11H24 + OH <=> S4XC11H23 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.12e+12, 'cm^3/(mol*s)'),
        n = 0.32,
        Ea = (846.5, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC11H24 + OH <=> S4XC11H23 + H2O""",
)

entry(
    index = 1249,
    label = "NC11H24 + OH <=> S5XC11H23 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.62e+11, 'cm^3/(mol*s)'),
        n = 0.32,
        Ea = (846.5, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC11H24 + OH <=> S5XC11H23 + H2O""",
)

entry(
    index = 1250,
    label = "NC11H24 + O2 <=> PXC11H23 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(50930, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC11H24 + O2 <=> PXC11H23 + HO2""",
)

entry(
    index = 1251,
    label = "NC11H24 + O2 <=> SXC11H23 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+13, 'cm^3/(mol*s)'), n=0, Ea=(47590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC11H24 + O2 <=> SXC11H23 + HO2""",
)

entry(
    index = 1252,
    label = "NC11H24 + O2 <=> S2XC11H23 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+13, 'cm^3/(mol*s)'), n=0, Ea=(47590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC11H24 + O2 <=> S2XC11H23 + HO2""",
)

entry(
    index = 1253,
    label = "NC11H24 + O2 <=> S3XC11H23 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+13, 'cm^3/(mol*s)'), n=0, Ea=(47590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC11H24 + O2 <=> S3XC11H23 + HO2""",
)

entry(
    index = 1254,
    label = "NC11H24 + O2 <=> S4XC11H23 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+13, 'cm^3/(mol*s)'), n=0, Ea=(47590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC11H24 + O2 <=> S4XC11H23 + HO2""",
)

entry(
    index = 1255,
    label = "NC11H24 + O2 <=> S5XC11H23 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(47590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC11H24 + O2 <=> S5XC11H23 + HO2""",
)

entry(
    index = 1256,
    label = "NC11H24 + HO2 <=> PXC11H23 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(61100, 'cm^3/(mol*s)'), n=2.65, Ea=(17496, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC11H24 + HO2 <=> PXC11H23 + H2O2""",
)

entry(
    index = 1257,
    label = "NC11H24 + HO2 <=> SXC11H23 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(14200, 'cm^3/(mol*s)'), n=2.77, Ea=(14913, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC11H24 + HO2 <=> SXC11H23 + H2O2""",
)

entry(
    index = 1258,
    label = "NC11H24 + HO2 <=> S2XC11H23 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(14200, 'cm^3/(mol*s)'), n=2.77, Ea=(14913, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC11H24 + HO2 <=> S2XC11H23 + H2O2""",
)

entry(
    index = 1259,
    label = "NC11H24 + HO2 <=> S3XC11H23 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(14200, 'cm^3/(mol*s)'), n=2.77, Ea=(14913, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC11H24 + HO2 <=> S3XC11H23 + H2O2""",
)

entry(
    index = 1260,
    label = "NC11H24 + HO2 <=> S4XC11H23 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(14200, 'cm^3/(mol*s)'), n=2.77, Ea=(14913, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC11H24 + HO2 <=> S4XC11H23 + H2O2""",
)

entry(
    index = 1261,
    label = "NC11H24 + HO2 <=> S5XC11H23 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7130, 'cm^3/(mol*s)'), n=2.77, Ea=(14913, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC11H24 + HO2 <=> S5XC11H23 + H2O2""",
)

entry(
    index = 1262,
    label = "NC11H24 + H <=> PXC11H23 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0645, 'cm^3/(mol*s)'), n=4.7, Ea=(3679, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC11H24 + H <=> PXC11H23 + H2""",
)

entry(
    index = 1263,
    label = "NC11H24 + H <=> SXC11H23 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0634, 'cm^3/(mol*s)'), n=4.65, Ea=(1340, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC11H24 + H <=> SXC11H23 + H2""",
)

entry(
    index = 1264,
    label = "NC11H24 + H <=> S2XC11H23 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0634, 'cm^3/(mol*s)'), n=4.65, Ea=(1340, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC11H24 + H <=> S2XC11H23 + H2""",
)

entry(
    index = 1265,
    label = "NC11H24 + H <=> S3XC11H23 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0634, 'cm^3/(mol*s)'), n=4.65, Ea=(1340, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC11H24 + H <=> S3XC11H23 + H2""",
)

entry(
    index = 1266,
    label = "NC11H24 + H <=> S4XC11H23 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0634, 'cm^3/(mol*s)'), n=4.65, Ea=(1340, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC11H24 + H <=> S4XC11H23 + H2""",
)

entry(
    index = 1267,
    label = "NC11H24 + H <=> S5XC11H23 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0317, 'cm^3/(mol*s)'), n=4.65, Ea=(1340, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC11H24 + H <=> S5XC11H23 + H2""",
)

entry(
    index = 1268,
    label = "NC11H24 + O <=> PXC11H23 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.49, 'cm^3/(mol*s)'), n=4.17, Ea=(2766, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC11H24 + O <=> PXC11H23 + OH""",
)

entry(
    index = 1269,
    label = "NC11H24 + O <=> SXC11H23 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(36.5, 'cm^3/(mol*s)'), n=3.75, Ea=(825, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC11H24 + O <=> SXC11H23 + OH""",
)

entry(
    index = 1270,
    label = "NC11H24 + O <=> S2XC11H23 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(36.5, 'cm^3/(mol*s)'), n=3.75, Ea=(825, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC11H24 + O <=> S2XC11H23 + OH""",
)

entry(
    index = 1271,
    label = "NC11H24 + O <=> S3XC11H23 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(36.5, 'cm^3/(mol*s)'), n=3.75, Ea=(825, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC11H24 + O <=> S3XC11H23 + OH""",
)

entry(
    index = 1272,
    label = "NC11H24 + O <=> S4XC11H23 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(36.5, 'cm^3/(mol*s)'), n=3.75, Ea=(825, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC11H24 + O <=> S4XC11H23 + OH""",
)

entry(
    index = 1273,
    label = "NC11H24 + O <=> S5XC11H23 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(18.2, 'cm^3/(mol*s)'), n=3.75, Ea=(825, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC11H24 + O <=> S5XC11H23 + OH""",
)

entry(
    index = 1274,
    label = "NC11H24 + CH3 <=> PXC11H23 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.16e-08, 'cm^3/(mol*s)'),
        n = 6.08,
        Ea = (6223, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC11H24 + CH3 <=> PXC11H23 + CH4""",
)

entry(
    index = 1275,
    label = "NC11H24 + CH3 <=> SXC11H23 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.82e-07, 'cm^3/(mol*s)'),
        n = 5.89,
        Ea = (4768, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC11H24 + CH3 <=> SXC11H23 + CH4""",
)

entry(
    index = 1276,
    label = "NC11H24 + CH3 <=> S2XC11H23 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.82e-07, 'cm^3/(mol*s)'),
        n = 5.89,
        Ea = (4768, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC11H24 + CH3 <=> S2XC11H23 + CH4""",
)

entry(
    index = 1277,
    label = "NC11H24 + CH3 <=> S3XC11H23 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.82e-07, 'cm^3/(mol*s)'),
        n = 5.89,
        Ea = (4768, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC11H24 + CH3 <=> S3XC11H23 + CH4""",
)

entry(
    index = 1278,
    label = "NC11H24 + CH3 <=> S4XC11H23 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.82e-07, 'cm^3/(mol*s)'),
        n = 5.89,
        Ea = (4768, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC11H24 + CH3 <=> S4XC11H23 + CH4""",
)

entry(
    index = 1279,
    label = "NC11H24 + CH3 <=> S5XC11H23 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.41e-07, 'cm^3/(mol*s)'),
        n = 5.89,
        Ea = (4768, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC11H24 + CH3 <=> S5XC11H23 + CH4""",
)

entry(
    index = 1280,
    label = "PXC12H23 + H <=> C12H24",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.01e+48, 'cm^6/(mol^2*s)'),
            n = -9.32,
            Ea = (5833.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.498,
        T3 = (1314, 'K'),
        T1 = (1314, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC12H23 + H <=> C12H24""",
)

entry(
    index = 1281,
    label = "PXC12H23 + H <=> CH3 + PXC11H21",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+21, 'cm^3/(mol*s)'), n=-2, Ea=(11000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC12H23 + H <=> CH3 + PXC11H21""",
)

entry(
    index = 1282,
    label = "PXC12H23 + HO2 <=> CH2O + OH + PXC11H21",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC12H23 + HO2 <=> CH2O + OH + PXC11H21""",
)

entry(
    index = 1283,
    label = "PXC12H23 + HCO <=> C12H24 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC12H23 + HCO <=> C12H24 + CO""",
)

entry(
    index = 1284,
    label = "C2H4 + PXC10H19 <=> PXC12H23",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+11, 'cm^3/(mol*s)'), n=0, Ea=(7300, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4 + PXC10H19 <=> PXC12H23""",
)

entry(
    index = 1285,
    label = "C12H24 <=> PXC9H19 + aC3H5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.07e+23, 's^-1'), n=-2.03, Ea=(74958, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C12H24 <=> PXC9H19 + aC3H5""",
)

entry(
    index = 1286,
    label = "C12H24 <=> C9H18 + C3H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.08e+06, 's^-1'), n=1.65, Ea=(53752, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C12H24 <=> C9H18 + C3H6""",
)

entry(
    index = 1287,
    label = "C12H24 + H <=> C2H4 + PXC10H21",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8e+21, 'cm^3/(mol*s)'),
        n = -2.39,
        Ea = (11180, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C12H24 + H <=> C2H4 + PXC10H21""",
)

entry(
    index = 1288,
    label = "C12H24 + H <=> C3H6 + PXC9H19",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.6e+22, 'cm^3/(mol*s)'),
        n = -2.39,
        Ea = (11180, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C12H24 + H <=> C3H6 + PXC9H19""",
)

entry(
    index = 1289,
    label = "C12H24 + H <=> PXC12H23 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0323, 'cm^3/(mol*s)'), n=4.7, Ea=(3679, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C12H24 + H <=> PXC12H23 + H2""",
)

entry(
    index = 1290,
    label = "C12H24 + O <=> PXC11H23 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.3e+08, 'cm^3/(mol*s)'),
        n = 1.45,
        Ea = (-402, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C12H24 + O <=> PXC11H23 + HCO""",
)

entry(
    index = 1291,
    label = "C12H24 + O <=> PXC12H23 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.745, 'cm^3/(mol*s)'), n=4.17, Ea=(2766, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C12H24 + O <=> PXC12H23 + OH""",
)

entry(
    index = 1292,
    label = "C12H24 + OH <=> PXC12H23 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(700, 'cm^3/(mol*s)'), n=2.66, Ea=(527, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C12H24 + OH <=> PXC12H23 + H2O""",
)

entry(
    index = 1293,
    label = "C12H24 + O2 <=> PXC12H23 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(50930, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C12H24 + O2 <=> PXC12H23 + HO2""",
)

entry(
    index = 1294,
    label = "C12H24 + HO2 <=> PXC12H23 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(14340, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C12H24 + HO2 <=> PXC12H23 + H2O2""",
)

entry(
    index = 1295,
    label = "C12H24 + CH3 <=> PXC12H23 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.58e-08, 'cm^3/(mol*s)'),
        n = 6.08,
        Ea = (6223, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C12H24 + CH3 <=> PXC12H23 + CH4""",
)

entry(
    index = 1296,
    label = "PXC12H25 <=> C2H4 + PXC10H21",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(9.12e+11, 's^-1'), n=0.31, Ea=(27237.8, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.8e-57, 'cm^3/(mol*s)'),
            n = 23.463,
            Ea = (-602.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -2.46,
        T3 = (206, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC12H25 <=> C2H4 + PXC10H21""",
)

entry(
    index = 1297,
    label = "SXC12H25 <=> C3H6 + PXC9H19",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(6.03e+10, 's^-1'), n=0.84, Ea=(27820, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1e-43, 'cm^3/(mol*s)'),
            n = 18.591,
            Ea = (-602.5, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -43.32,
        T3 = (200, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SXC12H25 <=> C3H6 + PXC9H19""",
)

entry(
    index = 1298,
    label = "S2XC12H25 <=> C4H81 + PXC8H17",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2.04e+13, 's^-1'), n=0.04, Ea=(28493.6, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3e-43, 'cm^3/(mol*s)'),
            n = 18.43,
            Ea = (-602.8, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -34.47,
        T3 = (208, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S2XC12H25 <=> C4H81 + PXC8H17""",
)

entry(
    index = 1299,
    label = "S2XC12H25 <=> C11H22 + CH3",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(9.55e+09, 's^-1'), n=1.08, Ea=(29387.7, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.3e-46, 'cm^3/(mol*s)'),
            n = 19.133,
            Ea = (-602.7, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -34.36,
        T3 = (210, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S2XC12H25 <=> C11H22 + CH3""",
)

entry(
    index = 1300,
    label = "S3XC12H25 <=> C5H10 + PXC7H15",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(5.5e+11, 's^-1'), n=0.55, Ea=(28084.3, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.1e-43, 'cm^3/(mol*s)'),
            n = 18.418,
            Ea = (-602.9, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -32.13,
        T3 = (207, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S3XC12H25 <=> C5H10 + PXC7H15""",
)

entry(
    index = 1301,
    label = "S3XC12H25 <=> C10H20 + C2H5",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(6.76e+09, 's^-1'), n=1.11, Ea=(27023.2, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (8.2e-43, 'cm^3/(mol*s)'),
            n = 18.276,
            Ea = (-602.9, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -30.04,
        T3 = (210, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S3XC12H25 <=> C10H20 + C2H5""",
)

entry(
    index = 1302,
    label = "S4XC12H25 <=> C6H12 + PXC6H13",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(5.5e+11, 's^-1'), n=0.55, Ea=(28084.3, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.1e-43, 'cm^3/(mol*s)'),
            n = 18.418,
            Ea = (-602.9, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -32.13,
        T3 = (207, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S4XC12H25 <=> C6H12 + PXC6H13""",
)

entry(
    index = 1303,
    label = "S4XC12H25 <=> C9H18 + nC3H7",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(5.5e+11, 's^-1'), n=0.55, Ea=(28084.3, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.1e-43, 'cm^3/(mol*s)'),
            n = 18.418,
            Ea = (-602.9, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -32.13,
        T3 = (207, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S4XC12H25 <=> C9H18 + nC3H7""",
)

entry(
    index = 1304,
    label = "S5XC12H25 <=> C7H14 + PXC5H11",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(5.5e+11, 's^-1'), n=0.55, Ea=(28084.3, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.1e-43, 'cm^3/(mol*s)'),
            n = 18.418,
            Ea = (-602.9, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -32.13,
        T3 = (207, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S5XC12H25 <=> C7H14 + PXC5H11""",
)

entry(
    index = 1305,
    label = "S5XC12H25 <=> C8H16 + pC4H9",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(5.5e+11, 's^-1'), n=0.55, Ea=(28084.3, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.1e-43, 'cm^3/(mol*s)'),
            n = 18.418,
            Ea = (-602.9, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -32.13,
        T3 = (207, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S5XC12H25 <=> C8H16 + pC4H9""",
)

entry(
    index = 1306,
    label = "PXC12H25 + H <=> NC12H26",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.01e+48, 'cm^6/(mol^2*s)'),
            n = -9.32,
            Ea = (5833.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.498,
        T3 = (1314, 'K'),
        T1 = (1314, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC12H25 + H <=> NC12H26""",
)

entry(
    index = 1307,
    label = "SXC12H25 + H <=> PXC10H21 + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.4e+28, 'cm^3/(mol*s)'),
        n = -3.94,
        Ea = (15916, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is SXC12H25 + H <=> PXC10H21 + C2H5""",
)

entry(
    index = 1308,
    label = "SXC12H25 + H <=> C12H24 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.2e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC12H25 + H <=> C12H24 + H2""",
)

entry(
    index = 1309,
    label = "SXC12H25 + O <=> CH3CHO + PXC10H21",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC12H25 + O <=> CH3CHO + PXC10H21""",
)

entry(
    index = 1310,
    label = "SXC12H25 + OH <=> C12H24 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC12H25 + OH <=> C12H24 + H2O""",
)

entry(
    index = 1311,
    label = "SXC12H25 + O2 <=> C12H24 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC12H25 + O2 <=> C12H24 + HO2""",
)

entry(
    index = 1312,
    label = "SXC12H25 + HO2 <=> CH3CHO + PXC10H21 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC12H25 + HO2 <=> CH3CHO + PXC10H21 + OH""",
)

entry(
    index = 1313,
    label = "SXC12H25 + HCO <=> NC12H26 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC12H25 + HCO <=> NC12H26 + CO""",
)

entry(
    index = 1314,
    label = "SXC12H25 + CH3 <=> CH4 + C12H24",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.2e+14, 'cm^3/(mol*s)'), n=-0.68, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC12H25 + CH3 <=> CH4 + C12H24""",
)

entry(
    index = 1315,
    label = "S2XC12H25 + O2 <=> C12H24 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is S2XC12H25 + O2 <=> C12H24 + HO2""",
)

entry(
    index = 1316,
    label = "S3XC12H25 + O2 <=> C12H24 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is S3XC12H25 + O2 <=> C12H24 + HO2""",
)

entry(
    index = 1317,
    label = "S4XC12H25 + O2 <=> C12H24 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is S4XC12H25 + O2 <=> C12H24 + HO2""",
)

entry(
    index = 1318,
    label = "S5XC12H25 + O2 <=> C12H24 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is S5XC12H25 + O2 <=> C12H24 + HO2""",
)

entry(
    index = 1319,
    label = "PXC12H25 <=> S3XC12H25",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(5.13, 's^-1'), n=3.23, Ea=(16847.8, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.1e-44, 'cm^3/(mol*s)'),
            n = 18.749,
            Ea = (-602.9, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -20.15,
        T3 = (205, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC12H25 <=> S3XC12H25""",
)

entry(
    index = 1320,
    label = "PXC11H23 <=> S3XC11H23",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(5.13, 's^-1'), n=3.23, Ea=(16847.8, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.1e-44, 'cm^3/(mol*s)'),
            n = 18.749,
            Ea = (-602.9, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -20.15,
        T3 = (205, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC11H23 <=> S3XC11H23""",
)

entry(
    index = 1321,
    label = "PXC10H21 <=> S3XC10H21",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(5.13, 's^-1'), n=3.23, Ea=(16847.8, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.1e-44, 'cm^3/(mol*s)'),
            n = 18.749,
            Ea = (-602.9, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -20.15,
        T3 = (205, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC10H21 <=> S3XC10H21""",
)

entry(
    index = 1322,
    label = "PXC9H19 <=> S3XC9H19",
    degeneracy = 1,
    duplicate = True,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(5.13, 's^-1'), n=3.23, Ea=(16847.8, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.1e-44, 'cm^3/(mol*s)'),
            n = 18.749,
            Ea = (-602.9, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -20.15,
        T3 = (205, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC9H19 <=> S3XC9H19""",
)

entry(
    index = 1323,
    label = "PXC8H17 <=> S3XC8H17",
    degeneracy = 1,
    duplicate = True,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(5.13, 's^-1'), n=3.23, Ea=(16847.8, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.1e-44, 'cm^3/(mol*s)'),
            n = 18.749,
            Ea = (-602.9, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -20.15,
        T3 = (205, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC8H17 <=> S3XC8H17""",
)

entry(
    index = 1324,
    label = "PXC7H15 <=> S3XC7H15",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(117, 's^-1'), n=2.85, Ea=(17247.2, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.9e-31, 'cm^3/(mol*s)'),
            n = 14.521,
            Ea = (-599, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -8.1,
        T3 = (259, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC7H15 <=> S3XC7H15""",
)

entry(
    index = 1325,
    label = "PXC6H13 <=> S2XC6H13",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(6.98, 's^-1'), n=3.2, Ea=(16557.7, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (2e-26, 'cm^3/(mol*s)'),
            n = 12.833,
            Ea = (-600.7, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -10.14,
        T3 = (307, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC6H13 <=> S2XC6H13""",
)

entry(
    index = 1326,
    label = "PXC5H11 <=> SXC5H11",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1e+12, 's^-1'), n=0, Ea=(22453.1, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (2e-26, 'cm^3/(mol*s)'),
            n = 12.833,
            Ea = (-600.7, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -10.14,
        T3 = (307, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC5H11 <=> SXC5H11""",
)

entry(
    index = 1327,
    label = "PXC12H25 <=> S4XC12H25",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(22.9, 's^-1'), n=2.82, Ea=(10755.6, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (9.9e-38, 'cm^3/(mol*s)'),
            n = 17.215,
            Ea = (-603, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -16.33,
        T3 = (200, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC12H25 <=> S4XC12H25""",
)

entry(
    index = 1328,
    label = "PXC11H23 <=> S4XC11H23",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(22.9, 's^-1'), n=2.82, Ea=(10755.6, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (9.9e-38, 'cm^3/(mol*s)'),
            n = 17.215,
            Ea = (-603, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -16.33,
        T3 = (200, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC11H23 <=> S4XC11H23""",
)

entry(
    index = 1329,
    label = "PXC10H21 <=> S4XC10H21",
    degeneracy = 1,
    duplicate = True,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(22.9, 's^-1'), n=2.82, Ea=(10755.6, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (9.9e-38, 'cm^3/(mol*s)'),
            n = 17.215,
            Ea = (-603, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -16.33,
        T3 = (200, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC10H21 <=> S4XC10H21""",
)

entry(
    index = 1330,
    label = "PXC9H19 <=> S4XC9H19",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(22.9, 's^-1'), n=2.82, Ea=(10755.6, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (9.9e-38, 'cm^3/(mol*s)'),
            n = 17.215,
            Ea = (-603, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -16.33,
        T3 = (200, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC9H19 <=> S4XC9H19""",
)

entry(
    index = 1331,
    label = "PXC8H17 <=> S3XC8H17",
    degeneracy = 1,
    duplicate = True,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(22.9, 's^-1'), n=2.82, Ea=(10755.6, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (9.9e-38, 'cm^3/(mol*s)'),
            n = 17.215,
            Ea = (-603, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -16.33,
        T3 = (200, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC8H17 <=> S3XC8H17""",
)

entry(
    index = 1332,
    label = "PXC7H15 <=> S2XC7H15",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(676, 's^-1'), n=2.39, Ea=(10405.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.4e-26, 'cm^3/(mol*s)'),
            n = 13.202,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -25.39,
        T3 = (215, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC7H15 <=> S2XC7H15""",
)

entry(
    index = 1333,
    label = "PXC6H13 <=> SXC6H13",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(183, 's^-1'), n=2.55, Ea=(10960.3, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.4e-26, 'cm^3/(mol*s)'),
            n = 13.087,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -25.38,
        T3 = (215, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC6H13 <=> SXC6H13""",
)

entry(
    index = 1334,
    label = "PXC12H25 <=> S5XC12H25",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2.95, 's^-1'), n=3.08, Ea=(11015.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.9e-34, 'cm^3/(mol*s)'),
            n = 15.855,
            Ea = (-606.2, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -15.24,
        T3 = (216, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC12H25 <=> S5XC12H25""",
)

entry(
    index = 1335,
    label = "PXC11H23 <=> S5XC11H23",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2.95, 's^-1'), n=3.08, Ea=(11015.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.9e-34, 'cm^3/(mol*s)'),
            n = 15.855,
            Ea = (-606.2, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -15.24,
        T3 = (216, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC11H23 <=> S5XC11H23""",
)

entry(
    index = 1336,
    label = "PXC10H21 <=> S4XC10H21",
    degeneracy = 1,
    duplicate = True,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2.95, 's^-1'), n=3.08, Ea=(11015.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.9e-34, 'cm^3/(mol*s)'),
            n = 15.855,
            Ea = (-606.2, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -15.24,
        T3 = (216, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC10H21 <=> S4XC10H21""",
)

entry(
    index = 1337,
    label = "PXC9H19 <=> S3XC9H19",
    degeneracy = 1,
    duplicate = True,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2.95, 's^-1'), n=3.08, Ea=(11015.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.9e-34, 'cm^3/(mol*s)'),
            n = 15.855,
            Ea = (-606.2, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -15.24,
        T3 = (216, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC9H19 <=> S3XC9H19""",
)

entry(
    index = 1338,
    label = "PXC8H17 <=> S2XC8H17",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2.95, 's^-1'), n=3.08, Ea=(11015.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.9e-34, 'cm^3/(mol*s)'),
            n = 15.855,
            Ea = (-606.2, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -15.24,
        T3 = (216, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC8H17 <=> S2XC8H17""",
)

entry(
    index = 1339,
    label = "PXC7H15 <=> SXC7H15",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(245, 's^-1'), n=2.51, Ea=(12502.2, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.3e-30, 'cm^3/(mol*s)'),
            n = 14.309,
            Ea = (-602.5, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -27.19,
        T3 = (220, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC7H15 <=> SXC7H15""",
)

entry(
    index = 1340,
    label = "S3XC12H25 <=> S5XC12H25",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.41, 's^-1'), n=3.32, Ea=(16144.4, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.2e-30, 'cm^3/(mol*s)'),
            n = 14.079,
            Ea = (-606.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -21.93,
        T3 = (219, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S3XC12H25 <=> S5XC12H25""",
)

entry(
    index = 1341,
    label = "S2XC12H25 <=> S5XC12H25",
    degeneracy = 1,
    duplicate = True,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.41, 's^-1'), n=3.32, Ea=(16144.4, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.2e-30, 'cm^3/(mol*s)'),
            n = 14.079,
            Ea = (-606.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -21.93,
        T3 = (219, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S2XC12H25 <=> S5XC12H25""",
)

entry(
    index = 1342,
    label = "SXC12H25 <=> S4XC12H25",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.41, 's^-1'), n=3.32, Ea=(16144.4, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.2e-30, 'cm^3/(mol*s)'),
            n = 14.079,
            Ea = (-606.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -21.93,
        T3 = (219, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SXC12H25 <=> S4XC12H25""",
)

entry(
    index = 1343,
    label = "S3XC11H23 <=> S4XC11H23",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.41, 's^-1'), n=3.32, Ea=(16144.4, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.2e-30, 'cm^3/(mol*s)'),
            n = 14.079,
            Ea = (-606.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -21.93,
        T3 = (219, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S3XC11H23 <=> S4XC11H23""",
)

entry(
    index = 1344,
    label = "S2XC11H23 <=> S5XC11H23",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.41, 's^-1'), n=3.32, Ea=(16144.4, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.2e-30, 'cm^3/(mol*s)'),
            n = 14.079,
            Ea = (-606.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -21.93,
        T3 = (219, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S2XC11H23 <=> S5XC11H23""",
)

entry(
    index = 1345,
    label = "SXC11H23 <=> S4XC11H23",
    degeneracy = 1,
    duplicate = True,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.41, 's^-1'), n=3.32, Ea=(16144.4, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.2e-30, 'cm^3/(mol*s)'),
            n = 14.079,
            Ea = (-606.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -21.93,
        T3 = (219, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SXC11H23 <=> S4XC11H23""",
)

entry(
    index = 1346,
    label = "S4XC10H21 <=> S2XC10H21",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.41, 's^-1'), n=3.32, Ea=(16144.4, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.2e-30, 'cm^3/(mol*s)'),
            n = 14.079,
            Ea = (-606.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -21.93,
        T3 = (219, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S4XC10H21 <=> S2XC10H21""",
)

entry(
    index = 1347,
    label = "S4XC10H21 <=> SXC10H21",
    degeneracy = 1,
    duplicate = True,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.41, 's^-1'), n=3.32, Ea=(16144.4, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.2e-30, 'cm^3/(mol*s)'),
            n = 14.079,
            Ea = (-606.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -21.93,
        T3 = (219, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S4XC10H21 <=> SXC10H21""",
)

entry(
    index = 1348,
    label = "SXC9H19 <=> S4XC9H19",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.41, 's^-1'), n=3.32, Ea=(16144.4, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.2e-30, 'cm^3/(mol*s)'),
            n = 14.079,
            Ea = (-606.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -21.93,
        T3 = (219, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SXC9H19 <=> S4XC9H19""",
)

entry(
    index = 1349,
    label = "S3XC9H19 <=> S2XC9H19",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.41, 's^-1'), n=3.32, Ea=(16144.4, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.2e-30, 'cm^3/(mol*s)'),
            n = 14.079,
            Ea = (-606.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -21.93,
        T3 = (219, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S3XC9H19 <=> S2XC9H19""",
)

entry(
    index = 1350,
    label = "SXC8H17 <=> S3XC8H17",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.41, 's^-1'), n=3.32, Ea=(16144.4, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.2e-30, 'cm^3/(mol*s)'),
            n = 14.079,
            Ea = (-606.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -21.93,
        T3 = (219, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SXC8H17 <=> S3XC8H17""",
)

entry(
    index = 1351,
    label = "SXC7H15 <=> S2XC7H15",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(24.5, 's^-1'), n=3.09, Ea=(18107.5, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.3e-32, 'cm^3/(mol*s)'),
            n = 14.834,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -28.41,
        T3 = (219, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SXC7H15 <=> S2XC7H15""",
)

entry(
    index = 1352,
    label = "SXC12H25 <=> S5XC12H25",
    degeneracy = 1,
    duplicate = True,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.86, 's^-1'), n=3.27, Ea=(13197.7, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3e-27, 'cm^3/(mol*s)'),
            n = 13.481,
            Ea = (-606.5, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -19.7,
        T3 = (215, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SXC12H25 <=> S5XC12H25""",
)

entry(
    index = 1353,
    label = "S2XC12H25 <=> S5XC12H25",
    degeneracy = 1,
    duplicate = True,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.86, 's^-1'), n=3.27, Ea=(13197.7, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3e-27, 'cm^3/(mol*s)'),
            n = 13.481,
            Ea = (-606.5, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -19.7,
        T3 = (215, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S2XC12H25 <=> S5XC12H25""",
)

entry(
    index = 1354,
    label = "S3XC12H25 <=> S4XC12H25",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.86, 's^-1'), n=3.27, Ea=(13197.7, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3e-27, 'cm^3/(mol*s)'),
            n = 13.481,
            Ea = (-606.5, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -19.7,
        T3 = (215, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S3XC12H25 <=> S4XC12H25""",
)

entry(
    index = 1355,
    label = "SXC11H23 <=> S5XC11H23",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.86, 's^-1'), n=3.27, Ea=(13197.7, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3e-27, 'cm^3/(mol*s)'),
            n = 13.481,
            Ea = (-606.5, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -19.7,
        T3 = (215, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SXC11H23 <=> S5XC11H23""",
)

entry(
    index = 1356,
    label = "S2XC11H23 <=> S4XC11H23",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.86, 's^-1'), n=3.27, Ea=(13197.7, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3e-27, 'cm^3/(mol*s)'),
            n = 13.481,
            Ea = (-606.5, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -19.7,
        T3 = (215, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S2XC11H23 <=> S4XC11H23""",
)

entry(
    index = 1357,
    label = "S2XC10H21 <=> S3XC10H21",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.86, 's^-1'), n=3.27, Ea=(13197.7, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3e-27, 'cm^3/(mol*s)'),
            n = 13.481,
            Ea = (-606.5, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -19.7,
        T3 = (215, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S2XC10H21 <=> S3XC10H21""",
)

entry(
    index = 1358,
    label = "SXC10H21 <=> S4XC10H21",
    degeneracy = 1,
    duplicate = True,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.86, 's^-1'), n=3.27, Ea=(13197.7, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3e-27, 'cm^3/(mol*s)'),
            n = 13.481,
            Ea = (-606.5, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -19.7,
        T3 = (215, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SXC10H21 <=> S4XC10H21""",
)

entry(
    index = 1359,
    label = "SXC9H19 <=> S2XC9H19",
    degeneracy = 1,
    duplicate = True,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.86, 's^-1'), n=3.27, Ea=(13197.7, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3e-27, 'cm^3/(mol*s)'),
            n = 13.481,
            Ea = (-606.5, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -19.7,
        T3 = (215, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SXC9H19 <=> S2XC9H19""",
)

entry(
    index = 1360,
    label = "SXC8H17 <=> S2XC8H17",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.86, 's^-1'), n=3.27, Ea=(13197.7, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3e-27, 'cm^3/(mol*s)'),
            n = 13.481,
            Ea = (-606.5, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -19.7,
        T3 = (215, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SXC8H17 <=> S2XC8H17""",
)

entry(
    index = 1361,
    label = "S2XC12H25 <=> S4XC12H25",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2.95, 's^-1'), n=3.08, Ea=(12865.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.9e-34, 'cm^3/(mol*s)'),
            n = 15.855,
            Ea = (1243.8, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -15.24,
        T3 = (216, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S2XC12H25 <=> S4XC12H25""",
)

entry(
    index = 1362,
    label = "SXC12H25 <=> S5XC12H25",
    degeneracy = 1,
    duplicate = True,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2.95, 's^-1'), n=3.08, Ea=(12865.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.9e-34, 'cm^3/(mol*s)'),
            n = 15.855,
            Ea = (1243.8, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -15.24,
        T3 = (216, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SXC12H25 <=> S5XC12H25""",
)

entry(
    index = 1363,
    label = "S2XC11H23 <=> S3XC11H23",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2.95, 's^-1'), n=3.08, Ea=(12865.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.9e-34, 'cm^3/(mol*s)'),
            n = 15.855,
            Ea = (1243.8, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -15.24,
        T3 = (216, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S2XC11H23 <=> S3XC11H23""",
)

entry(
    index = 1364,
    label = "SXC11H23 <=> S4XC11H23",
    degeneracy = 1,
    duplicate = True,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2.95, 's^-1'), n=3.08, Ea=(12865.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.9e-34, 'cm^3/(mol*s)'),
            n = 15.855,
            Ea = (1243.8, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -15.24,
        T3 = (216, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SXC11H23 <=> S4XC11H23""",
)

entry(
    index = 1365,
    label = "SXC10H21 <=> S3XC10H21",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2.95, 's^-1'), n=3.08, Ea=(12865.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.9e-34, 'cm^3/(mol*s)'),
            n = 15.855,
            Ea = (1243.8, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -15.24,
        T3 = (216, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SXC10H21 <=> S3XC10H21""",
)

entry(
    index = 1366,
    label = "SXC9H19 <=> S2XC9H19",
    degeneracy = 1,
    duplicate = True,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2.95, 's^-1'), n=3.08, Ea=(12865.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.9e-34, 'cm^3/(mol*s)'),
            n = 15.855,
            Ea = (1243.8, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -15.24,
        T3 = (216, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SXC9H19 <=> S2XC9H19""",
)

entry(
    index = 1367,
    label = "PXC11H23 + CH3 <=> NC12H26",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.93e+14, 'cm^3/(mol*s)'), n=-0.32, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC11H23 + CH3 <=> NC12H26""",
)

entry(
    index = 1368,
    label = "PXC10H21 + C2H5 <=> NC12H26",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.88e+14, 'cm^3/(mol*s)'), n=-0.5, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC10H21 + C2H5 <=> NC12H26""",
)

entry(
    index = 1369,
    label = "PXC9H19 + nC3H7 <=> NC12H26",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.88e+14, 'cm^3/(mol*s)'), n=-0.5, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC9H19 + nC3H7 <=> NC12H26""",
)

entry(
    index = 1370,
    label = "PXC8H17 + pC4H9 <=> NC12H26",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.88e+14, 'cm^3/(mol*s)'), n=-0.5, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC8H17 + pC4H9 <=> NC12H26""",
)

entry(
    index = 1371,
    label = "PXC7H15 + PXC5H11 <=> NC12H26",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.88e+14, 'cm^3/(mol*s)'), n=-0.5, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC7H15 + PXC5H11 <=> NC12H26""",
)

entry(
    index = 1372,
    label = "PXC6H13 + PXC6H13 <=> NC12H26",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.4e+14, 'cm^3/(mol*s)'), n=-0.5, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC6H13 + PXC6H13 <=> NC12H26""",
)

entry(
    index = 1373,
    label = "NC12H26 + OH <=> PXC12H25 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.73e+07, 'cm^3/(mol*s)'),
        n = 1.81,
        Ea = (868.3, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC12H26 + OH <=> PXC12H25 + H2O""",
)

entry(
    index = 1374,
    label = "NC12H26 + OH <=> SXC12H25 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.41e+10, 'cm^3/(mol*s)'),
        n = 0.94,
        Ea = (504.7, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC12H26 + OH <=> SXC12H25 + H2O""",
)

entry(
    index = 1375,
    label = "NC12H26 + OH <=> S2XC12H25 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.14e+07, 'cm^3/(mol*s)'),
        n = 1.81,
        Ea = (-1015.4, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC12H26 + OH <=> S2XC12H25 + H2O""",
)

entry(
    index = 1376,
    label = "NC12H26 + OH <=> S3XC12H25 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.12e+12, 'cm^3/(mol*s)'),
        n = 0.32,
        Ea = (846.5, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC12H26 + OH <=> S3XC12H25 + H2O""",
)

entry(
    index = 1377,
    label = "NC12H26 + OH <=> S4XC12H25 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.12e+12, 'cm^3/(mol*s)'),
        n = 0.32,
        Ea = (846.5, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC12H26 + OH <=> S4XC12H25 + H2O""",
)

entry(
    index = 1378,
    label = "NC12H26 + OH <=> S5XC12H25 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.12e+12, 'cm^3/(mol*s)'),
        n = 0.32,
        Ea = (846.5, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC12H26 + OH <=> S5XC12H25 + H2O""",
)

entry(
    index = 1379,
    label = "NC12H26 + O2 <=> PXC12H25 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(50930, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC12H26 + O2 <=> PXC12H25 + HO2""",
)

entry(
    index = 1380,
    label = "NC12H26 + O2 <=> SXC12H25 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+13, 'cm^3/(mol*s)'), n=0, Ea=(47590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC12H26 + O2 <=> SXC12H25 + HO2""",
)

entry(
    index = 1381,
    label = "NC12H26 + O2 <=> S2XC12H25 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+13, 'cm^3/(mol*s)'), n=0, Ea=(47590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC12H26 + O2 <=> S2XC12H25 + HO2""",
)

entry(
    index = 1382,
    label = "NC12H26 + O2 <=> S3XC12H25 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+13, 'cm^3/(mol*s)'), n=0, Ea=(47590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC12H26 + O2 <=> S3XC12H25 + HO2""",
)

entry(
    index = 1383,
    label = "NC12H26 + O2 <=> S4XC12H25 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+13, 'cm^3/(mol*s)'), n=0, Ea=(47590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC12H26 + O2 <=> S4XC12H25 + HO2""",
)

entry(
    index = 1384,
    label = "NC12H26 + O2 <=> S5XC12H25 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+13, 'cm^3/(mol*s)'), n=0, Ea=(47590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC12H26 + O2 <=> S5XC12H25 + HO2""",
)

entry(
    index = 1385,
    label = "NC12H26 + HO2 <=> PXC12H25 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(61100, 'cm^3/(mol*s)'), n=2.65, Ea=(17496, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC12H26 + HO2 <=> PXC12H25 + H2O2""",
)

entry(
    index = 1386,
    label = "NC12H26 + HO2 <=> SXC12H25 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(14200, 'cm^3/(mol*s)'), n=2.77, Ea=(14913, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC12H26 + HO2 <=> SXC12H25 + H2O2""",
)

entry(
    index = 1387,
    label = "NC12H26 + HO2 <=> S2XC12H25 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(14200, 'cm^3/(mol*s)'), n=2.77, Ea=(14913, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC12H26 + HO2 <=> S2XC12H25 + H2O2""",
)

entry(
    index = 1388,
    label = "NC12H26 + HO2 <=> S3XC12H25 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(14200, 'cm^3/(mol*s)'), n=2.77, Ea=(14913, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC12H26 + HO2 <=> S3XC12H25 + H2O2""",
)

entry(
    index = 1389,
    label = "NC12H26 + HO2 <=> S4XC12H25 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(14200, 'cm^3/(mol*s)'), n=2.77, Ea=(14913, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC12H26 + HO2 <=> S4XC12H25 + H2O2""",
)

entry(
    index = 1390,
    label = "NC12H26 + HO2 <=> S5XC12H25 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(14200, 'cm^3/(mol*s)'), n=2.77, Ea=(14913, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC12H26 + HO2 <=> S5XC12H25 + H2O2""",
)

entry(
    index = 1391,
    label = "NC12H26 + H <=> PXC12H25 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0645, 'cm^3/(mol*s)'), n=4.7, Ea=(3679, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC12H26 + H <=> PXC12H25 + H2""",
)

entry(
    index = 1392,
    label = "NC12H26 + H <=> SXC12H25 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0634, 'cm^3/(mol*s)'), n=4.65, Ea=(1340, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC12H26 + H <=> SXC12H25 + H2""",
)

entry(
    index = 1393,
    label = "NC12H26 + H <=> S2XC12H25 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0634, 'cm^3/(mol*s)'), n=4.65, Ea=(1340, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC12H26 + H <=> S2XC12H25 + H2""",
)

entry(
    index = 1394,
    label = "NC12H26 + H <=> S3XC12H25 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0634, 'cm^3/(mol*s)'), n=4.65, Ea=(1340, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC12H26 + H <=> S3XC12H25 + H2""",
)

entry(
    index = 1395,
    label = "NC12H26 + H <=> S4XC12H25 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0634, 'cm^3/(mol*s)'), n=4.65, Ea=(1340, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC12H26 + H <=> S4XC12H25 + H2""",
)

entry(
    index = 1396,
    label = "NC12H26 + H <=> S5XC12H25 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0634, 'cm^3/(mol*s)'), n=4.65, Ea=(1340, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC12H26 + H <=> S5XC12H25 + H2""",
)

entry(
    index = 1397,
    label = "NC12H26 + O <=> PXC12H25 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.49, 'cm^3/(mol*s)'), n=4.17, Ea=(2766, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC12H26 + O <=> PXC12H25 + OH""",
)

entry(
    index = 1398,
    label = "NC12H26 + O <=> SXC12H25 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(36.5, 'cm^3/(mol*s)'), n=3.75, Ea=(825, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC12H26 + O <=> SXC12H25 + OH""",
)

entry(
    index = 1399,
    label = "NC12H26 + O <=> S2XC12H25 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(36.5, 'cm^3/(mol*s)'), n=3.75, Ea=(825, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC12H26 + O <=> S2XC12H25 + OH""",
)

entry(
    index = 1400,
    label = "NC12H26 + O <=> S3XC12H25 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(36.5, 'cm^3/(mol*s)'), n=3.75, Ea=(825, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC12H26 + O <=> S3XC12H25 + OH""",
)

entry(
    index = 1401,
    label = "NC12H26 + O <=> S4XC12H25 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(36.5, 'cm^3/(mol*s)'), n=3.75, Ea=(825, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC12H26 + O <=> S4XC12H25 + OH""",
)

entry(
    index = 1402,
    label = "NC12H26 + O <=> S5XC12H25 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(36.5, 'cm^3/(mol*s)'), n=3.75, Ea=(825, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC12H26 + O <=> S5XC12H25 + OH""",
)

entry(
    index = 1403,
    label = "NC12H26 + CH3 <=> PXC12H25 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.16e-08, 'cm^3/(mol*s)'),
        n = 6.08,
        Ea = (6223, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC12H26 + CH3 <=> PXC12H25 + CH4""",
)

entry(
    index = 1404,
    label = "NC12H26 + CH3 <=> SXC12H25 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.82e-07, 'cm^3/(mol*s)'),
        n = 5.89,
        Ea = (4768, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC12H26 + CH3 <=> SXC12H25 + CH4""",
)

entry(
    index = 1405,
    label = "NC12H26 + CH3 <=> S2XC12H25 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.82e-07, 'cm^3/(mol*s)'),
        n = 5.89,
        Ea = (4768, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC12H26 + CH3 <=> S2XC12H25 + CH4""",
)

entry(
    index = 1406,
    label = "NC12H26 + CH3 <=> S3XC12H25 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.82e-07, 'cm^3/(mol*s)'),
        n = 5.89,
        Ea = (4768, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC12H26 + CH3 <=> S3XC12H25 + CH4""",
)

entry(
    index = 1407,
    label = "NC12H26 + CH3 <=> S4XC12H25 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.82e-07, 'cm^3/(mol*s)'),
        n = 5.89,
        Ea = (4768, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC12H26 + CH3 <=> S4XC12H25 + CH4""",
)

entry(
    index = 1408,
    label = "NC12H26 + CH3 <=> S5XC12H25 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.82e-07, 'cm^3/(mol*s)'),
        n = 5.89,
        Ea = (4768, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC12H26 + CH3 <=> S5XC12H25 + CH4""",
)

entry(
    index = 1409,
    label = "PXC12H25 + O2 <=> PC12H25O2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3e+48, 'cm^3/(mol*s)'),
        n = -11.66,
        Ea = (10000, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is PXC12H25 + O2 <=> PC12H25O2""",
)

entry(
    index = 1410,
    label = "SXC12H25 + O2 <=> PC12H25O2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3e+48, 'cm^3/(mol*s)'),
        n = -11.66,
        Ea = (10000, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is SXC12H25 + O2 <=> PC12H25O2""",
)

entry(
    index = 1411,
    label = "S2XC12H25 + O2 <=> PC12H25O2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3e+48, 'cm^3/(mol*s)'),
        n = -11.66,
        Ea = (10000, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is S2XC12H25 + O2 <=> PC12H25O2""",
)

entry(
    index = 1412,
    label = "S3XC12H25 + O2 <=> PC12H25O2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3e+48, 'cm^3/(mol*s)'),
        n = -11.66,
        Ea = (10000, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is S3XC12H25 + O2 <=> PC12H25O2""",
)

entry(
    index = 1413,
    label = "S4XC12H25 + O2 <=> PC12H25O2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3e+48, 'cm^3/(mol*s)'),
        n = -11.66,
        Ea = (10000, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is S4XC12H25 + O2 <=> PC12H25O2""",
)

entry(
    index = 1414,
    label = "S5XC12H25 + O2 <=> PC12H25O2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3e+48, 'cm^3/(mol*s)'),
        n = -11.66,
        Ea = (10000, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is S5XC12H25 + O2 <=> PC12H25O2""",
)

entry(
    index = 1415,
    label = "PC12H25O2 => P12OOHX2",
    degeneracy = 1,
    duplicate = True,
    reversible = False,
    kinetics = Arrhenius(A=(2e+12, 's^-1'), n=0, Ea=(17017.2, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PC12H25O2 => P12OOHX2""",
)

entry(
    index = 1416,
    label = "P12OOHX2 => PC12H25O2",
    degeneracy = 1,
    duplicate = True,
    reversible = False,
    kinetics = Arrhenius(A=(1e+11, 's^-1'), n=0, Ea=(12500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is P12OOHX2 => PC12H25O2""",
)

entry(
    index = 1417,
    label = "P12OOHX2 <=> C12H24 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.5e+12, 's^-1'), n=0, Ea=(25573.6, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is P12OOHX2 <=> C12H24 + HO2""",
)

entry(
    index = 1418,
    label = "P12OOHX2 + O2 => OC12OOH + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(4e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is P12OOHX2 + O2 => OC12OOH + OH""",
)

entry(
    index = 1419,
    label = "OC12OOH => CH2O + C2H4 + C2H4 + C2H4 + C2H4 + C2H5 + OH + CO",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(7e+14, 's^-1'), n=0, Ea=(42065, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is OC12OOH => CH2O + C2H4 + C2H4 + C2H4 + C2H4 + C2H5 + OH + CO""",
)

entry(
    index = 1420,
    label = "PXC3H6cC6H11 + CH3 <=> C4H9cC6H11",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.93e+14, 'cm^3/(mol*s)'), n=-0.32, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC3H6cC6H11 + CH3 <=> C4H9cC6H11""",
)

entry(
    index = 1421,
    label = "PXC2H4cC6H11 + C2H5 <=> C4H9cC6H11",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.88e+14, 'cm^3/(mol*s)'), n=-0.5, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC2H4cC6H11 + C2H5 <=> C4H9cC6H11""",
)

entry(
    index = 1422,
    label = "PXCH2cC6H11 + nC3H7 <=> C4H9cC6H11",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.88e+14, 'cm^3/(mol*s)'), n=-0.5, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXCH2cC6H11 + nC3H7 <=> C4H9cC6H11""",
)

entry(
    index = 1423,
    label = "cC6H11 + pC4H9 <=> C4H9cC6H11",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.29e+14, 'cm^3/(mol*s)'), n=-0.35, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC6H11 + pC4H9 <=> C4H9cC6H11""",
)

entry(
    index = 1424,
    label = "C4H9cC6H11 <=> C10H20-5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.17e+15, 's^-1'), n=0, Ea=(74000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H11 <=> C10H20-5""",
)

entry(
    index = 1425,
    label = "C4H9cC6H11 <=> C4H9-2-1C6H11",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.17e+15, 's^-1'), n=0, Ea=(74000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H11 <=> C4H9-2-1C6H11""",
)

entry(
    index = 1426,
    label = "C10H20-5 <=> SAXC7H13 + nC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.14e+23, 's^-1'), n=-2.03, Ea=(74958, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C10H20-5 <=> SAXC7H13 + nC3H7""",
)

entry(
    index = 1427,
    label = "C10H20-5 <=> C7H14 + C3H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.42e+07, 's^-1'), n=1.65, Ea=(53752, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C10H20-5 <=> C7H14 + C3H6""",
)

entry(
    index = 1428,
    label = "C10H20-5 + H <=> nC3H7 + C7H14",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8e+21, 'cm^3/(mol*s)'),
        n = -2.39,
        Ea = (11180, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C10H20-5 + H <=> nC3H7 + C7H14""",
)

entry(
    index = 1429,
    label = "C10H20-5 + H <=> pC4H9 + C6H12",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8e+21, 'cm^3/(mol*s)'),
        n = -2.39,
        Ea = (11180, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C10H20-5 + H <=> pC4H9 + C6H12""",
)

entry(
    index = 1430,
    label = "C4H9-2-1C6H11 <=> PAXCH2-2-1C6H11 + nC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.14e+23, 's^-1'), n=-2.03, Ea=(74958, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9-2-1C6H11 <=> PAXCH2-2-1C6H11 + nC3H7""",
)

entry(
    index = 1431,
    label = "C4H9-2-1C6H11 <=> CH3-2-1C6H11 + C3H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.42e+07, 's^-1'), n=1.65, Ea=(53752, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9-2-1C6H11 <=> CH3-2-1C6H11 + C3H6""",
)

entry(
    index = 1432,
    label = "C4H9-2-1C6H11 + H <=> nC3H7 + CH3-2-1C6H11",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.6e+22, 'cm^3/(mol*s)'),
        n = -2.39,
        Ea = (11180, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C4H9-2-1C6H11 + H <=> nC3H7 + CH3-2-1C6H11""",
)

entry(
    index = 1433,
    label = "C4H9-2-1C6H11 + H <=> C6H12 + pC4H9",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8e+21, 'cm^3/(mol*s)'),
        n = -2.39,
        Ea = (11180, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C4H9-2-1C6H11 + H <=> C6H12 + pC4H9""",
)

entry(
    index = 1434,
    label = "C4H9cC6H11 + H <=> S3XC4H8cC6H11 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(310000, 'cm^3/(mol*s)'), n=2.51, Ea=(4286, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H11 + H <=> S3XC4H8cC6H11 + H2""",
)

entry(
    index = 1435,
    label = "C4H9cC6H11 + H <=> S2XC4H8cC6H11 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(235000, 'cm^3/(mol*s)'), n=2.53, Ea=(4218, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H11 + H <=> S2XC4H8cC6H11 + H2""",
)

entry(
    index = 1436,
    label = "C4H9cC6H11 + H <=> SXC4H8cC6H11 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(380000, 'cm^3/(mol*s)'), n=2.5, Ea=(4467, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H11 + H <=> SXC4H8cC6H11 + H2""",
)

entry(
    index = 1437,
    label = "C4H9cC6H11 + H <=> PXC4H8cC6H11 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(306000, 'cm^3/(mol*s)'), n=2.59, Ea=(6603, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H11 + H <=> PXC4H8cC6H11 + H2""",
)

entry(
    index = 1438,
    label = "C4H9cC6H11 + H <=> C4H9TXcC6H10 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(383000, 'cm^3/(mol*s)'), n=2.41, Ea=(2653, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H11 + H <=> C4H9TXcC6H10 + H2""",
)

entry(
    index = 1439,
    label = "C4H9cC6H11 + H <=> C4H9S2XcC6H10 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(456000, 'cm^3/(mol*s)'), n=2.54, Ea=(5324, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H11 + H <=> C4H9S2XcC6H10 + H2""",
)

entry(
    index = 1440,
    label = "C4H9cC6H11 + H <=> C4H9S3XcC6H10 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(276000, 'cm^3/(mol*s)'), n=2.6, Ea=(4078, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H11 + H <=> C4H9S3XcC6H10 + H2""",
)

entry(
    index = 1441,
    label = "C4H9cC6H11 + H <=> C4H9S4XcC6H10 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(285000, 'cm^3/(mol*s)'), n=2.52, Ea=(4463, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H11 + H <=> C4H9S4XcC6H10 + H2""",
)

entry(
    index = 1442,
    label = "C4H9cC6H11 + O <=> PXC4H8cC6H11 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.2e+06, 'cm^3/(mol*s)'), n=2.4, Ea=(5504, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H11 + O <=> PXC4H8cC6H11 + OH""",
)

entry(
    index = 1443,
    label = "C4H9cC6H11 + O <=> SXC4H8cC6H11 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(310000, 'cm^3/(mol*s)'), n=2.5, Ea=(2225, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H11 + O <=> SXC4H8cC6H11 + OH""",
)

entry(
    index = 1444,
    label = "C4H9cC6H11 + O <=> S2XC4H8cC6H11 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(310000, 'cm^3/(mol*s)'), n=2.5, Ea=(2225, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H11 + O <=> S2XC4H8cC6H11 + OH""",
)

entry(
    index = 1445,
    label = "C4H9cC6H11 + O <=> S3XC4H8cC6H11 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(310000, 'cm^3/(mol*s)'), n=2.5, Ea=(2225, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H11 + O <=> S3XC4H8cC6H11 + OH""",
)

entry(
    index = 1446,
    label = "C4H9cC6H11 + O <=> C4H9TXcC6H10 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(157000, 'cm^3/(mol*s)'), n=2.5, Ea=(1110, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H11 + O <=> C4H9TXcC6H10 + OH""",
)

entry(
    index = 1447,
    label = "C4H9cC6H11 + O <=> C4H9S2XcC6H10 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(95200, 'cm^3/(mol*s)'), n=2.71, Ea=(2106, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H11 + O <=> C4H9S2XcC6H10 + OH""",
)

entry(
    index = 1448,
    label = "C4H9cC6H11 + O <=> C4H9S3XcC6H10 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(95200, 'cm^3/(mol*s)'), n=2.71, Ea=(2106, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H11 + O <=> C4H9S3XcC6H10 + OH""",
)

entry(
    index = 1449,
    label = "C4H9cC6H11 + O <=> C4H9S4XcC6H10 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(47600, 'cm^3/(mol*s)'), n=2.71, Ea=(2106, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H11 + O <=> C4H9S4XcC6H10 + OH""",
)

entry(
    index = 1450,
    label = "C4H9cC6H11 + OH <=> S3XC4H8cC6H11 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4640, 'cm^3/(mol*s)'), n=2.83, Ea=(-1746, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H11 + OH <=> S3XC4H8cC6H11 + H2O""",
)

entry(
    index = 1451,
    label = "C4H9cC6H11 + OH <=> S2XC4H8cC6H11 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3490, 'cm^3/(mol*s)'), n=2.82, Ea=(-1428, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H11 + OH <=> S2XC4H8cC6H11 + H2O""",
)

entry(
    index = 1452,
    label = "C4H9cC6H11 + OH <=> SXC4H8cC6H11 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3670, 'cm^3/(mol*s)'), n=2.87, Ea=(-1026, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H11 + OH <=> SXC4H8cC6H11 + H2O""",
)

entry(
    index = 1453,
    label = "C4H9cC6H11 + OH <=> PXC4H8cC6H11 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7530, 'cm^3/(mol*s)'), n=2.9, Ea=(605, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H11 + OH <=> PXC4H8cC6H11 + H2O""",
)

entry(
    index = 1454,
    label = "C4H9cC6H11 + OH <=> C4H9TXcC6H10 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1430, 'cm^3/(mol*s)'), n=2.92, Ea=(-2653, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H11 + OH <=> C4H9TXcC6H10 + H2O""",
)

entry(
    index = 1455,
    label = "C4H9cC6H11 + OH <=> C4H9S2XcC6H10 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(13000, 'cm^3/(mol*s)'), n=2.86, Ea=(-1846, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H11 + OH <=> C4H9S2XcC6H10 + H2O""",
)

entry(
    index = 1456,
    label = "C4H9cC6H11 + OH <=> C4H9S3XcC6H10 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(10300, 'cm^3/(mol*s)'), n=2.86, Ea=(-1197, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H11 + OH <=> C4H9S3XcC6H10 + H2O""",
)

entry(
    index = 1457,
    label = "C4H9cC6H11 + OH <=> C4H9S4XcC6H10 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8040, 'cm^3/(mol*s)'), n=2.88, Ea=(-1226, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H11 + OH <=> C4H9S4XcC6H10 + H2O""",
)

entry(
    index = 1458,
    label = "C4H9cC6H11 + O2 <=> PXC4H8cC6H11 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(50930, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H11 + O2 <=> PXC4H8cC6H11 + HO2""",
)

entry(
    index = 1459,
    label = "C4H9cC6H11 + O2 <=> SXC4H8cC6H11 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(47590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H11 + O2 <=> SXC4H8cC6H11 + HO2""",
)

entry(
    index = 1460,
    label = "C4H9cC6H11 + O2 <=> S2XC4H8cC6H11 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(47590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H11 + O2 <=> S2XC4H8cC6H11 + HO2""",
)

entry(
    index = 1461,
    label = "C4H9cC6H11 + O2 <=> S3XC4H8cC6H11 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(47590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H11 + O2 <=> S3XC4H8cC6H11 + HO2""",
)

entry(
    index = 1462,
    label = "C4H9cC6H11 + O2 <=> C4H9TXcC6H10 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(44000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H11 + O2 <=> C4H9TXcC6H10 + HO2""",
)

entry(
    index = 1463,
    label = "C4H9cC6H11 + O2 <=> C4H9S2XcC6H10 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+13, 'cm^3/(mol*s)'), n=0, Ea=(47590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H11 + O2 <=> C4H9S2XcC6H10 + HO2""",
)

entry(
    index = 1464,
    label = "C4H9cC6H11 + O2 <=> C4H9S3XcC6H10 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+13, 'cm^3/(mol*s)'), n=0, Ea=(47590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H11 + O2 <=> C4H9S3XcC6H10 + HO2""",
)

entry(
    index = 1465,
    label = "C4H9cC6H11 + O2 <=> C4H9S4XcC6H10 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(47590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H11 + O2 <=> C4H9S4XcC6H10 + HO2""",
)

entry(
    index = 1466,
    label = "C4H9cC6H11 + HO2 <=> PXC4H8cC6H11 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(30500, 'cm^3/(mol*s)'), n=2.65, Ea=(17496, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H11 + HO2 <=> PXC4H8cC6H11 + H2O2""",
)

entry(
    index = 1467,
    label = "C4H9cC6H11 + HO2 <=> SXC4H8cC6H11 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7130, 'cm^3/(mol*s)'), n=2.77, Ea=(14913, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H11 + HO2 <=> SXC4H8cC6H11 + H2O2""",
)

entry(
    index = 1468,
    label = "C4H9cC6H11 + HO2 <=> S2XC4H8cC6H11 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7130, 'cm^3/(mol*s)'), n=2.77, Ea=(14913, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H11 + HO2 <=> S2XC4H8cC6H11 + H2O2""",
)

entry(
    index = 1469,
    label = "C4H9cC6H11 + HO2 <=> S3XC4H8cC6H11 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7130, 'cm^3/(mol*s)'), n=2.77, Ea=(14913, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H11 + HO2 <=> S3XC4H8cC6H11 + H2O2""",
)

entry(
    index = 1470,
    label = "C4H9cC6H11 + HO2 <=> C4H9TXcC6H10 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1240, 'cm^3/(mol*s)'), n=2.77, Ea=(10500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H11 + HO2 <=> C4H9TXcC6H10 + H2O2""",
)

entry(
    index = 1471,
    label = "C4H9cC6H11 + HO2 <=> C4H9S2XcC6H10 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(281, 'cm^3/(mol*s)'), n=3.25, Ea=(14998, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H11 + HO2 <=> C4H9S2XcC6H10 + H2O2""",
)

entry(
    index = 1472,
    label = "C4H9cC6H11 + HO2 <=> C4H9S3XcC6H10 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(281, 'cm^3/(mol*s)'), n=3.25, Ea=(14998, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H11 + HO2 <=> C4H9S3XcC6H10 + H2O2""",
)

entry(
    index = 1473,
    label = "C4H9cC6H11 + HO2 <=> C4H9S4XcC6H10 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(140, 'cm^3/(mol*s)'), n=3.25, Ea=(14998, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H11 + HO2 <=> C4H9S4XcC6H10 + H2O2""",
)

entry(
    index = 1474,
    label = "C4H9cC6H11 + CH3 <=> S3XC4H8cC6H11 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(99.1, 'cm^3/(mol*s)'), n=3.26, Ea=(13774, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H11 + CH3 <=> S3XC4H8cC6H11 + CH4""",
)

entry(
    index = 1475,
    label = "C4H9cC6H11 + CH3 <=> S2XC4H8cC6H11 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(32.3, 'cm^3/(mol*s)'), n=3.29, Ea=(11053, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H11 + CH3 <=> S2XC4H8cC6H11 + CH4""",
)

entry(
    index = 1476,
    label = "C4H9cC6H11 + CH3 <=> SXC4H8cC6H11 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(178, 'cm^3/(mol*s)'), n=3.14, Ea=(11394, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H11 + CH3 <=> SXC4H8cC6H11 + CH4""",
)

entry(
    index = 1477,
    label = "C4H9cC6H11 + CH3 <=> PXC4H8cC6H11 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(27.5, 'cm^3/(mol*s)'), n=3.26, Ea=(10901, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H11 + CH3 <=> PXC4H8cC6H11 + CH4""",
)

entry(
    index = 1478,
    label = "C4H9cC6H11 + CH3 <=> C4H9TXcC6H10 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.99, 'cm^3/(mol*s)'), n=3.29, Ea=(8715, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H11 + CH3 <=> C4H9TXcC6H10 + CH4""",
)

entry(
    index = 1479,
    label = "C4H9cC6H11 + CH3 <=> C4H9S2XcC6H10 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(228, 'cm^3/(mol*s)'), n=3.2, Ea=(11161, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H11 + CH3 <=> C4H9S2XcC6H10 + CH4""",
)

entry(
    index = 1480,
    label = "C4H9cC6H11 + CH3 <=> C4H9S3XcC6H10 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(215, 'cm^3/(mol*s)'), n=3.17, Ea=(11468, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H11 + CH3 <=> C4H9S3XcC6H10 + CH4""",
)

entry(
    index = 1481,
    label = "C4H9cC6H11 + CH3 <=> C4H9S4XcC6H10 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(143, 'cm^3/(mol*s)'), n=3.24, Ea=(11347, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H11 + CH3 <=> C4H9S4XcC6H10 + CH4""",
)

entry(
    index = 1482,
    label = "PXC4H8cC6H11 <=> PXC2H4cC6H11 + C2H4",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(9.12e+11, 's^-1'), n=0.31, Ea=(27237.8, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.8e-57, 'cm^3/(mol*s)'),
            n = 23.463,
            Ea = (-602.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -2.46,
        T3 = (206, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC4H8cC6H11 <=> PXC2H4cC6H11 + C2H4""",
)

entry(
    index = 1483,
    label = "SXC4H8cC6H11 <=> PXCH2cC6H11 + C3H6",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(6.03e+10, 's^-1'), n=0.84, Ea=(27820, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1e-43, 'cm^3/(mol*s)'),
            n = 18.591,
            Ea = (-602.5, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -43.32,
        T3 = (200, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SXC4H8cC6H11 <=> PXCH2cC6H11 + C3H6""",
)

entry(
    index = 1484,
    label = "S2XC4H8cC6H11 <=> cC6H11 + C4H81",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(5.81e+10, 's^-1'), n=0.85, Ea=(23824.1, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (8.9e-37, 'cm^3/(mol*s)'),
            n = 16.685,
            Ea = (-661.7, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -64.15,
        T3 = (179, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S2XC4H8cC6H11 <=> cC6H11 + C4H81""",
)

entry(
    index = 1485,
    label = "S2XC4H8cC6H11 <=> C3H5cC6H11 + CH3",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(9.55e+09, 's^-1'), n=1.08, Ea=(29387.7, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.3e-46, 'cm^3/(mol*s)'),
            n = 19.133,
            Ea = (-602.7, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -34.36,
        T3 = (210, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S2XC4H8cC6H11 <=> C3H5cC6H11 + CH3""",
)

entry(
    index = 1486,
    label = "S3XC4H8cC6H11 <=> C2H3cC6H11 + C2H5",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(6.76e+09, 's^-1'), n=1.11, Ea=(27023.2, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (8.2e-43, 'cm^3/(mol*s)'),
            n = 18.276,
            Ea = (-602.9, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -30.04,
        T3 = (210, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S3XC4H8cC6H11 <=> C2H3cC6H11 + C2H5""",
)

entry(
    index = 1487,
    label = "S3XC4H8cC6H11 <=> PX10-4C10H19",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(6.03e+12, 's^-1'), n=0.07, Ea=(26982.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.1e-33, 'cm^3/(mol*s)'),
            n = 15.29,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -25.11,
        T3 = (225, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S3XC4H8cC6H11 <=> PX10-4C10H19""",
)

entry(
    index = 1488,
    label = "C4H9TXcC6H10 <=> CH2cC6H10 + nC3H7",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.63e+11, 's^-1'), n=0.68, Ea=(24579.2, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.1e-39, 'cm^3/(mol*s)'),
            n = 17.57,
            Ea = (-603.1, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -31.51,
        T3 = (219, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C4H9TXcC6H10 <=> CH2cC6H10 + nC3H7""",
)

entry(
    index = 1489,
    label = "C4H9TXcC6H10 <=> PXC4H8-2-1C6H11",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(6.03e+12, 's^-1'), n=0.07, Ea=(27982.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.1e-33, 'cm^3/(mol*s)'),
            n = 15.29,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -25.11,
        T3 = (225, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C4H9TXcC6H10 <=> PXC4H8-2-1C6H11""",
)

entry(
    index = 1490,
    label = "C4H9S2XcC6H10 <=> cC6H10 + pC4H9",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(5.47e+11, 's^-1'), n=0.5, Ea=(27798.1, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (8.5e-37, 'cm^3/(mol*s)'),
            n = 16.21,
            Ea = (-603.1, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -26.7,
        T3 = (216, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C4H9S2XcC6H10 <=> cC6H10 + pC4H9""",
)

entry(
    index = 1491,
    label = "C4H9S2XcC6H10 <=> PX10-5C10H19",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.01e+12, 's^-1'), n=0.07, Ea=(26982.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.1e-33, 'cm^3/(mol*s)'),
            n = 15.29,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -25.11,
        T3 = (225, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C4H9S2XcC6H10 <=> PX10-5C10H19""",
)

entry(
    index = 1492,
    label = "C4H9S2XcC6H10 <=> PXC3H6-3-1C7H13",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.01e+12, 's^-1'), n=0.07, Ea=(27982.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.1e-33, 'cm^3/(mol*s)'),
            n = 15.29,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -25.11,
        T3 = (225, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C4H9S2XcC6H10 <=> PXC3H6-3-1C7H13""",
)

entry(
    index = 1493,
    label = "C4H9S3XcC6H10 <=> S4XC10H19",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.01e+12, 's^-1'), n=0.07, Ea=(26982.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.1e-33, 'cm^3/(mol*s)'),
            n = 15.29,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -25.11,
        T3 = (225, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C4H9S3XcC6H10 <=> S4XC10H19""",
)

entry(
    index = 1494,
    label = "C4H9S3XcC6H10 <=> PXC2H4-4-1C8H15",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.01e+12, 's^-1'), n=0.07, Ea=(27982.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.1e-33, 'cm^3/(mol*s)'),
            n = 15.29,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -25.11,
        T3 = (225, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C4H9S3XcC6H10 <=> PXC2H4-4-1C8H15""",
)

entry(
    index = 1495,
    label = "C4H9S4XcC6H10 <=> PXCH2-5-1C9H17",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(6.03e+12, 's^-1'), n=0.07, Ea=(27982.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.1e-33, 'cm^3/(mol*s)'),
            n = 15.29,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -25.11,
        T3 = (225, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C4H9S4XcC6H10 <=> PXCH2-5-1C9H17""",
)

entry(
    index = 1496,
    label = "PXC4H8cC6H11 + H <=> C4H9cC6H11",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.5e+48, 'cm^6/(mol^2*s)'),
            n = -9.32,
            Ea = (5833.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.498,
        T3 = (1314, 'K'),
        T1 = (1314, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC4H8cC6H11 + H <=> C4H9cC6H11""",
)

entry(
    index = 1497,
    label = "PXC4H8cC6H11 + H <=> cC6H11 + pC4H9",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.85e+24, 'cm^3/(mol*s)'),
        n = -2.92,
        Ea = (12505, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is PXC4H8cC6H11 + H <=> cC6H11 + pC4H9""",
)

entry(
    index = 1498,
    label = "PXC4H8cC6H11 + O <=> PXC3H6cC6H11 + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC4H8cC6H11 + O <=> PXC3H6cC6H11 + CH2O""",
)

entry(
    index = 1499,
    label = "SXC4H8cC6H11 + H <=> C4H9cC6H11",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.7e+58, 'cm^6/(mol^2*s)'),
            n = -12.08,
            Ea = (11263.7, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.649,
        T3 = (1213.1, 'K'),
        T1 = (1213.1, 'K'),
        T2 = (13369.7, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SXC4H8cC6H11 + H <=> C4H9cC6H11""",
)

entry(
    index = 1500,
    label = "SXC4H8cC6H11 + H <=> cC6H11 + pC4H9",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.4e+28, 'cm^3/(mol*s)'),
        n = -3.94,
        Ea = (15916, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is SXC4H8cC6H11 + H <=> cC6H11 + pC4H9""",
)

entry(
    index = 1501,
    label = "S2XC4H8cC6H11 + H <=> C4H9cC6H11",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.7e+58, 'cm^6/(mol^2*s)'),
            n = -12.08,
            Ea = (11263.7, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.649,
        T3 = (1213.1, 'K'),
        T1 = (1213.1, 'K'),
        T2 = (13369.7, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S2XC4H8cC6H11 + H <=> C4H9cC6H11""",
)

entry(
    index = 1502,
    label = "S2XC4H8cC6H11 + H <=> cC6H11 + pC4H9",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.4e+28, 'cm^3/(mol*s)'),
        n = -3.94,
        Ea = (15916, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is S2XC4H8cC6H11 + H <=> cC6H11 + pC4H9""",
)

entry(
    index = 1503,
    label = "S3XC4H8cC6H11 + H <=> C4H9cC6H11",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.7e+58, 'cm^6/(mol^2*s)'),
            n = -12.08,
            Ea = (11263.7, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.649,
        T3 = (1213.1, 'K'),
        T1 = (1213.1, 'K'),
        T2 = (13369.7, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S3XC4H8cC6H11 + H <=> C4H9cC6H11""",
)

entry(
    index = 1504,
    label = "S3XC4H8cC6H11 + H <=> cC6H11 + pC4H9",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.4e+28, 'cm^3/(mol*s)'),
        n = -3.94,
        Ea = (15916, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is S3XC4H8cC6H11 + H <=> cC6H11 + pC4H9""",
)

entry(
    index = 1505,
    label = "C4H9TXcC6H10 + H <=> C4H9cC6H11",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.47e+61, 'cm^6/(mol^2*s)'),
            n = -12.94,
            Ea = (8000, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        T3 = (1456.4, 'K'),
        T1 = (1000, 'K'),
        T2 = (10000.5, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C4H9TXcC6H10 + H <=> C4H9cC6H11""",
)

entry(
    index = 1506,
    label = "C4H9TXcC6H10 + H <=> cC6H11 + pC4H9",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.6e+36, 'cm^3/(mol*s)'),
        n = -6.12,
        Ea = (25640, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C4H9TXcC6H10 + H <=> cC6H11 + pC4H9""",
)

entry(
    index = 1507,
    label = "C4H9S2XcC6H10 + H <=> C4H9cC6H11",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(4.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.7e+58, 'cm^6/(mol^2*s)'),
            n = -12.08,
            Ea = (11263.7, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.649,
        T3 = (1213.1, 'K'),
        T1 = (1213.1, 'K'),
        T2 = (13369.7, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C4H9S2XcC6H10 + H <=> C4H9cC6H11""",
)

entry(
    index = 1508,
    label = "C4H9S2XcC6H10 + H <=> cC6H11 + pC4H9",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.8e+28, 'cm^3/(mol*s)'),
        n = -3.94,
        Ea = (15916, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C4H9S2XcC6H10 + H <=> cC6H11 + pC4H9""",
)

entry(
    index = 1509,
    label = "C4H9S3XcC6H10 + H <=> C4H9cC6H11",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(4.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.7e+58, 'cm^6/(mol^2*s)'),
            n = -12.08,
            Ea = (11263.7, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.649,
        T3 = (1213.1, 'K'),
        T1 = (1213.1, 'K'),
        T2 = (13369.7, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C4H9S3XcC6H10 + H <=> C4H9cC6H11""",
)

entry(
    index = 1510,
    label = "C4H9S3XcC6H10 + H <=> cC6H11 + pC4H9",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.8e+28, 'cm^3/(mol*s)'),
        n = -3.94,
        Ea = (15916, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C4H9S3XcC6H10 + H <=> cC6H11 + pC4H9""",
)

entry(
    index = 1511,
    label = "C4H9S4XcC6H10 + H <=> C4H9cC6H11",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.7e+58, 'cm^6/(mol^2*s)'),
            n = -12.08,
            Ea = (11263.7, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.649,
        T3 = (1213.1, 'K'),
        T1 = (1213.1, 'K'),
        T2 = (13369.7, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C4H9S4XcC6H10 + H <=> C4H9cC6H11""",
)

entry(
    index = 1512,
    label = "C4H9S4XcC6H10 + H <=> cC6H11 + pC4H9",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.4e+28, 'cm^3/(mol*s)'),
        n = -3.94,
        Ea = (15916, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C4H9S4XcC6H10 + H <=> cC6H11 + pC4H9""",
)

entry(
    index = 1513,
    label = "PXC4H8cC6H11 <=> S3XC4H8cC6H11",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.14e+08, 's^-1'), n=0.93, Ea=(17627, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC4H8cC6H11 <=> S3XC4H8cC6H11""",
)

entry(
    index = 1514,
    label = "PXC4H8cC6H11 <=> C4H9TXcC6H10",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.1e+08, 's^-1'), n=0.73, Ea=(9950, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC4H8cC6H11 <=> C4H9TXcC6H10""",
)

entry(
    index = 1515,
    label = "PXC4H8cC6H11 <=> C4H9S2XcC6H10",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+09, 's^-1'), n=0.74, Ea=(11027, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC4H8cC6H11 <=> C4H9S2XcC6H10""",
)

entry(
    index = 1516,
    label = "PXC4H8cC6H11 <=> C4H9S3XcC6H10",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.76e+07, 's^-1'), n=0.8, Ea=(18020, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC4H8cC6H11 <=> C4H9S3XcC6H10""",
)

entry(
    index = 1517,
    label = "PXC4H8cC6H11 <=> C4H9S4XcC6H10",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.33e+07, 's^-1'), n=0.84, Ea=(26554, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC4H8cC6H11 <=> C4H9S4XcC6H10""",
)

entry(
    index = 1518,
    label = "SXC4H8cC6H11 <=> C4H9TXcC6H10",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.57e+08, 's^-1'), n=0.93, Ea=(17967, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC4H8cC6H11 <=> C4H9TXcC6H10""",
)

entry(
    index = 1519,
    label = "SXC4H8cC6H11 <=> C4H9S2XcC6H10",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.05e+09, 's^-1'), n=0.82, Ea=(15027, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC4H8cC6H11 <=> C4H9S2XcC6H10""",
)

entry(
    index = 1520,
    label = "SXC4H8cC6H11 <=> C4H9S3XcC6H10",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.95e+08, 's^-1'), n=0.79, Ea=(18372, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC4H8cC6H11 <=> C4H9S3XcC6H10""",
)

entry(
    index = 1521,
    label = "SXC4H8cC6H11 <=> C4H9S4XcC6H10",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.75e+08, 's^-1'), n=0.78, Ea=(24300, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC4H8cC6H11 <=> C4H9S4XcC6H10""",
)

entry(
    index = 1522,
    label = "S2XC4H8cC6H11 <=> C4H9S2XcC6H10",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.58e+08, 's^-1'), n=0.9, Ea=(22900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is S2XC4H8cC6H11 <=> C4H9S2XcC6H10""",
)

entry(
    index = 1523,
    label = "S2XC4H8cC6H11 <=> C4H9S3XcC6H10",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.94e+08, 's^-1'), n=0.78, Ea=(18644, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is S2XC4H8cC6H11 <=> C4H9S3XcC6H10""",
)

entry(
    index = 1524,
    label = "S2XC4H8cC6H11 <=> C4H9S4XcC6H10",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.97e+08, 's^-1'), n=0.78, Ea=(22740, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is S2XC4H8cC6H11 <=> C4H9S4XcC6H10""",
)

entry(
    index = 1525,
    label = "S3XC4H8cC6H11 <=> C4H9S3XcC6H10",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.65e+08, 's^-1'), n=1.02, Ea=(28687, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is S3XC4H8cC6H11 <=> C4H9S3XcC6H10""",
)

entry(
    index = 1526,
    label = "S3XC4H8cC6H11 <=> C4H9S4XcC6H10",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.19e+09, 's^-1'), n=0.92, Ea=(22700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is S3XC4H8cC6H11 <=> C4H9S4XcC6H10""",
)

entry(
    index = 1527,
    label = "PX10-4C10H19 <=> SAX6-4C10H19",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.46e+11, 's^-1'), n=0, Ea=(10516.8, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (9.9e-38, 'cm^3/(mol*s)'),
            n = 17.215,
            Ea = (-603, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -16.33,
        T3 = (200, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PX10-4C10H19 <=> SAX6-4C10H19""",
)

entry(
    index = 1528,
    label = "PX10-4C10H19 <=> PX1-4C8H15 + C2H4",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(9.12e+11, 's^-1'), n=0.31, Ea=(27237.8, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.8e-57, 'cm^3/(mol*s)'),
            n = 23.463,
            Ea = (-602.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -2.46,
        T3 = (206, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PX10-4C10H19 <=> PX1-4C8H15 + C2H4""",
)

entry(
    index = 1529,
    label = "PXC4H8-2-1C6H11 <=> SAXC4H8-2-1C6H11",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(155, 's^-1'), n=2.83, Ea=(15566.2, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.5e-30, 'cm^3/(mol*s)'),
            n = 14.56,
            Ea = (-602.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -13.59,
        T3 = (214, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC4H8-2-1C6H11 <=> SAXC4H8-2-1C6H11""",
)

entry(
    index = 1530,
    label = "PXC4H8-2-1C6H11 <=> PXC2H4-2-1C6H11 + C2H4",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(9.12e+11, 's^-1'), n=0.31, Ea=(27237.8, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.8e-57, 'cm^3/(mol*s)'),
            n = 23.463,
            Ea = (-602.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -2.46,
        T3 = (206, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC4H8-2-1C6H11 <=> PXC2H4-2-1C6H11 + C2H4""",
)

entry(
    index = 1531,
    label = "PX10-5C10H19 <=> SAX4-5C10H19",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(155, 's^-1'), n=2.83, Ea=(15566.2, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.5e-30, 'cm^3/(mol*s)'),
            n = 14.56,
            Ea = (-602.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -13.59,
        T3 = (214, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PX10-5C10H19 <=> SAX4-5C10H19""",
)

entry(
    index = 1532,
    label = "PX10-5C10H19 <=> PX1-3C8H15 + C2H4",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(9.12e+11, 's^-1'), n=0.31, Ea=(27237.8, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.8e-57, 'cm^3/(mol*s)'),
            n = 23.463,
            Ea = (-602.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -2.46,
        T3 = (206, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PX10-5C10H19 <=> PX1-3C8H15 + C2H4""",
)

entry(
    index = 1533,
    label = "PXC3H6-3-1C7H13 <=> C3H7-3-TAX1C7H13",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(86, 's^-1'), n=2.62, Ea=(8722.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (8.1e-33, 'cm^3/(mol*s)'),
            n = 15.214,
            Ea = (-677.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -30.39,
        T3 = (206, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC3H6-3-1C7H13 <=> C3H7-3-TAX1C7H13""",
)

entry(
    index = 1534,
    label = "PXC3H6-3-1C7H13 <=> PXCH2-3-1C7H13 + C2H4",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(9.12e+11, 's^-1'), n=0.31, Ea=(27237.8, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.8e-57, 'cm^3/(mol*s)'),
            n = 23.463,
            Ea = (-602.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -2.46,
        T3 = (206, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC3H6-3-1C7H13 <=> PXCH2-3-1C7H13 + C2H4""",
)

entry(
    index = 1535,
    label = "S4XC10H19 <=> SAXC10H19",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(155, 's^-1'), n=2.83, Ea=(15566.2, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.5e-30, 'cm^3/(mol*s)'),
            n = 14.56,
            Ea = (-602.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -13.59,
        T3 = (214, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S4XC10H19 <=> SAXC10H19""",
)

entry(
    index = 1536,
    label = "S4XC10H19 <=> C6H12 + C4H7",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(5.5e+11, 's^-1'), n=0.55, Ea=(28084.3, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.1e-43, 'cm^3/(mol*s)'),
            n = 18.418,
            Ea = (-602.9, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -32.13,
        T3 = (207, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S4XC10H19 <=> C6H12 + C4H7""",
)

entry(
    index = 1537,
    label = "S4XC10H19 <=> C7H12-16 + nC3H7",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(5.5e+11, 's^-1'), n=0.55, Ea=(28084.3, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.1e-43, 'cm^3/(mol*s)'),
            n = 18.418,
            Ea = (-602.9, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -32.13,
        T3 = (207, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S4XC10H19 <=> C7H12-16 + nC3H7""",
)

entry(
    index = 1538,
    label = "PXC2H4-4-1C8H15 <=> C2H5-4-SAX1C8H14",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(155, 's^-1'), n=2.83, Ea=(15566.2, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.5e-30, 'cm^3/(mol*s)'),
            n = 14.56,
            Ea = (-602.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -13.59,
        T3 = (214, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC2H4-4-1C8H15 <=> C2H5-4-SAX1C8H14""",
)

entry(
    index = 1539,
    label = "PXC2H4-4-1C8H15 <=> S4XC8H15 + C2H4",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(9.12e+11, 's^-1'), n=0.31, Ea=(27237.8, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.8e-57, 'cm^3/(mol*s)'),
            n = 23.463,
            Ea = (-602.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -2.46,
        T3 = (206, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC2H4-4-1C8H15 <=> S4XC8H15 + C2H4""",
)

entry(
    index = 1540,
    label = "PXCH2-5-1C9H17 <=> CH3-5-SAX1C9H16",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(155, 's^-1'), n=2.83, Ea=(15566.2, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.5e-30, 'cm^3/(mol*s)'),
            n = 14.56,
            Ea = (-602.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -13.59,
        T3 = (214, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXCH2-5-1C9H17 <=> CH3-5-SAX1C9H16""",
)

entry(
    index = 1541,
    label = "PXCH2-5-1C9H17 <=> C6H10-15 + pC4H9",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.76e+11, 's^-1'), n=0.57, Ea=(28791.6, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4.7e-39, 'cm^3/(mol*s)'),
            n = 16.77,
            Ea = (-603.1, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -27.89,
        T3 = (216, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXCH2-5-1C9H17 <=> C6H10-15 + pC4H9""",
)

entry(
    index = 1542,
    label = "PXCH2-5-1C9H17 <=> C4H7 + C6H12",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.76e+11, 's^-1'), n=0.57, Ea=(28791.6, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4.7e-39, 'cm^3/(mol*s)'),
            n = 16.77,
            Ea = (-603.1, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -27.89,
        T3 = (216, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXCH2-5-1C9H17 <=> C4H7 + C6H12""",
)

entry(
    index = 1543,
    label = "SAX6-4C10H19 <=> nC3H7 + C7H12-13",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.39e+11, 's^-1'), n=0.66, Ea=(32262.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4e-42, 'cm^3/(mol*s)'),
            n = 18.045,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -18.5,
        T3 = (246, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SAX6-4C10H19 <=> nC3H7 + C7H12-13""",
)

entry(
    index = 1544,
    label = "SAX6-4C10H19 <=> C2H5 + C8H14-13",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.39e+11, 's^-1'), n=0.66, Ea=(32262.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4e-42, 'cm^3/(mol*s)'),
            n = 18.045,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -18.5,
        T3 = (246, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SAX6-4C10H19 <=> C2H5 + C8H14-13""",
)

entry(
    index = 1545,
    label = "SAXC4H8-2-1C6H11 <=> CH2-3-1C7H12 + C2H5",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.39e+11, 's^-1'), n=0.66, Ea=(32262.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4e-42, 'cm^3/(mol*s)'),
            n = 18.045,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -18.5,
        T3 = (246, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SAXC4H8-2-1C6H11 <=> CH2-3-1C7H12 + C2H5""",
)

entry(
    index = 1546,
    label = "SAX4-5C10H19 <=> C2H5 + C8H14-13",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.39e+11, 's^-1'), n=0.66, Ea=(32262.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4e-42, 'cm^3/(mol*s)'),
            n = 18.045,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -18.5,
        T3 = (246, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SAX4-5C10H19 <=> C2H5 + C8H14-13""",
)

entry(
    index = 1547,
    label = "C3H7-3-TAX1C7H13 <=> CH2-3-1C7H12 + C2H5",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.39e+11, 's^-1'), n=0.66, Ea=(32262.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4e-42, 'cm^3/(mol*s)'),
            n = 18.045,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -18.5,
        T3 = (246, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C3H7-3-TAX1C7H13 <=> CH2-3-1C7H12 + C2H5""",
)

entry(
    index = 1548,
    label = "C3H7-3-TAX1C7H13 <=> CH2-3-1C6H10 + nC3H7",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.39e+11, 's^-1'), n=0.66, Ea=(32262.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4e-42, 'cm^3/(mol*s)'),
            n = 18.045,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -18.5,
        T3 = (246, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C3H7-3-TAX1C7H13 <=> CH2-3-1C6H10 + nC3H7""",
)

entry(
    index = 1549,
    label = "SAXC10H19 <=> C4H6 + PXC6H13",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.39e+11, 's^-1'), n=0.66, Ea=(32262.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4e-42, 'cm^3/(mol*s)'),
            n = 18.045,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -18.5,
        T3 = (246, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SAXC10H19 <=> C4H6 + PXC6H13""",
)

entry(
    index = 1550,
    label = "C2H5-4-SAX1C8H14 <=> C8H14-13 + C2H5",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.39e+11, 's^-1'), n=0.66, Ea=(32262.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4e-42, 'cm^3/(mol*s)'),
            n = 18.045,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -18.5,
        T3 = (246, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C2H5-4-SAX1C8H14 <=> C8H14-13 + C2H5""",
)

entry(
    index = 1551,
    label = "C2H5-4-SAX1C8H14 <=> C6H10-13 + pC4H9",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.39e+11, 's^-1'), n=0.66, Ea=(31262.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4e-42, 'cm^3/(mol*s)'),
            n = 18.045,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -18.5,
        T3 = (246, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C2H5-4-SAX1C8H14 <=> C6H10-13 + pC4H9""",
)

entry(
    index = 1552,
    label = "CH3-5-SAX1C9H16 <=> C4H6 + SXC6H13",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.39e+11, 's^-1'), n=0.66, Ea=(31262.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4e-42, 'cm^3/(mol*s)'),
            n = 18.045,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -18.5,
        T3 = (246, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is CH3-5-SAX1C9H16 <=> C4H6 + SXC6H13""",
)

entry(
    index = 1553,
    label = "PXC2H4cC6H11 + CH3 <=> C3H7cC6H11",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.93e+14, 'cm^3/(mol*s)'), n=-0.32, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC2H4cC6H11 + CH3 <=> C3H7cC6H11""",
)

entry(
    index = 1554,
    label = "PXCH2cC6H11 + C2H5 <=> C3H7cC6H11",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.88e+14, 'cm^3/(mol*s)'), n=-0.5, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXCH2cC6H11 + C2H5 <=> C3H7cC6H11""",
)

entry(
    index = 1555,
    label = "cC6H11 + nC3H7 <=> C3H7cC6H11",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.29e+14, 'cm^3/(mol*s)'), n=-0.35, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC6H11 + nC3H7 <=> C3H7cC6H11""",
)

entry(
    index = 1556,
    label = "C3H7cC6H11 <=> C9H18-4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.17e+15, 's^-1'), n=0, Ea=(74000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7cC6H11 <=> C9H18-4""",
)

entry(
    index = 1557,
    label = "C3H7cC6H11 <=> C3H7-2-1C6H11",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.17e+15, 's^-1'), n=0, Ea=(74000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7cC6H11 <=> C3H7-2-1C6H11""",
)

entry(
    index = 1558,
    label = "C9H18-4 <=> SAXC7H13 + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.07e+23, 's^-1'), n=-2.03, Ea=(74958, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C9H18-4 <=> SAXC7H13 + C2H5""",
)

entry(
    index = 1559,
    label = "C9H18-4 <=> SAXC6H11 + nC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.07e+23, 's^-1'), n=-2.03, Ea=(74958, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C9H18-4 <=> SAXC6H11 + nC3H7""",
)

entry(
    index = 1560,
    label = "C9H18-4 <=> C6H12 + C3H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.42e+07, 's^-1'), n=1.65, Ea=(53752, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C9H18-4 <=> C6H12 + C3H6""",
)

entry(
    index = 1561,
    label = "C9H18-4 + H <=> pC4H9 + C5H10",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8e+21, 'cm^3/(mol*s)'),
        n = -2.39,
        Ea = (11180, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C9H18-4 + H <=> pC4H9 + C5H10""",
)

entry(
    index = 1562,
    label = "C9H18-4 + H <=> C6H12 + nC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8e+21, 'cm^3/(mol*s)'),
        n = -2.39,
        Ea = (11180, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C9H18-4 + H <=> C6H12 + nC3H7""",
)

entry(
    index = 1563,
    label = "C3H7-2-1C6H11 <=> PAXCH2-2-1C6H11 + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.07e+23, 's^-1'), n=-2.03, Ea=(74958, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7-2-1C6H11 <=> PAXCH2-2-1C6H11 + C2H5""",
)

entry(
    index = 1564,
    label = "C3H7-2-1C6H11 <=> PAXCH2-2-1C5H9 + nC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.07e+23, 's^-1'), n=-2.03, Ea=(74958, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7-2-1C6H11 <=> PAXCH2-2-1C5H9 + nC3H7""",
)

entry(
    index = 1565,
    label = "C3H7-2-1C6H11 <=> CH3-2-1C5H9 + C3H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.42e+07, 's^-1'), n=1.65, Ea=(53752, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7-2-1C6H11 <=> CH3-2-1C5H9 + C3H6""",
)

entry(
    index = 1566,
    label = "C3H7-2-1C6H11 + H <=> CH3-2-1C5H9 + nC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8e+21, 'cm^3/(mol*s)'),
        n = -2.39,
        Ea = (11180, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C3H7-2-1C6H11 + H <=> CH3-2-1C5H9 + nC3H7""",
)

entry(
    index = 1567,
    label = "C3H7-2-1C6H11 + H <=> C6H12 + nC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8e+21, 'cm^3/(mol*s)'),
        n = -2.39,
        Ea = (11180, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C3H7-2-1C6H11 + H <=> C6H12 + nC3H7""",
)

entry(
    index = 1568,
    label = "C3H7cC6H11 + H <=> S2XC3H6cC6H11 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(380000, 'cm^3/(mol*s)'), n=2.44, Ea=(4688, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7cC6H11 + H <=> S2XC3H6cC6H11 + H2""",
)

entry(
    index = 1569,
    label = "C3H7cC6H11 + H <=> SXC3H6cC6H11 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(417000, 'cm^3/(mol*s)'), n=2.49, Ea=(4470, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7cC6H11 + H <=> SXC3H6cC6H11 + H2""",
)

entry(
    index = 1570,
    label = "C3H7cC6H11 + H <=> PXC3H6cC6H11 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(261000, 'cm^3/(mol*s)'), n=2.57, Ea=(6628, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7cC6H11 + H <=> PXC3H6cC6H11 + H2""",
)

entry(
    index = 1571,
    label = "C3H7cC6H11 + H <=> C3H7TXcC6H10 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(413000, 'cm^3/(mol*s)'), n=2.36, Ea=(2835, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7cC6H11 + H <=> C3H7TXcC6H10 + H2""",
)

entry(
    index = 1572,
    label = "C3H7cC6H11 + H <=> C3H7S2XcC6H10 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(392000, 'cm^3/(mol*s)'), n=2.53, Ea=(4353, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7cC6H11 + H <=> C3H7S2XcC6H10 + H2""",
)

entry(
    index = 1573,
    label = "C3H7cC6H11 + H <=> C3H7S3XcC6H10 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(303000, 'cm^3/(mol*s)'), n=2.58, Ea=(4116, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7cC6H11 + H <=> C3H7S3XcC6H10 + H2""",
)

entry(
    index = 1574,
    label = "C3H7cC6H11 + H <=> C3H7S4XcC6H10 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(252000, 'cm^3/(mol*s)'), n=2.53, Ea=(4406, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7cC6H11 + H <=> C3H7S4XcC6H10 + H2""",
)

entry(
    index = 1575,
    label = "C3H7cC6H11 + O <=> PXC3H6cC6H11 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.2e+06, 'cm^3/(mol*s)'), n=2.4, Ea=(5504, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7cC6H11 + O <=> PXC3H6cC6H11 + OH""",
)

entry(
    index = 1576,
    label = "C3H7cC6H11 + O <=> SXC3H6cC6H11 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(310000, 'cm^3/(mol*s)'), n=2.5, Ea=(2225, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7cC6H11 + O <=> SXC3H6cC6H11 + OH""",
)

entry(
    index = 1577,
    label = "C3H7cC6H11 + O <=> S2XC3H6cC6H11 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(310000, 'cm^3/(mol*s)'), n=2.5, Ea=(2225, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7cC6H11 + O <=> S2XC3H6cC6H11 + OH""",
)

entry(
    index = 1578,
    label = "C3H7cC6H11 + O <=> C3H7TXcC6H10 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(157000, 'cm^3/(mol*s)'), n=2.5, Ea=(1110, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7cC6H11 + O <=> C3H7TXcC6H10 + OH""",
)

entry(
    index = 1579,
    label = "C3H7cC6H11 + O <=> C3H7S2XcC6H10 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(95200, 'cm^3/(mol*s)'), n=2.71, Ea=(2106, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7cC6H11 + O <=> C3H7S2XcC6H10 + OH""",
)

entry(
    index = 1580,
    label = "C3H7cC6H11 + O <=> C3H7S3XcC6H10 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(95200, 'cm^3/(mol*s)'), n=2.71, Ea=(2106, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7cC6H11 + O <=> C3H7S3XcC6H10 + OH""",
)

entry(
    index = 1581,
    label = "C3H7cC6H11 + O <=> C3H7S4XcC6H10 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(47600, 'cm^3/(mol*s)'), n=2.71, Ea=(2106, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7cC6H11 + O <=> C3H7S4XcC6H10 + OH""",
)

entry(
    index = 1582,
    label = "C3H7cC6H11 + OH <=> S2XC3H6cC6H11 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3630, 'cm^3/(mol*s)'), n=2.81, Ea=(-1585, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7cC6H11 + OH <=> S2XC3H6cC6H11 + H2O""",
)

entry(
    index = 1583,
    label = "C3H7cC6H11 + OH <=> SXC3H6cC6H11 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4140, 'cm^3/(mol*s)'), n=2.85, Ea=(-1084, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7cC6H11 + OH <=> SXC3H6cC6H11 + H2O""",
)

entry(
    index = 1584,
    label = "C3H7cC6H11 + OH <=> PXC3H6cC6H11 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6280, 'cm^3/(mol*s)'), n=2.9, Ea=(426, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7cC6H11 + OH <=> PXC3H6cC6H11 + H2O""",
)

entry(
    index = 1585,
    label = "C3H7cC6H11 + OH <=> C3H7TXcC6H10 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1430, 'cm^3/(mol*s)'), n=2.92, Ea=(-2601, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7cC6H11 + OH <=> C3H7TXcC6H10 + H2O""",
)

entry(
    index = 1586,
    label = "C3H7cC6H11 + OH <=> C3H7S2XcC6H10 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7790, 'cm^3/(mol*s)'), n=2.85, Ea=(-1764, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7cC6H11 + OH <=> C3H7S2XcC6H10 + H2O""",
)

entry(
    index = 1587,
    label = "C3H7cC6H11 + OH <=> C3H7S3XcC6H10 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9940, 'cm^3/(mol*s)'), n=2.86, Ea=(-1139, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7cC6H11 + OH <=> C3H7S3XcC6H10 + H2O""",
)

entry(
    index = 1588,
    label = "C3H7cC6H11 + OH <=> C3H7S4XcC6H10 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6700, 'cm^3/(mol*s)'), n=2.86, Ea=(-1149, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7cC6H11 + OH <=> C3H7S4XcC6H10 + H2O""",
)

entry(
    index = 1589,
    label = "C3H7cC6H11 + O2 <=> PXC3H6cC6H11 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(50930, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7cC6H11 + O2 <=> PXC3H6cC6H11 + HO2""",
)

entry(
    index = 1590,
    label = "C3H7cC6H11 + O2 <=> SXC3H6cC6H11 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(47590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7cC6H11 + O2 <=> SXC3H6cC6H11 + HO2""",
)

entry(
    index = 1591,
    label = "C3H7cC6H11 + O2 <=> S2XC3H6cC6H11 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(47590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7cC6H11 + O2 <=> S2XC3H6cC6H11 + HO2""",
)

entry(
    index = 1592,
    label = "C3H7cC6H11 + O2 <=> C3H7TXcC6H10 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(44000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7cC6H11 + O2 <=> C3H7TXcC6H10 + HO2""",
)

entry(
    index = 1593,
    label = "C3H7cC6H11 + O2 <=> C3H7S2XcC6H10 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+13, 'cm^3/(mol*s)'), n=0, Ea=(47590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7cC6H11 + O2 <=> C3H7S2XcC6H10 + HO2""",
)

entry(
    index = 1594,
    label = "C3H7cC6H11 + O2 <=> C3H7S3XcC6H10 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+13, 'cm^3/(mol*s)'), n=0, Ea=(47590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7cC6H11 + O2 <=> C3H7S3XcC6H10 + HO2""",
)

entry(
    index = 1595,
    label = "C3H7cC6H11 + O2 <=> C3H7S4XcC6H10 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(47590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7cC6H11 + O2 <=> C3H7S4XcC6H10 + HO2""",
)

entry(
    index = 1596,
    label = "C3H7cC6H11 + HO2 <=> PXC3H6cC6H11 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(30500, 'cm^3/(mol*s)'), n=2.65, Ea=(17496, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7cC6H11 + HO2 <=> PXC3H6cC6H11 + H2O2""",
)

entry(
    index = 1597,
    label = "C3H7cC6H11 + HO2 <=> SXC3H6cC6H11 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7130, 'cm^3/(mol*s)'), n=2.77, Ea=(14913, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7cC6H11 + HO2 <=> SXC3H6cC6H11 + H2O2""",
)

entry(
    index = 1598,
    label = "C3H7cC6H11 + HO2 <=> S2XC3H6cC6H11 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7130, 'cm^3/(mol*s)'), n=2.77, Ea=(14913, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7cC6H11 + HO2 <=> S2XC3H6cC6H11 + H2O2""",
)

entry(
    index = 1599,
    label = "C3H7cC6H11 + HO2 <=> C3H7TXcC6H10 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1240, 'cm^3/(mol*s)'), n=2.77, Ea=(10500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7cC6H11 + HO2 <=> C3H7TXcC6H10 + H2O2""",
)

entry(
    index = 1600,
    label = "C3H7cC6H11 + HO2 <=> C3H7S2XcC6H10 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(281, 'cm^3/(mol*s)'), n=3.25, Ea=(14998, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7cC6H11 + HO2 <=> C3H7S2XcC6H10 + H2O2""",
)

entry(
    index = 1601,
    label = "C3H7cC6H11 + HO2 <=> C3H7S3XcC6H10 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(281, 'cm^3/(mol*s)'), n=3.25, Ea=(14998, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7cC6H11 + HO2 <=> C3H7S3XcC6H10 + H2O2""",
)

entry(
    index = 1602,
    label = "C3H7cC6H11 + HO2 <=> C3H7S4XcC6H10 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(140, 'cm^3/(mol*s)'), n=3.25, Ea=(14998, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7cC6H11 + HO2 <=> C3H7S4XcC6H10 + H2O2""",
)

entry(
    index = 1603,
    label = "C3H7cC6H11 + CH3 <=> S2XC3H6cC6H11 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(73.8, 'cm^3/(mol*s)'), n=3.3, Ea=(13603, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7cC6H11 + CH3 <=> S2XC3H6cC6H11 + CH4""",
)

entry(
    index = 1604,
    label = "C3H7cC6H11 + CH3 <=> SXC3H6cC6H11 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(49, 'cm^3/(mol*s)'), n=3.25, Ea=(11163, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7cC6H11 + CH3 <=> SXC3H6cC6H11 + CH4""",
)

entry(
    index = 1605,
    label = "C3H7cC6H11 + CH3 <=> PXC3H6cC6H11 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(53.2, 'cm^3/(mol*s)'), n=3.1, Ea=(11666, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7cC6H11 + CH3 <=> PXC3H6cC6H11 + CH4""",
)

entry(
    index = 1606,
    label = "C3H7cC6H11 + CH3 <=> C3H7TXcC6H10 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.34, 'cm^3/(mol*s)'), n=3.27, Ea=(8906, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7cC6H11 + CH3 <=> C3H7TXcC6H10 + CH4""",
)

entry(
    index = 1607,
    label = "C3H7cC6H11 + CH3 <=> C3H7S2XcC6H10 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(33.7, 'cm^3/(mol*s)'), n=3.36, Ea=(10569, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7cC6H11 + CH3 <=> C3H7S2XcC6H10 + CH4""",
)

entry(
    index = 1608,
    label = "C3H7cC6H11 + CH3 <=> C3H7S3XcC6H10 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(144, 'cm^3/(mol*s)'), n=3.21, Ea=(11387, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7cC6H11 + CH3 <=> C3H7S3XcC6H10 + CH4""",
)

entry(
    index = 1609,
    label = "C3H7cC6H11 + CH3 <=> C3H7S4XcC6H10 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(30, 'cm^3/(mol*s)'), n=3.28, Ea=(11233, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7cC6H11 + CH3 <=> C3H7S4XcC6H10 + CH4""",
)

entry(
    index = 1610,
    label = "PXC3H6cC6H11 <=> PXCH2cC6H11 + C2H4",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(9.12e+11, 's^-1'), n=0.31, Ea=(27237.8, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.8e-57, 'cm^3/(mol*s)'),
            n = 23.463,
            Ea = (-602.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -2.46,
        T3 = (206, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC3H6cC6H11 <=> PXCH2cC6H11 + C2H4""",
)

entry(
    index = 1611,
    label = "SXC3H6cC6H11 <=> cC6H11 + C3H6",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(6.03e+10, 's^-1'), n=0.84, Ea=(27820, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1e-43, 'cm^3/(mol*s)'),
            n = 18.591,
            Ea = (-602.5, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -43.32,
        T3 = (200, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SXC3H6cC6H11 <=> cC6H11 + C3H6""",
)

entry(
    index = 1612,
    label = "S2XC3H6cC6H11 <=> PX9-3C9H17",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(6.03e+12, 's^-1'), n=0.07, Ea=(26982.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.1e-33, 'cm^3/(mol*s)'),
            n = 15.29,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -25.11,
        T3 = (225, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S2XC3H6cC6H11 <=> PX9-3C9H17""",
)

entry(
    index = 1613,
    label = "S2XC3H6cC6H11 <=> C2H3cC6H11 + CH3",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(9.55e+09, 's^-1'), n=1.08, Ea=(29387.7, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.3e-46, 'cm^3/(mol*s)'),
            n = 19.133,
            Ea = (-602.7, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -34.36,
        T3 = (210, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S2XC3H6cC6H11 <=> C2H3cC6H11 + CH3""",
)

entry(
    index = 1614,
    label = "C3H7S2XcC6H10 <=> cC6H10 + nC3H7",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(5.47e+11, 's^-1'), n=0.5, Ea=(27798.1, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (8.5e-37, 'cm^3/(mol*s)'),
            n = 16.21,
            Ea = (-603.1, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -26.7,
        T3 = (216, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C3H7S2XcC6H10 <=> cC6H10 + nC3H7""",
)

entry(
    index = 1615,
    label = "C3H7S2XcC6H10 <=> PX9-4C9H17",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.01e+12, 's^-1'), n=0.07, Ea=(26982.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.1e-33, 'cm^3/(mol*s)'),
            n = 15.29,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -25.11,
        T3 = (225, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C3H7S2XcC6H10 <=> PX9-4C9H17""",
)

entry(
    index = 1616,
    label = "C3H7S2XcC6H10 <=> PXC3H6-3-1C6H11",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.01e+12, 's^-1'), n=0.07, Ea=(27982.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.1e-33, 'cm^3/(mol*s)'),
            n = 15.29,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -25.11,
        T3 = (225, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C3H7S2XcC6H10 <=> PXC3H6-3-1C6H11""",
)

entry(
    index = 1617,
    label = "C3H7TXcC6H10 <=> CH2cC6H10 + C2H5",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.63e+11, 's^-1'), n=0.68, Ea=(24579.2, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.1e-39, 'cm^3/(mol*s)'),
            n = 17.57,
            Ea = (-603.1, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -31.51,
        T3 = (219, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C3H7TXcC6H10 <=> CH2cC6H10 + C2H5""",
)

entry(
    index = 1618,
    label = "C3H7TXcC6H10 <=> C3H7-2-PXC6H10",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(6.03e+12, 's^-1'), n=0.07, Ea=(27982.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.1e-33, 'cm^3/(mol*s)'),
            n = 15.29,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -25.11,
        T3 = (225, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C3H7TXcC6H10 <=> C3H7-2-PXC6H10""",
)

entry(
    index = 1619,
    label = "C3H7S3XcC6H10 <=> S3XC9H17",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.01e+12, 's^-1'), n=0.07, Ea=(26982.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.1e-33, 'cm^3/(mol*s)'),
            n = 15.29,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -25.11,
        T3 = (225, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C3H7S3XcC6H10 <=> S3XC9H17""",
)

entry(
    index = 1620,
    label = "C3H7S3XcC6H10 <=> PXC2H4-4-1C7H13",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.01e+12, 's^-1'), n=0.07, Ea=(27982.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.1e-33, 'cm^3/(mol*s)'),
            n = 15.29,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -25.11,
        T3 = (225, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C3H7S3XcC6H10 <=> PXC2H4-4-1C7H13""",
)

entry(
    index = 1621,
    label = "C3H7S4XcC6H10 <=> PXCH2-5-1C8H15",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(6.03e+12, 's^-1'), n=0.07, Ea=(27982.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.1e-33, 'cm^3/(mol*s)'),
            n = 15.29,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -25.11,
        T3 = (225, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C3H7S4XcC6H10 <=> PXCH2-5-1C8H15""",
)

entry(
    index = 1622,
    label = "PXC3H6cC6H11 <=> C3H7S2XcC6H10",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.05e+09, 's^-1'), n=0.82, Ea=(15027, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC3H6cC6H11 <=> C3H7S2XcC6H10""",
)

entry(
    index = 1623,
    label = "PXC3H6cC6H11 <=> C3H7TXcC6H10",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.57e+08, 's^-1'), n=0.93, Ea=(17967, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC3H6cC6H11 <=> C3H7TXcC6H10""",
)

entry(
    index = 1624,
    label = "PXC3H6cC6H11 <=> C3H7S3XcC6H10",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.95e+08, 's^-1'), n=0.79, Ea=(18372, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC3H6cC6H11 <=> C3H7S3XcC6H10""",
)

entry(
    index = 1625,
    label = "PXC3H6cC6H11 <=> C3H7S4XcC6H10",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.75e+08, 's^-1'), n=0.78, Ea=(24300, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC3H6cC6H11 <=> C3H7S4XcC6H10""",
)

entry(
    index = 1626,
    label = "SXC3H6cC6H11 <=> C3H7S2XcC6H10",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.58e+08, 's^-1'), n=0.9, Ea=(22900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC3H6cC6H11 <=> C3H7S2XcC6H10""",
)

entry(
    index = 1627,
    label = "SXC3H6cC6H11 <=> C3H7S3XcC6H10",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.94e+08, 's^-1'), n=0.78, Ea=(18644, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC3H6cC6H11 <=> C3H7S3XcC6H10""",
)

entry(
    index = 1628,
    label = "SXC3H6cC6H11 <=> C3H7S4XcC6H10",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.97e+08, 's^-1'), n=0.78, Ea=(22740, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC3H6cC6H11 <=> C3H7S4XcC6H10""",
)

entry(
    index = 1629,
    label = "S2XC3H6cC6H11 <=> C3H7S3XcC6H10",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.65e+08, 's^-1'), n=1.02, Ea=(28687, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is S2XC3H6cC6H11 <=> C3H7S3XcC6H10""",
)

entry(
    index = 1630,
    label = "S2XC3H6cC6H11 <=> C3H7S4XcC6H10",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.19e+09, 's^-1'), n=0.92, Ea=(22700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is S2XC3H6cC6H11 <=> C3H7S4XcC6H10""",
)

entry(
    index = 1631,
    label = "PXC3H6cC6H11 + H <=> C3H7cC6H11",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.5e+48, 'cm^6/(mol^2*s)'),
            n = -9.32,
            Ea = (5833.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.498,
        T3 = (1314, 'K'),
        T1 = (1314, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC3H6cC6H11 + H <=> C3H7cC6H11""",
)

entry(
    index = 1632,
    label = "PXC3H6cC6H11 + H <=> cC6H11 + nC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.85e+24, 'cm^3/(mol*s)'),
        n = -2.92,
        Ea = (12505, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is PXC3H6cC6H11 + H <=> cC6H11 + nC3H7""",
)

entry(
    index = 1633,
    label = "PXC3H6cC6H11 + O <=> PXC2H4cC6H11 + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC3H6cC6H11 + O <=> PXC2H4cC6H11 + CH2O""",
)

entry(
    index = 1634,
    label = "SXC3H6cC6H11 + H <=> C3H7cC6H11",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(4.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.7e+58, 'cm^6/(mol^2*s)'),
            n = -12.08,
            Ea = (11263.7, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.649,
        T3 = (1213.1, 'K'),
        T1 = (1213.1, 'K'),
        T2 = (13369.7, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SXC3H6cC6H11 + H <=> C3H7cC6H11""",
)

entry(
    index = 1635,
    label = "SXC3H6cC6H11 + H <=> cC6H11 + nC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.8e+28, 'cm^3/(mol*s)'),
        n = -3.94,
        Ea = (15916, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is SXC3H6cC6H11 + H <=> cC6H11 + nC3H7""",
)

entry(
    index = 1636,
    label = "S2XC3H6cC6H11 + H <=> C3H7cC6H11",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(4.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.7e+58, 'cm^6/(mol^2*s)'),
            n = -12.08,
            Ea = (11263.7, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.649,
        T3 = (1213.1, 'K'),
        T1 = (1213.1, 'K'),
        T2 = (13369.7, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S2XC3H6cC6H11 + H <=> C3H7cC6H11""",
)

entry(
    index = 1637,
    label = "S2XC3H6cC6H11 + H <=> cC6H11 + nC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.8e+28, 'cm^3/(mol*s)'),
        n = -3.94,
        Ea = (15916, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is S2XC3H6cC6H11 + H <=> cC6H11 + nC3H7""",
)

entry(
    index = 1638,
    label = "C3H7S2XcC6H10 + H <=> C3H7cC6H11",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(4.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.7e+58, 'cm^6/(mol^2*s)'),
            n = -12.08,
            Ea = (11263.7, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.649,
        T3 = (1213.1, 'K'),
        T1 = (1213.1, 'K'),
        T2 = (13369.7, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C3H7S2XcC6H10 + H <=> C3H7cC6H11""",
)

entry(
    index = 1639,
    label = "C3H7S2XcC6H10 + H <=> cC6H11 + nC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.8e+28, 'cm^3/(mol*s)'),
        n = -3.94,
        Ea = (15916, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C3H7S2XcC6H10 + H <=> cC6H11 + nC3H7""",
)

entry(
    index = 1640,
    label = "C3H7TXcC6H10 + H <=> C3H7cC6H11",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.47e+61, 'cm^6/(mol^2*s)'),
            n = -12.94,
            Ea = (8000, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        T3 = (1456.4, 'K'),
        T1 = (1000, 'K'),
        T2 = (10000.5, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C3H7TXcC6H10 + H <=> C3H7cC6H11""",
)

entry(
    index = 1641,
    label = "C3H7TXcC6H10 + H <=> cC6H11 + nC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.6e+36, 'cm^3/(mol*s)'),
        n = -6.12,
        Ea = (25640, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C3H7TXcC6H10 + H <=> cC6H11 + nC3H7""",
)

entry(
    index = 1642,
    label = "C3H7S3XcC6H10 + H <=> C3H7cC6H11",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(4.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.7e+58, 'cm^6/(mol^2*s)'),
            n = -12.08,
            Ea = (11263.7, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.649,
        T3 = (1213.1, 'K'),
        T1 = (1213.1, 'K'),
        T2 = (13369.7, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C3H7S3XcC6H10 + H <=> C3H7cC6H11""",
)

entry(
    index = 1643,
    label = "C3H7S3XcC6H10 + H <=> cC6H11 + nC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.8e+28, 'cm^3/(mol*s)'),
        n = -3.94,
        Ea = (15916, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C3H7S3XcC6H10 + H <=> cC6H11 + nC3H7""",
)

entry(
    index = 1644,
    label = "C3H7S4XcC6H10 + H <=> C3H7cC6H11",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(4.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.7e+58, 'cm^6/(mol^2*s)'),
            n = -12.08,
            Ea = (11263.7, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.649,
        T3 = (1213.1, 'K'),
        T1 = (1213.1, 'K'),
        T2 = (13369.7, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C3H7S4XcC6H10 + H <=> C3H7cC6H11""",
)

entry(
    index = 1645,
    label = "C3H7S4XcC6H10 + H <=> cC6H11 + nC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.8e+28, 'cm^3/(mol*s)'),
        n = -3.94,
        Ea = (15916, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C3H7S4XcC6H10 + H <=> cC6H11 + nC3H7""",
)

entry(
    index = 1646,
    label = "PX9-3C9H17 <=> PX1-3C7H13 + C2H4",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(9.12e+11, 's^-1'), n=0.31, Ea=(27237.8, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.8e-57, 'cm^3/(mol*s)'),
            n = 23.463,
            Ea = (-602.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -2.46,
        T3 = (206, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PX9-3C9H17 <=> PX1-3C7H13 + C2H4""",
)

entry(
    index = 1647,
    label = "PX9-3C9H17 <=> SAX5-3C9H17",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.46e+11, 's^-1'), n=0, Ea=(10516.8, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (9.9e-38, 'cm^3/(mol*s)'),
            n = 17.215,
            Ea = (-603, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -16.33,
        T3 = (200, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PX9-3C9H17 <=> SAX5-3C9H17""",
)

entry(
    index = 1648,
    label = "PX9-4C9H17 <=> PX1-3C7H13 + C2H4",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(9.12e+11, 's^-1'), n=0.31, Ea=(27237.8, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.8e-57, 'cm^3/(mol*s)'),
            n = 23.463,
            Ea = (-602.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -2.46,
        T3 = (206, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PX9-4C9H17 <=> PX1-3C7H13 + C2H4""",
)

entry(
    index = 1649,
    label = "PX9-4C9H17 <=> SAX6-4C9H17",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(155, 's^-1'), n=2.83, Ea=(15566.2, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.5e-30, 'cm^3/(mol*s)'),
            n = 14.56,
            Ea = (-602.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -13.59,
        T3 = (214, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PX9-4C9H17 <=> SAX6-4C9H17""",
)

entry(
    index = 1650,
    label = "PXC3H6-3-1C6H11 <=> PXCH2-3-1C6H11 + C2H4",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(9.12e+11, 's^-1'), n=0.31, Ea=(27237.8, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.8e-57, 'cm^3/(mol*s)'),
            n = 23.463,
            Ea = (-602.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -2.46,
        T3 = (206, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC3H6-3-1C6H11 <=> PXCH2-3-1C6H11 + C2H4""",
)

entry(
    index = 1651,
    label = "PXC3H6-3-1C6H11 <=> C3H7-3-TAX1C6H10",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(86, 's^-1'), n=2.62, Ea=(8722.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (8.1e-33, 'cm^3/(mol*s)'),
            n = 15.214,
            Ea = (-677.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -30.39,
        T3 = (206, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC3H6-3-1C6H11 <=> C3H7-3-TAX1C6H10""",
)

entry(
    index = 1652,
    label = "C3H7-2-PXC6H10 <=> PXC2H4-2-1C5H9 + C2H4",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(9.12e+11, 's^-1'), n=0.31, Ea=(27237.8, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.8e-57, 'cm^3/(mol*s)'),
            n = 23.463,
            Ea = (-602.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -2.46,
        T3 = (206, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C3H7-2-PXC6H10 <=> PXC2H4-2-1C5H9 + C2H4""",
)

entry(
    index = 1653,
    label = "C3H7-2-PXC6H10 <=> C3H7-2-SAXC6H10",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.46e+11, 's^-1'), n=0, Ea=(10516.8, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (9.9e-38, 'cm^3/(mol*s)'),
            n = 17.215,
            Ea = (-603, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -16.33,
        T3 = (200, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C3H7-2-PXC6H10 <=> C3H7-2-SAXC6H10""",
)

entry(
    index = 1654,
    label = "S3XC9H17 <=> C7H12-16 + C2H5",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(5.5e+11, 's^-1'), n=0.55, Ea=(28084.3, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.1e-43, 'cm^3/(mol*s)'),
            n = 18.418,
            Ea = (-602.9, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -32.13,
        T3 = (207, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S3XC9H17 <=> C7H12-16 + C2H5""",
)

entry(
    index = 1655,
    label = "S3XC9H17 <=> C5H10 + C4H7",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(5.5e+11, 's^-1'), n=0.55, Ea=(28084.3, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.1e-43, 'cm^3/(mol*s)'),
            n = 18.418,
            Ea = (-602.9, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -32.13,
        T3 = (207, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S3XC9H17 <=> C5H10 + C4H7""",
)

entry(
    index = 1656,
    label = "S3XC9H17 <=> SAXC9H17",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(155, 's^-1'), n=2.83, Ea=(15566.2, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.5e-30, 'cm^3/(mol*s)'),
            n = 14.56,
            Ea = (-602.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -13.59,
        T3 = (214, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S3XC9H17 <=> SAXC9H17""",
)

entry(
    index = 1657,
    label = "PXC2H4-4-1C7H13 <=> S3XC7H13 + C2H4",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(9.12e+11, 's^-1'), n=0.31, Ea=(27237.8, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.8e-57, 'cm^3/(mol*s)'),
            n = 23.463,
            Ea = (-602.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -2.46,
        T3 = (206, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC2H4-4-1C7H13 <=> S3XC7H13 + C2H4""",
)

entry(
    index = 1658,
    label = "PXC2H4-4-1C7H13 <=> C2H5-4-SAX1C7H12",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(155, 's^-1'), n=2.83, Ea=(15566.2, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.5e-30, 'cm^3/(mol*s)'),
            n = 14.56,
            Ea = (-602.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -13.59,
        T3 = (214, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC2H4-4-1C7H13 <=> C2H5-4-SAX1C7H12""",
)

entry(
    index = 1659,
    label = "PXCH2-5-1C8H15 <=> CH3-5-SAX1C8H14",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(155, 's^-1'), n=2.83, Ea=(15566.2, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.5e-30, 'cm^3/(mol*s)'),
            n = 14.56,
            Ea = (-602.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -13.59,
        T3 = (214, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXCH2-5-1C8H15 <=> CH3-5-SAX1C8H14""",
)

entry(
    index = 1660,
    label = "PXCH2-5-1C8H15 <=> C5H10 + C4H7",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.76e+11, 's^-1'), n=0.57, Ea=(28791.6, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4.7e-39, 'cm^3/(mol*s)'),
            n = 16.77,
            Ea = (-603.1, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -27.89,
        T3 = (216, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXCH2-5-1C8H15 <=> C5H10 + C4H7""",
)

entry(
    index = 1661,
    label = "PXCH2-5-1C8H15 <=> C6H10-15 + nC3H7",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.76e+11, 's^-1'), n=0.57, Ea=(28791.6, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4.7e-39, 'cm^3/(mol*s)'),
            n = 16.77,
            Ea = (-603.1, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -27.89,
        T3 = (216, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXCH2-5-1C8H15 <=> C6H10-15 + nC3H7""",
)

entry(
    index = 1662,
    label = "SAX5-3C9H17 <=> C6H10-13 + nC3H7",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.39e+11, 's^-1'), n=0.66, Ea=(32262.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4e-42, 'cm^3/(mol*s)'),
            n = 18.045,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -18.5,
        T3 = (246, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SAX5-3C9H17 <=> C6H10-13 + nC3H7""",
)

entry(
    index = 1663,
    label = "SAX6-4C9H17 <=> C7H12-13 + C2H5",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.39e+11, 's^-1'), n=0.66, Ea=(32262.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4e-42, 'cm^3/(mol*s)'),
            n = 18.045,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -18.5,
        T3 = (246, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SAX6-4C9H17 <=> C7H12-13 + C2H5""",
)

entry(
    index = 1664,
    label = "C3H7-3-TAX1C6H10 <=> CH2-3-1C6H10 + C2H5",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.39e+11, 's^-1'), n=0.66, Ea=(32262.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4e-42, 'cm^3/(mol*s)'),
            n = 18.045,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -18.5,
        T3 = (246, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C3H7-3-TAX1C6H10 <=> CH2-3-1C6H10 + C2H5""",
)

entry(
    index = 1665,
    label = "C3H7-2-SAXC6H10 <=> CH2-3-1C6H10 + C2H5",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.39e+11, 's^-1'), n=0.66, Ea=(32262.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4e-42, 'cm^3/(mol*s)'),
            n = 18.045,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -18.5,
        T3 = (246, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C3H7-2-SAXC6H10 <=> CH2-3-1C6H10 + C2H5""",
)

entry(
    index = 1666,
    label = "SAXC9H17 <=> PXC5H11 + C4H6",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.39e+11, 's^-1'), n=0.66, Ea=(32262.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4e-42, 'cm^3/(mol*s)'),
            n = 18.045,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -18.5,
        T3 = (246, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SAXC9H17 <=> PXC5H11 + C4H6""",
)

entry(
    index = 1667,
    label = "C2H5-4-SAX1C7H12 <=> C7H12-13 + C2H5",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.39e+11, 's^-1'), n=0.66, Ea=(32262.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4e-42, 'cm^3/(mol*s)'),
            n = 18.045,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -18.5,
        T3 = (246, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C2H5-4-SAX1C7H12 <=> C7H12-13 + C2H5""",
)

entry(
    index = 1668,
    label = "C2H5-4-SAX1C7H12 <=> C6H10-12 + nC3H7",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.39e+11, 's^-1'), n=0.66, Ea=(32262.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4e-42, 'cm^3/(mol*s)'),
            n = 18.045,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -18.5,
        T3 = (246, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C2H5-4-SAX1C7H12 <=> C6H10-12 + nC3H7""",
)

entry(
    index = 1669,
    label = "CH3-5-SAX1C8H14 <=> PXC5H11 + C4H6",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.39e+11, 's^-1'), n=0.66, Ea=(32262.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4e-42, 'cm^3/(mol*s)'),
            n = 18.045,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -18.5,
        T3 = (246, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is CH3-5-SAX1C8H14 <=> PXC5H11 + C4H6""",
)

entry(
    index = 1670,
    label = "C3H5cC6H11 + H <=> C2H4 + PXCH2cC6H11",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8e+21, 'cm^3/(mol*s)'),
        n = -2.39,
        Ea = (11180, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C3H5cC6H11 + H <=> C2H4 + PXCH2cC6H11""",
)

entry(
    index = 1671,
    label = "C3H5cC6H11 + O <=> PXC2H4cC6H11 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.3e+08, 'cm^3/(mol*s)'),
        n = 1.45,
        Ea = (-402, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C3H5cC6H11 + O <=> PXC2H4cC6H11 + HCO""",
)

entry(
    index = 1672,
    label = "aC3H5 + cC6H11 <=> C3H5cC6H11",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.15e+14, 'cm^3/(mol*s)'), n=-0.35, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is aC3H5 + cC6H11 <=> C3H5cC6H11""",
)

entry(
    index = 1673,
    label = "PXCH2cC6H11 + CH3 <=> C2H5cC6H11",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.93e+14, 'cm^3/(mol*s)'), n=-0.32, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXCH2cC6H11 + CH3 <=> C2H5cC6H11""",
)

entry(
    index = 1674,
    label = "cC6H11 + C2H5 <=> C2H5cC6H11",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.95e+14, 'cm^3/(mol*s)'), n=-0.35, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC6H11 + C2H5 <=> C2H5cC6H11""",
)

entry(
    index = 1675,
    label = "C2H5cC6H11 <=> C8H16-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.17e+15, 's^-1'), n=0, Ea=(74000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5cC6H11 <=> C8H16-3""",
)

entry(
    index = 1676,
    label = "C2H5cC6H11 <=> C2H5-2-1C6H11",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.17e+15, 's^-1'), n=0, Ea=(74000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5cC6H11 <=> C2H5-2-1C6H11""",
)

entry(
    index = 1677,
    label = "C8H16-3 <=> SAXC5H9 + nC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.07e+23, 's^-1'), n=-2.03, Ea=(74958, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C8H16-3 <=> SAXC5H9 + nC3H7""",
)

entry(
    index = 1678,
    label = "C8H16-3 <=> C5H10 + C3H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.42e+07, 's^-1'), n=1.65, Ea=(53752, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C8H16-3 <=> C5H10 + C3H6""",
)

entry(
    index = 1679,
    label = "C8H16-3 + H <=> C5H10 + nC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8e+21, 'cm^3/(mol*s)'),
        n = -2.39,
        Ea = (11180, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C8H16-3 + H <=> C5H10 + nC3H7""",
)

entry(
    index = 1680,
    label = "C8H16-3 + H <=> C4H81 + pC4H9",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8e+21, 'cm^3/(mol*s)'),
        n = -2.39,
        Ea = (11180, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C8H16-3 + H <=> C4H81 + pC4H9""",
)

entry(
    index = 1681,
    label = "C2H5-2-1C6H11 <=> PAXCH2-2-1C4H7 + nC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.07e+23, 's^-1'), n=-2.03, Ea=(74958, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5-2-1C6H11 <=> PAXCH2-2-1C4H7 + nC3H7""",
)

entry(
    index = 1682,
    label = "C2H5-2-1C6H11 <=> C3H6 + CH3-2-1C4H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.54e+06, 's^-1'), n=1.65, Ea=(53752, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5-2-1C6H11 <=> C3H6 + CH3-2-1C4H7""",
)

entry(
    index = 1683,
    label = "C2H5-2-1C6H11 + H <=> C6H12 + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8e+21, 'cm^3/(mol*s)'),
        n = -2.39,
        Ea = (11180, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H5-2-1C6H11 + H <=> C6H12 + C2H5""",
)

entry(
    index = 1684,
    label = "C2H5-2-1C6H11 + H <=> nC3H7 + CH3-2-1C4H7",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.6e+22, 'cm^3/(mol*s)'),
        n = -2.39,
        Ea = (11180, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H5-2-1C6H11 + H <=> nC3H7 + CH3-2-1C4H7""",
)

entry(
    index = 1685,
    label = "C2H5cC6H11 + H <=> SXC2H4cC6H11 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(255000, 'cm^3/(mol*s)'), n=2.53, Ea=(4479, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5cC6H11 + H <=> SXC2H4cC6H11 + H2""",
)

entry(
    index = 1686,
    label = "C2H5cC6H11 + H <=> PXC2H4cC6H11 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(389000, 'cm^3/(mol*s)'), n=2.51, Ea=(6827, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5cC6H11 + H <=> PXC2H4cC6H11 + H2""",
)

entry(
    index = 1687,
    label = "C2H5cC6H11 + H <=> C2H5TXcC6H10 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(522000, 'cm^3/(mol*s)'), n=2.35, Ea=(2912, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5cC6H11 + H <=> C2H5TXcC6H10 + H2""",
)

entry(
    index = 1688,
    label = "C2H5cC6H11 + H <=> C2H5S2XcC6H10 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(310000, 'cm^3/(mol*s)'), n=2.58, Ea=(4216, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5cC6H11 + H <=> C2H5S2XcC6H10 + H2""",
)

entry(
    index = 1689,
    label = "C2H5cC6H11 + H <=> C2H5S3XcC6H10 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(207000, 'cm^3/(mol*s)'), n=2.63, Ea=(3929, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5cC6H11 + H <=> C2H5S3XcC6H10 + H2""",
)

entry(
    index = 1690,
    label = "C2H5cC6H11 + H <=> C2H5S4XcC6H10 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(207000, 'cm^3/(mol*s)'), n=2.55, Ea=(4330, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5cC6H11 + H <=> C2H5S4XcC6H10 + H2""",
)

entry(
    index = 1691,
    label = "C2H5cC6H11 + O <=> PXC2H4cC6H11 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.2e+06, 'cm^3/(mol*s)'), n=2.4, Ea=(5504, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5cC6H11 + O <=> PXC2H4cC6H11 + OH""",
)

entry(
    index = 1692,
    label = "C2H5cC6H11 + O <=> SXC2H4cC6H11 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(310000, 'cm^3/(mol*s)'), n=2.5, Ea=(2225, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5cC6H11 + O <=> SXC2H4cC6H11 + OH""",
)

entry(
    index = 1693,
    label = "C2H5cC6H11 + O <=> C2H5TXcC6H10 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(157000, 'cm^3/(mol*s)'), n=2.5, Ea=(1110, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5cC6H11 + O <=> C2H5TXcC6H10 + OH""",
)

entry(
    index = 1694,
    label = "C2H5cC6H11 + O <=> C2H5S2XcC6H10 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(95200, 'cm^3/(mol*s)'), n=2.71, Ea=(2106, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5cC6H11 + O <=> C2H5S2XcC6H10 + OH""",
)

entry(
    index = 1695,
    label = "C2H5cC6H11 + O <=> C2H5S3XcC6H10 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(95200, 'cm^3/(mol*s)'), n=2.71, Ea=(2106, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5cC6H11 + O <=> C2H5S3XcC6H10 + OH""",
)

entry(
    index = 1696,
    label = "C2H5cC6H11 + O <=> C2H5S4XcC6H10 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(47600, 'cm^3/(mol*s)'), n=2.71, Ea=(2106, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5cC6H11 + O <=> C2H5S4XcC6H10 + OH""",
)

entry(
    index = 1697,
    label = "C2H5cC6H11 + OH <=> SXC2H4cC6H11 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3050, 'cm^3/(mol*s)'), n=2.85, Ea=(-1229, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5cC6H11 + OH <=> SXC2H4cC6H11 + H2O""",
)

entry(
    index = 1698,
    label = "C2H5cC6H11 + OH <=> PXC2H4cC6H11 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8820, 'cm^3/(mol*s)'), n=2.9, Ea=(273, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5cC6H11 + OH <=> PXC2H4cC6H11 + H2O""",
)

entry(
    index = 1699,
    label = "C2H5cC6H11 + OH <=> C2H5TXcC6H10 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1410, 'cm^3/(mol*s)'), n=2.91, Ea=(-2538, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5cC6H11 + OH <=> C2H5TXcC6H10 + H2O""",
)

entry(
    index = 1700,
    label = "C2H5cC6H11 + OH <=> C2H5S2XcC6H10 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9230, 'cm^3/(mol*s)'), n=2.84, Ea=(-1603, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5cC6H11 + OH <=> C2H5S2XcC6H10 + H2O""",
)

entry(
    index = 1701,
    label = "C2H5cC6H11 + OH <=> C2H5S3XcC6H10 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(11100, 'cm^3/(mol*s)'), n=2.84, Ea=(-1055, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5cC6H11 + OH <=> C2H5S3XcC6H10 + H2O""",
)

entry(
    index = 1702,
    label = "C2H5cC6H11 + OH <=> C2H5S4XcC6H10 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4540, 'cm^3/(mol*s)'), n=2.84, Ea=(-1012, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5cC6H11 + OH <=> C2H5S4XcC6H10 + H2O""",
)

entry(
    index = 1703,
    label = "C2H5cC6H11 + O2 <=> PXC2H4cC6H11 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(50930, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5cC6H11 + O2 <=> PXC2H4cC6H11 + HO2""",
)

entry(
    index = 1704,
    label = "C2H5cC6H11 + O2 <=> SXC2H4cC6H11 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(47590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5cC6H11 + O2 <=> SXC2H4cC6H11 + HO2""",
)

entry(
    index = 1705,
    label = "C2H5cC6H11 + O2 <=> C2H5TXcC6H10 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(44000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5cC6H11 + O2 <=> C2H5TXcC6H10 + HO2""",
)

entry(
    index = 1706,
    label = "C2H5cC6H11 + O2 <=> C2H5S2XcC6H10 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+13, 'cm^3/(mol*s)'), n=0, Ea=(47590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5cC6H11 + O2 <=> C2H5S2XcC6H10 + HO2""",
)

entry(
    index = 1707,
    label = "C2H5cC6H11 + O2 <=> C2H5S3XcC6H10 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+13, 'cm^3/(mol*s)'), n=0, Ea=(47590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5cC6H11 + O2 <=> C2H5S3XcC6H10 + HO2""",
)

entry(
    index = 1708,
    label = "C2H5cC6H11 + O2 <=> C2H5S4XcC6H10 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(47590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5cC6H11 + O2 <=> C2H5S4XcC6H10 + HO2""",
)

entry(
    index = 1709,
    label = "C2H5cC6H11 + HO2 <=> PXC2H4cC6H11 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(30500, 'cm^3/(mol*s)'), n=2.65, Ea=(17496, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5cC6H11 + HO2 <=> PXC2H4cC6H11 + H2O2""",
)

entry(
    index = 1710,
    label = "C2H5cC6H11 + HO2 <=> SXC2H4cC6H11 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7130, 'cm^3/(mol*s)'), n=2.77, Ea=(14913, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5cC6H11 + HO2 <=> SXC2H4cC6H11 + H2O2""",
)

entry(
    index = 1711,
    label = "C2H5cC6H11 + HO2 <=> C2H5TXcC6H10 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1240, 'cm^3/(mol*s)'), n=2.77, Ea=(10500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5cC6H11 + HO2 <=> C2H5TXcC6H10 + H2O2""",
)

entry(
    index = 1712,
    label = "C2H5cC6H11 + HO2 <=> C2H5S2XcC6H10 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(281, 'cm^3/(mol*s)'), n=3.25, Ea=(14998, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5cC6H11 + HO2 <=> C2H5S2XcC6H10 + H2O2""",
)

entry(
    index = 1713,
    label = "C2H5cC6H11 + HO2 <=> C2H5S3XcC6H10 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(281, 'cm^3/(mol*s)'), n=3.25, Ea=(14998, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5cC6H11 + HO2 <=> C2H5S3XcC6H10 + H2O2""",
)

entry(
    index = 1714,
    label = "C2H5cC6H11 + HO2 <=> C2H5S4XcC6H10 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(140, 'cm^3/(mol*s)'), n=3.25, Ea=(14998, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5cC6H11 + HO2 <=> C2H5S4XcC6H10 + H2O2""",
)

entry(
    index = 1715,
    label = "C2H5cC6H11 + CH3 <=> SXC2H4cC6H11 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(71.3, 'cm^3/(mol*s)'), n=3.25, Ea=(13556, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5cC6H11 + CH3 <=> SXC2H4cC6H11 + CH4""",
)

entry(
    index = 1716,
    label = "C2H5cC6H11 + CH3 <=> PXC2H4cC6H11 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(25.9, 'cm^3/(mol*s)'), n=3.32, Ea=(11110, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5cC6H11 + CH3 <=> PXC2H4cC6H11 + CH4""",
)

entry(
    index = 1717,
    label = "C2H5cC6H11 + CH3 <=> C2H5TXcC6H10 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(31.3, 'cm^3/(mol*s)'), n=3.33, Ea=(8788, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5cC6H11 + CH3 <=> C2H5TXcC6H10 + CH4""",
)

entry(
    index = 1718,
    label = "C2H5cC6H11 + CH3 <=> C2H5S2XcC6H10 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(10.6, 'cm^3/(mol*s)'), n=3.47, Ea=(10333, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5cC6H11 + CH3 <=> C2H5S2XcC6H10 + CH4""",
)

entry(
    index = 1719,
    label = "C2H5cC6H11 + CH3 <=> C2H5S3XcC6H10 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(298, 'cm^3/(mol*s)'), n=3.12, Ea=(11716, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5cC6H11 + CH3 <=> C2H5S3XcC6H10 + CH4""",
)

entry(
    index = 1720,
    label = "C2H5cC6H11 + CH3 <=> C2H5S4XcC6H10 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(36, 'cm^3/(mol*s)'), n=3.25, Ea=(11369, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5cC6H11 + CH3 <=> C2H5S4XcC6H10 + CH4""",
)

entry(
    index = 1721,
    label = "PXC2H4cC6H11 <=> cC6H11 + C2H4",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(6.03e+10, 's^-1'), n=0.84, Ea=(27820, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1e-43, 'cm^3/(mol*s)'),
            n = 18.591,
            Ea = (-602.5, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -43.32,
        T3 = (200, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC2H4cC6H11 <=> cC6H11 + C2H4""",
)

entry(
    index = 1722,
    label = "SXC2H4cC6H11 <=> PX8-2C8H15",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(6.03e+12, 's^-1'), n=0.07, Ea=(26982.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.1e-33, 'cm^3/(mol*s)'),
            n = 15.29,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -25.11,
        T3 = (225, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SXC2H4cC6H11 <=> PX8-2C8H15""",
)

entry(
    index = 1723,
    label = "C2H5TXcC6H10 <=> CH2cC6H10 + CH3",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.63e+11, 's^-1'), n=0.68, Ea=(24579.2, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.1e-39, 'cm^3/(mol*s)'),
            n = 17.57,
            Ea = (-603.1, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -31.51,
        T3 = (219, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C2H5TXcC6H10 <=> CH2cC6H10 + CH3""",
)

entry(
    index = 1724,
    label = "C2H5TXcC6H10 <=> C2H5-2-PXC6H10",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(6.03e+12, 's^-1'), n=0.07, Ea=(27982.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.1e-33, 'cm^3/(mol*s)'),
            n = 15.29,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -25.11,
        T3 = (225, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C2H5TXcC6H10 <=> C2H5-2-PXC6H10""",
)

entry(
    index = 1725,
    label = "C2H5S2XcC6H10 <=> cC6H10 + C2H5",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(5.47e+11, 's^-1'), n=0.5, Ea=(27798.1, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (8.5e-37, 'cm^3/(mol*s)'),
            n = 16.21,
            Ea = (-603.1, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -26.7,
        T3 = (216, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C2H5S2XcC6H10 <=> cC6H10 + C2H5""",
)

entry(
    index = 1726,
    label = "C2H5S2XcC6H10 <=> PX8-3C8H15",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.01e+12, 's^-1'), n=0.07, Ea=(26982.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.1e-33, 'cm^3/(mol*s)'),
            n = 15.29,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -25.11,
        T3 = (225, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C2H5S2XcC6H10 <=> PX8-3C8H15""",
)

entry(
    index = 1727,
    label = "C2H5S2XcC6H10 <=> C2H5-3-PXC6H10",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.01e+12, 's^-1'), n=0.07, Ea=(27982.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.1e-33, 'cm^3/(mol*s)'),
            n = 15.29,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -25.11,
        T3 = (225, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C2H5S2XcC6H10 <=> C2H5-3-PXC6H10""",
)

entry(
    index = 1728,
    label = "C2H5S3XcC6H10 <=> S2XC8H15",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(6.03e+12, 's^-1'), n=0.07, Ea=(26982.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.1e-33, 'cm^3/(mol*s)'),
            n = 15.29,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -25.11,
        T3 = (225, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C2H5S3XcC6H10 <=> S2XC8H15""",
)

entry(
    index = 1729,
    label = "C2H5S3XcC6H10 <=> PXC2H4-4-1C6H11",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.01e+12, 's^-1'), n=0.07, Ea=(27982.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.1e-33, 'cm^3/(mol*s)'),
            n = 15.29,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -25.11,
        T3 = (225, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C2H5S3XcC6H10 <=> PXC2H4-4-1C6H11""",
)

entry(
    index = 1730,
    label = "C2H5S4XcC6H10 <=> PXCH2-5-1C7H13",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(6.03e+12, 's^-1'), n=0.07, Ea=(27982.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.1e-33, 'cm^3/(mol*s)'),
            n = 15.29,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -25.11,
        T3 = (225, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C2H5S4XcC6H10 <=> PXCH2-5-1C7H13""",
)

entry(
    index = 1731,
    label = "PXC2H4cC6H11 <=> C2H5S2XcC6H10",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.58e+08, 's^-1'), n=0.9, Ea=(22900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC2H4cC6H11 <=> C2H5S2XcC6H10""",
)

entry(
    index = 1732,
    label = "PXC2H4cC6H11 <=> C2H5S3XcC6H10",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.94e+08, 's^-1'), n=0.78, Ea=(18644, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC2H4cC6H11 <=> C2H5S3XcC6H10""",
)

entry(
    index = 1733,
    label = "PXC2H4cC6H11 <=> C2H5S4XcC6H10",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.97e+08, 's^-1'), n=0.78, Ea=(22740, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC2H4cC6H11 <=> C2H5S4XcC6H10""",
)

entry(
    index = 1734,
    label = "SXC2H4cC6H11 <=> C2H5S3XcC6H10",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.65e+08, 's^-1'), n=1.02, Ea=(28687, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC2H4cC6H11 <=> C2H5S3XcC6H10""",
)

entry(
    index = 1735,
    label = "SXC2H4cC6H11 <=> C2H5S4XcC6H10",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.19e+09, 's^-1'), n=0.92, Ea=(22700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC2H4cC6H11 <=> C2H5S4XcC6H10""",
)

entry(
    index = 1736,
    label = "PXC2H4cC6H11 + H <=> C2H5cC6H11",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.5e+48, 'cm^6/(mol^2*s)'),
            n = -9.32,
            Ea = (5833.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.498,
        T3 = (1314, 'K'),
        T1 = (1314, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC2H4cC6H11 + H <=> C2H5cC6H11""",
)

entry(
    index = 1737,
    label = "PXC2H4cC6H11 + H <=> cC6H11 + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.85e+24, 'cm^3/(mol*s)'),
        n = -2.92,
        Ea = (12505, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is PXC2H4cC6H11 + H <=> cC6H11 + C2H5""",
)

entry(
    index = 1738,
    label = "PXC2H4cC6H11 + O <=> PXCH2cC6H11 + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC2H4cC6H11 + O <=> PXCH2cC6H11 + CH2O""",
)

entry(
    index = 1739,
    label = "SXC2H4cC6H11 + H <=> C2H5cC6H11",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(4.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.7e+58, 'cm^6/(mol^2*s)'),
            n = -12.08,
            Ea = (11263.7, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.649,
        T3 = (1213.1, 'K'),
        T1 = (1213.1, 'K'),
        T2 = (13369.7, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SXC2H4cC6H11 + H <=> C2H5cC6H11""",
)

entry(
    index = 1740,
    label = "SXC2H4cC6H11 + H <=> cC6H11 + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.8e+28, 'cm^3/(mol*s)'),
        n = -3.94,
        Ea = (15916, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is SXC2H4cC6H11 + H <=> cC6H11 + C2H5""",
)

entry(
    index = 1741,
    label = "C2H5TXcC6H10 + H <=> C2H5cC6H11",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.47e+61, 'cm^6/(mol^2*s)'),
            n = -12.94,
            Ea = (8000, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        T3 = (1456.4, 'K'),
        T1 = (1000, 'K'),
        T2 = (10000.5, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C2H5TXcC6H10 + H <=> C2H5cC6H11""",
)

entry(
    index = 1742,
    label = "C2H5TXcC6H10 + H <=> cC6H11 + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.6e+36, 'cm^3/(mol*s)'),
        n = -6.12,
        Ea = (25640, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H5TXcC6H10 + H <=> cC6H11 + C2H5""",
)

entry(
    index = 1743,
    label = "C2H5S2XcC6H10 + H <=> C2H5cC6H11",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(4.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.7e+58, 'cm^6/(mol^2*s)'),
            n = -12.08,
            Ea = (11263.7, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.649,
        T3 = (1213.1, 'K'),
        T1 = (1213.1, 'K'),
        T2 = (13369.7, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C2H5S2XcC6H10 + H <=> C2H5cC6H11""",
)

entry(
    index = 1744,
    label = "C2H5S2XcC6H10 + H <=> cC6H11 + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.8e+28, 'cm^3/(mol*s)'),
        n = -3.94,
        Ea = (15916, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H5S2XcC6H10 + H <=> cC6H11 + C2H5""",
)

entry(
    index = 1745,
    label = "C2H5S3XcC6H10 + H <=> C2H5cC6H11",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(4.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.7e+58, 'cm^6/(mol^2*s)'),
            n = -12.08,
            Ea = (11263.7, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.649,
        T3 = (1213.1, 'K'),
        T1 = (1213.1, 'K'),
        T2 = (13369.7, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C2H5S3XcC6H10 + H <=> C2H5cC6H11""",
)

entry(
    index = 1746,
    label = "C2H5S3XcC6H10 + H <=> cC6H11 + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.8e+28, 'cm^3/(mol*s)'),
        n = -3.94,
        Ea = (15916, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H5S3XcC6H10 + H <=> cC6H11 + C2H5""",
)

entry(
    index = 1747,
    label = "C2H5S4XcC6H10 + H <=> C2H5cC6H11",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(4.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.7e+58, 'cm^6/(mol^2*s)'),
            n = -12.08,
            Ea = (11263.7, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.649,
        T3 = (1213.1, 'K'),
        T1 = (1213.1, 'K'),
        T2 = (13369.7, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C2H5S4XcC6H10 + H <=> C2H5cC6H11""",
)

entry(
    index = 1748,
    label = "C2H5S4XcC6H10 + H <=> cC6H11 + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.8e+28, 'cm^3/(mol*s)'),
        n = -3.94,
        Ea = (15916, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H5S4XcC6H10 + H <=> cC6H11 + C2H5""",
)

entry(
    index = 1749,
    label = "PX8-2C8H15 <=> SAX4-2C8H15",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.46e+11, 's^-1'), n=0, Ea=(10516.8, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (9.9e-38, 'cm^3/(mol*s)'),
            n = 17.215,
            Ea = (-603, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -16.33,
        T3 = (200, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PX8-2C8H15 <=> SAX4-2C8H15""",
)

entry(
    index = 1750,
    label = "PX8-2C8H15 <=> PX6-2C6H11 + C2H4",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(9.12e+11, 's^-1'), n=0.31, Ea=(27237.8, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.8e-57, 'cm^3/(mol*s)'),
            n = 23.463,
            Ea = (-602.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -2.46,
        T3 = (206, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PX8-2C8H15 <=> PX6-2C6H11 + C2H4""",
)

entry(
    index = 1751,
    label = "C2H5-2-PXC6H10 <=> PXC2H4-2-1C4H7 + C2H4",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(9.12e+11, 's^-1'), n=0.31, Ea=(27237.8, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.8e-57, 'cm^3/(mol*s)'),
            n = 23.463,
            Ea = (-602.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -2.46,
        T3 = (206, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C2H5-2-PXC6H10 <=> PXC2H4-2-1C4H7 + C2H4""",
)

entry(
    index = 1752,
    label = "PX8-3C8H15 <=> PX1-3C6H11 + C2H4",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(9.12e+11, 's^-1'), n=0.31, Ea=(27237.8, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.8e-57, 'cm^3/(mol*s)'),
            n = 23.463,
            Ea = (-602.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -2.46,
        T3 = (206, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PX8-3C8H15 <=> PX1-3C6H11 + C2H4""",
)

entry(
    index = 1753,
    label = "PX8-3C8H15 <=> SAX5-3C8H15",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.46e+11, 's^-1'), n=0, Ea=(10516.8, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (9.9e-38, 'cm^3/(mol*s)'),
            n = 17.215,
            Ea = (-603, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -16.33,
        T3 = (200, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PX8-3C8H15 <=> SAX5-3C8H15""",
)

entry(
    index = 1754,
    label = "C2H5-3-PXC6H10 <=> PXCH2-3-1C5H9 + C2H4",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(9.12e+11, 's^-1'), n=0.31, Ea=(27237.8, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.8e-57, 'cm^3/(mol*s)'),
            n = 23.463,
            Ea = (-602.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -2.46,
        T3 = (206, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C2H5-3-PXC6H10 <=> PXCH2-3-1C5H9 + C2H4""",
)

entry(
    index = 1755,
    label = "S2XC8H15 <=> C4H7 + C4H81",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(5.5e+11, 's^-1'), n=0.55, Ea=(28084.3, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.1e-43, 'cm^3/(mol*s)'),
            n = 18.418,
            Ea = (-602.9, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -32.13,
        T3 = (207, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S2XC8H15 <=> C4H7 + C4H81""",
)

entry(
    index = 1756,
    label = "S2XC8H15 <=> C7H12-16 + CH3",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(5.5e+11, 's^-1'), n=0.55, Ea=(28084.3, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.1e-43, 'cm^3/(mol*s)'),
            n = 18.418,
            Ea = (-602.9, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -32.13,
        T3 = (207, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S2XC8H15 <=> C7H12-16 + CH3""",
)

entry(
    index = 1757,
    label = "S2XC8H15 <=> SAXC8H15",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(155, 's^-1'), n=2.83, Ea=(15566.2, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.5e-30, 'cm^3/(mol*s)'),
            n = 14.56,
            Ea = (-602.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -13.59,
        T3 = (214, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S2XC8H15 <=> SAXC8H15""",
)

entry(
    index = 1758,
    label = "PXC2H4-4-1C6H11 <=> S2XC6H11 + C2H4",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(9.12e+11, 's^-1'), n=0.31, Ea=(27237.8, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.8e-57, 'cm^3/(mol*s)'),
            n = 23.463,
            Ea = (-602.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -2.46,
        T3 = (206, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC2H4-4-1C6H11 <=> S2XC6H11 + C2H4""",
)

entry(
    index = 1759,
    label = "PXC2H4-4-1C6H11 <=> C2H5-4-SAX1C6H10",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(155, 's^-1'), n=2.83, Ea=(15566.2, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.5e-30, 'cm^3/(mol*s)'),
            n = 14.56,
            Ea = (-602.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -13.59,
        T3 = (214, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC2H4-4-1C6H11 <=> C2H5-4-SAX1C6H10""",
)

entry(
    index = 1760,
    label = "PXCH2-5-1C7H13 <=> C4H7 + C4H81",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.76e+11, 's^-1'), n=0.57, Ea=(28791.6, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4.7e-39, 'cm^3/(mol*s)'),
            n = 16.77,
            Ea = (-603.1, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -27.89,
        T3 = (216, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXCH2-5-1C7H13 <=> C4H7 + C4H81""",
)

entry(
    index = 1761,
    label = "PXCH2-5-1C7H13 <=> C6H10-15 + C2H5",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.76e+11, 's^-1'), n=0.57, Ea=(28791.6, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4.7e-39, 'cm^3/(mol*s)'),
            n = 16.77,
            Ea = (-603.1, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -27.89,
        T3 = (216, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXCH2-5-1C7H13 <=> C6H10-15 + C2H5""",
)

entry(
    index = 1762,
    label = "PXCH2-5-1C7H13 <=> CH3-5-SAX1C7H12",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(155, 's^-1'), n=2.83, Ea=(15566.2, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.5e-30, 'cm^3/(mol*s)'),
            n = 14.56,
            Ea = (-602.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -13.59,
        T3 = (214, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXCH2-5-1C7H13 <=> CH3-5-SAX1C7H12""",
)

entry(
    index = 1763,
    label = "SAX4-2C8H15 <=> C5H8-13 + nC3H7",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.39e+11, 's^-1'), n=0.66, Ea=(32262.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4e-42, 'cm^3/(mol*s)'),
            n = 18.045,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -18.5,
        T3 = (246, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SAX4-2C8H15 <=> C5H8-13 + nC3H7""",
)

entry(
    index = 1764,
    label = "C2H5-2-SAX1C6H10 <=> CH2-3-1C5H8 + C2H5",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.39e+11, 's^-1'), n=0.66, Ea=(32262.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4e-42, 'cm^3/(mol*s)'),
            n = 18.045,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -18.5,
        T3 = (246, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C2H5-2-SAX1C6H10 <=> CH2-3-1C5H8 + C2H5""",
)

entry(
    index = 1765,
    label = "C6H10-12 + C2H5 <=> C2H5-2-SAX1C6H10",
    degeneracy = 1,
    duplicate = True,
    kinetics = Arrhenius(A=(2e+11, 'cm^3/(mol*s)'), n=0, Ea=(7500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H10-12 + C2H5 <=> C2H5-2-SAX1C6H10""",
)

entry(
    index = 1766,
    label = "SAX5-3C8H15 <=> C6H10-13 + C2H5",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.39e+11, 's^-1'), n=0.66, Ea=(32262.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4e-42, 'cm^3/(mol*s)'),
            n = 18.045,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -18.5,
        T3 = (246, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SAX5-3C8H15 <=> C6H10-13 + C2H5""",
)

entry(
    index = 1767,
    label = "C2H5-3-TAX1C6H10 <=> CH2-3-1C6H10 + CH3",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.39e+11, 's^-1'), n=0.66, Ea=(32262.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4e-42, 'cm^3/(mol*s)'),
            n = 18.045,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -18.5,
        T3 = (246, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C2H5-3-TAX1C6H10 <=> CH2-3-1C6H10 + CH3""",
)

entry(
    index = 1768,
    label = "C2H5-3-TAX1C6H10 <=> CH2-3-1C5H8 + C2H5",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.39e+11, 's^-1'), n=0.66, Ea=(32262.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4e-42, 'cm^3/(mol*s)'),
            n = 18.045,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -18.5,
        T3 = (246, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C2H5-3-TAX1C6H10 <=> CH2-3-1C5H8 + C2H5""",
)

entry(
    index = 1769,
    label = "SAXC8H15 <=> C4H6 + pC4H9",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.39e+11, 's^-1'), n=0.66, Ea=(32262.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4e-42, 'cm^3/(mol*s)'),
            n = 18.045,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -18.5,
        T3 = (246, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SAXC8H15 <=> C4H6 + pC4H9""",
)

entry(
    index = 1770,
    label = "C2H5-4-SAX1C6H10 <=> C6H10-13 + C2H5",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.39e+11, 's^-1'), n=0.66, Ea=(32262.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4e-42, 'cm^3/(mol*s)'),
            n = 18.045,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -18.5,
        T3 = (246, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C2H5-4-SAX1C6H10 <=> C6H10-13 + C2H5""",
)

entry(
    index = 1771,
    label = "CH3-5-SAX1C7H12 <=> C4H6 + sC4H9",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.39e+11, 's^-1'), n=0.66, Ea=(31262.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4e-42, 'cm^3/(mol*s)'),
            n = 18.045,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -18.5,
        T3 = (246, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is CH3-5-SAX1C7H12 <=> C4H6 + sC4H9""",
)

entry(
    index = 1772,
    label = "PX1-4C8H15 <=> C2H4 + SAXC6H11",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.04e+12, 's^-1'), n=-0.37, Ea=(25124.2, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.9e-31, 'cm^3/(mol*s)'),
            n = 13.982,
            Ea = (-603.1, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -14.78,
        T3 = (229, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PX1-4C8H15 <=> C2H4 + SAXC6H11""",
)

entry(
    index = 1773,
    label = "PXC2H4-2-1C6H11 <=> C2H5-2-SAX1C6H10",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(155, 's^-1'), n=2.83, Ea=(15566.2, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.5e-30, 'cm^3/(mol*s)'),
            n = 14.56,
            Ea = (-602.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -13.59,
        T3 = (214, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC2H4-2-1C6H11 <=> C2H5-2-SAX1C6H10""",
)

entry(
    index = 1774,
    label = "C2H5-2-SAX1C6H10 <=> C6H10-12 + C2H5",
    degeneracy = 1,
    duplicate = True,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.39e+11, 's^-1'), n=0.66, Ea=(32262.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4e-42, 'cm^3/(mol*s)'),
            n = 18.045,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -18.5,
        T3 = (246, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C2H5-2-SAX1C6H10 <=> C6H10-12 + C2H5""",
)

entry(
    index = 1775,
    label = "PX1-3C8H15 <=> C8H14-13 + H",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.05e+08, 's^-1'), n=-1.35, Ea=(32487.5, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (2.2e-38, 'cm^3/(mol*s)'),
            n = 13.17,
            Ea = (-603.1, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -58.42,
        T3 = (256, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PX1-3C8H15 <=> C8H14-13 + H""",
)

entry(
    index = 1776,
    label = "PXCH2-3-1C7H13 <=> C4H6 + pC4H9",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(4.6e+11, 's^-1'), n=-0.41, Ea=(31253.7, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.8e-28, 'cm^3/(mol*s)'),
            n = 12.34,
            Ea = (-603.2, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -19.24,
        T3 = (231, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXCH2-3-1C7H13 <=> C4H6 + pC4H9""",
)

entry(
    index = 1777,
    label = "S4XC8H15 <=> C5H8-14 + nC3H7",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.05e+13, 's^-1'), n=0.17, Ea=(13408.3, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (9.1e-19, 'cm^3/(mol*s)'),
            n = 11.346,
            Ea = (-602.9, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -56.53,
        T3 = (201, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S4XC8H15 <=> C5H8-14 + nC3H7""",
)

entry(
    index = 1778,
    label = "C2H3cC6H11 + H <=> C2H4 + cC6H11",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8e+21, 'cm^3/(mol*s)'),
        n = -2.39,
        Ea = (11180, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H3cC6H11 + H <=> C2H4 + cC6H11""",
)

entry(
    index = 1779,
    label = "C2H3cC6H11 + O <=> PXCH2cC6H11 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.3e+08, 'cm^3/(mol*s)'),
        n = 1.45,
        Ea = (-402, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H3cC6H11 + O <=> PXCH2cC6H11 + HCO""",
)

entry(
    index = 1780,
    label = "cC6H11 + C2H3 <=> C2H3cC6H11",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.25e+14, 'cm^3/(mol*s)'), n=-0.7, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC6H11 + C2H3 <=> C2H3cC6H11""",
)

entry(
    index = 1781,
    label = "C8H14-13 + H <=> C6H12 + C2H3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.46e+30, 'cm^3/(mol*s)'),
        n = -4.34,
        Ea = (21647, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C8H14-13 + H <=> C6H12 + C2H3""",
)

entry(
    index = 1782,
    label = "CH2-3-1C7H12 + H <=> C6H12 + C2H3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.46e+30, 'cm^3/(mol*s)'),
        n = -4.34,
        Ea = (21647, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2-3-1C7H12 + H <=> C6H12 + C2H3""",
)

entry(
    index = 1783,
    label = "PXCH2-2-C6H13 <=> C6H12 + CH3",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.76e+11, 's^-1'), n=0.57, Ea=(28791.6, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4.7e-39, 'cm^3/(mol*s)'),
            n = 16.77,
            Ea = (-603.1, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -27.89,
        T3 = (216, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXCH2-2-C6H13 <=> C6H12 + CH3""",
)

entry(
    index = 1784,
    label = "PXCH2-2-C6H13 <=> C3H6 + pC4H9",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(5.47e+11, 's^-1'), n=0.5, Ea=(27798.1, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (8.5e-37, 'cm^3/(mol*s)'),
            n = 16.21,
            Ea = (-603.1, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -26.7,
        T3 = (216, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXCH2-2-C6H13 <=> C3H6 + pC4H9""",
)

entry(
    index = 1785,
    label = "PXCH2-2-C6H13 <=> CH3-2-SXC6H12",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(183, 's^-1'), n=2.55, Ea=(10960.3, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.4e-26, 'cm^3/(mol*s)'),
            n = 13.087,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -25.38,
        T3 = (215, 'K'),
        T1 = (28, 'K'),
        T2 = (5e+06, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXCH2-2-C6H13 <=> CH3-2-SXC6H12""",
)

entry(
    index = 1786,
    label = "CH3-2-SXC6H12 <=> iC4H9 + C3H6",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(4.47e+11, 's^-1'), n=0.57, Ea=(28044.5, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.7e-33, 'cm^3/(mol*s)'),
            n = 14.91,
            Ea = (-600, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -6.53,
        T3 = (333, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is CH3-2-SXC6H12 <=> iC4H9 + C3H6""",
)

entry(
    index = 1787,
    label = "cC6H11 + CH3 <=> CH3cC6H11",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.4e+15, 'cm^3/(mol*s)'), n=-0.68, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC6H11 + CH3 <=> CH3cC6H11""",
)

entry(
    index = 1788,
    label = "CH3cC6H11 <=> C7H14-2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.17e+15, 's^-1'), n=0, Ea=(71000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3cC6H11 <=> C7H14-2""",
)

entry(
    index = 1789,
    label = "CH3cC6H11 <=> C7H14",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.17e+15, 's^-1'), n=0, Ea=(71000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3cC6H11 <=> C7H14""",
)

entry(
    index = 1790,
    label = "C7H14-2 <=> SAXC4H7 + nC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.07e+23, 's^-1'), n=-2.03, Ea=(74958, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-2 <=> SAXC4H7 + nC3H7""",
)

entry(
    index = 1791,
    label = "C7H14-2 <=> C4H81 + C3H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.42e+07, 's^-1'), n=1.65, Ea=(53752, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-2 <=> C4H81 + C3H6""",
)

entry(
    index = 1792,
    label = "C7H14-2 + H <=> C3H6 + pC4H9",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8e+21, 'cm^3/(mol*s)'),
        n = -2.39,
        Ea = (11180, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C7H14-2 + H <=> C3H6 + pC4H9""",
)

entry(
    index = 1793,
    label = "C7H14-2 + H <=> C4H81 + nC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8e+21, 'cm^3/(mol*s)'),
        n = -2.39,
        Ea = (11180, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C7H14-2 + H <=> C4H81 + nC3H7""",
)

entry(
    index = 1794,
    label = "C7H14-2 + O <=> C2H3CHO + pC4H9 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+08, 'cm^3/(mol*s)'), n=1.65, Ea=(327, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-2 + O <=> C2H3CHO + pC4H9 + H""",
)

entry(
    index = 1795,
    label = "C7H14-2 + OH <=> aC3H4 + H2O + pC4H9",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.1e+06, 'cm^3/(mol*s)'), n=2, Ea=(-298, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-2 + OH <=> aC3H4 + H2O + pC4H9""",
)

entry(
    index = 1796,
    label = "CH3cC6H11 + H <=> PXCH2cC6H11 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(168000, 'cm^3/(mol*s)'), n=2.61, Ea=(6535, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3cC6H11 + H <=> PXCH2cC6H11 + H2""",
)

entry(
    index = 1797,
    label = "CH3cC6H11 + H <=> CH3TXcC6H10 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(543000, 'cm^3/(mol*s)'), n=2.36, Ea=(2887, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3cC6H11 + H <=> CH3TXcC6H10 + H2""",
)

entry(
    index = 1798,
    label = "CH3cC6H11 + H <=> CH3S2XcC6H10 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(311000, 'cm^3/(mol*s)'), n=2.57, Ea=(4359, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3cC6H11 + H <=> CH3S2XcC6H10 + H2""",
)

entry(
    index = 1799,
    label = "CH3cC6H11 + H <=> CH3S3XcC6H10 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(297000, 'cm^3/(mol*s)'), n=2.59, Ea=(4138, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3cC6H11 + H <=> CH3S3XcC6H10 + H2""",
)

entry(
    index = 1800,
    label = "CH3cC6H11 + H <=> CH3S4XcC6H10 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(266000, 'cm^3/(mol*s)'), n=2.52, Ea=(4452, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3cC6H11 + H <=> CH3S4XcC6H10 + H2""",
)

entry(
    index = 1801,
    label = "CH3cC6H11 + O <=> PXCH2cC6H11 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.2e+06, 'cm^3/(mol*s)'), n=2.4, Ea=(5504, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3cC6H11 + O <=> PXCH2cC6H11 + OH""",
)

entry(
    index = 1802,
    label = "CH3cC6H11 + O <=> CH3TXcC6H10 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(157000, 'cm^3/(mol*s)'), n=2.5, Ea=(1110, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3cC6H11 + O <=> CH3TXcC6H10 + OH""",
)

entry(
    index = 1803,
    label = "CH3cC6H11 + O <=> CH3S2XcC6H10 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(95200, 'cm^3/(mol*s)'), n=2.71, Ea=(2106, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3cC6H11 + O <=> CH3S2XcC6H10 + OH""",
)

entry(
    index = 1804,
    label = "CH3cC6H11 + O <=> CH3S3XcC6H10 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(95200, 'cm^3/(mol*s)'), n=2.71, Ea=(2106, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3cC6H11 + O <=> CH3S3XcC6H10 + OH""",
)

entry(
    index = 1805,
    label = "CH3cC6H11 + O <=> CH3S4XcC6H10 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(47600, 'cm^3/(mol*s)'), n=2.71, Ea=(2106, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3cC6H11 + O <=> CH3S4XcC6H10 + OH""",
)

entry(
    index = 1806,
    label = "CH3cC6H11 + OH <=> PXCH2cC6H11 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5330, 'cm^3/(mol*s)'), n=2.9, Ea=(504, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3cC6H11 + OH <=> PXCH2cC6H11 + H2O""",
)

entry(
    index = 1807,
    label = "CH3cC6H11 + OH <=> CH3TXcC6H10 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6690, 'cm^3/(mol*s)'), n=2.82, Ea=(2966, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3cC6H11 + OH <=> CH3TXcC6H10 + H2O""",
)

entry(
    index = 1808,
    label = "CH3cC6H11 + OH <=> CH3S2XcC6H10 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8130, 'cm^3/(mol*s)'), n=2.84, Ea=(-1479, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3cC6H11 + OH <=> CH3S2XcC6H10 + H2O""",
)

entry(
    index = 1809,
    label = "CH3cC6H11 + OH <=> CH3S3XcC6H10 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9500, 'cm^3/(mol*s)'), n=2.85, Ea=(-1023, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3cC6H11 + OH <=> CH3S3XcC6H10 + H2O""",
)

entry(
    index = 1810,
    label = "CH3cC6H11 + OH <=> CH3S4XcC6H10 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5820, 'cm^3/(mol*s)'), n=2.85, Ea=(-1011, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3cC6H11 + OH <=> CH3S4XcC6H10 + H2O""",
)

entry(
    index = 1811,
    label = "CH3cC6H11 + O2 <=> PXCH2cC6H11 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(50930, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3cC6H11 + O2 <=> PXCH2cC6H11 + HO2""",
)

entry(
    index = 1812,
    label = "CH3cC6H11 + O2 <=> CH3TXcC6H10 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(44000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3cC6H11 + O2 <=> CH3TXcC6H10 + HO2""",
)

entry(
    index = 1813,
    label = "CH3cC6H11 + O2 <=> CH3S2XcC6H10 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+13, 'cm^3/(mol*s)'), n=0, Ea=(47590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3cC6H11 + O2 <=> CH3S2XcC6H10 + HO2""",
)

entry(
    index = 1814,
    label = "CH3cC6H11 + O2 <=> CH3S3XcC6H10 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+13, 'cm^3/(mol*s)'), n=0, Ea=(47590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3cC6H11 + O2 <=> CH3S3XcC6H10 + HO2""",
)

entry(
    index = 1815,
    label = "CH3cC6H11 + O2 <=> CH3S4XcC6H10 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(47590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3cC6H11 + O2 <=> CH3S4XcC6H10 + HO2""",
)

entry(
    index = 1816,
    label = "CH3cC6H11 + HO2 <=> PXCH2cC6H11 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(30500, 'cm^3/(mol*s)'), n=2.65, Ea=(17496, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3cC6H11 + HO2 <=> PXCH2cC6H11 + H2O2""",
)

entry(
    index = 1817,
    label = "CH3cC6H11 + HO2 <=> CH3TXcC6H10 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1240, 'cm^3/(mol*s)'), n=2.77, Ea=(10500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3cC6H11 + HO2 <=> CH3TXcC6H10 + H2O2""",
)

entry(
    index = 1818,
    label = "CH3cC6H11 + HO2 <=> CH3S2XcC6H10 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(281, 'cm^3/(mol*s)'), n=3.25, Ea=(14998, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3cC6H11 + HO2 <=> CH3S2XcC6H10 + H2O2""",
)

entry(
    index = 1819,
    label = "CH3cC6H11 + HO2 <=> CH3S3XcC6H10 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(281, 'cm^3/(mol*s)'), n=3.25, Ea=(14998, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3cC6H11 + HO2 <=> CH3S3XcC6H10 + H2O2""",
)

entry(
    index = 1820,
    label = "CH3cC6H11 + HO2 <=> CH3S4XcC6H10 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(140, 'cm^3/(mol*s)'), n=3.25, Ea=(14998, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3cC6H11 + HO2 <=> CH3S4XcC6H10 + H2O2""",
)

entry(
    index = 1821,
    label = "CH3cC6H11 + CH3 <=> PXCH2cC6H11 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(37.5, 'cm^3/(mol*s)'), n=3.27, Ea=(13516, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3cC6H11 + CH3 <=> PXCH2cC6H11 + CH4""",
)

entry(
    index = 1822,
    label = "CH3cC6H11 + CH3 <=> CH3TXcC6H10 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(19.1, 'cm^3/(mol*s)'), n=3.27, Ea=(9022, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3cC6H11 + CH3 <=> CH3TXcC6H10 + CH4""",
)

entry(
    index = 1823,
    label = "CH3cC6H11 + CH3 <=> CH3S2XcC6H10 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(66.8, 'cm^3/(mol*s)'), n=3.21, Ea=(11418, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3cC6H11 + CH3 <=> CH3S2XcC6H10 + CH4""",
)

entry(
    index = 1824,
    label = "CH3cC6H11 + CH3 <=> CH3S3XcC6H10 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(71.8, 'cm^3/(mol*s)'), n=3.26, Ea=(11303, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3cC6H11 + CH3 <=> CH3S3XcC6H10 + CH4""",
)

entry(
    index = 1825,
    label = "CH3cC6H11 + CH3 <=> CH3S4XcC6H10 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(18.2, 'cm^3/(mol*s)'), n=3.36, Ea=(10931, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3cC6H11 + CH3 <=> CH3S4XcC6H10 + CH4""",
)

entry(
    index = 1826,
    label = "PXCH2cC6H11 <=> PXC7H13",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(5.47e+11, 's^-1'), n=0.5, Ea=(27798.1, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (8.5e-37, 'cm^3/(mol*s)'),
            n = 16.21,
            Ea = (-603.1, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -26.7,
        T3 = (216, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXCH2cC6H11 <=> PXC7H13""",
)

entry(
    index = 1827,
    label = "CH3TXcC6H10 <=> CH3-2-PXC6H10",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(6.03e+12, 's^-1'), n=0.07, Ea=(27982.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.1e-33, 'cm^3/(mol*s)'),
            n = 15.29,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -25.11,
        T3 = (225, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is CH3TXcC6H10 <=> CH3-2-PXC6H10""",
)

entry(
    index = 1828,
    label = "CH3S2XcC6H10 <=> PX7-2C7H13",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.01e+12, 's^-1'), n=0.07, Ea=(26982.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.1e-33, 'cm^3/(mol*s)'),
            n = 15.29,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -25.11,
        T3 = (225, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is CH3S2XcC6H10 <=> PX7-2C7H13""",
)

entry(
    index = 1829,
    label = "CH3S2XcC6H10 <=> cC6H10 + CH3",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.76e+11, 's^-1'), n=0.57, Ea=(28791.6, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4.7e-39, 'cm^3/(mol*s)'),
            n = 16.77,
            Ea = (-603.1, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -27.89,
        T3 = (216, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is CH3S2XcC6H10 <=> cC6H10 + CH3""",
)

entry(
    index = 1830,
    label = "CH3S2XcC6H10 <=> CH3-3-PXC6H10",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.01e+12, 's^-1'), n=0.07, Ea=(27982.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.1e-33, 'cm^3/(mol*s)'),
            n = 15.29,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -25.11,
        T3 = (225, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is CH3S2XcC6H10 <=> CH3-3-PXC6H10""",
)

entry(
    index = 1831,
    label = "CH3S3XcC6H10 <=> SXC7H13",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.01e+12, 's^-1'), n=0.07, Ea=(26982.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.1e-33, 'cm^3/(mol*s)'),
            n = 15.29,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -25.11,
        T3 = (225, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is CH3S3XcC6H10 <=> SXC7H13""",
)

entry(
    index = 1832,
    label = "CH3S3XcC6H10 <=> CH3-4-PXC6H10",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.01e+12, 's^-1'), n=0.07, Ea=(27982.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.1e-33, 'cm^3/(mol*s)'),
            n = 15.29,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -25.11,
        T3 = (225, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is CH3S3XcC6H10 <=> CH3-4-PXC6H10""",
)

entry(
    index = 1833,
    label = "CH3S4XcC6H10 <=> PXCH2-5-1C6H11",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(6.03e+12, 's^-1'), n=0.07, Ea=(27982.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.1e-33, 'cm^3/(mol*s)'),
            n = 15.29,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -25.11,
        T3 = (225, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is CH3S4XcC6H10 <=> PXCH2-5-1C6H11""",
)

entry(
    index = 1834,
    label = "PXCH2cC6H11 <=> CH3S3XcC6H10",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.65e+08, 's^-1'), n=1.02, Ea=(28687, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXCH2cC6H11 <=> CH3S3XcC6H10""",
)

entry(
    index = 1835,
    label = "PXCH2cC6H11 <=> CH3S4XcC6H10",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.19e+09, 's^-1'), n=0.92, Ea=(22700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXCH2cC6H11 <=> CH3S4XcC6H10""",
)

entry(
    index = 1836,
    label = "PXCH2cC6H11 + H <=> CH3cC6H11",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.5e+48, 'cm^6/(mol^2*s)'),
            n = -9.32,
            Ea = (5833.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.498,
        T3 = (1314, 'K'),
        T1 = (1314, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXCH2cC6H11 + H <=> CH3cC6H11""",
)

entry(
    index = 1837,
    label = "PXCH2cC6H11 + H <=> cC6H11 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.85e+24, 'cm^3/(mol*s)'),
        n = -2.92,
        Ea = (12505, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is PXCH2cC6H11 + H <=> cC6H11 + CH3""",
)

entry(
    index = 1838,
    label = "PXCH2cC6H11 + O <=> cC6H11 + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXCH2cC6H11 + O <=> cC6H11 + CH2O""",
)

entry(
    index = 1839,
    label = "CH3TXcC6H10 + H <=> CH3cC6H11",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.47e+61, 'cm^6/(mol^2*s)'),
            n = -12.94,
            Ea = (8000, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        T3 = (1456.4, 'K'),
        T1 = (1000, 'K'),
        T2 = (10000.5, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is CH3TXcC6H10 + H <=> CH3cC6H11""",
)

entry(
    index = 1840,
    label = "CH3TXcC6H10 + H <=> cC6H11 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.6e+36, 'cm^3/(mol*s)'),
        n = -6.12,
        Ea = (25640, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3TXcC6H10 + H <=> cC6H11 + CH3""",
)

entry(
    index = 1841,
    label = "CH3S2XcC6H10 + H <=> CH3cC6H11",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(4.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.7e+58, 'cm^6/(mol^2*s)'),
            n = -12.08,
            Ea = (11263.7, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.649,
        T3 = (1213.1, 'K'),
        T1 = (1213.1, 'K'),
        T2 = (13369.7, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is CH3S2XcC6H10 + H <=> CH3cC6H11""",
)

entry(
    index = 1842,
    label = "CH3S2XcC6H10 + H <=> cC6H11 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.8e+28, 'cm^3/(mol*s)'),
        n = -3.94,
        Ea = (15916, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3S2XcC6H10 + H <=> cC6H11 + CH3""",
)

entry(
    index = 1843,
    label = "CH3S3XcC6H10 + H <=> CH3cC6H11",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(4.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.7e+58, 'cm^6/(mol^2*s)'),
            n = -12.08,
            Ea = (11263.7, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.649,
        T3 = (1213.1, 'K'),
        T1 = (1213.1, 'K'),
        T2 = (13369.7, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is CH3S3XcC6H10 + H <=> CH3cC6H11""",
)

entry(
    index = 1844,
    label = "CH3S3XcC6H10 + H <=> cC6H11 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.8e+28, 'cm^3/(mol*s)'),
        n = -3.94,
        Ea = (15916, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3S3XcC6H10 + H <=> cC6H11 + CH3""",
)

entry(
    index = 1845,
    label = "CH3S4XcC6H10 + H <=> CH3cC6H11",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(4.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.7e+58, 'cm^6/(mol^2*s)'),
            n = -12.08,
            Ea = (11263.7, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.649,
        T3 = (1213.1, 'K'),
        T1 = (1213.1, 'K'),
        T2 = (13369.7, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is CH3S4XcC6H10 + H <=> CH3cC6H11""",
)

entry(
    index = 1846,
    label = "CH3S4XcC6H10 + H <=> cC6H11 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.8e+28, 'cm^3/(mol*s)'),
        n = -3.94,
        Ea = (15916, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3S4XcC6H10 + H <=> cC6H11 + CH3""",
)

entry(
    index = 1847,
    label = "SXC7H13 <=> SAXC7H13",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(155, 's^-1'), n=2.83, Ea=(15566.2, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.5e-30, 'cm^3/(mol*s)'),
            n = 14.56,
            Ea = (-602.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -13.59,
        T3 = (214, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SXC7H13 <=> SAXC7H13""",
)

entry(
    index = 1848,
    label = "SXC7H13 <=> C4H7 + C3H6",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(5.5e+11, 's^-1'), n=0.55, Ea=(28084.3, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.1e-43, 'cm^3/(mol*s)'),
            n = 18.418,
            Ea = (-602.9, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -32.13,
        T3 = (207, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SXC7H13 <=> C4H7 + C3H6""",
)

entry(
    index = 1849,
    label = "CH3-4-PXC6H10 <=> SXC5H9 + C2H4",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.76e+11, 's^-1'), n=0.57, Ea=(28791, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (2.8e-44, 'cm^3/(mol*s)'),
            n = 18.729,
            Ea = (-602.5, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -14.66,
        T3 = (219, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is CH3-4-PXC6H10 <=> SXC5H9 + C2H4""",
)

entry(
    index = 1850,
    label = "CH3-4-PXC6H10 <=> CH3-4-SAXC6H10",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(155, 's^-1'), n=2.83, Ea=(15566.2, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.5e-30, 'cm^3/(mol*s)'),
            n = 14.56,
            Ea = (-602.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -13.59,
        T3 = (214, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is CH3-4-PXC6H10 <=> CH3-4-SAXC6H10""",
)

entry(
    index = 1851,
    label = "PXCH2-5-1C6H11 <=> C4H7 + C3H6",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.76e+11, 's^-1'), n=0.57, Ea=(28791.6, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4.7e-39, 'cm^3/(mol*s)'),
            n = 16.77,
            Ea = (-603.1, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -27.89,
        T3 = (216, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXCH2-5-1C6H11 <=> C4H7 + C3H6""",
)

entry(
    index = 1852,
    label = "PXCH2-5-1C6H11 <=> C6H10-15 + CH3",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.76e+11, 's^-1'), n=0.57, Ea=(28791.6, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4.7e-39, 'cm^3/(mol*s)'),
            n = 16.77,
            Ea = (-603.1, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -27.89,
        T3 = (216, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXCH2-5-1C6H11 <=> C6H10-15 + CH3""",
)

entry(
    index = 1853,
    label = "PXCH2-5-1C6H11 <=> CH3-5-SAXC6H10",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(155, 's^-1'), n=2.83, Ea=(15566.2, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.5e-30, 'cm^3/(mol*s)'),
            n = 14.56,
            Ea = (-602.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -13.59,
        T3 = (214, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXCH2-5-1C6H11 <=> CH3-5-SAXC6H10""",
)

entry(
    index = 1854,
    label = "PX7-2C7H13 <=> SAXC5H9 + C2H4",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(9.12e+11, 's^-1'), n=0.31, Ea=(27237.8, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.8e-57, 'cm^3/(mol*s)'),
            n = 23.463,
            Ea = (-602.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -2.46,
        T3 = (206, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PX7-2C7H13 <=> SAXC5H9 + C2H4""",
)

entry(
    index = 1855,
    label = "PX7-2C7H13 <=> SAX4-2C7H13",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(155, 's^-1'), n=2.83, Ea=(15566.2, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.5e-30, 'cm^3/(mol*s)'),
            n = 14.56,
            Ea = (-602.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -13.59,
        T3 = (214, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PX7-2C7H13 <=> SAX4-2C7H13""",
)

entry(
    index = 1856,
    label = "SAX4-2C7H13 <=> C5H8-13 + C2H5",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.05e+08, 's^-1'), n=-1.35, Ea=(32487.5, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (2.2e-38, 'cm^3/(mol*s)'),
            n = 13.17,
            Ea = (-603.1, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -58.42,
        T3 = (256, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SAX4-2C7H13 <=> C5H8-13 + C2H5""",
)

entry(
    index = 1857,
    label = "CH3-2-PXC6H10 <=> CH3-2-PXC4H6 + C2H4",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.98e+12, 's^-1'), n=0.12, Ea=(27571.6, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.3e-43, 'cm^3/(mol*s)'),
            n = 18.35,
            Ea = (-602.5, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -13.87,
        T3 = (227, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is CH3-2-PXC6H10 <=> CH3-2-PXC4H6 + C2H4""",
)

entry(
    index = 1858,
    label = "CH3-2-PXC6H10 + H <=> CH3-2-1C6H11",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.01e+48, 'cm^6/(mol^2*s)'),
            n = -9.32,
            Ea = (5833.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.498,
        T3 = (1314, 'K'),
        T1 = (1314, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is CH3-2-PXC6H10 + H <=> CH3-2-1C6H11""",
)

entry(
    index = 1859,
    label = "CH3-2-PXC6H10 <=> CH3-2-SAXC6H10",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(155, 's^-1'), n=2.83, Ea=(15566.2, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.5e-30, 'cm^3/(mol*s)'),
            n = 14.56,
            Ea = (-602.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -13.59,
        T3 = (214, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is CH3-2-PXC6H10 <=> CH3-2-SAXC6H10""",
)

entry(
    index = 1860,
    label = "CH3-2-SAXC6H10 <=> CH3-2-C4H5-13 + C2H5",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.05e+08, 's^-1'), n=-1.35, Ea=(32487.5, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (2.2e-38, 'cm^3/(mol*s)'),
            n = 13.17,
            Ea = (-603.1, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -58.42,
        T3 = (256, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is CH3-2-SAXC6H10 <=> CH3-2-C4H5-13 + C2H5""",
)

entry(
    index = 1861,
    label = "CH3-3-PXC6H10 <=> PXCH2-3-1C4H7 + C2H4",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2.38e+11, 's^-1'), n=0.51, Ea=(27281.5, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.8e-57, 'cm^3/(mol*s)'),
            n = 23.463,
            Ea = (-602.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -2.46,
        T3 = (206, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is CH3-3-PXC6H10 <=> PXCH2-3-1C4H7 + C2H4""",
)

entry(
    index = 1862,
    label = "CH3-3-PXC6H10 <=> CH3-3-TAXC6H10",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(86, 's^-1'), n=2.62, Ea=(8722.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (8.1e-33, 'cm^3/(mol*s)'),
            n = 15.214,
            Ea = (-677.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -30.39,
        T3 = (206, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is CH3-3-PXC6H10 <=> CH3-3-TAXC6H10""",
)

entry(
    index = 1863,
    label = "CH3-3-TAXC6H10 <=> CH3-2-C4H5-13 + C2H5",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(4.6e+11, 's^-1'), n=-0.41, Ea=(31253.7, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.8e-28, 'cm^3/(mol*s)'),
            n = 12.34,
            Ea = (-603.2, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -19.24,
        T3 = (231, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is CH3-3-TAXC6H10 <=> CH3-2-C4H5-13 + C2H5""",
)

entry(
    index = 1864,
    label = "CH3-4-SAXC6H10 <=> C5H8-13 + C2H5",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.39e+11, 's^-1'), n=0.66, Ea=(31262.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4e-42, 'cm^3/(mol*s)'),
            n = 18.045,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -18.5,
        T3 = (246, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is CH3-4-SAXC6H10 <=> C5H8-13 + C2H5""",
)

entry(
    index = 1865,
    label = "CH3-4-SAXC6H10 <=> C6H10-13 + CH3",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(4.6e+11, 's^-1'), n=-0.41, Ea=(31253.7, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.8e-28, 'cm^3/(mol*s)'),
            n = 12.34,
            Ea = (-603.2, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -19.24,
        T3 = (231, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is CH3-4-SAXC6H10 <=> C6H10-13 + CH3""",
)

entry(
    index = 1866,
    label = "SAXC7H13 <=> C4H6 + nC3H7",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.39e+11, 's^-1'), n=0.66, Ea=(32262.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4e-42, 'cm^3/(mol*s)'),
            n = 18.045,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -18.5,
        T3 = (246, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SAXC7H13 <=> C4H6 + nC3H7""",
)

entry(
    index = 1867,
    label = "CH3-5-SAXC6H10 <=> C4H6 + iC3H7",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.39e+11, 's^-1'), n=0.66, Ea=(31262.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4e-42, 'cm^3/(mol*s)'),
            n = 18.045,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -18.5,
        T3 = (246, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is CH3-5-SAXC6H10 <=> C4H6 + iC3H7""",
)

entry(
    index = 1868,
    label = "PAXCH2-2-1C6H11 <=> aC3H4 + pC4H9",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 's^-1'), n=0, Ea=(50078, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PAXCH2-2-1C6H11 <=> aC3H4 + pC4H9""",
)

entry(
    index = 1869,
    label = "PX1-3C7H13 <=> C7H12-13 + H",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(4.27e+09, 's^-1'), n=-0.96, Ea=(31962.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3e-45, 'cm^3/(mol*s)'),
            n = 13.55,
            Ea = (-603.1, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -35.89,
        T3 = (249, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PX1-3C7H13 <=> C7H12-13 + H""",
)

entry(
    index = 1870,
    label = "PXCH2-3-1C6H11 <=> nC3H7 + C4H6",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(4.6e+11, 's^-1'), n=-0.41, Ea=(31253.7, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.8e-28, 'cm^3/(mol*s)'),
            n = 12.34,
            Ea = (-603.2, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -19.24,
        T3 = (231, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXCH2-3-1C6H11 <=> nC3H7 + C4H6""",
)

entry(
    index = 1871,
    label = "PXC2H4-2-1C5H9 <=> C2H5-2-SAX1C5H9",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(155, 's^-1'), n=2.83, Ea=(15566.2, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.5e-30, 'cm^3/(mol*s)'),
            n = 14.56,
            Ea = (-602.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -13.59,
        T3 = (214, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC2H4-2-1C5H9 <=> C2H5-2-SAX1C5H9""",
)

entry(
    index = 1872,
    label = "C5H8-12 + C2H5 <=> C2H5-2-SAX1C5H9",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+11, 'cm^3/(mol*s)'), n=0, Ea=(7500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H8-12 + C2H5 <=> C2H5-2-SAX1C5H9""",
)

entry(
    index = 1873,
    label = "C2H5-2-C4H513 + CH3 <=> C2H5-2-SAX1C5H9",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+11, 'cm^3/(mol*s)'), n=0, Ea=(7500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5-2-C4H513 + CH3 <=> C2H5-2-SAX1C5H9""",
)

entry(
    index = 1874,
    label = "S3XC7H13 <=> C5H8-14 + C2H5",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(5.5e+11, 's^-1'), n=0.55, Ea=(28084.3, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.1e-43, 'cm^3/(mol*s)'),
            n = 18.418,
            Ea = (-602.9, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -32.13,
        T3 = (207, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S3XC7H13 <=> C5H8-14 + C2H5""",
)

entry(
    index = 1875,
    label = "S3XC7H13 <=> C2H3 + C5H10",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.39e+12, 's^-1'), n=-0.58, Ea=(37797.3, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (8e-40, 'cm^3/(mol*s)'),
            n = 15.349,
            Ea = (-603.1, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -38.64,
        T3 = (246, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S3XC7H13 <=> C2H3 + C5H10""",
)

entry(
    index = 1876,
    label = "C7H12-13 + H <=> C5H10 + C2H3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.46e+30, 'cm^3/(mol*s)'),
        n = -4.34,
        Ea = (21647, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C7H12-13 + H <=> C5H10 + C2H3""",
)

entry(
    index = 1877,
    label = "C2H5 + lC5H7 <=> C7H12-13",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.4e+14, 'cm^3/(mol*s)'),
        n = -0.38,
        Ea = (513, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H5 + lC5H7 <=> C7H12-13""",
)

entry(
    index = 1878,
    label = "C7H12-16 <=> aC3H5 + C4H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.14e+23, 's^-1'), n=-2.03, Ea=(74958, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H12-16 <=> aC3H5 + C4H7""",
)

entry(
    index = 1879,
    label = "C7H12-16 + H <=> C2H4 + PXC5H9",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.6e+22, 'cm^3/(mol*s)'),
        n = -2.39,
        Ea = (11180, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C7H12-16 + H <=> C2H4 + PXC5H9""",
)

entry(
    index = 1880,
    label = "CH2-3-1C6H10 + H <=> C5H10 + C2H3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.46e+30, 'cm^3/(mol*s)'),
        n = -4.34,
        Ea = (21647, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2-3-1C6H10 + H <=> C5H10 + C2H3""",
)

entry(
    index = 1881,
    label = "CH2cC6H10 + H <=> PXC7H13",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8e+21, 'cm^3/(mol*s)'),
        n = -2.39,
        Ea = (11180, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2cC6H10 + H <=> PXC7H13""",
)

entry(
    index = 1882,
    label = "cC6H12 <=> C6H12",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.01e+16, 's^-1'), n=0, Ea=(88232, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC6H12 <=> C6H12""",
)

entry(
    index = 1883,
    label = "cC6H12 + H <=> cC6H11 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(394000, 'cm^3/(mol*s)'), n=2.69, Ea=(3837, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC6H12 + H <=> cC6H11 + H2""",
)

entry(
    index = 1884,
    label = "cC6H12 + O <=> cC6H11 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(143000, 'cm^3/(mol*s)'), n=2.71, Ea=(2106, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC6H12 + O <=> cC6H11 + OH""",
)

entry(
    index = 1885,
    label = "cC6H12 + OH <=> cC6H11 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (283000, 'cm^3/(mol*s)'),
        n = 2.88,
        Ea = (-1015, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is cC6H12 + OH <=> cC6H11 + H2O""",
)

entry(
    index = 1886,
    label = "cC6H12 + O2 <=> cC6H11 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+14, 'cm^3/(mol*s)'), n=0, Ea=(47590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC6H12 + O2 <=> cC6H11 + HO2""",
)

entry(
    index = 1887,
    label = "cC6H12 + HO2 <=> cC6H11 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(844, 'cm^3/(mol*s)'), n=3.25, Ea=(14998, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC6H12 + HO2 <=> cC6H11 + H2O2""",
)

entry(
    index = 1888,
    label = "cC6H12 + CH3 <=> cC6H11 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(382, 'cm^3/(mol*s)'), n=3.21, Ea=(11633.6, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC6H12 + CH3 <=> cC6H11 + CH4""",
)

entry(
    index = 1889,
    label = "cC6H11 <=> PXC6H11",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(6.03e+12, 's^-1'), n=0.07, Ea=(27982.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.1e-33, 'cm^3/(mol*s)'),
            n = 15.29,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -25.11,
        T3 = (225, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is cC6H11 <=> PXC6H11""",
)

entry(
    index = 1890,
    label = "cC6H11 <=> cC6H10 + H",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.34e+11, 's^-1'), n=0.69, Ea=(33947.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3e-40, 'cm^3/(mol*s)'),
            n = 17.33,
            Ea = (-602.7, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -19.22,
        T3 = (230, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is cC6H11 <=> cC6H10 + H""",
)

entry(
    index = 1891,
    label = "cC6H11 + H <=> cC6H12",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(4.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.7e+58, 'cm^6/(mol^2*s)'),
            n = -12.08,
            Ea = (11263.7, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.649,
        T3 = (1213.1, 'K'),
        T1 = (1213.1, 'K'),
        T2 = (13369.7, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is cC6H11 + H <=> cC6H12""",
)

entry(
    index = 1892,
    label = "cC6H11 + O2 <=> OH + CO + C5H10",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.19e+13, 'cm^3/(mol*s)'), n=0, Ea=(18056, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC6H11 + O2 <=> OH + CO + C5H10""",
)

entry(
    index = 1893,
    label = "cC6H11 + O2 <=> cC6H10O2H-2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.95e+13, 'cm^3/(mol*s)'), n=0, Ea=(12080, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC6H11 + O2 <=> cC6H10O2H-2""",
)

entry(
    index = 1894,
    label = "cC6H11 + O2 <=> cC6H11O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.17e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1494, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC6H11 + O2 <=> cC6H11O2""",
)

entry(
    index = 1895,
    label = "cC6H11O2 <=> cC6H10O2H-2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.11e+09, 's^-1'), n=0, Ea=(20567, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC6H11O2 <=> cC6H10O2H-2""",
)

entry(
    index = 1896,
    label = "cC6H11O2 <=> cC6H10 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+12, 's^-1'), n=0, Ea=(33000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC6H11O2 <=> cC6H10 + HO2""",
)

entry(
    index = 1897,
    label = "cC6H10O2H-2 <=> cC6H10 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+11, 's^-1'), n=0, Ea=(20000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC6H10O2H-2 <=> cC6H10 + HO2""",
)

entry(
    index = 1898,
    label = "cC6H10O2H-2 => OH + CO + C5H10",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.45e+12, 's^-1'), n=0, Ea=(18060, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC6H10O2H-2 => OH + CO + C5H10""",
)

entry(
    index = 1899,
    label = "cC6H10O2H-2 + O2 <=> SOOcC6O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC6H10O2H-2 + O2 <=> SOOcC6O2H""",
)

entry(
    index = 1900,
    label = "SOOcC6O2H <=> SOOcC6O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 's^-1'), n=0, Ea=(18000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SOOcC6O2H <=> SOOcC6O + OH""",
)

entry(
    index = 1901,
    label = "SOOcC6O => CO + OH + C4H81 + HCO",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(3e+13, 's^-1'), n=0, Ea=(20000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SOOcC6O => CO + OH + C4H81 + HCO""",
)

entry(
    index = 1902,
    label = "PXC6H11 <=> C4H7 + C2H4",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.98e+12, 's^-1'), n=0.12, Ea=(27571.6, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.3e-43, 'cm^3/(mol*s)'),
            n = 18.35,
            Ea = (-602.5, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -13.87,
        T3 = (227, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC6H11 <=> C4H7 + C2H4""",
)

entry(
    index = 1903,
    label = "PXC6H11 <=> SAXC6H11",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(155, 's^-1'), n=2.83, Ea=(15566.2, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.5e-30, 'cm^3/(mol*s)'),
            n = 14.56,
            Ea = (-602.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -13.59,
        T3 = (214, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC6H11 <=> SAXC6H11""",
)

entry(
    index = 1904,
    label = "SAXC6H11 + H <=> C6H12",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.33e+60, 'cm^6/(mol^2*s)'),
            n = -12,
            Ea = (5967.8, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.02,
        T3 = (1096.6, 'K'),
        T1 = (1096.6, 'K'),
        T2 = (6859.5, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SAXC6H11 + H <=> C6H12""",
)

entry(
    index = 1905,
    label = "SAXC6H11 <=> C4H6 + C2H5",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.39e+11, 's^-1'), n=0.66, Ea=(32262.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4e-42, 'cm^3/(mol*s)'),
            n = 18.05,
            Ea = (-602.6, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -18.5,
        T3 = (246, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SAXC6H11 <=> C4H6 + C2H5""",
)

entry(
    index = 1906,
    label = "PXC6H11 <=> PXCH2cC5H9",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(9.55e+08, 's^-1'), n=0.36, Ea=(10704, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (2.3e-28, 'cm^3/(mol*s)'),
            n = 14.28,
            Ea = (-602.5, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -18.98,
        T3 = (214, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC6H11 <=> PXCH2cC5H9""",
)

entry(
    index = 1907,
    label = "PXCH2cC5H9 <=> CH3S3XcC5H8",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(66.1, 's^-1'), n=2.85, Ea=(21082.1, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.1e-36, 'cm^3/(mol*s)'),
            n = 16.12,
            Ea = (-602.5, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -21.57,
        T3 = (225, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXCH2cC5H9 <=> CH3S3XcC5H8""",
)

entry(
    index = 1908,
    label = "CH3S3XcC5H8 <=> SXC6H11",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(7.59e+12, 's^-1'), n=0.15, Ea=(32940.5, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.1e-36, 'cm^3/(mol*s)'),
            n = 16.19,
            Ea = (-602.5, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -21.6,
        T3 = (225, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is CH3S3XcC5H8 <=> SXC6H11""",
)

entry(
    index = 1909,
    label = "CH3S3XcC5H8 <=> PXCH2-4-1C5H9",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(3.55e+12, 's^-1'), n=0.15, Ea=(32404, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.9e-36, 'cm^3/(mol*s)'),
            n = 15.96,
            Ea = (-578.8, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -11.13,
        T3 = (269, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is CH3S3XcC5H8 <=> PXCH2-4-1C5H9""",
)

entry(
    index = 1910,
    label = "PXCH2-4-1C5H9 <=> aC3H5 + C3H6",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(5.37e+12, 's^-1'), n=0.12, Ea=(23947.3, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4.9e-31, 'cm^3/(mol*s)'),
            n = 14.54,
            Ea = (-578.8, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -11.9,
        T3 = (267, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXCH2-4-1C5H9 <=> aC3H5 + C3H6""",
)

entry(
    index = 1911,
    label = "SXC6H11 <=> C3H6 + aC3H5",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(4.57e+12, 's^-1'), n=0.13, Ea=(24386.4, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (2.5e-31, 'cm^3/(mol*s)'),
            n = 14.57,
            Ea = (-578.8, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -13.17,
        T3 = (268, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SXC6H11 <=> C3H6 + aC3H5""",
)

entry(
    index = 1912,
    label = "S2XC6H11 <=> C5H8-14 + CH3",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(8.13e+10, 's^-1'), n=0.78, Ea=(29648, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4e-39, 'cm^3/(mol*s)'),
            n = 16.782,
            Ea = (-600.4, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -7.03,
        T3 = (314, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S2XC6H11 <=> C5H8-14 + CH3""",
)

entry(
    index = 1913,
    label = "S2XC6H11 <=> PXCH2-3-1C5H9",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(7.59e+06, 's^-1'), n=1.81, Ea=(6447.8, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (9.3e-19, 'cm^3/(mol*s)'),
            n = 11.698,
            Ea = (-602.9, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -58.54,
        T3 = (201, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is S2XC6H11 <=> PXCH2-3-1C5H9""",
)

entry(
    index = 1914,
    label = "PXCH2-3-1C5H9 <=> C4H6 + C2H5",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(4.9e+12, 's^-1'), n=0.15, Ea=(11139.1, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.7e-14, 'cm^3/(mol*s)'),
            n = 9.652,
            Ea = (-603.1, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -76.6,
        T3 = (221, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXCH2-3-1C5H9 <=> C4H6 + C2H5""",
)

entry(
    index = 1915,
    label = "PX1-3C6H11 <=> C6H10-13 + H",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(4.27e+09, 's^-1'), n=-0.96, Ea=(31962.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3e-45, 'cm^3/(mol*s)'),
            n = 13.55,
            Ea = (-603.1, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -35.89,
        T3 = (249, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PX1-3C6H11 <=> C6H10-13 + H""",
)

entry(
    index = 1916,
    label = "PAXCH2-2-1C5H9 <=> aC3H4 + nC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 's^-1'), n=0, Ea=(50078, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PAXCH2-2-1C5H9 <=> aC3H4 + nC3H7""",
)

entry(
    index = 1917,
    label = "PX6-2C6H11 <=> SAXC4H7 + C2H4",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.04e+12, 's^-1'), n=-0.37, Ea=(25124.2, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.9e-31, 'cm^3/(mol*s)'),
            n = 13.982,
            Ea = (-603.1, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -14.78,
        T3 = (229, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PX6-2C6H11 <=> SAXC4H7 + C2H4""",
)

entry(
    index = 1918,
    label = "PXC2H4-2-1C4H7 <=> C2H3-2-1C4H7 + H",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(4.27e+09, 's^-1'), n=-0.96, Ea=(31962.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3e-45, 'cm^3/(mol*s)'),
            n = 13.55,
            Ea = (-603.1, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -35.89,
        T3 = (249, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXC2H4-2-1C4H7 <=> C2H3-2-1C4H7 + H""",
)

entry(
    index = 1919,
    label = "C6H10-13 + H <=> C4H81 + C2H3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.46e+30, 'cm^3/(mol*s)'),
        n = -4.34,
        Ea = (21647, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C6H10-13 + H <=> C4H81 + C2H3""",
)

entry(
    index = 1920,
    label = "lC5H7 + CH3 <=> C6H10-13",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(
            A = (1e+14, 'cm^3/(mol*s)'),
            n = -0.32,
            Ea = (-262.3, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        arrheniusLow = Arrhenius(
            A = (3.91e+60, 'cm^6/(mol^2*s)'),
            n = -12.81,
            Ea = (6250, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.104,
        T3 = (1606, 'K'),
        T1 = (60000, 'K'),
        T2 = (6118.4, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is lC5H7 + CH3 <=> C6H10-13""",
)

entry(
    index = 1921,
    label = "C6H10-15 + H <=> C4H7 + C2H4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.6e+22, 'cm^3/(mol*s)'),
        n = -2.39,
        Ea = (11180, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C6H10-15 + H <=> C4H7 + C2H4""",
)

entry(
    index = 1922,
    label = "aC3H5 + aC3H5 <=> C6H10-15",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.02e+13, 'cm^3/(mol*s)'), n=0, Ea=(262, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is aC3H5 + aC3H5 <=> C6H10-15""",
)

entry(
    index = 1923,
    label = "C6H10-15 + H <=> SAXC6H9-15 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(330000, 'cm^3/(mol*s)'), n=2.5, Ea=(2490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H10-15 + H <=> SAXC6H9-15 + H2""",
)

entry(
    index = 1924,
    label = "C6H10-15 + OH <=> SAXC6H9-15 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.14e+06, 'cm^3/(mol*s)'), n=2, Ea=(-298, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H10-15 + OH <=> SAXC6H9-15 + H2O""",
)

entry(
    index = 1925,
    label = "C6H10-15 + CH3 <=> SAXC6H9-15 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.94, 'cm^3/(mol*s)'), n=3.5, Ea=(5675, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H10-15 + CH3 <=> SAXC6H9-15 + CH4""",
)

entry(
    index = 1926,
    label = "C4H6 + C2H3 <=> SAXC6H9-15",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.2e+13, 'cm^3/(mol*s)'), n=0, Ea=(1300, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H6 + C2H3 <=> SAXC6H9-15""",
)

entry(
    index = 1927,
    label = "C6H10-12 + H <=> C4H6 + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(2000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H10-12 + H <=> C4H6 + C2H5""",
)

entry(
    index = 1928,
    label = "cC6H10 <=> C4H6 + C2H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.5e+12, 's^-1'), n=0.76, Ea=(62450, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC6H10 <=> C4H6 + C2H4""",
)

entry(
    index = 1929,
    label = "cC6H10 <=> cC6H8-13 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.08e+09, 's^-1'), n=1.12, Ea=(59560, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC6H10 <=> cC6H8-13 + H2""",
)

entry(
    index = 1930,
    label = "cC6H10 + H <=> SAXcC6H9 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(230000, 'cm^3/(mol*s)'), n=2.5, Ea=(2490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC6H10 + H <=> SAXcC6H9 + H2""",
)

entry(
    index = 1931,
    label = "cC6H10 + O <=> SAXcC6H9 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+11, 'cm^3/(mol*s)'), n=0.7, Ea=(5880, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC6H10 + O <=> SAXcC6H9 + OH""",
)

entry(
    index = 1932,
    label = "cC6H10 + OH <=> SAXcC6H9 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.14e+06, 'cm^3/(mol*s)'), n=2, Ea=(-298, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC6H10 + OH <=> SAXcC6H9 + H2O""",
)

entry(
    index = 1933,
    label = "SAXcC6H9 + HO2 <=> cC6H10 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.66e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SAXcC6H9 + HO2 <=> cC6H10 + O2""",
)

entry(
    index = 1934,
    label = "cC6H10 + CH3 <=> SAXcC6H9 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.94, 'cm^3/(mol*s)'), n=3.5, Ea=(5675, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC6H10 + CH3 <=> SAXcC6H9 + CH4""",
)

entry(
    index = 1935,
    label = "SAXcC6H9 + H <=> cC6H10",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.33e+60, 'cm^6/(mol^2*s)'),
            n = -12,
            Ea = (5967.8, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.02,
        T3 = (1096.6, 'K'),
        T1 = (1096.6, 'K'),
        T2 = (6859.5, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SAXcC6H9 + H <=> cC6H10""",
)

entry(
    index = 1936,
    label = "SAXcC6H9 <=> cC6H8-13 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.67e+12, 's^-1'), n=0.71, Ea=(49792.2, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SAXcC6H9 <=> cC6H8-13 + H""",
)

entry(
    index = 1937,
    label = "SAXcC6H9 + H <=> cC6H8-13 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SAXcC6H9 + H <=> cC6H8-13 + H2""",
)

entry(
    index = 1938,
    label = "SAXcC6H9 + O2 <=> cC6H8-13 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SAXcC6H9 + O2 <=> cC6H8-13 + HO2""",
)

entry(
    index = 1939,
    label = "SAXcC6H9 + CH3 <=> cC6H8-13 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SAXcC6H9 + CH3 <=> cC6H8-13 + CH4""",
)

entry(
    index = 1940,
    label = "cC6H8-13 <=> C6H6 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.51e+13, 's^-1'), n=0, Ea=(59020, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC6H8-13 <=> C6H6 + H2""",
)

entry(
    index = 1941,
    label = "cC6H8-13 + H <=> SAXcC6H7 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(230000, 'cm^3/(mol*s)'), n=2.5, Ea=(2490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC6H8-13 + H <=> SAXcC6H7 + H2""",
)

entry(
    index = 1942,
    label = "cC6H8-13 + O <=> SAXcC6H7 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+11, 'cm^3/(mol*s)'), n=0.7, Ea=(5880, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC6H8-13 + O <=> SAXcC6H7 + OH""",
)

entry(
    index = 1943,
    label = "cC6H8-13 + OH <=> SAXcC6H7 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.14e+06, 'cm^3/(mol*s)'), n=2, Ea=(-298, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC6H8-13 + OH <=> SAXcC6H7 + H2O""",
)

entry(
    index = 1944,
    label = "SAXcC6H7 + HO2 <=> cC6H8-13 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.66e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SAXcC6H7 + HO2 <=> cC6H8-13 + O2""",
)

entry(
    index = 1945,
    label = "cC6H8-13 + CH3 <=> SAXcC6H7 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.94, 'cm^3/(mol*s)'), n=3.5, Ea=(5675, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is cC6H8-13 + CH3 <=> SAXcC6H7 + CH4""",
)

entry(
    index = 1946,
    label = "SAXcC6H7 + H <=> cC6H8-13",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.33e+60, 'cm^6/(mol^2*s)'),
            n = -12,
            Ea = (5967.8, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.02,
        T3 = (1096.6, 'K'),
        T1 = (1096.6, 'K'),
        T2 = (6859.5, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is SAXcC6H7 + H <=> cC6H8-13""",
)

entry(
    index = 1947,
    label = "SAXcC6H7 <=> C6H6 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.76e+11, 's^-1'), n=0.78, Ea=(30230, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SAXcC6H7 <=> C6H6 + H""",
)

entry(
    index = 1948,
    label = "SAXcC6H7 + H <=> C6H6 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SAXcC6H7 + H <=> C6H6 + H2""",
)

entry(
    index = 1949,
    label = "SAXcC6H7 + O2 <=> C6H6 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SAXcC6H7 + O2 <=> C6H6 + HO2""",
)

entry(
    index = 1950,
    label = "SAXcC6H7 + CH3 <=> C6H6 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SAXcC6H7 + CH3 <=> C6H6 + CH4""",
)

entry(
    index = 1951,
    label = "PXCH2-2-C4H9 <=> C3H6 + C2H5",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(5.47e+11, 's^-1'), n=0.5, Ea=(27798.1, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (8.5e-37, 'cm^3/(mol*s)'),
            n = 16.21,
            Ea = (-603.1, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -26.7,
        T3 = (216, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is PXCH2-2-C4H9 <=> C3H6 + C2H5""",
)

entry(
    index = 1952,
    label = "PAXCH2-2-1C4H7 <=> aC3H4 + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 's^-1'), n=0, Ea=(50078, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PAXCH2-2-1C4H7 <=> aC3H4 + C2H5""",
)

entry(
    index = 1953,
    label = "CH3-2-PXC4H6 <=> CH3-2-C4H5-13 + H",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(4.27e+09, 's^-1'), n=-0.96, Ea=(31962.9, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3e-45, 'cm^3/(mol*s)'),
            n = 13.55,
            Ea = (-603.1, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = -35.89,
        T3 = (249, 'K'),
        T1 = (28, 'K'),
        T2 = (50000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is CH3-2-PXC4H6 <=> CH3-2-C4H5-13 + H""",
)

entry(
    index = 1954,
    label = "C2H4 + CH3CCH2 <=> CH3-2-PXC4H6",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (7.93e+38, 'cm^3/(mol*s)'),
        n = -8.47,
        Ea = (14220, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H4 + CH3CCH2 <=> CH3-2-PXC4H6""",
)

entry(
    index = 1955,
    label = "C5H8-13 + H <=> C2H4 + CH3CHCH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.46e+30, 'cm^3/(mol*s)'),
        n = -4.34,
        Ea = (21647, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C5H8-13 + H <=> C2H4 + CH3CHCH""",
)

entry(
    index = 1956,
    label = "C5H8-13 + H <=> C4H6-2 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+12, 'cm^3/(mol*s)'), n=0, Ea=(7000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H8-13 + H <=> C4H6-2 + CH3""",
)

entry(
    index = 1957,
    label = "C5H8-13 + H <=> C4H612 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+12, 'cm^3/(mol*s)'), n=0, Ea=(7000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H8-13 + H <=> C4H612 + CH3""",
)

entry(
    index = 1958,
    label = "C5H8-13 + H <=> lC5H7 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(173000, 'cm^3/(mol*s)'), n=2.5, Ea=(2490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H8-13 + H <=> lC5H7 + H2""",
)

entry(
    index = 1959,
    label = "C5H8-13 + OH <=> lC5H7 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.1e+06, 'cm^3/(mol*s)'), n=2, Ea=(-298, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H8-13 + OH <=> lC5H7 + H2O""",
)

entry(
    index = 1960,
    label = "C5H8-13 + CH3 <=> lC5H7 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.2, 'cm^3/(mol*s)'), n=3.5, Ea=(5675, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H8-13 + CH3 <=> lC5H7 + CH4""",
)

entry(
    index = 1961,
    label = "C5H8-12 <=> C5H8-13",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.2e+14, 's^-1'), n=0, Ea=(67000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H8-12 <=> C5H8-13""",
)

entry(
    index = 1962,
    label = "C5H8-12 <=> C2H4 + pC3H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.6e+12, 's^-1'), n=0, Ea=(58100, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H8-12 <=> C2H4 + pC3H4""",
)

entry(
    index = 1963,
    label = "nC4H5 + CH3 <=> C5H8-13",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.23e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is nC4H5 + CH3 <=> C5H8-13""",
)

entry(
    index = 1964,
    label = "C5H8-14 + H <=> aC3H5 + C2H4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.6e+22, 'cm^3/(mol*s)'),
        n = -2.39,
        Ea = (11180, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C5H8-14 + H <=> aC3H5 + C2H4""",
)

entry(
    index = 1965,
    label = "C5H8-14 + H <=> lC5H7 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(115000, 'cm^3/(mol*s)'), n=2.5, Ea=(2490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H8-14 + H <=> lC5H7 + H2""",
)

entry(
    index = 1966,
    label = "C5H8-14 + H <=> C4H6 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+12, 'cm^3/(mol*s)'), n=0, Ea=(7000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H8-14 + H <=> C4H6 + CH3""",
)

entry(
    index = 1967,
    label = "C5H8-14 + OH <=> lC5H7 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.07e+06, 'cm^3/(mol*s)'), n=2, Ea=(-298, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H8-14 + OH <=> lC5H7 + H2O""",
)

entry(
    index = 1968,
    label = "C5H8-14 + CH3 <=> lC5H7 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.47, 'cm^3/(mol*s)'), n=3.5, Ea=(5675, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H8-14 + CH3 <=> lC5H7 + CH4""",
)

entry(
    index = 1969,
    label = "CH3-2-C4H5-13 + H <=> C3H6 + C2H3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.46e+30, 'cm^3/(mol*s)'),
        n = -4.34,
        Ea = (21647, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3-2-C4H5-13 + H <=> C3H6 + C2H3""",
)

entry(
    index = 1970,
    label = "CH3-2-C4H5-13 + H <=> C4H6 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8e+21, 'cm^3/(mol*s)'),
        n = -2.39,
        Ea = (11180, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3-2-C4H5-13 + H <=> C4H6 + CH3""",
)

entry(
    index = 1971,
    label = "CH3-2-C4H5-13 + H <=> PAXCH2-2-C4H5 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.2e+06, 'cm^3/(mol*s)'),
        n = 2.54,
        Ea = (6760, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3-2-C4H5-13 + H <=> PAXCH2-2-C4H5 + H2""",
)

entry(
    index = 1972,
    label = "CH3-2-C4H5-13 + OH <=> PAXCH2-2-C4H5 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.1e+06, 'cm^3/(mol*s)'), n=2, Ea=(-298, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3-2-C4H5-13 + OH <=> PAXCH2-2-C4H5 + H2O""",
)

entry(
    index = 1973,
    label = "CH3-2-C4H5-13 + CH3 <=> PAXCH2-2-C4H5 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.2, 'cm^3/(mol*s)'), n=3.5, Ea=(5675, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3-2-C4H5-13 + CH3 <=> PAXCH2-2-C4H5 + CH4""",
)

entry(
    index = 1974,
    label = "CH3 + iC4H5 <=> CH3-2-C4H5-13",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + iC4H5 <=> CH3-2-C4H5-13""",
)

entry(
    index = 1975,
    label = "PAXCH2-2-C4H5 <=> aC3H4 + C2H3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 's^-1'), n=0, Ea=(50078, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PAXCH2-2-C4H5 <=> aC3H4 + C2H3""",
)

entry(
    index = 1976,
    label = "SAXC4H7 <=> C4H6 + H",
    degeneracy = 1,
    kinetics = Lindemann(
        arrheniusHigh = Arrhenius(A=(4.7e+08, 's^-1'), n=1.32, Ea=(44697.6, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4.6e-37, 'cm^3/(mol*s)'),
            n = 15.37,
            Ea = (-603.1, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5},
    ),
    shortDesc = u"""The chemkin file reaction is SAXC4H7 <=> C4H6 + H""",
)

entry(
    index = 1977,
    label = "SAXC4H7 + H <=> C4H81",
    degeneracy = 1,
    kinetics = Lindemann(
        arrheniusHigh = Arrhenius(A=(2e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (2.66e+60, 'cm^6/(mol^2*s)'),
            n = -12,
            Ea = (5967.8, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5},
    ),
    shortDesc = u"""The chemkin file reaction is SAXC4H7 + H <=> C4H81""",
)

entry(
    index = 1978,
    label = "SAXC4H7 + H <=> C4H82",
    degeneracy = 1,
    kinetics = Lindemann(
        arrheniusHigh = Arrhenius(A=(2e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (2.66e+60, 'cm^6/(mol^2*s)'),
            n = -12,
            Ea = (5967.8, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C-]#[O+]': 1.5},
    ),
    shortDesc = u"""The chemkin file reaction is SAXC4H7 + H <=> C4H82""",
)

entry(
    index = 1979,
    label = "SAXC4H7 + H <=> C4H6 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SAXC4H7 + H <=> C4H6 + H2""",
)

entry(
    index = 1980,
    label = "SAXC4H7 + H <=> C4H612 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SAXC4H7 + H <=> C4H612 + H2""",
)

entry(
    index = 1981,
    label = "C4H7 + H <=> SAXC4H7 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H7 + H <=> SAXC4H7 + H""",
)

entry(
    index = 1982,
    label = "SAXC4H7 + O <=> C2H3 + CH3CHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SAXC4H7 + O <=> C2H3 + CH3CHO""",
)

entry(
    index = 1983,
    label = "SAXC4H7 + OH <=> C4H6 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SAXC4H7 + OH <=> C4H6 + H2O""",
)

entry(
    index = 1984,
    label = "SAXC4H7 + OH <=> C4H612 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SAXC4H7 + OH <=> C4H612 + H2O""",
)

entry(
    index = 1985,
    label = "SAXC4H7 + HO2 <=> C4H81 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SAXC4H7 + HO2 <=> C4H81 + O2""",
)

entry(
    index = 1986,
    label = "SAXC4H7 + HO2 <=> C4H82 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SAXC4H7 + HO2 <=> C4H82 + O2""",
)

entry(
    index = 1987,
    label = "SAXC4H7 + HCO <=> C4H81 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SAXC4H7 + HCO <=> C4H81 + CO""",
)

entry(
    index = 1988,
    label = "SAXC4H7 + HCO <=> C4H82 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SAXC4H7 + HCO <=> C4H82 + CO""",
)

entry(
    index = 1989,
    label = "SAXC4H7 + CH3 <=> C4H6 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SAXC4H7 + CH3 <=> C4H6 + CH4""",
)

entry(
    index = 1990,
    label = "SAXC4H7 + CH3 <=> C4H612 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SAXC4H7 + CH3 <=> C4H612 + CH4""",
)

entry(
    index = 1991,
    label = "SAXC4H7 + C2H3 <=> C4H6 + C2H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SAXC4H7 + C2H3 <=> C4H6 + C2H4""",
)

entry(
    index = 1992,
    label = "SAXC4H7 + C2H3 <=> C4H612 + C2H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SAXC4H7 + C2H3 <=> C4H612 + C2H4""",
)

entry(
    index = 1993,
    label = "CH3 + aC3H4 <=> SAXC4H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(7500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + aC3H4 <=> SAXC4H7""",
)

entry(
    index = 1994,
    label = "CH3 + pC3H4 <=> SAXC4H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(7500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + pC3H4 <=> SAXC4H7""",
)

entry(
    index = 1995,
    label = "C4H81 + H <=> SAXC4H7 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+06, 'cm^3/(mol*s)'), n=2.4, Ea=(4471, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H81 + H <=> SAXC4H7 + H2""",
)

entry(
    index = 1996,
    label = "C4H81 + OH <=> SAXC4H7 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(27000, 'cm^3/(mol*s)'), n=2.39, Ea=(393, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H81 + OH <=> SAXC4H7 + H2O""",
)

entry(
    index = 1997,
    label = "C4H81 + O <=> SAXC4H7 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(47600, 'cm^3/(mol*s)'), n=2.71, Ea=(2106, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H81 + O <=> SAXC4H7 + OH""",
)

entry(
    index = 1998,
    label = "C4H81 + HO2 <=> SAXC4H7 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7130, 'cm^3/(mol*s)'), n=2.77, Ea=(14913, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H81 + HO2 <=> SAXC4H7 + H2O2""",
)

entry(
    index = 1999,
    label = "C4H81 + CH3 <=> SAXC4H7 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(7300, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H81 + CH3 <=> SAXC4H7 + CH4""",
)

entry(
    index = 2000,
    label = "C4H82 + H <=> SAXC4H7 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(346000, 'cm^3/(mol*s)'), n=2.5, Ea=(2490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H82 + H <=> SAXC4H7 + H2""",
)

entry(
    index = 2001,
    label = "C4H82 + O <=> SAXC4H7 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.6e+11, 'cm^3/(mol*s)'), n=0.7, Ea=(5880, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H82 + O <=> SAXC4H7 + OH""",
)

entry(
    index = 2002,
    label = "C4H82 + OH <=> SAXC4H7 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.2e+06, 'cm^3/(mol*s)'), n=2, Ea=(-298, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H82 + OH <=> SAXC4H7 + H2O""",
)

entry(
    index = 2003,
    label = "C4H82 + HO2 <=> SAXC4H7 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(14200, 'cm^3/(mol*s)'), n=2.77, Ea=(14913, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H82 + HO2 <=> SAXC4H7 + H2O2""",
)

entry(
    index = 2004,
    label = "C4H82 + CH3 <=> SAXC4H7 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.4, 'cm^3/(mol*s)'), n=3.5, Ea=(5675, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H82 + CH3 <=> SAXC4H7 + CH4""",
)

entry(
    index = 2005,
    label = "PXCH2cC6H11 + HO2 => OH + CH2O + cC6H11",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXCH2cC6H11 + HO2 => OH + CH2O + cC6H11""",
)

entry(
    index = 2006,
    label = "CH3TXcC6H10 + HO2 => OH + CH3CHO + cC5H9",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3TXcC6H10 + HO2 => OH + CH3CHO + cC5H9""",
)

entry(
    index = 2007,
    label = "CH3S2XcC6H10 + HO2 => OH + CH3CH2CHO + C4H7",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3S2XcC6H10 + HO2 => OH + CH3CH2CHO + C4H7""",
)

entry(
    index = 2008,
    label = "CH3S3XcC6H10 + HO2 => OH + CH3CH2CHO + C4H7",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3S3XcC6H10 + HO2 => OH + CH3CH2CHO + C4H7""",
)

entry(
    index = 2009,
    label = "CH3S4XcC6H10 + HO2 => OH + CH3CH2CHO + C4H7",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3S4XcC6H10 + HO2 => OH + CH3CH2CHO + C4H7""",
)

entry(
    index = 2010,
    label = "PXCH2cC6H11 + O2 => HCO + C3H6 + aC3H5 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(7.19e+13, 'cm^3/(mol*s)'), n=0, Ea=(18056, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXCH2cC6H11 + O2 => HCO + C3H6 + aC3H5 + OH""",
)

entry(
    index = 2011,
    label = "CH3TXcC6H10 + O2 => HCO + C3H6 + aC3H5 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(7.19e+13, 'cm^3/(mol*s)'), n=0, Ea=(18056, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3TXcC6H10 + O2 => HCO + C3H6 + aC3H5 + OH""",
)

entry(
    index = 2012,
    label = "CH3S2XcC6H10 + O2 => HCO + C3H6 + aC3H5 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(7.19e+13, 'cm^3/(mol*s)'), n=0, Ea=(18056, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3S2XcC6H10 + O2 => HCO + C3H6 + aC3H5 + OH""",
)

entry(
    index = 2013,
    label = "CH3S3XcC6H10 + O2 => HCO + C3H6 + aC3H5 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(7.19e+13, 'cm^3/(mol*s)'), n=0, Ea=(18056, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3S3XcC6H10 + O2 => HCO + C3H6 + aC3H5 + OH""",
)

entry(
    index = 2014,
    label = "CH3S4XcC6H10 + O2 => HCO + C3H6 + aC3H5 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(7.19e+13, 'cm^3/(mol*s)'), n=0, Ea=(18056, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3S4XcC6H10 + O2 => HCO + C3H6 + aC3H5 + OH""",
)

entry(
    index = 2015,
    label = "PXCH2cC6H11 + O2 <=> CH3cC6H10OO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.17e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1494, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXCH2cC6H11 + O2 <=> CH3cC6H10OO""",
)

entry(
    index = 2016,
    label = "CH3TXcC6H10 + O2 <=> CH3cC6H10OO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.17e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1494, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3TXcC6H10 + O2 <=> CH3cC6H10OO""",
)

entry(
    index = 2017,
    label = "CH3S2XcC6H10 + O2 <=> CH3cC6H10OO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.17e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1494, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3S2XcC6H10 + O2 <=> CH3cC6H10OO""",
)

entry(
    index = 2018,
    label = "CH3S3XcC6H10 + O2 <=> CH3cC6H10OO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.17e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1494, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3S3XcC6H10 + O2 <=> CH3cC6H10OO""",
)

entry(
    index = 2019,
    label = "CH3S4XcC6H10 + O2 <=> CH3cC6H10OO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.17e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1494, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3S4XcC6H10 + O2 <=> CH3cC6H10OO""",
)

entry(
    index = 2020,
    label = "PXCH2cC6H11 + O2 <=> CH3cC6H9OOH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.95e+13, 'cm^3/(mol*s)'), n=0, Ea=(12080, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXCH2cC6H11 + O2 <=> CH3cC6H9OOH""",
)

entry(
    index = 2021,
    label = "CH3TXcC6H10 + O2 <=> CH3cC6H9OOH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.95e+13, 'cm^3/(mol*s)'), n=0, Ea=(12080, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3TXcC6H10 + O2 <=> CH3cC6H9OOH""",
)

entry(
    index = 2022,
    label = "CH3S2XcC6H10 + O2 <=> CH3cC6H9OOH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.95e+13, 'cm^3/(mol*s)'), n=0, Ea=(12080, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3S2XcC6H10 + O2 <=> CH3cC6H9OOH""",
)

entry(
    index = 2023,
    label = "CH3S3XcC6H10 + O2 <=> CH3cC6H9OOH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.95e+13, 'cm^3/(mol*s)'), n=0, Ea=(12080, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3S3XcC6H10 + O2 <=> CH3cC6H9OOH""",
)

entry(
    index = 2024,
    label = "CH3S4XcC6H10 + O2 <=> CH3cC6H9OOH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.95e+13, 'cm^3/(mol*s)'), n=0, Ea=(12080, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3S4XcC6H10 + O2 <=> CH3cC6H9OOH""",
)

entry(
    index = 2025,
    label = "CH3cC6H10OO <=> CH3cC6H9OOH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.11e+09, 's^-1'), n=0, Ea=(20567, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3cC6H10OO <=> CH3cC6H9OOH""",
)

entry(
    index = 2026,
    label = "CH3cC6H10OO <=> CH2cC6H10 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+12, 's^-1'), n=0, Ea=(33000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3cC6H10OO <=> CH2cC6H10 + HO2""",
)

entry(
    index = 2027,
    label = "CH3cC6H9OOH <=> CH2cC6H10 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+11, 's^-1'), n=0, Ea=(13000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3cC6H9OOH <=> CH2cC6H10 + HO2""",
)

entry(
    index = 2028,
    label = "CH3cC6H9OOH => HCO + C3H6 + aC3H5 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.45e+12, 's^-1'), n=0, Ea=(18060, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3cC6H9OOH => HCO + C3H6 + aC3H5 + OH""",
)

entry(
    index = 2029,
    label = "CH3cC6H9OOH + O2 => CH3cC6H9O3 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(4e+12, 'cm^3/(mol*s)'), n=0, Ea=(10000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3cC6H9OOH + O2 => CH3cC6H9O3 + OH""",
)

entry(
    index = 2030,
    label = "CH3cC6H9O3 => OH + CH2CHO + CH2CO + C3H6",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(3e+13, 's^-1'), n=0, Ea=(20000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3cC6H9O3 => OH + CH2CHO + CH2CO + C3H6""",
)

entry(
    index = 2031,
    label = "PXC2H4cC6H11 + HO2 => OH + CH3CHO + cC6H11",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC2H4cC6H11 + HO2 => OH + CH3CHO + cC6H11""",
)

entry(
    index = 2032,
    label = "SXC2H4cC6H11 + HO2 => OH + CH3CHO + cC6H11",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC2H4cC6H11 + HO2 => OH + CH3CHO + cC6H11""",
)

entry(
    index = 2033,
    label = "C2H5TXcC6H10 + HO2 => OH + CH3CH2CHO + cC5H9",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5TXcC6H10 + HO2 => OH + CH3CH2CHO + cC5H9""",
)

entry(
    index = 2034,
    label = "C2H5S2XcC6H10 + HO2 => OH + C2H4 + CH3CHO + C4H7",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5S2XcC6H10 + HO2 => OH + C2H4 + CH3CHO + C4H7""",
)

entry(
    index = 2035,
    label = "C2H5S3XcC6H10 + HO2 => OH + C2H4 + CH3CHO + C4H7",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5S3XcC6H10 + HO2 => OH + C2H4 + CH3CHO + C4H7""",
)

entry(
    index = 2036,
    label = "C2H5S4XcC6H10 + HO2 => OH + C2H4 + CH3CHO + C4H7",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5S4XcC6H10 + HO2 => OH + C2H4 + CH3CHO + C4H7""",
)

entry(
    index = 2037,
    label = "PXC2H4cC6H11 + O2 => OH + HCO + C4H81 + aC3H5",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(7.19e+13, 'cm^3/(mol*s)'), n=0, Ea=(18056, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC2H4cC6H11 + O2 => OH + HCO + C4H81 + aC3H5""",
)

entry(
    index = 2038,
    label = "SXC2H4cC6H11 + O2 => OH + HCO + C4H81 + aC3H5",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(7.19e+13, 'cm^3/(mol*s)'), n=0, Ea=(18056, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC2H4cC6H11 + O2 => OH + HCO + C4H81 + aC3H5""",
)

entry(
    index = 2039,
    label = "C2H5TXcC6H10 + O2 => HCO + C4H81 + aC3H5 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(7.19e+13, 'cm^3/(mol*s)'), n=0, Ea=(18056, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5TXcC6H10 + O2 => HCO + C4H81 + aC3H5 + OH""",
)

entry(
    index = 2040,
    label = "C2H5S2XcC6H10 + O2 => HCO + C4H81 + aC3H5 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(7.19e+13, 'cm^3/(mol*s)'), n=0, Ea=(18056, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5S2XcC6H10 + O2 => HCO + C4H81 + aC3H5 + OH""",
)

entry(
    index = 2041,
    label = "C2H5S3XcC6H10 + O2 => HCO + C4H81 + aC3H5 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(7.19e+13, 'cm^3/(mol*s)'), n=0, Ea=(18056, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5S3XcC6H10 + O2 => HCO + C4H81 + aC3H5 + OH""",
)

entry(
    index = 2042,
    label = "C2H5S4XcC6H10 + O2 => HCO + C4H81 + aC3H5 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(7.19e+13, 'cm^3/(mol*s)'), n=0, Ea=(18056, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5S4XcC6H10 + O2 => HCO + C4H81 + aC3H5 + OH""",
)

entry(
    index = 2043,
    label = "PXC2H4cC6H11 + O2 <=> C2H5cC6H10OO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.17e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1494, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC2H4cC6H11 + O2 <=> C2H5cC6H10OO""",
)

entry(
    index = 2044,
    label = "SXC2H4cC6H11 + O2 <=> C2H5cC6H10OO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.17e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1494, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC2H4cC6H11 + O2 <=> C2H5cC6H10OO""",
)

entry(
    index = 2045,
    label = "C2H5TXcC6H10 + O2 <=> C2H5cC6H10OO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.17e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1494, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5TXcC6H10 + O2 <=> C2H5cC6H10OO""",
)

entry(
    index = 2046,
    label = "C2H5S2XcC6H10 + O2 <=> C2H5cC6H10OO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.17e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1494, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5S2XcC6H10 + O2 <=> C2H5cC6H10OO""",
)

entry(
    index = 2047,
    label = "C2H5S3XcC6H10 + O2 <=> C2H5cC6H10OO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.17e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1494, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5S3XcC6H10 + O2 <=> C2H5cC6H10OO""",
)

entry(
    index = 2048,
    label = "C2H5S4XcC6H10 + O2 <=> C2H5cC6H10OO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.17e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1494, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5S4XcC6H10 + O2 <=> C2H5cC6H10OO""",
)

entry(
    index = 2049,
    label = "PXC2H4cC6H11 + O2 <=> C2H5cC6H9OOH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.95e+13, 'cm^3/(mol*s)'), n=0, Ea=(12080, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC2H4cC6H11 + O2 <=> C2H5cC6H9OOH""",
)

entry(
    index = 2050,
    label = "SXC2H4cC6H11 + O2 <=> C2H5cC6H9OOH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.95e+13, 'cm^3/(mol*s)'), n=0, Ea=(12080, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC2H4cC6H11 + O2 <=> C2H5cC6H9OOH""",
)

entry(
    index = 2051,
    label = "C2H5TXcC6H10 + O2 <=> C2H5cC6H9OOH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.95e+13, 'cm^3/(mol*s)'), n=0, Ea=(12080, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5TXcC6H10 + O2 <=> C2H5cC6H9OOH""",
)

entry(
    index = 2052,
    label = "C2H5S2XcC6H10 + O2 <=> C2H5cC6H9OOH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.95e+13, 'cm^3/(mol*s)'), n=0, Ea=(12080, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5S2XcC6H10 + O2 <=> C2H5cC6H9OOH""",
)

entry(
    index = 2053,
    label = "C2H5S3XcC6H10 + O2 <=> C2H5cC6H9OOH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.95e+13, 'cm^3/(mol*s)'), n=0, Ea=(12080, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5S3XcC6H10 + O2 <=> C2H5cC6H9OOH""",
)

entry(
    index = 2054,
    label = "C2H5S4XcC6H10 + O2 <=> C2H5cC6H9OOH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.95e+13, 'cm^3/(mol*s)'), n=0, Ea=(12080, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5S4XcC6H10 + O2 <=> C2H5cC6H9OOH""",
)

entry(
    index = 2055,
    label = "C2H5cC6H10OO <=> C2H5cC6H9OOH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.11e+09, 's^-1'), n=0, Ea=(20567, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5cC6H10OO <=> C2H5cC6H9OOH""",
)

entry(
    index = 2056,
    label = "C2H5cC6H10OO <=> C2H4 + cC6H10 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+12, 's^-1'), n=0, Ea=(33000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5cC6H10OO <=> C2H4 + cC6H10 + HO2""",
)

entry(
    index = 2057,
    label = "C2H5cC6H9OOH <=> C2H4 + cC6H10 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+11, 's^-1'), n=0, Ea=(13000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5cC6H9OOH <=> C2H4 + cC6H10 + HO2""",
)

entry(
    index = 2058,
    label = "C2H5cC6H9OOH => HCO + C4H81 + aC3H5 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.45e+12, 's^-1'), n=0, Ea=(18060, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5cC6H9OOH => HCO + C4H81 + aC3H5 + OH""",
)

entry(
    index = 2059,
    label = "C2H5cC6H9OOH + O2 => C2H5cC6H9O3 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(4e+12, 'cm^3/(mol*s)'), n=0, Ea=(10000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5cC6H9OOH + O2 => C2H5cC6H9O3 + OH""",
)

entry(
    index = 2060,
    label = "C2H5cC6H9O3 => OH + CH2CHO + CH2CO + C4H81",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(3e+13, 's^-1'), n=0, Ea=(20000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5cC6H9O3 => OH + CH2CHO + CH2CO + C4H81""",
)

entry(
    index = 2061,
    label = "PXC3H6cC6H11 + HO2 => OH + CH3CH2CHO + cC6H11",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC3H6cC6H11 + HO2 => OH + CH3CH2CHO + cC6H11""",
)

entry(
    index = 2062,
    label = "SXC3H6cC6H11 + HO2 => OH + CH3CH2CHO + cC6H11",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC3H6cC6H11 + HO2 => OH + CH3CH2CHO + cC6H11""",
)

entry(
    index = 2063,
    label = "S2XC3H6cC6H11 + HO2 => OH + CH3CH2CHO + cC6H11",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is S2XC3H6cC6H11 + HO2 => OH + CH3CH2CHO + cC6H11""",
)

entry(
    index = 2064,
    label = "C3H7TXcC6H10 + HO2 => OH + CH3CH2CHO + cC6H11",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7TXcC6H10 + HO2 => OH + CH3CH2CHO + cC6H11""",
)

entry(
    index = 2065,
    label = "C3H7S2XcC6H10 + HO2 => OH + CH3CH2CHO + cC6H11",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7S2XcC6H10 + HO2 => OH + CH3CH2CHO + cC6H11""",
)

entry(
    index = 2066,
    label = "C3H7S3XcC6H10 + HO2 => OH + CH3CH2CHO + cC6H11",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7S3XcC6H10 + HO2 => OH + CH3CH2CHO + cC6H11""",
)

entry(
    index = 2067,
    label = "C3H7S4XcC6H10 + HO2 => OH + CH3CH2CHO + cC6H11",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7S4XcC6H10 + HO2 => OH + CH3CH2CHO + cC6H11""",
)

entry(
    index = 2068,
    label = "PXC3H6cC6H11 + O2 => HCO + C5H10 + aC3H5 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(7.19e+13, 'cm^3/(mol*s)'), n=0, Ea=(18056, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC3H6cC6H11 + O2 => HCO + C5H10 + aC3H5 + OH""",
)

entry(
    index = 2069,
    label = "SXC3H6cC6H11 + O2 => HCO + C5H10 + aC3H5 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(7.19e+13, 'cm^3/(mol*s)'), n=0, Ea=(18056, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC3H6cC6H11 + O2 => HCO + C5H10 + aC3H5 + OH""",
)

entry(
    index = 2070,
    label = "S2XC3H6cC6H11 + O2 => HCO + C5H10 + aC3H5 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(7.19e+13, 'cm^3/(mol*s)'), n=0, Ea=(18056, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is S2XC3H6cC6H11 + O2 => HCO + C5H10 + aC3H5 + OH""",
)

entry(
    index = 2071,
    label = "C3H7TXcC6H10 + O2 => HCO + C5H10 + aC3H5 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(7.19e+13, 'cm^3/(mol*s)'), n=0, Ea=(18056, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7TXcC6H10 + O2 => HCO + C5H10 + aC3H5 + OH""",
)

entry(
    index = 2072,
    label = "C3H7S2XcC6H10 + O2 => HCO + C5H10 + aC3H5 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(7.19e+13, 'cm^3/(mol*s)'), n=0, Ea=(18056, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7S2XcC6H10 + O2 => HCO + C5H10 + aC3H5 + OH""",
)

entry(
    index = 2073,
    label = "C3H7S3XcC6H10 + O2 => HCO + C5H10 + aC3H5 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(7.19e+13, 'cm^3/(mol*s)'), n=0, Ea=(18056, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7S3XcC6H10 + O2 => HCO + C5H10 + aC3H5 + OH""",
)

entry(
    index = 2074,
    label = "C3H7S4XcC6H10 + O2 => HCO + C5H10 + aC3H5 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(7.19e+13, 'cm^3/(mol*s)'), n=0, Ea=(18056, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7S4XcC6H10 + O2 => HCO + C5H10 + aC3H5 + OH""",
)

entry(
    index = 2075,
    label = "PXC3H6cC6H11 + O2 <=> C3H7cC6H10OO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.17e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1494, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC3H6cC6H11 + O2 <=> C3H7cC6H10OO""",
)

entry(
    index = 2076,
    label = "SXC3H6cC6H11 + O2 <=> C3H7cC6H10OO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.17e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1494, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC3H6cC6H11 + O2 <=> C3H7cC6H10OO""",
)

entry(
    index = 2077,
    label = "S2XC3H6cC6H11 + O2 <=> C3H7cC6H10OO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.17e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1494, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is S2XC3H6cC6H11 + O2 <=> C3H7cC6H10OO""",
)

entry(
    index = 2078,
    label = "C3H7TXcC6H10 + O2 <=> C3H7cC6H10OO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.17e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1494, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7TXcC6H10 + O2 <=> C3H7cC6H10OO""",
)

entry(
    index = 2079,
    label = "C3H7S2XcC6H10 + O2 <=> C3H7cC6H10OO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.17e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1494, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7S2XcC6H10 + O2 <=> C3H7cC6H10OO""",
)

entry(
    index = 2080,
    label = "C3H7S3XcC6H10 + O2 <=> C3H7cC6H10OO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.17e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1494, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7S3XcC6H10 + O2 <=> C3H7cC6H10OO""",
)

entry(
    index = 2081,
    label = "C3H7S4XcC6H10 + O2 <=> C3H7cC6H10OO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.17e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1494, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7S4XcC6H10 + O2 <=> C3H7cC6H10OO""",
)

entry(
    index = 2082,
    label = "PXC3H6cC6H11 + O2 <=> C3H7cC6H9OOH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.95e+13, 'cm^3/(mol*s)'), n=0, Ea=(12080, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC3H6cC6H11 + O2 <=> C3H7cC6H9OOH""",
)

entry(
    index = 2083,
    label = "SXC3H6cC6H11 + O2 <=> C3H7cC6H9OOH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.95e+13, 'cm^3/(mol*s)'), n=0, Ea=(12080, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC3H6cC6H11 + O2 <=> C3H7cC6H9OOH""",
)

entry(
    index = 2084,
    label = "S2XC3H6cC6H11 + O2 <=> C3H7cC6H9OOH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.95e+13, 'cm^3/(mol*s)'), n=0, Ea=(12080, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is S2XC3H6cC6H11 + O2 <=> C3H7cC6H9OOH""",
)

entry(
    index = 2085,
    label = "C3H7TXcC6H10 + O2 <=> C3H7cC6H9OOH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.95e+13, 'cm^3/(mol*s)'), n=0, Ea=(12080, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7TXcC6H10 + O2 <=> C3H7cC6H9OOH""",
)

entry(
    index = 2086,
    label = "C3H7S2XcC6H10 + O2 <=> C3H7cC6H9OOH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.95e+13, 'cm^3/(mol*s)'), n=0, Ea=(12080, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7S2XcC6H10 + O2 <=> C3H7cC6H9OOH""",
)

entry(
    index = 2087,
    label = "C3H7S3XcC6H10 + O2 <=> C3H7cC6H9OOH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.95e+13, 'cm^3/(mol*s)'), n=0, Ea=(12080, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7S3XcC6H10 + O2 <=> C3H7cC6H9OOH""",
)

entry(
    index = 2088,
    label = "C3H7S4XcC6H10 + O2 <=> C3H7cC6H9OOH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.95e+13, 'cm^3/(mol*s)'), n=0, Ea=(12080, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7S4XcC6H10 + O2 <=> C3H7cC6H9OOH""",
)

entry(
    index = 2089,
    label = "C3H7cC6H10OO <=> C3H7cC6H9OOH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.11e+09, 's^-1'), n=0, Ea=(20567, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7cC6H10OO <=> C3H7cC6H9OOH""",
)

entry(
    index = 2090,
    label = "C3H7cC6H10OO <=> C3H5cC6H11 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+12, 's^-1'), n=0, Ea=(33000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7cC6H10OO <=> C3H5cC6H11 + HO2""",
)

entry(
    index = 2091,
    label = "C3H7cC6H9OOH <=> C3H5cC6H11 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+11, 's^-1'), n=0, Ea=(13000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7cC6H9OOH <=> C3H5cC6H11 + HO2""",
)

entry(
    index = 2092,
    label = "C3H7cC6H9OOH => HCO + C5H10 + aC3H5 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.45e+12, 's^-1'), n=0, Ea=(18060, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7cC6H9OOH => HCO + C5H10 + aC3H5 + OH""",
)

entry(
    index = 2093,
    label = "C3H7cC6H9OOH + O2 => C3H7cC6H9O3 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(4e+12, 'cm^3/(mol*s)'), n=0, Ea=(10000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7cC6H9OOH + O2 => C3H7cC6H9O3 + OH""",
)

entry(
    index = 2094,
    label = "C3H7cC6H9O3 => OH + CH2CHO + CH2CO + C5H10",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(3e+13, 's^-1'), n=0, Ea=(20000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H7cC6H9O3 => OH + CH2CHO + CH2CO + C5H10""",
)

entry(
    index = 2095,
    label = "PXC4H8cC6H11 + HO2 => OH + C2H4 + CH3CHO + cC6H11",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC4H8cC6H11 + HO2 => OH + C2H4 + CH3CHO + cC6H11""",
)

entry(
    index = 2096,
    label = "SXC4H8cC6H11 + HO2 => OH + C2H4 + CH3CHO + cC6H11",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC4H8cC6H11 + HO2 => OH + C2H4 + CH3CHO + cC6H11""",
)

entry(
    index = 2097,
    label = "S2XC4H8cC6H11 + HO2 => OH + C2H4 + CH3CHO + cC6H11",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is S2XC4H8cC6H11 + HO2 => OH + C2H4 + CH3CHO + cC6H11""",
)

entry(
    index = 2098,
    label = "C4H9TXcC6H10 + HO2 => OH + C2H4 + CH3CHO + cC6H11",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9TXcC6H10 + HO2 => OH + C2H4 + CH3CHO + cC6H11""",
)

entry(
    index = 2099,
    label = "C4H9S2XcC6H10 + HO2 => OH + C2H4 + CH3CHO + cC6H11",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9S2XcC6H10 + HO2 => OH + C2H4 + CH3CHO + cC6H11""",
)

entry(
    index = 2100,
    label = "C4H9S3XcC6H10 + HO2 => OH + C2H4 + CH3CHO + cC6H11",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9S3XcC6H10 + HO2 => OH + C2H4 + CH3CHO + cC6H11""",
)

entry(
    index = 2101,
    label = "C4H9S4XcC6H10 + HO2 => OH + C2H4 + CH3CHO + cC6H11",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9S4XcC6H10 + HO2 => OH + C2H4 + CH3CHO + cC6H11""",
)

entry(
    index = 2102,
    label = "PXC4H8cC6H11 + O2 => HCO + C6H12 + aC3H5 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(7.19e+13, 'cm^3/(mol*s)'), n=0, Ea=(18056, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC4H8cC6H11 + O2 => HCO + C6H12 + aC3H5 + OH""",
)

entry(
    index = 2103,
    label = "SXC4H8cC6H11 + O2 => HCO + C6H12 + aC3H5 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(7.19e+13, 'cm^3/(mol*s)'), n=0, Ea=(18056, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC4H8cC6H11 + O2 => HCO + C6H12 + aC3H5 + OH""",
)

entry(
    index = 2104,
    label = "S2XC4H8cC6H11 + O2 => HCO + C6H12 + aC3H5 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(7.19e+13, 'cm^3/(mol*s)'), n=0, Ea=(18056, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is S2XC4H8cC6H11 + O2 => HCO + C6H12 + aC3H5 + OH""",
)

entry(
    index = 2105,
    label = "C4H9TXcC6H10 + O2 => HCO + C6H12 + aC3H5 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(7.19e+13, 'cm^3/(mol*s)'), n=0, Ea=(18056, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9TXcC6H10 + O2 => HCO + C6H12 + aC3H5 + OH""",
)

entry(
    index = 2106,
    label = "C4H9S2XcC6H10 + O2 => HCO + C6H12 + aC3H5 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(7.19e+13, 'cm^3/(mol*s)'), n=0, Ea=(18056, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9S2XcC6H10 + O2 => HCO + C6H12 + aC3H5 + OH""",
)

entry(
    index = 2107,
    label = "C4H9S3XcC6H10 + O2 => HCO + C6H12 + aC3H5 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(7.19e+13, 'cm^3/(mol*s)'), n=0, Ea=(18056, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9S3XcC6H10 + O2 => HCO + C6H12 + aC3H5 + OH""",
)

entry(
    index = 2108,
    label = "C4H9S4XcC6H10 + O2 => HCO + C6H12 + aC3H5 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(7.19e+13, 'cm^3/(mol*s)'), n=0, Ea=(18056, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9S4XcC6H10 + O2 => HCO + C6H12 + aC3H5 + OH""",
)

entry(
    index = 2109,
    label = "PXC4H8cC6H11 + O2 <=> C4H9cC6H10OO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.17e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1494, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC4H8cC6H11 + O2 <=> C4H9cC6H10OO""",
)

entry(
    index = 2110,
    label = "SXC4H8cC6H11 + O2 <=> C4H9cC6H10OO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.17e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1494, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC4H8cC6H11 + O2 <=> C4H9cC6H10OO""",
)

entry(
    index = 2111,
    label = "S2XC4H8cC6H11 + O2 <=> C4H9cC6H10OO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.17e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1494, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is S2XC4H8cC6H11 + O2 <=> C4H9cC6H10OO""",
)

entry(
    index = 2112,
    label = "C4H9TXcC6H10 + O2 <=> C4H9cC6H10OO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.17e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1494, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9TXcC6H10 + O2 <=> C4H9cC6H10OO""",
)

entry(
    index = 2113,
    label = "C4H9S2XcC6H10 + O2 <=> C4H9cC6H10OO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.17e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1494, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9S2XcC6H10 + O2 <=> C4H9cC6H10OO""",
)

entry(
    index = 2114,
    label = "C4H9S3XcC6H10 + O2 <=> C4H9cC6H10OO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.17e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1494, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9S3XcC6H10 + O2 <=> C4H9cC6H10OO""",
)

entry(
    index = 2115,
    label = "C4H9S4XcC6H10 + O2 <=> C4H9cC6H10OO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.17e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1494, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9S4XcC6H10 + O2 <=> C4H9cC6H10OO""",
)

entry(
    index = 2116,
    label = "PXC4H8cC6H11 + O2 <=> C4H9cC6H9OOH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.95e+13, 'cm^3/(mol*s)'), n=0, Ea=(12080, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PXC4H8cC6H11 + O2 <=> C4H9cC6H9OOH""",
)

entry(
    index = 2117,
    label = "SXC4H8cC6H11 + O2 <=> C4H9cC6H9OOH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.95e+13, 'cm^3/(mol*s)'), n=0, Ea=(12080, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SXC4H8cC6H11 + O2 <=> C4H9cC6H9OOH""",
)

entry(
    index = 2118,
    label = "S2XC4H8cC6H11 + O2 <=> C4H9cC6H9OOH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.95e+13, 'cm^3/(mol*s)'), n=0, Ea=(12080, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is S2XC4H8cC6H11 + O2 <=> C4H9cC6H9OOH""",
)

entry(
    index = 2119,
    label = "C4H9TXcC6H10 + O2 <=> C4H9cC6H9OOH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.95e+13, 'cm^3/(mol*s)'), n=0, Ea=(12080, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9TXcC6H10 + O2 <=> C4H9cC6H9OOH""",
)

entry(
    index = 2120,
    label = "C4H9S2XcC6H10 + O2 <=> C4H9cC6H9OOH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.95e+13, 'cm^3/(mol*s)'), n=0, Ea=(12080, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9S2XcC6H10 + O2 <=> C4H9cC6H9OOH""",
)

entry(
    index = 2121,
    label = "C4H9S3XcC6H10 + O2 <=> C4H9cC6H9OOH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.95e+13, 'cm^3/(mol*s)'), n=0, Ea=(12080, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9S3XcC6H10 + O2 <=> C4H9cC6H9OOH""",
)

entry(
    index = 2122,
    label = "C4H9S4XcC6H10 + O2 <=> C4H9cC6H9OOH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.95e+13, 'cm^3/(mol*s)'), n=0, Ea=(12080, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9S4XcC6H10 + O2 <=> C4H9cC6H9OOH""",
)

entry(
    index = 2123,
    label = "C4H9cC6H10OO <=> C4H9cC6H9OOH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.11e+09, 's^-1'), n=0, Ea=(20567, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H10OO <=> C4H9cC6H9OOH""",
)

entry(
    index = 2124,
    label = "C4H9cC6H10OO => C2H3cC6H11 + HO2 + C2H4",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 's^-1'), n=0, Ea=(33000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H10OO => C2H3cC6H11 + HO2 + C2H4""",
)

entry(
    index = 2125,
    label = "C4H9cC6H9OOH => C2H3cC6H11 + HO2 + C2H4",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(3e+11, 's^-1'), n=0, Ea=(13000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H9OOH => C2H3cC6H11 + HO2 + C2H4""",
)

entry(
    index = 2126,
    label = "C4H9cC6H9OOH => HCO + C6H12 + aC3H5 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.45e+12, 's^-1'), n=0, Ea=(18060, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H9OOH => HCO + C6H12 + aC3H5 + OH""",
)

entry(
    index = 2127,
    label = "C4H9cC6H9OOH + O2 => C4H9cC6H9O3 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(4e+12, 'cm^3/(mol*s)'), n=0, Ea=(10000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H9OOH + O2 => C4H9cC6H9O3 + OH""",
)

entry(
    index = 2128,
    label = "C4H9cC6H9O3 => OH + CH2CHO + CH2CO + C6H12",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(3e+13, 's^-1'), n=0, Ea=(20000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H9O3 => OH + CH2CHO + CH2CO + C6H12""",
)

entry(
    index = 2129,
    label = "CH3cC6H11 + O2 <=> CH2CHO + CH2O + pC4H9",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.3e+22, 'cm^3/(mol*s)'),
        n = -1.16,
        Ea = (72000, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3cC6H11 + O2 <=> CH2CHO + CH2O + pC4H9""",
)

entry(
    index = 2130,
    label = "C2H5cC6H11 + O2 <=> CH2CHO + CH2O + C2H4 + C2H4 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.5e+24, 'cm^3/(mol*s)'),
        n = -1.77,
        Ea = (76000, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H5cC6H11 + O2 <=> CH2CHO + CH2O + C2H4 + C2H4 + CH3""",
)

entry(
    index = 2131,
    label = "C3H7cC6H11 + O2 <=> CH2CHO + CH2O + C2H4 + C2H4 + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.5e+24, 'cm^3/(mol*s)'),
        n = -1.77,
        Ea = (76000, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C3H7cC6H11 + O2 <=> CH2CHO + CH2O + C2H4 + C2H4 + C2H5""",
)

entry(
    index = 2132,
    label = "C4H9cC6H11 + O2 <=> CH2CHO + CH2O + C2H4 + C2H4 + C2H4 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.5e+24, 'cm^3/(mol*s)'),
        n = -1.77,
        Ea = (76000, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C4H9cC6H11 + O2 <=> CH2CHO + CH2O + C2H4 + C2H4 + C2H4 + CH3""",
)

entry(
    index = 2133,
    label = "CH* => CH",
    degeneracy = 1,
    duplicate = True,
    reversible = False,
    kinetics = Arrhenius(A=(1.9e+06, 's^-1'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH* => CH""",
)

entry(
    index = 2134,
    label = "CH* => CH",
    degeneracy = 1,
    duplicate = True,
    reversible = False,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(A=(4e+10, 'cm^3/(mol*s)'), n=0.5, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is CH* => CH""",
)

entry(
    index = 2135,
    label = "CH* + O2 => CH + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.4e+12, 'cm^3/(mol*s)'), n=0.5, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH* + O2 => CH + O2""",
)

entry(
    index = 2136,
    label = "C2H + O2 => CH* + CO2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(4.5e+15, 'cm^3/(mol*s)'), n=0, Ea=(25000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H + O2 => CH* + CO2""",
)

entry(
    index = 2137,
    label = "H + O <=> OH*",
    degeneracy = 1,
    duplicate = True,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(
            A = (3.1e+14, 'cm^6/(mol^2*s)'),
            n = 0,
            Ea = (10000, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is H + O <=> OH*""",
)

entry(
    index = 2138,
    label = "OH* + AR <=> OH + AR",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.17e+10, 'cm^3/(mol*s)'),
        n = 0.5,
        Ea = (2060, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is OH* + AR <=> OH + AR""",
)

entry(
    index = 2139,
    label = "OH* + H2O <=> OH + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.92e+12, 'cm^3/(mol*s)'),
        n = 0.5,
        Ea = (-861, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is OH* + H2O <=> OH + H2O""",
)

entry(
    index = 2140,
    label = "OH* + CO2 <=> OH + CO2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.75e+12, 'cm^3/(mol*s)'),
        n = 0.5,
        Ea = (-968, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is OH* + CO2 <=> OH + CO2""",
)

entry(
    index = 2141,
    label = "OH* + CO <=> OH + CO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.23e+12, 'cm^3/(mol*s)'),
        n = 0.5,
        Ea = (-787, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is OH* + CO <=> OH + CO""",
)

entry(
    index = 2142,
    label = "OH* + H2 <=> OH + H2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.95e+12, 'cm^3/(mol*s)'),
        n = 0.5,
        Ea = (-444, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is OH* + H2 <=> OH + H2""",
)

entry(
    index = 2143,
    label = "OH* + O2 <=> OH + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.1e+12, 'cm^3/(mol*s)'), n=0.5, Ea=(-482, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is OH* + O2 <=> OH + O2""",
)

entry(
    index = 2144,
    label = "OH* + OH <=> OH + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.5e+12, 'cm^3/(mol*s)'), n=0.5, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is OH* + OH <=> OH + OH""",
)

entry(
    index = 2145,
    label = "OH* + H <=> OH + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.5e+12, 'cm^3/(mol*s)'), n=0.5, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is OH* + H <=> OH + H""",
)

entry(
    index = 2146,
    label = "OH* + O <=> OH + O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.5e+12, 'cm^3/(mol*s)'), n=0.5, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is OH* + O <=> OH + O""",
)

entry(
    index = 2147,
    label = "OH* + CH4 <=> OH + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.36e+12, 'cm^3/(mol*s)'),
        n = 0.5,
        Ea = (-635, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is OH* + CH4 <=> OH + CH4""",
)

entry(
    index = 2148,
    label = "OH* <=> OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.4e+06, 's^-1'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is OH* <=> OH""",
)

