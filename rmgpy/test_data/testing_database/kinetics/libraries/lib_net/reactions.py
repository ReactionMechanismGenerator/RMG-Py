#!/usr/bin/env python
# encoding: utf-8

name = "lib_net"
shortDesc = u""
longDesc = u"""
This is a fake library for testing PDep networks exploration of library reactions
"""

entry(
    index = 1,
    label = "CH3 + CH3 <=> ethane",
    degeneracy = 1.0,
    elementary_high_p = True,
    kinetics = Chebyshev(
        coeffs = [
            [12.8791, 0.560009, -0.115177, 0.0105093],
            [-1.04561, 0.857712, -0.116582, -0.0134747],
            [-0.599769, 0.449793, 0.0096248, -0.0271013],
            [-0.314421, 0.177204, 0.0458462, -0.0107188],
            [-0.155503, 0.0517415, 0.0323201, 0.00343069],
            [-0.0727158, 0.00622718, 0.0138549, 0.00630263],
        ],
        kunits = 'cm^3/(mol*s)',
        Tmin = (300, 'K'),
        Tmax = (3000, 'K'),
        Pmin = (0.001, 'bar'),
        Pmax = (100, 'bar'),
    ),
)

entry(
    index = 2,
    label = "H + CH3 <=> CH4",
    degeneracy = 1.0,
    elementary_high_p = True,
    kinetics = Arrhenius(
        A = (2.94005e-11, 'm^3/(mol*s)'),
        n = 5.135,
        Ea = (33.0118, 'kJ/mol'),
        T0 = (1, 'K')),
)

entry(
    index = 3,
    label = "O2 + H2 <=> H + [O]O",
    degeneracy = 4.0,
    kinetics = Arrhenius(
        A = (2.9e+14, 'cm^3/(mol*s)', '*|/', 5),
        n = 0,
        Ea = (236.982, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (800, 'K'),
    ),
)

entry(
    index = 4,
    label = "CO + OH <=> CO2 + H",
    degeneracy = 1,
    kinetics = PDepArrhenius(
        pressures = ([0.01315, 0.1315, 1.315, 13.158, 131.58], 'atm'),
        arrhenius = [
            Arrhenius(A=(210000, 'cm^3/(mol*s)'), n=1.9, Ea=(-1064, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(250000, 'cm^3/(mol*s)'), n=1.88, Ea=(-1043, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(870000, 'cm^3/(mol*s)'), n=1.73, Ea=(-685, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(6.8e+06, 'cm^3/(mol*s)'), n=1.48, Ea=(48, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2.3e+07, 'cm^3/(mol*s)'), n=1.35, Ea=(974, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CO + OH <=> CO2 + H""",
)

entry(
    index = 5,
    label = "C2H2 + CH3 <=> C3H5",
    degeneracy = 1,
    elementary_high_p = True,
    kinetics = PDepArrhenius(
        pressures = ([1, 2, 5], 'atm'),
        arrhenius = [
            Arrhenius(A=(4.99e+22, 'cm^3/(mol*s)'), n=-4.39, Ea = (18850, 'cal/mol'), T0 = (1, 'K')),
            Arrhenius(A=(6e+23, 'cm^3/(mol*s)'), n=-4.6, Ea=(19571, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(7.31e+25, 'cm^3/(mol*s)'), n=-5.06, Ea=(21150, 'cal/mol'), T0 = (1, 'K')),
        ],
    ),
)

entry(
    index = 6,
    label = "H + CH2 <=> CH3",
    degeneracy = 1,
    elementary_high_p = True,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(6e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.04e+26, 'cm^6/(mol^2*s)'),
            n = -2.76,
            Ea = (1600, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.562,
        T3 = (91, 'K'),
        T1 = (5836, 'K'),
        T2 = (8552, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C]=O': 1.5, '[Ar]': 0.7},
    ),
)

entry(
    index = 7,
    label = "O + CO <=> CO2",
    degeneracy = 1,
    elementary_high_p = True,
    kinetics = Lindemann(
        arrheniusHigh = Arrhenius(A=(1.8e+10, 'cm^3/(mol*s)'), n=0, Ea=(2385, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (6.02e+14, 'cm^6/(mol^2*s)'),
            n = 0,
            Ea = (3000, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        efficiencies = {'C': 2, 'O=C=O': 3.5, 'CC': 3, 'O': 6, '[H][H]': 2, '[O][O]': 6, '[C]=O': 1.5, '[Ar]': 0.5},
    ),
)
