#!/usr/bin/env python
# encoding: utf-8

name = "Surface"
shortDesc = u""
longDesc = u"""
test surface mechanism: based upon Olaf Deutschmann's work:
"Surface Reaction Kinetics of Steam- and CO2-Reforming as well as Oxidation of Methane over Nickel-Based Catalysts"
Delgado et al
Catalysts, 2015, 5, 871-904
"""

entry(
    index = 1,
    label = "O2 + Ni + Ni <=> OX + OX",
    kinetics = StickingCoefficient(
        A = 4.36E-2,
        n = -0.206,
        Ea=(1.5E3, 'J/mol'),
        Tmin = (200, 'K'),
        Tmax = (3000, 'K'),
    ),
    shortDesc = u"""Default""",
    longDesc = u"""Made up""",
	metal = "Ni",
)


entry(
    index = 2,
    label = "O2 + Ni <=> O2X",
    kinetics = StickingCoefficient(
        A = 0.0,
        n = 0,
        Ea=(0, 'J/mol'),
        Tmin = (200, 'K'),
        Tmax = (3000, 'K'),
    ),
    shortDesc = u"""Default""",
    longDesc = u"""Made up""",
	metal = "Ni",
)

entry(
    index = 3,
    label = "CH4 + Ni + Ni <=> CH3X + HX",
    kinetics = StickingCoefficient(
        A = 8.0E-3,
        n = 0,
        Ea=(0, 'J/mol'),
        Tmin = (200, 'K'),
        Tmax = (3000, 'K'),
    ),
    shortDesc = u"""Default""",
    longDesc = u"""Made up""",
	metal = "Ni",
)

entry(
    index = 4,
    label = "H2 + Ni + Ni <=> HX + HX",
    kinetics = StickingCoefficient(
        A = 3.2E-2,
        n = 0,
        Ea=(0, 'J/mol'),
        Tmin = (200, 'K'),
        Tmax = (3000, 'K'),
    ),
    shortDesc = u"""Default""",
    longDesc = u"""Made up""",
	metal = "Ni",
)


entry(
    index = 5,
    label = "HX + HOX <=> H2O + Ni + Ni",
    kinetics = SurfaceArrhenius(
        A=(1.85E16, 'm^2/(mol*s)'),
        n = 0.086,
        Ea=(41500.0, 'J/mol'),
        Tmin = (200, 'K'),
        Tmax = (3000, 'K'),
    ),
    shortDesc = u"""Default""",
    longDesc = u"""Made up""",
	metal = "Ni",
)


entry(
    index = 6,
    label = "CO + Ni <=> OCX",
    kinetics = StickingCoefficient(
        A = 5.0E-1,
        n = 0,
        Ea=(0, 'J/mol'),
        Tmin = (200, 'K'),
        Tmax = (3000, 'K'),
    ),
    shortDesc = u"""Default""",
    longDesc = u"""Made up""",
	metal = "Ni",
)



entry(
    index = 7,
    label = "OCX + OX <=> CO2 + Ni + Ni",
    kinetics = SurfaceArrhenius(
        A=(2.00E15, 'm^2/(mol*s)'),
        n = 0.0,
        Ea=(123600.0, 'J/mol'),
        Tmin = (200, 'K'),
        Tmax = (3000, 'K'),
    ),
    shortDesc = u"""Default""",
    longDesc = u"""Made up""",
	metal = "Ni",
)
