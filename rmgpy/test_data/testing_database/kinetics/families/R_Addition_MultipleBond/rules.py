#!/usr/bin/env python
# encoding: utf-8

name = "R_Addition_MultipleBond/rules"
shortDesc = u""
longDesc = u"""

"""
entry(
    index = 1,
    label = "Cd_R;CsJ",
    kinetics = ArrheniusEP(
        A = (20900, 'cm^3/(mol*s)'),
        n = 2.41,
        alpha = 0,
        E0 = (5.63, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (1500, 'K'),
    ),
    rank = 4,
    shortDesc = u"""Original label: Cds-HH_Cds-HH;CsJ-HHH""",
)

