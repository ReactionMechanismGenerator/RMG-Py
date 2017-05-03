#!/usr/bin/env python
# encoding: utf-8

name = "intra_H_migration/rules"
shortDesc = u""
longDesc = u"""

"""
entry(
    index = 614,
    label = "RnH;Y_rad_out;XH_out",
    kinetics = ArrheniusEP(
        A = (1e+10, 's^-1'),
        n = 0,
        alpha = 0,
        E0 = (25, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (1500, 'K'),
    ),
    rank = 0,
    shortDesc = u"""default""",
)

entry(
    index = 615,
    label = "R6H;C_rad_out_single;Cs_H_out",
    kinetics = ArrheniusEP(
        A = (5.48e+08, 's^-1'),
        n = 1.62,
        alpha = 0,
        E0 = (38.76, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (1500, 'K'),
    ),
    rank = 5,
    shortDesc = u"""Currans's estimation [5] in his reaction type 5.""",
    longDesc = 
u"""
Test case with modified group from database. This data is definitly not
scientifically accurate, so do not use it in model generation!!!"
""",
)

