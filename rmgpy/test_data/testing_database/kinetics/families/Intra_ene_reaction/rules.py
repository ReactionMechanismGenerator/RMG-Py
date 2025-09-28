#!/usr/bin/env python
# encoding: utf-8

name = "Intra_ene_reaction/rules"
shortDesc = u"Intramolecular H-shift from an allylic to an unsaturated endgroup (like in cyclopentadiene)"
longDesc = u"""

"""
entry(
    index = 1,
    label = "cyclopentadiene;CH_end;unsaturated_end",
    kinetics = ArrheniusEP(
        A = (5.06e+07, 's^-1'),
        n = 1.74,
        alpha = 0,
        E0 = (24.3, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (1500, 'K'),
    ),
    rank = 4,
    shortDesc = u"""AG Vandeputte, CBS-QB3""",
    longDesc = """Rate taken from H shift in ethyleneCPD""",
)

