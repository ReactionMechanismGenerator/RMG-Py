#!/usr/bin/env python
# encoding: utf-8

name = "custom_library"
shortDesc = u"dummy external kinetic library"
longDesc = u"""
This is a dummy kinetic library,
showcasing how to use a kinetic library that is external to the RMG-database.
"""

entry(
    index = 1,
    label = "H + O2 <=> O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (9.841e+13, 'cm^3/(mol*s)'),
        n = 0,
        Ea = (15310, 'cal/mol'),
        T0 = (1, 'K'),
    ),
)

entry(
    index = 2,
    label = "CH2(T) + CH <=> H + C2H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
)

