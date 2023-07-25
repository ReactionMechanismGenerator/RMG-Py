#!/usr/bin/env python
# encoding: utf-8

name = "kinetics"
shortDesc = ""
longDesc = """
Calculated using Arkane v3.1.0 using LevelOfTheory(method='cbsqb3').
"""
autoGenerated=True
entry(
    index = 0,
    label = "H + C2H4 <=> C2H5",
    degeneracy = 1.0,
    elementary_high_p = True,
    kinetics = Arrhenius(A=(4.00615e+08,'cm^3/(mol*s)'), n=1.66398, Ea=(5.95613,'kJ/mol'), T0=(1,'K'), Tmin=(400,'K'), Tmax=(1200,'K'), comment="""Fitted to 6 data points; dA = *|/ 1.0977, dn = +|- 0.0123114, dEa = +|- 0.0678956 kJ/mol"""),
)
