#!/usr/bin/env python
# encoding: utf-8

name = "Surface_Adsorption_Dissociative/rules"
shortDesc = u""
longDesc = u"""
Dissociative adsorption of a gas-phase species forming two adsorbates, each with a single bond to a surface site
"""
entry(
    index = 1,
    label = "Adsorbate;VacantSite1;VacantSite2",
    kinetics = StickingCoefficientBEP(
        A = 0.01,
        n = 0,
        alpha = 0,
        E0 = (10, 'kcal/mol'),
        Tmin = (200, 'K'),
        Tmax = (3000, 'K'),
    ),
    rank = 0,
    shortDesc = u"""Default""",
    longDesc = u"""Made up"""
)



