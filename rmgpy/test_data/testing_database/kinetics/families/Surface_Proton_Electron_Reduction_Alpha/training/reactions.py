#!/usr/bin/env python
# encoding: utf-8

name = "Surface_Proton_Electron_Reduction_Alpha/training"
shortDesc = u"Reaction kinetics used to generate rate rules"
longDesc = u"""
Put kinetic parameters for specific reactions in this file to use as a
training set for generating rate rules to populate this kinetics family.
"""

# entry(
#     index = 1,
#     label = "OX + H + e <=> HOX",
#     degeneracy = 1,
#     kinetics = SurfaceChargeTransfer(
#         a = 0.5, # charge transfer coeff
#         A = (1.0e10, 'm^3/(mol*s)'), # pre-exponential factor
#         n = 0, # temperature coeff
#         V0 = (0, 'V') # reference potential
#         Ea = (40, 'kJ/mol') # activation energy
#         Tmin = (200, 'K'),
#         Tmax = (3000, 'K'),
#         ne = -1, # electron stochiometric coeff 
#     ),
#     shortDesc = u"""Default""",
#     longDesc = u"""
# David made it up as example
# """,
#     metal = "Pt",
#     facet = "111",
#     site = ""
# )
