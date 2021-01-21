#!/usr/bin/env python
# encoding: utf-8

name = "Surface_Proton_Electron_Reduction_Alpha/rules"
shortDesc = u""
longDesc = u"""
Surface adsorption of a single radical forming a single bond to the surface site
"""

entry(
    index = 1,
    label = "Adsorbate;Proton;Electron",
    kinetics = SurfaceChargeTransfer(
        A = (2.483E21, 'cm^3/(mol*s)'), # pre-exponential factor 1E9/s / surface_site_density / [H+]0 (1M/L)
        n = 0, # temperature coeff, 0 default
        V0 = None, # Reversible potential (determined based on free energy of reaction)
        Ea = (1, 'kJ/mol'), # activation energy at the reversible potential
        Tmin = (200, 'K'),
        Tmax = (3000, 'K'),
        ne = -1, # electron stochiometric coeff
    ),
    rank = 0,
    shortDesc = u"""Default""",
    longDesc = u"""https://doi.org/10.1021/jp4100608"""
)

