#!/usr/bin/env python
# encoding: utf-8

name = "Surface_Dissociation/rules"
shortDesc = u""
longDesc = u"""
"""
entry(
    index = 1,
    label = "Combined;VacantSite",
    kinetics = SurfaceArrheniusBEP(
        A = (1.0e17, 'm^2/(mol*s)'),
        n = 0,
        alpha = 0.84,
        E0 = (44.25, 'kcal/mol'),
        Tmin = (200, 'K'),
        Tmax = (3000, 'K'),
    ),
    rank = 0,
    shortDesc = u"""Default""",
    longDesc = u"""
"Arrhenius preexponential values for surface recombination...reactions
are, in the SI system,... 10^13 - 10^14 m2/mol/s ...for bimolecular reactions"
from page 54 of "Silicon epitaxy"
Author:	Danilo Crippa; Daniel L Rode; Maurizio Masi
Publisher:	San Diego : Academic Press, 2001.
Series:	Semiconductors and semimetals, v. 72.

E0 and alpha are taken from:
"Universal Brønsted-Evans-Polanyi Relations for C–C, C–O, C–N, N–O, N–N, and O–O Dissociation Reactions" by Wang, ..., Norskov/ Catal. Lett (2011) 141:370-373.
DOI 10.1007/s10562-010-0477-y
(actual value for E0 was 1.92 eV.)
    """
)



