#!/usr/bin/env python
# encoding: utf-8

name = "Singlet_Val6_to_triplet/rules"
shortDesc = ""
longDesc = """
"""

entry(
    index = 1,
    label = "singlet",
    kinetics = Arrhenius(A=(4.5E+10, 's^-1'), n=0, Ea=(397, 'cal/mol')),
    rank = 1,
    shortDesc = """Default""",
    longDesc =
"""
taken from:
R. Atkinson, D.L. Baulch, R.A. Cox, R.F. Hampson, J.A. Kerr, J. Troe,
Evaluated Kinetic and Photochemical Data for Atmospheric Chemistry: Supplement IV.
IUPAC Subcommittee on Gas Kinetic Data Evaluation for Atmospheric Chemistry
Journal of Physical and Chemical Reference Data 21, 1125 (1992)
doi: 10.1063/1.555918

Adjusted to a first order reaction at 1 atm by alongd:
n/V = P/RT = 1 bar / (83 cm^3 bar K^-1 mol^-1 * 300 K) = 4E-05 mol cm^-3
1.81E+06 mol cm^-3 S^-1 / 4E-05 mol cm^-3 = 4.5E+10 s^-1
Original reaction is O2(1D) + M => O2 + M
""",
)
