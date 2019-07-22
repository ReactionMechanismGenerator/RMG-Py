#!/usr/bin/env python
# encoding: utf-8

"""
This input file generates a  .yml file in the ArkaneSpecies folder
(If a `thermo()` function exists and the species `structure` is specified,
a respective .yml file is automatically generated)
"""

modelChemistry = "CBS-QB3"
useHinderedRotors = True
useBondCorrections = False
author = 'I.B. Modeling'
# Example of how to manually set the frequencyScaleFactor
# Factor of 0.99 for CBS-QB3 ZPE scaling (dx.doi.org/10.1063/1.477924)
# Factor of 1.014 for converting between ZPE and frequency scaling factors (dx.doi.org/10.1021/ct100326h)
frequencyScaleFactor = 0.99*1.014

species('C2H6', 'C2H6.py',
        structure=SMILES('CC'))

thermo('C2H6', 'NASA')
