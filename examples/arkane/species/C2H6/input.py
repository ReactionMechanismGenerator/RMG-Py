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

species('C2H6', 'C2H6.py',
        structure=SMILES('CC'))

thermo('C2H6', 'NASA')
