#!/usr/bin/env python
# encoding: utf-8

"""
This example loads all species data from the C2H6.yml file in this folder.
The C2H4 and C2H6 examples demonstare how to generate such .yml files.
"""

modelChemistry = "CBS-QB3"
useHinderedRotors = True
useBondCorrections = False

species('C2H6', 'C2H6.yml')

thermo('C2H6', 'NASA')
