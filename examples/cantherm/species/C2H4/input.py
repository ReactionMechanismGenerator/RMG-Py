#!/usr/bin/env python
# encoding: utf-8

modelChemistry = "Klip_2"
frequencyScaleFactor = 0.99
useHinderedRotors = True
useBondCorrections = False

species('C2H4', 'ethene.py')

statmech('C2H4')
thermo('C2H4', 'NASA')
