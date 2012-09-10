#!/usr/bin/env python
# encoding: utf-8

modelChemistry = "Klip_2"
frequencyScaleFactor = 0.99
useHinderedRotors = True
useBondCorrections = False

species('C2H6', 'C2H6.py')

statmech('C2H6')
thermo('C2H6', 'NASA')
