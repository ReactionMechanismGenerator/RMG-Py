#!/usr/bin/env python
# encoding: utf-8

modelChemistry = "CBS-QB3"
frequencyScaleFactor = 0.99
useHinderedRotors = True
useBondCorrections = True

species('Toluene', 'toluene_HinderedRotor.py')

statmech('Toluene')
thermo('Toluene', 'NASA')
