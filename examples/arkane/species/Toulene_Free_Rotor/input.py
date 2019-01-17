#!/usr/bin/env python
# encoding: utf-8

modelChemistry = 'CBS-QB3'
useHinderedRotors = True
useBondCorrections = True

species('Toluene', 'toluene_FreeRotor.py')

statmech('Toluene')
thermo('Toluene', 'NASA')
