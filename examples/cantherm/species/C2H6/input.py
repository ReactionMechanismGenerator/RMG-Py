#!/usr/bin/env python
# encoding: utf-8

modelChemistry = "CBS-QB3"
useHinderedRotors = True
useBondCorrections = False

species('C2H6', 'C2H6.py')

statmech('C2H6')
thermo('C2H6', 'NASA')
