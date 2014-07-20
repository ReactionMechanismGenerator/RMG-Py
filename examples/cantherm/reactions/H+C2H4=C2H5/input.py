#!/usr/bin/env python
# -*- coding: utf-8 -*-

modelChemistry = "CBS-QB3"
frequencyScaleFactor = 0.99
useHinderedRotors = False
useBondCorrections = True

species('H', '../../species/H/H.py')
species('C2H4', '../../species/C2H4/ethene.py')
species('C2H5', '../../species/C2H5/ethyl.py')
transitionState('TS', 'TS.py')

reaction(
    label = 'H + C2H4 <=> C2H5',
    reactants = ['H', 'C2H4'],
    products = ['C2H5'],
    transitionState = 'TS',
    tunneling='Eckart',
)

statmech('TS')
kinetics('H + C2H4 <=> C2H5')
