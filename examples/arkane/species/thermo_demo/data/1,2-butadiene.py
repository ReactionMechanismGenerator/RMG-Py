#!/usr/bin/env python3
# encoding: utf-8

bonds = {'C-C': 1, 'C=C': 2, 'C-H': 6}

externalSymmetry = 1

spinMultiplicity = 1

opticalIsomers = 1

energy = Log('1,2-butadiene/sp.out')

geometry = Log('1,2-butadiene/freq.out')

frequencies = Log('1,2-butadiene/freq.out')

rotors = [HinderedRotor(scanLog=Log('1,2-butadiene/scan_1_2.out'), pivots=[1, 2], top=[1, 5, 6, 7], symmetry=3, fit='fourier')]

