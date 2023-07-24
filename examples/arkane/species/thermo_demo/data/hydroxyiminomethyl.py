#!/usr/bin/env python3
# encoding: utf-8

bonds = {'C-H': 1, 'H-O': 1, 'C=N': 1, 'N-O': 1}

externalSymmetry = 1

spinMultiplicity = 2

opticalIsomers = 1

energy = Log('hydroxyiminomethyl/sp.out')

geometry = Log('hydroxyiminomethyl/freq.out')

frequencies = Log('hydroxyiminomethyl/freq.out')

rotors = [HinderedRotor(scanLog=Log('hydroxyiminomethyl/scan_2_3.out'), pivots=[2, 3], top=[3, 5], symmetry=1, fit='fourier')]

