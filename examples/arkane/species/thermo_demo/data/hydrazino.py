#!/usr/bin/env python3
# encoding: utf-8

bonds = {'N-N': 1, 'H-N': 3}

externalSymmetry = 1

spinMultiplicity = 2

opticalIsomers = 1

energy = Log('hydrazino/sp.out')

geometry = Log('hydrazino/freq.out')

frequencies = Log('hydrazino/freq.out')

rotors = [HinderedRotor(scanLog=Log('hydrazino/scan_1_2.out'), pivots=[1, 2], top=[2, 5], symmetry=1, fit='fourier')]

