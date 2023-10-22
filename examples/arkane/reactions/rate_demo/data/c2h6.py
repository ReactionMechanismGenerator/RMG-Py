#!/usr/bin/env python3
# encoding: utf-8

externalSymmetry = 6

spinMultiplicity = 1

opticalIsomers = 1

energy = Log('c2h6/sp.out')

geometry = Log('c2h6/freq.out')

frequencies = Log('c2h6/freq.out')

rotors = [HinderedRotor(scanLog=Log('c2h6/scan_1_2.out'), pivots=[1, 2], top=[1, 3, 4, 5], symmetry=3, fit='fourier')]

