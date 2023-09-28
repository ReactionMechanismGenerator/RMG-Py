#!/usr/bin/env python3
# encoding: utf-8

externalSymmetry = 1

spinMultiplicity = 2

opticalIsomers = 1

energy = Log('ts2/sp.out')

geometry = Log('ts2/freq.out')

frequencies = Log('ts2/freq.out')

rotors = [HinderedRotor(scanLog=Log('ts2/scan_4_5.out'), pivots=[4, 5], top=[4, 6, 7, 8], symmetry=3, fit='fourier')]

