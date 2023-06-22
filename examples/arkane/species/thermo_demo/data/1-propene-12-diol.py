#!/usr/bin/env python3
# encoding: utf-8

bonds = {'C-C': 1, 'C-O': 2, 'C-H': 4, 'H-O': 2, 'C=C': 1}

externalSymmetry = 1

spinMultiplicity = 1

opticalIsomers = 2

energy = Log('1-propene-12-diol/sp.out')

geometry = Log('1-propene-12-diol/freq.out')

frequencies = Log('1-propene-12-diol/freq.out')

rotors = [HinderedRotor(scanLog=Log('1-propene-12-diol/scan_1_2.out'), pivots=[1, 2], top=[1, 6, 7, 8], symmetry=3, fit='fourier'),
          HinderedRotor(scanLog=Log('1-propene-12-diol/scan_2_3.out'), pivots=[2, 3], top=[3, 9], symmetry=1, fit='fourier'),
          HinderedRotor(scanLog=Log('1-propene-12-diol/scan_4_5.out'), pivots=[4, 5], top=[5, 11], symmetry=1, fit='fourier'),
]

