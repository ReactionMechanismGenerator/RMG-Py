#!/usr/bin/env python3
# encoding: utf-8

bonds = {'C-C': 3, 'H-N': 2, 'C-H': 9, 'C-N': 1}

externalSymmetry = 1

spinMultiplicity = 1

opticalIsomers = 1

energy = Log('2-methyl-2-propanamine/sp.out')

geometry = Log('2-methyl-2-propanamine/freq.out')

frequencies = Log('2-methyl-2-propanamine/freq.out')

rotors = [HinderedRotor(scanLog=Log('2-methyl-2-propanamine/scan_1_2.out'), pivots=[1, 2], top=[1, 6, 7, 8], symmetry=3, fit='fourier'),
          HinderedRotor(scanLog=Log('2-methyl-2-propanamine/scan_2_3.out'), pivots=[2, 3], top=[3, 9, 10, 11], symmetry=3, fit='fourier'),
          HinderedRotor(scanLog=Log('2-methyl-2-propanamine/scan_2_4.out'), pivots=[2, 4], top=[4, 12, 13, 14], symmetry=3, fit='fourier'),
          HinderedRotor(scanLog=Log('2-methyl-2-propanamine/scan_2_5.out'), pivots=[2, 5], top=[5, 15, 16], symmetry=3, fit='fourier')]

