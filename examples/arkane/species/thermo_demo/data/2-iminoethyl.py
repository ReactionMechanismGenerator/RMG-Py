#!/usr/bin/env python3
# encoding: utf-8

bonds = {'C-C': 1, 'C=N': 1, 'C-H': 3, 'N-H': 1}

externalSymmetry = 1

spinMultiplicity = 2

opticalIsomers = 1

energy = Log('2-iminoethyl/sp.out')

geometry = Log('2-iminoethyl/freq.out')

frequencies = Log('2-iminoethyl/freq.out')

rotors = [HinderedRotor(scanLog=Log('2-iminoethyl/scan_1_2.out'), pivots=[1, 2], top=[1, 4, 5], symmetry=2, fit='fourier'),
          HinderedRotor(scanLog=Log('2-iminoethyl/scan_2_3.out'), pivots=[2, 3], top=[3, 7], symmetry=1, fit='fourier')]

