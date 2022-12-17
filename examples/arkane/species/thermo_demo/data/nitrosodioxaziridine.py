#!/usr/bin/env python3
# encoding: utf-8

bonds = {'O-O': 1, 'N-O': 2, 'N-N': 1, 'N=O': 1}

externalSymmetry = 1

spinMultiplicity = 1

opticalIsomers = 1

energy = Log('nitrosodioxaziridine/sp.out')

geometry = Log('nitrosodioxaziridine/freq.out')

frequencies = Log('nitrosodioxaziridine/freq.out')

rotors = [HinderedRotor(scanLog=Log('nitrosodioxaziridine/scan_2_3.out'),
          pivots=[2, 3], top=[2, 1], symmetry=1, fit='fourier')]

