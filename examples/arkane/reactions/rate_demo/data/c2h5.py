#!/usr/bin/env python3
# encoding: utf-8

externalSymmetry = 1

spinMultiplicity = 2

opticalIsomers = 1

energy = Log('c2h5/sp.out')

geometry = Log('c2h5/freq.out')

frequencies = Log('c2h5/freq.out')

rotors = [FreeRotor(pivots=[1, 2], top=[1, 3, 4], symmetry=6)]

