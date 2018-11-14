#!/usr/bin/env python
# -*- coding: utf-8 -*-

bonds = {
    'C=C': 1, 
    'C-H': 4,
}

linear = False

externalSymmetry = 4

spinMultiplicity = 1

opticalIsomers = 1

energy = {
    'CBS-QB3': Log('ethene.log'),
    'Klip_2': -78.42735579,
    'CCSD(T)-F12/cc-pVTZ-F12': Log('ethene_f12.out'),
}

geometry = Log('ethene.log')

frequencies = Log('ethene_freq.log')

rotors = []
