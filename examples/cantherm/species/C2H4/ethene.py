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
    'CBS-QB3': GaussianLog('ethene.log'),
    'Klip_2': -78.42735579,
    'CCSD(T)-F12/cc-pVTZ-F12': MolproLog('ethene_f12.out'),
}

geometry = GaussianLog('ethene.log')

frequencies = GaussianLog('ethene_freq.log')

rotors = []
