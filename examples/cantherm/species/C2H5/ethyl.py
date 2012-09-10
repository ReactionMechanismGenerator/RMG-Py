#!/usr/bin/env python
# -*- coding: utf-8 -*-

atoms = {
    'C': 2,
    'H': 5,
}

bonds = {
    'C-C': 1, 
    'C-H': 5,
}

linear = False

externalSymmetry = 1

spinMultiplicity = 2

opticalIsomers = 1

energy = {
    'CBS-QB3': GaussianLog('ethyl_cbsqb3.log'),
    'Klip_2': -78.98344186,
}

geometry = GaussianLog('ethyl_cbsqb3.log')

frequencies = GaussianLog('ethyl_b3lyp.log')

rotors = [
    HinderedRotor(scanLog=GaussianLog('ethyl_scan_72.log'), pivots=[1,2], top=[1,3,4], symmetry=6)
]
