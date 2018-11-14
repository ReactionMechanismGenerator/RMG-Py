#!/usr/bin/env python
# -*- coding: utf-8 -*-

bonds = {
    'C-C': 4, 
    'C-H': 8,
    'C=C': 3,
}

linear = False

externalSymmetry = 1

spinMultiplicity = 1

opticalIsomers = 1

energy = {
    'CBS-QB3': Log('TolueneEnergy.log')
}

geometry = Log('TolueneFreq.log')

frequencies = Log('TolueneFreq.log')

rotors = [
HinderedRotor(scanLog=Log('TolueneRot1.log'),pivots=[3,12],top=[12,13,14,15],symmetry=6, fit='best'),
]
