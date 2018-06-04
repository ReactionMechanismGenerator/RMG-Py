#!/usr/bin/env python
# -*- coding: utf-8 -*-

bonds = {
	'C=C': 3,
    'C-C': 4, 
    'C-H': 7,
}

linear = False

externalSymmetry = 2

spinMultiplicity = 2

opticalIsomers = 1

energy = {
    'CBS-QB3': Log('BenzylEnergy.log'),
}

geometry = Log('BenzylFreq.log')

frequencies = Log('BenzylFreq.log')

rotors = [ 
    HinderedRotor(scanLog=Log('BenzylRot1.log'), pivots=[12,4], top=[12,13,14], symmetry=2, fit='best'),
    ]
