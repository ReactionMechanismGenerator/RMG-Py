#!/usr/bin/env python
# -*- coding: utf-8 -*-

atoms = {
    'C': 7,
    'H': 7,
}

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
    'CBS-QB3': GaussianLog('BenzylEnergy.log'),
}

geometry = GaussianLog('BenzylEnergy.log')

frequencies = GaussianLog('BenzylFreq.log')

rotors = [ 
    HinderedRotor(scanLog=GaussianLog('BenzylRot1.log'), pivots=[12,4], top=[12,13,14], symmetry=2, fit='best'), 
    ]
