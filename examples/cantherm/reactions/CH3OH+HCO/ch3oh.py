#!/usr/bin/env python
# -*- coding: utf-8 -*-

atoms = {
    'C': 1,
    'H': 4,
    'O': 1,
}

bonds = {
    'C-O': 1, 
    'C-H': 3,
    'O-H': 1,
}

linear = False

externalSymmetry = 1

spinMultiplicity = 2

opticalIsomers = 1 #confirmed

energy = {
    'M08SO': QchemLog('ch3oh.out'),
}

geometry = QchemLog('ch3oh.out')

frequencies = QchemLog('ch3oh.out')


