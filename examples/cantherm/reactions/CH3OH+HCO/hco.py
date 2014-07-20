#!/usr/bin/env python
# -*- coding: utf-8 -*-

atoms = {
    'C': 1,
    'H': 1,
    'O': 1,
}

bonds = {
    'C-O': 1, 
    'C-H': 1,
}

linear = False

externalSymmetry = 1

spinMultiplicity = 2

opticalIsomers = 1

energy = {
    'M08SO': QchemLog('hco.out'),
}

geometry = QchemLog('hco.out')

frequencies = QchemLog('hco.out')


