#!/usr/bin/env python
# -*- coding: utf-8 -*-

bonds = {}

linear = False

externalSymmetry = 1

spinMultiplicity = 2

opticalIsomers = 2 #confirmed

energy = {
    'M08SO': QchemLog('ts1.out'),
}

geometry = QchemLog('ts1.out')

frequencies = QchemLog('ts1.out')

rotors = []
