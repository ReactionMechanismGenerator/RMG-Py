#!/usr/bin/env python
# -*- coding: utf-8 -*-

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
    'm08so/mg3s*': Log('ch3oh.out'),
}

geometry = Log('ch3oh.out')

frequencies = Log('ch3oh.out')


