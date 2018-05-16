#!/usr/bin/env python
# -*- coding: utf-8 -*-

bonds = {
    'C-O': 1, 
    'C-H': 1,
}

linear = False

externalSymmetry = 1

spinMultiplicity = 2

opticalIsomers = 1

energy = {
    'm08so/mg3s*': Log('hco.out'),
}

geometry = Log('hco.out')

frequencies = Log('hco.out')


