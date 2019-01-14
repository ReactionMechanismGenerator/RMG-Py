#!/usr/bin/env python
# -*- coding: utf-8 -*-

bonds = {}

linear = False

externalSymmetry = 1

spinMultiplicity = 2

opticalIsomers = 1

energy = {
    'CBS-QB3': Log('TS_energy.log'),
}

geometry = Log('TS_freq.log')

frequencies = Log('TS_freq.log')

rotors = []
