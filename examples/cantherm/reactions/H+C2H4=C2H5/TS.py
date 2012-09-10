#!/usr/bin/env python
# -*- coding: utf-8 -*-

atoms = {
    'C': 2,
    'H': 5,
}

bonds = {}

linear = False

externalSymmetry = 1

spinMultiplicity = 2

opticalIsomers = 1

energy = {
    'CBS-QB3': GaussianLog('TS_energy.log'),
}

geometry = GaussianLog('TS_freq.log')

frequencies = GaussianLog('TS_freq.log')

rotors = []
