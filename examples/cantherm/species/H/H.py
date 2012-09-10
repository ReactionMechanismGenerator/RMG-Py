#!/usr/bin/env python
# -*- coding: utf-8 -*-

atoms = {
    'H': 1,
}

bonds = {}

linear = False

externalSymmetry = 1

spinMultiplicity = 2

opticalIsomers = 1

energy = {
    'CBS-QB3': GaussianLog('H_cbsqb3.log'),
    'Klip_2': -0.50003976,
}

geometry = GaussianLog('H_cbsqb3.log')

frequencies = GaussianLog('H_cbsqb3.log')

rotors = []
