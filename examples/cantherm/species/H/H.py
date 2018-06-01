#!/usr/bin/env python
# -*- coding: utf-8 -*-

bonds = {}

linear = False

externalSymmetry = 1

spinMultiplicity = 2

opticalIsomers = 1

energy = {
    'CBS-QB3': Log('H_cbsqb3.log'),
    'Klip_2': -0.50003976,
}

geometry = Log('H_cbsqb3.log')

frequencies = Log('H_cbsqb3.log')

rotors = []
