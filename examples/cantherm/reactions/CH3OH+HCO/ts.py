#!/usr/bin/env python
# -*- coding: utf-8 -*-

bonds = {}

linear = False

externalSymmetry = 1

spinMultiplicity = 2

opticalIsomers = 2 #confirmed

energy = {
    'm08so/mg3s*': Log('ts1.out'),
}

geometry = Log('ts1.out')

frequencies = Log('ts1.out')

rotors = []
