#!/usr/bin/env python
# -*- coding: utf-8 -*-

bonds = {'N-N': 3, 'N-H': 6}

externalSymmetry = 2

spinMultiplicity = 1

opticalIsomers = 2

energy = {'CCSD(T)-F12/cc-pVTZ-f12': Log('sp.out')}

geometry = Log('freq.log')

frequencies = Log('freq.log')
