#!/usr/bin/env python
# -*- coding: utf-8 -*-

modelChemistry = "CBS-QB3"
frequencyScaleFactor = 0.99
useHinderedRotors = True
useBondCorrections = True

species('H', 'H.py')

statmech('H')
thermo('H', 'Wilhoit')
