#!/usr/bin/env python
# -*- coding: utf-8 -*-

modelChemistry = "Klip_2"
frequencyScaleFactor = 0.99
useHinderedRotors = True
useBondCorrections = True

species('C2H5', 'ethyl.py')

statmech('C2H5')
thermo('C2H5', 'Wilhoit')
