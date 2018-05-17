#!/usr/bin/env python
# -*- coding: utf-8 -*-

modelChemistry = "CBS-QB3"
useHinderedRotors = True
useBondCorrections = True

species('C2H5', 'ethyl.py')

statmech('C2H5')
thermo('C2H5', 'Wilhoit')
