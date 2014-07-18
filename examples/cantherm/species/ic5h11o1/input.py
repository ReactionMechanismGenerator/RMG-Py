#!/usr/bin/env python
# -*- coding: utf-8 -*-

modelChemistry = "CBS-QB3"
frequencyScaleFactor = 0.983
useHinderedRotors = True
useBondCorrections = False

species('ic5h11o1', 'ic5h11o1.py')

statmech('ic5h11o1')
thermo('ic5h11o1', 'Wilhoit')
