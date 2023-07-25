#!/usr/bin/env python
# -*- coding: utf-8 -*-

modelChemistry = "CBS-QB3"
useHinderedRotors = True
useBondCorrections = False

species('CH2CHOOH', 'CH2CHOOH.py')

statmech('CH2CHOOH')
thermo('CH2CHOOH', 'Wilhoit')
