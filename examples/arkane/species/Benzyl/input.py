#!/usr/bin/env python
# -*- coding: utf-8 -*-

modelChemistry = "CBS-QB3"
useHinderedRotors = True
useBondCorrections = True

species('C7H7', 'benzyl.py')

statmech('C7H7')
thermo('C7H7', 'NASA')