#!/usr/bin/env python
# -*- coding: utf-8 -*-

title = 'N4H6 thermo'

modelChemistry = "CCSD(T)-F12/cc-pVTZ-f12//B3LYP/6-311+G(3df,2p)"

useHinderedRotors = True
useBondCorrections = True

species('N4H6', 'N4H6.py')

thermo('N4H6', 'NASA')
