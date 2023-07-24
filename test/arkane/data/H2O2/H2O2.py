#!/usr/bin/env python
# -*- coding: utf-8 -*-

bonds = {'H-O': 2, 'O-O': 1}

externalSymmetry = 2

spinMultiplicity = 1

opticalIsomers = 1

energy = {'b3lyp/6-311+g(3df,2p)': Log('/home/jackson/rmg/RMG-Py/test/arkane/data/H2O2/sp_a19032.out')}

geometry = Log('/home/jackson/rmg/RMG-Py/test/arkane/data/H2O2/freq_a19031.out')

frequencies = Log('/home/jackson/rmg/RMG-Py/test/arkane/data/H2O2/freq_a19031.out')

rotors = [HinderedRotor(scanLog=ScanLog('/home/jackson/rmg/RMG-Py/test/arkane/data/H2O2/scan.txt'), pivots=[1, 2], top=[1, 3], symmetry=1, fit='fourier')]

