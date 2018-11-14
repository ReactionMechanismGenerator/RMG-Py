#!/usr/bin/env python
# -*- coding: utf-8 -*-

bonds = {
    'C-C': 1, 
    'C-H': 5,
}

linear = False

externalSymmetry = 1

spinMultiplicity = 2

opticalIsomers = 1

energy = {
    'CBS-QB3': Log('ethyl_cbsqb3.log'),
    'Klip_2': -78.98344186,
}

geometry = Log('ethyl_cbsqb3.log')

frequencies = Log('ethyl_cbsqb3.log')

"""pivot are the two atoms that are attached to the rotor
top contains the atoms that are being rotated including one of the atoms from pivots
symmetry is the symmetry number of the scan
fit is fit of the scan data. It defaults to 'best', but can also be assigned as 'cosine' or 'fourier'"""
rotors = [
    HinderedRotor(scanLog=Log('ethyl_scan_72.log'), pivots=[1,2], top=[1,3,4], symmetry=6, fit='best')
]
