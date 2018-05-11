#!/usr/bin/env python
# encoding: utf-8

bonds = {
    'C-C': 1,
    'C-H': 6,
}

linear = False

externalSymmetry = 6

spinMultiplicity = 1

opticalIsomers = 1

energy = {
    'CBS-QB3': Log('ethane_cbsqb3.log'),
    'Klip_2': -79.64199436,
}

geometry = Log('ethane_cbsqb3.log')

frequencies = Log('ethane_cbsqb3.log')

frequencyScaleFactor = 0.99

"""pivot are the two atoms that are attached to the rotor
top contains the atoms that are being rotated including one of the atoms from pivots
symmetry is the symmetry number of the scan
fit is fit of the scan data. It defaults to 'best', but can also be assigned as 'cosine' or 'fourier'"""
rotors = [
    HinderedRotor(scanLog=Log('ethane_scan_1.log'), pivots=[1,5], top=[1,2,3,4], symmetry=3, fit='best'),
]
