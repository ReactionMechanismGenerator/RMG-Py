#!/usr/bin/env python
# encoding: utf-8

modelChemistry = "CBS-QB3"
useHinderedRotors = True
useBondCorrections = False

species('C2H4', 'ethene.py',
        structure=SMILES('C=C'))

thermo('C2H4', 'NASA')
