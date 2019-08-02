#!/usr/bin/env python
# encoding: utf-8

"""
This input file shows how to calculate the thermochemistry of a species using isodesmic reactions in Arkane
"""
title = '[CH2]CCCCCCCC Isodesmic Example'
modelChemistry = 'wB97x_D3BJ/def2-TZVP'
frequencyScaleFactor = 1

# Set the `useIsodesmicReactions` flag to True, which will also force useAtomCorrections and useBondCorrections to False
useIsodesmicReactions = True
useHinderedRotors = False
species('Nonyl', './Nonyl.py',
        structure=SMILES('[CH2]CCCCCCCC')
        )
thermo('Nonyl', 'NASA')
