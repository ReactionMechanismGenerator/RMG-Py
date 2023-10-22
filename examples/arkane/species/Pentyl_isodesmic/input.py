#!/usr/bin/env python
# encoding: utf-8

"""
This input file shows how to calculate the thermochemistry of a species using isodesmic reactions in Arkane
"""
title = 'Isodesmic Example'
modelChemistry = LevelOfTheory(
    method='wb97m-v',
    basis='def2-tzvpd',
    software='qchem'
)

# Setting `useIsodesmicReactions` to True will also force useAtomCorrections to True and useBondCorrections to False
useIsodesmicReactions = True
useHinderedRotors = False
species('Pentyl', './Pentyl.py',
        structure=SMILES('[CH2]CCCC')
        )
thermo('Pentyl', 'NASA')
