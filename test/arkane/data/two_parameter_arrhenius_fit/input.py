#!/usr/bin/env python
# -*- coding: utf-8 -*-

title = 'imipramine_rate_h5'

frequencyScaleFactor = 0.99
useHinderedRotors = False
useBondCorrections = False
useAtomCorrections = False

species('imipramine', 'imipramine.py', structure=SMILES('c1cc3c(cc1)CCc2c(cccc2)N3CCCN(C)C'))
species('imipraminerad', 'imipraminerad.py', structure=SMILES('c1cc3c(cc1)CCc2c(cccc2)N3[CH]CCN(C)C'))
species('CH3OO','CH3OO.py', structure=SMILES('CO[O]'))
species('CH3OOH','CH3OOH.py', structure=SMILES('COO'))
transitionState('TS','TS.py')

reaction(
    label = 'imipramine + CH3OO = imipraminerad + CH3OOH',
    reactants = ['imipramine', 'CH3OO'],
    products = ['imipraminerad', 'CH3OOH'],
    transitionState = 'TS',
    tunneling = 'Eckart',
)

kinetics(
label = 'imipramine + CH3OO = imipraminerad + CH3OOH',
Tmin = (275,'K'), Tmax = (350,'K'), Tcount = 100, 
three_params = False,
)
