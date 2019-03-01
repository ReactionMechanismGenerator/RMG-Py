#!/usr/bin/env python
# -*- coding: utf-8 -*-

title = 'NH2 + N2H3 = NH + N2H4'
description = """
This examples shows how to use species (and TS) YAML files to run Arkane
"""

modelChemistry = "CCSD(T)-F12/aug-cc-pVTZ"

useHinderedRotors = True
useBondCorrections = False

species('NH', 'NH.yml')
species('NH2','NH2.yml')
species('N2H3','N2H3.yml')
species('N2H4','N2H4.yml')
transitionState('TS1','TS1.yml')

reaction(
    label = 'NH2 + N2H3 = NH + N2H4',
    reactants = ['NH2','N2H3'],
    products = ['NH','N2H4'],
    transitionState = 'TS1',
    tunneling = 'Eckart',
)

kinetics(
label = 'NH2 + N2H3 = NH + N2H4',
Tmin = (500,'K'), Tmax = (3000,'K'), Tcount = 30,
)
