#!/usr/bin/env python
# -*- coding: utf-8 -*-

title = 'dimetpropoxy intra-H migrations and decomposition pathways'

description = \
"""
NONE
"""

modelChemistry = "M08SO"
frequencyScaleFactor = 0.983
useHinderedRotors = True
useBondCorrections = False

species('dimetpropoxy', './dimetpropoxy.py')
transitionState('dimetpropoxy_betasci', './dimetpropoxy_betasci.py')

reaction(
    label = 'dimetpropoxy = dimetpropoxy_betasci',
    reactants = ['dimetpropoxy'],
    products = ['dimetpropoxy'],
    transitionState = 'dimetpropoxy_betasci',
#    tunneling='Eckart',
)

kinetics('dimetpropoxy = dimetpropoxy_betasci') 

statmech('dimetpropoxy_betasci')
statmech('dimetpropoxy')