#!/usr/bin/env python
# -*- coding: utf-8 -*-

title = 'Example calculation of HCO + CH3OH = H2CO + CH2OH'

description = \
"""
This example is for a rigid rotor harmonic oscillator TST calculation calclated at the 
M08SO/MG3S level of DFT with a tight grid.  Tunneling is not included in this calculation.
Last, energy corrections for M08SO have not been added to cantherm yet, so you will see a 
warning: Unknown model chemistry "M08SO"; not applying energy corrections. Because this is
TST calculation and not a thermochemistry calculation, this can be disregarded
"""

modelChemistry = "M08SO"
frequencyScaleFactor = 0.983
useHinderedRotors = False
useBondCorrections = False


species('ch3oh', './ch3oh.py')
transitionState('ts', './ts.py')
species('hco', './hco.py')

reaction(
    label = 'hco + ch3oh <=> hco + ch3oh',
    reactants = ['hco', 'ch3oh'],
    products = ['hco', 'ch3oh'], #product channels are only important if tunneling is computed
    transitionState = 'ts',
#    tunneling='Eckart', #we dont want to comput tunneling for this calculation
)

statmech('ts')
kinetics('hco + ch3oh <=> hco + ch3oh')
