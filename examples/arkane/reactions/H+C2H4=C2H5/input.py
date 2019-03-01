#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This example also generates YAML files for all species and TS
# A YAML file is generated for species if they're structure is defined and Thermo() is called,
# and for TS if all the respective reactant/s and product/s have structures

modelChemistry = "CBS-QB3"
frequencyScaleFactor = 0.99
useHinderedRotors = False
useBondCorrections = True

species('H', '../../species/H/H.py',
       structure=SMILES('[H]'))
species('C2H4', '../../species/C2H4/ethene.py',
       structure=SMILES('C=C'))
species('C2H5', '../../species/C2H5/ethyl.py',
       structure=SMILES('[CH2]C'))
transitionState('TS', 'TS.py')

thermo('H','NASA')
thermo('C2H4','NASA')
thermo('C2H5','NASA')

reaction(
    label = 'H + C2H4 <=> C2H5',
    reactants = ['H', 'C2H4'],
    products = ['C2H5'],
    transitionState = 'TS',
    tunneling='Eckart',
)

statmech('TS')
kinetics(    
	label = 'H + C2H4 <=> C2H5',
    Tmin = (400,'K'), Tmax = (1200,'K'), Tcount = 6, # this can be changed to any desired temperature range with any number of temperatures
    Tlist = ([400,500,700,900,1100,1200],'K'),
    )
