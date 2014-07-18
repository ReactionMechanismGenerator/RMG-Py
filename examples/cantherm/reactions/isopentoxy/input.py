#!/usr/bin/env python
# -*- coding: utf-8 -*-

title = 'isopentoxy intra-H migrations and decomposition pathways'

description = \
"""
This is for an initial comparison of a single channel P-dep reaction using only the
RRHO assumption. 
NOTE: for all species and transition states in this example, the number of optical isomers is 2. 
- Relative energies should contain ZPE contributions.  

To Do / Questions:
- make sure that the PES scan shows the conformation you started with is the lowest energy conformer. 
- bond corections for M08-SO/MG3S? 
"""

modelChemistry = "CBS-QB3"
frequencyScaleFactor = 0.983
useHinderedRotors = True
useBondCorrections = False

species('ic5h11o1', '../../species/ic5h11o1/ic5h11o1.py')
species('ic5h11o2geo2', '../../species/ic5h11o2/ic5h11o2geo2.py')
species('ic5h11o2', '../../species/ic5h11o2/ic5h11o2.py')
species('ic5h11o3', '../../species/ic5h11o3/ic5h11o3.py')
species('ic5h11o4', '../../species/ic5h11o4/ic5h11o4.py')
species('ic5h11o5', '../../species/ic5h11o5/ic5h11o5.py')

transitionState('ts12', './ts12.py')
transitionState('ts13', './ts13.py')
transitionState('ts14', './ts14.py')
transitionState('ts15', './ts15.py')
transitionState('ts23', './ts23.py')
transitionState('ts24', './ts24.py')
transitionState('ts25', './ts25.py')
#transitionState('ts34', './ts34.py') # FIX THIS!!!!!Wrong geo
transitionState('ts35', './ts35.py')
transitionState('ts45', './ts45.py')
transitionState('tsC5H11O-1_betasc_CH2O', './iC5H11O-1_betasc_CH2O.py')
transitionState('tsC5H11O-2_betasc_CH2CHOH', './iC5H11O-2_betasc_CH2CHOH.py')
transitionState('tsC5H11O-2_betasc_CH2CHOH2', './iC5H11O-2_betasc_CH2CHOH2.py')
#transitionState('tsC5H11O-3_betasc_CH3', './iC5H11O-3_betasc_CH3.py') #qchem job failed
#transitionState('tsC5H11O-3_betasc_OH', './iC5H11O-3_betasc_OH.py') #not accurate
transitionState('tsC5H11O-4_betasc_CH2OH', './iC5H11O-4_betasc_CH2OH.py')
transitionState('tsC5H11O-5_betasc_C3H6', './iC5H11O-5_betasc_C3H6.py')
#transitionState('iC5H11O-5_betasc_CH3', './iC5H11O-5_betasc_CH3.py')

reaction(
    label = 'ic5h11o1 = iC5H11O-1_betasc_CH2O',
    reactants = ['ic5h11o1'],
    products = ['ic5h11o2'],
    transitionState = 'tsC5H11O-1_betasc_CH2O',
    tunneling='Eckart',
)
reaction(
    label = 'ic5h11o2 = iC5H11O-2_betasc_CH2CHOH',
    reactants = ['ic5h11o2'],
    products = ['ic5h11o3'],
    transitionState = 'tsC5H11O-2_betasc_CH2CHOH',
    tunneling='Eckart',
)
reaction(
    label = 'ic5h11o2geo2 = iC5H11O-2_betasc_CH2CHOH2',
    reactants = ['ic5h11o2geo2'],
    products = ['ic5h11o3'],
    transitionState = 'tsC5H11O-2_betasc_CH2CHOH2',
    tunneling='Eckart',
)
'''
reaction(
    label = 'ic5h11o3 = iC5H11O-3_betasc_CH3',
    reactants = ['ic5h11o3'],
    products = ['ic5h11o5'],
    transitionState = 'tsC5H11O-3_betasc_CH3',
    tunneling='Eckart',
)
reaction(
    label = 'ic5h11o3 = iC5H11O-3_betasc_OH',
    reactants = ['ic5h11o3'],
    products = ['ic5h11o5'],
    transitionState = 'tsC5H11O-3_betasc_OH',
    tunneling='Eckart',
)
'''
reaction(
    label = 'ic5h11o4 = iC5H11O-4_betasc_CH2OH',
    reactants = ['ic5h11o4'],
    products = ['ic5h11o5'],
    transitionState = 'tsC5H11O-4_betasc_CH2OH',
    tunneling='Eckart',
)
reaction(
    label = 'ic5h11o5 = iC5H11O-5_betasc_C3H6',
    reactants = ['ic5h11o5'],
    products = ['ic5h11o2'],
    transitionState = 'tsC5H11O-5_betasc_C3H6',
    tunneling='Eckart',
) 
reaction(
    label = 'ic5h11o1 = ic5h11o2',
    reactants = ['ic5h11o1'],
    products = ['ic5h11o2'],
    transitionState = 'ts12',
    tunneling='Eckart',
)
reaction(
    label = 'ic5h11o1 = ic5h11o3',
    reactants = ['ic5h11o1'],
    products = ['ic5h11o3'],
    transitionState = 'ts13',
    tunneling='Eckart',
)
reaction(
    label = 'ic5h11o1 = ic5h11o4',
    reactants = ['ic5h11o1'],
    products = ['ic5h11o4'],
    transitionState = 'ts14',
    tunneling='Eckart',
)
reaction(
    label = 'ic5h11o1 = ic5h11o5',
    reactants = ['ic5h11o1'],
    products = ['ic5h11o5'],
    transitionState = 'ts15',
    tunneling='Eckart',
)
reaction(
    label = 'ic5h11o2 = ic5h11o3',
    reactants = ['ic5h11o2'],
    products = ['ic5h11o3'],
    transitionState = 'ts23',
    tunneling='Eckart',
)
reaction(
    label = 'ic5h11o2 = ic5h11o4',
    reactants = ['ic5h11o2'],
    products = ['ic5h11o4'],
    transitionState = 'ts24',
    tunneling='Eckart',
)
reaction(
    label = 'ic5h11o2 = ic5h11o5',
    reactants = ['ic5h11o2'],
    products = ['ic5h11o5'],
    transitionState = 'ts25',
    tunneling='Eckart',
)
#reaction(
#    label = 'ic5h11o3 = ic5h11o4',
#    reactants = ['ic5h11o3'],
#    products = ['ic5h11o4'],
#    transitionState = 'ts34',
#    tunneling='Eckart',
#)
reaction(
    label = 'ic5h11o3 = ic5h11o5',
    reactants = ['ic5h11o3'],
    products = ['ic5h11o5'],
    transitionState = 'ts35',
    tunneling='Eckart',
)
reaction(
    label = 'ic5h11o4 = ic5h11o5',
    reactants = ['ic5h11o4'],
    products = ['ic5h11o5'],
    transitionState = 'ts45',
    tunneling='Eckart',
)
kinetics('ic5h11o1 = ic5h11o5')
kinetics('ic5h11o1 = ic5h11o4')
kinetics('ic5h11o1 = ic5h11o2')
kinetics('ic5h11o1 = ic5h11o3')
kinetics('ic5h11o2 = ic5h11o3')
kinetics('ic5h11o2 = ic5h11o4')
kinetics('ic5h11o2 = ic5h11o5')
#kinetics('ic5h11o3 = ic5h11o4')
kinetics('ic5h11o3 = ic5h11o5')
kinetics('ic5h11o4 = ic5h11o5')
kinetics('ic5h11o1 = iC5H11O-1_betasc_CH2O')
kinetics('ic5h11o2 = iC5H11O-2_betasc_CH2CHOH')
kinetics('ic5h11o2geo2 = iC5H11O-2_betasc_CH2CHOH2')
#kinetics('ic5h11o3 = iC5H11O-3_betasc_CH3')
#kinetics('ic5h11o3 = iC5H11O-3_betasc_OH')
kinetics('ic5h11o4 = iC5H11O-4_betasc_CH2OH')
kinetics('ic5h11o5 = iC5H11O-5_betasc_C3H6')



statmech('ts14')
statmech('ts12')
statmech('ts13')
statmech('ts15')
statmech('ts23')
statmech('ts24')
statmech('ts25')
#statmech('ts34')
statmech('ts35')
statmech('ts45')
statmech('tsC5H11O-1_betasc_CH2O')
statmech('tsC5H11O-2_betasc_CH2CHOH')
statmech('tsC5H11O-2_betasc_CH2CHOH2')
#statmech('tsC5H11O-3_betasc_CH3')
#statmech('tsC5H11O-3_betasc_OH')
statmech('tsC5H11O-4_betasc_CH2OH')
statmech('tsC5H11O-5_betasc_C3H6')