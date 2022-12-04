#!/usr/bin/env python
# encoding: utf-8

# Define the level of theory ("model chemistry"):
modelChemistry = LevelOfTheory(method='CCSD(T)-F12', basis='cc-pVTZ-F12', software='molpro')
useHinderedRotors = True

# Define the species and TSs:
species('CH4', 'data/ch4.py', structure=SMILES('C'))
species('NH2', 'data/nh2.py', structure=SMILES('[NH2]'))
species('CH3', 'data/ch3.py', structure=SMILES('[CH3]'))
species('NH3', 'data/nh3.py', structure=SMILES('N'))
species('C2H6', 'data/c2h6.py', structure=SMILES('CC'))
species('C2H5', 'data/c2h5.py', structure=SMILES('C[CH2]'))
transitionState('TS1', 'data/ts1.py')
transitionState('TS2', 'data/ts2.py')


# Define the reactions:
reaction(
    label = 'CH4 + NH2 <=> CH3 + NH3',
    reactants = ['CH4', 'NH2'],
    products = ['CH3', 'NH3'],
    transitionState = 'TS1',
    tunneling='Eckart',
)
reaction(
    label = 'C2H6 + NH2 <=> C2H5 + NH3',
    reactants = ['C2H6', 'NH2'],
    products = ['C2H5', 'NH3'],
    transitionState = 'TS2',
    tunneling='Eckart',
)


# Request rate coefficient computations:
kinetics(label='CH4 + NH2 <=> CH3 + NH3',
         Tmin=(300,'K'), Tmax=(2500,'K'), Tcount = 25)
kinetics(label='C2H6 + NH2 <=> C2H5 + NH3',
         Tmin=(300,'K'), Tmax=(2500,'K'), Tcount = 25)

