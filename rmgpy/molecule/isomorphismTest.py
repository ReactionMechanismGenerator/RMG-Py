#!/usr/bin/python
# -*- coding: utf-8 -*-

from nose.tools import assert_equal
try:
    from nose_parameterized import parameterized
except:
    print 'Install nose-parameterized via: "pip install nose-parameterized" !'
    
import itertools
import unittest
from rmgpy.molecule.molecule import Molecule

################################################################################
lp = {'C': 0, 'O': 2, 'N': 1, 'S': 2}

def load_test_cases():
    charges  = list(itertools.combinations_with_replacement(['0', '+1', '-1'], 2))
    elements = list(itertools.combinations(['C', 'O', 'N', 'S'], 1)) 
    input_test = []
    for el in elements:
        for charge_combo in charges:#prepend element to tuple
            input_test.append(el+charge_combo)
    
    return input_test

class TestIsomorphism(unittest.TestCase):
    def create(self, adjlist):
        return Molecule().fromAdjacencyList(adjlist, saturateH=True)
    
    '''
    requires:
    
    https://github.com/wolever/nose-parameterized
    
    Install via:
    
    pip install nose-parameterized
    
    '''
    @parameterized.expand(load_test_cases)
    def testIsomorphism_with_charge(self, element, charge_mol_1, charge_mol_2):
        '''
        Check whether isomorphism sees difference in charge
        '''
        global lp

        base = "1 "+element+" u0 p"+str(lp[element])+" c"
        mol1 = self.create(base+charge_mol_1)
        mol2 = self.create(base+charge_mol_2)
        err = charge_mol_1+charge_mol_2
        exp = charge_mol_1 == charge_mol_2#string comparison will give us expected value!
        assert_equal(mol1.isIsomorphic(mol2), exp, err)
################################################################################
        