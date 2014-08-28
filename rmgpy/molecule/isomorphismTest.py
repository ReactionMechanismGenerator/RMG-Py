#!/usr/bin/python
# -*- coding: utf-8 -*-

from nose.tools import assert_equal
from nose_parameterized import parameterized
import itertools
import unittest
from rmgpy.molecule.molecule import Molecule

################################################################################

def load_test_cases():
    allowed_charges = ['0', '+1', '-1']
    return list(itertools.combinations_with_replacement(allowed_charges, 2))

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
    def testIsomorphism_with_charge(self, charge_mol_1, charge_mol_2):
        '''
        Check whether isomorphism sees difference in charge
        '''
        base = "1 C u0 p0 c"
        mol1 = self.create(base+charge_mol_1)
        mol2 = self.create(base+charge_mol_2)
        err = charge_mol_1+charge_mol_2
        exp = charge_mol_1 == charge_mol_2#string comparison will give us expected value!
        assert_equal(mol1.isIsomorphic(mol2), exp, err)
################################################################################
