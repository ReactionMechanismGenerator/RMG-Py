#!/usr/bin/python
# -*- coding: utf-8 -*-


from nose.tools import assert_equal
try:
    '''
    Requires:
    
    https://github.com/wolever/nose-parameterized
    '''
    from nose_parameterized import parameterized
except:
    print 'Install nose-parameterized via: "pip install nose-parameterized" !'
    
import itertools
import unittest
from rmgpy.molecule.molecule import Molecule
from rmgpy.molecule.group import Group
################################################################################
lp = {'C': 0, 'O': 2, 'N': 1, 'S': 2}
elements = list(itertools.combinations(['C', 'O', 'N', 'S'], 1))

def createMolecule(adjlist):
    return Molecule().fromAdjacencyList(adjlist, saturateH=True)

def createGroup(adjlist):
    return Group().fromAdjacencyList(adjlist)
    
def load_test_cases_charge():
    charges  = list(itertools.combinations_with_replacement(['0', '+1', '-1'], 2))
     
    input_test = []
    global elements
    for el in elements:
        for charge_combo in charges:#prepend element to tuple
            input_test.append(el+charge_combo)
    
    return input_test

class TestIsomorphism(unittest.TestCase):
    '''
    @parameterized.expand(load_test_cases_charge)
    def testIsomorphism_with_charge(self, element, charge_mol_1, charge_mol_2):
        """
        Check whether isomorphism between 2 molecules consisting of each 1 atom
        perceives the difference in charge
        """
        global lp

        base = "1 "+element+" u0 p"+str(lp[element])+" c"
        mol1 = createMolecule(base+charge_mol_1)
        mol2 = createMolecule(base+charge_mol_2)
        err = charge_mol_1+charge_mol_2
        exp = charge_mol_1 == charge_mol_2#string comparison will give us expected value!
        calc = mol1.isIsomorphic(mol2)
        assert_equal(calc, exp, err)
    '''
    @parameterized.expand(load_test_cases_charge)
    def testSubgraphIsomorphism_with_charge(self, element, charge_mol_1, charge_mol_2):
        global lp
        
        base = "1 "+element+" u0 p"+str(lp[element])+" c"
        mol1 = createMolecule(base+charge_mol_1)
        group1 = createGroup(base+charge_mol_2)
        
        err = charge_mol_1+charge_mol_2
        exp = charge_mol_1 == charge_mol_2#string comparison will give us expected value!
        calc = len(mol1.findSubgraphIsomorphisms(group1)) > 0
        assert_equal(calc, exp, err)