#!/usr/bin/python
# -*- coding: utf-8 -*-


from nose.tools import assert_equal
import logging
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

lone_pairs         = {'H': 0, 'C': 0, 'N': 1, 'O': 2, 'Si':0, 'S': 2, 'Cl':3 }
elements           = [ 'C', 'O', 'N', 'S', 'Si', 'Cl']#'H',
unpaired_electrons = list(itertools.product(range(2), repeat=2))#, '2'
valency             = {'H': 1, 'C': 4, 'N': 5, 'O': 6, 'Si':4, 'S': 6, 'Cl':7 }

def get_multiplicity(unpaired_electrons):
    '''
    2*s + 1, with s = 1/2
    
    = unpaired electrons + 1 
    '''
    return unpaired_electrons+1

def get_molecule_string(element, unpaired_electrons, charge):
    lp = lone_pairs[element]
    
    charge = '+1' if charge == 1 else charge
    
    adjList = "1 {0} u{1} p{2} c{3}".format(element,unpaired_electrons,lp,charge)
    return adjList

def get_group_string(element, unpaired_electrons, charge):
    mult = get_multiplicity(unpaired_electrons)
    s = 'multiplicity [{0}]\n'.format(mult)
    return s+get_molecule_string(element, unpaired_electrons, charge)

def createMolecule(element, u1, c1):
    adjlist = get_molecule_string(element, u1, c1)
    logging.info('Creating molecule: {0}'.format(adjlist))
    return Molecule().fromAdjacencyList(adjlist, saturateH=True), adjlist

def createGroup(element, u1, c1):
    adjlist = get_group_string(element, u1, c1)
    return Group().fromAdjacencyList(adjlist), adjlist


def retrieve_unspecified_valency(element, unpaired_electrons):
    return valency[element] - 2*lone_pairs[element] - unpaired_electrons

def load_test_cases():
    output = []
    cross_element_unpaired = list(itertools.product(elements,unpaired_electrons))
    for item in cross_element_unpaired:
        charges = []#list containing tuples of charge for graph 1 and graph 2 [(0,0), (0,1), ...]
        for unp in item[1]:#unpaired electrons
                el = item[0]
                val = retrieve_unspecified_valency(el, unp)
                charges.append(range(min(val,1)+1))
                
        charge_combos = list(itertools.product(charges[0],charges[1]))
        for charge_combo in charge_combos:
            output.append((item[0],)+item[1]+tuple(charge_combo))
          
    return output

class TestIsomorphism(unittest.TestCase):
    
    @parameterized.expand(load_test_cases)
    def testIsIsomorphic(self, element, u1, u2, c1, c2):
        """
        Check whether isomorphism between 2 molecules consisting of each 1 atom
        perceives the difference in charge
        """
        mol1, adjList1 = createMolecule(element, u1, c1)
        mol2, adjList2 = createMolecule(element, u2, c2)

        exp = (c1 == c2) and (u1 == u2)#string comparison will give us expected value!        
        err = "\nGraph 1: {0},\nGraph 2: {1}. \nExpected: {2}".format(adjList1, adjList2, exp)

        calc = mol1.isIsomorphic(mol2)
        assert_equal(calc, exp, err)
    
    @parameterized.expand(load_test_cases)
    def testFindIsomorphisms(self, element, u1, u2, c1, c2):
        """
        Check whether isomorphism between 2 molecules consisting of each 1 atom
        perceives the difference in charge
        """
        
        mol1, adjList1 = createMolecule(element, u1, c1)
        mol2, adjList2 = createMolecule(element, u2, c2)

        exp = (c1 == c2) and (u1 == u2)#string comparison will give us expected value!        
        err = "\nGraph 1: {0},\nGraph 2: {1}. \nExpected: {2}".format(adjList1, adjList2, exp)

        calc = len(mol1.findIsomorphism(mol2)) > 0
        assert_equal(calc, exp, err)
        
    @parameterized.expand(load_test_cases)
    def testIsSubgraphIsomorphic(self, element, u1, u2, c1, c2):
        
        mol1, adjList1 = createMolecule(element, u1, c1)
        group1, adjList2 = createGroup(element, u2, c2)

        exp = (c1 == c2) and (u1 == u2)#string comparison will give us expected value! 
        
        err = "\nGraph 1: {0},\nGraph 2: {1}. \nExpected: {2}".format(adjList1, adjList2, exp)
        
        calc = mol1.isSubgraphIsomorphic(group1)
        assert_equal(calc, exp, err)
        
    @parameterized.expand(load_test_cases)
    def testFindSubgraphIsomorphisms(self, element, u1, u2, c1, c2):

        mol1, adjList1 = createMolecule(element, u1, c1)
        group1, adjList2 = createGroup(element, u2, c2)

        exp = (c1 == c2) and (u1 == u2)#string comparison will give us expected value! 
        
        err = "\nGraph 1: {0},\nGraph 2: {1}. \nExpected: {2}".format(adjList1, adjList2, exp)
        
        calc = len(mol1.findSubgraphIsomorphisms(group1)) > 0
        assert_equal(calc, exp, err)
    