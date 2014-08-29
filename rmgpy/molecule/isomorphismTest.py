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

lone_pairs         = {'H': 0, 'C': 0, 'N': 1, 'O': 2, 'Si':0, 'S': 2, 'Cl':3 }
elements           = list(itertools.combinations(['H', 'C', 'O', 'N', 'S', 'Si', 'Cl'], 1))
unpaired_electrons = list(itertools.product(['0', '1', '2'], repeat=2))
charges            = list(itertools.product(['0', '+1', '-1'], repeat=2))

def createMolecule(adjlist):
    return Molecule().fromAdjacencyList(adjlist, saturateH=True)

def createGroup(adjlist):
    return Group().fromAdjacencyList(adjlist)
    
def load_test_cases():
    
    output1 = []
    for el in elements:
        for charge_combo in charges:#prepend element to tuple
            output1.append(el+charge_combo)
    
    output = []
    for out in output1:
        for unp in unpaired_electrons:
            output.append(out+unp)
    
    return output

class TestIsomorphism(unittest.TestCase):
    
    @parameterized.expand(load_test_cases)
    def testIsIsomorphic(self, element, c1, c2, u1, u2):
        """
        Check whether isomorphism between 2 molecules consisting of each 1 atom
        perceives the difference in charge
        """
        lp = str(lone_pairs[element])
        
        adjList1 = "1 "+element+" u"+u1+" p"+lp+" c"+c1
        mol1 = createMolecule(adjList1)
        
        adjList2 = "1 "+element+" u"+u2+" p"+lp+" c"+c2
        mol2 = createMolecule(adjList2)

        exp = (c1 == c2) and (u1 == u2)#string comparison will give us expected value!        
        err = "\nGraph 1: {0},\nGraph 2: {1}. \nExpected: {2}".format(adjList1, adjList2, exp)

        calc = mol1.isIsomorphic(mol2)
        assert_equal(calc, exp, err)
        
    @parameterized.expand(load_test_cases)
    def testFindIsomorphisms(self, element, c1, c2, u1, u2):
        """
        Check whether isomorphism between 2 molecules consisting of each 1 atom
        perceives the difference in charge
        """
        lp = str(lone_pairs[element])
        
        adjList1 = "1 "+element+" u"+u1+" p"+lp+" c"+c1
        mol1 = createMolecule(adjList1)
        
        adjList2 = "1 "+element+" u"+u2+" p"+lp+" c"+c2
        mol2 = createMolecule(adjList2)

        exp = (c1 == c2) and (u1 == u2)#string comparison will give us expected value!        
        err = "\nGraph 1: {0},\nGraph 2: {1}. \nExpected: {2}".format(adjList1, adjList2, exp)

        calc = len(mol1.findIsomorphisms(mol2)) > 0
        assert_equal(calc, exp, err)
        
    @parameterized.expand(load_test_cases)
    def testIsSubgraphIsomorphic(self, element, c1, c2, u1, u2):
        
        lp = str(lone_pairs[element])
        
        adjList1 = "1 "+element+" u"+u1+" p"+lp+" c"+c1
        mol1 = createMolecule(adjList1)
        
        adjList2 = "1 "+element+" u"+u2+" p"+lp+" c"+c2
        group1 = createGroup(adjList2)

        exp = (c1 == c2) and (u1 == u2)#string comparison will give us expected value! 
        
        err = "\nGraph 1: {0},\nGraph 2: {1}. \nExpected: {2}".format(adjList1, adjList2, exp)
        
        calc = mol1.isSubgraphIsomorphic(group1)
        assert_equal(calc, exp, err)
        
    @parameterized.expand(load_test_cases)
    def testFindSubgraphIsomorphisms(self, element, c1, c2, u1, u2):
        
        lp = str(lone_pairs[element])
        
        adjList1 = "1 "+element+" u"+u1+" p"+lp+" c"+c1
        mol1 = createMolecule(adjList1)
        
        adjList2 = "1 "+element+" u"+u2+" p"+lp+" c"+c2
        group1 = createGroup(adjList2)

        exp = (c1 == c2) and (u1 == u2)#string comparison will give us expected value! 
        
        err = "\nGraph 1: {0},\nGraph 2: {1}. \nExpected: {2}".format(adjList1, adjList2, exp)
        
        calc = len(mol1.findSubgraphIsomorphisms(group1)) > 0
        assert_equal(calc, exp, err)