#!/usr/bin/python
# -*- coding: utf-8 -*-


from nose.tools import assert_equal, assert_true, assert_false
import logging
from rmgpy.molecule.adjlist import PeriodicSystem
from rmgpy.molecule.atomtype import atomTypes
from rmgpy.molecule.atomtypedatabase import create_atom_types 
  
import itertools
import unittest
from external.wip import work_in_progress
from rmgpy.molecule.molecule import Molecule
from rmgpy.molecule.group import Group

molecule_atom_types           = [ 'C', 'O', 'N', 'S', 'Si', 'Cl']
group_atomtypes = {}

for item in create_atom_types() :
    group_atomtypes[item.label] = item
    
def get_multiplicity(unpaired_electrons):
    '''
    2*s + 1, with s = 1/2
    
    = unpaired electrons + 1 
    '''
    return unpaired_electrons+1

def get_lone_pairs(atom_type):
    if atom_type in group_atomtypes:
        return group_atomtypes[atom_type].lp#
    return PeriodicSystem.lone_pairs[atom_type]

def get_molecule_string(element, unpaired_electrons, charge):
    
    lp = get_lone_pairs(element)
    charge = '+1' if charge == 1 else charge
    
    adjList = "1 {0} u{1} p{2} c{3}".format(element,unpaired_electrons,lp,charge)
    return adjList

def get_group_string(element, unpaired_electrons, charge):
    mult = get_multiplicity(unpaired_electrons)
    s = 'multiplicity [{0}]\n'.format(mult)
    return s+get_molecule_string(element, unpaired_electrons, charge)

def createMolecule(element, u1, c1):
    from rmgpy.molecule import InvalidAdjacencyListError
    adjlist = get_molecule_string(element, u1, c1)
    logging.info('Creating molecule: {0}'.format(adjlist))
    mol = None
    try: 
        mol = Molecule().fromAdjacencyList(adjlist)
    except (ValueError, InvalidAdjacencyListError):
        # Either a non-neutral molecule, incompatible multiplcity, or an improper valence 
        # was assigned with this combination.
        pass
    return mol, adjlist

def createGroup(element, u1, c1):
    adjlist = get_group_string(element, u1, c1)
    return Group().fromAdjacencyList(adjlist), adjlist

def retrieve_unspecified_valency(atom_type, unpaired_electrons):
    if not atom_type in molecule_atom_types:
        t = group_atomtypes[atom_type]
        order = 2*t.double + 3*t.triple + int(3*t.benzene/2)
        return PeriodicSystem.valence_electrons[t.element] - 2*t.lp - unpaired_electrons - order
    return PeriodicSystem.valence_electrons[atom_type] - 2*PeriodicSystem.lone_pairs[atom_type] - unpaired_electrons

def load_cases_group_atom_types():
    '''
    creates test cases in which the atom types of the 1st graph are 'molecule' atom types, like 'C', 'N', ...
    and are 'group' atom types like 'Cs', 'Od', ... for the 2nd graph.
    '''
    
    output = []
    
    a_types = list(itertools.product(molecule_atom_types, group_atomtypes.keys()))
    
    molecule_unpaired_electrons = [0, 1, 2]#0, 1, 2
    group_unpaired_electrons = [0]#0
    u_e = list(itertools.product(molecule_unpaired_electrons, group_unpaired_electrons))
    
    cross_element_unpaired = list(itertools.product(a_types,u_e))
    
    for item in cross_element_unpaired:
        charges = []#list containing tuples of charge for graph 1 and graph 2 [(0,0), (0,1), ...]
        
        '''
        for each atom we need to determine the unspecified valency, and generate
        a list of possible charges that go along with that unspecified valency. 
        ''' 
        for el, unp in zip(item[0], item[1]):#elements, unpaired electrons
                val = retrieve_unspecified_valency(el, unp)
                '''
                for now, only allow charges up to +1, not +2, +3, even
                if the unspecified valency allows for that.
                '''
                charges.append(range(min(val,1)+1))
                
        charge_combos = list(itertools.product(charges[0],charges[1]))#cross product for both graphs
        for charge_combo in charge_combos:#combine charge tuple with the cross product of element and unpaired
            output.append(item[0]+item[1]+tuple(charge_combo))
          
    return output

def load_cases_molecule_atom_types():
    '''
    creates test cases in which the atom types are 'molecule' atom types, like 'C', 'N', ...
    for both graphs.
    '''
    output = []
    a_types           = list(itertools.product(molecule_atom_types, repeat=2))
    unpaired_electrons = list(itertools.product(range(3), repeat=2))
    cross_element_unpaired = list(itertools.product(a_types,unpaired_electrons))
    for item in cross_element_unpaired:
        charges = []#list containing tuples of charge for graph 1 and graph 2 [(0,0), (0,1), ...]
        
        '''
        for each atom we need to determine the unspecified valency, and generate
        a list of possible charges that go along with that unspecified valency. 
        ''' 
        for el, unp in zip(item[0], item[1]):#elements, unpaired electrons
                val = retrieve_unspecified_valency(el, unp)
                '''
                for now, only allow charges up to +1, not +2, +3, even
                if the unspecified valency allows for that.
                '''
                charges.append(range(min(val,1)+1))
                
        charge_combos = list(itertools.product(charges[0],charges[1]))#cross product for both graphs
        for charge_combo in charge_combos:#combine charge tuple with the cross product of element and unpaired
            output.append(item[0]+item[1]+tuple(charge_combo))
            
    return output

def mol_atom_type_comparison(e1, e2, u1, u2, c1, c2):
    return  (e1 == e2) and (c1 == c2) and (u1 == u2)

def group_atom_type_comparison(a1, a2, u1, u2, c1, c2):
    return a1.equivalent(a2) and (c1 == c2) and (u1 == u2)
   
def run_parameter_tests():
    
    def failed(*args):
        raise AssertionError

    def exception(exc):
        raise exc

    def success():
        assert_equal(True, True)
        
    def isIsomorphic_mol_atom_types(e1, e2, u1, u2, c1, c2):
        """
        Check whether isomorphism between 2 molecules consisting of each 1 atom
        perceives the difference in charge
        """
        mol1, adjList1 = createMolecule(e1, u1, c1)
        mol2, adjList2 = createMolecule(e2, u2, c2)
        
        exp = mol_atom_type_comparison(e1, e2, u1, u2, c1, c2)    
        err = "\nGraph 1: {0},\nGraph 2: {1}. \nExpected: {2}".format(adjList1, adjList2, exp)
        
        if mol1 is not None and mol2 is not None:
            calc = mol1.isIsomorphic(mol2)
            assert_equal(calc, exp, err)
    
    def findIsomorphisms_mol_atom_types(e1, e2, u1, u2, c1, c2):
        """
        Check whether isomorphism between 2 molecules consisting of each 1 atom
        perceives the difference in charge
        """
        
        mol1, adjList1 = createMolecule(e1, u1, c1)
        mol2, adjList2 = createMolecule(e2, u2, c2)

        exp = mol_atom_type_comparison(e1, e2, u1, u2, c1, c2)        
        err = "\nGraph 1: {0},\nGraph 2: {1}. \nExpected: {2}".format(adjList1, adjList2, exp)
        
        if mol1 is not None and mol2 is not None:
            calc = len(mol1.findIsomorphism(mol2)) > 0
            assert_equal(calc, exp, err)
    
    def isSubgraphIsomorphic_mol_atom_types(e1, e2, u1, u2, c1, c2):
        mol1, adjList1 = createMolecule(e1, u1, c1)
        group1, adjList2 = createGroup(e2, u2, c2)

        exp = mol_atom_type_comparison(e1, e2, u1, u2, c1, c2)#string comparison will give us expected value! 
        
        err = "\nGraph 1: {0},\nGraph 2: {1}. \nExpected: {2}".format(adjList1, adjList2, exp)
        
        if mol1 is not None and group1 is not None:
            calc = mol1.isSubgraphIsomorphic(group1)
            assert_equal(calc, exp, err)
            
    def findSubgraphIsomorphisms_mol_atom_types(e1, e2, u1, u2, c1, c2):

        mol1, adjList1 = createMolecule(e1, u1, c1)
        group1, adjList2 = createGroup(e2, u2, c2)

        exp = mol_atom_type_comparison(e1, e2, u1, u2, c1, c2)#string comparison will give us expected value! 
        
        err = "\nGraph 1: {0},\nGraph 2: {1}. \nExpected: {2}".format(adjList1, adjList2, exp)
        
        if mol1 is not None and group1 is not None:
            calc = len(mol1.findSubgraphIsomorphisms(group1)) > 0
            assert_equal(calc, exp, err)
    
    output = load_cases_molecule_atom_types()
    for args in output:
        try:
            isIsomorphic_mol_atom_types(*args)
            findIsomorphisms_mol_atom_types(*args)
            isSubgraphIsomorphic_mol_atom_types(*args)
            
        except AssertionError:
            yield (failed, args)
            
        except Exception as e:
            yield (exception, e)
    
    
    def isIsomorphic_mol_group_atom_types(e1, e2, u1, u2, c1, c2):
        """
        Check whether isomorphism between 2 molecules consisting of each 1 atom
        perceives the difference in charge
        """
        mol1, adjList1 = createMolecule(e1, u1, c1)
        group1, adjList2 = createGroup(e2, u2, c2)
        if mol1 is not None and group1 is not None:
            a1 = mol1.atoms[0].atomType
            a2 = group1.atoms[0].atomType[0]
            exp = group_atom_type_comparison(a1, a2, u1, u2, c1, c2)#string comparison will give us expected value!        
            err = "\nGraph 1: {0},\nGraph 2: {1}. \nExpected: {2}".format(adjList1, adjList2, exp)

            calc = mol1.isSubgraphIsomorphic(group1)
            assert_equal(calc, exp, err)
    
    def findSubgraphIsomorphisms_mol_group_atom_types(e1, e2, u1, u2, c1, c2):
        mol1, adjList1 = createMolecule(e1, u1, c1)
        group1, adjList2 = createGroup(e2, u2, c2)
        if mol1 is not None and group1 is not None:
            a1 = mol1.atoms[0].atomType
            a2 = group1.atoms[0].atomType[0]
            exp = group_atom_type_comparison(a1, a2, u1, u2, c1, c2)        
            err = "\nGraph 1: {0},\nGraph 2: {1}. \nExpected: {2}".format(adjList1, adjList2, exp)

            calc = len(mol1.findSubgraphIsomorphisms(group1)) > 0
            assert_equal(calc, exp, err)
            
    output = load_cases_group_atom_types()
    for args in output:
        try:
            isIsomorphic_mol_group_atom_types(*args)
            findSubgraphIsomorphisms_mol_group_atom_types(*args)
        except AssertionError:
            yield (failed, args)
        except Exception as e:
            yield (exception, e)
                    
    # Make sure that one test is always returned
    yield (success, )

def testMultiplicity_mol_mol_distinct_multiplicity():
    '''
    distinct multiplicity for both molecules set by user.
    '''
    mol = Molecule().fromAdjacencyList("""
    multiplicity 1
    1 C u1 p0 c0 {2,S}
    2 C u0 p0 c0 {1,S} {3,S}
    3 C u1 p0 c0 {2,S}
    """, saturateH=True)
    
    mol2 = Molecule().fromAdjacencyList("""
    multiplicity 3
    1 C u1 p0 c0 {2,S}
    2 C u0 p0 c0 {1,S} {3,S}
    3 C u1 p0 c0 {2,S}
    """, saturateH=True)
    
    assert_false(mol.isIsomorphic(mol2))
    assert_false(len(mol.findIsomorphism(mol2)) > 0)
    
def testMultiplicity_mol_mol_identical_multiplicity():
    '''
    identical multiplicity for both molecules set by user.
    '''
    mol = Molecule().fromAdjacencyList("""
    multiplicity 3
    1 C u2 p0 c0
    """, saturateH=True)
    
    mol2 = Molecule().fromAdjacencyList("""
    multiplicity 3
    1 C u2 p0 c0
    """, saturateH=True)
    
    assert_true(mol.isIsomorphic(mol2))
    assert_true(len(mol.findIsomorphism(mol2)) > 0)
    
def testMultiplicity_mol_not_specified_mol_specified():
    '''
    Multiplicity not set for one of two molecules
    '''
    mol = Molecule().fromAdjacencyList("""
    1 C u2 p0 c0
    """, saturateH=True)
    
    mol2 = Molecule().fromAdjacencyList("""
    multiplicity 3
    1 C u2 p0 c0
    """, saturateH=True)
    
    assert_true(mol.isIsomorphic(mol2))
    assert_true(len(mol.findIsomorphism(mol2)) > 0)
    
def testMultiplicity_mol_not_specified_mol_not_specified():
    '''
    Both multiplicities not set.
    '''
    mol = Molecule().fromAdjacencyList("""
    1 C u2 p0 c0
    """, saturateH=True)
    
    mol2 = Molecule().fromAdjacencyList("""
    1 C u2 p0 c0
    """, saturateH=True)
    
    assert_true(mol.isIsomorphic(mol2))
    assert_true(len(mol.findIsomorphism(mol2)) > 0)

def test_isomorphism_R():
    mol = Molecule().fromAdjacencyList("""
    1 C u0 p0 c0
    """, saturateH=True)
    
    gp = Group().fromAdjacencyList("""
    1 R u0 p0 c0
    """)
    
    assert_true(len(mol.findSubgraphIsomorphisms(gp)) > 0)

def test_isomorphism_mol_group_not_identical():
    """
    Testing multiplicities in mol and group that don't match
    """
    mol = Molecule().fromAdjacencyList("""
    1 C u0 p0 c0
    """, saturateH=True)
    
    gp = Group().fromAdjacencyList("""
    multiplicity [2]
    1 R u0 p0 c0
    """)
    assert_false(mol.isSubgraphIsomorphic(gp))
    assert_false(len(mol.findSubgraphIsomorphisms(gp)) > 0)

def test_isomorphism_group_group():
    """
    Testing multiplicities in group vs. group 
    """
    gp1 = Group().fromAdjacencyList("""
    1 R u0 p0 c0
    """)
    
    gp2 = Group().fromAdjacencyList("""
    multiplicity [2]
    1 R u0 p0 c0
    """)
    
    gp3 = Group().fromAdjacencyList("""
    1 C u0 p0 c0
    """)
    
    assert_false(gp1.isSubgraphIsomorphic(gp2))
    assert_false(len(gp1.findSubgraphIsomorphisms(gp2)) > 0)

    assert_true(gp2.isSubgraphIsomorphic(gp1))
    assert_true(len(gp2.findSubgraphIsomorphisms(gp1)) > 0)

    assert_false(gp2.isIdentical(gp1))
    
    assert_true(gp3.isSubgraphIsomorphic(gp1))
    assert_true(len(gp3.findSubgraphIsomorphisms(gp1)) > 0)
    assert_false(gp3.isSubgraphIsomorphic(gp2))    
    assert_false(len(gp3.findSubgraphIsomorphisms(gp2)) > 0)

@work_in_progress
def test_isomorphism_sulfurGroup_sulfurMolecule():
    """
    Test isormophism check of a CS group vs. a sulfur containing molecule
    """
    gp1 = Group().fromAdjacencyList("""
1 S  u0 {2,S} {3,S}
2 H  u0 {1,S}
3 C u0 {1,S} {4,D}
4 S  u0 {3,D}
""")
    
    gp2 = Group().fromAdjacencyList("""
1 S  u0 {2,S} {3,S}
2 H  u0 {1,S}
3 CS u0 {1,S} {4,D}
4 S  u0 {3,D}
""")
    
    mol = Molecule().fromAdjacencyList("""
1 S  0 {2,S} {3,S}
2 H  0 {1,S}
3 C 0 {1,S} {4,D} {5,S}
4 S 0 {3,D}
5 H 0 {3,S}
""")
    assert_true(mol.isSubgraphIsomorphic(gp1))
    
    assert_true(mol.isSubgraphIsomorphic(gp2))