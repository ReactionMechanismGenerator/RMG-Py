#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2020 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

import itertools
import logging

from nose.tools import assert_equal, assert_true, assert_false

from rmgpy.molecule.atomtypedatabase import create_atom_types
from rmgpy.molecule.element import PeriodicSystem
from rmgpy.molecule.group import Group
from rmgpy.molecule.molecule import Molecule

molecule_atom_types = ['C', 'O', 'N', 'S', 'P', 'Si', 'Cl', 'Br', 'I', 'F']
group_atomtypes = {}

for item in create_atom_types():
    group_atomtypes[item.label] = item


def get_multiplicity(unpaired_electrons):
    """
    2*s + 1, with s = 1/2

    = unpaired electrons + 1
    """
    return unpaired_electrons + 1


def get_lone_pairs(atom_type):
    if atom_type in group_atomtypes:
        return group_atomtypes[atom_type].lp
    return PeriodicSystem.lone_pairs[atom_type]


def get_molecule_string(element, unpaired_electrons, charge):
    lp = get_lone_pairs(element)
    charge = '+1' if charge == 1 else charge

    adjList = "1 {0} u{1} p{2} c{3}".format(element, unpaired_electrons, lp, charge)
    return adjList


def get_group_string(element, unpaired_electrons, charge):
    mult = get_multiplicity(unpaired_electrons)
    s = 'multiplicity [{0}]\n'.format(mult)
    return s + get_molecule_string(element, unpaired_electrons, charge)


def create_molecule(element, u1, c1):
    from rmgpy.molecule.adjlist import InvalidAdjacencyListError
    adjlist = get_molecule_string(element, u1, c1)
    logging.info('Creating molecule: {0}'.format(adjlist))
    mol = None
    if element == 'P' and c1 == 1 and get_lone_pairs(element) == 1:
        # AtomTypeError is raised when phosphorus has +1 partial charge and 1 lone pair.
        # Although this will not raise InvalidAdjacencyListError, these are chemically infeasible.
        pass
    else:
        try:
            mol = Molecule().from_adjacency_list(adjlist)
        except (ValueError, InvalidAdjacencyListError):
            # Either a non-neutral molecule, incompatible multiplcity, or an improper valence
            # was assigned with this combination.
            pass
    return mol, adjlist


def create_group(element, u1, c1):
    adjlist = get_group_string(element, u1, c1)
    return Group().from_adjacency_list(adjlist), adjlist


def retrieve_unspecified_valency(atom_type, unpaired_electrons):
    if atom_type not in molecule_atom_types:
        t = group_atomtypes[atom_type]
        order = 2 * t.double + 3 * t.triple + int(3 * t.benzene / 2)
        return PeriodicSystem.valence_electrons[t.element] - 2 * t.lp - unpaired_electrons - order
    return PeriodicSystem.valence_electrons[atom_type] - 2 * PeriodicSystem.lone_pairs[atom_type] - unpaired_electrons


def load_cases_group_atom_types():
    """
    creates test cases in which the atom types of the 1st graph are 'molecule' atom types, like 'C', 'N', ...
    and are 'group' atom types like 'Cs', 'O2d', ... for the 2nd graph.
    """

    output = []

    a_types = list(itertools.product(molecule_atom_types, list(group_atomtypes.keys())))

    molecule_unpaired_electrons = [0, 1, 2]
    group_unpaired_electrons = [0]
    u_e = list(itertools.product(molecule_unpaired_electrons, group_unpaired_electrons))

    cross_element_unpaired = list(itertools.product(a_types, u_e))

    for item in cross_element_unpaired:
        charges = []  # list containing tuples of charge for graph 1 and graph 2 [(0,0), (0,1), ...]

        # for each atom we need to determine the unspecified valency, and generate
        # a list of possible charges that go along with that unspecified valency.
        for el, unp in zip(item[0], item[1]):  # elements, unpaired electrons
            val = retrieve_unspecified_valency(el, unp)
            # for now, only allow charges up to +1, not +2, +3, even
            # if the unspecified valency allows for that.
            charges.append(list(range(min(val, 1) + 1)))

        charge_combos = list(itertools.product(charges[0], charges[1]))  # cross product for both graphs
        for charge_combo in charge_combos:  # combine charge tuple with the cross product of element and unpaired
            output.append(item[0] + item[1] + tuple(charge_combo))

    return output


def load_cases_molecule_atom_types():
    """
    creates test cases in which the atom types are 'molecule' atom types, like 'C', 'N', ...
    for both graphs.
    """
    output = []
    a_types = list(itertools.product(molecule_atom_types, repeat=2))
    uncharged_a_types = ['Cl', 'Br', 'I', 'F']
    unpaired_electrons = list(itertools.product(list(range(3)), repeat=2))
    cross_element_unpaired = list(itertools.product(a_types, unpaired_electrons))
    for item in cross_element_unpaired:
        charges = []  # list containing tuples of charge for graph 1 and graph 2 [(0,0), (0,1), ...]

        # for each atom we need to determine the unspecified valency, and generate
        # a list of possible charges that go along with that unspecified valency.
        for el, unp in zip(item[0], item[1]):  # elements, unpaired electrons
            val = retrieve_unspecified_valency(el, unp)
            # for now, only allow charges up to +1, not +2, +3, even
            # if the unspecified valency allows for that.
            if el not in uncharged_a_types:
                charges.append(list(range(min(val, 1) + 1)))
            else:
                charges.append((0, 0))

        charge_combos = list(itertools.product(charges[0], charges[1]))  # cross product for both graphs
        for charge_combo in charge_combos:  # combine charge tuple with the cross product of element and unpaired
            output.append(item[0] + item[1] + tuple(charge_combo))

    return output


def mol_atom_type_comparison(e1, e2, u1, u2, c1, c2):
    return (e1 == e2) and (c1 == c2) and (u1 == u2)


def group_atom_type_comparison(a1, a2, u1, u2, c1, c2):
    return a1.equivalent(a2) and (c1 == c2) and (u1 == u2)


def run_parameter_tests():
    def failed(*args):
        raise AssertionError

    def exception(exc):
        raise exc

    def success():
        assert_equal(True, True)

    def is_isomorphic_mol_atom_types(e1, e2, u1, u2, c1, c2):
        """
        Check whether isomorphism between 2 molecules consisting of each 1 atom
        perceives the difference in charge
        """
        mol1, adjlist1 = create_molecule(e1, u1, c1)
        mol2, adjlist2 = create_molecule(e2, u2, c2)

        exp = mol_atom_type_comparison(e1, e2, u1, u2, c1, c2)
        err = "\nGraph 1: {0},\nGraph 2: {1}. \nExpected: {2}".format(adjlist1, adjlist2, exp)

        if mol1 is not None and mol2 is not None:
            calc = mol1.is_isomorphic(mol2)
            assert_equal(calc, exp, err)

    def find_isomorphisms_mol_atom_types(e1, e2, u1, u2, c1, c2):
        """
        Check whether isomorphism between 2 molecules consisting of each 1 atom
        perceives the difference in charge
        """

        mol1, adjlist1 = create_molecule(e1, u1, c1)
        mol2, adjlist2 = create_molecule(e2, u2, c2)

        exp = mol_atom_type_comparison(e1, e2, u1, u2, c1, c2)
        err = "\nGraph 1: {0},\nGraph 2: {1}. \nExpected: {2}".format(adjlist1, adjlist2, exp)

        if mol1 is not None and mol2 is not None:
            calc = len(mol1.find_isomorphism(mol2)) > 0
            assert_equal(calc, exp, err)

    def is_subgraph_isomorphic_mol_atom_types(e1, e2, u1, u2, c1, c2):
        mol1, adjlist1 = create_molecule(e1, u1, c1)
        group1, adjlist2 = create_group(e2, u2, c2)

        exp = mol_atom_type_comparison(e1, e2, u1, u2, c1, c2)  # string comparison will give us expected value!

        err = "\nGraph 1: {0},\nGraph 2: {1}. \nExpected: {2}".format(adjlist1, adjlist2, exp)

        if mol1 is not None and group1 is not None:
            calc = mol1.is_subgraph_isomorphic(group1)
            assert_equal(calc, exp, err)

    def find_subgraph_isomorphisms_mol_atom_types(e1, e2, u1, u2, c1, c2):

        mol1, adjlist1 = create_molecule(e1, u1, c1)
        group1, adjlist2 = create_group(e2, u2, c2)

        exp = mol_atom_type_comparison(e1, e2, u1, u2, c1, c2)  # string comparison will give us expected value!

        err = "\nGraph 1: {0},\nGraph 2: {1}. \nExpected: {2}".format(adjlist1, adjlist2, exp)

        if mol1 is not None and group1 is not None:
            calc = len(mol1.find_subgraph_isomorphisms(group1)) > 0
            assert_equal(calc, exp, err)

    output = load_cases_molecule_atom_types()
    for args in output:
        try:
            is_isomorphic_mol_atom_types(*args)
            find_isomorphisms_mol_atom_types(*args)
            is_subgraph_isomorphic_mol_atom_types(*args)

        except AssertionError:
            yield (failed, args)

        except Exception as e:
            yield (exception, e)

    def is_isomorphic_mol_group_atom_types(e1, e2, u1, u2, c1, c2):
        """
        Check whether isomorphism between 2 molecules consisting of each 1 atom
        perceives the difference in charge
        """
        mol1, adjlist1 = create_molecule(e1, u1, c1)
        group1, adjlist2 = create_group(e2, u2, c2)
        if mol1 is not None and group1 is not None:
            a1 = mol1.atoms[0].atomtype
            a2 = group1.atoms[0].atomtype[0]
            exp = group_atom_type_comparison(a1, a2, u1, u2, c1,
                                             c2)  # string comparison will give us expected value!
            err = "\nGraph 1: {0},\nGraph 2: {1}. \nExpected: {2}".format(adjlist1, adjlist2, exp)

            calc = mol1.is_subgraph_isomorphic(group1)
            assert_equal(calc, exp, err)

    def find_subgraph_isomorphisms_mol_group_atom_types(e1, e2, u1, u2, c1, c2):
        mol1, adjlist1 = create_molecule(e1, u1, c1)
        group1, adjlist2 = create_group(e2, u2, c2)
        if mol1 is not None and group1 is not None:
            a1 = mol1.atoms[0].atomtype
            a2 = group1.atoms[0].atomtype[0]
            exp = group_atom_type_comparison(a1, a2, u1, u2, c1, c2)
            err = "\nGraph 1: {0},\nGraph 2: {1}. \nExpected: {2}".format(adjlist1, adjlist2, exp)

            calc = len(mol1.find_subgraph_isomorphisms(group1)) > 0
            assert_equal(calc, exp, err)

    output = load_cases_group_atom_types()
    for args in output:
        try:
            is_isomorphic_mol_group_atom_types(*args)
            find_subgraph_isomorphisms_mol_group_atom_types(*args)
        except AssertionError:
            yield (failed, args)
        except Exception as e:
            yield (exception, e)

    # Make sure that one test is always returned
    yield (success,)


def test_multiplicity_mol_mol_distinct_multiplicity():
    """
    distinct multiplicity for both molecules set by user.
    """
    mol = Molecule().from_adjacency_list("""
    multiplicity 1
    1 C u1 p0 c0 {2,S}
    2 C u0 p0 c0 {1,S} {3,S}
    3 C u1 p0 c0 {2,S}
    """, saturate_h=True)

    mol2 = Molecule().from_adjacency_list("""
    multiplicity 3
    1 C u1 p0 c0 {2,S}
    2 C u0 p0 c0 {1,S} {3,S}
    3 C u1 p0 c0 {2,S}
    """, saturate_h=True)

    assert_false(mol.is_isomorphic(mol2))
    assert_false(len(mol.find_isomorphism(mol2)) > 0)


def test_multiplicity_mol_mol_identical_multiplicity():
    """
    identical multiplicity for both molecules set by user.
    """
    mol = Molecule().from_adjacency_list("""
    multiplicity 3
    1 C u2 p0 c0
    """, saturate_h=True)

    mol2 = Molecule().from_adjacency_list("""
    multiplicity 3
    1 C u2 p0 c0
    """, saturate_h=True)

    assert_true(mol.is_isomorphic(mol2))
    assert_true(len(mol.find_isomorphism(mol2)) > 0)


def test_multiplicity_mol_not_specified_mol_specified():
    """
    Multiplicity not set for one of two molecules
    """
    mol = Molecule().from_adjacency_list("""
    1 C u2 p0 c0
    """, saturate_h=True)

    mol2 = Molecule().from_adjacency_list("""
    multiplicity 3
    1 C u2 p0 c0
    """, saturate_h=True)

    assert_true(mol.is_isomorphic(mol2))
    assert_true(len(mol.find_isomorphism(mol2)) > 0)


def test_multiplicity_mol_not_specified_mol_not_specified():
    """
    Both multiplicities not set.
    """
    mol = Molecule().from_adjacency_list("""
    1 C u2 p0 c0
    """, saturate_h=True)

    mol2 = Molecule().from_adjacency_list("""
    1 C u2 p0 c0
    """, saturate_h=True)

    assert_true(mol.is_isomorphic(mol2))
    assert_true(len(mol.find_isomorphism(mol2)) > 0)


def test_isomorphism__r():
    mol = Molecule().from_adjacency_list("""
    1 C u0 p0 c0
    """, saturate_h=True)

    gp = Group().from_adjacency_list("""
    1 R u0 p0 c0
    """)

    assert_true(len(mol.find_subgraph_isomorphisms(gp)) > 0)


def test_isomorphism_mol_group_not_identical():
    """
    Testing multiplicities in mol and group that don't match
    """
    mol = Molecule().from_adjacency_list("""
    1 C u0 p0 c0
    """, saturate_h=True)

    gp = Group().from_adjacency_list("""
    multiplicity [2]
    1 R u0 p0 c0
    """)
    assert_false(mol.is_subgraph_isomorphic(gp))
    assert_false(len(mol.find_subgraph_isomorphisms(gp)) > 0)


def test_isomorphism_group_group():
    """
    Testing multiplicities in group vs. group 
    """
    gp1 = Group().from_adjacency_list("""
    1 R u0 p0 c0
    """)

    gp2 = Group().from_adjacency_list("""
    multiplicity [2]
    1 R u0 p0 c0
    """)

    gp3 = Group().from_adjacency_list("""
    1 C u0 p0 c0
    """)

    assert_false(gp1.is_subgraph_isomorphic(gp2))
    assert_false(len(gp1.find_subgraph_isomorphisms(gp2)) > 0)

    assert_true(gp2.is_subgraph_isomorphic(gp1))
    assert_true(len(gp2.find_subgraph_isomorphisms(gp1)) > 0)

    assert_false(gp2.is_identical(gp1))

    assert_true(gp3.is_subgraph_isomorphic(gp1))
    assert_true(len(gp3.find_subgraph_isomorphisms(gp1)) > 0)
    assert_false(gp3.is_subgraph_isomorphic(gp2))
    assert_false(len(gp3.find_subgraph_isomorphisms(gp2)) > 0)


def test_isomorphism_sulfur_group_sulfur_molecule():
    """
    Test isormophism check of a CS group vs. a sulfur containing molecule
    """
    gp1 = Group().from_adjacency_list("""
1 S  u0 {2,S} {3,S}
2 H  u0 {1,S}
3 C u0 {1,S} {4,D}
4 S  u0 {3,D}
""")

    gp2 = Group().from_adjacency_list("""
1 S  u0 {2,S} {3,S}
2 H  u0 {1,S}
3 CS u0 {1,S} {4,D}
4 S  u0 {3,D}
""")

    mol = Molecule().from_adjacency_list("""
1 S  0 {2,S} {3,S}
2 H  0 {1,S}
3 C 0 {1,S} {4,D} {5,S}
4 S 0 {3,D}
5 H 0 {3,S}
""")
    assert_true(mol.is_subgraph_isomorphic(gp1))

    assert_true(mol.is_subgraph_isomorphic(gp2))


def test_isotope_subgraph_isomorphism_molecule_and_group():
    """
    Checks that subgraph isomorphism works with enriched molecules and groups
    """
    methanei = Molecule().from_adjacency_list("""
    1 C u0 p0 c0 i13 {2,S} {3,S} {4,S} {5,S}
    2 H u0 p0 c0 {1,S}
    3 H u0 p0 c0 {1,S}
    4 H u0 p0 c0 {1,S}
    5 H u0 p0 c0 {1,S}
    """)
    methane = Molecule().from_adjacency_list("""
    1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
    2 H u0 p0 c0 {1,S}
    3 H u0 p0 c0 {1,S}
    4 H u0 p0 c0 {1,S}
    5 H u0 p0 c0 {1,S}
    """)
    group_methane = Group().from_adjacency_list("""
    1 C  u0 {2,S} {3,S} {4,S} {5,S}
    2 H  u0 {1,S}
    3 H  u0 {1,S}
    4 H  u0 {1,S}
    5 H  u0 {1,S}
    """)
    assert_true(methanei.is_subgraph_isomorphic(group_methane))
    assert_true(methane.is_subgraph_isomorphic(group_methane))
    assert_false(methanei.is_isomorphic(methane))

def test_isomorphism_wrong_mapping():
    """
    Checks isomorphism finds things not isomorphic if given a wrong mapping.
    """
    # These molecules are not isomorphic
    n1butane = Molecule().from_adjacency_list("""
        1 *1 X u0 p0 c0 {2,S}
        2 *2 C u0 p0 c0 {1,S} {3,S} {5,S} {6,S}
        3 *3 C u0 p0 c0 {2,S} {4,S} {7,S} {8,S}
        4 *4 H u0 p0 c0 {3,S}
        5    H u0 p0 c0 {2,S}
        6    H u0 p0 c0 {2,S}
        7    H u0 p0 c0 {3,S}
        8    C u0 p0 c0 {3,S} {9,S} {10,S} {11,S}
        9    H u0 p0 c0 {8,S}
        10   H u0 p0 c0 {8,S}
        11   C u0 p0 c0 {8,S} {12,S} {13,S} {14,S}
        12   H u0 p0 c0 {11,S}
        13   H u0 p0 c0 {11,S}
        14   H u0 p0 c0 {11,S}
        """)

    n2butane = Molecule().from_adjacency_list("""
        1 *1 X u0 p0 c0 {3,S}
        2 *2 C u0 p0 c0 {3,S} {5,S} {6,S} {4,S}
        3 *3 C u0 p0 c0 {2,S} {1,S} {7,S} {8,S}
        4 *4 H u0 p0 c0 {2,S}
        5    H u0 p0 c0 {2,S}
        6    H u0 p0 c0 {2,S}
        7    H u0 p0 c0 {3,S}
        8    C u0 p0 c0 {3,S} {9,S} {10,S} {11,S}
        9    H u0 p0 c0 {8,S}
        10   H u0 p0 c0 {8,S}
        11   C u0 p0 c0 {8,S} {12,S} {13,S} {14,S}
        12   H u0 p0 c0 {11,S}
        13   H u0 p0 c0 {11,S}
        14   H u0 p0 c0 {11,S}
        """)
    assert_false(n1butane.is_isomorphic(n2butane))
    assert_false(n1butane.is_isomorphic(n2butane, generate_initial_map=True))
    
    mapping = {}
    for label in ['*1', '*2', '*3', '*4']:
        mapping[n1butane.get_labeled_atoms(label)[0]] = n2butane.get_labeled_atoms(label)[0]
    assert_false(n1butane.is_isomorphic(n2butane, initial_map=mapping))
